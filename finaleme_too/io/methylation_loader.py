"""Multi-format methylation data loader.

Supported formats (architecture §3.2):
  - finaleme_bed     : FinaleMe prediction.bed.gz output
  - bissnp_6plus2    : Bis-SNP 6-column BED + 2 extra (methylation_pct, total_count)
  - wgbstools_beta   : Binary .beta file (NR_SITES x 2 uint8)
  - custom_bed       : User-defined column mapping via meth_col / total_col
"""

from __future__ import annotations

import gzip
import logging
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import pandas as pd

from finaleme_too.config import MeasurementMode
from finaleme_too.exceptions import InvalidInputFormatError
from finaleme_too.io.marker_regions import MarkerRegions

log = logging.getLogger(__name__)

_TABIX_FETCH_MAX_MERGE_GAP_BP = 250_000
_TABIX_FETCH_DENSE_COVERAGE_RATIO = 0.35


@dataclass(frozen=True)
class MarkerObservations:
    """Per-sample methylation observations aligned to marker regions.

    The chrom/start/end arrays represent the marker regions used for this
    observation. k and n are the per-marker methylated and total read counts.
    For FINALEME mode, predicted_beta optionally holds the raw FinaleMe
    prediction (uncalibrated) for each marker.

    ``called_state`` and ``context_bin`` are v3 FinaleMe binarization
    outputs populated by ``preprocessing.binarization.apply_binarization``.
    They default to ``None`` for WGBS mode and for any FinaleMe sample that
    has not yet been through the binarizer. The convention for
    ``called_state`` is a uint8 array where ``0 = U``, ``1 = M``,
    ``2 = Ambiguous``, ``3 = Excluded`` (bin not usable). ``context_bin``
    is an int32 array of per-marker bin indices in ``[0, n_bins)``.
    """

    sample_id: str
    chrom: np.ndarray
    start: np.ndarray
    end: np.ndarray
    k: np.ndarray  # int32
    n: np.ndarray  # int32
    predicted_beta: np.ndarray | None  # float32 or None
    mode: MeasurementMode
    called_state: np.ndarray | None = None  # uint8 (U=0, M=1, Ambiguous=2, Excluded=3)
    context_bin: np.ndarray | None = None  # int32 bin index

    def __len__(self) -> int:
        return len(self.k)

    @property
    def n_markers(self) -> int:
        return len(self.k)

    @property
    def total_reads(self) -> int:
        return int(np.sum(self.n))

    def with_counts(self, k: np.ndarray, n: np.ndarray) -> "MarkerObservations":
        return replace(self, k=k.astype(np.int32, copy=False), n=n.astype(np.int32, copy=False))

    def with_binarization(
        self,
        called_state: np.ndarray,
        context_bin: np.ndarray,
    ) -> "MarkerObservations":
        """Return a new ``MarkerObservations`` with binarization outputs filled in.

        v3 FinaleMe path: ``apply_binarization`` calls this to attach the
        ``called_state`` and ``context_bin`` arrays to a sample. All other
        fields (k, n, chrom, start, end, predicted_beta, mode, sample_id)
        are preserved bit-for-bit.
        """
        return replace(
            self,
            called_state=np.asarray(called_state, dtype=np.uint8),
            context_bin=np.asarray(context_bin, dtype=np.int32),
        )


class MethylationLoader:
    """Static loader: dispatches to format-specific parsers."""

    @staticmethod
    def load(
        filepath: str | Path,
        sample_id: str,
        mode: MeasurementMode,
        marker_regions: MarkerRegions,
        input_format: str | None = None,
        meth_col: int | None = None,
        total_col: int | None = None,
        cpg_index: dict | None = None,  # for wgbstools_beta
    ) -> MarkerObservations:
        path = Path(filepath)
        if not path.exists():
            raise InvalidInputFormatError(f"Methylation file not found: {path}")

        if input_format is None or input_format == "auto":
            input_format = MethylationLoader.auto_detect_format(path, mode)

        if input_format == "finaleme_bed":
            return MethylationLoader._parse_finaleme_bed(path, sample_id, marker_regions, mode)
        if input_format == "bissnp_6plus2":
            return MethylationLoader._parse_bissnp(path, sample_id, marker_regions, mode)
        if input_format == "wgbstools_beta":
            if cpg_index is None:
                raise InvalidInputFormatError(
                    "wgbstools_beta format requires --cpg-index"
                )
            return MethylationLoader._parse_wgbstools_beta(
                path, sample_id, marker_regions, cpg_index, mode
            )
        if input_format == "custom_bed":
            if meth_col is None or total_col is None:
                raise InvalidInputFormatError(
                    "custom_bed format requires --meth-col and --total-col (1-indexed)"
                )
            return MethylationLoader._parse_custom_bed(
                path, sample_id, marker_regions, meth_col, total_col, mode
            )
        raise InvalidInputFormatError(f"Unknown input_format: {input_format}")

    # ------------------------------------------------------------------
    # Format detection
    # ------------------------------------------------------------------

    @staticmethod
    def auto_detect_format(path: Path, mode: MeasurementMode) -> str:
        name = path.name.lower()
        if name.endswith(".beta") or name.endswith(".lbeta"):
            return "wgbstools_beta"
        if mode == MeasurementMode.FINALEME:
            return "finaleme_bed"

        # Peek at the first non-comment line to count columns
        opener = gzip.open if name.endswith(".gz") else open
        try:
            with opener(path, "rt") as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#") or line.startswith("track"):
                        continue
                    parts = line.split("\t")
                    if len(parts) >= 8:
                        try:
                            pct = float(parts[6])
                            if 0.0 <= pct <= 100.0:
                                return "bissnp_6plus2"
                        except ValueError:
                            pass
                    break
        except OSError:
            pass

        raise InvalidInputFormatError(
            f"Cannot auto-detect format for {path}. "
            "Specify --input-format or --meth-col/--total-col."
        )

    # ------------------------------------------------------------------
    # FinaleMe prediction.bed.gz parser
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_finaleme_bed(
        path: Path, sample_id: str, marker_regions: MarkerRegions, mode: MeasurementMode
    ) -> MarkerObservations:
        """Columns:
            #chr start end methy_perc_predict methy_count_predict total_count_predict
            methy_perc_obs methy_count_obs total_count_obs
        Per-CpG records are aggregated into the supplied marker regions.
        """
        # Fast path: tabix region queries. If unavailable, fall back to the
        # legacy full-file scan below.
        tabix_obs = _parse_finaleme_bed_with_tabix(path, sample_id, marker_regions, mode)
        if tabix_obs is not None:
            return tabix_obs

        opener = gzip.open if path.name.endswith(".gz") else open
        with opener(path, "rt") as fh:
            df = pd.read_csv(
                fh,
                sep="\t",
                comment="#",
                header=None,
                usecols=[0, 1, 2, 3, 4, 5],
                names=["chrom", "start", "end", "methy_pct_pred", "methy_count_pred", "total_count_pred"],
                dtype={
                    "chrom": str,
                    "start": np.int64,
                    "end": np.int64,
                    "methy_pct_pred": np.float64,
                    "methy_count_pred": np.float64,
                    "total_count_pred": np.float64,
                },
            )
        if df.empty:
            raise InvalidInputFormatError(f"Empty FinaleMe prediction file: {path}")
        return _aggregate_per_cpg_to_markers_first_match(
            sample_id=sample_id,
            cpg_chrom=df["chrom"].to_numpy(),
            cpg_start=df["start"].to_numpy(),
            cpg_methy=df["methy_count_pred"].to_numpy(),
            cpg_total=df["total_count_pred"].to_numpy(),
            marker_regions=marker_regions,
            mode=mode,
        )

    # ------------------------------------------------------------------
    # Bis-SNP 6+2 BED parser
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_bissnp(
        path: Path, sample_id: str, marker_regions: MarkerRegions, mode: MeasurementMode
    ) -> MarkerObservations:
        """Standard 6-column BED + col 7 = methylation_pct (0-100), col 8 = total_count.
        k = round(methylation_pct/100 * total_count)
        """
        opener = gzip.open if path.name.endswith(".gz") else open
        with opener(path, "rt") as fh:
            df = pd.read_csv(
                fh,
                sep="\t",
                comment="#",
                header=None,
                usecols=[0, 1, 2, 6, 7],
                names=["chrom", "start", "end", "methy_pct", "total"],
                dtype={
                    "chrom": str,
                    "start": np.int64,
                    "end": np.int64,
                    "methy_pct": np.float64,
                    "total": np.float64,
                },
            )
        if df.empty:
            raise InvalidInputFormatError(f"Empty Bis-SNP file: {path}")
        methy_count = np.round(df["methy_pct"].to_numpy() / 100.0 * df["total"].to_numpy())
        return _aggregate_per_cpg_to_markers(
            sample_id=sample_id,
            cpg_chrom=df["chrom"].to_numpy(),
            cpg_start=df["start"].to_numpy(),
            cpg_methy=methy_count,
            cpg_total=df["total"].to_numpy(),
            cpg_pct=None,
            marker_regions=marker_regions,
            mode=mode,
            keep_pct=False,
        )

    # ------------------------------------------------------------------
    # wgbstools .beta binary parser
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_wgbstools_beta(
        path: Path,
        sample_id: str,
        marker_regions: MarkerRegions,
        cpg_index: dict,
        mode: MeasurementMode,
    ) -> MarkerObservations:
        """Read a binary .beta file (NR_SITES x 2 uint8) and aggregate to marker regions.

        cpg_index has the structure produced by io.reference_panel._load_cpg_index:
            {
              "chr_positions": dict[str -> np.ndarray of sorted positions],
              "chr_offsets":   dict[str -> int],
              "total_sites":   int,
            }
        """
        with open(path, "rb") as fh:
            data = np.frombuffer(fh.read(), dtype=np.uint8)
        if data.size % 2 != 0:
            raise InvalidInputFormatError(
                f"Beta file size not a multiple of 2: {path}"
            )
        per_cpg = data.reshape((-1, 2)).astype(np.int32)  # cols: methy, total

        chr_positions = cpg_index["chr_positions"]
        chr_offsets = cpg_index["chr_offsets"]

        n_markers = marker_regions.n_markers
        k_arr = np.zeros(n_markers, dtype=np.int32)
        n_arr = np.zeros(n_markers, dtype=np.int32)

        for mi in range(n_markers):
            chrom = str(marker_regions.chrom[mi])
            start = int(marker_regions.start[mi])
            end = int(marker_regions.end[mi])
            positions = chr_positions.get(chrom)
            offset = chr_offsets.get(chrom)
            if positions is None or offset is None:
                continue
            lo = int(np.searchsorted(positions, start, side="left"))
            hi = int(np.searchsorted(positions, end, side="left"))
            if hi <= lo:
                continue
            global_lo = offset + lo
            global_hi = offset + hi
            if global_hi > per_cpg.shape[0]:
                global_hi = per_cpg.shape[0]
            if global_lo >= global_hi:
                continue
            block = per_cpg[global_lo:global_hi]
            k_arr[mi] = int(block[:, 0].sum())
            n_arr[mi] = int(block[:, 1].sum())

        return MarkerObservations(
            sample_id=sample_id,
            chrom=marker_regions.chrom,
            start=marker_regions.start,
            end=marker_regions.end,
            k=k_arr,
            n=n_arr,
            predicted_beta=None,
            mode=mode,
        )

    # ------------------------------------------------------------------
    # Custom BED parser
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_custom_bed(
        path: Path,
        sample_id: str,
        marker_regions: MarkerRegions,
        meth_col: int,
        total_col: int,
        mode: MeasurementMode,
    ) -> MarkerObservations:
        """meth_col and total_col are 1-indexed (matches CLI semantics)."""
        opener = gzip.open if path.name.endswith(".gz") else open
        meth_idx = meth_col - 1
        total_idx = total_col - 1
        with opener(path, "rt") as fh:
            df = pd.read_csv(
                fh,
                sep="\t",
                comment="#",
                header=None,
                usecols=[0, 1, 2, meth_idx, total_idx],
                names=["chrom", "start", "end", "methy", "total"],
            )
        return _aggregate_per_cpg_to_markers(
            sample_id=sample_id,
            cpg_chrom=df["chrom"].astype(str).to_numpy(),
            cpg_start=df["start"].to_numpy(),
            cpg_methy=df["methy"].to_numpy(),
            cpg_total=df["total"].to_numpy(),
            cpg_pct=None,
            marker_regions=marker_regions,
            mode=mode,
            keep_pct=False,
        )


def _parse_finaleme_bed_with_tabix(
    path: Path,
    sample_id: str,
    marker_regions: MarkerRegions,
    mode: MeasurementMode,
) -> MarkerObservations | None:
    """Attempt tabix-backed region loading for FinaleMe BED files.

    Returns ``None`` when tabix is unavailable or setup fails.
    """
    try:
        import pysam
    except Exception:
        return None

    tabix_path = _ensure_tabix_ready(path)
    if tabix_path is None:
        return None

    try:
        tbx = pysam.TabixFile(str(tabix_path))
    except Exception as exc:
        log.warning(
            "Failed to open tabix file %s (%s). Falling back to full-file parser.",
            tabix_path,
            exc,
        )
        return None

    try:
        contigs = set(tbx.contigs)
        n_markers = marker_regions.n_markers
        k_arr = np.zeros(n_markers, dtype=np.int64)
        n_arr = np.zeros(n_markers, dtype=np.int64)
        # Batch markers by chromosome, then do one tabix fetch per chromosome
        # (instead of one fetch per marker). This is much faster on network
        # filesystems than thousands of random tabix seeks.
        chrom_to_marker_idx: dict[str, list[int]] = {}
        for mi in range(n_markers):
            chrom = str(marker_regions.chrom[mi])
            chrom_to_marker_idx.setdefault(chrom, []).append(mi)

        for chrom, marker_ids_list in chrom_to_marker_idx.items():
            query_chrom = _resolve_query_chrom(chrom, contigs)
            if query_chrom is None:
                continue

            # Preserve original marker order for Java parity in overlap mode.
            marker_ids_orig = np.asarray(marker_ids_list, dtype=np.int64)
            starts_orig = np.asarray(marker_regions.start[marker_ids_orig], dtype=np.int64)
            ends_orig = np.asarray(marker_regions.end[marker_ids_orig], dtype=np.int64)
            if starts_orig.size == 0:
                continue
            valid_marker = ends_orig > starts_orig
            if not np.any(valid_marker):
                continue
            marker_ids_orig = marker_ids_orig[valid_marker]
            starts_orig = starts_orig[valid_marker]
            ends_orig = ends_orig[valid_marker]

            # Coordinate-sorted copies for fetch-window construction.
            order = np.argsort(starts_orig, kind="mergesort")
            marker_ids = marker_ids_orig[order]
            starts = starts_orig[order]
            ends = ends_orig[order]
            if starts.size <= 1:
                has_overlap = False
            else:
                prev_max_end = np.maximum.accumulate(ends[:-1])
                has_overlap = bool(np.any(starts[1:] < prev_max_end))

            # Java parity for overlapping markers:
            # assign each CpG to the first matching marker in original order.
            if has_overlap:
                fetch_windows = _build_tabix_fetch_windows(starts, ends)
                for _, _, fetch_start, fetch_end in fetch_windows:
                    if fetch_end <= fetch_start:
                        continue
                    cand_mask = (ends_orig > fetch_start) & (starts_orig < fetch_end)
                    if not np.any(cand_mask):
                        continue
                    cand_ids = marker_ids_orig[cand_mask]
                    cand_starts = starts_orig[cand_mask]
                    cand_ends = ends_orig[cand_mask]
                    try:
                        for line in tbx.fetch(query_chrom, int(fetch_start), int(fetch_end)):
                            parts = line.split("\t")
                            if len(parts) < 6:
                                continue
                            try:
                                cpg_start = int(parts[1])
                                meth = int(float(parts[4]))
                                total = int(float(parts[5]))
                            except ValueError:
                                continue
                            for idx in range(cand_ids.size):
                                if cand_starts[idx] <= cpg_start < cand_ends[idx]:
                                    mid = int(cand_ids[idx])
                                    k_arr[mid] += meth
                                    n_arr[mid] += total
                                    break
                    except ValueError:
                        pass
                continue

            # Non-overlapping markers: fast unique-assignment path.

            fetch_windows = _build_tabix_fetch_windows(starts, ends)
            for win_lo, win_hi, fetch_start, fetch_end in fetch_windows:
                marker_ids_w = marker_ids[win_lo:win_hi]
                starts_w = starts[win_lo:win_hi]
                ends_w = ends[win_lo:win_hi]
                if starts_w.size == 0 or fetch_end <= fetch_start:
                    continue

                ptr = 0
                n_local = starts_w.size
                try:
                    for line in tbx.fetch(query_chrom, fetch_start, fetch_end):
                        parts = line.split("\t")
                        if len(parts) < 6:
                            continue
                        try:
                            cpg_start = int(parts[1])
                            meth = float(parts[4])
                            total = float(parts[5])
                        except ValueError:
                            continue
                        # prediction.bed.gz is BED-like (0-based, half-open)
                        while ptr < n_local and ends_w[ptr] <= cpg_start:
                            ptr += 1
                        if ptr >= n_local:
                            break
                        if starts_w[ptr] <= cpg_start < ends_w[ptr]:
                            mid = int(marker_ids_w[ptr])
                            k_arr[mid] += int(meth)
                            n_arr[mid] += int(total)
                except ValueError:
                    pass
    finally:
        tbx.close()

    with np.errstate(invalid="ignore", divide="ignore"):
        predicted_beta = np.where(
            n_arr > 0, k_arr / np.maximum(n_arr, 1), np.nan
        ).astype(np.float32)

    return MarkerObservations(
        sample_id=sample_id,
        chrom=marker_regions.chrom,
        start=marker_regions.start,
        end=marker_regions.end,
        k=k_arr.astype(np.int32),
        n=n_arr.astype(np.int32),
        predicted_beta=predicted_beta,
        mode=mode,
    )


def _resolve_query_chrom(chrom: str, contigs: set[str]) -> str | None:
    """Resolve chr/no-chr naming differences for tabix queries."""
    if chrom in contigs:
        return chrom
    if chrom.startswith("chr"):
        alt = chrom[3:]
        if alt in contigs:
            return alt
    else:
        alt = f"chr{chrom}"
        if alt in contigs:
            return alt
    return None


def _build_tabix_fetch_windows(
    starts: np.ndarray,
    ends: np.ndarray,
    max_merge_gap_bp: int = _TABIX_FETCH_MAX_MERGE_GAP_BP,
    dense_coverage_ratio: float = _TABIX_FETCH_DENSE_COVERAGE_RATIO,
) -> list[tuple[int, int, int, int]]:
    """Build fetch windows from sorted marker intervals for one chromosome.

    Returns a list of tuples:
        (lo, hi, fetch_start, fetch_end)
    where ``lo:hi`` slice the marker arrays for that window.

    Strategy:
      * dense marker layout -> one broad fetch (fewer tabix seeks)
      * sparse marker layout -> clustered fetch windows (avoid scanning
        large chromosome spans with no markers)
    """
    if starts.size == 0:
        return []

    # Dense layout: one broad fetch is cheaper than many random seeks.
    chrom_span = int(max(1, int(np.max(ends)) - int(starts[0])))
    marker_bp = int(np.sum(np.maximum(0, ends - starts)))
    if starts.size <= 2 or (marker_bp / chrom_span) >= dense_coverage_ratio:
        fetch_start = int(max(0, starts[0]))
        fetch_end = int(max(fetch_start, np.max(ends)))
        return [(0, int(starts.size), fetch_start, fetch_end)]

    # Sparse layout: cluster nearby markers into local windows.
    windows: list[tuple[int, int, int, int]] = []
    lo = 0
    cur_end = int(ends[0])
    for i in range(1, int(starts.size)):
        s = int(starts[i])
        e = int(ends[i])
        if s <= cur_end + int(max_merge_gap_bp):
            if e > cur_end:
                cur_end = e
            continue
        fetch_start = int(max(0, starts[lo]))
        fetch_end = int(max(fetch_start, cur_end))
        windows.append((lo, i, fetch_start, fetch_end))
        lo = i
        cur_end = e
    fetch_start = int(max(0, starts[lo]))
    fetch_end = int(max(fetch_start, cur_end))
    windows.append((lo, int(starts.size), fetch_start, fetch_end))
    return windows


def _has_tabix_index(path: Path) -> bool:
    return path.with_suffix(path.suffix + ".tbi").exists() or path.with_suffix(
        path.suffix + ".csi"
    ).exists()


def _ensure_tabix_ready(path: Path) -> Path | None:
    """Return a tabix-ready path, creating a cached bgzip+tabix copy if needed."""
    try:
        import pysam
    except Exception:
        return None

    # Already indexed and readable.
    if _has_tabix_index(path):
        try:
            with pysam.TabixFile(str(path)):
                return path
        except Exception:
            pass

    # Try indexing in place (works when file is already sorted bgzip).
    try:
        pysam.tabix_index(str(path), preset="bed", force=True, keep_original=True)
        with pysam.TabixFile(str(path)):
            return path
    except Exception:
        pass

    # Otherwise build a sorted bgzip+tabix cache.
    cache_path = _build_tabix_cache(path)
    if cache_path is None:
        return None
    if _has_tabix_index(cache_path):
        try:
            with pysam.TabixFile(str(cache_path)):
                return cache_path
        except Exception:
            return None
    return None


def _build_tabix_cache(path: Path) -> Path | None:
    """Build a sorted bgzip/tabix cache copy from an arbitrary FinaleMe BED source."""
    try:
        import hashlib
        import pysam
    except Exception:
        return None

    stat = path.stat()
    sig = hashlib.sha1(
        f"{path.resolve()}::{stat.st_size}::{stat.st_mtime_ns}".encode("utf-8")
    ).hexdigest()[:12]
    cache_dir = path.parent / ".finaleme_too_tabix_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    out_bgz = cache_dir / f"{path.stem}.{sig}.prediction.bed.gz"

    if _has_tabix_index(out_bgz):
        return out_bgz

    tmp_plain = cache_dir / f"{path.stem}.{sig}.prediction.bed"
    opener = gzip.open if path.name.endswith(".gz") else open
    rows: list[tuple[str, int, int, float, float, float]] = []

    try:
        with opener(path, "rt") as fh:
            for line in fh:
                s = line.strip()
                if not s or s.startswith("#") or s.startswith("track"):
                    continue
                parts = s.split("\t")
                if len(parts) < 6:
                    continue
                try:
                    rows.append(
                        (
                            str(parts[0]),
                            int(parts[1]),
                            int(parts[2]),
                            float(parts[3]),
                            float(parts[4]),
                            float(parts[5]),
                        )
                    )
                except ValueError:
                    continue
    except OSError as exc:
        log.warning("Failed reading %s for tabix cache build (%s).", path, exc)
        return None

    if not rows:
        return None

    rows.sort(key=lambda x: (x[0], x[1], x[2]))
    with open(tmp_plain, "wt") as out:
        out.write(
            "#chr\tstart\tend\tmethy_perc_predict\tmethy_count_predict\ttotal_count_predict\n"
        )
        for chrom, start, end, pct, meth, total in rows:
            out.write(f"{chrom}\t{start}\t{end}\t{pct}\t{meth}\t{total}\n")

    try:
        pysam.tabix_compress(str(tmp_plain), str(out_bgz), force=True)
        pysam.tabix_index(str(out_bgz), preset="bed", force=True)
        log.info("Built tabix cache for %s at %s", path, out_bgz)
        return out_bgz
    except Exception as exc:
        log.warning("Failed building tabix cache for %s (%s).", path, exc)
        return None
    finally:
        try:
            tmp_plain.unlink(missing_ok=True)
        except OSError:
            pass


def _aggregate_per_cpg_to_markers(
    sample_id: str,
    cpg_chrom: np.ndarray,
    cpg_start: np.ndarray,
    cpg_methy: np.ndarray,
    cpg_total: np.ndarray,
    cpg_pct: np.ndarray | None,
    marker_regions: MarkerRegions,
    mode: MeasurementMode,
    keep_pct: bool,
) -> MarkerObservations:
    """Aggregate per-CpG methylation records into the supplied marker regions.

    Sorts CpGs by (chrom, start) once, then for each marker region uses
    binary search to find the CpG range. The same approach is used in the Java
    BetaValueDeconvolution._loadQueryFromPrediction method.
    """
    n_markers = marker_regions.n_markers
    k_arr = np.zeros(n_markers, dtype=np.int64)
    n_arr = np.zeros(n_markers, dtype=np.int64)

    cpg_chrom_arr = np.asarray(cpg_chrom)
    cpg_start_arr = np.asarray(cpg_start, dtype=np.int64)
    cpg_methy_arr = np.asarray(cpg_methy, dtype=np.float64)
    cpg_total_arr = np.asarray(cpg_total, dtype=np.float64)
    # cpg_pct (per-CpG percentage) is no longer used: predicted_beta is now
    # derived from the count-weighted ratio so it stays consistent with k/n.

    # Group CpGs by chromosome
    unique_chroms = np.unique(cpg_chrom_arr)
    chrom_to_indices: dict[str, np.ndarray] = {}
    for c in unique_chroms:
        idx = np.where(cpg_chrom_arr == c)[0]
        # Sort by start within chromosome
        idx_sorted = idx[np.argsort(cpg_start_arr[idx], kind="stable")]
        chrom_to_indices[str(c)] = idx_sorted

    for mi in range(n_markers):
        chrom = str(marker_regions.chrom[mi])
        idx_for_chrom = chrom_to_indices.get(chrom)
        if idx_for_chrom is None or idx_for_chrom.size == 0:
            continue
        start = int(marker_regions.start[mi])
        end = int(marker_regions.end[mi])
        positions = cpg_start_arr[idx_for_chrom]
        lo = int(np.searchsorted(positions, start, side="left"))
        hi = int(np.searchsorted(positions, end, side="left"))
        if hi <= lo:
            continue
        sel = idx_for_chrom[lo:hi]
        k_arr[mi] = int(np.sum(cpg_methy_arr[sel]))
        n_arr[mi] = int(np.sum(cpg_total_arr[sel]))

    if keep_pct:
        # Derive predicted_beta from the count-weighted ratio so it stays
        # consistent with k/n. An unweighted mean of per-CpG percentages
        # would diverge from k/n whenever CpGs have different coverage.
        with np.errstate(invalid="ignore", divide="ignore"):
            predicted_beta = np.where(
                n_arr > 0, k_arr / np.maximum(n_arr, 1), np.nan
            ).astype(np.float32)
    else:
        predicted_beta = None

    return MarkerObservations(
        sample_id=sample_id,
        chrom=marker_regions.chrom,
        start=marker_regions.start,
        end=marker_regions.end,
        k=k_arr.astype(np.int32),
        n=n_arr.astype(np.int32),
        predicted_beta=predicted_beta,
        mode=mode,
    )


def _aggregate_per_cpg_to_markers_first_match(
    sample_id: str,
    cpg_chrom: np.ndarray,
    cpg_start: np.ndarray,
    cpg_methy: np.ndarray,
    cpg_total: np.ndarray,
    marker_regions: MarkerRegions,
    mode: MeasurementMode,
) -> MarkerObservations:
    """Java-compatible aggregation for overlapping markers (first-hit assignment).

    For each CpG, assign counts to the first matching marker interval in the
    original marker order for that chromosome, then stop.
    """
    n_markers = marker_regions.n_markers
    k_arr = np.zeros(n_markers, dtype=np.int64)
    n_arr = np.zeros(n_markers, dtype=np.int64)

    cpg_chrom_arr = np.asarray(cpg_chrom, dtype=object)
    cpg_start_arr = np.asarray(cpg_start, dtype=np.int64)
    cpg_methy_arr = np.asarray(cpg_methy, dtype=np.float64)
    cpg_total_arr = np.asarray(cpg_total, dtype=np.float64)

    chrom_to_marker_idx: dict[str, list[int]] = {}
    for mi in range(n_markers):
        chrom = str(marker_regions.chrom[mi])
        chrom_to_marker_idx.setdefault(chrom, []).append(mi)

    for i in range(cpg_start_arr.size):
        chrom = str(cpg_chrom_arr[i])
        marker_ids = chrom_to_marker_idx.get(chrom)
        if not marker_ids:
            continue
        pos = int(cpg_start_arr[i])
        for mi in marker_ids:
            start = int(marker_regions.start[mi])
            end = int(marker_regions.end[mi])
            if start <= pos < end:
                k_arr[mi] += int(cpg_methy_arr[i])
                n_arr[mi] += int(cpg_total_arr[i])
                break

    with np.errstate(invalid="ignore", divide="ignore"):
        predicted_beta = np.where(
            n_arr > 0, k_arr / np.maximum(n_arr, 1), np.nan
        ).astype(np.float32)

    return MarkerObservations(
        sample_id=sample_id,
        chrom=marker_regions.chrom,
        start=marker_regions.start,
        end=marker_regions.end,
        k=k_arr.astype(np.int32),
        n=n_arr.astype(np.int32),
        predicted_beta=predicted_beta,
        mode=mode,
    )


__all__ = ["MarkerObservations", "MethylationLoader"]
