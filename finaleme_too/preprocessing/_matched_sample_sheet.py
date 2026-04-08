"""Shared parser helpers for matched WGBS / FinaleMe training data.

Used by both the v2 calibration path (``preprocessing/calibration.py``) and
the v3 binarization path (``preprocessing/binarization.py``). The module
contains no calibration- or binarization-specific logic — only:

  * Chromosome-name normalization
  * Sample-sheet + legacy-joined-TSV loading (thread-parallel)
  * Bis-SNP 6+2 BED and FinaleMe prediction.bed.gz parsers
  * CpG index loading
  * Per-row CpG density computation
  * Per-row region class computation (CGI / shore / shelf / open_sea)
  * Region-annotation TSV builder

During Phase A of the v2→v3 migration these functions live as COPIES of the
original calibration.py implementations. In Phase E the original definitions
are deleted and calibration.py (if still present) imports from here. In the
meantime both modules independently expose the same behavior so neither path
regresses.

The new ``_classify_region`` helper and the ``region_class`` column on
``compute_region_annotation`` are v3 additions that do not exist in
calibration.py. They are gated behind a default-off parameter so any v2
caller importing this module sees identical output to the original.
"""

from __future__ import annotations

import gzip
import io
from pathlib import Path

import numpy as np
import pandas as pd

from finaleme_too.exceptions import InvalidCalibrationError, InvalidMatchedDataError

_CANONICAL_COLUMNS = [
    "sample_id",
    "chrom",
    "start",
    "end",
    "methylated_count",
    "total_count",
]


# ---------------------------------------------------------------------------
# Chromosome normalization + matched-table loading
# ---------------------------------------------------------------------------


def _normalize_chrom(df: pd.DataFrame) -> pd.DataFrame:
    """Strip any leading ``chr`` from the chromosome column.

    Bis-SNP often writes chromosomes in the GRCh37 convention without a
    ``chr`` prefix (``1``, ``2``, ``X``), while FinaleMe output always has
    the prefix (``chr1``, ``chr2``, ``chrX``). Strip it on both sides so
    the (sample_id, chrom, start, end) join across modalities works.
    """
    df = df.copy()
    df["chrom"] = (
        df["chrom"].astype(str).str.replace(r"^chr", "", regex=True)
    )
    return df


def _load_matched_table(
    path: str | Path,
    modality: str = "wgbs",
    threads: int = 1,
) -> pd.DataFrame:
    """Load a matched WGBS or FinaleMe training table.

    Three formats are supported (auto-detected from the header row):

    1. **Legacy joined TSV** — columns ``sample_id chrom start end
       methylated_count total_count``.
    2. **Sample sheet → Bis-SNP 6+2 BEDs** — columns ``sample_id
       methylation_file [format]``. ``format`` defaults to ``bissnp_6plus2``
       when ``modality="wgbs"``.
    3. **Sample sheet → FinaleMe prediction BEDs** — same columns; ``format``
       defaults to ``finaleme_bed`` when ``modality="finaleme"``.

    Chromosome prefixes are always stripped before returning. When
    ``threads > 1`` and the input is a sample sheet, per-sample files are
    parsed in parallel via joblib.
    """
    df = pd.read_csv(path, sep="\t", comment="#")
    cols = set(df.columns)

    # Legacy format
    legacy_required = set(_CANONICAL_COLUMNS)
    if legacy_required.issubset(cols):
        return _normalize_chrom(df)

    # Sample sheet format
    sheet_required = {"sample_id", "methylation_file"}
    if sheet_required.issubset(cols):
        return _load_sample_sheet(df, Path(path).parent, modality, threads=threads)

    raise InvalidMatchedDataError(
        f"{path}: table has neither legacy columns {sorted(legacy_required)} "
        f"nor sample sheet columns {sorted(sheet_required)}; got {sorted(cols)}"
    )


def _parse_one_sample_row(
    row_dict: dict,
    base_dir: Path,
    default_format: str,
    has_format_col: bool,
) -> pd.DataFrame:
    """Worker: parse a single sample-sheet row into the canonical layout."""
    sample_id = str(row_dict["sample_id"])
    file_path = Path(str(row_dict["methylation_file"]))
    if not file_path.is_absolute():
        file_path = base_dir / file_path
    if not file_path.exists():
        raise InvalidMatchedDataError(
            f"Sample {sample_id}: file not found: {file_path}"
        )

    file_format = None
    if has_format_col:
        raw_fmt = row_dict.get("format")
        if raw_fmt is not None and pd.notna(raw_fmt) and str(raw_fmt).strip():
            file_format = str(raw_fmt).strip().lower()
    if file_format is None:
        file_format = default_format

    if file_format == "bissnp_6plus2":
        return _parse_bissnp_6plus2(file_path, sample_id)
    if file_format == "finaleme_bed":
        return _parse_finaleme_prediction(file_path, sample_id)
    raise InvalidMatchedDataError(
        f"Sample {sample_id}: unknown file format '{file_format}'. "
        "Supported: bissnp_6plus2, finaleme_bed"
    )


def _load_sample_sheet(
    df: pd.DataFrame,
    base_dir: Path,
    modality: str,
    threads: int = 1,
) -> pd.DataFrame:
    """Iterate a sample sheet and concatenate per-sample parsed tables.

    Parses per-sample files in parallel when ``threads > 1``. Parsing is
    the dominant cost on large cohorts (millions of CpGs per file through
    pandas.read_csv), so this is where parallelism pays off the most.
    """
    from finaleme_too.utils.parallel import parallel_map

    if df.empty:
        raise InvalidMatchedDataError("Sample sheet had no rows")

    default_format = _default_format_for_modality(modality)
    has_format_col = "format" in df.columns
    row_dicts = df.to_dict("records")

    parts = parallel_map(
        lambda row_dict: _parse_one_sample_row(
            row_dict, base_dir, default_format, has_format_col
        ),
        row_dicts,
        n_jobs=max(1, int(threads)),
        backend="threading",
    )
    combined = pd.concat(parts, ignore_index=True)
    return _normalize_chrom(combined)


def _default_format_for_modality(modality: str) -> str:
    mod = (modality or "").lower()
    if mod == "finaleme":
        return "finaleme_bed"
    return "bissnp_6plus2"


def _open_text(path: Path):
    """Open a file for text reading, decompressing on-the-fly if gzipped."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


# ---------------------------------------------------------------------------
# Per-sample file parsers
# ---------------------------------------------------------------------------


def _parse_bissnp_6plus2(path: Path, sample_id: str) -> pd.DataFrame:
    """Parse a Bis-SNP 6+2 BED file into the canonical six-column layout.

    Bis-SNP 6+2 format:
        chrom  start  end  name  score  strand  methylation_pct  total_count

    ``methylated_count = round(methylation_pct/100 * total_count)``, clamped
    to ``[0, total_count]``.
    """
    data_lines: list[str] = []
    with _open_text(path) as fh:
        for line in fh:
            if (
                not line.strip()
                or line.startswith("track")
                or line.startswith("browser")
                or line.startswith("#")
            ):
                continue
            data_lines.append(line)
    if not data_lines:
        return pd.DataFrame(columns=_CANONICAL_COLUMNS)

    df = pd.read_csv(
        io.StringIO("".join(data_lines)),
        sep="\t",
        header=None,
        usecols=[0, 1, 2, 6, 7],
        names=["chrom", "start", "end", "methy_pct", "total_count"],
        dtype={
            "chrom": str,
            "start": np.int64,
            "end": np.int64,
            "methy_pct": np.float64,
            "total_count": np.int64,
        },
    )
    methy_count = np.round(
        df["methy_pct"].to_numpy() / 100.0 * df["total_count"].to_numpy()
    ).astype(np.int64)
    methy_count = np.clip(methy_count, 0, df["total_count"].to_numpy())
    out = pd.DataFrame(
        {
            "sample_id": sample_id,
            "chrom": df["chrom"],
            "start": df["start"],
            "end": df["end"],
            "methylated_count": methy_count,
            "total_count": df["total_count"],
        }
    )
    return out


def _parse_finaleme_prediction(path: Path, sample_id: str) -> pd.DataFrame:
    """Parse a FinaleMe ``prediction.bed.gz`` file into the canonical layout.

    Expected header:
        #chr start end methy_perc_predict methy_count_predict total_count_predict
            methy_perc_obs methy_count_obs total_count_obs

    Clamps ``0 <= methy <= total`` defensively — malformed rows with
    ``methy_count_predict > total_count_predict`` would otherwise produce
    beta > 1 and poison downstream training.
    """
    with _open_text(path) as fh:
        df = pd.read_csv(
            fh,
            sep="\t",
            comment="#",
            header=None,
            usecols=[0, 1, 2, 4, 5],
            names=[
                "chrom", "start", "end",
                "methy_count_predict", "total_count_predict",
            ],
            dtype={
                "chrom": str,
                "start": np.int64,
                "end": np.int64,
                "methy_count_predict": np.float64,
                "total_count_predict": np.float64,
            },
        )
    methy = df["methy_count_predict"].round().astype(np.int64).to_numpy()
    total = df["total_count_predict"].round().astype(np.int64).to_numpy()
    total = np.maximum(total, 0)
    methy = np.clip(methy, 0, total)
    out = pd.DataFrame(
        {
            "sample_id": sample_id,
            "chrom": df["chrom"],
            "start": df["start"],
            "end": df["end"],
            "methylated_count": methy,
            "total_count": total,
        }
    )
    return out


# ---------------------------------------------------------------------------
# CpG index + region annotation
# ---------------------------------------------------------------------------


def _load_cpg_index_by_chrom(cpg_index_path: Path) -> dict[str, np.ndarray]:
    """Return ``{chrom (no prefix) -> sorted np.int64 array of CpG starts}``.

    Accepts the same CpG index BED format used elsewhere in the package:
    tab-separated ``chrom start end`` with optional gzip. Chromosome names
    are stripped of any leading ``chr``.
    """
    chr_to_pos: dict[str, list[int]] = {}
    with _open_text(cpg_index_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            chrom = parts[0]
            if chrom.startswith("chr"):
                chrom = chrom[3:]
            chr_to_pos.setdefault(chrom, []).append(pos)
    out: dict[str, np.ndarray] = {}
    for chrom, positions in chr_to_pos.items():
        arr = np.asarray(sorted(set(positions)), dtype=np.int64)
        out[chrom] = arr
    return out


# Default CpG density thresholds for the four standard region classes.
# These are approximate and configurable via BinarizationConfig at training
# time; they encode the "automatic" classification described in
# architecture doc §3.5.
DEFAULT_REGION_CLASS_THRESHOLDS: dict[str, float] = {
    # Class name   -> minimum cpg_density to qualify for THIS class
    "CGI": 0.08,       # dense CpG clusters → CGI
    "shore": 0.04,     # 2–8% density → shore
    "shelf": 0.015,    # 1.5–4% density → shelf
    # open_sea is the implicit bottom class (density < shelf threshold)
}

# Canonical ordering used when assigning a bin index or when a downstream
# consumer needs a stable column order. Keep CGI first so region_class == 0
# corresponds to the highest density class, matching how v2 bin_edges were
# ordered (dense → sparse).
REGION_CLASS_ORDER: tuple[str, ...] = ("CGI", "shore", "shelf", "open_sea")


def _classify_region(
    cpg_density: np.ndarray,
    thresholds: dict[str, float] | None = None,
) -> np.ndarray:
    """Classify each marker into ``{CGI, shore, shelf, open_sea}`` by density.

    v3.0 uses CpG density as the sole predictor of region class. Markers
    with density ≥ ``thresholds["CGI"]`` are CGI; markers below the lowest
    threshold are open_sea; the two bands in between are shore and shelf.
    See the architecture doc §3.5 "region_class ∈ {CGI, shore, shelf,
    open_sea}".

    NaN densities are labeled "open_sea" deterministically so downstream
    bin assignment never fails on missing annotations.
    """
    if thresholds is None:
        thresholds = DEFAULT_REGION_CLASS_THRESHOLDS
    density = np.asarray(cpg_density, dtype=np.float64)
    # Build output array of object dtype so we can write string labels.
    out = np.full(density.shape, "open_sea", dtype=object)
    # Treat NaN as open_sea (handled by the default fill above).
    finite = np.isfinite(density)
    out[finite & (density >= thresholds.get("shelf", 0.015))] = "shelf"
    out[finite & (density >= thresholds.get("shore", 0.04))] = "shore"
    out[finite & (density >= thresholds.get("CGI", 0.08))] = "CGI"
    return out


def compute_region_annotation(
    regions: pd.DataFrame,
    cpg_index_path: str | Path,
    window: int = 1000,
    region_class_thresholds: dict[str, float] | None = None,
) -> pd.DataFrame:
    """Compute per-row CpG density and region class from a CpG index BED.

    For each ``(chrom, start, end)`` row in ``regions``, density is
    ``count of CpGs in [center - window/2, center + window/2) / window``.
    The ``region_class`` column is derived from the density via
    :func:`_classify_region` (CGI / shore / shelf / open_sea).

    Parameters
    ----------
    regions
        Any DataFrame with ``chrom``, ``start``, ``end`` columns.
    cpg_index_path
        Path to a CpG index BED file (tab-separated, optionally gzipped).
    window
        Window size in base pairs centered on each row. Default 1000.
    region_class_thresholds
        Optional override for the default CpG-density thresholds
        ``{"CGI": 0.08, "shore": 0.04, "shelf": 0.015}``.

    Returns
    -------
    DataFrame with columns ``chrom start end cpg_density region_class`` —
    same row order and same row count as the input, chromosome prefix
    stripped. The v2 consumers of this function read only
    ``chrom, start, end, cpg_density`` and ignore the extra column, so
    adding ``region_class`` is backwards-compatible.
    """
    if window <= 0:
        raise InvalidCalibrationError(f"window must be positive, got {window}")
    if not {"chrom", "start", "end"}.issubset(regions.columns):
        raise InvalidCalibrationError(
            "regions DataFrame must have chrom, start, end columns"
        )

    cpg_by_chrom = _load_cpg_index_by_chrom(Path(cpg_index_path))
    if not cpg_by_chrom:
        raise InvalidCalibrationError(f"empty CpG index: {cpg_index_path}")

    chroms = regions["chrom"].astype(str).str.replace(r"^chr", "", regex=True)
    starts = regions["start"].astype(np.int64).to_numpy()
    ends = regions["end"].astype(np.int64).to_numpy()
    centers = (starts + ends) // 2
    half = window // 2
    lows = centers - half
    highs = centers + (window - half)

    density = np.zeros(len(regions), dtype=np.float64)
    unique_chroms = chroms.unique()
    for chrom in unique_chroms:
        mask = (chroms == chrom).to_numpy()
        cpg_arr = cpg_by_chrom.get(chrom)
        if cpg_arr is None or cpg_arr.size == 0:
            continue
        lo_sub = lows[mask]
        hi_sub = highs[mask]
        lo_idx = np.searchsorted(cpg_arr, lo_sub, side="left")
        hi_idx = np.searchsorted(cpg_arr, hi_sub, side="left")
        density[mask] = (hi_idx - lo_idx).astype(np.float64) / float(window)

    region_class = _classify_region(density, region_class_thresholds)

    return pd.DataFrame(
        {
            "chrom": chroms.to_numpy(),
            "start": starts,
            "end": ends,
            "cpg_density": density,
            "region_class": region_class,
        }
    )


def make_region_annotation_from_regions_file(
    regions_path: str | Path,
    cpg_index_path: str | Path,
    output_path: str | Path,
    window: int = 1000,
    region_class_thresholds: dict[str, float] | None = None,
) -> None:
    """Build a ``region_annotation.tsv`` for a BED/TSV of regions.

    The output TSV has columns ``chrom start end cpg_density region_class``
    and is directly consumable by ``finaleme-too train-binarization
    --region-annotation`` (and, for backwards compatibility, by
    ``finaleme-too train-calibration --region-annotation``, which reads
    only the ``cpg_density`` column).
    """
    path = Path(regions_path)
    first_data_line: str | None = None
    with _open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            first_data_line = line
            break
    has_header = False
    if first_data_line is not None:
        try:
            int(first_data_line.split("\t")[1])
        except (ValueError, IndexError):
            has_header = True

    if has_header:
        df = pd.read_csv(path, sep="\t", comment="#")
    else:
        df = pd.read_csv(
            path,
            sep="\t",
            comment="#",
            header=None,
            usecols=[0, 1, 2],
            names=["chrom", "start", "end"],
            dtype={"chrom": str, "start": np.int64, "end": np.int64},
        )
    if df.empty:
        raise InvalidCalibrationError(f"no regions found in {regions_path}")

    annotation = compute_region_annotation(
        df,
        cpg_index_path=cpg_index_path,
        window=window,
        region_class_thresholds=region_class_thresholds,
    )
    annotation.to_csv(output_path, sep="\t", index=False)


__all__ = [
    "_CANONICAL_COLUMNS",
    "DEFAULT_REGION_CLASS_THRESHOLDS",
    "REGION_CLASS_ORDER",
    "_classify_region",
    "_default_format_for_modality",
    "_load_cpg_index_by_chrom",
    "_load_matched_table",
    "_load_sample_sheet",
    "_normalize_chrom",
    "_open_text",
    "_parse_bissnp_6plus2",
    "_parse_finaleme_prediction",
    "_parse_one_sample_row",
    "compute_region_annotation",
    "make_region_annotation_from_regions_file",
]
