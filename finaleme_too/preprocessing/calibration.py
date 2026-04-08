"""FinaleMe calibration: apply (Phase B) + train (Phase C).

Math doc §6.1:
    logit(mu_calibrated) = a_b * logit(mu_FinaleMe) + c_b
    phi_b = exp(d_b)

Per CpG-density bin b. Parameters loaded from a JSON config file (or shipped
default). Phase C implements the training pipeline that produces such files.
"""

from __future__ import annotations

import gzip
import io
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from finaleme_too.exceptions import InvalidCalibrationError
from finaleme_too.io.methylation_loader import MarkerObservations


@dataclass
class CalibrationParams:
    """Per-bin calibration parameters."""

    n_bins: int
    bin_edges: np.ndarray  # shape (n_bins+1,) — CpG-density bin boundaries
    a: np.ndarray  # shape (n_bins,) — slopes
    c: np.ndarray  # shape (n_bins,) — intercepts
    log_dispersion: np.ndarray  # shape (n_bins,) — log phi per bin
    training_metadata: dict = field(default_factory=dict)

    @classmethod
    def from_dict(cls, raw: dict) -> "CalibrationParams":
        try:
            return cls(
                n_bins=int(raw["n_bins"]),
                bin_edges=np.asarray(raw["bin_edges"], dtype=np.float64),
                a=np.asarray(raw["a"], dtype=np.float64),
                c=np.asarray(raw["c"], dtype=np.float64),
                log_dispersion=np.asarray(raw["log_dispersion"], dtype=np.float64),
                training_metadata=dict(raw.get("training_metadata", {})),
            )
        except KeyError as exc:
            raise InvalidCalibrationError(
                f"Calibration JSON missing required key: {exc}"
            ) from exc

    @classmethod
    def load(cls, path: str | Path) -> "CalibrationParams":
        p = Path(path)
        if not p.exists():
            raise InvalidCalibrationError(f"Calibration file not found: {p}")
        with open(p) as fh:
            raw = json.load(fh)
        return cls.from_dict(raw)

    def save(self, path: str | Path) -> None:
        out = {
            "n_bins": self.n_bins,
            "bin_edges": self.bin_edges.tolist(),
            "a": self.a.tolist(),
            "c": self.c.tolist(),
            "log_dispersion": self.log_dispersion.tolist(),
            "training_metadata": self.training_metadata,
        }
        Path(path).write_text(json.dumps(out, indent=2))

    def assign_bin(self, density: np.ndarray) -> np.ndarray:
        """Map per-marker CpG density to a bin index in [0, n_bins-1].

        NaN densities are deterministically assigned to ``fallback_bin``.
        """
        density_arr = np.asarray(density, dtype=np.float64)
        nan_mask = ~np.isfinite(density_arr)
        clean = np.where(nan_mask, 0.0, density_arr)  # placeholder, overwritten
        idx = np.clip(
            np.searchsorted(self.bin_edges, clean, side="right") - 1,
            0,
            self.n_bins - 1,
        )
        idx = np.where(nan_mask, self.fallback_bin, idx)
        return idx.astype(np.int64)

    @property
    def fallback_bin(self) -> int:
        """Bin index used for markers without a known CpG density.

        Picks the bin whose midpoint (between its finite edges) is closest
        to the median of the *finite* training densities. When all edges
        are infinite (degenerate B=1 case), falls back to bin 0.
        """
        finite_edges = self.bin_edges[np.isfinite(self.bin_edges)]
        if finite_edges.size == 0:
            return 0
        # Use the median of the finite interior edges as a deterministic
        # representative density and look up its bin.
        rep = float(np.median(finite_edges))
        idx = int(
            np.clip(
                np.searchsorted(self.bin_edges, rep, side="right") - 1,
                0,
                self.n_bins - 1,
            )
        )
        return idx


# ----------------------------------------------------------------------------
# Apply path
# ----------------------------------------------------------------------------


def _logit(p: np.ndarray) -> np.ndarray:
    p_clip = np.clip(p, 1e-6, 1.0 - 1e-6)
    return np.log(p_clip / (1.0 - p_clip))


def _expit(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-x))


def apply_calibration(
    obs: MarkerObservations,
    params: CalibrationParams,
    region_annotations: pd.DataFrame | None,
) -> MarkerObservations:
    """Apply per-bin calibration to a sample's predicted methylation.

    Inputs
    ------
    obs : MarkerObservations
        FinaleMe-mode observations with predicted_beta populated.
    params : CalibrationParams
        Pre-trained per-bin parameters.
    region_annotations : DataFrame
        Must contain columns ``chrom, start, end, cpg_density``.

    Returns
    -------
    MarkerObservations with k recomputed from the calibrated beta value
    (k_calibrated = round(beta_calibrated * n)). n is preserved.
    """
    if obs.predicted_beta is None:
        # Nothing to calibrate
        return obs

    n = np.asarray(obs.n, dtype=np.int64)
    pred = np.asarray(obs.predicted_beta, dtype=np.float64)

    # Build a per-marker density vector by joining on (chrom, start, end).
    # Markers without a known density become NaN — assign_bin maps those
    # deterministically to params.fallback_bin (no NaN-mean bug).
    if region_annotations is not None and not region_annotations.empty:
        ann = region_annotations.set_index(["chrom", "start", "end"])["cpg_density"]
        keys = list(zip(obs.chrom.tolist(), obs.start.tolist(), obs.end.tolist()))
        density = np.array([float(ann.get(k, np.nan)) for k in keys], dtype=np.float64)
    else:
        density = np.full_like(pred, np.nan, dtype=np.float64)

    bin_idx = params.assign_bin(density)
    a = params.a[bin_idx]
    c = params.c[bin_idx]

    calibrated_logit = a * _logit(pred) + c
    calibrated = _expit(calibrated_logit)

    new_k = np.round(calibrated * n).astype(np.int32)
    new_k = np.clip(new_k, 0, n.astype(np.int32))
    return obs.with_counts(new_k, n.astype(np.int32))


def shipped_default_calibration_path() -> Path:
    """Path to the default calibration JSON shipped with the package."""
    from finaleme_too import data as data_pkg

    pkg_dir = Path(data_pkg.__file__).parent
    return pkg_dir / "default_calibration.json"


def load_default_calibration() -> CalibrationParams:
    return CalibrationParams.load(shipped_default_calibration_path())


# ----------------------------------------------------------------------------
# Phase C: training pipeline
# ----------------------------------------------------------------------------


_CANONICAL_COLUMNS = [
    "sample_id", "chrom", "start", "end", "methylated_count", "total_count",
]


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
) -> pd.DataFrame:
    """Load a matched WGBS or FinaleMe table for calibration training.

    Three formats are supported (auto-detected from the header row):

    1. **Legacy joined TSV** — columns ``sample_id chrom start end
       methylated_count total_count``. Used when the caller has already
       joined per-sample records into one table.

    2. **Sample sheet → Bis-SNP 6+2 BEDs** — columns ``sample_id
       methylation_file [format]``. Each row points at a Bis-SNP
       ``*.6plus2.bed`` (optionally gzipped) file, format defaults to
       ``bissnp_6plus2``. Used for ``--matched-wgbs``.

    3. **Sample sheet → FinaleMe prediction BEDs** — same columns, but
       ``format`` defaults to ``finaleme_bed`` and each file is a
       ``*.prediction.bed.gz`` from the FinaleMe decode step. Used for
       ``--matched-finaleme``.

    Chromosome prefixes are always stripped before returning so downstream
    join keys are uniform.
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
        return _load_sample_sheet(df, Path(path).parent, modality)

    raise InvalidCalibrationError(
        f"{path}: table has neither legacy columns {sorted(legacy_required)} "
        f"nor sample sheet columns {sorted(sheet_required)}; got {sorted(cols)}"
    )


def _load_sample_sheet(
    df: pd.DataFrame,
    base_dir: Path,
    modality: str,
) -> pd.DataFrame:
    """Iterate a sample sheet and concatenate per-sample parsed tables."""
    parts: list[pd.DataFrame] = []
    for _, row in df.iterrows():
        sample_id = str(row["sample_id"])
        file_path = Path(str(row["methylation_file"]))
        if not file_path.is_absolute():
            file_path = base_dir / file_path
        if not file_path.exists():
            raise InvalidCalibrationError(
                f"Sample {sample_id}: file not found: {file_path}"
            )

        file_format = None
        if "format" in df.columns:
            raw_fmt = row.get("format")
            if pd.notna(raw_fmt) and str(raw_fmt).strip():
                file_format = str(raw_fmt).strip().lower()
        if file_format is None:
            file_format = _default_format_for_modality(modality)

        if file_format == "bissnp_6plus2":
            parts.append(_parse_bissnp_6plus2(file_path, sample_id))
        elif file_format == "finaleme_bed":
            parts.append(_parse_finaleme_prediction(file_path, sample_id))
        else:
            raise InvalidCalibrationError(
                f"Sample {sample_id}: unknown file format '{file_format}'. "
                "Supported: bissnp_6plus2, finaleme_bed"
            )

    if not parts:
        raise InvalidCalibrationError("Sample sheet had no rows")
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


def _parse_bissnp_6plus2(path: Path, sample_id: str) -> pd.DataFrame:
    """Parse a Bis-SNP 6+2 BED file into the canonical six-column layout.

    Bis-SNP 6+2 format:
        chrom  start  end  name  score  strand  methylation_pct  total_count

    Example lines (note that a ``track name=...`` header line is valid):
        track name=HD_45 ...
        1  10469  10470  .  500  -  50.00  4
        1  10471  10472  .  750  -  75.00  4

    ``methylated_count`` is derived as ``round(methylation_pct/100 * total_count)``.
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
    out = pd.DataFrame(
        {
            "sample_id": sample_id,
            "chrom": df["chrom"],
            "start": df["start"],
            "end": df["end"],
            "methylated_count": df["methy_count_predict"].round().astype(np.int64),
            "total_count": df["total_count_predict"].round().astype(np.int64),
        }
    )
    return out


# ---------------------------------------------------------------------------
# Region annotation: local CpG density from a CpG index (architecture §3.5)
# ---------------------------------------------------------------------------


def _load_cpg_index_by_chrom(cpg_index_path: Path) -> dict[str, np.ndarray]:
    """Return ``{chrom (no prefix) -> sorted np.int64 array of CpG starts}``.

    Accepts the same CpG index BED format used elsewhere in the package:
    tab-separated ``chrom start end`` with optional gzip. Chromosome names
    are stripped of any leading ``chr`` so the lookup works regardless of
    whether the downstream data uses GRCh37 or UCSC conventions.
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


def compute_region_annotation(
    regions: pd.DataFrame,
    cpg_index_path: str | Path,
    window: int = 1000,
) -> pd.DataFrame:
    """Compute per-row CpG density from a CpG index BED.

    For each ``(chrom, start, end)`` row in ``regions``, the density is
    ``count of CpGs in [center - window/2, center + window/2) / window``
    where ``center = (start + end) / 2``. Rows whose chromosome is not in
    the CpG index get density = 0 (they won't join anywhere downstream
    anyway).

    Parameters
    ----------
    regions
        Any DataFrame with ``chrom``, ``start``, ``end`` columns. Any
        leading ``chr`` prefix is stripped before lookup.
    cpg_index_path
        Path to a CpG index BED file (tab-separated, optionally gzipped).
        Same file that FinaleMe ships as ``data/CpG_index.hg19.bed.gz``
        and that ``BetaValueDeconvolution -cpgIndex`` consumes.
    window
        Window size in base pairs centered on each row. Default 1000.

    Returns
    -------
    New DataFrame with columns ``chrom start end cpg_density`` — same row
    order and same row count as the input, with the chromosome prefix
    already stripped.
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
    highs = centers + (window - half)  # upper bound exclusive, matches window size

    density = np.zeros(len(regions), dtype=np.float64)
    # Group by chromosome for a single searchsorted per chrom
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

    return pd.DataFrame(
        {
            "chrom": chroms.to_numpy(),
            "start": starts,
            "end": ends,
            "cpg_density": density,
        }
    )


def make_region_annotation_from_regions_file(
    regions_path: str | Path,
    cpg_index_path: str | Path,
    output_path: str | Path,
    window: int = 1000,
) -> None:
    """Build a ``region_annotation.tsv`` for a BED/TSV of regions.

    The input can be:
      * a BED file with 3+ columns (``chrom start end``), or
      * a TSV with a header that includes the same columns.

    The output TSV has columns ``chrom start end cpg_density`` and is
    directly consumable by ``finaleme-too train-calibration --region-annotation``.
    """
    path = Path(regions_path)
    # Sniff the header — if the first non-comment line has non-numeric
    # columns after col 0, treat as a headered TSV; otherwise treat as
    # a headerless BED.
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
        df, cpg_index_path=cpg_index_path, window=window
    )
    annotation.to_csv(output_path, sep="\t", index=False)


def train_calibration(
    matched_wgbs: str | Path,
    matched_finaleme: str | Path,
    region_annotation: str | Path | None,
    n_bins_candidates: list[int],
    out_params: str | Path,
    out_report: str | Path,
    cpg_index: str | Path | None = None,
    region_annotation_window: int = 1000,
) -> CalibrationParams:
    """Train per-bin calibration parameters from matched WGBS / FinaleMe samples.

    Each of ``matched_wgbs`` / ``matched_finaleme`` can be:
      * a **legacy joined TSV** with columns ``sample_id chrom start end
        methylated_count total_count`` (pre-joined by the caller), or
      * a **sample sheet** with columns ``sample_id methylation_file [format]``
        pointing at one file per sample — Bis-SNP 6+2 BEDs for the WGBS
        side (``format=bissnp_6plus2``, the default for this modality) or
        FinaleMe ``prediction.bed.gz`` files for the FinaleMe side
        (``format=finaleme_bed``, the default for this modality).

    Chromosome-name prefixes (``chr``) are normalized across both tables
    and the optional region annotation so that GRCh37-style WGBS data
    (``1``, ``2``) joins cleanly with FinaleMe output (``chr1``, ``chr2``).

    CpG density source (priority order):
      1. ``region_annotation`` TSV if explicitly provided
      2. ``cpg_index`` → density is auto-computed via
         :func:`compute_region_annotation` over a ``region_annotation_window``-bp
         window centered on each (chrom, start, end) row
      3. Neither → all rows get density 0 and the calibration collapses
         to a single bin

    Workflow (math doc §6.4):
        1. Load + normalize WGBS and FinaleMe tables
        2. Join on (sample_id, chrom, start, end)
        3. Compute beta = methylated / total on each side
        4. Resolve CpG density per the priority above
        5. tune_n_bins via leave-one-sample-out CV
        6. Fit final calibration and write JSON params + report
    """
    from finaleme_too.preprocessing.calibration_eval import fit_calibration, tune_n_bins
    from finaleme_too.io.output_writer import write_calibration_report

    wgbs_df = _load_matched_table(matched_wgbs, modality="wgbs")
    fme_df = _load_matched_table(matched_finaleme, modality="finaleme")

    join_keys = ["sample_id", "chrom", "start", "end"]
    merged = wgbs_df.merge(
        fme_df, on=join_keys, suffixes=("_wgbs", "_fme")
    )
    if merged.empty:
        raise InvalidCalibrationError(
            "No overlapping (sample_id, chrom, start, end) between WGBS and "
            "FinaleMe tables. Check sample IDs, genomic coordinates, and "
            "chromosome naming — the loader already strips any 'chr' prefix."
        )

    # CpG density: prefer an explicit --region-annotation file; otherwise
    # auto-compute from --cpg-index so the user does not have to pre-build
    # the annotation themselves.
    if region_annotation is not None:
        ann = pd.read_csv(region_annotation, sep="\t", comment="#")
        ann = _normalize_chrom(ann)
        merged = merged.merge(
            ann[["chrom", "start", "end", "cpg_density"]],
            on=["chrom", "start", "end"],
            how="left",
        )
    elif cpg_index is not None:
        # Compute density on the fly for the distinct (chrom, start, end)
        # rows we actually need, rather than building a genome-wide file.
        unique_regions = merged[["chrom", "start", "end"]].drop_duplicates()
        ann = compute_region_annotation(
            unique_regions,
            cpg_index_path=cpg_index,
            window=region_annotation_window,
        )
        merged = merged.merge(
            ann[["chrom", "start", "end", "cpg_density"]],
            on=["chrom", "start", "end"],
            how="left",
        )
    if "cpg_density" not in merged.columns:
        merged["cpg_density"] = 0.0
    merged["cpg_density"] = merged["cpg_density"].fillna(0.0)

    wgbs_beta = (
        merged["methylated_count_wgbs"] / merged["total_count_wgbs"].clip(lower=1)
    ).to_numpy(dtype=np.float64)
    fme_beta = (
        merged["methylated_count_fme"] / merged["total_count_fme"].clip(lower=1)
    ).to_numpy(dtype=np.float64)
    density = merged["cpg_density"].to_numpy(dtype=np.float64)
    sample_ids = merged["sample_id"].astype(str).to_numpy()

    tune_result = tune_n_bins(
        finaleme_beta=fme_beta,
        wgbs_beta=wgbs_beta,
        cpg_density=density,
        sample_ids=sample_ids,
        n_bins_candidates=n_bins_candidates,
    )
    best_B = int(tune_result["selected_n_bins"])

    final = fit_calibration(
        finaleme_beta=fme_beta,
        wgbs_beta=wgbs_beta,
        cpg_density=density,
        n_bins=best_B,
    )

    params = CalibrationParams(
        n_bins=final.n_bins,
        bin_edges=final.bin_edges,
        a=final.a,
        c=final.c,
        log_dispersion=final.log_dispersion,
        training_metadata={
            "n_training_samples": int(len(np.unique(sample_ids))),
            "n_observations": int(len(merged)),
            "n_bins_candidates": list(n_bins_candidates),
            "tune_result": tune_result,
        },
    )
    params.save(out_params)

    report = {
        "calibration_version": "1.0",
        "n_training_samples": int(len(np.unique(sample_ids))),
        "n_bins": best_B,
        "overall_metrics": final.overall,
        "per_bin_metrics": final.per_bin_metrics,
        "candidates": tune_result["candidates"],
    }
    write_calibration_report(report, out_report)
    return params


__all__ = [
    "CalibrationParams",
    "apply_calibration",
    "compute_region_annotation",
    "load_default_calibration",
    "make_region_annotation_from_regions_file",
    "shipped_default_calibration_path",
    "train_calibration",
]
