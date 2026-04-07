"""FinaleMe calibration: apply (Phase B) + train (Phase C).

Math doc §6.1:
    logit(mu_calibrated) = a_b * logit(mu_FinaleMe) + c_b
    phi_b = exp(d_b)

Per CpG-density bin b. Parameters loaded from a JSON config file (or shipped
default). Phase C implements the training pipeline that produces such files.
"""

from __future__ import annotations

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


def _load_matched_table(path: str | Path) -> pd.DataFrame:
    """Load a matched WGBS or FinaleMe table for training.

    Expected columns: sample_id, chrom, start, end, methylated_count, total_count
    Optional: cpg_density (joined from region_annotation if absent).
    """
    df = pd.read_csv(path, sep="\t", comment="#")
    required = {"sample_id", "chrom", "start", "end", "methylated_count", "total_count"}
    missing = required - set(df.columns)
    if missing:
        raise InvalidCalibrationError(f"Matched table missing columns: {sorted(missing)}")
    return df


def train_calibration(
    matched_wgbs: str | Path,
    matched_finaleme: str | Path,
    region_annotation: str | Path | None,
    n_bins_candidates: list[int],
    out_params: str | Path,
    out_report: str | Path,
) -> CalibrationParams:
    """Train per-bin calibration parameters from matched WGBS / FinaleMe samples.

    Workflow (math doc §6.4):
        1. Join WGBS + FinaleMe per-marker observations on (sample_id, chrom, start, end)
        2. Compute beta = methylated / total for each
        3. Join CpG density from region_annotation (or fall back to bin index 0)
        4. tune_n_bins() over the candidate B values via leave-one-sample-out CV
        5. Fit final calibration on all data with the selected B
        6. Write JSON params + JSON report
    """
    from finaleme_too.preprocessing.calibration_eval import fit_calibration, tune_n_bins
    from finaleme_too.io.output_writer import write_calibration_report

    wgbs_df = _load_matched_table(matched_wgbs)
    fme_df = _load_matched_table(matched_finaleme)

    join_keys = ["sample_id", "chrom", "start", "end"]
    merged = wgbs_df.merge(
        fme_df, on=join_keys, suffixes=("_wgbs", "_fme")
    )
    if merged.empty:
        raise InvalidCalibrationError(
            "No overlapping (sample_id, chrom, start, end) between WGBS and FinaleMe tables"
        )

    # CpG density: join from region_annotation if available
    if region_annotation is not None:
        ann = pd.read_csv(region_annotation, sep="\t", comment="#")
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
    "load_default_calibration",
    "shipped_default_calibration_path",
    "train_calibration",
]
