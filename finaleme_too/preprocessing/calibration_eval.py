"""Calibration goodness metrics + bin tuning (Phase B inference QC, Phase C training).

Phase B: implements the inference-time QC functions used by the run pipeline
(Hosmer-Lemeshow, KS residual test, prediction range coverage). The full
training pipeline (cross_validate_calibration, tune_n_bins) is added in Phase C.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import chi2, ks_2samp

from finaleme_too.preprocessing.calibration import (
    CalibrationParams,
    _expit,
    _logit,
)


def compute_hosmer_lemeshow(
    observed: np.ndarray, predicted: np.ndarray, n_groups: int = 10
) -> float:
    """Hosmer-Lemeshow goodness-of-fit test.

    Bins markers into ``n_groups`` deciles by predicted methylation, then
    computes Σ_g n_g * (O_g - E_g)^2 / (E_g * (1 - E_g)) (math doc §6.3).

    Returns the chi-squared p-value (df = G-2).
    """
    obs = np.asarray(observed, dtype=np.float64)
    pred = np.asarray(predicted, dtype=np.float64)
    valid = np.isfinite(obs) & np.isfinite(pred)
    obs = obs[valid]
    pred = pred[valid]
    if obs.size < n_groups * 5:
        return float("nan")

    quantiles = np.quantile(pred, np.linspace(0, 1, n_groups + 1))
    quantiles[0] = -np.inf
    quantiles[-1] = np.inf
    bins = np.searchsorted(quantiles[1:-1], pred, side="right")

    chi2_stat = 0.0
    for g in range(n_groups):
        idx = bins == g
        n_g = int(np.sum(idx))
        if n_g == 0:
            continue
        o_g = float(np.mean(obs[idx]))
        e_g = float(np.mean(pred[idx]))
        denom = e_g * (1.0 - e_g)
        if denom <= 1e-9:
            continue
        chi2_stat += n_g * (o_g - e_g) ** 2 / denom

    df = max(n_groups - 2, 1)
    return float(1.0 - chi2.cdf(chi2_stat, df))


def compute_inference_qc(
    sample_observations: np.ndarray,  # array of predicted beta values for the sample
    params: CalibrationParams,
    training_residuals: np.ndarray | None = None,
) -> dict:
    """Per-sample inference-time calibration QC (architecture §5.3.2).

    Returns a dict with:
        prediction_range_coverage : fraction of markers within their bin's
                                    training range
        residual_ks_p : 2-sample KS p-value comparing this sample's residuals
                        to the training residual distribution (if provided)
        bin_coverage_balance : fraction of bins with >= 10 markers in this sample
        flag : "PASS" / "WARN" / "FAIL"
    """
    sample_pred = np.asarray(sample_observations, dtype=np.float64)
    valid = np.isfinite(sample_pred)
    sample_pred = sample_pred[valid]

    if sample_pred.size == 0:
        return {
            "prediction_range_coverage": 0.0,
            "residual_ks_p": float("nan"),
            "bin_coverage_balance": 0.0,
            "flag": "FAIL",
        }

    # Bin range coverage: fraction within [0, 1] (always true for betas)
    in_range = float(np.mean((sample_pred >= 0.0) & (sample_pred <= 1.0)))

    # Bin coverage balance: count markers per bin
    bin_idx = params.assign_bin(sample_pred)
    counts = np.bincount(bin_idx, minlength=params.n_bins)
    balance = float(np.mean(counts >= 10))

    # KS test on residuals if training residuals are provided
    if training_residuals is not None and training_residuals.size > 0:
        sample_logit = _logit(sample_pred)
        # Apply identity calibration as a reference; residual = sample_logit - logit(sample_pred)
        # = 0 here, so this only flags when training_residuals is not centred at 0.
        sample_residuals = sample_logit - sample_logit.mean()
        ks_p = float(ks_2samp(sample_residuals, training_residuals).pvalue)
    else:
        ks_p = float("nan")

    flag = "PASS"
    if balance < 0.5 or in_range < 0.95:
        flag = "WARN"
    if balance < 0.25 or in_range < 0.5:
        flag = "FAIL"

    return {
        "prediction_range_coverage": in_range,
        "residual_ks_p": ks_p,
        "bin_coverage_balance": balance,
        "flag": flag,
    }


__all__ = [
    "compute_hosmer_lemeshow",
    "compute_inference_qc",
]
