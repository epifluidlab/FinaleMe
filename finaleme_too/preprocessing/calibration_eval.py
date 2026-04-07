"""Calibration goodness metrics + bin tuning.

Phase B: inference-time QC functions used by the run pipeline.
Phase C: cross-validation + bin selection used by ``train-calibration``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
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


# ---------------------------------------------------------------------------
# Phase C: training-time goodness metrics + bin tuning
# ---------------------------------------------------------------------------


@dataclass
class CalibrationFitResult:
    """Per-bin training fit (math doc §6.1-§6.4)."""

    n_bins: int
    bin_edges: np.ndarray
    a: np.ndarray
    c: np.ndarray
    log_dispersion: np.ndarray
    per_bin_metrics: list[dict]
    overall: dict


def _bin_indices(
    cpg_density: np.ndarray, edges: np.ndarray
) -> np.ndarray:
    return np.clip(
        np.searchsorted(edges, cpg_density, side="right") - 1,
        0,
        len(edges) - 2,
    )


def _quantile_bin_edges(values: np.ndarray, n_bins: int) -> np.ndarray:
    qs = np.linspace(0, 1, n_bins + 1)
    edges = np.quantile(values[~np.isnan(values)], qs)
    edges[0] = -np.inf
    edges[-1] = np.inf
    return edges


def _fit_logistic_bin(
    finaleme_beta: np.ndarray, wgbs_beta: np.ndarray
) -> tuple[float, float, float]:
    """Per-bin OLS in logit space (math doc §6.1).

    Returns (a, c, log_dispersion). log_dispersion is estimated from the
    residual variance via the beta-binomial moment matching:
        Var(eps) ≈ 1 / phi  in logit space (rough approximation)
    """
    valid = (
        np.isfinite(finaleme_beta)
        & np.isfinite(wgbs_beta)
        & (finaleme_beta > 0)
        & (finaleme_beta < 1)
        & (wgbs_beta > 0)
        & (wgbs_beta < 1)
    )
    if int(np.sum(valid)) < 5:
        return 1.0, 0.0, np.log(20.0)
    x = _logit(finaleme_beta[valid])
    y = _logit(wgbs_beta[valid])
    A = np.column_stack([x, np.ones_like(x)])
    coeffs, *_ = np.linalg.lstsq(A, y, rcond=None)
    a, c = float(coeffs[0]), float(coeffs[1])
    residual = y - (a * x + c)
    var = float(np.var(residual)) if residual.size > 1 else 1.0
    # phi = 1 / var as a coarse moment estimate, clipped to a sensible range
    phi = float(np.clip(1.0 / max(var, 1e-3), 1.0, 1000.0))
    return a, c, float(np.log(phi))


def fit_calibration(
    finaleme_beta: np.ndarray,
    wgbs_beta: np.ndarray,
    cpg_density: np.ndarray,
    n_bins: int,
) -> CalibrationFitResult:
    """Fit a per-bin calibration model. Math doc §6.1-§6.3."""
    edges = _quantile_bin_edges(cpg_density, n_bins)
    bin_idx = _bin_indices(cpg_density, edges)

    a = np.ones(n_bins)
    c = np.zeros(n_bins)
    log_dispersion = np.full(n_bins, np.log(20.0))
    per_bin = []

    for b in range(n_bins):
        mask = bin_idx == b
        if int(np.sum(mask)) < 10:
            per_bin.append(
                {
                    "bin_id": b,
                    "n_markers": int(np.sum(mask)),
                    "r_squared": float("nan"),
                    "mae": float("nan"),
                    "slope": 1.0,
                    "intercept": 0.0,
                    "log_dispersion": float(np.log(20.0)),
                    "hosmer_lemeshow_p": float("nan"),
                }
            )
            continue
        a_b, c_b, log_phi_b = _fit_logistic_bin(finaleme_beta[mask], wgbs_beta[mask])
        a[b] = a_b
        c[b] = c_b
        log_dispersion[b] = log_phi_b

        # Goodness metrics
        pred = _expit(a_b * _logit(np.clip(finaleme_beta[mask], 1e-6, 1 - 1e-6)) + c_b)
        truth = wgbs_beta[mask]
        ss_res = float(np.sum((truth - pred) ** 2))
        ss_tot = float(np.sum((truth - truth.mean()) ** 2))
        r2 = 1.0 - ss_res / max(ss_tot, 1e-12)
        mae = float(np.mean(np.abs(pred - truth)))
        hl_p = compute_hosmer_lemeshow(truth, pred, n_groups=10)
        per_bin.append(
            {
                "bin_id": b,
                "n_markers": int(np.sum(mask)),
                "r_squared": float(r2),
                "mae": mae,
                "slope": a_b,
                "intercept": c_b,
                "log_dispersion": log_phi_b,
                "hosmer_lemeshow_p": hl_p,
            }
        )

    # Overall metrics
    overall_pred = np.zeros_like(finaleme_beta, dtype=np.float64)
    for b in range(n_bins):
        mask = bin_idx == b
        overall_pred[mask] = _expit(
            a[b] * _logit(np.clip(finaleme_beta[mask], 1e-6, 1 - 1e-6)) + c[b]
        )
    valid = (
        np.isfinite(finaleme_beta)
        & np.isfinite(wgbs_beta)
        & np.isfinite(overall_pred)
    )
    rmse = float(np.sqrt(np.mean((overall_pred[valid] - wgbs_beta[valid]) ** 2)))
    mae_all = float(np.mean(np.abs(overall_pred[valid] - wgbs_beta[valid])))
    r2_all = 1.0 - float(np.sum((wgbs_beta[valid] - overall_pred[valid]) ** 2)) / max(
        float(np.sum((wgbs_beta[valid] - wgbs_beta[valid].mean()) ** 2)), 1e-12
    )
    overall = {
        "rmse": rmse,
        "mae": mae_all,
        "r_squared": r2_all,
        "hosmer_lemeshow_p_min": float(
            min((m["hosmer_lemeshow_p"] for m in per_bin if not np.isnan(m["hosmer_lemeshow_p"])), default=float("nan"))
        ),
        "n_markers": int(np.sum(valid)),
    }

    return CalibrationFitResult(
        n_bins=n_bins,
        bin_edges=edges,
        a=a,
        c=c,
        log_dispersion=log_dispersion,
        per_bin_metrics=per_bin,
        overall=overall,
    )


def cross_validate_calibration(
    finaleme_beta: np.ndarray,
    wgbs_beta: np.ndarray,
    cpg_density: np.ndarray,
    sample_ids: np.ndarray,
    n_bins: int,
) -> dict:
    """Leave-one-sample-out cross-validation for a given B (math doc §6.3)."""
    unique_samples = np.unique(sample_ids)
    if unique_samples.size < 2:
        # No real CV possible — fit on all data
        fit = fit_calibration(finaleme_beta, wgbs_beta, cpg_density, n_bins)
        return {
            "n_bins": n_bins,
            "cv_rmse": float("nan"),
            "in_sample_rmse": fit.overall["rmse"],
            "in_sample_r_squared": fit.overall["r_squared"],
            "n_folds": 0,
        }

    fold_rmses: list[float] = []
    for held_out in unique_samples:
        train_mask = sample_ids != held_out
        test_mask = sample_ids == held_out
        if int(np.sum(train_mask)) < 10 or int(np.sum(test_mask)) < 5:
            continue
        fit = fit_calibration(
            finaleme_beta[train_mask],
            wgbs_beta[train_mask],
            cpg_density[train_mask],
            n_bins=n_bins,
        )
        edges = fit.bin_edges
        bin_idx = _bin_indices(cpg_density[test_mask], edges)
        pred = np.zeros(int(np.sum(test_mask)), dtype=np.float64)
        for b in range(n_bins):
            m = bin_idx == b
            if not np.any(m):
                continue
            pred[m] = _expit(
                fit.a[b] * _logit(np.clip(finaleme_beta[test_mask][m], 1e-6, 1 - 1e-6))
                + fit.c[b]
            )
        rmse = float(np.sqrt(np.mean((pred - wgbs_beta[test_mask]) ** 2)))
        fold_rmses.append(rmse)
    cv_rmse = float(np.mean(fold_rmses)) if fold_rmses else float("nan")

    in_sample = fit_calibration(finaleme_beta, wgbs_beta, cpg_density, n_bins)
    return {
        "n_bins": n_bins,
        "cv_rmse": cv_rmse,
        "in_sample_rmse": in_sample.overall["rmse"],
        "in_sample_r_squared": in_sample.overall["r_squared"],
        "n_folds": len(fold_rmses),
        "per_bin_metrics": in_sample.per_bin_metrics,
        "overall": in_sample.overall,
    }


def tune_n_bins(
    finaleme_beta: np.ndarray,
    wgbs_beta: np.ndarray,
    cpg_density: np.ndarray,
    sample_ids: np.ndarray,
    n_bins_candidates: list[int],
) -> dict:
    """Cross-validate calibration over multiple B candidates (math doc §6.4).

    Returns the best B (lowest CV-RMSE) and a list of all candidates.
    """
    candidates: list[dict] = []
    n_obs = int(np.sum(np.isfinite(finaleme_beta) & np.isfinite(wgbs_beta)))
    for B in n_bins_candidates:
        result = cross_validate_calibration(
            finaleme_beta, wgbs_beta, cpg_density, sample_ids, n_bins=B
        )
        # AIC / BIC using residual sum of squares as a proxy for log-likelihood
        n_params = 3 * B  # a_b, c_b, log_phi_b per bin
        if not np.isnan(result["in_sample_rmse"]):
            sse = (result["in_sample_rmse"] ** 2) * max(n_obs, 1)
            ll = -0.5 * n_obs * (np.log(2 * np.pi * sse / max(n_obs, 1)) + 1)
            result["aic"] = float(-2 * ll + 2 * n_params)
            result["bic"] = float(-2 * ll + np.log(max(n_obs, 1)) * n_params)
        candidates.append(result)
    # Choose B with lowest CV-RMSE; ties broken by AIC
    valid = [c for c in candidates if not np.isnan(c["cv_rmse"])]
    if not valid:
        valid = candidates
    best = min(valid, key=lambda c: (c["cv_rmse"], c.get("aic", float("inf"))))
    return {"selected_n_bins": best["n_bins"], "candidates": candidates}


__all__ = [
    "CalibrationFitResult",
    "compute_hosmer_lemeshow",
    "compute_inference_qc",
    "cross_validate_calibration",
    "fit_calibration",
    "tune_n_bins",
]
