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
    sample_predicted_beta: np.ndarray,  # raw FinaleMe predictions for the sample
    params: CalibrationParams,
    cpg_density: np.ndarray | None = None,
    training_residuals: np.ndarray | None = None,
    training_pred_range: tuple[float, float] | None = None,
) -> dict:
    """Per-sample inference-time calibration QC (architecture §5.3.2).

    Bins are assigned by **CpG density** (matching how the calibration model
    was trained), not by predicted beta. Prediction range coverage measures
    how many predictions fall inside the training prediction range of the
    bin they were assigned to. KS residuals compare this sample's per-marker
    calibration residuals (calibrated - raw, in logit space) against the
    training residual distribution.

    Returns
    -------
    Dict with:
        prediction_range_coverage : fraction of markers whose predicted beta
            lies within the training prediction range (defaults to [0, 1] if
            ``training_pred_range`` is omitted)
        bin_coverage_balance : fraction of bins with >= 10 markers in this sample
        residual_ks_p : two-sample KS p-value comparing the sample's per-marker
            residuals against ``training_residuals`` (if provided)
        flag : "PASS" / "WARN" / "FAIL"
    """
    sample_pred = np.asarray(sample_predicted_beta, dtype=np.float64)
    valid_mask = np.isfinite(sample_pred)
    sample_pred_valid = sample_pred[valid_mask]

    if sample_pred_valid.size == 0:
        return {
            "prediction_range_coverage": 0.0,
            "residual_ks_p": float("nan"),
            "bin_coverage_balance": 0.0,
            "flag": "FAIL",
        }

    # 1) Prediction range coverage: fraction of markers within the training
    # prediction range. If no explicit range is provided, default to (0, 1)
    # exclusive (a more meaningful check than [0, 1] inclusive).
    if training_pred_range is not None:
        lo, hi = training_pred_range
    else:
        lo, hi = 1e-6, 1.0 - 1e-6
    in_range = float(np.mean((sample_pred_valid > lo) & (sample_pred_valid < hi)))

    # 2) Bin coverage balance: bins are assigned by CpG density, NOT by
    # predicted beta. Markers without a known density are routed by
    # ``CalibrationParams.assign_bin`` to ``fallback_bin`` (deterministic).
    if cpg_density is not None:
        density = np.asarray(cpg_density, dtype=np.float64)
        density_valid = density[valid_mask]
        bin_idx = params.assign_bin(density_valid)
    else:
        bin_idx = np.full(sample_pred_valid.size, params.fallback_bin, dtype=np.int64)
    counts = np.bincount(bin_idx, minlength=params.n_bins)
    balance = float(np.mean(counts >= 10))

    # 3) KS test on per-marker calibration residuals (logit space).
    # Residual = logit(calibrated) - logit(raw). With the params' a/c the
    # calibrated value is logit_cal = a*logit(raw) + c, so the residual is
    # (a-1)*logit(raw) + c.
    ks_p = float("nan")
    if training_residuals is not None and np.asarray(training_residuals).size > 0:
        a_per = params.a[bin_idx]
        c_per = params.c[bin_idx]
        raw_logit = _logit(np.clip(sample_pred_valid, 1e-6, 1.0 - 1e-6))
        sample_residuals = (a_per - 1.0) * raw_logit + c_per
        ks_p = float(ks_2samp(sample_residuals, np.asarray(training_residuals)).pvalue)

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


def _cv_fold_task(
    finaleme_beta: np.ndarray,
    wgbs_beta: np.ndarray,
    cpg_density: np.ndarray,
    train_mask: np.ndarray,
    test_mask: np.ndarray,
    n_bins: int,
) -> float | None:
    """Fit one CV fold and return its test-set RMSE (or None if skipped).

    Factored out so ``joblib.Parallel`` can dispatch it across workers.
    """
    if int(np.sum(train_mask)) < 10 or int(np.sum(test_mask)) < 5:
        return None
    fit = fit_calibration(
        finaleme_beta[train_mask],
        wgbs_beta[train_mask],
        cpg_density[train_mask],
        n_bins=n_bins,
    )
    bin_idx = _bin_indices(cpg_density[test_mask], fit.bin_edges)
    pred = np.zeros(int(np.sum(test_mask)), dtype=np.float64)
    for b in range(n_bins):
        m = bin_idx == b
        if not np.any(m):
            continue
        pred[m] = _expit(
            fit.a[b] * _logit(np.clip(finaleme_beta[test_mask][m], 1e-6, 1 - 1e-6))
            + fit.c[b]
        )
    return float(np.sqrt(np.mean((pred - wgbs_beta[test_mask]) ** 2)))


def cross_validate_calibration(
    finaleme_beta: np.ndarray,
    wgbs_beta: np.ndarray,
    cpg_density: np.ndarray,
    sample_ids: np.ndarray,
    n_bins: int,
    threads: int = 1,
) -> dict:
    """Leave-one-sample-out cross-validation for a given B (math doc §6.3).

    ``threads > 1`` parallelizes the fold loop with joblib (architecture
    §14.2 says calibration training is embarrassingly parallel).
    """
    from finaleme_too.utils.parallel import parallel_map

    unique_samples = np.unique(sample_ids)
    if unique_samples.size < 2:
        fit = fit_calibration(finaleme_beta, wgbs_beta, cpg_density, n_bins)
        return {
            "n_bins": n_bins,
            "cv_rmse": float("nan"),
            "in_sample_rmse": fit.overall["rmse"],
            "in_sample_r_squared": fit.overall["r_squared"],
            "n_folds": 0,
        }

    fold_args = [
        (
            finaleme_beta,
            wgbs_beta,
            cpg_density,
            sample_ids != held_out,
            sample_ids == held_out,
            n_bins,
        )
        for held_out in unique_samples
    ]
    # Use the threading backend: _cv_fold_task is pure numpy / scipy, which
    # releases the GIL, and the shared arrays are huge so we avoid the
    # loky pickling overhead.
    fold_rmses_raw = parallel_map(
        lambda args: _cv_fold_task(*args),
        fold_args,
        n_jobs=max(1, threads),
        backend="threading",
    )
    fold_rmses = [r for r in fold_rmses_raw if r is not None]
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
    threads: int = 1,
) -> dict:
    """Cross-validate calibration over multiple B candidates (math doc §6.4).

    The fold × candidate cross-product is flattened and dispatched to joblib
    so a single pool of workers handles the entire grid. This maximizes CPU
    utilization when there are few candidates but many CV folds (or vice
    versa); nested parallelism would leave cores idle.
    """
    from finaleme_too.utils.parallel import parallel_map

    threads = max(1, int(threads))
    unique_samples = np.unique(sample_ids)
    n_obs = int(np.sum(np.isfinite(finaleme_beta) & np.isfinite(wgbs_beta)))

    # Fast path: if we have < 2 distinct samples the CV can't run — just
    # fit each candidate on all data. Keep this serial, it's cheap.
    if unique_samples.size < 2 or threads <= 1:
        candidates: list[dict] = []
        for B in n_bins_candidates:
            result = cross_validate_calibration(
                finaleme_beta, wgbs_beta, cpg_density, sample_ids, n_bins=B,
                threads=1,
            )
            _annotate_aic_bic(result, B, n_obs)
            candidates.append(result)
        return _select_best_candidate(candidates)

    # Parallel path: build the flat (B, held_out) task grid, run all folds
    # through a single joblib pool, then reduce back to per-candidate
    # results. This avoids nested parallelism.
    fold_jobs: list[tuple[int, object, tuple]] = []
    per_candidate_folds: dict[int, list[object]] = {B: [] for B in n_bins_candidates}
    for B in n_bins_candidates:
        for held_out in unique_samples:
            args = (
                finaleme_beta,
                wgbs_beta,
                cpg_density,
                sample_ids != held_out,
                sample_ids == held_out,
                B,
            )
            per_candidate_folds[B].append(held_out)
            fold_jobs.append((B, held_out, args))

    # Use threading backend for the same reason as cross_validate_calibration:
    # the inner work is pure numpy (releases the GIL) and the shared arrays
    # would be pickled once per task with loky.
    fold_rmses_all = parallel_map(
        lambda job: (job[0], job[1], _cv_fold_task(*job[2])),
        fold_jobs,
        n_jobs=threads,
        backend="threading",
    )
    # Group fold RMSEs by B
    rmses_by_B: dict[int, list[float]] = {B: [] for B in n_bins_candidates}
    for B, _held, rmse in fold_rmses_all:
        if rmse is not None:
            rmses_by_B[B].append(rmse)

    # Also run the in-sample fits in parallel (cheap but non-trivial)
    in_sample_fits = parallel_map(
        lambda B: fit_calibration(
            finaleme_beta, wgbs_beta, cpg_density, n_bins=B
        ),
        list(n_bins_candidates),
        n_jobs=threads,
        backend="threading",
    )

    candidates = []
    for B, in_sample in zip(n_bins_candidates, in_sample_fits):
        rmses = rmses_by_B[B]
        result = {
            "n_bins": B,
            "cv_rmse": float(np.mean(rmses)) if rmses else float("nan"),
            "in_sample_rmse": in_sample.overall["rmse"],
            "in_sample_r_squared": in_sample.overall["r_squared"],
            "n_folds": len(rmses),
            "per_bin_metrics": in_sample.per_bin_metrics,
            "overall": in_sample.overall,
        }
        _annotate_aic_bic(result, B, n_obs)
        candidates.append(result)
    return _select_best_candidate(candidates)


def _annotate_aic_bic(result: dict, B: int, n_obs: int) -> None:
    """Attach AIC / BIC columns to a candidate result in place."""
    n_params = 3 * B  # a_b, c_b, log_phi_b per bin
    if not np.isnan(result["in_sample_rmse"]):
        sse = (result["in_sample_rmse"] ** 2) * max(n_obs, 1)
        ll = -0.5 * n_obs * (np.log(2 * np.pi * sse / max(n_obs, 1)) + 1)
        result["aic"] = float(-2 * ll + 2 * n_params)
        result["bic"] = float(-2 * ll + np.log(max(n_obs, 1)) * n_params)


def _select_best_candidate(candidates: list[dict]) -> dict:
    """Pick the B with lowest CV-RMSE, ties broken by AIC."""
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
