"""Calibration goodness metrics + bin tuning.

Phase B: inference-time QC functions used by the run pipeline.
Phase C: cross-validation + bin selection used by ``train-calibration``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import chi2, ks_2samp, pearsonr, spearmanr

from finaleme_too.preprocessing.calibration import (
    CalibrationParams,
    _expit,
    _logit,
)


def _pearson_spearman(truth: np.ndarray, pred: np.ndarray) -> tuple[float, float]:
    """Return ``(pearson_r, spearman_r)`` on the paired finite entries.

    Returns ``(nan, nan)`` when fewer than 2 finite pairs remain or when
    either side is constant (correlation is undefined). Uses ``[0]``
    indexing so the same code works on scipy versions that return a
    plain tuple as well as the newer ``PearsonRResult`` / ``SpearmanrResult``
    namedtuples.
    """
    truth = np.asarray(truth, dtype=np.float64)
    pred = np.asarray(pred, dtype=np.float64)
    mask = np.isfinite(truth) & np.isfinite(pred)
    t = truth[mask]
    p = pred[mask]
    if t.size < 2 or float(np.std(t)) == 0.0 or float(np.std(p)) == 0.0:
        return float("nan"), float("nan")
    try:
        pr_r = float(pearsonr(t, p)[0])
    except Exception:
        pr_r = float("nan")
    try:
        sr_r = float(spearmanr(t, p)[0])
    except Exception:
        sr_r = float("nan")
    return pr_r, sr_r


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
    finaleme_beta: np.ndarray,
    wgbs_beta: np.ndarray,
    wgbs_k: np.ndarray | None = None,
    wgbs_n: np.ndarray | None = None,
) -> tuple[float, float, float]:
    """Per-bin OLS in logit space + beta-binomial dispersion MLE (§6.1, §6.2).

    Returns ``(a, c, log_dispersion)``.

    ``a`` and ``c`` are ordinary least squares coefficients of
    ``logit(wgbs_beta)`` on ``logit(finaleme_beta)`` (math doc §6.1).

    ``log_dispersion`` is the log of the beta-binomial MLE of φ on the
    calibration residuals (math doc §6.2), computed with the raw WGBS
    ``(k, n)`` counts against the calibrated mean ``μ̂ = expit(a·logit(β_fme) + c)``.
    When ``wgbs_k``/``wgbs_n`` are not provided (library users who still
    pass betas directly), fall back to the old residual-variance moment
    estimate so the function keeps working — but the training pipeline
    always passes the counts so production calibration uses MLE.
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

    if wgbs_k is not None and wgbs_n is not None:
        from finaleme_too.utils.beta_binomial import estimate_dispersion_mle

        k_bin = np.asarray(wgbs_k, dtype=np.float64)[valid]
        n_bin = np.asarray(wgbs_n, dtype=np.float64)[valid]
        # Calibrated mean on the original beta scale — this is what the
        # dispersion MLE treats as the mu parameter of the beta-binomial.
        mu_bin = _expit(a * x + c)
        try:
            phi = estimate_dispersion_mle(
                k=k_bin, n=n_bin, mu=mu_bin,
                phi_init=50.0, bounds=(1.0, 1000.0),
            )
        except Exception:
            # Extremely degenerate bins (e.g. all n=1) can still fail; fall
            # back to the residual-variance heuristic.
            residual = y - (a * x + c)
            var = float(np.var(residual)) if residual.size > 1 else 1.0
            phi = float(np.clip(1.0 / max(var, 1e-3), 1.0, 1000.0))
    else:
        # Backward-compat path: no raw counts available, moment estimate.
        residual = y - (a * x + c)
        var = float(np.var(residual)) if residual.size > 1 else 1.0
        phi = float(np.clip(1.0 / max(var, 1e-3), 1.0, 1000.0))

    return a, c, float(np.log(max(phi, 1e-6)))


def fit_calibration(
    finaleme_beta: np.ndarray,
    wgbs_beta: np.ndarray,
    cpg_density: np.ndarray,
    n_bins: int,
    wgbs_k: np.ndarray | None = None,
    wgbs_n: np.ndarray | None = None,
    compute_correlations: bool = True,
) -> CalibrationFitResult:
    """Fit a per-bin calibration model. Math doc §6.1-§6.3.

    When ``wgbs_k`` and ``wgbs_n`` (the raw WGBS methylated/total counts)
    are provided, per-bin dispersion is estimated by beta-binomial MLE on
    the calibration residuals (math doc §6.2). When they are omitted, the
    old residual-variance moment estimate is used for backward compatibility.

    ``compute_correlations`` toggles the per-bin and overall Pearson /
    Spearman computation. Set to ``False`` when the fit is a throwaway CV
    training-fold fit — the per-bin correlations aren't used downstream in
    that case (the fold's reported test-set correlations come from
    ``_cv_fold_task``, not from here) and ``spearmanr`` alone was ~25% of
    total runtime. The final ``train_calibration`` call keeps the default
    ``True`` so the report still carries the per-bin correlations.
    """
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
                    "pearson_raw": float("nan"),
                    "pearson_calibrated": float("nan"),
                    "spearman_raw": float("nan"),
                    "spearman_calibrated": float("nan"),
                }
            )
            continue
        k_slice = wgbs_k[mask] if wgbs_k is not None else None
        n_slice = wgbs_n[mask] if wgbs_n is not None else None
        a_b, c_b, log_phi_b = _fit_logistic_bin(
            finaleme_beta[mask],
            wgbs_beta[mask],
            wgbs_k=k_slice,
            wgbs_n=n_slice,
        )
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
        # Per-bin correlations + Hosmer-Lemeshow are only worth computing for
        # the final user-facing fit — spearmanr alone was ~25% of CV fold
        # runtime and the output isn't consumed during CV.
        if compute_correlations:
            raw_in_bin = finaleme_beta[mask]
            hl_p = compute_hosmer_lemeshow(truth, pred, n_groups=10)
            # Per-bin correlations before/after calibration. Note: within a
            # single bin, Spearman is invariant to the transform when
            # a_b > 0 (monotonic maps preserve ranks), so spearman_raw ≈
            # spearman_calibrated per-bin; it's the across-bin Spearman
            # that can change.
            pearson_raw_b, spearman_raw_b = _pearson_spearman(truth, raw_in_bin)
            pearson_cal_b, spearman_cal_b = _pearson_spearman(truth, pred)
        else:
            hl_p = float("nan")
            pearson_raw_b = float("nan")
            pearson_cal_b = float("nan")
            spearman_raw_b = float("nan")
            spearman_cal_b = float("nan")
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
                "pearson_raw": pearson_raw_b,
                "pearson_calibrated": pearson_cal_b,
                "spearman_raw": spearman_raw_b,
                "spearman_calibrated": spearman_cal_b,
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
    # Overall correlations before/after calibration (math doc §6.3 extension):
    #   - ``*_raw``         — Pearson/Spearman between WGBS truth and the
    #                         uncalibrated FinaleMe prediction. Invariant
    #                         to n_bins; this is a property of the input
    #                         data and serves as the "baseline".
    #   - ``*_calibrated``  — same correlation after applying the per-bin
    #                         calibration model. Should be >= the raw
    #                         correlation on the full dataset when the
    #                         calibration is adding predictive value.
    # Pearson captures linear agreement (changes noticeably after the
    # logit-space transform), Spearman captures rank/monotonic agreement
    # (can change across bins because different bins apply different maps).
    if compute_correlations:
        pearson_raw, spearman_raw = _pearson_spearman(
            wgbs_beta[valid], finaleme_beta[valid]
        )
        pearson_cal, spearman_cal = _pearson_spearman(
            wgbs_beta[valid], overall_pred[valid]
        )
    else:
        pearson_raw = float("nan")
        pearson_cal = float("nan")
        spearman_raw = float("nan")
        spearman_cal = float("nan")
    overall = {
        "rmse": rmse,
        "mae": mae_all,
        "r_squared": r2_all,
        "hosmer_lemeshow_p_min": float(
            min((m["hosmer_lemeshow_p"] for m in per_bin if not np.isnan(m["hosmer_lemeshow_p"])), default=float("nan"))
        ),
        "n_markers": int(np.sum(valid)),
        "pearson_raw": pearson_raw,
        "pearson_calibrated": pearson_cal,
        "spearman_raw": spearman_raw,
        "spearman_calibrated": spearman_cal,
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
    wgbs_k: np.ndarray | None = None,
    wgbs_n: np.ndarray | None = None,
) -> dict | None:
    """Fit one CV fold and return per-fold test-set metrics.

    Returns ``None`` when the fold is too small to fit, otherwise a dict
    with:

      * ``rmse`` — test-set RMSE between calibrated prediction and WGBS
      * ``pearson_raw``, ``spearman_raw`` — correlations of WGBS truth vs
        the uncalibrated FinaleMe prediction on the held-out rows
      * ``pearson_calibrated``, ``spearman_calibrated`` — same, after
        applying the per-bin calibration fit on the training rows
      * ``n`` — number of held-out rows contributing to the metrics

    Factored out so ``joblib.Parallel`` can dispatch it across workers.
    """
    if int(np.sum(train_mask)) < 10 or int(np.sum(test_mask)) < 5:
        return None
    # Skip the training-fold's own per-bin/overall correlations — they're
    # discarded. The only test-set correlations that make it into the
    # report are computed below on ``fme_test`` / ``wgbs_test`` / ``pred``.
    fit = fit_calibration(
        finaleme_beta[train_mask],
        wgbs_beta[train_mask],
        cpg_density[train_mask],
        n_bins=n_bins,
        wgbs_k=wgbs_k[train_mask] if wgbs_k is not None else None,
        wgbs_n=wgbs_n[train_mask] if wgbs_n is not None else None,
        compute_correlations=False,
    )
    bin_idx = _bin_indices(cpg_density[test_mask], fit.bin_edges)
    fme_test = finaleme_beta[test_mask]
    wgbs_test = wgbs_beta[test_mask]
    pred = np.zeros(fme_test.size, dtype=np.float64)
    for b in range(n_bins):
        m = bin_idx == b
        if not np.any(m):
            continue
        pred[m] = _expit(
            fit.a[b] * _logit(np.clip(fme_test[m], 1e-6, 1 - 1e-6))
            + fit.c[b]
        )
    rmse = float(np.sqrt(np.mean((pred - wgbs_test) ** 2)))
    pearson_raw, spearman_raw = _pearson_spearman(wgbs_test, fme_test)
    pearson_cal, spearman_cal = _pearson_spearman(wgbs_test, pred)
    return {
        "rmse": rmse,
        "pearson_raw": pearson_raw,
        "pearson_calibrated": pearson_cal,
        "spearman_raw": spearman_raw,
        "spearman_calibrated": spearman_cal,
        "n": int(fme_test.size),
    }


def _resolve_cv_strategy(strategy: str, n_unique_samples: int) -> str:
    """Resolve ``"auto"`` to a concrete CV strategy.

    With ``>= 2`` distinct samples, ``auto`` chooses leave-one-sample-out,
    which is the gold-standard estimator of generalization error to a new
    sample (it accounts for between-sample variability). With only one
    sample, leave-one-sample-out is undefined, so ``auto`` falls back to
    random region K-fold — this works with as few as 1 sample but can
    only estimate generalization to *new regions within the same samples*,
    not to a completely new sample. The trade-off is intentional: a
    biased but finite CV estimate is strictly better than NaN for picking
    the best ``n_bins``.
    """
    s = (strategy or "auto").lower()
    if s == "auto":
        return "sample" if n_unique_samples >= 2 else "region"
    if s not in {"sample", "region", "none"}:
        raise ValueError(
            f"Unknown cv_strategy {strategy!r}; expected one of "
            "'auto', 'sample', 'region', 'none'."
        )
    return s


def _build_region_fold_masks(
    n_rows: int,
    n_folds: int,
    seed: int | None,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Build ``n_folds`` random (train_mask, test_mask) pairs over ``n_rows``.

    The fold assignment is a deterministic permutation of ``range(n_rows)``
    cut into ``n_folds`` contiguous slices. With ``seed`` fixed, the same
    row index always lands in the same fold so the CV estimate is
    reproducible between runs.
    """
    if n_folds < 2:
        raise ValueError(f"cv_n_folds must be >= 2; got {n_folds}")
    n_folds = min(n_folds, n_rows)
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n_rows)
    # np.array_split handles non-divisible splits gracefully.
    test_chunks = np.array_split(perm, n_folds)
    masks: list[tuple[np.ndarray, np.ndarray]] = []
    for test_idx in test_chunks:
        if test_idx.size == 0:
            continue
        test_mask = np.zeros(n_rows, dtype=bool)
        test_mask[test_idx] = True
        masks.append((~test_mask, test_mask))
    return masks


def cross_validate_calibration(
    finaleme_beta: np.ndarray,
    wgbs_beta: np.ndarray,
    cpg_density: np.ndarray,
    sample_ids: np.ndarray,
    n_bins: int,
    threads: int = 1,
    wgbs_k: np.ndarray | None = None,
    wgbs_n: np.ndarray | None = None,
    cv_strategy: str = "auto",
    cv_n_folds: int = 10,
    cv_seed: int | None = None,
) -> dict:
    """Cross-validate a calibration fit for a given ``n_bins`` (math doc §6.3).

    Parameters
    ----------
    cv_strategy : {"auto", "sample", "region", "none"}
        * ``"auto"`` — leave-one-sample-out when ``>=2`` samples are
          available, otherwise random region K-fold. This is the default
          and the recommended setting.
        * ``"sample"`` — always leave-one-sample-out; returns
          ``cv_rmse=NaN`` when there are fewer than 2 samples.
        * ``"region"`` — random K-fold CV on row indices, regardless of
          how many samples are present. Use this when you want a finite
          CV estimate from a 1-sample run (with the caveat that the test
          and train folds share sample-level effects, so the CV error is
          biased downward compared to genuine leave-one-sample-out).
        * ``"none"`` — skip CV, return NaN ``cv_rmse`` plus the in-sample
          metrics; ``tune_n_bins`` will then pick ``n_bins`` by AIC.
    cv_n_folds : int
        Number of folds for ``"region"`` / auto-region mode. Default 10
        (each fold holds out ~10% of the rows).
    cv_seed : int | None
        Seed for the region-level permutation; ``None`` uses fresh
        entropy each call.
    ``threads > 1`` parallelizes the fold loop with joblib (architecture
    §14.2 says calibration training is embarrassingly parallel).

    ``wgbs_k``/``wgbs_n`` enable per-bin beta-binomial dispersion MLE
    (math doc §6.2) when the raw WGBS counts are available.
    """
    from finaleme_too.utils.parallel import parallel_map

    unique_samples = np.unique(sample_ids)
    strategy = _resolve_cv_strategy(cv_strategy, unique_samples.size)

    # Compute the in-sample fit once — it's reported by every branch.
    in_sample = fit_calibration(
        finaleme_beta, wgbs_beta, cpg_density, n_bins,
        wgbs_k=wgbs_k, wgbs_n=wgbs_n,
    )
    in_sample_result = {
        "n_bins": n_bins,
        "in_sample_rmse": in_sample.overall["rmse"],
        "in_sample_r_squared": in_sample.overall["r_squared"],
        "per_bin_metrics": in_sample.per_bin_metrics,
        "overall": in_sample.overall,
        "cv_strategy": strategy,
    }

    # "none": skip CV entirely.
    if strategy == "none":
        return {
            **in_sample_result,
            "cv_rmse": float("nan"),
            "n_folds": 0,
        }

    # "sample": leave-one-sample-out. Returns NaN if only one sample was
    # provided AND the user explicitly asked for sample mode (rather than
    # letting auto fall through to region).
    if strategy == "sample":
        if unique_samples.size < 2:
            return {
                **in_sample_result,
                "cv_rmse": float("nan"),
                "n_folds": 0,
            }
        fold_masks = [
            (sample_ids != held_out, sample_ids == held_out)
            for held_out in unique_samples
        ]
    else:
        # "region": random K-fold on row indices.
        fold_masks = _build_region_fold_masks(
            n_rows=finaleme_beta.size,
            n_folds=cv_n_folds,
            seed=cv_seed,
        )

    fold_args = [
        (
            finaleme_beta,
            wgbs_beta,
            cpg_density,
            train_mask,
            test_mask,
            n_bins,
            wgbs_k,
            wgbs_n,
        )
        for train_mask, test_mask in fold_masks
    ]
    # Use the threading backend: _cv_fold_task is pure numpy / scipy, which
    # releases the GIL, and the shared arrays are huge so we avoid the
    # loky pickling overhead.
    fold_results_raw = parallel_map(
        lambda args: _cv_fold_task(*args),
        fold_args,
        n_jobs=max(1, threads),
        backend="threading",
    )
    fold_results = [r for r in fold_results_raw if r is not None]
    aggregated = _aggregate_fold_metrics(fold_results)

    return {
        **in_sample_result,
        **aggregated,
    }


def _aggregate_fold_metrics(fold_results: list[dict]) -> dict:
    """Aggregate per-fold dicts into CV-level mean metrics.

    Returns a dict with:

      * ``cv_rmse`` — mean of per-fold RMSEs (NaN if no folds ran)
      * ``cv_rmse_std`` — standard deviation across folds (0 for a single
        fold, NaN for zero folds)
      * ``cv_pearson_raw``, ``cv_pearson_calibrated``,
        ``cv_spearman_raw``, ``cv_spearman_calibrated`` — same mean/std
        aggregation for the four test-set correlations
      * ``n_folds`` — number of folds that actually produced metrics
      * ``fold_metrics`` — the raw per-fold list for inspection

    NaN fold values are skipped in the mean so a single pathological fold
    doesn't poison the whole aggregate.
    """
    if not fold_results:
        return {
            "cv_rmse": float("nan"),
            "cv_rmse_std": float("nan"),
            "cv_pearson_raw": float("nan"),
            "cv_pearson_calibrated": float("nan"),
            "cv_spearman_raw": float("nan"),
            "cv_spearman_calibrated": float("nan"),
            "n_folds": 0,
            "fold_metrics": [],
        }

    def _mean(key: str) -> float:
        vals = [r[key] for r in fold_results if np.isfinite(r[key])]
        return float(np.mean(vals)) if vals else float("nan")

    def _std(key: str) -> float:
        vals = [r[key] for r in fold_results if np.isfinite(r[key])]
        if len(vals) < 2:
            return 0.0 if vals else float("nan")
        return float(np.std(vals, ddof=1))

    return {
        "cv_rmse": _mean("rmse"),
        "cv_rmse_std": _std("rmse"),
        "cv_pearson_raw": _mean("pearson_raw"),
        "cv_pearson_calibrated": _mean("pearson_calibrated"),
        "cv_spearman_raw": _mean("spearman_raw"),
        "cv_spearman_calibrated": _mean("spearman_calibrated"),
        "n_folds": len(fold_results),
        "fold_metrics": fold_results,
    }


def tune_n_bins(
    finaleme_beta: np.ndarray,
    wgbs_beta: np.ndarray,
    cpg_density: np.ndarray,
    sample_ids: np.ndarray,
    n_bins_candidates: list[int],
    threads: int = 1,
    wgbs_k: np.ndarray | None = None,
    wgbs_n: np.ndarray | None = None,
    cv_strategy: str = "auto",
    cv_n_folds: int = 10,
    cv_seed: int | None = None,
) -> dict:
    """Cross-validate calibration over multiple B candidates (math doc §6.4).

    The fold × candidate cross-product is flattened and dispatched to joblib
    so a single pool of workers handles the entire grid. This maximizes CPU
    utilization when there are few candidates but many CV folds (or vice
    versa); nested parallelism would leave cores idle.

    ``cv_strategy`` / ``cv_n_folds`` / ``cv_seed`` select between
    leave-one-sample-out and random region K-fold CV — see
    :func:`cross_validate_calibration` for the full semantics. The region
    folds are built **once** up front using ``cv_seed`` and reused across
    every ``n_bins`` candidate so the comparison between candidates is
    like-for-like.

    ``wgbs_k``/``wgbs_n`` enable per-bin beta-binomial dispersion MLE
    (math doc §6.2) when the raw WGBS counts are available. Without them
    the code falls back to the residual-variance moment estimate.
    """
    from finaleme_too.utils.parallel import parallel_map

    threads = max(1, int(threads))
    unique_samples = np.unique(sample_ids)
    n_obs = int(np.sum(np.isfinite(finaleme_beta) & np.isfinite(wgbs_beta)))
    strategy = _resolve_cv_strategy(cv_strategy, unique_samples.size)

    # Build the fold mask list ONCE — used for every n_bins candidate.
    # Region folds are seeded so tune results are reproducible across runs
    # with the same cv_seed.
    if strategy == "none":
        fold_masks: list[tuple[np.ndarray, np.ndarray]] = []
    elif strategy == "sample":
        if unique_samples.size < 2:
            fold_masks = []
        else:
            fold_masks = [
                (sample_ids != held_out, sample_ids == held_out)
                for held_out in unique_samples
            ]
    else:  # "region"
        fold_masks = _build_region_fold_masks(
            n_rows=finaleme_beta.size,
            n_folds=cv_n_folds,
            seed=cv_seed,
        )

    # Fast path: serial execution — keep it simple for single-thread runs.
    if threads <= 1 or not fold_masks:
        candidates: list[dict] = []
        for B in n_bins_candidates:
            # Reuse cross_validate_calibration for the CV-skipped branches
            # (no fold_masks) so "none" / single-sample-sample paths share
            # the same bookkeeping.
            if not fold_masks:
                result = cross_validate_calibration(
                    finaleme_beta, wgbs_beta, cpg_density, sample_ids,
                    n_bins=B, threads=1, wgbs_k=wgbs_k, wgbs_n=wgbs_n,
                    cv_strategy=cv_strategy, cv_n_folds=cv_n_folds,
                    cv_seed=cv_seed,
                )
            else:
                fold_results: list[dict] = []
                for train_mask, test_mask in fold_masks:
                    r = _cv_fold_task(
                        finaleme_beta, wgbs_beta, cpg_density,
                        train_mask, test_mask, B, wgbs_k, wgbs_n,
                    )
                    if r is not None:
                        fold_results.append(r)
                in_sample = fit_calibration(
                    finaleme_beta, wgbs_beta, cpg_density, n_bins=B,
                    wgbs_k=wgbs_k, wgbs_n=wgbs_n,
                )
                aggregated = _aggregate_fold_metrics(fold_results)
                result = {
                    "n_bins": B,
                    "in_sample_rmse": in_sample.overall["rmse"],
                    "in_sample_r_squared": in_sample.overall["r_squared"],
                    "per_bin_metrics": in_sample.per_bin_metrics,
                    "overall": in_sample.overall,
                    "cv_strategy": strategy,
                    **aggregated,
                }
            _annotate_aic_bic(result, B, n_obs)
            candidates.append(result)
        return _select_best_candidate(candidates)

    # Parallel path: build the flat (B, fold_id) task grid, run all folds
    # through a single joblib pool, then reduce back to per-candidate
    # results. This avoids nested parallelism.
    fold_jobs: list[tuple[int, int, tuple]] = []
    for B in n_bins_candidates:
        for fold_id, (train_mask, test_mask) in enumerate(fold_masks):
            args = (
                finaleme_beta,
                wgbs_beta,
                cpg_density,
                train_mask,
                test_mask,
                B,
                wgbs_k,
                wgbs_n,
            )
            fold_jobs.append((B, fold_id, args))

    # Use threading backend for the same reason as cross_validate_calibration:
    # the inner work is pure numpy (releases the GIL) and the shared arrays
    # would be pickled once per task with loky.
    fold_results_all = parallel_map(
        lambda job: (job[0], job[1], _cv_fold_task(*job[2])),
        fold_jobs,
        n_jobs=threads,
        backend="threading",
    )
    # Group fold result dicts by B
    fold_results_by_B: dict[int, list[dict]] = {B: [] for B in n_bins_candidates}
    for B, _fold_id, fold_result in fold_results_all:
        if fold_result is not None:
            fold_results_by_B[B].append(fold_result)

    # Also run the in-sample fits in parallel (cheap but non-trivial)
    in_sample_fits = parallel_map(
        lambda B: fit_calibration(
            finaleme_beta, wgbs_beta, cpg_density, n_bins=B,
            wgbs_k=wgbs_k, wgbs_n=wgbs_n,
        ),
        list(n_bins_candidates),
        n_jobs=threads,
        backend="threading",
    )

    candidates = []
    for B, in_sample in zip(n_bins_candidates, in_sample_fits):
        aggregated = _aggregate_fold_metrics(fold_results_by_B[B])
        result = {
            "n_bins": B,
            "in_sample_rmse": in_sample.overall["rmse"],
            "in_sample_r_squared": in_sample.overall["r_squared"],
            "per_bin_metrics": in_sample.per_bin_metrics,
            "overall": in_sample.overall,
            "cv_strategy": strategy,
            **aggregated,
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
    """Pick the B with lowest CV-RMSE, ties broken by AIC.

    When every candidate has ``cv_rmse=NaN`` (e.g. only one training
    sample so leave-one-sample-out CV can't run), fall back to lowest
    AIC rather than silently picking the first candidate — NaN sort
    keys in Python return the first item unchanged, which gave a
    meaningless "best" choice before. When AIC is also unavailable,
    fall through to lowest in-sample RMSE.
    """
    valid_cv = [c for c in candidates if not np.isnan(c.get("cv_rmse", np.nan))]
    if valid_cv:
        best = min(
            valid_cv, key=lambda c: (c["cv_rmse"], c.get("aic", float("inf")))
        )
    else:
        valid_aic = [
            c for c in candidates if not np.isnan(c.get("aic", float("nan")))
        ]
        if valid_aic:
            best = min(valid_aic, key=lambda c: c["aic"])
        else:
            valid_rmse = [
                c
                for c in candidates
                if not np.isnan(c.get("in_sample_rmse", float("nan")))
            ]
            if valid_rmse:
                best = min(valid_rmse, key=lambda c: c["in_sample_rmse"])
            else:
                best = candidates[0]
    return {"selected_n_bins": best["n_bins"], "candidates": candidates}


__all__ = [
    "CalibrationFitResult",
    "compute_hosmer_lemeshow",
    "compute_inference_qc",
    "cross_validate_calibration",
    "fit_calibration",
    "tune_n_bins",
]
