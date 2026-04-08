"""Binarization training, evaluation, and bin tuning (v3).

Replaces ``calibration_eval.py`` for FinaleMe mode. The key differences
from the v2 continuous calibration:

  * Per-bin training learns ``(τ_low_b, τ_high_b, ε_U_b, ε_M_b, usable_b)``
    by grid-searching thresholds that maximize classification accuracy on
    WGBS ground truth (math doc §6.2). No OLS, no logit regression.
  * CV uses **chromosome-blocked K-fold** splitting (math doc §6.5) — the
    natural fallback for few-sample training sets because chromosomes are
    genuine replicates of the marker population. Reusing v2's random
    region K-fold would leak the same sample's batch effects across folds.
  * Candidate selection optimizes
    ``score(B) = cv_accuracy × n_usable_markers / n_total_markers``
    (math doc §6.6) — not v2's ``min(cv_rmse)``.
  * Inference-time QC (``compute_inference_qc``) reports fraction-called,
    bin balance, and state-distribution KL shift — not v2's Hosmer-Lemeshow
    or residual KS.

The parallel fold-dispatch pattern is borrowed verbatim from
``calibration_eval.tune_n_bins``: flatten ``(B, fold_id)`` into a single
joblib task grid so we don't nest parallelism.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from finaleme_too.exceptions import InvalidBinarizationError
from finaleme_too.preprocessing._matched_sample_sheet import (
    DEFAULT_REGION_CLASS_THRESHOLDS,
    REGION_CLASS_ORDER,
    _classify_region,
)
from finaleme_too.preprocessing.binarization import (
    STATE_AMBIGUOUS,
    STATE_EXCLUDED,
    STATE_M,
    STATE_U,
    BinarizationParams,
)

log = logging.getLogger(__name__)


# Grid-search ranges for threshold learning (math doc §6.2).
DEFAULT_TAU_LOW_GRID = tuple(np.round(np.arange(0.05, 0.30 + 1e-9, 0.01), 2))
DEFAULT_TAU_HIGH_GRID = tuple(np.round(np.arange(0.35, 0.70 + 1e-9, 0.01), 2))
# Minimum per-bin called sample counts for a usable bin (arch §6.4).
_MIN_MARKERS_PER_STATE = 20


# ---------------------------------------------------------------------------
# BinarizationFitResult dataclass (training-time output)
# ---------------------------------------------------------------------------


@dataclass
class BinarizationFitResult:
    """Training-time output of ``fit_binarization`` for a single B value."""

    params: BinarizationParams
    per_bin_metrics: list[dict]
    overall: dict


# ---------------------------------------------------------------------------
# Single-bin threshold + error-rate fit (math doc §6.2, §6.3)
# ---------------------------------------------------------------------------


def fit_binarization_bin(
    predicted: np.ndarray,
    truth_beta: np.ndarray,
    tau_low_grid: tuple[float, ...] = DEFAULT_TAU_LOW_GRID,
    tau_high_grid: tuple[float, ...] = DEFAULT_TAU_HIGH_GRID,
    max_error_rate: float = 0.15,
    min_markers_per_state: int = _MIN_MARKERS_PER_STATE,
) -> dict:
    """Learn optimal thresholds + error rates for a single context bin.

    Implements math doc §6.2 / §6.3:

        (τ_low*, τ_high*) = argmax_{τ_l, τ_h} accuracy(τ_l, τ_h)
        ε_U = 1 - P(true U | called U)
        ε_M = 1 - P(true M | called M)

    Returns a dict with:
        tau_low, tau_high, eps_U, eps_M, accuracy,
        n_markers_U, n_markers_M, fraction_called, usable

    The per-bin ``usable`` flag (arch §6.4) is True iff
    ``eps_U < max_error_rate`` AND ``eps_M < max_error_rate`` AND
    ``n_markers_U ≥ min_markers_per_state`` AND
    ``n_markers_M ≥ min_markers_per_state``.

    When the bin has < ``min_markers_per_state`` total or no valid
    ground-truth states, returns a fallback dict marked ``usable=False``
    with placeholder thresholds.
    """
    predicted = np.asarray(predicted, dtype=np.float64)
    truth_beta = np.asarray(truth_beta, dtype=np.float64)
    valid = np.isfinite(predicted) & np.isfinite(truth_beta)
    predicted = predicted[valid]
    truth_beta = truth_beta[valid]

    # Binarize ground truth (WGBS) at the canonical 0.2 / 0.8 thresholds
    # (math doc §6.2 / §2B.2).
    true_U = truth_beta < 0.2
    true_M = truth_beta > 0.8

    if predicted.size < min_markers_per_state or (int(np.sum(true_U)) == 0 and int(np.sum(true_M)) == 0):
        return {
            "tau_low": 0.2,
            "tau_high": 0.8,
            "eps_U": 0.5,
            "eps_M": 0.5,
            "accuracy": float("nan"),
            "n_markers_U": 0,
            "n_markers_M": 0,
            "n_total": int(predicted.size),
            "fraction_called": 0.0,
            "usable": False,
        }

    # Grid-search the (τ_low, τ_high) pair maximizing classification accuracy
    # on non-ambiguous calls (math doc §6.2). We require
    # ``tau_low < tau_high`` so Ambiguous is a real category.
    best_acc = -1.0
    best_tau_l = float(tau_low_grid[0])
    best_tau_h = float(tau_high_grid[-1])
    for tau_l in tau_low_grid:
        for tau_h in tau_high_grid:
            if tau_l >= tau_h:
                continue
            called_U = predicted < tau_l
            called_M = predicted > tau_h
            n_called = int(np.sum(called_U) + np.sum(called_M))
            if n_called == 0:
                continue
            correct = int(np.sum(called_U & true_U) + np.sum(called_M & true_M))
            acc = correct / n_called
            if acc > best_acc:
                best_acc = acc
                best_tau_l = float(tau_l)
                best_tau_h = float(tau_h)

    # At optimal thresholds, compute error rates and marker counts.
    called_U = predicted < best_tau_l
    called_M = predicted > best_tau_h
    n_U = int(np.sum(called_U))
    n_M = int(np.sum(called_M))
    n_called = n_U + n_M

    # ε_U = P(true M | called U) = fraction of called-U markers where truth is M
    # ε_M = P(true U | called M) = fraction of called-M markers where truth is U
    if n_U > 0:
        eps_U = float(np.sum(true_M[called_U]) / n_U)
    else:
        eps_U = 0.5
    if n_M > 0:
        eps_M = float(np.sum(true_U[called_M]) / n_M)
    else:
        eps_M = 0.5

    fraction_called = n_called / predicted.size if predicted.size > 0 else 0.0

    usable = (
        eps_U < max_error_rate
        and eps_M < max_error_rate
        and n_U >= min_markers_per_state
        and n_M >= min_markers_per_state
    )

    return {
        "tau_low": best_tau_l,
        "tau_high": best_tau_h,
        "eps_U": eps_U,
        "eps_M": eps_M,
        "accuracy": float(best_acc),
        "n_markers_U": n_U,
        "n_markers_M": n_M,
        "n_total": int(predicted.size),
        "fraction_called": float(fraction_called),
        "usable": bool(usable),
    }


# ---------------------------------------------------------------------------
# Multi-bin fit (fit_binarization)
# ---------------------------------------------------------------------------


def _build_density_edges_per_class(
    density: np.ndarray,
    region_class: np.ndarray,
    class_order: list[str],
    density_sub_bins_per_class: int,
) -> np.ndarray:
    """Compute per-class quantile edges for density sub-bins.

    For each region class, split its density distribution into
    ``density_sub_bins_per_class`` quantile bins. Outer edges are fixed at
    ±inf. Classes with < ``density_sub_bins_per_class`` distinct values get
    a degenerate edge array (still valid for np.searchsorted).

    Shape: ``(len(class_order), density_sub_bins_per_class + 1)``.
    """
    n_classes = len(class_order)
    edges = np.zeros((n_classes, density_sub_bins_per_class + 1), dtype=np.float64)
    edges[:, 0] = -np.inf
    edges[:, -1] = np.inf
    if density_sub_bins_per_class == 1:
        return edges

    quantiles = np.linspace(0, 1, density_sub_bins_per_class + 1)[1:-1]
    for c_idx, cls in enumerate(class_order):
        mask = region_class == cls
        values = density[mask]
        values = values[np.isfinite(values)]
        if values.size < density_sub_bins_per_class:
            # Not enough data to split — use the class's midpoint heuristic
            # (the entire class lands in sub-bin 0).
            edges[c_idx, 1:-1] = 1.0
            continue
        quant_edges = np.quantile(values, quantiles)
        edges[c_idx, 1:-1] = quant_edges
    return edges


def fit_binarization(
    predicted: np.ndarray,
    truth_beta: np.ndarray,
    cpg_density: np.ndarray,
    region_class: np.ndarray | None = None,
    density_sub_bins_per_class: int = 2,
    region_class_thresholds: dict[str, float] | None = None,
    max_error_rate: float = 0.15,
    min_markers_per_state: int = _MIN_MARKERS_PER_STATE,
) -> BinarizationFitResult:
    """Fit a full binarization model across all (class × density sub-bin) bins.

    Parameters
    ----------
    predicted
        FinaleMe predicted_beta values, shape (M,).
    truth_beta
        WGBS ground-truth beta values, shape (M,).
    cpg_density
        Per-marker CpG density, shape (M,). Used to derive region class
        (when ``region_class=None``) and to build per-class density sub-bins.
    region_class
        Optional precomputed per-marker region class ('CGI', 'shore',
        'shelf', 'open_sea'). When None, derived from ``cpg_density`` via
        ``_classify_region`` with ``region_class_thresholds``.
    density_sub_bins_per_class
        Number of density quantile sub-bins within each region class. With
        the default 4 classes and ``density_sub_bins_per_class=2`` this
        yields B=8 total bins (architecture default).
    region_class_thresholds
        Optional override for the default CpG-density thresholds used when
        deriving region class from density.
    max_error_rate
        Bins with ``eps_U`` or ``eps_M`` ≥ this are marked unusable
        (architecture §6.4). Default 0.15.
    min_markers_per_state
        Bins with fewer than this many called-U or called-M markers are
        marked unusable.
    """
    if region_class_thresholds is None:
        region_class_thresholds = dict(DEFAULT_REGION_CLASS_THRESHOLDS)

    cpg_density = np.asarray(cpg_density, dtype=np.float64)
    if region_class is None:
        region_class = _classify_region(cpg_density, region_class_thresholds)
    else:
        region_class = np.asarray(region_class, dtype=object)

    class_order = list(REGION_CLASS_ORDER)
    n_classes = len(class_order)
    n_bins = n_classes * density_sub_bins_per_class

    # Build per-class density sub-bin edges from the training data quantiles.
    density_edges = _build_density_edges_per_class(
        cpg_density, region_class, class_order, density_sub_bins_per_class
    )

    # Assign a bin index to every marker. We build a scratch BinarizationParams
    # with placeholder thresholds just for assign_bin, then fill in the real
    # thresholds after fitting each bin.
    scratch = BinarizationParams(
        n_bins=n_bins,
        n_region_classes=n_classes,
        density_sub_bins_per_class=density_sub_bins_per_class,
        region_class_order=class_order,
        region_class_thresholds=region_class_thresholds,
        density_edges=density_edges,
        tau_low=np.full(n_bins, 0.2, dtype=np.float64),
        tau_high=np.full(n_bins, 0.8, dtype=np.float64),
        eps_U=np.full(n_bins, 0.5, dtype=np.float64),
        eps_M=np.full(n_bins, 0.5, dtype=np.float64),
        usable=np.ones(n_bins, dtype=bool),
        n_markers_U=np.zeros(n_bins, dtype=np.int64),
        n_markers_M=np.zeros(n_bins, dtype=np.int64),
        train_fraction_U=np.full(n_bins, 0.5, dtype=np.float64),
        train_fraction_M=np.full(n_bins, 0.5, dtype=np.float64),
    )
    bin_assignments = scratch.assign_bin(cpg_density, region_class)

    # Fit each bin independently
    tau_low = np.full(n_bins, 0.2, dtype=np.float64)
    tau_high = np.full(n_bins, 0.8, dtype=np.float64)
    eps_U = np.full(n_bins, 0.5, dtype=np.float64)
    eps_M = np.full(n_bins, 0.5, dtype=np.float64)
    usable = np.zeros(n_bins, dtype=bool)
    n_markers_U = np.zeros(n_bins, dtype=np.int64)
    n_markers_M = np.zeros(n_bins, dtype=np.int64)
    per_bin_metrics: list[dict] = []

    for b in range(n_bins):
        mask = bin_assignments == b
        fit = fit_binarization_bin(
            predicted=predicted[mask],
            truth_beta=truth_beta[mask],
            max_error_rate=max_error_rate,
            min_markers_per_state=min_markers_per_state,
        )
        tau_low[b] = fit["tau_low"]
        tau_high[b] = fit["tau_high"]
        eps_U[b] = fit["eps_U"]
        eps_M[b] = fit["eps_M"]
        usable[b] = fit["usable"]
        n_markers_U[b] = fit["n_markers_U"]
        n_markers_M[b] = fit["n_markers_M"]
        fit["bin_id"] = b
        # Also record the region class of this bin for debugging
        class_idx_for_bin = b // density_sub_bins_per_class
        fit["region_class"] = class_order[class_idx_for_bin]
        fit["density_sub_bin"] = b % density_sub_bins_per_class
        per_bin_metrics.append(fit)

    # Per-bin training fraction of U vs M calls (used for inference-time KL)
    total_called = n_markers_U + n_markers_M
    with np.errstate(divide="ignore", invalid="ignore"):
        train_fraction_U = np.where(
            total_called > 0, n_markers_U / np.maximum(total_called, 1), 0.5
        )
        train_fraction_M = np.where(
            total_called > 0, n_markers_M / np.maximum(total_called, 1), 0.5
        )

    params = BinarizationParams(
        n_bins=n_bins,
        n_region_classes=n_classes,
        density_sub_bins_per_class=density_sub_bins_per_class,
        region_class_order=class_order,
        region_class_thresholds=region_class_thresholds,
        density_edges=density_edges,
        tau_low=tau_low,
        tau_high=tau_high,
        eps_U=eps_U,
        eps_M=eps_M,
        usable=usable,
        n_markers_U=n_markers_U,
        n_markers_M=n_markers_M,
        train_fraction_U=train_fraction_U,
        train_fraction_M=train_fraction_M,
    )

    # Overall metrics: accuracy-weighted by called marker count
    total_called_all = int(np.sum(n_markers_U + n_markers_M))
    if total_called_all > 0:
        accuracies = np.array([m["accuracy"] for m in per_bin_metrics])
        counts = np.array([m["n_markers_U"] + m["n_markers_M"] for m in per_bin_metrics])
        finite_mask = np.isfinite(accuracies)
        if np.sum(counts[finite_mask]) > 0:
            overall_accuracy = float(
                np.sum(accuracies[finite_mask] * counts[finite_mask])
                / np.sum(counts[finite_mask])
            )
        else:
            overall_accuracy = float("nan")
    else:
        overall_accuracy = float("nan")

    n_usable_markers = int(
        np.sum([m["n_markers_U"] + m["n_markers_M"] for m, u in zip(per_bin_metrics, usable) if u])
    )
    overall = {
        "accuracy": overall_accuracy,
        "n_usable_bins": int(np.sum(usable)),
        "n_total_bins": int(n_bins),
        "n_usable_markers": n_usable_markers,
        "n_total_markers": int(predicted.size),
        "fraction_usable": float(n_usable_markers / max(predicted.size, 1)),
    }

    return BinarizationFitResult(
        params=params,
        per_bin_metrics=per_bin_metrics,
        overall=overall,
    )


# ---------------------------------------------------------------------------
# Chromosome-blocked K-fold CV (math doc §6.5)
# ---------------------------------------------------------------------------


def _chromosome_fold_masks(
    chrom: np.ndarray,
    n_folds: int,
    seed: int | None = None,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Partition markers into ``n_folds`` folds by chromosome group.

    Each chromosome is assigned to exactly one fold. ``seed`` controls the
    shuffle when there are more chromosomes than folds (so chromosome order
    does not leak into the split).

    Returns a list of (train_mask, test_mask) pairs of length ≤ n_folds.
    Folds with zero markers are dropped.
    """
    chrom = np.asarray(chrom, dtype=object)
    unique = np.unique(chrom)
    if unique.size == 0:
        return []
    if n_folds < 2:
        raise InvalidBinarizationError(f"cv_n_folds must be >= 2, got {n_folds}")

    rng = np.random.default_rng(seed)
    order = np.arange(unique.size)
    rng.shuffle(order)
    shuffled_chroms = unique[order]

    # Assign each chromosome to a fold in round-robin order.
    fold_assign: dict[str, int] = {}
    for i, c in enumerate(shuffled_chroms):
        fold_assign[str(c)] = i % n_folds

    # Per-fold masks
    masks: list[tuple[np.ndarray, np.ndarray]] = []
    for f in range(n_folds):
        test_mask = np.array(
            [fold_assign.get(str(c), -1) == f for c in chrom], dtype=bool
        )
        if int(np.sum(test_mask)) == 0:
            continue
        train_mask = ~test_mask
        if int(np.sum(train_mask)) == 0:
            continue
        masks.append((train_mask, test_mask))
    return masks


def _cv_fold_task(
    predicted: np.ndarray,
    truth_beta: np.ndarray,
    cpg_density: np.ndarray,
    region_class: np.ndarray | None,
    train_mask: np.ndarray,
    test_mask: np.ndarray,
    density_sub_bins_per_class: int,
    region_class_thresholds: dict[str, float],
    max_error_rate: float,
    min_markers_per_state: int,
) -> dict | None:
    """Fit one CV fold on the train mask, evaluate on the test mask.

    Returns a per-fold dict of metrics:
        accuracy  — test-set classification accuracy on called markers
        n_called  — number of test markers receiving a U or M call
        n_usable_markers — training-set count (for score(B) computation)
        n_total_markers  — training-set count
    """
    n_train = int(np.sum(train_mask))
    n_test = int(np.sum(test_mask))
    if n_train < min_markers_per_state or n_test < 5:
        return None

    train_rc = (
        None if region_class is None else np.asarray(region_class)[train_mask]
    )
    fit_result = fit_binarization(
        predicted=predicted[train_mask],
        truth_beta=truth_beta[train_mask],
        cpg_density=cpg_density[train_mask],
        region_class=train_rc,
        density_sub_bins_per_class=density_sub_bins_per_class,
        region_class_thresholds=region_class_thresholds,
        max_error_rate=max_error_rate,
        min_markers_per_state=min_markers_per_state,
    )
    params = fit_result.params

    # Classify the test set using the trained params
    test_rc = (
        None if region_class is None else np.asarray(region_class)[test_mask]
    )
    test_density = cpg_density[test_mask]
    test_pred = predicted[test_mask]
    test_truth = truth_beta[test_mask]

    bin_assignments = params.assign_bin(test_density, test_rc)
    tau_low_per = params.tau_low[bin_assignments]
    tau_high_per = params.tau_high[bin_assignments]
    usable_per = params.usable[bin_assignments]

    finite = np.isfinite(test_pred) & np.isfinite(test_truth)
    valid_mask = usable_per & finite
    if int(np.sum(valid_mask)) == 0:
        return {
            "accuracy": float("nan"),
            "n_called": 0,
            "n_usable_markers": fit_result.overall["n_usable_markers"],
            "n_total_markers": fit_result.overall["n_total_markers"],
        }

    pred_valid = test_pred[valid_mask]
    truth_valid = test_truth[valid_mask]
    tau_l = tau_low_per[valid_mask]
    tau_h = tau_high_per[valid_mask]

    called_U = pred_valid < tau_l
    called_M = pred_valid > tau_h
    truth_U = truth_valid < 0.2
    truth_M = truth_valid > 0.8

    n_called = int(np.sum(called_U) + np.sum(called_M))
    if n_called == 0:
        return {
            "accuracy": float("nan"),
            "n_called": 0,
            "n_usable_markers": fit_result.overall["n_usable_markers"],
            "n_total_markers": fit_result.overall["n_total_markers"],
        }
    correct = int(np.sum(called_U & truth_U) + np.sum(called_M & truth_M))
    accuracy = correct / n_called

    return {
        "accuracy": float(accuracy),
        "n_called": int(n_called),
        "n_usable_markers": int(fit_result.overall["n_usable_markers"]),
        "n_total_markers": int(fit_result.overall["n_total_markers"]),
    }


def cross_validate_binarization(
    predicted: np.ndarray,
    truth_beta: np.ndarray,
    cpg_density: np.ndarray,
    chrom: np.ndarray,
    region_class: np.ndarray | None = None,
    density_sub_bins_per_class: int = 2,
    region_class_thresholds: dict[str, float] | None = None,
    n_folds: int = 10,
    seed: int | None = None,
    max_error_rate: float = 0.15,
    min_markers_per_state: int = _MIN_MARKERS_PER_STATE,
    threads: int = 1,
) -> dict:
    """Chromosome-blocked K-fold CV for a single candidate model.

    Returns a dict with:
        cv_accuracy      — mean test-set accuracy across folds
        cv_accuracy_std  — std across folds
        n_folds          — number of folds actually run
        in_sample_accuracy — accuracy from a fit on the full dataset
        n_usable_markers — from the full-dataset fit (for score(B))
        n_total_markers  — from the full-dataset fit
        fold_metrics     — list of per-fold dicts
    """
    from finaleme_too.utils.parallel import parallel_map

    if region_class_thresholds is None:
        region_class_thresholds = dict(DEFAULT_REGION_CLASS_THRESHOLDS)

    # In-sample fit on the full dataset (used for score(B) + as a baseline).
    in_sample = fit_binarization(
        predicted=predicted,
        truth_beta=truth_beta,
        cpg_density=cpg_density,
        region_class=region_class,
        density_sub_bins_per_class=density_sub_bins_per_class,
        region_class_thresholds=region_class_thresholds,
        max_error_rate=max_error_rate,
        min_markers_per_state=min_markers_per_state,
    )

    fold_masks = _chromosome_fold_masks(chrom, n_folds=n_folds, seed=seed)
    if not fold_masks:
        return {
            "cv_accuracy": float("nan"),
            "cv_accuracy_std": float("nan"),
            "n_folds": 0,
            "in_sample_accuracy": in_sample.overall["accuracy"],
            "in_sample_fit": in_sample,
            "n_usable_markers": in_sample.overall["n_usable_markers"],
            "n_total_markers": in_sample.overall["n_total_markers"],
            "fold_metrics": [],
        }

    fold_args = [
        (
            predicted,
            truth_beta,
            cpg_density,
            region_class,
            train_mask,
            test_mask,
            density_sub_bins_per_class,
            region_class_thresholds,
            max_error_rate,
            min_markers_per_state,
        )
        for train_mask, test_mask in fold_masks
    ]
    fold_results_raw = parallel_map(
        lambda args: _cv_fold_task(*args),
        fold_args,
        n_jobs=max(1, int(threads)),
        backend="threading",
    )
    fold_metrics = [r for r in fold_results_raw if r is not None]
    return _aggregate_fold_metrics(
        fold_metrics=fold_metrics,
        in_sample=in_sample,
    )


def _aggregate_fold_metrics(
    fold_metrics: list[dict],
    in_sample: BinarizationFitResult,
) -> dict:
    """Reduce per-fold dicts into a single aggregated CV result."""
    if not fold_metrics:
        return {
            "cv_accuracy": float("nan"),
            "cv_accuracy_std": float("nan"),
            "n_folds": 0,
            "in_sample_accuracy": in_sample.overall["accuracy"],
            "in_sample_fit": in_sample,
            "n_usable_markers": in_sample.overall["n_usable_markers"],
            "n_total_markers": in_sample.overall["n_total_markers"],
            "fold_metrics": [],
        }

    accuracies = np.array([m["accuracy"] for m in fold_metrics])
    finite_mask = np.isfinite(accuracies)
    if np.sum(finite_mask) == 0:
        cv_acc = float("nan")
        cv_std = float("nan")
    else:
        cv_acc = float(np.mean(accuracies[finite_mask]))
        cv_std = (
            float(np.std(accuracies[finite_mask], ddof=1))
            if np.sum(finite_mask) > 1
            else 0.0
        )

    return {
        "cv_accuracy": cv_acc,
        "cv_accuracy_std": cv_std,
        "n_folds": len(fold_metrics),
        "in_sample_accuracy": in_sample.overall["accuracy"],
        "in_sample_fit": in_sample,
        "n_usable_markers": in_sample.overall["n_usable_markers"],
        "n_total_markers": in_sample.overall["n_total_markers"],
        "fold_metrics": fold_metrics,
    }


# ---------------------------------------------------------------------------
# Bin-count tuning (math doc §6.6)
# ---------------------------------------------------------------------------


def tune_n_bins(
    predicted: np.ndarray,
    truth_beta: np.ndarray,
    cpg_density: np.ndarray,
    chrom: np.ndarray,
    n_bins_candidates: list[int],
    region_class: np.ndarray | None = None,
    region_class_thresholds: dict[str, float] | None = None,
    n_folds: int = 10,
    seed: int | None = None,
    max_error_rate: float = 0.15,
    min_markers_per_state: int = _MIN_MARKERS_PER_STATE,
    threads: int = 1,
) -> dict:
    """Select the best total bin count B via chromosome-blocked CV.

    ``n_bins_candidates`` is the list of total bin counts to try. Each B is
    converted to a per-class sub-bin count ``B // n_region_classes`` (so B
    must be divisible by the number of region classes, 4 in v3.0). If B is
    not divisible, it's rounded up to the next multiple of 4 and a warning
    is logged.

    Scoring (math doc §6.6):
        score(B) = cv_accuracy × n_usable_markers / n_total_markers

    Returns a dict with ``selected_n_bins`` and a ``candidates`` list of
    per-B CV summaries.
    """
    n_classes = len(REGION_CLASS_ORDER)

    candidates: list[dict] = []
    for B in n_bins_candidates:
        # Round B up to the nearest multiple of n_classes
        sub_bins = max(1, int(np.ceil(B / n_classes)))
        effective_B = sub_bins * n_classes
        if effective_B != B:
            log.info(
                "tune_n_bins: rounding B=%d up to %d so it's divisible by %d region classes",
                B, effective_B, n_classes,
            )

        cv_result = cross_validate_binarization(
            predicted=predicted,
            truth_beta=truth_beta,
            cpg_density=cpg_density,
            chrom=chrom,
            region_class=region_class,
            density_sub_bins_per_class=sub_bins,
            region_class_thresholds=region_class_thresholds,
            n_folds=n_folds,
            seed=seed,
            max_error_rate=max_error_rate,
            min_markers_per_state=min_markers_per_state,
            threads=threads,
        )
        # score(B) — higher is better
        n_usable = cv_result["n_usable_markers"]
        n_total = max(cv_result["n_total_markers"], 1)
        if np.isfinite(cv_result["cv_accuracy"]):
            score = cv_result["cv_accuracy"] * (n_usable / n_total)
        else:
            score = float("nan")

        candidates.append(
            {
                "n_bins": effective_B,
                "density_sub_bins_per_class": sub_bins,
                "score": score,
                **cv_result,
            }
        )

    return _select_best_binarization_candidate(candidates)


def _select_best_binarization_candidate(candidates: list[dict]) -> dict:
    """Pick the candidate with the highest ``score(B)``.

    Fallback ladder (when all scores are NaN):
      1. Highest ``cv_accuracy`` (ignoring n_usable)
      2. Highest ``in_sample_accuracy``
      3. First candidate

    This mirrors the v2 ``_select_best_candidate`` fallback structure but
    with the v3 scoring rule.
    """
    if not candidates:
        raise InvalidBinarizationError("No candidates passed to _select_best_binarization_candidate")

    valid_score = [c for c in candidates if np.isfinite(c.get("score", float("nan")))]
    if valid_score:
        best = max(valid_score, key=lambda c: c["score"])
    else:
        valid_cv = [
            c for c in candidates if np.isfinite(c.get("cv_accuracy", float("nan")))
        ]
        if valid_cv:
            best = max(valid_cv, key=lambda c: c["cv_accuracy"])
        else:
            valid_is = [
                c
                for c in candidates
                if np.isfinite(c.get("in_sample_accuracy", float("nan")))
            ]
            if valid_is:
                best = max(valid_is, key=lambda c: c["in_sample_accuracy"])
            else:
                best = candidates[0]

    return {"selected_n_bins": int(best["n_bins"]), "candidates": candidates}


# ---------------------------------------------------------------------------
# Inference-time QC (math doc §6.7, arch §5.3.4)
# ---------------------------------------------------------------------------


def compute_inference_qc(
    called_state: np.ndarray,
    context_bin: np.ndarray,
    params: BinarizationParams,
    min_fraction_called: float = 0.5,
    min_bin_balance: float = 0.5,
    min_bin_size: int = 10,
    kl_warn_threshold: float = 0.5,
    kl_fail_threshold: float = 1.5,
) -> dict:
    """Per-sample inference-time binarization QC.

    Reports:

        fraction_called
            Fraction of markers receiving a U or M call (not Ambiguous/
            Excluded). Warns when below ``min_fraction_called``.
        bin_balance
            Fraction of usable bins with ≥ ``min_bin_size`` called markers.
            Warns when below ``min_bin_balance``.
        state_distribution_kl
            KL divergence of the sample's (U, M) state distribution versus
            the training distribution, computed per bin and averaged.
            Large values indicate distributional shift (e.g. a hyper-
            methylated disease sample vs a healthy training set).
        flag
            Overall ``"PASS"`` / ``"WARN"`` / ``"FAIL"``.
    """
    called_state = np.asarray(called_state, dtype=np.int64)
    context_bin = np.asarray(context_bin, dtype=np.int64)
    n = called_state.size
    if n == 0:
        return {
            "fraction_called": 0.0,
            "bin_balance": 0.0,
            "state_distribution_kl": float("nan"),
            "flag": "FAIL",
        }

    # Fraction called
    called_mask = (called_state == STATE_U) | (called_state == STATE_M)
    fraction_called = float(np.sum(called_mask)) / n

    # Per-bin called counts and U-fractions
    usable_bins = np.where(params.usable)[0]
    if usable_bins.size == 0:
        bin_balance = 0.0
        kl = float("nan")
    else:
        n_balanced = 0
        kls = []
        for b in usable_bins:
            bin_mask = context_bin == b
            if not np.any(bin_mask):
                continue
            states_in_bin = called_state[bin_mask]
            n_u = int(np.sum(states_in_bin == STATE_U))
            n_m = int(np.sum(states_in_bin == STATE_M))
            total = n_u + n_m
            if total >= min_bin_size:
                n_balanced += 1
            if total > 0:
                p_u = n_u / total
                p_m = n_m / total
                q_u = float(params.train_fraction_U[b])
                q_m = float(params.train_fraction_M[b])
                # KL(P || Q) with safe clipping
                eps = 1e-9
                p_u = max(p_u, eps)
                p_m = max(p_m, eps)
                q_u = max(q_u, eps)
                q_m = max(q_m, eps)
                kl_bin = p_u * np.log(p_u / q_u) + p_m * np.log(p_m / q_m)
                kls.append(kl_bin)
        bin_balance = float(n_balanced) / usable_bins.size if usable_bins.size else 0.0
        kl = float(np.mean(kls)) if kls else float("nan")

    # Decide the flag
    flag = "PASS"
    if fraction_called < min_fraction_called or bin_balance < min_bin_balance:
        flag = "WARN"
    if np.isfinite(kl) and kl >= kl_warn_threshold:
        flag = "WARN"
    if (
        fraction_called < min_fraction_called / 2.0
        or bin_balance < min_bin_balance / 2.0
    ):
        flag = "FAIL"
    if np.isfinite(kl) and kl >= kl_fail_threshold:
        flag = "FAIL"

    return {
        "fraction_called": fraction_called,
        "bin_balance": bin_balance,
        "state_distribution_kl": kl,
        "flag": flag,
    }


__all__ = [
    "BinarizationFitResult",
    "compute_inference_qc",
    "cross_validate_binarization",
    "fit_binarization",
    "fit_binarization_bin",
    "tune_n_bins",
]
