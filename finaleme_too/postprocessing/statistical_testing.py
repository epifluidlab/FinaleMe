"""Cross-sample statistical testing in ILR space + Wilcoxon legacy fallback.

Math doc §7: ILR per-cell-type regression with optional bootstrap SEs.
Math doc §8: Bayesian posterior probability of difference (Phase E).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from finaleme_too.utils.transforms import ilr_transform


@dataclass
class TestResult:
    cell_type: str
    contrast: str
    test_type: str
    mean_a: float
    mean_b: float
    effect_size: float
    se: float
    statistic: float
    p_value: float
    adjusted_p_value: float = float("nan")
    significant: bool = False
    extra: dict | None = None


def compositional_regression_test(
    proportions: np.ndarray,
    sample_ids: list[str],
    group_labels: list[str | None],
    cell_type_names: list[str],
    contrasts: list[tuple[str, str]],
    covariates: pd.DataFrame | None = None,
    fdr_alpha: float = 0.05,
) -> list[TestResult]:
    """Per-cell-type ILR regression Wald test (math doc §7.2-§7.3).

    proportions: (S, K+1) — last column is the unknown component.
    contrasts: list of (group_a, group_b) pairs.
    """
    # Drop samples with no group label
    valid = [i for i, g in enumerate(group_labels) if g is not None]
    if len(valid) < 3:
        return []
    P = proportions[valid]
    groups = [group_labels[i] for i in valid]

    # ILR transform — sum_w_unknown is included so that the simplex is K+1 parts
    # The result has K coords; coord j is the contrast that most directly
    # involves cell type j (because of the Helmert basis).
    ilr_coords = ilr_transform(P)  # (S, K)

    n_cells = len(cell_type_names)
    results: list[TestResult] = []

    for j in range(n_cells):
        y = ilr_coords[:, j]
        for group_a, group_b in contrasts:
            a_mask = np.array([g == group_a for g in groups])
            b_mask = np.array([g == group_b for g in groups])
            if int(a_mask.sum()) < 2 or int(b_mask.sum()) < 2:
                results.append(_skipped_result(cell_type_names[j], group_a, group_b))
                continue
            ya = y[a_mask]
            yb = y[b_mask]
            mean_a = float(np.mean(ya))
            mean_b = float(np.mean(yb))
            n_a = int(a_mask.sum())
            n_b = int(b_mask.sum())
            var_a = float(np.var(ya, ddof=1))
            var_b = float(np.var(yb, ddof=1))
            se = float(np.sqrt(var_a / n_a + var_b / n_b))
            if se < 1e-12:
                results.append(_skipped_result(cell_type_names[j], group_a, group_b))
                continue
            effect = mean_a - mean_b
            z = effect / se
            from scipy.stats import norm

            p = float(2.0 * (1.0 - norm.cdf(abs(z))))
            results.append(
                TestResult(
                    cell_type=cell_type_names[j],
                    contrast=f"{group_a}_vs_{group_b}",
                    test_type="ilr_regression",
                    mean_a=float(np.mean(P[a_mask, j])),
                    mean_b=float(np.mean(P[b_mask, j])),
                    effect_size=effect,
                    se=se,
                    statistic=z,
                    p_value=p,
                )
            )

    _apply_fdr(results, fdr_alpha)
    return results


def wilcoxon_test(
    proportions: np.ndarray,
    group_labels: list[str | None],
    cell_type_names: list[str],
    contrasts: list[tuple[str, str]],
    fdr_alpha: float = 0.05,
) -> list[TestResult]:
    """Mann-Whitney U test on raw proportions (legacy)."""
    from scipy.stats import mannwhitneyu

    results: list[TestResult] = []
    for j, ct in enumerate(cell_type_names):
        for group_a, group_b in contrasts:
            a_idx = [i for i, g in enumerate(group_labels) if g == group_a]
            b_idx = [i for i, g in enumerate(group_labels) if g == group_b]
            if len(a_idx) < 2 or len(b_idx) < 2:
                results.append(_skipped_result(ct, group_a, group_b, test="wilcoxon"))
                continue
            ya = proportions[a_idx, j]
            yb = proportions[b_idx, j]
            try:
                stat, p = mannwhitneyu(ya, yb, alternative="two-sided")
            except ValueError:
                results.append(_skipped_result(ct, group_a, group_b, test="wilcoxon"))
                continue
            results.append(
                TestResult(
                    cell_type=ct,
                    contrast=f"{group_a}_vs_{group_b}",
                    test_type="wilcoxon",
                    mean_a=float(np.mean(ya)),
                    mean_b=float(np.mean(yb)),
                    effect_size=float(np.mean(ya) - np.mean(yb)),
                    se=float("nan"),
                    statistic=float(stat),
                    p_value=float(p),
                )
            )

    _apply_fdr(results, fdr_alpha)
    return results


def _skipped_result(cell_type: str, group_a: str, group_b: str, test: str = "ilr_regression") -> TestResult:
    return TestResult(
        cell_type=cell_type,
        contrast=f"{group_a}_vs_{group_b}",
        test_type=test,
        mean_a=float("nan"),
        mean_b=float("nan"),
        effect_size=float("nan"),
        se=float("nan"),
        statistic=float("nan"),
        p_value=float("nan"),
    )


def apply_fdr(results: list[TestResult], alpha: float, method: str = "fdr_bh") -> None:
    """Benjamini-Hochberg FDR correction in place across the supplied results.

    Modifies each TestResult's ``adjusted_p_value`` and ``significant`` fields.
    Results with NaN p-values are skipped (their adjusted_p_value stays NaN).

    The ``method`` argument is forwarded to statsmodels.stats.multitest.multipletests
    so callers can pick a different correction (e.g. fdr_by, holm) via config.
    """
    pvals = [r.p_value for r in results if not np.isnan(r.p_value)]
    if not pvals:
        return
    from statsmodels.stats.multitest import multipletests

    method_map = {
        "BH": "fdr_bh",
        "bh": "fdr_bh",
        "fdr_bh": "fdr_bh",
        "BY": "fdr_by",
        "by": "fdr_by",
        "fdr_by": "fdr_by",
        "holm": "holm",
        "bonferroni": "bonferroni",
    }
    sm_method = method_map.get(method, method)
    _, padj, _, _ = multipletests(pvals, alpha=alpha, method=sm_method)
    j = 0
    for r in results:
        if np.isnan(r.p_value):
            continue
        r.adjusted_p_value = float(padj[j])
        r.significant = bool(padj[j] <= alpha)
        j += 1


# Backward-compatible alias
_apply_fdr = apply_fdr


def bayesian_group_comparison(
    posterior_samples_by_sample: dict[str, np.ndarray],
    sample_groups: dict[str, str | None],
    cell_type_names: list[str],
    contrasts: list[tuple[str, str]],
    fdr_alpha: float = 0.05,
) -> list[TestResult]:
    """Posterior probability of difference (math doc §8.1).

    posterior_samples_by_sample: maps sample_id → (T, K+1) array of MCMC draws.
    For each contrast (A, B) and cell type j we draw paired group means and
    compute P(mean_A > mean_B). The pseudo p-value 2 * min(P, 1-P) is used
    for FDR-compatible reporting.
    """
    if not posterior_samples_by_sample or not contrasts:
        return []

    # Determine the common number of draws (truncate if mismatched)
    sample_ids = list(posterior_samples_by_sample.keys())
    min_T = min(arr.shape[0] for arr in posterior_samples_by_sample.values())
    K = len(cell_type_names)

    # Stack into (S, T, K+1)
    stacked = np.stack(
        [posterior_samples_by_sample[sid][:min_T] for sid in sample_ids], axis=0
    )

    results: list[TestResult] = []
    for j in range(K):
        for group_a, group_b in contrasts:
            a_mask = np.array(
                [sample_groups.get(sid) == group_a for sid in sample_ids]
            )
            b_mask = np.array(
                [sample_groups.get(sid) == group_b for sid in sample_ids]
            )
            if not (np.any(a_mask) and np.any(b_mask)):
                continue
            # Group-level mean per draw
            mean_a_draws = stacked[a_mask, :, j].mean(axis=0)  # (T,)
            mean_b_draws = stacked[b_mask, :, j].mean(axis=0)  # (T,)
            prob = float(np.mean(mean_a_draws > mean_b_draws))
            pseudo_p = float(2.0 * min(prob, 1.0 - prob))
            results.append(
                TestResult(
                    cell_type=cell_type_names[j],
                    contrast=f"{group_a}_vs_{group_b}",
                    test_type="bayesian",
                    mean_a=float(np.mean(mean_a_draws)),
                    mean_b=float(np.mean(mean_b_draws)),
                    effect_size=float(np.mean(mean_a_draws) - np.mean(mean_b_draws)),
                    se=float(
                        np.std(mean_a_draws - mean_b_draws, ddof=1)
                        if min_T > 1
                        else 0.0
                    ),
                    statistic=float(prob),
                    p_value=pseudo_p,
                )
            )

    _apply_fdr(results, fdr_alpha)
    return results


__all__ = [
    "TestResult",
    "apply_fdr",
    "bayesian_group_comparison",
    "compositional_regression_test",
    "wilcoxon_test",
]
