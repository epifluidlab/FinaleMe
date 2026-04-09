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

    Fits a per-cell-type OLS regression in ILR space:

        y_{j,s} = β_0 + Σ_g β_g · I(group_s = g) + Σ_c β_c · covariate_{c,s} + ε_s

    and tests the contrast ``β_A - β_B`` via a Wald t-test. The first
    group alphabetically is the reference (absorbed into the intercept),
    so group coefficients are "difference from reference". The contrast
    is recovered as ``β_A - β_B`` regardless of which group is the
    reference; when the reference is one of the two contrast groups,
    ``β_ref = 0`` and the contrast collapses to ``±β_other``.

    Covariates flow through a pandas DataFrame indexed by ``sample_id``
    (or keyed by the ``sample_id`` column if present). Numeric columns
    are used directly; categorical columns are one-hot-encoded. Rows
    with any NaN covariate value are dropped for THAT cell type's fit
    (so a partially-populated covariate doesn't kill the whole analysis).

    When no covariates are supplied the regression reduces to the
    intercept + group indicator design — which for the two-group case
    gives the same p-value as Welch's t-test on the ILR coordinate,
    matching the previous implementation as a special case.

    proportions: (S, K+1) — last column is the unknown component.
    contrasts: list of (group_a, group_b) pairs.
    """
    # Drop samples with no group label
    valid = [i for i, g in enumerate(group_labels) if g is not None]
    if len(valid) < 3:
        return []
    P = proportions[valid]
    groups = [group_labels[i] for i in valid]
    valid_sample_ids = [sample_ids[i] for i in valid]

    # ILR transform — sum_w_unknown is included so that the simplex is K+1 parts
    # The result has K coords; coord j is the contrast that most directly
    # involves cell type j (because of the Helmert basis).
    ilr_coords = ilr_transform(P)  # (S, K)
    n_cells = len(cell_type_names)
    S = ilr_coords.shape[0]

    # Pre-build per-sample covariate row vector (float64, with one-hot
    # encoding for non-numeric columns). None if no covariates were supplied.
    cov_matrix, cov_col_names = _build_covariate_design(
        covariates, valid_sample_ids
    )
    # cov_matrix is (S, C) or None. NaN entries are preserved here; per-cell
    # dropna happens inside the fitting loop so cell types with different
    # NaN patterns don't drag each other down.

    # Build the group encoding. Reference group is the first in sorted order.
    unique_groups = sorted({g for g in groups})
    reference_group = unique_groups[0] if unique_groups else None
    non_ref_groups = [g for g in unique_groups if g != reference_group]
    # group_design is (S, n_non_ref_groups), 1 if sample is in that group
    group_design = np.zeros((S, len(non_ref_groups)), dtype=np.float64)
    group_col_to_idx: dict[str, int] = {}
    for col_idx, g in enumerate(non_ref_groups):
        group_col_to_idx[g] = col_idx
        group_design[:, col_idx] = np.array(
            [1.0 if gs == g else 0.0 for gs in groups], dtype=np.float64
        )

    results: list[TestResult] = []

    for j in range(n_cells):
        y = ilr_coords[:, j]
        for group_a, group_b in contrasts:
            a_mask = np.array([g == group_a for g in groups])
            b_mask = np.array([g == group_b for g in groups])
            if int(a_mask.sum()) < 2 or int(b_mask.sum()) < 2:
                results.append(_skipped_result(cell_type_names[j], group_a, group_b))
                continue

            # Build the design matrix for THIS fit:
            #   columns = [intercept] + non_ref_group_dummies + covariates
            # and drop rows where any covariate is NaN.
            X_parts = [np.ones((S, 1), dtype=np.float64), group_design]
            if cov_matrix is not None:
                X_parts.append(cov_matrix)
            X_full = np.hstack(X_parts)
            n_cols = X_full.shape[1]

            row_valid = np.all(np.isfinite(X_full), axis=1) & np.isfinite(y)
            # Require enough samples for the regression to be identified
            # (more rows than columns plus some slack).
            if int(np.sum(row_valid)) < n_cols + 1:
                results.append(_skipped_result(cell_type_names[j], group_a, group_b))
                continue

            # Post-drop group-size re-check. The initial pre-drop guard
            # above confirms each contrast group has ≥2 samples in the
            # full cohort, but covariate NaN row-dropping can silently
            # shrink either group below that threshold. A 1-vs-N contrast
            # post-drop yields a statistically meaningless p-value (the
            # group coefficient and its standard error are both driven
            # by a single observation), so we re-check per-group sizes
            # on the filtered mask and skip any contrast where either
            # group has <2 samples remaining.
            a_mask_filtered = a_mask & row_valid
            b_mask_filtered = b_mask & row_valid
            n_a_after = int(a_mask_filtered.sum())
            n_b_after = int(b_mask_filtered.sum())
            if n_a_after < 2 or n_b_after < 2:
                results.append(_skipped_result(cell_type_names[j], group_a, group_b))
                continue

            X = X_full[row_valid]
            y_v = y[row_valid]
            n_obs = X.shape[0]
            dof = n_obs - n_cols
            if dof < 1:
                results.append(_skipped_result(cell_type_names[j], group_a, group_b))
                continue

            # OLS closed form: β̂ = (Xᵀ X)⁻¹ Xᵀ y
            try:
                XtX = X.T @ X
                XtX_inv = np.linalg.pinv(XtX)
                beta = XtX_inv @ (X.T @ y_v)
            except np.linalg.LinAlgError:
                results.append(_skipped_result(cell_type_names[j], group_a, group_b))
                continue

            residuals = y_v - X @ beta
            rss = float(np.sum(residuals ** 2))
            sigma2 = rss / dof if dof > 0 else float("inf")
            cov_beta = sigma2 * XtX_inv

            # Contrast vector: difference between β_A and β_B where β_*
            # refers to the one-hot column for that group (0 for the
            # reference group since it's absorbed into the intercept).
            contrast_vec = np.zeros(n_cols, dtype=np.float64)
            # First column is intercept (coefficient index 0);
            # non-reference group dummies start at index 1.
            group_col_offset = 1
            if group_a in group_col_to_idx:
                contrast_vec[group_col_offset + group_col_to_idx[group_a]] = 1.0
            if group_b in group_col_to_idx:
                contrast_vec[group_col_offset + group_col_to_idx[group_b]] = -1.0

            effect = float(contrast_vec @ beta)
            var_effect = float(contrast_vec @ cov_beta @ contrast_vec)
            if var_effect <= 0 or not np.isfinite(var_effect):
                results.append(_skipped_result(cell_type_names[j], group_a, group_b))
                continue
            se = float(np.sqrt(var_effect))
            t_stat = effect / se

            from scipy.stats import t as student_t

            p = float(2.0 * student_t.sf(abs(t_stat), df=dof))

            # Raw proportion means use the unfiltered (pre-NaN-drop) mask
            # so the reported mean is the cohort-average even if some
            # samples were excluded from THIS fit due to missing covariates.
            results.append(
                TestResult(
                    cell_type=cell_type_names[j],
                    contrast=f"{group_a}_vs_{group_b}",
                    test_type="ilr_regression",
                    mean_a=float(np.mean(P[a_mask, j])),
                    mean_b=float(np.mean(P[b_mask, j])),
                    effect_size=effect,
                    se=se,
                    statistic=t_stat,
                    p_value=p,
                    extra={
                        "n_obs": n_obs,
                        "n_a_after_dropna": n_a_after,
                        "n_b_after_dropna": n_b_after,
                        "n_covariates": (
                            len(cov_col_names) if cov_col_names else 0
                        ),
                        "dof": dof,
                        "covariate_columns": list(cov_col_names or []),
                    },
                )
            )

    _apply_fdr(results, fdr_alpha)
    return results


def _build_covariate_design(
    covariates: pd.DataFrame | None,
    sample_ids: list[str],
) -> tuple[np.ndarray | None, list[str] | None]:
    """Build an (S, C) float64 covariate matrix aligned to ``sample_ids``.

    Numeric columns pass through unchanged. Categorical columns are
    one-hot-encoded with the first category dropped (treatment encoding)
    so the design matrix stays full-rank alongside the intercept. The
    ``sample_id`` column (if present) is consumed for alignment and
    dropped from the output matrix.

    NaN values are preserved in the output — the caller handles per-fit
    row dropping so one missing covariate doesn't cascade across all
    tests. Returns ``(None, None)`` when no covariates are supplied.
    """
    if covariates is None or covariates.empty:
        return None, None

    cov = covariates.copy()
    if "sample_id" in cov.columns:
        cov = cov.set_index("sample_id")
    # Reindex onto our sample order; missing rows become all-NaN.
    try:
        cov = cov.reindex(sample_ids)
    except Exception:
        # Index mismatch (non-string, duplicates, ...) — give up gracefully.
        return None, None

    # Drop columns that are entirely missing
    cov = cov.dropna(axis=1, how="all")
    if cov.empty:
        return None, None

    # One-hot-encode non-numeric columns (drop_first so the intercept
    # captures the reference category and the design is full-rank).
    cov_encoded = pd.get_dummies(cov, drop_first=True, dummy_na=False)
    if cov_encoded.empty:
        return None, None

    col_names = list(cov_encoded.columns)
    try:
        M = cov_encoded.to_numpy(dtype=np.float64)
    except (ValueError, TypeError):
        # Some column still isn't numeric after encoding (e.g. object dtype)
        # — drop those and retry.
        numeric_cols = [
            c for c in cov_encoded.columns
            if pd.api.types.is_numeric_dtype(cov_encoded[c])
        ]
        if not numeric_cols:
            return None, None
        cov_encoded = cov_encoded[numeric_cols]
        col_names = numeric_cols
        M = cov_encoded.to_numpy(dtype=np.float64)
    return M, col_names


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
