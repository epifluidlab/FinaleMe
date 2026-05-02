#!/usr/bin/env python3
"""too_diff_analysis.py — Cohort-level differential cell-type analysis from
BetaValueDeconvolution output.

Takes per-sample TOO deconvolution proportions (with optional bootstrap CIs)
and produces a per-cell-type group comparison adjusted for clinical and
technical covariates.

Workflow
--------
1. Load per-sample point estimates from one combined long-format TSV or
   from multiple per-sample TSVs (auto-concatenated).
2. Pivot to wide (samples x cell_types) and optionally refine NNLS-floored
   zeros using the bootstrap CI upper bound.
3. Centered log-ratio (CLR) transform — converts proportions on the
   simplex to unconstrained real values, suitable for standard linear
   models.
4. For each cell type, fit a linear model:
     clr(w_c) ~ group + clinical_covariates + technical_covariates
   Either omnibus (F-test across all group levels) or pairwise vs a
   reference group, configurable via --test.
5. BH-correct group-effect p-values across cell types.

Why not just compare per-sample "YES"/"NO" significance flags
------------------------------------------------------------
The per-sample q-values from BetaValueDeconvolution test
H0: w_c = 0 in THIS sample — they're QC diagnostics for individual
samples. The cohort question is "do groups differ in cell type
composition?" — a different inference. Per-sample binary thresholding
(a) loses magnitude, (b) is sensitive to noise near the cutoff, and
(c) doesn't model between-sample variance. The CLR + linear model
approach uses the actual point-estimate magnitudes as the data and
correctly handles the simplex constraint.

Example
-------
    python scripts/too_diff_analysis.py \\
        --deconv-tsv deconv_stats.tsv \\
        --metadata samples.tsv \\
        --group-col disease_status \\
        --reference-group Control \\
        --covariates age,sex,library_batch \\
        --output diff_results.tsv

    # Or merge per-sample outputs:
    python scripts/too_diff_analysis.py \\
        --deconv-files deconv_run/*.tsv \\
        --metadata samples.tsv \\
        --group-col disease_status \\
        --output diff_results.tsv
"""
import argparse
import glob
import sys
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

try:
    import statsmodels.api as sm  # noqa: F401
    import statsmodels.formula.api as smf
    from statsmodels.stats.anova import anova_lm
except ImportError as exc:  # pragma: no cover
    sys.exit(
        f"too_diff_analysis.py requires statsmodels (pip install statsmodels): {exc}"
    )

try:
    from scipy import stats as _scipy_stats
    from scipy.special import digamma as _digamma, polygamma as _polygamma
    from scipy.optimize import brentq as _brentq
except ImportError as exc:  # pragma: no cover
    sys.exit(
        f"too_diff_analysis.py requires scipy (pip install scipy): {exc}"
    )


# -------- IO ---------------------------------------------------------------

REQUIRED_COLS = {"sample", "cell_type", "proportion"}
OPTIONAL_COLS = {"CI_lower", "CI_upper", "p_value", "q_value", "significant"}


def load_deconv(deconv_tsv: Optional[str], deconv_files: Optional[Iterable[str]]) -> pd.DataFrame:
    """Load BetaValueDeconvolution long-format output.

    Accepts either a single combined TSV (when -bootstrap was run on
    multiple query files in one invocation) or a list of per-sample TSVs
    that get concatenated. If a per-sample file lacks a `sample` column,
    the basename (sans .tsv) is used as the sample id.
    """
    if deconv_tsv and deconv_files:
        raise ValueError("Provide either --deconv-tsv or --deconv-files, not both.")

    if deconv_tsv:
        df = pd.read_csv(deconv_tsv, sep="\t", comment="#")
        return df

    frames = []
    for f in deconv_files:
        df_i = pd.read_csv(f, sep="\t", comment="#")
        if "sample" not in df_i.columns:
            stem = Path(f).name
            for ext in (".tsv.gz", ".tsv"):
                if stem.endswith(ext):
                    stem = stem[: -len(ext)]
                    break
            df_i["sample"] = stem
        frames.append(df_i)
    if not frames:
        raise ValueError("No input files found.")
    return pd.concat(frames, ignore_index=True)


def pivot_long_to_wide(df: pd.DataFrame, value_col: str) -> pd.DataFrame:
    """Pivot long-format (sample, cell_type, value) into samples x cell_types."""
    if value_col not in df.columns:
        raise KeyError(f"Column {value_col!r} not in deconvolution output.")
    return df.pivot(index="sample", columns="cell_type", values=value_col)


# -------- Compositional preprocessing -------------------------------------

def refine_zeros_with_ci(
    w: pd.DataFrame, ci_upper: pd.DataFrame, min_ci: float = 1e-4
) -> pd.DataFrame:
    """Replace exact-zero point estimates with CI_upper/2 where CI_upper > min_ci.

    NNLS can drive small genuine fractions to exactly 0 because of the
    non-negative constraint. The bootstrap CI catches this: if the point
    estimate is 0 but the 97.5th percentile of bootstrap is non-trivial,
    the cell type is "below detection" rather than "absent". Replacing
    with CI_upper/2 keeps the magnitude information for the log-ratio
    transform without inventing precision.
    """
    out = w.copy()
    aligned = ci_upper.reindex(index=w.index, columns=w.columns)
    mask = (out == 0) & (aligned > min_ci)
    out = out.where(~mask, aligned / 2.0)
    return out


def clr_transform(w: pd.DataFrame, eps: float = 1e-6) -> pd.DataFrame:
    """Centered log-ratio transform.

        clr(w_i, c) = log(w_i,c) - (1/K) * sum_j log(w_i,j)

    Adds eps to zeros and renormalizes rows to sum to 1 before logging,
    so the result is finite for any non-negative input. Output is
    unconstrained real-valued; standard linear models apply.
    """
    if w.shape[1] < 2:
        raise ValueError("CLR requires at least 2 cell types.")
    w_safe = w.copy().astype(float)
    w_safe = w_safe.where(w_safe > eps, eps)
    w_safe = w_safe.div(w_safe.sum(axis=1), axis=0)
    log_w = np.log(w_safe)
    return log_w.sub(log_w.mean(axis=1), axis=0)


# -------- Modeling --------------------------------------------------------

def _quote_term(name: str, df: pd.DataFrame, force_categorical: bool = False) -> str:
    """Wrap a column name for a Patsy formula, marking categoricals."""
    safe = f'Q("{name}")'
    if force_categorical or df[name].dtype == object or str(df[name].dtype) == "category":
        return f"C({safe})"
    return safe


def build_formula(
    response: str, group_col: str, covariates: list[str], df: pd.DataFrame
) -> str:
    """Build a Patsy formula for `response ~ group + covariates`."""
    rhs_terms = [_quote_term(group_col, df, force_categorical=True)]
    for cov in covariates:
        if cov not in df.columns:
            raise KeyError(f"Covariate {cov!r} not found in metadata columns.")
        rhs_terms.append(_quote_term(cov, df))
    return f'Q("{response}") ~ ' + " + ".join(rhs_terms)


def fit_one_celltype(
    cell_type: str,
    df: pd.DataFrame,
    group_col: str,
    covariates: list[str],
    test: str,
    weights: Optional[pd.Series] = None,
) -> dict:
    """Fit a linear model and extract the group-effect statistics.

    test == 'omnibus': F-test for the group factor (any group differs).
        Returns the omnibus p-value and the maximum-magnitude pairwise
        effect for diagnostic purposes.
    test == 'pairwise': returns the per-level (vs reference) coefficients
        and their p-values for each non-reference group.
    """
    formula = build_formula(cell_type, group_col, covariates, df)
    fit_kwargs = {"data": df}
    if weights is not None:
        # WLS with provided weights. statsmodels expects 1D weights aligned
        # with the row order of `df`.
        fit_kwargs["weights"] = weights.reindex(df.index).values
        model = smf.wls(formula, **fit_kwargs)
    else:
        model = smf.ols(formula, **fit_kwargs)
    fit = model.fit()

    # Identify the group's coefficients (Patsy names them like
    # 'C(Q("group"))[T.<level>]' for the non-reference levels)
    group_pattern = f'C(Q("{group_col}"))'
    group_params = [t for t in fit.params.index if t.startswith(group_pattern)]

    # Residual diagnostics needed by empirical-Bayes shrinkage and permutation.
    # mse_resid = SSR / df_resid (statsmodels handles WLS weights correctly).
    base = {
        "cell_type": cell_type,
        "n_samples": int(fit.nobs),
        "r_squared": float(fit.rsquared),
        "s2_resid": float(fit.mse_resid) if fit.df_resid > 0 else np.nan,
        "df_resid": float(fit.df_resid),
    }

    if test == "omnibus":
        # F-test on the group factor only (Type II via anova_lm with typ=2)
        try:
            anova = anova_lm(fit, typ=2)
            row_label = next(
                (idx for idx in anova.index if idx.startswith(group_pattern)), None
            )
            if row_label is None:
                base["p_value"] = np.nan
                base["F_stat"] = np.nan
                base["df_group"] = np.nan
            else:
                base["p_value"] = float(anova.loc[row_label, "PR(>F)"])
                base["F_stat"] = float(anova.loc[row_label, "F"])
                base["df_group"] = float(anova.loc[row_label, "df"])
        except Exception as exc:  # pragma: no cover
            print(f"[{cell_type}] anova failed: {exc}", file=sys.stderr)
            base["p_value"] = np.nan
            base["F_stat"] = np.nan
            base["df_group"] = np.nan

        # Reporting effect: maximum-magnitude per-level coefficient
        if group_params:
            mags = [(t, float(fit.params[t])) for t in group_params]
            top_term, top_coef = max(mags, key=lambda x: abs(x[1]))
            level = _strip_level(top_term, group_pattern)
            base["max_effect_clr"] = top_coef
            base["max_effect_level"] = level
        else:
            base["max_effect_clr"] = np.nan
            base["max_effect_level"] = None
        return base

    # pairwise
    rows = []
    for term in group_params:
        level = _strip_level(term, group_pattern)
        rows.append(
            {
                **base,
                "level": level,
                "effect_clr": float(fit.params[term]),
                "std_err": float(fit.bse[term]),
                "p_value": float(fit.pvalues[term]),
            }
        )
    if not rows:
        rows.append({**base, "level": None, "effect_clr": np.nan,
                     "std_err": np.nan, "p_value": np.nan})
    return rows  # list of dicts (one per non-reference group level)


def _strip_level(term: str, group_pattern: str) -> str:
    """Extract the level name from a Patsy term like 'C(Q("group"))[T.Disease]'."""
    if "[T." in term:
        return term.split("[T.", 1)[1].rstrip("]")
    return term[len(group_pattern):].lstrip(".[]")


# -------- Empirical-Bayes variance shrinkage (limma / Smyth 2004) ----------

def smyth_eb_prior(s2: np.ndarray, df: np.ndarray) -> tuple[float, float]:
    """Estimate the prior variance s0^2 and prior degrees of freedom d0
    from a vector of per-cell-type residual variances and df.

    Implements Smyth 2004 (limma) — fits a scaled chi-square prior to the
    observed residual variances on the log scale, using the trigamma
    moment-matching equation
        var(log s2) - mean(trigamma(df/2)) = trigamma(d0/2)
    Returns (s0_sq, d0). Falls back to (mean(s2), inf) — i.e. infinite
    shrinkage to the pooled mean — when the moment estimator is degenerate
    (very few cell types or near-uniform variances).
    """
    s2 = np.asarray(s2, dtype=float)
    df = np.asarray(df, dtype=float)
    valid = (s2 > 0) & np.isfinite(s2) & np.isfinite(df) & (df > 0)
    if valid.sum() < 3:
        fallback = float(np.nanmean(s2[valid])) if valid.any() else 1.0
        return fallback, float("inf")
    s2v = s2[valid]
    dfv = df[valid]
    z = np.log(s2v)
    e = z - _digamma(dfv / 2.0) + np.log(dfv / 2.0)
    e_mean = float(e.mean())
    e_var = float(e.var(ddof=1))
    target = e_var - float(_polygamma(1, dfv / 2.0).mean())
    if target <= 0:
        # Observed between-cell-type variance is fully explained by within-
        # group sampling noise → prior is infinitely informative.
        return float(np.exp(e_mean)), float("inf")
    # Solve trigamma(x) = target for x = d0/2. Trigamma is strictly
    # decreasing on (0, inf), so we can bracket the root.
    def f(x):
        return float(_polygamma(1, x)) - target
    try:
        if f(1e-3) * f(1e6) > 0:
            return float(np.exp(e_mean)), float("inf")
        x = _brentq(f, 1e-3, 1e6, xtol=1e-6)
    except Exception:
        return float(np.exp(e_mean)), float("inf")
    d0 = 2.0 * float(x)
    log_s0_sq = e_mean + float(_digamma(x)) - float(np.log(x))
    s0_sq = float(np.exp(log_s0_sq))
    return s0_sq, d0


def moderated_t_pvalue(
    beta: float, se: float, df_resid: float, s2: float, d0: float, s0_sq: float
) -> tuple[float, float]:
    """Two-sided moderated t p-value (limma).

    Returns (moderated_t, p_value). df_total = df_resid + d0.
    """
    if not (np.isfinite(beta) and np.isfinite(se) and se > 0
            and np.isfinite(s2) and s2 > 0):
        return float("nan"), float("nan")
    if not np.isfinite(d0):
        # Prior dominates entirely → use pooled prior variance.
        s2_tilde = s0_sq
        df_total = 1e10
    else:
        s2_tilde = (d0 * s0_sq + df_resid * s2) / (d0 + df_resid)
        df_total = df_resid + d0
    se_mod = se * np.sqrt(s2_tilde / s2)
    t_mod = beta / se_mod
    p = 2.0 * float(_scipy_stats.t.sf(abs(t_mod), df_total))
    return float(t_mod), p


def moderated_F_pvalue(
    F_obs: float, df_num: float, df_resid: float,
    s2: float, d0: float, s0_sq: float,
) -> tuple[float, float]:
    """Moderated F-test p-value: F_mod = F_obs * s2 / s2_tilde, df_denom = df_resid + d0."""
    if not (np.isfinite(F_obs) and np.isfinite(df_num)
            and np.isfinite(s2) and s2 > 0):
        return float("nan"), float("nan")
    if not np.isfinite(d0):
        s2_tilde = s0_sq
        df_total = 1e10
    else:
        s2_tilde = (d0 * s0_sq + df_resid * s2) / (d0 + df_resid)
        df_total = df_resid + d0
    F_mod = F_obs * s2 / s2_tilde
    p = float(_scipy_stats.f.sf(F_mod, df_num, df_total))
    return float(F_mod), p


# -------- Permutation-based group-label null ------------------------------

def permutation_pvalues(
    cell_types: list[str],
    df: pd.DataFrame,
    group_col: str,
    covariates: list[str],
    test: str,
    weights: Optional[pd.Series],
    obs_results: list,
    n_perm: int,
    seed: int,
    verbose: bool = False,
) -> dict:
    """Compute empirical p-values by permuting the group-label column.

    Strategy: simple Manly-style permutation of the group column. The
    response (CLR-transformed proportions) and covariates stay aligned.
    For each permutation b in [1..B]:
       - shuffle df[group_col]
       - refit each cell-type model
       - record the test stat (F for omnibus, |t| for pairwise per level)
    Empirical p = (#{|stat_perm| >= |stat_obs|} + 1) / (B + 1).

    Caveat: simple group-label permutation is exact under the strong null
    (no covariate effect). With strong covariates the test is slightly
    anti-conservative; if you have non-trivial covariates use parametric
    or empirical-Bayes p-values as the primary, and treat permutation
    as a robustness check.

    Returns dict mapping cell_type -> p-value (omnibus) or
    cell_type -> {level: p-value} (pairwise).
    """
    rng = np.random.default_rng(seed)

    # Build observed-statistic lookup from the already-fit results.
    # obs_results matches the row layout of `results` dataframe.
    obs_stats: dict = {}
    for res in obs_results:
        if test == "omnibus":
            obs_stats[res["cell_type"]] = res.get("F_stat", float("nan"))
        else:
            ct = res["cell_type"]
            lvl = res.get("level")
            beta = res.get("effect_clr", float("nan"))
            se = res.get("std_err", float("nan"))
            t_obs = abs(beta / se) if (se and np.isfinite(se) and se > 0
                                       and np.isfinite(beta)) else float("nan")
            obs_stats.setdefault(ct, {})[lvl] = t_obs

    counts: dict = {}
    for ct in cell_types:
        counts[ct] = 0 if test == "omnibus" else {}

    # Pre-seed pairwise level counts to 0 for every observed (cell_type, level)
    # pair, so cell types where NO permutation beats |t_obs| still get
    # p = 1/(N+1) rather than NaN from a missing dict key.
    if test != "omnibus":
        for res in obs_results:
            ct = res["cell_type"]
            lvl = res.get("level")
            if lvl is not None and ct in counts:
                counts[ct].setdefault(lvl, 0)

    # Preserve the Categorical dtype (with the user's --reference-group
    # ordering) across permutations. A plain numpy assignment would coerce
    # to object dtype, after which Patsy reverts to alphabetical ordering
    # and the reference level changes — making the per-level counts diverge
    # from the observed-statistic keys and yielding NaN p-values.
    group_series = df[group_col]
    is_categorical = isinstance(group_series.dtype, pd.CategoricalDtype)
    group_categories = group_series.cat.categories if is_categorical else None
    group_vals = group_series.to_numpy(copy=True)

    df_perm = df.copy()

    for b in range(n_perm):
        perm_arr = rng.permutation(group_vals)
        if is_categorical:
            df_perm[group_col] = pd.Categorical(
                perm_arr, categories=group_categories
            )
        else:
            df_perm[group_col] = perm_arr
        for ct in cell_types:
            try:
                res = fit_one_celltype(
                    ct, df_perm, group_col, covariates, test, weights
                )
            except Exception:
                continue
            if test == "omnibus":
                f_perm = res.get("F_stat", float("nan"))
                f_obs = obs_stats.get(ct, float("nan"))
                if (np.isfinite(f_perm) and np.isfinite(f_obs)
                        and f_perm >= f_obs):
                    counts[ct] += 1
            else:
                for r in res:
                    lvl = r.get("level")
                    se = r.get("std_err", float("nan"))
                    beta = r.get("effect_clr", float("nan"))
                    if not (se and np.isfinite(se) and se > 0
                            and np.isfinite(beta)):
                        continue
                    t_perm = abs(beta / se)
                    t_obs = obs_stats.get(ct, {}).get(lvl, float("nan"))
                    if np.isfinite(t_perm) and np.isfinite(t_obs) and t_perm >= t_obs:
                        counts[ct][lvl] = counts[ct].get(lvl, 0) + 1
        if verbose and (b + 1) % max(1, n_perm // 10) == 0:
            print(f"  permutation {b + 1}/{n_perm}", file=sys.stderr)

    pvals: dict = {}
    for ct in cell_types:
        if test == "omnibus":
            pvals[ct] = (counts[ct] + 1) / (n_perm + 1)
        else:
            pvals[ct] = {
                lvl: (c + 1) / (n_perm + 1) for lvl, c in counts[ct].items()
            }
    return pvals


# -------- BH FDR correction (no scipy.stats.false_discovery_control dep) ---

def bh_correct(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values from a 1D array of p-values.

    NaNs are passed through (q = NaN). Result has the same shape as input.
    """
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return out
    pv = p[valid]
    n = pv.size
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * n / (np.arange(1, n + 1))
    # Step-up: cumulative minimum from the right
    q_min = np.minimum.accumulate(q[::-1])[::-1]
    q_min = np.clip(q_min, 0.0, 1.0)
    q_orig = np.empty_like(q_min)
    q_orig[order] = q_min
    out[valid] = q_orig
    return out


# -------- Main ------------------------------------------------------------

def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument(
        "--deconv-tsv",
        help="Single combined long-format TSV from BetaValueDeconvolution -bootstrap.",
    )
    src.add_argument(
        "--deconv-files",
        nargs="+",
        help="Per-sample long-format TSVs (one per sample). Glob patterns OK.",
    )

    p.add_argument(
        "--metadata",
        required=True,
        help="Metadata TSV with at minimum 'sample' and the --group-col column.",
    )
    p.add_argument(
        "--group-col",
        default="group",
        help="Metadata column holding the group label (default: 'group').",
    )
    p.add_argument(
        "--reference-group",
        default=None,
        help="Group level to use as the reference. If omitted, the "
             "alphabetically first level is used. Required for --test pairwise.",
    )
    p.add_argument(
        "--covariates",
        default="",
        help="Comma-separated metadata columns to adjust for (e.g. "
             "'age,sex,library_batch,coverage_tier').",
    )
    p.add_argument(
        "--test",
        choices=["omnibus", "pairwise"],
        default="omnibus",
        help="omnibus (default): one F-test per cell type for any group "
             "difference; suitable for >=2 groups. pairwise: each "
             "non-reference group vs the reference, one row per (cell_type, level).",
    )
    p.add_argument(
        "--exclude-cell-types",
        default="",
        help="Comma-separated cell types to skip (e.g. 'Unknown' if input "
             "comes from a tool that estimates Unknown).",
    )
    p.add_argument(
        "--min-prevalence",
        type=float,
        default=0.0,
        help="Drop cell types whose proportion exceeds --min-prevalence-thresh "
             "in fewer than this fraction of samples (e.g. 0.2 keeps cell "
             "types non-zero in >=20%% of samples). Reduces multiple-testing "
             "burden by removing cell types with no signal. Default: 0 (off).",
    )
    p.add_argument(
        "--min-prevalence-thresh",
        type=float,
        default=1e-4,
        help="Threshold above which a sample is counted as having the cell "
             "type 'present' for --min-prevalence filtering. Default: 1e-4.",
    )
    p.add_argument(
        "--refine-zeros-with-ci",
        action="store_true",
        help="Replace NNLS-floored zeros with CI_upper/2 where CI_upper "
             "exceeds --refine-zeros-ci-min. Requires CI_upper in input.",
    )
    p.add_argument(
        "--refine-zeros-ci-min",
        type=float,
        default=1e-4,
        help="Threshold for CI_upper above which a zero point estimate is "
             "treated as below-detection (refined to CI_upper/2). Default: 1e-4.",
    )
    p.add_argument(
        "--clr-eps",
        type=float,
        default=1e-6,
        help="Pseudocount added to zero proportions before log. Default: 1e-6.",
    )
    p.add_argument(
        "--weighted",
        action="store_true",
        help="Inverse-variance-weight samples by mean CI width across cell "
             "types. Down-weights low-quality samples. Requires CI columns.",
    )
    p.add_argument(
        "--fdr-alpha",
        type=float,
        default=0.05,
        help="q-value threshold for the 'significant' flag (default: 0.05).",
    )
    p.add_argument(
        "--empirical-bayes",
        action="store_true",
        help="Apply limma-style empirical-Bayes shrinkage (Smyth 2004) to "
             "per-cell-type residual variances before computing p-values. "
             "Pools variance information across cell types via a scaled "
             "chi-square prior, yielding a moderated t (pairwise) or F "
             "(omnibus) statistic with df = df_resid + d0. Substantially "
             "improves power with small N. The original parametric p-values "
             "are kept in 'p_value_unshrunk'.",
    )
    p.add_argument(
        "--permutations",
        type=int,
        default=0,
        help="If >0, also compute empirical p-values by permuting the group "
             "label that many times. Adds 'p_value_perm' and 'q_value_perm' "
             "columns. Robust to non-Gaussian residuals at the cost of "
             "B*K refits (B=permutations, K=cell types). Recommend "
             ">=1000 for q<=0.05 resolution.",
    )
    p.add_argument(
        "--permutation-seed",
        type=int,
        default=0,
        help="RNG seed for --permutations (default: 0, deterministic).",
    )
    p.add_argument(
        "--output",
        required=True,
        help="Output TSV path. One row per cell type (omnibus) or per "
             "(cell_type, level) pair (pairwise).",
    )
    p.add_argument(
        "--verbose",
        action="store_true",
        help="Print loaded sample counts, fit warnings, and a summary table.",
    )

    args = p.parse_args()

    # ---- Load deconv ---------------------------------------------------
    if args.deconv_files:
        # Expand any globs the shell didn't already
        expanded = []
        for f in args.deconv_files:
            matched = glob.glob(f)
            expanded.extend(matched if matched else [f])
        deconv = load_deconv(None, expanded)
    else:
        deconv = load_deconv(args.deconv_tsv, None)

    missing = REQUIRED_COLS - set(deconv.columns)
    if missing:
        sys.exit(
            f"ERROR: deconv input is missing required columns: {sorted(missing)}\n"
            f"Got columns: {sorted(deconv.columns)}\n"
            f"Was BetaValueDeconvolution run with -bootstrap (long format)?"
        )

    # ---- Pivot ---------------------------------------------------------
    w = pivot_long_to_wide(deconv, "proportion")

    # Optional zero refinement using CI_upper
    if args.refine_zeros_with_ci:
        if "CI_upper" not in deconv.columns:
            print(
                "WARNING: --refine-zeros-with-ci requested but CI_upper missing; "
                "skipping refinement.",
                file=sys.stderr,
            )
        else:
            ci_upper = pivot_long_to_wide(deconv, "CI_upper")
            w = refine_zeros_with_ci(w, ci_upper, min_ci=args.refine_zeros_ci_min)

    # Drop excluded cell types
    excluded = [c.strip() for c in args.exclude_cell_types.split(",") if c.strip()]
    if excluded:
        keep = [c for c in w.columns if c not in excluded]
        dropped = sorted(set(w.columns) - set(keep))
        if dropped and args.verbose:
            print(f"Excluded cell types: {dropped}", file=sys.stderr)
        w = w[keep]

    # Drop cell types with insufficient prevalence across samples
    if args.min_prevalence > 0:
        prev = (w > args.min_prevalence_thresh).mean(axis=0)
        keep = prev[prev >= args.min_prevalence].index.tolist()
        dropped = sorted(set(w.columns) - set(keep))
        if dropped and args.verbose:
            print(
                f"Dropped {len(dropped)} cell type(s) with prevalence < "
                f"{args.min_prevalence:.2f} (proportion > "
                f"{args.min_prevalence_thresh:g}): {dropped}",
                file=sys.stderr,
            )
        if not keep:
            sys.exit(
                f"ERROR: --min-prevalence {args.min_prevalence} dropped all "
                f"cell types. Lower the threshold or check the input."
            )
        w = w[keep]

    # ---- CLR ----------------------------------------------------------
    w_clr = clr_transform(w, eps=args.clr_eps)

    # ---- Load metadata ------------------------------------------------
    metadata = pd.read_csv(args.metadata, sep="\t", comment="#")
    if "sample" not in metadata.columns:
        sys.exit("ERROR: metadata must have a 'sample' column.")
    metadata = metadata.set_index("sample")
    if args.group_col not in metadata.columns:
        sys.exit(
            f"ERROR: metadata is missing --group-col {args.group_col!r}. "
            f"Available columns: {sorted(metadata.columns)}"
        )

    # Set reference group if requested
    if args.reference_group is not None:
        levels = list(metadata[args.group_col].dropna().unique())
        if args.reference_group not in levels:
            sys.exit(
                f"ERROR: --reference-group {args.reference_group!r} not found in "
                f"{args.group_col} levels: {sorted(levels)}"
            )
        ordered = [args.reference_group] + sorted(
            [g for g in levels if g != args.reference_group]
        )
        metadata = metadata.copy()
        metadata[args.group_col] = pd.Categorical(
            metadata[args.group_col], categories=ordered
        )

    # Restrict to samples present in both
    common = w_clr.index.intersection(metadata.index)
    if len(common) == 0:
        sys.exit(
            "ERROR: no overlap between deconv sample IDs and metadata sample IDs. "
            "Make sure 'sample' values match exactly."
        )
    if args.verbose:
        print(
            f"Samples in deconv: {w_clr.shape[0]:,}; "
            f"in metadata: {metadata.shape[0]:,}; "
            f"intersection used: {len(common):,}",
            file=sys.stderr,
        )
        group_counts = metadata.loc[common, args.group_col].value_counts()
        print(f"Group sizes:\n{group_counts.to_string()}", file=sys.stderr)
        print(f"Cell types: {w_clr.shape[1]:,}", file=sys.stderr)

    df = w_clr.loc[common].join(metadata.loc[common])

    # Drop samples with NaN in any covariate or group label
    needed = [args.group_col] + [
        c.strip() for c in args.covariates.split(",") if c.strip()
    ]
    before = len(df)
    df = df.dropna(subset=needed)
    if args.verbose and before != len(df):
        print(
            f"Dropped {before - len(df)} sample(s) with NaN in group/covariates.",
            file=sys.stderr,
        )

    # Compute per-sample weights if requested
    weights = None
    if args.weighted:
        if "CI_lower" not in deconv.columns or "CI_upper" not in deconv.columns:
            print(
                "WARNING: --weighted requires CI_lower/CI_upper; falling back "
                "to unweighted.",
                file=sys.stderr,
            )
        else:
            lo = pivot_long_to_wide(deconv, "CI_lower").reindex(df.index)
            hi = pivot_long_to_wide(deconv, "CI_upper").reindex(df.index)
            mean_width = (hi - lo).mean(axis=1)
            weights = 1.0 / (mean_width.pow(2) + 1e-9)

    # ---- Fit per cell type --------------------------------------------
    covariates = [c.strip() for c in args.covariates.split(",") if c.strip()]
    cell_types = list(w_clr.columns)

    rows: list[dict] = []
    for ct in cell_types:
        try:
            res = fit_one_celltype(
                ct, df, args.group_col, covariates, args.test, weights
            )
        except Exception as exc:  # pragma: no cover
            if args.verbose:
                print(f"WARNING: fit failed for {ct}: {exc}", file=sys.stderr)
            res = {
                "cell_type": ct,
                "p_value": np.nan,
                "n_samples": 0,
            }
        if isinstance(res, list):
            rows.extend(res)
        else:
            rows.append(res)

    results = pd.DataFrame(rows)

    # ---- Empirical-Bayes variance shrinkage (limma-style) ------------
    # Recomputes p_value using moderated t (pairwise) / moderated F (omnibus)
    # with a scaled chi-square variance prior pooled across cell types.
    # The original parametric p-value is preserved as p_value_unshrunk.
    if args.empirical_bayes:
        # For pairwise output we have one row per (cell_type, level) but
        # s2_resid/df_resid are properties of the cell-type fit, so dedupe
        # on cell_type for prior estimation.
        prior_table = (
            results[["cell_type", "s2_resid", "df_resid"]]
            .drop_duplicates("cell_type")
        )
        s0_sq, d0 = smyth_eb_prior(
            prior_table["s2_resid"].to_numpy(),
            prior_table["df_resid"].to_numpy(),
        )
        if args.verbose:
            print(
                f"EB prior (Smyth): s0^2={s0_sq:.4g}, d0={d0:.4g}, "
                f"n_celltypes={len(prior_table)}",
                file=sys.stderr,
            )
        results["p_value_unshrunk"] = results["p_value"].copy()
        new_p: list[float] = []
        for _, row in results.iterrows():
            s2 = row["s2_resid"]
            df_r = row["df_resid"]
            if args.test == "omnibus":
                _, p_mod = moderated_F_pvalue(
                    row.get("F_stat", float("nan")),
                    row.get("df_group", float("nan")),
                    df_r, s2, d0, s0_sq,
                )
            else:
                _, p_mod = moderated_t_pvalue(
                    row.get("effect_clr", float("nan")),
                    row.get("std_err", float("nan")),
                    df_r, s2, d0, s0_sq,
                )
            new_p.append(p_mod)
        results["p_value"] = new_p
        results["eb_s0_sq"] = s0_sq
        results["eb_d0"] = d0

    # ---- Permutation-based group-label null --------------------------
    if args.permutations > 0:
        if args.verbose:
            print(
                f"Running {args.permutations} permutations across "
                f"{len(cell_types)} cell types...",
                file=sys.stderr,
            )
        perm_p = permutation_pvalues(
            cell_types=cell_types,
            df=df,
            group_col=args.group_col,
            covariates=covariates,
            test=args.test,
            weights=weights,
            obs_results=rows,
            n_perm=args.permutations,
            seed=args.permutation_seed,
            verbose=args.verbose,
        )
        if args.test == "omnibus":
            results["p_value_perm"] = results["cell_type"].map(perm_p)
        else:
            results["p_value_perm"] = results.apply(
                lambda r: perm_p.get(r["cell_type"], {}).get(
                    r.get("level"), float("nan")
                ),
                axis=1,
            )
        results["q_value_perm"] = bh_correct(
            results["p_value_perm"].to_numpy()
        )

    # ---- BH FDR -------------------------------------------------------
    results["q_value"] = bh_correct(results["p_value"].to_numpy())
    results["significant"] = results["q_value"] <= args.fdr_alpha

    # Sort by q-value (significant cell types first)
    results = results.sort_values("q_value", na_position="last").reset_index(drop=True)

    # ---- Output -------------------------------------------------------
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.output, sep="\t", index=False, na_rep="NA")

    n_sig = int(results["significant"].sum())
    print(
        f"Wrote {args.output}: {len(results)} row(s), {n_sig} significant "
        f"at q<={args.fdr_alpha}.",
        file=sys.stderr,
    )
    if args.verbose and n_sig > 0:
        cols = [c for c in [
            "cell_type", "level", "effect_clr", "max_effect_clr",
            "max_effect_level", "p_value", "q_value", "n_samples",
        ] if c in results.columns]
        print("\nTop significant cell types:", file=sys.stderr)
        print(
            results[results["significant"]][cols].head(15).to_string(index=False),
            file=sys.stderr,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
