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

# matplotlib is only needed when --plot-pdf is set; defer import + error
# until the user actually requests plotting so the rest of the script
# stays usable on systems without matplotlib.
def _import_matplotlib():
    try:
        import matplotlib
        matplotlib.use("Agg")  # non-interactive backend (writes files only)
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        return matplotlib, plt, PdfPages
    except ImportError as exc:
        sys.exit(
            f"--plot-pdf requires matplotlib (pip install matplotlib): {exc}"
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

def _winsorize(arr: np.ndarray, k: float = 1.5) -> np.ndarray:
    """Symmetric IQR-based Winsorization: clip values outside median +- k*IQR."""
    arr = np.asarray(arr, dtype=float)
    if arr.size < 4:
        return arr
    q1, q3 = np.percentile(arr, [25, 75])
    iqr = q3 - q1
    if not np.isfinite(iqr) or iqr <= 0:
        return arr
    med = float(np.median(arr))
    lo = med - k * iqr
    hi = med + k * iqr
    return np.clip(arr, lo, hi)


def _solve_d0_from_target(target: float) -> float:
    """Solve trigamma(d0/2) = target for d0. Returns inf when degenerate."""
    if target <= 0:
        return float("inf")

    def f(x):
        return float(_polygamma(1, x)) - target

    try:
        if f(1e-3) * f(1e6) > 0:
            return float("inf")
        x = _brentq(f, 1e-3, 1e6, xtol=1e-6)
    except Exception:
        return float("inf")
    return 2.0 * float(x)


def smyth_eb_prior(
    s2: np.ndarray, df: np.ndarray, robust: bool = False
) -> tuple[float, float]:
    """Estimate the pooled prior variance s0^2 and prior df d0 from a vector
    of per-cell-type residual variances and df (Smyth 2004 limma).

    Fits a scaled chi-square prior to the observed residual variances on
    the log scale via the trigamma moment-matching equation
        var(log s2) - mean(trigamma(df/2)) = trigamma(d0/2)

    When `robust=True`, the prior moments (e_mean, e_var) are estimated on
    the IQR-Winsorized e-vector — so a few extreme variances do not pull
    the prior (Phipson, Lee, Majewski, Alexander, Smyth 2016). Useful when
    one or two cell types have anomalously large/small residual variances
    that would otherwise inflate d0 toward infinity (over-shrinkage) or
    pull d0 toward zero (no shrinkage).

    Falls back to (mean(s2), inf) — i.e. infinite shrinkage to the pooled
    mean — when the moment estimator is degenerate (fewer than 3 cell
    types or near-uniform variances).
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
    e_for_moments = _winsorize(e) if robust else e
    e_mean = float(e_for_moments.mean())
    e_var = float(e_for_moments.var(ddof=1))
    target = e_var - float(_polygamma(1, dfv / 2.0).mean())
    d0 = _solve_d0_from_target(target)
    if not np.isfinite(d0):
        return float(np.exp(e_mean)), float("inf")
    x = d0 / 2.0
    log_s0_sq = e_mean + float(_digamma(x)) - float(np.log(x))
    return float(np.exp(log_s0_sq)), d0


def smyth_eb_prior_trended(
    s2: np.ndarray,
    df: np.ndarray,
    abundance: np.ndarray,
    lowess_frac: float = 0.5,
    robust: bool = False,
) -> tuple[np.ndarray, float]:
    """Trended Smyth prior (limma `eBayes(trend=TRUE)`, Smyth 2009): fits
    a smooth function f(abundance) to the per-cell-type log residual
    variance, so each cell type gets its own prior mean s0^2_c instead
    of pooling toward a single global mean.

    Use this when residual variances are systematically heterogeneous
    across cell types — e.g. when rare cell types live near the CLR eps
    floor and have structurally larger variance than dominant lineages.
    Pooled EB then over-shrinks the high-variance cell types (false
    positives) and under-shrinks the low-variance ones (false negatives);
    the trend absorbs that systematic structure.

    Implementation: e_c = log(s2_c) - digamma(df_c/2) + log(df_c/2) is
    the bias-corrected log-variance. We fit lowess(e ~ abundance) to get
    a smoothed mean; d0 is estimated from the residuals around the trend
    via the same trigamma identity. Per-cell-type s0^2_c is then backed
    out from the smoothed e-values using the d0 correction:
        log(s0^2_c) = smoothed_c - digamma(d0/2) + log(d0/2).

    `robust=True` Winsorizes the post-trend residuals before fitting d0,
    matching the rationale of the pooled robust path.

    Returns (s0_sq_per_celltype_array, d0_scalar). The array is aligned
    with the input `s2` order; positions where validation failed get
    s0^2 = NaN and the caller should fall back to the pooled value or
    skip moderation for that row.
    """
    from statsmodels.nonparametric.smoothers_lowess import lowess

    s2 = np.asarray(s2, dtype=float)
    df = np.asarray(df, dtype=float)
    abundance = np.asarray(abundance, dtype=float)
    valid = (
        (s2 > 0) & np.isfinite(s2) & np.isfinite(df) & (df > 0)
        & np.isfinite(abundance)
    )
    if valid.sum() < 5:
        # Too few cell types for a meaningful trend. Fall back to pooled.
        s0_pool, d0 = smyth_eb_prior(s2, df, robust=robust)
        return np.full(s2.shape, s0_pool, dtype=float), d0

    s2v, dfv, abv = s2[valid], df[valid], abundance[valid]
    z = np.log(s2v)
    e = z - _digamma(dfv / 2.0) + np.log(dfv / 2.0)
    # lowess returns smoothed values aligned with the input order when
    # return_sorted=False.
    smoothed = lowess(e, abv, frac=lowess_frac, return_sorted=False)
    smoothed = np.asarray(smoothed, dtype=float)
    if not np.all(np.isfinite(smoothed)):
        # Lowess failed (e.g. too few unique abundance values). Fall back.
        s0_pool, d0 = smyth_eb_prior(s2, df, robust=robust)
        return np.full(s2.shape, s0_pool, dtype=float), d0

    # Estimate d0 from the post-detrending residual variance, optionally
    # Winsorized for robustness against outlier cell types.
    resid = e - smoothed
    resid_for_moments = _winsorize(resid) if robust else resid
    target = float(resid_for_moments.var(ddof=1)) - float(
        _polygamma(1, dfv / 2.0).mean()
    )
    d0 = _solve_d0_from_target(target)

    # Back out per-cell-type s0^2 from the smoothed e-values.
    if np.isfinite(d0):
        x = d0 / 2.0
        log_s0_sq = smoothed - float(_digamma(x)) + float(np.log(x))
    else:
        # Infinite-prior limit: correction term -> 0, s0^2_c = exp(smoothed_c).
        log_s0_sq = smoothed

    s0_sq_full = np.full(s2.shape, np.nan, dtype=float)
    s0_sq_full[valid] = np.exp(log_s0_sq)
    return s0_sq_full, d0


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
    n_jobs: int = 1,
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

    Determinism: results are reproducible for a fixed combination of
    (--permutations, --permutation-seed, --threads). Changing --threads
    re-shards the B permutations across workers and therefore samples a
    different set of permuted label vectors — output p-values will agree
    statistically but not bit-identically across thread counts. Differences
    are bounded by Monte Carlo SE ≈ sqrt(p(1-p)/B).

    Returns a dict with three keys, each mapping cell_type -> p (omnibus)
    or cell_type -> {level: p} (pairwise):
      - "perm_p":       raw permutation p, P(stat_perm >= stat_obs)
      - "perm_p_wy":    Westfall-Young max-T FWER-adjusted p (joint null)
      - "perm_q_pooled": SAM-style FDR (Tusher Tibshirani Chu 2001) using
                         the pooled empirical null across cell types

    The three corrections trade off differently:
      - perm_p alone is per-test exact but ignores multiplicity
      - BH(perm_p) (computed by the caller) is the standard FDR control
        and is conservative under positive dependence
      - perm_p_wy controls family-wise error using the joint null and is
        more powerful than Bonferroni when cell types are correlated
      - perm_q_pooled assumes test statistics are exchangeable across
        cell types under the null (best when stats are on a common scale,
        e.g. EB-moderated; can be anti-conservative when raw t-stats have
        very heterogeneous scales)
    """
    # Build observed-statistic lookup from the already-fit results.
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

    # Canonical (cell_type, level) ordering for the per-permutation
    # statistic matrix that workers will fill column-wise. For omnibus,
    # level is None.
    pair_keys: list = []
    obs_stat_vec: list = []
    if test == "omnibus":
        for ct in cell_types:
            pair_keys.append((ct, None))
            obs_stat_vec.append(obs_stats.get(ct, float("nan")))
    else:
        for r in obs_results:
            ct = r["cell_type"]
            lvl = r.get("level")
            if lvl is None:
                continue
            pair_keys.append((ct, lvl))
            obs_stat_vec.append(
                obs_stats.get(ct, {}).get(lvl, float("nan"))
            )
    obs_stat_arr = np.asarray(obs_stat_vec, dtype=float)

    # Preserve the Categorical dtype (with the user's --reference-group
    # ordering) across permutations. A plain numpy assignment would coerce
    # to object dtype, after which Patsy reverts to alphabetical ordering
    # and the reference level changes — making the per-level counts diverge
    # from the observed-statistic keys and yielding NaN p-values.
    group_series = df[group_col]
    is_categorical = isinstance(group_series.dtype, pd.CategoricalDtype)
    group_categories = (
        list(group_series.cat.categories) if is_categorical else None
    )

    # Initialize merged-counts skeleton (also pre-seeds every observed
    # (cell_type, level) pair to 0 so cell types with zero hits still get
    # p = 1/(N+1) rather than a missing dict key).
    counts: dict = {}
    for ct in cell_types:
        counts[ct] = 0 if test == "omnibus" else {}
    if test != "omnibus":
        for res in obs_results:
            ct = res["cell_type"]
            lvl = res.get("level")
            if lvl is not None and ct in counts:
                counts[ct].setdefault(lvl, 0)

    # Shard permutations across workers. Each worker gets an independent
    # PCG64 stream seeded from `seed + i` so the union over chunks is a
    # reproducible sequence (same total seed → same merged counts) and
    # different chunks can't alias the same permutations.
    n_jobs = max(1, int(n_jobs))
    base_chunk = n_perm // n_jobs
    remainder = n_perm % n_jobs
    chunk_sizes = [base_chunk + (1 if i < remainder else 0) for i in range(n_jobs)]
    chunk_sizes = [c for c in chunk_sizes if c > 0]

    worker_args = [
        (df, group_col, covariates, test, weights, cell_types,
         obs_stats, group_categories, chunk_n, seed + 100003 * i,
         pair_keys)
        for i, chunk_n in enumerate(chunk_sizes)
    ]

    # Aggregator state: counts dict (BH-source p) + concatenated
    # per-permutation statistic matrices (for WY and SAM).
    stat_chunks: list = []
    max_chunks: list = []

    def _merge(partial: dict) -> None:
        for ct, val in partial["counts"].items():
            if test == "omnibus":
                counts[ct] = counts.get(ct, 0) + val
            else:
                bucket = counts.setdefault(ct, {})
                for lvl, c in val.items():
                    bucket[lvl] = bucket.get(lvl, 0) + c
        stat_chunks.append(partial["stat_matrix"])
        max_chunks.append(partial["max_stats"])

    if n_jobs == 1 or len(worker_args) == 1:
        # In-process path — no pickling overhead.
        for i, args in enumerate(worker_args):
            partial = _permutation_chunk_worker(args)
            _merge(partial)
            if verbose:
                print(
                    f"  permutation chunk {i + 1}/{len(worker_args)} "
                    f"({chunk_sizes[i]} perms) done",
                    file=sys.stderr,
                )
    else:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        if verbose:
            print(
                f"  dispatching {n_perm} permutations across {n_jobs} workers "
                f"({len(worker_args)} chunks)",
                file=sys.stderr,
            )
        with ProcessPoolExecutor(max_workers=n_jobs) as ex:
            futures = {ex.submit(_permutation_chunk_worker, a): i
                       for i, a in enumerate(worker_args)}
            done = 0
            for fut in as_completed(futures):
                partial = fut.result()
                _merge(partial)
                done += 1
                if verbose:
                    print(
                        f"  chunk {done}/{len(worker_args)} done",
                        file=sys.stderr,
                    )

    # ---- BH-source p (existing) -----------------------------------
    pvals: dict = {}
    for ct in cell_types:
        if test == "omnibus":
            pvals[ct] = (counts[ct] + 1) / (n_perm + 1)
        else:
            pvals[ct] = {
                lvl: (c + 1) / (n_perm + 1) for lvl, c in counts[ct].items()
            }

    # ---- Aggregate per-permutation statistic matrices --------------
    full_stat = (
        np.vstack(stat_chunks) if stat_chunks
        else np.zeros((0, len(pair_keys)), dtype=float)
    )
    full_max = (
        np.concatenate(max_chunks) if max_chunks
        else np.zeros(0, dtype=float)
    )

    # ---- Westfall-Young max-T FWER -----------------------------------
    # For each cell type c, p_WY(c) = P(max_c' |t_perm,c'| >= |t_obs,c|)
    # under the joint null. This is exact-FWER and tends to be much more
    # powerful than Bonferroni when cell types are correlated.
    valid_max = full_max[np.isfinite(full_max)]
    B_eff = int(valid_max.size)
    wy: dict = {}
    for i, (ct, lvl) in enumerate(pair_keys):
        t_o = obs_stat_arr[i]
        if not np.isfinite(t_o) or B_eff == 0:
            p_wy = float("nan")
        else:
            p_wy = (int(np.sum(valid_max >= t_o)) + 1) / (B_eff + 1)
        if test == "omnibus":
            wy[ct] = p_wy
        else:
            wy.setdefault(ct, {})[lvl] = p_wy

    # ---- SAM-style pooled-null FDR (Tusher Tibshirani Chu 2001) -----
    # Builds the empirical null by pooling permutation stats across both
    # permutations and cell types. For threshold tau,
    #   R(tau) = #{|t_obs| >= tau across cell types}
    #   V(tau) = (1/B) * #{|t_perm| >= tau across (perms, cell types)}
    #   FDR(tau) = V(tau) / R(tau), monotonized so smaller |t_obs| -> larger q.
    # Caveat: assumes test stats are exchangeable across cell types under
    # the null. Best when used with --empirical-bayes (so stats live on a
    # common moderated scale) or when residual variances are similar.
    sam_q_arr = _sam_pooled_qvalues(obs_stat_arr, full_stat)
    sam: dict = {}
    for i, (ct, lvl) in enumerate(pair_keys):
        if test == "omnibus":
            sam[ct] = float(sam_q_arr[i])
        else:
            sam.setdefault(ct, {})[lvl] = float(sam_q_arr[i])

    return {
        "perm_p": pvals,
        "perm_p_wy": wy,
        "perm_q_pooled": sam,
    }


def _sam_pooled_qvalues(
    obs: np.ndarray, perm_matrix: np.ndarray
) -> np.ndarray:
    """SAM-style FDR using the pooled-across-cell-types empirical null.

    obs: shape (K,) — observed test statistics (F or |t|; non-negative).
    perm_matrix: shape (B, K) — permutation stats with matched columns.
                 NaNs allowed for failed fits.
    Returns q-values shape (K,).
    """
    K = int(obs.size)
    out = np.full(K, float("nan"), dtype=float)
    if K == 0 or perm_matrix.size == 0:
        return out
    B = int(perm_matrix.shape[0])
    finite_perm = perm_matrix[np.isfinite(perm_matrix)]
    if finite_perm.size == 0:
        return out
    # Order cell types by observed statistic descending (largest first).
    order = np.argsort(-obs)
    for rank_idx, c in enumerate(order):
        tau = obs[c]
        if not np.isfinite(tau):
            continue
        # Observed positives at threshold tau.
        R = int(np.sum(obs[np.isfinite(obs)] >= tau))
        if R == 0:
            out[c] = 1.0
            continue
        # Expected number of null positives at tau, per permutation.
        V = float(np.sum(finite_perm >= tau)) / B
        out[c] = float(min(V / R, 1.0))
    # Monotonize: q at decreasing |t_obs| must be non-decreasing.
    last_q = float("nan")
    for c in order:
        if not np.isfinite(out[c]):
            continue
        if not np.isfinite(last_q):
            last_q = out[c]
            continue
        if out[c] < last_q:
            out[c] = last_q
        else:
            last_q = out[c]
    return out


def _permutation_chunk_worker(args: tuple) -> dict:
    """Module-level worker for ProcessPoolExecutor: run `n_chunk` permutations
    and return per-(cell_type, level) hit counts plus the full per-perm
    statistic matrix (used downstream for Westfall-Young and SAM FDR).

    Defined at module level (not as a closure) so it pickles cleanly under
    spawn-based multiprocessing on macOS / Windows.
    """
    (df, group_col, covariates, test, weights, cell_types,
     obs_stats, group_categories, n_chunk, seed, pair_keys) = args

    rng = np.random.default_rng(seed)
    is_categorical = group_categories is not None
    # Permute on the underlying labels; we re-wrap as Categorical below.
    group_vals = df[group_col].astype(object).to_numpy(copy=True)

    counts: dict = {ct: (0 if test == "omnibus" else {}) for ct in cell_types}
    if test != "omnibus":
        for ct, level_map in obs_stats.items():
            if ct in counts:
                for lvl in level_map:
                    counts[ct].setdefault(lvl, 0)

    # Statistic matrix layout: rows = permutations, cols = pair_keys order.
    n_pairs = len(pair_keys)
    pair_idx = {pk: i for i, pk in enumerate(pair_keys)}
    stat_matrix = np.full((n_chunk, n_pairs), float("nan"), dtype=float)
    max_stats = np.full(n_chunk, float("nan"), dtype=float)

    df_perm = df.copy()
    for b in range(n_chunk):
        perm_arr = rng.permutation(group_vals)
        if is_categorical:
            df_perm[group_col] = pd.Categorical(
                perm_arr, categories=group_categories
            )
        else:
            df_perm[group_col] = perm_arr
        running_max = float("-inf")
        for ct in cell_types:
            df_ct, w_ct = _resolve_celltype_weights(weights, df_perm, ct)
            if df_ct.shape[0] < 3:
                continue
            try:
                res = fit_one_celltype(
                    ct, df_ct, group_col, covariates, test, w_ct
                )
            except Exception:
                continue
            if test == "omnibus":
                f_perm = res.get("F_stat", float("nan"))
                f_obs = obs_stats.get(ct, float("nan"))
                if np.isfinite(f_perm):
                    idx = pair_idx.get((ct, None))
                    if idx is not None:
                        stat_matrix[b, idx] = f_perm
                    if f_perm > running_max:
                        running_max = f_perm
                    if np.isfinite(f_obs) and f_perm >= f_obs:
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
                    idx = pair_idx.get((ct, lvl))
                    if idx is not None:
                        stat_matrix[b, idx] = t_perm
                    if t_perm > running_max:
                        running_max = t_perm
                    t_obs = obs_stats.get(ct, {}).get(lvl, float("nan"))
                    if np.isfinite(t_obs) and t_perm >= t_obs:
                        counts[ct][lvl] = counts[ct].get(lvl, 0) + 1
        if running_max != float("-inf"):
            max_stats[b] = running_max

    return {
        "counts": counts,
        "stat_matrix": stat_matrix,
        "max_stats": max_stats,
    }


# -------- Non-parametric companion test (Mann-Whitney / Kruskal-Wallis) ----

def compute_nonparametric_pvalues(
    w: pd.DataFrame,
    metadata: pd.DataFrame,
    group_col: str,
    test: str,
    results: pd.DataFrame,
) -> pd.Series:
    """Per-cell-type two-sided Mann-Whitney p-values, aligned with the rows
    of `results`. Auto-extends to Kruskal-Wallis when the analysis has
    three or more groups in omnibus mode (i.e. when MW does not apply).

    The non-parametric test is rank-based and therefore invariant to the
    monotone CLR transform — testing the proportions directly gives
    identical p-values to testing the CLR-transformed values, so we
    compute on the raw proportions for clarity.

    Returns a Series of p-values indexed identically to `results`. NaN
    where the row's level / cell type isn't sample-bearing in both groups.
    """
    if isinstance(metadata[group_col].dtype, pd.CategoricalDtype):
        group_levels = list(metadata[group_col].cat.categories)
    else:
        group_levels = sorted(
            metadata[group_col].dropna().unique().tolist()
        )

    common = w.index.intersection(metadata.index)
    meta = metadata.loc[common]
    group_samples = {
        g: meta[meta[group_col] == g].index for g in group_levels
    }

    pvals: list = []
    for _, row in results.iterrows():
        ct = row.get("cell_type")
        if ct is None or ct not in w.columns:
            pvals.append(float("nan"))
            continue
        w_ct = w[ct]
        if test == "pairwise":
            lvl = row.get("level")
            if (
                lvl is None or lvl not in group_levels
                or lvl == group_levels[0]  # reference level itself
            ):
                pvals.append(float("nan"))
                continue
            ref_g = group_levels[0]
            v_ref = (
                w_ct.reindex(group_samples[ref_g]).dropna().to_numpy()
            )
            v_lvl = (
                w_ct.reindex(group_samples[lvl]).dropna().to_numpy()
            )
            if len(v_ref) < 1 or len(v_lvl) < 1:
                pvals.append(float("nan"))
                continue
            try:
                p = float(
                    _scipy_stats.mannwhitneyu(
                        v_ref, v_lvl, alternative="two-sided"
                    ).pvalue
                )
            except Exception:
                p = float("nan")
            pvals.append(p)
        else:
            # omnibus: K=2 -> Mann-Whitney, K>=3 -> Kruskal-Wallis
            arrs = [
                w_ct.reindex(group_samples[g]).dropna().to_numpy()
                for g in group_levels
            ]
            arrs = [a for a in arrs if len(a) >= 1]
            if len(arrs) < 2:
                pvals.append(float("nan"))
                continue
            try:
                if len(arrs) == 2:
                    p = float(
                        _scipy_stats.mannwhitneyu(
                            arrs[0], arrs[1], alternative="two-sided"
                        ).pvalue
                    )
                else:
                    p = float(_scipy_stats.kruskal(*arrs).pvalue)
            except Exception:
                p = float("nan")
            pvals.append(p)

    return pd.Series(pvals, index=results.index)


# -------- Per-cell-type weights and detection masking --------------------

def build_per_celltype_weights(
    deconv: pd.DataFrame,
    sample_index: pd.Index,
    cell_types: list,
    weight_mode: str,
    mask_mode: str,
    mask_source: str,
    mask_thresh: float,
    eps: float = 1e-9,
    verbose: bool = False,
) -> Optional[dict]:
    """Build a dict mapping cell_type -> per-sample weight Series, combining a
    CI-derived inverse-variance precision weight with an optional per-sample
    detection-based soft/hard mask.

    weight_mode:
      'none'         all-ones (precision factor disabled)
      'ci-width'     w = 1 / (CI_upper - CI_lower)^2 on the proportion scale
      'ci-width-rel' w = 1 / ((CI_upper - CI_lower) / proportion)^2  -- the
                     delta-method-correct GLS weight on the log/CLR scale,
                     recommended.

    mask_mode (uses the per-sample p_value or q_value column from
    BetaValueDeconvolution -bootstrap output as a detection signal):
      'none' no masking
      'soft' multiply weight by (1 - detection); samples where the cell
             type isn't reliably detected get smoothly down-weighted.
      'hard' set weight to 0 (effectively excluding the sample) when
             detection > mask_thresh. Costs df_resid for that cell type's
             regression.

    Returns None when both weight_mode and mask_mode are 'none' (signaling
    that no per-CT weight spec is needed). Otherwise returns a dict whose
    keys are cell_types and values are pd.Series indexed by sample.
    """
    if weight_mode == "none" and mask_mode == "none":
        return None

    prop = pivot_long_to_wide(deconv, "proportion").reindex(
        index=sample_index, columns=cell_types
    )

    if weight_mode == "ci-width":
        if "CI_lower" not in deconv.columns or "CI_upper" not in deconv.columns:
            raise ValueError(
                "--per-celltype-weights ci-width requires CI_lower / "
                "CI_upper. Run BetaValueDeconvolution with -bootstrap."
            )
        ci_lo = pivot_long_to_wide(deconv, "CI_lower").reindex(
            index=sample_index, columns=cell_types
        )
        ci_hi = pivot_long_to_wide(deconv, "CI_upper").reindex(
            index=sample_index, columns=cell_types
        )
        width = (ci_hi - ci_lo).clip(lower=eps)
        w_mat = 1.0 / (width ** 2)
    elif weight_mode == "ci-width-rel":
        if "CI_lower" not in deconv.columns or "CI_upper" not in deconv.columns:
            raise ValueError(
                "--per-celltype-weights ci-width-rel requires CI_lower / "
                "CI_upper. Run BetaValueDeconvolution with -bootstrap."
            )
        ci_lo = pivot_long_to_wide(deconv, "CI_lower").reindex(
            index=sample_index, columns=cell_types
        )
        ci_hi = pivot_long_to_wide(deconv, "CI_upper").reindex(
            index=sample_index, columns=cell_types
        )
        # Relative width on the proportion scale; floor proportion at eps
        # so a zero-proportion sample doesn't blow up.
        rel_width = ((ci_hi - ci_lo) / prop.clip(lower=1e-6)).clip(lower=eps)
        w_mat = 1.0 / (rel_width ** 2)
    else:  # 'none' precision factor: start from ones, only mask
        w_mat = pd.DataFrame(1.0, index=sample_index, columns=cell_types)

    # Detection mask
    if mask_mode != "none":
        if mask_source not in deconv.columns:
            raise ValueError(
                f"--per-celltype-mask-source {mask_source!r} not in deconv "
                f"input. Run BetaValueDeconvolution with -bootstrap and "
                f"-permutationTest to produce per-sample p_value/q_value."
            )
        det = pivot_long_to_wide(deconv, mask_source).reindex(
            index=sample_index, columns=cell_types
        )
        if mask_mode == "soft":
            # (1 - q) factor: samples where cell type is fully detected
            # (q ≈ 0) get full weight; non-detected (q ≈ 1) get zero.
            factor = (1.0 - det).clip(lower=0.0, upper=1.0)
            w_mat = w_mat * factor
        elif mask_mode == "hard":
            # Drop samples where detection > threshold by setting weight=0.
            keep = det <= mask_thresh
            w_mat = w_mat.where(keep, 0.0)

    # Replace NaN with 0 for safe handling in downstream WLS / hard-filter.
    w_mat = w_mat.fillna(0.0)

    out = {ct: w_mat[ct].astype(float) for ct in cell_types}

    if verbose:
        print(
            f"Per-CT weights (mode={weight_mode}, mask={mask_mode}"
            + (f", source={mask_source}" if mask_mode != "none" else "")
            + (f", thresh={mask_thresh}" if mask_mode == "hard" else "")
            + "):",
            file=sys.stderr,
        )
        # Diagnostic: count of nonzero samples per cell type
        for i, ct in enumerate(cell_types[:5]):
            ws = out[ct]
            nz = int((ws > 0).sum())
            wmin = float(ws[ws > 0].min()) if nz > 0 else 0.0
            wmax = float(ws.max())
            print(
                f"  {ct}: nonzero_samples={nz}/{len(ws)}, "
                f"weight_range=[{wmin:.3g}, {wmax:.3g}]",
                file=sys.stderr,
            )
        if len(cell_types) > 5:
            print(
                f"  ... and {len(cell_types) - 5} more cell types",
                file=sys.stderr,
            )

    return out


def _resolve_celltype_weights(
    weights_spec, df: pd.DataFrame, cell_type: str
) -> tuple:
    """Return (df_for_ct, weights_for_ct) for one cell type's regression.

    - weights_spec is None: no weighting -> (df unchanged, None).
    - weights_spec is a pd.Series: per-sample weights -> (df unchanged,
      aligned series).
    - weights_spec is a dict: per-cell-type weights. Look up the Series
      for `cell_type`, drop df rows whose weight is 0/NaN/<=0 (hard-mask
      semantics), and return (filtered df, aligned positive weights).
    """
    if weights_spec is None:
        return df, None
    if isinstance(weights_spec, dict):
        w_ct = weights_spec.get(cell_type)
        if w_ct is None:
            return df, None
        aligned = w_ct.reindex(df.index)
        mask = np.isfinite(aligned) & (aligned > 0)
        if not mask.all():
            return df.loc[mask], aligned.loc[mask]
        return df, aligned
    # pd.Series — per-sample
    return df, weights_spec


# -------- ALDEx2-style Monte Carlo over bootstrap CIs ---------------------

def _apply_eb_inline(
    rows: list, w: pd.DataFrame, test: str, eb_config: dict
) -> list:
    """Apply EB variance shrinkage in-place to a list of fit-row dicts.

    Mirrors the EB block in main() but as a pure function callable from
    workers. Recomputes p_value using moderated_t / moderated_F and
    annotates rows with eb_s0_sq, eb_d0, p_value_unshrunk.
    """
    if not rows:
        return rows
    table = pd.DataFrame(rows).drop_duplicates("cell_type").reset_index(drop=True)
    s2_arr = table["s2_resid"].to_numpy()
    df_arr = table["df_resid"].to_numpy()
    if eb_config.get("trend"):
        mean_prop = w.reindex(columns=table["cell_type"]).mean(axis=0)
        abundance = np.log(mean_prop.to_numpy() + 1e-8)
        s0_arr, d0 = smyth_eb_prior_trended(
            s2_arr, df_arr, abundance,
            lowess_frac=eb_config.get("trend_frac", 0.5),
            robust=eb_config.get("robust", False),
        )
        if np.isnan(s0_arr).any():
            pool, _ = smyth_eb_prior(
                s2_arr, df_arr, robust=eb_config.get("robust", False)
            )
            s0_arr = np.where(np.isnan(s0_arr), pool, s0_arr)
    else:
        pool, d0 = smyth_eb_prior(
            s2_arr, df_arr, robust=eb_config.get("robust", False)
        )
        s0_arr = np.full_like(s2_arr, pool, dtype=float)
    s0_lookup = dict(zip(table["cell_type"], s0_arr))
    for r in rows:
        s2 = r.get("s2_resid", float("nan"))
        dfr = r.get("df_resid", float("nan"))
        s0 = s0_lookup.get(r["cell_type"], float("nan"))
        r["eb_s0_sq"] = float(s0) if np.isfinite(s0) else float("nan")
        r["eb_d0"] = float(d0)
        r["p_value_unshrunk"] = r.get("p_value", float("nan"))
        if not np.isfinite(s0):
            r["p_value"] = float("nan")
            continue
        if test == "omnibus":
            _, p_mod = moderated_F_pvalue(
                r.get("F_stat", float("nan")),
                r.get("df_group", float("nan")),
                dfr, s2, d0, float(s0),
            )
        else:
            _, p_mod = moderated_t_pvalue(
                r.get("effect_clr", float("nan")),
                r.get("std_err", float("nan")),
                dfr, s2, d0, float(s0),
            )
        r["p_value"] = p_mod
    return rows


def _aldex_mc_worker(args: tuple) -> list:
    """Worker for ALDEx2 MC: redo CLR + per-cell-type OLS (+ optional EB) on
    one Monte Carlo replicate of the (samples x cell_types) proportion
    matrix. Module-level so it pickles cleanly under spawn-based MP.
    """
    (w_array, sample_idx, cell_types, metadata, group_col,
     covariates, test, weights_arr, group_categories,
     clr_eps, eb_config) = args

    w_rep = pd.DataFrame(w_array, index=sample_idx, columns=cell_types)
    w_clr_rep = clr_transform(w_rep, eps=clr_eps)

    meta = metadata.copy()
    if group_categories is not None and group_col in meta.columns:
        meta[group_col] = pd.Categorical(
            meta[group_col].astype(object), categories=group_categories
        )

    common = w_clr_rep.index.intersection(meta.index)
    df = w_clr_rep.loc[common].join(meta.loc[common])
    needed = [group_col] + covariates
    df = df.dropna(subset=needed)

    # Reconstitute weight spec. weights_arr can be:
    #   None                      -> no weighting
    #   np.ndarray (n_samples,)   -> per-sample WLS
    #   dict[cell_type -> array]  -> per-(sample, cell-type) WLS / mask
    weights_spec = None
    if isinstance(weights_arr, dict):
        weights_spec = {
            ct: pd.Series(arr, index=sample_idx)
            for ct, arr in weights_arr.items()
        }
    elif weights_arr is not None:
        weights_spec = pd.Series(weights_arr, index=sample_idx).reindex(df.index)

    rows: list = []
    for ct in cell_types:
        if ct not in df.columns:
            continue
        df_ct, w_ct = _resolve_celltype_weights(weights_spec, df, ct)
        if df_ct.shape[0] < 3:
            continue
        try:
            res = fit_one_celltype(
                ct, df_ct, group_col, covariates, test, w_ct
            )
        except Exception:
            continue
        if isinstance(res, dict):
            rows.append(res)
        else:
            rows.extend(res)

    if eb_config is not None and rows:
        try:
            rows = _apply_eb_inline(rows, w_rep, test, eb_config)
        except Exception:
            # If EB fails on a particular replicate (e.g. lowess failure),
            # fall through with raw OLS p-values rather than dropping the
            # replicate entirely.
            pass

    return rows


def run_aldex_mc(
    w_point: pd.DataFrame,
    ci_lower_df: pd.DataFrame,
    ci_upper_df: pd.DataFrame,
    metadata: pd.DataFrame,
    group_col: str,
    covariates: list,
    test: str,
    weights: Optional[pd.Series],
    group_categories: Optional[list],
    clr_eps: float,
    eb_config: Optional[dict],
    n_mc: int,
    seed: int,
    n_jobs: int,
    fdr_alpha: float,
    eps_floor: float = 1e-6,
    verbose: bool = False,
) -> dict:
    """Run M Monte Carlo replicates over the bootstrap CIs and aggregate.

    Sampling model: per (sample, cell_type), draw
        proportion_replicate ~ Normal(point, SD)
    where SD = (CI_upper - CI_lower) / 3.92 (Wald approximation to a 95%
    CI). Replicates are clipped to [eps_floor, inf) and renormalized to
    sum to 1 per sample, so they live on the simplex. This is the same
    Gaussian-on-CLR-implied-by-Wald-CI approximation that ALDEx2 uses
    internally when it lacks raw counts.

    Each replicate is passed through the same CLR + per-cell-type OLS/WLS
    + (optional) EB shrinkage pipeline as the point-estimate baseline.
    Permutation is NOT redone per replicate (would be B*M fits) — use
    the point-estimate permutation null instead.

    Returns a dict mapping aggregator name -> {(cell_type, level): value}:
        effect_median, effect_iqr, p_median, p_max, p_stability
    """
    rng = np.random.default_rng(seed)
    sample_idx = list(w_point.index)
    cell_types = list(w_point.columns)
    n_samples = len(sample_idx)
    n_cells = len(cell_types)

    ci_lo = ci_lower_df.reindex(index=sample_idx, columns=cell_types).to_numpy()
    ci_hi = ci_upper_df.reindex(index=sample_idx, columns=cell_types).to_numpy()
    sd_arr = (ci_hi - ci_lo) / 3.92
    sd_arr = np.where(np.isfinite(sd_arr) & (sd_arr > 0), sd_arr, 0.0)
    w_arr = w_point.to_numpy()

    # Serialize weights for the worker:
    #   pd.Series  -> numpy array (per-sample)
    #   dict       -> dict of numpy arrays (per-cell-type)
    if weights is None:
        weights_arr = None
    elif isinstance(weights, dict):
        weights_arr = {
            ct: ser.reindex(sample_idx).to_numpy()
            for ct, ser in weights.items()
        }
    else:
        weights_arr = weights.reindex(sample_idx).to_numpy()

    if verbose:
        nz_frac = float((sd_arr > 0).mean())
        print(
            f"  MC sampling: SD>0 in {nz_frac:.0%} of (sample, cell-type) "
            f"cells; mean SD={float(sd_arr[sd_arr>0].mean()) if nz_frac > 0 else 0:.3g}",
            file=sys.stderr,
        )

    # Pre-generate the M w-replicates upfront. Memory at M=128, n=10, K=25:
    # 128 * 250 * 8 = 256 KB. Trivial.
    rep_seeds = rng.integers(0, 2**31 - 1, size=n_mc)
    w_replicates: list = []
    for m in range(n_mc):
        rng_m = np.random.default_rng(int(rep_seeds[m]))
        noise = rng_m.normal(0.0, sd_arr, size=(n_samples, n_cells))
        rep = np.clip(w_arr + noise, eps_floor, None)
        rep = rep / rep.sum(axis=1, keepdims=True)
        w_replicates.append(rep)

    worker_args = [
        (w_replicates[m], sample_idx, cell_types, metadata, group_col,
         covariates, test, weights_arr, group_categories,
         clr_eps, eb_config)
        for m in range(n_mc)
    ]

    rep_results: list = [None] * n_mc
    n_jobs = max(1, int(n_jobs))
    if n_jobs == 1 or n_mc <= 1:
        for i, a in enumerate(worker_args):
            rep_results[i] = _aldex_mc_worker(a)
            if verbose and (i + 1) % max(1, n_mc // 10) == 0:
                print(f"  MC replicate {i + 1}/{n_mc} done", file=sys.stderr)
    else:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        if verbose:
            print(
                f"  dispatching {n_mc} MC replicates across {n_jobs} workers",
                file=sys.stderr,
            )
        with ProcessPoolExecutor(max_workers=n_jobs) as ex:
            futures = {
                ex.submit(_aldex_mc_worker, a): i
                for i, a in enumerate(worker_args)
            }
            count = 0
            for fut in as_completed(futures):
                i = futures[fut]
                try:
                    rep_results[i] = fut.result()
                except Exception as exc:
                    if verbose:
                        print(
                            f"  MC replicate {i} failed: {exc}",
                            file=sys.stderr,
                        )
                    rep_results[i] = []
                count += 1
                if verbose and count % max(1, n_mc // 10) == 0:
                    print(
                        f"  MC replicate {count}/{n_mc} done",
                        file=sys.stderr,
                    )

    # ---- Aggregate per-(cell_type, level) across replicates ---------
    keys: list = []
    effects: dict = {}
    pvals: dict = {}
    for rep in rep_results:
        if not rep:
            continue
        for r in rep:
            ct = r.get("cell_type")
            lvl = r.get("level") if test != "omnibus" else None
            key = (ct, lvl)
            if key not in effects:
                keys.append(key)
                effects[key] = []
                pvals[key] = []
            if test == "omnibus":
                eff = r.get("F_stat", float("nan"))
            else:
                eff = r.get("effect_clr", float("nan"))
            effects[key].append(float(eff) if eff is not None else float("nan"))
            p = r.get("p_value", float("nan"))
            pvals[key].append(float(p) if p is not None else float("nan"))

    eff_med, eff_iqr, p_med, p_max, p_stab = {}, {}, {}, {}, {}
    for key in keys:
        e_arr = np.asarray(effects[key], dtype=float)
        p_arr = np.asarray(pvals[key], dtype=float)
        e_fin = e_arr[np.isfinite(e_arr)]
        p_fin = p_arr[np.isfinite(p_arr)]
        if e_fin.size:
            eff_med[key] = float(np.median(e_fin))
            eff_iqr[key] = float(
                np.percentile(e_fin, 75) - np.percentile(e_fin, 25)
            )
        else:
            eff_med[key] = float("nan")
            eff_iqr[key] = float("nan")
        if p_fin.size:
            p_med[key] = float(np.median(p_fin))
            p_max[key] = float(np.max(p_fin))
            p_stab[key] = float((p_fin < fdr_alpha).mean())
        else:
            p_med[key] = float("nan")
            p_max[key] = float("nan")
            p_stab[key] = float("nan")

    return {
        "effect_median": eff_med,
        "effect_iqr": eff_iqr,
        "p_median": p_med,
        "p_max": p_max,
        "p_stability": p_stab,
    }


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

def plot_per_celltype_pdf(
    w: pd.DataFrame,
    metadata: pd.DataFrame,
    results: pd.DataFrame,
    group_col: str,
    output_pdf: str,
    plot_scale: str = "proportion",
    pvalue_source: str = "p_value",
    test: str = "pairwise",
    sort_by: str = "p_value",
    fdr_alpha: float = 0.05,
    verbose: bool = False,
) -> None:
    """Write one violin+jitter plot per cell type to a multi-page PDF.

    Each page: cell-type proportion (or log10-proportion / CLR) on the y-axis,
    one violin + jittered points per group on the x-axis. A comparison
    bracket with the user-selected p-value (3 significant figures) is drawn
    between groups; a dashed line connects the per-group medians.

    Cell types are ordered by `sort_by` (typically p_value, so the most
    interesting plots appear first in the PDF).
    """
    _, plt, PdfPages = _import_matplotlib()

    # Determine cell-type ordering (most-significant first if sort_by exists).
    if sort_by in results.columns:
        ordered_celltypes = (
            results.sort_values(sort_by, na_position="last")["cell_type"]
            .drop_duplicates()
            .tolist()
        )
    else:
        ordered_celltypes = list(w.columns)

    # Build p-value lookup from the results DataFrame.
    if pvalue_source not in results.columns:
        sys.exit(
            f"ERROR: --plot-pvalue-source {pvalue_source!r} not found in "
            f"results columns: {sorted(results.columns)}"
        )
    has_mw = "p_value_mw" in results.columns
    if test == "pairwise" and "level" in results.columns:
        p_lookup: dict = {}
        mw_lookup: dict = {}
        for _, r in results.iterrows():
            p_lookup.setdefault(r["cell_type"], {})[
                r.get("level")
            ] = r.get(pvalue_source, float("nan"))
            if has_mw:
                mw_lookup.setdefault(r["cell_type"], {})[
                    r.get("level")
                ] = r.get("p_value_mw", float("nan"))
    else:
        p_lookup = dict(zip(results["cell_type"], results[pvalue_source]))
        mw_lookup = (
            dict(zip(results["cell_type"], results["p_value_mw"]))
            if has_mw else {}
        )

    # Group ordering: respect the Categorical category order when set
    # (so the user's --reference-group becomes the leftmost violin).
    if isinstance(metadata[group_col].dtype, pd.CategoricalDtype):
        group_levels = list(metadata[group_col].cat.categories)
    else:
        group_levels = sorted(metadata[group_col].dropna().unique().tolist())

    common_idx = w.index.intersection(metadata.index)
    meta_aligned = metadata.loc[common_idx]

    # Pre-compute CLR matrix once if the user wants the CLR scale.
    if plot_scale == "clr":
        try:
            w_clr_full = clr_transform(w)
        except Exception:
            w_clr_full = None
    else:
        w_clr_full = None

    rng = np.random.default_rng(0)  # reproducible jitter across runs

    n_pages_written = 0
    Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(output_pdf) as pdf:
        for ct in ordered_celltypes:
            if ct not in w.columns:
                continue

            # y values per scale
            if plot_scale == "log-proportion":
                yvals_full = np.log10(w[ct].clip(lower=1e-7))
                ylabel = f"log10(proportion) — {ct}"
            elif plot_scale == "clr" and w_clr_full is not None:
                yvals_full = w_clr_full[ct]
                ylabel = f"CLR(proportion) — {ct}"
            else:
                yvals_full = w[ct]
                ylabel = f"proportion — {ct}"

            data_per_group = []
            n_per_group = []
            for g in group_levels:
                samples_g = meta_aligned[
                    meta_aligned[group_col] == g
                ].index
                vals_g = (
                    yvals_full.reindex(samples_g)
                    .dropna()
                    .to_numpy()
                )
                data_per_group.append(vals_g)
                n_per_group.append(len(vals_g))

            if not any(len(v) > 0 for v in data_per_group):
                continue  # nothing to plot for this cell type

            fig, ax = plt.subplots(figsize=(5.5, 4.5))
            positions = list(range(len(group_levels)))

            # Violin (only for groups with >=2 points; matplotlib refuses
            # otherwise). Sparse-group cell types still get jittered points.
            violin_data = [v for v in data_per_group if len(v) >= 2]
            violin_pos = [
                p for p, v in zip(positions, data_per_group) if len(v) >= 2
            ]
            if violin_data:
                parts = ax.violinplot(
                    violin_data,
                    positions=violin_pos,
                    showmeans=False,
                    showmedians=True,
                    showextrema=False,
                    widths=0.7,
                )
                for body in parts["bodies"]:
                    body.set_alpha(0.35)
                    body.set_edgecolor("black")
                    body.set_linewidth(0.5)
                if "cmedians" in parts:
                    parts["cmedians"].set_color("black")
                    parts["cmedians"].set_linewidth(1.0)

            # Jittered points
            for i, vals in enumerate(data_per_group):
                if len(vals) == 0:
                    continue
                xj = i + rng.normal(0, 0.05, size=len(vals))
                ax.scatter(
                    xj, vals, s=22, alpha=0.75,
                    edgecolor="black", linewidth=0.4, zorder=3,
                )

            # Dashed line connecting per-group medians
            medians = [
                float(np.median(v)) if len(v) else float("nan")
                for v in data_per_group
            ]
            finite = [
                (p, m) for p, m in zip(positions, medians)
                if np.isfinite(m)
            ]
            if len(finite) >= 2:
                xs, ys = zip(*finite)
                ax.plot(
                    xs, ys, color="black", linestyle="--",
                    linewidth=0.9, alpha=0.6, zorder=2,
                )

            # Compute y-range for placing the comparison brackets.
            all_vals = np.concatenate(
                [v for v in data_per_group if len(v) > 0]
            )
            y_min, y_max = float(all_vals.min()), float(all_vals.max())
            y_range = y_max - y_min if y_max > y_min else max(abs(y_max), 1.0)

            # Choose the non-parametric label: MW for K=2 group comparisons,
            # KW when omnibus on K>=3 groups (the only case where MW
            # doesn't apply directly).
            np_label = (
                "KW"
                if (test != "pairwise" and len(group_levels) >= 3)
                else "MW"
            )

            # Comparison brackets.
            if test == "pairwise" and len(group_levels) >= 2:
                # Pairwise: ref vs each non-ref level. One bracket per pair.
                p_dict = (
                    p_lookup.get(ct, {})
                    if isinstance(p_lookup.get(ct), dict)
                    else {None: p_lookup.get(ct, float("nan"))}
                )
                mw_dict = (
                    mw_lookup.get(ct, {})
                    if isinstance(mw_lookup.get(ct), dict)
                    else {None: mw_lookup.get(ct, float("nan"))}
                )
                bracket_y = y_max + 0.04 * y_range
                tick = 0.015 * y_range
                for j, g in enumerate(group_levels[1:], start=1):
                    p = p_dict.get(g, float("nan"))
                    if not np.isfinite(p):
                        continue
                    mw_p = mw_dict.get(g, float("nan")) if has_mw else float("nan")
                    ax.plot(
                        [0, 0, j, j],
                        [bracket_y, bracket_y + tick,
                         bracket_y + tick, bracket_y],
                        color="black", linewidth=0.8,
                    )
                    if np.isfinite(mw_p):
                        label = f"p = {p:.3g}\n{np_label} p = {mw_p:.3g}"
                    else:
                        label = f"p = {p:.3g}"
                    ax.text(
                        (0 + j) / 2.0,
                        bracket_y + tick * 1.4,
                        label,
                        ha="center", va="bottom", fontsize=9,
                    )
                    # Allocate enough vertical room for the two-line label.
                    bracket_y += y_range * (0.14 if np.isfinite(mw_p) else 0.10)
            else:
                # Omnibus: single p-value annotated at top.
                p = p_lookup.get(ct, float("nan"))
                mw_p = mw_lookup.get(ct, float("nan")) if has_mw else float("nan")
                if np.isfinite(p):
                    if np.isfinite(mw_p):
                        label = f"p = {p:.3g}\n{np_label} p = {mw_p:.3g}"
                    else:
                        label = f"p = {p:.3g}"
                    ax.text(
                        0.5, 0.96,
                        label,
                        transform=ax.transAxes,
                        ha="center", va="top", fontsize=10,
                        bbox=dict(boxstyle="round,pad=0.3",
                                  facecolor="white",
                                  edgecolor="lightgray"),
                    )

            ax.set_xticks(positions)
            ax.set_xticklabels(
                [f"{g}\n(n={n})"
                 for g, n in zip(group_levels, n_per_group)],
                rotation=0,
            )
            ax.set_ylabel(ylabel)
            ax.set_title(ct, fontsize=11)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            # Make room for the bracket annotations at the top.
            cur_lo, cur_hi = ax.get_ylim()
            ax.set_ylim(cur_lo, cur_hi + 0.05 * y_range)

            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
            n_pages_written += 1

    if verbose:
        print(
            f"Wrote {n_pages_written} cell-type plot(s) to {output_pdf}",
            file=sys.stderr,
        )


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
             "types (single weight per sample, applied to every cell-type "
             "regression). Requires CI columns. Overridden by "
             "--per-celltype-weights when set.",
    )
    p.add_argument(
        "--per-celltype-weights",
        choices=["none", "ci-width", "ci-width-rel"],
        default="none",
        help="Per-(sample, cell-type) inverse-variance weighting from the "
             "bootstrap CI width. Strictly more powerful than --weighted "
             "because deconvolution precision varies across cell types "
             "within a sample. "
             "ci-width: w = 1/(CI_upper - CI_lower)^2 on the proportion "
             "scale (simple). "
             "ci-width-rel: w = 1/((CI_upper-CI_lower)/proportion)^2 -- "
             "the delta-method-correct GLS weight on log/CLR scale "
             "(recommended). Overrides --weighted. Default: none.",
    )
    p.add_argument(
        "--per-celltype-mask-by-detection",
        choices=["none", "soft", "hard"],
        default="none",
        help="Mask per-(sample, cell-type) data points using the per-sample "
             "permutation p_value / q_value from BetaValueDeconvolution as "
             "a detection signal -- complementary to CI-based precision "
             "weighting. "
             "soft: multiply per-CT weight by (1 - detection); samples "
             "where the cell type isn't reliably detected get smoothly "
             "down-weighted, df_resid preserved. "
             "hard: drop samples where detection > "
             "--per-celltype-mask-thresh; clean exclusion but costs "
             "df_resid. Default: none.",
    )
    p.add_argument(
        "--per-celltype-mask-source",
        choices=["p_value", "q_value"],
        default="q_value",
        help="Which detection column to use for "
             "--per-celltype-mask-by-detection. q_value is "
             "BH-corrected within-sample and is safer for hard masking; "
             "p_value is raw and is more sensitive for soft masking. "
             "Default: q_value.",
    )
    p.add_argument(
        "--per-celltype-mask-thresh",
        type=float,
        default=0.5,
        help="Detection threshold for --per-celltype-mask-by-detection "
             "hard. Samples with detection > thresh are dropped from that "
             "cell type's regression. Default: 0.5.",
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
        "--eb-trend",
        action="store_true",
        help="Trended EB (limma 'trend=TRUE', Smyth 2009): instead of "
             "pooling toward a single global prior mean, fit a lowess "
             "smoother of log(s^2) vs log(mean proportion) and shrink each "
             "cell type toward its trended prior. Strongly recommended when "
             "rare and dominant cell types coexist (rare ones have "
             "structurally larger CLR variance, so a pooled prior over- "
             "shrinks them). No effect without --empirical-bayes.",
    )
    p.add_argument(
        "--eb-trend-frac",
        type=float,
        default=0.5,
        help="lowess span for --eb-trend (fraction of cell types in each "
             "local fit). Default 0.5; raise toward 1.0 for very smooth "
             "trends with few cell types. No effect without --eb-trend.",
    )
    p.add_argument(
        "--eb-robust",
        action="store_true",
        help="Robust prior fitting (limma 'robust=TRUE', Phipson et al. "
             "2016): IQR-Winsorize the e-vector (or post-trend residuals "
             "when combined with --eb-trend) before estimating d0, so a "
             "few cell types with anomalously large/small variance do not "
             "drag the prior. No effect without --empirical-bayes.",
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
        "--aldex-mc",
        type=int,
        default=0,
        help="If >0, run ALDEx2-style Monte Carlo aggregation: sample M "
             "replicates of the (samples x cell-types) proportion matrix "
             "from the bootstrap CIs of BetaValueDeconvolution, redo the "
             "regression on each, and aggregate effect/p-value across "
             "replicates. Adds columns effect_clr_mc_median, "
             "effect_clr_mc_iqr, p_value_mc (= median p across replicates, "
             "the 'expected p'), p_value_mc_max, p_value_mc_stability "
             "(fraction of replicates with p<fdr_alpha), q_value_mc. "
             "Composes with --empirical-bayes (each replicate gets its own "
             "EB prior). Recommended M=128. Requires CI_lower/CI_upper in "
             "the deconv input. Default: 0 (off).",
    )
    p.add_argument(
        "--aldex-mc-seed",
        type=int,
        default=0,
        help="RNG seed for --aldex-mc (default: 0, deterministic).",
    )
    p.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Parallelize the permutation null across N processes "
             "(ProcessPoolExecutor). Each worker runs an independent shard "
             "of the B permutations and counts are merged at the end. "
             "Speedup is roughly linear up to physical core count. "
             "Only used when --permutations > 0. Default: 1 (single-process).",
    )
    p.add_argument(
        "--output",
        required=True,
        help="Output TSV path. One row per cell type (omnibus) or per "
             "(cell_type, level) pair (pairwise).",
    )
    p.add_argument(
        "--plot-pdf",
        default=None,
        help="If set, write a multi-page PDF with one violin + jittered-"
             "points plot per cell type comparing the groups. The chosen "
             "p-value (see --plot-pvalue-source) is annotated as a "
             "comparison bracket between groups, formatted to 3 "
             "significant figures. Cell types are ordered most-significant "
             "first. Requires matplotlib.",
    )
    p.add_argument(
        "--plot-scale",
        choices=["proportion", "log-proportion", "clr"],
        default="proportion",
        help="Y-axis scale for --plot-pdf. proportion: raw 0-1 fraction. "
             "log-proportion: log10(proportion). clr: centered log-ratio "
             "(matches the regression scale). Default: proportion.",
    )
    p.add_argument(
        "--plot-pvalue-source",
        default="p_value",
        help="Which column from the results to display on the plots "
             "(e.g. p_value, p_value_unshrunk, p_value_perm, p_value_mc). "
             "Default: p_value.",
    )
    p.add_argument(
        "--plot-show-refined",
        action="store_true",
        help="Plot the refined proportion matrix (the values the regression "
             "actually sees) rather than the raw deconvolution output. "
             "Only meaningful with --refine-zeros-with-ci, which silently "
             "replaces zero proportions with CI_upper/2 for the regression. "
             "Default behavior is to plot raw values so the violins reflect "
             "what BetaValueDeconvolution actually reported.",
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
    # Keep an unmodified copy of the raw proportion matrix. The plot path
    # uses this so per-cell-type violins reflect what the deconvolution
    # actually reported, not the refined values that the regression uses
    # to avoid log(0) issues in CLR. (Opt back into refined view via
    # --plot-show-refined if you want to see what the regression sees.)
    w_raw = w.copy()

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
            if args.verbose:
                # Count how many (sample, cell-type) cells changed.
                changed_mask = (w_raw == 0) & (w > 0)
                n_changed = int(changed_mask.sum().sum())
                n_zero_total = int((w_raw == 0).sum().sum())
                # Per-cell-type breakdown for cell types where the change
                # is large (>= 25% of samples).
                per_ct = changed_mask.sum(axis=0)
                heavy = per_ct[per_ct >= max(1, int(0.25 * len(w_raw)))].sort_values(ascending=False)
                print(
                    f"--refine-zeros-with-ci: refined {n_changed:,} of "
                    f"{n_zero_total:,} zero proportion cells to CI_upper/2 "
                    f"(threshold: {args.refine_zeros_ci_min:g}).",
                    file=sys.stderr,
                )
                if len(heavy) > 0:
                    print(
                        f"  Cell types with >=25% of samples refined: "
                        f"{', '.join(f'{ct}({int(n)})' for ct, n in heavy.items())}",
                        file=sys.stderr,
                    )

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

    covariates = [c.strip() for c in args.covariates.split(",") if c.strip()]
    cell_types = list(w_clr.columns)

    # Compute weight specification.
    # Three mutually exclusive modes:
    #   weights_spec = None        -> OLS per cell type
    #   weights_spec = pd.Series   -> per-sample WLS (current --weighted)
    #   weights_spec = dict        -> per-(sample, cell-type) WLS with
    #                                  optional detection mask
    weights_spec = None
    use_per_ct = (
        args.per_celltype_weights != "none"
        or args.per_celltype_mask_by_detection != "none"
    )
    if use_per_ct:
        if args.weighted and args.verbose:
            print(
                "Note: --per-celltype-weights / --per-celltype-mask-by-"
                "detection set; --weighted is overridden.",
                file=sys.stderr,
            )
        try:
            weights_spec = build_per_celltype_weights(
                deconv=deconv,
                sample_index=df.index,
                cell_types=cell_types,
                weight_mode=args.per_celltype_weights,
                mask_mode=args.per_celltype_mask_by_detection,
                mask_source=args.per_celltype_mask_source,
                mask_thresh=args.per_celltype_mask_thresh,
                verbose=args.verbose,
            )
        except ValueError as exc:
            sys.exit(f"ERROR: {exc}")
    elif args.weighted:
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
            weights_spec = 1.0 / (mean_width.pow(2) + 1e-9)

    # `weights` retained as the legacy variable name throughout the rest
    # of main() (some downstream call sites reuse it directly).
    weights = weights_spec

    # ---- Fit per cell type --------------------------------------------
    rows: list[dict] = []
    for ct in cell_types:
        df_ct, w_ct = _resolve_celltype_weights(weights_spec, df, ct)
        if df_ct.shape[0] < 3:
            # Hard mask removed nearly all samples for this cell type.
            if args.verbose:
                print(
                    f"WARNING: per-CT mask left {df_ct.shape[0]} samples "
                    f"for {ct}; skipping.",
                    file=sys.stderr,
                )
            res = {
                "cell_type": ct,
                "p_value": np.nan,
                "n_samples": int(df_ct.shape[0]),
            }
        else:
            try:
                res = fit_one_celltype(
                    ct, df_ct, args.group_col, covariates, args.test, w_ct
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

    # ---- Non-parametric companion test (Mann-Whitney / Kruskal-Wallis)
    # Computed on the RAW (unrefined) proportion matrix, so it serves as
    # a refinement-bias detector: if --refine-zeros-with-ci is doing the
    # heavy lifting (replacing many zeros with CI_upper/2), the
    # parametric p_value will diverge from p_value_mw because MW sees
    # mostly-zero ranks. Convergence -> robust signal; divergence ->
    # the parametric significance leans on the imputed values.
    # Two-sided MW for binary group comparisons; auto-extends to
    # Kruskal-Wallis when omnibus mode has K>=3 groups.
    w_raw_filtered = w_raw.reindex(index=df.index, columns=w.columns)
    results["p_value_mw"] = compute_nonparametric_pvalues(
        w=w_raw_filtered,
        metadata=metadata.loc[df.index],
        group_col=args.group_col,
        test=args.test,
        results=results,
    )
    results["q_value_mw"] = bh_correct(results["p_value_mw"].to_numpy())

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
            .reset_index(drop=True)
        )
        s2_arr = prior_table["s2_resid"].to_numpy()
        df_arr = prior_table["df_resid"].to_numpy()

        if args.eb_trend:
            # Per-cell-type abundance covariate for the variance trend:
            # log mean proportion across samples (post-prevalence-filter).
            # Add a small floor to avoid -inf when a cell type is at the
            # eps boundary in every sample.
            mean_prop = w.reindex(columns=prior_table["cell_type"]).mean(
                axis=0
            )
            abundance = np.log(mean_prop.to_numpy() + 1e-8)
            s0_sq_arr, d0 = smyth_eb_prior_trended(
                s2_arr, df_arr, abundance,
                lowess_frac=args.eb_trend_frac,
                robust=args.eb_robust,
            )
            # Fall back to pooled value for any NaN positions (e.g. the
            # cell type was filtered as invalid by the trend estimator).
            if np.isnan(s0_sq_arr).any():
                pool_s0, _ = smyth_eb_prior(
                    s2_arr, df_arr, robust=args.eb_robust
                )
                s0_sq_arr = np.where(
                    np.isnan(s0_sq_arr), pool_s0, s0_sq_arr
                )
            mode_label = (
                f"trended (frac={args.eb_trend_frac:.2g})"
                + (", robust" if args.eb_robust else "")
            )
            if args.verbose:
                lo, hi = float(np.min(s0_sq_arr)), float(np.max(s0_sq_arr))
                print(
                    f"EB prior {mode_label}: s0^2 range=[{lo:.4g}, {hi:.4g}], "
                    f"d0={d0:.4g}, n_celltypes={len(prior_table)}",
                    file=sys.stderr,
                )
        else:
            s0_pool, d0 = smyth_eb_prior(
                s2_arr, df_arr, robust=args.eb_robust
            )
            s0_sq_arr = np.full_like(s2_arr, s0_pool, dtype=float)
            mode_label = "pooled" + (", robust" if args.eb_robust else "")
            if args.verbose:
                print(
                    f"EB prior {mode_label}: s0^2={s0_pool:.4g}, d0={d0:.4g}, "
                    f"n_celltypes={len(prior_table)}",
                    file=sys.stderr,
                )

        # Map cell_type -> per-cell-type s0^2 so pairwise rows (multiple
        # rows per cell type) all see the same value.
        s0_lookup = dict(zip(prior_table["cell_type"], s0_sq_arr))

        results["p_value_unshrunk"] = results["p_value"].copy()
        new_p: list[float] = []
        eb_s0_col: list[float] = []
        for _, row in results.iterrows():
            s2 = row["s2_resid"]
            df_r = row["df_resid"]
            s0_sq_c = s0_lookup.get(row["cell_type"], float("nan"))
            eb_s0_col.append(s0_sq_c)
            if not np.isfinite(s0_sq_c):
                new_p.append(float("nan"))
                continue
            if args.test == "omnibus":
                _, p_mod = moderated_F_pvalue(
                    row.get("F_stat", float("nan")),
                    row.get("df_group", float("nan")),
                    df_r, s2, d0, s0_sq_c,
                )
            else:
                _, p_mod = moderated_t_pvalue(
                    row.get("effect_clr", float("nan")),
                    row.get("std_err", float("nan")),
                    df_r, s2, d0, s0_sq_c,
                )
            new_p.append(p_mod)
        results["p_value"] = new_p
        results["eb_s0_sq"] = eb_s0_col
        results["eb_d0"] = d0

    # ---- ALDEx2-style Monte Carlo over bootstrap CIs -----------------
    # Re-runs the CLR + per-cell-type OLS (+ optional EB) pipeline on M
    # Monte Carlo replicates of the (samples x cell-types) proportion
    # matrix sampled from the bootstrap CIs. The "expected p-value"
    # (median across replicates) folds the per-(sample, cell-type)
    # deconvolution uncertainty into the inference, and tends to be
    # smaller for cell types whose effect is consistent across replicates
    # while being larger for those whose effect is fragile.
    if args.aldex_mc > 0:
        if "CI_lower" not in deconv.columns or "CI_upper" not in deconv.columns:
            sys.exit(
                "ERROR: --aldex-mc requires CI_lower and CI_upper columns. "
                "Run BetaValueDeconvolution with -bootstrap to produce them."
            )
        if args.verbose:
            print(
                f"Running ALDEx2 Monte Carlo with {args.aldex_mc} replicates "
                f"on {len(cell_types)} cell types...",
                file=sys.stderr,
            )
        ci_lo_full = pivot_long_to_wide(deconv, "CI_lower")
        ci_hi_full = pivot_long_to_wide(deconv, "CI_upper")
        # Align to the post-prevalence-filter w (which the baseline used),
        # restricted to samples that survived the metadata join.
        common_idx = df.index
        w_for_mc = w.reindex(index=common_idx, columns=cell_types)
        ci_lo_for_mc = ci_lo_full.reindex(
            index=common_idx, columns=cell_types
        ).fillna(0.0)
        ci_hi_for_mc = ci_hi_full.reindex(
            index=common_idx, columns=cell_types
        ).fillna(0.0)

        eb_cfg = None
        if args.empirical_bayes:
            eb_cfg = {
                "trend": args.eb_trend,
                "trend_frac": args.eb_trend_frac,
                "robust": args.eb_robust,
            }

        meta_for_mc = metadata.loc[common_idx]
        group_categories_mc = (
            list(meta_for_mc[args.group_col].cat.categories)
            if isinstance(
                meta_for_mc[args.group_col].dtype, pd.CategoricalDtype
            )
            else None
        )

        mc_summary = run_aldex_mc(
            w_point=w_for_mc,
            ci_lower_df=ci_lo_for_mc,
            ci_upper_df=ci_hi_for_mc,
            metadata=meta_for_mc,
            group_col=args.group_col,
            covariates=covariates,
            test=args.test,
            weights=weights,
            group_categories=group_categories_mc,
            clr_eps=args.clr_eps,
            eb_config=eb_cfg,
            n_mc=args.aldex_mc,
            seed=args.aldex_mc_seed,
            n_jobs=args.threads,
            fdr_alpha=args.fdr_alpha,
            verbose=args.verbose,
        )

        def _mc_key(r):
            return (
                r["cell_type"],
                r.get("level") if args.test != "omnibus" else None,
            )

        results["effect_clr_mc_median"] = results.apply(
            lambda r: mc_summary["effect_median"].get(
                _mc_key(r), float("nan")
            ),
            axis=1,
        )
        results["effect_clr_mc_iqr"] = results.apply(
            lambda r: mc_summary["effect_iqr"].get(
                _mc_key(r), float("nan")
            ),
            axis=1,
        )
        results["p_value_mc"] = results.apply(
            lambda r: mc_summary["p_median"].get(_mc_key(r), float("nan")),
            axis=1,
        )
        results["p_value_mc_max"] = results.apply(
            lambda r: mc_summary["p_max"].get(_mc_key(r), float("nan")),
            axis=1,
        )
        results["p_value_mc_stability"] = results.apply(
            lambda r: mc_summary["p_stability"].get(
                _mc_key(r), float("nan")
            ),
            axis=1,
        )
        results["q_value_mc"] = bh_correct(
            results["p_value_mc"].to_numpy()
        )

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
            n_jobs=args.threads,
        )
        def _map_perm_dict(d: dict) -> "pd.Series":
            if args.test == "omnibus":
                return results["cell_type"].map(d)
            return results.apply(
                lambda r: d.get(r["cell_type"], {}).get(
                    r.get("level"), float("nan")
                ),
                axis=1,
            )

        results["p_value_perm"] = _map_perm_dict(perm_p["perm_p"])
        results["p_value_perm_wy"] = _map_perm_dict(perm_p["perm_p_wy"])
        results["q_value_perm_pooled"] = _map_perm_dict(
            perm_p["perm_q_pooled"]
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

    # ---- Per-cell-type violin/jitter PDF -----------------------------
    if args.plot_pdf is not None:
        # By default, plot the RAW proportion matrix (what BetaValueDeconvolution
        # actually reported), not the post-refinement values the regression sees.
        # This is so the violins are an honest diagnostic of the input.
        # Use --plot-show-refined to switch back to the refined view.
        w_for_plot = w if args.plot_show_refined else w_raw
        plot_data_label = (
            "refined-for-regression" if args.plot_show_refined else "raw"
        )
        if args.verbose:
            print(
                f"Writing per-cell-type plots to {args.plot_pdf} "
                f"(scale={args.plot_scale}, "
                f"data={plot_data_label}, "
                f"p-value source={args.plot_pvalue_source})...",
                file=sys.stderr,
            )
        plot_per_celltype_pdf(
            w=w_for_plot.reindex(df.index),
            metadata=metadata.loc[df.index],
            results=results,
            group_col=args.group_col,
            output_pdf=args.plot_pdf,
            plot_scale=args.plot_scale,
            pvalue_source=args.plot_pvalue_source,
            test=args.test,
            sort_by="q_value" if "q_value" in results.columns else "p_value",
            fdr_alpha=args.fdr_alpha,
            verbose=args.verbose,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
