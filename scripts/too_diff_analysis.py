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

    base = {
        "cell_type": cell_type,
        "n_samples": int(fit.nobs),
        "r_squared": float(fit.rsquared),
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
            else:
                base["p_value"] = float(anova.loc[row_label, "PR(>F)"])
                base["F_stat"] = float(anova.loc[row_label, "F"])
        except Exception as exc:  # pragma: no cover
            print(f"[{cell_type}] anova failed: {exc}", file=sys.stderr)
            base["p_value"] = np.nan
            base["F_stat"] = np.nan

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
