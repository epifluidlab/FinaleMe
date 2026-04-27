#!/usr/bin/env python3
"""too_diff_power_simulation.py — Power simulation for too_diff_analysis.py.

Estimates the statistical power of the cohort-level differential cell-type
analysis (CLR + linear model + BH; see scripts/too_diff_analysis.py) across
a 3D grid of:

  - Sample sizes (n_pairs of cases vs controls)
  - Sequencing coverage (fragments per sample, after QC filters)
  - Effect sizes (proportion-scale shift in the target cell type)

What "power" means here
-----------------------
For each (effect_size, n_pairs, n_fragments) cell, we run --n-iter replicate
cohorts. Power is the fraction of replicates where the *target* cell type's
BH-corrected q-value <= --fdr-alpha.

The companion column `any_other_flagged_rate` is the fraction of replicates
where any *non-target* cell type was also flagged significant. Its
interpretation depends on the effect size:
  - effect_size == 0 (true global null): this is the family-wise false
    positive rate. Should be near --fdr-alpha when BH is well-calibrated.
  - effect_size > 0: a mix of (a) genuine compositional consequences
    (when one component of a composition shifts hard, the others
    *necessarily* shift in the opposite direction; CLR correctly detects
    this) and (b) noise-driven spurious co-detections. Not a "Type I"
    error in the classical sense.

Noise model
-----------
Per-sample observed proportions combine two independent noise sources:

1. Biological between-sample variability (coverage-independent).
   Drawn from Dirichlet(baseline * --bio-alpha). Higher alpha = tighter
   clustering. Default --bio-alpha=200 gives a CV of roughly 8% at
   baseline=0.4 and 22% at baseline=0.1, matching healthy-plasma atlas
   variability.

2. Technical measurement noise from BetaValueDeconvolution
   (coverage-dependent). Modeled as Gaussian on the proportion scale with
   SD = --tech-sd-anchor * sqrt(--tech-anchor-n / n_fragments). Defaults
   give SD=0.005 at 100M fragments and SD=0.05 at 1M, matching bootstrap
   CIs reported in Liu et al. 2024 (Nat Commun). After noise the
   proportions are clipped to (eps, 1) and renormalized to sum to 1.

Statistical test
----------------
Mirrors too_diff_analysis.py exactly when --test pairwise --reference-group
Control with no covariates: per-cell-type Welch t-test on CLR-transformed
values (algebraically equivalent to OLS for the disease-vs-control coefficient),
then BH-correct p-values across cell types. With covariates, expect 10-20%
less power per added independent covariate (rule of thumb), so this script
gives a slightly optimistic upper bound for covariate-adjusted analyses.

Output
------
A tidy TSV with one row per (effect_size, n_pairs, n_fragments) cell and
columns: power, any_other_flagged_rate, expected_other_per_iter,
median_q_target, plus the simulation parameters as constants. Optionally:
  - --mde-output: a (n_pairs x n_fragments) summary of the smallest effect
    size with power >= --target-power.
  - --plot: heatmap PNGs (one panel per effect_size).

Example
-------
    python scripts/too_diff_power_simulation.py \\
        --n-pairs 2,5,10,20,50,100 \\
        --n-fragments 1e6,5e6,1e7,2e7,5e7,1e8 \\
        --effect-sizes 0.005,0.01,0.02,0.05,0.10,0.20 \\
        --target-celltype Hepatocyte \\
        --n-iter 500 \\
        --output results/too_diff_power.tsv \\
        --mde-output results/too_diff_mde.tsv \\
        --plot results/too_diff_power.png \\
        --seed 42 --verbose
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as _scipy_stats


# --- Defaults -------------------------------------------------------------

# Realistic healthy-plasma cell-type composition. Sums to 1 after
# normalize. Atlas-like, 13 cell types covering blood + common solid
# organs. Adjust via --baseline-tsv to use a different reference.
DEFAULT_BASELINE = {
    "Blood-T":          0.30,
    "Blood-B":          0.04,
    "Blood-NK":         0.02,
    "Blood-Mono+Macro": 0.10,
    "Blood-Granul":     0.40,
    "Endothelium":      0.03,
    "Hepatocyte":       0.03,
    "Lung":             0.01,
    "Colon":            0.01,
    "Pancreas":         0.005,
    "Prostate":         0.005,
    "Breast":           0.005,
    "Other":            0.045,
}

DEFAULT_N_PAIRS_STR = "2,5,10,20,50,100"
DEFAULT_N_FRAGMENTS_STR = "1e6,5e6,1e7,2e7,5e7,1e8"
DEFAULT_EFFECT_SIZES_STR = "0.005,0.01,0.02,0.05,0.10,0.20"

DEFAULT_TECH_SD_ANCHOR = 0.005     # SD of proportion at the anchor coverage
DEFAULT_TECH_ANCHOR_N = 1e8        # 100M fragments
DEFAULT_BIO_ALPHA = 200.0
DEFAULT_FDR_ALPHA = 0.05
DEFAULT_N_ITER = 500
DEFAULT_TARGET_CELLTYPE = "Hepatocyte"
DEFAULT_SEED = 42


# --- Reuse BH from too_diff_analysis.py if importable ----------------------

_THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(_THIS_DIR))
try:
    # too_diff_analysis.py imports statsmodels at module load. The
    # power simulation only needs bh_correct; tolerate missing statsmodels.
    from too_diff_analysis import bh_correct as _ext_bh_correct
except Exception:
    _ext_bh_correct = None
finally:
    if str(_THIS_DIR) in sys.path:
        sys.path.remove(str(_THIS_DIR))


def bh_correct(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values (delegates to too_diff_analysis if available)."""
    if _ext_bh_correct is not None:
        return _ext_bh_correct(p)
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
    q_min = np.minimum.accumulate(q[::-1])[::-1]
    q_min = np.clip(q_min, 0.0, 1.0)
    q_orig = np.empty_like(q_min)
    q_orig[order] = q_min
    out[valid] = q_orig
    return out


# --- Compositional core ---------------------------------------------------

def normalized_baseline(baseline: dict[str, float]) -> tuple[np.ndarray, list[str]]:
    cell_types = list(baseline.keys())
    p = np.array([baseline[c] for c in cell_types], dtype=float)
    p = p / p.sum()
    return p, cell_types


def load_baseline_tsv(path: str) -> dict[str, float]:
    """Load a 2-column TSV: cell_type, baseline_proportion."""
    df = pd.read_csv(path, sep="\t", comment="#")
    if df.shape[1] < 2:
        raise ValueError(f"--baseline-tsv {path} must have >= 2 columns")
    name_col, prop_col = df.columns[:2]
    out = {}
    for _, row in df.iterrows():
        out[str(row[name_col])] = float(row[prop_col])
    return out


def apply_effect(baseline: np.ndarray, target_idx: int, delta: float) -> np.ndarray:
    """Shift target cell type by +delta on the proportion scale.

    Other cell types are scaled proportionally so the result still sums to 1.
    delta can be negative (depletion).
    """
    out = baseline.copy()
    new_target = baseline[target_idx] + delta
    new_target = float(np.clip(new_target, 1e-4, 0.999))
    out[target_idx] = new_target
    others = np.ones(len(out), dtype=bool)
    others[target_idx] = False
    others_sum = baseline[others].sum()
    if others_sum > 0:
        out[others] = baseline[others] * ((1.0 - new_target) / others_sum)
    return out


def sample_bio(baseline: np.ndarray, n: int, alpha: float,
               rng: np.random.Generator) -> np.ndarray:
    """Sample n biologically variable per-sample proportions from Dirichlet."""
    return rng.dirichlet(baseline * alpha, size=n)


def add_tech_noise(props: np.ndarray, n_fragments: float,
                   sd_anchor: float, anchor_n: float,
                   rng: np.random.Generator, eps: float = 1e-6) -> np.ndarray:
    """Add Gaussian noise on the proportion scale; renormalize.

    SD scales as sd_anchor * sqrt(anchor_n / n_fragments).
    """
    sd = sd_anchor * math.sqrt(anchor_n / n_fragments)
    if sd <= 0:
        return props.copy()
    noise = rng.normal(0, sd, props.shape)
    noisy = props + noise
    noisy = np.clip(noisy, eps, None)
    noisy = noisy / noisy.sum(axis=-1, keepdims=True)
    return noisy


def simulate_cohort(n_pairs: int, n_fragments: float,
                    baseline_ctrl: np.ndarray, baseline_dis: np.ndarray,
                    bio_alpha: float, sd_anchor: float, anchor_n: float,
                    rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    bio_c = sample_bio(baseline_ctrl, n_pairs, bio_alpha, rng)
    bio_d = sample_bio(baseline_dis, n_pairs, bio_alpha, rng)
    obs_c = add_tech_noise(bio_c, n_fragments, sd_anchor, anchor_n, rng)
    obs_d = add_tech_noise(bio_d, n_fragments, sd_anchor, anchor_n, rng)
    return obs_c, obs_d


def clr_transform(w: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    w = np.clip(w, eps, None)
    w = w / w.sum(axis=1, keepdims=True)
    log_w = np.log(w)
    return log_w - log_w.mean(axis=1, keepdims=True)


def two_sample_test(obs_ctrl: np.ndarray, obs_dis: np.ndarray,
                    eps: float = 1e-6) -> np.ndarray:
    """Pooled-variance t-test per cell type on CLR values; returns p-values.

    For a 2-group / no-covariate model `clr(w_c) = b0 + b1 * I(group=Disease)`,
    the OLS t-statistic for b1 is exactly the pooled-variance two-sample
    t-test. So with `equal_var=True` this matches too_diff_analysis.py
    (--test pairwise --reference-group Control, no covariates) bit-for-bit.
    """
    all_obs = np.vstack([obs_ctrl, obs_dis])
    clr = clr_transform(all_obs, eps=eps)
    n_c = obs_ctrl.shape[0]
    ctrl_clr = clr[:n_c]
    dis_clr = clr[n_c:]
    _, p_arr = _scipy_stats.ttest_ind(dis_clr, ctrl_clr, axis=0, equal_var=True)
    return np.asarray(p_arr, dtype=float)


# --- Power loop -----------------------------------------------------------

def estimate_power(
    n_pairs: int, n_fragments: float, delta: float,
    baseline_ctrl: np.ndarray, target_idx: int,
    bio_alpha: float, sd_anchor: float, anchor_n: float,
    n_iter: int, fdr_alpha: float,
    rng: np.random.Generator,
) -> dict:
    """Estimate power and family-wise FPR for one (n_pairs, n_fragments, delta) cell."""
    baseline_dis = apply_effect(baseline_ctrl, target_idx, delta)
    K = len(baseline_ctrl)
    other_mask = np.ones(K, dtype=bool)
    other_mask[target_idx] = False

    detections = 0
    any_other = 0
    n_other_flagged = 0
    qs_target = np.empty(n_iter)
    for it in range(n_iter):
        obs_c, obs_d = simulate_cohort(
            n_pairs, n_fragments, baseline_ctrl, baseline_dis,
            bio_alpha, sd_anchor, anchor_n, rng,
        )
        p = two_sample_test(obs_c, obs_d)
        q = bh_correct(p)
        qs_target[it] = q[target_idx]
        if np.isfinite(q[target_idx]) and q[target_idx] <= fdr_alpha:
            detections += 1
        flagged_other = int(np.nansum(q[other_mask] <= fdr_alpha))
        n_other_flagged += flagged_other
        if flagged_other > 0:
            any_other += 1

    return {
        "effect_size": delta,
        "n_pairs": n_pairs,
        "n_fragments": int(n_fragments),
        "power": detections / n_iter,
        "any_other_flagged_rate": any_other / n_iter,
        "expected_other_per_iter": n_other_flagged / n_iter,
        "median_q_target": float(np.nanmedian(qs_target)),
    }


def run_grid(
    n_pairs_list: list[int],
    n_fragments_list: list[float],
    deltas: list[float],
    baseline_ctrl: np.ndarray,
    target_idx: int,
    bio_alpha: float,
    sd_anchor: float,
    anchor_n: float,
    n_iter: int,
    fdr_alpha: float,
    seed: int,
    verbose: bool,
) -> pd.DataFrame:
    rng_master = np.random.default_rng(seed)
    rows = []
    total = len(n_pairs_list) * len(n_fragments_list) * len(deltas)
    done = 0
    for delta in deltas:
        for npairs in n_pairs_list:
            for nf in n_fragments_list:
                # Independent substream per cell so tweaking the grid doesn't
                # change individual-cell results when the same seed is reused.
                child = int(rng_master.integers(0, 2**31 - 1))
                rng = np.random.default_rng(child)
                row = estimate_power(
                    npairs, nf, delta, baseline_ctrl, target_idx,
                    bio_alpha, sd_anchor, anchor_n,
                    n_iter, fdr_alpha, rng,
                )
                rows.append(row)
                done += 1
                if verbose:
                    print(
                        f"  [{done:>3}/{total}] delta={delta:6.3f} "
                        f"n_pairs={npairs:>3} n_frags={nf:8.0e} "
                        f"power={row['power']:.3f} "
                        f"any_other_flagged={row['any_other_flagged_rate']:.3f}",
                        file=sys.stderr,
                    )
    return pd.DataFrame(rows)


# --- MDE summary ----------------------------------------------------------

def minimum_detectable_effect(df: pd.DataFrame, target_power: float) -> pd.DataFrame:
    """For each (n_pairs, n_fragments), report the smallest tested effect_size with
    power >= target_power. NaN if no tested effect reaches the target.
    """
    out = []
    for (n_pairs, n_fragments), grp in df.groupby(["n_pairs", "n_fragments"], sort=True):
        grp = grp.sort_values("effect_size")
        feasible = grp[grp["power"] >= target_power]
        if len(feasible):
            mde = float(feasible["effect_size"].min())
            mde_power = float(feasible.iloc[0]["power"])
        else:
            mde = float("nan")
            mde_power = float("nan")
        out.append({
            "n_pairs": int(n_pairs),
            "n_fragments": int(n_fragments),
            "min_detectable_effect": mde,
            "power_at_mde": mde_power,
        })
    return pd.DataFrame(out)


# --- Plotting -------------------------------------------------------------

def plot_heatmaps(df: pd.DataFrame, out_path: str,
                  n_pairs_list: list[int], n_fragments_list: list[float],
                  deltas: list[float], target_name: str) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed; skipping --plot.", file=sys.stderr)
        return

    n_panels = len(deltas)
    ncols = min(3, n_panels)
    nrows = math.ceil(n_panels / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.6 * ncols, 3.8 * nrows),
                             squeeze=False)
    axes_flat = axes.flatten()

    n_pairs_sorted = sorted(set(n_pairs_list))
    n_frags_sorted = sorted(set(n_fragments_list))

    for ax, delta in zip(axes_flat, deltas):
        sub = df[df["effect_size"] == delta]
        mat = (
            sub.pivot(index="n_pairs", columns="n_fragments", values="power")
            .reindex(index=n_pairs_sorted, columns=n_frags_sorted)
        )
        im = ax.imshow(mat.values, vmin=0, vmax=1, aspect="auto", origin="lower",
                       cmap="viridis")
        ax.set_xticks(range(len(n_frags_sorted)))
        ax.set_xticklabels([f"{x:.0e}" for x in n_frags_sorted],
                           rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(len(n_pairs_sorted)))
        ax.set_yticklabels(n_pairs_sorted, fontsize=8)
        ax.set_xlabel("n_fragments / sample")
        ax.set_ylabel("n_pairs")
        ax.set_title(f"effect = {delta:g}")
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat.values[i, j]
                if np.isfinite(v):
                    ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                            color="white" if v < 0.5 else "black", fontsize=7)
        plt.colorbar(im, ax=ax, fraction=0.04, pad=0.04, label="power")

    for ax in axes_flat[n_panels:]:
        ax.set_visible(False)

    fig.suptitle(f"too_diff_analysis power - target: {target_name}", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}", file=sys.stderr)


# --- CLI ------------------------------------------------------------------

def parse_csv_floats(s: str) -> list[float]:
    return [float(x) for x in s.split(",") if x.strip()]


def parse_csv_ints(s: str) -> list[int]:
    return [int(x) for x in s.split(",") if x.strip()]


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--n-pairs", default=DEFAULT_N_PAIRS_STR,
                   help=f"Comma-separated sample sizes (pairs of cases vs controls). "
                        f"Default: {DEFAULT_N_PAIRS_STR}")
    p.add_argument("--n-fragments", default=DEFAULT_N_FRAGMENTS_STR,
                   help=f"Comma-separated total filtered fragments per sample. "
                        f"Default: {DEFAULT_N_FRAGMENTS_STR}")
    p.add_argument("--effect-sizes", default=DEFAULT_EFFECT_SIZES_STR,
                   help=f"Comma-separated proportion-scale shifts in the target "
                        f"cell type. Default: {DEFAULT_EFFECT_SIZES_STR}")
    p.add_argument("--target-celltype", default=DEFAULT_TARGET_CELLTYPE,
                   help=f"Cell type to apply the effect to. Default: "
                        f"{DEFAULT_TARGET_CELLTYPE}")
    p.add_argument("--baseline-tsv", default=None,
                   help="Optional TSV with columns (cell_type, baseline_proportion). "
                        "If omitted, a healthy-plasma default is used.")
    p.add_argument("--bio-alpha", type=float, default=DEFAULT_BIO_ALPHA,
                   help=f"Dirichlet concentration for biological variability "
                        f"(higher = tighter). Default: {DEFAULT_BIO_ALPHA}")
    p.add_argument("--tech-sd-anchor", type=float, default=DEFAULT_TECH_SD_ANCHOR,
                   help=f"Technical SD of proportion at the anchor coverage. "
                        f"Default: {DEFAULT_TECH_SD_ANCHOR}")
    p.add_argument("--tech-anchor-n", type=float, default=DEFAULT_TECH_ANCHOR_N,
                   help=f"Anchor n_fragments for tech-sd-anchor. "
                        f"Default: {DEFAULT_TECH_ANCHOR_N:.0e}")
    p.add_argument("--n-iter", type=int, default=DEFAULT_N_ITER,
                   help=f"Replicates per (n_pairs, n_fragments, effect_size). "
                        f"Default: {DEFAULT_N_ITER}")
    p.add_argument("--fdr-alpha", type=float, default=DEFAULT_FDR_ALPHA,
                   help=f"BH q-value threshold. Default: {DEFAULT_FDR_ALPHA}")
    p.add_argument("--target-power", type=float, default=0.8,
                   help="Power threshold for the MDE summary. Default: 0.8")
    p.add_argument("--output", required=True,
                   help="Output TSV path for the full power table.")
    p.add_argument("--mde-output", default=None,
                   help="Optional output TSV path for the (n_pairs x n_fragments) "
                        "minimum-detectable-effect summary.")
    p.add_argument("--plot", default=None,
                   help="Optional PNG path for power heatmaps faceted by effect_size.")
    p.add_argument("--seed", type=int, default=DEFAULT_SEED,
                   help=f"Master RNG seed. Default: {DEFAULT_SEED}")
    p.add_argument("--verbose", action="store_true")

    args = p.parse_args()

    # ---- Baseline ----
    if args.baseline_tsv:
        baseline_dict = load_baseline_tsv(args.baseline_tsv)
    else:
        baseline_dict = DEFAULT_BASELINE
    baseline_ctrl, cell_types = normalized_baseline(baseline_dict)
    if args.target_celltype not in cell_types:
        sys.exit(
            f"ERROR: --target-celltype {args.target_celltype!r} not in baseline. "
            f"Available: {cell_types}"
        )
    target_idx = cell_types.index(args.target_celltype)

    n_pairs_list = parse_csv_ints(args.n_pairs)
    n_fragments_list = parse_csv_floats(args.n_fragments)
    deltas = parse_csv_floats(args.effect_sizes)

    if args.verbose:
        print(f"Baseline ({len(cell_types)} cell types):", file=sys.stderr)
        for c, p_val in zip(cell_types, baseline_ctrl):
            marker = "  <-- target" if c == args.target_celltype else ""
            print(f"  {c:<20s}{p_val:7.4f}{marker}", file=sys.stderr)
        print(
            f"\nTarget cell type:  {args.target_celltype} "
            f"(baseline={baseline_ctrl[target_idx]:.4f})\n"
            f"n_pairs:           {n_pairs_list}\n"
            f"n_fragments:       {[f'{x:.0e}' for x in n_fragments_list]}\n"
            f"effect_sizes:      {deltas}\n"
            f"n_iter / cell:     {args.n_iter}\n"
            f"FDR alpha:         {args.fdr_alpha}\n"
            f"Bio Dirichlet a:   {args.bio_alpha}\n"
            f"Tech SD anchor:    {args.tech_sd_anchor} at "
            f"{args.tech_anchor_n:.0e} fragments\n"
            f"Seed:              {args.seed}\n",
            file=sys.stderr,
        )

    # ---- Run grid ----
    df = run_grid(
        n_pairs_list, n_fragments_list, deltas,
        baseline_ctrl, target_idx,
        args.bio_alpha, args.tech_sd_anchor, args.tech_anchor_n,
        args.n_iter, args.fdr_alpha, args.seed, args.verbose,
    )
    # Stamp parameters for provenance
    df["target_celltype"] = args.target_celltype
    df["target_baseline"] = float(baseline_ctrl[target_idx])
    df["bio_alpha"] = args.bio_alpha
    df["tech_sd_anchor"] = args.tech_sd_anchor
    df["tech_anchor_n"] = args.tech_anchor_n
    df["n_iter"] = args.n_iter
    df["fdr_alpha"] = args.fdr_alpha

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False, float_format="%.4f")
    print(f"Wrote {args.output} ({len(df)} rows)", file=sys.stderr)

    # ---- MDE table ----
    mde_df = minimum_detectable_effect(df, args.target_power)
    if args.mde_output:
        Path(args.mde_output).parent.mkdir(parents=True, exist_ok=True)
        mde_df.to_csv(args.mde_output, sep="\t", index=False, float_format="%.4f")
        print(f"Wrote {args.mde_output} ({len(mde_df)} rows)", file=sys.stderr)

    if args.verbose:
        print(
            f"\nMinimum detectable effect at power >= {args.target_power} "
            f"(target: {args.target_celltype}, baseline="
            f"{baseline_ctrl[target_idx]:.4f}):",
            file=sys.stderr,
        )
        wide = mde_df.pivot(index="n_pairs", columns="n_fragments",
                            values="min_detectable_effect")

        def _fmt(v: float) -> str:
            return f"{v:.3f}" if np.isfinite(v) else "  -  "

        print(wide.to_string(float_format=_fmt), file=sys.stderr)

    # ---- Plot ----
    if args.plot:
        plot_heatmaps(df, args.plot, n_pairs_list, n_fragments_list, deltas,
                      args.target_celltype)

    return 0


if __name__ == "__main__":
    sys.exit(main())
