#!/usr/bin/env python3
"""
NNLS recovery simulation: how often do low-fraction cell types get clipped
to exactly zero by the non-negativity constraint?

Empirical question motivating this script: for the cfDNA deconvolution
reference panel, when we know a cell type is truly present at a small
fraction (0.1%, 0.5%, 1%, ...), what does NNLS recover?

  * If recovery is "0% with high probability below some threshold, then
    quantitatively above": NNLS-zeros are LEFT-CENSORED measurements --
    a true detection-floor effect. Tobit-style left-censored regression
    is statistically appropriate downstream.

  * If recovery is "small noisy positive estimates that occasionally
    happen to land at zero by chance": NNLS-zeros are noise artifacts
    of a constrained optimizer hitting its boundary. Hurdle's two-part
    decomposition or simple OLS-on-refined are more appropriate.

The simulation:
  1. Load the FinaleMe reference atlas (markers x cell-types of methylation
     beta values).
  2. For each TARGET cell type c and each TRUE fraction f in a log-spaced
     grid:
       - Construct a typical cfDNA-like background (dominated by blood
         lineages) and INJECT cell type c at fraction f.
       - Compute the expected mixture beta values y = X @ w_true.
       - Add per-marker observation noise representing finite sequencing
         coverage (binomial-style).
       - Solve constrained NNLS: w_hat = argmin ||X @ w - y||^2 s.t. w >= 0.
       - Record w_hat[c].
  3. Aggregate across replicates -> per-condition zero-rate, recovery
     distribution, calibration vs true fraction.

Outputs:
  --output-tsv  long-format table (one row per replicate)
  --output-pdf  multi-panel diagnostic plots
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

try:
    from scipy.optimize import nnls
except ImportError as exc:  # pragma: no cover
    sys.exit(f"requires scipy: pip install scipy ({exc})")

try:
    import matplotlib  # noqa: F401
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
except ImportError as exc:  # pragma: no cover
    sys.exit(f"requires matplotlib: pip install matplotlib ({exc})")


# ---- Defaults ------------------------------------------------------------

DEFAULT_TARGETS = [
    # Diverse selection: rare-but-detectable, common, and edge cases.
    "Pancreas-Alpha",
    "Breast-Luminal-Ep",
    "Liver-Hep",
    "Blood-T",
    "Adipocytes",
    "Bone-Osteob",
    "Heart-Cardio",
    "Lung-Ep-Alveo",
]

DEFAULT_FRACTIONS = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]


def load_atlas(path: str) -> tuple[np.ndarray, list[str], np.ndarray]:
    """Load the FinaleMe reference atlas TSV.

    Returns (X, cell_types, target_per_marker) where:
      X.shape = (n_markers, n_cell_types) of float beta values in [0, 1]
      cell_types is the column header order
      target_per_marker is the `target` column (primary cell type) per row,
        used for stratified marker subsampling if requested
    """
    df = pd.read_csv(path, sep="\t")
    metadata_cols = {
        "chr", "start", "end", "startCpG", "endCpG",
        "target", "name", "direction",
    }
    cell_cols = [c for c in df.columns if c not in metadata_cols]

    def _parse_cell(v):
        """Convert 'methylated/total' counts string to a beta value.
        Falls back to float() for already-numeric atlases."""
        if isinstance(v, (int, float)):
            return float(v)
        s = str(v)
        if "/" in s:
            num, den = s.split("/", 1)
            try:
                d = float(den)
                return float(num) / d if d > 0 else 0.0
            except ValueError:
                return float("nan")
        try:
            return float(s)
        except ValueError:
            return float("nan")

    sample0 = df[cell_cols[0]].iloc[0]
    if isinstance(sample0, str) and "/" in sample0:
        # Counts atlas; vectorize the parse for speed.
        X_df = df[cell_cols].applymap(_parse_cell)
        X = X_df.to_numpy(dtype=float)
    else:
        X = df[cell_cols].to_numpy(dtype=float)

    target_per_marker = df["target"].to_numpy() if "target" in df.columns else None
    # Sanity checks
    if not np.all((X >= 0) & (X <= 1)):
        bad = ((X < 0) | (X > 1)).sum()
        print(
            f"warning: {bad} atlas entries outside [0, 1]; clipping",
            file=sys.stderr,
        )
        X = np.clip(X, 0.0, 1.0)
    return X, cell_cols, target_per_marker


def build_background_composition(
    cell_types: list[str], rng: np.random.Generator
) -> np.ndarray:
    """Build a cfDNA-like background composition (sums to 1).

    Dominated by blood lineages (granulocytes, monocytes, T, NK), with
    smaller contributions from common organ sources (liver, endothel).
    Other cell types get small Dirichlet noise so the background isn't
    perfectly sparse.
    """
    base_targets = {
        "Blood-Granul": 0.55,
        "Blood-Mono+Macro": 0.10,
        "Blood-T": 0.18,
        "Blood-NK": 0.04,
        "Blood-B": 0.03,
        "Endothel": 0.03,
        "Liver-Hep": 0.04,
    }
    w = np.zeros(len(cell_types))
    leftover = 0.03  # tiny mass distributed evenly across "other" CTs
    for ct, p in base_targets.items():
        if ct in cell_types:
            w[cell_types.index(ct)] = p
    # Distribute leftover across remaining cell types.
    rest_idx = [i for i, ct in enumerate(cell_types) if ct not in base_targets]
    if rest_idx:
        # Use a flat split so the background isn't dominated by one
        # arbitrary cell type.
        per_rest = leftover / len(rest_idx)
        for i in rest_idx:
            w[i] = per_rest
    # Renormalize defensively (cell types may be missing from the panel).
    w = w / w.sum()
    return w


def inject_target(
    background: np.ndarray, target_idx: int, true_fraction: float
) -> np.ndarray:
    """Replace the target's slot with `true_fraction` and rescale the
    rest so the composition still sums to 1."""
    w = background.copy()
    rest_total = 1.0 - true_fraction
    # Zero out the target slot then scale the rest.
    target_old = w[target_idx]
    w[target_idx] = 0.0
    rest_sum = w.sum()
    if rest_sum > 0:
        w = w * (rest_total / rest_sum)
    w[target_idx] = true_fraction
    return w


def add_marker_noise(
    y_true: np.ndarray, coverage: int, rng: np.random.Generator
) -> np.ndarray:
    """Simulate per-marker observation noise for finite sequencing depth.

    Binomial sampling: at marker i with true methylation p_i, observe
    m_i = Binomial(coverage, p_i) methylated reads. Returned beta value
    is m_i / coverage. Adds the heteroscedasticity that real cfDNA WGS
    shows (more noise near 0.5, less near 0/1).
    """
    p = np.clip(y_true, 0.0, 1.0)
    m = rng.binomial(coverage, p)
    return m / float(coverage)


def run_nnls_simulation(
    X: np.ndarray,
    cell_types: list[str],
    targets: list[str],
    fractions: list[float],
    n_replicates: int,
    coverage: int,
    seed: int,
    n_marker_subsample: Optional[int] = None,
    target_per_marker: Optional[np.ndarray] = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """Run the per-(target, fraction) NNLS recovery experiment.

    Returns a long-format DataFrame with columns:
        target_celltype, true_fraction, replicate, recovered_fraction,
        is_zero, n_markers, coverage
    """
    rng = np.random.default_rng(seed)
    n_markers_full = X.shape[0]

    rows: list = []
    total_conditions = len(targets) * len(fractions)
    cond_idx = 0

    for target in targets:
        if target not in cell_types:
            print(
                f"  skipping {target!r}: not in atlas cell types",
                file=sys.stderr,
            )
            continue
        target_idx = cell_types.index(target)
        background = build_background_composition(cell_types, rng)

        for true_f in fractions:
            cond_idx += 1
            if verbose:
                print(
                    f"  [{cond_idx}/{total_conditions}] "
                    f"{target} @ true_fraction={true_f:g} "
                    f"({n_replicates} replicates)",
                    file=sys.stderr,
                )

            for rep in range(n_replicates):
                # Optionally subsample markers per replicate for speed.
                if (n_marker_subsample is not None
                        and n_marker_subsample < n_markers_full):
                    if (target_per_marker is not None
                            and "stratified" not in ""):
                        # Plain random subsample. We could do stratified
                        # by `target` but uniform is OK for this purpose.
                        marker_idx = rng.choice(
                            n_markers_full, n_marker_subsample, replace=False,
                        )
                    else:
                        marker_idx = rng.choice(
                            n_markers_full, n_marker_subsample, replace=False,
                        )
                    X_sub = X[marker_idx]
                else:
                    X_sub = X

                w_true = inject_target(background, target_idx, true_f)
                y_true = X_sub @ w_true
                y_obs = add_marker_noise(y_true, coverage, rng)
                # NNLS: returns (w, residual_norm)
                w_hat, _ = nnls(X_sub, y_obs)
                recovered = float(w_hat[target_idx])
                rows.append({
                    "target_celltype": target,
                    "true_fraction": float(true_f),
                    "replicate": rep,
                    "recovered_fraction": recovered,
                    "is_zero": bool(recovered == 0.0),
                    "n_markers": int(X_sub.shape[0]),
                    "coverage": coverage,
                })

    return pd.DataFrame(rows)


def plot_results(results: pd.DataFrame, output_pdf: str) -> None:
    """Write a multi-panel diagnostic PDF.

    Each cell type gets two pages (or stacked sub-panels):
      * Recovery distribution (boxplot/violin) at each true fraction
      * Zero-rate vs true fraction
    Plus a summary calibration plot across all cell types.
    """
    targets = sorted(results["target_celltype"].unique())
    fractions = sorted(results["true_fraction"].unique())
    Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(output_pdf) as pdf:
        # ---- Page 1: zero-rate sweep across all cell types --------
        fig, ax = plt.subplots(figsize=(9, 6))
        for target in targets:
            sub = results[results["target_celltype"] == target]
            zero_rate = sub.groupby("true_fraction")["is_zero"].mean()
            ax.plot(
                zero_rate.index, zero_rate.values, "o-",
                label=target, linewidth=1.5, markersize=6,
            )
        ax.set_xscale("log")
        ax.set_xlabel("True fraction injected", fontsize=11)
        ax.set_ylabel("P(NNLS recovered = 0)", fontsize=11)
        ax.set_title(
            "NNLS zero-rate vs true fraction\n"
            "(higher zero-rate at low fraction -> censoring-like behavior; "
            "Tobit appropriate)",
            fontsize=11,
        )
        ax.axhline(0.5, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
        ax.legend(fontsize=8, loc="best", ncol=2)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-0.05, 1.05)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # ---- Page 2: median recovery vs true fraction (calibration) -
        fig, ax = plt.subplots(figsize=(9, 6))
        for target in targets:
            sub = results[results["target_celltype"] == target]
            med = sub.groupby("true_fraction")["recovered_fraction"].median()
            ax.plot(
                med.index, med.values, "o-",
                label=target, linewidth=1.5, markersize=6,
            )
        # Diagonal reference: perfect calibration.
        lim = [min(fractions) * 0.5, max(fractions) * 1.2]
        ax.plot(lim, lim, "k--", linewidth=0.8, alpha=0.6, label="y=x (perfect)")
        ax.set_xscale("log")
        ax.set_yscale("symlog", linthresh=1e-4)
        ax.set_xlabel("True fraction injected", fontsize=11)
        ax.set_ylabel("Median recovered fraction", fontsize=11)
        ax.set_title(
            "Calibration: median recovered vs true fraction",
            fontsize=11,
        )
        ax.legend(fontsize=8, loc="best", ncol=2)
        ax.grid(True, alpha=0.3, which="both")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # ---- Per-cell-type detail pages ---------------------------
        for target in targets:
            sub = results[results["target_celltype"] == target]
            fig, axes = plt.subplots(1, 2, figsize=(13, 5))

            # Left: recovery distribution per true fraction
            data_per_f = []
            labels = []
            for f in fractions:
                vals = sub[sub["true_fraction"] == f][
                    "recovered_fraction"
                ].to_numpy()
                if len(vals) > 0:
                    data_per_f.append(vals)
                    labels.append(f"{f:g}")

            if data_per_f:
                # Violin showing distribution
                parts = axes[0].violinplot(
                    data_per_f, positions=range(len(data_per_f)),
                    showmedians=True, showextrema=False, widths=0.7,
                )
                for body in parts["bodies"]:
                    body.set_alpha(0.4)
                    body.set_edgecolor("black")
                    body.set_linewidth(0.5)
                if "cmedians" in parts:
                    parts["cmedians"].set_color("black")
                # Diagonal reference at the "true" mark for each f
                for j, f in enumerate(fractions):
                    if j < len(data_per_f):
                        axes[0].plot(
                            [j - 0.35, j + 0.35], [f, f],
                            color="red", linewidth=1.5, alpha=0.8,
                        )
                axes[0].set_xticks(range(len(labels)))
                axes[0].set_xticklabels(labels, rotation=45)
                axes[0].set_xlabel("True fraction")
                axes[0].set_ylabel("NNLS recovered fraction")
                axes[0].set_title(
                    f"{target}: recovery distribution\n"
                    f"(red bar = true value)",
                    fontsize=10,
                )
                axes[0].grid(True, alpha=0.3, axis="y")

            # Right: zero-rate per fraction with confidence
            zero_rate = sub.groupby("true_fraction")["is_zero"].agg(
                ["mean", "count"]
            )
            zero_rate["se"] = np.sqrt(
                zero_rate["mean"] * (1 - zero_rate["mean"])
                / zero_rate["count"]
            )
            axes[1].errorbar(
                zero_rate.index, zero_rate["mean"], yerr=1.96 * zero_rate["se"],
                fmt="o-", color="C0", linewidth=1.5, markersize=7, capsize=3,
            )
            axes[1].axhline(0.5, color="gray", linestyle="--", linewidth=0.8)
            axes[1].set_xscale("log")
            axes[1].set_xlabel("True fraction")
            axes[1].set_ylabel("P(recovered = 0)")
            axes[1].set_title(
                f"{target}: zero-rate (95% CI)",
                fontsize=10,
            )
            axes[1].set_ylim(-0.05, 1.05)
            axes[1].grid(True, alpha=0.3)

            fig.suptitle(target, fontsize=13, fontweight="bold")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--atlas",
        default="data/Atlas.CGI_shore.U250.l3.hg19.tsv",
        help="Path to the FinaleMe reference panel TSV. "
             "Default: data/Atlas.CGI_shore.U250.l3.hg19.tsv",
    )
    parser.add_argument(
        "--targets",
        default=",".join(DEFAULT_TARGETS),
        help=f"Comma-separated cell types to inject. Default: "
             f"{','.join(DEFAULT_TARGETS)}",
    )
    parser.add_argument(
        "--fractions",
        default=",".join(str(f) for f in DEFAULT_FRACTIONS),
        help=f"Comma-separated true fractions to test. Default: "
             f"{','.join(str(f) for f in DEFAULT_FRACTIONS)}",
    )
    parser.add_argument(
        "--replicates",
        type=int,
        default=300,
        help="Replicates per (target, fraction) condition. Default: 300.",
    )
    parser.add_argument(
        "--coverage",
        type=int,
        default=10,
        help="Effective per-marker sequencing coverage (binomial draws). "
             "Default: 10.",
    )
    parser.add_argument(
        "--n-marker-subsample",
        type=int,
        default=2000,
        help="If set, subsample this many markers per replicate to speed "
             "up NNLS. The full panel has ~15,000 markers; 2000 is fast "
             "and preserves NNLS qualitative behavior. Default: 2000.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed. Default: 42.",
    )
    parser.add_argument(
        "--output-tsv",
        default="results/nnls_recovery_simulation.tsv",
        help="Output path for the long-format results TSV.",
    )
    parser.add_argument(
        "--output-pdf",
        default="results/nnls_recovery_simulation.pdf",
        help="Output path for the diagnostic plots PDF.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Per-condition progress logging on stderr.",
    )
    args = parser.parse_args()

    print(f"loading atlas: {args.atlas}", file=sys.stderr)
    X, cell_types, target_per_marker = load_atlas(args.atlas)
    print(
        f"  {X.shape[0]:,} markers x {X.shape[1]} cell types",
        file=sys.stderr,
    )

    targets = [t.strip() for t in args.targets.split(",") if t.strip()]
    fractions = [float(f) for f in args.fractions.split(",") if f.strip()]
    print(
        f"running simulation: {len(targets)} target cell types x "
        f"{len(fractions)} fractions x {args.replicates} replicates "
        f"= {len(targets) * len(fractions) * args.replicates:,} NNLS calls",
        file=sys.stderr,
    )

    results = run_nnls_simulation(
        X=X,
        cell_types=cell_types,
        targets=targets,
        fractions=fractions,
        n_replicates=args.replicates,
        coverage=args.coverage,
        seed=args.seed,
        n_marker_subsample=args.n_marker_subsample,
        target_per_marker=target_per_marker,
        verbose=args.verbose,
    )

    Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.output_tsv, sep="\t", index=False)
    print(
        f"wrote {len(results):,} rows to {args.output_tsv}",
        file=sys.stderr,
    )

    plot_results(results, args.output_pdf)
    print(f"wrote {args.output_pdf}", file=sys.stderr)

    # ---- Summary table to stderr ---------------------------------
    print("\n=== Zero-rate by (target, true_fraction) ===", file=sys.stderr)
    summary = (
        results.groupby(["target_celltype", "true_fraction"])
        .agg(
            zero_rate=("is_zero", "mean"),
            median_recovered=("recovered_fraction", "median"),
            mean_recovered=("recovered_fraction", "mean"),
            sd_recovered=("recovered_fraction", "std"),
            n=("replicate", "count"),
        )
        .reset_index()
    )
    pd.set_option("display.max_rows", None)
    pd.set_option("display.float_format", lambda x: f"{x:.4f}")
    print(summary.to_string(index=False), file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
