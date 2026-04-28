#!/usr/bin/env python3
"""plot_roc_curve.py -- render the FinaleMe -aucMode ROC/PR curve TSV to PDF.

Reads the TSV produced by `FinaleMe -aucMode -aucCurveTsv ...` (one row per
threshold, with FPR / TPR / Precision / Recall columns and AUROC/AUPRC in
commented header lines) and produces a 2-panel figure:

    +-------------------------+-------------------------+
    | ROC curve               | Precision-Recall curve  |
    | AUROC = ...             | AUPRC = ...             |
    +-------------------------+-------------------------+

Defaults to PDF; pass -o ... .png for raster output.

Example
-------
    python3 scripts/plot_roc_curve.py results/auc_curve.tsv -o results/roc.pdf
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def parse_header_metadata(path: str) -> dict:
    """Read leading '# key\\tvalue' header lines from the TSV."""
    meta = {}
    with open(path, "r") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            parts = line.lstrip("#").strip().split("\t", 1)
            if len(parts) == 2:
                meta[parts[0].strip()] = parts[1].strip()
    return meta


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("tsv", help="TSV produced by FinaleMe -aucMode -aucCurveTsv")
    p.add_argument("-o", "--output", required=True,
                   help="Output PDF (or PNG) path.")
    p.add_argument("--title", default=None,
                   help="Optional figure title (e.g. sample name).")
    p.add_argument("--dpi", type=int, default=150,
                   help="Output DPI (relevant for PNG/raster). Default: 150")
    args = p.parse_args()

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit("matplotlib is required: pip install matplotlib")

    meta = parse_header_metadata(args.tsv)
    df = pd.read_csv(args.tsv, sep="\t", comment="#")

    expected = {"threshold", "FPR", "TPR", "Precision", "Recall"}
    missing = expected - set(df.columns)
    if missing:
        sys.exit(f"ERROR: TSV is missing columns {sorted(missing)}; got {sorted(df.columns)}")

    auroc_str = meta.get("AUROC", "?")
    auprc_str = meta.get("AUPRC", "?")
    n_meth = meta.get("n_methylated", "?")
    n_unmeth = meta.get("n_unmethylated", "?")

    try:
        auroc_f = float(auroc_str)
        auprc_f = float(auprc_str)
    except ValueError:
        auroc_f = auprc_f = float("nan")

    # Sort points along each curve so the trapezoidal area visualization
    # matches the line we draw.
    roc = df.sort_values("FPR").reset_index(drop=True)
    pr = df.sort_values("Recall").reset_index(drop=True)
    # PR curve typically prepends (recall=0, precision=1) to anchor the
    # leftmost point.  Same convention as the Java AUPRC computation.
    if pr["Recall"].min() > 0:
        pr = pd.concat([
            pd.DataFrame({"Recall": [0.0], "Precision": [1.0]}),
            pr[["Recall", "Precision"]],
        ], ignore_index=True)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

    ax = axes[0]
    ax.plot(roc["FPR"], roc["TPR"], "-o", color="#1f77b4", markersize=3,
            linewidth=1.5, label=f"AUROC = {auroc_f:.4f}")
    ax.plot([0, 1], [0, 1], "--", color="grey", linewidth=0.8, label="chance")
    ax.fill_between(roc["FPR"], roc["TPR"], alpha=0.10, color="#1f77b4")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("ROC curve")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.plot(pr["Recall"], pr["Precision"], "-o", color="#d62728", markersize=3,
            linewidth=1.5, label=f"AUPRC = {auprc_f:.4f}")
    ax.fill_between(pr["Recall"], pr["Precision"], alpha=0.10, color="#d62728")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Precision-Recall curve")
    ax.legend(loc="lower left", fontsize=9)
    ax.grid(True, alpha=0.3)

    title = args.title or "FinaleMe AUC evaluation"
    fig.suptitle(
        f"{title}\nN methylated = {n_meth}, N unmethylated = {n_unmeth}",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {args.output} (AUROC={auroc_str}, AUPRC={auprc_str})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
