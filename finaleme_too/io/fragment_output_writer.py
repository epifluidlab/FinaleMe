"""Output writers for fragment-level TOO deconvolution results (design §3.5.9).

Three output files:
  - {prefix}_proportions.tsv     : per cell type proportions with CI, p/q values
  - {prefix}_fragment_responsibilities.tsv.gz : per-fragment gamma matrix
  - {prefix}_diagnostics.json    : EM convergence and filtering statistics
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd


def write_fragment_proportions(
    cell_types: list[str],
    w: np.ndarray,
    ci_lower: np.ndarray | None,
    ci_upper: np.ndarray | None,
    p_values: np.ndarray | None,
    q_values: np.ndarray | None,
    fdr_threshold: float,
    path: str | Path,
) -> None:
    """Write {prefix}_proportions.tsv per design §3.5.8.

    Columns: cell_type, proportion, proportion_renorm, CI_lower, CI_upper,
             p_value, q_value, significant.

    ``proportion_renorm`` is rescaled to sum to 1.0 across known cell types
    only (unknown excluded). Unknown row gets NA for proportion_renorm.
    """
    K = len(cell_types)
    K_total = K + 1  # includes unknown
    names = list(cell_types) + ["Unknown"]

    unknown_frac = float(w[K]) if K < len(w) else 0.0
    renorm_denom = max(1.0 - unknown_frac, 1e-12)

    rows = []
    for i, ct in enumerate(names):
        prop = float(w[i])
        if i < K:
            prop_renorm = prop / renorm_denom
        else:
            prop_renorm = float("nan")

        ci_lo = float(ci_lower[i]) if ci_lower is not None else float("nan")
        ci_hi = float(ci_upper[i]) if ci_upper is not None else float("nan")
        pval = float(p_values[i]) if p_values is not None else float("nan")
        qval = float(q_values[i]) if q_values is not None else float("nan")

        sig = "yes" if (not np.isnan(qval) and qval <= fdr_threshold) else "no"

        rows.append({
            "cell_type": ct,
            "proportion": f"{prop:.4f}",
            "proportion_renorm": f"{prop_renorm:.4f}" if not np.isnan(prop_renorm) else "",
            "CI_lower": f"{ci_lo:.4f}" if not np.isnan(ci_lo) else "",
            "CI_upper": f"{ci_hi:.4f}" if not np.isnan(ci_hi) else "",
            "p_value": f"{pval:.4g}" if not np.isnan(pval) else "",
            "q_value": f"{qval:.4g}" if not np.isnan(qval) else "",
            "significant": sig,
        })

    df = pd.DataFrame(rows)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def write_fragment_responsibilities(
    gamma: np.ndarray,
    cell_types: list[str],
    path: str | Path,
) -> None:
    """Write {prefix}_fragment_responsibilities.tsv.gz.

    Gzipped TSV with one row per fragment and columns for each cell type
    plus Unknown. Values are the per-fragment responsibilities gamma[f, j].
    """
    names = list(cell_types) + ["Unknown"]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.write("\t".join(names) + "\n")
        for row in gamma:
            fh.write("\t".join(f"{v:.6f}" for v in row) + "\n")


def write_fragment_diagnostics(
    diagnostics: dict,
    path: str | Path,
) -> None:
    """Write {prefix}_diagnostics.json.

    Expected keys:
        n_fragments_total, n_fragments_filtered, n_fragments_used,
        em_iterations, em_log_likelihood, em_converged, pi_unknown,
        informativeness_distribution (dict of percentile -> value),
        filter_breakdown (dict of reason -> count).
    """
    Path(path).parent.mkdir(parents=True, exist_ok=True)

    # Ensure numpy types are JSON-serializable
    def _serialize(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    clean = {k: _serialize(v) if not isinstance(v, dict)
             else {kk: _serialize(vv) for kk, vv in v.items()}
             for k, v in diagnostics.items()}

    with open(path, "w") as fh:
        json.dump(clean, fh, indent=2)


__all__ = [
    "write_fragment_proportions",
    "write_fragment_responsibilities",
    "write_fragment_diagnostics",
]
