"""Output writers for per-sample TSVs, cohort summaries, and reports."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from finaleme_too.core.deconvolution import DeconvolutionResult


def _format_proportion(x: float) -> float:
    return float(round(x, 4))


def write_per_sample_too(result: DeconvolutionResult, path: str | Path) -> None:
    """Write a per-sample TSV with one row per cell type (incl. unknown).

    Output columns:
        cell_type proportion ci_lower ci_upper p_goodness p_detection
        reliability n_markers
    """
    rows = []
    cell_type_labels = list(result.cell_types) + ["Unknown"]
    K = len(result.cell_types)
    for i, ct in enumerate(cell_type_labels):
        if i < K:
            p_good = float(result.p_goodness[i]) if result.p_goodness is not None else float("nan")
            n_mark = int(result.n_markers[i]) if result.n_markers is not None else 0
        else:
            p_good = float("nan")
            n_mark = 0
        rows.append(
            {
                "cell_type": ct,
                "proportion": _format_proportion(result.proportions[i]),
                "ci_lower": _format_proportion(result.ci_lower[i]),
                "ci_upper": _format_proportion(result.ci_upper[i]),
                "p_goodness": p_good,
                "p_detection": float(result.p_detection[i]),
                "reliability": str(result.reliability[i]),
                "n_markers": n_mark,
            }
        )
    df = pd.DataFrame(rows)
    df.to_csv(path, sep="\t", index=False, float_format="%.4f")


def write_cohort_proportions(
    results: list[DeconvolutionResult],
    sample_groups: dict[str, str | None],
    path: str | Path,
) -> None:
    """Write cohort_proportions.tsv with one row per sample.

    Each cell type contributes proportion + ci_lo + ci_hi + p_goodness + p_detection
    + reliability columns; the unknown component is appended at the end.
    """
    if not results:
        pd.DataFrame().to_csv(path, sep="\t", index=False)
        return

    base = results[0]
    cell_type_labels = list(base.cell_types) + ["Unknown"]
    K = len(base.cell_types)

    columns = ["sample_id", "group", "coverage_tier"]
    for ct in cell_type_labels:
        columns.append(ct)
        columns.append(f"{ct}_ci_lo")
        columns.append(f"{ct}_ci_hi")
        columns.append(f"{ct}_p_goodness")
        columns.append(f"{ct}_p_detection")
        columns.append(f"{ct}_reliability")
    columns.append("qc_flags")

    rows = []
    for r in results:
        row = {
            "sample_id": r.sample_id,
            "group": sample_groups.get(r.sample_id),
            "coverage_tier": r.coverage_tier.value if hasattr(r.coverage_tier, "value") else str(r.coverage_tier),
        }
        for i, ct in enumerate(cell_type_labels):
            row[ct] = _format_proportion(r.proportions[i])
            row[f"{ct}_ci_lo"] = _format_proportion(r.ci_lower[i])
            row[f"{ct}_ci_hi"] = _format_proportion(r.ci_upper[i])
            if i < K and r.p_goodness is not None:
                row[f"{ct}_p_goodness"] = float(r.p_goodness[i])
            else:
                row[f"{ct}_p_goodness"] = float("nan")
            row[f"{ct}_p_detection"] = float(r.p_detection[i])
            row[f"{ct}_reliability"] = str(r.reliability[i])
        row["qc_flags"] = ",".join(r.qc_flags) if r.qc_flags else "NONE"
        rows.append(row)

    pd.DataFrame(rows, columns=columns).to_csv(
        path, sep="\t", index=False, float_format="%.4f"
    )


def write_qc_summary(
    results: list[DeconvolutionResult],
    sample_groups: dict[str, str | None],
    path: str | Path,
) -> None:
    """Write a per-sample QC summary."""
    rows = []
    for r in results:
        unknown = float(r.proportions[-1])
        rows.append(
            {
                "sample_id": r.sample_id,
                "group": sample_groups.get(r.sample_id),
                "coverage_tier": r.coverage_tier.value if hasattr(r.coverage_tier, "value") else str(r.coverage_tier),
                "unknown_fraction": _format_proportion(unknown),
                "n_qc_flags": len(r.qc_flags),
                "qc_flags": ",".join(r.qc_flags) if r.qc_flags else "NONE",
            }
        )
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, float_format="%.4f")


def write_calibration_report(report: dict, path: str | Path) -> None:
    Path(path).write_text(json.dumps(report, indent=2, default=_json_default))


def _json_default(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    raise TypeError(f"Type {type(obj)} not JSON serializable")


__all__ = [
    "write_calibration_report",
    "write_cohort_proportions",
    "write_per_sample_too",
    "write_qc_summary",
]
