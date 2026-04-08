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

    Output columns (architecture §10.1):
        cell_type proportion ci_lower ci_upper p_goodness p_detection
        reliability n_markers mean_dispersion
    """
    rows = []
    cell_type_labels = list(result.cell_types) + ["Unknown"]
    K = len(result.cell_types)
    for i, ct in enumerate(cell_type_labels):
        if i < K:
            p_good = float(result.p_goodness[i]) if result.p_goodness is not None else float("nan")
            n_mark = int(result.n_markers[i]) if result.n_markers is not None else 0
            mean_disp = (
                float(result.mean_dispersion[i])
                if result.mean_dispersion is not None
                else float("nan")
            )
        else:
            p_good = float("nan")
            n_mark = 0
            mean_disp = float("nan")
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
                "mean_dispersion": mean_disp,
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

    columns = [
        "sample_id",
        "group",
        "coverage_tier",
        "mean_coverage",
        "n_markers_used",
        "pct_imputed",
    ]
    for ct in cell_type_labels:
        columns.append(ct)
        columns.append(f"{ct}_ci_lo")
        columns.append(f"{ct}_ci_hi")
        columns.append(f"{ct}_p_goodness")
        columns.append(f"{ct}_p_detection")
        columns.append(f"{ct}_reliability")
    columns.extend(
        [
            "calibration_flag",
            "hemolysis",
            "qc_flags",
            "overall_qc",
        ]
    )

    rows = []
    for r in results:
        row = {
            "sample_id": r.sample_id,
            "group": sample_groups.get(r.sample_id),
            "coverage_tier": r.coverage_tier.value if hasattr(r.coverage_tier, "value") else str(r.coverage_tier),
            "mean_coverage": float(r.mean_coverage),
            "n_markers_used": int(r.n_markers_used),
            "pct_imputed": _format_proportion(r.pct_imputed),
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
        row["calibration_flag"] = r.calibration_flag or "NA"
        row["hemolysis"] = (
            "TRUE" if r.hemolysis_flag is True
            else ("FALSE" if r.hemolysis_flag is False else "NA")
        )
        row["qc_flags"] = ",".join(r.qc_flags) if r.qc_flags else "NONE"
        row["overall_qc"] = r.overall_qc
        rows.append(row)

    pd.DataFrame(rows, columns=columns).to_csv(
        path, sep="\t", index=False, float_format="%.4f"
    )


def write_qc_summary(
    results: list[DeconvolutionResult],
    sample_groups: dict[str, str | None],
    path: str | Path,
) -> None:
    """Write a per-sample QC summary (architecture §10.6).

    Columns:
        sample_id group coverage_tier mean_coverage n_markers_used
        pct_imputed wbc_fraction unknown_fraction calibration_flag
        hemolysis qc_flags overall_qc
    """
    rows = []
    for r in results:
        unknown = float(r.proportions[-1])
        # WBC fraction = sum of proportions for cell types whose name starts
        # with "Blood" (matches the qc.py heuristic)
        blood_idx = [
            i for i, ct in enumerate(r.cell_types) if ct.lower().startswith("blood")
        ]
        wbc = float(np.sum(r.proportions[blood_idx])) if blood_idx else float("nan")
        rows.append(
            {
                "sample_id": r.sample_id,
                "group": sample_groups.get(r.sample_id),
                "coverage_tier": r.coverage_tier.value if hasattr(r.coverage_tier, "value") else str(r.coverage_tier),
                "mean_coverage": float(r.mean_coverage),
                "n_markers_used": int(r.n_markers_used),
                "pct_imputed": _format_proportion(r.pct_imputed),
                "wbc_fraction": _format_proportion(wbc) if not np.isnan(wbc) else float("nan"),
                "unknown_fraction": _format_proportion(unknown),
                "calibration_flag": r.calibration_flag or "NA",
                "hemolysis": (
                    "TRUE" if r.hemolysis_flag is True
                    else ("FALSE" if r.hemolysis_flag is False else "NA")
                ),
                "qc_flags": ",".join(r.qc_flags) if r.qc_flags else "NONE",
                "overall_qc": r.overall_qc,
            }
        )
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, float_format="%.4f")


def write_group_comparison(
    test_results: list,  # list[TestResult]
    path: str | Path,
) -> None:
    """Write group_comparison.tsv (architecture §10.3)."""
    rows = []
    for r in test_results:
        rows.append(
            {
                "test_type": r.test_type,
                "contrast": r.contrast,
                "cell_type": r.cell_type,
                "mean_A": r.mean_a,
                "mean_B": r.mean_b,
                "effect_size": r.effect_size,
                "se": r.se,
                "statistic": r.statistic,
                "p_value": r.p_value,
                "adjusted_p": r.adjusted_p_value,
                "significant": r.significant,
            }
        )
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, float_format="%.4g")


def write_residual_analysis(
    results: list[DeconvolutionResult],
    sample_groups: dict[str, str | None],
    path: str | Path,
    nmf_summary: dict | None = None,
) -> None:
    """Write residual_analysis.tsv (architecture §9.4 / §10.5).

    Columns:
        sample_id group coverage_tier unexplained_fraction mean_residual
        residual_sd nmf_top_component_loading qc_flag
    """
    rows = []
    loading_by_sample: dict[str, float] = {}
    if nmf_summary is not None:
        loadings = nmf_summary.get("loadings")
        sample_order = nmf_summary.get("sample_order")
        if loadings is not None and sample_order is not None and len(sample_order) > 0:
            # Use the per-sample max loading as a "top component" signature
            loadings_arr = np.asarray(loadings, dtype=np.float64)
            if loadings_arr.ndim == 2 and loadings_arr.shape[1] > 0:
                for sid, row in zip(sample_order, loadings_arr):
                    loading_by_sample[sid] = float(np.max(row))

    for r in results:
        if r.residuals is None or r.residuals.size == 0:
            unexplained = float(r.proportions[-1])
            mean_res = float("nan")
            sd_res = float("nan")
        else:
            finite = r.residuals[np.isfinite(r.residuals)]
            unexplained = float(r.proportions[-1])
            mean_res = float(np.mean(finite)) if finite.size else float("nan")
            sd_res = float(np.std(finite)) if finite.size else float("nan")
        rows.append(
            {
                "sample_id": r.sample_id,
                "group": sample_groups.get(r.sample_id),
                "coverage_tier": r.coverage_tier.value
                if hasattr(r.coverage_tier, "value")
                else str(r.coverage_tier),
                "unexplained_fraction": _format_proportion(unexplained),
                "mean_residual": mean_res,
                "residual_sd": sd_res,
                "nmf_top_component_loading": loading_by_sample.get(
                    r.sample_id, float("nan")
                ),
                "qc_flag": r.overall_qc,
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
    "write_group_comparison",
    "write_per_sample_too",
    "write_qc_summary",
    "write_residual_analysis",
]
