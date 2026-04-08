"""Regression tests for design gaps 5, 6, 7.

Gap 5 — per-tier min-read filtering now applied in _deconvolve_sample.
Gap 6a — FragmentLevelDeconvolver wired for ULTRALOW tier via .pat.gz.
Gap 6b — discover_residual_components wired + residual_analysis.tsv writer.
Gap 7  — enriched DeconvolutionResult + cohort/QC schema columns.
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from finaleme_too.config import CoverageTier, MeasurementMode
from finaleme_too.core.deconvolution import DeconvolutionResult, MLEDeconvolver
from finaleme_too.core.fragment_likelihood import Fragment
from finaleme_too.core.observation_model import BetaBinomialModel
from finaleme_too.io.marker_regions import MarkerRegions
from finaleme_too.io.methylation_loader import MarkerObservations
from finaleme_too.io.output_writer import (
    write_cohort_proportions,
    write_per_sample_too,
    write_qc_summary,
    write_residual_analysis,
)
from finaleme_too.io.pat_loader import load_fragments_from_pat
from finaleme_too.preprocessing.coverage import per_marker_min_reads


# ---------------------------------------------------------------------------
# Gap 5 — per-tier min-read filtering is actually applied
# ---------------------------------------------------------------------------


def test_gap5_per_marker_min_reads_table():
    assert per_marker_min_reads(CoverageTier.HIGH) == 3
    assert per_marker_min_reads(CoverageTier.LOW) == 2
    assert per_marker_min_reads(CoverageTier.ULTRALOW) == 1


def test_gap5_pipeline_filters_low_coverage_markers(
    tmp_path: Path, synthetic_marker_regions, synthetic_reference
):
    """Samples with some markers well below the mean should still keep
    those markers via effective-coverage down-tiering (architecture §4.2),
    while markers with absolutely zero coverage get filtered."""
    from finaleme_too.config import MeasurementMode, TOOConfig
    from finaleme_too.io.sample_sheet import Sample, SampleSheet
    from finaleme_too.pipeline import TOOPipeline

    # Half markers high-coverage, half zero-coverage.
    # Under per-marker min_reads_vector, zero-coverage markers fail even
    # in ULTRALOW (min=1), so they are filtered. Other markers pass.
    n_markers = synthetic_marker_regions.n_markers
    rng = np.random.default_rng(0)
    p = synthetic_reference.methylation[:, 0].astype(np.float64)
    n_arr = np.where(np.arange(n_markers) < n_markers // 2, 50, 0).astype(np.int32)
    k_arr = rng.binomial(n_arr.astype(np.int64), p).astype(np.int32)

    # Write a FinaleMe prediction.bed.gz with these counts (one CpG per marker)
    pred_path = tmp_path / "s1.prediction.bed.gz"
    with gzip.open(pred_path, "wt") as fh:
        fh.write(
            "#chr\tstart\tend\tmethy_perc_predict\tmethy_count_predict\ttotal_count_predict"
            "\to1\to2\to3\n"
        )
        for i in range(n_markers):
            mi_chr = str(synthetic_marker_regions.chrom[i])
            mi_start = int(synthetic_marker_regions.start[i]) + 100
            fh.write(
                f"{mi_chr}\t{mi_start}\t{mi_start + 1}\t"
                f"{100 * int(k_arr[i]) / max(int(n_arr[i]), 1):.1f}\t{int(k_arr[i])}\t{int(n_arr[i])}\t0\t0\t0\n"
            )

    # Build sample sheet + reference + marker files on disk.
    # Use FINALEME mode so the loader picks the finaleme_bed parser that
    # reads the predicted-count columns (4, 5), not the bissnp parser.
    sheet_path = tmp_path / "sheet.tsv"
    sheet_path.write_text(
        f"sample_id\tmethylation_file\tmode\tgroup\ns1\t{pred_path}\tFINALEME\tA\n"
    )
    ref_path = tmp_path / "ref.tsv"
    with open(ref_path, "w") as fh:
        fh.write("chrom\tstart\tend\t" + "\t".join(synthetic_reference.cell_types) + "\n")
        for i in range(n_markers):
            row = synthetic_reference.methylation[i]
            fh.write(
                f"{synthetic_marker_regions.chrom[i]}\t"
                f"{int(synthetic_marker_regions.start[i])}\t"
                f"{int(synthetic_marker_regions.end[i])}\t"
                + "\t".join(f"{x:.4f}" for x in row)
                + "\n"
            )
    mr_path = tmp_path / "markers.bed"
    with open(mr_path, "w") as fh:
        for i in range(n_markers):
            fh.write(
                f"{synthetic_marker_regions.chrom[i]}\t"
                f"{int(synthetic_marker_regions.start[i])}\t"
                f"{int(synthetic_marker_regions.end[i])}\tm{i}\n"
            )

    # Run the pipeline end-to-end. Force HIGH tier via thresholds so that
    # the per-tier filter actually runs on the synthetic sample (the
    # default heuristic classifies tiny synthetic samples as ULTRALOW).
    config = TOOConfig()
    config.uncertainty.n_bootstrap = 20
    config.coverage.tier_high = 0.0
    config.coverage.tier_low = -1.0
    pipeline = TOOPipeline(config)
    from finaleme_too.io.reference_panel import ReferencePanelLoader
    from finaleme_too.io.marker_regions import MarkerRegionsLoader

    reference = ReferencePanelLoader.load_matrix(ref_path)
    marker_regions = MarkerRegionsLoader.load(mr_path)
    sample_sheet = SampleSheet.from_tsv(sheet_path)
    cohort = pipeline.run(sample_sheet, reference, marker_regions, tmp_path / "out")
    assert len(cohort.samples) == 1
    r = cohort.samples[0]
    # Half the markers have n=50 (always pass), half have n=0 (always fail
    # even under ULTRALOW down-tiering). So exactly half pass the filter.
    assert r.n_markers_used == n_markers // 2
    # Proportions are still valid
    assert abs(np.sum(r.proportions) - 1.0) < 1e-6


# ---------------------------------------------------------------------------
# Gap 6a — pat_loader and FragmentLevelDeconvolver integration
# ---------------------------------------------------------------------------


def test_gap6a_pat_loader_parses_fragments(tmp_path: Path):
    """Build a trivial .pat.gz file + CpG index and verify fragments load."""
    # 10 CpGs on chr1 at positions 100, 110, 120, ...
    cpg_path = tmp_path / "cpg.bed"
    cpg_path.write_text(
        "\n".join(f"chr1\t{100 + i * 10}\t{101 + i * 10}" for i in range(10)) + "\n"
    )
    from finaleme_too.io.reference_panel import load_cpg_index

    cpg_index = load_cpg_index(cpg_path)
    assert cpg_index["total_sites"] == 10

    # Two marker regions: [100, 140) covers CpGs 0-3, [150, 190) covers 5-8
    mr = MarkerRegions(
        chrom=np.array(["chr1", "chr1"], dtype=object),
        start=np.array([100, 150], dtype=np.int64),
        end=np.array([140, 190], dtype=np.int64),
        marker_name=None,
    )

    pat_path = tmp_path / "sample.pat.gz"
    with gzip.open(pat_path, "wt") as fh:
        # Pattern "CCTT" starting at CpG index 1 covers CpGs at positions
        # 100, 110, 120, 130 — all inside marker 0 ([100, 140))
        fh.write("chr1\t1\tCCTT\t3\n")
        # Pattern "CCCT" starting at CpG index 6 covers CpGs at positions
        # 150, 160, 170, 180 — all inside marker 1 ([150, 190))
        fh.write("chr1\t6\tCCCT\t2\n")

    fragments = load_fragments_from_pat(pat_path, mr, cpg_index)
    # 3 + 2 = 5 fragments expanded from the count column
    assert len(fragments) == 5
    for f in fragments[:3]:
        assert np.all(f.cpg_indices == 0)
        np.testing.assert_array_equal(f.methylated, [1, 1, 0, 0])  # C,C,T,T
    for f in fragments[3:]:
        assert np.all(f.cpg_indices == 1)
        np.testing.assert_array_equal(f.methylated, [1, 1, 1, 0])  # C,C,C,T


def test_gap6a_resolved_pat_file_fallback_naming(tmp_path: Path):
    """Sample.resolved_pat_file() should find a sibling .prediction.pat.gz
    when no explicit pat_file column is set."""
    from finaleme_too.config import MeasurementMode
    from finaleme_too.io.sample_sheet import Sample

    bed_path = tmp_path / "s1.prediction.bed.gz"
    bed_path.write_text("")
    pat_path = tmp_path / "s1.prediction.pat.gz"
    pat_path.write_text("")

    s = Sample(
        sample_id="s1",
        methylation_file=bed_path,
        mode=MeasurementMode.FINALEME,
    )
    assert s.resolved_pat_file() == pat_path


# ---------------------------------------------------------------------------
# Gap 6b — residual_analysis.tsv writer
# ---------------------------------------------------------------------------


def _mk_decon_result(sample_id: str, unknown: float, residuals: np.ndarray) -> DeconvolutionResult:
    K = 3
    prop = np.array([0.3, 0.3, 0.4 - unknown, unknown])
    return DeconvolutionResult(
        sample_id=sample_id,
        cell_types=["CT1", "CT2", "CT3"],
        proportions=prop,
        ci_lower=prop - 0.05,
        ci_upper=prop + 0.05,
        p_goodness=np.full(K, 0.5),
        p_detection=np.full(K + 1, 0.9),
        reliability=np.array(["HIGH"] * (K + 1), dtype=object),
        n_markers=np.full(K, 20, dtype=np.int32),
        coverage_tier=CoverageTier.HIGH,
        qc_flags=[],
        mean_dispersion=np.full(K, 50.0),
        mean_coverage=25.0,
        n_markers_used=20,
        pct_imputed=0.0,
        calibration_flag="PASS",
        hemolysis_flag=False,
        overall_qc="PASS",
        residuals=residuals,
    )


def test_gap6b_residual_analysis_writer(tmp_path: Path):
    """write_residual_analysis should emit one row per sample with the
    expected columns and use any provided NMF summary."""
    rng = np.random.default_rng(0)
    residuals_a = rng.normal(0, 0.05, size=30)
    residuals_b = rng.normal(0.1, 0.05, size=30)
    results = [
        _mk_decon_result("s1", unknown=0.05, residuals=residuals_a),
        _mk_decon_result("s2", unknown=0.35, residuals=residuals_b),
    ]
    sample_groups = {"s1": "A", "s2": "A"}
    nmf_summary = {
        "loadings": np.array([[0.2, 0.1], [0.8, 0.05]]),
        "sample_order": ["s1", "s2"],
        "components": np.zeros((2, 30)),
        "explained_variance_ratio": np.array([0.7, 0.2]),
    }
    out_path = tmp_path / "residual_analysis.tsv"
    write_residual_analysis(results, sample_groups, out_path, nmf_summary=nmf_summary)
    df = pd.read_csv(out_path, sep="\t")
    assert list(df["sample_id"]) == ["s1", "s2"]
    assert "unexplained_fraction" in df.columns
    assert "mean_residual" in df.columns
    assert "residual_sd" in df.columns
    assert "nmf_top_component_loading" in df.columns
    assert "qc_flag" in df.columns
    # s2 should have the higher NMF loading
    s2_loading = float(df.loc[df["sample_id"] == "s2", "nmf_top_component_loading"].iloc[0])
    s1_loading = float(df.loc[df["sample_id"] == "s1", "nmf_top_component_loading"].iloc[0])
    assert s2_loading > s1_loading


def test_gap6b_residual_analysis_without_nmf(tmp_path: Path):
    """Residual TSV should still be emitted with NaN loadings when NMF is skipped."""
    results = [
        _mk_decon_result("s1", unknown=0.05, residuals=np.zeros(10)),
    ]
    out_path = tmp_path / "residual_analysis.tsv"
    write_residual_analysis(results, {"s1": "A"}, out_path, nmf_summary=None)
    df = pd.read_csv(out_path, sep="\t")
    assert len(df) == 1
    assert np.isnan(float(df["nmf_top_component_loading"].iloc[0]))


# ---------------------------------------------------------------------------
# Gap 7 — enriched DeconvolutionResult + output schema
# ---------------------------------------------------------------------------


def test_gap7_per_sample_tsv_has_mean_dispersion(tmp_path: Path):
    r = _mk_decon_result("s1", unknown=0.05, residuals=np.zeros(5))
    out = tmp_path / "s1.too.tsv"
    write_per_sample_too(r, out)
    df = pd.read_csv(out, sep="\t")
    # mean_dispersion must be present AND populated with real values for
    # the real cell types (NaN is fine for the Unknown row)
    assert "mean_dispersion" in df.columns
    ct_rows = df[df["cell_type"] != "Unknown"]
    assert (ct_rows["mean_dispersion"] == 50.0).all()


def test_gap7_cohort_tsv_has_enriched_columns(tmp_path: Path):
    results = [
        _mk_decon_result("s1", unknown=0.05, residuals=np.zeros(5)),
        _mk_decon_result("s2", unknown=0.10, residuals=np.zeros(5)),
    ]
    out = tmp_path / "cohort_proportions.tsv"
    write_cohort_proportions(results, {"s1": "A", "s2": "B"}, out)
    df = pd.read_csv(out, sep="\t")
    for col in (
        "mean_coverage",
        "n_markers_used",
        "pct_imputed",
        "calibration_flag",
        "hemolysis",
        "overall_qc",
    ):
        assert col in df.columns, f"missing column: {col}"
    assert list(df["overall_qc"]) == ["PASS", "PASS"]


def test_gap7_qc_summary_has_architecture_columns(tmp_path: Path):
    results = [
        _mk_decon_result("s1", unknown=0.05, residuals=np.zeros(5)),
        _mk_decon_result("s2", unknown=0.10, residuals=np.zeros(5)),
    ]
    out = tmp_path / "qc.tsv"
    write_qc_summary(results, {"s1": "A", "s2": "B"}, out)
    df = pd.read_csv(out, sep="\t")
    for col in (
        "sample_id",
        "group",
        "coverage_tier",
        "mean_coverage",
        "n_markers_used",
        "pct_imputed",
        "wbc_fraction",  # NaN is fine when no Blood-* cell types
        "unknown_fraction",
        "calibration_flag",
        "hemolysis",
        "qc_flags",
        "overall_qc",
    ):
        assert col in df.columns, f"missing column: {col}"
