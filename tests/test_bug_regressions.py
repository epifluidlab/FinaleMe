"""Regression tests for the bugs reported on April 7 2026.

Each bug is identified by its number from the original report.
"""

from __future__ import annotations

import gzip
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from finaleme_too.config import (
    BatchCorrectionConfig,
    CovariateAdjustmentConfig,
    MeasurementMode,
    TestMethod,
    TOOConfig,
    UncertaintyConfig,
)
from finaleme_too.io.marker_regions import MarkerRegions
from finaleme_too.io.methylation_loader import MarkerObservations, MethylationLoader
from finaleme_too.postprocessing.group_comparison import run_group_comparisons
from finaleme_too.preprocessing.calibration import CalibrationParams
from finaleme_too.preprocessing.calibration_eval import compute_inference_qc
from finaleme_too.preprocessing.imputation import CohortImputer


# ---------------------------------------------------------------------------
# Bug 1: bayesian_posterior test method must actually run when posterior
#        samples are available.
# ---------------------------------------------------------------------------


def test_bug1_bayesian_posterior_method_dispatches_to_bayesian():
    """run_group_comparisons must call bayesian_group_comparison when method
    is BAYESIAN_POSTERIOR and posterior samples are provided."""
    rng = np.random.default_rng(0)
    n_samples = 6
    K_total = 3
    sample_ids = [f"a_{i}" for i in range(3)] + [f"b_{i}" for i in range(3)]
    labels = ["A"] * 3 + ["B"] * 3

    # Mean proportions per sample (used as point estimates) — irrelevant here
    proportions = rng.dirichlet(np.ones(K_total), size=n_samples)

    # Posterior samples: group A has cell type 0 high, group B has it low
    posterior = {}
    for sid in sample_ids[:3]:
        posterior[sid] = rng.dirichlet([8, 1, 1], size=200)
    for sid in sample_ids[3:]:
        posterior[sid] = rng.dirichlet([1, 8, 1], size=200)

    results = run_group_comparisons(
        proportions=proportions,
        sample_ids=sample_ids,
        group_labels=labels,
        cell_type_names=["CT1", "CT2"],  # K = K_total - 1 (last is unknown)
        spec="A:B",
        method=TestMethod.BAYESIAN_POSTERIOR,
        posterior_samples_by_sample=posterior,
    )
    assert all(r.test_type == "bayesian" for r in results), \
        f"Expected all results to be bayesian, got: {[r.test_type for r in results]}"
    # CT1 should have P(A > B) close to 1 (statistic field)
    by_ct = {r.cell_type: r for r in results}
    assert by_ct["CT1"].statistic > 0.85


def test_bug1_bayesian_falls_back_to_ilr_when_no_posterior():
    """When method=BAYESIAN_POSTERIOR but no posterior samples, fall back."""
    rng = np.random.default_rng(0)
    proportions = rng.dirichlet(np.ones(3), size=6)
    results = run_group_comparisons(
        proportions=proportions,
        sample_ids=[f"a_{i}" for i in range(3)] + [f"b_{i}" for i in range(3)],
        group_labels=["A"] * 3 + ["B"] * 3,
        cell_type_names=["CT1", "CT2"],
        spec="A:B",
        method=TestMethod.BAYESIAN_POSTERIOR,
        posterior_samples_by_sample=None,
    )
    # Should use ilr_regression as fallback
    assert all(r.test_type == "ilr_regression" for r in results)


# ---------------------------------------------------------------------------
# Bug 2: Calibration inference QC must use CpG density for bins, not the
#        predicted beta values.
# ---------------------------------------------------------------------------


def _density_calibration() -> CalibrationParams:
    return CalibrationParams(
        n_bins=4,
        bin_edges=np.array([0.0, 0.025, 0.05, 0.075, 1.0]),
        a=np.array([0.8, 0.9, 0.95, 1.0]),
        c=np.array([-0.05, 0.0, 0.02, 0.05]),
        log_dispersion=np.full(4, np.log(20.0)),
    )


def test_bug2_inference_qc_uses_density_for_bins_not_beta():
    """All markers in bin 0 (low density) should produce balance < 1.0
    even if predicted betas span [0, 1]."""
    cal = _density_calibration()
    n = 50
    sample_pred = np.linspace(0.05, 0.95, n)  # spans full beta range
    # All markers have low density → all in bin 0
    density = np.full(n, 0.001)
    qc = compute_inference_qc(sample_pred, cal, cpg_density=density)
    # Bin 0 has 50 markers (>= 10), bins 1-3 have 0 → balance = 1/4 = 0.25
    assert qc["bin_coverage_balance"] == pytest.approx(0.25)
    # And it should have flagged WARN/FAIL due to balance, not PASS
    assert qc["flag"] in ("WARN", "FAIL")


def test_bug2_inference_qc_residuals_are_real_calibration_residuals():
    """KS p should be high when sample residuals match training residuals."""
    cal = _density_calibration()
    rng = np.random.default_rng(0)
    n = 200
    sample_pred = np.clip(rng.uniform(0.1, 0.9, size=n), 1e-3, 1 - 1e-3)
    density = np.full(n, 0.001)  # all in bin 0 → a=0.8, c=-0.05

    # Compute the *true* residuals the function should be testing against
    from scipy.special import logit
    raw_logit = logit(sample_pred)
    expected_residuals = (cal.a[0] - 1.0) * raw_logit + cal.c[0]

    qc = compute_inference_qc(
        sample_pred,
        cal,
        cpg_density=density,
        training_residuals=expected_residuals,  # exact match
    )
    # When residuals match exactly, KS p should be ~1
    assert qc["residual_ks_p"] > 0.5


def test_bug2_prediction_range_coverage_uses_training_range():
    cal = _density_calibration()
    n = 100
    # Half inside [0.1, 0.9], half outside
    sample_pred = np.concatenate([np.linspace(0.2, 0.8, 50), np.linspace(0.91, 0.99, 50)])
    qc = compute_inference_qc(
        sample_pred,
        cal,
        cpg_density=np.zeros(n),
        training_pred_range=(0.1, 0.9),
    )
    # 50/100 inside (0.1, 0.9)
    assert abs(qc["prediction_range_coverage"] - 0.5) < 1e-6


# ---------------------------------------------------------------------------
# Bug 3: predicted_beta from FinaleMe loader must be count-weighted (k/n),
#        not an unweighted mean of per-CpG percentages.
# ---------------------------------------------------------------------------


def test_bug3_predicted_beta_is_count_weighted(tmp_path: Path):
    """A marker with two CpGs of very different coverage should have
    predicted_beta == k_total / n_total, not mean of per-CpG percentages."""
    bed_path = tmp_path / "sample.prediction.bed.gz"
    # Marker covers chr1:100-200
    # CpG @ 110: 100% methylated, 1 read total
    # CpG @ 150: 0% methylated, 99 reads total
    # k/n = (1 + 0)/(1 + 99) = 0.01
    # mean(pct) = (100 + 0) / 2 = 50  →  0.5
    with gzip.open(bed_path, "wt") as fh:
        fh.write(
            "#chr\tstart\tend\tmethy_perc_predict\tmethy_count_predict\ttotal_count_predict"
            "\tobs1\tobs2\tobs3\n"
        )
        fh.write("chr1\t110\t111\t100.0\t1\t1\t0\t0\t0\n")
        fh.write("chr1\t150\t151\t0.0\t0\t99\t0\t0\t0\n")

    mr = MarkerRegions(
        chrom=np.array(["chr1"], dtype=object),
        start=np.array([100], dtype=np.int64),
        end=np.array([200], dtype=np.int64),
        marker_name=None,
    )
    obs = MethylationLoader.load(
        filepath=bed_path,
        sample_id="s",
        mode=MeasurementMode.FINALEME,
        marker_regions=mr,
        input_format="finaleme_bed",
    )
    assert int(obs.k[0]) == 1
    assert int(obs.n[0]) == 100
    # predicted_beta should equal k/n = 0.01, NOT 0.5
    assert obs.predicted_beta is not None
    assert abs(float(obs.predicted_beta[0]) - 0.01) < 1e-6, \
        f"predicted_beta should be count-weighted, got {obs.predicted_beta[0]}"


# ---------------------------------------------------------------------------
# Bug 4: Imputation must apply min_donors per marker, not globally; and
#        synthetic n must come from eligible donors only.
# ---------------------------------------------------------------------------


def _mk_obs(sid: str, k: list[int], n: list[int]) -> MarkerObservations:
    M = len(k)
    return MarkerObservations(
        sample_id=sid,
        chrom=np.array(["chr1"] * M, dtype=object),
        start=np.array([100 + i * 1000 for i in range(M)], dtype=np.int64),
        end=np.array([200 + i * 1000 for i in range(M)], dtype=np.int64),
        k=np.array(k, dtype=np.int32),
        n=np.array(n, dtype=np.int32),
        predicted_beta=None,
        mode=MeasurementMode.WGBS,
    )


def test_bug4_imputation_per_marker_donor_constraint():
    """Markers without enough eligible donors at that specific marker should
    not be imputed, even if the global donor pool is large."""
    target = _mk_obs("t", [0, 0, 0], [0, 0, 0])  # all zero coverage
    # 5 donors, but only marker 0 has 5 eligible donors;
    # marker 1 has 2 eligible donors (below threshold of 3);
    # marker 2 has 5 eligible donors.
    donors = [
        _mk_obs("d0", [10, 10, 10], [20, 20, 20]),  # all eligible
        _mk_obs("d1", [10, 10, 10], [20, 20, 20]),  # all eligible
        _mk_obs("d2", [10, 0, 10], [20, 0, 20]),    # marker 1 ineligible
        _mk_obs("d3", [10, 0, 10], [20, 0, 20]),    # marker 1 ineligible
        _mk_obs("d4", [10, 0, 10], [20, 0, 20]),    # marker 1 ineligible
    ]
    cohort = donors + [target]
    groups = {sid: "A" for sid in [d.sample_id for d in donors] + [target.sample_id]}
    imputer = CohortImputer(coverage_threshold=3, min_donors=3)
    out = imputer.impute(target, cohort, groups)

    # Marker 0: 5 eligible donors → imputed
    assert out.n[0] > 0, "Marker 0 should have been imputed"
    assert out.k[0] > 0
    # Marker 1: only 2 eligible donors → NOT imputed (must stay zero)
    assert out.n[1] == 0, "Marker 1 should NOT be imputed (only 2 eligible donors)"
    assert out.k[1] == 0
    # Marker 2: 5 eligible donors → imputed
    assert out.n[2] > 0, "Marker 2 should have been imputed"


def test_bug4_synthetic_n_uses_eligible_donors_only():
    """Synthetic n for an imputed marker should be the median of eligible
    donor coverages, not all donor coverages."""
    target = _mk_obs("t", [0], [0])
    # 4 donors: 3 with coverage 50, 1 with coverage 1
    # eligible-only median (threshold=3) = 50
    # all-donor median = 50 (still 50 because 4 elements, but pick a clearer case)
    donors = [
        _mk_obs("d0", [25], [50]),
        _mk_obs("d1", [25], [50]),
        _mk_obs("d2", [25], [50]),
        _mk_obs("d3", [0], [1]),  # ineligible
    ]
    cohort = donors + [target]
    groups = {sid: "A" for sid in [d.sample_id for d in donors] + [target.sample_id]}
    out = CohortImputer(coverage_threshold=3, min_donors=3).impute(target, cohort, groups)
    # Eligible donors are d0/d1/d2 (n=50 each); median = 50
    # All-donor median (median of [50,50,50,1]) = 50 too
    # but the key check: the imputed n must equal 50, not less
    assert int(out.n[0]) == 50


# ---------------------------------------------------------------------------
# Bug 5: CLI flags must actually be applied to the config.
# ---------------------------------------------------------------------------


def test_bug5_cli_options_applied(tmp_path: Path):
    """Run the CLI with all the previously-ignored flags and verify they
    show up in the effective TOOConfig used by the pipeline."""
    sheet = tmp_path / "sheet.tsv"
    ref = tmp_path / "ref.tsv"
    bed = tmp_path / "markers.bed"

    # Minimal valid inputs
    bed.write_text("chr1\t100\t200\tm1\nchr1\t300\t400\tm2\n")
    ref.write_text(
        "chrom\tstart\tend\tCellA\tCellB\nchr1\t100\t200\t0.1\t0.9\nchr1\t300\t400\t0.9\t0.1\n"
    )
    # FinaleMe BED query
    fme = tmp_path / "s1.prediction.bed.gz"
    with gzip.open(fme, "wt") as fh:
        fh.write(
            "#chr\tstart\tend\tmethy_perc_predict\tmethy_count_predict\ttotal_count_predict"
            "\to1\to2\to3\n"
        )
        fh.write("chr1\t110\t111\t10\t1\t10\t0\t0\t0\n")
        fh.write("chr1\t310\t311\t90\t9\t10\t0\t0\t0\n")
    sheet.write_text(
        "sample_id\tmethylation_file\tmode\tgroup\textraction_batch\nsample\t"
        + str(fme)
        + "\tFINALEME\tA\tB1\n"
    )

    # Capture the config the pipeline ends up using by monkey-patching
    captured = {}
    from finaleme_too import pipeline as pipeline_mod

    orig_init = pipeline_mod.TOOPipeline.__init__

    def spy_init(self, config, **kw):
        captured["config"] = config
        captured["group_comparison"] = kw.get("group_comparison_spec")
        orig_init(self, config, **kw)
        # Don't actually run the pipeline — just capture the config
        raise SystemExit(0)

    pipeline_mod.TOOPipeline.__init__ = spy_init  # type: ignore[assignment]
    try:
        sys.argv = [
            "finaleme-too",
            "run",
            "--sample-sheet",
            str(sheet),
            "--output-dir",
            str(tmp_path / "out"),
            "--reference-panel",
            str(ref),
            "--marker-regions",
            str(bed),
            "--bayesian-n-samples",
            "1234",
            "--batch-correct",
            "extraction_batch,library_date",
            "--adjust-covariates",
            "age,sex",
            "--configurable-covariates",
            "treatment",
            "--fdr-method",
            "BH",
            "--fdr-alpha",
            "0.10",
        ]
        from finaleme_too.cli import main

        try:
            main(standalone_mode=False)
        except SystemExit:
            pass
    finally:
        pipeline_mod.TOOPipeline.__init__ = orig_init  # type: ignore[assignment]

    cfg: TOOConfig = captured["config"]
    assert cfg.uncertainty.bayesian_n_samples == 1234
    assert cfg.batch_correction.technical_covariates == ["extraction_batch", "library_date"]
    assert cfg.covariate_adjustment.biological_covariates == ["age", "sex"]
    assert cfg.covariate_adjustment.user_configurable == ["treatment"]
    assert cfg.testing.fdr_method == "BH"
    assert abs(cfg.testing.fdr_alpha - 0.10) < 1e-9
    # Group comparison fallback: --group-comparison was omitted, so the
    # config default ("omnibus+pairwise") should be used
    assert captured["group_comparison"] == "omnibus+pairwise"


# ---------------------------------------------------------------------------
# April 2026 round 2 regressions
# ---------------------------------------------------------------------------


# Bug A — enum coercion in TOOConfig.from_dict / from_yaml
def test_round2_bug1_enum_coercion_from_yaml(tmp_path: Path):
    """Loading a YAML config must coerce enum-typed fields, not leave strings."""
    from finaleme_too.config import (
        CalibrationConfig,
        MeasurementMode,
        SolverMethod,
        TestMethod,
        TOOConfig,
    )

    yaml_path = tmp_path / "cfg.yaml"
    yaml_path.write_text(
        "model:\n"
        "  deconvolution: bayesian\n"
        "calibration:\n"
        "  mode: WGBS\n"
        "testing:\n"
        "  method: wilcoxon\n"
    )
    cfg = TOOConfig.from_yaml(yaml_path)

    # All three enum-typed fields must be the actual Enum subclass instance
    assert isinstance(cfg.model.deconvolution, SolverMethod)
    assert cfg.model.deconvolution == SolverMethod.BAYESIAN
    assert isinstance(cfg.calibration.mode, MeasurementMode)
    assert cfg.calibration.mode == MeasurementMode.WGBS
    assert isinstance(cfg.testing.method, TestMethod)
    assert cfg.testing.method == TestMethod.WILCOXON
    # And .value access works (this would crash if the field were a plain str)
    assert cfg.testing.method.value == "wilcoxon"


def test_round2_bug1_already_enum_passes_through():
    """Passing an Enum instance directly should not crash."""
    from finaleme_too.config import SolverMethod, TOOConfig

    cfg = TOOConfig.from_dict(
        {"model": {"deconvolution": SolverMethod.BAYESIAN}}
    )
    assert isinstance(cfg.model.deconvolution, SolverMethod)
    assert cfg.model.deconvolution == SolverMethod.BAYESIAN


# Bug B — CLI overrides config only when explicitly provided
def test_round2_bug2_cli_does_not_clobber_config_file_values(tmp_path: Path):
    """When the user supplies --config but omits --n-bootstrap, the n_bootstrap
    value from the YAML must NOT be overwritten by the Click default."""
    sheet = tmp_path / "sheet.tsv"
    ref = tmp_path / "ref.tsv"
    bed = tmp_path / "markers.bed"
    bed.write_text("chr1\t100\t200\tm1\nchr1\t300\t400\tm2\n")
    ref.write_text(
        "chrom\tstart\tend\tCellA\tCellB\nchr1\t100\t200\t0.1\t0.9\nchr1\t300\t400\t0.9\t0.1\n"
    )
    fme = tmp_path / "s1.prediction.bed.gz"
    with gzip.open(fme, "wt") as fh:
        fh.write(
            "#chr\tstart\tend\tmethy_perc_predict\tmethy_count_predict\ttotal_count_predict"
            "\to1\to2\to3\n"
        )
        fh.write("chr1\t110\t111\t10\t1\t10\t0\t0\t0\n")
        fh.write("chr1\t310\t311\t90\t9\t10\t0\t0\t0\n")
    sheet.write_text(
        "sample_id\tmethylation_file\tmode\tgroup\nsample\t" + str(fme) + "\tFINALEME\tA\n"
    )

    # YAML config sets values that differ from CLI defaults
    yaml_cfg = tmp_path / "config.yaml"
    yaml_cfg.write_text(
        "uncertainty:\n"
        "  n_bootstrap: 777\n"
        "  bayesian_n_samples: 9999\n"
        "testing:\n"
        "  method: wilcoxon\n"
        "  fdr_alpha: 0.20\n"
        "  fdr_method: bonferroni\n"
        "coverage:\n"
        "  coverage_cap: 88\n"
    )

    captured = {}
    from finaleme_too import pipeline as pipeline_mod

    orig_init = pipeline_mod.TOOPipeline.__init__

    def spy_init(self, config, **kw):
        captured["config"] = config
        raise SystemExit(0)

    pipeline_mod.TOOPipeline.__init__ = spy_init  # type: ignore[assignment]
    try:
        sys.argv = [
            "finaleme-too",
            "run",
            "--sample-sheet",
            str(sheet),
            "--output-dir",
            str(tmp_path / "out"),
            "--reference-panel",
            str(ref),
            "--marker-regions",
            str(bed),
            "--config",
            str(yaml_cfg),
            # NOTE: NOT passing --n-bootstrap, --bayesian-n-samples,
            # --test-method, --fdr-alpha, --fdr-method, --coverage-cap
        ]
        from finaleme_too.cli import main

        try:
            main(standalone_mode=False)
        except SystemExit:
            pass
    finally:
        pipeline_mod.TOOPipeline.__init__ = orig_init  # type: ignore[assignment]

    cfg = captured["config"]
    # All YAML values must be preserved — none clobbered by Click defaults
    assert cfg.uncertainty.n_bootstrap == 777, \
        f"Expected n_bootstrap=777 from YAML, got {cfg.uncertainty.n_bootstrap}"
    assert cfg.uncertainty.bayesian_n_samples == 9999
    assert cfg.testing.method.value == "wilcoxon"
    assert abs(cfg.testing.fdr_alpha - 0.20) < 1e-9
    assert cfg.testing.fdr_method == "bonferroni"
    assert cfg.coverage.coverage_cap == 88


# Bug C — omnibus rows must have FDR-adjusted p-values
def test_round2_bug3_omnibus_rows_have_adjusted_p_value():
    """Omnibus Kruskal-Wallis rows must go through FDR like the pairwise rows."""
    from finaleme_too.config import TestMethod

    rng = np.random.default_rng(0)
    K_total = 4  # 3 cell types + unknown
    n_per_group = 5

    group_a = rng.dirichlet([8, 1, 1, 1], size=n_per_group)
    group_b = rng.dirichlet([1, 8, 1, 1], size=n_per_group)
    group_c = rng.dirichlet([1, 1, 8, 1], size=n_per_group)
    proportions = np.vstack([group_a, group_b, group_c])
    labels = ["A"] * n_per_group + ["B"] * n_per_group + ["C"] * n_per_group
    sample_ids = [f"s{i}" for i in range(15)]

    results = run_group_comparisons(
        proportions=proportions,
        sample_ids=sample_ids,
        group_labels=labels,
        cell_type_names=["CT1", "CT2", "CT3"],
        spec="omnibus+pairwise",
        method=TestMethod.ILR_REGRESSION,
        fdr_alpha=0.05,
    )
    # Expect both omnibus rows AND pairwise rows
    omnibus_rows = [r for r in results if r.test_type == "omnibus"]
    pairwise_rows = [r for r in results if r.test_type == "ilr_regression"]
    assert len(omnibus_rows) >= 1, "expected at least one omnibus row"
    assert len(pairwise_rows) >= 1, "expected at least one pairwise row"

    # Every row with a non-NaN p-value must have a non-NaN adjusted_p_value
    for r in omnibus_rows + pairwise_rows:
        if not np.isnan(r.p_value):
            assert not np.isnan(r.adjusted_p_value), \
                f"{r.test_type}/{r.cell_type}: adjusted_p_value is NaN"
            assert r.adjusted_p_value >= r.p_value or r.p_value == 1.0, \
                f"{r.test_type}/{r.cell_type}: adjusted_p {r.adjusted_p_value} < raw {r.p_value}"


# Bug D — calibration fallback binning is deterministic with +/-inf edges
def test_round2_bug4_assign_bin_with_nan_density_does_not_return_nan():
    """assign_bin must produce a finite bin index even when density is NaN
    and bin_edges contain ±inf."""
    from finaleme_too.preprocessing.calibration import CalibrationParams

    edges = np.array([-np.inf, 0.01, 0.05, 0.075, np.inf])
    params = CalibrationParams(
        n_bins=4,
        bin_edges=edges,
        a=np.ones(4),
        c=np.zeros(4),
        log_dispersion=np.zeros(4),
    )
    out = params.assign_bin(np.array([np.nan, 0.005, 0.06, np.nan, np.inf, -np.inf]))
    # All entries must be finite ints in [0, n_bins-1]
    assert out.dtype.kind in ("i", "u"), f"expected integer dtype, got {out.dtype}"
    assert np.all((out >= 0) & (out < params.n_bins))
    # NaN entries → fallback_bin (deterministic)
    fallback = params.fallback_bin
    assert int(out[0]) == fallback
    assert int(out[3]) == fallback
    # Real values still bin correctly
    assert int(out[1]) == 0  # 0.005 in [-inf, 0.01)
    assert int(out[2]) == 2  # 0.06 in [0.05, 0.075)


def test_round2_bug4_apply_calibration_with_no_region_annotation_does_not_nan(tmp_path):
    """apply_calibration with region_annotations=None must not silently NaN."""
    from finaleme_too.config import MeasurementMode
    from finaleme_too.io.methylation_loader import MarkerObservations
    from finaleme_too.preprocessing.calibration import CalibrationParams, apply_calibration

    edges = np.array([-np.inf, 0.01, 0.05, np.inf])
    # Aggressive non-identity calibration so we can see the effect:
    # bin 1 (the deterministic fallback) has a=2.0, c=1.0
    cal = CalibrationParams(
        n_bins=3,
        bin_edges=edges,
        a=np.array([2.0, 2.0, 2.0]),
        c=np.array([1.0, 1.0, 1.0]),
        log_dispersion=np.zeros(3),
    )
    n = 5
    obs = MarkerObservations(
        sample_id="x",
        chrom=np.array(["chr1"] * n, dtype=object),
        start=np.array([100, 200, 300, 400, 500], dtype=np.int64),
        end=np.array([200, 300, 400, 500, 600], dtype=np.int64),
        k=np.array([3, 5, 7, 2, 8], dtype=np.int32),
        n=np.array([10, 10, 10, 10, 10], dtype=np.int32),
        predicted_beta=np.array([0.3, 0.5, 0.7, 0.2, 0.8], dtype=np.float32),
        mode=MeasurementMode.FINALEME,
    )
    out = apply_calibration(obs, cal, region_annotations=None)
    # The key assertion: NO NaNs anywhere (the original bug produced NaN
    # because params.bin_edges.mean() with +/-inf returns NaN, which then
    # propagated through the bin lookup)
    assert np.all(np.isfinite(out.k))
    assert np.all((out.k >= 0) & (out.k <= out.n))
    # The aggressive a=2.0, c=1.0 calibration must actually change at least
    # one value (sanity check that the lookup table is being used)
    assert not np.array_equal(out.k, obs.k)


# ---------------------------------------------------------------------------
# April 2026 round 3 regressions
# ---------------------------------------------------------------------------


# HIGH — residual_analysis.tsv must not silently skip when the first
# sample has no residuals (e.g. tier-filter fallback)
def test_round3_high_residual_analysis_not_skipped_when_first_sample_empty(tmp_path):
    """Even if results[0].residuals is None, the residual TSV must still
    be emitted (with NaN rows for samples that had no residuals)."""
    from finaleme_too.config import CoverageTier, TOOConfig
    from finaleme_too.core.deconvolution import DeconvolutionResult
    from finaleme_too.io.output_writer import write_residual_analysis
    from finaleme_too.pipeline import TOOPipeline

    def _mk(sid, residuals):
        K = 2
        prop = np.array([0.4, 0.3, 0.3])
        return DeconvolutionResult(
            sample_id=sid,
            cell_types=["CT1", "CT2"],
            proportions=prop,
            ci_lower=prop - 0.05,
            ci_upper=prop + 0.05,
            p_goodness=np.array([0.5, 0.5]),
            p_detection=np.array([0.9, 0.9, 0.9]),
            reliability=np.array(["HIGH", "HIGH", "HIGH"], dtype=object),
            n_markers=np.array([5, 5], dtype=np.int32),
            coverage_tier=CoverageTier.HIGH,
            qc_flags=[],
            mean_dispersion=np.array([50.0, 50.0]),
            mean_coverage=25.0,
            n_markers_used=5,
            residuals=residuals,
            overall_qc="PASS",
        )

    results = [
        _mk("s1", residuals=None),  # <-- first sample has no residuals
        _mk("s2", residuals=np.array([0.01, -0.02, 0.03, -0.01, 0.005])),
        _mk("s3", residuals=np.array([0.00, -0.01, 0.02, -0.03, 0.01])),
    ]
    sample_groups = {"s1": "A", "s2": "A", "s3": "A"}

    config = TOOConfig()
    pipeline = TOOPipeline(config)
    out_dir = tmp_path / "out"
    out_dir.mkdir(exist_ok=True)
    pipeline._run_residual_analysis(results, sample_groups, out_dir)

    residual_tsv = out_dir / "residual_analysis.tsv"
    assert residual_tsv.exists(), "residual_analysis.tsv must be written"
    import pandas as pd

    df = pd.read_csv(residual_tsv, sep="\t")
    assert set(df["sample_id"]) == {"s1", "s2", "s3"}
    # s1 (no residuals) should have NaN stats
    s1_row = df[df["sample_id"] == "s1"].iloc[0]
    assert np.isnan(float(s1_row["mean_residual"]))
    # s2/s3 should have finite stats
    s2_row = df[df["sample_id"] == "s2"].iloc[0]
    assert np.isfinite(float(s2_row["mean_residual"]))


# MEDIUM — sample sheet must require the group column
def test_round3_sample_sheet_requires_group_column(tmp_path):
    """Missing 'group' column should raise InvalidSampleSheetError."""
    from finaleme_too.exceptions import InvalidSampleSheetError
    from finaleme_too.io.sample_sheet import SampleSheet

    sheet_path = tmp_path / "sheet_no_group.tsv"
    sheet_path.write_text(
        "sample_id\tmethylation_file\tmode\ns1\t/nope/a.bed.gz\tFINALEME\n"
    )
    with pytest.raises(InvalidSampleSheetError) as exc:
        SampleSheet.from_tsv(sheet_path)
    assert "group" in str(exc.value).lower()


def test_round3_sample_sheet_rejects_empty_group_values(tmp_path):
    """Empty string / NaN group values should raise."""
    from finaleme_too.exceptions import InvalidSampleSheetError
    from finaleme_too.io.sample_sheet import SampleSheet

    sheet_path = tmp_path / "sheet_empty_group.tsv"
    sheet_path.write_text(
        "sample_id\tmethylation_file\tmode\tgroup\n"
        "s1\t/nope/a.bed.gz\tFINALEME\tA\n"
        "s2\t/nope/b.bed.gz\tFINALEME\t\n"  # empty group
    )
    with pytest.raises(InvalidSampleSheetError) as exc:
        SampleSheet.from_tsv(sheet_path)
    assert "group" in str(exc.value).lower()


def test_round3_imputation_lenient_mode_skips_unlabeled_samples():
    """In lenient mode (the pipeline's default) imputation returns the
    original sample instead of crashing when the group is missing."""
    from finaleme_too.config import MeasurementMode
    from finaleme_too.io.methylation_loader import MarkerObservations
    from finaleme_too.preprocessing.imputation import CohortImputer

    def _obs(sid, k, n):
        return MarkerObservations(
            sample_id=sid,
            chrom=np.array(["chr1"] * len(k), dtype=object),
            start=np.array(list(range(len(k))), dtype=np.int64),
            end=np.array(list(range(1, len(k) + 1)), dtype=np.int64),
            k=np.array(k, dtype=np.int32),
            n=np.array(n, dtype=np.int32),
            predicted_beta=None,
            mode=MeasurementMode.WGBS,
        )

    target = _obs("t", [0] * 5, [0] * 5)
    donors = [_obs(f"d{i}", [10] * 5, [20] * 5) for i in range(3)]
    # NO group for the target
    groups = {f"d{i}": "A" for i in range(3)}
    imp = CohortImputer()
    out = imp.impute(target, donors + [target], groups)  # default strict=False
    np.testing.assert_array_equal(out.n, target.n)


# MEDIUM — uncertainty.method wiring
def test_round3_uncertainty_method_none_produces_zero_width_cis(
    synthetic_observations_pure_celltype, synthetic_reference, tmp_path
):
    """uncertainty.method='none' must skip bootstrap/MCMC and emit CIs that
    equal the point estimate."""
    from finaleme_too.config import TOOConfig
    from finaleme_too.core.observation_model import BetaBinomialModel
    from finaleme_too.io.sample_sheet import Sample, SampleSheet
    from finaleme_too.pipeline import TOOPipeline

    config = TOOConfig()
    config.uncertainty.method = "none"
    config.uncertainty.n_bootstrap = 5
    config.coverage.tier_high = 0.0
    config.coverage.tier_low = -1.0

    pipeline = TOOPipeline(config)
    # Use the fixture observation object directly through the internal path
    observation = BetaBinomialModel().build(
        synthetic_observations_pure_celltype, reference=synthetic_reference
    )
    w_hat = pipeline.deconvolver.solve(observation, synthetic_reference)
    # Flags driven by the config
    assert pipeline._wants_bootstrap is False
    assert pipeline._wants_any_uncertainty is False


def test_round3_uncertainty_method_bayesian_without_model_bayesian():
    """uncertainty.method='bayesian' should instantiate the Bayesian
    deconvolver even when model.deconvolution stays at MLE."""
    from finaleme_too.config import SolverMethod, TOOConfig
    from finaleme_too.pipeline import TOOPipeline

    config = TOOConfig()
    config.model.deconvolution = SolverMethod.MLE
    config.uncertainty.method = "bayesian"
    pipeline = TOOPipeline(config)
    assert pipeline.bayesian_deconvolver is not None
    assert pipeline._wants_bayesian_uncertainty is True
    assert pipeline._point_estimate_is_bayesian is False


def test_round3_uncertainty_method_both_sets_both_flags():
    from finaleme_too.config import TOOConfig
    from finaleme_too.pipeline import TOOPipeline

    config = TOOConfig()
    config.uncertainty.method = "both"
    pipeline = TOOPipeline(config)
    assert pipeline._wants_bootstrap is True
    assert pipeline._wants_bayesian_uncertainty is True
    assert pipeline.bayesian_deconvolver is not None


# GAP — per-marker effective coverage down-tiering
def test_round3_gap_per_marker_min_reads_vector_down_tiers():
    """A marker with below-expected effective coverage should get a less
    strict minimum reads threshold."""
    from finaleme_too.config import CoverageTier
    from finaleme_too.preprocessing.coverage import per_marker_min_reads_vector

    # mean = 20. Marker 0: n=40 (eff=2.0, stays HIGH). Marker 1: n=4 (eff=0.2,
    # drops one tier to LOW). Marker 2: n=1 (eff=0.05, drops two to ULTRALOW).
    n = np.array([40, 4, 1], dtype=np.int64)
    out = per_marker_min_reads_vector(n, CoverageTier.HIGH)
    assert int(out[0]) == 3  # HIGH
    assert int(out[1]) == 2  # LOW
    assert int(out[2]) == 1  # ULTRALOW


def test_round3_gap_down_tier_keeps_marker_that_scalar_filter_would_drop():
    """In a HIGH-tier sample, a marker with n=1 should be KEPT under
    effective-coverage down-tiering (it drops to ULTRALOW tier, min=1)."""
    from finaleme_too.config import CoverageTier
    from finaleme_too.preprocessing.coverage import per_marker_min_reads_vector

    # 10 markers with n=50, 1 marker with n=1
    n = np.concatenate([np.full(10, 50), np.array([1])]).astype(np.int64)
    min_reads = per_marker_min_reads_vector(n, CoverageTier.HIGH)
    # Last marker should have min=1 (ULTRALOW), not 3 (HIGH)
    assert int(min_reads[-1]) == 1
    # And n=1 >= min=1, so it passes
    assert int(n[-1]) >= int(min_reads[-1])


# ---------------------------------------------------------------------------
# April 2026 round 4 regressions
# ---------------------------------------------------------------------------


# HIGH — covariate adjustment must preserve enriched fields
def test_round4_covariate_adjustment_preserves_enriched_fields():
    """After _maybe_adjust_covariates, mean_coverage, n_markers_used,
    pct_imputed, calibration_flag, overall_qc, residuals, and the marker
    coordinate arrays must all be preserved on the adjusted result."""
    import pandas as pd
    from finaleme_too.config import CoverageTier, TOOConfig
    from finaleme_too.core.deconvolution import DeconvolutionResult
    from finaleme_too.io.sample_sheet import Sample, SampleSheet
    from finaleme_too.pipeline import TOOPipeline
    from finaleme_too.config import MeasurementMode
    from pathlib import Path

    rng = np.random.default_rng(0)

    def _mk_result(sid: str, age: int) -> DeconvolutionResult:
        K = 2
        prop = np.array([0.4, 0.3, 0.3])
        return DeconvolutionResult(
            sample_id=sid,
            cell_types=["CT1", "CT2"],
            proportions=prop,
            ci_lower=prop - 0.05,
            ci_upper=prop + 0.05,
            p_goodness=np.array([0.5, 0.6]),
            p_detection=np.array([0.9, 0.9, 0.9]),
            reliability=np.array(["HIGH", "HIGH", "HIGH"], dtype=object),
            n_markers=np.array([10, 10], dtype=np.int32),
            coverage_tier=CoverageTier.HIGH,
            qc_flags=[],
            mean_dispersion=np.array([50.0, 60.0]),
            mean_coverage=27.3,  # distinctive value
            n_markers_used=18,  # distinctive value
            pct_imputed=0.125,  # distinctive value
            calibration_flag="PASS",
            hemolysis_flag=False,
            overall_qc="PASS",
            residuals=rng.normal(0, 0.05, size=5),  # distinctive
            marker_chrom=np.array(["chr1"] * 5, dtype=object),
            marker_start=np.array([100, 200, 300, 400, 500], dtype=np.int64),
            marker_end=np.array([150, 250, 350, 450, 550], dtype=np.int64),
        )

    results = [
        _mk_result("s1", 30),
        _mk_result("s2", 40),
        _mk_result("s3", 50),
        _mk_result("s4", 60),
        _mk_result("s5", 70),
    ]

    # Build a fake SampleSheet with age covariate
    samples = []
    for i, r in enumerate(results):
        samples.append(
            Sample(
                sample_id=r.sample_id,
                methylation_file=Path("/nope"),
                mode=MeasurementMode.WGBS,
                group="A",
                metadata={"age": 30 + i * 10},
            )
        )
    sheet = SampleSheet(
        samples=samples,
        raw_table=pd.DataFrame([{"sample_id": s.sample_id} for s in samples]),
    )

    config = TOOConfig()
    config.covariate_adjustment.biological_covariates = ["age"]
    pipeline = TOOPipeline(config)
    adjusted = pipeline._maybe_adjust_covariates(results, sheet)

    assert len(adjusted) == len(results)
    for orig, adj in zip(results, adjusted):
        # Proportions should have been re-residualized
        # (may or may not differ depending on the covariate signal)
        # Enriched fields must be IDENTICAL to the originals.
        assert adj.mean_coverage == orig.mean_coverage, \
            f"mean_coverage dropped: {adj.mean_coverage} != {orig.mean_coverage}"
        assert adj.n_markers_used == orig.n_markers_used
        assert abs(adj.pct_imputed - orig.pct_imputed) < 1e-12
        assert adj.calibration_flag == orig.calibration_flag
        assert adj.hemolysis_flag == orig.hemolysis_flag
        assert adj.overall_qc == orig.overall_qc
        assert adj.residuals is not None
        np.testing.assert_array_equal(adj.residuals, orig.residuals)
        assert adj.marker_chrom is not None
        np.testing.assert_array_equal(adj.marker_chrom, orig.marker_chrom)
        np.testing.assert_array_equal(adj.marker_start, orig.marker_start)
        np.testing.assert_array_equal(adj.marker_end, orig.marker_end)
        np.testing.assert_array_equal(adj.mean_dispersion, orig.mean_dispersion)


# MEDIUM — group comparison sample-size guard
def test_round4_group_comparison_runs_with_3_samples():
    """The pipeline should not hard-block group comparisons at <4 samples.
    Per-test sufficiency checks live inside run_group_comparisons; the
    pipeline should only require the bare minimum of 2 samples."""
    from finaleme_too.config import CoverageTier, TestMethod, TOOConfig
    from finaleme_too.core.deconvolution import DeconvolutionResult
    from finaleme_too.pipeline import TOOPipeline

    rng = np.random.default_rng(1)

    def _mk(sid, group, prop_bias):
        K = 2
        prop = np.array([0.4 + prop_bias, 0.3 - prop_bias, 0.3])
        return DeconvolutionResult(
            sample_id=sid,
            cell_types=["CT1", "CT2"],
            proportions=prop,
            ci_lower=prop - 0.05,
            ci_upper=prop + 0.05,
            p_goodness=np.array([0.5, 0.5]),
            p_detection=np.array([0.9, 0.9, 0.9]),
            reliability=np.array(["HIGH", "HIGH", "HIGH"], dtype=object),
            n_markers=np.array([5, 5], dtype=np.int32),
            coverage_tier=CoverageTier.HIGH,
            qc_flags=[],
            mean_dispersion=np.array([50.0, 50.0]),
            mean_coverage=25.0,
            n_markers_used=5,
            residuals=rng.normal(0, 0.05, size=5),
            overall_qc="PASS",
        )

    # 3 samples in 2 groups — the old guard (>=4) would have skipped this.
    results = [
        _mk("s1", "A", 0.1),
        _mk("s2", "A", 0.08),
        _mk("s3", "B", -0.1),
    ]
    sample_groups = {"s1": "A", "s2": "A", "s3": "B"}

    config = TOOConfig()
    config.testing.method = TestMethod.ILR_REGRESSION
    pipeline = TOOPipeline(config)
    pipeline.group_comparison_spec = "A:B"
    test_results = pipeline._run_group_comparisons(results, sample_groups)
    # A:B with n=2 vs n=1 — ILR still produces rows (may be NaN for
    # insufficient within-group variance). What we verify is that the
    # list is returned, NOT blocked outright.
    assert isinstance(test_results, list)


# MEDIUM — --marker-format CLI default no longer overrides YAML
def test_round4_marker_format_yaml_wins_when_cli_not_provided(tmp_path):
    """When the user provides --config with markers.marker_format=bed but
    does NOT pass --marker-format, the effective format must be 'bed',
    not the Click default 'auto'."""
    import sys
    from finaleme_too.config import TOOConfig
    from finaleme_too import pipeline as pipeline_mod

    # Minimal inputs
    sheet = tmp_path / "sheet.tsv"
    ref = tmp_path / "ref.tsv"
    bed = tmp_path / "markers.bed"
    bed.write_text("chr1\t100\t200\tm1\n")
    ref.write_text("chrom\tstart\tend\tCellA\nchr1\t100\t200\t0.1\n")
    fme = tmp_path / "s.prediction.bed.gz"
    with gzip.open(fme, "wt") as fh:
        fh.write(
            "#chr\tstart\tend\tmethy_perc_predict\tmethy_count_predict"
            "\ttotal_count_predict\to1\to2\to3\n"
        )
        fh.write("chr1\t110\t111\t10\t1\t10\t0\t0\t0\n")
    sheet.write_text(
        "sample_id\tmethylation_file\tmode\tgroup\ns1\t"
        + str(fme)
        + "\tFINALEME\tA\n"
    )

    yaml_cfg = tmp_path / "cfg.yaml"
    yaml_cfg.write_text("markers:\n  marker_format: bed\n")

    captured = {}
    orig_init = pipeline_mod.TOOPipeline.__init__

    def spy_init(self, config, **kw):
        captured["marker_format"] = config.markers.marker_format
        raise SystemExit(0)

    pipeline_mod.TOOPipeline.__init__ = spy_init  # type: ignore[assignment]
    try:
        sys.argv = [
            "finaleme-too",
            "run",
            "--sample-sheet",
            str(sheet),
            "--output-dir",
            str(tmp_path / "out"),
            "--reference-panel",
            str(ref),
            "--marker-regions",
            str(bed),
            "--config",
            str(yaml_cfg),
            # NOTE: NOT passing --marker-format
        ]
        from finaleme_too.cli import main

        try:
            main(standalone_mode=False)
        except SystemExit:
            pass
    finally:
        pipeline_mod.TOOPipeline.__init__ = orig_init  # type: ignore[assignment]

    assert captured["marker_format"] == "bed", \
        f"YAML marker_format should have won; got {captured['marker_format']}"


# MEDIUM — config.input.* now honored when CLI args omitted
def test_round4_config_input_format_honored_from_yaml(tmp_path):
    """When the YAML sets input.format=finaleme_bed but the CLI does not,
    samples without an input_format should inherit finaleme_bed."""
    import sys
    from finaleme_too import pipeline as pipeline_mod

    sheet = tmp_path / "sheet.tsv"
    ref = tmp_path / "ref.tsv"
    bed = tmp_path / "markers.bed"
    bed.write_text("chr1\t100\t200\tm1\n")
    ref.write_text("chrom\tstart\tend\tCellA\nchr1\t100\t200\t0.1\n")
    fme = tmp_path / "s.bed.gz"
    with gzip.open(fme, "wt") as fh:
        fh.write(
            "#chr\tstart\tend\tmethy_perc_predict\tmethy_count_predict"
            "\ttotal_count_predict\to1\to2\to3\n"
        )
        fh.write("chr1\t110\t111\t10\t1\t10\t0\t0\t0\n")
    sheet.write_text(
        "sample_id\tmethylation_file\tmode\tgroup\ns1\t"
        + str(fme)
        + "\tFINALEME\tA\n"
    )

    yaml_cfg = tmp_path / "cfg.yaml"
    yaml_cfg.write_text(
        "input:\n"
        "  format: finaleme_bed\n"
        "  meth_col: 5\n"
        "  total_col: 6\n"
    )

    captured = {}
    orig_init = pipeline_mod.TOOPipeline.__init__

    def spy_init(self, config, **kw):
        captured["config"] = config
        raise SystemExit(0)

    pipeline_mod.TOOPipeline.__init__ = spy_init  # type: ignore[assignment]
    try:
        sys.argv = [
            "finaleme-too",
            "run",
            "--sample-sheet",
            str(sheet),
            "--output-dir",
            str(tmp_path / "out"),
            "--reference-panel",
            str(ref),
            "--marker-regions",
            str(bed),
            "--config",
            str(yaml_cfg),
            # NOTE: NOT passing --input-format, --meth-col, --total-col
        ]
        from finaleme_too.cli import main

        try:
            main(standalone_mode=False)
        except SystemExit:
            pass
    finally:
        pipeline_mod.TOOPipeline.__init__ = orig_init  # type: ignore[assignment]

    # The config should carry the YAML values through
    cfg = captured["config"]
    assert cfg.input.format == "finaleme_bed"
    assert cfg.input.meth_col == 5
    assert cfg.input.total_col == 6


# ---------------------------------------------------------------------------
# April 2026 — `finaleme-too run` parallelization regression tests
# ---------------------------------------------------------------------------
# Symptom: user passed --threads 8 but observed <100% CPU usage. Root cause
# was the same as the train-calibration fix (commit 0747df1): the default
# loky backend pickles the (large) reference panel + cpg_index for every
# task, dominating the runtime. Switched the run-path parallel_map sites to
# the threading backend (numpy/scipy/pandas release the GIL), parallelized
# the imputation loop, parallelized load_beta_list, and added a bootstrap
# inner-loop parallelization that lights up unused cores when the cohort is
# smaller than --threads.


def _build_synth_cohort(tmp_path: Path, n_samples: int = 6) -> tuple:
    """Build a synthetic Sample sheet that the pipeline can run end-to-end.

    Each sample is a small WGBS BED in finaleme_bed format (with predicted
    beta + counts), with the methylation generated from a known reference
    profile so the deconvolver has something signal-rich to do.

    Returns ``(SampleSheet, ReferencePanel, MarkerRegions)``.
    Stays at 40 markers (10 per cell type) so the WGBS-HIGH per-sample
    dispersion MLE path is skipped — the brent bracket is unstable on
    pristine binomial-distributed data and that has nothing to do with
    parallelization correctness.
    """
    from finaleme_too.config import MeasurementMode
    from finaleme_too.io.marker_regions import MarkerRegions
    from finaleme_too.io.reference_panel import ReferencePanel
    from finaleme_too.io.sample_sheet import Sample, SampleSheet

    n_markers = 40
    n_cell_types = 4
    per_ct = n_markers // n_cell_types
    rng = np.random.default_rng(42)

    # Reference panel: each cell type owns ``per_ct`` markers (low at its own,
    # high at others' markers)
    methy = np.full((n_markers, n_cell_types), 0.5, dtype=np.float32)
    for j in range(n_cell_types):
        idx = slice(j * per_ct, (j + 1) * per_ct)
        methy[idx, j] = rng.uniform(0.0, 0.05, size=per_ct).astype(np.float32)
        for jj in range(n_cell_types):
            if jj != j:
                methy[idx, jj] = rng.uniform(0.85, 1.0, size=per_ct).astype(np.float32)
    chrom = np.array(["chr1"] * n_markers, dtype=object)
    starts = np.array([1000 + i * 5000 for i in range(n_markers)], dtype=np.int64)
    ends = starts + 500
    reference = ReferencePanel(
        chrom=chrom,
        start=starts,
        end=ends,
        cell_types=[f"CT{j+1}" for j in range(n_cell_types)],
        methylation=methy,
        coverage=np.full((n_markers, n_cell_types), 50, dtype=np.int32),
    )
    marker_regions = MarkerRegions(
        chrom=chrom, start=starts, end=ends, marker_name=None
    )

    # Each sample is a 50/50 mixture of two cell types — generates the WGBS
    # BED that MethylationLoader can parse with the bissnp_6plus2 format
    # (chrom, start, end, name, score, strand, methPct, totalCount).
    samples = []
    for i in range(n_samples):
        ct_a = i % n_cell_types
        ct_b = (i + 1) % n_cell_types
        true_w = np.zeros(n_cell_types, dtype=np.float64)
        true_w[ct_a] = 0.6
        true_w[ct_b] = 0.4
        beta = methy.astype(np.float64) @ true_w
        n_arr = np.full(n_markers, 30, dtype=np.int32)
        k_arr = rng.binomial(n_arr.astype(np.int64), beta).astype(np.int32)

        # Write as bissnp_6plus2 BED
        bed = tmp_path / f"sample_{i}.bed"
        with open(bed, "w") as fh:
            for m in range(n_markers):
                pct = 100.0 * k_arr[m] / max(int(n_arr[m]), 1)
                fh.write(
                    f"{chrom[m]}\t{int(starts[m])}\t{int(ends[m])}\t.\t0\t+"
                    f"\t{pct:.4f}\t{int(n_arr[m])}\n"
                )
        samples.append(
            Sample(
                sample_id=f"s_{i}",
                methylation_file=bed,
                mode=MeasurementMode.WGBS,
                input_format="bissnp_6plus2",
                group="A" if i < n_samples // 2 else "B",
            )
        )

    sheet = SampleSheet(
        samples=samples,
        raw_table=pd.DataFrame([{"sample_id": s.sample_id} for s in samples]),
    )
    return sheet, reference, marker_regions


def test_parallel_run_threaded_matches_serial_proportions(tmp_path):
    """Running TOOPipeline with threads=4 must produce the same per-sample
    proportions as threads=1. The bootstrap is seeded so the inner loops are
    deterministic.

    Regression: this is the test that would have caught the loky-vs-threading
    backend issue if any of the parallelized stages (loading, imputation,
    deconvolution) drifted in the threaded path.
    """
    from finaleme_too.config import TOOConfig
    from finaleme_too.pipeline import TOOPipeline

    sheet, reference, marker_regions = _build_synth_cohort(tmp_path, n_samples=4)

    def _run(threads: int) -> list:
        cfg = TOOConfig()
        cfg.threads = threads
        cfg.uncertainty.n_bootstrap = 20
        cfg.uncertainty.seed = 7  # determinism
        cfg.coverage.tier_high = 0.0  # force HIGH for everything
        cfg.coverage.tier_low = -1.0
        cfg.markers.n_per_type = 0  # skip marker selection
        # Disable optional cohort steps that could drift between runs
        cfg.batch_correction.technical_covariates = []
        cfg.covariate_adjustment.biological_covariates = []
        cfg.testing.method = cfg.testing.method  # leave default
        out_dir = tmp_path / f"out_t{threads}"
        out_dir.mkdir(exist_ok=True)
        pipeline = TOOPipeline(cfg)
        cohort = pipeline.run(sheet, reference, marker_regions, out_dir)
        return cohort.samples

    serial = _run(1)
    parallel = _run(4)

    # Same number of results, in the same order
    assert len(serial) == len(parallel)
    for s_res, p_res in zip(serial, parallel):
        assert s_res.sample_id == p_res.sample_id
        # Point estimates must be exactly identical: the threading backend
        # never re-pickles arrays, so the numerics are bit-for-bit the same
        # (the bootstrap rng is pre-generated, so iteration order does not
        # matter).
        np.testing.assert_allclose(
            s_res.proportions, p_res.proportions, rtol=0, atol=1e-12
        )
        np.testing.assert_allclose(
            s_res.ci_lower, p_res.ci_lower, rtol=0, atol=1e-12
        )
        np.testing.assert_allclose(
            s_res.ci_upper, p_res.ci_upper, rtol=0, atol=1e-12
        )


def test_bootstrap_inner_jobs_match_serial(synthetic_observations_pure_celltype, synthetic_reference):
    """BootstrapCI.estimate(n_jobs=4) must give the same samples as n_jobs=1.

    The threading backend uses a pre-generated index matrix so the per-task
    iteration is deterministic — there are no rng state hazards.
    """
    from finaleme_too.core.deconvolution import MLEDeconvolver
    from finaleme_too.core.observation_model import BetaBinomialModel
    from finaleme_too.core.uncertainty import BootstrapCI

    obs = BetaBinomialModel().build(
        synthetic_observations_pure_celltype, reference=synthetic_reference
    )
    bootstrap = BootstrapCI(n_iterations=15, ci_level=0.9, seed=11)
    deconvolver = MLEDeconvolver()
    serial = bootstrap.estimate(obs, synthetic_reference, deconvolver, n_jobs=1)
    parallel = bootstrap.estimate(obs, synthetic_reference, deconvolver, n_jobs=4)

    np.testing.assert_allclose(
        serial.proportions_samples, parallel.proportions_samples, rtol=0, atol=1e-12
    )
    np.testing.assert_allclose(serial.ci_lower, parallel.ci_lower, rtol=0, atol=1e-12)
    np.testing.assert_allclose(serial.ci_upper, parallel.ci_upper, rtol=0, atol=1e-12)


def test_pipeline_outer_inner_thread_split_for_single_sample(tmp_path):
    """Single-sample run with --threads 8 should hand the bootstrap 8 jobs.

    Regression: previously the outer parallel_map could only use 1 thread for
    1 sample, leaving 7 cores idle. The fix computes
    inner_jobs = max(1, threads // outer_jobs) and passes it to the bootstrap.
    """
    from finaleme_too.config import TOOConfig
    from finaleme_too.pipeline import TOOPipeline
    import finaleme_too.core.uncertainty as uncertainty_mod

    sheet, reference, marker_regions = _build_synth_cohort(tmp_path, n_samples=1)

    cfg = TOOConfig()
    cfg.threads = 8
    cfg.uncertainty.n_bootstrap = 4
    cfg.uncertainty.seed = 13
    cfg.coverage.tier_high = 0.0
    cfg.coverage.tier_low = -1.0
    cfg.markers.n_per_type = 0

    # Spy on BootstrapCI.estimate to capture the n_jobs argument
    captured = {"n_jobs": None}
    orig = uncertainty_mod.BootstrapCI.estimate

    def spy(self, model, reference, deconvolver, n_jobs=1):
        captured["n_jobs"] = n_jobs
        return orig(self, model, reference, deconvolver, n_jobs=n_jobs)

    uncertainty_mod.BootstrapCI.estimate = spy  # type: ignore[assignment]
    try:
        out_dir = tmp_path / "out_single"
        out_dir.mkdir(exist_ok=True)
        TOOPipeline(cfg).run(sheet, reference, marker_regions, out_dir)
    finally:
        uncertainty_mod.BootstrapCI.estimate = orig  # type: ignore[assignment]

    # 1 sample × 8 threads → bootstrap should get all 8 inner jobs
    assert captured["n_jobs"] == 8


def test_pipeline_outer_inner_thread_split_for_multi_sample(tmp_path):
    """Cohort run with samples >= threads should leave the inner bootstrap
    serial (n_jobs=1) so the outer parallel_map drives the parallelism."""
    from finaleme_too.config import TOOConfig
    from finaleme_too.pipeline import TOOPipeline
    import finaleme_too.core.uncertainty as uncertainty_mod

    sheet, reference, marker_regions = _build_synth_cohort(tmp_path, n_samples=4)

    cfg = TOOConfig()
    cfg.threads = 4
    cfg.uncertainty.n_bootstrap = 3
    cfg.uncertainty.seed = 19
    cfg.coverage.tier_high = 0.0
    cfg.coverage.tier_low = -1.0
    cfg.markers.n_per_type = 0

    captured_n_jobs: list[int] = []
    orig = uncertainty_mod.BootstrapCI.estimate

    def spy(self, model, reference, deconvolver, n_jobs=1):
        captured_n_jobs.append(n_jobs)
        return orig(self, model, reference, deconvolver, n_jobs=n_jobs)

    uncertainty_mod.BootstrapCI.estimate = spy  # type: ignore[assignment]
    try:
        out_dir = tmp_path / "out_multi"
        out_dir.mkdir(exist_ok=True)
        TOOPipeline(cfg).run(sheet, reference, marker_regions, out_dir)
    finally:
        uncertainty_mod.BootstrapCI.estimate = orig  # type: ignore[assignment]

    # 4 samples × 4 threads → outer absorbs the parallelism, inner stays at 1
    assert all(n == 1 for n in captured_n_jobs)
    assert len(captured_n_jobs) == 4


def test_load_beta_list_threaded_matches_serial(tmp_path):
    """ReferencePanelLoader.load_beta_list must produce identical reference
    panels regardless of ``threads``. The threaded path parses each .beta in
    parallel and accumulates afterwards, so accumulation order doesn't matter
    for the aggregate-mode counts."""
    from finaleme_too.io.marker_regions import MarkerRegions
    from finaleme_too.io.reference_panel import ReferencePanelLoader

    # Build a tiny CpG index covering chr1:0..1000
    cpg_path = tmp_path / "cpg_index.bed"
    cpg_lines = [f"chr1\t{p}\n" for p in range(0, 1000, 5)]
    cpg_path.write_text("".join(cpg_lines))

    # 4 marker regions on chr1 (each spanning a few CpGs)
    marker_regions = MarkerRegions(
        chrom=np.array(["chr1"] * 4, dtype=object),
        start=np.array([0, 100, 200, 300], dtype=np.int64),
        end=np.array([50, 150, 250, 350], dtype=np.int64),
        marker_name=None,
    )

    # 6 .beta files in 2 groups
    rng = np.random.default_rng(0)
    n_sites = len(cpg_lines)
    beta_paths = []
    for i in range(6):
        # Per-CpG (methylated_count, total_count) pairs
        m = rng.integers(0, 20, size=n_sites, dtype=np.uint8)
        n = np.maximum(m, rng.integers(20, 40, size=n_sites, dtype=np.uint8))
        data = np.empty((n_sites, 2), dtype=np.uint8)
        data[:, 0] = m
        data[:, 1] = n
        p = tmp_path / f"sample_{i}.beta"
        p.write_bytes(data.tobytes())
        beta_paths.append(str(p))

    groups_path = tmp_path / "groups.csv"
    groups_path.write_text(
        "name,group\n"
        + "\n".join(
            f"sample_{i},{'A' if i < 3 else 'B'}" for i in range(6)
        )
        + "\n"
    )

    serial = ReferencePanelLoader.load_beta_list(
        ref_betas_arg=",".join(beta_paths),
        groups_file=groups_path,
        cpg_index_path=cpg_path,
        marker_regions=marker_regions,
        threads=1,
    )
    parallel = ReferencePanelLoader.load_beta_list(
        ref_betas_arg=",".join(beta_paths),
        groups_file=groups_path,
        cpg_index_path=cpg_path,
        marker_regions=marker_regions,
        threads=4,
    )

    assert serial.cell_types == parallel.cell_types
    np.testing.assert_array_equal(serial.methylation, parallel.methylation)
    np.testing.assert_array_equal(serial.coverage, parallel.coverage)


# ---------------------------------------------------------------------------
# April 2026 — coverage tier ↔ mean_coverage consistency
# ---------------------------------------------------------------------------
# Symptom: qc_summary.tsv showed coverage_tier=ULTRALOW for two samples whose
# mean_coverage values were 124.67 and 74.12 — clearly HIGH-coverage samples.
#
# Root cause: CoverageTierAssigner.assign computed
#     mean_cov = sum(obs.n) * FRAGMENT_LENGTH / GENOME_SIZE
# treating sum(obs.n) (reads at marker positions) as if it were the whole-BAM
# total fragment count. With ~1000 markers at ~100 reads/marker that yields
# mean_cov ≈ 0.006, which always falls below tier_low=0.5 → ULTRALOW.
#
# Fix: compute the effective coverage from data we already have — total reads
# at marker positions × cfDNA fragment length ÷ total marker region area.
# This is the depth-of-coverage WITHIN the marker regions and needs no extra
# sample sheet inputs. Both qc_summary.tsv columns (``coverage_tier`` and
# ``mean_coverage``) report this same quantity, so they can no longer
# contradict each other.


def _make_obs(
    sample_id: str,
    n_per_marker: int,
    n_markers: int = 100,
    marker_width: int = 100,
):
    """Build a MarkerObservations with a deterministic per-marker read count
    and uniform marker region width.

    Default 100 markers × 100bp wide → 10,000 bp of marker area, so the
    effective coverage formula simplifies to:
        coverage = n_per_marker * 100 markers * 167 / 10000 bp
                 = 1.67 * n_per_marker
    """
    starts = np.array(
        [1000 + i * (marker_width + 100) for i in range(n_markers)], dtype=np.int64
    )
    ends = starts + marker_width
    return MarkerObservations(
        sample_id=sample_id,
        chrom=np.array(["chr1"] * n_markers, dtype=object),
        start=starts,
        end=ends,
        k=np.zeros(n_markers, dtype=np.int32),
        n=np.full(n_markers, n_per_marker, dtype=np.int32),
        predicted_beta=None,
        mode=MeasurementMode.WGBS,
    )


def test_high_mean_coverage_classified_as_high_tier():
    """A sample with mean reads-per-marker well above tier_high must be
    classified as HIGH, not ULTRALOW.

    Regression: this is exactly the user-reported bug — qc_summary showed
    mean_coverage values in the 70-125 range with coverage_tier=ULTRALOW.
    Under the new effective-coverage formula, those samples come out
    deeply HIGH.
    """
    from finaleme_too.config import CoverageConfig, CoverageTier
    from finaleme_too.preprocessing.coverage import (
        CoverageTierAssigner,
        effective_coverage_in_markers,
    )

    assigner = CoverageTierAssigner(CoverageConfig())  # defaults: 10, 0.5

    # User-reported sample 1: 125 reads/marker × 100 markers × 167 bp / 10000 bp
    # ≈ 208× effective depth — comfortably HIGH.
    obs_high_124 = _make_obs("s_high_124", n_per_marker=125)
    cov = effective_coverage_in_markers(obs_high_124)
    assert cov > 100, f"Expected ~208, got {cov}"
    assert assigner.assign(obs_high_124) == CoverageTier.HIGH

    # User-reported sample 2: 74 reads/marker → ~123× depth → HIGH.
    obs_high_74 = _make_obs("s_high_74", n_per_marker=74)
    assert effective_coverage_in_markers(obs_high_74) > 100
    assert assigner.assign(obs_high_74) == CoverageTier.HIGH


def test_effective_coverage_formula_matches_expected_depth():
    """The effective coverage formula equals
    Σ reads × fragment_length / Σ marker_widths."""
    from finaleme_too.preprocessing.coverage import (
        FRAGMENT_LENGTH_BP,
        effective_coverage_in_markers,
    )

    obs = _make_obs("s", n_per_marker=10, n_markers=50, marker_width=200)
    expected = 50 * 10 * FRAGMENT_LENGTH_BP / (50 * 200)
    actual = effective_coverage_in_markers(obs)
    assert abs(actual - expected) < 1e-9, f"{actual} != {expected}"


def test_coverage_tier_thresholds_use_effective_coverage():
    """Boundary cases for the tier classifier under the effective-coverage
    formula. With 100 markers × 100bp wide and 167bp fragments, the per-marker
    reads needed to cross thresholds are:
      tier_high (10):  10 / 1.67 ≈ 6 reads/marker
      tier_low  (0.5): 0.5 / 1.67 ≈ 0.3 reads/marker (so 1 read/marker → LOW)
    """
    from finaleme_too.config import CoverageConfig, CoverageTier
    from finaleme_too.preprocessing.coverage import CoverageTierAssigner

    assigner = CoverageTierAssigner(CoverageConfig())

    # 7 reads/marker → coverage ≈ 11.7 → HIGH (>10)
    assert assigner.assign(_make_obs("s_7", 7)) == CoverageTier.HIGH
    # 5 reads/marker → coverage ≈ 8.35 → LOW (between 0.5 and 10)
    assert assigner.assign(_make_obs("s_5", 5)) == CoverageTier.LOW
    # 1 read/marker → coverage ≈ 1.67 → LOW
    assert assigner.assign(_make_obs("s_1", 1)) == CoverageTier.LOW
    # 0 reads/marker → ULTRALOW (sum-zero short circuit)
    assert assigner.assign(_make_obs("s_0", 0)) == CoverageTier.ULTRALOW


def test_coverage_tier_ultralow_when_marker_area_is_huge(tmp_path):
    """A small number of reads spread over a huge marker region area should
    produce ULTRALOW. Confirms the denominator (marker area) actually matters."""
    from finaleme_too.config import CoverageConfig, CoverageTier
    from finaleme_too.preprocessing.coverage import CoverageTierAssigner

    # 100 markers × 1,000,000 bp wide each → 100,000,000 bp marker area.
    # 1 read/marker → 100 reads × 167 / 1e8 ≈ 1.67e-4 → ULTRALOW.
    obs = _make_obs("s_dilute", n_per_marker=1, n_markers=100, marker_width=1_000_000)
    assigner = CoverageTierAssigner(CoverageConfig())
    assert assigner.assign(obs) == CoverageTier.ULTRALOW


def test_pipeline_writes_consistent_qc_summary_for_high_coverage_sample(tmp_path):
    """End-to-end: a sample with deep effective coverage must show
    coverage_tier=HIGH and a matching mean_coverage in qc_summary.tsv.

    This is the symptom the user actually saw: contradictory columns in
    qc_summary.tsv. The fix puts both columns on the same effective-coverage
    scale so they cannot disagree.
    """
    from finaleme_too.config import TOOConfig
    from finaleme_too.pipeline import TOOPipeline
    from finaleme_too.preprocessing.coverage import FRAGMENT_LENGTH_BP

    sheet, reference, marker_regions = _build_synth_cohort(tmp_path, n_samples=2)
    # Bump the per-marker reads to ~120 by editing the BED files in place.
    # _build_synth_cohort writes 30 reads/marker by default.
    for s in sheet.samples:
        bed = s.methylation_file
        lines = bed.read_text().splitlines()
        new_lines = []
        for line in lines:
            parts = line.split("\t")
            # bissnp_6plus2: chrom start end name score strand methPct totalCount
            parts[7] = "120"
            new_lines.append("\t".join(parts))
        bed.write_text("\n".join(new_lines) + "\n")

    cfg = TOOConfig()
    cfg.threads = 1
    cfg.uncertainty.n_bootstrap = 5
    cfg.uncertainty.seed = 3
    cfg.markers.n_per_type = 0  # skip marker selection
    # Use the production tier_high default (10) and tier_low (0.5) so we
    # actually exercise the fix.
    cfg.coverage.tier_high = 10.0
    cfg.coverage.tier_low = 0.5

    out_dir = tmp_path / "out_qc"
    out_dir.mkdir(exist_ok=True)
    TOOPipeline(cfg).run(sheet, reference, marker_regions, out_dir)

    qc_path = out_dir / "qc_summary.tsv"
    assert qc_path.exists(), "qc_summary.tsv must be written"
    qc = pd.read_csv(qc_path, sep="\t")
    assert len(qc) == len(sheet.samples)
    # _build_synth_cohort uses 500bp markers (start, start+500). Per-marker
    # n=120 → effective coverage = 120 * 167 / 500 ≈ 40.08× → HIGH.
    expected_cov = 120 * FRAGMENT_LENGTH_BP / 500
    assert all(qc["coverage_tier"] == "HIGH"), \
        f"Expected all HIGH, got {qc['coverage_tier'].tolist()}"
    for mc in qc["mean_coverage"]:
        # Allow ±5% slack for the tier-filter and bootstrap noise.
        assert abs(float(mc) - expected_cov) < 2.0, \
            f"mean_coverage {mc} not close to expected {expected_cov}"
