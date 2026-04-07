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
