"""Tests for calibration apply path and calibration_eval QC (Phase B)."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from finaleme_too.config import MeasurementMode
from finaleme_too.io.methylation_loader import MarkerObservations
from finaleme_too.preprocessing.calibration import (
    CalibrationParams,
    apply_calibration,
    load_default_calibration,
)
from finaleme_too.preprocessing.calibration_eval import (
    compute_hosmer_lemeshow,
    compute_inference_qc,
)


def _identity_calibration(n_bins: int = 4) -> CalibrationParams:
    return CalibrationParams(
        n_bins=n_bins,
        bin_edges=np.linspace(0.0, 1.0, n_bins + 1),
        a=np.ones(n_bins),
        c=np.zeros(n_bins),
        log_dispersion=np.full(n_bins, np.log(20.0)),
    )


def test_calibration_save_load_round_trip(tmp_path):
    p = _identity_calibration(6)
    out = tmp_path / "cal.json"
    p.save(out)
    loaded = CalibrationParams.load(out)
    assert loaded.n_bins == 6
    np.testing.assert_allclose(loaded.bin_edges, p.bin_edges)
    np.testing.assert_allclose(loaded.a, p.a)
    np.testing.assert_allclose(loaded.c, p.c)


def test_apply_identity_calibration_is_no_op(synthetic_marker_regions):
    n = synthetic_marker_regions.n_markers
    pred = np.full(n, 0.4, dtype=np.float32)
    obs = MarkerObservations(
        sample_id="x",
        chrom=synthetic_marker_regions.chrom,
        start=synthetic_marker_regions.start,
        end=synthetic_marker_regions.end,
        k=np.full(n, 12, dtype=np.int32),
        n=np.full(n, 30, dtype=np.int32),
        predicted_beta=pred,
        mode=MeasurementMode.FINALEME,
    )
    cal = _identity_calibration()
    # No region annotations → defaults to middle bin (still identity)
    new_obs = apply_calibration(obs, cal, region_annotations=None)
    # Identity calibration should preserve k roughly equal to round(0.4 * 30) = 12
    assert int(new_obs.k[0]) == 12
    np.testing.assert_array_equal(new_obs.n, obs.n)


def test_apply_calibration_changes_methylation():
    """Calibration with a > 1 inflates extreme predictions."""
    pred = np.array([0.2, 0.8], dtype=np.float32)
    obs = MarkerObservations(
        sample_id="x",
        chrom=np.array(["chr1", "chr1"], dtype=object),
        start=np.array([100, 200], dtype=np.int64),
        end=np.array([200, 300], dtype=np.int64),
        k=np.array([6, 24], dtype=np.int32),
        n=np.array([30, 30], dtype=np.int32),
        predicted_beta=pred,
        mode=MeasurementMode.FINALEME,
    )
    # a = 1.5 stretches logit; c = 0
    cal = CalibrationParams(
        n_bins=1,
        bin_edges=np.array([0.0, 1.0]),
        a=np.array([1.5]),
        c=np.array([0.0]),
        log_dispersion=np.array([np.log(20.0)]),
    )
    new_obs = apply_calibration(obs, cal, region_annotations=None)
    # 0.2 should move toward 0; 0.8 should move toward 1
    assert int(new_obs.k[0]) < 6
    assert int(new_obs.k[1]) > 24


def test_default_calibration_is_loadable():
    cal = load_default_calibration()
    assert cal.n_bins == 8
    assert cal.bin_edges.size == 9
    assert cal.a.size == 8
    assert cal.c.size == 8
    assert cal.log_dispersion.size == 8


def test_hosmer_lemeshow_perfect_calibration():
    """Perfectly calibrated predictions should give a high p-value."""
    rng = np.random.default_rng(0)
    pred = np.linspace(0.05, 0.95, 200)
    obs = rng.binomial(50, pred) / 50
    p = compute_hosmer_lemeshow(obs, pred, n_groups=10)
    assert not np.isnan(p)
    assert p > 0.01  # not strictly significant under perfect calibration


def test_hosmer_lemeshow_miscalibrated():
    """Severely miscalibrated predictions should give a low p-value."""
    rng = np.random.default_rng(1)
    pred = np.full(500, 0.5)
    obs = rng.binomial(50, 0.9, size=500) / 50  # truly 0.9, predicted 0.5
    p = compute_hosmer_lemeshow(obs, pred, n_groups=10)
    assert p < 0.001


def test_inference_qc_pass_for_normal_sample():
    cal = _identity_calibration(8)
    sample_pred = np.linspace(0.05, 0.95, 200)
    qc = compute_inference_qc(sample_pred, cal)
    assert qc["flag"] in ("PASS", "WARN", "FAIL")
    assert 0.0 <= qc["bin_coverage_balance"] <= 1.0
    assert 0.0 <= qc["prediction_range_coverage"] <= 1.0


# ---------------------------------------------------------------------------
# Phase C: training-time fit + tuning
# ---------------------------------------------------------------------------


def test_fit_calibration_recovers_identity():
    from finaleme_too.preprocessing.calibration_eval import fit_calibration

    rng = np.random.default_rng(0)
    n = 500
    truth = rng.beta(2, 5, size=n)
    fme = truth + rng.normal(0, 0.02, size=n)
    fme = np.clip(fme, 0.01, 0.99)
    density = rng.uniform(0, 0.1, size=n)
    fit = fit_calibration(
        finaleme_beta=fme,
        wgbs_beta=truth,
        cpg_density=density,
        n_bins=4,
    )
    # Identity-like → slope ~ 1, intercept ~ 0
    np.testing.assert_allclose(fit.a, np.ones(4), atol=0.2)
    np.testing.assert_allclose(fit.c, np.zeros(4), atol=0.2)
    assert fit.overall["r_squared"] > 0.9


def test_tune_n_bins_selects_finite_value():
    from finaleme_too.preprocessing.calibration_eval import tune_n_bins

    rng = np.random.default_rng(1)
    n = 600
    sample_ids = np.array([f"S{i // 100}" for i in range(n)])
    truth = rng.beta(2, 5, size=n)
    fme = np.clip(truth + rng.normal(0, 0.05, size=n), 0.01, 0.99)
    density = rng.uniform(0, 0.1, size=n)
    out = tune_n_bins(
        finaleme_beta=fme,
        wgbs_beta=truth,
        cpg_density=density,
        sample_ids=sample_ids,
        n_bins_candidates=[2, 4, 8],
    )
    assert out["selected_n_bins"] in (2, 4, 8)
    assert len(out["candidates"]) == 3


def test_train_calibration_end_to_end(tmp_path):
    from scipy.special import expit, logit

    from finaleme_too.preprocessing.calibration import train_calibration

    rng = np.random.default_rng(2)
    n_samples = 4
    n_markers = 150
    rows_w = []
    rows_f = []
    rows_a = []
    for s in range(n_samples):
        sid = f"S{s}"
        for m in range(n_markers):
            true_b = float(rng.beta(2, 5))
            n_w = 30
            k_w = int(rng.binomial(n_w, true_b))
            density = float(rng.uniform(0, 0.1))
            # Miscalibration: a=0.8, c=0.1
            true_logit = logit(np.clip(true_b, 1e-6, 1 - 1e-6))
            fme_logit = (true_logit - 0.1) / 0.8 + rng.normal(0, 0.2)
            fme_b = float(expit(fme_logit))
            n_f = 30
            k_f = int(rng.binomial(n_f, fme_b))
            rows_w.append((sid, "chr1", m * 100, m * 100 + 50, k_w, n_w))
            rows_f.append((sid, "chr1", m * 100, m * 100 + 50, k_f, n_f))
            if s == 0:
                rows_a.append(("chr1", m * 100, m * 100 + 50, density))

    import pandas as pd

    pd.DataFrame(
        rows_w, columns=["sample_id", "chrom", "start", "end", "methylated_count", "total_count"]
    ).to_csv(tmp_path / "wgbs.tsv", sep="\t", index=False)
    pd.DataFrame(
        rows_f, columns=["sample_id", "chrom", "start", "end", "methylated_count", "total_count"]
    ).to_csv(tmp_path / "fme.tsv", sep="\t", index=False)
    pd.DataFrame(rows_a, columns=["chrom", "start", "end", "cpg_density"]).to_csv(
        tmp_path / "ann.tsv", sep="\t", index=False
    )

    params = train_calibration(
        matched_wgbs=tmp_path / "wgbs.tsv",
        matched_finaleme=tmp_path / "fme.tsv",
        region_annotation=tmp_path / "ann.tsv",
        n_bins_candidates=[2, 4],
        out_params=tmp_path / "cal.json",
        out_report=tmp_path / "report.json",
    )
    assert params.n_bins in (2, 4)
    assert (tmp_path / "cal.json").exists()
    assert (tmp_path / "report.json").exists()
    # Round-trip
    loaded = CalibrationParams.load(tmp_path / "cal.json")
    np.testing.assert_allclose(loaded.a, params.a)
