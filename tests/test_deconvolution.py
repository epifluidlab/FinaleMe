"""Tests for the MLE deconvolver and reliability p-values."""

from __future__ import annotations

import numpy as np

from finaleme_too.config import CoverageTier, MeasurementMode
from finaleme_too.core.deconvolution import MLEDeconvolver
from finaleme_too.core.observation_model import BetaBinomialModel
from finaleme_too.core.reliability import (
    assign_reliability,
    compute_fit_metrics,
    compute_p_detection,
    compute_p_goodness,
)
from finaleme_too.core.uncertainty import BootstrapCI
from finaleme_too.io.methylation_loader import MarkerObservations


def test_self_deconvolution_recovers_celltype1(
    synthetic_observations_pure_celltype, synthetic_reference
):
    """Pure CellType1 observations should yield w_CellType1 close to 1.0."""
    model = BetaBinomialModel().build(
        synthetic_observations_pure_celltype, tier=CoverageTier.HIGH
    )
    w_hat = MLEDeconvolver().solve(model, synthetic_reference)
    assert w_hat.shape == (synthetic_reference.n_cell_types + 1,)
    assert np.isclose(np.sum(w_hat), 1.0, atol=1e-6)
    assert np.all(w_hat >= -1e-9)
    # CellType1 should dominate (allow some leakage)
    assert w_hat[0] > 0.85, f"Expected w_CellType1 > 0.85, got {w_hat[0]:.3f}"


def test_unknown_component_dominates_for_uniform_methylation(
    synthetic_marker_regions, synthetic_reference
):
    """All-0.5 methylation should be explained by the unknown component."""
    n = synthetic_marker_regions.n_markers
    obs = MarkerObservations(
        sample_id="all_unknown",
        chrom=synthetic_marker_regions.chrom,
        start=synthetic_marker_regions.start,
        end=synthetic_marker_regions.end,
        k=np.full(n, 25, dtype=np.int32),  # 50%
        n=np.full(n, 50, dtype=np.int32),
        predicted_beta=None,
        mode=MeasurementMode.WGBS,
    )
    model = BetaBinomialModel().build(obs, tier=CoverageTier.HIGH)
    w_hat = MLEDeconvolver().solve(model, synthetic_reference)
    # Unknown is the last component
    assert w_hat[-1] > 0.4, f"Expected unknown > 0.4, got {w_hat[-1]:.3f}"


def test_two_celltype_mixture_recovery(synthetic_marker_regions, synthetic_reference):
    """50/50 mixture of CellType1 and CellType3 should recover both proportions."""
    n_markers = synthetic_marker_regions.n_markers
    rng = np.random.default_rng(0)
    p = 0.5 * synthetic_reference.methylation[:, 0] + 0.5 * synthetic_reference.methylation[:, 2]
    n = np.full(n_markers, 60, dtype=np.int32)
    k = rng.binomial(n.astype(int), p).astype(np.int32)
    obs = MarkerObservations(
        sample_id="mix",
        chrom=synthetic_marker_regions.chrom,
        start=synthetic_marker_regions.start,
        end=synthetic_marker_regions.end,
        k=k,
        n=n,
        predicted_beta=None,
        mode=MeasurementMode.WGBS,
    )
    model = BetaBinomialModel().build(obs, tier=CoverageTier.HIGH)
    w_hat = MLEDeconvolver().solve(model, synthetic_reference)
    # CellType1 (idx 0) and CellType3 (idx 2) should each be near 0.5
    assert abs(w_hat[0] - 0.5) < 0.15, f"CellType1: {w_hat[0]:.3f}"
    assert abs(w_hat[2] - 0.5) < 0.15, f"CellType3: {w_hat[2]:.3f}"
    # Other cell types and unknown should be small
    other = w_hat[1] + w_hat[3] + w_hat[4] + w_hat[-1]
    assert other < 0.20, f"Other components should be small, got {other:.3f}"


def test_constraints_satisfied(synthetic_observations_pure_celltype, synthetic_reference):
    model = BetaBinomialModel().build(
        synthetic_observations_pure_celltype, tier=CoverageTier.HIGH
    )
    w_hat = MLEDeconvolver().solve(model, synthetic_reference)
    assert np.all(w_hat >= -1e-9)
    assert np.isclose(np.sum(w_hat), 1.0, atol=1e-6)


def test_bootstrap_ci_contains_point_estimate(
    synthetic_observations_pure_celltype, synthetic_reference
):
    model = BetaBinomialModel().build(
        synthetic_observations_pure_celltype, tier=CoverageTier.HIGH
    )
    boot = BootstrapCI(n_iterations=50, ci_level=0.95, seed=0).estimate(
        model, synthetic_reference, MLEDeconvolver()
    )
    # CellType1 should be reliably high
    assert boot.point_estimate[0] > 0.8
    # CIs should be ordered
    assert np.all(boot.ci_lower <= boot.ci_upper + 1e-9)


def test_reliability_assignment_table():
    # New reliability uses p_detection + likelihood_score + p_likelihood.
    assert assign_reliability(0.99, 0.05, 1e-6) == "HIGH"
    assert assign_reliability(0.70, 0.01, 1e-2) == "MODERATE"
    assert assign_reliability(0.70, 0.0, 0.5) == "LOW"
    assert assign_reliability(0.10, -0.01, 0.9) == "UNRELIABLE"
    # Unknown-component style call (no fit metrics): detection-only fallback.
    assert assign_reliability(0.99, None, None) == "HIGH"
    assert assign_reliability(0.70, None, None) == "MODERATE"
    assert assign_reliability(0.05, None, None) == "UNRELIABLE"


def test_p_detection_above_noise_floor():
    boot = np.array([0.0, 0.0001, 0.5, 0.4, 0.3])
    p = compute_p_detection(boot, noise_floor=0.001)
    # 3 out of 5 are >= 0.001
    assert abs(p - 0.6) < 1e-9


def test_fit_metrics_positive_for_correct_model(
    synthetic_observations_pure_celltype, synthetic_reference
):
    model = BetaBinomialModel().build(
        synthetic_observations_pure_celltype, tier=CoverageTier.HIGH
    )
    w_hat = MLEDeconvolver().solve(model, synthetic_reference)
    lik, plik = compute_fit_metrics(
        w_hat=w_hat,
        reference_methylation=synthetic_reference.methylation,
        observation=model,
        cell_type_index=0,
        top_n=50,
        binarizer=None,
    )
    assert np.isfinite(lik)
    assert np.isfinite(plik)
    assert lik > 0.0
    assert 0.0 <= plik <= 1.0
