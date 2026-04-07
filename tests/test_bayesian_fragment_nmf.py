"""Tests for Bayesian deconvolver, fragment-level EM, and NMF residual discovery (Phase E)."""

from __future__ import annotations

import numpy as np
import pytest

from finaleme_too.config import CoverageTier, MeasurementMode
from finaleme_too.core.deconvolution import BayesianDeconvolver, MLEDeconvolver
from finaleme_too.core.fragment_likelihood import Fragment, FragmentLevelDeconvolver
from finaleme_too.core.observation_model import BetaBinomialModel
from finaleme_too.io.methylation_loader import MarkerObservations
from finaleme_too.postprocessing.residual_analysis import (
    compute_residuals_per_sample,
    discover_residual_components,
)
from finaleme_too.postprocessing.statistical_testing import bayesian_group_comparison


emcee = pytest.importorskip("emcee")


def test_bayesian_deconvolver_recovers_pure_celltype(
    synthetic_observations_pure_celltype, synthetic_reference
):
    """Bayesian posterior mean should match the MLE within ~0.1 on a clean signal."""
    model = BetaBinomialModel().build(
        synthetic_observations_pure_celltype, tier=CoverageTier.HIGH
    )
    mle = MLEDeconvolver().solve(model, synthetic_reference)
    bd = BayesianDeconvolver(n_walkers=32, n_steps=300, burn_in=100, prior_alpha=1.0, seed=0)
    samples = bd.solve(model, synthetic_reference)
    posterior_mean = samples.mean(axis=0)
    assert samples.ndim == 2
    assert samples.shape[1] == synthetic_reference.n_cell_types + 1
    # CellType1 should be dominant in both
    assert posterior_mean[0] > 0.7
    assert abs(posterior_mean[0] - mle[0]) < 0.15


def test_fragment_level_em_recovers_dominant_celltype():
    """Synthesize fragments where every CpG matches CellType1's pattern."""
    rng = np.random.default_rng(0)
    n_markers = 30
    K = 3
    # All non-owner cell types are highly methylated; the owner is unmethylated.
    R = np.full((n_markers, K), 0.95, dtype=np.float64)
    R[:10, 0] = 0.05  # CellType1 owns markers 0-9 (low at those markers)
    R[10:20, 1] = 0.05  # CellType2 owns markers 10-19
    R[20:, 2] = 0.05  # CellType3 owns markers 20-29

    # Generate 100 fragments drawn from CellType1's true profile.
    # P(methylated_i | CellType1) = R[i, 0]
    fragments = []
    for _ in range(100):
        idx = rng.choice(np.arange(n_markers), size=5, replace=False)
        m = (rng.uniform(0, 1, size=5) < R[idx, 0]).astype(np.uint8)
        fragments.append(Fragment(cpg_indices=idx, methylated=m))

    em = FragmentLevelDeconvolver(max_iter=100)
    w = em.solve(fragments, R)
    assert w.shape == (K + 1,)
    assert np.isclose(w.sum(), 1.0)
    # CellType1 should dominate
    assert w[0] > 0.4, f"Expected w[0] > 0.4, got {w[0]:.3f}"


def test_bayesian_group_comparison_high_when_groups_differ():
    """When group A has higher mean for cell type 0, P(A > B) should be near 1."""
    rng = np.random.default_rng(1)
    n_T = 200
    K_total = 4

    posterior_a = {}
    for s in range(4):
        # Group A: cell type 0 around 0.6
        draws = rng.dirichlet([6, 2, 1, 1], size=n_T)
        posterior_a[f"a_{s}"] = draws
    posterior_b = {}
    for s in range(4):
        draws = rng.dirichlet([2, 6, 1, 1], size=n_T)
        posterior_b[f"b_{s}"] = draws

    posterior = {**posterior_a, **posterior_b}
    sample_groups = {**{f"a_{s}": "A" for s in range(4)}, **{f"b_{s}": "B" for s in range(4)}}

    results = bayesian_group_comparison(
        posterior_samples_by_sample=posterior,
        sample_groups=sample_groups,
        cell_type_names=["CT1", "CT2", "CT3"],
        contrasts=[("A", "B")],
    )
    by_ct = {r.cell_type: r for r in results}
    # CT1 should have P(A > B) close to 1 (statistic field)
    assert by_ct["CT1"].statistic > 0.9
    assert by_ct["CT2"].statistic < 0.1


def test_compute_residuals_per_sample():
    R = np.array([[0.1, 0.9], [0.8, 0.2]])
    obs_meth = np.array([0.15, 0.78])
    w = np.array([0.5, 0.4, 0.1])  # K=2 + unknown
    out = compute_residuals_per_sample(w, R, obs_meth)
    assert "unexplained_fraction" in out
    assert abs(out["unexplained_fraction"] - 0.1) < 1e-9
    assert "mean_residual" in out
    assert "residual_sd" in out


def test_nmf_residual_discovery():
    rng = np.random.default_rng(2)
    # Synthesize a residual matrix with one strong rank-1 component
    profile = rng.uniform(0, 1, size=20)
    loadings = rng.uniform(0.3, 1.0, size=10)
    R = np.outer(loadings, profile) + rng.normal(0, 0.01, size=(10, 20))
    out = discover_residual_components(R, n_components=2)
    assert out["components"].shape[1] == 20
    assert out["loadings"].shape[0] == 10
    assert out["explained_variance_ratio"].size == 2
