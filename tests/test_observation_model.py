"""Tests for the beta-binomial observation model and primitives."""

from __future__ import annotations

import numpy as np
from scipy.stats import betabinom

from finaleme_too.config import CoverageTier, MeasurementMode
from finaleme_too.core.observation_model import BetaBinomialModel
from finaleme_too.io.methylation_loader import MarkerObservations
from finaleme_too.utils.beta_binomial import (
    gradient_w,
    log_likelihood_per_marker,
    log_pmf,
    variance,
)


def test_log_pmf_matches_scipy():
    n = np.array([20, 50, 100], dtype=np.float64)
    k = np.array([5, 25, 80], dtype=np.float64)
    alpha = np.array([5.0, 10.0, 20.0])
    beta = np.array([5.0, 10.0, 5.0])
    expected = betabinom.logpmf(k.astype(int), n.astype(int), alpha, beta)
    actual = log_pmf(k, n, alpha, beta)
    np.testing.assert_allclose(actual, expected, rtol=1e-10)


def test_log_likelihood_per_marker_matches_scipy():
    n = np.array([20, 50, 100], dtype=np.float64)
    k = np.array([5, 25, 80], dtype=np.float64)
    mu = np.array([0.25, 0.5, 0.8])
    phi = np.array([10.0, 20.0, 50.0])
    expected = betabinom.logpmf(k.astype(int), n.astype(int), mu * phi, (1 - mu) * phi)
    actual = log_likelihood_per_marker(k, n, mu, phi)
    np.testing.assert_allclose(actual, expected, rtol=1e-10)


def test_variance_matches_scipy():
    n = np.array([20, 50, 100], dtype=np.float64)
    mu = np.array([0.25, 0.5, 0.8])
    phi = np.array([10.0, 20.0, 50.0])
    expected = betabinom.var(n.astype(int), mu * phi, (1 - mu) * phi)
    actual = variance(n, mu, phi)
    np.testing.assert_allclose(actual, expected, rtol=1e-10)


def test_gradient_matches_numerical_gradient():
    """Compare analytical gradient to a finite-difference numerical gradient."""
    rng = np.random.default_rng(0)
    M, K = 30, 4
    R = rng.uniform(0, 1, size=(M, K + 1))  # +1 for unknown component
    R[:, -1] = 0.5
    n = rng.integers(20, 60, size=M).astype(np.float64)
    true_w = rng.dirichlet(np.ones(K + 1))
    true_mu = R @ true_w
    k = rng.binomial(n.astype(int), true_mu).astype(np.float64)
    phi = np.full(M, 30.0)

    w0 = rng.dirichlet(np.ones(K + 1))
    mu0 = np.clip(R @ w0, 1e-9, 1 - 1e-9)

    analytic = gradient_w(k, n, mu0, phi, R, weights=None)

    eps = 1e-5
    numeric = np.zeros_like(w0)
    for j in range(K + 1):
        wp = w0.copy()
        wm = w0.copy()
        wp[j] += eps
        wm[j] -= eps
        ll_p = float(np.sum(log_likelihood_per_marker(k, n, np.clip(R @ wp, 1e-9, 1 - 1e-9), phi)))
        ll_m = float(np.sum(log_likelihood_per_marker(k, n, np.clip(R @ wm, 1e-9, 1 - 1e-9), phi)))
        numeric[j] = (ll_p - ll_m) / (2 * eps)

    np.testing.assert_allclose(analytic, numeric, atol=1e-3, rtol=1e-3)


def test_observation_model_builds_for_wgbs_high(synthetic_marker_regions):
    n = 50
    obs = MarkerObservations(
        sample_id="x",
        chrom=synthetic_marker_regions.chrom,
        start=synthetic_marker_regions.start,
        end=synthetic_marker_regions.end,
        k=np.full(n, 25, dtype=np.int32),
        n=np.full(n, 50, dtype=np.int32),
        predicted_beta=None,
        mode=MeasurementMode.WGBS,
    )
    model = BetaBinomialModel().build(obs, tier=CoverageTier.HIGH)
    assert model.n_markers == n
    assert np.all(model.dispersion > 0)
    # WGBS HIGH MLE should produce a finite (positive) shared dispersion
    assert np.all(np.isfinite(model.dispersion))
    # Weights should equal min(n_i,cap)/cap = 50/50 = 1.0 (cap default = 50)
    np.testing.assert_allclose(model.weights, np.ones(n), rtol=1e-10)
