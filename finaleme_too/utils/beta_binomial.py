"""Vectorized beta-binomial distribution helpers.

Implements the formulas from TOO_MATH_FORMULATION_v2.md §2.2 (log-PMF) and
§2.3 (variance). All functions accept and return NumPy arrays and broadcast
naturally over the marker dimension.
"""

from __future__ import annotations

import numpy as np
from scipy import special


def log_pmf(
    k: np.ndarray,
    n: np.ndarray,
    alpha: np.ndarray,
    beta: np.ndarray,
) -> np.ndarray:
    """Log probability mass of the beta-binomial distribution.

    P(k | n, alpha, beta) = C(n, k) * B(k+alpha, n-k+beta) / B(alpha, beta)

    The log form (math doc §2.2) avoids loss of precision:

        log P = log Γ(n+1) - log Γ(k+1) - log Γ(n-k+1)
              + log Γ(k+alpha) + log Γ(n-k+beta) - log Γ(n+alpha+beta)
              - log Γ(alpha) - log Γ(beta) + log Γ(alpha+beta)
    """
    k_arr = np.asarray(k, dtype=np.float64)
    n_arr = np.asarray(n, dtype=np.float64)
    a_arr = np.asarray(alpha, dtype=np.float64)
    b_arr = np.asarray(beta, dtype=np.float64)

    log_binom = (
        special.gammaln(n_arr + 1.0)
        - special.gammaln(k_arr + 1.0)
        - special.gammaln(n_arr - k_arr + 1.0)
    )
    log_beta_num = (
        special.gammaln(k_arr + a_arr)
        + special.gammaln(n_arr - k_arr + b_arr)
        - special.gammaln(n_arr + a_arr + b_arr)
    )
    log_beta_den = (
        special.gammaln(a_arr) + special.gammaln(b_arr) - special.gammaln(a_arr + b_arr)
    )
    return log_binom + log_beta_num - log_beta_den


def log_likelihood_per_marker(
    k: np.ndarray,
    n: np.ndarray,
    mu: np.ndarray,
    phi: np.ndarray,
) -> np.ndarray:
    """Per-marker log-likelihood ℓ_i(w) given expected mean mu_i and dispersion phi_i.

    Parameterizes alpha = mu*phi, beta = (1-mu)*phi (math doc §2.2).
    """
    mu_clip = np.clip(mu, 1e-9, 1.0 - 1e-9)
    alpha = mu_clip * phi
    beta = (1.0 - mu_clip) * phi
    return log_pmf(k, n, alpha, beta)


def variance(
    n: np.ndarray, mu: np.ndarray, phi: np.ndarray
) -> np.ndarray:
    """Beta-binomial variance per math doc §2.3.

        Var(k_i) = n_i * mu_i * (1 - mu_i) * (n_i + phi_i) / (1 + phi_i)
    """
    n_arr = np.asarray(n, dtype=np.float64)
    mu_arr = np.asarray(mu, dtype=np.float64)
    phi_arr = np.asarray(phi, dtype=np.float64)
    return n_arr * mu_arr * (1.0 - mu_arr) * (n_arr + phi_arr) / (1.0 + phi_arr)


def gradient_w(
    k: np.ndarray,
    n: np.ndarray,
    mu: np.ndarray,
    phi: np.ndarray,
    R: np.ndarray,
    weights: np.ndarray | None = None,
) -> np.ndarray:
    """Gradient ∂L/∂w_j of the weighted beta-binomial log-likelihood.

    Math doc §3.3:

        ∂L/∂w_j = Σ_i ω_i * φ_i * r_{ij}
                  * [ψ(k_i + α_i) - ψ(n_i - k_i + β_i) - ψ(α_i) + ψ(β_i)]

    Parameters
    ----------
    k, n : shape (M,)
        Methylated count and total count per marker.
    mu : shape (M,)
        Expected mean methylation under current weights.
    phi : shape (M,)
        Per-marker dispersion.
    R : shape (M, K_total)
        Reference matrix including the unknown component column.
    weights : shape (M,) or None
        Per-marker objective weights ω_i. If None, uses ones.
    """
    mu_clip = np.clip(mu, 1e-9, 1.0 - 1e-9)
    alpha = mu_clip * phi
    beta = (1.0 - mu_clip) * phi
    psi_diff = (
        special.digamma(k + alpha)
        - special.digamma(n - k + beta)
        - special.digamma(alpha)
        + special.digamma(beta)
    )
    if weights is None:
        weights = np.ones_like(k, dtype=np.float64)
    factor = weights * phi * psi_diff  # shape (M,)
    return factor @ R  # shape (K_total,)


def estimate_dispersion_mle(
    k: np.ndarray,
    n: np.ndarray,
    mu: np.ndarray,
    phi_init: float = 50.0,
    bounds: tuple[float, float] = (1.0, 1000.0),
) -> float:
    """One-shot MLE of a single shared dispersion phi over all markers.

    Returns the scalar dispersion that maximizes
    ``Σ_i log P(k_i | n_i, mu_i, phi)`` under the beta-binomial model.

    Uses a bounded Brent search in log-phi space rather than the
    triple-bracket Brent method — the bracketed method fails on
    monotone objectives (which happens whenever the data is
    essentially binomial and the likelihood keeps climbing toward
    the upper bound), whereas ``method="bounded"`` handles the
    corner gracefully.
    """
    from scipy.optimize import minimize_scalar

    k_arr = np.asarray(k, dtype=np.float64).ravel()
    n_arr = np.asarray(n, dtype=np.float64).ravel()
    mu_arr = np.asarray(mu, dtype=np.float64).ravel()
    # Drop rows where coverage is zero — they contribute nothing to the
    # beta-binomial likelihood and can produce NaN in the digamma terms.
    valid = n_arr > 0
    if int(np.sum(valid)) == 0:
        return float(phi_init)
    k_v = k_arr[valid]
    n_v = n_arr[valid]
    mu_v = mu_arr[valid]

    def neg_ll(log_phi: float) -> float:
        phi = float(np.exp(log_phi))
        phi_arr = np.full_like(k_v, phi, dtype=np.float64)
        return -float(
            np.sum(log_likelihood_per_marker(k_v, n_v, mu_v, phi_arr))
        )

    log_lo = float(np.log(bounds[0]))
    log_hi = float(np.log(bounds[1]))
    res = minimize_scalar(
        neg_ll,
        bounds=(log_lo, log_hi),
        method="bounded",
        options={"xatol": 1e-4},
    )
    return float(np.exp(res.x))


__all__ = [
    "log_pmf",
    "log_likelihood_per_marker",
    "variance",
    "gradient_w",
    "estimate_dispersion_mle",
]
