"""Constrained MLE deconvolver and (lazy) Bayesian deconvolver.

Math doc §3.1: maximize Σ_i ω_i ℓ_i(w) subject to w_j ≥ 0 and Σ w_j = 1.
The unknown component is always included: an extra column of 0.5 is appended
to the reference matrix and an extra weight w_0 is fitted.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from scipy.optimize import minimize

from finaleme_too.config import CoverageTier
from finaleme_too.core.observation_model import ObservationModel
from finaleme_too.exceptions import DeconvolutionFailedError
from finaleme_too.utils.beta_binomial import gradient_w, log_likelihood_per_marker

if TYPE_CHECKING:
    from finaleme_too.io.reference_panel import ReferencePanel

UNKNOWN_PROFILE = 0.5  # flat methylation for the unknown component


@dataclass
class DeconvolutionResult:
    """Per-sample deconvolution result.

    proportions has length K+1 where index K is the unknown component.
    cell_types is the K cell-type names (without unknown). The unknown
    component is named "Unknown" in output writers.
    """

    sample_id: str
    cell_types: list[str]
    proportions: np.ndarray  # (K+1,)
    ci_lower: np.ndarray  # (K+1,)
    ci_upper: np.ndarray  # (K+1,)
    p_goodness: np.ndarray  # (K,) — no p_goodness for unknown
    p_detection: np.ndarray  # (K+1,)
    reliability: np.ndarray  # str array (K+1,)
    n_markers: np.ndarray  # (K,)
    bootstrap_proportions: np.ndarray | None = None  # (B, K+1)
    posterior_samples: np.ndarray | None = None  # (T, K+1)
    coverage_tier: CoverageTier = CoverageTier.HIGH
    qc_flags: list[str] = field(default_factory=list)


class MLEDeconvolver:
    """Constrained MLE deconvolution via SLSQP."""

    def __init__(self, max_iter: int = 200, tol: float = 1e-8):
        self.max_iter = max_iter
        self.tol = tol

    def solve(
        self,
        model: ObservationModel,
        reference: "ReferencePanel",
        marker_subset: np.ndarray | None = None,
    ) -> np.ndarray:
        """Return ŵ of length K+1 (with unknown last)."""
        R_full = self._augmented_reference(reference)
        k_arr = model.k
        n_arr = model.n
        phi = model.dispersion
        weights = model.weights

        if marker_subset is not None:
            R_full = R_full[marker_subset]
            k_arr = k_arr[marker_subset]
            n_arr = n_arr[marker_subset]
            phi = phi[marker_subset]
            weights = weights[marker_subset]

        # Drop markers with zero coverage or NaN reference rows
        valid = (
            (n_arr > 0)
            & np.all(np.isfinite(R_full), axis=1)
        )
        if int(np.sum(valid)) < R_full.shape[1]:
            # Insufficient markers — fall back to uniform with all weight on unknown
            uniform = np.zeros(R_full.shape[1], dtype=np.float64)
            uniform[-1] = 1.0
            return uniform

        R = R_full[valid]
        k = k_arr[valid].astype(np.float64)
        n = n_arr[valid].astype(np.float64)
        phi_v = phi[valid]
        w_obj = weights[valid]
        K_total = R.shape[1]

        # Objective: negative weighted log-likelihood
        def neg_ll(w: np.ndarray) -> float:
            mu = R @ w
            mu = np.clip(mu, 1e-9, 1.0 - 1e-9)
            ll = log_likelihood_per_marker(k, n, mu, phi_v)
            return -float(np.sum(w_obj * ll))

        def neg_grad(w: np.ndarray) -> np.ndarray:
            mu = R @ w
            mu = np.clip(mu, 1e-9, 1.0 - 1e-9)
            return -gradient_w(k, n, mu, phi_v, R, weights=w_obj)

        # Initial guess: uniform over all components
        w0 = np.full(K_total, 1.0 / K_total, dtype=np.float64)
        constraints = [{"type": "eq", "fun": lambda w: np.sum(w) - 1.0,
                        "jac": lambda w: np.ones_like(w)}]
        bounds = [(0.0, 1.0)] * K_total

        try:
            res = minimize(
                neg_ll,
                w0,
                jac=neg_grad,
                method="SLSQP",
                bounds=bounds,
                constraints=constraints,
                options={"maxiter": self.max_iter, "ftol": self.tol},
            )
        except Exception as exc:  # noqa: BLE001
            raise DeconvolutionFailedError(
                f"SLSQP failed for sample {model.sample_id}: {exc}"
            ) from exc

        if not res.success:
            # SLSQP often returns success=False but a usable result; check feasibility.
            w_hat = np.clip(res.x, 0.0, 1.0)
            s = float(np.sum(w_hat))
            if s <= 0:
                w_hat = np.full(K_total, 1.0 / K_total, dtype=np.float64)
            else:
                w_hat = w_hat / s
            return w_hat

        w_hat = np.clip(res.x, 0.0, 1.0)
        s = float(np.sum(w_hat))
        if s <= 0:
            return np.full(K_total, 1.0 / K_total, dtype=np.float64)
        return w_hat / s

    @staticmethod
    def _augmented_reference(reference: "ReferencePanel") -> np.ndarray:
        """Append a column of 0.5 for the always-on unknown component."""
        ref = np.asarray(reference.methylation, dtype=np.float64)
        unknown_col = np.full((ref.shape[0], 1), UNKNOWN_PROFILE, dtype=np.float64)
        return np.hstack([ref, unknown_col])


class BayesianDeconvolver:
    """Bayesian Dirichlet-BetaBinomial deconvolver via emcee.

    Lazy import: only imports emcee when ``solve`` is called. Skips
    gracefully (raises a clear error) if emcee is not installed.
    """

    def __init__(
        self,
        n_walkers: int = 64,
        n_steps: int = 2000,
        burn_in: int = 500,
        prior_alpha: float = 1.0,
    ):
        self.n_walkers = n_walkers
        self.n_steps = n_steps
        self.burn_in = burn_in
        self.prior_alpha = prior_alpha

    def solve(
        self,
        model: ObservationModel,
        reference: "ReferencePanel",
        marker_subset: np.ndarray | None = None,
    ) -> np.ndarray:
        try:
            import emcee  # noqa: F401  (lazy)
        except ImportError as exc:  # pragma: no cover
            raise DeconvolutionFailedError(
                "emcee not installed; install with `pip install finaleme-too[bayesian]`"
            ) from exc
        # Implementation deferred to Phase E (P3).
        raise NotImplementedError("BayesianDeconvolver.solve is implemented in Phase E")


__all__ = ["BayesianDeconvolver", "DeconvolutionResult", "MLEDeconvolver", "UNKNOWN_PROFILE"]
