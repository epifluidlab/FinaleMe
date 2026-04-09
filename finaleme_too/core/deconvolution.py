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

    Enriched output fields (architecture §10.1 / §10.2 / §10.6):
        mean_dispersion   — (K,) per-cell-type mean phi over that cell type's
                            discriminative markers (for the per-sample TSV)
        mean_coverage     — mean total reads per marker actually used
        n_markers_used    — total markers passing the tier coverage filter
        pct_imputed       — fraction of markers filled in by CohortImputer
        binarization_flag — PASS / WARN / FAIL from FinaleMe inference QC
                            (FinaleMe mode only; None otherwise). Renamed
                            from calibration_flag in v3 — under v3 the
                            FinaleMe path uses context-dependent
                            binarization with error rates instead of
                            continuous calibration, but the field name
                            stays semantically the same: it's the
                            sample-level inference QC for the FinaleMe
                            preprocessing step (whichever path was used).
        hemolysis_flag    — from sample metadata (None if missing)
        overall_qc        — PASS / WARN / FAIL derived from qc_flags
        residuals         — (M_used,) per-marker residual (mu_obs - mu_pred)
                            with marker coordinates; populated only when the
                            cohort-level NMF residual analysis is enabled
        marker_chrom      — (M_used,) chromosome per marker actually used
        marker_start      — (M_used,) start position per marker actually used
        marker_end        — (M_used,) end position per marker actually used
    """

    sample_id: str
    cell_types: list[str]
    proportions: np.ndarray  # (K+1,)
    ci_lower: np.ndarray  # (K+1,)
    ci_upper: np.ndarray  # (K+1,)
    # Legacy metric kept only for backward compatibility in in-memory
    # objects. It is no longer used for reliability assignment and is not
    # emitted in the default per-sample/cohort TSV outputs.
    p_goodness: np.ndarray | None  # (K,) — deprecated
    p_detection: np.ndarray  # (K+1,)
    reliability: np.ndarray  # str array (K+1,)
    n_markers: np.ndarray  # (K,)
    # New reliability metrics (K+1; unknown included as last entry):
    #   likelihood_score — weighted per-marker log-likelihood gain vs null
    #   p_likelihood     — raw LRT p-value vs ablated-null (smaller is better)
    #   q_likelihood     — BH-adjusted p-value across known cell types
    likelihood_score: np.ndarray | None = None  # (K+1,)
    p_likelihood: np.ndarray | None = None  # (K+1,)
    q_likelihood: np.ndarray | None = None  # (K+1,)
    # Deprecated legacy metric (kept for backward compatibility only).
    effect_size: np.ndarray | None = None  # (K+1,)
    bootstrap_proportions: np.ndarray | None = None  # (B, K+1)
    posterior_samples: np.ndarray | None = None  # (T, K+1)
    coverage_tier: CoverageTier = CoverageTier.HIGH
    qc_flags: list[str] = field(default_factory=list)
    # ---- enriched output fields ----
    mean_dispersion: np.ndarray | None = None
    mean_coverage: float = 0.0
    n_markers_used: int = 0
    pct_imputed: float = 0.0
    binarization_flag: str | None = None
    hemolysis_flag: bool | None = None
    overall_qc: str = "PASS"
    residuals: np.ndarray | None = None
    marker_chrom: np.ndarray | None = None
    marker_start: np.ndarray | None = None
    marker_end: np.ndarray | None = None
    # Optional per-sample FinaleMe binarization diagnostics. Populated only
    # for FinaleMe runs with binarization enabled; used by
    # io.output_writer.write_binarization_debug for troubleshooting.
    binarization_debug: dict | None = None


class MLEDeconvolver:
    """Constrained MLE deconvolution via SLSQP."""

    def __init__(self, max_iter: int = 200, tol: float = 1e-8):
        self.max_iter = max_iter
        self.tol = tol

    def solve(
        self,
        model,
        reference: "ReferencePanel",
        marker_subset: np.ndarray | None = None,
    ) -> np.ndarray:
        """Return ŵ of length K+1 (with unknown last).

        Dispatches on the observation model type:

        * ``ObservationModel`` (beta-binomial, WGBS or v2 FinaleMe path)
          → ``_solve_betabinom`` with the per-marker counts + dispersion.
        * ``BinarizationObservationModel`` (v3 FinaleMe path) →
          ``_solve_binarization`` with the precomputed per-marker linear
          coefficient matrix.

        The simplex constraint, bounds, SLSQP options, failure-fallback
        scaffolding, and result post-processing are mode-independent and
        live in ``_run_slsqp`` so both branches share them.
        """
        # Local import to avoid a top-level circular import (the
        # binarization module imports MarkerObservations from io which
        # depends on config, while deconvolution.py is imported by
        # observation_model.py which is imported by binarization).
        from finaleme_too.core.observation_model_binarization import (
            BinarizationObservationModel,
        )

        if isinstance(model, BinarizationObservationModel):
            return self._solve_binarization(model, marker_subset)
        return self._solve_betabinom(model, reference, marker_subset)

    def _solve_betabinom(
        self,
        model: ObservationModel,
        reference: "ReferencePanel",
        marker_subset: np.ndarray | None = None,
    ) -> np.ndarray:
        """Beta-binomial MLE solver (WGBS mode and v2 FinaleMe fallback)."""
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

        return self._run_slsqp(
            neg_ll, neg_grad, K_total, sample_id=model.sample_id
        )

    def _solve_binarization(
        self,
        model,  # BinarizationObservationModel
        marker_subset: np.ndarray | None = None,
    ) -> np.ndarray:
        """v3 FinaleMe MLE solver: binomial likelihood with error rates.

        The per-marker likelihood is

            P(observed_state_i | w) = coef[i] @ w

        which is linear in ``w``. The negative log-likelihood and gradient
        are therefore much cheaper than the beta-binomial path (no digamma):

            neg_ll(w)   = -Σ_i ω_i · log(coef[i] @ w)
            neg_grad(w) = -(ω_i / (coef[i] @ w)) @ coef
        """
        coef = model.coef          # (M_valid, K+1) — already augmented with unknown
        weights = model.weights    # (M_valid,)

        if marker_subset is not None:
            coef = coef[marker_subset]
            weights = weights[marker_subset]

        # Filter rows with NaN coefficients defensively (shouldn't happen
        # under normal binarization but matches the betabinom guard).
        valid = np.all(np.isfinite(coef), axis=1)
        if int(np.sum(valid)) < coef.shape[1]:
            uniform = np.zeros(coef.shape[1], dtype=np.float64)
            uniform[-1] = 1.0
            return uniform

        coef_v = coef[valid]
        w_obj = weights[valid]
        K_total = coef_v.shape[1]

        def neg_ll(w: np.ndarray) -> float:
            p_obs = coef_v @ w
            p_obs = np.clip(p_obs, 1e-15, 1.0)
            return -float(np.sum(w_obj * np.log(p_obs)))

        def neg_grad(w: np.ndarray) -> np.ndarray:
            p_obs = coef_v @ w
            p_obs = np.clip(p_obs, 1e-15, 1.0)
            # ∂(-log p_i)/∂w_j = -coef[i,j] / p_i
            # weighted: -ω_i · coef[i,j] / p_i
            # sum over i: -((ω/p) @ coef)
            return -((w_obj / p_obs) @ coef_v)

        return self._run_slsqp(
            neg_ll, neg_grad, K_total, sample_id=model.sample_id
        )

    def _run_slsqp(
        self,
        neg_ll,
        neg_grad,
        K_total: int,
        sample_id: str,
    ) -> np.ndarray:
        """Mode-independent SLSQP runner used by both solver paths.

        Encapsulates the simplex constraint, bounds, initial guess,
        scipy.optimize.minimize call with failure handling, and the
        post-processing (clip to [0, 1], renormalize). Both
        ``_solve_betabinom`` and ``_solve_binarization`` reuse this so
        the simplex / convergence behavior is identical.
        """
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
                f"SLSQP failed for sample {sample_id}: {exc}"
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
    """Bayesian Dirichlet-BetaBinomial deconvolver via emcee (math doc §4).

    Posterior: P(w | data) ∝ Dir(w | α0) * Π_i BetaBinom(k_i | n_i, μ_i(w), φ_i)

    Sampling strategy: parameterize w as a stick-breaking ALR (logit) vector
    of K free dimensions in unconstrained space, transform to a (K+1)-simplex
    via softmax inside the log-posterior. This avoids manual constraint
    enforcement and works well with emcee's affine-invariant sampler.
    """

    def __init__(
        self,
        n_walkers: int = 64,
        n_steps: int = 2000,
        burn_in: int = 500,
        prior_alpha: float = 1.0,
        seed: int | None = None,
        binarization_count_weight: float = 1.0,
    ):
        self.n_walkers = n_walkers
        self.n_steps = n_steps
        self.burn_in = burn_in
        self.prior_alpha = prior_alpha
        self.seed = seed
        self.binarization_count_weight = max(float(binarization_count_weight), 0.0)

    def solve(
        self,
        model: ObservationModel,
        reference: "ReferencePanel",
        marker_subset: np.ndarray | None = None,
    ) -> np.ndarray:
        """Return posterior samples shape (n_keep, K+1)."""
        try:
            import emcee
        except ImportError as exc:  # pragma: no cover
            raise DeconvolutionFailedError(
                "emcee not installed; install with `pip install finaleme-too[bayesian]`"
            ) from exc

        # Local import to avoid top-level circular imports.
        from finaleme_too.core.observation_model_binarization import (
            BinarizationObservationModel,
        )

        if isinstance(model, BinarizationObservationModel):
            return self._solve_hybrid_binarization(
                model=model,
                reference=reference,
                marker_subset=marker_subset,
                emcee_module=emcee,
            )

        R_full = MLEDeconvolver._augmented_reference(reference)
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

        valid = (n_arr > 0) & np.all(np.isfinite(R_full), axis=1)
        R = R_full[valid]
        k = k_arr[valid].astype(np.float64)
        n = n_arr[valid].astype(np.float64)
        phi_v = phi[valid]
        w_obj = weights[valid]
        K_total = R.shape[1]
        K_free = K_total - 1  # softmax-parameterized

        log_prior_const = -np.sum(np.log(np.maximum(self.prior_alpha, 1e-9))) * 0  # constant

        def softmax(z: np.ndarray) -> np.ndarray:
            z_max = np.max(z)
            ez = np.exp(z - z_max)
            return ez / ez.sum()

        def log_posterior(z: np.ndarray) -> float:
            # Last component is implicit zero in z; soft-max over [z; 0]
            full_z = np.concatenate([z, np.array([0.0])])
            w = softmax(full_z)
            mu = R @ w
            mu = np.clip(mu, 1e-9, 1.0 - 1e-9)
            ll = log_likelihood_per_marker(k, n, mu, phi_v)
            log_lik = float(np.sum(w_obj * ll))
            # Symmetric Dirichlet log prior on w (constant terms dropped).
            # Guard against alpha=1 (uniform Dirichlet) where the term is 0
            # but the explicit 0 * log(0) = NaN.
            if abs(self.prior_alpha - 1.0) < 1e-12:
                log_prior = 0.0
            else:
                log_prior = float(
                    (self.prior_alpha - 1.0) * np.sum(np.log(np.maximum(w, 1e-300)))
                )
            # Jacobian of softmax (log |det J|) — improper, but cancels for paired comparison
            return log_lik + log_prior

        rng = np.random.default_rng(self.seed)
        # Initial walkers: small noise around uniform (zero in unconstrained space)
        p0 = rng.normal(0, 0.1, size=(self.n_walkers, K_free))

        sampler = emcee.EnsembleSampler(self.n_walkers, K_free, log_posterior)
        sampler.run_mcmc(p0, self.n_steps, progress=False)

        chain = sampler.get_chain(discard=self.burn_in, flat=True)  # (n_keep, K_free)
        # Convert to (K+1)-simplex
        n_keep = chain.shape[0]
        out = np.empty((n_keep, K_total), dtype=np.float64)
        for i in range(n_keep):
            full_z = np.concatenate([chain[i], np.array([0.0])])
            out[i] = softmax(full_z)
        return out

    def _solve_hybrid_binarization(
        self,
        model,  # BinarizationObservationModel
        reference: "ReferencePanel",
        marker_subset: np.ndarray | None,
        emcee_module,
    ) -> np.ndarray:
        """Hybrid Bayesian inference for v3 binarization mode.

        Posterior combines:
          1) state likelihood: log P(called_state_i | w) from coef @ w
          2) count likelihood: beta-binomial log L(k_i, n_i | mu_i(w), phi_i)

        The count term is normalized per read and weighted by
        ``self.binarization_count_weight`` to keep scales compatible.
        """
        coef = np.asarray(model.coef, dtype=np.float64)
        weights = np.asarray(model.weights, dtype=np.float64)
        if model.k is None or model.n is None or model.dispersion is None:
            raise DeconvolutionFailedError(
                "BinarizationObservationModel is missing k/n/dispersion for "
                "hybrid Bayesian inference."
            )
        k = np.asarray(model.k, dtype=np.float64)
        n = np.asarray(model.n, dtype=np.float64)
        phi = np.asarray(model.dispersion, dtype=np.float64)

        # Continuous reference (not binarized) for the count-likelihood term.
        R_full = MLEDeconvolver._augmented_reference(reference)
        valid_mask = np.asarray(model.valid_mask, dtype=bool)
        R = R_full[valid_mask]

        if marker_subset is not None:
            coef = coef[marker_subset]
            weights = weights[marker_subset]
            k = k[marker_subset]
            n = n[marker_subset]
            phi = phi[marker_subset]
            R = R[marker_subset]

        valid = (
            np.all(np.isfinite(coef), axis=1)
            & np.all(np.isfinite(R), axis=1)
            & np.isfinite(k)
            & np.isfinite(n)
            & np.isfinite(phi)
            & (n > 0)
        )
        if int(np.sum(valid)) < R.shape[1]:
            raise DeconvolutionFailedError(
                "Not enough valid markers for hybrid Bayesian binarization solve."
            )

        coef_v = coef[valid]
        weights_v = weights[valid]
        k_v = k[valid]
        n_v = n[valid]
        phi_v = phi[valid]
        R_v = R[valid]

        K_total = coef_v.shape[1]
        K_free = K_total - 1

        def softmax(z: np.ndarray) -> np.ndarray:
            z_max = np.max(z)
            ez = np.exp(z - z_max)
            return ez / ez.sum()

        def log_posterior(z: np.ndarray) -> float:
            full_z = np.concatenate([z, np.array([0.0])])
            w = softmax(full_z)

            # Binarization state likelihood
            p_obs = np.clip(coef_v @ w, 1e-15, 1.0)
            ll_state = float(np.sum(weights_v * np.log(p_obs)))

            # Count likelihood (beta-binomial), normalized per read to
            # avoid overwhelming the state term at deep coverage.
            mu = np.clip(R_v @ w, 1e-9, 1.0 - 1e-9)
            ll_count_vec = log_likelihood_per_marker(k_v, n_v, mu, phi_v)
            ll_count_per_read = ll_count_vec / np.maximum(n_v, 1.0)
            ll_count = float(np.sum(weights_v * ll_count_per_read))

            log_lik = ll_state + self.binarization_count_weight * ll_count

            if abs(self.prior_alpha - 1.0) < 1e-12:
                log_prior = 0.0
            else:
                log_prior = float(
                    (self.prior_alpha - 1.0) * np.sum(np.log(np.maximum(w, 1e-300)))
                )
            return log_lik + log_prior

        rng = np.random.default_rng(self.seed)
        p0 = rng.normal(0, 0.1, size=(self.n_walkers, K_free))

        sampler = emcee_module.EnsembleSampler(self.n_walkers, K_free, log_posterior)
        sampler.run_mcmc(p0, self.n_steps, progress=False)

        chain = sampler.get_chain(discard=self.burn_in, flat=True)
        n_keep = chain.shape[0]
        out = np.empty((n_keep, K_total), dtype=np.float64)
        for i in range(n_keep):
            full_z = np.concatenate([chain[i], np.array([0.0])])
            out[i] = softmax(full_z)
        return out


__all__ = ["BayesianDeconvolver", "DeconvolutionResult", "MLEDeconvolver", "UNKNOWN_PROFILE"]
