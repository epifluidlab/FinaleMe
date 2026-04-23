"""Constrained MLE deconvolver and (lazy) Bayesian deconvolver.

Math doc §3.1: maximize Σ_i ω_i ℓ_i(w) subject to w_j ≥ 0 and Σ w_j = 1.
The unknown component is always included: an extra column of 0.5 is appended
to the reference matrix and an extra weight w_0 is fitted.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from scipy import special
from scipy.optimize import minimize

from finaleme_too.config import CoverageTier
from finaleme_too.core.observation_model import ObservationModel
from finaleme_too.exceptions import DeconvolutionFailedError
from finaleme_too.utils.beta_binomial import gradient_w, log_likelihood_per_marker

if TYPE_CHECKING:
    from finaleme_too.io.reference_panel import ReferencePanel

UNKNOWN_PROFILE = 0.5  # flat methylation for the unknown component
log = logging.getLogger(__name__)


def hierarchical_marker_log_likelihood(
    *,
    w: np.ndarray,
    reference_continuous: np.ndarray,
    k: np.ndarray,
    n: np.ndarray,
    phi: np.ndarray,
    tau_low: np.ndarray,
    tau_high: np.ndarray,
    call_zone_prob: np.ndarray,
    call_weight: float,
) -> np.ndarray:
    """Per-marker hierarchical log-likelihood for FinaleMe binarization.

    The joint marker term is:
      log P(k | n, mu(w), phi) + log E_posterior[q(call | zone(p))^call_weight]

    where the expectation is under p ~ Beta(k+alpha, n-k+beta), with
    alpha=mu*phi and beta=(1-mu)*phi.
    """
    mu = np.clip(reference_continuous @ w, 1e-9, 1.0 - 1e-9)
    ll_count = log_likelihood_per_marker(k, n, mu, phi)

    alpha = mu * phi
    beta = (1.0 - mu) * phi
    post_a = k + alpha
    post_b = (n - k) + beta

    t_low = np.clip(np.minimum(tau_low, tau_high), 0.0, 1.0)
    t_high = np.clip(np.maximum(tau_low, tau_high), 0.0, 1.0)

    p_low = special.betainc(post_a, post_b, t_low)
    p_high = 1.0 - special.betainc(post_a, post_b, t_high)
    p_mid = np.clip(1.0 - p_low - p_high, 0.0, 1.0)

    cz = np.clip(call_zone_prob, 1e-12, 1.0)
    cw = float(np.clip(call_weight, 0.0, 1.0))
    if abs(cw - 1.0) < 1e-12:
        g_low = cz[:, 0]
        g_mid = cz[:, 1]
        g_high = cz[:, 2]
    elif cw <= 0.0:
        g_low = np.ones_like(p_low)
        g_mid = np.ones_like(p_mid)
        g_high = np.ones_like(p_high)
    else:
        g_low = np.power(cz[:, 0], cw)
        g_mid = np.power(cz[:, 1], cw)
        g_high = np.power(cz[:, 2], cw)

    q_eff = np.clip(g_low * p_low + g_mid * p_mid + g_high * p_high, 1e-15, 1.0)
    ll_call = np.log(q_eff)
    return ll_count + ll_call


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
    """Constrained MLE / NNLS deconvolution.

    Supports two solver modes:

    * ``solver_method="mle"`` (default): SLSQP-based constrained MLE with
      beta-binomial or binarization likelihood. Uses the distributional
      assumptions encoded in the observation model.

    * ``solver_method="nnls"``: Lawson-Hanson non-negative least squares
      with a post-hoc sum=1 renormalization. Matches the behavior of the
      Java ``BetaValueDeconvolution`` tool (L2 residual minimization over
      per-marker predicted methylation). More robust than SLSQP when the
      reference/sample signal has been degraded by hard binarization,
      at the cost of dropping the likelihood-based uncertainty model.

    The ``disable_unknown`` flag is orthogonal to the solver choice: when
    True the 0.5 Unknown column is not appended to the reference, so all
    weight is forced onto the K known cell types (matches Java default
    behavior). The returned ŵ is always padded to length K+1 for
    compatibility (with w[-1]=0 when disabled).
    """

    def __init__(
        self,
        max_iter: int = 200,
        tol: float = 1e-8,
        binarization_count_weight: float = 1.0,
        binarization_model: str = "legacy",
        hierarchical_quadrature_points: int = 24,
        unknown_prior_weight: float = 0.0,
        unknown_prior_weight_auto: bool = False,
        unknown_prior_weight_auto_alpha: float = 0.01,
        solver_method: str = "mle",
        disable_unknown: bool = False,
    ):
        self.max_iter = max_iter
        self.tol = tol
        # Hybrid FinaleMe-v3 objective weight:
        #   logL = logL_state + binarization_count_weight * logL_count
        # Set to 0.0 to recover the legacy state-only solver.
        self.binarization_count_weight = max(float(binarization_count_weight), 0.0)
        self.binarization_model = str(binarization_model or "legacy").lower()
        self.hierarchical_quadrature_points = max(
            int(hierarchical_quadrature_points), 2
        )
        # Point-estimate solver: "mle" (SLSQP) or "nnls" (Lawson-Hanson).
        sm = str(solver_method or "mle").lower()
        if sm not in ("mle", "nnls"):
            raise ValueError(
                f"solver_method must be 'mle' or 'nnls', got {sm!r}"
            )
        self.solver_method = sm
        # When True, skip appending the 0.5 Unknown column to the reference.
        # Affects both MLE and NNLS paths.
        self.disable_unknown = bool(disable_unknown)
        # Prior on the Unknown component w_u: adds -lambda * log(1 - w_u + eps)
        # to the negative log-likelihood, equivalent to a Beta(1, lambda+1)
        # prior marginally on w_u. See design §X.Y and docstring below.
        self.unknown_prior_weight = max(float(unknown_prior_weight), 0.0)
        self.unknown_prior_weight_auto = bool(unknown_prior_weight_auto)
        self.unknown_prior_weight_auto_alpha = max(
            float(unknown_prior_weight_auto_alpha), 0.0
        )
        # Track which sample_ids we have already logged the effective
        # Unknown-prior lambda for, so bootstrap's 100s of solve() calls
        # do not flood the log.
        self._logged_prior_samples: set[str] = set()

    def _effective_unknown_prior_weight(
        self,
        m_valid: int,
        sample_id: str | None = None,
    ) -> float:
        """Resolve the Unknown-prior lambda for a given valid-marker count.

        When ``unknown_prior_weight_auto`` is enabled, returns
        ``alpha * M_valid`` so the prior:likelihood balance stays stable as
        the number of usable markers varies between samples. Otherwise
        returns the raw ``unknown_prior_weight``.

        Logs the resolved lambda the first time this method is called for
        a given ``sample_id`` (subsequent calls, e.g. from bootstrap
        resampling, are silenced to avoid log spam).
        """
        if self.unknown_prior_weight_auto:
            lam = self.unknown_prior_weight_auto_alpha * float(max(m_valid, 0))
        else:
            lam = self.unknown_prior_weight
        if (
            lam > 0.0
            and sample_id is not None
            and sample_id not in self._logged_prior_samples
        ):
            if self.unknown_prior_weight_auto:
                log.info(
                    "[%s] Unknown prior (auto): M_valid=%d, alpha=%.4f -> lambda=%.3f",
                    sample_id,
                    int(m_valid),
                    self.unknown_prior_weight_auto_alpha,
                    lam,
                )
            else:
                log.info(
                    "[%s] Unknown prior (fixed): lambda=%.3f (M_valid=%d)",
                    sample_id,
                    lam,
                    int(m_valid),
                )
            self._logged_prior_samples.add(sample_id)
        return lam

    def solve(
        self,
        model,
        reference: "ReferencePanel",
        marker_subset: np.ndarray | None = None,
    ) -> np.ndarray:
        """Return ŵ of length K+1 (with unknown last; 0 when disabled).

        Dispatches on the observation model type and solver method:

        * ``solver_method="nnls"`` → ``_solve_nnls`` regardless of model type.
          Uses k/n as the sample observation and solves constrained NNLS
          (L2 residual minimization) over the reference methylation matrix.
          Matches Java ``BetaValueDeconvolution`` when ``disable_unknown=True``.
        * ``solver_method="mle"`` (default) with ``ObservationModel`` →
          ``_solve_betabinom`` (beta-binomial MLE via SLSQP).
        * ``solver_method="mle"`` with ``BinarizationObservationModel`` →
          ``_solve_binarization`` or ``_solve_hierarchical_binarization``.

        When ``disable_unknown`` is True the solver works over K columns
        only; the returned vector is padded to K+1 with 0 for Unknown so
        downstream code can always index ``w[-1]``.
        """
        # Local import to avoid a top-level circular import (the
        # binarization module imports MarkerObservations from io which
        # depends on config, while deconvolution.py is imported by
        # observation_model.py which is imported by binarization).
        from finaleme_too.core.observation_model_binarization import (
            BinarizationObservationModel,
        )

        if self.solver_method == "nnls":
            return self._solve_nnls(model, reference, marker_subset)

        if isinstance(model, BinarizationObservationModel):
            obs_model = str(getattr(model, "binarization_model", "legacy")).lower()
            if self.binarization_model == "hierarchical" or obs_model == "hierarchical":
                return self._solve_hierarchical_binarization(model, marker_subset)
            return self._solve_binarization(model, reference, marker_subset)
        return self._solve_betabinom(model, reference, marker_subset)

    def _solve_nnls(
        self,
        model,
        reference: "ReferencePanel",
        marker_subset: np.ndarray | None = None,
    ) -> np.ndarray:
        """Non-negative least squares (Lawson-Hanson).

        Minimizes ``||y - X w||^2`` subject to ``w >= 0``, then
        renormalizes so ``sum(w) = 1``. Matches Java
        ``BetaValueDeconvolution`` when combined with ``--disable-unknown``.

        Two input paths, chosen from the observation model:

        1. ``BinarizationObservationModel`` with ``hard_threshold`` set
           (``--binarizeThreshold`` active): uses the SAME binarized inputs
           as Java — ``y = called_state`` (0 or 1 per marker) and ``X =
           reference_binary`` (0/1 after hard threshold; with 0.5 Unknown
           column when enabled). This is the exact-parity path.

        2. Otherwise (continuous / no hard threshold): uses ``y = k/n``
           (continuous observed methylation) and ``X = raw reference
           methylation``. This is the beta-binomial-MLE analogue.

        Optional Tikhonov regularization on the Unknown weight via
        ``unknown_prior_weight`` (or auto-scaled variant): adds a synthetic
        penalty row ``[0, ..., 0, sqrt(lambda)]`` with target 0 so the
        augmented-NNLS objective becomes ``||y - X w||^2 + lambda * w_u^2``.
        Has no effect when ``disable_unknown`` is True.

        Returns a weight vector of length K+1 (Unknown last). When
        ``disable_unknown`` is True, w[-1] is always 0.0.
        """
        from scipy.optimize import nnls

        from finaleme_too.core.observation_model_binarization import (
            BinarizationObservationModel,
        )

        K_tissues = int(reference.methylation.shape[1])

        # Path selection: binarized inputs (Java-parity) vs continuous.
        use_binarized = (
            isinstance(model, BinarizationObservationModel)
            and getattr(model, "hard_threshold", None) is not None
            and getattr(model, "reference_binary", None) is not None
            and getattr(model, "called_state", None) is not None
        )

        if use_binarized:
            # Java-parity path: y = called_state (0/1), X = reference_binary.
            # called_state and reference_binary are already subsetted to
            # valid markers (M_valid rows); no need for valid_mask gymnastics.
            called_state = np.asarray(model.called_state, dtype=np.float64)
            ref_bin = np.asarray(model.reference_binary, dtype=np.float64)

            # reference_binary has shape (M_valid, K+1) with Unknown as last
            # column. Drop it when disable_unknown is set.
            if self.disable_unknown and ref_bin.shape[1] == K_tissues + 1:
                ref_bin = ref_bin[:, :K_tissues]

            # Java NNLS applies uniform per-marker weights (no coverage-aware
            # weighting). Use uniform weights here to preserve exact parity.
            # Users can run with --solver mle --unknown-prior-weight if they
            # want the coverage-weighted likelihood objective instead.
            X = ref_bin
            y = called_state
            w_obj = np.ones(ref_bin.shape[0], dtype=np.float64)
        else:
            # Continuous / non-binarized path: use k/n and raw reference.
            R_full = self._augmented_reference(
                reference, disable_unknown=self.disable_unknown
            )

            if isinstance(model, BinarizationObservationModel):
                k_arr = (
                    np.asarray(model.k, dtype=np.float64)
                    if model.k is not None
                    else None
                )
                n_arr = (
                    np.asarray(model.n, dtype=np.float64)
                    if model.n is not None
                    else None
                )
                if k_arr is None or n_arr is None:
                    raise DeconvolutionFailedError(
                        f"Sample {model.sample_id}: NNLS (continuous path) "
                        "needs k and n per marker, but the binarization "
                        "observation model is missing them."
                    )
                valid_mask = getattr(model, "valid_mask", None)
                if valid_mask is not None and (
                    np.asarray(valid_mask).shape[0] == R_full.shape[0]
                ):
                    y_full = np.full(R_full.shape[0], np.nan, dtype=np.float64)
                    vm = np.asarray(valid_mask, dtype=bool)
                    y_valid = np.where(
                        n_arr > 0, k_arr / np.maximum(n_arr, 1), np.nan
                    )
                    y_full[vm] = y_valid
                    weights_full = np.zeros(R_full.shape[0], dtype=np.float64)
                    weights_full[vm] = np.asarray(model.weights, dtype=np.float64)
                else:
                    y_full = np.where(
                        n_arr > 0, k_arr / np.maximum(n_arr, 1), np.nan
                    )
                    weights_full = np.asarray(model.weights, dtype=np.float64)
            else:
                k_arr = np.asarray(model.k, dtype=np.float64)
                n_arr = np.asarray(model.n, dtype=np.float64)
                y_full = np.where(
                    n_arr > 0, k_arr / np.maximum(n_arr, 1), np.nan
                )
                weights_full = np.asarray(model.weights, dtype=np.float64)

            X = R_full
            y = y_full
            w_obj = weights_full

        K_total = X.shape[1]

        if marker_subset is not None:
            X = X[marker_subset]
            y = y[marker_subset]
            w_obj = w_obj[marker_subset]

        # Drop rows with NaN in reference or y
        valid = (
            np.isfinite(y)
            & np.all(np.isfinite(X), axis=1)
            & (w_obj > 0)
        )
        if int(np.sum(valid)) < K_total:
            uniform = np.zeros(K_tissues + 1, dtype=np.float64)
            if self.disable_unknown:
                uniform[:K_tissues] = 1.0 / K_tissues
            else:
                uniform[-1] = 1.0
            return uniform

        X_v = X[valid]
        y_v = y[valid]
        w_obj_v = w_obj[valid]

        # Apply marker weights as a row-wise scaling: NNLS does not
        # natively accept weights, but scaling rows by sqrt(w) yields
        # an equivalent weighted least-squares objective.
        sqrt_w = np.sqrt(np.maximum(w_obj_v, 0.0))
        X_w = X_v * sqrt_w[:, None]
        y_w = y_v * sqrt_w

        # Optional Tikhonov penalty on Unknown weight: augment the system
        # with a synthetic row that penalizes w_unknown.
        # NOTE: only meaningful when Unknown column exists (not disabled).
        if (not self.disable_unknown) and K_total > K_tissues:
            lam = self._effective_unknown_prior_weight(
                int(np.sum(valid)), sample_id=model.sample_id
            )
            if lam > 0.0:
                pen_row = np.zeros((1, K_total), dtype=np.float64)
                pen_row[0, -1] = float(np.sqrt(lam))
                X_w = np.vstack([X_w, pen_row])
                y_w = np.concatenate([y_w, np.array([0.0])])

        try:
            w_hat, _ = nnls(X_w, y_w, maxiter=max(self.max_iter, 3 * K_total))
        except Exception as exc:  # noqa: BLE001
            raise DeconvolutionFailedError(
                f"NNLS failed for sample {model.sample_id}: {exc}"
            ) from exc

        # Renormalize to the simplex. NNLS is unconstrained in sum, so
        # renormalize post-hoc (matches Java BetaValueDeconvolution).
        w_hat = np.clip(w_hat, 0.0, None)
        s = float(np.sum(w_hat))
        if s <= 0:
            # All-zero solution — fall back to uniform over active columns.
            w_hat = np.full(K_total, 1.0 / K_total, dtype=np.float64)
        else:
            w_hat = w_hat / s

        return self._pad_unknown(w_hat, K_tissues)

    def _solve_betabinom(
        self,
        model: ObservationModel,
        reference: "ReferencePanel",
        marker_subset: np.ndarray | None = None,
    ) -> np.ndarray:
        """Beta-binomial MLE solver (WGBS mode and v2 FinaleMe fallback)."""
        R_full = self._augmented_reference(
            reference, disable_unknown=self.disable_unknown
        )
        K_tissues = int(reference.methylation.shape[1])
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
            # Insufficient markers — fall back: all weight on Unknown
            # (if available) or uniform over tissues.
            if self.disable_unknown:
                uniform = np.zeros(K_tissues + 1, dtype=np.float64)
                uniform[:K_tissues] = 1.0 / K_tissues
            else:
                uniform = np.zeros(R_full.shape[1], dtype=np.float64)
                uniform[-1] = 1.0
            return uniform

        R = R_full[valid]
        k = k_arr[valid].astype(np.float64)
        n = n_arr[valid].astype(np.float64)
        phi_v = phi[valid]
        w_obj = weights[valid]
        K_total = R.shape[1]

        # Unknown prior is only meaningful when the Unknown column exists.
        lam = (
            self._effective_unknown_prior_weight(
                int(np.sum(valid)), sample_id=model.sample_id
            )
            if not self.disable_unknown
            else 0.0
        )
        eps_prior = 1e-9

        # Objective: negative weighted log-likelihood + optional Unknown prior
        #   prior term:  -lambda * log(1 - w_unknown + eps)
        # Unknown is the last component (index K_total - 1).
        def neg_ll(w: np.ndarray) -> float:
            mu = R @ w
            mu = np.clip(mu, 1e-9, 1.0 - 1e-9)
            ll = log_likelihood_per_marker(k, n, mu, phi_v)
            nll = -float(np.sum(w_obj * ll))
            if lam > 0.0:
                nll += -lam * float(np.log(1.0 - w[-1] + eps_prior))
            return nll

        def neg_grad(w: np.ndarray) -> np.ndarray:
            mu = R @ w
            mu = np.clip(mu, 1e-9, 1.0 - 1e-9)
            g = -gradient_w(k, n, mu, phi_v, R, weights=w_obj)
            if lam > 0.0:
                # d/dw_u [-lambda * log(1 - w_u + eps)] = lambda / (1 - w_u + eps)
                g[-1] += lam / (1.0 - w[-1] + eps_prior)
            return g

        w_hat = self._run_slsqp(
            neg_ll, neg_grad, K_total, sample_id=model.sample_id
        )
        return self._pad_unknown(w_hat, K_tissues)

    def _solve_binarization(
        self,
        model,  # BinarizationObservationModel
        reference: "ReferencePanel",
        marker_subset: np.ndarray | None = None,
    ) -> np.ndarray:
        """v3 FinaleMe hybrid MLE solver.

        Objective combines:

          1) state likelihood from called U/M states:
               P(observed_state_i | w) = coef[i] @ w

          2) count likelihood from (k_i, n_i) under beta-binomial with
             per-marker means mu_i(w) built from the continuous reference.

        The count term is normalized per read so deep-coverage markers do
        not dominate purely due to larger ``n_i``:

            ll_count = Σ_i ω_i * [ log P(k_i | n_i, mu_i(w), phi_i) / n_i ]

        Full objective:

            ll_total = ll_state + binarization_count_weight * ll_count
        """
        coef = np.asarray(model.coef, dtype=np.float64)       # (M_valid, K+1)
        weights = np.asarray(model.weights, dtype=np.float64)  # (M_valid,)
        K_tissues = int(reference.methylation.shape[1])
        # When disable_unknown is set, drop the trailing Unknown column from
        # the pre-computed coef matrix (built by BinarizationModel with K+1
        # columns by default). Renormalize each row so sum(coef_row)=1 still
        # holds over the K remaining columns (state-likelihood invariant).
        if self.disable_unknown and coef.shape[1] == K_tissues + 1:
            coef_wo = coef[:, :K_tissues]
            row_sums = coef_wo.sum(axis=1, keepdims=True)
            row_sums = np.where(row_sums > 1e-12, row_sums, 1.0)
            coef = coef_wo / row_sums
        K_total = coef.shape[1]

        use_count_term = (
            self.binarization_count_weight > 0.0
            and model.k is not None
            and model.n is not None
            and model.dispersion is not None
        )

        if use_count_term:
            k = np.asarray(model.k, dtype=np.float64)
            n = np.asarray(model.n, dtype=np.float64)
            phi = np.asarray(model.dispersion, dtype=np.float64)
            R_full = self._augmented_reference(
                reference, disable_unknown=self.disable_unknown
            )
            valid_mask = getattr(model, "valid_mask", None)
            if valid_mask is None:
                # Defensive fallback for legacy in-memory models that may
                # predate valid_mask; if shapes already align, consume full R.
                if R_full.shape[0] == coef.shape[0]:
                    R = R_full
                else:
                    log.warning(
                        "Sample %s: missing valid_mask for binarization count term; "
                        "falling back to state-only MLE.",
                        model.sample_id,
                    )
                    use_count_term = False
                    R = np.zeros((0, K_total), dtype=np.float64)
            else:
                vm = np.asarray(valid_mask, dtype=bool)
                if vm.shape[0] != R_full.shape[0]:
                    log.warning(
                        "Sample %s: valid_mask length (%d) != reference markers (%d); "
                        "falling back to state-only MLE.",
                        model.sample_id,
                        vm.shape[0],
                        R_full.shape[0],
                    )
                    use_count_term = False
                    R = np.zeros((0, K_total), dtype=np.float64)
                else:
                    R = R_full[vm]

        if marker_subset is not None:
            coef = coef[marker_subset]
            weights = weights[marker_subset]
            if use_count_term:
                k = k[marker_subset]
                n = n[marker_subset]
                phi = phi[marker_subset]
                R = R[marker_subset]

        # Filter invalid rows defensively; requires non-zero coverage for
        # the count term when enabled.
        valid = np.all(np.isfinite(coef), axis=1)
        if use_count_term:
            valid = (
                valid
                & np.all(np.isfinite(R), axis=1)
                & np.isfinite(k)
                & np.isfinite(n)
                & np.isfinite(phi)
                & (n > 0)
            )

        if int(np.sum(valid)) < K_total:
            if self.disable_unknown:
                uniform = np.zeros(K_tissues + 1, dtype=np.float64)
                uniform[:K_tissues] = 1.0 / K_tissues
            else:
                uniform = np.zeros(K_total, dtype=np.float64)
                uniform[-1] = 1.0
            return uniform

        coef_v = coef[valid]
        w_obj = weights[valid]
        if use_count_term:
            k_v = k[valid]
            n_v = n[valid]
            phi_v = phi[valid]
            R_v = R[valid]

        lam = (
            self._effective_unknown_prior_weight(
                int(np.sum(valid)), sample_id=model.sample_id
            )
            if not self.disable_unknown
            else 0.0
        )
        eps_prior = 1e-9

        def neg_ll(w: np.ndarray) -> float:
            p_obs = coef_v @ w
            p_obs = np.clip(p_obs, 1e-15, 1.0)
            ll_state = float(np.sum(w_obj * np.log(p_obs)))
            if use_count_term:
                mu = R_v @ w
                mu = np.clip(mu, 1e-9, 1.0 - 1e-9)
                ll_count_vec = log_likelihood_per_marker(k_v, n_v, mu, phi_v)
                ll_count_per_read = ll_count_vec / np.maximum(n_v, 1.0)
                ll_count = float(np.sum(w_obj * ll_count_per_read))
                nll = -(ll_state + self.binarization_count_weight * ll_count)
            else:
                nll = -ll_state
            if lam > 0.0:
                nll += -lam * float(np.log(1.0 - w[-1] + eps_prior))
            return nll

        def neg_grad(w: np.ndarray) -> np.ndarray:
            p_obs = coef_v @ w
            p_obs = np.clip(p_obs, 1e-15, 1.0)
            grad_state = (w_obj / p_obs) @ coef_v
            if use_count_term:
                mu = R_v @ w
                mu = np.clip(mu, 1e-9, 1.0 - 1e-9)
                grad_count = gradient_w(
                    k_v,
                    n_v,
                    mu,
                    phi_v,
                    R_v,
                    weights=(w_obj / np.maximum(n_v, 1.0)),
                )
                g = -(grad_state + self.binarization_count_weight * grad_count)
            else:
                g = -grad_state
            if lam > 0.0:
                g[-1] += lam / (1.0 - w[-1] + eps_prior)
            return g

        w_hat = self._run_slsqp(
            neg_ll, neg_grad, K_total, sample_id=model.sample_id
        )
        return self._pad_unknown(w_hat, K_tissues)

    def _solve_hierarchical_binarization(
        self,
        model,  # BinarizationObservationModel
        marker_subset: np.ndarray | None = None,
    ) -> np.ndarray:
        """Hierarchical FinaleMe MLE solve with joint count+state likelihood."""
        if (
            model.reference_continuous is None
            or model.k is None
            or model.n is None
            or model.dispersion is None
            or model.tau_low_per is None
            or model.tau_high_per is None
            or model.call_zone_prob is None
        ):
            raise DeconvolutionFailedError(
                f"Sample {model.sample_id}: hierarchical binarization model is missing "
                "reference_continuous/k/n/dispersion/tau/call_zone_prob."
            )

        R = np.asarray(model.reference_continuous, dtype=np.float64)
        k = np.asarray(model.k, dtype=np.float64)
        n = np.asarray(model.n, dtype=np.float64)
        phi = np.asarray(model.dispersion, dtype=np.float64)
        tau_low = np.asarray(model.tau_low_per, dtype=np.float64)
        tau_high = np.asarray(model.tau_high_per, dtype=np.float64)
        call_zone_prob = np.asarray(model.call_zone_prob, dtype=np.float64)
        weights = np.asarray(model.weights, dtype=np.float64)
        call_weight = float(np.clip(getattr(model, "call_weight", 1.0), 0.0, 1.0))

        # reference_continuous is (M, K+1) by construction; drop the trailing
        # Unknown column when disable_unknown is True so the solver works
        # over K tissues only. Output is padded back to K+1 before returning.
        K_tissues = R.shape[1] - 1
        if self.disable_unknown:
            R = R[:, :K_tissues]

        if marker_subset is not None:
            R = R[marker_subset]
            k = k[marker_subset]
            n = n[marker_subset]
            phi = phi[marker_subset]
            tau_low = tau_low[marker_subset]
            tau_high = tau_high[marker_subset]
            call_zone_prob = call_zone_prob[marker_subset]
            weights = weights[marker_subset]

        valid = (
            np.all(np.isfinite(R), axis=1)
            & np.isfinite(k)
            & np.isfinite(n)
            & np.isfinite(phi)
            & np.isfinite(tau_low)
            & np.isfinite(tau_high)
            & np.all(np.isfinite(call_zone_prob), axis=1)
            & (n > 0)
        )
        K_total = R.shape[1]
        if int(np.sum(valid)) < K_total:
            if self.disable_unknown:
                uniform = np.zeros(K_tissues + 1, dtype=np.float64)
                uniform[:K_tissues] = 1.0 / K_tissues
            else:
                uniform = np.zeros(K_total, dtype=np.float64)
                uniform[-1] = 1.0
            return uniform

        R_v = R[valid]
        k_v = k[valid]
        n_v = n[valid]
        phi_v = phi[valid]
        tl_v = tau_low[valid]
        th_v = tau_high[valid]
        cz_v = call_zone_prob[valid]
        w_obj = weights[valid]

        lam = (
            self._effective_unknown_prior_weight(
                int(np.sum(valid)), sample_id=model.sample_id
            )
            if not self.disable_unknown
            else 0.0
        )
        eps_prior = 1e-9

        def neg_ll(w: np.ndarray) -> float:
            ll_vec = hierarchical_marker_log_likelihood(
                w=w,
                reference_continuous=R_v,
                k=k_v,
                n=n_v,
                phi=phi_v,
                tau_low=tl_v,
                tau_high=th_v,
                call_zone_prob=cz_v,
                call_weight=call_weight,
            )
            nll = -float(np.sum(w_obj * ll_vec))
            if lam > 0.0:
                nll += -lam * float(np.log(1.0 - w[-1] + eps_prior))
            return nll

        # Start with scipy's internal finite-difference Jacobian for the
        # hierarchical model; this is robust and keeps implementation simple.
        w_hat = self._run_slsqp(
            neg_ll, None, K_total, sample_id=model.sample_id
        )
        return self._pad_unknown(w_hat, K_tissues)

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
            kwargs = {
                "method": "SLSQP",
                "bounds": bounds,
                "constraints": constraints,
                "options": {"maxiter": self.max_iter, "ftol": self.tol},
            }
            if neg_grad is not None:
                kwargs["jac"] = neg_grad
            res = minimize(
                neg_ll,
                w0,
                **kwargs,
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
    def _augmented_reference(
        reference: "ReferencePanel",
        disable_unknown: bool = False,
    ) -> np.ndarray:
        """Append a column of 0.5 for the always-on unknown component.

        When ``disable_unknown=True``, returns the raw reference methylation
        matrix without appending the Unknown column; the solver then works
        over K cell types only.
        """
        ref = np.asarray(reference.methylation, dtype=np.float64)
        if disable_unknown:
            return ref.copy()
        unknown_col = np.full((ref.shape[0], 1), UNKNOWN_PROFILE, dtype=np.float64)
        return np.hstack([ref, unknown_col])

    def _pad_unknown(self, w: np.ndarray, n_cell_types: int) -> np.ndarray:
        """Pad a K-length weight vector to K+1 by appending 0 for Unknown.

        When ``disable_unknown`` is True the solver returns a length-K
        vector; callers expect K+1 (Unknown is the last slot). Pad with 0
        so downstream code can always index ``w[-1]`` as the Unknown
        fraction without branching.
        """
        if not self.disable_unknown:
            return w
        if w.shape[0] == n_cell_types + 1:
            # Already K+1 — happens when callers pre-emptively padded.
            return w
        if w.shape[0] != n_cell_types:
            raise ValueError(
                f"expected weight vector of length {n_cell_types} or "
                f"{n_cell_types + 1}, got {w.shape[0]}"
            )
        out = np.zeros(n_cell_types + 1, dtype=np.float64)
        out[:n_cell_types] = w
        return out


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
        binarization_model: str = "legacy",
        hierarchical_quadrature_points: int = 24,
        unknown_prior_weight: float = 0.0,
        unknown_prior_weight_auto: bool = False,
        unknown_prior_weight_auto_alpha: float = 0.01,
    ):
        self.n_walkers = n_walkers
        self.n_steps = n_steps
        self.burn_in = burn_in
        self.prior_alpha = prior_alpha
        self.seed = seed
        self.binarization_count_weight = max(float(binarization_count_weight), 0.0)
        self.binarization_model = str(binarization_model or "legacy").lower()
        self.hierarchical_quadrature_points = max(
            int(hierarchical_quadrature_points), 2
        )
        # Unknown-prior regularization (same parameterization as MLE).
        # Adds log-prior  +lambda * log(1 - w_u + eps) to the log-posterior,
        # equivalent to a Beta(1, lambda+1) prior on the Unknown component.
        self.unknown_prior_weight = max(float(unknown_prior_weight), 0.0)
        self.unknown_prior_weight_auto = bool(unknown_prior_weight_auto)
        self.unknown_prior_weight_auto_alpha = max(
            float(unknown_prior_weight_auto_alpha), 0.0
        )
        # Suppress duplicate log lines when solve() is called repeatedly
        # for the same sample (e.g. from uncertainty re-estimation paths).
        self._logged_prior_samples: set[str] = set()

    def _effective_unknown_prior_weight(
        self,
        m_valid: int,
        sample_id: str | None = None,
    ) -> float:
        """Resolve the Unknown-prior lambda for a given valid-marker count.

        Mirrors ``MLEDeconvolver._effective_unknown_prior_weight`` so the
        MLE point estimate and the Bayesian posterior use a consistent prior.
        Logs the resolved lambda once per ``sample_id``.
        """
        if self.unknown_prior_weight_auto:
            lam = self.unknown_prior_weight_auto_alpha * float(max(m_valid, 0))
        else:
            lam = self.unknown_prior_weight
        if (
            lam > 0.0
            and sample_id is not None
            and sample_id not in self._logged_prior_samples
        ):
            if self.unknown_prior_weight_auto:
                log.info(
                    "[%s] Unknown prior (auto, Bayesian): M_valid=%d, "
                    "alpha=%.4f -> lambda=%.3f",
                    sample_id,
                    int(m_valid),
                    self.unknown_prior_weight_auto_alpha,
                    lam,
                )
            else:
                log.info(
                    "[%s] Unknown prior (fixed, Bayesian): lambda=%.3f "
                    "(M_valid=%d)",
                    sample_id,
                    lam,
                    int(m_valid),
                )
            self._logged_prior_samples.add(sample_id)
        return lam

    def _effective_n_walkers(self, k_free: int) -> int:
        """Return a sampler-safe walker count for the current dimensionality.

        emcee's default red-blue move requires at least ``2 * ndim`` walkers
        and behaves best with an even walker count. We enforce:
            n_walkers >= 2 * k_free + 2, and n_walkers is even.
        """
        required = max(2 * int(k_free) + 2, 2)
        n_walkers = max(int(self.n_walkers), required)
        if n_walkers % 2 != 0:
            n_walkers += 1
        return n_walkers

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
            obs_model = str(getattr(model, "binarization_model", "legacy")).lower()
            if self.binarization_model == "hierarchical" or obs_model == "hierarchical":
                return self._solve_hierarchical_binarization(
                    model=model,
                    marker_subset=marker_subset,
                    emcee_module=emcee,
                )
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
        lam_unk = self._effective_unknown_prior_weight(
            int(np.sum(valid)), sample_id=model.sample_id
        )
        eps_prior = 1e-9

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
            # Unknown-component prior: lambda * log(1 - w_unknown + eps).
            # Same parameterization as MLE; adds log-prior contribution.
            if lam_unk > 0.0:
                log_prior += lam_unk * float(np.log(1.0 - w[-1] + eps_prior))
            # Jacobian of softmax (log |det J|) — improper, but cancels for paired comparison
            return log_lik + log_prior

        n_walkers = self._effective_n_walkers(K_free)
        if n_walkers != self.n_walkers:
            log.info(
                "BayesianDeconvolver: increasing n_walkers from %d to %d for ndim=%d",
                self.n_walkers,
                n_walkers,
                K_free,
            )

        rng = np.random.default_rng(self.seed)
        # Initial walkers: small noise around uniform (zero in unconstrained space)
        p0 = rng.normal(0, 0.1, size=(n_walkers, K_free))

        sampler = emcee.EnsembleSampler(n_walkers, K_free, log_posterior)
        sampler.run_mcmc(p0, self.n_steps, progress=False)

        chain = sampler.get_chain(discard=self.burn_in, flat=True)  # (n_keep, K_free)
        # Convert to (K+1)-simplex
        n_keep = chain.shape[0]
        out = np.empty((n_keep, K_total), dtype=np.float64)
        for i in range(n_keep):
            full_z = np.concatenate([chain[i], np.array([0.0])])
            out[i] = softmax(full_z)
        return out

    def _solve_hierarchical_binarization(
        self,
        model,  # BinarizationObservationModel
        marker_subset: np.ndarray | None,
        emcee_module,
    ) -> np.ndarray:
        """Bayesian solve for the hierarchical FinaleMe binarization model."""
        if (
            model.reference_continuous is None
            or model.k is None
            or model.n is None
            or model.dispersion is None
            or model.tau_low_per is None
            or model.tau_high_per is None
            or model.call_zone_prob is None
        ):
            raise DeconvolutionFailedError(
                "BinarizationObservationModel is missing hierarchical fields "
                "(reference_continuous/k/n/dispersion/tau/call_zone_prob)."
            )

        R = np.asarray(model.reference_continuous, dtype=np.float64)
        weights = np.asarray(model.weights, dtype=np.float64)
        k = np.asarray(model.k, dtype=np.float64)
        n = np.asarray(model.n, dtype=np.float64)
        phi = np.asarray(model.dispersion, dtype=np.float64)
        tau_low = np.asarray(model.tau_low_per, dtype=np.float64)
        tau_high = np.asarray(model.tau_high_per, dtype=np.float64)
        call_zone_prob = np.asarray(model.call_zone_prob, dtype=np.float64)
        call_weight = float(np.clip(getattr(model, "call_weight", 1.0), 0.0, 1.0))

        if marker_subset is not None:
            R = R[marker_subset]
            weights = weights[marker_subset]
            k = k[marker_subset]
            n = n[marker_subset]
            phi = phi[marker_subset]
            tau_low = tau_low[marker_subset]
            tau_high = tau_high[marker_subset]
            call_zone_prob = call_zone_prob[marker_subset]

        valid = (
            np.all(np.isfinite(R), axis=1)
            & np.isfinite(k)
            & np.isfinite(n)
            & np.isfinite(phi)
            & np.isfinite(tau_low)
            & np.isfinite(tau_high)
            & np.all(np.isfinite(call_zone_prob), axis=1)
            & (n > 0)
        )
        R_v = R[valid]
        weights_v = weights[valid]
        k_v = k[valid]
        n_v = n[valid]
        phi_v = phi[valid]
        tl_v = tau_low[valid]
        th_v = tau_high[valid]
        cz_v = call_zone_prob[valid]
        if int(R_v.shape[0]) < int(R_v.shape[1]):
            raise DeconvolutionFailedError(
                "Not enough valid markers for hierarchical Bayesian binarization solve."
            )

        K_total = R_v.shape[1]
        K_free = K_total - 1
        lam_unk = self._effective_unknown_prior_weight(
            int(R_v.shape[0]), sample_id=getattr(model, "sample_id", None)
        )
        eps_prior = 1e-9

        def softmax(z: np.ndarray) -> np.ndarray:
            z_max = np.max(z)
            ez = np.exp(z - z_max)
            return ez / ez.sum()

        def log_posterior(z: np.ndarray) -> float:
            full_z = np.concatenate([z, np.array([0.0])])
            w = softmax(full_z)
            ll_vec = hierarchical_marker_log_likelihood(
                w=w,
                reference_continuous=R_v,
                k=k_v,
                n=n_v,
                phi=phi_v,
                tau_low=tl_v,
                tau_high=th_v,
                call_zone_prob=cz_v,
                call_weight=call_weight,
            )
            log_lik = float(np.sum(weights_v * ll_vec))
            if abs(self.prior_alpha - 1.0) < 1e-12:
                log_prior = 0.0
            else:
                log_prior = float(
                    (self.prior_alpha - 1.0) * np.sum(np.log(np.maximum(w, 1e-300)))
                )
            if lam_unk > 0.0:
                log_prior += lam_unk * float(np.log(1.0 - w[-1] + eps_prior))
            return log_lik + log_prior

        n_walkers = self._effective_n_walkers(K_free)
        if n_walkers != self.n_walkers:
            log.info(
                "BayesianDeconvolver (hierarchical): increasing n_walkers from %d to %d for ndim=%d",
                self.n_walkers,
                n_walkers,
                K_free,
            )

        rng = np.random.default_rng(self.seed)
        p0 = rng.normal(0, 0.1, size=(n_walkers, K_free))

        sampler = emcee_module.EnsembleSampler(n_walkers, K_free, log_posterior)
        sampler.run_mcmc(p0, self.n_steps, progress=False)

        chain = sampler.get_chain(discard=self.burn_in, flat=True)
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
        lam_unk = self._effective_unknown_prior_weight(
            int(coef_v.shape[0]), sample_id=getattr(model, "sample_id", None)
        )
        eps_prior = 1e-9

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
            if lam_unk > 0.0:
                log_prior += lam_unk * float(np.log(1.0 - w[-1] + eps_prior))
            return log_lik + log_prior

        n_walkers = self._effective_n_walkers(K_free)
        if n_walkers != self.n_walkers:
            log.info(
                "BayesianDeconvolver (binarization): increasing n_walkers from %d to %d for ndim=%d",
                self.n_walkers,
                n_walkers,
                K_free,
            )

        rng = np.random.default_rng(self.seed)
        p0 = rng.normal(0, 0.1, size=(n_walkers, K_free))

        sampler = emcee_module.EnsembleSampler(n_walkers, K_free, log_posterior)
        sampler.run_mcmc(p0, self.n_steps, progress=False)

        chain = sampler.get_chain(discard=self.burn_in, flat=True)
        n_keep = chain.shape[0]
        out = np.empty((n_keep, K_total), dtype=np.float64)
        for i in range(n_keep):
            full_z = np.concatenate([chain[i], np.array([0.0])])
            out[i] = softmax(full_z)
        return out


__all__ = ["BayesianDeconvolver", "DeconvolutionResult", "MLEDeconvolver", "UNKNOWN_PROFILE"]
