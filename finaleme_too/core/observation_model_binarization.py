"""Observation model for FinaleMe binarization mode (v3).

Sibling of ``BetaBinomialModel`` used when ``sample.mode == FINALEME`` and
a trained ``BinarizationParams`` is available. Produces a
``BinarizationObservationModel`` that the refactored MLE solver knows how
to consume.

The key difference from the beta-binomial model: instead of per-marker
``(k, n, dispersion)`` we carry per-marker linear coefficient vectors
``coef[:, i] = b_i ∈ R^(K+1)`` such that the likelihood is

    P(observed_state_i | w) = b_i @ w

Then the negative log-likelihood (for SLSQP) is

    neg_ll(w) = -Σ_i ω_i * log(b_i @ w)

and the gradient is

    neg_grad(w) = -(ω_i / (b_i @ w)) @ coef

Both quantities are computed via numpy matmul in the solver. See
``finaleme_too/core/deconvolution.py`` for the dispatch.

Math reference: ``design/TOO_MATH_FORMULATION_v3.md`` §2B.4 (observation
likelihood) and §3.3 (gradient).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from finaleme_too.config import CoverageTier, MeasurementMode
from finaleme_too.io.methylation_loader import MarkerObservations
from finaleme_too.preprocessing.binarization import (
    STATE_AMBIGUOUS,
    STATE_EXCLUDED,
    STATE_M,
    STATE_U,
    BinarizationParams,
)

if TYPE_CHECKING:
    from finaleme_too.io.reference_panel import ReferencePanel


UNKNOWN_STATE_PRIOR = 0.5  # w_0 contributes P(true_U) = P(true_M) = 0.5


@dataclass(frozen=True)
class BinarizationObservationModel:
    """Per-sample FinaleMe binarization model (v3).

    Attributes
    ----------
    sample_id
        Sample identifier.
    called_state
        uint8 array of per-marker state calls. Only markers where
        ``called_state ∈ {STATE_U, STATE_M}`` contribute to the likelihood.
        Markers where the state is ``STATE_AMBIGUOUS`` or ``STATE_EXCLUDED``
        are filtered out at construction time via ``valid_mask``.
    context_bin
        int32 array of per-marker bin indices into ``BinarizationParams``.
    reference_binary
        Shape ``(M, K+1)``. Reference panel after binarization with a soft
        intermediate region (``r < 0.2 → 0``, ``r > 0.8 → 1``, otherwise
        continuous). The last column (``K``) is the unknown component with
        constant 0.5. Only rows where ``valid_mask[i]`` is True contribute.
    coef
        Shape ``(M, K+1)``. Per-marker linear coefficient vectors
        ``b_i`` such that ``P(called_state_i | w) = coef[i] @ w``.
        Derived from ``reference_binary`` and the per-marker error rates
        ``(ε_U_b, ε_M_b)`` for the marker's bin.

            If called_state_i = U:
                b_i[j] = (1 - ε_U_b) * (1 - r_binary[i, j])
                       + ε_U_b * r_binary[i, j]

            If called_state_i = M:
                b_i[j] = (1 - ε_M_b) * r_binary[i, j]
                       + ε_M_b * (1 - r_binary[i, j])

        The unknown column (``j = K``) uses ``r_binary[i, K] = 0.5``,
        which gives ``b_i[K] = 0.5`` regardless of error rate (the unknown
        component contributes equally to P(true U) and P(true M)).
    weights
        Shape ``(M,)``. Per-marker objective weights ω_i (same formula as
        beta-binomial: ``min(n_i, n_cap) / n_cap``).
    valid_mask
        Shape ``(M_original,)`` bool. True for markers whose ``called_state``
        is U or M and whose context bin is usable. The other arrays on
        this dataclass are already filtered to the subset where
        ``valid_mask`` is True. ``valid_mask`` is kept for downstream
        callers that need to map back to the original marker indices.
    mode, coverage_tier, coverage_cap, kind
        Bookkeeping fields.
    """

    sample_id: str
    called_state: np.ndarray  # uint8 (M_valid,) — U or M only
    context_bin: np.ndarray  # int32 (M_valid,)
    reference_binary: np.ndarray  # float64 (M_valid, K+1)
    coef: np.ndarray  # float64 (M_valid, K+1)
    weights: np.ndarray  # float64 (M_valid,)
    valid_mask: np.ndarray  # bool (M_original,) — trail to original marker index
    mode: MeasurementMode
    coverage_tier: CoverageTier
    coverage_cap: int
    hard_threshold: float | None = None
    # Observation model flavor for FinaleMe:
    #   legacy       -> linear coef path
    #   hierarchical -> joint count+state path
    binarization_model: str = "legacy"
    # Per-marker continuous reference row for hierarchical likelihood.
    # Shape (M_valid, K+1), includes unknown column.
    reference_continuous: np.ndarray | None = None
    # Per-marker thresholds and call-zone probabilities for hierarchical path.
    # tau arrays have shape (M_valid,), call_zone_prob has shape (M_valid, 3)
    # with zones ordered [low, mid, high].
    tau_low_per: np.ndarray | None = None
    tau_high_per: np.ndarray | None = None
    call_zone_prob: np.ndarray | None = None
    # Effective call-channel tempering used in the hierarchical model.
    call_weight: float = 1.0
    # Optional count-level fields for hybrid Bayesian inference. These are
    # aligned to the same valid-marker subset as coef/called_state.
    k: np.ndarray | None = None  # int64 (M_valid,)
    n: np.ndarray | None = None  # int64 (M_valid,)
    dispersion: np.ndarray | None = None  # float64 (M_valid,)
    kind: str = "finaleme"

    @property
    def n_markers(self) -> int:
        return int(self.called_state.size)

    @property
    def n_cell_types_total(self) -> int:
        return int(self.coef.shape[1])

    def log_likelihood(self, w: np.ndarray) -> np.ndarray:
        """Vectorized per-marker log-likelihood given a weight vector w.

        ``w`` has shape ``(K+1,)``. Returns an ``(M_valid,)`` array of
        per-marker log-likelihood contributions (already per-marker, NOT
        weighted by ω_i — multiply by ``self.weights`` for the objective).
        """
        if (
            self.binarization_model == "hierarchical"
            and self.reference_continuous is not None
            and self.k is not None
            and self.n is not None
            and self.dispersion is not None
            and self.tau_low_per is not None
            and self.tau_high_per is not None
            and self.call_zone_prob is not None
        ):
            from finaleme_too.core.deconvolution import (
                hierarchical_marker_log_likelihood,
            )

            return hierarchical_marker_log_likelihood(
                w=np.asarray(w, dtype=np.float64),
                reference_continuous=np.asarray(self.reference_continuous, dtype=np.float64),
                k=np.asarray(self.k, dtype=np.float64),
                n=np.asarray(self.n, dtype=np.float64),
                phi=np.asarray(self.dispersion, dtype=np.float64),
                tau_low=np.asarray(self.tau_low_per, dtype=np.float64),
                tau_high=np.asarray(self.tau_high_per, dtype=np.float64),
                call_zone_prob=np.asarray(self.call_zone_prob, dtype=np.float64),
                call_weight=float(np.clip(self.call_weight, 0.0, 1.0)),
            )
        p_obs = self.coef @ w  # (M_valid,)
        p_obs = np.clip(p_obs, 1e-15, 1.0)
        return np.log(p_obs)

    def total_log_likelihood(self, w: np.ndarray) -> float:
        return float(np.sum(self.weights * self.log_likelihood(w)))


def _binarize_reference_panel(
    reference_methylation: np.ndarray,
    low_threshold: float = 0.2,
    high_threshold: float = 0.8,
    hard_threshold: float | None = None,
) -> np.ndarray:
    """Apply the soft reference binarization rule (math doc §2B.2).

    Default (soft) mode:
      * entries below ``low_threshold`` are set to 0 (U)
      * entries above ``high_threshold`` are set to 1 (M)
      * intermediate values are kept as-is to encode partial support

    Hard-threshold mode (when ``hard_threshold`` is set):
      * entries ``< hard_threshold`` become 0 (U)
      * entries ``>= hard_threshold`` become 1 (M)
      * no intermediate values remain for finite entries

    Returns a new ``float64`` array of the same shape.
    """
    ref = np.asarray(reference_methylation, dtype=np.float64)
    if hard_threshold is not None:
        thr = float(hard_threshold)
        out = ref.copy()
        finite = np.isfinite(ref)
        out[finite & (ref < thr)] = 0.0
        out[finite & (ref >= thr)] = 1.0
        return out

    out = ref.copy()
    out[ref < low_threshold] = 0.0
    out[ref > high_threshold] = 1.0
    return out


class BinarizationModel:
    """Build a ``BinarizationObservationModel`` from raw FinaleMe observations.

    The caller must have already run ``apply_binarization`` on the
    ``MarkerObservations`` (so ``obs.called_state`` and ``obs.context_bin``
    are populated). This builder then:

      1. Filters markers where the state is ``Ambiguous`` or ``Excluded``
         (those don't contribute to the likelihood).
      2. Binarizes the reference panel with the soft intermediate rule.
      3. Precomputes the per-marker linear coefficient matrix ``coef``
         that the SLSQP solver's objective/gradient consumes.
      4. Computes the per-marker coverage weight ``ω_i``.

    Mode / tier / coverage_cap plumbing matches ``BetaBinomialModel.build``
    so the pipeline dispatcher can treat the two builders interchangeably.
    """

    def __init__(
        self,
        hard_threshold: float | None = None,
        binarization_model: str = "legacy",
        call_weight_override: float | None = None,
        learned_threshold_from_params: bool = False,
        learned_threshold_use_all_bins: bool = False,
    ):
        # Optional universal hard threshold used by
        # ``finaleme-too run --binarizeThreshold``.
        self.hard_threshold = hard_threshold
        self.binarization_model = str(binarization_model or "legacy").lower()
        self.call_weight_override = call_weight_override
        self.learned_threshold_from_params = bool(learned_threshold_from_params)
        self.learned_threshold_use_all_bins = bool(learned_threshold_use_all_bins)

    def build(
        self,
        obs: MarkerObservations,
        binarization: BinarizationParams,
        reference: "ReferencePanel",
        tier: CoverageTier = CoverageTier.HIGH,
        coverage_cap: int = 50,
        balance_weights: np.ndarray | None = None,
    ) -> BinarizationObservationModel:
        if obs.called_state is None or obs.context_bin is None:
            raise ValueError(
                f"Sample {obs.sample_id} has no called_state/context_bin — "
                "did you forget to call apply_binarization() before build()?"
            )

        called_state = np.asarray(obs.called_state, dtype=np.uint8)
        context_bin = np.asarray(obs.context_bin, dtype=np.int32)
        n_original = called_state.size

        # Valid markers:
        #   legacy       -> called U/M and usable
        #   hierarchical -> any non-excluded call (U/M/A) and usable
        usable_mask = binarization.usable[context_bin]
        if self.learned_threshold_from_params and self.learned_threshold_use_all_bins:
            # Optional threshold-from-params "all bins" mode.
            usable_mask = np.ones_like(usable_mask, dtype=bool)
        if self.binarization_model == "hierarchical":
            called_valid = called_state != STATE_EXCLUDED
        else:
            called_valid = (called_state == STATE_U) | (called_state == STATE_M)
        valid_mask = usable_mask & called_valid

        # If absolutely nothing is valid, return an empty model — the solver
        # will see zero markers and fall back to the uniform (all-unknown)
        # result via its existing empty-input guard.
        if int(np.sum(valid_mask)) == 0:
            empty_k = 0
            K = int(reference.n_cell_types)
            return BinarizationObservationModel(
                sample_id=obs.sample_id,
                called_state=np.zeros(0, dtype=np.uint8),
                context_bin=np.zeros(0, dtype=np.int32),
                reference_binary=np.zeros((0, K + 1), dtype=np.float64),
                coef=np.zeros((0, K + 1), dtype=np.float64),
                weights=np.zeros(0, dtype=np.float64),
                valid_mask=valid_mask,
                mode=obs.mode,
                coverage_tier=tier,
                coverage_cap=coverage_cap,
                binarization_model=self.binarization_model,
                reference_continuous=np.zeros((0, K + 1), dtype=np.float64),
                tau_low_per=np.zeros(0, dtype=np.float64),
                tau_high_per=np.zeros(0, dtype=np.float64),
                call_zone_prob=np.zeros((0, 3), dtype=np.float64),
                call_weight=binarization.resolve_call_weight(self.call_weight_override),
                k=np.zeros(empty_k, dtype=np.int64),
                n=np.zeros(empty_k, dtype=np.int64),
                dispersion=np.zeros(empty_k, dtype=np.float64),
            )

        # Binarize the reference panel (full shape M_original × K) and
        # append the 0.5 unknown column.
        if self.hard_threshold is not None:
            ref_binary_full = _binarize_reference_panel(
                reference.methylation,
                hard_threshold=self.hard_threshold,
            )
        elif self.learned_threshold_from_params:
            # Per-marker binary reference thresholding using learned tau_high
            # from the marker's context bin. This keeps reference binarization
            # consistent with learned-threshold sample calling.
            thr_full = np.asarray(binarization.tau_high, dtype=np.float64)[context_bin]
            ref_arr = np.asarray(reference.methylation, dtype=np.float64)
            ref_binary_full = np.where(
                np.isfinite(ref_arr) & (ref_arr >= thr_full[:, None]),
                1.0,
                0.0,
            ).astype(np.float64)
        else:
            ref_binary_full = _binarize_reference_panel(
                reference.methylation,
                hard_threshold=None,
            )
        K = ref_binary_full.shape[1]
        unknown_col = np.full((ref_binary_full.shape[0], 1), 0.5, dtype=np.float64)
        ref_full = np.hstack([ref_binary_full, unknown_col])  # (M_original, K+1)
        ref_cont_full = np.hstack(
            [np.asarray(reference.methylation, dtype=np.float64), unknown_col]
        )

        # Subset to valid markers
        ref_valid = ref_full[valid_mask]  # (M_valid, K+1)
        ref_cont_valid = ref_cont_full[valid_mask]  # (M_valid, K+1)
        state_valid = called_state[valid_mask]
        bin_valid = context_bin[valid_mask]

        # Precompute per-marker linear coefficients b_i (math doc §2B.4).
        # For marker i with called state U:
        #     b_i[j] = (1 - ε_U) * (1 - r_binary[i,j]) + ε_U * r_binary[i,j]
        # For called state M:
        #     b_i[j] = (1 - ε_M) * r_binary[i,j] + ε_M * (1 - r_binary[i,j])
        eps_U_per = binarization.eps_U[bin_valid][:, None]  # (M_valid, 1)
        eps_M_per = binarization.eps_M[bin_valid][:, None]

        # coef for every marker assumes "called U" first, then overwrite
        # rows where called_state == M with the M formula. This is cleaner
        # than a Python-level branch per marker.
        coef = (1.0 - eps_U_per) * (1.0 - ref_valid) + eps_U_per * ref_valid
        m_rows = state_valid == STATE_M
        coef[m_rows] = (
            (1.0 - eps_M_per[m_rows]) * ref_valid[m_rows]
            + eps_M_per[m_rows] * (1.0 - ref_valid[m_rows])
        )
        # Ambiguous/Excluded markers are represented as residual probability
        # so legacy reliability calculations remain finite when hierarchical
        # mode includes Ambiguous calls.
        amb_or_exc_rows = (state_valid == STATE_AMBIGUOUS) | (state_valid == STATE_EXCLUDED)
        if np.any(amb_or_exc_rows):
            p_u = (1.0 - eps_U_per[amb_or_exc_rows]) * (1.0 - ref_valid[amb_or_exc_rows]) + (
                eps_U_per[amb_or_exc_rows] * ref_valid[amb_or_exc_rows]
            )
            p_m = (1.0 - eps_M_per[amb_or_exc_rows]) * ref_valid[amb_or_exc_rows] + (
                eps_M_per[amb_or_exc_rows] * (1.0 - ref_valid[amb_or_exc_rows])
            )
            p_a = np.clip(1.0 - p_u - p_m, 1e-12, 1.0)
            coef[amb_or_exc_rows] = p_a

        # Per-marker coverage weight ω_i. For FinaleMe mode we use the
        # per-marker total read count n_i (the number of fragments used by
        # FinaleMe to compute the prediction at this marker), capped at
        # coverage_cap. Matches math doc §3.2.
        n_valid = np.asarray(obs.n, dtype=np.float64)[valid_mask]
        capped = np.minimum(n_valid, float(coverage_cap)) / float(coverage_cap)
        if balance_weights is None:
            balance = np.ones_like(n_valid)
        else:
            balance = np.asarray(balance_weights, dtype=np.float64)[valid_mask]
        weights = capped * balance

        # Count-level view for optional hybrid Bayesian inference. We keep
        # the same FinaleMe tier defaults as the beta-binomial path.
        finaleme_phi_defaults = {
            CoverageTier.HIGH: 15.0,
            CoverageTier.LOW: 8.0,
            CoverageTier.ULTRALOW: 3.0,
        }
        phi_default = finaleme_phi_defaults.get(tier, 8.0)
        k_valid = np.asarray(obs.k, dtype=np.int64)[valid_mask]
        n_valid_i64 = np.asarray(obs.n, dtype=np.int64)[valid_mask]
        phi_valid = np.full(k_valid.size, float(phi_default), dtype=np.float64)
        tau_low_per = np.asarray(binarization.tau_low, dtype=np.float64)[bin_valid]
        tau_high_per = np.asarray(binarization.tau_high, dtype=np.float64)[bin_valid]
        call_zone_prob = binarization.call_zone_probabilities(bin_valid, state_valid)
        call_weight = binarization.resolve_call_weight(self.call_weight_override)

        return BinarizationObservationModel(
            sample_id=obs.sample_id,
            called_state=state_valid.astype(np.uint8),
            context_bin=bin_valid.astype(np.int32),
            reference_binary=ref_valid,
            coef=coef,
            weights=weights,
            valid_mask=valid_mask,
            mode=obs.mode,
            coverage_tier=tier,
            coverage_cap=coverage_cap,
            hard_threshold=self.hard_threshold,
            binarization_model=self.binarization_model,
            reference_continuous=ref_cont_valid,
            tau_low_per=tau_low_per,
            tau_high_per=tau_high_per,
            call_zone_prob=call_zone_prob,
            call_weight=call_weight,
            k=k_valid,
            n=n_valid_i64,
            dispersion=phi_valid,
        )


__all__ = [
    "BinarizationModel",
    "BinarizationObservationModel",
    "UNKNOWN_STATE_PRIOR",
]
