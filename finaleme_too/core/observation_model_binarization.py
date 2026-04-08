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
        p_obs = self.coef @ w  # (M_valid,)
        p_obs = np.clip(p_obs, 1e-15, 1.0)
        return np.log(p_obs)

    def total_log_likelihood(self, w: np.ndarray) -> float:
        return float(np.sum(self.weights * self.log_likelihood(w)))


def _binarize_reference_panel(
    reference_methylation: np.ndarray,
    low_threshold: float = 0.2,
    high_threshold: float = 0.8,
) -> np.ndarray:
    """Apply the soft reference binarization rule (math doc §2B.2).

    Entries below ``low_threshold`` are set to 0 (U), entries above
    ``high_threshold`` are set to 1 (M), and intermediate values are kept
    as-is to encode partial support for both states.

    Returns a new ``float64`` array of the same shape.
    """
    ref = np.asarray(reference_methylation, dtype=np.float64)
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

        # Valid markers: called U or M AND in a usable bin.
        usable_mask = binarization.usable[context_bin]
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
            )

        # Binarize the reference panel (full shape M_original × K) and
        # append the 0.5 unknown column.
        ref_binary_full = _binarize_reference_panel(reference.methylation)
        K = ref_binary_full.shape[1]
        unknown_col = np.full((ref_binary_full.shape[0], 1), 0.5, dtype=np.float64)
        ref_full = np.hstack([ref_binary_full, unknown_col])  # (M_original, K+1)

        # Subset to valid markers
        ref_valid = ref_full[valid_mask]  # (M_valid, K+1)
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
        )


__all__ = [
    "BinarizationModel",
    "BinarizationObservationModel",
    "UNKNOWN_STATE_PRIOR",
]
