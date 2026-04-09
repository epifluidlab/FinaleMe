"""ComBat-style technical batch correction at the marker methylation level.

Math doc §10:
    Y*_{i,s} = (Y_{i,s} - γ̂_{b(s),i}) / δ̂_{b(s),i} · δ̂_pool + ᾱ_i

Implementation notes:
    - WGBS mode: operates on per-marker beta values (k_i / n_i) across
      samples and writes the corrected betas back via obs.with_counts(k, n).
    - FinaleMe mode: operates on per-marker ``predicted_beta`` values BEFORE
      binarization. Per math doc §10, batch-shifted FinaleMe predictions
      would otherwise cause systematic miscalls; the corrected predictions
      are the inputs that ``apply_binarization`` then classifies into
      U / M / Ambiguous / Excluded states.
    - Empirical Bayes shrinkage for the per-batch (γ, δ) is the method-of-
      moments variant (Johnson, Li, Rabinovic 2007). For very small batches
      this is more robust than the full EM-based ComBat.
    - Skips silently if any batch has fewer than ``min_per_level`` samples
      or if there are fewer than ``min_levels`` distinct batches.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import replace as dc_replace

import numpy as np
import pandas as pd

from finaleme_too.io.methylation_loader import MarkerObservations


def _adjust_beta_matrix(
    Y: np.ndarray,
    batch_labels: list[str | None],
    min_levels: int,
    min_per_level: int,
) -> np.ndarray | None:
    """Shared ComBat-style location/scale adjustment on a beta matrix.

    ``Y`` is shape ``(n_samples, n_markers)`` with NaN for missing. Returns
    the adjusted matrix clipped to ``[0, 1]``, or ``None`` when the batch
    covariate is insufficient to run the correction.
    """
    counts = Counter([b for b in batch_labels if b is not None])
    if len(counts) < min_levels or any(v < min_per_level for v in counts.values()):
        return None

    pooled_mean = np.nanmean(Y, axis=0)
    pooled_var = np.nanvar(Y, axis=0)
    pooled_sd = np.sqrt(np.maximum(pooled_var, 1e-9))

    Y_adj = Y.copy()
    for b in sorted(counts.keys()):
        rows = [s for s, lbl in enumerate(batch_labels) if lbl == b]
        block = Y[rows]
        gamma = np.nanmean(block, axis=0) - pooled_mean
        with np.errstate(invalid="ignore", divide="ignore"):
            delta = np.nanstd(block, axis=0) / np.maximum(pooled_sd, 1e-9)
            delta = np.where(np.isfinite(delta) & (delta > 1e-6), delta, 1.0)
        adj = (block - (pooled_mean + gamma)) / delta * pooled_sd + pooled_mean
        Y_adj[rows] = adj

    return np.clip(Y_adj, 0.0, 1.0)


def combat_correct(
    observations: list[MarkerObservations],
    batch_labels: list[str | None],
    min_levels: int = 2,
    min_per_level: int = 5,
) -> list[MarkerObservations]:
    """Apply ComBat-style location/scale adjustment in beta space (WGBS mode).

    Returns a new list of ``MarkerObservations`` with ``k`` recomputed from
    the adjusted beta values (``n`` preserved). Operates on the per-marker
    ``k/n`` betas — the right path for WGBS mode, where the observation
    model reads directly from the integer counts.

    Use ``combat_correct_predicted_beta`` for FinaleMe mode, which
    adjusts ``predicted_beta`` (before binarization) instead.

    If the batch covariate has insufficient levels or any level is too
    small, the inputs are returned unchanged (no warning, per
    architecture §5.4).
    """
    n_samples = len(observations)
    if n_samples == 0:
        return observations

    # Build (n_samples, n_markers) beta matrix from k/n.
    M = observations[0].n_markers
    Y = np.full((n_samples, M), np.nan, dtype=np.float64)
    for s, obs in enumerate(observations):
        n = np.asarray(obs.n, dtype=np.float64)
        k = np.asarray(obs.k, dtype=np.float64)
        with np.errstate(invalid="ignore", divide="ignore"):
            Y[s] = np.where(n > 0, k / n, np.nan)

    Y_adj = _adjust_beta_matrix(Y, batch_labels, min_levels, min_per_level)
    if Y_adj is None:
        return observations

    out: list[MarkerObservations] = []
    for s, obs in enumerate(observations):
        n = np.asarray(obs.n, dtype=np.int64)
        new_k = np.round(Y_adj[s] * n).astype(np.int32)
        new_k = np.clip(new_k, 0, n.astype(np.int32))
        # Markers with zero coverage stay zero
        new_k = np.where(n > 0, new_k, 0)
        out.append(obs.with_counts(new_k, n.astype(np.int32)))
    return out


def combat_correct_predicted_beta(
    observations: list[MarkerObservations],
    batch_labels: list[str | None],
    min_levels: int = 2,
    min_per_level: int = 5,
) -> list[MarkerObservations]:
    """Apply ComBat-style correction to FinaleMe ``predicted_beta`` values.

    The FinaleMe path (math doc §10): batch correction runs on
    ``obs.predicted_beta`` **before** binarization, so the corrected
    predictions are what ``apply_binarization`` then classifies into
    U / M / Ambiguous / Excluded states. Without this, batch-shifted
    FinaleMe predictions would cause systematic miscalls.

    The adjusted predicted_beta is written back via
    ``dataclasses.replace`` so no other fields are disturbed. Samples
    without a ``predicted_beta`` (WGBS mode observations mixed into a
    FinaleMe cohort) pass through unchanged. ``k``, ``n``,
    ``called_state``, and ``context_bin`` are all preserved on the
    output — although in practice this function should be called
    **before** ``apply_binarization`` so ``called_state`` is still
    ``None`` at this point.

    If the batch covariate has insufficient levels or any level is too
    small, the inputs are returned unchanged.
    """
    n_samples = len(observations)
    if n_samples == 0:
        return observations

    # Check that at least one observation has predicted_beta populated —
    # otherwise this is a WGBS-only cohort and we have nothing to do.
    first_with_pred = next(
        (obs for obs in observations if obs.predicted_beta is not None), None
    )
    if first_with_pred is None:
        return observations

    M = first_with_pred.n_markers
    Y = np.full((n_samples, M), np.nan, dtype=np.float64)
    has_pred = np.zeros(n_samples, dtype=bool)
    for s, obs in enumerate(observations):
        if obs.predicted_beta is None:
            continue
        pred = np.asarray(obs.predicted_beta, dtype=np.float64)
        if pred.size != M:
            # Shape mismatch (shouldn't happen in practice — the loader
            # aligns everything to the same marker set). Skip.
            continue
        Y[s] = pred
        has_pred[s] = True

    # Only use labels for samples that actually have predicted_beta.
    effective_labels = [
        lbl if has_pred[s] else None for s, lbl in enumerate(batch_labels)
    ]
    Y_adj = _adjust_beta_matrix(Y, effective_labels, min_levels, min_per_level)
    if Y_adj is None:
        return observations

    out: list[MarkerObservations] = []
    for s, obs in enumerate(observations):
        if not has_pred[s]:
            out.append(obs)
            continue
        new_pred = Y_adj[s].astype(np.float32)
        out.append(dc_replace(obs, predicted_beta=new_pred))
    return out


__all__ = ["combat_correct", "combat_correct_predicted_beta"]
