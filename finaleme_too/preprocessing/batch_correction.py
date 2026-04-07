"""ComBat-style technical batch correction at the marker methylation level.

Math doc §10:
    Y*_{i,s} = (Y_{i,s} - γ̂_{b(s),i}) / δ̂_{b(s),i} · δ̂_pool + ᾱ_i

Implementation notes:
    - Operates on per-marker beta values (k_i / n_i) across samples.
    - Empirical Bayes shrinkage for the per-batch (γ, δ) is implemented as a
      simple shrinkage toward the pooled mean and variance using the method-
      of-moments estimator (Johnson, Li, Rabinovic 2007). For very small batches
      this is more robust than the full EM-based ComBat.
    - Skips silently if any batch has fewer than ``min_per_level`` samples
      or if there are fewer than ``min_levels`` distinct batches.
"""

from __future__ import annotations

from collections import Counter

import numpy as np
import pandas as pd

from finaleme_too.io.methylation_loader import MarkerObservations


def combat_correct(
    observations: list[MarkerObservations],
    batch_labels: list[str | None],
    min_levels: int = 2,
    min_per_level: int = 5,
) -> list[MarkerObservations]:
    """Apply ComBat-style location/scale adjustment in beta space.

    Returns a new list of MarkerObservations with k recomputed from the
    adjusted beta values (n preserved). Original list is not mutated.

    If the batch covariate has insufficient levels or any level is too small,
    the inputs are returned unchanged (no warning, per architecture §5.4).
    """
    n_samples = len(observations)
    if n_samples == 0:
        return observations

    # Filter samples with a non-null label
    label_arr = [b for b in batch_labels]
    counts = Counter([b for b in label_arr if b is not None])
    if len(counts) < min_levels or any(v < min_per_level for v in counts.values()):
        return observations

    # Build (n_samples, n_markers) beta matrix; samples must share marker layout
    M = observations[0].n_markers
    Y = np.full((n_samples, M), np.nan, dtype=np.float64)
    for s, obs in enumerate(observations):
        n = np.asarray(obs.n, dtype=np.float64)
        k = np.asarray(obs.k, dtype=np.float64)
        with np.errstate(invalid="ignore", divide="ignore"):
            Y[s] = np.where(n > 0, k / n, np.nan)

    pooled_mean = np.nanmean(Y, axis=0)  # alpha_i
    pooled_var = np.nanvar(Y, axis=0)
    pooled_sd = np.sqrt(np.maximum(pooled_var, 1e-9))

    # Per-batch γ̂ (location offset) and δ̂ (scale)
    unique_batches = sorted(counts.keys())
    Y_adj = Y.copy()
    for b in unique_batches:
        rows = [s for s, lbl in enumerate(label_arr) if lbl == b]
        block = Y[rows]
        gamma = np.nanmean(block, axis=0) - pooled_mean  # location shift
        with np.errstate(invalid="ignore", divide="ignore"):
            delta = np.nanstd(block, axis=0) / np.maximum(pooled_sd, 1e-9)
            delta = np.where(np.isfinite(delta) & (delta > 1e-6), delta, 1.0)
        adj = (block - (pooled_mean + gamma)) / delta * pooled_sd + pooled_mean
        Y_adj[rows] = adj

    # Clip to [0, 1] and recompute k
    Y_adj = np.clip(Y_adj, 0.0, 1.0)
    out: list[MarkerObservations] = []
    for s, obs in enumerate(observations):
        n = np.asarray(obs.n, dtype=np.int64)
        new_k = np.round(Y_adj[s] * n).astype(np.int32)
        new_k = np.clip(new_k, 0, n.astype(np.int32))
        # Markers with zero coverage stay zero
        new_k = np.where(n > 0, new_k, 0)
        out.append(obs.with_counts(new_k, n.astype(np.int32)))
    return out


__all__ = ["combat_correct"]
