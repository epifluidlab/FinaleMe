"""Fragment-level EM deconvolver for ULTRALOW coverage (math doc §11).

Each fragment f covers a set of CpGs {c_1, ..., c_L} with binary methylation
calls m_l ∈ {0, 1}. The likelihood under cell type j is

    P(f | j) = Π_l r_{c_l, j}^{m_l} * (1 - r_{c_l, j})^{1 - m_l}

The EM iteration is

    E:  γ_{f,j} = w_j P(f | j) / Σ_{j'} w_{j'} P(f | j')
    M:  w_j = (1/F) Σ_f γ_{f,j}

This module is intentionally lightweight: it accepts a list of fragments
(each represented as ``(cpg_indices, methylated_calls)``) and a reference
matrix indexed by CpG. Optional numba JIT acceleration is loaded lazily.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

log = logging.getLogger(__name__)


@dataclass
class Fragment:
    """Single fragment over a set of CpGs."""

    cpg_indices: np.ndarray  # int64 indices into the reference rows
    methylated: np.ndarray  # uint8 (0/1) calls aligned to cpg_indices


class FragmentLevelDeconvolver:
    """EM deconvolver for ultra-low coverage data."""

    def __init__(
        self,
        max_iter: int = 200,
        tol: float = 1e-6,
        unknown_profile: float = 0.5,
    ):
        self.max_iter = max_iter
        self.tol = tol
        self.unknown_profile = unknown_profile

    def solve(
        self,
        fragments: list[Fragment],
        reference_methylation: np.ndarray,
    ) -> np.ndarray:
        """Return (K+1,) mixture proportions including the unknown component."""
        R = np.asarray(reference_methylation, dtype=np.float64)
        # Augment with unknown column at 0.5
        R_aug = np.hstack(
            [R, np.full((R.shape[0], 1), self.unknown_profile, dtype=np.float64)]
        )
        K_total = R_aug.shape[1]

        F = len(fragments)
        if F == 0:
            uniform = np.zeros(K_total, dtype=np.float64)
            uniform[-1] = 1.0
            return uniform

        # Pre-compute log P(f | j) for all fragments and cell types
        log_p = np.zeros((F, K_total), dtype=np.float64)
        for f, frag in enumerate(fragments):
            idx = np.asarray(frag.cpg_indices, dtype=np.int64)
            mask = (idx >= 0) & (idx < R_aug.shape[0])
            if not np.any(mask):
                continue
            idx = idx[mask]
            m = np.asarray(frag.methylated, dtype=np.float64)[mask]
            r = np.clip(R_aug[idx], 1e-9, 1.0 - 1e-9)  # (L, K_total)
            ll = m[:, None] * np.log(r) + (1.0 - m[:, None]) * np.log(1.0 - r)
            log_p[f] = ll.sum(axis=0)

        # Initialize w uniform
        w = np.full(K_total, 1.0 / K_total, dtype=np.float64)

        for it in range(self.max_iter):
            # E-step: compute responsibilities (numerically stable softmax of log p + log w)
            with np.errstate(divide="ignore"):
                log_unnorm = log_p + np.log(np.maximum(w, 1e-12))[None, :]
            row_max = np.max(log_unnorm, axis=1, keepdims=True)
            log_unnorm = log_unnorm - row_max
            unnorm = np.exp(log_unnorm)
            denom = unnorm.sum(axis=1, keepdims=True)
            denom = np.where(denom > 0, denom, 1.0)
            gamma = unnorm / denom

            # M-step
            w_new = gamma.mean(axis=0)
            w_new = w_new / w_new.sum()

            if np.max(np.abs(w_new - w)) < self.tol:
                w = w_new
                break
            w = w_new

        return w


__all__ = ["Fragment", "FragmentLevelDeconvolver"]
