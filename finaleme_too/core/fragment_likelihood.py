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


def _em_from_log_p(
    log_p: np.ndarray,
    max_iter: int = 200,
    tol: float = 1e-6,
) -> tuple[np.ndarray, np.ndarray, float, bool]:
    """Run EM on a precomputed log-likelihood matrix.

    This is the shared building block used by ``solve_full()``, the fragment
    bootstrap, and the LRT/permutation tests.

    Parameters
    ----------
    log_p : ndarray, shape (F, K_total)
        Pre-computed log P(f | j) for all fragments and cell types.
    max_iter : int
        Maximum EM iterations.
    tol : float
        Convergence threshold on max absolute change in proportions.

    Returns
    -------
    w : ndarray, shape (K_total,)
        Estimated mixture proportions.
    gamma : ndarray, shape (F, K_total)
        Per-fragment responsibilities from the final E-step.
    ll : float
        Final log-likelihood Σ_f log Σ_j w_j P(f|j).
    converged : bool
        True if EM reached tolerance before max_iter.
    """
    K_total = log_p.shape[1]
    w = np.full(K_total, 1.0 / K_total, dtype=np.float64)

    gamma = np.empty_like(log_p)
    ll = -np.inf
    converged = False

    for _it in range(max_iter):
        # E-step: numerically stable softmax of log_p + log_w
        with np.errstate(divide="ignore"):
            log_unnorm = log_p + np.log(np.maximum(w, 1e-12))[None, :]
        row_max = np.max(log_unnorm, axis=1, keepdims=True)
        log_unnorm_shifted = log_unnorm - row_max
        unnorm = np.exp(log_unnorm_shifted)
        denom = unnorm.sum(axis=1, keepdims=True)
        denom = np.where(denom > 0, denom, 1.0)
        gamma = unnorm / denom

        # Log-likelihood: Σ_f log(Σ_j w_j P(f|j))
        ll = float(np.sum(np.log(denom.ravel()) + row_max.ravel()))

        # M-step
        w_new = gamma.mean(axis=0)
        w_new = w_new / w_new.sum()

        if np.max(np.abs(w_new - w)) < tol:
            w = w_new
            converged = True
            break
        w = w_new

    return w, gamma, ll, converged


class FragmentLevelDeconvolver:
    """EM deconvolver for ultra-low coverage data."""

    def __init__(
        self,
        max_iter: int = 200,
        tol: float = 1e-6,
        unknown_profile: float = 0.5,
        include_unknown: bool = True,
    ):
        self.max_iter = max_iter
        self.tol = tol
        self.unknown_profile = unknown_profile
        self.include_unknown = include_unknown

    @staticmethod
    def _compute_log_p(
        fragments: list[Fragment],
        reference_augmented: np.ndarray,
    ) -> np.ndarray:
        """Pre-compute log P(f | j) for all fragments and cell types.

        Parameters
        ----------
        fragments : list[Fragment]
            Each fragment has cpg_indices and methylated arrays.
        reference_augmented : ndarray, shape (M, K_total)
            Reference methylation matrix already augmented with the
            unknown column.

        Returns
        -------
        log_p : ndarray, shape (F, K_total)
        """
        R_aug = reference_augmented
        K_total = R_aug.shape[1]
        F = len(fragments)
        log_p = np.zeros((F, K_total), dtype=np.float64)

        for f, frag in enumerate(fragments):
            idx = np.asarray(frag.cpg_indices, dtype=np.int64)
            mask = (idx >= 0) & (idx < R_aug.shape[0])
            if not np.any(mask):
                continue
            idx = idx[mask]
            m = np.asarray(frag.methylated, dtype=np.float64)[mask]
            r = R_aug[idx]  # (L, K_total)
            # Replace NaN values with 0.5 so they contribute zero log-likelihood
            nan_mask = np.isnan(r)
            if nan_mask.any():
                r = r.copy()
                r[nan_mask] = 0.5
            r = np.clip(r, 1e-9, 1.0 - 1e-9)
            ll = m[:, None] * np.log(r) + (1.0 - m[:, None]) * np.log(1.0 - r)
            # Zero out contributions from NaN reference positions
            if nan_mask.any():
                ll[nan_mask] = 0.0
            log_p[f] = ll.sum(axis=0)

        return log_p

    def _augment_reference(self, reference_methylation: np.ndarray) -> np.ndarray:
        """Optionally add the unknown column (uniform at ``self.unknown_profile``)."""
        R = np.asarray(reference_methylation, dtype=np.float64)
        if self.include_unknown:
            return np.hstack(
                [R, np.full((R.shape[0], 1), self.unknown_profile, dtype=np.float64)]
            )
        return R

    def solve_full(
        self,
        fragments: list[Fragment],
        reference_methylation: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, float, np.ndarray, bool]:
        """Run fragment EM and return full results.

        Returns
        -------
        w : ndarray, shape (K+1,) or (K,)
            Estimated mixture proportions. Includes unknown component only
            when ``include_unknown=True``.
        gamma : ndarray, shape (N, K_total)
            Per-fragment responsibilities from the final E-step.
        ll : float
            Final log-likelihood.
        log_p : ndarray, shape (N, K_total)
            Pre-computed log P(f | j) matrix (for bootstrap / LRT reuse).
        converged : bool
            True if EM reached tolerance before max_iter.
        """
        R_aug = self._augment_reference(reference_methylation)
        K_total = R_aug.shape[1]
        F = len(fragments)

        if F == 0:
            uniform = np.zeros(K_total, dtype=np.float64)
            if self.include_unknown:
                uniform[-1] = 1.0
            else:
                uniform[:] = 1.0 / K_total
            empty_gamma = np.zeros((0, K_total), dtype=np.float64)
            empty_log_p = np.zeros((0, K_total), dtype=np.float64)
            return uniform, empty_gamma, 0.0, empty_log_p, False

        log_p = self._compute_log_p(fragments, R_aug)
        w, gamma, ll, converged = _em_from_log_p(log_p, self.max_iter, self.tol)
        return w, gamma, ll, log_p, converged

    def solve(
        self,
        fragments: list[Fragment],
        reference_methylation: np.ndarray,
    ) -> np.ndarray:
        """Return mixture proportions array.

        Backward-compatible wrapper around ``solve_full()``.
        """
        w, _gamma, _ll, _log_p, _converged = self.solve_full(fragments, reference_methylation)
        return w


__all__ = ["Fragment", "FragmentLevelDeconvolver", "_em_from_log_p"]
