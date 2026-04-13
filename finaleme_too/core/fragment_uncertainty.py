"""Fragment-level statistical testing: bootstrap CI, LRT, permutation (design §3.5.7-3.5.8).

All functions operate on a precomputed ``log_p`` matrix (N fragments x K+1
cell types) so that the expensive per-fragment likelihood computation is done
only once by ``FragmentLevelDeconvolver.solve_full()``. Bootstrap/LRT/
permutation replicates resample or modify rows/columns of this matrix and
re-run the lightweight EM via ``_em_from_log_p()``.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
from scipy.stats import chi2

from finaleme_too.core.fragment_likelihood import _em_from_log_p
from finaleme_too.utils.parallel import parallel_map

log = logging.getLogger(__name__)


@dataclass
class FragmentBootstrapResult:
    """Result of fragment-resampling bootstrap."""

    proportions_samples: np.ndarray  # (B, K+1)
    ci_lower: np.ndarray  # (K+1,)
    ci_upper: np.ndarray  # (K+1,)
    point_estimate: np.ndarray  # (K+1,)


class FragmentBootstrapCI:
    """Bootstrap CI by resampling fragment rows of the log-likelihood matrix.

    Design doc §3.5.7: for each replicate, draw N fragment indices with
    replacement, re-run EM on the resampled log_p rows.
    """

    def __init__(
        self,
        n_iterations: int = 1000,
        ci_level: float = 0.95,
        seed: int | None = None,
    ):
        self.n_iterations = n_iterations
        self.ci_level = ci_level
        self.seed = seed

    def estimate(
        self,
        log_p_matrix: np.ndarray,
        max_iter: int = 200,
        tol: float = 1e-6,
        n_jobs: int = 1,
    ) -> FragmentBootstrapResult:
        """Run bootstrap on precomputed log-likelihood matrix.

        Parameters
        ----------
        log_p_matrix : ndarray, shape (N, K+1)
        max_iter : int
            Max EM iterations per replicate.
        tol : float
            EM convergence threshold.
        n_jobs : int
            Parallel workers (threading backend recommended).

        Returns
        -------
        FragmentBootstrapResult
        """
        N, K_total = log_p_matrix.shape
        if N == 0:
            zero = np.zeros(K_total, dtype=np.float64)
            return FragmentBootstrapResult(
                proportions_samples=np.tile(zero, (1, 1)),
                ci_lower=zero.copy(),
                ci_upper=zero.copy(),
                point_estimate=zero.copy(),
            )

        rng = np.random.default_rng(self.seed)
        # Pre-generate all index arrays for reproducibility
        all_indices = [
            rng.integers(0, N, size=N) for _ in range(self.n_iterations)
        ]

        def _run_one(indices: np.ndarray) -> np.ndarray:
            resampled = log_p_matrix[indices]
            w, _gamma, _ll = _em_from_log_p(resampled, max_iter, tol)
            return w

        samples = parallel_map(
            _run_one,
            all_indices,
            n_jobs=n_jobs,
            backend="threading",
        )
        proportions_samples = np.array(samples, dtype=np.float64)

        alpha = 1.0 - self.ci_level
        ci_lower = np.quantile(proportions_samples, alpha / 2, axis=0)
        ci_upper = np.quantile(proportions_samples, 1.0 - alpha / 2, axis=0)
        point_estimate = proportions_samples.mean(axis=0)

        return FragmentBootstrapResult(
            proportions_samples=proportions_samples,
            ci_lower=ci_lower,
            ci_upper=ci_upper,
            point_estimate=point_estimate,
        )


def fragment_lrt_pvalues(
    log_p_matrix: np.ndarray,
    ll_full: float,
    max_iter: int = 200,
    tol: float = 1e-6,
    n_jobs: int = 1,
) -> np.ndarray:
    """Likelihood-ratio test p-value per cell type (design §3.5.8).

    For each cell type t, remove column t from log_p, re-run EM on the
    remaining K columns, and compute

        Λ = 2 * max(0, ll_full - ll_reduced)
        p = chi2.sf(Λ, df=1)

    Parameters
    ----------
    log_p_matrix : ndarray, shape (N, K+1)
        Includes the unknown component as the last column.
    ll_full : float
        Log-likelihood from the full (K+1)-component model.
    max_iter, tol : EM parameters.
    n_jobs : int
        Parallel workers.

    Returns
    -------
    p_values : ndarray, shape (K+1,)
    """
    K_total = log_p_matrix.shape[1]

    def _lrt_one(t: int) -> float:
        # Remove column t
        cols = list(range(K_total))
        cols.pop(t)
        log_p_reduced = log_p_matrix[:, cols]
        _w, _gamma, ll_reduced = _em_from_log_p(log_p_reduced, max_iter, tol)
        lr_stat = 2.0 * max(0.0, ll_full - ll_reduced)
        return float(chi2.sf(lr_stat, df=1))

    pvals = parallel_map(
        _lrt_one,
        list(range(K_total)),
        n_jobs=n_jobs,
        backend="threading",
    )
    return np.array(pvals, dtype=np.float64)


def fragment_permutation_pvalues(
    log_p_matrix: np.ndarray,
    w_observed: np.ndarray,
    n_permutations: int = 1000,
    max_iter: int = 200,
    tol: float = 1e-6,
    seed: int | None = None,
    n_jobs: int = 1,
) -> np.ndarray:
    """Permutation-based p-value per cell type (design §3.5.8 fallback).

    For each cell type t, independently permute the rows of column t in
    log_p, re-run EM, and record pi_t_null. The p-value is the fraction
    of null estimates >= the observed estimate.

    Parameters
    ----------
    log_p_matrix : ndarray, shape (N, K+1)
    w_observed : ndarray, shape (K+1,)
    n_permutations : int
    max_iter, tol : EM parameters.
    seed : int or None
    n_jobs : int

    Returns
    -------
    p_values : ndarray, shape (K+1,)
    """
    K_total = log_p_matrix.shape[1]
    N = log_p_matrix.shape[0]
    rng = np.random.default_rng(seed)

    p_values = np.full(K_total, np.nan, dtype=np.float64)

    for t in range(K_total):
        # Pre-generate permutation indices
        perm_indices = [rng.permutation(N) for _ in range(n_permutations)]

        def _perm_one(indices: np.ndarray) -> float:
            log_p_perm = log_p_matrix.copy()
            log_p_perm[:, t] = log_p_perm[indices, t]
            w, _gamma, _ll = _em_from_log_p(log_p_perm, max_iter, tol)
            return w[t]

        null_pi_t = parallel_map(
            _perm_one,
            perm_indices,
            n_jobs=n_jobs,
            backend="threading",
        )
        null_arr = np.array(null_pi_t, dtype=np.float64)
        p_values[t] = (np.sum(null_arr >= w_observed[t]) + 1) / (n_permutations + 1)

    return p_values


def fragment_bh_correction(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-value correction, preserving NaN entries.

    Applies BH across all T+1 components including the Unknown component
    (design §3.5.8: "Apply BH procedure across all T+1 components").
    """
    pvals = np.asarray(pvals, dtype=np.float64)
    out = np.full(pvals.shape, np.nan, dtype=np.float64)
    finite = np.isfinite(pvals)
    if not np.any(finite):
        return out
    p = np.clip(pvals[finite], 0.0, 1.0)
    m = p.size
    order = np.argsort(p)
    ranked = p[order]
    ranks = np.arange(1, m + 1, dtype=np.float64)
    adj_ranked = ranked * (m / ranks)
    adj_ranked = np.minimum.accumulate(adj_ranked[::-1])[::-1]
    adj_ranked = np.clip(adj_ranked, 0.0, 1.0)
    adj = np.empty_like(adj_ranked)
    adj[order] = adj_ranked
    out[finite] = adj
    return out


__all__ = [
    "FragmentBootstrapResult",
    "FragmentBootstrapCI",
    "fragment_lrt_pvalues",
    "fragment_permutation_pvalues",
    "fragment_bh_correction",
]
