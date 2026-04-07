"""Per-sample residual variance + cohort-level NMF residual discovery (architecture §9.4)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass
class ResidualAnalysisResult:
    sample_id: str
    unexplained_fraction: float
    mean_residual: float
    residual_sd: float
    qc_flag: str


def compute_residuals_per_sample(
    proportions: np.ndarray,
    reference_methylation: np.ndarray,
    observation_methylation: np.ndarray,
) -> dict:
    """Compute per-sample residual statistics for the unknown component.

    Returns a dict with keys: unexplained_fraction, mean_residual, residual_sd.
    The unexplained fraction is the unknown component (last entry of
    ``proportions``); residuals are computed against the model prediction
    excluding the unknown column.
    """
    K = reference_methylation.shape[1]
    w = np.asarray(proportions, dtype=np.float64)
    if w.size != K + 1:
        raise ValueError("proportions must have length K+1 (with unknown)")
    pred_no_unknown = reference_methylation @ w[:K]
    residuals = observation_methylation - pred_no_unknown
    valid = np.isfinite(residuals)
    return {
        "unexplained_fraction": float(w[-1]),
        "mean_residual": float(np.nanmean(residuals[valid])) if np.any(valid) else float("nan"),
        "residual_sd": float(np.nanstd(residuals[valid])) if np.any(valid) else float("nan"),
    }


def discover_residual_components(
    residual_matrix: np.ndarray,
    n_components: int = 3,
    random_state: int = 0,
) -> dict:
    """Run NMF on a (S, M) matrix of per-sample marker residuals.

    Returns a dict with:
        components: (n_components, M) array — the discovered profiles
        loadings:   (S, n_components) array — per-sample loadings
        explained_variance_ratio: 1D array of length n_components
    """
    from sklearn.decomposition import NMF

    R = np.asarray(residual_matrix, dtype=np.float64)
    # NMF requires non-negative input — clip negative residuals to 0
    R_pos = np.maximum(R, 0.0)
    if R_pos.size == 0 or R_pos.shape[0] < 2:
        return {
            "components": np.zeros((0, R.shape[1])),
            "loadings": np.zeros((R.shape[0], 0)),
            "explained_variance_ratio": np.array([]),
        }
    n_components = max(1, min(n_components, R_pos.shape[0] - 1, R_pos.shape[1]))
    model = NMF(n_components=n_components, init="random", random_state=random_state, max_iter=500)
    W = model.fit_transform(R_pos)
    H = model.components_
    total_var = float(np.var(R_pos))
    if total_var > 0:
        per_comp_var = []
        for c in range(n_components):
            recon_c = np.outer(W[:, c], H[c])
            per_comp_var.append(float(np.var(recon_c)) / total_var)
        ratio = np.array(per_comp_var)
    else:
        ratio = np.zeros(n_components)
    return {
        "components": H,
        "loadings": W,
        "explained_variance_ratio": ratio,
    }


__all__ = [
    "ResidualAnalysisResult",
    "compute_residuals_per_sample",
    "discover_residual_components",
]
