"""Per-cell-type reliability p-values (math doc §5).

Two complementary p-values per cell type:
  - p_goodness: chi-squared goodness-of-fit on the cell type's most
    discriminative markers
  - p_detection: bootstrap stability above a noise floor
"""

from __future__ import annotations

import numpy as np
from scipy.stats import chi2

from finaleme_too.core.deconvolution import UNKNOWN_PROFILE
from finaleme_too.core.observation_model import ObservationModel


def compute_p_goodness(
    w_hat: np.ndarray,
    reference_methylation: np.ndarray,
    observation: ObservationModel,
    cell_type_index: int,
    top_n: int = 50,
) -> float:
    """χ² goodness-of-fit on the top-N most discriminative markers for one cell type.

    Discriminativeness for cell type j is measured as
        |r_{ij} - mean_{j' != j} r_{ij'}|
    """
    R = np.asarray(reference_methylation, dtype=np.float64)
    K = R.shape[1]
    if cell_type_index >= K:
        return float("nan")

    # Discrimination score
    target = R[:, cell_type_index]
    others = np.delete(R, cell_type_index, axis=1)
    bg_mean = np.mean(others, axis=1)
    score = np.abs(target - bg_mean)

    valid = (observation.n > 0) & np.isfinite(score)
    score_valid = np.where(valid, score, -np.inf)
    top_n = min(top_n, int(np.sum(valid)))
    if top_n < 5:
        return float("nan")
    top_idx = np.argpartition(-score_valid, top_n - 1)[:top_n]

    R_full = np.hstack([R, np.full((R.shape[0], 1), UNKNOWN_PROFILE)])
    mu_pred = (R_full @ w_hat)[top_idx]
    mu_pred = np.clip(mu_pred, 1e-9, 1.0 - 1e-9)
    k = observation.k[top_idx].astype(np.float64)
    n = observation.n[top_idx].astype(np.float64)
    phi = observation.dispersion[top_idx]

    expected_k = mu_pred * n
    var = n * mu_pred * (1.0 - mu_pred) * (n + phi) / (1.0 + phi)
    var = np.maximum(var, 1e-10)
    chi2_stat = float(np.sum((k - expected_k) ** 2 / var))
    df = max(top_n - 1, 1)
    return float(1.0 - chi2.cdf(chi2_stat, df))


def compute_p_detection(
    bootstrap_proportions_j: np.ndarray, noise_floor: float = 0.001
) -> float:
    """Fraction of bootstrap replicates above a noise floor."""
    arr = np.asarray(bootstrap_proportions_j, dtype=np.float64)
    if arr.size == 0:
        return float("nan")
    above = float(np.mean(arr >= noise_floor))
    return above


def assign_reliability(p_goodness: float, p_detection: float) -> str:
    """Assign HIGH/MODERATE/LOW/UNRELIABLE per math doc §5.3 table."""
    if np.isnan(p_goodness) and np.isnan(p_detection):
        return "UNRELIABLE"
    if np.isnan(p_goodness):
        p_goodness = 1.0
    if np.isnan(p_detection):
        p_detection = 0.0
    if p_goodness > 0.05 and p_detection > 0.95:
        return "HIGH"
    if p_goodness > 0.05 and 0.50 <= p_detection <= 0.95:
        return "MODERATE"
    if p_goodness < 0.01 and p_detection < 0.50:
        return "UNRELIABLE"
    return "LOW"


__all__ = ["assign_reliability", "compute_p_detection", "compute_p_goodness"]
