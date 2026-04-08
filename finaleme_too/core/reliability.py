"""Per-cell-type reliability p-values (math doc §5).

Two complementary p-values per cell type:
  - p_goodness: model fit quality at the cell type's most discriminative markers.
    WGBS mode uses a chi-squared goodness-of-fit on the read counts; FinaleMe
    mode (v3 binarization) uses a one-sided binomial concordance test on the
    discrete state calls vs the expected states under the estimated mixture.
  - p_detection: bootstrap stability above a noise floor.

Both follow the convention HIGH = GOOD: a high p_goodness means the data is
consistent with the estimated mixture, NOT that some null hypothesis is
rejected.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import chi2

from finaleme_too.core.deconvolution import UNKNOWN_PROFILE
from finaleme_too.core.observation_model import ObservationModel


def compute_p_goodness(
    w_hat: np.ndarray,
    reference_methylation: np.ndarray,
    observation,
    cell_type_index: int,
    top_n: int = 50,
    binarizer=None,
) -> float:
    """Per-cell-type goodness-of-fit p-value at the top-N discriminative markers.

    Discriminativeness for cell type j is measured as
        |r_{ij} - mean_{j' != j} r_{ij'}|

    Parameters
    ----------
    w_hat
        Estimated mixture of length K+1 (with unknown last).
    reference_methylation
        Full (M, K) reference methylation matrix.
    observation
        Either an ``ObservationModel`` (WGBS / v2 path) or a
        ``BinarizationObservationModel`` (v3 FinaleMe path). The dispatch
        chooses the appropriate p-value formula.
    cell_type_index
        Cell type to score (in [0, K)).
    top_n
        Number of discriminative markers to consider per cell type.
    binarizer
        Optional ``BinarizationParams`` (only used by the FinaleMe path,
        for ``p_expected = 1 - mean(eps_U_b, eps_M_b)``). Ignored for WGBS.

    Returns
    -------
    A p-value in ``[0, 1]``. **High = good fit.** NaN when there are
    fewer than 5 valid markers for the cell type.
    """
    # Local import to avoid a circular dependency at module load.
    from finaleme_too.core.observation_model_binarization import (
        BinarizationObservationModel,
    )

    if isinstance(observation, BinarizationObservationModel):
        return _p_goodness_binarization(
            w_hat=w_hat,
            reference_methylation=reference_methylation,
            observation=observation,
            cell_type_index=cell_type_index,
            top_n=top_n,
            binarizer=binarizer,
        )
    return _p_goodness_betabinom(
        w_hat=w_hat,
        reference_methylation=reference_methylation,
        observation=observation,
        cell_type_index=cell_type_index,
        top_n=top_n,
    )


def _p_goodness_betabinom(
    w_hat: np.ndarray,
    reference_methylation: np.ndarray,
    observation: ObservationModel,
    cell_type_index: int,
    top_n: int = 50,
) -> float:
    """χ² goodness-of-fit on read counts (WGBS / v2 FinaleMe path)."""
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


def _p_goodness_binarization(
    w_hat: np.ndarray,
    reference_methylation: np.ndarray,
    observation,
    cell_type_index: int,
    top_n: int,
    binarizer,
) -> float:
    """Binomial concordance test on state calls (v3 FinaleMe path).

    Math doc §5.1 (FinaleMe variant):

      For cell type j, take its top-D_j discriminative markers (the same
      ranking by |r_{ij} - mean_{j' != j} r_{ij'}|). Of those that
      survived binarization (i.e. ``observation.valid_mask`` is True at
      that marker index), count how many called states agree with the
      expected state under the estimated mixture w_hat:

          expected_state_i = U if P(true_U | w_hat) > 0.5 else M

      Then run a one-sided binomial test:

          p_expected = 1 - mean(eps_U_b, eps_M_b)
          p_goodness = binomtest(n_concordant, n_tested, p_expected,
                                 alternative='less').pvalue

      A high p-value means we cannot reject the null "the model fits at
      least as well as the calibrated error rates", i.e. GOOD fit. A low
      p-value means too few calls agree with the expected state — the
      model is wrong.
    """
    from scipy.stats import binomtest

    from finaleme_too.preprocessing.binarization import STATE_M, STATE_U

    R = np.asarray(reference_methylation, dtype=np.float64)
    M_original, K = R.shape
    if cell_type_index >= K:
        return float("nan")

    # Use the soft binarized reference for both the discriminative score
    # and the expected-state computation, so the score reflects what the
    # binarization model actually sees.
    R_binary = R.copy()
    R_binary[R < 0.2] = 0.0
    R_binary[R > 0.8] = 1.0
    target = R_binary[:, cell_type_index]
    others = np.delete(R_binary, cell_type_index, axis=1)
    bg_mean = np.mean(others, axis=1)
    score = np.abs(target - bg_mean)

    # observation.valid_mask is shape (M_original,). Only markers in
    # valid_mask have entries in observation.coef / called_state.
    valid_mask = observation.valid_mask
    finite_score = np.isfinite(score)
    candidates = valid_mask & finite_score
    if int(np.sum(candidates)) < 5:
        return float("nan")

    # Pick top-N by discrimination score among the candidates.
    score_eligible = np.where(candidates, score, -np.inf)
    n_top = min(top_n, int(np.sum(candidates)))
    top_idx_original = np.argpartition(-score_eligible, n_top - 1)[:n_top]

    # Translate from original marker indices to indices in the filtered
    # observation arrays. valid_mask is bool over M_original; the cumulative
    # count gives the position in the filtered space minus 1 for True
    # entries.
    valid_to_filtered = np.cumsum(valid_mask) - 1
    filtered_idx = valid_to_filtered[top_idx_original]

    # Per-marker called states + bin indices
    called = np.asarray(observation.called_state, dtype=np.uint8)[filtered_idx]
    bins = np.asarray(observation.context_bin, dtype=np.int64)[filtered_idx]

    # Expected state under w_hat for each top marker. P(true_U | w) =
    # Σ_j w_j * (1 - r_binary[i,j]) + 0.5 * w_0. The augmented R already
    # has the unknown column with value 0.5 for "true_U(unknown)" because
    # the math doc says w_0 contributes 0.5 to both states.
    K_total = K + 1
    if w_hat.size != K_total:
        return float("nan")
    R_aug_U = np.hstack(
        [1.0 - R_binary, np.full((M_original, 1), 0.5, dtype=np.float64)]
    )
    p_true_U = (R_aug_U @ w_hat)[top_idx_original]
    expected_state = np.where(p_true_U > 0.5, STATE_U, STATE_M).astype(np.uint8)
    concordant = called == expected_state
    n_tested = int(concordant.size)
    n_concordant = int(np.sum(concordant))

    if n_tested == 0:
        return float("nan")

    # Per-bin error rates: average over the bins of these top markers.
    if binarizer is not None:
        bin_eps = 0.5 * (
            np.asarray(binarizer.eps_U)[bins] + np.asarray(binarizer.eps_M)[bins]
        )
        mean_eps = float(np.clip(np.mean(bin_eps), 0.0, 0.49))
        p_expected = 1.0 - mean_eps
    else:
        # Without a binarizer we have no error-rate estimate; fall back to
        # 0.5 (chance), which makes the test almost always pass.
        p_expected = 0.5

    # Clip to (0, 1) for binomtest stability
    p_expected = float(np.clip(p_expected, 1e-9, 1.0 - 1e-9))

    try:
        result = binomtest(
            n_concordant, n_tested, p=p_expected, alternative="less"
        )
        return float(result.pvalue)
    except Exception:  # noqa: BLE001
        return float("nan")


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
