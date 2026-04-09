"""Per-cell-type reliability metrics (math doc §5 + v3 reliability update).

Reliability combines:
  - p_detection: bootstrap stability above a noise floor.
  - likelihood_score: weighted per-marker log-likelihood gain over an
    ablated-null model.
  - p_likelihood: likelihood-ratio p-value vs the same ablated-null model.

``p_goodness`` remains available as a legacy helper function, but it is no
longer used for reliability assignment or default output tables.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import chi2

from finaleme_too.core.deconvolution import UNKNOWN_PROFILE
from finaleme_too.core.observation_model import ObservationModel
from finaleme_too.utils.beta_binomial import log_likelihood_per_marker


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
    from finaleme_too.core.observation_model_binarization import (
        _binarize_reference_panel,
    )

    R = np.asarray(reference_methylation, dtype=np.float64)
    M_original, K = R.shape
    if cell_type_index >= K:
        return float("nan")

    # Use the SAME reference-binarization rule as the observation model.
    # In hard-threshold mode this avoids a solver-vs-goodness mismatch
    # (solver used threshold t, goodness used 0.2/0.8).
    hard_thr = getattr(observation, "hard_threshold", None)
    R_binary = _binarize_reference_panel(
        R,
        hard_threshold=hard_thr,
    )
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
        # Optional floor so hard-threshold mode can tolerate a controlled
        # amount of mismatches in p_goodness (without changing the solver).
        tol = 0.0
        metadata = getattr(binarizer, "training_metadata", None)
        if isinstance(metadata, dict):
            tol_raw = metadata.get("p_goodness_mismatch_tolerance")
            if tol_raw is not None:
                try:
                    tol = float(tol_raw)
                except (TypeError, ValueError):
                    tol = 0.0
        tol = float(np.clip(tol, 0.0, 0.49))
        mean_eps = max(mean_eps, tol)
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


def _ablate_component_and_renormalize(w_hat: np.ndarray, cell_type_index: int) -> np.ndarray:
    """Return a null mixture with one cell type removed (w_j=0)."""
    w = np.asarray(w_hat, dtype=np.float64).copy()
    if cell_type_index < 0 or cell_type_index >= w.size:
        return w
    w[cell_type_index] = 0.0
    s = float(np.sum(w))
    if s <= 0.0:
        out = np.zeros_like(w)
        out[-1] = 1.0
        return out
    return w / s


def _top_discriminative_indices_betabinom(
    reference_methylation: np.ndarray,
    observation: ObservationModel,
    cell_type_index: int,
    top_n: int,
) -> np.ndarray:
    R = np.asarray(reference_methylation, dtype=np.float64)
    K = R.shape[1]
    if cell_type_index >= K:
        return np.zeros(0, dtype=np.int64)
    target = R[:, cell_type_index]
    others = np.delete(R, cell_type_index, axis=1)
    bg_mean = np.mean(others, axis=1)
    score = np.abs(target - bg_mean)
    valid = (observation.n > 0) & np.isfinite(score)
    n_top = min(top_n, int(np.sum(valid)))
    if n_top < 5:
        return np.zeros(0, dtype=np.int64)
    score_valid = np.where(valid, score, -np.inf)
    return np.argpartition(-score_valid, n_top - 1)[:n_top].astype(np.int64)


def _top_discriminative_indices_binarization(
    reference_methylation: np.ndarray,
    observation,
    cell_type_index: int,
    top_n: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    from finaleme_too.core.observation_model_binarization import (
        _binarize_reference_panel,
    )

    R = np.asarray(reference_methylation, dtype=np.float64)
    M_original, K = R.shape
    if cell_type_index >= K:
        return (
            np.zeros(0, dtype=np.int64),
            np.zeros(0, dtype=np.int64),
            np.zeros((M_original, K), dtype=np.float64),
        )
    hard_thr = getattr(observation, "hard_threshold", None)
    R_binary = _binarize_reference_panel(R, hard_threshold=hard_thr)
    target = R_binary[:, cell_type_index]
    others = np.delete(R_binary, cell_type_index, axis=1)
    bg_mean = np.mean(others, axis=1)
    score = np.abs(target - bg_mean)
    valid_mask = observation.valid_mask
    candidates = valid_mask & np.isfinite(score)
    n_top = min(top_n, int(np.sum(candidates)))
    if n_top < 5:
        return (
            np.zeros(0, dtype=np.int64),
            np.zeros(0, dtype=np.int64),
            R_binary,
        )
    score_eligible = np.where(candidates, score, -np.inf)
    top_idx_original = np.argpartition(-score_eligible, n_top - 1)[:n_top].astype(
        np.int64
    )
    valid_to_filtered = np.cumsum(valid_mask) - 1
    filtered_idx = valid_to_filtered[top_idx_original].astype(np.int64)
    return top_idx_original, filtered_idx, R_binary


def _binarization_mismatch_tolerance(binarizer) -> float:
    tol = 0.0
    if binarizer is not None:
        metadata = getattr(binarizer, "training_metadata", None)
        if isinstance(metadata, dict):
            raw = metadata.get("p_goodness_mismatch_tolerance")
            if raw is not None:
                try:
                    tol = float(raw)
                except (TypeError, ValueError):
                    tol = 0.0
    return float(np.clip(tol, 0.0, 0.49))


def _likelihood_ratio_pvalue(
    ll_full: float,
    ll_null: float,
    df: int = 1,
) -> float:
    """Likelihood-ratio p-value for nested full vs ablated-null models.

    Statistic:
        LR = 2 * max(0, ll_full - ll_null)
    p-value:
        P(ChiSq_df >= LR)
    """
    if not np.isfinite(ll_full) or not np.isfinite(ll_null):
        return float("nan")
    lr = float(max(0.0, 2.0 * (ll_full - ll_null)))
    return float(chi2.sf(lr, max(int(df), 1)))


def compute_fit_metrics(
    w_hat: np.ndarray,
    reference_methylation: np.ndarray,
    observation,
    cell_type_index: int,
    top_n: int = 50,
    binarizer=None,
) -> tuple[float, float]:
    """Return (likelihood_score, p_likelihood) for one cell type.

    likelihood_score:
      Weighted mean log-likelihood gain (nats per effective marker) over the
      same ablated-null model. Positive means better fit than null.

    p_likelihood:
      Likelihood-ratio p-value for full vs ablated-null fit.
      Lower = stronger evidence the cell type improves fit.
    """
    from finaleme_too.core.observation_model_binarization import (
        BinarizationObservationModel,
    )
    _ = binarizer  # retained for API compatibility

    if isinstance(observation, BinarizationObservationModel):
        _top_idx_original, filtered_idx, _ = _top_discriminative_indices_binarization(
            reference_methylation=reference_methylation,
            observation=observation,
            cell_type_index=cell_type_index,
            top_n=top_n,
        )
        if filtered_idx.size < 5:
            return float("nan"), float("nan")

        w_full = np.asarray(w_hat, dtype=np.float64)
        w_null = _ablate_component_and_renormalize(w_full, cell_type_index)

        coef = np.asarray(observation.coef, dtype=np.float64)[filtered_idx]
        weights = np.asarray(observation.weights, dtype=np.float64)[filtered_idx]
        wsum = float(np.sum(weights))
        if wsum <= 0:
            wsum = float(weights.size)
            weights = np.ones_like(weights, dtype=np.float64)
        p_full = np.clip(coef @ w_full, 1e-15, 1.0)
        p_null = np.clip(coef @ w_null, 1e-15, 1.0)
        ll_full = float(np.sum(weights * np.log(p_full)))
        ll_null = float(np.sum(weights * np.log(p_null)))
        likelihood_score = (ll_full - ll_null) / wsum
        p_likelihood = _likelihood_ratio_pvalue(ll_full, ll_null, df=1)
        return float(likelihood_score), float(p_likelihood)

    # WGBS / beta-binomial path
    top_idx = _top_discriminative_indices_betabinom(
        reference_methylation=reference_methylation,
        observation=observation,
        cell_type_index=cell_type_index,
        top_n=top_n,
    )
    if top_idx.size < 5:
        return float("nan"), float("nan")

    R = np.asarray(reference_methylation, dtype=np.float64)
    R_full = np.hstack([R, np.full((R.shape[0], 1), UNKNOWN_PROFILE, dtype=np.float64)])
    w_full = np.asarray(w_hat, dtype=np.float64)
    w_null = _ablate_component_and_renormalize(w_full, cell_type_index)

    k = np.asarray(observation.k, dtype=np.float64)[top_idx]
    n = np.asarray(observation.n, dtype=np.float64)[top_idx]
    phi = np.asarray(observation.dispersion, dtype=np.float64)[top_idx]
    weights = np.asarray(observation.weights, dtype=np.float64)[top_idx]
    wsum = float(np.sum(weights))
    if wsum <= 0:
        wsum = float(weights.size)
        weights = np.ones_like(weights, dtype=np.float64)

    mu_full = np.clip((R_full @ w_full)[top_idx], 1e-9, 1.0 - 1e-9)
    mu_null = np.clip((R_full @ w_null)[top_idx], 1e-9, 1.0 - 1e-9)
    with np.errstate(invalid="ignore", divide="ignore"):
        mu_obs = np.where(n > 0, k / np.maximum(n, 1), np.nan)
    valid = np.isfinite(mu_obs)
    if int(np.sum(valid)) < 5:
        return float("nan"), float("nan")
    mu_obs_v = mu_obs[valid]
    mu_full_v = mu_full[valid]
    mu_null_v = mu_null[valid]
    k_v = k[valid]
    n_v = n[valid]
    phi_v = phi[valid]
    w_v = weights[valid]
    wsum_v = float(np.sum(w_v))
    if wsum_v <= 0:
        wsum_v = float(w_v.size)
        w_v = np.ones_like(w_v, dtype=np.float64)

    ll_full_vec = log_likelihood_per_marker(k_v, n_v, mu_full_v, phi_v)
    ll_null_vec = log_likelihood_per_marker(k_v, n_v, mu_null_v, phi_v)
    ll_full = float(np.sum(w_v * ll_full_vec))
    ll_null = float(np.sum(w_v * ll_null_vec))
    likelihood_score = (ll_full - ll_null) / wsum_v
    p_likelihood = _likelihood_ratio_pvalue(ll_full, ll_null, df=1)
    return float(likelihood_score), float(p_likelihood)


def compute_unknown_fit_metrics(
    w_hat: np.ndarray,
    reference_methylation: np.ndarray,
    observation,
) -> tuple[float, float]:
    """Return (likelihood_score, p_likelihood) for the Unknown component.

    The null model removes the Unknown weight and renormalizes known
    cell types. Metrics are computed over all valid markers.
    """
    from finaleme_too.core.observation_model_binarization import (
        BinarizationObservationModel,
    )

    w_full = np.asarray(w_hat, dtype=np.float64)
    unknown_idx = int(w_full.size - 1)
    w_null = _ablate_component_and_renormalize(w_full, unknown_idx)

    if isinstance(observation, BinarizationObservationModel):
        coef = np.asarray(observation.coef, dtype=np.float64)
        weights = np.asarray(observation.weights, dtype=np.float64)
        if coef.shape[0] < 5:
            return float("nan"), float("nan")
        wsum = float(np.sum(weights))
        if wsum <= 0:
            wsum = float(weights.size)
            weights = np.ones_like(weights, dtype=np.float64)
        p_full = np.clip(coef @ w_full, 1e-15, 1.0)
        p_null = np.clip(coef @ w_null, 1e-15, 1.0)
        ll_full = float(np.sum(weights * np.log(p_full)))
        ll_null = float(np.sum(weights * np.log(p_null)))
        likelihood_score = float(np.sum(weights * (np.log(p_full) - np.log(p_null))) / wsum)
        p_likelihood = _likelihood_ratio_pvalue(ll_full, ll_null, df=1)
        return likelihood_score, p_likelihood

    # WGBS / beta-binomial path
    R = np.asarray(reference_methylation, dtype=np.float64)
    R_full = np.hstack([R, np.full((R.shape[0], 1), UNKNOWN_PROFILE, dtype=np.float64)])
    n = np.asarray(observation.n, dtype=np.float64)
    k = np.asarray(observation.k, dtype=np.float64)
    phi = np.asarray(observation.dispersion, dtype=np.float64)
    weights = np.asarray(observation.weights, dtype=np.float64)
    valid = (n > 0) & np.all(np.isfinite(R_full), axis=1)
    if int(np.sum(valid)) < 5:
        return float("nan"), float("nan")

    n_v = n[valid]
    k_v = k[valid]
    phi_v = phi[valid]
    weights_v = weights[valid]
    mu_full = np.clip((R_full @ w_full)[valid], 1e-9, 1.0 - 1e-9)
    mu_null = np.clip((R_full @ w_null)[valid], 1e-9, 1.0 - 1e-9)
    with np.errstate(invalid="ignore", divide="ignore"):
        mu_obs = np.where(n_v > 0, k_v / np.maximum(n_v, 1), np.nan)
    finite = np.isfinite(mu_obs)
    if int(np.sum(finite)) < 5:
        return float("nan"), float("nan")
    mu_obs_v = mu_obs[finite]
    mu_full_v = mu_full[finite]
    mu_null_v = mu_null[finite]
    k_f = k_v[finite]
    n_f = n_v[finite]
    phi_f = phi_v[finite]
    w_f = weights_v[finite]
    wsum = float(np.sum(w_f))
    if wsum <= 0:
        wsum = float(w_f.size)
        w_f = np.ones_like(w_f, dtype=np.float64)

    ll_full = np.sum(w_f * log_likelihood_per_marker(k_f, n_f, mu_full_v, phi_f))
    ll_null = np.sum(w_f * log_likelihood_per_marker(k_f, n_f, mu_null_v, phi_f))
    likelihood_score = float((ll_full - ll_null) / wsum)
    p_likelihood = _likelihood_ratio_pvalue(float(ll_full), float(ll_null), df=1)
    return float(likelihood_score), float(p_likelihood)


def assign_reliability(
    p_detection: float,
    likelihood_score: float | None = None,
    p_likelihood: float | None = None,
    effect_size: float | None = None,
) -> str:
    """Assign HIGH/MODERATE/LOW/UNRELIABLE from stability + likelihood metrics.

    Thresholds are intentionally conservative:
      HIGH:      p_detection >= 0.95, likelihood > 0, p_likelihood <= 1e-3
      MODERATE:  p_detection >= 0.50, likelihood > 0, p_likelihood <= 0.05
      UNRELIABLE (when fit metrics unavailable): p_detection < 0.10
      UNRELIABLE (when fit metrics available):   p_detection < 0.50 and
                                                 (likelihood <= 0 or p_likelihood > 0.20)
      else LOW.
    """
    _ = effect_size  # deprecated, ignored
    if np.isnan(p_detection):
        p_detection = 0.0

    lik_missing = likelihood_score is None or np.isnan(likelihood_score)
    pl_missing = p_likelihood is None or np.isnan(p_likelihood)
    if lik_missing or pl_missing:
        if p_detection >= 0.95:
            return "HIGH"
        if p_detection >= 0.50:
            return "MODERATE"
        if p_detection < 0.10:
            return "UNRELIABLE"
        return "LOW"

    lik = float(likelihood_score)
    pl = float(p_likelihood)

    if p_detection >= 0.95 and lik > 0.0 and pl <= 1e-3:
        return "HIGH"
    if p_detection >= 0.50 and lik > 0.0 and pl <= 0.05:
        return "MODERATE"
    if p_detection < 0.50 and (lik <= 0.0 or pl > 0.20):
        return "UNRELIABLE"
    return "LOW"


__all__ = [
    "assign_reliability",
    "compute_fit_metrics",
    "compute_unknown_fit_metrics",
    "compute_p_detection",
    "compute_p_goodness",
]
