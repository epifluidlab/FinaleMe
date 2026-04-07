"""ILR-space covariate residualization (architecture §9.1, math doc §7).

Applied AFTER deconvolution. Builds an OLS model in ILR coordinates and
returns the residuals projected back to the simplex. The user supplies the
list of biological covariates to adjust for; configurable covariates
(treatment, treatment_efficacy, mutation_status) can be opt-in.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from finaleme_too.utils.transforms import ilr_inverse, ilr_transform


def adjust_covariates(
    proportions: np.ndarray,
    sample_ids: list[str],
    covariates: pd.DataFrame,
    columns: list[str],
) -> np.ndarray:
    """Residualize ILR-transformed proportions on the chosen covariates.

    Parameters
    ----------
    proportions : (S, K+1)
        Per-sample proportions including the unknown component.
    sample_ids : list of str
        Used to align proportions with covariate rows by index.
    covariates : DataFrame
        Indexed (or column-keyed) by sample_id; must contain ``columns``.
    columns : list of str
        Covariate columns to adjust for. Skipped if missing or all-NaN.

    Returns
    -------
    Adjusted proportions, same shape as input.
    """
    if not columns or covariates is None or covariates.empty:
        return proportions
    df = covariates.copy()
    if "sample_id" in df.columns:
        df = df.set_index("sample_id")
    keep_cols = [c for c in columns if c in df.columns]
    if not keep_cols:
        return proportions
    aligned = df.reindex(sample_ids)[keep_cols]

    # Build numeric design matrix (one-hot encode categoricals)
    design = pd.get_dummies(aligned, drop_first=True, dummy_na=False)
    design = design.dropna(how="any")
    if design.shape[0] < 3 or design.shape[1] == 0:
        return proportions
    X_design = design.to_numpy(dtype=np.float64)
    X_design = np.column_stack([np.ones(X_design.shape[0]), X_design])  # intercept

    valid_idx = [sample_ids.index(s) for s in design.index.tolist()]
    Y = ilr_transform(proportions[valid_idx])

    coeffs, *_ = np.linalg.lstsq(X_design, Y, rcond=None)
    Y_pred = X_design @ coeffs
    Y_resid = Y - Y_pred + Y.mean(axis=0)
    adjusted_subset = ilr_inverse(Y_resid, n_parts=proportions.shape[1])

    # Write back into a copy of the original proportions
    out = proportions.copy()
    for local_i, sid in enumerate(design.index.tolist()):
        s_idx = sample_ids.index(sid)
        out[s_idx] = adjusted_subset[local_i]
    return out


__all__ = ["adjust_covariates"]
