"""Same-group cohort imputation for low-coverage markers (math doc §9.1)."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from finaleme_too.exceptions import IllegalImputationError
from finaleme_too.io.methylation_loader import MarkerObservations


@dataclass
class CohortImputer:
    """Impute missing-coverage markers from same-group cohort samples."""

    coverage_threshold: int = 3
    min_donors: int = 3

    def impute(
        self,
        sample: MarkerObservations,
        cohort: list[MarkerObservations],
        sample_groups: dict[str, str | None],
    ) -> MarkerObservations:
        """Return a new MarkerObservations with imputed counts where needed.

        Per math doc §9.1:
            μ̂_{i,s} = Σ_{s' ∈ G(s)} v_{s'} * μ_{i,s'} / Σ_{s'} v_{s'}
        where v_{s'} = n_{i,s'} * I(n_{i,s'} ≥ threshold).

        ``sample_groups`` maps sample_id → group label. Imputation only uses
        donors in the same group as ``sample``. Raises IllegalImputationError
        if no group label is set on ``sample``.
        """
        sample_group = sample_groups.get(sample.sample_id)
        if sample_group is None:
            raise IllegalImputationError(
                f"Sample {sample.sample_id} has no group label; cannot impute"
            )

        donors = [
            obs
            for obs in cohort
            if obs.sample_id != sample.sample_id
            and sample_groups.get(obs.sample_id) == sample_group
        ]
        if len(donors) < self.min_donors:
            return sample  # not enough donors — leave as-is

        n = np.asarray(sample.n, dtype=np.int64)
        k = np.asarray(sample.k, dtype=np.int64)
        below = n < self.coverage_threshold
        if not np.any(below):
            return sample

        donor_n = np.stack(
            [np.asarray(d.n, dtype=np.float64) for d in donors], axis=0
        )
        donor_k = np.stack(
            [np.asarray(d.k, dtype=np.float64) for d in donors], axis=0
        )
        eligible = donor_n >= self.coverage_threshold
        weights = donor_n * eligible

        weight_sum = weights.sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            beta_hat = np.where(
                weight_sum > 0,
                (donor_k * eligible).sum(axis=0) / np.maximum(weight_sum, 1),
                np.nan,
            )

        # Median donor coverage gives the synthetic n for imputed markers
        median_donor_n = np.median(donor_n, axis=0)
        synthetic_n = np.where(
            np.isfinite(beta_hat), np.maximum(median_donor_n, 1).astype(np.int64), n
        )
        synthetic_k = np.where(
            np.isfinite(beta_hat),
            np.round(beta_hat * synthetic_n).astype(np.int64),
            k,
        )

        new_k = np.where(below & np.isfinite(beta_hat), synthetic_k, k).astype(np.int32)
        new_n = np.where(below & np.isfinite(beta_hat), synthetic_n, n).astype(np.int32)
        return sample.with_counts(new_k, new_n)


__all__ = ["CohortImputer"]
