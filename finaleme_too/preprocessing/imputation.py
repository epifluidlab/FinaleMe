"""Same-group cohort imputation for low-coverage markers (math doc §9.1, §9.2)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from finaleme_too.exceptions import IllegalImputationError
from finaleme_too.io.methylation_loader import MarkerObservations


@dataclass
class CohortImputer:
    """Impute missing-coverage markers from same-group cohort samples.

    Math doc §9:
        μ̂_{i,s} = Σ_{s' ∈ G(s)} v_{s'} * μ_{i,s'} / Σ_{s'} v_{s'}
        v_{s'} = n_{i,s'} * I(n_{i,s'} ≥ threshold)

    Constraints (§9.2):
        - Never impute across comparison groups.
        - Require **≥ min_donors eligible donors *per marker*** (not globally).
          Markers without enough eligible donors are left unchanged.
    """

    coverage_threshold: int = 3
    min_donors: int = 3

    def impute(
        self,
        sample: MarkerObservations,
        cohort: list[MarkerObservations],
        sample_groups: dict[str, str | None],
        strict: bool = False,
    ) -> MarkerObservations:
        """Impute missing-coverage markers in ``sample`` from same-group cohort.

        When ``strict`` is True (library use), missing group labels raise
        ``IllegalImputationError``. When False (pipeline use, the default),
        the sample is returned unchanged so LOW/ULTRALOW runs do not crash
        on unlabeled samples. The sample sheet parser enforces non-empty
        groups at the entry point, so this is defense-in-depth.
        """
        sample_group = sample_groups.get(sample.sample_id)
        if sample_group is None:
            if strict:
                raise IllegalImputationError(
                    f"Sample {sample.sample_id} has no group label; cannot impute"
                )
            return sample

        donors = [
            obs
            for obs in cohort
            if obs.sample_id != sample.sample_id
            and sample_groups.get(obs.sample_id) == sample_group
        ]
        if not donors:
            return sample

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
        # Per-marker eligibility: only donors with enough coverage at THIS marker
        eligible = donor_n >= self.coverage_threshold  # (n_donors, n_markers)
        eligible_count = eligible.sum(axis=0)  # (n_markers,)
        # Per-marker constraint: >= min_donors eligible donors
        passes_donor_count = eligible_count >= self.min_donors

        # Weighted mean of donor methylation, but only over eligible donors
        weights = donor_n * eligible
        weight_sum = weights.sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            beta_hat = np.where(
                weight_sum > 0,
                (donor_k * eligible).sum(axis=0) / np.maximum(weight_sum, 1),
                np.nan,
            )

        # Synthetic n: median over **eligible** donors only (not all)
        # so that medians aren't dragged down by low-coverage donors at
        # this particular marker.
        donor_n_masked = np.where(eligible, donor_n, np.nan)
        with np.errstate(invalid="ignore"):
            median_eligible_n = np.nanmedian(donor_n_masked, axis=0)
        # When no eligible donor exists the median is NaN; fall back to 1.
        median_eligible_n = np.where(
            np.isfinite(median_eligible_n), median_eligible_n, 1.0
        )
        synthetic_n = np.where(
            np.isfinite(beta_hat) & passes_donor_count,
            np.maximum(median_eligible_n, 1).astype(np.int64),
            n,
        )
        synthetic_k = np.where(
            np.isfinite(beta_hat) & passes_donor_count,
            np.round(beta_hat * synthetic_n).astype(np.int64),
            k,
        )

        do_impute = below & np.isfinite(beta_hat) & passes_donor_count
        new_k = np.where(do_impute, synthetic_k, k).astype(np.int32)
        new_n = np.where(do_impute, synthetic_n, n).astype(np.int32)
        return sample.with_counts(new_k, new_n)


__all__ = ["CohortImputer"]
