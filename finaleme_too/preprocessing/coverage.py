"""Coverage tier assignment and effective coverage computation (architecture §4)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from finaleme_too.config import CoverageConfig, CoverageTier

if TYPE_CHECKING:
    from finaleme_too.io.methylation_loader import MarkerObservations


# Average human cfDNA fragment length and a typical genome size used to convert
# total reads to a rough genome-wide depth. These match the architecture doc
# defaults; users can override via config when needed.
FRAGMENT_LENGTH_BP = 167  # median cfDNA fragment length
HG_GENOME_BP = 3_000_000_000


@dataclass
class CoverageTierAssigner:
    """Assign per-sample coverage tier (HIGH / LOW / ULTRALOW)."""

    config: CoverageConfig

    def assign(self, obs: "MarkerObservations") -> CoverageTier:
        # Use total reads × fragment length / genome size as the genome-wide
        # mean coverage estimate (architecture §4.2).
        total_reads = int(np.sum(obs.n))
        if total_reads == 0:
            return CoverageTier.ULTRALOW
        mean_cov = total_reads * FRAGMENT_LENGTH_BP / HG_GENOME_BP
        if mean_cov > self.config.tier_high:
            return CoverageTier.HIGH
        if mean_cov >= self.config.tier_low:
            return CoverageTier.LOW
        return CoverageTier.ULTRALOW


def per_marker_effective_coverage(
    obs: "MarkerObservations", expected_per_marker: float | None = None
) -> np.ndarray:
    """Effective coverage per marker (architecture §4.2).

    effective(i) = observed_fragments(i) / expected_fragments(i)

    If ``expected_per_marker`` is None we use the genome-wide mean as a proxy
    so that the values are roughly centered on 1.0.
    """
    n = np.asarray(obs.n, dtype=np.float64)
    if expected_per_marker is None:
        if n.size == 0 or float(np.mean(n)) == 0:
            return np.zeros_like(n)
        expected_per_marker = float(np.mean(n))
    return n / max(expected_per_marker, 1.0)


def per_marker_min_reads(tier: CoverageTier) -> int:
    """Per-tier minimum read threshold (architecture §4.3)."""
    if tier == CoverageTier.HIGH:
        return 3
    if tier == CoverageTier.LOW:
        return 2
    return 1


# Thresholds (on the relative effective-coverage ratio) used to down-tier
# individual markers whose local coverage is below the sample-wide expected
# value (architecture §4.2). A marker with ratio >= 1.0 stays at the sample's
# global tier; ratio in [0.1, 1.0) drops by one tier; ratio < 0.1 drops by two.
_EFFECTIVE_DROP_ONE = 1.0
_EFFECTIVE_DROP_TWO = 0.1


_TIER_ORDER = [CoverageTier.HIGH, CoverageTier.LOW, CoverageTier.ULTRALOW]


def per_marker_min_reads_vector(
    obs_n: np.ndarray,
    global_tier: CoverageTier,
) -> np.ndarray:
    """Per-marker minimum read threshold with effective-coverage down-tiering.

    Implements architecture §4.2: a marker with effective coverage below a
    tier-specific threshold is treated as belonging to the next lower tier
    for that sample.

    The effective coverage ratio is ``n_i / mean_n`` where ``mean_n`` is the
    sample-wide mean reads-per-marker. Markers with ratio >= 1.0 keep the
    sample's global tier; ratio in [0.1, 1.0) drops by one tier (HIGH→LOW,
    LOW→ULTRALOW); ratio < 0.1 drops by two tiers. The down-tiered tier
    determines the per-marker min-reads threshold via ``per_marker_min_reads``.

    The *less strict* direction is intentional: down-tiering lets a marker
    with below-average coverage still pass the filter at a lower
    reliability bar, rather than being thrown out entirely.
    """
    n = np.asarray(obs_n, dtype=np.float64)
    if n.size == 0:
        return np.zeros(0, dtype=np.int64)
    mean_n = float(np.mean(n))
    if mean_n <= 0:
        return np.full(n.size, per_marker_min_reads(global_tier), dtype=np.int64)
    effective = n / mean_n

    try:
        global_idx = _TIER_ORDER.index(global_tier)
    except ValueError:
        global_idx = 0  # default to HIGH

    def _tier_for(eff: float) -> CoverageTier:
        if eff >= _EFFECTIVE_DROP_ONE:
            return _TIER_ORDER[global_idx]
        if eff >= _EFFECTIVE_DROP_TWO:
            return _TIER_ORDER[min(global_idx + 1, 2)]
        return _TIER_ORDER[min(global_idx + 2, 2)]

    out = np.empty(n.size, dtype=np.int64)
    for i in range(n.size):
        out[i] = per_marker_min_reads(_tier_for(float(effective[i])))
    return out


__all__ = [
    "CoverageTierAssigner",
    "FRAGMENT_LENGTH_BP",
    "HG_GENOME_BP",
    "per_marker_effective_coverage",
    "per_marker_min_reads",
    "per_marker_min_reads_vector",
]
