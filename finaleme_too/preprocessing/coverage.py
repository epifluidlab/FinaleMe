"""Coverage tier assignment and effective coverage computation (architecture §4)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from finaleme_too.config import CoverageConfig, CoverageTier
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

    def assign(self, obs: MarkerObservations) -> CoverageTier:
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
    obs: MarkerObservations, expected_per_marker: float | None = None
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


__all__ = [
    "CoverageTierAssigner",
    "FRAGMENT_LENGTH_BP",
    "HG_GENOME_BP",
    "per_marker_effective_coverage",
    "per_marker_min_reads",
]
