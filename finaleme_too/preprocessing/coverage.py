"""Coverage tier assignment and effective coverage computation (architecture §4)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from finaleme_too.config import CoverageConfig, CoverageTier

if TYPE_CHECKING:
    from finaleme_too.io.methylation_loader import MarkerObservations


# Median cfDNA fragment length used to convert total read counts at marker
# positions into a depth-of-coverage estimate over the marker region area.
FRAGMENT_LENGTH_BP = 167
# Whole human-genome size in bp. Only used as a sanity ceiling for the
# marker-region area (a marker bed should never exceed the genome size).
HG_GENOME_BP = 3_000_000_000


def effective_coverage_in_markers(obs: "MarkerObservations") -> float:
    """Mean depth-of-coverage **within the marker regions**.

    The numerator ``sum(obs.n)`` is the total read count observed at all
    CpG sites that fall within marker regions (these come from
    ``total_count_predict`` in FinaleMe prediction.bed.gz, or column 8 in
    Bis-SNP 6+2 BED, etc.). The denominator is the total marker-region
    area in base pairs, computed directly from the marker BED file
    (``end - start`` summed across all markers). Multiplying by the
    median cfDNA fragment length converts the read count into a depth-of-
    coverage figure that is comparable to whole-genome depth values:

        coverage = Σ obs.n * fragment_length / Σ (obs.end - obs.start)

    Returns 0.0 when the marker area is empty (defensive — should never
    happen in practice because the loader populates start/end from the
    marker BED).
    """
    n = np.asarray(obs.n, dtype=np.int64)
    if n.size == 0:
        return 0.0
    starts = np.asarray(obs.start, dtype=np.int64)
    ends = np.asarray(obs.end, dtype=np.int64)
    width = int(np.sum(np.maximum(ends - starts, 0)))
    if width <= 0:
        return 0.0
    total_reads = int(np.sum(n))
    return float(total_reads) * FRAGMENT_LENGTH_BP / float(width)


@dataclass
class CoverageTierAssigner:
    """Assign per-sample coverage tier (HIGH / LOW / ULTRALOW).

    Tier classification thresholds the **effective coverage in the marker
    regions** (architecture §4.2):

        coverage = Σ reads_at_markers * fragment_length / Σ marker_widths

    See ``effective_coverage_in_markers`` for the formula. This is the
    same number reported as ``mean_coverage`` in ``qc_summary.tsv``, so
    ``coverage_tier`` and ``mean_coverage`` are always on the same scale
    and self-consistent.

    The previous version of this routine implemented the architecture
    doc's whole-genome formula ``total_fragments * fragment_length /
    genome_size``, but substituted ``np.sum(obs.n)`` (reads at marker
    positions only) for ``total_fragments`` (whole-BAM total reads).
    That always under-counted dramatically — a sample with mean 100
    reads-per-marker over 1000 markers was misclassified as ULTRALOW
    (mean_cov ≈ 0.006). The marker-area denominator is the right
    denominator: it matches the numerator's scope and needs no extra
    sample-sheet inputs.

    Thresholds:
      * ``tier_high`` (default 10) — depth above this is HIGH
      * ``tier_low`` (default 0.5) — depth above this (but ≤ tier_high) is LOW
      * Below ``tier_low`` → ULTRALOW
    """

    config: CoverageConfig

    def assign(self, obs: "MarkerObservations") -> CoverageTier:
        n = np.asarray(obs.n, dtype=np.int64)
        if n.size == 0 or int(np.sum(n)) == 0:
            return CoverageTier.ULTRALOW

        coverage = effective_coverage_in_markers(obs)
        if coverage > self.config.tier_high:
            return CoverageTier.HIGH
        if coverage >= self.config.tier_low:
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
    "effective_coverage_in_markers",
    "per_marker_effective_coverage",
    "per_marker_min_reads",
    "per_marker_min_reads_vector",
]
