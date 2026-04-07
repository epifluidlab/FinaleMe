"""TOOPipeline orchestrator.

Phase A scope: single-sample MLE deconvolution end-to-end with bootstrap CIs
and per-cell-type reliability p-values. Coverage tiers, calibration, batch
correction, ILR testing, and Bayesian/fragment paths are added in later phases.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from finaleme_too.config import (
    CoverageTier,
    MeasurementMode,
    TOOConfig,
)
from finaleme_too.core.deconvolution import (
    DeconvolutionResult,
    MLEDeconvolver,
)
from finaleme_too.core.observation_model import BetaBinomialModel, ObservationModel
from finaleme_too.core.reliability import (
    assign_reliability,
    compute_p_detection,
    compute_p_goodness,
)
from finaleme_too.core.uncertainty import BootstrapCI
from finaleme_too.io.marker_regions import MarkerRegions, MarkerRegionsLoader
from finaleme_too.io.methylation_loader import MarkerObservations, MethylationLoader
from finaleme_too.io.output_writer import (
    write_cohort_proportions,
    write_per_sample_too,
    write_qc_summary,
)
from finaleme_too.io.reference_panel import (
    ReferencePanel,
    ReferencePanelLoader,
    load_cpg_index,
)
from finaleme_too.io.sample_sheet import Sample, SampleSheet
from finaleme_too.utils.parallel import parallel_map

log = logging.getLogger(__name__)


@dataclass
class CohortResult:
    samples: list[DeconvolutionResult]
    sample_groups: dict[str, str | None]


class TOOPipeline:
    """End-to-end TOO pipeline (Phase A scope: P0)."""

    def __init__(self, config: TOOConfig):
        self.config = config
        self.deconvolver = MLEDeconvolver()
        self.bootstrap = BootstrapCI(
            n_iterations=config.uncertainty.n_bootstrap,
            ci_level=config.uncertainty.ci_level,
        )
        self.observation_builder = BetaBinomialModel()

    # ------------------------------------------------------------------
    # Top-level
    # ------------------------------------------------------------------

    def run(
        self,
        sample_sheet: SampleSheet,
        reference: ReferencePanel,
        marker_regions: MarkerRegions,
        output_dir: str | Path,
        cpg_index: dict | None = None,
    ) -> CohortResult:
        """Run the full pipeline and write outputs to ``output_dir``."""
        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        # Per-sample work in parallel
        def _process(sample: Sample) -> DeconvolutionResult:
            return self._process_sample(sample, reference, marker_regions, cpg_index, out_dir)

        results = parallel_map(_process, sample_sheet.samples, n_jobs=self.config.threads)

        # Cohort outputs
        sample_groups = {s.sample_id: s.group for s in sample_sheet.samples}
        write_cohort_proportions(results, sample_groups, out_dir / "cohort_proportions.tsv")
        write_qc_summary(results, sample_groups, out_dir / "qc_summary.tsv")
        return CohortResult(samples=results, sample_groups=sample_groups)

    # ------------------------------------------------------------------
    # Per-sample
    # ------------------------------------------------------------------

    def _process_sample(
        self,
        sample: Sample,
        reference: ReferencePanel,
        marker_regions: MarkerRegions,
        cpg_index: dict | None,
        out_dir: Path,
    ) -> DeconvolutionResult:
        log.info("Processing sample %s (%s)", sample.sample_id, sample.mode.value)
        obs = MethylationLoader.load(
            filepath=sample.methylation_file,
            sample_id=sample.sample_id,
            mode=sample.mode,
            marker_regions=marker_regions,
            input_format=sample.input_format,
            meth_col=sample.meth_col,
            total_col=sample.total_col,
            cpg_index=cpg_index,
        )

        tier = self._assign_tier_simple(obs)
        observation = self.observation_builder.build(
            obs=obs,
            reference=reference,
            tier=tier,
            coverage_cap=self.config.coverage.coverage_cap,
        )

        # Point-estimate proportions
        w_hat = self.deconvolver.solve(observation, reference)

        # Bootstrap CIs
        boot = self.bootstrap.estimate(observation, reference, self.deconvolver)

        # Reliability p-values per cell type
        K = reference.n_cell_types
        p_goodness = np.full(K, np.nan, dtype=np.float64)
        p_detection = np.zeros(K + 1, dtype=np.float64)
        for j in range(K):
            p_goodness[j] = compute_p_goodness(
                w_hat=w_hat,
                reference_methylation=reference.methylation,
                observation=observation,
                cell_type_index=j,
            )
            p_detection[j] = compute_p_detection(
                boot.proportions_samples[:, j],
                noise_floor=self.config.uncertainty.noise_floor,
            )
        # Unknown component detection (no goodness)
        p_detection[K] = compute_p_detection(
            boot.proportions_samples[:, K],
            noise_floor=self.config.uncertainty.noise_floor,
        )

        reliability = np.empty(K + 1, dtype=object)
        for j in range(K):
            reliability[j] = assign_reliability(p_goodness[j], p_detection[j])
        reliability[K] = assign_reliability(np.nan, p_detection[K])

        # n_markers per cell type — count valid markers contributing to that cell type
        valid_n = np.sum(observation.n > 0)
        n_markers = np.full(K, int(valid_n), dtype=np.int32)

        result = DeconvolutionResult(
            sample_id=sample.sample_id,
            cell_types=list(reference.cell_types),
            proportions=w_hat,
            ci_lower=boot.ci_lower,
            ci_upper=boot.ci_upper,
            p_goodness=p_goodness,
            p_detection=p_detection,
            reliability=reliability,
            n_markers=n_markers,
            bootstrap_proportions=boot.proportions_samples,
            posterior_samples=None,
            coverage_tier=tier,
            qc_flags=self._compute_qc_flags(w_hat, observation, tier),
        )

        write_per_sample_too(result, out_dir / f"{sample.sample_id}.too.tsv")
        return result

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _assign_tier_simple(self, obs: MarkerObservations) -> CoverageTier:
        """Phase A: simple total-reads-per-marker heuristic.

        The full effective-coverage logic from §4 is added in Phase B.
        """
        if obs.n.size == 0:
            return CoverageTier.ULTRALOW
        median_cov = float(np.median(obs.n))
        if median_cov >= self.config.coverage.tier_high:
            return CoverageTier.HIGH
        if median_cov >= self.config.coverage.tier_low:
            return CoverageTier.LOW
        return CoverageTier.ULTRALOW

    def _compute_qc_flags(
        self,
        w_hat: np.ndarray,
        observation: ObservationModel,
        tier: CoverageTier,
    ) -> list[str]:
        flags: list[str] = []
        unknown = float(w_hat[-1])
        if unknown > self.config.qc.max_unknown_fraction:
            flags.append("HIGH_UNKNOWN")
        if tier == CoverageTier.ULTRALOW:
            flags.append("ULTRALOW_COVERAGE")
        return flags


# ---------------------------------------------------------------------------
# Convenience builder used by the CLI
# ---------------------------------------------------------------------------


def build_reference_and_markers(
    config: TOOConfig,
    explicit_reference: str | None = None,
    explicit_markers: str | None = None,
    explicit_ref_betas: str | None = None,
    explicit_ref_groups: str | None = None,
    explicit_cpg_index: str | None = None,
    explicit_marker_format: str | None = None,
) -> tuple[ReferencePanel, MarkerRegions, dict | None]:
    """Resolve reference panel and marker regions from CLI / config."""
    cpg_index = None
    cpg_index_path = explicit_cpg_index or config.reference.cpg_index
    if cpg_index_path is not None:
        cpg_index = load_cpg_index(cpg_index_path)

    # Marker regions
    marker_regions_path = explicit_markers or config.markers.marker_regions
    marker_format = explicit_marker_format or config.markers.marker_format

    marker_regions: MarkerRegions
    if marker_regions_path:
        marker_regions = MarkerRegionsLoader.load(marker_regions_path, marker_format=marker_format)
    elif explicit_reference or config.reference.reference_panel:
        # Coordinates can come from the matrix reference panel itself
        path = explicit_reference or config.reference.reference_panel
        ref_for_coords = ReferencePanelLoader.load_matrix(path)
        marker_regions = ref_for_coords.to_marker_regions()
    else:
        raise ValueError(
            "No marker regions: provide --marker-regions or --reference-panel"
        )

    # Reference panel
    if explicit_reference or config.reference.reference_panel:
        reference = ReferencePanelLoader.load_matrix(
            explicit_reference or config.reference.reference_panel
        )
    elif explicit_ref_betas or config.reference.ref_betas:
        ref_betas_arg = explicit_ref_betas or config.reference.ref_betas
        ref_groups = explicit_ref_groups or config.reference.ref_groups
        if cpg_index_path is None:
            raise ValueError(
                "--ref-betas requires --cpg-index for binary .beta parsing"
            )
        if ref_groups is None:
            raise ValueError("--ref-betas requires --ref-groups")
        reference = ReferencePanelLoader.load_beta_list(
            ref_betas_arg=ref_betas_arg,
            groups_file=ref_groups,
            cpg_index_path=cpg_index_path,
            marker_regions=marker_regions,
        )
    else:
        raise ValueError("No reference panel: provide --reference-panel or --ref-betas")

    return reference, marker_regions, cpg_index


__all__ = ["CohortResult", "TOOPipeline", "build_reference_and_markers"]
