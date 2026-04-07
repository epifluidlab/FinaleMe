"""TOOPipeline orchestrator.

Phase A scope (P0): single-sample MLE deconvolution end-to-end with bootstrap
CIs and per-cell-type reliability p-values.

Phase B scope (P1): coverage tiers, marker selection, FinaleMe calibration
apply path, ILR statistical testing, multi-group comparisons, QC summary.

Phases C-E (calibration training, batch correction, Bayesian, fragment-level)
extend further.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from finaleme_too.config import (
    CoverageTier,
    MeasurementMode,
    TOOConfig,
)
from finaleme_too.core.deconvolution import (
    BayesianDeconvolver,
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
    write_group_comparison,
    write_per_sample_too,
    write_qc_summary,
)
from finaleme_too.io.reference_panel import (
    ReferencePanel,
    ReferencePanelLoader,
    load_cpg_index,
)
from finaleme_too.io.sample_sheet import Sample, SampleSheet
from finaleme_too.preprocessing.batch_correction import combat_correct
from finaleme_too.preprocessing.calibration import (
    CalibrationParams,
    apply_calibration,
    load_default_calibration,
)
from finaleme_too.preprocessing.calibration_eval import compute_inference_qc
from finaleme_too.preprocessing.coverage import CoverageTierAssigner
from finaleme_too.preprocessing.imputation import CohortImputer
from finaleme_too.preprocessing.marker_selection import BalancedMarkerSelector
from finaleme_too.postprocessing.covariate_adjustment import adjust_covariates
from finaleme_too.postprocessing.group_comparison import run_group_comparisons
from finaleme_too.postprocessing.qc import compute_qc_flags
from finaleme_too.utils.parallel import parallel_map

log = logging.getLogger(__name__)


@dataclass
class CohortResult:
    samples: list[DeconvolutionResult]
    sample_groups: dict[str, str | None]


class TOOPipeline:
    """End-to-end TOO pipeline (Phase A + B scope: P0 + P1)."""

    def __init__(
        self,
        config: TOOConfig,
        calibration: CalibrationParams | None = None,
        region_annotations: pd.DataFrame | None = None,
        group_comparison_spec: str | None = None,
    ):
        from finaleme_too.config import SolverMethod

        self.config = config
        self.calibration = calibration
        self.region_annotations = region_annotations
        self.group_comparison_spec = group_comparison_spec
        self.deconvolver = MLEDeconvolver()
        self.bayesian_deconvolver = (
            BayesianDeconvolver(
                n_walkers=64,
                n_steps=config.uncertainty.bayesian_n_samples // 64 + 100,
                burn_in=100,
                prior_alpha=config.uncertainty.bayesian_prior_alpha,
            )
            if config.model.deconvolution == SolverMethod.BAYESIAN
            else None
        )
        self.bootstrap = BootstrapCI(
            n_iterations=config.uncertainty.n_bootstrap,
            ci_level=config.uncertainty.ci_level,
        )
        self.observation_builder = BetaBinomialModel()
        self.tier_assigner = CoverageTierAssigner(config.coverage)
        self.marker_selector = (
            BalancedMarkerSelector(
                n_per_type=config.markers.n_per_type,
                method=config.markers.specificity_method,
                strict_regions=config.markers.strict_regions,
            )
            if config.markers.n_per_type and config.markers.n_per_type > 0
            else None
        )

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

        # Cell-type-balanced marker selection (P1)
        marker_subset_indices: np.ndarray | None = None
        if self.marker_selector is not None and reference.n_markers > self.marker_selector.n_per_type:
            marker_subset_indices = self.marker_selector.select(
                reference,
                region_annotations=self.region_annotations,
            )
            log.info(
                "Marker selection: kept %d / %d markers",
                marker_subset_indices.size,
                reference.n_markers,
            )
            reference = _subset_reference(reference, marker_subset_indices)
            marker_regions = _subset_marker_regions(marker_regions, marker_subset_indices)

        # Phase D: load all samples first so we can do cohort-level operations
        # (batch correction, imputation) before deconvolution.
        log.info("Loading + calibrating %d samples", len(sample_sheet.samples))
        loaded_obs: list[MarkerObservations] = parallel_map(
            lambda s: self._load_and_calibrate(s, marker_regions, cpg_index),
            sample_sheet.samples,
            n_jobs=self.config.threads,
        )

        # Optional: ComBat-style technical batch correction
        if self.config.batch_correction.technical_covariates:
            batch_var = self.config.batch_correction.technical_covariates[0]
            batch_labels = [
                str(s.metadata.get(batch_var)) if batch_var in s.metadata else None
                for s in sample_sheet.samples
            ]
            log.info("Applying batch correction on %s", batch_var)
            loaded_obs = combat_correct(
                loaded_obs,
                batch_labels=batch_labels,
                min_levels=self.config.batch_correction.min_levels,
                min_per_level=self.config.batch_correction.min_samples_per_level,
            )

        # Optional: cohort imputation for low-coverage markers (LOW/ULTRALOW)
        sample_groups_map = {s.sample_id: s.group for s in sample_sheet.samples}
        if self._should_impute(loaded_obs):
            imputer = CohortImputer(coverage_threshold=self.config.coverage.min_reads)
            loaded_obs = [
                imputer.impute(obs, loaded_obs, sample_groups_map) for obs in loaded_obs
            ]

        # Per-sample deconvolution (parallel)
        def _decon(args: tuple[Sample, MarkerObservations]) -> DeconvolutionResult:
            sample, obs = args
            return self._deconvolve_sample(sample, obs, reference, out_dir)

        results = parallel_map(
            _decon,
            list(zip(sample_sheet.samples, loaded_obs)),
            n_jobs=self.config.threads,
        )

        # Phase D: optional ILR-space biological covariate adjustment
        adjusted_results = self._maybe_adjust_covariates(results, sample_sheet)

        # Cohort outputs
        sample_groups = {s.sample_id: s.group for s in sample_sheet.samples}
        write_cohort_proportions(
            adjusted_results, sample_groups, out_dir / "cohort_proportions.tsv"
        )
        write_qc_summary(adjusted_results, sample_groups, out_dir / "qc_summary.tsv")

        # Group comparison (P1)
        if self.group_comparison_spec and len(results) >= 4:
            test_results = self._run_group_comparisons(adjusted_results, sample_groups)
            if test_results:
                write_group_comparison(test_results, out_dir / "group_comparison.tsv")

        return CohortResult(samples=adjusted_results, sample_groups=sample_groups)

    def _should_impute(self, observations: list[MarkerObservations]) -> bool:
        """Phase D: enable imputation if any sample is in LOW/ULTRALOW tier."""
        for obs in observations:
            tier = self.tier_assigner.assign(obs)
            if tier in (CoverageTier.LOW, CoverageTier.ULTRALOW):
                return True
        return False

    def _maybe_adjust_covariates(
        self,
        results: list[DeconvolutionResult],
        sample_sheet: SampleSheet,
    ) -> list[DeconvolutionResult]:
        """Apply ILR-space biological covariate adjustment if configured."""
        bio = self.config.covariate_adjustment.biological_covariates
        if not bio or len(results) < 3:
            return results

        # Build a (S, K+1) proportions matrix and a covariates DataFrame
        sample_ids = [r.sample_id for r in results]
        prop = np.array([r.proportions for r in results], dtype=np.float64)
        rows = []
        for s in sample_sheet.samples:
            row = {"sample_id": s.sample_id}
            row.update(s.biological_covariates)
            rows.append(row)
        cov_df = pd.DataFrame(rows)
        adjusted = adjust_covariates(
            proportions=prop,
            sample_ids=sample_ids,
            covariates=cov_df,
            columns=bio,
        )
        out: list[DeconvolutionResult] = []
        for r, new_p in zip(results, adjusted):
            out.append(
                DeconvolutionResult(
                    sample_id=r.sample_id,
                    cell_types=r.cell_types,
                    proportions=new_p,
                    ci_lower=r.ci_lower,
                    ci_upper=r.ci_upper,
                    p_goodness=r.p_goodness,
                    p_detection=r.p_detection,
                    reliability=r.reliability,
                    n_markers=r.n_markers,
                    bootstrap_proportions=r.bootstrap_proportions,
                    posterior_samples=r.posterior_samples,
                    coverage_tier=r.coverage_tier,
                    qc_flags=list(r.qc_flags),
                )
            )
        return out

    def _run_group_comparisons(
        self,
        results: list[DeconvolutionResult],
        sample_groups: dict[str, str | None],
    ) -> list:
        K = len(results[0].cell_types)
        prop = np.array([r.proportions for r in results], dtype=np.float64)  # (S, K+1)
        sample_ids = [r.sample_id for r in results]
        labels = [sample_groups.get(sid) for sid in sample_ids]

        # Bayesian comparison needs per-sample posterior draws (Phase E only)
        posterior_by_sample: dict[str, np.ndarray] | None = None
        if (
            self.config.testing.method.value == "bayesian_posterior"
            and any(r.posterior_samples is not None for r in results)
        ):
            posterior_by_sample = {
                r.sample_id: r.posterior_samples
                for r in results
                if r.posterior_samples is not None
            }

        return run_group_comparisons(
            proportions=prop,
            sample_ids=sample_ids,
            group_labels=labels,
            cell_type_names=results[0].cell_types,
            spec=self.group_comparison_spec,
            method=self.config.testing.method,
            fdr_alpha=self.config.testing.fdr_alpha,
            fdr_method=self.config.testing.fdr_method,
            posterior_samples_by_sample=posterior_by_sample,
        )

    # ------------------------------------------------------------------
    # Per-sample
    # ------------------------------------------------------------------

    def _marker_cpg_density(self, obs: MarkerObservations) -> np.ndarray | None:
        """Look up per-marker CpG density from region_annotations, if loaded."""
        if self.region_annotations is None or self.region_annotations.empty:
            return None
        try:
            ann = self.region_annotations.set_index(["chrom", "start", "end"])["cpg_density"]
        except KeyError:
            return None
        keys = list(zip(obs.chrom.tolist(), obs.start.tolist(), obs.end.tolist()))
        return np.array([float(ann.get(k, np.nan)) for k in keys], dtype=np.float64)

    def _load_and_calibrate(
        self,
        sample: Sample,
        marker_regions: MarkerRegions,
        cpg_index: dict | None,
    ) -> MarkerObservations:
        """Phase D: load methylation data and apply FinaleMe calibration only.

        Coverage tier assignment, observation model construction, and
        deconvolution are deferred to ``_deconvolve_sample`` so that batch
        correction and imputation can run in between.
        """
        log.info("Loading sample %s (%s)", sample.sample_id, sample.mode.value)
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
        if sample.mode == MeasurementMode.FINALEME and self.calibration is not None:
            obs = apply_calibration(obs, self.calibration, self.region_annotations)
        return obs

    def _deconvolve_sample(
        self,
        sample: Sample,
        obs: MarkerObservations,
        reference: ReferencePanel,
        out_dir: Path,
    ) -> DeconvolutionResult:
        """Build observation model, deconvolve, bootstrap, write per-sample TSV."""
        log.info("Deconvolving sample %s", sample.sample_id)

        calibration_flag: str | None = None
        if (
            sample.mode == MeasurementMode.FINALEME
            and self.calibration is not None
            and obs.predicted_beta is not None
        ):
            density_vec = self._marker_cpg_density(obs)
            qc = compute_inference_qc(
                obs.predicted_beta,
                self.calibration,
                cpg_density=density_vec,
            )
            calibration_flag = qc.get("flag")

        tier = self.tier_assigner.assign(obs)
        observation = self.observation_builder.build(
            obs=obs,
            reference=reference,
            calibration=self.calibration,
            region_annotations=self.region_annotations,
            tier=tier,
            coverage_cap=self.config.coverage.coverage_cap,
        )

        # Point-estimate proportions (always run MLE for the point estimate)
        w_hat = self.deconvolver.solve(observation, reference)

        # Uncertainty: Bayesian posterior overrides bootstrap when configured
        posterior_samples = None
        if self.bayesian_deconvolver is not None:
            try:
                posterior_samples = self.bayesian_deconvolver.solve(
                    observation, reference
                )
                w_hat = posterior_samples.mean(axis=0)
                ci_lo = np.quantile(
                    posterior_samples,
                    (1.0 - self.config.uncertainty.ci_level) / 2.0,
                    axis=0,
                )
                ci_hi = np.quantile(
                    posterior_samples,
                    1.0 - (1.0 - self.config.uncertainty.ci_level) / 2.0,
                    axis=0,
                )

                @dataclass
                class _BootShim:
                    proportions_samples: np.ndarray
                    ci_lower: np.ndarray
                    ci_upper: np.ndarray
                    point_estimate: np.ndarray

                boot = _BootShim(
                    proportions_samples=posterior_samples,
                    ci_lower=ci_lo,
                    ci_upper=ci_hi,
                    point_estimate=w_hat,
                )
            except Exception as exc:  # noqa: BLE001
                log.warning("Bayesian solve failed for %s, falling back to bootstrap: %s",
                            sample.sample_id, exc)
                boot = self.bootstrap.estimate(observation, reference, self.deconvolver)
        else:
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
            bootstrap_proportions=boot.proportions_samples if posterior_samples is None else None,
            posterior_samples=posterior_samples,
            coverage_tier=tier,
            qc_flags=[],  # filled below
        )
        result.qc_flags = compute_qc_flags(
            result=result,
            observation=observation,
            qc_config=self.config.qc,
            calibration_flag=calibration_flag,
            hemolysis=sample.metadata.get("hemolysis_flag")
            if hasattr(sample, "metadata") else None,
        )

        write_per_sample_too(result, out_dir / f"{sample.sample_id}.too.tsv")
        return result



# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _subset_reference(reference: ReferencePanel, indices: np.ndarray) -> ReferencePanel:
    return ReferencePanel(
        chrom=reference.chrom[indices],
        start=reference.start[indices],
        end=reference.end[indices],
        cell_types=list(reference.cell_types),
        methylation=reference.methylation[indices],
        coverage=reference.coverage[indices] if reference.coverage is not None else None,
    )


def _subset_marker_regions(mr: MarkerRegions, indices: np.ndarray) -> MarkerRegions:
    return MarkerRegions(
        chrom=mr.chrom[indices],
        start=mr.start[indices],
        end=mr.end[indices],
        marker_name=mr.marker_name[indices] if mr.marker_name is not None else None,
    )


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


def load_optional_calibration(
    config: TOOConfig,
    explicit_path: str | None,
    use_default: bool = True,
) -> CalibrationParams | None:
    """Load calibration parameters JSON if available, else default, else None."""
    path = explicit_path or config.calibration.calibration_file
    if path is not None:
        return CalibrationParams.load(path)
    if use_default and config.calibration.use_default:
        try:
            return load_default_calibration()
        except Exception as exc:  # pragma: no cover
            log.warning("Could not load default calibration: %s", exc)
            return None
    return None


def load_optional_region_annotations(
    config: TOOConfig, explicit_path: str | None
) -> pd.DataFrame | None:
    path = explicit_path or config.markers.region_annotation
    if path is None:
        return None
    return pd.read_csv(path, sep="\t", comment="#")


__all__ = [
    "CohortResult",
    "TOOPipeline",
    "build_reference_and_markers",
    "load_optional_calibration",
    "load_optional_region_annotations",
]
