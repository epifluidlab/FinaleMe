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
from finaleme_too.core.fragment_likelihood import FragmentLevelDeconvolver
from finaleme_too.core.observation_model import BetaBinomialModel, ObservationModel
from finaleme_too.core.reliability import (
    assign_reliability,
    compute_p_detection,
    compute_p_goodness,
)
from finaleme_too.core.uncertainty import BootstrapCI
from finaleme_too.io.pat_loader import load_fragments_from_pat
from finaleme_too.io.marker_regions import MarkerRegions, MarkerRegionsLoader
from finaleme_too.io.methylation_loader import MarkerObservations, MethylationLoader
from finaleme_too.io.output_writer import (
    write_cohort_proportions,
    write_group_comparison,
    write_per_sample_too,
    write_qc_summary,
    write_residual_analysis,
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
from finaleme_too.preprocessing.coverage import (
    CoverageTierAssigner,
    per_marker_min_reads,
    per_marker_min_reads_vector,
)
from finaleme_too.preprocessing.imputation import CohortImputer
from finaleme_too.preprocessing.marker_selection import BalancedMarkerSelector
from finaleme_too.postprocessing.covariate_adjustment import adjust_covariates
from finaleme_too.postprocessing.group_comparison import run_group_comparisons
from finaleme_too.postprocessing.qc import compute_qc_flags
from finaleme_too.postprocessing.residual_analysis import (
    compute_residuals_per_sample,
    discover_residual_components,
)
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
        cpg_index: dict | None = None,
    ):
        from finaleme_too.config import SolverMethod

        self.config = config
        self.calibration = calibration
        self.region_annotations = region_annotations
        self.group_comparison_spec = group_comparison_spec
        self.cpg_index = cpg_index
        self.deconvolver = MLEDeconvolver()
        self.fragment_deconvolver = FragmentLevelDeconvolver()

        # Decouple uncertainty source from the point-estimate solver.
        # ``uncertainty.method`` selects bootstrap / bayesian / both / none;
        # ``model.deconvolution`` selects MLE / Bayesian for the point estimate.
        uncertainty_method = (config.uncertainty.method or "bootstrap").lower()
        self._wants_bootstrap = uncertainty_method in ("bootstrap", "both")
        self._wants_bayesian_uncertainty = uncertainty_method in (
            "bayesian", "both",
        )
        self._wants_any_uncertainty = uncertainty_method != "none"

        # Instantiate BayesianDeconvolver if needed either for the point
        # estimate (model.deconvolution==BAYESIAN) or for uncertainty
        # (uncertainty.method in {bayesian, both}).
        needs_bayesian = (
            config.model.deconvolution == SolverMethod.BAYESIAN
            or self._wants_bayesian_uncertainty
        )
        self.bayesian_deconvolver = (
            BayesianDeconvolver(
                n_walkers=64,
                n_steps=max(config.uncertainty.bayesian_n_samples // 64, 1) + 100,
                burn_in=100,
                prior_alpha=config.uncertainty.bayesian_prior_alpha,
            )
            if needs_bayesian
            else None
        )
        self.bootstrap = BootstrapCI(
            n_iterations=config.uncertainty.n_bootstrap,
            ci_level=config.uncertainty.ci_level,
        )
        # Track whether model.deconvolution is Bayesian so the solver
        # dispatcher can prefer the posterior mean as the point estimate.
        self._point_estimate_is_bayesian = (
            config.model.deconvolution == SolverMethod.BAYESIAN
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
        # Cache for fragment-level path (ULTRALOW tier) + residual analysis
        if cpg_index is not None:
            self.cpg_index = cpg_index
        self._final_marker_regions = marker_regions

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

        # Optional: cohort imputation for low-coverage markers (LOW/ULTRALOW).
        # Snapshot pre-imputation counts so we can compute pct_imputed later.
        pre_impute_n: list[np.ndarray] = [
            np.asarray(obs.n, dtype=np.int64).copy() for obs in loaded_obs
        ]
        sample_groups_map = {s.sample_id: s.group for s in sample_sheet.samples}
        if self._should_impute(loaded_obs):
            imputer = CohortImputer(coverage_threshold=self.config.coverage.min_reads)
            loaded_obs = [
                imputer.impute(obs, loaded_obs, sample_groups_map) for obs in loaded_obs
            ]

        # Per-sample deconvolution (parallel)
        def _decon(
            args: tuple[Sample, MarkerObservations, np.ndarray],
        ) -> DeconvolutionResult:
            sample, obs, pre_n = args
            return self._deconvolve_sample(sample, obs, reference, out_dir, pre_n)

        results = parallel_map(
            _decon,
            list(zip(sample_sheet.samples, loaded_obs, pre_impute_n)),
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

        # Group comparison (P1): delegate sample-sufficiency checks to the
        # per-test helpers in run_group_comparisons. Bayesian comparisons
        # can produce valid rows with as few as 3 samples; ILR regression
        # needs >=2 per group but not an overall cohort size threshold.
        if self.group_comparison_spec and len(results) >= 2:
            test_results = self._run_group_comparisons(adjusted_results, sample_groups)
            if test_results:
                write_group_comparison(test_results, out_dir / "group_comparison.tsv")

        # Residual analysis + NMF discovery (Gap 6b, architecture §9.4)
        self._run_residual_analysis(adjusted_results, sample_groups, out_dir)

        return CohortResult(samples=adjusted_results, sample_groups=sample_groups)

    def _run_residual_analysis(
        self,
        results: list[DeconvolutionResult],
        sample_groups: dict[str, str | None],
        out_dir: Path,
    ) -> None:
        """Cohort-level NMF residual discovery + per-sample residual TSV.

        Always emits ``residual_analysis.tsv`` as long as ``results`` is not
        empty (architecture §10.5). The NMF discovery step is optional and
        only runs if enough samples carry residuals + high unknown fraction.
        """
        if not results:
            return

        # Pick the reference residual length from the FIRST sample that
        # actually has residuals. A per-sample fallback (e.g. the sample
        # hit _emit_fallback_result) used to abort the whole analysis.
        key_len = 0
        for r in results:
            if r.residuals is not None and r.residuals.size > 0:
                key_len = int(r.residuals.size)
                break

        nmf_summary: dict | None = None
        if key_len > 0:
            residual_matrix = np.zeros((len(results), key_len), dtype=np.float64)
            sample_order: list[str] = []
            for i, r in enumerate(results):
                sample_order.append(r.sample_id)
                if r.residuals is not None and r.residuals.size == key_len:
                    vec = r.residuals.copy()
                    vec[~np.isfinite(vec)] = 0.0
                    residual_matrix[i] = vec

            high_unknown = sum(
                1
                for r in results
                if r.proportions[-1] > self.config.qc.max_unknown_fraction
            )
            if high_unknown >= 3 and residual_matrix.shape[0] >= 3:
                try:
                    nmf = discover_residual_components(residual_matrix, n_components=3)
                    nmf["sample_order"] = sample_order
                    nmf_summary = nmf
                    log.info(
                        "NMF residual discovery: %d components, explained variance %s",
                        nmf["components"].shape[0],
                        np.round(nmf["explained_variance_ratio"], 3).tolist(),
                    )
                except Exception as exc:  # noqa: BLE001
                    log.warning("NMF residual discovery failed: %s", exc)

        # Always write the per-sample residual TSV. Samples that had no
        # residuals (e.g. tier-filter fallback) show up with NaN stats.
        write_residual_analysis(
            results,
            sample_groups,
            out_dir / "residual_analysis.tsv",
            nmf_summary=nmf_summary,
        )

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
        # Use dataclasses.replace so every enriched field on the original
        # result (mean_dispersion, mean_coverage, n_markers_used, pct_imputed,
        # calibration_flag, hemolysis_flag, overall_qc, residuals, marker_*)
        # is preserved after covariate adjustment. Only the fields we
        # actually want to change are overridden.
        from dataclasses import replace as dc_replace

        out: list[DeconvolutionResult] = []
        for r, new_p in zip(results, adjusted):
            out.append(
                dc_replace(
                    r,
                    proportions=new_p,
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
        pre_impute_n: np.ndarray | None = None,
    ) -> DeconvolutionResult:
        """Build observation model, deconvolve, bootstrap, write per-sample TSV.

        ``pre_impute_n`` is the pre-imputation per-marker read count (n) for
        this sample. It is used to compute ``pct_imputed`` in the result.
        """
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

        # -------- Per-marker effective-coverage filtering --------
        # A marker with below-expected coverage is down-tiered to the next
        # lower tier for that sample (architecture §4.2), so its minimum
        # read threshold becomes less strict. Markers that still fall below
        # their effective tier's minimum are filtered before the observation
        # model is built so the solver, bootstrap, and reliability
        # computations all see the filtered marker set uniformly.
        min_reads_vec = per_marker_min_reads_vector(
            np.asarray(obs.n, dtype=np.int64), tier
        )
        valid_mask = np.asarray(obs.n, dtype=np.int64) >= min_reads_vec
        if not np.any(valid_mask):
            # No usable markers at this tier — return an all-unknown result
            return self._emit_fallback_result(
                sample=sample,
                reference=reference,
                tier=tier,
                calibration_flag=calibration_flag,
                pre_impute_n=pre_impute_n,
                post_impute_n=np.asarray(obs.n, dtype=np.int64),
                out_dir=out_dir,
            )
        obs_filtered = _subset_observations(obs, valid_mask)
        reference_filtered = _subset_reference_rows(reference, valid_mask)
        if pre_impute_n is not None:
            pre_impute_n_filtered = pre_impute_n[valid_mask]
        else:
            pre_impute_n_filtered = None

        observation = self.observation_builder.build(
            obs=obs_filtered,
            reference=reference_filtered,
            calibration=self.calibration,
            region_annotations=self.region_annotations,
            tier=tier,
            coverage_cap=self.config.coverage.coverage_cap,
        )
        reference = reference_filtered  # downstream code operates on the filtered reference

        # -------- Gap 6a: fragment-level dispatch for ULTRALOW tier --------
        # When the sample is ULTRALOW, attempt to load a .pat.gz companion
        # and run the fragment-level EM. Falls back silently to MLE if the
        # prerequisites are missing.
        w_hat_fragment: np.ndarray | None = None
        fragment_mode = (self.config.model.fragment_level or "auto").lower()
        should_try_fragment = (
            fragment_mode == "always"
            or (fragment_mode == "auto" and tier == CoverageTier.ULTRALOW)
        )
        if should_try_fragment and self.cpg_index is not None:
            pat_path = sample.resolved_pat_file() if hasattr(sample, "resolved_pat_file") else None
            if pat_path is not None and Path(pat_path).exists():
                try:
                    frag_marker_regions = MarkerRegions(
                        chrom=reference.chrom,
                        start=reference.start,
                        end=reference.end,
                        marker_name=None,
                    )
                    fragments = load_fragments_from_pat(
                        pat_path,
                        marker_regions=frag_marker_regions,
                        cpg_index=self.cpg_index,
                    )
                    if fragments:
                        w_hat_fragment = self.fragment_deconvolver.solve(
                            fragments, reference.methylation
                        )
                        log.info(
                            "Sample %s: used fragment-level EM on %d fragments (tier=%s)",
                            sample.sample_id, len(fragments), tier.value,
                        )
                except Exception as exc:  # noqa: BLE001
                    log.warning(
                        "Fragment-level EM failed for %s (%s); falling back to MLE",
                        sample.sample_id, exc,
                    )

        # Point-estimate proportions (always run MLE for the point estimate).
        # For ULTRALOW + fragment path, we use the fragment w as the point
        # estimate but still run MLE on the marker-aggregated data as a
        # reference comparison (the bootstrap needs something to iterate over).
        w_hat = self.deconvolver.solve(observation, reference)
        if w_hat_fragment is not None:
            w_hat = w_hat_fragment

        # Uncertainty dispatch honors ``config.uncertainty.method``:
        #   - bootstrap : marker-resampling bootstrap on the MLE
        #   - bayesian  : MCMC posterior quantiles
        #   - both      : run bootstrap AND MCMC (bootstrap is primary,
        #                 posterior samples are stored alongside)
        #   - none      : skip uncertainty, CIs = point estimate
        #
        # Separately, ``config.model.deconvolution == BAYESIAN`` says the
        # point estimate should come from the posterior mean (regardless
        # of which uncertainty source drives the CIs).
        posterior_samples: np.ndarray | None = None
        boot = None

        @dataclass
        class _BootShim:
            proportions_samples: np.ndarray
            ci_lower: np.ndarray
            ci_upper: np.ndarray
            point_estimate: np.ndarray

        def _posterior_to_boot(samples: np.ndarray) -> _BootShim:
            lo = np.quantile(
                samples, (1.0 - self.config.uncertainty.ci_level) / 2.0, axis=0
            )
            hi = np.quantile(
                samples, 1.0 - (1.0 - self.config.uncertainty.ci_level) / 2.0, axis=0
            )
            return _BootShim(
                proportions_samples=samples,
                ci_lower=lo,
                ci_upper=hi,
                point_estimate=samples.mean(axis=0),
            )

        # 1) Run Bayesian posterior if we need it (for point or uncertainty)
        if self.bayesian_deconvolver is not None:
            try:
                posterior_samples = self.bayesian_deconvolver.solve(
                    observation, reference
                )
                if self._point_estimate_is_bayesian:
                    w_hat = posterior_samples.mean(axis=0)
            except Exception as exc:  # noqa: BLE001
                log.warning(
                    "Bayesian solve failed for %s, falling back to bootstrap: %s",
                    sample.sample_id, exc,
                )
                posterior_samples = None

        # 2) Pick the "primary" boot struct that drives reliability + CIs.
        if self._wants_bootstrap:
            boot = self.bootstrap.estimate(observation, reference, self.deconvolver)
        elif self._wants_bayesian_uncertainty and posterior_samples is not None:
            boot = _posterior_to_boot(posterior_samples)
        elif not self._wants_any_uncertainty:
            # uncertainty.method == "none": zero-width CIs around the point
            zero_samples = np.tile(w_hat, (1, 1))
            boot = _BootShim(
                proportions_samples=zero_samples,
                ci_lower=w_hat.copy(),
                ci_upper=w_hat.copy(),
                point_estimate=w_hat.copy(),
            )
        else:
            # Requested bayesian/both but Bayesian failed → fall back to bootstrap
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
        valid_n = int(np.sum(observation.n > 0))
        n_markers = np.full(K, valid_n, dtype=np.int32)

        # Enriched output fields (Gap 7)
        mean_dispersion = _per_celltype_mean_dispersion(
            observation=observation,
            reference_methylation=reference.methylation,
            top_n=50,
        )
        mean_coverage = float(np.mean(observation.n)) if observation.n.size > 0 else 0.0

        # pct_imputed: fraction of markers where the pre-imputation n was below
        # the per-marker tier threshold but the post-imputation n is at or
        # above it.
        if pre_impute_n_filtered is not None and pre_impute_n_filtered.size:
            post = np.asarray(observation.n, dtype=np.int64)
            # Recompute the per-marker threshold vector on the filtered subset
            post_min_reads = per_marker_min_reads_vector(
                np.asarray(observation.n, dtype=np.int64), tier
            )
            was_low = pre_impute_n_filtered < post_min_reads
            was_lifted = was_low & (post >= post_min_reads)
            pct_imputed = float(np.mean(was_lifted))
        else:
            pct_imputed = 0.0

        # Per-marker residual (mu_obs - mu_predicted) for NMF later
        with np.errstate(invalid="ignore", divide="ignore"):
            mu_obs = np.where(
                observation.n > 0,
                observation.k / np.maximum(observation.n, 1),
                np.nan,
            )
        mu_pred = reference.methylation @ w_hat[:K]
        residuals = mu_obs - mu_pred

        hemolysis_flag = None
        if hasattr(sample, "metadata"):
            hv = sample.metadata.get("hemolysis_flag")
            if hv is not None:
                hemolysis_flag = bool(hv) if not isinstance(hv, str) else hv.lower() in ("true", "1", "yes")

        # Bootstrap vs posterior samples are stored independently: bootstrap
        # samples are populated iff we ran the bootstrap; posterior samples
        # are populated iff we ran the Bayesian deconvolver. Both can coexist
        # when uncertainty.method == "both".
        bootstrap_samples_field: np.ndarray | None = None
        if self._wants_bootstrap and boot is not None:
            bootstrap_samples_field = boot.proportions_samples

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
            bootstrap_proportions=bootstrap_samples_field,
            posterior_samples=posterior_samples,
            coverage_tier=tier,
            qc_flags=[],  # filled below
            mean_dispersion=mean_dispersion,
            mean_coverage=mean_coverage,
            n_markers_used=valid_n,
            pct_imputed=pct_imputed,
            calibration_flag=calibration_flag,
            hemolysis_flag=hemolysis_flag,
            residuals=residuals,
            marker_chrom=np.asarray(obs_filtered.chrom, dtype=object),
            marker_start=np.asarray(obs_filtered.start, dtype=np.int64),
            marker_end=np.asarray(obs_filtered.end, dtype=np.int64),
        )
        result.qc_flags = compute_qc_flags(
            result=result,
            observation=observation,
            qc_config=self.config.qc,
            calibration_flag=calibration_flag,
            hemolysis=hemolysis_flag,
        )
        # overall_qc derived from qc_flags
        if any(f.endswith("_FAIL") for f in result.qc_flags):
            result.overall_qc = "FAIL"
        elif result.qc_flags:
            result.overall_qc = "WARN"
        else:
            result.overall_qc = "PASS"

        write_per_sample_too(result, out_dir / f"{sample.sample_id}.too.tsv")
        return result

    def _emit_fallback_result(
        self,
        sample: Sample,
        reference: ReferencePanel,
        tier: CoverageTier,
        calibration_flag: str | None,
        pre_impute_n: np.ndarray | None,
        post_impute_n: np.ndarray,
        out_dir: Path,
    ) -> DeconvolutionResult:
        """Return an all-unknown result when no markers pass the tier filter."""
        K = reference.n_cell_types
        w = np.zeros(K + 1, dtype=np.float64)
        w[-1] = 1.0
        reliability = np.array(["UNRELIABLE"] * (K + 1), dtype=object)
        result = DeconvolutionResult(
            sample_id=sample.sample_id,
            cell_types=list(reference.cell_types),
            proportions=w,
            ci_lower=w.copy(),
            ci_upper=w.copy(),
            p_goodness=np.full(K, np.nan, dtype=np.float64),
            p_detection=np.zeros(K + 1, dtype=np.float64),
            reliability=reliability,
            n_markers=np.zeros(K, dtype=np.int32),
            coverage_tier=tier,
            qc_flags=["NO_MARKERS_PASS_TIER_FILTER"],
            mean_dispersion=np.zeros(K, dtype=np.float64),
            mean_coverage=float(np.mean(post_impute_n)) if post_impute_n.size else 0.0,
            n_markers_used=0,
            pct_imputed=0.0,
            calibration_flag=calibration_flag,
            hemolysis_flag=None,
            overall_qc="FAIL",
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


def _subset_observations(
    obs: MarkerObservations, mask: np.ndarray
) -> MarkerObservations:
    """Return a new MarkerObservations keeping only the masked markers."""
    return MarkerObservations(
        sample_id=obs.sample_id,
        chrom=obs.chrom[mask],
        start=obs.start[mask],
        end=obs.end[mask],
        k=np.asarray(obs.k, dtype=np.int32)[mask],
        n=np.asarray(obs.n, dtype=np.int32)[mask],
        predicted_beta=(
            obs.predicted_beta[mask] if obs.predicted_beta is not None else None
        ),
        mode=obs.mode,
    )


def _subset_reference_rows(
    reference: ReferencePanel, mask: np.ndarray
) -> ReferencePanel:
    """Return a new ReferencePanel keeping only the masked rows (markers)."""
    return ReferencePanel(
        chrom=reference.chrom[mask],
        start=reference.start[mask],
        end=reference.end[mask],
        cell_types=list(reference.cell_types),
        methylation=reference.methylation[mask],
        coverage=reference.coverage[mask] if reference.coverage is not None else None,
    )


def _per_celltype_mean_dispersion(
    observation,
    reference_methylation: np.ndarray,
    top_n: int = 50,
) -> np.ndarray:
    """Per-cell-type mean dispersion phi over the cell type's top discriminative markers.

    The dispersion scalar reported in the per-sample TSV (architecture §10.1)
    summarizes the model uncertainty attached to each cell type's best markers.
    """
    R = np.asarray(reference_methylation, dtype=np.float64)
    K = R.shape[1]
    out = np.zeros(K, dtype=np.float64)
    if observation.dispersion.size == 0:
        return out
    valid = observation.n > 0
    for j in range(K):
        target = R[:, j]
        others = np.delete(R, j, axis=1)
        bg_mean = np.mean(others, axis=1)
        score = np.where(valid, np.abs(target - bg_mean), -np.inf)
        if np.all(np.isneginf(score)):
            out[j] = float(np.mean(observation.dispersion)) if observation.dispersion.size else 0.0
            continue
        take = min(top_n, int(np.sum(valid)))
        top_idx = np.argpartition(-score, take - 1)[:take]
        out[j] = float(np.mean(observation.dispersion[top_idx]))
    return out


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
