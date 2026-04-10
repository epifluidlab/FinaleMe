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
from time import perf_counter
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
    compute_fit_metrics,
    compute_p_detection,
    compute_unknown_fit_metrics,
)
from finaleme_too.core.uncertainty import BootstrapCI
from finaleme_too.io.pat_loader import load_fragments_from_pat
from finaleme_too.io.marker_regions import MarkerRegions, MarkerRegionsLoader
from finaleme_too.io.methylation_loader import MarkerObservations, MethylationLoader
from finaleme_too.io.output_writer import (
    write_binarization_debug,
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
from finaleme_too.preprocessing.batch_correction import (
    combat_correct,
    combat_correct_predicted_beta,
)
from finaleme_too.preprocessing.coverage import (
    CoverageTierAssigner,
    effective_coverage_in_markers,
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
    """End-to-end TOO pipeline.

    FinaleMe-mode samples are processed through the v3 context-dependent
    binarization model: ``apply_binarization`` classifies markers into
    U / M / Ambiguous / Excluded, then ``BinarizationModel.build`` produces
    the observation model that ``MLEDeconvolver.solve`` consumes via its
    binomial-with-error-rates SLSQP path.

    WGBS-mode samples continue to use the beta-binomial path regardless of
    whether a ``binarization`` object is provided.

    The legacy ``calibration`` parameter was removed in v3. The old v2
    continuous-calibration pipeline is no longer available; use
    ``finaleme-too train-binarization`` to produce a ``BinarizationParams``
    JSON and pass it via ``binarization=`` or ``--binarization``.
    """

    def __init__(
        self,
        config: TOOConfig,
        region_annotations: pd.DataFrame | None = None,
        group_comparison_spec: str | None = None,
        cpg_index: dict | None = None,
        binarization=None,  # BinarizationParams | None — v3 FinaleMe path
        binarize_threshold: float | None = None,  # universal hard threshold mode
        binarize_threshold_from_params: bool = False,  # per-bin tau_high mode
        binarize_use_all_bins: bool = False,  # include unusable bins in tau_high mode
    ):
        from finaleme_too.config import SolverMethod
        from finaleme_too.core.observation_model_binarization import BinarizationModel

        self.config = config
        self.binarization = binarization
        if binarize_threshold is not None:
            thr = float(binarize_threshold)
            if thr < 0.0 or thr > 1.0:
                raise ValueError(
                    f"binarize_threshold must be within [0, 1], got {thr}"
                )
            self.binarize_threshold = thr
        else:
            self.binarize_threshold = None
        self.binarize_threshold_from_params = bool(binarize_threshold_from_params)
        self.binarize_use_all_bins = bool(binarize_use_all_bins)
        self.region_annotations = region_annotations
        self._region_annotation_lookup: pd.DataFrame | None = None
        self.group_comparison_spec = group_comparison_spec
        self.cpg_index = cpg_index
        self.deconvolver = MLEDeconvolver(
            binarization_count_weight=config.model.binarization_count_weight,
            binarization_model=config.model.binarization_model,
            hierarchical_quadrature_points=config.model.hierarchical_quadrature_points,
        )
        self.fragment_deconvolver = FragmentLevelDeconvolver()
        # v3 binarization observation-model builder. Cheap to construct
        # whether or not binarization is provided; the actual dispatch in
        # _deconvolve_sample only uses it when self.binarization is not None.
        self._binarization_builder = BinarizationModel(
            hard_threshold=self.binarize_threshold,
            binarization_model=config.model.binarization_model,
            call_weight_override=config.model.hierarchical_call_weight,
            learned_threshold_from_params=self.binarize_threshold_from_params,
            learned_threshold_use_all_bins=self.binarize_use_all_bins,
        )

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
                seed=config.uncertainty.seed,
                binarization_count_weight=config.model.binarization_count_weight,
                binarization_model=config.model.binarization_model,
                hierarchical_quadrature_points=config.model.hierarchical_quadrature_points,
            )
            if needs_bayesian
            else None
        )
        self.bootstrap = BootstrapCI(
            n_iterations=config.uncertainty.n_bootstrap,
            ci_level=config.uncertainty.ci_level,
            seed=config.uncertainty.seed,
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
        # Region annotation lookup can dominate runtime when rebuilt for
        # every sample. Pre-index once and reuse across all binarization calls.
        if region_annotations is not None and not region_annotations.empty:
            try:
                from finaleme_too.preprocessing.binarization import (
                    prepare_region_annotation_lookup,
                )

                self._region_annotation_lookup = prepare_region_annotation_lookup(
                    region_annotations
                )
            except Exception as exc:  # noqa: BLE001
                log.warning(
                    "Failed to pre-index region annotations (%s); falling back to per-sample path.",
                    exc,
                )
                self._region_annotation_lookup = region_annotations

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
        t_run_start = perf_counter()
        stage_seconds: dict[str, float] = {}
        skipped_stages: set[str] = set()

        def _mark_stage(name: str, t0: float, skipped: bool = False) -> None:
            stage_seconds[name] = perf_counter() - t0
            if skipped:
                skipped_stages.add(name)

        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        # Cache for fragment-level path (ULTRALOW tier) + residual analysis
        if cpg_index is not None:
            self.cpg_index = cpg_index
        self._final_marker_regions = marker_regions

        # Cell-type-balanced marker selection (P1)
        t0 = perf_counter()
        marker_subset_indices: np.ndarray | None = None
        marker_selection_skipped = True
        if self.marker_selector is not None and reference.n_markers > self.marker_selector.n_per_type:
            marker_selection_skipped = False
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
        _mark_stage("marker_selection", t0, skipped=marker_selection_skipped)

        # Phase D: load raw samples first — no preprocessing yet. We need
        # cohort-level operations (batch correction, imputation, and for
        # FinaleMe mode binarization itself) to run AFTER loading, because
        # batch correction of FinaleMe predicted_beta must happen BEFORE
        # binarization (math doc §10). Otherwise the per-marker state calls
        # would be derived from pre-correction betas and then the corrected
        # k/n would never feed the FinaleMe likelihood.
        #
        # Use the threading backend: pandas.read_csv / numpy releases the GIL
        # in the inner loops, and the marker_regions / cpg_index arguments
        # are large — pickling them per task (the loky default) would
        # dominate the runtime.
        t0 = perf_counter()
        log.info("Loading %d samples", len(sample_sheet.samples))
        loaded_obs: list[MarkerObservations] = parallel_map(
            lambda s: self._load_sample(s, marker_regions, cpg_index),
            sample_sheet.samples,
            n_jobs=self.config.threads,
            backend="threading",
        )
        _mark_stage("load_samples", t0)

        # Optional: ComBat-style technical batch correction. Mode-aware:
        #   * WGBS samples → combat_correct on k/n (operates on beta = k/n)
        #   * FinaleMe samples → combat_correct_predicted_beta on
        #     predicted_beta (BEFORE binarization, so the corrected
        #     predictions drive the U / M / Ambiguous state calls).
        # We run both correctors because a cohort can mix the two modes;
        # each only touches samples it applies to.
        t0 = perf_counter()
        batch_correction_skipped = True
        if self.config.batch_correction.technical_covariates:
            batch_correction_skipped = False
            batch_var = self.config.batch_correction.technical_covariates[0]
            batch_labels = [
                str(s.metadata.get(batch_var)) if batch_var in s.metadata else None
                for s in sample_sheet.samples
            ]
            log.info("Applying batch correction on %s", batch_var)
            # WGBS path: correct k/n. The FinaleMe-mode samples in the list
            # have predicted_beta set, so combat_correct (which reads k/n)
            # would also mutate them, but their likelihood goes through the
            # binarization path — k/n only feed the coverage weight, not
            # the state call. To keep things clean we run each corrector
            # over only its own mode's samples.
            wgbs_labels = [
                lbl if obs.mode == MeasurementMode.WGBS else None
                for obs, lbl in zip(loaded_obs, batch_labels)
            ]
            finaleme_labels = [
                lbl if obs.mode == MeasurementMode.FINALEME else None
                for obs, lbl in zip(loaded_obs, batch_labels)
            ]
            loaded_obs = combat_correct(
                loaded_obs,
                batch_labels=wgbs_labels,
                min_levels=self.config.batch_correction.min_levels,
                min_per_level=self.config.batch_correction.min_samples_per_level,
            )
            loaded_obs = combat_correct_predicted_beta(
                loaded_obs,
                batch_labels=finaleme_labels,
                min_levels=self.config.batch_correction.min_levels,
                min_per_level=self.config.batch_correction.min_samples_per_level,
            )
        _mark_stage("batch_correction", t0, skipped=batch_correction_skipped)

        # Now that FinaleMe samples have their predicted_beta corrected,
        # apply binarization. This is where called_state + context_bin are
        # populated — downstream code in _deconvolve_sample reads these.
        t0 = perf_counter()
        binarization_skipped = True
        if self.binarization is not None:
            binarization_skipped = False
            loaded_obs = parallel_map(
                self._apply_finaleme_binarization,
                loaded_obs,
                n_jobs=self.config.threads,
                backend="threading",
            )
        _mark_stage("apply_binarization", t0, skipped=binarization_skipped)

        # Optional: cohort imputation for low-coverage markers (LOW/ULTRALOW).
        # Snapshot pre-imputation counts so we can compute pct_imputed later.
        pre_impute_n: list[np.ndarray] = [
            np.asarray(obs.n, dtype=np.int64).copy() for obs in loaded_obs
        ]
        sample_groups_map = {s.sample_id: s.group for s in sample_sheet.samples}
        t0 = perf_counter()
        imputation_skipped = True
        if self._should_impute(loaded_obs):
            imputation_skipped = False
            imputer = CohortImputer(coverage_threshold=self.config.coverage.min_reads)
            # Imputation is pure numpy on per-sample arrays — releases the GIL
            # and shares the cohort donor list by reference, so the threading
            # backend wins by avoiding the per-task pickle of all donors.
            loaded_obs = parallel_map(
                lambda obs: imputer.impute(obs, loaded_obs, sample_groups_map),
                loaded_obs,
                n_jobs=self.config.threads,
                backend="threading",
            )
        _mark_stage("imputation", t0, skipped=imputation_skipped)

        # Compute the per-sample inner-thread count. When the cohort is
        # smaller than --threads (e.g. one sample run with --threads 8), the
        # outer parallel_map can only saturate `len(samples)` cores. Hand the
        # remaining budget to the bootstrap inner loop so all of --threads
        # are actually used.
        n_samples = len(sample_sheet.samples)
        outer_jobs = max(1, min(self.config.threads, n_samples))
        inner_jobs = max(1, self.config.threads // outer_jobs)

        # Per-sample deconvolution (parallel)
        def _decon(
            args: tuple[Sample, MarkerObservations, np.ndarray],
        ) -> DeconvolutionResult:
            sample, obs, pre_n = args
            return self._deconvolve_sample(
                sample, obs, reference, out_dir, pre_n, bootstrap_jobs=inner_jobs
            )

        # Threading backend: scipy.optimize SLSQP, the bootstrap, and the
        # observation-model construction are pure numpy/scipy and release the
        # GIL. Loky here would re-pickle the (large) reference panel for every
        # sample, which is what was capping CPU usage well below the user's
        # `--threads` setting.
        t0 = perf_counter()
        results = parallel_map(
            _decon,
            list(zip(sample_sheet.samples, loaded_obs, pre_impute_n)),
            n_jobs=outer_jobs,
            backend="threading",
        )
        _mark_stage("deconvolution", t0)

        # Phase D: optional ILR-space biological covariate adjustment
        t0 = perf_counter()
        adjusted_results = self._maybe_adjust_covariates(results, sample_sheet)
        _mark_stage(
            "covariate_adjustment",
            t0,
            skipped=adjusted_results is results,
        )

        # Cohort outputs
        t0 = perf_counter()
        sample_groups = {s.sample_id: s.group for s in sample_sheet.samples}
        write_cohort_proportions(
            adjusted_results, sample_groups, out_dir / "cohort_proportions.tsv"
        )
        write_qc_summary(adjusted_results, sample_groups, out_dir / "qc_summary.tsv")
        write_binarization_debug(
            adjusted_results, sample_groups, out_dir / "binarization_debug.tsv"
        )
        _mark_stage("write_primary_outputs", t0)

        # Group comparison (P1): delegate sample-sufficiency checks to the
        # per-test helpers in run_group_comparisons. Bayesian comparisons
        # can produce valid rows with as few as 3 samples; ILR regression
        # needs >=2 per group but not an overall cohort size threshold.
        t0 = perf_counter()
        group_comparison_skipped = True
        if self.group_comparison_spec and len(results) >= 2:
            group_comparison_skipped = False
            test_results = self._run_group_comparisons(
                adjusted_results, sample_groups, sample_sheet
            )
            if test_results:
                write_group_comparison(test_results, out_dir / "group_comparison.tsv")
        _mark_stage("group_comparison", t0, skipped=group_comparison_skipped)

        # Residual analysis + NMF discovery (Gap 6b, architecture §9.4)
        t0 = perf_counter()
        self._run_residual_analysis(adjusted_results, sample_groups, out_dir)
        _mark_stage("residual_analysis", t0)

        stage_order = [
            "marker_selection",
            "load_samples",
            "batch_correction",
            "apply_binarization",
            "imputation",
            "deconvolution",
            "covariate_adjustment",
            "write_primary_outputs",
            "group_comparison",
            "residual_analysis",
        ]
        stage_parts: list[str] = []
        for name in stage_order:
            sec = stage_seconds.get(name, 0.0)
            if name in skipped_stages:
                stage_parts.append(f"{name}=skipped({sec:.2f}s)")
            else:
                stage_parts.append(f"{name}={sec:.2f}s")
        log.info(
            "Pipeline runtime summary: %s | total=%.2fs",
            ", ".join(stage_parts),
            perf_counter() - t_run_start,
        )

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

    def _resolve_biological_covariates(
        self,
        sample_sheet: SampleSheet,
    ) -> list[str]:
        """Resolve the active biological covariate list for this run.

        Per architecture §9.1: biological covariates that are populated on
        at least one sample are adjusted by default. The user-configurable
        covariates (``treatment``, ``treatment_efficacy``, ``mutation_status``)
        stay opt-in — they're only included when the user explicitly
        configures them via ``config.covariate_adjustment.biological_covariates``.

        Resolution order:
          1. If ``config.covariate_adjustment.biological_covariates`` is
             non-empty, use it verbatim (explicit user config wins).
          2. Otherwise, auto-populate from the sample sheet: any key on
             any sample's ``biological_covariates`` dict that is NOT in
             ``config.covariate_adjustment.user_configurable`` is included.

        Returns the resolved list (possibly empty if no samples carry any
        biological-covariate metadata).
        """
        configured = list(self.config.covariate_adjustment.biological_covariates)
        if configured:
            return configured

        user_configurable = set(self.config.covariate_adjustment.user_configurable)
        present: set[str] = set()
        for s in sample_sheet.samples:
            for k, v in s.biological_covariates.items():
                if k in user_configurable:
                    continue
                if v is None:
                    continue
                # Skip obviously-empty values so "auto" doesn't enable a
                # column that's NaN on every sample.
                if isinstance(v, float) and np.isnan(v):
                    continue
                if isinstance(v, str) and not v.strip():
                    continue
                present.add(k)
        # Return in a stable (alphabetical) order so different runs on the
        # same cohort produce identical regression designs.
        return sorted(present)

    def _maybe_adjust_covariates(
        self,
        results: list[DeconvolutionResult],
        sample_sheet: SampleSheet,
    ) -> list[DeconvolutionResult]:
        """Apply ILR-space biological covariate adjustment if configured."""
        bio = self._resolve_biological_covariates(sample_sheet)
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
        # binarization_flag, hemolysis_flag, overall_qc, residuals, marker_*)
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
        sample_sheet: SampleSheet,
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

        # Build a covariate DataFrame from the sample sheet for the ILR
        # regression. We include the same biological covariates the
        # post-deconvolution covariate adjustment step used so the group
        # comparison is testing "residual group effect after adjusting for
        # confounders". ``compositional_regression_test`` drops any
        # sample with NaN covariates on a per-fit basis, so partially-
        # populated covariates don't kill the whole analysis.
        covariate_df = self._build_covariate_dataframe(sample_sheet, sample_ids)

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
            covariates=covariate_df,
        )

    def _build_covariate_dataframe(
        self,
        sample_sheet: SampleSheet,
        sample_ids: list[str],
    ) -> pd.DataFrame | None:
        """Build a per-sample biological-covariate DataFrame for the ILR
        regression, using the same auto-populated list as
        ``_maybe_adjust_covariates`` so the covariate adjustment and the
        group comparison see the same design matrix.

        Returns ``None`` when neither config nor sample metadata provide
        any biological covariates, in which case
        ``compositional_regression_test`` skips the covariate columns
        entirely.
        """
        bio_cols = self._resolve_biological_covariates(sample_sheet)
        if not bio_cols:
            return None
        sample_map = {s.sample_id: s for s in sample_sheet.samples}
        rows = []
        for sid in sample_ids:
            sample = sample_map.get(sid)
            if sample is None:
                rows.append({c: np.nan for c in bio_cols})
                continue
            cov = sample.biological_covariates
            rows.append({c: cov.get(c, np.nan) for c in bio_cols})
        cov_df = pd.DataFrame(rows, index=sample_ids)
        cov_df.index.name = "sample_id"
        return cov_df

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

    def _load_sample(
        self,
        sample: Sample,
        marker_regions: MarkerRegions,
        cpg_index: dict | None,
    ) -> MarkerObservations:
        """Load raw methylation data for one sample (no preprocessing).

        Binarization (for FinaleMe mode) is deliberately NOT applied here —
        it must run **after** batch correction (math doc §10) so the
        corrected ``predicted_beta`` is what the binarizer classifies into
        U / M / Ambiguous / Excluded states. See ``_apply_finaleme_binarization``
        which runs per-sample after the cohort-level batch correction step.
        """
        log.info("Loading sample %s (%s)", sample.sample_id, sample.mode.value)
        return MethylationLoader.load(
            filepath=sample.methylation_file,
            sample_id=sample.sample_id,
            mode=sample.mode,
            marker_regions=marker_regions,
            input_format=sample.input_format,
            meth_col=sample.meth_col,
            total_col=sample.total_col,
            cpg_index=cpg_index,
        )

    def _apply_finaleme_binarization(
        self, obs: MarkerObservations
    ) -> MarkerObservations:
        """Apply the trained binarization to a single (possibly batch-corrected)
        FinaleMe sample. WGBS and un-binarized FinaleMe samples pass through
        unchanged.
        """
        if (
            obs.mode != MeasurementMode.FINALEME
            or self.binarization is None
            or obs.predicted_beta is None
        ):
            return obs
        from finaleme_too.preprocessing.binarization import apply_binarization

        return apply_binarization(
            obs,
            params=self.binarization,
            region_annotations=(
                self._region_annotation_lookup
                if self._region_annotation_lookup is not None
                else self.region_annotations
            ),
            hard_threshold=self.binarize_threshold,
            learned_threshold_from_params=self.binarize_threshold_from_params,
            learned_threshold_use_all_bins=self.binarize_use_all_bins,
        )

    # Backwards-compat aliases so any stale caller finding the old method
    # names on TOOPipeline still works. The v2 name was _load_and_calibrate;
    # the post-Phase-E name was _load_and_preprocess_sample. Both now route
    # to _load_sample which no longer binarizes (binarization moved to
    # _apply_finaleme_binarization called after batch correction).
    _load_and_preprocess_sample = _load_sample
    _load_and_calibrate = _load_sample

    def _deconvolve_sample(
        self,
        sample: Sample,
        obs: MarkerObservations,
        reference: ReferencePanel,
        out_dir: Path,
        pre_impute_n: np.ndarray | None = None,
        bootstrap_jobs: int = 1,
    ) -> DeconvolutionResult:
        """Build observation model, deconvolve, bootstrap, write per-sample TSV.

        ``pre_impute_n`` is the pre-imputation per-marker read count (n) for
        this sample. It is used to compute ``pct_imputed`` in the result.

        ``bootstrap_jobs`` controls inner-loop parallelism over bootstrap
        iterations. The pipeline computes this from
        ``threads // num_samples`` so a single-sample run with ``--threads 8``
        actually uses all 8 cores instead of just the one outer thread.
        """
        log.info("Deconvolving sample %s", sample.sample_id)

        # Whether the v3 binarization path is active for THIS sample. Both
        # ``self.binarization`` must be set AND the sample must be FinaleMe
        # mode AND apply_binarization must already have run on the obs (which
        # populates called_state).
        use_binarization = (
            sample.mode == MeasurementMode.FINALEME
            and self.binarization is not None
            and obs.called_state is not None
            and obs.context_bin is not None
        )
        tier = self.tier_assigner.assign(obs)

        binarization_flag: str | None = None
        binarization_debug: dict | None = None
        if use_binarization:
            from finaleme_too.preprocessing.binarization_eval import (
                compute_inference_qc as _binarization_inference_qc,
            )
            from finaleme_too.preprocessing.binarization import (
                STATE_AMBIGUOUS,
                STATE_EXCLUDED,
                STATE_M,
                STATE_U,
            )
            # In universal hard-threshold mode (--binarizeThreshold),
            # binarization inference-QC calibration diagnostics are not
            # meaningful; treat the binarization stage as bypassed for QC
            # flagging so users do not get BINARIZATION_FAIL from this mode.
            if self.binarize_threshold is not None or self.binarize_threshold_from_params:
                qc = {
                    "flag": "HARD_THRESHOLD",
                    "fraction_called": float("nan"),
                    "bin_balance": float("nan"),
                    "state_distribution_kl": float("nan"),
                }
                binarization_flag = "HARD_THRESHOLD"
            else:
                qc = _binarization_inference_qc(
                    called_state=obs.called_state,
                    context_bin=obs.context_bin,
                    params=self.binarization,
                )
                binarization_flag = qc.get("flag")

            # Per-sample debug scaffold for troubleshooting cases where the
            # deconvolution falls back to Unknown=1 due to having no callable
            # markers after filtering/binarization.
            raw_state = np.asarray(obs.called_state, dtype=np.int64)
            raw_bin = np.asarray(obs.context_bin, dtype=np.int64)
            usable_bins_mask = np.asarray(self.binarization.usable, dtype=bool)
            raw_usable_per_marker = (
                usable_bins_mask[raw_bin]
                if raw_bin.size
                else np.zeros(0, dtype=bool)
            )
            raw_bins_hit = np.unique(raw_bin) if raw_bin.size else np.zeros(0, dtype=np.int64)
            raw_usable_bins_hit = np.unique(raw_bin[raw_usable_per_marker]) if raw_bin.size else np.zeros(0, dtype=np.int64)
            binarization_debug = {
                "coverage_tier": tier.value if hasattr(tier, "value") else str(tier),
                "binarization_flag": binarization_flag or "NA",
                "fraction_called_qc": float(qc.get("fraction_called", float("nan"))),
                "bin_balance_qc": float(qc.get("bin_balance", float("nan"))),
                "state_distribution_kl_qc": float(qc.get("state_distribution_kl", float("nan"))),
                "n_bins_total": int(self.binarization.n_bins),
                "n_usable_bins_total": int(np.sum(usable_bins_mask)),
                "n_markers_before_coverage_filter": int(raw_state.size),
                "called_u_before_filter": int(np.sum(raw_state == STATE_U)),
                "called_m_before_filter": int(np.sum(raw_state == STATE_M)),
                "ambiguous_before_filter": int(np.sum(raw_state == STATE_AMBIGUOUS)),
                "excluded_before_filter": int(np.sum(raw_state == STATE_EXCLUDED)),
                "n_bins_hit_before_filter": int(raw_bins_hit.size),
                "n_usable_bins_hit_before_filter": int(raw_usable_bins_hit.size),
                # Fields below are updated after the coverage + observation-model
                # filtering steps.
                "n_markers_after_coverage_filter": 0,
                "called_u_after_filter": 0,
                "called_m_after_filter": 0,
                "ambiguous_after_filter": 0,
                "excluded_after_filter": 0,
                "n_bins_hit_after_filter": 0,
                "n_usable_bins_hit_after_filter": 0,
                "n_unusable_bins_hit_after_filter": 0,
                "n_markers_used_for_likelihood": 0,
                "fail_reason": "",
            }

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
            if binarization_debug is not None:
                binarization_debug["coverage_tier"] = (
                    tier.value if hasattr(tier, "value") else str(tier)
                )
                binarization_debug["fail_reason"] = "no_markers_pass_coverage_filter"
            return self._emit_fallback_result(
                sample=sample,
                obs=obs,
                reference=reference,
                tier=tier,
                binarization_flag=binarization_flag,
                pre_impute_n=pre_impute_n,
                out_dir=out_dir,
                binarization_debug=binarization_debug,
            )
        obs_filtered = _subset_observations(obs, valid_mask)
        reference_filtered = _subset_reference_rows(reference, valid_mask)
        if pre_impute_n is not None:
            pre_impute_n_filtered = pre_impute_n[valid_mask]
        else:
            pre_impute_n_filtered = None

        if use_binarization and binarization_debug is not None:
            from finaleme_too.preprocessing.binarization import (
                STATE_AMBIGUOUS,
                STATE_EXCLUDED,
                STATE_M,
                STATE_U,
            )

            st = np.asarray(obs_filtered.called_state, dtype=np.int64)
            bn = np.asarray(obs_filtered.context_bin, dtype=np.int64)
            usable_bins_mask = np.asarray(self.binarization.usable, dtype=bool)
            usable_per_marker = (
                usable_bins_mask[bn] if bn.size else np.zeros(0, dtype=bool)
            )
            bins_hit = np.unique(bn) if bn.size else np.zeros(0, dtype=np.int64)
            usable_bins_hit = np.unique(bn[usable_per_marker]) if bn.size else np.zeros(0, dtype=np.int64)
            unusable_bins_hit = np.unique(bn[~usable_per_marker]) if bn.size else np.zeros(0, dtype=np.int64)
            binarization_debug.update(
                {
                    "coverage_tier": tier.value if hasattr(tier, "value") else str(tier),
                    "n_markers_after_coverage_filter": int(st.size),
                    "called_u_after_filter": int(np.sum(st == STATE_U)),
                    "called_m_after_filter": int(np.sum(st == STATE_M)),
                    "ambiguous_after_filter": int(np.sum(st == STATE_AMBIGUOUS)),
                    "excluded_after_filter": int(np.sum(st == STATE_EXCLUDED)),
                    "n_bins_hit_after_filter": int(bins_hit.size),
                    "n_usable_bins_hit_after_filter": int(usable_bins_hit.size),
                    "n_unusable_bins_hit_after_filter": int(unusable_bins_hit.size),
                }
            )

        if use_binarization:
            # v3 path: build a BinarizationObservationModel from the filtered
            # MarkerObservations + the filtered reference. The builder
            # internally drops Ambiguous / Excluded markers and precomputes
            # the per-marker linear coefficient matrix the SLSQP solver needs.
            observation = self._binarization_builder.build(
                obs=obs_filtered,
                binarization=self.binarization,
                reference=reference_filtered,
                tier=tier,
                coverage_cap=self.config.coverage.coverage_cap,
            )
        else:
            observation = self.observation_builder.build(
                obs=obs_filtered,
                reference=reference_filtered,
                region_annotations=self.region_annotations,
                tier=tier,
                coverage_cap=self.config.coverage.coverage_cap,
            )
        reference = reference_filtered  # downstream code operates on the filtered reference
        if use_binarization and binarization_debug is not None:
            n_used = int(observation.n_markers)
            binarization_debug["n_markers_used_for_likelihood"] = n_used
            if n_used == 0 and not binarization_debug.get("fail_reason"):
                binarization_debug["fail_reason"] = "no_callable_markers_after_binarization"
                log.warning(
                    "Sample %s: no callable markers after binarization "
                    "(all ambiguous/excluded or unusable bins).",
                    sample.sample_id,
                )

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
            boot = self.bootstrap.estimate(
                observation, reference, self.deconvolver, n_jobs=bootstrap_jobs
            )
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
            boot = self.bootstrap.estimate(
                observation, reference, self.deconvolver, n_jobs=bootstrap_jobs
            )

        # Reliability metrics per cell type.
        #
        # Reliability assignment uses:
        #   (1) p_detection
        #   (2) likelihood_score
        #   (3) q_likelihood (BH-adjusted LRT p-value vs ablated-null)
        #       with raw p_likelihood retained for transparency.
        K = reference.n_cell_types
        likelihood_score = np.full(K + 1, np.nan, dtype=np.float64)
        p_likelihood = np.full(K + 1, np.nan, dtype=np.float64)
        q_likelihood = np.full(K + 1, np.nan, dtype=np.float64)
        p_detection = np.zeros(K + 1, dtype=np.float64)
        for j in range(K):
            p_detection[j] = compute_p_detection(
                boot.proportions_samples[:, j],
                noise_floor=self.config.uncertainty.noise_floor,
            )
            lik_j, plik_j = compute_fit_metrics(
                w_hat=w_hat,
                reference_methylation=reference.methylation,
                observation=observation,
                cell_type_index=j,
                top_n=50,
                binarizer=self.binarization if use_binarization else None,
            )
            likelihood_score[j] = lik_j
            p_likelihood[j] = plik_j
        # Unknown component detection + fit metrics (ablation of unknown)
        p_detection[K] = compute_p_detection(
            boot.proportions_samples[:, K],
            noise_floor=self.config.uncertainty.noise_floor,
        )
        lik_u, plik_u = compute_unknown_fit_metrics(
            w_hat=w_hat,
            reference_methylation=reference.methylation,
            observation=observation,
        )
        likelihood_score[K] = lik_u
        p_likelihood[K] = plik_u
        # Multiple-testing correction across known cell types only.
        q_likelihood[:K] = _benjamini_hochberg(p_likelihood[:K])

        reliability = np.empty(K + 1, dtype=object)
        for j in range(K):
            reliability[j] = assign_reliability(
                p_detection=p_detection[j],
                likelihood_score=likelihood_score[j],
                p_likelihood=p_likelihood[j],
                q_likelihood=q_likelihood[j],
            )
        reliability[K] = assign_reliability(
            p_detection=p_detection[K],
            likelihood_score=likelihood_score[K],
            p_likelihood=p_likelihood[K],
            q_likelihood=q_likelihood[K],
        )

        # n_markers per cell type — count valid markers contributing to that
        # cell type. For WGBS / v2 paths "valid" means n > 0; for v3
        # binarization paths it means the marker survived binarization (i.e.
        # was called U or M and is in a usable bin).
        if use_binarization:
            valid_n = int(observation.n_markers)
        else:
            valid_n = int(np.sum(observation.n > 0))
        n_markers = np.full(K, valid_n, dtype=np.int32)

        # Enriched output fields (Gap 7). For the v3 binarization path
        # "mean_dispersion" doesn't exist as such — the analogous quantity
        # is the per-bin error rate. We report the mean ε across the cell
        # type's discriminative markers as a stand-in so the qc_summary
        # column stays populated.
        if use_binarization:
            mean_dispersion = _per_celltype_mean_error_rate(
                observation=observation,
                reference_methylation=reference.methylation,
                binarizer=self.binarization,
                top_n=50,
            )
        else:
            mean_dispersion = _per_celltype_mean_dispersion(
                observation=observation,
                reference_methylation=reference.methylation,
                top_n=50,
            )
        # Mean coverage = effective depth of coverage in marker regions
        # (Σ reads * fragment_length / Σ marker_widths). Same scale as the
        # tier classification thresholds, so qc_summary.tsv's coverage_tier
        # and mean_coverage columns are always self-consistent. Both paths
        # compute this from obs_filtered which always has k/n.
        mean_coverage = effective_coverage_in_markers(obs_filtered)

        # pct_imputed: fraction of markers where the pre-imputation n was below
        # the per-marker tier threshold but the post-imputation n is at or
        # above it. Always computed from obs_filtered (the raw counts) so the
        # binarization path produces the same value.
        if pre_impute_n_filtered is not None and pre_impute_n_filtered.size:
            post = np.asarray(obs_filtered.n, dtype=np.int64)
            # Recompute the per-marker threshold vector on the filtered subset
            post_min_reads = per_marker_min_reads_vector(
                np.asarray(obs_filtered.n, dtype=np.int64), tier
            )
            was_low = pre_impute_n_filtered < post_min_reads
            was_lifted = was_low & (post >= post_min_reads)
            pct_imputed = float(np.mean(was_lifted))
        else:
            pct_imputed = 0.0

        # Per-marker residual (mu_obs - mu_predicted) for NMF later. Computed
        # from obs_filtered.k / obs_filtered.n so the residual array is
        # the same shape regardless of which observation model is in use.
        with np.errstate(invalid="ignore", divide="ignore"):
            mu_obs = np.where(
                obs_filtered.n > 0,
                obs_filtered.k / np.maximum(obs_filtered.n, 1),
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
            p_goodness=None,
            p_detection=p_detection,
            likelihood_score=likelihood_score,
            p_likelihood=p_likelihood,
            q_likelihood=q_likelihood,
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
            binarization_flag=binarization_flag,
            hemolysis_flag=hemolysis_flag,
            residuals=residuals,
            marker_chrom=np.asarray(obs_filtered.chrom, dtype=object),
            marker_start=np.asarray(obs_filtered.start, dtype=np.int64),
            marker_end=np.asarray(obs_filtered.end, dtype=np.int64),
            binarization_debug=binarization_debug,
        )
        result.qc_flags = compute_qc_flags(
            result=result,
            observation=observation,
            qc_config=self.config.qc,
            binarization_flag=binarization_flag,
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
        obs: MarkerObservations,
        reference: ReferencePanel,
        tier: CoverageTier,
        binarization_flag: str | None,
        pre_impute_n: np.ndarray | None,
        out_dir: Path,
        binarization_debug: dict | None = None,
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
            p_goodness=None,
            p_detection=np.zeros(K + 1, dtype=np.float64),
            likelihood_score=np.full(K + 1, np.nan, dtype=np.float64),
            p_likelihood=np.full(K + 1, np.nan, dtype=np.float64),
            q_likelihood=np.full(K + 1, np.nan, dtype=np.float64),
            reliability=reliability,
            n_markers=np.zeros(K, dtype=np.int32),
            coverage_tier=tier,
            qc_flags=["NO_MARKERS_PASS_TIER_FILTER"],
            mean_dispersion=np.zeros(K, dtype=np.float64),
            # Same effective-coverage formula used everywhere else, so the
            # NO_MARKERS_PASS_TIER_FILTER row in qc_summary still has a
            # comparable mean_coverage value.
            mean_coverage=effective_coverage_in_markers(obs),
            n_markers_used=0,
            pct_imputed=0.0,
            binarization_flag=binarization_flag,
            hemolysis_flag=None,
            overall_qc="FAIL",
            binarization_debug=binarization_debug,
        )
        write_per_sample_too(result, out_dir / f"{sample.sample_id}.too.tsv")
        return result



# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """BH-adjust p-values in [0, 1], preserving NaN entries."""
    pvals = np.asarray(pvals, dtype=np.float64)
    out = np.full(pvals.shape, np.nan, dtype=np.float64)
    finite = np.isfinite(pvals)
    if not np.any(finite):
        return out
    p = np.clip(pvals[finite], 0.0, 1.0)
    m = p.size
    order = np.argsort(p)
    ranked = p[order]
    ranks = np.arange(1, m + 1, dtype=np.float64)
    adj_ranked = ranked * (m / ranks)
    adj_ranked = np.minimum.accumulate(adj_ranked[::-1])[::-1]
    adj_ranked = np.clip(adj_ranked, 0.0, 1.0)
    adj = np.empty_like(adj_ranked)
    adj[order] = adj_ranked
    out[finite] = adj
    return out


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
    """Return a new MarkerObservations keeping only the masked markers.

    Preserves the v3 binarization fields (``called_state`` and
    ``context_bin``) when present so the per-marker tier filter doesn't
    silently drop the FinaleMe binarization output.
    """
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
        called_state=(
            np.asarray(obs.called_state, dtype=np.uint8)[mask]
            if obs.called_state is not None
            else None
        ),
        context_bin=(
            np.asarray(obs.context_bin, dtype=np.int32)[mask]
            if obs.context_bin is not None
            else None
        ),
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


def _per_celltype_mean_error_rate(
    observation,
    reference_methylation: np.ndarray,
    binarizer,
    top_n: int = 50,
) -> np.ndarray:
    """Per-cell-type mean ε over the cell type's top discriminative markers.

    Analog of ``_per_celltype_mean_dispersion`` for the v3 binarization
    path. Reports the average bin error rate ``mean(ε_U_b, ε_M_b)`` over
    the top-N discriminative markers for each cell type — small values
    mean low classification noise at this cell type's most informative
    markers, larger values mean noisy bins. Reported in the
    ``mean_dispersion`` column of the per-sample TSV so the column has a
    consistent meaning regardless of mode (smaller = better).
    """
    R = np.asarray(reference_methylation, dtype=np.float64)
    M_orig, K = R.shape
    out = np.zeros(K, dtype=np.float64)

    # Soft-binarized reference for the discrimination score (matches
    # what compute_p_goodness uses).
    R_binary = R.copy()
    R_binary[R < 0.2] = 0.0
    R_binary[R > 0.8] = 1.0

    # observation.valid_mask maps the original M markers to the filtered
    # observation arrays. We score discriminativeness over the original
    # markers but only consider those that survived binarization.
    if observation.n_markers == 0:
        return out
    valid_mask = observation.valid_mask
    valid_to_filtered = np.cumsum(valid_mask) - 1

    bin_idx_per_marker = np.asarray(observation.context_bin, dtype=np.int64)
    eps_U_arr = np.asarray(binarizer.eps_U, dtype=np.float64)
    eps_M_arr = np.asarray(binarizer.eps_M, dtype=np.float64)
    mean_eps_per_marker = 0.5 * (
        eps_U_arr[bin_idx_per_marker] + eps_M_arr[bin_idx_per_marker]
    )

    for j in range(K):
        target = R_binary[:, j]
        others = np.delete(R_binary, j, axis=1)
        bg_mean = np.mean(others, axis=1)
        score = np.where(valid_mask, np.abs(target - bg_mean), -np.inf)
        if np.all(np.isneginf(score)):
            out[j] = float(np.mean(mean_eps_per_marker)) if mean_eps_per_marker.size else 0.0
            continue
        take = min(top_n, int(np.sum(valid_mask)))
        top_idx_orig = np.argpartition(-score, take - 1)[:take]
        # Translate to filtered indices
        top_idx_filtered = valid_to_filtered[top_idx_orig]
        out[j] = float(np.mean(mean_eps_per_marker[top_idx_filtered]))
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
    threads: int = 1,
) -> tuple[ReferencePanel, MarkerRegions, dict | None]:
    """Resolve reference panel and marker regions from CLI / config.

    ``threads`` is forwarded to the binary .beta loader so the per-sample
    .beta parsing in ``load_beta_list`` runs in parallel for large reference
    panels.
    """
    t_build_start = perf_counter()
    cpg_index = None
    cpg_index_path = explicit_cpg_index or config.reference.cpg_index
    if cpg_index_path is not None:
        t0 = perf_counter()
        log.info("Reference build: loading CpG index from %s", cpg_index_path)
        cpg_index = load_cpg_index(cpg_index_path)
        chr_positions = cpg_index.get("chr_positions", {})
        total_sites = int(cpg_index.get("total_sites", 0))
        log.info(
            "Reference build: loaded CpG index (%d chromosomes, %d CpGs) in %.2fs",
            len(chr_positions),
            total_sites,
            perf_counter() - t0,
        )

    # Marker regions
    marker_regions_path = explicit_markers or config.markers.marker_regions
    marker_format = explicit_marker_format or config.markers.marker_format

    marker_regions: MarkerRegions
    if marker_regions_path:
        t0 = perf_counter()
        log.info(
            "Reference build: loading marker regions from %s (format=%s)",
            marker_regions_path,
            marker_format,
        )
        marker_regions = MarkerRegionsLoader.load(marker_regions_path, marker_format=marker_format)
        log.info(
            "Reference build: loaded %d marker regions in %.2fs",
            marker_regions.n_markers,
            perf_counter() - t0,
        )
    elif explicit_reference or config.reference.reference_panel:
        # Coordinates can come from the matrix reference panel itself
        path = explicit_reference or config.reference.reference_panel
        t0 = perf_counter()
        log.info(
            "Reference build: deriving marker regions from reference panel coordinates in %s",
            path,
        )
        ref_for_coords = ReferencePanelLoader.load_matrix(path)
        marker_regions = ref_for_coords.to_marker_regions()
        log.info(
            "Reference build: derived %d marker regions from %s in %.2fs",
            marker_regions.n_markers,
            path,
            perf_counter() - t0,
        )
    else:
        raise ValueError(
            "No marker regions: provide --marker-regions or --reference-panel"
        )

    # Reference panel
    if explicit_reference or config.reference.reference_panel:
        reference_panel_path = explicit_reference or config.reference.reference_panel
        t0 = perf_counter()
        log.info(
            "Reference build: loading reference panel matrix from %s",
            reference_panel_path,
        )
        reference = ReferencePanelLoader.load_matrix(
            reference_panel_path
        )
        log.info(
            "Reference build: loaded matrix reference (%d markers x %d cell types) in %.2fs",
            reference.n_markers,
            reference.n_cell_types,
            perf_counter() - t0,
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
        t0 = perf_counter()
        log.info(
            "Reference build: building reference from .beta list (%s) + groups (%s) with %d thread(s)",
            ref_betas_arg,
            ref_groups,
            threads,
        )
        reference = ReferencePanelLoader.load_beta_list(
            ref_betas_arg=ref_betas_arg,
            groups_file=ref_groups,
            cpg_index_path=cpg_index_path,
            marker_regions=marker_regions,
            threads=threads,
        )
        log.info(
            "Reference build: built beta-list reference (%d markers x %d cell types) in %.2fs",
            reference.n_markers,
            reference.n_cell_types,
            perf_counter() - t0,
        )
    else:
        raise ValueError("No reference panel: provide --reference-panel or --ref-betas")

    log.info("Reference build: completed in %.2fs", perf_counter() - t_build_start)
    return reference, marker_regions, cpg_index


def load_optional_binarization(
    config: TOOConfig,
    explicit_path: str | None,
    use_default: bool = True,
):
    """Load v3 binarization parameters JSON if available, else default, else None.

    Reads from ``config.binarization.binarization_file`` when no explicit
    path is given. Falls back to the shipped default placeholder if
    ``use_default`` is True. Returns ``None`` when neither is available.
    """
    from finaleme_too.preprocessing.binarization import (
        BinarizationParams,
        load_default_binarization,
    )

    path = explicit_path or config.binarization.binarization_file
    if path is not None:
        return BinarizationParams.load(path)
    if use_default and config.binarization.use_default:
        try:
            return load_default_binarization()
        except Exception as exc:  # pragma: no cover
            log.warning("Could not load default binarization: %s", exc)
            return None
    return None


def load_optional_region_annotations(
    config: TOOConfig, explicit_path: str | None
) -> pd.DataFrame | None:
    """Load a region_annotation TSV and normalize chromosome naming.

    The output TSV from ``make-region-annotation`` / ``compute_region_annotation``
    strips any ``chr`` prefix (so rows look like ``1`` / ``2`` / ``X``),
    but the downstream consumer (``apply_binarization``) joins on the raw
    ``obs.chrom`` which may or may not carry the prefix depending on how
    the marker regions file was written. Without normalization the join
    silently misses and every marker falls back to the ``open_sea`` bin,
    producing wrong context-bin assignments at inference time.

    Normalize both sides by stripping the ``chr`` prefix here; the
    inference-time apply_binarization applies the same strip to its
    lookup keys for defense-in-depth.
    """
    path = explicit_path or config.markers.region_annotation
    if path is None:
        log.info("Region annotations: none provided; skipping.")
        return None
    from finaleme_too.preprocessing._matched_sample_sheet import _normalize_chrom

    t0 = perf_counter()
    log.info("Region annotations: loading from %s", path)
    df = pd.read_csv(path, sep="\t", comment="#")
    normalized = False
    if "chrom" in df.columns:
        df = _normalize_chrom(df)
        normalized = True
    log.info(
        "Region annotations: loaded %d rows x %d columns in %.2fs%s",
        int(df.shape[0]),
        int(df.shape[1]),
        perf_counter() - t0,
        " (normalized chromosome labels)" if normalized else "",
    )
    return df


__all__ = [
    "CohortResult",
    "TOOPipeline",
    "build_reference_and_markers",
    "load_optional_binarization",
    "load_optional_region_annotations",
]
