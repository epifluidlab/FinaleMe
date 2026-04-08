"""Click-based CLI for finaleme-too.

Phase A scope: ``finaleme-too run`` only. Phase C adds ``train-calibration``.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import click

from finaleme_too._version import __version__
from finaleme_too.config import (
    MeasurementMode,
    SolverMethod,
    TestMethod,
    TOOConfig,
)
from finaleme_too.io.sample_sheet import SampleSheet
from finaleme_too.pipeline import (
    TOOPipeline,
    build_reference_and_markers,
    load_optional_calibration,
    load_optional_region_annotations,
)


@click.group()
@click.version_option(__version__, prog_name="finaleme-too")
def main() -> None:
    """Tissue-of-origin deconvolution for cfDNA methylation data."""


@main.command("run")
@click.option("--sample-sheet", "sample_sheet_path", required=True, type=click.Path(exists=True),
              help="TSV with required columns: sample_id, methylation_file, mode.")
@click.option("--output-dir", required=True, type=click.Path(),
              help="Output directory (created if it does not exist).")
# Reference panel — choose one of:
@click.option("--reference-panel", default=None, type=click.Path(),
              help="TSV reference panel: chrom start end CellType1 ...")
@click.option("--ref-betas", default=None,
              help="Comma-separated reference .beta files OR .txt file listing one path per line. "
                   "Requires --ref-groups and --cpg-index.")
@click.option("--ref-groups", default=None, type=click.Path(),
              help="CSV file mapping sample names to cell types (columns: name, group).")
@click.option("--cpg-index", default=None, type=click.Path(),
              help="CpG index BED file (e.g. data/CpG_index.hg19.bed.gz).")
# Marker regions
@click.option("--marker-regions", default=None, type=click.Path(),
              help="BED or UXM atlas file specifying marker region coordinates.")
@click.option("--marker-format", default="auto", type=click.Choice(["auto", "bed", "uxm_atlas"]),
              help="Force a specific marker file format.")
# Mode + format
@click.option("--mode", type=click.Choice(["WGBS", "FINALEME"]), default=None,
              help="Measurement mode (overrides per-sample mode in sample sheet).")
@click.option("--input-format", default="auto",
              type=click.Choice(["auto", "finaleme_bed", "bissnp_6plus2", "wgbstools_beta", "custom_bed"]),
              help="Methylation file format (overrides per-sample input_format).")
@click.option("--meth-col", default=None, type=int, help="(custom_bed) 1-indexed methylated count column")
@click.option("--total-col", default=None, type=int, help="(custom_bed) 1-indexed total count column")
# Calibration (P1; ignored in Phase A)
@click.option("--calibration", default=None, type=click.Path(),
              help="(P1) Calibration parameters JSON for FinaleMe mode.")
@click.option("--region-annotation", default=None, type=click.Path(),
              help="(P1) Pre-computed CpG density / region class annotations.")
# Marker selection (P1)
@click.option("--strict-regions", default=None, help='(P1) e.g. "CGI+shore"')
@click.option("--n-markers-per-type", default=None, type=int)
# Model
@click.option("--n-bootstrap", default=200, type=int, help="Bootstrap iterations (default: 200).")
@click.option("--bayesian", is_flag=True, default=False, help="(P3) Use Bayesian deconvolution.")
@click.option("--bayesian-n-samples", default=5000, type=int)
@click.option(
    "--uncertainty-method",
    default="bootstrap",
    type=click.Choice(["bootstrap", "bayesian", "both", "none"]),
    help="Uncertainty source: bootstrap (default), bayesian (MCMC posterior), "
    "both (bootstrap is primary, posterior stored alongside), or none (CIs = point).",
)
# Covariates (P2)
@click.option("--batch-correct", default=None, help="(P2) comma-separated batch covariates.")
@click.option("--adjust-covariates", default=None, help="(P2) comma-separated biological covariates.")
@click.option("--configurable-covariates", default=None,
              help="(P2) comma-separated user-configurable covariates.")
# Testing (P1)
@click.option("--test-method", default="ilr_regression",
              type=click.Choice(["ilr_regression", "bayesian_posterior", "wilcoxon"]))
@click.option("--group-comparison", default=None,
              help="(P1) e.g. all, A:B,C:D, A:rest, omnibus+pairwise")
@click.option("--fdr-method", default="BH")
@click.option("--fdr-alpha", default=0.05, type=float)
# Performance
@click.option("--threads", default=1, type=int)
@click.option("--min-coverage", default=3, type=int)
@click.option("--coverage-cap", default=50, type=int)
@click.option("--config", "config_path", default=None, type=click.Path(),
              help="YAML config file. CLI options override config values.")
@click.option("--verbose", is_flag=True, default=False)
def run_cmd(
    sample_sheet_path: str,
    output_dir: str,
    reference_panel: str | None,
    ref_betas: str | None,
    ref_groups: str | None,
    cpg_index: str | None,
    marker_regions: str | None,
    marker_format: str,
    mode: str | None,
    input_format: str,
    meth_col: int | None,
    total_col: int | None,
    calibration: str | None,
    region_annotation: str | None,
    strict_regions: str | None,
    n_markers_per_type: int | None,
    n_bootstrap: int,
    bayesian: bool,
    bayesian_n_samples: int,
    uncertainty_method: str,
    batch_correct: str | None,
    adjust_covariates: str | None,
    configurable_covariates: str | None,
    test_method: str,
    group_comparison: str | None,
    fdr_method: str,
    fdr_alpha: float,
    threads: int,
    min_coverage: int,
    coverage_cap: int,
    config_path: str | None,
    verbose: bool,
) -> None:
    """Run TOO deconvolution end-to-end on a cohort."""
    _setup_logging(verbose)

    # Build / load config
    if config_path:
        config = TOOConfig.from_yaml(config_path)
    else:
        config = TOOConfig()

    # Apply CLI overrides ONLY for options the user explicitly passed.
    # Click's get_parameter_source returns COMMANDLINE for explicit args
    # and DEFAULT for the click option default. This way YAML values are
    # not silently clobbered by Click defaults.
    ctx = click.get_current_context()

    def _was_provided(name: str) -> bool:
        try:
            from click.core import ParameterSource

            return ctx.get_parameter_source(name) == ParameterSource.COMMANDLINE
        except Exception:  # pragma: no cover
            return False

    if _was_provided("threads"):
        config.threads = threads
    if _was_provided("coverage_cap"):
        config.coverage.coverage_cap = coverage_cap
    if _was_provided("min_coverage"):
        config.coverage.min_reads = min_coverage
    if _was_provided("n_bootstrap"):
        config.uncertainty.n_bootstrap = n_bootstrap
    if _was_provided("bayesian_n_samples"):
        config.uncertainty.bayesian_n_samples = bayesian_n_samples
    if _was_provided("uncertainty_method"):
        config.uncertainty.method = uncertainty_method
    if bayesian:  # is_flag — only True when explicitly set
        config.model.deconvolution = SolverMethod.BAYESIAN
        # --bayesian implies Bayesian uncertainty unless the user picked
        # something different via --uncertainty-method
        if not _was_provided("uncertainty_method"):
            config.uncertainty.method = "bayesian"
    if _was_provided("test_method"):
        config.testing.method = TestMethod(test_method)
    if _was_provided("fdr_alpha"):
        config.testing.fdr_alpha = fdr_alpha
    if _was_provided("fdr_method"):
        config.testing.fdr_method = fdr_method
    # Lists: set them only when the user explicitly passed them
    if batch_correct is not None:
        config.batch_correction.technical_covariates = [
            c.strip() for c in batch_correct.split(",") if c.strip()
        ]
    if adjust_covariates is not None:
        config.covariate_adjustment.biological_covariates = [
            c.strip() for c in adjust_covariates.split(",") if c.strip()
        ]
    if configurable_covariates is not None:
        config.covariate_adjustment.user_configurable = [
            c.strip() for c in configurable_covariates.split(",") if c.strip()
        ]

    sample_sheet = SampleSheet.from_tsv(sample_sheet_path)
    sample_sheet.validate_files_exist()

    # Optional global mode override
    if mode is not None:
        forced_mode = MeasurementMode(mode)
        for s in sample_sheet.samples:
            s.mode = forced_mode  # type: ignore[misc]

    # Optional global input format override:
    #   1. CLI flag wins if the user explicitly provided one
    #   2. Otherwise config.input.format wins (when not "auto")
    #   3. Otherwise the per-sample sheet input_format is kept
    if _was_provided("input_format"):
        effective_input_format = input_format
    else:
        effective_input_format = config.input.format
    if effective_input_format and effective_input_format != "auto":
        for s in sample_sheet.samples:
            s.input_format = effective_input_format  # type: ignore[misc]

    # meth_col / total_col follow the same precedence as input_format
    if _was_provided("meth_col"):
        effective_meth_col = meth_col
    else:
        effective_meth_col = config.input.meth_col
    if effective_meth_col is not None:
        for s in sample_sheet.samples:
            s.meth_col = int(effective_meth_col)  # type: ignore[misc]

    if _was_provided("total_col"):
        effective_total_col = total_col
    else:
        effective_total_col = config.input.total_col
    if effective_total_col is not None:
        for s in sample_sheet.samples:
            s.total_col = int(effective_total_col)  # type: ignore[misc]

    # Apply marker selection / strict regions overrides from CLI
    if strict_regions is not None:
        config.markers.strict_regions = strict_regions
    if n_markers_per_type is not None:
        config.markers.n_per_type = n_markers_per_type
    # --marker-format overrides the config value only when the user
    # explicitly passed it on the command line (otherwise YAML wins).
    if _was_provided("marker_format"):
        config.markers.marker_format = marker_format
    effective_marker_format = config.markers.marker_format

    reference, markers, cpg_idx = build_reference_and_markers(
        config=config,
        explicit_reference=reference_panel,
        explicit_markers=marker_regions,
        explicit_ref_betas=ref_betas,
        explicit_ref_groups=ref_groups,
        explicit_cpg_index=cpg_index,
        explicit_marker_format=effective_marker_format,
    )

    # Phase B: load calibration + region annotations (FinaleMe mode)
    calibration_params = None
    region_annotations = None
    has_finaleme = any(s.mode.value == "FINALEME" for s in sample_sheet.samples)
    if has_finaleme:
        calibration_params = load_optional_calibration(config, calibration)
        region_annotations = load_optional_region_annotations(config, region_annotation)

    # Group comparison: fall back to config default if CLI is omitted
    effective_group_comparison = (
        group_comparison if group_comparison is not None else config.testing.group_comparison
    )

    pipeline = TOOPipeline(
        config=config,
        calibration=calibration_params,
        region_annotations=region_annotations,
        group_comparison_spec=effective_group_comparison,
        cpg_index=cpg_idx,
    )
    pipeline.run(sample_sheet, reference, markers, output_dir, cpg_index=cpg_idx)
    click.echo(f"finaleme-too: wrote outputs to {output_dir}")


@main.command("train-calibration")
@click.option("--matched-wgbs", required=True, type=click.Path(exists=True),
              help="Matched WGBS counts table (TSV). Columns: sample_id, chrom, start, "
                   "end, methylated_count, total_count.")
@click.option("--matched-finaleme", required=True, type=click.Path(exists=True),
              help="Matched FinaleMe counts table (TSV). Same schema as --matched-wgbs.")
@click.option("--region-annotation", default=None, type=click.Path(),
              help="Per-marker CpG density annotation TSV (chrom, start, end, cpg_density). "
                   "If omitted, all markers are placed in bin 0.")
@click.option("--n-bins-candidates", default="4,6,8,10,12,16",
              help="Comma-separated list of bin counts to evaluate via CV.")
@click.option("--output", required=True, type=click.Path(),
              help="Output JSON path for the trained CalibrationParams.")
@click.option("--report", required=True, type=click.Path(),
              help="Output JSON path for the calibration training report.")
@click.option("--threads", default=1, type=int)
@click.option("--verbose", is_flag=True, default=False)
def train_calibration_cmd(
    matched_wgbs: str,
    matched_finaleme: str,
    region_annotation: str | None,
    n_bins_candidates: str,
    output: str,
    report: str,
    threads: int,
    verbose: bool,
) -> None:
    """Train per-bin FinaleMe calibration parameters."""
    _setup_logging(verbose)
    from finaleme_too.preprocessing.calibration import train_calibration

    candidates = [int(x.strip()) for x in n_bins_candidates.split(",") if x.strip()]
    params = train_calibration(
        matched_wgbs=matched_wgbs,
        matched_finaleme=matched_finaleme,
        region_annotation=region_annotation,
        n_bins_candidates=candidates,
        out_params=output,
        out_report=report,
    )
    click.echo(
        f"finaleme-too: trained calibration with B={params.n_bins} → {output} (report: {report})"
    )


def _setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="[%(asctime)s] %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
    )


if __name__ == "__main__":
    main()
