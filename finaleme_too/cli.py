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
from finaleme_too.pipeline import TOOPipeline, build_reference_and_markers


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

    # Apply CLI overrides
    config.threads = threads
    config.coverage.coverage_cap = coverage_cap
    config.coverage.min_reads = min_coverage
    config.uncertainty.n_bootstrap = n_bootstrap
    if bayesian:
        config.model.deconvolution = SolverMethod.BAYESIAN
    if test_method:
        config.testing.method = TestMethod(test_method)

    sample_sheet = SampleSheet.from_tsv(sample_sheet_path)
    sample_sheet.validate_files_exist()

    # Optional global mode override
    if mode is not None:
        forced_mode = MeasurementMode(mode)
        for s in sample_sheet.samples:
            s.mode = forced_mode  # type: ignore[misc]

    # Optional global input format override (set on every sample)
    if input_format != "auto":
        for s in sample_sheet.samples:
            s.input_format = input_format  # type: ignore[misc]
    if meth_col is not None:
        for s in sample_sheet.samples:
            s.meth_col = meth_col  # type: ignore[misc]
    if total_col is not None:
        for s in sample_sheet.samples:
            s.total_col = total_col  # type: ignore[misc]

    reference, markers, cpg_idx = build_reference_and_markers(
        config=config,
        explicit_reference=reference_panel,
        explicit_markers=marker_regions,
        explicit_ref_betas=ref_betas,
        explicit_ref_groups=ref_groups,
        explicit_cpg_index=cpg_index,
        explicit_marker_format=marker_format,
    )

    pipeline = TOOPipeline(config)
    pipeline.run(sample_sheet, reference, markers, output_dir, cpg_index=cpg_idx)
    click.echo(f"finaleme-too: wrote outputs to {output_dir}")


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
