#!/usr/bin/env python3
"""Fragment-level TOO deconvolution with bootstrap CI and statistical testing (design §3.5).

Usage:
    python fragment_too_deconvolution.py \\
        --pat-gz decoded_sample.pat.gz \\
        --reference-panel reference_panel.tsv \\
        --cpg-positions CG_motif.hg19.pos_only.bedgraph \\
        --min-cpgs 3 \\
        --informativeness-threshold 0.2 \\
        --estimate-unknown \\
        --bootstrap-replicates 1000 \\
        --test-method lrt \\
        --fdr-threshold 0.05 \\
        --output-prefix sample_001_too \\
        --n-jobs 8
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
_SCRIPT_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SCRIPT_DIR.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from finaleme_too.core.fragment_likelihood import FragmentLevelDeconvolver
from finaleme_too.core.fragment_uncertainty import (
    FragmentBootstrapCI,
    fragment_bh_correction,
    fragment_lrt_pvalues,
    fragment_permutation_pvalues,
)
from finaleme_too.io.fragment_output_writer import (
    write_fragment_diagnostics,
    write_fragment_proportions,
    write_fragment_responsibilities,
)
from finaleme_too.io.marker_regions import MarkerRegions
from finaleme_too.io.pat_loader import load_fragments_from_pat
from finaleme_too.io.reference_panel import ReferencePanelLoader, load_cpg_index

log = logging.getLogger(__name__)


def run_fragment_deconvolution(
    pat_path: str | Path,
    reference_panel_path: str | Path,
    cpg_index_path: str | Path,
    output_prefix: str,
    min_cpgs: int = 3,
    informativeness_threshold: float = 0.2,
    estimate_unknown: bool = True,
    bootstrap_replicates: int = 1000,
    test_method: str = "lrt",
    permutation_replicates: int = 1000,
    fdr_threshold: float = 0.05,
    n_jobs: int = 1,
    seed: int | None = None,
) -> dict:
    """Run the complete fragment-level TOO deconvolution pipeline.

    Returns diagnostics dict.
    """
    # 1. Load reference panel
    log.info("Loading reference panel from %s", reference_panel_path)
    reference = ReferencePanelLoader.load_matrix(reference_panel_path)
    cell_types = reference.cell_types
    K = reference.n_cell_types
    log.info("Reference panel: %d markers x %d cell types", reference.n_markers, K)

    # 2. Load CpG index
    log.info("Loading CpG index from %s", cpg_index_path)
    cpg_index = load_cpg_index(cpg_index_path)

    # 3. Build marker regions from reference panel coordinates
    marker_regions = reference.to_marker_regions()

    # 4. Load and filter fragments
    log.info("Loading fragments from %s (min_cpgs=%d, info_thresh=%.2f) ...",
             pat_path, min_cpgs, informativeness_threshold)
    fragments = load_fragments_from_pat(
        pat_path,
        marker_regions=marker_regions,
        cpg_index=cpg_index,
        reference_methylation=reference.methylation,
        min_cpgs=min_cpgs,
        informativeness_threshold=informativeness_threshold,
    )
    n_fragments_used = len(fragments)
    log.info("Fragments after filtering: %d", n_fragments_used)

    if n_fragments_used == 0:
        log.warning("No informative fragments found. Writing empty results.")
        K_total = (K + 1) if estimate_unknown else K
        empty_w = np.zeros(K_total)
        if estimate_unknown:
            empty_w[-1] = 1.0
        empty_gamma = np.zeros((0, K_total), dtype=np.float64)
        write_fragment_proportions(
            cell_types if estimate_unknown else cell_types,
            empty_w, None, None, None, None, fdr_threshold,
            f"{output_prefix}_proportions.tsv",
        )
        write_fragment_responsibilities(
            empty_gamma, cell_types,
            f"{output_prefix}_fragment_responsibilities.tsv.gz",
        )
        write_fragment_diagnostics(
            {"n_fragments_total": 0, "n_fragments_used": 0, "em_converged": False},
            f"{output_prefix}_diagnostics.json",
        )
        return {"n_fragments_used": 0}

    # 5. Run EM (full solve)
    log.info("Running fragment-level EM (unknown_component=%s) ...", estimate_unknown)
    deconvolver = FragmentLevelDeconvolver(
        max_iter=200,
        tol=1e-6,
        unknown_profile=0.5,
        include_unknown=estimate_unknown,
    )
    w, gamma, ll, log_p, em_converged = deconvolver.solve_full(fragments, reference.methylation)
    if em_converged:
        log.info("EM converged. Log-likelihood: %.2f", ll)
    else:
        log.warning("EM did not converge within max_iter. Log-likelihood: %.2f", ll)

    # 6. Bootstrap CI
    ci_lower = None
    ci_upper = None
    if bootstrap_replicates > 0:
        log.info("Running bootstrap (%d replicates, %d jobs) ...", bootstrap_replicates, n_jobs)
        bootstrap = FragmentBootstrapCI(
            n_iterations=bootstrap_replicates,
            ci_level=0.95,
            seed=seed,
        )
        boot_result = bootstrap.estimate(log_p, max_iter=200, tol=1e-6, n_jobs=n_jobs)
        ci_lower = boot_result.ci_lower
        ci_upper = boot_result.ci_upper
        log.info("Bootstrap done.")

    # 7. Statistical testing (p-values)
    p_values = None
    if test_method == "lrt":
        log.info("Computing LRT p-values (%d jobs) ...", n_jobs)
        p_values = fragment_lrt_pvalues(log_p, ll, max_iter=200, tol=1e-6, n_jobs=n_jobs)
    elif test_method == "permutation":
        log.info("Computing permutation p-values (%d permutations, %d jobs) ...",
                 permutation_replicates, n_jobs)
        p_values = fragment_permutation_pvalues(
            log_p, w, n_permutations=permutation_replicates,
            max_iter=200, tol=1e-6, seed=seed, n_jobs=n_jobs,
        )

    # 8. BH correction
    q_values = None
    if p_values is not None:
        q_values = fragment_bh_correction(p_values)

    # 9. Write output files
    log.info("Writing output files with prefix '%s' ...", output_prefix)

    write_fragment_proportions(
        cell_types, w, ci_lower, ci_upper, p_values, q_values, fdr_threshold,
        f"{output_prefix}_proportions.tsv",
    )

    write_fragment_responsibilities(
        gamma, cell_types,
        f"{output_prefix}_fragment_responsibilities.tsv.gz",
    )

    diagnostics = {
        "n_fragments_used": n_fragments_used,
        "em_log_likelihood": float(ll),
        "em_converged": em_converged,
        "pi_unknown": float(w[-1]),
        "bootstrap_replicates": bootstrap_replicates,
        "test_method": test_method,
        "fdr_threshold": fdr_threshold,
        "cell_types": cell_types,
        "proportions": w.tolist(),
    }
    write_fragment_diagnostics(diagnostics, f"{output_prefix}_diagnostics.json")

    log.info("Done. Results written to %s_*.tsv/.json", output_prefix)
    return diagnostics


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fragment-level TOO deconvolution (design §3.5).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--pat-gz", required=True, help="Input .pat.gz file from FinaleMe Viterbi decoding")
    parser.add_argument("--reference-panel", required=True, help="Reference panel TSV (from generate_reference_panel.py)")
    parser.add_argument("--cpg-positions", required=True, help="CG_motif bedgraph file with all CpG positions")
    parser.add_argument("--min-cpgs", type=int, default=3, help="Minimum CpGs per fragment (default: 3)")
    parser.add_argument("--informativeness-threshold", type=float, default=0.2,
                        help="Min informativeness score I(f) (default: 0.2)")
    parser.add_argument("--estimate-unknown", action=argparse.BooleanOptionalAction, default=True,
                        help="Include unknown cell type component (default: True; use --no-estimate-unknown to disable)")
    parser.add_argument("--bootstrap-replicates", type=int, default=1000,
                        help="Number of bootstrap replicates for CI (default: 1000; 0 to skip)")
    parser.add_argument("--test-method", choices=["lrt", "permutation", "none"], default="lrt",
                        help="Statistical test method (default: lrt)")
    parser.add_argument("--permutation-replicates", type=int, default=1000,
                        help="Number of permutations if --test-method=permutation (default: 1000)")
    parser.add_argument("--fdr-threshold", type=float, default=0.05,
                        help="FDR threshold for significance (default: 0.05)")
    parser.add_argument("--output-prefix", required=True, help="Output file prefix")
    parser.add_argument("--n-jobs", type=int, default=1, help="Parallel workers (default: 1)")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    run_fragment_deconvolution(
        pat_path=args.pat_gz,
        reference_panel_path=args.reference_panel,
        cpg_index_path=args.cpg_positions,
        output_prefix=args.output_prefix,
        min_cpgs=args.min_cpgs,
        informativeness_threshold=args.informativeness_threshold,
        estimate_unknown=args.estimate_unknown,
        bootstrap_replicates=args.bootstrap_replicates,
        test_method=args.test_method,
        permutation_replicates=args.permutation_replicates,
        fdr_threshold=args.fdr_threshold,
        n_jobs=args.n_jobs,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
