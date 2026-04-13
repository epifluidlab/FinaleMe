#!/usr/bin/env python3
"""Build a per-CpG reference methylation panel from WGBS data (design §3.4).

Supports four input formats via a unified manifest TSV:
  - pat.gz      : wgbstools/FinaleMe fragment patterns
  - .beta       : wgbstools binary (2 bytes per CpG: methy, total)
  - bissnp      : Bis-SNP 6+2 BED (col7=methy_pct, col8=total)
  - bigwig      : paired methylation + coverage bigWig files

Usage:
    python generate_reference_panel.py \\
        --manifest manifest.tsv \\
        --marker-regions cgi_shore_markers.bed \\
        --cpg-positions CG_motif.hg19.pos_only.bedgraph \\
        --min-coverage 5 \\
        --variance-percentile 0.01 \\
        --output reference_panel.tsv
"""

from __future__ import annotations

import argparse
import gzip
import logging
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Add project root to path so finaleme_too is importable when running as script
# ---------------------------------------------------------------------------
_SCRIPT_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SCRIPT_DIR.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from finaleme_too.io.reference_panel import ReferencePanel, ReferencePanelLoader, load_cpg_index
from finaleme_too.io.marker_regions import MarkerRegionsLoader

log = logging.getLogger(__name__)


# ============================================================================
# Manifest parsing
# ============================================================================

def load_manifest(path: str | Path) -> pd.DataFrame:
    """Load a TSV manifest and auto-detect format.

    3-column manifests (name, group, file_path) are used for pat/beta/bissnp.
    4-column manifests (name, group, methy_bw, cov_bw) are used for bigwig.

    The design doc uses ``# name  group  file_path`` as the header (the
    leading ``#`` is a comment marker in the spec, not part of the column
    name). We strip the ``#`` prefix from the header so that the column
    names parse correctly regardless of whether the ``#`` is present.
    """
    df = pd.read_csv(path, sep="\t", comment=None)
    # Strip leading '#' and whitespace from column names
    df.columns = [c.lstrip("# ").strip() for c in df.columns]
    # Drop any rows that are pure comment lines (first cell starts with #)
    first_col = df.columns[0]
    comment_mask = df[first_col].astype(str).str.startswith("#")
    df = df[~comment_mask].reset_index(drop=True)
    return df


def detect_input_format(manifest_df: pd.DataFrame, user_format: str) -> str:
    """Determine input format from manifest columns or user override."""
    if user_format != "auto":
        return user_format

    if manifest_df.empty:
        raise ValueError("Manifest is empty (no data rows after header)")

    # 4-column with methy_bw => bigwig
    if "methy_bw" in manifest_df.columns and "cov_bw" in manifest_df.columns:
        return "bigwig"

    # 3-column: detect from file extension of first entry
    if "file_path" not in manifest_df.columns:
        raise ValueError(
            f"Manifest must have 'file_path' column (or 'methy_bw'+'cov_bw' for bigwig). "
            f"Found columns: {list(manifest_df.columns)}"
        )

    first_file = str(manifest_df["file_path"].iloc[0]).lower()
    if first_file.endswith(".pat.gz") or first_file.endswith(".pat"):
        return "pat"
    elif first_file.endswith(".beta") or first_file.endswith(".lbeta"):
        return "beta"
    elif "6plus2" in first_file or first_file.endswith(".bed") or first_file.endswith(".bed.gz"):
        return "bissnp"
    else:
        raise ValueError(
            f"Cannot auto-detect format from '{first_file}'. "
            "Use --input-format to specify: pat, beta, bissnp, bigwig"
        )


# ============================================================================
# Per-format aggregation to CpG-level counts
# ============================================================================

def _aggregate_pat_files(
    file_paths: list[Path],
    cpg_index: dict,
    total_sites: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Aggregate pat.gz files into per-CpG (methylated, total) count arrays.

    Returns (methy_counts, total_counts) each of shape (total_sites,).
    """
    methy = np.zeros(total_sites, dtype=np.int64)
    total = np.zeros(total_sites, dtype=np.int64)

    for pat_path in file_paths:
        opener = gzip.open if str(pat_path).endswith(".gz") else open
        with opener(pat_path, "rt") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                try:
                    start_cpg = int(parts[1])  # 1-based global CpG index
                    pattern = parts[2]
                    count = int(parts[3]) if len(parts) >= 4 else 1
                except (ValueError, IndexError):
                    continue
                for offset, ch in enumerate(pattern):
                    if ch not in ("C", "T"):
                        continue
                    cpg_idx = start_cpg + offset - 1  # convert to 0-based
                    if 0 <= cpg_idx < total_sites:
                        total[cpg_idx] += count
                        if ch == "C":
                            methy[cpg_idx] += count

    return methy, total


def _aggregate_beta_files(
    file_paths: list[Path],
    total_sites: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Aggregate .beta binary files into per-CpG (methylated, total) count arrays."""
    methy = np.zeros(total_sites, dtype=np.int64)
    total = np.zeros(total_sites, dtype=np.int64)

    for beta_path in file_paths:
        with open(beta_path, "rb") as fh:
            data = np.frombuffer(fh.read(), dtype=np.uint8)
        if data.size % 2 != 0:
            log.warning("Beta file size not a multiple of 2: %s (skipping)", beta_path)
            continue
        per_cpg = data.reshape((-1, 2)).astype(np.int64)
        n = min(per_cpg.shape[0], total_sites)
        methy[:n] += per_cpg[:n, 0]
        total[:n] += per_cpg[:n, 1]

    return methy, total


def _aggregate_bissnp_files(
    file_paths: list[Path],
    cpg_index: dict,
) -> tuple[np.ndarray, np.ndarray]:
    """Aggregate Bis-SNP 6+2 BED files into per-CpG counts.

    Format: chr start end name score strand methy_pct total_count
    k = round(methy_pct / 100 * total_count)
    """
    chr_positions = cpg_index["chr_positions"]
    chr_offsets = cpg_index["chr_offsets"]
    total_sites = cpg_index["total_sites"]

    methy = np.zeros(total_sites, dtype=np.int64)
    total = np.zeros(total_sites, dtype=np.int64)

    for bed_path in file_paths:
        opener = gzip.open if str(bed_path).endswith(".gz") else open
        with opener(bed_path, "rt") as fh:
            df = pd.read_csv(
                fh, sep="\t", comment="#", header=None,
                usecols=[0, 1, 2, 6, 7],
                names=["chrom", "start", "end", "methy_pct", "total_count"],
                dtype={"chrom": str, "start": np.int64, "end": np.int64,
                       "methy_pct": np.float64, "total_count": np.float64},
            )
        for _, row in df.iterrows():
            chrom = str(row["chrom"])
            positions = chr_positions.get(chrom)
            offset = chr_offsets.get(chrom)
            if positions is None or offset is None:
                continue
            # Find the CpG index for this position
            pos = int(row["start"])
            idx = int(np.searchsorted(positions, pos, side="left"))
            if idx < len(positions) and positions[idx] == pos:
                global_idx = offset + idx
                if 0 <= global_idx < total_sites:
                    tc = int(row["total_count"])
                    mk = int(round(row["methy_pct"] / 100.0 * tc))
                    methy[global_idx] += mk
                    total[global_idx] += tc

    return methy, total


def _aggregate_bigwig_files(
    methy_bw_paths: list[Path],
    cov_bw_paths: list[Path],
    cpg_index: dict,
) -> tuple[np.ndarray, np.ndarray]:
    """Aggregate paired bigWig files into per-CpG counts."""
    try:
        import pyBigWig
    except ImportError:
        raise ImportError(
            "pyBigWig is required for bigwig input format. "
            "Install with: pip install pyBigWig"
        )

    chr_positions = cpg_index["chr_positions"]
    chr_offsets = cpg_index["chr_offsets"]
    total_sites = cpg_index["total_sites"]

    methy = np.zeros(total_sites, dtype=np.float64)
    total = np.zeros(total_sites, dtype=np.float64)

    for methy_bw_path, cov_bw_path in zip(methy_bw_paths, cov_bw_paths):
        bw_methy = pyBigWig.open(str(methy_bw_path))
        bw_cov = pyBigWig.open(str(cov_bw_path))

        for chrom, positions in chr_positions.items():
            offset = chr_offsets[chrom]
            for i, pos in enumerate(positions):
                global_idx = offset + i
                try:
                    mv = bw_methy.values(chrom, int(pos), int(pos) + 1)
                    cv = bw_cov.values(chrom, int(pos), int(pos) + 1)
                    if mv and cv and mv[0] is not None and cv[0] is not None:
                        methy[global_idx] += mv[0]
                        total[global_idx] += cv[0]
                except RuntimeError:
                    continue

        bw_methy.close()
        bw_cov.close()

    return methy.astype(np.int64), total.astype(np.int64)


# ============================================================================
# Aggregation dispatcher
# ============================================================================

def aggregate_to_cpg_matrix(
    manifest_df: pd.DataFrame,
    cpg_index: dict,
    input_format: str,
    min_coverage: int = 5,
    n_jobs: int = 1,
) -> tuple[np.ndarray, np.ndarray, list[str], np.ndarray, np.ndarray]:
    """Aggregate all input files into a CpG x cell_type matrix.

    Returns
    -------
    methylation : ndarray, shape (n_cpgs_used, K) — beta values
    coverage : ndarray, shape (n_cpgs_used, K) — total read counts
    cell_types : list[str]
    chrom_arr : ndarray, shape (n_cpgs_used,) — chromosome per CpG
    position_arr : ndarray, shape (n_cpgs_used,) — position per CpG
    """
    total_sites = cpg_index["total_sites"]
    chr_positions = cpg_index["chr_positions"]
    chr_offsets = cpg_index["chr_offsets"]

    # Group by cell type
    groups = sorted(manifest_df["group"].unique())
    cell_types = groups

    methy_matrix = np.zeros((total_sites, len(cell_types)), dtype=np.int64)
    total_matrix = np.zeros((total_sites, len(cell_types)), dtype=np.int64)

    for ci, ct in enumerate(cell_types):
        ct_rows = manifest_df[manifest_df["group"] == ct]

        if input_format == "pat":
            paths = [Path(p) for p in ct_rows["file_path"]]
            mk, tt = _aggregate_pat_files(paths, cpg_index, total_sites)
        elif input_format == "beta":
            paths = [Path(p) for p in ct_rows["file_path"]]
            mk, tt = _aggregate_beta_files(paths, total_sites)
        elif input_format == "bissnp":
            paths = [Path(p) for p in ct_rows["file_path"]]
            mk, tt = _aggregate_bissnp_files(paths, cpg_index)
        elif input_format == "bigwig":
            methy_paths = [Path(p) for p in ct_rows["methy_bw"]]
            cov_paths = [Path(p) for p in ct_rows["cov_bw"]]
            mk, tt = _aggregate_bigwig_files(methy_paths, cov_paths, cpg_index)
        else:
            raise ValueError(f"Unknown input format: {input_format}")

        methy_matrix[:, ci] = mk
        total_matrix[:, ci] = tt
        log.info("Cell type '%s': aggregated %d files", ct, len(ct_rows))

    # Compute beta values; CpGs below min_coverage in ANY cell type get NaN
    with np.errstate(invalid="ignore", divide="ignore"):
        beta = np.where(total_matrix > 0, methy_matrix / total_matrix, np.nan)

    low_cov_mask = np.any(total_matrix < min_coverage, axis=1)
    beta[low_cov_mask] = np.nan

    # Filter to CpGs with valid data in all cell types
    valid_mask = ~np.any(np.isnan(beta), axis=1)
    valid_indices = np.where(valid_mask)[0]

    # Build chrom/position arrays for valid CpGs
    chrom_arr = np.empty(len(valid_indices), dtype=object)
    position_arr = np.empty(len(valid_indices), dtype=np.int64)
    for chrom, positions in chr_positions.items():
        offset = chr_offsets[chrom]
        for i, pos in enumerate(positions):
            global_idx = offset + i
            local = np.searchsorted(valid_indices, global_idx)
            if local < len(valid_indices) and valid_indices[local] == global_idx:
                chrom_arr[local] = chrom
                position_arr[local] = pos

    return (
        beta[valid_indices],
        total_matrix[valid_indices],
        cell_types,
        chrom_arr,
        position_arr,
    )


# ============================================================================
# Filtering and panel construction
# ============================================================================

def filter_to_marker_regions(
    methylation: np.ndarray,
    coverage: np.ndarray,
    chrom_arr: np.ndarray,
    position_arr: np.ndarray,
    marker_regions,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Restrict CpG-level matrix to CpGs within marker regions."""
    keep = np.zeros(len(chrom_arr), dtype=bool)

    for mi in range(marker_regions.n_markers):
        m_chrom = str(marker_regions.chrom[mi])
        m_start = int(marker_regions.start[mi])
        m_end = int(marker_regions.end[mi])

        chrom_match = chrom_arr == m_chrom
        pos_match = (position_arr >= m_start) & (position_arr < m_end)
        keep |= (chrom_match & pos_match)

    return methylation[keep], coverage[keep], chrom_arr[keep], position_arr[keep]


def apply_variance_filter(
    methylation: np.ndarray,
    variance_percentile: float = 0.01,
) -> np.ndarray:
    """Return boolean mask of CpGs whose cross-cell-type variance exceeds threshold."""
    variances = np.nanvar(methylation, axis=1)
    threshold = np.nanquantile(variances, 1.0 - variance_percentile)
    return variances >= threshold


def build_reference_panel(
    methylation: np.ndarray,
    coverage: np.ndarray,
    cell_types: list[str],
    chrom_arr: np.ndarray,
    position_arr: np.ndarray,
    clip_range: tuple[float, float] = (0.01, 0.99),
) -> ReferencePanel:
    """Assemble final ReferencePanel from filtered matrices."""
    clipped = np.clip(methylation, clip_range[0], clip_range[1])

    # Build start/end arrays (end = position + 2 for a CpG dinucleotide)
    start = position_arr.copy()
    end = position_arr + 2

    return ReferencePanel(
        chrom=chrom_arr,
        start=start,
        end=end,
        cell_types=list(cell_types),
        methylation=clipped.astype(np.float32),
        coverage=coverage.astype(np.int32) if coverage is not None else None,
    )


# ============================================================================
# CLI
# ============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a per-CpG reference panel from WGBS data (design §3.4).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--manifest", required=True,
                        help="TSV manifest file (name, group, file_path) or (name, group, methy_bw, cov_bw)")
    parser.add_argument("--input-format", default="auto",
                        choices=["auto", "pat", "beta", "bissnp", "bigwig"],
                        help="Input file format (default: auto-detect)")
    parser.add_argument("--marker-regions", required=True,
                        help="BED file with marker regions")
    parser.add_argument("--cpg-positions", required=True,
                        help="CG_motif bedgraph file with all CpG positions")
    parser.add_argument("--min-coverage", type=int, default=5,
                        help="Minimum coverage per CpG per cell type (default: 5)")
    parser.add_argument("--variance-percentile", type=float, default=0.01,
                        help="Retain top X%% most variable CpGs (default: 0.01 = top 1%%)")
    parser.add_argument("--output", required=True,
                        help="Output reference panel TSV path")
    parser.add_argument("--threads", type=int, default=1,
                        help="Number of parallel threads (default: 1)")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    log.info("Loading manifest from %s", args.manifest)
    manifest_df = load_manifest(args.manifest)
    input_format = detect_input_format(manifest_df, args.input_format)
    log.info("Detected input format: %s", input_format)

    log.info("Loading CpG index from %s", args.cpg_positions)
    cpg_index = load_cpg_index(args.cpg_positions)
    log.info("Total CpG sites: %d", cpg_index["total_sites"])

    log.info("Loading marker regions from %s", args.marker_regions)
    marker_regions = MarkerRegionsLoader.load(args.marker_regions)
    log.info("Marker regions: %d", marker_regions.n_markers)

    log.info("Aggregating per-CpG methylation from %d files (%d threads) ...",
             len(manifest_df), args.threads)
    methylation, coverage, cell_types, chrom_arr, position_arr = aggregate_to_cpg_matrix(
        manifest_df, cpg_index, input_format, min_coverage=args.min_coverage,
        n_jobs=args.threads,
    )
    log.info("CpGs after min-coverage filter: %d", methylation.shape[0])

    log.info("Filtering to marker regions ...")
    methylation, coverage, chrom_arr, position_arr = filter_to_marker_regions(
        methylation, coverage, chrom_arr, position_arr, marker_regions
    )
    log.info("CpGs within marker regions: %d", methylation.shape[0])

    if args.variance_percentile > 0:
        log.info("Applying variance filter (top %.1f%%) ...", args.variance_percentile * 100)
        var_mask = apply_variance_filter(methylation, args.variance_percentile)
        methylation = methylation[var_mask]
        coverage = coverage[var_mask]
        chrom_arr = chrom_arr[var_mask]
        position_arr = position_arr[var_mask]
        log.info("CpGs after variance filter: %d", methylation.shape[0])

    log.info("Building reference panel with %d CpGs x %d cell types ...",
             methylation.shape[0], len(cell_types))
    panel = build_reference_panel(
        methylation, coverage, cell_types, chrom_arr, position_arr
    )

    log.info("Writing reference panel to %s", args.output)
    ReferencePanelLoader.write_matrix(panel, args.output)
    log.info("Done. Reference panel: %d CpGs x %d cell types", panel.n_markers, panel.n_cell_types)


if __name__ == "__main__":
    main()
