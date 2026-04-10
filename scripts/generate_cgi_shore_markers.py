#!/usr/bin/env python3
"""
Generate CGI+Shore restricted tissue-specific methylation markers and a
pre-aggregated reference panel for cfDNA tissue-of-origin deconvolution.

This script creates tissue-specific methylation markers restricted to CpG Island
(CGI) and CGI shore regions (±2kb around CGI), then builds a reference panel TSV
that ``finaleme-too run --reference-panel`` and ``BetaValueDeconvolution`` can
load directly without re-parsing every .beta file at runtime.

Pipeline stages:
  Stage 0: Download reference WGBS data (optional, if not user-provided)
  Stage 1: Generate CGI+shore BED regions
  Stage 2: Generate candidate blocks within CGI+shore
  Stage 3: Find tissue-specific markers with adaptive thresholds
  Stage 4: Build count-based panel outputs:
           - reference panel TSV (chrom/start/end + per-cell-type k/n)
           - UXM-atlas-compatible TSV (same marker rows as Stage 3 + per-cell-type k/n)
  Stage 5: Generate validation report

The Stage 4 reference-panel output is a TSV with format:

    chrom  start  end  CellType1     CellType2     ...
    chr1   1000   1500 125/250       30/250        ...

Each cell value is ``methylated_count/total_count`` aggregated across all
reference samples in that cell type's group. The downstream loader detects
this format and produces both a methylation matrix (k/n ratios) and a
coverage matrix (n values) for reference uncertainty weighting.

Stage 4 also writes an atlas-style file with the original marker metadata
columns (chr/start/end/startCpG/endCpG/target/name/direction) plus the same
per-cell-type ``k/n`` values, preserving marker row identity/order from
Stage 3.

Usage:
  python generate_cgi_shore_markers.py --genome hg19 --num-markers 250 \\
    --out-dir results/ --wgbstools-path /path/to/wgbs_tools \\
    --uxm-path /path/to/UXM_deconv

Requirements:
  - Python 3.8+ with pandas, numpy, scipy
  - wgbs_tools (https://github.com/nloyfer/wgbs_tools)
  - UXM_deconv (https://github.com/nloyfer/UXM_deconv) [for marker selection only]
  - bedtools (for BED intersection)
  - finaleme_too (optional; provides faster .beta parsing for Stage 4)
"""

import argparse
import glob
import gzip
import os
import os.path as op
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd


###############################################################################
#                                                                             #
#   Utility functions                                                         #
#                                                                             #
###############################################################################

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def run_cmd(cmd, verbose=False, check=True):
    """Run a shell command and return stdout."""
    if verbose:
        eprint(f'[CMD] {cmd}')
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if check and result.returncode != 0:
        eprint(f'Command failed: {cmd}')
        eprint(f'stderr: {result.stderr}')
        raise RuntimeError(f'Command failed with exit code {result.returncode}')
    return result.stdout.strip()


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path


def check_executable(name):
    """Check if an executable is on PATH."""
    for p in os.environ.get('PATH', '').split(':'):
        if os.access(op.join(p, name), os.X_OK):
            return True
    return False


###############################################################################
#                                                                             #
#   Stage 0: Download Reference WGBS Data                                    #
#                                                                             #
###############################################################################

# Known download URLs for Loyfer et al. reference data
# Users should update these if URLs change
REFERENCE_URLS = {
    'hg19': {
        'groups': None,  # Groups file URL (to be filled with actual URL)
        'beta_index': None,  # Index of beta files
        'pat_index': None,  # Index of pat files
    },
    'hg38': {
        'groups': None,
        'beta_index': None,
        'pat_index': None,
    },
}

# UCSC CGI download URLs
CGI_URLS = {
    'hg19': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cpgIslandExt.txt.gz',
    'hg38': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz',
}


def download_file(url, dest, verbose=False):
    """Download a file if it doesn't exist."""
    if op.isfile(dest):
        if verbose:
            eprint(f'  File exists, skipping: {dest}')
        return dest
    if verbose:
        eprint(f'  Downloading {url} -> {dest}')
    ensure_dir(op.dirname(dest))
    urllib.request.urlretrieve(url, dest)
    return dest


def download_cgi_bed(genome, out_dir, verbose=False):
    """Download UCSC CpG Island annotation and convert to BED format."""
    url = CGI_URLS.get(genome)
    if not url:
        raise ValueError(f'No CGI URL for genome {genome}')

    raw_path = op.join(out_dir, f'cpgIslandExt.{genome}.txt.gz')
    bed_path = op.join(out_dir, f'cpgIslandExt.{genome}.bed')

    if op.isfile(bed_path):
        if verbose:
            eprint(f'  CGI BED exists: {bed_path}')
        return bed_path

    download_file(url, raw_path, verbose)

    # Parse UCSC cpgIslandExt format: bin, chrom, chromStart, chromEnd, name, ...
    # We need columns 1 (chrom), 2 (start), 3 (end)
    eprint(f'  Converting to BED: {bed_path}')
    with gzip.open(raw_path, 'rt') as fin, open(bed_path, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            chrom, start, end = parts[1], parts[2], parts[3]
            # Filter to standard chromosomes
            if chrom.startswith('chr') and '_' not in chrom:
                fout.write(f'{chrom}\t{start}\t{end}\n')

    return bed_path


def stage0_download_references(args):
    """Download reference WGBS data if not user-provided."""
    eprint('\n=== Stage 0: Check/Download Reference Data ===')

    ref_dir = ensure_dir(op.join(args.out_dir, 'reference_data'))

    # Check if user provided data
    if args.betas and args.pats and args.groups:
        eprint('  User provided --betas, --pats, --groups. Skipping download.')
        return args.betas, args.pats, args.groups

    eprint('  WARNING: Automatic download of reference WGBS data is not yet implemented.')
    eprint('  Please provide reference data manually using:')
    eprint('    --betas /path/to/ref_betas/*.beta')
    eprint('    --pats /path/to/ref_pats/*.pat.gz')
    eprint('    --groups /path/to/groups.csv')
    eprint('')
    eprint('  Reference data sources:')
    eprint('    - Loyfer et al. Nature 2023 (GSE186458)')
    eprint('    - Download WGBS beta/pat files from GEO or the wgbs_tools atlas repository')
    eprint('    - Groups file: CSV with columns "name" and "group"')
    raise SystemExit(1)


###############################################################################
#                                                                             #
#   Stage 1: Generate CGI+Shore BED Regions                                  #
#                                                                             #
###############################################################################

def merge_intervals(intervals):
    """Merge overlapping intervals. Input: list of (chrom, start, end)."""
    if not intervals:
        return []

    # Sort by chrom, then start
    intervals.sort(key=lambda x: (x[0], x[1]))

    merged = [intervals[0]]
    for chrom, start, end in intervals[1:]:
        prev_chrom, prev_start, prev_end = merged[-1]
        if chrom == prev_chrom and start <= prev_end:
            merged[-1] = (chrom, prev_start, max(end, prev_end))
        else:
            merged.append((chrom, start, end))

    return merged


def stage1_generate_cgi_shore_bed(args):
    """Generate CGI + CGI shore BED file."""
    eprint('\n=== Stage 1: Generate CGI+Shore BED Regions ===')

    out_path = op.join(args.out_dir, f'cgi_plus_shore.{args.genome}.bed')
    if op.isfile(out_path) and not args.force:
        eprint(f'  Output exists: {out_path}')
        return out_path

    # Get CGI BED
    if args.cgi_bed:
        cgi_bed = args.cgi_bed
        eprint(f'  Using provided CGI BED: {cgi_bed}')
    else:
        eprint(f'  Downloading UCSC CpG Island annotation for {args.genome}...')
        cgi_bed = download_cgi_bed(args.genome, args.out_dir, args.verbose)

    # Load CGI regions
    eprint(f'  Loading CGI regions from {cgi_bed}')
    cgi_df = pd.read_csv(cgi_bed, sep='\t', header=None,
                         names=['chrom', 'start', 'end'],
                         usecols=[0, 1, 2])

    # Get chromosome sizes for boundary checking
    chrom_sizes = {}
    if args.chrom_sizes:
        with open(args.chrom_sizes) as f:
            for line in f:
                parts = line.strip().split('\t')
                chrom_sizes[parts[0]] = int(parts[1])

    # Extend by shore_size
    shore = args.shore_size
    eprint(f'  Extending CGI by ±{shore}bp (shore regions)...')
    intervals = []
    for _, row in cgi_df.iterrows():
        chrom = row['chrom']
        start = max(0, int(row['start']) - shore)
        end = int(row['end']) + shore
        if chrom in chrom_sizes:
            end = min(end, chrom_sizes[chrom])
        intervals.append((chrom, start, end))

    # Merge overlapping intervals
    merged = merge_intervals(intervals)
    eprint(f'  {len(cgi_df)} CGI regions -> {len(merged)} merged CGI+shore regions')

    # Calculate total coverage
    total_bp = sum(end - start for _, start, end in merged)
    eprint(f'  Total CGI+shore coverage: {total_bp:,} bp ({total_bp / 1e6:.1f} Mb)')

    # Write output
    with open(out_path, 'w') as f:
        for chrom, start, end in merged:
            f.write(f'{chrom}\t{start}\t{end}\n')

    eprint(f'  Written: {out_path}')
    return out_path


###############################################################################
#                                                                             #
#   Stage 2: Generate Candidate Blocks within CGI+Shore                      #
#                                                                             #
###############################################################################

def stage2_generate_blocks(args, cgi_shore_bed):
    """Generate candidate blocks within CGI+shore regions."""
    eprint('\n=== Stage 2: Generate Candidate Blocks within CGI+Shore ===')

    out_path = op.join(args.out_dir, f'cgi_shore_blocks.{args.genome}.bed')
    if op.isfile(out_path) and not args.force:
        eprint(f'  Output exists: {out_path}')
        return out_path

    wgbstools = op.join(args.wgbstools_path, 'wgbstools')
    if not op.isfile(wgbstools):
        # Try alternative path
        wgbstools = 'wgbstools'
        if not check_executable('wgbstools'):
            raise RuntimeError(
                f'wgbstools not found. Provide --wgbstools-path or add to PATH.')

    if args.blocks:
        # User provided blocks - just intersect with CGI+shore
        eprint(f'  Using provided blocks file: {args.blocks}')
        blocks_path = args.blocks
    else:
        # Segment genome using wgbstools
        eprint('  Running wgbstools segment on reference betas...')
        seg_dir = ensure_dir(op.join(args.out_dir, 'segmentation'))
        blocks_path = op.join(seg_dir, f'segments.{args.genome}.bed')

        if not op.isfile(blocks_path) or args.force:
            beta_files = ' '.join(args.betas[:10])  # Use subset for segmentation
            cmd = f'{wgbstools} segment --betas {beta_files} -o {blocks_path}'
            if args.threads:
                cmd += f' -@ {args.threads}'
            run_cmd(cmd, verbose=args.verbose)
        else:
            eprint(f'  Segments exist: {blocks_path}')

    # Intersect blocks with CGI+shore regions
    eprint('  Intersecting blocks with CGI+shore regions...')
    if not check_executable('bedtools'):
        raise RuntimeError('bedtools not found. Please install bedtools.')

    intersect_path = op.join(args.out_dir, f'cgi_shore_blocks_raw.{args.genome}.bed')
    cmd = (f'bedtools intersect -a {blocks_path} -b {cgi_shore_bed} '
           f'-wa -f 0.5 | sort -k1,1 -k2,2n | uniq > {intersect_path}')
    run_cmd(cmd, verbose=args.verbose)

    # Add CpG indices if not present
    eprint('  Adding CpG indices...')
    df = pd.read_csv(intersect_path, sep='\t', header=None, comment='#')

    # Drop header rows that bedtools may pass through (e.g., rows where
    # column 1 is non-numeric like "start" or "chromStart")
    df = df[pd.to_numeric(df.iloc[:, 1], errors='coerce').notna()].reset_index(drop=True)

    if df.shape[1] >= 5:
        # Already has CpG indices (startCpG, endCpG in cols 3,4)
        df.columns = ['chr', 'start', 'end', 'startCpG', 'endCpG'] + \
                     [f'col{i}' for i in range(5, df.shape[1])]
    else:
        # Need to add CpG indices via wgbstools
        cmd = f'{wgbstools} convert --add_cpg_count -L {intersect_path}'
        result = run_cmd(cmd, verbose=args.verbose)
        # Parse output
        lines = result.strip().split('\n')
        rows = [line.split('\t') for line in lines if line and not line.startswith('#')]
        df = pd.DataFrame(rows)
        df.columns = ['chr', 'start', 'end', 'startCpG', 'endCpG'] + \
                     [f'col{i}' for i in range(5, df.shape[1])]

    eprint(f'  Loaded {len(df)} blocks from intersection')

    # Convert types and ensure integer CpG indices
    # CRITICAL: pd.to_numeric with errors='coerce' produces float64.
    # wgbstools load_blocks_file expects integer CpG indices (e.g. "12345"
    # not "12345.0"). Writing floats causes all methylation lookups to fail.
    for col in ['start', 'end', 'startCpG', 'endCpG']:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df = df.dropna(subset=['start', 'end', 'startCpG', 'endCpG'])
    for col in ['start', 'end', 'startCpG', 'endCpG']:
        df[col] = df[col].astype(int)

    # Filter by CpG count and bp length
    df['lenCpG'] = df['endCpG'] - df['startCpG']
    df['lenBp'] = df['end'] - df['start']

    orig_count = len(df)
    df = df[(df['lenCpG'] >= args.min_cpg) & (df['lenCpG'] <= args.max_cpg)]
    df = df[(df['lenBp'] >= args.min_bp) & (df['lenBp'] <= args.max_bp)]
    df = df.drop_duplicates(subset=['chr', 'start', 'end'])

    # CRITICAL: wgbstools find_markers requires blocks sorted by startCpG
    # (sort -k4,4n). Sorting by chr,start instead causes silent failures
    # in beta file lookups.
    df = df.sort_values('startCpG').reset_index(drop=True)

    eprint(f'  {orig_count} raw blocks -> {len(df)} filtered blocks '
           f'({args.min_cpg}-{args.max_cpg} CpGs, {args.min_bp}-{args.max_bp} bp)')

    if len(df) == 0:
        eprint('  ERROR: No blocks remain after filtering! Check your input blocks file.')
        raise RuntimeError('No candidate blocks found in CGI+shore regions.')

    # Write output (5 columns: chr, start, end, startCpG, endCpG)
    # Use integer format explicitly to avoid "12345.0" float representation
    df[['chr', 'start', 'end', 'startCpG', 'endCpG']].to_csv(
        out_path, sep='\t', header=False, index=False)

    # Diagnostic: print first few lines of output
    if args.verbose:
        eprint(f'  First 3 lines of {out_path}:')
        with open(out_path) as f:
            for i, line in enumerate(f):
                if i >= 3:
                    break
                eprint(f'    {line.rstrip()}')

    eprint(f'  Written: {out_path}')
    return out_path


###############################################################################
#                                                                             #
#   Stage 3: Find Tissue-Specific Markers                                    #
#                                                                             #
###############################################################################

def run_find_markers(wgbstools, blocks_path, groups_file, betas, out_dir,
                     delta_means, top, threads, verbose):
    """Run wgbstools find_markers with specified parameters."""
    ensure_dir(out_dir)
    beta_str = ' '.join(betas)

    cmd = (f'{wgbstools} find_markers '
           f'--blocks_path {blocks_path} '
           f'--groups_file {groups_file} '
           f'--betas {beta_str} '
           f'--delta_means {delta_means} '
           f'--min_cpg {3} '
           f'--max_cpg {50} '
           f'--min_bp {50} '
           f'--max_bp {5000} '
           f'--sort_by delta_means '
           f'--out_dir {out_dir} '
           f'--pval 0.05')
    if top:
        cmd += f' --top {top}'
    if threads:
        cmd += f' -@ {threads}'
    if verbose:
        cmd += ' -v'

    run_cmd(cmd, verbose=verbose)


def read_marker_file(fpath):
    """Read a wgbstools find_markers output file, skipping #>/#< comment lines
    but preserving the #chr header."""
    header_lines = 0
    with open(fpath) as fh:
        for line in fh:
            if line.startswith('#>') or line.startswith('#<'):
                header_lines += 1
            else:
                break
    return pd.read_csv(fpath, sep='\t', skiprows=header_lines)


def count_markers_per_group(markers_dir):
    """Count markers in each per-cell-type marker file."""
    counts = {}
    for f in glob.glob(op.join(markers_dir, 'Markers.*.bed')):
        group = op.basename(f).replace('Markers.', '').replace('.bed', '')
        df = read_marker_file(f)
        counts[group] = len(df)
    return counts


def merge_marker_files(markers_dir, out_path):
    """Merge per-cell-type marker files into a single TSV."""
    all_dfs = []
    for f in sorted(glob.glob(op.join(markers_dir, 'Markers.*.bed'))):
        df = read_marker_file(f)
        if not df.empty:
            all_dfs.append(df)

    if not all_dfs:
        eprint('  WARNING: No markers found!')
        return pd.DataFrame()

    merged = pd.concat(all_dfs, ignore_index=True)
    # wgbstools find_markers uses '#chr' as column name
    chr_col = '#chr' if '#chr' in merged.columns else 'chr'
    merged.sort_values(by=[chr_col, 'start'], inplace=True)
    merged.to_csv(out_path, sep='\t', index=False)
    eprint(f'  Merged {len(merged)} markers -> {out_path}')
    return merged


def _run_legacy_atlas_build(args, markers_path, atlas_path):
    """Run the original UXM atlas builder to preserve legacy row behavior."""
    build_py = op.join(args.uxm_path, 'src', 'build.py')
    if not op.isfile(build_py):
        raise RuntimeError(f'UXM_deconv build.py not found at {build_py}')

    pat_str = ' '.join(args.pats)
    cmd = (f'python3 {build_py} '
           f'--markers {markers_path} '
           f'--groups {args.groups} '
           f'--pats {pat_str} '
           f'--output {atlas_path} '
           f'-l {args.rlen} '
           f'--use_um')
    if args.threads:
        cmd += f' -@ {args.threads}'
    if args.verbose:
        cmd += ' -v'

    eprint('  Building legacy atlas row template (UXM build.py)...')
    run_cmd(cmd, verbose=args.verbose)
    if not op.isfile(atlas_path):
        raise RuntimeError(f'Legacy atlas was not created: {atlas_path}')


def _select_legacy_coordinate_columns(df):
    """Detect chromosome/start/end columns in an atlas dataframe."""
    chr_col = 'chr' if 'chr' in df.columns else '#chr' if '#chr' in df.columns else 'chrom' if 'chrom' in df.columns else None
    if chr_col is None:
        raise RuntimeError(f'Cannot find chromosome column in atlas: {list(df.columns)}')
    if 'start' not in df.columns or 'end' not in df.columns:
        raise RuntimeError(f'Cannot find start/end columns in atlas: {list(df.columns)}')
    return chr_col, 'start', 'end'


def _select_legacy_cpg_columns(df):
    """Detect startCpG/endCpG columns in an atlas dataframe."""
    start_cpg_col = None
    end_cpg_col = None
    for col in df.columns:
        key = str(col).strip().lower().replace('_', '')
        if key == 'startcpg':
            start_cpg_col = col
        elif key == 'endcpg':
            end_cpg_col = col
    if start_cpg_col is None or end_cpg_col is None:
        raise RuntimeError(
            f'Cannot find startCpG/endCpG columns in legacy atlas: {list(df.columns)}'
        )
    return start_cpg_col, end_cpg_col


def _load_groups_for_stage4(groups_path):
    """Load groups CSV into a sample->group map (matching Java stripping rules)."""
    groups_df = pd.read_csv(groups_path, comment='#')
    cols_lower = {str(c).lower(): c for c in groups_df.columns}
    if 'name' in cols_lower and 'group' in cols_lower:
        name_col = cols_lower['name']
        group_col = cols_lower['group']
    elif groups_df.shape[1] >= 2:
        name_col = groups_df.columns[0]
        group_col = groups_df.columns[1]
    else:
        raise RuntimeError(f'Groups file must contain at least 2 columns: {groups_path}')

    names = (
        groups_df[name_col].astype(str)
        .str.replace(r'\.(pat\.gz|beta|pat)$', '', regex=True)
        .tolist()
    )
    groups = groups_df[group_col].astype(str).tolist()
    sample_to_group = dict(zip(names, groups))
    return sample_to_group


def _beta_basename_to_sample_name_local(beta_path):
    """Mirror Java beta filename stripping (.beta/.lbeta)."""
    name = Path(beta_path).name
    if name.endswith('.beta'):
        return name[:-5]
    if name.endswith('.lbeta'):
        return name[:-6]
    return name


def _load_beta_file_to_markers_by_cpg_index(beta_path, start_cpg, end_cpg):
    """Aggregate per-marker (methy,total) using Java marker-mode CpG index logic.

    Java equivalent:
        startIdx = startCpG - 1
        endIdx = endCpG - 1   // exclusive
        for i in [startIdx, endIdx): sum beta[i]
    """
    with open(beta_path, 'rb') as fh:
        data = np.frombuffer(fh.read(), dtype=np.uint8)
    if data.size % 2 != 0:
        raise RuntimeError(f'Beta file size not divisible by 2: {beta_path}')
    per_cpg = data.reshape((-1, 2)).astype(np.int64)
    num_sites = per_cpg.shape[0]

    n_markers = len(start_cpg)
    out = np.zeros((n_markers, 2), dtype=np.int64)
    for mi in range(n_markers):
        s = int(start_cpg[mi])
        e = int(end_cpg[mi])
        # Keep behavior aligned with Java assumptions (1-based positive indices).
        if s <= 0 or e <= 0:
            continue
        lo = s - 1
        hi = e - 1  # exclusive
        if hi <= lo or lo >= num_sites:
            continue
        hi = min(hi, num_sites)
        block = per_cpg[lo:hi]
        out[mi, 0] = int(block[:, 0].sum())
        out[mi, 1] = int(block[:, 1].sum())
    return out


def _build_count_output_frames(legacy_atlas_df, cell_types, agg_methy, agg_total):
    """Build Stage-4 outputs while preserving legacy atlas row identity/order.

    Returns:
        panel_df: chrom/start/end + per-cell-type k/n values
        atlas_df: legacy first 8 columns + per-cell-type k/n values
    """
    if legacy_atlas_df.shape[1] < 8:
        raise RuntimeError(
            f'Legacy atlas must have at least 8 metadata columns, got {legacy_atlas_df.shape[1]}'
        )

    chr_col, start_col, end_col = _select_legacy_coordinate_columns(legacy_atlas_df)
    start_vals = pd.to_numeric(legacy_atlas_df[start_col], errors='coerce')
    end_vals = pd.to_numeric(legacy_atlas_df[end_col], errors='coerce')
    if start_vals.isna().any() or end_vals.isna().any():
        raise RuntimeError('Legacy atlas has non-numeric start/end values.')

    marker_chroms = legacy_atlas_df[chr_col].astype(str).to_numpy()
    marker_starts = start_vals.astype(np.int64).to_numpy()
    marker_ends = end_vals.astype(np.int64).to_numpy()
    n_markers = len(legacy_atlas_df)

    def _kn_col(ci):
        return [f'{agg_methy[mi, ci]}/{agg_total[mi, ci]}' for mi in range(n_markers)]

    # Reference panel for -refPanel mode.
    panel_df = pd.DataFrame({
        'chrom': marker_chroms,
        'start': marker_starts,
        'end': marker_ends,
    })
    for ci, ct in enumerate(cell_types):
        panel_df[ct] = _kn_col(ci)

    # Atlas output: preserve first 8 metadata columns exactly as legacy.
    atlas_df = legacy_atlas_df.iloc[:, :8].copy()
    for ci, ct in enumerate(cell_types):
        atlas_df[ct] = _kn_col(ci)

    return panel_df, atlas_df


def stage3_find_markers(args, blocks_path):
    """Find tissue-specific markers with adaptive thresholds."""
    eprint('\n=== Stage 3: Find Tissue-Specific Markers ===')

    wgbstools = op.join(args.wgbstools_path, 'wgbstools')
    if not op.isfile(wgbstools):
        wgbstools = 'wgbstools'

    markers_dir = ensure_dir(op.join(args.out_dir, 'markers'))
    merged_path = op.join(args.out_dir,
                          f'Markers.CGI_shore.U{args.num_markers}.{args.genome}.tsv')

    if op.isfile(merged_path) and not args.force:
        eprint(f'  Output exists: {merged_path}')
        return merged_path

    # Pre-flight diagnostics
    eprint(f'  Blocks file: {blocks_path}')
    blocks_check = pd.read_csv(blocks_path, sep='\t', header=None, nrows=3)
    eprint(f'  Blocks columns: {blocks_check.shape[1]}, rows (sample): {blocks_check.shape[0]}')
    eprint(f'  Blocks dtypes: {dict(blocks_check.dtypes)}')
    eprint(f'  Blocks head:\n{blocks_check.to_string(index=False, header=False)}')

    eprint(f'  Beta files: {len(args.betas)} files')
    if args.verbose and args.betas:
        eprint(f'    First: {args.betas[0]}')
        eprint(f'    Last:  {args.betas[-1]}')

    # Load groups to know expected cell types
    groups_df = pd.read_csv(args.groups)
    expected_groups = sorted(groups_df['group'].unique())
    eprint(f'  Expected cell types: {len(expected_groups)}')
    eprint(f'  Groups file samples: {len(groups_df)}')

    # CRITICAL FIX: wgbstools find_markers expects the groups file 'name'
    # column to contain sample prefixes WITHOUT file extensions. It appends
    # '.beta' itself when matching (via match_prefix_to_bin). If the groups
    # file has names like "sample.pat.gz", it will look for "sample.pat.gz.beta"
    # which never matches, causing a silent failure with 0 markers.
    name_col = 'name' if 'name' in groups_df.columns else groups_df.columns[0]
    needs_fix = False
    for ext in ['.pat.gz', '.beta', '.pat']:
        if groups_df[name_col].str.endswith(ext).any():
            needs_fix = True
            eprint(f'  WARNING: Groups file names contain "{ext}" extension.')
            eprint(f'    wgbstools find_markers expects names WITHOUT extensions.')
            break

    if needs_fix:
        eprint(f'  Auto-fixing: stripping extensions from groups file names...')
        fixed_groups_path = op.join(args.out_dir, 'groups_fixed.csv')
        fixed_df = groups_df.copy()
        # Strip common extensions: .pat.gz, .beta, .pat
        fixed_df[name_col] = fixed_df[name_col].str.replace(r'\.(pat\.gz|beta|pat)$', '', regex=True)
        fixed_df.to_csv(fixed_groups_path, index=False)
        args.groups = fixed_groups_path
        eprint(f'  Fixed groups file written: {fixed_groups_path}')
        eprint(f'  Sample names now: {fixed_df[name_col].iloc[0]}')
        groups_df = fixed_df

    # Verify group names match beta file names.
    #
    # CRITICAL: wgbstools find_markers strips ONLY the `.beta` extension
    # from the beta file basename to get the sample name, then looks for
    # that name in the groups file. So the groups file names must match
    # the beta basename minus `.beta` EXACTLY.
    #
    # Common naming conventions that cause mismatches:
    #   hg19: beta = "Sample1.beta"       → wgbstools name = "Sample1"       ✓ matches groups
    #   hg38: beta = "Sample1.hg38.beta"  → wgbstools name = "Sample1.hg38" ✗ groups says "Sample1"
    #
    # We must match using the SINGLE-split basename (what wgbstools sees),
    # not the double-split one (which would hide the `.hg38` mismatch and
    # let the user think everything is fine while find_markers silently
    # produces 0 results).
    group_names = set(groups_df[name_col])
    # Single-split: strip .beta only, keeping any genome suffix
    beta_basenames_wgbstools = {
        op.splitext(op.basename(b))[0]  # "Sample1.hg38.beta" → "Sample1.hg38"
        for b in args.betas
    }
    matched = group_names & beta_basenames_wgbstools

    if not matched:
        # Check if appending the genome suffix (.hg38 / .hg19) to the group
        # names produces a match — this is the most common hg38 naming issue.
        for suffix in [f'.{args.genome}', '.hg38', '.hg19']:
            test_names = {n + suffix for n in group_names}
            test_matched = test_names & beta_basenames_wgbstools
            if test_matched:
                eprint(f'  Auto-fixing: appending "{suffix}" to group names to match beta files...')
                eprint(f'    Group name:    {sorted(group_names)[0]}')
                eprint(f'    Beta basename: {sorted(beta_basenames_wgbstools)[0]}')
                eprint(f'    After fix:     {sorted(test_names)[0]}')
                fixed_groups_path = op.join(args.out_dir, 'groups_genome_suffix_fixed.csv')
                fixed_df = groups_df.copy()
                fixed_df[name_col] = fixed_df[name_col] + suffix
                fixed_df.to_csv(fixed_groups_path, index=False)
                args.groups = fixed_groups_path
                groups_df = fixed_df
                group_names = set(fixed_df[name_col])
                matched = group_names & beta_basenames_wgbstools
                eprint(f'  Fixed groups file written: {fixed_groups_path}')
                break

    unmatched = group_names - beta_basenames_wgbstools
    if unmatched:
        eprint(f'  WARNING: {len(unmatched)} group names not found in beta files:')
        for name in sorted(list(unmatched)[:5]):
            eprint(f'    {name}')
        if len(unmatched) > 5:
            eprint(f'    ... and {len(unmatched) - 5} more')
    if not matched:
        eprint(f'  ERROR: No group names match any beta file!')
        eprint(f'    Group names (first 3): {sorted(list(group_names))[:3]}')
        eprint(f'    Beta basenames (first 3): {sorted(list(beta_basenames_wgbstools))[:3]}')
        raise RuntimeError(
            'Groups file names do not match beta file names. '
            f'wgbstools strips only ".beta" from the filename — if your '
            f'beta files are named like "Sample.{args.genome}.beta", the '
            f'groups file "name" column must be "Sample.{args.genome}" '
            f'(not just "Sample").')
    eprint(f'  Matched samples: {len(matched)}/{len(group_names)}')

    # Adaptive threshold relaxation
    delta_thresholds = [args.delta_means, 0.25, 0.2, 0.15]
    shortfall_groups = set(expected_groups)

    # Use absolute paths for beta files to avoid working directory issues
    abs_betas = [op.abspath(b) for b in args.betas]
    abs_groups = op.abspath(args.groups)
    abs_blocks = op.abspath(blocks_path)

    for i, delta in enumerate(delta_thresholds):
        if not shortfall_groups:
            break

        pass_label = f'Pass {i + 1}' if i > 0 else 'Initial pass'
        eprint(f'\n  {pass_label}: delta_means={delta}')

        pass_dir = ensure_dir(op.join(args.out_dir, 'markers', f'pass_{i}'))

        # Run find_markers for groups that still need more markers
        targets_str = ' '.join(shortfall_groups) if i > 0 else None
        target_flag = f'--targets {targets_str}' if targets_str else ''

        beta_str = ' '.join(abs_betas)
        cmd = (f'{wgbstools} find_markers '
               f'--blocks_path {abs_blocks} '
               f'--groups_file {abs_groups} '
               f'--betas {beta_str} '
               f'--delta_means {delta} '
               f'--min_cpg {args.min_cpg} '
               f'--max_cpg {args.max_cpg} '
               f'--min_bp {args.min_bp} '
               f'--max_bp {args.max_bp} '
               f'--sort_by delta_means '
               f'--top {args.num_markers} '
               f'--out_dir {pass_dir} '
               f'--pval 0.05')
        if args.unmeth_mean_thresh is not None:
            cmd += f' --unmeth_mean_thresh {args.unmeth_mean_thresh}'
        if args.meth_mean_thresh is not None:
            cmd += f' --meth_mean_thresh {args.meth_mean_thresh}'
        if target_flag:
            cmd += f' {target_flag}'
        if args.threads:
            cmd += f' -@ {args.threads}'
        if args.verbose:
            cmd += ' -v'

        try:
            run_cmd(cmd, verbose=args.verbose)
        except RuntimeError as e:
            err_msg = str(e)
            if 'mismatch' in err_msg.lower() or 'not found' in err_msg.lower():
                eprint(f'  ERROR: find_markers failed - likely a groups/beta file name mismatch!')
                eprint(f'  {e}')
                raise
            eprint(f'  WARNING: find_markers failed for pass {i}: {e}')
            continue

        # Check which groups still have shortfall
        counts = count_markers_per_group(pass_dir)
        new_shortfall = set()
        for group in shortfall_groups:
            n = counts.get(group, 0)
            if n < 5:
                new_shortfall.add(group)
                eprint(f'    {group}: {n} markers (shortfall)')
            else:
                eprint(f'    {group}: {n} markers')

        # Copy results for groups that passed threshold
        for group in shortfall_groups - new_shortfall:
            src = op.join(pass_dir, f'Markers.{group}.bed')
            dst = op.join(markers_dir, f'Markers.{group}.bed')
            if op.isfile(src):
                import shutil
                shutil.copy2(src, dst)

        # For groups that had shortfall and got relaxed results, also copy
        if i == len(delta_thresholds) - 1:
            # Last pass - accept whatever we got
            for group in new_shortfall:
                src = op.join(pass_dir, f'Markers.{group}.bed')
                dst = op.join(markers_dir, f'Markers.{group}.bed')
                if op.isfile(src):
                    import shutil
                    shutil.copy2(src, dst)
                    n = counts.get(group, 0)
                    eprint(f'    WARNING: {group} only has {n} markers '
                           f'(even after relaxation)')

        shortfall_groups = new_shortfall

    # Merge all per-cell-type markers
    eprint('\n  Merging markers across cell types...')
    merged_df = merge_marker_files(markers_dir, merged_path)

    # Print summary
    final_counts = count_markers_per_group(markers_dir)
    eprint(f'\n  Marker summary:')
    for group in sorted(final_counts):
        eprint(f'    {group}: {final_counts[group]}')
    eprint(f'  Total: {sum(final_counts.values())} markers '
           f'across {len(final_counts)} cell types')

    return merged_path


###############################################################################
#                                                                             #
#   Stage 4: Build Reference Panel                                           #
#                                                                             #
###############################################################################

def stage4_build_reference_panel(args, markers_path):
    """Build count-based outputs with legacy atlas row behavior.

    Row identity/order is produced by the original UXM ``build.py`` pipeline
    (legacy behavior), then per-cell-type value columns are replaced by
    ``methylated_count/total_count`` derived from reference .beta files.
    """
    eprint('\n=== Stage 4: Build Reference Panel ===')

    panel_path = op.join(
        args.out_dir, f'ReferencePanel.CGI_shore.U{args.num_markers}.{args.genome}.tsv'
    )
    atlas_path = op.join(
        args.out_dir, f'Atlas.CGI_shore.U{args.num_markers}.l{args.rlen}.{args.genome}.tsv'
    )

    if op.isfile(panel_path) and op.isfile(atlas_path) and not args.force:
        eprint(f'  Outputs exist: {panel_path} and {atlas_path}')
        return panel_path, atlas_path

    # Build or reuse legacy atlas row template.
    if op.isfile(atlas_path) and not args.force:
        eprint(f'  Reusing existing legacy atlas template: {atlas_path}')
    else:
        _run_legacy_atlas_build(args, markers_path, atlas_path)

    legacy_atlas_df = pd.read_csv(atlas_path, sep='\t')
    if legacy_atlas_df.empty:
        raise RuntimeError(f'Legacy atlas has no rows: {atlas_path}')
    if legacy_atlas_df.shape[1] < 8:
        raise RuntimeError(
            f'Legacy atlas must have >=8 columns, got {legacy_atlas_df.shape[1]} in {atlas_path}'
        )

    legacy_chr_col, legacy_start_col, legacy_end_col = _select_legacy_coordinate_columns(legacy_atlas_df)
    marker_chroms = legacy_atlas_df[legacy_chr_col].astype(str).to_numpy()
    marker_starts = pd.to_numeric(legacy_atlas_df[legacy_start_col], errors='coerce').astype(np.int64).to_numpy()
    marker_ends = pd.to_numeric(legacy_atlas_df[legacy_end_col], errors='coerce').astype(np.int64).to_numpy()
    n_markers = len(legacy_atlas_df)
    eprint(f'  Legacy atlas markers: {n_markers}')

    # Use startCpG/endCpG ranges from legacy atlas so Stage 4 k/n exactly matches
    # Java BetaValueDeconvolution -markerRegions reference extraction behavior.
    start_cpg_col, end_cpg_col = _select_legacy_cpg_columns(legacy_atlas_df)
    start_cpg = pd.to_numeric(legacy_atlas_df[start_cpg_col], errors='coerce').astype(np.int64).to_numpy()
    end_cpg = pd.to_numeric(legacy_atlas_df[end_cpg_col], errors='coerce').astype(np.int64).to_numpy()
    eprint('  Stage 4 k/n extraction mode: Java-compatible startCpG/endCpG index slicing')

    # Load groups
    sample_to_group = _load_groups_for_stage4(args.groups)
    # Preserve legacy atlas cell-type column order when available.
    atlas_cell_types = [str(c) for c in legacy_atlas_df.columns[8:]]
    if atlas_cell_types:
        cell_types = atlas_cell_types
    else:
        cell_types = sorted(set(sample_to_group.values()))
    ct_to_index = {c: i for i, c in enumerate(cell_types)}
    n_ct = len(cell_types)
    eprint(f'  Cell types: {n_ct} — {cell_types}')

    # Aggregate per-marker (methylated, total) counts per cell type
    agg_methy = np.zeros((n_markers, n_ct), dtype=np.int64)
    agg_total = np.zeros((n_markers, n_ct), dtype=np.int64)
    n_matched = 0
    skipped_group_mismatch = 0

    for beta_path_str in args.betas:
        beta_path = Path(beta_path_str)
        sample_name = _beta_basename_to_sample_name_local(beta_path)
        group = sample_to_group.get(sample_name)
        if group is None:
            continue
        if group not in ct_to_index:
            skipped_group_mismatch += 1
            continue
        ci = ct_to_index[group]
        counts = _load_beta_file_to_markers_by_cpg_index(beta_path, start_cpg, end_cpg)
        agg_methy[:, ci] += counts[:, 0]
        agg_total[:, ci] += counts[:, 1]
        n_matched += 1

    eprint(f'  Matched beta files: {n_matched}/{len(args.betas)}')
    if skipped_group_mismatch > 0:
        eprint(f'  WARNING: skipped {skipped_group_mismatch} matched betas due to group not present in legacy atlas columns.')

    if n_matched == 0:
        eprint('  ERROR: No beta files matched the groups file!')
        raise RuntimeError('No beta files matched — check groups file sample names.')

    panel_df, atlas_df = _build_count_output_frames(
        legacy_atlas_df=legacy_atlas_df,
        cell_types=cell_types,
        agg_methy=agg_methy,
        agg_total=agg_total,
    )
    panel_df.to_csv(panel_path, sep='\t', index=False)
    atlas_df.to_csv(atlas_path, sep='\t', index=False)

    eprint(f'  Reference panel written: {panel_path}')
    eprint(f'  Atlas-style k/n file written: {atlas_path} (legacy row structure preserved)')
    eprint(f'  Shape: {n_markers} markers x {n_ct} cell types')
    eprint(f'  Format: methylated_count/total_count per cell')
    eprint(f'  Usage:')
    eprint(f'    finaleme-too run --reference-panel {panel_path} ...')
    eprint(f'    BetaValueDeconvolution -refPanel {panel_path} ...')
    eprint(f'    BetaValueDeconvolution -markerRegions {atlas_path} ...')

    return panel_path, atlas_path


def _stage4_fallback_slow(args, legacy_atlas_df, chr_col, start_col, end_col, panel_path, atlas_path):
    """Fallback .beta parser when finaleme_too is not installed.

    Reads each .beta file as raw binary (uint8 pairs), looks up CpG
    positions from the wgbstools CpG index, and aggregates to marker
    regions with a double binary search. Much slower than the finaleme_too
    version (no numpy vectorization of the CpG index) but works without
    the package installed.
    """
    eprint('  [fallback] Loading CpG index...')
    cpg_index_path = op.join(args.wgbstools_path, 'references', args.genome,
                             'CpG.bed.gz')
    if not op.isfile(cpg_index_path):
        cpg_index_path = op.join(args.wgbstools_path, 'references', args.genome,
                                 'CpG.bed')
    if not op.isfile(cpg_index_path):
        raise RuntimeError(f'CpG index not found at {cpg_index_path}')

    # In fallback mode we still use Java-compatible CpG index slicing from
    # startCpG/endCpG in the legacy atlas; no coordinate overlap logic.
    start_cpg_col, end_cpg_col = _select_legacy_cpg_columns(legacy_atlas_df)
    start_cpg = pd.to_numeric(legacy_atlas_df[start_cpg_col], errors='coerce').astype(np.int64).to_numpy()
    end_cpg = pd.to_numeric(legacy_atlas_df[end_cpg_col], errors='coerce').astype(np.int64).to_numpy()
    eprint('  [fallback] Stage 4 k/n extraction mode: Java-compatible startCpG/endCpG index slicing')

    # Load groups
    groups_df = pd.read_csv(args.groups)
    name_col = 'name' if 'name' in groups_df.columns else groups_df.columns[0]
    groups_df['name_stripped'] = (
        groups_df[name_col].astype(str)
        .str.replace(r'\.(pat\.gz|beta|pat)$', '', regex=True)
    )
    sample_to_group = dict(zip(groups_df['name_stripped'], groups_df['group']))
    atlas_cell_types = [str(c) for c in legacy_atlas_df.columns[8:]]
    if atlas_cell_types:
        cell_types = atlas_cell_types
    else:
        cell_types = sorted(set(groups_df['group']))
    ct_to_index = {c: i for i, c in enumerate(cell_types)}

    n_markers = len(legacy_atlas_df)
    n_ct = len(cell_types)
    agg_methy = np.zeros((n_markers, n_ct), dtype=np.int64)
    agg_total = np.zeros((n_markers, n_ct), dtype=np.int64)
    n_matched = 0
    skipped_group_mismatch = 0

    for beta_path_str in args.betas:
        beta_path = Path(beta_path_str)
        name = beta_path.name
        for ext in ('.beta', '.lbeta'):
            if name.endswith(ext):
                name = name[:-len(ext)]
                break
        group = sample_to_group.get(name)
        if group is None:
            continue
        if group not in ct_to_index:
            skipped_group_mismatch += 1
            continue
        ci = ct_to_index[group]

        counts = _load_beta_file_to_markers_by_cpg_index(beta_path, start_cpg, end_cpg)
        agg_methy[:, ci] += counts[:, 0]
        agg_total[:, ci] += counts[:, 1]
        n_matched += 1

    eprint(f'  [fallback] Matched beta files: {n_matched}/{len(args.betas)}')
    if skipped_group_mismatch > 0:
        eprint(f'  [fallback] WARNING: skipped {skipped_group_mismatch} matched betas due to group not present in legacy atlas columns.')

    if n_matched == 0:
        raise RuntimeError('No beta files matched — check groups file names and atlas cell-type columns.')

    panel_df, atlas_df = _build_count_output_frames(
        legacy_atlas_df=legacy_atlas_df,
        cell_types=cell_types,
        agg_methy=agg_methy,
        agg_total=agg_total,
    )
    panel_df.to_csv(panel_path, sep='\t', index=False)
    atlas_df.to_csv(atlas_path, sep='\t', index=False)
    eprint(f'  [fallback] Reference panel written: {panel_path}')
    eprint(f'  [fallback] Atlas-style k/n file written: {atlas_path}')
    return panel_path, atlas_path


###############################################################################
#                                                                             #
#   Stage 5: Validation Report                                               #
#                                                                             #
###############################################################################

def stage5_validation_report(args, cgi_shore_bed, blocks_path, markers_path,
                             panel_path, atlas_path):
    """Generate validation report."""
    eprint('\n=== Stage 5: Validation Report ===')

    report_path = op.join(args.out_dir, 'report.txt')
    lines = []
    lines.append('=' * 70)
    lines.append('CGI+Shore Marker Generation Report')
    lines.append('=' * 70)
    lines.append('')

    # CGI+shore coverage
    lines.append('--- CGI+Shore Regions ---')
    cgi_df = pd.read_csv(cgi_shore_bed, sep='\t', header=None,
                         names=['chrom', 'start', 'end'])
    total_regions = len(cgi_df)
    total_bp = (cgi_df['end'] - cgi_df['start']).sum()
    lines.append(f'Total CGI+shore regions: {total_regions:,}')
    lines.append(f'Total coverage: {total_bp:,} bp ({total_bp / 1e6:.1f} Mb)')
    lines.append('')

    # Candidate blocks
    lines.append('--- Candidate Blocks ---')
    blocks_df = pd.read_csv(blocks_path, sep='\t', header=None,
                            names=['chr', 'start', 'end', 'startCpG', 'endCpG'])
    lines.append(f'Total candidate blocks: {len(blocks_df):,}')
    lines.append('')

    # Markers per cell type
    lines.append('--- Markers per Cell Type ---')
    if op.isfile(markers_path):
        markers_df = pd.read_csv(markers_path, sep='\t')
        target_col = 'target' if 'target' in markers_df.columns else '#target'
        if target_col in markers_df.columns:
            counts = markers_df[target_col].value_counts().sort_index()
            for group, count in counts.items():
                status = '' if count >= args.num_markers * 0.5 else ' [SHORTFALL]'
                lines.append(f'  {group}: {count}{status}')
            lines.append(f'  Total: {len(markers_df)} markers '
                         f'across {len(counts)} cell types')

            # Statistics
            lines.append('')
            lines.append('--- Marker Quality Statistics ---')
            if 'delta_means' in markers_df.columns:
                lines.append(f'Mean delta_means: {markers_df["delta_means"].mean():.3f}')
                lines.append(f'Min delta_means: {markers_df["delta_means"].min():.3f}')
            if 'ttest' in markers_df.columns:
                pvals = pd.to_numeric(markers_df['ttest'], errors='coerce')
                lines.append(f'Median p-value: {pvals.median():.2e}')
    else:
        lines.append('  Markers file not found!')
    lines.append('')

    # Overlap with original UXM markers
    lines.append('--- Overlap with Original UXM U250 Markers ---')
    orig_markers_path = op.join(
        args.uxm_path, 'supplemental', f'Markers.U250.{args.genome}.tsv')
    if op.isfile(orig_markers_path):
        orig_df = pd.read_csv(orig_markers_path, sep='\t')
        chr_col = '#chr' if '#chr' in orig_df.columns else 'chr'

        # Count how many original markers fall in CGI+shore
        overlap_count = 0
        for _, cgi_row in cgi_df.iterrows():
            mask = ((orig_df[chr_col] == cgi_row['chrom']) &
                    (orig_df['start'] >= cgi_row['start']) &
                    (orig_df['end'] <= cgi_row['end']))
            overlap_count += mask.sum()

        lines.append(f'Original U250 markers: {len(orig_df)}')
        lines.append(f'Overlapping CGI+shore: {overlap_count} '
                     f'({100 * overlap_count / len(orig_df):.1f}%)')
    else:
        lines.append(f'  Original markers file not found at {orig_markers_path}')
    lines.append('')

    # Reference panel info
    lines.append('--- Reference Panel ---')
    if op.isfile(panel_path):
        panel_df = pd.read_csv(panel_path, sep='\t', nrows=5)
        n_ct = panel_df.shape[1] - 3  # subtract chrom, start, end
        lines.append(f'Reference panel: {panel_path}')
        lines.append(f'Format: methylated_count/total_count per cell type')
        lines.append(f'Cell types ({n_ct}): {list(panel_df.columns[3:])}')
        # Count total markers by reading just the row count
        full_panel = pd.read_csv(panel_path, sep='\t', usecols=[0])
        lines.append(f'Markers: {len(full_panel)}')
    else:
        lines.append('  Reference panel file not found!')
    lines.append('')

    # Atlas-style k/n file info
    lines.append('--- Atlas-Style k/n Output ---')
    if op.isfile(atlas_path):
        atlas_df = pd.read_csv(atlas_path, sep='\t', nrows=5)
        lines.append(f'Atlas-style file: {atlas_path}')
        lines.append('Format: UXM-atlas columns + methylated_count/total_count per cell type')
        lines.append(f'Columns: {list(atlas_df.columns[:8])} + {max(0, atlas_df.shape[1] - 8)} cell-type columns')
        full_atlas = pd.read_csv(atlas_path, sep='\t', usecols=[0])
        lines.append(f'Rows: {len(full_atlas)}')
    else:
        lines.append('  Atlas-style file not found!')
    lines.append('')

    lines.append('=' * 70)

    report_text = '\n'.join(lines)
    with open(report_path, 'w') as f:
        f.write(report_text)

    eprint(report_text)
    eprint(f'\nReport written: {report_path}')
    return report_path


###############################################################################
#                                                                             #
#   Main                                                                     #
#                                                                             #
###############################################################################

def parse_args():
    parser = argparse.ArgumentParser(
        description='Generate CGI+Shore restricted tissue-specific markers '
                    'for UXM deconvolution.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)

    # Required
    parser.add_argument('--genome', required=True, choices=['hg19', 'hg38'],
                        help='Genome build (hg19 or hg38)')
    parser.add_argument('--out-dir', required=True,
                        help='Output directory')
    parser.add_argument('--wgbstools-path', required=True,
                        help='Path to wgbs_tools installation directory')
    parser.add_argument('--uxm-path', required=True,
                        help='Path to UXM_deconv installation directory')

    # Reference data (optional if auto-download)
    parser.add_argument('--betas', nargs='+',
                        help='Reference WGBS beta files')
    parser.add_argument('--pats', nargs='+',
                        help='Reference WGBS pat.gz files')
    parser.add_argument('--groups',
                        help='Groups CSV file mapping samples to cell types')

    # CGI annotation
    parser.add_argument('--cgi-bed',
                        help='CpG Island BED file. Auto-downloaded if not provided.')
    parser.add_argument('--shore-size', type=int, default=2000,
                        help='Shore size in bp around CGI (default: 2000)')
    parser.add_argument('--chrom-sizes',
                        help='Chromosome sizes file for boundary checking')

    # Block generation
    parser.add_argument('--blocks',
                        help='Pre-existing blocks file (skip segmentation)')

    # Marker selection
    parser.add_argument('--num-markers', type=int, default=250,
                        help='Target number of markers per cell type '
                             '(default: 250)')
    parser.add_argument('--delta-means', type=float, default=0.4,
                        help='Initial delta_means threshold (default: 0.4)')
    parser.add_argument('--unmeth-mean-thresh', type=float, default=0.1,
                        help='Max mean methylation for the unmethylated group. '
                             'For binarized beta deconvolution (threshold 0.1), '
                             'use 0.1 so U-markers have target mean < 0.1. '
                             '(default: 0.1)')
    parser.add_argument('--meth-mean-thresh', type=float, default=0.5,
                        help='Min mean methylation for the methylated group. '
                             'For binarized beta deconvolution (threshold 0.1), '
                             'use 0.5 so background is well above the binarization '
                             'boundary. (default: 0.5)')
    parser.add_argument('--min-cpg', type=int, default=1,
                        help='Minimum CpGs per block (default: 1)')
    parser.add_argument('--max-cpg', type=int, default=1000,
                        help='Maximum CpGs per block (default: 1000)')
    parser.add_argument('--min-bp', type=int, default=50,
                        help='Minimum block length in bp (default: 50)')
    parser.add_argument('--max-bp', type=int, default=5000,
                        help='Maximum block length in bp (default: 5000)')

    # Atlas building
    parser.add_argument('--rlen', type=int, default=3,
                        help='Minimum CpGs per read for homog (default: 3)')

    # Runtime
    parser.add_argument('--threads', type=int, default=10,
                        help='Number of threads (default: 10)')
    parser.add_argument('--force', action='store_true',
                        help='Overwrite existing output files')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Verbose output')

    args = parser.parse_args()

    # Validate paths
    if not op.isdir(args.wgbstools_path):
        parser.error(f'wgbstools path not found: {args.wgbstools_path}')
    if not op.isdir(args.uxm_path):
        parser.error(f'UXM_deconv path not found: {args.uxm_path}')

    # Add wgbstools to PATH for subprocess calls
    wgbs_bin = args.wgbstools_path
    os.environ['PATH'] = wgbs_bin + ':' + os.environ.get('PATH', '')

    # Also add wgbs_tools Python path for potential imports
    wgbs_python = op.join(args.wgbstools_path, 'src', 'python')
    if op.isdir(wgbs_python):
        sys.path.insert(0, wgbs_python)

    ensure_dir(args.out_dir)

    return args


def main():
    args = parse_args()

    eprint('=' * 70)
    eprint('CGI+Shore Marker Generation Pipeline')
    eprint(f'  Genome: {args.genome}')
    eprint(f'  Target markers per cell type: {args.num_markers}')
    eprint(f'  Shore size: ±{args.shore_size}bp')
    eprint(f'  Output: {args.out_dir}')
    eprint('=' * 70)

    # Stage 0: Download reference data if needed
    if not args.betas or not args.pats or not args.groups:
        stage0_download_references(args)

    # Stage 1: Generate CGI+shore BED
    cgi_shore_bed = stage1_generate_cgi_shore_bed(args)

    # Stage 2: Generate candidate blocks
    blocks_path = stage2_generate_blocks(args, cgi_shore_bed)

    # Stage 3: Find tissue-specific markers
    markers_path = stage3_find_markers(args, blocks_path)

    # Stage 4: Build reference panel
    panel_path, atlas_path = stage4_build_reference_panel(args, markers_path)

    # Stage 5: Validation report
    stage5_validation_report(args, cgi_shore_bed, blocks_path,
                            markers_path, panel_path, atlas_path)

    eprint('\n' + '=' * 70)
    eprint('Pipeline complete!')
    eprint(f'  Markers: {markers_path}')
    eprint(f'  Reference panel: {panel_path}')
    eprint(f'  Atlas-style k/n file: {atlas_path}')
    eprint(f'  Use with:')
    eprint(f'    finaleme-too run --reference-panel {panel_path} \\')
    eprint(f'      --marker-regions {markers_path} ...')
    eprint('=' * 70)


if __name__ == '__main__':
    main()
