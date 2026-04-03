#!/usr/bin/env python3
"""
Generate CGI+Shore restricted tissue-specific methylation markers for UXM deconvolution.

This script creates tissue-specific methylation markers restricted to CpG Island (CGI)
and CGI shore regions (±2kb around CGI), then builds a UXM-compatible atlas for
tissues-of-origin deconvolution from FinaleMe cfDNA predictions.

Pipeline stages:
  Stage 0: Download reference WGBS data (optional, if not user-provided)
  Stage 1: Generate CGI+shore BED regions
  Stage 2: Generate candidate blocks within CGI+shore
  Stage 3: Find tissue-specific markers with adaptive thresholds
  Stage 4: Build UXM atlas
  Stage 5: Generate validation report

Usage:
  python generate_cgi_shore_markers.py --genome hg19 --num-markers 25 \
    --out-dir results/ --wgbstools-path /path/to/wgbs_tools \
    --uxm-path /path/to/UXM_deconv

Requirements:
  - Python 3.8+ with pandas, numpy, scipy
  - wgbs_tools (https://github.com/nloyfer/wgbs_tools)
  - UXM_deconv (https://github.com/nloyfer/UXM_deconv)
  - bedtools (for BED intersection)
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
from collections import defaultdict
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
    df = pd.read_csv(intersect_path, sep='\t', header=None)

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

    # Convert types
    for col in ['start', 'end', 'startCpG', 'endCpG']:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    # Filter by CpG count and bp length
    df['lenCpG'] = df['endCpG'] - df['startCpG']
    df['lenBp'] = df['end'] - df['start']

    orig_count = len(df)
    df = df[(df['lenCpG'] >= args.min_cpg) & (df['lenCpG'] <= args.max_cpg)]
    df = df[(df['lenBp'] >= args.min_bp) & (df['lenBp'] <= args.max_bp)]
    df = df.drop_duplicates(subset=['chr', 'start', 'end'])
    df = df.sort_values(['chr', 'start']).reset_index(drop=True)

    eprint(f'  {orig_count} raw blocks -> {len(df)} filtered blocks '
           f'({args.min_cpg}-{args.max_cpg} CpGs, {args.min_bp}-{args.max_bp} bp)')

    # Write output (5 columns: chr, start, end, startCpG, endCpG)
    df[['chr', 'start', 'end', 'startCpG', 'endCpG']].to_csv(
        out_path, sep='\t', header=False, index=False)

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


def count_markers_per_group(markers_dir):
    """Count markers in each per-cell-type marker file."""
    counts = {}
    for f in glob.glob(op.join(markers_dir, 'Markers.*.bed')):
        group = op.basename(f).replace('Markers.', '').replace('.bed', '')
        df = pd.read_csv(f, sep='\t', comment='#')
        counts[group] = len(df)
    return counts


def merge_marker_files(markers_dir, out_path):
    """Merge per-cell-type marker files into a single TSV."""
    all_dfs = []
    for f in sorted(glob.glob(op.join(markers_dir, 'Markers.*.bed'))):
        df = pd.read_csv(f, sep='\t', comment='#')
        if not df.empty:
            all_dfs.append(df)

    if not all_dfs:
        eprint('  WARNING: No markers found!')
        return pd.DataFrame()

    merged = pd.concat(all_dfs, ignore_index=True)
    merged.sort_values(by=['#chr', 'start'], inplace=True)
    merged.to_csv(out_path, sep='\t', index=False)
    eprint(f'  Merged {len(merged)} markers -> {out_path}')
    return merged


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

    # Load groups to know expected cell types
    groups_df = pd.read_csv(args.groups)
    expected_groups = sorted(groups_df['group'].unique())
    eprint(f'  Expected cell types: {len(expected_groups)}')

    # Adaptive threshold relaxation
    delta_thresholds = [args.delta_means, 0.25, 0.2, 0.15]
    shortfall_groups = set(expected_groups)

    for i, delta in enumerate(delta_thresholds):
        if not shortfall_groups:
            break

        pass_label = f'Pass {i + 1}' if i > 0 else 'Initial pass'
        eprint(f'\n  {pass_label}: delta_means={delta}')

        pass_dir = ensure_dir(op.join(markers_dir, f'pass_{i}'))

        # Run find_markers for groups that still need more markers
        targets_str = ' '.join(shortfall_groups) if i > 0 else None
        target_flag = f'--targets {targets_str}' if targets_str else ''

        beta_str = ' '.join(args.betas)
        cmd = (f'{wgbstools} find_markers '
               f'--blocks_path {blocks_path} '
               f'--groups_file {args.groups} '
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
        if target_flag:
            cmd += f' {target_flag}'
        if args.threads:
            cmd += f' -@ {args.threads}'
        if args.verbose:
            cmd += ' -v'

        try:
            run_cmd(cmd, verbose=args.verbose)
        except RuntimeError as e:
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
#   Stage 4: Build UXM Atlas                                                 #
#                                                                             #
###############################################################################

def stage4_build_atlas(args, markers_path):
    """Build UXM atlas from markers and reference pat files."""
    eprint('\n=== Stage 4: Build UXM Atlas ===')

    atlas_path = op.join(args.out_dir,
                         f'Atlas.CGI_shore.U{args.num_markers}.l{args.rlen}.{args.genome}.tsv')

    if op.isfile(atlas_path) and not args.force:
        eprint(f'  Output exists: {atlas_path}')
        return atlas_path

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
           f'--use_um')  # Use both U and M directions
    if args.threads:
        cmd += f' -@ {args.threads}'
    if args.verbose:
        cmd += ' -v'

    eprint(f'  Building atlas...')
    run_cmd(cmd, verbose=args.verbose)

    # Validate atlas format
    if op.isfile(atlas_path):
        atlas_df = pd.read_csv(atlas_path, sep='\t', nrows=5)
        eprint(f'  Atlas created: {atlas_path}')
        eprint(f'  Shape: {atlas_df.shape[0]}+ markers x {atlas_df.shape[1]} columns')
        eprint(f'  Cell types: {list(atlas_df.columns[8:])}')
    else:
        eprint(f'  WARNING: Atlas file not created!')

    return atlas_path


###############################################################################
#                                                                             #
#   Stage 5: Validation Report                                               #
#                                                                             #
###############################################################################

def stage5_validation_report(args, cgi_shore_bed, blocks_path, markers_path,
                             atlas_path):
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

    # Atlas info
    lines.append('--- Atlas ---')
    if op.isfile(atlas_path):
        atlas_df = pd.read_csv(atlas_path, sep='\t')
        lines.append(f'Atlas shape: {atlas_df.shape[0]} markers x '
                     f'{atlas_df.shape[1] - 8} cell types')
        lines.append(f'Atlas file: {atlas_path}')
    else:
        lines.append('  Atlas file not found!')
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
    parser.add_argument('--num-markers', type=int, default=25,
                        help='Target number of markers per cell type '
                             '(default: 25, use 250 for high-resolution)')
    parser.add_argument('--delta-means', type=float, default=0.35,
                        help='Initial delta_means threshold (default: 0.35)')
    parser.add_argument('--min-cpg', type=int, default=3,
                        help='Minimum CpGs per block (default: 3)')
    parser.add_argument('--max-cpg', type=int, default=50,
                        help='Maximum CpGs per block (default: 50)')
    parser.add_argument('--min-bp', type=int, default=50,
                        help='Minimum block length in bp (default: 50)')
    parser.add_argument('--max-bp', type=int, default=5000,
                        help='Maximum block length in bp (default: 5000)')

    # Atlas building
    parser.add_argument('--rlen', type=int, default=4,
                        help='Minimum CpGs per read for homog (default: 4)')

    # Runtime
    parser.add_argument('--threads', type=int, default=4,
                        help='Number of threads (default: 4)')
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

    # Stage 4: Build UXM atlas
    atlas_path = stage4_build_atlas(args, markers_path)

    # Stage 5: Validation report
    stage5_validation_report(args, cgi_shore_bed, blocks_path,
                            markers_path, atlas_path)

    eprint('\n' + '=' * 70)
    eprint('Pipeline complete!')
    eprint(f'  Markers: {markers_path}')
    eprint(f'  Atlas: {atlas_path}')
    eprint(f'  Use with: uxm deconv -a {atlas_path} sample.pat.gz')
    eprint('=' * 70)


if __name__ == '__main__':
    main()
