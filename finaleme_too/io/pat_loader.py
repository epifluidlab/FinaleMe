"""Loader for wgbstools/FinaleMe ``.pat.gz`` fragment pattern files.

Pat format (bgzipped / plain TSV):
    chr  startCpG  pattern  count

where ``startCpG`` is the 1-based global CpG index of the first CpG in the
fragment pattern, ``pattern`` is a string of ``C`` (methylated) and ``T``
(unmethylated) calls — one character per consecutive CpG in the fragment —
and ``count`` is the multiplicity of that (chrom, startCpG, pattern) tuple.

This loader translates a .pat.gz file into a list of ``Fragment`` objects
that the :class:`~finaleme_too.core.fragment_likelihood.FragmentLevelDeconvolver`
can consume. Marker regions are provided to restrict the fragments to the
subset we care about — typically the top-N most discriminative markers
selected by :class:`~finaleme_too.preprocessing.marker_selection.BalancedMarkerSelector`.

The mapping from *marker region* to *reference-matrix row index* is a
straightforward lookup: each marker region becomes a block of consecutive
CpG positions, and the fragment-level deconvolver's reference matrix is
indexed *per marker* (not per CpG). So a fragment that covers CpGs inside
marker ``i`` contributes a single entry with ``cpg_indices[..] == i``.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np

from finaleme_too.core.fragment_likelihood import Fragment
from finaleme_too.exceptions import InvalidInputFormatError
from finaleme_too.io.marker_regions import MarkerRegions


def load_fragments_from_pat(
    pat_path: str | Path,
    marker_regions: MarkerRegions,
    cpg_index: dict,
    max_fragments: int | None = None,
) -> list[Fragment]:
    """Parse a .pat.gz file into fragments aligned to marker regions.

    Parameters
    ----------
    pat_path : path-like
        Tab-separated file ``chr startCpG pattern count`` (optionally gzipped).
    marker_regions : MarkerRegions
        Marker regions used by the deconvolver. Each marker becomes a
        single index into the reference matrix.
    cpg_index : dict
        Output of :func:`finaleme_too.io.reference_panel.load_cpg_index`.
        Maps chrom → sorted CpG positions + per-chrom 1-based offsets.
    max_fragments : int or None
        If set, stop parsing after this many fragments (useful for tests).

    Returns
    -------
    list[Fragment]
        Each Fragment's ``cpg_indices`` hold *marker* indices (not CpG
        indices), and ``methylated`` contains the per-CpG call aggregated
        into that marker.
    """
    path = Path(pat_path)
    if not path.exists():
        raise InvalidInputFormatError(f".pat.gz file not found: {path}")

    # Build a per-marker CpG range lookup: for each marker i, determine the
    # global 1-based CpG index range [lo, hi). A fragment's CpG indices are
    # matched against this range, then collapsed into marker indices.
    marker_ranges: list[tuple[str, int, int, int]] = []  # (chr, lo, hi, marker_idx)
    chr_positions = cpg_index["chr_positions"]
    chr_offsets = cpg_index["chr_offsets"]
    for mi in range(marker_regions.n_markers):
        chrom = str(marker_regions.chrom[mi])
        start = int(marker_regions.start[mi])
        end = int(marker_regions.end[mi])
        positions = chr_positions.get(chrom)
        offset = chr_offsets.get(chrom)
        if positions is None or offset is None:
            continue
        lo = int(np.searchsorted(positions, start, side="left"))
        hi = int(np.searchsorted(positions, end, side="left"))
        if hi <= lo:
            continue
        marker_ranges.append((chrom, offset + lo + 1, offset + hi + 1, mi))  # 1-based

    # Group ranges by chromosome for fast lookup during fragment scanning
    ranges_by_chrom: dict[str, list[tuple[int, int, int]]] = {}
    for chrom, lo, hi, mi in marker_ranges:
        ranges_by_chrom.setdefault(chrom, []).append((lo, hi, mi))
    for chrom in ranges_by_chrom:
        ranges_by_chrom[chrom].sort()

    def _find_marker(chrom: str, cpg_pos: int) -> int | None:
        entries = ranges_by_chrom.get(chrom)
        if not entries:
            return None
        # Binary search on sorted (lo, hi, marker_idx)
        los = np.asarray([e[0] for e in entries], dtype=np.int64)
        idx = int(np.searchsorted(los, cpg_pos, side="right") - 1)
        if idx < 0:
            return None
        lo, hi, mi = entries[idx]
        if lo <= cpg_pos < hi:
            return mi
        return None

    fragments: list[Fragment] = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                start_cpg = int(parts[1])
                pattern = parts[2]
                count = int(parts[3]) if len(parts) >= 4 else 1
            except ValueError:
                continue
            if not pattern:
                continue
            # Map each CpG in the pattern to a marker index; drop CpGs
            # that do not fall inside any marker.
            marker_idx_list: list[int] = []
            meth_calls: list[int] = []
            for offset, ch in enumerate(pattern):
                if ch not in ("C", "T"):
                    continue
                cpg_pos = start_cpg + offset
                mi = _find_marker(chrom, cpg_pos)
                if mi is None:
                    continue
                marker_idx_list.append(mi)
                meth_calls.append(1 if ch == "C" else 0)
            if not marker_idx_list:
                continue
            # Emit one Fragment per unique (marker, count) repetition. We
            # expand ``count`` by appending the same Fragment ``count`` times
            # because the EM algorithm sums responsibilities per-fragment.
            frag = Fragment(
                cpg_indices=np.asarray(marker_idx_list, dtype=np.int64),
                methylated=np.asarray(meth_calls, dtype=np.uint8),
            )
            for _ in range(count):
                fragments.append(frag)
                if max_fragments is not None and len(fragments) >= max_fragments:
                    return fragments
    return fragments


__all__ = ["load_fragments_from_pat"]
