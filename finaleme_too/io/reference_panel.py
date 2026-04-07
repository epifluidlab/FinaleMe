"""Reference panel loader.

Two supported input modes (architecture §3.3):
  - Mode A: matrix file (TSV with chrom/start/end + one column per cell type)
  - Mode B: beta file list + groups CSV (compatible with BetaValueDeconvolution)
"""

from __future__ import annotations

import gzip
from collections import defaultdict
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import pandas as pd

from finaleme_too.exceptions import InvalidReferencePanelError
from finaleme_too.io.marker_regions import MarkerRegions


@dataclass(frozen=True)
class ReferencePanel:
    """Reference methylation matrix.

    methylation has shape (M, K) where M = number of marker regions,
    K = number of cell types. Values are methylation density (0..1).
    Optional coverage matrix has the same shape; integer counts of total
    reads per (marker, cell_type) used for reference uncertainty weighting.
    """

    chrom: np.ndarray
    start: np.ndarray
    end: np.ndarray
    cell_types: list[str]
    methylation: np.ndarray  # float32, (M, K)
    coverage: np.ndarray | None = None  # int32, (M, K) — optional

    @property
    def n_markers(self) -> int:
        return self.methylation.shape[0]

    @property
    def n_cell_types(self) -> int:
        return self.methylation.shape[1]

    def to_marker_regions(self) -> MarkerRegions:
        return MarkerRegions(
            chrom=self.chrom,
            start=self.start,
            end=self.end,
            marker_name=None,
        )


# ---------------------------------------------------------------------------
# CpG index loading (mirrors Java BetaValueDeconvolution.CpgIndex)
# ---------------------------------------------------------------------------


def load_cpg_index(path: str | Path) -> dict:
    """Load a CpG index BED file (e.g. data/CpG_index.hg19.bed.gz).

    Returns a dict with keys:
        chr_positions: dict[str -> np.ndarray] (sorted positions per chromosome)
        chr_offsets:   dict[str -> int]        (cumulative position offset)
        total_sites:   int

    1-based global CpG index of position p on chromosome c is computed as:

        offset[c] + bisect_left(positions[c], p) + 1
    """
    p = Path(path)
    if not p.exists():
        raise InvalidReferencePanelError(f"CpG index not found: {p}")

    chr_to_positions: dict[str, list[int]] = defaultdict(list)
    opener = gzip.open if p.name.endswith(".gz") else open
    with opener(p, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            try:
                start = int(parts[1])
            except ValueError:
                continue
            chr_to_positions[parts[0]].append(start)

    if not chr_to_positions:
        raise InvalidReferencePanelError(f"Empty CpG index: {p}")

    chr_positions: dict[str, np.ndarray] = {}
    chr_offsets: dict[str, int] = {}
    total = 0
    for chrom, positions in chr_to_positions.items():
        arr = np.unique(np.asarray(positions, dtype=np.int64))
        chr_offsets[chrom] = total
        chr_positions[chrom] = arr
        total += int(arr.size)
    return {
        "chr_positions": chr_positions,
        "chr_offsets": chr_offsets,
        "total_sites": total,
    }


# ---------------------------------------------------------------------------
# Loader
# ---------------------------------------------------------------------------


class ReferencePanelLoader:
    """Static loader for reference panels."""

    @staticmethod
    def load_matrix(filepath: str | Path) -> ReferencePanel:
        """Load a TSV reference panel: chrom start end CellType1 CellType2 ..."""
        path = Path(filepath)
        if not path.exists():
            raise InvalidReferencePanelError(f"Reference panel not found: {path}")
        opener = gzip.open if path.name.endswith(".gz") else open
        with opener(path, "rt") as fh:
            df = pd.read_csv(fh, sep="\t", comment="#")
        if df.shape[1] < 4:
            raise InvalidReferencePanelError(
                f"Reference panel must have >=4 columns; got {df.shape[1]}"
            )
        # Rename leading columns case-insensitively
        rename = {}
        for col, target in zip(df.columns[:3], ["chrom", "start", "end"]):
            rename[col] = target
        df = df.rename(columns=rename)
        cell_types = [str(c) for c in df.columns[3:]]
        return ReferencePanel(
            chrom=df["chrom"].astype(str).to_numpy(),
            start=df["start"].astype(np.int64).to_numpy(),
            end=df["end"].astype(np.int64).to_numpy(),
            cell_types=cell_types,
            methylation=df[cell_types].to_numpy(dtype=np.float32),
            coverage=None,
        )

    @staticmethod
    def load_beta_list(
        ref_betas_arg: str,
        groups_file: str | Path,
        cpg_index_path: str | Path,
        marker_regions: MarkerRegions,
        replicate_mode: str = "aggregate",
    ) -> ReferencePanel:
        """Load reference panel from a list of .beta files + groups CSV.

        Mirrors BetaValueDeconvolution._loadReferenceMethylation. The
        ``ref_betas_arg`` may be either:
            * a comma-separated list of .beta paths, OR
            * a path to a .txt file with one .beta path per line.

        Groups CSV must have a header with at least ``name`` and ``group``
        columns. ``name`` is matched to beta file basenames (stripping the
        ``.beta`` extension and any ``.pat.gz`` artifact).
        """
        beta_paths = _parse_ref_betas_arg(ref_betas_arg)
        groups_df = _load_groups_file(groups_file)
        cpg_index = load_cpg_index(cpg_index_path)

        # Strip extensions on the groups names
        groups_df = groups_df.copy()
        groups_df["name_stripped"] = (
            groups_df["name"].astype(str).str.replace(r"\.(pat\.gz|beta|pat)$", "", regex=True)
        )
        sample_to_group = dict(zip(groups_df["name_stripped"], groups_df["group"]))

        cell_types = sorted(set(groups_df["group"]))
        ct_to_index = {c: i for i, c in enumerate(cell_types)}
        n_markers = marker_regions.n_markers
        n_ct = len(cell_types)

        # For aggregate mode we sum methylated and total counts per group
        agg_methy = np.zeros((n_markers, n_ct), dtype=np.int64)
        agg_total = np.zeros((n_markers, n_ct), dtype=np.int64)
        # For average mode we collect per-sample ratios
        ratio_sum = np.zeros((n_markers, n_ct), dtype=np.float64)
        ratio_count = np.zeros((n_markers, n_ct), dtype=np.int64)

        for beta_path in beta_paths:
            sample_name = _beta_basename_to_sample_name(beta_path)
            group = sample_to_group.get(sample_name)
            if group is None:
                # Try without _ to - replacement etc — fall back: skip silently
                continue
            ci = ct_to_index[group]
            counts = _load_beta_file_to_markers(beta_path, marker_regions, cpg_index)
            mk = counts[:, 0]
            mn = counts[:, 1]
            if replicate_mode == "aggregate":
                agg_methy[:, ci] += mk
                agg_total[:, ci] += mn
            else:  # "average"
                with np.errstate(invalid="ignore", divide="ignore"):
                    ratios = np.where(mn > 0, mk / np.maximum(mn, 1), np.nan)
                valid = mn > 0
                ratio_sum[valid, ci] += ratios[valid]
                ratio_count[valid, ci] += 1

        if replicate_mode == "aggregate":
            methylation = np.full((n_markers, n_ct), np.nan, dtype=np.float32)
            with np.errstate(invalid="ignore", divide="ignore"):
                ok = agg_total > 0
                methylation[ok] = (agg_methy[ok] / agg_total[ok]).astype(np.float32)
            coverage = agg_total.astype(np.int32)
        else:
            with np.errstate(invalid="ignore", divide="ignore"):
                methylation = np.where(
                    ratio_count > 0, ratio_sum / np.maximum(ratio_count, 1), np.nan
                ).astype(np.float32)
            coverage = ratio_count.astype(np.int32)

        return ReferencePanel(
            chrom=marker_regions.chrom,
            start=marker_regions.start,
            end=marker_regions.end,
            cell_types=cell_types,
            methylation=methylation,
            coverage=coverage,
        )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _parse_ref_betas_arg(arg: str) -> list[Path]:
    """Accept a comma-separated string OR a .txt file with one path per line."""
    if "," in arg or not arg.endswith(".txt"):
        return [Path(p.strip()) for p in arg.split(",") if p.strip()]
    paths: list[Path] = []
    with open(arg) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                paths.append(Path(line))
    return paths


def _load_groups_file(path: str | Path) -> pd.DataFrame:
    """Load groups CSV and validate columns."""
    df = pd.read_csv(path, comment="#")
    cols_lower = {c.lower(): c for c in df.columns}
    if "name" not in cols_lower or "group" not in cols_lower:
        raise InvalidReferencePanelError(
            f"Groups file must have 'name' and 'group' columns; got {list(df.columns)}"
        )
    df = df.rename(columns={cols_lower["name"]: "name", cols_lower["group"]: "group"})
    return df[["name", "group"]]


def _beta_basename_to_sample_name(path: Path) -> str:
    name = path.name
    for ext in (".beta", ".lbeta"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def _load_beta_file_to_markers(
    beta_path: Path, marker_regions: MarkerRegions, cpg_index: dict
) -> np.ndarray:
    """Read a binary .beta file and aggregate to per-marker (methy, total) counts.

    Returns int64 array of shape (n_markers, 2).
    """
    with open(beta_path, "rb") as fh:
        data = np.frombuffer(fh.read(), dtype=np.uint8)
    if data.size % 2 != 0:
        raise InvalidReferencePanelError(
            f"Beta file size not a multiple of 2: {beta_path}"
        )
    per_cpg = data.reshape((-1, 2)).astype(np.int64)
    chr_positions = cpg_index["chr_positions"]
    chr_offsets = cpg_index["chr_offsets"]
    total_sites = cpg_index["total_sites"]

    n_markers = marker_regions.n_markers
    out = np.zeros((n_markers, 2), dtype=np.int64)
    for mi in range(n_markers):
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
        global_lo = offset + lo
        global_hi = min(offset + hi, total_sites, per_cpg.shape[0])
        if global_lo >= global_hi:
            continue
        block = per_cpg[global_lo:global_hi]
        out[mi, 0] = int(block[:, 0].sum())
        out[mi, 1] = int(block[:, 1].sum())
    return out


__all__ = ["ReferencePanel", "ReferencePanelLoader", "load_cpg_index"]
