"""Marker regions loader: BED files and UXM atlas TSVs.

UXM atlas files contain U/M ratio patterns, NOT methylation levels. Only the
marker region coordinates are extracted; a separate reference panel with
actual methylation beta values is always required for TOO deconvolution.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from finaleme_too.exceptions import InvalidMarkerRegionsError


@dataclass(frozen=True)
class MarkerRegions:
    """Genomic coordinates of M marker regions used for deconvolution."""

    chrom: np.ndarray  # shape (M,) dtype=object (str)
    start: np.ndarray  # shape (M,) dtype=int64
    end: np.ndarray  # shape (M,) dtype=int64
    marker_name: np.ndarray | None = None  # shape (M,) dtype=object or None

    def __len__(self) -> int:
        return len(self.chrom)

    @property
    def n_markers(self) -> int:
        return len(self.chrom)


class MarkerRegionsLoader:
    """Static loader for marker region files."""

    @staticmethod
    def load(filepath: str | Path, marker_format: str = "auto") -> MarkerRegions:
        path = Path(filepath)
        if not path.exists():
            raise InvalidMarkerRegionsError(f"Marker regions file not found: {path}")

        if marker_format == "auto":
            marker_format = MarkerRegionsLoader._detect_format(path)

        if marker_format == "uxm_atlas":
            return MarkerRegionsLoader._parse_uxm_atlas(path)
        if marker_format == "bed":
            return MarkerRegionsLoader._parse_bed(path)
        raise InvalidMarkerRegionsError(f"Unknown marker_format: {marker_format}")

    @staticmethod
    def _detect_format(path: Path) -> str:
        name = path.name.lower()
        if name.endswith(".atlas.gz") or name.endswith(".atlas") or name.endswith(".atlas.tsv"):
            return "uxm_atlas"
        # If the first non-comment line has 8+ columns and column 4 is an integer
        # (startCpG) we assume it's a UXM-style atlas. Otherwise BED.
        try:
            opener = gzip.open if name.endswith(".gz") else open
            with opener(path, "rt") as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split("\t")
                    if len(parts) >= 8:
                        try:
                            int(parts[3])
                            int(parts[4])
                            return "uxm_atlas"
                        except (ValueError, IndexError):
                            pass
                    return "bed"
        except OSError:
            pass
        return "bed"

    @staticmethod
    def _parse_bed(path: Path) -> MarkerRegions:
        chroms: list[str] = []
        starts: list[int] = []
        ends: list[int] = []
        names: list[str] = []
        opener = gzip.open if path.name.endswith(".gz") else open
        with opener(path, "rt") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("track"):
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    continue
                chroms.append(parts[0])
                starts.append(start)
                ends.append(end)
                names.append(parts[3] if len(parts) >= 4 else "")
        if not chroms:
            raise InvalidMarkerRegionsError(f"No valid BED records found in {path}")
        marker_name_arr: np.ndarray | None
        if any(n != "" for n in names):
            marker_name_arr = np.array(names, dtype=object)
        else:
            marker_name_arr = None
        return MarkerRegions(
            chrom=np.array(chroms, dtype=object),
            start=np.array(starts, dtype=np.int64),
            end=np.array(ends, dtype=np.int64),
            marker_name=marker_name_arr,
        )

    @staticmethod
    def _parse_uxm_atlas(path: Path) -> MarkerRegions:
        """Extract only chrom/start/end (and optionally name) from a UXM atlas TSV.

        UXM atlas columns (typical):
            chr  start  end  startCpG  endCpG  target  region  direction  CellType1 ...
        """
        chroms: list[str] = []
        starts: list[int] = []
        ends: list[int] = []
        names: list[str] = []
        opener = gzip.open if path.name.endswith(".gz") else open
        with opener(path, "rt") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("#") or line.lower().startswith("chr\t"):
                    continue
                parts = line.split("\t")
                if len(parts) < 5:
                    continue
                # Skip header row written without leading '#' (e.g. "chr\tstart\t..."):
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    continue
                chroms.append(parts[0])
                starts.append(start)
                ends.append(end)
                if len(parts) >= 7:
                    names.append(parts[6])
                else:
                    names.append("")
        if not chroms:
            raise InvalidMarkerRegionsError(f"No valid UXM atlas records found in {path}")
        marker_name_arr: np.ndarray | None
        if any(n != "" for n in names):
            marker_name_arr = np.array(names, dtype=object)
        else:
            marker_name_arr = None
        return MarkerRegions(
            chrom=np.array(chroms, dtype=object),
            start=np.array(starts, dtype=np.int64),
            end=np.array(ends, dtype=np.int64),
            marker_name=marker_name_arr,
        )


__all__ = ["MarkerRegions", "MarkerRegionsLoader"]
