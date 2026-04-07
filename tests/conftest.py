"""Shared pytest fixtures."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from finaleme_too.config import MeasurementMode
from finaleme_too.io.marker_regions import MarkerRegions
from finaleme_too.io.methylation_loader import MarkerObservations
from finaleme_too.io.reference_panel import ReferencePanel


# Real test data lives outside the repo; only used by integration tests.
REAL_TEST_DATA_DIR = Path("/Users/yaping/Documents/workspace/projects/ccinference/test_bam/test_uxm")


@pytest.fixture
def real_test_data_dir() -> Path:
    if not REAL_TEST_DATA_DIR.exists():
        pytest.skip(f"Real test data not available at {REAL_TEST_DATA_DIR}")
    return REAL_TEST_DATA_DIR


@pytest.fixture
def synthetic_marker_regions() -> MarkerRegions:
    """50 marker regions on chr1, 1kb apart."""
    n = 50
    return MarkerRegions(
        chrom=np.array(["chr1"] * n, dtype=object),
        start=np.array([1000 + i * 2000 for i in range(n)], dtype=np.int64),
        end=np.array([1500 + i * 2000 for i in range(n)], dtype=np.int64),
        marker_name=None,
    )


@pytest.fixture
def synthetic_reference(synthetic_marker_regions) -> ReferencePanel:
    """5 cell types where each owns 10 markers (rest are 0.5)."""
    n_markers = 50
    n_cell_types = 5
    rng = np.random.default_rng(42)
    methy = np.full((n_markers, n_cell_types), 0.5, dtype=np.float32)
    for j in range(n_cell_types):
        idx = slice(j * 10, (j + 1) * 10)
        # Cell type j: U pattern at its 10 markers (low methylation),
        # all other cell types have high methylation at those markers.
        methy[idx, j] = rng.uniform(0.0, 0.05, size=10).astype(np.float32)
        for jj in range(n_cell_types):
            if jj == j:
                continue
            methy[idx, jj] = rng.uniform(0.85, 1.0, size=10).astype(np.float32)
    return ReferencePanel(
        chrom=synthetic_marker_regions.chrom,
        start=synthetic_marker_regions.start,
        end=synthetic_marker_regions.end,
        cell_types=[f"CellType{i+1}" for i in range(n_cell_types)],
        methylation=methy,
        coverage=np.full((n_markers, n_cell_types), 50, dtype=np.int32),
    )


@pytest.fixture
def synthetic_observations_pure_celltype(
    synthetic_marker_regions, synthetic_reference
) -> MarkerObservations:
    """Self-deconv: methylation matches CellType1's reference profile exactly.

    With 30 reads per marker, the deconvolver should recover w_CellType1 ≈ 1.0.
    """
    n_markers = synthetic_marker_regions.n_markers
    rng = np.random.default_rng(123)
    n = np.full(n_markers, 30, dtype=np.int32)
    p = synthetic_reference.methylation[:, 0].astype(np.float64)
    k = rng.binomial(n.astype(np.int64), p).astype(np.int32)
    return MarkerObservations(
        sample_id="self_celltype1",
        chrom=synthetic_marker_regions.chrom,
        start=synthetic_marker_regions.start,
        end=synthetic_marker_regions.end,
        k=k,
        n=n,
        predicted_beta=None,
        mode=MeasurementMode.WGBS,
    )
