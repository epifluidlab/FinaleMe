"""Tests for I/O format parsers."""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from finaleme_too.config import MeasurementMode
from finaleme_too.io.marker_regions import MarkerRegions, MarkerRegionsLoader
from finaleme_too.io.methylation_loader import MethylationLoader
from finaleme_too.io.reference_panel import ReferencePanelLoader, load_cpg_index
from finaleme_too.io.sample_sheet import SampleSheet


def test_marker_regions_bed_loader(tmp_path: Path):
    bed_path = tmp_path / "markers.bed"
    bed_path.write_text("chr1\t100\t200\tm1\nchr1\t500\t600\tm2\nchr2\t1000\t1100\tm3\n")
    mr = MarkerRegionsLoader.load(bed_path)
    assert mr.n_markers == 3
    assert list(mr.chrom) == ["chr1", "chr1", "chr2"]
    assert list(mr.start) == [100, 500, 1000]
    assert list(mr.end) == [200, 600, 1100]
    assert mr.marker_name is not None


def test_marker_regions_uxm_atlas_loader(tmp_path: Path):
    atlas_path = tmp_path / "test.atlas"
    atlas_path.write_text(
        "chr\tstart\tend\tstartCpG\tendCpG\ttarget\tname\tdirection\tCellA\tCellB\n"
        "chr1\t100\t200\t1\t5\tCellA\tchr1:100-200\tU\t0.1\t0.9\n"
        "chr1\t500\t600\t10\t15\tCellB\tchr1:500-600\tM\t0.95\t0.05\n"
    )
    mr = MarkerRegionsLoader.load(atlas_path)
    assert mr.n_markers == 2
    assert list(mr.start) == [100, 500]


def test_finaleme_bed_loader(tmp_path: Path):
    """Per-CpG records aggregated to marker regions."""
    bed_path = tmp_path / "sample.prediction.bed.gz"
    rows = [
        # 3 CpGs inside marker 1 (chr1:100-200), 2 inside marker 2 (chr1:500-600)
        ("chr1", 110, 111, 50.0, 5, 10),
        ("chr1", 130, 131, 100.0, 10, 10),
        ("chr1", 150, 151, 0.0, 0, 10),
        ("chr1", 510, 511, 100.0, 10, 10),
        ("chr1", 550, 551, 50.0, 5, 10),
    ]
    with gzip.open(bed_path, "wt") as fh:
        fh.write(
            "#chr\tstart\tend\tmethy_perc_predict\tmethy_count_predict\ttotal_count_predict"
            "\tmethy_perc_obs\tmethy_count_obs\ttotal_count_obs\n"
        )
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\t0\t0\t0\n")

    mr = MarkerRegions(
        chrom=np.array(["chr1", "chr1"], dtype=object),
        start=np.array([100, 500], dtype=np.int64),
        end=np.array([200, 600], dtype=np.int64),
        marker_name=None,
    )
    obs = MethylationLoader.load(
        filepath=bed_path,
        sample_id="s1",
        mode=MeasurementMode.FINALEME,
        marker_regions=mr,
        input_format="finaleme_bed",
    )
    assert int(obs.k[0]) == 5 + 10 + 0
    assert int(obs.n[0]) == 30
    assert int(obs.k[1]) == 10 + 5
    assert int(obs.n[1]) == 20


def test_sample_sheet_validation(tmp_path: Path):
    sheet_path = tmp_path / "sheet.tsv"
    sheet_path.write_text(
        "sample_id\tmethylation_file\tmode\tgroup\n"
        "s1\t/nope/a.bed.gz\tFINALEME\tA\n"
        "s2\t/nope/b.bed.gz\tWGBS\tB\n"
    )
    sheet = SampleSheet.from_tsv(sheet_path)
    assert len(sheet) == 2
    assert sheet.samples[0].mode == MeasurementMode.FINALEME
    assert sheet.samples[1].group == "B"
    assert sheet.groups() == ["A", "B"]


def test_sample_sheet_missing_required_column(tmp_path: Path):
    sheet_path = tmp_path / "bad.tsv"
    sheet_path.write_text("sample_id\tmode\ns1\tWGBS\n")
    with pytest.raises(Exception):
        SampleSheet.from_tsv(sheet_path)


def test_reference_panel_matrix_loader(tmp_path: Path):
    path = tmp_path / "ref.tsv"
    path.write_text(
        "chrom\tstart\tend\tCellA\tCellB\n"
        "chr1\t100\t200\t0.1\t0.9\n"
        "chr1\t500\t600\t0.95\t0.05\n"
    )
    ref = ReferencePanelLoader.load_matrix(path)
    assert ref.n_markers == 2
    assert ref.cell_types == ["CellA", "CellB"]
    assert ref.methylation.shape == (2, 2)
    np.testing.assert_allclose(ref.methylation[0], [0.1, 0.9], atol=1e-6)


def test_load_cpg_index(tmp_path: Path):
    path = tmp_path / "cpg.bed"
    path.write_text("chr1\t100\t101\nchr1\t200\t201\nchr2\t50\t51\n")
    idx = load_cpg_index(path)
    assert idx["total_sites"] == 3
    assert "chr1" in idx["chr_positions"]
    assert "chr2" in idx["chr_positions"]
    assert idx["chr_offsets"]["chr1"] == 0
    assert idx["chr_offsets"]["chr2"] == 2
    np.testing.assert_array_equal(idx["chr_positions"]["chr1"], [100, 200])
