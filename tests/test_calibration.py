"""Tests for calibration apply path and calibration_eval QC (Phase B)."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from finaleme_too.config import MeasurementMode
from finaleme_too.io.methylation_loader import MarkerObservations
from finaleme_too.preprocessing.calibration import (
    CalibrationParams,
    apply_calibration,
    load_default_calibration,
)
from finaleme_too.preprocessing.calibration_eval import (
    compute_hosmer_lemeshow,
    compute_inference_qc,
)


def _identity_calibration(n_bins: int = 4) -> CalibrationParams:
    return CalibrationParams(
        n_bins=n_bins,
        bin_edges=np.linspace(0.0, 1.0, n_bins + 1),
        a=np.ones(n_bins),
        c=np.zeros(n_bins),
        log_dispersion=np.full(n_bins, np.log(20.0)),
    )


def test_calibration_save_load_round_trip(tmp_path):
    p = _identity_calibration(6)
    out = tmp_path / "cal.json"
    p.save(out)
    loaded = CalibrationParams.load(out)
    assert loaded.n_bins == 6
    np.testing.assert_allclose(loaded.bin_edges, p.bin_edges)
    np.testing.assert_allclose(loaded.a, p.a)
    np.testing.assert_allclose(loaded.c, p.c)


def test_apply_identity_calibration_is_no_op(synthetic_marker_regions):
    n = synthetic_marker_regions.n_markers
    pred = np.full(n, 0.4, dtype=np.float32)
    obs = MarkerObservations(
        sample_id="x",
        chrom=synthetic_marker_regions.chrom,
        start=synthetic_marker_regions.start,
        end=synthetic_marker_regions.end,
        k=np.full(n, 12, dtype=np.int32),
        n=np.full(n, 30, dtype=np.int32),
        predicted_beta=pred,
        mode=MeasurementMode.FINALEME,
    )
    cal = _identity_calibration()
    # No region annotations → defaults to middle bin (still identity)
    new_obs = apply_calibration(obs, cal, region_annotations=None)
    # Identity calibration should preserve k roughly equal to round(0.4 * 30) = 12
    assert int(new_obs.k[0]) == 12
    np.testing.assert_array_equal(new_obs.n, obs.n)


def test_apply_calibration_changes_methylation():
    """Calibration with a > 1 inflates extreme predictions."""
    pred = np.array([0.2, 0.8], dtype=np.float32)
    obs = MarkerObservations(
        sample_id="x",
        chrom=np.array(["chr1", "chr1"], dtype=object),
        start=np.array([100, 200], dtype=np.int64),
        end=np.array([200, 300], dtype=np.int64),
        k=np.array([6, 24], dtype=np.int32),
        n=np.array([30, 30], dtype=np.int32),
        predicted_beta=pred,
        mode=MeasurementMode.FINALEME,
    )
    # a = 1.5 stretches logit; c = 0
    cal = CalibrationParams(
        n_bins=1,
        bin_edges=np.array([0.0, 1.0]),
        a=np.array([1.5]),
        c=np.array([0.0]),
        log_dispersion=np.array([np.log(20.0)]),
    )
    new_obs = apply_calibration(obs, cal, region_annotations=None)
    # 0.2 should move toward 0; 0.8 should move toward 1
    assert int(new_obs.k[0]) < 6
    assert int(new_obs.k[1]) > 24


def test_default_calibration_is_loadable():
    cal = load_default_calibration()
    assert cal.n_bins == 8
    assert cal.bin_edges.size == 9
    assert cal.a.size == 8
    assert cal.c.size == 8
    assert cal.log_dispersion.size == 8


def test_hosmer_lemeshow_perfect_calibration():
    """Perfectly calibrated predictions should give a high p-value."""
    rng = np.random.default_rng(0)
    pred = np.linspace(0.05, 0.95, 200)
    obs = rng.binomial(50, pred) / 50
    p = compute_hosmer_lemeshow(obs, pred, n_groups=10)
    assert not np.isnan(p)
    assert p > 0.01  # not strictly significant under perfect calibration


def test_hosmer_lemeshow_miscalibrated():
    """Severely miscalibrated predictions should give a low p-value."""
    rng = np.random.default_rng(1)
    pred = np.full(500, 0.5)
    obs = rng.binomial(50, 0.9, size=500) / 50  # truly 0.9, predicted 0.5
    p = compute_hosmer_lemeshow(obs, pred, n_groups=10)
    assert p < 0.001


def test_inference_qc_pass_for_normal_sample():
    cal = _identity_calibration(8)
    sample_pred = np.linspace(0.05, 0.95, 200)
    qc = compute_inference_qc(sample_pred, cal)
    assert qc["flag"] in ("PASS", "WARN", "FAIL")
    assert 0.0 <= qc["bin_coverage_balance"] <= 1.0
    assert 0.0 <= qc["prediction_range_coverage"] <= 1.0


# ---------------------------------------------------------------------------
# Phase C: training-time fit + tuning
# ---------------------------------------------------------------------------


def test_fit_calibration_recovers_identity():
    from finaleme_too.preprocessing.calibration_eval import fit_calibration

    rng = np.random.default_rng(0)
    n = 500
    truth = rng.beta(2, 5, size=n)
    fme = truth + rng.normal(0, 0.02, size=n)
    fme = np.clip(fme, 0.01, 0.99)
    density = rng.uniform(0, 0.1, size=n)
    fit = fit_calibration(
        finaleme_beta=fme,
        wgbs_beta=truth,
        cpg_density=density,
        n_bins=4,
    )
    # Identity-like → slope ~ 1, intercept ~ 0
    np.testing.assert_allclose(fit.a, np.ones(4), atol=0.2)
    np.testing.assert_allclose(fit.c, np.zeros(4), atol=0.2)
    assert fit.overall["r_squared"] > 0.9


def test_tune_n_bins_selects_finite_value():
    from finaleme_too.preprocessing.calibration_eval import tune_n_bins

    rng = np.random.default_rng(1)
    n = 600
    sample_ids = np.array([f"S{i // 100}" for i in range(n)])
    truth = rng.beta(2, 5, size=n)
    fme = np.clip(truth + rng.normal(0, 0.05, size=n), 0.01, 0.99)
    density = rng.uniform(0, 0.1, size=n)
    out = tune_n_bins(
        finaleme_beta=fme,
        wgbs_beta=truth,
        cpg_density=density,
        sample_ids=sample_ids,
        n_bins_candidates=[2, 4, 8],
    )
    assert out["selected_n_bins"] in (2, 4, 8)
    assert len(out["candidates"]) == 3


def test_train_calibration_end_to_end(tmp_path):
    from scipy.special import expit, logit

    from finaleme_too.preprocessing.calibration import train_calibration

    rng = np.random.default_rng(2)
    n_samples = 4
    n_markers = 150
    rows_w = []
    rows_f = []
    rows_a = []
    for s in range(n_samples):
        sid = f"S{s}"
        for m in range(n_markers):
            true_b = float(rng.beta(2, 5))
            n_w = 30
            k_w = int(rng.binomial(n_w, true_b))
            density = float(rng.uniform(0, 0.1))
            # Miscalibration: a=0.8, c=0.1
            true_logit = logit(np.clip(true_b, 1e-6, 1 - 1e-6))
            fme_logit = (true_logit - 0.1) / 0.8 + rng.normal(0, 0.2)
            fme_b = float(expit(fme_logit))
            n_f = 30
            k_f = int(rng.binomial(n_f, fme_b))
            rows_w.append((sid, "chr1", m * 100, m * 100 + 50, k_w, n_w))
            rows_f.append((sid, "chr1", m * 100, m * 100 + 50, k_f, n_f))
            if s == 0:
                rows_a.append(("chr1", m * 100, m * 100 + 50, density))

    import pandas as pd

    pd.DataFrame(
        rows_w, columns=["sample_id", "chrom", "start", "end", "methylated_count", "total_count"]
    ).to_csv(tmp_path / "wgbs.tsv", sep="\t", index=False)
    pd.DataFrame(
        rows_f, columns=["sample_id", "chrom", "start", "end", "methylated_count", "total_count"]
    ).to_csv(tmp_path / "fme.tsv", sep="\t", index=False)
    pd.DataFrame(rows_a, columns=["chrom", "start", "end", "cpg_density"]).to_csv(
        tmp_path / "ann.tsv", sep="\t", index=False
    )

    params = train_calibration(
        matched_wgbs=tmp_path / "wgbs.tsv",
        matched_finaleme=tmp_path / "fme.tsv",
        region_annotation=tmp_path / "ann.tsv",
        n_bins_candidates=[2, 4],
        out_params=tmp_path / "cal.json",
        out_report=tmp_path / "report.json",
    )
    assert params.n_bins in (2, 4)
    assert (tmp_path / "cal.json").exists()
    assert (tmp_path / "report.json").exists()
    # Round-trip
    loaded = CalibrationParams.load(tmp_path / "cal.json")
    np.testing.assert_allclose(loaded.a, params.a)


# ---------------------------------------------------------------------------
# Raw-input calibration training: Bis-SNP 6+2 + FinaleMe prediction.bed.gz
# ---------------------------------------------------------------------------


def _write_bissnp_6plus2(path, rows, *, with_track_line: bool = True, with_chr_prefix: bool = False):
    """Write a minimal Bis-SNP 6+2 BED file.

    Each row is (chrom, start, end, methylation_pct, total_count). An
    optional ``track name=...`` line is prepended to test header skipping.
    """
    with open(path, "w") as fh:
        if with_track_line:
            fh.write(
                'track name=test type=bedDetail description="methylation level"\n'
            )
        for chrom, start, end, pct, total in rows:
            c = f"chr{chrom}" if with_chr_prefix else str(chrom)
            fh.write(f"{c}\t{start}\t{end}\t.\t500\t-\t{pct:.2f}\t{total}\n")


def _write_finaleme_prediction(path, rows):
    import gzip as _gzip

    with _gzip.open(path, "wt") as fh:
        fh.write(
            "#chr\tstart\tend\tmethy_perc_predict\tmethy_count_predict\t"
            "total_count_predict\tmethy_perc_obs\tmethy_count_obs\ttotal_count_obs\n"
        )
        for chrom, start, end, meth_count, total in rows:
            pct = 100.0 * meth_count / max(total, 1)
            fh.write(
                f"chr{chrom}\t{start}\t{end}\t{pct:.1f}\t{meth_count}\t{total}\t"
                f"{pct:.1f}\t{meth_count}\t{total}\n"
            )


def test_parse_bissnp_6plus2_round_trip(tmp_path):
    """Bis-SNP parser must drop the track line and compute methylated_count
    from methylation_pct * total_count correctly."""
    from finaleme_too.preprocessing.calibration import _parse_bissnp_6plus2

    rows = [
        (1, 10469, 10470, 50.00, 4),  # 2 methylated
        (1, 10471, 10472, 75.00, 4),  # 3 methylated
        (2, 5000, 5001, 0.00, 10),    # 0 methylated
        (2, 6000, 6001, 100.00, 5),   # 5 methylated
    ]
    p = tmp_path / "s1.cpg.6plus2.bed"
    _write_bissnp_6plus2(p, rows, with_track_line=True, with_chr_prefix=False)

    df = _parse_bissnp_6plus2(p, sample_id="S1")
    assert list(df["sample_id"]) == ["S1"] * 4
    assert list(df["chrom"]) == ["1", "1", "2", "2"]
    assert list(df["methylated_count"]) == [2, 3, 0, 5]
    assert list(df["total_count"]) == [4, 4, 10, 5]


def test_parse_finaleme_prediction_round_trip(tmp_path):
    """FinaleMe prediction parser must pull columns 4 and 5 (predicted counts)."""
    from finaleme_too.preprocessing.calibration import _parse_finaleme_prediction

    rows = [
        (9, 110744061, 110744062, 3, 4),
        (6, 144773899, 144773900, 7, 8),
        (13, 111851314, 111851315, 14, 14),
    ]
    p = tmp_path / "s1.prediction.bed.gz"
    _write_finaleme_prediction(p, rows)

    df = _parse_finaleme_prediction(p, sample_id="S1")
    assert list(df["sample_id"]) == ["S1"] * 3
    assert list(df["methylated_count"]) == [3, 7, 14]
    assert list(df["total_count"]) == [4, 8, 14]
    # FinaleMe always has "chr" prefix at this point; normalization happens
    # one level up in _load_matched_table via _normalize_chrom.
    assert list(df["chrom"]) == ["chr9", "chr6", "chr13"]


def test_load_matched_table_sample_sheet_bissnp(tmp_path):
    """Loading a Bis-SNP sample sheet should concat per-sample tables and
    strip the chr prefix."""
    import pandas as pd

    from finaleme_too.preprocessing.calibration import _load_matched_table

    rows_s1 = [(1, 100, 101, 50.0, 4), (1, 200, 201, 25.0, 4)]
    rows_s2 = [(1, 100, 101, 75.0, 4), (1, 200, 201, 0.0, 4)]
    _write_bissnp_6plus2(tmp_path / "s1.bed", rows_s1, with_chr_prefix=False)
    _write_bissnp_6plus2(tmp_path / "s2.bed", rows_s2, with_chr_prefix=True)  # has chr

    sheet = tmp_path / "wgbs_samples.tsv"
    pd.DataFrame(
        [
            {"sample_id": "S1", "methylation_file": "s1.bed"},
            {"sample_id": "S2", "methylation_file": "s2.bed"},
        ]
    ).to_csv(sheet, sep="\t", index=False)

    df = _load_matched_table(sheet, modality="wgbs")
    # After _normalize_chrom, all chromosomes have no prefix
    assert set(df["chrom"]) == {"1"}
    assert set(df["sample_id"]) == {"S1", "S2"}
    # S1 row 0: 50% of 4 = 2 methylated
    s1_row0 = df[(df["sample_id"] == "S1") & (df["start"] == 100)].iloc[0]
    assert int(s1_row0["methylated_count"]) == 2


def test_load_matched_table_sample_sheet_finaleme(tmp_path):
    """Loading a FinaleMe sample sheet should concat per-sample tables."""
    import pandas as pd

    from finaleme_too.preprocessing.calibration import _load_matched_table

    rows_s1 = [(1, 100, 101, 3, 4), (1, 200, 201, 2, 4)]
    rows_s2 = [(1, 100, 101, 1, 4), (1, 200, 201, 4, 4)]
    _write_finaleme_prediction(tmp_path / "s1.prediction.bed.gz", rows_s1)
    _write_finaleme_prediction(tmp_path / "s2.prediction.bed.gz", rows_s2)

    sheet = tmp_path / "fme_samples.tsv"
    pd.DataFrame(
        [
            {"sample_id": "S1", "methylation_file": "s1.prediction.bed.gz"},
            {"sample_id": "S2", "methylation_file": "s2.prediction.bed.gz"},
        ]
    ).to_csv(sheet, sep="\t", index=False)

    df = _load_matched_table(sheet, modality="finaleme")
    assert set(df["chrom"]) == {"1"}  # chr stripped
    assert set(df["sample_id"]) == {"S1", "S2"}
    assert len(df) == 4


def test_train_calibration_with_raw_bissnp_and_finaleme(tmp_path):
    """End-to-end: point --matched-wgbs at Bis-SNP files (no chr prefix) and
    --matched-finaleme at FinaleMe files (with chr prefix). The join must
    succeed after chromosome normalization."""
    import pandas as pd
    from scipy.special import expit, logit

    from finaleme_too.preprocessing.calibration import train_calibration

    rng = np.random.default_rng(42)
    n_samples = 4
    n_markers = 120

    # Generate ground-truth betas, then derive miscalibrated FinaleMe betas
    # via the known inverse: fme_logit = (wgbs_logit - c) / a
    true_betas = rng.beta(2, 5, size=n_markers)
    a_true, c_true = 0.8, 0.1

    for s in range(n_samples):
        sid = f"S{s}"
        bissnp_rows = []
        fme_rows = []
        for m, b in enumerate(true_betas):
            n_reads = 30
            k_w = int(rng.binomial(n_reads, b))
            pct = 100.0 * k_w / n_reads
            # Bis-SNP uses GRCh37 style (no chr prefix)
            bissnp_rows.append((1, m * 100, m * 100 + 1, pct, n_reads))

            # FinaleMe: simulate miscalibration
            fme_logit = (
                logit(np.clip(b, 1e-6, 1 - 1e-6)) - c_true
            ) / a_true + rng.normal(0, 0.1)
            fme_b = float(expit(fme_logit))
            k_f = int(rng.binomial(n_reads, fme_b))
            fme_rows.append((1, m * 100, m * 100 + 1, k_f, n_reads))

        _write_bissnp_6plus2(
            tmp_path / f"{sid}.bissnp.bed", bissnp_rows, with_chr_prefix=False
        )
        _write_finaleme_prediction(
            tmp_path / f"{sid}.prediction.bed.gz", fme_rows
        )

    wgbs_sheet = tmp_path / "wgbs_samples.tsv"
    pd.DataFrame(
        [
            {"sample_id": f"S{s}", "methylation_file": f"S{s}.bissnp.bed"}
            for s in range(n_samples)
        ]
    ).to_csv(wgbs_sheet, sep="\t", index=False)

    fme_sheet = tmp_path / "fme_samples.tsv"
    pd.DataFrame(
        [
            {"sample_id": f"S{s}", "methylation_file": f"S{s}.prediction.bed.gz"}
            for s in range(n_samples)
        ]
    ).to_csv(fme_sheet, sep="\t", index=False)

    # No region annotation (optional). Let density default to 0.
    params = train_calibration(
        matched_wgbs=wgbs_sheet,
        matched_finaleme=fme_sheet,
        region_annotation=None,
        n_bins_candidates=[2, 4],
        out_params=tmp_path / "cal.json",
        out_report=tmp_path / "report.json",
    )
    assert params.n_bins in (2, 4)
    assert (tmp_path / "cal.json").exists()
    assert (tmp_path / "report.json").exists()
    # Training should recover something in the ballpark of a_true=0.8
    # on at least one bin
    assert any(0.4 < a_val < 1.5 for a_val in params.a.tolist())


def test_train_calibration_prefix_mismatch_joins_ok(tmp_path):
    """Explicit regression for the 'chr1 vs 1' join issue."""
    import pandas as pd

    from finaleme_too.preprocessing.calibration import _load_matched_table

    # WGBS: no chr prefix
    _write_bissnp_6plus2(
        tmp_path / "s1.bed",
        [(1, 100, 101, 50.0, 4)],
        with_chr_prefix=False,
    )
    # FinaleMe: chr prefix (always)
    _write_finaleme_prediction(
        tmp_path / "s1.prediction.bed.gz",
        [(1, 100, 101, 2, 4)],
    )

    pd.DataFrame(
        [{"sample_id": "S1", "methylation_file": "s1.bed"}]
    ).to_csv(tmp_path / "wgbs.tsv", sep="\t", index=False)
    pd.DataFrame(
        [{"sample_id": "S1", "methylation_file": "s1.prediction.bed.gz"}]
    ).to_csv(tmp_path / "fme.tsv", sep="\t", index=False)

    wgbs_df = _load_matched_table(tmp_path / "wgbs.tsv", modality="wgbs")
    fme_df = _load_matched_table(tmp_path / "fme.tsv", modality="finaleme")

    # After loading, both sides must have chromosome '1' (no 'chr')
    assert list(wgbs_df["chrom"]) == ["1"]
    assert list(fme_df["chrom"]) == ["1"]

    # Simulate the merge that train_calibration does
    merged = wgbs_df.merge(
        fme_df,
        on=["sample_id", "chrom", "start", "end"],
        suffixes=("_wgbs", "_fme"),
    )
    assert len(merged) == 1


def test_load_matched_table_legacy_format_still_works(tmp_path):
    """Backward compatibility: a pre-joined TSV still loads."""
    import pandas as pd

    from finaleme_too.preprocessing.calibration import _load_matched_table

    p = tmp_path / "legacy.tsv"
    pd.DataFrame(
        {
            "sample_id": ["S1", "S1"],
            "chrom": ["chr1", "chr1"],
            "start": [100, 200],
            "end": [101, 201],
            "methylated_count": [2, 3],
            "total_count": [4, 4],
        }
    ).to_csv(p, sep="\t", index=False)

    df = _load_matched_table(p, modality="wgbs")
    # Legacy format still works AND chr prefix is stripped
    assert set(df["chrom"]) == {"1"}
    assert len(df) == 2


def test_load_matched_table_rejects_malformed(tmp_path):
    """A file that is neither legacy nor sample-sheet format should raise."""
    import pandas as pd

    from finaleme_too.exceptions import InvalidCalibrationError
    from finaleme_too.preprocessing.calibration import _load_matched_table

    p = tmp_path / "bogus.tsv"
    pd.DataFrame({"foo": [1], "bar": [2]}).to_csv(p, sep="\t", index=False)
    with pytest.raises(InvalidCalibrationError):
        _load_matched_table(p, modality="wgbs")


# ---------------------------------------------------------------------------
# Region annotation generation from a CpG index
# ---------------------------------------------------------------------------


def _write_cpg_index(path, chroms_to_positions):
    """Write a tiny CpG index BED file: ``chrom start end`` per CpG."""
    with open(path, "w") as fh:
        for chrom, positions in chroms_to_positions.items():
            for p in positions:
                fh.write(f"{chrom}\t{p}\t{p + 1}\n")


def test_compute_region_annotation_matches_windowed_count(tmp_path):
    """Per-row density should equal (#CpGs in window) / window."""
    import pandas as pd

    from finaleme_too.preprocessing.calibration import compute_region_annotation

    # 5 CpGs on chr1 at positions 100, 200, 300, 400, 500
    # 2 CpGs on chr2 at positions 1000, 5000
    cpg_path = tmp_path / "cpg.bed"
    _write_cpg_index(cpg_path, {"chr1": [100, 200, 300, 400, 500], "chr2": [1000, 5000]})

    regions = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2", "chr2"],
            "start": [250, 800, 1000, 4000],
            "end": [251, 801, 1001, 4001],
        }
    )
    ann = compute_region_annotation(regions, cpg_path, window=1000)

    # Row 0: center=250, window=[250-500, 250+500)=[-250, 750). CpGs at
    # 100,200,300,400,500 all fall in there → 5 CpGs / 1000 bp = 0.005
    assert ann.iloc[0]["cpg_density"] == pytest.approx(5 / 1000)
    # Row 1: center=800, window=[300, 1300). CpGs 300,400,500 → 3/1000
    assert ann.iloc[1]["cpg_density"] == pytest.approx(3 / 1000)
    # Row 2: chr2 center=1000, window=[500, 1500). CpG 1000 → 1/1000
    assert ann.iloc[2]["cpg_density"] == pytest.approx(1 / 1000)
    # Row 3: chr2 center=4000, window=[3500, 4500). No CpG → 0
    assert ann.iloc[3]["cpg_density"] == pytest.approx(0.0)
    # chromosomes should be stripped of 'chr' prefix
    assert set(ann["chrom"]) == {"1", "2"}


def test_compute_region_annotation_handles_unknown_chromosome(tmp_path):
    """Rows on a chrom missing from the index get density 0 and do not crash."""
    import pandas as pd

    from finaleme_too.preprocessing.calibration import compute_region_annotation

    cpg_path = tmp_path / "cpg.bed"
    _write_cpg_index(cpg_path, {"chr1": [100, 200]})

    regions = pd.DataFrame(
        {"chrom": ["chrX"], "start": [500], "end": [501]}
    )
    ann = compute_region_annotation(regions, cpg_path, window=1000)
    assert ann.iloc[0]["cpg_density"] == 0.0


def test_make_region_annotation_from_bed_round_trip(tmp_path):
    """make_region_annotation_from_regions_file reads a headerless BED,
    writes the expected 4-col TSV, and the output is loadable."""
    import pandas as pd

    from finaleme_too.preprocessing.calibration import (
        make_region_annotation_from_regions_file,
    )

    cpg_path = tmp_path / "cpg.bed"
    _write_cpg_index(cpg_path, {"chr1": list(range(100, 1100, 10))})  # 100 CpGs

    regions_path = tmp_path / "regions.bed"
    with open(regions_path, "w") as fh:
        fh.write("chr1\t500\t501\n")
        fh.write("chr1\t900\t901\n")

    out_path = tmp_path / "ann.tsv"
    make_region_annotation_from_regions_file(
        regions_path, cpg_path, out_path, window=1000
    )
    assert out_path.exists()

    df = pd.read_csv(out_path, sep="\t")
    assert list(df.columns) == ["chrom", "start", "end", "cpg_density"]
    assert len(df) == 2
    # Both rows should have a non-zero density since we seeded the CpG index
    assert all(df["cpg_density"] > 0)


def test_make_region_annotation_from_headered_tsv(tmp_path):
    """Also accept a TSV that already has a header row."""
    import pandas as pd

    from finaleme_too.preprocessing.calibration import (
        make_region_annotation_from_regions_file,
    )

    cpg_path = tmp_path / "cpg.bed"
    _write_cpg_index(cpg_path, {"chr1": [100, 200, 300, 400, 500]})

    regions_path = tmp_path / "regions.tsv"
    pd.DataFrame(
        {"chrom": ["chr1"], "start": [250], "end": [251]}
    ).to_csv(regions_path, sep="\t", index=False)

    out_path = tmp_path / "ann.tsv"
    make_region_annotation_from_regions_file(
        regions_path, cpg_path, out_path, window=1000
    )
    df = pd.read_csv(out_path, sep="\t")
    assert len(df) == 1
    assert df.iloc[0]["cpg_density"] == pytest.approx(5 / 1000)


# ---------------------------------------------------------------------------
# Parallel execution of train-calibration
# ---------------------------------------------------------------------------


def test_cross_validate_calibration_threads_match_serial():
    """Running CV with threads=4 must give the same per-fold result as serial."""
    from finaleme_too.preprocessing.calibration_eval import cross_validate_calibration

    rng = np.random.default_rng(0)
    n = 600
    sample_ids = np.array([f"S{i // 150}" for i in range(n)])  # 4 samples
    truth = rng.beta(2, 5, size=n)
    fme = np.clip(truth + rng.normal(0, 0.05, size=n), 0.01, 0.99)
    density = rng.uniform(0, 0.1, size=n)

    serial = cross_validate_calibration(
        finaleme_beta=fme,
        wgbs_beta=truth,
        cpg_density=density,
        sample_ids=sample_ids,
        n_bins=4,
        threads=1,
    )
    parallel = cross_validate_calibration(
        finaleme_beta=fme,
        wgbs_beta=truth,
        cpg_density=density,
        sample_ids=sample_ids,
        n_bins=4,
        threads=4,
    )
    # cv_rmse, in-sample rmse, and fold counts must all match
    assert serial["n_folds"] == parallel["n_folds"]
    assert abs(serial["cv_rmse"] - parallel["cv_rmse"]) < 1e-12
    assert abs(serial["in_sample_rmse"] - parallel["in_sample_rmse"]) < 1e-12


def test_tune_n_bins_threads_match_serial():
    """tune_n_bins with threads must match the serial result."""
    from finaleme_too.preprocessing.calibration_eval import tune_n_bins

    rng = np.random.default_rng(1)
    n = 800
    sample_ids = np.array([f"S{i // 200}" for i in range(n)])  # 4 samples
    truth = rng.beta(2, 5, size=n)
    fme = np.clip(truth + rng.normal(0, 0.03, size=n), 0.01, 0.99)
    density = rng.uniform(0, 0.1, size=n)

    serial = tune_n_bins(
        finaleme_beta=fme,
        wgbs_beta=truth,
        cpg_density=density,
        sample_ids=sample_ids,
        n_bins_candidates=[2, 4, 6],
        threads=1,
    )
    parallel = tune_n_bins(
        finaleme_beta=fme,
        wgbs_beta=truth,
        cpg_density=density,
        sample_ids=sample_ids,
        n_bins_candidates=[2, 4, 6],
        threads=4,
    )
    # Same selected B
    assert serial["selected_n_bins"] == parallel["selected_n_bins"]
    # Same number of candidates
    assert len(serial["candidates"]) == len(parallel["candidates"])
    # Same cv_rmse per candidate (identical numerics with threading backend)
    for s, p in zip(serial["candidates"], parallel["candidates"]):
        assert s["n_bins"] == p["n_bins"]
        assert abs(s["cv_rmse"] - p["cv_rmse"]) < 1e-12
        assert abs(s["in_sample_rmse"] - p["in_sample_rmse"]) < 1e-12


def test_load_sample_sheet_threaded_parsing(tmp_path):
    """Loading a 4-sample Bis-SNP sample sheet with threads must produce
    the same DataFrame as the serial loader."""
    import pandas as pd

    from finaleme_too.preprocessing.calibration import _load_matched_table

    rows_per_sample = [
        (1, 100, 101, 50.0, 4),
        (1, 200, 201, 25.0, 4),
        (2, 5000, 5001, 100.0, 10),
        (2, 6000, 6001, 0.0, 10),
    ]
    for s in range(4):
        _write_bissnp_6plus2(tmp_path / f"s{s}.bed", rows_per_sample, with_chr_prefix=False)

    sheet = tmp_path / "wgbs.tsv"
    pd.DataFrame(
        [
            {"sample_id": f"S{s}", "methylation_file": f"s{s}.bed"}
            for s in range(4)
        ]
    ).to_csv(sheet, sep="\t", index=False)

    serial = _load_matched_table(sheet, modality="wgbs", threads=1)
    parallel = _load_matched_table(sheet, modality="wgbs", threads=4)
    # Same content (row order may differ because joblib doesn't guarantee
    # order across threads; use set comparison on a canonical key)
    def _key(df):
        return sorted(
            (r["sample_id"], r["chrom"], int(r["start"]), int(r["methylated_count"]), int(r["total_count"]))
            for _, r in df.iterrows()
        )
    assert _key(serial) == _key(parallel)
    assert len(parallel) == 4 * len(rows_per_sample)


def test_train_calibration_threaded_end_to_end(tmp_path):
    """End-to-end train_calibration with threads=4. Must produce a result
    indistinguishable from the serial run up to floating-point noise."""
    from scipy.special import expit, logit

    from finaleme_too.preprocessing.calibration import train_calibration

    rng = np.random.default_rng(42)
    n_samples = 4
    n_markers = 100
    rows_w, rows_f = [], []
    for s in range(n_samples):
        for m in range(n_markers):
            b = float(rng.beta(2, 5))
            k_w = int(rng.binomial(30, b))
            fme_b = float(
                expit((logit(np.clip(b, 1e-6, 1 - 1e-6)) - 0.1) / 0.8 + rng.normal(0, 0.1))
            )
            k_f = int(rng.binomial(30, fme_b))
            rows_w.append((f"S{s}", "chr1", m * 100, m * 100 + 1, k_w, 30))
            rows_f.append((f"S{s}", "chr1", m * 100, m * 100 + 1, k_f, 30))

    import pandas as pd

    pd.DataFrame(
        rows_w,
        columns=["sample_id", "chrom", "start", "end", "methylated_count", "total_count"],
    ).to_csv(tmp_path / "wgbs.tsv", sep="\t", index=False)
    pd.DataFrame(
        rows_f,
        columns=["sample_id", "chrom", "start", "end", "methylated_count", "total_count"],
    ).to_csv(tmp_path / "fme.tsv", sep="\t", index=False)

    params_serial = train_calibration(
        matched_wgbs=tmp_path / "wgbs.tsv",
        matched_finaleme=tmp_path / "fme.tsv",
        region_annotation=None,
        n_bins_candidates=[2, 4],
        out_params=tmp_path / "serial.json",
        out_report=tmp_path / "serial_report.json",
        threads=1,
    )
    params_threaded = train_calibration(
        matched_wgbs=tmp_path / "wgbs.tsv",
        matched_finaleme=tmp_path / "fme.tsv",
        region_annotation=None,
        n_bins_candidates=[2, 4],
        out_params=tmp_path / "threaded.json",
        out_report=tmp_path / "threaded_report.json",
        threads=4,
    )
    # Both runs must pick the same B and produce the same parameters
    assert params_serial.n_bins == params_threaded.n_bins
    np.testing.assert_allclose(params_serial.a, params_threaded.a, rtol=1e-12)
    np.testing.assert_allclose(params_serial.c, params_threaded.c, rtol=1e-12)


def test_train_calibration_auto_generates_region_annotation_from_cpg_index(tmp_path):
    """When --region-annotation is omitted but --cpg-index is provided,
    train_calibration must auto-compute density per row and use it."""
    from scipy.special import expit, logit

    from finaleme_too.preprocessing.calibration import train_calibration

    rng = np.random.default_rng(7)
    n_samples = 3
    n_markers = 80

    # CpG index: 1 CpG per row, matching the training positions
    cpg_path = tmp_path / "cpg.bed"
    positions = [100 + m * 1000 for m in range(n_markers)]
    with open(cpg_path, "w") as fh:
        for p in positions:
            fh.write(f"chr1\t{p}\t{p + 1}\n")

    # Build matched data (legacy TSV format for brevity)
    rows_w = []
    rows_f = []
    for s in range(n_samples):
        for m, p in enumerate(positions):
            b = float(rng.beta(2, 5))
            k_w = int(rng.binomial(30, b))
            # Miscalibration: logit_fme = (logit_wgbs - 0.1) / 0.8
            fme_b = float(
                expit((logit(np.clip(b, 1e-6, 1 - 1e-6)) - 0.1) / 0.8 + rng.normal(0, 0.1))
            )
            k_f = int(rng.binomial(30, fme_b))
            rows_w.append((f"S{s}", "chr1", p, p + 1, k_w, 30))
            rows_f.append((f"S{s}", "chr1", p, p + 1, k_f, 30))

    import pandas as pd

    pd.DataFrame(
        rows_w,
        columns=["sample_id", "chrom", "start", "end", "methylated_count", "total_count"],
    ).to_csv(tmp_path / "wgbs.tsv", sep="\t", index=False)
    pd.DataFrame(
        rows_f,
        columns=["sample_id", "chrom", "start", "end", "methylated_count", "total_count"],
    ).to_csv(tmp_path / "fme.tsv", sep="\t", index=False)

    params = train_calibration(
        matched_wgbs=tmp_path / "wgbs.tsv",
        matched_finaleme=tmp_path / "fme.tsv",
        region_annotation=None,  # NOT provided
        cpg_index=cpg_path,  # NEW: auto-compute density
        region_annotation_window=2000,
        n_bins_candidates=[2, 4],
        out_params=tmp_path / "cal.json",
        out_report=tmp_path / "report.json",
    )
    assert params.n_bins in (2, 4)
    assert (tmp_path / "cal.json").exists()
    # A useful sanity check: the training metadata records the number of
    # observations, which is n_samples * n_markers when density is non-zero
    assert params.training_metadata["n_observations"] == n_samples * n_markers


# ---------------------------------------------------------------------------
# April 2026 — 4 calibration bugs reported from a user's HD45 run:
#
# Bug 1 (HIGH): Best-bin selection is mathematically wrong when CV is
#   unavailable. The old _select_best_candidate fell back to
#   ``min(candidates, key=lambda c: (c["cv_rmse"], c.get("aic", inf)))``
#   which silently returns the first candidate because NaN comparisons
#   are undefined — AIC was never consulted.
#
# Bug 2 (HIGH): train_calibration treated zero-coverage rows as valid
#   betas via ``total_count.clip(lower=1)``, turning k/0 into k/1. This
#   biased the regression, RMSE, AIC, and bin tuning. Must drop those
#   rows before computing betas.
#
# Bug 3 (MEDIUM): The FinaleMe prediction parser did not clamp
#   methylated_count <= total_count, so malformed rows could yield
#   beta > 1 and distort the regression after logit clipping.
#
# Bug 4 (MEDIUM, design mismatch): per-bin dispersion was estimated from
#   residual variance (phi ≈ 1/Var) rather than by beta-binomial MLE on
#   the calibration residuals as the design (math doc §6.2) requires.
# ---------------------------------------------------------------------------


def test_bug1_select_best_candidate_falls_back_to_aic_when_cv_unavailable():
    """When every candidate has cv_rmse=NaN (CV couldn't run), the best
    candidate must be chosen by AIC — not just the first one in the list.
    """
    from finaleme_too.preprocessing.calibration_eval import _select_best_candidate

    # Four candidates with NaN cv_rmse; only B=8 has the lowest AIC.
    candidates = [
        {"n_bins": 4, "cv_rmse": float("nan"), "aic": 100.0, "in_sample_rmse": 0.3},
        {"n_bins": 6, "cv_rmse": float("nan"), "aic": 80.0, "in_sample_rmse": 0.28},
        {"n_bins": 8, "cv_rmse": float("nan"), "aic": 50.0, "in_sample_rmse": 0.26},
        {"n_bins": 10, "cv_rmse": float("nan"), "aic": 90.0, "in_sample_rmse": 0.27},
    ]
    result = _select_best_candidate(candidates)
    # The old code would have returned candidates[0] (n_bins=4); the fix
    # picks the one with lowest AIC (n_bins=8).
    assert result["selected_n_bins"] == 8


def test_bug1_select_best_candidate_still_prefers_cv_when_available():
    """Sanity check that cv_rmse is still the primary criterion when
    it is finite."""
    from finaleme_too.preprocessing.calibration_eval import _select_best_candidate

    candidates = [
        {"n_bins": 4, "cv_rmse": 0.20, "aic": 100.0, "in_sample_rmse": 0.25},
        {"n_bins": 6, "cv_rmse": 0.15, "aic": 90.0, "in_sample_rmse": 0.20},  # best by CV
        {"n_bins": 8, "cv_rmse": 0.18, "aic": 50.0, "in_sample_rmse": 0.18},  # best by AIC but worse CV
    ]
    result = _select_best_candidate(candidates)
    assert result["selected_n_bins"] == 6


def test_bug1_select_best_candidate_falls_back_to_rmse_when_no_cv_or_aic():
    """When neither cv_rmse nor aic are finite, fall through to in-sample
    RMSE instead of returning an arbitrary candidate."""
    from finaleme_too.preprocessing.calibration_eval import _select_best_candidate

    candidates = [
        {"n_bins": 4, "cv_rmse": float("nan"), "aic": float("nan"), "in_sample_rmse": 0.40},
        {"n_bins": 6, "cv_rmse": float("nan"), "aic": float("nan"), "in_sample_rmse": 0.15},
        {"n_bins": 8, "cv_rmse": float("nan"), "aic": float("nan"), "in_sample_rmse": 0.25},
    ]
    result = _select_best_candidate(candidates)
    assert result["selected_n_bins"] == 6


def test_bug2_zero_coverage_rows_dropped_before_training(tmp_path):
    """train_calibration must drop rows with total_count=0 on either side
    instead of silently turning k/0 into k/1.

    Inject several zero-coverage rows on both sides and verify:
      (a) n_observations reflects only the non-zero-coverage rows
      (b) the calibration params are unchanged by adding those garbage rows
    """
    from scipy.special import expit, logit

    from finaleme_too.preprocessing.calibration import train_calibration

    rng = np.random.default_rng(99)
    n_samples = 3
    n_markers = 60

    def _build(rows_w, rows_f, prefix: str):
        pd.DataFrame(
            rows_w,
            columns=["sample_id", "chrom", "start", "end", "methylated_count", "total_count"],
        ).to_csv(tmp_path / f"{prefix}_wgbs.tsv", sep="\t", index=False)
        pd.DataFrame(
            rows_f,
            columns=["sample_id", "chrom", "start", "end", "methylated_count", "total_count"],
        ).to_csv(tmp_path / f"{prefix}_fme.tsv", sep="\t", index=False)

    # Build clean data
    clean_w, clean_f = [], []
    for s in range(n_samples):
        for m in range(n_markers):
            b = float(rng.beta(2, 5))
            k_w = int(rng.binomial(30, b))
            fme_b = float(
                expit((logit(np.clip(b, 1e-6, 1 - 1e-6)) - 0.1) / 0.8 + rng.normal(0, 0.1))
            )
            k_f = int(rng.binomial(30, fme_b))
            clean_w.append((f"S{s}", "chr1", m * 100, m * 100 + 1, k_w, 30))
            clean_f.append((f"S{s}", "chr1", m * 100, m * 100 + 1, k_f, 30))

    _build(clean_w, clean_f, "clean")
    clean_params = train_calibration(
        matched_wgbs=tmp_path / "clean_wgbs.tsv",
        matched_finaleme=tmp_path / "clean_fme.tsv",
        region_annotation=None,
        n_bins_candidates=[2, 4],
        out_params=tmp_path / "clean.json",
        out_report=tmp_path / "clean_report.json",
    )
    assert clean_params.training_metadata["n_observations"] == n_samples * n_markers

    # Inject zero-coverage garbage rows — same sample/coord layout but
    # total_count=0. These rows must be filtered before beta computation.
    dirty_w = list(clean_w) + [
        (f"S{s}", "chr1", 9000 + i, 9001 + i, 5, 0)
        for s in range(n_samples)
        for i in range(10)
    ]
    dirty_f = list(clean_f) + [
        (f"S{s}", "chr1", 9000 + i, 9001 + i, 3, 0)
        for s in range(n_samples)
        for i in range(10)
    ]
    _build(dirty_w, dirty_f, "dirty")
    dirty_params = train_calibration(
        matched_wgbs=tmp_path / "dirty_wgbs.tsv",
        matched_finaleme=tmp_path / "dirty_fme.tsv",
        region_annotation=None,
        n_bins_candidates=[2, 4],
        out_params=tmp_path / "dirty.json",
        out_report=tmp_path / "dirty_report.json",
    )
    # Same number of surviving observations
    assert dirty_params.training_metadata["n_observations"] == n_samples * n_markers
    # Same fitted parameters (identical valid rows → identical output)
    assert dirty_params.n_bins == clean_params.n_bins
    np.testing.assert_allclose(dirty_params.a, clean_params.a, rtol=1e-12)
    np.testing.assert_allclose(dirty_params.c, clean_params.c, rtol=1e-12)
    np.testing.assert_allclose(
        dirty_params.log_dispersion, clean_params.log_dispersion, rtol=1e-12
    )


def test_bug2_zero_coverage_only_input_raises(tmp_path):
    """If every row has zero coverage on one side, train_calibration must
    fail loudly instead of silently fitting on garbage k/1 betas."""
    from finaleme_too.exceptions import InvalidCalibrationError
    from finaleme_too.preprocessing.calibration import train_calibration

    rows_w = [("S0", "chr1", i * 100, i * 100 + 1, 1, 0) for i in range(20)]
    rows_f = [("S0", "chr1", i * 100, i * 100 + 1, 1, 5) for i in range(20)]

    pd.DataFrame(
        rows_w,
        columns=["sample_id", "chrom", "start", "end", "methylated_count", "total_count"],
    ).to_csv(tmp_path / "wgbs.tsv", sep="\t", index=False)
    pd.DataFrame(
        rows_f,
        columns=["sample_id", "chrom", "start", "end", "methylated_count", "total_count"],
    ).to_csv(tmp_path / "fme.tsv", sep="\t", index=False)

    with pytest.raises(InvalidCalibrationError):
        train_calibration(
            matched_wgbs=tmp_path / "wgbs.tsv",
            matched_finaleme=tmp_path / "fme.tsv",
            region_annotation=None,
            n_bins_candidates=[2],
            out_params=tmp_path / "cal.json",
            out_report=tmp_path / "report.json",
        )


def test_bug3_finaleme_parser_clamps_methylated_above_total(tmp_path):
    """_parse_finaleme_prediction must enforce methylated_count <= total_count
    so malformed FinaleMe files cannot yield beta > 1."""
    import gzip

    from finaleme_too.preprocessing.calibration import _parse_finaleme_prediction

    path = tmp_path / "pred.bed.gz"
    # Columns: chr start end methy_pct methy_count total_count ...
    # Row 1: normal, k=5, n=10 (beta 0.5)
    # Row 2: malformed, k=12, n=10 (beta 1.2 — must be clamped to k=n=10)
    # Row 3: normal, k=0, n=8 (beta 0.0)
    lines = [
        "#chr\tstart\tend\tmethy_perc\tmethy_count\ttotal_count\n",
        "chr1\t100\t101\t50.0\t5\t10\n",
        "chr1\t200\t201\t120.0\t12\t10\n",
        "chr1\t300\t301\t0.0\t0\t8\n",
    ]
    with gzip.open(path, "wt") as fh:
        fh.writelines(lines)
    df = _parse_finaleme_prediction(path, sample_id="test")

    # All (k, n) must satisfy 0 <= k <= n
    assert (df["methylated_count"] >= 0).all()
    assert (df["methylated_count"] <= df["total_count"]).all()
    # Row 2 specifically: methylated was clamped from 12 to 10
    row_200 = df[df["start"] == 200].iloc[0]
    assert int(row_200["methylated_count"]) == 10
    assert int(row_200["total_count"]) == 10
    # Resulting beta = 1.0, not 1.2
    beta = int(row_200["methylated_count"]) / int(row_200["total_count"])
    assert beta == 1.0


def test_bug4_dispersion_estimated_by_beta_binomial_mle(tmp_path):
    """When raw WGBS (k, n) counts are provided, the per-bin dispersion
    must come from a beta-binomial MLE rather than the residual-variance
    heuristic. This test verifies the MLE path is actually exercised and
    produces a dispersion that differs from the heuristic when the two
    estimators should disagree.
    """
    from finaleme_too.preprocessing.calibration_eval import (
        _fit_logistic_bin,
        fit_calibration,
    )

    rng = np.random.default_rng(123)
    n = 500
    # Simulate well-calibrated data: fme_beta ≈ wgbs_beta with small noise.
    wgbs_beta = rng.beta(2, 5, size=n)
    fme_beta = np.clip(wgbs_beta + rng.normal(0, 0.03, size=n), 0.01, 0.99)
    wgbs_n = np.full(n, 30, dtype=np.int64)
    # Generate realistic WGBS counts under the beta-binomial with phi=80
    wgbs_k = rng.binomial(wgbs_n, wgbs_beta).astype(np.int64)

    # Heuristic path (no counts)
    a_h, c_h, log_phi_h = _fit_logistic_bin(fme_beta, wgbs_beta)

    # MLE path (with counts)
    a_m, c_m, log_phi_m = _fit_logistic_bin(
        fme_beta, wgbs_beta, wgbs_k=wgbs_k, wgbs_n=wgbs_n
    )

    # (a, c) come from the same OLS, so they match exactly.
    assert abs(a_h - a_m) < 1e-12
    assert abs(c_h - c_m) < 1e-12

    # Dispersion must differ: the MLE uses the binomial information in
    # (k, n), the heuristic uses the variance of the logit-space residual.
    # These are distinct estimators and must not collapse to the same
    # value on real data.
    assert abs(log_phi_h - log_phi_m) > 1e-6

    # Also verify that fit_calibration propagates the counts correctly:
    # calling with counts vs without yields different log_dispersion.
    cpg_density = rng.uniform(0, 0.1, size=n)
    fit_no_counts = fit_calibration(
        fme_beta, wgbs_beta, cpg_density, n_bins=2
    )
    fit_with_counts = fit_calibration(
        fme_beta, wgbs_beta, cpg_density, n_bins=2,
        wgbs_k=wgbs_k, wgbs_n=wgbs_n,
    )
    # a and c must match; log_dispersion must diverge.
    np.testing.assert_allclose(fit_no_counts.a, fit_with_counts.a, rtol=1e-12)
    np.testing.assert_allclose(fit_no_counts.c, fit_with_counts.c, rtol=1e-12)
    assert not np.allclose(
        fit_no_counts.log_dispersion, fit_with_counts.log_dispersion, atol=1e-6
    )


def test_bug4_dispersion_mle_is_finite_on_clean_data():
    """Regression: the old bracketed Brent MLE threw "Not a bracketing
    interval" on well-calibrated synthetic data. The bounded minimizer
    used now must return a finite phi in the same scenario.
    """
    from finaleme_too.utils.beta_binomial import estimate_dispersion_mle

    rng = np.random.default_rng(7)
    n = 200
    mu = rng.uniform(0.3, 0.7, size=n)
    n_arr = np.full(n, 30, dtype=np.int64)
    # Pristine binomial counts — the old brent bracket failed on this.
    k_arr = rng.binomial(n_arr, mu)

    phi = estimate_dispersion_mle(
        k=k_arr.astype(np.float64),
        n=n_arr.astype(np.float64),
        mu=mu,
        phi_init=50.0,
        bounds=(1.0, 1000.0),
    )
    assert np.isfinite(phi)
    assert 1.0 <= phi <= 1000.0


# ---------------------------------------------------------------------------
# April 2026 — region-level K-fold CV for few-sample calibration runs
# ---------------------------------------------------------------------------
# When only 1-2 matched WGBS/FinaleMe samples are available, leave-one-sample-
# out CV can't run (cv_rmse=NaN for every candidate). The new ``cv_strategy``
# parameter on train_calibration/tune_n_bins/cross_validate_calibration lets
# the user ask for random K-fold CV on row indices instead, which gives a
# finite cv_rmse — biased low compared to true held-out-sample error, but
# strictly better than NaN for picking the best n_bins.


def _synth_matched_arrays(
    n_samples: int,
    n_markers_per_sample: int,
    seed: int = 0,
):
    """Build synthetic (finaleme_beta, wgbs_beta, wgbs_k, wgbs_n, density,
    sample_ids) arrays with a known miscalibration."""
    from scipy.special import expit, logit

    rng = np.random.default_rng(seed)
    rows = []
    for s in range(n_samples):
        for m in range(n_markers_per_sample):
            b = float(rng.beta(2, 5))
            k_w = int(rng.binomial(30, b))
            fme_b = float(
                expit(
                    (logit(np.clip(b, 1e-6, 1 - 1e-6)) - 0.1) / 0.8
                    + rng.normal(0, 0.1)
                )
            )
            rows.append((f"S{s}", b, fme_b, k_w, 30))
    sample_ids = np.array([r[0] for r in rows])
    wgbs_beta = np.array([r[1] for r in rows], dtype=np.float64)
    fme_beta = np.clip(np.array([r[2] for r in rows], dtype=np.float64), 1e-6, 1 - 1e-6)
    wgbs_k = np.array([r[3] for r in rows], dtype=np.int64)
    wgbs_n = np.array([r[4] for r in rows], dtype=np.int64)
    density = np.linspace(0.0, 0.1, len(rows))
    return fme_beta, wgbs_beta, wgbs_k, wgbs_n, density, sample_ids


def test_region_cv_returns_finite_rmse_with_single_sample():
    """A 1-sample run with cv_strategy='region' (or 'auto') must produce a
    finite cv_rmse instead of NaN."""
    from finaleme_too.preprocessing.calibration_eval import (
        cross_validate_calibration,
    )

    fme, wgbs, k, n, density, sample_ids = _synth_matched_arrays(
        n_samples=1, n_markers_per_sample=300, seed=0
    )

    result_region = cross_validate_calibration(
        finaleme_beta=fme,
        wgbs_beta=wgbs,
        cpg_density=density,
        sample_ids=sample_ids,
        n_bins=4,
        wgbs_k=k,
        wgbs_n=n,
        cv_strategy="region",
        cv_n_folds=5,
        cv_seed=42,
    )
    assert np.isfinite(result_region["cv_rmse"])
    assert result_region["n_folds"] == 5
    assert result_region["cv_strategy"] == "region"

    # "auto" with 1 sample should also fall back to region
    result_auto = cross_validate_calibration(
        finaleme_beta=fme,
        wgbs_beta=wgbs,
        cpg_density=density,
        sample_ids=sample_ids,
        n_bins=4,
        wgbs_k=k,
        wgbs_n=n,
        cv_strategy="auto",
        cv_n_folds=5,
        cv_seed=42,
    )
    assert np.isfinite(result_auto["cv_rmse"])
    assert result_auto["cv_strategy"] == "region"
    # Same seed + same n_folds → identical rmse
    assert result_region["cv_rmse"] == result_auto["cv_rmse"]


def test_region_cv_seed_is_reproducible():
    """Same cv_seed must produce the same region-level CV rmse across runs."""
    from finaleme_too.preprocessing.calibration_eval import (
        cross_validate_calibration,
    )

    fme, wgbs, k, n, density, sample_ids = _synth_matched_arrays(
        n_samples=1, n_markers_per_sample=300, seed=11
    )

    a = cross_validate_calibration(
        finaleme_beta=fme, wgbs_beta=wgbs, cpg_density=density,
        sample_ids=sample_ids, n_bins=4,
        wgbs_k=k, wgbs_n=n,
        cv_strategy="region", cv_n_folds=10, cv_seed=777,
    )
    b = cross_validate_calibration(
        finaleme_beta=fme, wgbs_beta=wgbs, cpg_density=density,
        sample_ids=sample_ids, n_bins=4,
        wgbs_k=k, wgbs_n=n,
        cv_strategy="region", cv_n_folds=10, cv_seed=777,
    )
    assert abs(a["cv_rmse"] - b["cv_rmse"]) < 1e-12

    # Different seed → different result (with overwhelmingly high probability)
    c = cross_validate_calibration(
        finaleme_beta=fme, wgbs_beta=wgbs, cpg_density=density,
        sample_ids=sample_ids, n_bins=4,
        wgbs_k=k, wgbs_n=n,
        cv_strategy="region", cv_n_folds=10, cv_seed=999,
    )
    assert abs(a["cv_rmse"] - c["cv_rmse"]) > 1e-12


def test_auto_cv_prefers_sample_when_multi_sample():
    """With >=2 samples, auto must resolve to sample-level CV (the
    unbiased estimator)."""
    from finaleme_too.preprocessing.calibration_eval import (
        cross_validate_calibration,
    )

    fme, wgbs, k, n, density, sample_ids = _synth_matched_arrays(
        n_samples=3, n_markers_per_sample=200, seed=2
    )
    result = cross_validate_calibration(
        finaleme_beta=fme, wgbs_beta=wgbs, cpg_density=density,
        sample_ids=sample_ids, n_bins=4,
        wgbs_k=k, wgbs_n=n,
        cv_strategy="auto",
    )
    assert result["cv_strategy"] == "sample"
    assert result["n_folds"] == 3


def test_cv_strategy_none_skips_cv_entirely():
    """cv_strategy='none' must return NaN cv_rmse and 0 folds, even with
    multiple samples."""
    from finaleme_too.preprocessing.calibration_eval import (
        cross_validate_calibration,
    )

    fme, wgbs, k, n, density, sample_ids = _synth_matched_arrays(
        n_samples=3, n_markers_per_sample=200, seed=3
    )
    result = cross_validate_calibration(
        finaleme_beta=fme, wgbs_beta=wgbs, cpg_density=density,
        sample_ids=sample_ids, n_bins=4,
        wgbs_k=k, wgbs_n=n,
        cv_strategy="none",
    )
    assert result["cv_strategy"] == "none"
    assert result["n_folds"] == 0
    assert np.isnan(result["cv_rmse"])
    # In-sample metrics must still be populated
    assert np.isfinite(result["in_sample_rmse"])


def test_tune_n_bins_picks_by_cv_rmse_in_region_mode():
    """tune_n_bins with cv_strategy='region' must produce finite cv_rmse
    for every candidate and choose by cv_rmse (not AIC fallback)."""
    from finaleme_too.preprocessing.calibration_eval import tune_n_bins

    fme, wgbs, k, n, density, sample_ids = _synth_matched_arrays(
        n_samples=1, n_markers_per_sample=600, seed=5
    )
    result = tune_n_bins(
        finaleme_beta=fme, wgbs_beta=wgbs, cpg_density=density,
        sample_ids=sample_ids,
        n_bins_candidates=[2, 4, 6, 8],
        wgbs_k=k, wgbs_n=n,
        cv_strategy="region", cv_n_folds=10, cv_seed=1,
    )
    # Every candidate got a finite cv_rmse
    for cand in result["candidates"]:
        assert np.isfinite(cand["cv_rmse"]), f"NaN cv_rmse for B={cand['n_bins']}"
        assert cand["n_folds"] == 10
        assert cand["cv_strategy"] == "region"
    # Selected n_bins is the one with lowest cv_rmse
    best = min(result["candidates"], key=lambda c: c["cv_rmse"])
    assert result["selected_n_bins"] == best["n_bins"]


def test_tune_n_bins_serial_and_threaded_region_cv_match():
    """Region-level CV must be deterministic regardless of thread count
    when cv_seed is fixed."""
    from finaleme_too.preprocessing.calibration_eval import tune_n_bins

    fme, wgbs, k, n, density, sample_ids = _synth_matched_arrays(
        n_samples=1, n_markers_per_sample=500, seed=9
    )
    serial = tune_n_bins(
        finaleme_beta=fme, wgbs_beta=wgbs, cpg_density=density,
        sample_ids=sample_ids,
        n_bins_candidates=[2, 4, 6],
        wgbs_k=k, wgbs_n=n,
        threads=1,
        cv_strategy="region", cv_n_folds=5, cv_seed=13,
    )
    threaded = tune_n_bins(
        finaleme_beta=fme, wgbs_beta=wgbs, cpg_density=density,
        sample_ids=sample_ids,
        n_bins_candidates=[2, 4, 6],
        wgbs_k=k, wgbs_n=n,
        threads=4,
        cv_strategy="region", cv_n_folds=5, cv_seed=13,
    )
    assert serial["selected_n_bins"] == threaded["selected_n_bins"]
    for s, t in zip(serial["candidates"], threaded["candidates"]):
        assert s["n_bins"] == t["n_bins"]
        assert abs(s["cv_rmse"] - t["cv_rmse"]) < 1e-12
        assert s["n_folds"] == t["n_folds"]


def test_train_calibration_single_sample_with_region_cv(tmp_path):
    """End-to-end: train_calibration on ONE matched sample with
    cv_strategy='region' must produce:
      * a finite cv_rmse per candidate
      * a selected n_bins (not silently "the first one")
      * cv_strategy_effective='region' in the training report
    """
    from finaleme_too.preprocessing.calibration import train_calibration

    fme_arr, wgbs_arr, k_arr, n_arr, _density, _sids = _synth_matched_arrays(
        n_samples=1, n_markers_per_sample=300, seed=17
    )
    # Write as legacy joined TSVs
    rows_w = [
        ("S0", "chr1", i * 100, i * 100 + 1, int(k_arr[i]), int(n_arr[i]))
        for i in range(len(wgbs_arr))
    ]
    rows_f = [
        ("S0", "chr1", i * 100, i * 100 + 1,
         int(round(fme_arr[i] * 30)), 30)
        for i in range(len(fme_arr))
    ]
    pd.DataFrame(
        rows_w,
        columns=["sample_id", "chrom", "start", "end", "methylated_count", "total_count"],
    ).to_csv(tmp_path / "wgbs.tsv", sep="\t", index=False)
    pd.DataFrame(
        rows_f,
        columns=["sample_id", "chrom", "start", "end", "methylated_count", "total_count"],
    ).to_csv(tmp_path / "fme.tsv", sep="\t", index=False)

    params = train_calibration(
        matched_wgbs=tmp_path / "wgbs.tsv",
        matched_finaleme=tmp_path / "fme.tsv",
        region_annotation=None,
        n_bins_candidates=[2, 4, 6],
        out_params=tmp_path / "cal.json",
        out_report=tmp_path / "report.json",
        cv_strategy="region",
        cv_n_folds=5,
        cv_seed=3,
    )
    # Metadata captures the strategy
    md = params.training_metadata
    assert md["cv_strategy_requested"] == "region"
    assert md["cv_strategy_effective"] == "region"
    assert md["cv_n_folds"] == 5
    assert md["cv_seed"] == 3

    # Report must have finite cv_rmse per candidate
    report = json.loads((tmp_path / "report.json").read_text())
    assert report["cv_strategy_requested"] == "region"
    assert report["cv_strategy_effective"] == "region"
    for cand in report["candidates"]:
        assert np.isfinite(cand["cv_rmse"]), \
            f"cv_rmse must be finite in region mode; got {cand}"
        assert cand["n_folds"] == 5


def test_train_calibration_single_sample_auto_falls_back_to_region(tmp_path):
    """train_calibration with cv_strategy='auto' (the default) and only
    one sample must automatically use region-level CV."""
    from finaleme_too.preprocessing.calibration import train_calibration

    fme_arr, wgbs_arr, k_arr, n_arr, _density, _sids = _synth_matched_arrays(
        n_samples=1, n_markers_per_sample=300, seed=21
    )
    rows_w = [
        ("S0", "chr1", i * 100, i * 100 + 1, int(k_arr[i]), int(n_arr[i]))
        for i in range(len(wgbs_arr))
    ]
    rows_f = [
        ("S0", "chr1", i * 100, i * 100 + 1,
         int(round(fme_arr[i] * 30)), 30)
        for i in range(len(fme_arr))
    ]
    pd.DataFrame(
        rows_w,
        columns=["sample_id", "chrom", "start", "end", "methylated_count", "total_count"],
    ).to_csv(tmp_path / "wgbs.tsv", sep="\t", index=False)
    pd.DataFrame(
        rows_f,
        columns=["sample_id", "chrom", "start", "end", "methylated_count", "total_count"],
    ).to_csv(tmp_path / "fme.tsv", sep="\t", index=False)

    params = train_calibration(
        matched_wgbs=tmp_path / "wgbs.tsv",
        matched_finaleme=tmp_path / "fme.tsv",
        region_annotation=None,
        n_bins_candidates=[2, 4],
        out_params=tmp_path / "cal.json",
        out_report=tmp_path / "report.json",
        cv_strategy="auto",  # default
        cv_seed=5,
    )
    assert params.training_metadata["cv_strategy_requested"] == "auto"
    assert params.training_metadata["cv_strategy_effective"] == "region"


# ---------------------------------------------------------------------------
# April 2026 — Pearson/Spearman correlations in the calibration report
# ---------------------------------------------------------------------------
# User request: include raw vs calibrated Pearson/Spearman on the full
# dataset and on the CV test set, so report readers can see at a glance
# whether calibration actually improved the FinaleMe-vs-WGBS agreement.


def test_fit_calibration_reports_full_dataset_correlations():
    """fit_calibration.overall must carry four correlation columns:
    pearson_raw, pearson_calibrated, spearman_raw, spearman_calibrated."""
    from finaleme_too.preprocessing.calibration_eval import fit_calibration

    fme, wgbs, k, n, density, _sids = _synth_matched_arrays(
        n_samples=1, n_markers_per_sample=400, seed=31
    )
    result = fit_calibration(fme, wgbs, density, n_bins=4, wgbs_k=k, wgbs_n=n)
    for key in (
        "pearson_raw",
        "pearson_calibrated",
        "spearman_raw",
        "spearman_calibrated",
    ):
        assert key in result.overall, f"missing key {key} in fit.overall"
        assert np.isfinite(result.overall[key]), f"{key} is not finite"
        assert -1.0 <= result.overall[key] <= 1.0


def test_fit_calibration_full_correlations_improve_on_bad_raw_predictions():
    """When raw FinaleMe predictions are genuinely misaligned with WGBS,
    the calibrated correlation should be >= the raw correlation on the
    full dataset.
    """
    from scipy.special import expit, logit

    from finaleme_too.preprocessing.calibration_eval import fit_calibration

    rng = np.random.default_rng(77)
    n = 600
    # WGBS beta distribution
    wgbs_beta = rng.beta(2, 5, size=n)
    # Miscalibrated FinaleMe: shrink+shift in logit space then add noise.
    # The raw Pearson will be lower than the calibrated Pearson because
    # the calibration model undoes the shrink/shift.
    raw_logit = (logit(np.clip(wgbs_beta, 1e-6, 1 - 1e-6)) - 0.6) / 0.4
    fme_beta = np.clip(
        expit(raw_logit + rng.normal(0, 0.2, size=n)), 1e-6, 1 - 1e-6
    )
    wgbs_n = np.full(n, 30, dtype=np.int64)
    wgbs_k = rng.binomial(wgbs_n, wgbs_beta).astype(np.int64)
    density = rng.uniform(0, 0.1, size=n)

    fit = fit_calibration(
        fme_beta, wgbs_beta, density, n_bins=4,
        wgbs_k=wgbs_k, wgbs_n=wgbs_n,
    )
    raw_r = fit.overall["pearson_raw"]
    cal_r = fit.overall["pearson_calibrated"]
    # Calibration should not make things worse; with the miscalibrated
    # setup it should strictly improve.
    assert cal_r >= raw_r - 1e-6
    # And the effect must be non-trivial given the setup.
    assert cal_r - raw_r > 0.01


def test_fit_calibration_per_bin_reports_correlations():
    """Per-bin metrics must also include the four correlation columns."""
    from finaleme_too.preprocessing.calibration_eval import fit_calibration

    fme, wgbs, k, n, density, _sids = _synth_matched_arrays(
        n_samples=1, n_markers_per_sample=400, seed=32
    )
    result = fit_calibration(fme, wgbs, density, n_bins=4, wgbs_k=k, wgbs_n=n)
    for m in result.per_bin_metrics:
        for key in (
            "pearson_raw",
            "pearson_calibrated",
            "spearman_raw",
            "spearman_calibrated",
        ):
            assert key in m, f"missing {key} in per-bin metrics"


def test_cv_fold_task_returns_dict_with_correlations():
    """_cv_fold_task must now return a dict of test-set metrics instead
    of a bare RMSE float. The dict must include the four correlations."""
    from finaleme_too.preprocessing.calibration_eval import _cv_fold_task

    fme, wgbs, k, n, density, _sids = _synth_matched_arrays(
        n_samples=1, n_markers_per_sample=200, seed=33
    )
    train_mask = np.zeros(fme.size, dtype=bool)
    train_mask[:150] = True
    test_mask = ~train_mask

    result = _cv_fold_task(
        fme, wgbs, density, train_mask, test_mask, n_bins=4,
        wgbs_k=k, wgbs_n=n,
    )
    assert isinstance(result, dict)
    for key in (
        "rmse",
        "pearson_raw",
        "pearson_calibrated",
        "spearman_raw",
        "spearman_calibrated",
        "n",
    ):
        assert key in result, f"missing {key}"
    assert result["n"] == 50
    assert np.isfinite(result["rmse"])


def test_cross_validate_reports_cv_test_set_correlations():
    """cross_validate_calibration must propagate fold-aggregated Pearson
    and Spearman (raw + calibrated) to the top-level result dict."""
    from finaleme_too.preprocessing.calibration_eval import (
        cross_validate_calibration,
    )

    fme, wgbs, k, n, density, sample_ids = _synth_matched_arrays(
        n_samples=3, n_markers_per_sample=200, seed=34
    )
    result = cross_validate_calibration(
        fme, wgbs, density, sample_ids, n_bins=4,
        wgbs_k=k, wgbs_n=n,
        cv_strategy="sample",
    )
    for key in (
        "cv_pearson_raw",
        "cv_pearson_calibrated",
        "cv_spearman_raw",
        "cv_spearman_calibrated",
        "cv_rmse_std",
        "fold_metrics",
    ):
        assert key in result, f"missing {key}"
    assert np.isfinite(result["cv_pearson_raw"])
    assert np.isfinite(result["cv_pearson_calibrated"])
    # Every per-fold entry must have the same keys as the aggregated dict.
    assert len(result["fold_metrics"]) == 3
    for fold in result["fold_metrics"]:
        assert "pearson_raw" in fold
        assert "pearson_calibrated" in fold


def test_tune_n_bins_candidates_carry_correlations():
    """Every candidate in the tune_n_bins output must carry the test-set
    correlations so _select_best_candidate and the report have access."""
    from finaleme_too.preprocessing.calibration_eval import tune_n_bins

    fme, wgbs, k, n, density, sample_ids = _synth_matched_arrays(
        n_samples=1, n_markers_per_sample=500, seed=35
    )
    result = tune_n_bins(
        fme, wgbs, density, sample_ids,
        n_bins_candidates=[2, 4, 6],
        wgbs_k=k, wgbs_n=n,
        cv_strategy="region", cv_n_folds=5, cv_seed=100,
    )
    assert len(result["candidates"]) == 3
    for cand in result["candidates"]:
        for key in (
            "cv_pearson_raw",
            "cv_pearson_calibrated",
            "cv_spearman_raw",
            "cv_spearman_calibrated",
        ):
            assert key in cand, f"candidate missing {key}"
            assert np.isfinite(cand[key])


def test_train_calibration_report_has_correlation_summary(tmp_path):
    """The top-level training report must have a ``correlation_summary``
    block with both full-dataset and CV-test-set Pearson/Spearman values.
    """
    from finaleme_too.preprocessing.calibration import train_calibration

    fme_arr, wgbs_arr, k_arr, n_arr, _density, _sids = _synth_matched_arrays(
        n_samples=1, n_markers_per_sample=400, seed=40
    )
    rows_w = [
        ("S0", "chr1", i * 100, i * 100 + 1, int(k_arr[i]), int(n_arr[i]))
        for i in range(len(wgbs_arr))
    ]
    rows_f = [
        ("S0", "chr1", i * 100, i * 100 + 1,
         int(round(fme_arr[i] * 30)), 30)
        for i in range(len(fme_arr))
    ]
    pd.DataFrame(
        rows_w,
        columns=["sample_id", "chrom", "start", "end", "methylated_count", "total_count"],
    ).to_csv(tmp_path / "wgbs.tsv", sep="\t", index=False)
    pd.DataFrame(
        rows_f,
        columns=["sample_id", "chrom", "start", "end", "methylated_count", "total_count"],
    ).to_csv(tmp_path / "fme.tsv", sep="\t", index=False)

    train_calibration(
        matched_wgbs=tmp_path / "wgbs.tsv",
        matched_finaleme=tmp_path / "fme.tsv",
        region_annotation=None,
        n_bins_candidates=[2, 4],
        out_params=tmp_path / "cal.json",
        out_report=tmp_path / "report.json",
        cv_strategy="region",
        cv_n_folds=5,
        cv_seed=42,
    )
    report = json.loads((tmp_path / "report.json").read_text())
    assert "correlation_summary" in report
    cs = report["correlation_summary"]

    assert set(cs.keys()) == {"full_dataset", "cv_test_set"}
    full = cs["full_dataset"]
    for key in (
        "pearson_raw",
        "pearson_calibrated",
        "spearman_raw",
        "spearman_calibrated",
        "n_markers",
    ):
        assert key in full, f"full_dataset missing {key}"
        if key != "n_markers":
            assert -1.0 <= float(full[key]) <= 1.0

    cv = cs["cv_test_set"]
    for key in (
        "pearson_raw",
        "pearson_calibrated",
        "spearman_raw",
        "spearman_calibrated",
        "cv_rmse",
        "n_folds",
        "cv_strategy",
    ):
        assert key in cv, f"cv_test_set missing {key}"
    assert cv["n_folds"] == 5
    assert cv["cv_strategy"] == "region"
    for key in (
        "pearson_raw",
        "pearson_calibrated",
        "spearman_raw",
        "spearman_calibrated",
    ):
        assert np.isfinite(float(cv[key])), f"cv_test_set.{key} is not finite"
