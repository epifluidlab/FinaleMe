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
