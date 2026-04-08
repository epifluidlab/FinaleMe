"""Phase A unit tests for the v3 context-dependent binarization modules.

These tests cover the new modules in isolation:

  * ``finaleme_too/preprocessing/_matched_sample_sheet.py`` — region class
    classification + extended compute_region_annotation
  * ``finaleme_too/preprocessing/binarization.py`` — BinarizationParams
    dataclass, apply_binarization, save/load round-trip, 2-level bin
    assignment with NaN fallback
  * ``finaleme_too/preprocessing/binarization_eval.py`` — per-bin fit,
    multi-bin fit, chromosome-blocked CV, tune_n_bins, inference QC
  * ``finaleme_too/core/observation_model_binarization.py`` —
    BinarizationModel builder, reference panel binarization, log-likelihood
    discrimination

Phase A does not touch pipeline.py, cli.py, or any existing test. The
tests in this file exercise only the new modules.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from finaleme_too.config import CoverageTier, MeasurementMode
from finaleme_too.exceptions import InvalidBinarizationError
from finaleme_too.io.methylation_loader import MarkerObservations
from finaleme_too.io.reference_panel import ReferencePanel
from finaleme_too.preprocessing._matched_sample_sheet import (
    DEFAULT_REGION_CLASS_THRESHOLDS,
    REGION_CLASS_ORDER,
    _classify_region,
    compute_region_annotation,
)
from finaleme_too.preprocessing.binarization import (
    STATE_AMBIGUOUS,
    STATE_EXCLUDED,
    STATE_M,
    STATE_U,
    BinarizationParams,
    apply_binarization,
    build_identity_placeholder_params,
    load_default_binarization,
    shipped_default_binarization_path,
)
from finaleme_too.preprocessing.binarization_eval import (
    _chromosome_fold_masks,
    _select_best_binarization_candidate,
    compute_inference_qc,
    cross_validate_binarization,
    fit_binarization,
    fit_binarization_bin,
    tune_n_bins,
)
from finaleme_too.core.observation_model_binarization import (
    BinarizationModel,
    _binarize_reference_panel,
)


# ---------------------------------------------------------------------------
# Region class classification
# ---------------------------------------------------------------------------


def test_classify_region_cgi_shore_shelf_open_sea():
    """Density thresholds assign markers to the right region classes."""
    # Default thresholds: CGI >= 0.08, shore >= 0.04, shelf >= 0.015, else open_sea
    density = np.array([0.10, 0.08, 0.05, 0.04, 0.02, 0.015, 0.005, 0.0])
    classes = _classify_region(density)
    assert classes[0] == "CGI"       # 0.10 >= 0.08
    assert classes[1] == "CGI"       # 0.08 >= 0.08 (boundary inclusive)
    assert classes[2] == "shore"     # 0.05 in [0.04, 0.08)
    assert classes[3] == "shore"     # 0.04 boundary inclusive
    assert classes[4] == "shelf"     # 0.02 in [0.015, 0.04)
    assert classes[5] == "shelf"     # 0.015 boundary inclusive
    assert classes[6] == "open_sea"  # 0.005 < 0.015
    assert classes[7] == "open_sea"  # 0.0 < 0.015


def test_classify_region_nan_density_goes_to_open_sea():
    """Non-finite densities must deterministically route to open_sea
    (np.isfinite treats ±inf as non-finite)."""
    density = np.array([np.nan, -np.inf, np.inf, 0.5])
    classes = _classify_region(density)
    assert classes[0] == "open_sea"  # NaN → open_sea (non-finite)
    assert classes[1] == "open_sea"  # -inf → open_sea (non-finite)
    assert classes[2] == "open_sea"  # +inf → open_sea (non-finite)
    assert classes[3] == "CGI"       # 0.5 >= 0.08 (finite, CGI)


def test_classify_region_custom_thresholds():
    """Custom thresholds override the defaults."""
    density = np.array([0.30, 0.20, 0.10, 0.05])
    custom = {"CGI": 0.25, "shore": 0.15, "shelf": 0.08}
    classes = _classify_region(density, custom)
    assert classes[0] == "CGI"       # 0.30 >= 0.25
    assert classes[1] == "shore"     # 0.20 in [0.15, 0.25)
    assert classes[2] == "shelf"     # 0.10 in [0.08, 0.15)
    assert classes[3] == "open_sea"  # 0.05 < 0.08


def test_compute_region_annotation_emits_region_class_column(tmp_path):
    """compute_region_annotation must now include a region_class column."""
    # Build a tiny CpG index: 100 CpGs every 10bp on chr1
    cpg_index_path = tmp_path / "cpg_index.bed"
    lines = [f"chr1\t{i * 10}\t{i * 10 + 1}" for i in range(100)]
    cpg_index_path.write_text("\n".join(lines) + "\n")

    # Regions: a few markers at different positions
    regions = pd.DataFrame({
        "chrom": ["chr1", "chr1", "chr1"],
        "start": [0, 500, 950],
        "end": [100, 600, 1050],
    })

    annotation = compute_region_annotation(regions, cpg_index_path, window=1000)
    assert "cpg_density" in annotation.columns
    assert "region_class" in annotation.columns
    assert len(annotation) == 3
    # Every row should have a valid class label
    assert set(annotation["region_class"]).issubset(set(REGION_CLASS_ORDER))


# ---------------------------------------------------------------------------
# BinarizationParams save/load + 2-level bin assignment
# ---------------------------------------------------------------------------


def test_binarization_params_save_load_round_trip(tmp_path):
    """Saving then loading BinarizationParams must preserve all fields."""
    params = build_identity_placeholder_params()
    path = tmp_path / "binarization.json"
    params.save(path)

    loaded = BinarizationParams.load(path)
    assert loaded.n_bins == params.n_bins
    assert loaded.n_region_classes == params.n_region_classes
    assert loaded.density_sub_bins_per_class == params.density_sub_bins_per_class
    assert loaded.region_class_order == params.region_class_order
    np.testing.assert_array_equal(loaded.density_edges, params.density_edges)
    np.testing.assert_array_equal(loaded.tau_low, params.tau_low)
    np.testing.assert_array_equal(loaded.tau_high, params.tau_high)
    np.testing.assert_array_equal(loaded.eps_U, params.eps_U)
    np.testing.assert_array_equal(loaded.eps_M, params.eps_M)
    np.testing.assert_array_equal(loaded.usable, params.usable)


def test_binarization_params_load_rejects_missing_keys(tmp_path):
    path = tmp_path / "broken.json"
    path.write_text(json.dumps({"n_bins": 8}))  # missing everything else
    with pytest.raises(InvalidBinarizationError):
        BinarizationParams.load(path)


def test_binarization_params_load_missing_file(tmp_path):
    with pytest.raises(InvalidBinarizationError):
        BinarizationParams.load(tmp_path / "nonexistent.json")


def test_binarization_params_assign_bin_2d_lookup():
    """assign_bin must compute (class_idx × density_sub_bins_per_class + sub_bin)."""
    params = build_identity_placeholder_params()
    assert params.n_bins == 8  # 4 classes × 2 sub-bins
    # Placeholder has interior edges at 1.0 so all realistic densities go to sub-bin 0
    densities = np.array([0.15, 0.05, 0.02, 0.005])
    bin_idx = params.assign_bin(densities)
    # CGI = class 0 → bin 0, shore = class 1 → bin 2, shelf = class 2 → bin 4, open_sea = class 3 → bin 6
    assert bin_idx.tolist() == [0, 2, 4, 6]


def test_binarization_params_assign_bin_nan_density_falls_back_to_open_sea():
    """NaN density → open_sea class, sub-bin 0. This is the v2 bug-4 fix
    pattern re-implemented on the v3 dataclass."""
    params = build_identity_placeholder_params()
    bin_idx = params.assign_bin(np.array([np.nan, np.nan, np.nan]))
    # All NaN → open_sea (class 3) → bin 6 (3*2+0)
    assert bin_idx.tolist() == [6, 6, 6]


def test_binarization_params_assign_bin_unknown_region_class_falls_back():
    """An explicit region_class label not in region_class_order falls back
    to the open_sea class index."""
    params = build_identity_placeholder_params()
    bin_idx = params.assign_bin(
        np.array([0.15, 0.05]),
        region_class=np.array(["unknown_class", "alien"], dtype=object),
    )
    # Fallback to open_sea (class 3) → bin 6
    assert bin_idx.tolist() == [6, 6]


# ---------------------------------------------------------------------------
# apply_binarization
# ---------------------------------------------------------------------------


def _mk_obs(
    sample_id: str,
    predicted_beta: np.ndarray,
    mode: MeasurementMode = MeasurementMode.FINALEME,
) -> MarkerObservations:
    n = predicted_beta.size
    return MarkerObservations(
        sample_id=sample_id,
        chrom=np.array(["chr1"] * n, dtype=object),
        start=np.array([1000 + i * 1000 for i in range(n)], dtype=np.int64),
        end=np.array([1100 + i * 1000 for i in range(n)], dtype=np.int64),
        k=np.array([10] * n, dtype=np.int32),
        n=np.array([20] * n, dtype=np.int32),
        predicted_beta=predicted_beta.astype(np.float32),
        mode=mode,
    )


def test_apply_binarization_classifies_U_M_ambiguous_excluded():
    """The four call states are assigned correctly."""
    params = build_identity_placeholder_params()  # τ_low=0.2, τ_high=0.8, usable=True
    pred = np.array([0.05, 0.95, 0.5, 0.1, 0.9, np.nan], dtype=np.float32)
    obs = _mk_obs("s1", pred)
    binarized = apply_binarization(obs, params, region_annotations=None)

    assert binarized.called_state[0] == STATE_U          # 0.05 < 0.2
    assert binarized.called_state[1] == STATE_M          # 0.95 > 0.8
    assert binarized.called_state[2] == STATE_AMBIGUOUS  # 0.5 in [0.2, 0.8]
    assert binarized.called_state[3] == STATE_U          # 0.1 < 0.2
    assert binarized.called_state[4] == STATE_M          # 0.9 > 0.8
    assert binarized.called_state[5] == STATE_EXCLUDED   # NaN


def test_apply_binarization_respects_usable_flag():
    """Markers in bins with usable=False must be Excluded even when
    predicted_beta is in the U/M range."""
    params = build_identity_placeholder_params()
    # Mark bin 6 (open_sea) as unusable
    params.usable[6] = False
    # 3 markers, all NaN density → all in bin 6 (unusable)
    pred = np.array([0.05, 0.5, 0.95], dtype=np.float32)
    obs = _mk_obs("s1", pred)
    binarized = apply_binarization(obs, params, region_annotations=None)
    # All should be Excluded because bin 6 is unusable
    assert binarized.called_state[0] == STATE_EXCLUDED
    assert binarized.called_state[1] == STATE_EXCLUDED
    assert binarized.called_state[2] == STATE_EXCLUDED


def test_apply_binarization_with_no_region_annotation_does_not_crash():
    """Regression: apply_binarization with region_annotations=None must
    produce a valid output (NaN density → open_sea fallback)."""
    params = build_identity_placeholder_params()
    pred = np.array([0.05, 0.95], dtype=np.float32)
    obs = _mk_obs("s1", pred)
    binarized = apply_binarization(obs, params, region_annotations=None)
    assert binarized.called_state.shape == (2,)
    assert binarized.context_bin.shape == (2,)
    # NaN density → open_sea (class 3) → bin 6
    assert set(binarized.context_bin.tolist()) == {6}


def test_apply_binarization_preserves_k_n_predicted_beta():
    """apply_binarization must NOT modify k, n, or predicted_beta."""
    params = build_identity_placeholder_params()
    pred = np.array([0.05, 0.95, 0.5], dtype=np.float32)
    obs = _mk_obs("s1", pred)
    binarized = apply_binarization(obs, params, region_annotations=None)
    np.testing.assert_array_equal(binarized.k, obs.k)
    np.testing.assert_array_equal(binarized.n, obs.n)
    np.testing.assert_array_equal(binarized.predicted_beta, obs.predicted_beta)


def test_apply_binarization_noop_when_no_predicted_beta():
    """If predicted_beta is None (WGBS mode), apply_binarization is a no-op."""
    params = build_identity_placeholder_params()
    n = 3
    obs = MarkerObservations(
        sample_id="wgbs_s1",
        chrom=np.array(["chr1"] * n, dtype=object),
        start=np.array([1000, 2000, 3000], dtype=np.int64),
        end=np.array([1100, 2100, 3100], dtype=np.int64),
        k=np.array([5, 10, 15], dtype=np.int32),
        n=np.array([20, 20, 20], dtype=np.int32),
        predicted_beta=None,
        mode=MeasurementMode.WGBS,
    )
    result = apply_binarization(obs, params, region_annotations=None)
    assert result is obs  # exact same object — no work done


def test_default_binarization_is_loadable():
    """The shipped default_binarization.json must load as valid params."""
    params = load_default_binarization()
    assert params.n_bins >= 1
    assert params.n_region_classes == 4
    assert len(params.region_class_order) == 4
    assert shipped_default_binarization_path().exists()


# ---------------------------------------------------------------------------
# fit_binarization_bin — grid search recovers known thresholds
# ---------------------------------------------------------------------------


def test_fit_binarization_bin_recovers_good_thresholds():
    """With clean U/M training data, the bin fit must find thresholds that
    classify correctly with low error rates and marks the bin usable."""
    rng = np.random.default_rng(42)
    n = 1000
    # 40% clearly U, 40% clearly M, 20% intermediate (the training-time
    # "ambiguous" band that the fit should learn to exclude).
    truth = np.concatenate([
        rng.uniform(0.0, 0.15, size=400),   # U (< 0.2)
        rng.uniform(0.85, 1.0, size=400),   # M (> 0.8)
        rng.uniform(0.3, 0.7, size=200),    # intermediate
    ])
    # FinaleMe predictions: slightly noisy version of truth
    pred = np.clip(truth + rng.normal(0, 0.05, size=n), 0.0, 1.0)

    fit = fit_binarization_bin(pred, truth)
    assert fit["usable"] is True
    assert fit["accuracy"] > 0.90
    assert fit["eps_U"] < 0.15
    assert fit["eps_M"] < 0.15
    assert fit["n_markers_U"] >= 20
    assert fit["n_markers_M"] >= 20
    # Learned thresholds should be somewhere in the grid range
    assert 0.05 <= fit["tau_low"] <= 0.30
    assert 0.35 <= fit["tau_high"] <= 0.70


def test_fit_binarization_bin_marks_unusable_when_signal_is_poor():
    """When ground truth is uniform (no discrimination possible), the bin
    must be marked unusable regardless of threshold."""
    rng = np.random.default_rng(0)
    n = 500
    truth = rng.uniform(0.3, 0.7, size=n)  # all intermediate
    pred = rng.uniform(0.0, 1.0, size=n)    # random
    fit = fit_binarization_bin(pred, truth)
    assert fit["usable"] is False


def test_fit_binarization_bin_handles_empty_input():
    """Empty input shouldn't crash — returns fallback marked unusable."""
    fit = fit_binarization_bin(np.array([]), np.array([]))
    assert fit["usable"] is False
    assert fit["n_markers_U"] == 0


def test_fit_binarization_bin_respects_max_error_rate():
    """Raising max_error_rate allows noisier bins to be marked usable."""
    rng = np.random.default_rng(1)
    n = 400
    # Noisier dataset: enough U/M but with 15% noise → ε ≈ 0.15
    truth = np.concatenate([
        np.full(200, 0.05),  # U
        np.full(200, 0.95),  # M
    ])
    # Flip labels on 15% of markers
    pred = np.where(
        rng.random(n) < 0.15,
        1.0 - truth,  # flip
        truth,
    )
    pred = pred + rng.normal(0, 0.02, size=n)
    pred = np.clip(pred, 0.01, 0.99)

    # Tight error cap: should be unusable
    fit_tight = fit_binarization_bin(pred, truth, max_error_rate=0.05)
    # Loose error cap: should be usable
    fit_loose = fit_binarization_bin(pred, truth, max_error_rate=0.30)
    assert fit_loose["usable"] is True
    assert fit_tight["usable"] is False


# ---------------------------------------------------------------------------
# fit_binarization (multi-bin)
# ---------------------------------------------------------------------------


def test_fit_binarization_multi_class_distinguishes_usability():
    """Classes with clean training data should be usable; classes with all
    intermediate values should be unusable."""
    rng = np.random.default_rng(7)
    # 4 classes, 400 markers each. CGI and shore have clean U/M data;
    # shelf and open_sea have intermediate values only (hard to classify).
    density_per_class = {"CGI": 0.15, "shore": 0.05, "shelf": 0.02, "open_sea": 0.005}
    all_pred = []
    all_truth = []
    all_density = []
    for cls, d in density_per_class.items():
        if cls in ("CGI", "shore"):
            truth = np.concatenate([
                np.full(200, 0.05),  # U
                np.full(200, 0.95),  # M
            ])
            pred = truth + rng.normal(0, 0.03, size=400)
        else:
            truth = rng.uniform(0.3, 0.7, size=400)
            pred = rng.uniform(0.0, 1.0, size=400)
        all_pred.append(np.clip(pred, 0.01, 0.99))
        all_truth.append(truth)
        all_density.append(np.full(400, d))
    pred = np.concatenate(all_pred)
    truth = np.concatenate(all_truth)
    density = np.concatenate(all_density)

    result = fit_binarization(
        predicted=pred,
        truth_beta=truth,
        cpg_density=density,
        density_sub_bins_per_class=1,  # 1 sub-bin per class → 4 bins total
    )
    assert result.params.n_bins == 4
    # CGI (bin 0) and shore (bin 1) should be usable
    cgi_bin = result.per_bin_metrics[0]
    shore_bin = result.per_bin_metrics[1]
    shelf_bin = result.per_bin_metrics[2]
    open_sea_bin = result.per_bin_metrics[3]
    assert cgi_bin["usable"] is True, f"CGI should be usable: {cgi_bin}"
    assert shore_bin["usable"] is True, f"shore should be usable: {shore_bin}"
    assert shelf_bin["usable"] is False, f"shelf should be unusable: {shelf_bin}"
    assert open_sea_bin["usable"] is False, f"open_sea should be unusable: {open_sea_bin}"


def test_fit_binarization_records_train_fractions():
    """train_fraction_U and train_fraction_M should sum to 1 per bin (when
    the bin has any called markers)."""
    rng = np.random.default_rng(3)
    n = 800
    density = np.full(n, 0.15)  # all CGI
    # 30% U, 70% M → after threshold learning, the relative ordering should
    # still hold: train_fraction_M > train_fraction_U since M is more common.
    truth = np.concatenate([np.full(240, 0.05), np.full(560, 0.95)])
    pred = truth + rng.normal(0, 0.03, size=n)
    pred = np.clip(pred, 0.01, 0.99)

    result = fit_binarization(
        predicted=pred, truth_beta=truth, cpg_density=density,
        density_sub_bins_per_class=1,
    )
    params = result.params
    # CGI bin (class 0, sub-bin 0) = bin 0
    assert abs(params.train_fraction_U[0] + params.train_fraction_M[0] - 1.0) < 1e-9
    # M should dominate (we trained with 70% M). The exact U fraction depends
    # on the grid-search outcome — it may be less than 30% if the fit tightens
    # τ_low to reject noisy U calls — but M must still dominate.
    assert params.train_fraction_M[0] > params.train_fraction_U[0]
    assert params.train_fraction_M[0] > 0.6  # still roughly dominant


# ---------------------------------------------------------------------------
# Chromosome-blocked CV
# ---------------------------------------------------------------------------


def test_chromosome_blocked_cv_no_leakage():
    """Train and test folds must never share a chromosome."""
    chrom = np.concatenate([
        np.repeat("chr1", 100),
        np.repeat("chr2", 100),
        np.repeat("chr3", 100),
        np.repeat("chrX", 100),
    ])
    masks = _chromosome_fold_masks(chrom, n_folds=4, seed=0)
    assert len(masks) == 4
    for train_mask, test_mask in masks:
        train_chroms = set(chrom[train_mask].tolist())
        test_chroms = set(chrom[test_mask].tolist())
        assert len(train_chroms & test_chroms) == 0, "chromosome leakage!"


def test_chromosome_blocked_cv_seed_reproducibility():
    """Same seed → same fold assignment."""
    chrom = np.array([f"chr{i}" for i in range(22)] * 20, dtype=object)
    masks_a = _chromosome_fold_masks(chrom, n_folds=5, seed=42)
    masks_b = _chromosome_fold_masks(chrom, n_folds=5, seed=42)
    assert len(masks_a) == len(masks_b)
    for (ta, tea), (tb, teb) in zip(masks_a, masks_b):
        np.testing.assert_array_equal(ta, tb)
        np.testing.assert_array_equal(tea, teb)


def test_chromosome_blocked_cv_rejects_single_fold():
    with pytest.raises(InvalidBinarizationError):
        _chromosome_fold_masks(np.array(["chr1"] * 10), n_folds=1)


def test_cross_validate_binarization_returns_finite_accuracy():
    """End-to-end CV on a well-separated dataset must produce finite metrics."""
    rng = np.random.default_rng(5)
    n_per_chrom = 500
    chrom = np.repeat([f"chr{i}" for i in range(4)], n_per_chrom)
    density = np.full(n_per_chrom * 4, 0.15)  # all CGI
    truth = rng.choice([0.05, 0.95], size=n_per_chrom * 4)
    pred = np.clip(truth + rng.normal(0, 0.04, size=n_per_chrom * 4), 0, 1)

    cv_result = cross_validate_binarization(
        predicted=pred,
        truth_beta=truth,
        cpg_density=density,
        chrom=chrom,
        density_sub_bins_per_class=1,
        n_folds=4,
        seed=0,
    )
    assert cv_result["n_folds"] == 4
    assert np.isfinite(cv_result["cv_accuracy"])
    assert cv_result["cv_accuracy"] > 0.85
    assert len(cv_result["fold_metrics"]) == 4


# ---------------------------------------------------------------------------
# tune_n_bins + _select_best_binarization_candidate
# ---------------------------------------------------------------------------


def test_tune_n_bins_picks_by_score_B():
    """tune_n_bins must pick the candidate with the highest score(B)."""
    rng = np.random.default_rng(11)
    n_per_chrom = 500
    chrom = np.repeat([f"chr{i}" for i in range(4)], n_per_chrom)
    density = np.full(n_per_chrom * 4, 0.15)
    truth = rng.choice([0.05, 0.95], size=n_per_chrom * 4)
    pred = np.clip(truth + rng.normal(0, 0.05, size=n_per_chrom * 4), 0, 1)

    result = tune_n_bins(
        predicted=pred,
        truth_beta=truth,
        cpg_density=density,
        chrom=chrom,
        n_bins_candidates=[4, 8],  # 1 or 2 sub-bins per class
        n_folds=4,
        seed=0,
    )
    assert result["selected_n_bins"] in (4, 8)
    # The selected candidate must have the highest score among the candidates
    best_score = max(
        c.get("score", -np.inf) for c in result["candidates"]
        if np.isfinite(c.get("score", np.nan))
    )
    selected = next(
        c for c in result["candidates"]
        if c["n_bins"] == result["selected_n_bins"]
    )
    assert selected["score"] == best_score


def test_tune_n_bins_rounds_B_up_to_multiple_of_n_classes():
    """B that isn't divisible by n_region_classes (4) is rounded up and logged."""
    rng = np.random.default_rng(2)
    n = 1000
    chrom = np.repeat([f"chr{i}" for i in range(4)], 250)
    density = np.full(n, 0.15)
    truth = rng.choice([0.05, 0.95], size=n)
    pred = np.clip(truth + rng.normal(0, 0.04, size=n), 0, 1)

    result = tune_n_bins(
        predicted=pred,
        truth_beta=truth,
        cpg_density=density,
        chrom=chrom,
        n_bins_candidates=[5, 7],  # both NOT divisible by 4
        n_folds=4,
        seed=0,
    )
    # 5 → 8 (ceil(5/4)*4), 7 → 8 (ceil(7/4)*4)
    candidate_Bs = sorted([c["n_bins"] for c in result["candidates"]])
    assert candidate_Bs == [8, 8]


def test_select_best_falls_back_to_cv_accuracy_when_scores_nan():
    """When all score(B) are NaN, fall back to the highest cv_accuracy."""
    candidates = [
        {"n_bins": 4, "score": np.nan, "cv_accuracy": 0.9, "in_sample_accuracy": 0.95},
        {"n_bins": 8, "score": np.nan, "cv_accuracy": 0.7, "in_sample_accuracy": 0.85},
    ]
    result = _select_best_binarization_candidate(candidates)
    assert result["selected_n_bins"] == 4  # higher cv_accuracy


def test_select_best_falls_back_to_in_sample_when_all_else_nan():
    candidates = [
        {"n_bins": 4, "score": np.nan, "cv_accuracy": np.nan, "in_sample_accuracy": 0.80},
        {"n_bins": 8, "score": np.nan, "cv_accuracy": np.nan, "in_sample_accuracy": 0.95},
    ]
    result = _select_best_binarization_candidate(candidates)
    assert result["selected_n_bins"] == 8


def test_select_best_returns_first_when_all_metrics_nan():
    candidates = [
        {"n_bins": 4, "score": np.nan, "cv_accuracy": np.nan, "in_sample_accuracy": np.nan},
        {"n_bins": 8, "score": np.nan, "cv_accuracy": np.nan, "in_sample_accuracy": np.nan},
    ]
    result = _select_best_binarization_candidate(candidates)
    assert result["selected_n_bins"] == 4


# ---------------------------------------------------------------------------
# compute_inference_qc
# ---------------------------------------------------------------------------


def test_inference_qc_passes_well_balanced_sample():
    """A sample with markers spread across all 8 bins and a reasonable U/M
    ratio must get PASS."""
    params = build_identity_placeholder_params()
    n_per_bin = 20
    n_bins = 8
    # 20 markers per bin, 50/50 U/M each
    called_state = np.concatenate([
        np.array([STATE_U] * 10 + [STATE_M] * 10, dtype=np.uint8)
        for _ in range(n_bins)
    ])
    context_bin = np.concatenate([
        np.full(n_per_bin, b, dtype=np.int64) for b in range(n_bins)
    ])
    qc = compute_inference_qc(called_state, context_bin, params)
    assert qc["flag"] == "PASS"
    assert qc["fraction_called"] == 1.0
    assert qc["bin_balance"] == 1.0


def test_inference_qc_flags_all_ambiguous_as_fail():
    """Every marker Ambiguous → fraction_called = 0 → FAIL."""
    params = build_identity_placeholder_params()
    n = 100
    called_state = np.full(n, STATE_AMBIGUOUS, dtype=np.uint8)
    context_bin = np.zeros(n, dtype=np.int64)
    qc = compute_inference_qc(called_state, context_bin, params)
    assert qc["fraction_called"] == 0.0
    assert qc["flag"] == "FAIL"


def test_inference_qc_flags_distribution_shift():
    """A sample that's all-U when training was 50/50 must fire a KL warning."""
    params = build_identity_placeholder_params()
    n_per_bin = 20
    n_bins = 8
    # All U everywhere → strongly shifted from the 50/50 training fractions
    called_state = np.full(n_per_bin * n_bins, STATE_U, dtype=np.uint8)
    context_bin = np.concatenate([
        np.full(n_per_bin, b, dtype=np.int64) for b in range(n_bins)
    ])
    qc = compute_inference_qc(called_state, context_bin, params)
    # KL divergence of (1, 0) vs (0.5, 0.5) ≈ log(2) ≈ 0.693 → WARN
    assert qc["state_distribution_kl"] > 0.5
    assert qc["flag"] in ("WARN", "FAIL")


def test_inference_qc_empty_input_returns_fail():
    params = build_identity_placeholder_params()
    qc = compute_inference_qc(
        np.array([], dtype=np.uint8), np.array([], dtype=np.int64), params
    )
    assert qc["flag"] == "FAIL"
    assert qc["fraction_called"] == 0.0


# ---------------------------------------------------------------------------
# BinarizationModel (observation model builder)
# ---------------------------------------------------------------------------


def test_binarize_reference_panel_soft_intermediate():
    """The soft rule: r < 0.2 → 0, r > 0.8 → 1, else unchanged."""
    ref = np.array([[0.05, 0.5, 0.95], [0.15, 0.3, 0.85]], dtype=np.float32)
    result = _binarize_reference_panel(ref)
    # The clamped values are exact
    assert result[0, 0] == 0.0     # 0.05 < 0.2 → 0
    assert result[0, 2] == 1.0     # 0.95 > 0.8 → 1
    assert result[1, 0] == 0.0     # 0.15 < 0.2 → 0
    assert result[1, 2] == 1.0     # 0.85 > 0.8 → 1
    # Intermediate values are preserved with float32→float64 rounding tolerance
    assert abs(result[0, 1] - 0.5) < 1e-6
    assert abs(result[1, 1] - 0.3) < 1e-6


def _mk_reference(K: int, M: int) -> ReferencePanel:
    """Build a reference where each cell type owns one U marker."""
    methy = np.full((M, K), 0.9, dtype=np.float32)  # default: all M
    for j in range(min(K, M)):
        methy[j, j] = 0.05  # marker j is U for cell type j, M for others
    return ReferencePanel(
        chrom=np.array(["chr1"] * M, dtype=object),
        start=np.array([1000 + i * 1000 for i in range(M)], dtype=np.int64),
        end=np.array([1100 + i * 1000 for i in range(M)], dtype=np.int64),
        cell_types=[f"CT{i}" for i in range(K)],
        methylation=methy,
        coverage=None,
    )


def test_binarization_model_log_likelihood_discriminates_correct_w():
    """The correct mixture (all weight on the generating cell type) must
    have a higher log-likelihood than any other mixture."""
    K = 4
    M = 6
    reference = _mk_reference(K, M)
    params = build_identity_placeholder_params()

    # Sample looks like "pure CT0": U at marker 0, M at markers 1-3, rest ambiguous/M
    pred = np.array([0.05, 0.95, 0.95, 0.95, 0.5, 0.95], dtype=np.float32)
    obs = _mk_obs("pure_ct0", pred)
    binarized = apply_binarization(obs, params, region_annotations=None)

    model = BinarizationModel().build(binarized, params, reference, tier=CoverageTier.HIGH)
    assert model.n_markers == 5  # only the ambiguous marker is filtered out

    w_ct0 = np.array([1.0, 0.0, 0.0, 0.0, 0.0])
    w_ct1 = np.array([0.0, 1.0, 0.0, 0.0, 0.0])
    w_uniform = np.full(5, 0.2)

    ll_ct0 = model.total_log_likelihood(w_ct0)
    ll_ct1 = model.total_log_likelihood(w_ct1)
    ll_uniform = model.total_log_likelihood(w_uniform)

    assert ll_ct0 > ll_ct1, f"correct w={w_ct0} should beat wrong w={w_ct1}: {ll_ct0} vs {ll_ct1}"
    assert ll_ct0 > ll_uniform


def test_binarization_model_filters_ambiguous_and_excluded():
    """BinarizationModel.build must drop markers whose state is Ambiguous
    or whose bin is unusable."""
    K = 3
    M = 4
    reference = _mk_reference(K, M)
    params = build_identity_placeholder_params()

    # 2 called (U/M), 1 ambiguous, 1 NaN (→ Excluded)
    pred = np.array([0.05, 0.95, 0.5, np.nan], dtype=np.float32)
    obs = _mk_obs("s1", pred)
    binarized = apply_binarization(obs, params, region_annotations=None)

    model = BinarizationModel().build(binarized, params, reference)
    assert model.n_markers == 2  # only the 2 called markers pass


def test_binarization_model_raises_when_no_binarization_applied():
    """Trying to build a BinarizationObservationModel before calling
    apply_binarization must raise a clear error."""
    K = 3
    M = 3
    reference = _mk_reference(K, M)
    params = build_identity_placeholder_params()

    # No apply_binarization step → obs.called_state is None
    obs = _mk_obs("raw", np.array([0.05, 0.5, 0.95], dtype=np.float32))
    with pytest.raises(ValueError, match="called_state"):
        BinarizationModel().build(obs, params, reference)


def test_binarization_model_empty_when_nothing_valid():
    """When every marker is Ambiguous, the model is empty (not an error)."""
    K = 3
    M = 3
    reference = _mk_reference(K, M)
    params = build_identity_placeholder_params()

    pred = np.array([0.5, 0.5, 0.5], dtype=np.float32)  # all ambiguous
    obs = _mk_obs("all_amb", pred)
    binarized = apply_binarization(obs, params, region_annotations=None)

    model = BinarizationModel().build(binarized, params, reference)
    assert model.n_markers == 0
    assert model.coef.shape == (0, K + 1)
