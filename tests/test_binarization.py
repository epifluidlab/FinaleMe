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


def test_apply_binarization_hard_threshold_has_no_ambiguous():
    """Hard-threshold mode uses <t => U and >=t => M with no Ambiguous calls."""
    params = build_identity_placeholder_params()  # all bins usable
    pred = np.array([0.09, 0.10, 0.11, np.nan], dtype=np.float32)
    obs = _mk_obs("s1", pred)
    binarized = apply_binarization(
        obs,
        params,
        region_annotations=None,
        hard_threshold=0.10,
    )

    assert binarized.called_state[0] == STATE_U
    assert binarized.called_state[1] == STATE_M  # boundary: >= 0.10
    assert binarized.called_state[2] == STATE_M
    assert binarized.called_state[3] == STATE_EXCLUDED
    assert int(np.sum(binarized.called_state == STATE_AMBIGUOUS)) == 0


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


def test_binarize_reference_panel_hard_threshold():
    """Hard-threshold mode uses <t => 0 and >=t => 1 on finite entries."""
    ref = np.array([[0.05, 0.1, 0.95], [0.099, 0.101, np.nan]], dtype=np.float32)
    result = _binarize_reference_panel(ref, hard_threshold=0.1)
    assert result[0, 0] == 0.0
    assert result[0, 1] == 1.0  # boundary: >= 0.1
    assert result[0, 2] == 1.0
    assert result[1, 0] == 0.0
    assert result[1, 1] == 1.0
    assert np.isnan(result[1, 2])


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


# ---------------------------------------------------------------------------
# Phase B integration tests: MLE solver dispatch + p_goodness + pipeline
# ---------------------------------------------------------------------------


def _mk_richer_reference(K: int, n_per_ct: int) -> ReferencePanel:
    """Reference with ``n_per_ct`` own-U markers per cell type.

    Total M = K * n_per_ct. Markers ``[j*n_per_ct, (j+1)*n_per_ct)`` are
    CT_j's own-U markers (r=0.05 for CT_j, r=0.9 for all others). This
    gives the binarization solver enough discriminative evidence per cell
    type to recover known mixtures cleanly.
    """
    M = K * n_per_ct
    methy = np.full((M, K), 0.9, dtype=np.float32)
    for j in range(K):
        for k in range(n_per_ct):
            methy[j * n_per_ct + k, j] = 0.05
    return ReferencePanel(
        chrom=np.array(["chr1"] * M, dtype=object),
        start=np.array([1000 + i * 1000 for i in range(M)], dtype=np.int64),
        end=np.array([1100 + i * 1000 for i in range(M)], dtype=np.int64),
        cell_types=[f"CT{j}" for j in range(K)],
        methylation=methy,
        coverage=None,
    )


def _mk_pure_celltype_pred(K: int, n_per_ct: int, target_ct: int) -> np.ndarray:
    """Build a predicted_beta vector for a pure-target_ct sample.

    Markers 0..(K*n_per_ct - 1). At target_ct's own-U markers we predict
    U (0.05); at every other cell type's own-U markers we predict M (0.95).
    """
    M = K * n_per_ct
    pred = np.full(M, 0.95, dtype=np.float32)
    for k in range(n_per_ct):
        pred[target_ct * n_per_ct + k] = 0.05
    return pred


def test_mle_solver_dispatches_on_model_kind():
    """MLEDeconvolver.solve picks the right solver branch based on the
    observation model type — BinarizationObservationModel goes to
    _solve_binarization, ObservationModel goes to _solve_betabinom.
    """
    from finaleme_too.core.deconvolution import MLEDeconvolver
    from finaleme_too.core.observation_model_binarization import (
        BinarizationObservationModel,
    )

    K = 3
    n_per_ct = 4
    reference = _mk_richer_reference(K, n_per_ct)
    params = build_identity_placeholder_params()

    pred = _mk_pure_celltype_pred(K, n_per_ct, target_ct=0)
    obs = _mk_obs("pure_ct0", pred)
    binarized = apply_binarization(obs, params, region_annotations=None)

    model = BinarizationModel().build(binarized, params, reference)
    assert isinstance(model, BinarizationObservationModel)

    solver = MLEDeconvolver()
    w_hat = solver.solve(model, reference)
    assert w_hat.shape == (K + 1,)
    assert abs(w_hat.sum() - 1.0) < 1e-6
    # CT0 should be the dominant cell type
    assert w_hat.argmax() == 0
    assert w_hat[0] > 0.7


def test_mle_solver_binarization_recovers_50_50_mixture():
    """A 50/50 mixture of two cell types should produce a roughly balanced
    proportion estimate via the binarization solver."""
    from finaleme_too.core.deconvolution import MLEDeconvolver

    K = 4
    n_per_ct = 4
    reference = _mk_richer_reference(K, n_per_ct)
    params = build_identity_placeholder_params()

    # Mixture: half U at CT0's markers, half U at CT1's markers, M everywhere else.
    M = K * n_per_ct
    pred = np.full(M, 0.95, dtype=np.float32)
    # CT0's own-U markers → U (calling CT0)
    for k in range(n_per_ct):
        pred[0 * n_per_ct + k] = 0.05
    # CT1's own-U markers → U (calling CT1)
    for k in range(n_per_ct):
        pred[1 * n_per_ct + k] = 0.05

    obs = _mk_obs("mix50_50", pred)
    binarized = apply_binarization(obs, params, region_annotations=None)
    model = BinarizationModel().build(binarized, params, reference)
    solver = MLEDeconvolver()
    w_hat = solver.solve(model, reference)
    # CT0 and CT1 should both be elevated
    assert w_hat[0] > 0.2
    assert w_hat[1] > 0.2
    # CT2/CT3 should be near zero
    assert w_hat[2] < 0.1
    assert w_hat[3] < 0.1


def test_mle_solver_binarization_handles_marker_subset():
    """The bootstrap path passes a marker_subset to the solver. Verify the
    binarization solver subset-indexes coef + weights correctly.
    """
    from finaleme_too.core.deconvolution import MLEDeconvolver

    K = 3
    n_per_ct = 4
    reference = _mk_richer_reference(K, n_per_ct)
    params = build_identity_placeholder_params()
    pred = _mk_pure_celltype_pred(K, n_per_ct, target_ct=0)
    obs = _mk_obs("pure_ct0", pred)
    binarized = apply_binarization(obs, params, region_annotations=None)
    model = BinarizationModel().build(binarized, params, reference)

    solver = MLEDeconvolver()
    # Subset to just CT0's own-U markers (still pure CT0 evidence)
    subset = np.arange(n_per_ct)  # markers 0..3 (all CT0-U)
    w_subset = solver.solve(model, reference, marker_subset=subset)
    assert abs(w_subset.sum() - 1.0) < 1e-6
    assert w_subset[0] > 0.7  # CT0 still dominant on the subset


def test_compute_p_goodness_finaleme_binomial_concordance():
    """compute_p_goodness on a BinarizationObservationModel uses the
    binomial concordance test against eps_U/eps_M, not chi-square. With
    a perfect-fit sample the p-value should be high (>= 0.05 by the
    'high = good' convention)."""
    from finaleme_too.core.deconvolution import MLEDeconvolver
    from finaleme_too.core.reliability import compute_p_goodness

    K = 3
    n_per_ct = 4
    reference = _mk_richer_reference(K, n_per_ct)
    params = build_identity_placeholder_params()

    pred = _mk_pure_celltype_pred(K, n_per_ct, target_ct=0)
    obs = _mk_obs("pure_ct0", pred)
    binarized = apply_binarization(obs, params, region_annotations=None)
    model = BinarizationModel().build(binarized, params, reference)

    solver = MLEDeconvolver()
    w_hat = solver.solve(model, reference)

    # P_goodness for the dominant cell type should be HIGH (good fit)
    p_ct0 = compute_p_goodness(
        w_hat=w_hat,
        reference_methylation=reference.methylation,
        observation=model,
        cell_type_index=0,
        binarizer=params,
    )
    # All discriminative markers for CT0 should agree with the expected U
    # state under the estimated mixture, so the binomial test should NOT
    # reject (high p-value).
    assert p_ct0 > 0.05


def test_compute_p_goodness_finaleme_low_for_wrong_mixture():
    """If we lie to compute_p_goodness about the mixture (claim it's
    pure CT2 when the data really came from pure CT0), the binomial test
    should reject (low p-value)."""
    from finaleme_too.core.reliability import compute_p_goodness

    K = 3
    n_per_ct = 4
    reference = _mk_richer_reference(K, n_per_ct)
    params = build_identity_placeholder_params()

    pred = _mk_pure_celltype_pred(K, n_per_ct, target_ct=0)
    obs = _mk_obs("pure_ct0", pred)
    binarized = apply_binarization(obs, params, region_annotations=None)
    model = BinarizationModel().build(binarized, params, reference)

    # Wrong w_hat: claim it's pure CT2
    w_wrong = np.array([0.0, 0.0, 1.0, 0.0])  # K+1 = 4

    p_wrong = compute_p_goodness(
        w_hat=w_wrong,
        reference_methylation=reference.methylation,
        observation=model,
        cell_type_index=2,
        binarizer=params,
    )
    # Many of CT2's discriminative markers will disagree with the actual
    # state (since this is really CT0 data), so the test should reject.
    # The exact value depends on the test setup but it should be < 0.5
    # which is much lower than the perfect-fit p ≈ 1.0.
    assert p_wrong < 0.5


def test_pipeline_end_to_end_binarization_recovers_known_mixture(tmp_path):
    """End-to-end TOOPipeline run with binarization=BinarizationParams
    must recover the correct cell-type proportions for synthetic
    pure-cell-type samples."""
    from finaleme_too.config import TOOConfig
    from finaleme_too.io.marker_regions import MarkerRegions
    from finaleme_too.io.sample_sheet import Sample, SampleSheet
    from finaleme_too.pipeline import TOOPipeline

    K = 4
    M = 8
    methy = np.full((M, K), 0.9, dtype=np.float32)
    for j in range(K):
        methy[2 * j, j] = 0.05
        methy[2 * j + 1, j] = 0.05
    chrom = np.array(["chr1"] * M, dtype=object)
    starts = np.array([1000 + i * 1000 for i in range(M)], dtype=np.int64)
    ends = starts + 100
    reference = ReferencePanel(
        chrom=chrom,
        start=starts,
        end=ends,
        cell_types=[f"CT{i}" for i in range(K)],
        methylation=methy,
        coverage=None,
    )
    marker_regions = MarkerRegions(
        chrom=chrom, start=starts, end=ends, marker_name=None,
    )

    # Build 2 FinaleMe samples by writing prediction.bed files
    samples = []
    for ct_idx in range(2):
        bed_path = tmp_path / f"sample_ct{ct_idx}.prediction.bed"
        rows = []
        for i in range(M):
            pred_beta = 0.05 if i in (2 * ct_idx, 2 * ct_idx + 1) else 0.95
            methy_count = round(pred_beta * 20)
            total_count = 20
            methy_pct = pred_beta * 100
            rows.append(
                f"{chrom[i]}\t{starts[i]}\t{ends[i]}\t{methy_pct:.4f}"
                f"\t{methy_count}\t{total_count}\t0\t0\t0"
            )
        bed_path.write_text("\n".join(rows) + "\n")
        samples.append(
            Sample(
                sample_id=f"S{ct_idx}",
                methylation_file=bed_path,
                mode=MeasurementMode.FINALEME,
                input_format="finaleme_bed",
                group="A",
            )
        )
    sheet = SampleSheet(
        samples=samples,
        raw_table=pd.DataFrame([{"sample_id": s.sample_id} for s in samples]),
    )

    cfg = TOOConfig()
    cfg.threads = 1
    cfg.uncertainty.n_bootstrap = 5
    cfg.uncertainty.seed = 0
    cfg.coverage.tier_high = 0.0
    cfg.coverage.tier_low = -1.0
    cfg.markers.n_per_type = 0

    binarization = build_identity_placeholder_params()

    pipeline = TOOPipeline(
        config=cfg,
        binarization=binarization,
    )

    out_dir = tmp_path / "out"
    out_dir.mkdir(exist_ok=True)
    cohort = pipeline.run(sheet, reference, marker_regions, out_dir)

    assert len(cohort.samples) == 2
    s0, s1 = cohort.samples
    # Pure CT0 sample → CT0 should be ≥0.7
    assert s0.proportions[0] > 0.7
    # Pure CT1 sample → CT1 should be ≥0.7
    assert s1.proportions[1] > 0.7


def test_pipeline_v2_calibration_path_still_works_when_binarization_none(tmp_path):
    """The v2 calibration path must continue to work when binarization=None.
    This is the safety guarantee that the migration is non-breaking."""
    from finaleme_too.config import TOOConfig
    from finaleme_too.io.marker_regions import MarkerRegions
    from finaleme_too.io.sample_sheet import Sample, SampleSheet
    from finaleme_too.pipeline import TOOPipeline

    K = 3
    M = 6
    methy = np.full((M, K), 0.9, dtype=np.float32)
    for j in range(K):
        methy[2 * j, j] = 0.05
        methy[2 * j + 1, j] = 0.05
    chrom = np.array(["chr1"] * M, dtype=object)
    starts = np.array([1000 + i * 1000 for i in range(M)], dtype=np.int64)
    ends = starts + 100
    reference = ReferencePanel(
        chrom=chrom,
        start=starts,
        end=ends,
        cell_types=[f"CT{i}" for i in range(K)],
        methylation=methy,
        coverage=None,
    )
    marker_regions = MarkerRegions(
        chrom=chrom, start=starts, end=ends, marker_name=None,
    )

    # Build 1 WGBS sample (the simplest path that doesn't need calibration
    # or binarization at all)
    bed_path = tmp_path / "wgbs_sample.bed"
    rows = []
    for i in range(M):
        pct = 5.0 if i in (0, 1) else 95.0  # CT0 pattern
        total = 20
        rows.append(
            f"{chrom[i]}\t{starts[i]}\t{ends[i]}\t.\t0\t+\t{pct:.4f}\t{total}"
        )
    bed_path.write_text("\n".join(rows) + "\n")

    sample = Sample(
        sample_id="wgbs_s0",
        methylation_file=bed_path,
        mode=MeasurementMode.WGBS,
        input_format="bissnp_6plus2",
        group="A",
    )
    sheet = SampleSheet(
        samples=[sample],
        raw_table=pd.DataFrame([{"sample_id": "wgbs_s0"}]),
    )

    cfg = TOOConfig()
    cfg.threads = 1
    cfg.uncertainty.n_bootstrap = 5
    cfg.uncertainty.seed = 0
    cfg.coverage.tier_high = 0.0
    cfg.coverage.tier_low = -1.0
    cfg.markers.n_per_type = 0

    # binarization=None — v2 path should still run
    pipeline = TOOPipeline(config=cfg, binarization=None)
    out_dir = tmp_path / "out"
    out_dir.mkdir(exist_ok=True)
    cohort = pipeline.run(sheet, reference, marker_regions, out_dir)

    assert len(cohort.samples) == 1
    # Beta-binomial path on this synthetic data should still recover CT0
    assert cohort.samples[0].proportions[0] > 0.7


def test_load_optional_binarization_returns_default_when_no_path(tmp_path):
    """load_optional_binarization should return the shipped default
    placeholder when use_default=True and no explicit path is given."""
    from finaleme_too.config import TOOConfig
    from finaleme_too.pipeline import load_optional_binarization

    cfg = TOOConfig()
    params = load_optional_binarization(cfg, explicit_path=None, use_default=True)
    assert params is not None
    assert params.n_bins == 8


def test_load_optional_binarization_returns_none_when_default_disabled(tmp_path):
    """When use_default=False and no explicit path, return None so callers
    can fall through to a different path."""
    from finaleme_too.config import TOOConfig
    from finaleme_too.pipeline import load_optional_binarization

    cfg = TOOConfig()
    params = load_optional_binarization(cfg, explicit_path=None, use_default=False)
    assert params is None


def test_load_optional_binarization_loads_from_explicit_path(tmp_path):
    """An explicit path takes precedence and loads the params from it."""
    from finaleme_too.pipeline import load_optional_binarization
    from finaleme_too.config import TOOConfig

    placeholder = build_identity_placeholder_params()
    custom_path = tmp_path / "custom.json"
    placeholder.save(custom_path)

    cfg = TOOConfig()
    loaded = load_optional_binarization(
        cfg, explicit_path=str(custom_path), use_default=False
    )
    assert loaded is not None
    assert loaded.n_bins == placeholder.n_bins
    np.testing.assert_array_equal(loaded.tau_low, placeholder.tau_low)


# ---------------------------------------------------------------------------
# Phase D integration tests: CLI + config migration
# ---------------------------------------------------------------------------


def test_tooconfig_has_binarization_section():
    """TOOConfig gains a ``binarization`` subsection alongside the legacy
    ``calibration`` subsection. Both default-construct without errors."""
    from finaleme_too.config import BinarizationConfig, TOOConfig

    cfg = TOOConfig()
    assert isinstance(cfg.binarization, BinarizationConfig)
    assert cfg.binarization.n_context_bins == 8
    assert cfg.binarization.max_error_rate == 0.15
    assert cfg.binarization.cv_method == "chromosome_blocked"
    assert cfg.binarization.use_default is True


def test_tooconfig_from_yaml_accepts_binarization_section(tmp_path):
    """Loading a YAML file with a ``binarization:`` section populates the
    new subsection."""
    from finaleme_too.config import TOOConfig

    yaml_path = tmp_path / "v3_config.yaml"
    yaml_path.write_text(
        "binarization:\n"
        "  binarization_file: /tmp/custom_binarization.json\n"
        "  n_context_bins: 12\n"
        "  max_error_rate: 0.10\n"
        "  cv_n_folds: 5\n"
    )
    cfg = TOOConfig.from_yaml(yaml_path)
    assert cfg.binarization.binarization_file == "/tmp/custom_binarization.json"
    assert cfg.binarization.n_context_bins == 12
    assert cfg.binarization.max_error_rate == 0.10
    assert cfg.binarization.cv_n_folds == 5


def test_tooconfig_from_yaml_warns_on_legacy_calibration_section(tmp_path):
    """Loading a YAML file with a legacy ``calibration:`` section emits a
    DeprecationWarning and maps the known v2 keys into the new
    ``binarization:`` subsection."""
    import warnings
    from finaleme_too.config import TOOConfig

    yaml_path = tmp_path / "v2_config.yaml"
    yaml_path.write_text(
        "calibration:\n"
        "  calibration_file: /tmp/legacy_cal.json\n"
        "  n_density_bins: 6\n"
    )
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        cfg = TOOConfig.from_yaml(yaml_path)
        deprecations = [w for w in captured if issubclass(w.category, DeprecationWarning)]
        assert len(deprecations) == 1
        assert "calibration" in str(deprecations[0].message).lower()
    # The v2 calibration_file is remapped to binarization_file;
    # n_density_bins -> n_context_bins.
    assert cfg.binarization.binarization_file == "/tmp/legacy_cal.json"
    assert cfg.binarization.n_context_bins == 6


def test_tooconfig_from_yaml_both_sections_binarization_wins(tmp_path):
    """When a YAML file has BOTH sections, the explicit ``binarization:``
    section wins and the legacy ``calibration:`` section is silently
    dropped (the user has already migrated; no need to remap)."""
    import warnings
    from finaleme_too.config import TOOConfig

    yaml_path = tmp_path / "mixed_config.yaml"
    yaml_path.write_text(
        "binarization:\n"
        "  n_context_bins: 8\n"
        "  max_error_rate: 0.12\n"
        "calibration:\n"
        "  n_density_bins: 999  # ignored because binarization: wins\n"
    )
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        cfg = TOOConfig.from_yaml(yaml_path)
    # A DeprecationWarning still fires because ``calibration:`` is a
    # v2-only key; the migration notice is useful even when the user
    # already has a ``binarization:`` section alongside it.
    deprecations = [w for w in captured if issubclass(w.category, DeprecationWarning)]
    assert len(deprecations) == 1
    # The explicit binarization: section takes precedence for the values
    assert cfg.binarization.n_context_bins == 8
    assert cfg.binarization.max_error_rate == 0.12


def test_cli_train_calibration_hard_breaks_with_exit_2():
    """The v2 ``train-calibration`` CLI command must exit non-zero with a
    migration error message pointing at the new train-binarization command."""
    from click.testing import CliRunner
    from finaleme_too.cli import main

    runner = CliRunner()
    # With mix_stderr=True (default), click.echo(..., err=True) lands in
    # result.output alongside any stdout content.
    result = runner.invoke(main, ["train-calibration"])
    assert result.exit_code == 2
    assert "train-binarization" in result.output
    assert (
        "binarization" in result.output.lower()
        or "v3" in result.output
    )


def test_cli_train_binarization_has_help():
    """The new ``train-binarization`` command exists and documents its flags."""
    from click.testing import CliRunner
    from finaleme_too.cli import main

    runner = CliRunner()
    result = runner.invoke(main, ["train-binarization", "--help"])
    assert result.exit_code == 0
    assert "--matched-wgbs" in result.output
    assert "--matched-finaleme" in result.output
    assert "--n-bins-candidates" in result.output
    assert "--max-error-rate" in result.output
    assert "--cv-method" in result.output
    assert "--output" in result.output


def test_cli_run_has_binarization_flag():
    """The ``run`` command exposes the ``--binarization`` flag (Phase E
    removed the legacy ``--calibration`` flag entirely)."""
    from click.testing import CliRunner
    from finaleme_too.cli import main

    runner = CliRunner()
    result = runner.invoke(main, ["run", "--help"])
    assert result.exit_code == 0
    assert "--binarization" in result.output
    assert "--binarizeThreshold" in result.output
    # The legacy v2 --calibration flag should be gone by Phase E
    assert "--calibration" not in result.output


def test_train_binarization_end_to_end(tmp_path):
    """End-to-end smoke test of train_binarization() on a synthetic
    single-sample matched dataset."""
    from finaleme_too.preprocessing.binarization import (
        BinarizationParams,
        train_binarization,
    )

    # Build a tiny matched dataset: 100 markers across chr1..chr4, with
    # predictably-aligned WGBS and FinaleMe calls.
    rng = np.random.default_rng(42)
    n_chroms = 4
    n_per_chrom = 100
    n_total = n_chroms * n_per_chrom

    # Tile sample_id / chrom / start / end
    rows = []
    for ci in range(n_chroms):
        for i in range(n_per_chrom):
            start = 1000 + i * 200
            # WGBS methylation: half near 0, half near 1 (clean U/M)
            if i < 50:
                wgbs_meth_ct = 1
                fme_meth_ct = 1
            else:
                wgbs_meth_ct = 19
                fme_meth_ct = 19
            rows.append({
                "sample_id": "S1",
                "chrom": f"chr{ci + 1}",
                "start": start,
                "end": start + 100,
                "wgbs_methylated_count": wgbs_meth_ct,
                "wgbs_total_count": 20,
                "fme_methylated_count": fme_meth_ct,
                "fme_total_count": 20,
            })
    df = pd.DataFrame(rows)

    # Write WGBS and FinaleMe matched tables in the legacy joined format
    wgbs_path = tmp_path / "wgbs_matched.tsv"
    wgbs_df = df[["sample_id", "chrom", "start", "end"]].copy()
    wgbs_df["methylated_count"] = df["wgbs_methylated_count"]
    wgbs_df["total_count"] = df["wgbs_total_count"]
    wgbs_df.to_csv(wgbs_path, sep="\t", index=False)

    fme_path = tmp_path / "fme_matched.tsv"
    fme_df = df[["sample_id", "chrom", "start", "end"]].copy()
    fme_df["methylated_count"] = df["fme_methylated_count"]
    fme_df["total_count"] = df["fme_total_count"]
    fme_df.to_csv(fme_path, sep="\t", index=False)

    out_params = tmp_path / "bin_params.json"
    out_report = tmp_path / "bin_report.json"

    params = train_binarization(
        matched_wgbs=wgbs_path,
        matched_finaleme=fme_path,
        region_annotation=None,  # no CpG index → single density bin per class
        n_bins_candidates=[4, 8],  # 1 or 2 sub-bins per class
        out_params=out_params,
        out_report=out_report,
        cv_n_folds=3,
        cv_seed=0,
    )
    assert isinstance(params, BinarizationParams)
    assert out_params.exists()
    assert out_report.exists()

    # Round-trip the saved params to confirm the file is valid JSON
    loaded = BinarizationParams.load(out_params)
    assert loaded.n_bins == params.n_bins
    np.testing.assert_array_equal(loaded.tau_low, params.tau_low)

    # Report contains training metadata
    report = json.loads(out_report.read_text())
    assert report["binarization_version"] == "1.0"
    assert report["n_training_samples"] == 1
    assert report["cv_method"] == "chromosome_blocked"
