"""Tests for batch correction, imputation, and covariate adjustment (Phase D)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from finaleme_too.config import MeasurementMode
from finaleme_too.exceptions import IllegalImputationError
from finaleme_too.io.methylation_loader import MarkerObservations
from finaleme_too.postprocessing.covariate_adjustment import adjust_covariates
from finaleme_too.preprocessing.batch_correction import combat_correct
from finaleme_too.preprocessing.imputation import CohortImputer


def _mk_obs(sample_id: str, k: list[int], n: list[int]) -> MarkerObservations:
    M = len(k)
    return MarkerObservations(
        sample_id=sample_id,
        chrom=np.array(["chr1"] * M, dtype=object),
        start=np.array([100 + i * 1000 for i in range(M)], dtype=np.int64),
        end=np.array([200 + i * 1000 for i in range(M)], dtype=np.int64),
        k=np.array(k, dtype=np.int32),
        n=np.array(n, dtype=np.int32),
        predicted_beta=None,
        mode=MeasurementMode.WGBS,
    )


def test_combat_removes_synthetic_batch_shift():
    rng = np.random.default_rng(0)
    n_markers = 50
    truth = rng.uniform(0.2, 0.8, size=n_markers)
    samples = []
    labels = []
    n_per_batch = 10
    for s in range(n_per_batch):
        # Batch A: methylation shifted by +0.1
        beta = np.clip(truth + 0.1 + rng.normal(0, 0.02, size=n_markers), 0, 1)
        k = (beta * 100).astype(int).tolist()
        samples.append(_mk_obs(f"a_{s}", k, [100] * n_markers))
        labels.append("A")
    for s in range(n_per_batch):
        beta = np.clip(truth - 0.1 + rng.normal(0, 0.02, size=n_markers), 0, 1)
        k = (beta * 100).astype(int).tolist()
        samples.append(_mk_obs(f"b_{s}", k, [100] * n_markers))
        labels.append("B")

    # Pre-correction: per-batch means differ
    pre_a = np.mean([np.mean(o.k / o.n) for o in samples[:n_per_batch]])
    pre_b = np.mean([np.mean(o.k / o.n) for o in samples[n_per_batch:]])
    assert abs(pre_a - pre_b) > 0.1

    corrected = combat_correct(samples, labels, min_levels=2, min_per_level=5)
    post_a = np.mean([np.mean(o.k / o.n) for o in corrected[:n_per_batch]])
    post_b = np.mean([np.mean(o.k / o.n) for o in corrected[n_per_batch:]])
    # Post-correction: per-batch means should be much closer
    assert abs(post_a - post_b) < 0.05


def test_combat_skips_when_too_few_levels():
    samples = [_mk_obs(f"s{i}", [50] * 5, [100] * 5) for i in range(5)]
    out = combat_correct(samples, batch_labels=["A"] * 5, min_levels=2)
    # Identical: should be returned unchanged
    for a, b in zip(samples, out):
        np.testing.assert_array_equal(a.k, b.k)


def test_imputer_refuses_when_no_group():
    """In strict mode the imputer must raise when the sample has no group."""
    sample = _mk_obs("s1", [1] * 5, [1] * 5)
    cohort = [_mk_obs(f"c{i}", [10] * 5, [20] * 5) for i in range(3)]
    with pytest.raises(IllegalImputationError):
        CohortImputer().impute(sample, cohort, sample_groups={}, strict=True)


def test_imputer_lenient_returns_sample_when_no_group():
    """In lenient mode (pipeline default) the imputer returns the input
    unchanged instead of raising, so LOW/ULTRALOW runs do not crash on
    unlabeled samples."""
    sample = _mk_obs("s1", [1] * 5, [1] * 5)
    cohort = [_mk_obs(f"c{i}", [10] * 5, [20] * 5) for i in range(3)]
    out = CohortImputer().impute(sample, cohort, sample_groups={}, strict=False)
    np.testing.assert_array_equal(out.k, sample.k)
    np.testing.assert_array_equal(out.n, sample.n)


def test_imputer_only_uses_same_group_donors():
    target = _mk_obs("t1", [0] * 5, [0] * 5)  # zero coverage everywhere
    same = [_mk_obs(f"s{i}", [10] * 5, [20] * 5) for i in range(3)]
    other = [_mk_obs(f"o{i}", [18] * 5, [20] * 5) for i in range(3)]
    cohort = same + other + [target]
    groups = {**{f"s{i}": "A" for i in range(3)}, **{f"o{i}": "B" for i in range(3)}, "t1": "A"}
    imputed = CohortImputer(coverage_threshold=1, min_donors=3).impute(target, cohort, groups)
    # All beta values should be ~0.5 (10/20 from same-group donors), not ~0.9
    beta = imputed.k / np.maximum(imputed.n, 1)
    np.testing.assert_allclose(beta, np.full(5, 0.5), atol=0.1)


def test_imputer_skips_when_too_few_donors():
    target = _mk_obs("t1", [0] * 5, [0] * 5)
    cohort = [_mk_obs("d1", [10] * 5, [20] * 5)]
    groups = {"t1": "A", "d1": "A"}
    out = CohortImputer(min_donors=3).impute(target, cohort, groups)
    np.testing.assert_array_equal(out.k, target.k)
    np.testing.assert_array_equal(out.n, target.n)


def test_adjust_covariates_no_op_with_no_columns():
    rng = np.random.default_rng(0)
    p = rng.dirichlet(np.ones(4), size=5)
    out = adjust_covariates(
        proportions=p,
        sample_ids=[f"s{i}" for i in range(5)],
        covariates=pd.DataFrame({"sample_id": [f"s{i}" for i in range(5)], "age": [30, 40, 50, 60, 70]}),
        columns=[],
    )
    np.testing.assert_array_equal(p, out)


def test_adjust_covariates_residualizes_age():
    rng = np.random.default_rng(1)
    n = 12
    age = np.linspace(30, 70, n)
    # Inject linear age dependence on cell type 0 in ILR space
    base = rng.dirichlet(np.ones(4), size=n)
    base[:, 0] += 0.05 * (age - age.mean())
    base = np.clip(base, 1e-3, None)
    base /= base.sum(axis=1, keepdims=True)

    cov = pd.DataFrame({"sample_id": [f"s{i}" for i in range(n)], "age": age})
    adjusted = adjust_covariates(
        proportions=base,
        sample_ids=[f"s{i}" for i in range(n)],
        covariates=cov,
        columns=["age"],
    )
    # Adjusted proportions should be ~uniform across age (residualized)
    # Correlation between age and cell type 0 should be much smaller after adjustment.
    pre_corr = np.corrcoef(age, base[:, 0])[0, 1]
    post_corr = np.corrcoef(age, adjusted[:, 0])[0, 1]
    assert abs(post_corr) < abs(pre_corr)
