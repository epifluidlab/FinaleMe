"""Tests for ILR transform and compositional regression (Phase B)."""

from __future__ import annotations

import numpy as np

from finaleme_too.config import TestMethod
from finaleme_too.postprocessing.group_comparison import (
    parse_group_comparison,
    run_group_comparisons,
)
from finaleme_too.postprocessing.statistical_testing import (
    compositional_regression_test,
    wilcoxon_test,
)
from finaleme_too.utils.transforms import (
    clr_inverse,
    clr_transform,
    helmert_basis,
    ilr_inverse,
    ilr_transform,
)


def test_clr_round_trip():
    rng = np.random.default_rng(0)
    p = rng.dirichlet(np.ones(5), size=10)
    y = clr_transform(p)
    p2 = clr_inverse(y)
    np.testing.assert_allclose(p, p2, atol=1e-9)


def test_ilr_round_trip():
    rng = np.random.default_rng(1)
    p = rng.dirichlet(np.ones(5), size=10)
    y = ilr_transform(p)
    assert y.shape == (10, 4)
    p2 = ilr_inverse(y, n_parts=5)
    np.testing.assert_allclose(p, p2, atol=1e-9)


def test_helmert_basis_orthonormal():
    V = helmert_basis(5)
    # V.T @ V should be identity in (d-1)-dim
    G = V.T @ V
    np.testing.assert_allclose(G, np.eye(4), atol=1e-12)


def test_compositional_regression_detects_clear_difference():
    """Two groups with different proportions in cell type 1 should be detected."""
    rng = np.random.default_rng(42)
    n_per_group = 8
    K_total = 4  # 3 cell types + unknown

    group_a = []
    for _ in range(n_per_group):
        # Group A: cell type 0 around 0.6
        group_a.append(rng.dirichlet([0.6, 0.15, 0.15, 0.1] * np.array([10, 10, 10, 10])))
    group_b = []
    for _ in range(n_per_group):
        # Group B: cell type 0 around 0.2
        group_b.append(rng.dirichlet([0.2, 0.35, 0.35, 0.1] * np.array([10, 10, 10, 10])))

    proportions = np.vstack(group_a + group_b)
    sample_ids = [f"a{i}" for i in range(n_per_group)] + [f"b{i}" for i in range(n_per_group)]
    labels = ["A"] * n_per_group + ["B"] * n_per_group
    cell_types = ["CT1", "CT2", "CT3"]  # without unknown

    results = compositional_regression_test(
        proportions=proportions,
        sample_ids=sample_ids,
        group_labels=labels,
        cell_type_names=cell_types,
        contrasts=[("A", "B")],
        fdr_alpha=0.05,
    )
    assert len(results) == len(cell_types)
    assert all(r.p_value <= 1.0 and r.p_value >= 0.0 for r in results)
    # CT1 should be the most significant (the cell type that differs)
    ct1_p = next(r.adjusted_p_value for r in results if r.cell_type == "CT1")
    assert ct1_p < 0.1


def test_wilcoxon_matches_scipy():
    from scipy.stats import mannwhitneyu

    proportions = np.array([
        [0.1, 0.9],
        [0.15, 0.85],
        [0.2, 0.8],
        [0.5, 0.5],
        [0.55, 0.45],
        [0.6, 0.4],
    ])
    labels = ["A", "A", "A", "B", "B", "B"]
    results = wilcoxon_test(
        proportions, labels, ["CT1", "CT2"], [("A", "B")], fdr_alpha=0.05
    )
    expected_p = mannwhitneyu(
        proportions[:3, 0], proportions[3:, 0], alternative="two-sided"
    ).pvalue
    actual_p = next(r.p_value for r in results if r.cell_type == "CT1")
    assert abs(actual_p - expected_p) < 1e-9


def test_parse_group_comparison_all():
    do_omnibus, contrasts = parse_group_comparison("all", ["A", "B", "C"])
    assert not do_omnibus
    assert set(contrasts) == {("A", "B"), ("A", "C"), ("B", "C")}


def test_parse_group_comparison_omnibus_pairwise():
    do_omnibus, contrasts = parse_group_comparison("omnibus+pairwise", ["A", "B"])
    assert do_omnibus
    assert contrasts == [("A", "B")]


def test_parse_group_comparison_specific():
    _, contrasts = parse_group_comparison("X:Y,Z:W", ["X", "Y", "Z", "W"])
    assert contrasts == [("X", "Y"), ("Z", "W")]


def test_parse_group_comparison_one_vs_rest():
    _, contrasts = parse_group_comparison("X:rest", ["X", "Y", "Z"])
    assert contrasts == [("X", "Y"), ("X", "Z")]


def test_run_group_comparisons_dispatches_method():
    rng = np.random.default_rng(7)
    P = rng.dirichlet(np.ones(3), size=10)
    labels = ["A"] * 5 + ["B"] * 5
    out = run_group_comparisons(
        proportions=P,
        sample_ids=[f"s{i}" for i in range(10)],
        group_labels=labels,
        cell_type_names=["CT1", "CT2"],
        spec="A:B",
        method=TestMethod.WILCOXON,
    )
    assert all(r.test_type == "wilcoxon" for r in out)
