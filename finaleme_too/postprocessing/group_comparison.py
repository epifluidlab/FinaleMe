"""Parse --group-comparison syntax and run omnibus + pairwise tests."""

from __future__ import annotations

from itertools import combinations

import numpy as np

from finaleme_too.config import TestMethod
from finaleme_too.postprocessing.statistical_testing import (
    TestResult,
    compositional_regression_test,
    wilcoxon_test,
)
from finaleme_too.utils.transforms import ilr_transform


def parse_group_comparison(
    spec: str | None, available_groups: list[str]
) -> tuple[bool, list[tuple[str, str]]]:
    """Parse --group-comparison spec.

    Supported syntaxes:
        all                          → all pairs of groups, no omnibus
        omnibus+pairwise             → all pairs + omnibus
        A:B,C:D                      → specific contrasts
        A:rest                       → A vs each other group
    """
    if spec is None or spec.strip() == "":
        return False, []

    spec_clean = spec.strip().lower()
    do_omnibus = False

    if spec_clean.startswith("omnibus"):
        do_omnibus = True
        # If just "omnibus", run pairwise too by default
        if spec_clean in ("omnibus", "omnibus+pairwise"):
            return do_omnibus, list(combinations(available_groups, 2))

    if spec_clean == "all":
        return do_omnibus, list(combinations(available_groups, 2))

    # A:rest pattern
    if ":rest" in spec_clean:
        target = spec.split(":")[0]
        return do_omnibus, [(target, g) for g in available_groups if g != target]

    # Comma-separated A:B contrasts
    contrasts: list[tuple[str, str]] = []
    for piece in spec.split(","):
        piece = piece.strip()
        if not piece or ":" not in piece:
            continue
        a, b = piece.split(":", 1)
        contrasts.append((a.strip(), b.strip()))
    return do_omnibus, contrasts


def omnibus_kruskal(
    proportions: np.ndarray,
    group_labels: list[str | None],
    cell_type_names: list[str],
) -> list[TestResult]:
    """Per-cell-type Kruskal-Wallis omnibus test (architecture §9.3)."""
    from scipy.stats import kruskal

    groups = sorted({g for g in group_labels if g is not None})
    results: list[TestResult] = []
    for j, ct in enumerate(cell_type_names):
        samples_per_group = [
            proportions[[i for i, g in enumerate(group_labels) if g == grp], j]
            for grp in groups
        ]
        # Filter out empty groups
        samples_per_group = [s for s in samples_per_group if len(s) >= 2]
        if len(samples_per_group) < 2:
            continue
        try:
            stat, p = kruskal(*samples_per_group)
        except ValueError:
            continue
        results.append(
            TestResult(
                cell_type=ct,
                contrast="all_groups",
                test_type="omnibus",
                mean_a=float("nan"),
                mean_b=float("nan"),
                effect_size=float("nan"),
                se=float("nan"),
                statistic=float(stat),
                p_value=float(p),
            )
        )
    return results


def run_group_comparisons(
    proportions: np.ndarray,
    sample_ids: list[str],
    group_labels: list[str | None],
    cell_type_names: list[str],
    spec: str | None,
    method: TestMethod = TestMethod.ILR_REGRESSION,
    fdr_alpha: float = 0.05,
) -> list[TestResult]:
    """Top-level dispatcher."""
    available = sorted({g for g in group_labels if g is not None})
    if len(available) < 2 or spec is None:
        return []

    do_omnibus, contrasts = parse_group_comparison(spec, available)

    results: list[TestResult] = []
    if do_omnibus:
        results.extend(omnibus_kruskal(proportions, group_labels, cell_type_names))

    if contrasts:
        if method == TestMethod.WILCOXON:
            pairwise = wilcoxon_test(
                proportions, group_labels, cell_type_names, contrasts, fdr_alpha
            )
        else:
            pairwise = compositional_regression_test(
                proportions, sample_ids, group_labels, cell_type_names, contrasts,
                fdr_alpha=fdr_alpha,
            )
        results.extend(pairwise)

    return results


__all__ = [
    "omnibus_kruskal",
    "parse_group_comparison",
    "run_group_comparisons",
]
