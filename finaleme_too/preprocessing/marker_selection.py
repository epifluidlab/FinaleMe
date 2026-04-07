"""Cell-type-balanced marker selection (architecture §5.5).

Picks the top-N most discriminative markers per cell type using one of three
specificity scores: entropy (default), t-statistic, or delta-of-means.
"""

from __future__ import annotations

import numpy as np

from finaleme_too.io.reference_panel import ReferencePanel


def specificity_entropy(reference_methylation: np.ndarray) -> np.ndarray:
    """For each marker, compute Shannon entropy of normalized cell-type levels.

    Lower entropy = more cell-type-specific. Returns (M,) negative entropy
    so that higher = more specific.
    """
    R = np.asarray(reference_methylation, dtype=np.float64)
    eps = 1e-9
    p = (R + eps) / (R.sum(axis=1, keepdims=True) + eps * R.shape[1])
    H = -np.sum(p * np.log(p), axis=1)
    return -H  # negate so larger values = more specific


def specificity_t(reference_methylation: np.ndarray, target: int) -> np.ndarray:
    """Per-marker absolute t-statistic for ``target`` vs all other cell types."""
    R = np.asarray(reference_methylation, dtype=np.float64)
    target_vals = R[:, target]
    others = np.delete(R, target, axis=1)
    o_mean = others.mean(axis=1)
    o_var = others.var(axis=1, ddof=1) if others.shape[1] > 1 else np.full(R.shape[0], 1e-6)
    se = np.sqrt(o_var / max(others.shape[1], 1))
    return np.abs(target_vals - o_mean) / np.maximum(se, 1e-9)


def specificity_delta_means(reference_methylation: np.ndarray, target: int) -> np.ndarray:
    R = np.asarray(reference_methylation, dtype=np.float64)
    target_vals = R[:, target]
    others = np.delete(R, target, axis=1)
    return np.abs(target_vals - others.mean(axis=1))


class BalancedMarkerSelector:
    """Top-N markers per cell type, optionally restricted to a region class."""

    def __init__(
        self,
        n_per_type: int = 500,
        method: str = "entropy",
        strict_regions: str | None = None,
    ):
        self.n_per_type = n_per_type
        self.method = method
        self.strict_regions = strict_regions

    def select(
        self,
        reference: ReferencePanel,
        region_annotations: object | None = None,
    ) -> np.ndarray:
        """Return sorted unique indices of selected markers."""
        n_markers = reference.n_markers
        n_ct = reference.n_cell_types
        eligible = np.ones(n_markers, dtype=bool)

        # Strict region filter (e.g. CGI+shore)
        if self.strict_regions and region_annotations is not None:
            eligible &= self._filter_by_region_class(
                reference, region_annotations, self.strict_regions
            )

        chosen: set[int] = set()
        R = reference.methylation
        if self.method == "entropy":
            ent = specificity_entropy(R)
            for j in range(n_ct):
                # Combine global entropy with per-cell-type delta to break ties
                delta = specificity_delta_means(R, j)
                score = ent + delta
                self._take_top(score, eligible, chosen)
        else:
            scorer = specificity_t if self.method == "t_statistic" else specificity_delta_means
            for j in range(n_ct):
                self._take_top(scorer(R, j), eligible, chosen)

        if not chosen:
            return np.arange(n_markers)
        return np.array(sorted(chosen), dtype=np.int64)

    def _take_top(
        self, score: np.ndarray, eligible: np.ndarray, chosen: set[int]
    ) -> None:
        masked = np.where(eligible, score, -np.inf)
        n_to_take = min(self.n_per_type, int(np.sum(eligible)))
        if n_to_take <= 0:
            return
        idx = np.argpartition(-masked, n_to_take - 1)[:n_to_take]
        for i in idx:
            chosen.add(int(i))

    @staticmethod
    def _filter_by_region_class(
        reference: ReferencePanel,
        region_annotations,
        strict: str,
    ) -> np.ndarray:
        """Match marker (chrom, start, end) against region_annotations[region_class]."""
        import pandas as pd

        wanted = {x.strip().lower() for x in strict.split("+")}
        try:
            ann = pd.DataFrame(region_annotations)
        except Exception:
            return np.ones(reference.n_markers, dtype=bool)
        if "region_class" not in ann.columns:
            return np.ones(reference.n_markers, dtype=bool)
        key = list(zip(ann["chrom"], ann["start"], ann["end"], ann["region_class"]))
        lookup = {(c, s, e): rc.lower() for c, s, e, rc in key}
        mask = np.zeros(reference.n_markers, dtype=bool)
        for i in range(reference.n_markers):
            k = (str(reference.chrom[i]), int(reference.start[i]), int(reference.end[i]))
            cls = lookup.get(k)
            if cls is not None and cls in wanted:
                mask[i] = True
        return mask


__all__ = [
    "BalancedMarkerSelector",
    "specificity_delta_means",
    "specificity_entropy",
    "specificity_t",
]
