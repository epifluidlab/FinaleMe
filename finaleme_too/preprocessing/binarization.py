"""FinaleMe context-dependent binarization — apply + train (v3).

Replaces the v2 continuous calibration model. See:

  * ``design/TOO_ARCHITECTURE_v3.md`` §5.3 (ContextBinarizer)
  * ``design/TOO_MATH_FORMULATION_v3.md`` §2B (binarization likelihood) and §6

The model learns per-context-bin thresholds ``(τ_low_b, τ_high_b)`` and
error rates ``(ε_U_b, ε_M_b)`` for classifying FinaleMe predicted methylation
into ``{U, M, Ambiguous}`` states. Context bins are a 2-level grid over:

  1. Region class (CGI / shore / shelf / open_sea, derived from CpG density
     via ``_matched_sample_sheet._classify_region``)
  2. CpG density sub-bin within each class (quantile split)

The total bin count ``n_bins = n_region_classes × density_sub_bins_per_class``.
With the default 4 classes and ``density_sub_bins_per_class=2`` this gives
the architecture-doc default ``B=8``.

Calling state values (matches ``MarkerObservations.called_state``):

  0 = U          (called unmethylated)
  1 = M          (called methylated)
  2 = Ambiguous  (τ_low < predicted_beta < τ_high)
  3 = Excluded   (bin not usable, or no annotation)
"""

from __future__ import annotations

import json
import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from finaleme_too.exceptions import InvalidBinarizationError
from finaleme_too.io.methylation_loader import MarkerObservations
from finaleme_too.preprocessing._matched_sample_sheet import (
    DEFAULT_REGION_CLASS_THRESHOLDS,
    REGION_CLASS_ORDER,
    _classify_region,
)

log = logging.getLogger(__name__)

# Called-state integer codes. Kept at module scope so ``binarization_eval``
# and ``observation_model_binarization`` can import and reuse them.
STATE_U = 0
STATE_M = 1
STATE_AMBIGUOUS = 2
STATE_EXCLUDED = 3

STATE_LABELS = {
    STATE_U: "U",
    STATE_M: "M",
    STATE_AMBIGUOUS: "Ambiguous",
    STATE_EXCLUDED: "Excluded",
}


# ---------------------------------------------------------------------------
# BinarizationParams dataclass (JSON-serializable; shipped as default_binarization.json)
# ---------------------------------------------------------------------------


@dataclass
class BinarizationParams:
    """Per-context-bin binarization parameters (math doc §6, arch §5.3).

    Layout (v3.0): ``n_bins = n_region_classes × density_sub_bins_per_class``.
    Bin index is computed as:

        bin_idx = class_idx * density_sub_bins_per_class + density_sub_bin

    where ``class_idx`` is the index of the marker's region_class in
    ``REGION_CLASS_ORDER`` and ``density_sub_bin`` is found by
    ``np.searchsorted`` into the class's per-class density edges.

    Attributes
    ----------
    n_bins
        Total number of context bins (must equal
        ``n_region_classes * density_sub_bins_per_class``).
    n_region_classes
        Number of region classes. v3.0 uses 4 (CGI, shore, shelf, open_sea).
    density_sub_bins_per_class
        Number of density quantile sub-bins within each region class.
    region_class_order
        Ordered list of region class names. Index in this list determines
        ``class_idx`` for bin assignment.
    region_class_thresholds
        CpG density cutoffs used by ``_classify_region`` to map markers to
        region classes. Persisted so inference uses the same thresholds as
        training.
    density_edges
        Shape ``(n_region_classes, density_sub_bins_per_class + 1)``.
        Per-class density quantile edges. ``density_edges[c, 0] = -inf``,
        ``density_edges[c, -1] = +inf`` to catch out-of-range values.
    tau_low, tau_high
        Shape ``(n_bins,)``. Per-bin thresholds for calling U (< tau_low)
        and M (> tau_high). Markers with predicted_beta in ``[tau_low,
        tau_high]`` are Ambiguous.
    eps_U, eps_M
        Shape ``(n_bins,)``. Per-bin misclassification rates:
        ``eps_U = P(true M | called U)``, ``eps_M = P(true U | called M)``.
    usable
        Shape ``(n_bins,)``. Boolean mask. Markers in unusable bins are
        excluded from the likelihood (state = Excluded).
    n_markers_U, n_markers_M
        Shape ``(n_bins,)``. Training-time marker counts per bin per state,
        used to derive error rates + usability.
    train_fraction_U, train_fraction_M
        Shape ``(n_bins,)``. Training-time fraction of called markers in
        each state. Used by ``compute_inference_qc`` to detect state
        distribution shift (KL divergence vs. training U:M ratio).
    training_metadata
        Arbitrary dict of training provenance (n_training_samples,
        cv_accuracy, candidate scores, etc.).
    """

    n_bins: int
    n_region_classes: int
    density_sub_bins_per_class: int
    region_class_order: list[str]
    region_class_thresholds: dict[str, float]
    density_edges: np.ndarray  # shape (n_region_classes, density_sub_bins_per_class + 1)
    tau_low: np.ndarray  # shape (n_bins,)
    tau_high: np.ndarray  # shape (n_bins,)
    eps_U: np.ndarray  # shape (n_bins,)
    eps_M: np.ndarray  # shape (n_bins,)
    usable: np.ndarray  # shape (n_bins,) bool
    n_markers_U: np.ndarray  # shape (n_bins,) int
    n_markers_M: np.ndarray  # shape (n_bins,) int
    train_fraction_U: np.ndarray  # shape (n_bins,) float
    train_fraction_M: np.ndarray  # shape (n_bins,) float
    training_metadata: dict = field(default_factory=dict)

    # ------------------------------------------------------------------
    # Round-trip serialization
    # ------------------------------------------------------------------

    @classmethod
    def from_dict(cls, raw: dict) -> "BinarizationParams":
        try:
            return cls(
                n_bins=int(raw["n_bins"]),
                n_region_classes=int(raw["n_region_classes"]),
                density_sub_bins_per_class=int(raw["density_sub_bins_per_class"]),
                region_class_order=list(raw["region_class_order"]),
                region_class_thresholds=dict(raw["region_class_thresholds"]),
                density_edges=np.asarray(raw["density_edges"], dtype=np.float64),
                tau_low=np.asarray(raw["tau_low"], dtype=np.float64),
                tau_high=np.asarray(raw["tau_high"], dtype=np.float64),
                eps_U=np.asarray(raw["eps_U"], dtype=np.float64),
                eps_M=np.asarray(raw["eps_M"], dtype=np.float64),
                usable=np.asarray(raw["usable"], dtype=bool),
                n_markers_U=np.asarray(raw["n_markers_U"], dtype=np.int64),
                n_markers_M=np.asarray(raw["n_markers_M"], dtype=np.int64),
                train_fraction_U=np.asarray(raw["train_fraction_U"], dtype=np.float64),
                train_fraction_M=np.asarray(raw["train_fraction_M"], dtype=np.float64),
                training_metadata=dict(raw.get("training_metadata", {})),
            )
        except KeyError as exc:
            raise InvalidBinarizationError(
                f"Binarization JSON missing required key: {exc}"
            ) from exc

    @classmethod
    def load(cls, path: str | Path) -> "BinarizationParams":
        p = Path(path)
        if not p.exists():
            raise InvalidBinarizationError(f"Binarization file not found: {p}")
        with open(p) as fh:
            raw = json.load(fh)
        return cls.from_dict(raw)

    def save(self, path: str | Path) -> None:
        out = {
            "n_bins": self.n_bins,
            "n_region_classes": self.n_region_classes,
            "density_sub_bins_per_class": self.density_sub_bins_per_class,
            "region_class_order": list(self.region_class_order),
            "region_class_thresholds": dict(self.region_class_thresholds),
            "density_edges": self.density_edges.tolist(),
            "tau_low": self.tau_low.tolist(),
            "tau_high": self.tau_high.tolist(),
            "eps_U": self.eps_U.tolist(),
            "eps_M": self.eps_M.tolist(),
            "usable": self.usable.astype(bool).tolist(),
            "n_markers_U": self.n_markers_U.astype(np.int64).tolist(),
            "n_markers_M": self.n_markers_M.astype(np.int64).tolist(),
            "train_fraction_U": self.train_fraction_U.tolist(),
            "train_fraction_M": self.train_fraction_M.tolist(),
            "training_metadata": _to_jsonable(self.training_metadata),
        }
        Path(path).write_text(json.dumps(out, indent=2))

    # ------------------------------------------------------------------
    # Bin assignment
    # ------------------------------------------------------------------

    def assign_bin(
        self,
        cpg_density: np.ndarray,
        region_class: np.ndarray | None = None,
    ) -> np.ndarray:
        """Map per-marker (cpg_density, region_class) to a bin index.

        If ``region_class`` is not provided, it is derived from
        ``cpg_density`` using ``self.region_class_thresholds``. This is the
        automatic path that works even when the caller doesn't have a
        precomputed ``region_class`` column.

        NaN densities are routed deterministically to the ``open_sea``
        class and the lowest density sub-bin within that class, matching
        the v2 calibration bug-4 fix (no NaN propagation).
        """
        density_arr = np.asarray(cpg_density, dtype=np.float64)
        n = density_arr.size

        if region_class is None:
            region_class = _classify_region(density_arr, self.region_class_thresholds)
        else:
            region_class = np.asarray(region_class, dtype=object)

        # Map class name → class index
        class_to_idx = {c: i for i, c in enumerate(self.region_class_order)}
        # Unknown class names (e.g. typos, or a brand-new class not in the
        # training set) fall back to open_sea.
        fallback_class_idx = class_to_idx.get("open_sea", 0)
        class_idx = np.full(n, fallback_class_idx, dtype=np.int64)
        for c, idx in class_to_idx.items():
            class_idx[region_class == c] = idx

        # Per-class density sub-bin lookup
        nan_mask = ~np.isfinite(density_arr)
        clean_density = np.where(nan_mask, -np.inf, density_arr)

        sub_bin = np.zeros(n, dtype=np.int64)
        for c_idx in range(self.n_region_classes):
            mask = class_idx == c_idx
            if not np.any(mask):
                continue
            edges = self.density_edges[c_idx]
            sub_idx = np.clip(
                np.searchsorted(edges, clean_density[mask], side="right") - 1,
                0,
                self.density_sub_bins_per_class - 1,
            )
            sub_bin[mask] = sub_idx

        # NaN density → sub-bin 0 (deterministic fallback)
        sub_bin[nan_mask] = 0

        return (class_idx * self.density_sub_bins_per_class + sub_bin).astype(np.int64)


# ---------------------------------------------------------------------------
# Apply path: classify each marker into U / M / Ambiguous / Excluded
# ---------------------------------------------------------------------------


def apply_binarization(
    obs: MarkerObservations,
    params: BinarizationParams,
    region_annotations: pd.DataFrame | None = None,
) -> MarkerObservations:
    """Apply a trained binarization model to one sample's observations.

    Inputs
    ------
    obs
        FinaleMe-mode observations with ``predicted_beta`` populated.
    params
        Pre-trained binarization parameters.
    region_annotations
        Optional DataFrame with columns ``chrom, start, end, cpg_density``
        (and optionally ``region_class``). When provided, the per-marker
        density/region_class is looked up by joining on (chrom, start, end).
        When ``None`` or when a marker has no matching row, the marker
        falls back to ``cpg_density=NaN`` which routes it to the
        open_sea class + lowest sub-bin (see ``assign_bin``).

    Returns
    -------
    A new ``MarkerObservations`` with ``called_state`` and ``context_bin``
    populated. ``k``, ``n``, ``predicted_beta`` are unchanged.
    """
    if obs.predicted_beta is None:
        # No predictions to binarize — leave obs unchanged.
        return obs

    pred = np.asarray(obs.predicted_beta, dtype=np.float64)
    n = pred.size

    # Look up per-marker density + region_class from the annotation table.
    if region_annotations is not None and not region_annotations.empty:
        keys = list(zip(obs.chrom.tolist(), obs.start.tolist(), obs.end.tolist()))
        ann_indexed = region_annotations.set_index(["chrom", "start", "end"])
        density_lookup = ann_indexed["cpg_density"]
        density = np.array(
            [float(density_lookup.get(k, np.nan)) for k in keys], dtype=np.float64
        )
        if "region_class" in ann_indexed.columns:
            class_lookup = ann_indexed["region_class"]
            region_class = np.array(
                [
                    class_lookup.get(k, "open_sea") if k in class_lookup.index else "open_sea"
                    for k in keys
                ],
                dtype=object,
            )
        else:
            region_class = _classify_region(density, params.region_class_thresholds)
    else:
        density = np.full(n, np.nan, dtype=np.float64)
        region_class = _classify_region(density, params.region_class_thresholds)

    bin_idx = params.assign_bin(density, region_class)

    # Classify each marker using its bin's thresholds
    called_state = np.full(n, STATE_EXCLUDED, dtype=np.uint8)
    # Vectorized comparison: fetch per-marker τ_low, τ_high, usable
    tau_low_per = params.tau_low[bin_idx]
    tau_high_per = params.tau_high[bin_idx]
    usable_per = params.usable[bin_idx]

    # Start: everything defaults to Excluded
    # Usable bins: classify by predicted_beta
    finite = np.isfinite(pred)
    # U call
    u_mask = usable_per & finite & (pred < tau_low_per)
    called_state[u_mask] = STATE_U
    # M call
    m_mask = usable_per & finite & (pred > tau_high_per)
    called_state[m_mask] = STATE_M
    # Ambiguous: usable + finite + in [τ_low, τ_high]
    amb_mask = usable_per & finite & ~u_mask & ~m_mask
    called_state[amb_mask] = STATE_AMBIGUOUS
    # Markers in unusable bins or with NaN predicted_beta stay EXCLUDED.

    return obs.with_binarization(called_state=called_state, context_bin=bin_idx)


# ---------------------------------------------------------------------------
# Shipped default
# ---------------------------------------------------------------------------


def shipped_default_binarization_path() -> Path:
    """Path to the default binarization JSON shipped with the package."""
    from finaleme_too import data as data_pkg

    pkg_dir = Path(data_pkg.__file__).parent
    return pkg_dir / "default_binarization.json"


def load_default_binarization() -> BinarizationParams:
    return BinarizationParams.load(shipped_default_binarization_path())


def build_identity_placeholder_params(
    region_class_thresholds: dict[str, float] | None = None,
    density_sub_bins_per_class: int = 2,
    tau_low: float = 0.2,
    tau_high: float = 0.8,
    eps: float = 0.1,
) -> BinarizationParams:
    """Construct a placeholder ``BinarizationParams`` with no training.

    Used by Phase A to ship ``data/default_binarization.json`` before a
    real trained default is available. Every bin has the same τ_low, τ_high,
    ε_U, ε_M values and is marked usable. The density_edges are set to
    ``[-inf, 0, +inf]`` so every marker within a class lands in sub-bin 0
    (the sub-bin structure is preserved but effectively flat).

    Matches the Phase A plan: 8 bins (4 classes × 2 density sub-bins each),
    τ_low=0.2, τ_high=0.8, ε=0.1, usable=true, default density thresholds.
    """
    if region_class_thresholds is None:
        region_class_thresholds = dict(DEFAULT_REGION_CLASS_THRESHOLDS)
    classes = list(REGION_CLASS_ORDER)
    n_classes = len(classes)
    n_bins = n_classes * density_sub_bins_per_class

    # Per-class density edges. For the placeholder, interior edges are set
    # well above the realistic methylation density range so every realistic
    # marker lands in sub-bin 0. Real training replaces these with
    # class-specific quantile edges learned from data.
    density_edges = np.zeros(
        (n_classes, density_sub_bins_per_class + 1), dtype=np.float64
    )
    density_edges[:, 0] = -np.inf
    density_edges[:, -1] = np.inf
    # Interior edges are placed at 1.0 — far above any realistic CpG density
    # (which is typically < 0.2), so np.searchsorted assigns every marker to
    # sub-bin 0 within its region class.
    for j in range(1, density_sub_bins_per_class):
        density_edges[:, j] = 1.0

    return BinarizationParams(
        n_bins=n_bins,
        n_region_classes=n_classes,
        density_sub_bins_per_class=density_sub_bins_per_class,
        region_class_order=classes,
        region_class_thresholds=region_class_thresholds,
        density_edges=density_edges,
        tau_low=np.full(n_bins, tau_low, dtype=np.float64),
        tau_high=np.full(n_bins, tau_high, dtype=np.float64),
        eps_U=np.full(n_bins, eps, dtype=np.float64),
        eps_M=np.full(n_bins, eps, dtype=np.float64),
        usable=np.ones(n_bins, dtype=bool),
        n_markers_U=np.zeros(n_bins, dtype=np.int64),
        n_markers_M=np.zeros(n_bins, dtype=np.int64),
        train_fraction_U=np.full(n_bins, 0.5, dtype=np.float64),
        train_fraction_M=np.full(n_bins, 0.5, dtype=np.float64),
        training_metadata={
            "kind": "identity_placeholder",
            "note": (
                "Placeholder default shipped in v3 Phase A. Replace with a "
                "real trained binarization via `finaleme-too train-binarization`."
            ),
        },
    )


# ---------------------------------------------------------------------------
# Internals
# ---------------------------------------------------------------------------


def _to_jsonable(obj: Any) -> Any:
    """Recursively convert numpy types in the training_metadata dict to
    JSON-serializable Python natives."""
    if isinstance(obj, dict):
        return {k: _to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    return obj


__all__ = [
    "STATE_U",
    "STATE_M",
    "STATE_AMBIGUOUS",
    "STATE_EXCLUDED",
    "STATE_LABELS",
    "BinarizationParams",
    "apply_binarization",
    "build_identity_placeholder_params",
    "load_default_binarization",
    "shipped_default_binarization_path",
]
