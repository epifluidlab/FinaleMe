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
    # Defense-in-depth chromosome normalization: strip any leading ``chr``
    # from BOTH the obs chroms and the annotation chroms before joining.
    # The annotation writer (``compute_region_annotation``) strips the
    # prefix, but a hand-written or externally-sourced annotation file may
    # still carry it, and conversely the loaded obs may or may not.
    # Without this the join silently misses and every marker defaults to
    # the open_sea fallback bin — the exact issue reported by the user.
    if region_annotations is not None and not region_annotations.empty:
        obs_chrom_norm = [
            c[3:] if isinstance(c, str) and c.startswith("chr") else str(c)
            for c in obs.chrom.tolist()
        ]
        keys = list(zip(obs_chrom_norm, obs.start.tolist(), obs.end.tolist()))

        ann = region_annotations
        if "chrom" in ann.columns:
            ann = ann.copy()
            ann["chrom"] = (
                ann["chrom"].astype(str).str.replace(r"^chr", "", regex=True)
            )
        ann_indexed = ann.set_index(["chrom", "start", "end"])
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
# Phase D: training pipeline orchestrator
# ---------------------------------------------------------------------------


def train_binarization(
    matched_wgbs: str | Path,
    matched_finaleme: str | Path,
    region_annotation: str | Path | None,
    n_bins_candidates: list[int],
    out_params: str | Path,
    out_report: str | Path,
    cpg_index: str | Path | None = None,
    region_annotation_window: int = 1000,
    threads: int = 1,
    max_error_rate: float = 0.15,
    cv_method: str = "chromosome_blocked",
    cv_n_folds: int = 10,
    cv_seed: int | None = None,
) -> BinarizationParams:
    """Train per-bin binarization parameters from matched WGBS / FinaleMe samples.

    v3 replacement for ``train_calibration``. The matched-data loading,
    chromosome normalization, and CpG-density computation are reused from
    ``_matched_sample_sheet`` (the shared parser module). The new
    ``binarization_eval.tune_n_bins`` performs chromosome-blocked K-fold
    CV to pick the best total bin count, and ``fit_binarization`` produces
    the final per-bin (τ_low, τ_high, ε_U, ε_M, usable) parameters.

    Inputs
    ------
    matched_wgbs, matched_finaleme
        Paths to either legacy joined TSVs or sample-sheet TSVs (see
        ``_load_matched_table`` for the auto-detected formats).
    region_annotation
        Optional pre-computed CpG-density + region_class TSV. If omitted
        and ``cpg_index`` is provided, density + region_class are
        auto-computed via ``compute_region_annotation`` over a
        ``region_annotation_window``-bp window. If neither is provided,
        all rows go to the open_sea / sub-bin 0 fallback (degenerate).
    n_bins_candidates
        List of total bin counts to evaluate. v3.0 rounds each value up
        to the nearest multiple of n_region_classes (4) so the per-class
        sub-bin count is uniform. Default candidates: ``[4, 8, 12, 16]``.
    cv_method
        Currently always ``"chromosome_blocked"``. Kept as a parameter
        so future strategies can be added without changing the API.
    cv_n_folds, cv_seed
        Number of chromosome-blocked CV folds and the seed for the
        chromosome shuffle. ``cv_seed`` makes the split reproducible.

    Outputs
    -------
    Returns the trained ``BinarizationParams`` and writes them to
    ``out_params``. Also writes a JSON training report to ``out_report``
    with per-bin metrics, the chosen B, and CV summary.
    """
    from finaleme_too.io.output_writer import write_calibration_report
    from finaleme_too.preprocessing._matched_sample_sheet import (
        _load_matched_table,
        _normalize_chrom,
        compute_region_annotation,
    )
    from finaleme_too.preprocessing.binarization_eval import (
        fit_binarization,
        tune_n_bins,
    )

    threads = max(1, int(threads))
    wgbs_df = _load_matched_table(matched_wgbs, modality="wgbs", threads=threads)
    fme_df = _load_matched_table(matched_finaleme, modality="finaleme", threads=threads)

    join_keys = ["sample_id", "chrom", "start", "end"]
    merged = wgbs_df.merge(fme_df, on=join_keys, suffixes=("_wgbs", "_fme"))
    if merged.empty:
        raise InvalidBinarizationError(
            "No overlapping (sample_id, chrom, start, end) between WGBS and "
            "FinaleMe tables. Check sample IDs, genomic coordinates, and "
            "chromosome naming — the loader already strips any 'chr' prefix."
        )

    # Resolve per-row CpG density + region_class
    if region_annotation is not None:
        ann = pd.read_csv(region_annotation, sep="\t", comment="#")
        ann = _normalize_chrom(ann)
        merged = merged.merge(
            ann[
                [c for c in ("chrom", "start", "end", "cpg_density", "region_class") if c in ann.columns]
            ],
            on=["chrom", "start", "end"],
            how="left",
        )
    elif cpg_index is not None:
        unique_regions = merged[["chrom", "start", "end"]].drop_duplicates()
        ann = compute_region_annotation(
            unique_regions,
            cpg_index_path=cpg_index,
            window=region_annotation_window,
        )
        merged = merged.merge(
            ann[["chrom", "start", "end", "cpg_density", "region_class"]],
            on=["chrom", "start", "end"],
            how="left",
        )
    if "cpg_density" not in merged.columns:
        merged["cpg_density"] = 0.0
    merged["cpg_density"] = merged["cpg_density"].fillna(0.0)

    # Drop zero-coverage rows on either side (matches the v2 fix from bug 2)
    valid_rows = (
        merged["total_count_wgbs"].fillna(0).to_numpy() > 0
    ) & (
        merged["total_count_fme"].fillna(0).to_numpy() > 0
    ) & merged["methylated_count_wgbs"].notna().to_numpy() & merged[
        "methylated_count_fme"
    ].notna().to_numpy()
    n_dropped = int(valid_rows.size - valid_rows.sum())
    if n_dropped > 0:
        log.info(
            "train_binarization: dropped %d / %d rows with zero coverage on either side",
            n_dropped,
            int(valid_rows.size),
        )
    merged = merged.loc[valid_rows].reset_index(drop=True)
    if merged.empty:
        raise InvalidBinarizationError(
            "No rows with non-zero coverage on both WGBS and FinaleMe sides "
            "after the join. Check coverage and overlap of the input files."
        )

    # Compute beta values from k/n on each side
    wgbs_k = merged["methylated_count_wgbs"].to_numpy(dtype=np.float64)
    wgbs_n = merged["total_count_wgbs"].to_numpy(dtype=np.float64)
    fme_k = merged["methylated_count_fme"].to_numpy(dtype=np.float64)
    fme_n = merged["total_count_fme"].to_numpy(dtype=np.float64)
    wgbs_beta = (wgbs_k / wgbs_n).astype(np.float64)
    fme_beta = (fme_k / fme_n).astype(np.float64)
    density = merged["cpg_density"].to_numpy(dtype=np.float64)
    if "region_class" in merged.columns:
        region_class = merged["region_class"].astype(object).to_numpy()
    else:
        region_class = None
    chrom = merged["chrom"].astype(str).to_numpy()
    sample_ids = merged["sample_id"].astype(str).to_numpy()

    # Bin tuning via chromosome-blocked K-fold CV
    tune_result = tune_n_bins(
        predicted=fme_beta,
        truth_beta=wgbs_beta,
        cpg_density=density,
        chrom=chrom,
        n_bins_candidates=n_bins_candidates,
        region_class=region_class,
        n_folds=cv_n_folds,
        seed=cv_seed,
        max_error_rate=max_error_rate,
        threads=threads,
    )
    best_B = int(tune_result["selected_n_bins"])

    # Final fit at the chosen B (using the in_sample_fit from the chosen
    # candidate so we don't redo the work)
    final_candidate = next(
        c for c in tune_result["candidates"] if c["n_bins"] == best_B
    )
    final_fit = final_candidate["in_sample_fit"]
    params = final_fit.params

    # Attach training provenance
    params.training_metadata = {
        "n_training_samples": int(len(np.unique(sample_ids))),
        "n_observations": int(len(merged)),
        "n_bins_candidates": list(n_bins_candidates),
        "selected_n_bins": best_B,
        "cv_method": cv_method,
        "cv_n_folds": int(cv_n_folds),
        "cv_seed": cv_seed,
        "max_error_rate": float(max_error_rate),
        "tune_result": {
            "selected_n_bins": tune_result["selected_n_bins"],
            "candidates": [
                {
                    "n_bins": c["n_bins"],
                    "score": c.get("score"),
                    "cv_accuracy": c.get("cv_accuracy"),
                    "cv_accuracy_std": c.get("cv_accuracy_std"),
                    "in_sample_accuracy": c.get("in_sample_accuracy"),
                    "n_usable_markers": c.get("n_usable_markers"),
                    "n_total_markers": c.get("n_total_markers"),
                    "n_folds": c.get("n_folds"),
                }
                for c in tune_result["candidates"]
            ],
        },
    }
    params.save(out_params)

    # Write a JSON training report. Reuse the calibration report writer —
    # it's just a JSON file, no schema-specific behavior.
    report = {
        "binarization_version": "1.0",
        "n_training_samples": int(len(np.unique(sample_ids))),
        "n_observations": int(len(merged)),
        "n_bins": best_B,
        "n_bins_candidates": list(n_bins_candidates),
        "cv_method": cv_method,
        "cv_n_folds": int(cv_n_folds),
        "cv_seed": cv_seed,
        "max_error_rate": float(max_error_rate),
        "selected_score": final_candidate.get("score"),
        "selected_cv_accuracy": final_candidate.get("cv_accuracy"),
        "overall_metrics": final_fit.overall,
        "per_bin_metrics": final_fit.per_bin_metrics,
        "candidates": params.training_metadata["tune_result"]["candidates"],
    }
    write_calibration_report(report, out_report)
    return params


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
    "train_binarization",
]
