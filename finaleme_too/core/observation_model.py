"""Beta-binomial observation model.

Wraps per-marker counts (k_i, n_i) with per-marker dispersion phi_i and
weight omega_i, ready to be passed to the deconvolution solver. Implements
the formulas in TOO_MATH_FORMULATION_v2.md §2 and §3.2.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from finaleme_too.config import CoverageTier, MeasurementMode
from finaleme_too.io.methylation_loader import MarkerObservations
from finaleme_too.utils.beta_binomial import (
    estimate_dispersion_mle,
    log_likelihood_per_marker,
)

if TYPE_CHECKING:
    from finaleme_too.io.reference_panel import ReferencePanel


@dataclass(frozen=True)
class ObservationModel:
    """Per-sample data + per-marker dispersion + per-marker weight."""

    sample_id: str
    k: np.ndarray  # int64, (M,)
    n: np.ndarray  # int64, (M,)
    mu_obs: np.ndarray  # float64, (M,) — k/n where defined, NaN otherwise
    dispersion: np.ndarray  # float64, (M,)
    weights: np.ndarray  # float64, (M,)
    mode: MeasurementMode
    coverage_tier: CoverageTier
    coverage_cap: int

    @property
    def n_markers(self) -> int:
        return len(self.k)

    def log_likelihood(self, mu: np.ndarray) -> np.ndarray:
        """Vectorized per-marker log-likelihood given expected mu."""
        return log_likelihood_per_marker(self.k, self.n, mu, self.dispersion)

    def total_log_likelihood(self, mu: np.ndarray) -> float:
        return float(np.sum(self.weights * self.log_likelihood(mu)))


# Default dispersion ranges per (mode, tier) — architecture §6.2
_DEFAULT_DISPERSION = {
    (MeasurementMode.WGBS, CoverageTier.HIGH): 100.0,
    (MeasurementMode.WGBS, CoverageTier.LOW): 50.0,
    (MeasurementMode.WGBS, CoverageTier.ULTRALOW): 10.0,
    (MeasurementMode.FINALEME, CoverageTier.HIGH): 15.0,
    (MeasurementMode.FINALEME, CoverageTier.LOW): 8.0,
    (MeasurementMode.FINALEME, CoverageTier.ULTRALOW): 3.0,
}


class BetaBinomialModel:
    """Build an ObservationModel from raw MarkerObservations."""

    def build(
        self,
        obs: MarkerObservations,
        reference: "ReferencePanel | None" = None,
        calibration: object | None = None,  # CalibrationParams (P1)
        region_annotations: pd.DataFrame | None = None,
        tier: CoverageTier = CoverageTier.HIGH,
        coverage_cap: int = 50,
        balance_weights: np.ndarray | None = None,
    ) -> ObservationModel:
        k = np.asarray(obs.k, dtype=np.int64)
        n = np.asarray(obs.n, dtype=np.int64)
        with np.errstate(invalid="ignore", divide="ignore"):
            mu_obs = np.where(n > 0, k / np.maximum(n, 1), np.nan).astype(np.float64)

        # Dispersion: per-marker, defaults from mode+tier table.
        # WGBS HIGH gets a one-shot MLE of a single shared phi from the data.
        phi = self._compute_dispersion(
            obs=obs,
            mode=obs.mode,
            tier=tier,
            calibration=calibration,
            region_annotations=region_annotations,
            mu_obs=mu_obs,
        )

        # Marker weights ω_i = min(n_i, n_cap)/n_cap * balance_weight_i  (§3.2)
        capped = np.minimum(n.astype(np.float64), float(coverage_cap)) / float(coverage_cap)
        if balance_weights is None:
            balance = np.ones_like(k, dtype=np.float64)
        else:
            balance = np.asarray(balance_weights, dtype=np.float64)
        weights = capped * balance

        # Reference uncertainty (§2.5): inflate phi if reference coverage is low.
        if reference is not None and reference.coverage is not None:
            ref_min = np.min(reference.coverage, axis=1).astype(np.float64)
            n_ref_cap = max(coverage_cap, 1)
            scale = np.minimum(1.0, ref_min / float(n_ref_cap))
            phi = phi * np.maximum(scale, 1e-3)

        return ObservationModel(
            sample_id=obs.sample_id,
            k=k,
            n=n,
            mu_obs=mu_obs,
            dispersion=phi,
            weights=weights,
            mode=obs.mode,
            coverage_tier=tier,
            coverage_cap=coverage_cap,
        )

    def _compute_dispersion(
        self,
        obs: MarkerObservations,
        mode: MeasurementMode,
        tier: CoverageTier,
        calibration,
        region_annotations: pd.DataFrame | None,
        mu_obs: np.ndarray,
    ) -> np.ndarray:
        n_markers = obs.n_markers
        # FinaleMe mode + calibration: use per-bin log-dispersion when available.
        if mode == MeasurementMode.FINALEME and calibration is not None and region_annotations is not None:
            try:
                # Lazy attribute access to avoid hard import here
                bin_edges = np.asarray(calibration.bin_edges, dtype=np.float64)
                log_phi = np.asarray(calibration.log_dispersion, dtype=np.float64)
                # Map each marker to a bin via cpg_density
                density = region_annotations.set_index(["chrom", "start", "end"])["cpg_density"]
                key = list(zip(obs.chrom.tolist(), obs.start.tolist(), obs.end.tolist()))
                marker_density = np.asarray(
                    [float(density.get(k, np.nan)) for k in key], dtype=np.float64
                )
                bin_idx = np.clip(
                    np.searchsorted(bin_edges, marker_density, side="right") - 1,
                    0,
                    len(log_phi) - 1,
                )
                phi = np.exp(log_phi[bin_idx])
                # §2.4: scale by g(n_i) = min(n_i, n_cap)/n_cap with n_cap = 50
                cap = 50.0
                g_n = np.minimum(np.asarray(obs.n, dtype=np.float64), cap) / cap
                phi = phi * np.where(g_n > 0, g_n, 1e-3)
                return phi.astype(np.float64)
            except Exception:
                pass

        # WGBS HIGH tier: estimate a single shared phi via MLE on the data.
        if mode == MeasurementMode.WGBS and tier == CoverageTier.HIGH:
            valid = (~np.isnan(mu_obs)) & (obs.n > 0)
            if int(np.sum(valid)) > 50:
                shared_phi = estimate_dispersion_mle(
                    k=obs.k[valid].astype(np.float64),
                    n=obs.n[valid].astype(np.float64),
                    mu=mu_obs[valid],
                    phi_init=100.0,
                )
                return np.full(n_markers, shared_phi, dtype=np.float64)

        # Otherwise: use the default per-(mode,tier) value.
        default_phi = _DEFAULT_DISPERSION.get((mode, tier), 50.0)
        return np.full(n_markers, default_phi, dtype=np.float64)


__all__ = ["BetaBinomialModel", "ObservationModel"]
