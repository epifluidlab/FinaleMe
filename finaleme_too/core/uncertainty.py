"""Bootstrap CI estimation for deconvolution proportions."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from finaleme_too.core.deconvolution import MLEDeconvolver
from finaleme_too.core.observation_model import ObservationModel

if TYPE_CHECKING:
    from finaleme_too.io.reference_panel import ReferencePanel


@dataclass
class BootstrapResult:
    proportions_samples: np.ndarray  # (B, K+1)
    ci_lower: np.ndarray  # (K+1,)
    ci_upper: np.ndarray  # (K+1,)
    point_estimate: np.ndarray  # (K+1,) — mean over bootstrap


class BootstrapCI:
    """Marker-resampling bootstrap for proportion CIs."""

    def __init__(
        self,
        n_iterations: int = 1000,
        ci_level: float = 0.95,
        seed: int | None = None,
    ):
        self.n_iterations = n_iterations
        self.ci_level = ci_level
        self.seed = seed

    def estimate(
        self,
        model: ObservationModel,
        reference: "ReferencePanel",
        deconvolver: MLEDeconvolver,
    ) -> BootstrapResult:
        rng = np.random.default_rng(self.seed)
        n_markers = model.n_markers
        # K_total includes the unknown component
        K_total = reference.n_cell_types + 1
        samples = np.empty((self.n_iterations, K_total), dtype=np.float64)

        for b in range(self.n_iterations):
            idx = rng.integers(0, n_markers, size=n_markers)
            w_b = deconvolver.solve(model, reference, marker_subset=idx)
            samples[b] = w_b

        alpha = (1.0 - self.ci_level) / 2.0
        ci_lower = np.quantile(samples, alpha, axis=0)
        ci_upper = np.quantile(samples, 1.0 - alpha, axis=0)
        point = np.mean(samples, axis=0)
        return BootstrapResult(
            proportions_samples=samples,
            ci_lower=ci_lower,
            ci_upper=ci_upper,
            point_estimate=point,
        )


__all__ = ["BootstrapCI", "BootstrapResult"]
