"""QC flag generation for per-sample TOO results (architecture §10.6)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from finaleme_too.config import CoverageTier, QCConfig

if TYPE_CHECKING:
    from finaleme_too.core.deconvolution import DeconvolutionResult
    from finaleme_too.core.observation_model import ObservationModel


def compute_qc_flags(
    result: "DeconvolutionResult",
    observation: "ObservationModel",
    qc_config: QCConfig,
    binarization_flag: str | None = None,
    hemolysis: bool | None = None,
) -> list[str]:
    """Generate the qc_flags list for a single sample.

    The ``binarization_flag`` parameter was renamed from ``calibration_flag``
    in v3 to reflect the FinaleMe path now using context-dependent
    binarization. The flag values (PASS / WARN / FAIL) and their semantics
    are unchanged — they describe the per-sample inference QC for the
    FinaleMe preprocessing step.
    """
    flags: list[str] = []

    # Unknown component too high
    unknown = float(result.proportions[-1])
    if unknown > qc_config.max_unknown_fraction:
        flags.append("HIGH_UNKNOWN")

    # Coverage tier
    if result.coverage_tier == CoverageTier.ULTRALOW:
        flags.append("ULTRALOW_COVERAGE")
    elif result.coverage_tier == CoverageTier.LOW:
        flags.append("LOW_COVERAGE")

    # WBC fraction (sum of blood-related cell types named "Blood-*")
    blood_cts = [i for i, ct in enumerate(result.cell_types) if ct.lower().startswith("blood")]
    if blood_cts:
        wbc = float(np.sum(result.proportions[blood_cts]))
        if wbc > qc_config.max_wbc_fraction:
            flags.append("WBC_DOMINANT")

    # Binarization flag (FinaleMe mode). Renamed from CALIBRATION_* in v3.
    if binarization_flag in ("WARN", "FAIL"):
        flags.append(f"BINARIZATION_{binarization_flag}")

    # Hemolysis (from sample sheet)
    if hemolysis is True:
        flags.append("HEMOLYSIS")

    # No reliable cell types — everything is LOW or UNRELIABLE
    reliable = sum(
        1 for r in result.reliability[:-1] if r in ("HIGH", "MODERATE")
    )
    if reliable == 0:
        flags.append("NO_RELIABLE_CELLTYPES")

    return flags


__all__ = ["compute_qc_flags"]
