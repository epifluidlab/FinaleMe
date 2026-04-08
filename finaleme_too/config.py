"""Configuration dataclasses and YAML loading for finaleme_too.

Mirrors the YAML schema in design/TOO_ARCHITECTURE_v2.md §11.2.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields, asdict
from enum import Enum
from pathlib import Path
from typing import Any

import yaml


class MeasurementMode(str, Enum):
    """Measurement modality of the methylation data."""

    WGBS = "WGBS"
    FINALEME = "FINALEME"


class CoverageTier(str, Enum):
    """Per-sample coverage tier (architecture §4)."""

    HIGH = "HIGH"
    LOW = "LOW"
    ULTRALOW = "ULTRALOW"


class SolverMethod(str, Enum):
    MLE = "mle"
    BAYESIAN = "bayesian"


class TestMethod(str, Enum):
    ILR_REGRESSION = "ilr_regression"
    BAYESIAN_POSTERIOR = "bayesian_posterior"
    WILCOXON = "wilcoxon"


# ----------------------------------------------------------------------------
# Sub-configs (one per top-level YAML key)
# ----------------------------------------------------------------------------


@dataclass
class ModelConfig:
    observation: str = "beta_binomial"
    deconvolution: SolverMethod = SolverMethod.MLE
    unknown_component: bool = True  # always on; flag kept for documentation
    fragment_level: str = "auto"  # "auto", "always", "never"


@dataclass
class CoverageConfig:
    tier_high: float = 10.0
    tier_low: float = 0.5
    min_reads: int = 3
    coverage_cap: int = 50


@dataclass
class CalibrationConfig:
    mode: MeasurementMode = MeasurementMode.FINALEME
    calibration_file: str | None = None
    n_density_bins: int = 8
    use_default: bool = True  # use shipped default if calibration_file is None
    # CV strategy used during bin tuning (train-calibration only).
    #   "auto"   — leave-one-sample-out when >=2 samples, else random region K-fold
    #   "sample" — always leave-one-sample-out (returns NaN cv_rmse if <2 samples)
    #   "region" — always random K-fold on row indices (works with 1 sample)
    #   "none"   — skip CV, select by AIC / in-sample RMSE only
    cv_strategy: str = "auto"
    cv_n_folds: int = 10  # K for region-level CV (ignored for sample mode)
    cv_seed: int | None = None  # reproducible region-level shuffle


@dataclass
class MarkersConfig:
    marker_regions: str | None = None
    marker_format: str = "auto"  # "auto", "bed", "uxm_atlas"
    region_annotation: str | None = None
    strict_regions: str | None = None  # e.g. "CGI+shore"
    n_per_type: int = 500
    specificity_method: str = "entropy"  # "entropy", "t_statistic", "delta_means"


@dataclass
class UncertaintyConfig:
    method: str = "bootstrap"  # "bootstrap", "bayesian", "both", "none"
    n_bootstrap: int = 1000
    ci_level: float = 0.95
    noise_floor: float = 0.001
    bayesian_prior_alpha: float = 1.0
    bayesian_n_samples: int = 5000
    # Optional integer seed for the bootstrap rng. ``None`` (the default)
    # uses fresh entropy on every run; setting this enables reproducible
    # bootstrap CIs across both serial and threaded executions.
    seed: int | None = None


@dataclass
class BatchCorrectionConfig:
    technical_covariates: list[str] = field(default_factory=list)
    min_levels: int = 2
    min_samples_per_level: int = 5


@dataclass
class CovariateAdjustmentConfig:
    biological_covariates: list[str] = field(default_factory=list)
    user_configurable: list[str] = field(
        default_factory=lambda: ["treatment", "treatment_efficacy", "mutation_status"]
    )
    transform: str = "ILR"


@dataclass
class TestingConfig:
    method: TestMethod = TestMethod.ILR_REGRESSION
    group_comparison: str = "omnibus+pairwise"
    fdr_method: str = "BH"
    fdr_alpha: float = 0.05


@dataclass
class InputConfig:
    format: str = "auto"  # "auto", "finaleme_bed", "bissnp_6plus2", "wgbstools_beta", "custom_bed"
    meth_col: int | None = None
    total_col: int | None = None


@dataclass
class ReferenceConfig:
    format: str = "matrix"  # "matrix" or "beta_list"
    reference_panel: str | None = None
    ref_betas: str | None = None
    ref_groups: str | None = None
    cpg_index: str | None = None
    coverage_matrix: str | None = None


@dataclass
class QCConfig:
    max_wbc_fraction: float = 0.95
    max_unknown_fraction: float = 0.30
    max_residual_variance: float = 0.40


@dataclass
class TOOConfig:
    """Top-level finaleme-too configuration."""

    model: ModelConfig = field(default_factory=ModelConfig)
    coverage: CoverageConfig = field(default_factory=CoverageConfig)
    calibration: CalibrationConfig = field(default_factory=CalibrationConfig)
    markers: MarkersConfig = field(default_factory=MarkersConfig)
    uncertainty: UncertaintyConfig = field(default_factory=UncertaintyConfig)
    batch_correction: BatchCorrectionConfig = field(default_factory=BatchCorrectionConfig)
    covariate_adjustment: CovariateAdjustmentConfig = field(
        default_factory=CovariateAdjustmentConfig
    )
    testing: TestingConfig = field(default_factory=TestingConfig)
    input: InputConfig = field(default_factory=InputConfig)
    reference: ReferenceConfig = field(default_factory=ReferenceConfig)
    qc: QCConfig = field(default_factory=QCConfig)
    threads: int = 1

    @classmethod
    def from_yaml(cls, path: str | Path) -> TOOConfig:
        """Load a TOOConfig from a YAML file. Missing sections fall back to defaults."""
        with open(path) as fh:
            raw = yaml.safe_load(fh) or {}
        return cls.from_dict(raw)

    @classmethod
    def from_dict(cls, raw: dict[str, Any]) -> TOOConfig:
        """Build a TOOConfig from a nested dict (e.g. parsed YAML)."""
        cfg = cls()
        for f in fields(cls):
            if f.name not in raw:
                continue
            section_value = raw[f.name]
            if isinstance(section_value, dict):
                section_cls = type(getattr(cfg, f.name))
                if hasattr(section_cls, "__dataclass_fields__"):
                    setattr(cfg, f.name, _build_subconfig(section_cls, section_value))
                    continue
            setattr(cfg, f.name, section_value)
        return cfg

    def to_dict(self) -> dict[str, Any]:
        return _to_serializable(asdict(self))


def _resolve_field_enum_type(cls, field_obj) -> type | None:
    """Return the Enum subclass for a dataclass field, or None.

    With ``from __future__ import annotations`` everywhere, ``field.type``
    is a string and may use PEP 604 ``X | None`` syntax that fails to
    evaluate on Python 3.9-. We try ``typing.get_type_hints`` first and
    fall back to inspecting the default value's actual type, which is
    always the enum class for enum-typed fields.
    """
    import typing
    from dataclasses import MISSING

    # 1) Try the proper type-hint resolution path
    try:
        hints = typing.get_type_hints(cls)
        resolved = hints.get(field_obj.name)
        if isinstance(resolved, type) and issubclass(resolved, Enum):
            return resolved
        if resolved is not None:
            origin = typing.get_origin(resolved)
            if origin is typing.Union:
                for arg in typing.get_args(resolved):
                    if isinstance(arg, type) and issubclass(arg, Enum):
                        return arg
    except Exception:
        pass

    # 2) Fall back to inspecting the default value
    if field_obj.default is not MISSING and isinstance(field_obj.default, Enum):
        return type(field_obj.default)
    if field_obj.default_factory is not MISSING:  # type: ignore[misc]
        try:
            sample = field_obj.default_factory()  # type: ignore[misc]
        except Exception:
            sample = None
        if isinstance(sample, Enum):
            return type(sample)

    return None


def _build_subconfig(cls, raw: dict[str, Any]):
    """Build a sub-config dataclass, coercing enum-typed fields."""
    kwargs: dict[str, Any] = {}
    valid_fields = {f.name: f for f in fields(cls)}
    for key, value in raw.items():
        if key not in valid_fields:
            continue
        f = valid_fields[key]
        enum_cls = _resolve_field_enum_type(cls, f)
        if enum_cls is not None:
            if value is None:
                kwargs[key] = None
            elif isinstance(value, enum_cls) and not isinstance(value, str):
                # Genuine enum instance — pass through
                kwargs[key] = value
            else:
                kwargs[key] = enum_cls(value)
            continue
        kwargs[key] = value
    return cls(**kwargs)


def _to_serializable(obj: Any) -> Any:
    if isinstance(obj, Enum):
        return obj.value
    if isinstance(obj, dict):
        return {k: _to_serializable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_to_serializable(v) for v in obj]
    return obj


__all__ = [
    "MeasurementMode",
    "CoverageTier",
    "SolverMethod",
    "TestMethod",
    "ModelConfig",
    "CoverageConfig",
    "CalibrationConfig",
    "MarkersConfig",
    "UncertaintyConfig",
    "BatchCorrectionConfig",
    "CovariateAdjustmentConfig",
    "TestingConfig",
    "InputConfig",
    "ReferenceConfig",
    "QCConfig",
    "TOOConfig",
]
