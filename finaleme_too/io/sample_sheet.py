"""Sample sheet parser and validator (architecture §3.1)."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from finaleme_too.config import MeasurementMode
from finaleme_too.exceptions import InvalidSampleSheetError

REQUIRED_COLUMNS = {"sample_id", "methylation_file", "mode"}
OPTIONAL_COLUMNS = {
    "input_format",
    "group",
    "pat_file",  # optional .pat.gz path for fragment-level (ULTRALOW) mode
    "extraction_batch",
    "library_date",
    "sequencing_date",
    "institute",
    "plasma_draw_date",
    "age",
    "sex",
    "race",
    "ethnicity",
    "treatment",
    "treatment_efficacy",
    "mutation_status",
    "disease_duration",
    "follow_up_time",
    "smoking_status",
    "alcohol_usage",
    "bmi",
    "cfdna_concentration_ng_ml",
    "draw_time",
    "hemolysis_flag",
    "fasting_status",
    "meth_col",
    "total_col",
}


@dataclass
class Sample:
    sample_id: str
    methylation_file: Path
    mode: MeasurementMode
    input_format: str | None = None
    group: str | None = None
    meth_col: int | None = None
    total_col: int | None = None
    pat_file: Path | None = None
    metadata: dict = field(default_factory=dict)

    def resolved_pat_file(self) -> Path | None:
        """Return the .pat.gz companion path for this sample.

        Prefers the explicit ``pat_file`` sample-sheet column; otherwise
        falls back to the conventional naming where FinaleMe writes the
        fragment patterns next to the prediction BED (swapping
        ``.prediction.bed.gz`` → ``.prediction.pat.gz``).
        """
        if self.pat_file is not None:
            return self.pat_file
        p = self.methylation_file
        name = p.name
        for suffix in (".prediction.bed.gz", ".decode.prediction.bed.gz"):
            if name.endswith(suffix):
                candidate = p.with_name(name[: -len(suffix)] + suffix.replace(".bed.gz", ".pat.gz"))
                if candidate.exists():
                    return candidate
        return None

    @property
    def technical_covariates(self) -> dict[str, str]:
        keys = ("extraction_batch", "library_date", "sequencing_date", "institute", "plasma_draw_date")
        return {k: self.metadata[k] for k in keys if k in self.metadata}

    @property
    def biological_covariates(self) -> dict:
        biological_keys = (
            "age",
            "sex",
            "race",
            "ethnicity",
            "treatment",
            "treatment_efficacy",
            "mutation_status",
            "disease_duration",
            "follow_up_time",
            "smoking_status",
            "alcohol_usage",
            "bmi",
            "draw_time",
            "fasting_status",
        )
        return {k: self.metadata[k] for k in biological_keys if k in self.metadata}


@dataclass
class SampleSheet:
    samples: list[Sample]
    raw_table: pd.DataFrame

    def __len__(self) -> int:
        return len(self.samples)

    def __iter__(self):
        return iter(self.samples)

    @classmethod
    def from_tsv(cls, path: str | Path) -> SampleSheet:
        path = Path(path)
        if not path.exists():
            raise InvalidSampleSheetError(f"Sample sheet not found: {path}")
        df = pd.read_csv(path, sep="\t", comment="#")
        return cls.from_dataframe(df)

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame) -> SampleSheet:
        missing = REQUIRED_COLUMNS - set(df.columns)
        if missing:
            raise InvalidSampleSheetError(
                f"Sample sheet missing required columns: {sorted(missing)}"
            )
        # Validate uniqueness of sample_id
        if df["sample_id"].duplicated().any():
            dups = df[df["sample_id"].duplicated()]["sample_id"].tolist()
            raise InvalidSampleSheetError(f"Duplicate sample_id values: {dups}")

        samples: list[Sample] = []
        for _, row in df.iterrows():
            mode_str = str(row["mode"]).strip().upper()
            try:
                mode = MeasurementMode(mode_str)
            except ValueError as exc:
                raise InvalidSampleSheetError(
                    f"Invalid mode for sample {row['sample_id']}: {mode_str}"
                ) from exc
            metadata = {
                k: row[k]
                for k in df.columns
                if k not in REQUIRED_COLUMNS and k not in {"input_format", "meth_col", "total_col"}
                and pd.notna(row[k])
            }
            samples.append(
                Sample(
                    sample_id=str(row["sample_id"]),
                    methylation_file=Path(str(row["methylation_file"])),
                    mode=mode,
                    input_format=str(row["input_format"]) if "input_format" in df.columns and pd.notna(row.get("input_format")) else None,
                    group=str(row["group"]) if "group" in df.columns and pd.notna(row.get("group")) else None,
                    meth_col=int(row["meth_col"]) if "meth_col" in df.columns and pd.notna(row.get("meth_col")) else None,
                    total_col=int(row["total_col"]) if "total_col" in df.columns and pd.notna(row.get("total_col")) else None,
                    pat_file=Path(str(row["pat_file"])) if "pat_file" in df.columns and pd.notna(row.get("pat_file")) else None,
                    metadata=metadata,
                )
            )
        return cls(samples=samples, raw_table=df)

    def validate_files_exist(self) -> None:
        missing = [s.sample_id for s in self.samples if not s.methylation_file.exists()]
        if missing:
            raise InvalidSampleSheetError(
                f"Methylation files not found for samples: {missing}"
            )

    def groups(self) -> list[str]:
        return sorted({s.group for s in self.samples if s.group is not None})


__all__ = ["Sample", "SampleSheet"]
