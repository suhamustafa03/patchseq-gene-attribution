"""Schema validation and basic QC filters, applied before any modeling."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

CELLS_REQUIRED_COLUMNS: dict[str, type] = {
    "cell_id": str,
    "species": str,
    "subclass": str,
    "donor_id": str,
    "has_transcriptome": bool,
}

EPHYS_REQUIRED_COLUMNS: dict[str, type] = {
    "cell_id": str,
    "rheobase_pA": float,
    "ap_width_ms": float,
    "input_resistance_MOhm": float,
    "adaptation_ratio": float,
}

VALID_SPECIES = {"mouse", "human"}
# "unknown" added for the real-data run (2026-11-23): DANDI:000636 (human
# Patch-seq ephys) does not carry a confirmed transcriptomic subclass label,
# so real human cells are tagged subclass="unknown" rather than fabricating
# a label. Demo-fixture validation (Pvalb/Sst/Vip/Lamp5/Sncg only) is
# unaffected since demo_fixture/ never uses "unknown".
VALID_SUBCLASSES = {"Pvalb", "Sst", "Vip", "Lamp5", "Sncg", "unknown"}

EPHYS_FEATURE_COLUMNS = [
    "rheobase_pA",
    "ap_width_ms",
    "input_resistance_MOhm",
    "adaptation_ratio",
]


class SchemaError(ValueError):
    """Raised when a table does not match the expected pipeline schema."""


@dataclass
class ValidationReport:
    table_name: str
    n_rows_in: int
    n_rows_out: int = 0
    dropped_missing_key_features: int = 0
    errors: list[str] = field(default_factory=list)

    @property
    def is_valid(self) -> bool:
        return not self.errors


def _check_required_columns(df: pd.DataFrame, required: dict[str, type], table_name: str) -> list[str]:
    errors = []
    missing = [c for c in required if c not in df.columns]
    if missing:
        errors.append(f"{table_name}: missing required columns {missing}")
    return errors


def _check_no_duplicate_ids(df: pd.DataFrame, table_name: str) -> list[str]:
    if "cell_id" not in df.columns:
        return []
    dupes = df["cell_id"][df["cell_id"].duplicated()].unique().tolist()
    if dupes:
        return [f"{table_name}: duplicate cell_id values found: {dupes[:10]}"]
    return []


def validate_cells(cells: pd.DataFrame) -> ValidationReport:
    report = ValidationReport(table_name="cells", n_rows_in=len(cells))
    report.errors.extend(_check_required_columns(cells, CELLS_REQUIRED_COLUMNS, "cells"))
    report.errors.extend(_check_no_duplicate_ids(cells, "cells"))

    if not report.errors:
        bad_species = set(cells["species"].unique()) - VALID_SPECIES
        if bad_species:
            report.errors.append(f"cells: unexpected species values {bad_species}")
        bad_subclass = set(cells["subclass"].unique()) - VALID_SUBCLASSES
        if bad_subclass:
            report.errors.append(f"cells: unexpected subclass values {bad_subclass}")

    report.n_rows_out = len(cells) if report.is_valid else 0
    return report


def validate_ephys(ephys: pd.DataFrame) -> ValidationReport:
    report = ValidationReport(table_name="ephys_features", n_rows_in=len(ephys))
    report.errors.extend(_check_required_columns(ephys, EPHYS_REQUIRED_COLUMNS, "ephys_features"))
    report.errors.extend(_check_no_duplicate_ids(ephys, "ephys_features"))
    report.n_rows_out = len(ephys) if report.is_valid else 0
    return report


def require_valid(report: ValidationReport) -> None:
    if not report.is_valid:
        raise SchemaError("; ".join(report.errors))


def drop_missing_key_features(
    ephys: pd.DataFrame, feature_columns: list[str] = EPHYS_FEATURE_COLUMNS
) -> tuple[pd.DataFrame, int]:
    before = len(ephys)
    mask = ephys[feature_columns].notna().all(axis=1) & np.isfinite(
        ephys[feature_columns].to_numpy(dtype=float)
    ).all(axis=1)
    cleaned = ephys.loc[mask].reset_index(drop=True)
    return cleaned, before - len(cleaned)


def run_qc(cells: pd.DataFrame, ephys: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, list[ValidationReport]]:
    cells_report = validate_cells(cells)
    ephys_report = validate_ephys(ephys)
    require_valid(cells_report)
    require_valid(ephys_report)

    ephys_clean, n_dropped = drop_missing_key_features(ephys)
    ephys_report.n_rows_out = len(ephys_clean)
    ephys_report.dropped_missing_key_features = n_dropped

    valid_cell_ids = set(cells["cell_id"])
    ephys_clean = ephys_clean[ephys_clean["cell_id"].isin(valid_cell_ids)].reset_index(drop=True)

    return cells, ephys_clean, [cells_report, ephys_report]
