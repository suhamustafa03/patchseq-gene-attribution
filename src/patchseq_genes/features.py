"""Load/merge cells + ephys + expression tables, keeping the two analysis
arms (attribution vs. ephys-only cross-species comparison) separate.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .qc import EPHYS_FEATURE_COLUMNS


def load_tables(input_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    input_dir = Path(input_dir)
    cells = pd.read_csv(input_dir / "cells.csv")
    ephys = pd.read_csv(input_dir / "ephys_features.csv")
    expression = pd.read_csv(input_dir / "gene_expression.csv")
    cells["has_transcriptome"] = cells["has_transcriptome"].astype(bool)
    return cells, ephys, expression


def get_gene_columns(expression: pd.DataFrame) -> list[str]:
    return [c for c in expression.columns if c != "cell_id"]


def build_attribution_table(
    cells: pd.DataFrame, ephys: pd.DataFrame, expression: pd.DataFrame
) -> pd.DataFrame:
    """Mouse-only table for gene -> ephys attribution.

    Restricted to has_transcriptome=True because expression is only ever
    available for mouse cells in this project's scope (no human
    transcriptome): merging first and filtering after would silently
    include human rows with all-NaN gene columns.
    """
    mouse = cells[cells["has_transcriptome"]]
    merged = mouse.merge(ephys, on="cell_id", how="inner").merge(
        expression, on="cell_id", how="inner"
    )
    return merged.reset_index(drop=True)


def build_cross_species_ephys_table(cells: pd.DataFrame, ephys: pd.DataFrame) -> pd.DataFrame:
    """All cells, both species, ephys features only -- never expression.

    This is the ephys-only phenotypic comparison arm; expression columns
    must never appear here, since DANDI:000636 (human) has no
    transcriptome and mixing arms would implicitly claim cross-species
    gene attribution that is out of scope.
    """
    merged = cells.merge(ephys, on="cell_id", how="inner")
    return merged[["cell_id", "species", "subclass", "donor_id"] + EPHYS_FEATURE_COLUMNS].reset_index(
        drop=True
    )
