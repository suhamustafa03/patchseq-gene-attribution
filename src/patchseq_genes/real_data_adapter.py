"""Adapter: real mouse (Allen Cell Types API) + real human (DANDI:000636)
ephys tables -> this pipeline's QC + cross-species-ephys-only arm.

This does NOT run the gene-attribution arm: no real, honestly-obtained
linked mouse transcriptome was available (see docs/DATA_MANIFEST.md), so
there is no gene_expression table here and build_attribution_table() is
never called in this script.

Usage:
    python -m patchseq_genes.real_data_adapter \
        --mouse-cells data/raw/mouse_gouwens/cells_real.csv \
        --mouse-ephys data/raw/mouse_gouwens/ephys_features_real.csv \
        --human-cells data/raw/human_dandi000636/cells_real.csv \
        --human-ephys data/raw/human_dandi000636/ephys_features_real.csv \
        --output results/real
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from . import qc
from .confounders import build_confounder_matrix, donor_groups
from .cv import assert_no_donor_leakage, donor_grouped_train_test_split
from .features import build_cross_species_ephys_table


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mouse-cells", type=Path, required=True)
    p.add_argument("--mouse-ephys", type=Path, required=True)
    p.add_argument("--human-cells", type=Path, required=True)
    p.add_argument("--human-ephys", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--random-state", type=int, default=0)
    return p.parse_args()


def run(
    mouse_cells_path: Path,
    mouse_ephys_path: Path,
    human_cells_path: Path,
    human_ephys_path: Path,
    output_dir: Path,
    random_state: int = 0,
) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)

    mouse_cells = pd.read_csv(mouse_cells_path)
    mouse_ephys = pd.read_csv(mouse_ephys_path)
    human_cells = pd.read_csv(human_cells_path)
    human_ephys = pd.read_csv(human_ephys_path)

    cells = pd.concat([mouse_cells, human_cells], ignore_index=True)
    ephys = pd.concat([mouse_ephys, human_ephys], ignore_index=True)
    cells["has_transcriptome"] = cells["has_transcriptome"].astype(bool)

    cells, ephys_clean, qc_reports = qc.run_qc(cells, ephys)
    for report in qc_reports:
        print(
            f"[qc] {report.table_name}: in={report.n_rows_in} out={report.n_rows_out} "
            f"dropped_missing={report.dropped_missing_key_features}"
        )

    cross_species_table = build_cross_species_ephys_table(cells, ephys_clean)
    cross_species_table.to_csv(output_dir / "cross_species_ephys_table.csv", index=False)

    # Donor-grouped split check, exercised on REAL donor ids (not synthetic).
    groups = donor_groups(cross_species_table)
    donor_leakage_ok = False
    donor_leakage_error = None
    n_train_donors = n_test_donors = 0
    try:
        train_idx, test_idx = donor_grouped_train_test_split(groups, random_state=random_state)
        assert_no_donor_leakage(groups, train_idx, test_idx)
        donor_leakage_ok = True
        n_train_donors = len(set(groups[train_idx]))
        n_test_donors = len(set(groups[test_idx]))
    except Exception as exc:  # pragma: no cover - reported, not swallowed
        donor_leakage_error = str(exc)

    confounders = build_confounder_matrix(cross_species_table, include_species=True)
    confounders.to_csv(output_dir / "cross_species_confounders.csv", index=False)

    n_mouse = int((cells["species"] == "mouse").sum())
    n_human = int((cells["species"] == "human").sum())
    n_donors_mouse = cells.loc[cells["species"] == "mouse", "donor_id"].nunique()
    n_donors_human = cells.loc[cells["species"] == "human", "donor_id"].nunique()
    subclass_counts = cells["subclass"].value_counts().to_dict()

    per_feature_species_means = (
        cross_species_table.groupby("species")[qc.EPHYS_FEATURE_COLUMNS].mean().round(4).to_dict()
    )

    summary = {
        "status": "real_data_ephys_only_arm",
        "arms_run": ["human_ephys_only_cross_species_comparison"],
        "arms_blocked": {
            "mouse_gene_attribution": (
                "No real, honestly-licensed linked mouse transcriptome was obtained: "
                "NeMO dat-m2kd56v (transcriptome) now shows 'No license listed' on "
                "re-check (2026-11-23), and no confirmed cell-linked accession was "
                "found for this cohort during prior project scoping. Gene expression "
                "input is therefore absent and build_attribution_table()/stability "
                "selection/attribution models were never invoked on real data."
            )
        },
        "n_cells_total": int(len(cells)),
        "n_cells_mouse": n_mouse,
        "n_cells_human": n_human,
        "n_donors_mouse": int(n_donors_mouse),
        "n_donors_human": int(n_donors_human),
        "n_cells_cross_species_arm": int(len(cross_species_table)),
        "subclass_counts": subclass_counts,
        "per_species_mean_ephys_features": per_feature_species_means,
        "donor_grouped_split_check": {
            "ran_on_real_donor_ids": True,
            "no_leakage_confirmed": donor_leakage_ok,
            "error": donor_leakage_error,
            "n_train_donors": n_train_donors,
            "n_test_donors": n_test_donors,
        },
        "note": (
            "This summary describes a REAL-DATA run: mouse ephys+metadata from the "
            "Allen Institute Cell Types Database REST API (Gouwens 2020 Cre-line "
            "overlapping cohort), human ephys from a real downloaded subset of "
            "DANDI:000636 (CC-BY-4.0), features extracted with ipfx. No synthetic "
            "or substituted data is included. This is an ephys-ONLY phenotypic "
            "comparison across species; it makes no causal or gene-attribution "
            "claims, and the mouse gene-attribution arm was not run (see "
            "arms_blocked above) due to lack of accessible, properly licensed "
            "linked transcriptome data."
        ),
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return summary


def main() -> int:
    args = parse_args()
    run(
        args.mouse_cells,
        args.mouse_ephys,
        args.human_cells,
        args.human_ephys,
        args.output,
        args.random_state,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
