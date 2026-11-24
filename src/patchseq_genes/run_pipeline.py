"""CLI entry point: QC -> features -> baseline attribution -> stability
selection (with negative control) -> optional simulation plausibility
screen. Writes JSON/CSV results and a summary.json.

Usage:
    python -m patchseq_genes.run_pipeline --input demo_fixture --output results/demo
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from . import qc
from .attribution import lasso_attribution, reduced_rank_attribution, ridge_attribution
from .confounders import build_confounder_matrix, donor_groups
from .cv import donor_grouped_train_test_split, assert_no_donor_leakage
from .features import build_attribution_table, build_cross_species_ephys_table, get_gene_columns, load_tables
from .stability import stability_selection, stable_genes

# Pre-specified (before seeing any results) minimum count of stable genes
# required to call a target's attribution result anything other than
# null/inconclusive. Chosen once here, not tuned post hoc.
STABLE_GENE_MIN_COUNT = 2

ATTRIBUTION_TARGETS = ["rheobase_pA", "ap_width_ms"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--run-simulation",
        action="store_true",
        help="Also run the Brian2 biophysical plausibility screen (requires brian2).",
    )
    parser.add_argument("--n-resamples", type=int, default=100)
    parser.add_argument("--random-state", type=int, default=0)
    return parser.parse_args()


def run(input_dir: Path, output_dir: Path, run_simulation: bool, n_resamples: int, random_state: int) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)

    cells, ephys, expression = load_tables(input_dir)
    cells, ephys_clean, qc_reports = qc.run_qc(cells, ephys)
    for report in qc_reports:
        print(f"[qc] {report.table_name}: in={report.n_rows_in} out={report.n_rows_out} "
              f"dropped_missing={report.dropped_missing_key_features}")

    attribution_table = build_attribution_table(cells, ephys_clean, expression)
    cross_species_table = build_cross_species_ephys_table(cells, ephys_clean)
    gene_columns = get_gene_columns(expression)

    attribution_table.to_csv(output_dir / "attribution_table.csv", index=False)
    cross_species_table.to_csv(output_dir / "cross_species_ephys_table.csv", index=False)

    # Sanity-check donor-grouped CV leakage once per run, on the
    # attribution table, so a regression in cv.py fails loudly here too.
    groups = donor_groups(attribution_table)
    train_idx, test_idx = donor_grouped_train_test_split(groups, random_state=random_state)
    assert_no_donor_leakage(groups, train_idx, test_idx)

    confounders = build_confounder_matrix(attribution_table)
    X_genes = attribution_table[gene_columns].astype(float)
    X_full = pd.concat([X_genes, confounders], axis=1)

    attribution_results = []
    stability_by_target = {}
    for target in ATTRIBUTION_TARGETS:
        y = attribution_table[target]
        ridge_res = ridge_attribution(X_full, y, groups, target)
        lasso_res = lasso_attribution(X_full, y, groups, target)
        rrr_res = reduced_rank_attribution(X_full, y, groups, target)
        attribution_results.extend([ridge_res, lasso_res, rrr_res])

        stab = stability_selection(X_genes, y, n_resamples=n_resamples, random_state=random_state)
        stable = stable_genes(stab)
        stability_by_target[target] = {
            "selection_frequency": stab.selection_frequency.round(4).to_dict(),
            "null_selection_frequency": stab.null_selection_frequency.round(4).to_dict(),
            "stable_genes": stable,
            "n_stable_genes": len(stable),
            "is_null_or_inconclusive": len(stable) < STABLE_GENE_MIN_COUNT,
        }

    attribution_rows = [
        {
            "model": r.model_name,
            "target": r.target,
            "mean_r2": round(r.mean_r2, 4),
            "fold_r2": [round(v, 4) for v in r.fold_r2],
            "n_folds": r.n_folds,
        }
        for r in attribution_results
    ]
    (output_dir / "attribution_results.json").write_text(json.dumps(attribution_rows, indent=2))
    (output_dir / "stability_results.json").write_text(json.dumps(stability_by_target, indent=2))

    simulation_summary = None
    if run_simulation:
        from .simulation import GENE_CHANNEL_MAP, run_plausibility_screen

        sim_results = []
        for gene in GENE_CHANNEL_MAP:
            result = run_plausibility_screen(gene, conductance_multiplier=1.5)
            sim_results.append(result.__dict__)
        (output_dir / "simulation_results.json").write_text(json.dumps(sim_results, indent=2))
        simulation_summary = {"n_genes_screened": len(sim_results)}

    n_stable_total = sum(v["n_stable_genes"] for v in stability_by_target.values())
    any_null = any(v["is_null_or_inconclusive"] for v in stability_by_target.values())

    readme_path = input_dir / "README.md"
    is_synthetic = "demo_fixture" in str(input_dir).lower() or (
        readme_path.exists() and "SYNTHETIC" in readme_path.read_text()
    )

    summary = {
        "status": "demo_fixture_only" if is_synthetic else "unknown_input",
        "input_dir": str(input_dir),
        "n_cells_total": int(len(cells)),
        "n_cells_mouse": int((cells["species"] == "mouse").sum()),
        "n_cells_human": int((cells["species"] == "human").sum()),
        "n_cells_attribution_arm": int(len(attribution_table)),
        "n_cells_cross_species_arm": int(len(cross_species_table)),
        "n_genes": len(gene_columns),
        "attribution_targets": ATTRIBUTION_TARGETS,
        "stable_gene_min_count_threshold": STABLE_GENE_MIN_COUNT,
        "n_stable_genes_total": n_stable_total,
        "result_is_null_or_inconclusive_for_any_target": any_null,
        "stability_by_target": {
            t: {"n_stable_genes": v["n_stable_genes"], "stable_genes": v["stable_genes"],
                "is_null_or_inconclusive": v["is_null_or_inconclusive"]}
            for t, v in stability_by_target.items()
        },
        "simulation": simulation_summary,
        "note": (
            "This summary describes whatever data was in --input. If that "
            "was demo_fixture/, these numbers are 100% synthetic and "
            "validate the pipeline code only, not biology."
        ),
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return summary


def main() -> int:
    args = parse_args()
    run(args.input, args.output, args.run_simulation, args.n_resamples, args.random_state)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
