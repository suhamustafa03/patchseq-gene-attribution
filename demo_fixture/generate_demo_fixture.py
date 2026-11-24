"""Generate a small, 100% synthetic fixture matching the real pipeline's schema.

See demo_fixture/README.md: this is NOT real data and must never be
presented as a scientific result. It exists only so the downstream code
(QC, feature joins, CV, attribution, stability selection, simulation) can
be exercised and tested end-to-end without a multi-GB real download.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

SUBCLASSES = ["Pvalb", "Sst", "Vip", "Lamp5", "Sncg"]

# A few genes with a crude, illustrative association baked in so stability
# selection has *something* non-random to (sometimes) find; this is a code
# test fixture, not a biological claim.
SIGNAL_GENES = ["Kcnc1", "Scn1a", "Pvalb", "Gad1"]
NOISE_GENES = [
    "Sst", "Vip", "Lamp5", "Cck", "Npy", "Calb1", "Calb2", "Reln",
    "Cplx3", "Tac1", "Penk", "Crh", "Nos1", "Chat", "Th", "Gad2",
    "Slc17a6", "Slc17a7", "Htr3a", "Cux2",
]
ALL_GENES = SIGNAL_GENES + NOISE_GENES

N_MOUSE = 300
N_HUMAN = 100
SEED = 20261123


def make_cells(rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    cid = 0
    for species, n_cells, n_donors in (("mouse", N_MOUSE, 30), ("human", N_HUMAN, 20)):
        donor_ids = [f"{species}_donor_{i:03d}" for i in range(n_donors)]
        for _ in range(n_cells):
            cid += 1
            donor = rng.choice(donor_ids)
            subclass = rng.choice(SUBCLASSES, p=[0.30, 0.25, 0.20, 0.15, 0.10])
            rows.append(
                {
                    "cell_id": f"{species}_{cid:05d}",
                    "species": species,
                    "subclass": subclass,
                    "donor_id": donor,
                    "has_transcriptome": species == "mouse",
                }
            )
    return pd.DataFrame(rows)


def make_ephys(cells: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    # Species/subclass offsets are arbitrary but fixed, giving downstream
    # code real (synthetic) signal to recover rather than pure noise.
    species_offset = {"mouse": 0.0, "human": 15.0}
    subclass_offset = {
        "Pvalb": -20.0, "Sst": 10.0, "Vip": 25.0, "Lamp5": 5.0, "Sncg": 0.0,
    }
    rows = []
    for _, cell in cells.iterrows():
        base_rheobase = 80.0 + species_offset[cell["species"]] + subclass_offset[cell["subclass"]]
        rheobase = max(5.0, rng.normal(base_rheobase, 15.0))
        ap_width = max(0.1, rng.normal(0.6 + 0.01 * subclass_offset[cell["subclass"]], 0.12))
        input_resistance = max(10.0, rng.normal(200.0 - species_offset[cell["species"]], 40.0))
        adaptation_ratio = float(np.clip(rng.normal(1.2, 0.3), 0.1, 5.0))
        rows.append(
            {
                "cell_id": cell["cell_id"],
                "rheobase_pA": round(rheobase, 3),
                "ap_width_ms": round(ap_width, 4),
                "input_resistance_MOhm": round(input_resistance, 3),
                "adaptation_ratio": round(adaptation_ratio, 4),
            }
        )
    return pd.DataFrame(rows)


def make_expression(cells: pd.DataFrame, ephys: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    mouse_cells = cells[cells["has_transcriptome"]].merge(ephys, on="cell_id")
    rows = []
    for _, cell in mouse_cells.iterrows():
        counts = {}
        # Kcnc1 counts track rheobase (crude linear + noise) so a baseline
        # ridge/lasso attribution has a real, recoverable (if weak) signal.
        kcnc1 = max(0, rng.poisson(5 + 0.08 * cell["rheobase_pA"]))
        scn1a = max(0, rng.poisson(20 - 0.05 * cell["ap_width_ms"] * 10))
        counts["Kcnc1"] = kcnc1
        counts["Scn1a"] = scn1a
        counts["Pvalb"] = max(0, rng.poisson(30 if cell["subclass"] == "Pvalb" else 5))
        counts["Gad1"] = max(0, rng.poisson(40))
        for gene in NOISE_GENES:
            counts[gene] = max(0, rng.poisson(rng.integers(1, 25)))
        rows.append({"cell_id": cell["cell_id"], **counts})
    return pd.DataFrame(rows)[["cell_id"] + ALL_GENES]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).parent)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    args.output.mkdir(parents=True, exist_ok=True)

    cells = make_cells(rng)
    ephys = make_ephys(cells, rng)
    expression = make_expression(cells, ephys, rng)

    cells.to_csv(args.output / "cells.csv", index=False)
    ephys.to_csv(args.output / "ephys_features.csv", index=False)
    expression.to_csv(args.output / "gene_expression.csv", index=False)

    print(f"Wrote {len(cells)} cells ({(cells.species == 'mouse').sum()} mouse, "
          f"{(cells.species == 'human').sum()} human) to {args.output}")
    print("SYNTHETIC DATA ONLY -- see demo_fixture/README.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
