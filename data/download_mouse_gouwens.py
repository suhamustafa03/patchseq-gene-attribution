"""Download mouse Patch-seq ephys + metadata from the Gouwens et al. 2020 cohort.

Source: Allen Institute Cell Types Database (ephys/metadata) and the NeMO
Archive accession dat-m2kd56v (https://assets.nemoarchive.org/dat-m2kd56v)
for linked transcriptome counts. Parent NeMO collection is licensed
CC-BY-4.0; see data/README.md for the exact license text and citation
(Gouwens et al. 2020, Cell, doi:10.1016/j.cell.2020.09.057).

This script is NOT executed as part of building this repository (no
network data pulls of GBs of files here) -- it is provided so a user with
their own compute/storage can actually fetch the real cohort.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fetch mouse Patch-seq ephys + metadata (Gouwens et al. 2020) "
            "via the AllenSDK Cell Types API."
        )
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/raw/mouse_gouwens"),
        help="Directory to write downloaded ephys/metadata files into.",
    )
    parser.add_argument(
        "--transcriptome-dir",
        type=Path,
        default=None,
        help=(
            "Optional path to a manually-downloaded NeMO dat-m2kd56v "
            "transcriptome bundle to link against (not fetched by this "
            "script; see data/README.md)."
        ),
    )
    parser.add_argument(
        "--manifest-file",
        type=Path,
        default=None,
        help="Where AllenSDK should cache its manifest JSON (defaults "
        "to <output>/manifest.json).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    try:
        from allensdk.core.cell_types_cache import CellTypesCache
    except ImportError:
        print(
            "ERROR: 'allensdk' is not installed. Install it with:\n"
            "    pip install allensdk\n"
            "then re-run this script. AllenSDK provides CellTypesCache, "
            "which wraps the Cell Types API used to fetch ephys sweeps "
            "and per-cell metadata for the Gouwens et al. 2020 cohort.",
            file=sys.stderr,
        )
        return 1

    args.output.mkdir(parents=True, exist_ok=True)
    manifest_file = args.manifest_file or (args.output / "manifest.json")

    ctc = CellTypesCache(manifest_file=str(manifest_file))

    # species=['Mus musculus'] restricts to the mouse cohort; the Patch-seq
    # subset is identified downstream by cross-referencing transcriptome
    # linkage, which AllenSDK's ephys-only metadata call does not carry.
    cells = ctc.get_cells(species=["Mus musculus"])
    print(f"Found {len(cells)} mouse specimens via CellTypesApi.")

    for cell in cells:
        specimen_id = cell["id"]
        try:
            ctc.get_ephys_data(specimen_id)
        except Exception as exc:  # pragma: no cover - network-dependent
            print(f"  skip {specimen_id}: {exc}", file=sys.stderr)
            continue

    print(f"Ephys NWB files and feature tables cached under {args.output}")

    if args.transcriptome_dir is not None:
        if not args.transcriptome_dir.exists():
            print(
                f"NOTE: --transcriptome-dir {args.transcriptome_dir} does "
                "not exist yet. Download the NeMO dat-m2kd56v bundle "
                "manually (see data/README.md) and place it there to "
                "link transcriptome counts to these ephys records.",
                file=sys.stderr,
            )
    else:
        print(
            "NOTE: linked transcriptome counts are NOT fetched by this "
            "script. See data/README.md for the manual NeMO "
            "dat-m2kd56v download step required to build the full "
            "ephys+transcriptome linked dataset.",
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
