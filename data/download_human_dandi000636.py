"""Download human Patch-seq electrophysiology from DANDI:000636.

Source: DANDI Archive, dandiset 000636
(https://dandiarchive.org/dandiset/000636). License is CC-BY-4.0,
confirmed via the raw DANDI API dandiset metadata (`license` field reports
`spdx:CC-BY-4.0`). This dandiset is electrophysiology ONLY: there is no
linked transcriptome and no confirmed cell-subclass label, which is why
this project treats the human arm as an ephys-only phenotypic comparison,
not a gene-attribution arm.

This script is NOT executed as part of building this repository (this is
a multi-GB download) -- it documents the real command a user would run.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

DANDISET_ID = "000636"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=f"Fetch DANDI:{DANDISET_ID} (human Patch-seq ephys, CC-BY-4.0) "
        "via the dandi Python client."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/raw/human_dandi000636"),
        help="Directory to download NWB files into.",
    )
    parser.add_argument(
        "--dandiset",
        type=str,
        default=DANDISET_ID,
        help="Dandiset identifier (default: 000636).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    try:
        from dandi.download import download as dandi_download
    except ImportError:
        print(
            "ERROR: 'dandi' is not installed. Install it with:\n"
            "    pip install dandi\n"
            "then re-run this script, or run the CLI directly:\n"
            f"    dandi download DANDI:{args.dandiset} -o {args.output}",
            file=sys.stderr,
        )
        return 1

    args.output.mkdir(parents=True, exist_ok=True)
    url = f"https://dandiarchive.org/dandiset/{args.dandiset}"
    print(f"Downloading {url} -> {args.output} (this is a large, multi-GB download)")

    dandi_download([url], str(args.output))

    print(f"Done. NWB files for DANDI:{args.dandiset} are under {args.output}")
    print(
        "Reminder: this dandiset is ephys-only (no transcriptome, no "
        "confirmed subclass labels) -- do not use it for gene attribution."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
