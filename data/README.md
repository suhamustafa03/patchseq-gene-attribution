# Data

No raw data is committed to this repository. `data/raw/` and
`data/processed/` are listed in `.gitignore`. This directory only holds
download scripts and this license/provenance documentation.

## 1. Mouse: Gouwens et al. 2020 Patch-seq cohort (Allen Institute)

- **Source**: NeMO Archive, accessions `dat-7uhwm74` (morphology/ephys) and
  `dat-m2kd56v` (https://assets.nemoarchive.org/dat-m2kd56v); also
  accessible through the Allen Institute Cell Types Database via
  `allensdk.api.queries.cell_types_api.CellTypesApi` /
  `allensdk.core.cell_types_cache.CellTypesCache` for ephys sweeps and
  per-cell metadata.
- **License**: parent NeMO collection is licensed **CC-BY-4.0**. Quote
  from the NeMO Archive terms: "Data in the NeMO Archive are made
  available under the Creative Commons Attribution 4.0 International
  (CC-BY-4.0) license unless otherwise noted on a per-dataset basis."
  Re-confirm the per-dataset license tag on the NeMO page before use, as
  license metadata is sometimes set at the sub-collection level.
- **Citation**: Gouwens NW, Sorensen SA, Baftizadeh F, et al. (2020)
  "Integrated Morphoelectric and Transcriptomic Classification of Cortical
  GABAergic Cells." *Cell* 183(4):935-953.e19.
  doi:10.1016/j.cell.2020.09.057.
- **What `download_mouse_gouwens.py` does**: uses `allensdk`'s
  `CellTypesApi` to list Patch-seq-eligible specimens, pull per-cell ephys
  feature summaries and sweep data, and write them to
  `data/raw/mouse_gouwens/`. It does NOT fetch the linked transcriptome
  counts (those are distributed via NeMO separately); a second,
  documented manual step (download the NeMO `dat-m2kd56v` bundle and
  point `--transcriptome-dir` at it) is required for the full linked
  dataset. This script has not been run as part of this repository build.
- **Real download commands** (not executed here):
  ```bash
  pip install allensdk
  python data/download_mouse_gouwens.py --output data/raw/mouse_gouwens
  # Separately, for linked transcriptome counts:
  #   fetch https://assets.nemoarchive.org/dat-m2kd56v via the NeMO Archive
  #   browser/CLI and place the extracted matrix under
  #   data/raw/mouse_gouwens/transcriptome/
  ```

## 2. Human: DANDI:000636 (human interneuron Patch-seq electrophysiology)

- **Source**: DANDI Archive, dandiset `000636`
  (https://dandiarchive.org/dandiset/000636).
- **License**: **CC-BY-4.0**, confirmed by querying the raw DANDI API
  metadata for this dandiset (`license` field reports
  `spdx:CC-BY-4.0`). Ephys only (NWB `.nwb` files with intracellular
  electrophysiology); there is no linked transcriptome and no confirmed
  cell-subclass label in this dandiset as distributed.
- **Scale**: 109 subjects / 616 cells (per dandiset metadata as of
  project scoping; re-check the live dandiset page, which may be updated).
- **What `download_human_dandi000636.py` does**: uses the `dandi` Python
  client (`dandi.download.download`) to pull the dandiset's NWB files into
  `data/raw/human_dandi000636/`. This is a large (multi-GB) download and
  is NOT executed as part of building this repository.
- **Real download command** (not executed here):
  ```bash
  pip install dandi
  dandi download DANDI:000636 -o data/raw/human_dandi000636
  # or, via the wrapper script:
  python data/download_human_dandi000636.py --output data/raw/human_dandi000636
  ```

## Not in scope (yet)

**Lee, Dalley, et al. 2023** (*Science* 382:eadf6484, n=778 human
Patch-seq neurons with linked transcriptome) would be the natural upgrade
path to a true cross-species gene-attribution arm. Its exact data
accession and data-specific license have not been confirmed as of this
writing, so it is **not used** anywhere in this codebase or website. See
`docs/PROGRESS_REPORT.md` for the open items blocking that upgrade.
