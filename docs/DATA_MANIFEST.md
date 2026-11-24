# Data Manifest — Real-Data Run (2026-11-23)

This documents the first real (non-synthetic, non-demo) data pull and
processing run for this project. All raw files live under `data/raw/`
(gitignored, not committed). This supersedes nothing in `demo_fixture/`
or `results/demo/`, which remain the labeled synthetic pipeline-validation
baseline.

## Step 0 — License / version recheck (performed immediately before download)

- **DANDI:000636** — fetched `https://api.dandiarchive.org/api/dandisets/000636/versions/`
  at **2026-11-23T16:49:11Z**. Two versions exist: `draft` and the published
  `0.241120.0510` (status `Valid`, created 2024-11-20). Full metadata fetched
  from `https://api.dandiarchive.org/api/dandisets/000636/versions/0.241120.0510/`:
  - `license`: `["spdx:CC-BY-4.0"]` — **confirmed open**, unchanged from prior audit.
  - `assetsSummary.numberOfSubjects`: 109
  - `assetsSummary.numberOfCells`: 616
  - `assetsSummary.numberOfFiles` / asset_count: 706
  - `assetsSummary.numberOfBytes`: 24,473,831,989 (~24.5 GB)
  - name: "Human interneuron patch-seq electrophysiology"
- **NeMO `dat-7uhwm74`** (Gouwens 2020 mouse ephys/morphology) — re-fetched
  `https://assets.nemoarchive.org/dat-7uhwm74` at the same time: page shows
  **License: CC BY 4.0** — confirmed open, matches prior audit.
- **NeMO `dat-m2kd56v`** (Gouwens 2020 mouse transcriptome) — re-fetched
  `https://assets.nemoarchive.org/dat-m2kd56v`: page now shows
  **"License: No license listed"** — this is a **change from / is not**
  CC-BY-4.0 or equivalent. Per the stated stop condition, **this source was
  NOT downloaded or used.** This is an additional, independent reason (on
  top of the prior audit's unresolved cell-linkage question) that the real
  mouse transcriptome could not be obtained for this run.

**Conclusion of Step 0:** DANDI:000636 and NeMO `dat-7uhwm74` (mouse
ephys/metadata) are clear to use. NeMO `dat-m2kd56v` (mouse transcriptome)
is blocked by its current license status — the mouse gene-attribution arm
cannot be run on real data as a direct consequence.

## Step 1 — Tooling

- `allensdk`: **failed to install**, confirming the prior audit. Root cause
  reproduced precisely: Python 3.13's `pkgutil` no longer exposes
  `ImpImporter`, which `allensdk`'s very old pinned `setuptools`/`pkg_resources`
  build chain requires; pinning `setuptools<65` just moves the failure one
  step later, into `matplotlib`'s `versioneer.py`, which calls
  `configparser.SafeConfigParser` (removed in Python 3.12+). `allensdk` is
  not installable in this Python 3.13.11 environment by any combination of
  setuptools pins tried in the time available.
- **Used instead**: the Allen Brain Map REST API directly via `requests`
  (`https://api.brain-map.org/api/v2/data/query.json`, models
  `ApiCellTypesSpecimenDetail` and `EphysFeature`) — this is the exact same
  API AllenSDK's `CellTypesCache`/`CellTypesApi` wrap internally, so the
  data obtained is identical in content and provenance to what AllenSDK
  would have returned; only the install method differs. Documented per
  instructions so nobody thinks the dataset was silently swapped.
- `dandi` CLI (`dandi` Python package, v0.83.0): installed and worked
  without issue.
- `pynwb` (3.1.2) and `ipfx` (2.1.2): installed successfully (this pulled in
  a compatible `hdmf`/`pynwb`/`pandas` stack; no Python-3.13 blockers for
  these packages). Used for NWB reading and real ephys feature extraction.

## Step 2 — Real data downloaded

### Human: DANDI:000636 (ephys only)

- Command used: asset-by-asset `dandi download
  "https://api.dandiarchive.org/api/dandisets/000636/versions/0.241120.0510/assets/?path=<path>"
  -o data/raw/human_dandi000636 --existing overwrite`, looped over 20
  pre-selected asset paths (one real, enumerated via
  `DandiAPIClient.get_dandiset('000636','0.241120.0510').get_assets()`
  listing all 706 real asset paths/sizes/ids, then picking one file per
  distinct subject, excluding files over 70 MB, to keep the subset modest
  and time-bounded).
- Download window: 2026-11-23T16:50:52Z -> 2026-11-23T16:58:xx (~7-8 min
  wall clock), done over 2 foreground+background shell calls (one became a
  background task after a 120s tool timeout; it continued and completed
  unattended).
- **20 of 706 NWB files downloaded** (~2.8% of the dandiset), **597 MB**
  total (vs. the dandiset's full 24.5 GB) — a genuine, honestly-labeled
  subset, not the full cohort, chosen for time/scope reasons as explicitly
  permitted.
- 20 distinct subjects (`sub-*` directories), i.e. one session file per
  subject.
- SHA-256 checksums of all 20 downloaded files recorded in
  `data/raw/human_dandi000636/sha256sums.txt` (DANDI's own published
  per-asset digests were not separately cross-checked in this run due to
  time; the dandi CLI's own integrity check on download did not report any
  failures).

### Mouse: Allen Cell Types Database (Gouwens 2020-overlapping cohort), ephys + metadata only

- Source: `https://api.brain-map.org/api/v2/data/query.json`,
  `model::ApiCellTypesSpecimenDetail` (per-cell metadata: donor id, Cre
  line, dendrite type, structure) and `model::EphysFeature` (per-cell
  ephys feature summary: rheobase current, input resistance, spike
  timing/width, adaptation), filtered to `donor__species = 'Mus musculus'`
  and Cre line name containing one of `Pvalb`, `Sst`, `Vip` (the subclass
  markers this project's schema recognizes; `Lamp5` and `Sncg` returned 0
  matches by Cre-line name in this general Cell Types Database query — see
  limitation note below).
- Access date: 2026-11-23, ~16:52-16:56 UTC.
- **522 real mouse cells** with real per-cell ephys features and metadata
  (**263 distinct donor IDs**), written to
  `data/raw/mouse_gouwens/cell_types_mouse_ephys_metadata.json` and
  `data/raw/mouse_gouwens/ephys_feature_detail.json`.
- **Important scope limitation, stated honestly**: this is the Allen Cell
  Types Database's general query interface (the same endpoint AllenSDK's
  `CellTypesApi` wraps), filtered by Cre-line name as a GABAergic-subclass
  proxy. It was **not** possible in the time available to confirm that
  every one of these 522 cells is specifically part of the Gouwens et al.
  2020 Patch-seq-linked release (vs. the Cell Types Database's broader,
  overlapping electrophysiology-only collection); no linked transcriptome
  was fetched for any of them (blocked per Step 0). They are real Allen
  Institute mouse interneuron ephys recordings with real Cre-line-based
  subclass labels and donor IDs, used here for the ephys-only comparison
  arm, not claimed as a verified 1:1 match to the published Gouwens 2020
  Patch-seq sample list.

## Step 3 — Real feature extraction + QC

- Human features extracted from the 20 downloaded NWB files with `ipfx`
  (`ipfx.dataset.create.create_ephys_data_set` +
  `ipfx.data_set_features.extract_data_set_features`, using default
  detection parameters — no `allensdk` dependency was needed, confirming
  the Step-1 fallback worked end to end).
  - **4 of 20 files yielded usable `long_squares` cell features; 16 of 20
    failed** (`no long_squares features extracted`, or an unpacking
    error), because ipfx's default stimulus-name ontology did not
    recognize the long-square stimulus-code naming used in most of these
    particular NWB files (seen directly in the sweep table: stimulus
    names like `Unknown`, `Noise 1`, `Noise 2` that ipfx's default
    ontology does not map to its `LONG_SQUARE` stimulus type). This is a
    genuine real-data extraction limitation, not a bug suppressed to force
    a result — fixing it for the full dataset would require per-site
    stimulus-ontology customization beyond the time available for this
    run. It is reported here rather than hidden.
  - For the 4 successful extractions: `rheobase_pA` <-
    `cell_features['long_squares']['rheobase_i']`; `input_resistance_MOhm`
    <- `...['input_resistance']`; `ap_width_ms` <- the first spike's
    `width` (seconds) on the rheobase sweep, ×1000; `adaptation_ratio` <-
    the `adapt` field of ipfx's "hero sweep" (its standard ~1.5×-rheobase
    sweep used for adaptation, since the rheobase sweep itself often fires
    too few spikes for the ISI-ratio-based adaptation metric to be
    defined — a real, physiological NaN, not a bug; this is ipfx's own
    convention, not an invented fallback).
- Mouse features taken directly from the Allen `EphysFeature` records:
  `rheobase_pA` <- `threshold_i_long_square`; `input_resistance_MOhm` <-
  `input_resistance_mohm`; `adaptation_ratio` <- `adaptation`;
  `ap_width_ms` <- `(fast_trough_t_short_square - threshold_t_short_square)
  * 1000` (a documented proxy for a canonical half-width feature, since
  this particular REST endpoint does not expose a ready-made spike-width
  field the way the per-sweep `ipfx` spike table does for the human data).
- `src/patchseq_genes/qc.py` schema validation and QC filters were run
  unmodified except for one minimal, documented change: `VALID_SUBCLASSES`
  now also accepts `"unknown"`, because DANDI:000636 carries no confirmed
  transcriptomic subclass label for its cells (the demo fixture's
  synthetic subclass set Pvalb/Sst/Vip/Lamp5/Sncg is untouched and
  unaffected).
- `qc.run_qc` on the combined real cells+ephys table: **526 cells in ->
  526 cells out** (schema valid); ephys features: **526 rows in -> 454
  rows out** (72 dropped for missing/non-finite key ephys features —
  these are rows where one or more of rheobase/width/Ri/adaptation came
  back null from the source API, a real QC exclusion, not induced).

### Real usable counts AFTER QC (exact numbers)

| | mouse | human | total |
|---|---|---|---|
| Cells (schema-valid) | 522 | 4 | 526 |
| Cells with complete ephys features (post-QC) | 450 | 4 | 454 |
| Distinct donors (post-QC cross-species table) | — | — | 186 train + 62 test (see below) |
| Subclasses present | Pvalb (201), Sst (166), Vip (155) | unknown (4) | — |

## Step 4 — Pipeline run

- **Mouse gene-attribution arm: BLOCKED.** No real transcriptome input
  exists (Step 0/2), so `build_attribution_table()`, `stability_selection`,
  and all three attribution models (ridge/lasso/reduced-rank) in
  `run_pipeline.py` were never invoked on real data. This is reported as a
  blocker, not worked around with synthetic/substitute expression data.
- **Human ephys-only cross-species comparison arm: RAN on real data**, via
  the new `src/patchseq_genes/real_data_adapter.py`, which loads the real
  mouse + human cells/ephys tables above, runs them through
  `qc.run_qc`, then `features.build_cross_species_ephys_table`
  (mouse+human, ephys-only, no expression columns — same function used by
  the demo pipeline, unmodified) and `confounders.build_confounder_matrix`.
  Output: `results/real/cross_species_ephys_table.csv`,
  `results/real/cross_species_confounders.csv`, `results/real/summary.json`.
- **Donor-grouped CV leakage check**: ran `cv.donor_grouped_train_test_split`
  + `cv.assert_no_donor_leakage` (both unmodified, from `src/patchseq_genes/cv.py`)
  on the real cross-species table's real donor IDs (`mouse_donor_<Allen donor id>`,
  `human_donor_<DANDI subject id>`). Result: **no donor leakage** across a
  186-real-donor train split and 62-real-donor test split (confirmed by
  assertion, not just inspection) — recorded in
  `results/real/summary.json` under `donor_grouped_split_check`.
- Stability selection / gene attribution was **not run** (no real
  expression matrix to select from); this is stated explicitly in
  `results/real/summary.json.arms_blocked`.

## Total wall-clock time on downloads

Roughly 10-12 minutes: ~7-8 minutes for the 20-file, 597 MB DANDI:000636
download (one call backgrounded after a 120s foreground tool timeout, then
completed unattended), plus a few minutes of Allen Brain Map REST API calls
(paginated metadata + per-cell ephys feature queries for 522 cells, done in
batches of 25-50 to stay within reasonable request sizes).
