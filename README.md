# Candidate-Gene Attribution for Mouse-vs-Human Interneuron Excitability Differences

**Status: IN PROGRESS.** This repository implements the data pipeline, modeling,
attribution, and plausibility-screening code for the project described in
[`docs/stage2_research_design.md`](docs/stage2_research_design.md). No results
in this repository should be read as final; see
[`docs/PROGRESS_REPORT.md`](docs/PROGRESS_REPORT.md) for exactly what has and
has not been run on real data, and why.

## Scope (conditional-narrow, approved design)

- **Mouse-only candidate gene attribution** using the Gouwens et al. 2020
  (*Cell*) Patch-seq cohort (ephys + transcriptome + partial morphology,
  linked per-cell; NeMO `dat-7uhwm74`/`dat-m2kd56v`, parent collection
  licensed CC-BY-4.0).
- **Human electrophysiology-only phenotypic comparison** using DANDI 000636
  ("Human interneuron patch-seq electrophysiology," CC-BY-4.0, confirmed via
  raw DANDI API metadata, 109 subjects / 616 cells, ephys only — no
  transcriptome, no confirmed subclass labels).
- **No cross-species gene attribution is claimed.** The project does not use
  human transcriptomic data. A linked human Patch-seq dataset (Lee, Dalley,
  et al. 2023, *Science* 382:eadf6484, n=778) exists in principle, but its
  exact data accession and data-specific license have not been confirmed as
  of this writing — see `docs/PROGRESS_REPORT.md` §4. If and when that is
  resolved, the scope may be upgraded; until then, every human-side claim in
  this repo and its website is ephys-only.

## What this project does NOT claim

- No causal gene-to-channel-to-excitability proof from observational data.
- No cross-species gene attribution (no linked human transcriptome in scope).
- No simulation result is treated as validation — conductance-based models
  are a plausibility/sensitivity check only.
- Null, inconclusive, or unstable results are reported, not hidden.

## Repository layout

```
data/            Download scripts + data-use/license notes (no raw data committed)
src/patchseq_genes/   Modular pipeline: qc, features, attribution, stability, simulation
tests/           Schema, leakage-prevention, and pipeline tests (run on demo fixture)
demo_fixture/    Small, clearly synthetic fixture for testing code without downloading real data
website/         Static site (built only from data actually present in results/)
docs/            Research design, progress report, references
```

## Setup

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

See `data/README.md` for how to fetch the real datasets (not done automatically —
requires tens of GB of download and is left to the user/compute environment,
per the data-use terms discussed in `docs/stage2_research_design.md`).

## Running on the demo fixture (works now, no download needed)

```bash
pytest
python -m patchseq_genes.run_pipeline --input demo_fixture --output results/demo
```

This exercises the full pipeline (QC -> features -> baseline model ->
stability selection -> simulation plausibility check) end-to-end on
synthetic data, and is what the test suite checks. It is NOT a scientific
result — see `demo_fixture/README.md`.
