# Stage 2 Research Design: Candidate-Gene Attribution with Biophysical Plausibility Screening

## Status

IN PROGRESS. This document describes the approved "conditional-narrow" scope
and the analysis plan. No real-dataset results exist yet; see
`docs/PROGRESS_REPORT.md`.

## Research question

Among genes expressed in mouse Patch-seq interneurons, which show the
strongest statistical association with electrophysiological phenotype
(after controlling for subclass and donor), and are those associations
biophysically plausible given known ion-channel gene function? This is
framed throughout as **candidate-gene attribution with biophysical
plausibility screening**, never as a causal claim.

## Datasets in scope

1. **Gouwens et al. 2020** (*Cell*), mouse Patch-seq cohort, n=4,270
   interneurons with linked ephys + transcriptome + partial morphology.
   NeMO accessions `dat-7uhwm74` / `dat-m2kd56v`; parent collection
   CC-BY-4.0.
2. **DANDI:000636**, human interneuron Patch-seq electrophysiology
   (ephys only, no transcriptome, no confirmed subclass labels), 109
   subjects / 616 cells, CC-BY-4.0 (confirmed via raw DANDI API metadata).

No cross-species gene attribution is claimed: there is no linked human
transcriptome in scope. The human arm is restricted to an ephys-only
cross-species phenotypic comparison.

## Workflow

1. **QC and schema validation** (`qc.py`) — enforce required columns,
   types, no duplicate cell IDs, and drop rows with missing key features
   before any modeling, so silent NaNs never reach regressions.
2. **Feature construction** (`features.py`) — two analysis arms are kept
   strictly separate: (a) the attribution arm restricts to mouse cells
   with `has_transcriptome=True` and uses gene expression + ephys; (b) the
   ephys-only cross-species comparison arm uses all cells (both species)
   but never touches expression data, since human transcriptome is out of
   scope.
3. **Confounder control** (`confounders.py`) — subclass is one-hot encoded
   and included as a covariate; donor identity is retained as a grouping
   key (never as a model feature) so that cross-validation can respect it.
4. **Donor-grouped cross-validation** (`cv.py`) — splits are grouped by
   `donor_id`, not `cell_id`, because multiple cells from the same donor
   share donor-level technical and biological variation; a cell-level
   split would leak donor identity between train and test and inflate
   apparent predictability.
5. **Baseline attribution models** (`attribution.py`) — ridge and lasso
   regression mapping gene expression to a target ephys feature, plus a
   reduced-rank baseline (truncated SVD of the expression matrix followed
   by linear regression on the top components).
6. **Stability selection** (`stability.py`) — repeated subsampling + lasso
   across random subsets, reporting per-gene selection frequency. A
   negative control repeats the identical procedure on label-shuffled
   data; if real selection frequencies are not well above the null
   (shuffled) distribution, the result is reported as null/inconclusive,
   not hidden.
7. **Simulation plausibility screen** (`simulation.py`) — a minimal
   single-compartment conductance-based neuron model (Brian2). Candidate
   genes with a documented, simplified ion-channel mapping are perturbed
   to see whether the direction of the firing-rate / AP-width shift is
   at least qualitatively consistent with the statistical association.
   This is explicitly a plausibility screen, not independent validation.

## Pre-specified null/inconclusive criterion

A gene-level result is called **null/inconclusive** when fewer than a
pre-specified minimum number of genes exceed the negative-control
selection-frequency distribution; this threshold is written once into
`run_pipeline.py` (`STABLE_GENE_MIN_COUNT`) so it cannot be chosen after
seeing results. The pipeline reports this status explicitly in
`summary.json` rather than suppressing it.

## What this project does not claim

- No causal gene -> channel -> excitability claim.
- No cross-species gene attribution (no linked human transcriptome).
- Simulation results are a plausibility screen, not validation.
