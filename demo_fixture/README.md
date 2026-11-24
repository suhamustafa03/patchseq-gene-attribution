# Demo fixture

**THIS IS 100% SYNTHETIC DATA FOR TESTING ONLY. It is NOT real measurements
of any mouse or human cell, and it must NEVER be presented as a scientific
result.**

It exists only to let the pipeline code (`src/patchseq_genes/`) and its
tests run end-to-end without downloading the real, multi-GB Gouwens 2020
and DANDI:000636 datasets. Generate it with:

```bash
python3 demo_fixture/generate_demo_fixture.py
```

This writes:

- `cells.csv` -- `cell_id, species, subclass, donor_id, has_transcriptome`
  for ~300 synthetic mouse cells and ~100 synthetic human cells.
- `ephys_features.csv` -- `cell_id, rheobase_pA, ap_width_ms,
  input_resistance_MOhm, adaptation_ratio`, drawn from random
  distributions with small, arbitrary species/subclass offsets baked in
  so downstream code has some structure to recover.
- `gene_expression.csv` -- `cell_id` + gene-symbol columns of random
  Poisson counts, generated **only** for mouse cells
  (`has_transcriptome=True`), consistent with this project's scope (no
  human transcriptome).

Any numbers that come out of the pipeline when run on this fixture (see
`results/demo/summary.json`) describe the fixture and the code, not
biology.
