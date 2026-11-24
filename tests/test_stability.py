import numpy as np
import pandas as pd

from patchseq_genes.stability import stability_selection, stable_genes


def _signal_dataset(n=120, n_noise_genes=15, seed=0):
    rng = np.random.default_rng(seed)
    signal_gene = rng.normal(0, 1, n)
    noise = rng.normal(0, 1, (n, n_noise_genes))
    y = 3.0 * signal_gene + rng.normal(0, 0.5, n)

    X = pd.DataFrame(noise, columns=[f"noise_{i}" for i in range(n_noise_genes)])
    X.insert(0, "signal_gene", signal_gene)
    return X, pd.Series(y)


def test_negative_control_null_selection_frequency_is_low():
    X, y = _signal_dataset()
    result = stability_selection(X, y, n_resamples=60, alpha=0.5, random_state=42)
    # Sanity check: on label-shuffled data the "true" effect is destroyed,
    # so the null selection frequency for the real signal gene should not
    # be systematically high -- the method shouldn't find spurious signal
    # by construction.
    assert result.null_selection_frequency["signal_gene"] < 0.5
    assert result.null_selection_frequency.mean() < 0.4


def test_real_signal_gene_is_recovered_more_than_null():
    X, y = _signal_dataset()
    result = stability_selection(X, y, n_resamples=60, alpha=0.2, random_state=42)
    assert (
        result.selection_frequency["signal_gene"]
        > result.null_selection_frequency["signal_gene"]
    )


def test_stable_genes_uses_margin_threshold():
    X, y = _signal_dataset()
    result = stability_selection(X, y, n_resamples=60, alpha=0.2, random_state=42)
    stable = stable_genes(result, margin=0.2)
    assert "signal_gene" in stable
