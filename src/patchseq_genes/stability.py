"""Stability selection: repeated subsampling + lasso, with a negative
control that shuffles the target label to estimate a null selection-
frequency distribution.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.linear_model import Lasso
from sklearn.preprocessing import StandardScaler


@dataclass
class StabilityResult:
    selection_frequency: pd.Series  # gene -> fraction of resamples selected
    null_selection_frequency: pd.Series  # same, on label-shuffled data
    n_resamples: int
    subsample_fraction: float
    alpha: float


def _one_resample(
    X: np.ndarray, y: np.ndarray, rng: np.random.Generator, subsample_fraction: float, alpha: float
) -> np.ndarray:
    n = X.shape[0]
    n_sub = max(2, int(round(n * subsample_fraction)))
    idx = rng.choice(n, size=n_sub, replace=False)
    scaler = StandardScaler()
    X_sub = scaler.fit_transform(X[idx])
    model = Lasso(alpha=alpha, max_iter=10000)
    model.fit(X_sub, y[idx])
    return model.coef_ != 0


def stability_selection(
    X: pd.DataFrame,
    y: pd.Series,
    n_resamples: int = 100,
    subsample_fraction: float = 0.7,
    alpha: float = 0.1,
    random_state: int = 0,
) -> StabilityResult:
    rng = np.random.default_rng(random_state)
    genes = list(X.columns)
    X_arr = X.to_numpy(dtype=float)
    y_arr = y.to_numpy(dtype=float)

    selected_counts = np.zeros(len(genes), dtype=int)
    for _ in range(n_resamples):
        selected_counts += _one_resample(X_arr, y_arr, rng, subsample_fraction, alpha).astype(int)
    real_freq = pd.Series(selected_counts / n_resamples, index=genes)

    # Negative control: identical procedure, labels permuted. Fresh RNG
    # stream (not reusing `rng`) keeps this reproducible independent of
    # how many draws the real pass above consumed.
    null_rng = np.random.default_rng(random_state + 1)
    y_shuffled = y_arr.copy()
    null_rng.shuffle(y_shuffled)
    null_counts = np.zeros(len(genes), dtype=int)
    for _ in range(n_resamples):
        null_counts += _one_resample(X_arr, y_shuffled, null_rng, subsample_fraction, alpha).astype(int)
    null_freq = pd.Series(null_counts / n_resamples, index=genes)

    return StabilityResult(
        selection_frequency=real_freq,
        null_selection_frequency=null_freq,
        n_resamples=n_resamples,
        subsample_fraction=subsample_fraction,
        alpha=alpha,
    )


def stable_genes(result: StabilityResult, margin: float = 0.2) -> list[str]:
    """Genes whose real selection frequency exceeds the null (shuffled)
    frequency by at least `margin`; this is the pre-specified criterion
    used by run_pipeline to decide null/inconclusive vs. candidate genes.
    """
    diff = result.selection_frequency - result.null_selection_frequency
    return sorted(diff[diff > margin].index.tolist())
