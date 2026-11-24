"""Baseline regression models mapping gene expression -> an ephys feature.

Mouse-only, since only mouse cells carry a linked transcriptome in this
project's scope.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.linear_model import Lasso, LinearRegression, Ridge
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import StandardScaler

from .cv import donor_grouped_kfold_splits


@dataclass
class AttributionResult:
    model_name: str
    target: str
    mean_r2: float
    fold_r2: list[float]
    n_folds: int


def _cv_r2(model_factory, X: np.ndarray, y: np.ndarray, groups: np.ndarray, n_splits: int) -> list[float]:
    scores = []
    for train_idx, test_idx in donor_grouped_kfold_splits(groups, n_splits=n_splits):
        model = model_factory()
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X[train_idx])
        X_test = scaler.transform(X[test_idx])
        model.fit(X_train, y[train_idx])
        scores.append(model.score(X_test, y[test_idx]))
    return scores


def ridge_attribution(
    X: pd.DataFrame, y: pd.Series, groups: np.ndarray, target: str, alpha: float = 1.0, n_splits: int = 5
) -> AttributionResult:
    fold_r2 = _cv_r2(lambda: Ridge(alpha=alpha), X.to_numpy(dtype=float), y.to_numpy(dtype=float), groups, n_splits)
    return AttributionResult("ridge", target, float(np.mean(fold_r2)), fold_r2, len(fold_r2))


def lasso_attribution(
    X: pd.DataFrame, y: pd.Series, groups: np.ndarray, target: str, alpha: float = 0.1, n_splits: int = 5
) -> AttributionResult:
    fold_r2 = _cv_r2(
        lambda: Lasso(alpha=alpha, max_iter=10000), X.to_numpy(dtype=float), y.to_numpy(dtype=float), groups, n_splits
    )
    return AttributionResult("lasso", target, float(np.mean(fold_r2)), fold_r2, len(fold_r2))


def reduced_rank_attribution(
    X: pd.DataFrame, y: pd.Series, groups: np.ndarray, target: str, n_components: int = 5, n_splits: int = 5
) -> AttributionResult:
    """Truncated-SVD of the expression matrix + linear regression on the
    top components, used as a reduced-rank baseline when no dedicated
    reduced-rank-regression library is available.
    """
    n_components = min(n_components, X.shape[1] - 1, X.shape[0] - 1)
    n_components = max(1, n_components)

    def factory():
        return _SVDLinearRegression(n_components)

    fold_r2 = _cv_r2(factory, X.to_numpy(dtype=float), y.to_numpy(dtype=float), groups, n_splits)
    return AttributionResult("reduced_rank_svd", target, float(np.mean(fold_r2)), fold_r2, len(fold_r2))


class _SVDLinearRegression:
    def __init__(self, n_components: int) -> None:
        self.svd = TruncatedSVD(n_components=n_components, random_state=0)
        self.reg = LinearRegression()

    def fit(self, X, y):
        Z = self.svd.fit_transform(X)
        self.reg.fit(Z, y)
        return self

    def predict(self, X):
        Z = self.svd.transform(X)
        return self.reg.predict(Z)

    def score(self, X, y):
        from sklearn.metrics import r2_score

        return r2_score(y, self.predict(X))
