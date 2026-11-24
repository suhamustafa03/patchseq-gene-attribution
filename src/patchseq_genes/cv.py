"""Donor-grouped, leakage-safe cross-validation split utilities.

Splits are always grouped by donor_id, never by cell_id: multiple cells
from one donor share donor-level biological/technical variation, so a
cell-level split would let the model see a donor's cells in both train
and test and overstate generalization.
"""

from __future__ import annotations

import numpy as np
from sklearn.model_selection import GroupKFold, GroupShuffleSplit


def donor_grouped_kfold_splits(groups: np.ndarray, n_splits: int = 5):
    n_splits = min(n_splits, len(np.unique(groups)))
    if n_splits < 2:
        raise ValueError("Need at least 2 distinct donors for grouped CV.")
    gkf = GroupKFold(n_splits=n_splits)
    dummy_X = np.zeros((len(groups), 1))
    return list(gkf.split(dummy_X, groups=groups))


def donor_grouped_train_test_split(
    groups: np.ndarray, test_size: float = 0.25, random_state: int | None = 0
):
    gss = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=random_state)
    dummy_X = np.zeros((len(groups), 1))
    train_idx, test_idx = next(gss.split(dummy_X, groups=groups))
    return train_idx, test_idx


def assert_no_donor_leakage(groups: np.ndarray, train_idx: np.ndarray, test_idx: np.ndarray) -> None:
    train_donors = set(groups[train_idx])
    test_donors = set(groups[test_idx])
    overlap = train_donors & test_donors
    if overlap:
        raise AssertionError(f"Donor leakage across train/test: {overlap}")
