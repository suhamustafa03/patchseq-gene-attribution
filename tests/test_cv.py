import numpy as np
import pytest

from patchseq_genes import cv


def _groups(n_donors=6, per_donor=5):
    groups = []
    for d in range(n_donors):
        groups.extend([f"donor_{d}"] * per_donor)
    return np.array(groups)


def test_grouped_kfold_never_leaks_donor():
    groups = _groups()
    splits = cv.donor_grouped_kfold_splits(groups, n_splits=3)
    for train_idx, test_idx in splits:
        cv.assert_no_donor_leakage(groups, train_idx, test_idx)


def test_grouped_train_test_split_never_leaks_donor():
    groups = _groups()
    train_idx, test_idx = cv.donor_grouped_train_test_split(groups, test_size=0.3, random_state=1)
    cv.assert_no_donor_leakage(groups, train_idx, test_idx)
    assert len(train_idx) + len(test_idx) == len(groups)


def test_assert_no_donor_leakage_detects_violation():
    groups = _groups()
    train_idx = np.array([0, 1, 2])  # donor_0
    test_idx = np.array([2, 3])  # overlaps donor_0 via index 2
    with pytest.raises(AssertionError):
        cv.assert_no_donor_leakage(groups, train_idx, test_idx)


def test_too_few_donors_raises():
    groups = np.array(["donor_0"] * 10)
    with pytest.raises(ValueError):
        cv.donor_grouped_kfold_splits(groups, n_splits=3)
