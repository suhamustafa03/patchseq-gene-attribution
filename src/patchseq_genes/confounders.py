"""Design matrices that separate subclass/species covariates from the
donor grouping key, so donor never leaks into the feature matrix itself.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def one_hot_subclass(df: pd.DataFrame) -> pd.DataFrame:
    return pd.get_dummies(df["subclass"], prefix="subclass")


def one_hot_species(df: pd.DataFrame) -> pd.DataFrame:
    return pd.get_dummies(df["species"], prefix="species")


def build_confounder_matrix(df: pd.DataFrame, include_species: bool = False) -> pd.DataFrame:
    """Covariates to regress out / include as controls: subclass (always),
    species (only for the cross-species ephys arm). donor_id is
    deliberately excluded here -- it is a CV grouping key, not a model
    feature, because a model that is handed donor identity directly could
    memorize per-donor baselines instead of a generalizable relationship.
    """
    parts = [one_hot_subclass(df)]
    if include_species:
        parts.append(one_hot_species(df))
    return pd.concat(parts, axis=1)


def donor_groups(df: pd.DataFrame) -> np.ndarray:
    return df["donor_id"].to_numpy()
