import pandas as pd

from patchseq_genes.features import build_attribution_table, build_cross_species_ephys_table


def _cells():
    return pd.DataFrame(
        {
            "cell_id": ["m_1", "m_2", "h_1"],
            "species": ["mouse", "mouse", "human"],
            "subclass": ["Pvalb", "Sst", "Vip"],
            "donor_id": ["d1", "d1", "d2"],
            "has_transcriptome": [True, True, False],
        }
    )


def _ephys():
    return pd.DataFrame(
        {
            "cell_id": ["m_1", "m_2", "h_1"],
            "rheobase_pA": [80.0, 90.0, 100.0],
            "ap_width_ms": [0.5, 0.6, 0.7],
            "input_resistance_MOhm": [150.0, 160.0, 170.0],
            "adaptation_ratio": [1.1, 1.2, 1.3],
        }
    )


def _expression():
    return pd.DataFrame({"cell_id": ["m_1", "m_2"], "Kcnc1": [5, 7], "Scn1a": [10, 12]})


def test_attribution_table_is_mouse_only():
    table = build_attribution_table(_cells(), _ephys(), _expression())
    assert set(table["species"]) == {"mouse"}
    assert len(table) == 2


def test_cross_species_table_has_both_species_and_no_expression_columns():
    table = build_cross_species_ephys_table(_cells(), _ephys())
    assert set(table["species"]) == {"mouse", "human"}
    assert "Kcnc1" not in table.columns
    assert "Scn1a" not in table.columns
