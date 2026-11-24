import pandas as pd
import pytest

from patchseq_genes import qc


def _valid_cells():
    return pd.DataFrame(
        {
            "cell_id": ["m_1", "m_2", "h_1"],
            "species": ["mouse", "mouse", "human"],
            "subclass": ["Pvalb", "Sst", "Vip"],
            "donor_id": ["d1", "d1", "d2"],
            "has_transcriptome": [True, True, False],
        }
    )


def _valid_ephys():
    return pd.DataFrame(
        {
            "cell_id": ["m_1", "m_2", "h_1"],
            "rheobase_pA": [80.0, 90.0, 100.0],
            "ap_width_ms": [0.5, 0.6, 0.7],
            "input_resistance_MOhm": [150.0, 160.0, 170.0],
            "adaptation_ratio": [1.1, 1.2, 1.3],
        }
    )


def test_validate_cells_valid():
    report = qc.validate_cells(_valid_cells())
    assert report.is_valid
    assert report.n_rows_out == 3


def test_validate_cells_missing_column():
    bad = _valid_cells().drop(columns=["donor_id"])
    report = qc.validate_cells(bad)
    assert not report.is_valid
    assert any("missing required columns" in e for e in report.errors)


def test_validate_cells_bad_species():
    bad = _valid_cells()
    bad.loc[0, "species"] = "elephant"
    report = qc.validate_cells(bad)
    assert not report.is_valid


def test_no_duplicate_cell_ids_enforced():
    dup = _valid_cells()
    dup.loc[2, "cell_id"] = "m_1"
    report = qc.validate_cells(dup)
    assert not report.is_valid
    assert any("duplicate cell_id" in e for e in report.errors)


def test_require_valid_raises_schema_error():
    bad = _valid_cells().drop(columns=["donor_id"])
    report = qc.validate_cells(bad)
    with pytest.raises(qc.SchemaError):
        qc.require_valid(report)


def test_drop_missing_key_features():
    ephys = _valid_ephys()
    ephys.loc[1, "ap_width_ms"] = float("nan")
    cleaned, n_dropped = qc.drop_missing_key_features(ephys)
    assert n_dropped == 1
    assert len(cleaned) == 2


def test_run_qc_end_to_end():
    cells, ephys_clean, reports = qc.run_qc(_valid_cells(), _valid_ephys())
    assert all(r.is_valid for r in reports)
    assert len(ephys_clean) == 3
