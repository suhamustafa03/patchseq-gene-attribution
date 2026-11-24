import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_generate_demo_fixture_and_run_pipeline_end_to_end(tmp_path):
    fixture_dir = tmp_path / "fixture"
    output_dir = tmp_path / "out"

    gen_script = ROOT / "demo_fixture" / "generate_demo_fixture.py"
    subprocess.run(
        [sys.executable, str(gen_script), "--output", str(fixture_dir), "--seed", "1"],
        check=True,
        cwd=ROOT,
    )
    assert (fixture_dir / "cells.csv").exists()
    assert (fixture_dir / "ephys_features.csv").exists()
    assert (fixture_dir / "gene_expression.csv").exists()

    env_pythonpath = str(ROOT / "src")
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "patchseq_genes.run_pipeline",
            "--input",
            str(fixture_dir),
            "--output",
            str(output_dir),
            "--n-resamples",
            "15",
        ],
        cwd=ROOT,
        env={**__import__("os").environ, "PYTHONPATH": env_pythonpath},
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr

    summary_path = output_dir / "summary.json"
    assert summary_path.exists()
    summary = json.loads(summary_path.read_text())

    for key in (
        "n_cells_total",
        "n_cells_mouse",
        "n_cells_human",
        "n_genes",
        "stable_gene_min_count_threshold",
        "n_stable_genes_total",
        "result_is_null_or_inconclusive_for_any_target",
        "stability_by_target",
    ):
        assert key in summary

    assert summary["n_cells_mouse"] > 0
    assert summary["n_cells_human"] > 0
    assert summary["n_genes"] > 0
