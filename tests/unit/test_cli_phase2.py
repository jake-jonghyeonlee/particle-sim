from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from particle_sim.cli import main
from particle_sim.validation.synthetic_cases import build_case, write_case, write_fab_table

SCENARIO = """
schema_version: "1.0"
assumption_set: nominal
evidence: assumed
overlay_um:
  M2: {mean: 0.0, sigma: 0.3}
"""


@pytest.fixture
def case(tmp_path: Path) -> tuple[Path, Path]:
    manifest_path = write_case(build_case("parallel_lines"), tmp_path / "case")
    scenario_path = tmp_path / "scenario.yaml"
    scenario_path.write_text(SCENARIO, encoding="utf-8")
    return manifest_path, scenario_path


def test_simulate_probability_cli(case: tuple[Path, Path], tmp_path: Path,
                                  capsys: pytest.CaptureFixture[str]) -> None:
    manifest_path, scenario_path = case
    argv = [
        "simulate-probability",
        "--product-config", str(manifest_path),
        "--scenario-config", str(scenario_path),
        "--sobol-m", "6",
        "--seed", "11",
        "--output", str(tmp_path / "prob"),
    ]
    assert main(argv) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["scenario_count"] == 64
    assert Path(payload["outputs"]["convergence"]).is_file()
    # never overwrite an existing run
    assert main(argv) == 2
    assert "입력 오류" in capsys.readouterr().err


def test_make_synthetic_writes_fab_table(tmp_path: Path,
                                         capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["make-synthetic", "--case", "pilot", "--output", str(tmp_path / "out")]) == 0
    payload = json.loads(capsys.readouterr().out)
    fab = pd.read_csv(payload["fab_table"])
    assert {"FAB_SYN", "FAB_REF"} <= set(fab["fab_id"])
    assert (fab["density_per_mm2"] == fab["particle_count"] / fab["inspection_area_mm2"]).all()
    assert "합성" in payload["note"]


def test_apply_fab_cli(case: tuple[Path, Path], tmp_path: Path,
                       capsys: pytest.CaptureFixture[str]) -> None:
    manifest_path, scenario_path = case
    run_dir = tmp_path / "prob"
    assert main([
        "simulate-probability",
        "--product-config", str(manifest_path),
        "--scenario-config", str(scenario_path),
        "--sobol-m", "6", "--seed", "3", "--output", str(run_dir),
    ]) == 0
    capsys.readouterr()
    fab_path = write_fab_table(tmp_path / "fab")
    output = run_dir / "fab_risk.parquet"
    assert main([
        "apply-fab",
        "--ca", str(run_dir / "ca.parquet"),
        "--fab-data", str(fab_path),
        "--product-config", str(manifest_path),
        "--bin-method", "piecewise_constant",
        "--mode", "both",
        "--reference-fab", "FAB_REF",
        "--output", str(output),
    ]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["modes"] == ["common", "actual"]
    weighted = pd.read_parquet(output)
    assert set(weighted["weighting_mode"]) == {"common", "actual"}
    assert (weighted["expected_defect_count"] >= 0).all()
    # the two fabs have different distributions, so the totals must differ
    totals = weighted.groupby("weighting_mode")["expected_defect_count"].sum()
    assert totals["common"] != totals["actual"]


def test_apply_fab_requires_bin_method(case: tuple[Path, Path], tmp_path: Path) -> None:
    manifest_path, _ = case
    fab_path = write_fab_table(tmp_path / "fab")
    with pytest.raises(SystemExit):
        main([
            "apply-fab",
            "--ca", str(tmp_path / "missing.parquet"),
            "--fab-data", str(fab_path),
            "--product-config", str(manifest_path),
            "--mode", "common",
            "--output", str(tmp_path / "out.parquet"),
        ])
