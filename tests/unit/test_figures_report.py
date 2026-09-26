from __future__ import annotations

import json
from pathlib import Path

import pytest
from PIL import Image

from particle_sim.errors import InputContractError
from particle_sim.pipeline import run_nominal
from particle_sim.validation.synthetic_cases import build_case, write_case


@pytest.fixture(scope="module")
def pilot(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    base = tmp_path_factory.mktemp("report")
    manifest_path = write_case(build_case("pilot"), base / "case")
    run_nominal(manifest_path, base / "run")
    return manifest_path, base / "run"


@pytest.fixture(autouse=True)
def _mpl_config(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("MPLCONFIGDIR", str(tmp_path / "mpl"))


def test_make_report(pilot: tuple[Path, Path]) -> None:
    from particle_sim.visualization.reports import make_report

    manifest_path, run_dir = pilot
    report = make_report(run_dir, manifest_path)
    text = report.read_text(encoding="utf-8")
    assert "TODO_PHYSICS" in text
    assert "geometric contact risk" in text
    assert "실제 불량확률이 아님" in text
    assert "failure probability" not in text
    pngs = sorted(report.parent.glob("*.png"))
    assert len(pngs) == 16
    for png in pngs:
        with Image.open(png) as image:
            assert image.size[0] > 0
        assert f"({png.name})" in text
    assert make_report(run_dir, manifest_path) == report


def test_foreign_file_blocks_replacement(pilot: tuple[Path, Path]) -> None:
    from particle_sim.visualization.reports import make_report

    manifest_path, run_dir = pilot
    make_report(run_dir, manifest_path)
    (run_dir / "report" / "notes.txt").write_text("keep", encoding="utf-8")
    try:
        with pytest.raises(InputContractError):
            make_report(run_dir, manifest_path)
    finally:
        (run_dir / "report" / "notes.txt").unlink()


def test_cli_make_report(pilot: tuple[Path, Path], capsys: pytest.CaptureFixture[str]) -> None:
    from particle_sim.cli import main

    manifest_path, run_dir = pilot
    assert main(["make-report", "--run", str(run_dir), "--product-config",
                 str(manifest_path)]) == 0
    assert json.loads(capsys.readouterr().out)["status"] == "ok"


@pytest.mark.parametrize("legacy", [False, True])
def test_probability_report_supports_legacy_tables(tmp_path: Path, legacy: bool) -> None:
    import pandas as pd

    from particle_sim.pipeline import run_probability
    from particle_sim.visualization.reports import make_report

    manifest_path = write_case(build_case("parallel_lines"), tmp_path / "case")
    scenario_path = tmp_path / "scenario.yaml"
    scenario_path.write_text(
        'schema_version: "1.0"\nassumption_set: strict\nevidence: assumed\n'
    )
    run_dir = tmp_path / "run"
    run_probability(manifest_path, scenario_path, run_dir, sobol_exponent_m=1, seed=1)
    if legacy:
        convergence = pd.read_csv(run_dir / "convergence.csv")
        convergence.drop(columns=["scramble_seed"]).to_csv(
            run_dir / "convergence.csv", index=False
        )
        union = pd.read_parquet(run_dir / "ca_union.parquet")
        # Nonzero, distinct legacy values expose swapped labels even on a zero-CA fixture.
        union["ca_union_nominal_um2"] = 12.5
        union["ca_pair_sum_nominal_um2"] = 23.75
        union.drop(columns=["failure_unit_ids", "ca_union_mean_um2", "ca_union_p05_um2",
                            "ca_union_p95_um2", "ca_raw_unit_sum_mean_um2"]).to_parquet(
            run_dir / "ca_union.parquet", index=False
        )
    report = make_report(run_dir, manifest_path)
    text = report.read_text()
    assert "합집합 확률 CA" in text
    assert ("해당 없음" in text) == legacy
    if legacy:
        section = text.split("## 합집합 확률 CA", 1)[1].split("## 수렴", 1)[0]
        rows = [[cell.strip() for cell in line.strip("|").split("|")] for line in section.splitlines()
                if line.startswith("| S_")]
        assert rows
        for row in rows:
            assert row[3:] == ["해당 없음", "12.5", "해당 없음", "해당 없음", "23.75"]
    assert (report.parent / "convergence_S_M1.png").is_file()
