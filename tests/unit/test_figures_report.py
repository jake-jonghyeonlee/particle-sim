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
