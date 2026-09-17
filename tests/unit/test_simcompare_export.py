from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from particle_sim.cli import main
from particle_sim.errors import InputContractError
from particle_sim.export.simcompare import (
    BUNDLE_COLUMNS,
    CA_RECORD_COLUMNS,
    HOTSPOT_COLUMNS,
    LAYER_STACK_COLUMNS,
    LAYOUT_COLUMNS,
    export_simcompare,
)
from particle_sim.pipeline import run_nominal
from particle_sim.validation.synthetic_cases import build_case, write_case


def _rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


@pytest.fixture(scope="module")
def exported(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path, Path]:
    base = tmp_path_factory.mktemp("export")
    manifest_path = write_case(build_case("pilot"), base / "case")
    run_nominal(manifest_path, base / "run")
    export_simcompare(base / "run", manifest_path, base / "simcompare")
    return manifest_path, base / "run", base / "simcompare"


def test_headers_and_counts(exported: tuple[Path, Path, Path]) -> None:
    _, _, out = exported
    for name, columns, count in (
        ("bundles.csv", BUNDLE_COLUMNS, 9),
        ("layouts.csv", LAYOUT_COLUMNS, 18),
        ("hotspots.csv", HOTSPOT_COLUMNS, 45),
        ("ca_records.csv", CA_RECORD_COLUMNS, 45),
        ("layer_stack.csv", LAYER_STACK_COLUMNS, 5),
    ):
        header, rows = _rows(out / name)
        assert tuple(header) == columns
        assert len(rows) == count


def test_arrays_are_2d_float32(exported: tuple[Path, Path, Path]) -> None:
    _, _, out = exported
    for name in ("layouts.csv", "hotspots.csv"):
        _, rows = _rows(out / name)
        for row in rows[:6]:
            arr = np.load(out / row["path"], allow_pickle=False)
            assert arr.dtype == np.float32
            assert arr.shape == (200, 200)
            assert set(np.unique(arr).tolist()) <= {0.0, 1.0}


def test_ca_value_is_normalized_ca(exported: tuple[Path, Path, Path]) -> None:
    _, run_dir, out = exported
    ca = pd.read_parquet(run_dir / "ca.parquet")
    _, rows = _rows(out / "ca_records.csv")
    row = next(r for r in rows if r["sample_id"].endswith("S_M1-A__B")
               and float(r["particle_size"]) == 10.0)
    expected = ca[(ca.source_step == "S_M1") & (ca.pair_id == "A__B")
                  & (ca.particle_size_um == 10.0)]["ca_nominal_normalized"].iloc[0]
    assert float(row["ca_value"]) == pytest.approx(expected)
    assert row["ca_unit"] == "1"


def test_vertical_bundles_leave_fab_layer_empty(exported: tuple[Path, Path, Path]) -> None:
    _, _, out = exported
    _, bundles = _rows(out / "bundles.csv")
    warnings = json.loads((out / "export_warnings.json").read_text(encoding="utf-8"))
    assert any("fab.json" in w for w in warnings)
    for bundle in bundles:
        if bundle["failure_type"] == "vertical_short":
            assert bundle["assigned_process_layer"] == ""
            assert any(bundle["sample_id"] in w for w in warnings)
        else:
            assert bundle["assigned_process_layer"] != ""
    provenance = json.loads((out / "provenance.json").read_text(encoding="utf-8"))
    assert provenance["is_synthetic"] is True


def test_non_empty_output_rejected(exported: tuple[Path, Path, Path]) -> None:
    manifest_path, run_dir, out = exported
    with pytest.raises(InputContractError):
        export_simcompare(run_dir, manifest_path, out)


def test_cli_export(exported: tuple[Path, Path, Path], tmp_path: Path,
                    capsys: pytest.CaptureFixture[str]) -> None:
    manifest_path, run_dir, _ = exported
    argv = ["export-simcompare", "--run", str(run_dir), "--product-config",
            str(manifest_path), "--output", str(tmp_path / "out")]
    assert main(argv) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["bundles"] == 9
    assert main(argv) == 2
