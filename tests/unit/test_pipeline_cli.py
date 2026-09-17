from __future__ import annotations

import csv
import json
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from particle_sim.cli import main
from particle_sim.errors import InputContractError
from particle_sim.io.outputs import CA_COLUMNS, open_hotspot_store
from particle_sim.pipeline import run_nominal
from particle_sim.validation.synthetic_cases import build_case, write_case


@pytest.fixture(scope="module")
def pilot_run(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    base = tmp_path_factory.mktemp("pilot")
    manifest_path = write_case(build_case("pilot"), base / "case")
    run_nominal(manifest_path, base / "run", now=datetime(2026, 9, 18, tzinfo=UTC))
    return manifest_path, base / "run"


def test_run_outputs_shapes(pilot_run: tuple[Path, Path]) -> None:
    _, run_dir = pilot_run
    store = open_hotspot_store(run_dir / "hotspots.zarr")
    assert store["nominal_failure"].shape == (3, 5, 3, 200, 200)
    assert store["legacy_failure"].shape == (3, 5, 3, 200, 200)
    assert store["union_failure"].shape == (3, 5, 3, 200, 200)
    assert list(store.attrs["union_keys"]) == ["lateral_short", "vertical_short", "any"]


def test_ca_table_matches_hotspots(pilot_run: tuple[Path, Path]) -> None:
    _, run_dir = pilot_run
    ca = pd.read_parquet(run_dir / "ca.parquet")
    assert tuple(ca.columns) == CA_COLUMNS
    assert len(ca) == 45
    store = open_hotspot_store(run_dir / "hotspots.zarr")
    steps = list(store.attrs["source_steps"])
    sizes = [float(s) for s in store.attrs["particle_sizes_um"]]
    pairs = list(store.attrs["pair_ids"])
    for row in ca.itertuples(index=False):
        idx = (steps.index(row.source_step), sizes.index(row.particle_size_um),
               pairs.index(row.pair_id))
        assert row.ca_nominal_um2 == pytest.approx(
            float(np.asarray(store["nominal_failure"][idx]).sum()) * 0.25
        )
        assert row.ca_legacy_um2 == pytest.approx(
            float(np.asarray(store["legacy_failure"][idx]).sum()) * 0.25
        )
    assert ca["ca_mean_um2"].isna().all()
    assert (ca["scenario_count"] == 0).all()
    assert ca["ca_nominal_um2"].max() > 0


def test_union_never_exceeds_pair_sum(pilot_run: tuple[Path, Path]) -> None:
    _, run_dir = pilot_run
    union = pd.read_parquet(run_dir / "ca_union.parquet")
    assert len(union) == 45
    assert (union["ca_union_nominal_um2"] <= union["ca_pair_sum_nominal_um2"] + 1e-9).all()


def test_legacy_comparison_fully_explained(pilot_run: tuple[Path, Path]) -> None:
    _, run_dir = pilot_run
    with (run_dir / "legacy_vs_edt.csv").open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 15
    assert all(int(r["unexplained"]) == 0 for r in rows)


def test_run_manifest(pilot_run: tuple[Path, Path]) -> None:
    _, run_dir = pilot_run
    data = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    assert data["run_id"] == "20260918T000000Z_SYN_PILOT"
    assert data["config_hash"]
    assert "manifest" in data["input_hashes"]
    assert any("TODO_PHYSICS" in w for w in data["warnings"])


def test_refuses_non_empty_output(pilot_run: tuple[Path, Path]) -> None:
    manifest_path, run_dir = pilot_run
    with pytest.raises(InputContractError):
        run_nominal(manifest_path, run_dir)


def test_cli_simulate_nominal(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    manifest_path = write_case(build_case("parallel_lines"), tmp_path / "case")
    out = tmp_path / "run"
    assert main(["simulate-nominal", "--product-config", str(manifest_path),
                 "--output", str(out)]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["status"] == "ok"
    assert main(["simulate-nominal", "--product-config", str(manifest_path),
                 "--output", str(out)]) == 2
    assert "입력 오류" in capsys.readouterr().err
