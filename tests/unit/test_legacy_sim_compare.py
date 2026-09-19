from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml
import zarr

from particle_sim.comparison.legacy_sim import compare_legacy_sim
from particle_sim.errors import InputContractError
from particle_sim.pipeline import run_nominal
from particle_sim.validation.synthetic_cases import build_case, write_case


def _reference(tmp_path: Path) -> tuple[Path, Path]:
    manifest_path = write_case(build_case("pilot"), tmp_path / "case")
    run = run_nominal(manifest_path, tmp_path / "run")
    store = zarr.open_group(run.hotspots, mode="r")
    records = []
    for step_index, step in enumerate(store.attrs["source_steps"]):
        for size_index, size in enumerate(store.attrs["particle_sizes_um"]):
            name = f"reference-{step_index}-{size_index}.npy"
            array = np.asarray(
                store["nominal_failure"][step_index, size_index, 0], dtype=np.uint8
            )
            np.save(tmp_path / name, array)
            records.append(
                {
                    "source_step": step,
                    "particle_size_um": float(size),
                    "failure_unit_id": store.attrs["pair_ids"][0],
                    "path": name,
                }
            )
    config = {
        "schema_version": "1.0",
        "axes": "yx",
        "origin": "upper_left",
        "pixel_size_x_um": 0.5,
        "pixel_size_y_um": 0.5,
        "particle_size_definition": "diameter",
        "records": records,
    }
    path = tmp_path / "reference.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return path, run.output_dir


def test_roundtrip_reference_has_zero_xor_and_marks_new_pairs(tmp_path: Path) -> None:
    config, run = _reference(tmp_path)
    output = compare_legacy_sim(config, run, tmp_path / "comparison")
    frame = pd.read_parquet(output / "comparison.parquet")
    common = frame[frame["status"] == "common"]
    assert len(common) == 15
    assert not common["xor_pixels"].any()
    assert set(common["difference_classification"]) == {"equal"}
    assert (frame["status"] == "new_pair").any()
    manifest = json.loads((output / "comparison_manifest.json").read_text(encoding="utf-8"))
    assert manifest["unclassified_differences"] == 0


def test_axis_and_pixel_contract_fail_closed(tmp_path: Path) -> None:
    config, run = _reference(tmp_path)
    data = yaml.safe_load(config.read_text(encoding="utf-8"))
    data["axes"] = "xy"
    config.write_text(yaml.safe_dump(data), encoding="utf-8")
    with pytest.raises(InputContractError, match="axes=yx"):
        compare_legacy_sim(config, run, tmp_path / "bad-axis")

    data["axes"] = "yx"
    data["pixel_size_x_um"] = 9.0
    config.write_text(yaml.safe_dump(data), encoding="utf-8")
    with pytest.raises(InputContractError, match="pixel 크기"):
        compare_legacy_sim(config, run, tmp_path / "bad-pixel")


def test_reference_schema_version_fails_closed(tmp_path: Path) -> None:
    config, run = _reference(tmp_path)
    data = yaml.safe_load(config.read_text(encoding="utf-8"))
    data["schema_version"] = "9.0"
    config.write_text(yaml.safe_dump(data), encoding="utf-8")
    with pytest.raises(InputContractError, match="schema_version"):
        compare_legacy_sim(config, run, tmp_path / "bad-schema")
