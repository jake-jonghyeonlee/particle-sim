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
    run_manifest = json.loads(run.run_manifest.read_text(encoding="utf-8"))
    config = {
        "schema_version": "2.0",
        "product_id": run_manifest["product_id"],
        "design_version": run_manifest["design_version"],
        "input_hashes": {
            name: digest
            for name, digest in run_manifest["input_hashes"].items()
            if name.startswith("masks/")
        },
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


def test_reference_from_another_product_or_mask_is_rejected(tmp_path: Path) -> None:
    config, run = _reference(tmp_path)
    original = yaml.safe_load(config.read_text(encoding="utf-8"))

    data = dict(original, product_id="OTHER")
    config.write_text(yaml.safe_dump(data), encoding="utf-8")
    with pytest.raises(InputContractError, match="product_id"):
        compare_legacy_sim(config, run, tmp_path / "other-product")

    hashes = dict(original["input_hashes"])
    first = sorted(hashes)[0]
    hashes[first] = "0" * 64
    data = dict(original, input_hashes=hashes)
    config.write_text(yaml.safe_dump(data), encoding="utf-8")
    with pytest.raises(InputContractError, match="입력 mask가 다릅니다"):
        compare_legacy_sim(config, run, tmp_path / "other-mask")

    data = dict(original, input_hashes={})
    config.write_text(yaml.safe_dump(data), encoding="utf-8")
    with pytest.raises(InputContractError, match="input_hashes"):
        compare_legacy_sim(config, run, tmp_path / "no-hash")


def test_run_without_pixel_size_is_rejected(tmp_path: Path) -> None:
    config, run = _reference(tmp_path)
    manifest_path = run / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest.pop("pixel_size_x_um")
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(InputContractError, match="pixel 크기가 없어"):
        compare_legacy_sim(config, run, tmp_path / "no-pixel")


def test_shape_mismatch_leaves_no_partial_output(tmp_path: Path) -> None:
    config, run = _reference(tmp_path)
    data = yaml.safe_load(config.read_text(encoding="utf-8"))
    last = data["records"][-1]["path"]
    np.save(config.parent / last, np.zeros((5, 5), dtype=np.uint8))
    output = tmp_path / "bad-shape"
    with pytest.raises(InputContractError, match="shape"):
        compare_legacy_sim(config, run, output)
    assert not output.exists()
