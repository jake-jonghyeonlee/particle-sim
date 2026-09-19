from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from PIL import Image
from scipy.ndimage import label

from particle_sim.cli import main
from particle_sim.validation.synthetic_cases import CASE_NAMES, build_case, tile_case, write_case


@pytest.mark.parametrize("case", CASE_NAMES)
def test_cases_validate(case: str, tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    manifest_path = write_case(build_case(case), tmp_path / case)
    assert main(["validate-input", "--product-config", str(manifest_path)]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["status"] == "ok"
    assert payload["shape_yx"] == [200, 200]
    assert payload["source_steps"] == ["S_M1", "S_M2", "S_PLN"]
    assert all(
        row["min_xy_distance_um"] is not None
        and row["min_z_gap_um"] is not None
        and row["min_3d_distance_um"] is not None
        for row in payload["electrode_pair_audit"]
    )
    assert any("TODO_PHYSICS" in w for w in payload["warnings"])


def test_case_geometry_matches_contract() -> None:
    case = build_case("parallel_lines")
    a = case.electrode_masks["A"]
    b = case.electrode_masks["B"]
    assert a[:, 80:90].all() and not a[:, 90:96].any()
    assert b[:, 96:106].all()
    assert np.array_equal(case.layer_masks["M2"], a | b)
    assert not case.layer_masks["M1"].any()
    assert case.manifest["product_id"] == "SYN_PARALLEL_LINES"


def test_complex_multilayer_case_separates_electrodes_and_covers_layers() -> None:
    case = build_case("complex_multilayer")
    assert len(case.electrode_masks) == 10
    assert len(case.manifest["electrode_pairs"]) == 16
    electrodes = {row["electrode_id"]: row for row in case.manifest["electrodes"]}
    for layer_name, layer_mask in case.layer_masks.items():
        members = [
            case.electrode_masks[electrode_id]
            for electrode_id, row in electrodes.items()
            if row["layer_name"] == layer_name
        ]
        accumulated = np.zeros(layer_mask.shape, dtype=bool)
        for member in members:
            assert not np.any(accumulated & member)
            accumulated |= member
        assert np.array_equal(accumulated, layer_mask)
    assert label(case.electrode_masks["G3"])[1] == 2


def test_complex_case_tiles_to_1000_without_changing_pixel_size() -> None:
    large = tile_case(build_case("complex_multilayer"), (5, 5))
    assert {mask.shape for mask in large.layer_masks.values()} == {(1000, 1000)}
    assert {mask.shape for mask in large.electrode_masks.values()} == {(1000, 1000)}
    coordinate = large.manifest["coordinate_system"]
    assert coordinate["pixel_size_x_um"] == 0.5
    assert coordinate["pixel_size_y_um"] == 0.5
    assert coordinate["roi_width_um"] == 500.0
    assert coordinate["roi_height_um"] == 500.0


def test_unknown_case() -> None:
    with pytest.raises(ValueError):
        build_case("nope")


def test_make_synthetic_cli(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["make-synthetic", "--case", "pilot", "--output", str(tmp_path / "out")]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert Path(payload["manifest"]).is_file()


def test_corrupt_mask_is_reported(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    manifest_path = write_case(build_case("parallel_lines"), tmp_path / "case")
    bad = np.full((200, 200), 128, dtype=np.uint8)
    Image.fromarray(bad, mode="L").save(tmp_path / "case" / "masks" / "electrode_A.png")
    assert main(["validate-input", "--product-config", str(manifest_path)]) == 2
    err = capsys.readouterr().err
    assert "입력 오류" in err
    assert "field=" in err
