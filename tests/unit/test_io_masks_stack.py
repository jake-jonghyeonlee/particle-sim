"""io/masks.py 와 io/stack.py 단위 테스트."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import tifffile
import yaml
from PIL import Image

from particle_sim.config import load_manifest
from particle_sim.errors import InputContractError
from particle_sim.io.masks import GridSpec, MaskSet, load_binary_mask, load_product_masks
from particle_sim.io.stack import (
    build_stack_z,
    landing_surface_z_um,
    layers_after,
    layers_up_to,
)
from particle_sim.schema import ProductManifest

# ── helpers ──────────────────────────────────────────────────────────────────


def _write_png(path: Path, arr: np.ndarray) -> None:
    """bool/0-1/0-255 배열을 PNG로 저장."""
    if arr.dtype == bool:
        arr = arr.astype(np.uint8) * 255
    elif arr.max() <= 1:
        arr = arr.astype(np.uint8) * 255
    Image.fromarray(arr, mode="L").save(path)


def _write_tif(path: Path, arr: np.ndarray) -> None:
    """0/255 uint8 배열을 TIFF로 저장."""
    if arr.dtype == bool:
        arr = arr.astype(np.uint8) * 255
    elif arr.max() <= 1:
        arr = arr.astype(np.uint8) * 255
    tifffile.imwrite(str(path), arr)


def _base_manifest_dict(tmp_path: Path, *, n_electrodes: int = 1) -> dict:
    """최소 유효 manifest dict를 생성한다."""
    layers = [
        {
            "layer_name": "L0",
            "process_order": 0,
            "material_type": "conductor",
            "thickness_um": 1.0,
            "patterned": True,
        }
    ]
    layer_masks = [{"layer_name": "L0", "path": "masks/L0.png"}]
    electrodes = []
    pairs = []
    if n_electrodes >= 1:
        electrodes.append(
            {
                "electrode_id": "E0",
                "layer_name": "L0",
                "net_id": "N0",
                "electrode_role": "anode",
                "path": "masks/E0.png",
            }
        )
    if n_electrodes >= 2:
        electrodes.append(
            {
                "electrode_id": "E1",
                "layer_name": "L0",
                "net_id": "N1",
                "electrode_role": "cathode",
                "path": "masks/E1.png",
            }
        )
        pairs.append(
            {
                "pair_id": "P0",
                "electrode_a": "E0",
                "electrode_b": "E1",
                "failure_mode": "lateral_short",
            }
        )
    return {
        "schema_version": "1.0",
        "product_id": "P1",
        "product_family": "fam",
        "design_version": "v1",
        "design_lineage_id": "dl1",
        "fab_id": "fab1",
        "is_synthetic": True,
        "coordinate_system": {
            "origin": "upper_left",
            "x_direction": "right",
            "y_direction": "down",
            "pixel_size_x_um": 1.0,
            "pixel_size_y_um": 1.0,
            "roi_width_um": 10.0,
            "roi_height_um": 10.0,
        },
        "stack": layers,
        "layer_masks": layer_masks,
        "electrodes": electrodes,
        "particle_source_steps": [
            {
                "step_id": "S0",
                "lands_after_layer": "L0",
                "surface_model": "flat",
                "fab_distribution_layer": "L0",
            }
        ],
        "electrode_pairs": pairs,
        "embedding_depth_um": 0.0,
        "contact_margin_um": 0.0,
        "post_step_policy": {"mode": "exclude"},
    }


def _write_manifest(tmp_path: Path, mdict: dict) -> Path:
    mpath = tmp_path / "manifest.yaml"
    with mpath.open("w") as f:
        yaml.dump(mdict, f)
    return mpath


def _setup_valid_masks(tmp_path: Path) -> Path:
    """유효한 마스크 파일들을 tmp_path에 쓰고 manifest yaml 경로를 반환한다.

    layer mask == electrode mask (union 일치).
    """
    d = tmp_path / "masks"
    d.mkdir(exist_ok=True)
    mask = np.zeros((10, 10), dtype=bool)
    mask[2:8, 2:8] = True
    _write_png(d / "L0.png", mask)
    _write_png(d / "E0.png", mask)
    mdict = _base_manifest_dict(tmp_path)
    return _write_manifest(tmp_path, mdict)


# ── GridSpec ─────────────────────────────────────────────────────────────────


class TestGridSpec:
    def test_pixel_area(self) -> None:
        g = GridSpec(shape_yx=(10, 20), pixel_size_y_um=0.5, pixel_size_x_um=1.0)
        assert g.pixel_area_um2 == pytest.approx(0.5)

    def test_roi_area(self) -> None:
        g = GridSpec(shape_yx=(10, 20), pixel_size_y_um=0.5, pixel_size_x_um=1.0)
        assert g.roi_area_um2 == pytest.approx(10 * 0.5 * 20 * 1.0)

    @pytest.mark.parametrize(
        "py,px,expected",
        [(1.0, 1.0, True), (1.0, 1.0000000001, True), (1.0, 2.0, False)],
    )
    def test_is_isotropic(self, py: float, px: float, expected: bool) -> None:
        g = GridSpec(shape_yx=(10, 10), pixel_size_y_um=py, pixel_size_x_um=px)
        assert g.is_isotropic is expected


# ── load_binary_mask ─────────────────────────────────────────────────────────


class TestLoadBinaryMask:
    @pytest.mark.parametrize(
        "arr,fmt",
        [
            (np.zeros((5, 5), dtype=np.uint8), "png"),
            (np.ones((5, 5), dtype=np.uint8) * 255, "png"),
            (np.array([[0, 1], [1, 0]], dtype=np.uint8), "png"),
            (np.zeros((5, 5), dtype=np.uint8), "tif"),
            (np.ones((5, 5), dtype=np.uint8) * 255, "tif"),
        ],
    )
    def test_valid_binary(self, tmp_path: Path, arr: np.ndarray, fmt: str) -> None:
        path = tmp_path / f"m.{fmt}"
        if fmt == "png":
            _write_png(path, arr)
        else:
            _write_tif(path, arr)
        result = load_binary_mask(path, product_id="P1", field="test")
        assert result.dtype == bool
        assert result.shape == arr.shape

    def test_gray_rejected(self, tmp_path: Path) -> None:
        arr = np.full((5, 5), 128, dtype=np.uint8)
        path = tmp_path / "gray.png"
        _write_png(path, arr)
        with pytest.raises(InputContractError, match="이진 마스크가 아님"):
            load_binary_mask(path, product_id="P1", field="test")

    def test_rgb_rejected(self, tmp_path: Path) -> None:
        arr = np.zeros((5, 5, 3), dtype=np.uint8)
        path = tmp_path / "rgb.png"
        Image.fromarray(arr, mode="RGB").save(path)
        with pytest.raises(InputContractError, match="2-D 단일 채널"):
            load_binary_mask(path, product_id="P1", field="test")

    def test_bad_extension(self, tmp_path: Path) -> None:
        path = tmp_path / "m.bmp"
        path.write_bytes(b"fake")
        with pytest.raises(InputContractError, match="지원하지 않는"):
            load_binary_mask(path, product_id="P1", field="test")

    def test_3d_tif_rejected(self, tmp_path: Path) -> None:
        arr = np.zeros((5, 5, 2), dtype=np.uint8)
        path = tmp_path / "3d.tif"
        tifffile.imwrite(str(path), arr)
        with pytest.raises(InputContractError, match="2-D"):
            load_binary_mask(path, product_id="P1", field="test")


# ── load_product_masks ───────────────────────────────────────────────────────


class TestLoadProductMasks:
    def test_valid(self, tmp_path: Path) -> None:
        mpath = _setup_valid_masks(tmp_path)
        manifest = load_manifest(mpath)
        ms = load_product_masks(manifest, mpath)
        assert isinstance(ms, MaskSet)
        assert ms.grid.shape_yx == (10, 10)
        assert "L0" in ms.layers
        assert "E0" in ms.electrodes
        assert ms.layers["L0"].dtype == bool
        assert ms.electrodes["E0"].dtype == bool
        assert len(ms.input_hashes) == 2

    def test_shape_mismatch(self, tmp_path: Path) -> None:
        d = tmp_path / "masks"
        d.mkdir(exist_ok=True)
        mask = np.zeros((10, 10), dtype=bool)
        mask[2:8, 2:8] = True
        _write_png(d / "L0.png", mask)
        _write_png(d / "E0.png", np.ones((8, 8), dtype=bool))
        mpath = _write_manifest(tmp_path, _base_manifest_dict(tmp_path))
        manifest = load_manifest(mpath)
        with pytest.raises(InputContractError, match="shape 불일치"):
            load_product_masks(manifest, mpath)

    def test_roi_mismatch(self, tmp_path: Path) -> None:
        d = tmp_path / "masks"
        d.mkdir(exist_ok=True)
        mask = np.zeros((10, 10), dtype=bool)
        mask[2:8, 2:8] = True
        _write_png(d / "L0.png", mask)
        _write_png(d / "E0.png", mask)
        mdict = _base_manifest_dict(tmp_path)
        mdict["coordinate_system"]["roi_width_um"] = 99.0
        mpath = _write_manifest(tmp_path, mdict)
        manifest = load_manifest(mpath)
        with pytest.raises(InputContractError, match="ROI 폭 불일치"):
            load_product_masks(manifest, mpath)

    def test_empty_electrode(self, tmp_path: Path) -> None:
        d = tmp_path / "masks"
        d.mkdir(exist_ok=True)
        mask = np.zeros((10, 10), dtype=bool)
        mask[2:8, 2:8] = True
        _write_png(d / "L0.png", mask)
        _write_png(d / "E0.png", np.zeros((10, 10), dtype=bool))
        mpath = _write_manifest(tmp_path, _base_manifest_dict(tmp_path))
        manifest = load_manifest(mpath)
        with pytest.raises(InputContractError, match="비어있습니다"):
            load_product_masks(manifest, mpath)

    def test_electrode_overlap(self, tmp_path: Path) -> None:
        d = tmp_path / "masks"
        d.mkdir(exist_ok=True)
        # layer mask = union of both electrodes (no overlap between them)
        e0 = np.zeros((10, 10), dtype=bool)
        e0[2:5, 2:5] = True
        e1 = np.zeros((10, 10), dtype=bool)
        e1[5:8, 5:8] = True
        layer = e0 | e1
        _write_png(d / "L0.png", layer)
        _write_png(d / "E0.png", e0)
        _write_png(d / "E1.png", e1)
        mdict = _base_manifest_dict(tmp_path, n_electrodes=2)
        mpath = _write_manifest(tmp_path, mdict)
        manifest = load_manifest(mpath)
        # no overlap -> should pass
        ms = load_product_masks(manifest, mpath)
        assert "E0" in ms.electrodes
        assert "E1" in ms.electrodes

    def test_electrode_overlap_fail(self, tmp_path: Path) -> None:
        d = tmp_path / "masks"
        d.mkdir(exist_ok=True)
        e0 = np.zeros((10, 10), dtype=bool)
        e0[2:6, 2:6] = True
        e1 = np.zeros((10, 10), dtype=bool)
        e1[4:8, 4:8] = True  # overlaps with e0 at [4:6, 4:6]
        layer = e0 | e1
        _write_png(d / "L0.png", layer)
        _write_png(d / "E0.png", e0)
        _write_png(d / "E1.png", e1)
        mdict = _base_manifest_dict(tmp_path, n_electrodes=2)
        mpath = _write_manifest(tmp_path, mdict)
        manifest = load_manifest(mpath)
        with pytest.raises(InputContractError, match="겹칩니다"):
            load_product_masks(manifest, mpath)

    def test_union_mismatch(self, tmp_path: Path) -> None:
        d = tmp_path / "masks"
        d.mkdir(exist_ok=True)
        # layer mask is larger than electrode -> union != layer
        layer = np.ones((10, 10), dtype=bool)
        em = np.zeros((10, 10), dtype=bool)
        em[2:6, 2:6] = True
        _write_png(d / "L0.png", layer)
        _write_png(d / "E0.png", em)
        mpath = _write_manifest(tmp_path, _base_manifest_dict(tmp_path))
        manifest = load_manifest(mpath)
        with pytest.raises(InputContractError, match="합집합 불일치"):
            load_product_masks(manifest, mpath)


# ── build_stack_z ────────────────────────────────────────────────────────────


def _stack_manifest(tmp_path: Path, *, z_bottom: list[float | None] | None = None) -> ProductManifest:
    """3-layer stack manifest를 생성한다."""
    layers = [
        {
            "layer_name": "L0",
            "process_order": 0,
            "material_type": "conductor",
            "thickness_um": 1.0,
            "patterned": True,
        },
        {
            "layer_name": "L1",
            "process_order": 1,
            "material_type": "inorganic",
            "thickness_um": 2.0,
            "patterned": False,
        },
        {
            "layer_name": "L2",
            "process_order": 2,
            "material_type": "conductor",
            "thickness_um": 3.0,
            "patterned": True,
        },
    ]
    if z_bottom is not None:
        for i, zb in enumerate(z_bottom):
            if zb is not None:
                layers[i]["z_bottom_um"] = zb

    mdict = {
        "schema_version": "1.0",
        "product_id": "P1",
        "product_family": "fam",
        "design_version": "v1",
        "design_lineage_id": "dl1",
        "fab_id": "fab1",
        "is_synthetic": True,
        "coordinate_system": {
            "origin": "upper_left",
            "x_direction": "right",
            "y_direction": "down",
            "pixel_size_x_um": 1.0,
            "pixel_size_y_um": 1.0,
            "roi_width_um": 10.0,
            "roi_height_um": 10.0,
        },
        "stack": layers,
        "layer_masks": [
            {"layer_name": "L0", "path": "masks/L0.png"},
            {"layer_name": "L2", "path": "L2.png"},
        ],
        "electrodes": [],
        "particle_source_steps": [
            {
                "step_id": "S0",
                "lands_after_layer": "L1",
                "surface_model": "flat",
                "fab_distribution_layer": "L1",
            }
        ],
        "electrode_pairs": [],
        "embedding_depth_um": 0.0,
        "contact_margin_um": 0.0,
        "post_step_policy": {"mode": "exclude"},
    }
    mpath = tmp_path / "manifest.yaml"
    with mpath.open("w") as f:
        yaml.dump(mdict, f)
    return load_manifest(mpath)


class TestBuildStackZ:
    def test_cumulative_z(self, tmp_path: Path) -> None:
        manifest = _stack_manifest(tmp_path)
        sz = build_stack_z(manifest)
        assert list(sz.keys()) == ["L0", "L1", "L2"]
        assert sz["L0"].z_bottom_um == pytest.approx(0.0)
        assert sz["L0"].z_top_um == pytest.approx(1.0)
        assert sz["L1"].z_bottom_um == pytest.approx(1.0)
        assert sz["L1"].z_top_um == pytest.approx(3.0)
        assert sz["L2"].z_bottom_um == pytest.approx(3.0)
        assert sz["L2"].z_top_um == pytest.approx(6.0)

    def test_z_bottom_crosscheck_ok(self, tmp_path: Path) -> None:
        manifest = _stack_manifest(tmp_path, z_bottom=[0.0, 1.0, 3.0])
        sz = build_stack_z(manifest)
        assert sz["L1"].z_bottom_um == pytest.approx(1.0)

    def test_z_bottom_crosscheck_fail(self, tmp_path: Path) -> None:
        manifest = _stack_manifest(tmp_path, z_bottom=[0.0, 99.0, None])
        with pytest.raises(InputContractError, match="z_bottom_um 불일치"):
            build_stack_z(manifest)


# ── landing_surface_z_um / layers_up_to / layers_after ──────────────────────


class TestStepLayers:
    @pytest.fixture()
    def manifest_and_stack(self, tmp_path: Path):
        manifest = _stack_manifest(tmp_path)
        sz = build_stack_z(manifest)
        return manifest, sz

    def test_landing_z(self, manifest_and_stack) -> None:
        manifest, sz = manifest_and_stack
        z = landing_surface_z_um(manifest, "S0", sz)
        # S0 lands after L1, z_top of L1 = 3.0
        assert z == pytest.approx(3.0)

    def test_landing_unknown_step(self, manifest_and_stack) -> None:
        manifest, sz = manifest_and_stack
        with pytest.raises(InputContractError, match="알 수 없는 step_id"):
            landing_surface_z_um(manifest, "NOPE", sz)

    def test_layers_up_to(self, manifest_and_stack) -> None:
        manifest, sz = manifest_and_stack
        step = manifest.particle_source_steps[0]  # lands after L1 (order=1)
        up = layers_up_to(step, sz)
        assert [lz.layer_name for lz in up] == ["L0", "L1"]

    def test_layers_after(self, manifest_and_stack) -> None:
        manifest, sz = manifest_and_stack
        step = manifest.particle_source_steps[0]  # lands after L1 (order=1)
        after = layers_after(step, sz)
        assert [lz.layer_name for lz in after] == ["L2"]


def test_missing_mask_file_is_contract_error(tmp_path: Path) -> None:
    from particle_sim.validation.synthetic_cases import build_case, write_case

    manifest_path = write_case(build_case("parallel_lines"), tmp_path)
    (tmp_path / "masks" / "electrode_A.png").unlink()
    manifest = load_manifest(manifest_path)
    with pytest.raises(InputContractError) as exc:
        load_product_masks(manifest, manifest_path)
    message = str(exc.value)
    assert "file=" in message
    assert "electrodes[0].path" in message


def test_corrupt_mask_file_is_contract_error(tmp_path: Path) -> None:
    from particle_sim.validation.synthetic_cases import build_case, write_case

    manifest_path = write_case(build_case("parallel_lines"), tmp_path)
    (tmp_path / "masks" / "layer_M2.png").write_bytes(b"not a png")
    manifest = load_manifest(manifest_path)
    with pytest.raises(InputContractError) as exc:
        load_product_masks(manifest, manifest_path)
    assert "layer_masks[" in str(exc.value)
