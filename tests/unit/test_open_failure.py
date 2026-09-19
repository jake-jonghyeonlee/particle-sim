from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import yaml
from PIL import Image
from pydantic import ValidationError

from particle_sim.config import load_manifest
from particle_sim.errors import InputContractError
from particle_sim.failure.open_failure import (
    CONFIRMED,
    PROXY,
    blocked_mask,
    causes_open,
    causes_open_regions,
    open_candidates,
    particle_footprint,
    resolve_terminal_regions,
    resolve_terminals,
    simulate_open,
)
from particle_sim.geometry.contact import DistanceCache
from particle_sim.io.masks import GridSpec, load_product_masks
from particle_sim.io.stack import build_stack_z
from particle_sim.schema import OpenPath
from particle_sim.units import evaluate_units
from particle_sim.validation.synthetic_cases import build_case, write_case


def _grid(shape=(40, 60), py=0.5, px=0.5) -> GridSpec:
    return GridSpec(shape_yx=shape, pixel_size_y_um=py, pixel_size_x_um=px)


def _wire(shape=(40, 60)) -> np.ndarray:
    wire = np.zeros(shape, dtype=bool)
    wire[18:23, 5:55] = True
    return wire


def _brute_force(wire: np.ndarray, grid: GridSpec, radius_um: float, terminals) -> np.ndarray:
    yy, xx = np.mgrid[0 : wire.shape[0], 0 : wire.shape[1]]
    result = np.zeros(wire.shape, dtype=bool)
    for y in range(wire.shape[0]):
        for x in range(wire.shape[1]):
            distance = np.hypot(
                (yy - y) * grid.pixel_size_y_um, (xx - x) * grid.pixel_size_x_um
            )
            result[y, x] = causes_open(
                wire, blocked_mask(distance, radius_um, 0.0), *terminals
            )
    return result


def _case(tmp_path: Path, name: str):
    manifest_path = write_case(build_case(name), tmp_path / name)
    manifest = load_manifest(manifest_path)
    return manifest, load_product_masks(manifest, manifest_path)


def test_causes_open_needs_terminal_separation() -> None:
    wire = _wire()
    blocked = np.zeros(wire.shape, dtype=bool)
    assert not causes_open(wire, blocked, (20, 6), (20, 53))
    blocked[:, 30:33] = True
    assert causes_open(wire, blocked, (20, 6), (20, 53))
    # a blocked terminal is an open path as well
    only_terminal = np.zeros(wire.shape, dtype=bool)
    only_terminal[20, 6] = True
    assert causes_open(wire, only_terminal, (20, 6), (20, 53))


def test_terminal_regions_allow_partial_survival_but_detect_cut() -> None:
    wire = _wire()
    terminal_a = np.zeros(wire.shape, dtype=bool)
    terminal_b = np.zeros(wire.shape, dtype=bool)
    terminal_a[18:23, 5:8] = True
    terminal_b[18:23, 52:55] = True
    partial = np.zeros(wire.shape, dtype=bool)
    partial[18:20, 5:8] = True
    assert not causes_open_regions(wire, partial, terminal_a, terminal_b)
    covered = np.zeros(wire.shape, dtype=bool)
    covered[terminal_a] = True
    assert causes_open_regions(wire, covered, terminal_a, terminal_b)
    cut = np.zeros(wire.shape, dtype=bool)
    cut[:, 29:32] = True
    assert causes_open_regions(wire, cut, terminal_a, terminal_b)


@pytest.mark.parametrize("diameter", [1.0, 2.0, 3.0, 4.0, 6.0])
def test_candidate_reduction_never_misses_a_brute_force_hit(diameter: float) -> None:
    wire = _wire()
    grid = _grid()
    terminals = ((20, 6), (20, 53))
    brute = _brute_force(wire, grid, diameter / 2.0, terminals)
    candidates = open_candidates(wire, grid, diameter / 2.0, 0.0, terminals)
    assert not (brute & ~candidates).any()


def test_candidate_reduction_on_anisotropic_grid() -> None:
    wire = _wire()
    grid = _grid(py=1.0, px=0.5)
    terminals = ((20, 6), (20, 53))
    brute = _brute_force(wire, grid, 2.0, terminals)
    candidates = open_candidates(wire, grid, 2.0, 0.0, terminals)
    assert not (brute & ~candidates).any()


def test_straight_wire_is_cut_by_a_large_particle(tmp_path: Path) -> None:
    manifest, masks = _case(tmp_path, "open_straight")
    small = simulate_open(manifest, masks, manifest.open_paths[0], 2.0)
    large = simulate_open(manifest, masks, manifest.open_paths[0], 10.0)
    assert large.failure.sum() > small.failure.sum()
    assert large.classification == CONFIRMED
    # the cut positions sit on the wire, not somewhere in the empty field
    ys = np.nonzero(large.failure)[0]
    assert ys.min() >= 90 and ys.max() <= 110


def test_cropped_connectivity_matches_full_grid_reference(tmp_path: Path) -> None:
    manifest, masks = _case(tmp_path, "open_straight")
    open_path = manifest.open_paths[0]
    electrode = masks.electrodes[open_path.electrode]
    terminal_a, terminal_b, _classification, _kind = resolve_terminal_regions(
        open_path, masks, electrode, manifest.product_id
    )
    diameter = 6.0
    radius = diameter / 2.0
    candidates = open_candidates(electrode, masks.grid, radius, manifest.contact_margin_um)
    stamp = particle_footprint(masks.grid, radius, manifest.contact_margin_um)
    half_y, half_x = stamp.shape[0] // 2, stamp.shape[1] // 2
    height, width = electrode.shape
    blocked = np.zeros(electrode.shape, dtype=bool)
    expected = np.zeros(electrode.shape, dtype=bool)
    ys, xs = np.nonzero(candidates)
    for y, x in zip(ys.tolist(), xs.tolist(), strict=True):
        top, bottom = max(0, y - half_y), min(height, y + half_y + 1)
        left, right = max(0, x - half_x), min(width, x + half_x + 1)
        window = stamp[
            top - (y - half_y) : stamp.shape[0] - ((y + half_y + 1) - bottom),
            left - (x - half_x) : stamp.shape[1] - ((x + half_x + 1) - right),
        ]
        blocked[top:bottom, left:right] = window
        expected[y, x] = causes_open_regions(electrode, blocked, terminal_a, terminal_b)
        blocked[top:bottom, left:right] = False

    actual = simulate_open(manifest, masks, open_path, diameter)
    assert np.array_equal(actual.failure, expected)


def test_detour_keeps_terminals_connected(tmp_path: Path) -> None:
    straight_manifest, straight_masks = _case(tmp_path, "open_straight")
    detour_manifest, detour_masks = _case(tmp_path, "open_detour")
    straight = simulate_open(straight_manifest, straight_masks,
                             straight_manifest.open_paths[0], 6.0)
    detour = simulate_open(detour_manifest, detour_masks, detour_manifest.open_paths[0], 6.0)
    assert detour.failure.sum() < straight.failure.sum()
    # inside the bypassed span nothing opens any more
    assert not detour.failure[:, 80:120].any()
    assert straight.failure[:, 80:120].any()


def test_branch_cut_is_not_an_open(tmp_path: Path) -> None:
    straight_manifest, straight_masks = _case(tmp_path, "open_straight")
    branch_manifest, branch_masks = _case(tmp_path, "open_branch")
    straight = simulate_open(straight_manifest, straight_masks,
                             straight_manifest.open_paths[0], 6.0)
    branch = simulate_open(branch_manifest, branch_masks, branch_manifest.open_paths[0], 6.0)
    # cutting the dead-end branch (rows above the main line) never opens the path
    assert not branch.failure[:90, :].any()
    assert abs(int(branch.failure.sum()) - int(straight.failure.sum())) <= 10


def test_missing_terminals_use_proxy_or_fail(tmp_path: Path) -> None:
    manifest, masks = _case(tmp_path, "open_straight")
    wire = masks.electrodes["W"]
    proxy = OpenPath(path_id="P", electrode="W", fallback="roi_boundary")
    terminal_a, terminal_b, classification = resolve_terminals(proxy, wire, "SYN")
    assert classification == PROXY
    assert wire[terminal_a] and wire[terminal_b]
    assert terminal_a != terminal_b
    # a path without terminals must name a fallback; "error" keeps it fail-closed
    with pytest.raises(ValidationError):
        OpenPath(path_id="P", electrode="W")
    with pytest.raises(ValidationError):
        OpenPath(path_id="P", electrode="W", fallback="skeleton_cut")


def test_terminal_outside_electrode_is_rejected(tmp_path: Path) -> None:
    manifest, masks = _case(tmp_path, "open_straight")
    off_wire = OpenPath(
        path_id="P", electrode="W", terminal_a_yx_px=(0, 0), terminal_b_yx_px=(100, 178)
    )
    with pytest.raises(InputContractError, match="전극 위에 없습니다"):
        resolve_terminals(off_wire, masks.electrodes["W"], "SYN")
    outside = OpenPath(
        path_id="P", electrode="W", terminal_a_yx_px=(999, 0), terminal_b_yx_px=(100, 178)
    )
    with pytest.raises(InputContractError, match="이미지 밖"):
        resolve_terminals(outside, masks.electrodes["W"], "SYN")


def test_unknown_electrode(tmp_path: Path) -> None:
    manifest, masks = _case(tmp_path, "open_straight")
    unknown = OpenPath(
        path_id="P", electrode="NOPE", terminal_a_yx_px=(100, 22), terminal_b_yx_px=(100, 178)
    )
    with pytest.raises(InputContractError):
        simulate_open(manifest, masks, unknown, 4.0)


def test_candidate_reduction_covers_a_corner_wire() -> None:
    """The Codex counter-example: skeletonize used to drop corner pixels."""
    grid = _grid(shape=(12, 12), py=1.0, px=1.0)
    wire = np.zeros((12, 12), dtype=bool)
    wire[2:8, 3] = True
    wire[7, 3:9] = True
    terminals = ((2, 3), (7, 8))
    for diameter in (1.0, 2.0, 4.0, 8.0):
        brute = _brute_force(wire, grid, diameter / 2.0, terminals)
        candidates = open_candidates(wire, grid, diameter / 2.0, 0.0, terminals)
        assert not (brute & ~candidates).any()


@pytest.mark.parametrize("margin", [0.25, 1.0])
def test_candidate_reduction_with_contact_margin(margin: float) -> None:
    wire = _wire()
    grid = _grid()
    terminals = ((20, 6), (20, 53))
    yy, xx = np.mgrid[0 : wire.shape[0], 0 : wire.shape[1]]
    brute = np.zeros(wire.shape, dtype=bool)
    for y in range(wire.shape[0]):
        for x in range(wire.shape[1]):
            distance = np.hypot((yy - y) * grid.pixel_size_y_um,
                                (xx - x) * grid.pixel_size_x_um)
            brute[y, x] = causes_open(wire, blocked_mask(distance, 1.0, margin), *terminals)
    candidates = open_candidates(wire, grid, 1.0, margin, terminals)
    assert not (brute & ~candidates).any()


def test_terminal_region_masks_load_and_drive_open(tmp_path: Path) -> None:
    manifest_path = write_case(build_case("open_straight"), tmp_path / "regions")
    data = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    shape = (200, 200)
    terminal_a = np.zeros(shape, dtype=np.uint8)
    terminal_b = np.zeros(shape, dtype=np.uint8)
    terminal_a[98:103, 20:25] = 255
    terminal_b[98:103, 175:180] = 255
    Image.fromarray(terminal_a, mode="L").save(manifest_path.parent / "masks/terminal_a.png")
    Image.fromarray(terminal_b, mode="L").save(manifest_path.parent / "masks/terminal_b.png")
    path = data["open_paths"][0]
    path.pop("terminal_a_yx_px")
    path.pop("terminal_b_yx_px")
    path["terminal_a_mask_path"] = "masks/terminal_a.png"
    path["terminal_b_mask_path"] = "masks/terminal_b.png"
    path["applicable_particle_source_steps"] = [data["particle_source_steps"][0]["step_id"]]
    manifest_path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    manifest = load_manifest(manifest_path)
    masks = load_product_masks(manifest, manifest_path)
    result = simulate_open(manifest, masks, manifest.open_paths[0], 6.0)
    assert result.terminal_kind == "region_mask"
    assert result.classification == CONFIRMED
    assert result.failure.any()


def test_open_result_cache_reuses_identical_geometry_across_source_steps(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest, masks = _case(tmp_path, "open_straight")
    stack_z = build_stack_z(manifest)
    distance_cache = DistanceCache(masks)
    open_cache = {}
    calls = 0

    original = simulate_open

    def counted_simulate_open(*args, **kwargs):
        nonlocal calls
        calls += 1
        return original(*args, **kwargs)

    monkeypatch.setattr("particle_sim.units.simulate_open", counted_simulate_open)
    results = [
        evaluate_units(
            manifest,
            masks,
            stack_z,
            step.step_id,
            6.0,
            distance_cache,
            open_cache,
        )[0]
        for step in manifest.particle_source_steps
    ]

    assert calls == 1
    assert len(open_cache) == 1
    assert all(np.array_equal(result.failure, results[0].failure) for result in results[1:])

    evaluate_units(
        manifest,
        masks,
        stack_z,
        manifest.particle_source_steps[0].step_id,
        8.0,
        distance_cache,
        open_cache,
    )
    assert calls == 2
