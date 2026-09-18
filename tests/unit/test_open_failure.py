from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from pydantic import ValidationError

from particle_sim.config import load_manifest
from particle_sim.errors import InputContractError
from particle_sim.failure.open_failure import (
    CONFIRMED,
    PROXY,
    blocked_mask,
    causes_open,
    open_candidates,
    resolve_terminals,
    simulate_open,
)
from particle_sim.io.masks import GridSpec, load_product_masks
from particle_sim.schema import OpenPath
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
