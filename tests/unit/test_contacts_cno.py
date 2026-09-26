"""Contact-derived open terminals, CNO and the open baseline checks (2026-09-19 review)."""

from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml
from PIL import Image
from pydantic import ValidationError
from scipy.ndimage import correlate

from particle_sim.aggregation.risk_allocation import allocate_risk
from particle_sim.aggregation.risk_tables import build_risk_tables
from particle_sim.comparison.legacy_sim import compare_legacy_sim
from particle_sim.config import load_manifest
from particle_sim.errors import InputContractError
from particle_sim.failure.auto_open import (
    AUTO_CONFIRMED,
    INSUFFICIENT_TERMINALS,
    simulate_auto_open,
)
from particle_sim.failure.contact_not_open import (
    CNO_MODE,
    CONTACT_LOST,
    hole_blocked_centres,
    simulate_cno,
)
from particle_sim.failure.open_failure import (
    BASELINE_OPEN,
    TERMINAL_LOST,
    check_open_baselines,
    particle_footprint,
    simulate_open,
)
from particle_sim.failure.union import ANY_MODE, union_by_mode
from particle_sim.geometry.contact import DistanceCache
from particle_sim.geometry.contacts import build_contact_topology
from particle_sim.geometry.transforms import apply_scenario_to_masks, transform_mask
from particle_sim.io.fab import load_fab_table
from particle_sim.io.masks import GridSpec, MaskSet, load_product_masks
from particle_sim.io.stack import build_stack_z
from particle_sim.pipeline import run_nominal, run_probability
from particle_sim.schema import OpenPath
from particle_sim.uncertainty.sampler import Scenario
from particle_sim.units import evaluate_units, unit_specs

SHAPE = (60, 80)
PIXEL_UM = 0.5
STACK = [
    ("M1", "conductor", 0.3, None, True),
    ("GI", "inorganic", 0.4, 0.0, False),
    ("M2", "conductor", 0.3, None, True),
    ("PLN", "organic", 2.0, 0.8, False),
    ("M3", "conductor", 0.1, None, True),
]


def _rect(y0: int, y1: int, x0: int, x1: int) -> np.ndarray:
    mask = np.zeros(SHAPE, dtype=bool)
    mask[y0:y1, x0:x1] = True
    return mask


def _save(path: Path, mask: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(mask.astype(np.uint8) * 255).save(path)


def _contact_case(
    tmp_path: Path,
    *,
    wire: np.ndarray | None = None,
    holes: np.ndarray | None = None,
    pad_net: str = "N1",
    extra: dict | None = None,
) -> Path:
    """M2 wire W joined to two M1 pads (same net) through GI contact holes."""
    wire = _rect(28, 32, 10, 70) if wire is None else wire
    pad_a = _rect(24, 36, 8, 18)
    pad_b = _rect(24, 36, 62, 72)
    holes = (_rect(28, 32, 11, 15) | _rect(28, 32, 65, 69)) if holes is None else holes
    electrodes = {"W": ("M2", "N1", wire), "PA": ("M1", pad_net, pad_a), "PB": ("M1", "N1", pad_b)}
    for name, (_layer, _net, mask) in electrodes.items():
        _save(tmp_path / f"masks/electrode_{name}.png", mask)
    _save(tmp_path / "masks/layer_M1.png", pad_a | pad_b)
    _save(tmp_path / "masks/layer_M2.png", wire)
    _save(tmp_path / "masks/layer_M3.png", np.zeros(SHAPE, dtype=bool))
    _save(tmp_path / "masks/contact_GI.png", holes)
    manifest = {
        "schema_version": "2.0",
        "product_id": "SYN_CONTACT",
        "product_family": "SYN",
        "design_version": "V01",
        "design_lineage_id": "SYN_L1",
        "fab_id": "FAB_SYN",
        "is_synthetic": True,
        "coordinate_system": {
            "origin": "upper_left",
            "x_direction": "right",
            "y_direction": "down",
            "pixel_size_x_um": PIXEL_UM,
            "pixel_size_y_um": PIXEL_UM,
            "roi_width_um": SHAPE[1] * PIXEL_UM,
            "roi_height_um": SHAPE[0] * PIXEL_UM,
        },
        "stack": [
            {
                "layer_name": layer,
                "process_order": order,
                "material_type": material,
                "thickness_um": thickness,
                **({} if factor is None else {"planarization_factor": factor}),
                "patterned": patterned,
            }
            for order, (layer, material, thickness, factor, patterned) in enumerate(STACK)
        ],
        "layer_masks": [
            {"layer_name": layer, "path": f"masks/layer_{layer}.png"}
            for layer in ("M1", "M2", "M3")
        ],
        "electrodes": [
            {
                "electrode_id": name,
                "layer_name": layer,
                "net_id": net,
                "electrode_role": "DATA",
                "path": f"masks/electrode_{name}.png",
            }
            for name, (layer, net, _mask) in electrodes.items()
        ],
        "particle_source_steps": [
            {
                "step_id": step,
                "lands_after_layer": layer,
                "surface_model": model,
                "fab_distribution_layer": step,
            }
            for step, layer, model in (
                ("S_M1", "M1", "conformal"),
                ("S_M2", "M2", "conformal"),
                ("S_PLN", "PLN", "flat"),
            )
        ],
        "electrode_pairs": [],
        "contacts": [
            {
                "contact_id": "CNT",
                "path": "masks/contact_GI.png",
                "hole_layer": "GI",
                "applicable_particle_source_steps": ["S_M1"],
            }
        ],
        "open_applicable_steps_by_layer": {"M2": ["S_M2"], "M1": []},
        "particle_sizes_um": [2.0, 4.0],
        "embedding_depth_um": 0.0,
        "contact_margin_um": 0.0,
        "post_step_policy": {"mode": "protrusion", "coverage_break_factor": 1.5},
    }
    if extra:
        manifest.update(extra)
    path = tmp_path / "manifest.yaml"
    path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    return path


def _load(path: Path):
    manifest = load_manifest(path)
    return manifest, load_product_masks(manifest, path)


# ---------------------------------------------------------------- review defects 1 and 2


def _grid() -> GridSpec:
    return GridSpec(shape_yx=SHAPE, pixel_size_y_um=PIXEL_UM, pixel_size_x_um=PIXEL_UM)


def _single_masks(electrode: np.ndarray, terminal_a: np.ndarray, terminal_b: np.ndarray):
    return MaskSet(
        grid=_grid(),
        layers={},
        electrodes={"E": electrode},
        input_hashes={},
        terminal_regions={"P": (terminal_a, terminal_b)},
    )


class _Manifest:
    product_id = "SYN"
    contact_margin_um = 0.0
    open_paths = [
        OpenPath(path_id="P", electrode="E", terminal_a_mask_path="a", terminal_b_mask_path="b")
    ]


def test_already_disconnected_path_is_not_blamed_on_particles() -> None:
    islands = _rect(28, 32, 5, 25) | _rect(28, 32, 35, 55)
    masks = _single_masks(islands, _rect(28, 32, 5, 7), _rect(28, 32, 53, 55))
    result = simulate_open(_Manifest(), masks, _Manifest.open_paths[0], 4.0)
    assert result.classification == BASELINE_OPEN
    assert not result.failure.any()
    with pytest.raises(InputContractError, match="이물 없이도"):
        check_open_baselines(_Manifest(), masks)


def test_terminal_lost_under_cd_shrink_is_classified_not_counted() -> None:
    wire = _rect(28, 32, 5, 55)
    terminal = np.zeros(SHAPE, dtype=bool)
    terminal[28:32, 5] = True  # 1-px terminal on the electrode edge
    shrunk = transform_mask(wire, _grid(), overlay_xy_um=(0.1, 0.0), cd_bias_um=-0.3)
    moved = transform_mask(terminal, _grid(), overlay_xy_um=(0.1, 0.0))
    masks = _single_masks(shrunk, moved, _rect(28, 32, 53, 55))
    result = simulate_open(_Manifest(), masks, _Manifest.open_paths[0], 4.0)
    assert result.classification == TERMINAL_LOST
    assert not result.failure.any()


# ---------------------------------------------------------------- schema and topology


def test_contact_layers_resolve_from_stack_and_units_are_listed(tmp_path: Path) -> None:
    manifest, masks = _load(_contact_case(tmp_path))
    topology = masks.contact_topology
    assert topology is not None and topology.lost_holes == 0
    assert list(topology.links) == ["cno:CNT:PA__W", "cno:CNT:PB__W"]
    specs = unit_specs(manifest, masks)
    assert [(s.unit_id, s.unit_kind) for s in specs] == [
        ("auto_open:W", "open_electrode"),
        ("cno:CNT:PA__W", "contact"),
        ("cno:CNT:PB__W", "contact"),
    ]


def test_contact_joining_different_nets_is_rejected(tmp_path: Path) -> None:
    excluded = [
        {"electrode_a": a, "electrode_b": b, "reason": "contact net test"}
        for a, b in (("PA", "PB"), ("PA", "W"))
    ]
    path = _contact_case(tmp_path, pad_net="N9", extra={"excluded_electrode_pairs": excluded})
    manifest = load_manifest(path)
    with pytest.raises(InputContractError, match="서로 다른 net"):
        load_product_masks(manifest, path)


def test_contact_hole_touching_no_electrode_is_rejected(tmp_path: Path) -> None:
    holes = _rect(28, 32, 11, 15) | _rect(5, 8, 40, 44)
    path = _contact_case(tmp_path, holes=holes)
    manifest = load_manifest(path)
    with pytest.raises(InputContractError, match="아래 전극 0개"):
        load_product_masks(manifest, path)


def test_contact_schema_requirements(tmp_path: Path) -> None:
    path = _contact_case(tmp_path)
    data = yaml.safe_load(path.read_text(encoding="utf-8"))

    missing_layer = deepcopy(data)
    missing_layer["open_applicable_steps_by_layer"] = {"M2": ["S_M2"]}
    with pytest.raises(ValidationError, match="open_applicable_steps_by_layer에 선언"):
        type(load_manifest(path)).model_validate(missing_layer)

    no_steps = deepcopy(data)
    no_steps["contacts"][0]["applicable_particle_source_steps"] = []
    with pytest.raises(ValidationError, match="applicable_particle_source_steps가 필요"):
        type(load_manifest(path)).model_validate(no_steps)

    v1 = deepcopy(data)
    v1["schema_version"] = "1.0"
    with pytest.raises(ValidationError, match="schema 2.0에서만"):
        type(load_manifest(path)).model_validate(v1)

    reserved = deepcopy(data)
    reserved["open_paths"] = [
        {
            "path_id": "cno:x",
            "electrode": "W",
            "fallback": "roi_boundary",
            "applicable_particle_source_steps": ["S_M2"],
        }
    ]
    with pytest.raises(ValidationError, match="예약 접두사"):
        type(load_manifest(path)).model_validate(reserved)


# ---------------------------------------------------------------- automatic-terminal open


def test_auto_open_matches_explicit_terminals_except_blocked_contacts(tmp_path: Path) -> None:
    manifest, masks = _load(_contact_case(tmp_path))
    topology = masks.contact_topology
    assert topology is not None
    hole_regions = []
    for hole in topology.holes[:2]:
        region = np.zeros(SHAPE, dtype=bool)
        region[hole.ys, hole.xs] = True
        hole_regions.append(region)
    explicit = MaskSet(
        grid=masks.grid,
        layers=masks.layers,
        electrodes=masks.electrodes,
        input_hashes={},
        terminal_regions={"P": (hole_regions[0], hole_regions[1])},
    )
    path = OpenPath(path_id="P", electrode="W", terminal_a_mask_path="a", terminal_b_mask_path="b")
    for diameter in (2.0, 4.0, 6.0):
        auto = simulate_auto_open(masks, "W", diameter, 0.0)
        reference = simulate_open(manifest, explicit, path, diameter)
        assert auto.classification == AUTO_CONFIRMED
        cno = np.zeros(SHAPE, dtype=bool)
        for link in topology.links.values():
            cno |= simulate_cno(link, link.unit_id, masks.grid, diameter, 0.0).failure
        # Auto-open uses the original effective-terminal policy; CNO is a separate map.
        assert not (auto.failure & ~reference.failure).any()
        assert not ((reference.failure & ~auto.failure) & ~cno).any()
        if diameter >= 4.0:
            assert auto.failure.any()


def test_auto_open_ignores_fragments_without_terminals(tmp_path: Path) -> None:
    # a stub branch without a contact hangs off the wire; cutting it off is not an open
    wire = _rect(28, 32, 10, 70) | _rect(8, 28, 38, 42)
    manifest, masks = _load(_contact_case(tmp_path, wire=wire))
    result = simulate_auto_open(masks, "W", 4.0, 0.0)
    stub_cut = (18, 40)  # a 4 um disc here severs the stub only
    assert not result.failure[stub_cut]
    assert result.failure[30, 40]  # the main wire between the two contacts is cut


def test_roi_boundary_is_a_terminal(tmp_path: Path) -> None:
    # the wire leaves the ROI on the left; one contact on the right is enough for an open unit
    wire = _rect(28, 32, 0, 70)
    holes = _rect(28, 32, 65, 69)
    manifest, masks = _load(_contact_case(tmp_path, wire=wire, holes=holes))
    result = simulate_auto_open(masks, "W", 4.0, 0.0)
    assert result.classification == AUTO_CONFIRMED
    assert set(result.terminals) == {"CNT#1", "roi_boundary#1"}
    assert result.failure[30, 30]


def test_electrode_with_one_terminal_is_reported_not_evaluated(tmp_path: Path) -> None:
    holes = _rect(28, 32, 65, 69)
    manifest, masks = _load(_contact_case(tmp_path, holes=holes))
    result = simulate_auto_open(masks, "W", 4.0, 0.0)
    assert result.classification == INSUFFICIENT_TERMINALS
    assert not result.failure.any()


# ---------------------------------------------------------------- CNO


def test_sparse_cno_matches_correlate_randomized_and_clipped_edges() -> None:
    rng = np.random.default_rng(20260925)
    cases = [
        ((0, 0), (3, 5), (8, 9)),
        ((7, 8), (4, 2), (8, 9)),
        ((0, 8), (2, 4), (8, 9)),
    ]
    for _ in range(100):
        height, width = (int(value) for value in rng.integers(1, 9, size=2))
        stamp_shape = (height, width)
        shape = (int(rng.integers(1, 12)), int(rng.integers(1, 13)))
        n_hole = int(rng.integers(1, min(10, shape[0] * shape[1]) + 1))
        points = rng.choice(shape[0] * shape[1], size=n_hole, replace=False)
        ys, xs = np.divmod(points, shape[1])
        stamp = rng.integers(0, 2, size=stamp_shape, dtype=np.int8).astype(bool)
        if not stamp.any():
            stamp[height // 2, width // 2] = True
        rows, cols, actual = hole_blocked_centres(ys, xs, shape, stamp)
        top, bottom = rows.start, rows.stop
        left, right = cols.start, cols.stop
        hole = np.zeros((bottom - top, right - left), dtype=np.int32)
        hole[ys - top, xs - left] = 1
        expected = 2 * correlate(hole, stamp.astype(np.int32), mode="constant", cval=0) >= len(ys)
        np.testing.assert_array_equal(actual, expected)

    # Explicit odd/even stamp and ROI corner placements exercise the clipped
    # translated runs on every side of the candidate-centre window.
    for (y, x), stamp_shape, shape in cases[:3]:
        ys = np.array([y], dtype=int)
        xs = np.array([x], dtype=int)
        stamp = rng.integers(0, 2, size=stamp_shape, dtype=np.int8).astype(bool)
        rows, cols, actual = hole_blocked_centres(ys, xs, shape, stamp)
        top, bottom = rows.start, rows.stop
        left, right = cols.start, cols.stop
        hole = np.zeros((bottom - top, right - left), dtype=np.int32)
        hole[ys - top, xs - left] = 1
        expected = 2 * correlate(hole, stamp.astype(np.int32), mode="constant", cval=0) >= len(ys)
        np.testing.assert_array_equal(actual, expected)


def test_diagonal_cut_pixels_remain_separate_four_connected_holes(tmp_path: Path) -> None:
    holes = np.zeros(SHAPE, dtype=bool)
    holes[28, 11] = True
    holes[29, 12] = True
    _manifest, masks = _load(_contact_case(tmp_path, holes=holes))
    link = masks.contact_topology.links["cno:CNT:PA__W"]
    assert sorted(hole.cut_pixels for hole in link.holes) == [1, 1]


def test_cno_half_cut_hole_matches_brute_force(tmp_path: Path) -> None:
    manifest, masks = _load(_contact_case(tmp_path))
    link = masks.contact_topology.links["cno:CNT:PA__W"]
    hole = link.holes[0]
    grid = masks.grid
    for diameter in (1.0, 2.0, 3.0, 4.0):
        result = simulate_cno(link, link.unit_id, grid, diameter, 0.0)
        radius = diameter / 2.0
        brute = np.zeros(SHAPE, dtype=bool)
        stamp = particle_footprint(grid, radius, 0.0)
        half = stamp.shape[0] // 2
        for y in range(SHAPE[0]):
            for x in range(SHAPE[1]):
                dy = (hole.cut_ys - y) + half
                dx = (hole.cut_xs - x) + half
                inside = (dy >= 0) & (dy < stamp.shape[0]) & (dx >= 0) & (dx < stamp.shape[1])
                brute[y, x] = 2 * int(np.count_nonzero(stamp[dy[inside], dx[inside]])) >= hole.cut_pixels
        np.testing.assert_array_equal(result.failure, brute)
    # A 2x2 um hole cannot be half-covered by a 1 um particle.
    assert not simulate_cno(link, link.unit_id, grid, 1.0, 0.0).failure.any()
    assert simulate_cno(link, link.unit_id, grid, 4.0, 0.0).failure.any()


@pytest.mark.parametrize("covered,expected", [(49, False), (50, True), (51, True)])
def test_even_cut_hole_half_threshold(covered: int, expected: bool) -> None:
    ys, xs = np.nonzero(np.ones((10, 10), dtype=bool))
    ys, xs = ys + 2, xs + 2
    stamp = np.zeros((10, 10), dtype=bool)
    stamp.flat[:covered] = True
    rows, cols, blocked = hole_blocked_centres(ys, xs, (14, 14), stamp)
    assert bool(blocked[7 - rows.start, 7 - cols.start]) is expected


@pytest.mark.parametrize("covered,expected", [(112, False), (113, True)])
def test_odd_225_pixel_cut_hole_threshold(covered: int, expected: bool) -> None:
    ys, xs = np.nonzero(np.ones((15, 15), dtype=bool))
    ys, xs = ys + 2, xs + 2
    stamp = np.zeros((15, 15), dtype=bool)
    stamp.flat[:covered] = True
    rows, cols, blocked = hole_blocked_centres(ys, xs, (19, 19), stamp)
    assert bool(blocked[9 - rows.start, 9 - cols.start]) is expected


def test_auto_open_bytes_match_aae492a_baseline(tmp_path: Path) -> None:
    _manifest, masks = _load(_contact_case(tmp_path))
    payload = b"".join(
        simulate_auto_open(masks, "W", diameter, 0.5).failure.tobytes()
        for diameter in (2.0, 4.0, 6.0)
    )
    assert hashlib.sha256(payload).hexdigest() == (
        "1ef5579d0ea6dd445511be26e47951de70edbd4eb64f7afb9f8b5e39c382c335"
    )


def test_cno_uses_whole_cut_component_not_effective_endpoint(tmp_path: Path) -> None:
    wire = _rect(28, 30, 11, 13)
    hole_mask = _rect(28, 32, 11, 15)
    _manifest, masks = _load(_contact_case(tmp_path, wire=wire, holes=hole_mask))
    link = masks.contact_topology.links["cno:CNT:PA__W"]
    hole = link.holes[0]
    assert (hole.raw_pixels, hole.cut_pixels, hole.effective_pixels) == (None, 16, 4)
    stamp = particle_footprint(masks.grid, 0.5, 0.0)
    rows, cols, effective_blocked = hole_blocked_centres(hole.ys, hole.xs, SHAPE, stamp)
    assert effective_blocked[29 - rows.start, 12 - cols.start]
    assert not simulate_cno(link, link.unit_id, masks.grid, 1.0, 0.0).failure[29, 12]


def test_two_hole_link_ors_per_hole_results(tmp_path: Path) -> None:
    holes = _rect(28, 30, 10, 12) | _rect(28, 30, 15, 17)
    _manifest, masks = _load(_contact_case(tmp_path, holes=holes))
    link = masks.contact_topology.links["cno:CNT:PA__W"]
    assert len(link.holes) == 2
    stamp = particle_footprint(masks.grid, 0.5, 0.0)
    expected = np.zeros(SHAPE, dtype=bool)
    for hole in link.holes:
        rows, cols, blocked = hole_blocked_centres(hole.cut_ys, hole.cut_xs, SHAPE, stamp)
        expected[rows, cols] |= blocked
    result = simulate_cno(link, link.unit_id, masks.grid, 1.0, 0.0)
    np.testing.assert_array_equal(result.failure, expected)
    assert result.failure[29, 11] and result.failure[29, 16]


def test_raw_cut_counts_partial_cut_and_roi_edge(tmp_path: Path) -> None:
    holes = _rect(28, 32, 11, 15)
    manifest, masks = _load(_contact_case(tmp_path, holes=holes))
    raw = holes | _rect(27, 28, 11, 15)
    topology = build_contact_topology(
        manifest, {"CNT": holes}, masks.electrodes, strict=True,
        raw_contact_masks={"CNT": raw},
    )
    hole = topology.holes[0]
    assert (hole.raw_pixels, hole.cut_pixels, hole.effective_pixels) == (20, 16, 16)
    edge_ys = np.array([0, 0, 1, 1])
    edge_xs = np.array([0, 1, 0, 1])
    rows, cols, blocked = hole_blocked_centres(edge_ys, edge_xs, (3, 3), np.ones((3, 3), dtype=bool))
    assert blocked[0 - rows.start, 0 - cols.start]
    with pytest.raises(ValueError, match="nonempty"):
        hole_blocked_centres(np.array([], dtype=int), np.array([], dtype=int), (3, 3), np.ones((1, 1)))
    with pytest.raises(InputContractError, match="mask가 비어"):
        build_contact_topology(manifest, {"CNT": np.zeros(SHAPE, dtype=bool)}, masks.electrodes, strict=True)
    with pytest.raises(InputContractError, match="ROI edge"):
        _load(_contact_case(tmp_path / "edge", holes=_rect(0, 2, 0, 2)))


def test_cno_margin_is_zero_even_when_manifest_margin_is_half_um(tmp_path: Path) -> None:
    manifest, masks = _load(_contact_case(tmp_path, extra={"contact_margin_um": 0.5}))
    link = masks.contact_topology.links["cno:CNT:PA__W"]
    zero = simulate_cno(link, link.unit_id, masks.grid, 1.0, 0.0)
    nominal = simulate_cno(link, link.unit_id, masks.grid, 1.0, manifest.contact_margin_um)
    np.testing.assert_array_equal(nominal.failure, zero.failure)
    evaluated = {r.pair_id: r for r in evaluate_units(
        manifest, masks, build_stack_z(manifest), "S_M1", 1.0, DistanceCache(masks)
    )}
    np.testing.assert_array_equal(evaluated[link.unit_id].failure, zero.failure)
    stamp_with_margin = particle_footprint(masks.grid, 0.5, 0.5)
    rows, cols, inflated = hole_blocked_centres(link.holes[0].cut_ys, link.holes[0].cut_xs, SHAPE, stamp_with_margin)
    assert np.any(inflated & ~zero.failure[rows, cols])


def test_cno_link_missing_in_scenario_is_contact_lost(tmp_path: Path) -> None:
    _manifest, masks = _load(_contact_case(tmp_path))
    result = simulate_cno(None, "cno:CNT:PA__W", masks.grid, 4.0, 0.0)
    assert result.classification == CONTACT_LOST
    assert not result.failure.any()


def test_evaluate_units_applies_source_steps_and_allocation_closes(tmp_path: Path) -> None:
    manifest, masks = _load(_contact_case(tmp_path))
    stack_z = build_stack_z(manifest)
    cache = DistanceCache(masks)
    by_step = {
        step: {r.pair_id: r for r in evaluate_units(manifest, masks, stack_z, step, 4.0, cache)}
        for step in ("S_M1", "S_M2")
    }
    assert by_step["S_M1"]["auto_open:W"].classification == "not_applicable_source_step"
    assert by_step["S_M1"]["cno:CNT:PA__W"].failure.any()
    assert by_step["S_M2"]["cno:CNT:PA__W"].classification == "not_applicable_source_step"
    assert by_step["S_M2"]["auto_open:W"].failure.any()
    for results in by_step.values():
        listed = list(results.values())
        allocation = allocate_risk(listed, masks.grid)
        union = union_by_mode(listed)[ANY_MODE].any_failure
        assert sum(allocation.allocated_ca_um2.values()) == pytest.approx(
            float(union.sum()) * masks.grid.pixel_area_um2
        )
        assert CNO_MODE in {r.failure_mode for r in listed}


# ---------------------------------------------------------------- pipeline and outputs


def test_nominal_run_writes_cno_and_auto_open_units(tmp_path: Path) -> None:
    path = _contact_case(tmp_path)
    run = run_nominal(path, tmp_path / "run")
    ca = pd.read_parquet(run.ca_table)
    assert set(ca["unit_kind"]) == {"open_electrode", "contact"}
    cno = ca[ca["unit_kind"] == "contact"]
    assert set(cno["contact_id"]) == {"cno:CNT:PA__W", "cno:CNT:PB__W"}
    assert set(cno["electrode_b"]) == {"W"}
    union = pd.read_parquet(run.union_table)
    assert set(union["failure_mode"]) == {"lateral_open", CNO_MODE, ANY_MODE}
    for (step, size), group in ca.groupby(["source_step", "particle_size_um"]):
        any_row = union[
            (union["source_step"] == step)
            & (union["particle_size_um"] == size)
            & (union["failure_mode"] == ANY_MODE)
        ]
        assert group["ca_allocated_nominal_um2"].sum() == pytest.approx(
            float(any_row["ca_union_nominal_um2"].iloc[0])
        )
    run_manifest = json.loads(run.run_manifest.read_text(encoding="utf-8"))
    assert run_manifest["cno_coverage_rule"] == "at_least_half_cut_hole_any_hole"
    assert any("CNO" in warning for warning in run_manifest["warnings"])

    fab = tmp_path / "fab.csv"
    fab.write_text(
        "fab_id,product_family,fab_distribution_layer,size_bin_left_um,size_bin_right_um,"
        "size_representative_um,mean_particle_count,inspection_area_per_sample_mm2,"
        "measurement_period,particle_type\n"
        + "".join(
            f"FAB_SYN,SYN,{step},{left},{left + 2.0},{left + 1.0},2.0,10.0,2026Q3,any\n"
            for step in ("S_M1", "S_M2", "S_PLN")
            for left in (1.0, 3.0)
        ),
        encoding="utf-8",
    )
    tables = build_risk_tables(
        ca,
        union,
        load_fab_table(fab, bin_method="exact_point"),
        load_manifest(path),
        fab_id="FAB_SYN",
        weighting_mode="actual",
    )
    assert set(tables.by_contact["contact_link_id"]) == {"cno:CNT:PA__W", "cno:CNT:PB__W"}
    electrode_total = tables.by_electrode["allocated_electrode_contribution"].sum()
    union_total = tables.by_layer[tables.by_layer["layer_scope"] == "particle_source_layer"][
        "geometric_hit_expectation"
    ].sum()
    assert electrode_total == pytest.approx(union_total)
    assert "contact" in set(tables.summary["aggregation_level"])


def test_legacy_compare_marks_new_units_and_not_applicable(tmp_path: Path) -> None:
    run = run_nominal(_contact_case(tmp_path), tmp_path / "run")
    config = tmp_path / "legacy" / "reference.yaml"
    config.parent.mkdir()
    config.write_text(
        yaml.safe_dump(
            {
                "schema_version": "2.0",
                "product_id": "SYN_CONTACT",
                "design_version": "V01",
                "input_hashes": {
                    name: digest
                    for name, digest in json.loads(
                        run.run_manifest.read_text(encoding="utf-8")
                    )["input_hashes"].items()
                    if name.startswith("masks/")
                },
                "axes": "yx",
                "origin": "upper_left",
                "pixel_size_x_um": PIXEL_UM,
                "pixel_size_y_um": PIXEL_UM,
                "particle_size_definition": "diameter",
                "records": [],
            }
        ),
        encoding="utf-8",
    )
    output = compare_legacy_sim(config, run.output_dir, tmp_path / "cmp")
    frame = pd.read_parquet(output / "comparison.parquet")
    assert set(frame["status"]) == {"new_unit", "not_applicable"}
    assert "new_pair" not in set(frame["status"])


def test_probability_run_moves_contacts_with_hole_layer(tmp_path: Path) -> None:
    path = _contact_case(tmp_path)
    scenario = tmp_path / "scenario.yaml"
    scenario.write_text(
        """
schema_version: "2.0"
assumption_set: contact
evidence: measured
overlay_reference_layer: M1
overlay_sign_convention: positive_x_right_positive_y_down
overlay_um:
  GI:
    x: {mean: 0.0, sigma: 0.2, evidence: assumed}
    y: {mean: 0.0, sigma: 0.2, evidence: assumed}
cd_bias_um:
  M2: {mean: -0.1, sigma: 0.1, evidence: assumed}
""",
        encoding="utf-8",
    )
    result = run_probability(path, scenario, tmp_path / "prob", sobol_exponent_m=3, seed=5)
    ca = pd.read_parquet(result.ca_table)
    assert {"open_electrode", "contact"} == set(ca["unit_kind"])
    applicable = ca[(ca["unit_kind"] == "contact") & (ca["source_step"] == "S_M1")]
    assert (applicable["ca_mean_um2"] > 0).any()
    run_manifest = json.loads(result.run_manifest.read_text(encoding="utf-8"))
    assert "scenario_exclusions" in run_manifest
    assert "scenario_contact_lost_holes" in run_manifest


def test_scenario_topology_counts_lost_holes(tmp_path: Path) -> None:
    manifest, masks = _load(_contact_case(tmp_path))
    scenario = Scenario(
        scenario_id=0,
        thickness_delta_um={},
        overlay_xy_um={"GI": (5.0, 0.0)},  # both holes move off their M1 pads
        cd_bias_um={},
        particle_diameter_delta_um=0.0,
        center_offset_xy_um=(0.0, 0.0),
        embedding_depth_um=None,
        planarization_factor={},
        contact_margin_um=None,
    )
    moved = apply_scenario_to_masks(manifest, masks, scenario)
    assert moved.contact_topology is not None
    assert moved.contact_topology.lost_holes == 2
    stack_z = build_stack_z(manifest)
    results = evaluate_units(
        manifest,
        moved,
        stack_z,
        "S_M1",
        4.0,
        DistanceCache(moved),
        specs=unit_specs(manifest, masks),
    )
    cno = [r for r in results if r.failure_mode == CNO_MODE]
    assert {r.classification for r in cno} == {CONTACT_LOST}
    assert not any(r.failure.any() for r in cno)


def test_packed_scenario_masks_keep_contacts(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from particle_sim import pipeline

    path = _contact_case(tmp_path)
    scenario = tmp_path / "scenario.yaml"
    scenario.write_text(
        """
schema_version: "2.0"
assumption_set: contact
evidence: measured
overlay_reference_layer: M1
overlay_sign_convention: positive_x_right_positive_y_down
overlay_um:
  GI:
    x: {mean: 0.0, sigma: 0.3, evidence: assumed}
    y: {mean: 0.0, sigma: 0.3, evidence: assumed}
""",
        encoding="utf-8",
    )
    memory = run_probability(path, scenario, tmp_path / "memory", sobol_exponent_m=2, seed=9)
    monkeypatch.setattr(pipeline._ScenarioMaskProvider, "MEMORY_LIMIT_BYTES", 0)
    packed = run_probability(path, scenario, tmp_path / "packed", sobol_exponent_m=2, seed=9)
    manifest = json.loads(packed.run_manifest.read_text(encoding="utf-8"))
    assert manifest["scenario_mask_storage"] == "packed_temp"
    columns = ["failure_unit_id", "source_step", "particle_size_um", "ca_mean_um2"]
    pd.testing.assert_frame_equal(
        pd.read_parquet(memory.ca_table)[columns], pd.read_parquet(packed.ca_table)[columns]
    )


def test_structurally_excluded_unit_is_not_reported_as_scenario_exclusion(tmp_path: Path) -> None:
    # W has one contact only, so auto_open:W is insufficient in the nominal geometry already
    path = _contact_case(tmp_path, holes=_rect(28, 32, 65, 69))
    scenario = tmp_path / "scenario.yaml"
    scenario.write_text(
        """
schema_version: "2.0"
assumption_set: contact
evidence: measured
overlay_reference_layer: M1
overlay_sign_convention: positive_x_right_positive_y_down
overlay_um:
  M2:
    x: {mean: 0.0, sigma: 0.2, evidence: assumed}
    y: {mean: 0.0, sigma: 0.2, evidence: assumed}
""",
        encoding="utf-8",
    )
    result = run_probability(path, scenario, tmp_path / "prob", sobol_exponent_m=2, seed=3)
    run_manifest = json.loads(result.run_manifest.read_text(encoding="utf-8"))
    assert run_manifest["scenario_exclusions"] == []
    assert not any("SCENARIO_EXCLUSION" in w for w in run_manifest["warnings"])


@pytest.mark.parametrize("strict", [True, False])
def test_cno_unit_id_collision_is_rejected(tmp_path: Path, strict: bool) -> None:
    path = _contact_case(tmp_path)
    data = yaml.safe_load(path.read_text())
    data["electrodes"] = [
        {"electrode_id": name, "layer_name": layer, "net_id": "N1",
         "electrode_role": "DATA", "path": f"masks/electrode_{name}.png"}
        for name, layer in [("P__Q", "M1"), ("P", "M1"),
                            ("W", "M2"), ("Q__W", "M2")]
    ]
    path.write_text(yaml.safe_dump(data))
    left = _rect(24, 36, 8, 18)
    right = _rect(24, 36, 62, 72)
    electrodes = {"P__Q": left, "P": right, "W": left, "Q__W": right}
    for name, mask in electrodes.items():
        _save(tmp_path / f"masks/electrode_{name}.png", mask)
    _save(tmp_path / "masks/layer_M2.png", left | right)
    manifest = load_manifest(path)
    with pytest.raises(InputContractError, match="CNO unit_id .*중복됩니다") as caught:
        build_contact_topology(
            manifest, {"CNT": _rect(28, 32, 11, 15) | _rect(28, 32, 65, 69)},
            electrodes, strict=strict, manifest_path=str(path),
        )
    assert caught.value.product_id == manifest.product_id
    assert caught.value.file == str(path)
    assert caught.value.field == "contacts[0].path"
    if strict:
        from particle_sim.cli import main

        assert main(["validate-input", "--product-config", str(path)]) == 2


def _scenario_collision_case(tmp_path: Path, *, ambiguous: bool = True) -> Path:
    path = _contact_case(tmp_path, holes=_rect(28, 32, 9, 13) | _rect(28, 32, 38, 42))
    data = yaml.safe_load(path.read_text())
    groups = {
        "M1": {"P__Q": _rect(24, 36, 5, 17), "P": _rect(24, 36, 34, 62)},
        "M2": {"W": _rect(24, 36, 5, 17), "W2": _rect(24, 36, 34, 46),
               "Q__W" if ambiguous else "Q__X": _rect(24, 36, 48, 62)},
    }
    data["electrodes"] = [
        {"electrode_id": name, "layer_name": layer, "net_id": "N1",
         "electrode_role": "DATA", "path": f"masks/electrode_{name}.png"}
        for layer, group in groups.items() for name in group
    ]
    for layer, group in groups.items():
        _save(tmp_path / f"masks/layer_{layer}.png", np.logical_or.reduce(list(group.values())))
        for name, mask in group.items():
            _save(tmp_path / f"masks/electrode_{name}.png", mask)
    path.write_text(yaml.safe_dump(data))
    return path


def test_scenario_only_cno_collision_is_rejected_at_load(tmp_path: Path) -> None:
    from particle_sim.cli import main

    path = _scenario_collision_case(tmp_path)
    manifest = load_manifest(path)
    with pytest.raises(InputContractError, match="CNO unit_id .*중복됩니다") as caught:
        load_product_masks(manifest, path)
    assert caught.value.product_id == manifest.product_id
    assert caught.value.file == str(path)
    assert caught.value.field == "contacts[0].path"
    assert main(["validate-input", "--product-config", str(path)]) == 2


def test_unambiguous_double_underscore_ids_load_and_transform(tmp_path: Path) -> None:
    manifest, masks = _load(_scenario_collision_case(tmp_path, ambiguous=False))
    assert set(masks.contact_topology.links) == {"cno:CNT:P__Q__W", "cno:CNT:P__W2"}
    scenario = Scenario(
        scenario_id=0, thickness_delta_um={}, overlay_xy_um={"GI": (6.0, 0.0)},
        cd_bias_um={}, particle_diameter_delta_um=0.0, center_offset_xy_um=(0.0, 0.0),
        embedding_depth_um=None, planarization_factor={}, contact_margin_um=None,
    )
    moved = apply_scenario_to_masks(manifest, masks, scenario)
    assert set(moved.contact_topology.links) == {"cno:CNT:P__Q__X"}
    assert moved.contact_topology.lost_holes == 1
