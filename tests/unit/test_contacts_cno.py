"""Contact-derived open terminals, CNO and the open baseline checks (2026-09-19 review)."""

from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml
from PIL import Image
from pydantic import ValidationError

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
from particle_sim.failure.contact_not_open import CNO_MODE, CONTACT_LOST, simulate_cno
from particle_sim.failure.open_failure import (
    BASELINE_OPEN,
    TERMINAL_LOST,
    check_open_baselines,
    particle_footprint,
    simulate_open,
)
from particle_sim.failure.union import ANY_MODE, union_by_mode
from particle_sim.geometry.contact import DistanceCache
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
        # only positions where a contact is fully blocked move from open to CNO
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


def test_cno_requires_full_cover_and_matches_brute_force(tmp_path: Path) -> None:
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
                dy = (hole.ys - y) + half
                dx = (hole.xs - x) + half
                inside = (dy >= 0) & (dy < stamp.shape[0]) & (dx >= 0) & (dx < stamp.shape[1])
                brute[y, x] = bool(inside.all() and stamp[dy, dx].all())
        np.testing.assert_array_equal(result.failure, brute)
    # a 2x2 um hole cannot be fully covered by a 1 um particle
    assert not simulate_cno(link, link.unit_id, grid, 1.0, 0.0).failure.any()
    assert simulate_cno(link, link.unit_id, grid, 4.0, 0.0).failure.any()


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
    assert run_manifest["cno_coverage_rule"] == "full_cover_any_hole"
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
