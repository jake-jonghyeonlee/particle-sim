"""Ownership and compatibility checks for nominal evaluation placeholders."""

from __future__ import annotations

import gc
import weakref
from copy import deepcopy
from datetime import UTC, datetime
from pathlib import Path

import pytest

import particle_sim.pipeline as pipeline
from particle_sim.config import load_manifest
from particle_sim.geometry.contact import DistanceCache
from particle_sim.geometry.transforms import apply_scenario_to_masks
from particle_sim.io.masks import load_product_masks
from particle_sim.io.stack import build_stack_z
from particle_sim.uncertainty.sampler import Scenario
from particle_sim.units import (
    CNO_MODE,
    KIND_CONTACT,
    KIND_OPEN_ELECTRODE,
    NOT_APPLICABLE,
    UnitSpec,
    evaluate_units,
    unit_specs,
)
from particle_sim.validation.synthetic_cases import build_case, write_case


def _case(tmp_path: Path):
    case = build_case("open_straight")
    original = case.manifest["open_paths"][0]
    other = deepcopy(original)
    other["path_id"] = "INAPPLICABLE"
    other["applicable_particle_source_steps"] = [
        case.manifest["particle_source_steps"][-1]["step_id"]
    ]
    case.manifest["open_paths"].append(other)
    path = write_case(case, tmp_path / "input")
    manifest = load_manifest(path)
    masks = load_product_masks(manifest, path)
    return path, manifest, masks


def test_placeholder_bytes_aliases_and_default_mutation(tmp_path: Path) -> None:
    _, manifest, masks = _case(tmp_path)
    args = (manifest, masks, build_stack_z(manifest),
            manifest.particle_source_steps[0].step_id, 6.0, DistanceCache(masks))
    default = evaluate_units(*args)
    pooled = evaluate_units(*args, shared_readonly_placeholders=True)
    assert len(default) == len(pooled) == 2
    for old, new in zip(default, pooled, strict=True):
        for field in ("pair_id", "failure_mode", "electrode_a", "electrode_b", "classification"):
            assert getattr(old, field) == getattr(new, field)
        for field in ("failure", "contact_a", "contact_b", "pair_margin_um"):
            before, after = getattr(old, field), getattr(new, field)
            assert before.shape == after.shape and before.dtype == after.dtype
            assert before.flags.c_contiguous and after.flags.c_contiguous
            assert before.tobytes() == after.tobytes()
    assert pooled[0].contact_b is pooled[1].failure
    assert pooled[0].pair_margin_um is pooled[1].pair_margin_um
    assert pooled[1].contact_a is pooled[1].failure
    for array in (pooled[0].contact_b, pooled[0].pair_margin_um):
        with pytest.raises(ValueError):
            array[0, 0] = 1
        with pytest.raises(ValueError):
            array.setflags(write=True)
    assert default[0].contact_b is not default[1].failure
    default[0].contact_b[0, 0] = True
    assert not default[1].failure[0, 0]
    assert pooled[0].failure.flags.writeable  # actual open geometry is still independent


def test_pool_is_local_and_cache_keeps_geometry(tmp_path: Path) -> None:
    _, manifest, masks = _case(tmp_path)
    cache = DistanceCache(masks)
    open_cache = {}
    args = (manifest, masks, build_stack_z(manifest),
            manifest.particle_source_steps[0].step_id, 6.0, cache, open_cache)
    first = evaluate_units(*args, shared_readonly_placeholders=True)
    second = evaluate_units(*args, shared_readonly_placeholders=True)
    assert first[0].failure is second[0].failure
    assert first[0].contact_b is not second[0].contact_b
    assert second[0].contact_b is second[1].failure
    default_after_opt_in = evaluate_units(*args)
    assert default_after_opt_in[0].contact_b.flags.writeable
    assert default_after_opt_in[0].pair_margin_um.flags.writeable
    assert default_after_opt_in[0].contact_b is not default_after_opt_in[1].failure
    other_case = build_case("open_straight")
    other_shape = (193, 197)
    for maps in (other_case.layer_masks, other_case.electrode_masks):
        for key, value in maps.items():
            maps[key] = value[:other_shape[0], :other_shape[1]].copy()
    coordinates = other_case.manifest["coordinate_system"]
    coordinates["roi_height_um"] = other_shape[0] * coordinates["pixel_size_y_um"]
    coordinates["roi_width_um"] = other_shape[1] * coordinates["pixel_size_x_um"]
    other_shape_path = write_case(other_case, tmp_path / "other")
    other_manifest = load_manifest(other_shape_path)
    other_masks = load_product_masks(other_manifest, other_shape_path)
    other = evaluate_units(other_manifest, other_masks, build_stack_z(other_manifest),
                           other_manifest.particle_source_steps[0].step_id, 6.0,
                           DistanceCache(other_masks), shared_readonly_placeholders=True)
    assert other[0].contact_b is not second[0].contact_b
    assert other[0].contact_b.shape == other_shape


def test_pair_auto_open_and_missing_contact_keep_exact_results(tmp_path: Path) -> None:
    for name in ("pilot", "open_straight"):
        path = write_case(build_case(name), tmp_path / name)
        manifest = load_manifest(path)
        masks = load_product_masks(manifest, path)
        specs = unit_specs(manifest, masks)
        if name == "pilot":
            specs.append(UnitSpec("missing-contact", CNO_MODE, KIND_CONTACT,
                                  "A", "B", None))
        else:
            electrode = manifest.electrodes[0].electrode_id
            specs.append(UnitSpec("auto-test", "lateral_open", KIND_OPEN_ELECTRODE,
                                  electrode, NOT_APPLICABLE, None))
        args = (manifest, masks, build_stack_z(manifest),
                manifest.particle_source_steps[0].step_id, 4.0, DistanceCache(masks))
        ordinary = evaluate_units(*args, specs=specs)
        pooled = evaluate_units(*args, specs=specs, shared_readonly_placeholders=True)
        for old, new in zip(ordinary, pooled, strict=True):
            assert old.classification == new.classification
            for field in ("failure", "contact_a", "contact_b", "pair_margin_um"):
                assert getattr(old, field).tobytes() == getattr(new, field).tobytes()
        if name == "pilot":
            assert pooled[0].failure.flags.writeable
            assert pooled[-1].classification == "contact_lost"


def test_zero_applicable_nonpair_units(tmp_path: Path) -> None:
    case = build_case("open_straight")
    case.manifest["open_paths"][0]["applicable_particle_source_steps"] = [
        case.manifest["particle_source_steps"][-1]["step_id"]
    ]
    path = write_case(case, tmp_path / "input")
    manifest = load_manifest(path)
    masks = load_product_masks(manifest, path)
    result = evaluate_units(manifest, masks, build_stack_z(manifest),
                            manifest.particle_source_steps[0].step_id, 4.0,
                            DistanceCache(masks), shared_readonly_placeholders=True)
    assert len(result) == 1
    assert result[0].classification == "not_applicable_source_step"
    assert result[0].failure is result[0].contact_a is result[0].contact_b


def test_transformed_masks_get_separate_placeholders(tmp_path: Path) -> None:
    _, manifest, masks = _case(tmp_path)
    layer = manifest.electrodes[0].layer_name
    scenario = Scenario(
        scenario_id=0, thickness_delta_um={}, overlay_xy_um={layer: (1.0, 0.0)},
        cd_bias_um={}, particle_diameter_delta_um=0.0,
        center_offset_xy_um=(0.0, 0.0), embedding_depth_um=None,
        planarization_factor={}, contact_margin_um=None,
    )
    moved = apply_scenario_to_masks(manifest, masks, scenario)
    step = manifest.particle_source_steps[0].step_id
    stack = build_stack_z(manifest)
    nominal = evaluate_units(manifest, masks, stack, step, 4.0, DistanceCache(masks),
                             shared_readonly_placeholders=True)
    ordinary = evaluate_units(manifest, moved, stack, step, 4.0, DistanceCache(moved))
    pooled = evaluate_units(manifest, moved, stack, step, 4.0, DistanceCache(moved),
                            shared_readonly_placeholders=True)
    assert nominal[0].contact_b is not pooled[0].contact_b
    for old, new in zip(ordinary, pooled, strict=True):
        assert old.classification == new.classification
        for field in ("failure", "contact_a", "contact_b", "pair_margin_um"):
            assert getattr(old, field).tobytes() == getattr(new, field).tobytes()


def test_nominal_releases_uncached_results_before_next_evaluation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path, _, _ = _case(tmp_path)
    original = pipeline.evaluate_units
    previous: weakref.ReferenceType | None = None
    seen = 0

    def observed(*args, **kwargs):
        nonlocal previous, seen
        gc.collect()
        if previous is not None:
            assert previous() is None
        result = original(*args, **kwargs)
        previous = weakref.ref(result[1])
        seen += 1
        return result

    monkeypatch.setattr(pipeline, "evaluate_units", observed)
    pipeline.run_nominal(path, tmp_path / "run", now=datetime(2026, 9, 29, tzinfo=UTC))
    assert seen > 1
