"""schema-config 단위 테스트."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest
import yaml
from pydantic import ValidationError

from particle_sim.config import canonical_hash, file_sha256, load_manifest, resolve_input_path
from particle_sim.errors import InputContractError
from particle_sim.schema import (
    PARTICLE_SIZES_UM,
    ProductManifest,
    todo_physics_warnings,
)


def base_manifest() -> dict:
    return {
        "schema_version": "1.0",
        "product_id": "P1",
        "product_family": "fam",
        "design_version": "v1",
        "design_lineage_id": "lin1",
        "fab_id": "fab1",
        "is_synthetic": True,
        "coordinate_system": {
            "origin": "upper_left",
            "x_direction": "right",
            "y_direction": "down",
            "pixel_size_x_um": 0.1,
            "pixel_size_y_um": 0.1,
            "roi_width_um": 100.0,
            "roi_height_um": 100.0,
        },
        "stack": [
            {"layer_name": "L1", "process_order": 0, "material_type": "conductor",
             "thickness_um": 1.0, "patterned": True},
            {"layer_name": "L2", "process_order": 1, "material_type": "inorganic",
             "thickness_um": 2.0, "planarization_factor": 0.5, "patterned": False},
            {"layer_name": "L3", "process_order": 2, "material_type": "conductor",
             "thickness_um": 1.5, "patterned": True},
        ],
        "layer_masks": [
            {"layer_name": "L1", "path": "masks/L1.png"},
            {"layer_name": "L3", "path": "masks/L3.png"},
        ],
        "electrodes": [
            {"electrode_id": "E1", "layer_name": "L1", "net_id": "A",
             "electrode_role": "source", "path": "masks/E1.png"},
            {"electrode_id": "E2", "layer_name": "L1", "net_id": "B",
             "electrode_role": "drain", "path": "masks/E2.png"},
            {"electrode_id": "E3", "layer_name": "L3", "net_id": "C",
             "electrode_role": "gate", "path": "masks/E3.png"},
        ],
        "particle_source_steps": [
            {"step_id": "S1", "lands_after_layer": "L2", "surface_model": "conformal",
             "fab_distribution_layer": "L1"},
        ],
        "electrode_pairs": [
            {"pair_id": "P1", "electrode_a": "E1", "electrode_b": "E2",
             "failure_mode": "lateral_short"},
            {"pair_id": "P2", "electrode_a": "E1", "electrode_b": "E3",
             "failure_mode": "vertical_short"},
        ],
        "embedding_depth_um": 0.5,
        "contact_margin_um": 0.1,
        "post_step_policy": {"mode": "exclude"},
    }


def write_yaml(tmp_path: Path, data: dict, name: str = "manifest.yaml") -> Path:
    p = tmp_path / name
    p.write_text(yaml.safe_dump(data), encoding="utf-8")
    return p


def test_valid_manifest_roundtrip() -> None:
    m = ProductManifest.model_validate(base_manifest())
    assert m.product_id == "P1"
    assert len(m.stack) == 3
    assert m.particle_sizes_um == list(PARTICLE_SIZES_UM)
    assert m.legacy_kernel.size == 5


def test_particle_sizes() -> None:
    assert PARTICLE_SIZES_UM[:3] == (1.0, 2.0, 3.0)
    assert PARTICLE_SIZES_UM[-2:] == (40.0, 50.0)
    assert len(PARTICLE_SIZES_UM) == 32


@pytest.mark.parametrize(
    "name,mutate",
    [
        ("dup_layer_name", lambda d: d["stack"].append(
            {"layer_name": "L1", "process_order": 9, "material_type": "other",
             "thickness_um": 1.0, "patterned": False})),
        ("dup_process_order", lambda d: d["stack"].append(
            {"layer_name": "L9", "process_order": 0, "material_type": "other",
             "thickness_um": 1.0, "patterned": False})),
        ("dup_electrode_id", lambda d: d["electrodes"].append(
            {"electrode_id": "E1", "layer_name": "L1", "net_id": "Z",
             "electrode_role": "x", "path": "masks/E9.png"})),
        ("dup_step_id", lambda d: d["particle_source_steps"].append(
            {"step_id": "S1", "lands_after_layer": "L2", "surface_model": "flat",
             "fab_distribution_layer": "L1"})),
        ("dup_pair_id", lambda d: d["electrode_pairs"].append(
            {"pair_id": "P1", "electrode_a": "E2", "electrode_b": "E3",
             "failure_mode": "vertical_short"})),
        ("mask_unknown_layer", lambda d: d["layer_masks"].append(
            {"layer_name": "LX", "path": "masks/LX.png"})),
        ("electrode_unknown_layer", lambda d: d["electrodes"][0].update(layer_name="LX")),
        ("patterned_no_mask", lambda d: d["layer_masks"].remove(
            next(x for x in d["layer_masks"] if x["layer_name"] == "L1"))),
        ("patterned_two_masks", lambda d: d["layer_masks"].append(
            {"layer_name": "L1", "path": "masks/L1b.png"})),
        ("non_patterned_mask", lambda d: d["layer_masks"].append(
            {"layer_name": "L2", "path": "masks/L2.png"})),
        ("electrode_non_patterned", lambda d: d["electrodes"][0].update(layer_name="L2")),
        ("electrode_non_conductor", lambda d: d["stack"][0].update(material_type="inorganic")),
        ("step_lands_unknown", lambda d: d["particle_source_steps"][0].update(
            lands_after_layer="LX")),
        ("pair_a_missing", lambda d: d["electrode_pairs"][0].update(electrode_a="EX")),
        ("pair_b_missing", lambda d: d["electrode_pairs"][0].update(electrode_b="EX")),
        ("pair_same_electrode", lambda d: d["electrode_pairs"][0].update(electrode_b="E1")),
        ("pair_same_net", lambda d: d["electrodes"][1].update(net_id="A")),
        ("lateral_diff_layer", lambda d: d["electrode_pairs"][0].update(electrode_b="E3")),
        ("vertical_same_layer", lambda d: d["electrode_pairs"][1].update(
            electrode_b="E2", failure_mode="vertical_short")),
    ],
)
def test_cross_reference_violations(name: str, mutate) -> None:
    d = base_manifest()
    mutate(d)
    with pytest.raises(ValidationError):
        ProductManifest.model_validate(d)


def test_legacy_kernel_even_size() -> None:
    d = base_manifest()
    d["legacy_kernel"] = {"size": 4}
    with pytest.raises(ValidationError):
        ProductManifest.model_validate(d)


def test_protrusion_without_coverage() -> None:
    d = base_manifest()
    d["post_step_policy"] = {"mode": "protrusion"}
    with pytest.raises(ValidationError):
        ProductManifest.model_validate(d)


def test_partial_without_planarization() -> None:
    d = base_manifest()
    d["stack"][1]["planarization_factor"] = None
    d["particle_source_steps"][0]["surface_model"] = "partial"
    with pytest.raises(ValidationError):
        ProductManifest.model_validate(d)


def test_protrusion_without_planarization() -> None:
    d = base_manifest()
    d["stack"][1]["planarization_factor"] = None
    d["post_step_policy"] = {"mode": "protrusion", "coverage_break_factor": 0.5}
    with pytest.raises(ValidationError):
        ProductManifest.model_validate(d)


def test_load_manifest_happy(tmp_path: Path) -> None:
    p = write_yaml(tmp_path, base_manifest())
    m = load_manifest(p)
    assert m.product_id == "P1"


def test_load_manifest_yaml_error(tmp_path: Path) -> None:
    p = tmp_path / "bad.yaml"
    p.write_text("a: [unclosed", encoding="utf-8")
    with pytest.raises(InputContractError) as ei:
        load_manifest(p)
    assert "YAML" in str(ei.value)


def test_load_manifest_validation_error(tmp_path: Path) -> None:
    d = base_manifest()
    d["stack"][0]["layer_name"] = "L1"
    d["stack"].append(
        {"layer_name": "L1", "process_order": 9, "material_type": "other",
         "thickness_um": 1.0, "patterned": False})
    p = write_yaml(tmp_path, d)
    with pytest.raises(InputContractError) as ei:
        load_manifest(p)
    err = ei.value
    assert err.product_id == "P1"
    assert err.file == str(p)
    assert err.field is not None
    assert "[product=P1]" in str(err)
    assert "[file=" in str(err)
    assert "[field=" in str(err)


def test_resolve_input_path_valid(tmp_path: Path) -> None:
    mp = tmp_path / "manifest.yaml"
    mp.write_text("", encoding="utf-8")
    p = resolve_input_path(mp, "masks/L1.png", product_id="P1", field="layer_masks[0].path")
    assert p == (tmp_path / "masks" / "L1.png").resolve()


def test_resolve_input_path_absolute(tmp_path: Path) -> None:
    mp = tmp_path / "manifest.yaml"
    mp.write_text("", encoding="utf-8")
    with pytest.raises(InputContractError) as ei:
        resolve_input_path(mp, "/etc/passwd", product_id="P1", field="x")
    assert "[product=P1]" in str(ei.value)


def test_resolve_input_path_escape(tmp_path: Path) -> None:
    mp = tmp_path / "manifest.yaml"
    mp.write_text("", encoding="utf-8")
    with pytest.raises(InputContractError) as ei:
        resolve_input_path(mp, "../outside.png", product_id="P1", field="x")
    assert "[product=P1]" in str(ei.value)


def test_canonical_hash_deterministic() -> None:
    a = canonical_hash({"b": 1, "a": 2})
    b = canonical_hash({"a": 2, "b": 1})
    assert a == b
    assert a == hashlib.sha256(
        json.dumps({"a": 2, "b": 1}, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def test_file_sha256(tmp_path: Path) -> None:
    p = tmp_path / "f.bin"
    p.write_bytes(b"hello")
    assert file_sha256(p) == hashlib.sha256(b"hello").hexdigest()


def test_todo_physics_warnings() -> None:
    m = ProductManifest.model_validate(base_manifest())
    ws = todo_physics_warnings(m)
    assert any("particle_conductivity" in w for w in ws)
    assert any("embedding_depth_um=0.5" in w for w in ws)
    assert any("contact_margin_um=0.1" in w for w in ws)
    assert any("post_step_policy.mode=exclude" in w for w in ws)
    assert any("L2.planarization_factor=0.5" in w for w in ws)
    assert all(w.startswith("TODO_PHYSICS:") for w in ws)


def test_todo_physics_protrusion_coverage() -> None:
    d = base_manifest()
    d["post_step_policy"] = {"mode": "protrusion", "coverage_break_factor": 0.5}
    m = ProductManifest.model_validate(d)
    ws = todo_physics_warnings(m)
    assert any("coverage_break_factor=0.5" in w for w in ws)


def test_error_str_format() -> None:
    e = InputContractError("msg", product_id="P", file="F", field="X")
    assert str(e) == "[product=P][file=F][field=X] msg"
    assert str(InputContractError("msg")) == "msg"
    assert str(InputContractError("msg", product_id="P")) == "[product=P] msg"


def test_fab_distribution_layer_is_free_identifier() -> None:
    data = base_manifest()
    data["particle_source_steps"][0]["fab_distribution_layer"] = "S_FAB_ONLY"
    m = ProductManifest.model_validate(data)
    assert m.particle_source_steps[0].fab_distribution_layer == "S_FAB_ONLY"


def test_missing_manifest_file(tmp_path: Path) -> None:
    with pytest.raises(InputContractError) as exc:
        load_manifest(tmp_path / "missing.yaml")
    assert "file=" in str(exc.value)
