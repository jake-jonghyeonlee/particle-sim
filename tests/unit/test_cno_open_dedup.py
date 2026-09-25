"""CNO/open relation and per-condition accounting counterexamples."""

from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import zarr

from particle_sim.config import canonical_hash
from particle_sim.errors import InputContractError
from particle_sim.failure.contact_not_open import simulate_cno
from particle_sim.geometry.contacts import ContactHole, ContactLink, ContactTopology
from particle_sim.io.masks import GridSpec, MaskSet
from particle_sim.io.outputs import HotspotStore, open_hotspot_store
from particle_sim.schema import ProductManifest
from particle_sim.units import KIND_CONTACT, KIND_OPEN_ELECTRODE, KIND_OPEN_PATH, UnitSpec
from particle_sim.validation.cno_open_dedup import _relations, analyze_cno_open, split_open
from particle_sim.validation.synthetic_cases import build_case


def _hole(index: int, x: int) -> ContactHole:
    ys, xs = np.mgrid[40:55, x : x + 15]
    return ContactHole(
        "CNT",
        index,
        "A",
        "B",
        ys.ravel(),
        xs.ravel(),
        ys.ravel(),
        xs.ravel(),
        None,
        225,
        225,
    )


def test_parallel_holes_use_any_hole_link_and_same_particle_condition() -> None:
    h1, h2 = _hole(1, 30), _hole(2, 120)
    link = ContactLink("cno:CNT:A__B", "CNT", "A", "B", (h1, h2))
    cno = simulate_cno(link, link.unit_id, GridSpec((180, 180), 0.1, 0.1), 2.0, 0.5)
    assert cno.failure[47, 37] and cno.failure[47, 127]
    raw = np.zeros((180, 180), dtype=bool)
    raw[47, 37] = raw[47, 90] = True
    final, removed = split_open(raw, [cno.failure])
    assert removed[47, 37] and not removed[47, 90]
    assert final[47, 90] and not final[47, 37]
    np.testing.assert_array_equal(final | removed, raw)
    np.testing.assert_array_equal(raw | cno.failure, final | removed | cno.failure)
    # Same centre at a different source step or diameter has no related CNO event.
    other_final, other_removed = split_open(raw, [])
    np.testing.assert_array_equal(other_final, raw)
    assert not other_removed.any()


def _contract():
    link = ContactLink("cno:CNT:A__B", "CNT", "A", "B", (_hole(1, 30),))
    manifest = SimpleNamespace(
        open_paths=[
            SimpleNamespace(
                path_id="OPEN_1",
                electrode="A",
                fallback="error",
                terminal_a_yx_px=(1, 1),
                terminal_b_yx_px=(2, 2),
                terminal_a_mask_path=None,
                terminal_b_mask_path=None,
            )
        ],
        electrodes=[
            SimpleNamespace(electrode_id="A", net_id="NET_1"),
            SimpleNamespace(electrode_id="B", net_id="NET_1"),
        ],
        particle_source_steps=[
            SimpleNamespace(step_id="CNT_STEP"),
            SimpleNamespace(step_id="OTHER_STEP"),
        ],
    )
    masks = SimpleNamespace(
        contact_topology=ContactTopology((link.holes[0],), {link.unit_id: link})
    )
    specs = [
        UnitSpec("OPEN_1", "lateral_open", KIND_OPEN_PATH, "A", "N/A", ("CNT_STEP", "OTHER_STEP")),
        UnitSpec(link.unit_id, "contact_not_open", KIND_CONTACT, "A", "B", ("CNT_STEP",)),
    ]
    sidecar = {
        "circuit": {
            "net_status": "declared",
            "auto_open": "disabled",
            "open_path_circuit_metadata": [
                {
                    "path_id": "OPEN_1",
                    "circuit_connection_id": "CONN_1",
                    "circuit_net_id": "NET_1",
                }
            ],
            "cno_open_relations": [
                {
                    "cno_link_id": link.unit_id,
                    "open_path_id": "OPEN_1",
                    "circuit_connection_id": "CONN_1",
                }
            ],
        }
    }
    return manifest, masks, specs, sidecar


def _prepare_dedup_case(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    data = build_case("parallel_lines").manifest
    data["schema_version"] = "2.0"
    data["coordinate_system"].update(
        pixel_size_x_um=0.1,
        pixel_size_y_um=0.1,
        roi_width_um=120.0,
        roi_height_um=120.0,
    )
    data["particle_source_steps"] = data["particle_source_steps"][:2]
    data["particle_sizes_um"] = [1.0]
    data["electrodes"][0]["layer_name"] = "M1"
    data["electrode_pairs"][0]["failure_mode"] = "vertical_short"
    data["electrodes"].append(
        {
            "electrode_id": "C",
            "layer_name": "M2",
            "net_id": "N1",
            "electrode_role": "PAD",
            "path": "masks/electrode_C.png",
        }
    )
    data["contacts"] = [
        {
            "contact_id": "CNT",
            "path": "masks/contact_CNT.png",
            "hole_layer": "GI",
            "applicable_particle_source_steps": ["S_M1"],
        }
    ]
    data["open_applicable_steps_by_layer"] = {"M1": [], "M2": []}
    data["open_paths"] = [
        {
            "path_id": path_id,
            "electrode": electrode,
            "terminal_a_yx_px": [400, 400],
            "terminal_b_yx_px": [800, 800],
            "applicable_particle_source_steps": ["S_M1", "S_M2"],
            "fallback": "error",
        }
        for path_id, electrode in (("OPEN_1", "A"), ("OPEN_2", "B"))
    ]
    manifest = ProductManifest.model_validate(data)
    product_path = tmp_path / "product.yaml"
    product_path.write_text(json.dumps(manifest.model_dump(mode="json"), sort_keys=True))
    h1, h2 = _hole(1, 30), _hole(2, 120)
    link = ContactLink("cno:CNT:A__C", "CNT", "A", "C", (h1, h2))
    masks = MaskSet(
        GridSpec((1200, 1200), 0.1, 0.1),
        {},
        {},
        {},
        contact_topology=ContactTopology((h1, h2), {link.unit_id: link}),
    )
    monkeypatch.setattr("particle_sim.validation.cno_open_dedup.load_manifest", lambda _: manifest)
    monkeypatch.setattr(
        "particle_sim.validation.cno_open_dedup.load_product_masks", lambda *_: masks
    )
    sidecar = {
        "product_yaml_sha256": hashlib.sha256(product_path.read_bytes()).hexdigest(),
        "circuit": {
            "net_status": "declared",
            "auto_open": "disabled",
            "open_path_circuit_metadata": [
                {"path_id": "OPEN_1", "circuit_connection_id": "CONN_1", "circuit_net_id": "N1"},
                {"path_id": "OPEN_2", "circuit_connection_id": "CONN_2", "circuit_net_id": "N2"},
            ],
            "cno_open_relations": [
                {
                    "cno_link_id": link.unit_id,
                    "open_path_id": "OPEN_1",
                    "circuit_connection_id": "CONN_1",
                }
            ],
        },
    }
    sidecar_path = tmp_path / "adapter-provenance.json"
    sidecar_path.write_text(json.dumps(sidecar))
    run = tmp_path / "raw"
    run.mkdir()
    ids = ["A__B", "OPEN_1", "OPEN_2", link.unit_id]
    HotspotStore(
        run / "hotspots.zarr",
        source_steps=["S_M1", "S_M2"],
        particle_sizes_um=[1.0],
        pair_ids=ids,
        pair_modes=["vertical_short", "lateral_open", "lateral_open", "contact_not_open"],
        unit_kinds=["electrode_pair", "open_path", "open_path", "contact"],
        union_keys=["vertical_short", "lateral_open", "contact_not_open", "any"],
        shape_yx=(1200, 1200),
    )
    group = zarr.open_group(run / "hotspots.zarr", mode="r+")
    short = np.zeros((1200, 1200), dtype=np.uint8)
    short[350, 350] = 1
    open_related = short.copy()
    open_related[450, 450] = 1
    open_unrelated = np.zeros_like(short)
    open_unrelated[550, 550] = 1
    cno = np.zeros_like(short)
    cno[350, 350] = 1
    for si in (0, 1):
        for ui, item in enumerate(
            (short, open_related, open_unrelated, cno if si == 0 else cno * 0)
        ):
            group["nominal_failure"][si, 0, ui] = item
        group["union_failure"][si, 0, 3] = short | open_related | open_unrelated
    group = zarr.open_group(run / "hotspots.zarr", mode="r+")
    (run / "run_manifest.json").write_text(
        json.dumps(
            {
                "product_id": manifest.product_id,
                "config_hash": canonical_hash(manifest.model_dump(mode="json")),
                "random_seed": None,
            }
        )
    )
    pd.DataFrame([{"ca_allocated_nominal_um2": 9, "ca_unique_nominal_um2": 7}]).to_parquet(
        run / "ca.parquet"
    )
    pd.DataFrame([{"ca_allocated_nominal_um2": 8, "ca_unique_nominal_um2": 6}]).to_parquet(
        run / "ca_union.parquet"
    )
    return manifest, sidecar, product_path, sidecar_path, run, group, open_related, cno


def test_relation_requires_exact_link_path_connection_and_shared_step() -> None:
    manifest, masks, specs, sidecar = _contract()
    assert _relations(manifest, masks, sidecar, specs) == {"OPEN_1": ("cno:CNT:A__B",)}
    changes = [
        ("cno_link_id", "cno:unknown:A__B", "unknown"),
        ("open_path_id", "OPEN_MISSING", "unknown"),
        ("circuit_connection_id", "CONN_OTHER", "identity"),
    ]
    for field, value, message in changes:
        bad = deepcopy(sidecar)
        bad["circuit"]["cno_open_relations"][0][field] = value
        with pytest.raises(InputContractError, match=message):
            _relations(manifest, masks, bad, specs)
    bad = deepcopy(sidecar)
    bad["circuit"]["cno_open_relations"].append(deepcopy(bad["circuit"]["cno_open_relations"][0]))
    with pytest.raises(InputContractError, match="duplicate"):
        _relations(manifest, masks, bad, specs)
    with pytest.raises(InputContractError, match="common source step"):
        _relations(
            manifest,
            masks,
            sidecar,
            [
                specs[0].__class__(
                    "OPEN_1", "lateral_open", KIND_OPEN_PATH, "A", "N/A", ("OTHER_STEP",)
                ),
                specs[1],
            ],
        )
    bad = deepcopy(sidecar)
    del bad["circuit"]["open_path_circuit_metadata"][0]["circuit_connection_id"]
    with pytest.raises(InputContractError, match="identity"):
        _relations(manifest, masks, bad, specs)
    bad = deepcopy(sidecar)
    bad["circuit"]["open_path_circuit_metadata"].append(
        deepcopy(bad["circuit"]["open_path_circuit_metadata"][0])
    )
    with pytest.raises(InputContractError, match="duplicate"):
        _relations(manifest, masks, bad, specs)
    with pytest.raises(InputContractError, match="auto_open"):
        _relations(
            manifest,
            masks,
            sidecar,
            [
                *specs,
                UnitSpec(
                    "auto_open:A", "lateral_open", KIND_OPEN_ELECTRODE, "A", "N/A", ("CNT_STEP",)
                ),
            ],
        )


def test_derived_maps_keep_raw_tables_shorts_and_unrelated_open(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    data = build_case("parallel_lines").manifest
    data["schema_version"] = "2.0"
    data["coordinate_system"].update(
        pixel_size_x_um=0.1,
        pixel_size_y_um=0.1,
        roi_width_um=120.0,
        roi_height_um=120.0,
    )
    data["particle_source_steps"] = data["particle_source_steps"][:2]
    data["particle_sizes_um"] = [1.0]
    data["electrodes"][0]["layer_name"] = "M1"
    data["electrode_pairs"][0]["failure_mode"] = "vertical_short"
    data["electrodes"].append(
        {
            "electrode_id": "C",
            "layer_name": "M2",
            "net_id": "N1",
            "electrode_role": "PAD",
            "path": "masks/electrode_C.png",
        }
    )
    data["contacts"] = [
        {
            "contact_id": "CNT",
            "path": "masks/contact_CNT.png",
            "hole_layer": "GI",
            "applicable_particle_source_steps": ["S_M1"],
        }
    ]
    data["open_applicable_steps_by_layer"] = {"M1": [], "M2": []}
    data["open_paths"] = [
        {
            "path_id": path_id,
            "electrode": electrode,
            "terminal_a_yx_px": [400, 400],
            "terminal_b_yx_px": [800, 800],
            "applicable_particle_source_steps": ["S_M1", "S_M2"],
            "fallback": "error",
        }
        for path_id, electrode in (("OPEN_1", "A"), ("OPEN_2", "B"))
    ]
    manifest = ProductManifest.model_validate(data)
    h1, h2 = _hole(1, 30), _hole(2, 120)
    link = ContactLink("cno:CNT:A__C", "CNT", "A", "C", (h1, h2))
    masks = MaskSet(
        GridSpec((1200, 1200), 0.1, 0.1),
        {},
        {},
        {},
        contact_topology=ContactTopology((h1, h2), {link.unit_id: link}),
    )
    monkeypatch.setattr("particle_sim.validation.cno_open_dedup.load_manifest", lambda _: manifest)
    monkeypatch.setattr(
        "particle_sim.validation.cno_open_dedup.load_product_masks", lambda *_: masks
    )
    product_path = tmp_path / "product.yaml"
    product_path.write_text(json.dumps(manifest.model_dump(mode="json"), sort_keys=True))
    sidecar = {
        "product_yaml_sha256": hashlib.sha256(product_path.read_bytes()).hexdigest(),
        "circuit": {
            "net_status": "declared",
            "auto_open": "disabled",
            "open_path_circuit_metadata": [
                {"path_id": "OPEN_1", "circuit_connection_id": "CONN_1", "circuit_net_id": "N1"},
                {"path_id": "OPEN_2", "circuit_connection_id": "CONN_2", "circuit_net_id": "N2"},
            ],
            "cno_open_relations": [
                {
                    "cno_link_id": link.unit_id,
                    "open_path_id": "OPEN_1",
                    "circuit_connection_id": "CONN_1",
                }
            ],
        },
    }
    sidecar_path = tmp_path / "adapter-provenance.json"
    sidecar_path.write_text(json.dumps(sidecar))
    run = tmp_path / "raw"
    run.mkdir()
    ids = ["A__B", "OPEN_1", "OPEN_2", link.unit_id]
    HotspotStore(
        run / "hotspots.zarr",
        source_steps=["S_M1", "S_M2"],
        particle_sizes_um=[1.0],
        pair_ids=ids,
        pair_modes=["vertical_short", "lateral_open", "lateral_open", "contact_not_open"],
        unit_kinds=["electrode_pair", "open_path", "open_path", "contact"],
        union_keys=["vertical_short", "lateral_open", "contact_not_open", "any"],
        shape_yx=(1200, 1200),
    )
    group = zarr.open_group(run / "hotspots.zarr", mode="r+")
    short = np.zeros((1200, 1200), dtype=np.uint8)
    short[350, 350] = 1
    open_related = short.copy()
    open_related[450, 450] = 1
    open_unrelated = np.zeros_like(short)
    open_unrelated[550, 550] = 1
    cno = np.zeros_like(short)
    cno[350, 350] = 1
    for si in (0, 1):
        for ui, item in enumerate(
            (short, open_related, open_unrelated, cno if si == 0 else cno * 0)
        ):
            group["nominal_failure"][si, 0, ui] = item
        group["union_failure"][si, 0, 3] = short | open_related | open_unrelated
    group = open_hotspot_store(run / "hotspots.zarr")
    (run / "run_manifest.json").write_text(
        json.dumps(
            {
                "product_id": manifest.product_id,
                "config_hash": canonical_hash(manifest.model_dump(mode="json")),
                "random_seed": None,
            }
        )
    )
    raw_table = pd.DataFrame([{"ca_allocated_nominal_um2": 9, "ca_unique_nominal_um2": 7}])
    raw_table.to_parquet(run / "ca.parquet")
    raw_table.to_parquet(run / "ca_union.parquet")
    raw_bytes = {name: (run / name).read_bytes() for name in ("ca.parquet", "ca_union.parquet")}
    unrelated_bytes = np.asarray(group["nominal_failure"][0, 0, 2]).tobytes()
    short_bytes = np.asarray(group["nominal_failure"][0, 0, 0]).tobytes()
    any_bytes = np.asarray(group["union_failure"][0, 0, 3]).tobytes()

    result = analyze_cno_open(product_path, sidecar_path, run, tmp_path / "derived")
    assert result["central_denominator_pixels"] == 360_000
    table = pd.read_parquet(tmp_path / "derived" / "cno_open_dedup.parquet")
    rows = {(row.source_step, row.open_path_id): row for row in table.itertuples()}
    assert rows["S_M1", "OPEN_1"].removed_pixels == 1
    assert rows["S_M1", "OPEN_1"].final_pixels == 1
    assert rows["S_M1", "OPEN_1"].central_removed_pixels == 1
    assert rows["S_M2", "OPEN_1"].removed_pixels == 0
    assert rows["S_M1", "OPEN_2"].removed_pixels == 0
    derived = open_hotspot_store(tmp_path / "derived" / "dedup_open.zarr")
    final = np.asarray(derived["open_final"][0, 0, 0], dtype=bool)
    removed = np.asarray(derived["removed"][0, 0, 0], dtype=bool)
    np.testing.assert_array_equal(final | removed, open_related)
    np.testing.assert_array_equal(
        final | removed | cno | short | open_unrelated,
        np.frombuffer(any_bytes, dtype=np.uint8).reshape(1200, 1200),
    )
    assert np.asarray(group["nominal_failure"][0, 0, 0]).tobytes() == short_bytes
    assert np.asarray(group["nominal_failure"][0, 0, 2]).tobytes() == unrelated_bytes
    assert np.asarray(group["union_failure"][0, 0, 3]).tobytes() == any_bytes
    assert all((run / name).read_bytes() == value for name, value in raw_bytes.items())
    assert hashlib.sha256(raw_bytes["ca.parquet"]).hexdigest() == result["raw_ca_sha256"]


@pytest.mark.parametrize(
    "probe",
    [
        "non_binary_open_map",
        "any_union_omits_open",
        "any_union_omits_cno",
        "output_inside_run",
        "nonempty_output",
        "random_seed_set",
        "config_hash_mismatch",
        "missing_ca_union",
        "sidecar_not_mapping",
        "invalid_sidecar_json",
        "net_status_undeclared",
        "link_net_mismatch",
        "path_net_mismatch",
        "fallback_not_error",
        "tampered_sidecar_hash",
        "missing_sidecar_hash",
        "product_hash_mismatch",
    ],
)
def test_malformed_e3_inputs_never_leave_derived_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, probe: str
) -> None:
    manifest, sidecar, product_path, sidecar_path, run, group, _, _ = _prepare_dedup_case(
        tmp_path, monkeypatch
    )
    output = tmp_path / "derived"
    if probe == "non_binary_open_map":
        group["nominal_failure"][0, 0, 1, 10, 10] = 2
    elif probe == "any_union_omits_open":
        group["union_failure"][0, 0, 3, 450, 450] = 0
    elif probe == "any_union_omits_cno":
        group["nominal_failure"][0, 0, 3, 370, 370] = 1
    elif probe == "output_inside_run":
        output = run / "derived"
    elif probe == "nonempty_output":
        output.mkdir()
        (output / "keep.txt").write_text("existing")
    elif probe == "random_seed_set":
        run_manifest = json.loads((run / "run_manifest.json").read_text())
        run_manifest["random_seed"] = 17
        (run / "run_manifest.json").write_text(json.dumps(run_manifest))
    elif probe == "config_hash_mismatch":
        run_manifest = json.loads((run / "run_manifest.json").read_text())
        run_manifest["config_hash"] = "mismatch"
        (run / "run_manifest.json").write_text(json.dumps(run_manifest))
    elif probe == "missing_ca_union":
        (run / "ca_union.parquet").unlink()
    elif probe == "sidecar_not_mapping":
        sidecar_path.write_text("[]")
    elif probe == "invalid_sidecar_json":
        sidecar_path.write_text("{")
    elif probe == "net_status_undeclared":
        sidecar["circuit"]["net_status"] = "inferred"
        sidecar_path.write_text(json.dumps(sidecar))
    elif probe == "link_net_mismatch":
        electrodes = list(manifest.electrodes)
        electrodes[2] = electrodes[2].model_copy(update={"net_id": "OTHER_NET"})
        manifest = manifest.model_copy(update={"electrodes": electrodes})
        monkeypatch.setattr(
            "particle_sim.validation.cno_open_dedup.load_manifest", lambda _: manifest
        )
    elif probe == "path_net_mismatch":
        sidecar["circuit"]["open_path_circuit_metadata"][0]["circuit_net_id"] = "OTHER_NET"
        sidecar_path.write_text(json.dumps(sidecar))
    elif probe == "fallback_not_error":
        paths = list(manifest.open_paths)
        paths[0] = paths[0].model_copy(update={"fallback": "roi_boundary"})
        manifest = manifest.model_copy(update={"open_paths": paths})
        monkeypatch.setattr(
            "particle_sim.validation.cno_open_dedup.load_manifest", lambda _: manifest
        )
    elif probe == "tampered_sidecar_hash":
        sidecar["product_yaml_sha256"] = "0" * 64
        sidecar_path.write_text(json.dumps(sidecar))
    elif probe == "missing_sidecar_hash":
        sidecar.pop("product_yaml_sha256")
        sidecar_path.write_text(json.dumps(sidecar))
    elif probe == "product_hash_mismatch":
        product_path.write_text(product_path.read_text() + "\n# tampered\n")

    expected_error = (
        "product.yaml SHA-256"
        if probe in {"tampered_sidecar_hash", "missing_sidecar_hash", "product_hash_mismatch"}
        else None
    )
    with pytest.raises(InputContractError, match=expected_error):
        analyze_cno_open(product_path, sidecar_path, run, output)

    if probe == "nonempty_output":
        assert sorted(path.name for path in output.iterdir()) == ["keep.txt"]
    else:
        assert not output.exists()
    assert not list(output.parent.glob(f".{output.name}.*"))
