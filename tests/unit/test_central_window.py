from __future__ import annotations

import dataclasses
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from particle_sim.config import canonical_hash, file_sha256
from particle_sim.errors import InputContractError
from particle_sim.geometry.contacts import ContactHole, ContactLink, ContactTopology
from particle_sim.io.masks import GridSpec, MaskSet
from particle_sim.io.outputs import HotspotStore, ProbabilityStore, open_hotspot_store
from particle_sim.schema import PARTICLE_SIZES_UM, Contact, OpenPath, ProductManifest
from particle_sim.validation.central_window import (
    NONLOCAL_OPEN,
    analyze_central,
    central_area_um2,
    support_rows,
)
from particle_sim.validation.synthetic_cases import build_case


def _inputs(margin_um: float = 0.0) -> tuple[ProductManifest, MaskSet]:
    data = build_case("parallel_lines").manifest
    data["coordinate_system"].update(
        pixel_size_x_um=0.1, pixel_size_y_um=0.1,
        roi_width_um=120.0, roi_height_um=120.0,
    )
    data["particle_source_steps"] = data["particle_source_steps"][1:2]
    data["particle_sizes_um"] = list(PARTICLE_SIZES_UM)
    data["contact_margin_um"] = margin_um
    manifest = ProductManifest.model_validate(data)
    masks = MaskSet(
        GridSpec((1200, 1200), 0.1, 0.1),
        layers={}, electrodes={}, input_hashes={},
    )
    return manifest, masks


def test_exact_bins_window_and_margin_boundary() -> None:
    manifest, masks = _inputs(6.0)
    rows = support_rows(manifest, masks)
    assert len(rows) == 32
    assert [row["particle_size_um"] for row in rows] == list(PARTICLE_SIZES_UM)
    by_size = {row["particle_size_um"]: row for row in rows}
    assert [by_size[d]["radius_px"] for d in (1.0, 30.0, 40.0, 50.0)] == [
        5.0, 150.0, 200.0, 250.0,
    ]
    assert all(by_size[d]["supported"] for d in (1.0, 30.0, 40.0))
    assert not by_size[50.0]["supported"]
    assert "short_geometry_reach_exceeds_300px" in by_size[50.0]["reasons"]
    exact_manifest, _ = _inputs(5.0)
    assert support_rows(exact_manifest, masks)[31]["reaches_px"]["short_geometry"] == 300.0
    assert support_rows(exact_manifest, masks)[31]["supported"]

    full = np.zeros((1200, 1200), dtype=bool)
    full[300:900, 300:900] = True
    full[0, 0] = True
    assert central_area_um2(full, 0.01) == pytest.approx(3600.0)
    assert int(full[300:900, 300:900].sum()) == 360_000
    bad = manifest.model_copy(update={"particle_sizes_um": [1.0, 30.0, 40.0, 50.0]})
    with pytest.raises(InputContractError, match="32 bins"):
        support_rows(bad, masks)
    maskless_stage = build_case("parallel_lines").manifest["particle_source_steps"][2]
    with pytest.raises(InputContractError, match="maskless dielectric"):
        support_rows(manifest.model_copy(update={"particle_source_steps": [
            manifest.particle_source_steps[0].model_copy(update={
                "step_id": maskless_stage["step_id"],
                "lands_after_layer": maskless_stage["lands_after_layer"],
            }),
        ]}), masks)


def test_full_hole_extent_and_nonlocal_open_are_reported() -> None:
    manifest, masks = _inputs()
    small_y, small_x = np.mgrid[500:515, 500:515]
    small_hole = ContactHole(
        "CNT", 1, "A", "B", small_y.ravel(), small_x.ravel(),
        small_y.ravel(), small_x.ravel(), None, small_y.size, small_y.size,
    )
    ys, xs = np.mgrid[500:725, 500:725]
    hole = ContactHole(
        "CNT", 1, "A", "B", ys.ravel(), xs.ravel(), ys.ravel(), xs.ravel(),
        None, ys.size, ys.size,
    )
    link = ContactLink("cno:CNT:A__B", "CNT", "A", "B", (hole,))
    topology = ContactTopology((hole,), {link.unit_id: link})
    masks = MaskSet(masks.grid, {}, {}, {}, contact_topology=topology)
    contact = Contact(
        contact_id="CNT", path="unused", hole_layer="M2",
        applicable_particle_source_steps=[manifest.particle_source_steps[0].step_id],
    )
    open_path = OpenPath(
        path_id="OP", electrode="A", terminal_a_yx_px=(1, 1),
        terminal_b_yx_px=(2, 2),
    )
    manifest = manifest.model_copy(update={"contacts": [contact], "open_paths": [open_path]})
    small_link = ContactLink(link.unit_id, "CNT", "A", "B", (small_hole,))
    small_masks = MaskSet(
        masks.grid, {}, {}, {},
        contact_topology=ContactTopology((small_hole,), {link.unit_id: small_link}),
    )
    small_rows = support_rows(manifest, small_masks)
    assert small_hole.cut_pixels == 225
    assert small_rows[31]["full_cut_hole_half_span_px"] == 7.0
    assert small_rows[31]["reaches_px"]["full_cut_hole"] == 257.0
    rows = support_rows(manifest, masks)
    assert len(rows) == 32
    d30 = rows[29]
    d40 = rows[30]
    assert d30["full_cut_hole_half_span_px"] == 112.0
    assert d30["reaches_px"]["full_cut_hole"] == 262.0
    assert d30["local_geometry_supported"]
    assert not d30["supported"] and NONLOCAL_OPEN in d30["reasons"]
    assert "full_cut_hole_reach_exceeds_300px" in d40["reasons"]
    assert "flat_film_approximation" in d40["source_stage_semantics"]
    assert support_rows(manifest.model_copy(update={"particle_source_steps": []}), masks) == []


def test_analysis_keeps_full_map_and_writes_central_sidecars(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    manifest, masks = _inputs()
    run = tmp_path / "run"
    run.mkdir()
    store = HotspotStore(
        run / "hotspots.zarr", source_steps=[manifest.particle_source_steps[0].step_id],
        particle_sizes_um=list(PARTICLE_SIZES_UM), pair_ids=["A__B"],
        pair_modes=["lateral_short"], unit_kinds=["electrode_pair"],
        union_keys=["lateral_short", "any"], shape_yx=(1200, 1200),
    )
    full = np.zeros((1200, 1200), dtype=bool)
    full[300, 300] = True
    full[0, 0] = True
    store.write_pair(0, 0, 0, legacy=full, nominal=full,
                     margin_um=np.zeros((1200, 1200), dtype=np.float32))
    store.write_union(0, 0, 0, any_failure=full, pair_count=full.astype(np.uint16))
    (run / "run_manifest.json").write_text(json.dumps({
        "product_id": manifest.product_id,
        "config_hash": canonical_hash(manifest.model_dump(mode="json")),
    }), encoding="utf-8")
    monkeypatch.setattr("particle_sim.validation.central_window.load_manifest", lambda _: manifest)
    monkeypatch.setattr(
        "particle_sim.validation.central_window.load_product_masks", lambda *_: masks
    )
    result = analyze_central(tmp_path / "unused.yaml", tmp_path / "central", run_dir=run)
    assert result == {"support_rows": 32, "unit_rows": 32, "union_rows": 64,
                      "unsupported": 0}
    audit = json.loads((tmp_path / "central" / "central_support.json").read_text())
    assert audit["central_pixels"] == 360_000
    assert audit["halo_px"] == 300
    assert audit["central_area_um2"] == pytest.approx(3600.0)
    table = pd.read_parquet(tmp_path / "central" / "central_ca.parquet")
    assert table.iloc[0].central_ca_um2 == pytest.approx(0.01)
    assert table.iloc[0].central_normalized == pytest.approx(1 / 360_000)
    assert bool(table.iloc[0].supported)
    raw = open_hotspot_store(run / "hotspots.zarr")
    assert int(np.asarray(raw["nominal_failure"][0, 0, 0]).sum()) == 2


def test_probability_run_is_matched_by_manifest_hash(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    manifest, masks = _inputs()
    masks = dataclasses.replace(masks, input_hashes={"masks/M2.png": "a" * 64})
    manifest_path = tmp_path / "product.yaml"
    manifest_path.write_text("product", encoding="utf-8")
    run = tmp_path / "run"
    run.mkdir()
    store = ProbabilityStore(
        run / "probability.zarr", source_steps=[manifest.particle_source_steps[0].step_id],
        particle_sizes_um=list(PARTICLE_SIZES_UM), unit_ids=["A__B"],
        unit_modes=["lateral_short"], unit_kinds=["electrode_pair"],
        union_keys=["lateral_short", "any"], shape_yx=(1200, 1200),
    )
    zeros = np.zeros((1200, 1200), dtype=np.float32)
    for size_index in range(len(PARTICLE_SIZES_UM)):
        store.write_unit(0, size_index, 0, probability=zeros, margin_mean_um=zeros)
        for union_index in range(2):
            store.write_union(0, size_index, union_index, probability=zeros)
    run_manifest = {
        "product_id": manifest.product_id,
        "design_version": manifest.design_version,
        # Covers manifest and scenario, so it never equals the manifest-only hash.
        "config_hash": canonical_hash({"manifest": manifest.model_dump(mode="json"),
                                       "scenario": {}}),
        "input_hashes": {"manifest": file_sha256(manifest_path), "masks/M2.png": "a" * 64},
    }
    (run / "run_manifest.json").write_text(json.dumps(run_manifest), encoding="utf-8")
    monkeypatch.setattr("particle_sim.validation.central_window.load_manifest", lambda _: manifest)
    monkeypatch.setattr(
        "particle_sim.validation.central_window.load_product_masks", lambda *_: masks
    )
    result = analyze_central(manifest_path, tmp_path / "central", run_dir=run)
    assert result["unit_rows"] == 32 and result["unsupported"] == 32
    audit = json.loads((tmp_path / "central" / "central_support.json").read_text())
    assert all(
        "scenario_variation_influence_not_bounded" in row["reasons"] for row in audit["support"]
    )

    other_mask = dataclasses.replace(masks, input_hashes={"masks/M2.png": "b" * 64})
    monkeypatch.setattr(
        "particle_sim.validation.central_window.load_product_masks", lambda *_: other_mask
    )
    with pytest.raises(InputContractError, match="does not match product config"):
        analyze_central(manifest_path, tmp_path / "central-other-mask", run_dir=run)

    manifest_path.write_text("other product", encoding="utf-8")
    with pytest.raises(InputContractError, match="does not match product config"):
        analyze_central(manifest_path, tmp_path / "central-other", run_dir=run)
