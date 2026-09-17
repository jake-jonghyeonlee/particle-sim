from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import tifffile

from particle_sim.aggregation.critical_area import critical_area_um2, normalized_ca, um2_to_mm2
from particle_sim.io.masks import GridSpec
from particle_sim.io.outputs import (
    CA_COLUMNS,
    HotspotStore,
    build_run_manifest,
    make_run_id,
    open_hotspot_store,
    write_float_tiff,
    write_json,
    write_table,
)


def _grid(shape=(10, 10), py=0.5, px=0.5) -> GridSpec:
    return GridSpec(shape_yx=shape, pixel_size_y_um=py, pixel_size_x_um=px)


def test_critical_area_values() -> None:
    grid = _grid()
    assert critical_area_um2(np.ones((10, 10), dtype=bool), grid) == pytest.approx(25.0)
    assert critical_area_um2(np.full((10, 10), 0.5), grid) == pytest.approx(12.5)
    assert um2_to_mm2(25.0) == pytest.approx(2.5e-5)
    assert normalized_ca(25.0, grid) == pytest.approx(1.0)
    assert normalized_ca(12.5, grid) == pytest.approx(0.5)


def test_critical_area_rejects_bad_maps() -> None:
    grid = _grid()
    with pytest.raises(ValueError):
        critical_area_um2(np.ones((4, 4)), grid)
    with pytest.raises(ValueError):
        critical_area_um2(np.full((10, 10), np.nan), grid)
    with pytest.raises(ValueError):
        critical_area_um2(np.full((10, 10), 1.5), grid)


def _ca_row() -> dict:
    return {name: 0 for name in CA_COLUMNS}


def test_write_table_enforces_columns(tmp_path: Path) -> None:
    path = tmp_path / "ca.parquet"
    write_table([_ca_row()], CA_COLUMNS, path)
    assert tuple(pd.read_parquet(path).columns) == CA_COLUMNS
    missing = _ca_row()
    del missing["run_id"]
    with pytest.raises(ValueError, match="누락"):
        write_table([missing], CA_COLUMNS, path)
    extra = _ca_row() | {"surprise": 1}
    with pytest.raises(ValueError, match="추가"):
        write_table([extra], CA_COLUMNS, path)


def test_hotspot_store_round_trip(tmp_path: Path) -> None:
    store = HotspotStore(
        tmp_path / "h.zarr",
        source_steps=["S1"],
        particle_sizes_um=[2.0, 4.0],
        pair_ids=["A__B"],
        pair_modes=["lateral_short"],
        union_keys=["lateral_short", "any"],
        shape_yx=(6, 6),
    )
    failure = np.zeros((6, 6), dtype=bool)
    failure[1, 2] = True
    margin = np.full((6, 6), -np.inf, dtype=np.float32)
    margin[1, 2] = 0.25
    store.write_pair(0, 1, 0, legacy=failure, nominal=failure, margin_um=margin)
    store.write_union(0, 1, 1, any_failure=failure, pair_count=failure.astype(np.uint16))
    with pytest.raises(ValueError):
        store.write_pair(0, 0, 0, legacy=np.zeros((2, 2), dtype=bool), nominal=failure,
                         margin_um=margin)

    opened = open_hotspot_store(tmp_path / "h.zarr")
    assert opened["nominal_failure"].shape == (1, 2, 1, 6, 6)
    assert opened["union_failure"].shape == (1, 2, 2, 6, 6)
    assert opened["pair_margin_um"].dtype == np.float32
    assert bool(np.asarray(opened["nominal_failure"][0, 1, 0])[1, 2])
    assert np.asarray(opened["pair_margin_um"][0, 1, 0])[1, 2] == pytest.approx(0.25)
    assert list(opened.attrs["union_keys"]) == ["lateral_short", "any"]
    assert "not failure probability" in opened.attrs["value_semantics"]
    with pytest.raises((FileExistsError, ValueError, OSError)):
        HotspotStore(
            tmp_path / "h.zarr", source_steps=["S1"], particle_sizes_um=[2.0],
            pair_ids=["A__B"], pair_modes=["lateral_short"], union_keys=["any"], shape_yx=(6, 6),
        )


def test_run_manifest_and_helpers(tmp_path: Path) -> None:
    created = datetime(2026, 9, 18, 1, 2, 3, tzinfo=UTC)
    assert make_run_id("P1", created) == "20260918T010203Z_P1"
    manifest = build_run_manifest(
        run_id="R", product_id="P1", design_version="V1", config_hash="abc",
        input_hashes={"manifest": "d"}, particle_sizes_um=[1.0], source_steps=["S1"],
        assumption_set="protrusion", warnings=["TODO_PHYSICS: x"], legacy_kernel={"size": 5},
        post_step_policy={"mode": "protrusion"}, created_at=created,
    )
    assert manifest["units"] == {"length": "um", "area": "um2"}
    assert manifest["random_seed"] is None
    assert manifest["simulator_version"]
    path = tmp_path / "run_manifest.json"
    write_json(manifest, path)
    assert "TODO_PHYSICS" in path.read_text(encoding="utf-8")

    tiff_path = tmp_path / "m.tif"
    write_float_tiff(np.array([[0.5, 1.5]]), tiff_path)
    assert tifffile.imread(str(tiff_path)).dtype == np.float32


def test_hotspot_store_rejects_out_of_range_index(tmp_path: Path) -> None:
    store = HotspotStore(
        tmp_path / "bounds.zarr",
        source_steps=["S1"],
        particle_sizes_um=[2.0],
        pair_ids=["A__B"],
        pair_modes=["lateral_short"],
        union_keys=["lateral_short", "any"],
        shape_yx=(4, 4),
    )
    failure = np.zeros((4, 4), dtype=bool)
    margin = np.zeros((4, 4), dtype=np.float32)
    with pytest.raises(ValueError, match="step_index"):
        store.write_pair(5, 0, 0, legacy=failure, nominal=failure, margin_um=margin)
    with pytest.raises(ValueError, match="size_index"):
        store.write_pair(0, 3, 0, legacy=failure, nominal=failure, margin_um=margin)
    with pytest.raises(ValueError, match="pair_or_union_index"):
        store.write_pair(0, 0, 2, legacy=failure, nominal=failure, margin_um=margin)
    with pytest.raises(ValueError, match="pair_or_union_index"):
        store.write_union(0, 0, 9, any_failure=failure, pair_count=failure.astype(np.uint16))
