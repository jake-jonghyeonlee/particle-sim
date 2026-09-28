from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from particle_sim.aggregation.critical_area import critical_area_um2
from particle_sim.config import load_manifest
from particle_sim.geometry.legacy import legacy_short_hotspot
from particle_sim.io.masks import load_product_masks
from particle_sim.io.outputs import open_hotspot_store
from particle_sim.pipeline import _unpack_legacy_map, run_nominal
from particle_sim.validation.synthetic_cases import build_case, write_case


@pytest.mark.parametrize("shape", [(1, 1), (3, 5), (7, 9)])
@pytest.mark.parametrize("fill", [False, True, None])
def test_legacy_packbits_round_trip(shape: tuple[int, int], fill: bool | None) -> None:
    if fill is None:
        source = np.arange(shape[0] * shape[1]).reshape(shape) % 3 == 0
    else:
        source = np.full(shape, fill, dtype=bool)

    packed = np.packbits(source)
    restored = _unpack_legacy_map(packed, shape)

    assert packed.nbytes == (source.size + 7) // 8
    assert restored.dtype == np.dtype(bool)
    np.testing.assert_array_equal(restored, source)


def test_nominal_legacy_maps_and_unions_match_direct_calculation(tmp_path: Path) -> None:
    case = build_case("pilot")
    case.manifest["particle_sizes_um"] = [2.0, 5.0]
    manifest_path = write_case(case, tmp_path / "input")
    manifest = load_manifest(manifest_path)
    masks = load_product_masks(manifest, manifest_path)
    result = run_nominal(
        manifest_path, tmp_path / "output", now=datetime(2026, 9, 29, tzinfo=UTC)
    )
    store = open_hotspot_store(result.hotspots)
    unions = pd.read_parquet(result.union_table)

    for size_index, diameter in enumerate(manifest.particle_sizes_um):
        maps = {}
        for pair_index, pair in enumerate(manifest.electrode_pairs):
            legacy = legacy_short_hotspot(
                masks.electrodes[pair.electrode_a],
                masks.electrodes[pair.electrode_b],
                diameter,
                masks.grid,
                manifest.legacy_kernel,
            )
            maps[pair.pair_id] = legacy
            for step_index in range(len(manifest.particle_source_steps)):
                np.testing.assert_array_equal(
                    store["legacy_failure"][step_index, size_index, pair_index], legacy
                )

        for _, row in unions[unions["particle_size_um"] == diameter].iterrows():
            expected = np.zeros(masks.grid.shape_yx, dtype=bool)
            for pair_id in row["pair_ids"].split(";"):
                if pair_id in maps:
                    expected |= maps[pair_id]
            assert row["ca_union_legacy_um2"] == critical_area_um2(expected, masks.grid)


def test_nominal_missing_legacy_key_uses_empty_map(tmp_path: Path) -> None:
    manifest_path = write_case(build_case("open_straight"), tmp_path / "input")
    result = run_nominal(
        manifest_path, tmp_path / "output", now=datetime(2026, 9, 29, tzinfo=UTC)
    )
    store = open_hotspot_store(result.hotspots)
    assert not np.any(store["legacy_failure"][:])
    assert pd.read_parquet(result.union_table)["ca_union_legacy_um2"].eq(0).all()
