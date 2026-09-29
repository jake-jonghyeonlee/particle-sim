from __future__ import annotations

from copy import deepcopy
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import particle_sim.pipeline as pipeline
import particle_sim.units as units
from particle_sim.io.outputs import open_hotspot_store
from particle_sim.validation.synthetic_cases import build_case, write_case


def test_nominal_open_cache_expires_after_last_applicable_step(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    case = build_case("open_straight")
    template = case.manifest["open_paths"][0]
    case.manifest["open_paths"] = []
    steps = [step["step_id"] for step in case.manifest["particle_source_steps"]]
    for path_id, applicable in (
        ("ONE", [steps[0]]),
        ("MULTI", [steps[2], steps[0]]),
        ("ALL", []),
    ):
        path = deepcopy(template)
        path["path_id"] = path_id
        path["applicable_particle_source_steps"] = applicable
        case.manifest["open_paths"].append(path)
    case.manifest["particle_sizes_um"] = [6.0, 8.0]
    case.manifest["contact_margin_um"] = 0.5
    manifest_path = write_case(case, tmp_path / "input")
    now = datetime(2026, 9, 29, tzinfo=UTC)

    original_evaluate = pipeline.evaluate_units
    unbounded_cache: units.OpenResultCache = {}

    def unbounded(*args, **kwargs):
        return original_evaluate(*args[:6], unbounded_cache, **kwargs)

    monkeypatch.setattr(pipeline, "evaluate_units", unbounded)
    baseline = pipeline.run_nominal(manifest_path, tmp_path / "baseline", now=now)

    observations = []
    active_cache: units.OpenResultCache | None = None
    calls: dict[str, int] = {}
    original_open = units.simulate_open

    def counted_open(*args, **kwargs):
        path_id = args[2].path_id
        calls[path_id] = calls.get(path_id, 0) + 1
        return original_open(*args, **kwargs)

    def observed(*args, **kwargs):
        nonlocal active_cache
        active_cache = args[6]
        before = {key[0] for key in active_cache}
        results = original_evaluate(*args, **kwargs)
        after = {key[0] for key in active_cache}
        observations.append((args[3], args[4], before, after))
        return results

    monkeypatch.setattr(units, "simulate_open", counted_open)
    monkeypatch.setattr(pipeline, "evaluate_units", observed)
    candidate = pipeline.run_nominal(manifest_path, tmp_path / "candidate", now=now)

    assert observations == [
        (steps[0], 6.0, set(), {"ONE", "MULTI", "ALL"}),
        (steps[0], 8.0, {"MULTI", "ALL"}, {"ONE", "MULTI", "ALL"}),
        (steps[1], 6.0, {"MULTI", "ALL"}, {"MULTI", "ALL"}),
        (steps[1], 8.0, {"MULTI", "ALL"}, {"MULTI", "ALL"}),
        (steps[2], 6.0, {"MULTI", "ALL"}, {"MULTI", "ALL"}),
        (steps[2], 8.0, {"MULTI", "ALL"}, {"MULTI", "ALL"}),
    ]
    assert active_cache == {}
    assert calls == {"ONE": 2, "MULTI": 2, "ALL": 2}
    assert len(unbounded_cache) == 6

    for name in ("ca.parquet", "ca_union.parquet"):
        pd.testing.assert_frame_equal(
            pd.read_parquet(baseline.output_dir / name),
            pd.read_parquet(candidate.output_dir / name),
            check_exact=True,
        )
    assert baseline.legacy_comparison.read_bytes() == candidate.legacy_comparison.read_bytes()
    base_store = open_hotspot_store(baseline.hotspots)
    candidate_store = open_hotspot_store(candidate.hotspots)
    assert dict(base_store.attrs) == dict(candidate_store.attrs)
    for name in base_store.array_keys():
        np.testing.assert_array_equal(base_store[name][:], candidate_store[name][:])


def test_nominal_zero_source_steps_with_open_path(tmp_path: Path) -> None:
    case = build_case("open_straight")
    case.manifest["particle_source_steps"] = []
    manifest_path = write_case(case, tmp_path / "input")

    result = pipeline.run_nominal(
        manifest_path,
        tmp_path / "output",
        now=datetime(2026, 9, 29, tzinfo=UTC),
    )

    assert pd.read_parquet(result.ca_table).empty
    assert pd.read_parquet(result.union_table).empty
    assert open_hotspot_store(result.hotspots).attrs["source_steps"] == []
