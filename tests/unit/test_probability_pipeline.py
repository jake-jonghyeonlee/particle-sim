from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import zarr

from particle_sim.aggregation.probability import ProbabilityAccumulator, probability_ca
from particle_sim.config import load_manifest, load_scenario_config
from particle_sim.errors import InputContractError
from particle_sim.geometry.transforms import (
    apply_scenario_to_masks,
    apply_scenario_to_stack,
    transform_mask,
)
from particle_sim.io.masks import GridSpec, load_product_masks
from particle_sim.io.outputs import CA_COLUMNS
from particle_sim.io.stack import build_stack_z
from particle_sim.pipeline import run_nominal, run_probability
from particle_sim.uncertainty.sampler import Scenario
from particle_sim.validation.convergence import convergence_series, is_converged
from particle_sim.validation.synthetic_cases import build_case, write_case

NOMINAL_SCENARIO = """
schema_version: "1.0"
assumption_set: nominal
evidence: assumed
overlay_um:
  M2: {mean: 0.0, sigma: 0.3}
particle_diameter_delta_um: {mean: 0.0, sigma: 0.2}
"""
STRICT_SCENARIO = """
schema_version: "1.0"
assumption_set: strict
evidence: assumed
"""


def _scenario_file(tmp_path: Path, text: str, name: str) -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def _grid() -> GridSpec:
    return GridSpec(shape_yx=(20, 20), pixel_size_y_um=0.5, pixel_size_x_um=0.5)


def _empty_scenario(**kwargs) -> Scenario:
    base = {
        "scenario_id": 0,
        "thickness_delta_um": {},
        "overlay_xy_um": {},
        "cd_bias_um": {},
        "particle_diameter_delta_um": 0.0,
        "center_offset_xy_um": (0.0, 0.0),
        "embedding_depth_um": None,
        "planarization_factor": {},
        "contact_margin_um": None,
    }
    base.update(kwargs)
    return Scenario(**base)


def test_transform_mask_shifts_sub_pixel() -> None:
    grid = _grid()
    mask = np.zeros((20, 20), dtype=bool)
    mask[:, 8:12] = True
    same = transform_mask(mask, grid)
    assert np.array_equal(same, mask)
    shifted = transform_mask(mask, grid, overlay_xy_um=(1.0, 0.0))  # +2 px in x
    assert np.array_equal(shifted[:, 10:14], mask[:, 8:12])
    # the signed distance field reproduces the mask exactly at threshold 0
    assert np.array_equal(transform_mask(mask, grid, cd_bias_um=0.0), mask)
    grown = transform_mask(mask, grid, cd_bias_um=0.5)  # +1 px on each side
    assert np.where(grown[0])[0].tolist() == list(range(7, 13))
    shrunk = transform_mask(mask, grid, cd_bias_um=-0.5)  # -1 px on each side
    assert np.where(shrunk[0])[0].tolist() == [9, 10]


def test_apply_scenario_to_stack_reaccumulates_z(tmp_path: Path) -> None:
    manifest_path = write_case(build_case("pilot"), tmp_path)
    manifest = load_manifest(manifest_path)
    stack_z = build_stack_z(manifest)
    scenario = _empty_scenario(thickness_delta_um={"M1": 0.1})
    updated = apply_scenario_to_stack(stack_z, scenario)
    assert updated["M1"].thickness_um == pytest.approx(0.4)
    assert updated["GI"].z_bottom_um == pytest.approx(0.4)
    assert updated["M3"].z_top_um == pytest.approx(stack_z["M3"].z_top_um + 0.1)


def test_apply_scenario_to_masks_moves_only_named_layers(tmp_path: Path) -> None:
    manifest_path = write_case(build_case("pilot"), tmp_path)
    manifest = load_manifest(manifest_path)
    masks = load_product_masks(manifest, manifest_path)
    scenario = _empty_scenario(overlay_xy_um={"M2": (1.0, 0.0)})
    moved = apply_scenario_to_masks(manifest, masks, scenario)
    assert not np.array_equal(moved.electrodes["A"], masks.electrodes["A"])  # M2 electrode
    assert np.array_equal(moved.electrodes["G"], masks.electrodes["G"])  # M1 electrode
    assert moved.electrodes["A"].sum() == pytest.approx(masks.electrodes["A"].sum(), rel=0.05)


def test_probability_accumulator_counts_and_percentiles() -> None:
    grid = _grid()
    accumulator = ProbabilityAccumulator(grid.shape_yx, grid)
    always = np.zeros(grid.shape_yx, dtype=bool)
    always[0, 0] = True
    sometimes = always.copy()
    sometimes[0, 1] = True
    accumulator.add(always)
    accumulator.add(sometimes)
    probability = accumulator.probability()
    assert probability[0, 0] == pytest.approx(1.0)
    assert probability[0, 1] == pytest.approx(0.5)
    summary = probability_ca(accumulator)
    assert summary.scenario_count == 2
    assert summary.ca_probability_mean_um2 == pytest.approx(0.375)  # (1 + 0.5) * 0.25
    assert summary.ca_scenario_mean_um2 == pytest.approx(0.375)
    assert summary.ca_p05_um2 <= summary.ca_p50_um2 <= summary.ca_p95_um2
    with pytest.raises(ValueError):
        ProbabilityAccumulator(grid.shape_yx, grid).probability()


def test_convergence_series_and_tolerance() -> None:
    steps = convergence_series({64: 100.0, 128: 102.0, 256: 102.5})
    assert [s.scenario_count for s in steps] == [64, 128, 256]
    assert np.isnan(steps[0].delta_previous)
    assert steps[1].delta_previous == pytest.approx(2.0 / 102.0)
    assert is_converged(steps, 0.01)
    assert not is_converged(steps, 0.001)
    with pytest.raises(ValueError):
        convergence_series({})


def test_zero_variation_matches_nominal(tmp_path: Path) -> None:
    manifest_path = write_case(build_case("pilot"), tmp_path / "case")
    scenario_path = _scenario_file(tmp_path, STRICT_SCENARIO, "strict.yaml")
    nominal = run_nominal(manifest_path, tmp_path / "nominal")
    probability = run_probability(
        manifest_path, scenario_path, tmp_path / "prob", sobol_exponent_m=6, seed=1
    )
    assert probability.scenario_count == 1
    key = ["source_step", "particle_size_um", "pair_id"]
    left = pd.read_parquet(nominal.ca_table).set_index(key)["ca_nominal_um2"]
    right = pd.read_parquet(probability.ca_table).set_index(key)["ca_mean_um2"]
    assert float((left - right).abs().max()) == 0.0


def test_probability_run_outputs(tmp_path: Path) -> None:
    manifest_path = write_case(build_case("parallel_lines"), tmp_path / "case")
    scenario_path = _scenario_file(tmp_path, NOMINAL_SCENARIO, "nominal.yaml")
    result = run_probability(
        manifest_path, scenario_path, tmp_path / "run", sobol_exponent_m=6, seed=42
    )
    assert result.scenario_count == 64
    store = zarr.open_group(result.probability, mode="r")
    probability = np.asarray(store["probability"])
    assert probability.shape == (3, 5, 1, 200, 200)
    assert np.nanmin(probability) >= 0.0 and np.nanmax(probability) <= 1.0
    assert 0.0 < float(np.nanmax(probability)) <= 1.0
    assert "not an actual failure probability" in store.attrs["value_semantics"]

    ca = pd.read_parquet(result.ca_table)
    assert tuple(ca.columns) == CA_COLUMNS
    assert (ca["scenario_count"] == 64).all()
    assert (ca["ca_p05_um2"] <= ca["ca_p50_um2"]).all()
    assert (ca["ca_p50_um2"] <= ca["ca_p95_um2"]).all()

    with result.convergence.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows and all(int(r["scenario_count"]) == 64 for r in rows)

    manifest = json.loads(result.run_manifest.read_text(encoding="utf-8"))
    assert manifest["random_seed"] == 42
    assert manifest["sobol_exponent_m"] == 6
    assert any("evidence=assumed" in w for w in manifest["warnings"])
    assert "scenario" in manifest["input_hashes"]


def test_probability_run_is_reproducible(tmp_path: Path) -> None:
    manifest_path = write_case(build_case("parallel_lines"), tmp_path / "case")
    scenario_path = _scenario_file(tmp_path, NOMINAL_SCENARIO, "nominal.yaml")
    first = run_probability(
        manifest_path, scenario_path, tmp_path / "a", sobol_exponent_m=6, seed=7
    )
    second = run_probability(
        manifest_path, scenario_path, tmp_path / "b", sobol_exponent_m=6, seed=7
    )
    left = pd.read_parquet(first.ca_table)["ca_mean_um2"].to_numpy()
    right = pd.read_parquet(second.ca_table)["ca_mean_um2"].to_numpy()
    assert np.array_equal(left, right)


def test_convergence_prefixes_for_large_run(tmp_path: Path) -> None:
    manifest_path = write_case(build_case("parallel_lines"), tmp_path / "case")
    scenario_path = _scenario_file(tmp_path, NOMINAL_SCENARIO, "nominal.yaml")
    result = run_probability(
        manifest_path, scenario_path, tmp_path / "run", sobol_exponent_m=7, seed=3
    )
    with result.convergence.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    counts = sorted({int(r["scenario_count"]) for r in rows})
    assert counts == [64, 128]
    deltas = [float(r["delta_previous"]) for r in rows if int(r["scenario_count"]) == 128]
    assert all(np.isfinite(d) for d in deltas)


def test_scenario_config_validation(tmp_path: Path) -> None:
    good = _scenario_file(tmp_path, NOMINAL_SCENARIO, "good.yaml")
    assert load_scenario_config(good).evidence == "assumed"
    bad = _scenario_file(
        tmp_path, 'schema_version: "1.0"\nassumption_set: x\n', "bad.yaml"
    )
    with pytest.raises(InputContractError) as exc:
        load_scenario_config(bad)
    assert "evidence" in str(exc.value)
    unknown = _scenario_file(
        tmp_path,
        'schema_version: "1.0"\nassumption_set: x\nevidence: assumed\nsurprise: 1\n',
        "unknown.yaml",
    )
    with pytest.raises(InputContractError):
        load_scenario_config(unknown)
