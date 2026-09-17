"""Nominal (no process variation) simulation pipeline."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

from particle_sim.aggregation.critical_area import critical_area_um2, normalized_ca, um2_to_mm2
from particle_sim.aggregation.probability import ProbabilityAccumulator, probability_ca
from particle_sim.config import canonical_hash, file_sha256, load_manifest, load_scenario_config
from particle_sim.errors import InputContractError
from particle_sim.failure.union import ANY_MODE, union_by_mode
from particle_sim.geometry.contact import DistanceCache
from particle_sim.geometry.legacy import legacy_short_hotspot
from particle_sim.geometry.transforms import (
    apply_scenario_to_masks,
    apply_scenario_to_stack,
    scenario_manifest,
)
from particle_sim.io.masks import load_product_masks
from particle_sim.io.outputs import (
    CA_COLUMNS,
    SCHEMA_VERSION,
    UNION_CA_COLUMNS,
    HotspotStore,
    ProbabilityStore,
    build_run_manifest,
    make_run_id,
    write_json,
    write_table,
)
from particle_sim.io.stack import build_stack_z
from particle_sim.schema import scenario_warnings, todo_physics_warnings
from particle_sim.uncertainty.sampler import build_scenarios
from particle_sim.units import NOT_APPLICABLE, OPEN_MODE, evaluate_units, unit_ids, unit_modes
from particle_sim.validation.convergence import CONVERGENCE_COLUMNS, convergence_series
from particle_sim.validation.regression import explain_legacy_vs_edt

LEGACY_COMPARISON_COLUMNS = (
    "particle_size_um",
    "pair_id",
    "equal",
    "n_diff",
    "only_in_reference",
    "only_in_candidate",
    "kernel_quantization",
    "kernel_shape",
    "unexplained",
)


@dataclass(frozen=True)
class NominalRunResult:
    run_id: str
    output_dir: Path
    hotspots: Path
    ca_table: Path
    union_table: Path
    legacy_comparison: Path
    run_manifest: Path


def _require_empty_output(output_dir: Path, product_id: str) -> None:
    if output_dir.exists() and (not output_dir.is_dir() or any(output_dir.iterdir())):
        raise InputContractError(
            f"출력 경로가 비어 있지 않습니다(덮어쓰기 금지): {output_dir}",
            product_id=product_id,
            field="output",
        )


def run_nominal(
    manifest_path: Path, output_dir: Path, *, now: datetime | None = None
) -> NominalRunResult:
    """Run legacy + nominal 2.5D short simulation for every step, size and pair."""
    manifest_path = Path(manifest_path)
    output_dir = Path(output_dir)
    manifest = load_manifest(manifest_path)
    _require_empty_output(output_dir, manifest.product_id)
    masks = load_product_masks(manifest, manifest_path)
    stack_z = build_stack_z(manifest)
    created_at = now or datetime.now(UTC)
    run_id = make_run_id(manifest.product_id, created_at)

    grid = masks.grid
    electrodes = {e.electrode_id: e for e in manifest.electrodes}
    pairs = list(manifest.electrode_pairs)
    unit_id_list = unit_ids(manifest)
    unit_mode_list = unit_modes(manifest)
    steps = list(manifest.particle_source_steps)
    sizes = list(manifest.particle_sizes_um)
    modes: list[str] = []
    for mode in unit_mode_list:
        if mode not in modes:
            modes.append(mode)
    union_keys = [*modes, ANY_MODE]
    assumption_set = manifest.post_step_policy.mode

    from particle_sim import __version__

    output_dir.mkdir(parents=True, exist_ok=True)
    result = NominalRunResult(
        run_id=run_id,
        output_dir=output_dir,
        hotspots=output_dir / "hotspots.zarr",
        ca_table=output_dir / "ca.parquet",
        union_table=output_dir / "ca_union.parquet",
        legacy_comparison=output_dir / "legacy_vs_edt.csv",
        run_manifest=output_dir / "run_manifest.json",
    )
    store = HotspotStore(
        result.hotspots,
        source_steps=[s.step_id for s in steps],
        particle_sizes_um=sizes,
        pair_ids=unit_id_list,
        pair_modes=unit_mode_list,
        union_keys=union_keys,
        shape_yx=grid.shape_yx,
    )

    cache = DistanceCache(masks)
    ca_rows: list[dict] = []
    union_rows: list[dict] = []
    legacy_cache: dict[tuple[int, int], np.ndarray] = {}
    comparison_rows: list[dict] = []

    for size_index, diameter in enumerate(sizes):
        for pair_index, pair in enumerate(pairs):
            legacy = legacy_short_hotspot(
                masks.electrodes[pair.electrode_a],
                masks.electrodes[pair.electrode_b],
                diameter,
                grid,
                manifest.legacy_kernel,
            )
            legacy_cache[(size_index, pair_index)] = legacy
            cmp = explain_legacy_vs_edt(
                masks.electrodes[pair.electrode_a],
                masks.electrodes[pair.electrode_b],
                diameter,
                grid,
                manifest.legacy_kernel,
            )
            comparison_rows.append(
                {
                    "particle_size_um": diameter,
                    "pair_id": pair.pair_id,
                    "equal": cmp.equal,
                    "n_diff": cmp.n_diff,
                    "only_in_reference": cmp.only_in_reference,
                    "only_in_candidate": cmp.only_in_candidate,
                    "kernel_quantization": cmp.cause_counts.get("kernel_quantization", 0),
                    "kernel_shape": cmp.cause_counts.get("kernel_shape", 0),
                    "unexplained": cmp.cause_counts.get("unexplained", 0),
                }
            )

    pair_index_by_id = {unit_id: i for i, unit_id in enumerate(unit_id_list)}
    empty_legacy = np.zeros(grid.shape_yx, dtype=bool)
    for step_index, step in enumerate(steps):
        for size_index, diameter in enumerate(sizes):
            results = evaluate_units(manifest, masks, stack_z, step.step_id, diameter, cache)
            by_id = {r.pair_id: r for r in results}
            nominal_ca: dict[str, float] = {}
            for pair_index, unit_id in enumerate(unit_id_list):
                pair_result = by_id[unit_id]
                legacy = legacy_cache.get((size_index, pair_index), empty_legacy)
                is_open = pair_result.failure_mode == OPEN_MODE
                store.write_pair(
                    step_index,
                    size_index,
                    pair_index,
                    legacy=legacy,
                    nominal=pair_result.failure,
                    margin_um=pair_result.pair_margin_um,
                )
                ca_nominal = critical_area_um2(pair_result.failure, grid)
                nominal_ca[unit_id] = ca_nominal
                ea = electrodes[pair_result.electrode_a]
                eb = None if is_open else electrodes[pair_result.electrode_b]
                ca_rows.append(
                    {
                        "schema_version": SCHEMA_VERSION,
                        "run_id": run_id,
                        "product_id": manifest.product_id,
                        "design_version": manifest.design_version,
                        "source_step": step.step_id,
                        "fab_distribution_layer": step.fab_distribution_layer,
                        "particle_size_um": float(diameter),
                        "pair_id": unit_id,
                        "electrode_a": pair_result.electrode_a,
                        "electrode_b": pair_result.electrode_b,
                        "net_a": ea.net_id,
                        "net_b": NOT_APPLICABLE if eb is None else eb.net_id,
                        "failure_mode": pair_result.failure_mode,
                        "classification": pair_result.classification,
                        "ca_legacy_um2": (
                            float("nan") if is_open else critical_area_um2(legacy, grid)
                        ),
                        "ca_nominal_um2": ca_nominal,
                        "ca_nominal_mm2": um2_to_mm2(ca_nominal),
                        "ca_nominal_normalized": normalized_ca(ca_nominal, grid),
                        "ca_mean_um2": float("nan"),
                        "ca_p05_um2": float("nan"),
                        "ca_p50_um2": float("nan"),
                        "ca_p95_um2": float("nan"),
                        "scenario_count": 0,
                        "convergence_delta": float("nan"),
                        "assumption_set": assumption_set,
                        "simulator_version": __version__,
                    }
                )

            unions = union_by_mode(results)
            for union_index, key in enumerate(union_keys):
                union = unions[key]
                store.write_union(
                    step_index,
                    size_index,
                    union_index,
                    any_failure=union.any_failure,
                    pair_count=union.pair_count,
                )
                legacy_union = np.zeros(grid.shape_yx, dtype=bool)
                for pair_id in union.pair_ids:
                    legacy_key = (size_index, pair_index_by_id[pair_id])
                    if legacy_key in legacy_cache:
                        legacy_union |= legacy_cache[legacy_key]
                ca_union = critical_area_um2(union.any_failure, grid)
                union_rows.append(
                    {
                        "schema_version": SCHEMA_VERSION,
                        "run_id": run_id,
                        "product_id": manifest.product_id,
                        "design_version": manifest.design_version,
                        "source_step": step.step_id,
                        "fab_distribution_layer": step.fab_distribution_layer,
                        "particle_size_um": float(diameter),
                        "failure_mode": key,
                        "pair_ids": ";".join(union.pair_ids),
                        "ca_union_legacy_um2": critical_area_um2(legacy_union, grid),
                        "ca_union_nominal_um2": ca_union,
                        "ca_union_nominal_normalized": normalized_ca(ca_union, grid),
                        "ca_pair_sum_nominal_um2": float(
                            sum(nominal_ca[p] for p in union.pair_ids)
                        ),
                        "assumption_set": assumption_set,
                        "simulator_version": __version__,
                    }
                )

    write_table(ca_rows, CA_COLUMNS, result.ca_table)
    write_table(union_rows, UNION_CA_COLUMNS, result.union_table)
    with result.legacy_comparison.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(LEGACY_COMPARISON_COLUMNS))
        writer.writeheader()
        writer.writerows(comparison_rows)

    input_hashes = dict(masks.input_hashes)
    input_hashes["manifest"] = file_sha256(manifest_path)
    run_manifest = build_run_manifest(
        run_id=run_id,
        product_id=manifest.product_id,
        design_version=manifest.design_version,
        config_hash=canonical_hash(manifest.model_dump(mode="json")),
        input_hashes=input_hashes,
        particle_sizes_um=[float(s) for s in sizes],
        source_steps=[s.step_id for s in steps],
        assumption_set=assumption_set,
        warnings=todo_physics_warnings(manifest),
        legacy_kernel=manifest.legacy_kernel.model_dump(mode="json"),
        post_step_policy=manifest.post_step_policy.model_dump(mode="json"),
        created_at=created_at,
        git_commit=None,
    )
    write_json(run_manifest, result.run_manifest)
    return result


CONVERGENCE_PREFIXES: tuple[int, ...] = (64, 128, 256)


@dataclass(frozen=True)
class ProbabilityRunResult:
    run_id: str
    output_dir: Path
    probability: Path
    ca_table: Path
    union_table: Path
    convergence: Path
    run_manifest: Path
    scenario_count: int


def _prefix_means(values: list[float]) -> dict[int, float]:
    """Mean CA over the first n scenarios for the standard convergence prefixes."""
    result: dict[int, float] = {}
    array = np.asarray(values, dtype=np.float64)
    for prefix in CONVERGENCE_PREFIXES:
        if array.size >= prefix:
            result[prefix] = float(array[:prefix].mean())
    if array.size and array.size not in result:
        result[int(array.size)] = float(array.mean())
    return result


def run_probability(
    manifest_path: Path,
    scenario_path: Path,
    output_dir: Path,
    *,
    sobol_exponent_m: int,
    seed: int,
    now: datetime | None = None,
) -> ProbabilityRunResult:
    """Accumulate scenario hotspots into probability maps, CA percentiles and convergence."""
    manifest_path = Path(manifest_path)
    output_dir = Path(output_dir)
    manifest = load_manifest(manifest_path)
    _require_empty_output(output_dir, manifest.product_id)
    scenario_config = load_scenario_config(Path(scenario_path))
    masks = load_product_masks(manifest, manifest_path)
    stack_z = build_stack_z(manifest)
    created_at = now or datetime.now(UTC)
    run_id = make_run_id(manifest.product_id, created_at)
    scenarios = build_scenarios(scenario_config, sobol_exponent_m, seed)

    grid = masks.grid
    electrodes = {e.electrode_id: e for e in manifest.electrodes}
    unit_id_list = unit_ids(manifest)
    unit_mode_list = unit_modes(manifest)
    steps = list(manifest.particle_source_steps)
    sizes = list(manifest.particle_sizes_um)
    modes: list[str] = []
    for mode in unit_mode_list:
        if mode not in modes:
            modes.append(mode)
    union_keys = [*modes, ANY_MODE]

    from particle_sim import __version__

    output_dir.mkdir(parents=True, exist_ok=True)
    result = ProbabilityRunResult(
        run_id=run_id,
        output_dir=output_dir,
        probability=output_dir / "probability.zarr",
        ca_table=output_dir / "ca.parquet",
        union_table=output_dir / "ca_union.parquet",
        convergence=output_dir / "convergence.csv",
        run_manifest=output_dir / "run_manifest.json",
        scenario_count=len(scenarios),
    )
    store = ProbabilityStore(
        result.probability,
        source_steps=[s.step_id for s in steps],
        particle_sizes_um=sizes,
        unit_ids=unit_id_list,
        unit_modes=unit_mode_list,
        union_keys=union_keys,
        shape_yx=grid.shape_yx,
    )

    scenario_masks = [apply_scenario_to_masks(manifest, masks, s) for s in scenarios]
    scenario_stacks = [apply_scenario_to_stack(stack_z, s) for s in scenarios]
    scenario_manifests = [scenario_manifest(manifest, s) for s in scenarios]
    nominal_cache = DistanceCache(masks)
    ca_rows: list[dict] = []
    union_rows: list[dict] = []
    convergence_rows: list[dict] = []
    assumption_set = scenario_config.assumption_set

    for step_index, step in enumerate(steps):
        for size_index, diameter in enumerate(sizes):
            unit_accumulators = {
                unit_id: ProbabilityAccumulator(grid.shape_yx, grid) for unit_id in unit_id_list
            }
            union_accumulators = {
                key: ProbabilityAccumulator(grid.shape_yx, grid) for key in union_keys
            }
            for scenario, s_masks, s_stack, s_manifest in zip(
                scenarios, scenario_masks, scenario_stacks, scenario_manifests, strict=True
            ):
                diameter_um = max(diameter + scenario.particle_diameter_delta_um, 1e-9)
                results = evaluate_units(
                    s_manifest,
                    s_masks,
                    s_stack,
                    step.step_id,
                    diameter_um,
                    DistanceCache(s_masks),
                )
                for unit in results:
                    unit_accumulators[unit.pair_id].add(unit.failure, unit.pair_margin_um)
                for key, union in union_by_mode(results).items():
                    union_accumulators[key].add(union.any_failure)

            nominal_results = {
                r.pair_id: r
                for r in evaluate_units(
                    manifest, masks, stack_z, step.step_id, diameter, nominal_cache
                )
            }
            for unit_index, unit_id in enumerate(unit_id_list):
                accumulator = unit_accumulators[unit_id]
                probability_map = accumulator.probability()
                store.write_unit(
                    step_index,
                    size_index,
                    unit_index,
                    probability=probability_map,
                    margin_mean_um=accumulator.mean_margin_um(),
                )
                summary = probability_ca(accumulator)
                nominal = nominal_results[unit_id]
                is_open = nominal.failure_mode == OPEN_MODE
                ea = electrodes[nominal.electrode_a]
                eb = None if is_open else electrodes[nominal.electrode_b]
                ca_nominal = critical_area_um2(nominal.failure, grid)
                prefixes = _prefix_means(accumulator.scenario_ca_um2)
                steps_series = convergence_series(prefixes)
                for entry in steps_series:
                    convergence_rows.append(
                        {
                            "source_step": step.step_id,
                            "particle_size_um": float(diameter),
                            "pair_id": unit_id,
                            "failure_mode": nominal.failure_mode,
                            "scenario_count": entry.scenario_count,
                            "ca_um2": entry.ca_um2,
                            "delta_previous": entry.delta_previous,
                        }
                    )
                ca_rows.append(
                    {
                        "schema_version": SCHEMA_VERSION,
                        "run_id": run_id,
                        "product_id": manifest.product_id,
                        "design_version": manifest.design_version,
                        "source_step": step.step_id,
                        "fab_distribution_layer": step.fab_distribution_layer,
                        "particle_size_um": float(diameter),
                        "pair_id": unit_id,
                        "electrode_a": nominal.electrode_a,
                        "electrode_b": nominal.electrode_b,
                        "net_a": ea.net_id,
                        "net_b": NOT_APPLICABLE if eb is None else eb.net_id,
                        "failure_mode": nominal.failure_mode,
                        "classification": nominal.classification,
                        "ca_legacy_um2": float("nan"),
                        "ca_nominal_um2": ca_nominal,
                        "ca_nominal_mm2": um2_to_mm2(ca_nominal),
                        "ca_nominal_normalized": normalized_ca(ca_nominal, grid),
                        "ca_mean_um2": summary.ca_probability_mean_um2,
                        "ca_p05_um2": summary.ca_p05_um2,
                        "ca_p50_um2": summary.ca_p50_um2,
                        "ca_p95_um2": summary.ca_p95_um2,
                        "scenario_count": summary.scenario_count,
                        "convergence_delta": (
                            steps_series[-1].delta_previous if steps_series else float("nan")
                        ),
                        "assumption_set": assumption_set,
                        "simulator_version": __version__,
                    }
                )

            for union_index, key in enumerate(union_keys):
                accumulator = union_accumulators[key]
                union_probability = accumulator.probability()
                store.write_union(
                    step_index, size_index, union_index, probability=union_probability
                )
                summary = probability_ca(accumulator)
                member_ids = [
                    unit_id
                    for unit_id, mode in zip(unit_id_list, unit_mode_list, strict=True)
                    if key in (ANY_MODE, mode)
                ]
                union_rows.append(
                    {
                        "schema_version": SCHEMA_VERSION,
                        "run_id": run_id,
                        "product_id": manifest.product_id,
                        "design_version": manifest.design_version,
                        "source_step": step.step_id,
                        "fab_distribution_layer": step.fab_distribution_layer,
                        "particle_size_um": float(diameter),
                        "failure_mode": key,
                        "pair_ids": ";".join(member_ids),
                        "ca_union_legacy_um2": float("nan"),
                        "ca_union_nominal_um2": summary.ca_probability_mean_um2,
                        "ca_union_nominal_normalized": normalized_ca(
                            summary.ca_probability_mean_um2, grid
                        ),
                        "ca_pair_sum_nominal_um2": float(
                            sum(
                                probability_ca(unit_accumulators[unit_id])
                                .ca_probability_mean_um2
                                for unit_id in member_ids
                            )
                        ),
                        "assumption_set": assumption_set,
                        "simulator_version": __version__,
                    }
                )

    write_table(ca_rows, CA_COLUMNS, result.ca_table)
    write_table(union_rows, UNION_CA_COLUMNS, result.union_table)
    with result.convergence.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(CONVERGENCE_COLUMNS))
        writer.writeheader()
        writer.writerows(convergence_rows)

    input_hashes = dict(masks.input_hashes)
    input_hashes["manifest"] = file_sha256(manifest_path)
    input_hashes["scenario"] = file_sha256(Path(scenario_path))
    run_manifest = build_run_manifest(
        run_id=run_id,
        product_id=manifest.product_id,
        design_version=manifest.design_version,
        config_hash=canonical_hash(
            {
                "manifest": manifest.model_dump(mode="json"),
                "scenario": scenario_config.model_dump(mode="json"),
            }
        ),
        input_hashes=input_hashes,
        particle_sizes_um=[float(s) for s in sizes],
        source_steps=[s.step_id for s in steps],
        assumption_set=assumption_set,
        warnings=[*todo_physics_warnings(manifest), *scenario_warnings(scenario_config)],
        legacy_kernel=manifest.legacy_kernel.model_dump(mode="json"),
        post_step_policy=manifest.post_step_policy.model_dump(mode="json"),
        created_at=created_at,
        git_commit=None,
    )
    run_manifest["random_seed"] = seed
    run_manifest["sobol_exponent_m"] = sobol_exponent_m
    run_manifest["scenario_count"] = len(scenarios)
    write_json(run_manifest, result.run_manifest)
    return result
