"""Nominal (no process variation) simulation pipeline."""

from __future__ import annotations

import csv
import tempfile
from dataclasses import dataclass, replace
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

from particle_sim.aggregation.critical_area import critical_area_um2, normalized_ca, um2_to_mm2
from particle_sim.aggregation.probability import ProbabilityAccumulator, probability_ca
from particle_sim.aggregation.risk_allocation import allocate_risk
from particle_sim.config import (
    canonical_hash,
    file_sha256,
    load_manifest,
    load_scenario_config,
    validate_scenario_for_manifest,
)
from particle_sim.errors import InputContractError
from particle_sim.failure.auto_open import INSUFFICIENT_TERMINALS
from particle_sim.failure.contact_not_open import CNO_COVERAGE_RULE, CONTACT_LOST
from particle_sim.failure.open_failure import EXCLUDED_CLASSIFICATIONS, check_open_baselines
from particle_sim.failure.union import ANY_MODE, union_by_mode
from particle_sim.geometry.contact import DistanceCache
from particle_sim.geometry.contacts import build_contact_topology
from particle_sim.geometry.legacy import legacy_short_hotspot
from particle_sim.geometry.transforms import (
    apply_scenario_to_masks,
    apply_scenario_to_stack,
    scenario_manifest,
)
from particle_sim.io.masks import MaskSet, load_product_masks
from particle_sim.io.outputs import (
    CA_COLUMNS,
    PROBABILITY_SEMANTICS,
    SCHEMA_VERSION,
    UNION_CA_COLUMNS,
    HotspotStore,
    ProbabilityStore,
    artifact_records,
    build_run_manifest,
    make_run_id,
    write_json,
    write_table,
)
from particle_sim.io.stack import build_stack_z
from particle_sim.schema import scenario_warnings, todo_physics_warnings
from particle_sim.uncertainty.sampler import build_scenarios
from particle_sim.units import (
    KIND_CONTACT,
    KIND_OPEN_ELECTRODE,
    KIND_OPEN_PATH,
    KIND_PAIR,
    NOT_APPLICABLE,
    OpenResultCache,
    UnitSpec,
    evaluate_units,
    unit_specs,
)
from particle_sim.validation.convergence import (
    CONVERGENCE_COLUMNS,
    convergence_series,
    is_converged,
)
from particle_sim.validation.regression import explain_legacy_vs_edt
from particle_sim.validation.short_feasibility import structural_zero_short_warnings

# scenario classifications that are excluded from particle-caused CA and must be reported
SCENARIO_EXCLUSIONS = (*EXCLUDED_CLASSIFICATIONS, CONTACT_LOST, INSUFFICIENT_TERMINALS)


def _unit_identity(spec: UnitSpec, electrodes: dict) -> dict:
    """ID and electrode columns shared by the nominal and probability CA rows."""
    ea = electrodes[spec.electrode_a]
    eb = electrodes.get(spec.electrode_b)
    is_open = spec.unit_kind in (KIND_OPEN_PATH, KIND_OPEN_ELECTRODE)
    return {
        "failure_unit_id": spec.unit_id,
        "unit_kind": spec.unit_kind,
        "pair_id": spec.unit_id,
        "electrode_pair_id": spec.unit_id if spec.unit_kind == KIND_PAIR else NOT_APPLICABLE,
        "open_path_id": spec.unit_id if is_open else NOT_APPLICABLE,
        "contact_id": spec.unit_id if spec.unit_kind == KIND_CONTACT else NOT_APPLICABLE,
        "electrode_a": spec.electrode_a,
        "electrode_b": spec.electrode_b,
        "electrode_a_layer": ea.layer_name,
        "electrode_b_layer": NOT_APPLICABLE if eb is None else eb.layer_name,
        "net_a": ea.net_id,
        "net_b": NOT_APPLICABLE if eb is None else eb.net_id,
        "failure_mode": spec.failure_mode,
    }

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


def _unpack_legacy_map(packed: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    """Restore a cached legacy map, excluding padding bits in the final byte."""
    height, width = shape
    return np.unpackbits(packed, count=height * width).reshape(height, width).astype(bool)


def _require_empty_output(output_dir: Path, product_id: str) -> None:
    if output_dir.exists() and (not output_dir.is_dir() or any(output_dir.iterdir())):
        raise InputContractError(
            f"출력 경로가 비어 있지 않습니다(덮어쓰기 금지): {output_dir}",
            product_id=product_id,
            field="output",
        )


def final_convergence_delta(series_by_seed: list[list]) -> float:
    """Worst last comparison over the seeds, matching ``is_converged`` (last prefix step).

    Each seed contributes the relative change of its final prefix step (e.g. 128->256). An
    infinite change (CA dropped to zero) is kept, not discarded, so it cannot look converged.
    NaN means no seed had two prefixes to compare.
    """
    finals = [
        series[-1].delta_previous
        for series in series_by_seed
        if len(series) >= 2 and not np.isnan(series[-1].delta_previous)
    ]
    return float(max(finals)) if finals else float("nan")


def _distribution(values: list[float]) -> tuple[float, float, float, float]:
    array = np.asarray(values, dtype=np.float64)
    if array.size == 0:
        return (float("nan"),) * 4
    return (
        float(array.mean()),
        float(np.percentile(array, 5)),
        float(np.percentile(array, 50)),
        float(np.percentile(array, 95)),
    )


class _ScenarioMaskProvider:
    """작은 실행은 메모리, 큰 실행은 packbits 임시 파일로 변형 mask를 재사용한다."""

    MEMORY_LIMIT_BYTES = 256 * 1024 * 1024

    def __init__(self, manifest, masks: MaskSet, scenarios, output_dir: Path) -> None:
        self._base = masks
        self._manifest = manifest
        self._contact_names = list(masks.contact_masks or {})
        self._layer_names = list(masks.layers)
        self._electrode_names = list(masks.electrodes)
        self._terminal_names = list(masks.terminal_regions or {})
        mask_count = (
            len(self._layer_names)
            + len(self._electrode_names)
            + 2 * len(self._terminal_names)
            + len(self._contact_names)
        )
        estimated = len(scenarios) * mask_count * int(np.prod(masks.grid.shape_yx))
        self.storage = "memory" if estimated <= self.MEMORY_LIMIT_BYTES else "packed_temp"
        self.estimated_eager_bytes = estimated
        self._memory: list[MaskSet] = []
        self._temporary: tempfile.TemporaryDirectory[str] | None = None
        self._paths: list[Path] = []
        if self.storage == "memory":
            self._memory = [apply_scenario_to_masks(manifest, masks, item) for item in scenarios]
            return
        self._temporary = tempfile.TemporaryDirectory(prefix="particle-sim-masks-", dir=output_dir)
        try:
            self._pack(manifest, masks, scenarios, Path(self._temporary.name))
        except BaseException:
            self.close()
            raise

    def _pack(self, manifest, masks: MaskSet, scenarios, directory: Path) -> None:
        for index, scenario in enumerate(scenarios):
            transformed = apply_scenario_to_masks(manifest, masks, scenario)
            arrays = [transformed.layers[name] for name in self._layer_names]
            arrays.extend(transformed.electrodes[name] for name in self._electrode_names)
            for name in self._terminal_names:
                arrays.extend((transformed.terminal_regions or {})[name])
            arrays.extend((transformed.contact_masks or {})[name] for name in self._contact_names)
            packed = np.stack([np.packbits(array.reshape(-1)) for array in arrays], axis=0)
            path = directory / f"{index:06d}.npy"
            np.save(path, packed, allow_pickle=False)
            self._paths.append(path)

    def get(self, index: int) -> MaskSet:
        if self.storage == "memory":
            return self._memory[index]
        packed = np.load(self._paths[index], allow_pickle=False)
        shape = self._base.grid.shape_yx
        size = int(np.prod(shape))
        arrays = [np.unpackbits(row, count=size).reshape(shape).astype(bool) for row in packed]
        position = 0
        layers = {}
        for name in self._layer_names:
            layers[name] = arrays[position]
            position += 1
        electrodes = {}
        for name in self._electrode_names:
            electrodes[name] = arrays[position]
            position += 1
        terminal_regions = {}
        for name in self._terminal_names:
            terminal_regions[name] = (arrays[position], arrays[position + 1])
            position += 2
        contact_masks = {}
        for name in self._contact_names:
            contact_masks[name] = arrays[position]
            position += 1
        return MaskSet(
            grid=self._base.grid,
            layers=layers,
            electrodes=electrodes,
            input_hashes=self._base.input_hashes,
            terminal_regions=terminal_regions,
            layer_unassigned_pixels=self._base.layer_unassigned_pixels,
            contact_masks=contact_masks,
            contact_topology=build_contact_topology(
                self._manifest, contact_masks, electrodes, strict=False
            ),
        )

    def close(self) -> None:
        if self._temporary is not None:
            self._temporary.cleanup()


def run_nominal(
    manifest_path: Path, output_dir: Path, *, now: datetime | None = None
) -> NominalRunResult:
    """Run legacy + nominal 2.5D short simulation for every step, size and pair."""
    manifest_path = Path(manifest_path)
    output_dir = Path(output_dir)
    manifest = load_manifest(manifest_path)
    _require_empty_output(output_dir, manifest.product_id)
    masks = load_product_masks(manifest, manifest_path)
    check_open_baselines(manifest, masks)
    stack_z = build_stack_z(manifest)
    created_at = now or datetime.now(UTC)
    run_id = make_run_id(manifest.product_id, created_at)

    grid = masks.grid
    electrodes = {e.electrode_id: e for e in manifest.electrodes}
    pairs = list(manifest.electrode_pairs)
    specs = unit_specs(manifest, masks)
    unit_id_list = [spec.unit_id for spec in specs]
    unit_mode_list = [spec.failure_mode for spec in specs]
    spec_by_id = {spec.unit_id: spec for spec in specs}
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
        unit_kinds=[spec.unit_kind for spec in specs],
        union_keys=union_keys,
        shape_yx=grid.shape_yx,
    )

    cache = DistanceCache(masks)
    open_cache: OpenResultCache = {}
    ca_rows: list[dict] = []
    union_rows: list[dict] = []
    legacy_cache: dict[tuple[int, int], tuple[np.ndarray, tuple[int, int]]] = {}
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
            legacy_cache[(size_index, pair_index)] = (np.packbits(legacy), legacy.shape)
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
    last_step_by_unit = {
        spec.unit_id: max(
            (
                index
                for index, step in enumerate(steps)
                if spec.applicable_steps is None or step.step_id in spec.applicable_steps
            ),
            default=-1,
        )
        for spec in specs
        if spec.unit_kind != KIND_PAIR
    }
    contact_margin = float(manifest.contact_margin_um)
    for step_index, step in enumerate(steps):
        for size_index, diameter in enumerate(sizes):
            results = evaluate_units(
                manifest,
                masks,
                stack_z,
                step.step_id,
                diameter,
                cache,
                open_cache,
            )
            by_id = {r.pair_id: r for r in results}
            allocation = allocate_risk(results, grid)
            nominal_ca: dict[str, float] = {}
            for pair_index, unit_id in enumerate(unit_id_list):
                pair_result = by_id[unit_id]
                packed_legacy = legacy_cache.get((size_index, pair_index))
                legacy = (
                    _unpack_legacy_map(*packed_legacy)
                    if packed_legacy is not None
                    else empty_legacy
                )
                spec = spec_by_id[unit_id]
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
                ca_rows.append(
                    {
                        "schema_version": SCHEMA_VERSION,
                        "run_id": run_id,
                        "product_id": manifest.product_id,
                        "design_version": manifest.design_version,
                        "source_step": step.step_id,
                        "fab_distribution_layer": step.fab_distribution_layer,
                        "particle_size_um": float(diameter),
                        **_unit_identity(spec, electrodes),
                        "classification": pair_result.classification,
                        "ca_legacy_um2": (
                            critical_area_um2(legacy, grid)
                            if spec.unit_kind == KIND_PAIR
                            else float("nan")
                        ),
                        "ca_nominal_um2": ca_nominal,
                        "ca_allocated_nominal_um2": allocation.allocated_ca_um2[unit_id],
                        "ca_unique_nominal_um2": allocation.unique_ca_um2[unit_id],
                        "ca_nominal_mm2": um2_to_mm2(ca_nominal),
                        "ca_nominal_normalized": normalized_ca(ca_nominal, grid),
                        "ca_mean_um2": float("nan"),
                        "ca_p05_um2": float("nan"),
                        "ca_p50_um2": float("nan"),
                        "ca_p95_um2": float("nan"),
                        "ca_allocated_mean_um2": float("nan"),
                        "ca_allocated_p05_um2": float("nan"),
                        "ca_allocated_p50_um2": float("nan"),
                        "ca_allocated_p95_um2": float("nan"),
                        "ca_unique_mean_um2": float("nan"),
                        "ca_unique_p05_um2": float("nan"),
                        "ca_unique_p50_um2": float("nan"),
                        "ca_unique_p95_um2": float("nan"),
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
                    packed_legacy = legacy_cache.get(legacy_key)
                    if packed_legacy is not None:
                        legacy_union |= _unpack_legacy_map(*packed_legacy)
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
                        "union_scope": "all_failure_modes" if key == ANY_MODE else "failure_mode",
                        "pair_ids": ";".join(union.pair_ids),
                        "failure_unit_ids": ";".join(union.pair_ids),
                        "ca_union_legacy_um2": critical_area_um2(legacy_union, grid),
                        "ca_union_nominal_um2": ca_union,
                        "ca_union_nominal_normalized": normalized_ca(ca_union, grid),
                        "ca_pair_sum_nominal_um2": float(
                            sum(nominal_ca[p] for p in union.pair_ids)
                        ),
                        "ca_union_mean_um2": float("nan"),
                        "ca_union_p05_um2": float("nan"),
                        "ca_union_p50_um2": float("nan"),
                        "ca_union_p95_um2": float("nan"),
                        "ca_raw_unit_sum_mean_um2": float("nan"),
                        "scenario_count": 0,
                        "convergence_delta": float("nan"),
                        "assumption_set": assumption_set,
                        "simulator_version": __version__,
                    }
                )

            # Results remain live in ``results`` through writes and unions. Once this
            # step has consumed a unit's final applicable result at this diameter,
            # its cache entry cannot be reused by the remaining step-major traversal.
            for unit_id, last_step_index in last_step_by_unit.items():
                if step_index == last_step_index:
                    open_cache.pop((unit_id, float(diameter), contact_margin), None)

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
        warnings=[*todo_physics_warnings(manifest),
                  *structural_zero_short_warnings(manifest, stack_z)],
        legacy_kernel=manifest.legacy_kernel.model_dump(mode="json"),
        post_step_policy=manifest.post_step_policy.model_dump(mode="json"),
        created_at=created_at,
        git_commit=None,
    )
    run_manifest["layer_unassigned_pixels"] = dict(masks.layer_unassigned_pixels or {})
    if manifest.contacts:
        run_manifest["cno_coverage_rule"] = CNO_COVERAGE_RULE
    run_manifest["input_schema_versions"] = {"product_manifest": manifest.schema_version}
    run_manifest["shape_yx"] = list(grid.shape_yx)
    run_manifest["pixel_size_x_um"] = grid.pixel_size_x_um
    run_manifest["pixel_size_y_um"] = grid.pixel_size_y_um
    run_manifest["artifacts"] = artifact_records(
        [result.hotspots, result.ca_table, result.union_table, result.legacy_comparison],
        base=output_dir,
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
    additional_seeds: list[int] | None = None,
    convergence_tolerance: float | None = None,
    now: datetime | None = None,
) -> ProbabilityRunResult:
    """Accumulate scenario hotspots into probability maps, CA percentiles and convergence."""
    manifest_path = Path(manifest_path)
    output_dir = Path(output_dir)
    manifest = load_manifest(manifest_path)
    _require_empty_output(output_dir, manifest.product_id)
    scenario_config = load_scenario_config(Path(scenario_path))
    validate_scenario_for_manifest(
        scenario_config, manifest, scenario_path=Path(scenario_path)
    )
    masks = load_product_masks(manifest, manifest_path)
    check_open_baselines(manifest, masks)
    stack_z = build_stack_z(manifest)
    created_at = now or datetime.now(UTC)
    run_id = make_run_id(manifest.product_id, created_at)
    seed_values = list(dict.fromkeys([seed, *(additional_seeds or [])]))
    scenarios = []
    for seed_value in seed_values:
        for scenario in build_scenarios(scenario_config, sobol_exponent_m, seed_value):
            scenarios.append(replace(scenario, scenario_id=len(scenarios)))

    grid = masks.grid
    electrodes = {e.electrode_id: e for e in manifest.electrodes}
    specs = unit_specs(manifest, masks)
    unit_id_list = [spec.unit_id for spec in specs]
    unit_mode_list = [spec.failure_mode for spec in specs]
    spec_by_id = {spec.unit_id: spec for spec in specs}
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
        unit_kinds=[spec.unit_kind for spec in specs],
        union_keys=union_keys,
        shape_yx=grid.shape_yx,
    )

    mask_provider = _ScenarioMaskProvider(manifest, masks, scenarios, output_dir)
    try:
        scenario_stacks = [apply_scenario_to_stack(stack_z, scenario) for scenario in scenarios]
        scenario_manifests = [scenario_manifest(manifest, scenario) for scenario in scenarios]
        nominal_cache = DistanceCache(masks)
        nominal_open_cache: OpenResultCache = {}
        scenario_exclusions: dict[str, dict[str, int]] = {}
        lost_holes_by_scenario: dict[int, int] = {}
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
                nominal_results = {
                    r.pair_id: r
                    for r in evaluate_units(
                        manifest,
                        masks,
                        stack_z,
                        step.step_id,
                        diameter,
                        nominal_cache,
                        nominal_open_cache,
                    )
                }
                allocated_samples = {unit_id: [] for unit_id in unit_id_list}
                unique_samples = {unit_id: [] for unit_id in unit_id_list}
                for scenario_index, scenario in enumerate(scenarios):
                    s_masks = mask_provider.get(scenario_index)
                    if s_masks.contact_topology is not None:
                        lost_holes_by_scenario[scenario_index] = s_masks.contact_topology.lost_holes
                    s_stack = scenario_stacks[scenario_index]
                    s_manifest = scenario_manifests[scenario_index]
                    diameter_um = max(diameter + scenario.particle_diameter_delta_um, 1e-9)
                    results = evaluate_units(
                        s_manifest,
                        s_masks,
                        s_stack,
                        step.step_id,
                        diameter_um,
                        DistanceCache(s_masks),
                        specs=specs,
                    )
                    for unit in results:
                        unit_accumulators[unit.pair_id].add(unit.failure, unit.pair_margin_um)
                        # a unit already excluded in the nominal geometry (e.g. an electrode
                        # without two terminals) is structural, not a scenario effect
                        if (
                            unit.classification in SCENARIO_EXCLUSIONS
                            and nominal_results[unit.pair_id].classification
                            not in SCENARIO_EXCLUSIONS
                        ):
                            key = f"{step.step_id}|{float(diameter)}|{unit.pair_id}"
                            counts = scenario_exclusions.setdefault(key, {})
                            counts[unit.classification] = counts.get(unit.classification, 0) + 1
                    allocation = allocate_risk(results, grid)
                    for unit_id in unit_id_list:
                        allocated_samples[unit_id].append(allocation.allocated_ca_um2[unit_id])
                        unique_samples[unit_id].append(allocation.unique_ca_um2[unit_id])
                    for key, union in union_by_mode(results).items():
                        union_accumulators[key].add(union.any_failure)

                nominal_allocation = allocate_risk(list(nominal_results.values()), grid)
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
                    spec = spec_by_id[unit_id]
                    ca_nominal = critical_area_um2(nominal.failure, grid)
                    allocated_mean, allocated_p05, allocated_p50, allocated_p95 = _distribution(
                        allocated_samples[unit_id]
                    )
                    unique_mean, unique_p05, unique_p50, unique_p95 = _distribution(
                        unique_samples[unit_id]
                    )
                    offset = 0
                    unit_steps_by_seed: list = []
                    for seed_value in seed_values:
                        seed_count = sum(item.scramble_seed == seed_value for item in scenarios)
                        seed_samples = accumulator.scenario_ca_um2[offset : offset + seed_count]
                        offset += seed_count
                        seed_steps = convergence_series(_prefix_means(seed_samples))
                        unit_steps_by_seed.append(seed_steps)
                        for entry in seed_steps:
                            convergence_rows.append(
                                {
                                    "scramble_seed": seed_value,
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
                            **_unit_identity(spec, electrodes),
                            "classification": nominal.classification,
                            "ca_legacy_um2": float("nan"),
                            "ca_nominal_um2": ca_nominal,
                            "ca_allocated_nominal_um2": nominal_allocation.allocated_ca_um2[unit_id],
                            "ca_unique_nominal_um2": nominal_allocation.unique_ca_um2[unit_id],
                            "ca_nominal_mm2": um2_to_mm2(ca_nominal),
                            "ca_nominal_normalized": normalized_ca(ca_nominal, grid),
                            "ca_mean_um2": summary.ca_probability_mean_um2,
                            "ca_p05_um2": summary.ca_p05_um2,
                            "ca_p50_um2": summary.ca_p50_um2,
                            "ca_p95_um2": summary.ca_p95_um2,
                            "ca_allocated_mean_um2": allocated_mean,
                            "ca_allocated_p05_um2": allocated_p05,
                            "ca_allocated_p50_um2": allocated_p50,
                            "ca_allocated_p95_um2": allocated_p95,
                            "ca_unique_mean_um2": unique_mean,
                            "ca_unique_p05_um2": unique_p05,
                            "ca_unique_p50_um2": unique_p50,
                            "ca_unique_p95_um2": unique_p95,
                            "scenario_count": summary.scenario_count,
                            "convergence_delta": final_convergence_delta(unit_steps_by_seed),
                            "assumption_set": assumption_set,
                            "simulator_version": __version__,
                        }
                    )

                for union_index, key in enumerate(union_keys):
                    accumulator = union_accumulators[key]
                    union_probability = accumulator.probability()
                    offset = 0
                    union_steps_by_seed: list = []
                    for seed_value in seed_values:
                        seed_count = sum(item.scramble_seed == seed_value for item in scenarios)
                        seed_samples = accumulator.scenario_ca_um2[offset : offset + seed_count]
                        offset += seed_count
                        seed_steps = convergence_series(_prefix_means(seed_samples))
                        union_steps_by_seed.append(seed_steps)
                        for entry in seed_steps:
                            convergence_rows.append(
                                {
                                    "scramble_seed": seed_value,
                                    "source_step": step.step_id,
                                    "particle_size_um": float(diameter),
                                    "pair_id": f"union:{key}",
                                    "failure_mode": key,
                                    "scenario_count": entry.scenario_count,
                                    "ca_um2": entry.ca_um2,
                                    "delta_previous": entry.delta_previous,
                                }
                            )
                    store.write_union(
                        step_index, size_index, union_index, probability=union_probability
                    )
                    summary = probability_ca(accumulator)
                    nominal_union = union_by_mode(list(nominal_results.values()))[key]
                    ca_union_nominal = critical_area_um2(nominal_union.any_failure, grid)
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
                            "union_scope": "all_failure_modes" if key == ANY_MODE else "failure_mode",
                            "pair_ids": ";".join(member_ids),
                            "failure_unit_ids": ";".join(member_ids),
                            "ca_union_legacy_um2": float("nan"),
                            "ca_union_nominal_um2": ca_union_nominal,
                            "ca_union_nominal_normalized": normalized_ca(
                                ca_union_nominal, grid
                            ),
                            "ca_pair_sum_nominal_um2": float(
                                sum(
                                    critical_area_um2(nominal_results[unit_id].failure, grid)
                                    for unit_id in member_ids
                                )
                            ),
                            "ca_union_mean_um2": summary.ca_probability_mean_um2,
                            "ca_union_p05_um2": summary.ca_p05_um2,
                            "ca_union_p50_um2": summary.ca_p50_um2,
                            "ca_union_p95_um2": summary.ca_p95_um2,
                            "ca_raw_unit_sum_mean_um2": float(
                                sum(
                                    probability_ca(unit_accumulators[unit_id])
                                    .ca_probability_mean_um2
                                    for unit_id in member_ids
                                )
                            ),
                            "scenario_count": summary.scenario_count,
                            "convergence_delta": final_convergence_delta(union_steps_by_seed),
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
            warnings=[*todo_physics_warnings(manifest),
                      *scenario_warnings(scenario_config),
                      *structural_zero_short_warnings(manifest, stack_z)],
            legacy_kernel=manifest.legacy_kernel.model_dump(mode="json"),
            post_step_policy=manifest.post_step_policy.model_dump(mode="json"),
            created_at=created_at,
            git_commit=None,
        )
        run_manifest["random_seed"] = seed
        run_manifest["random_seeds"] = seed_values
        run_manifest["sobol_exponent_m"] = sobol_exponent_m
        run_manifest["scenario_count"] = len(scenarios)
        run_manifest["value_semantics"] = PROBABILITY_SEMANTICS
        run_manifest["scenario_mask_storage"] = mask_provider.storage
        run_manifest["scenario_exclusions"] = [
            {
                "source_step": key.split("|")[0],
                "particle_size_um": float(key.split("|")[1]),
                "failure_unit_id": key.split("|", 2)[2],
                "classification_counts": counts,
            }
            for key, counts in sorted(scenario_exclusions.items())
        ]
        run_manifest["scenario_contact_lost_holes"] = {
            "scenarios_with_lost_holes": sum(1 for n in lost_holes_by_scenario.values() if n),
            "max_lost_holes": max(lost_holes_by_scenario.values(), default=0),
        }
        if scenario_exclusions:
            run_manifest["warnings"] = [
                *run_manifest.get("warnings", []),
                "SCENARIO_EXCLUSION: 기준 단선·단자 소실·contact 소실 조건은 이물 기인 CA에서 "
                "제외했습니다. scenario_exclusions를 확인하세요",
            ]
        run_manifest["estimated_eager_mask_bytes"] = mask_provider.estimated_eager_bytes
        run_manifest["convergence_tolerance"] = convergence_tolerance
        if convergence_tolerance is not None:
            by_series: dict[tuple, dict[int, float]] = {}
            for row in convergence_rows:
                key = (
                    row["scramble_seed"],
                    row["source_step"],
                    row["particle_size_um"],
                    row["pair_id"],
                )
                by_series.setdefault(key, {})[int(row["scenario_count"])] = float(row["ca_um2"])
            not_converged = []
            undecidable = 0
            for key, samples in by_series.items():
                steps_for_key = convergence_series(samples)
                comparable = [s for s in steps_for_key if np.isfinite(s.delta_previous)]
                if not comparable:
                    # a single scenario count cannot be compared with anything
                    undecidable += 1
                    continue
                if not is_converged(steps_for_key, convergence_tolerance):
                    not_converged.append(
                        {
                            "scramble_seed": key[0],
                            "source_step": key[1],
                            "particle_size_um": key[2],
                            "pair_id": key[3],
                        }
                    )
            run_manifest["convergence_not_met"] = not_converged
            run_manifest["convergence_undecidable"] = undecidable
            if undecidable:
                run_manifest["warnings"] = [
                    *run_manifest["warnings"],
                    f"CONVERGENCE: 비교할 시나리오 수가 하나뿐이라 판정할 수 없는 계열 "
                    f"{undecidable}개가 있습니다(sobol-m을 7 이상으로 올리십시오)",
                ]
            if not_converged:
                run_manifest["warnings"] = [
                    *run_manifest["warnings"],
                    f"CONVERGENCE: 허용 오차 {convergence_tolerance}를 만족하지 못한 계열 "
                    f"{len(not_converged)}개가 있습니다",
                ]
        run_manifest["layer_unassigned_pixels"] = dict(masks.layer_unassigned_pixels or {})
        if manifest.contacts:
            run_manifest["cno_coverage_rule"] = CNO_COVERAGE_RULE
        run_manifest["input_schema_versions"] = {
            "product_manifest": manifest.schema_version,
            "scenario": scenario_config.schema_version,
        }
        run_manifest["scenario_evidence"] = scenario_config.model_dump(mode="json")
        run_manifest["shape_yx"] = list(grid.shape_yx)
        run_manifest["pixel_size_x_um"] = grid.pixel_size_x_um
        run_manifest["pixel_size_y_um"] = grid.pixel_size_y_um
        run_manifest["artifacts"] = artifact_records(
            [result.probability, result.ca_table, result.union_table, result.convergence],
            base=output_dir,
        )
        write_json(run_manifest, result.run_manifest)
    finally:
        # the packed masks live under output_dir; never leave them behind on failure
        mask_provider.close()
    return result
