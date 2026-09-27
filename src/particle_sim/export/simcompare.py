"""Export a nominal run into the simcompare example adapter data_dir format."""

from __future__ import annotations

import csv
import json
import re
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from particle_sim.config import load_manifest
from particle_sim.errors import InputContractError
from particle_sim.failure.auto_open import INSUFFICIENT_TERMINALS
from particle_sim.failure.open_failure import BASELINE_OPEN, TERMINAL_LOST
from particle_sim.io.masks import load_product_masks
from particle_sim.io.outputs import check_run_manifest, open_hotspot_store
from particle_sim.io.stack import build_stack_z
from particle_sim.units import NOT_APPLICABLE_STEP, unit_specs

BUNDLE_COLUMNS = (
    "sample_id", "product_id", "product_family", "design_version", "design_lineage_id",
    "electrode_a_role", "electrode_b_role", "layer_a", "layer_b", "failure_type", "fab_id",
    "fab_distribution_id", "assigned_process_layer", "fab_pair_key", "unit_width", "unit_height",
    "physical_length_unit", "parent_design_id", "change_description", "calc_version",
)
LAYOUT_COLUMNS = ("sample_id", "slot", "path", "pixel_size", "origin_x", "origin_y")
HOTSPOT_COLUMNS = (
    "sample_id", "particle_size", "path", "particle_size_unit", "pixel_size", "origin_x",
    "origin_y",
)
CA_RECORD_COLUMNS = ("sample_id", "particle_size", "particle_size_unit", "ca_value", "ca_unit")
LAYER_STACK_COLUMNS = (
    "product_id", "design_version", "layer_name", "stack_order", "fab_process_layer",
)
_SAMPLE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")
FAB_INPUT_WARNING = (
    "Fab 입력(fab.json 또는 fab_measurements.csv)은 생성하지 않았습니다. "
    "simcompare 실행 전 실제 Fab 자료를 data_dir에 추가해야 합니다."
)


def _write_csv(path: Path, columns: tuple[str, ...], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns))
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "" if row[key] is None else str(row[key]) for key in columns})


def _fab_pair_key(raw_unit: str, layer_a: str, layer_b: str,
                  registry: dict[str, tuple[str, str, str]]) -> str:
    raw = (raw_unit, layer_a, layer_b)
    if any(not isinstance(c, str) or not c or c != c.strip() for c in raw):
        raise InputContractError("fab_pair_key components must be nonblank without edged whitespace",
                                 field="fab_pair_key")
    key = f"{raw_unit}:{layer_a}-{layer_b}"
    previous = registry.get(key)
    if previous is not None and previous != raw:
        raise InputContractError(f"fab_pair_key collision {key!r}: {previous!r} and {raw!r}",
                                 field="fab_pair_key")
    registry[key] = raw
    return key


def export_simcompare(run_dir: Path, manifest_path: Path, out_dir: Path) -> Path:
    """Write simcompare bundles/layouts/hotspots/ca_records for every applicable supported (step, unit)."""
    run_dir = Path(run_dir)
    manifest_path = Path(manifest_path)
    out_dir = Path(out_dir)
    manifest = load_manifest(manifest_path)
    product_id = manifest.product_id
    if out_dir.exists() and (not out_dir.is_dir() or any(out_dir.iterdir())):
        raise InputContractError(
            f"출력 경로가 비어 있지 않습니다(덮어쓰기 금지): {out_dir}",
            product_id=product_id,
            field="output",
        )
    masks = load_product_masks(manifest, manifest_path)
    grid = masks.grid
    if not grid.is_isotropic:
        raise InputContractError(
            "simcompare는 단일 pixel_size만 지원하므로 x/y pixel 크기가 같아야 합니다",
            product_id=product_id,
            file=str(manifest_path),
            field="coordinate_system.pixel_size",
        )
    stack_z = build_stack_z(manifest)

    run_manifest_path = run_dir / "run_manifest.json"
    try:
        run_manifest = json.loads(run_manifest_path.read_text(encoding="utf-8"))
        ca = pd.read_parquet(run_dir / "ca.parquet")
        store = open_hotspot_store(run_dir / "hotspots.zarr")
    except (OSError, ValueError) as exc:
        raise InputContractError(
            f"실행 결과를 읽을 수 없습니다: {exc}",
            product_id=product_id,
            file=str(run_dir),
            field="run",
        ) from exc
    missing = check_run_manifest(run_manifest)
    if missing:
        raise InputContractError(
            f"run manifest에 필수 필드가 없습니다: {missing}",
            product_id=product_id,
            file=str(run_manifest_path),
            field=missing[0],
        )
    ca_id = "failure_unit_id" if "failure_unit_id" in ca.columns else "pair_id"
    missing_columns = [c for c in ("source_step", ca_id, "particle_size_um",
                                   "ca_nominal_normalized") if c not in ca.columns]
    if missing_columns:
        raise InputContractError(
            f"ca.parquet에 필수 열이 없습니다: {missing_columns}",
            product_id=product_id,
            file=str(run_dir / "ca.parquet"),
            field=missing_columns[0],
        )
    missing_attrs = [a for a in ("source_steps", "particle_sizes_um")
                     if a not in store.attrs]
    if missing_attrs or "nominal_failure" not in store:
        raise InputContractError(
            f"hotspots.zarr에 필수 항목이 없습니다: {missing_attrs or ['nominal_failure']}",
            product_id=product_id,
            file=str(run_dir / "hotspots.zarr"),
            field=(missing_attrs or ["nominal_failure"])[0],
        )
    if run_manifest.get("product_id") != product_id:
        raise InputContractError(
            f"실행 결과의 product_id({run_manifest.get('product_id')})가 manifest와 다릅니다",
            product_id=product_id,
            file=str(run_manifest_path),
            field="product_id",
        )

    step_ids = list(store.attrs["source_steps"])
    sizes = [float(s) for s in store.attrs["particle_sizes_um"]]

    def contract(message: str, field: str = "failure_unit_ids") -> None:
        raise InputContractError(message, product_id=product_id, field=field)

    raw_ids = store.attrs.get("failure_unit_ids", store.attrs.get("pair_ids"))
    if raw_ids is None:
        contract("hotspots.zarr에 failure_unit_ids 또는 pair_ids가 필요합니다")
    pair_ids = list(raw_ids)
    if "pair_ids" in store.attrs and list(store.attrs["pair_ids"]) != pair_ids:
        contract("failure_unit_ids와 pair_ids가 다릅니다")
    if len(set(pair_ids)) != len(pair_ids):
        contract("hotspots.zarr에 중복 unit ID가 있습니다")
    nominal = store["nominal_failure"]
    steps = {s.step_id: s for s in manifest.particle_source_steps}
    specs = {u.unit_id: u for u in unit_specs(manifest, masks)}
    electrodes = {e.electrode_id: e for e in manifest.electrodes}
    pixel_size = grid.pixel_size_x_um
    calc_version = f"particle-sim-{run_manifest['simulator_version']}"
    warnings: list[str] = [FAB_INPUT_WARNING]
    for attr, member in (("unit_kinds", "unit_kind"), ("pair_modes", "failure_mode")):
        values = store.attrs.get(attr)
        if values is not None and len(values) != len(pair_ids):
            contract(f"{attr} 길이가 unit 축과 다릅니다", attr)
        for index, unit_id in enumerate(pair_ids):
            if unit_id not in specs:
                contract(f"알 수 없는 failure unit ID: {unit_id!r}")
            if values is not None and values[index] != getattr(specs[unit_id], member):
                contract(f"{unit_id!r}: {attr}가 manifest와 다릅니다", attr)
    if nominal.shape != (len(step_ids), len(sizes), len(pair_ids), *grid.shape_yx):
        contract("nominal_failure shape가 축/마스크와 다릅니다", "nominal_failure")
    if "failure_unit_id" in ca and "pair_id" in ca and not ca.failure_unit_id.equals(ca.pair_id):
        contract("ca.parquet failure_unit_id와 pair_id가 다릅니다")

    # Validate the complete sample set before creating any output files.
    samples = []
    identities: dict[str, tuple[str, str, str, str]] = {}
    sample_id_mapping = []
    pair_keys: dict[str, tuple[str, str, str]] = {}
    array_paths: dict[str, tuple[str, str]] = {}
    layout_paths: dict[str, str] = {}
    for unit_id in pair_ids:
        spec = specs[unit_id]
        if spec.unit_kind == "contact" and spec.failure_mode == "contact_not_open":
            warnings.append(f"{unit_id}: kind={spec.unit_kind}, mode={spec.failure_mode}; "
                            "excluded: simcompare does not support contact_not_open")
        elif (spec.unit_kind, spec.failure_mode) not in {
            ("electrode_pair", "lateral_short"), ("electrode_pair", "vertical_short"),
            ("open_path", "lateral_open"), ("open_electrode", "lateral_open"),
        }:
            contract(f"지원하지 않는 unit metadata: {unit_id!r}")
    for step_index, step_id in enumerate(step_ids):
        if step_id not in steps:
            contract(f"알 수 없는 source step: {step_id!r}", "source_steps")
        assigned = steps[step_id].fab_distribution_layer
        if not isinstance(assigned, str) or not assigned.strip():
            contract(f"{step_id!r}: source Fab layer must be nonblank", "fab_distribution_layer")
        if "fab_distribution_layer" in ca:
            labels = ca.loc[ca.source_step == step_id, "fab_distribution_layer"]
            if labels.empty or not labels.eq(assigned).all():
                contract(f"{step_id!r}: ca.parquet fab_distribution_layer must be uniform "
                         "and equal the manifest source Fab layer", "fab_distribution_layer")
        for pair_index, pair_id in enumerate(pair_ids):
            pair = specs[pair_id]
            if pair.unit_kind == "contact":
                continue
            rows = ca[(ca["source_step"] == step_id) & (ca[ca_id] == pair_id)]
            for column, expected in (("unit_kind", pair.unit_kind),
                                     ("failure_mode", pair.failure_mode)):
                if column in rows and not rows[column].eq(expected).all():
                    contract(f"{pair_id!r}: ca.parquet {column}가 manifest와 다릅니다", column)
            values = []
            classifications = []
            for size in sizes:
                match = rows[np.isclose(rows["particle_size_um"].astype(float), size)]
                if len(match) != 1:
                    contract(f"ca.parquet에서 ({step_id}, {pair_id}, {size:g} um) 행이 "
                             f"{len(match)}개입니다", "ca_nominal_normalized")
                values.append(float(match["ca_nominal_normalized"].iloc[0]))
                if "classification" in match:
                    classifications.append(match["classification"].iloc[0])
            if len(set(classifications)) > 1:
                contract(f"({step_id}, {pair_id}): mixed classifications", "classification")
            skipped = classifications == [NOT_APPLICABLE_STEP] * len(sizes)
            applicable = pair.applicable_steps is None or step_id in pair.applicable_steps
            if classifications and skipped == applicable:
                contract(f"({step_id}, {pair_id}): applicability contradicts manifest",
                         "classification")
            if skipped:
                continue
            if not all(np.isfinite(values)):
                contract(f"({step_id}, {pair_id}): ca_value must be finite", "ca_nominal_normalized")
            layer_a = electrodes[pair.electrode_a].layer_name
            layer_b = ("N/A" if pair.failure_mode == "lateral_open" else
                       electrodes[pair.electrode_b].layer_name)
            key = _fab_pair_key(pair_id, layer_a, layer_b, pair_keys)
            raw_tuple = (product_id, manifest.design_version, step_id, pair_id)
            sample_id = re.sub(r"[^A-Za-z0-9_.-]", "_", "-".join(raw_tuple))
            previous = identities.get(sample_id)
            if previous is not None and previous != raw_tuple:
                contract(f"sample_id collision {sample_id!r}: {previous!r} and {raw_tuple!r}",
                         "sample_id")
            identities[sample_id] = raw_tuple
            slots = ("a",) if pair.failure_mode == "lateral_open" else ("a", "b")
            paths = []
            for slot in slots:
                electrode_id = pair.electrode_a if slot == "a" else pair.electrode_b
                if electrode_id not in layout_paths:
                    sanitized = re.sub(r"[^A-Za-z0-9_.-]", "_", electrode_id)
                    rel = f"arrays/layout-{sanitized}.npy"
                    layout_paths[electrode_id] = rel
                    paths.append(rel)
            paths.extend(f"arrays/{sample_id}-hot{size:g}.npy" for size in sizes)
            for path in paths:
                previous_path = array_paths.get(path.casefold())
                if previous_path is not None:
                    contract(f"array path collision {path!r}: sample {previous_path[0]!r} "
                             f"({previous_path[1]!r}) and sample {sample_id!r} ({path!r})",
                             "array_path")
                array_paths[path.casefold()] = (sample_id, path)
            if (pair.failure_mode == "lateral_open" and classifications and
                    classifications[0] in (BASELINE_OPEN, TERMINAL_LOST, INSUFFICIENT_TERMINALS)):
                warnings.append(
                    f"{pair_id!r}: step={step_id!r}, classification={classifications[0]}; "
                    "CA 0 reflects an excluded baseline condition, not low risk"
                )
            sample_id_mapping.append({"raw_tuple": list(raw_tuple), "sample_id": sample_id,
                                      "fab_pair_key": key, "assigned_process_layer": assigned})
            samples.append((step_index, step_id, pair_index, pair, sample_id, values))

    hotspot_writes = []
    bundles: list[dict] = []
    layouts: list[dict] = []
    hotspots: list[dict] = []
    ca_records: list[dict] = []
    geometry: dict[str, dict] = {}
    thickness = {name: layer.thickness_um for name, layer in stack_z.items()}

    for step_index, step_id, pair_index, pair, sample_id, values in samples:
        step = steps[step_id]
        ea = electrodes[pair.electrode_a]
        eb = None if pair.failure_mode == "lateral_open" else electrodes[pair.electrode_b]
        assigned = step.fab_distribution_layer
        bundles.append(
            {
                "sample_id": sample_id,
                "product_id": product_id,
                "product_family": manifest.product_family,
                "design_version": manifest.design_version,
                "design_lineage_id": manifest.design_lineage_id,
                "electrode_a_role": ea.electrode_role,
                "electrode_b_role": eb.electrode_role if eb else "N/A",
                "layer_a": ea.layer_name,
                "layer_b": eb.layer_name if eb else "N/A",
                "failure_type": pair.failure_mode,
                "fab_id": manifest.fab_id,
                "fab_distribution_id": None,
                "assigned_process_layer": assigned,
                "fab_pair_key": _fab_pair_key(pair.unit_id, ea.layer_name,
                                             eb.layer_name if eb else "N/A", pair_keys),
                "unit_width": manifest.coordinate_system.roi_width_um,
                "unit_height": manifest.coordinate_system.roi_height_um,
                "physical_length_unit": "um",
                "parent_design_id": None,
                "change_description": None,
                "calc_version": calc_version,
            }
        )
        slots = [("a", pair.electrode_a)]
        if eb is not None:
            slots.append(("b", pair.electrode_b))
        for slot, electrode_id in slots:
            rel = layout_paths[electrode_id]
            layouts.append(
                {
                    "sample_id": sample_id,
                    "slot": slot,
                    "path": rel,
                    "pixel_size": pixel_size,
                    "origin_x": 0.0,
                    "origin_y": 0.0,
                }
            )
        for size_index, size in enumerate(sizes):
            rel = f"arrays/{sample_id}-hot{size:g}.npy"
            hotspot_writes.append((rel, step_index, size_index, pair_index))
            hotspots.append(
                {
                    "sample_id": sample_id,
                    "particle_size": size,
                    "path": rel,
                    "particle_size_unit": "um",
                    "pixel_size": pixel_size,
                    "origin_x": 0.0,
                    "origin_y": 0.0,
                }
            )
            ca_records.append(
                {
                    "sample_id": sample_id,
                    "particle_size": size,
                    "particle_size_unit": "um",
                    "ca_value": values[size_index],
                    "ca_unit": "1",
                }
            )
        geometry[sample_id] = {
            "thickness": thickness,
            "shape_metrics": {},
            "length_unit": "um",
        }

    layer_stack = [
        {
            "product_id": product_id,
            "design_version": manifest.design_version,
            "layer_name": layer.layer_name,
            "stack_order": layer.process_order,
            "fab_process_layer": layer.layer_name,
        }
        for layer in stack_z.values()
    ]
    provenance = {
        "mapping_contract_version": "1.2",
        "calc_settings": {
            "run_id": run_manifest["run_id"],
            "config_hash": run_manifest["config_hash"],
            "assumption_set": run_manifest["assumption_set"],
            "legacy_kernel": run_manifest["legacy_kernel"],
            "post_step_policy": run_manifest["post_step_policy"],
            "hotspot_semantics": "nominal binary geometric contact risk",
        },
        "calc_version": calc_version,
        "is_synthetic": manifest.is_synthetic,
        "sample_id_mapping": sample_id_mapping,
        "array_storage": {
            "format": "npy",
            "layout_dtype": "uint8",
            "hotspot_dtype": "uint8 if binary, otherwise float32",
            "layout_deduplication": "one file per electrode_id across all samples and steps",
        },
    }
    json_outputs = {
        name: json.dumps(payload, indent=2, ensure_ascii=False) + "\n"
        for name, payload in (
            ("geometry.json", geometry),
            ("provenance.json", provenance),
            ("export_warnings.json", warnings),
        )
    }
    csv_outputs = (
        ("bundles.csv", BUNDLE_COLUMNS, bundles),
        ("layouts.csv", LAYOUT_COLUMNS, layouts),
        ("hotspots.csv", HOTSPOT_COLUMNS, hotspots),
        ("ca_records.csv", CA_RECORD_COLUMNS, ca_records),
        ("layer_stack.csv", LAYER_STACK_COLUMNS, layer_stack),
    )
    # Inspect each plane to size the actual binary uint8 / nonbinary float32 output.
    # Headers, CSV quoting and filesystem allocation retain conservative slack.
    hotspot_dtypes = []
    for _, si, di, ui in hotspot_writes:
        plane = np.asarray(nominal[si, di, ui], dtype=np.float32)
        hotspot_dtypes.append(np.dtype(np.uint8 if np.all((plane == 0) | (plane == 1))
                                       else np.float32))
    pixels = int(np.prod(grid.shape_yx))
    estimated_bytes = sum(pixels + 128 + 4096 for _ in layout_paths)
    estimated_bytes += sum(pixels * dtype.itemsize + 128 + 4096
                           for dtype in hotspot_dtypes)
    estimated_bytes += sum(len(value.encode("utf-8")) + 4096
                           for value in json_outputs.values())
    for _, columns, rows in csv_outputs:
        estimated_bytes += 4096 + sum(len(c.encode("utf-8")) + 3 for c in columns)
        estimated_bytes += sum(
            sum(2 * len(str(row[c]).encode("utf-8")) + 3 for c in columns)
            for row in rows
        )
    parent = out_dir.absolute().parent
    capacity_parent = parent
    while not capacity_parent.exists():
        capacity_parent = capacity_parent.parent
    free_bytes = shutil.disk_usage(capacity_parent).free
    if estimated_bytes > free_bytes:
        contract(f"Insufficient disk space: estimated {estimated_bytes} bytes, "
                 f"free {free_bytes} bytes", "output")

    parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{out_dir.name}-", dir=parent))
    try:
        (staging / "arrays").mkdir()
        # mkdir applies the process umask without temporarily changing global state.
        staging.chmod((staging / "arrays").stat().st_mode & 0o777)
        for electrode_id, rel in layout_paths.items():
            np.save(staging / rel, masks.electrodes[electrode_id].astype(np.uint8))
        for (rel, step_index, size_index, pair_index), dtype in zip(
                hotspot_writes, hotspot_dtypes, strict=True):
            hotspot = np.asarray(nominal[step_index, size_index, pair_index], dtype=np.float32)
            hotspot = hotspot.astype(dtype)
            np.save(staging / rel, hotspot)
        for name, columns, rows in csv_outputs:
            _write_csv(staging / name, columns, rows)
        for name, content in json_outputs.items():
            (staging / name).write_text(content, encoding="utf-8")
        # Recheck the overwrite guard in case another process populated the target.
        if out_dir.exists() and (not out_dir.is_dir() or any(out_dir.iterdir())):
            contract(f"출력 경로가 비어 있지 않습니다(덮어쓰기 금지): {out_dir}", "output")
        staging.rename(out_dir)
    except BaseException:
        shutil.rmtree(staging)
        raise
    return out_dir
