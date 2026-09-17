"""Export a nominal run into the simcompare example adapter data_dir format."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from particle_sim.config import load_manifest
from particle_sim.errors import InputContractError
from particle_sim.io.masks import load_product_masks
from particle_sim.io.outputs import open_hotspot_store
from particle_sim.io.stack import build_stack_z

BUNDLE_COLUMNS = (
    "sample_id", "product_id", "product_family", "design_version", "design_lineage_id",
    "electrode_a_role", "electrode_b_role", "layer_a", "layer_b", "failure_type", "fab_id",
    "fab_distribution_id", "assigned_process_layer", "unit_width", "unit_height",
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


def export_simcompare(run_dir: Path, manifest_path: Path, out_dir: Path) -> Path:
    """Write simcompare bundles/layouts/hotspots/ca_records for every (step, pair)."""
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
    if run_manifest.get("product_id") != product_id:
        raise InputContractError(
            f"실행 결과의 product_id({run_manifest.get('product_id')})가 manifest와 다릅니다",
            product_id=product_id,
            file=str(run_manifest_path),
            field="product_id",
        )

    step_ids = list(store.attrs["source_steps"])
    sizes = [float(s) for s in store.attrs["particle_sizes_um"]]
    pair_ids = list(store.attrs["pair_ids"])
    nominal = store["nominal_failure"]
    steps = {s.step_id: s for s in manifest.particle_source_steps}
    pairs = {p.pair_id: p for p in manifest.electrode_pairs}
    electrodes = {e.electrode_id: e for e in manifest.electrodes}
    pixel_size = grid.pixel_size_x_um
    calc_version = f"particle-sim-{run_manifest['simulator_version']}"

    arrays_dir = out_dir / "arrays"
    arrays_dir.mkdir(parents=True, exist_ok=True)
    bundles: list[dict] = []
    layouts: list[dict] = []
    hotspots: list[dict] = []
    ca_records: list[dict] = []
    geometry: dict[str, dict] = {}
    warnings: list[str] = [FAB_INPUT_WARNING]
    thickness = {name: layer.thickness_um for name, layer in stack_z.items()}

    for step_index, step_id in enumerate(step_ids):
        step = steps[step_id]
        for pair_index, pair_id in enumerate(pair_ids):
            pair = pairs[pair_id]
            sample_id = f"{product_id}-{manifest.design_version}-{step_id}-{pair_id}"
            if not _SAMPLE_ID.match(sample_id):
                raise InputContractError(
                    f"sample_id에 허용되지 않는 문자가 있습니다: {sample_id!r}",
                    product_id=product_id,
                    field="sample_id",
                )
            ea = electrodes[pair.electrode_a]
            eb = electrodes[pair.electrode_b]
            if pair.failure_mode == "vertical_short":
                assigned = ""
                warnings.append(
                    f"{sample_id}: simcompare는 vertical_short의 Fab layer를 상부 전극 layer로 "
                    f"강제하지만 simulator는 이물 발생 공정({step_id}) 기준이므로 "
                    "assigned_process_layer를 비워 두었습니다"
                )
            else:
                assigned = step.fab_distribution_layer
            bundles.append(
                {
                    "sample_id": sample_id,
                    "product_id": product_id,
                    "product_family": manifest.product_family,
                    "design_version": manifest.design_version,
                    "design_lineage_id": manifest.design_lineage_id,
                    "electrode_a_role": ea.electrode_role,
                    "electrode_b_role": eb.electrode_role,
                    "layer_a": ea.layer_name,
                    "layer_b": eb.layer_name,
                    "failure_type": pair.failure_mode,
                    "fab_id": manifest.fab_id,
                    "fab_distribution_id": None,
                    "assigned_process_layer": assigned,
                    "unit_width": manifest.coordinate_system.roi_width_um,
                    "unit_height": manifest.coordinate_system.roi_height_um,
                    "physical_length_unit": "um",
                    "parent_design_id": None,
                    "change_description": None,
                    "calc_version": calc_version,
                }
            )
            for slot, electrode_id in (("a", pair.electrode_a), ("b", pair.electrode_b)):
                rel = f"arrays/{sample_id}-{slot}.npy"
                np.save(out_dir / rel, masks.electrodes[electrode_id].astype(np.float32))
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
                hotspot = np.asarray(nominal[step_index, size_index, pair_index], dtype=np.float32)
                np.save(out_dir / rel, hotspot)
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
                match = ca[
                    (ca["source_step"] == step_id)
                    & (ca["pair_id"] == pair_id)
                    & np.isclose(ca["particle_size_um"].astype(float), size)
                ]
                if len(match) != 1:
                    raise InputContractError(
                        f"ca.parquet에서 ({step_id}, {pair_id}, {size:g} um) 행이 "
                        f"{len(match)}개입니다",
                        product_id=product_id,
                        file=str(run_dir / "ca.parquet"),
                        field="ca_nominal_normalized",
                    )
                ca_records.append(
                    {
                        "sample_id": sample_id,
                        "particle_size": size,
                        "particle_size_unit": "um",
                        "ca_value": float(match["ca_nominal_normalized"].iloc[0]),
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
    _write_csv(out_dir / "bundles.csv", BUNDLE_COLUMNS, bundles)
    _write_csv(out_dir / "layouts.csv", LAYOUT_COLUMNS, layouts)
    _write_csv(out_dir / "hotspots.csv", HOTSPOT_COLUMNS, hotspots)
    _write_csv(out_dir / "ca_records.csv", CA_RECORD_COLUMNS, ca_records)
    _write_csv(out_dir / "layer_stack.csv", LAYER_STACK_COLUMNS, layer_stack)
    provenance = {
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
    }
    for name, payload in (
        ("geometry.json", geometry),
        ("provenance.json", provenance),
        ("export_warnings.json", warnings),
    ):
        (out_dir / name).write_text(
            json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
    return out_dir
