"""Derived nominal CNO/open accounting from independently evaluated raw maps.

The adapter's explicit circuit relation is the only authority for removing an
open event. Source-run maps and CA tables are opened read-only and never changed.
"""

from __future__ import annotations

import json
import math
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import zarr

from particle_sim.config import canonical_hash, file_sha256, load_manifest
from particle_sim.errors import InputContractError
from particle_sim.io.masks import load_product_masks
from particle_sim.io.outputs import open_hotspot_store, write_json
from particle_sim.units import KIND_CONTACT, KIND_OPEN_ELECTRODE, KIND_OPEN_PATH, unit_specs
from particle_sim.validation.central_window import (
    CENTRE_PIXELS,
    CENTRE_START,
    CENTRE_STOP,
    ROI_SHAPE,
)


def split_open(raw: np.ndarray, related_cno: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """Classify one particle condition, preserving the independent raw result."""
    cno = np.zeros(raw.shape, dtype=bool)
    for item in related_cno:
        if item.shape != raw.shape:
            raise ValueError("CNO and open map shapes disagree")
        cno |= item
    removed = raw & cno
    final = raw & ~cno
    return final, removed


def _binary_map(array) -> np.ndarray:
    values = np.asarray(array)
    if values.shape != ROI_SHAPE or np.any((values != 0) & (values != 1)):
        raise InputContractError("source map must be a binary 1200x1200 map", field="run")
    return values.astype(bool)


def _relations(manifest, masks, sidecar: dict, specs) -> dict[str, tuple[str, ...]]:
    if not isinstance(sidecar, dict):
        raise InputContractError("adapter sidecar must be a mapping", field="circuit")
    circuit = sidecar.get("circuit")
    if not isinstance(circuit, dict) or circuit.get("net_status") != "declared":
        raise InputContractError("dedup requires declared circuit identities", field="circuit")
    if circuit.get("auto_open") != "disabled" or any(
        spec.unit_kind == KIND_OPEN_ELECTRODE for spec in specs
    ):
        raise InputContractError("auto_open units are unsupported by dedup", field="auto_open")
    ids = [spec.unit_id for spec in specs]
    if len(ids) != len(set(ids)):
        raise InputContractError("duplicate failure unit ID", field="failure_units")
    by_id = {spec.unit_id: spec for spec in specs}
    for path in manifest.open_paths:
        if path.fallback != "error" or not (
            (path.terminal_a_yx_px is not None and path.terminal_b_yx_px is not None)
            or (path.terminal_a_mask_path is not None and path.terminal_b_mask_path is not None)
        ):
            raise InputContractError(
                "explicit open terminals and fallback=error required", field="open_paths"
            )
    details = circuit.get("open_path_circuit_metadata")
    rows = circuit.get("cno_open_relations")
    if not isinstance(details, list) or not isinstance(rows, list):
        raise InputContractError("missing explicit path identities or relations", field="circuit")
    path_details = {}
    for detail in details:
        if not isinstance(detail, dict) or not all(
            isinstance(detail.get(key), str) and detail[key].strip()
            for key in ("path_id", "circuit_connection_id", "circuit_net_id")
        ):
            raise InputContractError("missing circuit connection identity", field="circuit")
        path_id = detail["path_id"]
        if (
            path_id in path_details
            or path_id not in by_id
            or by_id[path_id].unit_kind != KIND_OPEN_PATH
        ):
            raise InputContractError("duplicate or unknown open path ID", field="circuit")
        path_details[path_id] = detail
    if set(path_details) != {path.path_id for path in manifest.open_paths}:
        raise InputContractError("open path circuit identities incomplete", field="circuit")
    electrodes = {item.electrode_id: item for item in manifest.electrodes}
    for path in manifest.open_paths:
        if electrodes[path.electrode].net_id != path_details[path.path_id]["circuit_net_id"]:
            raise InputContractError("open path circuit net disagrees", field="circuit")
    links = {} if masks.contact_topology is None else masks.contact_topology.links
    for link_id in links:
        if link_id not in by_id or by_id[link_id].unit_kind != KIND_CONTACT:
            raise InputContractError("unknown CNO link unit", field="circuit")
    related: dict[str, list[str]] = {path.path_id: [] for path in manifest.open_paths}
    seen = set()
    for row in rows:
        if not isinstance(row, dict) or not all(
            isinstance(row.get(key), str) and row[key].strip()
            for key in ("cno_link_id", "open_path_id", "circuit_connection_id")
        ):
            raise InputContractError(
                "missing relation circuit connection identity", field="circuit"
            )
        link_id, path_id = row["cno_link_id"], row["open_path_id"]
        if link_id not in links or path_id not in path_details:
            raise InputContractError("unknown CNO link or open path ID", field="circuit")
        if (link_id, path_id) in seen:
            raise InputContractError("duplicate CNO/open relation", field="circuit")
        seen.add((link_id, path_id))
        if row["circuit_connection_id"] != path_details[path_id]["circuit_connection_id"]:
            raise InputContractError(
                "relation circuit connection identity disagrees", field="circuit"
            )
        link = links[link_id]
        net = path_details[path_id]["circuit_net_id"]
        if (
            electrodes[link.lower_electrode].net_id != net
            or electrodes[link.upper_electrode].net_id != net
        ):
            raise InputContractError("relation link circuit net disagrees", field="circuit")
        link_steps = by_id[link_id].applicable_steps
        path_steps = by_id[path_id].applicable_steps
        all_steps = {step.step_id for step in manifest.particle_source_steps}
        if not (set(link_steps or all_steps) & set(path_steps or all_steps)):
            raise InputContractError("relation has no common source step", field="circuit")
        related[path_id].append(link_id)
    return {path_id: tuple(link_ids) for path_id, link_ids in related.items()}


def analyze_cno_open(
    manifest_path: Path, sidecar_path: Path, run_dir: Path, output_dir: Path
) -> dict:
    """Write derived nominal maps and full/central statistics outside the raw run."""
    manifest_path, sidecar_path = Path(manifest_path), Path(sidecar_path)
    run_dir, output_dir = Path(run_dir), Path(output_dir)
    manifest = load_manifest(manifest_path)
    masks = load_product_masks(manifest, manifest_path)
    if (
        masks.grid.shape_yx != ROI_SHAPE
        or not math.isclose(masks.grid.pixel_size_x_um, 0.1, abs_tol=1e-9)
        or not math.isclose(masks.grid.pixel_size_y_um, 0.1, abs_tol=1e-9)
    ):
        raise InputContractError("dedup requires 1200x1200 at 0.1 um/px", field="coordinate_system")
    specs = unit_specs(manifest, masks)
    try:
        sidecar = json.loads(sidecar_path.read_text(encoding="utf-8"))
        run_manifest = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise InputContractError(f"dedup input cannot be read: {exc}", field="dedup") from exc
    if not isinstance(sidecar, dict):
        raise InputContractError("adapter sidecar must be a mapping", field="circuit")
    if sidecar.get("product_yaml_sha256") != file_sha256(manifest_path):
        raise InputContractError(
            "adapter sidecar product.yaml SHA-256 is missing or mismatched", field="sidecar"
        )
    related = _relations(manifest, masks, sidecar, specs)
    if (
        output_dir.resolve().is_relative_to(run_dir.resolve())
        or output_dir.resolve() == run_dir.resolve()
    ):
        raise InputContractError("dedup output must be outside the raw run", field="output")
    if output_dir.exists() and (not output_dir.is_dir() or any(output_dir.iterdir())):
        raise InputContractError("dedup output must be empty", field="output")
    if (
        run_manifest.get("product_id") != manifest.product_id
        or run_manifest.get("config_hash") != canonical_hash(manifest.model_dump(mode="json"))
        or run_manifest.get("random_seed") is not None
    ):
        raise InputContractError("nominal source run does not match product config", field="run")
    if not (run_dir / "hotspots.zarr").exists():
        raise InputContractError("independent nominal raw maps required", field="run")
    group = open_hotspot_store(run_dir / "hotspots.zarr")
    steps = [step.step_id for step in manifest.particle_source_steps]
    sizes = list(manifest.particle_sizes_um)
    ids = [spec.unit_id for spec in specs]
    if (
        list(group.attrs["source_steps"]) != steps
        or list(group.attrs["particle_sizes_um"]) != sizes
        or list(group.attrs["failure_unit_ids"]) != ids
        or list(group.attrs["unit_kinds"]) != [spec.unit_kind for spec in specs]
        or tuple(group["nominal_failure"].shape[-2:]) != ROI_SHAPE
    ):
        raise InputContractError("raw map axes disagree with manifest", field="run")
    if not (run_dir / "ca.parquet").exists() or not (run_dir / "ca_union.parquet").exists():
        raise InputContractError("raw CA tables required", field="run")
    open_ids = [path.path_id for path in manifest.open_paths]
    parent_existed = output_dir.parent.exists()
    stage_dir = None
    try:
        output_dir.parent.mkdir(parents=True, exist_ok=True)
        stage_dir = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
        result = _write_cno_open_output(
            manifest,
            sidecar_path,
            run_dir,
            stage_dir,
            group,
            steps,
            sizes,
            specs,
            ids,
            open_ids,
            related,
            sidecar,
            masks,
        )
        if output_dir.exists():
            output_dir.rmdir()
        stage_dir.rename(output_dir)
        return result
    except BaseException:
        if stage_dir is not None:
            shutil.rmtree(stage_dir, ignore_errors=True)
        if not parent_existed:
            try:
                output_dir.parent.rmdir()
            except OSError:
                pass
        raise


def _write_cno_open_output(
    manifest,
    sidecar_path: Path,
    run_dir: Path,
    output_dir: Path,
    group,
    steps: list[str],
    sizes: list[float],
    specs,
    ids: list[str],
    open_ids: list[str],
    related: dict[str, tuple[str, ...]],
    sidecar: dict,
    masks,
) -> dict:
    out = zarr.open_group(output_dir / "dedup_open.zarr", mode="w-")
    shape = (len(steps), len(sizes), len(open_ids), *ROI_SHAPE)
    chunks = (1, 1, 1, 512, 512)
    for key in ("open_final", "removed"):
        out.create_array(key, shape=shape, chunks=chunks, dtype="uint8", fill_value=0)
    out.attrs.update(
        {
            "source_steps": steps,
            "particle_sizes_um": sizes,
            "open_path_ids": open_ids,
            "particle_shape": "disc",
            "source_run": str(run_dir),
        }
    )
    id_index = {unit_id: index for index, unit_id in enumerate(ids)}
    connection_by_path = {
        item["path_id"]: item["circuit_connection_id"]
        for item in sidecar["circuit"]["open_path_circuit_metadata"]
    }
    any_index = list(group.attrs["union_keys"]).index("any")
    area = masks.grid.pixel_area_um2
    rows = []
    for si, step in enumerate(steps):
        for di, diameter in enumerate(sizes):
            raw_any = _binary_map(group["union_failure"][si, di, any_index])
            cno_maps = {}
            for path_index, path_id in enumerate(open_ids):
                raw = _binary_map(group["nominal_failure"][si, di, id_index[path_id]])
                if np.any(raw & ~raw_any):
                    raise InputContractError("raw any union omits an open pixel", field="run")
                applicable = set(specs[id_index[path_id]].applicable_steps or steps)
                links = [
                    link_id
                    for link_id in related[path_id]
                    if step in set(specs[id_index[link_id]].applicable_steps or steps)
                ]
                if step not in applicable:
                    links = []
                for link_id in links:
                    if link_id not in cno_maps:
                        cno_maps[link_id] = _binary_map(
                            group["nominal_failure"][si, di, id_index[link_id]]
                        )
                    if np.any(cno_maps[link_id] & ~raw_any):
                        raise InputContractError("raw any union omits a CNO pixel", field="run")
                final, removed = split_open(raw, [cno_maps[link_id] for link_id in links])
                if not np.array_equal(final | removed, raw):
                    raise AssertionError("open partition failed")
                out["open_final"][si, di, path_index] = final.astype(np.uint8)
                out["removed"][si, di, path_index] = removed.astype(np.uint8)
                centre = np.s_[CENTRE_START:CENTRE_STOP, CENTRE_START:CENTRE_STOP]
                rows.append(
                    {
                        "source_step": step,
                        "particle_size_um": diameter,
                        "particle_shape": "disc",
                        "sample_id": None,
                        "open_path_id": path_id,
                        "circuit_connection_id": connection_by_path[path_id],
                        "related_cno_link_ids": list(links),
                        "raw_pixels": int(raw.sum()),
                        "removed_pixels": int(removed.sum()),
                        "final_pixels": int(final.sum()),
                        "raw_ca_um2": float(raw.sum() * area),
                        "removed_ca_um2": float(removed.sum() * area),
                        "final_ca_um2": float(final.sum() * area),
                        "central_denominator_pixels": CENTRE_PIXELS,
                        "central_raw_pixels": int(raw[centre].sum()),
                        "central_removed_pixels": int(removed[centre].sum()),
                        "central_final_pixels": int(final[centre].sum()),
                        "central_raw_ca_um2": float(raw[centre].sum() * area),
                        "central_removed_ca_um2": float(removed[centre].sum() * area),
                        "central_final_ca_um2": float(final[centre].sum() * area),
                        "central_raw_normalized": float(raw[centre].sum() / CENTRE_PIXELS),
                        "central_removed_normalized": float(removed[centre].sum() / CENTRE_PIXELS),
                        "central_final_normalized": float(final[centre].sum() / CENTRE_PIXELS),
                    }
                )
    pd.DataFrame(rows).to_parquet(output_dir / "cno_open_dedup.parquet", index=False)
    result = {
        "contract": "derived geometric CNO/open accounting; no physical process correctness claim",
        "source_run": str(run_dir),
        "sidecar": str(sidecar_path),
        "source_run_manifest_sha256": file_sha256(run_dir / "run_manifest.json"),
        "sidecar_sha256": file_sha256(sidecar_path),
        "relations": sidecar["circuit"]["cno_open_relations"],
        "raw_ca_sha256": file_sha256(run_dir / "ca.parquet"),
        "raw_ca_union_sha256": file_sha256(run_dir / "ca_union.parquet"),
        "central_window_yx_half_open": [CENTRE_START, CENTRE_STOP],
        "central_denominator_pixels": CENTRE_PIXELS,
        "any_union": "pixel_identical_by_CNO_subset_and_unchanged_raw_units",
        "unrelated_open_and_short_maps": "unchanged_in_read_only_source_run",
        "rows": len(rows),
    }
    write_json(result, output_dir / "dedup_manifest.json")
    return result
