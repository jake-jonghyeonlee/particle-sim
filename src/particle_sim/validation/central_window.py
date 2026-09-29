"""E2 central-window statistics and conservative local-support audit.

This is a derived analysis of existing full-ROI maps. A local reach check cannot
certify the nonlocal connectivity used to decide circuit opens.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from particle_sim.config import canonical_hash, file_sha256, load_manifest
from particle_sim.errors import InputContractError
from particle_sim.geometry.legacy import legacy_iterations
from particle_sim.io.masks import MaskSet, load_product_masks
from particle_sim.io.outputs import open_hotspot_store, write_json
from particle_sim.schema import PARTICLE_SIZES_UM, ProductManifest
from particle_sim.units import KIND_OPEN_ELECTRODE, KIND_OPEN_PATH, unit_specs

ROI_SHAPE = (1200, 1200)
CENTRE_START = 300
CENTRE_STOP = 900
HALO_PX = 300
CENTRE_PIXELS = 360_000
PIXEL_UM = 0.1
NONLOCAL_OPEN = "nonlocal_open_connectivity_not_certified"


def _check_geometry(manifest: ProductManifest, masks: MaskSet) -> None:
    grid = masks.grid
    if grid.shape_yx != ROI_SHAPE:
        raise InputContractError("E2 requires a 1200x1200 ROI", field="coordinate_system")
    if not (math.isclose(grid.pixel_size_x_um, PIXEL_UM, abs_tol=1e-9)
            and math.isclose(grid.pixel_size_y_um, PIXEL_UM, abs_tol=1e-9)):
        raise InputContractError("E2 requires 0.1 um/px in both axes", field="coordinate_system")
    if tuple(manifest.particle_sizes_um) != PARTICLE_SIZES_UM:
        raise InputContractError(
            "E2 requires exactly 1..30,40,50 um (32 bins)", field="particle_sizes_um"
        )
    mask_bearing_layers = {layer.layer_name for layer in manifest.layer_masks} | {
        contact.hole_layer for contact in manifest.contacts
    }
    for step in manifest.particle_source_steps:
        if step.lands_after_layer not in mask_bearing_layers:
            raise InputContractError(
                "maskless dielectric cannot have a particle source stage",
                field=f"particle_source_steps.{step.step_id}",
            )
    if (CENTRE_STOP - CENTRE_START) ** 2 != CENTRE_PIXELS:
        raise AssertionError("central window area changed")


def central_area_um2(values: np.ndarray, pixel_area_um2: float) -> float:
    """Integrate only [300,900)^2 from a full-size map."""
    if values.shape != ROI_SHAPE:
        raise ValueError(f"expected {ROI_SHAPE} full map, got {values.shape}")
    centre = np.asarray(values[CENTRE_START:CENTRE_STOP, CENTRE_START:CENTRE_STOP])
    if not np.isfinite(centre).all() or np.any((centre < 0) | (centre > 1)):
        raise ValueError("central hotspot values must be finite in [0,1]")
    return float(centre.sum(dtype=np.float64) * pixel_area_um2)


def _hole_half_span_px(masks: MaskSet, step_id: str, manifest: ProductManifest) -> float:
    topology = masks.contact_topology
    if topology is None:
        return 0.0
    applicable = {
        contact.contact_id
        for contact in manifest.contacts
        if step_id in contact.applicable_particle_source_steps
    }
    spans = [
        max((int(hole.cut_ys.max()) - int(hole.cut_ys.min())) / 2,
            (int(hole.cut_xs.max()) - int(hole.cut_xs.min())) / 2)
        for hole in topology.holes if hole.contact_id in applicable
    ]
    return max(spans, default=0.0)


def support_rows(manifest: ProductManifest, masks: MaskSet) -> list[dict]:
    """Audit each size and source stage without evaluating any pair hotspot."""
    _check_geometry(manifest, masks)
    grid = masks.grid
    specs = unit_specs(manifest, masks)
    rows: list[dict] = []
    for step in manifest.particle_source_steps:
        step_id = step.step_id
        hole_half_span = _hole_half_span_px(masks, step_id, manifest)
        open_units = [
            spec.unit_id for spec in specs
            if spec.unit_kind in (KIND_OPEN_PATH, KIND_OPEN_ELECTRODE)
            and (spec.applicable_steps is None or step_id in spec.applicable_steps)
        ]
        flat_film = any(
            step_id in contact.applicable_particle_source_steps
            and step.lands_after_layer == contact.hole_layer
            for contact in manifest.contacts
        )
        for diameter in PARTICLE_SIZES_UM:
            radius = diameter / (2 * PIXEL_UM)
            explicit_margin = max(0.0, manifest.contact_margin_um / PIXEL_UM)
            footprint_reach = radius
            contact_reach = radius + explicit_margin
            short_reach = contact_reach  # both contacts are queried at the same centre
            kernel_reach = legacy_iterations(diameter, grid, manifest.legacy_kernel) * (
                manifest.legacy_kernel.size // 2
            )
            hole_reach = radius + hole_half_span if hole_half_span else radius
            reaches = {
                "footprint": footprint_reach,
                "contact_margin": contact_reach,
                "legacy_short_kernel": kernel_reach,
                "short_geometry": short_reach,
                "full_cut_hole": hole_reach,
            }
            reasons = [
                f"{name}_reach_exceeds_300px"
                for name, reach in reaches.items() if reach > HALO_PX + 1e-9
            ]
            if open_units:
                reasons.append(NONLOCAL_OPEN)
            rows.append({
                "source_step": step_id,
                "particle_size_um": diameter,
                "supported": not reasons,
                "local_geometry_supported": not any(
                    reach > HALO_PX + 1e-9 for reach in reaches.values()
                ),
                "reasons": reasons,
                "radius_px": radius,
                "explicit_contact_margin_px": explicit_margin,
                "full_cut_hole_half_span_px": hole_half_span,
                "reaches_px": reaches,
                "open_unit_ids": open_units,
                "source_stage_semantics": (
                    "flat_film_approximation_not_physical_CNT_open_before_fill"
                    if flat_film else "configured_surface_model"
                ),
            })
    return rows


def analyze_central(
    manifest_path: Path, output_dir: Path, *, run_dir: Path | None = None
) -> dict:
    """Write support audit and optional central CA sidecars; never alter the source run."""
    manifest_path, output_dir = Path(manifest_path), Path(output_dir)
    manifest = load_manifest(manifest_path)
    masks = load_product_masks(manifest, manifest_path)
    rows = support_rows(manifest, masks)
    if output_dir.exists() and (not output_dir.is_dir() or any(output_dir.iterdir())):
        raise InputContractError("central analysis output must be empty", field="output")

    unit_rows: list[dict] = []
    union_rows: list[dict] = []
    source_info = None
    if run_dir is not None:
        run_dir = Path(run_dir)
        if output_dir.resolve().is_relative_to(run_dir.resolve()):
            raise InputContractError("central output must be outside the source run", field="output")
        source_info = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
        nominal = (run_dir / "hotspots.zarr").exists()
        if nominal:
            same_config = source_info.get("config_hash") == canonical_hash(
                manifest.model_dump(mode="json")
            )
        else:
            # A probability config_hash also covers the scenario, which the run does not keep;
            # match the product through the recorded manifest and mask file hashes instead.
            source_hashes = source_info.get("input_hashes") or {}
            same_config = (
                source_info.get("design_version") == manifest.design_version
                and source_hashes.get("manifest") == file_sha256(manifest_path)
                and all(
                    source_hashes.get(name) == digest
                    for name, digest in masks.input_hashes.items()
                )
            )
        if source_info.get("product_id") != manifest.product_id or not same_config:
            raise InputContractError("source run does not match product config", field="run")
        if not nominal:
            # The manifest's nominal bins do not bound scenario diameter/CD/overlay changes.
            for row in rows:
                row["supported"] = False
                row["reasons"].append("scenario_variation_influence_not_bounded")
        store_path = run_dir / ("hotspots.zarr" if nominal else "probability.zarr")
        group = open_hotspot_store(store_path)
        steps = [step.step_id for step in manifest.particle_source_steps]
        if (list(group.attrs["source_steps"]) != steps or
                tuple(group.attrs["particle_sizes_um"]) != PARTICLE_SIZES_UM or
                tuple(group["nominal_failure" if nominal else "probability"].shape[-2:])
                != ROI_SHAPE):
            raise InputContractError("source map axes do not match E2 contract", field="run")
        unit_ids = list(group.attrs["failure_unit_ids"])
        unit_modes = list(group.attrs["pair_modes" if nominal else "unit_modes"])
        unit_kinds = list(group.attrs["unit_kinds"])
        specs = unit_specs(manifest, masks)
        if (unit_ids != [spec.unit_id for spec in specs] or
                unit_modes != [spec.failure_mode for spec in specs] or
                unit_kinds != [spec.unit_kind for spec in specs]):
            raise InputContractError("source map failure units do not match config", field="run")
        spec_by_id = {spec.unit_id: spec for spec in specs}
        union_keys = list(group.attrs["union_keys"])
        map_key = "nominal_failure" if nominal else "probability"
        union_key = "union_failure" if nominal else "union_probability"
        area = masks.grid.pixel_area_um2
        for step_index, step_id in enumerate(steps):
            for size_index, diameter in enumerate(PARTICLE_SIZES_UM):
                support = rows[step_index * len(PARTICLE_SIZES_UM) + size_index]
                for unit_index, (unit_id, mode, kind) in enumerate(
                    zip(unit_ids, unit_modes, unit_kinds, strict=True)
                ):
                    # Read only the 600x600 central slice; the full-map Zarr stays untouched.
                    centre = np.asarray(group[map_key][
                        step_index, size_index, unit_index,
                        CENTRE_START:CENTRE_STOP, CENTRE_START:CENTRE_STOP,
                    ])
                    if not np.isfinite(centre).all() or np.any((centre < 0) | (centre > 1)):
                        raise InputContractError("invalid central map values", field="run")
                    spec = spec_by_id.get(unit_id)
                    is_open = kind in (KIND_OPEN_PATH, KIND_OPEN_ELECTRODE) and (
                        spec is None or spec.applicable_steps is None
                        or step_id in spec.applicable_steps
                    )
                    reasons = [reason for reason in support["reasons"] if reason != NONLOCAL_OPEN]
                    if is_open:
                        reasons.append(NONLOCAL_OPEN)
                    unit_rows.append({
                        "source_step": step_id, "particle_size_um": diameter,
                        "failure_unit_id": unit_id, "failure_mode": mode, "unit_kind": kind,
                        "central_ca_um2": float(centre.sum(dtype=np.float64) * area),
                        "central_normalized": float(centre.mean(dtype=np.float64)),
                        "supported": not reasons,
                        "unsupported_reasons": ";".join(reasons),
                    })
                for union_index, mode in enumerate(union_keys):
                    centre = np.asarray(group[union_key][
                        step_index, size_index, union_index,
                        CENTRE_START:CENTRE_STOP, CENTRE_START:CENTRE_STOP,
                    ])
                    if not np.isfinite(centre).all() or np.any((centre < 0) | (centre > 1)):
                        raise InputContractError("invalid central union values", field="run")
                    union_reasons = [
                        reason for reason in support["reasons"]
                        if reason != NONLOCAL_OPEN or mode in ("lateral_open", "any")
                    ]
                    union_rows.append({
                        "source_step": step_id, "particle_size_um": diameter,
                        "failure_mode": mode,
                        "central_ca_um2": float(centre.sum(dtype=np.float64) * area),
                        "central_normalized": float(centre.mean(dtype=np.float64)),
                        "supported": not union_reasons,
                        "unsupported_reasons": ";".join(union_reasons),
                    })

    output_dir.mkdir(parents=True, exist_ok=True)
    audit = {
        "contract": "E2 central geometric statistic, not physical simulation",
        "source_run": str(run_dir) if run_dir is not None else None,
        "source_run_manifest_sha256": (
            file_sha256(run_dir / "run_manifest.json") if run_dir is not None else None
        ),
        "roi_shape_yx": list(ROI_SHAPE),
        "pixel_size_um": PIXEL_UM,
        "central_window_yx_half_open": {
            "y": [CENTRE_START, CENTRE_STOP],
            "x": [CENTRE_START, CENTRE_STOP],
        },
        "central_pixels": CENTRE_PIXELS,
        "central_area_um2": CENTRE_PIXELS * masks.grid.pixel_area_um2,
        "halo_px": HALO_PX,
        "particle_sizes_um": list(PARTICLE_SIZES_UM),
        "source_steps": [step.step_id for step in manifest.particle_source_steps],
        "support": rows,
        "limitation": "Local reach does not certify nonlocal circuit-open connectivity.",
        "raw_full_map_outputs": "preserved in the source run; central tables are derived sidecars",
    }
    write_json(audit, output_dir / "central_support.json")
    if run_dir is not None:
        pd.DataFrame(unit_rows).to_parquet(output_dir / "central_ca.parquet", index=False)
        pd.DataFrame(union_rows).to_parquet(output_dir / "central_ca_union.parquet", index=False)
    return {"support_rows": len(rows), "unit_rows": len(unit_rows),
            "union_rows": len(union_rows), "unsupported": sum(not r["supported"] for r in rows)}
