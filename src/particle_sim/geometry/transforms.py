"""Apply one process-variation scenario to the masks and the layer stack."""

from __future__ import annotations

import numpy as np
from scipy.ndimage import distance_transform_edt, map_coordinates

from particle_sim.geometry.contacts import build_contact_topology
from particle_sim.io.masks import GridSpec, MaskSet
from particle_sim.io.stack import LayerZ
from particle_sim.schema import ProductManifest
from particle_sim.uncertainty.sampler import Scenario


def _half_step_um(mask: np.ndarray, grid: GridSpec) -> np.ndarray:
    """Half a pixel measured along the direction to the nearest opposite-phase pixel.

    On an anisotropic grid a vertical edge sits half of `pixel_size_x_um` away while a
    horizontal edge sits half of `pixel_size_y_um` away, so a single scalar cannot serve both.
    """
    sampling = (grid.pixel_size_y_um, grid.pixel_size_x_um)
    rows, cols = np.mgrid[0 : mask.shape[0], 0 : mask.shape[1]]
    half = np.zeros(mask.shape, dtype=np.float64)
    for phase in (True, False):
        target = ~mask if phase else mask
        if not target.any():
            continue
        _, indices = distance_transform_edt(~target, sampling=sampling, return_indices=True)
        dy_px = np.abs(indices[0] - rows).astype(np.float64)
        dx_px = np.abs(indices[1] - cols).astype(np.float64)
        length = np.hypot(dy_px, dx_px)
        with np.errstate(invalid="ignore", divide="ignore"):
            unit_y = np.where(length > 0, dy_px / length, 0.0)
            unit_x = np.where(length > 0, dx_px / length, 0.0)
        step = 0.5 * (unit_y * grid.pixel_size_y_um + unit_x * grid.pixel_size_x_um)
        selection = mask if phase else ~mask
        half[selection] = step[selection]
    return half


def signed_distance_um(mask: np.ndarray, grid: GridSpec) -> np.ndarray:
    """Distance from each pixel centre to the mask edge: negative inside, positive outside.

    `distance_transform_edt` measures the distance to the nearest pixel centre of the other
    phase, while the drawn edge sits half a pixel before it, so half a step along that
    direction is removed. `signed_distance_um(mask) <= 0` reproduces the original mask.
    """
    sampling = (grid.pixel_size_y_um, grid.pixel_size_x_um)
    if not mask.any():
        return np.full(mask.shape, np.inf, dtype=np.float64)
    if mask.all():
        return np.full(mask.shape, -np.inf, dtype=np.float64)
    half = _half_step_um(mask, grid)
    outside = distance_transform_edt(~mask, sampling=sampling) - half
    inside = distance_transform_edt(mask, sampling=sampling) - half
    return np.where(mask, -inside, outside)


def transform_mask(
    mask: np.ndarray,
    grid: GridSpec,
    *,
    overlay_xy_um: tuple[float, float] = (0.0, 0.0),
    cd_bias_um: float = 0.0,
) -> np.ndarray:
    """Sub-pixel overlay shift and CD bias applied through the signed distance field.

    Repeated nearest-neighbour shifts would stair-step the edges, so the distance field is
    interpolated instead and re-thresholded (spec 7.3). A positive cd_bias grows the mask.
    """
    dx_um, dy_um = overlay_xy_um
    if dx_um == 0.0 and dy_um == 0.0 and cd_bias_um == 0.0:
        return mask.copy()
    if not mask.any() or mask.all():
        return mask.copy()

    field = signed_distance_um(mask, grid)
    dy_px = dy_um / grid.pixel_size_y_um
    dx_px = dx_um / grid.pixel_size_x_um
    rows, cols = np.mgrid[0 : mask.shape[0], 0 : mask.shape[1]].astype(np.float64)
    # sampling the source position of each target pixel shifts the pattern by (+dy, +dx)
    shifted = map_coordinates(
        field,
        np.stack([rows - dy_px, cols - dx_px]),
        order=1,
        mode="nearest",
    )
    # signed distance <= 0 is the original mask; a positive bias grows it outwards
    return shifted <= cd_bias_um


def apply_scenario_to_masks(
    manifest: ProductManifest, masks: MaskSet, scenario: Scenario
) -> MaskSet:
    """Overlay, CD bias and the particle centre offset applied to every mask."""
    offset_x, offset_y = scenario.center_offset_xy_um
    layer_of = {e.electrode_id: e.layer_name for e in manifest.electrodes}
    grid = masks.grid

    def variation(layer_name: str) -> tuple[tuple[float, float], float]:
        overlay = scenario.overlay_xy_um.get(layer_name, (0.0, 0.0))
        # a particle centre offset is equivalent to shifting every mask the other way
        shifted = (overlay[0] - offset_x, overlay[1] - offset_y)
        return shifted, scenario.cd_bias_um.get(layer_name, 0.0)

    layers = {}
    for layer_name, mask in masks.layers.items():
        overlay, cd_bias = variation(layer_name)
        layers[layer_name] = transform_mask(
            mask, grid, overlay_xy_um=overlay, cd_bias_um=cd_bias
        )
    electrodes = {}
    for electrode_id, mask in masks.electrodes.items():
        overlay, cd_bias = variation(layer_of[electrode_id])
        electrodes[electrode_id] = transform_mask(
            mask, grid, overlay_xy_um=overlay, cd_bias_um=cd_bias
        )
    terminal_regions: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    open_paths = {item.path_id: item for item in manifest.open_paths}
    for path_id, regions in (masks.terminal_regions or {}).items():
        layer_name = layer_of[open_paths[path_id].electrode]
        overlay, _cd_bias = variation(layer_name)
        terminal_regions[path_id] = tuple(
            transform_mask(region, grid, overlay_xy_um=overlay, cd_bias_um=0.0)
            for region in regions
        )
    contact_masks: dict[str, np.ndarray] = {}
    for contact in manifest.contacts:
        overlay, cd_bias = variation(contact.hole_layer)
        contact_masks[contact.contact_id] = transform_mask(
            (masks.contact_masks or {})[contact.contact_id],
            grid,
            overlay_xy_um=overlay,
            cd_bias_um=cd_bias,
        )
    return MaskSet(
        grid=grid,
        layers=layers,
        electrodes=electrodes,
        input_hashes=masks.input_hashes,
        terminal_regions=terminal_regions,
        layer_unassigned_pixels=masks.layer_unassigned_pixels,
        contact_masks=contact_masks,
        contact_topology=build_contact_topology(
            manifest, contact_masks, electrodes, strict=False
        ),
    )


def apply_scenario_to_stack(
    stack_z: dict[str, LayerZ], scenario: Scenario
) -> dict[str, LayerZ]:
    """Re-accumulate z after adding the per-layer thickness deltas."""
    if not scenario.thickness_delta_um and not scenario.planarization_factor:
        return stack_z
    ordered = sorted(stack_z.values(), key=lambda layer: layer.process_order)
    result: dict[str, LayerZ] = {}
    z_bottom = 0.0
    for layer in ordered:
        thickness = layer.thickness_um + scenario.thickness_delta_um.get(layer.layer_name, 0.0)
        thickness = max(thickness, 0.0)
        z_top = z_bottom + thickness
        result[layer.layer_name] = LayerZ(
            layer_name=layer.layer_name,
            process_order=layer.process_order,
            material_type=layer.material_type,
            thickness_um=thickness,
            z_bottom_um=z_bottom,
            z_top_um=z_top,
            planarization_factor=scenario.planarization_factor.get(
                layer.layer_name, layer.planarization_factor
            ),
            patterned=layer.patterned,
        )
        z_bottom = z_top
    return result


def scenario_manifest(manifest: ProductManifest, scenario: Scenario) -> ProductManifest:
    """Manifest copy whose scalar physics values follow the scenario."""
    update: dict[str, object] = {}
    if scenario.embedding_depth_um is not None:
        update["embedding_depth_um"] = scenario.embedding_depth_um
    if scenario.contact_margin_um is not None:
        update["contact_margin_um"] = scenario.contact_margin_um
    if not update:
        return manifest
    return manifest.model_copy(update=update)
