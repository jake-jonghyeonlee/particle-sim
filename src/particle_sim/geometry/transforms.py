"""Apply one process-variation scenario to the masks and the layer stack."""

from __future__ import annotations

import numpy as np
from scipy.ndimage import distance_transform_edt, map_coordinates

from particle_sim.io.masks import GridSpec, MaskSet
from particle_sim.io.stack import LayerZ
from particle_sim.schema import ProductManifest
from particle_sim.uncertainty.sampler import Scenario


def signed_distance_um(mask: np.ndarray, grid: GridSpec) -> np.ndarray:
    """Negative inside the mask, positive outside, in micrometres."""
    sampling = (grid.pixel_size_y_um, grid.pixel_size_x_um)
    if not mask.any():
        return np.full(mask.shape, np.inf, dtype=np.float64)
    if mask.all():
        return np.full(mask.shape, -np.inf, dtype=np.float64)
    outside = distance_transform_edt(~mask, sampling=sampling)
    inside = distance_transform_edt(mask, sampling=sampling)
    return outside - inside


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
    return shifted <= -cd_bias_um


def apply_scenario_to_masks(
    manifest: ProductManifest, masks: MaskSet, scenario: Scenario
) -> MaskSet:
    """Overlay, CD bias and the particle centre offset applied to every mask."""
    offset_y, offset_x = scenario.center_offset_xy_um
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
    return MaskSet(
        grid=grid, layers=layers, electrodes=electrodes, input_hashes=masks.input_hashes
    )


def apply_scenario_to_stack(
    stack_z: dict[str, LayerZ], scenario: Scenario
) -> dict[str, LayerZ]:
    """Re-accumulate z after adding the per-layer thickness deltas."""
    if not scenario.thickness_delta_um:
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
