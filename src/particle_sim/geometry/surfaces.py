"""Particle landing surface height (flat / conformal / partial planarization)."""

from __future__ import annotations

import numpy as np

from particle_sim.errors import InputContractError
from particle_sim.io.masks import MaskSet
from particle_sim.io.stack import LayerZ, layers_up_to
from particle_sim.schema import ProductManifest


def surface_height_um(
    manifest: ProductManifest,
    step_id: str,
    stack_z: dict[str, LayerZ],
    masks: MaskSet,
) -> np.ndarray:
    """Return the landing surface height (um) of a particle source step."""
    step = next((s for s in manifest.particle_source_steps if s.step_id == step_id), None)
    if step is None:
        raise InputContractError(
            f"알 수 없는 step_id: '{step_id}'",
            product_id=manifest.product_id,
            field="step_id",
        )
    landing = stack_z[step.lands_after_layer]
    surface = np.full(masks.grid.shape_yx, landing.z_top_um, dtype=np.float64)
    if step.surface_model == "flat":
        return surface

    below = layers_up_to(step, stack_z)
    for layer in below:
        if not layer.patterned:
            continue
        retention = 1.0
        if step.surface_model == "partial":
            for upper in below:
                if upper.process_order <= layer.process_order:
                    continue
                if upper.material_type == "conductor":
                    continue
                if upper.planarization_factor is None:
                    raise InputContractError(
                        f"partial 표면 모델에는 레이어 '{upper.layer_name}'의 "
                        "planarization_factor가 필요합니다",
                        product_id=manifest.product_id,
                        field=f"stack[{upper.layer_name}].planarization_factor",
                    )
                retention *= 1.0 - upper.planarization_factor
        missing = ~masks.layers[layer.layer_name]
        surface -= layer.thickness_um * missing * retention
    return surface
