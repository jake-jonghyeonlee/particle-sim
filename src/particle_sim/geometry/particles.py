"""Spherical particle height model and post-step layer contact eligibility."""

from __future__ import annotations

import numpy as np

from particle_sim.errors import InputContractError
from particle_sim.io.stack import LayerZ, layers_after
from particle_sim.schema import ProductManifest


def particle_center_z_um(
    surface_z_um: np.ndarray, diameter_um: float, embedding_depth_um: float
) -> np.ndarray:
    """Sphere centre height: surface + d/2 - embedding depth."""
    if diameter_um <= 0:
        raise InputContractError(
            f"이물 직경은 0보다 커야 합니다: {diameter_um}", field="particle_sizes_um"
        )
    if embedding_depth_um < 0 or embedding_depth_um > diameter_um:
        raise InputContractError(
            f"매립 깊이는 0 이상 직경 이하여야 합니다: {embedding_depth_um} (d={diameter_um})",
            field="embedding_depth_um",
        )
    return np.asarray(surface_z_um, dtype=np.float64) + diameter_um / 2.0 - embedding_depth_um


def post_step_eligibility(
    manifest: ProductManifest,
    step_id: str,
    stack_z: dict[str, LayerZ],
    diameter_um: float,
) -> dict[str, bool]:
    """Can an electrode layer formed after the source step touch the particle?"""
    step = next((s for s in manifest.particle_source_steps if s.step_id == step_id), None)
    if step is None:
        raise InputContractError(
            f"알 수 없는 step_id: '{step_id}'",
            product_id=manifest.product_id,
            field="step_id",
        )
    policy = manifest.post_step_policy
    after = layers_after(step, stack_z)
    result: dict[str, bool] = {}
    height = diameter_um - manifest.embedding_depth_um
    covered = False
    for layer in after:
        if layer.material_type == "conductor":
            if policy.mode == "exclude":
                result[layer.layer_name] = False
            elif policy.mode == "wrap_contact":
                result[layer.layer_name] = True
            else:
                result[layer.layer_name] = height > 0.0 and not covered
            continue
        if policy.mode != "protrusion":
            continue
        if policy.coverage_break_factor is None:
            raise InputContractError(
                "protrusion 모드에는 coverage_break_factor가 필요합니다",
                product_id=manifest.product_id,
                field="post_step_policy.coverage_break_factor",
            )
        if layer.planarization_factor is None:
            raise InputContractError(
                f"protrusion 모드에는 레이어 '{layer.layer_name}'의 planarization_factor가 "
                "필요합니다",
                product_id=manifest.product_id,
                field=f"stack[{layer.layer_name}].planarization_factor",
            )
        if height <= policy.coverage_break_factor * layer.thickness_um:
            covered = True
        height = max(0.0, height - layer.thickness_um * layer.planarization_factor)
    return result
