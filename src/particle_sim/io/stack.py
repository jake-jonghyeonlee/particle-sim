"""누적 z stack 계산과 source step 기준 layer 분할."""

from __future__ import annotations

from dataclasses import dataclass

from particle_sim.errors import InputContractError
from particle_sim.schema import ParticleSourceStep, ProductManifest


@dataclass(frozen=True)
class LayerZ:
    """단일 레이어의 z 범위."""

    layer_name: str
    process_order: int
    material_type: str
    thickness_um: float
    z_bottom_um: float
    z_top_um: float
    planarization_factor: float | None
    patterned: bool


def build_stack_z(manifest: ProductManifest) -> dict[str, LayerZ]:
    """process_order 순서로 누적 z를 계산하여 LayerZ dict를 반환한다.

    첫 레이어 z_bottom = 0.0, 이후 각 z_bottom = 이전 z_top.
    z_bottom_um이 명시된 경우 1e-6 이내 일치 여부를 검증한다.
    """
    sorted_layers = sorted(manifest.stack, key=lambda layer: layer.process_order)
    result: dict[str, LayerZ] = {}
    z = 0.0
    for layer in sorted_layers:
        z_top = z + layer.thickness_um
        if layer.z_bottom_um is not None:
            if abs(layer.z_bottom_um - z) > 1e-6:
                raise InputContractError(
                    f"z_bottom_um 불일치: given={layer.z_bottom_um}, computed={z}",
                    product_id=manifest.product_id,
                    field=f"stack[{layer.layer_name}].z_bottom_um",
                )
        result[layer.layer_name] = LayerZ(
            layer_name=layer.layer_name,
            process_order=layer.process_order,
            material_type=layer.material_type,
            thickness_um=layer.thickness_um,
            z_bottom_um=z,
            z_top_um=z_top,
            planarization_factor=layer.planarization_factor,
            patterned=layer.patterned,
        )
        z = z_top
    return result


def landing_surface_z_um(
    manifest: ProductManifest, step_id: str, stack_z: dict[str, LayerZ]
) -> float:
    """step의 landing surface z (lands_after_layer의 z_top)를 반환한다."""
    for step in manifest.particle_source_steps:
        if step.step_id == step_id:
            landing_layer = step.lands_after_layer
            if landing_layer not in stack_z:
                raise InputContractError(
                    f"알 수 없는 layer: {landing_layer}",
                    product_id=manifest.product_id,
                    field="step_id",
                )
            return stack_z[landing_layer].z_top_um
    raise InputContractError(
        f"알 수 없는 step_id: {step_id}",
        product_id=manifest.product_id,
        field="step_id",
    )


def layers_up_to(step: ParticleSourceStep, stack_z: dict[str, LayerZ]) -> list[LayerZ]:
    """landing layer 포함, process_order <= landing order인 레이어 목록."""
    landing_order = stack_z[step.lands_after_layer].process_order
    return [lz for lz in stack_z.values() if lz.process_order <= landing_order]


def layers_after(step: ParticleSourceStep, stack_z: dict[str, LayerZ]) -> list[LayerZ]:
    """landing layer 이후, process_order > landing order인 레이어 목록."""
    landing_order = stack_z[step.lands_after_layer].process_order
    return [lz for lz in stack_z.values() if lz.process_order > landing_order]
