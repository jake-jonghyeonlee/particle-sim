"""3-D and post-step contact maps between a particle and electrodes."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np

from particle_sim.errors import InputContractError
from particle_sim.geometry.distance import distance_to_mask_um, distance_to_z_interval_um
from particle_sim.geometry.particles import particle_center_z_um, post_step_eligibility
from particle_sim.geometry.surfaces import surface_height_um
from particle_sim.io.masks import MaskSet
from particle_sim.io.stack import LayerZ
from particle_sim.schema import ProductManifest


@dataclass(frozen=True)
class ContactResult:
    contact: np.ndarray
    margin_um: np.ndarray
    kind: str


def contact_3d(
    distance_xy_um: np.ndarray,
    center_z_um: np.ndarray,
    z_bottom_um: float,
    z_top_um: float,
    radius_um: float,
    contact_margin_um: float,
) -> ContactResult:
    """Contact when r + margin - sqrt(Dxy^2 + Dz^2) >= 0 (Q == R counts as contact)."""
    distance_z = distance_to_z_interval_um(center_z_um, z_bottom_um, z_top_um)
    with np.errstate(invalid="ignore"):
        distance_3d = np.hypot(np.asarray(distance_xy_um, dtype=np.float64), distance_z)
        margin = radius_um + contact_margin_um - distance_3d
    margin = np.where(np.isfinite(distance_3d), margin, -np.inf)
    return ContactResult(
        contact=margin >= 0.0, margin_um=margin.astype(np.float32), kind="3d"
    )


def contact_post_step(
    distance_xy_um: np.ndarray,
    eligible: bool,
    radius_um: float,
    contact_margin_um: float,
) -> ContactResult:
    """A later-formed electrode wraps the protruding particle: only Dxy matters."""
    distance_xy = np.asarray(distance_xy_um, dtype=np.float64)
    if not eligible:
        margin = np.full(distance_xy.shape, -np.inf, dtype=np.float32)
        return ContactResult(
            contact=np.zeros(distance_xy.shape, dtype=bool), margin_um=margin, kind="excluded"
        )
    margin = radius_um + contact_margin_um - distance_xy
    margin = np.where(np.isfinite(distance_xy), margin, -np.inf)
    return ContactResult(
        contact=margin >= 0.0, margin_um=margin.astype(np.float32), kind="post_step"
    )


class DistanceCache:
    """Per-electrode XY distance maps, computed once per product."""

    def __init__(self, masks: MaskSet) -> None:
        self._masks = masks
        self._cache: dict[str, np.ndarray] = {}

    def get(self, electrode_id: str) -> np.ndarray:
        if electrode_id not in self._cache:
            if electrode_id not in self._masks.electrodes:
                raise InputContractError(
                    f"알 수 없는 electrode_id: '{electrode_id}'", field="electrodes"
                )
            self._cache[electrode_id] = distance_to_mask_um(
                self._masks.electrodes[electrode_id], self._masks.grid
            )
        return self._cache[electrode_id]

    @property
    def computed(self) -> int:
        return len(self._cache)


def compute_electrode_contacts(
    manifest: ProductManifest,
    masks: MaskSet,
    stack_z: dict[str, LayerZ],
    step_id: str,
    diameter_um: float,
    cache: DistanceCache,
    electrode_ids: Iterable[str] | None = None,
) -> dict[str, ContactResult]:
    """Contact map of every requested electrode for one source step and particle size."""
    step = next((s for s in manifest.particle_source_steps if s.step_id == step_id), None)
    if step is None:
        raise InputContractError(
            f"알 수 없는 step_id: '{step_id}'",
            product_id=manifest.product_id,
            field="step_id",
        )
    surface = surface_height_um(manifest, step_id, stack_z, masks)
    center_z = particle_center_z_um(surface, diameter_um, manifest.embedding_depth_um)
    eligibility = post_step_eligibility(manifest, step_id, stack_z, diameter_um)
    landing_order = stack_z[step.lands_after_layer].process_order
    radius = diameter_um / 2.0
    electrodes = {e.electrode_id: e for e in manifest.electrodes}
    wanted = list(electrodes) if electrode_ids is None else list(electrode_ids)

    results: dict[str, ContactResult] = {}
    for electrode_id in wanted:
        electrode = electrodes.get(electrode_id)
        if electrode is None:
            raise InputContractError(
                f"알 수 없는 electrode_id: '{electrode_id}'",
                product_id=manifest.product_id,
                field="electrodes",
            )
        layer = stack_z[electrode.layer_name]
        distance_xy = cache.get(electrode_id)
        if layer.process_order <= landing_order:
            results[electrode_id] = contact_3d(
                distance_xy,
                center_z,
                layer.z_bottom_um,
                layer.z_top_um,
                radius,
                manifest.contact_margin_um,
            )
        else:
            results[electrode_id] = contact_post_step(
                distance_xy,
                eligibility.get(layer.layer_name, False),
                radius,
                manifest.contact_margin_um,
            )
    return results
