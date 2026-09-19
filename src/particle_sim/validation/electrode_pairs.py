"""다른 net 전극 조합이 계산 또는 명시적 제외로 모두 해소됐는지 검사한다."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import TYPE_CHECKING

import numpy as np

from particle_sim.schema import ProductManifest

if TYPE_CHECKING:
    from particle_sim.io.masks import MaskSet
    from particle_sim.io.stack import LayerZ


@dataclass(frozen=True)
class ElectrodePairAuditRow:
    electrode_a: str
    electrode_b: str
    layer_a: str
    layer_b: str
    failure_mode: str
    status: str
    pair_id: str | None
    exclusion_reason: str | None
    min_xy_distance_um: float | None = None
    min_z_gap_um: float | None = None
    min_3d_distance_um: float | None = None


def audit_electrode_pairs(
    manifest: ProductManifest,
    masks: MaskSet | None = None,
    stack_z: dict[str, LayerZ] | None = None,
) -> list[ElectrodePairAuditRow]:
    """모든 다른 net 조합을 결정적 순서로 선언/제외/미결정으로 분류한다."""
    electrodes = {item.electrode_id: item for item in manifest.electrodes}
    declared = {
        tuple(sorted((item.electrode_a, item.electrode_b))): item
        for item in manifest.electrode_pairs
    }
    excluded = {
        tuple(sorted((item.electrode_a, item.electrode_b))): item
        for item in manifest.excluded_electrode_pairs
    }
    distance_maps: dict[str, np.ndarray] = {}
    rows: list[ElectrodePairAuditRow] = []
    for electrode_a, electrode_b in combinations(sorted(electrodes), 2):
        ea = electrodes[electrode_a]
        eb = electrodes[electrode_b]
        if ea.net_id == eb.net_id:
            continue
        key = (electrode_a, electrode_b)
        pair = declared.get(key)
        exclusion = excluded.get(key)
        mode = "lateral_short" if ea.layer_name == eb.layer_name else "vertical_short"
        min_xy: float | None = None
        min_z: float | None = None
        min_3d: float | None = None
        if masks is not None:
            from particle_sim.geometry.distance import distance_to_mask_um

            if electrode_a not in distance_maps:
                distance_maps[electrode_a] = distance_to_mask_um(
                    masks.electrodes[electrode_a], masks.grid
                )
            distance = distance_maps[electrode_a]
            min_xy = float(distance[masks.electrodes[electrode_b]].min())
        if stack_z is not None:
            layer_a = stack_z[ea.layer_name]
            layer_b = stack_z[eb.layer_name]
            min_z = max(
                layer_a.z_bottom_um - layer_b.z_top_um,
                layer_b.z_bottom_um - layer_a.z_top_um,
                0.0,
            )
        if min_xy is not None and min_z is not None:
            min_3d = float((min_xy**2 + min_z**2) ** 0.5)
        rows.append(
            ElectrodePairAuditRow(
                electrode_a=electrode_a,
                electrode_b=electrode_b,
                layer_a=ea.layer_name,
                layer_b=eb.layer_name,
                failure_mode=mode,
                status="declared" if pair else "excluded" if exclusion else "unresolved",
                pair_id=None if pair is None else pair.pair_id,
                exclusion_reason=None if exclusion is None else exclusion.reason,
                min_xy_distance_um=min_xy,
                min_z_gap_um=min_z,
                min_3d_distance_um=min_3d,
            )
        )
    return rows


def unresolved_electrode_pairs(manifest: ProductManifest) -> list[ElectrodePairAuditRow]:
    return [row for row in audit_electrode_pairs(manifest) if row.status == "unresolved"]
