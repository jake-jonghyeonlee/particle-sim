"""Common view over the evaluated failure units (electrode pairs and open paths)."""

from __future__ import annotations

import numpy as np

from particle_sim.failure.lateral_short import PairResult
from particle_sim.failure.open_failure import simulate_open
from particle_sim.failure.union import simulate_pairs
from particle_sim.geometry.contact import DistanceCache
from particle_sim.io.masks import MaskSet
from particle_sim.io.stack import LayerZ
from particle_sim.schema import ProductManifest

OPEN_MODE = "lateral_open"
NOT_APPLICABLE = "N/A"


def open_as_pair_result(result, shape_yx: tuple[int, int]) -> PairResult:
    """Wrap an OpenResult so unions and CA tables treat every unit the same way."""
    return PairResult(
        pair_id=result.path_id,
        failure_mode=OPEN_MODE,
        electrode_a=result.electrode,
        electrode_b=NOT_APPLICABLE,
        failure=result.failure,
        contact_a=result.failure,
        contact_b=np.zeros(shape_yx, dtype=bool),
        pair_margin_um=np.full(shape_yx, np.nan, dtype=np.float32),
        classification=result.classification,
    )


def unit_ids(manifest: ProductManifest) -> list[str]:
    """Pair ids followed by open path ids, in manifest order."""
    return [p.pair_id for p in manifest.electrode_pairs] + [
        o.path_id for o in manifest.open_paths
    ]


def unit_modes(manifest: ProductManifest) -> list[str]:
    return [p.failure_mode for p in manifest.electrode_pairs] + [
        OPEN_MODE for _ in manifest.open_paths
    ]


def evaluate_units(
    manifest: ProductManifest,
    masks: MaskSet,
    stack_z: dict[str, LayerZ],
    step_id: str,
    diameter_um: float,
    cache: DistanceCache,
) -> list[PairResult]:
    """Short pairs and open paths for one source step and particle size."""
    results = simulate_pairs(manifest, masks, stack_z, step_id, diameter_um, cache)
    for open_path in manifest.open_paths:
        open_result = simulate_open(
            manifest, masks, open_path, diameter_um, contact_margin_um=manifest.contact_margin_um
        )
        results.append(open_as_pair_result(open_result, masks.grid.shape_yx))
    return results
