"""Common view over the evaluated failure units (electrode pairs and open paths)."""

from __future__ import annotations

import numpy as np

from particle_sim.failure.lateral_short import PairResult
from particle_sim.failure.open_failure import simulate_open_sizes
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


class OpenCache:
    """Open results reused across source steps.

    An open only depends on the electrode mask, the particle size and the contact margin, so
    the same geometry must not be evaluated once per source step.
    """

    def __init__(self, manifest: ProductManifest, masks: MaskSet, sizes: list[float]) -> None:
        self._manifest = manifest
        self._masks = masks
        self._sizes = sorted(sizes)
        self._results: dict[tuple[str, float], object] = {}

    def get(self, open_path, diameter_um: float):
        key = (open_path.path_id, float(diameter_um))
        if key not in self._results:
            # every size is computed once, in increasing order, reusing the carry-forward
            batch = simulate_open_sizes(
                self._manifest,
                self._masks,
                open_path,
                self._sizes if float(diameter_um) in self._sizes else [float(diameter_um)],
                contact_margin_um=self._manifest.contact_margin_um,
            )
            for size, result in batch.items():
                self._results[(open_path.path_id, float(size))] = result
        return self._results[key]


def evaluate_units(
    manifest: ProductManifest,
    masks: MaskSet,
    stack_z: dict[str, LayerZ],
    step_id: str,
    diameter_um: float,
    cache: DistanceCache,
    open_cache: OpenCache | None = None,
) -> list[PairResult]:
    """Short pairs and open paths for one source step and particle size."""
    results = simulate_pairs(manifest, masks, stack_z, step_id, diameter_um, cache)
    if manifest.open_paths and open_cache is None:
        open_cache = OpenCache(manifest, masks, [diameter_um])
    for open_path in manifest.open_paths:
        open_result = open_cache.get(open_path, diameter_um)
        results.append(open_as_pair_result(open_result, masks.grid.shape_yx))
    return results
