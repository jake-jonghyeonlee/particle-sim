"""Per-scenario union of pair hotspots (never a sum of CA numbers)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from particle_sim.failure.lateral_short import PairResult, lateral_short
from particle_sim.failure.vertical_short import vertical_short
from particle_sim.geometry.contact import DistanceCache, compute_electrode_contacts
from particle_sim.io.masks import MaskSet
from particle_sim.io.stack import LayerZ
from particle_sim.schema import ProductManifest

ANY_MODE = "any"


@dataclass(frozen=True)
class UnionResult:
    failure_mode: str
    pair_ids: tuple[str, ...]
    any_failure: np.ndarray
    pair_count: np.ndarray


def simulate_pairs(
    manifest: ProductManifest,
    masks: MaskSet,
    stack_z: dict[str, LayerZ],
    step_id: str,
    diameter_um: float,
    cache: DistanceCache,
) -> list[PairResult]:
    """Evaluate every electrode pair of the manifest for one step and particle size."""
    needed = sorted(
        {e for pair in manifest.electrode_pairs for e in (pair.electrode_a, pair.electrode_b)}
    )
    contacts = compute_electrode_contacts(
        manifest, masks, stack_z, step_id, diameter_um, cache, needed
    )
    results: list[PairResult] = []
    for pair in manifest.electrode_pairs:
        if pair.failure_mode == "lateral_short":
            results.append(lateral_short(manifest, pair, contacts))
        else:
            results.append(vertical_short(manifest, pair, contacts))
    return results


def union_by_mode(results: list[PairResult]) -> dict[str, UnionResult]:
    """OR pair maps inside each failure mode and over all modes ("any")."""
    if not results:
        raise ValueError("union 계산에는 최소 한 개의 PairResult가 필요합니다")
    shape = results[0].failure.shape
    if any(r.failure.shape != shape for r in results):
        raise ValueError("PairResult의 shape가 서로 다릅니다")

    groups: dict[str, list[PairResult]] = {}
    for result in results:
        groups.setdefault(result.failure_mode, []).append(result)
    groups[ANY_MODE] = list(results)

    unions: dict[str, UnionResult] = {}
    for mode, members in groups.items():
        any_failure = np.zeros(shape, dtype=bool)
        pair_count = np.zeros(shape, dtype=np.uint16)
        for member in members:
            any_failure |= member.failure
            pair_count += member.failure.astype(np.uint16)
        unions[mode] = UnionResult(
            failure_mode=mode,
            pair_ids=tuple(m.pair_id for m in members),
            any_failure=any_failure,
            pair_count=pair_count,
        )
    return unions
