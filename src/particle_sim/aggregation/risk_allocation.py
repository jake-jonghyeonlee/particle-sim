"""겹치는 hotspot을 union 합계와 일치하도록 실패 단위에 배분한다."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from particle_sim.failure.lateral_short import PairResult
from particle_sim.io.masks import GridSpec


@dataclass(frozen=True)
class RiskAllocation:
    union: np.ndarray
    active_count: np.ndarray
    raw_ca_um2: dict[str, float]
    allocated_ca_um2: dict[str, float]
    unique_ca_um2: dict[str, float]


def allocate_risk(results: list[PairResult], grid: GridSpec) -> RiskAllocation:
    """활성 단위 k개가 겹친 픽셀을 각 단위에 1/k씩 배분한다."""
    if not results:
        raise ValueError("위험 배분에는 최소 한 개의 실패 단위가 필요합니다")
    shape = results[0].failure.shape
    if any(result.failure.shape != shape for result in results):
        raise ValueError("실패 단위 hotspot shape가 서로 다릅니다")
    if all(
        type(result.failure) is np.ndarray
        and result.failure.dtype == np.bool_
        and result.failure.flags.c_contiguous
        for result in results
    ):
        active_count = np.zeros(shape, dtype=np.uint32)
        for result in results:
            active_count += result.failure
        failures = (result.failure for result in results)
    else:
        stack = np.stack([result.failure.astype(bool) for result in results], axis=0)
        active_count = stack.sum(axis=0, dtype=np.uint32)
        failures = iter(stack)
    union = active_count > 0
    denominator = np.where(union, active_count, 1).astype(np.float64)
    pixel_area = grid.pixel_area_um2
    raw: dict[str, float] = {}
    allocated: dict[str, float] = {}
    unique: dict[str, float] = {}
    for result, failure in zip(results, failures, strict=True):
        raw[result.pair_id] = float(failure.sum()) * pixel_area
        allocated[result.pair_id] = float(np.sum(failure / denominator)) * pixel_area
        unique[result.pair_id] = float(np.sum(failure & (active_count == 1))) * pixel_area
    return RiskAllocation(
        union=union,
        active_count=active_count,
        raw_ca_um2=raw,
        allocated_ca_um2=allocated,
        unique_ca_um2=unique,
    )
