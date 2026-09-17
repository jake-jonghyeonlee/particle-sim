"""legacy vs EDT binary map 비교 및 차이 원인 분류."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from particle_sim.geometry.legacy import (
    edt_short_hotspot,
    legacy_effective_radius_um,
    legacy_short_hotspot,
)
from particle_sim.io.masks import GridSpec
from particle_sim.schema import LegacyKernel

CAUSE_CODES: dict[str, int] = {
    "same": 0,
    "kernel_quantization": 1,
    "kernel_shape": 2,
    "unexplained": 3,
}


@dataclass(frozen=True)
class PixelComparison:
    """두 binary map의 픽셀 단위 비교 결과."""

    equal: bool
    n_diff: int
    only_in_reference: int
    only_in_candidate: int
    cause_counts: dict[str, int]
    cause_map: np.ndarray


def compare_binary_maps(reference: np.ndarray, candidate: np.ndarray) -> PixelComparison:
    """두 binary map을 비교하여 PixelComparison을 반환한다.

    shape 불일치 -> ValueError.
    차이 픽셀은 모두 "unexplained"로 분류한다.
    """
    if reference.shape != candidate.shape:
        raise ValueError(
            f"shape 불일치: reference {reference.shape} != candidate {candidate.shape}"
        )
    ref = reference.astype(bool)
    cand = candidate.astype(bool)
    diff = ref != cand
    n_diff = int(diff.sum())
    only_ref = int((ref & ~cand).sum())
    only_cand = int((~ref & cand).sum())
    equal = n_diff == 0

    cause_map = np.zeros(reference.shape, dtype=np.uint8)
    if n_diff > 0:
        cause_map[diff] = CAUSE_CODES["unexplained"]

    cause_counts: dict[str, int] = {
        "kernel_quantization": 0,
        "kernel_shape": 0,
        "unexplained": n_diff,
    }
    return PixelComparison(
        equal=equal,
        n_diff=n_diff,
        only_in_reference=only_ref,
        only_in_candidate=only_cand,
        cause_counts=cause_counts,
        cause_map=cause_map,
    )


def explain_legacy_vs_edt(
    mask_a: np.ndarray,
    mask_b: np.ndarray,
    diameter_um: float,
    grid: GridSpec,
    kernel: LegacyKernel,
) -> PixelComparison:
    """legacy_short_hotspot와 edt_short_hotspot를 비교하고 차이 원인을 분류한다.

    ref  = legacy_short_hotspot(diameter_um)
    cand = edt_short_hotspot(r = diameter_um / 2)
    q    = edt_short_hotspot(r_eff = legacy_effective_radius_um)

    ref != cand인 픽셀에 대해:
      cand != q -> kernel_quantization
      ref  != q -> kernel_shape
      그 외     -> unexplained
    """
    ref = legacy_short_hotspot(mask_a, mask_b, diameter_um, grid, kernel)
    r = diameter_um / 2.0
    cand = edt_short_hotspot(mask_a, mask_b, r, grid)
    r_eff = legacy_effective_radius_um(diameter_um, grid, kernel)
    q = edt_short_hotspot(mask_a, mask_b, r_eff, grid)

    diff = ref != cand
    n_diff = int(diff.sum())
    only_ref = int((ref & ~cand).sum())
    only_cand = int((~ref & cand).sum())
    equal = n_diff == 0

    cause_map = np.zeros(ref.shape, dtype=np.uint8)
    n_quant = 0
    n_shape = 0
    n_unexplained = 0

    if n_diff > 0:
        diff_idx = np.where(diff)
        cand_ne_q = cand[diff_idx] != q[diff_idx]
        ref_ne_q = ref[diff_idx] != q[diff_idx]

        quant_mask = cand_ne_q
        shape_mask = ~cand_ne_q & ref_ne_q
        unexp_mask = ~cand_ne_q & ~ref_ne_q

        cause_map[diff_idx] = CAUSE_CODES["unexplained"]
        if quant_mask.any():
            q_idx = tuple(diff_idx[i][quant_mask] for i in range(len(diff_idx)))
            cause_map[q_idx] = CAUSE_CODES["kernel_quantization"]
            n_quant = int(quant_mask.sum())
        if shape_mask.any():
            s_idx = tuple(diff_idx[i][shape_mask] for i in range(len(diff_idx)))
            cause_map[s_idx] = CAUSE_CODES["kernel_shape"]
            n_shape = int(shape_mask.sum())
        n_unexplained = int(unexp_mask.sum())

    cause_counts: dict[str, int] = {
        "kernel_quantization": n_quant,
        "kernel_shape": n_shape,
        "unexplained": n_unexplained,
    }
    return PixelComparison(
        equal=equal,
        n_diff=n_diff,
        only_in_reference=only_ref,
        only_in_candidate=only_cand,
        cause_counts=cause_counts,
        cause_map=cause_map,
    )
