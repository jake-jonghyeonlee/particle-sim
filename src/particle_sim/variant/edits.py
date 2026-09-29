"""회피안 검토용 국소 마스크 편집.

전역 CD bias나 overlay 이동이 아니라, 지정한 ROI 안의 전극 픽셀만 제거한다.
금속을 없애는 편집이므로 호출자는 같은 픽셀을 해당 layer 마스크에서도 지운다.

- remove_region: ROI 안의 전극 픽셀을 모두 제거한다.
- trim_near: ROI 안에서 상대 전극과의 픽셀 중심 거리가 (ROI 안 현재 최소 거리 + amount)보다
  가까운 전극 픽셀만 제거해 그 구간의 간격만 amount 이상 넓힌다.
- remove_overlap: ROI 안에서 상대 전극과 평면 투영이 겹치거나 amount 이내인 전극 픽셀을 제거한다.

거리는 scipy EDT로 픽셀 중심 사이를 um 단위로 잰다. 결과는 픽셀 격자로 양자화된다.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import distance_transform_edt

from particle_sim.io.masks import GridSpec

OPS = ("remove_region", "trim_near", "remove_overlap")


@dataclass(frozen=True)
class EditResult:
    """편집 뒤 전극 마스크와 제거 픽셀, 전후 측정값."""

    electrode: np.ndarray
    removed: np.ndarray
    min_distance_before_um: float | None
    min_distance_after_um: float | None
    overlap_px_before: int | None
    overlap_px_after: int | None


def roi_mask(grid: GridSpec, roi_um: tuple[float, float, float, float]) -> np.ndarray:
    """(x0, y0, x1, y1) um ROI를 픽셀 중심이 안에 드는 bool 마스크로 바꾼다(좌상단 원점, y 아래)."""
    x0, y0, x1, y1 = roi_um
    if not (x1 > x0 and y1 > y0):
        raise ValueError(f"ROI는 x1 > x0, y1 > y0이어야 합니다: {roi_um}")
    h, w = grid.shape_yx
    xs = (np.arange(w) + 0.5) * grid.pixel_size_x_um
    ys = (np.arange(h) + 0.5) * grid.pixel_size_y_um
    inside_x = (xs >= x0) & (xs <= x1)
    inside_y = (ys >= y0) & (ys <= y1)
    return inside_y[:, None] & inside_x[None, :]


def distance_to_um(target: np.ndarray, grid: GridSpec) -> np.ndarray:
    """모든 픽셀에서 target 픽셀 중심까지의 거리(um). target 픽셀은 0."""
    if not target.any():
        return np.full(target.shape, np.inf)
    return distance_transform_edt(
        ~target, sampling=(grid.pixel_size_y_um, grid.pixel_size_x_um)
    )


def _min_in(values: np.ndarray, where: np.ndarray) -> float | None:
    return float(values[where].min()) if where.any() else None


def apply_edit(op: str, electrode: np.ndarray, grid: GridSpec, roi: np.ndarray, *,
               target: np.ndarray | None = None, amount_um: float = 0.0) -> EditResult:
    """전극 마스크에 국소 편집 하나를 적용한다. 원본 배열은 바꾸지 않는다."""
    if op not in OPS:
        raise ValueError(f"알 수 없는 편집 {op!r}: {OPS} 중 하나여야 합니다")
    if amount_um < 0:
        raise ValueError("amount_um은 0 이상이어야 합니다")
    electrode = np.asarray(electrode, dtype=bool)
    local = electrode & roi
    if not local.any():
        raise ValueError("ROI 안에 편집할 전극 픽셀이 없습니다")
    distance = None
    if op == "remove_region":
        removed = local
    else:
        if target is None or not target.any():
            raise ValueError(f"{op}에는 비어 있지 않은 상대 전극 마스크가 필요합니다")
        distance = distance_to_um(np.asarray(target, dtype=bool), grid)
        if op == "trim_near":
            if amount_um <= 0:
                raise ValueError("trim_near의 amount_um은 0보다 커야 합니다")
            nearest = float(distance[local].min())
            removed = local & (distance < nearest + amount_um)
        else:
            removed = local & (distance <= amount_um)
            if not removed.any():
                raise ValueError("ROI 안에 상대 전극과 겹치는 전극 픽셀이 없습니다")
    edited = electrode & ~removed
    before = after = None
    overlap_before = overlap_after = None
    if distance is not None:
        before = _min_in(distance, local)
        after = _min_in(distance, edited & roi)
        overlap_before = int((electrode & (distance == 0)).sum())
        overlap_after = int((edited & (distance == 0)).sum())
    return EditResult(edited, removed, before, after, overlap_before, overlap_after)
