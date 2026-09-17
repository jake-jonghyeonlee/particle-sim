"""µm 단위 EDT 거리맵과 z 구간 거리 계산."""

from __future__ import annotations

import numpy as np
from scipy.ndimage import distance_transform_edt

from particle_sim.io.masks import GridSpec


def distance_to_mask_um(mask: np.ndarray, grid: GridSpec) -> np.ndarray:
    """mask 경계까지의 µm 단위 유클리드 거리맵을 반환한다.

    mask 내부 픽셀은 0.0, 외부 픽셀은 경계까지의 거리를 가진다.
    빈 mask(모두 False)는 전체 np.inf를 반환한다.
    shape 불일치는 ValueError를 발생시킨다.
    """
    if mask.shape != grid.shape_yx:
        raise ValueError(
            f"mask shape {mask.shape}와 grid.shape_yx {grid.shape_yx}이(가) 일치하지 않습니다"
        )
    if not mask.any():
        return np.full(mask.shape, np.inf, dtype=np.float64)
    sampling = (grid.pixel_size_y_um, grid.pixel_size_x_um)
    return distance_transform_edt(~mask, sampling=sampling).astype(np.float64)


def distance_to_z_interval_um(
    z_center: np.ndarray | float, z_bottom_um: float, z_top_um: float
) -> np.ndarray:
    """z_center가 [z_bottom_um, z_top_um] 구간 밖일 때의 거리를 반환한다.

    max(z_bottom - z, z - z_top, 0)을 elementwise로 계산한다.
    z_top_um <= z_bottom_um이면 ValueError를 발생시킨다.
    """
    if z_top_um <= z_bottom_um:
        raise ValueError(
            f"z_top_um({z_top_um})은 z_bottom_um({z_bottom_um})보다 커야 합니다"
        )
    z = np.asarray(z_center, dtype=np.float64)
    below = z_bottom_um - z
    above = z - z_top_um
    return np.maximum(np.maximum(below, above), 0.0)
