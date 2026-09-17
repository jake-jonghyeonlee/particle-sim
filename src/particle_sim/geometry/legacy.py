"""기존 5x5 Erode/Dilate 기반 hotspot 재현."""

from __future__ import annotations

import math

import numpy as np
from scipy.ndimage import binary_dilation, binary_erosion

from particle_sim.errors import InputContractError
from particle_sim.geometry.distance import distance_to_mask_um
from particle_sim.io.masks import GridSpec
from particle_sim.schema import LegacyKernel


def legacy_structure(kernel: LegacyKernel) -> np.ndarray:
    """kernel shape에 따른 bool 구조를 반환한다.

    square: np.ones((size, size), bool)
    disk: (dy, dx) in [-h, h] where dy*dy + dx*dx <= h*(h+1), h = size // 2
    """
    size = kernel.size
    h = size // 2
    if kernel.shape == "square":
        return np.ones((size, size), dtype=bool)
    # disk
    ys, xs = np.mgrid[-h : h + 1, -h : h + 1]
    return ys * ys + xs * xs <= h * (h + 1)


def legacy_iterations(diameter_um: float, grid: GridSpec, kernel: LegacyKernel) -> int:
    """기존 방식의 반복 횟수를 계산한다.

    grid.is_isotropic이 False이면 InputContractError(field="coordinate_system.pixel_size").
    radius_px = diameter_um / 2 / pixel_size
    x = radius_px / (size // 2)
    round -> floor(x + 0.5) (half-up, banker's 아님)
    ceil  -> ceil(x - 1e-9)
    floor -> floor(x + 1e-9)
    """
    if not grid.is_isotropic:
        raise InputContractError(
            "legacy 반복 규칙은 등방성(isotropic) grid만 지원합니다",
            field="coordinate_system.pixel_size",
        )
    pixel_size = grid.pixel_size_y_um
    radius_px = diameter_um / 2.0 / pixel_size
    h = kernel.size // 2
    x = radius_px / h
    if kernel.iterations == "round":
        return int(math.floor(x + 0.5))
    if kernel.iterations == "ceil":
        return int(math.ceil(x - 1e-9))
    # floor
    return int(math.floor(x + 1e-9))


def legacy_effective_radius_um(
    diameter_um: float, grid: GridSpec, kernel: LegacyKernel
) -> float:
    """기존 방식의 유효 반경(µm)을 반환한다.

    iterations * (size // 2) * pixel_size
    """
    n = legacy_iterations(diameter_um, grid, kernel)
    h = kernel.size // 2
    return float(n * h * grid.pixel_size_y_um)


def legacy_dilate(
    mask: np.ndarray, diameter_um: float, grid: GridSpec, kernel: LegacyKernel
) -> np.ndarray:
    """기존 binary_dilation을 적용한다.

    n == 0이면 mask.copy(), 아니면 scipy.ndimage.binary_dilation(
    mask, structure=legacy_structure(kernel), iterations=n, border_value=0).
    """
    n = legacy_iterations(diameter_um, grid, kernel)
    if n == 0:
        return mask.copy()
    structure = legacy_structure(kernel)
    return binary_dilation(mask, structure=structure, iterations=n, border_value=0)


def legacy_erode(
    mask: np.ndarray, diameter_um: float, grid: GridSpec, kernel: LegacyKernel
) -> np.ndarray:
    """기존 binary_erosion을 적용한다.

    n == 0이면 mask.copy(), 아니면 scipy.ndimage.binary_erosion(
    mask, structure=legacy_structure(kernel), iterations=n, border_value=0).
    """
    n = legacy_iterations(diameter_um, grid, kernel)
    if n == 0:
        return mask.copy()
    structure = legacy_structure(kernel)
    return binary_erosion(mask, structure=structure, iterations=n, border_value=0)


def legacy_short_hotspot(
    mask_a: np.ndarray,
    mask_b: np.ndarray,
    diameter_um: float,
    grid: GridSpec,
    kernel: LegacyKernel,
) -> np.ndarray:
    """legacy_dilate(a) & legacy_dilate(b)를 반환한다.

    vertical pair도 동일한 2-D 투영을 사용한다.
    """
    dilated_a = legacy_dilate(mask_a, diameter_um, grid, kernel)
    dilated_b = legacy_dilate(mask_b, diameter_um, grid, kernel)
    return dilated_a & dilated_b


def edt_short_hotspot(
    mask_a: np.ndarray, mask_b: np.ndarray, radius_um: float, grid: GridSpec
) -> np.ndarray:
    """EDT 기반 hotspot: dist_a <= r AND dist_b <= r.

    (distance_to_mask_um(a) <= radius_um) & (distance_to_mask_um(b) <= radius_um)
    """
    dist_a = distance_to_mask_um(mask_a, grid)
    dist_b = distance_to_mask_um(mask_b, grid)
    return (dist_a <= radius_um) & (dist_b <= radius_um)
