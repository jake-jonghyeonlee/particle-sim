"""Critical area integration over hotspot maps."""

from __future__ import annotations

import numpy as np

from particle_sim.io.masks import GridSpec


def critical_area_um2(hotspot: np.ndarray, grid: GridSpec) -> float:
    """Spatial integral of a binary or [0, 1] hotspot map, in um^2."""
    values = np.asarray(hotspot)
    if values.shape != grid.shape_yx:
        raise ValueError(f"hotspot shape {values.shape} != grid shape {grid.shape_yx}")
    values = values.astype(np.float64)
    if not np.isfinite(values).all():
        raise ValueError("hotspot에 NaN 또는 무한대가 있습니다")
    if values.min() < 0.0 or values.max() > 1.0:
        raise ValueError("hotspot 값은 0과 1 사이여야 합니다")
    return float(values.sum() * grid.pixel_area_um2)


def um2_to_mm2(value_um2: float) -> float:
    """um^2 -> mm^2."""
    return float(value_um2) * 1e-6


def normalized_ca(ca_um2: float, grid: GridSpec) -> float:
    """CA divided by the analysed ROI area (dimensionless, unit "1")."""
    return float(ca_um2) / grid.roi_area_um2
