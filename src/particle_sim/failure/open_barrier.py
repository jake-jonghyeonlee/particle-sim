"""Certified dual-background judgement for explicit open paths.

The certificate requires a connected baseline, one effective electrode pixel in
each terminal, opposite non-corner image sides, usable positive finite spacing,
and exactly two distinct 8-connected background components on the two frame
arcs. A one-pixel background frame surrounds the *original* image. Its only
foreground pixels are the two outward 4-neighbor ports, each attached to exactly
one terminal. No second exterior ring is present.

Foreground 4-connectivity and background 8-connectivity are complementary in a
rectangle: the ports disconnect exactly when the intervening frame arcs join in
background. The removed particle stamp is 4-connected because walking toward
its center never increases distance. It can join the initial distinct frame
components exactly when it intersects both of their closed 8-neighborhoods.
Other initial background components, including holes, have no edges to either
frame component; any new connection must pass through that connected stamp.
Covering either singleton terminal disconnects its fixed port as well. Thus the
geometric critical reach is min(max(distance to both sides), distance to A,
distance to B). The caller repairs floating ties against the actual footprint
and unchanged local connectivity judgement, including footprint-box rounding.

The EDT guarantee is limited to a conservative finite arithmetic envelope.
SciPy's feature transform can form products beyond squared distances, so both
physical spacing and image extent stay well inside the floating-point range.
Unsupported spacing or extent falls back to local judgement. No result is ever
inferred from a failed certificate.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import binary_dilation, distance_transform_edt, label


@dataclass(frozen=True)
class Ineligible:
    reason: str


@dataclass(frozen=True)
class PreparedBarrier:
    critical: np.ndarray
    pixel_size_y_um: float
    pixel_size_x_um: float

    def ambiguous(self, reach_um: float) -> np.ndarray:
        height, width = self.critical.shape
        tolerance = 32 * np.finfo(float).eps * max(
            1.0,
            abs(reach_um),
            (height + 2) * self.pixel_size_y_um,
            (width + 2) * self.pixel_size_x_um,
        )
        return np.abs(self.critical - reach_um) <= tolerance


_EIGHT = np.ones((3, 3), dtype=bool)


def _port(point: tuple[int, int], height: int, width: int) -> tuple[int, int] | None:
    y, x = point
    if y == 0 and 0 < x < width - 1:
        return 0, x + 1
    if y == height - 1 and 0 < x < width - 1:
        return height + 1, x + 1
    if x == 0 and 0 < y < height - 1:
        return y + 1, 0
    if x == width - 1 and 0 < y < height - 1:
        return y + 1, width + 1
    return None


def _opposite(a: tuple[int, int], b: tuple[int, int], h: int, w: int) -> bool:
    return (a[0] == 0 and b[0] == h - 1) or (
        b[0] == 0 and a[0] == h - 1
    ) or (a[1] == 0 and b[1] == w - 1) or (
        b[1] == 0 and a[1] == w - 1
    )


def _distance_to_sites(sites: np.ndarray, sy: float, sx: float) -> np.ndarray:
    iy, ix = distance_transform_edt(
        ~sites, sampling=(sy, sx), return_distances=False, return_indices=True
    )
    yy, xx = np.ogrid[: sites.shape[0], : sites.shape[1]]
    return np.hypot((iy - yy) * sy, (ix - xx) * sx)


def prepare_barrier(
    electrode: np.ndarray,
    terminal_a: np.ndarray,
    terminal_b: np.ndarray,
    pixel_size_y_um: float,
    pixel_size_x_um: float,
) -> PreparedBarrier | Ineligible:
    """Prepare one call's certified barrier or explain why local judgement is needed.

    The caller has already resolved terminals and checked the connected baseline.
    """
    height, width = electrode.shape
    if height < 3 or width < 3:
        return Ineligible("small_image")
    a = np.argwhere(electrode & terminal_a)
    b = np.argwhere(electrode & terminal_b)
    if len(a) != 1 or len(b) != 1:
        return Ineligible("effective_terminal_not_singleton")
    point_a = tuple(int(v) for v in a[0])
    point_b = tuple(int(v) for v in b[0])
    port_a = _port(point_a, height, width)
    port_b = _port(point_b, height, width)
    if port_a is None or port_b is None or not _opposite(
        point_a, point_b, height, width
    ):
        return Ineligible("terminals_not_opposite_noncorner_sides")

    sy, sx = float(pixel_size_y_um), float(pixel_size_x_um)
    extent_y, extent_x = (height + 2) * sy, (width + 2) * sx
    # Feature-transform site selection can involve length times squared distance.
    # Keep both ends far from that cubic overflow/underflow regime.
    lower, upper = 1e-30, 1e30
    if not all(np.isfinite(v) for v in (sy, sx, extent_y, extent_x)) or not (
        lower < min(sy, sx) <= max(sy, sx) < upper
        and max(extent_y, extent_x) < upper
        and max(sy, sx) / min(sy, sx) <= 1e12
    ):
        return Ineligible("unsafe_pixel_arithmetic")

    padded = np.zeros((height + 2, width + 2), dtype=bool)
    padded[1:-1, 1:-1] = electrode
    padded[port_a] = True
    padded[port_b] = True
    labels, _ = label(~padded, structure=_EIGHT)

    # Clockwise perimeter. The two fixed ports split it into exactly two arcs.
    perimeter = (
        [(0, x) for x in range(width + 2)]
        + [(y, width + 1) for y in range(1, height + 2)]
        + [(height + 1, x) for x in range(width, -1, -1)]
        + [(y, 0) for y in range(height, 0, -1)]
    )
    cut = [i for i, point in enumerate(perimeter) if point in (port_a, port_b)]
    if len(cut) != 2:
        return Ineligible("invalid_frame_ports")
    first, second = cut
    arcs = (perimeter[first + 1 : second], perimeter[second + 1 :] + perimeter[:first])
    arc_labels = [{int(labels[point]) for point in arc} for arc in arcs]
    if any(len(values) != 1 or 0 in values for values in arc_labels):
        return Ineligible("frame_arc_not_single_component")
    side_a = next(iter(arc_labels[0]))
    side_b = next(iter(arc_labels[1]))
    if side_a == side_b:
        return Ineligible("frame_arcs_connected")

    distances = []
    for side in (side_a, side_b):
        sites = binary_dilation(labels == side, structure=_EIGHT)[1:-1, 1:-1]
        if not sites.any():
            return Ineligible("empty_side_sites")
        distances.append(_distance_to_sites(sites, sy, sx))
    yy, xx = np.ogrid[:height, :width]
    terminal_distances = [
        np.hypot((yy - y) * sy, (xx - x) * sx)
        for y, x in (point_a, point_b)
    ]
    critical = np.minimum(
        np.maximum(distances[0], distances[1]),
        np.minimum(terminal_distances[0], terminal_distances[1]),
    )
    if not np.isfinite(critical).all():
        return Ineligible("nonfinite_critical_distance")
    critical.setflags(write=False)
    return PreparedBarrier(critical, sy, sx)
