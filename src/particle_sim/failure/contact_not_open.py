"""CNO (Contact Not Open): a particle blocks a contact·via hole between two electrodes.

Synthetic criterion: a hole fails when a margin-free particle disc covers at least half
of its entire TR-cut component. A link fails when any one of its holes fails.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from particle_sim.failure.open_failure import particle_footprint
from particle_sim.geometry.contacts import ContactLink
from particle_sim.io.masks import GridSpec

CNO_MODE = "contact_not_open"
CNO_CONFIRMED = "cno_at_least_half_cut_hole"
CONTACT_LOST = "contact_lost"
CNO_COVERAGE_RULE = "at_least_half_cut_hole_any_hole"


@dataclass(frozen=True)
class CnoResult:
    unit_id: str
    failure: np.ndarray
    classification: str


def hole_blocked_centres(
    ys: np.ndarray, xs: np.ndarray, shape_yx: tuple[int, int], stamp: np.ndarray
) -> tuple[slice, slice, np.ndarray]:
    """Centres covering at least half the physical cut hole, in a local window."""
    if len(ys) == 0 or len(ys) != len(xs):
        raise ValueError("physical cut hole must contain matching, nonempty coordinates")
    half_y, half_x = stamp.shape[0] // 2, stamp.shape[1] // 2
    top, bottom = max(0, int(ys.min()) - half_y), min(shape_yx[0], int(ys.max()) + half_y + 1)
    left, right = max(0, int(xs.min()) - half_x), min(shape_yx[1], int(xs.max()) + half_x + 1)
    # Each nonempty stamp row is represented by one or more contiguous runs.
    # Translating a run for a hole pixel adds one to a horizontal range of
    # candidate centres. A difference image plus an exact integer prefix sum
    # avoids materializing the hole-by-footprint correlation tensor.
    width = right - left
    diff = np.zeros((bottom - top, width + 1), dtype=np.int32)
    stamp_y, stamp_x = np.nonzero(stamp)
    half_x = stamp.shape[1] // 2
    row_runs: dict[int, list[tuple[int, int]]] = {}
    for sy in np.unique(stamp_y):
        run_x = stamp_x[stamp_y == sy]
        breaks = np.flatnonzero(np.diff(run_x) != 1) + 1
        for run in np.split(run_x, breaks):
            row_runs[int(sy)] = row_runs.get(int(sy), []) + [(int(run[0]), int(run[-1]))]

    half_y = stamp.shape[0] // 2
    for hy, hx in zip(ys, xs, strict=True):
        for sy, runs in row_runs.items():
            cy = int(hy) - (sy - half_y)
            if cy < top or cy >= bottom:
                continue
            diff_row = diff[cy - top]
            for sx0, sx1 in runs:
                # Correlation placement: p - c = stamp_index - stamp_center.
                x0 = max(left, int(hx) - (sx1 - half_x))
                x1 = min(right, int(hx) - (sx0 - half_x) + 1)
                if x0 < x1:
                    diff_row[x0 - left] += 1
                    diff_row[x1 - left] -= 1

    covered = np.cumsum(diff[:, :width], axis=1, dtype=np.int32)
    return slice(top, bottom), slice(left, right), 2 * covered >= len(ys)


def simulate_cno(
    link: ContactLink | None,
    unit_id: str,
    grid: GridSpec,
    diameter_um: float,
    contact_margin_um: float,
) -> CnoResult:
    failure = np.zeros(grid.shape_yx, dtype=bool)
    if link is None or not link.holes:
        return CnoResult(unit_id, failure, CONTACT_LOST)
    # Keep the argument for callers; CNO deliberately ignores the contact/open margin.
    stamp = particle_footprint(grid, diameter_um / 2.0, 0.0)
    for hole in link.holes:
        rows, cols, blocked = hole_blocked_centres(
            hole.cut_ys, hole.cut_xs, grid.shape_yx, stamp
        )
        failure[rows, cols] |= blocked
    return CnoResult(unit_id, failure, CNO_CONFIRMED)
