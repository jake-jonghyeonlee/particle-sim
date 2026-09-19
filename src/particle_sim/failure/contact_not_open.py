"""CNO (Contact Not Open): a particle blocks a contact·via hole between two electrodes.

Provisional criterion (TODO_PHYSICS, see docs/CODE-CHANGE-BACKLOG.md 11): a hole fails only
when the particle covers every effective contact pixel, and a link fails when any one of its
holes fails. Partial coverage and redundant parallel holes are open questions.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import correlate

from particle_sim.failure.open_failure import particle_footprint
from particle_sim.geometry.contacts import ContactLink
from particle_sim.io.masks import GridSpec

CNO_MODE = "contact_not_open"
CNO_CONFIRMED = "cno_full_cover"
CONTACT_LOST = "contact_lost"
CNO_COVERAGE_RULE = "full_cover_any_hole"


@dataclass(frozen=True)
class CnoResult:
    unit_id: str
    failure: np.ndarray
    classification: str


def hole_blocked_centres(
    ys: np.ndarray, xs: np.ndarray, shape_yx: tuple[int, int], stamp: np.ndarray
) -> tuple[slice, slice, np.ndarray]:
    """Centres whose particle disc covers every given hole pixel, in a local window."""
    half_y, half_x = stamp.shape[0] // 2, stamp.shape[1] // 2
    top, bottom = max(0, int(ys.min()) - half_y), min(shape_yx[0], int(ys.max()) + half_y + 1)
    left, right = max(0, int(xs.min()) - half_x), min(shape_yx[1], int(xs.max()) + half_x + 1)
    hole = np.zeros((bottom - top, right - left), dtype=np.int32)
    hole[ys - top, xs - left] = 1
    # covered(c) = sum_p hole[p] * disc[p - c]; the disc is symmetric, so this is a correlation
    covered = correlate(hole, stamp.astype(np.int32), mode="constant", cval=0)
    return slice(top, bottom), slice(left, right), covered == int(hole.sum())


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
    stamp = particle_footprint(grid, diameter_um / 2.0, contact_margin_um)
    for hole in link.holes:
        rows, cols, blocked = hole_blocked_centres(hole.ys, hole.xs, grid.shape_yx, stamp)
        failure[rows, cols] |= blocked
    return CnoResult(unit_id, failure, CNO_CONFIRMED)
