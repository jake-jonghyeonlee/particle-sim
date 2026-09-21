"""Window-local open judgement must answer exactly like the full labelling.

Every test walks all candidate centres of a shape and compares the two implementations pixel
by pixel. The reference is `causes_open_regions`, which labels the whole electrode each time.
"""

from __future__ import annotations

import numpy as np
import pytest

from particle_sim.failure.open_failure import causes_open_regions
from particle_sim.failure.open_fast import LocalConnectivity


def disc(reach_y: int, reach_x: int) -> np.ndarray:
    yy, xx = np.mgrid[-reach_y : reach_y + 1, -reach_x : reach_x + 1]
    return np.hypot(yy, xx) <= max(reach_y, reach_x)


def stamp(shape: tuple[int, int], footprint: np.ndarray, y: int, x: int) -> np.ndarray:
    half_y, half_x = footprint.shape[0] // 2, footprint.shape[1] // 2
    blocked = np.zeros(shape, dtype=bool)
    top, bottom = max(y - half_y, 0), min(y + half_y + 1, shape[0])
    left, right = max(x - half_x, 0), min(x + half_x + 1, shape[1])
    if top >= bottom or left >= right:
        return blocked
    blocked[top:bottom, left:right] = footprint[
        top - (y - half_y) : footprint.shape[0] - ((y + half_y + 1) - bottom),
        left - (x - half_x) : footprint.shape[1] - ((x + half_x + 1) - right),
    ]
    return blocked


def compare_all_centres(
    electrode: np.ndarray,
    terminal_a: np.ndarray,
    terminal_b: np.ndarray,
    reach: int,
    *,
    max_tiles: int = 32,
) -> tuple[int, int]:
    """Return (candidates, open count) after checking both implementations agree everywhere."""
    footprint = disc(reach, reach)
    fast = LocalConnectivity(
        electrode, [terminal_a, terminal_b], reach, reach, max_tiles=max_tiles
    )
    opened = 0
    checked = 0
    for y in range(electrode.shape[0]):
        for x in range(electrode.shape[1]):
            blocked = stamp(electrode.shape, footprint, y, x)
            if not (blocked & electrode).any():
                continue  # a disc that removes nothing cannot change connectivity
            expected = causes_open_regions(electrode, blocked, terminal_a, terminal_b)
            actual = fast.splits((y, x), footprint)
            assert actual == expected, f"({y},{x})에서 판정이 다릅니다: {actual} != {expected}"
            checked += 1
            opened += int(expected)
    return checked, opened


def ends(shape: tuple[int, int], a: tuple[slice, slice], b: tuple[slice, slice]):
    terminal_a = np.zeros(shape, dtype=bool)
    terminal_b = np.zeros(shape, dtype=bool)
    terminal_a[a] = True
    terminal_b[b] = True
    return terminal_a, terminal_b


def test_straight_wire_matches_reference() -> None:
    electrode = np.zeros((24, 40), dtype=bool)
    electrode[10:14, 2:38] = True
    terminal_a, terminal_b = ends(electrode.shape, (slice(10, 14), slice(2, 5)),
                                  (slice(10, 14), slice(35, 38)))
    checked, opened = compare_all_centres(electrode, terminal_a, terminal_b, reach=3)
    assert checked > 0 and opened > 0


def test_corner_wire_matches_reference() -> None:
    """ㄱ자. The first implementation lost candidates exactly on the corner."""
    electrode = np.zeros((30, 30), dtype=bool)
    electrode[4:8, 4:26] = True
    electrode[4:26, 22:26] = True
    terminal_a, terminal_b = ends(electrode.shape, (slice(4, 8), slice(4, 7)),
                                  (slice(22, 26), slice(22, 26)))
    checked, opened = compare_all_centres(electrode, terminal_a, terminal_b, reach=3)
    assert checked > 0 and opened > 0


def test_branch_is_ignored_when_it_carries_no_terminal() -> None:
    electrode = np.zeros((30, 40), dtype=bool)
    electrode[12:16, 2:38] = True
    electrode[2:12, 18:22] = True  # stub with no terminal on it
    terminal_a, terminal_b = ends(electrode.shape, (slice(12, 16), slice(2, 5)),
                                  (slice(12, 16), slice(35, 38)))
    checked, opened = compare_all_centres(electrode, terminal_a, terminal_b, reach=3)
    assert checked > 0 and opened > 0


def test_detour_keeps_the_path_connected() -> None:
    """Two parallel runs: cutting one never opens the path, so the window must see the other."""
    electrode = np.zeros((28, 40), dtype=bool)
    electrode[6:9, 4:36] = True
    electrode[19:22, 4:36] = True
    electrode[6:22, 4:7] = True
    electrode[6:22, 33:36] = True
    terminal_a, terminal_b = ends(electrode.shape, (slice(6, 9), slice(4, 7)),
                                  (slice(19, 22), slice(33, 36)))
    checked, opened = compare_all_centres(electrode, terminal_a, terminal_b, reach=2)
    assert checked > 0


def test_separate_islands_are_open_everywhere() -> None:
    """The electrode is already in two pieces: the reference says open, so must the fast path."""
    electrode = np.zeros((20, 30), dtype=bool)
    electrode[8:12, 2:12] = True
    electrode[8:12, 18:28] = True
    terminal_a, terminal_b = ends(electrode.shape, (slice(8, 12), slice(2, 5)),
                                  (slice(8, 12), slice(25, 28)))
    checked, opened = compare_all_centres(electrode, terminal_a, terminal_b, reach=2)
    assert checked == opened > 0


def test_terminal_covered_by_the_particle_is_open() -> None:
    electrode = np.zeros((16, 30), dtype=bool)
    electrode[7:10, 2:28] = True
    terminal_a, terminal_b = ends(electrode.shape, (slice(7, 10), slice(2, 4)),
                                  (slice(7, 10), slice(26, 28)))
    checked, opened = compare_all_centres(electrode, terminal_a, terminal_b, reach=4)
    assert checked > 0 and opened > 0


def test_anisotropic_reach_matches_reference() -> None:
    electrode = np.zeros((26, 40), dtype=bool)
    electrode[11:15, 2:38] = True
    terminal_a, terminal_b = ends(electrode.shape, (slice(11, 15), slice(2, 5)),
                                  (slice(11, 15), slice(35, 38)))
    footprint = disc(2, 5)
    fast = LocalConnectivity(electrode, [terminal_a, terminal_b], 2, 5)
    for y in range(electrode.shape[0]):
        for x in range(electrode.shape[1]):
            blocked = stamp(electrode.shape, footprint, y, x)
            if not (blocked & electrode).any():
                continue
            expected = causes_open_regions(electrode, blocked, terminal_a, terminal_b)
            assert fast.splits((y, x), footprint) == expected, f"({y},{x})"


def test_tile_eviction_does_not_change_the_answer() -> None:
    """A one-tile cache forces a rebuild on nearly every centre; answers must not move."""
    electrode = np.zeros((30, 40), dtype=bool)
    electrode[13:17, 2:38] = True
    terminal_a, terminal_b = ends(electrode.shape, (slice(13, 17), slice(2, 5)),
                                  (slice(13, 17), slice(35, 38)))
    checked, opened = compare_all_centres(
        electrode, terminal_a, terminal_b, reach=3, max_tiles=1
    )
    assert checked > 0 and opened > 0


def test_footprint_larger_than_reach_is_refused() -> None:
    electrode = np.zeros((12, 20), dtype=bool)
    electrode[5:8, 2:18] = True
    terminal_a, terminal_b = ends(electrode.shape, (slice(5, 8), slice(2, 4)),
                                  (slice(5, 8), slice(16, 18)))
    fast = LocalConnectivity(electrode, [terminal_a, terminal_b], 1, 1)
    with pytest.raises(ValueError, match="footprint"):
        fast.splits((6, 10), disc(4, 4))


def test_fewer_than_two_terminals_is_refused() -> None:
    electrode = np.zeros((10, 10), dtype=bool)
    electrode[4:6, 1:9] = True
    with pytest.raises(ValueError, match="단자 영역"):
        LocalConnectivity(electrode, [electrode], 1, 1)
