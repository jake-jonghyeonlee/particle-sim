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


# ---------------------------------------------------------------------------
# 1. 단자 영역이 서로 다른 연결 성분에 걸친 전극
# ---------------------------------------------------------------------------
def test_terminal_spanning_multiple_components() -> None:
    """One terminal has pixels on two separate connected components of the electrode.

    The electrode is an H-shape: two vertical bars joined by a thin bridge.
    Terminal A spans both bars (left + right), so it sits on two components
    when the bridge is cut.  This exercises the shared-intersection logic
    where transitivity cannot be assumed across components.
    """
    electrode = np.zeros((30, 40), dtype=bool)
    # left vertical bar
    electrode[4:26, 4:8] = True
    # right vertical bar
    electrode[4:26, 30:34] = True
    # thin bridge connecting them
    electrode[13:15, 8:30] = True

    # terminal A spans both bars (left bar + a piece on the right bar)
    terminal_a = np.zeros((30, 40), dtype=bool)
    terminal_a[4:8, 4:8] = True    # on the left bar
    terminal_a[4:8, 30:34] = True   # also on the right bar

    # terminal B is only on the right bar, at the bottom
    terminal_b = np.zeros((30, 40), dtype=bool)
    terminal_b[22:26, 30:34] = True

    checked, opened = compare_all_centres(electrode, terminal_a, terminal_b, reach=3)
    assert checked > 0 and opened > 0


# ---------------------------------------------------------------------------
# 2. 단자 3개 이상을 넘긴 경우
# ---------------------------------------------------------------------------
def test_three_terminals() -> None:
    """LocalConnectivity accepts ≥2 terminal regions; test with 3.

    The reference for 3 terminals is not causes_open_regions (which takes exactly 2).
    Instead we build a direct reference: all terminals must share a common component
    in the damaged electrode.
    """
    from scipy.ndimage import label as sp_label

    electrode = np.zeros((30, 50), dtype=bool)
    # T-shape: horizontal bar + vertical stub
    electrode[12:16, 2:48] = True
    electrode[2:12, 23:27] = True

    terminal_a = np.zeros((30, 50), dtype=bool)
    terminal_a[12:16, 2:5] = True
    terminal_b = np.zeros((30, 50), dtype=bool)
    terminal_b[12:16, 45:48] = True
    terminal_c = np.zeros((30, 50), dtype=bool)
    terminal_c[2:5, 23:27] = True

    def reference_three(elec: np.ndarray, blk: np.ndarray, *terms: np.ndarray) -> bool:
        """True when any terminal loses all pixels or terminals don't share a component."""
        damaged = elec & ~blk
        labels, _ = sp_label(damaged)
        shared: set[int] | None = None
        for region in terms:
            alive = damaged & region
            if not alive.any():
                return True
            classes = set(np.unique(labels[alive]).tolist()) - {0}
            shared = classes if shared is None else (shared & classes)
            if not shared:
                return True
        return False

    footprint = disc(3, 3)
    fast = LocalConnectivity(
        electrode, [terminal_a, terminal_b, terminal_c], 3, 3,
    )
    checked = 0
    for y in range(electrode.shape[0]):
        for x in range(electrode.shape[1]):
            blocked = stamp(electrode.shape, footprint, y, x)
            if not (blocked & electrode).any():
                continue
            expected = reference_three(electrode, blocked, terminal_a, terminal_b, terminal_c)
            actual = fast.splits((y, x), footprint)
            assert actual == expected, f"({y},{x})에서 판정이 다릅니다: {actual} != {expected}"
            checked += 1
    assert checked > 0


# ---------------------------------------------------------------------------
# 3. 단자 영역이 전극 밖으로 일부 삐져나온 경우
# ---------------------------------------------------------------------------
def test_terminal_partially_outside_electrode() -> None:
    """Terminal regions extend beyond the electrode mask.

    LocalConnectivity intersects terminal_regions with the electrode (region & component).
    The reference causes_open_regions also uses (damaged & terminal), so both must agree.
    """
    electrode = np.zeros((24, 40), dtype=bool)
    electrode[10:14, 4:36] = True

    # terminal A overflows above the electrode (rows 7-13 vs electrode rows 10-13)
    terminal_a = np.zeros((24, 40), dtype=bool)
    terminal_a[7:14, 4:8] = True

    # terminal B overflows below the electrode (rows 10-17 vs electrode rows 10-13)
    terminal_b = np.zeros((24, 40), dtype=bool)
    terminal_b[10:17, 33:38] = True

    checked, opened = compare_all_centres(electrode, terminal_a, terminal_b, reach=3)
    assert checked > 0 and opened > 0


# ---------------------------------------------------------------------------
# 4. 전극이 이미지 가장자리에 닿아 창이 max(0)/min(height)로 잘리는 경우
# ---------------------------------------------------------------------------
def test_electrode_touching_image_edge() -> None:
    """Electrode pixels sit on the very edge of the array.

    The window must clamp to max(0) / min(height) correctly.
    """
    electrode = np.zeros((20, 36), dtype=bool)
    # horizontal wire along the top edge (row 0)
    electrode[0:3, 0:36] = True

    terminal_a, terminal_b = ends(electrode.shape, (slice(0, 3), slice(0, 3)),
                                  (slice(0, 3), slice(33, 36)))
    checked, opened = compare_all_centres(electrode, terminal_a, terminal_b, reach=3)
    assert checked > 0 and opened > 0


# ---------------------------------------------------------------------------
# 5. 디스크가 전극과 전혀 겹치지 않아 early return False를 타는 경로
# ---------------------------------------------------------------------------
def test_disc_no_overlap_returns_false() -> None:
    """A disc placed entirely outside the electrode must return False (no open).

    compare_all_centres skips such centres (blocked & electrode is empty), so we
    call LocalConnectivity.splits directly.
    """
    electrode = np.zeros((24, 40), dtype=bool)
    electrode[10:14, 15:25] = True

    terminal_a, terminal_b = ends(electrode.shape, (slice(10, 14), slice(15, 18)),
                                  (slice(10, 14), slice(22, 25)))
    footprint = disc(2, 2)
    fast = LocalConnectivity(electrode, [terminal_a, terminal_b], 2, 2)

    # Pick centres far from the electrode where the disc does not touch it at all
    checked = 0
    for y, x in [(0, 0), (0, 20), (23, 0), (23, 39), (20, 5)]:
        blocked = stamp(electrode.shape, footprint, y, x)
        if (blocked & electrode).any():
            continue  # not an interesting centre for this test
        result = fast.splits((y, x), footprint)
        assert result is False, f"({y},{x})에서 겹침 없는데 open: {result}"
        checked += 1
    assert checked > 0, "겹침 없는 후보 중심이 하나도 없었습니다"
