"""Window-local open judgement must answer exactly like the full labelling.

Every test walks all candidate centres of a shape and compares the two implementations pixel
by pixel. The reference is `causes_open_regions`, which labels the whole electrode each time.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.ndimage import label

from particle_sim.failure.open_failure import causes_open_regions
from particle_sim.failure.open_fast import LocalConnectivity, _UnionFind


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


def _original_tile_terminals(fast: LocalConnectivity, center: tuple[int, int]):
    """Original whole-crop terminal preparation, independent of the new slice path."""
    tile = fast._tile_for(*center)
    outside = fast.component.copy()
    outside[tile.top:tile.bottom, tile.left:tile.right] = False
    labels, _ = label(outside)
    expected = []
    for region in fast.terminals:
        inside = np.zeros_like(region)
        inside[tile.top:tile.bottom, tile.left:tile.right] = region[
            tile.top:tile.bottom, tile.left:tile.right
        ]
        ys, xs = np.nonzero(inside)
        beyond = region & ~inside
        outer = frozenset(int(value) for value in np.unique(labels[beyond]) if value)
        expected.append(((ys - tile.top, xs - tile.left), outer))
    return tile, expected


def _original_terminal_classes(
    fast: LocalConnectivity, center: tuple[int, int], footprint: np.ndarray
):
    """Original four-strip merge loop, including every repeated border pair."""
    tile = fast._tile_for(*center)
    blocked = stamp(fast.component.shape, footprint, *center)[
        tile.top:tile.bottom, tile.left:tile.right
    ]
    damaged = tile.inside_component & ~blocked
    if np.array_equal(damaged, tile.inside_component):
        return None
    inside_labels, _ = label(damaged)
    union = _UnionFind()
    height, width = damaged.shape
    for inside, outside in (
        (inside_labels[0], tile.outside_above),
        (inside_labels[height - 1], tile.outside_below),
        (inside_labels[:, 0], tile.outside_left),
        (inside_labels[:, width - 1], tile.outside_right),
    ):
        for index in np.nonzero((inside > 0) & (outside > 0))[0]:
            union.union(("in", int(inside[index])), ("out", int(outside[index])))
    result = []
    for (ys, xs), outer in zip(tile.terminal_local, tile.terminal_outside, strict=True):
        classes = {union.find(("out", value)) for value in outer}
        if ys.size:
            alive = inside_labels[ys, xs]
            classes.update(union.find(("in", int(value))) for value in np.unique(alive) if value)
        result.append(classes)
    return result


@pytest.mark.parametrize("seed", [0, 1, 2, 3])
def test_tile_preparation_and_classes_match_original(seed: int) -> None:
    """Compare every centre, including diagonal contacts and all clipped ROI sides."""
    rng = np.random.default_rng(seed)
    electrode = rng.random((20, 24)) < 0.64
    electrode[9:11, :] = True
    electrode[:, 11:13] = True
    a = np.zeros_like(electrode)
    a[0:5, 0:6] = True
    a[14:20, 18:24] = True  # disconnected pieces of one terminal
    b = np.zeros_like(electrode)
    b[2:18, 4:19] = True  # overlaps A and many components
    c = np.zeros_like(electrode)
    c[0:2, 20:24] = True
    c[18:20, 0:4] = True
    fast = LocalConnectivity(electrode, [a, b, c], 2, 3, max_tiles=2)
    footprint = disc(2, 3)
    for y in range(electrode.shape[0]):
        for x in range(electrode.shape[1]):
            tile, reference = _original_tile_terminals(fast, (y, x))
            for (ys, xs), outer, (expected_local, expected_outer) in zip(
                tile.terminal_local, tile.terminal_outside, reference, strict=True
            ):
                np.testing.assert_array_equal(ys, expected_local[0])
                np.testing.assert_array_equal(xs, expected_local[1])
                assert outer == expected_outer
            assert fast.terminal_classes((y, x), footprint) == _original_terminal_classes(
                fast, (y, x), footprint
            )


def test_border_reduction_keeps_separated_repeated_pairs(monkeypatch) -> None:
    rng = np.random.default_rng(4)
    electrode = rng.random((20, 24)) < 0.64
    electrode[9:11, :] = True
    electrode[:, 11:13] = True
    fast = LocalConnectivity(electrode, [electrode, electrode], 2, 3)
    observed = []
    original_union = _UnionFind.union

    def record_union(self, left, right):
        observed.append((left, right))
        original_union(self, left, right)

    monkeypatch.setattr(_UnionFind, "union", record_union)
    fast.terminal_classes((0, 6), np.ones((5, 7), dtype=bool))
    pair = (("in", 1), ("out", 2))
    positions = [i for i, value in enumerate(observed) if value == pair]
    assert len(positions) >= 2
    assert all(right > left + 1 for left, right in zip(positions[:-1], positions[1:], strict=True))
    assert (("in", 1), ("out", 7)) in observed
    assert (("in", 1), ("out", 8)) in observed


def test_diagonal_contact_across_tile_border_stays_separate() -> None:
    electrode = np.zeros((10, 10), dtype=bool)
    electrode[2, 2] = True  # damage makes terminal_classes take the label path
    electrode[4, 4] = True  # inside the centre's tile
    electrode[5, 5] = True  # diagonal outside the tile, no 4-connected edge
    a = np.zeros_like(electrode)
    b = np.zeros_like(electrode)
    a[4, 4] = True
    b[5, 5] = True
    fast = LocalConnectivity(electrode, [a, b], 1, 1)
    footprint = np.ones((3, 3), dtype=bool)
    tile, reference = _original_tile_terminals(fast, (2, 2))
    assert (tile.bottom, tile.right) == (5, 5)
    for (ys, xs), outer, (expected_local, expected_outer) in zip(
        tile.terminal_local, tile.terminal_outside, reference, strict=True
    ):
        np.testing.assert_array_equal(ys, expected_local[0])
        np.testing.assert_array_equal(xs, expected_local[1])
        assert outer == expected_outer
    classes = fast.terminal_classes((2, 2), footprint)
    assert classes == _original_terminal_classes(fast, (2, 2), footprint)
    assert classes is not None and classes[0].isdisjoint(classes[1])
