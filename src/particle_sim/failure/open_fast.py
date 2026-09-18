"""Window-local connectivity for open judgement.

`causes_open` labels the whole image for every candidate centre, which costs ~1.7 ms on a
1000x1000 mask and does not scale to 32 particle sizes. A removed particle is local, so the
connectivity outside a window never changes: this module labels only the window and merges the
result with pre-computed outside components. The answer is identical to the full labelling,
which `tests/unit/test_open_failure.py` checks against the reference implementation.
"""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass

import numpy as np
from scipy.ndimage import label

Terminal = tuple[int, int]
_Class = tuple[str, int]


class _UnionFind:
    __slots__ = ("parent",)

    def __init__(self) -> None:
        self.parent: dict[_Class, _Class] = {}

    def find(self, item: _Class) -> _Class:
        parent = self.parent.setdefault(item, item)
        while parent != item:
            item, parent = parent, self.parent.setdefault(parent, parent)
        return item

    def union(self, left: _Class, right: _Class) -> None:
        root_left, root_right = self.find(left), self.find(right)
        if root_left != root_right:
            self.parent[root_left] = root_right


@dataclass(frozen=True)
class _Tile:
    """A window plus only what the outside connectivity contributes to it.

    The outside label image is not kept: a tile stores the four border strips and the outside
    class of each terminal, which is all the merge needs. Keeping the full label image per tile
    costs 4 MB each and exhausted memory on 1000x1000 masks.
    """

    top: int
    bottom: int
    left: int
    right: int
    inside_component: np.ndarray
    outside_above: np.ndarray
    outside_below: np.ndarray
    outside_left: np.ndarray
    outside_right: np.ndarray
    terminal_outside: tuple[int | None, int | None]


class LocalConnectivity:
    """Answers `causes_open` for many centres of one electrode and one particle size."""

    MAX_TILES = 32

    def __init__(
        self,
        component: np.ndarray,
        terminal_a: Terminal,
        terminal_b: Terminal,
        reach_px_y: int,
        reach_px_x: int,
    ) -> None:
        self.component = component
        self.terminal_a = terminal_a
        self.terminal_b = terminal_b
        self.reach_y = max(int(reach_px_y), 0)
        self.reach_x = max(int(reach_px_x), 0)
        # a tile of roughly the disc size keeps the window at ~3x the disc in each axis
        self.tile_y = max(self.reach_y, 1) * 2
        self.tile_x = max(self.reach_x, 1) * 2
        self._tiles: OrderedDict[tuple[int, int], _Tile] = OrderedDict()

    def _tile_for(self, y: int, x: int) -> _Tile:
        key = (y // self.tile_y, x // self.tile_x)
        tile = self._tiles.get(key)
        if tile is not None:
            return tile
        height, width = self.component.shape
        tile_top = key[0] * self.tile_y
        tile_left = key[1] * self.tile_x
        top = max(tile_top - self.reach_y, 0)
        bottom = min(tile_top + self.tile_y + self.reach_y, height)
        left = max(tile_left - self.reach_x, 0)
        right = min(tile_left + self.tile_x + self.reach_x, width)

        outside = self.component.copy()
        outside[top:bottom, left:right] = False
        outside_labels, _ = label(outside)
        empty_row = np.zeros(right - left, dtype=outside_labels.dtype)
        empty_col = np.zeros(bottom - top, dtype=outside_labels.dtype)
        tile = _Tile(
            top=top,
            bottom=bottom,
            left=left,
            right=right,
            inside_component=self.component[top:bottom, left:right].copy(),
            outside_above=(
                outside_labels[top - 1, left:right].copy() if top > 0 else empty_row
            ),
            outside_below=(
                outside_labels[bottom, left:right].copy() if bottom < height else empty_row
            ),
            outside_left=(
                outside_labels[top:bottom, left - 1].copy() if left > 0 else empty_col
            ),
            outside_right=(
                outside_labels[top:bottom, right].copy() if right < width else empty_col
            ),
            terminal_outside=(
                self._outside_label(outside_labels, self.terminal_a, top, bottom, left, right),
                self._outside_label(outside_labels, self.terminal_b, top, bottom, left, right),
            ),
        )
        if len(self._tiles) >= self.MAX_TILES:
            self._tiles.popitem(last=False)
        self._tiles[key] = tile
        return tile

    @staticmethod
    def _outside_label(
        outside_labels: np.ndarray,
        terminal: Terminal,
        top: int,
        bottom: int,
        left: int,
        right: int,
    ) -> int | None:
        y, x = terminal
        if top <= y < bottom and left <= x < right:
            return None
        return int(outside_labels[y, x])

    def _terminal_class(
        self,
        terminal: Terminal,
        index: int,
        tile: _Tile,
        inside_labels: np.ndarray,
        blocked: np.ndarray,
    ) -> _Class | None:
        y, x = terminal
        if tile.top <= y < tile.bottom and tile.left <= x < tile.right:
            local_y, local_x = y - tile.top, x - tile.left
            if blocked[local_y, local_x]:
                return None  # the particle covers the terminal itself
            return ("in", int(inside_labels[local_y, local_x]))
        outside = tile.terminal_outside[index]
        return ("out", int(outside if outside is not None else 0))

    def causes_open(self, center: Terminal, footprint: np.ndarray) -> bool:
        """True when a particle centred at `center` disconnects the two terminals."""
        y, x = center
        tile = self._tile_for(y, x)
        height, width = tile.bottom - tile.top, tile.right - tile.left
        half_y, half_x = footprint.shape[0] // 2, footprint.shape[1] // 2

        blocked = np.zeros((height, width), dtype=bool)
        top = max(y - half_y, tile.top)
        bottom = min(y + half_y + 1, tile.bottom)
        left = max(x - half_x, tile.left)
        right = min(x + half_x + 1, tile.right)
        if top < bottom and left < right:
            window = footprint[
                top - (y - half_y) : footprint.shape[0] - ((y + half_y + 1) - bottom),
                left - (x - half_x) : footprint.shape[1] - ((x + half_x + 1) - right),
            ]
            blocked[
                top - tile.top : bottom - tile.top, left - tile.left : right - tile.left
            ] = window

        damaged = tile.inside_component & ~blocked
        if np.array_equal(damaged, tile.inside_component):
            return False  # nothing was removed inside the window
        inside_labels, _ = label(damaged)

        union = _UnionFind()
        # the window border ties the local components to the outside components
        for inside, outside in (
            (inside_labels[0], tile.outside_above),
            (inside_labels[height - 1], tile.outside_below),
            (inside_labels[:, 0], tile.outside_left),
            (inside_labels[:, width - 1], tile.outside_right),
        ):
            for index in np.nonzero((inside > 0) & (outside > 0))[0]:
                union.union(("in", int(inside[index])), ("out", int(outside[index])))

        class_a = self._terminal_class(self.terminal_a, 0, tile, inside_labels, blocked)
        class_b = self._terminal_class(self.terminal_b, 1, tile, inside_labels, blocked)
        if class_a is None or class_b is None:
            return True
        if class_a[1] == 0 or class_b[1] == 0:
            return True
        return union.find(class_a) != union.find(class_b)
