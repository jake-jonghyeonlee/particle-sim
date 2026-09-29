"""Window-local connectivity for open judgement.

The reference judgement in `open_failure.causes_open_regions` labels the whole electrode for
every candidate centre. That is fine on a small crop and hopeless on a long wire: a 1000x1000
electrode costs minutes per particle size, because the electrode bounding box is the whole ROI
and the crop buys nothing.

A particle only removes pixels under its own disc, so connectivity outside a window around it
never changes. This module labels the window alone and merges the result with the components
outside it, which it pre-computes once per tile and reuses for every centre in that tile. The
answer is identical to the full labelling; `tests/unit/test_open_fast.py` checks that on every
candidate of several shapes.

The outside label image is not kept per tile. It costs 4 MB on a 1000x1000 mask and exhausted
memory in the first implementation. A tile stores the four border strips of the window and, for
each terminal region, the outside labels it touches — all the merge needs.

Adapted from `open_fast.py` on ops/offline-packaging (f0e5cca), which handled two point
terminals on a single component. Terminals here are regions and the electrode may hold several
components, which is what the schema 2.0 open paths and contact-derived terminals need.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from scipy.ndimage import label

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
    """A window plus only what the outside connectivity contributes to it."""

    top: int
    bottom: int
    left: int
    right: int
    inside_component: np.ndarray
    outside_above: np.ndarray
    outside_below: np.ndarray
    outside_left: np.ndarray
    outside_right: np.ndarray
    # per terminal: local (y, x) of its pixels inside the window, and the outside labels its
    # pixels outside the window sit on
    terminal_local: tuple[tuple[np.ndarray, np.ndarray], ...]
    terminal_outside: tuple[frozenset[int], ...]


class LocalConnectivity:
    """Answers the open question for many centres of one electrode and one particle size.

    `component` is the electrode mask and may hold several connected components. Each terminal
    region is intersected with it, so only electrode pixels count, exactly as the reference
    judgement does.
    """

    def __init__(
        self,
        component: np.ndarray,
        terminal_regions: Sequence[np.ndarray],
        reach_y: int,
        reach_x: int,
        *,
        max_tiles: int = 32,
    ) -> None:
        if len(terminal_regions) < 2:
            raise ValueError("단자 영역이 둘 이상이어야 open을 판정할 수 있습니다")
        self.component = component
        self.terminals = tuple(np.asarray(region) & component for region in terminal_regions)
        self.reach_y = max(int(reach_y), 0)
        self.reach_x = max(int(reach_x), 0)
        # a tile of roughly the disc size keeps the window at ~3x the disc in each axis
        self.tile_y = max(self.reach_y, 1) * 2
        self.tile_x = max(self.reach_x, 1) * 2
        self.max_tiles = max(int(max_tiles), 1)
        self._tiles: dict[tuple[int, int], _Tile] = {}
        self._order: list[tuple[int, int]] = []

    def tile_key(self, y: int, x: int) -> tuple[int, int]:
        """Tile a centre belongs to. Iterating candidates in this order avoids rebuilding."""
        return (y // self.tile_y, x // self.tile_x)

    def _tile_for(self, y: int, x: int) -> _Tile:
        key = self.tile_key(y, x)
        tile = self._tiles.get(key)
        if tile is not None:
            return tile

        height, width = self.component.shape
        top = max(key[0] * self.tile_y - self.reach_y, 0)
        bottom = min(key[0] * self.tile_y + self.tile_y + self.reach_y, height)
        left = max(key[1] * self.tile_x - self.reach_x, 0)
        right = min(key[1] * self.tile_x + self.tile_x + self.reach_x, width)

        outside = self.component.copy()
        outside[top:bottom, left:right] = False
        outside_labels, _ = label(outside)
        empty_row = np.zeros(right - left, dtype=outside_labels.dtype)
        empty_col = np.zeros(bottom - top, dtype=outside_labels.dtype)

        local: list[tuple[np.ndarray, np.ndarray]] = []
        outer: list[frozenset[int]] = []
        for region in self.terminals:
            local.append(np.nonzero(region[top:bottom, left:right]))
            # The outside image is zero in the tile window, so indexing by the
            # full terminal region gives exactly the labels beyond the window.
            labels = outside_labels[region]
            outer.append(frozenset(int(v) for v in np.unique(labels) if v))

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
            terminal_local=tuple(local),
            terminal_outside=tuple(outer),
        )
        if len(self._tiles) >= self.max_tiles:
            self._tiles.pop(self._order.pop(0), None)
        self._tiles[key] = tile
        self._order.append(key)
        return tile

    def terminal_classes(
        self, center: tuple[int, int], footprint: np.ndarray
    ) -> list[set[_Class]] | None:
        """Per-terminal union-find root class sets after placing a particle.

        Returns ``None`` when the particle does not remove any electrode pixel inside
        the window (the early-return path that ``splits`` also uses).  Otherwise returns
        a list with one ``set`` per terminal.  An empty set means the terminal lost every
        electrode pixel it had.
        """
        y, x = center
        half_y, half_x = footprint.shape[0] // 2, footprint.shape[1] // 2
        if half_y > self.reach_y or half_x > self.reach_x:
            raise ValueError("footprint가 타일 여유보다 큽니다. reach를 크게 잡아야 합니다")
        tile = self._tile_for(y, x)
        height, width = tile.bottom - tile.top, tile.right - tile.left

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
            return None  # nothing was removed inside the window
        inside_labels, _ = label(damaged)

        union = _UnionFind()
        # the window border ties the local components to the outside components
        for inside, outside in (
            (inside_labels[0], tile.outside_above),
            (inside_labels[height - 1], tile.outside_below),
            (inside_labels[:, 0], tile.outside_left),
            (inside_labels[:, width - 1], tile.outside_right),
        ):
            valid = (inside > 0) & (outside > 0)
            valid[1:] &= (inside[1:] != inside[:-1]) | (outside[1:] != outside[:-1])
            for index in np.nonzero(valid)[0]:
                union.union(("in", int(inside[index])), ("out", int(outside[index])))

        result: list[set[_Class]] = []
        for (ys, xs), outer in zip(tile.terminal_local, tile.terminal_outside, strict=True):
            classes = {union.find(("out", value)) for value in outer}
            if ys.size:
                alive = inside_labels[ys, xs]
                classes.update(union.find(("in", int(v))) for v in np.unique(alive) if v)
            result.append(classes)
        return result

    def splits(self, center: tuple[int, int], footprint: np.ndarray) -> bool:
        """True when a particle centred at `center` disconnects the terminal regions.

        Matches `open_failure.causes_open_regions`: a terminal with no surviving electrode
        pixel counts as open, and otherwise the regions must share a component.
        """
        classes = self.terminal_classes(center, footprint)
        if classes is None:
            return False  # nothing was removed inside the window
        shared: set[_Class] | None = None
        for terminal_set in classes:
            if not terminal_set:
                return True  # this terminal lost every electrode pixel it had
            shared = terminal_set if shared is None else (shared & terminal_set)
            if not shared:
                return True
        return False
