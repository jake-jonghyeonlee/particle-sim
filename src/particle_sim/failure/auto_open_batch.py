"""Exact footprint-envelope certificates for groups of auto-open centres.

A 9×9 block contains every candidate at its grid locations, while the two
envelopes include all 81 possible centres. This can only make a certificate
more conservative. The caller supplies the unchanged terminal group rule.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
from scipy.ndimage import maximum_filter, minimum_filter

from particle_sim.failure.open_fast import LocalConnectivity

_Judge = Callable[[list[set[tuple[str, int]]] | None, list[list[int]], tuple[str, ...]], bool]
_BLOCK_RADIUS = 4
_BLOCK_WIDTH = 2 * _BLOCK_RADIUS + 1


@dataclass
class BatchStats:
    certified_false: int = 0
    certified_true: int = 0
    fallback_centres: int = 0
    upper_evaluations: int = 0
    lower_evaluations: int = 0
    blocks: int = 0


def footprint_envelopes(stamp: np.ndarray, radius: int = _BLOCK_RADIUS) -> tuple[np.ndarray, np.ndarray]:
    """OR and AND of the actual Boolean stamp at every centre in a block."""
    if stamp.ndim != 2 or not all(size % 2 for size in stamp.shape):
        raise ValueError("particle footprint must have odd two-dimensional extent")
    padded = np.pad(np.asarray(stamp, dtype=bool), radius, constant_values=False)
    size = 2 * radius + 1
    return (
        maximum_filter(padded, size=size, mode="constant", cval=False),
        minimum_filter(padded, size=size, mode="constant", cval=False),
    )


def evaluate_blocks(
    crop: np.ndarray,
    terminal_regions: Sequence[np.ndarray],
    groups: list[list[int]],
    kinds: tuple[str, ...],
    stamp: np.ndarray,
    ys: np.ndarray,
    xs: np.ndarray,
    y0: int,
    x0: int,
    failure: np.ndarray,
    judge: _Judge,
) -> BatchStats:
    """Write original candidate pixels only, falling back for uncertified blocks."""
    stats = BatchStats()
    if not ys.size:
        return stats
    upper, lower = footprint_envelopes(stamp)
    half_y, half_x = stamp.shape[0] // 2, stamp.shape[1] // 2
    lc = LocalConnectivity(crop, terminal_regions, half_y + _BLOCK_RADIUS,
                           half_x + _BLOCK_RADIUS)
    cy = ys.astype(np.int64, copy=False) - y0
    cx = xs.astype(np.int64, copy=False) - x0
    by = cy // _BLOCK_WIDTH
    bx = cx // _BLOCK_WIDTH
    anchor_y = by * _BLOCK_WIDTH + _BLOCK_RADIUS
    anchor_x = bx * _BLOCK_WIDTH + _BLOCK_RADIUS
    order = np.lexsort((bx, by, anchor_x // lc.tile_x, anchor_y // lc.tile_y))
    sorted_by, sorted_bx = by[order], bx[order]
    starts = np.flatnonzero(np.r_[True,
        (sorted_by[1:] != sorted_by[:-1]) | (sorted_bx[1:] != sorted_bx[:-1])])
    participating = {index for members in groups for index in members}
    for run, start in enumerate(starts):
        stop = int(starts[run + 1]) if run + 1 < len(starts) else len(order)
        members = order[start:stop]
        center = (int(sorted_by[start]) * _BLOCK_WIDTH + _BLOCK_RADIUS,
                  int(sorted_bx[start]) * _BLOCK_WIDTH + _BLOCK_RADIUS)
        stats.blocks += 1
        upper_classes = lc.terminal_classes(center, upper)
        stats.upper_evaluations += 1
        certified: bool | None = None
        if upper_classes is None:
            certified = False
        else:
            upper_open = judge(upper_classes, groups, kinds)
            if not upper_open and all(upper_classes[i] for i in participating):
                certified = False
            else:
                lower_classes = lc.terminal_classes(center, lower)
                stats.lower_evaluations += 1
                if lower_classes is not None and all(
                    bool(upper_classes[i]) == bool(lower_classes[i]) for i in participating
                ):
                    lower_open = judge(lower_classes, groups, kinds)
                    if upper_open == lower_open:
                        certified = upper_open
        if certified is None:
            stats.fallback_centres += len(members)
            for index in members:
                i = int(index)
                failure[int(ys[i]), int(xs[i])] = judge(
                    lc.terminal_classes((int(cy[i]), int(cx[i])), stamp), groups, kinds
                )
        else:
            if certified:
                stats.certified_true += len(members)
            else:
                stats.certified_false += len(members)
            failure[ys[members], xs[members]] = certified
    return stats
