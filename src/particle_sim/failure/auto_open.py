"""Open judgement with terminals derived from contacts and the ROI boundary.

Each contact hole touching an electrode and each run of electrode pixels on the ROI boundary
(the wire continues outside the ROI) is a terminal. Terminals that are connected before any
particle lands form a group; a particle opens the electrode when it splits a group. A fragment
without a terminal carries no function and is ignored. A fully covered contact is a CNO, not an
open, so it drops out of the comparison; a fully covered boundary run cuts the wire and is open.
"""

from __future__ import annotations

import os
from dataclasses import dataclass

import numpy as np
from scipy.ndimage import label

from particle_sim.failure.auto_open_batch import evaluate_blocks
from particle_sim.failure.open_failure import _mask_bounds, open_candidates, particle_footprint
from particle_sim.failure.open_fast import LocalConnectivity
from particle_sim.io.masks import MaskSet

AUTO_CONFIRMED = "open_auto_terminals"
INSUFFICIENT_TERMINALS = "insufficient_terminals"


@dataclass(frozen=True)
class Terminal:
    kind: str  # "contact" | "roi_boundary"
    name: str
    ys: np.ndarray
    xs: np.ndarray


@dataclass(frozen=True)
class AutoOpenResult:
    electrode: str
    failure: np.ndarray
    classification: str
    terminals: tuple[str, ...]
    candidates: int


def electrode_terminals(masks: MaskSet, electrode_id: str) -> list[Terminal]:
    """Contact holes on the electrode, then boundary runs, in a deterministic order."""
    electrode = masks.electrodes[electrode_id]
    terminals: list[Terminal] = []
    topology = masks.contact_topology
    for hole in () if topology is None else topology.holes:
        if electrode_id in (hole.lower_electrode, hole.upper_electrode):
            terminals.append(
                Terminal(
                    kind="contact",
                    name=f"{hole.contact_id}#{hole.hole_index}",
                    ys=hole.ys,
                    xs=hole.xs,
                )
            )
    border = np.zeros(electrode.shape, dtype=bool)
    border[0, :] = border[-1, :] = True
    border[:, 0] = border[:, -1] = True
    runs, count = label(electrode & border)
    for run in range(1, count + 1):
        ys, xs = np.nonzero(runs == run)
        terminals.append(
            Terminal(kind="roi_boundary", name=f"roi_boundary#{run}", ys=ys, xs=xs)
        )
    return terminals


def _union_groups(label_sets: list[set[int]]) -> list[list[int]]:
    """Group member indices that share a connected-component label (union-find)."""
    parent = list(range(len(label_sets)))

    def find(item: int) -> int:
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    owner: dict[int, int] = {}
    for index, labels in enumerate(label_sets):
        for value in labels:
            if value in owner:
                parent[find(index)] = find(owner[value])
            else:
                owner[value] = index
    groups: dict[int, list[int]] = {}
    for index in range(len(label_sets)):
        groups.setdefault(find(index), []).append(index)
    return list(groups.values())


def _judges_open(
    classes: list[set[tuple[str, int]]] | None,
    groups: list[list[int]],
    kinds: tuple[str, ...],
) -> bool:
    """Apply the original per-group terminal rule to local component classes."""
    if classes is None:
        return False
    for members in groups:
        alive: list[set[tuple[str, int]]] = []
        for index in members:
            if classes[index]:
                alive.append(classes[index])
            elif kinds[index] == "roi_boundary":
                return True
            # A covered contact is CNO and leaves the group comparison.
        if len(alive) >= 2:
            root_to_int: dict[tuple[str, int], int] = {}
            mapped: list[set[int]] = []
            for class_set in alive:
                mapped.append({root_to_int.setdefault(root, len(root_to_int)) for root in class_set})
            if len(_union_groups(mapped)) > 1:
                return True
    return False


def _reference_auto_open_judgement(
    crop: np.ndarray,
    local: list[tuple[str, np.ndarray, np.ndarray]],
    groups: list[list[int]],
    stamp: np.ndarray,
    candidate_ys: np.ndarray,
    candidate_xs: np.ndarray,
    y0: int,
    y1: int,
    x0: int,
    x1: int,
    height: int,
    width: int,
) -> np.ndarray:
    """Reference per-candidate brute-force labelling (original algorithm).

    Returns a bool array shaped ``(height, width)`` marking open centres.
    """
    half_y, half_x = stamp.shape[0] // 2, stamp.shape[1] // 2
    failure = np.zeros((height, width), dtype=bool)
    blocked = np.zeros(crop.shape, dtype=bool)
    for y, x in zip(candidate_ys.tolist(), candidate_xs.tolist(), strict=True):
        stamp_top, stamp_left = y - half_y, x - half_x
        top, bottom = max(0, y0, stamp_top), min(height, y1, y + half_y + 1)
        left, right = max(0, x0, stamp_left), min(width, x1, x + half_x + 1)
        if top >= bottom or left >= right:
            continue
        window = stamp[top - stamp_top : bottom - stamp_top, left - stamp_left : right - stamp_left]
        blocked[top - y0 : bottom - y0, left - x0 : right - x0] = window
        damaged = crop & ~blocked
        labels, _ = label(damaged)
        is_open = False
        for members in groups:
            alive: list[set[int]] = []
            for index in members:
                kind, tys, txs = local[index]
                values = set(labels[tys, txs].tolist()) - {0}
                if values:
                    alive.append(values)
                elif kind == "roi_boundary":
                    is_open = True  # the wire is cut where it leaves the ROI
                    break
                # a fully covered contact is counted as CNO, not as open
            if is_open or (len(alive) >= 2 and len(_union_groups(alive)) > 1):
                is_open = True
                break
        failure[y, x] = is_open
        blocked[top - y0 : bottom - y0, left - x0 : right - x0] = False
    return failure


def simulate_auto_open(
    masks: MaskSet, electrode_id: str, diameter_um: float, contact_margin_um: float
) -> AutoOpenResult:
    backend = os.environ.get("PARTICLE_SIM_OPEN_BACKEND", "auto")
    if backend not in ("auto", "local"):
        raise ValueError(
            "PARTICLE_SIM_OPEN_BACKEND must be 'auto' or 'local' "
            f"(got {backend!r})"
        )
    electrode = masks.electrodes[electrode_id]
    grid = masks.grid
    failure = np.zeros(electrode.shape, dtype=bool)
    terminals = electrode_terminals(masks, electrode_id)
    names = tuple(item.name for item in terminals)
    if len(terminals) < 2 or not electrode.any():
        return AutoOpenResult(electrode_id, failure, INSUFFICIENT_TERMINALS, names, 0)

    y0, y1, x0, x1 = _mask_bounds(electrode)
    crop = electrode[y0:y1, x0:x1]
    local = [(item.kind, item.ys - y0, item.xs - x0) for item in terminals]

    # particle-free grouping: terminals already apart are never blamed on a particle
    base_labels, _ = label(crop)
    base_sets = [set(base_labels[ys, xs].tolist()) - {0} for _kind, ys, xs in local]
    groups = [members for members in _union_groups(base_sets) if len(members) >= 2]
    if not groups:
        return AutoOpenResult(electrode_id, failure, INSUFFICIENT_TERMINALS, names, 0)

    radius_um = diameter_um / 2.0
    candidates = open_candidates(electrode, grid, radius_um, contact_margin_um)
    stamp = particle_footprint(grid, radius_um, contact_margin_um)
    half_y, half_x = stamp.shape[0] // 2, stamp.shape[1] // 2
    height, width = electrode.shape

    # Build terminal region masks in the cropped coordinate system for LocalConnectivity
    terminal_regions: list[np.ndarray] = []
    for _kind, tys, txs in local:
        region = np.zeros(crop.shape, dtype=bool)
        region[tys, txs] = True
        terminal_regions.append(region)

    ys, xs = np.nonzero(candidates)
    kinds = tuple(item[0] for item in local)
    if backend == "auto":
        evaluate_blocks(crop, terminal_regions, groups, kinds, stamp, ys, xs,
                        y0, x0, failure, _judges_open)
    else:
        lc = LocalConnectivity(crop, terminal_regions, half_y, half_x)
        centres = [(int(y) - y0, int(x) - x0) for y, x in zip(ys, xs, strict=True)]
        order = sorted(range(len(centres)), key=lambda i: lc.tile_key(*centres[i]))
        for i in order:
            failure[int(ys[i]), int(xs[i])] = _judges_open(
                lc.terminal_classes(centres[i], stamp), groups, kinds
            )
    return AutoOpenResult(electrode_id, failure, AUTO_CONFIRMED, names, int(candidates.sum()))
