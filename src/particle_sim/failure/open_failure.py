"""Open judgement: a particle cuts an electrode between its terminals."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import binary_dilation, label

from particle_sim.errors import InputContractError
from particle_sim.io.masks import GridSpec, MaskSet
from particle_sim.schema import OpenPath, ProductManifest

# A particle only blocks the electrode where it actually overlaps it; the blocked disc is
# built from the XY distance map so it stays consistent with the contact model.
CONFIRMED = "open_confirmed"
PROXY = "open_proxy"


@dataclass(frozen=True)
class OpenResult:
    path_id: str
    electrode: str
    failure: np.ndarray  # bool (y, x): particle centre positions that open the path
    classification: str  # CONFIRMED | PROXY
    terminal_a_yx: tuple[int, int]
    terminal_b_yx: tuple[int, int]
    candidates: int  # evaluated particle centres (skeleton-reduced)


def blocked_mask(
    distance_to_center_um: np.ndarray, radius_um: float, contact_margin_um: float
) -> np.ndarray:
    """Pixels covered by a particle of `radius_um` centred at distance 0."""
    return distance_to_center_um <= radius_um + contact_margin_um


def causes_open(
    electrode_mask: np.ndarray,
    blocked: np.ndarray,
    terminal_a_yx: tuple[int, int],
    terminal_b_yx: tuple[int, int],
) -> bool:
    """True when the two terminals stop being connected after removing `blocked`."""
    damaged = electrode_mask & ~blocked
    labels, _ = label(damaged)
    label_a = int(labels[terminal_a_yx])
    label_b = int(labels[terminal_b_yx])
    return label_a == 0 or label_b == 0 or label_a != label_b


def resolve_terminals(
    open_path: OpenPath, electrode_mask: np.ndarray, product_id: str
) -> tuple[tuple[int, int], tuple[int, int], str]:
    """Return both terminals and whether the result is confirmed or a proxy."""
    shape = electrode_mask.shape
    if open_path.terminal_a_yx_px is not None and open_path.terminal_b_yx_px is not None:
        terminals = (
            tuple(open_path.terminal_a_yx_px),
            tuple(open_path.terminal_b_yx_px),
        )
        for name, terminal in zip(("terminal_a_yx_px", "terminal_b_yx_px"), terminals,
                                  strict=True):
            if not (0 <= terminal[0] < shape[0] and 0 <= terminal[1] < shape[1]):
                raise InputContractError(
                    f"open_path '{open_path.path_id}'의 {name} {terminal}가 이미지 밖입니다",
                    product_id=product_id,
                    field=f"open_paths[{open_path.path_id}].{name}",
                )
            if not electrode_mask[terminal]:
                raise InputContractError(
                    f"open_path '{open_path.path_id}'의 {name} {terminal}가 전극 위에 없습니다",
                    product_id=product_id,
                    field=f"open_paths[{open_path.path_id}].{name}",
                )
        return terminals[0], terminals[1], CONFIRMED

    if open_path.fallback == "roi_boundary":
        labels, count = label(electrode_mask)
        if count == 0:
            raise InputContractError(
                f"open_path '{open_path.path_id}'의 전극이 비어 있습니다",
                product_id=product_id,
                field=f"open_paths[{open_path.path_id}].electrode",
            )
        sizes = np.bincount(labels.ravel())
        sizes[0] = 0
        main = int(np.argmax(sizes))
        ys, xs = np.nonzero(labels == main)
        first = int(np.argmin(ys * shape[1] + xs))
        last = int(np.argmax(ys * shape[1] + xs))
        return (int(ys[first]), int(xs[first])), (int(ys[last]), int(xs[last])), PROXY

    raise InputContractError(
        f"open_path '{open_path.path_id}'의 단자를 확인할 수 없습니다 "
        f"(fallback={open_path.fallback})",
        product_id=product_id,
        field=f"open_paths[{open_path.path_id}].fallback",
    )


def _reach_footprint(grid: GridSpec, reach_um: float) -> np.ndarray:
    """Disc of `reach_um` in pixels, rounded outwards so no candidate is lost."""
    reach_y = int(np.ceil(reach_um / grid.pixel_size_y_um))
    reach_x = int(np.ceil(reach_um / grid.pixel_size_x_um))
    yy, xx = np.mgrid[-reach_y : reach_y + 1, -reach_x : reach_x + 1]
    return np.hypot(yy * grid.pixel_size_y_um, xx * grid.pixel_size_x_um) <= reach_um


def particle_footprint(grid: GridSpec, radius_um: float, contact_margin_um: float
                      ) -> np.ndarray:
    """Disc stamp of a particle centred on the middle pixel."""
    reach_um = max(radius_um + contact_margin_um, 0.0)
    reach_y = int(np.floor(reach_um / grid.pixel_size_y_um))
    reach_x = int(np.floor(reach_um / grid.pixel_size_x_um))
    yy, xx = np.mgrid[-reach_y : reach_y + 1, -reach_x : reach_x + 1]
    distance = np.hypot(yy * grid.pixel_size_y_um, xx * grid.pixel_size_x_um)
    return blocked_mask(distance, radius_um, contact_margin_um)


def open_candidates(
    electrode_mask: np.ndarray,
    grid: GridSpec,
    radius_um: float,
    contact_margin_um: float,
    terminals: tuple[tuple[int, int], ...] = (),
) -> np.ndarray:
    """Particle centres that can possibly open the path.

    Only a disc that actually covers at least one electrode pixel can change the connectivity,
    so every centre within reach of the electrode (or of a terminal) is a candidate. This is a
    provable superset: a disc that touches nothing removes nothing. A skeleton/half-width
    filter was tried first and was dropped because it loses centres on corners and junctions,
    where `skeletonize` deletes pixels that still carry the connection.
    """
    reach_um = radius_um + contact_margin_um
    if reach_um <= 0:
        return np.zeros(electrode_mask.shape, dtype=bool)
    seeds = electrode_mask.copy()
    for terminal in terminals:
        seeds[terminal] = True
    return binary_dilation(seeds, structure=_reach_footprint(grid, reach_um))


def simulate_open(
    manifest: ProductManifest,
    masks: MaskSet,
    open_path: OpenPath,
    diameter_um: float,
    *,
    contact_margin_um: float | None = None,
) -> OpenResult:
    """Evaluate every candidate particle centre for one open path and particle size."""
    if open_path.electrode not in masks.electrodes:
        raise InputContractError(
            f"open_path '{open_path.path_id}'의 전극 '{open_path.electrode}'을(를) "
            "찾을 수 없습니다",
            product_id=manifest.product_id,
            field=f"open_paths[{open_path.path_id}].electrode",
        )
    electrode_mask = masks.electrodes[open_path.electrode]
    grid = masks.grid
    margin = manifest.contact_margin_um if contact_margin_um is None else contact_margin_um
    radius_um = diameter_um / 2.0
    terminal_a, terminal_b, classification = resolve_terminals(
        open_path, electrode_mask, manifest.product_id
    )

    failure = np.zeros(electrode_mask.shape, dtype=bool)
    candidates = open_candidates(
        electrode_mask, grid, radius_um, margin, (terminal_a, terminal_b)
    )
    if not candidates.any():
        return OpenResult(
            path_id=open_path.path_id,
            electrode=open_path.electrode,
            failure=failure,
            classification=classification,
            terminal_a_yx=terminal_a,
            terminal_b_yx=terminal_b,
            candidates=0,
        )

    stamp = particle_footprint(grid, radius_um, margin)
    half_y, half_x = stamp.shape[0] // 2, stamp.shape[1] // 2
    height, width = electrode_mask.shape
    blocked = np.zeros(electrode_mask.shape, dtype=bool)
    ys, xs = np.nonzero(candidates)
    for y, x in zip(ys.tolist(), xs.tolist(), strict=True):
        top, bottom = max(0, y - half_y), min(height, y + half_y + 1)
        left, right = max(0, x - half_x), min(width, x + half_x + 1)
        window = stamp[
            top - (y - half_y) : stamp.shape[0] - ((y + half_y + 1) - bottom),
            left - (x - half_x) : stamp.shape[1] - ((x + half_x + 1) - right),
        ]
        blocked[top:bottom, left:right] = window
        failure[y, x] = causes_open(electrode_mask, blocked, terminal_a, terminal_b)
        blocked[top:bottom, left:right] = False
    return OpenResult(
        path_id=open_path.path_id,
        electrode=open_path.electrode,
        failure=failure,
        classification=classification,
        terminal_a_yx=terminal_a,
        terminal_b_yx=terminal_b,
        candidates=int(candidates.sum()),
    )
