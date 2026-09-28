"""Open judgement: a particle cuts an electrode between its terminals."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import binary_dilation, distance_transform_edt, label
from scipy.spatial import cKDTree

from particle_sim.errors import InputContractError
from particle_sim.failure.open_fast import LocalConnectivity
from particle_sim.io.masks import GridSpec, MaskSet
from particle_sim.schema import OpenPath, ProductManifest

# A particle only blocks the electrode where it actually overlaps it; the blocked disc is
# built from the XY distance map so it stays consistent with the contact model.
CONFIRMED = "open_confirmed"
PROXY = "open_proxy"
# Excluded conditions: the path is already broken, or a terminal has no electrode left, before
# any particle lands. Counting every candidate as open there would blame particles for it.
BASELINE_OPEN = "baseline_open"
TERMINAL_LOST = "terminal_lost"
EXCLUDED_CLASSIFICATIONS = (BASELINE_OPEN, TERMINAL_LOST)


@dataclass(frozen=True)
class OpenResult:
    path_id: str
    electrode: str
    failure: np.ndarray  # bool (y, x): particle centre positions that open the path
    classification: str  # CONFIRMED | PROXY | BASELINE_OPEN | TERMINAL_LOST
    terminal_a_yx: tuple[int, int]
    terminal_b_yx: tuple[int, int]
    candidates: int  # evaluated particle centres (skeleton-reduced)
    terminal_kind: str = "point"


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


def causes_open_regions(
    electrode_mask: np.ndarray,
    blocked: np.ndarray,
    terminal_a: np.ndarray,
    terminal_b: np.ndarray,
) -> bool:
    """두 단자 영역에서 살아남은 전극 픽셀 사이의 경로가 사라지면 True."""
    damaged = electrode_mask & ~blocked
    alive_a = damaged & terminal_a
    alive_b = damaged & terminal_b
    if not alive_a.any() or not alive_b.any():
        return True
    labels, _ = label(damaged)
    labels_a = set(np.unique(labels[alive_a]).tolist()) - {0}
    labels_b = set(np.unique(labels[alive_b]).tolist()) - {0}
    return labels_a.isdisjoint(labels_b)


def resolve_terminals(
    open_path: OpenPath,
    electrode_mask: np.ndarray,
    product_id: str,
    *,
    require_on_electrode: bool = True,
) -> tuple[tuple[int, int], tuple[int, int], str]:
    """Return both terminals and whether the result is confirmed or a proxy.

    A transformed scenario may move the electrode away from a fixed point terminal; the caller
    then passes ``require_on_electrode=False`` and classifies the lost terminal instead.
    """
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
            if require_on_electrode and not electrode_mask[terminal]:
                raise InputContractError(
                    f"open_path '{open_path.path_id}'의 {name} {terminal}가 전극 위에 없습니다",
                    product_id=product_id,
                    field=f"open_paths[{open_path.path_id}].{name}",
                )
        return terminals[0], terminals[1], CONFIRMED

    if open_path.fallback == "roi_boundary":
        labels, count = label(electrode_mask)
        if count == 0 and not require_on_electrode:
            return (-1, -1), (-1, -1), PROXY
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


def resolve_terminal_regions(
    open_path: OpenPath,
    masks: MaskSet,
    electrode_mask: np.ndarray,
    product_id: str,
    *,
    require_on_electrode: bool = True,
) -> tuple[np.ndarray, np.ndarray, str, str]:
    """단자 mask를 우선 사용하고 점/fallback은 1픽셀 영역으로 변환한다."""
    regions = (masks.terminal_regions or {}).get(open_path.path_id)
    if regions is not None:
        terminal_a, terminal_b = regions
        return terminal_a, terminal_b, CONFIRMED, "region_mask"
    point_a, point_b, classification = resolve_terminals(
        open_path, electrode_mask, product_id, require_on_electrode=require_on_electrode
    )
    terminal_a = np.zeros(electrode_mask.shape, dtype=bool)
    terminal_b = np.zeros(electrode_mask.shape, dtype=bool)
    for region, point in ((terminal_a, point_a), (terminal_b, point_b)):
        if point[0] >= 0:
            region[point] = True
    return terminal_a, terminal_b, classification, "point"


def baseline_status(
    electrode_mask: np.ndarray, terminal_a: np.ndarray, terminal_b: np.ndarray
) -> str | None:
    """Classify the particle-free state: None when both terminals are connected."""
    if not (electrode_mask & terminal_a).any() or not (electrode_mask & terminal_b).any():
        return TERMINAL_LOST
    empty = np.zeros(electrode_mask.shape, dtype=bool)
    if causes_open_regions(electrode_mask, empty, terminal_a, terminal_b):
        return BASELINE_OPEN
    return None


def check_open_baselines(manifest: ProductManifest, masks: MaskSet) -> None:
    """Nominal inputs must connect every open path's terminals before any particle lands."""
    for index, open_path in enumerate(manifest.open_paths):
        electrode_mask = masks.electrodes[open_path.electrode]
        terminal_a, terminal_b, _classification, _kind = resolve_terminal_regions(
            open_path, masks, electrode_mask, manifest.product_id
        )
        status = baseline_status(electrode_mask, terminal_a, terminal_b)
        if status == TERMINAL_LOST:
            raise InputContractError(
                f"open_path '{open_path.path_id}'의 단자가 전극과 겹치지 않습니다",
                product_id=manifest.product_id,
                field=f"open_paths[{index}]",
            )
        if status == BASELINE_OPEN:
            raise InputContractError(
                f"open_path '{open_path.path_id}'의 두 단자가 이물 없이도 연결돼 있지 않습니다",
                product_id=manifest.product_id,
                field=f"open_paths[{index}]",
            )


def _first_point(region: np.ndarray) -> tuple[int, int]:
    points = np.argwhere(region)
    if points.size == 0:
        return (-1, -1)
    return int(points[0, 0]), int(points[0, 1])


def _mask_bounds(mask: np.ndarray) -> tuple[int, int, int, int]:
    """Return the tight non-empty mask bounds as half-open y/x coordinates."""
    ys, xs = np.nonzero(mask)
    if ys.size == 0:
        raise ValueError("비어 있는 mask의 경계를 계산할 수 없습니다")
    return int(ys.min()), int(ys.max()) + 1, int(xs.min()), int(xs.max()) + 1


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
    if not np.isfinite(reach_um):
        _reach_footprint(grid, reach_um)
    seeds = np.asarray(electrode_mask, dtype=bool).copy()
    for terminal in terminals:
        seeds[terminal] = True
    if not seeds.any():
        return np.zeros(electrode_mask.shape, dtype=bool)

    nearest_y, nearest_x = distance_transform_edt(
        ~seeds,
        sampling=(grid.pixel_size_y_um, grid.pixel_size_x_um),
        return_distances=False,
        return_indices=True,
    )
    yy, xx = np.ogrid[: seeds.shape[0], : seeds.shape[1]]
    distance = np.hypot(
        (nearest_y - yy) * grid.pixel_size_y_um,
        (nearest_x - xx) * grid.pixel_size_x_um,
    )
    candidates = distance <= reach_um

    # The EDT chooses one of equally near seeds using squared distances. At the
    # boundary, another seed can round differently under the footprint's hypot.
    # Check only unresolved centres against nearby seeds with that exact predicate.
    tolerance = 32 * np.finfo(float).eps * max(
        1.0,
        abs(reach_um),
        seeds.shape[0] * grid.pixel_size_y_um,
        seeds.shape[1] * grid.pixel_size_x_um,
    )
    boundary = np.argwhere((~candidates) & (distance <= reach_um + tolerance))
    if boundary.size:
        seed_points = np.argwhere(seeds)
        seed_um = seed_points * (grid.pixel_size_y_um, grid.pixel_size_x_um)
        tree = cKDTree(seed_um)
        for y, x in boundary:
            nearby = tree.query_ball_point(
                (y * grid.pixel_size_y_um, x * grid.pixel_size_x_um),
                reach_um + tolerance,
            )
            offsets = seed_points[nearby] - (y, x)
            candidates[y, x] = np.any(
                np.hypot(
                    offsets[:, 0] * grid.pixel_size_y_um,
                    offsets[:, 1] * grid.pixel_size_x_um,
                ) <= reach_um
            )
    return candidates


def _open_candidates_reference(
    electrode_mask: np.ndarray,
    grid: GridSpec,
    radius_um: float,
    contact_margin_um: float,
    terminals: tuple[tuple[int, int], ...] = (),
) -> np.ndarray:
    """Original dilation, retained for exact-equivalence tests only."""
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
    terminal_a_region, terminal_b_region, classification, terminal_kind = (
        resolve_terminal_regions(
            open_path, masks, electrode_mask, manifest.product_id, require_on_electrode=False
        )
    )
    terminal_a = _first_point(terminal_a_region)
    terminal_b = _first_point(terminal_b_region)

    failure = np.zeros(electrode_mask.shape, dtype=bool)
    excluded = baseline_status(electrode_mask, terminal_a_region, terminal_b_region)
    if excluded is not None:
        return OpenResult(
            path_id=open_path.path_id,
            electrode=open_path.electrode,
            failure=failure,
            classification=excluded,
            terminal_a_yx=terminal_a,
            terminal_b_yx=terminal_b,
            candidates=0,
            terminal_kind=terminal_kind,
        )
    candidates = open_candidates(
        electrode_mask, grid, radius_um, margin
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
            terminal_kind=terminal_kind,
        )

    stamp = particle_footprint(grid, radius_um, margin)
    half_y, half_x = stamp.shape[0] // 2, stamp.shape[1] // 2
    height, width = electrode_mask.shape
    # Crop to the electrode, widened by the disc so every candidate centre still falls inside.
    # Candidates are the electrode dilated by the reach, so this bound holds by construction.
    crop_y0, crop_y1, crop_x0, crop_x1 = _mask_bounds(electrode_mask)
    crop_y0 = max(crop_y0 - half_y, 0)
    crop_y1 = min(crop_y1 + half_y, height)
    crop_x0 = max(crop_x0 - half_x, 0)
    crop_x1 = min(crop_x1 + half_x, width)
    cropped_electrode = electrode_mask[crop_y0:crop_y1, crop_x0:crop_x1]
    cropped_terminal_a = terminal_a_region[crop_y0:crop_y1, crop_x0:crop_x1]
    cropped_terminal_b = terminal_b_region[crop_y0:crop_y1, crop_x0:crop_x1]
    local = LocalConnectivity(
        cropped_electrode,
        [cropped_terminal_a, cropped_terminal_b],
        half_y,
        half_x,
    )
    ys, xs = np.nonzero(candidates)
    centres = [
        (int(y) - crop_y0, int(x) - crop_x0)
        for y, x in zip(ys.tolist(), xs.tolist(), strict=True)
    ]
    # Walking tile by tile keeps each tile's outside labelling built once instead of on every
    # row sweep, which otherwise evicts and rebuilds it across a wide electrode.
    order = sorted(range(len(centres)), key=lambda i: local.tile_key(*centres[i]))
    for i in order:
        cy, cx = centres[i]
        failure[int(ys[i]), int(xs[i])] = local.splits((cy, cx), stamp)
    return OpenResult(
        path_id=open_path.path_id,
        electrode=open_path.electrode,
        failure=failure,
        classification=classification,
        terminal_a_yx=terminal_a,
        terminal_b_yx=terminal_b,
        candidates=int(candidates.sum()),
        terminal_kind=terminal_kind,
    )
