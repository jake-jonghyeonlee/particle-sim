"""Probe whether open judgement is monotone in particle diameter.

Premise: if a centre is open at diameter d1, it is also open at every d2 > d1.
This script tests both simulate_open (explicit 2-terminal) and
simulate_auto_open (contact-derived N-terminal) on synthetic shapes designed
to stress the suspected counterexample mechanism: a contact terminal that is
alive and part of a split at a small diameter but fully covered (→ CNO, not
open) at a larger diameter.

A violation is a (centre, d1, d2) where judgement(centre, d1) = open and
judgement(centre, d2) != open with d1 < d2.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from particle_sim.failure.auto_open import simulate_auto_open
from particle_sim.failure.open_failure import (
    _mask_bounds,
    open_candidates,
    particle_footprint,
)
from particle_sim.io.masks import GridSpec, MaskSet

# ── helpers ──────────────────────────────────────────────────────────


class _FakeHole:
    def __init__(
        self,
        contact_id: str,
        hole_index: int,
        lower_electrode: str,
        upper_electrode: str,
        ys: np.ndarray,
        xs: np.ndarray,
    ) -> None:
        self.contact_id = contact_id
        self.hole_index = hole_index
        self.lower_electrode = lower_electrode
        self.upper_electrode = upper_electrode
        self.ys = ys
        self.xs = xs


class _FakeTopology:
    def __init__(self, holes: list[_FakeHole]) -> None:
        self.holes = holes
        self.links: dict = {}
        self.lost_holes = 0


def _rect(
    y0: int, y1: int, x0: int, x1: int, shape: tuple[int, int]
) -> np.ndarray:
    mask = np.zeros(shape, dtype=bool)
    mask[y0:y1, x0:x1] = True
    return mask


def _masks_with_topology(
    electrode: np.ndarray,
    electrode_id: str,
    holes: list[_FakeHole],
    shape: tuple[int, int],
) -> MaskSet:
    return MaskSet(
        grid=GridSpec(
            shape_yx=shape, pixel_size_y_um=1.0, pixel_size_x_um=1.0
        ),
        layers={},
        electrodes={electrode_id: electrode},
        input_hashes={},
        contact_topology=_FakeTopology(holes),
    )


# ── shape generators ─────────────────────────────────────────────────


def _shape_narrow_wire_small_contacts(
    canvas: int,
) -> tuple[str, np.ndarray, list[_FakeHole], tuple[int, int]]:
    """Narrow wire with tiny contacts that a moderate particle can fully cover.

    This is the prime suspect: at small d a contact is alive and the wire is
    split → open.  At large d the contact is fully covered → CNO → not open.
    """
    shape = (canvas, canvas)
    mid_y = canvas // 2
    electrode = _rect(mid_y, mid_y + 2, 5, canvas - 5, shape)
    h1_ys, h1_xs = np.nonzero(_rect(mid_y, mid_y + 1, 8, 10, shape))
    h2_ys, h2_xs = np.nonzero(_rect(mid_y, mid_y + 1, canvas - 10, canvas - 8, shape))
    holes = [
        _FakeHole("C", 0, "X", "W", h1_ys, h1_xs),
        _FakeHole("C", 1, "X", "W", h2_ys, h2_xs),
    ]
    return "narrow_wire_small_contacts", electrode, holes, shape


def _shape_t_junction_tiny_contact(
    canvas: int,
) -> tuple[str, np.ndarray, list[_FakeHole], tuple[int, int]]:
    """T-junction: the branch has a tiny contact that can be fully covered."""
    shape = (canvas, canvas)
    mid_y = canvas // 2
    electrode = _rect(mid_y, mid_y + 2, 5, canvas - 5, shape)
    branch_x = canvas // 2
    electrode |= _rect(mid_y - 10, mid_y, branch_x, branch_x + 2, shape)
    h1_ys, h1_xs = np.nonzero(_rect(mid_y, mid_y + 1, 8, 10, shape))
    h2_ys, h2_xs = np.nonzero(_rect(mid_y, mid_y + 1, canvas - 10, canvas - 8, shape))
    h3_ys, h3_xs = np.nonzero(
        _rect(mid_y - 10, mid_y - 9, branch_x, branch_x + 1, shape)
    )
    holes = [
        _FakeHole("C", 0, "X", "W", h1_ys, h1_xs),
        _FakeHole("C", 1, "X", "W", h2_ys, h2_xs),
        _FakeHole("C", 2, "X", "W", h3_ys, h3_xs),
    ]
    return "t_junction_tiny_contact", electrode, holes, shape


def _shape_three_contacts_asymmetric(
    canvas: int,
) -> tuple[str, np.ndarray, list[_FakeHole], tuple[int, int]]:
    """Three contacts: one tiny (1 px), two normal (4 px).

    The tiny one is the easiest to cover fully.  A particle near it at small d
    might split the tiny contact from the others (open), but at large d the
    tiny contact is fully covered (CNO) and the other two remain connected.
    """
    shape = (canvas, canvas)
    mid_y = canvas // 2
    electrode = _rect(mid_y, mid_y + 4, 3, canvas - 3, shape)
    h1_ys, h1_xs = np.nonzero(_rect(mid_y + 1, mid_y + 3, 6, 8, shape))
    h2_ys, h2_xs = np.nonzero(_rect(mid_y + 1, mid_y + 3, canvas - 8, canvas - 6, shape))
    h3_ys, h3_xs = np.nonzero(_rect(mid_y + 1, mid_y + 2, canvas // 2, canvas // 2 + 1, shape))
    holes = [
        _FakeHole("C", 0, "X", "W", h1_ys, h1_xs),
        _FakeHole("C", 1, "X", "W", h2_ys, h2_xs),
        _FakeHole("C", 2, "X", "W", h3_ys, h3_xs),
    ]
    return "three_contacts_asymmetric", electrode, holes, shape


def _shape_two_contacts_one_boundary(
    canvas: int,
) -> tuple[str, np.ndarray, list[_FakeHole], tuple[int, int]]:
    """Wire touching ROI boundary (auto roi_boundary terminal) + two tiny contacts."""
    shape = (canvas, canvas)
    mid_y = canvas // 2
    electrode = _rect(mid_y, mid_y + 2, 0, canvas - 5, shape)
    h1_ys, h1_xs = np.nonzero(_rect(mid_y, mid_y + 1, 15, 17, shape))
    h2_ys, h2_xs = np.nonzero(_rect(mid_y, mid_y + 1, canvas - 10, canvas - 8, shape))
    holes = [
        _FakeHole("C", 0, "X", "W", h1_ys, h1_xs),
        _FakeHole("C", 1, "X", "W", h2_ys, h2_xs),
    ]
    return "two_contacts_one_boundary", electrode, holes, shape


def _shape_single_pixel_contact(
    canvas: int,
) -> tuple[str, np.ndarray, list[_FakeHole], tuple[int, int]]:
    """Extreme case: 1-pixel contact on a 2-pixel-wide wire.

    Even a diameter-3 particle can cover a single pixel.
    """
    shape = (canvas, canvas)
    mid_y = canvas // 2
    electrode = _rect(mid_y, mid_y + 2, 5, canvas - 5, shape)
    h1_ys = np.array([mid_y])
    h1_xs = np.array([10])
    h2_ys, h2_xs = np.nonzero(_rect(mid_y, mid_y + 2, canvas - 10, canvas - 6, shape))
    holes = [
        _FakeHole("C", 0, "X", "W", h1_ys, h1_xs),
        _FakeHole("C", 1, "X", "W", h2_ys, h2_xs),
    ]
    return "single_pixel_contact", electrode, holes, shape


ALL_SHAPES = [
    _shape_narrow_wire_small_contacts,
    _shape_t_junction_tiny_contact,
    _shape_three_contacts_asymmetric,
    _shape_two_contacts_one_boundary,
    _shape_single_pixel_contact,
]


# ── monotonicity probe ───────────────────────────────────────────────


def _probe_auto_open(
    name: str,
    masks: MaskSet,
    electrode_id: str,
    diameters: list[float],
    contact_margin_um: float,
    max_violations: int = 20,
) -> dict:
    """Check monotonicity of simulate_auto_open across ascending diameters."""
    results_by_d: dict[float, np.ndarray] = {}
    for d in sorted(diameters):
        result = simulate_auto_open(masks, electrode_id, d, contact_margin_um)
        results_by_d[d] = result.failure.copy()

    sorted_ds = sorted(diameters)
    checked = 0
    violations: list[dict] = []
    for i in range(len(sorted_ds)):
        d1 = sorted_ds[i]
        f1 = results_by_d[d1]
        for j in range(i + 1, len(sorted_ds)):
            d2 = sorted_ds[j]
            f2 = results_by_d[d2]
            open_at_d1 = f1.astype(bool)
            open_at_d2 = f2.astype(bool)
            violation_mask = open_at_d1 & ~open_at_d2
            checked += int(open_at_d1.sum())
            vy, vx = np.nonzero(violation_mask)
            for k in range(min(len(vy), max_violations - len(violations))):
                violations.append({
                    "centre_y": int(vy[k]),
                    "centre_x": int(vx[k]),
                    "d1": d1,
                    "d2": d2,
                    "label_d1": "open",
                    "label_d2": "not_open",
                })
            if len(violations) >= max_violations:
                break
        if len(violations) >= max_violations:
            break

    return {
        "shape": name,
        "function": "simulate_auto_open",
        "diameters": sorted_ds,
        "open_counts_by_diameter": {
            str(d): int(results_by_d[d].sum()) for d in sorted_ds
        },
        "checked_pairs": checked,
        "violation_count": len(violations),
        "violations": violations[:max_violations],
    }


def _probe_simulate_open_monotonicity(
    name: str,
    electrode: np.ndarray,
    terminal_a_yx: tuple[int, int],
    terminal_b_yx: tuple[int, int],
    grid: GridSpec,
    diameters: list[float],
    contact_margin_um: float,
    max_violations: int = 20,
) -> dict:
    """Check monotonicity of simulate_open (explicit 2-terminal) across diameters."""
    from particle_sim.failure.open_fast import LocalConnectivity

    y0, y1, x0, x1 = _mask_bounds(electrode)
    sorted_ds = sorted(diameters)
    results_by_d: dict[float, np.ndarray] = {}
    for d in sorted_ds:
        radius_um = d / 2.0
        stamp = particle_footprint(grid, radius_um, contact_margin_um)
        half_y, half_x = stamp.shape[0] // 2, stamp.shape[1] // 2
        height, width = electrode.shape
        crop_y0 = max(y0 - half_y, 0)
        crop_y1 = min(y1 + half_y, height)
        crop_x0 = max(x0 - half_x, 0)
        crop_x1 = min(x1 + half_x, width)
        cropped = electrode[crop_y0:crop_y1, crop_x0:crop_x1]
        ta_region = np.zeros(electrode.shape, dtype=bool)
        tb_region = np.zeros(electrode.shape, dtype=bool)
        ta_region[terminal_a_yx] = True
        tb_region[terminal_b_yx] = True
        cta = ta_region[crop_y0:crop_y1, crop_x0:crop_x1]
        ctb = tb_region[crop_y0:crop_y1, crop_x0:crop_x1]
        lc = LocalConnectivity(cropped, [cta, ctb], half_y, half_x)
        candidates = open_candidates(electrode, grid, radius_um, contact_margin_um)
        ys, xs = np.nonzero(candidates)
        failure = np.zeros(electrode.shape, dtype=bool)
        centres = [
            (int(y) - crop_y0, int(x) - crop_x0)
            for y, x in zip(ys.tolist(), xs.tolist(), strict=True)
        ]
        order = sorted(range(len(centres)), key=lambda i: lc.tile_key(*centres[i]))
        for i in order:
            cy, cx = centres[i]
            failure[int(ys[i]), int(xs[i])] = lc.splits((cy, cx), stamp)
        results_by_d[d] = failure

    checked = 0
    violations: list[dict] = []
    for i in range(len(sorted_ds)):
        d1 = sorted_ds[i]
        f1 = results_by_d[d1]
        for j in range(i + 1, len(sorted_ds)):
            d2 = sorted_ds[j]
            f2 = results_by_d[d2]
            violation_mask = f1 & ~f2
            checked += int(f1.sum())
            vy, vx = np.nonzero(violation_mask)
            for k in range(min(len(vy), max_violations - len(violations))):
                violations.append({
                    "centre_y": int(vy[k]),
                    "centre_x": int(vx[k]),
                    "d1": d1,
                    "d2": d2,
                    "label_d1": "open",
                    "label_d2": "not_open",
                })
            if len(violations) >= max_violations:
                break
        if len(violations) >= max_violations:
            break

    return {
        "shape": name,
        "function": "simulate_open",
        "diameters": sorted_ds,
        "open_counts_by_diameter": {
            str(d): int(results_by_d[d].sum()) for d in sorted_ds
        },
        "checked_pairs": checked,
        "violation_count": len(violations),
        "violations": violations[:max_violations],
    }


# ── main ─────────────────────────────────────────────────────────────


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Probe open judgement monotonicity across particle sizes."
    )
    parser.add_argument(
        "--canvas-size", type=int, default=50, help="Canvas side length in pixels"
    )
    parser.add_argument(
        "--sizes",
        type=str,
        default="2,3,4,5,6,8,10,12",
        help="Comma-separated diameters in µm",
    )
    parser.add_argument(
        "--output-json",
        type=str,
        default=None,
        help="Path to write JSON results",
    )
    parser.add_argument(
        "--work-root",
        type=str,
        default=None,
        help="(unused, kept for CLI convention parity)",
    )
    args = parser.parse_args()

    canvas = args.canvas_size
    diameters = [float(s) for s in args.sizes.split(",")]
    contact_margin = 0.0

    all_results: list[dict] = []
    total_checked = 0
    total_violations = 0

    for shape_fn in ALL_SHAPES:
        name, electrode, holes, shape = shape_fn(canvas)
        masks = _masks_with_topology(electrode, "W", holes, shape)
        print(f"[auto_open] {name} ...", end=" ", flush=True)
        result = _probe_auto_open(name, masks, "W", diameters, contact_margin)
        total_checked += result["checked_pairs"]
        total_violations += result["violation_count"]
        print(
            f"checked={result['checked_pairs']}, violations={result['violation_count']}"
        )
        all_results.append(result)

        mid_y = canvas // 2
        ta = (mid_y, 8)
        tb = (mid_y, canvas - 9)
        if electrode[ta] and electrode[tb]:
            grid = GridSpec(
                shape_yx=shape, pixel_size_y_um=1.0, pixel_size_x_um=1.0
            )
            print(f"[sim_open]  {name} ...", end=" ", flush=True)
            result2 = _probe_simulate_open_monotonicity(
                name, electrode, ta, tb, grid, diameters, contact_margin
            )
            total_checked += result2["checked_pairs"]
            total_violations += result2["violation_count"]
            print(
                f"checked={result2['checked_pairs']}, violations={result2['violation_count']}"
            )
            all_results.append(result2)

    summary = {
        "canvas_size": canvas,
        "diameters": diameters,
        "contact_margin_um": contact_margin,
        "total_shapes": len(ALL_SHAPES),
        "total_checked_pairs": total_checked,
        "total_violations": total_violations,
        "results": all_results,
    }

    print()
    print(f"Total checked pairs: {total_checked}")
    print(f"Total violations:    {total_violations}")

    if args.output_json:
        out = Path(args.output_json)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(summary, indent=2))
        print(f"Results written to {out}")

    sys.exit(0 if total_violations == 0 else 1)


if __name__ == "__main__":
    main()
