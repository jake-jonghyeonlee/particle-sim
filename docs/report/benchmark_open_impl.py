"""Benchmark open_failure.simulate_open on a synthetic serpentine electrode.

Generates a single serpentine wire that crosses the ROI horizontally with turns,
two point terminals at opposite ends, and a schema_version '1.0' manifest.
Measures wall time for each particle diameter and writes a JSON report.
"""

from __future__ import annotations

import argparse
import json
import platform
import statistics
import sys
import time
from pathlib import Path

import numpy as np
import yaml
from PIL import Image

from particle_sim.failure.open_failure import simulate_open
from particle_sim.io.masks import load_product_masks
from particle_sim.schema import ProductManifest

PIXEL_UM = 0.5

STACK = [
    ("M1", "conductor", 0.3, None, True),
    ("GI", "inorganic", 0.4, 0.0, False),
    ("M2", "conductor", 0.3, None, True),
    ("PLN", "organic", 2.0, 0.8, False),
    ("M3", "conductor", 0.1, None, True),
]
STEPS = [
    ("S_M1", "M1", "conformal"),
    ("S_M2", "M2", "conformal"),
    ("S_PLN", "PLN", "flat"),
]


def _serpentine_electrode(
    canvas: int, width: int = 4, margin: int = 6
) -> tuple[np.ndarray, tuple[int, int], tuple[int, int]]:
    """Build a serpentine wire across the canvas.

    Returns (electrode_mask, terminal_a_yx, terminal_b_yx).
    """
    mask = np.zeros((canvas, canvas), dtype=bool)
    y = margin
    going_right = True
    terminal_a = (y, margin)
    last_point = terminal_a

    while y + width <= canvas - margin:
        if going_right:
            mask[y : y + width, margin : canvas - margin] = True
            last_point = (y + width // 2, canvas - margin - 1)
        else:
            mask[y : y + width, margin : canvas - margin] = True
            last_point = (y + width // 2, margin)

        next_y = y + width + margin
        if next_y + width > canvas - margin:
            break

        # vertical connector on the side
        if going_right:
            mask[y : next_y + width, canvas - margin - width : canvas - margin] = True
        else:
            mask[y : next_y + width, margin : margin + width] = True

        y = next_y
        going_right = not going_right

    terminal_b = last_point
    return mask, terminal_a, terminal_b


def _write_mask(path: Path, mask: np.ndarray) -> None:
    Image.fromarray(mask.astype(np.uint8) * 255, mode="L").save(path)


def prepare_case(
    work_root: Path, canvas: int
) -> tuple[ProductManifest, Path]:
    """Write a minimal synthetic case and return (manifest, manifest_path)."""
    case_dir = work_root / f"open-bench-{canvas}"
    masks_dir = case_dir / "masks"
    masks_dir.mkdir(parents=True, exist_ok=True)

    electrode, terminal_a, terminal_b = _serpentine_electrode(canvas)
    electrode_px = int(electrode.sum())

    # Write masks
    _write_mask(masks_dir / "electrode_W.png", electrode)

    layer_masks: dict[str, np.ndarray] = {}
    for layer, material, *_ in STACK:
        if material == "conductor":
            layer_mask = electrode.copy() if layer == "M2" else np.zeros_like(electrode)
            layer_masks[layer] = layer_mask
            _write_mask(masks_dir / f"layer_{layer}.png", layer_mask)

    manifest_dict: dict = {
        "schema_version": "1.0",
        "product_id": f"SYN_OPEN_BENCH_{canvas}",
        "product_family": "SYN",
        "design_version": "V01",
        "design_lineage_id": "SYN_L1",
        "fab_id": "FAB_SYN",
        "is_synthetic": True,
        "coordinate_system": {
            "origin": "upper_left",
            "x_direction": "right",
            "y_direction": "down",
            "pixel_size_x_um": PIXEL_UM,
            "pixel_size_y_um": PIXEL_UM,
            "roi_width_um": canvas * PIXEL_UM,
            "roi_height_um": canvas * PIXEL_UM,
        },
        "stack": [
            {
                "layer_name": layer,
                "process_order": order,
                "material_type": material,
                "thickness_um": thickness,
                **(
                    {}
                    if factor is None
                    else {"planarization_factor": factor}
                ),
                "patterned": patterned,
            }
            for order, (layer, material, thickness, factor, patterned) in enumerate(
                STACK
            )
        ],
        "layer_masks": [
            {"layer_name": layer, "path": f"masks/layer_{layer}.png"}
            for layer in layer_masks
        ],
        "electrodes": [
            {
                "electrode_id": "W",
                "layer_name": "M2",
                "net_id": "N1",
                "electrode_role": "DATA",
                "path": "masks/electrode_W.png",
            }
        ],
        "particle_source_steps": [
            {
                "step_id": sid,
                "lands_after_layer": layer,
                "surface_model": smodel,
                "fab_distribution_layer": sid,
            }
            for sid, layer, smodel in STEPS
        ],
        "electrode_pairs": [],
        "open_paths": [
            {
                "path_id": "W_PATH",
                "electrode": "W",
                "terminal_a_yx_px": list(terminal_a),
                "terminal_b_yx_px": list(terminal_b),
            }
        ],
        "particle_sizes_um": [2.0, 4.0, 6.0, 8.0, 10.0],
        "embedding_depth_um": 0.0,
        "contact_margin_um": 0.0,
        "post_step_policy": {
            "mode": "protrusion",
            "coverage_break_factor": 1.5,
        },
    }
    manifest_path = case_dir / "manifest.yaml"
    manifest_path.write_text(
        yaml.safe_dump(manifest_dict, sort_keys=False, allow_unicode=True),
        encoding="utf-8",
    )
    manifest = ProductManifest(**manifest_dict)
    print(f"  canvas={canvas}  electrode_px={electrode_px}")
    return manifest, manifest_path


def benchmark(
    manifest: ProductManifest,
    manifest_path: Path,
    sizes_um: list[float],
    repeats: int,
) -> dict:
    """Measure simulate_open wall time for each diameter.

    Returns a dict with per-size median, total, and open pixel counts.
    """
    masks = load_product_masks(manifest, manifest_path)
    open_path = manifest.open_paths[0]
    electrode = masks.electrodes[open_path.electrode]
    electrode_px = int(electrode.sum())

    per_size_times: dict[str, list[float]] = {str(s): [] for s in sizes_um}
    open_px: dict[str, int] = {}

    for repeat in range(repeats):
        for size in sizes_um:
            t0 = time.perf_counter()
            result = simulate_open(manifest, masks, open_path, size)
            elapsed = time.perf_counter() - t0
            per_size_times[str(size)].append(elapsed)
            if repeat == 0:
                open_px[str(size)] = int(result.failure.sum())

    per_size_median = {
        size: statistics.median(times) for size, times in per_size_times.items()
    }
    total = sum(per_size_median.values())

    return {
        "electrode_px": electrode_px,
        "per_size_median_s": per_size_median,
        "total_median_s": total,
        "open_px_per_size": open_px,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Benchmark open_failure.simulate_open on a serpentine electrode."
    )
    parser.add_argument("--canvas-size", type=int, default=200)
    parser.add_argument(
        "--sizes",
        type=str,
        default="2.0,4.0,6.0,8.0,10.0",
        help="Comma-separated particle diameters in um",
    )
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument(
        "--output-json",
        type=Path,
        default=None,
        help="Output JSON path (default: <work-root>/benchmark-open-impl-new.json)",
    )
    args = parser.parse_args()

    if args.repeats < 1:
        parser.error("--repeats must be at least 1")
    if args.canvas_size < 20:
        parser.error("--canvas-size must be at least 20")

    sizes_um = [float(s) for s in args.sizes.split(",")]
    work_root = args.work_root.resolve()
    work_root.mkdir(parents=True, exist_ok=True)
    output_json = (
        args.output_json
        if args.output_json is not None
        else work_root / "benchmark-open-impl-new.json"
    )

    print(f"Preparing case (canvas={args.canvas_size}) ...")
    manifest, manifest_path = prepare_case(work_root, args.canvas_size)

    print(f"Benchmarking ({args.repeats} repeats, sizes={sizes_um}) ...")
    result = benchmark(manifest, manifest_path, sizes_um, args.repeats)

    report = {
        "label": f"open-impl-bench-{args.canvas_size}x{args.canvas_size}",
        "canvas_size": args.canvas_size,
        "pixel_size_um": PIXEL_UM,
        "sizes_um": sizes_um,
        "repeats": args.repeats,
        "python": sys.version,
        "platform": platform.platform(),
        **result,
    }

    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(
        json.dumps(report, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(f"\nResult written to {output_json}")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
