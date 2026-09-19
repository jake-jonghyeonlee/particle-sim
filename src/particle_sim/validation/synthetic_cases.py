"""Synthetic fixture products (not physical simulation results)."""

from __future__ import annotations

import csv
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml
from PIL import Image

CASE_NAMES: tuple[str, ...] = (
    "parallel_lines",
    "vertical_aligned",
    "diagonal_offset",
    "triple_overlap",
    "pilot",
    "complex_multilayer",
    "open_straight",
    "open_detour",
    "open_branch",
)
SHAPE_YX = (200, 200)
PIXEL_UM = 0.5
PARTICLE_SIZES_UM = [2.0, 4.0, 6.0, 8.0, 10.0]

STACK = [
    ("M1", "conductor", 0.3, None, True),
    ("GI", "inorganic", 0.4, 0.0, False),
    ("M2", "conductor", 0.3, None, True),
    ("PLN", "organic", 2.0, 0.8, False),
    ("M3", "conductor", 0.1, None, True),
]
STEPS = [("S_M1", "M1", "conformal"), ("S_M2", "M2", "conformal"), ("S_PLN", "PLN", "flat")]

# case -> (electrodes, pairs); electrode = (id, layer, net, role, x0, x1)
_CASES: dict[str, tuple[list[tuple[str, str, str, str, int, int]], list[tuple[str, str, str]]]] = {
    "parallel_lines": (
        [("A", "M2", "N1", "DATA", 80, 90), ("B", "M2", "N2", "DATA", 96, 106)],
        [("A", "B", "lateral_short")],
    ),
    "vertical_aligned": (
        [("G", "M1", "N1", "GATE", 90, 110), ("P", "M3", "N2", "PIXEL", 90, 110)],
        [("G", "P", "vertical_short")],
    ),
    "diagonal_offset": (
        [("G", "M1", "N1", "GATE", 90, 110), ("P", "M3", "N2", "PIXEL", 116, 136)],
        [("G", "P", "vertical_short")],
    ),
    "triple_overlap": (
        [
            ("A", "M2", "N1", "DATA", 80, 90),
            ("B", "M2", "N2", "DATA", 92, 93),
            ("C", "M2", "N3", "DATA", 95, 105),
        ],
        [("A", "B", "lateral_short"), ("B", "C", "lateral_short"), ("A", "C", "lateral_short")],
    ),
    "pilot": (
        [
            ("G", "M1", "N1", "GATE", 40, 60),
            ("A", "M2", "N2", "DATA", 66, 76),
            ("B", "M2", "N3", "DATA", 82, 92),
            ("P", "M3", "N4", "PIXEL", 98, 118),
        ],
        [("A", "B", "lateral_short"), ("G", "A", "vertical_short"), ("B", "P", "vertical_short")],
    ),
}


# Open cases draw a wire on M2 as pixel rectangles instead of full-height stripes.
# (name -> list of (y0, y1, x0, x1) blocks, terminals in pixels)
_OPEN_SHAPES: dict[str, tuple[list[tuple[int, int, int, int]], tuple[int, int], tuple[int, int]]] = {
    # straight 6 px wide line; a large particle cuts it
    "open_straight": ([(97, 103, 20, 180)], (100, 22), (100, 178)),
    # same line plus a wide detour that keeps the terminals connected
    "open_detour": (
        [(97, 103, 20, 180), (60, 66, 60, 140), (60, 103, 60, 66), (60, 103, 134, 140)],
        (100, 22),
        (100, 178),
    ),
    # main line with a dead-end branch: cutting the branch must not count as open
    "open_branch": (
        [(97, 103, 20, 180), (40, 100, 100, 106)],
        (100, 22),
        (100, 178),
    ),
}


@dataclass(frozen=True)
class SyntheticCase:
    name: str
    manifest: dict
    layer_masks: dict[str, np.ndarray]
    electrode_masks: dict[str, np.ndarray]


def _blocks(blocks: list[tuple[int, int, int, int]]) -> np.ndarray:
    mask = np.zeros(SHAPE_YX, dtype=bool)
    for y0, y1, x0, x1 in blocks:
        mask[y0:y1, x0:x1] = True
    return mask


def _stripe(x0: int, x1: int) -> np.ndarray:
    mask = np.zeros(SHAPE_YX, dtype=bool)
    mask[:, x0:x1] = True
    return mask


def _complex_multilayer_masks() -> tuple[dict[str, np.ndarray], dict[str, str]]:
    """Ten non-overlapping electrodes with branches, bends, islands and layer crossings."""
    blocks_by_electrode = {
        "G1": [(18, 28, 15, 88), (18, 115, 15, 25), (70, 80, 25, 70)],
        "G2": [(42, 52, 70, 150), (42, 135, 105, 115), (125, 135, 75, 115)],
        "G3": [(145, 165, 28, 48), (145, 165, 72, 92)],
        "D1": [(35, 45, 30, 105), (35, 100, 95, 105), (90, 100, 55, 105)],
        "D2": [(105, 115, 30, 100), (105, 175, 30, 40), (150, 160, 40, 85)],
        "D3": [(45, 55, 125, 185), (45, 135, 150, 160), (85, 95, 120, 150)],
        "D4": [(145, 165, 115, 137), (145, 165, 168, 190), (155, 160, 137, 168)],
        "P1": [(25, 65, 42, 82), (55, 75, 60, 98)],
        "P2": [(80, 95, 105, 130), (95, 110, 120, 145), (110, 125, 135, 160)],
        "P3": [(135, 145, 100, 185), (135, 185, 175, 185), (175, 185, 125, 185)],
    }
    layer_by_electrode = {
        "G1": "M1", "G2": "M1", "G3": "M1",
        "D1": "M2", "D2": "M2", "D3": "M2", "D4": "M2",
        "P1": "M3", "P2": "M3", "P3": "M3",
    }
    return {key: _blocks(value) for key, value in blocks_by_electrode.items()}, layer_by_electrode


def _build_complex_multilayer_case() -> SyntheticCase:
    electrode_masks, layer_by_electrode = _complex_multilayer_masks()
    layer_masks = {
        layer: np.zeros(SHAPE_YX, dtype=bool)
        for layer, material, *_ in STACK
        if material == "conductor"
    }
    for electrode_id, mask in electrode_masks.items():
        layer_masks[layer_by_electrode[electrode_id]] |= mask

    roles = {"M1": "GATE", "M2": "DATA", "M3": "PIXEL"}
    pair_specs = [
        ("G1", "G2", "lateral_short"),
        ("G2", "G3", "lateral_short"),
        ("D1", "D2", "lateral_short"),
        ("D1", "D3", "lateral_short"),
        ("D2", "D4", "lateral_short"),
        ("D3", "D4", "lateral_short"),
        ("G1", "D1", "vertical_short"),
        ("G2", "D1", "vertical_short"),
        ("G2", "D3", "vertical_short"),
        ("G3", "D2", "vertical_short"),
        ("D1", "P1", "vertical_short"),
        ("D2", "P3", "vertical_short"),
        ("D3", "P2", "vertical_short"),
        ("D4", "P3", "vertical_short"),
        ("G1", "D4", "vertical_short"),
        ("G3", "D3", "vertical_short"),
    ]
    manifest = _base_manifest(
        "complex_multilayer",
        electrodes=[
            {
                "electrode_id": electrode_id,
                "layer_name": layer,
                "net_id": f"NET_{electrode_id}",
                "electrode_role": roles[layer],
                "path": f"masks/electrode_{electrode_id}.png",
            }
            for electrode_id, layer in layer_by_electrode.items()
        ],
        layer_masks=layer_masks,
        pairs=[
            {
                "pair_id": f"{a}__{b}",
                "electrode_a": a,
                "electrode_b": b,
                "failure_mode": mode,
            }
            for a, b, mode in pair_specs
        ],
        open_paths=[],
    )
    return SyntheticCase(
        name="complex_multilayer",
        manifest=manifest,
        layer_masks=layer_masks,
        electrode_masks=electrode_masks,
    )


def _base_manifest(
    name: str,
    *,
    electrodes: list[dict],
    layer_masks: dict[str, np.ndarray],
    pairs: list[dict],
    open_paths: list[dict],
) -> dict:
    """Manifest dict shared by every synthetic case."""
    return {
        "schema_version": "1.0",
        "product_id": f"SYN_{name.upper()}",
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
            "roi_width_um": SHAPE_YX[1] * PIXEL_UM,
            "roi_height_um": SHAPE_YX[0] * PIXEL_UM,
        },
        "stack": [
            {
                "layer_name": layer,
                "process_order": order,
                "material_type": material,
                "thickness_um": thickness,
                **({} if factor is None else {"planarization_factor": factor}),
                "patterned": patterned,
            }
            for order, (layer, material, thickness, factor, patterned) in enumerate(STACK)
        ],
        "layer_masks": [
            {"layer_name": layer, "path": f"masks/layer_{layer}.png"} for layer in layer_masks
        ],
        "electrodes": electrodes,
        "particle_source_steps": [
            {
                "step_id": step_id,
                "lands_after_layer": layer,
                "surface_model": surface_model,
                "fab_distribution_layer": step_id,
            }
            for step_id, layer, surface_model in STEPS
        ],
        "electrode_pairs": pairs,
        "open_paths": open_paths,
        "particle_sizes_um": list(PARTICLE_SIZES_UM),
        "embedding_depth_um": 0.0,
        "contact_margin_um": 0.0,
        "post_step_policy": {"mode": "protrusion", "coverage_break_factor": 1.5},
    }


def _build_open_case(name: str) -> SyntheticCase:
    """Single-electrode wire case used for open connectivity tests."""
    blocks, terminal_a, terminal_b = _OPEN_SHAPES[name]
    wire = _blocks(blocks)
    electrode_masks = {"W": wire}
    layer_masks = {
        layer: (wire.copy() if layer == "M2" else np.zeros(SHAPE_YX, dtype=bool))
        for layer, material, *_ in STACK
        if material == "conductor"
    }
    manifest = _base_manifest(
        name,
        electrodes=[
            {
                "electrode_id": "W",
                "layer_name": "M2",
                "net_id": "N1",
                "electrode_role": "DATA",
                "path": "masks/electrode_W.png",
            }
        ],
        layer_masks=layer_masks,
        pairs=[],
        open_paths=[
            {
                "path_id": "W_PATH",
                "electrode": "W",
                "terminal_a_yx_px": list(terminal_a),
                "terminal_b_yx_px": list(terminal_b),
            }
        ],
    )
    return SyntheticCase(
        name=name,
        manifest=manifest,
        layer_masks=layer_masks,
        electrode_masks=electrode_masks,
    )


def build_case(name: str) -> SyntheticCase:
    """Build one synthetic product manifest with its layer and electrode masks."""
    if name in _OPEN_SHAPES:
        return _build_open_case(name)
    if name == "complex_multilayer":
        return _build_complex_multilayer_case()
    if name not in _CASES:
        raise ValueError(f"알 수 없는 합성 case: {name!r} (허용: {', '.join(CASE_NAMES)})")
    electrode_specs, pair_specs = _CASES[name]
    electrode_masks = {spec[0]: _stripe(spec[4], spec[5]) for spec in electrode_specs}
    layer_masks = {
        layer: np.zeros(SHAPE_YX, dtype=bool)
        for layer, material, *_ in STACK
        if material == "conductor"
    }
    for spec in electrode_specs:
        layer_masks[spec[1]] |= electrode_masks[spec[0]]

    manifest = _base_manifest(
        name,
        electrodes=[
            {
                "electrode_id": spec[0],
                "layer_name": spec[1],
                "net_id": spec[2],
                "electrode_role": spec[3],
                "path": f"masks/electrode_{spec[0]}.png",
            }
            for spec in electrode_specs
        ],
        layer_masks=layer_masks,
        pairs=[
            {
                "pair_id": f"{a}__{b}",
                "electrode_a": a,
                "electrode_b": b,
                "failure_mode": mode,
            }
            for a, b, mode in pair_specs
        ],
        open_paths=[],
    )
    return SyntheticCase(
        name=name,
        manifest=manifest,
        layer_masks=layer_masks,
        electrode_masks=electrode_masks,
    )


def tile_case(case: SyntheticCase, repeats_yx: tuple[int, int]) -> SyntheticCase:
    """Tile a fixture without changing pixel size, for deterministic scale measurements."""
    repeat_y, repeat_x = repeats_yx
    if repeat_y < 1 or repeat_x < 1:
        raise ValueError(f"반복 수는 1 이상이어야 합니다: {repeats_yx}")
    manifest = deepcopy(case.manifest)
    shape_y = next(iter(case.layer_masks.values())).shape[0] * repeat_y
    shape_x = next(iter(case.layer_masks.values())).shape[1] * repeat_x
    manifest["product_id"] = f"{manifest['product_id']}_{shape_y}X{shape_x}"
    coordinate = manifest["coordinate_system"]
    coordinate["roi_width_um"] = shape_x * coordinate["pixel_size_x_um"]
    coordinate["roi_height_um"] = shape_y * coordinate["pixel_size_y_um"]
    return SyntheticCase(
        name=f"{case.name}_{shape_y}x{shape_x}",
        manifest=manifest,
        layer_masks={key: np.tile(mask, repeats_yx) for key, mask in case.layer_masks.items()},
        electrode_masks={
            key: np.tile(mask, repeats_yx) for key, mask in case.electrode_masks.items()
        },
    )


FAB_COLUMNS: tuple[str, ...] = (
    "fab_id", "product_family", "source_step", "size_bin_left_um", "size_bin_right_um",
    "size_representative_um", "particle_count", "inspection_area_mm2", "measurement_period",
    "density_per_mm2", "particle_type",
)
# (fab_id, counts per representative size) — synthetic, not measured data.
_FAB_COUNTS: dict[str, list[int]] = {
    "FAB_SYN": [120, 60, 30, 12, 4],
    "FAB_REF": [80, 40, 20, 10, 5],
}
FAB_INSPECTION_AREA_MM2 = 100.0


def write_fab_table(out_dir: Path, *, file_name: str = "fab_particles.csv") -> Path:
    """Write a synthetic Fab particle table covering every synthetic source step."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / file_name
    rows: list[dict[str, object]] = []
    for fab_id, counts in _FAB_COUNTS.items():
        for step_id, _layer, _model in STEPS:
            for size, count in zip(PARTICLE_SIZES_UM, counts, strict=True):
                rows.append(
                    {
                        "fab_id": fab_id,
                        "product_family": "SYN",
                        "source_step": step_id,
                        "size_bin_left_um": size - 1.0,
                        "size_bin_right_um": size + 1.0,
                        "size_representative_um": size,
                        "particle_count": count,
                        "inspection_area_mm2": FAB_INSPECTION_AREA_MM2,
                        "measurement_period": "2026-09 synthetic",
                        "density_per_mm2": count / FAB_INSPECTION_AREA_MM2,
                        "particle_type": "unknown",
                    }
                )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(FAB_COLUMNS))
        writer.writeheader()
        writer.writerows(rows)
    return path


def _write_mask(path: Path, mask: np.ndarray) -> None:
    Image.fromarray((mask.astype(np.uint8) * 255), mode="L").save(path)


def write_case(case: SyntheticCase, out_dir: Path) -> Path:
    """Write masks and manifest.yaml under out_dir and return the manifest path."""
    out_dir = Path(out_dir)
    (out_dir / "masks").mkdir(parents=True, exist_ok=True)
    for layer, mask in case.layer_masks.items():
        _write_mask(out_dir / "masks" / f"layer_{layer}.png", mask)
    for electrode_id, mask in case.electrode_masks.items():
        _write_mask(out_dir / "masks" / f"electrode_{electrode_id}.png", mask)
    manifest_path = out_dir / "manifest.yaml"
    manifest_path.write_text(
        yaml.safe_dump(case.manifest, sort_keys=False, allow_unicode=True), encoding="utf-8"
    )
    return manifest_path
