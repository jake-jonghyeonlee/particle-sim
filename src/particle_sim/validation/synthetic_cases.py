"""Synthetic fixture products (not physical simulation results)."""

from __future__ import annotations

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
