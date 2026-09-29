"""Full-crop contract cases for auto-open batching and future implementations."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.ndimage import label

from particle_sim.failure.auto_open import (
    _reference_auto_open_judgement,
    electrode_terminals,
    simulate_auto_open,
)
from particle_sim.failure.open_failure import _mask_bounds, open_candidates, particle_footprint
from particle_sim.io.masks import GridSpec, MaskSet


def _rect(y0: int, y1: int, x0: int, x1: int, shape=(32, 40)) -> np.ndarray:
    array = np.zeros(shape, dtype=bool)
    array[y0:y1, x0:x1] = True
    return array


def _case(name: str) -> MaskSet:
    wire = _rect(14, 18, 0, 40)
    branch = wire | _rect(4, 28, 18, 22)
    ring = _rect(4, 28, 4, 36)
    ring[9:23, 10:30] = False
    cases = {
        "branch": (branch, [_rect(5, 8, 18, 22), _rect(24, 27, 18, 22)]),
        "hole": (ring, [_rect(4, 8, 5, 9), _rect(24, 28, 31, 35)]),
        "interior": (_rect(10, 22, 3, 37), [_rect(14, 17, 6, 9), _rect(14, 17, 31, 34)]),
        "components": (wire | _rect(3, 6, 3, 9),
                       [_rect(14, 18, 4, 8), _rect(14, 18, 30, 34), _rect(3, 6, 3, 9)]),
        "corner": (_rect(0, 4, 0, 30) | _rect(0, 32, 0, 4), [_rect(0, 4, 20, 24)]),
        "edges": (wire, []),
        "bridge": (wire, [_rect(14, 18, 4, 9), _rect(14, 18, 7, 33),
                          _rect(14, 18, 31, 36)]),
        "solid": (np.ones((32, 40), bool), [_rect(6, 9, 8, 11), _rect(24, 27, 31, 34)]),
        "empty": (np.zeros((32, 40), bool), []),
        "disconnected_contact": (
            _rect(10, 12, 2, 38) | _rect(20, 22, 2, 38),
            [_rect(10, 12, 4, 8) | _rect(20, 22, 4, 8),
             _rect(10, 12, 30, 34), _rect(20, 22, 30, 34)],
        ),
        "tiny": (np.ones((2, 3), bool),
                 [_rect(0, 1, 0, 1, (2, 3)), _rect(1, 2, 2, 3, (2, 3))]),
    }
    if name.startswith("random"):
        seed = int(name.removeprefix("random"))
        rng = np.random.default_rng(seed)
        electrode = rng.random((32, 40)) < 0.68
        electrode[15, 1:39] = True
        regions = [_rect(14, 18, 4, 8), _rect(14, 18, 30, 34)]
    else:
        electrode, regions = cases[name]
    holes = []
    for index, region in enumerate(regions):
        ys, xs = np.nonzero(region & electrode)
        holes.append(SimpleNamespace(contact_id="C", hole_index=index,
                                     lower_electrode="X", upper_electrode="W", ys=ys, xs=xs))
    return MaskSet(
        grid=GridSpec(shape_yx=electrode.shape, pixel_size_y_um=1, pixel_size_x_um=1),
        layers={}, electrodes={"W": electrode}, input_hashes={},
        contact_topology=SimpleNamespace(holes=holes),
    )


def _groups(sets: list[set[int]]) -> list[list[int]]:
    """Independent transitive closure, separate from production _union_groups."""
    remaining = set(range(len(sets)))
    groups = []
    while remaining:
        component = {remaining.pop()}
        while True:
            touching = {i for i in remaining if any(sets[i] & sets[j] for j in component)}
            if not touching:
                break
            component |= touching
            remaining -= touching
        if len(component) >= 2:
            groups.append(sorted(component))
    return groups


def _compare_full_crop(masks: MaskSet, diameter: float, margin: float = 0.0) -> None:
    result = simulate_auto_open(masks, "W", diameter, margin)
    electrode = masks.electrodes["W"]
    terminals = electrode_terminals(masks, "W")
    assert result.electrode == "W"
    assert result.terminals == tuple(terminal.name for terminal in terminals)
    assert result.failure.dtype == np.bool_
    assert result.failure.shape == electrode.shape
    if len(terminals) < 2 or not electrode.any():
        assert result.classification == "insufficient_terminals"
        assert result.candidates == 0
        assert not result.failure.any()
        return
    y0, y1, x0, x1 = _mask_bounds(electrode)
    crop = electrode[y0:y1, x0:x1]
    local = [(t.kind, t.ys - y0, t.xs - x0) for t in terminals]
    base, _ = label(crop)
    groups = _groups([set(base[ys, xs].tolist()) - {0} for _, ys, xs in local])
    if not groups:
        assert result.classification == "insufficient_terminals"
        assert result.candidates == 0
        assert not result.failure.any()
        return
    candidates = open_candidates(electrode, masks.grid, diameter / 2, margin)
    ys, xs = np.nonzero(candidates)
    stamp = particle_footprint(masks.grid, diameter / 2, margin)
    reference = _reference_auto_open_judgement(
        crop, local, groups, stamp, ys, xs, y0, y1, x0, x1, *electrode.shape
    )
    different = np.argwhere(reference != result.failure)
    assert not different.size, different[:10].tolist()
    assert result.classification == "open_auto_terminals"
    assert result.candidates == int(candidates.sum())


@pytest.mark.parametrize("name", ["branch", "hole", "interior", "components", "corner",
                                  "edges", "bridge", "solid", "empty", "disconnected_contact",
                                  "tiny", "random0", "random1", "random2", "random3"])
@pytest.mark.parametrize("diameter", [1, 2, 5, 10, 20, 30, 50])
def test_full_crop_matrix(name: str, diameter: int) -> None:
    _compare_full_crop(_case(name), diameter)


@pytest.mark.parametrize("diameter,margin,spacing", [
    (0.0, 0.0, (1.0, 1.0)), (2.0, 0.03, (0.13, 0.07)),
    (0.2, -0.2, (0.1, 0.1)), (8.6, 0.0, (0.1, 0.1)),
    (8.6, 0.0, (0.13, 0.07)),
])
def test_numeric_footprints(diameter: float, margin: float, spacing: tuple[float, float]) -> None:
    masks = _case("branch")
    masks = replace(masks, grid=GridSpec(shape_yx=(32, 40), pixel_size_y_um=spacing[0],
                                        pixel_size_x_um=spacing[1]))
    _compare_full_crop(masks, diameter, margin)


def test_cno_reversal() -> None:
    electrode = np.zeros((14, 14), dtype=bool)
    electrode[7:9, 2:12] = True
    holes = [
        SimpleNamespace(contact_id="C", hole_index=0, lower_electrode="X", upper_electrode="W",
                        ys=np.array([7, 7, 8, 8]), xs=np.array([3, 4, 3, 4])),
        SimpleNamespace(contact_id="C", hole_index=1, lower_electrode="X", upper_electrode="W",
                        ys=np.array([7]), xs=np.array([10])),
    ]
    masks = MaskSet(grid=GridSpec(shape_yx=(14, 14), pixel_size_y_um=1,
                                   pixel_size_x_um=1), layers={}, electrodes={"W": electrode},
                    input_hashes={}, contact_topology=SimpleNamespace(holes=holes))
    _compare_full_crop(masks, 2)
    _compare_full_crop(masks, 3)
    assert simulate_auto_open(masks, "W", 2, 0).failure[7, 4]
    assert not simulate_auto_open(masks, "W", 3, 0).failure[7, 4]


@pytest.mark.parametrize("mutation,old,new", [
    ("cno_as_open", "# Rule 3: contact terminal fully covered → skip (CNO, not open)",
     "is_open = True; break  # mutant contact loss"),
    ("boundary_ignored", 'elif kind == "roi_boundary":', "elif False:"),
    ("candidate_count", "int(candidates.sum()))", "0)"),
])
def test_contract_rejects_disposable_mutant(
    tmp_path: Path, mutation: str, old: str, new: str
) -> None:
    root = Path(__file__).resolve().parents[2]
    disposable = tmp_path / mutation
    shutil.copytree(root / "src", disposable / "src")
    shutil.copy2(root / "pyproject.toml", disposable / "pyproject.toml")
    path = disposable / "src/particle_sim/failure/auto_open.py"
    code = path.read_text()
    assert old in code
    path.write_text(code.replace(old, new, 1))
    probe = disposable / "probe.py"
    probe.write_text(
        "from test_auto_open_batch import _case, _compare_full_crop, test_cno_reversal\n"
        + ("test_cno_reversal()\n" if mutation == "cno_as_open" else
           "_compare_full_crop(_case('edges'), 10)\n" if mutation == "boundary_ignored" else
           "_compare_full_crop(_case('branch'), 2)\n")
    )
    env = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1",
           "PYTHONPATH": os.pathsep.join((str(disposable / "src"), str(root / "tests/unit")))}
    completed = subprocess.run([sys.executable, str(probe)], cwd=disposable,
                               env=env, capture_output=True, text=True, timeout=120)
    assert completed.returncode != 0, f"mutant survived: {mutation}"
    assert "AssertionError" in completed.stderr, completed.stderr
