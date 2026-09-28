"""Exact candidate equivalence for EDT and the original disc dilation."""

from __future__ import annotations

import time

import numpy as np
import pytest

from particle_sim.failure.open_failure import _open_candidates_reference, open_candidates
from particle_sim.io.masks import GridSpec


def _shape(name: str) -> tuple[np.ndarray, tuple[tuple[int, int], ...]]:
    mask = np.zeros((71, 79), dtype=bool)
    terminals: tuple[tuple[int, int], ...] = ((0, 2), (70, 76))
    if name == "straight":
        mask[5:66, 39] = True
    elif name == "bent":
        mask[5:40, 20] = True
        mask[39, 20:60] = True
        mask[39:66, 59] = True
    elif name == "branched":
        mask[4:67, 39] = True
        mask[28, 10:69] = True
        mask[48, 20:60] = True
    elif name == "hole":
        mask[12:59, 14:65] = True
        mask[25:47, 27:52] = False
    elif name == "components":
        mask[7:22, 9:25] = True
        mask[49:65, 55:74] = True
    elif name == "border":
        mask[0, :] = True
        mask[:, 0] = True
        mask[-1, 35:79] = True
    return mask, terminals


@pytest.mark.parametrize(
    "shape", ["straight", "bent", "branched", "hole", "components", "border"]
)
@pytest.mark.parametrize("pixel_y,pixel_x", [(0.1, 0.1), (0.1, 0.2), (0.2, 0.1)])
@pytest.mark.parametrize("reach", [-0.1, 0.0, 0.1, 0.5, 1.3, 2.5, 5.0])
def test_open_candidates_equal_reference(
    shape: str, pixel_y: float, pixel_x: float, reach: float
) -> None:
    electrode, terminals = _shape(shape)
    grid = GridSpec(electrode.shape, pixel_y, pixel_x)
    actual = open_candidates(electrode, grid, reach - 0.07, 0.07, terminals)
    expected = _open_candidates_reference(electrode, grid, reach - 0.07, 0.07, terminals)
    np.testing.assert_array_equal(actual, expected)


def test_empty_seeds_and_terminal_only() -> None:
    electrode = np.zeros((21, 19), dtype=bool)
    grid = GridSpec(electrode.shape, 0.1, 0.2)
    np.testing.assert_array_equal(
        open_candidates(electrode, grid, 2.5, 0.0),
        _open_candidates_reference(electrode, grid, 2.5, 0.0),
    )
    np.testing.assert_array_equal(
        open_candidates(electrode, grid, 2.5, 0.0, ((0, 0),)),
        _open_candidates_reference(electrode, grid, 2.5, 0.0, ((0, 0),)),
    )


def test_lattice_radius_uses_reference_hypot_rounding() -> None:
    grid = GridSpec((61, 61), 0.1, 0.1)
    electrode = np.zeros(grid.shape_yx, dtype=bool)
    electrode[37, 54] = True  # (7, 24): hypot is one ULP above 2.5.
    assert not open_candidates(electrode, grid, 2.5, 0.0)[30, 30]
    electrode[45, 50] = True  # (15, 20): hypot is exactly 2.5.
    assert open_candidates(electrode, grid, 2.5, 0.0)[30, 30]
    np.testing.assert_array_equal(
        open_candidates(electrode, grid, 2.5, 0.0),
        _open_candidates_reference(electrode, grid, 2.5, 0.0),
    )


@pytest.mark.parametrize(
    "shape,pixel_y,pixel_x,centre,seeds,reach",
    [
        ((61, 61), 0.1, 0.1, (30, 30), ((39, 50), (45, 46)), np.hypot(0.9, 2.0)),
        ((61, 61), 0.1, 0.1, (30, 30), ((38, 55), (47, 50)), np.hypot(0.8, 2.5)),
        (
            (600, 600), 0.1, 0.1, (560, 560),
            ((568, 585), (577, 580)), np.hypot(0.8, 2.5),
        ),
        ((61, 61), 0.1, 0.2, (30, 30), ((39, 40), (45, 38)), np.hypot(0.9, 2.0)),
    ],
)
def test_equal_distance_seed_ties_match_reference(
    shape: tuple[int, int],
    pixel_y: float,
    pixel_x: float,
    centre: tuple[int, int],
    seeds: tuple[tuple[int, int], ...],
    reach: float,
) -> None:
    electrode = np.zeros(shape, dtype=bool)
    for seed in seeds:
        electrode[seed] = True
    grid = GridSpec(shape, pixel_y, pixel_x)
    actual = open_candidates(electrode, grid, reach, 0.0)
    expected = _open_candidates_reference(electrode, grid, reach, 0.0)
    assert expected[centre]
    assert actual[centre]
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize("dtype", [np.uint8, np.int64])
def test_integer_electrode_mask_matches_reference(dtype: np.dtype) -> None:
    electrode = np.zeros((21, 19), dtype=dtype)
    electrode[10, 9] = 1
    grid = GridSpec(electrode.shape, 0.1, 0.2)
    np.testing.assert_array_equal(
        open_candidates(electrode, grid, 0.5, 0.0),
        _open_candidates_reference(electrode, grid, 0.5, 0.0),
    )


@pytest.mark.parametrize("reach", [float("nan"), float("inf")])
def test_nonfinite_reach_raises_like_reference(reach: float) -> None:
    electrode = np.zeros((3, 3), dtype=bool)
    electrode[1, 1] = True
    grid = GridSpec(electrode.shape, 0.1, 0.1)
    with pytest.raises((ValueError, OverflowError)) as reference_error:
        _open_candidates_reference(electrode, grid, reach, 0.0)
    with pytest.raises(type(reference_error.value), match=str(reference_error.value)):
        open_candidates(electrode, grid, reach, 0.0)


def test_large_reach_finishes_quickly() -> None:
    electrode = np.zeros((1200, 1200), dtype=bool)
    electrode[:, 580:620] = True
    grid = GridSpec(electrode.shape, 0.1, 0.1)
    started = time.perf_counter()
    candidates = open_candidates(electrode, grid, 25.0, 0.0)
    elapsed = time.perf_counter() - started
    assert candidates.shape == electrode.shape
    assert candidates.dtype == bool
    assert candidates[:, 580:620].all()
    assert elapsed < 10.0
