"""legacy-distance (T4) 단위 테스트."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.ndimage import distance_transform_cdt

from particle_sim.errors import InputContractError
from particle_sim.geometry.distance import distance_to_mask_um, distance_to_z_interval_um
from particle_sim.geometry.legacy import (
    edt_short_hotspot,
    legacy_dilate,
    legacy_effective_radius_um,
    legacy_erode,
    legacy_iterations,
    legacy_short_hotspot,
    legacy_structure,
)
from particle_sim.io.masks import GridSpec
from particle_sim.schema import LegacyKernel
from particle_sim.validation.regression import (
    CAUSE_CODES,
    compare_binary_maps,
    explain_legacy_vs_edt,
)


def _grid(py: float = 1.0, px: float = 1.0, shape: tuple[int, int] = (10, 10)) -> GridSpec:
    return GridSpec(shape_yx=shape, pixel_size_y_um=py, pixel_size_x_um=px, origin="upper_left")


# ── distance_to_mask_um ──────────────────────────────────────────────────────


class TestDistanceToMaskUm:
    def test_inside_mask_is_zero(self) -> None:
        mask = np.zeros((5, 5), dtype=bool)
        mask[2, 2] = True
        dist = distance_to_mask_um(mask, _grid(shape=(5, 5)))
        assert dist[2, 2] == 0.0

    def test_3_4_5_distance(self) -> None:
        """(0,0)에 True, (3,4)까지 거리가 5.0 (pixel 1 um)."""
        mask = np.zeros((10, 10), dtype=bool)
        mask[0, 0] = True
        dist = distance_to_mask_um(mask, _grid(py=1.0, px=1.0))
        assert dist[3, 4] == pytest.approx(5.0)

    def test_anisotropic_sampling(self) -> None:
        """py=2, px=1: y 방향 2배, x 방향 1배."""
        mask = np.zeros((10, 10), dtype=bool)
        mask[0, 0] = True
        dist = distance_to_mask_um(mask, _grid(py=2.0, px=1.0))
        assert dist[5, 0] == pytest.approx(10.0)
        assert dist[0, 5] == pytest.approx(5.0)

    def test_empty_mask_returns_inf(self) -> None:
        mask = np.zeros((5, 5), dtype=bool)
        dist = distance_to_mask_um(mask, _grid(shape=(5, 5)))
        assert np.all(np.isinf(dist))

    def test_full_mask_returns_zeros(self) -> None:
        mask = np.ones((5, 5), dtype=bool)
        dist = distance_to_mask_um(mask, _grid(shape=(5, 5)))
        assert np.all(dist == 0.0)

    def test_shape_mismatch_raises(self) -> None:
        mask = np.zeros((5, 5), dtype=bool)
        with pytest.raises(ValueError):
            distance_to_mask_um(mask, _grid(shape=(3, 3)))


# ── distance_to_z_interval_um ────────────────────────────────────────────────


class TestDistanceToZIntervalUm:
    @pytest.mark.parametrize(
        ("z", "z_bottom", "z_top", "expected"),
        [
            (5.0, 0.0, 10.0, 0.0),
            (-3.0, 0.0, 10.0, 3.0),
            (15.0, 0.0, 10.0, 5.0),
        ],
    )
    def test_scalar(self, z: float, z_bottom: float, z_top: float, expected: float) -> None:
        assert distance_to_z_interval_um(z, z_bottom, z_top) == pytest.approx(expected)

    def test_array_input(self) -> None:
        z = np.array([-1.0, 5.0, 12.0])
        result = distance_to_z_interval_um(z, 0.0, 10.0)
        assert result[0] == pytest.approx(1.0)
        assert result[1] == 0.0
        assert result[2] == pytest.approx(2.0)

    def test_invalid_interval_raises(self) -> None:
        with pytest.raises(ValueError):
            distance_to_z_interval_um(5.0, 10.0, 0.0)


# ── legacy_structure ─────────────────────────────────────────────────────────


class TestLegacyStructure:
    def test_square_5x5(self) -> None:
        s = legacy_structure(LegacyKernel(size=5, shape="square"))
        assert s.shape == (5, 5)
        assert s.all()

    def test_disk_5x5_corner_less(self) -> None:
        """5x5 disk = square without 4 corners (21 True)."""
        s = legacy_structure(LegacyKernel(size=5, shape="disk"))
        assert s.shape == (5, 5)
        assert not s[0, 0] and not s[0, 4] and not s[4, 0] and not s[4, 4]
        assert s[2, 2] and s[0, 2] and s[2, 0]
        assert s.sum() == 21


# ── legacy_iterations ────────────────────────────────────────────────────────


class TestLegacyIterations:
    def test_d2_pixel05_size5_all_rules_give_1(self) -> None:
        """d=2, pixel 0.5, size 5: radius_px=2, x=1 -> n=1 for all rules."""
        grid = _grid(py=0.5, px=0.5)
        for mode in ("round", "ceil", "floor"):
            kernel = LegacyKernel(size=5, shape="square", iterations=mode)
            assert legacy_iterations(2.0, grid, kernel) == 1

    @pytest.mark.parametrize(
        ("diameter", "mode", "expected"),
        [
            (3.0, "round", 2),
            (3.0, "ceil", 2),
            (3.0, "floor", 1),
            (2.5, "round", 1),
            (2.5, "ceil", 2),
            (2.5, "floor", 1),
        ],
    )
    def test_iteration_rules(self, diameter: float, mode: str, expected: int) -> None:
        grid = _grid(py=0.5, px=0.5)
        kernel = LegacyKernel(size=5, shape="square", iterations=mode)
        assert legacy_iterations(diameter, grid, kernel) == expected

    def test_anisotropic_raises(self) -> None:
        grid = _grid(py=2.0, px=1.0)
        kernel = LegacyKernel(size=5, shape="square")
        with pytest.raises(InputContractError) as exc_info:
            legacy_iterations(2.0, grid, kernel)
        assert exc_info.value.field == "coordinate_system.pixel_size"


# ── legacy_effective_radius_um ───────────────────────────────────────────────


class TestLegacyEffectiveRadiusUm:
    def test_d2_pixel05_size5_round(self) -> None:
        """n=1, h=2, pixel=0.5 -> r_eff = 1.0."""
        grid = _grid(py=0.5, px=0.5)
        kernel = LegacyKernel(size=5, shape="square", iterations="round")
        assert legacy_effective_radius_um(2.0, grid, kernel) == pytest.approx(1.0)

    def test_d3_pixel05_size5_round(self) -> None:
        """n=2, h=2, pixel=0.5 -> r_eff = 2.0."""
        grid = _grid(py=0.5, px=0.5)
        kernel = LegacyKernel(size=5, shape="square", iterations="round")
        assert legacy_effective_radius_um(3.0, grid, kernel) == pytest.approx(2.0)


# ── legacy_dilate / legacy_erode ─────────────────────────────────────────────


class TestLegacyDilateErode:
    def test_dilate_zero_iterations_returns_copy(self) -> None:
        mask = np.zeros((5, 5), dtype=bool)
        mask[2, 2] = True
        grid = _grid()
        kernel = LegacyKernel(size=5, shape="square", iterations="floor")
        result = legacy_dilate(mask, 0.1, grid, kernel)
        assert result is not mask
        assert np.array_equal(result, mask)

    def test_erode_zero_iterations_returns_copy(self) -> None:
        mask = np.zeros((5, 5), dtype=bool)
        mask[2, 2] = True
        grid = _grid()
        kernel = LegacyKernel(size=5, shape="square", iterations="floor")
        result = legacy_erode(mask, 0.1, grid, kernel)
        assert result is not mask
        assert np.array_equal(result, mask)

    def test_dilate_expands(self) -> None:
        mask = np.zeros((10, 10), dtype=bool)
        mask[5, 5] = True
        grid = _grid()
        kernel = LegacyKernel(size=5, shape="square", iterations="round")
        result = legacy_dilate(mask, 2.0, grid, kernel)
        assert result[5, 5]
        assert result.sum() > 1

    def test_square_kernel_matches_chessboard_cdt(self) -> None:
        """square kernel n iterations == chessboard CDT <= 2n."""
        mask = np.zeros((10, 10), dtype=bool)
        mask[5, 5] = True
        grid = _grid()
        kernel = LegacyKernel(size=5, shape="square", iterations="round")
        n = legacy_iterations(2.0, grid, kernel)
        assert n == 1
        dilated = legacy_dilate(mask, 2.0, grid, kernel)
        cdt = distance_transform_cdt(~mask, metric="chessboard")
        assert np.all(cdt[dilated] <= 2 * n)
        assert np.array_equal(dilated, cdt <= 2 * n)


# ── legacy_short_hotspot / edt_short_hotspot ─────────────────────────────────


class TestShortHotspot:
    def test_legacy_hotspot_no_overlap(self) -> None:
        mask_a = np.zeros((10, 10), dtype=bool)
        mask_b = np.zeros((10, 10), dtype=bool)
        mask_a[5, 3] = True
        mask_b[5, 9] = True
        grid = _grid()
        kernel = LegacyKernel(size=5, shape="square", iterations="round")
        result = legacy_short_hotspot(mask_a, mask_b, 4.0, grid, kernel)
        assert not result.any()

    def test_edt_hotspot(self) -> None:
        mask_a = np.zeros((10, 10), dtype=bool)
        mask_b = np.zeros((10, 10), dtype=bool)
        mask_a[5, 3] = True
        mask_b[5, 7] = True
        grid = _grid()
        result = edt_short_hotspot(mask_a, mask_b, 2.0, grid)
        assert result[5, 5]


# ── compare_binary_maps ──────────────────────────────────────────────────────


class TestCompareBinaryMaps:
    def test_identical_maps(self) -> None:
        m = np.array([[True, False], [False, True]])
        result = compare_binary_maps(m, m)
        assert result.equal
        assert result.n_diff == 0
        assert result.only_in_reference == 0
        assert result.only_in_candidate == 0
        assert result.cause_counts["unexplained"] == 0

    @pytest.mark.parametrize(
        ("ref", "cand", "only_ref", "only_cand"),
        [
            (np.array([[True, True], [False, False]]), np.array([[True, False], [False, False]]), 1, 0),
            (np.array([[True, False], [False, False]]), np.array([[True, True], [False, False]]), 0, 1),
        ],
    )
    def test_directional_diffs(
        self, ref: np.ndarray, cand: np.ndarray, only_ref: int, only_cand: int
    ) -> None:
        result = compare_binary_maps(ref, cand)
        assert not result.equal
        assert result.n_diff == 1
        assert result.only_in_reference == only_ref
        assert result.only_in_candidate == only_cand

    def test_shape_mismatch_raises(self) -> None:
        with pytest.raises(ValueError):
            compare_binary_maps(np.zeros((3, 3), dtype=bool), np.zeros((2, 2), dtype=bool))


# ── explain_legacy_vs_edt ────────────────────────────────────────────────────


class TestExplainLegacyVsEdt:
    def test_parallel_stripes_equal_when_r_eff_equals_r(self) -> None:
        """pixel 0.5, d=2: r_eff = 1.0 = r -> equal, no unexplained."""
        mask_a = np.zeros((20, 20), dtype=bool)
        mask_b = np.zeros((20, 20), dtype=bool)
        # parallel stripes: a at x=5, b at x=9 (4 px apart = 2 um)
        mask_a[:, 5] = True
        mask_b[:, 9] = True
        grid = _grid(py=0.5, px=0.5, shape=(20, 20))
        kernel = LegacyKernel(size=5, shape="square", iterations="round")
        result = explain_legacy_vs_edt(mask_a, mask_b, 2.0, grid, kernel)
        assert result.equal
        assert result.cause_counts["unexplained"] == 0

    def test_cause_counts_sum_equals_n_diff(self) -> None:
        mask_a = np.zeros((10, 10), dtype=bool)
        mask_b = np.zeros((10, 10), dtype=bool)
        mask_a[5, 3] = True
        mask_b[5, 6] = True
        grid = _grid(py=0.5, px=0.5)
        kernel = LegacyKernel(size=5, shape="square", iterations="round")
        result = explain_legacy_vs_edt(mask_a, mask_b, 3.0, grid, kernel)
        assert sum(result.cause_counts.values()) == result.n_diff

    def test_cause_map_dtype_and_shape(self) -> None:
        mask_a = np.zeros((10, 10), dtype=bool)
        mask_b = np.zeros((10, 10), dtype=bool)
        mask_a[5, 3] = True
        mask_b[5, 6] = True
        grid = _grid(py=0.5, px=0.5)
        kernel = LegacyKernel(size=5, shape="square", iterations="round")
        result = explain_legacy_vs_edt(mask_a, mask_b, 3.0, grid, kernel)
        assert result.cause_map.shape == (10, 10)
        assert result.cause_map.dtype == np.uint8

    def test_cause_codes_valid(self) -> None:
        mask_a = np.zeros((10, 10), dtype=bool)
        mask_b = np.zeros((10, 10), dtype=bool)
        mask_a[5, 3] = True
        mask_b[5, 6] = True
        grid = _grid(py=0.5, px=0.5)
        kernel = LegacyKernel(size=5, shape="square", iterations="round")
        result = explain_legacy_vs_edt(mask_a, mask_b, 3.0, grid, kernel)
        valid_codes = set(CAUSE_CODES.values())
        assert set(np.unique(result.cause_map).tolist()) <= valid_codes
