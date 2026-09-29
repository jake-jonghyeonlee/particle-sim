"""Bitwise compatibility with the pre-MEM-3B stack allocator."""

from __future__ import annotations

import struct

import numpy as np
import pytest

from particle_sim.aggregation.risk_allocation import allocate_risk
from particle_sim.failure.lateral_short import PairResult
from particle_sim.io.masks import GridSpec


def _old_allocate(results: list[PairResult], grid: GridSpec):
    if not results:
        raise ValueError("위험 배분에는 최소 한 개의 실패 단위가 필요합니다")
    shape = results[0].failure.shape
    if any(result.failure.shape != shape for result in results):
        raise ValueError("실패 단위 hotspot shape가 서로 다릅니다")
    stack = np.stack([result.failure.astype(bool) for result in results], axis=0)
    active_count = stack.sum(axis=0, dtype=np.uint32)
    union = active_count > 0
    denominator = np.where(union, active_count, 1).astype(np.float64)
    pixel_area = grid.pixel_area_um2
    raw: dict[str, float] = {}
    allocated: dict[str, float] = {}
    unique: dict[str, float] = {}
    for result, failure in zip(results, stack, strict=True):
        raw[result.pair_id] = float(failure.sum()) * pixel_area
        allocated[result.pair_id] = float(np.sum(failure / denominator)) * pixel_area
        unique[result.pair_id] = float(np.sum(failure & (active_count == 1))) * pixel_area
    return union, active_count, raw, allocated, unique


def _pair(index: int, failure: np.ndarray, pair_id: str | None = None) -> PairResult:
    return PairResult(
        pair_id or str(index), "lateral_short", "A", "B", failure, failure,
        failure, np.zeros(failure.shape, dtype=np.float32), "confirmed",
    )


def _assert_exact(results: list[PairResult], grid: GridSpec) -> None:
    old = _old_allocate(results, grid)
    new = allocate_risk(results, grid)
    assert new.union.dtype == old[0].dtype
    assert new.active_count.dtype == old[1].dtype
    assert new.union.tobytes() == old[0].tobytes()
    assert new.active_count.tobytes() == old[1].tobytes()
    for current, expected in zip(
        (new.raw_ca_um2, new.allocated_ca_um2, new.unique_ca_um2),
        old[2:], strict=True,
    ):
        assert list(current) == list(expected)
        assert [struct.pack("=d", value) for value in current.values()] == [
            struct.pack("=d", value) for value in expected.values()
        ]


@pytest.mark.parametrize("shape", [(1, 1), (7, 11), (193, 197), (200, 200)])
@pytest.mark.parametrize("count", [1, 3, 17, 96])
@pytest.mark.parametrize("layout", ["C", "F", "strided"])
def test_layout_and_overlap_oracle(shape, count, layout) -> None:
    # Cover the same shape, count, and layout matrix as the numerical experiment.
    rng = np.random.default_rng(30)
    maps = []
    for index in range(count):
        failure = rng.random(shape) < (0.0 if index == 0 else 1.0 if index == 1 else 0.2)
        if layout == "F":
            failure = np.asfortranarray(failure)
        elif layout == "strided":
            failure = np.repeat(failure, 2, axis=1)[:, ::2]
        maps.append(_pair(index, failure))
    _assert_exact(maps, GridSpec(shape, 0.3, 0.7))


def test_observed_six_fortran_counterexamples() -> None:
    # Preserve the experiment's single RNG stream across all 48 cases.
    rng = np.random.default_rng(30)
    mismatches = []
    for shape in ((1, 1), (7, 11), (193, 197), (200, 200)):
        for count in (1, 3, 17, 96):
            for layout in ("C", "F", "strided"):
                results = []
                for index in range(count):
                    failure = rng.random(shape) < (
                        0.0 if index == 0 else 1.0 if index == 1 else 0.2
                    )
                    if layout == "F":
                        failure = np.asfortranarray(failure)
                    elif layout == "strided":
                        failure = np.repeat(failure, 2, axis=1)[:, ::2]
                    results.append(_pair(index, failure))
                grid = GridSpec(shape, 0.3, 0.7)
                old = _old_allocate(results, grid)
                if layout == "F" and count in (17, 96) and shape != (1, 1):
                    streamed_count = np.zeros(shape, dtype=np.uint32)
                    for result in results:
                        streamed_count += result.failure
                    denominator = np.where(streamed_count > 0, streamed_count, 1).astype(
                        np.float64
                    )
                    unguarded = (
                        [float(result.failure.sum()) * grid.pixel_area_um2 for result in results],
                        [float(np.sum(result.failure / denominator)) * grid.pixel_area_um2
                         for result in results],
                        [float(np.sum(result.failure & (streamed_count == 1))) * grid.pixel_area_um2
                         for result in results],
                    )
                    if any(
                        struct.pack("=d", value) != struct.pack("=d", expected)
                        for values, expected_dict in zip(unguarded, old[2:], strict=True)
                        for value, expected in zip(values, expected_dict.values(), strict=True)
                    ):
                        mismatches.append((shape, count))
                _assert_exact(results, grid)
    assert mismatches == [
        ((7, 11), 17), ((7, 11), 96),
        ((193, 197), 17), ((193, 197), 96),
        ((200, 200), 17), ((200, 200), 96),
    ]


@pytest.mark.parametrize("kind", ["all_false", "all_true", "disjoint", "overlap"])
def test_pattern_oracle(kind) -> None:
    shape = (7, 11)
    maps = [np.zeros(shape, dtype=bool) for _ in range(7)]
    if kind == "all_true":
        for array in maps:
            array[:] = True
    elif kind == "disjoint":
        for index, array in enumerate(maps):
            array.flat[index] = True
    elif kind == "overlap":
        for index, array in enumerate(maps):
            array.flat[: index + 1] = True
    _assert_exact([_pair(index, array) for index, array in enumerate(maps)], GridSpec(shape, 0.3, 0.7))


def test_mixed_numeric_readonly_subclass_duplicate_ids() -> None:
    class BoolSubclass(np.ndarray):
        pass

    rng = np.random.default_rng(4)
    plain = rng.random((7, 11)) < 0.5
    readonly = plain.copy()
    readonly.flags.writeable = False
    numeric = plain.astype(np.int16) * 3
    subclass = plain.view(BoolSubclass)
    for maps in (
        [plain, readonly], [plain, np.asfortranarray(plain)],
        [plain, numeric], [plain, subclass],
    ):
        _assert_exact([_pair(i, value, "duplicate") for i, value in enumerate(maps)], GridSpec((7, 11), 0.3, 0.7))


def test_eligible_path_never_stacks(monkeypatch) -> None:
    def forbidden(*_args, **_kwargs):
        raise AssertionError("eligible path called np.stack")

    monkeypatch.setattr(np, "stack", forbidden)
    arrays = [np.ones((7, 11), dtype=bool) for _ in range(1000)]
    arrays[0].flags.writeable = False
    result = allocate_risk([_pair(i, array) for i, array in enumerate(arrays)], GridSpec((7, 11), 0.3, 0.7))
    assert result.active_count.tobytes() == np.full((7, 11), 1000, dtype=np.uint32).tobytes()


def test_ineligible_path_uses_original_stack(monkeypatch) -> None:
    original_stack = np.stack
    calls = 0

    def counted(*args, **kwargs):
        nonlocal calls
        calls += 1
        return original_stack(*args, **kwargs)

    monkeypatch.setattr(np, "stack", counted)
    failure = np.asfortranarray(np.ones((7, 11), dtype=bool))
    allocate_risk([_pair(0, failure)], GridSpec((7, 11), 0.3, 0.7))
    assert calls == 1


def test_errors_match_reference() -> None:
    grid = GridSpec((7, 11), 0.3, 0.7)
    for values in ([], [_pair(0, np.zeros((7, 11), bool)), _pair(1, np.zeros((7, 12), bool))]):
        with pytest.raises(ValueError) as old:
            _old_allocate(values, grid)
        with pytest.raises(ValueError) as new:
            allocate_risk(values, grid)
        assert str(new.value) == str(old.value)
