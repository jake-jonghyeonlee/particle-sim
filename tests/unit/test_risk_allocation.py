from __future__ import annotations

import numpy as np
import pytest

from particle_sim.aggregation.risk_allocation import allocate_risk
from particle_sim.failure.lateral_short import PairResult
from particle_sim.io.masks import GridSpec


def _result(unit_id: str, failure: np.ndarray) -> PairResult:
    return PairResult(
        pair_id=unit_id,
        failure_mode="lateral_short",
        electrode_a=f"{unit_id}A",
        electrode_b=f"{unit_id}B",
        failure=failure,
        contact_a=failure,
        contact_b=failure,
        pair_margin_um=np.zeros(failure.shape, dtype=np.float32),
        classification="confirmed",
    )


def test_allocation_sums_to_union_for_partial_overlap() -> None:
    a = np.zeros((4, 4), dtype=bool)
    b = np.zeros((4, 4), dtype=bool)
    a[0, 0:3] = True
    b[0, 2:4] = True
    grid = GridSpec((4, 4), 0.5, 0.5)
    result = allocate_risk([_result("A", a), _result("B", b)], grid)
    union_ca = float(result.union.sum()) * grid.pixel_area_um2
    assert sum(result.allocated_ca_um2.values()) == pytest.approx(union_ca)
    assert result.raw_ca_um2 == pytest.approx({"A": 0.75, "B": 0.5})
    assert result.allocated_ca_um2 == pytest.approx({"A": 0.625, "B": 0.375})
    assert result.unique_ca_um2 == pytest.approx({"A": 0.5, "B": 0.25})
