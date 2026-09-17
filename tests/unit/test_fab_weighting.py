from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from particle_sim.aggregation.fab_weighting import (
    WEIGHTED_COLUMNS,
    apply_fab_weighting,
    expected_defect_count,
    match_sizes,
)
from particle_sim.errors import InputContractError
from particle_sim.io.fab import FAB_COLUMNS, fab_density_for, load_fab_table


def _fab_rows(**overrides) -> list[dict]:
    rows = []
    for index, size in enumerate((2.0, 4.0)):
        row = {
            "fab_id": "FAB_A",
            "product_family": "F1",
            "source_step": "S1",
            "size_bin_left_um": size - 1.0,
            "size_bin_right_um": size + 1.0,
            "size_representative_um": size,
            "particle_count": 500 if index == 0 else 100,
            "inspection_area_mm2": 100.0,
            "measurement_period": "2026-09",
            "density_per_mm2": 5.0 if index == 0 else 1.0,
            "particle_type": "unknown",
        }
        row.update(overrides)
        rows.append(row)
    return rows


def _write_fab(tmp_path: Path, rows: list[dict], name: str = "fab.csv") -> Path:
    path = tmp_path / name
    pd.DataFrame(rows, columns=list(FAB_COLUMNS)).to_csv(path, index=False)
    return path


def _ca_frame(sizes=(2.0, 4.0)) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "source_step": ["S1"] * len(sizes),
            "particle_size_um": list(sizes),
            "pair_id": ["A__B"] * len(sizes),
            "failure_mode": ["lateral_short"] * len(sizes),
            "ca_nominal_um2": [2000.0] * len(sizes),
        }
    )


def test_expected_defect_count_hand_calculation() -> None:
    # 2000 um^2 = 0.002 mm^2; 0.002 * 5 = 0.01
    assert expected_defect_count(np.array([0.002]), np.array([5.0])) == pytest.approx(0.01)
    with pytest.raises(ValueError):
        expected_defect_count(np.array([1.0, 2.0]), np.array([1.0]))


def test_load_fab_table_fills_and_checks_density(tmp_path: Path) -> None:
    rows = _fab_rows()
    rows[0]["density_per_mm2"] = np.nan
    table = load_fab_table(_write_fab(tmp_path, rows), bin_method="piecewise_constant")
    assert float(table.frame["density_per_mm2"].iloc[0]) == pytest.approx(5.0)

    wrong = _fab_rows()
    wrong[1]["density_per_mm2"] = 9.9
    with pytest.raises(InputContractError, match="density_per_mm2"):
        load_fab_table(_write_fab(tmp_path, wrong, "wrong.csv"), bin_method="exact_point")


@pytest.mark.parametrize(
    ("mutation", "match"),
    [
        ({"size_bin_left_um": 9.0, "size_bin_right_um": 3.0}, "size_bin_left_um"),
        ({"size_representative_um": 99.0}, "size_representative_um"),
        ({"particle_count": -1}, "particle_count"),
        ({"inspection_area_mm2": 0.0}, "inspection_area_mm2"),
        ({"measurement_period": ""}, "measurement_period"),
    ],
)
def test_fab_table_rejects_bad_rows(tmp_path: Path, mutation: dict, match: str) -> None:
    rows = _fab_rows()
    rows[0].update(mutation)
    rows[0]["density_per_mm2"] = np.nan
    rows[1]["density_per_mm2"] = np.nan
    with pytest.raises(InputContractError, match=match):
        load_fab_table(_write_fab(tmp_path, rows, "bad.csv"), bin_method="piecewise_constant")


def test_fab_table_rejects_overlap_and_duplicates(tmp_path: Path) -> None:
    rows = _fab_rows()
    rows[1]["size_bin_left_um"] = 1.5  # overlaps the first bin (1.0 .. 3.0)
    with pytest.raises(InputContractError, match="겹칩"):
        load_fab_table(_write_fab(tmp_path, rows, "overlap.csv"), bin_method="piecewise_constant")
    duplicated = _fab_rows()
    duplicated[1] = dict(duplicated[0])
    # identical rows overlap as well, so either fail-closed message is acceptable
    with pytest.raises(InputContractError, match="중복|겹칩"):
        load_fab_table(_write_fab(tmp_path, duplicated, "dup.csv"), bin_method="exact_point")


def test_fab_table_rejects_unknown_bin_method_and_extension(tmp_path: Path) -> None:
    path = _write_fab(tmp_path, _fab_rows())
    with pytest.raises(InputContractError, match="bin_method"):
        load_fab_table(path, bin_method="guess")
    other = tmp_path / "fab.txt"
    other.write_text("x", encoding="utf-8")
    with pytest.raises(InputContractError, match="확장자"):
        load_fab_table(other, bin_method="exact_point")


def test_fab_density_lookup_prefers_exact_family(tmp_path: Path) -> None:
    rows = _fab_rows()
    generic = dict(rows[0])
    generic.update({"product_family": "", "size_bin_left_um": 5.0, "size_bin_right_um": 7.0,
                    "size_representative_um": 6.0, "particle_count": 100,
                    "density_per_mm2": 1.0})
    table = load_fab_table(
        _write_fab(tmp_path, [*rows, generic], "mixed.csv"), bin_method="piecewise_constant"
    )
    exact = fab_density_for(table, fab_id="FAB_A", product_family="F1", source_step="S1")
    assert set(exact["product_family"]) == {"F1"}
    fallback = fab_density_for(table, fab_id="FAB_A", product_family="OTHER", source_step="S1")
    assert set(fallback["product_family"]) == {"*"}
    with pytest.raises(InputContractError, match="source_step"):
        fab_density_for(table, fab_id="FAB_A", product_family="F1", source_step="S9")


@pytest.mark.parametrize("bin_method", ["exact_point", "piecewise_constant"])
def test_matching_and_weighting_hand_calculation(tmp_path: Path, bin_method: str) -> None:
    table = load_fab_table(_write_fab(tmp_path, _fab_rows()), bin_method=bin_method)
    ca = _ca_frame()
    weighted = apply_fab_weighting(
        ca, table, weighting_mode="actual", fab_id="FAB_A", product_family="F1",
        run_id="R", product_id="P",
    )
    assert tuple(weighted.columns) == WEIGHTED_COLUMNS
    assert len(weighted) == 2
    # 2000 um^2 -> 0.002 mm^2; densities 5 and 1 per mm^2
    assert weighted["expected_defect_count"].tolist() == pytest.approx([0.01, 0.002])
    assert ca.equals(_ca_frame())  # geometry input untouched


def test_piecewise_linear_interpolates_inside_range(tmp_path: Path) -> None:
    table = load_fab_table(_write_fab(tmp_path, _fab_rows()), bin_method="piecewise_linear")
    ca = _ca_frame(sizes=(3.0,))
    weighted = apply_fab_weighting(
        ca, table, weighting_mode="actual", fab_id="FAB_A", product_family="F1",
        run_id="R", product_id="P", require_full_coverage=False,
    )
    assert float(weighted["density_per_mm2"].iloc[0]) == pytest.approx(3.0)  # between 5 and 1
    outside = _ca_frame(sizes=(9.0,))
    matched, unmatched_sizes, _ = match_sizes(
        outside,
        fab_density_for(table, fab_id="FAB_A", product_family="F1", source_step="S1"),
        bin_method="piecewise_linear",
    )
    assert matched.empty and unmatched_sizes == [9.0]


def test_missing_coverage_blocks_by_default(tmp_path: Path) -> None:
    table = load_fab_table(_write_fab(tmp_path, _fab_rows()), bin_method="exact_point")
    ca = _ca_frame(sizes=(2.0, 4.0, 8.0))
    with pytest.raises(InputContractError, match="연결되지 않은"):
        apply_fab_weighting(
            ca, table, weighting_mode="actual", fab_id="FAB_A", product_family="F1",
            run_id="R", product_id="P",
        )
    partial = apply_fab_weighting(
        ca, table, weighting_mode="actual", fab_id="FAB_A", product_family="F1",
        run_id="R", product_id="P", require_full_coverage=False,
    )
    assert len(partial) == 2


def test_common_and_actual_differ(tmp_path: Path) -> None:
    rows = _fab_rows()
    reference = [dict(row, fab_id="FAB_REF", particle_count=row["particle_count"] * 2,
                      density_per_mm2=row["density_per_mm2"] * 2) for row in rows]
    table = load_fab_table(
        _write_fab(tmp_path, [*rows, *reference], "two.csv"), bin_method="piecewise_constant"
    )
    ca = _ca_frame()
    actual = apply_fab_weighting(
        ca, table, weighting_mode="actual", fab_id="FAB_A", product_family="F1",
        run_id="R", product_id="P",
    )
    common = apply_fab_weighting(
        ca, table, weighting_mode="common", fab_id="FAB_REF", product_family="F1",
        run_id="R", product_id="P",
    )
    assert common["expected_defect_count"].sum() == pytest.approx(
        2 * actual["expected_defect_count"].sum()
    )


def test_unsupported_mode_and_method(tmp_path: Path) -> None:
    path = _write_fab(tmp_path, _fab_rows())
    # the unsupported bin method is rejected as early as loading
    with pytest.raises(InputContractError, match="empirical_samples"):
        load_fab_table(path, bin_method="empirical_samples")
    table = load_fab_table(path, bin_method="piecewise_constant")
    with pytest.raises(InputContractError, match="weighting_mode"):
        apply_fab_weighting(
            _ca_frame(), table, weighting_mode="guess", fab_id="FAB_A", product_family="F1",
            run_id="R", product_id="P",
        )
    with pytest.raises(InputContractError, match="empirical_samples"):
        match_sizes(_ca_frame(), table.frame, bin_method="empirical_samples")
