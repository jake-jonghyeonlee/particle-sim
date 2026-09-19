"""Fab 밀도를 적용한 layer, 전극 조합, open 경로와 개별 전극 위험도 표."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from particle_sim.aggregation.critical_area import um2_to_mm2
from particle_sim.aggregation.fab_weighting import match_sizes
from particle_sim.errors import InputContractError
from particle_sim.io.fab import FabTable, fab_density_for_layer
from particle_sim.schema import ProductManifest

COMMON_COLUMNS = [
    "run_id", "product_id", "fab_id", "weighting_mode", "particle_source_step",
    "fab_distribution_layer", "particle_size_um", "size_bin_left_um", "size_bin_right_um",
    "mean_density_per_mm2", "fab_input_semantics", "particle_count_value",
    "inspection_area_mm2", "measurement_period", "particle_type", "sample_count",
    "fab_source_id", "fab_source_hash", "failure_unit_id", "failure_mode", "raw_ca_um2",
    "allocated_ca_um2", "unique_ca_um2", "raw_geometric_hit_expectation",
    "allocated_geometric_hit_expectation", "unique_geometric_hit_expectation",
    "allocation_semantics",
]
PAIR_COLUMNS = [
    *COMMON_COLUMNS, "electrode_pair_id", "electrode_a", "electrode_b",
    "electrode_a_layer", "electrode_b_layer",
]
OPEN_COLUMNS = [
    *COMMON_COLUMNS, "open_path_id", "electrode_id", "electrode_layer",
]
CONTACT_COLUMNS = [
    *COMMON_COLUMNS, "contact_link_id", "lower_electrode", "upper_electrode",
    "lower_electrode_layer", "upper_electrode_layer",
]
ELECTRODE_COLUMNS = [
    *COMMON_COLUMNS, "electrode_id", "electrode_layer",
    "involvement_geometric_hit_expectation", "allocated_electrode_contribution",
    "electrode_share",
]
LAYER_COLUMNS = [
    "run_id", "product_id", "fab_id", "weighting_mode", "particle_source_step",
    "fab_distribution_layer", "particle_size_um", "size_bin_left_um", "size_bin_right_um",
    "mean_density_per_mm2", "fab_input_semantics", "particle_count_value",
    "inspection_area_mm2", "measurement_period", "particle_type", "sample_count",
    "fab_source_id", "fab_source_hash", "layer_scope", "layer_id", "failure_mode", "ca_um2",
    "geometric_hit_expectation", "allocation_semantics",
]


OPEN_KINDS = ("open_path", "open_electrode")
UNIT_KINDS = ("electrode_pair", *OPEN_KINDS, "contact")
SHARE_SEMANTICS = {
    "electrode_pair": "short_equal_share_between_electrodes",
    "open_path": "open_full_contribution",
    "open_electrode": "open_full_contribution",
    "contact": "cno_equal_share_between_electrodes",
}


@dataclass(frozen=True)
class RiskTables:
    by_layer: pd.DataFrame
    by_pair: pd.DataFrame
    by_open_path: pd.DataFrame
    by_contact: pd.DataFrame
    by_electrode: pd.DataFrame
    summary: pd.DataFrame


def _first_value(row: pd.Series, names: tuple[str, ...]) -> float:
    for name in names:
        if name in row and pd.notna(row[name]):
            return float(row[name])
    raise InputContractError(f"CA 값이 없습니다: {names}", field=names[0])


def _density_by_size(
    rows: pd.DataFrame,
    table: FabTable,
    *,
    fab_id: str,
    product_family: str,
    fab_distribution_layer: str,
) -> dict[float, tuple[float, pd.Series]]:
    density = fab_density_for_layer(
        table,
        fab_id=fab_id,
        product_family=product_family,
        fab_distribution_layer=fab_distribution_layer,
    )
    matched, missing_sizes, missing_bins = match_sizes(rows, density, bin_method=table.bin_method)
    if missing_sizes or missing_bins:
        raise InputContractError(
            "Fab 크기 bin 연결이 완전하지 않습니다: "
            f"연결되지 않은 크기={missing_sizes}, 연결되지 않은 bin={missing_bins}",
            field="particle_size_um",
        )
    result: dict[float, tuple[float, pd.Series]] = {}
    for item in matched.to_dict("records"):
        size = float(item["row"]["particle_size_um"])
        if size in result:
            continue
        result[size] = (float(item["density_per_mm2"]), item["bin"])
    return result


def _weighted(ca_um2: float, density_per_mm2: float) -> float:
    return um2_to_mm2(ca_um2) * density_per_mm2


def _fab_metadata(bin_row: pd.Series, table: FabTable) -> dict:
    """정규화된 Fab 행에서도 원래 측정 의미와 출처를 잃지 않는다."""
    def optional(name: str) -> object:
        value = bin_row.get(name)
        return None if pd.isna(value) else value

    return {
        "fab_input_semantics": table.input_semantics,
        "particle_count_value": float(bin_row["particle_count"]),
        "inspection_area_mm2": float(bin_row["inspection_area_mm2"]),
        "measurement_period": str(bin_row["measurement_period"]),
        "particle_type": str(bin_row["particle_type"]),
        "sample_count": optional("sample_count"),
        "fab_source_id": optional("source_id"),
        "fab_source_hash": optional("source_hash"),
    }


def build_risk_tables(
    ca: pd.DataFrame,
    union: pd.DataFrame,
    table: FabTable,
    manifest: ProductManifest,
    *,
    fab_id: str,
    weighting_mode: str,
) -> RiskTables:
    """한 Fab 분포에 대한 중복 없는 위험도 표를 만든다."""
    if table.input_semantics == "mean_per_sample" and table.bin_method != "exact_point":
        raise InputContractError(
            "평균 이물 개수 입력은 exact_point 연결만 지원합니다", field="bin_method"
        )
    required_ca = {
        "run_id",
        "source_step",
        "fab_distribution_layer",
        "particle_size_um",
        "failure_unit_id",
        "unit_kind",
        "electrode_pair_id",
        "open_path_id",
        "electrode_a",
        "electrode_b",
        "electrode_a_layer",
        "electrode_b_layer",
        "failure_mode",
    }
    missing_ca = sorted(required_ca - set(ca.columns))
    if missing_ca:
        raise InputContractError(
            f"schema 2.0 ca.parquet 열이 누락됐습니다: {missing_ca}", field=missing_ca[0]
        )
    unknown_kinds = sorted(set(ca["unit_kind"].dropna()) - set(UNIT_KINDS))
    if unknown_kinds:
        raise InputContractError(
            f"알 수 없는 unit_kind가 있습니다: {unknown_kinds}", field="unit_kind"
        )
    ca_key = ["source_step", "fab_distribution_layer", "particle_size_um", "failure_unit_id"]
    if ca.duplicated(ca_key).any():
        raise InputContractError("CA 표에 같은 불량 단위 조건이 중복됩니다", field="failure_unit_id")
    required_union = {
        "run_id",
        "source_step",
        "fab_distribution_layer",
        "particle_size_um",
        "failure_mode",
    }
    missing_union = sorted(required_union - set(union.columns))
    if missing_union:
        raise InputContractError(
            f"schema 2.0 ca_union.parquet 열이 누락됐습니다: {missing_union}",
            field=missing_union[0],
        )
    any_union = union[union["failure_mode"] == "any"]
    union_key = ["source_step", "fab_distribution_layer", "particle_size_um"]
    expected_union = ca[union_key].drop_duplicates()
    actual_union = any_union[union_key].drop_duplicates()
    if any_union.duplicated(union_key).any() or set(
        map(tuple, expected_union.to_numpy())
    ) != set(map(tuple, actual_union.to_numpy())):
        raise InputContractError(
            "ca_union.parquet의 failure_mode=any 조건이 CA 표와 1:1로 맞지 않습니다",
            field="failure_mode",
        )

    product_family = manifest.product_family
    pair_rows: list[dict] = []
    open_rows: list[dict] = []
    contact_rows: list[dict] = []
    electrode_rows: list[dict] = []
    affected_layer_rows: list[dict] = []
    layer_rows: list[dict] = []

    for (source_step, fab_layer), ca_group in ca.groupby(
        ["source_step", "fab_distribution_layer"], sort=False
    ):
        density_map = _density_by_size(
            ca_group,
            table,
            fab_id=fab_id,
            product_family=product_family,
            fab_distribution_layer=str(fab_layer),
        )
        for _, row in ca_group.iterrows():
            size = float(row["particle_size_um"])
            density_value, bin_row = density_map[size]
            raw_ca = _first_value(row, ("ca_mean_um2", "ca_nominal_um2"))
            allocated_ca = _first_value(
                row, ("ca_allocated_mean_um2", "ca_allocated_nominal_um2")
            )
            unique_ca = _first_value(row, ("ca_unique_mean_um2", "ca_unique_nominal_um2"))
            base = {
                "run_id": str(row["run_id"]),
                "product_id": manifest.product_id,
                "fab_id": fab_id,
                "weighting_mode": weighting_mode,
                "particle_source_step": str(source_step),
                "fab_distribution_layer": str(fab_layer),
                "particle_size_um": size,
                "size_bin_left_um": float(bin_row["size_bin_left_um"]),
                "size_bin_right_um": float(bin_row["size_bin_right_um"]),
                "mean_density_per_mm2": density_value,
                **_fab_metadata(bin_row, table),
                "failure_unit_id": str(row["failure_unit_id"]),
                "failure_mode": str(row["failure_mode"]),
                "raw_ca_um2": raw_ca,
                "allocated_ca_um2": allocated_ca,
                "unique_ca_um2": unique_ca,
                "raw_geometric_hit_expectation": _weighted(raw_ca, density_value),
                "allocated_geometric_hit_expectation": _weighted(
                    allocated_ca, density_value
                ),
                "unique_geometric_hit_expectation": _weighted(unique_ca, density_value),
                "allocation_semantics": "pixel_equal_share_1_over_k",
            }
            if row["unit_kind"] == "contact":
                contact_rows.append(
                    {
                        **base,
                        "contact_link_id": str(row["failure_unit_id"]),
                        "lower_electrode": str(row["electrode_a"]),
                        "upper_electrode": str(row["electrode_b"]),
                        "lower_electrode_layer": str(row["electrode_a_layer"]),
                        "upper_electrode_layer": str(row["electrode_b_layer"]),
                    }
                )
                participants = [
                    (str(row["electrode_a"]), str(row["electrode_a_layer"]), 0.5),
                    (str(row["electrode_b"]), str(row["electrode_b_layer"]), 0.5),
                ]
            elif row["unit_kind"] in OPEN_KINDS:
                open_rows.append(
                    {
                        **base,
                        "open_path_id": str(row["open_path_id"]),
                        "electrode_id": str(row["electrode_a"]),
                        "electrode_layer": str(row["electrode_a_layer"]),
                    }
                )
                participants = [(str(row["electrode_a"]), str(row["electrode_a_layer"]), 1.0)]
            else:
                pair_rows.append(
                    {
                        **base,
                        "electrode_pair_id": str(row["electrode_pair_id"]),
                        "electrode_a": str(row["electrode_a"]),
                        "electrode_b": str(row["electrode_b"]),
                        "electrode_a_layer": str(row["electrode_a_layer"]),
                        "electrode_b_layer": str(row["electrode_b_layer"]),
                    }
                )
                participants = [
                    (str(row["electrode_a"]), str(row["electrode_a_layer"]), 0.5),
                    (str(row["electrode_b"]), str(row["electrode_b_layer"]), 0.5),
                ]
            for electrode_id, layer_name, share in participants:
                electrode_rows.append(
                    {
                        **base,
                        "electrode_id": electrode_id,
                        "electrode_layer": layer_name,
                        "involvement_geometric_hit_expectation": base[
                            "allocated_geometric_hit_expectation"
                        ],
                        "allocated_electrode_contribution": base[
                            "allocated_geometric_hit_expectation"
                        ]
                        * share,
                        "electrode_share": share,
                    }
                )
                affected_layer_rows.append(
                    {
                        "run_id": base["run_id"],
                        "product_id": base["product_id"],
                        "fab_id": base["fab_id"],
                        "weighting_mode": base["weighting_mode"],
                        "particle_source_step": base["particle_source_step"],
                        "fab_distribution_layer": base["fab_distribution_layer"],
                        "particle_size_um": base["particle_size_um"],
                        "size_bin_left_um": base["size_bin_left_um"],
                        "size_bin_right_um": base["size_bin_right_um"],
                        "mean_density_per_mm2": base["mean_density_per_mm2"],
                        **{
                            name: base[name]
                            for name in (
                                "fab_input_semantics", "particle_count_value",
                                "inspection_area_mm2", "measurement_period", "particle_type",
                                "sample_count", "fab_source_id", "fab_source_hash",
                            )
                        },
                        "layer_scope": "affected_electrode_layer",
                        "layer_id": layer_name,
                        "failure_mode": base["failure_mode"],
                        "ca_um2": allocated_ca * share,
                        "geometric_hit_expectation": base[
                            "allocated_geometric_hit_expectation"
                        ]
                        * share,
                        "allocation_semantics": SHARE_SEMANTICS[str(row["unit_kind"])],
                    }
                )

        union_group = union[
            (union["source_step"] == source_step)
            & (union["fab_distribution_layer"] == fab_layer)
            & (union["failure_mode"] == "any")
        ]
        for _, row in union_group.iterrows():
            size = float(row["particle_size_um"])
            density_value, bin_row = density_map[size]
            union_ca = _first_value(row, ("ca_union_mean_um2", "ca_union_nominal_um2"))
            layer_rows.append(
                {
                    "run_id": str(row["run_id"]),
                    "product_id": manifest.product_id,
                    "fab_id": fab_id,
                    "weighting_mode": weighting_mode,
                    "particle_source_step": str(source_step),
                    "fab_distribution_layer": str(fab_layer),
                    "particle_size_um": size,
                    "size_bin_left_um": float(bin_row["size_bin_left_um"]),
                    "size_bin_right_um": float(bin_row["size_bin_right_um"]),
                    "mean_density_per_mm2": density_value,
                    **_fab_metadata(bin_row, table),
                    "layer_scope": "particle_source_layer",
                    "layer_id": str(fab_layer),
                    "failure_mode": "any",
                    "ca_um2": union_ca,
                    "geometric_hit_expectation": _weighted(union_ca, density_value),
                    "allocation_semantics": "union_any",
                }
            )

    pair_frame = pd.DataFrame(pair_rows, columns=PAIR_COLUMNS)
    open_frame = pd.DataFrame(open_rows, columns=OPEN_COLUMNS)
    contact_frame = pd.DataFrame(contact_rows, columns=CONTACT_COLUMNS)
    electrode_frame = pd.DataFrame(electrode_rows, columns=ELECTRODE_COLUMNS)
    layer_frame = pd.concat(
        [
            pd.DataFrame(layer_rows, columns=LAYER_COLUMNS),
            pd.DataFrame(affected_layer_rows, columns=LAYER_COLUMNS),
        ],
        ignore_index=True,
        sort=False,
    )
    summary_frames: list[pd.DataFrame] = []
    identity = ["run_id", "product_id", "fab_id", "weighting_mode"]
    definitions = (
        ("layer", layer_frame, ["layer_scope", "layer_id"], "geometric_hit_expectation"),
        (
            "electrode_pair",
            pair_frame,
            ["fab_distribution_layer", "electrode_pair_id", "failure_mode"],
            "allocated_geometric_hit_expectation",
        ),
        (
            "open_path",
            open_frame,
            ["fab_distribution_layer", "open_path_id"],
            "allocated_geometric_hit_expectation",
        ),
        (
            "contact",
            contact_frame,
            ["fab_distribution_layer", "contact_link_id"],
            "allocated_geometric_hit_expectation",
        ),
        (
            "electrode",
            electrode_frame,
            ["fab_distribution_layer", "electrode_id"],
            "allocated_electrode_contribution",
        ),
    )
    for level, frame, keys, value in definitions:
        if frame.empty:
            continue
        grouped = frame.groupby([*identity, *keys], dropna=False, as_index=False)[value].sum()
        grouped = grouped.sort_values(value, ascending=False, kind="stable").reset_index(drop=True)
        grouped.insert(0, "rank", np.arange(1, len(grouped) + 1))
        grouped.insert(0, "aggregation_level", level)
        grouped = grouped.rename(columns={value: "geometric_hit_expectation"})
        summary_frames.append(grouped)
    summary = pd.concat(summary_frames, ignore_index=True, sort=False) if summary_frames else pd.DataFrame()
    return RiskTables(
        by_layer=layer_frame,
        by_pair=pair_frame,
        by_open_path=open_frame,
        by_contact=contact_frame,
        by_electrode=electrode_frame,
        summary=summary,
    )
