"""Fab 이물 표 로더와 fail-closed 검증."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from particle_sim.errors import InputContractError

FAB_COLUMNS: tuple[str, ...] = (
    "fab_id",
    "product_family",
    "source_step",
    "size_bin_left_um",
    "size_bin_right_um",
    "size_representative_um",
    "particle_count",
    "inspection_area_mm2",
    "measurement_period",
    "density_per_mm2",
    "particle_type",
)

BIN_METHODS: tuple[str, ...] = (
    "exact_point",
    "piecewise_constant",
    "piecewise_linear",
    "empirical_samples",
)


@dataclass(frozen=True)
class FabTable:
    """검증된 Fab 이물 표."""

    frame: pd.DataFrame
    bin_method: str
    source: str


def load_fab_table(path: Path, *, bin_method: str) -> FabTable:
    """Fab 이물 표를 읽고 fail-closed 검증을 수행한다."""
    if bin_method == "empirical_samples":
        raise InputContractError(
            "empirical_samples는 아직 지원하지 않습니다(원본 표본 연결 규칙 미확정)",
            file=str(path),
            field="bin_method",
        )
    if bin_method not in BIN_METHODS:
        raise InputContractError(
            f"bin_method '{bin_method}'은(는) 지원되지 않습니다",
            file=str(path),
            field="bin_method",
        )

    suffix = path.suffix.lower()
    if suffix == ".parquet":
        try:
            frame = pd.read_parquet(path)
        except (OSError, ValueError) as exc:
            raise InputContractError(
                f"파일을 읽을 수 없습니다: {exc}", file=str(path), field="fab_data"
            ) from exc
    elif suffix == ".csv":
        try:
            frame = pd.read_csv(path)
        except (OSError, ValueError) as exc:
            raise InputContractError(
                f"파일을 읽을 수 없습니다: {exc}", file=str(path), field="fab_data"
            ) from exc
    else:
        raise InputContractError(
            f"지원하지 않는 확장자입니다: {suffix}", file=str(path), field="fab_data"
        )

    for col in FAB_COLUMNS:
        if col not in frame.columns:
            raise InputContractError(
                f"필수 열이 누락되었습니다: {col}", file=str(path), field=col
            )

    if frame.empty:
        raise InputContractError("표가 비어 있습니다", file=str(path), field="fab_data")

    frame = frame[list(FAB_COLUMNS)].copy()

    numeric_columns = (
        "size_bin_left_um", "size_bin_right_um", "size_representative_um",
        "particle_count", "inspection_area_mm2",
    )
    for column in numeric_columns:
        values = pd.to_numeric(frame[column], errors="coerce")
        if not np.isfinite(values).all():
            raise InputContractError(
                f"{column}에 비어 있거나 유한하지 않은 값이 있습니다",
                file=str(path),
                field=column,
            )
        frame[column] = values.astype(float)
    density = pd.to_numeric(frame["density_per_mm2"], errors="coerce")
    if np.isinf(density).any():
        raise InputContractError(
            "density_per_mm2에 유한하지 않은 값이 있습니다",
            file=str(path),
            field="density_per_mm2",
        )
    frame["density_per_mm2"] = density

    # product_family가 비어 있으면 "*"로 정규화
    # CSV의 빈 칸은 NaN으로 읽히므로 함께 정규화한다
    frame["product_family"] = (
        frame["product_family"].astype("object").where(frame["product_family"].notna(), "*")
    )
    frame["product_family"] = frame["product_family"].replace("", "*").astype(str)

    # density_per_mm2 자동 계산 또는 검증
    computed = frame["particle_count"] / frame["inspection_area_mm2"]
    missing = frame["density_per_mm2"].isna()
    # every value that IS given must match the count/area calculation, including mixed tables
    given = ~missing
    if given.any():
        rel_err = (frame.loc[given, "density_per_mm2"] - computed[given]).abs() / computed[
            given
        ].abs().clip(lower=1e-30)
        if (rel_err > 1e-6).any():
            raise InputContractError(
                "density_per_mm2 값이 particle_count / inspection_area_mm2와 일치하지 않습니다",
                file=str(path),
                field="density_per_mm2",
            )
    if missing.any():
        frame.loc[missing, "density_per_mm2"] = computed[missing]

    # 1. size_bin_left_um < size_bin_right_um, 둘 다 > 0
    if (frame["size_bin_left_um"] <= 0).any() or (frame["size_bin_right_um"] <= 0).any():
        raise InputContractError(
            "size_bin_left_um과 size_bin_right_um은 모두 0보다 커야 합니다",
            file=str(path),
            field="size_bin_left_um",
        )
    if (frame["size_bin_left_um"] >= frame["size_bin_right_um"]).any():
        raise InputContractError(
            "size_bin_left_um은 size_bin_right_um보다 작아야 합니다",
            file=str(path),
            field="size_bin_left_um",
        )

    # 2. size_bin_left_um <= size_representative_um <= size_bin_right_um
    if (
        (frame["size_representative_um"] < frame["size_bin_left_um"])
        | (frame["size_representative_um"] > frame["size_bin_right_um"])
    ).any():
        raise InputContractError(
            "size_representative_um은 bin 구간 내에 있어야 합니다",
            file=str(path),
            field="size_representative_um",
        )

    # 3. 같은 (fab_id, product_family, source_step) 안에서 bin 구간이 겹치지 않음
    for _, group in frame.groupby(["fab_id", "product_family", "source_step"]):
        intervals = sorted(
            zip(group["size_bin_left_um"], group["size_bin_right_um"], strict=True)
        )
        for i in range(len(intervals) - 1):
            # bins are half-open [left, right); touching bins are fine, overlaps are not
            if intervals[i][1] > intervals[i + 1][0]:
                raise InputContractError(
                    "같은 그룹 안에서 bin 구간이 겹칩니다",
                    file=str(path),
                    field="size_bin_left_um",
                )

    # 4. particle_count >= 0, inspection_area_mm2 > 0
    if (frame["particle_count"] < 0).any():
        raise InputContractError(
            "particle_count는 0 이상이어야 합니다",
            file=str(path),
            field="particle_count",
        )
    if (frame["inspection_area_mm2"] <= 0).any():
        raise InputContractError(
            "inspection_area_mm2는 0보다 커야 합니다",
            file=str(path),
            field="inspection_area_mm2",
        )

    # 5. measurement_period가 비어 있지 않음
    if frame["measurement_period"].isna().any() or (frame["measurement_period"] == "").any():
        raise InputContractError(
            "measurement_period가 비어 있을 수 없습니다",
            file=str(path),
            field="measurement_period",
        )

    # 7. 같은 (fab_id, product_family, source_step, size_bin_left_um) 행 중복
    dup_mask = frame.duplicated(
        subset=["fab_id", "product_family", "source_step", "size_bin_left_um"],
        keep=False,
    )
    if dup_mask.any():
        raise InputContractError(
            "같은 (fab_id, product_family, source_step, size_bin_left_um) 행이 중복됩니다",
            file=str(path),
            field="size_bin_left_um",
        )

    # deterministic order so that bin matching never depends on the input row order
    frame = frame.sort_values(
        ["fab_id", "product_family", "source_step", "size_bin_left_um"], kind="stable"
    ).reset_index(drop=True)
    return FabTable(frame=frame, bin_method=bin_method, source=str(path))


def fab_density_for(
    table: FabTable, *, fab_id: str, product_family: str, source_step: str
) -> pd.DataFrame:
    """fab_id와 source_step이 일치하는 행을 돌려준다.

    정확히 일치하는 product_family가 있으면 그 행만 쓰고 "*"는 무시한다.
    결과가 비면 InputContractError(field="source_step")를 발생시킨다.
    """
    base = table.frame[
        (table.frame["fab_id"] == fab_id) & (table.frame["source_step"] == source_step)
    ]
    exact = base[base["product_family"] == product_family]
    if not exact.empty:
        result = exact
    else:
        result = base[base["product_family"] == "*"]
    if result.empty:
        raise InputContractError(
            f"source_step '{source_step}'에 대한 Fab 분포가 없습니다",
            field="source_step",
        )
    return result
