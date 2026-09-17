"""Weight geometry CA with a Fab particle size distribution.

The result is an expected defect count derived from CA and measured particle density, not an
actual failure probability. Geometry results are never modified here.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from particle_sim.errors import InputContractError
from particle_sim.io.fab import FabTable, fab_density_for

WEIGHTED_COLUMNS: tuple[str, ...] = (
    "run_id", "product_id", "fab_id", "product_family", "weighting_mode", "source_step",
    "failure_mode", "pair_id", "size_bin_left_um", "size_bin_right_um",
    "size_representative_um", "particle_size_um", "ca_um2", "ca_mm2", "density_per_mm2",
    "expected_defect_count", "bin_method",
)
WEIGHTING_MODES: tuple[str, ...] = ("common", "actual")


def expected_defect_count(ca_mm2: np.ndarray, density_per_mm2: np.ndarray) -> float:
    """Sum of CA (mm^2) times particle density (1/mm^2)."""
    ca = np.asarray(ca_mm2, dtype=np.float64)
    density = np.asarray(density_per_mm2, dtype=np.float64)
    if ca.shape != density.shape:
        raise ValueError(f"길이가 다릅니다: {ca.shape} vs {density.shape}")
    return float(np.sum(ca * density))


def match_sizes(
    ca: pd.DataFrame, density: pd.DataFrame, *, bin_method: str
) -> tuple[pd.DataFrame, list[float], list[float]]:
    """Join CA rows with Fab bins; also report the sizes and bins that stayed unmatched."""
    if bin_method == "empirical_samples":
        raise InputContractError(
            "empirical_samples는 이 단계에서 지원하지 않습니다(원본 표본 연결 규칙 미확정)",
            field="bin_method",
        )
    sizes = sorted({float(v) for v in ca["particle_size_um"]})
    matched_rows: list[dict] = []
    matched_sizes: set[float] = set()
    matched_bins: set[float] = set()

    if bin_method == "piecewise_linear":
        points = density.sort_values("size_representative_um")
        xs = points["size_representative_um"].to_numpy(dtype=np.float64)
        ys = points["density_per_mm2"].to_numpy(dtype=np.float64)
        for _, row in ca.iterrows():
            size = float(row["particle_size_um"])
            if xs.size == 0 or size < xs.min() or size > xs.max():
                continue  # never extrapolate outside the observed range
            value = float(np.interp(size, xs, ys))
            nearest = points.iloc[int(np.argmin(np.abs(xs - size)))]
            matched_rows.append(
                {
                    "row": row,
                    "bin": nearest,
                    "density_per_mm2": value,
                }
            )
            matched_sizes.add(size)
            matched_bins.add(float(nearest["size_bin_left_um"]))
    else:
        for _, row in ca.iterrows():
            size = float(row["particle_size_um"])
            for _, bin_row in density.iterrows():
                if bin_method == "exact_point":
                    hit = np.isclose(size, float(bin_row["size_representative_um"]))
                else:  # piecewise_constant
                    hit = (
                        float(bin_row["size_bin_left_um"])
                        <= size
                        <= float(bin_row["size_bin_right_um"])
                    )
                if hit:
                    matched_rows.append(
                        {
                            "row": row,
                            "bin": bin_row,
                            "density_per_mm2": float(bin_row["density_per_mm2"]),
                        }
                    )
                    matched_sizes.add(size)
                    matched_bins.add(float(bin_row["size_bin_left_um"]))
                    break

    unmatched_sizes = sorted(set(sizes) - matched_sizes)
    unmatched_bins = sorted(
        {float(v) for v in density["size_bin_left_um"]} - matched_bins
    )
    return pd.DataFrame(matched_rows), unmatched_sizes, unmatched_bins


def _ca_column(ca: pd.DataFrame) -> str:
    for name in ("ca_mean_um2", "ca_nominal_um2"):
        if name in ca.columns and ca[name].notna().any():
            return name
    raise InputContractError(
        "CA 표에 ca_mean_um2 또는 ca_nominal_um2 값이 없습니다", field="ca_um2"
    )


def apply_fab_weighting(
    ca: pd.DataFrame,
    table: FabTable,
    *,
    weighting_mode: str,
    fab_id: str,
    product_family: str,
    run_id: str,
    product_id: str,
    require_full_coverage: bool = True,
) -> pd.DataFrame:
    """Multiply every CA row with the matching Fab density (common or actual weighting)."""
    if weighting_mode not in WEIGHTING_MODES:
        raise InputContractError(
            f"weighting_mode '{weighting_mode}'은(는) common 또는 actual이어야 합니다",
            field="weighting_mode",
        )
    for column in ("source_step", "particle_size_um", "pair_id", "failure_mode"):
        if column not in ca.columns:
            raise InputContractError(
                f"CA 표에 필수 열이 없습니다: {column}", field=column
            )
    value_column = _ca_column(ca)

    rows: list[dict] = []
    missing_sizes: dict[str, list[float]] = {}
    missing_bins: dict[str, list[float]] = {}
    for source_step, ca_step in ca.groupby("source_step", sort=False):
        density = fab_density_for(
            table, fab_id=fab_id, product_family=product_family, source_step=str(source_step)
        )
        matched, unmatched_sizes, unmatched_bins = match_sizes(
            ca_step, density, bin_method=table.bin_method
        )
        if unmatched_sizes:
            missing_sizes[str(source_step)] = unmatched_sizes
        if unmatched_bins:
            missing_bins[str(source_step)] = unmatched_bins
        for entry in matched.to_dict("records"):
            ca_row = entry["row"]
            bin_row = entry["bin"]
            ca_um2 = float(ca_row[value_column])
            ca_mm2 = ca_um2 * 1e-6
            density_value = float(entry["density_per_mm2"])
            rows.append(
                {
                    "run_id": run_id,
                    "product_id": product_id,
                    "fab_id": fab_id,
                    "product_family": product_family,
                    "weighting_mode": weighting_mode,
                    "source_step": str(source_step),
                    "failure_mode": str(ca_row["failure_mode"]),
                    "pair_id": str(ca_row["pair_id"]),
                    "size_bin_left_um": float(bin_row["size_bin_left_um"]),
                    "size_bin_right_um": float(bin_row["size_bin_right_um"]),
                    "size_representative_um": float(bin_row["size_representative_um"]),
                    "particle_size_um": float(ca_row["particle_size_um"]),
                    "ca_um2": ca_um2,
                    "ca_mm2": ca_mm2,
                    "density_per_mm2": density_value,
                    "expected_defect_count": ca_mm2 * density_value,
                    "bin_method": table.bin_method,
                }
            )

    if require_full_coverage and (missing_sizes or missing_bins):
        raise InputContractError(
            "Fab 크기 bin 연결이 완전하지 않습니다: "
            f"연결되지 않은 크기={missing_sizes}, 연결되지 않은 bin={missing_bins}",
            product_id=product_id,
            field="particle_size_um",
        )
    return pd.DataFrame(rows, columns=list(WEIGHTED_COLUMNS))
