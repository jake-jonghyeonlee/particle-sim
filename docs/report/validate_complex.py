"""Validate the persisted complex fixture run and write report evidence."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import zarr

ROOT = Path(__file__).resolve().parent
RUN = ROOT / "work/complex-run"
STORE = zarr.open_group(RUN / "hotspots.zarr", mode="r")

steps = list(STORE.attrs["source_steps"])
sizes = [float(value) for value in STORE.attrs["particle_sizes_um"]]
pairs = list(STORE.attrs["pair_ids"])
modes = list(STORE.attrs["pair_modes"])
union_keys = list(STORE.attrs["union_keys"])
nominal = np.asarray(STORE["nominal_failure"]).astype(bool)
unions = np.asarray(STORE["union_failure"]).astype(bool)
pair_counts = np.asarray(STORE["union_pair_count"])

monotonic_violations: list[dict[str, object]] = []
for step_index, step in enumerate(steps):
    for pair_index, pair in enumerate(pairs):
        for size_index in range(1, len(sizes)):
            lost = nominal[step_index, size_index - 1, pair_index] & ~nominal[
                step_index, size_index, pair_index
            ]
            if lost.any():
                monotonic_violations.append(
                    {
                        "source_step": step,
                        "pair_id": pair,
                        "smaller_um": sizes[size_index - 1],
                        "larger_um": sizes[size_index],
                        "lost_pixels": int(lost.sum()),
                    }
                )

union_mismatches: list[dict[str, object]] = []
for step_index, step in enumerate(steps):
    for size_index, size in enumerate(sizes):
        for union_index, mode in enumerate(union_keys):
            members = [
                nominal[step_index, size_index, pair_index]
                for pair_index, pair_mode in enumerate(modes)
                if mode == "any" or pair_mode == mode
            ]
            expected_union = np.logical_or.reduce(members)
            expected_count = sum(member.astype(np.uint16) for member in members)
            if not np.array_equal(unions[step_index, size_index, union_index], expected_union):
                union_mismatches.append({"source_step": step, "size_um": size, "mode": mode})
            if not np.array_equal(
                pair_counts[step_index, size_index, union_index], expected_count
            ):
                union_mismatches.append(
                    {"source_step": step, "size_um": size, "mode": f"{mode}:pair_count"}
                )

ca = pd.read_parquet(RUN / "ca.parquet")
pixel_area_um2 = 0.25
ca_max_error_um2 = 0.0
for row in ca.itertuples(index=False):
    index = (
        steps.index(row.source_step),
        sizes.index(float(row.particle_size_um)),
        pairs.index(row.pair_id),
    )
    expected = float(nominal[index].sum()) * pixel_area_um2
    ca_max_error_um2 = max(ca_max_error_um2, abs(float(row.ca_nominal_um2) - expected))

step_index = steps.index("S_M1")
size_index = sizes.index(10.0)
pixel_counts = {
    pair: int(nominal[step_index, size_index, pair_index].sum())
    for pair_index, pair in enumerate(pairs)
}
expected_positive = [
    "G1__D1", "G2__D1", "G2__D3", "G3__D2",
    "D1__P1", "D2__P3", "D3__P2", "D4__P3",
]
negative_controls = ["G1__D4", "G3__D3"]
legacy = pd.read_csv(RUN / "legacy_vs_edt.csv")

result = {
    "shape_yx": list(nominal.shape[-2:]),
    "electrode_count": 10,
    "pair_count": len(pairs),
    "source_steps": steps,
    "particle_sizes_um": sizes,
    "expected_positive_s_m1_d10": {
        pair: pixel_counts[pair] for pair in expected_positive
    },
    "negative_controls_s_m1_d10": {
        pair: pixel_counts[pair] for pair in negative_controls
    },
    "max_simultaneous_pairs": int(pair_counts.max()),
    "monotonic_violations": monotonic_violations,
    "union_mismatches": union_mismatches,
    "ca_max_error_um2": ca_max_error_um2,
    "legacy_difference_pixels": int(legacy["n_diff"].sum()),
    "legacy_unexplained_pixels": int(legacy["unexplained"].sum()),
    "external_existing_sim_compared": False,
    "external_comparison_reason": (
        "기존 현장 sim 실행 파일과 동일 입력 변환기가 저장소에 없어 내부 legacy/EDT만 대조"
    ),
}

assert all(result["expected_positive_s_m1_d10"].values())
assert not any(result["negative_controls_s_m1_d10"].values())
assert result["max_simultaneous_pairs"] >= 3
assert not monotonic_violations
assert not union_mismatches
assert ca_max_error_um2 == 0.0
assert result["legacy_unexplained_pixels"] == 0

output = ROOT / "evidence/complex-validation.json"
output.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(result, ensure_ascii=False, indent=2))
