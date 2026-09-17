# arch: ca-outputs (T7)

Files: `src/particle_sim/aggregation/critical_area.py`, `src/particle_sim/io/outputs.py`,
tests `tests/unit/test_ca_outputs.py`. Depends on T2 `GridSpec(shape_yx, pixel_size_y_um,
pixel_size_x_um, origin)` with `pixel_area_um2`, `roi_area_um2`. Pure numeric/IO module:
does not import failure/ or geometry/.

## aggregation/critical_area.py
```python
def critical_area_um2(hotspot: np.ndarray, grid: GridSpec) -> float
# bool or float map in [0, 1]; shape must equal grid.shape_yx; NaN or out of range -> ValueError
# float(hotspot.astype(float64).sum() * grid.pixel_area_um2)
def um2_to_mm2(value_um2: float) -> float          # * 1e-6
def normalized_ca(ca_um2: float, grid: GridSpec) -> float   # ca / grid.roi_area_um2 (unit "1")
```

## io/outputs.py
```python
SCHEMA_VERSION = "1.0"
CA_COLUMNS: tuple[str, ...] = (
    "schema_version", "run_id", "product_id", "design_version", "source_step",
    "fab_distribution_layer", "particle_size_um", "pair_id", "electrode_a", "electrode_b",
    "net_a", "net_b", "failure_mode", "classification", "ca_legacy_um2", "ca_nominal_um2",
    "ca_nominal_mm2", "ca_nominal_normalized", "ca_mean_um2", "ca_p05_um2", "ca_p50_um2",
    "ca_p95_um2", "scenario_count", "convergence_delta", "assumption_set", "simulator_version")
UNION_CA_COLUMNS: tuple[str, ...] = (
    "schema_version", "run_id", "product_id", "design_version", "source_step",
    "fab_distribution_layer", "particle_size_um", "failure_mode", "pair_ids",
    "ca_union_legacy_um2", "ca_union_nominal_um2", "ca_union_nominal_normalized",
    "ca_pair_sum_nominal_um2", "assumption_set", "simulator_version")
def write_table(rows: list[dict], columns: tuple[str, ...], path: Path) -> None
```
Each row keys must equal set(columns) exactly (else ValueError naming missing/extra keys);
DataFrame in column order, `to_parquet(path, index=False, engine="pyarrow")`.
Nominal runs put NaN in ca_mean/p05/p50/p95 and convergence_delta and 0 in scenario_count.

```python
class HotspotStore:
    def __init__(self, path: Path, *, source_steps: list[str], particle_sizes_um: list[float],
                 pair_ids: list[str], pair_modes: list[str], union_keys: list[str],
                 shape_yx: tuple[int, int]) -> None
```
Creates a NEW zarr v3 group (`zarr.open_group(path, mode="w-")`; existing path -> error).
Arrays (chunks `(1, 1, 1, min(512, Y), min(512, X))`):
`legacy_failure` uint8 [S, D, P, Y, X] fill 0; `nominal_failure` uint8 [S, D, P, Y, X] fill 0;
`pair_margin_um` float32 [S, D, P, Y, X] fill NaN; `union_failure` uint8 [S, D, U, Y, X];
`union_pair_count` uint16 [S, D, U, Y, X]. Group attrs: `axes` = ["source_step",
"particle_size_um", "pair_or_union", "y", "x"], `source_steps`, `particle_sizes_um`,
`pair_ids`, `pair_modes`, `union_keys`, `units` {"length": "um", "margin": "um"},
`value_semantics` "binary geometric contact risk (not failure probability)".
```python
    def write_pair(self, step_index: int, size_index: int, pair_index: int, *,
                   legacy: np.ndarray, nominal: np.ndarray, margin_um: np.ndarray) -> None
    def write_union(self, step_index: int, size_index: int, union_index: int, *,
                    any_failure: np.ndarray, pair_count: np.ndarray) -> None
def open_hotspot_store(path: Path) -> zarr.Group     # mode="r"
```
Shape mismatch or index out of range -> ValueError. Margin -inf is stored as-is.

```python
def make_run_id(product_id: str, now: datetime) -> str   # f"{now:%Y%m%dT%H%M%SZ}_{product_id}" (UTC)
def build_run_manifest(*, run_id: str, product_id: str, design_version: str,
                       config_hash: str, input_hashes: dict[str, str],
                       particle_sizes_um: list[float], source_steps: list[str],
                       assumption_set: str, warnings: list[str],
                       legacy_kernel: dict, post_step_policy: dict,
                       created_at: datetime, git_commit: str | None = None) -> dict
# adds "schema_version", "simulator_version" (particle_sim.__version__),
# "created_at" ISO-8601 UTC, "units": {"length": "um", "area": "um2"},
# "random_seed": None, "sobol_exponent_m": None
def write_json(data: dict, path: Path) -> None   # indent=2, ensure_ascii=False, sort_keys=True
def write_float_tiff(array: np.ndarray, path: Path) -> None   # tifffile.imwrite float32
```
Tests: CA of all-ones 10x10 (pixel 0.5 -> 25 um^2), 0.5 map, mm^2, normalized; NaN
rejected; write_table round trip + missing/extra key errors; HotspotStore shapes, dtypes,
attrs, write/read back, existing path error; run manifest keys; tiff round trip dtype.
