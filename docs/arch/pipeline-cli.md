# arch: pipeline-cli (T8)

Files: `src/particle_sim/pipeline.py`, `src/particle_sim/cli.py` (add one subcommand),
tests `tests/unit/test_pipeline_cli.py`. Depends on T1-T7 public APIs:
`load_manifest`, `canonical_hash`, `file_sha256` (config); `todo_physics_warnings` (schema);
`load_product_masks -> MaskSet(grid, layers, electrodes, input_hashes)` (io.masks);
`build_stack_z` (io.stack); `legacy_short_hotspot(mask_a, mask_b, d, grid, kernel)` (geometry.legacy);
`explain_legacy_vs_edt(mask_a, mask_b, d, grid, kernel) -> PixelComparison(equal, n_diff,
only_in_reference, only_in_candidate, cause_counts, cause_map)` (validation.regression);
`DistanceCache(masks)` (geometry.contact); `simulate_pairs(manifest, masks, stack_z, step_id,
d, cache) -> list[PairResult]`, `union_by_mode(results) -> dict[str, UnionResult]`
(failure.union; PairResult fields pair_id, failure_mode, electrode_a, electrode_b, failure,
contact_a, contact_b, pair_margin_um, classification; UnionResult fields failure_mode,
pair_ids, any_failure, pair_count); `critical_area_um2`, `um2_to_mm2`, `normalized_ca`
(aggregation.critical_area); `CA_COLUMNS`, `UNION_CA_COLUMNS`, `write_table`, `HotspotStore`,
`make_run_id`, `build_run_manifest`, `write_json`, `SCHEMA_VERSION` (io.outputs).
cli.py already has `build_parser()` with subparsers (dest "command", `func` handler) and
`validate-input`, `make-synthetic`; InputContractError -> stderr `입력 오류: ...`, exit 2.

## pipeline.py
```python
@dataclass(frozen=True)
class NominalRunResult:
    run_id: str
    output_dir: Path
    hotspots: Path          # output_dir / "hotspots.zarr"
    ca_table: Path          # output_dir / "ca.parquet"
    union_table: Path       # output_dir / "ca_union.parquet"
    legacy_comparison: Path # output_dir / "legacy_vs_edt.csv"
    run_manifest: Path      # output_dir / "run_manifest.json"

def run_nominal(manifest_path: Path, output_dir: Path, *,
                now: datetime | None = None) -> NominalRunResult
```
1. `output_dir` exists and is non-empty -> InputContractError(field="output") (never overwrite).
2. Load manifest, masks, stack_z. `now` default `datetime.now(UTC)`.
3. union_keys = modes in first-appearance order of electrode_pairs + ["any"].
4. For step s, size d (manifest order): `results = simulate_pairs(...)`; for each pair:
   legacy = legacy_short_hotspot(electrode masks, d); store legacy/nominal/margin;
   CA row (ca_legacy_um2 from legacy, nominal from failure, net ids from manifest,
   assumption_set = post_step_policy.mode, probability columns NaN, scenario_count 0).
   unions = union_by_mode(results); legacy union = OR of legacy maps of the same pair_ids;
   union row with pair_ids joined by ";" and ca_pair_sum_nominal_um2 = sum of pair CAs.
5. Legacy comparison is step independent: for each size and pair once, one CSV row
   `particle_size_um, pair_id, equal, n_diff, only_in_reference, only_in_candidate,
   kernel_quantization, kernel_shape, unexplained`.
6. Run manifest: config_hash = canonical_hash(manifest.model_dump(mode="json")),
   input_hashes = masks.input_hashes plus {"manifest": file_sha256(manifest_path)},
   warnings = todo_physics_warnings(manifest), git_commit None.

## cli.py addition
`simulate-nominal --product-config PATH --output DIR` -> run_nominal; prints JSON
`{"status": "ok", "run_id", "outputs": {name: path}}`; return 0.

Tests (synthetic `pilot` via `build_case` + `write_case` in tmp_path):
- run succeeds; zarr `nominal_failure` shape (3, 5, 3, 200, 200); union arrays U == 3
  (lateral_short, vertical_short, any).
- ca.parquet columns == CA_COLUMNS, 45 rows; ca_nominal_um2 == nominal hotspot sum * 0.25.
- union rows: ca_union_nominal_um2 <= ca_pair_sum_nominal_um2 for all rows.
- legacy_vs_edt.csv has 15 rows and unexplained all 0.
- run_manifest.json contains warnings with "TODO_PHYSICS" and config_hash.
- second run into same dir -> InputContractError; CLI simulate-nominal exit 0 then 2.
