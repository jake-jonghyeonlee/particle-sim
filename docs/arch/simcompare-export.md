# arch: simcompare-export (T9)

Files: `src/particle_sim/export/simcompare.py`, `src/particle_sim/cli.py` (add subcommand),
tests `tests/unit/test_simcompare_export.py`. Depends on T8 run directory layout:
`hotspots.zarr` (zarr group; array `nominal_failure` uint8 [S, D, P, Y, X]; attrs
`source_steps`, `particle_sizes_um`, `pair_ids`), `ca.parquet` (columns include
source_step, particle_size_um, pair_id, ca_nominal_normalized), `run_manifest.json`
(run_id, config_hash, simulator_version, assumption_set, legacy_kernel, post_step_policy).
Uses `load_manifest`, `load_product_masks`, `build_stack_z`, `open_hotspot_store`.

Target: the simcompare example adapter data_dir (read-only reference contract):
```text
bundles.csv      sample_id,product_id,product_family,design_version,design_lineage_id,
                 electrode_a_role,electrode_b_role,layer_a,layer_b,failure_type,fab_id,
                 fab_distribution_id,assigned_process_layer,unit_width,unit_height,
                 physical_length_unit,parent_design_id,change_description,calc_version
layouts.csv      sample_id,slot,path,pixel_size,origin_x,origin_y
hotspots.csv     sample_id,particle_size,path,particle_size_unit,pixel_size,origin_x,origin_y
ca_records.csv   sample_id,particle_size,particle_size_unit,ca_value,ca_unit
layer_stack.csv  product_id,design_version,layer_name,stack_order,fab_process_layer
geometry.json    {sample_id: {"thickness": {layer: um}, "shape_metrics": {}, "length_unit": "um"}}
provenance.json  {"calc_settings": {...}, "calc_version": str, "is_synthetic": bool}
arrays/*.npy     2-D float32, allow_pickle=False
```

```python
def export_simcompare(run_dir: Path, manifest_path: Path, out_dir: Path) -> Path
```
- out_dir exists and non-empty -> InputContractError(field="output"). Returns out_dir.
- Grid must be isotropic (simcompare has one pixel_size) else InputContractError.
- One bundle per (source_step, pair) in store order:
  `sample_id = f"{product_id}-{design_version}-{step_id}-{pair_id}"`; must match
  `^[A-Za-z0-9_.-]+$` else InputContractError(field="sample_id").
- bundles: roles/layers from electrodes a/b; failure_type = pair failure_mode; fab_id from
  manifest; fab_distribution_id, parent_design_id, change_description empty;
  assigned_process_layer = step.fab_distribution_layer for lateral_short, EMPTY for
  vertical_short; unit_width/unit_height = roi_width_um/roi_height_um;
  physical_length_unit "um"; calc_version `f"particle-sim-{simulator_version}"`.
- layouts: slot a/b rows -> `arrays/{sample_id}-a.npy` / `-b.npy` electrode masks as float32
  0/1; pixel_size = pixel size um; origin 0.0, 0.0.
- hotspots: one row per size -> `arrays/{sample_id}-hot{size:g}.npy` float32 from
  nominal_failure; particle_size_unit "um".
- ca_records: ca_value = ca_nominal_normalized (dimensionless), ca_unit "1".
- layer_stack: every stack layer, stack_order = process_order, fab_process_layer = layer_name,
  design_version = manifest design_version.
- geometry.json thickness = all stack layers; provenance calc_settings = {"run_id",
  "config_hash", "assumption_set", "legacy_kernel", "post_step_policy",
  "hotspot_semantics": "nominal binary geometric contact risk"}, is_synthetic from manifest.
- `export_warnings.json` list of Korean strings: always a note that Fab input
  (fab.json / fab_measurements.csv) is NOT generated; per vertical_short bundle a note that
  simcompare forces the upper electrode layer as Fab layer while the simulator uses the
  particle source step (`step_id`), so assigned_process_layer was left empty.
CSV written with the csv module, header exactly as listed, floats via `repr`-free `str()`.

## cli.py addition
`export-simcompare --run DIR --product-config PATH --output DIR`; prints JSON
`{"status": "ok", "output": str, "bundles": n, "warnings": n}`; return 0.

Tests (synthetic pilot -> `run_nominal` -> export in tmp_path):
- headers exactly as above; 9 bundles (3 steps x 3 pairs); 45 hotspot rows; 18 layout rows
  (9 bundles x 2 slots); npy arrays 2-D float32 shape (200, 200).
- ca_value equals ca.parquet ca_nominal_normalized for a sampled row.
- vertical bundles have empty assigned_process_layer and a matching warning.
- non-empty output dir -> error; CLI exit codes 0 / 2.
