# arch: Fab input and risk weighting

`io/fab.py` supports the legacy total-count table and the preferred mean-per-sample contract. The preferred
columns are `fab_id`, `product_family`, `fab_distribution_layer`, size-bin bounds and representative size,
`mean_particle_count`, `inspection_area_per_sample_mm2`, `measurement_period`, `particle_type`, with
optional `mean_density_per_mm2`, `sample_count`, `source_id`, and `source_hash`.

Density is mean count divided by per-sample inspection area; a supplied density must agree. Mean inputs use
`exact_point` only so one measured bin is not repeated over several simulated diameters.

`apply-fab --run` reads `ca.parquet` and `ca_union.parquet`. It emits:

- `risk_by_layer.parquet`: source-layer OR-union and allocated affected-layer views;
- `risk_by_pair.parquet`: raw, `1/k` allocated and unique pair risk;
- `risk_by_open_path.parquet`: the same metrics for open;
- `risk_by_electrode.parquet`: involvement and additive electrode contribution;
- `risk_summary.csv` and `risk_manifest.json`.

All rows retain Fab input semantics, count, inspection area, period/type and available source IDs/hashes.
The calculated value is a geometric hit expectation, not a confirmed electrical defect count.
Before weighting, the CLI checks schema/run/product IDs and verifies `ca.parquet` and `ca_union.parquet`
against their run-manifest hashes. The risk manifest stores hashes for the product config, Fab table and
source run manifest as well as every generated risk artifact.
