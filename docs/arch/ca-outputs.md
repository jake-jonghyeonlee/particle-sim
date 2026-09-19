# arch: CA and outputs

Output schema is `2.0`. `failure_unit_id` is the canonical axis ID and `unit_kind` distinguishes
`electrode_pair`, `open_path`, `open_electrode` (automatic terminals) and `contact` (CNO); compatibility `pair_id` fields remain in core CA/Zarr outputs.

`ca.parquet` preserves legacy, nominal, mean, p05/p50/p95, allocated and unique CA. It also records actual
electrode IDs/layers, failure mode, scenario count and convergence delta. `convergence_delta` is the worst
(largest) final prefix comparison over the scramble seeds, the same step `is_converged` checks; an infinite
change is kept and NaN means no seed had two prefixes. `ca_union.parquet` keeps nominal
union separate from variation mean and percentiles, plus the raw unit sum used to expose overlap.

Nominal `hotspots.zarr` stores legacy/nominal failure, margin, union failure and simultaneous unit count.
Probability Zarr stores unit probability, mean margin and union probability. Attributes name source steps,
particle sizes, failure-unit IDs, unit kinds, union keys, axes, units and value semantics.

CNO rows use `failure_mode=contact_not_open`, fill `contact_id`, and name the lower and upper electrode in
`electrode_a`/`electrode_b`. Unions add a `contact_not_open` key and `any` includes CNO.

`run_manifest.json` records schema/input versions, config and input hashes, git commit/dirty state, build ID,
shape and pixel sizes, physical assumptions, scenario provenance, storage strategy and SHA-256/size for each
output artifact. Runs with contacts add `cno_coverage_rule`; probability runs add `scenario_exclusions`
and `scenario_contact_lost_holes`.
