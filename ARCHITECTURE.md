# ARCHITECTURE: particle-sim (phase 0-3)

Module contracts are in `docs/arch/<task>.md`; SPEC.md holds the global rules.

| task | modules | depends |
|---|---|---|
| schema-config | errors, schema, config | - |
| io-masks-stack | io/masks, io/stack | schema-config |
| synthetic-validate | validation/synthetic_cases, cli, __main__ | io-masks-stack |
| legacy-distance | geometry/distance, geometry/legacy, validation/regression | schema-config |
| surface-contact | geometry/surfaces, geometry/particles, geometry/contact | io-masks-stack, legacy-distance |
| short-union | failure/lateral_short, failure/vertical_short, failure/union | surface-contact, synthetic-validate |
| ca-outputs | aggregation/critical_area, io/outputs | io-masks-stack |
| pipeline-cli | pipeline, cli | short-union, ca-outputs |
| simcompare-export | export/simcompare, cli | pipeline-cli |
| figures-report | visualization/hotspots, visualization/reports, cli | simcompare-export |

Data flow:
```text
manifest.yaml -> load_manifest -> load_product_masks (MaskSet) + build_stack_z (LayerZ)
  -> for source_step, particle size:
       DistanceCache (xy EDT per electrode, reused across steps/sizes)
       surface_height -> particle center z / post-step eligibility -> ContactResult per electrode
       -> PairResult (lateral / vertical) -> union_by_mode (OR, pair_count)
       legacy_short_hotspot (5x5 erode/dilate baseline) + explain_legacy_vs_edt
  -> HotspotStore (zarr) + ca.parquet + ca_union.parquet + legacy_vs_edt.csv + run_manifest.json
  -> export_simcompare (simcompare data_dir) / make_report (PNG + report.md)
```

Known approximations (reported as TODO_PHYSICS, not hidden):
- electrode z intervals are nominal flat cumulative stack values; only the landing surface
  has topography (flat / conformal / partial)
- particle is a conductive sphere; embedding depth and contact margin are manifest inputs
- layers formed after the source step: exclude / protrusion(planarization, coverage break k)
  / wrap_contact
- vertical results are potential_vertical_short until process calibration

## Phase 2 modules

| module | role |
|---|---|
| `failure/open_failure` | terminal connectivity, skeleton+terminal candidate reduction, proxy fallback |
| `units` | one view over pairs and open paths (`lateral_open`, electrode_b = "N/A") |
| `uncertainty/transforms`, `uncertainty/sampler` | inverse-CDF transforms, deterministic variable order, Sobol scenarios |
| `geometry/transforms` | scenario overlay / CD bias (signed distance field), thickness delta, scalar overrides |
| `aggregation/probability` | per-pixel failure counts, probability map, CA percentiles |
| `validation/convergence` | CA against scenario count, caller-provided tolerance |
| `io/fab`, `aggregation/fab_weighting` | Fab table validation and common / actual weighting |
| `pipeline.run_probability` | scenario loop, probability store, CA and convergence tables |
| `visualization` | probability map, convergence and Fab comparison figures |

```text
manifest + scenario.yaml -> build_scenarios (Sobol, deterministic order)
  -> per scenario: transform masks/stack/scalars -> evaluate_units -> accumulate counts
  -> probability map + CA percentiles + convergence prefixes (64/128/256)
  -> optional: apply_fab_weighting(CA, fab table) -> common / actual expected defect count
```

Still deferred: tile + halo for very large ROIs, empirical Fab size samples, electrical
calibration of vertical shorts.
