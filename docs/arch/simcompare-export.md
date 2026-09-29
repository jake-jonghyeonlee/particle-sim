# arch: simcompare-export

The nominal exporter writes mapping contract 1.2, one bundle per applicable supported
(source step, failure unit). Supported units are lateral/vertical electrode pairs and
explicit/automatic lateral opens. CNO is excluded with warnings; excluded baseline opens
remain exported with their classification warning. Simulation and CA values are unchanged.

CA exports retain raw area in `um^2`: `ca_nominal_um2` by default, or `ca_mean_um2`
via the Python API's `ca_mode="mean"`. Mean mode requires a probability run, identified
by the run manifest's `value_semantics = PROBABILITY_SEMANTICS` as written by
`run_probability`; nominal runs are rejected before writing. Provenance `calc_settings`
records `ca_mode` and `ca_source` (the selected raw CA column). Hotspot arrays always
remain the nominal binary footprint, including when CA is the probability mean.
Limitation: `simulate-probability` output directories contain `probability.zarr` but no
`hotspots.zarr`, so mean-mode export currently stops with a controlled "cannot read run"
error on pipeline probability runs; inspect probability-run mean CA with
particle-sim-compare's run checker (`--ca mean`) instead.

`bundles.csv` carries `assigned_process_layer = step.fab_distribution_layer` for every type,
including vertical shorts. Preflight rejects blank source labels and, when CA carries
`fab_distribution_layer`, any label that differs from the manifest or is not uniform within
its source step. Physical upper electrode layers do not determine Fab assignment.

`fab_pair_key = "<raw unit>:<layer_a>-<layer_b>"` preserves endpoint order, punctuation and
Unicode; opens use literal `N/A` for b and only layout slot a. Components must be nonblank
and have no edged whitespace. A registry rejects different raw tuples that compose the
same key; identical tuples may repeat across steps. Sample IDs retain their existing
sanitization and collision checks. Pair IDs are product-local; layer-wide Fab is the default
and pair-specific real measurements require explicit IDs and established identity scope.

`provenance.json` carries `mapping_contract_version: "1.2"`, existing calculation/storage
metadata, and each `sample_id_mapping` entry's raw tuple, sample ID, pair key and assignment.
`layer_stack.csv` retains its legacy columns and physical layer aliases, which do not
provide 1.2 source assignment. Consumers must support 1.2; keep old exports as legacy artifacts.

Layouts share one file per electrode. Each layout and hotspot plane is inspected during
preflight: binary arrays use lossless 8-bit grayscale PNG (0/255, zlib compression level 9),
and nonbinary arrays retain float32 NPY. The adapter's L-mode float32 division by 255
recovers exact binary 0/1 values; layout thresholding remains unchanged. Provenance records
each array path's format and both encodings. Casefold collision checks use selected extensions.
Capacity conservatively bounds PNG by uncompressed uint8 pixels plus row filter bytes,
worst-case DEFLATE expansion and header/allocation slack; it never assumes compression
savings. NPY uses float32 payload size plus slack; CSV quoting and metadata are included.
Planes are read again during writing to bound memory use; no full-run arrays are retained.

All validation and capacity checks precede writing. Outputs stage beneath the target parent
and rename atomically, with cleanup on failure and a second overwrite guard. Final directories
use 0777 masked by the process umask; files use 0666 masked by umask. The staging root gets
its mode from the normal arrays directory, without temporarily modifying global umask.
No Fab measurements are generated. Nonempty output paths are never overwritten.

CLI: `export-simcompare --run DIR --product-config PATH --output DIR`.
Tests cover source assignment/key/provenance, legacy layer stack, raw identity and collisions,
applicability, dtype preservation, capacity, permissions, atomic cleanup, and CLI results.
