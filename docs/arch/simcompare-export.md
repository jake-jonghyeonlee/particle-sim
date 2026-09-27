# arch: simcompare-export

The nominal exporter writes mapping contract 1.2, one bundle per applicable supported
(source step, failure unit). Supported units are lateral/vertical electrode pairs and
explicit/automatic lateral opens. CNO is excluded with warnings; excluded baseline opens
remain exported with their classification warning. Simulation and CA values are unchanged.

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

Layouts share one uint8 NPY per electrode. Each hotspot plane is inspected during preflight
and stored as uint8 if binary, otherwise float32. Capacity uses each actual storage dtype;
only NPY headers, CSV quoting and filesystem allocation have conservative slack. Planes
are read again during writing to bound memory use. No arrays are retained for the full run.

All validation and capacity checks precede writing. Outputs stage beneath the target parent
and rename atomically, with cleanup on failure and a second overwrite guard. Final directories
use 0777 masked by the process umask; files use 0666 masked by umask. The staging root gets
its mode from the normal arrays directory, without temporarily modifying global umask.
No Fab measurements are generated. Nonempty output paths are never overwritten.

CLI: `export-simcompare --run DIR --product-config PATH --output DIR`.
Tests cover source assignment/key/provenance, legacy layer stack, raw identity and collisions,
applicability, dtype preservation, capacity, permissions, atomic cleanup, and CLI results.
