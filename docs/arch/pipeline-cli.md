# arch: pipeline and CLI

`validate-input` loads the manifest/masks, reports residual layer regions and emits the complete electrode
pair audit. Schema 2.0 unresolved pairs fail before geometry work.

`simulate-nominal` evaluates every source step, particle size and failure unit, writes binary/margin maps,
mode unions, raw/allocated/unique CA and a reproducibility manifest. `simulate-probability` does the same
per Sobol scenario, then writes probabilities, CA distributions and seed-aware convergence.

`compare-sim` consumes an explicit reference YAML/JSON (schema 2.0) that names the product, design version
and the SHA-256 of the mask files the old sim used; all must match the candidate run manifest. Its records map source step, diameter and failure
unit to a 2-D `.npy` or TIFF binary map. It rejects axis, origin, pixel-size (including a run manifest without
pixel sizes), shape or incomplete-common-axis mismatches before creating any output. Outputs are comparison tables, intersection/reference-only/candidate-only/XOR Zarr maps and a
manifest. Candidate-only electrode pairs absent from the old sim are `new_pair`; open and CNO units are `new_unit`;
conditions whose source step does not apply are `not_applicable`. None of them is a zero-reference difference.

`apply-fab --run` creates the additive risk views. `--ca` remains a legacy compatibility mode.
`make-report` renders per-run diagnostic figures. `export-simcompare` remains the downstream bundle
compatibility path and is updated only after the local old-sim replacement decision.
