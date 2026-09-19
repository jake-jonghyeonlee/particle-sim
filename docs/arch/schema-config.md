# arch: schema and configuration

`schema.py` defines frozen Pydantic models with `extra="forbid"`; `config.py` converts YAML and
validation failures to `InputContractError`.

Product schema `1.0` remains readable for old fixtures. New product inputs use `2.0`:

- `layer_masks` describes the full patterned layer; `electrodes` describes evaluated conductor regions.
- `electrode_pairs` contains calculated different-net pairs.
- `excluded_electrode_pairs` contains every intentionally skipped different-net pair and a nonempty reason.
- loading schema 2.0 audits every possible different-net combination and rejects unresolved pairs.
- `contacts` lists contact·via masks with a `hole_layer`; lower/upper conductor layers default to the
  nearest patterned conductors in `process_order`. Each contact needs `applicable_particle_source_steps`.
  Every hole must join exactly one lower and one upper electrode of the same net (`geometry/contacts.py`).
- `open_applicable_steps_by_layer` enables automatic-terminal open per layer; every layer joined by a contact
  must be declared (an empty list disables automatic open there). Unit IDs may not start with the reserved
  `auto_open:` or `cno:` prefixes.
- `open_paths` accepts either two legacy `(y,x)` points or two terminal mask paths, never both. Schema 2.0
  requires `applicable_particle_source_steps`.

Scenario schema `2.0` uses `OverlayXYSpec(x: NormalSpec, y: NormalSpec)`,
`overlay_reference_layer`, `overlay_sign_convention`, and `independence_model: independent`.
A measured overlay or CD `NormalSpec` requires `source_id`; optional provenance fields are
`source_hash`, `measurement_period`, `fab_process`, and `product_family`.
Probability execution cross-checks overlay/CD names against `layer_masks` and contact hole layers, thickness/planarization names
against `stack`, and rejects a moving or unknown overlay reference layer.

All identifiers are unique; cross references, material/layer rules, lateral/vertical layer rules, source
steps, open electrodes and applicable steps are checked before any mask calculation. Paths are relative to
the config directory and cannot escape it. `canonical_hash` and `file_sha256` provide reproducibility IDs.
