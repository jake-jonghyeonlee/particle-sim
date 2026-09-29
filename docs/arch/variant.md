# variant (local mask edits for avoidance candidates)

`particle-sim make-variant --product-config BASE.yaml --edits SPEC --output DIR [--candidate ID ...]`

Global CD bias / overlay (scenario sampling) is not an avoidance edit. Variants remove electrode
pixels only inside an ROI; removed pixels are also cleared from the electrode's layer mask unless
another electrode of the same layer occupies them.

Edit spec (`variant/package.py`, pydantic, extra keys forbidden), JSON or YAML:

```yaml
schema_version: "1.0"
length_unit: um                 # ROI origin upper-left, x right, y down
candidates:
  - candidate_id: P1-v1-S_M2-A__B-e1   # [A-Za-z0-9_.-]+, output folder name
    sample_id: P1-v1-S_M2-A__B         # optional; needed for slot/target_slot
    design_version: null               # default: <base>-<last '-' token of candidate_id>
    change_description: "..."
    edits:
      - op: trim_near | remove_region | remove_overlap
        electrode: B        # or slot: a|b (resolved from sample_id with export-simcompare rules)
        target: A           # or target_slot; required except remove_region
        roi_um: [x0, y0, x1, y1]
        amount_um: 1.0      # trim_near > 0; remove_overlap distance allowance (default 0)
    basis: {...}            # free-form, recorded only
    info: {...}             # free-form, recorded only
```

Ops (`variant/edits.py`, scipy EDT on pixel centres in um):
- `remove_region`: remove every electrode pixel in the ROI.
- `trim_near`: remove ROI pixels closer to the target than (ROI minimum distance + amount).
- `remove_overlap`: remove ROI pixels within `amount_um` of the target (0 = planar overlap).

Per variant `DIR/<candidate_id>/`: manifest.yaml (only design_version changed), every referenced
mask (edited ones rewritten in the same format, others copied), `variant_record.json`
(parent_design_id = `<product>-<base version>`, edits with removed area and ROI distance /
overlap before-after, checks: electrode components before/after and `split`, contact masks
touched), `edit_preview.png` (grey original, red removed, blue target, yellow ROI).

Errors (no package written): empty ROI, missing target, unknown electrode/sample_id, emptied
electrode, broken open-path terminals (`check_open_baselines`), same design_version, non-empty
output. Other checks are reported only; tolerance and design rules are the designer's decision.

`export-simcompare` reads `variant_record.json` beside the manifest and fills `parent_design_id`
and `change_description`; a record whose product_id/design_version differ is an error.
