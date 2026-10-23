# arch: synthetic-validate (T3)

Files: `src/particle_sim/validation/synthetic_cases.py`, `src/particle_sim/cli.py`,
`src/particle_sim/__main__.py`, tests `tests/unit/test_synthetic_validate.py`.
Depends on schema-config (T1) and io-masks-stack (T2).

## validation/synthetic_cases.py
Synthetic data only (not physical results). All cases: shape (200, 200), pixel 0.5 um
(roi 100 x 100 um), product_family "SYN", design_version "V01", design_lineage_id "SYN_L1",
fab_id "FAB_SYN", is_synthetic True, embedding_depth_um 0.0, contact_margin_um 0.0,
particle_sizes_um [2.0, 4.0, 6.0, 8.0, 10.0], legacy_kernel default,
post_step_policy {mode: "protrusion", coverage_break_factor: 1.5}.

Stack (process_order, name, material, thickness_um, planarization_factor, patterned):
`0 M1 conductor 0.3 None True`, `1 GI inorganic 0.4 0.0 False`,
`2 M2 conductor 0.3 None True`, `3 PLN organic 2.0 0.8 False`, `4 M3 conductor 0.1 None True`.
Steps (step_id lands_after_layer surface_model; fab_distribution_layer omitted, derived as the
mask layer (patterned or contact hole layer) directly above: S_M1 -> M2, S_M2 -> M3, S_PLN -> M3; `write_fab_table` writes M2 and M3):
`S_M1 M1 conformal`, `S_M2 M2 conformal`, `S_PLN PLN flat`.

Electrodes are full-height vertical stripes: `mask[:, x0:x1] = True` (x1 exclusive).
Layer mask = union of that layer's electrodes (all-False layer mask if none).

| case | electrodes (id layer net role x0:x1) | pairs (pair_id mode) |
|---|---|---|
| parallel_lines | A M2 N1 DATA 80:90; B M2 N2 DATA 96:106 | A__B lateral_short |
| vertical_aligned | G M1 N1 GATE 90:110; P M3 N2 PIXEL 90:110 | G__P vertical_short |
| diagonal_offset | G M1 N1 GATE 90:110; P M3 N2 PIXEL 116:136 | G__P vertical_short |
| triple_overlap | A M2 N1 DATA 80:90; B M2 N2 DATA 92:93; C M2 N3 DATA 95:105 | A__B, B__C, A__C lateral_short |
| pilot | G M1 N1 GATE 40:60; A M2 N2 DATA 66:76; B M2 N3 DATA 82:92; P M3 N4 PIXEL 98:118 | A__B lateral_short; G__A vertical_short; B__P vertical_short |

pair_id is `f"{a}__{b}"`; electrode_a is the first id.

```python
CASE_NAMES: tuple[str, ...] = ("parallel_lines", "vertical_aligned", "diagonal_offset",
                               "triple_overlap", "pilot")
@dataclass(frozen=True)
class SyntheticCase:
    name: str
    manifest: dict                        # plain YAML-ready dict, paths relative
    layer_masks: dict[str, np.ndarray]    # bool
    electrode_masks: dict[str, np.ndarray]  # bool
def build_case(name: str) -> SyntheticCase   # unknown name -> ValueError
def write_case(case: SyntheticCase, out_dir: Path) -> Path
```
`write_case` writes `masks/layer_<name>.png` and `masks/electrode_<id>.png` as uint8 0/255
(mode "L"), sets manifest paths to those relative paths, product_id = `"SYN_" + name.upper()`,
writes `manifest.yaml` (yaml.safe_dump, sort_keys=False) and returns its path.

## cli.py
```python
def build_parser() -> argparse.ArgumentParser   # prog "particle-sim", Korean help
def main(argv: list[str] | None = None) -> int
```
Subparsers stored with `dest="command"`, each sets `func=<handler(args) -> int>` so later
tasks add commands without changing existing ones.
- `validate-input --product-config PATH`: load_manifest, load_product_masks, build_stack_z;
  print JSON (indent=2, ensure_ascii=False) to stdout:
  `{"status": "ok", "product_id", "shape_yx": [y, x], "pixel_size_um": [py, px],
  "layers": [...], "electrodes": n, "pairs": n, "source_steps": [...],
  "warnings": todo_physics_warnings(...)}`; return 0.
- `make-synthetic --case NAME --output DIR`: write_case, print manifest path, return 0.
- InputContractError anywhere -> print `입력 오류: {err}` to stderr, return 2.
`__main__.py`: `raise SystemExit(main())`.

Tests: every case builds and `write_case` output passes `validate-input` (exit 0, JSON has
warnings with "TODO_PHYSICS"); corrupt one electrode PNG to gray 128 -> exit 2 and stderr
contains "field="; pair stripes match the table exactly.
