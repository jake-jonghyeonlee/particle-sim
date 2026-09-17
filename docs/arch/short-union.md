# arch: short-union (T6)

Files: `src/particle_sim/failure/lateral_short.py`, `src/particle_sim/failure/vertical_short.py`,
`src/particle_sim/failure/union.py`, tests `tests/unit/test_short_union.py`.
Depends on T5: `ContactResult(contact: bool ndarray, margin_um: float32 ndarray, kind: str)`
and `compute_electrode_contacts(manifest, masks, stack_z, step_id, diameter_um, cache,
electrode_ids=None) -> dict[str, ContactResult]`, `DistanceCache(masks)`.
Manifest types from T1 (`ProductManifest`, `ElectrodePair`, `Electrode`).

## Shared result type (define in failure/lateral_short.py, import elsewhere)
```python
@dataclass(frozen=True)
class PairResult:
    pair_id: str
    failure_mode: str            # "lateral_short" | "vertical_short"
    electrode_a: str
    electrode_b: str
    failure: np.ndarray          # bool (y, x) = contact_a & contact_b
    contact_a: np.ndarray        # bool
    contact_b: np.ndarray        # bool
    pair_margin_um: np.ndarray   # float32 = minimum(margin_a, margin_b)
    classification: str
```

## failure/lateral_short.py
```python
def lateral_short(manifest: ProductManifest, pair: ElectrodePair,
                  contacts: dict[str, ContactResult]) -> PairResult
```
Re-check in this order (InputContractError, field=f"electrode_pairs[{pair.pair_id}]"): mode is
lateral_short, both electrodes exist, same layer_name, different net_id, contacts present.
classification: `"direct_geometric_contact"`.

## failure/vertical_short.py
```python
def vertical_short(manifest: ProductManifest, pair: ElectrodePair,
                   contacts: dict[str, ContactResult]) -> PairResult
```
Same checks but mode vertical_short and different layer_name.
classification: `"potential_vertical_short"` always in this phase (particle conductivity,
insulator continuity, step coverage are TODO_PHYSICS). If either ContactResult.kind is
"post_step" append nothing else — kind is kept in the per-electrode contact only.

```python
def simulate_pairs(manifest: ProductManifest, masks: MaskSet, stack_z: dict[str, LayerZ],
                   step_id: str, diameter_um: float, cache: DistanceCache) -> list[PairResult]
```
(define in failure/union.py) computes contacts only for electrodes used by pairs, then
dispatches each pair in manifest order to lateral_short / vertical_short.

## failure/union.py
```python
ANY_MODE = "any"
@dataclass(frozen=True)
class UnionResult:
    failure_mode: str         # a mode name or "any"
    pair_ids: tuple[str, ...]
    any_failure: np.ndarray   # bool (y, x) OR of member pair failures
    pair_count: np.ndarray    # uint16 (y, x) number of member pairs failing at the pixel

def union_by_mode(results: list[PairResult]) -> dict[str, UnionResult]
```
Keys: each failure_mode present (manifest order of first appearance) plus "any" over all
pairs. Empty list -> ValueError. Shapes must match -> ValueError. Never sums CA numbers.

Tests (use `particle_sim.validation.synthetic_cases.build_case` + `write_case` to tmp_path,
then `load_manifest`, `load_product_masks`, `build_stack_z`; synthetic contract:
cases parallel_lines, vertical_aligned, diagonal_offset, triple_overlap, pilot with steps
S_M1 / S_M2 / S_PLN and sizes 2..10 um, pixel 0.5 um, 200x200):
- parallel_lines, step S_M2, d=10: failure non-empty and lies between the stripes;
  d=2: failure empty. Hotspots are monotone in size: H(d1) subset H(d2) for d1 < d2.
- vertical_aligned vs diagonal_offset, step S_M1, d=10: both non-empty and the mean x of the
  diagonal hotspot is larger than the aligned one (moves toward the shifted upper stripe).
- post_step_policy exclude (model_copy with update) -> vertical failure empty from S_M1.
- triple_overlap S_M2: d=10 -> union "lateral_short" area < sum of pair areas and
  pair_count.max() == 3; d=6 -> A__C empty and union area == sum of pair areas.
- same net or wrong layer relation -> InputContractError.
- pair_margin_um == minimum of the two margins; classification strings.
