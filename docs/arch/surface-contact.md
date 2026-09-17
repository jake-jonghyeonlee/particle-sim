# arch: surface-contact (T5)

Files: `src/particle_sim/geometry/surfaces.py`, `src/particle_sim/geometry/particles.py`,
`src/particle_sim/geometry/contact.py`, tests `tests/unit/test_surface_contact.py`.
Depends on T2 (`MaskSet`, `GridSpec`, `LayerZ`, `build_stack_z`, `layers_up_to`,
`layers_after`) and T4 (`distance_to_mask_um`, `distance_to_z_interval_um`).
Electrode z intervals are NOMINAL (flat cumulative stack) everywhere; only the particle
landing surface has topography. This approximation is intentional for this phase.

## geometry/surfaces.py
```python
def surface_height_um(manifest: ProductManifest, step_id: str,
                      stack_z: dict[str, LayerZ], masks: MaskSet) -> np.ndarray
```
float64 (y, x). `z_land = stack_z[step.lands_after_layer].z_top_um`; `below = layers_up_to`.
- flat: `np.full(shape, z_land)`
- conformal: `z_land - sum(k.thickness_um * (~masks.layers[k]) for patterned k in below)`
- partial: same sum but each term multiplied by `retention_k = prod(1 - j.planarization_factor
  for j in below if j.process_order > k.process_order and j.material_type != "conductor")`
  (conductor films keep topography, factor 1). Missing factor -> InputContractError.
Unknown step -> InputContractError(field="step_id").

## geometry/particles.py
```python
def particle_center_z_um(surface_z_um: np.ndarray, diameter_um: float,
                         embedding_depth_um: float) -> np.ndarray
# surface + d/2 - e. d <= 0 or e < 0 or e > d -> InputContractError(field="embedding_depth_um")
def post_step_eligibility(manifest: ProductManifest, step_id: str,
                          stack_z: dict[str, LayerZ], diameter_um: float) -> dict[str, bool]
```
Keys: every conductor layer formed after the landing layer (`layers_after`).
- exclude: all False
- wrap_contact: all True
- protrusion: `h = d - e`, `covered = False`, `k = coverage_break_factor`; iterate layers_after
  in process order:
  - conductor: `result[name] = (h > 0) and not covered`
  - other: `broken = h > k * thickness`; if not broken: `covered = True`;
    `h = max(0.0, h - thickness * planarization_factor)`

## geometry/contact.py
```python
@dataclass(frozen=True)
class ContactResult:
    contact: np.ndarray      # bool (y, x)
    margin_um: np.ndarray    # float32 (y, x); -inf where contact impossible
    kind: str                # "3d" | "post_step" | "excluded"

def contact_3d(distance_xy_um: np.ndarray, center_z_um: np.ndarray, z_bottom_um: float,
               z_top_um: float, radius_um: float, contact_margin_um: float) -> ContactResult
# margin = r + cm - hypot(dxy, dz); contact = margin >= 0 (Q == R counts as contact)
def contact_post_step(distance_xy_um: np.ndarray, eligible: bool, radius_um: float,
                      contact_margin_um: float) -> ContactResult
# eligible: margin = r + cm - dxy, kind "post_step"; else margin -inf, contact False, "excluded"

class DistanceCache:
    def __init__(self, masks: MaskSet) -> None
    def get(self, electrode_id: str) -> np.ndarray   # distance_to_mask_um, computed once
    @property
    def computed(self) -> int                          # number of distance maps computed

def compute_electrode_contacts(manifest: ProductManifest, masks: MaskSet,
                               stack_z: dict[str, LayerZ], step_id: str,
                               diameter_um: float, cache: DistanceCache,
                               electrode_ids: Iterable[str] | None = None
                               ) -> dict[str, ContactResult]
```
surface -> center z -> for each electrode (all or `electrode_ids`, unknown id -> error):
if its layer process_order <= landing order use `contact_3d` with that layer's z interval,
else `contact_post_step` with `post_step_eligibility[layer]`. r = d / 2.

Tests (build manifests/masks in memory with small stripes, pixel 0.5 um):
- flat / conformal / partial surface values at chosen pixels (hand computed).
- center z formula and embedding error cases.
- protrusion: with GI 0.4 (p 0), PLN 2.0 (p 0.8), k 1.5 from landing M1: d=0.5 -> M2 False
  (0.5 <= 0.6 covered); d=4 -> M2 True, M3 True (4 > 3.0; h=4-1.6>0); d=2 -> M3 False.
- contact_3d: Dz > r -> no contact anywhere; margin sign; Q == R boundary is contact;
  diagonal hypot(3, 4) with r=5 -> margin 0 -> contact.
- flat surface with embedding_depth = d/2 on the landing layer electrodes gives contact maps
  identical to `distance_to_mask_um <= r` (legacy EDT equivalence when Dz = 0).
- DistanceCache computes each electrode once across repeated calls / sizes.
