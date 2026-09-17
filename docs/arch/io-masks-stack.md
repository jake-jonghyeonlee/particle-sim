# arch: io-masks-stack (T2)

Files: `src/particle_sim/io/masks.py`, `src/particle_sim/io/stack.py`,
tests `tests/unit/test_io_masks_stack.py`. Depends on schema-config (T1).

## io/masks.py
```python
@dataclass(frozen=True)
class GridSpec:
    shape_yx: tuple[int, int]
    pixel_size_y_um: float
    pixel_size_x_um: float
    origin: str = "upper_left"
    @property
    def pixel_area_um2(self) -> float
    @property
    def roi_area_um2(self) -> float          # shape_y*pixel_y * shape_x*pixel_x
    @property
    def is_isotropic(self) -> bool           # abs(py - px) <= 1e-9

@dataclass(frozen=True)
class MaskSet:
    grid: GridSpec
    layers: dict[str, np.ndarray]      # layer_name -> bool (y, x)
    electrodes: dict[str, np.ndarray]  # electrode_id -> bool (y, x)
    input_hashes: dict[str, str]       # manifest-relative path -> sha256

def load_binary_mask(path: Path, *, product_id: str, field: str) -> np.ndarray
```
- `.png` via Pillow, `.tif/.tiff` via tifffile; other suffix -> InputContractError.
- Accept 2-D only (mode "L", "1", or single-channel TIFF). RGB/RGBA -> error.
- Allowed unique values: subset of {0, 1} or subset of {0, 255} (bool also ok).
  Any other value -> error "이진 마스크가 아님(resize/antialias 의심)" with count of bad values.
- Returns `arr != 0` as bool. Never resize or threshold.

```python
def load_product_masks(manifest: ProductManifest, manifest_path: Path) -> MaskSet
```
Uses `config.resolve_input_path` + `file_sha256`. Checks, each raising InputContractError
with product_id, file and field (`layer_masks[i].path` / `electrodes[i].path`):
1. all masks share one shape (report first mismatching file and both shapes)
2. `abs(shape_x*pixel_size_x_um - roi_width_um) <= pixel_size_x_um/2` and same for y
3. no electrode is empty (all False)
4. electrodes of the same layer do not overlap (report overlapping pixel count)
5. for each layer that has at least one electrode: union of its electrodes == that layer
   mask exactly (report missing/extra counts). Layers without electrodes are topography only.
GridSpec built from shape + coordinate_system.

## io/stack.py
```python
@dataclass(frozen=True)
class LayerZ:
    layer_name: str
    process_order: int
    material_type: str
    thickness_um: float
    z_bottom_um: float
    z_top_um: float
    planarization_factor: float | None
    patterned: bool

def build_stack_z(manifest: ProductManifest) -> dict[str, LayerZ]
```
Sort stack by process_order; z_bottom of first = 0.0; each next z_bottom = previous z_top;
z_top = z_bottom + thickness. If a StackLayer has z_bottom_um, require
`abs(given - computed) <= 1e-6` else InputContractError(field=f"stack[{name}].z_bottom_um").
Return dict in process order (insertion ordered).

```python
def landing_surface_z_um(manifest: ProductManifest, step_id: str,
                         stack_z: dict[str, LayerZ]) -> float
# z_top of step.lands_after_layer; unknown step -> InputContractError(field="step_id")
def layers_up_to(step: ParticleSourceStep, stack_z: dict[str, LayerZ]) -> list[LayerZ]
# layers with process_order <= landing layer's order, process order
def layers_after(step: ParticleSourceStep, stack_z: dict[str, LayerZ]) -> list[LayerZ]
# layers with process_order > landing layer's order, process order
```
Tests: write PNG/TIFF masks in tmp_path with a manifest dict built in the test (same fields
as docs/arch/schema-config.md), cover every check above, 0/255 and 0/1 accepted, gray value
128 rejected, RGB rejected, cumulative z values, z_bottom cross-check failure.
