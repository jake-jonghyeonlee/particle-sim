# arch: legacy-distance (T4)

Files: `src/particle_sim/geometry/distance.py`, `src/particle_sim/geometry/legacy.py`,
`src/particle_sim/validation/regression.py`, tests `tests/unit/test_legacy_distance.py`.
Depends on schema-config (T1). Uses `GridSpec` from `particle_sim.io.masks` (T2 contract:
frozen dataclass `GridSpec(shape_yx, pixel_size_y_um, pixel_size_x_um, origin)` with
properties `pixel_area_um2`, `roi_area_um2`, `is_isotropic`).

## geometry/distance.py
```python
def distance_to_mask_um(mask: np.ndarray, grid: GridSpec) -> np.ndarray
# float64 (y, x). scipy.ndimage.distance_transform_edt(~mask, sampling=(py, px)).
# Inside mask 0. Empty mask -> np.full(shape, np.inf). Shape != grid.shape_yx -> ValueError.
def distance_to_z_interval_um(z_center: np.ndarray | float, z_bottom_um: float,
                              z_top_um: float) -> np.ndarray
# max(z_bottom - z, z - z_top, 0) elementwise, float64; z_top <= z_bottom -> ValueError.
```

## geometry/legacy.py (baseline to keep: binary Erode/Dilate distance measurement)
```python
def legacy_structure(kernel: LegacyKernel) -> np.ndarray
# square: np.ones((size, size), bool). disk: offsets (dy, dx) in [-h, h] with
# dy*dy + dx*dx <= h*(h+1), h = size // 2  (5x5 disk = square without 4 corners).
def legacy_iterations(diameter_um: float, grid: GridSpec, kernel: LegacyKernel) -> int
# requires grid.is_isotropic else InputContractError(field="coordinate_system.pixel_size").
# radius_px = diameter_um / 2 / pixel_size; x = radius_px / (size // 2)
# round -> floor(x + 0.5) (half up, NOT banker's); ceil -> ceil(x - 1e-9); floor -> floor(x + 1e-9)
def legacy_effective_radius_um(diameter_um: float, grid: GridSpec, kernel: LegacyKernel) -> float
# iterations * (size // 2) * pixel_size
def legacy_dilate(mask: np.ndarray, diameter_um: float, grid: GridSpec,
                  kernel: LegacyKernel) -> np.ndarray
# bool; n == 0 -> mask.copy(); else scipy.ndimage.binary_dilation(mask,
# structure=legacy_structure(kernel), iterations=n, border_value=0)
def legacy_erode(mask, diameter_um, grid, kernel) -> np.ndarray   # same with binary_erosion
def legacy_short_hotspot(mask_a: np.ndarray, mask_b: np.ndarray, diameter_um: float,
                         grid: GridSpec, kernel: LegacyKernel) -> np.ndarray
# legacy_dilate(a) & legacy_dilate(b)  (vertical pairs use the same 2-D projection)
def edt_short_hotspot(mask_a, mask_b, radius_um: float, grid: GridSpec) -> np.ndarray
# (distance_to_mask_um(a) <= radius_um) & (distance_to_mask_um(b) <= radius_um)
```

## validation/regression.py
```python
CAUSE_CODES = {"same": 0, "kernel_quantization": 1, "kernel_shape": 2, "unexplained": 3}
@dataclass(frozen=True)
class PixelComparison:
    equal: bool
    n_diff: int
    only_in_reference: int
    only_in_candidate: int
    cause_counts: dict[str, int]     # keys of CAUSE_CODES except "same"
    cause_map: np.ndarray            # uint8 (y, x) codes
def compare_binary_maps(reference: np.ndarray, candidate: np.ndarray) -> PixelComparison
# shape mismatch -> ValueError; differing pixels coded "unexplained".
def explain_legacy_vs_edt(mask_a, mask_b, diameter_um: float, grid: GridSpec,
                          kernel: LegacyKernel) -> PixelComparison
```
`explain_legacy_vs_edt`: ref = legacy_short_hotspot; cand = edt_short_hotspot(r = d/2);
q = edt_short_hotspot(r_eff = legacy_effective_radius_um). For each pixel with ref != cand:
if cand != q -> kernel_quantization, elif ref != q -> kernel_shape, else unexplained.

Tests (tests/unit/test_legacy_distance.py):
- EDT: inside 0; 3-4-5 distance (points (0,0) and (3,4) px, pixel 1 um -> 5.0);
  anisotropic sampling (py=2, px=1); empty mask inf; full mask zeros; z interval inside/below/above.
- iteration rules: pixel 0.5 um, size 5: d=2 -> radius_px 2 -> n=1 all rules;
  d=3 -> x=1.5 -> round 2, ceil 2, floor 1; d=2.5 -> x=1.25 -> round 1, ceil 2, floor 1.
- square kernel n iterations == `distance_transform_cdt(~m, metric="chessboard") <= 2n`.
- disk structure equals the 5x5 corner-less pattern.
- anisotropic grid -> InputContractError.
- parallel stripes: explain_legacy_vs_edt has cause_counts["unexplained"] == 0 and
  equal iff the diameter makes r_eff == r with straight stripes (e.g. pixel 0.5, d=2).
- compare_binary_maps counts only_in_reference / only_in_candidate correctly.
