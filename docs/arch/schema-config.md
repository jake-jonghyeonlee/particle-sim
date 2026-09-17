# arch: schema-config (T1)

Files: `src/particle_sim/__init__.py` (`__version__ = "0.1.0"`), `errors.py`, `schema.py`,
`config.py`, empty `io/__init__.py`, `geometry/__init__.py`, `failure/__init__.py`,
`aggregation/__init__.py`, `validation/__init__.py`, `visualization/__init__.py`, `export/__init__.py`.

## errors.py
```python
class InputContractError(ValueError):
    def __init__(self, message: str, *, product_id: str | None = None,
                 file: str | None = None, field: str | None = None) -> None
    # attributes: message, product_id, file, field
    # str(e) == "[product=P][file=F][field=X] message" (omit missing brackets)
PARTICLE_SIZES_UM: tuple[float, ...]  # in schema.py: 1..30 step 1, then 40, 50 (floats)
```

## schema.py (pydantic v2, every model `model_config = ConfigDict(extra="forbid", frozen=True)`)
```python
MaterialType = Literal["conductor", "inorganic", "organic", "other"]
SurfaceModel = Literal["flat", "conformal", "partial"]
FailureMode = Literal["lateral_short", "vertical_short"]

class CoordinateSystem(BaseModel):
    origin: Literal["upper_left"]; x_direction: Literal["right"]; y_direction: Literal["down"]
    pixel_size_x_um: float  # >0
    pixel_size_y_um: float  # >0
    roi_width_um: float     # >0
    roi_height_um: float    # >0
class StackLayer(BaseModel):
    layer_name: str; process_order: int  # >=0
    material_type: MaterialType
    thickness_um: float                   # >0
    thickness_sigma_um: float | None = None  # >=0, unused in this phase
    planarization_factor: float | None = None  # 0..1; 0=conformal keeps topography, 1=flat
    z_bottom_um: float | None = None      # optional cross-check only
    patterned: bool                        # True -> must appear in layer_masks
class LayerMask(BaseModel): layer_name: str; path: str
class Electrode(BaseModel):
    electrode_id: str; layer_name: str; net_id: str; electrode_role: str; path: str
class ParticleSourceStep(BaseModel):
    step_id: str
    lands_after_layer: str        # particle lands on the surface right after this layer
    surface_model: SurfaceModel
    fab_distribution_layer: str
class ElectrodePair(BaseModel):
    pair_id: str; electrode_a: str; electrode_b: str; failure_mode: FailureMode
class LegacyKernel(BaseModel):
    size: int = 5                                  # odd >= 3
    shape: Literal["square", "disk"] = "square"
    iterations: Literal["round", "ceil", "floor"] = "round"
class PostStepPolicy(BaseModel):
    mode: Literal["exclude", "protrusion", "wrap_contact"]   # REQUIRED, no default
    coverage_break_factor: float | None = None  # >0, REQUIRED when mode == "protrusion"
class ProductManifest(BaseModel):
    schema_version: Literal["1.0"]
    product_id: str; product_family: str; design_version: str; design_lineage_id: str
    fab_id: str; is_synthetic: bool
    coordinate_system: CoordinateSystem
    stack: list[StackLayer]; layer_masks: list[LayerMask]; electrodes: list[Electrode]
    particle_source_steps: list[ParticleSourceStep]; electrode_pairs: list[ElectrodePair]
    particle_sizes_um: list[float] = list(PARTICLE_SIZES_UM)  # >0, strictly increasing
    embedding_depth_um: float      # REQUIRED, >=0
    contact_margin_um: float       # REQUIRED (may be 0 or negative)
    legacy_kernel: LegacyKernel = LegacyKernel()
    post_step_policy: PostStepPolicy
```
`model_validator(mode="after")` on ProductManifest checks (raise ValueError, Korean text):
unique layer_name / process_order / electrode_id / step_id / pair_id; every layer_masks and
electrode layer_name exists in stack; patterned layers have exactly one LayerMask and
non-patterned have none; each electrode layer is patterned and material conductor;
steps reference existing layers; pairs reference existing electrodes, a != b, net_id differs;
lateral_short -> same layer_name; vertical_short -> different layer_name;
partial surface_model -> every non-conductor layer up to lands_after_layer has
planarization_factor; protrusion mode -> every non-conductor layer has planarization_factor.

```python
def todo_physics_warnings(m: ProductManifest) -> list[str]
# always includes Korean notes for: particle conductivity assumed conductive,
# embedding_depth_um, contact_margin_um, post_step_policy mode (and coverage_break_factor),
# every planarization_factor used; format "TODO_PHYSICS: <field>=<value> <note>".
```

## config.py
```python
def load_manifest(path: Path) -> ProductManifest
# yaml.safe_load; non-mapping/YAML error -> InputContractError(file=path).
# pydantic ValidationError -> InputContractError with product_id (if readable),
# file=str(path), field=dotted loc of the FIRST error, Korean message incl. pydantic msg.
def resolve_input_path(manifest_path: Path, relative: str, *, product_id: str, field: str) -> Path
# relative to manifest directory; absolute paths or resolving outside that dir -> error.
def canonical_hash(obj: object) -> str
# sha256 hex of json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str).
def file_sha256(path: Path) -> str
```
Tests: tests/unit/test_schema_config.py — valid manifest (write YAML in tmp_path), each
rule above failing, error string contains product/file/field, path escape rejected.
