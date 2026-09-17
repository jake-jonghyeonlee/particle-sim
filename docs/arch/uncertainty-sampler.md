# arch: uncertainty-sampler (P2-T1)

Files: `src/particle_sim/uncertainty/__init__.py` (empty),
`src/particle_sim/uncertainty/transforms.py`, `src/particle_sim/uncertainty/sampler.py`,
tests `tests/unit/test_uncertainty_sampler.py`.
Depends only on `particle_sim.errors.InputContractError` and these schema models
(already in `src/particle_sim/schema.py`, read them): `NormalSpec(mean, sigma)`,
`UniformSpec(low, high)`, `ScenarioConfig(schema_version, assumption_set, evidence,
thickness_delta_um: dict[str, NormalSpec], overlay_um: dict[str, NormalSpec],
cd_bias_um: dict[str, NormalSpec], particle_diameter_delta_um: NormalSpec | None,
center_offset_um: NormalSpec | None, embedding_depth_um: UniformSpec | None,
planarization_factor: dict[str, UniformSpec], contact_margin_um: UniformSpec | None)`
with `is_degenerate()`. Do NOT import io/, geometry/, failure/ or edit any other file.

## transforms.py
```python
def normal_from_unit(u: np.ndarray | float, spec: NormalSpec) -> np.ndarray
# inverse CDF: mean + sigma * sqrt(2) * erfinv(2u - 1) via scipy.stats.norm.ppf.
# sigma == 0 -> constant mean (no ppf call, so u == 0 or 1 is safe).
# u outside [0, 1] -> InputContractError(field="unit_sample").
def uniform_from_unit(u: np.ndarray | float, spec: UniformSpec) -> np.ndarray
# low + (high - low) * u; low == high -> constant.
def clamp(values, low: float | None, high: float | None) -> np.ndarray
```
All return float64 arrays (0-d array for scalar input is fine); never mutate the input.

## sampler.py
```python
@dataclass(frozen=True)
class Scenario:
    scenario_id: int
    thickness_delta_um: dict[str, float]      # layer_name -> delta
    overlay_xy_um: dict[str, tuple[float, float]]
    cd_bias_um: dict[str, float]
    particle_diameter_delta_um: float
    center_offset_xy_um: tuple[float, float]
    embedding_depth_um: float | None          # None -> use the manifest value
    planarization_factor: dict[str, float]
    contact_margin_um: float | None           # None -> use the manifest value

def variable_names(config: ScenarioConfig) -> list[str]
```
Deterministic dimension order, built from sorted dict keys so the same config always maps
to the same Sobol columns:
`[f"thickness_delta_um.{k}" ...], [f"overlay_um.{k}.x", f"overlay_um.{k}.y" ...],`
`[f"cd_bias_um.{k}" ...], ["particle_diameter_delta_um"], ["center_offset_um.x",`
`"center_offset_um.y"], ["embedding_depth_um"], [f"planarization_factor.{k}" ...],`
`["contact_margin_um"]` — only for the entries present in the config.

```python
def make_sobol_unit_samples(dimensions: int, exponent_m: int, seed: int) -> np.ndarray
# scipy.stats.qmc.Sobol(d=dimensions, scramble=True, seed=seed).random_base2(m=exponent_m)
# dimensions <= 0 or exponent_m < 0 -> InputContractError(field="sobol").
# dimensions == 0 is impossible; a degenerate config is handled by build_scenarios.
def build_scenarios(config: ScenarioConfig, exponent_m: int, seed: int) -> list[Scenario]
```
- `config.is_degenerate()` -> return exactly one Scenario (`scenario_id=0`) holding the fixed
  values (mean / low), regardless of `exponent_m`. This keeps sigma=0 identical to nominal.
- otherwise: one Scenario per Sobol row, `scenario_id` = row index, each column converted with
  `normal_from_unit` / `uniform_from_unit` in `variable_names` order.
- `embedding_depth_um` / `contact_margin_um` stay None when the config omits them.
- Same (config, exponent_m, seed) must produce identical numbers on every call.

Tests (tests/unit/test_uncertainty_sampler.py):
- `normal_from_unit`: u=0.5 -> mean; sigma=0 -> mean for u in {0.0, 0.5, 1.0};
  symmetric values around the mean; out-of-range u raises.
- `uniform_from_unit`: endpoints and midpoint; low == high constant.
- `variable_names` order is deterministic and matches a hand-written expected list.
- `make_sobol_unit_samples(3, 6, 42).shape == (64, 3)`, all values in [0, 1),
  same seed reproduces, different seed differs.
- `build_scenarios` with a degenerate config returns 1 scenario with the fixed values.
- `build_scenarios(config, 6, 7)` returns 64 scenarios, ids 0..63, reproducible, and the
  sample mean of a `NormalSpec(mean=0.0, sigma=1.0)` variable is within 0.2 of 0.0.
- scenario dicts only contain the layer keys named in the config.
