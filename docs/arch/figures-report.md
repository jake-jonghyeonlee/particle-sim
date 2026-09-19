# arch: figures-report (T10)

Files: `src/particle_sim/visualization/hotspots.py`, `src/particle_sim/visualization/reports.py`,
`src/particle_sim/cli.py` (add subcommand), tests `tests/unit/test_figures_report.py`.
Depends on T8 run directory: `hotspots.zarr` (arrays legacy_failure, nominal_failure uint8,
pair_margin_um float32 [S, D, P, Y, X]; attrs source_steps, particle_sizes_um, pair_ids,
pair_modes), `ca.parquet` (CA_COLUMNS incl. ca_legacy_um2, ca_nominal_um2),
`ca_union.parquet` (failure_mode, ca_union_nominal_um2, ca_pair_sum_nominal_um2),
`legacy_vs_edt.csv`, `run_manifest.json` (run_id, warnings, assumption_set,
simulator_version, config_hash). Uses `load_manifest`, `load_product_masks`,
`open_hotspot_store`.

matplotlib only with `matplotlib.use("Agg")` before pyplot import; every figure closed
after `savefig(dpi=100)`. Korean text is NOT required inside figures (font issues);
titles/labels in English.

## visualization/hotspots.py
```python
LAYER_COLORS: tuple[str, ...] = ("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")
def plot_layout_overlay(masks: MaskSet, layer_order: list[str], path: Path) -> Path
# layers drawn in order with alpha 0.5, legend by layer name, axes in um (extent from grid)
def plot_pair_panels(legacy: np.ndarray, nominal: np.ndarray, margin_um: np.ndarray,
                     grid: GridSpec, title: str, path: Path) -> Path
# 3 panels: legacy binary, nominal 2.5D binary (both cmap "gray_r", vmin 0, vmax 1),
# contact margin (cmap "RdBu", symmetric fixed range +-max(finite |margin|, 1e-9),
# non-finite shown as NaN), colorbars, extent in um
def plot_ca_curves(ca: pandas.DataFrame, step_id: str, path: Path) -> Path
# x particle_size_um, y ca_nominal_um2 (solid) and ca_legacy_um2 (dashed) per failure_unit_id
def plot_union_vs_sum(union: pandas.DataFrame, step_id: str, path: Path) -> Path
# per failure_mode: ca_union_nominal_um2 vs ca_pair_sum_nominal_um2 over particle size
```
Every function creates parent dirs and returns the path.

## visualization/reports.py
```python
def make_report(run_dir: Path, manifest_path: Path) -> Path
```
Writes into `run_dir / "report"`: `layout_overlay.png`; for each step:
`ca_curve_{step}.png`, `union_vs_sum_{step}.png`, and for each pair the largest size
`pair_{step}_{pair}_d{size:g}.png`; then `report.md` (Korean) containing:
title, run_id / simulator_version / config_hash / assumption_set; the notice
"결과값은 기하 접촉 위험(geometric contact risk)이며 실제 불량확률이 아님";
a "TODO_PHYSICS / 가정" section listing every run_manifest warning;
legacy vs EDT table (all rows of legacy_vs_edt.csv) with a sentence that differences are
classified as kernel_quantization / kernel_shape;
CA summary table per step x pair for the largest size (legacy, nominal, normalized);
union vs pair-sum table per step x mode for the largest size;
relative image links to every PNG. Existing `report/` directory is replaced only if it
contains nothing but files this function generates (`*.png`, `report.md`), otherwise
InputContractError(field="report"). Returns report.md path.

## cli.py addition
`make-report --run DIR --product-config PATH`; prints JSON `{"status": "ok", "report": str}`.

Tests (synthetic pilot -> run_nominal -> make_report, MPLCONFIGDIR set to tmp_path via
monkeypatch.setenv before import): report.md exists, contains "TODO_PHYSICS",
"geometric contact risk" and "실제 불량확률이 아님", and does not contain the English phrase
"failure probability"; expected PNG count = 1 + 3 steps x (2 + 3 pairs) = 16; each PNG opens
with Pillow; running twice succeeds; a foreign file in report/ -> error; CLI exit 0.
