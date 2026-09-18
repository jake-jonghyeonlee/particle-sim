# SPEC: particle-sim (worker summary)

Purpose: 2.5D particle defect GEOMETRY simulator for display layouts. From per-layer
electrode masks + layer stack thickness it computes, per particle source step and particle
diameter, binary short hotspots (lateral / vertical / diagonal), electrode-pair and union
critical area (CA). Legacy Erode/Dilate (5x5 kernel) hotspots are reproduced as the
regression baseline. Outputs are "geometric contact risk", never "failure probability".

Scope, phase 0-3: manifest + validation, legacy morphology, 2.5D nominal short, union, CA,
Zarr/Parquet/manifest outputs, simcompare export, figures.
Scope, phase 2: open connectivity (terminal based, `open_proxy` fallback), Sobol process
variation with probability hotspots and CA percentiles, convergence prefixes, Fab weighting
(common / actual). Probability values are conditioned on the scenario model: when the input
sigma is not a measured process value the result is not an actual Fab failure probability.
Still OUT: tile + halo for very large ROIs, empirical Fab size samples, electrical
calibration of vertical shorts, `skeleton_cut` open proxy.

Inputs per product (one YAML manifest, see docs/arch/schema-config.md):
- layer masks: original patterned mask of each layer (binary PNG/TIFF, same shape)
- electrode masks: electrode-separated masks per layer; their union MUST equal the layer mask
- stack: ordered films (conductor/inorganic/organic/other) with thickness; z is cumulative
- pixel size (um/px), electrode->layer/net, electrode pairs, particle source steps

Physics rules (never hardcode, never silently default):
- Units: um, um^2 (mm^2 only in tables). Image index (y, x). Physical (x_um, y_um, z_um).
- Particle = sphere diameter d, radius r=d/2. Center z = surface_z + r - embedding_depth.
- Layers at/below the landing surface: contact iff r + contact_margin - hypot(Dxy, Dz) >= 0.
- Layers formed AFTER the source step use post_step_policy.mode:
  exclude (never contact) | protrusion (buried-height model) | wrap_contact (conservative).
- Short pairs require different net_id. Lateral: same layer. Vertical: different layers.
- Vertical results are labelled potential_vertical_short (not confirmed shorts).
- Union: OR of pair binary maps inside the same (step, size, mode); never sum CA numbers.
- Unknown physics values are REQUIRED config fields; record them as TODO_PHYSICS warnings.

Hard code rules:
- Python 3.11+, numpy, scipy, pandas, pydantic v2, pyyaml, pillow, tifffile, zarr v3,
  matplotlib (Agg), scikit-image, pyarrow. No network. CPU only.
- Package: src/particle_sim/... Import style `from particle_sim.io.masks import ...`.
- Each module contract lives in docs/arch/<task>.md; follow names/signatures exactly.
- Invalid input: raise particle_sim.errors.InputContractError with product/file/field.
  Never auto-correct, resize, blur or fill missing values.
- User-facing text (error messages, CLI help, warnings) is Korean. Identifiers, JSON keys,
  CSV/Parquet columns and enum values are English.
- Type hints everywhere. Code must pass `ruff check` (line length 100).
- Tests: pytest under tests/unit/, small synthetic arrays (<=200x200), seconds on CPU,
  use tmp_path for files, set MPLCONFIGDIR to tmp_path when importing matplotlib.
- Do not modify files outside the task write_paths.
