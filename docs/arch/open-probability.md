# arch: open, probability and risk integration

`failure/open_failure.py` models open as loss of connectivity after a particle blocks an electrode.
`resolve_terminal_regions` selects two terminal masks, legacy terminal points, or an explicit ROI-boundary
proxy. A terminal remains connected when at least one pixel survives. Full coverage of either terminal or
loss of a 4-connected path between the surviving terminal regions is open.

Before any particle is placed, `baseline_status` checks the path. A terminal with no electrode pixel is
`terminal_lost`; terminals that are already disconnected are `baseline_open`. Nominal inputs reject both
(`check_open_baselines`); transformed scenarios record them with an all-zero map, so a pinched line or a
terminal moved off its electrode never turns every candidate into an open.

`failure/auto_open.py` derives terminals without user input: every effective contact hole on the electrode
and every run of electrode pixels on the ROI boundary. Terminals connected in the particle-free state form
groups; a particle opens the electrode when it splits a group or fully covers a boundary run. A fully covered
contact drops out because it is a CNO, and fragments without a terminal are ignored. The unit is
`auto_open:<electrode>` for layers listed in `open_applicable_steps_by_layer`; electrodes with an explicit
`open_path` keep that path. Fewer than two connected terminals gives `insufficient_terminals`.

Automatic open uses exact 9×9 candidate-block certificates from the union and intersection of the
actual particle footprint; ambiguous blocks use the original per-centre local judgement. Set
`PARTICLE_SIM_OPEN_BACKEND=local` to bypass block preparation and use the per-centre route for
automatic open. The default and `auto` setting use certificates; this setting also controls the
existing explicit-open backend route.

`failure/contact_not_open.py` evaluates CNO per contact link `cno:<contact>:<lower>__<upper>`. The
provisional rule (`full_cover_any_hole`) fails a hole when the particle disc covers every effective contact
pixel (hole ∩ lower ∩ upper electrode) and fails the link when any hole fails. A link that disappears in a
scenario is `contact_lost`.

`OpenResult` records `path_id`, electrode, binary failure map, classification, terminal kind and candidate
count. `units.unit_specs` lists pairs, explicit open paths, automatic open electrodes and CNO links;
`units.evaluate_units` evaluates them on one failure-unit axis. Open and CNO units are evaluated only for their
applicable source steps; other steps retain a zero row/map. Scenarios reuse the nominal unit list.
For nominal runs, the binary open result is cached by path, diameter and contact margin while geometry stays
identical, then reused across applicable source steps. The cache stores the wrapped unit result for explicit open,
automatic open and CNO units. Probability runs deliberately do not retain an unbounded cache across transformed scenarios;
bounded or disk-backed reuse is required before that optimization is enabled.

`simulate_open` keeps candidate centres and the particle footprint in the full input coordinate system, but
labels connectivity only inside the tight bounding box of the target electrode. Pixels outside that box are
non-conductive and cannot join the terminal regions. A full-grid reference test and preserved 2B artifacts
verify pixel-identical failure maps, CA tables and union maps.

`geometry/transforms.py` applies one layer overlay to the layer mask, all electrodes and terminal masks.
CD bias changes layer and electrode boundaries; terminal regions shift but are not expanded. Contact masks use
the overlay and CD of their `hole_layer`, and the contact topology is rebuilt non-strictly per scenario.

`pipeline.run_probability` creates Sobol scenarios for each requested scramble seed, transforms masks and
stack, evaluates units, and accumulates per-unit and OR-union maps. It computes nominal, mean, p05, p50,
p95, `1/k` allocated and unique CA. Scenario conditions classified as `baseline_open`, `terminal_lost`,
`contact_lost` or `insufficient_terminals` are counted in `run_manifest.json` `scenario_exclusions`. Convergence rows include the seed. If eager transformed masks are
estimated above 256 MB, bool masks are packbits-encoded in a temporary cache.

`aggregation/risk_tables.py` applies Fab density after simulation. It creates layer, pair, open-path, contact and
electrode tables plus a ranked CSV. CNO contributions are split equally between the lower and upper electrode.
Allocated pair/open/contact contributions sum to source-layer union risk, and
allocated electrode contributions sum to the same total. `risk_manifest.json` records units, semantics,
source paths and artifact hashes.
