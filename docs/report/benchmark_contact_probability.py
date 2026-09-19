"""Sobol probability pilot on the 200x200 contact case (assumed variation, synthetic).

Runs simulate-probability once per case (explicit terminals vs contact terminals) in a separate
process and records wall time, peak RSS, CA by unit kind and scenario exclusions. The overlay/CD
sigmas are assumptions for a TAT check, not measured Fab data.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
EXPLICIT = ROOT / "docs/report/work/failure-mode-pilot-200/combined-case/manifest.yaml"
CONTACT = ROOT / "docs/report/work/contact-pilot-200/contact-case/manifest.yaml"
SCENARIO = """schema_version: "2.0"
assumption_set: contact_probability_pilot
evidence: assumed
overlay_reference_layer: M1
overlay_sign_convention: positive_x_right_positive_y_down
overlay_um:
{overlay}
cd_bias_um:
{cd}
"""


def _scenario(path: Path, layers: list[str]) -> Path:
    overlay = "\n".join(
        f"  {layer}:\n    x: {{mean: 0.0, sigma: 0.2, evidence: assumed}}\n"
        f"    y: {{mean: 0.0, sigma: 0.2, evidence: assumed}}"
        for layer in layers
        if layer != "M1"
    )
    cd = "\n".join(
        f"  {layer}: {{mean: 0.0, sigma: 0.1, evidence: assumed}}" for layer in layers
    )
    path.write_text(SCENARIO.format(overlay=overlay, cd=cd), encoding="utf-8")
    return path


def _run(manifest: Path, scenario: Path, output: Path, sobol_m: int, seed: int) -> dict:
    if output.exists():
        shutil.rmtree(output)
    code = (
        "import sys,time,resource,json;from pathlib import Path;"
        "from particle_sim.pipeline import run_probability;t=time.perf_counter();"
        "run_probability(Path(sys.argv[1]),Path(sys.argv[2]),Path(sys.argv[3]),"
        "sobol_exponent_m=int(sys.argv[4]),seed=int(sys.argv[5]));"
        "print(json.dumps({'wall':time.perf_counter()-t,"
        "'rss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))"
    )
    done = subprocess.run(
        [sys.executable, "-c", code, str(manifest), str(scenario), str(output), str(sobol_m),
         str(seed)],
        check=True, capture_output=True, text=True,
    )
    measured = json.loads(done.stdout.strip().splitlines()[-1])
    ca = pd.read_parquet(output / "ca.parquet")
    run_manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
    rss_scale = 1 if sys.platform == "darwin" else 1024
    return {
        "wall_seconds": measured["wall"],
        "peak_rss_bytes": measured["rss"] * rss_scale,
        "scenario_count": run_manifest["scenario_count"],
        "scenario_mask_storage": run_manifest["scenario_mask_storage"],
        "ca_mean_um2_by_unit_kind": ca.groupby("unit_kind")["ca_mean_um2"].sum().to_dict(),
        "ca_nominal_um2_by_unit_kind": ca.groupby("unit_kind")["ca_nominal_um2"].sum().to_dict(),
        "scenario_exclusions": run_manifest.get("scenario_exclusions", []),
        "scenario_contact_lost_holes": run_manifest.get("scenario_contact_lost_holes"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-root", required=True)
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--sobol-m", type=int, default=6)
    parser.add_argument("--seed", type=int, default=20260920)
    args = parser.parse_args()
    work = Path(args.work_root)
    work.mkdir(parents=True, exist_ok=True)
    if not CONTACT.is_file():
        raise SystemExit("run docs/report/benchmark_contact_pilot.py first")
    results = {}
    for name, manifest, layers in (
        ("explicit_terminals", EXPLICIT, ["M1", "M2", "M3"]),
        ("contact_terminals", CONTACT, ["M1", "M2", "M3", "GI", "PLN"]),
    ):
        scenario = _scenario(work / f"{name}-scenario.yaml", layers)
        results[name] = _run(manifest, scenario, work / f"{name}-prob", args.sobol_m, args.seed)
        print(name, round(results[name]["wall_seconds"], 1), flush=True)
    payload = {
        "label": "contact-probability-pilot-200x200",
        "purpose": "TAT and exclusion bookkeeping of automatic open and CNO under assumed variation",
        "python": sys.version.split()[0],
        "sobol_exponent_m": args.sobol_m,
        "seed": args.seed,
        "repeats": 1,
        "variation": "assumed: overlay x/y sigma 0.2 um (non-reference layers), CD sigma 0.1 um",
        "results": results,
    }
    Path(args.output_json).write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
