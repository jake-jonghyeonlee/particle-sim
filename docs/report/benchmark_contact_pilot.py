"""Contact-derived open and CNO pilot on the 200x200 complex geometry.

The explicit terminal regions of G1_OPEN, D1_OPEN and D2_OPEN are turned into real contacts:
a same-net pad on the next conductor layer and a contact hole over the old terminal region.
The script checks that automatic-terminal open equals the explicit-terminal open except at
positions that fully block a contact (now CNO) or that involve an ROI-boundary terminal, and
measures TAT/RSS against the explicit-terminal run.
"""

from __future__ import annotations

import argparse
import json
import shutil
import statistics
import subprocess
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import zarr
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "docs/report/work/failure-mode-pilot-200/combined-case"
PADS = {  # open path -> (pad layer, hole layer)
    "G1_OPEN": ("M2", "GI"),
    "D1_OPEN": ("M3", "PLN"),
    "D2_OPEN": ("M3", "PLN"),
}


def _read(path: Path) -> np.ndarray:
    return np.asarray(Image.open(path).convert("L")) > 127


def _write(path: Path, mask: np.ndarray) -> None:
    Image.fromarray(mask.astype(np.uint8) * 255).save(path)


def build_case(target: Path) -> Path:
    if target.exists():
        shutil.rmtree(target)
    shutil.copytree(SOURCE, target)
    manifest = yaml.safe_load((SOURCE / "manifest.yaml").read_text(encoding="utf-8"))
    electrodes = {e["electrode_id"]: e for e in manifest["electrodes"]}
    layer_paths = {item["layer_name"]: item["path"] for item in manifest["layer_masks"]}
    holes: dict[str, np.ndarray] = {}
    for open_path in manifest["open_paths"]:
        pad_layer, hole_layer = PADS[open_path["path_id"]]
        wire = electrodes[open_path["electrode"]]
        for side in ("a", "b"):
            region = _read(target / open_path[f"terminal_{side}_mask_path"])
            others = [
                _read(target / e["path"]) for e in manifest["electrodes"]
                if e["layer_name"] == pad_layer
            ]
            if any((region & other).any() for other in others):
                raise SystemExit(f"pad for {open_path['path_id']}:{side} overlaps {pad_layer}")
            pad_id = f"{open_path['electrode']}_PAD_{side.upper()}"
            pad_path = f"masks/electrode_{pad_id}.png"
            _write(target / pad_path, region)
            manifest["electrodes"].append(
                {
                    "electrode_id": pad_id,
                    "layer_name": pad_layer,
                    "net_id": wire["net_id"],
                    "electrode_role": "PAD",
                    "path": pad_path,
                }
            )
            layer_mask = _read(target / layer_paths[pad_layer]) | region
            _write(target / layer_paths[pad_layer], layer_mask)
            holes[hole_layer] = holes.get(hole_layer, np.zeros(region.shape, bool)) | region
    manifest["contacts"] = []
    for hole_layer, mask in holes.items():
        path = f"masks/contact_{hole_layer}.png"
        _write(target / path, mask)
        manifest["contacts"].append(
            {
                "contact_id": f"CNT_{hole_layer}",
                "path": path,
                "hole_layer": hole_layer,
                "applicable_particle_source_steps": ["S_M2"] if hole_layer == "GI" else ["S_PLN"],
            }
        )
    manifest["open_applicable_steps_by_layer"] = {
        "M1": ["S_M1", "S_M2", "S_PLN"],
        "M2": ["S_M2", "S_PLN"],
        "M3": [],
    }
    manifest["open_paths"] = []
    declared = {tuple(sorted((p["electrode_a"], p["electrode_b"]))) for p in manifest["electrode_pairs"]}
    excluded = {
        tuple(sorted((p["electrode_a"], p["electrode_b"])))
        for p in manifest.get("excluded_electrode_pairs", [])
    }
    net = {e["electrode_id"]: e["net_id"] for e in manifest["electrodes"]}
    for a, b in combinations(sorted(net), 2):
        if net[a] != net[b] and (a, b) not in declared | excluded:
            # keep the short workload identical to the explicit-terminal v1 pilot
            reason = (
                "benchmark contact pad" if "_PAD_" in a or "_PAD_" in b
                else "not declared in the v1 pilot"
            )
            manifest.setdefault("excluded_electrode_pairs", []).append(
                {"electrode_a": a, "electrode_b": b, "reason": reason}
            )
    manifest["schema_version"] = "2.0"
    manifest["product_id"] = manifest["product_id"] + "_CONTACT"
    (target / "manifest.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    return target / "manifest.yaml"


def _timed_run(manifest: Path, output: Path) -> dict:
    if output.exists():
        shutil.rmtree(output)
    code = (
        "import sys,time,resource,json;from pathlib import Path;"
        "from particle_sim.pipeline import run_nominal;t=time.perf_counter();"
        "run_nominal(Path(sys.argv[1]),Path(sys.argv[2]));"
        "print(json.dumps({'wall':time.perf_counter()-t,"
        "'rss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))"
    )
    done = subprocess.run(
        [sys.executable, "-c", code, str(manifest), str(output)],
        check=True, capture_output=True, text=True,
    )
    return json.loads(done.stdout.strip().splitlines()[-1])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-root", required=True)
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--repeats", type=int, default=5)
    args = parser.parse_args()
    work = Path(args.work_root)
    work.mkdir(parents=True, exist_ok=True)
    contact_manifest = build_case(work / "contact-case")
    runs = {"explicit_terminals": SOURCE / "manifest.yaml", "contact_terminals": contact_manifest}
    timings: dict[str, list[dict]] = {name: [] for name in runs}
    for name, manifest in runs.items():
        _timed_run(manifest, work / f"{name}-warmup")
    for repeat in range(args.repeats):
        order = list(runs) if repeat % 2 == 0 else list(reversed(runs))
        for name in order:
            timings[name].append(_timed_run(runs[name], work / f"{name}-result"))
    rss_scale = 1 if sys.platform == "darwin" else 1024

    explicit = zarr.open_group(work / "explicit_terminals-result/hotspots.zarr", mode="r")
    contact = zarr.open_group(work / "contact_terminals-result/hotspots.zarr", mode="r")
    e_units = list(explicit.attrs["failure_unit_ids"])
    c_units = list(contact.attrs["failure_unit_ids"])
    e_map = np.asarray(explicit["nominal_failure"], dtype=bool)
    c_map = np.asarray(contact["nominal_failure"], dtype=bool)
    ca = pd.read_parquet(work / "contact_terminals-result/ca.parquet")
    equivalence = {}
    for path_id, electrode in (("G1_OPEN", "G1"), ("D1_OPEN", "D1"), ("D2_OPEN", "D2")):
        e = e_map[:, :, e_units.index(path_id)]
        c = c_map[:, :, c_units.index(f"auto_open:{electrode}")]
        cno = np.zeros_like(c)
        for index, unit in enumerate(c_units):
            if unit.startswith("cno:") and electrode in unit.split(":", 2)[2].split("__"):
                # CNO maps are step-dependent; the geometric blockage is the OR over steps
                cno |= c_map[:, :, index].any(axis=0, keepdims=True)
        equivalence[path_id] = {
            "explicit_open_px": int(e.sum()),
            "auto_open_px": int(c.sum()),
            "auto_only_px": int((c & ~e).sum()),
            "explicit_only_px": int((e & ~c).sum()),
            "explicit_only_outside_cno_px": int((e & ~c & ~cno).sum()),
        }
    summary = {
        name: {
            "repeats": len(items),
            "wall_median_seconds": statistics.median(item["wall"] for item in items),
            "wall_min_seconds": min(item["wall"] for item in items),
            "wall_max_seconds": max(item["wall"] for item in items),
            "peak_rss_max_bytes": max(item["rss"] for item in items) * rss_scale,
        }
        for name, items in timings.items()
    }
    payload = {
        "label": "contact-terminal-open-and-cno-pilot-200x200",
        "purpose": "automatic contact/ROI terminals and CNO vs explicit terminal regions",
        "python": sys.version.split()[0],
        "repeats": args.repeats,
        "tat_summary": summary,
        "open_equivalence": equivalence,
        "unit_counts": ca.drop_duplicates("failure_unit_id")["unit_kind"].value_counts().to_dict(),
        "classifications": ca.groupby("unit_kind")["classification"].value_counts().to_dict(),
        "cno_ca_um2_by_step": ca[ca["unit_kind"] == "contact"]
        .groupby("source_step")["ca_nominal_um2"].sum().to_dict(),
    }
    payload["classifications"] = {"|".join(k): int(v) for k, v in payload["classifications"].items()}
    Path(args.output_json).write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
