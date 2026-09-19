"""Benchmark failure modes on one shared complex synthetic wiring layout.

The benchmark creates four manifests from ``complex_multilayer``:

* lateral short pairs only
* vertical short pairs only
* three open paths only
* all pairs and open paths combined

It reports end-to-end CLI time, peak RSS, time per active failure-unit condition,
and the allocated-CA contribution of each mode from the combined run.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from PIL import Image

from particle_sim.validation.synthetic_cases import SyntheticCase, build_case, write_case

MODE_LABEL = {
    "lateral_short": "lateral short",
    "vertical_short": "vertical short",
    "lateral_open": "open",
}
ISOLATED_MODES = ("lateral_short", "vertical_short", "open")

# Terminal regions are deliberately placed at electrically distant ends of three
# connected, bent/branched electrodes. Coordinates are half-open [y0:y1, x0:x1].
OPEN_PATHS = (
    {
        "path_id": "G1_OPEN",
        "electrode": "G1",
        "terminal_a": (18, 28, 80, 88),
        "terminal_b": (107, 115, 15, 25),
        "applicable_particle_source_steps": ["S_M1", "S_M2", "S_PLN"],
    },
    {
        "path_id": "D1_OPEN",
        "electrode": "D1",
        "terminal_a": (35, 45, 30, 38),
        "terminal_b": (90, 100, 55, 63),
        "applicable_particle_source_steps": ["S_M2", "S_PLN"],
    },
    {
        "path_id": "D2_OPEN",
        "electrode": "D2",
        "terminal_a": (105, 115, 92, 100),
        "terminal_b": (150, 160, 77, 85),
        "applicable_particle_source_steps": ["S_M2", "S_PLN"],
    },
)


def _region(shape: tuple[int, int], box: tuple[int, int, int, int]) -> np.ndarray:
    y0, y1, x0, x1 = box
    result = np.zeros(shape, dtype=bool)
    result[y0:y1, x0:x1] = True
    return result


def _write_binary_mask(path: Path, mask: np.ndarray) -> None:
    Image.fromarray(mask.astype(np.uint8) * 255, mode="L").save(path)


def _open_manifest_rows(open_paths=OPEN_PATHS) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for spec in open_paths:
        path_id = str(spec["path_id"])
        rows.append(
            {
                "path_id": path_id,
                "electrode": spec["electrode"],
                "terminal_a_mask_path": f"masks/terminal_{path_id}_a.png",
                "terminal_b_mask_path": f"masks/terminal_{path_id}_b.png",
                "applicable_particle_source_steps": spec["applicable_particle_source_steps"],
                "fallback": "error",
            }
        )
    return rows


def _pad_case(
    case: SyntheticCase, canvas_size: int
) -> tuple[SyntheticCase, tuple[dict[str, object], ...]]:
    source_shape = next(iter(case.layer_masks.values())).shape
    if source_shape[0] != source_shape[1]:
        raise ValueError(f"padding은 정사각형 fixture만 지원합니다: {source_shape}")
    if canvas_size < source_shape[0]:
        raise ValueError(f"canvas size {canvas_size}가 fixture {source_shape[0]}보다 작습니다")
    if canvas_size == source_shape[0]:
        return case, OPEN_PATHS

    offset = (canvas_size - source_shape[0]) // 2

    def padded(mask: np.ndarray) -> np.ndarray:
        result = np.zeros((canvas_size, canvas_size), dtype=bool)
        result[offset : offset + source_shape[0], offset : offset + source_shape[1]] = mask
        return result

    manifest = deepcopy(case.manifest)
    coordinate = manifest["coordinate_system"]
    coordinate["roi_width_um"] = canvas_size * coordinate["pixel_size_x_um"]
    coordinate["roi_height_um"] = canvas_size * coordinate["pixel_size_y_um"]
    translated_paths = []
    for original in OPEN_PATHS:
        item = dict(original)
        for terminal in ("terminal_a", "terminal_b"):
            y0, y1, x0, x1 = original[terminal]
            item[terminal] = (y0 + offset, y1 + offset, x0 + offset, x1 + offset)
        translated_paths.append(item)
    return (
        SyntheticCase(
            name=f"{case.name}_padded_{canvas_size}",
            manifest=manifest,
            layer_masks={name: padded(mask) for name, mask in case.layer_masks.items()},
            electrode_masks={name: padded(mask) for name, mask in case.electrode_masks.items()},
        ),
        tuple(translated_paths),
    )


def _prepare_cases(work_root: Path, canvas_size: int) -> dict[str, Path]:
    base = build_case("complex_multilayer")
    base, open_path_specs = _pad_case(base, canvas_size)
    all_pairs = list(base.manifest["electrode_pairs"])
    variants = {
        "lateral_short": (
            [pair for pair in all_pairs if pair["failure_mode"] == "lateral_short"],
            [],
        ),
        "vertical_short": (
            [pair for pair in all_pairs if pair["failure_mode"] == "vertical_short"],
            [],
        ),
        "open": ([], _open_manifest_rows(open_path_specs)),
        "combined": (all_pairs, _open_manifest_rows(open_path_specs)),
    }

    manifests: dict[str, Path] = {}
    for name, (pairs, open_paths) in variants.items():
        manifest = deepcopy(base.manifest)
        size_suffix = "" if canvas_size == 200 else f"_{canvas_size}X{canvas_size}"
        manifest["product_id"] = f"SYN_COMPLEX_FAILURE_MODE_{name.upper()}{size_suffix}"
        manifest["electrode_pairs"] = pairs
        manifest["open_paths"] = open_paths
        case = SyntheticCase(
            name=f"complex_failure_mode_{name}",
            manifest=manifest,
            layer_masks=base.layer_masks,
            electrode_masks=base.electrode_masks,
        )
        case_dir = work_root / f"{name}-case"
        manifest_path = write_case(case, case_dir)
        for spec in open_path_specs if open_paths else ():
            electrode = base.electrode_masks[str(spec["electrode"])]
            terminal_a = _region(electrode.shape, spec["terminal_a"])
            terminal_b = _region(electrode.shape, spec["terminal_b"])
            if not np.all(terminal_a <= electrode) or not np.all(terminal_b <= electrode):
                raise RuntimeError(f"{spec['path_id']} 단자 영역이 전극 밖에 있습니다")
            if np.any(terminal_a & terminal_b):
                raise RuntimeError(f"{spec['path_id']} 단자 영역이 서로 겹칩니다")
            path_id = str(spec["path_id"])
            _write_binary_mask(case_dir / "masks" / f"terminal_{path_id}_a.png", terminal_a)
            _write_binary_mask(case_dir / "masks" / f"terminal_{path_id}_b.png", terminal_b)
        manifests[name] = manifest_path
    return manifests


def _rss_bytes(pid: int) -> int:
    result = subprocess.run(
        ["ps", "-o", "rss=", "-p", str(pid)],
        capture_output=True,
        check=False,
        text=True,
    )
    try:
        return int(result.stdout.strip()) * 1024
    except ValueError:
        return 0


def _directory_bytes(path: Path) -> int:
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def _run_cli(manifest: Path, output: Path) -> dict[str, object]:
    command = [
        sys.executable,
        "-m",
        "particle_sim",
        "simulate-nominal",
        "--product-config",
        str(manifest),
        "--output",
        str(output),
    ]
    started = time.perf_counter()
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    peak_rss = 0
    while process.poll() is None:
        peak_rss = max(peak_rss, _rss_bytes(process.pid))
        time.sleep(0.01)
    stdout, stderr = process.communicate()
    wall_seconds = time.perf_counter() - started
    if process.returncode != 0:
        raise RuntimeError(f"benchmark command failed ({process.returncode}): {stderr or stdout}")
    run_manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
    return {
        "wall_seconds": wall_seconds,
        "peak_rss_bytes": peak_rss,
        "output_bytes": _directory_bytes(output),
        "config_hash": run_manifest["config_hash"],
        "input_hashes": run_manifest["input_hashes"],
    }


def _temporary_run(manifest: Path) -> dict[str, object]:
    with tempfile.TemporaryDirectory(prefix="particle-sim-mode-benchmark-") as temporary:
        return _run_cli(manifest, Path(temporary) / "run")


def _validate(manifest: Path) -> None:
    command = [
        sys.executable,
        "-m",
        "particle_sim",
        "validate-input",
        "--product-config",
        str(manifest),
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"input validation failed: {result.stderr or result.stdout}")


def _active_counts(manifests: dict[str, Path]) -> dict[str, dict[str, int]]:
    result: dict[str, dict[str, int]] = {}
    for name, manifest_path in manifests.items():
        manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
        step_count = len(manifest["particle_source_steps"])
        size_count = len(manifest["particle_sizes_um"])
        pair_count = len(manifest["electrode_pairs"])
        open_count = len(manifest["open_paths"])
        active_open_steps = sum(
            len(item.get("applicable_particle_source_steps") or manifest["particle_source_steps"])
            for item in manifest["open_paths"]
        )
        active_evaluations = pair_count * step_count * size_count + active_open_steps * size_count
        result[name] = {
            "electrode_pair_count": pair_count,
            "open_path_count": open_count,
            "source_step_count": step_count,
            "particle_size_count": size_count,
            "active_failure_unit_conditions": active_evaluations,
        }
    return result


def _summarize_runs(
    runs: dict[str, list[dict[str, object]]],
    active: dict[str, dict[str, int]],
) -> dict[str, dict[str, float | int]]:
    summary: dict[str, dict[str, float | int]] = {}
    for name, mode_runs in runs.items():
        wall = [float(run["wall_seconds"]) for run in mode_runs]
        rss = [int(run["peak_rss_bytes"]) for run in mode_runs]
        sizes = [int(run["output_bytes"]) for run in mode_runs]
        median = statistics.median(wall)
        mean = statistics.mean(wall)
        count = active[name]["active_failure_unit_conditions"]
        summary[name] = {
            "repeats": len(mode_runs),
            "wall_median_seconds": median,
            "wall_min_seconds": min(wall),
            "wall_max_seconds": max(wall),
            "wall_mean_seconds": mean,
            "wall_cv": statistics.pstdev(wall) / mean if mean else 0.0,
            "peak_rss_max_bytes": max(rss),
            "output_bytes_median": int(statistics.median(sizes)),
            "active_failure_unit_conditions": count,
            "wall_seconds_per_active_condition": median / count,
        }

    isolated_total = sum(summary[name]["wall_median_seconds"] for name in ISOLATED_MODES)
    for name in ISOLATED_MODES:
        summary[name]["isolated_tat_share"] = summary[name]["wall_median_seconds"] / isolated_total
    summary["combined"]["vs_isolated_sum_ratio"] = (
        summary["combined"]["wall_median_seconds"] / isolated_total
    )
    return summary


def _analyze_combined(output: Path) -> dict[str, object]:
    ca = pd.read_parquet(output / "ca.parquet")
    union = pd.read_parquet(output / "ca_union.parquet")
    ca["report_mode"] = ca["failure_mode"].replace(MODE_LABEL)
    union["report_mode"] = union["failure_mode"].replace(MODE_LABEL)

    applicable = ca[ca["classification"] != "not_applicable_source_step"].copy()
    grouped = (
        applicable.groupby("report_mode", sort=False)
        .agg(
            row_count=("failure_unit_id", "size"),
            nonzero_row_count=("ca_nominal_um2", lambda values: int((values > 0).sum())),
            raw_ca_sum_um2=("ca_nominal_um2", "sum"),
            allocated_ca_sum_um2=("ca_allocated_nominal_um2", "sum"),
            unique_ca_sum_um2=("ca_unique_nominal_um2", "sum"),
        )
        .reset_index()
    )
    allocated_total = float(grouped["allocated_ca_sum_um2"].sum())
    grouped["allocated_ca_share"] = grouped["allocated_ca_sum_um2"] / allocated_total

    any_union = union[union["failure_mode"] == "any"]
    any_union_total = float(any_union["ca_union_nominal_um2"].sum())
    allocated_delta = allocated_total - any_union_total

    mode_union = (
        union[union["failure_mode"].isin(MODE_LABEL)]
        .groupby("report_mode", sort=False)["ca_union_nominal_um2"]
        .sum()
        .rename("mode_union_ca_sum_um2")
        .reset_index()
    )
    grouped = grouped.merge(mode_union, on="report_mode", how="left")

    by_size = (
        applicable.groupby(["particle_size_um", "report_mode"], sort=True)[
            "ca_allocated_nominal_um2"
        ]
        .sum()
        .rename("allocated_ca_um2")
        .reset_index()
    )
    size_totals = by_size.groupby("particle_size_um")["allocated_ca_um2"].transform("sum")
    by_size["allocated_ca_share"] = by_size["allocated_ca_um2"] / size_totals

    return {
        "mode_contribution": grouped.to_dict(orient="records"),
        "mode_contribution_by_particle_size": by_size.to_dict(orient="records"),
        "allocated_ca_total_um2": allocated_total,
        "any_union_ca_total_um2": any_union_total,
        "allocation_minus_union_um2": allocated_delta,
        "allocation_matches_union": abs(allocated_delta) <= 1e-8,
        "share_basis": (
            "combined run ca_allocated_nominal_um2; overlapping failure pixels are split 1/k"
        ),
    }


def _write_summary_csv(
    path: Path,
    summary: dict[str, dict[str, float | int]],
) -> None:
    rows = [{"mode": mode, **values} for mode, values in summary.items()]
    pd.DataFrame(rows).to_csv(path, index=False)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--canvas-size", type=int, default=200)
    parser.add_argument(
        "--reset",
        action="store_true",
        help="전용 work-root를 지우고 fixture와 보존 결과를 다시 만든다",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="기존 output JSON의 실행을 재사용하고 목표 반복 수까지 이어서 측정한다",
    )
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be at least 1")
    if args.canvas_size < 200:
        parser.error("--canvas-size must be at least 200")
    if args.reset and args.resume:
        parser.error("--reset과 --resume은 함께 쓸 수 없습니다")

    work_root = args.work_root.resolve()
    if args.reset and work_root.exists():
        shutil.rmtree(work_root)
    work_root.mkdir(parents=True, exist_ok=True)
    manifests = _prepare_cases(work_root, args.canvas_size)
    for manifest in manifests.values():
        _validate(manifest)

    active = _active_counts(manifests)
    previous: dict[str, object] | None = None
    if args.resume:
        if not args.output_json.exists():
            parser.error("--resume에는 기존 --output-json 파일이 필요합니다")
        previous = json.loads(args.output_json.read_text(encoding="utf-8"))
        for name, manifest in manifests.items():
            expected = hashlib.sha256(manifest.read_bytes()).hexdigest()
            actual = previous["manifests"][name]["sha256"]
            if actual != expected:
                parser.error(f"{name} manifest hash가 이전 측정과 다릅니다")

    if previous is None:
        warmups = {name: _temporary_run(manifest) for name, manifest in manifests.items()}
        runs: dict[str, list[dict[str, object]]] = {name: [] for name in manifests}
    else:
        warmups = previous["warmups"]
        runs = previous["runs"]
        existing_counts = {name: len(runs[name]) for name in manifests}
        if len(set(existing_counts.values())) != 1:
            parser.error(f"기존 mode별 반복 수가 서로 다릅니다: {existing_counts}")
        if next(iter(existing_counts.values())) > args.repeats:
            parser.error("기존 반복 수가 요청한 --repeats보다 큽니다")

    names = list(manifests)
    completed_repeats = len(runs[names[0]])
    for repeat_index in range(completed_repeats, args.repeats):
        order = names[repeat_index % len(names) :] + names[: repeat_index % len(names)]
        for name in order:
            runs[name].append(_temporary_run(manifests[name]))

    combined_output = work_root / "combined-result"
    if previous is not None and combined_output.exists():
        retained_combined_run = previous["combined_result"]["retained_run_measurement"]
    else:
        if combined_output.exists():
            shutil.rmtree(combined_output)
        retained_combined_run = _run_cli(manifests["combined"], combined_output)
    tat_summary = _summarize_runs(runs, active)
    result = {
        "label": f"complex-failure-mode-pilot-{args.canvas_size}x{args.canvas_size}",
        "canvas_size_yx": [args.canvas_size, args.canvas_size],
        "purpose": "relative TAT and geometric-risk mix; synthetic data, not Fab probability",
        "python": sys.version,
        "platform": platform.platform(),
        "repeats": args.repeats,
        "manifests": {
            name: {
                "path": str(path),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
            for name, path in manifests.items()
        },
        "active_counts": active,
        "warmups": warmups,
        "runs": runs,
        "tat_summary": tat_summary,
        "combined_result": {
            "path": str(combined_output),
            "retained_run_measurement": retained_combined_run,
            **_analyze_combined(combined_output),
        },
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    _write_summary_csv(args.output_csv, tat_summary)
    print(json.dumps({"tat_summary": tat_summary}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
