"""Repeat a nominal run and preserve timing/RSS evidence without keeping temporary outputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path


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


def _run_once(manifest: Path) -> dict[str, object]:
    with tempfile.TemporaryDirectory(prefix="particle-sim-benchmark-") as temporary:
        output = Path(temporary) / "run"
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
            raise RuntimeError(
                f"benchmark command failed ({process.returncode}): {stderr or stdout}"
            )
        run_manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
        return {
            "wall_seconds": wall_seconds,
            "peak_rss_bytes": peak_rss,
            "output_bytes": _directory_bytes(output),
            "config_hash": run_manifest["config_hash"],
            "input_hashes": run_manifest["input_hashes"],
        }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be at least 1")

    manifest = args.manifest.resolve()
    warmup = _run_once(manifest)
    runs = [_run_once(manifest) for _ in range(args.repeats)]
    wall = [float(run["wall_seconds"]) for run in runs]
    rss = [int(run["peak_rss_bytes"]) for run in runs]
    sizes = [int(run["output_bytes"]) for run in runs]
    mean = statistics.mean(wall)
    result = {
        "label": args.label,
        "manifest": str(manifest),
        "manifest_sha256": hashlib.sha256(manifest.read_bytes()).hexdigest(),
        "python": sys.version,
        "platform": platform.platform(),
        "warmup": warmup,
        "runs": runs,
        "summary": {
            "repeats": args.repeats,
            "wall_median_seconds": statistics.median(wall),
            "wall_min_seconds": min(wall),
            "wall_max_seconds": max(wall),
            "wall_mean_seconds": mean,
            "wall_cv": statistics.pstdev(wall) / mean if mean else 0.0,
            "peak_rss_max_bytes": max(rss),
            "output_bytes_median": int(statistics.median(sizes)),
        },
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(result["summary"], ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
