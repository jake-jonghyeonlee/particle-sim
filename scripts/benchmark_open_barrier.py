"""Opt-in explicit-open equality and timing harness.

Examples:
  python scripts/benchmark_open_barrier.py synthetic
  python scripts/benchmark_open_barrier.py equality --manifest /path/product.yaml
  python scripts/benchmark_open_barrier.py benchmark --manifest /path/product.yaml
  python scripts/benchmark_open_barrier.py explicit32 --manifest /path/product.yaml
  python scripts/benchmark_open_barrier.py full32 --manifest /path/product.yaml
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import resource
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import yaml
import zarr

from particle_sim.config import load_manifest
from particle_sim.failure import open_failure as of
from particle_sim.failure.open_barrier import PreparedBarrier, prepare_barrier
from particle_sim.geometry.contact import DistanceCache
from particle_sim.io.masks import load_product_masks
from particle_sim.io.stack import build_stack_z
from particle_sim.units import KIND_OPEN_PATH, evaluate_units, unit_specs

SIZES = [*range(1, 31), 40, 50]
COMPARE = [2, 5, 10, 20, 30, 50]


def _rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value * (1024 if sys.platform.startswith("linux") else 1))


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _record(**values) -> None:
    print(json.dumps(values, sort_keys=True), flush=True)


def _copy_manifest(source: Path, destination: Path) -> None:
    data = yaml.safe_load(source.read_text())

    def copy_inputs(value):
        if isinstance(value, dict):
            for key, item in value.items():
                if key == "path" and isinstance(item, str):
                    path = Path(item)
                    if path.is_absolute() or ".." in path.parts:
                        raise ValueError(f"input path is not relative: {item}")
                    target = destination.parent / path
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(source.parent / path, target)
                else:
                    copy_inputs(item)
        elif isinstance(value, list):
            for item in value:
                copy_inputs(item)

    copy_inputs(data)
    data["particle_sizes_um"] = SIZES
    destination.write_text(yaml.safe_dump(data, sort_keys=False))


def _meta(source: Path, masks) -> None:
    open_masks = {
        name: hashlib.sha256(masks.electrodes[name].tobytes()).hexdigest()
        for name in sorted({path.electrode for path in load_manifest(source).open_paths})
    }
    _record(kind="meta", manifest=str(source), manifest_sha256=_sha(source),
            open_electrode_hashes=open_masks, python=sys.version.split()[0],
            numpy=np.__version__, scipy=scipy.__version__, hardware=platform.platform(),
            cpu=platform.processor(), pid=os.getpid())


def _measure(manifest, masks, path, diameter, backend):
    os.environ["PARTICLE_SIM_OPEN_BACKEND"] = backend
    stats = {"preparation_calls": 0, "preparation_s": 0.0,
             "candidate_s": 0.0, "local_checks": 0, "local_s": 0.0}
    original_prepare = of.prepare_barrier
    original_candidates = of.open_candidates
    original_local = of._evaluate_local_centres

    def prepare(*args):
        start = time.perf_counter()
        answer = original_prepare(*args)
        stats["preparation_calls"] += 1
        stats["preparation_s"] += time.perf_counter() - start
        return answer

    def candidates(*args):
        start = time.perf_counter()
        answer = original_candidates(*args)
        stats["candidate_s"] += time.perf_counter() - start
        return answer

    def local(*args):
        start = time.perf_counter()
        stats["local_checks"] += int(args[3].sum())
        answer = original_local(*args)
        stats["local_s"] += time.perf_counter() - start
        return answer

    of.prepare_barrier = prepare
    of.open_candidates = candidates
    of._evaluate_local_centres = local
    start_cpu = time.process_time()
    start = time.perf_counter()
    try:
        result = of.simulate_open(manifest, masks, path, diameter)
        elapsed = time.perf_counter() - start
        cpu = time.process_time() - start_cpu
    finally:
        of.prepare_barrier = original_prepare
        of.open_candidates = original_candidates
        of._evaluate_local_centres = original_local
    return result, elapsed, cpu, stats


def _equality(manifest, masks, sizes):
    for path in manifest.open_paths:
        electrode = masks.electrodes[path.electrode]
        ta, tb, _, _ = of.resolve_terminal_regions(
            path, masks, electrode, manifest.product_id, require_on_electrode=False
        )
        prepared = prepare_barrier(electrode, ta, tb, masks.grid.pixel_size_y_um,
                                   masks.grid.pixel_size_x_um)
        _record(kind="certificate", path=path.path_id,
                eligible=isinstance(prepared, PreparedBarrier),
                reason=None if isinstance(prepared, PreparedBarrier) else prepared.reason)
        for diameter in sizes:
            local, local_s, local_cpu, local_stats = _measure(
                manifest, masks, path, diameter, "local"
            )
            fast, fast_s, fast_cpu, fast_stats = _measure(
                manifest, masks, path, diameter, "auto"
            )
            differences = np.argwhere(local.failure != fast.failure)
            stamp = of.particle_footprint(masks.grid, diameter / 2,
                                          manifest.contact_margin_um)
            _record(kind="equality", path=path.path_id, diameter_um=diameter,
                    candidates=fast.candidates, stamp_yx=list(stamp.shape),
                    window_yx=[min(electrode.shape[0], stamp.shape[0] * 3),
                               min(electrode.shape[1], stamp.shape[1] * 3)],
                    failure_pixels=int(fast.failure.sum()),
                    local_failure_pixels=int(local.failure.sum()),
                    mismatch_pixels=len(differences), first_mismatches=differences[:10].tolist(),
                    metadata_equal=(local.path_id, local.electrode, local.classification,
                                    local.terminal_a_yx, local.terminal_b_yx,
                                    local.terminal_kind, local.candidates) ==
                                   (fast.path_id, fast.electrode, fast.classification,
                                    fast.terminal_a_yx, fast.terminal_b_yx,
                                    fast.terminal_kind, fast.candidates),
                    local_s=local_s, local_cpu_s=local_cpu,
                    auto_s=fast_s, auto_cpu_s=fast_cpu,
                    local_breakdown=local_stats, auto_breakdown=fast_stats,
                    local_us_per_candidate=1e6 * local_s / max(local.candidates, 1),
                    auto_us_per_candidate=1e6 * fast_s / max(fast.candidates, 1),
                    peak_rss_bytes=_rss_bytes())
            if differences.size:
                raise SystemExit("pixel mismatch")


def _benchmark(manifest, masks):
    for path in manifest.open_paths:
        for diameter in COMPARE:
            _, elapsed, cpu, stats = _measure(manifest, masks, path, diameter, "auto")
            _record(kind="cold", path=path.path_id, diameter_um=diameter,
                    wall_s=elapsed, cpu_s=cpu, breakdown=stats,
                    peak_rss_bytes=_rss_bytes())


def _explicit32(manifest, masks):
    stack = build_stack_z(manifest)
    specs = [spec for spec in unit_specs(manifest, masks)
             if spec.unit_kind == KIND_OPEN_PATH]
    explicit_manifest = manifest.model_copy(update={"electrode_pairs": []})
    distance_cache = DistanceCache(masks)
    result_cache = {}
    count = 0
    peak_cache_entries = 0
    preparation_count = 0
    boundary_checks = 0
    candidate_count = 0
    original_prepare = of.prepare_barrier
    original_local = of._evaluate_local_centres
    original_candidates = of.open_candidates

    def prepare(*args):
        nonlocal preparation_count
        preparation_count += 1
        return original_prepare(*args)

    def local(*args):
        nonlocal boundary_checks
        boundary_checks += int(args[3].sum())
        return original_local(*args)

    def candidates(*args):
        nonlocal candidate_count
        result = original_candidates(*args)
        candidate_count += int(result.sum())
        return result

    of.prepare_barrier = prepare
    of._evaluate_local_centres = local
    of.open_candidates = candidates
    start = time.perf_counter()
    cpu_start = time.process_time()
    try:
        for diameter in SIZES:
            for step in manifest.particle_source_steps:
                results = evaluate_units(explicit_manifest, masks, stack, step.step_id, diameter,
                                         distance_cache, result_cache, specs)
                peak_cache_entries = max(peak_cache_entries, len(result_cache))
                count += sum(1 for r in results if r.pair_id in
                             {path.path_id for path in manifest.open_paths})
            # This isolated harness needs no earlier size once every step consumed it.
            result_cache.clear()
    finally:
        of.prepare_barrier = original_prepare
        of._evaluate_local_centres = original_local
        of.open_candidates = original_candidates
    _record(kind="explicit32", wall_s=time.perf_counter() - start,
            cpu_s=time.process_time() - cpu_start, evaluations=count,
            sizes=len(SIZES), steps=len(manifest.particle_source_steps),
            peak_result_cache_entries=peak_cache_entries, preparations=preparation_count,
            boundary_checks=boundary_checks, candidates=candidate_count,
            geometry_cache_bytes=0, geometry_cache_evictions=0,
            peak_rss_bytes=_rss_bytes())


def _full32(source: Path):
    with tempfile.TemporaryDirectory(prefix="open-barrier-full32-") as scratch:
        scratch = Path(scratch)
        copied = scratch / "product.yaml"
        output = scratch / "output"
        _copy_manifest(source, copied)
        command = [str(Path(sys.executable).parent / "particle-sim"), "simulate-nominal",
                   "--product-config", str(copied), "--output", str(output)]
        start = time.perf_counter()
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        try:
            result = subprocess.run(command, check=False, capture_output=True, text=True,
                                    timeout=1800)
        except subprocess.TimeoutExpired as exc:
            after = resource.getrusage(resource.RUSAGE_CHILDREN)
            _record(kind="full32", timed_out=True, limit_s=1800, command=command,
                    wall_s=time.perf_counter() - start,
                    cpu_s=after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime,
                    peak_rss_bytes=int(after.ru_maxrss *
                                       (1024 if sys.platform.startswith("linux") else 1)),
                    stdout_tail=str(exc.stdout)[-2000:], stderr_tail=str(exc.stderr)[-2000:])
            raise SystemExit(124) from exc
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        dimensions = None
        explicit_classifications = None
        if result.returncode == 0:
            group = zarr.open_group(str(output / "hotspots.zarr"), mode="r")
            dimensions = {"nominal_failure": list(group["nominal_failure"].shape),
                          "source_steps": len(group.attrs["source_steps"]),
                          "particle_sizes_um": len(group.attrs["particle_sizes_um"]),
                          "failure_unit_ids": len(group.attrs["failure_unit_ids"])}
            ca = pd.read_parquet(output / "ca.parquet")
            explicit_rows = ca[ca["unit_kind"] == KIND_OPEN_PATH]
            explicit_classifications = {
                str(key): int(value) for key, value in
                explicit_rows["classification"].value_counts().items()
            }
        _record(kind="full32", command=command, returncode=result.returncode,
                wall_s=time.perf_counter() - start, output_files=sorted(p.name for p in
                output.iterdir()) if output.exists() else [],
                cpu_s=after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime,
                peak_rss_bytes=int(after.ru_maxrss *
                                   (1024 if sys.platform.startswith("linux") else 1)),
                dimensions=dimensions, explicit_classifications=explicit_classifications,
                stdout_tail=result.stdout[-2000:], stderr_tail=result.stderr[-2000:])
        if result.returncode:
            raise SystemExit(result.returncode)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["synthetic", "equality", "benchmark",
                                         "explicit32", "full32"])
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--sizes", nargs="*", type=float)
    args = parser.parse_args()
    if args.mode == "synthetic":
        subprocess.run([sys.executable, "-m", "pytest", "-q",
                        "tests/unit/test_open_barrier.py"], check=True)
        return
    if args.manifest is None:
        parser.error("--manifest is required for real-data modes")
    source = args.manifest.resolve()
    if args.mode == "full32":
        _full32(source)
        return
    manifest = load_manifest(source)
    masks = load_product_masks(manifest, source)
    _meta(source, masks)
    if args.mode == "equality":
        _equality(manifest, masks, COMPARE if args.sizes is None else args.sizes)
    elif args.mode == "benchmark":
        _benchmark(manifest, masks)
    elif args.mode == "explicit32":
        _explicit32(manifest, masks)


if __name__ == "__main__":
    main()
