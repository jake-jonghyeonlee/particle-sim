"""Serial, resumable auto-open oracle and candidate comparison.

Example::

    python scripts/benchmark_auto_open.py prepare --oracle-dir /tmp/ao-oracle.XXXXXX
    python scripts/benchmark_auto_open.py oracle --oracle-dir /tmp/ao-oracle.XXXXXX
    python scripts/benchmark_auto_open.py compare --oracle-dir /tmp/ao-oracle.XXXXXX

The oracle is an extracted git archive, never a checkout or the current package.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import resource
import shutil
import subprocess
import sys
import tarfile
import time
from pathlib import Path

BASE = "8393f83"
DEFAULT_MANIFEST = Path("/tmp/particle-mem-rca.0kdEy2/product-D3-five.yaml")
JUDGED = ("E00022_METAL1", "E00023_METAL1", "E00084_METAL1", "E00085_METAL1")
SIZES = (1, 2, 5, 10, 20, 30, 50)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _write_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def _disk_guard(path: Path) -> None:
    if shutil.disk_usage(path).free < 15_000_000_000:
        raise RuntimeError("free disk below 15 GB; stopping before next case")


def _source_dir(root: Path, mode: str, candidate: Path) -> Path:
    return (root / "oracle-src" / "src" if mode == "oracle" else candidate / "src").resolve()


def _import_source(source: Path):
    if not (source / "particle_sim" / "failure" / "auto_open.py").is_file():
        raise FileNotFoundError(source)
    sys.path.insert(0, str(source))
    import numpy as np
    import scipy

    from particle_sim.config import load_manifest
    from particle_sim.failure import auto_open
    from particle_sim.io.masks import load_product_masks

    actual = Path(auto_open.__file__).resolve()
    expected = (source / "particle_sim" / "failure" / "auto_open.py").resolve()
    if actual != expected:
        raise RuntimeError(f"wrong source import: {actual} != {expected}")
    return np, scipy, load_manifest, load_product_masks, auto_open, actual


def prepare(root: Path) -> None:
    root.mkdir(parents=True, exist_ok=True)
    _disk_guard(root)
    source = root / "oracle-src"
    if not source.exists():
        source.mkdir()
        archive = subprocess.run(
            ["git", "archive", "--format=tar", BASE, "src"],
            check=True, capture_output=True,
        ).stdout
        with tarfile.open(fileobj=__import__("io").BytesIO(archive)) as tar:
            tar.extractall(source, filter="data")
    expected = source / "src" / "particle_sim" / "failure" / "auto_open.py"
    if not expected.is_file():
        raise RuntimeError(f"incomplete archive: {expected}")
    print(json.dumps({"archive": str(source), "auto_open_sha256": sha256(expected.read_bytes())}))


def _metadata(result, diameter: int, elapsed: float, np, rss_before: int) -> dict:
    array = result.failure
    if array.dtype != np.dtype("bool"):
        raise AssertionError(f"non-bool failure dtype: {array.dtype}")
    return {
        "electrode": result.electrode,
        "diameter_um": diameter,
        "classification": result.classification,
        "terminals": list(result.terminals),
        "candidates": int(result.candidates),
        "dtype": str(array.dtype),
        "shape": list(array.shape),
        "failure_sha256": sha256(array.tobytes(order="C")),
        "opens": int(np.count_nonzero(array)),
        "wall_s": elapsed,
        "peak_rss_bytes": _rss_bytes(),
        "rss_before_bytes": rss_before,
    }


def _rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value * (1024 if sys.platform.startswith("linux") else 1))


def _cases(masks, sizes: list[int], all_electrodes: bool) -> list[tuple[str, int]]:
    electrodes = (
        sorted(name for name in masks.electrodes if name.endswith("_METAL1"))
        if all_electrodes else JUDGED
    )
    missing = set(electrodes) - set(masks.electrodes)
    if missing:
        raise ValueError(f"missing electrodes: {sorted(missing)}")
    return [(name, diameter) for diameter in sizes for name in electrodes]


def _verify_existing(path: Path, meta: dict, np) -> None:
    with np.load(path, allow_pickle=False) as stored:
        array = stored["failure"]
    if (str(array.dtype) != meta["dtype"] or list(array.shape) != meta["shape"]
            or sha256(array.tobytes(order="C")) != meta["failure_sha256"]):
        raise RuntimeError(f"stored oracle corruption: {path}")


def run(args: argparse.Namespace) -> None:
    root = args.oracle_dir.resolve()
    if args.mode == "oracle":
        prepare(root)
    source = _source_dir(root, args.mode, args.candidate.resolve())
    np, scipy, load_manifest, load_product_masks, auto_open, actual = _import_source(source)
    manifest_path = args.manifest.resolve()
    manifest = load_manifest(manifest_path)
    masks = load_product_masks(manifest, manifest_path)
    inputs = {str(k): str(v) for k, v in masks.input_hashes.items()}
    provenance = {
        "base_commit": subprocess.check_output(["git", "rev-parse", BASE], text=True).strip(),
        "candidate_head": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "mode": args.mode,
        "source": str(source),
        "asserted_import_path": str(actual),
        "source_auto_open_sha256": sha256(actual.read_bytes()),
        "manifest": str(manifest_path),
        "manifest_sha256": sha256(manifest_path.read_bytes()),
        "input_hashes": inputs,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
    }
    _write_json(root / f"{args.mode}-provenance.json", provenance)
    print(json.dumps({"kind": "provenance", "mode": args.mode,
                      "asserted_import_path": str(actual),
                      "manifest_sha256": provenance["manifest_sha256"],
                      "input_count": len(inputs)}), flush=True)
    oracle_cases = root / "cases"
    oracle_cases.mkdir(exist_ok=True)
    mismatches = 0
    for electrode, diameter in _cases(masks, args.sizes, args.all_electrodes):
        _disk_guard(root)
        stem = f"{electrode}-d{diameter}"
        map_path = oracle_cases / f"{stem}.npz"
        meta_path = oracle_cases / f"{stem}.json"
        if args.mode == "oracle" and map_path.exists() and meta_path.exists():
            meta = json.loads(meta_path.read_text())
            _verify_existing(map_path, meta, np)
            if meta["manifest_sha256"] != provenance["manifest_sha256"] or meta["input_hashes"] != inputs:
                raise RuntimeError(f"input drift in {stem}")
            print(json.dumps({"case": stem, "status": "resumed", "opens": meta["opens"]}), flush=True)
            continue
        rss_before = _rss_bytes()
        start = time.perf_counter()
        result = auto_open.simulate_auto_open(
            masks, electrode, float(diameter), float(manifest.contact_margin_um)
        )
        meta = _metadata(result, diameter, time.perf_counter() - start, np, rss_before)
        meta.update(manifest_sha256=provenance["manifest_sha256"], input_hashes=inputs)
        if args.mode == "oracle":
            temporary = map_path.with_suffix(".tmp.npz")
            np.savez_compressed(temporary, failure=result.failure)
            temporary.replace(map_path)
            _write_json(meta_path, meta)
            print(json.dumps({"case": stem, "status": "stored", "opens": meta["opens"],
                              "candidates": meta["candidates"], "wall_s": meta["wall_s"],
                              "peak_rss_bytes": meta["peak_rss_bytes"]}), flush=True)
        else:
            if not map_path.exists() or not meta_path.exists():
                raise FileNotFoundError(f"oracle missing: {stem}")
            expected = json.loads(meta_path.read_text())
            _verify_existing(map_path, expected, np)
            if expected["manifest_sha256"] != provenance["manifest_sha256"] or expected["input_hashes"] != inputs:
                raise RuntimeError(f"input drift in {stem}")
            fields = ("electrode", "classification", "terminals", "candidates", "dtype", "shape", "failure_sha256")
            differences = {field: [expected[field], meta[field]] for field in fields if expected[field] != meta[field]}
            with np.load(map_path, allow_pickle=False) as stored:
                reference = stored["failure"]
            if reference.shape == result.failure.shape:
                coordinates = np.argwhere(reference != result.failure)
                count = len(coordinates)
                examples = coordinates[:10].tolist()
            else:
                count, examples = None, []
            mismatches += bool(differences)
            print(json.dumps({"case": stem, "status": "equal" if not differences else "mismatch",
                              "field_differences": differences, "mismatch_pixels": count,
                              "example_coordinates": examples, "candidate_wall_s": meta["wall_s"],
                              "candidate_peak_rss_bytes": meta["peak_rss_bytes"]}), flush=True)
    if mismatches:
        raise SystemExit(f"{mismatches} mismatched cases")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("prepare", "oracle", "compare"))
    parser.add_argument("--oracle-dir", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, default=Path.cwd())
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--sizes", type=int, nargs="+", default=list(SIZES))
    parser.add_argument("--all-electrodes", action="store_true")
    args = parser.parse_args()
    if args.mode == "prepare":
        prepare(args.oracle_dir.resolve())
    else:
        run(args)


if __name__ == "__main__":
    main()
