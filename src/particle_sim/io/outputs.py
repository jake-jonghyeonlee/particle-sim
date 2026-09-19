"""Run outputs: zarr hotspot store, parquet CA tables and the run manifest."""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
import zarr

SCHEMA_VERSION = "2.0"

CA_COLUMNS: tuple[str, ...] = (
    "schema_version", "run_id", "product_id", "design_version", "source_step",
    "fab_distribution_layer", "particle_size_um", "failure_unit_id", "unit_kind", "pair_id",
    "electrode_pair_id", "open_path_id", "contact_id", "electrode_a", "electrode_b",
    "electrode_a_layer", "electrode_b_layer", "net_a", "net_b", "failure_mode",
    "classification", "ca_legacy_um2",
    "ca_nominal_um2", "ca_allocated_nominal_um2", "ca_unique_nominal_um2",
    "ca_nominal_mm2", "ca_nominal_normalized", "ca_mean_um2", "ca_p05_um2", "ca_p50_um2",
    "ca_p95_um2", "ca_allocated_mean_um2", "ca_allocated_p05_um2",
    "ca_allocated_p50_um2", "ca_allocated_p95_um2", "ca_unique_mean_um2",
    "ca_unique_p05_um2", "ca_unique_p50_um2", "ca_unique_p95_um2", "scenario_count",
    "convergence_delta", "assumption_set", "simulator_version",
)
UNION_CA_COLUMNS: tuple[str, ...] = (
    "schema_version", "run_id", "product_id", "design_version", "source_step",
    "fab_distribution_layer", "particle_size_um", "failure_mode", "union_scope", "pair_ids",
    "failure_unit_ids",
    "ca_union_legacy_um2", "ca_union_nominal_um2", "ca_union_nominal_normalized",
    "ca_pair_sum_nominal_um2", "ca_union_mean_um2", "ca_union_p05_um2",
    "ca_union_p50_um2", "ca_union_p95_um2", "ca_raw_unit_sum_mean_um2", "scenario_count",
    "convergence_delta", "assumption_set", "simulator_version",
)
def _default_unit_kinds(modes: list[str]) -> list[str]:
    """Kinds for callers that only know the failure modes (pairs and explicit open paths)."""
    return ["open_path" if mode == "lateral_open" else "electrode_pair" for mode in modes]


VALUE_SEMANTICS = "binary geometric contact risk (not failure probability)"
PROBABILITY_SEMANTICS = (
    "process-variation-conditioned geometric contact risk in [0, 1] "
    "(not an actual failure probability)"
)


def write_table(rows: list[dict], columns: tuple[str, ...], path: Path) -> None:
    """Write rows as parquet with exactly `columns`, in that order."""
    expected = set(columns)
    for index, row in enumerate(rows):
        keys = set(row)
        if keys != expected:
            raise ValueError(
                f"행 {index}의 열이 계약과 다릅니다: 누락={sorted(expected - keys)}, "
                f"추가={sorted(keys - expected)}"
            )
    frame = pd.DataFrame(rows, columns=list(columns))
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(path, index=False, engine="pyarrow")


class HotspotStore:
    """Zarr v3 store for legacy / nominal hotspots, margins and unions."""

    def __init__(
        self,
        path: Path,
        *,
        source_steps: list[str],
        particle_sizes_um: list[float],
        pair_ids: list[str],
        pair_modes: list[str],
        union_keys: list[str],
        shape_yx: tuple[int, int],
        unit_kinds: list[str] | None = None,
    ) -> None:
        if len(pair_ids) != len(pair_modes):
            raise ValueError("pair_ids와 pair_modes의 길이가 다릅니다")
        self.path = Path(path)
        self.shape_yx = shape_yx
        self._counts = (len(source_steps), len(particle_sizes_um))
        self._pair_count = len(pair_ids)
        self._union_count = len(union_keys)
        chunks = (1, 1, 1, min(512, shape_yx[0]), min(512, shape_yx[1]))
        group = zarr.open_group(self.path, mode="w-")
        pair_shape = (len(source_steps), len(particle_sizes_um), len(pair_ids), *shape_yx)
        union_shape = (len(source_steps), len(particle_sizes_um), len(union_keys), *shape_yx)
        group.create_array("legacy_failure", shape=pair_shape, dtype="uint8", chunks=chunks,
                           fill_value=0)
        group.create_array("nominal_failure", shape=pair_shape, dtype="uint8", chunks=chunks,
                           fill_value=0)
        group.create_array("pair_margin_um", shape=pair_shape, dtype="float32", chunks=chunks,
                           fill_value=float("nan"))
        group.create_array("union_failure", shape=union_shape, dtype="uint8", chunks=chunks,
                           fill_value=0)
        group.create_array("union_pair_count", shape=union_shape, dtype="uint16", chunks=chunks,
                           fill_value=0)
        group.attrs["axes"] = ["source_step", "particle_size_um", "pair_or_union", "y", "x"]
        group.attrs["source_steps"] = list(source_steps)
        group.attrs["particle_sizes_um"] = [float(s) for s in particle_sizes_um]
        group.attrs["pair_ids"] = list(pair_ids)
        group.attrs["failure_unit_ids"] = list(pair_ids)
        group.attrs["unit_kinds"] = (
            list(unit_kinds) if unit_kinds is not None else _default_unit_kinds(pair_modes)
        )
        group.attrs["pair_modes"] = list(pair_modes)
        group.attrs["union_keys"] = list(union_keys)
        group.attrs["units"] = {"length": "um", "margin": "um"}
        group.attrs["value_semantics"] = VALUE_SEMANTICS
        self._group = group

    def _check(self, array: np.ndarray) -> None:
        if tuple(array.shape) != tuple(self.shape_yx):
            raise ValueError(f"map shape {array.shape} != {self.shape_yx}")

    def _check_index(self, step_index: int, size_index: int, third: int, limit: int) -> None:
        for name, value, bound in (
            ("step_index", step_index, self._counts[0]),
            ("size_index", size_index, self._counts[1]),
            ("pair_or_union_index", third, limit),
        ):
            if not 0 <= value < bound:
                raise ValueError(f"{name}={value}가 범위 [0, {bound})를 벗어났습니다")

    def write_pair(
        self,
        step_index: int,
        size_index: int,
        pair_index: int,
        *,
        legacy: np.ndarray,
        nominal: np.ndarray,
        margin_um: np.ndarray,
    ) -> None:
        self._check_index(step_index, size_index, pair_index, self._pair_count)
        for array in (legacy, nominal, margin_um):
            self._check(np.asarray(array))
        index = (step_index, size_index, pair_index)
        self._group["legacy_failure"][index] = np.asarray(legacy, dtype=np.uint8)
        self._group["nominal_failure"][index] = np.asarray(nominal, dtype=np.uint8)
        self._group["pair_margin_um"][index] = np.asarray(margin_um, dtype=np.float32)

    def write_union(
        self,
        step_index: int,
        size_index: int,
        union_index: int,
        *,
        any_failure: np.ndarray,
        pair_count: np.ndarray,
    ) -> None:
        self._check_index(step_index, size_index, union_index, self._union_count)
        for array in (any_failure, pair_count):
            self._check(np.asarray(array))
        index = (step_index, size_index, union_index)
        self._group["union_failure"][index] = np.asarray(any_failure, dtype=np.uint8)
        self._group["union_pair_count"][index] = np.asarray(pair_count, dtype=np.uint16)


class ProbabilityStore:
    """Zarr v3 store for probability hotspots, mean margins and union probabilities."""

    def __init__(
        self,
        path: Path,
        *,
        source_steps: list[str],
        particle_sizes_um: list[float],
        unit_ids: list[str],
        unit_modes: list[str],
        union_keys: list[str],
        shape_yx: tuple[int, int],
        unit_kinds: list[str] | None = None,
    ) -> None:
        if len(unit_ids) != len(unit_modes):
            raise ValueError("unit_ids와 unit_modes의 길이가 다릅니다")
        self.path = Path(path)
        self.shape_yx = shape_yx
        self._counts = (len(source_steps), len(particle_sizes_um))
        self._unit_count = len(unit_ids)
        self._union_count = len(union_keys)
        chunks = (1, 1, 1, min(512, shape_yx[0]), min(512, shape_yx[1]))
        group = zarr.open_group(self.path, mode="w-")
        unit_shape = (len(source_steps), len(particle_sizes_um), len(unit_ids), *shape_yx)
        union_shape = (len(source_steps), len(particle_sizes_um), len(union_keys), *shape_yx)
        group.create_array("probability", shape=unit_shape, dtype="float32", chunks=chunks,
                           fill_value=float("nan"))
        group.create_array("contact_margin_mean_um", shape=unit_shape, dtype="float32",
                           chunks=chunks, fill_value=float("nan"))
        group.create_array("union_probability", shape=union_shape, dtype="float32",
                           chunks=chunks, fill_value=float("nan"))
        group.attrs["axes"] = ["source_step", "particle_size_um", "unit_or_union", "y", "x"]
        group.attrs["source_steps"] = list(source_steps)
        group.attrs["particle_sizes_um"] = [float(s) for s in particle_sizes_um]
        group.attrs["unit_ids"] = list(unit_ids)
        group.attrs["failure_unit_ids"] = list(unit_ids)
        group.attrs["unit_kinds"] = (
            list(unit_kinds) if unit_kinds is not None else _default_unit_kinds(unit_modes)
        )
        group.attrs["unit_modes"] = list(unit_modes)
        group.attrs["union_keys"] = list(union_keys)
        group.attrs["units"] = {"length": "um", "margin": "um"}
        group.attrs["value_semantics"] = PROBABILITY_SEMANTICS
        self._group = group

    def _check(self, array: np.ndarray, step_index: int, size_index: int, third: int,
               limit: int) -> None:
        if tuple(array.shape) != tuple(self.shape_yx):
            raise ValueError(f"map shape {array.shape} != {self.shape_yx}")
        for name, value, bound in (
            ("step_index", step_index, self._counts[0]),
            ("size_index", size_index, self._counts[1]),
            ("unit_or_union_index", third, limit),
        ):
            if not 0 <= value < bound:
                raise ValueError(f"{name}={value}가 범위 [0, {bound})를 벗어났습니다")

    def write_unit(self, step_index: int, size_index: int, unit_index: int, *,
                   probability: np.ndarray, margin_mean_um: np.ndarray) -> None:
        for array in (probability, margin_mean_um):
            self._check(np.asarray(array), step_index, size_index, unit_index,
                        self._unit_count)
        index = (step_index, size_index, unit_index)
        self._group["probability"][index] = np.asarray(probability, dtype=np.float32)
        self._group["contact_margin_mean_um"][index] = np.asarray(
            margin_mean_um, dtype=np.float32
        )

    def write_union(self, step_index: int, size_index: int, union_index: int, *,
                    probability: np.ndarray) -> None:
        self._check(np.asarray(probability), step_index, size_index, union_index,
                    self._union_count)
        self._group["union_probability"][(step_index, size_index, union_index)] = np.asarray(
            probability, dtype=np.float32
        )


REQUIRED_RUN_MANIFEST_KEYS: tuple[str, ...] = (
    "run_id", "product_id", "design_version", "config_hash", "assumption_set",
    "simulator_version", "legacy_kernel", "post_step_policy", "warnings",
)


def check_run_manifest(manifest: dict) -> list[str]:
    """Return the required keys missing from a run manifest."""
    return [key for key in REQUIRED_RUN_MANIFEST_KEYS if key not in manifest]


def open_hotspot_store(path: Path) -> zarr.Group:
    """Open an existing hotspot store read-only."""
    return zarr.open_group(Path(path), mode="r")


def make_run_id(product_id: str, now: datetime) -> str:
    """UTC timestamped run id."""
    return f"{now.strftime('%Y%m%dT%H%M%SZ')}_{product_id}"


def build_run_manifest(
    *,
    run_id: str,
    product_id: str,
    design_version: str,
    config_hash: str,
    input_hashes: dict[str, str],
    particle_sizes_um: list[float],
    source_steps: list[str],
    assumption_set: str,
    warnings: list[str],
    legacy_kernel: dict,
    post_step_policy: dict,
    created_at: datetime,
    git_commit: str | None = None,
) -> dict:
    """Reproducibility record stored next to every run."""
    from particle_sim import __version__

    detected_commit, git_dirty = _git_identity()
    return {
        "schema_version": SCHEMA_VERSION,
        "run_id": run_id,
        "product_id": product_id,
        "design_version": design_version,
        "git_commit": git_commit or detected_commit,
        "git_dirty": git_dirty,
        "build_id": f"particle-sim-{__version__}",
        "config_hash": config_hash,
        "input_hashes": dict(input_hashes),
        "particle_sizes_um": [float(s) for s in particle_sizes_um],
        "source_steps": list(source_steps),
        "assumption_set": assumption_set,
        "legacy_kernel": legacy_kernel,
        "post_step_policy": post_step_policy,
        "warnings": list(warnings),
        "random_seed": None,
        "sobol_exponent_m": None,
        "simulator_version": __version__,
        "created_at": created_at.isoformat(),
        "units": {"length": "um", "area": "um2"},
        "value_semantics": VALUE_SEMANTICS,
        "artifacts": [],
    }


def _git_identity() -> tuple[str | None, bool | None]:
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "status", "--porcelain"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError):
        return None, None
    return commit or None, bool(status.strip())


def artifact_digest(path: Path) -> tuple[str, int]:
    """파일 또는 디렉터리 트리의 결정적 sha256과 총 바이트 수."""
    path = Path(path)
    if path.is_file():
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        return digest, path.stat().st_size
    hasher = hashlib.sha256()
    total = 0
    for child in sorted(item for item in path.rglob("*") if item.is_file()):
        relative = child.relative_to(path).as_posix().encode("utf-8")
        payload = child.read_bytes()
        hasher.update(len(relative).to_bytes(8, "big"))
        hasher.update(relative)
        hasher.update(len(payload).to_bytes(8, "big"))
        hasher.update(payload)
        total += len(payload)
    return hasher.hexdigest(), total


def artifact_records(paths: list[Path], *, base: Path) -> list[dict]:
    records = []
    for path in paths:
        digest, size = artifact_digest(path)
        records.append(
            {
                "path": path.relative_to(base).as_posix(),
                "sha256": digest,
                "size_bytes": size,
                "kind": "directory" if path.is_dir() else "file",
            }
        )
    return records


def write_json(data: dict, path: Path) -> None:
    """Deterministic JSON dump."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(
        json.dumps(data, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8"
    )


def write_float_tiff(array: np.ndarray, path: Path) -> None:
    """Exchange format for a single float32 map."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(str(path), np.asarray(array, dtype=np.float32))
