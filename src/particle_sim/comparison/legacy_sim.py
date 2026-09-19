"""외부 기존 sim hotspot을 명시적 좌표 계약으로 읽어 픽셀 단위 비교한다."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import tifffile
import yaml
import zarr

from particle_sim.config import file_sha256
from particle_sim.errors import InputContractError
from particle_sim.io.outputs import artifact_records, write_json

COMPARISON_COLUMNS = (
    "source_step", "particle_size_um", "failure_unit_id", "failure_mode", "status",
    "reference_ca_um2", "candidate_ca_um2", "ca_delta_um2", "ca_relative_delta",
    "intersection_pixels", "reference_only_pixels", "candidate_only_pixels", "xor_pixels",
    "difference_classification",
)


def _load_config(path: Path) -> dict:
    try:
        if path.suffix.lower() == ".json":
            data = json.loads(path.read_text(encoding="utf-8"))
        else:
            data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, yaml.YAMLError) as exc:
        raise InputContractError(
            f"기존 sim 비교 설정을 읽을 수 없습니다: {exc}", file=str(path)
        ) from exc
    if not isinstance(data, dict):
        raise InputContractError("기존 sim 비교 설정은 mapping이어야 합니다", file=str(path))
    return data


def _reference_array(config_path: Path, relative: str) -> np.ndarray:
    base = config_path.resolve().parent
    path = (base / relative).resolve()
    if not path.is_relative_to(base):
        raise InputContractError(
            "기존 sim 배열 경로가 설정 폴더 밖으로 나갑니다", file=str(config_path), field="path"
        )
    try:
        if path.suffix.lower() == ".npy":
            array = np.load(path, allow_pickle=False)
        elif path.suffix.lower() in (".tif", ".tiff"):
            array = tifffile.imread(path)
        else:
            raise InputContractError(
                f"지원하지 않는 기존 sim 배열 확장자입니다: {path.suffix}",
                file=str(path),
                field="path",
            )
    except (OSError, ValueError) as exc:
        raise InputContractError(
            f"기존 sim 배열을 읽을 수 없습니다: {exc}", file=str(path), field="path"
        ) from exc
    if array.ndim != 2:
        raise InputContractError(
            f"기존 sim hotspot은 2-D여야 합니다: shape={array.shape}",
            file=str(path),
            field="path",
        )
    values = set(np.unique(array).tolist())
    if not values <= {0, 1, 255}:
        raise InputContractError(
            f"기존 sim hotspot은 이진 배열이어야 합니다: values={sorted(values)[:10]}",
            file=str(path),
            field="path",
        )
    return array.astype(bool)


def _not_applicable_keys(run_dir: Path) -> set[tuple[str, float, str]]:
    """Conditions the candidate run did not evaluate because the source step does not apply."""
    path = run_dir / "ca.parquet"
    if not path.is_file():
        return set()
    import pandas as pd

    frame = pd.read_parquet(path)
    unit_column = "failure_unit_id" if "failure_unit_id" in frame.columns else "pair_id"
    if "classification" not in frame.columns or unit_column not in frame.columns:
        return set()
    rows = frame[frame["classification"] == "not_applicable_source_step"]
    return {
        (str(row.source_step), float(row.particle_size_um), str(getattr(row, unit_column)))
        for row in rows.itertuples()
    }


def compare_legacy_sim(reference_config: Path, run_dir: Path, output_dir: Path) -> Path:
    """기존 sim과 nominal hotspot을 공통 항목에서 비교하고 신규 조합을 분리한다."""
    reference_config = Path(reference_config)
    run_dir = Path(run_dir)
    output_dir = Path(output_dir)
    if output_dir.exists() and (not output_dir.is_dir() or any(output_dir.iterdir())):
        raise InputContractError(
            f"비교 출력 경로가 비어 있지 않습니다: {output_dir}", field="output"
        )
    config = _load_config(reference_config)
    required = {
        "schema_version", "axes", "origin", "pixel_size_x_um", "pixel_size_y_um",
        "particle_size_definition", "records",
    }
    missing = sorted(required - set(config))
    if missing:
        raise InputContractError(
            f"기존 sim 비교 설정 필드가 누락됐습니다: {missing}",
            file=str(reference_config),
            field=missing[0],
        )
    if str(config["schema_version"]) != "1.0":
        raise InputContractError(
            "기존 sim 비교 설정 schema_version은 1.0이어야 합니다",
            file=str(reference_config),
            field="schema_version",
        )
    if config["axes"] != "yx" or config["origin"] != "upper_left":
        raise InputContractError(
            "기존 sim 배열은 axes=yx, origin=upper_left로 명시해야 합니다",
            file=str(reference_config),
            field="axes",
        )
    if config["particle_size_definition"] != "diameter":
        raise InputContractError(
            "particle_size_definition은 diameter여야 합니다",
            file=str(reference_config),
            field="particle_size_definition",
        )
    try:
        run_manifest = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
        store = zarr.open_group(run_dir / "hotspots.zarr", mode="r")
    except (OSError, ValueError) as exc:
        raise InputContractError(f"비교할 실행 결과를 읽을 수 없습니다: {exc}", field="run") from exc
    pixel_x = float(run_manifest.get("pixel_size_x_um", config["pixel_size_x_um"]))
    pixel_y = float(run_manifest.get("pixel_size_y_um", config["pixel_size_y_um"]))
    # Older manifests did not store each axis, so derive it from CA normalization only when config
    # and stored array shape agree. New manifests always provide the explicit values.
    if "pixel_size_x_um" in run_manifest:
        if not np.isclose(pixel_x, float(config["pixel_size_x_um"])) or not np.isclose(
            pixel_y, float(config["pixel_size_y_um"])
        ):
            raise InputContractError("기존 sim과 현재 실행의 pixel 크기가 다릅니다", field="pixel_size")
    steps = list(store.attrs["source_steps"])
    sizes = [float(item) for item in store.attrs["particle_sizes_um"]]
    units = list(store.attrs.get("failure_unit_ids", store.attrs.get("pair_ids", [])))
    modes = list(store.attrs.get("unit_modes", store.attrs.get("pair_modes", [])))
    kinds = list(store.attrs.get("unit_kinds", ["electrode_pair"] * len(units)))
    not_applicable = _not_applicable_keys(run_dir)
    candidate = store["nominal_failure"]
    index = {
        (step, size, unit): (step_i, size_i, unit_i)
        for step_i, step in enumerate(steps)
        for size_i, size in enumerate(sizes)
        for unit_i, unit in enumerate(units)
    }
    records = config["records"]
    if not isinstance(records, list):
        raise InputContractError("records는 list여야 합니다", field="records")
    reference_by_key: dict[tuple[str, float, str], tuple[np.ndarray, str]] = {}
    reference_hashes = {"reference_config": file_sha256(reference_config)}
    for position, record in enumerate(records):
        if not isinstance(record, dict):
            raise InputContractError("record는 mapping이어야 합니다", field=f"records[{position}]")
        try:
            key = (
                str(record["source_step"]),
                float(record["particle_size_um"]),
                str(record["failure_unit_id"]),
            )
            relative = str(record["path"])
        except (KeyError, TypeError, ValueError) as exc:
            raise InputContractError(
                f"기존 sim record가 잘못됐습니다: {exc}", field=f"records[{position}]"
            ) from exc
        if key in reference_by_key:
            raise InputContractError(f"기존 sim record가 중복됩니다: {key}", field="records")
        array = _reference_array(reference_config, relative)
        reference_by_key[key] = (array, relative)
        reference_hashes[relative] = file_sha256(reference_config.parent / relative)

    unknown = sorted(set(reference_by_key) - set(index))
    if unknown:
        raise InputContractError(
            f"현재 실행에 없는 기존 sim 항목이 있습니다: {unknown}", field="records"
        )
    known_reference_units = {key[2] for key in reference_by_key}
    missing_common = sorted(
        key
        for key in index
        if key[2] in known_reference_units
        and key not in reference_by_key
        and key not in not_applicable
    )
    if missing_common:
        raise InputContractError(
            f"기존 sim 공통 전극 조합의 공정 또는 이물 크기 결과가 누락됐습니다: {missing_common}",
            field="records",
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    common_keys = sorted(reference_by_key)
    shape = tuple(candidate.shape[-2:])
    diff_store = zarr.open_group(output_dir / "comparison_maps.zarr", mode="w-")
    chunks = (1, min(512, shape[0]), min(512, shape[1]))
    for name in ("intersection", "reference_only", "candidate_only", "xor"):
        diff_store.create_array(name, shape=(len(common_keys), *shape), dtype="uint8", chunks=chunks)
    rows: list[dict] = []
    pixel_area = float(config["pixel_size_x_um"]) * float(config["pixel_size_y_um"])
    for row_index, key in enumerate(common_keys):
        reference, _relative = reference_by_key[key]
        step_i, size_i, unit_i = index[key]
        current = np.asarray(candidate[step_i, size_i, unit_i], dtype=bool)
        if reference.shape != current.shape:
            raise InputContractError(
                f"hotspot shape가 다릅니다: reference={reference.shape}, current={current.shape}",
                field="shape",
            )
        intersection = reference & current
        reference_only = reference & ~current
        candidate_only = current & ~reference
        xor = reference ^ current
        for name, array in (
            ("intersection", intersection),
            ("reference_only", reference_only),
            ("candidate_only", candidate_only),
            ("xor", xor),
        ):
            diff_store[name][row_index] = array.astype(np.uint8)
        reference_ca = float(reference.sum()) * pixel_area
        candidate_ca = float(current.sum()) * pixel_area
        relative_delta = (
            float("nan") if reference_ca == 0 else (candidate_ca - reference_ca) / reference_ca
        )
        rows.append(
            {
                "source_step": key[0],
                "particle_size_um": key[1],
                "failure_unit_id": key[2],
                "failure_mode": modes[unit_i],
                "status": "common",
                "reference_ca_um2": reference_ca,
                "candidate_ca_um2": candidate_ca,
                "ca_delta_um2": candidate_ca - reference_ca,
                "ca_relative_delta": relative_delta,
                "intersection_pixels": int(intersection.sum()),
                "reference_only_pixels": int(reference_only.sum()),
                "candidate_only_pixels": int(candidate_only.sum()),
                "xor_pixels": int(xor.sum()),
                "difference_classification": "equal" if not xor.any() else "unclassified",
            }
        )
    for key in sorted(
        item for item in set(index) - set(reference_by_key) if item[2] not in known_reference_units
    ):
        step_i, size_i, unit_i = index[key]
        current = np.asarray(candidate[step_i, size_i, unit_i], dtype=bool)
        candidate_ca = float(current.sum()) * pixel_area
        if key in not_applicable:
            status = "not_applicable"
        elif kinds[unit_i] == "electrode_pair":
            status = "new_pair"
        else:
            status = "new_unit"  # open path/electrode or CNO absent from the existing sim
        rows.append(
            {
                "source_step": key[0],
                "particle_size_um": key[1],
                "failure_unit_id": key[2],
                "failure_mode": modes[unit_i],
                "status": status,
                "reference_ca_um2": float("nan"),
                "candidate_ca_um2": candidate_ca,
                "ca_delta_um2": float("nan"),
                "ca_relative_delta": float("nan"),
                "intersection_pixels": 0,
                "reference_only_pixels": 0,
                "candidate_only_pixels": int(current.sum()),
                "xor_pixels": int(current.sum()),
                "difference_classification": status,
            }
        )
    diff_store.attrs["record_keys"] = [list(key) for key in common_keys]
    import pandas as pd

    frame = pd.DataFrame(rows, columns=list(COMPARISON_COLUMNS))
    frame.to_parquet(output_dir / "comparison.parquet", index=False)
    frame.to_csv(output_dir / "comparison.csv", index=False)
    manifest = {
        "schema_version": "1.0",
        "candidate_run_id": run_manifest["run_id"],
        "candidate_artifacts": run_manifest.get("artifacts", []),
        "reference_hashes": reference_hashes,
        "common_records": len(common_keys),
        "new_pairs": int((frame["status"] == "new_pair").sum()),
        "new_units": int((frame["status"] == "new_unit").sum()),
        "not_applicable": int((frame["status"] == "not_applicable").sum()),
        "unclassified_differences": int(
            (frame["difference_classification"] == "unclassified").sum()
        ),
    }
    manifest["artifacts"] = artifact_records(
        [
            output_dir / "comparison_maps.zarr",
            output_dir / "comparison.parquet",
            output_dir / "comparison.csv",
        ],
        base=output_dir,
    )
    write_json(manifest, output_dir / "comparison_manifest.json")
    return output_dir
