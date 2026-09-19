"""YAML manifest 로더와 경로/해시 도우미."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import yaml
from pydantic import ValidationError

from particle_sim.errors import InputContractError
from particle_sim.schema import ProductManifest, ScenarioConfig


def load_manifest(path: Path) -> ProductManifest:
    """YAML 파일에서 ProductManifest를 로드한다."""
    try:
        with path.open("r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
    except OSError as e:
        raise InputContractError(f"manifest 파일을 읽을 수 없습니다: {e}", file=str(path)) from e
    except yaml.YAMLError as e:
        raise InputContractError(f"YAML 파싱 실패: {e}", file=str(path)) from e
    if not isinstance(data, dict):
        raise InputContractError("manifest는 mapping이어야 합니다", file=str(path))
    try:
        manifest = ProductManifest.model_validate(data)
        if manifest.schema_version == "2.0":
            from particle_sim.validation.electrode_pairs import unresolved_electrode_pairs

            unresolved = unresolved_electrode_pairs(manifest)
            if unresolved:
                pairs = [f"{row.electrode_a}--{row.electrode_b}" for row in unresolved]
                raise InputContractError(
                    f"계산 또는 제외로 정하지 않은 전극 조합이 있습니다: {pairs}",
                    product_id=manifest.product_id,
                    file=str(path),
                    field="electrode_pairs",
                )
        return manifest
    except InputContractError:
        raise
    except ValidationError as e:
        product_id: str | None = None
        if isinstance(data.get("product_id"), str):
            product_id = data["product_id"]
        first = e.errors()[0]
        loc = ".".join(str(x) for x in first["loc"])
        raise InputContractError(
            f"manifest 검증 실패: {first['msg']}",
            product_id=product_id,
            file=str(path),
            field=loc,
        ) from e


def load_scenario_config(path: Path) -> ScenarioConfig:
    """YAML 파일에서 산포 시나리오 설정을 로드한다."""
    try:
        with path.open("r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
    except OSError as e:
        raise InputContractError(
            f"시나리오 파일을 읽을 수 없습니다: {e}", file=str(path)
        ) from e
    except yaml.YAMLError as e:
        raise InputContractError(f"YAML 파싱 실패: {e}", file=str(path)) from e
    if not isinstance(data, dict):
        raise InputContractError("시나리오 설정은 mapping이어야 합니다", file=str(path))
    try:
        return ScenarioConfig.model_validate(data)
    except ValidationError as e:
        first = e.errors()[0]
        raise InputContractError(
            f"시나리오 설정 검증 실패: {first['msg']}",
            file=str(path),
            field=".".join(str(x) for x in first["loc"]),
        ) from e


def validate_scenario_for_manifest(
    scenario: ScenarioConfig, manifest: ProductManifest, *, scenario_path: Path
) -> None:
    """산포 항목이 실제 manifest layer에 적용되는지 fail-closed로 확인한다."""
    stack_names = {layer.layer_name for layer in manifest.stack}
    # contact·via mask는 hole_layer 이름으로 overlay·CD를 받는다
    mask_names = {layer.layer_name for layer in manifest.layer_masks} | {
        contact.hole_layer for contact in manifest.contacts
    }
    checks = (
        ("overlay_um", set(scenario.overlay_um), mask_names),
        ("cd_bias_um", set(scenario.cd_bias_um), mask_names),
        ("thickness_delta_um", set(scenario.thickness_delta_um), stack_names),
        ("planarization_factor", set(scenario.planarization_factor), stack_names),
    )
    for field, provided, allowed in checks:
        unknown = sorted(provided - allowed)
        if unknown:
            raise InputContractError(
                f"{field}에 manifest에 없는 layer가 있습니다: {unknown}",
                product_id=manifest.product_id,
                file=str(scenario_path),
                field=field,
            )
    reference = scenario.overlay_reference_layer
    if reference is not None and reference not in mask_names:
        raise InputContractError(
            f"overlay 기준 layer '{reference}'가 layer_masks에 없습니다",
            product_id=manifest.product_id,
            file=str(scenario_path),
            field="overlay_reference_layer",
        )
    if reference is not None and reference in scenario.overlay_um:
        raise InputContractError(
            "overlay 기준 layer에는 이동 분포를 지정할 수 없습니다",
            product_id=manifest.product_id,
            file=str(scenario_path),
            field="overlay_reference_layer",
        )


def resolve_input_path(
    manifest_path: Path, relative: str, *, product_id: str, field: str
) -> Path:
    """manifest 디렉토리 기준 상대 경로를 해석한다.

    절대 경로 또는 manifest 디렉토리 밖으로 벗어나는 경로는 거부한다.
    """
    base = manifest_path.resolve().parent
    p = Path(relative)
    if p.is_absolute():
        raise InputContractError(
            f"경로 '{relative}'은(는) 절대 경로일 수 없습니다",
            product_id=product_id,
            file=str(manifest_path),
            field=field,
        )
    resolved = (base / p).resolve()
    if not resolved.is_relative_to(base):
        raise InputContractError(
            f"경로 '{relative}'은(는) manifest 디렉토리 밖으로 벗어날 수 없습니다",
            product_id=product_id,
            file=str(manifest_path),
            field=field,
        )
    return resolved


def canonical_hash(obj: object) -> str:
    """객체의 정준 JSON에 대한 sha256 hex를 반환한다."""
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    ).hexdigest()


def file_sha256(path: Path) -> str:
    """파일의 sha256 hex를 반환한다."""
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()
