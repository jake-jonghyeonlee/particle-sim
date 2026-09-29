"""편집 지정서로 제품 입력의 국소 편집 변경본 패키지를 만든다.

지정서(JSON/YAML)는 simcompare의 edit_candidates.json이나 설계자가 직접 쓴 파일이다.
변경본마다 출력 폴더 아래 <candidate_id>/에 다음을 쓴다.

- manifest.yaml: 원본과 같고 design_version만 바뀐 manifest
- 원본 manifest가 참조한 모든 마스크(편집한 전극·layer 마스크는 편집본, 나머지는 복사본)
- variant_record.json: 원안, 편집 내용, 전후 측정값, 점검 결과
- edit_preview.png: 원안 전극(회색), 제거된 픽셀(빨강), 상대 전극(파랑), ROI(노랑 테두리)

점검은 표시만 하고 판정하지 않는다. 다른 불량 증가 허용과 설계 규칙은 설계자가 판단한다.
open 경로 단자 연결이 끊기거나 전극이 비면 재시뮬레이션 입력이 되지 않으므로 오류로 멈춘다.
"""

from __future__ import annotations

import json
import re
import shutil
from dataclasses import replace
from pathlib import Path
from typing import Literal

import numpy as np
import yaml
from PIL import Image
from pydantic import BaseModel, ConfigDict, Field, ValidationError
from scipy.ndimage import label

from particle_sim.config import file_sha256, load_manifest, resolve_input_path
from particle_sim.errors import InputContractError
from particle_sim.failure.open_failure import check_open_baselines
from particle_sim.io.masks import MaskSet, load_product_masks
from particle_sim.schema import ProductManifest
from particle_sim.units import unit_specs
from particle_sim.variant.edits import apply_edit, roi_mask

_SAFE_ID = r"^[A-Za-z0-9_.-]+$"
_EIGHT = np.ones((3, 3), dtype=bool)


class EditOp(BaseModel):
    """국소 편집 하나. 전극은 electrode(ID) 또는 sample_id의 slot으로 지정한다."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    op: Literal["remove_region", "trim_near", "remove_overlap"]
    electrode: str | None = None
    slot: Literal["a", "b"] | None = None
    target: str | None = None
    target_slot: Literal["a", "b"] | None = None
    roi_um: tuple[float, float, float, float]
    amount_um: float = Field(default=0.0, ge=0)


class VariantSpec(BaseModel):
    """변경본 하나. info는 simcompare가 붙인 참고 정보이며 그대로 기록만 한다."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    candidate_id: str = Field(pattern=_SAFE_ID)
    sample_id: str | None = None
    design_version: str | None = Field(default=None, pattern=_SAFE_ID)
    change_description: str = Field(min_length=1)
    edits: list[EditOp] = Field(min_length=1)
    basis: dict | None = None
    info: dict | None = None


class EditSpecFile(BaseModel):
    """편집 지정서. 길이 단위는 um, 좌표는 좌상단 원점·x 오른쪽·y 아래다."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    schema_version: Literal["1.0"]
    length_unit: Literal["um"]
    note: str | None = None
    candidates: list[VariantSpec] = Field(min_length=1)


def load_edit_spec(path: Path) -> EditSpecFile:
    """JSON 또는 YAML 편집 지정서를 읽는다."""
    try:
        text = Path(path).read_text(encoding="utf-8")
        data = json.loads(text) if Path(path).suffix.lower() == ".json" else yaml.safe_load(text)
    except (OSError, json.JSONDecodeError, yaml.YAMLError) as exc:
        raise InputContractError(f"편집 지정서를 읽을 수 없습니다: {exc}", file=str(path)) from exc
    try:
        return EditSpecFile.model_validate(data)
    except ValidationError as exc:
        first = exc.errors()[0]
        raise InputContractError(
            f"편집 지정서 검증 실패: {first['msg']}", file=str(path),
            field=".".join(str(x) for x in first["loc"]),
        ) from exc


def _sample_units(manifest: ProductManifest, masks: MaskSet) -> dict[str, tuple[str, str]]:
    """export-simcompare와 같은 규칙의 sample_id -> (electrode_a, electrode_b)."""
    out = {}
    for step in manifest.particle_source_steps:
        for spec in unit_specs(manifest, masks):
            if spec.unit_kind == "contact":
                continue
            raw = "-".join((manifest.product_id, manifest.design_version, step.step_id, spec.unit_id))
            out[re.sub(r"[^A-Za-z0-9_.-]", "_", raw)] = (spec.electrode_a, spec.electrode_b)
    return out


def _resolve(edit: EditOp, spec: VariantSpec, units: dict, electrodes: dict) -> tuple[str, str | None]:
    """편집 대상 전극 ID와 상대 전극 ID."""
    def pick(explicit, slot, name):
        if explicit is not None:
            if slot is not None:
                raise InputContractError(f"{spec.candidate_id}: {name}와 {name}_slot을 함께 쓸 수 없습니다")
            return explicit
        if slot is None:
            return None
        if spec.sample_id is None or spec.sample_id not in units:
            raise InputContractError(
                f"{spec.candidate_id}: sample_id {spec.sample_id!r}을(를) 이 manifest에서 찾을 수 없습니다")
        value = units[spec.sample_id][0 if slot == "a" else 1]
        if value not in electrodes:
            raise InputContractError(f"{spec.candidate_id}: slot {slot}에 전극이 없습니다(open 불량)")
        return value

    electrode = pick(edit.electrode, edit.slot, "electrode")
    target = pick(edit.target, edit.target_slot, "target")
    if electrode is None:
        raise InputContractError(f"{spec.candidate_id}: electrode 또는 slot을 지정해야 합니다")
    for value in (electrode, target):
        if value is not None and value not in electrodes:
            raise InputContractError(f"{spec.candidate_id}: 알 수 없는 전극 {value!r}")
    if edit.op != "remove_region" and target is None:
        raise InputContractError(f"{spec.candidate_id}: {edit.op}에는 target 또는 target_slot이 필요합니다")
    if target == electrode:
        raise InputContractError(f"{spec.candidate_id}: 편집 전극과 상대 전극이 같습니다")
    return electrode, target


def _components(mask: np.ndarray) -> int:
    return int(label(mask, structure=_EIGHT)[1])


def _write_mask(path: Path, mask: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data = mask.astype(np.uint8) * 255
    if path.suffix.lower() in (".tif", ".tiff"):
        import tifffile

        tifffile.imwrite(path, data)
    else:
        Image.fromarray(data, mode="L").save(path)


def _preview(path: Path, base: np.ndarray, removed: np.ndarray, target: np.ndarray | None,
             roi: np.ndarray) -> None:
    rgb = np.full((*base.shape, 3), 250, dtype=np.uint8)
    rgb[base] = (170, 176, 186)
    if target is not None:
        rgb[target & ~base] = (43, 95, 199)
    rgb[removed] = (217, 60, 50)
    edge = roi & ~(np.roll(roi, 1, 0) & np.roll(roi, -1, 0) & np.roll(roi, 1, 1) & np.roll(roi, -1, 1))
    rgb[edge] = (230, 180, 30)
    image = Image.fromarray(rgb, mode="RGB")
    scale = max(1, -(-600 // max(base.shape)))  # 작은 격자는 600 px 이상으로 확대해 본다.
    if scale > 1:
        image = image.resize((base.shape[1] * scale, base.shape[0] * scale), Image.NEAREST)
    image.save(path)


def _referenced_paths(manifest: ProductManifest) -> list[tuple[str, str]]:
    """(manifest 필드 이름, 상대 경로) 목록."""
    paths = [(f"layer_masks[{i}]", m.path) for i, m in enumerate(manifest.layer_masks)]
    paths += [(f"electrodes[{i}]", e.path) for i, e in enumerate(manifest.electrodes)]
    for i, o in enumerate(manifest.open_paths):
        for side in ("a", "b"):
            value = getattr(o, f"terminal_{side}_mask_path")
            if value is not None:
                paths.append((f"open_paths[{i}].terminal_{side}_mask_path", value))
    paths += [(f"contacts[{i}]", c.path) for i, c in enumerate(manifest.contacts)]
    return paths


def make_variant(manifest_path: Path, spec: VariantSpec, out_dir: Path) -> dict:
    """변경본 패키지 하나를 out_dir/<candidate_id>/에 만들고 기록 dict를 반환한다."""
    manifest_path = Path(manifest_path)
    manifest = load_manifest(manifest_path)
    masks = load_product_masks(manifest, manifest_path)
    grid = masks.grid
    electrodes = {e.electrode_id: e for e in manifest.electrodes}
    units = _sample_units(manifest, masks)
    edited = {k: v.copy() for k, v in masks.electrodes.items()}
    layers = {k: v.copy() for k, v in masks.layers.items()}
    base_electrodes = {k: v.copy() for k, v in masks.electrodes.items()}
    records, touched = [], set()
    removed_total = {}
    for index, edit in enumerate(spec.edits):
        electrode_id, target_id = _resolve(edit, spec, units, electrodes)
        roi = roi_mask(grid, edit.roi_um)
        try:
            result = apply_edit(edit.op, edited[electrode_id], grid, roi,
                                target=None if target_id is None else edited[target_id],
                                amount_um=edit.amount_um)
        except ValueError as exc:
            raise InputContractError(f"{spec.candidate_id} edits[{index}]: {exc}") from exc
        edited[electrode_id] = result.electrode
        layer = electrodes[electrode_id].layer_name
        others = [edited[e] for e, info in electrodes.items() if info.layer_name == layer and e != electrode_id]
        keep = np.logical_or.reduce(others) if others else np.zeros_like(result.removed)
        layers[layer] &= ~(result.removed & ~keep)
        touched.add(electrode_id)
        removed_total[electrode_id] = removed_total.get(electrode_id, np.zeros_like(roi)) | result.removed
        records.append({
            "op": edit.op, "electrode": electrode_id, "target": target_id,
            "roi_um": list(edit.roi_um), "amount_um": edit.amount_um,
            "removed_px": int(result.removed.sum()),
            "removed_area_um2": float(result.removed.sum() * grid.pixel_area_um2),
            "min_center_distance_in_roi_um": {"before": result.min_distance_before_um,
                                               "after": result.min_distance_after_um},
            "overlap_area_um2": None if result.overlap_px_before is None else {
                "before": result.overlap_px_before * grid.pixel_area_um2,
                "after": result.overlap_px_after * grid.pixel_area_um2},
            "_roi": roi, "_target": target_id,
        })

    checks = {"electrodes": {}, "open_paths": "ok", "contact_overlap": []}
    for electrode_id in sorted(touched):
        before, after = base_electrodes[electrode_id], edited[electrode_id]
        if not after.any():
            raise InputContractError(f"{spec.candidate_id}: 편집 후 전극 {electrode_id!r}이(가) 비었습니다")
        checks["electrodes"][electrode_id] = {
            "components_before": _components(before), "components_after": _components(after),
            "split": _components(after) > _components(before),
            "area_um2_before": float(before.sum() * grid.pixel_area_um2),
            "area_um2_after": float(after.sum() * grid.pixel_area_um2),
        }
        for contact_id, contact in (masks.contact_masks or {}).items():
            if (removed_total[electrode_id] & contact).any():
                checks["contact_overlap"].append({"electrode": electrode_id, "contact": contact_id})
    edited_masks = replace(masks, electrodes=edited, layers=layers)
    try:
        check_open_baselines(manifest, edited_masks)
    except InputContractError as exc:
        raise InputContractError(f"{spec.candidate_id}: 편집으로 open 경로 단자 연결이 끊깁니다: {exc}") from exc

    version = spec.design_version or f"{manifest.design_version}-{spec.candidate_id.rsplit('-', 1)[-1]}"
    if version == manifest.design_version:
        raise InputContractError(f"{spec.candidate_id}: 변경본 design_version이 원안과 같습니다")
    target_dir = Path(out_dir) / spec.candidate_id
    if target_dir.exists() and any(target_dir.iterdir()):
        raise InputContractError(f"출력 폴더가 비어 있지 않습니다: {target_dir}")
    target_dir.mkdir(parents=True, exist_ok=True)
    written = {}
    layer_paths = {m.layer_name: m.path for m in manifest.layer_masks}
    for field, rel in _referenced_paths(manifest):
        source = resolve_input_path(manifest_path, rel, product_id=manifest.product_id, field=field)
        dest = target_dir / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        electrode_id = next((e for e, info in electrodes.items() if info.path == rel), None)
        layer_name = next((n for n, p in layer_paths.items() if p == rel), None)
        if electrode_id in touched:
            _write_mask(dest, edited[electrode_id])
        elif layer_name is not None and layer_name in {electrodes[e].layer_name for e in touched}:
            _write_mask(dest, layers[layer_name])
        else:
            shutil.copyfile(source, dest)
        written[rel] = file_sha256(dest)
    raw = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    raw["design_version"] = version
    new_manifest = target_dir / "manifest.yaml"
    new_manifest.write_text(yaml.safe_dump(raw, sort_keys=False, allow_unicode=True), encoding="utf-8")
    # 쓴 결과를 다시 읽어 입력 계약을 통과하는지 확인한다.
    load_product_masks(load_manifest(new_manifest), new_manifest)

    first = records[0]
    target_mask = None if first["_target"] is None else edited[first["_target"]]
    _preview(target_dir / "edit_preview.png", base_electrodes[first["electrode"]],
             removed_total[first["electrode"]], target_mask, first["_roi"])
    record = {
        "candidate_id": spec.candidate_id,
        "product_id": manifest.product_id,
        "parent_design_version": manifest.design_version,
        "parent_design_id": f"{manifest.product_id}-{manifest.design_version}",
        "design_version": version,
        "change_description": spec.change_description,
        "sample_id": spec.sample_id,
        "edits": [{k: v for k, v in r.items() if not k.startswith("_")} for r in records],
        "checks": checks,
        "basis": spec.basis,
        "info": spec.info,
        "parent_manifest_sha256": file_sha256(manifest_path),
        "files_sha256": written,
        "decision": "designer_judgement",
        "note": "국소 편집 후보이며 효과는 재시뮬레이션과 simcompare 재비교로 확인해야 합니다",
    }
    (target_dir / "variant_record.json").write_text(
        json.dumps(record, indent=2, ensure_ascii=False), encoding="utf-8")
    return {**record, "manifest": str(new_manifest)}


def make_variants(manifest_path: Path, spec_path: Path, out_dir: Path,
                  candidate_ids: list[str] | None = None) -> list[dict]:
    """지정서의 후보(또는 candidate_ids로 고른 후보)마다 변경본 패키지를 만든다."""
    spec = load_edit_spec(Path(spec_path))
    chosen = [c for c in spec.candidates if not candidate_ids or c.candidate_id in candidate_ids]
    missing = sorted(set(candidate_ids or []) - {c.candidate_id for c in spec.candidates})
    if missing:
        raise InputContractError(f"지정서에 없는 후보: {missing}", file=str(spec_path))
    ids = [c.candidate_id for c in chosen]
    if len(set(ids)) != len(ids):
        raise InputContractError("candidate_id가 중복됩니다", file=str(spec_path))
    return [make_variant(Path(manifest_path), c, Path(out_dir)) for c in chosen]
