"""이진 마스크 로더와 제품 MaskSet 검증."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import tifffile
from PIL import Image

from particle_sim.config import file_sha256, resolve_input_path
from particle_sim.errors import InputContractError
from particle_sim.schema import ProductManifest


@dataclass(frozen=True)
class GridSpec:
    """이미지 격자 사양."""

    shape_yx: tuple[int, int]
    pixel_size_y_um: float
    pixel_size_x_um: float
    origin: str = "upper_left"

    @property
    def pixel_area_um2(self) -> float:
        """1픽셀의 면적 (um^2)."""
        return self.pixel_size_y_um * self.pixel_size_x_um

    @property
    def roi_area_um2(self) -> float:
        """ROI 전체 면적 (um^2)."""
        return self.shape_yx[0] * self.pixel_size_y_um * self.shape_yx[1] * self.pixel_size_x_um

    @property
    def is_isotropic(self) -> bool:
        """x/y 픽셀 크기가 동일한지."""
        return abs(self.pixel_size_y_um - self.pixel_size_x_um) <= 1e-9


@dataclass(frozen=True)
class MaskSet:
    """제품 전체 마스크 집합."""

    grid: GridSpec
    layers: dict[str, np.ndarray]
    electrodes: dict[str, np.ndarray]
    input_hashes: dict[str, str]


def load_binary_mask(path: Path, *, product_id: str, field: str) -> np.ndarray:
    """이진 마스크 파일을 읽어 bool (y, x) 배열을 반환한다.

    .png: Pillow, .tif/.tiff: tifffile.
    2-D 단일 채널만 허용. 허용 값: {0,1} 또는 {0,255} (bool 포함).
    resize/threshold 금지.
    """
    suffix = path.suffix.lower()
    if suffix == ".png":
        try:
            img = Image.open(path)
            img.load()
        except (OSError, ValueError) as exc:
            raise InputContractError(
                f"PNG 마스크를 읽을 수 없습니다: {exc}",
                product_id=product_id,
                file=str(path),
                field=field,
            ) from exc
        if img.mode not in ("L", "1"):
            raise InputContractError(
                f"2-D 단일 채널 PNG가 아닙니다 (mode={img.mode})",
                product_id=product_id,
                file=str(path),
                field=field,
            )
        arr = np.array(img)
    elif suffix in (".tif", ".tiff"):
        try:
            arr = tifffile.imread(str(path))
        except (OSError, ValueError, KeyError, IndexError) as exc:
            raise InputContractError(
                f"TIFF 마스크를 읽을 수 없습니다: {exc}",
                product_id=product_id,
                file=str(path),
                field=field,
            ) from exc
        if arr.ndim != 2:
            raise InputContractError(
                f"2-D 단일 채널 TIFF가 아닙니다 (ndim={arr.ndim})",
                product_id=product_id,
                file=str(path),
                field=field,
            )
    else:
        raise InputContractError(
            f"지원하지 않는 확장자: {suffix} (허용: .png, .tif, .tiff)",
            product_id=product_id,
            file=str(path),
            field=field,
        )

    if arr.dtype == bool:
        return arr

    unique = np.unique(arr)
    allowed_01 = set(unique.tolist()) <= {0, 1}
    allowed_0255 = set(unique.tolist()) <= {0, 255}
    if not (allowed_01 or allowed_0255):
        bad_count = int(np.sum(~np.isin(arr, [0, 1, 255])))
        raise InputContractError(
            f"이진 마스크가 아님(resize/antialias 의심): {bad_count}개의 유효하지 않은 값",
            product_id=product_id,
            file=str(path),
            field=field,
        )

    return arr != 0


def load_product_masks(manifest: ProductManifest, manifest_path: Path) -> MaskSet:
    """manifest의 모든 layer/electrode 마스크를 로드하고 계약을 검증한다."""
    cs = manifest.coordinate_system
    product_id = manifest.product_id

    # 모든 마스크 경로와 field를 수집
    entries: list[tuple[str, str, str]] = []  # (relative_path, field, kind)
    for index, lm in enumerate(manifest.layer_masks):
        entries.append((lm.path, f"layer_masks[{index}].path", "layer"))
    for index, el in enumerate(manifest.electrodes):
        entries.append((el.path, f"electrodes[{index}].path", "electrode"))

    if not entries:
        raise InputContractError(
            "마스크가 없습니다",
            product_id=product_id,
            file=str(manifest_path),
            field="layer_masks",
        )

    electrode_index = {e.electrode_id: i for i, e in enumerate(manifest.electrodes)}
    layer_index = {lm.layer_name: i for i, lm in enumerate(manifest.layer_masks)}

    # 1. 모든 마스크 shape 일치
    ref_shape: tuple[int, int] | None = None
    ref_file: str | None = None
    masks: dict[str, np.ndarray] = {}
    input_hashes: dict[str, str] = {}

    for rel_path, field, _kind in entries:
        resolved = resolve_input_path(
            manifest_path, rel_path, product_id=product_id, field=field
        )
        arr = load_binary_mask(resolved, product_id=product_id, field=field)
        try:
            input_hashes[rel_path] = file_sha256(resolved)
        except OSError as exc:
            raise InputContractError(
                f"마스크 해시를 계산할 수 없습니다: {exc}",
                product_id=product_id,
                file=str(resolved),
                field=field,
            ) from exc
        shape = arr.shape
        if ref_shape is None:
            ref_shape = shape
            ref_file = rel_path
        elif shape != ref_shape:
            raise InputContractError(
                f"shape 불일치: '{rel_path}'={shape} vs '{ref_file}'={ref_shape}",
                product_id=product_id,
                file=str(resolved),
                field=field,
            )
        masks[rel_path] = arr

    assert ref_shape is not None
    shape_y, shape_x = ref_shape

    # 2. ROI 크기 검증
    roi_w = shape_x * cs.pixel_size_x_um
    roi_h = shape_y * cs.pixel_size_y_um
    if abs(roi_w - cs.roi_width_um) > cs.pixel_size_x_um / 2:
        raise InputContractError(
            f"ROI 폭 불일치: 계산={roi_w}um, manifest={cs.roi_width_um}um",
            product_id=product_id,
            file=str(manifest_path),
            field="coordinate_system.roi_width_um",
        )
    if abs(roi_h - cs.roi_height_um) > cs.pixel_size_y_um / 2:
        raise InputContractError(
            f"ROI 높이 불일치: 계산={roi_h}um, manifest={cs.roi_height_um}um",
            product_id=product_id,
            file=str(manifest_path),
            field="coordinate_system.roi_height_um",
        )

    # layer_name -> mask, electrode_id -> mask, electrode_id -> layer_name 매핑
    layers: dict[str, np.ndarray] = {}
    electrodes: dict[str, np.ndarray] = {}
    electrode_layer: dict[str, str] = {}

    for lm in manifest.layer_masks:
        layers[lm.layer_name] = masks[lm.path]
    for el in manifest.electrodes:
        electrodes[el.electrode_id] = masks[el.path]
        electrode_layer[el.electrode_id] = el.layer_name

    # 3. 빈 전극 검증
    for eid, emask in electrodes.items():
        if not emask.any():
            raise InputContractError(
                f"전극 '{eid}' 마스크가 비어있습니다",
                product_id=product_id,
                file=str(manifest_path),
                field=f"electrodes[{electrode_index[eid]}].path",
            )

    # 4. 동일 레이어 전극 겹침 검증
    layer_electrodes: dict[str, list[str]] = {}
    for eid, lname in electrode_layer.items():
        layer_electrodes.setdefault(lname, []).append(eid)

    for lname, eids in layer_electrodes.items():
        if len(eids) < 2:
            continue
        # 누적 OR으로 겹침 감지
        accumulated = np.zeros(shape, dtype=bool)
        for eid in eids:
            em = electrodes[eid]
            overlap = accumulated & em
            if overlap.any():
                raise InputContractError(
                    f"레이어 '{lname}'의 전극이 겹칩니다: {int(overlap.sum())}픽셀",
                    product_id=product_id,
                    file=str(manifest_path),
                    field=f"electrodes[{electrode_index[eid]}].path",
                )
            accumulated |= em

    # 5. 레이어 합집합 일치 검증
    for lname, lmask in layers.items():
        eids = layer_electrodes.get(lname, [])
        if not eids:
            continue  # topography only
        union = np.zeros(shape, dtype=bool)
        for eid in eids:
            union |= electrodes[eid]
        missing = int(np.sum(lmask & ~union))  # layer에 있지만 union에 없음
        extra = int(np.sum(union & ~lmask))  # union에 있지만 layer에 없음
        if missing or extra:
            raise InputContractError(
                f"레이어 '{lname}' 합집합 불일치: missing={missing}, extra={extra}",
                product_id=product_id,
                file=str(manifest_path),
                field=f"layer_masks[{layer_index[lname]}].path",
            )

    grid = GridSpec(
        shape_yx=(shape_y, shape_x),
        pixel_size_y_um=cs.pixel_size_y_um,
        pixel_size_x_um=cs.pixel_size_x_um,
        origin=cs.origin,
    )

    return MaskSet(
        grid=grid,
        layers=layers,
        electrodes=electrodes,
        input_hashes=input_hashes,
    )
