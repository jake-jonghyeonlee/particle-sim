"""Contact·via topology: which hole joins which lower and upper electrode."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import label

from particle_sim.errors import InputContractError
from particle_sim.schema import CNO_PREFIX, ProductManifest, resolve_contact_layers


@dataclass(frozen=True)
class ContactHole:
    contact_id: str
    hole_index: int  # 1-based connected component of the contact mask
    lower_electrode: str
    upper_electrode: str
    # effective contact pixels: hole ∩ lower electrode ∩ upper electrode
    ys: np.ndarray
    xs: np.ndarray
    # Physical TR-cut component, distinct from the effective endpoint above.
    cut_ys: np.ndarray
    cut_xs: np.ndarray
    raw_pixels: int | None  # None when the input supplies only a cut mask
    cut_pixels: int
    effective_pixels: int


@dataclass(frozen=True)
class ContactLink:
    """One CNO failure unit: every hole of one contact that joins the same two electrodes."""

    unit_id: str
    contact_id: str
    lower_electrode: str
    upper_electrode: str
    holes: tuple[ContactHole, ...]


@dataclass(frozen=True)
class ContactTopology:
    holes: tuple[ContactHole, ...]
    links: dict[str, ContactLink]
    lost_holes: int = 0  # holes that no longer join two electrodes (transformed scenarios)


def cno_unit_id(contact_id: str, lower_electrode: str, upper_electrode: str) -> str:
    return f"{CNO_PREFIX}{contact_id}:{lower_electrode}__{upper_electrode}"


def _electrode_index(
    manifest: ProductManifest, electrodes: dict[str, np.ndarray], layer_name: str, shape
) -> tuple[np.ndarray, list[str]]:
    ids = [e.electrode_id for e in manifest.electrodes if e.layer_name == layer_name]
    index = np.zeros(shape, dtype=np.int32)
    for position, electrode_id in enumerate(ids, start=1):
        index[electrodes[electrode_id]] = position
    return index, ids


def build_contact_topology(
    manifest: ProductManifest,
    contact_masks: dict[str, np.ndarray],
    electrodes: dict[str, np.ndarray],
    *,
    strict: bool,
    manifest_path: str | None = None,
    raw_contact_masks: dict[str, np.ndarray] | None = None,
) -> ContactTopology:
    """Map every contact hole to its lower and upper electrode.

    ``strict`` is used for the nominal input: a hole that does not join exactly one electrode
    below and one above, of the same net, is an input error. Transformed scenarios pass
    ``strict=False`` and count such holes as lost instead of failing the whole run.
    """
    net_of = {e.electrode_id: e.net_id for e in manifest.electrodes}
    holes: list[ContactHole] = []
    lost = 0
    index_cache: dict[str, tuple[np.ndarray, list[str]]] = {}
    for position, contact in enumerate(manifest.contacts):
        field = f"contacts[{position}]"
        mask = contact_masks[contact.contact_id]
        raw_labels = None
        if raw_contact_masks is not None:
            raw_mask = raw_contact_masks[contact.contact_id]
            if raw_mask.shape != mask.shape or np.any(mask & ~raw_mask):
                raise InputContractError(
                    f"contact '{contact.contact_id}'의 cut mask가 raw mask와 맞지 않습니다",
                    product_id=manifest.product_id,
                    file=manifest_path,
                    field=f"{field}.path",
                )
            raw_labels, _ = label(raw_mask)
        lower_layer, upper_layer = resolve_contact_layers(manifest.stack, contact)
        assert lower_layer is not None and upper_layer is not None
        for layer_name in (lower_layer, upper_layer):
            if layer_name not in index_cache:
                index_cache[layer_name] = _electrode_index(
                    manifest, electrodes, layer_name, mask.shape
                )
        lower_index, lower_ids = index_cache[lower_layer]
        upper_index, upper_ids = index_cache[upper_layer]
        labels, count = label(mask)
        if strict and count == 0:
            raise InputContractError(
                f"contact '{contact.contact_id}' mask가 비어 있습니다",
                product_id=manifest.product_id,
                file=manifest_path,
                field=f"{field}.path",
            )
        for hole_index in range(1, count + 1):
            ys, xs = np.nonzero(labels == hole_index)
            if (ys.min() == 0 or xs.min() == 0 or
                    ys.max() == mask.shape[0] - 1 or xs.max() == mask.shape[1] - 1):
                if strict:
                    raise InputContractError(
                        f"contact '{contact.contact_id}'의 hole {hole_index}가 ROI edge에 닿아 "
                        "전체 physical hole 면적을 확인할 수 없습니다",
                        product_id=manifest.product_id,
                        file=manifest_path,
                        field=f"{field}.path",
                    )
                lost += 1
                continue
            raw_pixels = None
            if raw_labels is not None:
                raw_hits = set(raw_labels[ys, xs].tolist()) - {0}
                if len(raw_hits) != 1:
                    raise InputContractError(
                        f"contact '{contact.contact_id}'의 cut hole {hole_index}가 여러 raw hole을 합칩니다",
                        product_id=manifest.product_id,
                        file=manifest_path,
                        field=f"{field}.path",
                    )
                raw_pixels = int(np.count_nonzero(raw_labels == raw_hits.pop()))
            lower_hits = sorted(set(lower_index[ys, xs].tolist()) - {0})
            upper_hits = sorted(set(upper_index[ys, xs].tolist()) - {0})
            problem: str | None = None
            if len(lower_hits) != 1 or len(upper_hits) != 1:
                problem = (
                    f"아래 전극 {len(lower_hits)}개, 위 전극 {len(upper_hits)}개와 겹칩니다"
                    " (각각 정확히 1개여야 함)"
                )
            else:
                lower_id = lower_ids[lower_hits[0] - 1]
                upper_id = upper_ids[upper_hits[0] - 1]
                if net_of[lower_id] != net_of[upper_id]:
                    problem = f"서로 다른 net 전극 {lower_id}, {upper_id}을(를) 잇습니다"
                else:
                    effective = (lower_index[ys, xs] == lower_hits[0]) & (
                        upper_index[ys, xs] == upper_hits[0]
                    )
                    if not effective.any():
                        problem = "위·아래 전극이 hole 안에서 함께 겹치는 픽셀이 없습니다"
            if problem is not None:
                if strict:
                    first = (int(ys[0]), int(xs[0]))
                    raise InputContractError(
                        f"contact '{contact.contact_id}'의 hole {hole_index}(첫 픽셀 {first})이 "
                        + problem,
                        product_id=manifest.product_id,
                        file=manifest_path,
                        field=f"{field}.path",
                    )
                lost += 1
                continue
            holes.append(
                ContactHole(
                    contact_id=contact.contact_id,
                    hole_index=hole_index,
                    lower_electrode=lower_id,
                    upper_electrode=upper_id,
                    ys=ys[effective],
                    xs=xs[effective],
                    cut_ys=ys,
                    cut_xs=xs,
                    raw_pixels=raw_pixels,
                    cut_pixels=len(ys),
                    effective_pixels=int(np.count_nonzero(effective)),
                )
            )
    links: dict[str, list[ContactHole]] = {}
    for hole in holes:
        unit_id = cno_unit_id(hole.contact_id, hole.lower_electrode, hole.upper_electrode)
        links.setdefault(unit_id, []).append(hole)
    return ContactTopology(
        holes=tuple(holes),
        links={
            unit_id: ContactLink(
                unit_id=unit_id,
                contact_id=members[0].contact_id,
                lower_electrode=members[0].lower_electrode,
                upper_electrode=members[0].upper_electrode,
                holes=tuple(members),
            )
            for unit_id, members in links.items()
        },
        lost_holes=lost,
    )
