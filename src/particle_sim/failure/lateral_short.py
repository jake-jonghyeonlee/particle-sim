"""Lateral short between two nets on the same layer."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from particle_sim.errors import InputContractError
from particle_sim.geometry.contact import ContactResult
from particle_sim.schema import ElectrodePair, ProductManifest


@dataclass(frozen=True)
class PairResult:
    pair_id: str
    failure_mode: str
    electrode_a: str
    electrode_b: str
    failure: np.ndarray
    contact_a: np.ndarray
    contact_b: np.ndarray
    pair_margin_um: np.ndarray
    classification: str


def build_pair_result(
    manifest: ProductManifest,
    pair: ElectrodePair,
    contacts: dict[str, ContactResult],
    *,
    expected_mode: str,
    same_layer: bool,
    classification: str,
) -> PairResult:
    """Shared checks and map algebra for both short modes."""
    field = f"electrode_pairs[{pair.pair_id}]"
    if pair.failure_mode != expected_mode:
        raise InputContractError(
            f"pair '{pair.pair_id}'의 failure_mode가 {expected_mode}가 아닙니다: "
            f"{pair.failure_mode}",
            product_id=manifest.product_id,
            field=field,
        )
    electrodes = {e.electrode_id: e for e in manifest.electrodes}
    for electrode_id in (pair.electrode_a, pair.electrode_b):
        if electrode_id not in electrodes:
            raise InputContractError(
                f"pair '{pair.pair_id}'의 전극 '{electrode_id}'이(가) 없습니다",
                product_id=manifest.product_id,
                field=field,
            )
    ea = electrodes[pair.electrode_a]
    eb = electrodes[pair.electrode_b]
    if same_layer and ea.layer_name != eb.layer_name:
        raise InputContractError(
            f"lateral_short pair '{pair.pair_id}'의 두 전극은 같은 layer여야 합니다",
            product_id=manifest.product_id,
            field=field,
        )
    if not same_layer and ea.layer_name == eb.layer_name:
        raise InputContractError(
            f"vertical_short pair '{pair.pair_id}'의 두 전극은 다른 layer여야 합니다",
            product_id=manifest.product_id,
            field=field,
        )
    if ea.net_id == eb.net_id:
        raise InputContractError(
            f"pair '{pair.pair_id}'의 두 전극 net_id가 같습니다: {ea.net_id}",
            product_id=manifest.product_id,
            field=field,
        )
    for electrode_id in (pair.electrode_a, pair.electrode_b):
        if electrode_id not in contacts:
            raise InputContractError(
                f"pair '{pair.pair_id}'의 전극 '{electrode_id}' 접촉 결과가 없습니다",
                product_id=manifest.product_id,
                field=field,
            )

    contact_a = contacts[pair.electrode_a]
    contact_b = contacts[pair.electrode_b]
    return PairResult(
        pair_id=pair.pair_id,
        failure_mode=pair.failure_mode,
        electrode_a=pair.electrode_a,
        electrode_b=pair.electrode_b,
        failure=contact_a.contact & contact_b.contact,
        contact_a=contact_a.contact,
        contact_b=contact_b.contact,
        pair_margin_um=np.minimum(contact_a.margin_um, contact_b.margin_um).astype(np.float32),
        classification=classification,
    )


def lateral_short(
    manifest: ProductManifest, pair: ElectrodePair, contacts: dict[str, ContactResult]
) -> PairResult:
    """Particle touches two different nets of the same layer."""
    return build_pair_result(
        manifest,
        pair,
        contacts,
        expected_mode="lateral_short",
        same_layer=True,
        classification="direct_geometric_contact",
    )
