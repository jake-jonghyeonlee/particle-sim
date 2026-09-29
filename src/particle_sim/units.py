"""Common view over the evaluated failure units (electrode pairs and open paths)."""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from particle_sim.failure.auto_open import simulate_auto_open
from particle_sim.failure.contact_not_open import CNO_MODE, simulate_cno
from particle_sim.failure.lateral_short import PairResult
from particle_sim.failure.open_failure import simulate_open
from particle_sim.failure.union import simulate_pairs
from particle_sim.geometry.contact import DistanceCache
from particle_sim.io.masks import MaskSet
from particle_sim.io.stack import LayerZ
from particle_sim.schema import AUTO_OPEN_PREFIX, ProductManifest

OPEN_MODE = "lateral_open"
NOT_APPLICABLE = "N/A"
NOT_APPLICABLE_STEP = "not_applicable_source_step"
# key: (unit id, particle diameter, contact margin). Values stay valid only while the manifest
# and masks are identical; see evaluate_units.
OpenResultCache = dict[tuple[str, float, float], PairResult]

KIND_PAIR = "electrode_pair"
KIND_OPEN_PATH = "open_path"
KIND_OPEN_ELECTRODE = "open_electrode"
KIND_CONTACT = "contact"


@dataclass(frozen=True)
class UnitSpec:
    """One failure unit: an electrode pair, an open path/electrode or a contact link."""

    unit_id: str
    failure_mode: str
    unit_kind: str
    electrode_a: str
    electrode_b: str  # NOT_APPLICABLE for open units
    applicable_steps: tuple[str, ...] | None  # None: every source step (v1 open paths, pairs)


def unit_specs(manifest: ProductManifest, masks: MaskSet | None = None) -> list[UnitSpec]:
    """Pairs, explicit open paths, automatic open electrodes, then CNO links.

    Automatic open and CNO units come from the nominal contact topology, so ``masks`` must be
    the nominal MaskSet whenever the manifest declares them.
    """
    specs = [
        UnitSpec(p.pair_id, p.failure_mode, KIND_PAIR, p.electrode_a, p.electrode_b, None)
        for p in manifest.electrode_pairs
    ]
    explicit = {o.electrode for o in manifest.open_paths}
    for open_path in manifest.open_paths:
        specs.append(
            UnitSpec(
                open_path.path_id,
                OPEN_MODE,
                KIND_OPEN_PATH,
                open_path.electrode,
                NOT_APPLICABLE,
                tuple(open_path.applicable_particle_source_steps) or None,
            )
        )
    for layer_name, steps in manifest.open_applicable_steps_by_layer.items():
        if not steps:
            continue
        for electrode in manifest.electrodes:
            if electrode.layer_name != layer_name or electrode.electrode_id in explicit:
                continue
            specs.append(
                UnitSpec(
                    f"{AUTO_OPEN_PREFIX}{electrode.electrode_id}",
                    OPEN_MODE,
                    KIND_OPEN_ELECTRODE,
                    electrode.electrode_id,
                    NOT_APPLICABLE,
                    tuple(steps),
                )
            )
    if manifest.contacts:
        if masks is None or masks.contact_topology is None:
            raise ValueError("contact가 있는 manifest의 불량 단위에는 명목 MaskSet이 필요합니다")
        steps_by_contact = {
            c.contact_id: tuple(c.applicable_particle_source_steps) for c in manifest.contacts
        }
        for link in masks.contact_topology.links.values():
            specs.append(
                UnitSpec(
                    link.unit_id,
                    CNO_MODE,
                    KIND_CONTACT,
                    link.lower_electrode,
                    link.upper_electrode,
                    steps_by_contact[link.contact_id],
                )
            )
    return specs


def unit_ids(manifest: ProductManifest, masks: MaskSet | None = None) -> list[str]:
    return [spec.unit_id for spec in unit_specs(manifest, masks)]


def unit_modes(manifest: ProductManifest, masks: MaskSet | None = None) -> list[str]:
    return [spec.failure_mode for spec in unit_specs(manifest, masks)]


def _empty_result(
    spec: UnitSpec, shape: tuple[int, int], classification: str,
    placeholders: tuple[np.ndarray, np.ndarray] | None = None,
) -> PairResult:
    failure = np.zeros(shape, dtype=bool) if placeholders is None else placeholders[0]
    return _single_electrode_result(spec, failure, classification, placeholders)


def _single_electrode_result(
    spec: UnitSpec, failure: np.ndarray, classification: str,
    placeholders: tuple[np.ndarray, np.ndarray] | None = None,
) -> PairResult:
    shape = failure.shape
    contact_b = np.zeros(shape, dtype=bool) if placeholders is None else placeholders[0]
    margin = (
        np.full(shape, np.nan, dtype=np.float32)
        if placeholders is None else placeholders[1]
    )
    return PairResult(
        pair_id=spec.unit_id,
        failure_mode=spec.failure_mode,
        electrode_a=spec.electrode_a,
        electrode_b=spec.electrode_b,
        failure=failure,
        contact_a=failure,
        contact_b=contact_b,
        pair_margin_um=margin,
        classification=classification,
    )


def evaluate_units(
    manifest: ProductManifest,
    masks: MaskSet,
    stack_z: dict[str, LayerZ],
    step_id: str,
    diameter_um: float,
    cache: DistanceCache,
    open_cache: OpenResultCache | None = None,
    specs: list[UnitSpec] | None = None,
    *,
    shared_readonly_placeholders: bool = False,
) -> list[PairResult]:
    """Every failure unit for one source step and particle size, in ``specs`` order.

    ``specs`` defaults to the units of ``masks``; transformed scenarios pass the nominal specs
    so a contact link that disappears is reported as lost instead of silently dropped.

    ``open_cache`` may be shared only while the manifest and masks stay identical. Open and CNO
    geometry does not depend on the particle source step, so applicable steps reuse the same
    unit/diameter result. A caller evaluating transformed scenarios must isolate caches by
    geometry and bound their memory use.

    ``shared_readonly_placeholders`` shares zero and NaN maps within this evaluation.
    These placeholder arrays are read-only; the default keeps independent mutable maps.
    """
    if specs is None:
        specs = unit_specs(manifest, masks)
    shape = masks.grid.shape_yx
    placeholders = None
    if shared_readonly_placeholders:
        # A bytes owner cannot be made writeable through an ndarray alias. Build the
        # same C-order bytes as the default constructors, scoped to this evaluation.
        zero = np.frombuffer(np.zeros(shape, dtype=bool).tobytes(), dtype=bool).reshape(shape)
        nan = np.frombuffer(
            np.full(shape, np.nan, dtype=np.float32).tobytes(), dtype=np.float32
        ).reshape(shape)
        placeholders = (zero, nan)
    margin = float(manifest.contact_margin_um)
    shorts = {r.pair_id: r for r in simulate_pairs(manifest, masks, stack_z, step_id, diameter_um, cache)}
    open_paths = {o.path_id: o for o in manifest.open_paths}
    links = {} if masks.contact_topology is None else masks.contact_topology.links
    results: list[PairResult] = []
    for spec in specs:
        if spec.unit_kind == KIND_PAIR:
            results.append(shorts[spec.unit_id])
            continue
        if spec.applicable_steps is not None and step_id not in spec.applicable_steps:
            results.append(_empty_result(spec, shape, NOT_APPLICABLE_STEP, placeholders))
            continue
        cache_key = (spec.unit_id, float(diameter_um), margin)
        cached = None if open_cache is None else open_cache.get(cache_key)
        if cached is None:
            if spec.unit_kind == KIND_OPEN_PATH:
                opened = simulate_open(
                    manifest, masks, open_paths[spec.unit_id], diameter_um, contact_margin_um=margin
                )
                cached = _single_electrode_result(
                    spec, opened.failure, opened.classification, placeholders
                )
            elif spec.unit_kind == KIND_OPEN_ELECTRODE:
                auto = simulate_auto_open(masks, spec.electrode_a, diameter_um, margin)
                cached = _single_electrode_result(
                    spec, auto.failure, auto.classification, placeholders
                )
            else:
                cno = simulate_cno(
                    links.get(spec.unit_id), spec.unit_id, masks.grid, diameter_um, margin
                )
                cached = _single_electrode_result(
                    spec, cno.failure, cno.classification, placeholders
                )
            if open_cache is not None:
                open_cache[cache_key] = cached
        elif placeholders is not None:
            # Reuse the cached geometry while moving its placeholders to this
            # evaluation's pool. Otherwise cached units retain older pools here.
            cached = replace(cached, contact_b=placeholders[0], pair_margin_um=placeholders[1])
            if open_cache is not None:
                open_cache[cache_key] = cached
        elif not cached.contact_b.flags.writeable or not cached.pair_margin_um.flags.writeable:
            # A caller may reuse a cache populated by an opt-in evaluation.
            # The default contract still returns independently mutable maps.
            cached = replace(
                cached,
                contact_b=np.zeros(shape, dtype=bool),
                pair_margin_um=np.full(shape, np.nan, dtype=np.float32),
            )
            if open_cache is not None:
                open_cache[cache_key] = cached
        results.append(cached)
    return results
