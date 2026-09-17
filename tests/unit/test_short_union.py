from __future__ import annotations

import numpy as np
import pytest

from particle_sim.config import load_manifest
from particle_sim.errors import InputContractError
from particle_sim.failure.lateral_short import lateral_short
from particle_sim.failure.union import ANY_MODE, simulate_pairs, union_by_mode
from particle_sim.failure.vertical_short import vertical_short
from particle_sim.geometry.contact import DistanceCache, compute_electrode_contacts
from particle_sim.io.masks import load_product_masks
from particle_sim.io.stack import build_stack_z
from particle_sim.validation.synthetic_cases import build_case, write_case


def _load(tmp_path_factory: pytest.TempPathFactory, case: str):
    path = write_case(build_case(case), tmp_path_factory.mktemp(case))
    manifest = load_manifest(path)
    masks = load_product_masks(manifest, path)
    return manifest, masks, build_stack_z(manifest), DistanceCache(masks)


def _failures(bundle, step_id: str, diameter: float) -> dict[str, np.ndarray]:
    manifest, masks, stack_z, cache = bundle
    results = simulate_pairs(manifest, masks, stack_z, step_id, diameter, cache)
    return {r.pair_id: r.failure for r in results}


def test_lateral_short_between_stripes(tmp_path_factory: pytest.TempPathFactory) -> None:
    bundle = _load(tmp_path_factory, "parallel_lines")
    small = _failures(bundle, "S_M2", 2.0)["A__B"]
    large = _failures(bundle, "S_M2", 10.0)["A__B"]
    assert not small.any()
    assert large.any()
    columns = np.where(large.any(axis=0))[0]
    assert 90 <= columns.min() and columns.max() < 96  # strictly inside the gap


def test_hotspot_is_monotone_in_size(tmp_path_factory: pytest.TempPathFactory) -> None:
    bundle = _load(tmp_path_factory, "parallel_lines")
    previous = None
    for diameter in (2.0, 4.0, 6.0, 8.0, 10.0):
        current = _failures(bundle, "S_M2", diameter)["A__B"]
        if previous is not None:
            assert np.all(previous <= current)
        previous = current


def test_diagonal_hotspot_moves_toward_upper_electrode(
    tmp_path_factory: pytest.TempPathFactory,
) -> None:
    aligned = _failures(_load(tmp_path_factory, "vertical_aligned"), "S_M1", 10.0)["G__P"]
    diagonal = _failures(_load(tmp_path_factory, "diagonal_offset"), "S_M1", 10.0)["G__P"]
    assert aligned.any() and diagonal.any()
    assert np.where(diagonal)[1].mean() > np.where(aligned)[1].mean()


def test_exclude_policy_removes_post_step_contacts(
    tmp_path_factory: pytest.TempPathFactory,
) -> None:
    manifest, masks, stack_z, cache = _load(tmp_path_factory, "vertical_aligned")
    excluded = manifest.model_copy(
        update={"post_step_policy": manifest.post_step_policy.model_copy(
            update={"mode": "exclude"})}
    )
    results = simulate_pairs(excluded, masks, stack_z, "S_M1", 10.0, cache)
    assert not results[0].failure.any()


def test_union_removes_overlap(tmp_path_factory: pytest.TempPathFactory) -> None:
    bundle = _load(tmp_path_factory, "triple_overlap")
    manifest, masks, stack_z, cache = bundle
    large = simulate_pairs(manifest, masks, stack_z, "S_M2", 10.0, cache)
    unions = union_by_mode(large)
    lateral = unions["lateral_short"]
    pair_sum = sum(int(r.failure.sum()) for r in large)
    assert int(lateral.any_failure.sum()) < pair_sum
    assert lateral.pair_count.max() == 3
    assert set(unions) == {"lateral_short", ANY_MODE}

    small = simulate_pairs(manifest, masks, stack_z, "S_M2", 6.0, cache)
    small_union = union_by_mode(small)["lateral_short"]
    assert not next(r for r in small if r.pair_id == "A__C").failure.any()
    assert int(small_union.any_failure.sum()) == sum(int(r.failure.sum()) for r in small)


def test_pair_margin_and_classification(tmp_path_factory: pytest.TempPathFactory) -> None:
    manifest, masks, stack_z, cache = _load(tmp_path_factory, "pilot")
    contacts = compute_electrode_contacts(
        manifest, masks, stack_z, "S_M1", 10.0, cache, ["G", "A", "B", "P"]
    )
    pairs = {p.pair_id: p for p in manifest.electrode_pairs}
    lateral = lateral_short(manifest, pairs["A__B"], contacts)
    vertical = vertical_short(manifest, pairs["G__A"], contacts)
    assert lateral.classification == "direct_geometric_contact"
    assert vertical.classification == "potential_vertical_short"
    expected = np.minimum(contacts["A"].margin_um, contacts["B"].margin_um)
    assert np.array_equal(lateral.pair_margin_um, expected)


def test_wrong_mode_or_relation_rejected(tmp_path_factory: pytest.TempPathFactory) -> None:
    manifest, masks, stack_z, cache = _load(tmp_path_factory, "pilot")
    contacts = compute_electrode_contacts(
        manifest, masks, stack_z, "S_M1", 4.0, cache, ["G", "A", "B", "P"]
    )
    pairs = {p.pair_id: p for p in manifest.electrode_pairs}
    with pytest.raises(InputContractError):
        vertical_short(manifest, pairs["A__B"], contacts)  # same layer pair
    with pytest.raises(InputContractError):
        lateral_short(manifest, pairs["G__A"], contacts)  # different layers
    same_net = manifest.model_copy(
        update={"electrodes": [e.model_copy(update={"net_id": "SAME"})
                               for e in manifest.electrodes]}
    )
    with pytest.raises(InputContractError):
        lateral_short(same_net, pairs["A__B"], contacts)


def test_union_requires_results() -> None:
    with pytest.raises(ValueError):
        union_by_mode([])
