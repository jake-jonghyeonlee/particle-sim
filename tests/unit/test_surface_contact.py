from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from particle_sim.config import load_manifest
from particle_sim.errors import InputContractError
from particle_sim.geometry.contact import (
    DistanceCache,
    compute_electrode_contacts,
    contact_3d,
    contact_post_step,
)
from particle_sim.geometry.distance import distance_to_mask_um
from particle_sim.geometry.particles import particle_center_z_um, post_step_eligibility
from particle_sim.geometry.surfaces import surface_height_um
from particle_sim.io.masks import load_product_masks
from particle_sim.io.stack import build_stack_z
from particle_sim.validation.synthetic_cases import build_case, write_case


@pytest.fixture(scope="module")
def pilot(tmp_path_factory: pytest.TempPathFactory):
    path = write_case(build_case("pilot"), tmp_path_factory.mktemp("surface"))
    manifest = load_manifest(path)
    masks = load_product_masks(manifest, path)
    return manifest, masks, build_stack_z(manifest)


def test_flat_surface_is_layer_top(pilot) -> None:
    manifest, masks, stack_z = pilot
    surface = surface_height_um(manifest, "S_PLN", stack_z, masks)
    assert np.unique(surface).tolist() == [pytest.approx(3.0)]


def test_conformal_surface_follows_missing_layers(pilot) -> None:
    manifest, masks, stack_z = pilot
    surface = surface_height_um(manifest, "S_M2", stack_z, masks)
    # over the M2 electrode A (x 66..76 px): M1 missing (-0.3), M2 present -> 1.0 - 0.3
    assert surface[0, 70] == pytest.approx(0.7)
    # bare region: M1 and M2 missing -> 1.0 - 0.3 - 0.3
    assert surface[0, 5] == pytest.approx(0.4)
    # over the M1 gate (x 40..60 px): only M2 missing
    assert surface[0, 50] == pytest.approx(0.7)


def test_partial_surface_uses_planarization(pilot) -> None:
    manifest, masks, stack_z = pilot
    partial = manifest.model_copy(
        update={
            "particle_source_steps": [
                s.model_copy(update={"surface_model": "partial"})
                for s in manifest.particle_source_steps
            ]
        }
    )
    surface = surface_height_um(partial, "S_M2", stack_z, masks)
    # M1 topography is damped by GI (planarization 0.0 -> retention 1.0); M2 has nothing above
    assert surface[0, 5] == pytest.approx(0.4)


def test_particle_center_z_and_embedding_limits() -> None:
    surface = np.zeros((2, 2))
    assert particle_center_z_um(surface, 4.0, 1.0)[0, 0] == pytest.approx(1.0)
    with pytest.raises(InputContractError):
        particle_center_z_um(surface, 4.0, 5.0)
    with pytest.raises(InputContractError):
        particle_center_z_um(surface, 0.0, 0.0)


@pytest.mark.parametrize(
    ("diameter", "expected"),
    [(0.5, {"M2": False, "M3": False}), (2.0, {"M2": True, "M3": False}),
     (4.0, {"M2": True, "M3": True})],
)
def test_protrusion_eligibility(pilot, diameter: float, expected: dict[str, bool]) -> None:
    manifest, _, stack_z = pilot
    assert post_step_eligibility(manifest, "S_M1", stack_z, diameter) == expected


def test_exclude_and_wrap_policies(pilot) -> None:
    manifest, _, stack_z = pilot
    excluded = manifest.model_copy(update={"post_step_policy": manifest.post_step_policy
                                           .model_copy(update={"mode": "exclude"})})
    assert post_step_eligibility(excluded, "S_M1", stack_z, 10.0) == {"M2": False, "M3": False}
    wrapped = manifest.model_copy(update={"post_step_policy": manifest.post_step_policy
                                          .model_copy(update={"mode": "wrap_contact"})})
    assert post_step_eligibility(wrapped, "S_M1", stack_z, 0.5) == {"M2": True, "M3": True}


def test_contact_3d_boundary_and_diagonal() -> None:
    distance_xy = np.array([[0.0, 3.0]])
    center_z = np.array([[5.0, 5.0]])
    result = contact_3d(distance_xy, center_z, 0.0, 1.0, 5.0, 0.0)
    # left pixel: Dz = 4 -> Q = 4 <= 5; right pixel: hypot(3, 4) == 5 -> exactly on boundary
    assert result.contact.tolist() == [[True, True]]
    assert result.margin_um[0, 1] == pytest.approx(0.0)
    far = contact_3d(distance_xy, np.array([[20.0, 20.0]]), 0.0, 1.0, 5.0, 0.0)
    assert not far.contact.any()


def test_contact_post_step_requires_eligibility() -> None:
    distance_xy = np.array([[0.0, 9.0]])
    eligible = contact_post_step(distance_xy, True, 5.0, 0.0)
    assert eligible.contact.tolist() == [[True, False]]
    assert eligible.kind == "post_step"
    blocked = contact_post_step(distance_xy, False, 5.0, 0.0)
    assert not blocked.contact.any()
    assert blocked.kind == "excluded"
    assert np.isneginf(blocked.margin_um).all()


def test_zero_dz_matches_edt_threshold(pilot) -> None:
    manifest, masks, stack_z = pilot
    diameter = 6.0
    # embedding = radius keeps the centre exactly on the landing surface (flat model)
    tuned = manifest.model_copy(update={"embedding_depth_um": diameter / 2.0})
    cache = DistanceCache(masks)
    contacts = compute_electrode_contacts(
        tuned, masks, stack_z, "S_PLN", diameter, cache, ["P"]
    )
    expected = distance_to_mask_um(masks.electrodes["P"], masks.grid) <= diameter / 2.0
    assert np.array_equal(contacts["P"].contact, expected)


def test_distance_cache_computes_once(pilot) -> None:
    manifest, masks, stack_z = pilot
    cache = DistanceCache(masks)
    for diameter in (2.0, 4.0):
        compute_electrode_contacts(manifest, masks, stack_z, "S_M1", diameter, cache, ["G", "A"])
    assert cache.computed == 2
    with pytest.raises(InputContractError):
        cache.get("NOPE")


def test_unknown_step(pilot) -> None:
    manifest, masks, stack_z = pilot
    with pytest.raises(InputContractError):
        surface_height_um(manifest, "S_NOPE", stack_z, masks)
    with pytest.raises(InputContractError):
        compute_electrode_contacts(
            manifest, masks, stack_z, "S_NOPE", 2.0, DistanceCache(masks), ["A"]
        )


def test_written_case_paths_are_relative(tmp_path: Path) -> None:
    manifest_path = write_case(build_case("parallel_lines"), tmp_path)
    manifest = load_manifest(manifest_path)
    assert all(m.path.startswith("masks/") for m in manifest.layer_masks)
    assert (tmp_path / "masks" / "electrode_A.png").is_file()
