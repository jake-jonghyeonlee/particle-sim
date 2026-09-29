from __future__ import annotations

import json

import pandas as pd
import pytest

from particle_sim.cli import main
from particle_sim.config import load_manifest
from particle_sim.io.stack import build_stack_z
from particle_sim.pipeline import run_nominal
from particle_sim.schema import ProductManifest
from particle_sim.validation.short_feasibility import structural_zero_short_warnings
from particle_sim.validation.synthetic_cases import build_case, write_case


def _manifest(**changes):
    payload = build_case("pilot").manifest
    for step in payload["particle_source_steps"]:
        step["surface_model"] = "flat"
    payload.update(changes)
    return ProductManifest.model_validate(payload)


@pytest.mark.parametrize("policy,counts", [
    ("protrusion", [2, 2]), ("wrap_contact", [2, 2]), ("exclude", [3, 3, 3]),
])
def test_flat_counts_and_post_step_policy(policy, counts):
    manifest = _manifest(post_step_policy={"mode": policy, "coverage_break_factor": 1.5})
    warnings = structural_zero_short_warnings(manifest, build_stack_z(manifest))
    assert len(warnings) == len(counts)
    expected_steps = ["S_M1", "S_M2", "S_PLN"] if policy == "exclude" else ["S_M2", "S_PLN"]
    for warning, step, count in zip(warnings, expected_steps, counts, strict=True):
        assert f"step {step}:" in warning
        assert f"{count} of 3 short pairs" in warning
        assert "zero contact margin and embedding" in warning
    assert warnings == structural_zero_short_warnings(manifest, build_stack_z(manifest))


@pytest.mark.parametrize("surface,margin,embedding", [
    ("conformal", 0, 0), ("partial", 0, 0), ("flat", 0.1, 0), ("flat", 0, 0.1),
])
def test_no_unproven_claim(surface, margin, embedding):
    manifest = _manifest(contact_margin_um=margin, embedding_depth_um=embedding)
    manifest = manifest.model_copy(update={"particle_source_steps": [
        step.model_copy(update={"surface_model": surface}) for step in manifest.particle_source_steps
    ]})
    assert structural_zero_short_warnings(manifest, build_stack_z(manifest)) == []


def test_open_and_contact_not_counted(tmp_path):
    from tests.unit.test_contacts_cno import _contact_case

    for path in (write_case(build_case("open_straight"), tmp_path / "open"),
                 _contact_case(tmp_path / "contact")):
        manifest = load_manifest(path)
        manifest = manifest.model_copy(update={"particle_source_steps": [
            step.model_copy(update={"surface_model": "flat"})
            for step in manifest.particle_source_steps
        ]})
        assert not manifest.electrode_pairs
        assert structural_zero_short_warnings(manifest, build_stack_z(manifest)) == []


@pytest.mark.parametrize("case_name", ["parallel_lines", "pilot"])
@pytest.mark.parametrize("surface,zero", [("flat", True), ("conformal", False)])
def test_nominal_warning_matches_short_ca(tmp_path, capsys, surface, zero, case_name):
    case = build_case(case_name)
    step = next(s for s in case.manifest["particle_source_steps"] if s["step_id"] == "S_M2")
    step["surface_model"] = surface
    case.manifest["particle_source_steps"] = [step]
    case.manifest["particle_sizes_um"] = [10.0]
    manifest_path = write_case(case, tmp_path / "case")
    assert main(["validate-input", "--product-config", str(manifest_path)]) == 0
    validation = json.loads(capsys.readouterr().out)
    warnings = [w for w in validation["warnings"] if w.startswith("STRUCTURAL_ZERO_SHORT:")]
    assert bool(warnings) == zero
    run_nominal(manifest_path, tmp_path / "run")
    run = json.loads((tmp_path / "run/run_manifest.json").read_text())
    assert [w for w in run["warnings"] if w.startswith("STRUCTURAL_ZERO_SHORT:")] == warnings
    ca = pd.read_parquet(tmp_path / "run/ca.parquet")
    ca = ca[ca.pair_id.isin(["A__B", "G__A"])]
    assert bool((ca.ca_nominal_um2 == 0).all()) == zero
    if zero:
        expected = "1 of 1" if case_name == "parallel_lines" else "2 of 3"
        assert f"{expected} short pairs" in warnings[0]
    else:
        assert ca.ca_nominal_um2.max() > 0
