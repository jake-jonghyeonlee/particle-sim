from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from particle_sim.cli import main
from particle_sim.validation.synthetic_cases import build_case, write_case, write_fab_table

SCENARIO = """
schema_version: "1.0"
assumption_set: nominal
evidence: assumed
overlay_um:
  M2: {mean: 0.0, sigma: 0.3}
"""


@pytest.fixture
def case(tmp_path: Path) -> tuple[Path, Path]:
    manifest_path = write_case(build_case("parallel_lines"), tmp_path / "case")
    scenario_path = tmp_path / "scenario.yaml"
    scenario_path.write_text(SCENARIO, encoding="utf-8")
    return manifest_path, scenario_path


def test_simulate_probability_cli(case: tuple[Path, Path], tmp_path: Path,
                                  capsys: pytest.CaptureFixture[str]) -> None:
    manifest_path, scenario_path = case
    argv = [
        "simulate-probability",
        "--product-config", str(manifest_path),
        "--scenario-config", str(scenario_path),
        "--sobol-m", "6",
        "--seed", "11",
        "--output", str(tmp_path / "prob"),
    ]
    assert main(argv) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["scenario_count"] == 64
    assert Path(payload["outputs"]["convergence"]).is_file()
    # never overwrite an existing run
    assert main(argv) == 2
    assert "입력 오류" in capsys.readouterr().err


def test_make_synthetic_writes_fab_table(tmp_path: Path,
                                         capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["make-synthetic", "--case", "pilot", "--output", str(tmp_path / "out")]) == 0
    payload = json.loads(capsys.readouterr().out)
    fab = pd.read_csv(payload["fab_table"])
    assert {"FAB_SYN", "FAB_REF"} <= set(fab["fab_id"])
    assert (fab["density_per_mm2"] == fab["particle_count"] / fab["inspection_area_mm2"]).all()
    assert "합성" in payload["note"]


def test_apply_fab_cli(case: tuple[Path, Path], tmp_path: Path,
                       capsys: pytest.CaptureFixture[str]) -> None:
    manifest_path, scenario_path = case
    run_dir = tmp_path / "prob"
    assert main([
        "simulate-probability",
        "--product-config", str(manifest_path),
        "--scenario-config", str(scenario_path),
        "--sobol-m", "6", "--seed", "3", "--output", str(run_dir),
    ]) == 0
    capsys.readouterr()
    fab_path = write_fab_table(tmp_path / "fab")
    output = run_dir / "fab_risk.parquet"
    assert main([
        "apply-fab",
        "--ca", str(run_dir / "ca.parquet"),
        "--fab-data", str(fab_path),
        "--product-config", str(manifest_path),
        "--bin-method", "piecewise_constant",
        "--mode", "both",
        "--reference-fab", "FAB_REF",
        "--output", str(output),
    ]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["modes"] == ["common", "actual"]
    weighted = pd.read_parquet(output)
    assert set(weighted["weighting_mode"]) == {"common", "actual"}
    assert (weighted["expected_defect_count"] >= 0).all()
    # the two fabs have different distributions, so the totals must differ
    totals = weighted.groupby("weighting_mode")["expected_defect_count"].sum()
    assert totals["common"] != totals["actual"]


def test_apply_fab_requires_bin_method(case: tuple[Path, Path], tmp_path: Path) -> None:
    manifest_path, _ = case
    fab_path = write_fab_table(tmp_path / "fab")
    with pytest.raises(SystemExit):
        main([
            "apply-fab",
            "--ca", str(tmp_path / "missing.parquet"),
            "--fab-data", str(fab_path),
            "--product-config", str(manifest_path),
            "--mode", "common",
            "--output", str(tmp_path / "out.parquet"),
        ])


def test_apply_fab_run_writes_additive_risk_views(
    case: tuple[Path, Path], tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    manifest_path, _ = case
    run_dir = tmp_path / "nominal"
    assert main([
        "simulate-nominal", "--product-config", str(manifest_path),
        "--output", str(run_dir),
    ]) == 0
    capsys.readouterr()
    fab_path = write_fab_table(tmp_path / "fab-run")
    risk_dir = tmp_path / "risk"
    assert main([
        "apply-fab", "--run", str(run_dir), "--fab-data", str(fab_path),
        "--product-config", str(manifest_path), "--bin-method", "exact_point",
        "--mode", "actual", "--output", str(risk_dir),
    ]) == 0
    capsys.readouterr()
    pair = pd.read_parquet(risk_dir / "risk_by_pair.parquet")
    layer = pd.read_parquet(risk_dir / "risk_by_layer.parquet")
    electrode = pd.read_parquet(risk_dir / "risk_by_electrode.parquet")
    source_total = layer.loc[
        layer["layer_scope"] == "particle_source_layer", "geometric_hit_expectation"
    ].sum()
    assert pair["allocated_geometric_hit_expectation"].sum() == pytest.approx(source_total)
    assert electrode["allocated_electrode_contribution"].sum() == pytest.approx(
        source_total
    )
    assert (
        electrode["involvement_geometric_hit_expectation"]
        >= electrode["allocated_electrode_contribution"]
    ).all()
    assert (risk_dir / "risk_by_open_path.parquet").is_file()
    assert (risk_dir / "risk_summary.csv").is_file()
    assert (risk_dir / "risk_manifest.json").is_file()
    assert {
        "fab_input_semantics", "particle_count_value", "inspection_area_mm2",
        "measurement_period", "particle_type", "fab_source_id", "fab_source_hash",
    } <= set(pair.columns)
    summary = pd.read_csv(risk_dir / "risk_summary.csv")
    assert {"fab_id", "weighting_mode"} <= set(summary.columns)
    risk_manifest = json.loads((risk_dir / "risk_manifest.json").read_text(encoding="utf-8"))
    assert risk_manifest["fab_input_semantics"] == "total_count"
    assert len(risk_manifest["fab_source_sha256"]) == 64
    assert len(risk_manifest["product_config_sha256"]) == 64
    assert len(risk_manifest["source_run_manifest_sha256"]) == 64
    assert {item["path"] for item in risk_manifest["artifacts"]} == {
        "risk_by_layer.parquet", "risk_by_pair.parquet", "risk_by_open_path.parquet",
        "risk_by_contact.parquet", "risk_by_electrode.parquet", "risk_summary.csv",
    }


def test_apply_fab_run_rejects_modified_ca_artifact(
    case: tuple[Path, Path], tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    manifest_path, _ = case
    run_dir = tmp_path / "nominal-modified"
    assert main([
        "simulate-nominal", "--product-config", str(manifest_path),
        "--output", str(run_dir),
    ]) == 0
    capsys.readouterr()
    ca_path = run_dir / "ca.parquet"
    ca = pd.read_parquet(ca_path)
    ca.loc[0, "ca_nominal_um2"] += 1.0
    ca.to_parquet(ca_path, index=False)
    fab_path = write_fab_table(tmp_path / "fab-modified")
    assert main([
        "apply-fab", "--run", str(run_dir), "--fab-data", str(fab_path),
        "--product-config", str(manifest_path), "--bin-method", "exact_point",
        "--mode", "actual", "--output", str(tmp_path / "risk-modified"),
    ]) == 2
    assert "run manifest 기록과 다릅니다" in capsys.readouterr().err


def test_risk_tables_reject_unknown_unit_and_missing_any_union(
    case: tuple[Path, Path], tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    from particle_sim.aggregation.risk_tables import build_risk_tables
    from particle_sim.config import load_manifest
    from particle_sim.errors import InputContractError
    from particle_sim.io.fab import load_fab_table

    manifest_path, _ = case
    run_dir = tmp_path / "nominal-risk-contract"
    assert main([
        "simulate-nominal", "--product-config", str(manifest_path),
        "--output", str(run_dir),
    ]) == 0
    capsys.readouterr()
    ca = pd.read_parquet(run_dir / "ca.parquet")
    union = pd.read_parquet(run_dir / "ca_union.parquet")
    table = load_fab_table(write_fab_table(tmp_path / "fab-contract"), bin_method="exact_point")
    manifest = load_manifest(manifest_path)

    bad_ca = ca.copy()
    bad_ca.loc[0, "unit_kind"] = "unknown"
    with pytest.raises(InputContractError, match="unit_kind"):
        build_risk_tables(
            bad_ca, union, table, manifest, fab_id=manifest.fab_id, weighting_mode="actual"
        )

    missing_union = union.drop(union[union["failure_mode"] == "any"].index[0])
    with pytest.raises(InputContractError, match="1:1"):
        build_risk_tables(
            ca, missing_union, table, manifest,
            fab_id=manifest.fab_id, weighting_mode="actual",
        )


def test_apply_fab_run_rejects_other_design_or_manifest(
    case: tuple[Path, Path], tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    import yaml

    manifest_path, _ = case
    run_dir = tmp_path / "nominal-identity"
    assert main([
        "simulate-nominal", "--product-config", str(manifest_path),
        "--output", str(run_dir),
    ]) == 0
    capsys.readouterr()
    fab_path = write_fab_table(tmp_path / "fab-identity")
    data = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    for name, change, message in (
        ("other-design", {"design_version": "V99_OTHER"}, "design_version"),
        ("other-family", {"product_family": "OTHER_FAMILY"}, "input_hashes.manifest"),
    ):
        other = manifest_path.with_name(f"{name}.yaml")
        other.write_text(yaml.safe_dump({**data, **change}), encoding="utf-8")
        assert main([
            "apply-fab", "--run", str(run_dir), "--fab-data", str(fab_path),
            "--product-config", str(other), "--bin-method", "exact_point",
            "--mode", "actual", "--output", str(tmp_path / f"risk-{name}"),
        ]) == 2
        assert message in capsys.readouterr().err
        assert not (tmp_path / f"risk-{name}").exists()
