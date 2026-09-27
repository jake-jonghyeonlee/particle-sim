from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from particle_sim.cli import main
from particle_sim.errors import InputContractError
from particle_sim.export.simcompare import (
    BUNDLE_COLUMNS,
    CA_RECORD_COLUMNS,
    HOTSPOT_COLUMNS,
    LAYER_STACK_COLUMNS,
    LAYOUT_COLUMNS,
    export_simcompare,
)
from particle_sim.failure.auto_open import INSUFFICIENT_TERMINALS
from particle_sim.failure.open_failure import BASELINE_OPEN, TERMINAL_LOST
from particle_sim.pipeline import run_nominal
from particle_sim.validation.synthetic_cases import build_case, write_case


def _rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


@pytest.fixture(scope="module")
def exported(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path, Path]:
    base = tmp_path_factory.mktemp("export")
    manifest_path = write_case(build_case("pilot"), base / "case")
    run_nominal(manifest_path, base / "run")
    export_simcompare(base / "run", manifest_path, base / "simcompare")
    return manifest_path, base / "run", base / "simcompare"


def test_headers_and_counts(exported: tuple[Path, Path, Path]) -> None:
    _, _, out = exported
    for name, columns, count in (
        ("bundles.csv", BUNDLE_COLUMNS, 9),
        ("layouts.csv", LAYOUT_COLUMNS, 18),
        ("hotspots.csv", HOTSPOT_COLUMNS, 45),
        ("ca_records.csv", CA_RECORD_COLUMNS, 45),
        ("layer_stack.csv", LAYER_STACK_COLUMNS, 5),
    ):
        header, rows = _rows(out / name)
        assert tuple(header) == columns
        assert len(rows) == count


def test_arrays_are_2d_uint8(exported: tuple[Path, Path, Path]) -> None:
    _, _, out = exported
    for name in ("layouts.csv", "hotspots.csv"):
        _, rows = _rows(out / name)
        for row in rows[:6]:
            arr = np.load(out / row["path"], allow_pickle=False)
            assert arr.dtype == np.uint8
            assert arr.shape == (200, 200)
            assert set(np.unique(arr).tolist()) <= {0.0, 1.0}


def test_ca_value_is_normalized_ca(exported: tuple[Path, Path, Path]) -> None:
    _, run_dir, out = exported
    ca = pd.read_parquet(run_dir / "ca.parquet")
    _, rows = _rows(out / "ca_records.csv")
    row = next(r for r in rows if r["sample_id"].endswith("S_M1-A__B")
               and float(r["particle_size"]) == 10.0)
    expected = ca[(ca.source_step == "S_M1") & (ca.pair_id == "A__B")
                  & (ca.particle_size_um == 10.0)]["ca_nominal_normalized"].iloc[0]
    assert float(row["ca_value"]) == pytest.approx(expected)
    assert row["ca_unit"] == "1"


def test_vertical_bundles_leave_fab_layer_empty(exported: tuple[Path, Path, Path]) -> None:
    _, _, out = exported
    _, bundles = _rows(out / "bundles.csv")
    warnings = json.loads((out / "export_warnings.json").read_text(encoding="utf-8"))
    assert any("fab.json" in w for w in warnings)
    for bundle in bundles:
        if bundle["failure_type"] == "vertical_short":
            assert bundle["assigned_process_layer"] == ""
            assert any(bundle["sample_id"] in w for w in warnings)
        else:
            assert bundle["assigned_process_layer"] != ""
    provenance = json.loads((out / "provenance.json").read_text(encoding="utf-8"))
    assert provenance["is_synthetic"] is True


def test_non_empty_output_rejected(exported: tuple[Path, Path, Path]) -> None:
    manifest_path, run_dir, out = exported
    with pytest.raises(InputContractError):
        export_simcompare(run_dir, manifest_path, out)


def test_cli_export(exported: tuple[Path, Path, Path], tmp_path: Path,
                    capsys: pytest.CaptureFixture[str]) -> None:
    manifest_path, run_dir, _ = exported
    argv = ["export-simcompare", "--run", str(run_dir), "--product-config",
            str(manifest_path), "--output", str(tmp_path / "out")]
    assert main(argv) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["bundles"] == 9
    assert main(argv) == 2


@pytest.fixture(scope="module")
def mapping_contract() -> dict:
    return json.loads((Path(__file__).parents[1] / "fixtures/particle_sim_mapping_v1.json")
                      .read_text(encoding="utf-8"))


def _generated_run(manifest_path: Path, run_dir: Path, template: Path,
                   *, contact_first: bool = False) -> tuple[list, object]:
    """Tiny deterministic run planes distinguish every step, unit and size without kernels."""
    import zarr

    from particle_sim.config import load_manifest
    from particle_sim.io.masks import load_product_masks
    from particle_sim.units import unit_specs

    manifest = load_manifest(manifest_path)
    masks = load_product_masks(manifest, manifest_path)
    specs = unit_specs(manifest, masks)
    if contact_first:
        specs.sort(key=lambda u: u.unit_kind != "contact")
    run_dir.mkdir()
    payload = json.loads((template / "run_manifest.json").read_text())
    payload["product_id"] = manifest.product_id
    (run_dir / "run_manifest.json").write_text(json.dumps(payload))
    store = zarr.open_group(run_dir / "hotspots.zarr", mode="w")
    steps = [s.step_id for s in manifest.particle_source_steps]
    sizes = [2.0, 4.0]
    ids = [u.unit_id for u in specs]
    store.attrs.update({"source_steps": steps, "particle_sizes_um": sizes,
                        "pair_ids": ids, "failure_unit_ids": ids,
                        "unit_kinds": [u.unit_kind for u in specs],
                        "pair_modes": [u.failure_mode for u in specs]})
    shape = (len(steps), len(sizes), len(specs), *masks.grid.shape_yx)
    data = np.zeros(shape, dtype=np.uint8)
    rows = []
    for si, step in enumerate(steps):
        for ui, unit in enumerate(specs):
            applicable = unit.applicable_steps is None or step in unit.applicable_steps
            for di, size in enumerate(sizes):
                data[si, di, ui, ui, :si + ui + di + 1] = 1
                rows.append({"source_step": step, "failure_unit_id": unit.unit_id,
                             "pair_id": unit.unit_id, "unit_kind": unit.unit_kind,
                             "failure_mode": unit.failure_mode, "particle_size_um": size,
                             "classification": "confirmed" if applicable else
                             "not_applicable_source_step",
                             "ca_nominal_normalized": float(data[si, di, ui].mean())})
    store.create_array("nominal_failure", data=data)
    pd.DataFrame(rows).to_parquet(run_dir / "ca.parquet", index=False)
    return specs, masks


def _open_manifest(tmp_path: Path, *, auto: bool = False, unit_id: str = "W_PATH") -> Path:
    case = build_case("open_straight")
    case.manifest["particle_sizes_um"] = [2.0, 4.0]
    if auto:
        case.manifest["schema_version"] = "2.0"
        case.manifest["open_paths"] = []
        case.manifest["open_applicable_steps_by_layer"] = {"M2": ["S_M2"]}
    else:
        case.manifest["open_paths"][0]["path_id"] = unit_id
    return write_case(case, tmp_path / "case")


@pytest.mark.parametrize("auto", [False, True])
def test_open_mapping_outputs(tmp_path, exported, mapping_contract, auto):
    manifest = _open_manifest(tmp_path, auto=auto)
    specs, masks = _generated_run(manifest, tmp_path / "run", exported[1])
    out = export_simcompare(tmp_path / "run", manifest, tmp_path / "out")
    _, bundles = _rows(out / "bundles.csv")
    _, layouts = _rows(out / "layouts.csv")
    row = next(r for r in mapping_contract["mapping_rows"]
               if r["unit_kind"] == specs[0].unit_kind)
    assert len(bundles) == (1 if auto else 3)
    assert len(layouts) == len(bundles)
    assert {r["slot"] for r in layouts} == set(row["layout_slots"])
    for b in bundles:
        assert b["failure_type"] == row["target"]
        assert b["electrode_b_role"] == row["electrode_b_role"]
        assert b["layer_b"] == row["layer_b"]
    for layout in layouts:
        np.testing.assert_array_equal(np.load(out / layout["path"]), masks.electrodes["W"])
    assert not list((out / "arrays").glob("*-b.npy"))


@pytest.mark.parametrize("vector_index", range(5))
def test_conformance_id_vectors(tmp_path, exported, mapping_contract, vector_index):
    import yaml

    vector = mapping_contract["sample_id"]["vectors"][vector_index]
    product, version, step, unit = vector["raw_tuple"]
    manifest = _open_manifest(tmp_path, unit_id=unit)
    payload = yaml.safe_load(manifest.read_text())
    payload.update(product_id=product, design_version=version)
    payload["particle_source_steps"] = [dict(payload["particle_source_steps"][0], step_id=step)]
    manifest.write_text(yaml.safe_dump(payload))
    _generated_run(manifest, tmp_path / "run", exported[1])
    out = export_simcompare(tmp_path / "run", manifest, tmp_path / "out")
    sample_id = vector["sample_id"]
    assert _rows(out / "bundles.csv")[1][0]["sample_id"] == sample_id
    for name in ("layouts.csv", "hotspots.csv", "ca_records.csv"):
        assert {r["sample_id"] for r in _rows(out / name)[1]} == {sample_id}
    assert list(json.loads((out / "geometry.json").read_text())) == [sample_id]
    provenance = json.loads((out / "provenance.json").read_text())
    assert provenance["sample_id_mapping"] == [vector]
    for row in _rows(out / "ca_records.csv")[1]:
        size_index = [2.0, 4.0].index(float(row["particle_size"]))
        assert float(row["ca_value"]) == pytest.approx((size_index + 1) / 40000)


@pytest.mark.parametrize("vector_index", range(3))
def test_collision_preflight(tmp_path, exported, mapping_contract, vector_index, capsys):
    import yaml

    vector = mapping_contract["sample_id"]["collision_vectors"][vector_index]
    t1, t2 = vector["raw_tuples"]
    case = build_case("pilot")
    case.manifest.update(product_id=t1[0], design_version=t1[1])
    originals = case.manifest["electrode_pairs"]
    case.manifest["electrode_pairs"] = [dict(originals[i], pair_id=u)
                                        for i, u in enumerate(dict.fromkeys([t1[3], t2[3]]))]
    step = case.manifest["particle_source_steps"][0]
    case.manifest["particle_source_steps"] = [dict(step, step_id=s)
                                              for s in dict.fromkeys([t1[2], t2[2]])]
    manifest = write_case(case, tmp_path / "case")
    assert yaml.safe_load(manifest.read_text())["product_id"] == t1[0]
    _generated_run(manifest, tmp_path / "run", exported[1])
    out = tmp_path / "out"
    with pytest.raises(InputContractError) as error:
        export_simcompare(tmp_path / "run", manifest, out)
    message = str(error.value)
    assert vector["sample_id"] in message
    assert repr(tuple(t1)) in message and repr(tuple(t2)) in message
    assert not out.exists()
    assert main(["export-simcompare", "--run", str(tmp_path / "run"),
                 "--product-config", str(manifest), "--output", str(out)]) == 2
    assert "collision" in capsys.readouterr().err
    assert not out.exists()


@pytest.mark.parametrize("all_cno", [False, True])
def test_contact_exclusion_and_original_plane(tmp_path, exported, mapping_contract, all_cno):
    import yaml
    from PIL import Image

    from tests.unit.test_contacts_cno import _contact_case

    extra = {"open_applicable_steps_by_layer": {"M1": [], "M2": []}} if all_cno else {}
    manifest = _contact_case(tmp_path, extra=extra)
    if not all_cno:
        payload = yaml.safe_load(manifest.read_text())
        mask = np.zeros((60, 80), dtype=np.uint8)
        mask[40:45, 10:20] = 255
        Image.fromarray(mask).save(tmp_path / "masks/electrode_X.png")
        layer_path = tmp_path / "masks/layer_M1.png"
        Image.fromarray(np.maximum(np.asarray(Image.open(layer_path)), mask)).save(layer_path)
        payload["electrodes"].append({"electrode_id": "X", "layer_name": "M1",
                                      "net_id": "N2", "electrode_role": "GATE",
                                      "path": "masks/electrode_X.png"})
        payload["electrode_pairs"] = [{"pair_id": "pair:X-PA", "electrode_a": "X",
                                       "electrode_b": "PA", "failure_mode": "lateral_short"}]
        payload["excluded_electrode_pairs"] = [
            {"electrode_a": "X", "electrode_b": other, "reason": "export fixture"}
            for other in ("PB", "W")]
        manifest.write_text(yaml.safe_dump(payload))
    specs, masks = _generated_run(manifest, tmp_path / "run", exported[1], contact_first=True)
    out = export_simcompare(tmp_path / "run", manifest, tmp_path / "out")
    _, bundles = _rows(out / "bundles.csv")
    warnings = json.loads((out / "export_warnings.json").read_text())
    excluded = [u for u in specs if u.unit_kind == "contact"]
    assert next(r for r in mapping_contract["mapping_rows"] if r["unit_kind"] == "contact")["target"] is None
    for unit in excluded:
        assert sum(unit.unit_id in w for w in warnings) == 1
        assert any(unit.unit_id in w and "kind=contact" in w and "mode=contact_not_open" in w
                   and "excluded" in w for w in warnings)
    if all_cno:
        for name in ("bundles.csv", "layouts.csv", "hotspots.csv", "ca_records.csv"):
            assert _rows(out / name)[1] == []
        assert list((out / "arrays").iterdir()) == []
    else:
        assert {b["failure_type"] for b in bundles} == {"lateral_short", "lateral_open"}
        mapping = json.loads((out / "provenance.json").read_text())["sample_id_mapping"]
        for layout in _rows(out / "layouts.csv")[1]:
            raw = next(m["raw_tuple"] for m in mapping if m["sample_id"] == layout["sample_id"])
            unit = next(u for u in specs if u.unit_id == raw[3])
            electrode = unit.electrode_a if layout["slot"] == "a" else unit.electrode_b
            np.testing.assert_array_equal(np.load(out / layout["path"]), masks.electrodes[electrode])
        for row in _rows(out / "hotspots.csv")[1]:
            mapping = next(m for m in json.loads((out / "provenance.json").read_text())
                           ["sample_id_mapping"] if m["sample_id"] == row["sample_id"])
            ui = next(i for i, u in enumerate(specs) if u.unit_id == mapping["raw_tuple"][3])
            expected = np.zeros(masks.grid.shape_yx)
            si = [s["step_id"] for s in yaml.safe_load(manifest.read_text())
                  ["particle_source_steps"]].index(mapping["raw_tuple"][2])
            expected[ui, :si + ui + [2., 4.].index(float(row["particle_size"])) + 1] = 1
            np.testing.assert_array_equal(np.load(out / row["path"]), expected)
            ca_row = next(r for r in _rows(out / "ca_records.csv")[1]
                          if r["sample_id"] == row["sample_id"] and
                          r["particle_size"] == row["particle_size"])
            assert float(ca_row["ca_value"]) == pytest.approx(expected.mean())


@pytest.mark.parametrize("defect", ["unknown", "aliases", "kind", "mode", "mixed", "nan"])
def test_controlled_preflight_errors(tmp_path, exported, defect):
    import zarr

    manifest = _open_manifest(tmp_path)
    _generated_run(manifest, tmp_path / "run", exported[1])
    store = zarr.open_group(tmp_path / "run/hotspots.zarr", mode="a")
    if defect == "unknown":
        store.attrs["failure_unit_ids"] = store.attrs["pair_ids"] = ["missing:unit"]
    elif defect == "aliases":
        store.attrs["pair_ids"] = ["other"]
    elif defect == "kind":
        store.attrs["unit_kinds"] = ["electrode_pair"]
    elif defect == "mode":
        store.attrs["pair_modes"] = ["vertical_short"]
    else:
        ca = pd.read_parquet(tmp_path / "run/ca.parquet")
        ca.loc[0, "classification" if defect == "mixed" else "ca_nominal_normalized"] = (
            "not_applicable_source_step" if defect == "mixed" else float("nan"))
        ca.to_parquet(tmp_path / "run/ca.parquet", index=False)
    with pytest.raises(InputContractError) as error:
        export_simcompare(tmp_path / "run", manifest, tmp_path / "out")
    if defect == "unknown":
        assert "missing:unit" in str(error.value)
    assert not (tmp_path / "out").exists()


def test_schema_one_and_repeat_determinism(tmp_path, exported):
    import zarr

    case = build_case("pilot")
    manifest = write_case(case, tmp_path / "case")
    _generated_run(manifest, tmp_path / "run", exported[1])
    store = zarr.open_group(tmp_path / "run/hotspots.zarr", mode="a")
    for key in ("failure_unit_ids", "unit_kinds"):
        del store.attrs[key]
    ca = pd.read_parquet(tmp_path / "run/ca.parquet")
    ca.drop(columns=["failure_unit_id", "unit_kind", "classification"]).to_parquet(
        tmp_path / "run/ca.parquet", index=False)
    first = export_simcompare(tmp_path / "run", manifest, tmp_path / "first")
    second = export_simcompare(tmp_path / "run", manifest, tmp_path / "second")
    assert len(_rows(first / "bundles.csv")[1]) == 9
    for path in first.rglob("*"):
        if path.is_file():
            assert path.read_bytes() == (second / path.relative_to(first)).read_bytes()


def test_stack_content_and_slot_order(tmp_path, exported):
    expected = [("M1", "0"), ("GI", "1"), ("M2", "2"), ("PLN", "3"), ("M3", "4")]
    rows = _rows(exported[2] / "layer_stack.csv")[1]
    assert [(r["layer_name"], r["stack_order"]) for r in rows] == expected
    assert {(r["product_id"], r["design_version"]) for r in rows} == {("SYN_PILOT", "V01")}
    assert all(r["fab_process_layer"] == r["layer_name"] for r in rows)
    case = build_case("pilot")
    for pair in case.manifest["electrode_pairs"]:
        pair["electrode_a"], pair["electrode_b"] = pair["electrode_b"], pair["electrode_a"]
    manifest = write_case(case, tmp_path / "case")
    _generated_run(manifest, tmp_path / "run", exported[1])
    out = export_simcompare(tmp_path / "run", manifest, tmp_path / "out")
    assert (out / "layer_stack.csv").read_bytes() == (exported[2] / "layer_stack.csv").read_bytes()
    order = {r["layer_name"]: int(r["stack_order"]) for r in rows}
    def upper_layers(path):
        return {b["sample_id"]: max([b["layer_a"], b["layer_b"]], key=order.get)
                for b in _rows(path / "bundles.csv")[1] if b["failure_type"] == "vertical_short"}
    assert upper_layers(out) == upper_layers(exported[2])
    warnings = json.loads((out / "export_warnings.json").read_text())
    for b in _rows(out / "bundles.csv")[1]:
        if b["failure_type"] == "vertical_short":
            assert b["assigned_process_layer"] == ""
            assert any(b["sample_id"] in w for w in warnings)


def test_anisotropic_rejected_before_files(tmp_path, exported):
    import yaml

    manifest = _open_manifest(tmp_path)
    payload = yaml.safe_load(manifest.read_text())
    payload["coordinate_system"]["pixel_size_y_um"] = 2.0
    payload["coordinate_system"]["roi_height_um"] = 400.0
    manifest.write_text(yaml.safe_dump(payload))
    with pytest.raises(InputContractError, match="pixel"):
        export_simcompare(exported[1], manifest, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_pair_mapping_conformance(exported, mapping_contract):
    bundles = _rows(exported[2] / "bundles.csv")[1]
    layouts = _rows(exported[2] / "layouts.csv")[1]
    for row in mapping_contract["mapping_rows"]:
        if row["unit_kind"] != "electrode_pair":
            continue
        matches = [b for b in bundles if b["failure_type"] == row["mode"]]
        assert matches
        for bundle in matches:
            assert bundle["failure_type"] == row["target"]
            assert [layout["slot"] for layout in layouts
                    if layout["sample_id"] == bundle["sample_id"]] == row["layout_slots"]


@pytest.mark.parametrize("vector_index", range(3))
def test_applicability_conformance(tmp_path, exported, mapping_contract, vector_index):
    vector = mapping_contract["applicability"]["vectors"][vector_index]
    manifest = _open_manifest(tmp_path, auto=True)
    _generated_run(manifest, tmp_path / "run", exported[1])
    ca_path = tmp_path / "run/ca.parquet"
    ca = pd.read_parquet(ca_path)
    step = "S_M1" if vector["result"] == "skip" else "S_M2"
    ca.loc[ca.source_step == step, "classification"] = vector["classifications"]
    ca.to_parquet(ca_path, index=False)
    if vector["result"] == "error":
        with pytest.raises(InputContractError, match="mixed classifications"):
            export_simcompare(tmp_path / "run", manifest, tmp_path / "out")
        assert not (tmp_path / "out").exists()
    else:
        out = export_simcompare(tmp_path / "run", manifest, tmp_path / "out")
        exported_steps = [m["raw_tuple"][2] for m in
                          json.loads((out / "provenance.json").read_text())["sample_id_mapping"]]
        assert (step in exported_steps) == (vector["result"] == "export")
        assert all(np.isfinite(float(r["ca_value"])) for r in _rows(out / "ca_records.csv")[1])


def test_failure_unit_ids_without_legacy_alias(tmp_path, exported):
    import zarr

    manifest = _open_manifest(tmp_path, unit_id="open:primary")
    _generated_run(manifest, tmp_path / "run", exported[1])
    store = zarr.open_group(tmp_path / "run/hotspots.zarr", mode="a")
    del store.attrs["pair_ids"]
    ca_path = tmp_path / "run/ca.parquet"
    pd.read_parquet(ca_path).drop(columns="pair_id").to_parquet(ca_path, index=False)
    out = export_simcompare(tmp_path / "run", manifest, tmp_path / "out")
    assert len(_rows(out / "bundles.csv")[1]) == 3
    assert all(b["sample_id"].endswith("open_primary") for b in _rows(out / "bundles.csv")[1])


def test_casefold_array_collision_preflight(tmp_path, exported):
    case = build_case("pilot")
    originals = case.manifest["electrode_pairs"]
    case.manifest["electrode_pairs"] = [
        dict(originals[0], pair_id="X"), dict(originals[1], pair_id="x")]
    manifest = write_case(case, tmp_path / "case")
    _generated_run(manifest, tmp_path / "run", exported[1])
    with pytest.raises(InputContractError, match="array path collision") as error:
        export_simcompare(tmp_path / "run", manifest, tmp_path / "out")
    assert "S_M1-X" in str(error.value) and "S_M1-x" in str(error.value)
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("auto", [False, True])
@pytest.mark.parametrize("classification", [BASELINE_OPEN, TERMINAL_LOST, INSUFFICIENT_TERMINALS])
def test_excluded_open_warning_conformance(tmp_path, exported, mapping_contract, auto,
                                          classification):
    assert mapping_contract["contract_version"] == "1.1"
    policy = mapping_contract["excluded_open_warnings"]
    assert policy["classifications"] == sorted(
        [BASELINE_OPEN, TERMINAL_LOST, INSUFFICIENT_TERMINALS])
    assert "CA 0 reflects an excluded baseline condition, not low risk" in policy["rule"]
    manifest = _open_manifest(tmp_path, auto=auto)
    specs, _ = _generated_run(manifest, tmp_path / "run", exported[1])
    ca_path = tmp_path / "run/ca.parquet"
    ca = pd.read_parquet(ca_path)
    applicable = ca.classification == "confirmed"
    ca.loc[applicable, "classification"] = classification
    ca.loc[applicable, "ca_nominal_normalized"] = 0.0
    ca.to_parquet(ca_path, index=False)
    out = export_simcompare(tmp_path / "run", manifest, tmp_path / "out")
    warnings = json.loads((out / "export_warnings.json").read_text())
    excluded = [w for w in warnings if "excluded baseline condition" in w]
    steps = ca.loc[applicable, "source_step"].unique()
    assert len(excluded) == len(steps)
    for step in steps:
        assert sum(repr(specs[0].unit_id) in w and repr(step) in w and
                   f"classification={classification}" in w and
                   "CA 0 reflects an excluded baseline condition, not low risk" in w
                   for w in excluded) == 1
    assert len(_rows(out / "bundles.csv")[1]) == len(steps)
    assert all(float(r["ca_value"]) == 0.0 for r in _rows(out / "ca_records.csv")[1])


def _adapter_read(path, kind):
    arr = np.asarray(np.load(path, allow_pickle=False), dtype=np.float32)
    return (arr > 0.5).astype(np.float32) if kind == "layout" else arr


def test_shared_layouts_and_exact_adapter_values(exported):
    import zarr

    from particle_sim.config import load_manifest
    from particle_sim.io.masks import load_product_masks
    from particle_sim.units import unit_specs

    manifest_path, run, out = exported
    manifest = load_manifest(manifest_path)
    masks = load_product_masks(manifest, manifest_path)
    specs = {u.unit_id: u for u in unit_specs(manifest, masks)}
    provenance = json.loads((out / "provenance.json").read_text())
    mapping = {m["sample_id"]: m["raw_tuple"] for m in provenance["sample_id_mapping"]}
    paths = {}
    for row in _rows(out / "layouts.csv")[1]:
        unit = specs[mapping[row["sample_id"]][3]]
        electrode = unit.electrode_a if row["slot"] == "a" else unit.electrode_b
        paths.setdefault(electrode, []).append(row["path"])
        np.testing.assert_array_equal(_adapter_read(out / row["path"], "layout"),
                                      masks.electrodes[electrode].astype(np.float32))
    assert all(len(set(values)) == 1 and len(values) > 1 for values in paths.values())
    assert len(list((out / "arrays").glob("layout-*.npy"))) == len(paths)
    store = zarr.open_group(run / "hotspots.zarr", mode="r")
    for row in _rows(out / "hotspots.csv")[1]:
        raw = mapping[row["sample_id"]]
        si = list(store.attrs["source_steps"]).index(raw[2])
        ui = list(store.attrs["pair_ids"]).index(raw[3])
        di = list(store.attrs["particle_sizes_um"]).index(float(row["particle_size"]))
        np.testing.assert_array_equal(_adapter_read(out / row["path"], "hotspot"),
                                      np.asarray(store["nominal_failure"][si, di, ui],
                                                 dtype=np.float32))
    assert provenance["array_storage"]["format"] == "npy"
    assert provenance["array_storage"]["layout_dtype"] == "uint8"


def test_nonbinary_hotspot_preserved(tmp_path, exported):
    import shutil

    import zarr

    manifest, run, _ = exported
    shutil.copytree(run, tmp_path / "run")
    store = zarr.open_group(tmp_path / "run/hotspots.zarr", mode="a")
    data = np.asarray(store["nominal_failure"][:], dtype=np.float32)
    data[0, 0, 0, 0, 0] = 0.125
    del store["nominal_failure"]
    store.create_array("nominal_failure", data=data)
    out = export_simcompare(tmp_path / "run", manifest, tmp_path / "out")
    row = _rows(out / "hotspots.csv")[1][0]
    assert np.load(out / row["path"]).dtype == np.float32
    np.testing.assert_array_equal(_adapter_read(out / row["path"], "hotspot"), data[0, 0, 0])


def test_capacity_preflight_writes_nothing(tmp_path, exported, monkeypatch):
    from types import SimpleNamespace

    monkeypatch.setattr("particle_sim.export.simcompare.shutil.disk_usage",
                        lambda path: SimpleNamespace(free=1))
    with pytest.raises(InputContractError, match=r"estimated \d+ bytes, free 1 bytes"):
        export_simcompare(exported[1], exported[0], tmp_path / "missing" / "out")
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("existing_empty", [False, True])
def test_atomic_cleanup_on_midway_failure(tmp_path, exported, monkeypatch, existing_empty):
    original = np.save
    count = 0

    def fail_midway(*args, **kwargs):
        nonlocal count
        count += 1
        if count == 3:
            raise OSError("injected write failure")
        return original(*args, **kwargs)

    out = tmp_path / "out"
    if existing_empty:
        out.mkdir()
    monkeypatch.setattr("particle_sim.export.simcompare.np.save", fail_midway)
    with pytest.raises(OSError, match="injected write failure"):
        export_simcompare(exported[1], exported[0], out)
    assert list(tmp_path.iterdir()) == ([out] if existing_empty else [])
    if existing_empty:
        assert list(out.iterdir()) == []
    monkeypatch.setattr("particle_sim.export.simcompare.np.save", original)
    export_simcompare(exported[1], exported[0], out)
    assert (out / "bundles.csv").exists()


def test_layout_casefold_collision_preflight(tmp_path, exported):
    import yaml

    case = build_case("pilot")
    manifest = write_case(case, tmp_path / "case")
    payload = yaml.safe_load(manifest.read_text())
    next(e for e in payload["electrodes"] if e["electrode_id"] == "B")["electrode_id"] = "a"
    for pair in payload["electrode_pairs"]:
        for slot in ("electrode_a", "electrode_b"):
            if pair[slot] == "B":
                pair[slot] = "a"
    manifest.write_text(yaml.safe_dump(payload))
    _generated_run(manifest, tmp_path / "run", exported[1])
    with pytest.raises(InputContractError, match="array path collision"):
        export_simcompare(tmp_path / "run", manifest, tmp_path / "out")
    assert not (tmp_path / "out").exists()
