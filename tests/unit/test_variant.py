from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
import yaml

from particle_sim.cli import main
from particle_sim.config import load_manifest
from particle_sim.errors import InputContractError
from particle_sim.io.masks import GridSpec, load_product_masks
from particle_sim.validation.synthetic_cases import build_case, write_case
from particle_sim.variant.edits import apply_edit, roi_mask
from particle_sim.variant.package import load_edit_spec, make_variants

GRID = GridSpec(shape_yx=(40, 40), pixel_size_y_um=0.5, pixel_size_x_um=0.5)


def _spec(tmp_path: Path, candidates: list[dict]) -> Path:
    path = tmp_path / "edits.yaml"
    path.write_text(yaml.safe_dump({"schema_version": "1.0", "length_unit": "um",
                                    "candidates": candidates}, allow_unicode=True), encoding="utf-8")
    return path


def test_roi_uses_pixel_centres_in_um() -> None:
    roi = roi_mask(GRID, (1.0, 2.0, 3.0, 4.0))
    rows, cols = np.nonzero(roi)
    # 픽셀 중심 (i + 0.5) * 0.5 um가 ROI 안: x 1.0..3.0 -> 열 2..5, y 2.0..4.0 -> 행 4..7
    assert set(cols) == {2, 3, 4, 5} and set(rows) == {4, 5, 6, 7}
    with pytest.raises(ValueError):
        roi_mask(GRID, (3.0, 0.0, 1.0, 1.0))


def test_trim_near_widens_only_the_local_gap() -> None:
    a = np.zeros((40, 40), bool)
    a[:, 5:9] = True
    b = np.zeros((40, 40), bool)
    b[:, 12:20] = True  # 중심 거리 4 px = 2.0 um
    roi = roi_mask(GRID, (0.0, 5.0, 20.0, 10.0))  # 10..19행
    result = apply_edit("trim_near", b, GRID, roi, target=a, amount_um=1.0)
    assert result.min_distance_before_um == pytest.approx(2.0)
    assert result.min_distance_after_um == pytest.approx(3.0)
    assert result.removed[10:20, 12:14].all() and result.removed.sum() == 20
    assert result.electrode[:10, 12].all() and result.electrode[20:, 12].all()
    assert b[10, 12]  # 원본은 바꾸지 않는다


def test_remove_overlap_and_errors() -> None:
    a = np.zeros((40, 40), bool)
    a[10:14, :] = True
    b = np.zeros((40, 40), bool)
    b[:, 20:24] = True
    roi = roi_mask(GRID, (0.0, 0.0, 20.0, 20.0))
    result = apply_edit("remove_overlap", b, GRID, roi, target=a)
    assert result.removed.sum() == 16 and result.overlap_px_after == 0
    with pytest.raises(ValueError, match="ROI 안에"):
        apply_edit("remove_region", b, GRID, roi_mask(GRID, (0.0, 0.0, 1.0, 1.0)))
    with pytest.raises(ValueError, match="상대 전극"):
        apply_edit("trim_near", b, GRID, roi, target=None, amount_um=1.0)
    with pytest.raises(ValueError, match="0보다"):
        apply_edit("trim_near", b, GRID, roi, target=a, amount_um=0.0)


def test_make_variant_writes_valid_package(tmp_path: Path) -> None:
    base = write_case(build_case("parallel_lines"), tmp_path / "base")
    spec = _spec(tmp_path, [{
        "candidate_id": "PL-e1", "change_description": "B 국소 후퇴 1 um",
        "edits": [{"op": "trim_near", "electrode": "B", "target": "A",
                   "roi_um": [40, 40, 60, 60], "amount_um": 1.0}],
    }])
    record = make_variants(base, spec, tmp_path / "out")[0]
    manifest_path = Path(record["manifest"])
    manifest = load_manifest(manifest_path)
    assert manifest.design_version == "V01-e1" and record["parent_design_id"] == "SYN_PARALLEL_LINES-V01"
    masks = load_product_masks(manifest, manifest_path)
    original = load_product_masks(load_manifest(base), base)
    removed = original.electrodes["B"] & ~masks.electrodes["B"]
    assert removed.sum() == 80  # 40행 x 2열(1 um)
    assert not (masks.layers["M2"] & removed).any()  # layer 마스크에서도 지운다
    assert np.array_equal(masks.electrodes["A"], original.electrodes["A"])
    edit = record["edits"][0]
    assert edit["removed_area_um2"] == pytest.approx(20.0)
    assert edit["min_center_distance_in_roi_um"] == {"before": 3.5, "after": 4.5}
    assert record["checks"]["electrodes"]["B"]["split"] is False
    assert record["decision"] == "designer_judgement"
    saved = json.loads((manifest_path.parent / "variant_record.json").read_text(encoding="utf-8"))
    assert saved["design_version"] == "V01-e1"
    assert (manifest_path.parent / "edit_preview.png").exists()
    with pytest.raises(InputContractError, match="비어 있지 않습니다"):
        make_variants(base, spec, tmp_path / "out")


def test_slots_resolve_from_simcompare_sample_id(tmp_path: Path) -> None:
    base = write_case(build_case("vertical_aligned"), tmp_path / "base")
    spec = _spec(tmp_path, [{
        "candidate_id": "VA-e1", "sample_id": "SYN_VERTICAL_ALIGNED-V01-S_M2-G__P",
        "change_description": "P 겹침 일부 제거",
        "edits": [{"op": "remove_overlap", "slot": "b", "target_slot": "a",
                   "roi_um": [40, 40, 60, 60]}],
    }])
    record = make_variants(base, spec, tmp_path / "out")[0]
    edit = record["edits"][0]
    assert edit["electrode"] == "P" and edit["target"] == "G"
    assert edit["overlap_area_um2"]["after"] < edit["overlap_area_um2"]["before"]
    # 선 전체 폭의 겹침을 지우면 P가 둘로 나뉜다: 판정은 하지 않고 표시만 한다.
    assert record["checks"]["electrodes"]["P"]["split"] is True
    bad = _spec(tmp_path, [{
        "candidate_id": "VA-e2", "sample_id": "SYN_VERTICAL_ALIGNED-V09-S_M2-G__P",
        "change_description": "x",
        "edits": [{"op": "remove_region", "slot": "a", "roi_um": [40, 40, 60, 60]}],
    }])
    with pytest.raises(InputContractError, match="찾을 수 없습니다"):
        make_variants(base, bad, tmp_path / "out2")


def test_open_path_break_is_rejected_but_detour_is_kept(tmp_path: Path) -> None:
    cut = [{"op": "remove_region", "electrode": "W", "roi_um": [45, 45, 55, 55]}]
    straight = write_case(build_case("open_straight"), tmp_path / "straight")
    spec = _spec(tmp_path, [{"candidate_id": "OS-e1", "change_description": "cut", "edits": cut}])
    with pytest.raises(InputContractError, match="open 경로 단자 연결"):
        make_variants(straight, spec, tmp_path / "out")
    detour = write_case(build_case("open_detour"), tmp_path / "detour")
    record = make_variants(detour, spec, tmp_path / "out_detour")[0]
    assert record["checks"]["open_paths"] == "ok"


def test_spec_rejects_unknown_keys_and_cli_runs(tmp_path: Path, capsys) -> None:
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"schema_version": "1.0", "length_unit": "um", "candidates": [{
        "candidate_id": "x", "change_description": "x", "typo": 1,
        "edits": [{"op": "remove_region", "electrode": "B", "roi_um": [0, 0, 1, 1]}]}]}))
    with pytest.raises(InputContractError, match="검증 실패"):
        load_edit_spec(bad)
    base = write_case(build_case("parallel_lines"), tmp_path / "base")
    spec = _spec(tmp_path, [
        {"candidate_id": "PL-e1", "change_description": "a",
         "edits": [{"op": "remove_region", "electrode": "B", "roi_um": [47, 90, 54, 100]}]},
        {"candidate_id": "PL-e2", "change_description": "b",
         "edits": [{"op": "remove_region", "electrode": "B", "roi_um": [47, 0, 54, 10]}]},
    ])
    assert main(["make-variant", "--product-config", str(base), "--edits", str(spec),
                 "--output", str(tmp_path / "cli"), "--candidate", "PL-e2"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert [v["candidate_id"] for v in payload["variants"]] == ["PL-e2"]
    assert main(["make-variant", "--product-config", str(base), "--edits", str(spec),
                 "--output", str(tmp_path / "cli2"), "--candidate", "nope"]) == 2


def test_export_carries_variant_parent(tmp_path: Path) -> None:
    import csv

    base = write_case(build_case("parallel_lines"), tmp_path / "base")
    spec = _spec(tmp_path, [{
        "candidate_id": "PL-e1", "change_description": "B 국소 후퇴 1 um",
        "edits": [{"op": "trim_near", "electrode": "B", "target": "A",
                   "roi_um": [40, 40, 60, 60], "amount_um": 1.0}],
    }])
    manifest = Path(make_variants(base, spec, tmp_path / "var")[0]["manifest"])
    assert main(["simulate-nominal", "--product-config", str(manifest),
                 "--output", str(tmp_path / "run")]) == 0
    assert main(["export-simcompare", "--run", str(tmp_path / "run"), "--product-config",
                 str(manifest), "--output", str(tmp_path / "exp")]) == 0
    with (tmp_path / "exp" / "bundles.csv").open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert {r["parent_design_id"] for r in rows} == {"SYN_PARALLEL_LINES-V01"}
    assert {r["change_description"] for r in rows} == {"B 국소 후퇴 1 um"}
    assert all(r["sample_id"].startswith("SYN_PARALLEL_LINES-V01-e1-") for r in rows)

    record_path = manifest.parent / "variant_record.json"
    record = json.loads(record_path.read_text(encoding="utf-8"))
    record_path.write_text(json.dumps({**record, "design_version": "other"}), encoding="utf-8")
    assert main(["export-simcompare", "--run", str(tmp_path / "run"), "--product-config",
                 str(manifest), "--output", str(tmp_path / "exp2")]) == 2
