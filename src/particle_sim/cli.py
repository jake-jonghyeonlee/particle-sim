"""particle-sim 명령줄 도구."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from particle_sim.config import load_manifest
from particle_sim.errors import InputContractError
from particle_sim.failure.open_failure import check_open_baselines
from particle_sim.io.masks import load_product_masks
from particle_sim.io.stack import build_stack_z
from particle_sim.schema import todo_physics_warnings
from particle_sim.units import unit_specs
from particle_sim.validation.electrode_pairs import audit_electrode_pairs
from particle_sim.validation.short_feasibility import structural_zero_short_warnings


def _print_json(payload: dict) -> None:
    print(json.dumps(payload, indent=2, ensure_ascii=False))


def _validate_input(args: argparse.Namespace) -> int:
    manifest_path = Path(args.product_config)
    manifest = load_manifest(manifest_path)
    masks = load_product_masks(manifest, manifest_path)
    check_open_baselines(manifest, masks)
    stack_z = build_stack_z(manifest)
    pair_audit = audit_electrode_pairs(manifest, masks, stack_z)
    topology = masks.contact_topology
    contacts = {
        "holes": 0 if topology is None else len(topology.holes),
        "links": [] if topology is None else sorted(topology.links),
    }
    total_layer_pixels = {
        name: int(mask.sum()) for name, mask in masks.layers.items()
    }
    unassigned = {
        name: {
            "pixels": int(count),
            "area_um2": float(count * masks.grid.pixel_area_um2),
            "layer_fraction": (
                float(count / total_layer_pixels[name]) if total_layer_pixels[name] else 0.0
            ),
        }
        for name, count in (masks.layer_unassigned_pixels or {}).items()
    }
    _print_json(
        {
            "status": "ok",
            "product_id": manifest.product_id,
            "shape_yx": list(masks.grid.shape_yx),
            "pixel_size_um": [masks.grid.pixel_size_y_um, masks.grid.pixel_size_x_um],
            "layers": list(stack_z),
            "electrodes": len(manifest.electrodes),
            "pairs": len(manifest.electrode_pairs),
            "source_steps": [s.step_id for s in manifest.particle_source_steps],
            "layer_unassigned_regions": unassigned,
            "electrode_pair_audit": [row.__dict__ for row in pair_audit],
            "contacts": contacts,
            "failure_units": [
                {"unit_id": spec.unit_id, "unit_kind": spec.unit_kind}
                for spec in unit_specs(manifest, masks)
            ],
            "unresolved_electrode_pairs": sum(
                row.status == "unresolved" for row in pair_audit
            ),
            "warnings": [*todo_physics_warnings(manifest),
                         *structural_zero_short_warnings(manifest, stack_z)],
        }
    )
    return 0


def _make_synthetic(args: argparse.Namespace) -> int:
    from particle_sim.validation.synthetic_cases import build_case, write_case, write_fab_table

    manifest_path = write_case(build_case(args.case), Path(args.output))
    fab_path = write_fab_table(Path(args.output))
    _print_json(
        {
            "status": "ok",
            "case": args.case,
            "manifest": str(manifest_path),
            "fab_table": str(fab_path),
            "note": "합성 자료이며 실제 측정값이 아닙니다",
        }
    )
    return 0


def _simulate_nominal(args: argparse.Namespace) -> int:
    from particle_sim.pipeline import run_nominal

    result = run_nominal(Path(args.product_config), Path(args.output))
    _print_json(
        {
            "status": "ok",
            "run_id": result.run_id,
            "outputs": {
                "hotspots": str(result.hotspots),
                "ca_table": str(result.ca_table),
                "union_table": str(result.union_table),
                "legacy_comparison": str(result.legacy_comparison),
                "run_manifest": str(result.run_manifest),
            },
        }
    )
    return 0


def _simulate_probability(args: argparse.Namespace) -> int:
    from particle_sim.pipeline import run_probability

    result = run_probability(
        Path(args.product_config),
        Path(args.scenario_config),
        Path(args.output),
        sobol_exponent_m=args.sobol_m,
        seed=args.seed,
        additional_seeds=args.additional_seed,
        convergence_tolerance=args.convergence_tolerance,
    )
    _print_json(
        {
            "status": "ok",
            "run_id": result.run_id,
            "scenario_count": result.scenario_count,
            "outputs": {
                "probability": str(result.probability),
                "ca_table": str(result.ca_table),
                "union_table": str(result.union_table),
                "convergence": str(result.convergence),
                "run_manifest": str(result.run_manifest),
            },
        }
    )
    return 0


def _analyze_central(args: argparse.Namespace) -> int:
    from particle_sim.validation.central_window import analyze_central

    result = analyze_central(
        Path(args.product_config), Path(args.output),
        run_dir=Path(args.run) if args.run else None,
    )
    _print_json({"status": "ok", "output": args.output, **result})
    return 0


def _analyze_cno_open(args: argparse.Namespace) -> int:
    from particle_sim.validation.cno_open_dedup import analyze_cno_open

    result = analyze_cno_open(
        Path(args.product_config), Path(args.adapter_sidecar),
        Path(args.run), Path(args.output),
    )
    _print_json({"status": "ok", "output": args.output, **result})
    return 0


def _apply_fab(args: argparse.Namespace) -> int:
    import pandas as pd

    from particle_sim.aggregation.fab_weighting import apply_fab_weighting
    from particle_sim.aggregation.risk_tables import build_risk_tables
    from particle_sim.config import file_sha256
    from particle_sim.io.fab import load_fab_table
    from particle_sim.io.outputs import (
        SCHEMA_VERSION,
        artifact_digest,
        artifact_records,
        write_json,
    )

    table = load_fab_table(Path(args.fab_data), bin_method=args.bin_method)
    manifest = load_manifest(Path(args.product_config))
    modes = ["common", "actual"] if args.mode == "both" else [args.mode]
    if args.run:
        if args.ca:
            raise InputContractError("--run과 --ca는 함께 쓸 수 없습니다", field="run")
        run_dir = Path(args.run)
        try:
            ca = pd.read_parquet(run_dir / "ca.parquet")
            union = pd.read_parquet(run_dir / "ca_union.parquet")
            source_manifest = json.loads(
                (run_dir / "run_manifest.json").read_text(encoding="utf-8")
            )
        except (OSError, ValueError) as exc:
            raise InputContractError(
                f"실행 결과 CA를 읽을 수 없습니다: {exc}", file=str(run_dir), field="run"
            ) from exc
        if not isinstance(source_manifest, dict):
            raise InputContractError(
                "run_manifest.json은 mapping이어야 합니다", file=str(run_dir), field="run"
            )
        for label, frame in (("ca.parquet", ca), ("ca_union.parquet", union)):
            for column in ("schema_version", "run_id", "product_id"):
                if column not in frame.columns or frame[column].nunique(dropna=False) != 1:
                    raise InputContractError(
                        f"{label}의 {column} 식별자가 하나로 고정되지 않았습니다",
                        field=column,
                    )
            if str(frame["schema_version"].iloc[0]) != SCHEMA_VERSION:
                raise InputContractError(
                    f"{label}는 schema {SCHEMA_VERSION} 결과여야 합니다",
                    field="schema_version",
                )
        run_ids = {str(ca["run_id"].iloc[0]), str(union["run_id"].iloc[0])}
        product_ids = {
            str(ca["product_id"].iloc[0]), str(union["product_id"].iloc[0])
        }
        if run_ids != {str(source_manifest.get("run_id"))}:
            raise InputContractError("CA 표와 run manifest의 run_id가 다릅니다", field="run_id")
        if product_ids != {manifest.product_id} or source_manifest.get(
            "product_id"
        ) != manifest.product_id:
            raise InputContractError(
                "CA 표, run manifest와 제품 YAML의 product_id가 다릅니다",
                field="product_id",
            )
        artifacts = {
            str(item.get("path")): item
            for item in source_manifest.get("artifacts", [])
            if isinstance(item, dict)
        }
        for name in ("ca.parquet", "ca_union.parquet"):
            record = artifacts.get(name)
            if record is None:
                raise InputContractError(
                    f"run manifest에 {name} artifact 기록이 없습니다", field="artifacts"
                )
            digest, size = artifact_digest(run_dir / name)
            if digest != record.get("sha256") or size != record.get("size_bytes"):
                raise InputContractError(
                    f"{name}이 run manifest 기록과 다릅니다", field="artifacts"
                )
        output_dir = Path(args.output)
        if output_dir.exists() and (not output_dir.is_dir() or any(output_dir.iterdir())):
            raise InputContractError(
                f"출력 경로가 비어 있지 않습니다: {output_dir}", field="output"
            )
        output_dir.mkdir(parents=True, exist_ok=True)
        combined: dict[str, list[pd.DataFrame]] = {
            "risk_by_layer.parquet": [],
            "risk_by_pair.parquet": [],
            "risk_by_open_path.parquet": [],
            "risk_by_contact.parquet": [],
            "risk_by_electrode.parquet": [],
            "risk_summary.csv": [],
        }
        for mode in modes:
            fab_id = args.reference_fab if mode == "common" else manifest.fab_id
            tables = build_risk_tables(
                ca,
                union,
                table,
                manifest,
                fab_id=fab_id,
                weighting_mode=mode,
            )
            combined["risk_by_layer.parquet"].append(tables.by_layer)
            combined["risk_by_pair.parquet"].append(tables.by_pair)
            combined["risk_by_open_path.parquet"].append(tables.by_open_path)
            combined["risk_by_contact.parquet"].append(tables.by_contact)
            combined["risk_by_electrode.parquet"].append(tables.by_electrode)
            combined["risk_summary.csv"].append(tables.summary)
        outputs = {}
        for name, frames in combined.items():
            frame = pd.concat(frames, ignore_index=True, sort=False)
            path = output_dir / name
            if path.suffix == ".csv":
                frame.to_csv(path, index=False)
            else:
                frame.to_parquet(path, index=False, engine="pyarrow")
            outputs[name] = str(path)
        artifact_paths = [output_dir / name for name in combined]
        risk_manifest = output_dir / "risk_manifest.json"
        run_id = str(ca["run_id"].iloc[0]) if len(ca) else ""
        write_json(
            {
                "schema_version": SCHEMA_VERSION,
                "run_id": run_id,
                "product_id": manifest.product_id,
                "source_run": str(run_dir),
                "fab_source": table.source,
                "source_run_manifest_sha256": file_sha256(run_dir / "run_manifest.json"),
                "product_config_sha256": file_sha256(Path(args.product_config)),
                "fab_source_sha256": file_sha256(Path(args.fab_data)),
                "source_run_artifacts": source_manifest.get("artifacts", []),
                "fab_input_semantics": table.input_semantics,
                "weighting_modes": modes,
                "units": {
                    "particle_size": "um",
                    "critical_area": "um2",
                    "inspection_area": "mm2",
                    "density": "count/mm2",
                    "geometric_hit_expectation": "expected geometric hits",
                },
                "allocation_semantics": {
                    "failure_unit": "overlap pixel split equally by 1/k active units",
                    "short_electrode": "allocated unit contribution split equally between electrodes",
                    "open_electrode": "all allocated unit contribution assigned to its electrode",
                    "cno_electrode": "allocated CNO contribution split equally between the "
                    "lower and upper electrode",
                },
                "value_semantics": (
                    "critical area multiplied by measured mean particle density; "
                    "geometric hit expectation, not confirmed electrical failure rate"
                ),
                "artifacts": artifact_records(artifact_paths, base=output_dir),
            },
            risk_manifest,
        )
        outputs[risk_manifest.name] = str(risk_manifest)
        _print_json({"status": "ok", "output": str(output_dir), "outputs": outputs})
        return 0
    if not args.ca:
        raise InputContractError("--run 또는 --ca 중 하나가 필요합니다", field="ca")
    ca = pd.read_parquet(Path(args.ca))
    frames = []
    for mode in modes:
        frames.append(
            apply_fab_weighting(
                ca,
                table,
                weighting_mode=mode,
                fab_id=args.reference_fab if mode == "common" else manifest.fab_id,
                product_family=manifest.product_family,
                run_id=str(ca["run_id"].iloc[0]) if len(ca) else "",
                product_id=manifest.product_id,
                require_full_coverage=not args.allow_missing_bins,
            )
        )
    weighted = pd.concat(frames, ignore_index=True)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    weighted.to_parquet(output, index=False, engine="pyarrow")
    _print_json(
        {
            "status": "ok",
            "output": str(output),
            "rows": int(len(weighted)),
            "modes": modes,
            "expected_defect_count": {
                mode: float(
                    weighted[weighted["weighting_mode"] == mode]["expected_defect_count"].sum()
                )
                for mode in modes
            },
        }
    )
    return 0


def _export_simcompare(args: argparse.Namespace) -> int:
    import csv

    from particle_sim.export.simcompare import export_simcompare

    out_dir = export_simcompare(Path(args.run), Path(args.product_config), Path(args.output))
    with (out_dir / "bundles.csv").open(encoding="utf-8") as handle:
        bundles = len(list(csv.DictReader(handle)))
    warnings = json.loads((out_dir / "export_warnings.json").read_text(encoding="utf-8"))
    _print_json(
        {"status": "ok", "output": str(out_dir), "bundles": bundles, "warnings": len(warnings)}
    )
    return 0


def _make_report(args: argparse.Namespace) -> int:
    from particle_sim.visualization.reports import make_report

    report = make_report(Path(args.run), Path(args.product_config))
    _print_json({"status": "ok", "report": str(report)})
    return 0


def _compare_sim(args: argparse.Namespace) -> int:
    from particle_sim.comparison.legacy_sim import compare_legacy_sim

    output = compare_legacy_sim(
        Path(args.reference_config), Path(args.run), Path(args.output)
    )
    manifest = json.loads((output / "comparison_manifest.json").read_text(encoding="utf-8"))
    _print_json({"status": "ok", "output": str(output), **manifest})
    return 0


def build_parser() -> argparse.ArgumentParser:
    """CLI 인자 구조를 만든다."""
    parser = argparse.ArgumentParser(
        prog="particle-sim",
        description="2.5D 이물 기하 접촉 위험 시뮬레이터 (실제 불량확률이 아님)",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate = subparsers.add_parser("validate-input", help="제품 manifest와 마스크를 검증한다")
    validate.add_argument("--product-config", required=True, help="제품 manifest YAML 경로")
    validate.set_defaults(func=_validate_input)

    synthetic = subparsers.add_parser("make-synthetic", help="합성 fixture 제품을 생성한다")
    synthetic.add_argument("--case", required=True, help="합성 case 이름")
    synthetic.add_argument("--output", required=True, help="출력 디렉터리")
    synthetic.set_defaults(func=_make_synthetic)

    nominal = subparsers.add_parser("simulate-nominal", help="명목 조건 hotspot과 CA를 계산한다")
    nominal.add_argument("--product-config", required=True, help="제품 manifest YAML 경로")
    nominal.add_argument("--output", required=True, help="실행 결과 디렉터리(비어 있어야 함)")
    nominal.set_defaults(func=_simulate_nominal)

    probability = subparsers.add_parser(
        "simulate-probability", help="공정 산포를 반영한 확률 hotspot과 CA를 계산한다"
    )
    probability.add_argument("--product-config", required=True, help="제품 manifest YAML 경로")
    probability.add_argument("--scenario-config", required=True, help="산포 시나리오 YAML 경로")
    probability.add_argument("--sobol-m", type=int, required=True, help="Sobol 지수 m (2^m 조건)")
    probability.add_argument("--seed", type=int, required=True, help="난수 seed")
    probability.add_argument(
        "--additional-seed",
        type=int,
        action="append",
        default=[],
        help="추가 Sobol scramble seed. 3회 검증이면 두 번 지정한다",
    )
    probability.add_argument("--output", required=True, help="실행 결과 디렉터리(비어 있어야 함)")
    probability.add_argument(
        "--convergence-tolerance", type=float, default=None,
        help="수렴 판정 허용 오차(예: 0.01). 지정하면 만족하지 못한 계열을 경고로 남긴다",
    )
    probability.set_defaults(func=_simulate_probability)

    central = subparsers.add_parser(
        "analyze-central", help="1200 ROI의 중앙 600x600 통계와 지름별 support 검사"
    )
    central.add_argument("--product-config", required=True, help="제품 manifest YAML 경로")
    central.add_argument("--run", help="기존 simulate 결과 폴더; 생략하면 support만 검사")
    central.add_argument("--output", required=True, help="비어 있는 별도 분석 출력 폴더")
    central.set_defaults(func=_analyze_central)

    dedup = subparsers.add_parser(
        "analyze-cno-open", help="명시 회로 관계에 따른 CNO/open 중복 회계"
    )
    dedup.add_argument("--product-config", required=True)
    dedup.add_argument("--adapter-sidecar", required=True)
    dedup.add_argument("--run", required=True, help="독립적으로 계산한 nominal 원시 결과")
    dedup.add_argument("--output", required=True, help="별도 빈 결과 디렉터리")
    dedup.set_defaults(func=_analyze_cno_open)

    fab = subparsers.add_parser("apply-fab", help="CA에 Fab 이물 분포를 곱해 기대 위험을 만든다")
    fab.add_argument("--ca", help="이전 방식: simulate-* 결과의 ca.parquet 경로")
    fab.add_argument("--run", help="권장 방식: ca.parquet와 ca_union.parquet가 있는 실행 폴더")
    fab.add_argument("--fab-data", required=True, help="Fab 이물 표(parquet/csv) 경로")
    fab.add_argument("--product-config", required=True, help="제품 manifest YAML 경로")
    fab.add_argument(
        "--bin-method", required=True,
        choices=["exact_point", "piecewise_constant", "piecewise_linear", "empirical_samples"],
        help="Fab 크기 bin 연결 방식(기본값 없음)",
    )
    fab.add_argument("--mode", default="both", choices=["common", "actual", "both"],
                     help="공통 기준 Fab, 실제 Fab 또는 둘 다")
    fab.add_argument("--reference-fab", default="", help="common 모드에서 쓸 기준 Fab ID")
    fab.add_argument("--allow-missing-bins", action="store_true",
                     help="이전 --ca 방식에서만 누락 크기 bin을 허용한다")
    fab.add_argument(
        "--output", required=True,
        help="--run이면 비어 있는 출력 폴더, --ca이면 가중 결과 parquet 경로",
    )
    fab.set_defaults(func=_apply_fab)

    export = subparsers.add_parser(
        "export-simcompare", help="실행 결과를 simcompare 입력 형식으로 내보낸다"
    )
    export.add_argument("--run", required=True, help="simulate-nominal 결과 디렉터리")
    export.add_argument("--product-config", required=True, help="제품 manifest YAML 경로")
    export.add_argument("--output", required=True, help="출력 디렉터리(비어 있어야 함)")
    export.set_defaults(func=_export_simcompare)

    compare = subparsers.add_parser(
        "compare-sim", help="기존 sim hotspot과 현재 기준 조건 결과를 픽셀 단위로 비교한다"
    )
    compare.add_argument("--reference-config", required=True, help="기존 sim 배열 매핑 YAML/JSON")
    compare.add_argument("--run", required=True, help="simulate-nominal 결과 디렉터리")
    compare.add_argument("--output", required=True, help="비어 있는 비교 결과 디렉터리")
    compare.set_defaults(func=_compare_sim)

    report = subparsers.add_parser("make-report", help="그림과 한국어 보고서를 생성한다")
    report.add_argument("--run", required=True, help="simulate-nominal 결과 디렉터리")
    report.add_argument("--product-config", required=True, help="제품 manifest YAML 경로")
    report.set_defaults(func=_make_report)

    return parser


def main(argv: list[str] | None = None) -> int:
    """CLI 진입점. 입력 계약 위반은 stderr로 알리고 2를 반환한다."""
    import sys

    args = build_parser().parse_args(argv)
    try:
        return int(args.func(args))
    except InputContractError as err:
        print(f"입력 오류: {err}", file=sys.stderr)
        return 2
    except ValueError as err:
        print(f"입력 오류: {err}", file=sys.stderr)
        return 2
