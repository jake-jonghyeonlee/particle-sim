"""particle-sim 명령줄 도구."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from particle_sim.config import load_manifest
from particle_sim.errors import InputContractError
from particle_sim.io.masks import load_product_masks
from particle_sim.io.stack import build_stack_z
from particle_sim.schema import todo_physics_warnings


def _print_json(payload: dict) -> None:
    print(json.dumps(payload, indent=2, ensure_ascii=False))


def _validate_input(args: argparse.Namespace) -> int:
    manifest_path = Path(args.product_config)
    manifest = load_manifest(manifest_path)
    masks = load_product_masks(manifest, manifest_path)
    stack_z = build_stack_z(manifest)
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
            "warnings": todo_physics_warnings(manifest),
        }
    )
    return 0


def _make_synthetic(args: argparse.Namespace) -> int:
    from particle_sim.validation.synthetic_cases import build_case, write_case

    manifest_path = write_case(build_case(args.case), Path(args.output))
    _print_json({"status": "ok", "case": args.case, "manifest": str(manifest_path)})
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

    export = subparsers.add_parser(
        "export-simcompare", help="실행 결과를 simcompare 입력 형식으로 내보낸다"
    )
    export.add_argument("--run", required=True, help="simulate-nominal 결과 디렉터리")
    export.add_argument("--product-config", required=True, help="제품 manifest YAML 경로")
    export.add_argument("--output", required=True, help="출력 디렉터리(비어 있어야 함)")
    export.set_defaults(func=_export_simcompare)

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
