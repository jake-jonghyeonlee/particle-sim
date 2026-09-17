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


def _apply_fab(args: argparse.Namespace) -> int:
    import pandas as pd

    from particle_sim.aggregation.fab_weighting import apply_fab_weighting
    from particle_sim.io.fab import load_fab_table

    ca = pd.read_parquet(Path(args.ca))
    table = load_fab_table(Path(args.fab_data), bin_method=args.bin_method)
    manifest = load_manifest(Path(args.product_config))
    modes = ["common", "actual"] if args.mode == "both" else [args.mode]
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
    probability.add_argument("--output", required=True, help="실행 결과 디렉터리(비어 있어야 함)")
    probability.add_argument(
        "--convergence-tolerance", type=float, default=None,
        help="수렴 판정 허용 오차(예: 0.01). 지정하면 만족하지 못한 계열을 경고로 남긴다",
    )
    probability.set_defaults(func=_simulate_probability)

    fab = subparsers.add_parser("apply-fab", help="CA에 Fab 이물 분포를 곱해 기대 위험을 만든다")
    fab.add_argument("--ca", required=True, help="simulate-* 결과의 ca.parquet 경로")
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
                     help="연결되지 않은 크기 bin이 있어도 계산을 계속한다")
    fab.add_argument("--output", required=True, help="가중 결과 parquet 경로")
    fab.set_defaults(func=_apply_fab)

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
