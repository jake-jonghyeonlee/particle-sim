"""Korean Markdown report with review figures for a nominal run."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from particle_sim.config import load_manifest
from particle_sim.errors import InputContractError
from particle_sim.io.masks import load_product_masks
from particle_sim.io.outputs import check_run_manifest, open_hotspot_store
from particle_sim.io.stack import build_stack_z
from particle_sim.visualization.hotspots import (
    plot_ca_curves,
    plot_layout_overlay,
    plot_pair_panels,
    plot_union_vs_sum,
)

NOTICE = "결과값은 기하 접촉 위험(geometric contact risk)이며 실제 불량확률이 아님"


def _table(headers: list[str], rows: list[list[object]]) -> list[str]:
    def cell(value: object) -> str:
        if isinstance(value, float):
            return f"{value:.6g}"
        return str(value).replace("|", "\\|")

    lines = ["| " + " | ".join(headers) + " |", "|" + "---|" * len(headers)]
    lines.extend("| " + " | ".join(cell(v) for v in row) + " |" for row in rows)
    return lines


def _prepare_report_dir(report_dir: Path, product_id: str) -> None:
    if not report_dir.exists():
        report_dir.mkdir(parents=True)
        return
    foreign = [
        p.name
        for p in report_dir.iterdir()
        if not (p.is_file() and (p.suffix == ".png" or p.name == "report.md"))
    ]
    if foreign:
        raise InputContractError(
            f"report 디렉터리에 생성물이 아닌 항목이 있어 교체하지 않습니다: {sorted(foreign)}",
            product_id=product_id,
            file=str(report_dir),
            field="report",
        )
    shutil.rmtree(report_dir)
    report_dir.mkdir()


def make_report(run_dir: Path, manifest_path: Path) -> Path:
    run_dir = Path(run_dir)
    manifest_path = Path(manifest_path)
    manifest = load_manifest(manifest_path)
    product_id = manifest.product_id
    try:
        run_manifest = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
        ca = pd.read_parquet(run_dir / "ca.parquet")
        union = pd.read_parquet(run_dir / "ca_union.parquet")
        legacy_cmp = pd.read_csv(run_dir / "legacy_vs_edt.csv")
        store = open_hotspot_store(run_dir / "hotspots.zarr")
    except (OSError, ValueError) as exc:
        raise InputContractError(
            f"실행 결과를 읽을 수 없습니다: {exc}",
            product_id=product_id,
            file=str(run_dir),
            field="run",
        ) from exc
    missing = check_run_manifest(run_manifest)
    if missing:
        raise InputContractError(
            f"run manifest에 필수 필드가 없습니다: {missing}",
            product_id=product_id,
            file=str(run_dir / "run_manifest.json"),
            field=missing[0],
        )
    missing_attrs = [a for a in ("source_steps", "particle_sizes_um", "pair_ids")
                     if a not in store.attrs]
    if missing_attrs:
        raise InputContractError(
            f"hotspots.zarr에 필수 항목이 없습니다: {missing_attrs}",
            product_id=product_id,
            file=str(run_dir / "hotspots.zarr"),
            field=missing_attrs[0],
        )
    masks = load_product_masks(manifest, manifest_path)
    stack_z = build_stack_z(manifest)
    report_dir = run_dir / "report"
    _prepare_report_dir(report_dir, product_id)

    step_ids = list(store.attrs["source_steps"])
    sizes = [float(s) for s in store.attrs["particle_sizes_um"]]
    pair_ids = list(store.attrs["pair_ids"])
    largest_index = int(np.argmax(sizes))
    largest = sizes[largest_index]
    figures: list[tuple[str, str]] = []

    patterned = [name for name, layer in stack_z.items() if layer.patterned]
    plot_layout_overlay(masks, patterned, report_dir / "layout_overlay.png")
    figures.append(("Layout overlay", "layout_overlay.png"))
    for step_index, step_id in enumerate(step_ids):
        name = f"ca_curve_{step_id}.png"
        plot_ca_curves(ca, step_id, report_dir / name)
        figures.append((f"{step_id} 크기별 CA", name))
        name = f"union_vs_sum_{step_id}.png"
        plot_union_vs_sum(union, step_id, report_dir / name)
        figures.append((f"{step_id} 합집합 CA vs pair 단순합", name))
        for pair_index, pair_id in enumerate(pair_ids):
            name = f"pair_{step_id}_{pair_id}_d{largest:g}.png"
            plot_pair_panels(
                store["legacy_failure"][step_index, largest_index, pair_index],
                store["nominal_failure"][step_index, largest_index, pair_index],
                store["pair_margin_um"][step_index, largest_index, pair_index],
                masks.grid,
                f"{step_id} / {pair_id} / d={largest:g} um",
                report_dir / name,
            )
            figures.append((f"{step_id} {pair_id} d={largest:g} um", name))

    lines = [
        f"# particle-sim nominal 결과 보고서: {product_id}",
        "",
        f"> {NOTICE}. 산포(Monte Carlo/Sobol)는 반영되지 않은 명목 조건 결과입니다.",
        "",
        f"- run_id: `{run_manifest['run_id']}`",
        f"- simulator_version: `{run_manifest['simulator_version']}`",
        f"- config_hash: `{run_manifest['config_hash']}`",
        f"- assumption_set: `{run_manifest['assumption_set']}`",
        f"- 합성 데이터 여부: `{manifest.is_synthetic}`",
        "",
        "## TODO_PHYSICS / 가정",
        "",
    ]
    lines.extend(f"- {warning}" for warning in run_manifest.get("warnings", []))
    lines += [
        "",
        "## 기존 방식(Erode/Dilate)과 EDT 비교",
        "",
        "기존 방식과 EDT 기준 hotspot의 차이 pixel은 kernel_quantization(커널 반복 단위로 반경이 "
        "양자화됨) 또는 kernel_shape(커널 모양과 Euclidean 거리 차이)로 분류합니다.",
        "",
    ]
    lines += _table(
        list(legacy_cmp.columns),
        [list(row) for row in legacy_cmp.itertuples(index=False)],
    )
    lines += ["", f"## 전극쌍별 CA (d={largest:g} um)", ""]
    ca_rows = ca[np.isclose(ca["particle_size_um"].astype(float), largest)]
    lines += _table(
        ["source_step", "pair_id", "failure_mode", "classification", "ca_legacy_um2",
         "ca_nominal_um2", "ca_nominal_normalized"],
        [
            [r.source_step, r.pair_id, r.failure_mode, r.classification,
             float(r.ca_legacy_um2), float(r.ca_nominal_um2), float(r.ca_nominal_normalized)]
            for r in ca_rows.itertuples(index=False)
        ],
    )
    lines += ["", f"## 합집합 CA vs pair 단순합 (d={largest:g} um)", ""]
    union_rows = union[np.isclose(union["particle_size_um"].astype(float), largest)]
    lines += _table(
        ["source_step", "failure_mode", "pair_ids", "ca_union_nominal_um2",
         "ca_pair_sum_nominal_um2"],
        [
            [r.source_step, r.failure_mode, r.pair_ids, float(r.ca_union_nominal_um2),
             float(r.ca_pair_sum_nominal_um2)]
            for r in union_rows.itertuples(index=False)
        ],
    )
    lines += ["", "## 그림", ""]
    for title, name in figures:
        lines += [f"### {title}", "", f"![{title}]({name})", ""]
    report_path = report_dir / "report.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")
    return report_path
