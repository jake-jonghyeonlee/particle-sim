"""Build the compact benchmark figure used on report page 11."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
# isort: off
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import zarr
# isort: on


REPORT_DIR = Path(__file__).resolve().parent
EVIDENCE_DIR = REPORT_DIR / "evidence"
FIGURE_DIR = REPORT_DIR / "figures"

CASES = (
    ("pilot-200", "Pilot\n200²", "#287A91"),
    ("complex-200", "복합\n200²", "#E8873A"),
    ("pilot-1000", "Pilot\n1000²", "#287A91"),
    ("complex-1000", "복합\n1000²", "#E8873A"),
)


def load_summary(name: str) -> dict[str, float]:
    path = EVIDENCE_DIR / f"benchmark-{name}.json"
    return json.loads(path.read_text(encoding="utf-8"))["summary"]


def label_bars(axis: plt.Axes, bars, suffix: str) -> None:
    for bar in bars:
        value = bar.get_height()
        axis.annotate(
            f"{value:.1f}{suffix}",
            (bar.get_x() + bar.get_width() / 2, value),
            xytext=(0, 5),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
            color="#193247",
        )


def make_ca_figures() -> None:
    ca = pd.read_parquet(REPORT_DIR / "work" / "run" / "ca.parquet")
    selected = ca[ca["source_step"] == "S_M1"]
    fig, axis = plt.subplots(figsize=(7.0, 5.0), constrained_layout=True)
    colors = {"A__B": "#287DB2", "G__A": "#F28A2B", "B__P": "#2D9A48"}
    for pair_id, rows in selected.groupby("pair_id", sort=False):
        rows = rows.sort_values("particle_size_um")
        axis.plot(
            rows["particle_size_um"],
            rows["ca_nominal_um2"],
            marker="o",
            color=colors[pair_id],
            label=f"{pair_id} · 기준 조건 2.5D",
        )
        axis.plot(
            rows["particle_size_um"],
            rows["ca_legacy_um2"],
            linestyle="--",
            color=colors[pair_id],
            label=f"{pair_id} · 기존 방식",
        )
    axis.set(
        xlabel="이물 직경 (µm)",
        ylabel="CA (µm²)",
        title="이물 유입 공정 S_M1의 크기별 CA",
    )
    axis.grid(color="#DCE6EA", linewidth=0.8)
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(fontsize=8)
    fig.savefig(FIGURE_DIR / "nominal_ca_curve_S_M1.png", dpi=180, facecolor="white")
    plt.close(fig)

    whole = pd.read_parquet(REPORT_DIR / "work" / "triple" / "ca_union.parquet")
    rows = whole[(whole["source_step"] == "S_M2") & (whole["failure_mode"] == "any")].sort_values("particle_size_um")
    fig, axis = plt.subplots(figsize=(7.0, 5.0), constrained_layout=True)
    axis.plot(rows["particle_size_um"], rows["ca_union_nominal_um2"], marker="o", linewidth=2, label="중복 제거 전체 CA")
    axis.plot(rows["particle_size_um"], rows["ca_pair_sum_nominal_um2"], linestyle=":", linewidth=2, label="전극 조합 CA 단순합")
    axis.set(
        xlabel="이물 직경 (µm)",
        ylabel="CA (µm²)",
        title="S_M2 · 중복 제거 CA와 단순합",
    )
    axis.grid(color="#DCE6EA", linewidth=0.8)
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(fontsize=8)
    fig.savefig(FIGURE_DIR / "triple_union_vs_sum_S_M2.png", dpi=180, facecolor="white")
    plt.close(fig)


def make_complex_contact_figure() -> None:
    store = zarr.open_group(REPORT_DIR / "work" / "complex-run" / "hotspots.zarr", mode="r")
    step = list(store.attrs["source_steps"]).index("S_M1")
    size = list(store.attrs["particle_sizes_um"]).index(10.0)
    pair = list(store.attrs["pair_ids"]).index("G2__D3")
    index = (step, size, pair)
    legacy = np.asarray(store["legacy_failure"][index], dtype=float)
    baseline = np.asarray(store["nominal_failure"][index], dtype=float)
    margin = np.asarray(store["pair_margin_um"][index], dtype=float)

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.2), constrained_layout=True)
    for axis, data, title in (
        (axes[0], legacy, "기존 방식 hotspot"),
        (axes[1], baseline, "기준 조건 2.5D hotspot"),
    ):
        image = axis.imshow(data, cmap="Greys", vmin=0, vmax=1, extent=(0, 100, 100, 0), interpolation="nearest")
        axis.set(title=title, xlabel="x (µm)", ylabel="y (µm)")
        fig.colorbar(image, ax=axis, fraction=0.046, pad=0.04)
    image = axes[2].imshow(margin, cmap="RdBu", extent=(0, 100, 100, 0), interpolation="nearest")
    axes[2].set(title="전극 조합 접촉 여유 (µm)", xlabel="x (µm)", ylabel="y (µm)")
    fig.colorbar(image, ax=axes[2], fraction=0.046, pad=0.04)
    fig.suptitle("S_M1 / G2__D3 / d=10 µm", fontsize=13)
    fig.savefig(FIGURE_DIR / "complex_pair_S_M1_G2__D3_d10.png", dpi=180, facecolor="white")
    plt.close(fig)


def main() -> None:
    summaries = [load_summary(name) for name, _, _ in CASES]
    labels = [label for _, label, _ in CASES]
    colors = [color for _, _, color in CASES]
    tat = [summary["wall_median_seconds"] for summary in summaries]
    rss = [summary["peak_rss_max_bytes"] / 1e6 for summary in summaries]

    plt.rcParams.update(
        {
            "font.family": ["Apple SD Gothic Neo", "Arial Unicode MS", "DejaVu Sans"],
            "axes.unicode_minus": False,
            "axes.titleweight": "bold",
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.0), constrained_layout=True)

    for axis, values, title, ylabel, suffix in (
        (axes[0], tat, "기준 조건 TAT", "실행 시간 중앙값 (s)", " s"),
        (axes[1], rss, "최대 메모리", "최대 RSS (MB)", " MB"),
    ):
        bars = axis.bar(labels, values, color=colors, width=0.62)
        axis.set_title(title, color="#133148", pad=10)
        axis.set_ylabel(ylabel, color="#304C5D")
        axis.set_axisbelow(True)
        axis.grid(axis="y", color="#DCE6EA", linewidth=0.8)
        axis.spines[["top", "right"]].set_visible(False)
        axis.spines[["left", "bottom"]].set_color("#9BB0BA")
        axis.set_ylim(0, max(values) * 1.22)
        label_bars(axis, bars, suffix)

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_DIR / "complex_performance.png", dpi=180, facecolor="white")
    plt.close(fig)
    make_ca_figures()
    make_complex_contact_figure()


if __name__ == "__main__":
    main()
