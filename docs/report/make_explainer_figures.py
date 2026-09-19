"""Build explanatory figures for report pages 3 through 6."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

# isort: off
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import zarr  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402
from PIL import Image  # noqa: E402

from particle_sim.config import load_manifest  # noqa: E402
from particle_sim.failure.open_failure import (  # noqa: E402
    blocked_mask,
    causes_open,
    simulate_open,
)
from particle_sim.io.masks import GridSpec, MaskSet  # noqa: E402
from particle_sim.schema import OpenPath  # noqa: E402
# isort: on


REPORT_DIR = Path(__file__).resolve().parent
FIGURE_DIR = REPORT_DIR / "figures"
CASE_DIR = REPORT_DIR / "evidence" / "case"
WORK_DIR = REPORT_DIR / "work"

plt.rcParams.update(
    {
        "font.family": ["Apple SD Gothic Neo", "Arial Unicode MS", "DejaVu Sans"],
        "axes.unicode_minus": False,
    }
)

ELECTRODES = (
    ("G", "M1", "N1", "#4C9BC4"),
    ("A", "M2", "N2", "#F2A65A"),
    ("B", "M2", "N3", "#E67E3F"),
    ("P", "M3", "N4", "#67B66A"),
)


def electrode_masks() -> dict[str, np.ndarray]:
    return {
        electrode: np.asarray(
            Image.open(CASE_DIR / "masks" / f"electrode_{electrode}.png")
        )
        > 0
        for electrode, _, _, _ in ELECTRODES
    }


def mask_span_um(mask: np.ndarray, pixel_size_um: float = 0.5) -> tuple[float, float]:
    _, xs = np.where(mask)
    return float(xs.min()) * pixel_size_um, float(xs.max() + 1) * pixel_size_um


def make_layout() -> None:
    masks = electrode_masks()
    fig, axis = plt.subplots(figsize=(5.8, 5.0), constrained_layout=True)
    for electrode, layer, net, color in ELECTRODES:
        rgba = np.zeros((*masks[electrode].shape, 4), dtype=float)
        rgba[..., :3] = matplotlib.colors.to_rgb(color)
        rgba[..., 3] = masks[electrode].astype(float) * 0.82
        axis.imshow(rgba, extent=(0, 100, 100, 0), interpolation="nearest")
        left, right = mask_span_um(masks[electrode])
        axis.text(
            (left + right) / 2,
            50,
            f"{electrode}\n{layer} · {net}",
            ha="center",
            va="center",
            rotation=90,
            fontsize=9,
            fontweight="bold",
            color="#173247",
            bbox={"boxstyle": "round,pad=0.25", "fc": "white", "ec": "none", "alpha": 0.82},
        )
    axis.set(xlim=(0, 100), ylim=(100, 0), xlabel="x (µm)", ylabel="y (µm)")
    axis.set_title("Pilot 전극 mask · 평면도", color="#173247", fontweight="bold")
    axis.text(
        98,
        96,
        "각 전극 띠는 y방향 전체로 이어짐",
        ha="right",
        va="bottom",
        fontsize=8,
        color="#526B78",
    )
    axis.grid(color="#DCE6EA", linewidth=0.6, alpha=0.75)
    fig.savefig(FIGURE_DIR / "explainer_layout.png", dpi=180, facecolor="white")
    plt.close(fig)


def make_layout_3d_draft() -> None:
    """Draft page-3 concept figure showing mask height and short direction."""
    masks = electrode_masks()
    z_by_mask = {"mask1": 0.0, "mask2": 0.7, "mask3": 3.0}
    mask_by_electrode = {"G": "mask1", "A": "mask2", "B": "mask2", "P": "mask3"}
    colors = {electrode: color for electrode, _, _, color in ELECTRODES}
    spans = {electrode: mask_span_um(masks[electrode]) for electrode, *_ in ELECTRODES}

    fig = plt.figure(figsize=(11.2, 5.4), constrained_layout=True)
    axis = fig.add_subplot(121, projection="3d")
    for mask_name, z in z_by_mask.items():
        plane = [[(0, 0, z), (100, 0, z), (100, 100, z), (0, 100, z)]]
        axis.add_collection3d(
            Poly3DCollection(
                plane,
                facecolors="#DDE8EC",
                edgecolors="#9AB1BB",
                linewidths=0.8,
                alpha=0.22,
            )
        )

    for electrode, _layer, _net, _color in ELECTRODES:
        left, right = spans[electrode]
        z = z_by_mask[mask_by_electrode[electrode]] + 0.035
        surface = [[(left, 0, z), (right, 0, z), (right, 100, z), (left, 100, z)]]
        axis.add_collection3d(
            Poly3DCollection(
                surface,
                facecolors=colors[electrode],
                edgecolors="#173247",
                linewidths=0.7,
                alpha=0.88,
            )
        )
        axis.text((left + right) / 2, 72, z + 0.11, f"{electrode}", ha="center", fontweight="bold")

    axis.set(
        xlim=(0, 100),
        ylim=(0, 100),
        zlim=(-0.1, 3.6),
        xlabel="x (µm)",
        ylabel="y (µm)",
        title="1. mask의 상하 배치",
    )
    axis.set_zticks([0.0, 0.7, 3.0], labels=["mask1", "mask2", "mask3"])
    axis.view_init(elev=25, azim=-62)
    axis.set_box_aspect((1.45, 1.0, 0.75))

    relation_axis = fig.add_subplot(122)
    height = 0.16
    centers: dict[str, tuple[float, float]] = {}
    for electrode, _layer, _net, _color in ELECTRODES:
        left, right = spans[electrode]
        z = z_by_mask[mask_by_electrode[electrode]]
        relation_axis.add_patch(
            Rectangle((left, z - height / 2), right - left, height, facecolor=colors[electrode], edgecolor="#173247", linewidth=1.0)
        )
        centers[electrode] = ((left + right) / 2, z)
        relation_axis.text((left + right) / 2, z, electrode, ha="center", va="center", fontweight="bold")

    def arrow(a: str, b: str, label: str, color: str, text_offset: tuple[float, float]) -> None:
        start, end = centers[a], centers[b]
        relation_axis.annotate(
            "",
            xy=end,
            xytext=start,
            arrowprops={"arrowstyle": "<->", "linestyle": "--", "linewidth": 2.0, "color": color},
        )
        relation_axis.text(
            (start[0] + end[0]) / 2 + text_offset[0],
            (start[1] + end[1]) / 2 + text_offset[1],
            label,
            color=color,
            ha="center",
            fontweight="bold",
            fontsize=9,
            bbox={"boxstyle": "round,pad=0.18", "fc": "white", "ec": "none", "alpha": 0.9},
        )

    arrow("A", "B", "A__B\nlateral short", "#D36F1D", (0, 0.28))
    arrow("G", "A", "G__A\nvertical short", "#6D4AB5", (-7, 0.05))
    arrow("B", "P", "B__P\nvertical short", "#6D4AB5", (8, 0.10))
    relation_axis.set(
        xlim=(15, 65),
        ylim=(-0.35, 3.35),
        xlabel="x (µm)",
        ylabel="높이 방향",
        title="2. 평가 전극 조합",
    )
    relation_axis.set_yticks([0.0, 0.7, 3.0], labels=["mask1", "mask2", "mask3"])
    relation_axis.grid(color="#DCE6EA", linewidth=0.8)
    relation_axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle("전극 mask의 높이와 lateral short·vertical short", fontsize=15, fontweight="bold", color="#173247")
    fig.savefig(FIGURE_DIR / "explainer_layout_3d_draft.png", dpi=180, facecolor="white")
    plt.close(fig)


def failure_span_um(failure: np.ndarray) -> tuple[float, float]:
    _, xs = np.where(np.asarray(failure, dtype=bool))
    return float(xs.min()) * 0.5, float(xs.max() + 1) * 0.5


def make_contact_explainer() -> None:
    masks = electrode_masks()
    store = zarr.open_group(WORK_DIR / "run" / "hotspots.zarr", mode="r")
    pair_index = list(store.attrs["pair_ids"]).index("G__A")
    size_index = list(store.attrs["particle_sizes_um"]).index(10.0)
    step_index = list(store.attrs["source_steps"]).index("S_M1")
    index = (step_index, size_index, pair_index)
    legacy = np.asarray(store["legacy_failure"][index], dtype=bool)
    nominal = np.asarray(store["nominal_failure"][index], dtype=bool)
    margin = np.asarray(store["pair_margin_um"][index], dtype=float)
    if not np.allclose(margin, margin[0], equal_nan=True):
        raise ValueError("The pilot contact explainer assumes a y-invariant margin.")

    fig, axes = plt.subplots(1, 3, figsize=(10.4, 3.7), constrained_layout=True)
    colors = {electrode: color for electrode, _, _, color in ELECTRODES}

    for y, electrode in ((1, "G"), (0, "A")):
        left, right = mask_span_um(masks[electrode])
        axes[0].barh(y, right - left, left=left, height=0.52, color=colors[electrode])
        axes[0].text((left + right) / 2, y, electrode, ha="center", va="center", fontweight="bold")
    axes[0].set(
        xlim=(18, 42),
        yticks=(0, 1),
        yticklabels=("A · M2", "G · M1"),
        xlabel="x (µm)",
        title="1. 전극 투영",
    )

    for y, (_name, values, color) in enumerate(
        (("기존 방식 2D", legacy, "#AAB7BE"), ("기준 조건 2.5D", nominal, "#F08A38"))
    ):
        left, right = failure_span_um(values)
        axes[1].barh(y, right - left, left=left, height=0.5, color=color)
        axes[1].text(
            (left + right) / 2,
            y,
            f"{left:g}–{right:g} µm",
            ha="center",
            va="center",
            fontsize=8,
            fontweight="bold",
        )
    axes[1].set(
        xlim=(18, 42),
        yticks=(0, 1),
        yticklabels=("기존 방식 2D", "기준 조건 2.5D"),
        xlabel="이물 중심 x (µm)",
        title="2. 전극 조합 hotspot 구간",
    )

    x_um = np.arange(margin.shape[1], dtype=float) * 0.5
    profile = margin[0]
    selected = (x_um >= 26) & (x_um <= 34)
    axes[2].plot(x_um[selected], profile[selected], color="#173247", linewidth=2)
    axes[2].fill_between(
        x_um[selected],
        profile[selected],
        0,
        where=profile[selected] >= 0,
        color="#F08A38",
        alpha=0.75,
        interpolate=True,
        label="접촉",
    )
    axes[2].fill_between(
        x_um[selected],
        profile[selected],
        0,
        where=profile[selected] < 0,
        color="#D8E4E9",
        alpha=0.9,
        interpolate=True,
    )
    axes[2].axhline(0, color="#8AA0AA", linewidth=1, linestyle="--")
    axes[2].set(
        xlim=(26, 34),
        ylim=(-2.2, 0.55),
        xlabel="이물 중심 x (µm)",
        ylabel="접촉 여유 (µm)",
        title="3. 2.5D 접촉 판정",
    )
    axes[2].annotate(
        "접촉: 여유 ≥ 0",
        xy=(30.0, 0.2),
        xytext=(27.0, 0.42),
        arrowprops={"arrowstyle": "->", "color": "#B66125"},
        fontsize=8,
        color="#8B4A20",
    )
    for axis in axes:
        axis.grid(axis="x", color="#DCE6EA", linewidth=0.7)
        axis.spines[["top", "right"]].set_visible(False)
        axis.title.set_fontweight("bold")
        axis.title.set_color("#173247")
    fig.savefig(FIGURE_DIR / "explainer_contact.png", dpi=180, facecolor="white")
    plt.close(fig)


def open_fixture(name: str) -> tuple[np.ndarray, tuple[tuple[int, int], tuple[int, int]]]:
    wire = np.zeros((24, 24), dtype=bool)
    wire[11:14, 3:21] = True
    terminals = ((12, 3), (12, 20))
    if name == "detour":
        wire[5:8, 6:18] = True
        wire[5:14, 6:9] = True
        wire[5:14, 15:18] = True
    elif name == "branch":
        wire[3:12, 11:14] = True
    elif name == "corner":
        wire[:] = False
        wire[3:15, 6:9] = True
        wire[12:15, 6:21] = True
        terminals = ((3, 7), (13, 20))
    return wire, terminals


def make_open_explainer() -> None:
    manifest = load_manifest(CASE_DIR / "manifest.yaml")
    names = ("straight", "detour", "branch", "corner")
    rows = ((0.5, 0.5), (1.0, 0.5))
    fig, axes = plt.subplots(2, 4, figsize=(10.4, 5.1))
    for row, (pixel_y, pixel_x) in enumerate(rows):
        for column, name in enumerate(names):
            wire, terminals = open_fixture(name)
            grid = GridSpec(
                shape_yx=wire.shape,
                pixel_size_y_um=pixel_y,
                pixel_size_x_um=pixel_x,
            )
            masks = MaskSet(
                grid=grid,
                layers={"M2": wire},
                electrodes={"W": wire},
                input_hashes={},
            )
            path = OpenPath(
                path_id="W",
                electrode="W",
                terminal_a_yx_px=terminals[0],
                terminal_b_yx_px=terminals[1],
            )
            predicted = simulate_open(manifest, masks, path, 4.0).failure

            yy, xx = np.mgrid[:24, :24]
            brute = np.zeros_like(wire)
            for y, x in np.ndindex(wire.shape):
                blocked = blocked_mask(
                    np.hypot((yy - y) * pixel_y, (xx - x) * pixel_x), 2.0, 0.0
                )
                brute[y, x] = causes_open(wire, blocked, *terminals)
            if not np.array_equal(predicted, brute):
                raise ValueError(f"Open explainer mismatch: {name}, {pixel_y}, {pixel_x}")

            axis = axes[row, column]
            axis.imshow(
                np.ma.masked_where(~wire, wire),
                cmap=ListedColormap(["#8798A2"]),
                vmin=0,
                vmax=1,
                interpolation="nearest",
            )
            axis.imshow(
                np.ma.masked_where(~predicted, predicted),
                cmap=ListedColormap(["#FF6B2C"]),
                vmin=0,
                vmax=1,
                alpha=0.72,
                interpolation="nearest",
            )
            axis.contour(wire.astype(float), levels=[0.5], colors=["#455B68"], linewidths=0.7)
            axis.scatter(
                [terminal[1] for terminal in terminals],
                [terminal[0] for terminal in terminals],
                s=24,
                c="#18C7BF",
                edgecolors="#173247",
                linewidths=0.5,
                zorder=3,
            )
            if row == 0:
                axis.set_title({"straight": "직선", "detour": "우회", "branch": "분기", "corner": "ㄱ자"}[name], fontsize=10, fontweight="bold", color="#173247")
            if column == 0:
                axis.set_ylabel(f"dy={pixel_y:g}\ndx={pixel_x:g}\nµm/px", fontsize=8)
            axis.set_xticks([])
            axis.set_yticks([])

    handles = (
        Patch(facecolor="#8798A2", edgecolor="#455B68", label="배선"),
        Patch(facecolor="#FF6B2C", alpha=0.72, label="open을 만드는 이물 중심"),
        Line2D([], [], marker="o", linestyle="", color="#18C7BF", markeredgecolor="#173247", label="단자"),
    )
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=9)
    fig.subplots_adjust(left=0.06, right=0.99, top=0.91, bottom=0.14, wspace=0.12, hspace=0.12)
    fig.savefig(FIGURE_DIR / "explainer_open.png", dpi=180, facecolor="white")
    plt.close(fig)


def make_probability_profile() -> None:
    store = zarr.open_group(WORK_DIR / "prob" / "probability.zarr", mode="r")
    step_index = list(store.attrs["source_steps"]).index("S_M1")
    size_index = list(store.attrs["particle_sizes_um"]).index(10.0)
    pair_ids = list(store.attrs["unit_ids"])
    union_index = list(store.attrs["union_keys"]).index("any")
    pair_maps = np.asarray(store["probability"][step_index, size_index], dtype=float)
    union_map = np.asarray(
        store["union_probability"][step_index, size_index, union_index], dtype=float
    )
    if np.ptp(pair_maps, axis=1).max() != 0 or np.ptp(union_map, axis=0).max() != 0:
        raise ValueError("The pilot probability profile assumes y-invariant results.")

    x_um = np.arange(pair_maps.shape[-1], dtype=float) * 0.5
    shown = (x_um >= 22) & (x_um <= 58)
    colors = {"A__B": "#E8873A", "G__A": "#2796A6", "B__P": "#69A85B"}

    fig, axis = plt.subplots(figsize=(10.4, 3.8), constrained_layout=True)
    union_profile = union_map[0]
    axis.fill_between(
        x_um[shown],
        union_profile[shown],
        step="post",
        color="#DCE8ED",
        alpha=0.9,
    )
    axis.step(
        x_um[shown],
        union_profile[shown],
        where="post",
        color="#173247",
        linewidth=3,
        label="중복 제거 전체",
    )
    for index, pair_id in enumerate(pair_ids):
        axis.step(
            x_um[shown],
            pair_maps[index, 0, shown],
            where="post",
            linewidth=1.8,
            color=colors[pair_id],
            label=pair_id,
        )
    axis.set(
        xlim=(22, 58),
        ylim=(-0.03, 1.08),
        xlabel="이물 중심 x (µm)",
        ylabel="조건부 기하 확률",
        title="S_M1 · d=10 µm · 산포 64조건",
    )
    axis.set_yticks((0, 0.25, 0.5, 0.75, 1.0))
    axis.grid(color="#DCE6EA", linewidth=0.75)
    axis.spines[["top", "right"]].set_visible(False)
    axis.title.set_fontweight("bold")
    axis.title.set_color("#173247")
    axis.legend(loc="upper right", ncol=2, frameon=True, fontsize=8.5)
    axis.annotate(
        "0.5 = 64조건 중 32조건에서 접촉",
        xy=(36.0, 0.5),
        xytext=(23.0, 0.68),
        arrowprops={"arrowstyle": "->", "color": "#6F818B"},
        fontsize=8.5,
        color="#526B78",
    )
    fig.savefig(FIGURE_DIR / "explainer_probability_profile.png", dpi=180, facecolor="white")
    plt.close(fig)


def make_fab_comparison() -> None:
    data = pd.read_parquet(REPORT_DIR / "evidence" / "prob" / "fab_risk.parquet")
    grouped = (
        data.groupby(["source_step", "weighting_mode"], sort=False)["expected_defect_count"]
        .sum()
        .unstack(fill_value=0)
        .reindex(index=["S_M1", "S_M2", "S_PLN"], columns=["actual", "common"])
    )
    x = np.arange(len(grouped.index), dtype=float)
    fig, axis = plt.subplots(figsize=(7.0, 5.0), constrained_layout=True)
    axis.bar(x - 0.2, grouped["actual"], width=0.4, label="실제 Fab 분포", color="#287DB2")
    axis.bar(x + 0.2, grouped["common"], width=0.4, label="공통 이물 분포", color="#F28A2B")
    axis.set(
        xticks=x,
        xticklabels=grouped.index,
        ylabel="모델상 기대 결함 수 (CA × 밀도)",
        title="공통 이물 분포와 실제 Fab 분포 적용 결과",
    )
    axis.grid(axis="y", color="#DCE6EA", linewidth=0.8)
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend()
    fig.savefig(FIGURE_DIR / "probability_fab_comparison.png", dpi=180, facecolor="white")
    plt.close(fig)


def main() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    make_layout()
    make_layout_3d_draft()
    make_contact_explainer()
    make_open_explainer()
    make_probability_profile()
    make_fab_comparison()


if __name__ == "__main__":
    main()
