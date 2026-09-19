"""Figures for the four failure categories: lateral short, vertical short, open and CNO.

Schematics explain each algorithm; result panels read the saved pilot runs under
docs/report/work (run benchmark_contact_pilot.py first). Algorithm panels call the real
judgement functions on a small synthetic wire so the pictures show computed output.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

# isort: off
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import zarr  # noqa: E402
from matplotlib.patches import Circle, Rectangle  # noqa: E402
from PIL import Image  # noqa: E402
from scipy.ndimage import label  # noqa: E402

from particle_sim.failure.auto_open import simulate_auto_open  # noqa: E402
from particle_sim.failure.contact_not_open import simulate_cno  # noqa: E402
from particle_sim.failure.open_failure import particle_footprint  # noqa: E402
from particle_sim.geometry.contacts import ContactHole, ContactLink, ContactTopology  # noqa: E402
from particle_sim.io.masks import GridSpec, MaskSet  # noqa: E402
# isort: on

REPORT = Path(__file__).resolve().parent
FIG = REPORT / "figures"
PILOT = REPORT / "work" / "contact-pilot-200"
EVIDENCE = REPORT / "evidence"

plt.rcParams.update(
    {
        "font.family": ["Apple SD Gothic Neo", "Arial Unicode MS", "DejaVu Sans"],
        "axes.unicode_minus": False,
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)
NAVY, ORANGE, BLUE, GREEN, PURPLE, GRAY = (
    "#193247", "#e8833a", "#3b82c4", "#3a9d6b", "#8e5bb5", "#b8c2cc",
)
MODE_COLORS = {
    "lateral short": BLUE,
    "vertical short": GREEN,
    "open": ORANGE,
    "CNO": PURPLE,
}
MODE_NAMES = {
    "lateral_short": "lateral short",
    "vertical_short": "vertical short",
    "lateral_open": "open",
    "contact_not_open": "CNO",
}


def _save(fig, name: str) -> None:
    fig.savefig(FIG / name, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _clean(ax) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def _mask(name: str) -> np.ndarray:
    return np.asarray(Image.open(PILOT / "contact-case" / "masks" / name).convert("L")) > 127


def _outline(ax, mask: np.ndarray, color: str, width: float = 1.2, label_text=None) -> None:
    ax.contour(mask.astype(float), levels=[0.5], colors=[color], linewidths=width)
    if label_text:
        ax.plot([], [], color=color, lw=width, label=label_text)


def _pilot():
    result = PILOT / "contact_terminals-result"
    store = zarr.open_group(result / "hotspots.zarr", mode="r")
    ca = pd.read_parquet(result / "ca.parquet")
    union = pd.read_parquet(result / "ca_union.parquet")
    explicit = zarr.open_group(PILOT / "explicit_terminals-result" / "hotspots.zarr", mode="r")
    return store, ca, union, explicit


def _map(store, unit: str, step: str, size: float) -> np.ndarray:
    steps = list(store.attrs["source_steps"])
    sizes = [float(s) for s in store.attrs["particle_sizes_um"]]
    units = list(store.attrs["failure_unit_ids"])
    return np.asarray(
        store["nominal_failure"][steps.index(step), sizes.index(size), units.index(unit)],
        dtype=bool,
    )


def _hotspot_panel(ax, hotspot, outlines, title, crop=None, color=ORANGE, legend="lower left"):
    view = np.zeros((*hotspot.shape, 4))
    view[hotspot] = matplotlib.colors.to_rgba(color, 0.85)
    ax.imshow(np.ones(hotspot.shape), cmap="gray", vmin=0, vmax=1)
    ax.imshow(view)
    for mask, color, text in outlines:
        _outline(ax, mask, color, label_text=text)
    if crop:
        (y0, y1), (x0, x1) = crop
        ax.set_xlim(x0, x1)
        ax.set_ylim(y1, y0)
    _clean(ax)
    ax.set_title(title)
    ax.legend(loc=legend, fontsize=7, frameon=True, framealpha=0.9)


def _ca_by_size(ax, ca, mode: str, title: str, color: str) -> None:
    rows = ca[ca["failure_mode"] == mode]
    for step, style in zip(("S_M1", "S_M2", "S_PLN"), ("-", "--", ":"), strict=True):
        series = rows[rows["source_step"] == step].groupby("particle_size_um")["ca_nominal_um2"].sum()
        ax.plot(series.index, series.values, style, color=color, marker="o", ms=3,
                label=f"{step} 이후 유입")
    ax.set_xlabel("이물 직경 d [µm]")
    ax.set_ylabel("조합 CA 합 [µm²]")
    ax.set_title(title)
    ax.legend(fontsize=7, frameon=False)
    ax.grid(alpha=0.25)


# ------------------------------------------------------------------ overview


def overview() -> None:
    fig, axes = plt.subplots(1, 4, figsize=(11, 2.9))
    # lateral short: top view
    ax = axes[0]
    ax.add_patch(Rectangle((0, 0), 1.2, 4, color=BLUE, alpha=0.8))
    ax.add_patch(Rectangle((2.2, 0), 1.2, 4, color="#7fb3de", alpha=0.9))
    ax.add_patch(Circle((1.7, 2), 0.62, color=ORANGE, alpha=0.9))
    ax.text(0.6, 4.3, "전극 A", ha="center")
    ax.text(2.8, 4.3, "전극 B", ha="center")
    ax.set_title("lateral short\n같은 layer 두 전극에 동시 접촉")
    ax.set_xlim(-0.4, 3.8)
    ax.set_ylim(-0.5, 4.9)
    # vertical short: section
    ax = axes[1]
    ax.add_patch(Rectangle((0, 0), 1.6, 0.35, color=GREEN))
    ax.add_patch(Rectangle((0, 0.35), 4, 1.0, color=GRAY, alpha=0.5))
    ax.add_patch(Rectangle((2.2, 1.35), 1.6, 0.35, color="#8fd0ad"))
    ax.add_patch(Circle((1.85, 1.0), 0.72, color=ORANGE, alpha=0.9))
    ax.text(0.8, -0.35, "아래 전극 (mask1)", ha="center", fontsize=8)
    ax.text(3.0, 2.0, "위 전극 (mask2)", ha="center", fontsize=8)
    ax.set_title("vertical short\n다른 높이 두 전극에 동시 접촉")
    ax.set_xlim(-0.2, 4.2)
    ax.set_ylim(-0.7, 2.5)
    # open: top view
    ax = axes[2]
    ax.add_patch(Rectangle((0, 1.6), 4, 0.8, color=NAVY, alpha=0.75))
    ax.add_patch(Rectangle((0, 1.6), 0.5, 0.8, color=BLUE))
    ax.add_patch(Rectangle((3.5, 1.6), 0.5, 0.8, color=BLUE))
    ax.add_patch(Circle((2, 2), 0.62, color=ORANGE, alpha=0.9))
    ax.text(0.25, 2.8, "단자", ha="center", fontsize=8)
    ax.text(3.75, 2.8, "단자", ha="center", fontsize=8)
    ax.set_title("open\n두 단자 사이 도체 경로가 끊김")
    ax.set_xlim(-0.3, 4.3)
    ax.set_ylim(0.5, 3.4)
    # CNO: section
    ax = axes[3]
    ax.add_patch(Rectangle((0, 0), 4, 0.35, color=NAVY, alpha=0.8))
    ax.add_patch(Rectangle((0, 0.35), 1.5, 0.9, color=GRAY, alpha=0.6))
    ax.add_patch(Rectangle((2.5, 0.35), 1.5, 0.9, color=GRAY, alpha=0.6))
    ax.add_patch(Rectangle((0, 1.25), 4, 0.3, color="#7f95a8"))
    ax.add_patch(Circle((2.0, 0.8), 0.46, color=PURPLE, alpha=0.9))
    ax.text(2.0, -0.35, "아래 전극", ha="center", fontsize=8)
    ax.text(2.0, 1.75, "위 전극", ha="center", fontsize=8)
    ax.text(3.3, 0.72, "절연막", ha="center", fontsize=8)
    ax.set_title("CNO (Contact Not Open)\ncontact hole이 막혀 층간 연결 실패")
    ax.set_xlim(-0.2, 4.2)
    ax.set_ylim(-0.7, 2.2)
    for ax in axes:
        ax.set_aspect("equal")
        _clean(ax)
    fig.tight_layout()
    _save(fig, "fm_overview.png")


# ------------------------------------------------------------------ short modes


def lateral(store, ca) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), gridspec_kw={"width_ratios": [1, 1.1, 1.2]})
    ax = axes[0]
    ax.add_patch(Rectangle((0, 0), 1.0, 4, color=BLUE, alpha=0.8))
    ax.add_patch(Rectangle((2.4, 0), 1.0, 4, color="#7fb3de"))
    ax.add_patch(Circle((1.7, 2.6), 0.72, color=ORANGE, alpha=0.9))
    ax.add_patch(Circle((1.7, 1.0), 0.45, fill=False, ec=NAVY, ls="--", lw=1.2))
    ax.annotate("", xy=(1.0, 2.6), xytext=(2.4, 2.6),
                arrowprops={"arrowstyle": "<->", "color": NAVY, "lw": 0.9})
    ax.text(1.7, 4.2, "큰 이물: 두 전극 모두 접촉 → hotspot", ha="center", fontsize=8)
    ax.text(1.7, -0.45, "작은 이물(점선): 접촉 부족 → 정상", ha="center", fontsize=8)
    ax.set_xlim(-0.2, 3.6)
    ax.set_ylim(-0.8, 4.5)
    ax.set_aspect("equal")
    _clean(ax)
    ax.set_title("① 판정: 원판이 A·B에 모두 닿는가")
    hotspot = _map(store, "D1__D2", "S_M1", 10.0)
    _hotspot_panel(
        axes[1], hotspot,
        [(_mask("electrode_D1.png"), BLUE, "D1 (mask2)"), (_mask("electrode_D2.png"), NAVY, "D2 (mask2)")],
        "② D1–D2, d=10 µm hotspot",
    )
    _ca_by_size(axes[2], ca, "lateral_short", "③ lateral short CA 합 (6개 조합)", BLUE)
    fig.tight_layout()
    _save(fig, "fm_lateral.png")


def vertical(store, ca) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), gridspec_kw={"width_ratios": [1.1, 1.1, 1.2]})
    ax = axes[0]
    ax.add_patch(Rectangle((0, 0), 1.6, 0.3, color=GREEN))
    ax.add_patch(Rectangle((2.6, 1.6), 1.6, 0.3, color="#8fd0ad"))
    ax.add_patch(Circle((1.9, 0.95), 0.9, color=ORANGE, alpha=0.35))
    ax.plot([1.9], [0.95], "o", color=NAVY, ms=3)
    ax.plot([1.9, 2.6], [0.95, 1.6], color=NAVY, lw=1.2)
    ax.plot([1.9, 2.6], [0.95, 0.95], ":", color=NAVY, lw=1)
    ax.plot([2.6, 2.6], [0.95, 1.6], ":", color=NAVY, lw=1)
    ax.text(2.25, 0.78, "Dxy", fontsize=8, ha="center")
    ax.text(2.72, 1.25, "Dz", fontsize=8)
    ax.text(2.0, 1.4, "Q", fontsize=9, color=NAVY, fontweight="bold")
    ax.text(0.8, -0.3, "아래 전극", ha="center", fontsize=8)
    ax.text(3.4, 2.1, "위 전극", ha="center", fontsize=8)
    ax.text(2.1, -0.75, "Q = √(Dxy² + Dz²) ≤ r 이면 접촉", ha="center", fontsize=8)
    ax.set_xlim(-0.2, 4.4)
    ax.set_ylim(-1.0, 2.5)
    ax.set_aspect("equal")
    _clean(ax)
    ax.set_title("① 판정: 두 전극까지 대각 거리 Q")
    hotspot = _map(store, "G2__D1", "S_M1", 10.0)
    _hotspot_panel(
        axes[1], hotspot,
        [(_mask("electrode_G2.png"), GREEN, "G2 (mask1)"), (_mask("electrode_D1.png"), NAVY, "D1 (mask2)")],
        "② G2–D1, d=10 µm hotspot",
    )
    _ca_by_size(axes[2], ca, "vertical_short", "③ vertical short CA 합 (10개 조합)", GREEN)
    fig.tight_layout()
    _save(fig, "fm_vertical.png")


# ------------------------------------------------------------------ open algorithm


def _toy():
    """A 60x90 wire with two contacts and a stub, evaluated with the real functions."""
    shape = (60, 90)
    grid = GridSpec(shape_yx=shape, pixel_size_y_um=0.5, pixel_size_x_um=0.5)
    wire = np.zeros(shape, dtype=bool)
    wire[28:34, 8:82] = True
    wire[8:28, 44:49] = True  # stub without a terminal
    holes = []
    for index, (x0, x1) in enumerate(((10, 16), (74, 80)), start=1):
        region = np.zeros(shape, dtype=bool)
        region[29:33, x0:x1] = True
        ys, xs = np.nonzero(region)
        holes.append(ContactHole("CNT", index, f"PAD{index}", "W", ys, xs))
    links = {
        f"cno:CNT:PAD{h.hole_index}__W": ContactLink(
            f"cno:CNT:PAD{h.hole_index}__W", "CNT", f"PAD{h.hole_index}", "W", (h,)
        )
        for h in holes
    }
    masks = MaskSet(
        grid=grid, layers={}, electrodes={"W": wire}, input_hashes={},
        contact_topology=ContactTopology(tuple(holes), links),
    )
    return grid, wire, holes, links, masks


def open_algorithm() -> None:
    grid, wire, holes, links, masks = _toy()
    diameter = 4.0
    stamp = particle_footprint(grid, diameter / 2, 0.0)
    half = stamp.shape[0] // 2
    result = simulate_auto_open(masks, "W", diameter, 0.0)
    cno = np.zeros(wire.shape, dtype=bool)
    for link in links.values():
        cno |= simulate_cno(link, link.unit_id, grid, 8.0, 0.0).failure

    def damaged_at(y, x):
        blocked = np.zeros(wire.shape, dtype=bool)
        blocked[y - half : y + half + 1, x - half : x + half + 1] = stamp
        return wire & ~blocked, blocked

    fig, grid_axes = plt.subplots(2, 2, figsize=(9.2, 4.6))
    axes = grid_axes.ravel()
    palette = ["#5b8fb9", "#6aa16f", "#c3a1dd", "#9aa5b1"]

    def draw_components(ax, damaged, blocked, title, verdict, color):
        labels, count = label(damaged)
        view = np.ones((*wire.shape, 3))
        for index in range(1, count + 1):
            view[labels == index] = matplotlib.colors.to_rgb(palette[(index - 1) % len(palette)])
        view[blocked & wire] = matplotlib.colors.to_rgb("#f3c7a6")
        ax.imshow(view)
        for hole in holes:
            ax.add_patch(Rectangle((hole.xs.min() - 0.5, hole.ys.min() - 0.5),
                                   np.ptp(hole.xs) + 1, np.ptp(hole.ys) + 1,
                                   fill=False, ec="white", lw=1.4))
        ax.set_title(title)
        ax.text(45, 57, verdict, ha="center", va="bottom", fontsize=9, color=color,
                fontweight="bold")
        _clean(ax)

    view = np.ones((*wire.shape, 3))
    view[wire] = matplotlib.colors.to_rgb(GRAY)
    ax = axes[0]
    ax.imshow(view)
    for hole in holes:
        ax.add_patch(Rectangle((hole.xs.min() - 0.5, hole.ys.min() - 0.5),
                               np.ptp(hole.xs) + 1, np.ptp(hole.ys) + 1, color=BLUE))
    ax.text(13, 26, "contact 단자", ha="center", fontsize=7, color=BLUE)
    ax.text(77, 26, "contact 단자", ha="center", fontsize=7, color=BLUE)
    ax.text(58, 14, "단자 없는 가지", fontsize=7, color=NAVY)
    ax.set_title("① 입력: 전극 + 자동 단자")
    _clean(ax)
    damaged, blocked = damaged_at(31, 30)
    draw_components(axes[1], damaged, blocked, "② 본선 절단", "단자 분리 → open", ORANGE)
    damaged, blocked = damaged_at(18, 46)
    draw_components(axes[2], damaged, blocked, "③ 가지 절단", "단자 없는 조각 → 무시", GREEN)
    ax = axes[3]
    overlay = np.ones((*wire.shape, 3))
    overlay[wire] = matplotlib.colors.to_rgb("#dfe5ea")
    overlay[result.failure] = matplotlib.colors.to_rgb(ORANGE)
    overlay[cno] = matplotlib.colors.to_rgb(PURPLE)
    ax.imshow(overlay)
    ax.plot([], [], "s", color=ORANGE, label="open (d=4 µm)")
    ax.plot([], [], "s", color=PURPLE, label="CNO (d=8 µm)")
    ax.legend(loc="lower center", fontsize=7, ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.02))
    ax.set_title("④ 계산 결과: 이물 중심 위치")
    _clean(ax)
    fig.tight_layout()
    _save(fig, "fm_open_algorithm.png")


def open_pilot(store, explicit) -> dict:
    step, size = "S_M2", 10.0
    auto = _map(store, "auto_open:G1", step, size)
    steps = list(explicit.attrs["source_steps"])
    sizes = [float(s) for s in explicit.attrs["particle_sizes_um"]]
    units = list(explicit.attrs["failure_unit_ids"])
    manual = np.asarray(
        explicit["nominal_failure"][steps.index(step), sizes.index(size), units.index("G1_OPEN")],
        dtype=bool,
    )
    moved = manual & ~auto
    electrode = _mask("electrode_G1.png")
    pads = _mask("electrode_G1_PAD_A.png") | _mask("electrode_G1_PAD_B.png")
    ys, xs = np.nonzero(electrode)
    crop = ((max(0, ys.min() - 12), min(200, ys.max() + 12)), (max(0, xs.min() - 12), min(200, xs.max() + 12)))
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.1))
    _hotspot_panel(axes[0], manual, [(electrode, NAVY, "G1 (mask1)"), (pads, BLUE, "단자 영역")],
                   f"① 명시 단자 mask: {int(manual.sum()):,} px", crop, legend="center right")
    view = np.zeros((*auto.shape, 4))
    view[auto] = matplotlib.colors.to_rgba(ORANGE, 0.85)
    view[moved] = matplotlib.colors.to_rgba(PURPLE, 0.95)
    ax = axes[1]
    ax.imshow(np.ones(auto.shape), cmap="gray", vmin=0, vmax=1)
    ax.imshow(view)
    _outline(ax, electrode, NAVY, label_text="G1 (mask1)")
    _outline(ax, pads, BLUE, label_text="contact hole")
    ax.plot([], [], "s", color=ORANGE, label="open")
    ax.plot([], [], "s", color=PURPLE, label="contact 막힘 → CNO")
    (y0, y1), (x0, x1) = crop
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    _clean(ax)
    ax.set_title(f"② contact 자동 단자: {int(auto.sum()):,} px + CNO 이동 {int(moved.sum()):,} px")
    ax.legend(loc="center right", fontsize=7, framealpha=0.9)
    fig.tight_layout()
    _save(fig, "fm_open_pilot.png")
    return {"explicit": int(manual.sum()), "auto": int(auto.sum()), "moved": int(moved.sum())}


# ------------------------------------------------------------------ CNO


def cno(store, ca) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.5), gridspec_kw={"width_ratios": [1.15, 1, 1.1]})
    ax = axes[0]
    for ox, radius, verdict, color in ((0.0, 1.35, "전부 가림 → CNO", PURPLE),
                                       (3.4, 0.8, "일부 가림 → 정상(잠정)", GRAY)):
        ax.add_patch(Rectangle((ox + 0.2, 0.6), 1.6, 1.2, fill=False, ec=NAVY, lw=1.4))
        ax.add_patch(Circle((ox + 1.0 + (0 if radius > 1 else 0.5), 1.2), radius,
                            color=color, alpha=0.55))
        ax.text(ox + 1.0, -0.55, verdict, ha="center", fontsize=8,
                color=PURPLE if radius > 1 else NAVY)
    ax.text(1.0, 2.75, "contact hole", ha="center", fontsize=8)
    ax.set_xlim(-0.6, 6.0)
    ax.set_ylim(-0.9, 3.0)
    ax.set_aspect("equal")
    _clean(ax)
    ax.set_title("① 판정: 원판이 hole 접촉 픽셀을 모두 덮는가")
    hotspot = _map(store, "cno:CNT_GI:G1__G1_PAD_A", "S_M2", 10.0)
    hole = _mask("electrode_G1_PAD_A.png")
    ys, xs = np.nonzero(hole)
    crop = ((ys.min() - 10, ys.max() + 11), (xs.min() - 10, xs.max() + 11))
    _hotspot_panel(axes[1], hotspot, [(hole, BLUE, "contact hole"), (_mask("electrode_G1.png"), NAVY, "G1")],
                   "② G1 contact, d=10 µm CNO 위치", crop, color=PURPLE)
    rows = ca[ca["failure_mode"] == "contact_not_open"]
    series = rows.groupby("particle_size_um")["ca_nominal_um2"].sum()
    axes[2].bar(series.index, series.values, width=1.2, color=PURPLE)
    for x, y in series.items():
        axes[2].text(x, y + 1, f"{y:g}", ha="center", fontsize=8)
    axes[2].set_xlabel("이물 직경 d [µm]")
    axes[2].set_ylabel("CNO CA 합 [µm²]")
    axes[2].set_title("③ CNO CA (hole 6개, 적용 공정 합)")
    axes[2].grid(axis="y", alpha=0.25)
    fig.tight_layout()
    _save(fig, "fm_cno.png")


# ------------------------------------------------------------------ mix and TAT


def mix(ca, union) -> None:
    probability = json.loads((EVIDENCE / "benchmark-contact-probability-200.json").read_text(encoding="utf-8"))
    nominal = json.loads((EVIDENCE / "benchmark-contact-pilot-200.json").read_text(encoding="utf-8"))
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.3), gridspec_kw={"width_ratios": [1.6, 1, 1]})
    ax = axes[0]
    pivot = ca.assign(mode=ca["failure_mode"].map(MODE_NAMES)).pivot_table(
        index="particle_size_um", columns="mode", values="ca_allocated_nominal_um2", aggfunc="sum"
    )
    bottom = np.zeros(len(pivot))
    for mode in ("lateral short", "vertical short", "open", "CNO"):
        values = pivot[mode].to_numpy()
        ax.bar(pivot.index, values, bottom=bottom, width=1.3, color=MODE_COLORS[mode], label=mode)
        bottom += values
    totals = union[union["failure_mode"] == "any"].groupby("particle_size_um")["ca_union_nominal_um2"].sum()
    ax.plot(totals.index, totals.values, "k_", ms=22, mew=2, label="union any")
    ax.set_xlabel("이물 직경 d [µm]")
    ax.set_ylabel("1/k 배분 CA [µm²] (3개 공정 합)")
    ax.set_title("① 4개 유형 배분 CA = union any")
    ax.legend(fontsize=7, frameon=False, ncol=2)
    ax.grid(axis="y", alpha=0.25)
    labels = ["명시 단자", "contact 단자"]
    keys = ["explicit_terminals", "contact_terminals"]
    values = [nominal["tat_summary"][k]["wall_median_seconds"] for k in keys]
    axes[1].bar(labels, values, color=[GRAY, PURPLE])
    for x, v in enumerate(values):
        axes[1].text(x, v * 1.02, f"{v:.2f} s", ha="center", fontsize=8)
    axes[1].set_title("② 기준 조건 TAT (5회 중앙값)")
    axes[1].set_ylabel("초")
    values = [probability["results"][k]["wall_seconds"] for k in keys]
    axes[2].bar(labels, values, color=[GRAY, PURPLE])
    for x, v in enumerate(values):
        axes[2].text(x, v * 1.02, f"{v:.0f} s", ha="center", fontsize=8)
    axes[2].set_title("③ 산포 64조건 TAT (1회)")
    axes[2].set_ylabel("초")
    for ax in axes[1:]:
        ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    _save(fig, "fm_mix_tat.png")


def main() -> None:
    FIG.mkdir(exist_ok=True)
    store, ca, union, explicit = _pilot()
    overview()
    lateral(store, ca)
    vertical(store, ca)
    open_algorithm()
    counts = open_pilot(store, explicit)
    cno(store, ca)
    mix(ca, union)
    print(json.dumps(counts))


if __name__ == "__main__":
    main()
