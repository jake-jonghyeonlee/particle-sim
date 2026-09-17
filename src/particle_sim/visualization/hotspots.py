"""Review figures for nominal runs (PNG, matplotlib Agg)."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from particle_sim.io.masks import GridSpec, MaskSet  # noqa: E402

LAYER_COLORS: tuple[str, ...] = ("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")


def _extent(grid: GridSpec) -> tuple[float, float, float, float]:
    height = grid.shape_yx[0] * grid.pixel_size_y_um
    width = grid.shape_yx[1] * grid.pixel_size_x_um
    return (0.0, width, height, 0.0)


def _save(fig: plt.Figure, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=100)
    plt.close(fig)
    return path


def plot_layout_overlay(masks: MaskSet, layer_order: list[str], path: Path) -> Path:
    grid = masks.grid
    fig, ax = plt.subplots(figsize=(6, 6))
    handles = []
    for index, name in enumerate(layer_order):
        color = LAYER_COLORS[index % len(LAYER_COLORS)]
        rgba = np.zeros((*grid.shape_yx, 4), dtype=float)
        rgba[..., :3] = matplotlib.colors.to_rgb(color)
        rgba[..., 3] = masks.layers[name].astype(float) * 0.5
        ax.imshow(rgba, extent=_extent(grid), interpolation="nearest")
        handles.append(Patch(color=color, alpha=0.5, label=name))
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title("Layout overlay")
    if handles:
        ax.legend(handles=handles, loc="upper right")
    return _save(fig, path)


def plot_pair_panels(
    legacy: np.ndarray,
    nominal: np.ndarray,
    margin_um: np.ndarray,
    grid: GridSpec,
    title: str,
    path: Path,
) -> Path:
    extent = _extent(grid)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, data, label in (
        (axes[0], legacy, "Legacy erode/dilate (binary)"),
        (axes[1], nominal, "Nominal 2.5D contact (binary)"),
    ):
        image = ax.imshow(
            np.asarray(data, dtype=float), cmap="gray_r", vmin=0, vmax=1,
            extent=extent, interpolation="nearest",
        )
        ax.set_title(label)
        ax.set_xlabel("x (um)")
        ax.set_ylabel("y (um)")
        fig.colorbar(image, ax=ax, fraction=0.046)
    margin = np.asarray(margin_um, dtype=float)
    finite = np.isfinite(margin)
    bound = float(np.abs(margin[finite]).max()) if finite.any() else 0.0
    bound = max(bound, 1e-9)
    shown = np.where(finite, margin, np.nan)
    image = axes[2].imshow(
        shown, cmap="RdBu", vmin=-bound, vmax=bound, extent=extent, interpolation="nearest"
    )
    axes[2].set_title("Pair contact margin (um)")
    axes[2].set_xlabel("x (um)")
    axes[2].set_ylabel("y (um)")
    fig.colorbar(image, ax=axes[2], fraction=0.046)
    fig.suptitle(title)
    fig.tight_layout()
    return _save(fig, path)


def plot_ca_curves(ca: pd.DataFrame, step_id: str, path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(7, 5))
    rows = ca[ca["source_step"] == step_id]
    for index, (pair_id, group) in enumerate(rows.groupby("pair_id", sort=False)):
        group = group.sort_values("particle_size_um")
        color = LAYER_COLORS[index % len(LAYER_COLORS)]
        ax.plot(group["particle_size_um"], group["ca_nominal_um2"], "-o", color=color,
                label=f"{pair_id} nominal 2.5D")
        ax.plot(group["particle_size_um"], group["ca_legacy_um2"], "--", color=color,
                label=f"{pair_id} legacy")
    ax.set_xlabel("particle diameter (um)")
    ax.set_ylabel("critical area (um^2)")
    ax.set_title(f"CA by particle size - source step {step_id}")
    if len(rows):
        ax.legend(fontsize="small")
    return _save(fig, path)


def plot_union_vs_sum(union: pd.DataFrame, step_id: str, path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(7, 5))
    rows = union[union["source_step"] == step_id]
    for index, (mode, group) in enumerate(rows.groupby("failure_mode", sort=False)):
        group = group.sort_values("particle_size_um")
        color = LAYER_COLORS[index % len(LAYER_COLORS)]
        ax.plot(group["particle_size_um"], group["ca_union_nominal_um2"], "-o", color=color,
                label=f"{mode} union")
        ax.plot(group["particle_size_um"], group["ca_pair_sum_nominal_um2"], ":", color=color,
                label=f"{mode} pair sum")
    ax.set_xlabel("particle diameter (um)")
    ax.set_ylabel("critical area (um^2)")
    ax.set_title(f"Union CA vs pair sum - source step {step_id}")
    if len(rows):
        ax.legend(fontsize="small")
    return _save(fig, path)
