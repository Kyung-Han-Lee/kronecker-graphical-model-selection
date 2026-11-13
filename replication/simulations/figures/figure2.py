from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from . import style

FILES = ("figure2_path.csv",)


LOGDET_MAX = -np.log(2.0)


LOGDET_FLOOR = -4.0

Y_MIN_B = -1.75


ARROWS = ((18, 24), (58, 66))
C_LMAX = "0.45"


def landscape(n: int = 500) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    theta = np.radians(np.linspace(0.5, 89.5, n))
    z = np.linspace(-0.75, 0.75, n)
    theta_grid, z_grid = np.meshgrid(theta, z)
    det = ((1.0 - z_grid**2) * np.sin(2.0 * theta_grid) - z_grid**2) / 2.0
    values = np.where(det > 0, np.log(np.maximum(det, 1e-12)), np.nan)
    return np.degrees(theta_grid), z_grid, values


def draw_panel_a(axis, figure, path: pd.DataFrame, lambda_max: float) -> None:
    theta_grid, z_grid, values = landscape()
    levels = np.linspace(LOGDET_FLOOR, LOGDET_MAX, 14)
    filled = axis.contourf(theta_grid, z_grid, values, levels=levels, cmap="viridis")

    axis.contour(
        theta_grid,
        z_grid,
        values,
        levels=levels[1:-1:2],
        colors="white",
        linewidths=0.5,
        alpha=0.6,
    )
    axis.axhline(0, color=style.color("empty_graph"), lw=style.line("main"), zorder=3)
    theta = path.theta.to_numpy()
    z = path.z.to_numpy()
    color = style.color("path")
    axis.plot(
        theta,
        z,
        color=color,
        lw=style.line("emphasis"),
        solid_capstyle="round",
        zorder=5,
    )
    arrow = {
        "arrowstyle": "-|>,head_width=0.4,head_length=0.8",
        "lw": style.line("main"),
        "shrinkA": 0,
        "shrinkB": 0,
        "mutation_scale": 7,
    }
    for a, b in ARROWS:
        axis.annotate(
            "",
            xy=(theta[b], z[b]),
            xytext=(theta[a], z[a]),
            arrowprops={**arrow, "color": color},
            zorder=7,
        )

    axis.scatter(
        [theta[0]],
        [z[0]],
        facecolor="white",
        edgecolor="k",
        s=30,
        lw=0.8,
        zorder=8,
    )
    k = int(np.argmax(path["lambda"].to_numpy() > lambda_max))
    axis.scatter([theta[k]], [z[k]], color=C_LMAX, marker="D", s=26, zorder=8)
    axis.scatter([45.0], [0.0], color="k", marker="*", s=70, zorder=9)
    axis.set_xlim(0, 90)
    axis.set_ylim(-0.75, 0.75)
    axis.set_xticks([0, 30, 45, 60, 90])
    axis.set_xticklabels(["0°", "30°", "45°", "60°", "90°"])

    axis.set_xlabel(r"$\theta$", labelpad=5)
    axis.set_ylabel(r"$z$", labelpad=6)
    style.panel_label(axis, "(a)")

    colorbar = figure.colorbar(filled, ax=axis, shrink=1.0, aspect=28, pad=0.02)
    colorbar.ax.set_title(
        r"$\log\det\boldsymbol{\Omega}$",
        fontsize=style.STYLE["font"]["title_pt"],
        pad=6,
    )
    colorbar.set_ticks([-4, -3, -2, -1, LOGDET_MAX])
    colorbar.set_ticklabels(["$-4$", "$-3$", "$-2$", "$-1$", r"$-\log 2$"])


def draw_panel_b(axis, path: pd.DataFrame, lambda_max: float) -> None:
    axis.plot(
        path["lambda"],
        path.log_det,
        color=style.color("path"),
        lw=style.line("emphasis"),
    )
    axis.axhline(LOGDET_MAX, color="k", lw=style.line("reference"), ls=":")
    axis.axvline(lambda_max, color="k", lw=style.line("reference"), ls="-.")

    axis.text(
        5.95,
        LOGDET_MAX + 0.012,
        "upper bound",
        ha="right",
        va="bottom",
        fontsize=style.STYLE["font"]["legend_pt"],
    )
    axis.set_xlim(0, 6.0)
    axis.set_ylim(Y_MIN_B, -0.6)
    axis.set_xticks([0, lambda_max, 2, 4, 6])
    axis.set_xticklabels(["$0$", r"$\lambda_{\max}$", "$2$", "$4$", "$6$"])
    ticks = [-1.6, -1.4, -1.2, -1.0, -0.8, LOGDET_MAX]
    axis.set_yticks(ticks)
    axis.set_yticklabels(
        ["$-1.6$", "$-1.4$", "$-1.2$", "$-1.0$", "$-0.8$", r"$-\log 2$"]
    )
    axis.set_xlabel(r"$\lambda$", labelpad=5)

    axis.set_ylabel(r"$\log\det\boldsymbol{\Omega}$", labelpad=-2)
    style.panel_label(axis, "(b)")
    axis.spines[["top", "right"]].set_visible(False)


def centre_horizontally(figure: Figure) -> None:

    dpi = figure.dpi
    figure.set_dpi(300)
    figure.canvas.draw()
    image = np.asarray(figure.canvas.buffer_rgba())
    columns = np.flatnonzero(np.any(image[..., :3] < 250, axis=(0, 2)))
    left = columns[0]
    right = image.shape[1] - 1 - columns[-1]
    shift = (right - left) / 2 / image.shape[1]
    figure.set_dpi(dpi)
    figure.set_layout_engine("none")
    for axis in figure.axes:
        position = axis.get_position()
        axis.set_position(
            [position.x0 + shift, position.y0, position.width, position.height]
        )


def draw(path: pd.DataFrame) -> Figure:
    style.apply()
    figure, (axis_a, axis_b) = plt.subplots(
        1,
        2,
        figsize=(style.width(), style.height("figure2")),
        layout="constrained",
        gridspec_kw={"width_ratios": [1.15, 1]},
    )
    lambda_max = float(path.lambda_max.iloc[0])
    draw_panel_a(axis_a, figure, path[path.panel == "a"], lambda_max)
    draw_panel_b(axis_b, path[path.panel == "b"], lambda_max)

    for axis in (axis_a, axis_b):
        axis.set_box_aspect(1)
    centre_horizontally(figure)
    return figure


def make(folder) -> Figure:
    path = pd.read_csv(folder / FILES[0]).sort_values(["panel", "point"])
    return draw(path)
