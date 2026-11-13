from __future__ import annotations

import math

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.ticker import LogFormatterMathtext, StrMethodFormatter

from . import style

CRITERIA = ("ebic", "scale", "relaxation")
LABELS = {
    "ebic": r"$\operatorname{EBIC}_{\gamma}$",
    "scale": r"$\operatorname{EBIC}^{\mathrm{scale}}_{\gamma}$",
    "relaxation": r"$\operatorname{EBIC}^{\mathrm{relax}}_{\gamma}$",
}
FILES = ("figure1_curves.csv", "figure1_medians.csv")


def x_ticks(eta: np.ndarray) -> list[float]:
    low = math.ceil(math.log10(eta.min()))
    high = math.floor(math.log10(eta.max()))
    return [10.0**e for e in range(low, high + 1) if e % 2 == 0]


def draw(curves: pd.DataFrame, medians: pd.DataFrame) -> Figure:
    style.apply()
    fig = style.figure("figure1")
    axes = fig.subplots(1, 3)

    centred = np.stack(
        [
            curves[curves.criterion == c]
            .pivot(index="rep", columns="k", values="centred")
            .to_numpy()
            for c in CRITERIA
        ]
    )
    rows = [medians[medians.criterion == c].sort_values("k") for c in CRITERIA]
    eta = rows[0].eta.to_numpy()
    median = np.stack([r["median"].to_numpy() for r in rows])
    start = float(rows[0].eta[rows[0].shading_start].iloc[0])

    reach = np.abs(centred).max(axis=(1, 2))
    ticks = (
        style.round_tick(reach[0]),
        *2 * [style.round_tick(reach[1:].max())],
    )
    for j, axis in enumerate(axes):
        axis.set_box_aspect(1)
        axis.set_xscale("log")
        axis.set_xlim(eta[0] / 1.25, eta[-1] * 1.25)
        axis.axvspan(start, eta[-1] * 1.25, color=style.color("shade"), zorder=0)
        axis.plot(
            eta,
            centred[j].T,
            color=style.color("black"),
            alpha=0.07,
            linewidth=0.24,
            zorder=1,
        )
        axis.plot(
            eta,
            median[j],
            color=style.color("median"),
            linewidth=style.line("main"),
            zorder=4,
        )
        axis.axvline(
            float(rows[j].eta[rows[j].median_minimum].iloc[0]),
            color=style.color("selection"),
            linestyle="--",
            linewidth=style.line("secondary"),
            zorder=3,
        )
        axis.set_xticks(x_ticks(eta))
        axis.xaxis.set_major_formatter(LogFormatterMathtext())
        axis.minorticks_off()
        style.set_symmetric_y(axis, ticks[j])
        if ticks[j] >= 1:
            axis.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
        axis.set_xlabel(r"$\eta$")
        axis.set_ylabel(LABELS[CRITERIA[j]])
        style.panel_label(axis, f"({'abc'[j]})")
    return fig


def make(folder) -> Figure:
    curves, medians = (pd.read_csv(folder / name) for name in FILES)
    return draw(curves, medians)
