from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure
from matplotlib.ticker import NullLocator

from . import style

FILES = ("figure3_mcc.csv",)
FAMILIES = {"er": "ER", "ba": "BA"}
METHODS = ("scale", "relaxation", "cv")


COLORS = {
    "scale": "#800080",
    "relaxation": "#0000FF",
    "cv": "#CC0000",
    "stars": "#008000",
    "stars_01": "#FFD700",
}
MARKERS = {"scale": "s", "relaxation": "o", "cv": "D", "stars": "^", "stars_01": "v"}
LINES = {"scale": "--", "relaxation": "-", "cv": ":", "stars": "-.", "stars_01": "--"}
LABELS = {
    "scale": "Scale",
    "relaxation": "Relaxation",
    "cv": "CV",
    "stars": "StARS (0.05)",
    "stars_01": "StARS (0.1)",
}


STARS_ZORDER = 1.5
VIOLIN_WIDTH = 0.45


CUT = 2.0
BANDWIDTH_ADJUST = 1.25
OUTLINE_POINTS = 400


FILL_ALPHA = 0.65
LEGEND_BAND_IN = 0.22


def violin_outline(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    bandwidth = BANDWIDTH_ADJUST * values.std(ddof=1) * values.size**-0.2
    grid = np.linspace(
        values.min() - CUT * bandwidth,
        min(values.max() + CUT * bandwidth, 1.0),
        OUTLINE_POINTS,
    )

    density = np.exp(-0.5 * ((grid[:, None] - values[None, :]) / bandwidth) ** 2).sum(
        axis=1
    )
    return grid, density / density.max() * VIOLIN_WIDTH / 2


def draw(mcc: pd.DataFrame) -> Figure:
    families = [f for f in FAMILIES if f in set(mcc.family)]
    sizes = sorted(mcc.n.unique())
    methods = [*METHODS, "stars", "stars_01"]
    arrays = {
        (f, n, m): g.sort_values("rep").mcc.to_numpy(dtype=float)
        for (f, n, m), g in mcc.groupby(["family", "n", "method"])
    }
    missing = [
        (f, n, m)
        for f in families
        for n in sizes
        for m in methods
        if (f, n, m) not in arrays
    ]
    if missing:
        raise ValueError(f"Figure 3: no MCC of {missing[0]}")
    style.apply()
    figure, axes = plt.subplots(
        len(families),
        1,
        sharex=True,
        sharey=True,
        squeeze=False,
        figsize=(style.width(), style.height("figure3")),
    )
    outlines = {key: violin_outline(v) for key, v in arrays.items() if np.ptp(v) > 0}

    lowest = min(
        [grid[0] for grid, _ in outlines.values()]
        + [v[0] for key, v in arrays.items() if key not in outlines]
    )
    lower = np.floor(lowest / 0.05) * 0.05
    positions = np.arange(len(sizes), dtype=float)
    for panel, family in enumerate(families):
        axis = axes[panel, 0]
        for method in methods:
            color = COLORS[method]
            zorder = STARS_ZORDER if method.startswith("stars") else 2
            for i, n in enumerate(sizes):
                key = (family, n, method)
                if key not in outlines:
                    axis.hlines(
                        arrays[key][0],
                        positions[i] - VIOLIN_WIDTH / 2,
                        positions[i] + VIOLIN_WIDTH / 2,
                        color=color,
                        linewidth=style.line("secondary"),
                        zorder=zorder,
                    )
                    continue
                grid, half = outlines[key]
                axis.fill_betweenx(
                    grid,
                    positions[i] - half,
                    positions[i] + half,
                    facecolor=to_rgba(color, FILL_ALPHA),
                    linewidth=0,
                    zorder=zorder,
                )
            axis.plot(
                positions,
                [float(np.mean(arrays[family, n, method])) for n in sizes],
                color=color,
                linestyle=LINES[method],
                marker=MARKERS[method],
                markerfacecolor="white",
                markeredgewidth=style.line("secondary"),
                linewidth=style.line("secondary"),
                markersize=style.STYLE["marker_pt"],
                label=LABELS[method],
                zorder=3,
            )
        axis.set_xlim(-0.55, len(sizes) - 0.45)
        axis.set_ylim(lower, 1.02)
        axis.set_yticks(np.arange(np.ceil(lower / 0.25) * 0.25, 1.01, 0.25))
        axis.set_axisbelow(True)
        axis.set_ylabel("MCC")
        axis.xaxis.set_minor_locator(NullLocator())
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(
            axis="y",
            color=style.color("shade"),
            linewidth=style.line("reference"),
            zorder=0,
        )
        style.panel_label(axis, f"({'ab'[panel]})")
        axis.set_title(FAMILIES[family])
    axes[-1, 0].set_xticks(range(len(sizes)), labels=[str(n) for n in sizes])
    axes[-1, 0].set_xlabel(r"$n$")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="lower center",
        ncol=len(handles),
        bbox_to_anchor=(0.5, -0.01),
    )

    figure.tight_layout(
        rect=(0, LEGEND_BAND_IN / figure.get_figheight(), 1, 1),
        pad=0.4,
        h_pad=0.5,
    )
    return figure


def make(folder) -> Figure:
    return draw(pd.read_csv(folder / FILES[0]))
