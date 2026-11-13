from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import numpy as np
import pandas as pd

from . import style as fig_style
from replication.io import output_directory

ROOT = Path(__file__).resolve().parents[1]

FIG_DIR = ROOT / "../outputs/adni/figures"
FIG_STEM = "FigureC1"
RESULT_DIR = ROOT / "../outputs/adni/joint_profile"
PATH_CSV = RESULT_DIR / "probe_path.csv"
VERDICT_CSV = RESULT_DIR / "probe_verdict.csv"

GAMMA = 0.5
GROUPS = ("CN", "MCI", "Dementia")
ARM = "projected"

XLIM_PAD_HI = 1.15
XLIM_PAD_LO = 0.8

C_BEST = fig_style.color("selection")
C_SHADE = fig_style.color("shade")
GROUP_STYLE = {g: dict(color=fig_style.color(f"group.{g}")) for g in GROUPS}
PANELS = (
    ("ebic_plug", "plug", "(a)", r"$\mathrm{EBIC}_{\gamma}$"),
    (
        "ebic_joint_uncounted",
        "joint_uncounted",
        "(b)",
        r"$\mathrm{EBIC}^{\mathrm{relax}}_{\gamma}$",
    ),
)


PANEL_ROW_SHARE = 0.77


PANEL_ROW_SHIFT = 0.023

Y_LABEL_PAD = 2

fig_style.apply()


def load_frames(path_csv=None, verdict_csv=None) -> tuple[pd.DataFrame, pd.DataFrame]:
    d = pd.read_csv(PATH_CSV if path_csv is None else path_csv)
    d = d[(d.arm == ARM) & np.isclose(d.gamma, GAMMA)]
    d = d.sort_values(["group", "eta_index"])
    v = pd.read_csv(VERDICT_CSV if verdict_csv is None else verdict_csv)
    v = v[(v.arm == ARM) & np.isclose(v.gamma, GAMMA)]

    missing = set(GROUPS) - set(d.group)
    if missing:
        raise RuntimeError(f"경로 CSV 에 군이 없다: {sorted(missing)}")
    dup = d.groupby(["group", "eta_index"]).size()
    if (dup != 1).any():
        raise RuntimeError(
            f"경로 CSV 에 (군, eta_index) 중복이 있다 (arm={ARM}, "
            f"gamma={GAMMA}): {dup[dup != 1].to_dict()}"
        )
    vdup = v.groupby(["group", "criterion"]).size()
    if (vdup != 1).any():
        raise RuntimeError(
            f"판정 CSV 에 (군, 기준) 중복이 있다: {vdup[vdup != 1].to_dict()}"
        )
    return d, v


def empty_region_start(d: pd.DataFrame) -> tuple[float, dict[str, float]]:
    first_empty = {}
    for g in GROUPS:
        empty = d[(d.group == g) & (d.df == 0)]
        if empty.empty:
            raise RuntimeError(
                f"{g} 는 격자 안에서 빈 그래프에 닿지 않는다. 음영 구간을 "
                f"정의할 수 없으니 격자를 넓히거나 그림 설계를 고쳐라."
            )
        first_empty[g] = float(empty.eta.min())

    start = max(first_empty.values())
    live = d[(d.eta >= start) & (d.df != 0)]
    if not live.empty:
        raise RuntimeError(
            f"eta >= {start:.4g} 인데 비지 않은 격자점이 있다 "
            f"({len(live)} 개). 음영이 캡션 주장과 다르다."
        )
    return start, first_empty


def center_criterion(raw: np.ndarray) -> np.ndarray:
    return raw - raw.mean()


def plot_curves(d, v) -> None:
    etas = np.sort(d.eta.unique())
    empty_from, first_empty = empty_region_start(d)
    fig, axes = plt.subplots(
        1,
        len(PANELS),
        figsize=(fig_style.width(), fig_style.height("adni_ebic_curve")),
        layout="constrained",
    )
    fig.get_layout_engine().set(
        rect=((1 - PANEL_ROW_SHARE) / 2 + PANEL_ROW_SHIFT, 0, PANEL_ROW_SHARE, 1)
    )
    reach = []
    for ax, (col, crit, title, ylabel) in zip(axes, PANELS):
        ax.set_box_aspect(1)
        ax.axvspan(empty_from, etas[-1] * XLIM_PAD_HI, color=C_SHADE, zorder=0)
        sel_etas = []
        panel_reach = 0.0
        for g in GROUPS:
            sub = d[d.group == g]
            raw = sub[col].values
            if not np.isfinite(raw).all():
                bad = sub.eta.values[~np.isfinite(raw)]
                raise RuntimeError(
                    f"{g}/{crit}: 기준값이 유한하지 않은 격자점이 있다 (eta={np.array2string(bad, precision=4)})"
                )
            y = center_criterion(raw)
            panel_reach = max(panel_reach, float(np.abs(y).max()))
            ax.plot(
                sub.eta.values, y, lw=fig_style.line("main"), label=g, **GROUP_STYLE[g]
            )
            k = int(np.argmin(y))
            want = v[(v.group == g) & (v.criterion == crit)].best_eta.iloc[0]
            got = sub.eta.values[k]
            if not np.isclose(got, want):
                raise RuntimeError(
                    f"{g}/{crit}: 그림 argmin eta={got:.4g} 이 판정 {want:.4g} 와 다르다"
                )
            sel_etas.append(got)
            ax.plot(
                got,
                y[k],
                "o",
                ms=fig_style.STYLE["marker_pt"],
                color=GROUP_STYLE[g]["color"],
            )
        if len(set(np.round(sel_etas, 10))) == 1:
            ax.axvline(
                sel_etas[0], color=C_BEST, ls="--", lw=fig_style.line("secondary")
            )
        ax.set_xscale("log")
        ax.set_xlim(etas[0] * XLIM_PAD_LO, etas[-1] * XLIM_PAD_HI)
        ax.set_xticks([0.001, 0.1, 10.0])
        fig_style.panel_label(ax, title)
        ax.set_xlabel("$\\eta$")
        ax.set_ylabel(ylabel, labelpad=Y_LABEL_PAD)
        ax.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
        reach.append(panel_reach)
    for ax, panel_reach in zip(axes, reach):
        fig_style.set_symmetric_y(ax, fig_style.round_tick(panel_reach))
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=len(GROUPS))
    return fig


def draw_curves(path_csv, verdict_csv, output):
    with output_directory(output) as directory:
        paths, verdicts = load_frames(path_csv, verdict_csv)
        figure = plot_curves(paths, verdicts)
        try:
            fig_style.save(figure, directory, FIG_STEM)
        finally:
            plt.close(figure)
    return Path(output)
