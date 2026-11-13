from __future__ import annotations

import json
import math
from pathlib import Path

from matplotlib import font_manager

import matplotlib.pyplot as plt

STYLE_PATH = Path(__file__).resolve().with_name("style.json")
STYLE = json.loads(STYLE_PATH.read_text())


def register_fonts() -> str:
    family = STYLE["font"]["family"]
    otf_dir = Path(STYLE["font"]["otf_dir"])
    for name in STYLE["font"]["otf_files"]:
        path = otf_dir / name
        if path.exists():
            font_manager.fontManager.addfont(str(path))
    return family


def apply() -> None:
    family = register_fonts()
    font = STYLE["font"]
    lines = STYLE["line_pt"]
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": [family, "DejaVu Serif"],
            "font.size": font["size_pt"],
            "axes.labelsize": font["size_pt"],
            "axes.titlesize": font["title_pt"],
            "xtick.labelsize": font["tick_pt"],
            "ytick.labelsize": font["tick_pt"],
            "legend.fontsize": font["legend_pt"],
            "legend.frameon": False,
            "mathtext.fontset": font["mathtext_fontset"],
            "mathtext.fallback": font["mathtext_fallback"],
            "mathtext.default": "it",
            "lines.linewidth": lines["main"],
            "lines.markersize": STYLE["marker_pt"],
            "axes.linewidth": lines["axes"],
            "xtick.major.width": lines["axes"],
            "ytick.major.width": lines["axes"],
            "savefig.dpi": STYLE["dpi"],
            "pdf.fonttype": 3,
            "ps.fonttype": 3,
        }
    )


def width() -> float:
    return float(STYLE["text_width_in"])


def height(figure_key: str) -> float:
    return float(STYLE["heights_in"][figure_key])


def raster_pdf_dpi() -> int:
    raster = STYLE["raster"]
    return math.ceil(
        raster["min_effective_ppi"] * raster["widest_inclusion_width_in"] / width()
    )


def line(kind: str) -> float:
    return float(STYLE["line_pt"][kind])


def color(role: str) -> str:
    value = STYLE["roles"].get(role, role)
    return STYLE["palette"].get(value, value)


def figure(figure_key: str, **kwargs) -> plt.Figure:
    kwargs.setdefault("layout", "constrained")
    return plt.figure(figsize=(width(), height(figure_key)), **kwargs)


def panel_label(axis, text: str):
    spec = STYLE["panel_label"]
    return axis.annotate(
        text,
        xy=(spec["x_axes"], spec["y_axes"]),
        xycoords="axes fraction",
        xytext=tuple(spec["offset_pt"]),
        textcoords="offset points",
        ha=spec["ha"],
        va=spec["va"],
        fontsize=STYLE["font"]["panel_label_pt"],
        fontweight=STYLE["font"]["panel_label_weight"],
    )


def round_tick(reach: float) -> float:
    spec = STYLE["y_ticks"]
    target = reach / spec["limit_factor"]
    scale = 10.0 ** math.floor(math.log10(target))
    return next(step * scale for step in spec["round_steps"] if step * scale >= target)


def set_symmetric_y(axis, tick: float) -> None:
    factor = STYLE["y_ticks"]["limit_factor"]
    axis.set_ylim(-factor * tick, factor * tick)
    axis.set_yticks([-tick, 0, tick])


def save(fig: plt.Figure, fig_dir: Path, stem: str) -> list[Path]:
    written = []
    for extension in ("pdf", "png"):
        path = Path(fig_dir) / f"{stem}.{extension}"
        dpi = raster_pdf_dpi() if extension == "pdf" else STYLE["dpi"]
        fig.savefig(path, dpi=dpi)
        written.append(path)
        print(f"saved {path}")
    return written
