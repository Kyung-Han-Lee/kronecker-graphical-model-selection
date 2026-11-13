from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from . import figure1
from . import figure2
from . import figure3
from . import style

FIGURES = {"figure1": figure1, "figure2": figure2, "figure3": figure3}


def draw_all(results: Path) -> dict[str, str]:
    folder = results / "figures"
    print(f"text font: {style.apply()}")
    skipped = {}
    for name, module in FIGURES.items():
        missing = [f for f in module.FILES if not (folder / f).is_file()]
        if missing:
            skipped[name] = "no " + ", ".join(missing)
            print(f"{name} skipped: {skipped[name]}")
            continue
        figure = module.make(folder)
        path = style.save(figure, folder / f"{name.capitalize()}.pdf")
        plt.close(figure)
        print(f"{name}: {path}")
    return skipped


def main() -> None:
    default = Path(__file__).resolve().parents[2] / "outputs" / "results"
    parser = argparse.ArgumentParser(description="Draw simulation Figures 1-3.")
    parser.add_argument("--results", type=Path, default=default)
    draw_all(parser.parse_args().results)


if __name__ == "__main__":
    main()
