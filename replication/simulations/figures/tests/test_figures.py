from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection, PolyCollection


from replication.simulations.figures import draw_figures
from replication.simulations.figures import figure1
from replication.simulations.figures import figure2
from replication.simulations.figures import figure3
from replication.simulations.figures import style


def figure1_values() -> tuple[pd.DataFrame, pd.DataFrame]:
    eta = 10.0 ** np.array([-4, -2, 0, 2])
    curves, medians = [], []
    for criterion, values in (
        ("ebic", [3.0, 1.0, -1.0, -3.0]),
        ("scale", [1.0, -3.0, -1.0, 3.0]),
        ("relaxation", [4.0, 2.0, -4.0, -2.0]),
    ):
        for rep in (1, 2):
            for k in range(4):
                curves.append(
                    {
                        "criterion": criterion,
                        "rep": rep,
                        "k": k + 1,
                        "eta": eta[k],
                        "centred": values[k] * rep / 2,
                    }
                )
        best = int(np.argmin(values))
        for k in range(4):
            medians.append(
                {
                    "criterion": criterion,
                    "k": k + 1,
                    "eta": eta[k],
                    "median": values[k] * 0.75,
                    "median_minimum": k == best,
                    "shading_start": k == 1,
                }
            )
    return pd.DataFrame(curves), pd.DataFrame(medians)


def figure2_values() -> pd.DataFrame:
    rows = []
    for panel, top, count in (("a", 60.0, 100), ("b", 6.0, 120)):
        lam = np.concatenate(
            [np.linspace(0, 0.882, 40), np.linspace(0.918, top, count - 40)]
        )
        for i, value in enumerate(lam):
            rows.append(
                {
                    "panel": panel,
                    "point": i + 1,
                    "lambda": value,
                    "theta": 26.0 + 19.0 * value / top,
                    "z": min(0.0, value - 0.9) / 2,
                    "log_det": -1.7 + value / top,
                    "lambda_max": 0.9,
                }
            )
    return pd.DataFrame(rows)


def figure3_values(families=("er", "ba")) -> pd.DataFrame:
    rng = np.random.default_rng(1)
    rows = []
    for family in families:
        for n in (10, 50):
            for method in ("scale", "relaxation", "cv", "stars", "stars_01"):
                values = rng.uniform(0.2, 0.9, 5)
                if method == "cv" and n == 10:
                    values[:] = 0.4
                for rep, mcc in enumerate(values, 1):
                    rows.append(
                        {
                            "family": family,
                            "n": n,
                            "method": method,
                            "rep": rep,
                            "mcc": mcc,
                        }
                    )
    return pd.DataFrame(rows)


class StyleTest(unittest.TestCase):
    def test_round_tick(self):
        self.assertEqual(style.round_tick(1000), 1000)
        self.assertEqual(style.round_tick(1200), 1500)
        self.assertAlmostEqual(style.round_tick(0.3), 0.3)

    def test_colours_by_role(self):
        self.assertEqual(style.color("median"), "#D55E00")
        self.assertEqual(style.color("shade"), "#E8E8E8")
        self.assertEqual(style.color("#123456"), "#123456")

    def test_text_family(self):
        family = style.apply()
        if style.font_files():
            self.assertEqual(family, style.STYLE["font"]["family"])
        else:
            self.assertEqual(family, style.STYLE["font"]["fallback"])


class Figure1Test(unittest.TestCase):
    def test_draw(self):
        curves, medians = figure1_values()
        fig = figure1.draw(curves, medians)
        self.assertEqual(len(fig.axes), 3)
        for axis, criterion, minimum in zip(
            fig.axes, ("ebic", "scale", "relaxation"), (100.0, 0.01, 1.0)
        ):
            lines = axis.get_lines()
            self.assertEqual(len(lines), 4)
            np.testing.assert_array_equal(
                lines[2].get_ydata(),
                medians[medians.criterion == criterion]["median"],
            )
            self.assertEqual(lines[3].get_xdata()[0], minimum)
            low, high = axis.get_ylim()
            self.assertAlmostEqual(low, -high)
            self.assertEqual(list(axis.get_xticks()), [1e-4, 1e-2, 1.0, 100.0])

            self.assertEqual(axis.patches[0].get_x(), 0.01)

    def test_profiled_panels_share_ticks(self):
        curves, medians = figure1_values()
        fig = figure1.draw(curves, medians)
        ticks = [tuple(axis.get_yticks()) for axis in fig.axes]
        self.assertEqual(ticks[0], (-3.0, 0.0, 3.0))
        self.assertEqual(ticks[1], (-4.0, 0.0, 4.0))
        self.assertEqual(ticks[2], ticks[1])


class Figure2Test(unittest.TestCase):
    def test_draw(self):
        path = figure2_values()
        fig = figure2.draw(path)
        axis_a, axis_b = fig.axes[:2]
        a = path[path.panel == "a"]

        line = axis_a.get_lines()[1]
        np.testing.assert_array_equal(line.get_xdata(), a.theta)
        np.testing.assert_array_equal(line.get_ydata(), a.z)

        offsets = [c.get_offsets()[0] for c in axis_a.collections[-3:]]
        np.testing.assert_allclose(offsets[0], a[["theta", "z"]].iloc[0])
        np.testing.assert_allclose(offsets[1], a[["theta", "z"]].iloc[40])
        np.testing.assert_allclose(offsets[2], [45.0, 0.0])
        b = path[path.panel == "b"]
        np.testing.assert_array_equal(axis_b.get_lines()[0].get_ydata(), b.log_det)
        self.assertEqual(axis_b.get_xticks()[1], 0.9)


class Figure3Test(unittest.TestCase):
    def test_draw(self):
        mcc = figure3_values()
        fig = figure3.draw(mcc)
        axes = [a for a in fig.axes]
        self.assertEqual([a.get_title() for a in axes], ["ER", "BA"])
        labels = [t.get_text() for t in fig.legends[0].get_texts()]
        self.assertEqual(
            labels, ["Scale", "Relaxation", "CV", "StARS (0.05)", "StARS (0.1)"]
        )
        means = mcc[mcc.family == "ba"].groupby(["method", "n"]).mcc.mean()
        for line, method in zip(
            axes[1].get_lines(), ("scale", "relaxation", "cv", "stars", "stars_01")
        ):
            np.testing.assert_allclose(line.get_ydata(), means[method].to_numpy())

        kinds = [type(c) for c in axes[0].collections]
        self.assertEqual(sum(issubclass(k, LineCollection) for k in kinds), 1)
        self.assertEqual(sum(issubclass(k, PolyCollection) for k in kinds), 9)

    def test_one_family(self):
        fig = figure3.draw(figure3_values(("er",)))
        self.assertEqual([a.get_title() for a in fig.axes], ["ER"])

    def test_missing_method(self):
        mcc = figure3_values()
        with self.assertRaisesRegex(ValueError, "no MCC"):
            figure3.draw(
                mcc[~((mcc.method == "cv") & (mcc.n == 50) & (mcc.family == "er"))]
            )


class DrawFiguresTest(unittest.TestCase):
    def test_draw_all(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp) / "figures"
            folder.mkdir()
            curves, medians = figure1_values()
            curves.to_csv(folder / "figure1_curves.csv", index=False)
            medians.to_csv(folder / "figure1_medians.csv", index=False)
            figure3_values().to_csv(folder / "figure3_mcc.csv", index=False)
            with contextlib.redirect_stdout(io.StringIO()):
                skipped = draw_figures.draw_all(Path(tmp))
            self.assertEqual(skipped, {"figure2": "no figure2_path.csv"})
            self.assertEqual(
                sorted(p.name for p in folder.glob("*.pdf")),
                ["Figure1.pdf", "Figure3.pdf"],
            )

            first = (folder / "Figure3.pdf").read_bytes()
            with contextlib.redirect_stdout(io.StringIO()):
                draw_figures.draw_all(Path(tmp))
            self.assertEqual((folder / "Figure3.pdf").read_bytes(), first)


if __name__ == "__main__":
    unittest.main()
