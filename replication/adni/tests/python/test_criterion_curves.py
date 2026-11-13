import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from replication.adni.figures import criterion_curves as curve


@pytest.fixture
def frames():
    rows, verdicts = [], []
    eta = np.array([0.001, 0.01, 0.1, 1.0, 10.0])
    for index, group in enumerate(curve.GROUPS):
        plug = 1e6 + index * 1000 + np.array([20, 15, 10, 5, 0]) * 1e5
        relax = 1e6 + index * 1000 + np.array([20, 10, 0, 5, 15]) * 1e4
        for point in range(len(eta)):
            rows.append(
                {
                    "group": group,
                    "eta": eta[point],
                    "df": 0 if point == 4 else 2,
                    "ebic_plug": plug[point],
                    "ebic_joint_uncounted": relax[point],
                }
            )
        verdicts.extend(
            [
                {"group": group, "criterion": "plug", "best_eta": eta[-1]},
                {"group": group, "criterion": "joint_uncounted", "best_eta": eta[2]},
            ]
        )
    path, verdict = pd.DataFrame(rows), pd.DataFrame(verdicts)
    return path, verdict


def test_curve_has_only_ebic_and_relaxation_panels():
    assert [(p[0], p[1]) for p in curve.PANELS] == [
        ("ebic_plug", "plug"),
        ("ebic_joint_uncounted", "joint_uncounted"),
    ]


def test_two_panels_plot_the_original_centered_values(frames):
    path, verdict = frames
    figure = curve.plot_curves(path, verdict)
    try:
        assert len(figure.axes) == 2
        for axis, column in zip(figure.axes, ("ebic_plug", "ebic_joint_uncounted")):
            for group in curve.GROUPS:
                line = next(line for line in axis.lines if line.get_label() == group)
                values = path[path.group == group]
                np.testing.assert_array_equal(line.get_xdata(), values.eta)
                np.testing.assert_allclose(
                    line.get_ydata(),
                    values[column] - values[column].mean(),
                )
    finally:
        plt.close(figure)


def test_panel_b_label_is_separate_and_panels_are_square(frames):
    path, verdict = frames
    figure = curve.plot_curves(path, verdict)
    try:
        assert len(figure.axes) == 2
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        width, height = figure.get_size_inches()
        left, right = figure.axes
        for axis in figure.axes:
            box = axis.get_position()
            assert box.width * width == pytest.approx(box.height * height, abs=1e-3)
        label = right.yaxis.label.get_window_extent(renderer)
        frame = left.get_window_extent(renderer).x1
        assert (label.x0 - frame) / figure.dpi > 0.15
        tight = figure.get_tightbbox(renderer)
        assert tight.x0 >= 0 and tight.x1 <= width
        assert tight.y0 >= 0 and tight.y1 <= height
    finally:
        plt.close(figure)


def test_draw_curves_writes_submission_names(frames, tmp_path):
    path, verdict = frames
    path = path.assign(arm="projected", gamma=0.5)
    path["eta_index"] = path.groupby("group").cumcount() + 1
    verdict = verdict.assign(arm="projected", gamma=0.5)
    path.to_csv(tmp_path / "path.csv", index=False)
    verdict.to_csv(tmp_path / "verdict.csv", index=False)
    output = tmp_path / "figures"
    curve.draw_curves(tmp_path / "path.csv", tmp_path / "verdict.csv", output)
    assert sorted(p.name for p in output.iterdir()) == ["FigureC1.pdf", "FigureC1.png"]
