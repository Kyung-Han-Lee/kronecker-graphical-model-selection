import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from replication.adni.figures import connectome as fig3


def _toy_scene():

    vertices = np.array(
        [
            [-40.0, -40.0, -30.0],
            [40.0, -40.0, -30.0],
            [0.0, 40.0, -30.0],
            [0.0, 0.0, 40.0],
        ]
    )
    faces = np.array([[0, 1, 2], [0, 1, 3], [1, 2, 3], [0, 2, 3]])
    coordinates = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [0.0, 10.0, 5.0]])
    shown = np.zeros((3, 3))
    shown[0, 1] = shown[1, 0] = -0.4
    return vertices, faces, coordinates, shown


def test_generate_creates_a_new_output_directory(tmp_path, monkeypatch):
    vertices, faces, coordinates, shown = _toy_scene()
    scene = vertices, faces, coordinates, {group: shown for group in fig3.common.GROUPS}
    monkeypatch.setattr(fig3, "load_scene", lambda *args: scene)
    monkeypatch.setattr(fig3, "plot_connectome", lambda **kwargs: plt.figure())
    output = tmp_path / "network"
    fig3.generate(output)
    assert (output / "Figure4.pdf").exists()
    assert (output / "Figure4.png").exists()
    assert len(list((output / "brainnet").glob("*.edge"))) == 3
    with pytest.raises(FileExistsError):
        fig3.generate(output)


def test_node_radii_grow_linearly_with_degree():
    degree = np.array([0, 1, 4])
    got = fig3.node_radii(degree)

    assert got[0] == fig3.NODE_RADIUS_BASE
    assert np.allclose(np.diff(got) / np.diff(degree), fig3.NODE_RADIUS_PER_DEGREE)
    assert np.all(np.diff(got) > 0)


def test_edge_pairs_list_each_undirected_edge_once():
    shown = np.zeros((4, 4))
    for i, j in ((0, 2), (1, 3), (2, 3)):
        shown[i, j] = shown[j, i] = 0.5
    got = fig3.edge_pairs(shown)

    assert got.tolist() == [[0, 2], [1, 3], [2, 3]]
    assert np.all(got[:, 0] < got[:, 1])


def test_union_bbox_covers_every_mask_and_clips_the_margin():
    first = np.zeros((10, 12), dtype=bool)
    first[2:4, 3:6] = True
    second = np.zeros((10, 12), dtype=bool)
    second[7:9, 8:11] = True

    assert fig3.union_bbox([first, second], margin=1) == (1, 10, 2, 12)
    assert fig3.union_bbox([first, second], margin=5) == (0, 10, 0, 12)


def test_views_are_left_lateral_then_dorsal():
    names = [view[0] for view in fig3.VIEWS]
    assert names == ["left-lateral", "dorsal"]
    lateral_eye, lateral_up = fig3.VIEWS[0][1], fig3.VIEWS[0][2]
    dorsal_eye, dorsal_up = fig3.VIEWS[1][1], fig3.VIEWS[1][2]

    assert lateral_eye == (-1.0, 0.0, 0.0) and lateral_up == (0.0, 0.0, 1.0)

    assert dorsal_eye == (0.0, 0.0, 1.0) and dorsal_up == (0.0, 1.0, 0.0)


def test_plot_connectome_places_two_rendered_views_in_each_group_column(
    monkeypatch,
):
    omega = np.eye(116)
    omega[0, 1] = omega[1, 0] = 0.3

    omega[2, 2] = omega[3, 3] = 1e-6
    omega[2, 3] = omega[3, 2] = 5e-7
    calls = []

    monkeypatch.setattr(fig3.common, "aal_centroids", lambda: np.zeros((116, 3)))
    monkeypatch.setattr(
        fig3.common,
        "brainnet_surface_mesh",
        lambda path=None: (np.zeros((4, 3)), np.array([[0, 1, 2], [0, 1, 3]])),
    )
    monkeypatch.setattr(
        fig3.common,
        "load_joint_group",
        lambda group: (omega, pd.Series({"edges_roi": 1, "n": 10})),
    )

    def fake_render(
        vertices,
        faces,
        coordinates,
        shown,
        view,
        size,
        *,
        common_mask,
        additional_color,
    ):
        assert len(fig3.edge_pairs(shown)) == 1
        assert shown[2, 3] == 0
        assert common_mask[0, 1] and not common_mask[2, 3]
        assert additional_color in fig3.ADDITIONAL_COLORS.values()
        calls.append(view[0])
        image = np.zeros((size[1], size[0], 4), dtype=np.uint8)
        image[10:30, 10:50] = (200, 190, 180, 255)
        return image

    figure = fig3.plot_connectome(render=fake_render)
    try:
        assert calls == ["left-lateral", "dorsal"] * 3
        positions = [
            (
                axis.get_subplotspec().rowspan.start,
                axis.get_subplotspec().colspan.start,
            )
            for axis in figure.axes
        ]
        assert positions == [(0, 0), (1, 0), (0, 1), (1, 1), (0, 2), (1, 2)]
        assert [figure.axes[index].get_title() for index in (0, 2, 4)] == [
            "CN",
            "MCI",
            "Dementia",
        ]
        labels = [
            text.get_text() for index in (0, 2, 4) for text in figure.axes[index].texts
        ]
        assert labels == ["(a)", "(b)", "(c)"]
        assert all(not axis.axison for axis in figure.axes)
    finally:
        plt.close(figure)


def test_embedded_panels_retain_600_ppi_at_article_width(tmp_path):
    groups = fig3.common.GROUPS
    image = np.full((200, 1600, 4), 255, dtype=np.uint8)
    panels = {(group, view[0]): image for group in groups for view in fig3.VIEWS}
    figure = fig3.compose_figure(panels, groups)
    try:
        pixels_per_panel = figure.axes[0].images[0].get_array().shape[1]
        article_width = 6.4
        assert pixels_per_panel / (article_width / len(groups)) >= 600
        pdf = next(
            path
            for path in fig3.fig_style.save(figure, tmp_path, "probe")
            if path.suffix == ".pdf"
        )
        images = subprocess.run(
            ["pdfimages", "-list", str(pdf)],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.splitlines()[2:]
        native_ppi = [
            int(row.split()[12]) for row in images if row.split()[2] == "image"
        ]
        assert native_ppi and min(native_ppi) * fig3.FIGURE_WIDTH / article_width >= 600
    finally:
        plt.close(figure)


def test_brainnet_files_follow_node_and_edge_layout(tmp_path):
    _, _, coordinates, shown = _toy_scene()
    node_path, edge_path = fig3.write_brainnet_files("CN", coordinates, shown, tmp_path)

    node = np.loadtxt(node_path, dtype=str)
    assert node.shape == (3, 6)
    assert np.allclose(node[:, :3].astype(float), coordinates)
    assert node[:, 5].tolist() == ["-", "-", "-"]
    edge = np.loadtxt(edge_path)
    assert edge.shape == (3, 3)
    assert np.allclose(edge, shown)


def test_composite_veils_only_the_network_behind_the_surface():
    surface = np.zeros((1, 4, 4), dtype=np.uint8)
    surface[..., :3] = 200
    surface[0, :3, 3] = 255
    surface_depth = np.array([[-500.0, -500.0, -500.0, np.nan]])
    network = np.zeros((1, 4, 4), dtype=np.uint8)
    network[0, [0, 1, 3], 3] = 255
    network_depth = np.array([[-400.0, -600.0, np.nan, -450.0]])

    got = fig3.composite_layers(
        surface,
        surface_depth,
        [(network, network_depth, 1.0)],
        behind_alpha=0.7,
    )

    assert got[0, 0].tolist() == [0, 0, 0, 255]
    assert got[0, 1].tolist() == [
        60,
        60,
        60,
        255,
    ]
    assert got[0, 2].tolist() == [200, 200, 200, 255]
    assert got[0, 3].tolist() == [0, 0, 0, 255]


def test_composite_layers_orders_fragments_by_depth_at_each_pixel():

    surface = np.zeros((1, 4, 4), dtype=np.uint8)
    surface_depth = np.full((1, 4), np.nan)
    nodes = np.zeros((1, 4, 4), dtype=np.uint8)
    nodes[0, :3, 3] = 255
    node_depth = np.array([[-450.0, -450.0, -450.0, np.nan]])
    edges = np.zeros((1, 4, 4), dtype=np.uint8)
    edges[..., 2] = 255
    edges[0, [0, 1, 3], 3] = 255
    edge_depth = np.array([[-470.0, -430.0, np.nan, -440.0]])
    layers = [(edges, edge_depth, 0.6), (nodes, node_depth, 1.0)]

    got = fig3.composite_layers(surface, surface_depth, layers)

    assert got[0, 0].tolist() == [0, 0, 0, 255]

    assert got[0, 1].tolist() == [0, 0, 153, 255]

    assert got[0, 2].tolist() == [0, 0, 0, 255]
    assert got[0, 3].tolist() == [0, 0, 255, 153]

    np.testing.assert_array_equal(
        got, fig3.composite_layers(surface, surface_depth, layers[::-1])
    )


def test_composite_layers_veils_fragments_behind_the_surface_toward_its_colour():
    surface = np.zeros((1, 3, 4), dtype=np.uint8)
    surface[..., :3] = 200
    surface[..., 3] = 255
    surface_depth = np.full((1, 3), -500.0)
    layer = np.zeros((1, 3, 4), dtype=np.uint8)
    layer[0, :2, 3] = 255
    depth = np.array([[-600.0, -400.0, np.nan]])

    opaque = fig3.composite_layers(
        surface, surface_depth, [(layer, depth, 1.0)], behind_alpha=0.85
    )

    assert opaque[0, 0].tolist() == [30, 30, 30, 255]
    assert opaque[0, 1].tolist() == [0, 0, 0, 255]
    assert opaque[0, 2].tolist() == [200, 200, 200, 255]
    translucent = fig3.composite_layers(
        surface, surface_depth, [(layer, depth, 0.7)], behind_alpha=0.85
    )

    assert translucent[0, 0].tolist() == [81, 81, 81, 255]
    assert translucent[0, 1].tolist() == [60, 60, 60, 255]


def test_surface_is_white_opaque_and_outlined_with_a_partial_veil():
    assert fig3.C_SURFACE == "#FFFFFF" and fig3.SURFACE_OPACITY == 1.0
    assert fig3.C_OUTLINE and fig3.OUTLINE_WIDTH > 0
    assert 0.0 < fig3.BEHIND_ALPHA < 1.0
    assert fig3.FIG_STEM == "Figure4"


def test_render_view_returns_rgba_with_transparent_background():
    pytest.importorskip("pyvista")
    vertices, faces, coordinates, shown = _toy_scene()
    image = fig3.render_view(
        vertices, faces, coordinates, shown, fig3.VIEWS[0], size=(160, 120)
    )

    assert image.shape == (120, 160, 4)
    assert image.dtype == np.uint8
    assert image[..., 3].max() == 255
    assert image[0, 0, 3] == 0


def test_edge_radii_grow_linearly_with_absolute_partial_correlation():
    got = fig3.edge_radii(np.array([0.0, -0.1, 0.5]))

    assert got[0] == fig3.EDGE_RADIUS_BASE
    assert np.isclose(got[1], fig3.EDGE_RADIUS_BASE + 0.1 * fig3.EDGE_RADIUS_SLOPE)
    assert np.isclose(got[2], fig3.EDGE_RADIUS_BASE + 0.5 * fig3.EDGE_RADIUS_SLOPE)


def test_edge_tubes_carry_one_radius_per_edge():
    pytest.importorskip("pyvista")
    _, _, coordinates, shown = _toy_scene()
    shown[0, 2] = shown[2, 0] = 0.6
    tubes = fig3.edge_tubes(coordinates, shown)

    assert tubes is not None and tubes.n_points > 0
    assert fig3.edge_tubes(coordinates, np.zeros((3, 3))) is None


def test_panels_span_the_central_panel_scale_of_the_text_width():

    assert fig3.PANEL_SCALE == 1.0
    groups = fig3.common.GROUPS
    image = np.full((300, 400, 4), 255, dtype=np.uint8)
    panels = {(group, view[0]): image for group in groups for view in fig3.VIEWS}
    figure = fig3.compose_figure(panels, groups)
    try:
        width, height = figure.get_size_inches()
        assert width == pytest.approx(fig3.FIGURE_WIDTH)
        boxes = [axis.get_position() for axis in figure.axes]
        margin = (1 - fig3.PANEL_SCALE) / 2
        assert min(box.x0 for box in boxes) == pytest.approx(margin)
        assert max(box.x1 for box in boxes) == pytest.approx(1 - margin)
        rows = [
            axis.images[0].get_array().shape[0]
            for axis in figure.axes[: len(fig3.VIEWS)]
        ]
        widest = max(axis.images[0].get_array().shape[1] for axis in figure.axes)
        expected = (
            sum(rows) * fig3.PANEL_SCALE * fig3.FIGURE_WIDTH / len(groups) / widest
        )
        assert height == pytest.approx(expected + fig3.TITLE_HEIGHT, rel=1e-6)
    finally:
        plt.close(figure)


def test_common_mask_excludes_edges_shared_by_only_two_groups():
    shown = {group: np.zeros((4, 4)) for group in fig3.common.GROUPS}
    for group, weight in zip(fig3.common.GROUPS, (0.25, -0.3, 0.4), strict=True):
        shown[group][0, 1] = shown[group][1, 0] = weight
    for group in ("CN", "MCI"):
        shown[group][0, 2] = shown[group][2, 0] = 0.27
    shown["Dementia"][2, 3] = shown["Dementia"][3, 2] = -0.5
    mask = fig3.common_edge_mask(shown)
    assert fig3.edge_pairs(mask).tolist() == [[0, 1]]
    for group, weights in shown.items():
        shared = np.where(mask, weights, 0)
        added = np.where(~mask, weights, 0)
        np.testing.assert_array_equal(shared + added, weights)
        assert len(fig3.edge_pairs(added)) == 1
    assert np.where(mask, shown["MCI"], 0)[0, 1] == -0.3


def test_common_and_group_colours_are_distinct_without_a_legend():
    groups = fig3.common.GROUPS
    image = np.full((300, 400, 4), 255, dtype=np.uint8)
    panels = {(group, view[0]): image for group in groups for view in fig3.VIEWS}
    figure = fig3.compose_figure(panels, groups)
    try:
        assert len({fig3.COMMON_COLOR, *fig3.ADDITIONAL_COLORS.values()}) == 4
        assert not figure.legends
    finally:
        plt.close(figure)


def test_render_view_accepts_common_and_additional_edge_layers():
    pytest.importorskip("pyvista")
    vertices, faces, coordinates, shown = _toy_scene()
    shown[0, 2] = shown[2, 0] = 0.5
    mask = np.zeros_like(shown, dtype=bool)
    mask[0, 1] = mask[1, 0] = True
    image = fig3.render_view(
        vertices,
        faces,
        coordinates,
        shown,
        fig3.VIEWS[0],
        size=(160, 120),
        common_mask=mask,
        additional_color=fig3.ADDITIONAL_COLORS["CN"],
    )
    assert image.shape == (120, 160, 4)
    visible = image[..., 3] > 0
    rgb = image[..., :3].astype(float)
    assert np.any(visible & (rgb[..., 2] > rgb[..., 0] + 30))


def test_edge_radii_scale_multiplies_the_radius():
    weights = np.array([0.2, -0.45, 0.7])
    np.testing.assert_allclose(
        fig3.edge_radii(weights, scale=4 / 3), 4 / 3 * fig3.edge_radii(weights)
    )


def test_edge_tubes_scale_the_tube_radius():
    pytest.importorskip("pyvista")
    _, _, coordinates, shown = _toy_scene()

    def width(scale):
        bounds = fig3.edge_tubes(coordinates, shown, scale=scale).bounds
        return bounds[3] - bounds[2]

    assert width(4 / 3) == pytest.approx(4 / 3 * width(1.0), rel=1e-6)


def test_render_view_draws_only_the_coloured_layer_wider(monkeypatch):
    pytest.importorskip("pyvista")
    vertices, faces, coordinates, shown = _toy_scene()
    shown[0, 2] = shown[2, 0] = 0.5
    mask = np.zeros_like(shown, dtype=bool)
    mask[0, 1] = mask[1, 0] = True
    original = fig3.edge_tubes
    scales = {}

    def recording(coordinates, layer, scale=1.0):
        pairs = tuple(tuple(pair) for pair in fig3.edge_pairs(layer).tolist())
        scales[pairs] = scale
        return original(coordinates, layer, scale=scale)

    monkeypatch.setattr(fig3, "edge_tubes", recording)
    fig3.render_view(
        vertices,
        faces,
        coordinates,
        shown,
        fig3.VIEWS[0],
        size=(80, 60),
        common_mask=mask,
        additional_color=fig3.ADDITIONAL_COLORS["CN"],
    )
    assert scales == {
        ((0, 1),): fig3.COMMON_EDGE_SCALE,
        ((0, 2),): fig3.ADDITIONAL_EDGE_SCALE,
    }
    assert fig3.ADDITIONAL_EDGE_SCALE == pytest.approx(4 / 3)


def test_render_view_composites_common_additional_and_node_layers(monkeypatch):
    pytest.importorskip("pyvista")
    vertices, faces, coordinates, shown = _toy_scene()
    shown[0, 2] = shown[2, 0] = 0.5
    mask = np.zeros_like(shown, dtype=bool)
    mask[0, 1] = mask[1, 0] = True
    original = fig3.composite_layers
    recorded = []

    def recording(surface, surface_depth, layers, behind_alpha=fig3.BEHIND_ALPHA):
        recorded.append([opacity for _, _, opacity in layers])
        return original(surface, surface_depth, layers, behind_alpha=behind_alpha)

    monkeypatch.setattr(fig3, "composite_layers", recording)
    fig3.render_view(
        vertices,
        faces,
        coordinates,
        shown,
        fig3.VIEWS[0],
        size=(80, 60),
        common_mask=mask,
        additional_color=fig3.ADDITIONAL_COLORS["CN"],
    )
    fig3.render_view(vertices, faces, coordinates, shown, fig3.VIEWS[0], size=(80, 60))

    assert recorded == [[fig3.COMMON_OPACITY, 1.0, 1.0], [1.0, 1.0]]


def test_common_edges_are_dark_grey_translucent_and_thinner():

    assert fig3.COMMON_COLOR == "#1A1A1A"
    assert fig3.COMMON_OPACITY == 0.7
    assert fig3.COMMON_EDGE_SCALE == 0.8
    assert fig3.BEHIND_ALPHA == 0.85
