from __future__ import annotations

from pathlib import Path
from replication.io import output_directory

from . import connectome_inputs as common
from . import style as fig_style
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from PIL import Image

FIG_STEM = "Figure4"
BRAINNET_DIR = common.ROOT / "../outputs/adni/brainnet"

COMMON_COLOR = "#1A1A1A"
COMMON_OPACITY = 0.7
ADDITIONAL_COLORS = {
    "CN": "#0000FF",
    "MCI": "#FFD700",
    "Dementia": "#FF0000",
}


C_SURFACE = "#FFFFFF"
C_OUTLINE = "#B8B8B8"

OUTLINE_WIDTH = 2.0
SURFACE_OPACITY = 1.0

BEHIND_ALPHA = 0.85

SURFACE_AMBIENT = 0.12
SURFACE_DIFFUSE = 1.0
SURFACE_SPECULAR = 0.2
SURFACE_SPECULAR_POWER = 30

SSAO_RADIUS = 8.0
SSAO_BIAS = 0.3
SSAO_KERNEL_SIZE = 256

NODE_RADIUS_BASE = 1.8
NODE_RADIUS_PER_DEGREE = 0.25


EDGE_RADIUS_BASE = 0.45
EDGE_RADIUS_SLOPE = 1.5

COMMON_EDGE_SCALE = 0.8
ADDITIONAL_EDGE_SCALE = 4 / 3

EDGE_AMBIENT = 0.9
EDGE_DIFFUSE = 0.2


RENDER_SIZE = (2000, 2000)
SUPERSAMPLE = 2

CROP_MARGIN = 24
FIGURE_WIDTH = fig_style.width()
TITLE_HEIGHT = 0.26


PANEL_SCALE = 1.0

CAMERA_DISTANCE = 600.0


VIEWS = (
    ("left-lateral", (-1.0, 0.0, 0.0), (0.0, 0.0, 1.0)),
    ("dorsal", (0.0, 0.0, 1.0), (0.0, 1.0, 0.0)),
)


def node_radii(degree: np.ndarray) -> np.ndarray:
    return NODE_RADIUS_BASE + NODE_RADIUS_PER_DEGREE * np.asarray(degree, dtype=float)


def edge_radii(weights: np.ndarray, scale: float = 1.0) -> np.ndarray:
    return scale * (
        EDGE_RADIUS_BASE + EDGE_RADIUS_SLOPE * np.abs(np.asarray(weights, dtype=float))
    )


def edge_tubes(coordinates: np.ndarray, shown: np.ndarray, scale: float = 1.0):
    import pyvista as pv

    pairs = edge_pairs(shown)
    if len(pairs) == 0:
        return None

    points = coordinates[pairs].reshape(-1, 3)
    lines = np.hstack(
        [np.full((len(pairs), 1), 2), np.arange(2 * len(pairs)).reshape(-1, 2)]
    )
    poly = pv.PolyData(points, lines=lines.ravel())
    poly["radius"] = np.repeat(edge_radii(shown[pairs[:, 0], pairs[:, 1]], scale), 2)
    return poly.tube(scalars="radius", absolute=True, n_sides=12)


def edge_pairs(shown: np.ndarray) -> np.ndarray:
    return np.transpose(np.nonzero(np.triu(shown, k=1)))


def common_edge_mask(shown: dict[str, np.ndarray]) -> np.ndarray:
    return np.logical_and.reduce([weights != 0 for weights in shown.values()])


def union_bbox(masks: list[np.ndarray], margin: int) -> tuple[int, int, int, int]:
    stacked = np.any(np.stack(masks), axis=0)
    rows = np.flatnonzero(stacked.any(axis=1))
    cols = np.flatnonzero(stacked.any(axis=0))
    if rows.size == 0 or cols.size == 0:
        raise ValueError("cannot crop: every rendered panel is empty")
    height, width = stacked.shape
    return (
        max(int(rows[0]) - margin, 0),
        min(int(rows[-1]) + 1 + margin, height),
        max(int(cols[0]) - margin, 0),
        min(int(cols[-1]) + 1 + margin, width),
    )


def crop(image: np.ndarray, box: tuple[int, int, int, int]) -> np.ndarray:
    top, bottom, left, right = box
    return image[top:bottom, left:right]


def _distance(depth: np.ndarray) -> np.ndarray:
    depth = np.asarray(depth, dtype=np.float32)
    return np.where(np.isfinite(depth), np.abs(depth), np.float32(np.inf))


def composite_layers(
    surface: np.ndarray,
    surface_depth: np.ndarray,
    layers: list[tuple[np.ndarray, np.ndarray, float]],
    behind_alpha: float = BEHIND_ALPHA,
) -> np.ndarray:
    surface = np.asarray(surface, dtype=np.float32)
    surface_rgb = surface[..., :3]
    surface_distance = _distance(surface_depth)
    prepared = []
    for image, depth, opacity in layers:
        image = np.asarray(image, dtype=np.float32)
        distance = _distance(depth)
        behind = (np.isfinite(distance) & (distance > surface_distance))[..., None]
        rgb = np.where(
            behind,
            behind_alpha * image[..., :3] + (1.0 - behind_alpha) * surface_rgb,
            image[..., :3],
        )
        prepared.append((rgb, image[..., 3:4] / 255.0 * opacity, distance))

    order = np.argsort(
        np.stack([-distance for _, _, distance in prepared]),
        axis=0,
        kind="stable",
    )
    alpha = surface[..., 3:4] / 255.0
    rgb = surface_rgb * alpha
    for rank in range(len(prepared)):
        chosen = order[rank]
        layer_rgb = np.zeros_like(rgb)
        layer_alpha = np.zeros_like(alpha)
        for index, (candidate_rgb, candidate_alpha, _) in enumerate(prepared):
            mask = chosen == index
            layer_rgb[mask] = candidate_rgb[mask]
            layer_alpha[mask] = candidate_alpha[mask]
        rgb = layer_rgb * layer_alpha + rgb * (1.0 - layer_alpha)
        alpha = layer_alpha + alpha * (1.0 - layer_alpha)
    rgb = np.where(alpha > 0, rgb / np.maximum(alpha, 1e-9), 255.0)
    out = np.concatenate([rgb, alpha * 255.0], axis=-1)
    return np.clip(np.rint(out), 0, 255).astype(np.uint8)


def _render_pass(add_actors, vertices, view, size, ssao: bool):
    import pyvista as pv

    pv.OFF_SCREEN = True
    _, eye, up = view
    width, height = size
    plotter = pv.Plotter(
        off_screen=True,
        window_size=(width * SUPERSAMPLE, height * SUPERSAMPLE),
        lighting="light_kit",
    )
    plotter.set_background("white")
    if ssao:
        plotter.enable_ssao(
            radius=SSAO_RADIUS,
            bias=SSAO_BIAS,
            kernel_size=SSAO_KERNEL_SIZE,
            blur=True,
        )
    add_actors(plotter)

    lower, upper = vertices.min(axis=0), vertices.max(axis=0)
    centre = (lower + upper) / 2
    camera = plotter.camera
    camera.focal_point = tuple(centre)
    camera.position = tuple(centre + CAMERA_DISTANCE * np.asarray(eye))
    camera.up = up
    camera.parallel_projection = True
    camera.parallel_scale = 0.55 * float((upper - lower).max())
    image = np.asarray(
        plotter.screenshot(None, transparent_background=True, return_img=True)
    )
    depth = np.asarray(plotter.get_image_depth(fill_value=np.nan), dtype=float)
    plotter.close()
    return image, depth


def render_view(
    vertices: np.ndarray,
    faces: np.ndarray,
    coordinates: np.ndarray,
    shown: np.ndarray,
    view: tuple[str, tuple[float, float, float], tuple[float, float, float]],
    size: tuple[int, int] = RENDER_SIZE,
    *,
    common_mask: np.ndarray | None = None,
    additional_color: str = common.C_EDGE,
) -> np.ndarray:
    try:
        import pyvista as pv
    except ImportError as error:
        raise ImportError(
            "Figure 4 needs pyvista with a headless VTK build; "
            "see the module docstring for the install lines"
        ) from error

    surface = pv.PolyData(
        vertices, np.hstack([np.full((len(faces), 1), 3), faces]).ravel()
    )
    degree = np.count_nonzero(shown, axis=1)
    points = pv.PolyData(coordinates)
    points["radius"] = node_radii(degree)
    spheres = points.glyph(
        geom=pv.Sphere(radius=1.0, theta_resolution=24, phi_resolution=24),
        scale="radius",
        orient=False,
    )
    if common_mask is None:
        edge_layers = [(edge_tubes(coordinates, shown), common.C_EDGE, 1.0)]
    else:
        edge_layers = [
            (
                edge_tubes(
                    coordinates,
                    np.where(common_mask, shown, 0),
                    scale=COMMON_EDGE_SCALE,
                ),
                COMMON_COLOR,
                COMMON_OPACITY,
            ),
            (
                edge_tubes(
                    coordinates,
                    np.where(~common_mask, shown, 0),
                    scale=ADDITIONAL_EDGE_SCALE,
                ),
                additional_color,
                1.0,
            ),
        ]

    def add_surface(plotter):
        plotter.add_mesh(
            surface,
            color=C_SURFACE,
            opacity=SURFACE_OPACITY,
            smooth_shading=True,
            ambient=SURFACE_AMBIENT,
            diffuse=SURFACE_DIFFUSE,
            specular=SURFACE_SPECULAR,
            specular_power=SURFACE_SPECULAR_POWER,
        )
        plotter.add_silhouette(surface, color=C_OUTLINE, line_width=OUTLINE_WIDTH)

    def add_tubes(tubes, color):
        def add(plotter):
            plotter.add_mesh(
                tubes,
                color=color,
                smooth_shading=True,
                ambient=EDGE_AMBIENT,
                diffuse=EDGE_DIFFUSE,
                specular=0.2,
            )

        return add

    def add_nodes(plotter):
        plotter.add_mesh(
            spheres, color=common.C_NODE, smooth_shading=True, specular=0.3
        )

    surface_image, surface_depth = _render_pass(
        add_surface, vertices, view, size, ssao=True
    )

    passes = [
        (add_tubes(tubes, color), opacity)
        for tubes, color, opacity in edge_layers
        if tubes is not None
    ] + [(add_nodes, 1.0)]
    layers = []
    for add_actors, opacity in passes:
        image, depth = _render_pass(add_actors, vertices, view, size, ssao=False)
        layers.append((image, depth, opacity))
    composed = composite_layers(surface_image, surface_depth, layers)
    width, height = size
    shrunk = Image.fromarray(composed).resize((width, height), Image.LANCZOS)
    return np.asarray(shrunk)


def load_scene(
    groups: tuple[str, ...],
    fits=None,
    atlas_nii=None,
    atlas_labels=None,
    surface=None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    vertices, faces = (
        common.brainnet_surface_mesh(surface)
        if surface is not None
        else common.brainnet_surface_mesh()
    )
    coordinates = (
        common.aal_centroids(atlas_nii, atlas_labels)
        if atlas_nii is not None or atlas_labels is not None
        else common.aal_centroids()
    )
    shown = {}
    for group in groups:
        omega, _ = (
            common.load_joint_group(group, fits)
            if fits is not None
            else common.load_joint_group(group)
        )
        graph = common.edges_above(omega, common.SUPPORT_TOL) != 0
        shown[group] = common.edges_above(
            np.where(graph, common.partial_correlation(omega), 0.0),
            common.PARTIAL_CORRELATION_THRESHOLD,
        )
    return vertices, faces, coordinates, shown


def compose_figure(panels: dict, groups: tuple[str, ...]) -> plt.Figure:
    cropped = {}
    for name, _, _ in VIEWS:
        box = union_bbox(
            [panels[group, name][..., 3] > 0 for group in groups], CROP_MARGIN
        )
        for group in groups:
            cropped[group, name] = crop(panels[group, name], box)

    raster = fig_style.STYLE["raster"]
    target_pdf_ppi = fig_style.raster_pdf_dpi()
    widest = max(cropped[groups[0], name].shape[1] for name, _, _ in VIEWS)
    scale = target_pdf_ppi * FIGURE_WIDTH / len(groups) / widest
    if scale < 0.999:
        for key, panel in list(cropped.items()):
            size = (
                max(1, round(panel.shape[1] * scale)),
                max(1, round(panel.shape[0] * scale)),
            )
            cropped[key] = np.asarray(
                Image.fromarray(panel).resize(size, Image.LANCZOS)
            )
    heights = [cropped[groups[0], name].shape[0] for name, _, _ in VIEWS]
    widest = max(cropped[groups[0], name].shape[1] for name, _, _ in VIEWS)
    inch_per_px = PANEL_SCALE * FIGURE_WIDTH / len(groups) / widest
    raster_dpi = 1.0 / inch_per_px
    effective_ppi = (
        widest * len(groups) / (PANEL_SCALE * raster["widest_inclusion_width_in"])
    )
    print(
        f"raster panels at {raster_dpi:.0f} ppi in PDF, "
        f"{effective_ppi:.0f} ppi at {raster['widest_inclusion_width_in']} in"
    )
    if effective_ppi < raster["min_effective_ppi"]:
        print(
            f"warning: below the {raster['min_effective_ppi']} ppi target; raise RENDER_SIZE"
        )
    figure_height = sum(heights) * inch_per_px + TITLE_HEIGHT
    figure = plt.figure(figsize=(FIGURE_WIDTH, figure_height))
    grid = GridSpec(
        len(VIEWS),
        len(groups),
        figure=figure,
        height_ratios=heights,
        hspace=0.0,
        wspace=0.0,
        left=(1.0 - PANEL_SCALE) / 2,
        right=1.0 - (1.0 - PANEL_SCALE) / 2,
        top=1.0 - TITLE_HEIGHT / figure_height,
        bottom=0.0,
    )
    for column, group in enumerate(groups):
        for row, (name, _, _) in enumerate(VIEWS):
            axis = figure.add_subplot(grid[row, column])
            axis.imshow(cropped[group, name])
            axis.set_axis_off()
            if row == 0:
                fig_style.panel_label(axis, f"({common.PANEL_LABELS[column]})")
                axis.set_title(group)
    return figure


def plot_connectome(
    groups: tuple[str, ...] = common.GROUPS, render=render_view, scene=None
) -> plt.Figure:
    vertices, faces, coordinates, shown = scene or load_scene(groups)
    common_mask = common_edge_mask(shown)
    panels = {
        (group, view[0]): render(
            vertices,
            faces,
            coordinates,
            shown[group],
            view,
            RENDER_SIZE,
            common_mask=common_mask,
            additional_color=ADDITIONAL_COLORS[group],
        )
        for group in groups
        for view in VIEWS
    }
    return compose_figure(panels, groups)


def write_brainnet_files(
    group: str, coordinates: np.ndarray, shown: np.ndarray, out_dir: Path
) -> tuple[Path, Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    node_path = out_dir / f"{group}.node"
    edge_path = out_dir / f"{group}.edge"
    radii = node_radii(np.count_nonzero(shown, axis=1))
    with node_path.open("w") as handle:
        for (x, y, z), radius in zip(coordinates, radii):
            handle.write(f"{x:.3f}\t{y:.3f}\t{z:.3f}\t1\t{radius:.2f}\t-\n")
    np.savetxt(edge_path, shown, fmt="%.6g", delimiter="\t")
    return node_path, edge_path


def generate(output, fits=None, atlas_nii=None, atlas_labels=None, surface=None):
    output = Path(output)
    with output_directory(output) as directory:
        scene = load_scene(common.GROUPS, fits, atlas_nii, atlas_labels, surface)
        figure = plot_connectome(scene=scene)
        try:
            fig_style.save(figure, directory, FIG_STEM)
        finally:
            plt.close(figure)
        for group in common.GROUPS:
            write_brainnet_files(
                group, scene[2], scene[3][group], directory / "brainnet"
            )
    return output
