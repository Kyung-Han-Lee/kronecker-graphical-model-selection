from __future__ import annotations

from pathlib import Path

from . import style as fig_style
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "../outputs/adni/figures"
FIT_DIR = ROOT / "../outputs/adni/selected_fits"
ATLAS_NII = ROOT / "deps/atlas/ROI_MNI_V4.nii"
ATLAS_TXT = ROOT / "deps/atlas/ROI_MNI_V4.txt"
BRAINNET_MESH = ROOT / "deps/brainnet-viewer/BrainMesh_ICBM152.nv"
GROUPS = ("CN", "MCI", "Dementia")
GAMMA_TAG = "0p5"
SUPPORT_TOL = 1e-6


PARTIAL_CORRELATION_THRESHOLD = 0.20
C_NODE = fig_style.color("node")
C_EDGE = fig_style.color("edge")
PANEL_LABELS = ("a", "b", "c")

fig_style.apply()


def partial_correlation(omega: np.ndarray) -> np.ndarray:
    omega = np.asarray(omega, dtype=float)
    if omega.ndim != 2 or omega.shape[0] != omega.shape[1]:
        raise ValueError("omega must be square")
    if not np.allclose(omega, omega.T, atol=1e-10):
        raise ValueError("omega must be symmetric")
    diagonal = np.diag(omega)
    if np.any(diagonal <= 0):
        raise ValueError("omega must have a positive diagonal")
    scale = np.sqrt(diagonal)
    out = -omega / np.outer(scale, scale)
    np.fill_diagonal(out, 0.0)
    return out


def edges_above(matrix: np.ndarray, threshold: float) -> np.ndarray:
    if threshold < 0:
        raise ValueError("threshold must be nonnegative")
    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("matrix must be square")
    out = np.where(np.abs(matrix) > threshold, matrix, 0.0)
    np.fill_diagonal(out, 0.0)
    return out


def connectome_panel_title(index: int, group: str) -> str:
    return f"({PANEL_LABELS[index]}) {group}"


def brainnet_surface_mesh(
    path: Path | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    path = BRAINNET_MESH if path is None else Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"missing BrainNet Viewer surface: {path}; see fig/scripts/README.md"
        )
    with path.open() as handle:
        n_vertices = int(handle.readline())
        vertices = np.loadtxt(handle, max_rows=n_vertices)
        n_faces = int(handle.readline())
        faces = np.loadtxt(handle, dtype=int, max_rows=n_faces) - 1
    if vertices.shape != (n_vertices, 3):
        raise ValueError(f"unexpected BrainNet vertex array: {vertices.shape}")
    if faces.shape != (n_faces, 3):
        raise ValueError(f"unexpected BrainNet face array: {faces.shape}")
    if faces.min() < 0 or faces.max() >= n_vertices:
        raise ValueError("BrainNet face indices fall outside the vertex array")
    return vertices, faces


def aal_centroids(atlas_nii=None, atlas_labels=None) -> np.ndarray:
    rows = [
        line.split("\t")
        for line in (ATLAS_TXT if atlas_labels is None else Path(atlas_labels))
        .read_text()
        .splitlines()
    ]
    codes = [int(row[2]) for row in rows]
    if codes != sorted(codes):
        raise RuntimeError("AAL label file is not ordered by label code")
    image = nib.load(str(ATLAS_NII if atlas_nii is None else atlas_nii))
    volume = np.asarray(image.dataobj)
    present = sorted(int(value) for value in np.unique(volume) if value != 0)
    if present != codes:
        raise RuntimeError("AAL image and label file contain different codes")
    coordinates = np.zeros((len(codes), 3))
    for index, code in enumerate(codes):
        voxel = np.argwhere(volume == code).mean(axis=0)
        coordinates[index] = nib.affines.apply_affine(image.affine, voxel)
    return coordinates


def load_joint_group(group: str, fit_dir=None) -> tuple[np.ndarray, pd.Series]:
    directory = FIT_DIR if fit_dir is None else Path(fit_dir)
    omega_path = directory / f"Omega1_{group}_projected_joint_g{GAMMA_TAG}.csv"
    summary_path = directory / "selected_fits.csv"
    if not omega_path.exists():
        raise FileNotFoundError(f"missing joint factor: {omega_path}")
    summary = pd.read_csv(summary_path)
    row = summary[
        (summary["group"] == group)
        & (summary["criterion"] == "joint")
        & np.isclose(summary["gamma"], 0.5)
    ]
    if len(row) != 1:
        raise ValueError(
            f"expected one joint summary row for {group}, found {len(row)}"
        )
    omega = pd.read_csv(omega_path).to_numpy(dtype=float)
    if omega.shape != (116, 116):
        raise ValueError(f"unexpected ROI factor shape for {group}: {omega.shape}")
    edge_count = np.count_nonzero(np.triu(np.abs(omega) > SUPPORT_TOL, k=1))
    expected_edges = int(row.iloc[0]["matrix_edges_roi"])
    if edge_count != expected_edges:
        raise ValueError(
            f"joint ROI edge count for {group} is {edge_count}, "
            f"but the summary records {expected_edges}"
        )
    return omega, row.iloc[0]


def save_figure(figure: plt.Figure, stem: str) -> None:
    fig_style.save(figure, FIG_DIR, stem)
