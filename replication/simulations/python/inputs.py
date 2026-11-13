from pathlib import Path

from replication.io import read_frame


def load_paths(folder):
    frame = read_frame(
        Path(folder) / "paths.csv",
        required=(
            "id",
            "experiment",
            "M",
            "dimensions",
            "family",
            "degree",
            "n",
            "rep",
        ),
    )
    if frame["id"].duplicated().any():
        raise ValueError("duplicate dataset ids")
    return frame


def paper_scope(paths, dimensions=(40, 60, 80), sample_sizes=(100, 200, 800, 3200)):
    parts = paths["dimensions"].astype(str).str.split("x")
    p1 = parts.str[0].astype(int)
    p2 = parts.str[1].astype(int)
    return paths[
        (paths["experiment"] == "matrix")
        & (paths["M"] == 2)
        & (paths["family"].str.lower() == "er")
        & p1.isin(dimensions)
        & (p1 == p2)
        & paths["n"].isin(sample_sizes)
    ].copy()
