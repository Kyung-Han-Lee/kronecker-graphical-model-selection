from replication.io import finite
from .inputs import paper_scope


def summarize(
    paths, covariances, dimensions=(40, 60, 80), sample_sizes=(100, 200, 800, 3200)
):
    scope = paper_scope(paths, dimensions, sample_sizes)
    rows = covariances[covariances["id"].isin(scope["id"])]
    if rows.empty or set(rows["id"]) != set(scope["id"]):
        raise ValueError("missing covariance paths")
    finite(rows, ("minimum_eigenvalue",))
    return {
        "paths": len(scope),
        "grid_points": len(rows[["id", "k"]].drop_duplicates()),
        "minimum_eigenvalue": float(rows["minimum_eigenvalue"].min()),
        "nonpositive_values": int((rows["minimum_eigenvalue"] <= 0).sum()),
    }
