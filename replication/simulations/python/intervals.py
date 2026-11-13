import pandas as pd


def summarize(paths, intervals):
    paths = pd.DataFrame(paths)
    intervals = pd.DataFrame(intervals).copy()
    paths = paths[paths["experiment"].isin(("matrix", "tensor"))]
    if paths.empty:
        raise ValueError("no matrix or tensor paths")
    intervals = intervals[intervals["id"].isin(paths["id"])]
    expected = paths.set_index("id")["M"].astype(int)
    violations = 0
    first = set()
    for name, rows in intervals.groupby("id", observed=True, sort=False):
        k = rows["k"].astype(int)
        q = rows["q"].astype(int)
        if k.duplicated().any() or set(k) != set(range(1, int(k.max()) + 1)):
            raise ValueError(f"incomplete interval path: {name}")
        violations += bool((q[k >= 2] != expected.loc[name]).any())
        first.update(q[k == 1].tolist())
    if set(intervals["id"]) != set(paths["id"]):
        raise ValueError("missing interval paths")
    return {
        "paths": len(paths),
        "violations": violations,
        "q_at_first_grid_point": sorted(first),
    }
