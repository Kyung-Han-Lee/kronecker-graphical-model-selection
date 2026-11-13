from replication.io import finite
from .inputs import paper_scope


def summarize(
    paths, conditions, dimensions=(40, 60, 80), sample_sizes=(100, 200, 800, 3200)
):
    scope = paper_scope(paths, dimensions, sample_sizes)
    rows = conditions[conditions["id"].isin(scope["id"])]
    if rows.empty or set(rows["id"]) != set(scope["id"]):
        raise ValueError("missing condition-number paths")
    finite(rows, ("cond_base", "cond_relaxed"))
    maxima = {}
    for mode, group in rows.groupby("mode", sort=True):
        for column, label in (
            ("cond_base", "cond_base"),
            ("cond_relaxed", "cond_joint"),
        ):
            row = group.loc[group[column].idxmax()]
            maxima[label + str(int(mode))] = {
                "value": float(row[column]),
                "id": str(row["id"]),
                "k": int(row["k"]),
            }
    return {
        "cells": len(scope),
        "grid_points": len(rows[["id", "k"]].drop_duplicates()),
        "maxima": maxima,
    }
