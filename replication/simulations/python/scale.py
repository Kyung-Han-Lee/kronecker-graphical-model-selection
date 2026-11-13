import math
from statistics import median

from replication.io import finite
from .inputs import paper_scope


def summarize(
    paths, scores, scales, dimensions=(40, 60, 80), sample_sizes=(100, 200, 800, 3200)
):
    scope = paper_scope(paths, dimensions, sample_sizes)
    selected = scores[
        (scores["criterion"] == "scale")
        & scores["id"].isin(scope["id"])
        & (scores["status"] == "ok")
    ].copy()
    if selected.empty:
        raise ValueError("no scale scores in the requested scope")
    finite(selected, ("ebic", "gamma", "k"))
    selected = selected.sort_values(["ebic", "k"]).drop_duplicates(["id", "gamma"])
    metadata = scope.set_index("id")
    scale_rows = scales.set_index(["id", "k"])
    margins = []
    ratios = {}
    for row in selected.itertuples(index=False):
        design = metadata.loc[row.id]
        value = scale_rows.loc[(row.id, row.k)]
        n = int(design.n)
        p = math.prod(int(x) for x in design.dimensions.split("x"))
        penalty = math.log(n) + 4 * float(row.gamma) * math.log(p)
        h_plus_p = (float(row.ebic) - penalty) / (n * p) + math.log(p) - 1
        margin = float(value.B) - h_plus_p
        if not math.isfinite(margin):
            raise ValueError(f"nonfinite certificate margin: {row.id}")
        margins.append(
            {
                "id": str(row.id),
                "gamma": float(row.gamma),
                "k": int(row.k),
                "margin": margin,
            }
        )
        if float(row.gamma) == 0.5:
            key = (design.dimensions, int(design.degree), n)
            ratios.setdefault(key, []).append(
                float(value.tau) / float(value.truth_scale)
            )
    if set(selected["id"]) != set(scope["id"]):
        raise ValueError("missing selected scale paths")
    by_n = {}
    for key, values in ratios.items():
        by_n.setdefault(key[2], []).append(median(values))
    if not by_n:
        raise ValueError("no gamma=0.5 scale ratios")
    return {
        "pairs": len(margins),
        "positive": sum(row["margin"] > 0 for row in margins),
        "minimum": min(margins, key=lambda row: row["margin"]),
        "tau_medians": {str(n): median(values) for n, values in sorted(by_n.items())},
    }
