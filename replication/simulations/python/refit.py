from collections import defaultdict
from dataclasses import asdict, dataclass
import math
from pathlib import Path
from statistics import mean, median

import numpy as np

from replication import io
from .inputs import load_paths
from .likelihood import compare_likelihoods, focus_indices, order_disagreement

METHODS = ("scale", "relaxation", "refit")


class AnalysisError(ValueError):
    pass


@dataclass(frozen=True)
class PathInfo:
    id: str
    family: str
    M: int
    dimensions: tuple[int, ...]
    degree: int
    n: int
    rep: int


@dataclass(frozen=True)
class Candidate:
    k: int
    eta: float
    df: int
    df_edges: int
    loglik: float
    ebic: float
    mcc: float


def summarize_comparisons(paths, manifests, supports, candidates) -> None:
    manifests = {row["id"]: row for row in manifests}
    supports = {(row["id"], int(row["k"])): row for row in supports}
    if not paths or len(supports) != len(paths) * 40:
        raise AnalysisError("incomplete analysis inputs")
    groups = defaultdict(list)
    for row in paths:
        chosen = supports[row["id"], int(row["k_relaxation"])]
        other = supports[row["id"], int(row["k_refit"])]
        row["same_support"] = str(chosen["pattern_hash"] == other["pattern_hash"])
        row["selected_exact_zero_match"] = str(int(chosen["tiny_edges"]) == 0)
        row["unique_supports"] = str(
            len({supports[row["id"], k]["pattern_hash"] for k in range(1, 41)})
        )
        groups[
            tuple((row[key] for key in ("family", "M", "dimensions", "degree", "n")))
        ].append(row)
    summaries = []
    metrics = (
        "selected_gain_fraction",
        "selected_gap_scale",
        "selected_gap_relaxation",
        "refit_score_regret",
        "refit_score_regret_normalized",
        "mcc_difference",
        "all_rank_disagreement",
        "focus_rank_disagreement",
        "all_gain_fraction_median",
        "focus_gain_fraction_median",
        "unique_supports",
    )
    for key, group in groups.items():
        if len({int(r["rep"]) for r in group}) != len(group):
            raise AnalysisError(f"incomplete setting: {key}")
        entry = dict(
            zip(("family", "M", "dimensions", "degree", "n"), key, strict=True)
        )
        entry["replications"] = len(group)
        for flag in (
            "same_index",
            "same_support",
            "selected_exact_zero_match",
            "selected_refit_below_relaxation",
        ):
            entry[flag + "_count"] = sum((r[flag] == "True" for r in group))
        for metric in metrics:
            values = [float(row[metric]) for row in group if row[metric] != ""]
            entry[metric + "_count"] = len(values)
            entry[metric + "_mean"] = mean(values) if values else None
            for q, name in (
                (0, "min"),
                (0.05, "q05"),
                (0.5, "median"),
                (0.95, "q95"),
                (1, "max"),
            ):
                entry[metric + "_" + name] = (
                    float(np.quantile(values, q)) if values else None
                )
        summaries.append(entry)
    candidate_groups = defaultdict(list)
    negative = []
    for row in candidates:
        support = supports[row["id"], int(row["k"])]
        if row["refit_below_relaxation"] == "True":
            negative.append(
                {
                    **row,
                    "tiny_edges": support["tiny_edges"],
                    "max_tiny": support["max_tiny"],
                }
            )
        key = (row["family"], row["M"], row["dimensions"], row["degree"], row["n"])
        candidate_groups[key].append((row, int(support["tiny_edges"]) == 0))
    candidate_summaries = []
    for key, group in candidate_groups.items():
        for region in ("all", "focus", "selected_relaxation"):
            for mask in ("all", "exact_zero_match"):
                selected = [
                    row
                    for row, exact in group
                    if (region == "all" or row[region] == "True")
                    and (mask == "all" or exact)
                ]
                entry = dict(
                    zip(("family", "M", "dimensions", "degree", "n"), key, strict=True)
                )
                entry.update(region=region, mask=mask, candidates=len(selected))
                entry["refit_below_relaxation"] = sum(
                    (r["refit_below_relaxation"] == "True" for r in selected)
                )
                entry["relaxation_below_scale"] = sum(
                    (r["relaxation_below_scale"] == "True" for r in selected)
                )
                for metric in ("gap_scale", "gap_relaxation", "gain_fraction"):
                    values = [float(r[metric]) for r in selected if r[metric]]
                    entry[metric + "_count"] = len(values)
                    entry[metric + "_median"] = (
                        float(np.median(values)) if values else None
                    )
                candidate_summaries.append(entry)
    diagnostic_jobs = {}
    for order in sorted({row["M"] for row in paths}):
        group = [row for row in paths if row["M"] == order]
        for reason, row in (
            (
                "largest_selection_regret",
                max(group, key=lambda r: float(r["refit_score_regret"])),
            ),
            (
                "lowest_selected_gain_fraction",
                min(group, key=lambda r: float(r["selected_gain_fraction"])),
            ),
        ):
            for method in ("relaxation", "refit"):
                k = int(row["k_" + method])
                diagnostic_jobs[row["id"], k] = {
                    **manifests[row["id"]],
                    "k": k,
                    "reason": reason,
                }
    for row in sorted(negative, key=lambda r: float(r["gap_relaxation"]))[:3]:
        diagnostic_jobs[row["id"], int(row["k"])] = {
            **manifests[row["id"]],
            "k": row["k"],
            "reason": "refit_below_relaxation",
        }
    outputs = (
        ("settings.csv", summaries),
        ("candidate_summary.csv", candidate_summaries),
        ("paths_with_support.csv", paths),
        ("negative_gaps.csv", negative),
        ("refit_jobs.csv", list(diagnostic_jobs.values())),
    )
    audit = {
        "settings": len(summaries),
        "paths": len(paths),
        "same_index": sum((row["same_index"] == "True" for row in paths)),
        "same_support": sum((row["same_support"] == "True" for row in paths)),
        "exact_zero_match_at_relaxation_choice": sum(
            (row["selected_exact_zero_match"] == "True" for row in paths)
        ),
        "negative_gaps": len(negative),
        "negative_gaps_exact_zero_match": sum(
            (int(row["tiny_edges"]) == 0 for row in negative)
        ),
        "diagnostic_jobs": len(diagnostic_jobs),
    }
    return (dict(outputs), audit)


def path_comparison(spec, scored):
    candidates = []
    scale, relaxation, refit = (scored[m] for m in METHODS)
    selected = tuple(
        (min(range(40), key=lambda j: (scored[m][j].ebic, j)) for m in METHODS)
    )
    focus = focus_indices(selected, 40)
    comparisons = []
    p = math.prod(spec.dimensions)
    common_unit = math.log(spec.n) + 2 * math.log(p)
    normalized_scores = {
        m: np.array([r.ebic - (r.df - r.df_edges) * common_unit for r in scored[m]])
        for m in METHODS
    }
    metadata = {
        "id": spec.id,
        "family": spec.family,
        "M": spec.M,
        "dimensions": "x".join(map(str, spec.dimensions)),
        "degree": spec.degree,
        "n": spec.n,
        "rep": spec.rep,
    }
    for j, (s, r, f) in enumerate(zip(scale, relaxation, refit, strict=True)):
        if not (
            s.k == r.k == f.k
            and s.eta == r.eta == f.eta
            and (s.df_edges == r.df_edges == f.df_edges)
        ):
            raise AnalysisError(f"unmatched candidate: {spec.id}, {j + 1}")
        if (s.df - s.df_edges, r.df - r.df_edges, f.df - f.df_edges) != (1, 1, 0):
            raise AnalysisError(f"unexpected common scale count: {spec.id}")
        comparison = compare_likelihoods((s.loglik, r.loglik, f.loglik), spec.n * p)
        comparisons.append(comparison)
        candidates.append(
            {
                **metadata,
                "k": s.k,
                "eta": s.eta,
                "edges": s.df_edges,
                "density": s.df_edges / sum((d * (d - 1) / 2 for d in spec.dimensions)),
                "ll_scale": s.loglik,
                "ll_relaxation": r.loglik,
                "ll_refit": f.loglik,
                "focus": j in focus,
                "selected_relaxation": j == selected[1],
                **asdict(comparison),
            }
        )
    at_relaxation = comparisons[selected[1]]
    record = {
        **metadata,
        **{f"k_{m}": selected[i] + 1 for i, m in enumerate(METHODS)},
        "same_index": selected[1] == selected[2],
        "refit_score_regret": refit[selected[1]].ebic - refit[selected[2]].ebic,
        "refit_score_regret_normalized": (
            refit[selected[1]].ebic - refit[selected[2]].ebic
        )
        / (spec.n * p),
        "mcc_difference": relaxation[selected[1]].mcc - refit[selected[2]].mcc,
        "selected_gap_scale": at_relaxation.gap_scale,
        "selected_gap_relaxation": at_relaxation.gap_relaxation,
        "selected_gain_fraction": at_relaxation.gain_fraction,
        "selected_ratio_status": at_relaxation.ratio_status,
        "selected_refit_below_relaxation": at_relaxation.refit_below_relaxation,
    }
    for name, indices in (("all", tuple(range(40))), ("focus", focus)):
        values = [comparisons[j] for j in indices]
        ratios = [v.gain_fraction for v in values if v.gain_fraction is not None]
        record.update(
            {
                f"{name}_candidates": len(indices),
                f"{name}_gap_scale_median": median((v.gap_scale for v in values)),
                f"{name}_gap_relaxation_median": median(
                    (v.gap_relaxation for v in values)
                ),
                f"{name}_gain_fraction_median": median(ratios) if ratios else None,
                f"{name}_gain_fraction_count": len(ratios),
                f"{name}_refit_below_relaxation": sum(
                    (v.refit_below_relaxation for v in values)
                ),
                f"{name}_rank_disagreement": order_disagreement(
                    normalized_scores["relaxation"][list(indices)],
                    normalized_scores["refit"][list(indices)],
                ),
            }
        )
    return (record, candidates)


def normalised_deficits(paths):
    groups = defaultdict(list)
    for row in paths:
        if int(row["M"]) != 2:
            continue
        n = int(row["n"])
        p = math.prod(int(value) for value in row["dimensions"].split("x"))
        key = (row["family"], row["dimensions"], row["degree"], row["n"])
        groups[key].append(
            float(row["selected_gap_relaxation"])
            * n
            * p
            / (math.log(n) + 2 * math.log(p))
        )
    values = [median(group) for group in groups.values()]
    return {
        "settings": len(groups),
        "minimum_setting_median": min(values) if values else None,
        "maximum_setting_median": max(values) if values else None,
    }


def string_records(records):
    return [
        {key: "" if value is None else str(value) for key, value in row.items()}
        for row in records
    ]


def compare_inputs(metadata, scores):
    paths, candidates = [], []
    metadata = metadata[metadata["experiment"] == "refit_comparison"]
    if metadata.empty:
        raise ValueError("no Refit-comparison inputs")
    info = metadata.set_index("id")
    for name, frame in scores[scores["id"].isin(info.index)].groupby(
        "id", observed=True, sort=False
    ):
        row = info.loc[name]
        spec = PathInfo(
            str(name),
            row.family,
            int(row.M),
            tuple(map(int, row.dimensions.split("x"))),
            int(row.degree),
            int(row.n),
            int(row.rep),
        )
        scored = {}
        for method in METHODS:
            selected = frame[
                (frame["criterion"] == method) & (frame["gamma"] == 0.5)
            ].sort_values("k")
            if (
                selected["k"].tolist() != list(range(1, 41))
                or (selected["status"] != "ok").any()
            ):
                raise ValueError(f"incomplete scored path: {name} {method}")
            io.finite(selected, ("eta", "df", "edges", "loglik", "ebic", "mcc"))
            scored[method] = [
                Candidate(
                    int(r.k),
                    float(r.eta),
                    int(r.df),
                    int(r.edges),
                    float(r.loglik),
                    float(r.ebic),
                    float(r.mcc),
                )
                for r in selected.itertuples(index=False)
            ]
        path, records = path_comparison(spec, scored)
        paths.append(path)
        candidates.extend(records)
    if len(paths) != len(metadata):
        raise ValueError("missing Refit-comparison paths")
    return paths, candidates, metadata.to_dict("records")


def run(folder, output):
    io.check_output(output)
    folder = Path(folder)
    metadata = load_paths(folder)
    scores = io.read_frame(
        folder / "scores.csv",
        required=(
            "id",
            "criterion",
            "gamma",
            "k",
            "ebic",
            "loglik",
            "edges",
            "df",
            "mcc",
            "eta",
            "status",
        ),
    )
    supports = io.read_rows(
        folder / "supports.csv", required=("id", "k", "pattern_hash", "tiny_edges")
    )
    paths, candidates, manifests = compare_inputs(metadata, scores)
    scope = {row["id"] for row in paths}
    supports = [row for row in supports if row["id"] in scope]
    outputs, summary = summarize_comparisons(
        string_records(paths),
        string_records(manifests),
        supports,
        string_records(candidates),
    )
    summary["normalised_deficit"] = normalised_deficits(paths)
    with io.output_directory(output) as target:
        io.write_rows(target / "paths.csv", paths)
        io.write_rows(target / "candidates.csv", candidates)
        for name, rows in outputs.items():
            if rows:
                io.write_rows(target / name, rows)
        io.write_json(target / "summary.json", summary)
    return summary
