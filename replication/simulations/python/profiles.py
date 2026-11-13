import json
import math
import time
from .profile_bounds import Profile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from replication import io


def check_candidate(task):
    candidate, fraction, budget = task
    start = time.monotonic()
    profile = Profile(
        candidate["mu"],
        candidate["coefficients"],
        candidate["upper"],
        candidate["p_vec"],
    )
    if len(profile.active) != candidate["q"]:
        raise ValueError("active-coordinate count changed")
    p = math.prod(candidate["p_vec"])
    multiplier = candidate["n"] * p
    penalty = math.log(candidate["n"]) + 2 * math.log(p)
    result = profile.maximise(
        fraction * penalty / multiplier, max_boxes=budget, seed=candidate["u_replay"]
    )
    old_value = 2 * (candidate["loglik_saved"] - candidate["loglik_zero"]) / multiplier
    lower_gain = multiplier * (result["lower"] - old_value)
    upper_gain = multiplier * (result["upper"] - old_value)
    if lower_gain < -0.0001:
        raise ArithmeticError("refined likelihood is below the saved feasible value")
    lower_gain = max(0.0, lower_gain)
    upper_gain = max(lower_gain, upper_gain)
    return dict(
        id=candidate["id"],
        n=candidate["n"],
        k=candidate["k"],
        q=candidate["q"],
        fp=candidate["fp"],
        fn=candidate["fn"],
        support_hash=candidate["support_hash"],
        ebic_saved=candidate["ebic_saved"],
        ebic_feasible=candidate["ebic_saved"] - lower_gain,
        ebic_lower_bound=candidate["ebic_saved"] - upper_gain,
        improvement_over_edge_penalty=lower_gain / penalty,
        old_loss_upper_over_edge_penalty=upper_gain / penalty,
        remaining_gap_over_edge_penalty=(upper_gain - lower_gain) / penalty,
        resolved=result["resolved"],
        boxes=result["boxes"],
        starts=result["starts"],
        converged_starts=result["successes"],
        u=json.dumps(result["point"], separators=(",", ":")),
        spectral_probe_max_error=candidate["spectral_probe_max_error"],
        replay_difference=candidate["replay_difference"],
        seconds=time.monotonic() - start,
    )


def refine_path(task):
    case, candidates, fraction, budget = task
    expected = int(case.get("grid_points", 40))
    if (
        len(candidates) != expected
        or sorted(row["k"] for row in candidates) != list(range(1, expected + 1))
        or any(
            row["id"] != case["id"] or row["n"] != int(case["n"]) for row in candidates
        )
    ):
        raise ValueError("inconsistent candidate path")
    rows = [check_candidate((candidate, fraction, budget)) for candidate in candidates]
    summary = dict(
        summarise_path(rows, case),
        family=case["family"],
        dimensions=case["dimensions"],
        degree=int(case["degree"]),
    )
    return rows, summary


def run(folder, output, workers=1, fraction=0.01, max_boxes=20000):
    io.check_output(output)
    if not 1 <= workers <= 8 or fraction <= 0 or max_boxes < 1:
        raise ValueError("invalid worker count or optimisation budget")
    folder = Path(folder)
    cases = io.read_json(folder / "cases.json")
    if not cases or len({case["id"] for case in cases}) != len(cases):
        raise ValueError("empty or duplicated profile cases")
    tasks = [
        (
            case,
            io.read_json(folder / "profiles" / (case["id"] + ".json"))["candidates"],
            fraction,
            max_boxes,
        )
        for case in cases
    ]
    with io.output_directory(output) as target:
        if workers == 1:
            results = [refine_path(task) for task in tasks]
        else:
            with ProcessPoolExecutor(max_workers=workers) as pool:
                results = list(pool.map(refine_path, tasks))
        rows = [row for candidates, summary in results for row in candidates]
        paths = [summary for candidates, summary in results]
        chosen = {(case["id"], int(case["selected_k"])) for case in cases}
        selected = [row for row in rows if (row["id"], row["k"]) in chosen]
        true_support = [row for row in rows if row["fp"] == row["fn"] == 0]
        summary = {
            "paths": len(paths),
            "candidates": len(rows),
            "index_changes": sum(row["old_k"] != row["refined_k"] for row in paths),
            "support_changes": sum(row["support_changed"] for row in paths),
            "maximum_improvement": max(
                row["improvement_over_edge_penalty"] for row in rows
            ),
            "unresolved_candidates": sum(not row["resolved"] for row in rows),
        }
        for name, group in (("selected", selected), ("true_support", true_support)):
            summary[name + "_maximum_improvement"] = max(
                (row["improvement_over_edge_penalty"] for row in group), default=None
            )
            summary[name + "_maximum_loss_upper"] = max(
                (row["old_loss_upper_over_edge_penalty"] for row in group), default=None
            )
        io.write_rows(target / "candidates.csv", rows)
        io.write_rows(target / "paths.csv", paths)
        io.write_json(target / "summary.json", summary)
    return summary


def summarise_path(rows, original_case):
    ordered_old = sorted(rows, key=lambda x: (x["ebic_saved"], x["k"]))
    ordered_new = sorted(rows, key=lambda x: (x["ebic_feasible"], x["k"]))
    old, new = (ordered_old[0], ordered_new[0])
    if old["k"] != int(original_case["selected_k"]):
        raise ValueError("saved path selection no longer matches the input audit")
    others = [row for row in rows if row["k"] != new["k"]]
    other_supports = [row for row in rows if row["support_hash"] != new["support_hash"]]
    true = [row for row in rows if row["fp"] == row["fn"] == 0]
    if len(true) != int(original_case["exact_candidate_count"]):
        raise ValueError("true-candidate count changed")
    index_margin = (
        min((row["ebic_lower_bound"] for row in others)) - new["ebic_feasible"]
    )
    support_margin = (
        min((row["ebic_lower_bound"] for row in other_supports), default=math.inf)
        - new["ebic_feasible"]
    )
    best_true = min(true, key=lambda x: x["ebic_feasible"]) if true else None
    truth_margin = (
        min((row["ebic_lower_bound"] for row in true)) - new["ebic_feasible"]
        if true
        else None
    )
    return dict(
        id=old["id"],
        n=old["n"],
        rep=int(original_case["rep"]),
        role=original_case["outcome"],
        true_candidates=len(true),
        old_k=old["k"],
        refined_k=new["k"],
        old_exact=old["fp"] == old["fn"] == 0,
        refined_exact=new["fp"] == new["fn"] == 0,
        support_changed=old["support_hash"] != new["support_hash"],
        index_margin_lower=index_margin,
        support_margin_lower=support_margin,
        nominal_index_separated=index_margin > 1e-05,
        nominal_support_separated=support_margin > 1e-05,
        best_true_k=best_true["k"] if best_true else None,
        best_true_old_loss_upper=best_true["old_loss_upper_over_edge_penalty"]
        if best_true
        else None,
        truth_cannot_win_margin=truth_margin,
        selected_improvement=old["improvement_over_edge_penalty"],
        selected_old_loss_upper=old["old_loss_upper_over_edge_penalty"],
        max_candidate_improvement=max(
            (row["improvement_over_edge_penalty"] for row in rows)
        ),
        unresolved=sum((not row["resolved"] for row in rows)),
        max_remaining_gap=max((row["remaining_gap_over_edge_penalty"] for row in rows)),
        old_runner_up_gap=ordered_old[1]["ebic_saved"] - old["ebic_saved"],
    )
