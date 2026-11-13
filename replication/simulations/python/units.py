from collections import defaultdict
import math
from pathlib import Path
from statistics import median

import numpy as np

from replication import io
from .inputs import load_paths

GAMMAS = (0.0, 0.25, 0.5, 1.0)
GAMMA_INDEX = {g: i for i, g in enumerate(GAMMAS)}
N_GRID = 40
C_GRID = (0.01, 0.1, 1.0, 10**0.5, 10.0, 10**1.5, 100.0, 1000.0)


def invert_t(phi: np.ndarray) -> np.ndarray:
    phi = np.asarray(phi, dtype=float)
    if np.any(phi < -1e-10):
        raise ValueError("phi must be nonnegative; got %g" % float(phi.min()))
    phi = np.maximum(phi, 0.0)
    out = np.ones_like(phi)
    pos = phi > 0.0
    if pos.any():
        p = phi[pos]
        s = np.sqrt(2.0 * p)
        u = -s * (1.0 + s)
        for _ in range(200):
            step = (np.expm1(u) - u - p) / np.expm1(u)
            u = u - step
            if np.all(np.abs(step) <= 4e-16 * np.maximum(1.0, np.abs(u))):
                break
        out[pos] = np.exp(u)
    return out


def first_min(values: np.ndarray) -> np.ndarray | int:
    idx = np.argmin(values, axis=-1)
    return int(idx) if np.ndim(idx) == 0 else idx


def plugin_ebic_at(ebic1: np.ndarray, t: np.ndarray, n: int, p: int, c) -> np.ndarray:
    c = np.asarray(c, dtype=float)
    shift = float(n) * float(p) * (c**2 - 1.0)
    return ebic1 + np.multiply.outer(shift, t) if c.ndim else ebic1 + shift * t


def agreement_window(ebic1, t, n, p, target):
    e = np.asarray(ebic1, dtype=float)
    t = np.asarray(t, dtype=float)
    e_t, t_t = (e[target], t[target])
    others = np.arange(len(e)) != target
    if np.any(others & (t == t_t) & (e < e_t)):
        return (math.nan, math.nan)
    steeper = others & (t > t_t)
    flatter = others & (t < t_t)
    a_lo = (
        np.max((e_t - e[steeper]) / (t[steeper] - t_t)) if steeper.any() else -math.inf
    )
    a_hi = (
        np.min((e[flatter] - e_t) / (t_t - t[flatter])) if flatter.any() else math.inf
    )
    if a_lo > a_hi:
        return (math.nan, math.nan)
    scale = float(n) * float(p)
    lo = math.sqrt(1.0 + a_lo / scale) if a_lo > -scale else 0.0
    hi = math.sqrt(1.0 + a_hi / scale) if a_hi < math.inf else math.inf
    return (lo, hi)


def _pick(sc, k) -> dict:
    return {
        "k": int(k) + 1,
        "mcc": sc["mcc"][k],
        "empty": int(sc["df_edges"][k] == 0),
        "exact": int(sc["fp"][k] == 0 and sc["fn"][k] == 0),
        "df_edges": int(sc["df_edges"][k]),
    }


def analyse(row: dict, sc) -> list[dict]:
    rows = []
    for g, gamma in enumerate(GAMMAS):
        k_scale = first_min(sc["scale"][g])
        k_relax = first_min(sc["relax"][g])
        t_scale = float(sc["t"][k_scale])
        c_star = 1.0 / math.sqrt(t_scale)
        lo, hi = agreement_window(sc["plug"][g], sc["t"], sc["n"], sc["p"], k_scale)
        rec = {
            key: row[key]
            for key in (
                "id",
                "experiment",
                "family",
                "degree",
                "p1",
                "p2",
                "n",
                "rep",
                "group",
            )
        }
        rec.update(
            {
                "gamma": gamma,
                "t_scale": t_scale,
                "c_star": c_star,
                "window_lo": lo,
                "window_hi": hi,
                "t_min": float(sc["t"].min()),
                "t_max": float(sc["t"].max()),
                "k_argmin_t": first_min(sc["t"]) + 1,
                "t_gap": float(np.max(np.abs(sc["t_inv"] / sc["t_tau"] - 1.0)))
                if sc["t_tau"] is not None
                else math.nan,
            }
        )
        for name, k in (("scale", k_scale), ("relax", k_relax)):
            for key, value in _pick(sc, k).items():
                rec[f"{name}_{key}"] = value
        for c in C_GRID:
            k = first_min(plugin_ebic_at(sc["plug"][g], sc["t"], sc["n"], sc["p"], c))
            tag = f"c{c:.6g}"
            for key, value in _pick(sc, k).items():
                rec[f"plug_{tag}_{key}"] = value
            rec[f"plug_{tag}_same_scale"] = int(k == k_scale)
            rec[f"plug_{tag}_same_relax"] = int(k == k_relax)
        rows.append(rec)
    return rows


def _tag(c) -> str:
    return f"c{c:.6g}"


def summarise_by_c(records: list[dict]) -> list[dict]:
    out = []
    for gamma in GAMMAS:
        sub = [r for r in records if r["gamma"] == gamma]
        total = len(sub)
        for c in C_GRID:
            t = _tag(c)
            ks = [r[f"plug_{t}_k"] for r in sub]
            out.append(
                {
                    "gamma": gamma,
                    "c": c,
                    "datasets": total,
                    "empty": sum((r[f"plug_{t}_empty"] for r in sub)),
                    "smallest_grid_value": sum((k == 1 for k in ks)),
                    "largest_grid_value": sum((k == N_GRID for k in ks)),
                    "interior": sum((1 < k < N_GRID for k in ks)),
                    "exact": sum((r[f"plug_{t}_exact"] for r in sub)),
                    "same_as_scale": sum((r[f"plug_{t}_same_scale"] for r in sub)),
                    "same_as_relaxation": sum((r[f"plug_{t}_same_relax"] for r in sub)),
                    "mean_mcc": float(np.mean([r[f"plug_{t}_mcc"] for r in sub])),
                }
            )
        for name in ("scale", "relax"):
            ks = [r[f"{name}_k"] for r in sub]
            out.append(
                {
                    "gamma": gamma,
                    "c": "any",
                    "method": name,
                    "datasets": total,
                    "empty": sum((r[f"{name}_empty"] for r in sub)),
                    "smallest_grid_value": sum((k == 1 for k in ks)),
                    "largest_grid_value": sum((k == N_GRID for k in ks)),
                    "interior": sum((1 < k < N_GRID for k in ks)),
                    "exact": sum((r[f"{name}_exact"] for r in sub)),
                    "mean_mcc": float(np.mean([r[f"{name}_mcc"] for r in sub])),
                }
            )
    return out


def summarise_by_setting(records: list[dict], gamma: float = 0.5) -> list[dict]:
    groups = defaultdict(list)
    for r in records:
        if r["gamma"] == gamma:
            groups[r["group"], r["p1"], r["p2"], r["n"]].append(r)
    out = []
    for (group, p1, p2, n), sub in sorted(groups.items()):
        cs = np.array([r["c_star"] for r in sub])
        lo = np.array([r["window_lo"] / r["c_star"] for r in sub])
        hi = np.array([r["window_hi"] / r["c_star"] for r in sub])
        rec = {
            "group": group,
            "p1": p1,
            "p2": p2,
            "n": n,
            "datasets": len(sub),
            "c_star_median": float(np.median(cs)),
            "c_star_q1": float(np.quantile(cs, 0.25)),
            "c_star_q3": float(np.quantile(cs, 0.75)),
            "c_star_min": float(cs.min()),
            "c_star_max": float(cs.max()),
            "window_lo_ratio_median": float(np.median(lo)),
            "window_hi_ratio_median": float(np.median(hi)),
            "scale_exact": sum((r["scale_exact"] for r in sub)),
            "relax_exact": sum((r["relax_exact"] for r in sub)),
        }
        for c in C_GRID:
            t = _tag(c)
            rec[f"{t}_empty"] = sum((r[f"plug_{t}_empty"] for r in sub))
            rec[f"{t}_smallest"] = sum((r[f"plug_{t}_k"] == 1 for r in sub))
            rec[f"{t}_same_scale"] = sum((r[f"plug_{t}_same_scale"] for r in sub))
            rec[f"{t}_exact"] = sum((r[f"plug_{t}_exact"] for r in sub))
        out.append(rec)
    return out


def summarize_multipliers(records):
    values = [row["c_star"] for row in records if row["gamma"] == 0.5]
    if not values:
        raise ValueError("no gamma=0.5 multipliers")
    return {
        "datasets": len(values),
        "median": median(values),
        "minimum": min(values),
        "maximum": max(values),
    }


def candidate_arrays(rows, n, p):
    parts = {}
    for criterion in ("ebic", "scale", "relaxation"):
        selected = rows[rows["criterion"] == criterion]
        if selected.empty or set(selected["gamma"]) != set(GAMMAS):
            raise ValueError(f"incomplete criterion: {criterion}")
        if (selected["status"] != "ok").any():
            raise ValueError(f"unscored candidates: {criterion}")
        io.finite(selected, ("ebic", "loglik", "edges", "fp", "fn", "mcc"))
        frame = selected.pivot(index="gamma", columns="k", values="ebic")
        if frame.shape != (4, 40) or list(frame.columns) != list(range(1, 41)):
            raise ValueError("incomplete gamma by candidate grid")
        parts[criterion] = (selected, frame.to_numpy())
    scale, scale_values = parts["scale"]
    scale = scale[scale["gamma"] == 0].sort_values("k")
    plug = parts["ebic"][0]
    plug = plug[plug["gamma"] == 0].sort_values("k")
    t = invert_t(2 * (scale["loglik"].to_numpy() - plug["loglik"].to_numpy()) / (n * p))
    return {
        "n": n,
        "p": p,
        "scale": scale_values,
        "relax": parts["relaxation"][1],
        "plug": parts["ebic"][1],
        "df_edges": scale["edges"].to_numpy(),
        "fp": scale["fp"].to_numpy(),
        "fn": scale["fn"].to_numpy(),
        "mcc": scale["mcc"].to_numpy(),
        "t_inv": t,
        "t_tau": None,
        "t": t,
    }


def run(folder, output):
    io.check_output(output)
    folder = Path(folder)
    paths = load_paths(folder)
    paths = paths[(paths["experiment"] == "matrix") & (paths["M"] == 2)].set_index("id")
    if paths.empty:
        raise ValueError("no two-way matrix inputs")
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
            "fp",
            "fn",
            "mcc",
            "status",
        ),
    )
    scores = scores[scores["id"].isin(paths.index)]
    records = []
    for name, rows in scores.groupby("id", observed=True, sort=False):
        row = paths.loc[name].to_dict()
        dims = list(map(int, row["dimensions"].split("x")))
        row.update(
            id=str(name),
            p1=dims[0],
            p2=dims[1],
            group=f"{row['family']}_{row['degree']}",
        )
        records.extend(
            analyse(row, candidate_arrays(rows, int(row["n"]), math.prod(dims)))
        )
    if len(records) != len(paths) * len(GAMMAS):
        raise ValueError("missing unit-multiplier inputs")
    summary = summarize_multipliers(records)
    with io.output_directory(output) as target:
        io.write_rows(target / "per_dataset.csv", records)
        io.write_rows(target / "summary_by_c.csv", summarise_by_c(records))
        io.write_rows(
            target / "summary_by_setting_gamma0.5.csv", summarise_by_setting(records)
        )
        io.write_json(target / "summary.json", summary)
    return summary
