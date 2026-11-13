import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = ROOT / "../outputs/adni/bootstrap"
SUMMARY_PATH = RESULT_DIR / "bootstrap_summary.csv"
RUN_PATH = RESULT_DIR / "RUN_bootstrap.json"
OUTPUT_PATH = ROOT / "../outputs/adni/tables/TableC1.csv"

GROUPS = ("CN", "MCI", "Dementia")
DIFFERENCE = "Dementia-CN"
PAIRS = ("posterior_cingulate", "precuneus")
GROUP_LABELS = {
    "CN": "Cognitively normal",
    "MCI": "Mild cognitive impairment",
    "Dementia": "Dementia",
    DIFFERENCE: "Dementia $-$ cognitively normal",
}
PAIR_LABELS = {
    "posterior_cingulate": "Posterior cingulate",
    "precuneus": "Precuneus",
}


def group_sizes(run_path: Path = RUN_PATH) -> dict:
    run = json.loads(run_path.read_text(encoding="utf-8"))
    return dict(zip(run["groups"], run["n_by_group"]))


def prepare_table_data(
    summary: pd.DataFrame, n_by_group: dict | None = None
) -> pd.DataFrame:
    required = {"group", "pair", "estimate", "lower", "upper", "n_boot"}
    missing = sorted(required - set(summary.columns))
    if missing:
        raise ValueError(f"summary is missing columns: {', '.join(missing)}")
    keys = set(summary[["group", "pair"]].itertuples(index=False, name=None))
    expected = {(group, pair) for group in (*GROUPS, DIFFERENCE) for pair in PAIRS}
    if summary.duplicated(["group", "pair"]).any() or keys != expected:
        raise ValueError(
            "expected one row per group and pair, including the difference"
        )
    if summary["n_boot"].nunique() != 1:
        raise ValueError("every row must use the same number of resamples")
    if n_by_group is None:
        n_by_group = group_sizes()

    indexed = summary.set_index(["group", "pair"])
    rows = []
    for group in (*GROUPS, DIFFERENCE):
        row = {
            "group": group,
            "n": n_by_group[group] if group in GROUPS else None,
        }
        for pair in PAIRS:
            for column in ("estimate", "lower", "upper"):
                row[f"{pair}_{column}"] = float(indexed.loc[(group, pair), column])
        rows.append(row)
    data = pd.DataFrame(rows)

    data["n"] = pd.Series([row["n"] for row in rows], dtype=object)
    return data


def _number(value: float) -> str:

    text = f"{value:.3f}"
    return f"${text}$" if value < 0 else text


def _interval(lower: float, upper: float) -> str:
    return f"$({lower:.3f},{upper:.3f})$"


def render_table(data: pd.DataFrame) -> str:
    lines = [
        r"\begin{table}[t]",
        r"\centering",
        r"\begin{minipage}{\textwidth}",
        r"\caption{Participant-bootstrap intervals for left--right partial correlations "
        r"in the ADNI3 ROI networks selected by Relaxation at $\gamma=0.5$\label{tab:adni-bootstrap}}",
        r"\begingroup",
        r"\renewcommand{\baselinestretch}{1}",
        r"\fontsize{10}{10}\selectfont",
        r"\renewcommand{\arraystretch}{1.6}",
        r"\tabcolsep=2pt",
        r"\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}llcccc@{\extracolsep{\fill}}}",
        r"\toprule",
        rf"& & \multicolumn{{2}}{{c}}{{{PAIR_LABELS[PAIRS[0]]}}} & \multicolumn{{2}}{{c}}{{{PAIR_LABELS[PAIRS[1]]}}} \\",
        r"\cline{3-4}\cline{5-6}",
        r"Group & $n$ & Estimate & $95\%$ interval & Estimate & $95\%$ interval \\",
        r"\midrule",
    ]
    for row in data.to_dict("records"):
        cells = [
            GROUP_LABELS[row["group"]],
            "" if row["n"] is None else str(row["n"]),
        ]
        for pair in PAIRS:
            cells.append(_number(row[f"{pair}_estimate"]))
            cells.append(_interval(row[f"{pair}_lower"], row[f"{pair}_upper"]))
        lines.append(" & ".join(cells) + r" \\")
    lines.extend(
        [
            r"\botrule",
            r"\end{tabular*}",
            r"\begin{tablenotes}",
            r"\item \textit{Note}. Estimate = partial correlation between the left and right regions in the relaxation-profiled ROI factor; "
            r"$95\%$ interval = $2.5$th and $97.5$th bootstrap percentiles.",
            r"The last row subtracts the cognitively normal value from the dementia value at each resample index.",
            r"\end{tablenotes}",
            r"\endgroup",
            r"\end{minipage}",
            r"\end{table}",
            "",
        ]
    )
    return "\n".join(lines)


def run(summary_path, metadata_path, output):
    summary = pd.read_csv(summary_path)
    data = prepare_table_data(summary, group_sizes(metadata_path))
    output = Path(output)
    if output.exists():
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.suffix == ".tex":
        output.write_text(render_table(data))
    else:
        data.to_csv(output, index=False)
    return output
