from pathlib import Path

import pandas as pd
import pytest

from replication.adni.python import bootstrap_table as tab_bootstrap


@pytest.fixture(autouse=True)
def counts(monkeypatch):
    monkeypatch.setattr(
        tab_bootstrap,
        "group_sizes",
        lambda *args: {"CN": 289, "MCI": 158, "Dementia": 73},
    )


def summary_frame() -> pd.DataFrame:
    rows = [
        ("Dementia", "posterior_cingulate", 0.5866034, 0.5589, 0.62402, 1.0),
        ("Dementia", "precuneus", 0.5480045, 0.5140, 0.59324, 1.0),
        ("MCI", "posterior_cingulate", 0.6161255, 0.6087, 0.64545, 1.0),
        ("MCI", "precuneus", 0.6003716, 0.5937, 0.63193, 1.0),
        ("CN", "posterior_cingulate", 0.6630227, 0.6472, 0.67475, 1.0),
        ("CN", "precuneus", 0.6379070, 0.6226, 0.64704, 1.0),
        (
            "Dementia-CN",
            "posterior_cingulate",
            -0.0764193,
            -0.1051,
            -0.03441,
            None,
        ),
        ("Dementia-CN", "precuneus", -0.0899025, -0.1234, -0.03937, None),
    ]
    frame = pd.DataFrame(
        rows,
        columns=[
            "group",
            "pair",
            "estimate",
            "lower",
            "upper",
            "selected_freq",
        ],
    )
    frame["n_boot"] = 500
    return frame


def test_output_is_an_appendix_table_file():
    assert tab_bootstrap.OUTPUT_PATH.name == "TableC1.csv"
    assert tab_bootstrap.SUMMARY_PATH.parts[-2:] == (
        "bootstrap",
        "bootstrap_summary.csv",
    )


def test_rows_follow_group_order_and_end_with_the_difference():
    data = tab_bootstrap.prepare_table_data(summary_frame())
    assert list(data["group"]) == ["CN", "MCI", "Dementia", "Dementia-CN"]
    assert list(data["n"]) == [289, 158, 73, None]


def test_render_prints_three_decimals_and_math_minus():
    tex = tab_bootstrap.render_table(tab_bootstrap.prepare_table_data(summary_frame()))
    assert (
        "Cognitively normal & 289 & 0.663 & $(0.647,0.675)$ & 0.638 & $(0.623,0.647)$ \\\\"
        in tex
    )
    assert (
        "Dementia $-$ cognitively normal &  & $-0.076$ & $(-0.105,-0.034)$"
        " & $-0.090$ & $(-0.123,-0.039)$ \\\\" in tex
    )
    assert r"\label{tab:adni-bootstrap}" in tex
    assert tex.count(r"\begin{table}") == 1 and tex.rstrip().endswith(r"\end{table}")


def test_render_leaves_the_number_of_resamples_to_appendix_c6():
    tex = tab_bootstrap.render_table(tab_bootstrap.prepare_table_data(summary_frame()))
    assert "resamples" not in tex


def test_missing_group_is_an_error():
    frame = summary_frame()
    frame = frame[frame["group"] != "MCI"]
    with pytest.raises(ValueError, match="one row per group and pair"):
        tab_bootstrap.prepare_table_data(frame)


def test_unequal_resample_counts_are_an_error():
    frame = summary_frame()
    frame.loc[0, "n_boot"] = 499
    with pytest.raises(ValueError, match="same number of resamples"):
        tab_bootstrap.prepare_table_data(frame)


def test_written_table_matches_the_supplied_summary(tmp_path: Path):
    source = tmp_path / "summary.csv"
    metadata = tmp_path / "run.json"
    output = tmp_path / "table.csv"
    summary_frame().to_csv(source, index=False)
    metadata.write_text('{"groups":["CN","MCI","Dementia"],"n_by_group":[289,158,73]}')
    tab_bootstrap.run(source, metadata, output)
    data = pd.read_csv(output)
    assert data.iloc[0]["posterior_cingulate_estimate"] == pytest.approx(0.6630227)


def test_title_ends_without_a_full_stop_and_the_note_has_two_sentences():
    tex = tab_bootstrap.render_table(tab_bootstrap.prepare_table_data(summary_frame()))
    assert ".\\label{" not in tex
    notes = tex[
        tex.index(r"\begin{tablenotes}") : tex.index(r"\end{tablenotes}")
    ].splitlines()[1:]
    assert notes == [
        r"\item \textit{Note}. Estimate = partial correlation between the left and right regions in the relaxation-profiled ROI factor; $95\%$ interval = $2.5$th and $97.5$th bootstrap percentiles.",
        r"The last row subtracts the cognitively normal value from the dementia value at each resample index.",
    ]


def test_title_names_matching_left_and_right_regions():
    tex = tab_bootstrap.render_table(tab_bootstrap.prepare_table_data(summary_frame()))
    assert (
        r"\caption{Participant-bootstrap intervals for left--right partial correlations "
        r"in the ADNI3 ROI networks "
        r"selected by Relaxation at $\gamma=0.5$\label{tab:adni-bootstrap}}"
    ) in tex
