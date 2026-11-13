import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from replication.adni.python.projection import passband_basis


def diagnostic():
    name = "replication.adni.python.band_diagnostics"
    assert importlib.util.find_spec(name) is not None, (
        "band diagnostic is missing"
    )
    return __import__(name, fromlist=["run"])


def inputs(folder, covariate_effect=False):
    rng = np.random.default_rng(17)
    covariates = rng.normal(size=(24, 3))
    block = rng.normal(size=(24, 4, 34))
    block[:, 1, :] += 0.7 * block[:, 0, :]
    block[:, 2, :] -= 0.4 * block[:, 0, :]
    projected = np.concatenate([block, block, block], axis=2)
    if covariate_effect:
        projected += np.einsum(
            "na,arc->nrc", covariates, rng.normal(size=(3, 4, 102))
        )
    basis = passband_basis(187, 3.0, 0.01, 0.1)
    full = np.einsum("nrc,tc->nrt", projected, basis)
    outside = np.cos(2 * np.pi * 70 * np.arange(187) / 187)
    outside *= np.sqrt(0.25 * np.sum(full**2) / (24 * 4 * np.sum(outside**2)))
    full += outside
    (folder / "combat").mkdir(parents=True)
    (folder / "combat_pb102").mkdir()
    np.save(folder / "combat/Y_CN_combat.npy", full)
    np.save(folder / "combat_pb102/Y_CN_combat.npy", projected)
    pd.DataFrame(covariates, columns=["Age", "Sex", "APOE4"]).to_csv(
        folder / "combat/covars_CN.csv", index=False
    )
    return projected


def test_identical_bands_and_orthogonal_variance(tmp_path):
    source = tmp_path / "input"
    inputs(source)
    output = tmp_path / "output"
    diagnostic().run(source, output, groups=("CN",), splits=5)
    summary = pd.read_csv(output / "summary.csv").iloc[0]
    assert summary["partial_correlation_min"] == pytest.approx(1.0)
    assert summary["partial_correlation_max"] == pytest.approx(1.0)
    assert summary["outside_variance_fraction"] == pytest.approx(0.2)
    bands = pd.read_csv(output / "bands.csv")
    assert len(bands) == 6
    assert set(bands["pair"]) == {"low-mid", "mid-high", "low-high"}
    assert bands["rms_difference"].max() < 1e-12


def test_covariate_adjustment_removes_known_effect(tmp_path):
    for name, effect in [("original", False), ("adjusted", True)]:
        inputs(tmp_path / name, covariate_effect=effect)
        diagnostic().run(
            tmp_path / name,
            tmp_path / (name + "-out"),
            groups=("CN",),
            splits=5,
        )
    a = pd.read_csv(tmp_path / "original-out/bands.csv")
    b = pd.read_csv(tmp_path / "adjusted-out/bands.csv")
    np.testing.assert_allclose(
        a[["correlation", "rms_difference"]],
        b[["correlation", "rms_difference"]],
        atol=1e-12,
    )


def test_random_splits_are_reproducible_and_respond_to_seed(tmp_path):
    inputs(tmp_path / "input")
    for name, seed in [
        ("one", 20261101),
        ("two", 20261101),
        ("three", 20261102),
    ]:
        diagnostic().run(
            tmp_path / "input",
            tmp_path / name,
            groups=("CN",),
            splits=5,
            seed=seed,
        )
    assert (tmp_path / "one/random_splits.csv").read_bytes() == (
        tmp_path / "two/random_splits.csv"
    ).read_bytes()
    assert (tmp_path / "one/random_splits.csv").read_bytes() != (
        tmp_path / "three/random_splits.csv"
    ).read_bytes()
    assert len(pd.read_csv(tmp_path / "one/random_splits.csv")) == 30
    meta = json.loads((tmp_path / "one/metadata.json").read_text())
    assert meta["seed"] == 20261101 and meta["splits"] == 5


def test_existing_output_is_preserved_before_input_loading(tmp_path):
    output = tmp_path / "output"
    output.mkdir()
    sentinel = output / "keep.txt"
    sentinel.write_text("keep")
    with pytest.raises(FileExistsError):
        diagnostic().run(tmp_path / "missing", output)
    assert sentinel.read_text() == "keep"


@pytest.mark.parametrize(
    "problem",
    [
        "missing_covariate",
        "row_count",
        "nonfinite",
        "component_count",
        "zero_variance",
    ],
)
def test_invalid_inputs_leave_no_output(tmp_path, problem):
    source = tmp_path / "input"
    projected = inputs(source)
    covars = source / "combat/covars_CN.csv"
    frame = pd.read_csv(covars)
    if problem == "missing_covariate":
        frame.drop(columns="APOE4").to_csv(covars, index=False)
    elif problem == "row_count":
        frame.iloc[:-1].to_csv(covars, index=False)
    elif problem == "nonfinite":
        projected[0, 0, 0] = np.nan
        np.save(source / "combat_pb102/Y_CN_combat.npy", projected)
    elif problem == "component_count":
        np.save(source / "combat_pb102/Y_CN_combat.npy", projected[:, :, :-1])
    else:
        np.save(source / "combat/Y_CN_combat.npy", np.zeros((24, 4, 187)))
    with pytest.raises(ValueError):
        diagnostic().run(source, tmp_path / "output", groups=("CN",), splits=5)
    assert not (tmp_path / "output").exists()


def test_command_runs_from_repository_root(tmp_path):
    diagnostic()
    inputs(tmp_path / "input")
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "replication.adni.scripts.check_bands",
            "--aggregated",
            str(tmp_path / "input"),
            "--output",
            str(tmp_path / "output"),
            "--groups",
            "CN",
            "--splits",
            "5",
        ],
        cwd=Path(__file__).resolve().parents[4],
        capture_output=True,
        text=True,
        env={**os.environ, "OPENBLAS_NUM_THREADS": "1"},
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "output/summary.csv").is_file()
