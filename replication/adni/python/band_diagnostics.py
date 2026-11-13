from pathlib import Path

import numpy as np

from replication.io import (
    check_output,
    output_directory,
    read_frame,
    write_json,
    write_rows,
)

GROUPS = ("CN", "MCI", "Dementia")
BANDS = {"low": range(0, 34), "mid": range(34, 68), "high": range(68, 102)}
PAIRS = (("low", "mid"), ("mid", "high"), ("low", "high"))


def residuals(values, covariates):
    design = np.column_stack([np.ones(len(covariates)), covariates])
    flat = values.reshape(len(values), -1)
    return (
        flat - design @ np.linalg.lstsq(design, flat, rcond=None)[0]
    ).reshape(values.shape)


def correlations(covariance, partial):
    matrix = np.linalg.inv(covariance) if partial else covariance
    diagonal = np.sqrt(np.diag(matrix))
    return (-matrix if partial else matrix) / np.outer(diagonal, diagonal)


def agreement(left, right):
    upper = np.triu_indices(len(left), 1)
    a, b = left[upper], right[upper]
    if np.std(a) == 0 or np.std(b) == 0:
        raise ValueError("off-diagonal correlations have zero variance")
    return float(np.corrcoef(a, b)[0, 1]), float(np.sqrt(np.mean((a - b) ** 2)))


def load_group(aggregated, group):
    full = np.load(
        aggregated / "combat" / f"Y_{group}_combat.npy", mmap_mode="r"
    )
    projected = np.load(aggregated / "combat_pb102" / f"Y_{group}_combat.npy")
    columns = ["Age", "Sex", "APOE4"]
    covariates = read_frame(
        aggregated / "combat" / f"covars_{group}.csv", required=columns
    )[columns].to_numpy(float)
    if full.ndim != 3 or projected.ndim != 3:
        raise ValueError("ROI tensors must have three dimensions")
    if (
        full.shape[:2] != projected.shape[:2]
        or full.shape[2] != 187
        or projected.shape[2] != 102
    ):
        raise ValueError(
            "expected matching ROI tensors with 187 time points and 102 components"
        )
    if len(covariates) != len(full) or len(full) <= 4 or full.shape[1] < 3:
        raise ValueError(
            "covariate rows must match the tensors, with more than four participants and at least three ROIs"
        )
    if any(
        not np.isfinite(value).all() for value in (full, projected, covariates)
    ):
        raise ValueError("inputs must contain finite values")
    total = float(np.sum(np.square(full)))
    if total <= 0:
        raise ValueError("input variance must be positive")
    share = float(np.sum(np.square(projected))) / total
    if share > 1 + 1e-10:
        raise ValueError("projected power exceeds input power")
    return residuals(projected, covariates), share


def compare_group(values, group, rng, splits):
    components = np.stack(
        [values[:, :, c].T @ values[:, :, c] / len(values) for c in range(102)]
    )
    covariances = {
        name: components[list(indices)].mean(axis=0)
        for name, indices in BANDS.items()
    }
    observed, random = [], []
    for partial, metric in [
        (False, "correlation"),
        (True, "partial_correlation"),
    ]:
        matrices = {
            name: correlations(value, partial)
            for name, value in covariances.items()
        }
        for left, right in PAIRS:
            corr, rms = agreement(matrices[left], matrices[right])
            observed.append(
                {
                    "group": group,
                    "metric": metric,
                    "pair": f"{left}-{right}",
                    "correlation": corr,
                    "rms_difference": rms,
                }
            )
    for split in range(splits):
        permutation = rng.permutation(102)
        parts = [permutation[:34], permutation[34:68], permutation[68:]]
        covariances = [components[part].mean(axis=0) for part in parts]
        for partial, metric in [
            (False, "correlation"),
            (True, "partial_correlation"),
        ]:
            matrices = [correlations(value, partial) for value in covariances]
            for left, right in [(0, 1), (1, 2), (0, 2)]:
                corr, rms = agreement(matrices[left], matrices[right])
                random.append(
                    {
                        "group": group,
                        "split": split + 1,
                        "metric": metric,
                        "pair": f"{left + 1}-{right + 1}",
                        "correlation": corr,
                        "rms_difference": rms,
                    }
                )
    if not all(
        np.isfinite(row[key])
        for rows in [observed, random]
        for row in rows
        for key in ["correlation", "rms_difference"]
    ):
        raise ValueError("nonfinite band diagnostic")
    return observed, random


def run(aggregated, output, groups=GROUPS, splits=100, seed=20261101):
    check_output(output)
    groups = tuple(groups)
    if (
        not groups
        or len(set(groups)) != len(groups)
        or not set(groups).issubset(GROUPS)
    ):
        raise ValueError("groups must be distinct ADNI group names")
    if isinstance(splits, bool) or not isinstance(splits, int) or splits < 1:
        raise ValueError("splits must be a positive integer")
    aggregated = Path(aggregated)
    rng = np.random.default_rng(seed)
    summaries, bands, random_splits = [], [], []
    for group in groups:
        values, share = load_group(aggregated, group)
        observed, random = compare_group(values, group, rng, splits)
        summary = {
            "group": group,
            "participants": len(values),
            "rois": values.shape[1],
            "inside_variance_fraction": share,
            "outside_variance_fraction": 1 - share,
        }
        for metric in ["correlation", "partial_correlation"]:
            ordered = [
                row["correlation"]
                for row in observed
                if row["metric"] == metric
            ]
            shuffled = [
                row["correlation"] for row in random if row["metric"] == metric
            ]
            summary.update(
                {
                    metric + "_min": min(ordered),
                    metric + "_max": max(ordered),
                    metric + "_random_mean": float(np.mean(shuffled)),
                }
            )
        summaries.append(summary)
        bands.extend(observed)
        random_splits.extend(random)
    frequencies = np.fft.rfftfreq(187, d=3.0)
    frequencies = frequencies[(frequencies >= 0.01) & (frequencies <= 0.1)]
    metadata = {
        "seed": seed,
        "splits": splits,
        "groups": list(groups),
        "covariates": ["Age", "Sex", "APOE4"],
        "band_hz": {
            name: [
                float(frequencies[indices[0] // 2]),
                float(frequencies[indices[-1] // 2]),
            ]
            for name, indices in BANDS.items()
        },
    }
    with output_directory(output) as folder:
        write_rows(folder / "summary.csv", summaries)
        write_rows(folder / "bands.csv", bands)
        write_rows(folder / "random_splits.csv", random_splits)
        write_json(folder / "metadata.json", metadata)
    return summaries
