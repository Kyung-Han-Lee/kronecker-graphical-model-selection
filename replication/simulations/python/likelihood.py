from __future__ import annotations
from dataclasses import dataclass
from typing import Literal
import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True, slots=True)
class Comparison:
    gap_scale: float
    gap_relaxation: float
    gain_fraction: float | None
    ratio_status: Literal["finite", "near_zero", "negative_refit_gain"]
    refit_below_scale: bool
    refit_below_relaxation: bool
    relaxation_below_scale: bool
    subtraction_tolerance: float


def compare_likelihoods(
    loglik: tuple[float, float, float], total_values: int
) -> Comparison:
    scale, relaxation, refit = loglik
    tolerance = 1e-07 + 1e-10 * max(map(abs, loglik))
    gain = refit - scale
    ratio_status: Literal["finite", "near_zero", "negative_refit_gain"] = "finite"
    ratio = None
    if abs(gain) <= tolerance:
        ratio_status = "near_zero"
    elif gain < 0:
        ratio_status = "negative_refit_gain"
    else:
        ratio = (relaxation - scale) / gain
    return Comparison(
        2 * gain / total_values,
        2 * (refit - relaxation) / total_values,
        ratio,
        ratio_status,
        gain < -tolerance,
        refit - relaxation < -tolerance,
        relaxation - scale < -tolerance,
        tolerance,
    )


def focus_indices(selected: tuple[int, ...], size: int) -> tuple[int, ...]:
    return tuple(
        sorted(
            {
                neighbor
                for index in selected
                for neighbor in (index - 1, index, index + 1)
                if 0 <= neighbor < size
            }
        )
    )


def order_disagreement(a: NDArray[np.float64], b: NDArray[np.float64]) -> float | None:
    left, right = np.triu_indices(len(a), k=1)
    da, db = (a[left] - a[right], b[left] - b[right])
    tolerance = 2e-07 + 1e-10 * max(float(np.max(np.abs(a))), float(np.max(np.abs(b))))
    comparable = (np.abs(da) > tolerance) & (np.abs(db) > tolerance)
    if not comparable.any():
        return None
    return float(np.mean(np.sign(da[comparable]) != np.sign(db[comparable])))
