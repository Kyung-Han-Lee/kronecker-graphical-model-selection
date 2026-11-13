import numpy as np


def summarize(frame):
    if frame.empty:
        raise ValueError("no convergence comparisons")
    differences = abs(frame["strict_refit_ll"] - frame["saved_refit_ll"])
    if not np.isfinite(differences).all():
        raise ValueError("nonfinite likelihood difference")
    return {
        "candidates": len(frame),
        "maximum_loglikelihood_change": float(differences.max()),
    }
