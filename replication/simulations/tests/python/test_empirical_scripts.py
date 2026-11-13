import math
import unittest

import numpy as np
import pandas as pd

from replication.simulations.python import (
    conditions,
    convergence,
    covariances,
    intervals,
    likelihood,
    profile_bounds,
    refit,
    scale,
    units,
)


class EmpiricalScriptTests(unittest.TestCase):
    def metadata(self):
        return pd.DataFrame(
            [
                {
                    "id": "a",
                    "experiment": "matrix",
                    "M": 2,
                    "dimensions": "40x40",
                    "family": "ER",
                    "degree": 1,
                    "n": 100,
                    "rep": 1,
                }
            ]
        )

    def test_interval_ignores_the_first_grid_point(self):
        rows = [{"id": "a", "k": 1, "q": 0}, {"id": "a", "k": 2, "q": 2}]
        self.assertEqual(intervals.summarize(self.metadata(), rows)["violations"], 0)
        rows[1]["q"] = 1
        self.assertEqual(intervals.summarize(self.metadata(), rows)["violations"], 1)

    def test_unit_multiplier_recovers_the_trace_scale(self):
        t = np.array([0.02, 0.1, 0.7, 1.0])
        np.testing.assert_allclose(units.invert_t(t - 1 - np.log(t)), t, rtol=1e-13)
        self.assertEqual(units.first_min(np.array([2.0, 1.0, 1.0])), 1)
        np.testing.assert_allclose(
            units.plugin_ebic_at(np.array([1.0, 2.0]), np.array([0.5, 0.25]), 2, 3, 2),
            [10.0, 6.5],
        )

    def test_scale_point_uses_the_selected_candidate(self):
        scores = pd.DataFrame(
            [
                {
                    "id": "a",
                    "criterion": "scale",
                    "status": "ok",
                    "gamma": 0.5,
                    "k": 2,
                    "ebic": 4,
                }
            ]
        )
        values = pd.DataFrame(
            [{"id": "a", "k": 2, "tau": 19, "truth_scale": 20, "B": 8}]
        )
        result = scale.summarize(self.metadata(), scores, values)
        self.assertEqual(result["tau_medians"]["100"], 0.95)
        expected = 8 - (
            (4 - math.log(100) - 2 * math.log(1600)) / 160000 + math.log(1600) - 1
        )
        self.assertAlmostEqual(result["minimum"]["margin"], expected)

    def test_partial_covariance_uses_the_reported_scope(self):
        metadata = pd.concat(
            [
                self.metadata(),
                self.metadata().assign(id="outside", dimensions="116x116"),
            ]
        )
        values = pd.DataFrame(
            [
                {"id": "a", "k": 1, "minimum_eigenvalue": 0.02},
                {"id": "a", "k": 2, "minimum_eigenvalue": 0.01},
                {"id": "outside", "k": 1, "minimum_eigenvalue": -5},
            ]
        )
        result = covariances.summarize(metadata, values)
        self.assertEqual(result["paths"], 1)
        self.assertEqual(result["grid_points"], 2)
        self.assertEqual(result["minimum_eigenvalue"], 0.01)
        self.assertEqual(result["nonpositive_values"], 0)

    def test_condition_numbers_scan_all_recorded_modes(self):
        values = pd.DataFrame(
            [
                {"id": "a", "k": 1, "mode": 1, "cond_base": 2, "cond_relaxed": 4},
                {"id": "a", "k": 1, "mode": 2, "cond_base": 3, "cond_relaxed": 5},
                {"id": "a", "k": 2, "mode": 1, "cond_base": 8, "cond_relaxed": 6},
                {"id": "a", "k": 2, "mode": 2, "cond_base": 7, "cond_relaxed": 4},
            ]
        )
        result = conditions.summarize(self.metadata(), values)
        self.assertEqual(result["grid_points"], 2)
        self.assertEqual(result["maxima"]["cond_base1"]["value"], 8)
        self.assertEqual(result["maxima"]["cond_joint2"]["value"], 5)

    def test_profile_objective_and_gradient_match_dense_calculation(self):
        profile = profile_bounds.Profile([[-0.2, 0.4]], [2.0, 0.3], [1.0], [2])
        point = np.array([0.3])
        expected = np.log1p(point[0] * np.array([-0.2, 0.4])).mean()
        expected -= np.log((2 + 0.3 * point[0]) / 2)
        self.assertAlmostEqual(profile.value(point), expected, places=13)
        h = 1e-6
        numerical = (profile.value(point + h) - profile.value(point - h)) / (2 * h)
        self.assertAlmostEqual(profile.gradient(point)[0], numerical, places=8)

    def test_refit_comparison_preserves_likelihood_units(self):
        result = likelihood.compare_likelihoods((10.0, 18.0, 20.0), 4)
        self.assertEqual(result.gain_fraction, 0.8)
        self.assertEqual(result.gap_scale, 5)
        self.assertEqual(result.gap_relaxation, 1)
        self.assertFalse(result.refit_below_relaxation)
        self.assertTrue(
            likelihood.compare_likelihoods((10.0, 21.0, 20.0), 4).refit_below_relaxation
        )

    def test_convergence_uses_strict_minus_saved_refit(self):
        result = convergence.summarize(
            pd.DataFrame(
                {
                    "saved_refit_ll": [100.0, 100.0],
                    "strict_refit_ll": [100.00000001, 99.99999998],
                    "fresh_refit_ll": [200, 300],
                }
            )
        )
        self.assertEqual(result["candidates"], 2)
        self.assertAlmostEqual(result["maximum_loglikelihood_change"], 2e-8, places=13)

    def test_unit_summary_filters_gamma(self):
        result = units.summarize_multipliers(
            [
                {"gamma": 0.5, "c_star": 2.0},
                {"gamma": 0.5, "c_star": 4.0},
                {"gamma": 0.0, "c_star": 1000.0},
            ]
        )
        self.assertEqual(
            result, {"datasets": 2, "median": 3.0, "minimum": 2.0, "maximum": 4.0}
        )

    def test_deficit_uses_per_setting_medians_and_two_way_paths(self):
        factor = (math.log(100) + 2 * math.log(6)) / 600
        rows = [
            {
                "M": "2",
                "n": "100",
                "dimensions": "2x3",
                "family": family,
                "degree": "1",
                "selected_gap_relaxation": str(value * factor),
            }
            for family, value in (("ER", 1), ("ER", 3), ("BA", 10))
        ]
        rows.append({"M": "3"})
        result = refit.normalised_deficits(rows)
        self.assertEqual(result["settings"], 2)
        self.assertAlmostEqual(result["minimum_setting_median"], 2)
        self.assertAlmostEqual(result["maximum_setting_median"], 10)


if __name__ == "__main__":
    unittest.main()
