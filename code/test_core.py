"""Code-only checks using in-memory synthetic arrays, not study records."""
import contextlib
import io
import unittest

import numpy as np
from scipy.stats import spearmanr

import GAM_methods
import fixed_GAM_methods
import principal_curve_methods
from compare_k_resolutions import historical_bins, interp_nan, select_directions


class NumericalTests(unittest.TestCase):
    def test_native_segment_projection(self):
        queries = np.array([[0.2, 2.0], [1.8, -1.0]])
        line = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
        arc, _, squared_distance, _ = principal_curve_methods.project_polyline(queries, line)
        np.testing.assert_allclose(arc, [0.2, 1.8])
        np.testing.assert_allclose(squared_distance, [4.0, 1.0])

    def test_curve_invariances(self):
        with contextlib.redirect_stdout(io.StringIO()):
            principal_curve_methods.test_algorithm()

    def test_fixed_spline(self):
        age = np.repeat(np.linspace(20, 90, 31), 4)
        rng = np.random.default_rng(17)
        score = ((age - 20) / 70) ** 2 + rng.normal(0, 0.02, len(age))
        curve, report = fixed_GAM_methods.fixed_gam(age, score, GAM_methods, total_edf=6)
        self.assertAlmostEqual(report['effective_degrees_of_freedom'], 6, places=7)
        self.assertTrue(np.isfinite(curve.fitted_pseudotime).all())
        self.assertTrue(report['no_monotonicity_convexity_or_concavity_constraint'])

    def test_historical_bins_and_fallback(self):
        bins, centers = historical_bins()
        self.assertEqual(len(bins), 63)
        self.assertIn((44, 46), bins)
        self.assertIn((45, 46), bins)
        self.assertFalse(any(lower <= 90 < upper for lower, upper in bins))
        np.testing.assert_allclose(interp_nan(np.array([np.nan, 1., np.nan, 3., np.nan]),
                                            np.arange(5)), [1, 1, 2, 3, 3])
        self.assertEqual(len(centers), len(bins))

    def test_directional_selection(self):
        import pandas as pd
        table = pd.DataFrame({'cluster': [0, 1, 2], 'n_analytes': [2, 3, 4],
                              'rho_module_age': [0.4, -0.8, 0.6]})
        up, down, value = select_directions(table)
        self.assertEqual(int(up.cluster), 2)
        self.assertEqual(int(down.cluster), 1)
        self.assertAlmostEqual(value, 0.7)
        missing = table.copy()
        missing['rho_module_age'] = abs(missing.rho_module_age)
        self.assertTrue(np.isnan(select_directions(missing)[2]))


if __name__ == '__main__':
    unittest.main(verbosity=2)
