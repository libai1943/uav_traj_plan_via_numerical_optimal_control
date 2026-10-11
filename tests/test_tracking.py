"""Nontrivial QP regression: known bound-active constrained optimum."""
import unittest
import numpy as np
from python.tracking import box_qp


class TrackingQPTest(unittest.TestCase):
    def test_active_and_interior_variables(self):
        H = np.diag([2., 4., 1.])
        f = np.array([-8., 2., -.25])
        lo, hi = -np.ones(3), np.ones(3)
        x, residual = box_qp(H, f, lo, hi, np.zeros(3))
        np.testing.assert_allclose(x, [1., -.5, .25], atol=2e-5)
        self.assertLess(residual, 2e-3)


if __name__ == '__main__':
    unittest.main()
