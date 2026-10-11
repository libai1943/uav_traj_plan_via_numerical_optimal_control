"""Geometry and model regression checks. Run: python -m unittest discover -s tests -v"""
import json
from pathlib import Path
import unittest
import numpy as np
from python.geometry import segment_clearance, box_clearance, astar_3d, spline_seed
from python.model import motion


class GeometryTests(unittest.TestCase):
    def test_edge_check_detects_collision_between_free_endpoints(self):
        obs = np.array([[0., 0., 1.]])
        self.assertGreater(segment_clearance([-2, 0, 100], [-2, 0, 100], obs), 0)
        self.assertGreater(segment_clearance([2, 0, 100], [2, 0, 100], obs), 0)
        self.assertLess(segment_clearance([-2, 0, 100], [2, 0, 100], obs), 0)

    def test_whole_box_check_detects_obstacle_interior(self):
        obs = np.array([[0., 0., 1.]])
        # All four XY corners are clear, but the rectangle contains the column axis.
        self.assertAlmostEqual(box_clearance(np.array([-2, -2, 80]), np.array([2, 2, 180]), obs), -1)
        self.assertAlmostEqual(box_clearance(np.array([2, -1, 80]), np.array([3, 1, 180]), obs), 1)

    def test_true_3d_search_and_checked_spline(self):
        cfg = json.loads((Path(__file__).resolve().parents[1]/'scenes.json').read_text())
        case = cfg['cases'][0]; obs = np.array(case['columns'])
        a, b = np.array(cfg['start_state'][:3]), np.array(case['goal_state'][:3])
        path, count = astar_3d(a, b, obs, cfg)
        np.testing.assert_array_equal(path[0], a); np.testing.assert_array_equal(path[-1], b)
        self.assertGreater(count, 0); self.assertGreater(np.ptp(path[:, 2]), 0)
        moves = np.diff(path, axis=0)/cfg['grid_step_m']
        self.assertTrue(np.all(np.abs(moves) <= 1))
        for p, q in zip(path[:-1], path[1:]):
            self.assertGreater(segment_clearance(p, q, obs), cfg['safety_margin_m']+cfg['astar_extra_clearance_m'])
        seed, dense, _, certificate = spline_seed(path, obs, cfg)
        self.assertEqual(seed.shape, (cfg['intervals']+1, 3))
        self.assertGreater(certificate, cfg['safety_margin_m'])
        actual_samples = min(np.min(np.linalg.norm(dense[:, :2]-o[:2], axis=1)-o[2]) for o in obs)
        self.assertLessEqual(certificate, actual_samples+1e-8)


class DynamicsTests(unittest.TestCase):
    def test_straight_level_trim(self):
        x = np.array([0., 0., 100., 20., 0., 0., 0.]); u = np.array([0., 1., 0.])
        np.testing.assert_allclose(motion(x, u), [20, 0, 0, 0, 0, 0, 0], atol=1e-12)

    def test_all_dynamics_jacobian_columns_off_trim(self):
        x = np.array([12., -8., 115., 19., .12, -.18, .25]); u = np.array([.2, 1.15, -.07])
        _, A, B = motion(x, u, True)
        for j in range(7):
            dx = np.zeros(7); dx[j] = 1e-6
            np.testing.assert_allclose((motion(x+dx, u)-motion(x-dx, u))/2e-6, A[:, j], atol=2e-8)
        for j in range(3):
            du = np.zeros(3); du[j] = 1e-6
            np.testing.assert_allclose((motion(x, u+du)-motion(x, u-du))/2e-6, B[:, j], atol=2e-8)


if __name__ == '__main__':
    unittest.main()
