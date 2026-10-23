"""体积检查及参数几何的独立回归；不替代实际内核检查。"""

from math import atan2, degrees
import unittest

import numpy as np

from blade_interference import overlap_from_losses
from pump_impeller import Parameters, point_at
from pump_meridional import Parameters as BlankParameters


class OverlapTests(unittest.TestCase):
    def test_disjoint_bodies_have_zero_loss_even_if_receivers_stay_nonempty(self):
        self.assertEqual(overlap_from_losses(12, 9, 12, 9, 1e-6), 0)

    def test_true_overlap_is_measured_in_both_directions(self):
        self.assertEqual(overlap_from_losses(12, 9, 10, 7, 1e-6), 2)

    def test_full_containment(self):
        self.assertEqual(overlap_from_losses(12, 9, 3, 0, 1e-6), 9)

    def test_inconsistent_or_invalid_results_fail_closed(self):
        for values in ((12, 9, 10, 9, 1e-6), (12, 9, 13, 9, 1e-6),
                       (12, 9, float("nan"), 9, 1e-6)):
            with self.subTest(values=values), self.assertRaises(ValueError):
                overlap_from_losses(*values)

    def test_nominal_twisted_sections_leave_clearance_at_72_degree_pitch(self):
        p, blank = Parameters(), BlankParameters()
        minimum_clearance = 360.0/p.blade_count
        # 对各高度剖面的左右边界按真实圆柱半径对齐，独立于切向法平面距离。
        # 将曲线延拓到整个高度包络是保守检查（部分点在工具体外）。
        for z in np.linspace(-p.axial_overcut, 18+p.axial_overcut, 41):
            sides = []
            for side in (-1, 1):
                points = [point_at(r, float(z), side, p, blank)
                          for r in np.linspace(p.inlet_radius, p.outlet_radius, 501)]
                radii = np.array([(x*x+y*y)**0.5 for x,y,_ in points])
                angles = np.unwrap([atan2(y,x) for x,y,_ in points])
                self.assertTrue(np.all(np.diff(radii) > 0))
                sides.append((radii, angles))
            lower = max(s[0][0] for s in sides)
            upper = min(s[0][-1] for s in sides)
            rho = np.linspace(lower, upper, 501)
            width = np.abs(np.interp(rho,*sides[0])-np.interp(rho,*sides[1]))
            clearance = 360.0/p.blade_count-degrees(float(np.max(width)))
            minimum_clearance = min(minimum_clearance, clearance)
        self.assertGreater(minimum_clearance, 50.0)


if __name__ == "__main__":
    unittest.main()
