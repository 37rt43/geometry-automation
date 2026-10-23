"""新 Bezier 毛坯的几何回归：独立切片体积、边界和精确 NURBS 表达。"""

import ast
from dataclasses import replace
from math import hypot, pi
import unittest

import numpy as np
from scipy.integrate import quad
from geomdl import BSpline

from pump_meridional import (Parameters, bezier, boundary_z, check_profile,
                             derivative, expected_volume, profile_data, spaceclaim_script)
from spaceclaim_diagnostics import instrument_script


class MeridionalTests(unittest.TestCase):
    def setUp(self):
        self.p = Parameters()

    def test_specified_controls_and_endpoints(self):
        self.assertEqual(self.p.hub_controls, ((12,60),(12,15),(70,0),(140,0)))
        self.assertEqual(self.p.front_controls, ((45,60),(45,25),(90,8.5),(140,8.5)))
        self.assertEqual(self.p.vertices(), ((12,75),(12,60),(140,0),(140,8.5),(45,60),(45,75)))
        self.assertEqual(bezier(self.p.hub_controls,.5), (49.75,13.125))
        self.assertEqual(bezier(self.p.front_controls,.5), (73.75,21.125))
        self.assertNotEqual(bezier(self.p.hub_controls,1/3),self.p.hub_controls[1])

    def test_tangency_and_strictly_positive_narrowing_width(self):
        report = check_profile(self.p)
        self.assertAlmostEqual(report["minimum_width_mm"],8.5)
        for controls in (self.p.hub_controls,self.p.front_controls):
            self.assertEqual(derivative(controls,0)[0],0)
            self.assertLess(derivative(controls,0)[1],0)
            self.assertGreater(derivative(controls,1)[0],0)
            self.assertEqual(derivative(controls,1)[1],0)

    def test_no_axis_fill_and_no_polygon_self_intersection(self):
        outline = np.asarray(profile_data(self.p)["outline"])
        self.assertTrue(np.array_equal(outline[0],outline[-1]))
        self.assertGreaterEqual(float(outline[:,0].min()),12)
        starts, ends = outline[:-1], outline[1:]
        def cross(a,b):
            return a[...,0]*b[...,1]-a[...,1]*b[...,0]
        for i in range(len(starts)):
            indices = np.arange(i+2,len(starts))
            if i == 0:
                indices = indices[indices != len(starts)-1]
            c,d = starts[indices],ends[indices]
            a,b = starts[i],ends[i]
            crossing = ((cross(b-a,c-a)*cross(b-a,d-a) < -1e-12)
                        & (cross(d-c,a-c)*cross(d-c,b-c) < -1e-12))
            self.assertFalse(np.any(crossing))

    def test_exact_clamped_cubic_nurbs_matches_bernstein(self):
        for controls in (self.p.hub_controls,self.p.front_controls):
            curve = BSpline.Curve()
            curve.degree = 3
            curve.ctrlpts = [[r,0,z] for r,z in controls]
            curve.knotvector = [0,0,0,0,1,1,1,1]
            for u in np.linspace(0,1,501):
                actual = curve.evaluate_single(float(u))
                target = bezier(controls,float(u))
                self.assertLess(hypot(actual[0]-target[0],actual[2]-target[1]),1e-10)

    def test_exact_volume_against_independent_horizontal_annular_slices(self):
        def radius_at_z(controls,z):
            lo,hi = 0.,1.
            for _ in range(60):
                u = (lo+hi)/2
                if bezier(controls,u)[1] > z:
                    lo=u
                else:
                    hi=u
            return bezier(controls,(lo+hi)/2)[0]
        def area(z):
            if z >= 60:
                return pi*(45**2-12**2)
            inner = radius_at_z(self.p.hub_controls,z)
            outer = 140 if z < 8.5 else radius_at_z(self.p.front_controls,z)
            return pi*(outer**2-inner**2)
        volume = sum(quad(area,a,b,epsabs=1e-4,epsrel=1e-10)[0]
                     for a,b in ((0,8.5),(8.5,60),(60,75)))
        self.assertAlmostEqual(volume,expected_volume(self.p),places=4)

    def test_bad_width_or_tangent_rejected(self):
        for p in (replace(self.p,front_outlet_z=-1),replace(self.p,hub_control1_z=60),
                  replace(self.p,front_control2_r=50),replace(self.p,hub_inlet_radius=0)):
            with self.subTest(parameters=p),self.assertRaises(ValueError):
                p.validate()

    def test_native_script_is_new_bezier_geometry_and_valid_syntax(self):
        source = spaceclaim_script(self.p)
        ast.parse(instrument_script(source))
        source.encode("ascii")
        self.assertIn("NurbsData(4,False,False,knots)",source)
        self.assertIn("ControlPoint(point(rz),1.0)",source)
        self.assertNotIn("CreateArc",source)
        self.assertNotIn("Subtract",source)


if __name__ == "__main__":
    unittest.main()
