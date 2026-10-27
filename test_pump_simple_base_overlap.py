"""纯参数/原生模板回归；不启动 SpaceClaim，也不代表 V241 CAD 验收。"""
import math
import unittest

import pump_simple as pump


class BaseOverlapTests(unittest.TestCase):
    def test_extension_preserves_nominal_edges_and_stays_in_feed_bridge(self):
        for amount in (0.5, 1.0, 2.0):
            with self.subTest(overrun=amount):
                p = dict(pump.PARAMS, volute_loft_radial_overrun=amount)
                pump.validate(p)
                nominal = pump.volute_profiles(p)
                construction = pump.volute_profiles(p, construction=True)
                for i, (target, raw) in enumerate(zip(nominal, construction)):
                    if i >= p['volute_section_count']:
                        self.assertEqual(target['points'], raw['points'])
                        continue
                    delta = raw['radial_overrun']
                    self.assertGreaterEqual(raw['construction_inner_height'], p['volute_min_positive_thickness'])
                    self.assertEqual(target['points'][1], raw['points'][1])
                    self.assertEqual(target['points'][2], raw['points'][2])
                    for inner, outer in ((0, 1), (3, 2)):
                        point = raw['points'][inner]
                        radius = math.hypot(point[0], point[1])
                        self.assertGreater(radius, p['impeller_diameter']/2)
                        self.assertLessEqual(radius, p['volute_base_diameter']/2+1e-10)
                        self.assertGreaterEqual(point[2], p['front_shroud_thickness']-1e-10)
                        self.assertLessEqual(point[2], p['front_shroud_thickness']+p['blade_height']+1e-10)
                        fraction = delta/(delta+target['width'])
                        restored = [(1-fraction)*x+fraction*y for x, y in zip(point, raw['points'][outer])]
                        for actual, expected in zip(restored, target['points'][inner]):
                            self.assertAlmostEqual(actual, expected, places=10)

    def test_invalid_or_zero_thickness_construction_is_rejected(self):
        for changes in ({'volute_loft_radial_overrun': 0},
                        {'volute_loft_radial_overrun': 6},
                        {'volute_boolean_axial_overrun': 0.1},
                        {'volute_min_positive_thickness': 11}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                pump.validate(dict(pump.PARAMS, **changes))

    def test_native_uses_positive_fusion_not_tangent_normalization(self):
        source = pump.spiral_native_source(pump.PARAMS)
        compile(source, '<V241 volute template>', 'exec')
        self.assertNotIn('VoluteInnerTrim', source)
        self.assertNotIn('volute.Subtract', source)
        self.assertIn('verify_base_overlap_and_fuse', source)
        self.assertIn('core_volume + loft_volume - combined.Volume', source)
        self.assertIn('base_overlap_mm3', source)


if __name__ == '__main__':
    unittest.main()
