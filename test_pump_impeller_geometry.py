"""入口相切布尔退化的离线几何回归；实际布尔仍须在 SpaceClaim 验证。"""

from dataclasses import replace
import ast
from math import cos, hypot, sin
import unittest

from pump_impeller import Parameters, blade_sections, loft_sections, native_script
from test_legacy_blank_fixture import LegacyBlankParameters as BlankParameters


class InletLeadinTests(unittest.TestCase):
    def test_new_bezier_blank_requires_new_blade_design(self):
        from pump_meridional import Parameters as NewBlankParameters
        with self.assertRaisesRegex(ValueError, "R140"):
            Parameters().validate(NewBlankParameters())

    def test_current_stage_is_single_blade_only(self):
        Parameters().validate(BlankParameters())
        with self.assertRaises(ValueError):
            replace(Parameters(), blade_count=5).validate(BlankParameters())

    def test_native_single_blade_path_has_no_array_or_pair_check(self):
        tree = ast.parse(native_script(Parameters(), BlankParameters()))
        calls = [node.func.attr for node in ast.walk(tree)
                 if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)]
        self.assertNotIn("Transform", calls)
        self.assertNotIn("Intersect", calls)
        names = {node.id for node in ast.walk(tree) if isinstance(node, ast.Name)}
        self.assertNotIn("overlap_from_losses", names)
        self.assertNotIn("pair_diagnostics", names)
        # 唯一放入主减法工具集合的是 first（眼保护裁切使用另一集合）。
        tools_add = [node for node in ast.walk(tree) if isinstance(node, ast.Call)
                     and isinstance(node.func, ast.Attribute) and node.func.attr == "Add"
                     and isinstance(node.func.value, ast.Name) and node.func.value.id == "tools"]
        self.assertEqual(len(tools_add), 1)
        self.assertEqual(tools_add[0].args[0].id, "first")

    def test_subtraction_uses_detached_blank_before_attaching_result(self):
        tree = ast.parse(native_script(Parameters(), BlankParameters()))
        subtraction = [node for node in ast.walk(tree) if isinstance(node, ast.Call)
                       and isinstance(node.func, ast.Attribute) and node.func.attr == "Subtract"]
        # Both operations target detached shapes, not any DesignBody.Shape.
        self.assertEqual({node.func.value.id for node in subtraction}, {"first", "fluid"})
        copies = [node for node in ast.walk(tree) if isinstance(node, ast.Assign)
                  and any(isinstance(t, ast.Name) and t.id == "fluid" for t in node.targets)]
        self.assertEqual(len(copies), 1)
        self.assertEqual(ast.unparse(copies[0].value), "blank.Shape.Copy()")
        attach = [node for node in ast.walk(tree) if isinstance(node, ast.Call)
                  and isinstance(node.func, ast.Attribute)
                  and isinstance(node.func.value, ast.Name) and node.func.value.id == "DesignBody"
                  and node.func.attr == "Create"]
        self.assertEqual(len(attach), 1)
        self.assertEqual(attach[0].args[2].id, "fluid")
        fluid_subtract = next(node for node in subtraction if node.func.value.id == "fluid")
        self.assertLess(fluid_subtract.lineno, attach[0].lineno)

    def test_preserves_all_31_design_sections(self):
        p, blank = Parameters(), BlankParameters()
        design, angles = blade_sections(p, blank)
        loft, loft_angles = loft_sections(p, blank)
        self.assertEqual(len(design), 31)
        self.assertEqual(len(loft), 32)
        self.assertEqual(loft[1:], design)
        self.assertEqual(loft_angles[1:], angles)

    def test_entire_sacrificial_cap_is_strictly_inside_eye(self):
        # 圆盘为凸集：端面的四个角均严格在圆内，则整个端面均在圆内。
        # 最小延伸量故意减小，验证自动加深仍会避免相切端面。
        for p in (Parameters(), replace(Parameters(), inlet_leadin=0.001),
                  replace(Parameters(), tangential_thickness=5, inlet_twist_degrees=12)):
            with self.subTest(parameters=p):
                sections, angles = loft_sections(p, BlankParameters())
                self.assertTrue(all(hypot(x,y) < p.inlet_radius-p.inlet_cap_clearance
                                    for x,y,z in sections[0]))
                er = cos(angles[0]), sin(angles[0])
                projection = [x*er[0]+y*er[1] for x,y,z in sections[0]]
                self.assertGreater(min(projection), 0)
                self.assertLess(max(projection), p.inlet_radius)

    def test_impossible_cap_is_rejected(self):
        with self.assertRaises(ValueError):
            loft_sections(replace(Parameters(), tangential_thickness=60), BlankParameters())


if __name__ == "__main__":
    unittest.main()
