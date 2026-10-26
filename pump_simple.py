"""阶段一：空心轴向进口管和带直叶片的封闭式叶轮，长度数值单位 mm。

Windows / SpaceClaim 2024R1: python pump_simple.py --output-dir output
离线参数检查: python pump_simple.py --check-only
修改 PARAMS 后整模重建；STEP 不保存参数特征树，以本脚本为参数源。
Z 为轴线，前盖板 [-tf, 0]，叶片 [0, h]，后盖板 [h, h+tb]。
进口管位于前盖板上游，轴向留装配间隙。流道为空腔，未另建流体实体。
"""

PARAMS = {
    "inlet_diameter": 50.0,       # 进口管内径，与叶轮进口眼相等
    "inlet_length": 80.0,         # 进口直管长度
    "inlet_wall_thickness": 4.0,  # 进口管壁厚，避免实心圆柱堵塞流道
    "inlet_axial_gap": 0.5,       # 进口管端面至前盖板上游端面间隙
    "hub_diameter": 30.0,         # 直筒形旋转轮毂外径
    "shaft_diameter": 16.0,       # 轮毂及后盖板贯通轴孔直径
    "front_shroud_diameter": 120.0,  # 前盖板外径
    "rear_shroud_diameter": 120.0,   # 后盖板外径
    "front_shroud_thickness": 3.0,  # 前盖板厚度
    "rear_shroud_thickness": 4.0,   # 后盖板厚度
    "eye_diameter": 50.0,         # 前盖板进口眼直径，保留轴向入口
    "blade_count": 6,             # 直叶片圆周阵列数量（整数）
    "blade_thickness": 3.0,       # 叶片切向厚度，矩形单截面
    "blade_height": 12.0,         # 前后盖板之间的轴向净流道高度
    "impeller_outlet_diameter": 120.0,  # 叶片出口包络直径，不是出口管径
    "volute_base_diameter": 132.0,  # 预留：蜗壳基圆直径，本阶段不建模
    "volute_height": 24.0,         # 预留：蜗壳轴向高度
    "outlet_pipe_diameter": 40.0,  # 预留：出口管内径
    "outlet_pipe_length": 100.0,   # 预留：出口直管长度
    "outlet_wall_thickness": 4.0,  # 预留：出口管壁厚
    "seal_chamber_diameter": 40.0, # 预留：密封腔腔径，轴径用 shaft_diameter
    "seal_chamber_length": 20.0,   # 预留：密封腔轴向长度
    "wear_ring_radial_clearance": 0.25, # 预留：口环单边径向间隙
    "wear_ring_length": 8.0,       # 预留：口环轴向长度
    "rear_chamber_gap": 2.0,       # 预留：后腔轴向间隙
    "volume_relative_tolerance": 1e-5, # CAD/解析体积相对容差（无量纲）
    "volume_absolute_tolerance": 0.01, # CAD/解析体积绝对容差（mm³）
}

import argparse
import json
from math import asin, cos, isfinite, pi, sin, sqrt
import os
from pathlib import Path
import tempfile

from main import connect  # 无 host 时调用 launch_modeler_with_spaceclaim(version=241, ...)


def validate(p):
    for key, value in p.items():
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not isfinite(value) or value <= 0:
            raise ValueError(key + " 必须为有限正数")
    n = p["blade_count"]
    if type(n) is not int or n < 2:
        raise ValueError("blade_count 必须为 >= 2 的整数")
    rs, rh, re, ro = (p[k]/2 for k in ("shaft_diameter", "hub_diameter", "eye_diameter", "impeller_outlet_diameter"))
    a = p["blade_thickness"]/2
    if not rs < rh < re < ro:
        raise ValueError("需满足 轴孔 < 轮毂 < 进口眼 < 叶轮出口直径")
    if p["inlet_diameter"] != p["eye_diameter"]:
        raise ValueError("简化直通入口要求 inlet_diameter == eye_diameter")
    if a >= rh*sin(pi/n) or a*a + rs*rs >= rh*rh:
        raise ValueError("叶片在轮毂外相互干涉，须减小厚度或叶片数")
    if p["inlet_diameter"] + 2*p["inlet_wall_thickness"] >= p["front_shroud_diameter"]:
        raise ValueError("进口管外径须小于前盖板外径")
    if min(p["front_shroud_diameter"], p["rear_shroud_diameter"]) < 2*ro:
        raise ValueError("前后盖板须覆盖叶片出口包络")
    if p["seal_chamber_diameter"] <= p["shaft_diameter"]:
        raise ValueError("预留密封腔径须大于轴径")
    return p


def expected_volumes(p):
    validate(p)
    ri, rh, rs = (p[k]/2 for k in ("inlet_diameter", "hub_diameter", "shaft_diameter"))
    a, ro = p["blade_thickness"]/2, p["impeller_outlet_diameter"]/2
    # 矩形叶片端部为内接弦；扣除与圆轮毂相交面积，避免并集重复计体积。
    tip_x = sqrt(ro*ro-a*a)
    hub_strip = a*sqrt(rh*rh-a*a) + rh*rh*asin(a/rh)
    blade_area = 2*a*tip_x - hub_strip
    return {
        "inlet": pi*((ri+p["inlet_wall_thickness"])**2-ri**2)*p["inlet_length"],
        "impeller": (pi*(rh**2-rs**2) + p["blade_count"]*blade_area)*p["blade_height"]
            + pi*((p["front_shroud_diameter"]/2)**2-(p["eye_diameter"]/2)**2)*p["front_shroud_thickness"]
            + pi*((p["rear_shroud_diameter"]/2)**2-rs**2)*p["rear_shroud_thickness"],
    }


def run_native(modeler, source):
    # V241 原生脚本支持旋转与拓扑检查；PyAnsys revolve_sketch/rotate 要求 24R2。
    with tempfile.TemporaryDirectory(prefix="pump_simple_") as folder:
        path = Path(folder)/"native.py"
        path.write_text(source, encoding="utf-8")
        return modeler.run_script_file(path, import_design=True, api_version="241")


def hub_script(p):
    return '''# Python Script, API Version = V241
from System.Collections.Generic import List
root = GetRootPart()
if root.Bodies.Count != 0:
    raise Exception("Hub requires an empty design")
coordinates = %r
curves = List[ITrimmedCurve]()
def pt(rz):
    return Point.Create(MM(rz[0]), 0, MM(rz[1]))
for i in range(len(coordinates)):
    curves.Add(CurveSegment.Create(pt(coordinates[i]), pt(coordinates[(i+1) %% len(coordinates)])))
plane = Plane.Create(Frame.Create(Point.Origin, Direction.DirX, Direction.DirZ))
profile = PlanarBody.Create(plane, curves).CreatedBody
options = RevolveFaceOptions()
options.ExtrudeType = ExtrudeType.Add
RevolveFaces.Execute(FaceSelection.Create(profile.Faces[0]),
                    Line.Create(Point.Origin, Direction.DirZ), DEG(360), options)
if root.Bodies.Count != 1:
    raise Exception("Expected one revolved hub")
root.Bodies[0].SetName("impeller")
result = {"hub_revolved": "True"}
''' % ([(p["shaft_diameter"]/2, 0), (p["hub_diameter"]/2, 0),
        (p["hub_diameter"]/2, p["blade_height"]), (p["shaft_diameter"]/2, p["blade_height"])],)


CLOSED_CHECK = '''# Python Script, API Version = V241
root = GetRootPart()
if root.Bodies.Count != 1 or root.Components.Count != 0:
    raise Exception("Expected exactly one solid, no leftover tools/components")
body = root.Bodies[0]
if not body.Shape.IsClosed or body.Shape.PieceCount != 1 or body.Shape.Volume <= 0:
    raise Exception("Body must be closed, connected and positive-volume")
result = {"closed": "True", "pieces": "1"}
'''


def build_part(modeler, name, p):
    from ansys.geometry.core.math import Plane, Point2D, Point3D
    from ansys.geometry.core.misc import UNITS
    from ansys.geometry.core.sketch import Sketch
    mm = UNITS.mm
    design = modeler.create_design("pump_simple_"+name)

    def cylinder(label, z, radius, length):
        sk = Sketch(Plane(origin=Point3D([0, 0, z], unit=mm)))
        sk.circle(Point2D([0, 0], unit=mm), radius*mm)
        return design.extrude_sketch(label, sk, distance=length*mm)

    def ring(label, z, outer, inner, length):
        body = cylinder(label, z, outer, length)
        body.subtract(cylinder(label+"_bore", z, inner, length))
        return body

    if name == "inlet":
        z = -p["front_shroud_thickness"]-p["inlet_axial_gap"]-p["inlet_length"]
        ring(name, z, p["inlet_diameter"]/2+p["inlet_wall_thickness"], p["inlet_diameter"]/2, p["inlet_length"])
    else:
        values, design = run_native(modeler, hub_script(p))
        if values.get("hub_revolved") != "True":
            raise RuntimeError("原生轮毂旋转未确认完成")
        hub = design.get_all_bodies()[0]
        h = p["blade_height"]
        # 后盖板先与轮毂合并，随后合并阵列叶片，最后加前盖板。
        hub.unite(ring("rear_shroud", h, p["rear_shroud_diameter"]/2,
                       p["shaft_diameter"]/2, p["rear_shroud_thickness"]))
        half = p["blade_thickness"]/2
        x0 = (p["shaft_diameter"]/2 + sqrt((p["hub_diameter"]/2)**2-half**2))/2
        # 根部完全落在轮毂内且在轴孔外，合并后轴孔保持贯通。
        x1 = sqrt((p["impeller_outlet_diameter"]/2)**2-half**2)
        for i in range(p["blade_count"]):
            theta = 2*pi*i/p["blade_count"]
            # 24.1 通过旋转草图坐标重建同一单截面，实现圆周阵列；不用新版 body.rotate。
            plane = Plane(direction_x=[cos(theta), sin(theta), 0],
                          direction_y=[-sin(theta), cos(theta), 0])
            sk = Sketch(plane)
            sk.box(Point2D([(x0+x1)/2, 0], unit=mm), (x1-x0)*mm, 2*half*mm)
            blade = design.extrude_sketch("blade_%02d" % (i+1), sk, distance=h*mm)
            hub.unite(blade)
        hub.unite(ring("front_shroud", -p["front_shroud_thickness"],
                       p["front_shroud_diameter"]/2, p["eye_diameter"]/2, p["front_shroud_thickness"]))
    values, design = run_native(modeler, CLOSED_CHECK)
    if values.get("closed") != "True" or values.get("pieces") != "1":
        raise RuntimeError(name+" 封闭性/连通性检查未通过")
    bodies = design.get_all_bodies()
    if len(bodies) != 1 or bodies[0].is_surface:
        raise RuntimeError(name+" 必须为一个实体")
    actual = bodies[0].volume.to("mm^3").magnitude
    target = expected_volumes(p)[name]
    allowed = max(p["volume_absolute_tolerance"], target*p["volume_relative_tolerance"])
    if not isfinite(actual) or abs(actual-target) > allowed:
        raise RuntimeError(f"{name} 体积异常：{actual}，预期 {target} mm³，容差 {allowed}")
    return design, {"closed": True, "pieces": 1, "volume_mm3": actual, "expected_mm3": target}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default=os.getenv("ANSRV_GEO_HOST"))
    parser.add_argument("--port", type=int, default=int(os.getenv("ANSRV_GEO_PORT", "50051")))
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent/"output")
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    p = validate(dict(PARAMS))
    targets = expected_volumes(p)
    if args.check_only:
        print(json.dumps({"status": "offline_only_no_CAD", "expected_volumes_mm3": targets,
                          "axial_part_gap_mm": p["inlet_axial_gap"]}, indent=2))
        return
    from ansys.geometry.core.connection.backend import BackendType
    from ansys.geometry.core.designer.design import DesignFileFormat
    modeler = connect(args.host, args.port)
    try:
        if not modeler.client.healthy or modeler.client.backend_type != BackendType.SPACECLAIM:
            raise RuntimeError("需要健康的 SPACECLAIM 后端")
        output = args.output_dir.resolve()
        output.mkdir(parents=True, exist_ok=True)
        # 先清理旧完成标记，失败的重建不会留下本轮成功报告。
        report_path = output/"pump_simple_checks.json"
        report_path.unlink(missing_ok=True)
        report = {"params": p, "parts": {}, "interference_check": {
            "method": "analytic axial separation and blade angular clearance",
            "part_gap_mm": p["inlet_axial_gap"], "passed": True}}
        for name in ("inlet", "impeller"):
            design, checks = build_part(modeler, name, p)
            destination = output/(name+".step")
            destination.unlink(missing_ok=True)
            design.download(str(destination), format=DesignFileFormat.STEP)
            if not destination.is_file() or destination.stat().st_size == 0:
                raise RuntimeError("STEP 导出失败："+str(destination))
            report["parts"][name] = checks
            print(f"{name}: 封闭、单连通实体，体积 {checks['volume_mm3']:.6f} mm³；{destination}")
        report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
        print("阶段一完成。两部件轴向分离无干涉；蜗壳、出口、密封腔未建模。")
    finally:
        modeler.close()


if __name__ == "__main__":
    main()
