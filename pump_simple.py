"""阶段一 CFD 流体域：实心进口圆柱 + 圆盘减叶片和轴（长度单位 mm）。

python pump_simple.py --check-only
python pump_simple.py --output-dir output
本机启动沿用 main.connect() / launch_modeler_with_spaceclaim，后端为 241。
只用草图拉伸、圆周重复拉伸和布尔减法；不调用原生脚本，不构造金属零件。
"""

PARAMS = {
    "eye_diameter": 50.0,          # 进口眼/进口管共用的唯一流体直径
    "inlet_length": 80.0,          # 进口流体圆柱长度，Z=-length..0（含穿过前盖板的段）
    "impeller_diameter": 120.0,    # 两平面盖板之间的流体圆盘外径
    "hub_diameter": 30.0,          # 叶片根部定位圆直径；此阶段不扣除轮毂实体
    "shaft_diameter": 16.0,        # 叶轮域贯通轴孔直径，轴所占空间不属于流体
    "blade_count": 5,              # 圆周等分的叶片刀具数量（整数）
    "blade_thickness": 3.0,        # 直叶片法向厚度，单一矩形截面
    "blade_height": 12.0,          # 叶片高度及两盖板之间的流体域高度
    "blade_installation_angle_deg": 0.0, # 叶片直线相对根部径向的偏角；0 为径向直叶片
    "front_shroud_thickness": 3.0, # 前盖板厚度，仅定位金属边界 Z=-t..0，不生成金属
    "rear_shroud_thickness": 4.0,  # 后盖板厚度，仅定位金属边界 Z=h..h+t，不生成金属
    "blade_tool_overrun": 2.0,     # 叶片刀具超出圆盘外缘和上下边界的长度
    "section_angle_deg": 0.0,      # 含 Z 轴的剖切平面方位角
    "section_grid_count": 181,    # CAD 剖面填色采样网格边长（整数，非几何尺寸）
    "mesh_weld_tolerance": 1e-5,   # 三角化表面合并重合顶点的绝对距离
    "probe_clearance": 0.3,        # 验证探针距解析边界的最小距离，避开三角化边界误差
    "volume_relative_tolerance": 1e-5, # CAD/解析体积相对容差（无量纲）
    "volume_absolute_tolerance": 0.01, # CAD/解析体积绝对容差，mm³
    "volute_base_diameter": 132.0, # 后续预留：蜗壳基圆直径，本阶段不建模
    "volute_height": 24.0,         # 后续预留：蜗壳轴向高度
    "outlet_pipe_diameter": 40.0,  # 后续预留：出口管流体直径
    "outlet_pipe_length": 100.0,   # 后续预留：出口管长度
    "seal_chamber_diameter": 40.0, # 后续预留：密封腔腔径，轴径使用 shaft_diameter
    "seal_chamber_length": 20.0,   # 后续预留：密封腔长度
    "wear_ring_radial_clearance": 0.25, # 后续预留：口环单边径向间隙
    "wear_ring_length": 8.0,       # 后续预留：口环轴向长度
    "rear_chamber_gap": 2.0,       # 后续预留：后腔轴向间隙
}

import argparse
import json
from math import asin, cos, isfinite, pi, radians, sin, sqrt
import os
from pathlib import Path

from main import connect


def blade_layout(p):
    """单截面接口：返回每片的局部轴 d/n、矩形起点 q、横向偏移 b、刀具长度。

    根部中心为 (hub_radius, 0)，沿相对径向 beta 的直线向外。
    未来曲线叶片可替换此接口和 blade_sketch，默认仍为直线拉伸。
    """
    rh = p["hub_diameter"]/2
    beta = radians(p["blade_installation_angle_deg"])
    q, b = rh*cos(beta), -rh*sin(beta)
    length = p["impeller_diameter"]/2 + rh + p["blade_tool_overrun"]
    return [(2*pi*i/p["blade_count"]+beta, q, b, length)
            for i in range(p["blade_count"])]


def validate(p):
    signed = {"blade_installation_angle_deg", "section_angle_deg"}
    for key, value in p.items():
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not isfinite(value):
            raise ValueError(key+" 必须为有限数")
        if key not in signed and value <= 0:
            raise ValueError(key+" 必须为正数")
    if type(p["blade_count"]) is not int or p["blade_count"] < 2:
        raise ValueError("blade_count 必须为 >=2 的整数")
    if type(p["section_grid_count"]) is not int or p["section_grid_count"] < 21:
        raise ValueError("section_grid_count 必须为 >=21 的整数")
    rs, rh, re, ro = (p[k]/2 for k in ("shaft_diameter", "hub_diameter", "eye_diameter", "impeller_diameter"))
    a, c = p["blade_thickness"]/2, p["probe_clearance"]
    beta = abs(radians(p["blade_installation_angle_deg"]))
    if not rs < rh < re < ro:
        raise ValueError("需满足 轴径 < 叶片根部定位圆直径 < 进口眼直径 < 叶轮外径")
    if rh+a+c >= re or rh-a-c <= rs or a <= c:
        raise ValueError("叶片根部须在轴孔外、进口眼内，并留足探针净空")
    # 每个矩形刀具在自己的径向扇区内，扇区间留正间隙；保守限制避免刀具相交。
    if beta + asin((a+c)/rh) >= pi/p["blade_count"]:
        raise ValueError("安放角/厚度/叶片数导致叶片扇区交叠，请减小参数")
    if p["inlet_length"] <= p["front_shroud_thickness"]:
        raise ValueError("进口长度须大于前盖板厚度")
    if c >= min(rs, p["blade_height"]/2, (rh-a-rs)/2, (ro-re)/2):
        raise ValueError("probe_clearance 过大，无法放置内部验证点")
    if p["seal_chamber_diameter"] <= p["shaft_diameter"]:
        raise ValueError("预留密封腔径须大于轴径")
    return p


def circle_primitive(radius, y):
    """积分 sqrt(R²-y²) dy，单位 mm²。"""
    return (y*sqrt(radius*radius-y*y) + radius*radius*asin(y/radius))/2


def expected_volumes(p):
    validate(p)
    r, rs, a = p["impeller_diameter"]/2, p["shaft_diameter"]/2, p["blade_thickness"]/2
    _, q, b, _ = blade_layout(p)[0]
    # 所有根部角点都在圆内且在轴外；刀具相互不交。精确计圆盘内的矩形刀具面积。
    blade_area = circle_primitive(r, b+a)-circle_primitive(r, b-a)-2*a*q
    return {
        "inlet": pi*(p["eye_diameter"]/2)**2*p["inlet_length"],
        "impeller": (pi*(r*r-rs*rs)-p["blade_count"]*blade_area)*p["blade_height"],
    }


def blade_sketch(p, index):
    """仅单矩形截面，旋转草图坐标实现 24.1 圆周阵列，无 body.rotate。"""
    from ansys.geometry.core.math import Plane, Point2D, Point3D
    from ansys.geometry.core.misc import UNITS
    from ansys.geometry.core.sketch import Sketch
    theta, q, b, length = blade_layout(p)[index]
    sk = Sketch(Plane(origin=Point3D([0, 0, -p["blade_tool_overrun"]], unit=UNITS.mm),
                      direction_x=[cos(theta), sin(theta), 0],
                      direction_y=[-sin(theta), cos(theta), 0]))
    return sk.box(Point2D([q+length/2, b], unit=UNITS.mm),
                  length*UNITS.mm, p["blade_thickness"]*UNITS.mm)


def build_part(modeler, name, p):
    from ansys.geometry.core.math import Plane, Point2D, Point3D
    from ansys.geometry.core.misc import UNITS
    from ansys.geometry.core.sketch import Sketch
    design = modeler.create_design("pump_simple_"+name)

    def cylinder(label, z, diameter, length):
        sk = Sketch(Plane(origin=Point3D([0, 0, z], unit=UNITS.mm)))
        sk.circle(Point2D([0, 0], unit=UNITS.mm), diameter/2*UNITS.mm)
        return design.extrude_sketch(label, sk, distance=length*UNITS.mm)

    if name == "inlet":
        cylinder("InletFluid", -p["inlet_length"], p["eye_diameter"], p["inlet_length"])
    elif name == "impeller":
        fluid = cylinder("ImpellerFluid", 0, p["impeller_diameter"], p["blade_height"])
        extra = p["blade_tool_overrun"]
        for i in range(p["blade_count"]):
            tool = design.extrude_sketch("BladeTool_%02d" % (i+1), blade_sketch(p, i),
                                        distance=(p["blade_height"]+2*extra)*UNITS.mm)
            fluid.subtract(tool)  # 叶片为非流体，刀具被消耗
        shaft = cylinder("ShaftTool", -extra, p["shaft_diameter"], p["blade_height"]+2*extra)
        fluid.subtract(shaft)
    else:
        raise ValueError("未知部件："+name)
    return design


def mesh_inside(mesh, points):
    import numpy as np
    import pyvista as pv
    # 全部转为 mm；VTK 的封闭表面检查是本机 CAD 运行时的必要检查。
    result = pv.PolyData(np.asarray(points)).select_enclosed_points(mesh, tolerance=0.0, check_surface=True)
    return np.asarray(result["SelectedPoints"], dtype=bool)


def check_part(design, name, p):
    """24.1 公共 API 的实体/体积检查，以及 CAD 三角化表面的封闭性和连通性检查。"""
    bodies = design.get_all_bodies()
    if len(bodies) != 1 or bodies[0].is_surface:
        raise RuntimeError(name+" 应为一个流体实体且没有残余刀具")
    body = bodies[0]
    actual = body.volume.to("mm^3").magnitude
    target = expected_volumes(p)[name]
    if not isfinite(actual) or actual <= 0 or abs(actual-target) > max(
            p["volume_absolute_tolerance"], target*p["volume_relative_tolerance"]):
        raise RuntimeError(f"{name} CAD 体积 {actual} 不符合解析体积 {target} mm³")
    # 24.1 支持基础 tessellate；不传仅新版本支持的三角化选项。
    mesh = body.tessellate(merge=True).triangulate()
    mesh.points *= 1000.0  # Geometry Service 返回 SI 米；截面工具统一使用 mm
    mesh = mesh.clean(tolerance=p["mesh_weld_tolerance"], absolute=True)
    bad_edges = mesh.extract_feature_edges(boundary_edges=True, non_manifold_edges=True,
                                           feature_edges=False, manifold_edges=False)
    regions = mesh.connectivity()
    count = len(set(regions.cell_data["RegionId"])) if regions.n_cells else 0
    if bad_edges.n_cells or count != 1:
        raise RuntimeError(f"{name} 三角化表面未通过封闭单连通检查：异常边 {bad_edges.n_cells}，区域 {count}")
    return mesh, {"volume_mm3": actual, "expected_mm3": target,
                  "closed_mesh": True, "surface_regions": count, "nonmanifold_or_boundary_edges": 0}


def check_channels(mesh, p):
    """在 CAD 表面内外分类：轴孔/五片叶片为空，根部内环及每个扇区通道为流体。"""
    import numpy as np
    rh, rs, r = p["hub_diameter"]/2, p["shaft_diameter"]/2, p["impeller_diameter"]/2
    h, a, c = p["blade_height"], p["blade_thickness"]/2, p["probe_clearance"]
    points, expected = [[0, 0, h/2]], [False]
    for theta, q, b, length in blade_layout(p):
        # 叶片中线内部多个径向点（非流体），并检查厚度两侧的流体。
        tip = sqrt(r*r-b*b)
        for fraction in (0.2, 0.5, 0.8):
            s = q+(tip-q)*fraction
            for offset, fluid in ((0, False), (a+c, True), (-a-c, True)):
                v = b+offset
                x, y = s*cos(theta)-v*sin(theta), s*sin(theta)+v*cos(theta)
                if x*x+y*y < (r-c)**2:
                    points.append([x, y, h/2]); expected.append(fluid)
    # 轴孔外、叶片根部内的完整流体环，把所有叶片间通道连接到进口眼。
    inner = (rs+rh-a)/2
    for i in range(p["blade_count"]*4):
        angle = 2*pi*i/(p["blade_count"]*4)
        points.append([inner*cos(angle), inner*sin(angle), h/2]); expected.append(True)
    # 每个扇区边界上，直接验证进口眼至外缘的径向流体路径。
    for i in range(p["blade_count"]):
        angle = (2*i+1)*pi/p["blade_count"]
        for radius in (inner, p["eye_diameter"]/2, r-2*c):
            points.append([radius*cos(angle), radius*sin(angle), h/2]); expected.append(True)
    actual = mesh_inside(mesh, points)
    if not np.array_equal(actual, expected):
        raise RuntimeError("CAD 叶片空区/轴孔/内环通道探针检查失败")
    # 解析证明：每条叶片根部在 shaft 外，留完整内环；刀具在互不相交的扇区内且向外开口。
    return {"probe_count": len(points), "passed": True,
            "connectivity": "continuous inner annulus joins all inter-blade passages; disjoint blade sectors"}


def save_sections(mesh, p, output):
    """使用 PyVista 截面工具处理 CAD 三角化表面，同时输出轴中心剖面和中高截面。"""
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    angle = radians(p["section_angle_deg"])
    r, h, n = p["impeller_diameter"]/2, p["blade_height"], p["section_grid_count"]
    axial = mesh.slice(normal=[-sin(angle), cos(angle), 0], origin=[0, 0, 0])
    mid = mesh.slice(normal=[0, 0, 1], origin=[0, 0, h/2])
    if not axial.n_cells or not mid.n_cells:
        raise RuntimeError("轴中心/中高截面为空")
    axial.save(output/"impeller_axial_section.vtp")
    mid.save(output/"impeller_mid_section.vtp")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for ax, section, vertical in ((axes[0], axial, True), (axes[1], mid, False)):
        u = np.linspace(-r, r, n)
        v = np.linspace(0, h, n) if vertical else u
        uu, vv = np.meshgrid(u, v)
        if vertical:
            pts = np.c_[uu.ravel()*cos(angle), uu.ravel()*sin(angle), vv.ravel()]
        else:
            pts = np.c_[uu.ravel(), vv.ravel(), np.full(uu.size, h/2)]
        selected = mesh_inside(mesh, pts).reshape(uu.shape)
        ax.pcolormesh(uu, vv, selected, cmap="Blues", vmin=0, vmax=1, shading="nearest")
        # 实际 CAD slice 线段覆盖在分类填色上。
        lines, cursor = section.lines, 0
        while cursor < len(lines):
            k = lines[cursor]; ids = lines[cursor+1:cursor+1+k]; xyz = section.points[ids]
            x = xyz[:, 0]*cos(angle)+xyz[:, 1]*sin(angle) if vertical else xyz[:, 0]
            y = xyz[:, 2] if vertical else xyz[:, 1]
            ax.plot(x, y, color="black", linewidth=0.7)
            cursor += k+1
        ax.set_aspect("equal"); ax.set_xlabel("radial coordinate (mm)" if vertical else "X (mm)")
        ax.set_ylabel("Z (mm)" if vertical else "Y (mm)")
        ax.set_title("CAD axial section through Z axis" if vertical else "CAD mid-height section: blade/shaft voids")
    fig.suptitle("CFD domain: blue = fluid, white = excluded material")
    fig.savefig(output/"impeller_sections.png", dpi=180)
    plt.close(fig)
    return {"source": "CAD tessellation / PyVista slice", "axial_cells": axial.n_cells, "mid_cells": mid.n_cells}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default=os.getenv("ANSRV_GEO_HOST"))
    parser.add_argument("--port", type=int, default=int(os.getenv("ANSRV_GEO_PORT", "50051")))
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent/"output")
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    p = validate(dict(PARAMS))
    if args.check_only:
        print(json.dumps({"status": "offline_only_no_CAD", "expected_volumes_mm3": expected_volumes(p),
                          "blade_count": p["blade_count"], "interface_z_mm": 0.0,
                          "inlet_diameter_parameter": "eye_diameter"}, indent=2))
        return
    from ansys.geometry.core.connection.backend import BackendType
    from ansys.geometry.core.designer.design import DesignFileFormat
    modeler = connect(args.host, args.port)
    try:
        if not modeler.client.healthy or modeler.client.backend_type != BackendType.SPACECLAIM:
            raise RuntimeError("需要健康的 SPACECLAIM 后端")
        output = args.output_dir.resolve()
        output.mkdir(parents=True, exist_ok=True)
        report_path = output/"pump_simple_checks.json"
        report_path.unlink(missing_ok=True)
        report = {"params": p, "parts": {}, "interface": {
            "z_mm": 0, "gap_mm": 0, "overlap_volume_mm3": 0,
            "method": "inlet z<=0 and impeller z>=0; shared eye_diameter, no positive-volume intersection",
            "note": "interface is inlet disk minus shaft/blade-root footprints; full disk includes a shaft-end wall"}}
        for name in ("inlet", "impeller"):
            design = build_part(modeler, name, p)
            mesh, checks = check_part(design, name, p)
            # 直接从 CAD 返回包围盒复核轴向范围，避免只依赖参数推断。
            expected_z = (-p["inlet_length"], 0) if name == "inlet" else (0, p["blade_height"])
            if any(abs(a-b) > p["mesh_weld_tolerance"] for a, b in zip(mesh.bounds[4:6], expected_z)):
                raise RuntimeError(name+" CAD 轴向包围盒与共面接口不符")
            if name == "impeller":
                checks["channels"] = check_channels(mesh, p)
                checks["sections"] = save_sections(mesh, p, output)
            destination = output/(name+".step")
            destination.unlink(missing_ok=True)
            design.download(str(destination), format=DesignFileFormat.STEP)
            if not destination.is_file() or destination.stat().st_size == 0:
                raise RuntimeError("STEP 导出失败："+str(destination))
            report["parts"][name] = checks
            print(f"{name}: CFD 流体实体，体积 {checks['volume_mm3']:.6f} mm³；{destination}")
        report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
        print("阶段一流体域检查及导出完成；请查看 impeller_sections.png 的实际 CAD 剖面。")
    finally:
        modeler.close()


if __name__ == "__main__":
    main()
