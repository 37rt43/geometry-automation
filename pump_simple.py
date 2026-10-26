"""阶段一 CFD 流体域：环形进口流体 + 圆盘减叶片和轴（长度单位 mm）。

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
    "shaft_diameter": 16.0,        # 进口和叶轮共用的贯通轴孔直径，轴所占空间不属于流体
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
import tempfile

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
        "inlet": pi*((p["eye_diameter"]/2)**2-rs*rs)*p["inlet_length"],
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


def build(modeler, p):
    """一次建出同一设计中的两个同轴域，保持独立实体供静止/旋转域分别导出。"""
    from ansys.geometry.core.math import Plane, Point2D, Point3D
    from ansys.geometry.core.misc import UNITS
    from ansys.geometry.core.sketch import Sketch
    validate(p)
    design = modeler.create_design("pump_simple_stage1")

    def cylinder(label, z, diameter, length):
        sk = Sketch(Plane(origin=Point3D([0, 0, z], unit=UNITS.mm)))
        sk.circle(Point2D([0, 0], unit=UNITS.mm), diameter/2*UNITS.mm)
        return design.extrude_sketch(label, sk, distance=length*UNITS.mm)

    # 统一的 XY 原点、Z 轴、接口 Z=0；进口外径/轴孔直接复用叶轮参数。
    extra = p["blade_tool_overrun"]
    inlet = cylinder("InletFluid", -p["inlet_length"], p["eye_diameter"], p["inlet_length"])
    inlet.subtract(cylinder("InletShaftTool", -p["inlet_length"]-extra,
                            p["shaft_diameter"], p["inlet_length"]+2*extra))

    # 以下叶轮圆盘、叶片刀具及轴孔的尺寸/位置/减法逻辑保持不变。
    fluid = cylinder("ImpellerFluid", 0, p["impeller_diameter"], p["blade_height"])
    for i in range(p["blade_count"]):
        tool = design.extrude_sketch("BladeTool_%02d" % (i+1), blade_sketch(p, i),
                                    distance=(p["blade_height"]+2*extra)*UNITS.mm)
        fluid.subtract(tool)
    shaft = cylinder("ShaftTool", -extra, p["shaft_diameter"], p["blade_height"]+2*extra)
    fluid.subtract(shaft)
    bodies = design.get_all_bodies()
    if len(bodies) != 2 or {body.name for body in bodies} != {"InletFluid", "ImpellerFluid"}:
        raise RuntimeError("联合建模后应只有进口、叶轮两个命名流体实体，无残余刀具")
    return design, {"inlet": inlet, "impeller": fluid}


def mesh_inside(mesh, points):
    import numpy as np
    import pyvista as pv
    # 全部转为 mm；VTK 的封闭表面检查是本机 CAD 运行时的必要检查。
    result = pv.PolyData(np.asarray(points)).select_enclosed_points(mesh, tolerance=0.0, check_surface=True)
    return np.asarray(result["SelectedPoints"], dtype=bool)


def check_part(body, name, p):
    """24.1 公共 API 的实体/体积检查，以及 CAD 三角化表面的封闭性和连通性检查。"""
    if body.is_surface:
        raise RuntimeError(name+" 应为流体实体而非曲面")
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


def save_sections(mesh, p, output, name="impeller"):
    """使用 PyVista 截面工具处理 CAD 三角化表面，同时输出轴中心剖面和中高截面。"""
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    angle = radians(p["section_angle_deg"])
    r = p["impeller_diameter"]/2 if name == "impeller" else p["eye_diameter"]/2
    zlo, zhi = (0, p["blade_height"]) if name == "impeller" else (-p["inlet_length"], 0)
    n = p["section_grid_count"]
    axial = mesh.slice(normal=[-sin(angle), cos(angle), 0], origin=[0, 0, 0])
    mid = mesh.slice(normal=[0, 0, 1], origin=[0, 0, (zlo+zhi)/2])
    if not axial.n_cells or not mid.n_cells:
        raise RuntimeError("轴中心/中高截面为空")
    axial.save(output/(name+"_axial_section.vtp"))
    mid.save(output/(name+"_mid_section.vtp"))
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for ax, section, vertical in ((axes[0], axial, True), (axes[1], mid, False)):
        u = np.linspace(-r, r, n)
        v = np.linspace(zlo, zhi, n) if vertical else u
        uu, vv = np.meshgrid(u, v)
        if vertical:
            pts = np.c_[uu.ravel()*cos(angle), uu.ravel()*sin(angle), vv.ravel()]
        else:
            pts = np.c_[uu.ravel(), vv.ravel(), np.full(uu.size, (zlo+zhi)/2)]
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
        ax.set_title(name+": "+("CAD axial section through Z axis" if vertical else "CAD transverse section"))
    fig.suptitle("CFD domain: blue = fluid, white = excluded material")
    fig.savefig(output/(name+"_sections.png"), dpi=180)
    plt.close(fig)
    return {"source": "CAD tessellation / PyVista slice", "axial_cells": axial.n_cells, "mid_cells": mid.n_cells}


def check_alignment(meshes, p):
    """两域包围盒贴合；轴孔/外圆内外探针；Z=0 两侧通流探针。"""
    import numpy as np
    rs, re = p["shaft_diameter"]/2, p["eye_diameter"]/2
    c, length, height = p["probe_clearance"], p["inlet_length"], p["blade_height"]
    for name, mesh in meshes.items():
        expected_z = (-length, 0) if name == "inlet" else (0, height)
        if any(abs(a-b) > p["mesh_weld_tolerance"] for a, b in zip(mesh.bounds[4:6], expected_z)):
            raise RuntimeError(name+" CAD 轴向包围盒与共面接口不符")
    # 圆周多个方向、多处轴向站点，确认进口中心始终空、内外径同轴且无误挖/堵塞。
    points, expected = [], []
    for fraction in (0.25, 0.5, 0.75):
        z = -length*fraction
        points.append([0, 0, z]); expected.append(False)
        for i in range(p["blade_count"]*4):
            angle = 2*pi*i/(p["blade_count"]*4)
            for radius, fluid in ((rs-c, False), (rs+c, True), (re-c, True), (re+c, False)):
                points.append([radius*cos(angle), radius*sin(angle), z]); expected.append(fluid)
    if not np.array_equal(mesh_inside(meshes["inlet"], points), expected):
        raise RuntimeError("进口轴孔/环形流体截面/外径探针检查失败")
    # 检查同一 XY 点在接口两侧都位于流体内：内环及所有叶间通道。
    rh, a = p["hub_diameter"]/2, p["blade_thickness"]/2
    inner = (rs+rh-a)/2
    xy = []
    for i in range(p["blade_count"]):
        angle = (2*i+1)*pi/p["blade_count"]
        for radius in (inner, (inner+re)/2, re-c):
            xy.append([radius*cos(angle), radius*sin(angle)])
    dz = min(c, length/4, height/4)
    for name, z in (("inlet", -dz), ("impeller", dz)):
        probes = [[x, y, z] for x, y in xy]
        if not mesh_inside(meshes[name], probes).all():
            raise RuntimeError(name+" 在 Z=0 接口附近的通流探针失败")
        # 两侧同轴贯通轴孔，避免只检查进口轴心而漏掉叶轮孔径不匹配。
        for i in range(p["blade_count"]*4):
            angle = 2*pi*i/(p["blade_count"]*4)
            pair = [[radius*cos(angle), radius*sin(angle), z] for radius in (rs-c, rs+c)]
            if not np.array_equal(mesh_inside(meshes[name], pair), [False, True]):
                raise RuntimeError("接口两侧轴孔径/同轴探针失败")
    return {"passed": True, "inlet_probe_count": len(points), "throughflow_paths": len(xy),
            "interface_z_mm": 0, "axial_gap_mm": 0, "overlap_volume_mm3": 0,
            "diameter_parameters": ["eye_diameter", "shaft_diameter"],
            "method": "CAD bounds and probes on both sides; disjoint axial interiors; connected inner annulus",
            "note": "shared annular eye minus blade-root footprints; blade-root end faces remain walls"}


def save_alignment_section(meshes, p, output):
    """将两域实际 CAD 的轴中心截线绘在同一坐标系，核对贴合和贯通轴孔。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    theta = radians(p["section_angle_deg"])
    fig, ax = plt.subplots(figsize=(7, 8), constrained_layout=True)
    for name, color in (("inlet", "tab:blue"), ("impeller", "tab:orange")):
        section = meshes[name].slice(normal=[-sin(theta), cos(theta), 0], origin=[0, 0, 0])
        if not section.n_cells:
            raise RuntimeError(name+" 组合剖面的轴中心截线为空")
        lines, cursor = section.lines, 0
        while cursor < len(lines):
            k = lines[cursor]; ids = lines[cursor+1:cursor+1+k]; xyz = section.points[ids]
            ax.plot(xyz[:, 0]*cos(theta)+xyz[:, 1]*sin(theta), xyz[:, 2], color=color)
            cursor += k+1
        ax.plot([], [], color=color, label=name)
    ax.axhline(0, color="gray", linestyle="--", label="shared interface Z=0")
    ax.set(xlabel="radial coordinate (mm)", ylabel="Z (mm)", title="CAD axial section: annular inlet + impeller")
    ax.set_aspect("equal"); ax.legend()
    fig.savefig(output/"pump_alignment_section.png", dpi=180)
    plt.close(fig)


def export_parts(modeler, design, p, output):
    """从同一个联合模型快照分离导出；不重新建模，不跨设计复制，不依赖新版抑制 API。"""
    from ansys.geometry.core.designer.design import DesignFileFormat
    with tempfile.TemporaryDirectory(prefix="pump_export_") as folder:
        snapshot = Path(folder)/"stage1.scdocx"
        design.download(str(snapshot), format=DesignFileFormat.SCDOCX)
        if not snapshot.is_file() or not snapshot.stat().st_size:
            raise RuntimeError("联合模型快照保存失败")
        for name, label in (("inlet", "InletFluid"), ("impeller", "ImpellerFluid")):
            # open_file 会关闭当前设计；每次均从同一未删减的联合快照恢复。
            isolated = modeler.open_file(str(snapshot))
            bodies = isolated.get_all_bodies()
            if len(bodies) != 2 or {b.name for b in bodies} != {"InletFluid", "ImpellerFluid"}:
                raise RuntimeError("联合快照重载后未保留两个命名域")
            for body in bodies:
                if body.name != label:
                    isolated.delete_body(body)
            remaining = isolated.get_all_bodies()
            if len(remaining) != 1 or remaining[0].name != label:
                raise RuntimeError(name+" 导出隔离失败")
            # 快照重载后再次检查实体和体积，确保导出对应同一几何。
            check_part(remaining[0], name, p)
            destination = output/(name+".step")
            destination.unlink(missing_ok=True)
            isolated.download(str(destination), format=DesignFileFormat.STEP)
            if not destination.is_file() or not destination.stat().st_size:
                raise RuntimeError("STEP 导出失败："+str(destination))
            print("已独立导出："+str(destination))


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
    modeler = connect(args.host, args.port)
    try:
        if not modeler.client.healthy or modeler.client.backend_type != BackendType.SPACECLAIM:
            raise RuntimeError("需要健康的 SPACECLAIM 后端")
        output = args.output_dir.resolve()
        output.mkdir(parents=True, exist_ok=True)
        report_path = output/"pump_simple_checks.json"
        report_path.unlink(missing_ok=True)
        design, bodies = build(modeler, p)
        report, meshes = {"params": p, "parts": {}}, {}
        for name, body in bodies.items():
            mesh, checks = check_part(body, name, p)
            meshes[name] = mesh
            if name == "impeller":
                checks["channels"] = check_channels(mesh, p)
            checks["sections"] = save_sections(mesh, p, output, name)
            report["parts"][name] = checks
            print(f"{name}: CFD 流体实体，体积 {checks['volume_mm3']:.6f} mm³")
        report["interface"] = check_alignment(meshes, p)
        save_alignment_section(meshes, p, output)
        export_parts(modeler, design, p, output)
        report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
        print("阶段一流体域检查及导出完成；请查看 inlet_sections.png、impeller_sections.png 和 pump_alignment_section.png。")
    finally:
        modeler.close()


if __name__ == "__main__":
    main()
