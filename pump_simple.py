"""阶段一 CFD 流体域：环形进口 + 旋转包络减前后盖板、叶片及轴（长度单位 mm）。

python pump_simple.py --check-only
python pump_simple.py --output-dir output
本机启动沿用 main.connect() / launch_modeler_with_spaceclaim，后端为 241。
只用草图拉伸、圆周重复拉伸和布尔减法；不调用原生脚本；金属实体仅作布尔刀具，不导出。
"""

PARAMS = {
    "eye_diameter": 30.0,          # 进口眼/进口管共用的唯一流体直径
    "inlet_length": 80.0,          # 进口流体圆柱长度，Z=-length..0，止于前盖板外表面
    "impeller_diameter": 120.0,    # 两平面盖板之间的流体圆盘外径
    "hub_diameter": 30.0,          # 叶片根部定位圆直径；此阶段不扣除轮毂实体
    "shaft_diameter": 16.0,        # 进口和叶轮共用的贯通轴孔直径，轴所占空间不属于流体
    "blade_count": 5,              # 圆周等分的叶片刀具数量（整数）
    "blade_thickness": 3.0,        # 直叶片法向厚度，单一矩形截面
    "blade_height": 12.0,          # 叶片高度及两盖板之间的流体域高度
    "blade_installation_angle_deg": 0.0, # 叶片直线相对根部径向的偏角；0 为径向直叶片
    "front_shroud_thickness": 3.0, # 前盖板厚度，布尔刀具 Z=0..tf，保留眼部环形流体
    "rear_shroud_thickness": 4.0,  # 后盖板厚度，布尔刀具 Z=tf+h..tf+h+tb
    "blade_tool_overrun": 2.0,     # 刀具径向延伸、轴/开孔超切量；叶片轴向不跨盖板
    "section_angle_deg": 18.0,     # 默认避开 5 个径向叶片槽的轴截面角；0° 会沿 +X 叶片槽切开
    "section_grid_count": 181,    # CAD 剖面填色采样网格边长（整数，非几何尺寸）
    "mesh_weld_tolerance": 1e-5,   # 三角化表面合并重合顶点的绝对距离
    "alignment_tolerance": 0.01,   # 拟合轴心/半径和接口位置的绝对容差，mm
    "wall_normal_z_tolerance": 1e-6, # 提取圆柱侧壁时法向 Z 分量上限（无量纲）
    "mesh_inside_relative_tolerance": 1e-7, # VTK 内外分类相对于包围盒对角线的容差
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
    if not (rs < re < ro and rs < rh < ro):
        raise ValueError("轴径须小于 eye 与叶片根部定位圆直径，二者均须小于叶轮外径")
    # eye 可与根部定位圆同径；叶片从前盖板内表面开始，不会切坏眼部通道。
    if rh+a+c >= ro or rh-a-c <= rs or a <= c:
        raise ValueError("叶片根部须在轴孔外、叶轮外缘内，并留足探针净空")
    # 每个矩形刀具在自己的径向扇区内，扇区间留正间隙；保守限制避免刀具相交。
    if beta + asin((a+c)/rh) >= pi/p["blade_count"]:
        raise ValueError("安放角/厚度/叶片数导致叶片扇区交叠，请减小参数")
    if c >= min(rs, p["blade_height"]/2, (rh-a-rs)/2, (ro-re)/2, (re-rs)/4):
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
        "impeller": (pi*(r*r-rs*rs)-p["blade_count"]*blade_area)*p["blade_height"]
            + pi*((p["eye_diameter"]/2)**2-rs*rs)*p["front_shroud_thickness"],
    }


def blade_sketch(p, index):
    """仅单矩形截面，旋转草图坐标实现 24.1 圆周阵列，无 body.rotate。"""
    from ansys.geometry.core.math import Plane, Point2D, Point3D
    from ansys.geometry.core.misc import UNITS
    from ansys.geometry.core.sketch import Sketch
    theta, q, b, length = blade_layout(p)[index]
    sk = Sketch(Plane(origin=Point3D([0, 0, p["front_shroud_thickness"]], unit=UNITS.mm),
                      direction_x=[cos(theta), sin(theta), 0],
                      direction_y=[-sin(theta), cos(theta), 0]))
    return sk.box(Point2D([q+length/2, b], unit=UNITS.mm),
                  length*UNITS.mm, p["blade_thickness"]*UNITS.mm)


def build(modeler, p):
    """先生成唯一整体流体，再沿同一 Z=0 平面精确分区。"""
    from ansys.geometry.core.math import Plane, Point2D, Point3D
    from ansys.geometry.core.misc import UNITS
    from ansys.geometry.core.sketch import Sketch
    validate(p)
    design = modeler.create_design("pump_simple_stage1")

    def cylinder(label, z, diameter, length):
        sk = Sketch(Plane(origin=Point3D([0, 0, z], unit=UNITS.mm)))
        sk.circle(Point2D([0, 0], unit=UNITS.mm), diameter/2*UNITS.mm)
        return design.extrude_sketch(label, sk, distance=length*UNITS.mm)

    extra = p["blade_tool_overrun"]
    length = p["inlet_length"]
    tf, h, tb = p["front_shroud_thickness"], p["blade_height"], p["rear_shroud_thickness"]
    # 一个毛坯覆盖进口和完整旋转区域，没有进口/叶轮的拼接或并集。
    fluid = cylinder("InletFluid", -length, p["impeller_diameter"], length+tf+h+tb)
    # 将进口外侧非流体区与前盖板外环作为同一个连续套筒刀具，
    # 形成从 Z=-length 到 Z=tf 的同一解析 eye 圆柱侧壁，Z=0 不预留拼接面。
    sleeve = cylinder("InletExteriorAndFrontShroudTool", -length,
                      p["impeller_diameter"], length+tf)
    sleeve.subtract(cylinder("EyeOpeningTool", -length-extra,
                             p["eye_diameter"], length+tf+2*extra))
    fluid.subtract(sleeve)
    fluid.subtract(cylinder("RearShroudTool", tf+h, p["impeller_diameter"], tb))
    for i in range(p["blade_count"]):
        tool = design.extrude_sketch("BladeTool_%02d" % (i+1), blade_sketch(p, i),
                                    distance=h*UNITS.mm)
        fluid.subtract(tool)
    fluid.subtract(cylinder("ShaftTool", -length-extra, p["shaft_diameter"], length+tf+h+tb+2*extra))
    whole_bodies = design.get_all_bodies()
    if len(whole_bodies) != 1 or whole_bodies[0].is_surface:
        raise RuntimeError("分区前必须仅有一个整体流体实体，无残余刀具")
    targets = expected_volumes(p)
    whole_volume = fluid.volume.to("mm^3").magnitude
    expected = targets["inlet"]+targets["impeller"]
    allowed = max(p["volume_absolute_tolerance"], expected*p["volume_relative_tolerance"])
    if not isfinite(whole_volume) or abs(whole_volume-expected) > allowed:
        raise RuntimeError(f"整体流体体积异常：{whole_volume}，预期 {expected} mm³")
    inlet, impeller = partition_at_interface(design, fluid, p)
    total = sum(body.volume.to("mm^3").magnitude for body in (inlet, impeller))
    if abs(total-whole_volume) > allowed:
        raise RuntimeError("Z=0 分区前后体积不守恒")
    bodies = design.get_all_bodies()
    if len(bodies) != 2 or {body.name for body in bodies} != {"InletFluid", "ImpellerFluid"}:
        raise RuntimeError("平面分区后应只有进口和叶轮两个命名实体")
    print(f"整体流体 {whole_volume:.6f} mm³；Z=0 分区后体积守恒。")
    return design, {"inlet": inlet, "impeller": impeller}


def partition_at_interface(design, whole, p):
    """24.1 用布尔半空间求交实现精确 Z=0 分区，不调用要求 25R2 的 split_body。

    裁剪盒其他五个面均在整体流体外；只有共用 Z=0 平面切到流体。
    此方法得到几何重合的两张接口面，不承诺自动共享拓扑或共形网格。
    """
    from ansys.geometry.core.math import Plane, Point2D, Point3D
    from ansys.geometry.core.misc import UNITS
    from ansys.geometry.core.sketch import Sketch
    margin = p["blade_tool_overrun"]
    width = p["impeller_diameter"]+2*margin
    def halfspace(label, z, height):
        sk = Sketch(Plane(origin=Point3D([0, 0, z], unit=UNITS.mm)))
        sk.box(Point2D([0, 0], unit=UNITS.mm), width*UNITS.mm, width*UNITS.mm)
        return design.extrude_sketch(label, sk, distance=height*UNITS.mm)
    # 在任何裁剪前复制唯一的整体流体；不各自重建入口/叶轮。
    # 24.1 只在 extrude/copy 创建请求中指定名称；绝不调用 Body.name setter。
    impeller = whole.copy(design, "ImpellerFluid")
    whole.intersect(halfspace("BelowInterfaceTool", -p["inlet_length"]-margin,
                              p["inlet_length"]+margin))
    impeller.intersect(halfspace("AboveInterfaceTool", 0,
                                 p["front_shroud_thickness"]+p["blade_height"]+p["rear_shroud_thickness"]+margin))
    return whole, impeller


def mesh_inside(mesh, points):
    import numpy as np
    import pyvista as pv
    # 全部转为 mm；VTK 的封闭表面检查是本机 CAD 运行时的必要检查。
    result = pv.PolyData(np.asarray(points, dtype=float)).select_enclosed_points(
        mesh, tolerance=PARAMS["mesh_inside_relative_tolerance"], check_surface=True)
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
    zmid = p["front_shroud_thickness"]+h/2
    points, expected = [[0, 0, zmid]], [False]
    for theta, q, b, length in blade_layout(p):
        # 叶片中线内部多个径向点（非流体），并检查厚度两侧的流体。
        tip = sqrt(r*r-b*b)
        for fraction in (0.2, 0.5, 0.8):
            s = q+(tip-q)*fraction
            for offset, fluid in ((0, False), (a+c, True), (-a-c, True)):
                v = b+offset
                x, y = s*cos(theta)-v*sin(theta), s*sin(theta)+v*cos(theta)
                if x*x+y*y < (r-c)**2:
                    points.append([x, y, zmid]); expected.append(fluid)
    # 轴孔外、叶片根部内的完整流体环，把所有叶片间通道连接到进口眼。
    inner = (rs+rh-a)/2
    for i in range(p["blade_count"]*4):
        angle = 2*pi*i/(p["blade_count"]*4)
        points.append([inner*cos(angle), inner*sin(angle), zmid]); expected.append(True)
    # 每个扇区边界上，直接验证进口眼至外缘的径向流体路径。
    for i in range(p["blade_count"]):
        angle = (2*i+1)*pi/p["blade_count"]
        for radius in (inner, p["eye_diameter"]/2, r-2*c):
            points.append([radius*cos(angle), radius*sin(angle), zmid]); expected.append(True)
    actual = mesh_inside(mesh, points)
    if not np.array_equal(actual, expected):
        raise RuntimeError("CAD 叶片空区/轴孔/内环通道探针检查失败")
    # 解析证明：每条叶片根部在 shaft 外，留完整内环；刀具在互不相交的扇区内且向外开口。
    return {"probe_count": len(points), "passed": True,
            "connectivity": "continuous inner annulus joins all inter-blade passages; disjoint blade sectors"}


def save_png_safely(fig, destination):
    """先在内存渲染，再同目录临时写入并原子替换；路径/占用/渲染错误只记警告。"""
    from io import BytesIO
    from pathlib import PureWindowsPath
    import matplotlib.pyplot as plt
    destination = Path(destination)
    temporary = None
    try:
        if any(ord(ch) < 32 for ch in str(destination)) or destination.suffix.lower() != ".png":
            raise ValueError("PNG 路径含控制字符或扩展名错误")
        if os.name == "nt":
            windows_path = PureWindowsPath(destination)
            for part in windows_path.parts:
                if part == windows_path.anchor:
                    continue
                if (any(ch in '<>:"|?*' for ch in part) or part.endswith((" ", "."))
                        or PureWindowsPath(part).is_reserved()):
                    raise ValueError("Windows 路径组件不合法："+part)
        destination = destination.resolve()
        destination.parent.mkdir(parents=True, exist_ok=True)
        buffer = BytesIO()
        fig.savefig(buffer, format="png", dpi=180)
        # 不直接向目标 PNG 写入；查看器锁定目标时替换会给出明确的 OSError。
        with tempfile.NamedTemporaryFile(dir=destination.parent, prefix="pump_png_", suffix=".tmp", delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(buffer.getvalue())
        os.replace(temporary, destination)
        return {"saved": True, "path": str(destination)}
    except Exception as error:
        warning = {"saved": False, "path": str(destination), "error": str(error),
                   "errno": getattr(error, "errno", None),
                   "hint": "请检查路径、写权限和图片查看器占用；已有同名图片可能为旧版本。"}
        print("PNG 保存警告（不影响 STEP/几何检查）："+json.dumps(warning, ensure_ascii=False))
        return warning
    finally:
        if temporary is not None:
            try:
                temporary.unlink(missing_ok=True)
            except OSError:
                pass
        plt.close(fig)


def optional_plot(function, *args, **kwargs):
    """绘图/截线文件 IO 是可选输出，不混入几何验收的成功或失败判定。"""
    try:
        return function(*args, **kwargs)
    except Exception as error:
        import matplotlib.pyplot as plt
        plt.close("all")
        warning = {"saved": False, "error": str(error), "function": function.__name__}
        print("可选图像输出警告："+json.dumps(warning, ensure_ascii=False))
        return warning


def save_sections(mesh, p, output, name="impeller"):
    """使用 PyVista 截面工具处理 CAD 三角化表面，同时输出轴中心剖面和中高截面。"""
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    angle = radians(p["section_angle_deg"])
    r = p["impeller_diameter"]/2 if name == "impeller" else p["eye_diameter"]/2
    zlo, zhi = (0, p["front_shroud_thickness"]+p["blade_height"]) if name == "impeller" else (-p["inlet_length"], 0)
    n = p["section_grid_count"]
    axial = mesh.slice(normal=[-sin(angle), cos(angle), 0], origin=[0, 0, 0])
    mid_z = p["front_shroud_thickness"]+p["blade_height"]/2 if name == "impeller" else (zlo+zhi)/2
    mid = mesh.slice(normal=[0, 0, 1], origin=[0, 0, mid_z])
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
            pts = np.c_[uu.ravel(), vv.ravel(), np.full(uu.size, mid_z)]
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
    png = save_png_safely(fig, output/(name+"_sections.png"))
    return {"source": "CAD tessellation / PyVista slice", "axial_cells": axial.n_cells, "mid_cells": mid.n_cells,
            "png": png}


def fit_cylinder_section(mesh, radius_range, label, p, z_range=None):
    """用竖直侧壁的原始顶点拟合圆；避开三角面弦中点和封盖面内部顶点。

    这里只测量已生成的几何，不修改几何；圆柱的弦面片不应被误当成精确圆边界。
    """
    import numpy as np
    triangles = mesh.triangulate()
    ids = triangles.faces.reshape(-1, 4)[:, 1:]
    xyz = triangles.points[ids]
    normals = np.cross(xyz[:, 1]-xyz[:, 0], xyz[:, 2]-xyz[:, 0])
    norm = np.linalg.norm(normals, axis=1)
    vertical = (norm > 0) & (np.abs(normals[:, 2]) <= p["wall_normal_z_tolerance"]*norm)
    radii = np.linalg.norm(xyz[:, :, :2], axis=2)
    # 要求整个三角形都处于指定半径带，排除叶片壁面。
    selected = vertical & np.all((radii > radius_range[0]) & (radii < radius_range[1]), axis=1)
    if z_range is not None:
        selected &= np.all((xyz[:, :, 2] >= z_range[0]-p["alignment_tolerance"])
                           & (xyz[:, :, 2] <= z_range[1]+p["alignment_tolerance"]), axis=1)
    xy = np.unique(xyz[selected, :, :2].reshape(-1, 2), axis=0)
    if len(xy) < 3:
        raise RuntimeError(label+" 未找到足够圆柱侧壁顶点（可能缺孔、半径错误或明显偏轴）")
    origin = xy.mean(axis=0)
    local = xy-origin
    coefficients, _, rank, _ = np.linalg.lstsq(
        np.c_[2*local, np.ones(len(local))], np.sum(local*local, axis=1), rcond=None)
    if rank != 3:
        raise RuntimeError(label+" 圆柱顶点不足以确定圆")
    center = coefficients[:2]+origin
    radius = float(np.sqrt(coefficients[2]+np.dot(coefficients[:2], coefficients[:2])))
    residual = float(np.max(np.abs(np.linalg.norm(xy-center, axis=1)-radius)))
    if not np.isfinite(radius) or residual > p["alignment_tolerance"]:
        raise RuntimeError(f"{label} 非圆柱侧壁：拟合残差 {residual:.6g} mm")
    return {"center_mm": center.tolist(), "radius_mm": radius, "residual_mm": residual,
            "vertex_count": len(xy)}


def require_probes(mesh, points, expected, label):
    import numpy as np
    points = np.asarray(points, dtype=float).reshape(-1, 3)
    expected = np.asarray(expected, dtype=bool)
    actual = mesh_inside(mesh, points)
    failed = np.flatnonzero(actual != expected)
    if failed.size:
        details = [{"xyz_mm": points[i].tolist(), "expected_fluid": bool(expected[i]),
                    "actual_fluid": bool(actual[i])} for i in failed[:8]]
        raise RuntimeError(f"{label}: {failed.size}/{len(points)} 个探针失败；"+json.dumps(details))


def check_alignment(meshes, p):
    """圆柱拟合检查轴心/半径；深内部探针验证环形流体和跨接口通路。"""
    import numpy as np
    rs, re = p["shaft_diameter"]/2, p["eye_diameter"]/2
    length, height, tol = p["inlet_length"], p["front_shroud_thickness"]+p["blade_height"], p["alignment_tolerance"]
    rh, a = p["hub_diameter"]/2, p["blade_thickness"]/2
    for name, mesh in meshes.items():
        expected_z = (-length, 0) if name == "inlet" else (0, height)
        if any(abs(x-y) > tol for x, y in zip(mesh.bounds[4:6], expected_z)):
            raise RuntimeError(f"{name} CAD 轴向包围盒 {list(mesh.bounds[4:6])} 与接口 {expected_z} 不符")
    # 轴孔壁位于根部内侧；eye 壁在前盖板厚度范围内，可独立拟合校验。
    inner = (rs+rh-a)/2
    fits = {
        "inlet_shaft": fit_cylinder_section(meshes["inlet"], (0, inner), "进口轴孔", p),
        "impeller_shaft": fit_cylinder_section(meshes["impeller"], (0, inner), "叶轮轴孔", p),
        "inlet_eye": fit_cylinder_section(meshes["inlet"], (inner, re+(re-rs)/2), "进口外圆", p),
        "impeller_eye": fit_cylinder_section(meshes["impeller"], (inner, re+(re-rs)/2),
                                              "前盖板进口眼", p, (0, p["front_shroud_thickness"])),
    }
    for name, fit in fits.items():
        expected_radius = re if name.endswith("_eye") else rs
        if np.linalg.norm(fit["center_mm"]) > tol or abs(fit["radius_mm"]-expected_radius) > tol:
            raise RuntimeError(f"{name} 轴心/半径不符：{fit}，预期轴心 [0,0]、半径 {expected_radius} mm，容差 {tol} mm")
    left, right = fits["inlet_shaft"], fits["impeller_shaft"]
    if (np.linalg.norm(np.array(left["center_mm"])-right["center_mm"]) > tol
            or abs(left["radius_mm"]-right["radius_mm"]) > tol):
        raise RuntimeError("两域轴孔拟合轴心/半径不一致："+json.dumps(fits))
    # 不在 rs±0.3、re±0.3 处分类：粗三角化的弦面会使这些点跨越网格边界。
    # 半径/同轴由上面的拟合验证；内外分类仅用远离边界的代表性点。
    points, expected = [], []
    ring_middle = (rs+re)/2
    for fraction in (0.25, 0.5, 0.75):
        z = -length*fraction
        points.append([0.0, 0.0, z]); expected.append(False)
        for i in range(p["blade_count"]*4):
            angle = 2*pi*(i+0.5)/(p["blade_count"]*4)
            for radius, fluid in ((rs/2, False), (ring_middle, True), (re+(re-rs)/2, False)):
                points.append([radius*cos(angle), radius*sin(angle), z]); expected.append(fluid)
    require_probes(meshes["inlet"], points, expected, "进口环形流体/轴孔检查")
    # 采样位置由半径与角度在全局 XY 中直接构造，两个域用完全相同的 XY。
    xy = []
    for i in range(p["blade_count"]):
        angle = (2*i+1)*pi/p["blade_count"]
        for radius in (inner, (inner+re)/2, (inner+3*re)/4):
            xy.append([radius*cos(angle), radius*sin(angle)])
    dz = min(p["front_shroud_thickness"]/2, length/4)  # 离开接口三角面，不在 Z=±0.3 的近边界层采样
    for name, z in (("inlet", -dz), ("impeller", dz)):
        probes = [[x, y, z] for x, y in xy]
        require_probes(meshes[name], probes, [True]*len(probes), name+" 接口通流检查")
        require_probes(meshes[name], [[0, 0, z]], [False], name+" 贯通轴孔检查")
    return {"passed": True, "inlet_probe_count": len(points), "throughflow_paths": len(xy),
            "cylinder_fits": fits, "alignment_tolerance_mm": tol,
            "interface_z_mm": 0, "axial_gap_mm": float(meshes["impeller"].bounds[4]-meshes["inlet"].bounds[5]),
            "diameter_parameters": ["eye_diameter", "shaft_diameter"],
            "method": "cylindrical sidewall fits and deep-interior probes; shared XY, Z-axis and eye parameter",
            "note": "full annular interface at Z=0; actual front-shroud eye wall fitted over 0..tf"}


def check_interface_caps(meshes, p):
    """核对 Z=0 两个实际三角化环形端面的面积和内外半径。"""
    import numpy as np
    caps = {}
    for name, mesh in meshes.items():
        xyz = mesh.points[mesh.faces.reshape(-1, 4)[:, 1:]]
        chosen = np.all(np.abs(xyz[:, :, 2]) <= p["alignment_tolerance"], axis=1)
        faces = xyz[chosen]
        if not len(faces):
            raise RuntimeError(name+" 缺少 Z=0 接口面")
        normals = np.cross(faces[:, 1]-faces[:, 0], faces[:, 2]-faces[:, 0])
        area = float(np.linalg.norm(normals, axis=1).sum()/2)
        radii = np.linalg.norm(faces[:, :, :2], axis=2)
        caps[name] = {"area_mm2": area, "inner_radius_mm": float(radii.min()),
                      "outer_radius_mm": float(radii.max())}
        for actual, expected in ((radii.min(), p["shaft_diameter"]/2),
                                 (radii.max(), p["eye_diameter"]/2)):
            if abs(actual-expected) > p["alignment_tolerance"]:
                raise RuntimeError(name+" 接口内外圆不匹配")
    allowed = max(p["alignment_tolerance"]*pi*p["eye_diameter"], p["volume_absolute_tolerance"])
    if abs(caps["inlet"]["area_mm2"]-caps["impeller"]["area_mm2"]) > allowed:
        raise RuntimeError("两侧 Z=0 环形接口面积不一致："+json.dumps(caps))
    return {"passed": True, "caps": caps, "area_tolerance_mm2": allowed,
            "shared_plane": "Z=0; geometric coincidence, no automatic shared topology"}


def check_impeller_footprint(mesh, p):
    """中高平面按完整圆盘减轴/矩形槽的解析定义检查，不强求奇数叶片轴截面镜像对称。"""
    import numpy as np
    r, rs, a = p["impeller_diameter"]/2, p["shaft_diameter"]/2, p["blade_thickness"]/2
    u = np.linspace(-r, r, p["section_grid_count"])
    xx, yy = np.meshgrid(u, u)
    x, y = xx.ravel(), yy.ravel()
    radius = np.hypot(x, y)
    expected = (radius < r) & (radius > rs)
    margin = np.minimum(np.abs(radius-r), np.abs(radius-rs))
    for theta, q, b, length in blade_layout(p):
        along, across = x*cos(theta)+y*sin(theta), -x*sin(theta)+y*cos(theta)-b
        slot = (along >= q) & (along <= q+length) & (np.abs(across) <= a)
        expected &= ~slot
        outside = np.hypot(np.maximum(np.maximum(q-along, along-q-length), 0),
                           np.maximum(np.abs(across)-a, 0))
        inside = np.minimum(np.minimum(along-q, q+length-along), a-np.abs(across))
        margin = np.minimum(margin, np.where(slot, inside, outside))
    points = np.c_[x, y, np.full(x.size, p["front_shroud_thickness"]+p["blade_height"]/2)]
    # 同时避开解析边界和三角化边界，区别面片近似误差与大块流体缺失。
    _, closest = mesh.find_closest_cell(points, return_closest_point=True)
    safe = (margin > p["probe_clearance"]) & (np.linalg.norm(points-closest, axis=1) > p["probe_clearance"])
    if not safe.any():
        raise RuntimeError("叶轮平面网格不足以检查圆盘包络和叶片槽")
    require_probes(mesh, points[safe], expected[safe], "叶轮完整圆盘/轴孔/阵列槽平面检查")
    return {"passed": True, "probe_count": int(safe.sum()), "blade_count": p["blade_count"],
            "definition": "full circular disk minus shaft and equally spaced straight blade slots"}


def check_shrouds(mesh, p):
    """前盖板眼内通流、眼外封闭、后盖板封闭；眼部不应被叶片超切。"""
    from math import cos, sin, pi
    tf, h, tb = p["front_shroud_thickness"], p["blade_height"], p["rear_shroud_thickness"]
    rs, re, r = p["shaft_diameter"]/2, p["eye_diameter"]/2, p["impeller_diameter"]/2
    points, expected = [], []
    for i in range(p["blade_count"]*4):
        angle = 2*pi*i/(p["blade_count"]*4)
        for radius, fluid in ((rs/2, False), ((rs+re)/2, True), ((re+r)/2, False)):
            points.append([radius*cos(angle), radius*sin(angle), tf/2]); expected.append(fluid)
        # 后盖板原体积中没有流体。
        points.append([(rs+r)/2*cos(angle), (rs+r)/2*sin(angle), tf+h+tb/2]); expected.append(False)
    require_probes(mesh, points, expected, "前后盖板/进口眼壁面检查")
    return {"passed": True, "probe_count": len(points), "front_eye_z_mm": [0, tf],
            "blade_channel_z_mm": [tf, tf+h], "rear_shroud_z_mm": [tf+h, tf+h+tb]}


def save_alignment_section(meshes, p, output):
    """实际填色+截线，轴截面旁显示完整中高截面，明确槽引起的非对称。"""
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    theta = radians(p["section_angle_deg"])
    r, h, n = p["impeller_diameter"]/2, p["blade_height"], p["section_grid_count"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 8), constrained_layout=True)
    ax = axes[0]
    for name, color in (("inlet", "tab:blue"), ("impeller", "tab:orange")):
        lo, hi = (-p["inlet_length"], 0) if name == "inlet" else (0, p["front_shroud_thickness"]+h)
        uu, zz = np.meshgrid(np.linspace(-r, r, n), np.linspace(lo, hi, n))
        pts = np.c_[uu.ravel()*cos(theta), uu.ravel()*sin(theta), zz.ravel()]
        fluid = mesh_inside(meshes[name], pts).reshape(uu.shape)
        ax.pcolormesh(uu, zz, np.ma.masked_where(~fluid, fluid.astype(float)),
                      cmap=ListedColormap([color]), shading="nearest", vmin=0, vmax=1)
        section = meshes[name].slice(normal=[-sin(theta), cos(theta), 0], origin=[0, 0, 0])
        if not section.n_cells:
            raise RuntimeError(name+" 组合剖面的轴中心截线为空")
        section.save(output/(name+"_assembly_axial_section.vtp"))
        lines, cursor = section.lines, 0
        while cursor < len(lines):
            k = lines[cursor]; ids = lines[cursor+1:cursor+1+k]; xyz = section.points[ids]
            ax.plot(xyz[:, 0]*cos(theta)+xyz[:, 1]*sin(theta), xyz[:, 2], color="black", linewidth=0.6)
            cursor += k+1
        ax.plot([], [], color=color, label=name)
    ax.axhline(0, color="gray", linestyle="--", label="interface Z=0")
    ax.axhline(p["front_shroud_thickness"], color="gray", linestyle=":", label="front inner wall Z=tf")
    ax.set(xlim=(-r, r), xlabel="signed radial coordinate (mm)", ylabel="Z (mm)",
           title="Actual axial section at %.1f degrees" % p["section_angle_deg"])
    ax.set_aspect("equal"); ax.legend()
    ax = axes[1]
    xx, yy = np.meshgrid(np.linspace(-r, r, n), np.linspace(-r, r, n))
    pts = np.c_[xx.ravel(), yy.ravel(), np.full(xx.size, p["front_shroud_thickness"]+h/2)]
    fluid = mesh_inside(meshes["impeller"], pts).reshape(xx.shape)
    ax.pcolormesh(xx, yy, fluid, cmap="Blues", vmin=0, vmax=1, shading="nearest")
    t = np.linspace(0, 2*pi, n)
    ax.plot(r*np.cos(t), r*np.sin(t), "k--", linewidth=0.8, label="disk envelope")
    ax.plot([-r*cos(theta), r*cos(theta)], [-r*sin(theta), r*sin(theta)],
            color="red", linestyle="--", label="axial cut direction")
    ax.set(xlabel="X (mm)", ylabel="Y (mm)", title="Actual mid-height fluid / blade and shaft voids")
    ax.set_aspect("equal"); ax.legend()
    fig.suptitle("CFD assembly: actual sections, no mirroring or filling of blade voids")
    return save_png_safely(fig, output/"pump_alignment_section.png")


def export_parts(modeler, design, p, output):
    """从同一个联合模型快照分离导出；不重新建模，不跨设计复制，不依赖新版抑制 API。"""
    from ansys.geometry.core.designer.design import DesignFileFormat
    assembly = output/"pump_assembly.step"
    assembly.unlink(missing_ok=True)
    design.download(str(assembly), format=DesignFileFormat.STEP)
    if not assembly.is_file() or not assembly.stat().st_size:
        raise RuntimeError("装配 STEP 导出失败")
    print("已导出双流体域装配（不合并）："+str(assembly))
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
            # 这里只做必要的导出隔离；所有几何/网格/探针检查均放在两个 STEP 写出后。
            destination = output/(name+".step")
            destination.unlink(missing_ok=True)
            isolated.download(str(destination), format=DesignFileFormat.STEP)
            if not destination.is_file() or not destination.stat().st_size:
                raise RuntimeError("STEP 导出失败："+str(destination))
            print("已独立导出："+str(destination))
        # open_file 关闭过原设计；恢复联合快照后才进行后续检查。
        restored = modeler.open_file(str(snapshot))
        lookup = {body.name: body for body in restored.get_all_bodies()}
        return restored, {"inlet": lookup["InletFluid"], "impeller": lookup["ImpellerFluid"]}


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
        design, bodies = export_parts(modeler, design, p, output)
        report, meshes = {"params": p, "status": "exported_checks_pending",
                          "step_exports": ["pump_assembly.step", "inlet.step", "impeller.step"],
                          "construction": "single envelope minus solids; exact Z=0 halfspace Boolean partition (24.1)",
                          "parts": {}}, {}
        try:
            for name, body in bodies.items():
                mesh, checks = check_part(body, name, p)
                meshes[name] = mesh
                if name == "impeller":
                    checks["channels"] = check_channels(mesh, p)
                    checks["footprint"] = check_impeller_footprint(mesh, p)
                    checks["shrouds"] = check_shrouds(mesh, p)
                checks["sections"] = optional_plot(save_sections, mesh, p, output, name)
                report["parts"][name] = checks
                print(f"{name}: CFD 流体实体，体积 {checks['volume_mm3']:.6f} mm³")
            report["interface"] = check_alignment(meshes, p)
            report["interface"]["cap_match"] = check_interface_caps(meshes, p)
            report["alignment_plot"] = optional_plot(save_alignment_section, meshes, p, output)
            report["status"] = "exported_checks_passed"
        except Exception as error:
            report["status"] = "exported_checks_failed"
            report["error"] = str(error)
            print("三个 STEP 已导出并保留；后续检查失败，详情见 pump_simple_checks.json。")
            raise
        finally:
            report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
        print("STEP 导出及几何检查完成；图片保存状态见 pump_simple_checks.json。")
    finally:
        modeler.close()


if __name__ == "__main__":
    main()
