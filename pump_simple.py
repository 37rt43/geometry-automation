"""参数化离心泵 CFD：整体进口/叶轮 + 变截面梯形蜗壳与切向扩散段（长度单位 mm）。

python pump_simple.py --check-only
python pump_simple.py --output-dir output
本机启动沿用 main.connect() / launch_modeler_with_spaceclaim，后端为 241。
核心用基本特征；仅蜗壳按用户授权调用 V241 原生平滑放样。金属刀具不导出。
"""

PARAMS = {
    "eye_diameter": 30.0,          # 进口眼/进口管共用的唯一流体直径
    "inlet_length": 80.0,          # 进口流体圆柱长度，Z=-length..0，止于前盖板外表面
    "impeller_diameter": 120.0,    # 两平面盖板之间的流体圆盘外径
    "hub_diameter": 30.0,          # 叶片根部定位圆直径；此阶段不扣除轮毂实体
    "shaft_diameter": 16.0,        # 进口和叶轮共用的贯通轴孔直径，轴所占空间不属于流体
    "blade_count": 5,              # 圆周等分的叶片刀具数量（整数）
    "blade_thickness": 3.0,        # 外圆弧处的周向弧长；扇区角宽=此值/外半径
    "blade_height": 12.0,          # 叶片高度及两盖板之间的流体域高度
    "blade_installation_angle_deg": 0.0, # 环形扇区相对等分方位的整体偏角，两侧仍为径向直线
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
    "volute_base_diameter": 132.0, # 主蜗壳基圆；Ri..Rb 为均匀环形流体间隙
    "volute_section_shape": "trapezoid", # 本版固定使用非矩形梯形，四角对应一致
    "volute_start_width": 12.0,    # 隔舌邻近起始截面径向宽度
    "volute_end_width": 40.0,      # 扩散前终了截面径向宽度
    "volute_axial_width": 32.0,    # 梯形外侧轴向高度；内侧高度共用 blade_height
    "volute_wrap_deg": 330.0,      # 主螺旋包角；从 tongue_position_deg 起计
    "volute_section_count": 57,   # 平滑放样的主蜗壳截面数，非分段棱柱
    "volute_connection_overlap": 1.0, # 均匀间隙桥与核心/蜗壳的布尔重叠量
    "tongue_position_deg": 0.0,    # 隔舌和螺旋起点的绝对方位
    "tongue_orientation_deg": 0.0, # 隔舌长轴相对起点径向的偏角
    "tongue_nose_radius": 3.0,     # 隔舌鼻部精确半圆半径，不使用尖点
    "tongue_tail_length": 160.0,   # 隔舌刀具尾部延伸到整个蜗壳之外
    "diffuser_length": 30.0,       # 切向扩散段长度，不含任何外接直管
    "diffuser_angle_deg": 8.0,     # 径向总扩散角；两侧各 angle/2
    "diffuser_blend_length": 6.0,  # 入口壁斜率从零平滑增至目标角的长度
    "diffuser_section_count": 13,  # 扩散段放样截面数；与主蜗壳一次放样
    "volute_min_edge_length": 0.05, # CAD 最短边下限 mm，用于定位细碎尖刺特征
    "volute_min_face_area": 0.05,  # CAD 面积下限 mm²，小于此值需排查碎面
    "volute_area_relative_tolerance": 0.08, # 截面采样面积与设计面积相对容差
    "volute_check_grid_count": 65, # 实际截面填色/面积抽查网格边长
    "seal_chamber_diameter": 40.0, # 后续预留：密封腔腔径，轴径使用 shaft_diameter
    "seal_chamber_length": 20.0,   # 后续预留：密封腔长度
    "wear_ring_radial_clearance": 0.25, # 后续预留：口环单边径向间隙
    "wear_ring_length": 8.0,       # 后续预留：口环轴向长度
    "rear_chamber_gap": 2.0,       # 后续预留：后腔轴向间隙
}

import argparse
import json
from math import cos, isfinite, pi, radians, sin, tan
import os
from pathlib import Path
import tempfile

from main import connect


DOMAIN_LABELS = {"inlet": "InletFluid", "impeller": "ImpellerFluid",
                 "volute": "VoluteFluid"}


def smooth_growth(t):
    return t*t*t*(10+t*(-15+6*t))


def diffuser_gain(distance, p):
    """单侧增宽：入口增宽速度和加速度为零，过渡后为 tan(总角/2)。"""
    length=p["diffuser_blend_length"]; rate=tan(radians(p["diffuser_angle_deg"])/2)
    if distance >= length:
        return rate*(distance-length/2)
    t=distance/length
    return rate*length*(2.5*t**4-3*t**5+t**6)


def stationary_layout(p):
    rb=p["volute_base_diameter"]/2; ri=p["impeller_diameter"]/2
    z=p["front_shroud_thickness"]+p["blade_height"]/2
    end=radians(p["tongue_position_deg"]+p["volute_wrap_deg"])
    center_r=rb+p["volute_end_width"]/2
    scale=1+2*diffuser_gain(p["diffuser_length"],p)/p["volute_end_width"]
    bound=rb+p["volute_end_width"]+p["diffuser_length"]+p["volute_end_width"]*scale
    return {"ri":ri,"rb":rb,"z":z,"ro":bound,
            "zlo":z-p["volute_axial_width"]*scale/2,"zhi":z+p["volute_axial_width"]*scale/2,
            "end_theta":end,"end_center":[center_r*cos(end),center_r*sin(end),z],
            "direction":[-sin(end),cos(end),0],"diffuser_end_scale":scale}


def volute_profiles(p):
    """四角一致的梯形截面；中心轨迹半径 rb+w/2 单调增大形成螺旋。"""
    g=stationary_layout(p); sections=[]
    def section(theta,width,scale=1.0,advance=0.0):
        er=[cos(theta),sin(theta)];et=[-sin(theta),cos(theta)]
        center_r=g["rb"]+width/2
        center=[center_r*er[0]+advance*et[0],center_r*er[1]+advance*et[1],g["z"]]
        w=width*scale; hi=p["blade_height"]*scale; ho=p["volute_axial_width"]*scale
        corners=[]
        for u,v in ((-w/2,-hi/2),(w/2,-ho/2),(w/2,ho/2),(-w/2,hi/2)):
            corners.append([center[0]+u*er[0],center[1]+u*er[1],center[2]+v])
        return {"points":corners,"theta":theta,"width":w,"inner_height":hi,"outer_height":ho,
                "area":w*(hi+ho)/2,"center":center,"advance":advance}
    for i in range(p["volute_section_count"]):
        t=i/(p["volute_section_count"]-1)
        theta=radians(p["tongue_position_deg"]+p["volute_wrap_deg"]*t)
        width=p["volute_start_width"]+(p["volute_end_width"]-p["volute_start_width"])*smooth_growth(t)
        sections.append(section(theta,width))
    for i in range(1,p["diffuser_section_count"]):
        advance=p["diffuser_length"]*i/(p["diffuser_section_count"]-1)
        scale=1+2*diffuser_gain(advance,p)/p["volute_end_width"]
        sections.append(section(g["end_theta"],p["volute_end_width"],scale,advance))
    return sections


def blade_layout(p):
    """返回每片环形扇区的中心角、半角宽、内半径、外半径（rad/mm）。"""
    rh, ro = p["hub_diameter"]/2, p["impeller_diameter"]/2
    half = p["blade_thickness"]/(2*ro)
    beta = radians(p["blade_installation_angle_deg"])
    return [(beta+2*pi*i/p["blade_count"], half, rh, ro) for i in range(p["blade_count"])]


def validate(p):
    signed = {"blade_installation_angle_deg", "section_angle_deg", "tongue_position_deg", "tongue_orientation_deg"}
    for key, value in p.items():
        if key == "volute_section_shape":
            if value != "trapezoid":
                raise ValueError("本版蜗壳截面必须为 trapezoid")
            continue
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
    if not (rs < re < ro and rs < rh < ro):
        raise ValueError("轴径须小于 eye 与叶片根部定位圆直径，二者均须小于叶轮外径")
    half = p["blade_thickness"]/(2*ro)
    if not 0 < half < pi/p["blade_count"]:
        raise ValueError("扇区角宽须小于阵列节距，叶片之间必须有通道")
    if rh*sin(half) <= c or half+2*c/rh >= pi/p["blade_count"]:
        raise ValueError("扇区根部宽度或通道间隙不足以放置验证探针，请减小 probe_clearance")
    if c >= min(rs, p["blade_height"]/2, (rh-rs)/2, (ro-rh)/2, (ro-re)/2, (re-rs)/4):
        raise ValueError("probe_clearance 过大，无法放置内部验证点")
    if p["seal_chamber_diameter"] <= p["shaft_diameter"]:
        raise ValueError("预留密封腔径须大于轴径")
    g=stationary_layout(p)
    for key,minimum in (("volute_section_count",25),("diffuser_section_count",5),("volute_check_grid_count",21)):
        if type(p[key]) is not int or p[key] < minimum:
            raise ValueError(key+" 必须是至少 "+str(minimum)+" 的整数")
    if not (g["rb"]>g["ri"] and p["volute_connection_overlap"]<min(g["rb"]-g["ri"],p["volute_start_width"])):
        raise ValueError("基圆须在叶轮外，并为均匀间隙保留正宽度")
    if not p["volute_end_width"]>p["volute_start_width"] or p["volute_axial_width"]<=p["blade_height"]:
        raise ValueError("蜗壳必须渐扩且为非矩形梯形（外侧高度大于内侧）")
    if not 270 <= p["volute_wrap_deg"] <= 345 or not 0 < p["diffuser_angle_deg"] <= 10:
        raise ValueError("包角范围 270..345°，总扩散角范围 0..10°")
    if p["diffuser_blend_length"]>=p["diffuser_length"]:
        raise ValueError("扩散入口平滑段应短于总扩散段")
    if abs(p["tongue_orientation_deg"])>10 or g["rb"]+p["tongue_tail_length"]<g["ro"]:
        raise ValueError("隔舌角限 ±10°，刀具尾长须延伸到蜗壳包络外")
    if p["tongue_nose_radius"]>=p["volute_start_width"]/2:
        raise ValueError("隔舌鼻半径应小于最小截面宽度的一半")
    # 扩散段不得穿过隔舌长轴，否则可能被截断为孤立体；在角落和中心路径上检查。
    tongue=radians(p["tongue_position_deg"]+p["tongue_orientation_deg"])
    nose_center=(g["rb"]+p["tongue_nose_radius"])
    start=radians(p["tongue_position_deg"])
    # 起始端盖必须完全落在隔舌刀具内，避免留下未修圆的端盖尖角。
    for x,y,_ in volute_profiles(p)[0]["points"]:
        dx=x-nose_center*cos(start);dy=y-nose_center*sin(start)
        u=dx*cos(tongue)+dy*sin(tongue);v=-dx*sin(tongue)+dy*cos(tongue)
        radius=p["tongue_nose_radius"]
        inside=(u*u+v*v<=radius*radius*(1+1e-12)) if u<0 else (abs(v)<=radius*(1+1e-12) and u<=p["tongue_tail_length"])
        if not inside:
            raise ValueError("起始截面未被圆弧隔舌覆盖，请减小隔舌偏角或增大鼻半径")
    for section in volute_profiles(p)[p["volute_section_count"]:]:
        for x,y,_ in section["points"]:
            start=radians(p["tongue_position_deg"])
            lateral=-(x-nose_center*cos(start))*sin(tongue)+(y-nose_center*sin(start))*cos(tongue)
            if lateral >= -p["tongue_nose_radius"]-p["alignment_tolerance"]:
                raise ValueError("扩散段触及隔舌/起始平面，请缩短扩散段或减小包角/扩散角")
            if (x*x+y*y)**0.5 < g["rb"]-p["alignment_tolerance"]:
                raise ValueError("扩散段侵入基圆，请增加入口平滑长度或减小扩散角")
    return p


def expected_volumes(p):
    validate(p)
    ro, rs, rh = p["impeller_diameter"]/2, p["shaft_diameter"]/2, p["hub_diameter"]/2
    angle = p["blade_thickness"]/ro
    blade_area = angle*(ro*ro-rh*rh)/2
    return {
        **stationary_volumes(p),
        "inlet": pi*((p["eye_diameter"]/2)**2-rs*rs)*p["inlet_length"],
        "impeller": (pi*(ro*ro-rs*rs)-p["blade_count"]*blade_area)*p["blade_height"]
            + pi*((p["eye_diameter"]/2)**2-rs*rs)*p["front_shroud_thickness"],
    }


def stationary_volumes(p):
    # 放样、圆柱修剪及隔舌布尔后没有简单解析并集体积，不伪造精确值。
    return {"volute":None}


def blade_sketch(p, index):
    """两条精确同心圆弧和两条径向线闭合；只沿轴向拉伸，不拟合圆弧。"""
    from ansys.geometry.core.math import Plane, Point2D, Point3D
    from ansys.geometry.core.misc import UNITS
    from ansys.geometry.core.sketch import Sketch
    theta, half, rh, ro = blade_layout(p)[index]
    sk = Sketch(Plane(origin=Point3D([0, 0, p["front_shroud_thickness"]], unit=UNITS.mm),
                      direction_x=[cos(theta), sin(theta), 0],
                      direction_y=[-sin(theta), cos(theta), 0]))
    def point(r, t):
        return Point2D([r*cos(t), r*sin(t)], unit=UNITS.mm)
    il, ol, ou, iu = point(rh,-half), point(ro,-half), point(ro,half), point(rh,half)
    center = Point2D([0,0], unit=UNITS.mm)
    sk.segment(il, ol)
    sk.arc(ol, ou, center, clockwise=False)
    sk.segment(ou, iu)
    sk.arc(iu, il, center, clockwise=True)
    return sk


def build(modeler, p):
    """整体流体并入环腔、桥接和切向管，再作圆柱面及平面精确分区。"""
    from ansys.geometry.core.math import Plane, Point2D, Point3D
    from ansys.geometry.core.misc import UNITS
    from ansys.geometry.core.sketch import Sketch
    validate(p)
    design = modeler.create_design("pump_simple")

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
    design, fluid = add_stationary_flow(modeler, design, fluid, p)
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
    expected_core=targets["inlet"]+targets["impeller"]
    allowed=max(p["volume_absolute_tolerance"],whole_volume*p["volume_relative_tolerance"])
    if not isfinite(whole_volume) or whole_volume<=expected_core:
        raise RuntimeError("整体流体体积未包含有效蜗壳")
    # 先按 R=叶轮外径和切向 Y 接口分离静止域，再对核心整体执行原 Z=0 分区。
    stationary = partition_stationary(design, fluid, p)
    inlet, impeller = partition_at_interface(design, fluid, p)
    domains = {"inlet": inlet, "impeller": impeller, **stationary}
    total = sum(body.volume.to("mm^3").magnitude for body in domains.values())
    if abs(total-whole_volume) > allowed:
        raise RuntimeError("Z=0 分区前后体积不守恒")
    bodies = design.get_all_bodies()
    if len(bodies) != 3 or {body.name for body in bodies} != set(DOMAIN_LABELS.values()):
        raise RuntimeError("分区后应只有进口、叶轮、蜗壳三个命名实体")
    print(f"整体流体 {whole_volume:.6f} mm³；三域分区后体积守恒。")
    return design, domains


def spiral_native_source(p):
    sections=volute_profiles(p)
    return r'''# Python Script, API Version = V241
from System.Collections.Generic import List, ICollection
import math
stage = "read_volute_tools"
root = GetRootPart()
core = next(b for b in root.Bodies if b.Name == "InletFluid")
trim = next(b for b in root.Bodies if b.Name == "VoluteInnerTrim")
tongue = next(b for b in root.Bodies if b.Name == "TongueTool")
KernelBody = type(core.Shape)
profiles = List[ICollection[ITrimmedCurve]]()
planes = []
section_data = %r
for data in section_data:
    xyz, theta = data
    pts = [Point.Create(MM(x),MM(y),MM(z)) for x,y,z in xyz]
    curves = List[ITrimmedCurve]()
    for i in range(4):
        curves.Add(CurveSegment.Create(pts[i],pts[(i+1) %% 4]))
    profiles.Add(curves)
    planes.append(Plane.Create(Frame.Create(pts[0],Direction.Create(math.cos(theta),math.sin(theta),0),Direction.DirZ)))
stage = "smooth_spiral_and_diffuser_single_loft"
volute = KernelBody.LoftProfiles(profiles,False,False)
if not volute.IsClosed:
    caps = List[KernelBody]()
    caps.Add(KernelBody.CreatePlanarBody(planes[0],profiles[0]))
    caps.Add(KernelBody.CreatePlanarBody(planes[-1],profiles[profiles.Count-1]))
    volute.Fuse(caps,False,None)
if not volute.IsClosed or volute.PieceCount != 1 or volute.Volume <= 0:
    raise Exception("Spiral/diffuser loft is not a connected closed solid")
stage = "enforce_exact_base_cylinder"
tools = List[KernelBody]()
tools.Add(trim.Shape.Copy())
volute.Subtract(tools)
stage = "fuse_uniform_gap_and_volute"
combined = core.Shape.Copy()
parts = List[KernelBody]()
parts.Add(volute)
combined.Fuse(parts,False,None)
stage = "round_cutwater_nose"
tools = List[KernelBody]()
tools.Add(tongue.Shape.Copy())
combined.Subtract(tools)
if not combined.IsClosed or combined.PieceCount != 1 or combined.Volume <= core.Shape.Volume:
    raise Exception("Tongue or loft created an open/disconnected/empty fluid")
stage = "commit_volute"
core.Delete()
trim.Delete()
tongue.Delete()
DesignBody.Create(root,"InletFluid",combined)
result = {"closed":"True","pieces":"1","profile_count":str(len(section_data)),"volume_m3":str(combined.Volume)}
''' % ([(v["points"],v["theta"]) for v in sections],)


def add_stationary_flow(modeler,design,whole,p):
    """仅蜗壳调用已获授权的 V241 原生平滑放样；其余工具保持公共基础 API。"""
    from ansys.geometry.core.math import Plane,Point2D,Point3D
    from ansys.geometry.core.misc import UNITS
    from ansys.geometry.core.sketch import Sketch
    from spaceclaim_diagnostics import run_checked_script
    g=stationary_layout(p);mm=UNITS.mm;ov=p["volute_connection_overlap"]
    def cylinder(label, radius, z, height):
        sk=Sketch(Plane(origin=Point3D([0,0,z],unit=mm))).circle(Point2D([0,0],unit=mm),radius*mm)
        return design.extrude_sketch(label,sk,distance=height*mm)
    bridge=cylinder("UniformGap",g["rb"]+ov,p["front_shroud_thickness"],p["blade_height"])
    bridge.subtract(cylinder("UniformGapBore",g["ri"]-ov,p["front_shroud_thickness"],p["blade_height"]))
    whole.unite(bridge)
    extra=p["blade_tool_overrun"];zlo=g["zlo"]-extra;height=g["zhi"]-g["zlo"]+2*extra
    cylinder("VoluteInnerTrim",g["rb"],zlo,height)
    angle=radians(p["tongue_position_deg"]); orientation=angle+radians(p["tongue_orientation_deg"])
    radius=p["tongue_nose_radius"];length=p["tongue_tail_length"]
    center=(g["rb"]+radius)
    plane=Plane(origin=Point3D([center*cos(angle),center*sin(angle),zlo],unit=mm),
                direction_x=[cos(orientation),sin(orientation),0],direction_y=[-sin(orientation),cos(orientation),0])
    sk=Sketch(plane)
    point=lambda x,y:Point2D([x,y],unit=mm)
    a,b,c,d=point(0,-radius),point(length,-radius),point(length,radius),point(0,radius)
    sk.segment(a,b).arc(b,c,point(length,0)).segment(c,d).arc(d,a,point(0,0))
    design.extrude_sketch("TongueTool",sk,distance=height*mm)
    values,design=run_checked_script(modeler,spiral_native_source(p),"variable_volute")
    if values.get("closed")!="True" or values.get("pieces")!="1":
        raise RuntimeError("V241 蜗壳原生放样/隔舌检查失败")
    bodies=design.get_all_bodies()
    if len(bodies)!=1:
        raise RuntimeError("蜗壳提交后不应残留工具体")
    return design,bodies[0]


def partition_stationary(design,whole,p):
    from ansys.geometry.core.math import Plane,Point2D,Point3D
    from ansys.geometry.core.misc import UNITS
    from ansys.geometry.core.sketch import Sketch
    g=stationary_layout(p);mm=UNITS.mm;m=p["blade_tool_overrun"]
    zlo=min(-p["inlet_length"],g["zlo"])-m;zhi=max(g["zhi"],p["front_shroud_thickness"]+p["blade_height"]+p["rear_shroud_thickness"])+m
    def mask():
        sk=Sketch(Plane(origin=Point3D([0,0,zlo],unit=mm))).circle(Point2D([0,0],unit=mm),g["ri"]*mm)
        return design.extrude_sketch("RotorCylinderMask",sk,distance=(zhi-zlo)*mm)
    volute=whole.copy(design,"VoluteFluid")
    whole.intersect(mask());volute.subtract(mask())
    return {"volute":volute}


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
    if not isfinite(actual) or actual <= 0 or (target is not None and abs(actual-target) > max(
            p["volume_absolute_tolerance"], target*p["volume_relative_tolerance"])):
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
    """探针遵循实际环形扇区：叶片中线为空，角边之外及根部内环为流体。"""
    rs = p["shaft_diameter"]/2
    zmid = p["front_shroud_thickness"]+p["blade_height"]/2
    c = p["probe_clearance"]
    points, expected = [[0,0,zmid]], [False]
    for theta, half, rh, ro in blade_layout(p):
        for radius in (rh+(ro-rh)*0.2, (rh+ro)/2, ro-2*c):
            for offset, fluid in ((0,False), (half+2*c/radius,True), (-half-2*c/radius,True)):
                points.append([radius*cos(theta+offset),radius*sin(theta+offset),zmid]); expected.append(fluid)
        angle = theta+pi/p["blade_count"]
        for radius in ((rs+rh)/2, (rh+ro)/2, ro-2*c):
            points.append([radius*cos(angle),radius*sin(angle),zmid]); expected.append(True)
    for i in range(p["blade_count"]*4):
        angle = 2*pi*i/(p["blade_count"]*4)
        radius = (rs+p["hub_diameter"]/2)/2
        points.append([radius*cos(angle),radius*sin(angle),zmid]); expected.append(True)
    probes = require_probes(mesh,points,expected,"扇区叶片/叶间通道检查")
    return {"passed":probes["passed"],"probes":probes,"probe_count":len(points),"blade_width_definition":"outer arc length"}


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


def fit_cylinder_section(mesh, radius_range, label, p, z_range=None, radius_hint=None):
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
    if radius_hint is not None:
        selected &= np.all(np.abs(radii-radius_hint) <= p["alignment_tolerance"], axis=1)
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
    """诊断探针：失败只返回完整详情并警告，不改变几何，也不抛出异常。"""
    import numpy as np
    result = {"kind": "probe_check", "label": label, "passed": False, "failures": []}
    try:
        points = np.asarray(points, dtype=float).reshape(-1, 3)
        expected = np.asarray(expected, dtype=bool)
        actual = mesh_inside(mesh, points)
        failed = np.flatnonzero(actual != expected)
        result.update(probe_count=len(points), failed_count=int(failed.size), passed=not bool(failed.size))
        result["failures"] = [{"index": int(i), "xyz_mm": points[i].tolist(),
                               "expected_fluid": bool(expected[i]), "actual_fluid": bool(actual[i])}
                              for i in failed]
    except Exception as error:
        result["error"] = str(error)
    if not result["passed"]:
        print(f"探针警告（非致命）：{label}；失败数 {result.get('failed_count', '未完成分类')}，详情写入检查报告。")
    return result


def run_probe_check(function, *args):
    """仅隔离探针检查阶段的分类/采样异常；几何和导出异常保持原处理。"""
    try:
        return function(*args)
    except Exception as error:
        warning = {"kind": "probe_check", "label": function.__name__, "passed": False,
                   "failures": [], "error": str(error)}
        print("探针警告（非致命）："+function.__name__+"："+str(error))
        return warning


def collect_probe_warnings(value, path=""):
    warnings = []
    if isinstance(value, dict):
        if value.get("kind") == "probe_check" and not value.get("passed", False):
            warnings.append({"path": path, **value})
        else:
            for key, item in value.items():
                warnings.extend(collect_probe_warnings(item, path+"/"+str(key)))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            warnings.extend(collect_probe_warnings(item, path+"/"+str(index)))
    return warnings


def sector_fluid_expected(x, y, p):
    """主流道解析掩码：圆盘减轴孔及 n 个实际环形扇区；叶片位置始终为 False。"""
    import numpy as np
    radius = np.hypot(x, y)
    expected = (radius > p["shaft_diameter"]/2) & (radius < p["impeller_diameter"]/2)
    angle = np.arctan2(y, x)
    for theta, half, rh, ro in blade_layout(p):
        delta = np.arctan2(np.sin(angle-theta), np.cos(angle-theta))
        epsilon = 32*np.finfo(float).eps
        blade = (radius >= rh-epsilon*ro) & (radius <= ro+epsilon*ro) & (np.abs(delta) <= half+epsilon)
        expected &= ~blade
    return expected


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
    inner = (rs+rh)/2
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
    probe_results = [require_probes(meshes["inlet"], points, expected, "进口环形流体/轴孔检查")]
    # 采样位置由半径与角度在全局 XY 中直接构造，两个域用完全相同的 XY。
    xy = []
    for i in range(p["blade_count"]):
        angle = radians(p["blade_installation_angle_deg"])+(2*i+1)*pi/p["blade_count"]
        for radius in (inner, (inner+re)/2, (inner+3*re)/4):
            xy.append([radius*cos(angle), radius*sin(angle)])
    dz = min(p["front_shroud_thickness"]/2, length/4)  # 离开接口三角面，不在 Z=±0.3 的近边界层采样
    for name, z in (("inlet", -dz), ("impeller", dz)):
        probes = [[x, y, z] for x, y in xy]
        probe_results.append(require_probes(meshes[name], probes, [True]*len(probes), name+" 接口通流检查"))
        probe_results.append(require_probes(meshes[name], [[0, 0, z]], [False], name+" 贯通轴孔检查"))
    return {"passed": all(item["passed"] for item in probe_results), "probe_results": probe_results, "inlet_probe_count": len(points), "throughflow_paths": len(xy),
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
    """中高平面按完整圆盘减轴/环形扇区槽的解析定义检查，不强求奇数叶片轴截面镜像对称。"""
    import numpy as np
    r, rs, a = p["impeller_diameter"]/2, p["shaft_diameter"]/2, p["blade_thickness"]/2
    u = np.linspace(-r, r, p["section_grid_count"])
    xx, yy = np.meshgrid(u, u)
    x, y = xx.ravel(), yy.ravel()
    radius = np.hypot(x, y)
    expected = sector_fluid_expected(x, y, p)
    margin = np.minimum(np.abs(radius-r), np.abs(radius-rs))
    for theta, half, rh, ro in blade_layout(p):
        delta = np.arctan2(np.sin(np.arctan2(y,x)-theta),np.cos(np.arctan2(y,x)-theta))
        # 保守排除圆弧和两条角边附近的点，避免三角化近似引起误判。
        margin = np.minimum(margin, np.abs(radius-rh))
        margin = np.minimum(margin, radius*np.abs(np.sin(np.abs(delta)-half)))
    points = np.c_[x, y, np.full(x.size, p["front_shroud_thickness"]+p["blade_height"]/2)]
    # 同时避开解析边界和三角化边界，区别面片近似误差与大块流体缺失。
    _, closest = mesh.find_closest_cell(points, return_closest_point=True)
    safe = (margin > p["probe_clearance"]) & (np.linalg.norm(points-closest, axis=1) > p["probe_clearance"])
    if not safe.any():
        warning = {"kind": "probe_check", "label": "check_impeller_footprint", "passed": False,
                   "probe_count": 0, "failures": [], "error": "无远离边界的可用探针，跳过平面抽查"}
        print("探针警告（非致命）："+warning["error"])
        return warning
    probes = require_probes(mesh, points[safe], expected[safe], "叶轮完整圆盘/轴孔/阵列槽平面检查")
    return {"passed": probes["passed"], "probes": probes, "probe_count": int(safe.sum()), "blade_count": p["blade_count"],
            "definition": "full circular disk minus shaft and equally spaced annular sector slots"}


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
    probes = require_probes(mesh, points, expected, "前后盖板/进口眼壁面检查")
    return {"passed": probes["passed"], "probes": probes, "probe_count": len(points), "front_eye_z_mm": [0, tf],
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


def check_outer_rim(mesh, p):
    rh, ro = p["hub_diameter"]/2, p["impeller_diameter"]/2
    fit = fit_cylinder_section(mesh, (ro-(ro-rh)/4, ro+(ro-rh)/4), "叶轮外圆", p, radius_hint=ro)
    if abs(fit["radius_mm"]-ro) > p["alignment_tolerance"]:
        raise RuntimeError("叶轮外圆半径与扇区外弧不匹配")
    return {"passed": True, "outer_cylinder_fit": fit,
            "blade_tip_arc_radius_mm": ro, "blade_root_arc_radius_mm": rh}


def save_outer_rim_section(mesh, p, output):
    """放大实际 CAD 中高截面，并叠加精确外圆/扇区角边作为核对参考。"""
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    theta, half, rh, ro = blade_layout(p)[0]
    z = p["front_shroud_thickness"]+p["blade_height"]/2
    cut = mesh.slice(normal=[0,0,1],origin=[0,0,z])
    cut.save(output/"impeller_outer_rim_section.vtp")
    fig, ax = plt.subplots(figsize=(7,5),constrained_layout=True)
    cursor=0
    while cursor < len(cut.lines):
        k=cut.lines[cursor];ids=cut.lines[cursor+1:cursor+1+k];xyz=cut.points[ids]
        ax.plot(xyz[:,0],xyz[:,1],color="black",linewidth=1)
        cursor+=k+1
    t=np.linspace(theta-4*half,theta+4*half,301)
    ax.plot(ro*np.cos(t),ro*np.sin(t),"r--",label="exact OD reference")
    for a in (theta-half,theta+half):
        ax.plot([rh*cos(a),ro*cos(a)],[rh*sin(a),ro*sin(a)],"b--",linewidth=0.7)
    center=np.array([ro*cos(theta),ro*sin(theta)])
    extent=max(2*p["blade_thickness"],4*p["probe_clearance"])
    ax.set(xlim=(center[0]-extent,center[0]+extent),ylim=(center[1]-extent,center[1]+extent),
           xlabel="X (mm)",ylabel="Y (mm)",title="Actual outer-rim section vs analytic circle / sector sides")
    ax.set_aspect("equal");ax.legend()
    return save_png_safely(fig,output/"impeller_outer_rim_detail.png")


def check_stationary_paths(meshes,p):
    g=stationary_layout(p);points=[]
    for i in range(p["blade_count"]*3):
        t=2*pi*i/(p["blade_count"]*3)
        r=(g["ri"]+g["rb"])/2
        points.append([r*cos(t),r*sin(t),g["z"]])
    points.extend(v["center"] for v in volute_profiles(p)[3::4])
    probes=require_probes(meshes["volute"],points,[True]*len(points),"均匀间隙/渐扩蜗壳/切向扩散通流")
    return {"passed":probes["passed"],"probes":probes}


def check_stationary_geometry(meshes,p,body=None):
    import numpy as np
    mesh=meshes["volute"];profiles=volute_profiles(p);n=p["volute_check_grid_count"]
    errors=[];measurements=[]
    # 避开隔舌局部修剪区，核对主蜗壳及扩散段非矩形截面面积。
    indices=list(range(4,p["volute_section_count"]-2,max(1,p["volute_section_count"]//8)))
    indices+=list(range(p["volute_section_count"]+1,len(profiles),max(1,p["diffuser_section_count"]//4)))
    for index in indices:
        v=profiles[index];theta=v["theta"];w=v["width"];ho=v["outer_height"]
        uu,vv=np.meshgrid(np.linspace(-w/2,w/2,n),np.linspace(-ho/2,ho/2,n))
        points=np.c_[v["center"][0]+uu.ravel()*cos(theta),v["center"][1]+uu.ravel()*sin(theta),v["center"][2]+vv.ravel()]
        try:
            selected=mesh_inside(mesh,points)
            area=float(selected.mean()*w*ho)
            error=abs(area-v["area"])/v["area"]
            measurements.append({"profile":index,"measured_area_mm2":area,"design_area_mm2":v["area"],"relative_error":error})
            if error>p["volute_area_relative_tolerance"]:
                errors.append(index)
        except Exception as error:
            measurements.append({"profile":index,"error":str(error)})
            errors.append(index)
    small=[];short_edges=[]
    if body is not None:
        for i,face in enumerate(body.faces):
            area=float(face.area.to("mm^2").magnitude)
            if area<p["volute_min_face_area"]:
                small.append({"face_index":i,"area_mm2":area})
        for i,edge in enumerate(body.edges):
            length=float(edge.length.to("mm").magnitude)
            if length<p["volute_min_edge_length"]:
                short_edges.append({"edge_index":i,"length_mm":length})
    if small or short_edges:
        raise RuntimeError("蜗壳碎面/短边检查失败："+json.dumps({"small_faces":small,"short_edges":short_edges}))
    area_probes={"kind":"probe_check","label":"蜗壳/扩散截面面积抽查","passed":not bool(errors),
                 "failures":[item for item in measurements if item["profile"] in errors],"failed_profiles":errors}
    if errors:
        print("探针警告（非致命）：蜗壳截面面积抽查偏差，详情写入 JSON。")
    return {"passed":not bool(errors),"section_type":"trapezoid","measurements":measurements,"area_probes":area_probes,"small_faces":small,"short_edges":short_edges,
            "uniform_gap_mm":(p["volute_base_diameter"]-p["impeller_diameter"])/2,
            "note":"one smooth loft across spiral/diffuser; capsule tongue uses exact tangent arcs, CAD acceptance still required"}


def save_stationary_sections(mesh,p,output,name):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    g=stationary_layout(p);profiles=volute_profiles(p)
    chosen=[profiles[i] for i in (5,p["volute_section_count"]//2,p["volute_section_count"]-1,len(profiles)-1)]
    fig,axes=plt.subplots(1,5,figsize=(19,5),constrained_layout=True)
    xy=mesh.slice(normal=[0,0,1],origin=[0,0,g["z"]]);xy.save(output/"volute_mid_section.vtp")
    def lines(ax,cut,project):
        cursor=0
        while cursor<len(cut.lines):
            k=cut.lines[cursor];ids=cut.lines[cursor+1:cursor+1+k];xyz=cut.points[ids]
            a,b=project(xyz);ax.plot(a,b,color="black",linewidth=0.7);cursor+=k+1
    lines(axes[0],xy,lambda xyz:(xyz[:,0],xyz[:,1]));axes[0].set_title("Actual XY spiral + diffuser")
    for i,(ax,v) in enumerate(zip(axes[1:],chosen)):
        theta=v["theta"];center=np.array(v["center"])
        cut=mesh.slice(normal=[-sin(theta),cos(theta),0],origin=center)
        cut.save(output/("volute_profile_%d.vtp"%i))
        lines(ax,cut,lambda xyz:((xyz[:,0]-center[0])*cos(theta)+(xyz[:,1]-center[1])*sin(theta),xyz[:,2]-center[2]))
        ax.set(xlim=(-v["width"],v["width"]),ylim=(-p["volute_axial_width"],p["volute_axial_width"]),title="Section / area %.1f"%v["area"])
    for ax in axes:ax.set_aspect("equal")
    return save_png_safely(fig,output/"volute_sections.png")


def export_parts(modeler, design, p, output):
    """从同一个联合模型快照分离导出；不重新建模，不跨设计复制，不依赖新版抑制 API。"""
    from ansys.geometry.core.designer.design import DesignFileFormat
    assembly = output/"pump_assembly.step"
    assembly.unlink(missing_ok=True)
    design.download(str(assembly), format=DesignFileFormat.STEP)
    if not assembly.is_file() or not assembly.stat().st_size:
        raise RuntimeError("装配 STEP 导出失败")
    print("已导出三流体域装配（不合并）："+str(assembly))
    with tempfile.TemporaryDirectory(prefix="pump_export_") as folder:
        snapshot = Path(folder)/"stage1.scdocx"
        design.download(str(snapshot), format=DesignFileFormat.SCDOCX)
        if not snapshot.is_file() or not snapshot.stat().st_size:
            raise RuntimeError("联合模型快照保存失败")
        for name, label in DOMAIN_LABELS.items():
            # open_file 会关闭当前设计；每次均从同一未删减的联合快照恢复。
            isolated = modeler.open_file(str(snapshot))
            bodies = isolated.get_all_bodies()
            if len(bodies) != 3 or {b.name for b in bodies} != set(DOMAIN_LABELS.values()):
                raise RuntimeError("联合快照重载后未保留三个命名域")
            for body in bodies:
                if body.name != label:
                    isolated.delete_body(body)
            remaining = isolated.get_all_bodies()
            if len(remaining) != 1 or remaining[0].name != label:
                raise RuntimeError(name+" 导出隔离失败")
            # 这里只做必要的导出隔离；所有几何/网格/探针检查均放在四个 STEP 写出后。
            destination = output/(name+".step")
            destination.unlink(missing_ok=True)
            isolated.download(str(destination), format=DesignFileFormat.STEP)
            if not destination.is_file() or not destination.stat().st_size:
                raise RuntimeError("STEP 导出失败："+str(destination))
            print("已独立导出："+str(destination))
        # open_file 关闭过原设计；恢复联合快照后才进行后续检查。
        restored = modeler.open_file(str(snapshot))
        lookup = {body.name: body for body in restored.get_all_bodies()}
        return restored, {name:lookup[label] for name,label in DOMAIN_LABELS.items()}


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
        if (output/"outlet.step").exists():
            print("提示：旧 outlet.step 本次不更新，不属于当前三域装配；本阶段没有外接直管。")
        report_path = output/"pump_simple_checks.json"
        report_path.unlink(missing_ok=True)
        design, bodies = build(modeler, p)
        design, bodies = export_parts(modeler, design, p, output)
        report, meshes = {"params": p, "status": "exported_checks_pending",
                          "step_exports": ["pump_assembly.step", *[name+".step" for name in DOMAIN_LABELS]],
                          "construction": "V241 smooth variable trapezoid spiral and diffuser; rounded tongue; three exact Boolean partitions",
                          "volute_design":{"section_type":p["volute_section_shape"],"profile_count":len(volute_profiles(p)),
                                           "start_area_mm2":volute_profiles(p)[0]["area"],
                                           "spiral_end_area_mm2":volute_profiles(p)[p["volute_section_count"]-1]["area"],
                                           "diffuser_end_area_mm2":volute_profiles(p)[-1]["area"],
                                           "native_api":"V241 LoftProfiles(periodic=False, ruled=False)"},
                          "parts": {}}, {}
        try:
            for name, body in bodies.items():
                mesh, checks = check_part(body, name, p)
                meshes[name] = mesh
                if name == "impeller":
                    checks["channels"] = run_probe_check(check_channels, mesh, p)
                    checks["footprint"] = run_probe_check(check_impeller_footprint, mesh, p)
                    checks["shrouds"] = run_probe_check(check_shrouds, mesh, p)
                    checks["outer_rim"] = check_outer_rim(mesh, p)
                    checks["outer_rim_plot"] = optional_plot(save_outer_rim_section, mesh, p, output)
                checks["sections"] = optional_plot(save_sections if name in ("inlet","impeller") else save_stationary_sections, mesh, p, output, name)
                report["parts"][name] = checks
                print(f"{name}: CFD 流体实体，体积 {checks['volume_mm3']:.6f} mm³")
            core_meshes={name:meshes[name] for name in ("inlet","impeller")}
            report["interface"] = check_alignment(core_meshes, p)
            report["interface"]["cap_match"] = check_interface_caps(core_meshes, p)
            report["stationary_geometry"] = check_stationary_geometry(meshes, p, bodies["volute"])
            report["stationary_probes"] = run_probe_check(check_stationary_paths, meshes, p)
            report["alignment_plot"] = optional_plot(save_alignment_section, core_meshes, p, output)
            report["probe_warnings"] = collect_probe_warnings(report)
            report["status"] = "exported_with_probe_warnings" if report["probe_warnings"] else "exported_checks_passed"
            if report["probe_warnings"]:
                print(f"四个 STEP 导出成功；{len(report['probe_warnings'])} 项探针警告已记录，正常结束（退出码 0）。")
        except Exception as error:
            report["status"] = "exported_checks_failed"
            report["error"] = str(error)
            print("四个 STEP 已导出并保留；后续检查失败，详情见 pump_simple_checks.json。")
            raise
        finally:
            report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
        print("STEP 主流程完成；检查和警告详情见 pump_simple_checks.json。")
    finally:
        modeler.close()


if __name__ == "__main__":
    main()
