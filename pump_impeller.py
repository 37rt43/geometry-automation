"""五个后弯空间叶片扣除后的叶轮流体域；长度 mm，角度 deg。

python pump_impeller.py
python pump_impeller.py --host HOST --port 50051
python pump_impeller.py --preview-only

复用 pump_meridional 的 connect/build_and_check，不建立其他泵腔。
24.1 的放样、复制旋转和连通性检查通过 run_script_file 调用原生内核，
不使用要求 24R2 的 PyAnsys loft/rotate/get_collision 接口。
预览是根据截面参数计算的 Z=5 mm 剖面示意，不冒充已运行的 CAD 截图。
"""

import argparse
from dataclasses import dataclass
from math import cos, hypot, isfinite, log, radians, sin, sqrt, tan
import os
from pathlib import Path
from tempfile import TemporaryDirectory

from pump_meridional import Parameters as BlankParameters
from pump_meridional import build_and_check as build_blank, connect


@dataclass(frozen=True)
class Parameters:
    blade_count: int = 5
    inlet_radius: float = 27.5
    outlet_radius: float = 50.0
    wrap_degrees: float = 130.0
    tangential_thickness: float = 3.0
    inlet_twist_degrees: float = 8.0
    section_count: int = 31
    axial_overcut: float = 0.75
    preview_z: float = 5.0
    volume_tolerance_mm3: float = 1e-3

    def validate(self, blank):
        blank.validate()
        if not isinstance(self.blade_count, int) or self.blade_count < 2:
            raise ValueError("blade_count 必须为至少 2 的整数")
        if not isinstance(self.section_count, int) or self.section_count < 3:
            raise ValueError("section_count 必须为至少 3 的整数")
        values = (self.inlet_radius, self.outlet_radius, self.wrap_degrees,
                  self.tangential_thickness, self.axial_overcut, self.volume_tolerance_mm3)
        if any(not isfinite(v) or v <= 0 for v in values):
            raise ValueError("半径、包角、厚度、超切量和体积容差须为有限正数")
        if not isfinite(self.inlet_twist_degrees) or abs(self.inlet_twist_degrees) >= 45:
            raise ValueError("进口扭转偏移须为有限数且绝对值小于 45 度")
        if not (blank.hub_center[0] < self.inlet_radius < self.outlet_radius
                and self.inlet_radius >= blank.front_inlet_radius
                and self.outlet_radius == blank.outlet_radius):
            raise ValueError("叶片须从轮毂弯段外开始，出口半径须等于毛坯出口半径")
        if not blank.hub_outlet_z < self.preview_z < blank.front_outlet_z:
            raise ValueError("preview_z 须位于径向出口流道内部")
        return self


def front_z(radius, blank):
    """前盖板边界：精确四分之一圆弧，随后为水平出口边界。"""
    cr, cz = blank.front_center
    if radius < blank.front_inlet_radius - 1e-9:
        raise ValueError("查询位置不在前盖板径向段")
    if radius >= cr:
        return blank.front_outlet_z
    return cz - sqrt(max(0.0, blank.front_bend_radius**2 - (radius-cr)**2))


def section_data(radius, p, blank):
    fraction = (radius-p.inlet_radius)/(p.outlet_radius-p.inlet_radius)
    theta = -radians(p.wrap_degrees)*log(radius/p.inlet_radius)/log(p.outlet_radius/p.inlet_radius)
    twist = radians(p.inlet_twist_degrees)*(1-3*fraction**2+2*fraction**3)
    tangent_top = radius*tan(twist)
    # 顶部中心的真实圆柱半径略大于所在径向法平面的法向距离。
    physical_top_z = front_z(hypot(radius, tangent_top), blank)
    bottom_z = blank.hub_outlet_z  # 默认 r>=27.5>轮毂弧终点25，轮毂边界为Z=0。
    slope = tangent_top/(physical_top_z-bottom_z)
    return theta, bottom_z, physical_top_z, slope


def point_at(radius, z, side, p, blank):
    theta, bottom, top, slope = section_data(radius, p, blank)
    tangent = slope*(z-bottom) + side*p.tangential_thickness/2
    return (radius*cos(theta)-tangent*sin(theta),
            radius*sin(theta)+tangent*cos(theta), z)


def blade_sections(p, blank):
    p.validate(blank)
    sections, angles = [], []
    for i in range(p.section_count):
        radius = p.inlet_radius+(p.outlet_radius-p.inlet_radius)*i/(p.section_count-1)
        theta, bottom, top, slope = section_data(radius, p, blank)
        # 用该法平面内最小可能圆柱半径 r 的前盖板高度作保守上界。
        # 四边形上下均超切，扭转斜率外推；流体边界处保持指定偏移。
        low, high = bottom-p.axial_overcut, front_z(radius, blank)+p.axial_overcut
        quad = [point_at(radius, low, -1, p, blank), point_at(radius, low, 1, p, blank),
                point_at(radius, high, 1, p, blank), point_at(radius, high, -1, p, blank)]
        sections.append(quad)
        angles.append(theta)
    return sections, angles


def native_script(p, blank):
    sections, angles = blade_sections(p, blank)
    return '''# Python Script, API Version = V241
from System.Collections.Generic import List, ICollection
from SpaceClaim.Api.V241.Modeler import Body as KernelBody
sections = %r
angles = %r
blade_count = %r
eye_radius = %r
eye_z = %r
eye_height = %r
volume_tolerance = %r * 1e-9
root = GetRootPart()
if root.Bodies.Count != 1 or root.Bodies[0].Name != "ImpellerFluidBlank":
    raise Exception("Expected exactly one ImpellerFluidBlank")
blank = root.Bodies[0]
before = blank.Shape.Volume
profiles = List[ICollection[ITrimmedCurve]]()
planes = []
for corners, theta in zip(sections, angles):
    import math
    curves = List[ITrimmedCurve]()
    points = [Point.Create(MM(x), MM(y), MM(z)) for x,y,z in corners]
    for i in range(4):
        curves.Add(CurveSegment.Create(points[i], points[(i+1) %% 4]))
    profiles.Add(curves)
    planes.append(Plane.Create(Frame.Create(points[0],
                  Direction.Create(-math.sin(theta), math.cos(theta), 0), Direction.DirZ)))
# periodic=False, ruled=False: smooth loft, with identical corner correspondence.
first = KernelBody.LoftProfiles(profiles, False, False)
if not first.IsClosed:
    caps = List[KernelBody]()
    caps.Add(KernelBody.CreatePlanarBody(planes[0], profiles[0]))
    caps.Add(KernelBody.CreatePlanarBody(planes[-1], profiles[profiles.Count-1]))
    first.Fuse(caps, False, None)
if not first.IsClosed or first.PieceCount != 1 or first.Volume <= 0:
    raise Exception("Smooth loft did not produce one closed blade tool")
# Protect the inlet eye even against small spline overshoot between sections.
eye_plane = Plane.Create(Frame.Create(Point.Create(0, 0, MM(eye_z)), Direction.DirX, Direction.DirY))
eye_curves = List[ITrimmedCurve]()
eye_curves.Add(CurveSegment.Create(Circle.Create(eye_plane.Frame, MM(eye_radius))))
eye = KernelBody.ExtrudeProfile(eye_plane, eye_curves, MM(eye_height))
eye_tools = List[KernelBody]()
eye_tools.Add(eye)
first.Subtract(eye_tools)
tools = List[KernelBody]()
tools.Add(first)
for i in range(1, blade_count):
    copy = first.Copy()
    copy.Transform(Matrix.CreateRotation(Line.Create(Point.Origin, Direction.DirZ),
                                        DEG(i*360.0/blade_count)))
    tools.Add(copy)
for tool in tools:
    if not tool.IsClosed or tool.PieceCount != 1 or tool.Volume <= 0:
        raise Exception("Blade tool is not a connected solid")
# Check all pairs on transient copies, never consuming the actual tools.
for i in range(blade_count):
    for j in range(i+1, blade_count):
        probe = tools[i].Copy()
        others = List[KernelBody]()
        others.Add(tools[j].Copy())
        probe.Intersect(others)
        if probe.Faces.Count > 0 and probe.Volume > volume_tolerance:
            raise Exception("Blade tools intersect: %%s, %%s" %% (i,j))
blank.Shape.Subtract(tools)
after = blank.Shape.Volume
if not (0 < after < before-volume_tolerance):
    raise Exception("Boolean subtraction did not reduce fluid volume")
if not blank.Shape.IsClosed or blank.Shape.PieceCount != 1:
    raise Exception("Fluid must be one connected closed solid")
blank.SetName("ImpellerFluid")
# Tools/caps/eye are transient kernel bodies, never added to the document.
tools.Clear()
eye_tools.Clear()
if root.Bodies.Count != 1:
    raise Exception("Unexpected residual bodies")
result = {"pieces": "1", "pair_checks": str(blade_count*(blade_count-1)//2),
          "closed": "True", "before_mm3": str(before*1e9), "after_mm3": str(after*1e9)}
''' % (sections, angles, p.blade_count, p.inlet_radius,
       blank.hub_outlet_z-5*p.axial_overcut,
       blank.inlet_z-blank.hub_outlet_z+10*p.axial_overcut, p.volume_tolerance_mm3)


def save_preview(p, blank, destination):
    """参数计算剖面：蓝色流体、白色叶片扣除区域，可见五条叶道及连通眼。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle as PlotCircle, Polygon

    p.validate(blank)
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.add_patch(PlotCircle((0, 0), p.outlet_radius, color="#89cee8"))
    cr, cz = blank.hub_center
    inner = cr-sqrt(max(0, blank.hub_bend_radius**2-(p.preview_z-cz)**2))
    ax.add_patch(PlotCircle((0, 0), inner, color="white"))
    for i in range(p.blade_count):
        rotation = radians(i*360/p.blade_count)
        curves = []
        for side in (-1, 1):
            edge = []
            for j in range(501):
                r = p.inlet_radius+(p.outlet_radius-p.inlet_radius)*j/500
                x,y,z = point_at(r, p.preview_z, side, p, blank)
                edge.append((x*cos(rotation)-y*sin(rotation), x*sin(rotation)+y*cos(rotation)))
            curves.append(edge)
        ax.add_patch(Polygon(curves[0]+curves[1][::-1], closed=True,
                             facecolor="white", edgecolor="#354f60", linewidth=1))
    ax.add_patch(PlotCircle((0,0), p.inlet_radius, fill=False,
                            linestyle="--", edgecolor="#28707a", linewidth=1))
    ax.annotate("Connected inlet eye", (0, 22), ha="center", fontsize=10)
    ax.set(xlim=(-55,55), ylim=(-55,55), xlabel="X (mm)", ylabel="Y (mm)",
           title=f"ImpellerFluid: Z={p.preview_z:g} mm section\nParameter-derived preview (not a CAD screenshot)")
    ax.set_aspect("equal")
    ax.grid(alpha=0.15)
    fig.tight_layout()
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=180)
    plt.close(fig)


def build_impeller(modeler, p, blank):
    p.validate(blank)
    design = build_blank(modeler, blank)
    before = design.get_all_bodies()[0].volume.to("mm^3").magnitude
    with TemporaryDirectory(prefix="pump_impeller_") as directory:
        path = Path(directory)/"blades_v241.py"
        path.write_text(native_script(p, blank), encoding="ascii")
        values, design = modeler.run_script_file(path, import_design=True)
    if (values.get("pieces") != "1" or values.get("closed") != "True"
            or values.get("pair_checks") != str(p.blade_count*(p.blade_count-1)//2)):
        raise RuntimeError(f"原生叶片/连通性检查未通过：{values}")
    bodies = design.get_all_bodies()
    if len(bodies) != 1 or bodies[0].is_surface or bodies[0].name != "ImpellerFluid":
        raise RuntimeError("最终必须只有一个 ImpellerFluid 实体，无残留工具体")
    after = bodies[0].volume.to("mm^3").magnitude
    if not isfinite(after) or not 0 < after < before-p.volume_tolerance_mm3:
        raise RuntimeError("布尔减后的流体体积未减小")
    print(f"叶片两两检查通过；流体连通块=1；体积 {before:.6f} -> {after:.6f} mm³")
    return design


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default=os.getenv("ANSRV_GEO_HOST"))
    parser.add_argument("--port", type=int, default=int(os.getenv("ANSRV_GEO_PORT", "50051")))
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).resolve().with_name("pump_impeller.step"))
    parser.add_argument("--preview", type=Path,
                        default=Path(__file__).resolve().with_name("pump_impeller_section.png"))
    parser.add_argument("--preview-only", action="store_true")
    args = parser.parse_args()
    if args.output.suffix.lower() not in (".step", ".stp"):
        parser.error("--output 必须为 .step 或 .stp 文件")
    blank = BlankParameters()
    p = Parameters().validate(blank)
    save_preview(p, blank, args.preview.resolve())
    print(f"参数剖面预览：{args.preview.resolve()}")
    if args.preview_only:
        print("离线预览已生成；尚未进行 CAD 放样、布尔减或导出 STEP。")
        return
    from ansys.geometry.core.designer.design import DesignFileFormat
    modeler = connect(args.host, args.port)
    try:
        design = build_impeller(modeler, p, blank)
        destination = args.output.resolve()
        design.download(str(destination), format=DesignFileFormat.STEP)
        if not destination.is_file() or destination.stat().st_size == 0:
            raise RuntimeError("STEP 导出失败：文件不存在或为空")
        print(f"已导出：{destination}")
    finally:
        modeler.close()


if __name__ == "__main__":
    main()
