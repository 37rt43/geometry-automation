"""闭式离心泵无叶片流体回转毛坯：精确三次 Bezier 子午面，单位 mm。

python pump_meridional.py
python pump_meridional.py --host HOST --port 50051
python pump_meridional.py --preview-only

复用 main.connect 的 SpaceClaim 241 启动/连接方式。24.1 原生 NURBS
用四个单位权重控制点精确表达 Bezier，原生 RevolveFaces 回转；不调用
要求 24R2 的 PyAnsys revolve_sketch。不创建金属、叶片或其他泵腔。
build_and_check(modeler, Parameters()) 可供后续与新尺寸匹配的叶片脚本复用。
"""

import argparse
from dataclasses import dataclass, fields
from math import isclose, isfinite, pi
import os
from pathlib import Path

from main import connect
from spaceclaim_diagnostics import run_checked_script


@dataclass(frozen=True)
class Parameters:
    outlet_radius: float = 140.0
    hub_inlet_radius: float = 12.0
    front_inlet_radius: float = 45.0
    inlet_z: float = 75.0
    straight_end_z: float = 60.0
    hub_outlet_z: float = 0.0
    front_outlet_z: float = 8.5
    hub_control1_z: float = 15.0
    hub_control2_r: float = 70.0
    front_control1_z: float = 25.0
    front_control2_r: float = 90.0
    curve_tolerance_mm: float = 0.01
    preview_samples: int = 401

    @property
    def profile_kind(self):
        return "cubic_bezier"

    @property
    def outlet_width(self):
        return self.front_outlet_z-self.hub_outlet_z

    @property
    def hub_controls(self):
        return ((self.hub_inlet_radius, self.straight_end_z),
                (self.hub_inlet_radius, self.hub_control1_z),
                (self.hub_control2_r, self.hub_outlet_z),
                (self.outlet_radius, self.hub_outlet_z))

    @property
    def front_controls(self):
        return ((self.front_inlet_radius, self.straight_end_z),
                (self.front_inlet_radius, self.front_control1_z),
                (self.front_control2_r, self.front_outlet_z),
                (self.outlet_radius, self.front_outlet_z))

    def vertices(self):
        """A,B,C,F,E,D，闭合端点依次对应指定子午面。"""
        return ((self.hub_inlet_radius, self.inlet_z), self.hub_controls[0],
                self.hub_controls[-1], self.front_controls[-1],
                self.front_controls[0], (self.front_inlet_radius, self.inlet_z))

    def validate(self):
        if any(not isfinite(getattr(self, f.name)) for f in fields(self)):
            raise ValueError("所有尺寸和采样数必须有限")
        if type(self.preview_samples) is not int or self.preview_samples < 101:
            raise ValueError("preview_samples 必须为至少 101 的整数")
        if not (0 < self.hub_inlet_radius < self.front_inlet_radius < self.outlet_radius
                and self.inlet_z > self.straight_end_z > self.front_outlet_z > self.hub_outlet_z
                and 0 < self.curve_tolerance_mm <= 0.01):
            raise ValueError("进口半径、高度、正出口宽度或曲线误差限制不合法")
        for controls in (self.hub_controls, self.front_controls):
            if not (controls[0][0] == controls[1][0] < controls[2][0] < controls[3][0]
                    and controls[0][1] > controls[1][1] > controls[2][1] == controls[3][1]):
                raise ValueError("控制点必须保证轴向进口切向、径向出口切向及单调曲线")
        # Bezier 凸包/导数性质：两侧各自 r 严格递增、z 递减，无自交。
        # Q 的 r,z 控制分量均不小于 P；在同一 r 上 Q 比 P 更高，不会交叉。
        if any(q[0] < p[0] or q[1] < p[1]
               for p,q in zip(self.hub_controls, self.front_controls)):
            raise ValueError("两侧控制多边形不满足不交叉的充分条件")
        return self


def bezier(controls, u):
    if not 0 <= u <= 1:
        raise ValueError("Bezier 参数须在 [0,1] 内")
    weights = ((1-u)**3, 3*(1-u)**2*u, 3*(1-u)*u*u, u**3)
    return tuple(sum(w*p[k] for w,p in zip(weights, controls)) for k in (0,1))


def derivative(controls, u):
    return tuple(3*((1-u)**2*(controls[1][k]-controls[0][k])
                    + 2*(1-u)*u*(controls[2][k]-controls[1][k])
                    + u*u*(controls[3][k]-controls[2][k])) for k in (0,1))


def boundary_z(p, side, radius):
    """给定真实圆柱半径查盖板高度，供检查或后续叶片脚本复用。"""
    controls = p.hub_controls if side == "hub" else p.front_controls if side == "front" else None
    if controls is None:
        raise ValueError("side 必须是 hub 或 front")
    if radius < controls[0][0]-1e-9 or radius > controls[-1][0]+1e-9:
        raise ValueError("半径超出该曲线定义范围")
    if radius <= controls[0][0]:
        return controls[0][1]
    if radius >= controls[-1][0]:
        return controls[-1][1]
    low, high = 0.0, 1.0
    for _ in range(60):
        mid = (low+high)/2
        if bezier(controls, mid)[0] < radius:
            low = mid
        else:
            high = mid
    return bezier(controls, (low+high)/2)[1]


def profile_data(p):
    """绘图、解析体积和原生建模均取自这里的同一组控制点。"""
    p.validate()
    hub = [bezier(p.hub_controls, i/(p.preview_samples-1)) for i in range(p.preview_samples)]
    front = [bezier(p.front_controls, i/(p.preview_samples-1)) for i in range(p.preview_samples)]
    a,b,c,f,e,d = p.vertices()
    outline = [a]+hub+[f]+front[-2::-1]+[d,a]
    return {"hub_controls": p.hub_controls, "front_controls": p.front_controls,
            "hub": hub, "front": front, "outline": outline, "vertices": (a,b,c,f,e,d)}


def check_profile(p):
    p.validate()
    for controls in (p.hub_controls, p.front_controls):
        start, end = derivative(controls, 0), derivative(controls, 1)
        if not (start[0] == 0 and start[1] < 0 and end[0] > 0 and end[1] == 0):
            raise ValueError("进口或出口切向不满足指定方向")
    widths = [boundary_z(p,"front",r)-boundary_z(p,"hub",r)
              for r in (p.front_inlet_radius+(p.outlet_radius-p.front_inlet_radius)*i/1000
                        for i in range(1001))]
    if min(widths) <= 0:
        raise ValueError("径向流道宽度非正或两侧曲线相交")
    if any(b-a > 1e-6 for a,b in zip(widths, widths[1:])):
        raise ValueError("径向流道未向出口逐渐收窄")
    outline = profile_data(p)["outline"]
    if outline[0] != outline[-1] or min(r for r,z in outline) <= 0:
        raise ValueError("轮廓未闭合或误接到旋转轴")
    # validate 的单调性和控制多边形排序保证完整边界简单；此处另查同半径宽度。
    return {"minimum_width_mm": min(widths), "maximum_radial_width_mm": max(widths),
            "closed": True, "tangent_continuous": True}


def expected_volume(p):
    """精确多项式积分 V=pi*闭环积分(r² dz)，没有采样积分误差。"""
    from numpy.polynomial import Polynomial
    p.validate()

    def integral(controls):
        def coefficients(k):
            a,b,c,d = (point[k] for point in controls)
            return (a, 3*(b-a), 3*(a-2*b+c), -a+3*b-3*c+d)
        r, z = Polynomial(coefficients(0)), Polynomial(coefficients(1))
        antiderivative = (r*r*z.deriv()).integ()
        return antiderivative(1)-antiderivative(0)

    moment = (p.hub_inlet_radius**2*(p.straight_end_z-p.inlet_z)
              + integral(p.hub_controls) + p.outlet_radius**2*p.outlet_width
              - integral(p.front_controls)
              + p.front_inlet_radius**2*(p.inlet_z-p.straight_end_z))
    volume = pi*moment
    if not isfinite(volume) or volume <= 0:
        raise ValueError("流体回转体积必须为正")
    return volume


def save_preview(p, destination):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    data = profile_data(p)
    fig, ax = plt.subplots(figsize=(12,7))
    ax.fill(*zip(*data["outline"]), color="#a9ddee", alpha=0.85, label="Fluid section")
    for name,color in (("hub","#176684"),("front","#c66c27")):
        controls = data[name+"_controls"]
        ax.plot(*zip(*controls), "--o", color=color, alpha=0.65, label=name+" control polygon")
        ax.plot(*zip(*data[name]), color=color, linewidth=2.5, label=name+" cubic Bezier")
        for i,(r,z) in enumerate(controls):
            label = ("P" if name=="hub" else "Q")+str(i)+f" ({r:g},{z:g})"
            ax.annotate(label,(r,z),xytext=(7,12 if name=="front" else -17),
                        textcoords="offset points",fontsize=8)
    for label,(r,z) in zip(("A","B","C","F","E","D"),data["vertices"]):
        ax.plot(r,z,"o",color="#333333",markersize=4)
        ax.annotate(label+(f" ({r:g},{z:g})" if label in ("A","D") else ""),
                    (r,z),xytext=(-13,7),textcoords="offset points",fontweight="bold")
    ax.axvline(0,color="#b95353",linestyle="--",linewidth=1)
    ax.annotate("",(p.hub_inlet_radius,p.inlet_z+7),(p.front_inlet_radius,p.inlet_z+7),
                arrowprops={"arrowstyle":"<->"})
    ax.text((p.hub_inlet_radius+p.front_inlet_radius)/2,p.inlet_z+9,
            f"Rh={p.hub_inlet_radius:g}, Re={p.front_inlet_radius:g}",ha="center",fontsize=9)
    ax.annotate("",(p.outlet_radius+8,p.hub_outlet_z),(p.outlet_radius+8,p.front_outlet_z),
                arrowprops={"arrowstyle":"<->"})
    ax.text(p.outlet_radius+10,(p.hub_outlet_z+p.front_outlet_z)/2,
            f"b2={p.outlet_width:g}",va="center",fontsize=9)
    ax.text(62,70,f"R2={p.outlet_radius:g} / D2={2*p.outlet_radius:g}\n"
            f"Inlet Z={p.inlet_z:g}; straight length={p.inlet_z-p.straight_end_z:g} mm",fontsize=10)
    ax.set(xlabel="r = X (mm)",ylabel="Z (mm)",title="ImpellerFluidBlank: cubic Bezier meridional section",
           xlim=(-5,p.outlet_radius+27),ylim=(p.hub_outlet_z-10,p.inlet_z+17))
    ax.set_aspect("equal",adjustable="box")
    ax.grid(alpha=0.22)
    ax.legend(loc="upper right",fontsize=8)
    fig.tight_layout()
    destination.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(destination,dpi=180)
    plt.close(fig)


def save_3d_preview(p, destination, section=False):
    """由同一轮廓数据生成参数曲面预览；非后端截图。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    data = profile_data(p)
    theta = np.linspace(0,pi if section else 2*pi,101)
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111,projection="3d")
    faces,colors = [],[]

    def surface(rz,color):
        points = np.array(rz)[sorted(set(range(0,len(rz),4)) | {len(rz)-1})]
        r,z = points[:,0,None],points[:,1,None]
        xyz = np.stack(np.broadcast_arrays(r*np.cos(theta),r*np.sin(theta),z),axis=-1)
        for i in range(len(points)-1):
            for j in range(len(theta)-1):
                faces.append([xyz[i,j],xyz[i+1,j],xyz[i+1,j+1],xyz[i,j+1]])
                colors.append(color)
    surface(data["hub"],"#548a9d")
    surface(data["front"],"#98c9df")
    surface([(p.hub_inlet_radius,p.inlet_z),(p.hub_inlet_radius,p.straight_end_z)],"#467a91")
    surface([(p.front_inlet_radius,p.inlet_z),(p.front_inlet_radius,p.straight_end_z)],"#81b5d0")
    surface([(p.hub_inlet_radius,p.inlet_z),(p.front_inlet_radius,p.inlet_z)],"#a9dfee")
    surface([(p.outlet_radius,p.hub_outlet_z),(p.outlet_radius,p.front_outlet_z)],"#6b9eb8")
    if section:
        for angle in (0,pi):
            hub,front = data["hub"][::4],data["front"][::4]
            if hub[-1] != data["hub"][-1]:
                hub.append(data["hub"][-1]); front.append(data["front"][-1])
            strips = [[data["vertices"][0],hub[0],front[0],data["vertices"][-1]]]
            strips += [[hub[i],hub[i+1],front[i+1],front[i]] for i in range(len(hub)-1)]
            for strip in strips:
                faces.append([(r*np.cos(angle),r*np.sin(angle),z) for r,z in strip])
                colors.append("#e7b56d")
    # Sort individual faces across all surfaces, avoiding per-surface occlusion artifacts.
    collection = Poly3DCollection(faces,facecolors=colors,linewidths=0,shade=True,zsort="average")
    collection.set_edgecolor("none")
    ax.add_collection3d(collection)
    for r,z in ((p.hub_inlet_radius,p.inlet_z),(p.front_inlet_radius,p.inlet_z),
                (p.outlet_radius,p.front_outlet_z),(p.outlet_radius,p.hub_outlet_z)):
        ax.plot(r*np.cos(theta),r*np.sin(theta),np.full_like(theta,z),color="#34566a",linewidth=0.8)
    ax.set(xlabel="X (mm)",ylabel="Y (mm)",zlabel="Z (mm)",
           xlim=(-p.outlet_radius,p.outlet_radius),ylim=(-p.outlet_radius,p.outlet_radius),
           zlim=(p.hub_outlet_z,p.inlet_z))
    ax.set_box_aspect((2*p.outlet_radius,2*p.outlet_radius,p.inlet_z-p.hub_outlet_z))
    ax.view_init(elev=35.2643897,azim=-45)
    ax.set_title("ImpellerFluidBlank: "+("axial half-section" if section else "isometric")
                 +"\nProfile-derived preview (not a CAD screenshot)")
    fig.tight_layout()
    destination.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(destination,dpi=180)
    plt.close(fig)


def spaceclaim_script(p):
    data = profile_data(p)
    return '''# Python Script, API Version = V241
from System import Array
from System.Collections.Generic import List
import clr
hub_controls = %r
front_controls = %r
inlet_z = %r
tolerance = %r
stage = "create_exact_bezier_curves"
geometry_point_type = Point.Origin.GetType()
geometry_namespace = geometry_point_type.Namespace
geometry_assembly = geometry_point_type.Assembly
api_context = geometry_namespace + " :: " + geometry_assembly.FullName

def geometry_type(name):
    actual_type = geometry_assembly.GetType(geometry_namespace + "." + name)
    if actual_type is None:
        raise Exception("Loaded geometry API has no type: " + name)
    return clr.GetPythonType(actual_type)

# Resolve types in the actual installed API instead of hardcoding V241 CLR names.
NurbsData = geometry_type("NurbsData")
NurbsCurve = geometry_type("NurbsCurve")
Knot = geometry_type("Knot")
ControlPoint = geometry_type("ControlPoint")
root = GetRootPart()
if root.Bodies.Count != 0 or root.Components.Count != 0:
    raise Exception("Expected a new empty design")

def point(rz):
    return Point.Create(MM(rz[0]), 0, MM(rz[1]))

def make_bezier(controls):
    # Order 4 = degree 3. Clamped knot multiplicities 4,4 and unit weights
    # express the specified cubic Bezier exactly, not through-point fitting.
    knots = Array[Knot]([Knot(0.0,4),Knot(1.0,4)])
    data = NurbsData(4,False,False,knots)
    poles = Array[ControlPoint]([ControlPoint(point(rz),1.0) for rz in controls])
    creator = getattr(NurbsCurve,"CreateFromControlPoints",None)
    if creator is None:
        creator = NurbsCurve.Create
    curve = creator(data,poles)
    maximum_error = 0.0
    for i in range(501):
        u = i/500.0
        weights = ((1-u)**3,3*(1-u)**2*u,3*(1-u)*u*u,u**3)
        expected = [MM(sum(w*rz[k] for w,rz in zip(weights,controls))) for k in (0,1)]
        actual = curve.Evaluate(u).Point
        error = ((actual.X-expected[0])**2+actual.Y**2+(actual.Z-expected[1])**2)**0.5
        maximum_error = max(maximum_error,error)
    if maximum_error > MM(tolerance):
        raise Exception("Native Bezier differs from formula: %%s mm" %% (maximum_error*1000))
    return CurveSegment.Create(curve,Interval.Create(0,1)),maximum_error*1000

hub,error_hub = make_bezier(hub_controls)
front_reverse,error_front = make_bezier(tuple(reversed(front_controls)))
stage = "close_bezier_meridional_profile"
curves = List[ITrimmedCurve]()
a,b,c = (hub_controls[0][0],inlet_z),hub_controls[0],hub_controls[-1]
d,e,f = (front_controls[0][0],inlet_z),front_controls[0],front_controls[-1]
curves.Add(CurveSegment.Create(point(a),point(b)))
curves.Add(hub)
curves.Add(CurveSegment.Create(point(c),point(f)))
curves.Add(front_reverse)
curves.Add(CurveSegment.Create(point(e),point(d)))
curves.Add(CurveSegment.Create(point(d),point(a)))
plane = Plane.Create(Frame.Create(Point.Origin,Direction.DirX,Direction.DirZ))
profile = PlanarBody.Create(plane,curves).CreatedBody
if profile.Faces.Count != 1:
    raise Exception("Profile must create exactly one closed planar face")
stage = "revolve_bezier_meridional_profile"
options = RevolveFaceOptions()
options.ExtrudeType = ExtrudeType.Add
RevolveFaces.Execute(FaceSelection.Create(profile.Faces[0]),
                    Line.Create(Point.Origin,Direction.DirZ),DEG(360),options)
stage = "check_closed_bezier_blank"
if root.Bodies.Count != 1:
    raise Exception("Expected exactly one revolved body")
body = root.Bodies[0]
if not body.Shape.IsClosed or body.Shape.PieceCount != 1 or body.Shape.Volume <= 0:
    raise Exception("Blank must be one connected closed positive-volume solid")
body.SetName("ImpellerFluidBlank")
result = {"closed":"True","body_count":"1","pieces":"1",
          "curve_error_mm":str(max(error_hub,error_front))}
''' % (data["hub_controls"],data["front_controls"],p.inlet_z,p.curve_tolerance_mm)


def build_and_check(modeler, p=None):
    """新建并验证此 R140 Bezier 毛坯，返回可供后续叶片脚本操作的 Design。"""
    from ansys.geometry.core.connection.backend import BackendType
    p = Parameters() if p is None else p
    check_profile(p)
    if not modeler.client.healthy or modeler.client.backend_type != BackendType.SPACECLAIM:
        raise RuntimeError("需要 Healthy 的 SPACECLAIM 后端；不根据 Version=0.0.0 判错")
    modeler.create_design("impeller_blank")
    values,design = run_checked_script(modeler,spaceclaim_script(p),"bezier_meridional_blank")
    if (values.get("closed") != "True" or values.get("body_count") != "1"
            or values.get("pieces") != "1"
            or float(values.get("curve_error_mm","inf")) > p.curve_tolerance_mm):
        raise RuntimeError(f"原生曲线/封闭性检查未通过：{values}")
    bodies = design.get_all_bodies()
    if len(bodies) != 1 or bodies[0].name != "ImpellerFluidBlank" or bodies[0].is_surface:
        raise RuntimeError("应只有一个名为 ImpellerFluidBlank 的流体实体")
    actual = bodies[0].volume.to("mm^3").magnitude
    target = expected_volume(p)
    if not isfinite(actual) or not isclose(actual,target,rel_tol=1e-6,abs_tol=0.01):
        raise RuntimeError(f"CAD 体积 {actual} 与精确 Bezier 回转积分 {target} mm³ 不符")
    print(f"曲线误差 {values['curve_error_mm']} mm；封闭单实体，体积 {actual:.6f} mm³")
    return design


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host",default=os.getenv("ANSRV_GEO_HOST"))
    parser.add_argument("--port",type=int,default=int(os.getenv("ANSRV_GEO_PORT","50051")))
    parser.add_argument("--output",type=Path,default=Path(__file__).resolve().with_name("impeller_blank.step"))
    parser.add_argument("--preview",type=Path,default=Path(__file__).resolve().with_name("impeller_meridional.png"))
    parser.add_argument("--isometric",type=Path,default=Path(__file__).resolve().with_name("impeller_blank_isometric.png"))
    parser.add_argument("--section-preview",type=Path,default=Path(__file__).resolve().with_name("impeller_blank_section.png"))
    parser.add_argument("--preview-only",action="store_true")
    args = parser.parse_args()
    if args.output.suffix.lower() not in (".step",".stp"):
        parser.error("--output 必须使用 .step 或 .stp 扩展名")
    p = Parameters()
    checks = check_profile(p)
    save_preview(p,args.preview.resolve())
    save_3d_preview(p,args.isometric.resolve())
    save_3d_preview(p,args.section_preview.resolve(),section=True)
    print(f"轮廓闭合、无自交、正宽度、切向检查通过：{checks}")
    print(f"精确回转积分体积：{expected_volume(p):.6f} mm³")
    print(f"已生成子午面、等轴测和轴向剖切参数预览：{args.preview.resolve().parent}")
    if args.preview_only:
        print("仅离线预览；未运行 CAD 或生成新的 STEP。")
        return
    from ansys.geometry.core.designer.design import DesignFileFormat
    modeler = connect(args.host,args.port)
    try:
        design = build_and_check(modeler,p)
        destination = args.output.resolve()
        design.download(str(destination),format=DesignFileFormat.STEP)
        if not destination.is_file() or destination.stat().st_size == 0:
            raise RuntimeError("STEP 导出失败：文件不存在或为空")
        print(f"已导出新的 Bezier 毛坯：{destination}")
    finally:
        modeler.close()


if __name__ == "__main__":
    main()
