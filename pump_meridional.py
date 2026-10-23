"""精确子午面旋转得到无叶片 ImpellerFluidBlank；尺寸单位 mm。

python pump_meridional.py                  # 本机 SpaceClaim 241
python pump_meridional.py --host HOST --port 50051
python pump_meridional.py --preview-only   # 离线校验截面并生成预览

24.1 不支持 PyAnsys revolve_sketch（该接口要求 24R2）；这里通过支持
24.1 的 Modeler.run_script_file 执行 SpaceClaim 原生 PlanarBody.Create /
RevolveFaces.Execute。圆弧是精确四分之一圆，不以折线或圆柱替代。
只创建本毛坯；中心 r<10 mm 始终为非流体空间。
"""

import argparse
from dataclasses import dataclass, fields
from math import cos, isclose, isfinite, pi, sin
import os
from pathlib import Path

from main import connect
from spaceclaim_diagnostics import run_checked_script


@dataclass(frozen=True)
class Parameters:
    hub_inlet_radius: float = 10.0
    front_inlet_radius: float = 27.5
    inlet_z: float = 40.0
    hub_bend_radius: float = 15.0
    front_bend_radius: float = 8.0
    hub_outlet_z: float = 0.0
    front_outlet_z: float = 10.0
    outlet_radius: float = 50.0

    @property
    def hub_center(self):
        return (self.hub_inlet_radius + self.hub_bend_radius,
                self.hub_outlet_z + self.hub_bend_radius)

    @property
    def front_center(self):
        return (self.front_inlet_radius + self.front_bend_radius,
                self.front_outlet_z + self.front_bend_radius)

    def validate(self):
        for field in fields(self):
            value = getattr(self, field.name)
            if not isfinite(value):
                raise ValueError(f"{field.name} 必须为有限数值")
        if min(self.hub_inlet_radius, self.front_inlet_radius,
               self.hub_bend_radius, self.front_bend_radius) <= 0:
            raise ValueError("进口半径和转弯半径必须为正，截面不得接触旋转轴")
        # 这些充分条件保证两条单调侧边不相交，且各直线段均有正长度。
        if not (self.hub_center[0] < self.front_inlet_radius
                < self.front_center[0] < self.outlet_radius):
            raise ValueError("需满足轮毂弧圆心半径 < 前侧进口半径 < 前侧弧圆心半径 < 出口半径")
        if not (self.hub_outlet_z < self.front_outlet_z
                and self.inlet_z > max(self.hub_center[1], self.front_center[1])):
            raise ValueError("出口宽度须为正，进口高度须高于两个圆弧圆心")
        return self

    def vertices(self):
        """按闭环顺序：轮毂进口→轮毂侧→出口→前侧→前侧进口。"""
        self.validate()
        return [
            (self.hub_inlet_radius, self.inlet_z),
            (self.hub_inlet_radius, self.hub_center[1]),
            (self.hub_center[0], self.hub_outlet_z),
            (self.outlet_radius, self.hub_outlet_z),
            (self.outlet_radius, self.front_outlet_z),
            (self.front_center[0], self.front_outlet_z),
            (self.front_inlet_radius, self.front_center[1]),
            (self.front_inlet_radius, self.inlet_z),
        ]


def expected_volume(p):
    """Green 公式 V=pi*闭环积分(r² dz)，精确包含两段四分之一圆弧。"""
    p.validate()

    def arc_integral(center_r, radius):
        # t=pi..3pi/2，r=c+R*cos(t)，dz=R*cos(t)dt。
        return -radius*center_r**2 + center_r*radius**2*pi/2 - 2*radius**3/3

    moment = (
        p.outlet_radius**2 * (p.front_outlet_z - p.hub_outlet_z)
        - p.hub_inlet_radius**2 * (p.inlet_z - p.hub_center[1])
        + p.front_inlet_radius**2 * (p.inlet_z - p.front_center[1])
        + arc_integral(p.hub_center[0], p.hub_bend_radius)
        - arc_integral(p.front_center[0], p.front_bend_radius)
    )
    volume = pi * moment
    if not isfinite(volume) or volume <= 0:
        raise ValueError("截面旋转体积必须为正")
    return volume


def save_preview(p, destination):
    """仅预览使用采样；后端几何使用精确圆弧。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    vertices = p.vertices()
    outline = vertices[:2]
    for center, radius, start, end, tail in (
        (p.hub_center, p.hub_bend_radius, pi, 1.5*pi, vertices[3:6]),
        (p.front_center, p.front_bend_radius, 1.5*pi, pi, vertices[7:]),
    ):
        for i in range(1, 181):
            t = start + (end - start)*i/180
            outline.append((center[0] + radius*cos(t), center[1] + radius*sin(t)))
        outline.extend(tail)
    outline.append(vertices[0])
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.fill(*zip(*outline), facecolor="#bce6f8", edgecolor="#146c94", linewidth=2)
    ax.plot(*zip(*vertices), "o", color="#146c94", markersize=4)
    for i, (r, z) in enumerate(vertices):
        offset = (-8, 10) if i in (0, 6, 7) else (8, -15)
        ax.annotate(f"({r:g}, {z:g})", (r, z), xytext=offset,
                    textcoords="offset points", fontsize=9)
    for center, radius in ((p.hub_center, p.hub_bend_radius),
                           (p.front_center, p.front_bend_radius)):
        ax.plot(*center, "+", color="#777777")
        ax.annotate(f"R{radius:g}", center, xytext=(8, 5),
                    textcoords="offset points", color="#666666")
    ax.axvline(0, color="#d05a42", linestyle="--", label="Z rotation axis (r=0)")
    ax.set(xlabel="r = X (mm)", ylabel="Z (mm)", title="ImpellerFluidBlank: meridional section")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-5, p.outlet_radius + 10)
    ax.set_ylim(p.hub_outlet_z - 7, p.inlet_z + 8)
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right")
    fig.tight_layout()
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=180)
    plt.close(fig)


def spaceclaim_script(p):
    """生成在 SpaceClaim 的 V241 脚本环境内执行的 IronPython 代码。"""
    vertices = p.vertices()
    # 使用数值 repr 嵌入参数，避免依赖远程脚本参数字典或本机路径。
    return '''# Python Script, API Version = V241
from System.Collections.Generic import List
stage = "create_meridional_profile"
vertices = %r
hub_center = %r
front_center = %r

def point(rz):
    return Point.Create(MM(rz[0]), MM(0), MM(rz[1]))

root = GetRootPart()
if root.Bodies.Count != 0 or root.Components.Count != 0:
    raise Exception("Expected a new empty design")
curves = List[ITrimmedCurve]()
for i in range(8):
    start, end = point(vertices[i]), point(vertices[(i + 1) %% 8])
    if i == 1:
        # XZ local (r,z) normal is -Y: left to bottom, CCW quarter arc.
        curves.Add(CurveSegment.CreateArc(point(hub_center), start, end, -Direction.DirY))
    elif i == 5:
        # Traverse front-side arc in reverse: bottom to left.
        curves.Add(CurveSegment.CreateArc(point(front_center), start, end, Direction.DirY))
    else:
        curves.Add(CurveSegment.Create(start, end))
plane = Plane.Create(Frame.Create(Point.Origin, Direction.DirX, Direction.DirZ))
profile = PlanarBody.Create(plane, curves).CreatedBody
if profile.Faces.Count != 1:
    raise Exception("Meridional profile must form one closed planar face")
options = RevolveFaceOptions()
options.ExtrudeType = ExtrudeType.Add
stage = "revolve_meridional_profile"
RevolveFaces.Execute(FaceSelection.Create(profile.Faces[0]),
                    Line.Create(Point.Origin, Direction.DirZ), DEG(360), options)
if root.Bodies.Count != 1:
    raise Exception("Revolve must produce exactly one body, with no sheet/tool bodies")
body = root.Bodies[0]
body.SetName("ImpellerFluidBlank")
stage = "check_closed_blank"
if not body.Shape.IsClosed:
    raise Exception("Revolved body is not closed")
result = {"closed": "True", "body_count": "1"}
''' % (vertices, p.hub_center, p.front_center)


def build_and_check(modeler, p):
    from ansys.geometry.core.connection.backend import BackendType

    if not modeler.client.healthy or modeler.client.backend_type != BackendType.SPACECLAIM:
        raise RuntimeError("需要 Healthy 的 SPACECLAIM 后端；不根据 Version=0.0.0 判错")
    modeler.create_design("impeller_blank")
    # 新建文档后运行原生旋转，再通过 PyAnsys 刷新并检查实体。
    values, design = run_checked_script(modeler, spaceclaim_script(p), "meridional_blank")
    if values.get("closed") != "True" or values.get("body_count") != "1":
        raise RuntimeError(f"原生封闭性检查未通过：{values}")
    bodies = design.get_all_bodies()
    if len(bodies) != 1 or bodies[0].name != "ImpellerFluidBlank" or bodies[0].is_surface:
        raise RuntimeError("应只有一个名为 ImpellerFluidBlank 的实体")
    actual = bodies[0].volume.to("mm^3").magnitude
    expected = expected_volume(p)
    if not isfinite(actual) or not isclose(actual, expected, rel_tol=1e-6, abs_tol=1e-3):
        raise RuntimeError(f"体积异常：CAD={actual}，解析值={expected} mm³")
    print(f"封闭性和单实体检查通过；体积 {actual:.6f} mm³（解析值 {expected:.6f}）")
    return design


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default=os.getenv("ANSRV_GEO_HOST"))
    parser.add_argument("--port", type=int, default=int(os.getenv("ANSRV_GEO_PORT", "50051")))
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).resolve().with_name("impeller_blank.step"))
    parser.add_argument("--preview", type=Path,
                        default=Path(__file__).resolve().with_name("impeller_meridional.png"))
    parser.add_argument("--preview-only", action="store_true")
    args = parser.parse_args()
    if args.output.suffix.lower() not in (".step", ".stp"):
        parser.error("--output 必须使用 .step 或 .stp 扩展名")
    p = Parameters().validate()
    volume = expected_volume(p)
    save_preview(p, args.preview.resolve())
    print(f"子午面预览：{args.preview.resolve()}；解析体积：{volume:.6f} mm³")
    if args.preview_only:
        print("离线截面检查通过；未连接 CAD 后端，未导出 STEP。")
        return

    from ansys.geometry.core.designer.design import DesignFileFormat

    modeler = connect(args.host, args.port)
    try:
        design = build_and_check(modeler, p)
        destination = args.output.resolve()
        design.download(str(destination), format=DesignFileFormat.STEP)
        if not destination.is_file() or destination.stat().st_size == 0:
            raise RuntimeError("STEP 导出失败：文件不存在或为空")
        print(f"已导出：{destination}")
    finally:
        modeler.close()


if __name__ == "__main__":
    main()
