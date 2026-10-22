"""2024R1 简化离心泵流体域；所有长度参数的数值单位均为 mm。

python pump_stage1.py                     # 本机 SpaceClaim 241
python pump_stage1.py --host HOST --port 50051
python pump_stage1.py --check-only        # 无后端：参数和解析体积检查

Z 为旋转轴，进口流向 +Z，出口流向 +X。默认轴向区间：
进口 [-80, 0]，叶轮 [0, 10]，后盖板留空 [10, 13]，
后腔 [13, 15.5]，密封腔 [15.5, 27.5]。
轴/后盖板只通过流体域边界表达，不创建金属实体。

未指定的集流腔采用 R50~R65、Z=-11~21 的等截面圆环，外圈同时
充当通向后腔的旁通通道。出口连接段采用 32x32 mm 方截面直段，
在 X=85 平面接直径 32 mm 的圆管；方圆连接是本阶段的简化台阶。
连接段与圆环合并为 VoluteFluid；不合并六个命名流体域。

独立域之间仅共边界：进口/叶轮在 Z=0，叶轮/集流腔在 R=50，
集流腔/出口在 X=85，集流腔/后腔在 R=50，后腔/密封腔在 Z=15.5。
不使用 24R2 才有的 get_collision，也不调用新版共享拓扑等接口。
STEP 保留六个独立体；网格所需的共形共享拓扑需在后续阶段处理。
"""

import argparse
from dataclasses import dataclass, fields
from math import asin, isclose, isfinite, pi, sqrt
import os
from pathlib import Path

from main import connect


@dataclass(frozen=True)
class Parameters:
    inlet_diameter: float = 40.0
    inlet_length: float = 80.0
    impeller_diameter: float = 100.0
    shaft_diameter: float = 20.0
    impeller_width: float = 10.0
    collector_radial_width: float = 15.0  # 未指定尺寸的简化假设
    collector_height: float = 32.0
    connector_length: float = 20.0  # 集流腔外半径至出口管起点
    connector_embed: float = 10.0   # 方形直段嵌入圆环的长度
    outlet_diameter: float = 32.0
    outlet_length: float = 100.0
    rear_gap: float = 2.5
    backplate_thickness: float = 3.0
    seal_inner_diameter: float = 26.0
    seal_outer_diameter: float = 34.0
    seal_length: float = 12.0

    @property
    def rear_start(self):
        return self.impeller_width + self.backplate_thickness

    @property
    def rear_end(self):
        return self.rear_start + self.rear_gap

    @property
    def collector_radius(self):
        return self.impeller_diameter / 2 + self.collector_radial_width

    @property
    def outlet_start(self):
        return self.collector_radius + self.connector_length

    def validate(self):
        for field in fields(self):
            value = getattr(self, field.name)
            if not isfinite(value) or value <= 0:
                raise ValueError(f"{field.name} 必须是有限正数（mm）")
        if not self.shaft_diameter < self.inlet_diameter < self.impeller_diameter:
            raise ValueError("需满足：轴径 < 进口直径 < 叶轮直径")
        if not (self.shaft_diameter < self.seal_inner_diameter
                < self.seal_outer_diameter < self.impeller_diameter):
            raise ValueError("需满足：轴径 < 密封内径 < 密封外径 < 叶轮直径")
        center_z = self.impeller_width / 2
        if (center_z - self.collector_height / 2 > 0
                or center_z + self.collector_height / 2 < self.rear_end
                or self.collector_height < self.outlet_diameter):
            raise ValueError("集流腔轴向范围必须覆盖叶轮、后腔和整个出口截面")
        start = self.collector_radius - self.connector_embed
        half_side = self.outlet_diameter / 2
        if not (start > self.impeller_diameter / 2
                and sqrt(start**2 + half_side**2) < self.collector_radius):
            raise ValueError("连接段必须在叶轮半径外，并在整个宽度上嵌入集流圆环")
        return self


def expected_volumes(p):
    """解析体积，mm³；圆环与方连接段的相交体积只计一次。"""
    p.validate()
    ri, rs = p.impeller_diameter / 2, p.shaft_diameter / 2
    radius, a = p.collector_radius, p.outlet_diameter / 2
    x0 = radius - p.connector_embed
    # 积分 y=-a..a 的 sqrt(R²-y²)-x0；验证条件保证无需分段。
    intersection_area = a * sqrt(radius**2 - a**2) + radius**2 * asin(a / radius) - 2*a*x0
    overlap_height = p.outlet_diameter
    ring_volume = pi * (radius**2 - ri**2) * p.collector_height
    connector_volume = (p.connector_embed + p.connector_length) * p.outlet_diameter**2
    return {
        "InletFluid": pi * (p.inlet_diameter / 2)**2 * p.inlet_length,
        "ImpellerFluid": pi * (ri**2 - rs**2) * p.impeller_width,
        "VoluteFluid": ring_volume + connector_volume - intersection_area * overlap_height,
        "OutletFluid": pi * a**2 * p.outlet_length,
        "RearChamberFluid": pi * (ri**2 - rs**2) * p.rear_gap,
        "SealChamberFluid": pi * ((p.seal_outer_diameter / 2)**2
                                  - (p.seal_inner_diameter / 2)**2) * p.seal_length,
    }


def build_pump(modeler, p):
    """只用基本草图、拉伸、差集与并集；工具体由布尔操作消耗。"""
    from ansys.geometry.core.math import Plane, Point2D, Point3D
    from ansys.geometry.core.misc import UNITS
    from ansys.geometry.core.sketch import Sketch

    p.validate()
    mm = UNITS.mm
    design = modeler.create_design("pump_stage1")

    def plane_z(z):
        return Plane(origin=Point3D([0, 0, z], unit=mm))

    def plane_x(x):
        return Plane(origin=Point3D([x, 0, p.impeller_width / 2], unit=mm),
                     direction_x=[0, 1, 0], direction_y=[0, 0, 1])

    def cylinder(name, plane, radius, length):
        sketch = Sketch(plane)
        sketch.circle(Point2D([0, 0], unit=mm), radius * mm)
        return design.extrude_sketch(name, sketch, distance=length * mm)

    def annulus(name, z, outer_radius, inner_radius, length):
        body = cylinder(name, plane_z(z), outer_radius, length)
        tool = cylinder(name + "_BoreTool", plane_z(z), inner_radius, length)
        body.subtract(tool)  # keep_other=False，刀具体被消耗
        return body

    ri, rs = p.impeller_diameter / 2, p.shaft_diameter / 2
    cylinder("InletFluid", plane_z(-p.inlet_length), p.inlet_diameter / 2, p.inlet_length)
    annulus("ImpellerFluid", 0, ri, rs, p.impeller_width)
    volute = annulus("VoluteFluid", (p.impeller_width - p.collector_height) / 2,
                     p.collector_radius, ri, p.collector_height)
    # 仅合并属于 VoluteFluid 的连接段，保留与其他流体体的接口。
    x0 = p.collector_radius - p.connector_embed
    connector_sketch = Sketch(plane_x(x0))
    connector_sketch.box(Point2D([0, 0], unit=mm),
                         p.outlet_diameter * mm, p.outlet_diameter * mm)
    connector = design.extrude_sketch("CollectorConnectorTool", connector_sketch,
                                     distance=(p.outlet_start - x0) * mm)
    volute.unite(connector)
    cylinder("OutletFluid", plane_x(p.outlet_start), p.outlet_diameter / 2, p.outlet_length)
    annulus("RearChamberFluid", p.rear_start, ri, rs, p.rear_gap)
    annulus("SealChamberFluid", p.rear_end, p.seal_outer_diameter / 2,
             p.seal_inner_diameter / 2, p.seal_length)
    return design


def check_bodies(design, p):
    """导出前检查最终六体及逐体体积，拒绝残留工具体、零体积或错误尺寸。"""
    expected = expected_volumes(p)
    bodies = design.get_all_bodies()
    if len(bodies) != 6 or {b.name for b in bodies} != set(expected):
        raise RuntimeError(f"应有六个命名流体体，实际为 {[b.name for b in bodies]}")
    for body in bodies:
        actual = body.volume.to("mm^3").magnitude
        target = expected[body.name]
        if body.is_surface or not isfinite(actual) or actual <= 0:
            raise RuntimeError(f"{body.name} 不是有效的正体积实体")
        if not isclose(actual, target, rel_tol=1e-6, abs_tol=1e-3):
            raise RuntimeError(f"{body.name} 体积异常：{actual}，解析值 {target} mm³")
        print(f"{body.name}: {actual:.6f} mm³（解析值 {target:.6f}）")
    print("六体数量、名称及体积检查通过；接口位置由参数约束保证。")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default=os.getenv("ANSRV_GEO_HOST"))
    parser.add_argument("--port", type=int, default=int(os.getenv("ANSRV_GEO_PORT", "50051")))
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).resolve().with_name("pump_stage1.step"))
    parser.add_argument("--check-only", action="store_true", help="仅检查参数和解析体积，不连接后端")
    args = parser.parse_args()
    p = Parameters().validate()  # 全部尺寸集中于 Parameters；改尺寸后重新运行。
    if args.check_only:
        for name, volume in expected_volumes(p).items():
            print(f"{name}: {volume:.6f} mm³（解析值，未运行 CAD 后端）")
        return
    if args.output.suffix.lower() not in (".step", ".stp"):
        parser.error("--output 必须是 .step 或 .stp 文件")

    from ansys.geometry.core.designer.design import DesignFileFormat

    modeler = connect(args.host, args.port)
    try:
        # 远程方式也明确验证后端为用户要求的 24.1。
        version = modeler.client.backend_version
        if (version.major, version.minor) != (24, 1):
            raise RuntimeError(f"需要 2024R1（24.1）后端，实际为 {version}")
        design = build_pump(modeler, p)
        check_bodies(design, p)
        destination = args.output.resolve()
        design.download(str(destination), format=DesignFileFormat.STEP)
        if not destination.is_file() or destination.stat().st_size == 0:
            raise RuntimeError("STEP 导出失败：文件不存在或为空")
        print(f"已导出六体 STEP：{destination}")
    finally:
        modeler.close()


if __name__ == "__main__":
    main()
