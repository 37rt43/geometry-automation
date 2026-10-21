"""使用本机 SpaceClaim 2024R1 生成板件，导出脚本同目录的 plate.step。

运行：python plate_with_hole.py
依赖：ansys-geometry-core；实际建模需要 Windows 本机 SpaceClaim 和许可证。
24.1 没有 PyAnsys Geometry 圆角操作接口，因此使用四段精确的 R2 圆弧
构造 YZ 截面并沿 X 拉伸，几何上等价于对长方体的四条 X 向长边倒圆角。
不调用新版圆角、切割或脚本执行接口。
"""

from pathlib import Path


def connect():
    """启动并连接本机 241 产品，显式使用 24.1 API。"""
    from ansys.geometry.core.connection import launch_modeler_with_spaceclaim
    from ansys.geometry.core.connection.backend import ApiVersions

    return launch_modeler_with_spaceclaim(version=241, api_version=ApiVersions.V_241)


def build_plate(modeler, destination: Path):
    """外包尺寸 100×50×10 mm；中心 Z 向通孔 Ø10；四条长边 R2。"""
    from ansys.geometry.core.designer.design import DesignFileFormat
    from ansys.geometry.core.math import Plane, Point2D, Point3D
    from ansys.geometry.core.misc import UNITS
    from ansys.geometry.core.sketch import Sketch

    length, width, thickness = 100.0, 50.0, 10.0
    hole_diameter, fillet_radius = 10.0, 2.0
    overshoot = 1.0  # 刀具体两端各伸出 1 mm，保证通孔。
    mm = UNITS.mm
    half_w, half_t, radius = width / 2, thickness / 2, fillet_radius

    def point(y, z):
        return Point2D([y, z], unit=mm)

    design = modeler.create_design("plate_with_hole")
    # 截面的局部 x/y 分别对应全局 Y/Z，法向为全局 +X。
    section = Sketch(Plane(
        origin=Point3D([-length / 2, 0, 0], unit=mm),
        direction_x=[0, 1, 0],
        direction_y=[0, 0, 1],
    ))
    # 从下边开始逆时针构造封闭轮廓；四段圆弧都是 90 度。
    a = point(-half_w + radius, -half_t)
    b = point(half_w - radius, -half_t)
    c = point(half_w, -half_t + radius)
    d = point(half_w, half_t - radius)
    e = point(half_w - radius, half_t)
    f = point(-half_w + radius, half_t)
    g = point(-half_w, half_t - radius)
    h = point(-half_w, -half_t + radius)
    section.segment(a, b)
    section.arc(b, c, point(half_w - radius, -half_t + radius))
    section.segment(c, d)
    section.arc(d, e, point(half_w - radius, half_t - radius))
    section.segment(e, f)
    section.arc(f, g, point(-half_w + radius, half_t - radius))
    section.segment(g, h)
    section.arc(h, a, point(-half_w + radius, -half_t + radius))
    plate = design.extrude_sketch("plate", section, distance=length * mm)

    hole = Sketch(Plane(origin=Point3D([0, 0, -half_t - overshoot], unit=mm)))
    hole.circle(Point2D([0, 0], unit=mm), hole_diameter / 2 * mm)
    cutter = design.extrude_sketch(
        "hole_tool", hole, distance=(thickness + 2 * overshoot) * mm
    )
    plate.subtract(cutter)  # 24.1 布尔差集；默认消耗刀具体。
    design.download(str(destination), format=DesignFileFormat.STEP)
    if not destination.is_file() or destination.stat().st_size == 0:
        raise RuntimeError("STEP 导出失败：未生成非空文件。")
    return destination


def main():
    destination = Path(__file__).resolve().with_name("plate.step")
    modeler = connect()
    try:
        print("已导出：", build_plate(modeler, destination))
    finally:
        modeler.close()


if __name__ == "__main__":
    main()
