"""
PyAnsys Geometry 最小示例：连接 Geometry Service，画草图并拉伸成实体。

三种用法：
1. 本机 / 本地 Docker 自动启动（需要本机装好 Ansys 2025R2+，或有 Docker 镜像）：
   python main.py

2. 连接已经启动的远程 / Docker Geometry 服务（默认端口 50051）：
   python main.py --host 192.168.1.10 --port 50051

3. 手动启动 Docker 容器：
   docker run --name ans_geo -e LICENSE_SERVER=<license服务器IP> \
              -p 50051:50051 ghcr.io/ansys/geometry:core-linux-latest
"""

import argparse
import os


def connect(host: str | None = None, port: int = 50051):
    if host:
        # 连接已有的 Geometry 服务（远程机器 / Docker / 手动启动）
        from ansys.geometry.core import Modeler

        return Modeler(host=host, port=port)

    # 自动启动本地 Geometry 服务（本机为 Ansys 2024R1，显式指定 241）
    from ansys.geometry.core import launch_modeler_with_spaceclaim
    from ansys.geometry.core.connection.backend import ApiVersions

    return launch_modeler_with_spaceclaim(version=241, api_version=ApiVersions.V_241)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--host", default=os.getenv("ANSRV_GEO_HOST"),
                        help="Geometry 服务地址，不填则本机自动启动")
    parser.add_argument("--port", type=int,
                        default=int(os.getenv("ANSRV_GEO_PORT", "50051")))
    args = parser.parse_args()

    modeler = connect(args.host, args.port)
    print(modeler)  # 打印服务信息，确认连接成功

    # ---- 下面是一个最简单的建模示例：圆草图拉伸成圆柱 ----
    from ansys.geometry.core.math import Plane, Point3D, Point2D
    from ansys.geometry.core.misc import UNITS
    from ansys.geometry.core.sketch import Sketch

    design = modeler.create_design("demo_part")

    # 在 z=0 平面建草图
    plane = Plane(Point3D([0, 0, 0]))
    sketch = Sketch(plane)
    sketch.circle(Point2D([0, 0]), 30 * UNITS.mm)   # 半径 30mm 的圆

    # 拉伸 100mm 成实体
    body = design.extrude_sketch("cylinder", sketch, distance=100 * UNITS.mm)
    print("创建的 body：", body.name)

    # 导出为 SCDX（SpaceClaim 格式，供后续在 Mechanical 中使用）
    # design.upload("demo_part.scdx", body)

    modeler.close()


if __name__ == "__main__":
    main()
#（注：内容由AI生成）
