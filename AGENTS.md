# AGENTS.md —— 给 Codex / AI 编程助手的项目说明

## 项目目标
使用 PyAnsys Geometry（`ansys-geometry-core`，GitHub: ansys/pyansys-geometry）通过 Python
调用 Ansys Geometry Service 完成参数化三维建模（草图、拉伸、旋转、布尔运算、倒角、导入导出 CAD 等），
替代手工 SpaceClaim 建模。

## 环境
- Python 3.10+（最新版 ansys-geometry-core 要求 3.12+）
- 依赖见 `requirements.txt`（含 graphics 扩展，用于草图/模型预览）：`pip install -r requirements.txt`
- **本机后端版本固定为 Ansys 2024R1（v241，DMS API 24.1）**，Windows 本地 SpaceClaim/Discovery 服务。
  写代码时只能使用在 24.1 上可用的 API，禁止调用标注 "available starting on 24R2/25R1/25R2" 的
  方法和属性；如某个功能确实需要新版本，在代码注释中明确标注，不要静默使用。
- 连接方式：
  - 本机 Windows：`launch_modeler()`（product_version 按本机安装，如 `launch_modeler(product_version="241")`）。
  - 远程：`Modeler(host=..., port=...)` 连接。
  - **禁止使用 Linux Docker Core Service 方案**（仅支持 2025R2+，本机版本不满足）。
  - 云端 Codex 沙箱内没有 Geometry Service 和许可证，代码只能做静态检查；
    真正建模必须连接本机后端。

## 常用 API
- 连接：`from ansys.geometry.core import Modeler` / `from ansys.geometry.core.connection import launch_modeler`
  本机启动显式指定版本：`launch_modeler(product_version="241")``
- 草图：`ansys.geometry.core.sketch.Sketch`（circle、polygon、segment、arc、slot…）
- 平面与点：`ansys.geometry.core.math.Plane / Point3D / Point2D / Matrix`
- 单位：`ansys.geometry.core.misc.UNITS`，所有尺寸必须带单位（如 `10 * UNITS.mm`）
- 设计：`modeler.create_design(name)`，再 `design.extrude_sketch / revolve_sketch / add_…`
- 不要修改预置的 BGE-M3 embedding 模型配置（那是 AALI 相关功能，与建模无关）。

## 代码约定
- 连接逻辑统一放在 `main.py` 的 `connect()`。
- 每个建模任务一个函数或一个脚本，参数化变量集中放在函数开头，方便改尺寸。
- 主机/端口走环境变量 `ANSRV_GEO_HOST` / `ANSRV_GEO_PORT`，不要硬编码许可证服务器。
- 导出文件格式用 SCDX 或 STEP，便于导入 Mechanical / SpaceClaim。

## 验证
- 不依赖后端：`python -m py_compile main.py`
- 连接验证：`python main.py --host <服务器IP> --port 50051`
