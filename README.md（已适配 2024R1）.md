# geometry-automation

基于 [PyAnsys Geometry](https://github.com/ansys/pyansys-geometry) 的参数化三维建模脚本仓库，配合 ChatGPT Codex 云端环境使用。

## 1. 安装依赖

```bash
pip install -r requirements.txt
```

## 2. 准备 Geometry Service 后端

`ansys-geometry-core` 是客户端，必须连接 Geometry Service。

**本机环境为 Ansys 2024R1（v241），使用 Windows 本地服务模式**：
确保安装了 SpaceClaim（或 Discovery），直接运行：

```bash
python main.py
```

代码会通过 `launch_modeler(product_version="241")` 启动本机 24.1 服务。

> 注意：
> - Linux Docker Core Service 仅支持 2025R2+，2024R1 不要使用该方案。
> - 仅可使用在 24.1 上存在的 API，24R2/25R1 后新增的功能会报错（AGENTS.md 中已约束 Codex）。

连接远程已有服务：

```bash
python main.py --host <服务器IP> --port 50051
```

## 3. 接入 Codex

1. 在 GitHub 新建空仓库 `geometry-automation`（不要勾选自动生成 README）。
2. 把本目录所有文件上传到仓库根目录（网页 Add file → Upload files），Commit changes。
3. 在 ChatGPT → Codex → 新建云端环境 → 选择 `geometry-automation` 仓库 → 开始设置。
   Codex 会自动安装 `requirements.txt` 并读取 `AGENTS.md`。

> 注意：Codex 云端环境没有 Geometry Service 和许可证，只负责写代码；
> 实际生成几何要在装了 Ansys 的机器上运行，或让脚本通过 `--host` 连接可达的 Geometry 服务。
