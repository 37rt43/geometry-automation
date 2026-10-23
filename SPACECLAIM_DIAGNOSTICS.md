# SpaceClaim 24.1 脚本诊断

运行方式不变：

```bat
git pull origin main
python pump_impeller.py
```

毛坯和叶片两个阶段都通过 `spaceclaim_diagnostics.run_checked_script` 执行。
失败时控制台打印完整诊断 JSON，并保留实际执行的脚本及响应文件：

```text
work/spaceclaim_diagnostics/<时间_唯一编号>/meridional_blank.py
work/spaceclaim_diagnostics/<时间_唯一编号>/meridional_blank_response.json
work/spaceclaim_diagnostics/<时间_唯一编号>/impeller_blades.py
work/spaceclaim_diagnostics/<时间_唯一编号>/impeller_blades_response.json
```

每次调用使用独立目录；文件不再随临时目录删除，不加入 Git。
诊断包含实际收到的原始 protobuf 响应、服务端返回值、Python 堆栈、
可取得的 .NET 异常全文及客户端异常堆栈。叶片阶段还记录实际加载的
内核类型/程序集，而不是假设产品版本 241 对应某个 CLR 命名空间。

`stage` 区分轮廓构造、旋转、31 个放样截面、平滑放样、封盖缝合、
进口眼保护、复制旋转、各对叶片相交检查、布尔减和连通性检查。
`phase` 区分脚本执行与客户端模型刷新。不自动重试失败的建模操作。
若服务在执行包装器之前就失败（例如脚本引擎初始化失败），可能仍只收到
通用消息；此时记录收到的原始响应，不虚构不存在的详细堆栈。

## API 核对结果

- PyAnsys `Modeler.run_script_file` 文档注明支持 DMS API 24.1 及以后版本。
- `import_design` 是客户端参数，不在 `RunScriptFileRequest` 中。
  SDK 在服务端执行成功后才根据该参数决定是否调用 `read_existing_design()`。
  现使用 `import_design=False`，检查脚本结果后显式刷新，以区分两类失败。
- PyAnsys 的 `revolve_sketch`、`create_body_from_loft_profile` 和 `rotate`
  标注最低要求 24R2，不能将它们用作本项目 24.1 的直接调用替代方案。
- 原生 SpaceClaim 的 `Body.LoftProfiles`、`CreatePlanarBody`、`Fuse`、
  `Copy`、`Transform`、布尔运算及 `PieceCount` 可在早于
  24.1 的原生 API 示例中找到；这与 PyAnsys 对相应 gRPC 接口的版本限制不同。
- 已删除叶片脚本里未经本机确认的 `SpaceClaim.Api.V241.Modeler` 导入，改从
  当前毛坯的 `Shape` 取得实际加载的内核类型。产品版本和原生程序集命名
  不能直接等同。各重载在用户安装中的绑定仍须以实际运行结果确认。
- 本机已确认原生 `ExtrudeProfile` 不接受原先的三个参数。入口眼保护圆柱
  现改在本机 Python 中用 24.1 支持的 `Sketch.circle` + `extrude_sketch`
  创建，不再调用该原生重载；原生叶片裁切使用其形状副本，完成后删除工具体。
- 入口眼布尔减曾在相切端面处报 non-manifold。现在保留原有 31 个设计截面，
  额外增加一个完全落在 R27.5 保护圆柱内的牺牲入口截面，再按精确 R27.5
  裁回。延伸量至少 1 mm，并根据扭曲/厚度自动加深，使端面具备 0.05 mm
  内缩余量；这些参数集中于 `Parameters`。不缩小进口眼，也不忽略布尔异常。
- 叶片干涉检查改用双向差集体积损失：`V(A)-V(A-B)` 必须与
  `V(B)-V(B-A)` 一致。先用已知相离和半高重叠的圆柱校准该内核的差集行为，
  不再仅凭 `Intersect` 后接收体仍非空就判断干涉。每对叶片记录实测体积和
  重叠量；复制后逐顶点验证 72° 布置。真实重叠或校准/一致性失败仍阻止导出。
  叶片根部与盖板的 0.75 mm 超切不变，本模型仍只输出流体域，不建立金属轮毂。

离线验证：

```bat
python -m unittest test_spaceclaim_diagnostics -v
python -m py_compile pump_impeller.py pump_meridional.py spaceclaim_diagnostics.py
```

云端没有 Windows SpaceClaim/许可证；离线检查不能证明本机放样和布尔运算成功。
