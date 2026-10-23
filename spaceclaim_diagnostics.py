"""保留 SpaceClaim 脚本、原始 RPC 响应及服务端堆栈；不重试建模操作。"""

from datetime import datetime, timezone
import json
from pathlib import Path
import sys
import traceback
from uuid import uuid4


def instrument_script(source):
    """在服务端捕获异常并用 result 返回，避免服务只返回 Script run failed。"""
    return '''# Python Script, API Version = V241
import traceback
stage = "bootstrap"
result = {}
try:
''' + "\n".join("    " + line for line in source.splitlines()) + '''
except Exception as error:
    result = {"diagnostic_ok": "False", "stage": stage,
              "error_type": str(type(error)), "error_repr": repr(error),
              "traceback": traceback.format_exc()}
    if "api_context" in globals():
        result["api_context"] = str(api_context)
    try:
        result["dotnet_exception"] = error.clsException.ToString()
    except Exception:
        try:
            result["dotnet_exception"] = error.ToString()
        except Exception:
            pass
else:
    result = dict(result or {})
    result["diagnostic_ok"] = "True"
    result["stage"] = stage
    if "api_context" in globals():
        result["api_context"] = str(api_context)
'''


class _RecordingStub:
    """仅观察这一服务实例的 RunScriptFile 返回值，其他方法原样转发。"""

    def __init__(self, original, report):
        self.original = original
        self.report = report

    def __getattr__(self, name):
        return getattr(self.original, name)

    def RunScriptFile(self, *args, **kwargs):
        response = self.original.RunScriptFile(*args, **kwargs)
        self.report["response_received"] = True
        from google.protobuf.json_format import MessageToDict
        from google.protobuf.text_format import MessageToString

        self.report["raw_response"] = MessageToDict(response, preserving_proto_field_name=True)
        self.report["rpc_success"] = response.success
        self.report["raw_response_protobuf"] = MessageToString(response)
        return response


def run_checked_script(modeler, source, label, directory=None):
    """执行一次脚本，再单独刷新。失败打印完整已收到的响应并保留诊断文件。

    import_design 是 Modeler 客户端选项，不属于 RunScriptFileRequest。
    使用 False 将服务端执行错误与后续 read_existing_design 错误区分开。
    V0/V1 的 SDK 适配层目前均公开 stub；若未来结构改变，保留客户端异常
    和服务端 result，明确注明无法观察原始 protobuf，而不伪造响应。
    """
    if directory is None:
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S") + "_" + uuid4().hex[:8]
        directory = Path(__file__).resolve().parent / "work" / "spaceclaim_diagnostics" / stamp
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / (label + ".py")
    path.write_text(instrument_script(source), encoding="ascii")
    report_path = directory / (label + "_response.json")
    report = {"script": str(path), "phase": "script_execution", "response_received": False}
    service = modeler.client.services.commands_script
    original = getattr(service, "stub", None)
    if original is not None:
        service.stub = _RecordingStub(original, report)
    else:
        report["raw_capture_unavailable"] = "SDK service has no observable stub"
    try:
        values, _ = modeler.run_script_file(path, import_design=False)
        report["server_values"] = values
        if not values or values.get("diagnostic_ok") != "True":
            raise RuntimeError("SpaceClaim server script failed at stage: "
                               + str((values or {}).get("stage", "no diagnostic result")))
        report["phase"] = "design_refresh"
        design = modeler.read_existing_design()
        report["phase"] = "complete"
        return values, design
    except Exception as error:
        report["client_exception"] = repr(error)
        report["client_traceback"] = traceback.format_exc()
        # gRPC 异常可能在收到 protobuf 之前发生；保留实际可取得的传输诊断。
        for attribute in ("code", "details", "debug_error_string"):
            method = getattr(error, attribute, None)
            if callable(method):
                try:
                    report[attribute] = str(method())
                except Exception:
                    pass
        print(json.dumps(report, ensure_ascii=False, indent=2), file=sys.stderr)
        print(f"SpaceClaim 诊断文件：{report_path}", file=sys.stderr)
        raise RuntimeError(f"SpaceClaim 建模失败（{report['phase']}）；完整响应：{report_path}") from error
    finally:
        if original is not None:
            service.stub = original
        report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
