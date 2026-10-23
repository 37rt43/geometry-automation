"""无需后端的诊断回归测试：python -m unittest test_spaceclaim_diagnostics。"""

from contextlib import redirect_stderr
from io import StringIO
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from ansys.api.dbu.v0.dbuapplication_pb2 import RunScriptFileResponse

from spaceclaim_diagnostics import instrument_script, run_checked_script


class FakeModeler:
    def __init__(self, response, refresh_error=None):
        self.response = response
        self.refresh_error = refresh_error
        self.refresh_count = 0
        self.calls = 0
        self.original_stub = SimpleNamespace(RunScriptFile=lambda request: response)
        self.client = SimpleNamespace(services=SimpleNamespace(
            commands_script=SimpleNamespace(stub=self.original_stub)))

    def run_script_file(self, path, import_design):
        assert import_design is False
        self.calls += 1
        response = self.client.services.commands_script.stub.RunScriptFile(None)
        if not response.success:
            raise RuntimeError(response.message)
        return dict(response.values), None

    def read_existing_design(self):
        self.refresh_count += 1
        if self.refresh_error:
            raise self.refresh_error
        return "refreshed design"


class DiagnosticsTests(unittest.TestCase):
    def test_server_wrapper_preserves_traceback_and_stage(self):
        scope = {}
        exec(instrument_script('stage = "smooth_loft"\nraise AttributeError("missing API")'), scope)
        self.assertEqual(scope["result"]["diagnostic_ok"], "False")
        self.assertEqual(scope["result"]["stage"], "smooth_loft")
        self.assertIn("AttributeError: missing API", scope["result"]["traceback"])

    def test_wrapper_keeps_success_values(self):
        scope = {}
        exec(instrument_script('result = {"pieces": "1"}'), scope)
        self.assertEqual(scope["result"]["pieces"], "1")
        self.assertEqual(scope["result"]["diagnostic_ok"], "True")

    def test_raw_failure_response_is_printed_and_preserved_without_retry(self):
        modeler = FakeModeler(RunScriptFileResponse(
            success=False, message="Script run failed", values={"detail": "native import failed"}))
        with TemporaryDirectory() as directory, redirect_stderr(StringIO()) as output:
            with self.assertRaises(RuntimeError):
                run_checked_script(modeler, "pass", "test", directory)
            report = json.loads((Path(directory)/"test_response.json").read_text())
            self.assertIn("native import failed", output.getvalue())
            self.assertEqual(report["raw_response"]["values"]["detail"], "native import failed")
            self.assertFalse(report["rpc_success"])
            self.assertTrue((Path(directory)/"test.py").is_file())
        self.assertEqual(modeler.calls, 1)
        self.assertEqual(modeler.refresh_count, 0)
        self.assertIs(modeler.client.services.commands_script.stub, modeler.original_stub)

    def test_caught_server_error_stops_before_refresh(self):
        modeler = FakeModeler(RunScriptFileResponse(success=True, values={
            "diagnostic_ok": "False", "stage": "smooth_loft", "traceback": "full server stack"}))
        with TemporaryDirectory() as directory, redirect_stderr(StringIO()) as output:
            with self.assertRaises(RuntimeError):
                run_checked_script(modeler, "pass", "test", directory)
            self.assertIn("full server stack", output.getvalue())
        self.assertEqual(modeler.refresh_count, 0)

    def test_success_refreshes_separately_and_restores_stub(self):
        modeler = FakeModeler(RunScriptFileResponse(success=True, values={"diagnostic_ok": "True"}))
        with TemporaryDirectory() as directory:
            values, design = run_checked_script(modeler, "pass", "test", directory)
            self.assertEqual(design, "refreshed design")
            report = json.loads((Path(directory)/"test_response.json").read_text())
            self.assertEqual(report["phase"], "complete")
        self.assertEqual(modeler.refresh_count, 1)
        self.assertIs(modeler.client.services.commands_script.stub, modeler.original_stub)

    def test_refresh_failure_has_distinct_phase(self):
        modeler = FakeModeler(RunScriptFileResponse(success=True, values={"diagnostic_ok": "True"}),
                              refresh_error=ValueError("refresh failed"))
        with TemporaryDirectory() as directory, redirect_stderr(StringIO()):
            with self.assertRaises(RuntimeError):
                run_checked_script(modeler, "pass", "test", directory)
            report = json.loads((Path(directory)/"test_response.json").read_text())
            self.assertEqual(report["phase"], "design_refresh")
        self.assertIs(modeler.client.services.commands_script.stub, modeler.original_stub)


if __name__ == "__main__":
    unittest.main()
