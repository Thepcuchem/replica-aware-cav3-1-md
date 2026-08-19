#!/usr/bin/env python3
"""Small contract tests for the ReplicaLab analysis registry."""

from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("replicalab_server", ROOT / "gui/server.py")
SERVER = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(SERVER)

class RegistryTests(unittest.TestCase):
    def test_registered_scripts_and_output_argument(self):
        self.assertGreaterEqual(len(SERVER.ANALYSES), 10)
        for analysis in SERVER.ANALYSES.values():
            script = SERVER.SOURCE / analysis["script"]
            self.assertTrue(script.is_file(), script)
            self.assertIn("--output-dir", script.read_text())

    def test_command_uses_isolated_output_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            analysis, command, values = SERVER.build_command(
                "distance-determinants", {"feature-dir": temporary, "top": 7},
                Path(temporary) / "output")
            self.assertEqual(analysis["name"], "Reproducible distance determinants")
            self.assertEqual(values["top"], 7)
            self.assertEqual(command[-2:], ["--output-dir", str(Path(temporary) / "output")])

    def test_missing_path_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Directory not found"):
            SERVER.build_command("distance-validation",
                {"feature-dir": "/path/that/does/not/exist", "components": 10, "seed": 1},
                Path("output"))

if __name__ == "__main__":
    unittest.main()
