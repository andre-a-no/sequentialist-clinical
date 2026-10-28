# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Novikov Laboratories LLC
# Commercial licenses for use outside the AGPL: see COMMERCIAL.md
"""The `sequentialist clinical` command."""
import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

from sequentialist.clinical import Trial
from sequentialist.clinical.cli import register
from sequentialist.core import Engine, FileStorage

from .helpers import DAY, T0, references, spec, survival_records


def run(*argv):
    import argparse

    parser = argparse.ArgumentParser(prog="sequentialist")
    register(parser.add_subparsers(dest="command"))
    args = parser.parse_args(list(argv))
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()):
        code = args.handler(args)
    return code, out.getvalue()


class CliTest(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = Path(tmp.name)

    def test_design_check_plan_schema(self):
        source = self.dir / "trial.json"
        data = spec().to_json()
        data["designed"] = None
        source.write_text(json.dumps(data), encoding="utf-8")
        designed = self.dir / "designed.json"
        self.assertEqual(run("clinical", "design", str(source), "-o", str(designed))[0], 0)
        self.assertEqual(run("clinical", "check", str(designed))[0], 0)
        edited = json.loads(designed.read_text(encoding="utf-8"))
        edited["designed"]["maximum"] += 1
        designed.write_text(json.dumps(edited), encoding="utf-8")
        code, out = run("clinical", "check", str(designed))
        self.assertEqual(code, 1)
        self.assertIn("maximum", out)
        _, out = run("clinical", "plan", str(source.with_name("designed.json")), "--json")
        self.assertEqual(json.loads(out)["format"], "sequentialist.report/1")
        _, out = run("clinical", "schema")
        self.assertEqual(json.loads(out)["title"], "clinical.trial/1")

    def test_simulate(self):
        path = self.dir / "designed.json"
        path.write_text(json.dumps(spec(design={"analyses": [0.5, 1.0], "maximum": 40, "power": None}).to_json()),
                        encoding="utf-8")
        _, out = run("clinical", "simulate", str(path), "--runs", "3", "--json")
        self.assertEqual([r["under"] for r in json.loads(out)], ["null", "alternative"])

    def test_report_and_protocol_from_the_storage(self):
        storage = FileStorage(self.dir / "trials")
        with Trial.create(Engine(storage), "t-1", spec(design={"maximum": 900, "power": None})) as trial:
            trial.ingest(survival_records(references()["survival"]["rows"]))
            trial.analyse(T0 + 500 * DAY, actor="s", reason="first")
        name = f"files:{self.dir / 'trials'}"
        _, out = run("clinical", "report", name, "t-1")
        self.assertIn("open part", out)
        _, out = run("clinical", "report", name, "t-1", "--part", "closed", "--json")
        self.assertEqual(json.loads(out)["task_id"], "t-1")
        _, out = run("clinical", "protocol", name, "t-1", "1")
        self.assertIn("protocol of analysis 1", out)


if __name__ == "__main__":
    unittest.main()
