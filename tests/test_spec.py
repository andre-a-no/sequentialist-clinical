# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Novikov Laboratories LLC
# Commercial licenses for use outside the AGPL: see COMMERCIAL.md
"""The spec: defaults, every error at once with its pointer, the records a trial takes."""
import unittest

from sequentialist.clinical import TrialSpec
from sequentialist.clinical.spec import record_problem
from sequentialist.core import SpecError, json_schema, parse_spec

from .helpers import spec


def problems(data):
    try:
        parse_spec(data, TrialSpec)
    except SpecError as error:
        return dict(error.errors)
    return {}


def base(**design):
    return {"kind": "clinical.trial/1", "name": "t", "endpoint": {"type": "time-to-event", "alternative": "less"},
            "arms": {"treatment": "A", "control": "B"},
            "design": {"alpha": 0.025, "analyses": [0.5, 1.0], "maximum": 100, **design}}


class SpecTest(unittest.TestCase):
    def test_group_sequential_is_the_default(self):
        s = parse_spec(base(), TrialSpec)
        self.assertTrue(s.group_sequential)
        self.assertEqual(s.design["spending"], "obrien-fleming")
        self.assertEqual(s.data["small_cells"], 5)
        self.assertEqual(s.data["roles"], {"unblinded": [], "dsmb": [], "sponsor": []})

    def test_relations_between_members(self):
        found = problems({**base(analyses=[0.6, 0.5], spending="power", futility="pocock", two_sided=True),
                          "arms": {"treatment": "A", "control": "A"}})
        self.assertIn("/design/analyses", found)
        self.assertIn("/design/rho", found)
        self.assertIn("/design/futility", found)
        self.assertIn("/design/power", found)
        self.assertIn("/arms", found)

    def test_what_the_maximum_needs(self):
        data = base()
        del data["design"]["maximum"]
        found = problems(data)
        self.assertIn("/design/power", found)
        self.assertIn("/design/assumptions/hazard_ratio", found)

    def test_anytime_designs(self):
        found = problems({**base(type="anytime", futility="pocock"),
                          "endpoint": {"type": "continuous", "alternative": "greater"}})
        self.assertIn("/endpoint/type", found)
        self.assertIn("/design/futility", found)
        self.assertEqual(problems(base(type="anytime", method="mixture")), {})
        self.assertIn("/design/method", problems(base(type="anytime", method="universal")))

    def test_endpoint_members(self):
        found = problems({**base(), "endpoint": {"type": "binary", "alternative": "greater", "sigma": 1.0,
                                                 "strata": ["a", "a"]}})
        self.assertIn("/endpoint/sigma", found)
        self.assertIn("/endpoint/strata", found)

    def test_schema_errors_at_once(self):
        found = problems({"kind": "clinical.trial/1", "name": "", "endpoint": {"type": "x", "alternative": "up"},
                          "arms": {}, "design": {"alpha": 0.7}, "extra": 1})
        for pointer in ("/name", "/endpoint/type", "/endpoint/alternative", "/arms/treatment", "/design/alpha",
                        "/extra"):
            self.assertIn(pointer, found)

    def test_json_schema(self):
        schema = json_schema(TrialSpec)
        self.assertEqual(schema["title"], "clinical.trial/1")
        self.assertIn("design", schema["properties"])


class RecordTest(unittest.TestCase):
    def test_records(self):
        s = spec("binary", "greater")
        good = {"type": "patient", "patient": "p-1", "arm": "TRT", "entry": "2026-01-01T00:00:00Z"}
        self.assertIsNone(record_problem(s, good))
        self.assertIsNone(record_problem(s, {"type": "outcome", "patient": "p-1", "value": 1, "time": 5}))
        self.assertIn("0 or 1", record_problem(s, {"type": "outcome", "patient": "p-1", "value": 2, "time": 5}))
        self.assertIn("event record", record_problem(s, {"type": "event", "patient": "p-1", "time": 5,
                                                         "event": True}))
        self.assertIn("no strata", record_problem(s, {**good, "stratum": "x"}))

    def test_identifying_data_is_refused(self):
        s = spec()
        self.assertIn("does not match", record_problem(s, {"type": "patient", "patient": "Иванов Иван",
                                                           "arm": "TRT", "entry": 0}))
        self.assertIn("unknown member", record_problem(s, {"type": "event", "patient": "p-1", "time": 0,
                                                           "event": True, "birth_date": "1950-01-01"}))

    def test_strata(self):
        s = spec(more_endpoint={"strata": ["<65", ">=65"]})
        self.assertIsNone(record_problem(s, {"type": "patient", "patient": "p", "arm": "TRT", "entry": 0,
                                             "stratum": "<65"}))
        self.assertIn("not one of", record_problem(s, {"type": "patient", "patient": "p", "arm": "TRT", "entry": 0,
                                                       "stratum": "old"}))


if __name__ == "__main__":
    unittest.main()
