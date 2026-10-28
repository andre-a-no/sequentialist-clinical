# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Novikov Laboratories LLC
# Commercial licenses for use outside the AGPL: see COMMERCIAL.md
"""ADaM: Dataset-JSON, NDJSON and CSV give the same records; the CDISC pilot analysed against lifelines."""
import csv
import json
import tempfile
import unittest
from pathlib import Path

from sequentialist.clinical import Trial
from sequentialist.clinical.adam import Dataset, date_seconds, outcome_records, read_dataset, time_to_event_records
from sequentialist.core import Engine, MemoryStorage

from .helpers import DATA, references, spec

PILOT = DATA / "cdiscpilot01"
KEY = "a key of the data owner, 16+ bytes"


def pilot_spec():
    s = spec("time-to-event", "greater", design={"analyses": [1.0], "maximum": 300, "power": None})
    return s.__class__({**s.data, "arms": {"treatment": "Xanomeline High Dose", "control": "Placebo"}})


class PilotTest(unittest.TestCase):
    """Time to the first dermatologic event in the CDISC pilot (cdiscpilot01), high dose against placebo."""

    def records(self, folder=PILOT, suffix=".json"):
        adsl, adtte = read_dataset(folder / f"adsl{suffix}"), read_dataset(folder / f"adtte{suffix}")
        return time_to_event_records(adsl, adtte, paramcd="TTDE", key=KEY, cut="cut-1")

    def test_log_rank_on_the_pilot(self):
        ref = references()["pilot"]
        trial = Trial.create(Engine(MemoryStorage()), "pilot", pilot_spec())
        trial.ingest(self.records())
        result = trial.analyse("2016-01-01T00:00:00Z", actor="s", reason="r", final=True)
        self.assertEqual(result["patients"], ref["patients"])
        self.assertEqual(result["events"], ref["events"])
        self.assertEqual(result["left_out"]["other_arm"], 84)  # the low dose
        # lifelines uses the hypergeometric variance for tied event times, sequential-analysis takes tied events
        # one at a time; the pilot has ties (days), so the two agree closely, not exactly
        self.assertAlmostEqual(result["z"] ** 2 / ref["chi2"], 1.0, delta=0.005)

    def test_identifiers_are_pseudonymised(self):
        records = self.records()
        text = json.dumps([r.to_json() for r in records])
        self.assertNotIn("01-701-1015", text)
        self.assertTrue(all(len(r.data["patient"]) == 64 for r in records))  # HMAC-SHA-256, hex

    def test_the_same_cut_twice_changes_nothing(self):
        trial = Trial.create(Engine(MemoryStorage()), "pilot", pilot_spec())
        first = trial.ingest(self.records())
        again = trial.ingest(self.records())
        self.assertEqual((first.accepted, again.accepted, again.duplicates), (508, 0, 508))

    def test_csv_and_ndjson_give_the_same_records(self):
        expected = [r.to_json() for r in self.records()]
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            for name in ("adsl", "adtte"):
                meta = json.loads((PILOT / f"{name}.json").read_text(encoding="utf-8"))
                rows = meta.pop("rows")
                with (folder / f"{name}.csv").open("w", newline="", encoding="utf-8") as f:
                    writer = csv.writer(f)
                    writer.writerow([c["name"] for c in meta["columns"]])
                    writer.writerows(rows)
                lines = [json.dumps(meta)] + [json.dumps(row) for row in rows]
                (folder / f"{name}.ndjson").write_text("\n".join(lines) + "\n", encoding="utf-8")
            self.assertEqual([r.to_json() for r in self.records(folder, ".csv")], expected)
            self.assertEqual([r.to_json() for r in self.records(folder, ".ndjson")], expected)

    def test_other_versions_and_formats_are_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            old = Path(tmp) / "old.json"
            old.write_text(json.dumps({"datasetJSONVersion": "1.0.0", "columns": [], "rows": []}), encoding="utf-8")
            with self.assertRaises(ValueError):
                read_dataset(old)
            with self.assertRaises(ValueError):
                read_dataset(Path(tmp) / "adsl.xpt")


class OutcomeTest(unittest.TestCase):
    def test_outcomes_from_a_bds_dataset(self):
        adsl = Dataset("ADSL", ["USUBJID", "TRT01P", "ITTFL", "TRTSDT", "AGEGR1"], [
            {"USUBJID": "S-1", "TRT01P": "TRT", "ITTFL": "Y", "TRTSDT": "2026-01-05", "AGEGR1": "<65"},
            {"USUBJID": "S-2", "TRT01P": "CTL", "ITTFL": "Y", "TRTSDT": 24112, "AGEGR1": ">=65"},  # a SAS date
            {"USUBJID": "S-3", "TRT01P": "CTL", "ITTFL": "N", "TRTSDT": "2026-01-07", "AGEGR1": "<65"},
        ])
        adrs = Dataset("ADRS", ["USUBJID", "PARAMCD", "AVAL", "ADT"], [
            {"USUBJID": "S-1", "PARAMCD": "RESP", "AVAL": 1, "ADT": "2026-03-01"},
            {"USUBJID": "S-2", "PARAMCD": "RESP", "AVAL": "0", "ADT": "2026-03-02"},
            {"USUBJID": "S-3", "PARAMCD": "RESP", "AVAL": 1, "ADT": "2026-03-03"},
            {"USUBJID": "S-1", "PARAMCD": "OTHER", "AVAL": 7, "ADT": "2026-03-01"},
        ])
        records = outcome_records(adsl, adrs, paramcd="RESP", key=KEY, cut="c", stratum="AGEGR1")
        kinds = [(r.source, r.seq, r.data["type"]) for r in records]
        self.assertEqual(kinds, [("c/ADSL", 0, "patient"), ("c/ADSL", 1, "patient"), ("c/ADRS", 0, "outcome"),
                                 ("c/ADRS", 1, "outcome")])
        self.assertEqual(records[1].data["entry"], date_seconds("2026-01-06"))
        self.assertEqual(records[0].data["stratum"], "<65")
        self.assertEqual([r.data["value"] for r in records[2:]], [1.0, 0.0])

    def test_dates(self):
        self.assertEqual(date_seconds(0), -315619200.0)  # 1960-01-01, the SAS epoch
        self.assertEqual(date_seconds("2014-01-02"), 1388620800.0)
        self.assertEqual(date_seconds("2014-01-02T12:00:00"), 1388620800.0 + 43200)
        with self.assertRaises(ValueError):
            date_seconds("02.01.2014")


if __name__ == "__main__":
    unittest.main()
