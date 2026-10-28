# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Novikov Laboratories LLC
# Commercial licenses for use outside the AGPL: see COMMERCIAL.md
"""Trials: statistics at a cut-off against external references, commands, recovery, blinding of the open report."""
import json
import math
import tempfile
import unittest
from pathlib import Path

from sequential_analysis.group_sequential import group_sequential_design

from sequentialist.clinical import Trial, event, outcome, patient
from sequentialist.core import (
    Engine,
    FileStorage,
    InvalidRecord,
    MemoryStorage,
    SQLiteStorage,
    TaskError,
    verify_audit,
)

from .helpers import DAY, T0, outcome_records, references, spec, survival_records


def new_trial(s, storage=None, task_id="trial"):
    return Trial.create(Engine(storage or MemoryStorage(), clock="event"), task_id, s)


def big(endpoint="time-to-event", alternative="less", **design):
    """A design whose maximum is far away, so that one analysis is an interim one."""
    return spec(endpoint, alternative, design={"analyses": [0.5, 1.0], "maximum": 5000, "power": None, **design})


class ReferenceTest(unittest.TestCase):
    """The statistic of an analysis is the one external implementations compute on the same data."""

    def test_log_rank_against_lifelines(self):
        ref = references()["survival"]
        trial = new_trial(big())
        trial.ingest(survival_records(ref["rows"]))
        result = trial.analyse(T0 + ref["cutoff"] * DAY, actor="s", reason="r")
        self.assertEqual(result["patients"], ref["patients"])
        self.assertEqual(result["events"], ref["events"])
        self.assertEqual(result["information"], sum(ref["events"]))
        self.assertAlmostEqual(result["z"] ** 2, ref["chi2"], delta=1e-9 * ref["chi2"])
        self.assertGreater(result["z"], 0)  # fewer events on treatment: positive for alternative "less"

    def test_two_proportions_against_scipy(self):
        ref = references()["binary"]
        trial = new_trial(big("binary", "greater"))
        trial.ingest(outcome_records(ref["rows"]))
        result = trial.analyse(T0 + 1000 * DAY, actor="s", reason="r")
        self.assertAlmostEqual(result["z"] ** 2, ref["chi2"], delta=1e-9 * ref["chi2"])
        treated = [r for r in ref["rows"] if r["arm"] == 1]
        self.assertEqual(result["successes"][0], sum(r["value"] for r in treated))

    def test_two_means_against_scipy(self):
        ref = references()["continuous"]
        trial = new_trial(big("continuous", "greater", assumptions={}))
        trial.ingest(outcome_records(ref["rows"]))
        result = trial.analyse(T0 + 1000 * DAY, actor="s", reason="r")
        self.assertAlmostEqual(result["z"], ref["t"], delta=1e-9 * abs(ref["t"]))

    def test_two_means_with_a_known_sigma(self):
        rows = references()["continuous"]["rows"]
        trial = new_trial(big("continuous", "greater", assumptions={}).__class__(
            {**big("continuous", "greater", assumptions={}).data,
             "endpoint": {"type": "continuous", "alternative": "greater", "sigma": 1.3, "strata": []}}))
        trial.ingest(outcome_records(rows))
        result = trial.analyse(T0 + 1000 * DAY, actor="s", reason="r")
        t = [r["value"] for r in rows if r["arm"] == 1]
        c = [r["value"] for r in rows if r["arm"] == 0]
        z = (math.fsum(t) / len(t) - math.fsum(c) / len(c)) / (1.3 * math.sqrt(1 / len(t) + 1 / len(c)))
        self.assertAlmostEqual(result["z"], z, places=10)

    def test_boundaries_at_the_information_reached(self):
        ref = references()["survival"]
        s = spec(design={"analyses": [0.5, 1.0], "maximum": 200, "power": None})
        trial = new_trial(s)
        trial.ingest(survival_records(ref["rows"]))
        first = trial.analyse(T0 + 400 * DAY, actor="s", reason="r")
        second = trial.analyse(T0 + 600 * DAY, actor="s", reason="r")
        fractions = [first["information"] / 200, second["information"] / 200]
        realised = group_sequential_design(0.025, [0.5, 1.0]).at(fractions)
        # the second analysis is an interim one too: its boundary spends alpha(t2), not all that is left
        interim = group_sequential_design(0.025, [*fractions, 1.0]).bounds
        self.assertAlmostEqual(first["bound"], realised.bounds[0], places=9)
        self.assertAlmostEqual(second["bound"], interim[1], places=6)


class TrialTest(unittest.TestCase):
    def test_late_data_change_only_later_analyses(self):
        trial = new_trial(big())
        for i in range(40):
            arm = "TRT" if i % 2 else "CTL"
            trial.ingest([patient("cut-1", 2 * i, f"p{i}", arm, T0 + i * DAY),
                          event("cut-1", 2 * i + 1, f"p{i}", T0 + (i + 30 + (i % 7) * 11) * DAY, i % 3 == 0)])
        first = trial.analyse(T0 + 60 * DAY, actor="s", reason="r")
        # an event that happened before the first cut-off, reported later, and a correction of an arm
        trial.ingest([event("cut-2", 0, "p1", T0 + 40 * DAY, True), patient("cut-2", 1, "p2", "TRT", T0 + 2 * DAY)])
        second = trial.analyse(T0 + 90 * DAY, actor="s", reason="r")
        self.assertEqual(trial.analyses[0], first)
        self.assertEqual(trial.handler.counts["updates"], 1)
        self.assertGreaterEqual(sum(second["events"]), sum(first["events"]) + 1)

    def test_censoring_at_the_cut_off(self):
        trial = new_trial(big())
        trial.ingest([patient("c", 0, "a", "TRT", T0), event("c", 1, "a", T0 + 100 * DAY, True),
                      patient("c", 2, "b", "CTL", T0), event("c", 3, "b", T0 + 10 * DAY, True),
                      patient("c", 4, "c", "CTL", T0 + 50 * DAY), patient("c", 5, "d", "XYZ", T0)])
        result = trial.analyse(T0 + 20 * DAY, actor="s", reason="r")
        self.assertEqual(result["events"], [0, 1])  # a's event comes after the cut-off: censored at it
        self.assertEqual(result["left_out"], {"other_arm": 1, "incomplete": 0, "no_follow_up": 0,
                                              "inconsistent": 0})  # c enters after the cut-off

    def test_invalid_batches_are_refused_whole(self):
        trial = new_trial(big())
        with self.assertRaises(InvalidRecord):
            trial.ingest([patient("c", 0, "a", "TRT", T0), patient("c", 1, "Ivan Petrov", "TRT", T0)])
        self.assertEqual(trial.handler.patients, {})
        self.assertEqual(trial.task.counters["records"], 0)

    def test_commands_and_roles(self):
        s = big()
        s = s.__class__({**s.data, "roles": {"unblinded": ["stat"], "dsmb": ["chair"], "sponsor": ["sp"]}})
        trial = new_trial(s)
        trial.ingest(survival_records(references()["survival"]["rows"]))
        before = len(trial.task.audit())
        with self.assertRaises(TaskError):
            trial.analyse(T0 + 500 * DAY, actor="chair", reason="r")
        with self.assertRaises(TaskError):
            trial.recommend(1, "continue", actor="chair", reason="r")  # no analysis yet
        self.assertEqual(len(trial.task.audit()), before)  # refused before the journal
        trial.analyse(T0 + 500 * DAY, actor="stat", reason="first")
        with self.assertRaises(TaskError):
            trial.analyse(T0 + 400 * DAY, actor="stat", reason="before the previous cut-off")
        with self.assertRaises(TaskError):
            trial.analyse(T0 + 500 * DAY, actor="stat", reason="no new information")
        trial.recommend(1, "continue", actor="chair", reason="fine")
        with self.assertRaises(TaskError):
            trial.decide(1, True, actor="chair", reason="not the sponsor")
        trial.decide(1, True, actor="sp", reason="agreed")
        with self.assertRaises(TaskError):
            trial.decide(1, False, actor="sp", reason="twice")
        trial.close_trial(actor="sp", reason="end")
        trial.ingest([patient("late", 0, "z", "TRT", T0)])
        self.assertEqual(trial.handler.counts["after_close"], 1)
        with self.assertRaises(TaskError):
            trial.analyse(T0 + 900 * DAY, actor="stat", reason="closed")
        kinds = [(e.kind, e.decision) for e in trial.task.audit()]
        self.assertEqual(kinds[-4:], [("analysis", "continue"), ("action", "recommendation"),
                                      ("action", "accepted"), ("action", "closed")])

    def test_a_refused_command_from_the_task_is_recorded(self):
        trial = new_trial(big())
        trial.task.command("recommend", {"analysis": 3, "action": "continue"}, actor="x", reason="y")
        last = trial.task.audit()[-1]
        self.assertEqual((last.kind, last.decision), ("action", "refused"))

    def test_the_test_ends_at_a_boundary(self):
        ref = references()["survival"]
        s = spec(design={"analyses": [0.5, 1.0], "maximum": 120, "power": None})
        trial = new_trial(s)
        trial.ingest(survival_records(ref["rows"]))
        result = trial.analyse(T0 + 700 * DAY, actor="s", reason="r")  # 122 events: beyond the maximum, final
        self.assertTrue(result["final"])
        self.assertEqual(result["decision"], "reject_h0" if result["z"] >= result["bound"] else "accept_h0")
        self.assertEqual(trial.handler.finished, result["decision"])
        with self.assertRaises(TaskError):
            trial.analyse(T0 + 800 * DAY, actor="s", reason="r")

    def test_anytime_design(self):
        ref = references()["survival"]
        s = spec(design={"type": "anytime", "analyses": None, "power": None, "assumptions": {}})
        trial = new_trial(s)
        trial.ingest(survival_records(ref["rows"]))
        looks = [trial.analyse(T0 + d * DAY, actor="s", reason="r") for d in (300, 400, 500)]
        self.assertTrue(all("log_e_value" in a for a in looks))
        self.assertEqual(trial.handler.analyses[-1]["decision"], looks[-1]["decision"])
        closed = trial.report("closed")
        guarantees = {f.guarantee for f in closed.findings()}
        self.assertIn("exact", guarantees)

    def test_stratified_log_rank(self):
        s = spec(more_endpoint={"strata": ["young", "old"]}, design={"maximum": 5000, "power": None})
        trial = new_trial(s)
        records = []
        for i, r in enumerate(references()["survival"]["rows"]):
            arm = "TRT" if r["arm"] == 1 else "CTL"
            records.append(patient("c", 2 * i, r["id"], arm, T0 + r["entry"] * DAY, "young" if i % 3 else "old"))
            records.append(event("c", 2 * i + 1, r["id"], T0 + r["end"] * DAY, r["event"]))
        trial.ingest(records)
        result = trial.analyse(T0 + 700 * DAY, actor="s", reason="r")
        self.assertEqual(result["events"], references()["survival"]["events"])
        self.assertTrue(math.isfinite(result["z"]))


def snapshot(trial):
    """What must be identical between an uninterrupted and a recovered trial (task ids aside)."""
    audit = []
    for e in trial.task.audit():
        entry = e.to_json()
        for key in ("task_id", "hash", "prev"):
            entry.pop(key)
        entry["evidence"] = [{k: v for k, v in ev.items() if k != "run_id"} for ev in entry["evidence"]]
        audit.append(entry)
    return {"state": trial.handler.state(), "audit": audit}


class RecoveryTest(unittest.TestCase):
    """A trial that crashes and is reopened is the uninterrupted trial, bit for bit, on every storage."""

    def storages(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        made = [MemoryStorage(), FileStorage(Path(tmp.name) / "files", sync="flush"),
                SQLiteStorage(Path(tmp.name) / "trials.db", sync="flush")]
        for storage in made:
            self.addCleanup(storage.close)
        return made

    def run_trial(self, trial, crash_after=None):
        rows = references()["survival"]["rows"]
        records = survival_records(rows)
        engine = trial.task.engine
        steps = [("ingest", records[:200]), ("analyse", 300), ("ingest", records[200:]), ("analyse", 500),
                 ("recommend", 2), ("analyse", 700)]
        for i, (what, arg) in enumerate(steps):
            if what == "ingest":
                trial.ingest(arg)
            elif what == "analyse":
                trial.analyse(T0 + arg * DAY, actor="s", reason="r")
            else:
                trial.recommend(arg, "continue", actor="c", reason="r")
            if crash_after == i:
                trial.task.closed = True  # the process dies: nothing more is written, the lock goes with it
                trial.task._lock.release()
                trial = Trial.open(engine, trial.task.task_id)
        return trial

    def test_recovery(self):
        s = spec(design={"analyses": [0.5, 1.0], "maximum": 400, "power": None})
        for storage in self.storages():
            straight = self.run_trial(new_trial(s, storage, "a"))
            expected = snapshot(straight)
            straight.close()
            for crash_after in range(6):
                with self.subTest(storage=type(storage).__name__, crash_after=crash_after):
                    recovered = self.run_trial(new_trial(s, storage, f"b{crash_after}"), crash_after)
                    got = snapshot(recovered)
                    self.assertEqual(got["state"], expected["state"])
                    self.assertEqual(got["audit"], expected["audit"])
                    self.assertEqual(verify_audit([x.to_json() for x in recovered.task.audit()]),
                                     len(got["audit"]))
                    recovered.close()


class BlindingTest(unittest.TestCase):
    def test_the_open_report_does_not_reveal_the_comparison(self):
        s = spec(design={"analyses": [0.5, 1.0], "maximum": 400, "power": None})
        s = s.__class__({**s.data, "arms": {"treatment": "ARM-TREATED", "control": "ARM-CONTROL"}})
        trial = new_trial(s)
        trial.ingest(survival_records(references()["survival"]["rows"], arms=("ARM-TREATED", "ARM-CONTROL")))
        result = trial.analyse(T0 + 500 * DAY, actor="s", reason="r")
        trial.recommend(1, "continue", actor="c", reason="r")
        opened = json.dumps(trial.report("open").to_json())
        closed = json.dumps(trial.report("closed").to_json())
        for secret in ("ARM-TREATED", "ARM-CONTROL", repr(result["z"]), repr(result["bound"]), "reject_h0",
                       "accept_h0", "Z ="):
            self.assertNotIn(secret, opened)
        self.assertIn(repr(result["z"]), closed)
        self.assertIn("ARM-TREATED", closed)

    def test_small_cells(self):
        s = spec(design={"maximum": 400, "power": None})
        trial = new_trial(s)
        trial.ingest([patient("c", i, f"p{i}", "TRT", T0) for i in range(3)])
        table = trial.report("open").sections[0].tables[0]
        self.assertEqual(table.rows[0][1], "<5")

    def test_protocol(self):
        trial = new_trial(big())
        trial.ingest(survival_records(references()["survival"]["rows"]))
        trial.analyse(T0 + 500 * DAY, actor="s", reason="first look")
        trial.recommend(1, "continue", actor="c", reason="r")
        text = trial.protocol(1).to_markdown()
        self.assertIn("group_sequential_logrank/1", text)
        self.assertIn("group_sequential_monitor/1", text)
        self.assertIn("first look", text)
        with self.assertRaises(ValueError):
            trial.protocol(2)


class DeterminismTest(unittest.TestCase):
    def test_the_same_inputs_give_the_same_audit_log(self):
        hashes = []
        for _ in range(2):
            trial = new_trial(big("binary", "greater"))
            trial.ingest(outcome_records(references()["binary"]["rows"]))
            trial.analyse(T0 + 200 * DAY, actor="s", reason="r")
            trial.ingest([outcome("late", 0, "b0000", 1, T0 + 5 * DAY)])
            trial.analyse(T0 + 400 * DAY, actor="s", reason="r")
            hashes.append([e.hash for e in trial.task.audit()])
        self.assertEqual(hashes[0], hashes[1])


if __name__ == "__main__":
    unittest.main()
