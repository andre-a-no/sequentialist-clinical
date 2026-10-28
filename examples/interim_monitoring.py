# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Novikov Laboratories LLC
# Commercial licenses for use outside the AGPL: see COMMERCIAL.md
"""A time-to-event trial monitored at three data cuts, on synthetic data.

The design is examples/survival.designed.json: one-sided alpha 0.025,
O'Brien-Fleming alpha spending, non-binding futility boundaries, analyses
at 50 %, 75 % and 100 % of the events needed for 90 % power at a hazard
ratio of 0.7.  Patients are synthetic (a true hazard ratio of 0.75); at
each data cut the trial receives what is known by then: who entered, who
had the event, and for everyone else the last contact at the cut.  The
unblinded statistician runs the analysis, the committee recommends, the
sponsor decides.  The trial is stored in SQLite and opened anew for every
cut, as a service would.
"""
import math
import random
import tempfile
from pathlib import Path

from sequentialist.clinical import Trial, TrialSpec, event, patient
from sequentialist.core import Engine, SQLiteStorage, load_spec, verify_audit
from sequentialist.core.schema import format_timestamp

HERE = Path(__file__).resolve().parent
DAY = 86400.0
START = 1768176000.0  # 2026-01-12T00:00:00Z
SEED = 6  # with this seed the trial stops at the second analysis; other seeds stop later, or for futility


RECOMMENDATION = {
    "continue": ("continue", "no boundary crossed, no safety concern"),
    "reject_h0": ("stop-efficacy", "the efficacy boundary is crossed"),
    "accept_h0": ("stop-futility", "the futility boundary is crossed, or the final analysis did not reject"),
}


def synthetic_patients(rng, n, hazard_ratio):
    people = []
    for i in range(n):
        treated = rng.random() < 0.5
        entry = START + i * DAY / 1.4
        rate = math.log(2) / 365 * (hazard_ratio if treated else 1.0)
        people.append({"id": f"pt-{i:05d}", "arm": "A" if treated else "B", "entry": entry,
                       "event": entry + rng.expovariate(rate) * DAY})
    return people


def data_cut(people, cut, label):
    """The records of a data cut: entries by then, events by then, the last contact at the cut for the rest."""
    records = []
    for p in people:
        if p["entry"] > cut:
            continue
        records.append(patient(label, len(records), p["id"], p["arm"], p["entry"]))
        happened = p["event"] <= cut
        records.append(event(label, len(records), p["id"], p["event"] if happened else cut, happened))
    return records


def main():
    spec = load_spec(HERE / "survival.designed.json", TrialSpec)
    designed = spec.data["designed"]
    people = synthetic_patients(random.Random(SEED), math.ceil(1.5 * designed["maximum"]), 0.75)
    events = sorted(p["event"] for p in people)
    cuts = [math.ceil(events[math.ceil(t * designed["maximum"]) - 1]) for t in designed["fractions"]]
    with tempfile.TemporaryDirectory() as tmp:
        engine = Engine(SQLiteStorage(Path(tmp) / "trials.db"))
        Trial.create(engine, "rx-301", spec).close()
        for k, cut in enumerate(cuts, start=1):
            with Trial.open(engine, "rx-301") as trial:
                trial.ingest(data_cut(people, cut, f"cut-{k}"))
                result = trial.analyse(cut, actor="statistician-1", reason=f"planned analysis {k}",
                                       final=k == len(cuts))
                print(f"analysis {k} at {format_timestamp(cut)}: {result['information']} events, "
                      f"Z = {result['z']:.3f}, efficacy boundary {result['bound']:.3f}: {result['decision']}")
                action, why = RECOMMENDATION[result["decision"]]
                trial.recommend(k, action, actor="dsmb-chair", reason=why)
                trial.decide(k, True, actor="sponsor-1", reason="the committee's recommendation is followed")
                if trial.handler.finished is not None:
                    trial.unblind(actor="sponsor-1", reason="the test has decided")
                    trial.close_trial(actor="sponsor-1", reason="the trial is stopped")
                    break
        with Trial.open(engine, "rx-301") as trial:
            print()
            print(trial.report("open").to_markdown())
            print(trial.report("closed").to_markdown())
            print(trial.protocol(len(trial.analyses)).to_markdown())
            entries = [e.to_json() for e in trial.task.audit()]
            print(f"audit log: {verify_audit(entries, 'rx-301')} entries, chain verified")


if __name__ == "__main__":
    main()
