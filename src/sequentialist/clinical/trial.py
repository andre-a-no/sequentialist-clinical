# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Novikov Laboratories LLC
# Commercial licenses for use outside the AGPL: see COMMERCIAL.md
"""A trial as a task of the sequentialist engine: facts about patients in, analyses at data cut-offs out.

The task keeps what is known about every patient: arm, stratum, entry, and
the outcome or the follow-up for the event.  A later record about the
same fact replaces the earlier one, so new follow-up, a new data cut and a
correction all arrive the same way, and the order in which facts are
reported does not matter.

An analysis is a command, ``analyse``, with a data cut-off.  The data as
of the cut-off are put in order and run through a fresh run of a
sequential-analysis procedure (the stream ``data``):

* binary and continuous endpoints: every patient with an outcome observed
  by the cut-off, in the order of observation;
* time to event: every patient who entered by the cut-off, on the time
  scale of the study (from entry), with an event if it happened by the
  cut-off and censored at the last contact or at the cut-off otherwise.
  Patients without any follow-up record are left out (and counted).

A **group sequential** design then gives the standardised statistic and
the information reached (patients with an outcome, or events) to
`GroupSequentialMonitor` (the stream ``monitor``), which computes the
boundary of this analysis by alpha spending at the information actually
reached (Lan and DeMets): analyses need not happen when planned.  The
analysis is final when the information reaches the maximum or when the
command says so; a crossed boundary ends the test (the trial itself goes
on until the sponsor closes it).  An **anytime** design reads the decision
of its anytime-valid test (``data``) at the cut-off.

The audit log holds every analysis with its numbers and the evidence of
both streams (procedure, guarantee, history digest: the exact sequence of
data the decision rests on), every recommendation of the data monitoring
committee and every decision of the sponsor, with actor and reason.
Reaching a boundary is recorded; it stops nothing by itself.
"""
from __future__ import annotations

import math
from typing import Any, Iterable, Mapping

from sequential_analysis.procedures import (
    AnytimeLogRankTest,
    GroupSequentialLogRank,
    GroupSequentialMonitor,
    GroupSequentialProportions,
    GroupSequentialTwoMeans,
    TwoProportionTest,
)

from sequentialist.core import (
    AuditEntry,
    Command,
    Context,
    Engine,
    Handler,
    IngestResult,
    InvalidRecord,
    Record,
    Stream,
    Task,
    TaskError,
)
from sequentialist.core.schema import parse_timestamp

from .design import spending_of
from .spec import ACTIONS, KIND, TrialSpec, normal_record, record_problem

#: The arms as sequential-analysis labels them.
TREATMENT, CONTROL = 1, 0
ENTER, EVENT, CENSOR = 0, 1, 2
#: An analysis no stream reaches: the ``data`` stream of a group sequential design only accumulates.
_NEVER_ANALYSED = 2**53
_NO_BOUNDARY = 40.0
COMMANDS = ("analyse", "recommend", "decide", "unblind", "close")
_ROLE = {"analyse": "unblinded", "recommend": "dsmb", "decide": "sponsor", "unblind": "sponsor", "close": "sponsor"}


def _data_stream(spec: TrialSpec) -> Stream:
    endpoint = spec.data["endpoint"]
    design = spec.design
    alternative = endpoint["alternative"]
    arms = (TREATMENT, CONTROL)
    if spec.group_sequential:
        common: dict[str, Any] = {"alternative": alternative, "arms": arms, "label": "data at the cut-off"}
        if spec.endpoint == "binary":
            return Stream(GroupSequentialProportions([_NEVER_ANALYSED], [_NO_BOUNDARY], **common))
        if spec.endpoint == "continuous":
            return Stream(GroupSequentialTwoMeans(endpoint["sigma"], [_NEVER_ANALYSED], [_NO_BOUNDARY], **common))
        return Stream(GroupSequentialLogRank([_NEVER_ANALYSED], [_NO_BOUNDARY], strata=max(1, len(spec.strata())),
                                             **common))
    side = "two-sided" if design["two_sided"] else alternative
    method = {} if design["method"] is None else {"method": design["method"]}
    if spec.endpoint == "binary":
        return Stream(TwoProportionTest(design["alpha"], alternative=side, arms=arms, **method))
    return Stream(AnytimeLogRankTest(design["alpha"], alternative=side, arms=arms,
                                     strata=max(1, len(spec.strata())), **method))


def _monitor(spec: TrialSpec) -> GroupSequentialMonitor:
    design, designed = spec.design, spec.data["designed"]
    futility = spending_of(design["futility"], design["futility_rho"])
    extra: dict[str, Any] = {}
    if futility is not None:
        extra = {"futility": futility, "power": designed["power"], "drift": designed["drift"]}
    return GroupSequentialMonitor(design["alpha"], float(designed["maximum"]),
                                  spending_of(design["spending"], design["rho"]), two_sided=design["two_sided"],
                                  final_input="final", **extra)


class TrialHandler(Handler):
    """The domain logic of a trial (see the module docstring); deterministic, as the engine requires."""

    commands = COMMANDS

    def __init__(self, spec: TrialSpec | Mapping[str, Any]):
        self.spec = spec if isinstance(spec, TrialSpec) else TrialSpec(dict(spec))
        if self.spec.group_sequential and self.spec.data["designed"] is None:
            raise ValueError("a group sequential trial runs from a designed spec: run `sequentialist clinical design`")
        self.patients: dict[str, dict[str, Any]] = {}
        self.counts = {"records": 0, "updates": 0, "rejected": 0, "after_close": 0}
        self.analyses: list[dict[str, Any]] = []
        self.recommendations: list[dict[str, Any]] = []
        self.decisions: list[dict[str, Any]] = []
        self.finished: str | None = None  # the decision of the test, once it has made one
        self.unblinded = False
        self.closed = False
        self.fed = False  # the data stream holds the data of an earlier analysis

    def streams(self) -> dict[str, Stream]:
        streams = {"data": _data_stream(self.spec)}
        if self.spec.group_sequential:
            streams["monitor"] = Stream(_monitor(self.spec))
        return streams

    # ------------------------------------------------------------------ state

    def state(self) -> dict[str, Any]:
        return {"patients": self.patients, "counts": self.counts, "analyses": self.analyses,
                "recommendations": self.recommendations, "decisions": self.decisions, "finished": self.finished,
                "unblinded": self.unblinded, "closed": self.closed, "fed": self.fed}

    def restore(self, state: Any) -> None:
        self.patients = state["patients"]
        self.counts = state["counts"]
        self.analyses, self.recommendations = state["analyses"], state["recommendations"]
        self.decisions, self.finished = state["decisions"], state["finished"]
        self.unblinded, self.closed, self.fed = state["unblinded"], state["closed"], state["fed"]

    # ------------------------------------------------------------------ records

    def on_record(self, record: Record, ctx: Context) -> None:
        self.counts["records"] += 1
        if self.closed:
            self.counts["after_close"] += 1
            return
        if record_problem(self.spec, record.data) is not None:  # Trial.ingest refuses these before the journal
            self.counts["rejected"] += 1
            return
        data, _ = normal_record(record.data)
        assert data is not None
        facts = self.patients.setdefault(data["patient"], {"arm": None, "stratum": None, "entry": None,
                                                           "outcome": None, "follow": None})
        kind = data["type"]
        if kind == "patient":
            new: dict[str, Any] = {"arm": data["arm"], "stratum": data["stratum"], "entry": data["entry"]}
        elif kind == "outcome":
            new = {"outcome": [data["value"], data["time"]]}
        else:
            new = {"follow": [data["time"], data["event"]]}
        if kind != "event" and any(facts[k] is not None and facts[k] != v for k, v in new.items()):
            self.counts["updates"] += 1  # a correction; new follow-up for an event is not one
        facts.update(new)

    # ------------------------------------------------------------------ the data at a cut-off

    def snapshot(self, cutoff: float) -> dict[str, Any]:
        """The inputs of the data stream as of `cutoff`, in order, and what was left out."""
        arms = self.spec.data["arms"]
        label = {arms["treatment"]: TREATMENT, arms["control"]: CONTROL}
        strata = {name: i for i, name in enumerate(self.spec.strata())}
        out: dict[str, Any] = {"inputs": [], "patients": [0, 0], "units": [0, 0], "events": [0, 0],
                               "successes": [0, 0], "other_arm": 0, "incomplete": 0, "no_follow_up": 0,
                               "inconsistent": 0}
        included = []
        for pid in sorted(self.patients):
            facts = self.patients[pid]
            if facts["arm"] is None or facts["entry"] is None:
                out["incomplete"] += 1
                continue
            if facts["entry"] > cutoff:
                continue
            if facts["arm"] not in label:
                out["other_arm"] += 1
                continue
            arm = label[facts["arm"]]
            out["patients"][1 - arm] += 1  # [treatment, control]
            included.append((pid, arm, facts))
        if self.spec.endpoint == "time-to-event":
            self._survival(included, cutoff, strata, out)
        else:
            units = []
            for pid, arm, facts in included:
                if facts["outcome"] is not None and facts["outcome"][1] <= cutoff:
                    value, time = facts["outcome"]
                    units.append((time, pid, arm, value))
            for _, _, arm, value in sorted(units):
                out["inputs"].append({"arm": arm, "y": value})
                out["units"][1 - arm] += 1
                if self.spec.endpoint == "binary" and value == 1:
                    out["successes"][1 - arm] += 1
        out["information"] = sum(out["events"]) if self.spec.endpoint == "time-to-event" else sum(out["units"])
        return out

    def _survival(self, included: list, cutoff: float, strata: dict[str, int], out: dict[str, Any]) -> None:
        ends = []
        entered = []
        for pid, arm, facts in included:
            if facts["follow"] is None:
                out["no_follow_up"] += 1
                continue
            time, happened = facts["follow"]
            if time < facts["entry"]:
                out["inconsistent"] += 1
                continue
            if time <= cutoff:
                kind, end = (EVENT if happened else CENSOR), time
            else:
                kind, end = CENSOR, cutoff
            extra = {"stratum": strata[facts["stratum"]]} if strata else {}  # an unstratified test takes no stratum
            entered.append({"arm": arm, "kind": ENTER, **extra})
            ends.append((end - facts["entry"], kind, pid, {"arm": arm, "kind": kind, **extra}))
            if kind == EVENT:
                out["events"][1 - arm] += 1
            out["units"][1 - arm] += 1
        # on the study's time scale everybody enters at 0; at equal times events come before censorings
        out["inputs"] = entered + [x for _, _, _, x in sorted(ends, key=lambda e: (e[0], e[1], e[2]))]

    # ------------------------------------------------------------------ commands

    def problem(self, command: Command) -> str | None:
        """Why the command cannot be carried out now (None when it can)."""
        role = self.spec.data["roles"][_ROLE[command.name]]
        if role and command.actor not in role:
            return f"{command.actor!r} is not among those who may {command.name} ({_ROLE[command.name]})"
        args = command.args
        if command.name == "analyse":
            given, final = args.get("cutoff"), args.get("final", False)
            if set(args) - {"cutoff", "final"}:
                return f"analyse takes cutoff and final, got {sorted(args)}"
            if isinstance(given, bool) or not isinstance(given, (int, float)) or not math.isfinite(given):
                return "analyse needs the data cut-off, in seconds since the epoch"
            cutoff = float(given)
            if type(final) is not bool:
                return "final must be true or false"
            if self.closed:
                return "the trial is closed"
            if self.finished is not None:
                return f"the test has decided already ({self.finished}); no further analysis"
            if self.analyses and cutoff < self.analyses[-1]["cutoff"]:
                return "the cut-off is before that of the previous analysis"
            if self.spec.group_sequential:
                information = self.snapshot(cutoff)["information"]
                previous = self.analyses[-1]["information"] if self.analyses else 0
                if information <= previous:
                    return (f"no new information since the previous analysis ({information} "
                            f"{self.spec.unit} against {previous})")
            return None
        if command.name == "recommend":
            if set(args) != {"analysis", "action"}:
                return "recommend takes analysis (its number) and action"
            if args["action"] not in ACTIONS:
                return f"action must be one of {list(ACTIONS)}"
            if type(args["analysis"]) is not int or not 1 <= args["analysis"] <= len(self.analyses):
                return f"there is no analysis {args['analysis']!r}"
            return None
        if command.name == "decide":
            if set(args) != {"recommendation", "accept"}:
                return "decide takes recommendation (its number) and accept"
            if type(args["accept"]) is not bool:
                return "accept must be true or false"
            n = args["recommendation"]
            if type(n) is not int or not 1 <= n <= len(self.recommendations):
                return f"there is no recommendation {n!r}"
            if any(d["recommendation"] == n for d in self.decisions):
                return f"recommendation {n} has been decided already"
            return None
        if args:
            return f"{command.name} takes no arguments"
        if command.name == "unblind" and self.unblinded:
            return "the trial is unblinded already"
        if command.name == "close" and self.closed:
            return "the trial is closed already"
        return None

    def on_command(self, command: Command, ctx: Context) -> None:
        if command.name not in COMMANDS:
            super().on_command(command, ctx)
        problem = self.problem(command)
        if problem is not None:  # Trial refuses these before the journal; a direct task.command lands here
            ctx.record("action", "refused", basis={"command": command.name, "args": command.args,
                                                   "problem": problem},
                       actor=command.actor, reason=command.reason)
            return
        getattr(self, "_" + command.name)(command, ctx)

    def _analyse(self, command: Command, ctx: Context) -> None:
        cutoff, final = command.args["cutoff"], command.args.get("final", False)
        snap = self.snapshot(cutoff)
        if self.fed:
            ctx.restart("data", f"data cut-off of analysis {len(self.analyses) + 1}")
        for x in snap["inputs"]:
            ctx.step("data", x)
        self.fed = True
        run = ctx.run("data")
        if run.counters["invalid_skipped"] or run.t != len(snap["inputs"]):
            # a run skips invalid inputs; an analysis on part of its data must not happen
            raise TaskError(f"analysis {len(self.analyses) + 1}: the data stream took {run.t} of "
                            f"{len(snap['inputs'])} inputs")
        counts = {"patients": snap["patients"], "units": snap["units"], "information": snap["information"],
                  "left_out": {k: snap[k] for k in ("other_arm", "incomplete", "no_follow_up", "inconsistent")}}
        if self.spec.endpoint == "time-to-event":
            counts["events"] = snap["events"]
        if self.spec.endpoint == "binary":
            counts["successes"] = snap["successes"]
        k = len(self.analyses) + 1
        if self.spec.group_sequential:
            run = ctx.run("data")
            z = float(run.procedure.statistic(run.state))  # type: ignore[attr-defined]
            event = ctx.step("monitor", {"z": z, "information": float(snap["information"]), "final": int(final)})
            assert event is not None
            out = ctx.last_output("monitor")
            decision = out["decision"]
            result = {"z": z, "fraction": out["fraction"], "bound": out["bound"], "futility": out["futility"],
                      "spent": out["spent"], "final": out["final"]}
            evidence = ["data", "monitor"]
        else:
            out = ctx.last_output("data") or {}
            decision = out.get("decision")
            result = {key: out.get(key) for key in
                      ("log_e_value", "p_value", "hazard_ratio_estimate", "lower", "upper", "rates", "difference")
                      if key in out}
            result["final"] = bool(final)
            evidence = ["data"]
        if decision is None and final:
            decision = "accept_h0"  # the last look of an anytime design that did not reject (the monitor decides)
        self.finished = decision
        summary = {"analysis": k, "cutoff": cutoff, **counts, **result, "decision": decision or "continue"}
        self.analyses.append(summary)
        ctx.record("analysis", summary["decision"], basis=summary, evidence=evidence, actor=command.actor,
                   reason=command.reason)

    def _recommend(self, command: Command, ctx: Context) -> None:
        item = {"n": len(self.recommendations) + 1, **command.args, "actor": command.actor}
        self.recommendations.append(item)
        ctx.record("action", "recommendation", basis=item, actor=command.actor, reason=command.reason)

    def _decide(self, command: Command, ctx: Context) -> None:
        item = {"n": len(self.decisions) + 1, **command.args, "actor": command.actor}
        self.decisions.append(item)
        ctx.record("action", "accepted" if command.args["accept"] else "rejected", basis=item, actor=command.actor,
                   reason=command.reason)

    def _unblind(self, command: Command, ctx: Context) -> None:
        self.unblinded = True
        ctx.record("action", "unblinded", actor=command.actor, reason=command.reason)

    def _close(self, command: Command, ctx: Context) -> None:
        self.closed = True
        ctx.record("action", "closed", basis={"analyses": len(self.analyses), "finished": self.finished},
                   actor=command.actor, reason=command.reason)


# ---------------------------------------------------------------------- records made easy


def _seconds(value: float | str) -> float:
    return parse_timestamp(value) if isinstance(value, str) else float(value)


def patient(source: str, seq: int, pid: str, arm: str, entry: float | str, stratum: str | None = None) -> Record:
    """A patient's arm, entry (randomisation) and stratum."""
    t = _seconds(entry)
    return Record(source, seq, {"type": "patient", "patient": pid, "arm": arm, "entry": t, "stratum": stratum}, t)


def outcome(source: str, seq: int, pid: str, value: float, time: float | str) -> Record:
    """A patient's outcome (0 or 1, or a measurement), observed at `time`."""
    t = _seconds(time)
    return Record(source, seq, {"type": "outcome", "patient": pid, "value": value, "time": t}, t)


def event(source: str, seq: int, pid: str, time: float | str, happened: bool) -> Record:
    """The event at `time` (happened=True), or the last contact without it (False: censored then)."""
    t = _seconds(time)
    return Record(source, seq, {"type": "event", "patient": pid, "time": t, "event": happened}, t)


# ---------------------------------------------------------------------- the trial


class Trial:
    """A trial open in an engine: the only writer of its task until `close`."""

    def __init__(self, task: Task):
        if task.kind != KIND:
            raise TaskError(f"task {task.task_id!r} is a {task.kind!r}, not a {KIND!r}")
        self.task = task
        handler = task.handler
        assert isinstance(handler, TrialHandler)
        self.handler: TrialHandler = handler

    @classmethod
    def create(cls, engine: Engine, task_id: str, spec: TrialSpec) -> Trial:
        """A new trial from a spec (designed, for a group sequential design)."""
        return cls(engine.create(task_id, KIND, spec.data, TrialHandler(spec)))

    @classmethod
    def open(cls, engine: Engine, task_id: str) -> Trial:
        return cls(engine.open(task_id, TrialHandler))

    @property
    def spec(self) -> TrialSpec:
        return self.handler.spec

    @property
    def analyses(self) -> list[dict[str, Any]]:
        return list(self.handler.analyses)

    def ingest(self, records: Iterable[Record]) -> IngestResult:
        """Take a batch of records: all of it, or none if one of them is invalid (InvalidRecord)."""
        batch = list(records)
        for i, r in enumerate(batch):
            problem = record_problem(self.spec, r.data) if isinstance(r, Record) else "not a Record"
            if problem is not None:
                raise InvalidRecord(f"record {i} of the batch ({getattr(r, 'source', '?')}, "
                                    f"{getattr(r, 'seq', '?')}): {problem}; nothing of the batch was accepted")
        return self.task.ingest(batch)

    def _command(self, name: str, args: dict[str, Any], actor: str, reason: str) -> AuditEntry:
        problem = self.handler.problem(Command(name, args, actor, reason))
        if problem is not None:
            raise TaskError(problem)
        result = self.task.command(name, args, actor=actor, reason=reason)
        return result.entries[-1]

    def analyse(self, cutoff: float | str, *, actor: str, reason: str, final: bool = False) -> dict[str, Any]:
        """An interim (or, with final=True, the final) analysis of the data as of `cutoff`; its summary."""
        self._command("analyse", {"cutoff": _seconds(cutoff), "final": final}, actor, reason)
        return dict(self.handler.analyses[-1])

    def recommend(self, analysis: int, action: str, *, actor: str, reason: str) -> AuditEntry:
        """A recommendation of the data monitoring committee after an analysis."""
        return self._command("recommend", {"analysis": analysis, "action": action}, actor, reason)

    def decide(self, recommendation: int, accept: bool, *, actor: str, reason: str) -> AuditEntry:
        """The sponsor accepts or rejects a recommendation."""
        return self._command("decide", {"recommendation": recommendation, "accept": accept}, actor, reason)

    def unblind(self, *, actor: str, reason: str) -> AuditEntry:
        return self._command("unblind", {}, actor, reason)

    def close_trial(self, *, actor: str, reason: str) -> AuditEntry:
        """End the trial: later records are counted and ignored, no analysis follows."""
        return self._command("close", {}, actor, reason)

    def report(self, part: str = "open") -> Any:
        from .report import dsmb_report

        return dsmb_report(self.spec, self.task.task_id, self.handler.state(), self.task.audit(), part)

    def protocol(self, analysis: int) -> Any:
        from .report import analysis_protocol

        return analysis_protocol(self.spec, self.task.task_id, self.handler.state(), self.task.audit(), analysis)

    def close(self) -> None:
        """Close the task (not the trial: see `close_trial`); it can be opened again."""
        self.task.close()

    def __enter__(self) -> Trial:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()
