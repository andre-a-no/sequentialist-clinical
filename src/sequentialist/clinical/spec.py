# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Novikov Laboratories LLC
# Commercial licenses for use outside the AGPL: see COMMERCIAL.md
"""The spec of a trial (`clinical.trial/1`) and the records a trial takes.

A trial compares a treatment arm with a control arm on one primary
endpoint: binary (response, an adverse outcome), continuous (a measurement,
a change from baseline) or the time to an event (death, progression).  It
is monitored either by a group sequential design (the default: alpha
spending, efficacy and optional futility boundaries) or by an
anytime-valid test, whose type I error holds however often the data are
looked at.

    {
      "kind": "clinical.trial/1",
      "name": "rx-301",
      "endpoint": {"type": "time-to-event", "alternative": "less"},
      "arms": {"treatment": "A", "control": "B"},
      "design": {"alpha": 0.025, "analyses": [0.5, 0.75, 1.0], "spending": "obrien-fleming",
                 "power": 0.9, "assumptions": {"hazard_ratio": 0.7}}
    }

`alternative` is the direction of benefit: "greater" when the treatment
should raise the endpoint (response rate, mean), "less" when it should
lower it (hazard, a symptom score).  A two-sided design rejects in either
direction; the alternative then only fixes the sign of the statistic.

A group sequential spec is designed once (`design_trial`), which writes
the boundaries, the maximum number of units (or events) and the design's
operating characteristics into ``designed``; the trial runs from the
designed spec, and its spec id names the plan that was executed.

Records are facts about patients, identified by a pseudonym, never by a
name or a document number:

    {"type": "patient", "patient": "p-0001", "arm": "A", "entry": 1767225600, "stratum": "<65"}
    {"type": "outcome", "patient": "p-0001", "value": 1, "time": 1769904000}
    {"type": "event",   "patient": "p-0001", "time": 1772323200, "event": true}

An ``event`` record says when the patient had the event (``event``: true)
or was last known to be free of it (false: censored then).  A later record
about the same fact replaces the earlier one: new follow-up, a new data
cut, a correction.  Times are seconds since the epoch or RFC 3339
timestamps with a time zone.
"""
from __future__ import annotations

import math
from typing import Any

from sequentialist.core import Spec
from sequentialist.core.schema import (
    INVALID,
    Array,
    Boolean,
    Errors,
    Field,
    Nullable,
    Number,
    Object,
    OneOf,
    String,
    Timestamp,
)

KIND = "clinical.trial/1"
ENDPOINTS = ("binary", "continuous", "time-to-event")
DESIGNS = ("group-sequential", "anytime")
SPENDINGS = ("obrien-fleming", "pocock", "power")
METHODS = {"binary": ("grow", "universal"), "time-to-event": ("mixture", "grow")}
ACTIONS = ("continue", "stop-efficacy", "stop-futility", "stop-harm", "modify")
#: Pseudonyms: letters, digits and . _ : -, so that a name or a document number with spaces is refused.
PATIENT_PATTERN = r"[A-Za-z0-9._:-]{1,128}"

_label = String(min_length=1, max_length=64)
_people = Array(String(min_length=1, max_length=256), max_items=100)
_probability = Nullable(Number(exclusive_minimum=0, exclusive_maximum=1))

_endpoint = Object({
    "type": Field(String(enum=ENDPOINTS), doc="binary, continuous or time-to-event."),
    "alternative": Field(String(enum=["greater", "less"]),
                         doc="The direction of benefit: greater (the treatment raises the endpoint) or less "
                             "(it lowers it, a hazard for instance)."),
    "sigma": Field(Nullable(Number(exclusive_minimum=0)), default=None,
                   doc="Continuous endpoint: the known standard deviation, or null to estimate it."),
    "strata": Field(Array(_label, max_items=1000), default=[],
                    doc="Time to event: the labels of the strata of a stratified log-rank test (none: unstratified)."),
}, doc="The primary endpoint.")

_assumptions = Object({
    "p_control": Field(_probability, default=None, doc="Binary: the expected rate in the control arm."),
    "p_treatment": Field(_probability, default=None, doc="Binary: the rate the trial is powered for."),
    "effect": Field(Nullable(Number()), default=None,
                    doc="Continuous: the difference of means (treatment - control) the trial is powered for."),
    "sigma": Field(Nullable(Number(exclusive_minimum=0)), default=None,
                   doc="Continuous: the standard deviation assumed for the design when endpoint.sigma is null."),
    "hazard_ratio": Field(Nullable(Number(exclusive_minimum=0)), default=None,
                          doc="Time to event: the hazard ratio (treatment / control) the trial is powered for."),
    "allocation": Field(Number(exclusive_minimum=0, exclusive_maximum=1), default=0.5,
                        doc="The share of patients allocated to the treatment arm."),
}, doc="What the sample size is computed for.")

_design = Object({
    "type": Field(String(enum=DESIGNS), default="group-sequential",
                  doc="group-sequential (alpha spending, the default) or anytime (valid at any look)."),
    "alpha": Field(Number(exclusive_minimum=0, exclusive_maximum=0.5),
                   doc="The type I error: one-sided, or the total of both sides when two_sided."),
    "two_sided": Field(Boolean(), default=False),
    "spending": Field(String(enum=SPENDINGS), default="obrien-fleming",
                      doc="Alpha spending function (Lan-DeMets O'Brien-Fleming, Lan-DeMets Pocock, Kim-DeMets power)."),
    "rho": Field(Nullable(Number(exclusive_minimum=0)), default=None, doc="The exponent of the power family."),
    "analyses": Field(Nullable(Array(Number(exclusive_minimum=0, maximum=1), min_items=1, max_items=20)),
                      default=None,
                      doc="Planned analyses as fractions of the maximum information, increasing, the last 1."),
    "futility": Field(Nullable(String(enum=SPENDINGS)), default=None,
                      doc="Beta spending function of non-binding futility boundaries (one-sided designs)."),
    "futility_rho": Field(Nullable(Number(exclusive_minimum=0)), default=None),
    "power": Field(_probability, default=None, doc="The power the maximum information is computed for."),
    "maximum": Field(Nullable(Number(integer=True, minimum=1)), default=None,
                     doc="The maximum number of patients with an outcome (binary, continuous) or of events "
                         "(time to event); computed from power and assumptions when null."),
    "assumptions": Field(_assumptions, default=dict),
    "method": Field(Nullable(String(enum=sorted({m for ms in METHODS.values() for m in ms}))), default=None,
                    doc="Anytime designs: the method of the test of sequential-analysis (its default when null)."),
}, doc="How the trial is monitored.")

_designed = Object({
    "unit": Field(String(enum=["patients", "events"])),
    "maximum": Field(Number(integer=True, minimum=1)),
    "fractions": Field(Array(Number(exclusive_minimum=0, maximum=1), min_items=1)),
    "bounds": Field(Array(Nullable(Number()), min_items=1),
                    doc="Efficacy boundaries for Z at the planned analyses (null: no alpha spent there)."),
    "futility": Field(Nullable(Array(Number(), min_items=1)), default=None),
    "spent": Field(Array(Number(minimum=0), min_items=1), doc="Cumulative alpha (one side) by analysis."),
    "drift": Field(Nullable(Number()), default=None,
                   doc="theta * sqrt(maximum information) the design is powered for."),
    "power": Field(Nullable(Number()), default=None),
    "expected_fraction": Field(Nullable(Number()), default=None,
                               doc="Expected information at stopping over the maximum, under the alternative."),
    "expected_fraction_null": Field(Nullable(Number()), default=None, doc="The same under H0."),
    "sequential_analysis": Field(String(), doc="The version of sequential-analysis that computed the design."),
}, doc="Written by `sequentialist clinical design`; do not edit.")


class TrialSpec(Spec):
    kind = KIND
    schema = Object({
        "name": Field(String(min_length=1, max_length=200)),
        "endpoint": Field(_endpoint),
        "arms": Field(Object({
            "treatment": Field(_label, doc="The code of the treatment arm in the records."),
            "control": Field(_label, doc="The code of the control arm in the records."),
        }), doc="Arm codes as the records carry them; records of other arms are counted and ignored."),
        "design": Field(_design),
        "roles": Field(Object({
            "unblinded": Field(_people, default=[], doc="Who may run analyses (the unblinded statistician)."),
            "dsmb": Field(_people, default=[], doc="Members of the data monitoring committee."),
            "sponsor": Field(_people, default=[], doc="Who decides on recommendations, unblinds and closes."),
        }), default=dict, doc="Pseudonyms of the people who may give each command; an empty list lets anyone."),
        "small_cells": Field(Number(integer=True, minimum=1), default=5,
                             doc="Counts below this are shown as \"<n\" in the open report."),
        "designed": Field(Nullable(_designed), default=None),
    }, doc="A clinical trial monitored by interim analyses.")

    # ------------------------------------------------------------------ shortcuts

    @property
    def endpoint(self) -> str:
        return str(self.data["endpoint"]["type"])

    @property
    def design(self) -> dict[str, Any]:
        return dict(self.data["design"])

    @property
    def group_sequential(self) -> bool:
        return self.data["design"]["type"] == "group-sequential"

    @property
    def unit(self) -> str:
        return "events" if self.endpoint == "time-to-event" else "patients"

    @property
    def maximum(self) -> int | None:
        designed = self.data["designed"]
        if designed is not None:
            return int(designed["maximum"])
        given = self.data["design"]["maximum"]
        return None if given is None else int(given)

    def strata(self) -> list[str]:
        return list(self.data["endpoint"]["strata"])

    # ------------------------------------------------------------------ rules the schema cannot express

    def problems(self) -> list[tuple[str, str]]:
        out: list[tuple[str, str]] = []
        endpoint, design = self.data["endpoint"], self.data["design"]
        kind, assume = endpoint["type"], design["assumptions"]
        if endpoint["sigma"] is not None and kind != "continuous":
            out.append(("/endpoint/sigma", "only a continuous endpoint has a standard deviation"))
        if endpoint["strata"] and kind != "time-to-event":
            out.append(("/endpoint/strata", "strata are supported for the log-rank test only"))
        if len(set(endpoint["strata"])) != len(endpoint["strata"]):
            out.append(("/endpoint/strata", "the labels of the strata must differ"))
        if self.data["arms"]["treatment"] == self.data["arms"]["control"]:
            out.append(("/arms", "the treatment and the control arm need different codes"))
        for name, rho in (("spending", "rho"), ("futility", "futility_rho")):
            if design[name] == "power" and design[rho] is None:
                out.append((f"/design/{rho}", f"{name} \"power\" needs its exponent"))
            if design[name] != "power" and design[rho] is not None:
                out.append((f"/design/{rho}", f"only the power family has an exponent; {name} is {design[name]!r}"))
        analyses = design["analyses"]
        if analyses is not None:
            if any(b <= a for a, b in zip(analyses, analyses[1:])) or analyses[-1] != 1:
                out.append(("/design/analyses", "fractions must increase strictly and end at 1"))
        if design["type"] == "anytime":
            out += self._anytime_problems()
        else:
            out += self._group_sequential_problems()
        unknown = [k for k, v in assume.items() if v is not None and k != "allocation" and k not in _USES[kind]]
        for key in unknown:
            out.append((f"/design/assumptions/{key}", f"not used by a {kind} endpoint"))
        return out

    def _anytime_problems(self) -> list[tuple[str, str]]:
        out = []
        design, kind = self.data["design"], self.data["endpoint"]["type"]
        if kind == "continuous":
            out.append(("/endpoint/type", "sequential-analysis 0.4 has no anytime-valid two-sample test of means; "
                                          "use a group-sequential design"))
        for name in ("futility", "futility_rho", "rho"):
            if design[name] is not None:
                out.append((f"/design/{name}", "an anytime design has no spending functions"))
        if design["method"] is not None and design["method"] not in METHODS.get(kind, ()):
            out.append(("/design/method", f"a {kind} endpoint takes one of {list(METHODS.get(kind, ()))}"))
        if self.data["designed"] is not None:
            out.append(("/designed", "an anytime design is not designed in advance"))
        return out

    def _group_sequential_problems(self) -> list[tuple[str, str]]:
        out = []
        design = self.data["design"]
        endpoint = self.data["endpoint"]
        if design["analyses"] is None:
            out.append(("/design/analyses", "a group sequential design needs its planned analyses"))
        if design["method"] is not None:
            out.append(("/design/method", "only anytime designs take a method"))
        if design["futility"] is not None:
            if design["two_sided"]:
                out.append(("/design/futility", "futility boundaries are for one-sided designs"))
            if design["power"] is None:
                out.append(("/design/power", "futility boundaries are designed for a power"))
        if design["maximum"] is None:
            if design["power"] is None:
                out.append(("/design/power", "without a maximum, the power to compute it for is needed"))
            assume = design["assumptions"]
            missing = [k for k in _NEEDS[endpoint["type"]] if assume[k] is None]
            if endpoint["type"] == "continuous" and endpoint["sigma"] is None and assume["sigma"] is None:
                missing.append("sigma")
            for key in missing:
                out.append((f"/design/assumptions/{key}", "needed to compute the maximum (or give design.maximum)"))
            if (endpoint["type"] == "binary" and not missing
                    and assume["p_control"] == assume["p_treatment"]):
                out.append(("/design/assumptions", "the rates the trial is powered for must differ"))
            if endpoint["type"] == "time-to-event" and assume["hazard_ratio"] == 1:
                out.append(("/design/assumptions/hazard_ratio", "must differ from 1"))
            if endpoint["type"] == "continuous" and assume["effect"] == 0:
                out.append(("/design/assumptions/effect", "must differ from 0"))
        return out


_NEEDS = {"binary": ("p_control", "p_treatment"), "continuous": ("effect",), "time-to-event": ("hazard_ratio",)}
_USES = {"binary": ("p_control", "p_treatment"), "continuous": ("effect", "sigma"),
         "time-to-event": ("hazard_ratio",)}


# ---------------------------------------------------------------------- records

_patient_id = String(pattern=PATIENT_PATTERN, doc="A pseudonym (see sequentialist.core.pseudonymize).")

RECORD = OneOf("type", {
    "patient": Object({
        "patient": Field(_patient_id),
        "arm": Field(_label),
        "entry": Field(Timestamp(), doc="Randomisation (or start of follow-up)."),
        "stratum": Field(Nullable(_label), default=None),
    }),
    "outcome": Object({
        "patient": Field(_patient_id),
        "value": Field(Number(), doc="0 or 1 for a binary endpoint, a number for a continuous one."),
        "time": Field(Timestamp(), doc="When the outcome was observed."),
    }),
    "event": Object({
        "patient": Field(_patient_id),
        "time": Field(Timestamp(), doc="The time of the event, or of the last contact without it."),
        "event": Field(Boolean(), doc="true: the event happened then; false: censored then."),
    }),
}, doc="A fact about a patient.")


def normal_record(data: Any) -> tuple[dict[str, Any] | None, Errors]:
    """The normal form of a record's data and the problems found (an empty list when it is valid)."""
    errors: Errors = []
    normal = RECORD.check(data, "", errors)
    if errors or normal is INVALID:
        return None, errors
    for key in ("entry", "time"):
        if key in normal and not math.isfinite(normal[key]):  # pragma: no cover - Timestamp refuses these
            return None, [(f"/{key}", "must be finite")]
    return normal, []


def record_problem(spec: TrialSpec, data: Any) -> str | None:
    """Why a record cannot be taken by a trial of this spec (None when it can)."""
    normal, errors = normal_record(data)
    if normal is None:
        return "; ".join(f"{where or '/'}: {message}" for where, message in errors)
    if normal["type"] == "event" and spec.endpoint != "time-to-event":
        return f"an event record, but the endpoint is {spec.endpoint}"
    if normal["type"] == "outcome":
        if spec.endpoint == "time-to-event":
            return "an outcome record, but the endpoint is the time to an event"
        if spec.endpoint == "binary" and normal["value"] not in (0, 1):
            return f"a binary outcome must be 0 or 1, got {normal['value']!r}"
    if normal["type"] == "patient":
        strata = spec.strata()
        stratum = normal["stratum"]
        if strata and stratum not in strata:
            return f"stratum {stratum!r} is not one of {strata}"
        if not strata and stratum is not None:
            return "this trial has no strata"
    return None
