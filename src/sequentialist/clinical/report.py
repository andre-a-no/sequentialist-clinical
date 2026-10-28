# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Novikov Laboratories LLC
# Commercial licenses for use outside the AGPL: see COMMERCIAL.md
"""Reports of a trial: the plan, the report of the data monitoring committee, the protocol of an analysis.

The report of the committee has two parts.  The **open** part shows what
does not reveal the comparison: enrolment and follow-up pooled over the
arms (counts below `small_cells` shown as "<n"), when the analyses
happened and how much information they had, and the recommendations and
decisions.  The **closed** part, for the committee and the unblinded
statistician, adds the results by arm and the statistics, boundaries and
decisions of every analysis, each with its guarantee.

The package separates what each part shows; it does not control who reads
the storage or the reports.

Effect estimates at the end of a group sequential trial are biased; the
closed part gives the rates or means by arm as plain descriptions
(guarantee "none"), not adjusted estimates.
"""
from __future__ import annotations

import math
import time as _time
from typing import Any, Mapping, Sequence, Union

from sequentialist.core import AuditEntry, Finding, Report, Section, Table
from sequentialist.core.schema import format_timestamp

from .spec import TrialSpec

Entry = Union[AuditEntry, Mapping[str, Any]]


def _entries(audit: Sequence[Entry]) -> list[AuditEntry]:
    return [e if isinstance(e, AuditEntry) else AuditEntry.from_json(dict(e)) for e in audit]


def _now(generated_at: float | None) -> float:
    return float(_time.time()) if generated_at is None else float(generated_at)


def _describe(spec: TrialSpec) -> str:
    d, e = spec.design, spec.data["endpoint"]
    sides = "two-sided" if d["two_sided"] else "one-sided"
    benefit = "higher" if e["alternative"] == "greater" else "lower"
    text = (f"Endpoint: {spec.endpoint} ({benefit} is better), treatment arm `{spec.data['arms']['treatment']}` "
            f"against control `{spec.data['arms']['control']}`. ")
    if spec.group_sequential:
        spending = d["spending"] + (f" (rho = {d['rho']})" if d["rho"] is not None else "")
        text += f"Group sequential design, {sides} alpha = {d['alpha']}, alpha spending {spending}"
        if d["futility"] is not None:
            text += f", non-binding futility boundaries (beta spending {d['futility']})"
        return text + "."
    return text + f"Anytime-valid design, {sides} alpha = {d['alpha']}: valid at any number of looks."


def plan_report(spec: TrialSpec, generated_at: float | None = None) -> Report:
    """The plan: boundaries at the planned analyses, maximum information, power, expected information."""
    report = Report(f"Plan of trial {spec.data['name']}", spec.data["name"], spec.id, _now(generated_at),
                    summary=_describe(spec))
    designed = spec.data["designed"]
    if not spec.group_sequential:
        section = Section("Design")
        section.findings.append(Finding(
            "The type I error is at most alpha however often and whenever the data are analysed.",
            _anytime_guarantee(spec), {"alpha": spec.design["alpha"], "maximum": spec.maximum}))
        report.sections.append(section)
        return report
    if designed is None:
        raise ValueError("the spec is not designed: run `sequentialist clinical design`")
    maximum, unit = designed["maximum"], designed["unit"]
    rows = []
    futility = designed["futility"] or [None] * len(designed["fractions"])
    for k, (t, b, f, s) in enumerate(zip(designed["fractions"], designed["bounds"], futility, designed["spent"])):
        rows.append([k + 1, t, _count(t, maximum), b, f, s])
    section = Section("Planned analyses", text=f"Maximum information: {maximum} {unit}. Boundaries are for the "
                                                "standardised statistic Z; at the analyses actually done they are "
                                                "recomputed for the information reached (alpha spending).")
    section.tables.append(Table(["analysis", "fraction", unit, "efficacy Z", "futility Z", "cumulative alpha"], rows))
    values: dict[str, Any] = {"alpha": spec.design["alpha"], "maximum": maximum}
    for key in ("power", "drift", "expected_fraction", "expected_fraction_null"):
        if designed[key] is not None:
            values[key] = designed[key]
    if designed["expected_fraction"] is not None:
        values[f"expected {unit} under the alternative"] = designed["expected_fraction"] * maximum
    section.findings.append(Finding(
        f"The type I error is at most alpha = {spec.design['alpha']}, futility boundaries ignored or obeyed; "
        f"the power is computed for the assumptions of the spec.", "asymptotic", values,
        notes=["Guarantee of the normal approximation of the statistics (sequential-analysis: "
               "GroupSequentialMonitor and the two-arm group sequential procedures)."]))
    report.sections.append(section)
    return report


def _count(fraction: float, maximum: int) -> int:
    return math.ceil(fraction * maximum - 1e-9)


def _anytime_guarantee(spec: TrialSpec) -> str:
    from .trial import _data_stream

    return str(_data_stream(spec).procedure.guarantee)


def _small(n: int, threshold: int) -> Any:
    return f"<{threshold}" if 0 < n < threshold else n


def _arms(spec: TrialSpec) -> tuple[str, str]:
    return spec.data["arms"]["treatment"], spec.data["arms"]["control"]


def _open_sections(spec: TrialSpec, state: Mapping[str, Any]) -> list[Section]:
    small = spec.data["small_cells"]
    strata: list[str | None] = [*spec.strata()] or [None]
    patients = state["patients"]
    rows = []
    for stratum in strata:
        group = [p for p in patients.values() if p["arm"] is not None and (stratum is None or p["stratum"] == stratum)]
        observed = sum(1 for p in group if (p["outcome"] if spec.endpoint != "time-to-event" else p["follow"]))
        events = sum(1 for p in group if p["follow"] is not None and p["follow"][1])
        row: list[Any] = [stratum or "all", _small(len(group), small), _small(observed, small)]
        if spec.endpoint == "time-to-event":
            row.append(_small(events, small))
        rows.append(row)
    columns = ["stratum", "patients", "with outcome" if spec.endpoint != "time-to-event" else "with follow-up"]
    if spec.endpoint == "time-to-event":
        columns.append("events")
    counts = state["counts"]
    enrolment = Section("Enrolment and data, pooled over the arms",
                        text=f"Records taken: {counts['records']}; corrections: {counts['updates']}; "
                             f"refused: {counts['rejected']}; after the trial was closed: {counts['after_close']}.")
    enrolment.tables.append(Table(columns, rows))
    analyses = Section("Analyses", text="When the analyses were done and how much information they had.")
    maximum = spec.maximum
    analyses.tables.append(Table(
        ["analysis", "cut-off", f"information ({spec.unit})", "fraction of maximum"],
        [[a["analysis"], format_timestamp(a["cutoff"]), a["information"],
          None if not maximum else a["information"] / maximum] for a in state["analyses"]]))
    actions = Section("Recommendations and decisions")
    decided = {d["recommendation"]: d for d in state["decisions"]}
    actions.tables.append(Table(
        ["recommendation", "after analysis", "action", "by", "sponsor"],
        [[r["n"], r["analysis"], r["action"], r["actor"],
          ("accepted" if decided[r["n"]]["accept"] else "rejected") if r["n"] in decided else None]
         for r in state["recommendations"]]))
    status = Section("Status", text=("The trial is closed. " if state["closed"] else "")
                     + ("The trial is unblinded." if state["unblinded"] else "The trial is blinded."))
    return [enrolment, analyses, actions, status]


def _analysis_entry(entries: Sequence[AuditEntry], k: int) -> AuditEntry:
    for e in entries:
        if e.kind == "analysis" and e.basis.get("analysis") == k:
            return e
    raise ValueError(f"there is no analysis {k}")


def _findings(spec: TrialSpec, summary: Mapping[str, Any], entry: AuditEntry | None) -> list[Finding]:
    treatment, control = _arms(spec)
    k, decision = summary["analysis"], summary["decision"]
    evidence = {e["stream"]: e for e in (entry.evidence if entry is not None else ())}
    when = format_timestamp(summary["cutoff"])
    out = []
    if spec.group_sequential:
        guarantee = evidence["monitor"]["guarantee"] if "monitor" in evidence else "asymptotic"
        bound = summary["bound"]
        text = (f"Analysis {k} (cut-off {when}, {summary['fraction']:.1%} of the maximum information): "
                f"Z = {summary['z']:.3f} against the efficacy boundary "
                f"{'none' if bound is None else format(bound, '.3f')}")
        if summary["futility"] is not None:
            text += f" and the futility boundary {summary['futility']:.3f}"
        text += f": {_verdict(decision)}."
        values = {key: summary[key] for key in ("z", "bound", "futility", "spent", "fraction", "information")}
        out.append(Finding(text, guarantee, values, stream="monitor",
                           procedure=evidence.get("monitor", {}).get("procedure"),
                           history_digest=evidence.get("monitor", {}).get("history_digest")))
    else:
        guarantee = evidence["data"]["guarantee"] if "data" in evidence else "exact"
        text = (f"Analysis {k} (cut-off {when}): e-value {_e(summary.get('log_e_value'))}, "
                f"p-value {_fmt(summary.get('p_value'))}: {_verdict(decision)}.")
        values = {key: summary[key] for key in ("log_e_value", "p_value", "hazard_ratio_estimate", "lower", "upper",
                                                "difference") if summary.get(key) is not None}
        notes = []
        if summary.get("lower") is not None:
            notes.append(f"Confidence sequence for the hazard ratio ({treatment} / {control}): "
                         f"[{summary['lower']:.3f}, {summary['upper']:.3f}], valid at any look.")
        out.append(Finding(text, guarantee, values, stream="data", procedure=evidence.get("data", {}).get("procedure"),
                           history_digest=evidence.get("data", {}).get("history_digest"), notes=notes))
    out.append(Finding(_by_arm(spec, summary, evidence.get("data")), "none",
                       {"patients": summary["patients"], "units": summary["units"]},
                       notes=["A description of the data, not an estimate adjusted for the sequential design."]))
    return out


def _by_arm(spec: TrialSpec, summary: Mapping[str, Any], data: Mapping[str, Any] | None) -> str:
    treatment, control = _arms(spec)
    units = summary["units"]
    if spec.endpoint == "binary":
        s = summary["successes"]
        return (f"{treatment}: {s[0]} of {units[0]} ({_rate(s[0], units[0])}); "
                f"{control}: {s[1]} of {units[1]} ({_rate(s[1], units[1])}).")
    if spec.endpoint == "time-to-event":
        e = summary["events"]
        return f"{treatment}: {e[0]} events in {units[0]} patients; {control}: {e[1]} events in {units[1]} patients."
    output = (data or {}).get("output") or {}
    return (f"{treatment}: mean {_fmt(output.get('mean_treatment'))} in {units[0]} patients; "
            f"{control}: mean {_fmt(output.get('mean_control'))} in {units[1]} patients.")


def _rate(a: int, n: int) -> str:
    return f"{a / n:.1%}" if n else "—"


def _fmt(x: Any) -> str:
    return "—" if x is None else f"{x:.4g}"


def _e(log_e: Any) -> str:
    if log_e is None:
        return "—"
    return f"{math.exp(log_e):.3g}" if log_e < 700 else f"exp({log_e:.1f})"


def _verdict(decision: str) -> str:
    return {"reject_h0": "H0 rejected (efficacy)", "accept_h0": "H0 not rejected (futility or final analysis)",
            "continue": "no boundary reached"}.get(decision, decision)


def dsmb_report(spec: TrialSpec, task_id: str, state: Mapping[str, Any], audit: Sequence[Entry], part: str = "open",
                generated_at: float | None = None) -> Report:
    """The report of the data monitoring committee: part "open" (pooled) or "closed" (by arm)."""
    if part not in ("open", "closed"):
        raise ValueError("part must be 'open' or 'closed'")
    entries = _entries(audit)
    title = f"Trial {spec.data['name']}: report to the data monitoring committee ({part} part)"
    summary = _describe(spec) if part == "closed" else ("Pooled over the arms; the comparison is in the closed "
                                                         "part.")
    report = Report(title, task_id, spec.id, _now(generated_at), summary=summary)
    report.sections += _open_sections(spec, state)
    if part == "closed":
        results = Section("Results by analysis")
        for a in state["analyses"]:
            results.findings += _findings(spec, a, _analysis_entry(entries, a["analysis"]))
        report.sections.insert(2, results)
    return report


def analysis_protocol(spec: TrialSpec, task_id: str, state: Mapping[str, Any], audit: Sequence[Entry], analysis: int,
                      generated_at: float | None = None) -> Report:
    """Everything about one analysis: the data cut, the result, its evidence, what was recommended and decided."""
    entries = _entries(audit)
    entry = _analysis_entry(entries, analysis)
    summary = entry.basis
    report = Report(f"Trial {spec.data['name']}: protocol of analysis {analysis}", task_id, spec.id,
                    _now(generated_at), summary=_describe(spec))
    left = summary["left_out"]
    data = Section("Data", text=(f"Cut-off {format_timestamp(summary['cutoff'])}, analysis run by {entry.actor} "
                                 f"({entry.reason}). Left out: {left['other_arm']} patients of other arms, "
                                 f"{left['incomplete']} without arm or entry, {left['no_follow_up']} without "
                                 f"follow-up, {left['inconsistent']} with follow-up before entry."))
    treatment, control = _arms(spec)
    rows = [["patients", *summary["patients"]], [spec.unit if spec.unit != "patients" else "with outcome",
                                                 *(summary["events"] if spec.unit == "events" else summary["units"])]]
    data.tables.append(Table(["", treatment, control], rows))
    result = Section("Result")
    result.findings += _findings(spec, summary, entry)
    evidence = Section("Evidence", text=f"Audit entry {entry.n}, hash `{entry.hash}`.")
    evidence.tables.append(Table(
        ["stream", "procedure", "guarantee", "inputs", "history digest"],
        [[e["stream"], f"{e['procedure']['name']}/{e['procedure']['version']}", e["guarantee"], e["t"],
          e["history_digest"]] for e in entry.evidence]))
    decided = {d["recommendation"]: d for d in state["decisions"]}
    actions = Section("Recommendations and decisions")
    actions.tables.append(Table(
        ["recommendation", "action", "by", "sponsor", "decided by"],
        [[r["n"], r["action"], r["actor"],
          ("accepted" if decided[r["n"]]["accept"] else "rejected") if r["n"] in decided else None,
          decided[r["n"]]["actor"] if r["n"] in decided else None]
         for r in state["recommendations"] if r["analysis"] == analysis]))
    report.sections += [data, result, evidence, actions]
    return report
