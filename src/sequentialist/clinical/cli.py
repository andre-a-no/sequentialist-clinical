# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Novikov Laboratories LLC
# Commercial licenses for use outside the AGPL: see COMMERCIAL.md
"""`sequentialist clinical`: design a trial, check it, simulate it, report on it.

    sequentialist clinical design SPEC [-o DESIGNED]        compute the boundaries and the maximum, once
    sequentialist clinical check DESIGNED                   recompute and list what differs
    sequentialist clinical plan SPEC [--json]               the plan: boundaries, power, expected information
    sequentialist clinical simulate SPEC [--runs N]         type I error and power of the whole pipeline
    sequentialist clinical report STORAGE TRIAL [--part open|closed] [--json]
    sequentialist clinical protocol STORAGE TRIAL ANALYSIS [--json]
    sequentialist clinical schema                           JSON Schema of clinical.trial/1

STORAGE is ``files:<directory>`` or ``sqlite:<file>``.  Reports are made
from the trial's last checkpoint, without opening (and locking) the task.
"""
from __future__ import annotations

import argparse
import json
import sys
from typing import Any

from sequentialist.core import json_schema, load_spec, open_storage

from .design import check_trial, design_trial
from .report import analysis_protocol, dsmb_report, plan_report
from .simulate import simulate_trial
from .spec import KIND, TrialSpec


def _print(report: Any, as_json: bool) -> None:
    print(json.dumps(report.to_json(), indent=2, ensure_ascii=False) if as_json else report.to_markdown())


def _design(args: argparse.Namespace) -> int:
    spec = design_trial(load_spec(args.spec, TrialSpec))
    text = json.dumps(spec.to_json(), indent=2, ensure_ascii=False) + "\n"
    if args.output:
        with open(args.output, "w", encoding="utf-8") as f:
            f.write(text)
        designed = spec.data["designed"]
        print(f"{args.output}: maximum {designed['maximum']} {designed['unit']}, spec id {spec.id}", file=sys.stderr)
    else:
        sys.stdout.write(text)
    return 0


def _check(args: argparse.Namespace) -> int:
    problems = check_trial(load_spec(args.spec, TrialSpec))
    for p in problems:
        print(p)
    if not problems:
        print(f"{args.spec}: the design is what sequential-analysis computes", file=sys.stderr)
    return 1 if problems else 0


def _plan(args: argparse.Namespace) -> int:
    _print(plan_report(load_spec(args.spec, TrialSpec)), args.json)
    return 0


def _simulate(args: argparse.Namespace) -> int:
    spec = load_spec(args.spec, TrialSpec)
    results = [simulate_trial(spec, under=u, runs=args.runs, seed=args.seed)
               for u in (("null", "alternative") if args.under == "both" else (args.under,))]
    if args.json:
        print(json.dumps([r.to_json() for r in results], indent=2))
        return 0
    for r in results:
        lo, hi = r.rejected.interval
        early_lo, early_hi = r.stopped_early.interval
        print(f"under the {r.under}: H0 rejected in {r.rejected.hits} of {r.runs} trials, {r.rejected.rate:.4f} "
              f"(95% Monte Carlo interval {lo:.4f} to {hi:.4f}); decided early in {r.stopped_early.rate:.3f} "
              f"({early_lo:.3f} to {early_hi:.3f}); mean information at the decision {r.mean_fraction:.3f} "
              f"of the maximum")
    return 0


def _stored(storage_name: str, task_id: str) -> tuple[TrialSpec, dict[str, Any], list[dict[str, Any]]]:
    with open_storage(storage_name) as storage:
        loaded = storage.load(task_id)
        if loaded is None:
            raise ValueError(f"no task {task_id!r} in {storage_name}")
        state = loaded[1]
        audit = storage.audit_read(task_id)
    if state.get("kind") != KIND:
        raise ValueError(f"task {task_id!r} is a {state.get('kind')!r}, not a trial")
    return TrialSpec(state["spec"]), state["handler"], audit


def _report(args: argparse.Namespace) -> int:
    spec, state, audit = _stored(args.storage, args.trial)
    _print(dsmb_report(spec, args.trial, state, audit, args.part), args.json)
    return 0


def _protocol(args: argparse.Namespace) -> int:
    spec, state, audit = _stored(args.storage, args.trial)
    _print(analysis_protocol(spec, args.trial, state, audit, args.analysis), args.json)
    return 0


def _schema(args: argparse.Namespace) -> int:
    print(json.dumps(json_schema(TrialSpec), indent=2, ensure_ascii=False))
    return 0


def register(commands: Any) -> None:
    """Add `clinical` to the `sequentialist` command (entry point sequentialist.commands)."""
    parser = commands.add_parser("clinical", help="interim monitoring of clinical trials",
                                 description="Interim monitoring of clinical trials.")
    sub = parser.add_subparsers(dest="clinical_command")

    p = sub.add_parser("design", help="compute the boundaries and the maximum of a group sequential trial")
    p.add_argument("spec")
    p.add_argument("-o", "--output", help="write the designed spec here (default: standard output)")
    p.set_defaults(handler=_design)

    p = sub.add_parser("check", help="recompute a designed spec and list what differs")
    p.add_argument("spec")
    p.set_defaults(handler=_check)

    p = sub.add_parser("plan", help="the plan: boundaries, power, expected information")
    p.add_argument("spec")
    p.add_argument("--json", action="store_true")
    p.set_defaults(handler=_plan)

    p = sub.add_parser("simulate", help="type I error and power of the whole pipeline, by simulation")
    p.add_argument("spec")
    p.add_argument("--under", choices=["null", "alternative", "both"], default="both")
    p.add_argument("--runs", type=int, default=1000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--json", action="store_true")
    p.set_defaults(handler=_simulate)

    p = sub.add_parser("report", help="the report to the data monitoring committee")
    p.add_argument("storage")
    p.add_argument("trial")
    p.add_argument("--part", choices=["open", "closed"], default="open")
    p.add_argument("--json", action="store_true")
    p.set_defaults(handler=_report)

    p = sub.add_parser("protocol", help="the protocol of one analysis")
    p.add_argument("storage")
    p.add_argument("trial")
    p.add_argument("analysis", type=int)
    p.add_argument("--json", action="store_true")
    p.set_defaults(handler=_protocol)

    p = sub.add_parser("schema", help="JSON Schema of clinical.trial/1")
    p.set_defaults(handler=_schema)

    parser.set_defaults(handler=lambda args: (parser.print_help(), 2)[1])
