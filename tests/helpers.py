# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Novikov Laboratories LLC
# Commercial licenses for use outside the AGPL: see COMMERCIAL.md
"""Specs and data for the tests."""
import json
from pathlib import Path

from sequentialist.clinical import TrialSpec, design_trial, event, outcome, patient
from sequentialist.core import parse_spec

DATA = Path(__file__).resolve().parent / "data"
DAY = 86400.0
T0 = 1767225600.0  # 2026-01-01T00:00:00Z


def references():
    return json.loads((DATA / "references.json").read_text(encoding="utf-8"))


def spec(endpoint="time-to-event", alternative="less", design=None, more_endpoint=None, **extra):
    """A spec; `design` is merged into a default group sequential design for the endpoint."""
    assumptions = {"binary": {"p_control": 0.3, "p_treatment": 0.45}, "continuous": {"effect": 0.5, "sigma": 1.0},
                   "time-to-event": {"hazard_ratio": 0.7}}[endpoint]
    base = {"alpha": 0.025, "analyses": [0.5, 0.75, 1.0], "power": 0.9, "assumptions": assumptions}
    base.update(design or {})
    data = {"kind": "clinical.trial/1", "name": "test", "endpoint": {"type": endpoint, "alternative": alternative},
            "arms": {"treatment": "TRT", "control": "CTL"}, "design": base}
    data["endpoint"].update(more_endpoint or {})
    data.update(extra)
    s = parse_spec(data, TrialSpec)
    return design_trial(s) if s.group_sequential else s


def survival_records(rows, source="cut", arms=("TRT", "CTL")):
    """Records of the reference survival data: arm 1 is the treatment."""
    out = []
    for r in rows:
        out.append(patient(source, len(out), r["id"], arms[0] if r["arm"] == 1 else arms[1], T0 + r["entry"] * DAY))
        out.append(event(source, len(out), r["id"], T0 + r["end"] * DAY, r["event"]))
    return out


def outcome_records(rows, source="cut", arms=("TRT", "CTL")):
    out = []
    for i, r in enumerate(rows):
        entry = T0 + i * DAY
        out.append(patient(source, len(out), r["id"], arms[0] if r["arm"] == 1 else arms[1], entry))
        out.append(outcome(source, len(out), r["id"], r["value"], entry + 10 * DAY))
    return out
