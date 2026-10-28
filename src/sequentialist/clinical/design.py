# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Novikov Laboratories LLC
# Commercial licenses for use outside the AGPL: see COMMERCIAL.md
"""Designing a group sequential trial once, before the first patient.

`design_trial` computes, with `group_sequential_design` of
sequential-analysis, the efficacy (and futility) boundaries at the planned
analyses, the maximum number of patients or events, the power and the
expected information at stopping, and writes them into the spec's
``designed``.  The trial is run from that spec: the plan in the audit
trail is the plan that was executed.  `check_trial` recomputes the design
and lists what differs (another version of sequential-analysis, an edited
file).

The maximum is that of rpact and gsDesign: the fixed-sample size (or
number of events, by Schoenfeld's formula) multiplied by the inflation
factor of the design.  For a time-to-event trial it is a number of events;
how many patients and how much follow-up yield them depends on accrual and
on the baseline hazard, which are not modelled here.
"""
from __future__ import annotations

import math
from typing import Any

from sequential_analysis import __version__ as _sa_version
from sequential_analysis.group_sequential import GroupSequentialDesign, group_sequential_design

from .spec import TrialSpec


def spending_of(name: str | None, rho: float | None) -> Any:
    """A spending function as sequential-analysis takes it: a name or ("power", rho)."""
    if name is None:
        return None
    return ("power", float(rho)) if name == "power" else name  # type: ignore[arg-type]


def core_design(spec: TrialSpec) -> GroupSequentialDesign:
    """The design of sequential-analysis behind a group sequential spec."""
    if not spec.group_sequential:
        raise ValueError("an anytime design has no group sequential design")
    d = spec.design
    futility = spending_of(d["futility"], d["futility_rho"])
    return group_sequential_design(
        d["alpha"], d["analyses"], spending_of(d["spending"], d["rho"]), two_sided=d["two_sided"],
        futility=futility, power=d["power"] if futility is not None else None)


def _maximum(spec: TrialSpec, design: GroupSequentialDesign) -> int:
    d = spec.design
    if d["maximum"] is not None:
        return int(d["maximum"])
    assume = d["assumptions"]
    power = None if design.futility is not None else d["power"]  # a futility design carries its own power
    allocation = assume["allocation"]
    if spec.endpoint == "binary":
        return design.max_sample_size_two_proportions(assume["p_control"], assume["p_treatment"], allocation, power)
    if spec.endpoint == "continuous":
        sigma = spec.data["endpoint"]["sigma"] or assume["sigma"]
        return design.max_sample_size_two_means(assume["effect"], sigma, allocation, power)
    return design.max_events(assume["hazard_ratio"], allocation, power)


def _plain(x: float) -> float | None:
    return float(x) if math.isfinite(x) else None


def design_trial(spec: TrialSpec) -> TrialSpec:
    """The spec with ``designed`` filled in; ValueError for an anytime design."""
    design = core_design(spec)
    maximum = _maximum(spec, design)
    power = spec.design["power"]
    drift: float | None = None
    if design.drift is not None:
        drift = design.drift
    elif power is not None:
        drift = design.drift_for_power(power)
    designed = {
        "unit": spec.unit,
        "maximum": maximum,
        "fractions": [float(t) for t in design.fractions],
        "bounds": [_plain(b) for b in design.bounds],
        "futility": None if design.futility is None else [float(f) for f in design.futility],
        "spent": [float(s) for s in design.spent],
        "drift": drift,
        "power": None if drift is None else design.power(drift),
        "expected_fraction": None if drift is None else design.expected_fraction(drift),
        "expected_fraction_null": design.expected_fraction(0.0),
        "sequential_analysis": _sa_version,
    }
    return TrialSpec({**spec.data, "designed": designed})


def check_trial(spec: TrialSpec, tolerance: float = 1e-9) -> list[str]:
    """What differs between the stored design and a recomputed one (empty: nothing)."""
    stored = spec.data["designed"]
    if stored is None:
        return ["the spec is not designed: run `sequentialist clinical design`"]
    fresh = design_trial(TrialSpec({**spec.data, "designed": None})).data["designed"]
    out = []
    for key in ("unit", "maximum"):
        if stored[key] != fresh[key]:
            out.append(f"{key}: stored {stored[key]!r}, recomputed {fresh[key]!r}")
    for key in ("fractions", "bounds", "futility", "spent"):
        a, b = stored[key], fresh[key]
        if (a is None) != (b is None) or (a is not None and len(a) != len(b)):
            out.append(f"{key}: stored {a!r}, recomputed {b!r}")
            continue
        for i, (x, y) in enumerate(zip(a or [], b or [])):
            if (x is None) != (y is None) or (x is not None and abs(x - y) > tolerance * max(1.0, abs(y))):
                out.append(f"{key}[{i}]: stored {x!r}, recomputed {y!r}")
    for key in ("drift", "power", "expected_fraction", "expected_fraction_null"):
        x, y = stored[key], fresh[key]
        if (x is None) != (y is None) or (x is not None and abs(x - y) > 1e-6 * max(1.0, abs(y))):
            out.append(f"{key}: stored {x!r}, recomputed {y!r}")
    return out
