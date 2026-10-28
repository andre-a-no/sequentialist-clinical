# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Novikov Laboratories LLC
# Commercial licenses for use outside the AGPL: see COMMERCIAL.md
"""The type I error and the power of a whole trial pipeline, by simulation.

Each simulated trial runs through the same code as a real one: records of
synthetic patients go into a `Trial` (in memory, on the data's clock), and
analyses are done at the data cut-offs where the planned information
fractions are reached.  What is counted is the final decision of the
test; the Monte Carlo error is stated as a Wilson interval
(`sequentialist.core.proportion`).

The synthetic trial: patients enter one a day and are allocated to the
treatment arm with the probability of the design's allocation.

* binary: the outcome is observed 30 days after entry, with the rates of
  the assumptions (both at p_control under the null);
* continuous: normal outcomes 30 days after entry, the difference of means
  of the assumptions (0 under the null), standard deviation sigma;
* time to event: exponential times, a median of 365 days in the control
  arm and the hazard ratio of the assumptions (1 under the null); 1.5 times
  as many patients as the maximum number of events, followed until the
  event.

These are the conditions under which the design's numbers hold; a real
trial's accrual, drop-out and delays are not modelled.
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Any

from sequentialist.core import Engine, MemoryStorage, Proportion, proportion, run_simulations

from .spec import TrialSpec
from .trial import Trial, event, outcome, patient

DAY = 86400.0


@dataclass(frozen=True)
class SimulationResult:
    under: str  # "null" or "alternative"
    runs: int
    rejected: Proportion  # H0 rejected at some analysis
    stopped_early: Proportion  # the test decided before the last planned analysis
    mean_fraction: float  # mean information at the decision, over the maximum

    def to_json(self) -> dict[str, Any]:
        return {"under": self.under, "runs": self.runs, "rejected": self.rejected.to_json(),
                "stopped_early": self.stopped_early.to_json(), "mean_fraction": self.mean_fraction}


def _looks(spec: TrialSpec) -> list[float]:
    designed = spec.data["designed"]
    fractions = designed["fractions"] if designed is not None else spec.design["analyses"]
    if fractions is None or spec.maximum is None:
        raise ValueError("a simulation needs planned analyses and a maximum (design the trial, or give "
                         "design.analyses and design.maximum)")
    return [float(t) for t in fractions]


def _records(spec: TrialSpec, rng: random.Random, under: str) -> tuple[list, list[float]]:
    """The synthetic records and the time at which each unit of information arrives, in order."""
    assume, maximum = spec.design["assumptions"], int(spec.maximum)  # type: ignore[arg-type]
    arms = spec.data["arms"]
    strata = spec.strata()
    alt = under == "alternative"
    records: list = []
    arrivals: list[float] = []
    n = math.ceil(1.5 * maximum) if spec.endpoint == "time-to-event" else maximum
    for i in range(n):
        treated = rng.random() < assume["allocation"]
        pid = f"p{i:06d}"
        entry = i * DAY
        stratum = strata[rng.randrange(len(strata))] if strata else None
        records.append(patient("sim", len(records), pid, arms["treatment" if treated else "control"], entry, stratum))
        if spec.endpoint == "binary":
            p = assume["p_treatment"] if (treated and alt) else assume["p_control"]
            value = 1 if rng.random() < p else 0
            records.append(outcome("sim", len(records), pid, value, entry + 30 * DAY))
            arrivals.append(entry + 30 * DAY)
        elif spec.endpoint == "continuous":
            sigma = spec.data["endpoint"]["sigma"] or assume["sigma"]
            mean = assume["effect"] if (treated and alt) else 0.0
            records.append(outcome("sim", len(records), pid, rng.gauss(mean, sigma), entry + 30 * DAY))
            arrivals.append(entry + 30 * DAY)
        else:
            rate = math.log(2.0) / 365.0 * (assume["hazard_ratio"] if (treated and alt) else 1.0)
            happened = entry + rng.expovariate(rate) * DAY
            records.append(event("sim", len(records), pid, happened, True))
            arrivals.append(happened)
    return records, sorted(arrivals)


def simulate_trial(spec: TrialSpec, *, under: str = "alternative", runs: int = 200, seed: int = 0) -> SimulationResult:
    """Simulate `runs` trials under the null or the alternative of the assumptions; count the decisions."""
    if under not in ("null", "alternative"):
        raise ValueError("under must be 'null' or 'alternative'")
    looks = _looks(spec)
    maximum = int(spec.maximum)  # type: ignore[arg-type]
    free = TrialSpec({**spec.data, "roles": {"unblinded": [], "dsmb": [], "sponsor": []}})

    def one(rng: random.Random, i: int) -> tuple[bool, bool, float]:
        records, arrivals = _records(free, rng, under)
        engine = Engine(MemoryStorage(), clock="event", checkpoint_every=10**9)
        trial = Trial.create(engine, f"simulation-{i}", free)
        trial.ingest(records)
        summary: dict[str, Any] = {}
        for k, t in enumerate(looks):
            n = min(len(arrivals), math.ceil(t * maximum - 1e-9))
            last = k + 1 == len(looks)
            summary = trial.analyse(arrivals[n - 1], actor="simulation", reason="planned look", final=last)
            if trial.handler.finished is not None:
                break
        trial.close()
        return (summary["decision"] == "reject_h0", summary["analysis"] < len(looks),
                summary["information"] / maximum)

    results = run_simulations(one, runs, seed=seed)
    return SimulationResult(under, runs, proportion(sum(r[0] for r in results), runs),
                            proportion(sum(r[1] for r in results), runs),
                            math.fsum(r[2] for r in results) / runs)
