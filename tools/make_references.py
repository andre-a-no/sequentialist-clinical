# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Novikov Laboratories LLC
# Commercial licenses for use outside the AGPL: see COMMERCIAL.md
"""Write tests/data/references.json: data sets and the values external implementations compute for them.

Run by hand, in an environment with scipy and lifelines (neither is a
dependency of the package or of its tests):

    python -m venv /tmp/ref && /tmp/ref/bin/pip install scipy lifelines
    /tmp/ref/bin/python tools/make_references.py

References:

* boundaries of group sequential designs: computed here from the
  multivariate normal distribution of (Z_1, ..., Z_K), corr(Z_i, Z_j) =
  sqrt(t_i / t_j), with scipy's multivariate_normal.cdf (Genz's
  quasi-Monte Carlo integration): an algorithm independent of the
  recursive numerical integration of sequential-analysis;
* the log-rank statistic: lifelines.statistics.logrank_test, on survival
  times without ties (with ties, sequential-analysis takes tied events one
  at a time, lifelines and R's survdiff use the hypergeometric variance);
  and on the CDISC pilot data (tests/data/cdiscpilot01), which has ties;
* two proportions: scipy.stats.chi2_contingency without continuity
  correction (the square of the pooled Z);
* two means: scipy.stats.ttest_ind with equal variances (Student's t).
"""
import json
import math
import random
import sys
from pathlib import Path

import lifelines
import numpy as np
import scipy
from lifelines.statistics import logrank_test
from scipy.optimize import brentq
from scipy.stats import chi2_contingency, multivariate_normal, norm, ttest_ind

ROOT = Path(__file__).resolve().parent.parent
TARGET = ROOT / "tests" / "data" / "references.json"


def spending(kind, alpha, t, rho=None):
    if kind == "obrien-fleming":
        return 2 - 2 * norm.cdf(norm.ppf(1 - alpha / 2) / math.sqrt(t))
    if kind == "pocock":
        return alpha * math.log(1 + (math.e - 1) * t)
    return alpha * t**rho


def boundaries(alpha, fractions, kind, rho=None, two_sided=False):
    """Efficacy boundaries by the multivariate normal distribution, one analysis after the other."""
    side = alpha / 2 if two_sided else alpha
    ts = list(fractions)
    bounds = []
    previous = 0.0
    for k, t in enumerate(ts):
        target = (spending(kind, side, t, rho) if k + 1 < len(ts) else side) - previous
        target *= 2 if two_sided else 1
        cov = np.array([[math.sqrt(min(a, b) / max(a, b)) for b in ts[:k + 1]] for a in ts[:k + 1]])

        def crossing(c, bounds=bounds, cov=cov, k=k):
            # P(no crossing before k, crossing at k)
            if k == 0:
                return (2 if two_sided else 1) * norm.sf(c)
            lower = [-b if two_sided else -np.inf for b in bounds]
            continue_before = _box(cov[:k, :k], lower, bounds)
            continue_through = _box(cov, [*lower, -c if two_sided else -np.inf], [*bounds, c])
            return continue_before - continue_through

        bounds.append(brentq(lambda c, crossing=crossing, target=target: crossing(c) - target, 0.5, 8.0,
                             xtol=1e-10))
        previous += target / (2 if two_sided else 1)
    return bounds


def _box(cov, lower, upper):
    """P(lower < X < upper) for X ~ N(0, cov)."""
    dist = multivariate_normal(mean=np.zeros(len(upper)), cov=cov, maxpts=10_000_000, abseps=1e-10, releps=1e-10,
                               seed=1)
    return float(dist.cdf(np.array(upper, dtype=float), lower_limit=np.array(lower, dtype=float)))


def survival_data(rng, n, hazard_ratio):
    rows = []
    for i in range(n):
        arm = i % 2  # 1 treatment, 0 control
        entry = rng.uniform(0, 400)
        rate = math.log(2) / 300 * (hazard_ratio if arm == 1 else 1.0)
        time = rng.expovariate(rate)
        dropout = rng.expovariate(1 / 900)
        follow = min(time, dropout)
        rows.append({"id": f"s{i:04d}", "arm": arm, "entry": round(entry, 6), "end": round(entry + follow, 6),
                     "event": time <= dropout})
    return rows


def main() -> int:
    rng = random.Random(20261028)
    out = {"_generated_by": "tools/make_references.py", "_versions": {"scipy": scipy.__version__,
                                                                       "lifelines": lifelines.__version__}}
    designs = []
    for alpha, fractions, kind, rho, two_sided in [
        (0.025, [1 / 3, 2 / 3, 1.0], "obrien-fleming", None, False),
        (0.025, [0.25, 0.5, 0.75, 1.0], "pocock", None, False),
        (0.025, [0.4, 1.0], "power", 2.0, False),
        (0.05, [0.5, 0.75, 1.0], "obrien-fleming", None, True),
    ]:
        designs.append({"alpha": alpha, "analyses": fractions, "spending": kind, "rho": rho, "two_sided": two_sided,
                        "bounds": boundaries(alpha, fractions, kind, rho, two_sided)})
    out["boundaries"] = designs

    # survival without ties: the cut-off censors whoever is still followed then
    rows = survival_data(rng, 240, 0.65)
    cutoff = 700.0
    treated, control = [], []
    for r in rows:
        if r["entry"] > cutoff:
            continue
        end = min(r["end"], cutoff)
        event = r["event"] and r["end"] <= cutoff
        (treated if r["arm"] == 1 else control).append((end - r["entry"], event))
    result = logrank_test([d for d, _ in treated], [d for d, _ in control], [e for _, e in treated],
                          [e for _, e in control])
    observed = sum(e for _, e in treated)
    out["survival"] = {"rows": rows, "cutoff": cutoff, "chi2": float(result.test_statistic),
                       "events": [observed, sum(e for _, e in control)],
                       "patients": [len(treated), len(control)]}

    # the CDISC pilot: time to the first dermatologic event, high dose against placebo
    pilot = json.loads((ROOT / "tests" / "data" / "cdiscpilot01" / "adtte.json").read_text(encoding="utf-8"))
    names = [c["name"] for c in pilot["columns"]]
    tte = [dict(zip(names, row)) for row in pilot["rows"]]
    high = [r for r in tte if r["TRTP"] == "Xanomeline High Dose"]
    placebo = [r for r in tte if r["TRTP"] == "Placebo"]
    result = logrank_test([r["AVAL"] for r in high], [r["AVAL"] for r in placebo],
                          [r["CNSR"] == 0 for r in high], [r["CNSR"] == 0 for r in placebo])
    out["pilot"] = {"chi2": float(result.test_statistic), "patients": [len(high), len(placebo)],
                    "events": [sum(r["CNSR"] == 0 for r in high), sum(r["CNSR"] == 0 for r in placebo)]}

    # two proportions
    binary = [{"id": f"b{i:04d}", "arm": i % 3 % 2, "value": int(rng.random() < (0.42 if i % 3 % 2 else 0.3))}
              for i in range(301)]
    table = [[sum(r["value"] for r in binary if r["arm"] == a), sum(1 - r["value"] for r in binary if r["arm"] == a)]
             for a in (1, 0)]
    chi2 = chi2_contingency(table, correction=False)[0]
    out["binary"] = {"rows": binary, "chi2": float(chi2)}

    # two means
    means = [{"id": f"m{i:04d}", "arm": i % 2, "value": round(rng.gauss(0.4 if i % 2 else 0.0, 1.3), 6)}
             for i in range(121)]
    t = ttest_ind([r["value"] for r in means if r["arm"] == 1], [r["value"] for r in means if r["arm"] == 0])
    out["continuous"] = {"rows": means, "t": float(t.statistic)}

    TARGET.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {TARGET.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
