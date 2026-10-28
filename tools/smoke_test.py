# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Novikov Laboratories LLC
# Commercial licenses for use outside the AGPL: see COMMERCIAL.md
"""Smoke test of an installed sequentialist-clinical: `python tools/smoke_test.py VERSION`, outside the checkout.

Checks that the installed package is the expected version, installed next
to sequentialist-core in the shared namespace, that its command is
registered, and that a trial is designed, takes records, is analysed,
reported on and keeps a verifiable audit log.
"""
import sys
from importlib.metadata import entry_points, version


def main(expected: str) -> int:
    import sequentialist
    import sequentialist.clinical as clinical
    from sequentialist.core import Engine, MemoryStorage, parse_spec, verify_audit

    assert clinical.__version__ == expected == version("sequentialist-clinical"), (clinical.__version__, expected)
    assert getattr(sequentialist, "__file__", None) is None, "the sequentialist namespace must stay open"
    found = entry_points()
    commands = found.select(group="sequentialist.commands") if hasattr(found, "select") else \
        found.get("sequentialist.commands", [])
    assert "clinical" in {p.name for p in commands}

    spec = clinical.design_trial(parse_spec({
        "kind": "clinical.trial/1", "name": "smoke", "endpoint": {"type": "time-to-event", "alternative": "less"},
        "arms": {"treatment": "A", "control": "B"},
        "design": {"alpha": 0.025, "analyses": [0.5, 1.0], "power": 0.9, "assumptions": {"hazard_ratio": 0.7}},
    }, clinical.TrialSpec))
    day = 86400.0
    records = []
    for i in range(80):
        records.append(clinical.patient("cut-1", len(records), f"p{i}", "A" if i % 2 else "B", i * day))
        records.append(clinical.event("cut-1", len(records), f"p{i}", (i + 20 + i % 13 * 9) * day, i % 3 != 1))
    trial = clinical.Trial.create(Engine(MemoryStorage(), clock="event"), "smoke", spec)
    trial.ingest(records)
    result = trial.analyse(150 * day, actor="statistician", reason="smoke test")
    assert result["analysis"] == 1 and result["information"] > 0
    assert "open part" in trial.report("open").to_markdown()
    assert trial.report("closed").findings()
    entries = [e.to_json() for e in trial.task.audit()]
    trial.close()
    assert verify_audit(entries, "smoke") == len(entries)
    print(f"sequentialist-clinical {expected}: installed, a trial designed and analysed "
          f"({result['information']} events), {len(entries)} audit entries verified")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
