# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Novikov Laboratories LLC
# Commercial licenses for use outside the AGPL: see COMMERCIAL.md
"""Design the group sequential example trials: X.trial.json -> X.designed.json; with --check, compare."""
import json
import sys
from pathlib import Path

from sequentialist.clinical import TrialSpec, check_trial, design_trial
from sequentialist.core import load_spec

HERE = Path(__file__).resolve().parent


def main() -> int:
    check = "--check" in sys.argv[1:]
    failed = 0
    for source in sorted(HERE.glob("*.trial.json")):
        spec = load_spec(source, TrialSpec)
        if not spec.group_sequential:
            continue
        target = source.with_name(source.name.replace(".trial.json", ".designed.json"))
        if check:
            problems = check_trial(load_spec(target, TrialSpec)) if target.exists() else ["missing"]
            if not problems and load_spec(target, TrialSpec).data["design"] != spec.data["design"]:
                problems = ["the design part differs from " + source.name]
            for p in problems:
                print(f"{target.name}: {p}")
            failed += bool(problems)
            continue
        target.write_text(json.dumps(design_trial(spec).to_json(), indent=2, ensure_ascii=False) + "\n",
                          encoding="utf-8")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
