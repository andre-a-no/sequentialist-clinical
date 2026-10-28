# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Novikov Laboratories LLC
# Commercial licenses for use outside the AGPL: see COMMERCIAL.md
"""Write schemas/clinical.trial-1.schema.json from the code; with --check, fail if the file differs."""
import json
import sys
from pathlib import Path

from sequentialist.clinical import TrialSpec
from sequentialist.core import json_schema

ROOT = Path(__file__).resolve().parent.parent
TARGET = ROOT / "schemas" / "clinical.trial-1.schema.json"


def main() -> int:
    text = json.dumps(json_schema(TrialSpec), indent=2, ensure_ascii=False) + "\n"
    if "--check" in sys.argv[1:]:
        if not TARGET.exists() or TARGET.read_text(encoding="utf-8") != text:
            print(f"{TARGET.relative_to(ROOT)} is not what the code produces: run tools/make_schemas.py")
            return 1
        return 0
    TARGET.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
