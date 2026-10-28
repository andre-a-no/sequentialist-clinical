# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Novikov Laboratories LLC
# Commercial licenses for use outside the AGPL: see COMMERCIAL.md
"""Records of a trial from CDISC ADaM datasets: ADSL with ADTTE, or ADSL with a dataset of outcomes.

An interim analysis usually starts from the analysis datasets of a data
cut.  `read_dataset` reads one, from CSV or from CDISC Dataset-JSON 1.1
(``.json``, or ``.ndjson`` with the metadata on the first line), with the
standard library only.  SAS transport files (XPT) are not read.

`time_to_event_records` and `outcome_records` turn a subject-level dataset
and an analysis dataset into the records of a `Trial`:

* the patient: ``USUBJID``, replaced by its pseudonym under the data
  owner's key (`sequentialist.core.pseudonymize`) before anything is kept;
  the arm (``TRT01P``, planned treatment, by default) and the start of
  follow-up; only patients of the population (``ITTFL`` = "Y" by default);
* time to event (ADTTE, one parameter ``PARAMCD``): the start ``STARTDT``,
  the date of the event or censoring ``ADT`` and ``CNSR`` (0: event);
* an outcome (any BDS dataset, one parameter): ``AVAL`` at ``ADT``.

Only these variables are read; the package does not claim conformance to
the ADaM standard.  Dates are ISO 8601 (Dataset-JSON, most CSV exports),
taken at midnight UTC, or numbers of days since 1960-01-01 (SAS dates).

Every record of a data cut gets the source ``<cut>/<dataset>`` and its row
number, so ingesting the same cut twice changes nothing and a later cut
replaces what an earlier one said about the same patient.
"""
from __future__ import annotations

import csv
import datetime as _dt
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

from sequentialist.core import Record, pseudonymize

from .trial import event, outcome, patient

_DATE = re.compile(r"(\d{4})-(\d{2})-(\d{2})(?:[Tt ](\d{2}):(\d{2})(?::(\d{2})(?:\.\d+)?)?)?")
_SAS_EPOCH = _dt.datetime(1960, 1, 1, tzinfo=_dt.timezone.utc)


@dataclass(frozen=True)
class Dataset:
    """An analysis dataset: its name and its rows as {variable: value}."""

    name: str
    columns: Sequence[str]
    rows: Sequence[Mapping[str, Any]] = field(default_factory=list)

    def where(self, **equal: Any) -> list[Mapping[str, Any]]:
        return [r for r in self.rows if all(r.get(k) == v for k, v in equal.items())]


def read_dataset(path: str | Path, name: str | None = None) -> Dataset:
    """A dataset from CSV, Dataset-JSON (.json) or Dataset-JSON NDJSON (.ndjson)."""
    p = Path(path)
    suffix = p.suffix.lower()
    if suffix == ".csv":
        with p.open(newline="", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            rows = [dict(row) for row in reader]
            columns = list(reader.fieldnames or [])
        return Dataset(name or p.stem.upper(), columns, rows)
    if suffix == ".ndjson":
        lines = [line for line in p.read_text(encoding="utf-8").splitlines() if line.strip()]
        meta = json.loads(lines[0])
        return _dataset_json(meta, [json.loads(line) for line in lines[1:]], name)
    if suffix == ".json":
        meta = json.loads(p.read_text(encoding="utf-8"))
        return _dataset_json(meta, meta.get("rows", []), name)
    raise ValueError(f"{p}: read .csv, .json or .ndjson (Dataset-JSON); SAS transport files are not read")


def _dataset_json(meta: Mapping[str, Any], rows: Sequence[Sequence[Any]], name: str | None) -> Dataset:
    version = str(meta.get("datasetJSONVersion", ""))
    if not version.startswith("1.1"):
        raise ValueError(f"Dataset-JSON {version or '(no version)'}: version 1.1 is read")
    columns = [c["name"] for c in meta["columns"]]
    out = []
    for i, row in enumerate(rows):
        if len(row) != len(columns):
            raise ValueError(f"row {i} has {len(row)} values for {len(columns)} columns")
        out.append(dict(zip(columns, row)))
    return Dataset(name or str(meta.get("name", "")), columns, out)


def date_seconds(value: Any, what: str = "date") -> float:
    """Seconds since the epoch of an ISO 8601 date or datetime (UTC) or of a SAS date (days since 1960)."""
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return (_SAS_EPOCH + _dt.timedelta(days=float(value))).timestamp()
    text = str(value).strip()
    if re.fullmatch(r"-?\d+(\.\d+)?", text):
        return date_seconds(float(text), what)
    match = _DATE.fullmatch(text)
    if match is None:
        raise ValueError(f"{what}: {value!r} is not an ISO 8601 date or a SAS date")
    y, m, d, hh, mm, ss = (int(g) if g else 0 for g in match.groups())
    return _dt.datetime(y, m, d, hh, mm, ss, tzinfo=_dt.timezone.utc).timestamp()


def _number(value: Any, what: str) -> float:
    if isinstance(value, bool) or value is None or value == "":
        raise ValueError(f"{what}: missing or not a number ({value!r})")
    try:
        return float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{what}: {value!r} is not a number") from None


class _Patients:
    """The population of ADSL: pseudonym, arm and stratum by USUBJID."""

    def __init__(self, adsl: Dataset, key: bytes | str, arm: str, population: str | None, stratum: str | None):
        self.by_id: dict[str, dict[str, Any]] = {}
        for row in adsl.rows:
            if population is not None and row.get(population) != "Y":
                continue
            usubjid = row["USUBJID"]
            self.by_id[usubjid] = {"pid": pseudonymize(str(usubjid), key), "arm": str(row[arm]),
                                   "stratum": None if stratum is None else str(row[stratum]), "row": row}


def _numbered() -> Any:
    """Dense sequence numbers for one source, as the engine's intake expects."""
    counter = iter(range(2**53))
    return lambda: next(counter)


def time_to_event_records(adsl: Dataset, adtte: Dataset, *, paramcd: str, key: bytes | str, cut: str,
                          arm: str = "TRT01P", population: str | None = "ITTFL",
                          stratum: str | None = None) -> list[Record]:
    """Patient and event records of the population from ADSL and one parameter of ADTTE."""
    people = _Patients(adsl, key, arm, population, stratum)
    rows = adtte.where(PARAMCD=paramcd)
    if not rows:
        raise ValueError(f"{adtte.name}: no rows with PARAMCD = {paramcd!r}")
    source = f"{cut}/{adtte.name}"
    seq = _numbered()
    out = []
    for row in rows:
        who = people.by_id.get(row["USUBJID"])
        if who is None:
            continue
        where = f"{adtte.name} {row['USUBJID']}"
        start = date_seconds(row["STARTDT"], f"{where} STARTDT")
        out.append(patient(source, seq(), who["pid"], who["arm"], start, who["stratum"]))
        censored = _number(row["CNSR"], f"{where} CNSR")
        out.append(event(source, seq(), who["pid"], date_seconds(row["ADT"], f"{where} ADT"), censored == 0))
    return out


def outcome_records(adsl: Dataset, bds: Dataset, *, paramcd: str, key: bytes | str, cut: str,
                    arm: str = "TRT01P", population: str | None = "ITTFL", stratum: str | None = None,
                    entry: str = "TRTSDT", value: str = "AVAL", date: str = "ADT") -> list[Record]:
    """Patient records from ADSL and outcome records from one parameter of a BDS dataset."""
    people = _Patients(adsl, key, arm, population, stratum)
    out = []
    seq = _numbered()
    for usubjid, who in sorted(people.by_id.items()):
        start = date_seconds(who["row"][entry], f"{adsl.name} {usubjid} {entry}")
        out.append(patient(f"{cut}/{adsl.name}", seq(), who["pid"], who["arm"], start, who["stratum"]))
    source = f"{cut}/{bds.name}"
    seq = _numbered()
    rows = bds.where(PARAMCD=paramcd)
    if not rows:
        raise ValueError(f"{bds.name}: no rows with PARAMCD = {paramcd!r}")
    for row in rows:
        found = people.by_id.get(row["USUBJID"])
        if found is None:
            continue
        where = f"{bds.name} {row['USUBJID']}"
        out.append(outcome(source, seq(), found["pid"], _number(row[value], f"{where} {value}"),
                           date_seconds(row[date], f"{where} {date}")))
    return out
