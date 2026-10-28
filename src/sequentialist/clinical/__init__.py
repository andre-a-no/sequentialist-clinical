# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Novikov Laboratories LLC
# Commercial licenses for use outside the AGPL: see COMMERCIAL.md
"""sequentialist.clinical: interim monitoring of clinical trials.

spec       the spec of a trial (`TrialSpec`, kind "clinical.trial/1") and its records
design     group sequential designs, computed once (`design_trial`, `check_trial`)
trial      a trial as a task of the engine: facts about patients, analyses at data cut-offs (`Trial`)
report     the plan, the report of the data monitoring committee (open and closed parts), analysis protocols
adam       records from CDISC ADaM datasets (CSV, Dataset-JSON)
simulate   type I error and power of the whole pipeline, by simulation
"""
from .design import check_trial, design_trial
from .report import analysis_protocol, dsmb_report, plan_report
from .simulate import simulate_trial
from .spec import KIND, TrialSpec
from .trial import Trial, TrialHandler, event, outcome, patient
from .version import __version__

__all__ = [
    "KIND",
    "Trial",
    "TrialHandler",
    "TrialSpec",
    "__version__",
    "analysis_protocol",
    "check_trial",
    "design_trial",
    "dsmb_report",
    "event",
    "outcome",
    "patient",
    "plan_report",
    "simulate_trial",
]
