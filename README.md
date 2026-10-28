# sequentialist-clinical

Interim monitoring of clinical trials: group sequential designs with alpha
spending (the default) or anytime-valid tests, analyses at data cut-offs,
reports to the data monitoring committee with an open and a closed part,
and an audit trail of every analysis, recommendation and decision. The
statistics are those of
[sequential-analysis](https://github.com/andre-a-no/sequential-analysis);
every conclusion names its guarantee (`exact`, `asymptotic`, `none`) as the
procedure behind it states it.

Part of the [sequentialist](https://github.com/andre-a-no/sequentialist-core)
family; built on
[sequentialist-core](https://github.com/andre-a-no/sequentialist-core) (specs,
the checkpointed engine, storage, audit log, reports).

```
pip install sequentialist-clinical        # Python 3.8+ (not on PyPI yet)
```

Until the packages are on PyPI, install from the repositories:

```
pip install "sequentialist-core @ git+https://github.com/andre-a-no/sequentialist-core@27f71847ea39cedf0a8e41a5f05042428a533443#subdirectory=packages/core"
pip install "sequentialist-clinical @ git+https://github.com/andre-a-no/sequentialist-clinical"
```

This is a tool for the trial statistician. It is not a validated system
for regulatory submissions, and it does not claim compliance with 21 CFR
Part 11, ICH E9 or GCP.

## A trial

A spec (`clinical.trial/1`) names the endpoint, the arms and the design:

```json
{
  "kind": "clinical.trial/1",
  "name": "rx-301",
  "endpoint": {"type": "time-to-event", "alternative": "less"},
  "arms": {"treatment": "A", "control": "B"},
  "design": {"alpha": 0.025, "analyses": [0.5, 0.75, 1.0], "spending": "obrien-fleming",
             "futility": "obrien-fleming", "power": 0.9, "assumptions": {"hazard_ratio": 0.7}},
  "roles": {"unblinded": ["statistician-1"], "dsmb": ["dsmb-chair"], "sponsor": ["sponsor-1"]}
}
```

| endpoint | group sequential (default) | anytime |
|---|---|---|
| `binary` | `GroupSequentialProportions` statistic: asymptotic | `TwoProportionTest`: exact |
| `continuous` | `GroupSequentialTwoMeans` statistic (Z with a known sigma, Student's t otherwise) | — (none in sequential-analysis 0.4) |
| `time-to-event` | `GroupSequentialLogRank` statistic, stratified if `strata` are given | `AnytimeLogRankTest`: exact, with a confidence sequence for the hazard ratio |

`sequentialist clinical design trial.json -o trial.designed.json` computes
the design once: the efficacy (and non-binding futility) boundaries at the
planned analyses, the maximum number of patients with an outcome or of
events (that of rpact and gsDesign: the fixed-sample size, by Schoenfeld's
formula for events, times the inflation factor), the power and the
expected information at stopping. The trial runs from the designed spec,
so the plan in the audit trail is the plan that was executed;
`sequentialist clinical check` recomputes it and lists what differs.

## Records and analyses

Records are facts about patients, identified by pseudonyms
(`sequentialist.core.pseudonymize`); a name or a document number with
spaces is refused, and so is any field the format does not know:

```python
from sequentialist.clinical import Trial, TrialSpec, event, patient
from sequentialist.core import Engine, SQLiteStorage, load_spec

trial = Trial.create(Engine(SQLiteStorage("trials.db")), "rx-301", load_spec("trial.designed.json", TrialSpec))
trial.ingest([patient("cut-1", 0, "p-0001", "A", "2026-02-03T00:00:00Z"),
              event("cut-1", 1, "p-0001", "2026-09-30T00:00:00Z", False)])   # last contact, no event
result = trial.analyse("2026-10-01T00:00:00Z", actor="statistician-1", reason="planned analysis 1")
trial.recommend(1, "continue", actor="dsmb-chair", reason="no boundary crossed")
trial.decide(1, True, actor="sponsor-1", reason="recommendation followed")
print(trial.report("open").to_markdown())        # pooled over the arms
print(trial.report("closed").to_markdown())      # by arm, with the statistics
print(trial.protocol(1).to_markdown())           # the data cut, the result, its evidence
```

A later record about the same fact replaces the earlier one: new
follow-up, a new data cut and a correction arrive the same way, and the
order of reporting does not matter. An analysis takes the data as of its
cut-off: patients with an outcome observed by then, or, for the time to
an event, every patient who entered by then on the study's time scale,
with the event if it happened by the cut-off and censored at the last
contact or at the cut-off otherwise. Patients without any follow-up
record are left out and counted.

For a group sequential design the standardised statistic and the
information reached (patients, or events) go to `GroupSequentialMonitor`,
which computes the boundary of this analysis by alpha spending at the
information actually reached (Lan and DeMets): analyses need not happen
when planned. The analysis is final when the information reaches the
maximum or when it is called final. Crossing a boundary is recorded in the
audit log with the evidence of both streams (procedure, guarantee, history
digest of the exact data); it stops nothing by itself: the committee
recommends, the sponsor decides, each with actor and reason.

`examples/interim_monitoring.py` runs a time-to-event trial on synthetic
data through three data cuts.

## CDISC ADaM

`sequentialist.clinical.adam` reads ADSL and an analysis dataset from CSV
or CDISC Dataset-JSON 1.1 (JSON or NDJSON) with the standard library, and
turns them into records: `time_to_event_records` (ADTTE: `STARTDT`, `ADT`,
`CNSR` of one `PARAMCD`) and `outcome_records` (any BDS dataset: `AVAL` at
`ADT`). `USUBJID` is replaced by its pseudonym under the data owner's key
before anything is kept. Only these variables are read; SAS transport
files (XPT) are not. The tests read the ADSL and ADTTE of the CDISC pilot
study (`tests/data/cdiscpilot01`, from CDISC's Dataset-JSON repository,
MIT license).

## How it is checked

* Boundaries: against an independent computation of the multivariate
  normal probabilities (scipy, Genz's algorithm) for O'Brien-Fleming,
  Pocock, power-family and two-sided designs, and against the published
  Lan-DeMets O'Brien-Fleming values of gsDesign.
* Maximum information: the closed formulas of the fixed-sample tests
  (Schoenfeld, two means, two proportions) for one analysis.
* The statistic at a cut-off: the log-rank test against lifelines on data
  without ties (on the CDISC pilot, which has ties, within 0.5 %:
  sequential-analysis takes tied events one at a time, lifelines and R's
  survdiff use the hypergeometric variance); two proportions against
  scipy's chi-square test; two means against scipy's t-test.
* The whole pipeline: simulated type I error and power with Monte Carlo
  intervals (`sequentialist clinical simulate`).
* The engine: a trial that crashes at any step and is reopened is the
  uninterrupted trial, state and audit log alike, on every storage; the
  open report contains no arm code and no statistic.

The reference values are produced by `tools/make_references.py` and stored
in `tests/data/references.json`.

## Limits of 0.1

* Two arms; no adaptive randomisation, sample size re-estimation or
  conditional power.
* Effect estimates after a group sequential trial are biased; the closed
  report gives rates, means or events by arm as descriptions (guarantee
  "none"), not adjusted estimates, and no hazard ratio for the group
  sequential log-rank test.
* With an estimated standard deviation the t statistic of two means is
  compared with boundaries computed for Z: approximate for small samples.
* Tied event times are taken one at a time (see above).
* The open and closed parts of a report are separated by the package; who
  can read the storage and the reports is up to the environment. The spec
  says which arm code is the treatment.
* The whole state of a trial (every patient) is written at every
  checkpoint, as sequentialist-core 0.1 does.
* Drug safety surveillance (MaxSPRT) and quality of care (risk-adjusted
  CUSUM) wait for their procedures in sequential-analysis.

## License

GNU AGPL v3 or a commercial license ([COMMERCIAL.md](COMMERCIAL.md)).
Contributions need the [Contributor License Agreement](CLA.md); see
[CONTRIBUTING.md](CONTRIBUTING.md).

Contact: andrei@hobukob.ru, a.hobukob@gmail.com.
