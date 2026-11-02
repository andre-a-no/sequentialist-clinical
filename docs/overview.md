# Clinical trial monitoring with sequentialist-clinical

[English](overview.md) · [Русский](overview.ru.md)

`sequentialist-clinical` helps a trial statistician run and document interim
analyses of a two-arm clinical trial. It connects statistical procedures to a
declared trial plan, patient records, data cut-off dates and the decisions of
the data monitoring committee and sponsor.

In plain language: the package helps answer **what the accumulated trial data
show at this analysis, under the declared method and its assumptions**. A
crossed statistical boundary is recorded as evidence. It does not itself stop
the trial or administer a treatment.

## Who it is for

- Trial statisticians implementing a reproducible interim-analysis workflow.
- Data monitoring committees and sponsors who need understandable reports and
  a record of recommendations and decisions.
- Developers integrating trial analysis with a data platform.

It is not a patient-facing diagnostic tool or a treatment recommendation system.
Broader diagnosis and treatment-option search are possible future application
directions, not features of the current implementation.

## What goes in and what comes out

| Input | Meaning |
|---|---|
| Trial specification | Endpoint, treatment/control arms, analysis method, planned analyses and roles |
| Designed specification | The boundaries and maximum information computed for the plan |
| Pseudonymous patient records | Enrolment, outcomes or time-to-event follow-up/censoring data |
| Data cut-off | The date defining which observed information belongs to an analysis |
| Human commands | A committee recommendation or sponsor decision, with actor and reason |

Outputs include analysis results, boundary crossings, an analysis protocol,
open and closed committee reports, and an audit trail. Reports identify the
procedure and its declared guarantee. Descriptive group summaries after a
group sequential trial are not presented as adjusted effect estimates.

## A workflow without needing to read the implementation

1. A statistician defines the endpoint, arms, design and assumptions.
2. The design command computes the plan that will be executed.
3. Records are imported, using pseudonyms rather than direct patient identifiers.
4. An analysis is requested for a data cut-off.
5. The committee reviews the result and records a recommendation.
6. The sponsor records the decision. The software preserves the distinction
   between statistical evidence, recommendation and decision.

The [README](../README.md) contains a concrete trial specification and commands.
[The synthetic example](../examples/interim_monitoring.py) shows a time-to-event
trial at three data cut-offs.

## Current endpoint and data support

- Binary, continuous and time-to-event endpoints have the group sequential
  options listed in the README. Group sequential design is the default.
- The available anytime-valid options depend on the endpoint and upstream
  procedure. There is no general anytime-valid continuous-endpoint option in
  this package's documented 0.1 scope.
- ADaM ADSL plus supported analysis datasets can be read from CSV or CDISC
  Dataset-JSON 1.1, including NDJSON. XPT is not implemented.
- The package is application software. Statistical methods come from
  [sequential-analysis](https://github.com/andre-a-no/sequential-analysis),
  while storage, checkpoints and audit support come from
  [sequentialist-core](https://github.com/andre-a-no/sequentialist-core).

## Questions people commonly ask

### Does this find a treatment for an individual patient?

No. Current code monitors declared clinical trials. It does not choose a
patient's treatment, perform differential diagnosis or search a treatment graph.

### Does a boundary crossing automatically end the trial?

No. It is recorded. The committee recommends and the sponsor decides, each with
an actor and reason. The intended human workflow is explicit in
[trial.py](../src/sequentialist/clinical/trial.py).

### Can analyses happen at different information levels than planned?

For the group sequential workflow, the accumulated statistic and information
are passed to the upstream spending monitor, as described in the README. The
validity of the analysis still depends on the selected design and assumptions;
this is not a guarantee for arbitrary changes to an ongoing trial.

### Are the reports regulatory-ready or validated for submission?

The repository makes no claim of validation for regulatory submissions or
compliance with 21 CFR Part 11, ICH E9 or GCP. Open/closed report separation is
implemented, but access to storage and generated reports is the deployment's
responsibility.

### What are the main limitations?

Two arms; no adaptive randomisation, sample-size re-estimation or conditional
power. Group sequential post-stopping effect estimates are not adjusted. The
upstream treatment of tied event times differs from common grouped-tie
references. Each checkpoint writes the whole trial state. Drug safety
surveillance and quality-of-care procedures remain future work in this scope.
See [the full limits](../README.md#limits-of-01) before selecting a workflow.

### How is the implementation checked?

The repository includes boundary/reference comparisons, closed-form checks,
data cut-off tests, pipeline simulations and recovery/audit tests. The reference
generation script is [tools/make_references.py](../tools/make_references.py),
and tests and committed fixtures are in [tests](../tests). These checks do not
turn every output into an exact or universally applicable guarantee.

### Is this available as an installed package?

The repository contains a Python package, version 0.1.0 in
[version.py](../src/sequentialist/clinical/version.py), and declares Python 3.8+
in [pyproject.toml](../pyproject.toml). Follow the current checkout's installation
instructions and dependency metadata. A version in source does not establish
that a distribution has been published to a registry.

## Documentation for assistants and search systems

[llms.txt](../llms.txt) is a curated source index, and [AGENTS.md](../AGENTS.md)
contains development guidance. Claude, Codex, Qwen and Gemini clients have
different context-loading behavior; Grok, Perplexity and other retrieval tools
can use the public prose and linked sources when accessible. These files help
make claims attributable; they do not guarantee indexing or recommendation.
