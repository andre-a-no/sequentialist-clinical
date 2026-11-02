# Repository guidance

Read [CLAUDE.md](CLAUDE.md) for the existing owner instructions and
[CONTRIBUTING.md](CONTRIBUTING.md) for contribution rules. Preserve them.

## Purpose and boundaries

`sequentialist-clinical` is application software for interim monitoring of
clinical trials. The Python import is `sequentialist.clinical`; the CLI is
`sequentialist clinical`. Reusable statistics belong to `sequential-analysis`;
durable task execution belongs to `sequentialist-core`.

Current code handles two-arm trials, group sequential or supported anytime-valid
procedures, data cut-offs, ADaM imports and open/closed committee reports.
Treatment selection and patient-level diagnostic search are not implemented.
Do not claim regulatory validation, ICH E9/GCP/21 CFR Part 11 compliance,
adjusted post-stopping estimates or unsupported error guarantees.

## Source map

- [README.md](README.md): user workflow, endpoint support and limits.
- [docs/overview.md](docs/overview.md): plain-language overview and FAQ.
- [src/sequentialist/clinical/spec.py](src/sequentialist/clinical/spec.py): trial specs.
- [src/sequentialist/clinical/design.py](src/sequentialist/clinical/design.py): design integration.
- [src/sequentialist/clinical/trial.py](src/sequentialist/clinical/trial.py): analyses, recommendations and decisions.
- [src/sequentialist/clinical/adam.py](src/sequentialist/clinical/adam.py): ADaM import.
- [src/sequentialist/clinical/report.py](src/sequentialist/clinical/report.py): reports.
- [tests](tests): behavior and reference checks.
- [tools/make_references.py](tools/make_references.py): reference-value generation.

## Validation

Run commands from this repository after installing compatible core and analysis
dependencies. Existing installation and CI details are in CONTRIBUTING and
[the test workflow](.github/workflows/tests.yml).

```bash
python -m unittest discover -s tests -t .
ruff check .
mypy
vermin --no-tips --violations -t=3.8- src tests examples tools
```

For changes to generated schemas or example designs:

```bash
python tools/make_schemas.py --check
python examples/make_trials.py --check
```

For documentation-only changes, verify referenced files, links, actual command
names and current-versus-planned claims; `git diff --check` is an appropriate
format check. Do not regenerate references or run simulations solely to validate
unchanged explanatory prose. Report tests actually executed.

## Assistant entry points

Codex-compatible clients can read this file. Claude's existing context is in
CLAUDE.md; Qwen and Gemini context files point here. Loading these files depends
on the client. Grok, Perplexity and other retrieval tools can use the README,
overview and [llms.txt](llms.txt); no universal automatic ingestion is assumed.
Keep documentation factual, source-linked and readable without programming.
