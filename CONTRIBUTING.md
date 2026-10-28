# Contributing

Thank you for your interest in sequentialist.

## Contributor License Agreement

sequentialist is offered under the GNU AGPL v3 and under commercial licenses
(see [COMMERCIAL.md](COMMERCIAL.md)). To keep that possible, Novikov
Laboratories LLC must hold the rights needed to license every part of the
code under both. A contribution can therefore be merged only after its
author has signed the Contributor License Agreement ([CLA.md](CLA.md),
Russian and English; a draft until it has been reviewed by a lawyer), sent
to andrei@hobukob.ru. Small fixes of typos in documentation are welcome
without it.

If you write code as part of your job, your employer may own it; in that
case the employer has to sign the CLA as well.

## Development

```
python -m pip install "sequentialist-core @ git+https://github.com/andre-a-no/sequentialist-core@27f71847ea39cedf0a8e41a5f05042428a533443#subdirectory=packages/core"
python -m pip install -e . ruff mypy
python -m unittest discover -s tests -t .
ruff check . && mypy
```

- sequentialist-clinical runs on Python 3.8 and later and needs nothing but
  sequentialist-core and sequential-analysis (CI checks both, `vermin`
  checks the language level).
- Statistics belong in sequential-analysis; this package runs procedures,
  it does not compute tests.
- Every result is checked against an external reference or a closed
  formula: `tools/make_references.py` (scipy, lifelines) writes the
  reference values the tests compare with.
- New source files start with the SPDX header used throughout the project.
