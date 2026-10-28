# Changelog

## sequentialist-clinical 0.1.0 — unreleased

First version, on sequentialist-core 0.1.0 and sequential-analysis 0.4.

- Specs `clinical.trial/1`: a binary, continuous or time-to-event endpoint
  in two arms; group sequential designs (the default: alpha spending of
  the Lan-DeMets O'Brien-Fleming, Lan-DeMets Pocock or Kim-DeMets power
  family, non-binding futility boundaries) or anytime-valid tests; roles
  of the people who may run analyses, recommend and decide.
- Designs computed once (`sequentialist clinical design`) and stored in the
  spec: boundaries, maximum patients or events, power, expected
  information; `check` recomputes them.
- Trials on the engine of sequentialist-core: facts about patients by
  pseudonym, later facts replacing earlier ones; analyses at data
  cut-offs, boundaries at the information actually reached
  (`GroupSequentialMonitor`); recommendations of the data monitoring
  committee and decisions of the sponsor in the audit log.
- Reports: the plan; the report to the committee with an open part (pooled,
  small counts hidden) and a closed part (by arm, every finding with its
  guarantee); the protocol of an analysis with its evidence.
- Records from CDISC ADaM (ADSL with ADTTE or a BDS dataset), from CSV and
  Dataset-JSON 1.1.
- Simulations of the type I error and the power of the whole pipeline.
- Checked against scipy (multivariate normal boundaries, chi-square,
  t-test), lifelines (log-rank) and closed formulas (maximum information).
