# Numerical paper replication

## Study I

[Study I reproduction](study1-README.md) follows the same separate-machine and
reporting standards as Study II. `study1-eivGP.R` covers eivGP;
`study1-competitors.R` covers original and calibrated competitors. Default action
is plan; use `--action=run` explicitly. The former follow-up comparisons are
included automatically. `study1.R` remains a combined compatibility entry point.
`study1.Rmd` is report-only. Both studies use `convergence-report.R` / `.Rmd` for
failure rates and matched diagnostic-screened accuracy; keep Study I regimes separate.

## Shared conventions

`common.R` checks installed source and locked dependencies. `--cores` is the
total CPU budget; set BLAS libraries to one thread. Data and fits default to
`~/eivGP-results`, with `--data-root` and `--output` overrides. Changes to source,
priors, data or numerical budgets require compatible provenance or fresh roots;
never relabel historical caches to bypass checks. Scores retain failures and
diagnostic warnings. Report generation does not fit missing models.

`check-numerical.R` and `check-portability.R` are developer checks, not paper
replication steps; some checks fit small models and should not be run when fitting
is prohibited. The same package supplies all reusable modeling and reporting.

## Paper protocol (2 November 2026)

Use the separate study/method lane scripts above. They preserve all reported
primary, calibration and pressure-test cells while excluding superseded PI/CC
component ablations. Study II primary fitting is fixed-budget; its secondary CUT
comparison uses the reported 28/14 pairs. FULL/CC comparisons retain 50 datasets.
The scripts default to 12 cores with four chains per fit. See each study guide
for the distinct primary and 512-state evaluation protocols and archive selection.
