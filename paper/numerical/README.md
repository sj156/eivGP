# Numerical paper replication

## Study II

Use [Study II reproduction](study2-README.md). The active files are:

- `study2-eivGP.R`: complete EIV-GP paper workload, including CUT and independent-test scoring.
- `study2-competitors.R`: all original and complete-case competitors, independently.
- `study2-common.R`: thin lane selection; inference remains in the package.
- `study2.R`, `study2.Rmd`, `study2.html`: saved-result reporting.

The combined temporary Study II launcher is retired. No pilot or integration
sensitivity run is needed merely to reproduce the current paper. Use separate
machine-local output roots and the same source, settings, and frozen data.

## Study I

Install with `Rscript tools/setup_numerical.R`. The complete workflow is
`Rscript paper/numerical/study1.R --action=plan --profile=paper`; select
`--component=core` for the original study or `--component=complete-case` for the
published complete-case adaptations. `--action=run` executes deliberately;
`--action=report` reads saved outputs. Notebook: `study1.Rmd`.

The complete paper workflow includes all 600 calibrated competitor fits with
four measurement chains, 12000/4000 total/warmup iterations, 16 outer states,
128 latent draws for means and 64 for response scoring. A separate follow-up
launcher is unnecessary for fresh replication. Historical saved scores remain
available for reproducing the current manuscript without new fitting.

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
