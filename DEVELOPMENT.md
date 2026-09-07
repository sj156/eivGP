# Repository development experiments with eivGP 0.3.0

This document describes the separate numerical-experiment layer. Its study
cells, replication counts, evaluation budgets, and output paths are not
package defaults or installed package APIs. The current posterior engines
come from the canonical `codes/` modules used to build `eivGP` 0.3.0.
Existing older fits and caches require refitting under the new posterior.

Development estimates are scientifically informative, but replication
uncertainty and diagnostic warnings must accompany their interpretation.

## Dataset sizes and replications

Both development and publication configurations use 100 training and 100 test
observations per dataset in both studies. A replication is a newly simulated dataset under the same
setting, with a distinct prespecified seed; all applicable methods share it.
Development uses three replications in all seven settings. Publication uses
100 replications in every setting.

Previously frozen larger test datasets are not overwritten or silently truncated.
Their manifests may be incompatible with this revised design; select a new
data directory or explicitly prepare and audit the revised data before running.
Equal test sizes alone do not yet implement shared frozen data across modes.

## Core budget and iterations

Enter a core budget, not a dataset-worker count. On macOS/Linux, 12 cores
allow 3 datasets × 4 parallel chains; 16 cores allow 4 × 4. Workers are
capped by replications in the current setting: the current three-replication
profile uses at most three dataset workers even on a 16-core machine.
Settings are still processed sequentially. Windows falls back to serial.

Each chain runs **500 warmup + 1,250 sampling = 1,750 total iterations**.
Four chains retain all **5,000 post-warmup draws**, without thinning.
Continuation remains explicit via `continue_eivgp()`, never automatic.

From the GitHub repository root (use 16 for the MacBook):

```sh
export MIXEDGP_CORE_BUDGET=12
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
Rscript --vanilla experiments/run_development_study.R study1 plan
Rscript --vanilla experiments/run_development_study.R study2 plan
```

After reviewing settings, replace `plan` with `run`. Numerical-library
thread controls should be set before R starts; support depends on the linked
library. These controls do not switch the BLAS implementation.
The visible wrapper is `experiments/development_numerical_experiment.Rmd`.
Outputs stay under `reproduction/development/`, separate from publication.

Both `codes/run_study1_simulation.R` and `codes/run_study2_simulation.R`
also accept `MIXEDGP_RUN_MODE=development` and `MIXEDGP_CORE_BUDGET`.
A core budget takes precedence over the old worker setting. Without one,
the masters retain historical allocation.

## Package interface

```r
library(eivGP)
settings <- eivgp_run_settings(core_budget = 12L, pending_datasets = 3L)
settings
# Inside each dataset worker:
# fit <- do.call(fit_eivgp, c(list(X = X, y = y, C = C,
#   U_obs = U_obs, engine = "multivariate", latent_dim = 2L),
#   settings$fit_args))
# Optional continuation:
# fit <- continue_eivgp(fit, n_iter = 1000L)
```

This reusable helper returns an allocation and fitting arguments. The
environment variables below are read by repository experiment scripts, not
the package model API. `MIXEDGP_DEV_DRAWS` overrides retained iterations
per chain, excluding warmup; `MIXEDGP_DEV_BURN` overrides warmup;
`MIXEDGP_DEV_REPS` controls replications.

## Scope and remaining work

Fitting/evaluation and reporting now have separate entry points. Use the
development launcher's `fit` action to save cell-level `report_inputs.rds`
without figures or publication summaries. `run` performs reporting only after
the fitting stage. `experiments/report_study.R RUN_DIRECTORY summarize|plot|report`
regenerates outputs from saved data, without any model fitting. It accepts
older Study II raw bundles as well as the new checkpoints; missing cells are
listed explicitly. Reporting writes to a separate `reporting/` tree and does
not modify replication caches. See README.md for copy-ready commands.

The former monolithic drivers are split into `setup_study*_experiment.R`
(configuration and helpers), `02_study*_monte_carlo.R` (fitting/evaluation),
and `report_study*_results.R` (tables/figures). `experiment_reporting.R` restores
the reporting context, records provenance, and handles empty plots. Reporting
source files are tracked in provenance but excluded from the fit fingerprint,
so editing a figure does not by itself invalidate a fit. This is a repository
experiment-layer change, not a change to the installed posterior sampler.

The repository experiment layer retains iteration accounting and core-budgeted
concurrency. Study I crosses balanced/imbalanced categories with eta=0/1,
each with calibration 0/20/50. Study II uses primary q=2 (calibration 50),
primary q=4 (0/20/50/80), and logistic misspecification q=4 (50).
See NUMERICAL_DESIGN.md for the agreed design.

Remaining work includes shared frozen-data storage between modes, scheduling
across settings, method-level time caps/recovery, fit/evaluation cache
separation, and recovery from genuine fitting/evaluation exceptions.
Diagnostic handling is now shared: convergence and completeness checks
produce flags and warnings in both modes, never gate-triggered stops.
Missing competitors are recorded and skipped by default. Competitor fitting is
now a standalone experiment stage: run `experiments/install_eivgp_dependencies.R`
first, then `experiments/run_competitors.R study1 plan publication` (or `study2`).
Use `run` to prepare missing fits, `retry` to retry failed fits, and `export
development` to extract the first three publication datasets without fitting.
Set `EIVGP_OVERLEAF_ROOT` to the existing paper folder for run/export actions.

`codes/competitor_cache.R` owns the shared optimization protocol and per-method
cache. The two Monte Carlo drivers consume it read-only; neither contains an
independent competitor optimizer configuration. Calibration and mode are not
cache keys because these competitors use the same observed training data and
no calibration measurements. Inputs, method settings, seeds, package/R identity
and adapter source are checked before reuse. The cache is experiment-only and
does not change the installed MCMC package or migrate old simulation bundles.

Regression checks (run from repository root):

```sh
Rscript --vanilla codes/tests/test_competitor_cache.R
Rscript --vanilla codes/tests/test_competitor_export.R
Rscript --vanilla codes/tests/test_reporting_separation.R
```

The export test uses fixture fits and temporary data/paper directories, not
publication results. Real full-grid competitor runs must be launched explicitly.
No claim is made that the full brief's acceptance tests pass. Development
runs can still take substantial time and stop on failures.

## Interpreting diagnostic flags

Main MCMC tables retain pass/fail fields and add `diagnostic_warning` and
`diagnostic_advice`. Flagged fitted objects are saved before proceeding.
Study II measurement warnings no longer suppress PI-GP/CC-GP calculations;
their status and measurement diagnostic tables identify the warning.
Run-level `config/diagnostic_gates.csv` records unresolved checks and advice.

Low ESS or high MCSE with otherwise stable mixing can justify more sampling.
High R-hat, separated traces or unidentified coordinates need investigation;
longer runs are not a guaranteed remedy. Use `continue_eivgp()` only on a
compatible public `eivgp_fit`; raw simulation-engine fit bundles are not
automatically converted to that interface. Continuation remains opt-in.

Flags are not permission to present unreliable numbers as established results.
Retain them alongside estimates, and report missing-method and uncertainty
information. Corrupt input data, missing core runtime dependencies, disk errors
and unhandled fitting/evaluation exceptions still stop rather than invent outputs.
Changed cache schemas prevent silently reusing old gate-suppressed bundles;
existing cache files are preserved, but a new run can require recomputation.

No output is automatically copied into Overleaf or pushed to GitHub.
`experiments/run_publication_study.R` resolves its study constructors from the
same repository helper bundle. Package release validation does not run either
development or publication studies.
