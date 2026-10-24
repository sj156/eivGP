# Unreleased

- Production Study I now continues exact saved states until all monitored scientific
  targets have bulk and tail ESS strictly above 100; retain R-hat and MCSE for review.
- Add `continue_eivgp_until_ess()` and `continue_study1_results()` to extend chains
  and refresh only affected scores, with fixed prediction RNG blocks and provenance.
- Remove redundant Study I checkpoint files only after complete target diagnostics,
  final fit histories, scores and reporting have been verified.

- Make Study I EIV-GP saved-result reporting independent of fitting/resume checks;
  the existing study1.Rmd renders a complete report with scientific-target ESS > 100.
- Save small target traces for all new Study I fits and reconstruct selected
  historical traces from verified saved inputs without new MCMC.
- Render mixing figures serially to avoid macOS Objective-C fork crashes.

# eivGP 0.4.0.9012

- Ocean retries only missing baseline methods, preserving successful fitted
  objects and the original method-specific seeds.
- LVGP optimizer starts share the application CPU pool through a scoped parallel
  dispatcher; the author package likelihood, starts and optimizer are unchanged.
- Ocean paper fits retain finite optimizer iteration limits but no longer abort
  each LVGP attempt after 900 seconds. Recovery provenance remains separate from
  historical EIV fit identities.

# eivGP 0.4.0.9011

- Ocean accepts a single `--cores` allocation for the full reproduction workflow,
  exports a manuscript-linked table/figure inventory, and retains per-method failures.
- Failed native comparators no longer abort valid EIV analysis. Partial comparator
  summaries retain fold counts; missing curves are not fabricated.
- Prediction and functional stages use the shared pool with separate resident-fit
  limits, and CPU-slot events remain available for concurrency audits.

# eivGP 0.4.0.9010

- Add the opt-in `sampler="HYBRID"`: joint HMC with one complete ESS sweep
  after every `hybrid_every` iterations (default 20), targeting the same full
  GP posterior. ESS remains the default.
- Validate hybrid frequency, expose truthful sampler metadata and HMC
  diagnostics, and preserve the schedule across exact continuation.
- Preserve existing ESS/HMC trajectories and legacy checkpoint continuation.

# eivGP 0.4.0

- Consolidate reusable implementation into helpers, MCMC, and plotting/reporting modules.
- Add package APIs for numerical studies, application workflows, summaries and reports.
- Add lossless compact fit storage with legacy reads and exact continuation.
- Retain the 0.3.1 posterior sampler and scientific settings.
- Replace duplicated repository implementation with four thin R/Rmd/HTML analyses.

## Publication recovery

- UC-GP and EzGP now support documented, bounded optimizer rescue sequences;
  original successful fits return immediately and attempt histories are retained.
- Added internal, refinement-checked two-dimensional Gauss--Hermite oracle mean
  integration for independent Gaussian/logistic ordinal scores. Numerical truth
  diagnostics distinguish quadrature error checks from Monte Carlo standard errors.
- Reference-task failures can be recorded without aborting fitted-method scoring.
  Reference calculations preserve the caller's RNG stream.
- Repository experiments add immediate fit/predictive checkpoints, replication-ID
  selection, cache-lock isolation, and a separate resumable publication recovery
  runner. They do not alter the posterior, MCMC budgets, or frozen datasets.

# eivGP 0.3.1

This patch changes transition kernels, not the 0.3.0 posterior, priors, finite
dictionary, calibration rules, or inferential targets. No thinning or HMC.

- Deterministic-threshold models now update cutoffs one at a time by exact
  truncated-Beta draws in cumulative-probability space. Every observed and
  currently imputed U constrains the interval. No GP factorization is needed.
- The existing cutoff/input transport is retained. The previous cutoff ESS
  step remains available as `sampler_control = list(threshold_update = "ess")`
  for controlled comparisons. The new default is `"gibbs"`.
- Ordinal-probit joint moves reuse the exact Schur-complement likelihood when
  enabled, with the existing dense fallback, and evaluate only the union of
  affected ordinal factors. This improves evaluation cost, not the transition's
  ideal-arithmetic mixing per iteration. No new probit transport is introduced.
- Chain statistics record `threshold_gibbs_updates`. Existing ESS counters
  still count slice updates only. Diagnostic thresholds are not weakened.
- Checkpoints explicitly identify sampler version 0.3.1. Use the matching
  version for continuation. 0.3.0 draws still target the same posterior, but
  continuing them under the old algorithm requires 0.3.0; start fresh to use
  the new transitions. Earlier-than-0.3.0 fits target a different posterior.

The package remains reusable independently of the manuscript's experiment
scripts. No experimental designs, datasets, cached results, or run budgets are
changed by this patch. Numerical-study caches must distinguish code revisions.
Small correctness/performance checks do not guarantee mixing in general.

- `fit_eivgp()` now accepts `ess_min = NULL` and `ess_panel = NULL`.
  A positive ESS threshold requests adaptive saved-state extensions after the
  initial budget until all scientific bulk and tail ESS values strictly exceed
  it. NULL preserves fixed-budget behavior. Scientific history and diagnostics
  are retained in `fit$scientific_ess`; raw-coordinate diagnostics remain separate.
- Study II core FULL fits enable threshold 100 in the paper profile, including
  compatible saved fits. Extensions precede scoring and retain checkpoint state.
  Two full-history scientific trace figures replace aggregate mixing dashboards
  when scientific trace artifacts are available. CUT's nested conditional workflow
  retains its existing sampling controls; this change applies to joint EIV-GP fits.


- Ocean and ADNI method CLIs accept optional `--ess-min=100` (default NULL).
  Final scientific checks and adaptive saved-state extensions precede scoring;
  ADNI applies the rule only at the final nominal stage and records actual
  per-fit budgets. Both save full-history scientific traces, with no fitting
  during report-only execution.
- CUT retains fixed budgets and all completed comparisons. New conditional-GP
  and outer-integration scientific traces are descriptive, with no ESS stopping
  requirement. Conditional diagnostics freeze the correct outer measurement state;
  outer diagnostics explicitly retain finite inner Monte Carlo uncertainty.
