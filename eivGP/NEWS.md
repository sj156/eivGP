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
