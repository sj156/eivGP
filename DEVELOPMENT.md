# Development experiments — eivmixgp 0.2.1

Development estimates are scientifically informative, but replication
uncertainty and diagnostic warnings must accompany their interpretation.

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
library(eivmixgp)
settings <- eivgp_run_settings(core_budget = 12L, pending_datasets = 3L)
settings
# Inside each dataset worker:
# fit <- do.call(fit_eivgp, c(list(X = X, y = y, C = C,
#   U_obs = U_obs, engine = "multivariate", latent_dim = 2L),
#   settings$fit_args))
# Optional continuation:
# fit <- continue_eivgp(fit, n_iter = 1000L)
```

This helper returns settings, not a running experiment. It does not change
low-level sampler defaults. `MIXEDGP_DEV_DRAWS` overrides retained iterations
per chain, excluding warmup; `MIXEDGP_DEV_BURN` overrides warmup;
`MIXEDGP_DEV_REPS` controls replications.

## Scope and remaining work

This release implements iteration accounting and core-budgeted concurrency,
not the complete earlier development brief. Current cells remain Study I
eta0/eta1 with calibration 5/20 and Study II primary q2/q4 with calibration
12 and 6/24, three replications and 200 test observations.

Remaining work includes the ten-cell design and affine control, scheduling
across settings, method-level time caps/recovery, fit/evaluation cache
separation, and uniform nonfatal diagnostic handling. Publication still has
its existing stricter gates; modes do not yet differ only in experiment scale.
No claim is made that the full brief's acceptance tests pass. Development
runs can still take substantial time and stop on failures.

No output is automatically copied into Overleaf or pushed to GitHub.
The older `experiments/run_publication_study.R` is the legacy eivGP launcher,
not the eivmixgp development entry point.
