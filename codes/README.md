# Mixed-input EIV-GP computational workflow

September 6 status: optional exact no-HMC mixing improvements and stricter
diagnostic precision/window checks are implemented. Package checks pass, but
the small pilots still fail the convergence gates; package success is not
publication readiness. Transition-kernel defaults are unchanged; every
post-warm-up draw is now retained by the public package interface.
See `../COMPUTATION-REVISION-20260906.md` and run the bounded comparison:

```sh
Rscript 17_no_hmc_mixing_pilot.R
```

That script compares both studies, freezes its datasets and source snapshots,
and prints its output directory. Its header lists workstation size, iteration,
kernel, calibration, parallel-core, and output-directory controls. New moves
are available through `fit_eivgp()` and the low-level engines; the publication
masters continue to use the unchanged baseline until a strategy is validated
for the intended scenarios. No HMC is used.

For the publication simulations, start with
[`PUBLICATION_SIMULATIONS.md`](PUBLICATION_SIMULATIONS.md). The current
publication interface is one helper bundle (`simulation_helpers.R`) and one
master script per study (`run_study1_simulation.R` and
`run_study2_simulation.R`). Both masters default to a read-only dry run.

The posterior engines are shared with the planned R package, while publication
orchestration is intentionally narrower. Run the following from this directory:

```sh
Rscript run_study1_simulation.R
Rscript run_study2_simulation.R
```

Those commands are dry runs. Generate frozen data, smoke-test the complete
workflow, and launch the archival runs by setting `MIXEDGP_RUN_MODE` to
`data`, `smoke`, and `publication`, respectively. For example:

```sh
MIXEDGP_RUN_MODE=publication MIXEDGP_WORKERS=8 \
  Rscript run_study1_simulation.R
```

Set `MIXEDGP_R_LIBRARY=/absolute/path/to/R-library` when the audited package
library is not stored beside the project root. Optional `MIXEDGP_DATA_ROOT`
and `MIXEDGP_OUTPUT_ROOT` variables place large workstation artifacts outside
the synchronized Overleaf tree without changing the masters.

`simulation_helpers.R` is the single source of truth for design cells, seeds,
task eligibility, frozen-data checks, dispatch, gates, aggregation, and paired
comparisons. The masters contain only the user-editable run configuration and
one function call. The scientific posterior engines remain separate because
they are also the package source of truth; they are not duplicated in the
simulation layer.

## Canonical modules

- `00_public_api.R`: package-facing EIV-GP and competitor calls.
- `00_mcmc_workflow.R`: public `diagnose_eivgp()` and `continue_eivgp()`.
- `00_diagnostics.R`: shared full-window, target-aware diagnostic computation.
- `simulation_helpers.R`: publication design and orchestration bundle.
- `run_study1_simulation.R` and `run_study2_simulation.R`: publication masters.
- `00_experiment_runner.R`: package-facing and historical unified workflow.
- `00_parallel_utils.R`: deterministic `mclapply`/`mcmapply` wrappers.
- `00_synthetic_data.R`: versioned frozen datasets and checksum manifests.
- `00_study1_functions.R` and `00_study2_functions.R`: posterior engines.
- `03_study2_published_competitors.R`: audited public-package adapters shared
  by both studies; the historical filename is retained for compatibility.

The `run_study1_all.R` and `run_study2_all.R` files are compatibility aliases
for the audited masters. `run_mixedgp_experiment()` remains available for
package examples and historical workflows, but it is not the publication
simulation entry point.

## Install, diagnose, and extend

The installable release is `../package-release/eivmixgp_0.2.0.tar.gz`.
From the project root, install dependencies if needed, then the local archive:

```r
install.packages(c("posterior", "TruncatedNormal")) # once, if not installed
install.packages("package-release/eivmixgp_0.2.0.tar.gz", repos = NULL, type = "source")
library(eivmixgp)
```

After fitting at least four chains with `fit_eivgp()`, use a diagnostic panel
covering the scientific comparisons you intend to report:

```r
report <- diagnose_eivgp(fit, X = X_panel, C = C_panel)
report
report$table
saveRDS(fit, "fit-checkpoint.rds")
# Only after reviewing whether additional sampling is appropriate:
fit <- continue_eivgp(readRDS("fit-checkpoint.rds"), n_iter = 2000L)
report <- diagnose_eivgp(fit, X = X_panel, C = C_panel)
```

The extension adds 2,000 transitions/draws per chain, without reinitialization
or repeated warm-up. Terminal states, RNG streams, update schedules and frozen
tuning survive save/reload. Use the same package version and numerical
environment. Historical fits without checkpoints must be refitted; the package
cannot recover discarded draws from thinned fits. Checkpoints are returned at
the end of a fitting call, not periodically during an interrupted call.

`poor_exploration` calls for trace/chain/model investigation; simply extending
may not help. `insufficient_precision` means agreement passes the R-hat screen
but ESS/mean-MCSE criteria are not met. `insufficient_diagnostics` cannot issue
a pass. `screen_passed` applies only to the monitored functionals and is not a
convergence or identification guarantee. Supply `U` at explicit panel locations
for calibrated surface diagnostics, and `additional_series` for important
contrasts. Quantile MCSE and integration sensitivity require separate checks.
Default target panels use all retained draws with common integration randomness.
Increasing `n_latent`/changing the integration seed checks error not included
in the MCMC MCSE. No thinning is used or accepted by the public interface.

## Reproducibility gates

Before a publication run:

1. run each master without environment overrides and inspect the dry-run
   specification and package preflight;
2. run `08_published_competitor_validation.R` and require every prespecified
   competitor for the archival run;
3. run `09_experiment_design_validation.R`;
4. run both masters with `MIXEDGP_RUN_MODE=smoke` as end-to-end tests;
5. run both masters with `MIXEDGP_RUN_MODE=publication`; publication mode
   activates strict competitor and MCMC gates;
6. archive manifests, raw RDS results, diagnostics, package versions, and
   session information together.

Generated datasets live under `../data-synthetic/publication-v2/`. Their
CSV/RDS manifests record schema versions, seeds, sizes, and MD5 checksums.
Before fitting, the masters also verify the common-random-number and nested-
proxy contracts. Fitting stops when a required artifact is absent or altered;
it never silently regenerates data in a fit-only run.
