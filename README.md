# eivGP

`eivGP` implements errors-in-variables Gaussian process regression with
numeric predictors, ordinal proxies for latent continuous inputs, and sparse
calibration measurements. This repository is both the literate package source
and the reproducibility companion for the numerical studies.

## Install

Install the generated package subdirectory from GitHub:

```r
pak::pak("sj156/eivGP/eivGP", dependencies = TRUE)
```

From a checkout, regenerate and check the package with:

```r
litr::render("index.Rmd")
devtools::check("eivGP", document = FALSE)
```

## Stable package interface

```r
fit <- eivGP::fit_eivgp(
  X, y, C, U_obs,
  engine = "multivariate",
  latent_dim = 2L,
  kernel = "matern"
)

mean_draws <- eivGP::predict_eivgp(
  fit, new_X = X_new, new_C = C_new, target = "mean"
)
response_draws <- predict(
  fit, new_X = X_new, new_C = C_new, target = "response"
)
latent_draws <- eivGP::impute_eivgp(fit, new_C = C_new)
```

The published mixed-input competitors are available through
`fit_ucgp()`, `fit_lvgp()`, `fit_ezgp()`, and
`fit_mixedgp_competitor()` when their suggested packages are installed.

## Package organization

The installed package is organized by reusable responsibility, rather than by
paper study:

- `model_api.R` provides the unified user-facing fit, prediction, imputation,
  summary, and diagnostic interface.
- `model_univariate.R` implements the one-latent-input model with arbitrary
  numeric predictors.
- `model_multivariate.R` implements multivariate latent inputs and multiple
  ordinal proxies.
- `core_numerics.R` contains the single shared implementation of linear
  algebra, MCMC diagnostic, and predictive scoring helpers.
- `core_parallel.R` contains deterministic parallel and adaptive-MCMC control
  utilities.
- `application_interface.R` validates general real-data inputs.
- `competitors.R` contains optional external-method adapters.
- `reproduction_data.R`, `reproduction_workflows.R`, and
  `reproduction_compat.R` isolate paper-reproduction concerns from the model
  API.

Study I and Study II names therefore appear only where they describe frozen
paper designs, data generators, output tables, or compatibility entry points.
They are not separate packages or separate public model implementations.

## Reproduce the paper's numerical experiments

This section is the supported route for reproducing the numerical tables and
figures in the paper. Real-data analyses are intentionally not included yet.

### 1. Obtain the reproducibility companion and install its dependencies

Clone this repository, then run the following from the repository root:

```sh
Rscript --vanilla experiments/install_eivgp_dependencies.R
R CMD INSTALL eivGP
```

This installs eivGP's imported packages (`posterior`, `TruncatedNormal`) and
the packages needed to render the numerical-study reports. The second command
installs the local `eivGP` package. Repeat `R CMD INSTALL eivGP` after pulling
changes.

The documents use `library(eivGP)`, as an external user would. The package
installs the frozen simulation source used by its public functions. Optional
published-competitor packages are reported separately by the preflight and do
not block EIV-GP.

### 2. Generate frozen synthetic data once

Generate each study's deterministic input data before fitting. Open
[`experiments/study1_synthetic_data.Rmd`](experiments/study1_synthetic_data.Rmd)
or [`experiments/study2_synthetic_data.Rmd`](experiments/study2_synthetic_data.Rmd),
set `run_mode: "publication"`, and render it. Each creates the revised frozen
design below the chosen reproduction root and refuses to overwrite existing files.

From the repository root, the safest data-only command uses the built-in
publication and `reproduction/` defaults:

```sh
Rscript --vanilla experiments/run_publication_study.R study1-data
```

Use `study2-data` for Study II. The generated manifests and checksums are the
fixed inputs for every later fitting run. Set `EIVGP_ARTIFACT_ROOT` only when
you intentionally want outputs outside `reproduction/`.

### 3. Common manual publication procedure

Use this sequence for any ready study after its frozen inputs have been
generated. Run every command from the repository root. Do not regenerate data
when reproducing an existing study.

```sh
R CMD INSTALL eivGP
```

Select the study-specific input and output directories:

```sh
STUDY=study1
DATA_ROOT="reproduction/data/synthetic/$STUDY"
RESULT_ROOT="reproduction/results/$STUDY"
export STUDY DATA_ROOT RESULT_ROOT
```

For Study II, change the first line to `STUDY=study2`. Confirm that the frozen
cell manifests are present:

```sh
find "$DATA_ROOT" -name manifest.rds
```

Run the read-only preflight. It must show every design cell with
`complete_count TRUE` and all required packages as available:

```sh
Rscript --vanilla -e 'library(eivGP); study <- Sys.getenv("STUDY"); data_root <- Sys.getenv("DATA_ROOT"); result_root <- Sys.getenv("RESULT_ROOT"); config <- if (identical(study, "study1")) study1_simulation_config(mode = "dry_run", workers = 8L, data_root = data_root, output_root = result_root) else study2_simulation_config(mode = "dry_run", workers = 8L, data_root = data_root, output_root = result_root); if (identical(study, "study1")) run_study1_simulation(config) else run_study2_simulation(config)'
```

Start the formal fit and aggregation. With the default repository-local output
location, this avoids shell line-continuation mistakes:

```sh
Rscript --vanilla experiments/run_publication_study.R "$STUDY"
```

Keep the terminal open until the command exits successfully. If the machine
does not have enough memory for eight workers, restart the final command as
`EIVGP_WORKERS=4 Rscript --vanilla experiments/run_publication_study.R "$STUDY"`.
Do not treat smoke-mode output as publication output.

The real-data application will follow this same pattern—freeze an audited data
bundle, run a preflight, fit through `eivGP`, then write results and a report.
Its data bundle and run wrapper are not released yet, so
[`applications/real_data_application.Rmd`](applications/real_data_application.Rmd)
is intentionally not a publication command at this stage.

### Monitor a long-running fit

The runner writes a single, atomically updated progress record at
`config/progress.csv` inside its fingerprinted run directory. From a second
terminal, locate the newest Study I run and print its current phase:

```sh
RUN_DIR="$(find reproduction/results/study1 -maxdepth 1 -type d -name 'study1-publication-*' | sort | tail -n 1)"
cat "$RUN_DIR/config/progress.csv"
```

For a live terminal monitor that refreshes every 30 seconds:

```sh
while true; do clear; date; cat "$RUN_DIR/config/progress.csv"; sleep 30; done
```

`config/cell_status.csv` records completed cells and elapsed times;
`config/diagnostic_gates.csv` records their quality gates. Stop the monitor
with `Ctrl-C`; this does not stop the experiment.

After the first design cell finishes, `progress.csv` also contains elapsed
wall-clock time, an estimated remaining time, and an estimated finish time.
This estimate is a completed-cell average, so it is informative rather than a
guarantee: cells and optional competitor fits can have different costs.

For both studies, each EIV-GP fit also writes sampler-level ETA updates every
1,000 iterations under the cell's `mcmc_progress/` directory. For example,
monitor the first Study I fit with:

```sh
tail -f "$RUN_DIR/cells/eta0_balanced/results/study1_publication/mcmc_progress/rep001_calib005.log"
```

### Fast code-path demonstration (not paper reproduction)

Readers who want to inspect the complete workflow quickly should use smoke
mode. It uses separate frozen smoke data, one replication per design cell, and
short chains; it exercises the data, fitting, aggregation, table, and figure
code but must not be compared with or cited as the paper's results.

```sh
EIVGP_RUN_MODE=smoke Rscript --vanilla experiments/run_publication_study.R study1-data
EIVGP_RUN_MODE=smoke Rscript --vanilla experiments/run_publication_study.R study1
```

Substitute `study2-data` and `study2` for the Study II demonstration. Formal
publication mode retains the full frozen replication grid and strict MCMC
diagnostic gates, and can require substantially longer runtimes.

### Shared adaptive MCMC protocol and mixing diagnostics

Study I and Study II use the same continuation protocol for the main EIV-GP
posterior sampler. Warmup is separate. Each chain first retains 5,000
post-warmup draws, then the fit evaluates effective sample size for the key
model parameters. With four chains, each chain targets one fourth of the
prespecified total ESS of 200. If the target is missed, the required chain
length is projected from the observed ESS rate, rounded up to the next 1,000
retained draws, and capped at 15,000 retained draws per chain. The chain is
continued from its current parameter state and RNG state: warmup and completed
draws are not restarted. A nondefault cap can be supplied through the
simulation configuration or `adaptive_control` in `fit_eivgp()`.

The simulation drivers parallelize independent frozen data replications.
Chains within a replication are serial, preventing nested process trees and
oversubscription. The final diagnostic summary records actual iterations,
retained draws, minimum key-parameter ESS, R-hat, elapsed time, extension
count, termination reason, and parallel backend. The detailed decision history
is written as `study1_mcmc_schedule.csv` or
`study2_mcmc_schedule_<fingerprint>.csv`; sampler ETA logs are stored in
`mcmc_progress/`.

The measurement-model ablation samplers remain on their separately declared
fixed schedules; the adaptive protocol governs the main EIV-GP sampler used
for the reported Study I and Study II results.

### Study I mixing plots

Study I no longer runs a separate mandatory publication pilot. The former
pilot duplicated costly fits and gave Study I a different execution path from
Study II. Both studies now start their frozen replication grids directly and
apply the shared adaptive ESS protocol to every main EIV-GP fit.

The Study I publication sampler uses four chains, 1,000 warmup iterations and
the shared 5,000-to-15,000 retained-draw adaptive schedule, a thorough
latent-update schedule, joint elliptical-slice updates for correlated GP
hyperparameters, and centered/noncentered interweaving for ordinal thresholds
and latent inputs. The transition kernels preserve the specified posterior.
Fits that still miss the predeclared ESS or R-hat gate at the cap are retained
with an explicit diagnostic warning. The simulation continues through later
replications and design cells, and `config/diagnostic_gates.csv` provides the
audit record that must be reviewed before treating the aggregate as final.
External-competitor fit failures follow the same nonfatal, auditable policy.
Actual cell execution errors still stop immediately.

For the first actual Study I replication at the smallest and largest anchored
calibration sizes, the runner writes three mixing plots under
`cells/<design-cell>/figures/`:

- `*_trace.pdf` overlays post-warmup traces from all chains;
- `*_acf.pdf` compares within-chain autocorrelation;
- `*_rank.pdf` compares pooled-rank histograms across chains.

The same plots can be created for an independently fitted univariate model:

```r
diagnostics <- eivGP::plot_eivgp_mcmc_diagnostics(fit)
diagnostics$trace
diagnostics$autocorrelation
diagnostics$rank
```

### 4. Run another numerical study

Use the portable launcher from the repository root. Substitute any writable
directory for `/path/to/output-directory`; it is the only location that will
receive generated data, simulation results, tables, figures, and HTML reports.

```sh
EIVGP_ARTIFACT_ROOT=/path/to/output-directory \
  Rscript experiments/run_publication_study.R study1
```

Use `study2` in place of `study1` to reproduce Study II. The command defaults
to the frozen, reportable `publication` design. It verifies frozen inputs and
then performs these stages in order:

1. Runs the prespecified `eivGP` and competitor analyses.
2. Validates the paired design and the exact competitor implementations,
   aggregates results, and writes the paper-ready result artifacts.
3. Renders an HTML audit report.

Every run is fingerprinted. Re-run the same command after an interruption;
the frozen inputs are checksum-verified before work begins.

### 5. Locate and inspect the results

For Study I, the command writes these files beneath the chosen reproduction
directory:

```
data/synthetic/study1/                  # frozen synthetic datasets and manifest
results/study1/<fingerprinted-run>/
  cells/<design-cell>/figures/           # paper figures for each design cell
  cells/<design-cell>/tables/            # CSV and LaTeX tables for each design cell
  combined/                             # aggregated raw results and paired comparisons
  config/                               # validation records, hashes, versions, session details
reports/study1_numerical_experiment.html
```

Study II uses the analogous `study2` directories. The HTML report lists the
exact aggregated files produced for that run. Keep the full fingerprinted run
directory with any manuscript result: it is the audit trail linking tables and
figures to input data, code hashes, package versions, diagnostics, and
validation outcomes.

### Dry runs and computing controls

Before committing compute time, verify the local setup without modifying
outputs:

```sh
EIVGP_RUN_MODE=dry_run \
Rscript experiments/run_publication_study.R study1
```

For a small, explicitly non-reportable end-to-end test, use
`EIVGP_RUN_MODE=smoke`. Generate its separate smoke inputs before fitting:

```sh
EIVGP_ARTIFACT_ROOT=/path/to/output-directory EIVGP_RUN_MODE=smoke \
  Rscript experiments/run_publication_study.R study1-data
EIVGP_ARTIFACT_ROOT=/path/to/output-directory EIVGP_RUN_MODE=smoke \
  Rscript experiments/run_publication_study.R study1
```

Smoke data are stored under `data/smoke/` and are never paper results. Other
optional controls are:

- `EIVGP_WORKERS=8` sets the number of replication workers.
- `MIXEDGP_R_LIBRARY=/path/to/R-library` adds a custom R library before the
  preflight, useful on shared systems.
- `EIVGP_ARTIFACT_ROOT` is optional; without it, reproducibility files go to
  `reproduction/` under the repository checkout.

### Numerical linear algebra and macOS performance

The EIV-GP sampler repeatedly factors GP covariance matrices, so an optimized
BLAS/LAPACK library can improve runtime. Its benefit must be benchmarked on the
target design: the numerical studies parallelize independent replications, and
their relatively small per-fit matrices may benefit less than a large dense
matrix workload.

Check the numerical library used by the current R installation:

```sh
Rscript --vanilla -e 'cat("BLAS:  ", La_library(), "\n", sep = "")'
R_HOME="$(R RHOME)"
otool -L "$R_HOME/lib/libRblas.dylib"
```

On macOS, output that names `libRblas.dylib` without `Accelerate` or
`openblas` normally indicates R's bundled reference BLAS. Environment variables
such as `OPENBLAS_NUM_THREADS` select the number of threads only; they do not
switch R to an OpenBLAS backend.

For a safe comparison, install a **separate** R build linked to one backend,
then reinstall `eivGP` and its dependencies in that R library. R's official
administration manual documents the relevant configure options:

- Apple Accelerate: `--with-newAccelerate=lapack` on supported macOS/SDK
  combinations.
- OpenBLAS: install a compatible OpenBLAS library and configure R with an
  explicit `--with-blas` value, for example
  `--with-blas="-L/path/to/openblas/lib -lopenblas"`.

Do not replace `libRblas.dylib` inside a working R installation in place: it
can invalidate compiled packages and makes a reproducibility environment hard
to audit. Benchmark a separate R installation first, and keep the selected R
version and `La_library()` output with the run metadata.

The publication runners deliberately set `OMP_NUM_THREADS`,
`OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, and `VECLIB_MAXIMUM_THREADS` to
`1` while multiple replication workers are active. This prevents, for example,
eight R workers each starting eight BLAS threads. For these parallel numerical
studies, start by setting `EIVGP_WORKERS` to the desired number of independent
workers and leave BLAS single-threaded; evaluate a different allocation only
by benchmark.

No personal filesystem path is embedded in these commands. The standalone
documents—[`experiments/study1_numerical_experiment.Rmd`](experiments/study1_numerical_experiment.Rmd)
and [`experiments/study2_numerical_experiment.Rmd`](experiments/study2_numerical_experiment.Rmd)—may also be rendered from RStudio after setting their `run_mode` parameter to `publication`.

### Real-data application

The real-data workflow will be released separately after its inputs and
analysis specification have been audited. Its placeholder document,
[`applications/real_data_application.Rmd`](applications/real_data_application.Rmd),
is not part of the current paper-reproduction instructions.

## Source organization

The explicitly ordered root `.Rmd` chapters are the authoritative `litr`
package source. Their newly written and maintained R payload is under
`script/newly-written-and-modified-scripts/`; historical R code is isolated
under `script/legacy-code-scripts/`. Historical literate chapters are retained
under `legacy/` for provenance and excluded from the active package build.
