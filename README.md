# eivGP

Gaussian process regression with ordinal proxies and partially observed continuous
inputs. This repository provides the package and the paper's numerical studies,
Ocean application and ADNI application.

## Repository structure

```text
eivGP/R/               Shared model, sampler, prediction and reporting utilities
  helpers.R            Data, kernels, prediction, storage and workflow helpers
  mcmc.R               Fitting, continuation and sampler diagnostics
  hmc.R                Joint HMC coordinates, adaptation and transitions
  plotting.R           Shared summaries, figures and report utilities
eivGP/src/             Analytic gradients and leapfrog integration
tools/                 Application scripts, directly in this directory
  ocean-tools.R        All Ocean-specific functions
  adni-tools.R         All ADNI-specific functions
  ocean-eivGP.R, ocean-competitors.R, ocean-report.R
  adni-eivGP.R, adni-competitors.R, adni-report.R
  ocean-prepare-folds.R, adni-prepare-folds.R
paper/numerical/        Study I/II entrypoints, Rmd/HTML and reproduction guides
paper/applications/     Ocean/ADNI Rmd, HTML and reproduction guides
data/                  Supplied inputs, metadata and private-data placeholders
tests/application/     Application tests and workflow checks
```

The package contains reusable model and numerical implementations. Each
`tools/*-tools.R` contains its application's preparation, orchestration,
evaluation and reporting. Entry scripts load that file and select a workload.
Competitor adapters use `kergp`, `LVGP` and `EzGP`.

## Install and configure

Use R >= 4.1, a working package compiler and Pandoc. From the repository root,
install dependencies, refresh source provenance and install this checkout:

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "knitr", "patchwork", "posterior",
                   "jsonlite", "TruncatedNormal", "digest", "rmarkdown", "testthat",
                   "Rcpp", "RcppArmadillo", "kergp", "LVGP", "EzGP"))
```

```sh
Rscript tools/build_package.R
R CMD INSTALL eivGP
```

For numerical paper replication, use `Rscript tools/setup_numerical.R` to restore
`renv.lock` and install the matching checkout into `.r-library/`. Application
scripts also find `.r-library/`; set `EIVGP_R_LIBRARY` for a custom library.
Restart R after updating a package. Do not replace a library used by active jobs.

Set numerical-library threads to one **before starting R**:

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
```

`--cores=N` is a shared CPU budget for chains, prediction, diagnostics and figure
export, capped at detected CPU capacity. It does not alter chains, seeds or
numerical settings. Waiting coordinators release slots; resident fits have a
separate memory limit. macOS/Linux use forked workers; Windows uses the serial
fallback. HTML assembly remains serial.

## Model interfaces and samplers

Both engines support ESS, HMC and the HMC–ESS hybrid:

```r
fit <- fit_eivgp(X, y, C, U_obs=U_obs, engine="univariate", sampler="ESS")
fit_hmc <- fit_eivgp(X, y, C, U_obs=U_obs, engine="univariate", sampler="HMC")
fit_hybrid <- fit_eivgp(X, y, C, U_obs=U_obs, engine="univariate",
  sampler="HYBRID", sampler_control=list(hybrid_every=20L))
continued <- continue_eivgp(fit_hmc, n_iter=1000L)
# Multivariate inputs: engine="multivariate", latent_dim=2L.
```

ESS is the package default. HMC jointly updates missing inputs, observed-input
kernel coefficients and signal fraction; auxiliary ESS updates remain in the
algorithm. HMC supports conditional dictionary updates and SE or Matérn 3/2 and
5/2 kernels. Use ESS for other smoothness or marginal dictionary updates;
unsupported HMC/HYBRID requests raise an error. HYBRID adds an ESS sweep after
every 20 HMC iterations by default and retains one draw per combined transition.

Supply tuning through `sampler_control=list(hmc=list(step=.05,
leapfrog=c(8L,24L),target_accept=.8))`. Adaptation occurs during warmup only.
Continuation preserves the sampler, metric, step size and RNG state; it does not
repeat warmup or switch samplers. Check scientific-target R-hat, bulk/tail ESS,
MCSE and traces. ESS efficiency alone does not establish posterior mixing.
See the [package guide](eivGP/README.md) and function documentation for details.

## Reproduce the paper experiments

Choose an experiment below, then follow its numbered preparation, fitting and
reporting steps. Commands run from the repository root and use 16 cores as an
example; choose a budget appropriate for your machine. New runs require matching
source, dependencies and frozen inputs. Existing fits keep their original
provenance; use new output roots for incompatible sources or numerical settings.

- [1. Synthetic experiments: Studies I and II](#1-synthetic-experiments-studies-i-and-ii)
- [2. Ocean application](#2-ocean-application)
- [3. ADNI application](#3-adni-application)

### 1. Synthetic experiments: Studies I and II

**1.1. Prepare the environment and inputs.** Run `Rscript tools/setup_numerical.R`.
The study scripts generate fixed-seed datasets and masks locally. EIV-GP and
competitor lanes can run on separate machines with the same source and data
protocol. Inspect work with `--action=plan` before replacing it with `run`.

**1.2. Fit and analyze EIV-GP.** These commands include each study's complete
paper workload and calibrated comparisons:

```sh
Rscript paper/numerical/study1-eivGP.R --action=run --cores=16 --output=results/study1-eivGP --data-root=results/data/study1
Rscript paper/numerical/study2-eivGP.R --action=run --cores=16 --output=results/study2-eivGP --data-root=results/data/study2
```

Study II reuses calibrated four-proxy fits for independent-test FULL scoring and
includes the CUT feedback comparison. FULL is rescored, not fitted twice.

**1.3. Fit and analyze competitors.** Each lane includes original UC-GP, LVGP and
EzGP and their calibration-aware complete-case versions:

```sh
Rscript paper/numerical/study1-competitors.R --action=run --cores=16 --output=results/study1-competitors --data-root=results/data/study1
Rscript paper/numerical/study2-competitors.R --action=run --cores=16 --output=results/study2-competitors --data-root=results/data/study2
```

The competitor machine does not require EIV-GP response fits. Numerical scripts
default to planning; fitting requires `--action=run`. Repeat compatible commands
to resume saved work. Original grids and calibrated comparisons retain their
separate paper protocols and budgets.

**1.4. Read saved reports.** Study I EIV-GP runs automatically generate
`core/study1-report/study1.html`. Repeat a lane's command with `--action=report`
to read its saved results; the guides describe explicit run selection when
several runs exist. The `study1.Rmd` and `study2.Rmd` notebooks are report-only.
When collecting lanes from different machines, validate matching settings/data
and pair raw scores by replication; do not average machine-level summaries.

See [Study I reproduction](paper/numerical/study1-README.md),
[Study II reproduction](paper/numerical/study2-README.md), and the
[convergence/failure report](paper/numerical/convergence-report.Rmd).

### 2. Ocean application

**2.1. Select the inputs.** `ocean500` is the main paper cohort and the default.
Its frozen five-fold inputs are supplied in `data/ocean_silicate/n500/`.
`ocean200` is also supplied; prepare `ocean1000` separately with
`Rscript tools/ocean-prepare-folds.R ocean1000`.

**2.2. Fit and analyze EIV-GP.** Preview the workload, then run it:

```sh
Rscript tools/ocean-eivGP.R --plan=true --profile=paper --cores=16
Rscript tools/ocean-eivGP.R --profile=paper --cores=16
```

Paper sampling is fixed at four chains, initially 50,000 iterations per chain
including 5,000 burn-in. Optional `--ess-min=100` continues saved states until all
monitored scientific bulk and tail ESS exceed 100; its default is `NULL`.
For a separate explicit fixed-length extension, use
`tools/ocean.R --action=continue --profile=paper --dataset=ocean500 --continuation-iter=20000 --cores=16`.
Neither operation repeats warmup. Reporting never starts continuation.

**2.3. Fit and analyze competitors.** Run against the same input/output root:

```sh
Rscript tools/ocean-competitors.R --profile=paper --cores=16
```

This includes original categorical UC-GP/LVGP/EzGP and complete-case comparisons.
Rerunning a compatible command reuses successful tasks; failures and optimizer
nonconvergence remain visible.

**2.4. Assemble and inspect the paper report.**

```sh
Rscript tools/ocean-report.R --profile=paper --cores=16
```

The combined report reads saved fits and scores. `reports/paper/` contains the
manuscript figures, tables and `artifact_manifest.csv`; method reports are
`reports/eivGP/index.html` and `reports/competitors/index.html`. To rebuild method
reports without fitting, add `--reuse-only=true` to the corresponding method
command. Direct saved-result rendering through `paper/applications/ocean.Rmd`
is also available.

| Path | Contents |
| --- | --- |
| `data/ocean_silicate/n500/` | Frozen main cohort, five-fold inputs and masks |
| `results/ocean500/` | Paper fits, predictions, diagnostics and `ocean.html` |
| `results/ocean500_pilot/` | Separate short-chain workflow check |
| `paper/applications/ocean.html` | Versioned 200-profile example report |

New roots use `ocean-flat-v1`: `eiv_foldNN_O*/`, `baseline_foldNN/` and
`cc_foldNN_O*_METHOD/` are siblings, with `data/`, `logs/`, `reports/` and `_meta/`
for shared material. Existing unmarked roots keep their layout. Migration is
explicit and dry-run by default; see the
[Ocean reproduction and storage guide](paper/applications/README_ocean.md).

### 3. ADNI application

**3.1. Supply the private prepared inputs.** Place the flat fold files and `DATASET.json` in
`data/adni/` using the [data contract](data/adni/README.md). No ADNI
participant records are distributed in Git.

The default design uses ADNI composite memory, FDG-PET SUVR and hippocampal volume
normalized by intracranial volume; nine native study-partner ECog items
(MEMORY1–8 and DIVATT3); and continuous age/education within women. Coordinate
labels, dimensions, masks and seeds come from the supplied protocols. For the
three-item cohort, use `data/adni_3C/` and add `--dataset=3C` to all three
ADNI commands; its paper root is `results/3C/`.

**3.2. Fit and analyze EIV-GP.**

```sh
Rscript tools/adni-eivGP.R --plan=true --profile=paper --cores=16
Rscript tools/adni-eivGP.R --profile=paper --cores=16
```

Paper runs use five folds, four chains and initially 50,000 iterations per chain
including 5,000 burn-in. Defaults are HMC at O0/O20, ESS at O50, 250 prediction
states per chain and full competitor controls. Explicit numerical overrides are
recorded in the run identity. Optional `--ess-min=100` enables scientific-target
continuation after the final planned stage. O20 calibration is nested in O50
under the supplied C-stratified design.

**3.3. Fit and analyze competitors.** Run sequentially against the same root:

```sh
Rscript tools/adni-competitors.R --profile=paper --cores=16
```

This covers original ordinal-input UC-GP/LVGP/EzGP and matched complete-case
comparisons. `--fit-jobs` and `--evaluation-jobs` can limit resident fits
independently of the total shared CPU budget; see the ADNI guide.

**3.4. Assemble and inspect the paper report.**

```sh
Rscript tools/adni-report.R --profile=paper --cores=16
```

The report includes settings and C definitions, Y/U/X and C/U illustrations,
prediction, complete-case comparisons, scientific diagnostics, response functions
and U imputation. Method reports and paper exports use the same relative paths
as Ocean. `--reuse-only=true` rebuilds a method's saved report without fitting;
`adni-report.R` always reads saved results. The report source is
`paper/applications/adni.Rmd`.

| Path | Contents |
| --- | --- |
| `data/adni/` | Privately supplied metadata, five-fold inputs and masks |
| `results/adni/` | Paper fits, predictions, diagnostics and `adni.html` |
| `paper/applications/adni.html` | Versioned data-free plan page |

New roots use `adni-flat-v1`, matching Ocean's task folders and `data/`, `logs/`,
`reports/`, `_meta/` organization. Earlier stages stay inside their EIV task;
final fit, prediction and diagnostic files sit together. Existing unmarked roots
retain their original paths; select old folders explicitly with `--output`.
All participant inputs, frozen snapshots, fits and empirical reports stay local.
See the [ADNI reproduction and storage guide](paper/applications/README_adni.md).

## Check and interpret results

For both applications, `--output` overrides the result root and `--data-root`
overrides prepared inputs. Use the same dataset, profile and path flags for
fitting and reporting. Completion, optimizer status, scientific-target ESS,
R-hat, MCSE and traces are separate checks. ESS thresholds and a completed
workflow do not establish convergence. Exact historical numerical values require
the original saved fits and matching environment; new runs reproduce the design.

A synthetic ADNI check exercises the actual three CLI entrypoints, six comparator
scores, HTML, flat storage and saved-fit reuse without participant data:

```sh
Rscript tests/application/check-adni.R --cores=2
```

Additional application tests are under `tests/application/`. Run package tests
after installation:

```r
library(eivGP)
testthat::test_dir("eivGP/tests/testthat",env=new.env(parent=asNamespace("eivGP")),
                   stop_on_failure=TRUE)
```

After package-source edits, run `Rscript tools/build_package.R` and reinstall.
Keep application-specific functions in the corresponding `tools/*-tools.R` file;
shared model/numerical functions belong in the package. The
[paper overview](paper/README.md) links all entrypoints and report sources.
