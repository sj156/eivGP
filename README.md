# eivGP

Gaussian process regression with ordinal proxies and partially observed continuous
inputs. This repository provides the package and the paper's numerical studies,
Ocean, ADNI and OAI applications.

Follow the numbered sections to locate the code, set up the environment, review
model settings, reproduce an experiment and check its results. Experiment groups
in Section 4 use letters A–D; each application has numbered preparation, fitting and reporting steps.

1. [Repository structure](#1-repository-structure)
2. [Install and configure](#2-install-and-configure)
3. [Model interfaces and samplers](#3-model-interfaces-and-samplers)
4. [Reproduce the paper experiments](#4-reproduce-the-paper-experiments)
5. [Check and interpret results](#5-check-and-interpret-results)

## 1. Repository structure

```
eivGP/R/
  helpers.R     # data preparation, kernels, prediction, storage, workflow helpers
  mcmc.R        # posterior sampling, fitting, continuation, sampler diagnostics
  hmc.R         # joint HMC coordinates, warmup adaptation and transitions
eivGP/src/      # analytic GP/ordinal gradients and leapfrog integration
  plotting.R    # summaries, comparison tables, figures, report generation
data/ocean_silicate/
  n200/         # frozen 200-profile data and five-fold inputs
  n500/         # frozen 500-profile data and five-fold inputs
  n1000/        # raw selected 1000-profile cohort
data/adni/
  inputs/       # local prepared ADNI folds; only .gitkeep is versioned
data/oai/       # private flat fold files; only contract and example are versioned
paper/
  numerical/    # study1 and study2: one .R, .Rmd, and .html each
  applications/ # Ocean, ADNI and OAI Rmd reports, data-free HTML and guides
  manifests/    # selected numerical runs and frozen Ocean input identities
tools/          # flat application scripts; no application subdirectories
  adni-tools.R  # all ADNI application functions
  ocean-tools.R # all Ocean application functions
  oai-tools.R   # all OAI application functions
  adni-eivGP.R, adni-competitors.R, adni-report.R
  ocean-eivGP.R, ocean-competitors.R, ocean-report.R
  oai-eivGP.R, oai-competitors.R, oai-report.R
  adni-prepare-folds.R, ocean-prepare-folds.R, oai-prepare-folds.R
  build_package.R
  render_reports.R
```

The package contains reusable model and numerical implementations. Each
`tools/*-tools.R` contains its application's preparation, orchestration,
evaluation and reporting. Entry scripts load that file and select a workload.
Competitor adapters use `kergp`, `LVGP` and `EzGP`.

## 2. Install and configure

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
export, capped at detected CPU capacity. For ADNI/Ocean/OAI, new EIV fits use min(N,4) chains; it does not alter seeds or
numerical settings. Waiting coordinators release slots; resident fits have a
separate memory limit. macOS/Linux use forked workers; Windows uses the serial
fallback. HTML assembly remains serial.

## 3. Model interfaces and samplers

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

## 4. Reproduce the paper experiments

Choose an experiment below, then follow its numbered preparation, fitting and
reporting steps. Commands run from the repository root and use 12 cores as an
example; choose a budget appropriate for your machine. New runs require matching
source, dependencies and frozen inputs. Existing fits keep their original
provenance; use new output roots for incompatible sources or numerical settings.

- [A. Synthetic experiments: Studies I and II](#a-synthetic-experiments-studies-i-and-ii)
- [B. Ocean application](#b-ocean-application)
- [C. ADNI application](#c-adni-application)
- [D. OAI application](#d-oai-application)

### A. Synthetic experiments: Studies I and II

**A1. Prepare the environment and inputs.** Run `Rscript tools/setup_numerical.R`.
The study scripts generate fixed-seed datasets and masks locally. EIV-GP and
competitor lanes can run on separate machines with the same source and data
protocol. Inspect work with `--action=plan` before replacing it with `run`. All method launchers accept `--action=plan|run|report`, `--cores`, `--data-root`, and `--output`. Application `--plan=true` and `--reuse-only=true` remain compatible aliases; do not combine them with `--action`. Combined report scripts accept `plan` or `report` only.

**A2. Fit and analyze EIV-GP.** These commands include each study's complete
paper workload and calibrated comparisons:

```sh
Rscript tools/study1-eivGP.R --action=run --cores=12 --output=results/study1-eivGP --data-root=results/data/study1
Rscript tools/study2-eivGP.R --action=run --cores=12 --output=results/study2-eivGP --data-root=results/data/study2
```

Study II reuses calibrated four-proxy fits for independent-test FULL scoring and
includes the CUT feedback comparison. FULL is rescored, not fitted twice.

**A3. Fit and analyze competitors.** Each lane includes original UC-GP, LVGP and
EzGP and their calibration-aware complete-case versions:

```sh
Rscript tools/study1-competitors.R --action=run --cores=12 --output=results/study1-competitors --data-root=results/data/study1
Rscript tools/study2-competitors.R --action=run --cores=12 --output=results/study2-competitors --data-root=results/data/study2
```

The competitor machine does not require EIV-GP response fits. Numerical scripts
default to planning; fitting requires `--action=run`. Repeat compatible commands
to resume saved work. Original grids and calibrated comparisons retain their
separate paper protocols and budgets.

**A4. Read saved reports.** Study I EIV-GP runs automatically generate
`core/study1-report/study1.html`. Repeat a lane's command with `--action=report`
to read its saved results; the guides describe explicit run selection when
several runs exist. The `study1.Rmd` and `study2.Rmd` notebooks are report-only.
When collecting lanes from different machines, validate matching settings/data
and pair raw scores by replication; do not average machine-level summaries.

See [Study I reproduction](paper/study1-README.md),
[Study II reproduction](paper/study2-README.md), and the
[convergence/failure report](paper/convergence-report.Rmd).

### B. Ocean application

**B1. Select the inputs.** `ocean500` is the main paper cohort and the default.
The supplied five-fold inputs in `data/ocean_silicate/n500/` are ready to run;
`ocean200` is also included. Prepare the raw `ocean1000` cohort once with
`Rscript tools/ocean-prepare-folds.R ocean1000`. See the
[Ocean guide](paper/README_ocean.md) for data preparation and sampling details.

**B2. Fit and analyze EIV-GP.** Preview the workload, then run it:

```sh
Rscript tools/ocean-eivGP.R --action=plan --profile=paper --cores=12 --data-root=data/ocean_silicate/n500 --output=results/ocean500
Rscript tools/ocean-eivGP.R --action=run --profile=paper --cores=12 --data-root=data/ocean_silicate/n500 --output=results/ocean500
```

Paper sampling uses min(cores,4) chains, initially 50,000 iterations per chain
including 5,000 burn-in. Optional `--ess-min=100` continues saved states until all
monitored scientific bulk and tail ESS exceed 100; its default is `NULL`.
For a separate explicit fixed-length extension, use
`tools/ocean.R --action=continue --profile=paper --dataset=ocean500 --continuation-iter=20000 --cores=12`.
Neither operation repeats warmup. Reporting never starts continuation.
The EIV-GP command summarizes its own predictions, oxygen imputation and
scientific-target diagnostics without requiring completed competitors.

**B3. Fit and analyze competitors.** Run against the same input/output root:

```sh
Rscript tools/ocean-competitors.R --action=run --profile=paper --cores=12 --data-root=data/ocean_silicate/n500 --output=results/ocean500
```

This includes original categorical UC-GP/LVGP/EzGP and complete-case comparisons.
The competitor command summarizes its own results without requiring completed
EIV-GP fits. The two method commands may run concurrently against the same result
directory; either method group may be started first. `--cores` is a budget per
command, so two commands requesting 12 cores request 24 workers in total.
Rerunning a compatible command reuses successful tasks; failures and optimizer
nonconvergence remain visible.

**B4. Assemble and inspect the paper report.**

```sh
Rscript tools/ocean-report.R --action=report --profile=paper --cores=12 --data-root=data/ocean_silicate/n500 --output=results/ocean500
```

After both method groups finish, the combined report reads saved fits and scores
and assembles their comparisons. `reports/paper/` contains the
manuscript figures, tables and `artifact_manifest.csv`; method reports are
`reports/eivGP/index.html` and `reports/competitors/index.html`. To rebuild method
reports without fitting, set `--action=report` on the corresponding method
command. Direct saved-result rendering through `paper/ocean.Rmd`
is also available.

**B5. Short workflow check.** Use the supplied `ocean200` cohort and a separate
output root. These commands run the same fitting, analysis and reporting
entrypoints with short chains:

```sh
Rscript tools/ocean-eivGP.R --action=run --profile=quick --dataset=ocean200 --cores=12 --output=results/ocean200_check --data-root=data/ocean_silicate/n200
Rscript tools/ocean-competitors.R --action=run --profile=quick --dataset=ocean200 --cores=12 --output=results/ocean200_check --data-root=data/ocean_silicate/n200
Rscript tools/ocean-report.R --action=report --profile=quick --dataset=ocean200 --cores=12 --output=results/ocean200_check --data-root=data/ocean_silicate/n200
```

Quick checks cover one fold. Use `--profile=pilot` with a separate output root
for five short-chain folds and the m/f diagnostic plots.

| Path | Contents |
| --- | --- |
| `data/ocean_silicate/n500/` | Frozen main cohort, five-fold inputs and masks |
| `results/ocean500/` | Fits, predictions, diagnostics and `ocean.html` |
| `results/ocean500_pilot/` | Separate short-chain workflow check |
| `paper/ocean.html` | Analysis-plan HTML for the supplied calibration design |

Ocean results use `ocean-flat-v1`: `eiv_foldNN_O*/`, `baseline_foldNN/` and
`cc_foldNN_O*_METHOD/` are siblings, with `data/`, `logs/`, `reports/` and `_meta/`
for shared material. See the
[Ocean reproduction and storage guide](paper/README_ocean.md).

### C. ADNI application

**C1. Supply the private prepared inputs.** Place the flat fold files and `DATASET.json` in
`data/adni/` using the [data contract](data/adni/README.md). No ADNI
participant records are distributed in Git.

The default design uses ADNI composite memory, FDG-PET SUVR and hippocampal volume
normalized by intracranial volume; nine native study-partner ECog items
(MEMORY1–8 and DIVATT3); and continuous age/education within women. Coordinate
labels, dimensions, masks and seeds come from the supplied protocols. For the
three-item cohort, use `data/adni_3C/` and add `--dataset=3C` to all three
ADNI commands; its paper root is `results/3C/`.

**C2. Fit and analyze EIV-GP.**

```sh
Rscript tools/adni-eivGP.R --action=plan --profile=paper --cores=12 --data-root=data/adni --output=results/adni
Rscript tools/adni-eivGP.R --action=run --profile=paper --cores=12 --data-root=data/adni --output=results/adni
```

Paper runs use five folds, min(cores,4) chains and initially 50,000 iterations per chain
including 5,000 burn-in. Defaults are HMC at O0/O20, ESS at O50, 250 prediction
states per chain and full competitor controls. Explicit numerical overrides are
recorded in the run identity. Optional `--ess-min=100` enables scientific-target
continuation after the final planned stage. O20 calibration is nested in O50
under the supplied C-stratified design.

**C3. Fit and analyze competitors.** This command may run in a second Terminal
while C2 is running, using the same result root:

```sh
Rscript tools/adni-competitors.R --action=run --profile=paper --cores=12 --data-root=data/adni --output=results/adni
```

This covers original ordinal-input UC-GP/LVGP/EzGP and matched complete-case
comparisons. Specify only `--cores` for resource allocation. Fitting and evaluation
use complete chain-sized groups for eivGP and one core per competitor job; lower the budget if memory pressure is
high. Each command summarizes its own results independently. The CPU budget is
per command: two commands requesting 12 cores request 24 workers in total.
See the ADNI guide.

**C4. Assemble and inspect the paper report.** After both fitting commands finish:

```sh
Rscript tools/adni-report.R --action=report --profile=paper --cores=12 --data-root=data/adni --output=results/adni
```

The report includes settings and C definitions, Y/U/X and C/U illustrations,
prediction, complete-case comparisons, scientific diagnostics, response functions
and U imputation. Method reports and paper exports use the same relative paths
as Ocean. `--action=report` rebuilds a method's saved report without fitting;
`adni-report.R` always reads saved results. The report source is
`paper/adni.Rmd`.

| Path | Contents |
| --- | --- |
| `data/adni/` | Privately supplied metadata, five-fold inputs and masks |
| `results/adni/` | Paper fits, predictions, diagnostics and `adni.html` |
| `paper/adni.html` | Versioned data-free plan page |

ADNI results use `adni-flat-v1`, with task folders and `data/`, `logs/`,
`reports/`, `_meta/` directories. Intermediate sampling stages stay inside their
EIV task; final fit, prediction and diagnostic files sit together. `--output`
selects the result root.
All participant inputs, frozen snapshots, fits and empirical reports stay local.
See the [ADNI reproduction and storage guide](paper/README_adni.md).

### D. OAI application

**D1. Supply the private prepared inputs.** Place `DATASET.json` and the 45 flat
fold/protocol files directly in `data/oai/`, following the
[data contract](data/oai/README.md). Participant data are not distributed in Git.
The default application relates baseline MRI medial cartilage thickness to
baseline fixed-location radiographic joint-space width, using continuous age,
BMI and height and eight native osteophyte/sclerosis grades within women.
Dimensions, grade mappings, folds and observation masks come from supplied
protocols. Another prepared cohort can use `data/oai_NAME/` with `--dataset=NAME`.

**D2. Fit and analyze EIV-GP.**

```sh
Rscript tools/oai-eivGP.R --action=plan --profile=paper --cores=12 --data-root=data/oai --output=results/oai
Rscript tools/oai-eivGP.R --action=run --profile=paper --cores=12 --data-root=data/oai --output=results/oai
```

Paper sampling uses five folds, min(cores,4) chains and initially 50,000
iterations including 5,000 burn-in. OAI defaults are HMC at O0/O20 and ESS at O50. Optional
`--ess-min=100` enables scientific-target continuation. Explicit numerical
overrides are recorded in the fitting identity.

**D3. Fit and analyze competitors.** This command may run in a second Terminal
while D2 is running, using the same result root:

```sh
Rscript tools/oai-competitors.R --action=run --profile=paper --cores=12 --data-root=data/oai --output=results/oai
```

This fits UC-GP, LVGP and EzGP with ordinal inputs and matched complete-case
inputs. Use only `--cores` for resources: EIV jobs use complete chain-sized groups
and competitor jobs use one core each. The two commands have separate CPU
budgets: 12 + 12 requests 24 workers in total, so choose their sum to fit the
machine and available memory. Each prints job starts, completions and failures;
`results/oai/PROGRESS_eivGP.txt` and `PROGRESS_competitors.txt` retain the latest
group status. OAI application functions are self-contained in `tools/oai-tools.R`.
Each method report summarizes its own results without requiring the other group
to be complete.

**D4. Assemble and inspect the paper report.** After both fitting commands finish:

```sh
Rscript tools/oai-report.R --action=report --profile=paper --cores=12 --data-root=data/oai --output=results/oai
```

The report covers settings, data/C definitions and illustrations, numerical
diagnostics, OI/CC prediction, response functions,
training-U imputation and reproducibility. `oai-report.R` reads saved fits only;
`--reuse-only=true` also makes either method script saved-only. All methods read
the fixed observation masks supplied with the prepared inputs.

#### Where OAI results are saved

| Path | Contents |
| --- | --- |
| `data/oai/` | Privately supplied metadata, five-fold inputs and masks |
| `results/oai/` | Paper fits, predictions, diagnostics and `oai.html` |
| `paper/oai.html` | Versioned data-free plan page |

OAI results use flat task folders and `data/`, `logs/`, `reports/` and `_meta/`
directories. Participant inputs, fits and empirical reports
stay local. See the [OAI reproduction and storage guide](paper/README_oai.md).

## 5. Check and interpret results

For all three applications, `--output` overrides the result root and `--data-root`
overrides prepared inputs. Use the same dataset, profile and path flags for
fitting and reporting. Completion, optimizer status, scientific-target ESS,
R-hat, MCSE and traces are separate checks. ESS thresholds and a completed
workflow do not establish convergence. Exact numerical reproduction requires
matching seeds, settings, package versions and numerical environment. Saved fits
preserve the posterior histories used for reporting.

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

## Repository layout

- `eivGP/`: reusable model, prediction, diagnostic, storage and reporting functions.
- `tools/`: executable study and application workflows, shared script helpers and developer checks.
- `paper/`: report `.Rmd` sources and `.html` outputs, manuscript sources, figures and tables. Report documents live directly here; `manifests/` retains structured result provenance.
- `data/` and ignored `results/`: inputs and computation outputs; use external stores for large runs.

Rendering a report reads saved results. It does not start or resume an experiment.

## Paths and working directories

Run the documented commands from the checkout root. Study input/output paths,
archive-map paths, maintenance paths, and convergence-manifest source roots are
resolved from the directory where the command starts. Ocean, ADNI and OAI
input/output/config/report paths are resolved from their checkout root; Ocean
fold-preparation outputs follow the same rule. Absolute paths and `~` work in
both cases. New destinations are made absolute before any fitting or rendering
changes the working directory; resolving a path does not create it.

From another directory, call scripts by absolute path and supply absolute input
and output paths. Report rendering locates `paper/` from its script, not the shell's
working directory. Relative roots inside a convergence CSV remain relative to the
launch directory. Keep saved provenance intact; no path repair authorizes moving,
relabeling, deleting, or refitting historical results. Plan checks do not run fits.
