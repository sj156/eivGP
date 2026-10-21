# eivGP

Ocean workflows and the 200/500/1000-profile inputs follow the shared package and paper-analysis structure. See [Ocean instructions](paper/applications/README_ocean.md) for data preparation, fitting and reporting.

Gaussian process regression with ordinal measurements and potentially available
continuous measurements associated with those inputs. The package provides
fitting, prediction, imputation, diagnostics, continuation, and reporting.

## Code structure

```
eivGP/R/
  helpers.R     # data preparation, kernels, prediction, storage, workflow helpers
  mcmc.R        # posterior sampling, fitting, continuation, sampler diagnostics
  hmc.R         # joint HMC coordinates, warmup adaptation and transitions
eivGP/src/      # analytic GP/ordinal gradients and leapfrog integration
  plotting.R    # summaries, comparison tables, figures, report generation
data/ocean_silicate/
  n200/         # frozen 200-profile data and five-fold inputs
  n500/         # updated 500-profile data and five-fold inputs
  n1000/        # raw selected 1000-profile cohort
  prepare_folds.R # standalone fold preparation entrypoint
paper/
  numerical/    # study1 and study2: one .R, .Rmd, and .html each
  applications/ # ocean and adni: one .R, .Rmd, and .html each
  config/       # Ocean dataset/compute parameter templates
  manifests/    # selected numerical runs and frozen Ocean input identities
tools/
  build_package.R
  render_reports.R
```

`eivGP/R/` is the authoritative implementation. The R modules and native HMC implementation are installed together as one package. Metadata and static-analysis declarations remain separate. Analysis scripts call package APIs and select
settings, inputs, and output locations. Competitor adapters call `kergp`, `LVGP`,
and `EzGP`. Reporting functions also live in the package.

## Install

Use R >= 4.1. Install dependencies once, then install this checkout:

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "knitr", "patchwork", "posterior",
                   "TruncatedNormal", "digest", "rmarkdown", "testthat", "Rcpp", "RcppArmadillo"))
# Needed when fitting the competitor methods:
install.packages(c("kergp", "LVGP", "EzGP"))
```

```sh
Rscript tools/build_package.R
R CMD INSTALL eivGP
```

The build helper refreshes the module identity manifest and checks parsing. It
never generates or replaces the package source. Reports require Pandoc, such as
the version included with RStudio.

## Choose ESS, HMC or the HMC–ESS hybrid

The two engines support an explicit main continuous-update option:

```r
fit_ess <- fit_eivgp(X, y, C, U_obs = U_obs, engine = "univariate",
                    sampler = "ESS")
fit_hmc <- fit_eivgp(X, y, C, U_obs = U_obs, engine = "univariate",
                    sampler = "HMC")
fit_hybrid <- fit_eivgp(X, y, C, U_obs = U_obs, engine = "univariate",
                       sampler = "HYBRID",
                       sampler_control = list(hybrid_every = 20L))
# For Study II-style data, use engine="multivariate", latent_dim=2L.
continued <- continue_eivgp(fit_hmc, n_iter = 1000L)
fit_hmc$diagnostics$hmc  # chain-specific acceptance and energy-error diagnostics
```

ESS remains the default. In the representative matched-CPU comparison, the
every-20 hybrid improved efficiency on many targets but did not remove
persistent scientific-target separation in the difficult calibrated Study I
case. A favorable median effective sample size does not validate posterior
basin weights. Joint HMC moves all missing inputs, observed-input
kernel coefficients and the signal fraction together, retaining the same
full-GP posterior and existing measurement/dictionary updates. Thus auxiliary
ESS updates remain in the HMC-within-Gibbs algorithm. HMC supports conditional
dictionary updates and SE or Matérn 3/2 and 5/2 kernels; select ESS for other
kernel smoothness values or marginal dictionary updates. Unsupported explicit
HMC or HYBRID requests raise an error.

HYBRID performs the same HMC sweep at every iteration and adds a complete
ESS sweep after every 20 iterations by default, retaining one draw after the
combined transition. `hybrid_every` must be a positive integer and applies only
to HYBRID. This composes posterior-preserving kernels without changing the
full GP, priors, or calibration constraints. The extra sweep includes the
measurement and dictionary updates. Its count is saved in each checkpoint
state as `hybrid_extra_sweeps`.

HMC controls can be supplied through
`sampler_control=list(hmc=list(step=.05, leapfrog=c(8L,24L), target_accept=.8))`.
Step size and a diagonal mass matrix adapt during warmup only. Continuation
preserves the sampler, frozen metric, step size and RNG state; it cannot switch
methods, hybrid frequency or tuning. The default retains every post-warmup draw. Inspect
scientific-target R-hat, bulk/tail ESS and MCSE rather than assuming that a
higher median ESS per second certifies all posterior quantities.

## Reproduce and read the paper analyses

See [paper/README.md](paper/README.md) for the analysis workflows. Use `.R`
entrypoints for computation and `.Rmd` notebooks for reader-facing HTML reports.

- [Synthetic experiments: Studies I and II](#synthetic-experiments-studies-i-and-ii)
- [Ocean application: 500-profile cohort](#ocean-application-500-profile-cohort)

<a id="numerical-studies-mac-mini-or-linux"></a>

### Synthetic experiments: Studies I and II

#### Setup: Mac mini or Linux

Run commands from the repository root. For numerical reproduction, use the locked
setup below instead of the general package-install instructions above. Install
Git and R first (the lockfile records R 4.6.0); source packages may require
C/C++/Fortran compilers and CMake. See the [platform prerequisites](paper/numerical/README.md#macos-or-linux-setup).
No private data or archived fits are required.

```sh
# First checkout:
git clone https://github.com/sj156/eivGP.git
cd eivGP
# On an existing clean checkout, use this instead of cloning:
# git pull --ff-only

Rscript tools/setup_numerical.R
Rscript tools/setup_numerical.R --check
Rscript paper/numerical/check-portability.R
Rscript paper/numerical/check-numerical.R
```

The setup restores `renv.lock`, installs the package into `.r-library/`, and
checks its source identity. Freeze the checkout and dependencies during a run;
do not pull changes midway through an experiment you intend to resume.

#### Plan and run the complete studies

```sh
# Inspect the workload without fitting.
Rscript paper/numerical/study1.R --action=plan --profile=paper --cores=4
Rscript paper/numerical/study2.R --action=plan --profile=paper --cores=4

# Generate missing data, fit models, and aggregate results.
Rscript paper/numerical/study1.R --action=run --profile=paper --cores=4
Rscript paper/numerical/study2.R --action=run --profile=paper --cores=4
```

Study I runs the univariate grid. Study II runs the original multivariate grid
and the targeted CUT/complete-case comparisons. `--component=core` restricts
execution to an original grid; `--component=mechanism` selects the targeted
Study II comparisons. Four cores is a starting choice, not a runtime estimate.
Quick profiles test software only and must use a separate output root.

#### Run only missing or new experiments

To resume an interrupted run, reissue its original command with the **same
output root, profile and scientific settings**. Compatible saved work is reused;
missing work is computed. Mechanism jobs with saved scores are skipped. Original
grids reuse their available caches/checkpoints; interrupted work without a
checkpoint may have to be recomputed.

```sh
# Paired cost/mixing pilot: two replications, two calibration levels.
# Four FULL references plus 16 new CUT/complete-case jobs.
Rscript paper/numerical/study2.R --component=mechanism --action=run \
  --profile=pilot --cores=4 --output="$HOME/eivGP-results/study2-pilot"

# Inspect mechanism/report.html and stage_timing.csv before scaling up.
# Paper runs need a separate root because their sampling budgets differ.
# After review, a fresh complete Study II run uses the default two components.

# Only CUT and complete-case methods, if compatible baseline jobs already exist.
Rscript paper/numerical/study2.R --component=mechanism --action=run \
  --profile=paper --reps=1:50 --cores=4 \
  --methods=CUT,UC-GP-CC,LVGP-CC,EzGP-CC

# Retry one missing/failed job; this does not refit an already scored job.
Rscript paper/numerical/study2.R --component=mechanism --action=run \
  --profile=paper --reps=7 --calib=20 --methods=LVGP-CC --cores=4
```

The last two commands alone do not produce a complete comparison without the
other methods. Mechanism summaries include all jobs registered in that output
root, not just the latest selection. Historical/pre-refactor outputs are not
made compatible by moving or renaming them. See [optional fit imports and resume
rules](paper/numerical/README.md#resume-transfer-and-optional-imports).

#### Where results are saved

With no environment/path overrides, numerical outputs are saved under `~/eivGP-results/`, outside the checkout:

| Path | Contents |
| --- | --- |
| `~/eivGP-results/data/study1/`, `~/eivGP-results/data/study2/` | Fixed-seed synthetic datasets shared across methods |
| `~/eivGP-results/study1/core/<run-id>-<fingerprint>/` | Univariate run, caches, diagnostics and aggregate outputs |
| `~/eivGP-results/study2/core/<run-id>-<fingerprint>/` | Original multivariate run and aggregate outputs |
| `~/eivGP-results/study2/mechanism/` | Targeted comparison jobs, diagnostics and score tables |
| `.r-library/` | Local R installation and dependency environment |

Core runs put status files in `config/`, aggregate CSVs in `combined/`, and
combined competitor outputs in `reporting/comparison/`. Mechanism outputs include
`point_scores.csv`, `replicate_scores.csv`, `summary.csv`, and
`paired_differences.csv`, paired summaries and unseen-minus-seen interactions.
`report.html` includes completion, diagnostics, paired scores and stage wall times. Per-job directories retain diagnostic files and
predictions; shared directories retain measurement diagnostics.

Use `--output=/absolute/local/path/study2` to change the study's parent result
root; `core/` and `mechanism/` are appended automatically. Use
`--data-root=/absolute/local/path/data/study2` for that study's input directory.
Repeat these overrides on every run/report command. `EIVGP_RESULTS_ROOT` and
`EIVGP_DATA_ROOT` also override defaults; explicit flags take precedence. Keep
outputs outside Dropbox/Overleaf. Do not run two processes against the same
output root. The two studies have separate default roots.

#### Checkpoints and runtime forecasts

Sampling saves explicitly every 300 seconds or 1,000 iterations, whichever is
reached first, plus at completion. Time checks occur at iteration boundaries.
Timing and chain forecasts refresh at least every 60 seconds of sampling, also
at iteration boundaries; batch forecasts update as tasks start and finish.
Use `--checkpoint-seconds=120 --checkpoint-iterations=500` to change the save
interval. These operational settings do not change the statistical model.
See [restart and timing details](paper/numerical/README.md#checkpoints-and-runtime-forecasts)
for file locations, forecast scope, and recovery after a power outage.

#### Check completion and convergence

**Completion means outputs were computed; it does not mean convergence passed.**
Inspect both the status and diagnostic files:

```sh
# Locate core run status and gate summaries (use the printed run directory).
find "$HOME/eivGP-results/study1/core" "$HOME/eivGP-results/study2/core" -type f \
  \( -name cell_status.csv -o -name diagnostic_gates.csv -o -name statuses.csv \)

# Count mechanism jobs by completion and diagnostic status; show problems.
Rscript -e 's <- read.csv(path.expand("~/eivGP-results/study2/mechanism/status.csv")); print(with(s, table(status, diagnostics_pass, useNA="ifany"))); print(s[s$status != "complete" | is.na(s$diagnostics_pass) | !s$diagnostics_pass, ])'
```

Run these after the relevant workflow has created its outputs. For core runs,
check `config/cell_status.csv` for `failed`, `gate_failed`, or
`completed_with_gate_warnings`, and inspect `config/diagnostic_gates.csv` for
failed gates. Check competitor `statuses.csv` files too. An early failure may
leave later cells absent from the status table; compare with the planned grid.

For mechanism runs, `status.csv` lists registered jobs as `complete`, `failed`,
or `pending`; `diagnostics_pass` reports screening separately. A partial selection
can be entirely complete while the full paper design remains unfinished.
Inspect `full_diagnostics.csv`, `conditional-*.diagnostics.csv`, shared
`measurement_diagnostics.csv`, and `integration_stability.csv` for flagged jobs.
MCMC screens require R-hat <= 1.01 and bulk/tail ESS >= 400. Optimizer failures
are also flagged. Scores from flagged fits remain visible and must not be
silently treated as validated results.

A completed but diagnostically flagged job is not automatically refitted.
Investigate it first; changed sampling budgets/code require a **new output
root**, and all conditional-fit histories are retained for investigation. Passing screens
still require target-specific mixing and nested-integration precision review.
An abrupt kill can leave a `RUN.lock` or job `.lock`; inspect its PID and confirm
the process has stopped before removing that exact lock. Never remove a live lock.

#### Rebuild reports without fitting

```sh
Rscript paper/numerical/study1.R --action=report --profile=paper
Rscript paper/numerical/study2.R --action=report --profile=paper
# If only targeted Study II work has been run:
Rscript paper/numerical/study2.R --component=mechanism --action=report --profile=paper
```

Report commands require existing results and the same settings/path overrides;
they do not fill missing fits. For custom mechanism sampler budgets, repeat those
flags when reporting. The [numerical guide](paper/numerical/README.md) describes
statistical targets, provenance, and compact output retention. Publication-scale
storage has not been measured; all conditional fits retain chain histories.
The [manuscript section](paper/numerical/manuscript/section.tex) retains TBDs for
pending scientific validation and final displays.

### Ocean application: 500-profile cohort

The main Ocean dataset is `ocean500`; both the `.R` runner and `.Rmd` default to
it. Use the same installed checkout as the numerical studies above. For the
separate application installation and dataset guide, see
[Ocean instructions](paper/applications/README_ocean.md).

#### Fit and analyze the Ocean experiment

```sh
Rscript paper/applications/ocean.R --action=plan --profile=paper --dataset=ocean500 --cores=16

# Fit, save checkpoints, analyze with the shared core budget, and generate HTML.
Rscript paper/applications/ocean.R --action=fit --profile=paper --dataset=ocean500 --cores=16

# Explicitly extend flagged fits, then update analysis and HTML.
Rscript paper/applications/ocean.R --action=continue --profile=paper --dataset=ocean500 --continuation-iter=20000 --cores=16
```

Sampling settings are fixed by profile. `paper` uses four chains, 50,000 initial
iterations per chain including 5,000 burn-in; `pilot` uses four chains and 80/30;
`quick` uses one chain and 80/30. `--chains`, `--iter` and `--burn` are not CLI
options. Use `pilot` or `quick` for workflow checks.

`--cores` is one shared budget across MCMC, prediction, diagnostics, scoring,
and figure export. Four chains per paper fit are unchanged. Waiting coordinators
release their slots; active chain and analysis workers acquire slots from the
same pool. After MCMC, prediction and post-processing can borrow the full budget,
including when only one fit remains. The number of simultaneously resident fits
is bounded separately to limit memory use. HTML assembly remains serial.
Windows currently executes serially. Allocation logs under `logs/<run_id>/` in new flat Ocean roots (legacy `jobs/<run_id>/` remains readable)
count outer job coordinators, not instantaneous CPU utilization.

Set numerical-library threads to one **before starting R**, so each active worker
uses at most one numerical-library thread:

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
```

Numerical study entrypoints and the ADNI driver use the same pool. Random
prediction tasks use fixed blocks of eight selected posterior states and task
seeds independent of core count (prediction stream schema 1). Results therefore
reproduce across core counts, including joint predictions and Gaussian mixture
components. MCMC streams and chain counts are unchanged. Prediction Monte Carlo
realizations can differ from the older serial stream; source/cache identities
must match. Do not relabel old results or remove identity files to resume an
incompatible run. Existing saved fits and outputs are not migrated automatically.


Continuation keeps the original chains and adds 20,000 iterations per chain by
default in paper mode. `--continuation-iter=10000` requests another extension
length. Each fit supports one explicit extension; passing fits are skipped.
`run` aliases `fit`. Independent `analyze` and `report` actions remain available.

#### Minimal command-line workflow check

Run these commands after installation. They use ocean200, two folds/O20, three
CC comparators per fold and an 8-core fitting budget. A temporary result folder
keeps this check separate from the formal ocean500 analysis.

```sh
OCEAN_CHECK_OUTPUT="$(mktemp -d "${TMPDIR:-/tmp}/eivGP-ocean-check.XXXXXX")"
Rscript paper/applications/ocean.R --action=fit --profile=pilot --dataset=ocean200 --folds=1,2 --cells=O20 --modules=eiv,complete_case_baselines,continuation,prediction --cores=8 --output="$OCEAN_CHECK_OUTPUT"
Rscript paper/applications/ocean.R --action=continue --profile=pilot --dataset=ocean200 --folds=1,2 --cells=O20 --modules=eiv,complete_case_baselines,continuation,prediction --continuation-iter=40 --cores=8 --output="$OCEAN_CHECK_OUTPUT"
```

Both commands finish with single-worker analysis and generate or update
`$OCEAN_CHECK_OUTPUT/ocean.html`. Only the second requests continuation. The CLI
finds a custom library through `EIVGP_R_LIBRARY`, otherwise through `.r-library/`
or standard library paths. Pandoc must be available before fitting starts.

#### Where Ocean results are saved

| Path | Contents |
| --- | --- |
| `data/ocean_silicate/n500/` | Frozen 500-profile cohort, five-fold inputs and masks |
| `results/ocean500/` (or printed run folder) | Formal-budget fits, predictions, receipts, diagnostics and `ocean.html` |
| `results/ocean500_pilot/` | Separate short-chain workflow results |
| `paper/applications/ocean.html` | Reader-facing example report generated from the Rmd |

`--output=/absolute/path/ocean500-paper` overrides the result root;
`--data-root` overrides the prepared input directory. Both overrides must also
be supplied to report generation. Check `jobs/ocean500_paper_fit/task_status.csv`
for completion and the saved convergence/optimizer diagnostics separately.
The versioned HTML is the verified 200-profile short-chain example. The 500-profile
formal run remains to be completed; a finished short-chain run is not evidence
of convergence. `ocean200` is also available; prepare `ocean1000` folds separately
as described in the Ocean guide.

#### Rebuild the Ocean report without fitting

```sh
# Analyze saved chains and regenerate HTML; no model fitting.
Rscript paper/applications/ocean.R --action=analyze --profile=paper --dataset=ocean500
# Rebuild HTML from saved predictions and summaries.
Rscript paper/applications/ocean.R --action=report --profile=paper --dataset=ocean500
```

The CLI already generates HTML. Direct Rmd rendering remains available using the same dataset, profile and result root:

```r
rmarkdown::render("paper/applications/ocean.Rmd",
  params=list(dataset="ocean500",profile="paper",action="report",
    output_dir=""))
```

The `.R` script calls separate fitting and analysis stages in sequence; the `.Rmd` reads saved results and produces
HTML. The notebook rejects `action="run"`. See the
[Ocean guide](paper/applications/README_ocean.md) for installation, optional task
selectors, output imports and dataset preparation.

## Develop and test

Edit the package modules and native sampler sources, then refresh the manifest
and reinstall. API
exports and help pages are maintained in `NAMESPACE` and `man/`.

```r
library(eivGP)
testthat::test_dir("eivGP/tests/testthat", env=new.env(parent=asNamespace("eivGP")),
                   stop_on_failure=TRUE)
```

Version 0.4.0 consolidates the implementation and analysis/reporting workflows.
It retains the 0.3.1 sampler and scientific defaults. Package tests cover both
engines, prediction, diagnostics, continuation, and lossless compact storage.


Both numerical entry points now include the two published-competitor protocols: original `(X,C)` fits on all records and complete-case `(X,U,C)` fits with `U` continuous. `study1.R --component=complete-case --profile=pilot --action=run` starts a 12-job paired pilot; `study2.R --component=mechanism --profile=pilot --action=run` starts the CUT/CC pilot. Full Study I additions (600 CC fits) are separate from the existing 400 new Study II jobs. See [the numerical workflow](paper/numerical/README.md#two-versions-of-the-published-competitors-in-both-studies) for data access, diagnostics, reporting and the two-script paper reproduction commands. Pilots never launch the full grids automatically.

CPU requests are upper bounds: the effective shared budget is the smaller of
`--cores` and the detected logical CPU count. For example, `--cores=56` uses
16 workers on a 16-CPU host, while `--cores=8` uses eight. Failed CPU detection
falls back to one worker; Windows remains serial. Four-chain fits retain all
four chains even when fewer workers are available. Keep BLAS/OpenMP threads
limited to one per worker. Detection describes host CPU capacity, not current
idle CPUs or a guaranteed scheduler/container allocation.

New Ocean roots use storage schema `ocean-flat-v1`: EIV, continuation and
competitor task folders are siblings under the dataset root, with `data/`,
`logs/`, `reports/` and `_meta/` for shared material. Existing unmarked roots
retain their layout. See [Ocean storage and migration](paper/applications/README_ocean.md)
for the directory specification and the dry-run migration command. No automatic
migration runs during fitting, prediction, report generation or resume.
