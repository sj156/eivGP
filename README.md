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

`eivGP/R/` is the authoritative implementation. The three implementation modules are installed together as one package, following the upstream structure. Metadata and static-analysis declarations remain separate. Analysis scripts call package APIs and select
settings, inputs, and output locations. Competitor adapters call `kergp`, `LVGP`,
and `EzGP`. Reporting functions also live in the package.

## Install

Use R >= 4.1. Install dependencies once, then install this checkout:

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "knitr", "patchwork", "posterior",
                   "TruncatedNormal", "digest", "rmarkdown", "testthat"))
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
# Targeted comparisons only, including FULL and categorical-only baselines.
# Start with one replication to assess cost and diagnostics.
Rscript paper/numerical/study2.R --component=mechanism --action=run \
  --profile=paper --reps=1 --cores=4

# Expand or resume the same protocol; completed jobs are reused.
Rscript paper/numerical/study2.R --component=mechanism --action=run \
  --profile=paper --reps=1:50 --cores=4

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

With no environment/path overrides, paths are relative to this checkout:

| Path | Contents |
| --- | --- |
| `results/data/study1/`, `results/data/study2/` | Fixed-seed synthetic datasets shared across methods |
| `results/study1/core/<run-id>-<fingerprint>/` | Univariate run, caches, diagnostics and aggregate outputs |
| `results/study2/core/<run-id>-<fingerprint>/` | Original multivariate run and aggregate outputs |
| `results/study2/mechanism/` | Targeted comparison jobs, diagnostics and score tables |
| `.r-library/` | Local R installation and dependency environment |

Core runs put status files in `config/`, aggregate CSVs in `combined/`, and
combined competitor outputs in `reporting/comparison/`. Mechanism outputs include
`point_scores.csv`, `replicate_scores.csv`, `summary.csv`, and
`paired_differences.csv`. Per-job directories retain diagnostic files and
predictions; shared directories retain measurement diagnostics.

Use `--output=/absolute/local/path/study2` to change the study's parent result
root; `core/` and `mechanism/` are appended automatically. Use
`--data-root=/absolute/local/path/data/study2` for that study's input directory.
Repeat these overrides on every run/report command. `EIVGP_RESULTS_ROOT` and
`EIVGP_DATA_ROOT` also override defaults; explicit flags take precedence. Keep
outputs outside Dropbox/Overleaf. Do not run two processes against the same
output root. The two studies have separate default roots.

#### Check completion and convergence

**Completion means outputs were computed; it does not mean convergence passed.**
Inspect both the status and diagnostic files:

```sh
# Locate core run status and gate summaries (use the printed run directory).
find results/study1/core results/study2/core -type f \
  \( -name cell_status.csv -o -name diagnostic_gates.csv -o -name statuses.csv \)

# Count mechanism jobs by completion and diagnostic status; show problems.
Rscript -e 's <- read.csv("results/study2/mechanism/status.csv"); print(with(s, table(status, diagnostics_pass, useNA="ifany"))); print(s[s$status != "complete" | is.na(s$diagnostics_pass) | !s$diagnostics_pass, ])'
```

Run these after the relevant workflow has created its outputs. For core runs,
check `config/cell_status.csv` for `failed`, `gate_failed`, or
`completed_with_gate_warnings`, and inspect `config/diagnostic_gates.csv` for
failed gates. Check competitor `statuses.csv` files too. An early failure may
leave later cells absent from the status table; compare with the planned grid.

For mechanism runs, `status.csv` lists registered jobs as `complete`, `failed`,
or `pending`; `diagnostics_pass` reports screening separately. A partial selection
can be entirely complete while the full paper design remains unfinished.
Inspect `full_diagnostics.csv`, `conditional-*-diagnostics.csv`, shared
`measurement_diagnostics.csv`, and `integration_stability.csv` for flagged jobs.
MCMC screens require R-hat <= 1.01 and bulk/tail ESS >= 400. Optimizer failures
are also flagged. Scores from flagged fits remain visible and must not be
silently treated as validated results.

A completed but diagnostically flagged job is not automatically refitted.
Investigate it first; changed sampling budgets/code require a **new output
root**, and failed-fit histories are retained for investigation. Passing screens
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
storage has not been measured; failed fits can still retain substantial history.
The [manuscript section](paper/numerical/manuscript/section.tex) retains TBDs for
pending scientific validation and final displays.

### Ocean application: 500-profile cohort

The main Ocean dataset is `ocean500`; both the `.R` runner and `.Rmd` default to
it. Use the same installed checkout as the numerical studies above. For the
separate application installation and dataset guide, see
[Ocean instructions](paper/applications/README_ocean.md).

#### Fit and analyze the Ocean experiment

```sh
Rscript paper/applications/ocean.R --action=plan --profile=paper --dataset=ocean500

# Fit and save initial checkpoints. No automatic continuation or EIV postprocessing.
Rscript paper/applications/ocean.R --action=fit --profile=paper --dataset=ocean500 --cores=4

# Read saved chains and analyze using one worker.
Rscript paper/applications/ocean.R --action=analyze --profile=paper --dataset=ocean500

# Explicitly extend initially flagged chains if needed, then analyze again.
Rscript paper/applications/ocean.R --action=continue --profile=paper --dataset=ocean500 --continuation-iter=20000 --cores=4
Rscript paper/applications/ocean.R --action=analyze --profile=paper --dataset=ocean500
```

The `paper` profile uses four chains, 20,000 initial iterations and 5,000 burn-in;
`continue` requests one screen-selected extension, defaulting to 20,000 additional
iterations per chain. Set `--continuation-iter=10000` on the `continue` command
to choose another length after initial fitting and analysis. Saved lengths are
reused by subsequent analysis and reporting; an already-started extension cannot
be changed or extended a second time. Use `profile=pilot`
for short chains (80 initial, 30 burn-in, +80 extension). All stages are separate;
`run` aliases `fit`. Analysis always uses one worker and has no core-allocation
option. Repeating a stage reuses compatible results; changed inputs, source or
initial-fit numerical settings require a fresh output root. Each output root has one owner.

For a 10,000-iteration extension per chain, use this command at the continuation
step, then run `analyze` again:

```sh
Rscript paper/applications/ocean.R --action=continue --profile=paper --dataset=ocean500 --continuation-iter=10000 --cores=4
```

#### Minimal command-line workflow check

After installing this checkout, copy these commands into a Terminal at the
repository root. This uses the 200-profile inputs, one fold and one calibration
condition with short chains, plus the three CC comparators. It checks the staged
workflow; the main experiment above continues to use 500 profiles. Use a fresh
`--output` directory after changing source or settings.

```sh
Rscript paper/applications/ocean.R --action=fit --profile=pilot --dataset=ocean200 --folds=1 --cells=O20 --modules=eiv,complete_case_baselines,continuation,prediction --cores=2 --output=results/ocean/ocean200/terminal-check
Rscript paper/applications/ocean.R --action=analyze --profile=pilot --dataset=ocean200 --folds=1 --cells=O20 --modules=eiv,complete_case_baselines,continuation,prediction --output=results/ocean/ocean200/terminal-check
Rscript paper/applications/ocean.R --action=continue --profile=pilot --dataset=ocean200 --folds=1 --cells=O20 --modules=eiv,complete_case_baselines,continuation,prediction --continuation-iter=40 --cores=2 --output=results/ocean/ocean200/terminal-check
Rscript paper/applications/ocean.R --action=analyze --profile=pilot --dataset=ocean200 --folds=1 --cells=O20 --modules=eiv,complete_case_baselines,continuation,prediction --output=results/ocean/ocean200/terminal-check
Rscript paper/applications/ocean.R --action=report --profile=pilot --dataset=ocean200 --folds=1 --cells=O20 --modules=eiv,complete_case_baselines,continuation,prediction --output=results/ocean/ocean200/terminal-check
Rscript -e ".libPaths(c(Sys.getenv('EIVGP_R_LIBRARY','.r-library'),.libPaths())); rmarkdown::render('paper/applications/ocean.Rmd', output_file='ocean_check.html', output_dir='results/ocean/ocean200/terminal-check', params=list(dataset='ocean200', profile='pilot', action='report', folds=1, cells='O20', modules=c('eiv','complete_case_baselines','continuation','prediction'), output_dir='results/ocean/ocean200/terminal-check', appendix=FALSE))"
```

The first analysis does not extend chains. The explicit `continue` command
adds 40 iterations per chain to a flagged fit, and the next analysis reads the saved extension. The
last command renders `results/ocean/ocean200/terminal-check/ocean_check.html`.
The final command locates the installed package through `EIVGP_R_LIBRARY` when
set, otherwise through the checkout-local `.r-library/`.

#### Where Ocean results are saved

| Path | Contents |
| --- | --- |
| `data/ocean_silicate/n500/` | Frozen 500-profile cohort, five-fold inputs and masks |
| `results/ocean/ocean500/paper/` | Formal-budget fits, predictions, receipts and diagnostics |
| `results/ocean/ocean500/pilot/` | Separate short-chain workflow results |
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
# Rebuild saved-result summaries and figures; no model fitting.
Rscript paper/applications/ocean.R --action=report --profile=paper --dataset=ocean500 --cores=4
```

Render the reader-facing HTML in R using the same dataset, profile and result root:

```r
rmarkdown::render("paper/applications/ocean.Rmd",
  params=list(dataset="ocean500",profile="paper",action="report",
    output_dir="results/ocean/ocean500/paper"))
```

The `.R` script separates fitting, analysis and continuation; the `.Rmd` reads saved results and produces
HTML. The notebook rejects `action="run"`. See the
[Ocean guide](paper/applications/README_ocean.md) for installation, optional task
selectors, output imports and dataset preparation.

## Develop and test

Edit the three package modules, then refresh the manifest and reinstall. API
exports and help pages are maintained in `NAMESPACE` and `man/`.

```r
library(eivGP)
testthat::test_dir("eivGP/tests/testthat", env=new.env(parent=asNamespace("eivGP")),
                   stop_on_failure=TRUE)
```

Version 0.4.0 consolidates the implementation and analysis/reporting workflows.
It retains the 0.3.1 sampler and scientific defaults. Package tests cover both
engines, prediction, diagnostics, continuation, and lossless compact storage.
