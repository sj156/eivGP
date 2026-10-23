# eivGP

Ocean workflows and the 200/500/1000-profile inputs follow the shared package and paper-analysis structure. See [Ocean instructions](paper/applications/README_ocean.md) for data preparation, fitting and reporting.

The ADNI application uses the same three method-oriented entrypoints with
privately supplied prepared data. See [ADNI instructions](paper/applications/README_adni.md)
for the input contract, pilot/paper settings and saved-result HTML.

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
data/adni/
  inputs/       # local prepared ADNI folds; only .gitkeep is versioned
  prepare_folds.R # preparation from an already matched private cohort
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
install.packages(c("dplyr", "tidyr", "ggplot2", "knitr", "patchwork", "posterior", "jsonlite",
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
- [ADNI application: memory and partially observed imaging](#adni-application-memory-and-partially-observed-imaging)

<a id="numerical-studies-mac-mini-or-linux"></a>

### Synthetic experiments: Studies I and II

Install the locked dependencies and this checkout with `Rscript tools/setup_numerical.R`.
[Study II paper replication](paper/numerical/study2-README.md) separates EIV-GP
and competitors for different machines, with saved-result reporting kept separate.
Each command covers its complete paper workload, including the formerly separate
follow-up comparisons. Optional method selections support partial resumption.

```sh
Rscript paper/numerical/study2-eivGP.R --action=plan
Rscript paper/numerical/study2-competitors.R --action=plan
```

Study I uses `study1-eivGP.R` and `study1-competitors.R` for the same split,
including all calibrated comparisons. See [Study I reproduction](paper/numerical/study1-README.md).
The combined `paper/numerical/study1.R` remains available. See the
[numerical guide](paper/numerical/README.md). All execution requires an explicit
run action. Data, fits, and checkpoints belong in external non-synced result
folders; manuscript sources and figures are exported separately.

### Ocean application: 500-profile cohort

The main Ocean dataset is `ocean500`; both the `.R` runner and `.Rmd` default to
it. Use the same installed checkout as the numerical studies above. For the
separate application installation and dataset guide, see
[Ocean instructions](paper/applications/README_ocean.md).

#### Fit and analyze the Ocean experiment

```sh
Rscript paper/applications/ocean-eivGP.R --plan=true --cores=16

# Fit and analyze each method group, then assemble manuscript comparisons.
Rscript paper/applications/ocean-eivGP.R --cores=16
Rscript paper/applications/ocean-competitors.R --cores=16
Rscript paper/applications/ocean-report.R --cores=16

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

Both commands retain the allocated core budget for analysis and generate or update
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
The versioned example HTML is a 200-profile short-chain workflow check. The revised
manuscript uses the completed 500-profile fits, including all three competitors
and complete-case comparisons. Regenerate its report from those saved fits;
workflow completion alone does not establish convergence. `ocean200` is also available; prepare `ocean1000` folds separately
as described in the Ocean guide.

#### Rebuild the Ocean report without fitting

```sh
# Analyze completed saved fits and scores; fail rather than fit missing models.
Rscript paper/applications/ocean-eivGP.R --reuse-only=true --cores=16
Rscript paper/applications/ocean-competitors.R --reuse-only=true --cores=16
Rscript paper/applications/ocean-report.R --cores=16
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

### ADNI application: memory and partially observed imaging

The default metadata example uses ADNI composite memory (ADNI_MEM), two imaging
coordinates (FDG-PET SUVR and hippocampal volume normalized by intracranial
volume), nine native study-partner ECog items (MEMORY1--8 and DIVATT3), and
continuous age/education within women. Coordinate labels, dimensions, folds,
observation masks and seeds are read from the supplied protocol; other compatible
C sets use the same code. No participant data are included in this repository.

#### Fit and analyze the ADNI experiment

After obtaining the prepared data separately, place its `MC/` directory and
`DATASET.json` under `data/adni/inputs/`. See the [data contract](data/adni/README.md).
Run from the repository root:

```sh
Rscript paper/applications/adni-eivGP.R --plan=true --profile=paper --cores=16

# Fit and analyze each method group, then assemble manuscript comparisons.
Rscript paper/applications/adni-eivGP.R --profile=paper --cores=16
Rscript paper/applications/adni-competitors.R --profile=paper --cores=16
Rscript paper/applications/adni-report.R --profile=paper --cores=16
```

ADNI `paper` uses five folds, four chains and 50,000 total iterations per chain
including 5,000 burn-in. It uses HMC at O0/O20, ESS at O50, 250 prediction states
per chain and full public competitor controls. Explicit numerical overrides are
supported and frozen in the run identity. Use `--profile=paper` consistently
for fitting and reporting. Cores are independently configurable through the shared
public CPU pool; choose `--cores` for the available machine and use the
single-thread numerical-library settings described above. O20 calibration is
nested in O50 under the supplied C-stratified sampling design.

#### Where ADNI results are saved

| Path | Contents |
| --- | --- |
| `data/adni/inputs/` | Privately supplied cohort metadata, five-fold inputs and masks |
| `results/adni_paper/` (or supplied output folder) | Formal-budget fits, predictions, receipts, diagnostics and `adni.html` |
| `paper/applications/adni.html` | Versioned data-free plan page generated from the Rmd |

`--output` overrides the result root; `--data-root` overrides the prepared input
directory. Use the same flags for fitting and reporting. For the three-item
cohort, add `--dataset=3C` to all three entrypoints: inputs go in
`data/adni/3C/inputs/` and paper results in `results/3C_paper/`.

Participant inputs and empirical results stay outside Git; the versioned HTML
contains only the data-free plan. Check task completion and convergence
separately. Detailed output locations are in the
[ADNI guide](paper/applications/README_adni.md#where-adni-results-are-saved).

#### Rebuild the ADNI report without fitting

```sh
# Analyze completed saved fits and scores; fail rather than fit missing models.
Rscript paper/applications/adni-eivGP.R --reuse-only=true --profile=paper --cores=16
Rscript paper/applications/adni-competitors.R --reuse-only=true --profile=paper --cores=16
Rscript paper/applications/adni-report.R --profile=paper --cores=16
```

The report is generated through `paper/applications/adni.Rmd`. It includes the
complete C definition table, Y/U/X and C/U illustrations, response prediction,
matched complete-case comparisons, numerical diagnostics, response functions and
imputation. `adni-report.R` always reads saved results. See the
[ADNI guide](paper/applications/README_adni.md) for input preparation, profiles,
numerical overrides and report settings.

#### Quick workflow check without ADNI data

```sh
Rscript tests/application/check-adni.R --cores=2
```

This creates a small synthetic fixture under ignored `results/`, exercises the
three command-line entrypoints with reduced settings, checks reports and
core-count-only reuse, and writes `CHECK.json`. It does not fit participant data.
Completion is checked separately from convergence.

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


For current numerical commands, see [Study II paper replication](paper/numerical/study2-README.md) and [the Study I guide](paper/numerical/README.md).

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

### Ocean: one CPU allocation, complete paper exports

After installation, run `ocean-eivGP.R`, `ocean-competitors.R`, then
`ocean-report.R` from `paper/applications/`, each with `--cores=16`. The first two
combine fitting and analysis by method group; the third assembles paper tables,
figures and HTML. Use `--plan=true` to preview or `--reuse-only=true` on the method
scripts to prohibit fitting when reproducing completed saved results. See [Ocean reproduction](paper/applications/README_ocean.md).
`reports/paper/artifact_manifest.csv` maps the manuscript labels to generated files;
method failures remain explicit and no historical 200-profile values are relabeled.
