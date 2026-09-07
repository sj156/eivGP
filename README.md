# eivGP 0.3.1

`eivGP` fits Gaussian process regressions with numeric predictors, ordinal
proxies for latent continuous inputs, and optional calibration measurements.
The reusable package and the manuscript's numerical experiments are separate:
`eivGP/` is the single installable package; `codes/` is its canonical source;
`experiments/` and the repository simulation helpers define the studies.

The current source default is `r ~ Beta(13,3)`, with probability 0.9021 between
0.65 and 0.95 and mean 0.8125. Earlier 0.3.1 fits used Beta(32,8); version alone
does not identify the prior. Inspect saved `priors$signal_shape`. Continuation
preserves a fit's original prior; adopting the new prior requires a fresh fit.

The original version 0.3.1 release improved transitions without changing the 0.3.0 posterior:
deterministic-threshold cutoffs use exact coordinate Gibbs updates, while
ordinal-probit joint moves reuse exact block-GP calculations and evaluate
only affected ordinal factors. Existing cutoff/input transports are retained.
Set `sampler_control = list(threshold_update = "ess")` to compare with the
previous cutoff transition; the default is `"gibbs"`. Neither option thins.
Continuation requires a matching 0.3.1 checkpoint. Existing 0.3.0 draws target
the same posterior, but should be continued using 0.3.0, or refitted to use
the new transition schedule. Do not relabel old checkpoints.

Version 0.3.0 restored the package name `eivGP` and replaced the previous
posterior prior and sampler. Fits and checkpoints from `eivGP` 0.1.x or
`eivmixgp` 0.1–0.2.x must be refitted. Renaming or relabeling an old fit does
not migrate it to the new posterior.

## Repository layout

The working tree keeps the current package and paper-experiment implementation:

- `codes/`: canonical model modules, experiment helpers, and regression tests.
- `eivGP/`: the generated, installable package (keep this alongside its sources).
- `litr/`: the active literate package build and build-support tests.
- `experiments/`: dependency setup, frozen-data generation, competitor/MCMC
  launchers, and separate reporting.
- `applications/` and `codes/real-data/`: application work in progress.
- `reproduction/`: local frozen data, reusable caches, and saved results;
  Git-ignored and not distributed automatically with the repository.

Obsolete package/book sources, the old `artifacts/` tree, archived superseded
data, and selected older result runs have been removed. Current frozen data,
competitor caches, and the latest saved development run per study are retained.
Retained results must still be checked against the current design; being the
latest saved run does not make them final publication results. Previously
committed versions remain in Git history; cleanup does not rewrite history.

## Install

The current five-setting numerical design (100 training / 100 test observations,
50 publication versus 3 development replications) is documented in
[NUMERICAL_DESIGN.md](NUMERICAL_DESIGN.md). Publication MCMC budgets remain
to be finalized. Existing frozen datasets have not been replaced.

For package-only use, install its dependencies outside the package, then install
from a repository checkout:

```sh
Rscript --vanilla -e 'install.packages(c("posterior", "TruncatedNormal"), repos="https://cloud.r-project.org")'
R CMD INSTALL eivGP
```

The required non-base R dependencies are `posterior` and `TruncatedNormal`.
To regenerate the package from canonical code, install `litr`, `rmarkdown`,
`usethis`, `devtools`, and `roxygen2`, then run:

```sh
Rscript --vanilla litr/render-package.R
R CMD build eivGP
R CMD check --no-manual eivGP_0.3.1.tar.gz
```

## Reproduce the numerical experiments: a new-reader workflow

You can currently run the revised experiments in **development mode**. This is
not yet exact reproduction of final paper tables: the revised publication
results have not been generated, the publication MCMC budget is specified below,
and shared frozen-data storage between modes is not yet implemented.
Development uses the same five settings and sample sizes, but only three
replicated datasets per setting instead of 50. Real-data reproduction
instructions will be added separately when that workflow is ready.

### 1. Download the repository

Install Git, R >= 4.1, and Pandoc (available with RStudio). In Terminal:

```sh
git clone https://github.com/sj156/eivGP.git
cd eivGP
git rev-parse HEAD
```

Keep the commit identifier with your results. Run all subsequent commands from
this repository root. An existing checkout does not need to be cloned again.

### 2. Install dependencies and eivGP

Install the experiment dependencies outside the package:

```sh
Rscript --vanilla experiments/install_eivgp_dependencies.R
R CMD INSTALL eivGP
Rscript --vanilla -e 'library(eivGP); print(packageVersion("eivGP")); stopifnot(rmarkdown::pandoc_available())'
```

The package version for this workflow is 0.3.1. If dependency installation
fails, resolve the reported error before continuing. Missing optional
competitors (`kergp`, `LVGP`, or `EzGP`) are recorded and skipped, so their
comparisons will be incomplete. If Pandoc is not detected from Terminal,
make it available there before rendering the experiment documents.

The setup script installs only missing dependencies, including `kergp` (UC-GP),
`LVGP`, and `EzGP`, and prints their versions and library locations. It belongs
to `experiments/`, not the installed package. Package loading, fitting, and
reporting never install software automatically. Run setup on each machine,
using the same R installation as the experiments. If you set `MIXEDGP_R_LIBRARY`,
retain that setting for both setup and execution; an existing repository
`R-library/` is also recognized. These comparison packages remain optional
(`Suggests`) for package-only use. Installation does not backfill missing
comparisons in existing cached results.

### 3. Set the computing budget

On macOS/Linux, enter the number of cores you want to allocate; for example,
12 for a 12-core Mac mini, or replace 12 with 16 for a 16-core MacBook:

```sh
export MIXEDGP_CORE_BUDGET=12
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
```

Set these before R starts. They limit numerical-library threads where supported;
they do not switch the BLAS implementation. Development has three datasets per
setting, so at these budgets it uses at most three dataset workers, each with
four parallel chains. Settings are processed sequentially. The current Windows
execution path falls back to serial; the shell commands above target macOS/Linux.

### 4. Preview both studies without starting MCMC

```sh
Rscript --vanilla experiments/run_development_study.R study1 plan
Rscript --vanilla experiments/run_development_study.R study2 plan
```

Check the plans against [NUMERICAL_DESIGN.md](NUMERICAL_DESIGN.md):

- Both studies: 100 training and 100 test observations, three replications per setting.
- Study I: balanced categories with eta = 0/1; calibration 0, 20, 50 in each setting. Imbalanced settings are excluded from new runs.
- Study II: primary q = 2 at calibration 50; primary q = 4 at 0, 20, 50, 80; logistic misspecification q = 4 at 50.

A replication is a separately generated synthetic dataset, not an MCMC chain.
Calibration observations are part of the 100 training observations, not extra
training units.

### 4b. Prepare reusable competitor results (no MCMC)

After installing dependencies, run competitors as a separate experiment stage.
This stage loads frozen data, fits UC-GP/LVGP/EzGP, saves each method separately,
and exports paper-facing tables. It never installs packages or runs EIV-GP.

Use the publication frozen collection in `reproduction/data/synthetic/study1`
and `study2`. If you have not obtained those data, generate them once using:

```sh
Rscript --vanilla experiments/run_publication_study.R study1-data
Rscript --vanilla experiments/run_publication_study.R study2-data
```

Existing incompatible datasets must be inspected and preserved, not overwritten.
For exact reproduction, obtain the matching frozen collection and code revision;
generated data and caches under `reproduction/` are Git-ignored, not automatically
distributed by cloning the repository.

Set an existing paper project folder (your local Overleaf/Dropbox folder is fine):

```sh
export EIVGP_OVERLEAF_ROOT="/absolute/path/to/your/paper"
export EIVGP_WORKERS=10
Rscript --vanilla experiments/run_competitors.R study1 plan publication
Rscript --vanilla experiments/run_competitors.R study1 run publication
Rscript --vanilla experiments/run_competitors.R study2 plan publication
Rscript --vanilla experiments/run_competitors.R study2 run publication
```

`plan` only inspects the setup. `run publication` processes 50 datasets per
setting, one dataset per worker, with its methods run serially. Keep the BLAS
thread limits above; use fewer workers if memory is tight. The entire selected
frozen collection is verified before fitting. Progress identifies each dataset
and method; successful fits and failed attempts are checkpointed separately.

Both modes use the same competitor optimization settings and publication data
root. To extract the first three datasets per setting without fitting anything:

```sh
Rscript --vanilla experiments/run_competitors.R study1 export development
Rscript --vanilla experiments/run_competitors.R study2 export development
```

Alternatively, `run development` prepares only that subset using the same fit
protocol. A later publication run reuses those fits. `retry publication` retries
failed or missing fits, preserving successful fits and previous failed attempts.
Ordinary `run` does not repeatedly retry a saved failure.

Outputs:

- `reproduction/competitor-cache/`: fitted objects, predictive moments, status,
  settings and provenance, saved per method/dataset.
- `reproduction/competitor-reports/<study>/<mode>/`: status, prediction and
  per-replication metric CSVs, plus means/Monte Carlo standard errors and provenance.
- `$EIVGP_OVERLEAF_ROOT/tables/competitors/<study>/<mode>/`: exported tables,
  including `predictive_summary.tex`. Check `statuses.csv` and `n_success`
  before using any summary; unsuccessful fits are not silently scored as zero.

These are prediction comparisons for `Y* | X*, C*`, not physical latent-input
recovery. Competitors do not use calibration measurements, so a fitted predictor
is reused across calibration sizes. Cache identity checks training/prediction
inputs, seed, method settings, adapter source, R/platform and package version;
it does not depend on mode, calibration size or number of predictive draws.
Changed identities require a new fit. Current publication optimization budgets
are retained, not newly certified as final: LVGP's limit is **per attempt**
(1,800 seconds in Study I, 3,600 in Study II, up to three attempts).

`EIVGP_DATA_ROOT` can select another matching frozen collection;
`EIVGP_COMPETITOR_CACHE` can select a shared/copied cache directory;
`EIVGP_COMPETITOR_METHODS` can select methods, e.g. `UC-GP,EzGP`.
Cache reuse on another machine requires matching recorded runtime/package
identities. Unset a study-specific `EIVGP_DATA_ROOT` before switching studies.

The Study I/II MCMC experiment drivers neither fit competitors nor read their
caches. They save dataset identities alongside their own metrics and diagnostics.
Competitors run independently, using the same frozen datasets, on another machine
if desired. Only reporting needs access to both sets of outputs.

### Combine independent results for the paper

After fitting, copy/sync the MCMC run folder and competitor report folder to one
machine. No competitor model cache is needed for this step. For example:

```sh
Rscript --vanilla experiments/combine_results.R \\
  "/path/to/mcmc-run" \\
  "/path/to/competitor-reports/study1/publication" \\
  "$EIVGP_OVERLEAF_ROOT/tables/combined/study1"
```

Use the exact run directory containing `config/resolved_config.rds` and `cells/`.
The competitor directory must contain `metrics.csv`, `statuses.csv` and
`provenance.rds`. The merger checks study and per-dataset content checksums,
selects the MCMC run's settings/replications, and rejects mismatched data or
duplicate metric rows. Publication competitor results can therefore serve a
three-replication development MCMC run without refitting.

This command writes overall response-prediction summaries, matched-replication
EIV-minus-competitor differences and MCSEs, availability/status records, and
MCMC diagnostics. Baselines retain `n_calib=NA`: they are not independent fits
at each calibration size. Coverage/width differences are descriptive, not
universally lower-is-better. Missing results remain explicit, never zero-filled.
Other tasks (latent recovery, response surfaces and ablations) retain their
existing study-specific reporting; this merger does not fabricate competitor
results for those tasks or combine incompatible evaluation strata.

New MCMC runs save `dataset_identity.csv`. Older runs without that evidence
are deliberately rejected by the merger; they require a separate provenance
audit, not automatic relabeling or a forced MCMC restart. Reporting does not
change the original fits, and fitting jobs do not wait for synchronization.

### 5. Run one study at a time

Run Study I, then inspect its outputs before starting Study II:

```sh
Rscript --vanilla experiments/run_development_study.R study1 run
```

```sh
Rscript --vanilla experiments/run_development_study.R study2 run
```

The visible wrapper is
[`experiments/development_numerical_experiment.Rmd`](experiments/development_numerical_experiment.Rmd).
It generates and freezes development datasets, then fits and aggregates results.
Subsequent runs verify those data. If existing data fail verification, preserve
them and resolve the mismatch before proceeding; do not silently overwrite or
relabel older datasets.

Each chain runs 500 warmup plus 1,250 sampling iterations: 1,750 total per chain,
with 5,000 post-warmup draws retained across four chains. Continuation is opt-in,
not automatic. These runs can still take substantial time. The knitting progress
percentage measures document chunks, not completed datasets or remaining MCMC
time; a document can remain at its fitting chunk for a long time.

### Run fitting and reporting independently

For long experiments, you can finish fitting/evaluation without generating
figures or publication tables:

```sh
Rscript --vanilla experiments/run_development_study.R study1 fit
Rscript --vanilla experiments/run_development_study.R study2 fit
```

The existing `run` action fits/evaluates all cells first, then aggregates and
reports. Each completed cell saves `report_inputs.rds` before reporting begins.
The fitting stage still computes predictive/evaluation metrics; those costly
calculations are not repeated by the reporting commands below.

To regenerate tables and figures, use the exact saved run directory printed by
the fitting command. Replace the example path with yours:

```sh
RUN_DIR="reproduction/development/results/study2/study2-development-YOUR_RUN_ID"
Rscript --vanilla experiments/report_study.R "$RUN_DIR" summarize
Rscript --vanilla experiments/report_study.R "$RUN_DIR" plot
```

Use `report` instead of `summarize`/`plot` to generate both. These commands do
not fit models, generate datasets, or select a run automatically. `plot` reads
saved evaluation results and computes plotting summaries in memory; it does not
require the `summarize` command first. New files go under
`RUN_DIR/reporting/<cell>/tables/` and `figures/`, with study-wide summaries
under `RUN_DIR/reporting/combined/`, and status and provenance
files. Plot-only runs leave table files unchanged. Empty plots are skipped and
recorded in `plot_status.csv`; other reporting failures are recorded separately
from fitting failures. Incomplete reporting returns a nonzero CLI exit status,
but never triggers refitting or deletes cached results.

Previously saved Study II `study2_results_*.rds` bundles can also be reported
without rerunning MCMC, provided the run's saved configuration is present and
each cell has a unique bundle. Cells with no saved results are explicitly
marked `missing_inputs`, not silently treated as complete. Older Study I runs
without `report_inputs.rds` are not automatically migrated. Preserve all old
run directories: this refactor changes fitting code identity, so a fresh `fit`
or `run` may use a new directory rather than reuse an older run automatically.

### 6. Inspect results and diagnostics

Outputs stay inside the repository checkout:

```text
reproduction/development/
├── data/
├── results/
└── reports/
```

Start with the study HTML reports in `reports/`. A `fit` action writes
`study1_development_fit.html` or `study2_development_fit.html`. The report identifies its run
directory; inspect its combined outputs, missing-method records, and
`config/diagnostic_gates.csv`. Review R-hat, ESS, Monte Carlo errors, and diagnostic
advice alongside estimates and figures. Diagnostic warnings do not automatically
stop execution or erase completed estimates; genuine runtime errors can still
stop a run. Three replications check the workflow and offer preliminary insight,
but do not establish precise simulation-performance comparisons.

No output is automatically copied into Overleaf or pushed to GitHub. Paths are
relative to the checkout and require no author-specific local folders.

### Moving to publication reproduction

Publication uses 50 replications per setting. Both studies use the following
EIV-GP sampling budget:

| Setting | Chains | Warmup per chain | Retained per chain | Total iterations per chain |
|---|---:|---:|---:|---:|
| Study I and Study II | 4 | 5,000 | 15,000 | 20,000 |

All 60,000 post-warmup draws are kept (`thin = 1`). Development and smoke
budgets, measurement-only comparator budgets, and evaluation/integration
budgets are unchanged. Keeping every sampling draw does not imply evaluating
every draw in every prediction routine; those separate budgets remain recorded
in the resolved configuration.

```bash
Rscript --vanilla experiments/run_publication_study.R study1-data
Rscript --vanilla experiments/run_publication_study.R study2-data
Rscript --vanilla experiments/run_publication_study.R study1
Rscript --vanilla experiments/run_publication_study.R study2
```

The launcher defaults to publication mode. `EIVGP_WORKERS` controls replication
workers, not chain count. Use a worker count appropriate to memory and avoid
nested parallel oversubscription. Existing frozen data are verified rather than
silently replaced; changed fitting configurations receive distinct run identities.

The budget is a starting protocol, not a convergence guarantee. Runs finish
with diagnostic warnings rather than requiring a convergence gate, and do not
automatically extend. Review target-specific R-hat (reference 1.01), bulk/tail
ESS (reference 400), MCSE, and traces, particularly for variance summaries and
latent imputation. Extend selected fits explicitly when warranted. Numerical
failures remain errors, not mixing warnings.

#### Interruptions and restart

Results are saved atomically after each complete synthetic-data replication
within a design cell (including its configured calibration settings and
evaluation), not every 500 MCMC iterations. Rerun the same command with the
same source, configuration, and data to validate and reuse completed caches.
An unfinished replication restarts from its beginning; completed replications
do not. This is dataset-level recovery, not within-chain continuation.

Each cell's `tables/study1_replication_status.csv.tasks/` or
`tables/study2_replication_status.csv.tasks/` contains one live status file per
task: `pending`, `running`, `success`, or `failed`. The adjacent status CSV
summarizes the batch after workers return. After a killed process, `running`
may be stale; the validated result cache, not the status label, determines reuse.
Do not run two launchers concurrently against the same output directory.
Diagnostic warnings do not invalidate a completed result. Changes to priors,
code, or fitting configuration require a distinct compatible run identity;
old-prior results are never relabeled as new-prior fits.

This budget change does not change priors or promote experimental input--signal
moves into the default sampler. Code and manuscript now use Beta(13,3).
The saved continuation used Beta(32,8), and the earlier joint-move Study I pilots
used Beta(13.052446,2.794765); neither validates the current prior. First check one
representative fit per study at the publication budget for runtime and mixing.

Before final tables are reported, audit the frozen collection and evaluation
budgets. Development and publication still have separate default data directories.
Fits from before 0.3.0 must be refitted under the new posterior; 0.3.0 fits cannot
be relabeled as results from the 0.3.1 algorithm. See [DEVELOPMENT.md](DEVELOPMENT.md)
for remaining design limitations and diagnostic handling.

## Fit, predict, and diagnose

`X` is a numeric predictor matrix, `y` a response vector, and `C` an ordered
factor or a data frame of ordered factors. `U_obs` contains calibrated latent
inputs, with `NA` for unobserved entries.

```r
library(eivGP)
fit <- fit_eivgp(
  X, y, C, U_obs,
  engine = "multivariate", latent_dim = 2L, ident = "none",
  kernel = "matern", n_chains = 4L,
  n_iter = 1750L, burn = 500L, seed = 123L
)
summary(fit)
mean_draws <- predict(fit, new_X = X_new, new_C = C_new, target = "mean")
response_draws <- predict(fit, new_X = X_new, new_C = C_new, target = "response")
surface_draws <- predict(fit, new_X = X_grid, new_U = U_grid, target = "surface")
latent_draws <- impute_eivgp(fit, new_C = C_new)
report <- diagnose_eivgp(fit, X = X_panel, C = C_panel)
report$table
```

Choose `engine = "univariate"` for a single deterministically thresholded
latent input; omit `latent_dim` and `ident` in that case. The multivariate
engine uses a noisy ordinal-probit measurement model. Calibration and loading
constraints determine the interpretation of the latent scale; inspect
`summary(fit)` and `fit$model_specification` before interpreting physical units.

The default budget is four chains of 500 warmup and 1,250 retained transitions
per chain. These general fitting controls are tunable. Every post-warmup draw
is kept (`thin = 1`); no HMC is used. Fits finish their requested budget and
report diagnostic warnings. Review target-specific R-hat, ESS, and MCSE before
scientific interpretation; numerical failures remain errors.

Continuation is explicit and requires a compatible 0.3.1 checkpoint:

```r
saveRDS(fit, "fit-checkpoint.rds")
fit <- continue_eivgp(readRDS("fit-checkpoint.rds"), n_iter = 1000L)
```

This adds transitions to each chain without repeating warmup. Additional draws
can improve Monte Carlo precision; poor mixing also requires investigation.

## Posterior and repository experiments

The current sampler integrates out the common variance and uses
`V ~ IG(3,2)`, a signal fraction `r ~ Beta(13,3)`, continuous numeric-input
kernel coefficients, a finite dictionary of latent-input kernel coefficients,
and Dirichlet priors
on ordinal category probabilities. Conditional variance recovery preserves the
joint posterior used by prediction. These priors differ from earlier releases.
The fraction `r` is pointwise signal variance, not realized regression R-squared.
`X` and `y` are standardized using training statistics. Stored `samples_V` and
the variance prior use standardized response units; multiply by
`fit$data$y_scale^2` to recover variance in response units. Dictionary entries
are inverse squared-distance coefficients on model-scaled inputs, not lengths.
Set `priors` and `sampler_control` explicitly to change the model or algorithm;
for example, a two-dimensional coordinate dictionary can be supplied as
`priors = list(u_dictionary = list(c(.1, .5, 2), c(.1, .5, 2)))`.
`sampler_control = list(u_block_size = 8L, cross_every = 10L)` sets the reference
blocking schedule. The help for `fit_eivgp` lists all configurable fields.

Study designs, replications, dataset workers, evaluation budgets, and output
paths belong to the repository experiment layer. They are not installed as
package APIs. See `DEVELOPMENT.md` and the scripts in `experiments/`; this
release does not regenerate study results or establish publication readiness.
The package retains optional reusable `fit_ucgp()`, `fit_lvgp()`, and
`fit_ezgp()` competitor adapters. Run `citation("eivGP")` for the manuscript
citation.
# Unified diagnostic exports

After `d <- diagnose_eivgp(fit)`, use
`write_diagnostics_eivgp(d, "diagnostics")` to save summary, raw-parameter,
target/invariant, and complete-detail CSVs plus the full diagnostic report RDS.
Existing exports require explicit `overwrite = TRUE`. NA diagnostics remain NA;
exporting does not refit, thin, extend chains, or certify convergence.

Both study fitting drivers and saved-result reporting use the same table writer.
The stable names are `studyN_mcmc_diagnostics.csv`,
`studyN_mcmc_parameter_diagnostics.csv`, `studyN_mcmc_target_diagnostics.csv`, and
`studyN_mcmc_diagnostic_details.csv`. Parameter files now contain raw coordinates
only; the details file preserves the combined table formerly exported by Study I.
Existing Study II cache-tagged reports remain available for compatibility.
Study II `measurement_parameter_diagnostics` refers to the optional response-free
measurement comparator, not EIV-GP; an empty comparator table does not imply
missing EIV-GP diagnostics. Export normalization preserves study-specific gates
and recorded draw windows; it does not retrospectively apply package thresholds.
