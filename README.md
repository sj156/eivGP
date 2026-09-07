# eivGP 0.3.0

`eivGP` fits Gaussian process regressions with numeric predictors, ordinal
proxies for latent continuous inputs, and optional calibration measurements.
The reusable package and the manuscript's numerical experiments are separate:
`eivGP/` is the single installable package; `codes/` is its canonical source;
`experiments/` and the repository simulation helpers define the studies.

Version 0.3.0 restores the package name `eivGP` and replaces the previous
posterior prior and sampler. Fits and checkpoints from `eivGP` 0.1.x or
`eivmixgp` 0.1–0.2.x must be refitted. Renaming or relabeling an old fit does
not migrate it to the new posterior.

## Install

The current seven-setting numerical design (100 training / 100 test observations,
100 publication versus 3 development replications) is documented in
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
R CMD check --no-manual eivGP_0.3.0.tar.gz
```

## Reproduce the numerical experiments: a new-reader workflow

You can currently run the revised experiments in **development mode**. This is
not yet exact reproduction of final paper tables: the revised publication
results have not been generated, publication MCMC budgets remain undecided,
and shared frozen-data storage between modes is not yet implemented.
Development uses the same seven settings and sample sizes, but only three
replicated datasets per setting instead of 100. Real-data reproduction
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
Rscript --vanilla -e 'install.packages(c("posterior", "TruncatedNormal", "rmarkdown", "knitr", "ggplot2", "dplyr", "tidyr", "patchwork", "kergp", "LVGP", "EzGP"), repos="https://cloud.r-project.org")'
R CMD INSTALL eivGP
Rscript --vanilla -e 'library(eivGP); print(packageVersion("eivGP")); stopifnot(rmarkdown::pandoc_available())'
```

The package version for this workflow is 0.3.0. If dependency installation
fails, resolve the reported error before continuing. Missing optional
competitors (`kergp`, `LVGP`, or `EzGP`) are recorded and skipped, so their
comparisons will be incomplete. If Pandoc is not detected from Terminal,
make it available there before rendering the experiment documents.

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
- Study I: balanced/imbalanced categories crossed with eta = 0/1; calibration 0, 10, 50 in each setting.
- Study II: primary q = 2 at calibration 50; primary q = 4 at 0, 20, 50, 80; logistic misspecification q = 4 at 50.

A replication is a separately generated synthetic dataset, not an MCMC chain.
Calibration observations are part of the 100 training observations, not extra
training units.

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

### 6. Inspect results and diagnostics

Outputs stay inside the repository checkout:

```text
reproduction/development/
├── data/
├── results/
└── reports/
```

Start with the study HTML reports in `reports/`. The report identifies its run
directory; inspect its combined outputs, missing-method records, and
`config/diagnostic_gates.csv`. Review R-hat, ESS, Monte Carlo errors, and diagnostic
advice alongside estimates and figures. Diagnostic warnings do not automatically
stop execution or erase completed estimates; genuine runtime errors can still
stop a run. Three replications check the workflow and offer preliminary insight,
but do not establish precise simulation-performance comparisons.

No output is automatically copied into Overleaf or pushed to GitHub. Paths are
relative to the checkout and require no author-specific local folders.

### Moving to publication reproduction

Publication uses 100 replications per setting. Before final tables can be
reproduced, we must finalize the publication MCMC/evaluation budgets, prepare and
audit the revised frozen collection, and make development use a prescribed
subset of that same collection. Currently the two modes have separate default
data directories. Old fits must be refitted under 0.3.0; they cannot be relabeled
as new results. Do not interpret the existing publication launcher's numeric
defaults as a finalized paper protocol. See [DEVELOPMENT.md](DEVELOPMENT.md) for
current limitations and diagnostic handling.

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

Continuation is explicit and requires a compatible 0.3.0 checkpoint:

```r
saveRDS(fit, "fit-checkpoint.rds")
fit <- continue_eivgp(readRDS("fit-checkpoint.rds"), n_iter = 1000L)
```

This adds transitions to each chain without repeating warmup. Additional draws
can improve Monte Carlo precision; poor mixing also requires investigation.

## Posterior and repository experiments

The current sampler integrates out the common variance and uses
`V ~ IG(3,2)`, a signal fraction `r ~ Beta(32,8)`, continuous numeric-input
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
