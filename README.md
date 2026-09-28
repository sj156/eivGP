# eivGP 0.3.1

`eivGP` fits Gaussian process regression with numeric predictors, ordinal proxies
for latent continuous inputs, and optional observed calibration inputs. It supports
univariate threshold and multivariate ordinal-probit measurement models, prediction,
latent-input imputation, MCMC diagnostics, and explicit continuation.

The installable package is in `eivGP/`; `litr/` builds it from `codes/core/`.
Supporting R code is grouped by purpose under [`codes/`](codes/README.md).

## Read an example, then run the same workflow

There are **two analysis notebooks**. Each starts with a small exploratory
analysis and computes its displayed results automatically. The bundled HTML
uses reduced computation to demonstrate the workflow.

| Workflow | Notebook | Worked HTML example |
| --- | --- | --- |
| Numerical experiments, Studies I and II | [Notebook](replication/01_numerical_experiments.Rmd) | [HTML](docs/replication/01_numerical_experiments.html) |
| GLODAP: silicate and partially observed oxygen | [Notebook](replication/02_ocean.Rmd) | [HTML](docs/replication/02_ocean.html) |

After installing dependencies, open a notebook and Knit, or call
`rmarkdown::render("replication/01_numerical_experiments.Rmd")` from the repository
root. Defaults are `profile: "quick"` and `action: "run"`. Use `profile: "paper"`
in the **same notebook** for the full computation. The numerical notebook also
offers `development`, and `study` selects either or both studies. See the
[notebook guide](replication/README.md) for data requirements and saved-result reuse.

The paper numerical configuration is 50 datasets per setting, four chains,
20,000 total iterations per chain, and 5,000 warmup. The quick HTML illustrates
analysis; short-chain scores are not publication estimates or convergence evidence.
The ocean notebook lists its full budget and continuation protocol. The GLODAP
frozen subset is bundled. ADNI is deferred from this reader-facing collection.

`replication/` contains notebook sources; `docs/replication/` contains the selected
rendered examples; `reproduction/` and application output folders contain local
fits, logs, and caches. Preprocessing and operational scripts are supporting code,
not separate reader notebooks.

## Install and use the package

Run from a repository checkout. Requires R >= 4.1; experiment reports also need
Pandoc (available with RStudio).

```sh
git clone https://github.com/sj156/eivGP.git
cd eivGP
Rscript --vanilla -e 'install.packages(c("posterior", "TruncatedNormal"), repos="https://cloud.r-project.org")'
R CMD INSTALL eivGP
```

For your own data, `X` is a numeric matrix, `y` a response vector, `C` an ordered
factor (or data frame of ordered factors), and `U_obs` contains calibrated inputs
with `NA` for unobserved values. A univariate example:

```r
library(eivGP)
fit <- fit_eivgp(X, y, C, U_obs, engine = "univariate",
                n_chains = 4L, n_iter = 1750L, burn = 500L, seed = 123L)
pred <- predict(fit, new_X = X_new, new_C = C_new, target = "response")
diagnose_eivgp(fit)$table
```

Use `engine = "multivariate"` for noisy ordinal proxies of multivariate latent
inputs; see `?fit_eivgp` for calibration, identification, priors, and kernel options.
Use `impute_eivgp()` for latent inputs and `continue_eivgp()` for a compatible saved
fit. More iterations do not by themselves establish convergence.

## Reproduce the numerical experiments

Run all commands below from the repository root. On an existing checkout, skip
cloning. Keep the Git commit and saved configuration with your results: package
version alone does not identify the priors or experiment settings.

### 1. Set up each machine

```sh
Rscript --vanilla codes/cli/install_eivgp_dependencies.R
R CMD INSTALL eivGP
git rev-parse HEAD

export MIXEDGP_CORE_BUDGET=16
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
```

Use `12` for a 12-core Mac mini or a smaller budget if memory is limited.
Dependencies are installed explicitly, never during fitting. The shell examples
target macOS/Linux; Windows parallel execution falls back to serial.

### 2. Design and frozen data

Both studies use **100 training and 100 test observations** per dataset.

| Study | Settings | Calibration sizes |
|---|---|---|
| I | Balanced η = 0 and η = 1 | 0, 20, 50 |
| II | Primary q = 2 | 50 |
| II | Primary q = 4 | 0, 20, 50, 80 |
| II | Logistic misspecification q = 4 | 50 |

Here `q` is the number of ordinal proxies; Study II has two latent dimensions.
Calibration uses a nested subset of training observations, not additional data.

| Mode | Datasets per setting | Chains | Warmup + sampling per chain |
|---|---:|---:|---:|
| Development | 3 | 4 | 500 + 1,250 |
| Publication | 50 | 4 | 5,000 + 15,000 |

These are the current main EIV-GP budgets for both studies. Evaluation and
ablation budgets are separate; inspect the printed plan and saved configuration.

Generate publication data once, before publication MCMC or competitor runs:

```sh
Rscript --vanilla codes/cli/run_publication_study.R study1-data
Rscript --vanilla codes/cli/run_publication_study.R study2-data
```

Files go to `reproduction/data/synthetic/`. Development MCMC currently prepares
its own `reproduction/development/data/` collection using the same per-replication
generator settings and seeds. Regeneration requires unchanged code/design/seeds;
incompatible existing data must be resolved, not silently relabeled. Across
machines, use matching frozen data. `reproduction/` is Git-ignored, so cloning
does not download data or results.

### 3. Run MCMC independently

Start with development mode; run one study at a time:

```sh
Rscript --vanilla codes/cli/run_development_study.R study1 plan
Rscript --vanilla codes/cli/run_development_study.R study1 fit
Rscript --vanilla codes/cli/run_development_study.R study2 plan
Rscript --vanilla codes/cli/run_development_study.R study2 fit
```

`fit` saves MCMC and evaluation results without plotting. Use `run` instead to
also report. Neither action fits competitors or reads their caches.

For publication runs, after reviewing development diagnostics:

```sh
EIVGP_RUN_MODE=dry_run Rscript --vanilla codes/cli/run_publication_study.R study1
EIVGP_RUN_MODE=publication Rscript --vanilla codes/cli/run_publication_study.R study1
EIVGP_RUN_MODE=dry_run Rscript --vanilla codes/cli/run_publication_study.R study2
EIVGP_RUN_MODE=publication Rscript --vanilla codes/cli/run_publication_study.R study2
```

Keep the exact run-directory path printed at completion. Development results are
under `reproduction/development/results/`; publication results under
`reproduction/results/`. HTML reports identify the corresponding run.

### 4. Run competitors independently

This can run on another machine while MCMC is running. Set an **existing** paper
folder; for the author's setup, use the local Dropbox/Overleaf project folder.

```sh
export EIVGP_OVERLEAF_ROOT="/absolute/path/to/your/paper"
export EIVGP_WORKERS=10
Rscript --vanilla codes/cli/run_competitors.R study1 plan publication
Rscript --vanilla codes/cli/run_competitors.R study1 run publication
Rscript --vanilla codes/cli/run_competitors.R study2 plan publication
Rscript --vanilla codes/cli/run_competitors.R study2 run publication
```

Each worker handles one dataset, fitting UC-GP, LVGP, and EzGP sequentially.
Compatible fits are cached; `retry` retries failed fits, while `export development`
extracts the first three datasets without fitting. Competitors do not refit for
each calibration size. LVGP's time limit is per attempt: 30 minutes in Study I,
60 minutes in Study II, up to three attempts.

Reports are saved under `reproduction/competitor-reports/<study>/<mode>/` and
exported to `$EIVGP_OVERLEAF_ROOT/tables/competitors/<study>/<mode>/`.

### 5. Inspect and combine results

Generate study tables and figures from a saved MCMC run, without refitting:

```sh
RUN_DIR="/absolute/path/to/the/saved/mcmc-run"
Rscript --vanilla codes/cli/report_study.R "$RUN_DIR" report
```

Review R̂, bulk/tail ESS, MCSE, traces, and failure records before interpreting
tables. Diagnostic warnings preserve results; runtime errors can still stop a run.
Three development datasets give preliminary evidence, not precise comparisons.

To combine independent outputs, copy/sync the MCMC run folder and competitor
report folder onto one machine. No competitor model cache is needed:

```sh
COMPETITOR_REPORT="/absolute/path/to/competitor-reports/study1/publication"
Rscript --vanilla codes/cli/combine_results.R "$RUN_DIR" "$COMPETITOR_REPORT" "$EIVGP_OVERLEAF_ROOT/tables/combined/study1"
```

Use matching study paths. The merger verifies dataset checksums and writes
overall predictive summaries, paired differences, diagnostics and availability
records. Other tasks retain their study-specific reports. Older MCMC runs without
dataset identity records require a provenance audit before merging.

**Restarting:** rerun the same fitting command with unchanged code, settings,
and data to reuse validated completed replications. An unfinished replication
restarts from its beginning; this is not within-chain continuation. Do not launch
two jobs into the same output directory or remove locks while workers are active.

For implementation details, see [DEVELOPMENT.md](DEVELOPMENT.md) and the
[literate build instructions](litr/README.md).
