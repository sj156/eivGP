# eivGP

Gaussian process regression with ordinal proxies and partially observed continuous
inputs. Use the R package for your own data, or the scripts to reproduce Studies I
and II and the Ocean, ADNI and OAI applications.

The workflow is **install → prepare inputs → plan → run → inspect results**.
Run the commands below from the repository root. Use absolute input/output paths
when launching scripts from another directory.

## 1. Install and verify

Use a working R package compiler. The locked paper environment uses R 4.6.0;
the package supports R >= 4.1. HTML reports also require Pandoc.

Set numerical-library threads to one, then restore the locked dependencies and
install this checkout into `.r-library/`:

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
Rscript tools/setup_numerical.R
Rscript tools/setup_numerical.R --check
```

The first command restores `renv.lock` and installs the package. `--check`
verifies dependency versions and that the installed package matches the source.
The study and application scripts find `.r-library/` automatically. Set
`EIVGP_R_LIBRARY` if you use a different library.

After editing package sources, refresh their identities and reinstall using the
already-restored dependencies:

```sh
Rscript tools/build_package.R
Rscript tools/setup_numerical.R --install-local
```

Restart existing R sessions after reinstalling. Update the library when no active
jobs are using it.

Basic inference uses Rfast, Rcpp, posterior and TruncatedNormal; compilation also
uses RcppArmadillo. Plotting/report dependencies are optional for the package.
The setup script installs the full environment needed for paper reproduction.

## 2. Walk through the R API

`X` is an `n × p` numeric matrix, `y` is a length-`n` response, and `C` is an
`n × q` matrix of ordinal codes starting at 1. `U_obs` is an `n × d` matrix:
observed calibration rows contain measured inputs; missing rows contain `NA`.
Each latent-input row must be either fully observed or fully missing.

### Fit a small example

Start R in the repository root and load the installed checkout:

```r
.libPaths(c(normalizePath(".r-library"), .libPaths()))
library(eivGP)

set.seed(1)
n <- 40L
X <- cbind(x = runif(n))
U <- cbind(u = rnorm(n))
C <- cbind(proxy = findInterval(U[, 1], c(-Inf, -0.5, 0.5, Inf)))
y <- sin(2 * pi * X[, 1]) + U[, 1] + rnorm(n, sd = 0.2)
U_obs <- U
U_obs[11:n, ] <- NA_real_  # Only the first ten inputs are measured.

fit <- fit_eivgp(
  X, y, C, U_obs = U_obs, engine = "univariate", m_vec = 3L,
  sampler = "ESS", n_iter = 2000L, burn = 500L,
  n_chains = 4L, n_cores = 12L
)
summary(fit)
```

`fit_eivgp()` fits the response surface and latent-input model together. This
example demonstrates the interface; its iteration budget does not establish
convergence. For several latent inputs, use `engine="multivariate"` and set
`latent_dim=d`; specify the category counts in `m_vec`.

### Diagnose, predict and impute

```r
checks <- diagnose_eivgp(fit)

new_X <- cbind(x = c(0.25, 0.75))
new_C <- cbind(proxy = c(1L, 3L))
response_draws <- predict_eivgp(
  fit, new_X = new_X, new_C = new_C, target = "response", seed = 2L
)
imputed_inputs <- impute_eivgp(fit, rows = 11:15)
```

| Function | What to use it for |
| --- | --- |
| `fit_eivgp()` | Fit the model and retain chain histories and terminal states. |
| `diagnose_eivgp()` | Review R-hat, bulk/tail effective sample size, MCSE and diagnostic status. Supply `X`, `C` and, when available, `U` to monitor a scientific target panel. |
| `predict_eivgp()` | Draw responses (`target="response"`), observable conditional means (`"mean"`), or the latent-input surface (`"surface"`; supply `new_U`). |
| `impute_eivgp()` | Recover missing training inputs; supply `new_C` instead of `rows` for prospective inputs. |
| `save_fit()` / `read_fit()` | Save and restore a fit with its continuation state. |
| `continue_eivgp()` | Explicitly add iterations from saved states without repeating warmup. |

Inspect diagnostics for the quantities you intend to report. Predictive accuracy,
parameter mixing and physical-input recovery are different assessments.

### Save and, when justified, continue

```r
dir.create("results", showWarnings = FALSE)
save_fit(fit, "results/example-fit.rds")
restored <- read_fit("results/example-fit.rds")

# Optional: extend after reviewing diagnostics.
# extended <- continue_eivgp(restored, n_iter = 1000L, n_cores = 12L)
```

ESS means **elliptical slice sampling** when selecting a sampler; in diagnostic
reports, ESS means **effective sample size**. ESS is the package's default
sampler. HMC and the HMC–ESS hybrid are also available. The application scripts
choose their own documented samplers by calibration setting. See the
[package guide](eivGP/README.md) for sampler restrictions, tuning, scales and
prediction targets.

## 3. Reproduce Studies I and II

The paper workflows have separate EIV-GP and competitor scripts. They share
frozen inputs and fixed seeds, and can run on different machines. Choose one
study and use the same input path for its two method groups.

For Study I, work through these steps:

```sh
# 1. Inspect the complete EIV-GP workload; nothing is fitted.
Rscript tools/study1-eivGP.R --action=plan --cores=12 --data-root=results/data/study1 --output=results/study1-eivGP

# 2. Generate/verify the fixed-seed inputs; no MCMC is run.
Rscript tools/study1-eivGP.R --action=data --cores=12 --data-root=results/data/study1 --output=results/study1-eivGP

# 3. Fit EIV-GP, then fit the original and complete-case competitors.
Rscript tools/study1-eivGP.R --action=run --cores=12 --data-root=results/data/study1 --output=results/study1-eivGP
Rscript tools/study1-competitors.R --action=run --cores=12 --data-root=results/data/study1 --output=results/study1-competitors

# 4. Rebuild the EIV-GP report from saved results.
Rscript tools/study1-eivGP.R --action=report --cores=12 --data-root=results/data/study1 --output=results/study1-eivGP
```

For Study II, replace `study1` with `study2` in these commands and paths. Its EIV-GP
script includes the FULL/CUT feedback comparison; calibrated FULL scoring reuses
saved response fits. The detailed guides explain the separate original-grid and
calibrated-comparison budgets, appendix subsets and saved-run selection.

| Action | Effect |
| --- | --- |
| `plan` | List the requested workload and settings. This is the study-script default. |
| `data` | Prepare or verify frozen synthetic inputs. |
| `run` | Fit/evaluate the selected workflow, reusing compatible completed work. |
| `report` | Read saved results and regenerate reports without fitting. |

Keep compatible inputs/settings when resuming work. Changed sources or settings
need a separate output root; archived results retain their original provenance.
The notebooks in `paper/` read saved results.

Start with [Study I reproduction](paper/study1-README.md) or
[Study II reproduction](paper/study2-README.md). For custom experiments, the
lower-level package functions are `study_config()` and `run_numerical_study()`;
`read_study_results()` and `summarize_comparison()` load and summarize saved scores.

## 4. Run an application

Each application uses the same three scripts: `APP-eivGP.R` fits/evaluates EIV-GP,
`APP-competitors.R` fits/evaluates UC-GP, LVGP and EzGP comparisons, and
`APP-report.R` combines their saved results.

### Choose and prepare the inputs

| Application / script prefix | Prepared inputs | Example output root | Guide |
| --- | --- | --- | --- |
| Ocean / `ocean` | `data/ocean_silicate/n500/` | `results/ocean500/` | [Ocean](paper/README_ocean.md) |
| ADNI / `adni` | `data/adni/` | `results/adni/` | [ADNI](paper/README_adni.md) |
| OAI / `oai` | `data/oai/` | `results/oai/` | [OAI](paper/README_oai.md) |

Ocean500 is the default paper cohort; its frozen folds are supplied. Ocean200 is
also supplied. Prepare the raw Ocean1000 cohort separately with
`Rscript tools/ocean-prepare-folds.R ocean1000`.

ADNI and OAI require privately supplied prepared folds, metadata and observation
masks. Follow the [ADNI data contract](data/adni/README.md) or
[OAI data contract](data/oai/README.md). Participant records are not distributed
in Git. Keep participant inputs, fits and empirical reports local. The
`*-prepare-folds.R` scripts prepare new cohorts under those contracts;
existing frozen folds should be reused.

### Walk through Ocean

```sh
# 1. Confirm the workload and resource allocation.
Rscript tools/ocean-eivGP.R --action=plan --profile=paper --cores=12 --data-root=data/ocean_silicate/n500 --output=results/ocean500

# 2. Fit/evaluate EIV-GP and its comparisons, using the same inputs/output root.
Rscript tools/ocean-eivGP.R --action=run --profile=paper --cores=12 --data-root=data/ocean_silicate/n500 --output=results/ocean500
Rscript tools/ocean-competitors.R --action=run --profile=paper --cores=12 --data-root=data/ocean_silicate/n500 --output=results/ocean500

# 3. Assemble the combined report after both method groups finish.
Rscript tools/ocean-report.R --action=report --profile=paper --cores=12 --data-root=data/ocean_silicate/n500 --output=results/ocean500
```

For ADNI or OAI, use the same commands with the script prefix and input/output
paths from the table. Each method group produces its own report before the
combined report. Use `--action=report` on a method script to rebuild that group's
report from saved results. Combined report scripts accept `plan` or `report`.

`--cores` is a budget **per command**, capped by detected capacity. Run the two
method groups sequentially to
keep this example within 12 cores. Concurrent commands need budgets whose sum
fits the machine. Use at least four cores for paper workflows to retain four
chains; Windows executes serially. Thread counts should remain one per worker.

For a short Ocean check, use `--profile=quick`, a separate output root,
and Ocean200 with `--dataset=ocean200 --data-root=data/ocean_silicate/n200`.
A quick check exercises the workflow; it is not a convergence assessment.
Continuation is an explicit operation documented in each application guide;
reporting does not start it.

## 5. Find outputs and inspect results

Applications place group reports in `reports/eivGP/index.html` and
`reports/competitors/index.html` beneath the output root. Combined paper exports
are under `reports/paper/`; inspect `artifact_manifest.csv` for their identities.
Study I's EIV-GP workflow writes `core/study1-report/study1.html`. Study II report
locations and selections are described in its reproduction guide. The versioned
application HTML pages in `paper/` are data-free plans; empirical application
reports come from your output root.

Before reporting results, check completion/failure status, scientific-target
R-hat, bulk/tail ESS, MCSE and traces. A completed command or an ESS threshold
alone does not establish convergence. Use replication-level scores for
comparisons, preserving failed fits and diagnostic warnings.

| Location | What it contains |
| --- | --- |
| `eivGP/R/` | Model fitting, numerical helpers, prediction, diagnostics and plotting. |
| `eivGP/src/` | Native HMC target, gradient and leapfrog calculations. |
| `tools/` | Executable workflows; `APP-tools.R` holds application-specific functions. |
| `data/` | Related input collections grouped by dataset/cohort. |
| `paper/` | Report sources, readable HTML, manuscript assets and provenance manifests. |
| `results/` | Ignored local computation outputs; use external stores for large runs. |

For development checks, run `Rscript tools/check-portability.R` after installation.
Run the package tests from an R session that has loaded this checkout:

```r
testthat::test_dir("eivGP/tests/testthat", package = "eivGP",
                   stop_on_failure = TRUE)
```

Application checks are in `tests/application/`. The [paper overview](paper/README.md)
links report sources and saved-result guides.
