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

## Reproduce the numerical studies

The two standalone documents are:

- [`experiments/study1_numerical_experiment.Rmd`](experiments/study1_numerical_experiment.Rmd)
- [`experiments/study2_numerical_experiment.Rmd`](experiments/study2_numerical_experiment.Rmd)

Both default to a read-only dry run. Change their single `run_mode` parameter
to `publication` and render from top to bottom to execute the frozen paper
design. The runners store source hashes, package versions, frozen-data
manifests, diagnostics, paired comparisons, result checksums, and session
information in fingerprinted resumable directories. Before every smoke or
publication fit they automatically execute the exact published-competitor and
paired-design validators. Publication runs require all competitor validations
to succeed, and validator outputs are archived under the run's `config/`
directory.

The same workflows are callable directly:

```r
eivGP::run_study1_simulation(
  eivGP::study1_simulation_config(mode = "dry_run")
)
eivGP::run_study2_simulation(
  eivGP::study2_simulation_config(mode = "dry_run")
)
```

The developing real-data workflow has its own document at
[`applications/real_data_application.Rmd`](applications/real_data_application.Rmd).
It renders safely without data until the audited application bundle is ready.

## Source organization

The cleaned R scripts under `script/R-scripts/` are the computational source
of truth. The explicitly ordered `.Rmd` chapters form the `litr` package
skeleton and copy those scripts into the generated `eivGP/` package. Historical
chapters are retained under `legacy/` for provenance and are excluded from the
active package build.
