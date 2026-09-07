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

From a repository checkout:

```sh
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
