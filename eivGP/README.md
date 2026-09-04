# eivGP

`eivGP` implements errors-in-variables Gaussian process regression
with ordinal proxies for latent continuous inputs and sparse calibration.

## Stable model interface

```r
fit <- fit_eivgp(X, y, C, U_obs, engine = 'multivariate')
mean_draws <- predict(fit, new_X = X_new, new_C = C_new, target = 'mean')
response_draws <- predict(fit, new_X = X_new, new_C = C_new, target = 'response')
latent_draws <- impute_eivgp(fit, new_C = C_new)
```

## Package organization

Model code is organized by responsibility rather than by paper study:
`model_api.R`, `model_univariate.R`, `model_multivariate.R`,
`core_numerics.R`, `core_parallel.R`, and `application_interface.R`.
Study-specific names are confined to frozen reproduction workflows.

## Publication workflows

Both exported runners default to a read-only dry run:

```r
run_study1_simulation()
run_study2_simulation()
```

Use `study1_simulation_config()` or `study2_simulation_config()` to
select data generation, smoke testing, or the frozen publication design.
Smoke and publication fits automatically run the exact competitor and
paired-design validators before computation; publication requires every
competitor validation to succeed. Results are archived under `config/`.
See the repository's standalone experiment R Markdown files for the full
top-to-bottom reproduction workflow.
