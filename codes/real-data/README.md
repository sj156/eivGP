# Ocean application workflow

This directory contains the audited merge of the collaborator's ocean scripts
uploaded on 2026-08-30. The shared sampler is
`../00_study1_functions.R`.

## Data

Set `OCEAN_DATA_DIR` to a directory containing `study1_data.csv` and
`meta.json`. The CSV must contain the train/test split, response `y`, ordinal
code `c`, standardized latent input `u_log_std`, and every predictor named by
`meta.json`. The prepared data themselves are not stored in this repository.

The class-proportional calibration design is stored in
`data/class6_proportional_calibration_allocation.csv`. The loader verifies that
its class counts match the prepared training data before fitting.

## Run

From any working directory:

```sh
Rscript revision/codes/real-data/run_ocean_all.R
```

Edit `OCEAN_QUICK`, `OCEAN_KERNEL`, and `OCEAN_MATERN_NU` near the top of the
runner. Quick and paper runs, covariance families, and cache schema versions
use separate cache filenames.

## Comparison target

All primary predictive methods are scored for `y` given `x` and ordinal `c`.
The complete-case GP uses only calibrated training rows and integrates a new
latent input over its ordinal interval. The full-`u` oracle uses all training
latent inputs but also integrates the test latent input over its interval; it
never conditions on held-out test `u`.

Run `../07_study1_collaborator_merge_validation.R` for a synthetic regression
test of multivariate predictors, sparse calibration values, SE/Matérn kernels,
same-target prediction, and posterior inference for `f(x,u)`.
