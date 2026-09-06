# eivmixgp

`eivmixgp` fits Gaussian process regressions with numeric
predictors and ordinal proxies for latent continuous inputs. Sparse
calibration measurements of the latent inputs can be supplied without
discarding the larger sample observed only through the ordinal proxies.

## Installation

Build and install the package from the generated source directory:

```sh
R CMD build eivmixgp
R CMD INSTALL eivmixgp_0.2.1.tar.gz
```

## Statistical targets

The fitted posterior supports four related tasks:

- inference on the observed-input mean `m(x,c) = E{f(x,U) | C=c}`;
- inference on the latent response surface `f(x,u)`;
- prediction of a future response given `(x,c)` or an exact new `u`; and
- posterior imputation of latent inputs for training or prospective units.

## Stable interface

```r
fit <- fit_eivgp(
  X, y, C, U_obs,
  engine = "multivariate",
  latent_dim = 2,
  ident = "none", # unrestricted loadings; requires some calibration
  kernel = "matern",
  n_chains = 4
)

summary(fit)

m_draws <- predict(
  fit, new_X = X_new, new_C = C_new, target = "mean"
)
y_draws <- predict(
  fit, new_X = X_new, new_C = C_new, target = "response"
)
f_draws <- predict(
  fit, new_X = X_grid, new_U = U_grid, target = "surface"
)
u_draws <- impute_eivgp(fit, new_C = C_new)
```

A small self-contained univariate smoke test is below. Its one short chain
is for checking installation only, not for scientific inference:

```r
set.seed(1)
n <- 24L
u <- seq(-2, 2, length.out = n)
X <- cbind(exposure = stats::runif(n, -1, 1))
C <- ordered(
  cut(u, c(-Inf, -0.5, 0.6, Inf), labels = c('low', 'mid', 'high')),
  levels = c('low', 'mid', 'high')
)
y <- sin(X[, 1]) + 0.7 * tanh(u) + stats::rnorm(n, sd = 0.1)
U_obs <- rep(NA_real_, n)
U_obs[c(2L, 7L, 13L, 19L, 23L)] <- u[c(2L, 7L, 13L, 19L, 23L)]
fit <- fit_eivgp(
  X, y, C, U_obs, engine = 'univariate', parallel = FALSE,
  n_iter = 300L, burn = 150L, n_chains = 1L, # smoke test only
  preset = 'fast', seed = 11L, verbose = FALSE
)
new_X <- cbind(exposure = c(-0.5, 0.5))
new_C <- ordered(c('low', 'high'), levels = levels(C))
m_draws <- predict(
  fit, new_X = new_X, new_C = new_C, target = 'mean',
  n_latent = 64L, seed = 12L
)
dim(m_draws)
```

Select `engine = "univariate"` for one deterministically thresholded
latent input and `engine = "multivariate"` for the noisy ordinal-probit
measurement model. The multivariate latent dimension must be specified
when no calibration measurements are available; it is never guessed.

## Diagnose, extend, and recheck

Use at least four independently initialized chains for scientific work.
Every post-warm-up draw is kept; the public API rejects `thin > 1`.
Choose panel locations covering the scientific comparisons you will report:

```r
report <- diagnose_eivgp(fit, X = X_new, C = C_new)
report
report$table # R-hat, bulk/tail ESS, mean MCSE, and per-functional status
saveRDS(fit, 'fit-checkpoint.rds')
# After reviewing diagnostics, if additional sampling is appropriate:
fit <- continue_eivgp(readRDS('fit-checkpoint.rds'), n_iter = 2000L)
report <- diagnose_eivgp(fit, X = X_new, C = C_new)
```

`n_iter` in `continue_eivgp()` means additional transitions per chain.
It resumes terminal states and RNG streams, freezes warm-up tuning,
and appends draws within each chain. It does not restart or repeat warm-up.
Keep the same package version and numerical environment when resuming.
Old fits without checkpoints must be refitted. This is not periodic crash
recovery while a fitting call is still running.

Poor exploration requires investigation, not an automatic longer run.
Agreement with low ESS or excessive MCSE suggests extending and rechecking.
A diagnostic pass is only a screen for the monitored functionals, not proof
of convergence, identification, or model adequacy. Supply `U` with the panel
to additionally check calibrated `f(x,u)` summaries, and `additional_series`
for scientific contrasts aligned with every retained draw in every chain.
Check sensitivity to `n_latent` and the integration seed: integration error
is not included in reported MCMC error. Quantiles need separate MCSE checks.

## Calibration and latent scale

Raw or physical latent coordinates are reported only when calibrated
latent inputs have affine rank `d + 1`. This is a conservative software
reporting rule, not a necessary-condition theorem for identification.
Otherwise outputs use model coordinates. Without calibration, inference on
`m(x,c)` and response prediction given `(x,c)` remain model-defined, while
`f(x,u)` and latent-input imputations must be interpreted on the working
scale. `summary(fit)` reports the calibration count, affine rank, and active
scale; requesting `scale = "raw"` when the scale is unanchored is an error.

Both engines currently fix the latent population law to `N_d(0,I_d)`,
independent of X in model coordinates. Calibration standardization is
on by default in the public interface and uses calibrated training rows;
these plug-in centers/scales have no posterior uncertainty, and scaling
does not remove correlation. Inspect `fit$model_specification`.

With sparse/rank-deficient calibration, select the loading structure
explicitly: `ident = "none"` leaves loadings unrestricted, whereas
`ident = "lower_triangular"` imposes structural zeros and positive
diagonal entries. Triangular loadings with diagonal ARD are not merely
a coordinate convention. A richer `U | X` population model is not yet
implemented; assess these assumptions before physical-scale inference.

## Prospective latent simulation

For the multivariate engine, a minimax-tilted accept--reject sampler for
`U | C` is the default (`latent_sampler = "minimax_tilting"`) and is exact
on successful completion. Exact prior rejection
is available as `"rejection"`; finite Gibbs sampling is retained only as
an explicitly requested diagnostic approximation. The default uses
`TruncatedNormal::rtmvnorm`.

## Output contracts

`predict()` and `predict_eivgp()` return a numeric matrix whose columns
follow the evaluation rows and whose rows follow posterior `draw_ids`, with
the `n_per_draw` replicates for each posterior state contiguous. The `joint`
attribute distinguishes coherent multivariate draws from pointwise draws;
use `joint = TRUE` for simultaneous functionals. `impute_eivgp()` returns a
posterior-draw by unit by latent-coordinate array and labels the source,
coordinate scale, and calibration anchoring in attributes.

## Reproducible numerical experiments

Both simulation studies use the same validated runner:

```r
study1 <- run_mixedgp_experiment(
  "study1", config = "quick", parallel_level = "chains"
)
study2 <- run_mixedgp_experiment(
  "study2", config = "quick", parallel_level = "chains"
)
```

Use `config = "balanced"` or `"thorough"` for publication-scale runs.
Frozen synthetic datasets are collision-checked, and seeded serial and fork
backends preserve the caller's random-number stream. Previously archived
Study I fits are reused only through an explicit `study_options` opt-in.
Published competitor packages and reporting packages are optional
dependencies; the experiment preflight records availability and versions.

## Citation

Run `citation('eivmixgp')` for the manuscript citation bundled with the
package.
