# Study II reproducibility guide

This folder contains the frozen multivariate-ordinal experiment described in
the manuscript. Run the scripts from this directory.

## Files

- `00_study2_functions.R`: EIV-GP sampler, prediction, scoring, and the frozen
  Study II data-generating mechanisms.
- `01_study2_representative_figures.R`: representative-data figures and the
  coordinate-wise imputation table.
- `02_study2_monte_carlo.R`: resumable repeated experiment, manuscript tables,
  pattern-frequency analysis, and MCMC publication gate.
- `03_study2_published_competitors.R`: common Study I/II wrappers for the
  public `kergp`, `LVGP`, and `EzGP` implementations, plus preflight and
  availability/failure records. The historical filename is retained for
  compatibility.
- `04_study2_ablations.R`: response-free PI-GP, CC-GP, and Full-U GP appendix
  ablations. These are not literature competitors.
- `05_study2_sampler_validation.R`: small exact-sampler audit.
- `13_study2_prospective_latent_validation.R`: focused distributional and
  failure-mode checks for prospective latent draws used by EIV-GP, PI-GP, and
  CC-GP.
- `06_realdata_interface_validation.R`: end-to-end checks for multivariate
  quantitative inputs, ordinal proxies, latent variables, and both kernels.
- `08_published_competitor_validation.R`: interface checks and small real-fit
  tests for every installed public competitor package.
- `09_experiment_design_validation.R`: fast checks that the Study I negative
  control and the three Gaussian Study II scenarios retain their intended
  paired random quantities.
- `real-data.R`: application-facing fitting and prediction wrappers.
- `run_study2_all.R`: user-facing master configuration.

The design identifier is
`study2-manuscript-v10-tailstable_q4-d2_A-fixed_balanced-4-level_sigma0.12_exact-minimax-interwoven_common-random-numbers`.
Caches with another identifier must not be reused.

For a fixed replication, the primary, latent-additive, and high-uncertainty
Gaussian-score scenarios share `X`, `U`, Gaussian score innovations, and
response innovations. The runner reports paired changes in each
competitor-minus-EIV CRPS and interval-score gap. The logistic robustness
scenario shares `X` and `U` but uses a separate logistic score-error
transformation.

## Real-data interface

`fit_eivgp_real_data()` accepts numeric matrix/data-frame `X`, a univariate
numeric response, multiple ordinal columns in `C`, and a multivariate `U_obs`.
Rows of `U_obs` must be either completely observed or completely missing.
Ordered factors retain their declared level order; character or unordered-factor
proxies require an explicit `ordinal_levels` list.
When calibrated latent covariates anchor substantively named coordinates, the
real-data wrapper uses unconstrained loadings (`ident = "none"`). It uses the
lower-triangular convention only when the latent coordinates are completely
unobserved; users can override this choice explicitly.

Choose `kernel = "se"` for the squared-exponential covariance or
`kernel = "matern"` together with `matern_nu` (2.5 by default). The selected
kernel is used consistently in MCMC and prediction. Use
`predict_eivgp_y_given_xc()` for the primary target `Y* | X*, C*`,
`predict_eivgp_f_given_xu()` for the latent surface, and
`posterior_u_draws()` to recover latent-input draws on either the model or raw
scale. Run `Rscript 06_realdata_interface_validation.R` after changing the
application interface.

## Required workflow

1. Install and freeze all R dependencies before running. The scripts never
   install packages silently. Record the lockfile or library snapshot with the
   archive. The three public comparators can be installed with
   `install.packages(c("kergp", "LVGP", "EzGP"))`.
2. Run `run_study2_all.R` with `STUDY2_CONFIG <- "quick"` after any code
   change. This is a workflow check, not a numerical result.
3. Audit the competitor status CSV. An unavailable method remains unavailable;
   it is never replaced by one of the legacy embeddings in
   `00_study2_functions.R`.
   Run `Rscript 08_published_competitor_validation.R`; for an archival run set
   `MIXEDGP_REQUIRE_ALL_COMPETITORS=true` so any unavailable or failed package
   stops the validation.
4. Use a multi-chain development run to tune the sampler. Do not launch the
   final run until the prespecified MCMC gate is plausible at every calibration
   size, especially 0 and 10.
   Also run `Rscript 13_study2_prospective_latent_validation.R`. Prospective
   `U | C` draws use a minimax-tilted accept--reject sampler for the truncated
   Gaussian score law. On successful completion, the following exact Gaussian
   `U | S` update yields exact draws. Exact prior rejection
   remains a validation option in EIV-GP and both appendix ablations. The
   optional finite-Gibbs path is diagnostic only and must be requested with
   `STUDY2_PREDICTIVE_LATENT_SAMPLER <- "gibbs"`; its sweep count is recorded in
   the design manifest.
5. For the final run, use `STUDY2_CONFIG <- "thorough"`, start with
   `STUDY2_USE_CACHE <- FALSE`, and set
   `STUDY2_ENFORCE_MCMC_GATE <- TRUE`. Resume only from caches carrying the
   exact final cache tag.
6. Archive the raw RDS bundle, all status/diagnostic CSVs, the design manifest,
   session information, package lockfile, and generated tables and figures.

## Published competitor wrappers

- UC-GP uses `kergp::q1Symm` for each ordinal factor and eight likelihood
  starts in the publication run. Because `kergp` can return a nonconverged
  start with the largest objective, the wrapper selects the largest likelihood
  among starts whose package convergence flag is successful.
- LVGP calls the authors' `LVGP::LVGP_fit(..., noise = TRUE)` and
  `LVGP::LVGP_predict()` functions with the published two-dimensional latent
  representation and eight optimization starts. It first uses the package
  default iteration limits. A package-selected convergence code 1 triggers a
  deterministic rerun of the same starts and seed with limits 300/100; one
  further prespecified-seed rescue is allowed. Each complete package call has
  a prespecified elapsed-time cap (15 minutes in Study I and 30 minutes in
  Study II under the current publication configuration). If none converges,
  or every attempt times out or errors, the fit is reported as failed. The
  status file distinguishes default convergence, rescued convergence, timeout,
  and other failure.
- EzGP calls the full authors' `EzGP::EzGP_fit()` and `EzGP::EzGP_predict()`
  functions. Its nugget is selected by three-fold training-only
  cross-validation over the frozen variance-fraction grid in the manuscript.

Each wrapper converts the package output to predictive draws for the noisy
response `Y* | X*, C*`. It records the package version, elapsed time, fitted or
tuned noise variance, and failure message. The four Study II ordinal inputs
remain four factors; they are never collapsed into a 256-level joint factor.
`run_study2_all.R` performs a package preflight. Thorough mode stops before any
simulation work if a named package is missing.

LVGP 2.1.6's internal `parallel = TRUE` path was also tested. In the current
R environment its workers fail because the package does not make its internal
`to_latent` function available to them. We therefore retain the unmodified
serial author routine (`parallel = FALSE`) and parallelize independent Monte
Carlo replications at the job level. We do not patch the package or substitute
a home-built LVGP optimizer.

## Reusable API, frozen data, and parallel runs

`load_mixedgp.R` loads the canonical computational modules. New analysis code
should use one stable call per method: `fit_eivgp()`, `fit_ucgp()`,
`fit_lvgp()`, or `fit_ezgp()`. The generic `fit_mixedgp_competitor()` is
available for loops, while the method-specific functions keep the exact
package arguments visible.

Both numerical studies now use the same launcher. For example:

```r
study1 <- run_mixedgp_experiment(
  "study1", config = "quick", output_dir = ".."
)
study2 <- run_mixedgp_experiment(
  "study2", config = "quick", output_dir = ".."
)
```

The common arguments control stages, competitors, caching, ablations, worker
count, and chain- versus replication-level parallelism. Genuine statistical
differences remain behind the dispatcher. Less common controls are passed in a
named `study_options` list using the corresponding `STUDY1_` or `STUDY2_`
prefix. Use `dry_run = TRUE` to validate and print a specification without
generating data or fitting models.

`12_generate_synthetic_datasets.R` generates the publication train/test data
and calibration subsets before any methods are fitted. The resulting RDS files
live under `../data-synthetic/`; their manifests contain seeds, schema version,
file sizes, and MD5 checksums. Publication drivers load these artifacts and
stop if they are missing.

EIV-GP chains use deterministic `parallel::mclapply()` through the shared
parallel utility. Dataset and task grids use deterministic
`parallel::mcmapply()`. On Windows both calls fall back to serial execution.
The Monte Carlo drivers permit either chain-level or replication-level
parallelism, but not nested fork pools.

## Statistical roles

- Main fitted methods: EIV-GP, UC-GP, LVGP, and EzGP.
- Reference distribution: the data-generating oracle for `Y* | X*, C*`.
- Appendix ablations: PI-GP and CC-GP.
- Infeasible fitted surface benchmark: Full-U GP.
- Latent-imputation comparator: the response-free ordinal-probit measurement
  model, labeled `Ordinal model (no Y)`. It is scored alongside EIV-GP for
  missing training states and for prospective `U* | C*`, but it is not
  presented as a mixed-input regression competitor. It is run in the
  prespecified ablation scenarios (`primary` by default).

Raw-coordinate imputation error is reported only when the calibration design
has affine rank `d + 1`. In particular, calibration size zero receives an
explicit `not_identified` row in the imputation-status file rather than a
coordinate-wise RMSE. Training and prospective latent targets are stored in
separate `target` rows, because the former can use the unit's observed response
under EIV-GP whereas the latter cannot.

The main predictive table never pools methods into a composite “best” row. It
shows each successful prespecified method separately and reports the number and
reason for failed or unavailable fits. Scenario contrasts and sparse-pattern
tables additionally report paired competitor-minus-EIV score advantages;
positive values favor EIV-GP. If a method fails, its score summary is labeled
as conditional on successful fits and paired contrasts retain only joint
successes; failures are never assigned an artificial penalty score.
