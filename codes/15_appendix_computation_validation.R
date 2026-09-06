############################################################
## 15_appendix_computation_validation.R
##
## Focused release gates for the appendix GP baselines, the response-free
## measurement-model diagnostics, predictive-variance safeguards, and exact
## rejection-sampler telemetry.  This script is intentionally small and does
## not run a publication Monte Carlo replication.
############################################################

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
codes_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
} else {
  getwd()
}

source(file.path(codes_dir, "00_parallel_utils.R"))
source(file.path(codes_dir, "00_study2_functions.R"))
source(file.path(codes_dir, "03_study2_published_competitors.R"))
source(file.path(codes_dir, "04_study2_ablations.R"))

if (!requireNamespace("posterior", quietly = TRUE)) {
  stop("The appendix computation release gate requires package 'posterior'.")
}

## Multistart GP MLE must preserve caller RNG state and persist every start.
set.seed(191L)
X_gp <- cbind(
  x1 = seq(-1, 1, length.out = 20L),
  x2 = sin(seq(-1, 1, length.out = 20L))
)
y_gp <- sin(2 * X_gp[, 1L]) + rnorm(nrow(X_gp), 0, 0.05)
set.seed(777L)
rng_before <- .Random.seed
fit_gp <- gp_mle_fit(
  X_gp, y_gp, n_starts = 3L, seed = 92L, maxit = 200L
)
required_attempt_columns <- c(
  "start", "n_starts", "optimizer_seed", "maxit",
  "initial_log_sigma2", "initial_log_rho", "initial_log_theta1",
  "final_log_sigma2", "final_log_rho", "final_log_theta1",
  "convergence", "objective", "message"
)
stopifnot(
  fit_gp$convergence == 0L,
  nrow(fit_gp$optimizer_attempts) == 3L,
  fit_gp$selected_start %in% seq_len(3L),
  all(required_attempt_columns %in% names(fit_gp$optimizer_attempts)),
  identical(rng_before, .Random.seed)
)

## A deliberately zero-iteration optimizer budget must fail closed and retain
## the attempt table on the error object.
failed_gp <- tryCatch(
  gp_mle_fit(X_gp, y_gp, n_starts = 2L, seed = 92L, maxit = 0L),
  error = function(e) e
)
stopifnot(
  inherits(failed_gp, "mixedgp_gp_mle_error"),
  nrow(failed_gp$optimizer_attempts) == 2L,
  all(failed_gp$optimizer_attempts$convergence != 0L)
)

## Variance guards clamp only negative values within a roundoff tolerance.
variance_near_zero <- predictive_variance_diagonal(
  c(-.Machine$double.eps, 1), 2L, "validation adapter"
)
variance_failure <- tryCatch(
  predictive_variance_diagonal(c(-1e-4, 1), 2L, "validation adapter"),
  error = function(e) e
)
prepared <- list(
  X = matrix(0, 1L, 1L), alpha = 0,
  chol_train = matrix(sqrt(1 / (1 + 1e-15)), 1L, 1L),
  sigma2 = 1, rho = 1, theta = 1
)
appendix_near_zero <- gp_mle_predict_prepared(prepared, matrix(0, 1L, 1L))
prepared$chol_train <- matrix(0.5, 1L, 1L)
appendix_variance_failure <- tryCatch(
  gp_mle_predict_prepared(prepared, matrix(0, 1L, 1L)),
  error = function(e) e
)
stopifnot(
  variance_near_zero[[1L]] == 0,
  inherits(variance_failure, "error"),
  appendix_near_zero$var[[1L]] == 0,
  inherits(appendix_variance_failure, "error")
)

## Exact rejection records pattern-specific proposal probabilities and returns
## the same telemetry on a structured underfill error.
A_rejection <- matrix(c(1, 0.5), 2L, 1L)
tau_rejection <- list(0, 0)
C_rejection <- rbind(
  matrix(c(1L, 1L), 3L, 2L, byrow = TRUE),
  matrix(c(2L, 2L), 4L, 2L, byrow = TRUE)
)
set.seed(88L)
rejection_draw <- sample_u_given_c_ordprobit_rejection(
  C_rejection,
  A_rejection,
  tau_rejection,
  c(2L, 2L),
  batch_size = 1000L,
  max_batches = 2L
)
rejection_telemetry <- attr(rejection_draw, "pattern_telemetry")
stopifnot(
  nrow(rejection_telemetry) == 2L,
  all(rejection_telemetry$retained_draws ==
        rejection_telemetry$requested_draws),
  all(rejection_telemetry$proposal_hits >=
        rejection_telemetry$retained_draws),
  all(rejection_telemetry$complete)
)
rejection_failure <- tryCatch(
  sample_u_given_c_ordprobit_rejection(
    matrix(c(1L, 1L), 2L, 2L, byrow = TRUE),
    A_rejection,
    tau_rejection,
    c(2L, 2L),
    batch_size = 1L,
    max_batches = 1L
  ),
  error = function(e) e
)
stopifnot(
  inherits(rejection_failure, "mixedgp_rejection_sampling_error"),
  any(rejection_failure$pattern_telemetry$shortfall > 0L),
  rejection_failure$candidate_draws == 1
)

## Chain-aware rank-normalized R-hat and bulk/tail ESS must be finite for two
## chains.  A one-chain fit is never permitted to pass the gate.
set.seed(501L)
n_measurement <- 36L
U_measurement <- matrix(rnorm(n_measurement), n_measurement, 1L)
A_measurement <- matrix(c(1.4, 0.8), 2L, 1L)
S_measurement <- U_measurement %*% t(A_measurement) +
  matrix(rnorm(2L * n_measurement), n_measurement, 2L)
C_measurement <- apply(S_measurement, 2L, function(z) {
  findInterval(z, c(-Inf, -0.5, 0.5, Inf))
})
C_measurement <- matrix(as.integer(C_measurement), n_measurement, 2L)

measurement_fit <- fit_ordinalprobit_measurement_fb(
  C_ord = C_measurement,
  U_obs = U_measurement,
  calib_idx = seq_len(8L),
  d = 1L,
  m_vec = c(3L, 3L),
  n_iter = 100L,
  burn = 40L,
  thin = 1L,
  n_chains = 2L,
  seed = 73L,
  parallel_chains = FALSE,
  rhat_limit = 1e6,
  bulk_ess_limit = 1,
  tail_ess_limit = 1
)
key_diagnostics <- measurement_fit$diagnostic_parameters$key
latent_diagnostics <- measurement_fit$diagnostic_parameters$latent_rhat
stopifnot(
  identical(
    measurement_fit$diagnostics$diagnostic_backend,
    "posterior_rank_normalized"
  ),
  isTRUE(measurement_fit$diagnostics$convergence_pass),
  all(is.finite(key_diagnostics$rhat)),
  all(is.finite(key_diagnostics$ess_bulk)),
  all(is.finite(key_diagnostics$ess_tail)),
  all(is.finite(latent_diagnostics$rhat)),
  all(is.finite(latent_diagnostics$ess_bulk)),
  all(is.finite(latent_diagnostics$ess_tail))
)

## PI-GP and CC-GP must use the same multistart engine and retain the exact
## rejection telemetry when that independent validation sampler is requested.
X_measurement <- cbind(
  x1 = seq(-1, 1, length.out = n_measurement),
  x2 = cos(seq(-1, 1, length.out = n_measurement))
)
y_measurement <- sin(X_measurement[, 1L]) +
  0.4 * U_measurement[, 1L] + rnorm(n_measurement, 0, 0.05)
pi_fit <- fit_study2_pi_gp(
  X_measurement,
  y_measurement,
  measurement_fit,
  gp_n_starts = 2L,
  gp_seed = 300L,
  gp_maxit = 200L
)
cc_fit <- fit_study2_cc_gp(
  X_measurement,
  y_measurement,
  U_measurement,
  calib_idx = seq_len(8L),
  measurement_fit = measurement_fit,
  gp_n_starts = 2L,
  gp_seed = 301L,
  gp_maxit = 200L
)
pi_draw <- sample_study2_pi_gp(
  pi_fit,
  X_measurement[seq_len(5L), , drop = FALSE],
  C_measurement[seq_len(5L), , drop = FALSE],
  n_draw = 6L,
  n_measurement_draw = 3L,
  seed = 302L,
  latent_sampler = "rejection",
  rejection_batch_size = 2000L,
  rejection_max_batches = 4L
)
cc_draw <- sample_study2_cc_gp(
  cc_fit,
  X_measurement[seq_len(5L), , drop = FALSE],
  C_measurement[seq_len(5L), , drop = FALSE],
  n_draw = 6L,
  seed = 303L,
  latent_sampler = "rejection",
  rejection_batch_size = 2000L,
  rejection_max_batches = 4L
)
stopifnot(
  nrow(pi_fit$gp$optimizer_attempts) == 2L,
  nrow(cc_fit$gp$optimizer_attempts) == 2L,
  sum(pi_fit$gp$optimizer_attempts$convergence == 0L) >= 1L,
  sum(cc_fit$gp$optimizer_attempts$convergence == 0L) >= 1L,
  identical(dim(pi_draw), c(6L, 5L)),
  identical(dim(cc_draw), c(6L, 5L)),
  is.data.frame(attr(pi_draw, "pattern_telemetry")),
  is.data.frame(attr(cc_draw, "pattern_telemetry"))
)

one_chain_fit <- fit_ordinalprobit_measurement_fb(
  C_ord = C_measurement,
  U_obs = U_measurement,
  calib_idx = seq_len(8L),
  d = 1L,
  m_vec = c(3L, 3L),
  n_iter = 50L,
  burn = 20L,
  thin = 1L,
  n_chains = 1L,
  seed = 73L,
  parallel_chains = FALSE,
  rhat_limit = 1e6,
  bulk_ess_limit = 1,
  tail_ess_limit = 1
)
stopifnot(!isTRUE(one_chain_fit$diagnostics$convergence_pass))

## Exercise the identical guard in the ocean real-data prediction helper.
source(file.path(codes_dir, "00_study1_functions.R"))
source(file.path(codes_dir, "real-data", "ocean_data_helpers.R"))
ocean_prepared <- list(
  X = matrix(0, 1L, 1L, dimnames = list(NULL, "x1")),
  alpha = 0,
  chol_train = matrix(sqrt(1 / (1 + 1e-15)), 1L, 1L),
  sigma2 = 1,
  rho = 1,
  theta = 1,
  kernel = list(name = "se", matern_nu = 2.5)
)
Xstar_ocean <- matrix(0, 1L, 1L, dimnames = list(NULL, "x1"))
ocean_near_zero <- gp_mle_predict_prepared_1d(ocean_prepared, Xstar_ocean)
ocean_prepared$chol_train <- matrix(0.5, 1L, 1L)
ocean_variance_failure <- tryCatch(
  gp_mle_predict_prepared_1d(ocean_prepared, Xstar_ocean),
  error = function(e) e
)
stopifnot(
  ocean_near_zero$var[[1L]] == 0,
  inherits(ocean_variance_failure, "error")
)

cat(
  "Appendix computation validation passed: multistart/fail-closed GP MLE; ",
  "rank-normalized two-chain MCMC gate; tolerance-guarded variances; ",
  "pattern-wise rejection telemetry.\n",
  sep = ""
)
