############################################################
## Regression checks for the 2026-08-30 collaborator merge.
############################################################

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = TRUE)
} else {
  normalizePath(sys.frame(1)$ofile, mustWork = TRUE)
}
code_dir <- dirname(script_file)
source(file.path(code_dir, "00_study1_functions.R"))
OCEAN_REALDATA_DIR <- file.path(code_dir, "real-data")
source(file.path(OCEAN_REALDATA_DIR, "ocean_data_helpers.R"))

set.seed(20260830)
n <- 24L
X <- cbind(
  temperature = rnorm(n),
  salinity = runif(n),
  depth = rnorm(n)
)
u <- rnorm(n)
tau <- c(-0.45, 0.35)
c_ord <- make_class(u, tau)
while (any(tabulate(c_ord, nbins = 3L) < 3L)) {
  u <- rnorm(n)
  c_ord <- make_class(u, tau)
}
y <- sin(X[, 1L]) + 0.3 * X[, 2L] - 0.2 * X[, 3L] + cos(u) + rnorm(n, 0, 0.1)

calib_idx <- c(1L, 5L, 9L, 13L, 17L, 21L)
u_obs <- rep(NA_real_, n)
u_obs[calib_idx] <- u[calib_idx]

fit <- fit_eivgp_1d(
  x_raw = X,
  y_raw = y,
  c_ord = c_ord,
  calib_idx = calib_idx,
  m = 3L,
  tau_true = tau,
  n_iter = 30L,
  burn = 10L,
  thin = 1L,
  n_chains = 2L,
  preset = "fast",
  seed = 81L,
  kernel = "matern",
  matern_nu = 1.5,
  u_obs = u_obs
)

stopifnot(
  fit$data$p_x == 3L,
  fit$kernel$name == "matern",
  ncol(fit$mcmc$samples_logtheta) == 5L,
  all(c("rhat", "ess_bulk", "ess_tail", "ess") %in%
        names(fit$diagnostics$ess_key))
)

## A named p-vector denotes one prediction row; reversed data-frame columns
## are restored to the training order.
x_one <- c(temperature = 0, salinity = 0.5, depth = -0.2)
y_one <- sample_eiv_test_y(x_one, 2L, fit, draw_ids = 1:2)
stopifnot(identical(dim(y_one), c(2L, 1L)), all(is.finite(y_one)))

x_two <- data.frame(
  depth = c(-0.2, 0.1),
  salinity = c(0.5, 0.6),
  temperature = c(0, 0.2)
)
y_two <- sample_eiv_test_y(x_two, c(2L, 3L), fit, draw_ids = 1:2)
stopifnot(identical(dim(y_two), c(2L, 2L)), all(is.finite(y_two)))

f_draws <- sample_eiv_f_given_xu_1d(
  x_star_raw = x_two,
  u_star = c(0, 0.7),
  fit_obj = fit,
  draw_ids = 1:2
)
stopifnot(identical(dim(f_draws), c(2L, 2L)), all(is.finite(f_draws)))

baselines <- fit_embedding_baselines(
  X, y, c_ord, m = 3L, n_starts_learned = 2L,
  kernel = "matern", matern_nu = 1.5
)
baseline_draws <- predict_embedding_baseline_samples(
  baselines, x_one, 2L, m = 3L, n_draw = 3L
)
stopifnot(all(vapply(
  baseline_draws,
  function(z) identical(dim(z), c(3L, 1L)) && all(is.finite(z)),
  logical(1)
)))

cc_fit <- fit_ocean_latent_gp(
  X, y, u, train_idx = calib_idx,
  kernel = "matern", matern_nu = 1.5
)
cc_draws <- sample_ocean_latent_gp_y_given_xc(
  cc_fit, X[1:4, , drop = FALSE], c_ord[1:4], tau,
  n_draw = 4L, seed = 91L
)
stopifnot(identical(dim(cc_draws), c(4L, 4L)), all(is.finite(cc_draws)))

allocation <- data.frame(
  c = 1:3,
  n_calib_3 = rep(1L, 3L),
  n_calib_6 = rep(2L, 3L)
)
nested <- make_ocean_nested_calibration(c_ord, allocation, c(3L, 6L), seed = 19L)
stopifnot(
  length(nested[["3"]]) == 3L,
  length(nested[["6"]]) == 6L,
  all(nested[["3"]] %in% nested[["6"]])
)

## Preserve the original p_x = 1 interface used by Study I scripts.
dat_1d <- simulate_1d_data(n = 24L, n_test = 5L, m = 6L, seed = 32L)
fit_1d <- fit_eivgp_1d(
  dat_1d$train$x,
  dat_1d$train$y,
  dat_1d$train$c,
  u_true = dat_1d$train$u,
  calib_idx = seq_len(6L),
  m = 6L,
  tau_true = dat_1d$tau_true,
  n_iter = 24L,
  burn = 8L,
  thin = 2L,
  n_chains = 2L,
  preset = "fast",
  seed = 33L
)
stopifnot(ncol(fit_1d$mcmc$samples_logtheta) == 3L)

cat("Study I collaborator merge validation: OK\n")
