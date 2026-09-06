## Estimator-specific MCSE and complete-window regression tests.
args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", args[1L]), mustWork = TRUE)
code_dir <- dirname(dirname(script))
for (module in c("00_parallel_utils.R", "00_study2_functions.R", "simulation_helpers.R")) {
  source(file.path(code_dir, module))
}

set.seed(606L)
skewed <- list(skewed = lapply(1:4, function(i) {
  exp(as.numeric(stats::arima.sim(list(ar = 0.9), n = 4000L)) / 2)
}))
tab <- mixedgp_summarize_diagnostic_series(skewed)
mat <- do.call(cbind, skewed$skewed)
stopifnot(isTRUE(all.equal(tab$mcse, as.numeric(posterior::mcse_mean(mat)))),
  isTRUE(all.equal(tab$ess_mean, as.numeric(posterior::ess_mean(mat)))),
  abs(tab$mcse - tab$posterior_sd / sqrt(tab$ess_bulk)) > 1e-5)
iid <- list(mean = replicate(4L, rnorm(2000), simplify = FALSE))
good <- mixedgp_summarize_diagnostic_series(iid)
stopifnot(mixedgp_diagnostic_table_pass(good),
  mixedgp_diagnostic_table_pass(good, mcse_limit = good$mcse * 2),
  !mixedgp_diagnostic_table_pass(good, mcse_limit = good$mcse / 2),
  !mixedgp_diagnostic_table_pass(good, mcse_sd_ratio_limit = good$mcse_sd_ratio / 2))
iid <- mixedgp_mark_diagnostic_window(iid, FALSE)
short <- mixedgp_summarize_diagnostic_series(iid)
stopifnot(!short$target_pass, !mixedgp_diagnostic_table_pass(short),
  !mixedgp_summarize_diagnostic_series(c(iid))$target_pass,
  !mixedgp_summarize_diagnostic_series(iid["mean"])$target_pass,
  identical(mixedgp_diagnostic_draw_indices(1600L), seq_len(1600L)),
  identical(mixedgp_diagnostic_draw_indices(1600L, 1000L), 601:1600))
for (value in list(NA_real_, numeric(0), -Inf, 0, 1.5, c(1, 2))) {
  stopifnot(inherits(try(mixedgp_diagnostic_draw_indices(10L, value), silent = TRUE), "try-error"))
}

## Verify the Study II integrated-mean variance against an independent full
## posterior covariance average, and prospective moments on supplied U units.
n_saved <- 5L
X <- matrix(c(-1, 0, 1), 3L, 1L)
U <- matrix(c(0, -.5, .7, 0, .4, -.3), 3L, 2L)
A <- matrix(c(1, .4, 0, .9), 2L, 2L)
tau <- c(-.6, .6, -.6, .6)
tiny <- list(data = list(p = 1L, q = 2L, d = 2L, m_vec = c(3L, 3L),
  X = X, y = c(-.5, .1, .4), X_center = 0, X_scale = 1,
  y_center = 2, y_scale = 3, U_center = c(10, 20), U_scale = c(2, 4),
  calib_idx = 1L, miss_idx = 2:3, ident = "lower_triangular"),
  kernel = list(name = "se", matern_nu = 2.5), mcmc = list(samples_by_chain = list(
    A = replicate(4L, array(rep(as.numeric(A), each = n_saved), c(n_saved, 2L, 2L)), simplify = FALSE),
    U = replicate(4L, array(rep(as.numeric(U), each = n_saved), c(n_saved, 3L, 2L)), simplify = FALSE),
    tau = replicate(4L, matrix(rep(tau, each = n_saved), n_saved, 4L), simplify = FALSE),
    logtheta = replicate(4L, matrix(rep(log(c(2, .5, .5, .5)), each = n_saved), n_saved, 4L), simplify = FALSE),
    sigma2 = replicate(4L, rep(.1, n_saved), simplify = FALSE))))
xp <- matrix(.25, 1L, 1L)
cp <- matrix(c(2L, 3L), 1L, 2L)
for (kernel in c("se", "matern")) {
  tiny$kernel$name <- kernel
  set.seed(81L)
  before <- .Random.seed
  series <- mixedgp_study2_target_series(tiny, xp, cp, n_latent = 8L, seed = 91L)
  stopifnot(identical(before, .Random.seed), isTRUE(attr(series, "draw_window_complete")),
    all(lengths(unlist(series, recursive = FALSE)) == n_saved))
  set.seed(91L)
  um <- sample_u_given_c_ordprobit_minimax(cp[rep(1L, 8L), , drop = FALSE],
    A, list(tau[1:2], tau[3:4]), c(3L, 3L))
  pred <- gp_predict_draw_general(X, U, tiny$data$y,
    xp[rep(1L, 8L), , drop = FALSE], um, log(c(2, .5, .5, .5)), .1,
    noisy = FALSE, kernel = kernel, matern_nu = 2.5, return_cov = TRUE)
  stopifnot(isTRUE(all.equal(series[["m_conditional_variance[1]"]][[1L]][1L],
    9 * mean(pred$cov), tolerance = 1e-10)),
    isTRUE(all.equal(series[["prospective_U_mean[1,2]"]][[1L]][1L],
      mean(20 + 4 * um[, 2L]), tolerance = 1e-12)))
  truncated <- mixedgp_study2_target_series(tiny, xp, cp,
    max_draws_per_chain = 4L, n_latent = 8L, seed = 91L)
  stopifnot(identical(attr(truncated, "draw_window_complete"), FALSE))
  stopifnot(!any(mixedgp_summarize_diagnostic_series(c(truncated))$draw_window_complete))
}
cat("Diagnostic precision, reporting window, integrated variance, and scale tests passed.\n")
