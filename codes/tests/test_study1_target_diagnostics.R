# Target-panel regression checks, not a substantive convergence certification.
args <- commandArgs(FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
codes <- dirname(dirname(normalizePath(script)))
for (module in c("00_parallel_utils.R", "00_study1_functions.R", "simulation_helpers.R")) {
  source(file.path(codes, module))
}
set.seed(415L)
X <- cbind(a = seq(-1, 1, length.out = 12L), b = rnorm(12L))
u <- seq(-1.8, 1.8, length.out = 12L)
C <- as.integer(cut(u, c(-Inf, -0.5, 0.5, Inf)))
y <- sin(X[, 1L]) + 0.3 * u + rnorm(12L, sd = 0.1)
for (kernel in c("se", "matern")) {
  fit <- fit_eivgp_1d(X, y, C, u_true = u, calib_idx = c(1L, 6L, 12L),
    m = 3L, n_chains = 4L, n_iter = 80L, burn = 40L, thin = 1L,
    kernel = kernel, matern_nu = 2.5, parallel_chains = FALSE, seed = 21L)
  seed_before <- .Random.seed
  panel <- c(1L, 6L, 12L)
  series <- mixedgp_study1_target_series(fit, X[panel, ], C[panel], u[panel],
    max_draws_per_chain = 20L, n_latent = 16L, seed = 4L)
  stopifnot(identical(seed_before, .Random.seed), length(series) == 21L,
    all(vapply(series, length, integer(1)) == 4L),
    all(lengths(unlist(series, recursive = FALSE)) == 20L),
    all(is.finite(unlist(series))))
  repeat_series <- mixedgp_study1_target_series(fit, X[panel, ], C[panel], u[panel],
    max_draws_per_chain = 20L, n_latent = 16L, seed = 4L)
  stopifnot(identical(series, repeat_series))
  diagnostics <- mixedgp_summarize_diagnostic_series(series)
  stopifnot(!mixedgp_diagnostic_table_pass(diagnostics))
  # First signature is independently reproducible from the empirical-Q routine.
  s <- 21L
  tau <- fit$mcmc$samples_by_chain$tau[[1L]][s, ]
  set.seed(4L)
  um <- matrix(rtruncnorm_vec(0, 1,
    rep(c(-Inf, tau)[C[panel]], each = 16L),
    rep(c(tau, Inf)[C[panel]], each = 16L)), 16L, 3L)
  xs <- sweep(sweep(X[panel, ], 2L, fit$data$x_center, "-"), 2L, fit$data$x_scale, "/")
  direct <- gp_integrated_mean_state_1d(fit$data$x,
    fit$mcmc$samples_by_chain$u[[1L]][s, ], fit$data$y, xs, um,
    fit$mcmc$samples_by_chain$logtheta[[1L]][s, ],
    fit$mcmc$samples_by_chain$sigma2[[1L]][s], kernel = kernel, matern_nu = 2.5)
  stopifnot(isTRUE(all.equal(series[["m_conditional_mean[1]"]][[1L]][1L],
    fit$data$y_center + fit$data$y_scale * direct$mean[1L], tolerance = 1e-12)))
}
message("Study I target diagnostics: SE/Matérn shapes, RNG, repeatability, formula checks passed.")
