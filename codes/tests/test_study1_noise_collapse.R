## Run from any directory:
## Rscript /path/to/codes/tests/test_study1_noise_collapse.R
## Bounded target/order checks, not scientific convergence certification.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(script_arg) == 1L)
codes_dir <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
source(file.path(codes_dir, "00_study1_functions.R"))
expect_error <- function(expr, pattern) {
  message <- tryCatch(force(expr), error = conditionMessage)
  stopifnot(is.character(message), grepl(pattern, message, fixed = TRUE))
}
batch_mcse <- function(x, batch = 100L) {
  blocks <- length(x) %/% batch
  sd(colMeans(matrix(x[seq_len(blocks * batch)], batch))) / sqrt(blocks)
}
short_fit_warnings <- character(0)
short_fit <- function(expr) withCallingHandlers(force(expr), warning = function(w) {
  ## These intentionally short chains exercise storage/order, not convergence.
  ## Record the known short-chain ESS warning; unexpected warnings fail tests.
  if (!identical(conditionMessage(w), "The ESS has been capped to avoid unstable estimates.")) {
    stop(w)
  }
  short_fit_warnings <<- c(short_fit_warnings, conditionMessage(w))
  invokeRestart("muffleWarning")
})

## Independent integration over precision verifies the NORMALIZED marginal
## likelihood, including the inverse-gamma prior's parameterization.
X <- matrix(c(-1, 0, 1), ncol = 1L)
Dx <- list(pairwise_sqdist(X))
u <- c(-.7, .3, .9)
y <- c(-.8, .1, .7)
lt <- log(c(2, .7, .4))
spec <- make_theta_spec_multix(1L)
max_integral_error <- 0
for (kernel in c("se", "matern")) {
  for (perturbation in c(-.4, 0, .6)) {
    theta <- lt + perturbation
    state <- gp_state_1d(y, u, Dx, theta, .12, spec, kernel = kernel)
    A <- state$A
    determinant <- as.numeric(determinant(A, logarithm = TRUE)$modulus)
    quadratic <- drop(crossprod(y, solve(A, y)))
    integral <- integrate(function(precision) {
      exp(-length(y) / 2 * log(2 * pi) - determinant / 2 +
            length(y) / 2 * log(precision) - quadratic * precision / 2 +
            dgamma(precision, shape = 2, rate = .05, log = TRUE))
    }, 0, Inf, rel.tol = 1e-10)$value
    actual <- gp_update_loglik_1d(state, length(y), TRUE)
    max_integral_error <- max(max_integral_error, abs(actual - log(integral)))
    stopifnot(abs(actual - log(integral)) < 1e-8,
              identical(gp_update_loglik_1d(state, length(y)), state$loglik))
    stopifnot(abs(theta_logpost_1d(y, u, Dx, theta, NA_real_, spec,
                 kernel = kernel, collapse_noise = TRUE) - actual -
                 log_prior_logtheta(theta, spec)) < 1e-10)
  }
}
cat("Maximum integrated-noise log-density error:", max_integral_error, "\n")

## All proposal paths must ignore stale sigma in marginal mode. Compare the
## SAME seeded transition using two different stale noise values, with both
## dense and exact Schur likelihoods, interior and unbounded ordinal cells.
set.seed(731)
n <- 48L
X48 <- matrix(rnorm(2L * n), n, 2L)
u48 <- sort(rnorm(n))
y48 <- sin(u48) + .2 * X48[, 1L]
D48 <- lapply(1:2, function(j) pairwise_sqdist(X48[, j, drop = FALSE]))
spec48 <- make_theta_spec_multix(2L)
lt48 <- log(c(2, .4, .7, .6))
tau48 <- c(-.5, .5)
c48 <- as.integer(cut(u48, c(-Inf, tau48, Inf), labels = FALSE))
idx48 <- c(which(c48 == 1L)[1L], which(c48 == 2L)[1L], tail(which(c48 == 3L), 1L))
max_block_error <- 0
for (kernel in c("se", "matern")) for (use_blocks in c(FALSE, TRUE)) {
  call_u <- function(fun, sigma, indices, seed) {
    set.seed(seed)
    do.call(fun, c(list(y48, u48, D48, c48, tau48, lt48, sigma, spec48),
      list(indices), list(kernel = kernel, gp_cache = new_gp_cache_1d(use_blocks),
                          collapse_noise = TRUE)))$u
  }
  for (fun in list(update_u_ess_block_1d, update_u_local_z_slice_1d)) {
    stopifnot(identical(call_u(fun, .001, idx48, 732),
                        call_u(fun, NA_real_, idx48, 732)))
  }
  set.seed(733)
  th1 <- update_logtheta_slice_1d(y48, u48, D48, lt48, .001,
       rep(.4, 4L), spec48, kernel = kernel,
       gp_cache = new_gp_cache_1d(use_blocks), collapse_noise = TRUE)
  set.seed(733)
  th2 <- update_logtheta_slice_1d(y48, u48, D48, lt48, NA_real_,
       rep(.4, 4L), spec48, kernel = kernel,
       gp_cache = new_gp_cache_1d(use_blocks), collapse_noise = TRUE)
  stopifnot(identical(th1, th2))
  if (use_blocks) {
    cache <- new_gp_cache_1d(TRUE)
    prepared <- prepare_gp_block_1d(y48, u48, D48, lt48, spec48,
                                    idx48, kernel = kernel, gp_cache = cache)
    proposed <- u48
    proposed[idx48] <- proposed[idx48] + .03
    fast <- gp_block_state_1d(prepared, proposed[idx48], 1, cache)
    dense <- gp_state_1d(y48, proposed, D48, lt48, 1, spec48, kernel = kernel)
    err <- abs(gp_update_loglik_1d(fast, n, TRUE) -
                 gp_update_loglik_1d(dense, n, TRUE))
    max_block_error <- max(max_block_error, err)
    stopifnot(err < 1e-8, cache$block_evaluations > 0L)
    ## A failed Schur shortcut must retry the SAME marginal target densely.
    prepared$chol_rest[,] <- NA_real_
    proposed[idx48] <- proposed[idx48] + .02
    fallback <- gp_block_state_1d(prepared, proposed[idx48], 1, cache)
    direct <- gp_state_1d(y48, proposed, D48, lt48, 1, spec48, kernel = kernel)
    stopifnot(cache$block_fallbacks == 1L,
      abs(gp_update_loglik_1d(fallback, n, TRUE) -
            gp_update_loglik_1d(direct, n, TRUE)) < 1e-8)
  }
}
stopifnot(identical(gp_loglik_with_constraints_1d(y, c(.1, .3, .9), Dx,
  c(1L, 2L, 2L), 0, lt, NA_real_, spec, collapse_noise = TRUE), -Inf))
cat("Maximum collapsed Schur/dense log-density error:", max_block_error, "\n")

## One missing U with fixed thresholds/kernel parameters has a scalar,
## independently integrated posterior. Both sampler strategies must recover
## its U moments and the JOINT U/noise moments, not only a marginal score.
direct_components <- function(z) {
  proposed <- u
  proposed[1L] <- z
  A <- diag(3L) + 4 * exp(-.7 * outer(X[, 1L], X[, 1L], "-")^2 -
                           .4 * outer(proposed, proposed, "-")^2)
  quadratic <- drop(crossprod(y, solve(A, y)))
  log_density <- -.5 * as.numeric(determinant(A, logarithm = TRUE)$modulus) -
    (2 + 3 / 2) * log(.05 + quadratic / 2) + dnorm(z, log = TRUE)
  list(log_density = log_density, mean_sigma = (.05 + quadratic / 2) / (2 + 3 / 2 - 1))
}
integrand <- function(z, moment) vapply(z, function(value) {
  a <- direct_components(value)
  multiplier <- switch(moment, one = 1, u = value, u2 = value^2,
                       sigma = a$mean_sigma, u_sigma = value * a$mean_sigma)
  multiplier * exp(a$log_density)
}, numeric(1))
normalizer <- integrate(function(z) integrand(z, "one"), -Inf, 0,
                        rel.tol = 1e-9)$value
reference <- vapply(c("u", "u2", "sigma", "u_sigma"), function(moment) {
  integrate(function(z) integrand(z, moment), -Inf, 0,
            rel.tol = 1e-9)$value / normalizer
}, numeric(1))
kept <- 5000L
warmup <- 500L
target_results <- list()
for (strategy in c("conditional", "collapsed")) {
  set.seed(if (strategy == "collapsed") 735L else 734L)
  marginal <- strategy == "collapsed"
  state_u <- u
  sigma <- .1
  cache <- new_gp_cache_1d(FALSE)
  draws <- matrix(NA_real_, kept, 4L, dimnames = list(NULL, names(reference)))
  for (iter in seq_len(warmup + kept)) {
    if (!marginal) sigma <- sample_sigma2_eps_1d(y, state_u, Dx, lt, spec,
                                                gp_cache = cache)
    update <- if (iter %% 2L) update_u_local_z_slice_1d else update_u_ess_block_1d
    state_u <- update(y, state_u, Dx, c(1L, 2L, 2L), 0, lt, sigma,
                      spec, 1L, gp_cache = cache, collapse_noise = marginal)$u
    if (marginal) sigma <- sample_sigma2_eps_1d(y, state_u, Dx, lt, spec,
                                               gp_cache = cache)
    if (iter > warmup) draws[iter - warmup, ] <-
      c(state_u[1L], state_u[1L]^2, sigma, state_u[1L] * sigma)
  }
  se <- apply(draws, 2L, batch_mcse)
  error <- abs(colMeans(draws) - reference)
  stopifnot(all(error < 7 * se + 1e-5),
            all(is.finite(draws)), all(draws[, "sigma"] > 0))
  target_results[[strategy]] <- data.frame(strategy = strategy,
    target = names(reference), estimate = colMeans(draws), reference = reference,
    mcse = se, standardized_error = error / se, row.names = NULL)
}
print(do.call(rbind, target_results), row.names = FALSE)

## Instrument the ACTUAL fitter: every retained collapsed state must use the
## sigma draw made at its final U/theta, including all-calibrated data and a
## singleton missing row. Check SE/Matérn, multivariate X, and thinning.
fit_env <- new.env(parent = globalenv())
source(file.path(codes_dir, "00_study1_functions.R"), local = fit_env)
original_noise <- fit_env$sample_sigma2_eps_1d
noise_records <- list()
fit_env$sample_sigma2_eps_1d <- function(y, u, Dx_list, logtheta, ...) {
  value <- original_noise(y, u, Dx_list, logtheta, ...)
  noise_records[[length(noise_records) + 1L]] <<-
    list(u = u, logtheta = logtheta, sigma = value)
  value
}
u8 <- c(-1.5, -1, -.5, -.2, .2, .5, 1, 1.5)
x8 <- cbind(x1 = 1:8, x2 = c(0, 1, 3, 2, 6, 5, 4, 7))
y8 <- sin(u8) + c(0, .1, -.1, .05, 0, -.1, .1, 0)
c8 <- rep(1:2, each = 4L)
for (kernel in c("se", "matern")) for (missing in list(integer(0), 8L, seq_len(8))) {
  obs <- u8
  obs[missing] <- NA_real_
  noise_records <- list()
  fit <- short_fit(fit_env$fit_eivgp_1d(x8, y8, c8, u_obs = obs, m = 2L,
    n_iter = 40L, burn = 10L, thin = 2L, n_chains = 2L,
    parallel_chains = FALSE, seed = 736L, kernel = kernel,
    noise_strategy = "collapsed"))
  stopifnot(identical(fit$control$noise_strategy, "collapsed"),
            length(noise_records) == 80L)
  saved_indices <- c(seq(12L, 40L, by = 2L), 40L + seq(12L, 40L, by = 2L))
  stopifnot(nrow(fit$mcmc$samples_u) == length(saved_indices))
  for (i in seq_along(saved_indices)) {
    recorded <- noise_records[[saved_indices[i]]]
    stopifnot(identical(as.numeric(fit$mcmc$samples_u[i, ]), as.numeric(recorded$u)),
      identical(as.numeric(fit$mcmc$samples_logtheta[i, ]), as.numeric(recorded$logtheta)),
      identical(fit$mcmc$samples_sigma2[i], recorded$sigma),
      check_constraints_1d(fit$mcmc$samples_u[i, ], c8, fit$mcmc$samples_tau[i, ]))
    observed <- setdiff(seq_len(8L), missing)
    stopifnot(identical(as.numeric(fit$mcmc$samples_u[i, observed]), u8[observed]))
  }
}

## End-to-end collapsed fitting exercises the actual Schur-enabled path.
observed48 <- u48
observed48[idx48] <- NA_real_
fit48_args <- list(x_raw = X48, y_raw = y48, c_ord = c48, u_obs = observed48,
  m = 3L, n_iter = 40L, burn = 10L, n_chains = 1L, parallel_chains = FALSE,
  seed = 738L, kernel = "matern", noise_strategy = "collapsed")
fit48_dense <- short_fit(do.call(fit_eivgp_1d,
  c(fit48_args, list(gp_block_schur = FALSE))))
fit48_schur <- short_fit(do.call(fit_eivgp_1d,
  c(fit48_args, list(gp_block_schur = TRUE))))
stopifnot(sum(fit48_schur$mcmc$chain_stats$gp_block_evaluations) > 0L,
          sum(fit48_dense$mcmc$chain_stats$gp_block_evaluations) == 0L)
for (name in c("samples_u", "samples_tau", "samples_logtheta", "samples_sigma2")) {
  stopifnot(max(abs(fit48_dense$mcmc[[name]] - fit48_schur$mcmc[[name]])) < 1e-8)
}

## Default behavior remains the explicit baseline, with exact seeded paths.
args <- list(x_raw = x8, y_raw = y8, c_ord = c8, m = 2L, n_iter = 30L,
             burn = 10L, n_chains = 1L, parallel_chains = FALSE, seed = 737L)
default <- short_fit(do.call(fit_eivgp_1d, args))
baseline <- short_fit(do.call(fit_eivgp_1d, c(args, list(noise_strategy = "conditional"))))
stopifnot(identical(default$control$noise_strategy, "conditional"),
          identical(default$mcmc$samples_by_chain, baseline$mcmc$samples_by_chain))
for (bad in list(NULL, NA_character_, character(0), c("conditional", "collapsed"),
                 TRUE, 1, "marginal", list("collapsed"))) {
  expect_error(do.call(fit_eivgp_1d, c(args, list(noise_strategy = bad))), "noise_strategy")
}
cat("Recorded short-chain ESS-cap warnings:", length(short_fit_warnings), "\n")
cat("PASS: Study I exact noise collapse, target moments, proposal paths, and retained-state order.\n")
