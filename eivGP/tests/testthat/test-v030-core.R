## Bounded implementation checks for the 0.3.0 posterior target.
## Run with Rscript from any directory, or copy into package tests/testthat.
## The few seeded Monte Carlo checks test small known distributions; they do
## not certify convergence or mixing in a scientific application.
if (!exists("mixedgp_v030_collapsed", mode = "function")) {
  cli <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(cli)) stop("Load codes/load_mixedgp.R before sourcing this test.")
  code_dir <- dirname(dirname(normalizePath(sub("^--file=", "", cli[1L]))))
  local({
    old_dir <- setwd(dirname(code_dir))
    on.exit(setwd(old_dir))
    source(file.path(code_dir, "load_mixedgp.R"))
  })
}

testthat::test_that("reference priors and the product dictionary match the target", {
  prior <- mixedgp_v030_priors(list(), p = 2L, d = 2L, m_vec = c(2L, 4L))
  testthat::expect_equal(prior$variance_shape, 3)
  testthat::expect_equal(prior$variance_rate, 2)
  testthat::expect_equal(prior$signal_shape, c(32, 8))
  testthat::expect_equal(prior$log_theta_x_mean, rep(log(.5), 2L))
  testthat::expect_equal(prior$log_theta_x_sd, rep(1.5, 2L))
  testthat::expect_equal(prior$loading_sd, 2.5)
  testthat::expect_identical(prior$dictionary_type, "product")
  grid <- exp(log(.5) + 1.5 * qnorm((seq_len(5L) - .5) / 5))
  testthat::expect_length(prior$u_dictionary, 2L)
  testthat::expect_length(prior$u_weights, 2L)
  for (j in seq_len(2L)) {
    testthat::expect_equal(prior$u_dictionary[[j]], grid)
    testthat::expect_equal(prior$u_weights[[j]], rep(.2, 5L))
  }
  testthat::expect_equal(prior$category_alpha, list(rep(2, 2L), rep(2, 4L)))
  custom <- mixedgp_v030_priors(list(signal_shape = c(16, 4),
    u_dictionary = list(c(.1, 1), c(.2, 2, 4)),
    u_weights = list(c(.3, .7), c(.2, .3, .5))),
    p = 2L, d = 2L, m_vec = c(2L, 4L))
  testthat::expect_equal(custom$signal_shape, c(16, 4))
  testthat::expect_equal(custom$u_weights, list(c(.3, .7), c(.2, .3, .5)))
  dictionary <- mixedgp_v030_dictionary(custom)
  testthat::expect_equal(dim(dictionary$indices), c(6L, 2L))
  testthat::expect_equal(exp(dictionary$logweights),
    as.vector(outer(c(.3, .7), c(.2, .3, .5))), tolerance = 1e-12)
  testthat::expect_error(mixedgp_v030_dictionary(custom, limit = 5L))
  vector_grid <- rbind(c(.1, 1), c(1, .1), c(2, 3))
  vector_prior <- mixedgp_v030_priors(list(u_dictionary = vector_grid,
    u_weights = c(2, 3, 5)), p = 2L, d = 2L, m_vec = c(2L, 4L))
  testthat::expect_identical(vector_prior$dictionary_type, "vector")
  testthat::expect_equal(vector_prior$u_dictionary, vector_grid)
  testthat::expect_equal(vector_prior$u_weights, c(.2, .3, .5))
  testthat::expect_equal(exp(mixedgp_v030_dictionary(vector_prior)$logweights),
    c(.2, .3, .5))
  testthat::expect_error(mixedgp_v030_priors(list(variance_rate = -1), 1L, 1L, 2L))
  testthat::expect_error(mixedgp_v030_priors(list(signal_shape = c(0, 8)), 1L, 1L, 2L))
})

testthat::test_that("collapsed likelihood matches independent scale integration", {
  y <- c(-.8, .2, .7, -.1)
  z <- c(-1, -.2, .5, 1.3)
  R <- exp(-.7 * outer(z, z, "-")^2)
  n <- length(y)
  a <- 3; b <- 2
  constant <- -n / 2 * log(2 * pi) + a * log(b) - lgamma(a) + lgamma(a + n / 2)
  for (r in c(.12, .8, .995)) {
    actual <- mixedgp_v030_collapsed(y, R, r, a, b)
    B <- r * R + (1 - r) * diag(n)
    logdet <- as.numeric(determinant(B, logarithm = TRUE)$modulus)
    Q <- drop(crossprod(y, solve(B, y)))
    expected <- -.5 * logdet - (a + n / 2) * log(b + Q / 2)
    testthat::expect_equal(actual$B, B, tolerance = 1e-12)
    testthat::expect_equal(actual$logdet, logdet, tolerance = 1e-10)
    testthat::expect_equal(actual$Q, Q, tolerance = 1e-10)
    testthat::expect_equal(actual$loglik, expected, tolerance = 1e-10)
    integral <- integrate(function(precision) {
      exp(-n / 2 * log(2 * pi) - logdet / 2 + n / 2 * log(precision) -
        Q * precision / 2 + dgamma(precision, shape = a, rate = b, log = TRUE))
    }, 0, Inf, rel.tol = 1e-10)$value
    testthat::expect_equal(actual$loglik + constant, log(integral), tolerance = 1e-8)
  }
  one <- mixedgp_v030_collapsed(.3, matrix(1, 1L, 1L), .8)
  testthat::expect_equal(one$Q, .09)
  testthat::expect_equal(one$logdet, 0)
  testthat::expect_equal(dim(one$B), c(1L, 1L))
  testthat::expect_error(mixedgp_v030_collapsed(c(0, 0),
    matrix(c(1, 2, 2, 1), 2L), .8), "factorization failed")
  testthat::expect_error(mixedgp_v030_collapsed(y, R, 0))
  testthat::expect_error(mixedgp_v030_collapsed(y, R, 1))
})

testthat::test_that("Dirichlet Gaussian sticks preserve support and inverse coordinates", {
  for (alpha in list(c(2, 2), c(.5, 2, 4), rep(2, 5L))) {
    e <- seq(-2.5, 2.5, length.out = length(alpha) - 1L)
    probability <- mixedgp_v030_stick_forward(e, alpha)
    testthat::expect_length(probability, length(alpha))
    testthat::expect_true(all(is.finite(probability) & probability > 0))
    testthat::expect_equal(sum(probability), 1, tolerance = 1e-13)
    testthat::expect_equal(mixedgp_v030_stick_inverse(probability, alpha), e,
      tolerance = 1e-7)
    sticks <- qbeta(pnorm(e), alpha[-length(alpha)],
      rev(cumsum(rev(alpha)))[-1L])
    independent <- c(sticks, 1) * c(1, cumprod(1 - sticks))
    testthat::expect_equal(probability, independent, tolerance = 1e-12)
  }
  for (e in c(-8, 8)) {
    probability <- mixedgp_v030_stick_forward(e, c(2, 2))
    testthat::expect_true(all(is.finite(probability) & probability > 0))
    testthat::expect_equal(sum(probability), 1, tolerance = 1e-13)
  }
})

testthat::test_that("truncated Gaussian coordinates remain accurate in both tails", {
  z <- c(-3, -.4, 0, .7, 3)
  for (bounds in list(c(-Inf, Inf), c(-Inf, -12), c(12, Inf),
                      c(-35.1, -35), c(35, 35.1), c(-.2, .9))) {
    u <- mixedgp_v030_interval_forward(z, bounds[1L], bounds[2L])
    testthat::expect_true(all(is.finite(u)))
    testthat::expect_true(all(u > bounds[1L] & u < bounds[2L]))
    recovered <- mixedgp_v030_interval_inverse(u, bounds[1L], bounds[2L])
    testthat::expect_equal(recovered, z, tolerance = 2e-7)
    if (all(is.finite(bounds)) && max(abs(bounds)) < 3) {
      direct <- qnorm(pnorm(bounds[1L]) +
        (pnorm(bounds[2L]) - pnorm(bounds[1L])) * pnorm(z))
      testthat::expect_equal(u, direct, tolerance = 1e-12)
    }
  }
  testthat::expect_equal(mixedgp_v030_interval_forward(z, -Inf, Inf), z,
    tolerance = 1e-12)
  testthat::expect_error(mixedgp_v030_interval_forward(0, 1, 1))
  testthat::expect_error(mixedgp_v030_interval_inverse(1, 1, 2))
})

testthat::test_that("Gaussian block ESS recovers the Study I category-count target", {
  set.seed(30491)
  alpha <- c(2, 2, 2)
  counts <- c(12, 3, 0)
  residual <- function(e) sum(counts * log(mixedgp_v030_stick_forward(e, alpha)))
  e <- c(0, 0)
  draws <- matrix(NA_real_, 2600L, 3L)
  for (i in seq_len(nrow(draws))) {
    move <- mixedgp_v030_ess(e, residual)
    e <- move$value
    if (i == 1L) {
      testthat::expect_length(e, 2L)
      testthat::expect_true(is.finite(move$angle))
      testthat::expect_true(move$evaluations >= 1L)
      testthat::expect_equal(move$loglik, residual(e), tolerance = 1e-12)
    }
    draws[i, ] <- mixedgp_v030_stick_forward(e, alpha)
  }
  kept <- draws[-seq_len(400L), , drop = FALSE]
  posterior <- alpha + counts
  expected_mean <- posterior / sum(posterior)
  expected_second <- posterior * (posterior + 1) /
    (sum(posterior) * (sum(posterior) + 1))
  testthat::expect_lt(max(abs(colMeans(kept) - expected_mean)), .055)
  testthat::expect_lt(max(abs(colMeans(kept^2) - expected_second)), .055)
  testthat::expect_error(mixedgp_v030_ess(c(0, 0), residual, max_steps = 0L))
  testthat::expect_error(mixedgp_v030_ess(c(0, 0), function(e) NaN))
  testthat::expect_error(mixedgp_v030_ess(0,
    function(e) if (identical(e, 0)) 0 else -Inf, max_steps = 3L), "exhausted")
})

testthat::test_that("measurement decoding enforces structural loadings and prior scales", {
  prior <- mixedgp_v030_priors(p = 1L, d = 2L, m_vec = c(2L, 3L))
  ctx <- list(priors = prior, n = 3L, q = 2L, d = 2L,
    measurement = "probit", ident = "lower_triangular",
    C = rbind(c(1L, 1L), c(2L, 2L), c(1L, 3L)))
  state <- list(e_pi = list(.2, c(-.3, .5)),
    e_A = rbind(c(.25, 8), c(-.4, .6)), e_r = .3, J = c(1L, 4L),
    U = rbind(c(-.5, .2), c(.1, -.3), c(.8, .6)))
  decoded <- mixedgp_v030_decode(state, ctx)
  testthat::expect_equal(decoded$A[1L, 2L], 0)
  testthat::expect_true(all(diag(decoded$A) > 0))
  testthat::expect_equal(decoded$A[2L, 1L], -1)
  testthat::expect_equal(decoded$r, qbeta(pnorm(.3), 32, 8), tolerance = 1e-12)
  testthat::expect_equal(decoded$theta_u,
    c(prior$u_dictionary[[1L]][1L], prior$u_dictionary[[2L]][4L]))
  internal_tau <- vector("list", 2L)
  for (j in 1:2) {
    cuts <- decoded$tau[[j]]
    testthat::expect_equal(cuts[c(1L, length(cuts))], c(-Inf, Inf))
    internal_tau[[j]] <- cuts[-c(1L, length(cuts))]
    testthat::expect_equal(pnorm(internal_tau[[j]] /
      sqrt(1 + sum(decoded$A[j, ]^2))),
      head(cumsum(decoded$pi[[j]]), -1L), tolerance = 1e-12)
  }
  testthat::expect_equal(mixedgp_v030_measurement_loglik(state, decoded, ctx),
    ordinal_loglik_marginal(ctx$C, state$U, decoded$A, internal_tau), tolerance = 1e-12)
  testthat::expect_equal(mixedgp_v030_measurement_loglik(state, decoded, ctx,
    rows = 2L, items = 1L), ordinal_loglik_marginal(
      ctx$C[2L, 1L, drop = FALSE], state$U[2L, , drop = FALSE],
      decoded$A[1L, , drop = FALSE], internal_tau[1L]), tolerance = 1e-12)
})

testthat::test_that("Study I cutpoint transport uses missing counts and fixes calibration", {
  prior <- mixedgp_v030_priors(p = 1L, d = 1L, m_vec = 3L)
  ctx <- list(priors = prior, n = 4L, p = 1L, q = 1L, d = 1L,
    measurement = "threshold", ident = "none", observed = 1L, missing = 2:4,
    C = matrix(c(1L, 1L, 2L, 3L), ncol = 1L),
    control = list(u_block_size = 2L, cross_every = 1L, ess_max_steps = 100L,
      dictionary_mode = "conditional", threshold_update = "ess"))
  state <- list(e_pi = list(mixedgp_v030_stick_inverse(c(.25, .4, .35), rep(2, 3L))),
    e_A = matrix(0, 1L, 1L), e_r = .2, J = 2L, logtheta_x = log(.7),
    U = matrix(c(-1.5, 0, 0, 0), ncol = 1L))
  decoded <- mixedgp_v030_decode(state, ctx)
  classes <- ctx$C[ctx$missing, 1L]
  original_z <- c(-.3, .1, .6)
  state$U[ctx$missing, 1L] <- mixedgp_v030_interval_forward(original_z,
    decoded$tau[[1L]][classes], decoded$tau[[1L]][classes + 1L])
  ## Isolate the actual sweep's residual closures while holding all but its
  ## cross move fixed. The mocked GP is a known nonconstant function of U.
  isolated <- new.env(parent = environment(mixedgp_v030_sweep))
  sweep <- mixedgp_v030_sweep
  environment(sweep) <- isolated
  calls <- 0L
  cross_state <- NULL
  cross_value <- NA_real_
  last_gp_state <- NULL
  isolated$mixedgp_v030_prepare_block <- function(...) NULL
  isolated$mixedgp_v030_record_ess <- function(...) invisible(NULL)
  isolated$mixedgp_v030_update_dictionary <- function(state, ctx) state
  isolated$mixedgp_v030_loggp <- function(state, ctx, prepared = NULL) {
    last_gp_state <<- state
    -.1 * sum(state$U^2)
  }
  isolated$mixedgp_v030_ess <- function(v, loglik, max_steps) {
    calls <<- calls + 1L
    candidate <- v
    if (calls == 3L) {
      ## Fixed-physical-U cutpoint updates have indicator-only residuals.
      testthat::expect_equal(loglik(v), 0)
      testthat::expect_identical(loglik(c(-8, 0)), -Inf)
    }
    if (calls == 4L) {
      candidate <- v + c(.15, -.07)
      cross_value <<- loglik(candidate)
      cross_state <<- last_gp_state
    }
    list(value = candidate, loglik = loglik(candidate), evaluations = 1L, angle = .1)
  }
  set.seed(391L)
  result <- sweep(state, ctx, iteration = 1L)
  testthat::expect_equal(calls, 5L)
  testthat::expect_equal(result$U[1L, 1L], state$U[1L, 1L], tolerance = 0)
  changed <- mixedgp_v030_decode(cross_state, ctx)
  testthat::expect_true(all(abs(result$U[ctx$missing, 1L] -
    state$U[ctx$missing, 1L]) > 1e-6))
  testthat::expect_equal(mixedgp_v030_interval_inverse(result$U[ctx$missing, 1L],
    changed$tau[[1L]][classes], changed$tau[[1L]][classes + 1L]), original_z,
    tolerance = 1e-10)
  expected <- -.1 * sum(cross_state$U^2) + sum(log(changed$pi[[1L]][classes]))
  testthat::expect_equal(cross_value, expected, tolerance = 1e-12)
})

testthat::test_that("marginal probit category probabilities use the loading-dependent scale", {
  probability <- c(.12, .31, .57)
  for (loading in c(-2.2, 0, 1.3)) {
    cuts <- sqrt(1 + loading^2) * qnorm(cumsum(probability)[1:2])
    observed <- vapply(seq_along(probability), function(category) {
      integrate(function(u) vapply(u, function(value) {
        exp(ordinal_loglik_marginal(matrix(category, 1L, 1L),
          matrix(value, 1L, 1L), matrix(loading, 1L, 1L), list(cuts))) * dnorm(value)
      }, numeric(1L)), -Inf, Inf, rel.tol = 1e-9)$value
    }, numeric(1L))
    testthat::expect_equal(observed, probability, tolerance = 1e-8)
  }
  testthat::expect_true(all(is.finite(log_normal_interval_prob(
    c(35, -35.1), c(35.1, -35)))))
})

testthat::test_that("Schur likelihoods agree with the same dense covariance", {
  n <- 8L
  prior <- mixedgp_v030_priors(p = 1L, d = 2L, m_vec = c(2L, 3L))
  state <- list(e_pi = list(.2, c(-.3, .5)), e_A = matrix(.3, 2L, 2L),
    e_r = .3, J = c(1L, 4L), logtheta_x = log(.7),
    U = cbind(seq(-1, 1, length.out = n), sin(seq_len(n))))
  ctx <- list(priors = prior, n = n, p = 1L, q = 2L, d = 2L,
    measurement = "probit", ident = "none", kernel = "se", matern_nu = 2.5,
    y = cos(seq_len(n)) / 2, X = matrix(seq(-.7, .8, length.out = n), ncol = 1L),
    control = list(gp_block_schur = TRUE, dictionary_mode = "conditional"),
    cache = new.env(parent = emptyenv()), counts = list2env(list(
      gp_block_evaluations = 0L, gp_block_fallbacks = 0L,
      gp_full_factorizations = 0L, gp_block_setups = 0L)))
  for (kernel in c("se", "matern")) for (rows in list(2L, c(2L, 7L))) {
    ctx$kernel <- kernel
    ctx$cache$key <- NULL
    prepared <- mixedgp_v030_prepare_block(state, ctx, rows)
    proposal <- state
    proposal$U[rows, ] <- proposal$U[rows, , drop = FALSE] + .13
    fast <- mixedgp_v030_gp(proposal, ctx, prepared)
    decoded <- mixedgp_v030_decode(proposal, ctx)
    dense <- mixedgp_v030_collapsed(ctx$y,
      mixedgp_v030_corr(proposal, decoded, ctx), decoded$r)
    testthat::expect_equal(fast$loglik, dense$loglik, tolerance = 1e-10)
    testthat::expect_equal(fast$Q, dense$Q, tolerance = 1e-10)
    testthat::expect_equal(fast$logdet, dense$logdet, tolerance = 1e-10)
  }
  testthat::expect_equal(ctx$counts$gp_block_evaluations, 4L)
  prepared <- mixedgp_v030_prepare_block(state, ctx, c(2L, 7L))
  proposal <- state
  proposal$U[c(2L, 7L), ] <- proposal$U[c(2L, 7L), , drop = FALSE] + .13
  reference <- function(st) {
    dec <- mixedgp_v030_decode(st, ctx)
    mixedgp_v030_collapsed(ctx$y, mixedgp_v030_corr(st, dec, ctx), dec$r)
  }
  for (failure in c("singular_factor", "nonfinite_moment")) {
    broken <- prepared
    if (failure == "singular_factor") broken$chol[,] <- 0 else broken$Q <- NaN
    ctx$cache$key <- NULL
    before <- ctx$counts$gp_block_fallbacks
    actual <- mixedgp_v030_gp(proposal, ctx, broken)
    direct <- reference(proposal)
    testthat::expect_equal(ctx$counts$gp_block_fallbacks, before + 1L)
    testthat::expect_equal(actual$B, direct$B, tolerance = 0)
    testthat::expect_equal(actual$loglik, direct$loglik, tolerance = 1e-12)
  }
  for (changed in c("complement", "theta_x", "dictionary", "r")) {
    stale <- proposal
    if (changed == "complement") stale$U[1L, ] <- stale$U[1L, ] + .2
    if (changed == "theta_x") stale$logtheta_x <- stale$logtheta_x + .2
    if (changed == "dictionary") stale$J[1L] <- 2L
    if (changed == "r") stale$e_r <- stale$e_r + .2
    ctx$cache$key <- NULL
    before_dense <- ctx$counts$gp_full_factorizations
    before_block <- ctx$counts$gp_block_evaluations
    actual <- mixedgp_v030_gp(stale, ctx, prepared)
    testthat::expect_equal(ctx$counts$gp_full_factorizations, before_dense + 1L)
    testthat::expect_equal(ctx$counts$gp_block_evaluations, before_block)
    testthat::expect_equal(actual$loglik, reference(stale)$loglik, tolerance = 1e-12)
  }
  original <- mixedgp_v030_gp(state, ctx)
  measurement_changed <- state
  measurement_changed$e_A <- state$e_A + .3
  measurement_changed$e_pi <- lapply(state$e_pi, function(e) e + .2)
  before_dense <- ctx$counts$gp_full_factorizations
  testthat::expect_identical(mixedgp_v030_gp(measurement_changed, ctx), original)
  testthat::expect_equal(ctx$counts$gp_full_factorizations, before_dense)
})

testthat::test_that("prediction covariance agrees with recovered V and r directly", {
  X <- matrix(c(-1, .1, 1), ncol = 1L)
  U <- matrix(c(-.6, .2, .8), ncol = 1L)
  y <- c(-.5, .1, .7)
  Xstar <- matrix(c(-.2, .4), ncol = 1L)
  Ustar <- matrix(c(-.1, .6), ncol = 1L)
  V <- 1.7; r <- .76; theta_x <- .7; theta_u <- .9
  legacy_theta <- c(.5 * log(r / (1 - r)), log(theta_x), log(theta_u))
  noise <- (1 - r) * V
  correlation <- function(xa, ua, xb = xa, ub = ua) {
    exp(-theta_x * outer(xa[, 1L], xb[, 1L], "-")^2 -
      theta_u * outer(ua[, 1L], ub[, 1L], "-")^2)
  }
  C <- V * (r * correlation(X, U) + (1 - r) * diag(nrow(X)))
  Kstar <- V * r * correlation(Xstar, Ustar, X, U)
  expected_mean <- drop(Kstar %*% solve(C, y))
  expected_cov <- V * r * correlation(Xstar, Ustar) -
    Kstar %*% solve(C, t(Kstar))
  fixed <- gp_predict_draw_general(X, U, y, Xstar, Ustar,
    legacy_theta, noise, return_cov = TRUE)
  testthat::expect_equal(fixed$mean, expected_mean, tolerance = 1e-10)
  testthat::expect_equal(fixed$cov, expected_cov, tolerance = 1e-10)
  noisy <- gp_predict_draw_general(X, U, y, Xstar, Ustar,
    legacy_theta, noise, noisy = TRUE, return_cov = TRUE)
  testthat::expect_equal(noisy$cov - fixed$cov, noise * diag(2L), tolerance = 1e-10)
  univariate <- gp_predict_draw(X, drop(U), y, Xstar, drop(Ustar),
    legacy_theta, noise, return_cov = TRUE)
  testthat::expect_equal(univariate$mean, expected_mean, tolerance = 1e-10)
  testthat::expect_equal(univariate$cov, expected_cov, tolerance = 1e-10)

  Umc <- array(c(-.8, -.1, .5, -.3, .4, 1), c(3L, 2L, 1L))
  expanded_X <- Xstar[rep(seq_len(2L), each = 3L), , drop = FALSE]
  expanded_U <- matrix(Umc[, , 1L], ncol = 1L)
  weights <- kronecker(diag(2L), matrix(1 / 3, 1L, 3L))
  expanded_K <- V * r * correlation(expanded_X, expanded_U, X, U)
  expanded_cov <- V * r * correlation(expanded_X, expanded_U) -
    expanded_K %*% solve(C, t(expanded_K))
  integrated <- gp_integrated_mean_state_general(X, U, y, Xstar, Umc,
    legacy_theta, noise, return_cov = TRUE)
  testthat::expect_equal(integrated$mean,
    drop(weights %*% expanded_K %*% solve(C, y)), tolerance = 1e-10)
  testthat::expect_equal(integrated$cov,
    weights %*% expanded_cov %*% t(weights), tolerance = 1e-10)
})

testthat::test_that("tiny fits preserve calibration, absent levels and recovered quantities", {
  configurations <- list(
    c("threshold", "product", "conditional", "singleton"),
    c("threshold", "product", "marginal", "none"),
    c("threshold", "vector", "conditional", "all"),
    c("probit", "product", "conditional", "none"),
    c("probit", "product", "marginal", "singleton"),
    c("probit", "vector", "conditional", "all"),
    c("probit", "vector", "marginal", "none"))
  n <- 6L
  X <- matrix(seq(-1, 1, length.out = n), ncol = 1L)
  U <- cbind(seq(-1.2, 1.2, length.out = n), cos(seq_len(n)))
  for (config in configurations) {
    probit <- config[1L] == "probit"
    d <- if (probit) 2L else 1L
    C <- matrix(c(1L, 1L, 1L, 3L, 3L, 3L), ncol = 1L)
    if (probit) C <- cbind(C, c(2L, 1L, 2L, 1L, 2L, 1L))
    q <- ncol(C)
    observed <- switch(config[4L], none = integer(0), singleton = 1L, all = seq_len(n))
    Uobs <- U[, seq_len(d), drop = FALSE]
    Uobs[setdiff(seq_len(n), observed), ] <- NA_real_
    dictionary <- if (config[2L] == "product") rep(list(c(.2, .8)), d) else
      matrix(c(rep(.2, d), rep(.8, d)), nrow = 2L, byrow = TRUE)
    result <- mixedgp_v030_fit(X, sin(U[, 1L]) + X[, 1L] / 3, C, Uobs,
      m_vec = rep(3L, q), measurement = config[1L],
      ident = if (probit) "lower_triangular" else "none",
      n_iter = 12L, burn = 4L, n_chains = 1L, seed = 736L,
      parallel_chains = FALSE, priors = list(u_dictionary = dictionary),
      sampler_control = list(u_block_size = 2L, cross_every = 2L,
        dictionary_mode = config[3L]), store_scores = probit)
    chain <- result$chains[[1L]]
    testthat::expect_identical(dim(chain$samples_U), c(8L, n, d))
    testthat::expect_true(all(is.finite(chain$samples_U)))
    testthat::expect_true(all(is.finite(chain$samples_V) & chain$samples_V > 0))
    testthat::expect_true(all(chain$samples_r > 0 & chain$samples_r < 1))
    testthat::expect_equal(chain$samples_sigma2,
      chain$samples_V * (1 - chain$samples_r), tolerance = 1e-13)
    testthat::expect_equal(exp(2 * chain$samples_logtheta[, 1L]),
      chain$samples_r / (1 - chain$samples_r), tolerance = 1e-12)
    testthat::expect_false(any(c("S", "V", "sigma2") %in% names(chain$checkpoint$state)))
    for (i in observed) for (k in seq_len(d)) {
      testthat::expect_equal(chain$samples_U[, i, k], rep(Uobs[i, k], 8L), tolerance = 0)
    }
    for (k in seq_len(d)) {
      expected <- if (config[2L] == "product")
        dictionary[[k]][chain$samples_J[, k]] else dictionary[chain$samples_J[, 1L], k]
      testthat::expect_equal(unname(exp(chain$samples_logtheta[, 2L + k])),
        expected, tolerance = 1e-12)
    }
    for (j in seq_len(q)) {
      pi_columns <- (j - 1L) * 3L + seq_len(3L)
      tau_columns <- (j - 1L) * 2L + seq_len(2L)
      probabilities <- chain$samples_pi[, pi_columns, drop = FALSE]
      testthat::expect_equal(rowSums(probabilities), rep(1, 8L), tolerance = 1e-12)
      scale <- if (probit) sqrt(1 + rowSums(matrix(chain$samples_A[, j, ], 8L, d)^2)) else rep(1, 8L)
      cumulative <- t(apply(probabilities, 1L, cumsum))[, 1:2, drop = FALSE]
      cuts <- chain$samples_tau[, tau_columns, drop = FALSE]
      testthat::expect_equal(unname(cuts), qnorm(cumulative) * scale, tolerance = 1e-9)
      for (i in seq_len(n)) {
        lo <- if (C[i, j] == 1L) rep(-Inf, 8L) else cuts[, C[i, j] - 1L]
        hi <- if (C[i, j] == 3L) rep(Inf, 8L) else cuts[, C[i, j]]
        value <- if (probit) chain$samples_S[, i, j] else chain$samples_U[, i, 1L]
        testthat::expect_true(all(value > lo & value <= hi))
      }
    }
    if (probit) {
      testthat::expect_true(all(chain$samples_A[, 1L, 2L] == 0))
      testthat::expect_true(all(chain$samples_A[, 1L, 1L] > 0))
      testthat::expect_true(all(chain$samples_A[, 2L, 2L] > 0))
    }
  }
})
