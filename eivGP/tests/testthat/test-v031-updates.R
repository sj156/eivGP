## Exact conditional and regression checks for the 0.3.1 transition updates.
if (!exists("mixedgp_v031_threshold_gibbs", mode = "function")) {
  cli <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  code_dir <- dirname(dirname(normalizePath(sub("^--file=", "", cli[1L]))))
  source(file.path(code_dir, "load_mixedgp.R"), chdir = TRUE)
}

testthat::test_that("diagnostic exports share roles and preserve undefined diagnostics", {
  out <- tempfile("diagnostic-export-")
  detail <- data.frame(parameter = c("u", "invariant", "m", "contrast"),
    group = c("raw_and_training_imputation", "measurement_invariant", "scientific_panel", "user_functional"),
    rhat = c(NA, 1, 1.2, 1), ess_bulk = c(NA, 500, 10, 600))
  obj <- structure(list(status = "poor_exploration", recommendation = "Inspect chains",
    table = detail, settings = list(full_retained_window = TRUE)), class = "eivgp_diagnostics")
  paths <- write_diagnostics_eivgp(obj, out)
  testthat::expect_true(all(file.exists(paths)))
  testthat::expect_identical(readRDS(paths[["report"]]), obj)
  raw <- read.csv(paths[["mcmc_parameter_diagnostics"]])
  testthat::expect_equal(raw$parameter, "u")
  testthat::expect_true(is.na(raw$rhat))
  testthat::expect_equal(nrow(read.csv(paths[["mcmc_target_diagnostics"]])), 3L)
  testthat::expect_error(write_diagnostics_eivgp(obj, out), "exists")
  testthat::expect_error(write_diagnostics_eivgp(obj, out, prefix = "../bad"), "prefix")
  testthat::expect_silent(write_diagnostics_eivgp(obj, out, overwrite = TRUE))
  legacy <- data.frame(parameter = c("u", "m", "inv"),
    target = c("raw_coordinate", "conditional_mean", "ordinal_measurement_invariant"))
  p <- mixedgp_write_diagnostic_tables(data.frame(n = 1), legacy, out, "study2")
  testthat::expect_equal(read.csv(p[["mcmc_parameter_diagnostics"]])$parameter, "u")
  testthat::expect_equal(nrow(read.csv(p[["mcmc_target_diagnostics"]])), 2L)
  p <- mixedgp_write_diagnostic_tables(data.frame(), data.frame(), out, "empty")
  testthat::expect_equal(nrow(read.csv(p[["mcmc_diagnostic_details"]])), 0L)
})

testthat::test_that("truncated Beta quantiles match an independent conditional CDF", {
  for (ab in list(c(.7, 3), c(2, 2), c(5, .8))) {
    for (bounds in list(c(0, 1), c(.1, .7), c(1e-8, 2e-7), c(.99999, .999999))) {
      for (z in c(.1, .5, .9)) {
        x <- mixedgp_v031_truncated_beta(bounds[1], bounds[2], ab, z)
        testthat::expect_true(x > bounds[1] && x < bounds[2])
        ## Integrate on the interval's rescaled coordinate to avoid cancellation.
        width <- diff(bounds)
        center <- mean(bounds)
        density <- function(t) exp(dbeta(bounds[1] + width * t, ab[1], ab[2], log = TRUE) -
                                    dbeta(center, ab[1], ab[2], log = TRUE))
        actual <- integrate(density, 0, (x - bounds[1]) / width, rel.tol = 1e-9)$value /
          integrate(density, 0, 1, rel.tol = 1e-9)$value
        testthat::expect_equal(actual, z, tolerance = 2e-7)
      }
    }
  }
  testthat::expect_error(mixedgp_v031_truncated_beta(.2, .2, c(2, 2)))
  testthat::expect_error(mixedgp_v031_truncated_beta(0, 1, c(0, 2)))
})

v031_threshold_context <- function(pi, alpha, u, categories) {
  ctx <- list(priors = mixedgp_v030_priors(list(category_alpha = alpha),
    p = 1L, d = 1L, m_vec = length(alpha)), measurement = "threshold", ident = "none",
    n = length(u), q = 1L, d = 1L, C = matrix(categories, ncol = 1L),
    counts = new.env(parent = emptyenv()))
  state <- list(U = matrix(u, ncol = 1L), e_pi = list(mixedgp_v030_stick_inverse(pi, alpha)),
    e_A = matrix(0, 1L, 1L), e_r = 0, J = 1L, logtheta_x = 0)
  list(state = state, ctx = ctx)
}

testthat::test_that("cutoff sweep equals sequential truncated Dirichlet splits", {
  pi <- c(.2, .3, .1, .4); alpha <- c(.7, 2, 3, 1.2)
  ## Empty category 3, but all observations still constrain the cutoffs.
  obj <- v031_threshold_context(pi, alpha, c(-1.5, -.6, -.1, 1), c(1L, 2L, 2L, 4L))
  set.seed(31); uniforms <- runif(3)
  expected <- pi
  for (k in 1:3) {
    P <- c(0, cumsum(expected)); total <- expected[k] + expected[k + 1L]
    lo <- max(c(P[k], pnorm(obj$state$U[obj$ctx$C[, 1] <= k, 1])))
    hi <- min(c(P[k + 2L], pnorm(obj$state$U[obj$ctx$C[, 1] > k, 1])))
    l <- (lo - P[k]) / total; h <- (hi - P[k]) / total
    flip <- expected[k] > expected[k + 1L]
    z <- if (flip) 1 - uniforms[k] else uniforms[k]
    v <- qbeta(pbeta(l, alpha[k], alpha[k + 1L]) + z *
      (pbeta(h, alpha[k], alpha[k + 1L]) - pbeta(l, alpha[k], alpha[k + 1L])),
      alpha[k], alpha[k + 1L])
    expected[c(k, k + 1L)] <- total * c(v, 1 - v)
  }
  set.seed(31); updated <- mixedgp_v031_threshold_gibbs(obj$state, obj$ctx)
  dec <- mixedgp_v030_decode(updated, obj$ctx)
  testthat::expect_equal(dec$pi[[1L]], expected, tolerance = 1e-11)
  testthat::expect_identical(updated$U, obj$state$U)
  testthat::expect_identical(updated$J, obj$state$J)
  testthat::expect_equal(obj$ctx$counts$threshold_gibbs_updates, 3)
  testthat::expect_equal(mixedgp_v030_measurement_loglik(updated, dec, obj$ctx), 0)
})

testthat::test_that("binary cutoff Gibbs recovers its exact conditional distribution", {
  alpha <- c(2, 3)
  obj <- v031_threshold_context(c(.5, .5), alpha, c(-.8, 1.2), c(1L, 2L))
  set.seed(3191)
  P <- replicate(1600, {
    obj$state <- mixedgp_v031_threshold_gibbs(obj$state, obj$ctx)
    mixedgp_v030_decode(obj$state, obj$ctx)$pi[[1L]][1L]
  })
  lo <- pnorm(-.8); hi <- pnorm(1.2)
  pit <- (pbeta(P, 2, 3) - pbeta(lo, 2, 3)) / (pbeta(hi, 2, 3) - pbeta(lo, 2, 3))
  testthat::expect_true(all(P > lo & P < hi))
  testthat::expect_equal(mean(pit), .5, tolerance = .025)
  testthat::expect_equal(var(pit), 1 / 12, tolerance = .01)
  testthat::expect_equal(as.numeric(quantile(pit, c(.1, .5, .9))), c(.1, .5, .9), tolerance = .04)
  ## A tight calibrated interval must remain valid without GP calls.
  tight <- v031_threshold_context(c(.5, .5), c(2, 2), c(-1e-7, 1e-7), c(1L, 2L))
  draw <- mixedgp_v031_threshold_gibbs(tight$state, tight$ctx)
  testthat::expect_true(abs(mixedgp_v030_decode(draw, tight$ctx)$tau[[1]][2]) < 1e-7)
})

testthat::test_that("probit joint likelihood subset has the same candidate differences", {
  prior <- mixedgp_v030_priors(p = 1L, d = 2L, m_vec = c(3L, 4L, 2L))
  ctx <- list(priors = prior, n = 5L, q = 3L, d = 2L, measurement = "probit", ident = "none",
    C = cbind(c(1, 2, 3, 1, 2), c(1, 3, 2, 4, 2), c(1, 2, 1, 2, 1)))
  state <- list(e_pi = list(c(0, .2), c(-.3, .2, .6), -.2), e_A = matrix(.3, 3L, 2L),
    U = cbind(seq(-1, 1, length.out = 5), cos(1:5)), e_r = 0, J = c(1L, 2L))
  rows <- c(2L, 4L)
  for (j in 1:3) {
    proposed <- state
    proposed$e_A[j, ] <- proposed$e_A[j, ] + .2
    proposed$e_pi[[j]] <- proposed$e_pi[[j]] - .3
    proposed$U[rows, ] <- proposed$U[rows, ] + .4
    full <- function(st) mixedgp_v030_measurement_loglik(st, mixedgp_v030_decode(st, ctx), ctx)
    affected <- function(st) {
      dec <- mixedgp_v030_decode(st, ctx)
      mixedgp_v030_measurement_loglik(st, dec, ctx, items = j) +
        mixedgp_v030_measurement_loglik(st, dec, ctx, rows = rows, items = setdiff(1:3, j))
    }
    testthat::expect_equal(full(proposed) - full(state), affected(proposed) - affected(state), tolerance = 1e-12)
  }
})

testthat::test_that("real sweeps use cutoff Gibbs and prepare the probit joint shortcut", {
  X <- matrix(seq(-1, 1, length.out = 8), ncol = 1L)
  u <- as.numeric(X); C <- matrix(as.integer(u > 0) + 1L, ncol = 1L)
  observed <- matrix(u, ncol = 1L); observed[3:6, ] <- NA_real_
  for (measurement in c("threshold", "probit")) {
    fit <- mixedgp_v030_fit(X, sin(u), C, observed, 2L, measurement = measurement,
      n_iter = 8L, burn = 2L, n_chains = 1L, parallel_chains = FALSE,
      sampler_control = list(u_block_size = 2L, cross_every = 2L), seed = 903)
    stats <- fit$chains[[1]]$stats
    testthat::expect_equal(stats$threshold_gibbs_updates, if (measurement == "threshold") 8 else 0)
    testthat::expect_equal(stats$measurement_ess_total, if (measurement == "threshold") 0 else 8)
    testthat::expect_equal(stats$gp_block_setups, if (measurement == "threshold") 16 else 20)
    testthat::expect_equal(stats$cross_ess_total, 4)
    testthat::expect_true(all(fit$chains[[1]]$samples_U[, 1, 1] == u[1]))
  }
  testthat::expect_error(mixedgp_v030_control(list(threshold_update = "unknown"), 4))
})
