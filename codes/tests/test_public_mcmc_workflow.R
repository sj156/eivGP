## Run with Rscript from any directory, or copy into package tests/testthat.
if (!exists("continue_eivgp", mode = "function")) {
  cli <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  code_dir <- dirname(dirname(normalizePath(sub("^--file=", "", cli[1L]))))
  for (module in c("00_parallel_utils.R", "00_study1_functions.R",
                   "00_study2_functions.R", "00_public_api.R",
                   "00_diagnostics.R", "00_mcmc_workflow.R")) {
    source(file.path(code_dir, module))
  }
}

workflow_arguments <- function(engine) {
  X <- cbind(x1 = seq(-1, 1, length.out = 10L), x2 = cos(seq_len(10L)))
  u <- seq(-1.2, 1.2, length.out = 10L)
  C <- matrix(as.integer(u > 0) + 1L, ncol = 1L)
  U <- matrix(u, ncol = 1L)
  if (engine == "multivariate") {
    C <- cbind(C, rep(1:2, 5L))
    U <- cbind(U, sin(seq_len(10L)))
  }
  U[3:8, ] <- NA_real_
  out <- list(X = X, y = sin(u) + X[, 2L] / 5, C = C, U_obs = U,
    engine = engine, parallel = FALSE, n_chains = 2L, burn = 205L,
    preset = "fast", seed = 913L)
  if (engine == "multivariate") out$ident <- "none"
  out
}

testthat::test_that("core budgets allocate datasets and chains without oversubscription", {
  for (budget in c(1L, 2L, 3L, 4L, 8L, 12L, 16L)) {
    plan <- eivgp_run_settings(budget, pending_datasets = 10L)
    testthat::expect_lte(plan$workers * plan$chain_workers, budget)
    testthat::expect_identical(plan$fit_args$n_iter, 1750L)
    testthat::expect_identical(plan$fit_args$burn, 500L)
    testthat::expect_equal(plan$retained_draws, 5000)
    testthat::expect_false(plan$automatic_continuation)
  }
  if (.Platform$OS.type != "windows") {
    testthat::expect_identical(eivgp_run_settings(12L, 10L)$workers, 3L)
    testthat::expect_identical(eivgp_run_settings(16L, 10L)$workers, 4L)
  }
  testthat::expect_identical(eivgp_run_settings(16L, 1L)$workers, 1L)
  testthat::expect_identical(eivgp_run_settings(16L, 0L)$workers, 0L)
  for (bad in list(0, -1, 1.5, NA_real_, Inf, "12")) {
    testthat::expect_error(eivgp_run_settings(bad))
  }
})

testthat::test_that("both engines continue exactly, without thinning or re-warmup", {
  for (variant in c("conditional", "collapsed", "legacy", "interwoven")) {
    engine <- if (variant %in% c("conditional", "collapsed")) "univariate" else "multivariate"
    args <- workflow_arguments(engine)
    if (engine == "univariate") {
      args$noise_strategy <- variant
    } else {
      args$sampler_strategy <- variant
      args$store_scores <- variant == "legacy"
      if (variant == "interwoven") args$control_overrides <- list(
        joint_measurement_every = 2L, loading_transport_every = 5L)
    }
    if (variant %in% c("collapsed", "interwoven")) args$kernel <- "matern"
    first <- suppressWarnings(do.call(fit_eivgp, c(args, list(n_iter = 243L))))
    path <- tempfile(fileext = ".rds")
    saveRDS(first, path)
    first <- readRDS(path)
    unlink(path)
    set.seed(488)
    rng <- mixedgp_rng_state()
    longer <- suppressWarnings(continue_eivgp(first, 42L, parallel = FALSE))
    testthat::expect_identical(mixedgp_rng_state(), rng)
    direct <- suppressWarnings(do.call(fit_eivgp, c(args, list(n_iter = 285L))))
    testthat::expect_identical(longer$mcmc$samples_by_chain, direct$mcmc$samples_by_chain)
    testthat::expect_identical(longer$checkpoint$states, direct$checkpoint$states)
    testthat::expect_identical(longer$mcmc$chain_stats, direct$mcmc$chain_stats)
    testthat::expect_equal(longer$diagnostics$summary$saved_per_chain, 80L)
    testthat::expect_equal(longer$diagnostics$summary$thin, 1L)
    for (i in 1:2) testthat::expect_identical(
      first$checkpoint$states[[i]]$state$theta_slice_width,
      longer$checkpoint$states[[i]]$state$theta_slice_width)
    parallel_fit <- suppressWarnings(continue_eivgp(first, 42L, parallel = TRUE, n_cores = 2L))
    testthat::expect_identical(parallel_fit$mcmc$samples_by_chain, longer$mcmc$samples_by_chain)
    twice <- suppressWarnings(continue_eivgp(
      suppressWarnings(continue_eivgp(first, 1L, parallel = FALSE)), 41L, parallel = FALSE))
    testthat::expect_identical(twice$mcmc$samples_by_chain, longer$mcmc$samples_by_chain)
    testthat::expect_identical(twice$checkpoint$states, longer$checkpoint$states)
    testthat::expect_identical(first$data, longer$data)
  }
})

testthat::test_that("continuation rejects modified and unsupported fits", {
  args <- workflow_arguments("univariate")
  args$burn <- 10L
  first <- suppressWarnings(do.call(fit_eivgp, c(args, list(n_iter = 22L))))
  old <- first; old$checkpoint <- NULL
  testthat::expect_error(continue_eivgp(old, 5L), "no supported")
  modified <- first; modified$data$y[1L] <- 9
  testthat::expect_error(continue_eivgp(modified, 5L), "modified")
  modified <- first; modified$mcmc$samples_by_chain$u[[1L]][1L, 1L] <- 9
  testthat::expect_error(continue_eivgp(modified, 5L), "modified")
  testthat::expect_error(continue_eivgp(first, 1.5), "integer")
  testthat::expect_error(continue_eivgp(first, 5L, seed = 2L), "unused argument")
  testthat::expect_error(do.call(fit_eivgp, c(args, list(n_iter = 22L, thin = 2L))), "every post-warm-up")
  testthat::expect_error(do.call(fit_eivgp, c(args, list(n_iter = 22L, .resume = list()))), "internal")
})

testthat::test_that("diagnostic classification separates disagreement from precision", {
  set.seed(981)
  iid <- list(mu = replicate(4L, rnorm(2000L), simplify = FALSE))
  tab <- mixedgp_summarize_diagnostic_series(iid)
  testthat::expect_identical(mixedgp_diagnostic_status(tab, 1.01, 400, 400, .05), "screen_passed")
  testthat::expect_identical(mixedgp_diagnostic_status(tab, 1.01, 400, 400, 1e-7), "insufficient_precision")
  tab$rhat <- 1.2
  testthat::expect_identical(mixedgp_diagnostic_status(tab, 1.01, 400, 400, .05), "poor_exploration")
  tab$rhat <- NA_real_
  testthat::expect_identical(mixedgp_diagnostic_status(tab, 1.01, 400, 400, .05), "insufficient_diagnostics")
})

testthat::test_that("public diagnostics use all retained draws and actual functionals", {
  for (engine in c("univariate", "multivariate")) {
    args <- workflow_arguments(engine)
    args$burn <- 10L; args$n_chains <- 4L
    fit <- suppressWarnings(do.call(fit_eivgp, c(args, list(n_iter = 30L))))
    X <- args$X[1:2, , drop = FALSE]
    C <- args$C[1:2, , drop = FALSE]
    U <- args$U_obs[1:2, , drop = FALSE]
    custom <- list(noise_variance = fit$mcmc$samples_by_chain$sigma2)
    set.seed(533)
    before <- mixedgp_rng_state()
    out <- suppressWarnings(diagnose_eivgp(fit, X, C, U, n_latent = 4L,
      additional_series = custom))
    testthat::expect_identical(mixedgp_rng_state(), before)
    testthat::expect_true(all(out$table$min_draws_per_chain == 20L))
    testthat::expect_true(all(out$table$draw_window_complete))
    testthat::expect_true(any(startsWith(out$table$parameter, "f_conditional_mean")))
    testthat::expect_true(any(grepl("\\^2$", out$table$parameter)))
    row <- out$table[out$table$parameter == "noise_variance", ]
    testthat::expect_equal(row$mcse, as.numeric(posterior::mcse_mean(do.call(cbind, custom$noise_variance))))
    testthat::expect_error(diagnose_eivgp(fit, X, C, n_latent = 4L,
      additional_series = list(bad = list(1:2))), "all retained")
    testthat::expect_error(diagnose_eivgp(fit, X, NULL), "both X and C")
    testthat::expect_error(diagnose_eivgp(fit, n_latent = 1L), "integer")
    testthat::expect_output(print(out), "EIV-GP diagnostic screen")
  }
})
