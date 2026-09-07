## Canonical 0.3.0 public workflow tests, also copied into installed package tests.
if (!exists("mixedgp_v030_fit", mode = "function")) {
  cli <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  code_dir <- dirname(dirname(normalizePath(sub("^--file=", "", cli[1L]))))
  source(file.path(code_dir, "load_mixedgp.R"), chdir = TRUE)
}

v030_workflow_args <- function(engine, dictionary_mode = "conditional") {
  n <- 10L
  u <- seq(-1.2, 1.2, length.out = n)
  U <- matrix(u, ncol = 1L, dimnames = list(NULL, "latent1"))
  U[3:8, ] <- NA_real_
  C <- matrix(as.integer(u > 0) + 1L, ncol = 1L, dimnames = list(NULL, "ordinal1"))
  if (engine == "multivariate") C <- cbind(C, ordinal2 = rep(1:2, 5L))
  args <- list(X = cbind(exposure = seq(-1, 1, length.out = n),
                        context = cos(seq_len(n))),
    y = sin(u) + cos(seq_len(n)) / 5, C = C, U_obs = U,
    engine = engine, parallel = FALSE, n_chains = 2L, burn = 8L,
    seed = 913L, priors = list(u_dictionary = list(c(.2, .7, 1.8))),
    sampler_control = list(u_block_size = 3L, cross_every = 2L,
                           dictionary_mode = dictionary_mode))
  if (engine == "multivariate") {
    args$ident <- "none"
    args$store_scores <- TRUE
  }
  args
}

testthat::test_that("0.3.0 public operations preserve the sampler and draw contracts", {
  for (engine in c("univariate", "multivariate")) {
    args <- v030_workflow_args(engine)
    args$n_chains <- 4L
    set.seed(781L)
    rng <- mixedgp_rng_state()
    fit <- do.call(fit_eivgp, c(args, list(n_iter = 20L)))
    testthat::expect_identical(mixedgp_rng_state(), rng)
    testthat::expect_identical(fit$sampler_version, "0.3.0")
    testthat::expect_identical(fit$checkpoint$version, 2L)
    testthat::expect_equal(fit$mcmc$samples_sigma2,
                           fit$mcmc$samples_V * (1 - fit$mcmc$samples_r))
    testthat::expect_equal(exp(2 * fit$mcmc$samples_logtheta[, 1L]),
                           fit$mcmc$samples_r / (1 - fit$mcmc$samples_r))
    X <- args$X[1:2, , drop = FALSE]
    C <- args$C[1:2, , drop = FALSE]
    U <- args$U_obs[1:2, , drop = FALSE]
    imputed <- impute_eivgp(fit, rows = 1:2, draw_ids = 1:2)
    testthat::expect_identical(dim(imputed), c(2L, 2L, 1L))
    testthat::expect_equal(as.numeric(imputed[1L, , 1L]), as.numeric(U))
    for (target in c("mean", "surface", "response")) {
      extras <- if (target == "surface") list(new_U = U) else list(new_C = C)
      if (target == "mean") extras$n_latent <- 4L
      pred <- do.call(predict_eivgp, c(list(object = fit, new_X = X,
        target = target, draw_ids = 1:2, seed = 18L, joint = TRUE), extras))
      testthat::expect_identical(dim(pred), c(2L, 2L))
      testthat::expect_true(all(is.finite(pred)))
      testthat::expect_true(attr(pred, "joint"))
    }
    report <- suppressWarnings(diagnose_eivgp(fit, X, C, U, n_latent = 4L))
    testthat::expect_true(all(report$table$min_draws_per_chain == 12L))
    testthat::expect_true(all(report$table$draw_window_complete))
    testthat::expect_true(all(c("V", "r") %in% report$table$parameter))
    testthat::expect_true(any(startsWith(report$table$parameter, "dictionary_occupancy[")))
    testthat::expect_true(any(startsWith(report$table$parameter, "f_conditional_mean[")))
    testthat::expect_identical(mixedgp_rng_state(), rng)
  }
})

testthat::test_that("0.3.0 continuation appends exact same chains across serial and forked execution", {
  for (engine in c("univariate", "multivariate")) {
    for (mode in c("conditional", "marginal")) {
      args <- v030_workflow_args(engine, mode)
      if (mode == "marginal") args$kernel <- "matern"
      first <- do.call(fit_eivgp, c(args, list(n_iter = 16L)))
      path <- tempfile(fileext = ".rds")
      saveRDS(first, path)
      restored <- readRDS(path)
      unlink(path)
      set.seed(883L)
      rng <- mixedgp_rng_state()
      extended <- continue_eivgp(restored, 6L, parallel = FALSE)
      testthat::expect_identical(mixedgp_rng_state(), rng)
      direct <- do.call(fit_eivgp, c(args, list(n_iter = 22L)))
      testthat::expect_identical(extended$mcmc$samples_by_chain, direct$mcmc$samples_by_chain)
      testthat::expect_identical(extended$checkpoint$states, direct$checkpoint$states)
      testthat::expect_equal(extended$mcmc$chain_stats, direct$mcmc$chain_stats, tolerance = 1e-12)
      testthat::expect_identical(extended$data, first$data)
      testthat::expect_equal(extended$diagnostics$summary$saved_per_chain, 14L)
      twice <- continue_eivgp(continue_eivgp(first, 1L, parallel = FALSE), 5L, parallel = FALSE)
      testthat::expect_identical(twice$mcmc$samples_by_chain, extended$mcmc$samples_by_chain)
      parallel_fit <- continue_eivgp(first, 6L, parallel = TRUE, n_cores = 2L)
      testthat::expect_identical(parallel_fit$mcmc$samples_by_chain, extended$mcmc$samples_by_chain)
      testthat::expect_identical(parallel_fit$checkpoint$states, extended$checkpoint$states)
    }
  }
})

testthat::test_that("0.3.0 continuation rejects incompatible or modified checkpoints", {
  args <- v030_workflow_args("univariate")
  first <- do.call(fit_eivgp, c(args, list(n_iter = 16L)))
  old <- first; old$checkpoint <- NULL
  testthat::expect_error(continue_eivgp(old, 5L), "no supported")
  old <- first; old$checkpoint$version <- 1L
  testthat::expect_error(continue_eivgp(old, 5L), "Incompatible sampler checkpoint")
  old <- first; old$checkpoint$sampler_version <- "0.2.1"
  testthat::expect_error(continue_eivgp(old, 5L), "Incompatible sampler checkpoint")
  for (field in c("data", "draws", "pooled_draws", "dictionary", "state")) {
    modified <- first
    if (field == "data") modified$data$y[1L] <- 9
    if (field == "draws") modified$mcmc$samples_by_chain$V[[1L]][1L] <- 9
    if (field == "pooled_draws") modified$mcmc$samples_V[1L] <- 9
    if (field == "dictionary") modified$priors$u_dictionary[[1L]][1L] <- 9
    if (field == "state") modified$checkpoint$states[[1L]]$state$modified <- TRUE
    testthat::expect_error(continue_eivgp(modified, 5L), "modified")
  }
  testthat::expect_error(do.call(fit_eivgp, c(args, list(n_iter = 16L, thin = 2L))), "every post-warm-up")
  testthat::expect_error(do.call(fit_eivgp, c(args, list(n_iter = 16L, noise_strategy = "conditional"))), "Unsupported 0.3.0")
  testthat::expect_error(do.call(fit_eivgp, c(args, list(n_iter = 16L, .resume = list()))), "internal")
})

testthat::test_that("multivariate vector dictionaries support no or complete calibration", {
  n <- 10L
  X <- cbind(exposure = seq(-1, 1, length.out = n))
  U <- cbind(latent1 = sin(seq_len(n)), latent2 = cos(seq_len(n)))
  C <- cbind(ordinal1 = rep(1:2, 5L), ordinal2 = as.integer(U[, 2L] > 0) + 1L)
  dictionary <- rbind(c(.2, 1.2), c(.7, .7), c(1.2, .2))
  for (calibrated in c(FALSE, TRUE)) {
    fit <- suppressWarnings(fit_eivgp(X, sin(X[, 1L]) + U[, 1L], C,
      U_obs = if (calibrated) U else NULL, latent_dim = 2L,
      engine = "multivariate", ident = if (calibrated) "none" else "lower_triangular",
      n_iter = 14L, burn = 6L, n_chains = 2L, parallel = FALSE, seed = 791L,
      priors = list(u_dictionary = dictionary)))
    testthat::expect_identical(dim(fit$mcmc$samples_U), c(16L, n, 2L))
    testthat::expect_identical(ncol(fit$mcmc$samples_J), 1L)
    testthat::expect_equal(unname(exp(fit$mcmc$samples_logtheta[, 3:4])),
      unname(dictionary[fit$mcmc$samples_J[, 1L], , drop = FALSE]))
    imputed <- impute_eivgp(fit, draw_ids = 1L)
    if (calibrated) testthat::expect_equal(unname(imputed[1L, , ]), unname(U)) else {
      testthat::expect_true(all(fit$mcmc$samples_A[, 1L, 1L] > 0))
      testthat::expect_true(all(fit$mcmc$samples_A[, 2L, 2L] > 0))
      testthat::expect_true(all(fit$mcmc$samples_A[, 1L, 2L] == 0))
    }
  }
  selected <- mixedgp_v030_calibration(U, n, 2L, 1:3)
  testthat::expect_true(all(is.na(selected$U_obs[-(1:3), ])))
  testthat::expect_equal(unname(selected$U_obs[1:3, ]), unname(U[1:3, ]))
})
