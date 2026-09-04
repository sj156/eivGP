# Canonical package tests copied into the generated eivGP package.
mixedgp_test_univariate_fit <- function() {
  x <- matrix(c(-1, 0, 1), ncol = 1L, dimnames = list(NULL, "x1"))
  level_maps <- list(severity = c("low", "medium", "high"))
  structure(
    list(
      data = list(
        x = x, x_raw = x, y = c(-0.5, 0, 0.6),
        x_center = 0, x_scale = 1, predictor_names = "x1",
        y_center = 0, y_scale = 1, m = 3L, C_names = "severity",
        C_level_maps = level_maps,
        U_center = 0, U_scale = 1, U_names = "latent_severity",
        calib_idx = integer(0), miss_idx = 1:3,
        latent_scale_anchored = FALSE,
        latent_anchor_rank = 0L, latent_anchor_required_rank = 2L
      ),
      mcmc = list(
        samples_u = rbind(c(-0.8, 0, 0.8), c(-0.7, 0.1, 0.9)),
        samples_tau = rbind(c(-0.25, 0.4), c(-0.2, 0.45)),
        samples_logtheta = rbind(
          log(c(1.2, 0.7, 0.9)), log(c(1.1, 0.8, 1))
        ),
        samples_sigma2 = c(0.2, 0.25)
      ),
      kernel = list(name = "se", matern_nu = NA_real_),
      diagnostics = list(summary = data.frame(
        kernel = "se", matern_nu = NA_real_, n_chains = 2L,
        parallel_backend = "serial", parallel_cores = 1L,
        n_iter = 20L, burn = 10L, thin = 5L,
        saved_per_chain = 1L, total_saved_draws = 2L,
        max_rhat_hyper = 1.01, max_rhat_tau = 1.02,
        median_rhat_missing_u = 1.00, max_rhat_missing_u = 1.01,
        min_ess_key = 75, time_seconds = 0.1
      )),
      interface = list(
        engine = "univariate", ordinal_level_maps = level_maps
      )
    ),
    class = c("eivgp_fit", "list")
  )
}

mixedgp_test_multivariate_fit <- function() {
  x <- matrix(c(-1, 0, 1), ncol = 1L, dimnames = list(NULL, "x1"))
  samples_U <- array(NA_real_, dim = c(2L, 3L, 1L))
  samples_U[1L, , 1L] <- c(-0.7, 0, 0.7)
  samples_U[2L, , 1L] <- c(-0.6, 0.1, 0.8)
  samples_A <- array(NA_real_, dim = c(2L, 2L, 1L))
  samples_A[1L, , 1L] <- c(1, 0.7)
  samples_A[2L, , 1L] <- c(1.1, 0.8)
  level_maps <- list(
    proxy_a = c("low", "medium", "high"),
    proxy_b = c("none", "mild", "severe")
  )
  structure(
    list(
      data = list(
        X = x, X_raw = x, y = c(-0.5, 0, 0.6),
        X_center = 0, X_scale = 1, y_center = 0, y_scale = 1,
        p = 1L, d = 1L, q = 2L, m_vec = c(3L, 3L),
        C_names = names(level_maps), C_level_maps = level_maps,
        C_ord = cbind(
          proxy_a = c(1L, 2L, 3L), proxy_b = c(1L, 2L, 3L)
        ),
        U_center = 0, U_scale = 1, U_names = "latent1",
        calib_idx = integer(0), miss_idx = 1:3, ident = "lower_triangular",
        latent_scale_anchored = FALSE,
        latent_anchor_rank = 0L, latent_anchor_required_rank = 2L
      ),
      mcmc = list(
        samples_U = samples_U, samples_A = samples_A,
        samples_tau = rbind(
          c(-0.5, 0.5, -0.5, 0.5),
          c(-0.4, 0.6, -0.45, 0.55)
        ),
        samples_logtheta = rbind(
          log(c(1.2, 0.7, 0.9)), log(c(1.1, 0.8, 1))
        ),
        samples_sigma2 = c(0.2, 0.25)
      ),
      kernel = list(name = "se", matern_nu = NA_real_),
      diagnostics = list(summary = data.frame(
        sampler_strategy = "blocked", kernel = "se", matern_nu = NA_real_,
        n_chains = 2L, parallel_backend = "serial", parallel_cores = 1L,
        n_iter = 20L, burn = 10L, thin = 5L,
        saved_per_chain = 1L, total_saved_draws = 2L,
        max_rhat_hyper = 1.01, max_rhat_A = 1.02,
        max_rhat_tau = 1.01, median_rhat_missing_U = 1.00,
        max_rhat_missing_U = 1.02, min_ess_key = 80,
        mean_u_ess_accept = 1, mean_global_u_accept = 1,
        time_seconds = 0.2
      )),
      interface = list(
        engine = "multivariate", ordinal_level_maps = level_maps
      )
    ),
    class = c("eivgp_fit", "list")
  )
}

testthat::test_that("robust Study I sampler and mixing plots are available", {
  control <- make_default_control(100L, 95L, preset = "thorough", p_x = 1L)
  testthat::expect_equal(control$theta_joint_every, 2L)
  testthat::expect_equal(control$tau_noncentered_every, 10L)

  testthat::skip_if_not_installed("ggplot2")
  fit <- mixedgp_test_univariate_fit()
  make_chain <- function(offset) {
    draw <- seq_len(30L)
    list(
      logtheta = cbind(
        log(1.2 + 0.01 * sin(draw + offset)),
        log(0.7 + 0.01 * cos(draw + offset)),
        log(0.9 + 0.01 * sin(draw / 2 + offset))
      ),
      sigma2 = 0.2 + 0.005 * cos(draw + offset),
      tau = cbind(
        -0.25 + 0.01 * sin(draw + offset),
        0.4 + 0.01 * cos(draw + offset)
      ),
      u = cbind(
        -0.8 + 0.01 * sin(draw + offset),
        0.01 * cos(draw + offset),
        0.8 + 0.01 * sin(draw / 2 + offset)
      )
    )
  }
  chains <- list(make_chain(0), make_chain(0.7))
  fit$mcmc$samples_by_chain <- list(
    logtheta = lapply(chains, `[[`, "logtheta"),
    sigma2 = lapply(chains, `[[`, "sigma2"),
    tau = lapply(chains, `[[`, "tau"),
    u = lapply(chains, `[[`, "u")
  )
  fit$diagnostics$rhat_u <- data.frame(
    parameter = c("u[1]", "u[2]", "u[3]"),
    global_index = 1:3,
    rhat = c(1.02, 1.01, 1.03)
  )
  plots <- plot_eivgp_mcmc_diagnostics(
    fit, max_draws = 20L, max_lag = 10L, max_latent = 2L
  )
  testthat::expect_s3_class(plots$trace, "ggplot")
  testthat::expect_s3_class(plots$autocorrelation, "ggplot")
  testthat::expect_s3_class(plots$rank, "ggplot")
  testthat::expect_true(nrow(plots$draws) > 0L)
})

testthat::test_that("noncentered Study I threshold move preserves constraints", {
  set.seed(812L)
  x <- matrix(seq(-1, 1, length.out = 6L), ncol = 1L)
  y <- c(-0.6, -0.4, -0.1, 0.1, 0.4, 0.7)
  c_ord <- rep(1:3, each = 2L)
  tau <- c(-0.35, 0.35)
  u <- c(-0.9, -0.6, -0.2, 0.2, 0.6, 0.9)
  calib_idx <- c(1L, 6L)
  u_obs <- rep(NA_real_, 6L)
  u_obs[calib_idx] <- u[calib_idx]
  theta_spec <- make_theta_spec_multix(1L)
  answer <- update_tau_u_noncentered_1d(
    y = y, u = u,
    Dx_list = list(pairwise_sqdist(x)),
    c_ord = c_ord, tau = tau,
    logtheta = log(c(2, 0.7, 0.8)), sigma2_eps = 0.2,
    theta_spec = theta_spec, miss_idx = setdiff(1:6, calib_idx),
    calib_idx = calib_idx, u_obs = u_obs
  )
  testthat::expect_true(check_constraints_1d(answer$u, c_ord, answer$tau))
  testthat::expect_equal(answer$u[calib_idx], u_obs[calib_idx])
  testthat::expect_true(is.finite(answer$n_eval) && answer$n_eval > 0L)
})

testthat::test_that("parallel maps are reproducible", {
  f <- function(i) c(i = i, u = stats::runif(1))
  serial <- mixedgp_parallel_lapply(
    as.list(1:3), f, n_cores = 1, seeds = 101:103
  )
  parallel <- mixedgp_parallel_lapply(
    as.list(1:3), f, n_cores = 2, seeds = 101:103
  )
  testthat::expect_identical(serial, parallel)

  serial_m <- mixedgp_parallel_mapply(
    function(a, b) c(sum = a + b, u = stats::runif(1)),
    a = 1:3, b = 3:1, n_cores = 1, seeds = 201:203
  )
  parallel_m <- mixedgp_parallel_mapply(
    function(a, b) c(sum = a + b, u = stats::runif(1)),
    a = 1:3, b = 3:1, n_cores = 2, seeds = 201:203
  )
  testthat::expect_identical(serial_m, parallel_m)
})

testthat::test_that("adaptive MCMC projects bounded continuation lengths", {
  control <- mixedgp_adaptive_mcmc_control(
    initial_draws = 5000L,
    target_ess = 200,
    max_draws = 15000L,
    extension_draws = 1000L
  )
  extend <- mixedgp_adaptive_next_draws(
    current_draws = 5000L,
    observed_ess = 25,
    control = control,
    n_chains = 4L
  )
  stop_now <- mixedgp_adaptive_next_draws(
    current_draws = 5000L,
    observed_ess = 55,
    control = control,
    n_chains = 4L
  )
  testthat::expect_equal(extend$next_draws, 10000L)
  testthat::expect_equal(stop_now$next_draws, 5000L)
  testthat::expect_equal(stop_now$reason, "ess_target_met")
})

testthat::test_that("Study I adaptive MCMC continues without restarting", {
  set.seed(4113)
  n <- 18L
  x <- seq(-1, 1, length.out = n)
  u <- stats::rnorm(n)
  c_ord <- cut(u, c(-Inf, -0.4, 0.4, Inf), labels = FALSE)
  y <- sin(x) + u + stats::rnorm(n, sd = 0.1)
  fit <- fit_eivgp_1d(
    x, y, c_ord, u_true = u, calib_idx = 1:6, m = 3L,
    n_iter = 25L, burn = 5L, thin = 1L, n_chains = 1L,
    preset = "fast", parallel_chains = FALSE,
    adaptive_control = list(
      initial_draws = 10L, target_ess = 100,
      max_draws = 20L, extension_draws = 5L
    )
  )
  testthat::expect_equal(fit$mcmc$chain_stats$saved, 20L)
  testthat::expect_equal(fit$mcmc$chain_stats$iterations_run, 25L)
  testthat::expect_true(fit$mcmc$chain_stats$adaptive_extensions >= 1L)
  testthat::expect_true(nrow(fit$mcmc$adaptive_schedule) >= 2L)
  testthat::expect_true(length(fit$mcmc$chain_final[[1L]]$rng_state) > 1L)
})

testthat::test_that("Study II uses the shared adaptive MCMC protocol", {
  set.seed(4114)
  n <- 18L
  X <- matrix(seq(-1, 1, length.out = n), ncol = 1L)
  U <- cbind(stats::rnorm(n), stats::rnorm(n))
  C <- cbind(
    cut(U[, 1L], c(-Inf, -0.3, 0.3, Inf), labels = FALSE),
    cut(U[, 2L], c(-Inf, -0.3, 0.3, Inf), labels = FALSE)
  )
  y <- X[, 1L] + U[, 1L] - 0.5 * U[, 2L] + stats::rnorm(n, sd = 0.2)
  fit <- fit_eivgp_ordprobit_fb(
    X, y, C, U_obs = U, calib_idx = 1:8,
    d = 2L, m_vec = c(3L, 3L),
    n_iter = 15L, burn = 5L, thin = 1L, n_chains = 1L,
    preset = "fast", sampler_strategy = "legacy",
    parallel_chains = FALSE, store_scores = FALSE,
    adaptive_control = list(
      initial_draws = 6L, target_ess = 100,
      max_draws = 10L, extension_draws = 2L
    )
  )
  testthat::expect_equal(fit$mcmc$chain_stats$saved, 10L)
  testthat::expect_equal(fit$mcmc$chain_stats$iterations_run, 15L)
  testthat::expect_true(fit$mcmc$chain_stats$adaptive_extensions >= 1L)
  testthat::expect_true(nrow(fit$mcmc$adaptive_schedule) >= 2L)
  testthat::expect_true(length(fit$mcmc$chain_final[[1L]]$rng_state) > 1L)
})

testthat::test_that("isolated simulation engines retain adaptive RNG state", {
  testthat::skip_if(
    isTRUE(getOption("knitr.in.progress")),
    "requires installed system.file() semantics; exercised by R CMD check"
  )
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  original_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", original_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (had_seed) rm(".Random.seed", envir = .GlobalEnv)

  source_dir <- system.file("simulation-source", package = "eivGP")
  engine <- mixedgp_simulation_engine(source_dir)
  x <- seq(-1, 1, length.out = 12L)
  u <- seq(-0.8, 0.8, length.out = 12L)
  c_ord <- cut(u, c(-Inf, -0.3, 0.3, Inf), labels = FALSE)
  fit <- engine$fit_eivgp_1d(
    x, sin(x) + u, c_ord, u_true = u, calib_idx = 1:4, m = 3L,
    n_iter = 15L, burn = 5L, thin = 1L, n_chains = 1L,
    preset = "fast", parallel_chains = FALSE,
    adaptive_control = list(
      initial_draws = 5L, target_ess = 100,
      max_draws = 10L, extension_draws = 5L
    )
  )
  testthat::expect_true(length(fit$mcmc$chain_final[[1L]]$rng_state) > 1L)
})

testthat::test_that("seeded parallel maps preserve the caller RNG stream", {
  set.seed(4101)
  expected_lapply_continuation <- stats::runif(4)
  set.seed(4101)
  invisible(mixedgp_parallel_lapply(
    as.list(1:3), function(i) stats::runif(i),
    n_cores = 2, seeds = 501:503
  ))
  testthat::expect_identical(stats::runif(4), expected_lapply_continuation)

  set.seed(4102)
  expected_mapply_continuation <- stats::runif(4)
  set.seed(4102)
  invisible(mixedgp_parallel_mapply(
    function(a, b) c(a + b, stats::runif(1)),
    a = 1:3, b = 3:1, n_cores = 2, seeds = 601:603
  ))
  testthat::expect_identical(stats::runif(4), expected_mapply_continuation)
})

testthat::test_that("seeded public operations preserve the caller RNG stream", {
  univariate <- mixedgp_test_univariate_fit()
  new_X <- matrix(
    c(-0.25, 0.25), ncol = 1L, dimnames = list(NULL, "x1")
  )
  new_C <- data.frame(severity = ordered(
    c("high", "low"), levels = c("low", "high")
  ))

  set.seed(4110)
  expected_predict <- stats::runif(4)
  set.seed(4110)
  invisible(predict_eivgp(
    univariate, new_X, new_C, target = "mean", draw_ids = 1L,
    n_latent = 4L, include_process_uncertainty = FALSE, seed = 9001L
  ))
  testthat::expect_identical(stats::runif(4), expected_predict)

  set.seed(4111)
  expected_impute <- stats::runif(4)
  set.seed(4111)
  invisible(impute_eivgp(
    univariate, new_C = new_C, draw_ids = 1L, seed = 9002L
  ))
  testthat::expect_identical(stats::runif(4), expected_impute)

  X <- matrix(seq(-1, 1, length.out = 8L), ncol = 1L)
  u <- seq(-1.2, 1.2, length.out = 8L)
  C <- ordered(
    ifelse(u < 0, "low", "high"), levels = c("low", "high")
  )
  U_obs <- rep(NA_real_, 8L)
  U_obs[c(1L, 4L, 8L)] <- u[c(1L, 4L, 8L)]
  set.seed(4112)
  expected_fit <- stats::runif(4)
  set.seed(4112)
  invisible(fit_eivgp(
    X, sin(u), C, U_obs, engine = "univariate", parallel = FALSE,
    n_iter = 20L, burn = 10L, thin = 2L, n_chains = 1L,
    preset = "fast", seed = 9003L, verbose = FALSE
  ))
  testthat::expect_identical(stats::runif(4), expected_fit)

  adapter <- function(..., seed) {
    set.seed(seed)
    stats::runif(2L)
  }
  testthat::local_mocked_bindings(
    mixedgp_adapter_ucgp = adapter,
    mixedgp_adapter_lvgp = adapter,
    mixedgp_adapter_ezgp = adapter,
    .package = "eivGP"
  )
  set.seed(4113)
  expected_competitors <- stats::runif(4)
  set.seed(4113)
  invisible(fit_ucgp(1, 1, 1, 1, 1, seed = 9004L))
  invisible(fit_lvgp(1, 1, 1, 1, 1, seed = 9005L))
  invisible(fit_ezgp(1, 1, 1, 1, 1, seed = 9006L))
  testthat::expect_identical(stats::runif(4), expected_competitors)
})

testthat::test_that("public controls and input schemas fail closed", {
  X <- matrix(c(-1, 0, 1), ncol = 1L)
  y <- c(-0.5, 0, 0.5)
  C1 <- ordered(c("low", "middle", "high"))
  testthat::expect_error(
    fit_eivgp(X, y, C1, U_obs = factor(c(10, 20, 30)),
              engine = "univariate"),
    "must be numeric"
  )
  testthat::expect_error(
    fit_eivgp(X, y, C1, engine = "univariate", parallel = 1),
    "parallel must be TRUE or FALSE"
  )
  testthat::expect_error(
    fit_eivgp(X, y, C1, engine = "univariate", m_vec = 3.5),
    "integer-valued"
  )

  C2 <- data.frame(
    proxy_a = ordered(c("low", "middle", "high")),
    proxy_b = ordered(c("none", "mild", "severe"))
  )
  U_rank_deficient <- matrix(
    c(0, 0, 1, 1, NA, NA), nrow = 3L, byrow = TRUE
  )
  testthat::expect_error(
    fit_eivgp(
      X, y, C2, U_obs = U_rank_deficient, latent_dim = 2L,
      engine = "multivariate", ident = "none"
    ),
    "full affine rank"
  )
  testthat::expect_error(
    fit_eivgp(X, y, C2, latent_dim = 1.5, engine = "multivariate"),
    "integer-valued"
  )

  duplicate_X <- matrix(1:6, nrow = 3L)
  colnames(duplicate_X) <- c("x", "x")
  testthat::expect_error(
    as_numeric_matrix_strict(duplicate_X, "X"), "complete and unique"
  )
  duplicate_C <- data.frame(
    ordered(c("a", "b", "a")), ordered(c("c", "d", "c")),
    check.names = FALSE
  )
  names(duplicate_C) <- c("proxy", "proxy")
  testthat::expect_error(
    prepare_ordinal_matrix(duplicate_C, name = "C"), "complete and unique"
  )
  testthat::expect_error(
    mixedgp_do_call(identity, list(x = 1), list(2)), "must be named"
  )
})

testthat::test_that("package loading has no reporting or thread side effects", {
  package_path <- getNamespaceInfo(asNamespace("eivGP"), "path")
  check_script <- tempfile("eivGP-load-side-effects-", fileext = ".R")
  on.exit(unlink(check_script), add = TRUE)
  writeLines(
    c(
      "args <- commandArgs(trailingOnly = TRUE)",
      "package_path <- normalizePath(args[[1L]], winslash = '/', mustWork = TRUE)",
      "thread_vars <- c('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',",
      "  'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')",
      "reporting <- c('ggplot2', 'dplyr', 'tidyr')",
      "have_ggplot <- requireNamespace('ggplot2', quietly = TRUE)",
      "theme_before <- if (have_ggplot) ggplot2::theme_get() else NULL",
      "loaded_before <- reporting %in% loadedNamespaces()",
      "names(loaded_before) <- reporting",
      "search_before <- search()",
      "env_before <- Sys.getenv(thread_vars, unset = NA_character_)",
      "source_modules <- c(",
      "  'core_parallel.R', 'core_numerics.R', 'model_univariate.R',",
      "  'model_multivariate.R', 'competitors.R',",
      "  'reproduction_data.R', 'model_api.R',",
      "  'reproduction_compat.R')",
      "source_paths <- file.path(package_path, 'R', source_modules)",
      "if (all(file.exists(source_paths))) {",
      "  package_env <- new.env(parent = baseenv())",
      "  for (path in source_paths) sys.source(path, envir = package_env)",
      "} else {",
      "  suppressPackageStartupMessages(library(",
      "    'eivGP', character.only = TRUE, lib.loc = dirname(package_path)",
      "  ))",
      "}",
      "search_after <- search()",
      "allowed <- unique(c(search_before, 'package:eivGP'))",
      "stopifnot(length(setdiff(search_after, allowed)) == 0L)",
      "stopifnot(identical(",
      "  Sys.getenv(thread_vars, unset = NA_character_), env_before",
      "))",
      "loaded_after <- reporting %in% loadedNamespaces()",
      "names(loaded_after) <- reporting",
      "stopifnot(identical(loaded_after, loaded_before))",
      "if (have_ggplot) stopifnot(identical(ggplot2::theme_get(), theme_before))"
    ),
    check_script
  )
  output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(check_script), shQuote(package_path)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  testthat::expect_true(
    is.null(status) || identical(as.integer(status), 0L),
    info = paste(output, collapse = "\n")
  )
})

testthat::test_that("tail-stable truncated normals respect extreme support", {
  lower <- rep(c(8, -Inf, 40, -41), each = 64L)
  upper <- rep(c(Inf, -8, 41, -40), each = 64L)
  set.seed(4103)
  z <- rtruncnorm_vec(0, 1, lower, upper)
  testthat::expect_true(all(is.finite(z)))
  testthat::expect_true(all(z >= lower & z <= upper))
})

testthat::test_that("stored synthetic data round trips", {
  dat <- make_study1_synthetic_dataset(
    rep_id = 1, n = 30, n_test = 10, m = 6, min_class_count = 1,
    calib_grid = c(0L, 5L, 10L, 20L, 30L)
  )
  testthat::expect_s3_class(dat, "mixedgp_synthetic_dataset")
  testthat::expect_equal(nrow(dat$data$train), 30)
})

testthat::test_that("frozen-data filename collisions fail closed", {
  data_dir <- tempfile("eivGP-frozen-")
  dir.create(data_dir)
  on.exit(unlink(data_dir, recursive = TRUE), add = TRUE)

  original <- make_study1_synthetic_dataset(
    rep_id = 2, n = 30, n_test = 10, m = 6, min_class_count = 1,
    heterogeneity_eta = 0, calib_grid = c(0L, 5L, 10L, 20L, 30L)
  )
  collision <- make_study1_synthetic_dataset(
    rep_id = 2, n = 30, n_test = 10, m = 6, min_class_count = 1,
    heterogeneity_eta = 1, calib_grid = c(0L, 5L, 10L, 20L, 30L)
  )
  manifest <- store_mixedgp_synthetic_dataset(original, data_dir)
  saveRDS(manifest, file.path(data_dir, "manifest.rds"), version = 3L)
  testthat::expect_error(
    store_mixedgp_synthetic_dataset(collision, data_dir),
    "different synthetic dataset"
  )
  reloaded <- load_mixedgp_synthetic_dataset(
    file.path(data_dir, manifest$file)
  )
  testthat::expect_identical(reloaded, original)
  strict <- load_mixedgp_synthetic_dataset_strict(
    file.path(data_dir, manifest$file),
    expected = list(
      study = "study1", scenario = "active", rep_id = 2L,
      n = 30L, n_test = 10L, m = 6L,
      calib_grid = c(0L, 5L, 10L, 20L, 30L)
    )
  )
  testthat::expect_identical(strict$data, original$data)
  testthat::expect_identical(strict$design, original$design)
  testthat::expect_identical(attr(strict, "manifest_md5"), manifest$md5)
  testthat::expect_error(
    load_mixedgp_synthetic_dataset_strict(
      file.path(data_dir, manifest$file),
      expected = list(study = "study1", n = 31L)
    ),
    "expected 31"
  )
  testthat::expect_error(
    make_nested_calibration_sets(30L, c(0L, 31L), seed = 1L),
    "between zero and n"
  )
})

testthat::test_that("synthetic manifests merge atomically and fail on worker errors", {
  data_dir <- tempfile("eivGP-incremental-")
  dir.create(data_dir)
  on.exit(unlink(data_dir, recursive = TRUE), add = TRUE)

  first <- generate_study2_synthetic_datasets(
    n_rep = 1L, scenarios = "primary", directory = data_dir,
    n_cores = 1L, n = 12L, n_test = 5L, q = 2L, m = 3L,
    calib_grid = c(0L, 3L)
  )
  merged <- generate_study2_synthetic_datasets(
    n_rep = 1L, scenarios = "high_uncertainty", directory = data_dir,
    n_cores = 1L, n = 12L, n_test = 5L, q = 2L, m = 3L,
    calib_grid = c(0L, 3L)
  )
  testthat::expect_equal(nrow(first), 1L)
  testthat::expect_equal(nrow(merged), 2L)
  testthat::expect_identical(readRDS(file.path(data_dir, "manifest.rds")), merged)
  testthat::expect_setequal(
    utils::read.csv(file.path(data_dir, "manifest.csv"))$file, merged$file
  )
  testthat::expect_true(all(nzchar(merged$generator_tag)))
  testthat::expect_true(all(grepl("^[0-9a-f]{32}$", merged$design_fingerprint)))

  failed <- structure("worker boom", class = "try-error")
  testthat::expect_error(
    mixedgp_stop_on_worker_errors(list(failed), "generation"),
    "no manifest was written"
  )
  inconsistent <- merged[1L, , drop = FALSE]
  inconsistent$md5 <- paste(rep("0", 32L), collapse = "")
  testthat::expect_error(
    mixedgp_merge_manifests(merged, inconsistent), "metadata disagree"
  )
  testthat::expect_error(
    generate_study2_synthetic_datasets(
      n_rep = 0L, scenarios = "primary", directory = data_dir
    ),
    "n_rep"
  )
  testthat::expect_error(
    generate_study2_synthetic_datasets(
      n_rep = 1L, scenarios = character(), directory = data_dir
    ),
    "one or more"
  )
})

testthat::test_that("Study II frozen truth and fingerprints are validated", {
  artifact <- make_study2_synthetic_dataset(
    rep_id = 1L, scenario = "primary", n = 12L, n_test = 5L,
    q = 2L, m = 3L, calib_grid = c(0L, 3L)
  )
  testthat::expect_invisible(validate_mixedgp_synthetic_dataset(artifact))

  bad_U <- artifact
  bad_U$data$train$U[1L, 1L] <- NA_real_
  testthat::expect_error(
    validate_mixedgp_synthetic_dataset(bad_U), "Study II data"
  )
  bad_truth <- artifact
  bad_truth$data$test$f <- bad_truth$data$test$f[-1L]
  testthat::expect_error(
    validate_mixedgp_synthetic_dataset(bad_truth), "Study II data"
  )
  stale <- artifact
  stale$design$generator_tag <- "stale-generator"
  testthat::expect_error(
    validate_mixedgp_synthetic_dataset(stale), "generator tag"
  )
})

testthat::test_that("Study II oracle pools are tied to their generating model", {
  pars <- make_study2_true_params("primary", q = 4L, m = 4L)
  pool <- make_oracle_pool_2d(pars, n_pool = 200L, seed = 41025L)
  testthat::expect_invisible(validate_oracle_pool_2d(pool, pars))
  changed <- pars
  changed$A[1L, 1L] <- changed$A[1L, 1L] + 0.01
  testthat::expect_error(
    validate_oracle_pool_2d(pool, changed),
    "provenance"
  )
})

testthat::test_that("ordinal factor maps survive dropped prediction levels", {
  training <- data.frame(severity = ordered(
    c("low", "high", "medium"),
    levels = c("low", "medium", "high")
  ))
  encoded <- prepare_ordinal_matrix(training, name = "C")
  prediction <- data.frame(severity = ordered(
    c("high", "low"), levels = c("low", "high")
  ))
  encoded_new <- prepare_ordinal_matrix(
    prediction,
    m_vec = encoded$m_vec,
    level_maps = encoded$level_maps,
    expected_names = encoded$column_names,
    name = "new_C"
  )
  testthat::expect_identical(unname(encoded_new$C[, 1L]), c(3L, 1L))
  testthat::expect_error(
    prepare_ordinal_matrix(
      data.frame(severity = ordered("unseen")),
      m_vec = encoded$m_vec,
      level_maps = encoded$level_maps,
      expected_names = encoded$column_names,
      name = "new_C"
    ),
    "unknown level"
  )
})

testthat::test_that("multivariate latent dimension is never silently guessed", {
  X <- matrix(c(-1, 0, 1), ncol = 1L)
  y <- c(-0.5, 0, 0.5)
  C <- data.frame(
    proxy_a = ordered(c("low", "medium", "high")),
    proxy_b = ordered(c("none", "mild", "severe"))
  )
  testthat::expect_error(
    fit_eivgp(X, y, C, engine = "multivariate", U_obs = NULL),
    "latent_dim must be supplied explicitly"
  )
  testthat::expect_error(
    fit_eivgp(
      X, y, C, engine = "multivariate",
      U_obs = matrix(NA_real_, nrow = 3L, ncol = 2L)
    ),
    "U_obs contains no calibration values"
  )

  rank_deficient <- mixedgp_latent_anchor_status(
    matrix(c(0, 0, 1, 1), ncol = 2L, byrow = TRUE), 1:2, d = 2L
  )
  anchored <- mixedgp_latent_anchor_status(
    matrix(c(0, 0, 1, 0, 0, 1), ncol = 2L, byrow = TRUE), 1:3, d = 2L
  )
  testthat::expect_false(rank_deficient$anchored)
  testthat::expect_identical(rank_deficient$affine_rank, 2L)
  testthat::expect_true(anchored$anchored)
  testthat::expect_identical(anchored$affine_rank, 3L)
})

testthat::test_that("duplicate-point GP predictions remain positive semidefinite", {
  prediction <- gp_predict_draw_general(
    X_train = matrix(c(-1, 0, 1), ncol = 1L),
    U_train = matrix(c(-0.5, 0, 0.5), ncol = 1L),
    y_train = c(-0.7, 0, 0.8),
    X_star = matrix(c(0.25, 0.25), ncol = 1L),
    U_star = matrix(c(0.1, 0.1), ncol = 1L),
    logtheta = log(c(1, 1, 1)),
    sigma2_eps = 0.1,
    return_cov = TRUE
  )
  testthat::expect_equal(prediction$mean[1L], prediction$mean[2L])
  testthat::expect_true(
    min(eigen(prediction$cov, symmetric = TRUE, only.values = TRUE)$values) >=
      -1e-10
  )
  set.seed(4104)
  draws <- rmvnorm_psd(32L, prediction$mean, prediction$cov)
  testthat::expect_true(all(is.finite(draws)))
  testthat::expect_equal(draws[, 1L], draws[, 2L], tolerance = 1e-7)
})

testthat::test_that("sampler safety caps fail explicitly", {
  testthat::expect_error(
    bounded_slice_update_1d(
      x0 = 0, logf = function(x) stats::dnorm(x, log = TRUE),
      max_steps_out = 0L, max_iter = 0L, fail_on_limit = TRUE
    ),
    "exceeded max_iter"
  )

  X <- matrix(c(-0.5, 0.5), ncol = 1L)
  U <- matrix(c(-0.5, 0.5), ncol = 1L)
  C <- matrix(c(1L, 2L), ncol = 1L)
  prior <- make_gp_prior(p = 1L, d = 1L)
  testthat::expect_error(
    update_U_theta_ess_integrated_general(
      y = c(-0.2, 0.2), X = X, C = C, U_curr = U,
      logtheta_curr = prior$mean, gp_prior = prior,
      block_idx = 1:2, reference = "prior",
      A = matrix(1, nrow = 1L), tau = list(0), max_try = 0L
    ),
    "exceeded max_try"
  )
})

testthat::test_that("both engines expose mean inference and labeled imputation", {
  univariate <- mixedgp_test_univariate_fit()
  multivariate <- mixedgp_test_multivariate_fit()
  new_X <- matrix(
    c(-0.25, 0.25), ncol = 1L, dimnames = list(NULL, "x1")
  )
  new_C_1d <- data.frame(severity = ordered(
    c("high", "low"), levels = c("low", "high")
  ))
  new_C_multi <- data.frame(
    proxy_a = ordered(c("medium", "low"), levels = c("low", "medium")),
    proxy_b = ordered(c("mild", "none"), levels = c("none", "mild"))
  )

  mean_1d <- predict_eivgp(
    univariate, new_X, new_C_1d, target = "mean", draw_ids = 1:2,
    n_latent = 4L, include_process_uncertainty = FALSE, seed = 4201
  )
  mean_multi <- predict_eivgp(
    multivariate, new_X, new_C_multi, target = "mean", draw_ids = 1:2,
    n_latent = 4L, include_process_uncertainty = FALSE,
    seed = 4202
  )
  testthat::expect_identical(dim(mean_1d), c(2L, 2L))
  testthat::expect_identical(dim(mean_multi), c(2L, 2L))
  testthat::expect_true(all(is.finite(mean_1d)))
  testthat::expect_true(all(is.finite(mean_multi)))
  testthat::expect_identical(attr(mean_1d, "target"), "m(x,c)")
  testthat::expect_identical(attr(mean_multi, "target"), "m(x,c)")

  train_1d <- impute_eivgp(univariate, rows = 1:2, draw_ids = 1:2)
  train_multi <- impute_eivgp(multivariate, rows = 1:2, draw_ids = 1:2)
  prospective_1d <- impute_eivgp(
    univariate, new_C = new_C_1d, draw_ids = 1L,
    n_per_draw = 2L, seed = 4203
  )
  prospective_multi <- impute_eivgp(
    multivariate, new_C = new_C_multi, draw_ids = 1L,
    n_per_draw = 2L, seed = 4204
  )
  testthat::expect_identical(dim(train_1d), c(2L, 2L, 1L))
  testthat::expect_identical(dim(train_multi), c(2L, 2L, 1L))
  testthat::expect_identical(attr(train_1d, "source"), "training")
  testthat::expect_identical(attr(train_multi, "source"), "training")
  testthat::expect_identical(attr(prospective_1d, "source"), "prospective")
  testthat::expect_identical(attr(prospective_multi, "source"), "prospective")
  testthat::expect_identical(attr(prospective_1d, "scale"), "working")
  testthat::expect_identical(attr(prospective_multi, "scale"), "working")
  testthat::expect_false(attr(prospective_1d, "latent_scale_anchored"))
  testthat::expect_false(attr(prospective_multi, "latent_scale_anchored"))
  testthat::expect_identical(
    dimnames(prospective_1d)[[3L]], "latent_severity"
  )
  testthat::expect_identical(dimnames(prospective_multi)[[3L]], "latent1")
  testthat::expect_true(all(prospective_1d[, 1L, 1L] > 0.4))
  testthat::expect_true(all(prospective_1d[, 2L, 1L] < -0.25))
  testthat::expect_identical(
    attr(prospective_multi, "latent_sampler"), "minimax_tilting"
  )
})

testthat::test_that("fit S3 methods report identification and delegate prediction", {
  univariate <- mixedgp_test_univariate_fit()
  multivariate <- mixedgp_test_multivariate_fit()
  summary_1d <- summary.eivgp_fit(univariate)
  summary_multi <- summary.eivgp_fit(multivariate)

  testthat::expect_s3_class(summary_1d, "summary_eivgp_fit")
  testthat::expect_s3_class(summary_multi, "summary_eivgp_fit")
  testthat::expect_identical(summary_1d$engine, "univariate")
  testthat::expect_identical(summary_multi$engine, "multivariate")
  testthat::expect_identical(
    summary_1d$dimensions, c(n = 3L, p = 1L, q = 1L, d = 1L)
  )
  testthat::expect_identical(
    summary_multi$dimensions, c(n = 3L, p = 1L, q = 2L, d = 1L)
  )
  testthat::expect_identical(summary_1d$calibration$status, "none")
  testthat::expect_identical(summary_multi$calibration$status, "none")
  testthat::expect_false(summary_1d$calibration$anchored)
  testthat::expect_false(summary_multi$calibration$anchored)
  testthat::expect_match(
    summary_1d$calibration$interpretation, "No calibration data"
  )
  testthat::expect_equal(summary_1d$convergence$max_rhat_tau, 1.02)
  testthat::expect_equal(summary_multi$convergence$min_ess_key, 80)
  testthat::expect_false("mean_u_ess_accept" %in% names(summary_multi$convergence))
  printed_fit <- paste(
    utils::capture.output(print.eivgp_fit(univariate)), collapse = "\n"
  )
  printed_summary <- paste(
    utils::capture.output(print.summary_eivgp_fit(summary_multi)),
    collapse = "\n"
  )
  if (!grepl("working latent scale", printed_fit, fixed = TRUE)) {
    stop("The fitted-model print method omitted its latent-scale label.")
  }
  if (!grepl("No calibration data", printed_summary, fixed = TRUE)) {
    stop("The summary print method omitted its calibration interpretation.")
  }

  anchored_1d <- univariate
  anchored_1d$data$calib_idx <- 1:2
  anchored_1d$data$latent_scale_anchored <- TRUE
  anchored_1d$data$latent_anchor_rank <- 2L
  anchored_summary <- summary.eivgp_fit(anchored_1d)
  testthat::expect_identical(
    anchored_summary$calibration$status, "affine_anchored"
  )
  testthat::expect_identical(anchored_summary$calibration$output_scale, "raw")

  rank_deficient_multi <- multivariate
  rank_deficient_multi$data$calib_idx <- 1L
  rank_deficient_multi$data$latent_anchor_rank <- 1L
  rank_deficient_summary <- summary.eivgp_fit(rank_deficient_multi)
  testthat::expect_identical(
    rank_deficient_summary$calibration$status, "rank_deficient"
  )
  testthat::expect_match(
    rank_deficient_summary$calibration$interpretation,
    "affine rank 1 of 2"
  )

  new_X <- matrix(
    c(-0.25, 0.25), ncol = 1L, dimnames = list(NULL, "x1")
  )
  new_C <- data.frame(severity = ordered(
    c("high", "low"), levels = c("low", "high")
  ))
  via_s3 <- predict.eivgp_fit(
    univariate, new_X = new_X, new_C = new_C, target = "mean",
    draw_ids = 1L, n_latent = 4L,
    include_process_uncertainty = FALSE, seed = 4210
  )
  via_public <- predict_eivgp(
    univariate, new_X = new_X, new_C = new_C, target = "mean",
    draw_ids = 1L, n_latent = 4L,
    include_process_uncertainty = FALSE, seed = 4210
  )
  testthat::expect_identical(via_s3, via_public)
})

testthat::test_that("exact new U and unanchored scale semantics are explicit", {
  univariate <- mixedgp_test_univariate_fit()
  multivariate <- mixedgp_test_multivariate_fit()
  new_X <- matrix(
    c(-0.25, 0.25), ncol = 1L, dimnames = list(NULL, "x1")
  )
  new_U_multi <- matrix(
    c(-0.1, 0.2), ncol = 1L,
    dimnames = list(NULL, "latent1")
  )
  response_1d <- predict_eivgp(
    univariate, new_X, new_U = c(-0.1, 0.2),
    target = "response", draw_ids = 1:2, seed = 4205
  )
  response_multi <- predict_eivgp(
    multivariate, new_X, new_U = new_U_multi,
    target = "response", draw_ids = 1:2, seed = 4206
  )
  testthat::expect_true(all(is.finite(response_1d)))
  testthat::expect_true(all(is.finite(response_multi)))
  testthat::expect_identical(dim(response_1d), c(2L, 2L))
  testthat::expect_identical(dim(response_multi), c(2L, 2L))
  testthat::expect_identical(attr(response_1d, "latent_input_scale"), "working")
  testthat::expect_identical(
    attr(response_multi, "latent_input_scale"), "working"
  )
  testthat::expect_true(all(is.na(
    attr(response_multi, "rejection_acceptance")
  )))
  testthat::expect_false(attr(response_1d, "latent_scale_anchored"))
  testthat::expect_false(attr(response_multi, "latent_scale_anchored"))

  imputed_1d <- impute_eivgp(univariate, rows = 1:2, scale = "auto")
  imputed_multi <- impute_eivgp(multivariate, rows = 1:2, scale = "auto")
  testthat::expect_identical(attr(imputed_1d, "scale"), "working")
  testthat::expect_identical(attr(imputed_multi, "scale"), "working")
  testthat::expect_false(attr(imputed_1d, "latent_scale_anchored"))
  testthat::expect_false(attr(imputed_multi, "latent_scale_anchored"))
  testthat::expect_error(
    impute_eivgp(univariate, scale = "raw"), "requires calibration data"
  )
  testthat::expect_error(
    impute_eivgp(multivariate, scale = "raw"), "requires calibration data"
  )
  testthat::expect_error(
    predict_eivgp(
      univariate, new_X, new_U = c(-0.1, 0.2),
      new_U_scale = "raw", target = "response"
    ),
    "requires calibration data"
  )
  testthat::expect_error(
    predict_eivgp(
      multivariate, new_X, new_U = new_U_multi,
      new_U_scale = "raw", target = "response"
    ),
    "requires calibration data"
  )
})

testthat::test_that("real-data conveniences delegate to the stable API", {
  X <- matrix(seq(-1, 1, length.out = 4L), ncol = 1L)
  compact_U <- matrix(c(-1, 1), ncol = 1L)
  expanded <- mixedgp_expand_real_data_calibration(
    X, compact_U, calib_idx = c(1L, 4L), d = 1L
  )
  testthat::expect_identical(expanded$d, 1L)
  testthat::expect_equal(
    as.numeric(expanded$U_obs), c(-1, NA_real_, NA_real_, 1)
  )
  stable_fit_stub <- function(...) list(...)
  testthat::local_mocked_bindings(
    fit_eivgp = stable_fit_stub,
    .package = "eivGP"
  )
  delegated_fit <- fit_eivgp_real_data(
    X = X, y = seq_len(4L),
    C = data.frame(proxy = ordered(c("a", "b", "a", "b"))),
    U_obs = compact_U, calib_idx = c(1L, 4L), d = 1L,
    parallel = FALSE
  )
  testthat::expect_identical(delegated_fit$engine, "multivariate")
  testthat::expect_equal(
    as.numeric(delegated_fit$U_obs), c(-1, NA_real_, NA_real_, 1)
  )

  multivariate <- mixedgp_test_multivariate_fit()
  new_X <- matrix(
    c(-0.25, 0.25), ncol = 1L, dimnames = list(NULL, "x1")
  )
  new_C <- data.frame(
    proxy_a = ordered(c("medium", "low"), levels = c("low", "medium")),
    proxy_b = ordered(c("mild", "none"), levels = c("none", "mild"))
  )

  set.seed(4211)
  expected_continuation <- stats::runif(4)
  set.seed(4211)
  delegated_mean <- predict_eivgp_m_given_xc(
    multivariate, new_X, new_C, draw_ids = 1L, n_latent = 4L,
    include_process_uncertainty = FALSE, seed = 9211L
  )
  testthat::expect_identical(stats::runif(4), expected_continuation)
  testthat::expect_identical(attr(delegated_mean, "target"), "m(x,c)")
  testthat::expect_error(
    predict_eivgp_f_given_xu(
      multivariate, new_X, matrix(c(-0.1, 0.2), ncol = 1L),
      U_new_scale = "raw", draw_ids = 1L, seed = 9212L
    ),
    "requires calibration data"
  )
})

testthat::test_that("irrelevant prediction and imputation controls are rejected", {
  univariate <- mixedgp_test_univariate_fit()
  multivariate <- mixedgp_test_multivariate_fit()
  new_X <- matrix(
    c(-0.25, 0.25), ncol = 1L, dimnames = list(NULL, "x1")
  )
  new_C_1d <- data.frame(severity = ordered(
    c("high", "low"), levels = c("low", "high")
  ))
  new_C_multi <- data.frame(
    proxy_a = ordered(c("medium", "low"), levels = c("low", "medium")),
    proxy_b = ordered(c("mild", "none"), levels = c("none", "mild"))
  )

  testthat::expect_error(
    predict_eivgp(
      univariate, new_X, new_C_1d, new_U = c(0, 0), target = "mean"
    ),
    "new_U is not used"
  )
  testthat::expect_error(
    predict_eivgp(
      univariate, new_X, new_C_1d, new_U = c(0, 0), target = "surface"
    ),
    "new_C is not used"
  )
  testthat::expect_error(
    predict_eivgp(
      univariate, new_X, new_C_1d, target = "mean", new_U_scale = "model"
    ),
    "only when new_U"
  )
  testthat::expect_error(
    predict_eivgp(
      univariate, new_X, new_C_1d, target = "response",
      include_process_uncertainty = FALSE
    ),
    "applies only"
  )
  testthat::expect_error(
    predict_eivgp(
      univariate, new_X, new_C_1d, target = "response", n_latent = 4L
    ),
    "only to target='mean'"
  )
  testthat::expect_error(
    predict_eivgp(
      univariate, new_X, new_C_1d, target = "mean",
      latent_sampler = "rejection"
    ),
    "univariate engine"
  )
  testthat::expect_error(
    predict_eivgp(
      multivariate, new_X, new_C_multi, target = "mean",
      rejection_batch_size = 64L
    ),
    "only for latent_sampler='rejection'"
  )
  testthat::expect_error(
    impute_eivgp(univariate, n_per_draw = NA_integer_), "integer-valued"
  )
  testthat::expect_error(
    impute_eivgp(univariate, seed = 1L), "not used"
  )
  testthat::expect_error(
    impute_eivgp(
      univariate, new_C = new_C_1d, latent_sampler = "rejection"
    ),
    "univariate engine"
  )
  testthat::expect_error(
    impute_eivgp(
      multivariate, new_C = new_C_multi, rejection_max_batches = 2L
    ),
    "only for latent_sampler='rejection'"
  )
})

testthat::test_that("joint and pointwise prediction modes are labeled", {
  univariate <- mixedgp_test_univariate_fit()
  multivariate <- mixedgp_test_multivariate_fit()
  new_X <- matrix(
    c(-0.25, 0.25), ncol = 1L, dimnames = list(NULL, "x1")
  )
  new_C_1d <- data.frame(severity = ordered(
    c("high", "low"), levels = c("low", "high")
  ))
  new_C_multi <- data.frame(
    proxy_a = ordered(c("medium", "low"), levels = c("low", "medium")),
    proxy_b = ordered(c("mild", "none"), levels = c("none", "mild"))
  )

  pointwise_1d <- predict_eivgp(
    univariate, new_X, new_C_1d, target = "mean", draw_ids = 1L,
    n_latent = 4L, joint = FALSE, seed = 4207
  )
  joint_1d <- predict_eivgp(
    univariate, new_X, new_C_1d, target = "mean", draw_ids = 1L,
    n_latent = 4L, joint = TRUE, seed = 4207
  )
  pointwise_multi <- predict_eivgp(
    multivariate, new_X, new_C_multi, target = "mean", draw_ids = 1L,
    n_latent = 4L, joint = FALSE, seed = 4208
  )
  joint_multi <- predict_eivgp(
    multivariate, new_X, new_C_multi, target = "mean", draw_ids = 1L,
    n_latent = 4L, joint = TRUE, seed = 4208
  )
  testthat::expect_false(attr(pointwise_1d, "joint"))
  testthat::expect_true(attr(joint_1d, "joint"))
  testthat::expect_false(attr(pointwise_multi, "joint"))
  testthat::expect_true(attr(joint_multi, "joint"))
  testthat::expect_true(attr(joint_1d, "include_process_uncertainty"))
  testthat::expect_true(attr(joint_multi, "include_process_uncertainty"))
})

testthat::test_that("empirical-Q GP integration matches brute-force aggregation", {
  X_train <- matrix(c(-1, 0, 1), ncol = 1L)
  U_train <- matrix(c(-0.6, 0.1, 0.8), ncol = 1L)
  y_train <- c(-0.5, 0.1, 0.7)
  X_star <- matrix(c(-0.2, 0.4), ncol = 1L)
  U_mc <- array(NA_real_, dim = c(3L, 2L, 1L))
  U_mc[, 1L, 1L] <- c(-0.8, -0.1, 0.5)
  U_mc[, 2L, 1L] <- c(-0.3, 0.4, 1)
  logtheta <- log(c(1.2, 0.7, 0.9))
  sigma2_eps <- 0.2

  integrated <- gp_integrated_mean_state_general(
    X_train, U_train, y_train, X_star, U_mc,
    logtheta = logtheta, sigma2_eps = sigma2_eps, return_cov = TRUE
  )
  M <- dim(U_mc)[1L]
  N <- nrow(X_star)
  brute_mean <- vapply(seq_len(N), function(i) {
    mean(vapply(seq_len(M), function(a) {
      gp_predict_draw_general(
        X_train, U_train, y_train,
        X_star = X_star[i, , drop = FALSE],
        U_star = matrix(U_mc[a, i, ], nrow = 1L),
        logtheta = logtheta, sigma2_eps = sigma2_eps
      )$mean
    }, numeric(1L)))
  }, numeric(1L))
  brute_cov <- matrix(NA_real_, N, N)
  for (i in seq_len(N)) {
    for (j in seq_len(N)) {
      values <- numeric(M * M)
      pos <- 0L
      for (a in seq_len(M)) {
        for (b in seq_len(M)) {
          pos <- pos + 1L
          pair <- gp_predict_draw_general(
            X_train, U_train, y_train,
            X_star = rbind(
              X_star[i, , drop = FALSE], X_star[j, , drop = FALSE]
            ),
            U_star = matrix(
              c(U_mc[a, i, ], U_mc[b, j, ]), ncol = 1L
            ),
            logtheta = logtheta, sigma2_eps = sigma2_eps,
            return_cov = TRUE
          )
          values[pos] <- pair$cov[1L, 2L]
        }
      }
      brute_cov[i, j] <- mean(values)
    }
  }
  testthat::expect_equal(integrated$mean, brute_mean, tolerance = 1e-10)
  testthat::expect_equal(integrated$cov, brute_cov, tolerance = 1e-10)

  integrated_1d <- gp_integrated_mean_state_1d(
    X_train, as.numeric(U_train), y_train, X_star, U_mc[, , 1L],
    logtheta = logtheta, sigma2_eps = sigma2_eps, return_cov = TRUE
  )
  brute_mean_1d <- vapply(seq_len(N), function(i) {
    mean(vapply(seq_len(M), function(a) {
      gp_predict_draw(
        X_train, as.numeric(U_train), y_train,
        x_star = X_star[i, , drop = FALSE],
        u_star = U_mc[a, i, 1L],
        logtheta = logtheta, sigma2_eps = sigma2_eps
      )$mean
    }, numeric(1L)))
  }, numeric(1L))
  brute_cov_1d <- matrix(NA_real_, N, N)
  for (i in seq_len(N)) {
    for (j in seq_len(N)) {
      values <- numeric(M * M)
      pos <- 0L
      for (a in seq_len(M)) {
        for (b in seq_len(M)) {
          pos <- pos + 1L
          pair <- gp_predict_draw(
            X_train, as.numeric(U_train), y_train,
            x_star = rbind(
              X_star[i, , drop = FALSE], X_star[j, , drop = FALSE]
            ),
            u_star = c(U_mc[a, i, 1L], U_mc[b, j, 1L]),
            logtheta = logtheta, sigma2_eps = sigma2_eps,
            return_cov = TRUE
          )
          values[pos] <- pair$cov[1L, 2L]
        }
      }
      brute_cov_1d[i, j] <- mean(values)
    }
  }
  testthat::expect_equal(
    integrated_1d$mean, brute_mean_1d, tolerance = 1e-10
  )
  testthat::expect_equal(
    integrated_1d$cov, brute_cov_1d, tolerance = 1e-10
  )
})

testthat::test_that("both numerical studies use one validated interface", {
  s1 <- run_mixedgp_experiment("study1", dry_run = TRUE)
  s2 <- run_mixedgp_experiment("study2", dry_run = TRUE)
  testthat::expect_s3_class(s1, "mixedgp_experiment_spec")
  testthat::expect_s3_class(s2, "mixedgp_experiment_spec")
  testthat::expect_identical(s1$config, s2$config)
  testthat::expect_identical(s1$stages, s2$stages)

  thorough <- run_mixedgp_experiment(
    "study1", config = "thorough", stages = "data", dry_run = TRUE
  )
  run_env <- mixedgp_experiment_environment(thorough)
  mixedgp_assign_experiment_controls(thorough, run_env)
  testthat::expect_false(run_env$STUDY1_REUSE_LOCKED_EIV)

  abbreviated <- run_mixedgp_experiment(
    "study1", config = "t", stages = "data", dry_run = TRUE
  )
  quick <- run_mixedgp_experiment(
    "study1", config = "q", stages = "data", dry_run = TRUE
  )
  explicit <- run_mixedgp_experiment(
    "study1", config = "thorough", stages = "data",
    strict_competitors = FALSE, dry_run = TRUE
  )
  testthat::expect_identical(abbreviated$config, "thorough")
  testthat::expect_true(abbreviated$strict_competitors)
  testthat::expect_false(quick$strict_competitors)
  testthat::expect_false(explicit$strict_competitors)

  serial_spec <- run_mixedgp_experiment(
    "study1", stages = "data", parallel_level = "none", n_cores = 4L,
    dry_run = TRUE
  )
  testthat::expect_identical(mixedgp_experiment_data_cores(serial_spec), 1L)
  testthat::expect_false(mixedgp_parallel_chains_enabled("none", 4L))
  chain_spec <- run_mixedgp_experiment(
    "study1", stages = "data", parallel_level = "chains", n_cores = 4L,
    dry_run = TRUE
  )
  replication_spec <- run_mixedgp_experiment(
    "study1", stages = "data", parallel_level = "replications", n_cores = 2L,
    dry_run = TRUE
  )
  testthat::expect_identical(mixedgp_experiment_data_cores(chain_spec), 1L)
  testthat::expect_identical(
    mixedgp_experiment_data_cores(replication_spec),
    replication_spec$execution$data$cores
  )
  testthat::expect_identical(chain_spec$execution$data$backend, "serial")
  policy_scripts <- file.path(
    s1$code_dir,
    c(
      "01_study1_representative_figures.R",
      "02_study1_targeted_controls.R",
      "01_study2_representative_figures.R"
    )
  )
  for (path in policy_scripts) {
    text <- readLines(path, warn = FALSE)
    testthat::expect_true(any(grepl("mixedgp_parallel_chains_enabled", text)))
    testthat::expect_false(any(grepl("parallel::detectCores", text)))
  }

  testthat::expect_error(
    run_mixedgp_experiment(
      "study2", stages = "data", dry_run = TRUE,
      study_options = list(STUDY2_MC_N_REPP = 2L)
    ),
    "Unknown study_options"
  )
  testthat::expect_error(
    run_mixedgp_experiment(
      "study2", stages = "data", dry_run = TRUE,
      study_options = list(STUDY2_MC_N_REP = 2.5)
    ),
    "integer-valued"
  )
  testthat::expect_error(
    run_mixedgp_experiment(
      "study1", stages = "data", use_cache = 1, dry_run = TRUE
    ),
    "use_cache must be TRUE or FALSE"
  )
})

testthat::test_that("experiment execution preserves RNG and package search state", {
  testthat::skip_if(
    isTRUE(getOption("knitr.in.progress")),
    "requires the installed-style package namespace; exercised by R CMD check"
  )
  output_dir <- tempfile("eivGP-runner-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
  testthat::local_mocked_bindings(
    mixedgp_generate_experiment_data = function(spec) {
      set.seed(9010L)
      stats::runif(2L)
      data.frame()
    },
    .package = "eivGP"
  )
  set.seed(4114)
  expected <- stats::runif(4)
  set.seed(4114)
  invisible(run_mixedgp_experiment(
    "study1", stages = "data", output_dir = output_dir,
    parallel_level = "none"
  ))
  testthat::expect_identical(stats::runif(4), expected)

  testthat::skip_if_not_installed("ggplot2")
  spec <- run_mixedgp_experiment(
    "study1", stages = "data", output_dir = output_dir, dry_run = TRUE
  )
  run_env <- mixedgp_experiment_environment(spec)
  search_before <- search()
  eval(quote(suppressPackageStartupMessages(library(ggplot2))), run_env)
  testthat::expect_identical(search(), search_before)
  testthat::expect_true(exists("ggplot", envir = run_env, inherits = FALSE))
})

testthat::test_that("core-count environment variables are parsed strictly", {
  old <- Sys.getenv("MIXEDGP_CORES", unset = NA_character_)
  on.exit({
    if (is.na(old)) Sys.unsetenv("MIXEDGP_CORES") else
      Sys.setenv(MIXEDGP_CORES = old)
  }, add = TRUE)
  Sys.setenv(MIXEDGP_CORES = "2")
  testthat::expect_identical(mixedgp_resolve_cores(), 2L)
  Sys.setenv(MIXEDGP_CORES = "2.5")
  testthat::expect_error(mixedgp_resolve_cores(), "base-10 integer")
})

testthat::test_that("experiment defaults do not leak from the global workspace", {
  sentinel <- "MIXEDGP_GLOBAL_SENTINEL_FOR_TEST"
  assign(sentinel, TRUE, envir = .GlobalEnv)
  on.exit(rm(list = sentinel, envir = .GlobalEnv), add = TRUE)

  spec <- run_mixedgp_experiment(
    "study1", stages = "data", dry_run = TRUE
  )
  run_env <- mixedgp_experiment_environment(spec)
  testthat::expect_false(eval(substitute(exists(NAME), list(NAME = sentinel)), run_env))
  testthat::expect_error(eval(as.name(sentinel), run_env), "not found")
  assign(sentinel, TRUE, envir = run_env)
  testthat::expect_true(eval(substitute(exists(NAME), list(NAME = sentinel)), run_env))
  testthat::expect_true(eval(as.name(sentinel), run_env))
})

testthat::test_that("the installed experiment runner ignores cwd script copies", {
  testthat::skip_if(
    isTRUE(getOption("knitr.in.progress")),
    "requires installed system.file() semantics; exercised by R CMD check"
  )
  fake_dir <- tempfile("eivGP-fake-cwd-")
  dir.create(fake_dir)
  invisible(file.create(
    file.path(fake_dir, "01_study1_representative_figures.R")
  ))
  old_wd <- setwd(fake_dir)
  on.exit({
    setwd(old_wd)
    unlink(fake_dir, recursive = TRUE)
  }, add = TRUE)

  spec <- run_mixedgp_experiment("study1", stages = "data", dry_run = TRUE)
  installed_scripts <- normalizePath(
    system.file("experiments", package = "eivGP"),
    winslash = "/", mustWork = TRUE
  )
  testthat::expect_identical(spec$code_dir, installed_scripts)
})

testthat::test_that("publication simulation APIs resolve bundled sources", {
  testthat::skip_if(
    isTRUE(getOption("knitr.in.progress")),
    "requires installed system.file() semantics; exercised by R CMD check"
  )
  installed_sources <- normalizePath(
    system.file("simulation-source", package = "eivGP"),
    winslash = "/", mustWork = TRUE
  )
  study1 <- study1_simulation_config(mode = "dry_run", workers = 1L)
  study2 <- study2_simulation_config(mode = "dry_run", workers = 1L)
  publication1 <- study1_simulation_config(mode = "publication", workers = 1L)
  publication2 <- study2_simulation_config(mode = "publication", workers = 1L)
  testthat::expect_identical(study1$study, "study1")
  testthat::expect_identical(study2$study, "study2")
  testthat::expect_identical(study1$code_dir, installed_sources)
  testthat::expect_identical(study2$code_dir, installed_sources)
  testthat::expect_length(study1$stages, 0L)
  testthat::expect_length(study2$stages, 0L)
  testthat::expect_false(publication1$fail_closed)
  testthat::expect_false(publication2$fail_closed)
  testthat::expect_true(publication1$mcmc$require_gate)
  testthat::expect_true(publication2$mcmc$require_gate)

  printed <- utils::capture.output(run <- run_study1_simulation(study1))
  testthat::expect_true(length(printed) > 0L)
  testthat::expect_true(is.list(run$task_plan))
})

testthat::test_that("external competitor validators are archived and nonfatal", {
  validator_files <- c(
    "08_published_competitor_validation.R",
    "09_experiment_design_validation.R"
  )
  testthat::expect_true(all(
    validator_files %in% mixedgp_simulation_modules()
  ))

  engine <- new.env(parent = emptyenv())
  engine$run_published_competitor_validation <- function() {
    data.frame(
      method = c("UC-GP", "LVGP", "EzGP"),
      status = c("success", "success", "unavailable_or_failed"),
      stringsAsFactors = FALSE
    )
  }
  engine$run_experiment_design_validation <- function() {
    data.frame(
      validator = "paired_experiment_design",
      pass = TRUE,
      stringsAsFactors = FALSE
    )
  }

  smoke_dir <- tempfile("eivGP-validator-smoke-")
  publication_dir <- tempfile("eivGP-validator-publication-")
  on.exit(unlink(c(smoke_dir, publication_dir), recursive = TRUE), add = TRUE)

  smoke <- mixedgp_run_automatic_validators(
    list(mode = "smoke"), engine, smoke_dir
  )
  testthat::expect_true(all(smoke$status$pass))
  testthat::expect_true(all(file.exists(file.path(
    smoke_dir, "config",
    c(
      "published_competitor_validation.csv",
      "experiment_design_validation.csv",
      "automatic_validation_status.csv"
    )
  ))))

  publication <- mixedgp_run_automatic_validators(
    list(mode = "publication"), engine, publication_dir
  )
  testthat::expect_true(all(publication$status$pass))
  status <- utils::read.csv(file.path(
    publication_dir, "config", "automatic_validation_status.csv"
  ))
  testthat::expect_true(status$pass[status$validator == "published_competitors"])
  testthat::expect_match(
    status$message[status$validator == "published_competitors"],
    "External competitor validation did not succeed"
  )
})
