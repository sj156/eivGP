test_that('disk continuation preserves both engines, samplers and all posterior draws', {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  set.seed(617); n <- 12L; X <- matrix(rnorm(n), ncol = 1)
  U <- matrix(rnorm(n), ncol = 1); observed <- U; observed[5:n, ] <- NA
  for (engine in c('univariate', 'multivariate')) for (sampler in c('ESS', 'HMC', 'HYBRID')) {
    C <- if (engine == 'univariate') matrix(as.integer(U[, 1] > 0) + 1L) else
      cbind(rep(1:3, 4), rep(3:1, 4))
    args <- list(X = X, y = X[, 1] + U[, 1], C = C, U_obs = observed,
      engine = engine, m_vec = if (engine == 'univariate') 2L else c(3L, 3L),
      sampler = sampler, sampler_control = if (sampler == 'HYBRID') list(hybrid_every = 3L) else list(),
      n_iter = 20L, burn = 8L, n_chains = 2L, seed = 731L, parallel = FALSE)
    if (engine == 'multivariate') { args$ident <- 'none'; args$store_scores <- TRUE }
    fit <- suppressWarnings(do.call(fit_eivgp, args))
    expected <- suppressWarnings(continue_eivgp(fit, 8L, parallel = FALSE, n_cores = 1L))
    path <- file.path(root, paste0(engine, '-', sampler, '.rds')); save_fit(fit, path)
    old_hash <- unname(tools::md5sum(path))
    set.seed(909); caller_rng <- .Random.seed
    index <- continue_eivgp_disk(path, 8L, chunk_iter = 3L, parallel = FALSE, n_cores = 1L)
    expect_identical(.Random.seed, caller_rng)
    expect_identical(unname(tools::md5sum(file.path(paste0(path, '.draws'), 'legacy-fit.rds'))), old_hash)
    expect_null(index$mcmc$samples_by_chain)
    expect_true(all(index$disk_history$blocks$to - index$disk_history$blocks$from + 1L <= 3L))
    if (sampler != 'ESS') expect_true(all(vapply(index$checkpoint$states,
      function(x) nrow(x$state$hmc$history) == 0L, logical(1))))
    actual <- read_fit(path)
    expect_identical(actual$mcmc$samples_by_chain, expected$mcmc$samples_by_chain, info = paste(engine, sampler))
    expect_equal(actual$mcmc$chain_stats, expected$mcmc$chain_stats, tolerance = 1e-13)
    expect_identical(actual$checkpoint$states, expected$checkpoint$states)
    expect_identical(actual$control, fit$control)
    expect_identical(actual$priors, fit$priors)
    expect_identical(actual$data, fit$data)
    expect_identical(actual$diagnostics$ess_key, expected$diagnostics$ess_key)
    expect_error(save_fit(actual, path), 'Disk histories')
    later <- continue_eivgp(actual, 2L, parallel = FALSE, n_cores = 1L)
    continue_eivgp_disk(path, 2L, chunk_iter = 2L, parallel = FALSE, n_cores = 1L)
    expect_identical(read_fit(path)$checkpoint$states, later$checkpoint$states)
    if (engine == 'multivariate') {
      selected <- read_fit(path, draw_ids = c(1:4, 23:26))
      expect_identical(selected$mcmc$samples_by_chain$S[[1]], later$mcmc$samples_by_chain$S[[1]][1:4, , , drop = FALSE])
      expect_false(isTRUE(attr(selected, 'full_retained_history')))
      lean <- read_fit(path, fields = setdiff(index$disk_history$fields, 'S'))
      expect_true(all(vapply(lean$mcmc$samples_by_chain$S, is.null, logical(1))))
      expect_identical(lean$diagnostics$ess_key, later$diagnostics$ess_key)
      predict <- function(fit) predict_eivgp(fit, new_X = X[1:2, , drop = FALSE],
        new_C = C[1:2, , drop = FALSE], draw_ids = 1:4, seed = 51L)
      expect_identical(predict(lean), predict(later))
    }
  }
})

test_that('disk transactions resume after interruption and reject modified blocks', {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  X <- matrix(seq_len(10) / 10); C <- matrix(rep(1:2, 5))
  fit <- fit_eivgp(X, sin(X[, 1]), C, engine = 'univariate', sampler = 'HMC',
    n_iter = 16L, burn = 4L, n_chains = 2L, parallel = FALSE, seed = 19L)
  path <- file.path(root, 'fit.rds'); save_fit(fit, path)
  continue_eivgp_disk(path, 0L, chunk_iter = 4L)
  original_commit <- mixedgp_disk_commit
  expect_error(with_mocked_bindings(continue_eivgp_disk(path, 8L, chunk_iter = 4L, parallel = FALSE),
    mixedgp_disk_commit = function(index, path) stop('interrupted before commit'),
    .package = 'eivGP'), 'interrupted before commit')
  expect_identical(mixedgp_disk_index(path)$checkpoint$iteration, 16L)
  # Orphan draw blocks must neither advance RNG nor duplicate committed history.
  continue_eivgp_disk(path, 8L, chunk_iter = 4L, parallel = FALSE, n_cores = 1L)
  expected <- continue_eivgp(fit, 8L, parallel = FALSE, n_cores = 1L)
  expect_identical(read_fit(path)$mcmc$samples_by_chain, expected$mcmc$samples_by_chain)
  expect_identical(read_fit(path)$checkpoint$states, expected$checkpoint$states)
  # Interruption after committing the archive index, before the small pointer.
  write <- mixedgp_atomic_rds
  expect_error(with_mocked_bindings(continue_eivgp_disk(path, 4L, chunk_iter = 4L, parallel = FALSE),
    mixedgp_atomic_rds = function(x, destination) {
      if (identical(destination, path)) stop('interrupted pointer update')
      write(x, destination)
    }, .package = 'eivGP'), 'interrupted pointer update')
  expect_identical(read_fit(path)$checkpoint$iteration, 28L)
  expected <- continue_eivgp(expected, 4L, parallel = FALSE, n_cores = 1L)
  expect_identical(read_fit(path)$checkpoint$states, expected$checkpoint$states)
  block <- mixedgp_disk_index(path)$disk_history$blocks$file[1L]
  cat('changed', file = file.path(paste0(path, '.draws'), block), append = TRUE)
  expect_error(continue_eivgp_disk(path, 1L), 'missing or changed')
})

test_that('four-chain disk continuation is invariant to worker allocation', {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  X <- matrix(seq_len(12) / 12); C <- cbind(rep(1:3, 4), rep(3:1, 4))
  U <- X; U[5:12, ] <- NA
  fit <- fit_eivgp(X, cos(X[, 1]), C, U_obs = U, ident = 'none', engine = 'multivariate',
    store_scores = TRUE, sampler = 'ESS', n_iter = 16L, burn = 4L,
    n_chains = 4L, parallel = FALSE, seed = 99L)
  path <- file.path(root, 'fit.rds'); save_fit(fit, path)
  expected <- continue_eivgp(fit, 8L, parallel = FALSE, n_cores = 1L)
  continue_eivgp_disk(path, 8L, chunk_iter = 4L, parallel = TRUE, n_cores = 2L)
  actual <- read_fit(path)
  expect_identical(actual$mcmc$samples_by_chain, expected$mcmc$samples_by_chain)
  expect_equal(actual$mcmc$chain_stats, expected$mcmc$chain_stats, tolerance = 1e-13)
  expect_identical(actual$checkpoint$states, expected$checkpoint$states)
})
