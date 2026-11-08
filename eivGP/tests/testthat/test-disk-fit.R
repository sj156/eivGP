test_that('fresh disk fits reproduce every retained draw and terminal state', {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  set.seed(1701)
  X <- matrix(rnorm(12)); U <- matrix(rnorm(12)); observed <- U; observed[5:12, ] <- NA
  for (engine in c('univariate', 'multivariate')) for (sampler in c('ESS', 'HMC', 'HYBRID')) {
    C <- if (engine == 'univariate') matrix(as.integer(U[, 1] > 0) + 1L) else
      cbind(rep(1:3, 4), rep(3:1, 4))
    args <- list(X = X, y = cos(X[, 1]) + U[, 1], C = C, U_obs = observed,
      engine = engine, m_vec = if (engine == 'univariate') 2L else c(3L, 3L),
      sampler = sampler, n_iter = 25L, burn = 8L, n_chains = 2L, seed = 31L,
      parallel = FALSE, sampler_control = if (sampler == 'HYBRID') list(hybrid_every = 3L) else list())
    if (engine == 'multivariate') { args$ident <- 'none'; args$store_scores <- TRUE }
    expected <- suppressWarnings(do.call(fit_eivgp, args))
    for (chunk in c(3L, 7L)) {
      path <- file.path(root, paste(engine, sampler, chunk, 'rds', sep = '.'))
      set.seed(727); rng <- .Random.seed
      index <- suppressWarnings(do.call(fit_eivgp, c(args, list(storage = 'disk', path = path, chunk_iter = chunk))))
      expect_identical(.Random.seed, rng)
      expect_s3_class(index, 'eivgp_disk_fit')
      expect_null(index$mcmc$samples_by_chain)
      expect_true(all(index$disk_history$blocks$to - index$disk_history$blocks$from + 1L <= chunk))
      actual <- read_fit(path)
      expect_identical(actual$mcmc$samples_by_chain, expected$mcmc$samples_by_chain, info = paste(engine, sampler, chunk))
      expect_identical(actual$checkpoint$states, expected$checkpoint$states)
      expect_equal(actual$mcmc$chain_stats, expected$mcmc$chain_stats, tolerance = 1e-13)
      expect_identical(actual$control, expected$control)
      expect_identical(actual$priors, expected$priors)
      expect_identical(actual$data, expected$data)
      expect_identical(actual$diagnostics$ess_key, expected$diagnostics$ess_key)
      expect_equal(index$checkpoint$iteration, args$n_iter)
      expect_equal(index$checkpoint$burn, args$burn)
      expect_equal(length(index$checkpoint$states), args$n_chains)
    }
  }
})

test_that('fresh disk fits preserve RNG without an existing caller seed and refuse overwrite', {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  args <- list(X = matrix(1:10 / 10), y = sin(1:10), C = matrix(rep(1:2, 5)),
    engine = 'univariate', n_iter = 12L, burn = 4L, n_chains = 2L, parallel = FALSE)
  path <- file.path(root, 'fit.rds')
  rng <- mixedgp_rng_state(); on.exit(mixedgp_restore_rng_state(rng), add = TRUE)
  if (exists('.Random.seed', .GlobalEnv)) rm('.Random.seed', envir = .GlobalEnv)
  do.call(fit_eivgp, c(args, list(storage = 'disk', path = path, chunk_iter = 3L)))
  expect_false(exists('.Random.seed', .GlobalEnv, inherits = FALSE))
  hash <- tools::md5sum(path)
  expect_error(do.call(fit_eivgp, c(args, list(storage = 'disk', path = path))), 'already exists')
  expect_identical(tools::md5sum(path), hash)
  expect_error(do.call(fit_eivgp, c(args, list(storage = 'disk'))), 'checkpoint path')
  expect_error(do.call(fit_eivgp, c(args, list(path = file.path(root, 'unused.rds')))), 'requires storage')
})

test_that('fresh disk fits can resume an interrupted committed block exactly', {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  args <- list(X = matrix(1:12 / 12), y = cos(1:12), C = matrix(rep(1:3, 4)),
    engine = 'univariate', sampler = 'HYBRID', sampler_control = list(hybrid_every = 3L),
    n_iter = 26L, burn = 6L, n_chains = 4L, parallel = FALSE, seed = 117L)
  expected <- do.call(fit_eivgp, args)
  path <- file.path(root, 'fit.rds'); commit <- mixedgp_disk_commit
  expect_error(with_mocked_bindings(do.call(fit_eivgp, c(args, list(storage = 'disk', path = path, chunk_iter = 4L))),
    mixedgp_disk_commit = function(index, path) {
      if (index$checkpoint$iteration >= 16L) stop('interrupted fresh fit')
      commit(index, path)
    }, .package = 'eivGP'), 'interrupted fresh fit')
  last <- mixedgp_disk_index(path)$checkpoint$iteration
  continue_eivgp_disk(path, args$n_iter - last, chunk_iter = 4L, parallel = TRUE, n_cores = 2L)
  actual <- read_fit(path)
  expect_identical(actual$mcmc$samples_by_chain, expected$mcmc$samples_by_chain)
  expect_identical(actual$checkpoint$states, expected$checkpoint$states)
})

test_that('scientific ESS uses the same checks and extension schedule in disk mode', {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  args <- list(X = matrix(1:12 / 12), y = cos(1:12), C = matrix(rep(1:3, 4)),
    engine = 'univariate', n_iter = 20L, burn = 6L, n_chains = 4L,
    parallel = FALSE, seed = 119L, ess_min = 1L,
    ess_panel = list(X = matrix(c(.2, .7)), C = matrix(c(1L, 3L)), n_latent = 2L, seed = 11L))
  run <- function(extra = list()) {
    checks <- 0L
    with_mocked_bindings(do.call(fit_eivgp, c(args, extra)),
      mixedgp_scientific_ess_pass = function(...) { checks <<- checks + 1L; checks > 1L },
      mixedgp_ess_extension_size = function(...) 4L, .package = 'eivGP')
  }
  expected <- suppressWarnings(run())
  path <- file.path(root, 'fit.rds')
  index <- suppressWarnings(run(list(storage = 'disk', path = path, chunk_iter = 3L)))
  actual <- read_fit(path)
  expect_identical(actual$mcmc$samples_by_chain, expected$mcmc$samples_by_chain)
  expect_identical(actual$checkpoint$states, expected$checkpoint$states)
  expect_identical(actual$scientific_ess, expected$scientific_ess)
  expect_true(isTRUE(index$scientific_ess$success))
  expect_null(index$scientific_ess$target_series)
  with_mocked_bindings(mixedgp_application_disk_ess(path, 1L, args$ess_panel, 1L, 3L),
    mixedgp_finish_ess_disk = function(...) stop('matching science should be reused'), .package = 'eivGP')
  expect_identical(read_fit(path)$scientific_ess, expected$scientific_ess)
})

test_that('copying disk and legacy fits preserves their original archives', {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  args <- list(X = matrix(1:12 / 12), y = cos(1:12), C = matrix(rep(1:3, 4)),
    engine = 'univariate', n_iter = 18L, burn = 6L, n_chains = 2L, parallel = FALSE, seed = 117L)
  for (storage in c('memory', 'disk')) {
    source <- file.path(root, paste0(storage, '.rds'))
    if (storage == 'memory') save_fit(do.call(fit_eivgp, args), source) else
      do.call(fit_eivgp, c(args, list(storage = 'disk', path = source, chunk_iter = 3L)))
    expected <- continue_eivgp(read_fit(source), 4L, parallel = FALSE)
    hash <- tools::md5sum(source)
    destination <- file.path(root, paste0(storage, '-continued.rds'))
    mixedgp_disk_copy(source, destination)
    continue_eivgp_disk(destination, 4L, chunk_iter = 3L, parallel = FALSE)
    actual <- read_fit(destination)
    expect_identical(tools::md5sum(source), hash)
    expect_identical(mixedgp_saved_iteration(source), 18L)
    expect_identical(actual$mcmc$samples_by_chain, expected$mcmc$samples_by_chain)
    expect_identical(actual$checkpoint$states, expected$checkpoint$states)
  }
})

test_that('disk fitting bounds sampler buffers regardless of the total iteration budget', {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  sizes <- integer(); past <- logical(); sampler <- mixedgp_v030_fit
  index <- with_mocked_bindings(fit_eivgp(matrix(1:10 / 10), sin(1:10), matrix(rep(1:2, 5)),
    engine = 'univariate', n_iter = 95L, burn = 4L, n_chains = 1L, parallel = FALSE,
    storage = 'disk', path = file.path(root, 'fit.rds'), chunk_iter = 7L),
    mixedgp_v030_fit = function(...) {
      a <- list(...); start <- if (is.null(a$.resume)) 0L else a$.resume$iteration
      sizes <<- c(sizes, a$n_iter - max(a$burn, start))
      past <<- c(past, !is.null(a$.resume$old_chains))
      sampler(...)
    }, .package = 'eivGP')
  expect_gt(length(sizes), 10L)
  expect_true(all(sizes <= 7L))
  expect_false(any(past))
  expect_equal(index$checkpoint$iteration, 95L)
  expect_equal(summary(index)$mcmc$total_saved_draws, 91L)
  expect_output(print(index), 'retained draws stored on disk')
})
