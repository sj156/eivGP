# Parallel conversion must preserve histories and leave the parent history-free.
repo <- Sys.getenv('EIVGP_TEST_REPO')
stopifnot(nzchar(repo))
e <- new.env(parent = globalenv())
sys.source(file.path(repo, 'tools/resume-oai-tools.R'), e)

testthat::test_that('legacy import respects cores and estimated available RAM', {
  slots <- function(cores, sizes, available) e$resume_oai_migration_capacity(cores, sizes, available)$workers
  testthat::expect_identical(slots(48L, rep(1.5, 12), 113), 12L)
  testthat::expect_identical(slots(8L, rep(1.5, 12), 113), 8L)
  testthat::expect_identical(slots(4L, rep(1.5, 12), 113), 4L)
  testthat::expect_identical(slots(24L, rep(1.5, 12), 113), 12L)
  testthat::expect_identical(slots(48L, rep(1.5, 3), 113), 3L)
  testthat::expect_identical(slots(48L, rep(1.5, 12), 32), 2L)
  testthat::expect_identical(slots(48L, rep(1.5, 12), 23), 0L)
  testthat::expect_identical(slots(48L, rep(1.5, 12), NA_real_), 1L)
  testthat::expect_identical(slots(48L, rep(4, 12), 48), 2L)
  testthat::expect_identical(slots(48L, numeric(), 113), 0L)
})

testthat::test_that('isolated parallel imports preserve every draw, RNG and original file', {
  testthat::skip_on_os('windows')
  root <- tempfile(); dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  scoped <- new.env(parent = globalenv())
  sys.source(file.path(repo, 'tools/resume-oai-tools.R'), scoped)
  scoped$resume_oai_available_gib <- function() 113
  X <- matrix(seq_len(12) / 12); C <- matrix(rep(1:3, 4))
  fits <- lapply(c('ESS', 'HMC'), function(sampler) eivGP::fit_eivgp(X, cos(X[, 1]), C,
    engine = 'univariate', sampler = sampler, n_iter = 16L, burn = 4L,
    n_chains = 2L, parallel = FALSE, seed = 93L))
  paths <- file.path(root, paste0('fit', 1:6, '.rds'))
  for (i in seq_along(paths)) eivGP::save_fit(fits[[1L + (i %% 2L)]], paths[i])
  hashes <- unname(tools::md5sum(paths))
  set.seed(901); rng <- .Random.seed
  status <- scoped$resume_oai_migrate(paths, 5L, 3L, root)
  testthat::expect_identical(.Random.seed, rng)
  testthat::expect_true(all(status$state == 'DONE'))
  testthat::expect_length(unique(status$worker_pid), 6L)
  testthat::expect_false(Sys.getpid() %in% status$worker_pid)
  for (i in seq_along(paths)) {
    actual <- eivGP::read_fit(paths[i]); expected <- fits[[1L + (i %% 2L)]]
    testthat::expect_identical(actual$mcmc$samples_by_chain, expected$mcmc$samples_by_chain)
    testthat::expect_identical(actual$checkpoint$states, expected$checkpoint$states)
    testthat::expect_identical(unname(tools::md5sum(file.path(paste0(paths[i], '.draws'),
      'legacy-fit.rds'))), hashes[i])
    testthat::expect_true(file.info(paths[i])$size < file.info(file.path(paste0(paths[i],
      '.draws'), 'legacy-fit.rds'))$size)
  }
  logs <- list.files(file.path(root, 'logs'), 'fit_status[.]csv', recursive = TRUE, full.names = TRUE)
  testthat::expect_length(logs, 1L)
  testthat::expect_identical(read.csv(logs, stringsAsFactors = FALSE)$state, rep('DONE', 6L))
})

testthat::test_that('failed import records the error and preserves the legacy checkpoint', {
  testthat::skip_on_os('windows')
  root <- tempfile(); dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  scoped <- new.env(parent = globalenv())
  sys.source(file.path(repo, 'tools/resume-oai-tools.R'), scoped)
  scoped$resume_oai_available_gib <- function() 113
  path <- file.path(root, 'invalid.rds'); saveRDS(list(invalid = TRUE), path)
  hash <- unname(tools::md5sum(path))
  testthat::expect_error(scoped$resume_oai_migrate(path, 4L, 3L, root), 'migration failed')
  testthat::expect_identical(unname(tools::md5sum(path)), hash)
  logs <- list.files(file.path(root, 'logs'), 'fit_status[.]csv', recursive = TRUE, full.names = TRUE)
  tab <- read.csv(logs, stringsAsFactors = FALSE)
  testthat::expect_identical(tab$state, 'FAILED')
  testthat::expect_match(tab$error, 'legacy fit')
  scoped$resume_oai_available_gib <- function() 20
  testthat::expect_error(scoped$resume_oai_migrate(path, 4L, 3L, root), 'Insufficient available RAM')
  testthat::expect_identical(unname(tools::md5sum(path)), hash)
})
