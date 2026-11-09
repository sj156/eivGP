legacy_disk_library <- function() {
  library <- Sys.getenv('EIVGP_LEGACY_LIBRARY')
  skip_if_not(dir.exists(file.path(library, 'eivGP')), 'Original installed package fixture not supplied')
  normalizePath(library)
}

test_that('new disk indexes preserve the call source without a namespace closure', {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  path <- file.path(root, 'fit.rds')
  X <- matrix(seq_len(12) / 12); C <- matrix(rep(1:2, 6))
  fit <- do.call(fit_eivgp, list(X = X, y = sin(X[, 1]), C = C, engine = 'univariate',
    n_iter = 16L, burn = 4L, n_chains = 2L, parallel = FALSE, seed = 19L))
  expect_true(is.function(fit$call[[1L]]))
  expected <- fit$call
  save_fit(fit, path)
  index <- continue_eivgp_disk(path, 0L, chunk_iter = 4L)
  expect_true(is.call(index$call[[1L]]))
  expect_identical(index$call[[1L]][[1L]], as.name('function'))
  expect_identical(as.list(index$call)[-1L], as.list(expected)[-1L])
  expect_identical(index$call[[1L]][[2L]], formals(expected[[1L]]))
  expect_identical(index$call[[1L]][[3L]], body(expected[[1L]]))
  expect_identical(read_fit(path)$mcmc$samples_by_chain, fit$mcmc$samples_by_chain)
  # A different installed package must be able to validate the same bytes.
  library <- legacy_disk_library()
  code <- 'a<-commandArgs(TRUE);.libPaths(c(a[2],.libPaths()));getFromNamespace("mixedgp_disk_index","eivGP")(a[1],verify=TRUE);invisible(NULL)'
  expect_equal(system2(file.path(R.home('bin'), 'Rscript'), shQuote(c('--vanilla', '-e', code, path, library)),
    stdout = file.path(root, 'old-reader.log'), stderr = file.path(root, 'old-reader.log')), 0L)
})

test_that('legacy version-bound signatures are validated without modifying original files', {
  library <- legacy_disk_library()
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  path <- file.path(root, 'fit.rds')
  code <- paste0('a<-commandArgs(TRUE);.libPaths(c(a[2],.libPaths()));',
    'X<-matrix(seq_len(12)/12);C<-matrix(rep(1:2,6));',
    'fit<-do.call(eivGP::fit_eivgp,list(X=X,y=sin(X[,1]),C=C,engine="univariate",',
    'sampler="HMC",n_iter=16L,burn=4L,n_chains=2L,parallel=FALSE,seed=19L));',
    'eivGP::save_fit(fit,a[1]);eivGP::continue_eivgp_disk(a[1],4L,chunk_iter=4L,parallel=FALSE);',
    'fit<-eivGP::read_fit(a[1]);saveRDS(list(draws=fit$mcmc$samples_by_chain,states=fit$checkpoint$states),paste0(a[1],".expected"))')
  expect_equal(system2(file.path(R.home('bin'), 'Rscript'), shQuote(c('--vanilla', '-e', code, path, library)),
    stdout = file.path(root, 'fixture.log'), stderr = file.path(root, 'fixture.log')), 0L)
  files <- c(path, list.files(paste0(path, '.draws'), full.names = TRUE))
  hashes <- tools::md5sum(files)
  expect_error(mixedgp_disk_index(path), 'Invalid or modified')
  expect_true(mixedgp_disk_prepare_compatibility(path, library))
  expect_identical(tools::md5sum(files), hashes)
  expect_false(mixedgp_disk_prepare_compatibility(path, library))
  index <- mixedgp_disk_index(path, verify = TRUE)
  expect_identical(index$checkpoint$iteration, 20L)
  actual <- read_fit(path); expected <- readRDS(paste0(path, '.expected'))
  expect_identical(actual$mcmc$samples_by_chain, expected$draws)
  expect_identical(actual$checkpoint$states, expected$states)
  copy <- file.path(root, 'copy.rds')
  copied <- mixedgp_disk_copy(path, copy)
  expect_true(is.call(copied$call[[1L]]))
  expect_identical(read_fit(copy)$mcmc$samples_by_chain, expected$draws)
  # The compatibility certificate cannot accept a changed terminal state.
  changed <- index; changed$checkpoint$states[[1L]]$rng[2L] <- 99L
  index_path <- file.path(paste0(path, '.draws'), 'index.rds')
  backup <- file.path(root, 'index-original.rds'); file.copy(index_path, backup)
  saveRDS(changed, index_path)
  expect_error(mixedgp_disk_index(path), 'Invalid or modified')
  expect_error(mixedgp_disk_prepare_compatibility(path, library), 'Legacy disk validation failed')
  # Even the exact certified index must still pass immutable-block checksums.
  file.copy(backup, index_path, overwrite = TRUE)
  expect_identical(unname(tools::md5sum(index_path)), unname(hashes[index_path]))
  block <- file.path(paste0(path, '.draws'), index$disk_history$blocks$file[1L])
  block_backup <- file.path(root, 'block-original.rds'); file.copy(block, block_backup)
  con <- file(block, 'r+b'); seek(con, 40L); writeBin(as.raw(0:3), con); close(con)
  expect_error(mixedgp_disk_index(path, verify = TRUE), 'checksum differs')
  file.copy(block_backup, block, overwrite = TRUE)
  later <- continue_eivgp(actual, 4L, parallel = FALSE, n_cores = 1L)
  continue_eivgp_disk(path, 4L, chunk_iter = 4L, parallel = FALSE, n_cores = 1L)
  resumed <- read_fit(path)
  expect_identical(resumed$mcmc$samples_by_chain, later$mcmc$samples_by_chain)
  expect_identical(resumed$checkpoint$states, later$checkpoint$states)
  expect_true(is.call(readRDS(index_path)$call[[1L]]))
})
