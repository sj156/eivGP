# OAI evaluates up to eight independent histories within the same CPU pool.
repo <- Sys.getenv('EIVGP_TEST_REPO'); stopifnot(nzchar(repo))
source(file.path(repo, 'tools/oai.R'))
api <- oai_load_api(repo)

testthat::test_that('OAI evaluation uses chain-sized groups with an eight-job ceiling', {
  for (cores in c(4L, 8L, 12L, 24L, 32L, 36L, 48L, 64L)) {
    testthat::expect_identical(api$oai_resident_workers(cores, 'evaluate', 15L, 4L),
      as.integer(min(8L, cores %/% 4L)))
    testthat::expect_identical(api$oai_resident_workers(cores, 'fit', 30L, 4L), cores %/% 4L)
  }
  testthat::expect_identical(api$oai_resident_workers(32L, 'evaluate', 2L, 4L), 2L)
  testthat::expect_identical(api$oai_resident_workers(32L, 'evaluate', 0L, 4L), 0L)
  testthat::expect_error(api$oai_resident_workers(2.5, 'evaluate', 15L, 4L), 'integer')
})

testthat::test_that('eight evaluation jobs and 32 nested workers overlap without exceeding the shared budget', {
  testthat::skip_on_os('windows')
  testthat::skip_if(getFromNamespace('mixedgp_available_cores', 'eivGP')() < 32L,
    'This concurrency check needs 32 available cores')
  root <- tempfile(); dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  ready <- file.path(root, 'ready'); dir.create(ready)
  close <- api$mixedgp_budget_enter(32L); on.exit(close(), add = TRUE)
  tasks <- data.frame(plan = 'synthetic', fold = 1:8, cell = 'O20', method = 'EIV-GP',
    kind = 'eiv', stage = 'S20', id = paste0('task', 1:8))
  wait_for <- function(pattern, count) {
    start <- proc.time()[['elapsed']]
    while (length(list.files(ready, pattern)) < count) {
      if (proc.time()[['elapsed']] - start > 20) stop('Concurrent worker barrier timed out')
      Sys.sleep(.01)
    }
  }
  work <- function(task, cores) {
    stopifnot(cores == 4L)
    saveRDS(Sys.getpid(), file.path(ready, paste0('job-', task$fold, '.rds')))
    wait_for('^job-', 8L)
    leaves <- getFromNamespace('mixedgp_parallel_lapply', 'eivGP')(1:4, function(chain) {
      pool <- getOption('eivGP.core_pool')
      saveRDS(list(pid = Sys.getpid(), slots = pool$cores),
        file.path(ready, paste0('chain-', task$fold, '-', chain, '.rds')))
      wait_for('^chain-', 32L)
      # All 32 children must hold their pool slot at the same time.
      stopifnot(length(list.files(pool$path, '^[0-9]+$')) <= 32L)
      chain * chain
    }, n_cores = cores, .preserve_rng = TRUE)
    stopifnot(identical(as.numeric(unlist(leaves)), (1:4)^2))
    invisible('DONE')
  }
  set.seed(217L); rng <- .Random.seed
  status <- api$oai_screen_phase(tasks, 'evaluate', 32L, api, work,
    file.path(root, 'logs'), chains = 4L, root = root)
  testthat::expect_true(all(status$state == 'DONE'))
  testthat::expect_length(unique(status$worker_pid), 8L)
  testthat::expect_false(Sys.getpid() %in% status$worker_pid)
  testthat::expect_identical(.Random.seed, rng)
  children <- lapply(list.files(ready, '^chain-', full.names = TRUE), readRDS)
  testthat::expect_length(unique(vapply(children, `[[`, integer(1), 'pid')), 32L)
  testthat::expect_true(all(vapply(children, `[[`, integer(1), 'slots') == 32L))
  record <- jsonlite::fromJSON(file.path(root, '_meta/eivGP_status.json'))
  testthat::expect_identical(record$completed, 8L)
  testthat::expect_identical(record$failed, 0L)
  testthat::expect_length(record$active_jobs, 0L)
})
