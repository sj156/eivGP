# Recovery must preserve active, remote and newly acquired workflow locks.
repo <- Sys.getenv('EIVGP_TEST_REPO')
stopifnot(nzchar(repo))
e <- new.env(parent = globalenv())
sys.source(file.path(repo, 'tools/resume-oai-tools.R'), e)

testthat::test_that('resume archives only confirmed stale local locks', {
  testthat::skip_if_not(dir.exists('/proc'))
  root <- tempfile(); dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  lock <- file.path(root, '_meta/locks/eivGP.lock')
  owner_file <- file.path(lock, 'owner.json')
  dead <- 999999999L
  stopifnot(!dir.exists(file.path('/proc', dead)))
  owner <- list(pid = dead, host = unname(Sys.info()[['nodename']]), token = 'original-owner')
  write_owner <- function(value = owner) {
    dir.create(lock, recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(value, owner_file, auto_unbox = TRUE)
  }
  # Current process and a recycled/existing PID both require a refusal.
  live <- owner; live$pid <- Sys.getpid(); write_owner(live)
  before <- tools::md5sum(owner_file)
  testthat::expect_error(e$resume_oai_recover_locks(root, 'eivGP'), 'process IDs still exist')
  testthat::expect_identical(tools::md5sum(owner_file), before)
  remote <- owner; remote$host <- paste0(owner$host, '-other-host'); write_owner(remote)
  testthat::expect_error(e$resume_oai_recover_locks(root, 'eivGP'), 'lock host')
  testthat::expect_true(dir.exists(lock))
  invalid <- owner; invalid$token <- NULL; write_owner(invalid)
  testthat::expect_error(e$resume_oai_recover_locks(root, 'eivGP'), 'lock token')
  write_owner()
  phase <- file.path(root, 'logs/eivGP_fixture/S50000/fit_status.csv')
  dir.create(dirname(phase), recursive = TRUE)
  write.csv(data.frame(state = 'RUNNING', worker_pid = Sys.getpid()), phase, row.names = FALSE)
  testthat::expect_error(e$resume_oai_recover_locks(root, 'eivGP'), 'process IDs still exist')
  testthat::expect_true(dir.exists(lock))
  write.csv(data.frame(state = 'RUNNING', worker_pid = dead), phase, row.names = FALSE)
  # The other method group's lock is independent and must remain untouched.
  other <- file.path(root, '_meta/locks/competitors.lock'); dir.create(other)
  jsonlite::write_json(live, file.path(other, 'owner.json'), auto_unbox = TRUE)
  archived <- e$resume_oai_recover_locks(root, 'eivGP')
  testthat::expect_length(archived, 1L)
  testthat::expect_false(dir.exists(lock))
  testthat::expect_identical(jsonlite::fromJSON(file.path(archived, 'owner.json')), owner)
  testthat::expect_true(dir.exists(other))
  # Legacy and modern lock formats can be recovered together after a crash.
  write_owner()
  legacy <- file.path(root, '.application.lock'); dir.create(legacy)
  writeLines(as.character(dead), file.path(legacy, 'pid'))
  archived <- e$resume_oai_recover_locks(root, 'eivGP')
  testthat::expect_length(archived, 2L)
  testthat::expect_false(dir.exists(legacy))
  testthat::expect_false(dir.exists(lock))
  testthat::expect_true(dir.exists(other))
})

testthat::test_that('resume detects ownership changes during stale-lock checks', {
  testthat::skip_if_not(dir.exists('/proc'))
  root <- tempfile(); dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  lock <- file.path(root, '_meta/locks/eivGP.lock'); dir.create(lock, recursive = TRUE)
  path <- file.path(lock, 'owner.json')
  owner <- list(pid = 999999999L, host = unname(Sys.info()[['nodename']]), token = 'old')
  jsonlite::write_json(owner, path, auto_unbox = TRUE)
  scoped <- new.env(parent = globalenv())
  sys.source(file.path(repo, 'tools/resume-oai-tools.R'), scoped)
  scoped$system2 <- function(...) {
    owner$token <- 'new'
    jsonlite::write_json(owner, path, auto_unbox = TRUE)
    base::system2(...)
  }
  testthat::expect_error(scoped$resume_oai_recover_locks(root, 'eivGP'), 'ownership changed')
  testthat::expect_true(dir.exists(lock))
  testthat::expect_identical(jsonlite::fromJSON(path)$token, 'new')
})

testthat::test_that('unreaped zombie workers do not block checkpoint recovery', {
  testthat::skip_if_not(dir.exists('/proc'))
  python <- Sys.which('python3')
  testthat::skip_if_not(nzchar(python))
  root <- tempfile(); dir.create(root)
  stop_file <- file.path(root, 'reap')
  pid <- NA_integer_
  on.exit({
    file.create(stop_file)
    for (i in 1:100) {
      if (is.na(pid) || is.null(e$resume_oai_process(pid))) break
      Sys.sleep(.01)
    }
    unlink(root, recursive = TRUE)
  })
  fixture <- file.path(root, 'zombie.py')
  writeLines(c('import os, sys, time', 'from pathlib import Path',
    'root = Path(sys.argv[1])', 'pid = os.fork()',
    'if pid == 0: os._exit(0)',
    '(root / "pid").write_text(str(pid))',
    'for i in range(1000):', '    if (root / "reap").exists(): break',
    '    time.sleep(.01)', 'os.waitpid(pid, 0)'), fixture)
  system2(python, shQuote(c(fixture, root)), wait = FALSE)
  for (i in 1:200) {
    if (file.exists(file.path(root, 'pid'))) {
      pid <- suppressWarnings(as.integer(readLines(file.path(root, 'pid'), warn = FALSE)))
      if (length(pid) == 1L && !is.na(pid) && identical(e$resume_oai_process(pid)$state, 'Z')) break
    }
    Sys.sleep(.01)
  }
  testthat::expect_identical(e$resume_oai_process(pid)$state, 'Z')
  lock <- file.path(root, '_meta/locks/eivGP.lock'); dir.create(lock, recursive = TRUE)
  owner <- list(pid = pid, host = unname(Sys.info()[['nodename']]), token = 'zombie-owner')
  jsonlite::write_json(owner, file.path(lock, 'owner.json'), auto_unbox = TRUE)
  phase <- file.path(root, 'logs/eivGP_fixture/S50000/fit_status.csv')
  dir.create(dirname(phase), recursive = TRUE)
  write.csv(data.frame(state = 'RUNNING', worker_pid = pid), phase, row.names = FALSE)
  archived <- e$resume_oai_recover_locks(root, 'eivGP')
  testthat::expect_length(archived, 1L)
  testthat::expect_false(dir.exists(lock))
  testthat::expect_identical(jsonlite::fromJSON(file.path(archived, 'owner.json')), owner)
})
