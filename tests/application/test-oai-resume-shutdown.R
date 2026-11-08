# A failed coordinator must not leave nested chain workers blocking the next run.
repo <- Sys.getenv('EIVGP_TEST_REPO')
stopifnot(nzchar(repo))
e <- new.env(parent = globalenv())
sys.source(file.path(repo, 'tools/resume-oai-tools.R'), e)

testthat::test_that('launcher cleans nested workers after the session leader exits', {
  testthat::skip_if_not(dir.exists('/proc'))
  testthat::skip_if_not(nzchar(Sys.which('setsid')))
  root <- tempfile(); dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  fixture <- file.path(root, 'worker.R')
  code <- c(
    'a <- commandArgs(TRUE)',
    'source(a[1])',
    'resume_oai_record_session(Sys.getenv("EIVGP_OAI_SESSION_FILE"))',
    'job <- parallel::mcparallel({',
    '  inner <- parallel::mcparallel({',
    '    writeLines(as.character(Sys.getpid()), file.path(a[2], "leaf"))',
    '    Sys.sleep(30)',
    '  }, mc.set.seed = FALSE)',
    '  writeLines(as.character(Sys.getpid()), file.path(a[2], "parent"))',
    '  Sys.sleep(30)',
    '}, mc.set.seed = FALSE)',
    'for (i in 1:200) {',
    '  if (file.exists(file.path(a[2], "leaf"))) break',
    '  Sys.sleep(.01)',
    '}',
    'stopifnot(file.exists(file.path(a[2], "leaf")))',
    'quit(save = "no", status = 7L, runLast = FALSE)')
  writeLines(code, fixture)
  status <- suppressWarnings(e$resume_oai_launch_worker(fixture,
    c(file.path(repo, 'tools/resume-oai-tools.R'), root)))
  testthat::expect_identical(status, 7L)
  pids <- as.integer(vapply(file.path(root, c('parent', 'leaf')), readLines, character(1)))
  info <- lapply(pids, e$resume_oai_process)
  testthat::expect_true(all(vapply(info, function(x) is.null(x) || x$state %in% c('Z', 'X'), logical(1))),
    info = paste(capture.output(str(info)), collapse = '\n'))
  testthat::expect_identical(Sys.getenv('EIVGP_OAI_SESSION_FILE'), '')
})

testthat::test_that('cleanup refuses a changed session identity', {
  testthat::skip_if_not(dir.exists('/proc'))
  path <- tempfile(); on.exit(unlink(path))
  owner <- e$resume_oai_process(Sys.getpid())
  owner$group <- owner$session <- owner$pid
  owner$start <- 'changed-start-time'
  saveRDS(owner, path)
  testthat::expect_error(e$resume_oai_cleanup_session(path), 'identity changed')
  testthat::expect_identical(e$resume_oai_cleanup_session(paste0(path, '-absent')), integer())
})

testthat::test_that('isolated worker recognizes its own reparented resume launcher', {
  testthat::skip_if_not(dir.exists('/proc'))
  testthat::skip_if_not(nzchar(Sys.which('setsid')))
  root <- tempfile(); dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  worker <- file.path(root, 'worker.R')
  launcher <- file.path(root, 'resume-oai-eivGP.r')
  helper <- file.path(repo, 'tools/resume-oai-tools.R')
  writeLines(c('a <- commandArgs(TRUE); source(a[1])',
    'resume_oai_record_session(Sys.getenv("EIVGP_OAI_SESSION_FILE"))',
    'resume_oai_recover_locks(a[2], "eivGP")'), worker)
  writeLines(c('a <- commandArgs(TRUE); source(a[1])',
    'status <- resume_oai_launch_worker(a[2], c(a[1], a[3]))',
    'quit(save = "no", status = status)'), launcher)
  status <- suppressWarnings(system2(file.path(R.home('bin'), 'Rscript'),
    shQuote(c('--vanilla', launcher, helper, worker, root))))
  testthat::expect_identical(status, 0L)
})

testthat::test_that('interrupting the launcher also stops nested workers', {
  testthat::skip_if_not(dir.exists('/proc'))
  testthat::skip_if_not(nzchar(Sys.which('setsid')))
  root <- tempfile(); dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  worker <- file.path(root, 'worker.R'); launcher <- file.path(root, 'launcher.R')
  helper <- file.path(repo, 'tools/resume-oai-tools.R')
  writeLines(c(
    'a <- commandArgs(TRUE); source(a[1])',
    'resume_oai_record_session(Sys.getenv("EIVGP_OAI_SESSION_FILE"))',
    'job <- parallel::mcparallel({',
    '  inner <- parallel::mcparallel({',
    '    writeLines(as.character(Sys.getpid()), file.path(a[2], "leaf"))',
    '    Sys.sleep(30)',
    '  }, mc.set.seed = FALSE)',
    '  writeLines(as.character(Sys.getpid()), file.path(a[2], "parent"))',
    '  Sys.sleep(30)',
    '}, mc.set.seed = FALSE)',
    'for (i in 1:200) {',
    '  if (file.exists(file.path(a[2], "leaf"))) break',
    '  Sys.sleep(.01)',
    '}',
    'stopifnot(file.exists(file.path(a[2], "leaf")))',
    'tools::pskill(as.integer(a[3]), tools::SIGINT)',
    'Sys.sleep(30)'), worker)
  writeLines(c('a <- commandArgs(TRUE); source(a[1])',
    'status <- resume_oai_launch_worker(a[2], c(a[1], a[3], as.character(Sys.getpid())))',
    'quit(save = "no", status = status)'), launcher)
  status <- suppressWarnings(system2(file.path(R.home('bin'), 'Rscript'),
    shQuote(c('--vanilla', launcher, helper, worker, root))))
  testthat::expect_true(status != 0L)
  pids <- as.integer(vapply(file.path(root, c('parent', 'leaf')), readLines, character(1)))
  info <- lapply(pids, e$resume_oai_process)
  testthat::expect_true(all(vapply(info, function(x) is.null(x) || x$state %in% c('Z', 'X'), logical(1))),
    info = paste(capture.output(str(info)), collapse = '\n'))
})
