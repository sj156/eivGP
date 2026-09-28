f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
project <- dirname(dirname(normalizePath(sub("^--file=", "", f[1]))))
source(file.path(project, "applications/ADNI_parallel.R"))
for (n in 1:64) {
  p <- adni_core_plan(n)
  stopifnot(p$active_chain_limit <= n, p$chain_workers >= 1, p$chain_workers <= 4, p$fold_workers <= 5)
}
for (n in c(12L,16L,56L)) {
  p <- adni_core_plan(n)
  expected <- switch(as.character(n), `12`=c(3L,4L),`16`=c(4L,4L),`56`=c(5L,4L))
  stopifnot(identical(c(p$fold_workers,p$chain_workers),expected))
}
stopifnot(adni_core_plan(56L,os_type="windows")$active_chain_limit==1L)
work <- function() {
  root <- tempfile("adni-scheduler-"); dir.create(root)
  on.exit(unlink(root,recursive=TRUE),add=TRUE)
  jobs <- data.frame(fold=1:5,output_root=root)
  status <- file.path(root,'status.csv')
  ans <- adni_run_jobs(jobs,function(fold) {
    Sys.sleep(.05);c(fold=fold,pid=Sys.getpid())
  },adni_core_plan(12L),status)
  stopifnot(length(ans)==5L,all(read.csv(status)$success),
    identical(vapply(ans,function(x) as.integer(x['fold']),integer(1)),1:5))
  stopifnot(!any(file.exists(file.path(root,paste0('fold_',1:5),'.fold-lock'))))
  fail <- try(adni_run_jobs(jobs,function(fold) {
    if (fold==2L) stop('intentional failure'); fold
  },adni_core_plan(12L),status),silent=TRUE)
  stopifnot(inherits(fail,'try-error'),sum(read.csv(status)$success)==4L,
    file.exists(file.path(root,'fold_2/WORKER_FAILURE.txt')))
  locked <- file.path(root,'fold_1/.fold-lock');dir.create(locked)
  fail <- try(adni_run_jobs(jobs[1,,drop=FALSE],function(...) stop('must not execute'),
    adni_core_plan(4L,1L),status),silent=TRUE)
  stopifnot(inherits(fail,'try-error'),dir.exists(locked))
}
work()
cat('PASS: budgets 1-64; five-fold dispatch; failure aggregation; lock ownership.\n')
