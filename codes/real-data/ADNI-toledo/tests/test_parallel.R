f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
project <- dirname(dirname(normalizePath(sub("^--file=", "", f[1]))))
source(file.path(project, "adni_parallel.R"))
for (n in 1:64) {
  p <- adni_core_plan(n, 6L, 4L, 6L)
  stopifnot(p$active_chain_limit <= n, p$chain_workers >= 1, p$chain_workers <= 4,
            p$fold_workers <= 6)
}
for (n in c(12L, 16L, 56L)) {
  p <- adni_core_plan(n, 6L, 4L, 6L)
  expected <- switch(as.character(n), `12`=c(3L,4L),`16`=c(4L,4L),`56`=c(6L,4L))
  stopifnot(identical(c(p$fold_workers,p$chain_workers),expected))
}
stopifnot(adni_core_plan(56L,6L,4L,6L,os_type="windows")$active_chain_limit==1L)
root <- tempfile("adni-scheduler-"); dir.create(root)
work <- function() {
  on.exit(unlink(root,recursive=TRUE),add=TRUE)
  jobs <- expand.grid(fold=1:3,repeat_id=2:3)
  jobs$output_root <- file.path(root,paste0('repeat',jobs$repeat_id))
  for (r in unique(jobs$output_root)) dir.create(r)
  status <- file.path(root,'status.csv')
  ans <- adni_run_jobs(jobs,function(repeat_id,fold) {
    Sys.sleep(.15)
    c(repeat_id=repeat_id,fold=fold,pid=Sys.getpid())
  },adni_core_plan(16L,6L,4L,6L),status)
  stopifnot(length(ans)==6L, all(read.csv(status)$success),
    length(unique(vapply(ans,function(x) x['pid'],numeric(1))))>1L)
  stopifnot(!any(file.exists(file.path(jobs$output_root,paste0('fold_',jobs$fold),'.fold-lock'))))
  fail <- try(adni_run_jobs(jobs,function(repeat_id,fold) {
    if (repeat_id==2L && fold==2L) stop('intentional failure')
    paste(repeat_id,fold)
  },adni_core_plan(12L,6L,4L,6L),status),silent=TRUE)
  stopifnot(inherits(fail,'try-error'),sum(read.csv(status)$success)==5L,
    file.exists(file.path(root,'repeat2/fold_2/WORKER_FAILURE.txt')))
  # Existing fold lock is never taken over or removed.
  locked <- file.path(root,'repeat2/fold_1/.fold-lock');dir.create(locked)
  fail <- try(adni_run_jobs(jobs[1,,drop=FALSE],function(...) stop('must not execute'),
    adni_core_plan(4L,1L),status),silent=TRUE)
  stopifnot(inherits(fail,'try-error'),dir.exists(locked));unlink(locked,recursive=TRUE)
}
if (.Platform$OS.type!='windows') work() else unlink(root,recursive=TRUE)
cat('PASS: allocations for 1-64 cores; 12/16/56-core plans; cross-repeat dispatch; failure aggregation; lock ownership.\n')
