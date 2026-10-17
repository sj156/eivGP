test_that("disk restart preserves warmup, draws, counters and terminal RNG", {
  root <- tempfile("restart-"); dir.create(root)
  on.exit(unlink(root,recursive=TRUE))
  old <- options(eivGP.checkpoint_root=NULL,eivGP.checkpoint_seconds=300,
    eivGP.checkpoint_iterations=2L)
  on.exit(options(old),add=TRUE)
  env <- new.env(parent=environment(mixedgp_v030_fit))
  fit <- mixedgp_v030_fit; environment(fit) <- env
  args <- list(X=matrix(seq(-1,1,length.out=6),ncol=1),y=c(-.2,.1,.4,.8,.3,-.1),
    C=matrix(c(1,1,2,2,3,3),ncol=1),U_obs=matrix(NA_real_,6,1),m_vec=3L,
    measurement="threshold",n_iter=10L,burn=4L,n_chains=2L,seed=736L,
    parallel_chains=FALSE,priors=list(u_dictionary=list(c(.2,.8))))
  expected <- do.call(fit,args)
  writer <- mixedgp_checkpoint_write
  for(at in c(2L,6L)) {
    options(eivGP.checkpoint_root=file.path(root,paste0("at",at)))
    env$mixedgp_checkpoint_write <- function(object,path) {
      writer(object,path)
      if(!is.null(object$iteration) && object$iteration==at)stop("simulated outage")
    }
    expect_error(do.call(fit,args),"simulated outage")
    env$mixedgp_checkpoint_write <- writer
    actual <- do.call(fit,args)
    expect_identical(actual$chains,expected$chains)
    expect_identical(do.call(fit,args)$chains,expected$chains)
    csv <- list.files(getOption("eivGP.checkpoint_root"),pattern="timing.csv$",recursive=TRUE,full.names=TRUE)
    expect_length(csv,2L)
    z <- read.csv(csv[1]); expect_equal(z$iteration,10)
    expect_equal(z$checkpoint_iteration,10);expect_equal(z$remaining_seconds,0)
  }
  # Model changes get a different checkpoint identity, never stale draws.
  directories <- list.dirs(getOption("eivGP.checkpoint_root"),recursive=FALSE)
  args$y[1] <- args$y[1]+.01
  do.call(fit,args)
  expect_length(list.dirs(getOption("eivGP.checkpoint_root"),recursive=FALSE),length(directories)+1L)
})

test_that("batch forecast is grouped, conservative, and excludes cache-load time", {
  root<-tempfile();on.exit(unlink(root,recursive=TRUE))
  record<-mixedgp_batch_timing(root,c("a","b","c","d"),c("FULL","FULL","CUT","CUT"))
  record(1,"complete",10)
  latest<-function()read.csv(list.files(root,"^forecast-[0-9]+[.]csv$",full.names=TRUE)[1])
  expect_true(is.na(latest()$remaining_seconds))
  record(3,"complete",100)
  expect_equal(latest()$remaining_seconds,110)
  record(1,"running");record(1,"cached",.001)
  expect_equal(readRDS(file.path(root,"a.rds"))$duration_seconds,10)
  expect_equal(latest()$remaining_seconds,110)
  record(2,"complete",20);record(4,"failed",30)
  expect_equal(latest()$remaining_seconds,100)
  expect_equal(latest()$recorded_task_seconds,160)
  record(4,"complete",80)
  expect_equal(latest()$remaining_seconds,0)
  expect_equal(latest()$completed,4)
})

test_that("CUT measurement restarts without changing selected states or traces", {
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  old<-options(eivGP.checkpoint_root=NULL,eivGP.checkpoint_iterations=2L)
  on.exit(options(old),add=TRUE)
  env<-new.env(parent=environment(eivgp_numerical_api))
  maker<-eivgp_numerical_api;environment(maker)<-env
  api<-maker(tempdir());api$F<-function(name)get(name,envir=env,inherits=TRUE)
  C<-matrix(rep(1:4,8),8,4); U<-matrix(NA_real_,8,2)
  U[1:2,]<-cbind(c(-1,1),c(.5,-.5))
  cfg<-list(iter=10L,warmup=4L,outer=8L,cores=1L)
  reference<-suppressWarnings(api$measurement_cut(C,U,cfg,123L))
  writer<-mixedgp_checkpoint_write
  for(at in c(2L,6L)) {
    options(eivGP.checkpoint_root=file.path(root,paste0("at",at)))
    env$mixedgp_checkpoint_write<-function(object,path) {
      writer(object,path)
      if(!is.null(object$iteration)&&object$iteration==at)stop("simulated outage")
    }
    expect_error(api$measurement_cut(C,U,cfg,123L),"simulated outage")
    env$mixedgp_checkpoint_write<-writer
    actual<-suppressWarnings(api$measurement_cut(C,U,cfg,123L))
    expect_identical(actual,reference)
  }
})

test_that("default numerical paths are in home and overrides remain explicit", {
  api<-eivgp_numerical_api(tempdir())
  expect_identical(api$numerical_default_output("study1"),file.path(path.expand("~/eivGP-results"),"study1"))
  expect_identical(api$numerical_default_data("study2"),file.path(path.expand("~/eivGP-results"),"data","study2"))
})

test_that("forked chain checkpoints reproduce serial draws and can be reused", {
  skip_on_os("windows")
  root<-tempfile();on.exit(unlink(root,recursive=TRUE))
  old<-options(eivGP.checkpoint_root=NULL,eivGP.checkpoint_iterations=2L)
  on.exit(options(old),add=TRUE)
  args<-list(X=matrix(seq(-1,1,length.out=6),ncol=1),y=c(-.2,.1,.4,.8,.3,-.1),
    C=matrix(c(1,1,2,2,3,3),ncol=1),U_obs=matrix(NA_real_,6,1),m_vec=3L,
    measurement="threshold",n_iter=8L,burn=4L,n_chains=2L,seed=736L,
    parallel_chains=FALSE,priors=list(u_dictionary=list(c(.2,.8))))
  reference<-do.call(mixedgp_v030_fit,args)
  options(eivGP.checkpoint_root=root)
  args$parallel_chains<-TRUE;args$n_cores<-2L
  actual<-do.call(mixedgp_v030_fit,args)
  expect_identical(actual$chains,reference$chains)
  expect_identical(do.call(mixedgp_v030_fit,args)$chains,reference$chains)
})

test_that("replication timing preserves failures and cache-hit return values", {
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  fun<-function(x) { if(x==2)stop("expected failure"); z<-list(value=x);attr(z,"eivgp_cache_hit")<-TRUE;z }
  map<-function(X,FUN,...)lapply(X,FUN)
  expect_warning(out<-mixedgp_run_replications(as.list(1:3),fun,1L,1:3,
    file.path(root,"status.csv"),"test",map),"replications failed")
  expect_identical(out,list(list(value=1L),list(value=3L)))
  status<-read.csv(file.path(root,"status.csv"))
  expect_identical(status$status,c("success","failed","success"))
  expect_true(file.exists(file.path(root,"status.csv.timing","X2.rds")))
})
