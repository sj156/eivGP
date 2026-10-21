test_that("default GP policy never writes or reads partial chains and restarts a failed fit", {
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  old<-options(eivGP.checkpoint_root=root,eivGP.checkpoint_seconds=NULL,eivGP.checkpoint_iterations=NULL,eivGP.checkpoint_after_seconds=NULL,eivGP.checkpoint_script_start=NULL)
  on.exit(options(old),add=TRUE)
  env<-new.env(parent=environment(mixedgp_v030_fit));fit<-mixedgp_v030_fit;environment(fit)<-env
  args<-list(X=matrix(seq(-1,1,length.out=6),ncol=1),y=c(-.2,.1,.4,.8,.3,-.1),
    C=matrix(c(1,1,2,2,3,3),ncol=1),U_obs=matrix(NA_real_,6,1),m_vec=3L,
    measurement="threshold",n_iter=8L,burn=4L,n_chains=2L,seed=736L,
    parallel_chains=FALSE,priors=list(u_dictionary=list(c(.2,.8))))
  calls<-0L;sweep<-mixedgp_v030_sweep
  env$mixedgp_v030_sweep<-function(...) {calls<<-calls+1L;if(calls==11L)stop("outage");sweep(...)}
  expect_error(do.call(fit,args),"outage")
  expect_length(list.files(root,'^chain-[0-9]+[.]rds$',recursive=TRUE),0L)
  # A legacy checkpoint must not be read, overwritten or removed by the default.
  dir<-list.dirs(root,recursive=FALSE)[1];legacy<-file.path(dir,'chain-01.rds')
  writeLines('legacy checkpoint placeholder',legacy);hash<-tools::md5sum(legacy)
  calls<-0L;env$mixedgp_v030_sweep<-function(...) {calls<<-calls+1L;sweep(...)}
  actual<-do.call(fit,args)
  expect_identical(calls,16L)
  expect_identical(tools::md5sum(legacy),hash)
  expect_null(actual$restart_storage)
  timing<-list.files(root,'timing.csv$',recursive=TRUE,full.names=TRUE)
  expect_length(timing,2L)
  expect_true(all(vapply(timing,function(p)tail(read.csv(p)$iteration,1)==8,logical(1))))
  options(eivGP.checkpoint_root=NULL)
  expected<-do.call(fit,args)
  expect_identical(actual$chains,expected$chains)
  options(eivGP.checkpoint_root=file.path(root,'delayed'),eivGP.checkpoint_seconds=10800,
    eivGP.checkpoint_after_seconds=21600,eivGP.checkpoint_script_start=proc.time()[['elapsed']])
  delayed<-do.call(fit,args)
  expect_identical(delayed$chains,expected$chains)
  expect_null(delayed$restart_storage)
  expect_length(list.files(file.path(root,'delayed'),'^chain-[0-9]+[.]rds$',recursive=TRUE),0L)
})

test_that("default measurement policy saves timing but no partial-chain files", {
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  old<-options(eivGP.checkpoint_root=root,eivGP.checkpoint_seconds=NULL,eivGP.checkpoint_iterations=NULL,eivGP.checkpoint_after_seconds=NULL,eivGP.checkpoint_script_start=NULL)
  on.exit(options(old),add=TRUE)
  api<-eivgp_numerical_api(tempdir())
  C<-matrix(rep(1:4,8),8,4);U<-matrix(NA_real_,8,2);U[1:2,]<-cbind(c(-1,1),c(.5,-.5))
  cfg<-list(iter=8L,warmup=4L,outer=8L,cores=1L)
  actual<-suppressWarnings(api$measurement_cut(C,U,cfg,123L))
  expect_length(list.files(root,'^chain-[0-9]+[.]rds$',recursive=TRUE),0L)
  expect_length(list.files(root,'timing.csv$',recursive=TRUE),4L)
  options(eivGP.checkpoint_root=NULL)
  expect_identical(actual,suppressWarnings(api$measurement_cut(C,U,cfg,123L)))
})

test_that("automatic checkpoint schedule waits six hours and uses three-hour script boundaries", {
  old<-options(eivGP.checkpoint_seconds=10800,eivGP.checkpoint_iterations=Inf,
    eivGP.checkpoint_after_seconds=21600,eivGP.checkpoint_script_start=100)
  on.exit(options(old),add=TRUE)
  due<-function(t,last,force=FALSE)mixedgp_checkpoint_due(100+t,100+last,999999,force)
  expect_false(due(10800,0))
  expect_false(due(21599,0,TRUE))
  expect_true(due(21600,20000))
  expect_false(due(21601,21600))
  expect_false(due(32399,21600))
  expect_true(due(32400,30000))
  expect_true(due(43200,32400))
  expect_true(due(25000,21600,TRUE))
  options(eivGP.checkpoint_seconds=Inf)
  expect_false(due(43200,0,TRUE))
})
