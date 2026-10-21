test_that("Study I measurement resumes exactly during warmup and sampling", {
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  old<-options(eivGP.checkpoint_seconds=Inf,eivGP.checkpoint_iterations=2L,
    eivGP.checkpoint_after_seconds=0,eivGP.checkpoint_script_start=0)
  on.exit(options(old),add=TRUE)
  env<-new.env(parent=environment(fit_threshold_measurement_response_free))
  sampler<-fit_threshold_measurement_response_free;environment(sampler)<-env
  args<-list(c_ord=rep(1:6,2),u_obs=rep(NA_real_,12),calib_idx=integer(),m=6L,n_iter=12L,burn=4L,seed=321L)
  expected<-do.call(sampler,args)
  writer<-mixedgp_checkpoint_write
  for(at in c(2L,8L)) {
    path<-file.path(root,paste0('chain-',at,'.rds'))
    env$mixedgp_checkpoint_write<-function(object,path) {
      writer(object,path)
      if(!is.null(object$iteration)&&object$iteration==at)stop('simulated outage')
    }
    expect_error(do.call(sampler,c(args,list(.checkpoint_path=path))),'simulated outage')
    cp<-readRDS(path);expect_equal(nrow(cp$draws$samples_u),max(0L,at-args$burn))
    env$mixedgp_checkpoint_write<-writer
    actual<-do.call(sampler,c(args,list(.checkpoint_path=path)))
    expect_identical(actual,expected)
    expect_identical(do.call(sampler,c(args,list(.checkpoint_path=path))),expected)
    timing<-read.csv(sub('[.]rds$','-timing.csv',path))
    expect_equal(timing$iteration,12);expect_equal(timing$remaining_seconds,0)
  }
  # Even a forced chain-completion save is gated before six hours.
  options(eivGP.checkpoint_seconds=10800,eivGP.checkpoint_iterations=Inf,
    eivGP.checkpoint_after_seconds=21600,eivGP.checkpoint_script_start=proc.time()[['elapsed']])
  path<-file.path(root,'delayed.rds')
  expect_identical(do.call(sampler,c(args,list(.checkpoint_path=path))),expected)
  expect_false(file.exists(path))
  expect_true(file.exists(sub('[.]rds$','-timing.csv',path)))
  # Simulate a six-hour boundary using the same schedule function as production.
  due<-mixedgp_checkpoint_due
  env$mixedgp_checkpoint_due<-function(now,last_save,iterations_since,force=FALSE)
    due(21600,21599,iterations_since,force)
  options(eivGP.checkpoint_script_start=0)
  env$mixedgp_checkpoint_write<-function(object,path){writer(object,path);if(!is.null(object$iteration))stop('six hour save')}
  path<-file.path(root,'scheduled.rds')
  expect_error(do.call(sampler,c(args,list(.checkpoint_path=path))),'six hour save')
  expect_equal(readRDS(path)$iteration,1)
})

test_that("Study I completed measurement publication retires only exact restart copies", {
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  old<-options(eivGP.checkpoint_root=file.path(root,'checkpoints'),eivGP.checkpoint_seconds=Inf,
    eivGP.checkpoint_iterations=2L,eivGP.checkpoint_after_seconds=0)
  on.exit(options(old),add=TRUE)
  cfg<-list(iter=12L,warmup=4L,outer=8L,cores=1L)
  C<-rep(1:6,2);U<-rep(NA_real_,12)
  fit<-suppressWarnings(mixedgp_study1_cc_measurement(C,U,cfg,123L))
  disk<-attr(fit,'eivgp_study1_restart_storage')$directory
  files<-file.path(disk,sprintf('chain-%02d.rds',1:4))
  expect_true(all(file.exists(files)))
  api<-eivgp_numerical_api(tempdir());path<-file.path(root,'measurement.rds')
  good<-readRDS(files[1]);bad<-good;bad$draws$samples_u[1,1]<-999
  saveRDS(bad,files[1])
  expect_warning(api$atomic_rds(fit,path),'differs')
  expect_true(all(file.exists(files)))
  saveRDS(good,files[1]);api$atomic_rds(fit,path)
  expect_identical(readRDS(path),fit)
  expect_false(any(file.exists(files)))
  expect_true(file.exists(file.path(disk,'identity.rds')))
  expect_true(file.exists(file.path(disk,'chain-01-timing.csv')))
  expect_identical(readRDS(file.path(disk,'retired.rds'))$status,'removed')
  options(eivGP.checkpoint_root=NULL)
  reference<-suppressWarnings(mixedgp_study1_cc_measurement(C,U,cfg,123L))
  attr(fit,'eivgp_study1_restart_storage')<-NULL
  expect_identical(fit,reference)
})


test_that("Study I workflow restores the script clock between invocations", {
  old<-options(eivGP.checkpoint_script_start=NULL);on.exit(options(old),add=TRUE)
  env<-new.env(parent=environment(run_paper_study));runner<-run_paper_study;environment(runner)<-env
  seen<-NULL
  env$mixedgp_study1_cc_workflow<-function(...) {seen<<-getOption("eivGP.checkpoint_script_start");invisible(NULL)}
  suppressMessages(runner("study1",list(component="complete-case",action="plan"),tempdir()))
  expect_true(is.finite(seen))
  expect_null(getOption("eivGP.checkpoint_script_start"))
  options(eivGP.checkpoint_script_start=123)
  suppressMessages(runner("study1",list(component="complete-case",action="plan"),tempdir()))
  expect_equal(seen,123)
  expect_equal(getOption("eivGP.checkpoint_script_start"),123)
})
