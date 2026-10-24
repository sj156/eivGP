repo<-Sys.getenv('EIVGP_TEST_REPO');stopifnot(nzchar(repo))
source(file.path(repo,'tools/ocean-tools.R'),local=TRUE)
testthat::test_that("Ocean continuation accepts supported totals before held-out prediction", {
  # Use an existing checkpoint and stop at prediction; no MCMC or GP is fitted.
  prediction_calls<-0L
  testthat::local_mocked_bindings(predict_eivgp=function(...) {
    prediction_calls<<-prediction_calls+1L
    stop("held-out prediction reached")
  },.package="eivGP")
  root<-tempfile("ocean-continuation-");dir.create(root)
  on.exit(unlink(root,recursive=TRUE),add=TRUE)
  dir.create(file.path(root,"data"))
  output<-file.path(root,"continued");dir.create(output)
  cell<-file.path(output,"fold01/O20");dir.create(cell,recursive=TRUE)
  diagnostic<-data.frame(parameter="test",rhat=1,ess_bulk=800,
    ess_tail=800,target_pass=TRUE)
  fit<-list(diagnostics=list(rhat_hyper=diagnostic,rhat_tau=diagnostic,rhat_u=diagnostic),
    mcmc=list(samples_tau=matrix(1,20L,1L)))
  saveRDS(fit,file.path(cell,"fit_plus20k.rds"))
  env<-new.env(parent=environment(ocean_initial_functions))
  env$root<-root;env$out<-output;env$continuation_iter<-80L
  env$atomic_csv<-mixedgp_atomic_csv
  env$set_task<-env$write_progress<-function(...)invisible(NULL)
  ocean_continuation_functions(env)
  for(n in c(200L,500L,1000L)) {
    records<-data.frame(profile=paste0("profile",seq_len(n)),
      U=rep(c(25,75,125,165),each=n/4L),
      G2salinity=35,log_pressure=log(650))
    held<-unlist(lapply(split(seq_len(n),records$U),function(ids)head(ids,n/20L)))
    train<-records[-held,,drop=FALSE];test<-records[held,,drop=FALSE]
    train_file<-file.path(root,"data/fold01_train.csv")
    write.csv(train,train_file,row.names=FALSE)
    write.csv(test,file.path(root,"data/fold01_test.csv"),row.names=FALSE)
    testthat::expect_error(env$run_one(1L,"O20"),"held-out prediction reached")
    write.csv(train[-1,,drop=FALSE],train_file,row.names=FALSE)
    testthat::expect_error(env$run_one(1L,"O20"),"nrow")
  }
  testthat::expect_identical(prediction_calls,3L)
})


testthat::test_that("extension budgets are chosen after fitting and saved per task", {
  root<-tempfile("ocean-budget-");dir.create(root)
  on.exit(unlink(root,recursive=TRUE),add=TRUE)
  testthat::expect_identical(ocean_continuation_policy("paper")$n_iter,20000L)
  testthat::expect_identical(ocean_continuation_policy("pilot")$n_iter,80L)
  testthat::expect_null(ocean_settings()$continuation_iter)
  for(n in list(0,-1,1.5,NA_real_,Inf,"100",numeric(),2^31))
    testthat::expect_error(ocean_settings(continuation_iter=n),"positive integer")
  settings<-ocean_settings(profile="paper",action="continue",continuation_iter=12000)
  testthat::expect_identical(settings$continuation_iter,12000L)
  testthat::local_mocked_bindings(eivgp_module_identity=function()"test",.package="eivGP")
  id<-ocean_identity(root,root,ocean_settings())
  testthat::expect_identical(id,ocean_identity(root,root,settings))
  testthat::expect_identical(ocean_continuation_iterations(root,"paper",1,"O20",12000,TRUE),12000L)
  testthat::expect_identical(ocean_continuation_iterations(root,"paper",1,"O20"),12000L)
  testthat::expect_identical(ocean_continuation_label(root),"+12,000")
  testthat::expect_error(ocean_continuation_iterations(root,"paper",1,"O20",10000,TRUE),
    "differs from the saved task")
  # A separate task can choose another length; the report must display both.
  testthat::expect_identical(ocean_continuation_iterations(root,"paper",1,"O50",10000,TRUE),10000L)
  testthat::expect_identical(ocean_continuation_label(root),"+10,000 / +12,000")
  task<-data.frame(module="continuation",fold=1,cell="O20")
  testthat::expect_true(any(grepl("continuation_settings.rds",ocean_task_files(task,stage="continue"),fixed=TRUE)))
  testthat::expect_true(any(grepl("continuation_settings.rds",ocean_task_files(task,stage="analyze"),fixed=TRUE)))
  path<-file.path(root,"full5_class4/fold01/O20/convergence_parameters.csv")
  mixedgp_atomic_csv(data.frame(rhat=1.2,ess_bulk=20,ess_tail=10),path)
  decisions<-ocean_continuation_decisions(root)
  testthat::expect_identical(decisions$decision[decisions$fold==1 & decisions$cell=="O20"],"append_12000")
})
