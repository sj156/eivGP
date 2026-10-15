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
  env<-new.env(parent=asNamespace("eivGP"))
  env$root<-root;env$out<-output;env$continuation_iter<-80L
  env$atomic_csv<-ocean_atomic_csv
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
