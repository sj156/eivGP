repo<-Sys.getenv('EIVGP_TEST_REPO');stopifnot(nzchar(repo))
source(file.path(repo,'tools/ocean-tools.R'),local=TRUE)
test_that('Ocean ESS settings are validated',{
  expect_null(ocean_settings()$ess_min)
  expect_equal(ocean_settings(ess_min='100')$ess_min,100)
  expect_null(ocean_settings(ess_min='NULL')$ess_min)
  expect_error(ocean_settings(profile='quick',cores=1L,ess_min=100),'four chains')
  expect_equal(ocean_settings(profile='quick',cores=4L,ess_min=100)$ess_min,100)
})
test_that('Ocean fit phase attaches scientific diagnostics without changing a sufficient fit', {
  set.seed(29);X<-matrix(rnorm(24),12,2);u<-sort(rnorm(12));C<-rep(1:4,each=3)
  fit<-fit_eivgp(X,X[,1]+u+rnorm(12),C,U_obs=matrix(u,ncol=1),engine='univariate',m_vec=4L,
    n_iter=24L,burn=10L,n_chains=4L,parallel=FALSE,seed=8L)
  root<-tempfile();dir.create(root);dir.create(file.path(root,'O50'));on.exit(unlink(root,recursive=TRUE))
  saveRDS(fit,file.path(root,'O50','fit.rds'))
  env<-new.env(parent=environment(ocean_initial_functions));ocean_initial_functions(env)
  env$stage<-'fit';env$ess_min<-.01;env$cores<-1L;env$smoke<-FALSE
  env$atomic_rds<-mixedgp_atomic_rds;env$atomic_csv<-mixedgp_atomic_csv
  tr<-data.frame(G2salinity=X[,1],log_pressure=X[,2],Y=fit$data$y_raw,C=C,U=u,calib_o50=1L)
  suppressWarnings(env$run_eiv(1L,'O50',tr,tr,root))
  saved<-readRDS(file.path(root,'O50','fit.rds'))
  expect_true(saved$scientific_ess$success)
  expect_identical(saved$mcmc,fit$mcmc)
  expect_equal(read.csv(file.path(root,'O50','convergence_summary.csv'))$iterations,24L)
  expect_true(file.exists(file.path(root,'O50','scientific-trace-targets.rds')))
})

test_that('CPU selection preserves the Ocean chain count',{
  # CPU resolution never modifies the independently selected four-chain setting.
  cfg <- ocean_settings(profile="paper",cores=56L)
  expect_equal(cfg$chains,4L)
})
