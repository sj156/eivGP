test_that('application ESS settings are opt-in and validated', {
  expect_null(ocean_settings()$ess_min)
  expect_equal(ocean_settings(ess_min='100')$ess_min,100)
  expect_null(ocean_settings(ess_min='NULL')$ess_min)
  expect_error(ocean_settings(profile='quick',ess_min=100),'four chains')
  expect_null(adni_settings()$ess_min)
  expect_equal(adni_settings(ess_min='100')$ess_min,100)
  expect_equal(adni_cli_params('--ess-min=100')$ess_min,'100')
  expect_error(adni_settings(ess_min='nonsense'),'ess_min')
  m<-list(stages=c(100L,200L),burn=20L,folds=1L,cells='O50',ess_min=100)
  expect_null(adni_screen_settings(m,100L,1L)$ess_min)
  expect_equal(adni_screen_settings(m,200L,1L)$ess_min,100)
  n<-mixedgp_application_numerics(data.frame(iterations=c(50000L,60000L),burn=5000L,chains=4L))
  expect_equal(n$iterations,50000L);expect_equal(n$iterations_max,60000L)
  expect_error(mixedgp_application_numerics(data.frame(iterations=c(50,60),burn=c(5,6),chains=4)),'warmup')
})

test_that('CUT conditional targets freeze outer measurement parameters', {
  set.seed(19);X<-matrix(rnorm(24),12,2);U<-matrix(rnorm(24),12,2)
  C<-matrix(sample(1:4,48,TRUE),12,4);y<-X[,1]+U[,1]+rnorm(12)
  fit<-fit_eivgp(X,y,C,U_obs=U,engine='multivariate',latent_dim=2,
    m_vec=rep(4L,4),ident='lower_triangular',standardize_U=FALSE,
    n_iter=22,burn=10,n_chains=4,parallel=FALSE,seed=7)
  measurement<-list(U=U,A=fit$mcmc$samples_A[1,,],tau=unflatten_tau(fit$mcmc$samples_tau[1,],rep(4L,4)))
  copy<-fit
  copy$mcmc$samples_by_chain$A<-lapply(copy$mcmc$samples_by_chain$A,function(a)a+3)
  a<-mixedgp_cut_conditional_series(fit,measurement,X[1:2,],C[1:2,],U[1:2,],4L,99L)
  b<-mixedgp_cut_conditional_series(copy,measurement,X[1:2,],C[1:2,],U[1:2,],4L,99L)
  expect_identical(a,b)
  expect_true(all(lengths(a[[1]])==12L))
  expect_false(identical(copy$mcmc$samples_by_chain$A,fit$mcmc$samples_by_chain$A))
  # Actual ADNI task routes a saved fit through final-stage science checks.
  out<-tempfile();dir.create(out);on.exit(unlink(out,recursive=TRUE))
  task<-data.frame(fold=1L,cell='O50')
  path<-adni_plan_fit_path(out,task);adni_plan_save_fit(fit,path)
  te<-data.frame(X1=X[,1],X2=X[,2],C,U1=U[,1],U2=U[,2]);names(te)[3:6]<-paste0('C',1:4)
  z<-list(test=te,x=c('X1','X2'),c=paste0('C',1:4),u=c('U1','U2'),protocol=list(seeds=list(prediction=11L)))
  s<-list(iter=22L,burn=10L,chains=4L,cores=1L,final_stage=TRUE,ess_min=.01)
  suppressWarnings(adni_plan_fit_task(task,list(cells=list('1/O50'=z)),out,s))
  saved<-adni_plan_read_fit(path)
  expect_true(saved$scientific_ess$success)
  expect_identical(saved$mcmc,fit$mcmc)
  expect_true(file.exists(file.path(dirname(path),'scientific-trace-targets.rds')))
})

test_that('CUT outer diagnostics retain chronology and never fit or extend', {
  outer_fun<-mixedgp_cut_outer_trace
  env<-new.env(parent=environment(outer_fun));env$continue_eivgp<-function(...)stop('must not extend')
  env$fit_eivgp<-function(...)stop('must not fit');environment(outer_fun)<-env
  states<-list()
  for(it in c(20L,40L))for(ch in 1:4)for(gp in 1:2)
    states[[length(states)+1]]<-list(chain=ch,iteration=it)
  N<-6L;S<-length(states)
  pred<-list(m_mean=matrix(seq_len(N*S)/100,N,S),f_mean=matrix(seq_len(N*S)/50,N,S),
    y_mean=matrix(seq_len(N*S*2)/100,N,S*2),y_var=matrix(.5,N,S*2))
  test<-list(X=matrix(seq_len(N*2),N,2),C=matrix(rep(1:3,length.out=N*4),N,4))
  out<-tempfile();dir.create(out);on.exit(unlink(out,recursive=TRUE))
  outer_fun(pred,states,test,list(gp_response_latent=2L),out,'CUT test')
  z<-readRDS(file.path(out,'scientific-trace-cut-outer.rds'))
  expect_identical(z$iteration,c(20L,40L))
  expect_true(all(z$table$scope=='outer_integration'))
  expect_true(all(lengths(z$series[[1]])==2L))
  expect_true(all(is.na(z$table$ess_bulk))) # Too short; never invented or passed.
})

test_that('Ocean fit phase attaches scientific diagnostics without changing a sufficient fit', {
  set.seed(29);X<-matrix(rnorm(24),12,2);u<-sort(rnorm(12));C<-rep(1:4,each=3)
  fit<-fit_eivgp(X,X[,1]+u+rnorm(12),C,U_obs=matrix(u,ncol=1),engine='univariate',m_vec=4L,
    n_iter=24L,burn=10L,n_chains=4L,parallel=FALSE,seed=8L)
  root<-tempfile();dir.create(root);dir.create(file.path(root,'O50'));on.exit(unlink(root,recursive=TRUE))
  saveRDS(fit,file.path(root,'O50','fit.rds'))
  env<-new.env(parent=asNamespace('eivGP'));ocean_initial_functions(env)
  env$stage<-'fit';env$ess_min<-.01;env$cores<-1L;env$smoke<-FALSE
  env$atomic_rds<-ocean_atomic_rds;env$atomic_csv<-ocean_atomic_csv
  tr<-data.frame(G2salinity=X[,1],log_pressure=X[,2],Y=fit$data$y_raw,C=C,U=u,calib_o50=1L)
  suppressWarnings(env$run_eiv(1L,'O50',tr,tr,root))
  saved<-readRDS(file.path(root,'O50','fit.rds'))
  expect_true(saved$scientific_ess$success)
  expect_identical(saved$mcmc,fit$mcmc)
  expect_equal(read.csv(file.path(root,'O50','convergence_summary.csv'))$iterations,24L)
  expect_true(file.exists(file.path(root,'O50','scientific-trace-targets.rds')))
})
