test_that('application ESS settings are opt-in and validated', {
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
