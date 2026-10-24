test_that('internal multivariate extensions preserve the uninterrupted trajectory', {
  set.seed(14);X<-matrix(rnorm(12),12,1);U<-matrix(rnorm(12),12,1)
  C<-matrix(as.integer(cut(U,c(-Inf,0,Inf))),12,1);y<-X[,1]+U[,1]+rnorm(12)
  args<-list(X_raw=X,y_raw=y,C_ord=C,U_obs=U,calib_idx=1:6,d=1L,
    m_vec=2L,ident='none',burn=10L,n_chains=4L,seed=44L,parallel_chains=FALSE)
  a<-do.call(eivGP:::fit_eivgp_ordprobit_v030,c(args,list(n_iter=25L)))
  b<-eivGP:::mixedgp_continue_study2_fit(a,10L)
  z<-do.call(eivGP:::fit_eivgp_ordprobit_v030,c(args,list(n_iter=35L)))
  expect_identical(b$mcmc$samples_by_chain,z$mcmc$samples_by_chain)
  expect_identical(b$checkpoint$states,z$checkpoint$states)
})

test_that('public ESS option defaults to the fixed budget and retains target histories', {
  set.seed(17);X<-matrix(rnorm(12),12,1);U<-matrix(rnorm(12),12,1)
  C<-matrix(as.integer(cut(U,c(-Inf,0,Inf))),12,1);y<-X[,1]+U[,1]+rnorm(12)
  args<-list(X=X,y=y,C=C,U_obs=U,latent_dim=1L,engine='multivariate',
    m_vec=2L,ident='none',burn=10L,n_iter=30L,n_chains=4L,seed=44L,parallel=FALSE)
  a<-do.call(fit_eivgp,args)
  b<-do.call(fit_eivgp,c(args,list(ess_min=NULL)))
  expect_identical(a$mcmc,b$mcmc)
  panel<-list(X=X[1:2,,drop=FALSE],C=C[1:2,,drop=FALSE],U=U[1:2,,drop=FALSE],n_latent=8L)
  z<-suppressWarnings(do.call(fit_eivgp,c(args,list(ess_min=.01,ess_panel=panel))))
  expect_true(z$scientific_ess$success)
  expect_identical(z$mcmc,a$mcmc)
  expect_true(all(z$scientific_ess$diagnostics$draw_window_complete))
  expect_true(all(lengths(z$scientific_ess$target_series[[1]])==20L))
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  prefix<-file.path(root,'scientific-trace-test')
  eivGP:::mixedgp_save_scientific_trace(z,prefix,'test')
  tr<-readRDS(paste0(prefix,'.rds'))
  expect_identical(range(tr$iteration),c(11L,30L))
  expect_length(tr$series,2L)
})

test_that('adaptive ESS budgets are bounded and invalid requests fail', {
  f<-eivGP:::mixedgp_ess_extension_size
  expect_equal(f(data.frame(ess_bulk=50,ess_tail=25),15000,100),20000L)
  expect_equal(f(data.frame(ess_bulk=99,ess_tail=120),1000,100),5000L)
  expect_equal(f(data.frame(ess_bulk=NA,ess_tail=120),1000,100),5000L)
  for(v in list(0,-1,NA,Inf,c(100,200),'100'))
    expect_error(eivGP:::mixedgp_validate_ess_request(v),'ess_min')
})
