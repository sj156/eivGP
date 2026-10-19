test_that('hybrid crosses scheduled boundaries under exact continuation', {
 set.seed(617);n<-12L;X<-matrix(rnorm(n),n);U<-matrix(rnorm(n*2),n);obs<-U;obs[5:n,]<-NA
 for(engine in c('univariate','multivariate')){
  C<-if(engine=='univariate')matrix(findInterval(U[,1],qnorm(seq(0,1,.25)),rightmost.closed=TRUE),n)else matrix(sample(1:4,n*3,TRUE),n)
  args<-list(X=X,y=sin(U[,1])+rnorm(n,sd=.2),C=C,U_obs=if(engine=='univariate')obs[,1,drop=FALSE]else obs,engine=engine,m_vec=rep(4L,ncol(C)),standardize_U=FALSE,n_chains=2L,burn=12L,parallel=FALSE,seed=67L,sampler='HYBRID')
  if(engine=='multivariate'){args$latent_dim<-2L;args$ident<-'lower_triangular'}
  for(every in c(1L,5L,20L)){
   aa<-c(args,list(sampler_control=list(hybrid_every=every)))
   a<-suppressWarnings(do.call(fit_eivgp,c(aa,list(n_iter=19L))))
   b<-continue_eivgp(a,6L,parallel=FALSE)
   full<-suppressWarnings(do.call(fit_eivgp,c(aa,list(n_iter=25L))))
   expect_identical(b$mcmc$samples_by_chain,full$mcmc$samples_by_chain)
   expect_identical(b$sampler_strategy,'collapsed_dictionary_hybrid')
   expect_identical(b$control$hybrid_every,every)
   expect_equal(b$diagnostics$hmc$retained,rep(13L,2))
   for(ch in 1:2){
    st<-b$checkpoint$states[[ch]]$state
    expect_identical(st$hybrid_extra_sweeps,as.integer(25L%/%every))
    expect_equal(unname(st$U[1:4,,drop=FALSE]),unname(args$U_obs[1:4,,drop=FALSE]))
    expect_identical(st$hmc$inverse_mass,a$checkpoint$states[[ch]]$state$hmc$inverse_mass)
    expect_identical(st$hmc$step,a$checkpoint$states[[ch]]$state$hmc$step)
   }
   changed<-a;changed$control$hybrid_every<-every+1L
   expect_error(continue_eivgp(changed,2L),'control')
  }
 }
})

test_that('hybrid controls reject invalid schedules and incompatible targets', {
 for(x in list(0L,-1L,1.5,NA_real_,Inf,c(5L,20L)))
  expect_error(mixedgp_v030_control(list(sampler='HYBRID',hybrid_every=x),4L),'hybrid_every')
 expect_equal(mixedgp_v030_control(list(sampler='HYBRID'),4L)$hybrid_every,20L)
 expect_error(mixedgp_v030_control(list(sampler='ESS',hybrid_every=20L),4L),'only used')
 expect_error(mixedgp_v030_control(list(sampler='HYBRID',dictionary_mode='marginal'),4L),'conditional')
 expect_false('hybrid_every'%in%names(mixedgp_v030_control(list(sampler='ESS'),4L)))
 expect_false('hybrid_every'%in%names(mixedgp_v030_control(list(sampler='HMC'),4L)))
})
