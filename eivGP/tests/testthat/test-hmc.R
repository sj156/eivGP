test_that('ESS and HMC preserve continuation and fixed calibration', {
 set.seed(617);n<-12L;X<-matrix(rnorm(n),n);U<-matrix(rnorm(n*2),n);obs<-U;obs[5:n,]<-NA
 for(engine in c('univariate','multivariate')){
  C<-if(engine=='univariate')matrix(findInterval(U[,1],qnorm(seq(0,1,.25)),rightmost.closed=TRUE),n)else matrix(sample(1:4,n*3,TRUE),n)
  args<-list(X=X,y=sin(U[,1])+rnorm(n,sd=.2),C=C,U_obs=if(engine=='univariate')obs[,1,drop=FALSE]else obs,engine=engine,m_vec=rep(4L,ncol(C)),standardize_U=FALSE,n_chains=2L,burn=12L,parallel=FALSE,seed=67L)
  if(engine=='multivariate'){args$latent_dim<-2L;args$ident<-'lower_triangular'}
  for(sampler in c('ESS','HMC')){
   a<-suppressWarnings(do.call(fit_eivgp,c(args,list(sampler=sampler,n_iter=20L))))
   b<-continue_eivgp(a,4L,parallel=FALSE)
   full<-suppressWarnings(do.call(fit_eivgp,c(args,list(sampler=sampler,n_iter=24L))))
   expect_identical(b$mcmc$samples_by_chain,full$mcmc$samples_by_chain)
   expect_identical(a$control$sampler,sampler)
   if(sampler=="HMC"){expect_identical(a$sampler_strategy,"collapsed_dictionary_hmc");expect_equal(b$diagnostics$hmc$retained,rep(12L,2));expect_true(all(b$diagnostics$hmc$divergence_fraction>=0))}
   if(sampler=='HMC')for(cp in b$checkpoint$states){h<-cp$state$hmc;expect_true(all(h$history[13:24,'step']>=.9*h$step & h$history[13:24,'step']<=1.1*h$step))}
   for(cp in b$checkpoint$states)expect_equal(unname(cp$state$U[1:4,,drop=FALSE]),unname(args$U_obs[1:4,,drop=FALSE]))
  }
 }
})

test_that('unsupported HMC controls are rejected explicitly', {
 expect_error(mixedgp_v030_control(list(sampler='HMC',dictionary_mode='marginal'),4L),'conditional')
 expect_error(mixedgp_v030_control(list(sampler='HMC',hmc=list(step=-1)),4L),'step')
 expect_error(mixedgp_v030_control(list(sampler='HMC',hmc=list(leapfrog=c(24,8))),4L),'ordered')
 expect_error(mixedgp_v030_control(list(sampler='HMC',hmc=list(target_accept=1)),4L),'acceptance')
})

test_that('native GP HMC target matches independent gradients and reversibility', {
 set.seed(72);n<-8L;p<-2L;d<-2L
 data<-list(X=matrix(rnorm(n*p),n),U=matrix(rnorm(n*d),n),y=rnorm(n),mis=as.integer(c(0:3,8:11)),
  lo=numeric(),hi=numeric(),logmass=numeric(),A=matrix(c(1,0,.3,.8),2,2,byrow=TRUE),
  measurement_lo=matrix(-Inf,n,2),measurement_hi=matrix(.7,n,2),txmean=rep(log(.5),p),txsd=rep(1.5,p),tu=c(.5,.7),av=3,bv=2,ar=13,br=3,threshold=FALSE,se=TRUE,nu=2.5)
 q<-c(as.numeric(data$U)[data$mis+1L],c(.1,-.2),qlogis(.8))
 density<-function(q){
  u<-data$U;u[data$mis+1L]<-q[1:8];tx<-exp(data$txmean+data$txsd*q[9:10]);r<-plogis(q[11])
  dd<-weighted_sqdist_general(data$X,u,tx,data$tu);K<-kernel_from_weighted_sqdist(dd,if(data$se)'se'else'matern',data$nu);B<-r*K+(1-r)*diag(n);L<-chol(B)
  gp<--sum(log(diag(L)))-(3+n/2)*log(2+sum(forwardsolve(t(L),data$y)^2)/2)
  mu<-u%*%t(data$A);gp+sum(log_normal_interval_prob(data$measurement_lo-mu,data$measurement_hi-mu))-sum(q[1:10]^2)/2+13*log(r)+3*log1p(-r)
 }
 for(kernel in c('se','m32','m52')){
  data$se<-kernel=='se';data$nu<-if(kernel=='m32')1.5 else 2.5
  z<-mixedgp_hmc_target_cpp(q,data);expect_equal(z$logp,density(q),tolerance=1e-10)
  grad<-vapply(seq_along(q),function(j){a<-b<-q;a[j]<-a[j]+1e-5;b[j]<-b[j]-1e-5;(density(a)-density(b))/2e-5},0)
  expect_equal(as.numeric(z$gradient),grad,tolerance=1e-6)
  mom<-rnorm(length(q));a<-mixedgp_hmc_leapfrog_cpp(q,mom,data,.001,7L,rep(1,length(q)));b<-mixedgp_hmc_leapfrog_cpp(a$q,-a$p,data,.001,7L,rep(1,length(q)))
  expect_equal(as.numeric(b$q),q,tolerance=1e-10);expect_equal(as.numeric(b$p),-mom,tolerance=1e-10)
 }
})


test_that('a fixed HMC trajectory length is honored', {
 set.seed(34);X<-matrix(rnorm(10),10);C<-matrix(rep(1:2,5),10)
 f<-fit_eivgp(X,rnorm(10),C,engine="univariate",m_vec=2L,sampler="HMC",
  sampler_control=list(hmc=list(leapfrog=c(8L,8L))),n_iter=12L,burn=4L,n_chains=2L,parallel=FALSE)
 for(cp in f$checkpoint$states)expect_true(all(cp$state$hmc$history[,"leapfrog"]==8L))
 expect_match(f$model_specification$sampler,"joint HMC")
})
