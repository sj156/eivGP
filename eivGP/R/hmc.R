# Joint conditional HMC: all missing inputs, numeric-kernel coefficients and r.
# Measurement parameters and dictionary states retain their existing updates.
mixedgp_hmc_coordinates <- function(state,ctx) {
 d<-mixedgp_v030_decode(state,ctx);u<-as.numeric(state$U[ctx$missing,,drop=FALSE])
 if(ctx$measurement=='threshold'){
  cl<-ctx$C[ctx$missing,1];u<-mixedgp_v030_interval_inverse(u,d$tau[[1]][cl],d$tau[[1]][cl+1])
 }
 c(u,(state$logtheta_x-ctx$priors$log_theta_x_mean)/ctx$priors$log_theta_x_sd,qlogis(d$r))
}
mixedgp_hmc_data <- function(state,ctx) {
 d<-mixedgp_v030_decode(state,ctx);p<-ctx$priors;cl<-ctx$C[ctx$missing,1]
 threshold<-ctx$measurement=='threshold';lo<-hi<-lm<-numeric()
 if(threshold){lo<-d$tau[[1]][cl];hi<-d$tau[[1]][cl+1];lm<-log_normal_interval_prob(lo,hi)}
 list(X=ctx$X,U=state$U,y=ctx$y,mis=as.integer(as.vector(outer(ctx$missing-1L,(seq_len(ctx$d)-1L)*ctx$n,'+'))),
  lo=lo,hi=hi,logmass=lm,A=d$A,
  measurement_lo=matrix(unlist(lapply(seq_len(ctx$q),function(j)d$tau[[j]][ctx$C[,j]])),ctx$n),
  measurement_hi=matrix(unlist(lapply(seq_len(ctx$q),function(j)d$tau[[j]][ctx$C[,j]+1L])),ctx$n),
  txmean=rep_len(p$log_theta_x_mean,ctx$p),txsd=rep_len(p$log_theta_x_sd,ctx$p),tu=d$theta_u,
  av=p$variance_shape,bv=p$variance_rate,ar=p$signal_shape[1],br=p$signal_shape[2],
  threshold=threshold,se=ctx$kernel=='se',nu=ctx$matern_nu)
}
mixedgp_hmc_unpack <- function(q,state,ctx) {
 m<-length(ctx$missing)*ctx$d;u<-q[seq_len(m)]
 if(ctx$measurement=='threshold'){
  d<-mixedgp_v030_decode(state,ctx);cl<-ctx$C[ctx$missing,1]
  u<-mixedgp_v030_interval_forward(u,d$tau[[1]][cl],d$tau[[1]][cl+1])
 }
 state$U[ctx$missing,]<-matrix(u,length(ctx$missing),ctx$d)
 state$logtheta_x<-ctx$priors$log_theta_x_mean+ctx$priors$log_theta_x_sd*q[m+seq_len(ctx$p)]
 r<-plogis(q[m+ctx$p+1L]);sh<-ctx$priors$signal_shape
 lp<-pbeta(r,sh[1],sh[2],log.p=TRUE);lq<-pbeta(r,sh[1],sh[2],lower.tail=FALSE,log.p=TRUE)
 state$e_r<-if(lp<lq)qnorm(lp,log.p=TRUE)else qnorm(lq,lower.tail=FALSE,log.p=TRUE)
 state
}
mixedgp_hmc_da_init <- function(step)list(t=0L,hbar=0,mu=log(10*step),logbar=log(step))
mixedgp_hmc_joint <- function(state,ctx,iteration,burn) {
 settings<-ctx$control$hmc;q<-mixedgp_hmc_coordinates(state,ctx);d<-length(q);h<-state$hmc
 if(is.null(h))h<-list(step=settings$step,inverse_mass=rep(1,d),da=mixedgp_hmc_da_init(settings$step),wn=0L,mean=rep(0,d),m2=rep(0,d),
   history=matrix(numeric(),0,7,dimnames=list(NULL,c('iteration','accept_prob','accepted','delta_h','max_error','step','leapfrog'))))
 step<-h$step*runif(1,.9,1.1);choices<-seq.int(settings$leapfrog[1],settings$leapfrog[2]);L<-choices[sample.int(length(choices),1L)]
 proposal<-mixedgp_hmc_leapfrog_cpp(q,rnorm(d)/sqrt(h$inverse_mass),mixedgp_hmc_data(state,ctx),step,L,h$inverse_mass)
 ap<-if(is.finite(proposal$delta_h))exp(min(0,-proposal$delta_h))else 0;accept<-runif(1)<ap
 if(accept){state<-mixedgp_hmc_unpack(proposal$q,state,ctx);q<-proposal$q}
 h$history<-rbind(h$history,c(iteration,ap,accept,proposal$delta_h,proposal$max_error,step,L))
 if(iteration<=burn){
  h$wn<-h$wn+1L;delta<-q-h$mean;h$mean<-h$mean+delta/h$wn;h$m2<-h$m2+delta*(q-h$mean)
  a<-h$da;a$t<-a$t+1L;t<-a$t;eta<-1/(t+10);a$hbar<-(1-eta)*a$hbar+eta*(settings$target_accept-ap)
  ls<-a$mu-sqrt(t)/.05*a$hbar;a$logbar<-t^(-.75)*ls+(1-t^(-.75))*a$logbar;h$da<-a;h$step<-exp(max(-12,min(1,ls)))
  if(iteration %in% unique(as.integer(burn*c(.25,.5))) && iteration<burn){
   h$inverse_mass<-pmax(.001,pmin(100,(h$wn*h$m2/max(1,h$wn-1L)+5)/(h$wn+5)))
   h$step<-exp(h$da$logbar);h$da<-mixedgp_hmc_da_init(h$step);h$wn<-0L;h$mean[]<-0;h$m2[]<-0
  }
  if(iteration==burn)h$step<-exp(h$da$logbar)
 }
 state$hmc<-h;state
}


# Diagnostic-only recovery; does not draw randomness or modify sampler states.
mixedgp_hmc_attach_diagnostics <- function(fit) {
 if(!identical(fit$control$sampler,"HMC"))return(fit)
 tab<-do.call(rbind,lapply(seq_along(fit$checkpoint$states),function(ch){
  h<-fit$checkpoint$states[[ch]]$state$hmc
  v<-h$history[h$history[,"iteration"]>fit$checkpoint$burn,,drop=FALSE]
  divergent<-!is.finite(v[,"max_error"]) | v[,"max_error"]>1000
  data.frame(chain=ch,retained=nrow(v),mean_acceptance=mean(v[,"accept_prob"]),
   divergences=sum(divergent),divergence_fraction=mean(divergent),step_size=h$step)
 }))
 fit$diagnostics$hmc<-tab
 fit$diagnostics$summary$hmc_divergences<-sum(tab$divergences)
 fit$diagnostics$summary$hmc_max_chain_divergence_fraction<-max(tab$divergence_fraction)
 fit
}
