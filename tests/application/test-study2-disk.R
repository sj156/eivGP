# Exact Study II adapter regression: no study result files are used or changed.
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1L])
repo<-dirname(dirname(dirname(normalizePath(script))))
lib<-Sys.getenv('EIVGP_R_LIBRARY',file.path(repo,'.r-library'))
.libPaths(c(lib,file.path(repo,'.r-library'),.libPaths()))
if(utils::packageVersion('eivGP')<'0.4.0.9020')stop('Install/use eivGP >=0.4.0.9020')
source(file.path(repo,'tools/study2-disk-tools.R'));F<-function(n)getFromNamespace(n,'eivGP')
root<-tempfile('disk-adapter-');dir.create(root)
tryCatch({
 set.seed(819);n<-12L;X<-matrix(rnorm(n*2),n);C<-cbind(rep(1:3,4),rep(3:1,4),rep(c(1,3,2),4));y<-rnorm(n);U<-matrix(rnorm(n*2),n)
 for(scores in c(FALSE,TRUE))for(cal in c(0L,4L))for(method in c('ESS','HMC','HYBRID')) {
  obs<-matrix(NA_real_,n,2);if(cal>0)obs[1:cal,]<-U[1:cal,]
  fit<-suppressWarnings(eivGP::fit_eivgp(X,y,C,U_obs=obs,engine='multivariate',latent_dim=2L,m_vec=rep(3L,3),
   ident='lower_triangular',standardize_U=FALSE,store_scores=scores,sampler=method,
   sampler_control=if(method=='HYBRID')list(hybrid_every=3L)else list(),
   n_iter=16L,burn=4L,n_chains=4L,seed=611L,parallel=FALSE))
  original_public<-fit
  class(fit)<-NULL;fit$interface<-NULL;fit$model_specification<-NULL;fit$call<-NULL
  fit$checkpoint$arguments<-NULL;fit$checkpoint$signature<-NULL
  expected<-if(scores)eivGP::continue_eivgp(original_public,5L,parallel=FALSE,n_cores=1L)else F('mixedgp_continue_study2_fit')(fit,5L,FALSE,1L)
  public<-study2_disk_public(fit)
  path<-file.path(root,paste0(method,'-',cal,'-',scores,'.rds'));study2_disk_atomic(F('mixedgp_pack_artifact')(public),path)
  callback<-function(z)stopifnot(is.null(z$mcmc$samples_by_chain),all(vapply(z$checkpoint$states,function(x)is.null(x$state$hmc$history)||nrow(x$state$hmc$history)==0L,logical(1))))
  z<-eivGP::continue_eivgp_disk(path,5L,chunk_iter=2L,parallel=FALSE,n_cores=1L,progress=callback)
  actual<-eivGP::read_fit(path)
  stopifnot(identical(actual$mcmc$samples_by_chain,expected$mcmc$samples_by_chain),
   identical(actual$checkpoint$states,expected$checkpoint$states),
   identical(actual$priors,fit$priors),identical(actual$control,fit$control),identical(actual$data,fit$data),
   all(z$disk_history$blocks$to-z$disk_history$blocks$from+1L<=2L))
  cat('Exact Study II internal-to-disk equivalence:',method,'calibration',cal,'scores',scores,'\n')
 }
},finally=unlink(root,recursive=TRUE))
