args<-commandArgs(FALSE);here<-dirname(normalizePath(sub("^--file=","",grep("^--file=",args,value=TRUE)[1])))
source(file.path(here,"common.R"))
numerical_package_check()
api<-getFromNamespace("eivgp_numerical_api","eivGP")(NUMERICAL_REPO_ROOT)
invisible(list2env(as.list(api),envir=environment()))
cfg<-settings(list("data-root"="unused",output=tempdir(),profile="quick"))
z<-F("simulate_study2_data")(n=30,n_test=3,seed=31)
obs<-z$train$U;obs[21:30,]<-NA
m<-measurement_cut(z$train$C,obs,cfg,19)
stopifnot(!any(c("y","X","truth")%in%names(formals(measurement_cut))))
stopifnot(identical(m,measurement_cut(z$train$C,obs,cfg,19)))
for(s in m$states)stopifnot(identical(unname(s$U[1:20,]),unname(obs[1:20,])))
f<-fit_full(z$train,m$states[[1]]$U,cfg,21,TRUE)
states<-select_gp(f,1,m$states[[1]])
p<-gp_evaluate(states,z$train,z$test,cfg,22)
changed<-z$test;changed$y<-changed$y+100;changed$U<-changed$U+2
p2<-gp_evaluate(states,z$train,changed,cfg,22)
stopifnot(identical(p$m_mean,p2$m_mean),identical(p$y_mean,p2$y_mean),!identical(p$f_mean,p2$f_mean))
expected<-1/sqrt(pi)*(sqrt(2)-1)
s<-mixture_score(0,0,1);stopifnot(abs(s['crps']-expected)<1e-12,abs(s['nlpd']-log(sqrt(2*pi)))<1e-12)
stopifnot(isTRUE(all.equal(mixture_score(.3,c(0,1),c(.5,2)),mixture_score(.3,rep(c(0,1),2),rep(c(.5,2),2)),tolerance=1e-10)))
truth<-list(m=F("oracle_m0_quadrature_2d")(z$test$X,z$test$C,z$true_params),f=F("f0_2d")(z$test$X,z$test$U))
p$U<-f$mcmc$samples_U
scores<-score_result(p,z,1:20,truth,"CUT")
stopifnot(all(c("y","m","f","U")%in%scores$target),all(is.finite(scores$sq_error)))
# Real author-package adapters and prediction-after-reload for all CC competitors.
us<-prospective_u(m$states,z$test$C,cfg$latent,23)
for(method in c("UC-GP-CC","LVGP-CC","EzGP-CC")) {
  message("Checking ",method)
  path<-tempfile(fileext=".rds")
  pr<-cc_evaluate(method,z$train,z$test,1:30,us,cfg,24,path)
  pr2<-cc_evaluate(method,z$train,z$test,1:30,us,cfg,24,path)
  stopifnot(identical(pr,pr2),all(is.finite(pr$y_mean)),all(pr$y_var>=0));unlink(path)
  cat_path<-file.path(tempdir(),paste0("categorical-",sub("-CC$","",method)),"predictions.rds")
  cat_pred<-categorical_evaluate(sub("-CC$","",method),z$train,z$test,cfg,24,cat_path)
  stopifnot(all(is.finite(cat_pred$m_point)),all(cat_pred$y_var>=0),
    identical(cat_pred,categorical_evaluate(sub("-CC$","",method),z$train,z$test,cfg,24,cat_path)))
  unlink(dirname(cat_path),recursive=TRUE)
}
cat("PASS: response-free determinism, fixed calibration, no test truth leakage, mixture scores, all targets, categorical-only and CC adapters, and fit reload.\n")

# Compact fit roundtrip preserves predictions and resumable chain history.
fp<-tempfile();sp<-tempfile();atomic_rds(f,fp)
stopifnot(inherits(readRDS(fp),"mixedgp_compact_fit"),
          identical(f,F("mixedgp_read_artifact")(fp)))
atomic_rds(list(pass=FALSE,states=states),sp)
stopifnot(identical(retire_conditional_fit(fp,sp)$states,states),file.exists(fp))
atomic_rds(list(pass=TRUE,states=list()),sp)
invisible(retire_conditional_fit(fp,sp));stopifnot(file.exists(fp))
atomic_rds(list(pass=TRUE,states=states),sp)
retained<-retire_conditional_fit(fp,sp)
stopifnot(!file.exists(fp),identical(retained$states,states),
          identical(p2,gp_evaluate(retained$states,z$train,changed,cfg,22)))
unlink(sp)
cat("PASS: lossless compact fit roundtrip, failed-fit retention, and passing-fit retirement with exact prediction reconstruction.\n")
