#!/usr/bin/env Rscript
# Locate the shared bootstrap beside this script, including when source() is used.
.local_file <- if(sys.nframe()>0L && !is.null(sys.frame(1)$ofile))sys.frame(1)$ofile else sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(normalizePath(.local_file)),"common.R"))
# Targeted Study II additions. All functions are sourceable for testing.
# CUT targets q(U,psi | C,U_obs) p(theta | X,y,U); no response reweighting.
F <- function(name) getFromNamespace(name, "eivGP")
atomic_rds <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp<-tempfile(tmpdir=dirname(path));on.exit(unlink(tmp));saveRDS(F("mixedgp_pack_artifact")(x),tmp)
  if(!file.rename(tmp,path)) stop("Cannot publish checkpoint: ",path)
}
retire_conditional_fit <- function(path, states_path) {
  saved<-readRDS(states_path)
  if(isTRUE(saved$pass) && length(saved$states)>0L && file.exists(path))unlink(path)
  saved
}
hash_object <- function(x) {
  p<-tempfile();on.exit(unlink(p));saveRDS(x,p,version=2);unname(tools::md5sum(p))
}
parse_cli <- function(args) {
  out<-list(); for(a in args) {
    if(!grepl("^--[^=]+=.+$",a)) stop("Use --name=value: ",a)
    k<-sub("^--([^=]+)=.*$","\\1",a);out[[k]]<-sub("^--[^=]+=","",a)
  };out
}
integers <- function(x) {
  z<-unlist(lapply(strsplit(x,",",fixed=TRUE)[[1]],function(a) {
    b<-as.integer(strsplit(a,":",fixed=TRUE)[[1]]);if(length(b)==2)seq.int(b[1],b[2]) else b
  }));if(anyNA(z)||any(z<1)||anyDuplicated(z))stop("Invalid integer selection: ",x);z
}
settings <- function(a) {
  allowed<-c("action","profile","data-root","output","reps","calib","methods","outer","inner-draws",
    "latent","iter","warmup","inner-iter","inner-warmup","cores","seed","full-fit-map")
  if(length(setdiff(names(a),allowed)))stop("Unknown option: ",paste(setdiff(names(a),allowed),collapse=","))
  get<-function(k,default)if(is.null(a[[k]]))default else a[[k]]
  quick<-identical(get("profile","paper"),"quick")
  if(!get("profile","paper")%in%c("paper","quick"))stop("profile must be paper or quick")
  cfg<-list(full_fit_map=get("full-fit-map",""),profile=get("profile","paper"),data_root=get("data-root",Sys.getenv("EIVGP_DATA_ROOT",numerical_default_data("study2"))),
    output=get("output",Sys.getenv("EIVGP_SUPPLEMENT_ROOT",numerical_default_output("supplement"))),reps=integers(get("reps",if(quick)"1"else"1:50")),
    calib=integers(get("calib","20,50")),methods=strsplit(get("methods","FULL,CUT,UC-GP,LVGP,EzGP,UC-GP-CC,LVGP-CC,EzGP-CC"),",",fixed=TRUE)[[1]])
  defaults<-if(quick)c(outer=4,inner_draws=1,latent=4,iter=40,warmup=20,inner_iter=40,inner_warmup=20,cores=1,seed=77101) else
    c(outer=64,inner_draws=4,latent=64,iter=20000,warmup=5000,inner_iter=4000,inner_warmup=1000,cores=4,seed=77101)
  for(k in names(defaults))cfg[[k]]<-as.integer(get(gsub("_","-",k),defaults[[k]]))
  if(any(!is.finite(unlist(cfg[names(defaults)])))||any(unlist(cfg[names(defaults)])<1))stop("Numeric controls must be positive")
  if(cfg$iter<=cfg$warmup||cfg$inner_iter<=cfg$inner_warmup||cfg$outer%%4!=0)stop("iterations must exceed warmup; outer must be divisible by four")
  if(cfg$latent<2||cfg$inner_draws>cfg$inner_iter-cfg$inner_warmup)stop("Need at least two latent draws and enough retained conditional draws")
  if(cfg$outer/4>cfg$iter-cfg$warmup)stop("Too few measurement draws for outer selection")
  if(!all(cfg$calib%in%c(20,50)))stop("Minimal supplement uses calibration 20 and 50 only")
  if(!all(cfg$methods%in%c("CUT","FULL","UC-GP","LVGP","EzGP","UC-GP-CC","LVGP-CC","EzGP-CC")))stop("Unknown method")
  if(!nzchar(cfg$output)||!nzchar(cfg$data_root))stop("Supply --data-root and --output")
  if(grepl("Dropbox|Overleaf",cfg$output,ignore.case=TRUE))stop("Keep output outside the synced manuscript tree")
  cfg
}
# This API intentionally does not accept y, X, or evaluation truth.
measurement_cut <- function(C,U_obs,cfg,seed) {
  run<-function(chain) {
    set.seed(seed+10000L*chain)
    ctx<-list(C=C,U_obs=U_obs,n=nrow(C),p=2L,q=ncol(C),d=ncol(U_obs),m_vec=rep(4L,ncol(C)),
      observed=which(complete.cases(U_obs)),missing=which(!complete.cases(U_obs)),measurement="probit",ident="lower_triangular")
    ctx$priors<-F("mixedgp_v030_priors")(list(),ctx$p,ctx$d,ctx$m_vec)
    st<-F("mixedgp_v030_initial_state")(ctx)
    logmeas<-function(s)F("mixedgp_v030_measurement_loglik")(s,F("mixedgp_v030_decode")(s,ctx),ctx)
    K<-cfg$iter-cfg$warmup;states<-vector("list",cfg$outer/4)
    selected<-unique(round(seq(1,K,length.out=cfg$outer/4)))
    trace<-matrix(NA_real_,K,length(ctx$missing)*ctx$d+ctx$q*ctx$d+sum(ctx$m_vec))
    for(it in seq_len(cfg$iter)) {
      for(rows in split(ctx$missing,ceiling(seq_along(ctx$missing)/8))) {
        old<-as.vector(st$U[rows,,drop=FALSE]);base<-st
        proposal<-F("mixedgp_v030_ess")(old,function(v){s<-base;s$U[rows,]<-v;logmeas(s)})
        st$U[rows,]<-proposal$value
      }
      for(j in seq_len(ctx$q)) {
        active<-seq_len(min(j,ctx$d));np<-length(st$e_pi[[j]]);base<-st
        update<-function(v){s<-base;s$e_pi[[j]]<-v[seq_len(np)];s$e_A[j,active]<-v[np+seq_along(active)];s}
        v<-F("mixedgp_v030_ess")(c(st$e_pi[[j]],st$e_A[j,active]),function(v)logmeas(update(v)))$value
        st<-update(v)
      }
      if(it>cfg$warmup) {
        k<-it-cfg$warmup;dec<-F("mixedgp_v030_decode")(st,ctx)
        trace[k,]<-c(st$U[ctx$missing,,drop=FALSE],dec$A,unlist(dec$pi))
        if(k%in%selected)states[[match(k,selected)]]<-list(U=st$U,A=dec$A,tau=lapply(dec$tau,function(t)t[-c(1,length(t))]),chain=chain,iteration=it)
      }
    }
    list(states=states,trace=trace)
  }
  z<-if(.Platform$OS.type=="unix"&&cfg$cores>1)parallel::mclapply(1:4,run,mc.cores=min(4,cfg$cores),mc.set.seed=FALSE)else lapply(1:4,run)
  if(any(vapply(z,inherits,logical(1),"try-error")))stop("Measurement chain failed")
  mis<-which(!complete.cases(U_obs))
  labels<-c(unlist(lapply(1:2,function(k)sprintf("U[%d,%d]",mis,k))),
    unlist(lapply(1:2,function(k)sprintf("A[%d,%d]",seq_len(ncol(C)),k))),
    unlist(lapply(seq_len(ncol(C)),function(j)sprintf("pi[%d,%d]",j,1:4))))
  arr<-array(NA_real_,c(nrow(z[[1]]$trace),4,ncol(z[[1]]$trace)),
    dimnames=list(NULL,paste0("chain",1:4),labels))
  for(j in 1:4)arr[,j,]<-z[[j]]$trace
  # Constant constrained loading entries have no diagnostic interpretation.
  keep<-apply(arr,3,function(x)stats::sd(as.vector(x))>0)
  diagnostics<-posterior::summarise_draws(posterior::as_draws_array(arr[,,keep,drop=FALSE]))
  # Interleave chains, so first/second-half integration checks include every chain.
  states<-unlist(lapply(seq_len(cfg$outer/4),function(i)lapply(z,function(c)c$states[[i]])),recursive=FALSE)
  list(states=states,diagnostics=diagnostics,trace=arr,
    target="p(U,psi | C,U_obs); same category-probability and loading priors as FULL")
}
diagnostic_pass <- function(z) {
  nrow(z)>0 && all(is.finite(z$rhat) & z$rhat<=1.01 &
    is.finite(z$ess_bulk) & z$ess_bulk>=400 & is.finite(z$ess_tail) & z$ess_tail>=400)
}
record_gp_diagnostics <- function(fit,path,conditional=FALSE) {
  z<-fit$diagnostics$ess_key
  if(conditional)z<-z[grepl("^(V$|r$|theta_x|J\\[|dictionary_occupancy)",z$parameter),,drop=FALSE]
  else z<-z[is.finite(z$posterior_sd)&z$posterior_sd>0,,drop=FALSE]
  write.csv(z,path,row.names=FALSE)
  diagnostic_pass(z)
}
fit_full <- function(train,U_obs,cfg,seed,conditional=FALSE) {
  eivGP::fit_eivgp(train$X,train$y,train$C,U_obs=U_obs,latent_dim=2L,engine="multivariate",
    m_vec=rep(4L,4),ident="lower_triangular",standardize_U=FALSE,n_chains=4L,
    n_iter=if(conditional)cfg$inner_iter else cfg$iter,burn=if(conditional)cfg$inner_warmup else cfg$warmup,
    seed=seed,parallel=cfg$cores>1,n_cores=cfg$cores,verbose=FALSE)
}
read_matched_full <- function(path,train,U_obs) {
  fit<-readRDS(path)
  if(is.null(fit$mcmc)&&!is.null(fit$fit))fit<-fit$fit
  equal<-function(a,b)isTRUE(all.equal(unname(a),unname(b),check.attributes=FALSE,tolerance=1e-12))
  expected<-F("mixedgp_v030_priors")(list(),ncol(train$X),2L,rep(4L,4))
  if(!identical(fit$sampler_version,"0.3.1")||!equal(fit$priors,expected)||
     !equal(fit$data$X_raw,train$X)||!equal(fit$data$y_raw,train$y)||
     !equal(fit$data$C_ord,train$C)||!equal(fit$data$U_obs,U_obs)||
     !identical(fit$data$ident,"lower_triangular")||
     !all(fit$data$U_center==0)||!all(fit$data$U_scale==1)||
     !identical(fit$kernel$name,"se"))stop("Saved FULL fit does not match data, calibration, priors or model: ",path)
  fit
}
select_gp <- function(fit,number,measurement=NULL) {
  mm<-fit$mcmc;info<-mm$mcmc_draw_info
  ids<-unlist(lapply(split(seq_len(nrow(info)),info$chain),function(i)i[unique(round(seq(1,length(i),length.out=number)))]))
  lapply(ids,function(i) {
    st<-if(is.null(measurement))list(U=mm$samples_U[i,,],A=mm$samples_A[i,,],
      tau=F("unflatten_tau")(mm$samples_tau[i,],rep(4L,4)))else measurement
    c(st,list(logtheta=mm$samples_logtheta[i,],sigma2=mm$samples_sigma2[i],gp_chain=info$chain[i]))
  })
}
# Same integration points for all CC methods. Only C_test enters this function.
prospective_u <- function(states,C_test,L,seed) {
  set.seed(seed)
  lapply(states,function(s) {
    a<-F("sample_u_given_c_ordprobit_dispatch")(C_test[rep(seq_len(nrow(C_test)),each=L),,drop=FALSE],
      s$A,s$tau,rep(4L,4),latent_sampler="minimax_tilting")
    array(as.numeric(a),c(L,nrow(C_test),2L))
  })
}
# Components of the posterior distribution of m and f, plus predictive y.
# The covariance of the integrated GP is retained: mean(cov), not mean(diag(cov)).
gp_evaluate <- function(states,train,test,cfg,seed) {
  xy<-F("mixedgp_v030_standardize_xy")(train$X,train$y)
  Xt<-sweep(sweep(test$X,2,xy$X_center,"-"),2,xy$X_scale,"/")
  us<-prospective_u(states,test$C,cfg$latent,seed)
  N<-nrow(test$X);S<-length(states)
  out<-lapply(c("y_mean","y_var","m_mean","m_var","f_mean","f_var"),function(x)matrix(NA_real_,N,S))
  names(out)<-c("y_mean","y_var","m_mean","m_var","f_mean","f_var")
  for(s in seq_len(S)) {
    st<-states[[s]]
    predict<-function(x,u,cov=FALSE,noisy=FALSE)F("gp_predict_draw_general")(xy$X,st$U,xy$y,x,u,
      st$logtheta,st$sigma2,noisy=noisy,return_cov=cov)
    f<-predict(Xt,test$U);out$f_mean[,s]<-xy$y_center+xy$y_scale*f$mean;out$f_var[,s]<-xy$y_scale^2*f$var
    # One independent prospective U per posterior component for the y mixture.
    yp<-predict(Xt,matrix(us[[s]][1,,],N,2),noisy=TRUE)
    out$y_mean[,s]<-xy$y_center+xy$y_scale*yp$mean;out$y_var[,s]<-xy$y_scale^2*yp$var
    for(i in seq_len(N)) {
      mp<-predict(Xt[rep(i,cfg$latent),,drop=FALSE],matrix(us[[s]][,i,],cfg$latent,2),cov=TRUE)
      out$m_mean[i,s]<-xy$y_center+xy$y_scale*mean(mp$mean)
      out$m_var[i,s]<-xy$y_scale^2*max(0,mean(mp$cov))
    }
  }
  out
}
categorical_evaluate <- function(method,train,test,cfg,seed,path) {
  if(file.exists(path))return(readRDS(path))
  controls<-if(cfg$profile=="quick")switch(method,`UC-GP`=list(n_starts=2L),
    LVGP=list(n_starts=1L,max_retries=1L,max_iter_ini=10L,max_iter_lat=5L),
    EzGP=list(tau_fractions=c(.01,.04),cv_folds=2L,maxeval=10L))else list()
  model<-do.call(eivGP::fit_mixedgp_competitor,c(list(method=method,X=train$X,y=train$y,C=train$C,
    new_X=test$X,new_C=test$C,m_vec=rep(4L,4),n_draw=1L,seed=seed),controls))
  atomic_rds(model,sub("predictions.rds$","fit.rds",path))
  out<-list(y_mean=matrix(model$predictive_mean,ncol=1),y_var=matrix(model$predictive_variance,ncol=1),
    m_point=as.numeric(model$predictive_mean),optimization_status=model$optimization_status)
  atomic_rds(out,path);out
}
cc_predict <- function(method,model,X,C) {
  if(method=="UC-GP") {
    nd<-F("make_mixedgp_competitor_data")(X,C,rep(4L,4))
    p<-stats::predict(model$fit,newdata=nd,type="UK",seCompute=TRUE,covCompute=FALSE,lightReturn=TRUE)
    mu<-as.numeric(p$mean);v<-as.numeric(p$sd)^2
  } else {
    p<-if(method=="LVGP")LVGP::LVGP_predict(cbind(X,C),model$fit,MSE_on=1) else EzGP::EzGP_predict(cbind(X,C),model$fit,MSE_on=1)
    mu<-as.numeric(p$Y_hat);v<-F("predictive_variance_diagonal")(p$MSE,nrow(X),context=method)
  }
  list(mean=mu,var=v+model$fitted_noise_variance)
}
cc_evaluate <- function(method,train,test,cal,us,cfg,seed,path) {
  base<-sub("-CC$","",method)
  if(file.exists(path))model<-readRDS(path)else {
    controls<-if(cfg$profile=="quick")switch(base,`UC-GP`=list(n_starts=2L),LVGP=list(n_starts=1L,max_retries=1L,max_iter_ini=10L,max_iter_lat=5L),EzGP=list(tau_fractions=c(.01,.04),cv_folds=2L,maxeval=10L))else list()
    model<-do.call(eivGP::fit_mixedgp_competitor,c(list(method=base,X=cbind(train$X[cal,,drop=FALSE],train$U[cal,,drop=FALSE]),
      y=train$y[cal],C=train$C[cal,,drop=FALSE],new_X=cbind(test$X[rep(1L,2L),,drop=FALSE],matrix(rep(us[[1]][1,1,],each=2L),2,2)),
      new_C=test$C[rep(1L,2L),,drop=FALSE],m_vec=rep(4L,4),n_draw=1L,seed=seed),controls))
    atomic_rds(model,path)
  }
  N<-nrow(test$X);S<-length(us)*cfg$latent;mu<-v<-matrix(NA_real_,N,S)
  for(s in seq_along(us))for(i in seq_len(N)) {
    p<-cc_predict(base,model,cbind(test$X[rep(i,cfg$latent),,drop=FALSE],matrix(us[[s]][,i,],cfg$latent,2)),test$C[rep(i,cfg$latent),,drop=FALSE])
    ids<-(s-1)*cfg$latent+seq_len(cfg$latent);mu[i,ids]<-p$mean;v[i,ids]<-p$var
  }
  ids<-unlist(lapply(seq_along(us),function(s)(s-1)*cfg$latent+unique(round(seq(1,cfg$latent,length.out=min(cfg$latent,4L))))))
  list(y_mean=mu[,ids,drop=FALSE],y_var=v[,ids,drop=FALSE],m_point=rowMeans(mu),optimization_status=model$optimization_status)
}
normal_abs <- function(a,s) {s<-pmax(s,1e-12);2*s*dnorm(a/s)+a*(2*pnorm(a/s)-1)}
mixture_score <- function(y,mu,v) {
  if(any(!is.finite(c(y,mu,v)))||any(v< -1e-10))stop("Invalid predictive components")
  v<-pmax(v,1e-14);sd<-sqrt(v);lp<-dnorm(y,mu,sd,log=TRUE);mx<-max(lp)
  q<-function(p)uniroot(function(x)mean(pnorm(x,mu,sd))-p,c(min(mu-12*sd),max(mu+12*sd)),tol=1e-8)$root
  lo<-q(.025);hi<-q(.975)
  # Block pairwise calculation bounds memory for CC mixtures (64*64 components).
  pair<-0;for(ids in split(seq_along(mu),ceiling(seq_along(mu)/128)))
    pair<-pair+sum(normal_abs(outer(mu[ids],mu,"-"),sqrt(outer(v[ids],v,"+"))))
  c(sq_error=(mean(mu)-y)^2,crps=mean(normal_abs(y-mu,sd))-.5*pair/length(mu)^2,
    nlpd=-(mx+log(mean(exp(lp-mx)))),coverage=as.numeric(y>=lo&&y<=hi),width=hi-lo)
}
score_result <- function(pred,dat,cal,truth,method) {
  train<-dat$train;test<-dat$test;key<-function(C)apply(C,1,paste,collapse=":")
  strata<-ifelse(key(test$C)%in%key(train$C),"seen","unseen")
  cc_support<-key(test$C)%in%key(train$C[cal,,drop=FALSE])
  rows<-list()
  add<-function(target,values,group=strata,component=rep(NA_integer_,length(group))) {
    z<-as.data.frame(values);z$target<-target;z$id<-seq_len(nrow(z));z$stratum<-group;z$coordinate<-component
    z$seen_in_complete_cases<-if(nrow(z)==length(cc_support))cc_support else NA
    rows[[length(rows)+1L]]<<-z
  }
  add("y",t(vapply(seq_along(test$y),function(i)mixture_score(test$y[i],pred$y_mean[i,],pred$y_var[i,]),numeric(5))))
  mp<-if(!is.null(pred$m_point))pred$m_point else rowMeans(pred$m_mean)
  add("m",data.frame(sq_error=(mp-truth$m)^2,crps=NA_real_,nlpd=NA_real_,coverage=NA_real_,width=NA_real_))
  if(!is.null(pred$f_mean)) {
    add("f",t(vapply(seq_along(truth$f),function(i)mixture_score(truth$f[i],pred$f_mean[i,],pred$f_var[i,]),numeric(5))))
    # m intervals use GP covariance under integration, with finite-L approximation.
    rows[[2]]$coverage<-vapply(seq_along(truth$m),function(i)mixture_score(truth$m[i],pred$m_mean[i,],pred$m_var[i,])["coverage"],numeric(1))
    rows[[2]]$width<-vapply(seq_along(truth$m),function(i)mixture_score(truth$m[i],pred$m_mean[i,],pred$m_var[i,])["width"],numeric(1))
  }
  if(!is.null(pred$U)) {
    miss<-setdiff(seq_len(nrow(train$X)),cal)
    for(d in 1:2) {
      u<-pred$U[,miss,d,drop=FALSE];dim(u)<-c(dim(pred$U)[1],length(miss))
      lo<-apply(u,2,quantile,.025);hi<-apply(u,2,quantile,.975);t<-train$U[miss,d]
      add("U",data.frame(sq_error=(colMeans(u)-t)^2,crps=NA_real_,nlpd=NA_real_,coverage=as.numeric(t>=lo&t<=hi),width=hi-lo),
        rep("missing_training",length(miss)),rep(d,length(miss)))
      rows[[length(rows)]]$id<-miss
    }
  }
  out<-do.call(rbind,rows);out$method<-method;out
}
aggregate_scores <- function(jobs,root) {
  rows<-list();status<-list()
  for(j in seq_len(nrow(jobs))) {
    p<-file.path(root,jobs$job[j],"scores.rds");err<-file.path(root,jobs$job[j],"failure.rds")
    dg<-file.path(root,jobs$job[j],"diagnostics.rds")
    passed<-if(file.exists(dg))readRDS(dg)$pass else NA
    state<-if(file.exists(p))"complete"else if(file.exists(err))"failed"else"pending"
    status[[j]]<-data.frame(jobs[j,],status=state,diagnostics_pass=passed,message=if(state=="failed")readRDS(err)$message else "")
    if(state=="complete") {z<-readRDS(p);z$rep<-jobs$rep[j];z$calib<-jobs$calib[j];z$diagnostics_pass<-passed;rows[[length(rows)+1]]<-z}
  }
  write.csv(do.call(rbind,status),file.path(root,"status.csv"),row.names=FALSE)
  if(!length(rows))return(invisible(NULL))
  points<-do.call(rbind,rows);write.csv(points,file.path(root,"point_scores.csv"),row.names=FALSE)
  both<-rbind(points,transform(points,stratum="overall"))
  groups<-split(both,interaction(both$rep,both$calib,both$method,both$target,both$stratum,drop=TRUE))
  reps<-do.call(rbind,lapply(groups,function(z){
    out<-z[1,c("rep","calib","method","target","stratum")];out$n<-nrow(z);out$diagnostics_pass<-all(z$diagnostics_pass)
    for(k in c("sq_error","crps","nlpd","coverage","width"))out[[k]]<-if(all(is.na(z[[k]])))NA_real_ else mean(z[[k]],na.rm=TRUE)
    out$rmse<-sqrt(out$sq_error);out
  }))
  write.csv(reps,file.path(root,"replicate_scores.csv"),row.names=FALSE)
  groups<-split(reps,interaction(reps$calib,reps$method,reps$target,reps$stratum,drop=TRUE))
  sums<-do.call(rbind,lapply(groups,function(z){
    out<-z[1,c("calib","method","target","stratum")];out$n_rep<-nrow(z);out$n_diagnostics_pass<-sum(z$diagnostics_pass,na.rm=TRUE);out$n_points<-sum(z$n)
    for(k in c("rmse","crps","nlpd","coverage","width")) {
      v<-z[[k]];out[[k]]<-if(all(is.na(v)))NA_real_ else mean(v,na.rm=TRUE)
      out[[paste0(k,"_se")]]<-if(sum(is.finite(v))>1)sd(v,na.rm=TRUE)/sqrt(sum(is.finite(v)))else NA_real_
    };out
  }))
  write.csv(sums,file.path(root,"summary.csv"),row.names=FALSE)
  # Pair on replication, calibration, target and stratum; never pool test rows as independent replications.
  pairs<-list();methods<-unique(reps$method)
  if(length(methods)>1)for(a in seq_len(length(methods)-1))for(b in (a+1):length(methods)) {
    z<-merge(reps[reps$method==methods[a],],reps[reps$method==methods[b],],by=c("rep","calib","target","stratum"))
    if(nrow(z)) {z$method_a<-methods[a];z$method_b<-methods[b];for(k in c("rmse","crps","nlpd","coverage","width"))z[[paste0(k,"_difference")]]<-z[[paste0(k,".x")]]-z[[paste0(k,".y")]];pairs[[length(pairs)+1]]<-z}
  }
  if(length(pairs))write.csv(do.call(rbind,pairs),file.path(root,"paired_differences.csv"),row.names=FALSE)
  invisible(sums)
}
run_supplement <- function(cfg,action="plan",script_path=NULL) {
  if(!action%in%c("plan","data","run","report"))stop("action must be plan, data, run or report")
  numerical_package_check()
  jobs<-expand.grid(rep=cfg$reps,calib=cfg$calib,method=cfg$methods,stringsAsFactors=FALSE)
  jobs$job<-sprintf("rep%03d-cal%02d-%s",jobs$rep,jobs$calib,jobs$method)
  dir.create(cfg$output,recursive=TRUE,showWarnings=FALSE)
  files<-setNames(file.path(cfg$data_root,"primary_q4_calibration",sprintf("study2_primary_q4_n100_ntest100_m4_rep%03d.rds",cfg$reps)),cfg$reps)
  if(action=="plan") {
    print(jobs);cat(nrow(jobs),"jobs; missing datasets will be generated on run. No fits started.\n")
    return(invisible(jobs))
  }
  if(action %in% c("data","run")) {
    for(r in cfg$reps)if(!file.exists(files[as.character(r)])) {
      artifact<-F("make_study2_synthetic_dataset")(rep_id=r,scenario="primary",n=100L,n_test=100L,
        q=4L,m=4L,calib_grid=c(0L,20L,50L,80L),data_seed_base=1000000L)
      F("store_mixedgp_synthetic_dataset")(artifact,dirname(files[as.character(r)]),overwrite=FALSE)
    }
  }
  if(any(!file.exists(files)))stop("Missing datasets; run with --action=data first")
  for(r in cfg$reps) {
    artifact<-readRDS(files[as.character(r)])
    expected<-list(n=100L,n_test=100L,q=4L,m=4L,d=2L,calib_grid=c(0L,20L,50L,80L))
    if(!identical(artifact$study,"study2")||!identical(artifact$scenario,"primary")||
       !identical(artifact$rep_id,as.integer(r))||
       !isTRUE(all.equal(artifact$design[names(expected)],expected,check.attributes=FALSE)))stop("Frozen design mismatch: ",files[as.character(r)])
  }
  if(action=="data")return(invisible(files))
  pkgs<-c("eivGP","TruncatedNormal","posterior","kergp","LVGP","EzGP")
  versions<-setNames(vapply(pkgs,function(p)if(requireNamespace(p,quietly=TRUE))as.character(packageVersion(p))else"missing",character(1)),pkgs)
  # A root is one scientific configuration; job selection/cores can change on resume.
  scientific<-cfg[setdiff(names(cfg),c("reps","methods","calib","output","data_root","cores"))]
  identity<-list(settings=scientific,versions=versions,R=as.character(getRversion()),
    code=c(supplement=if(is.null(script_path))hash_object(body(run_supplement)) else unname(tools::md5sum(script_path)),
      common=unname(tools::md5sum(file.path(NUMERICAL_REPO_ROOT,"paper/numerical/common.R")))) ,
    package_manifest=F("eivgp_module_identity")(),
    dependency_lock=unname(tools::md5sum(file.path(NUMERICAL_REPO_ROOT,"renv.lock"))),
    full_map_hash=if(nzchar(cfg$full_fit_map))unname(tools::md5sum(cfg$full_fit_map))else NULL)
  manifest<-file.path(cfg$output,"protocol.rds")
  if(file.exists(manifest)&&!identical(readRDS(manifest),identity))stop("Configuration/code changed: use a new output directory")
  if(!file.exists(manifest))atomic_rds(identity,manifest)
  capture.output(sessionInfo(),file=file.path(cfg$output,"sessionInfo.txt"))
  revision<-tryCatch(system2("git",c("-C",shQuote(NUMERICAL_REPO_ROOT),"rev-parse","HEAD"),stdout=TRUE,stderr=FALSE),error=function(e)"unavailable")
  writeLines(revision,file.path(cfg$output,"git_revision.txt"))
  write.csv(jobs,file.path(cfg$output,"requested_jobs.csv"),row.names=FALSE)
  for(r in cfg$reps) {
    fp<-file.path(cfg$output,sprintf("data-rep%03d.rds",r));record<-list(file=basename(files[as.character(r)]),md5=unname(tools::md5sum(files[as.character(r)])))
    if(file.exists(fp)&&!identical(readRDS(fp),record))stop("Frozen dataset changed for replication ",r)
    if(!file.exists(fp))atomic_rds(record,fp)
  }
  if(action=="plan") {print(jobs);cat(nrow(jobs),"jobs; completed scores will be skipped. Profile:",cfg$profile,"\n");return(invisible(jobs))}
  registry<-file.path(cfg$output,"job_registry.rds")
  alljobs<-if(file.exists(registry))unique(rbind(readRDS(registry),jobs))else jobs
  if(action=="report")return(aggregate_scores(alljobs,cfg$output))
  runlock<-file.path(cfg$output,"RUN.lock")
  if(!dir.create(runlock,showWarnings=FALSE))stop("Run already locked: ",runlock)
  on.exit(unlink(runlock,recursive=TRUE),add=TRUE)
  writeLines(as.character(Sys.getpid()),file.path(runlock,"pid"))
  atomic_rds(alljobs,registry)
  failed<-0L
  for(j in seq_len(nrow(jobs))) {
    job<-jobs[j,];folder<-file.path(cfg$output,job$job);scorefile<-file.path(folder,"scores.rds")
    if(file.exists(scorefile))next
    dir.create(folder,recursive=TRUE,showWarnings=FALSE)
    # Atomic lock prevents two workers writing the same job. Remove only after checking no process owns it.
    lock<-paste0(folder,".lock");if(!dir.create(lock,showWarnings=FALSE))stop("Job locked: ",lock)
    writeLines(as.character(Sys.getpid()),file.path(lock,"pid"))
    message("Running ",job$job)
    tryCatch({
      artifact<-readRDS(files[as.character(job$rep)]);dat<-artifact$data
      stopifnot(artifact$study=="study2",artifact$scenario=="primary",artifact$design$q==4,artifact$design$d==2)
      cal<-artifact$calibration_sets[[as.character(job$calib)]];stopifnot(length(cal)==job$calib)
      if(qr(cbind(1,dat$train$U[cal,,drop=FALSE]))$rank<3L)stop("Calibration fails affine-rank anchoring")
      Uobs<-matrix(NA_real_,nrow(dat$train$U),2);Uobs[cal,]<-dat$train$U[cal,]
      base_seed<-cfg$seed+job$rep*1000L+job$calib
      shared<-file.path(cfg$output,sprintf("shared-rep%03d-cal%02d",job$rep,job$calib));dir.create(shared,showWarnings=FALSE)
      mp<-file.path(shared,"measurement.rds")
      if(job$method=="CUT"||grepl("-CC$",job$method)) {
        if(!file.exists(mp))atomic_rds(measurement_cut(dat$train$C,Uobs,cfg,base_seed),mp)
        meas<-readRDS(mp)
        write.csv(meas$diagnostics,file.path(shared,"measurement_diagnostics.csv"),row.names=FALSE)
      }
      tp<-file.path(shared,"truth.rds")
      if(!file.exists(tp))atomic_rds(list(m=F("oracle_m0_quadrature_2d")(dat$test$X,dat$test$C,dat$true_params),
        f=F("f0_2d")(dat$test$X,dat$test$U)),tp)
      truth<-readRDS(tp)
      if(job$method=="CUT") {
        passes<-diagnostic_pass(meas$diagnostics);states<-list()
        for(b in seq_along(meas$states)) {
          cp<-file.path(folder,sprintf("conditional-%03d.rds",b))
          sp<-sub("[.]rds$","-states.rds",cp)
          if(!file.exists(sp)) {
            if(!file.exists(cp))atomic_rds(fit_full(dat$train,meas$states[[b]]$U,cfg,base_seed+200000L+b,TRUE),cp)
            fit<-F("mixedgp_read_artifact")(cp)
            passed<-record_gp_diagnostics(fit,sub("rds$","diagnostics.csv",cp),TRUE)
            atomic_rds(list(pass=passed,states=select_gp(fit,cfg$inner_draws,meas$states[[b]])),sp)
          }
          # Publish and verify prediction states before retiring a passing fit.
          saved<-retire_conditional_fit(cp,sp)
          passes<-c(passes,saved$pass)
          states<-c(states,saved$states)
        }
        pred<-gp_evaluate(states,dat$train,dat$test,cfg,base_seed+300000L)
        # Keep every measurement draw for imputation, not just integration states.
        mis<-which(!complete.cases(Uobs));arr<-meas$trace[,,seq_len(length(mis)*2),drop=FALSE]
        pred$U<-array(NA_real_,c(dim(arr)[1]*4,nrow(Uobs),2));pred$U[,mis,]<-matrix(arr,ncol=length(mis)*2)
      } else if(job$method=="FULL") {
        fp<-file.path(folder,"fit.rds")
        if(nzchar(cfg$full_fit_map)) {
          map<-read.csv(cfg$full_fit_map,stringsAsFactors=FALSE)
          if(!all(c("rep","calib","path")%in%names(map)))stop("FULL map requires rep,calib,path columns")
          selected<-map[map$rep==job$rep & map$calib==job$calib,,drop=FALSE]
          if(nrow(selected)!=1||!file.exists(selected$path))stop("FULL map must supply exactly one existing fit for this job")
          provenance<-list(path=normalizePath(selected$path),md5=unname(tools::md5sum(selected$path)))
          sourcefile<-file.path(folder,"fit_source.rds")
          if(file.exists(sourcefile)&&!identical(readRDS(sourcefile),provenance))stop("Saved FULL source changed")
          fit<-read_matched_full(selected$path,dat$train,Uobs);atomic_rds(provenance,sourcefile)
        } else {
          if(!file.exists(fp))atomic_rds(fit_full(dat$train,Uobs,cfg,base_seed+400000L),fp)
          fit<-F("mixedgp_read_artifact")(fp)
        }
        passes<-record_gp_diagnostics(fit,file.path(folder,"full_diagnostics.csv"));states<-select_gp(fit,cfg$outer*cfg$inner_draws)
        pred<-gp_evaluate(states,dat$train,dat$test,cfg,base_seed+300000L);pred$U<-fit$mcmc$samples_U
      } else if(job$method%in%c("UC-GP","LVGP","EzGP")) {
        cp<-file.path(cfg$output,sprintf("categorical-rep%03d",job$rep),job$method,"predictions.rds")
        pred<-categorical_evaluate(job$method,dat$train,dat$test,cfg,cfg$seed+job$rep*1000L+600000L,cp)
        passes<-is.null(pred$optimization_status)||!grepl("nonconverged",pred$optimization_status)
      } else {
        up<-file.path(shared,"prospective.rds")
        if(!file.exists(up))atomic_rds(prospective_u(meas$states,dat$test$C,cfg$latent,base_seed+300000L),up)
        passes<-diagnostic_pass(meas$diagnostics)
        pred<-cc_evaluate(job$method,dat$train,dat$test,cal,readRDS(up),cfg,base_seed+500000L,file.path(folder,"fit.rds"))
        if(!is.null(pred$optimization_status)&&grepl("nonconverged",pred$optimization_status))passes<-c(passes,FALSE)
      }
      atomic_rds(list(pass=all(passes),checks=passes,profile=cfg$profile,
        note="MCMC screening only; also inspect target stability and nested integration error before publication"),file.path(folder,"diagnostics.rds"))
      if(!is.null(pred$m_mean)) {
        cols<-seq_len(ncol(pred$m_mean));half<-length(cols)%/%2
        stability<-data.frame(test_id=seq_len(nrow(pred$m_mean)),
          m_first=rowMeans(pred$m_mean[,cols<=half,drop=FALSE]),m_second=rowMeans(pred$m_mean[,cols>half,drop=FALSE]))
        stability$difference<-stability$m_second-stability$m_first
        write.csv(stability,file.path(folder,"integration_stability.csv"),row.names=FALSE)
      }
      scores<-score_result(pred,dat,cal,truth,job$method)
      # Imputation summaries are in scores; full draws remain in the source fit
      # or shared measurement trace, so do not duplicate them per prediction file.
      pred$U<-NULL
      atomic_rds(pred,file.path(folder,"predictions.rds"))
      atomic_rds(scores,scorefile);unlink(file.path(folder,"failure.rds"))
    },error=function(e){failed<<-failed+1L;atomic_rds(list(message=conditionMessage(e),time=Sys.time()),file.path(folder,"failure.rds"));message("FAILED: ",conditionMessage(e))},finally=unlink(lock,recursive=TRUE))
  }
  aggregate_scores(alljobs,cfg$output)
  if(failed)stop(failed," jobs failed; see status.csv. Rerun the identical command to retry.")
  invisible(jobs)
}
if(sys.nframe()==0L) {
  a<-parse_cli(commandArgs(TRUE));cfg<-settings(a)
  script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
  run_supplement(cfg,if(is.null(a$action))"plan"else a$action,normalizePath(script))
}
