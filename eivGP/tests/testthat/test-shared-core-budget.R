with_budget <- function(cores,expr) {
  close <- mixedgp_budget_enter(cores)
  on.exit(close(),add=TRUE)
  force(expr)
}

test_that("nested stages lend slots, respect the total budget and preserve seeds", {
  events<-tempfile();dir.create(events)
  on.exit(unlink(events,recursive=TRUE),add=TRUE)
  old<-options(eivGP.core_events=events);on.exit(options(old),add=TRUE)
  run<-function(cores)with_budget(cores,
    mixedgp_parallel_lapply(1:2,function(i) {
      chains<-mixedgp_parallel_lapply(1:4,function(j){Sys.sleep(.015);runif(3)},
        n_cores=4L,seeds=100L*i+1:4)
      prediction<-mixedgp_stage_lapply(1:8,function(j){Sys.sleep(.005);sum(chains[[1]])+j})
      list(chains=chains,prediction=prediction)
    },n_cores=2L,seeds=11:12))
  serial<-run(1L)
  unlink(list.files(events,full.names=TRUE))
  expect_identical(run(4L),serial)
  records<-do.call(rbind,lapply(list.files(events,full.names=TRUE),read.csv,header=FALSE))
  records<-records[order(records$V1),]
  active<-cumsum(ifelse(records$V3=="acquire",1L,-1L))
  expect_true(all(active>=0L & active<=4L))
  expect_equal(tail(active,1),0)
  if(.Platform$OS.type=="unix")expect_equal(max(active),4L)
  expect_null(mixedgp_budget())
  expect_false(mixedgp_budget_owned())
})

test_that("worker failures release slots and permit subsequent work", {
  with_budget(2L,{
    expect_error(suppressWarnings(mixedgp_stage_lapply(1:2,function(i)stop("expected failure"))),"failure|worker failed")
    expect_identical(unlist(mixedgp_stage_lapply(1:4,function(i)i*i)),c(1L,4L,9L,16L))
    expect_true(mixedgp_budget_owned())
  })
  expect_null(mixedgp_budget())
})

test_that("mapply shares the same budget and preserves names and RNG", {
  run<-function(b)with_budget(b,mixedgp_parallel_mapply(function(a,b)a+b+runif(1),
    a=c(first=1,second=2),b=3,n_cores=4,seeds=31:32,SIMPLIFY=TRUE))
  expect_identical(run(1),run(4))
  expect_identical(names(run(2)),c("first","second"))
})

test_that("four-chain fits and prediction targets reproduce across core counts", {
  args<-list(X=matrix(seq(-1,1,length.out=8),ncol=1),y=sin(seq(-1,1,length.out=8)),
    C=matrix(rep(1:2,each=4),ncol=1),engine="univariate",m_vec=2L,
    n_iter=20L,burn=8L,n_chains=4L,seed=736L,verbose=FALSE)
  old<-options(eivGP.checkpoint_root=NULL);on.exit(options(old),add=TRUE)
  fit<-with_budget(1L,do.call(fit_eivgp,c(args,list(n_cores=1L))))
  parallel_fit<-with_budget(4L,do.call(fit_eivgp,c(args,list(n_cores=4L))))
  expect_identical(fit$mcmc$samples_by_chain,parallel_fit$mcmc$samples_by_chain)
  expect_equal(length(fit$mcmc$samples_by_chain$u),4L)
  for(target in c("response","mean","surface"))for(joint in c(FALSE,TRUE)) {
    a<-list(object=fit,new_X=matrix(c(-.25,.3),ncol=1),target=target,
      draw_ids=c(1:9,17:24),n_per_draw=2L,joint=joint,seed=123L)
    if(target=="surface") {a$new_U<-matrix(c(-.2,.2),ncol=1);a$new_U_scale<-"model"}
    else a$new_C<-matrix(c(1L,2L),ncol=1)
    if(target=="mean"){a$n_latent<-8L;a$return_components<-TRUE}
    set.seed(88);rng<-.Random.seed
    reference<-with_budget(1L,do.call(predict_eivgp,a))
    expect_identical(.Random.seed,rng)
    for(b in c(2L,4L))expect_identical(with_budget(b,do.call(predict_eivgp,a)),reference)
    expect_equal(dim(reference),c(34L,2L))
    if(target!="surface")expect_equal(ncol(attr(reference,"conditional_means")),2L)
  }
})

test_that("plot workers produce the same PNG contents", {
  skip_if_not_installed("ggplot2")
  root<-tempfile();on.exit(unlink(root,recursive=TRUE),add=TRUE)
  p<-ggplot2::ggplot(data.frame(x=1:3,y=3:1),ggplot2::aes(x,y))+ggplot2::geom_line()
  for(b in c(1L,4L))with_budget(b,mixedgp_save_plot_jobs(list(a=p,b=p),file.path(root,b),function(n)c(3,2)))
  expect_identical(unname(tools::md5sum(file.path(root,1,c("a.png","b.png")))),
                   unname(tools::md5sum(file.path(root,4,c("a.png","b.png")))))
})

test_that("multivariate prediction metadata and joint draws survive parallel blocks", {
  n<-8L;u<-seq(-1,1,length.out=n);U<-matrix(u,ncol=1);U[3:6,]<-NA_real_
  X<-matrix(u,ncol=1);C<-cbind(as.integer(u>0)+1L,rep(1:2,4))
  fit<-fit_eivgp(X,sin(u),C,U_obs=U,engine="multivariate",ident="none",
    n_iter=16L,burn=4L,n_chains=4L,parallel=FALSE,seed=91L,verbose=FALSE)
  for(target in c("response","mean","surface")) {
    a<-list(object=fit,new_X=X[1:2,,drop=FALSE],draw_ids=1:17,
      n_per_draw=2L,seed=128L,target=target,joint=TRUE)
    if(target=="surface")a$new_U<-U[1:2,,drop=FALSE]else a$new_C<-C[1:2,,drop=FALSE]
    if(target=="mean"){a$n_latent<-4L;a$return_components<-TRUE}
    serial<-with_budget(1L,do.call(predict_eivgp,a))
    expect_identical(with_budget(4L,do.call(predict_eivgp,a)),serial)
    expect_equal(dim(serial),c(34L,2L))
    if(target=="response")expect_equal(length(attr(serial,"rejection_acceptance")),34L)
  }
})

test_that("numerical posterior integration and scoring reproduce in the pool", {
  api<-eivgp_numerical_api(tempdir())
  train<-list(X=matrix(seq(-1,1,length.out=6),ncol=1),U=cbind(seq(-.5,.5,length.out=6),seq(.3,-.3,length.out=6)),
    y=seq(-.2,.8,length.out=6),C=matrix(rep(1:4,6),6,4))
  test<-list(X=matrix(c(-.2,.4),ncol=1),U=matrix(c(-.1,.2,.2,-.2),2,2),y=c(.2,.7),C=matrix(c(1,2,3,4,2,3,4,1),2,4))
  state<-list(U=train$U,A=matrix(c(.5,.3,.2,.1,0,.6,.2,.4),4,2),tau=rep(list(c(-1,0,1)),4),
    logtheta=rep(0,4),sigma2=.2)
  states<-rep(list(state),9L);cfg<-list(latent=4L,gp_response_latent=2L)
  serial<-with_budget(1L,api$gp_evaluate(states,train,test,cfg,765L))
  parallel<-with_budget(4L,api$gp_evaluate(states,train,test,cfg,765L))
  expect_identical(serial,parallel)
  dat<-list(train=train,test=test);truth<-list(m=c(.1,.5),f=c(.2,.6))
  expect_identical(with_budget(1L,api$score_result(serial,dat,1:3,truth,"FULL")),
                   with_budget(4L,api$score_result(parallel,dat,1:3,truth,"FULL")))
  path<-tempfile();saveRDS(list(optimization_status="ok"),path);on.exit(unlink(path),add=TRUE)
  api$cc_predict<-function(method,model,X,C,m_vec)list(mean=rowSums(X),var=rep(.2,nrow(X)))
  us<-api$prospective_u(states,test$C,4L,456L)
  expect_identical(with_budget(1L,api$cc_evaluate("UC-GP-CC",train,test,1:3,us,cfg,45L,path)),
                   with_budget(4L,api$cc_evaluate("UC-GP-CC",train,test,1:3,us,cfg,45L,path)))
})

test_that("an abruptly terminated worker cannot strand the shared budget", {
  skip_if(.Platform$OS.type!="unix")
  with_budget(2L,{
    expect_error(suppressWarnings(mixedgp_stage_lapply(1:2,function(i) {
      if(i==1L)tools::pskill(Sys.getpid(),9L)
      Sys.sleep(.02);i
    })),"worker failed")
    expect_identical(unlist(mixedgp_stage_lapply(1:3,identity)),1:3)
  })
})
