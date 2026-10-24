test_that("internal Study I continuation matches uninterrupted chains exactly", {
  set.seed(44);x<-rnorm(12);u<-rnorm(12);y<-x+u+rnorm(12)
  c<-as.integer(cut(u,c(-Inf,-.5,.5,Inf)))
  args<-list(x_raw=x,y_raw=y,c_ord=c,burn=10L,n_chains=4L,seed=55L,parallel_chains=FALSE,m=3L)
  a<-do.call(eivGP:::fit_eivgp_1d_v030,c(args,list(n_iter=40L)))
  p<-tempfile(fileext=".rds");on.exit(unlink(p))
  eivGP:::mixedgp_save_artifact(a,p);a<-read_fit(p)
  b<-eivGP:::mixedgp_continue_study1_fit(a,20L)
  z<-do.call(eivGP:::fit_eivgp_1d_v030,c(args,list(n_iter=60L)))
  expect_identical(b$mcmc$samples_by_chain,z$mcmc$samples_by_chain)
  expect_identical(b$checkpoint$states,z$checkpoint$states)
  serial<-eivGP:::mixedgp_continuation_target_series(a,matrix(x[1:2],ncol=1),c[1:2],NULL,8L,123L,FALSE,1L)
  parallel<-eivGP:::mixedgp_continuation_target_series(a,matrix(x[1:2],ncol=1),c[1:2],NULL,8L,123L,TRUE,4L)
  expect_identical(serial,parallel)
})

test_that("scientific ESS is strict and requires a complete four-chain window", {
  tab<-data.frame(ess_bulk=101,ess_tail=101,n_chains=4L,draw_window_complete=TRUE)
  expect_true(eivGP:::mixedgp_scientific_ess_pass(tab))
  for(value in c(100,99,NA,Inf)) {
    z<-tab;z$ess_tail<-value;expect_false(eivGP:::mixedgp_scientific_ess_pass(z))
  }
  z<-tab;z$n_chains<-3L;expect_false(eivGP:::mixedgp_scientific_ess_pass(z))
  for(value in c(FALSE,NA)) {
    z<-tab;z$draw_window_complete<-value;expect_false(eivGP:::mixedgp_scientific_ess_pass(z))
  }
})

test_that("an operational limit preserves unsuccessful states without claiming success", {
  tab<-data.frame(ess_bulk=80,ess_tail=90,n_chains=4L,draw_window_complete=TRUE,rhat=1.1)
  local_mocked_bindings(mixedgp_study1_target_series=function(...)list(a=rep(list(1:20),4)),
    mixedgp_summarize_diagnostic_series=function(...)tab,
    mixedgp_fit_engine=function(...)"univariate",.package="eivGP")
  object<-list(mcmc=list(samples_by_chain=list(logtheta=rep(list(matrix(0,20,1)),4))),checkpoint=list(iteration=30L))
  p<-tempfile(fileext=".rds");on.exit(unlink(p))
  r<-continue_eivgp_until_ess(object,matrix(0,1,1),1,max_extensions=0,checkpoint_file=p)
  expect_false(r$success);expect_identical(read_fit(p)$fit,object)
})

test_that("cleanup retains checkpoint files when any target is deficient", {
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  dir.create(file.path(root,"config"));dir.create(file.path(root,"cells/c"),recursive=TRUE)
  dir.create(file.path(root,"checkpoints"));p<-file.path(root,"checkpoints","chain.rds");saveRDS(list(state=1),p)
  cfg<-list(study="study1",cells=list(list(id="c",n_rep=1L,calibration_grid=0L)))
  saveRDS(cfg,file.path(root,"config","resolved_config.rds"))
  tab<-expand.grid(parameter=paste0(rep(c("m_conditional_mean","m_conditional_variance","predictive_second_moment"),each=5),"[",rep(1:5,3),"]"))
  tab$rep<-1L;tab$n_calib<-0L;tab$diagnostic_role<-"target";tab$ess_bulk<-101
  tab$ess_tail<-101;tab$ess_tail[1]<-100;tab$n_chains<-4L;tab$draw_window_complete<-TRUE
  ctx<-list(mcmc_parameter_diagnostics=tab,STUDY1_HETEROGENEITY_ETA=0)
  saveRDS(list(context=ctx),file.path(root,"cells/c/report_inputs.rds"))
  expect_error(clean_study1_checkpoints(root),"successful scientific-target")
  expect_true(file.exists(p))
})

test_that("continued targets stay visible while the report keeps three trace figures", {
  root<-tempfile();dir.create(root);dir.create(file.path(root,"combined"));on.exit(unlink(root,recursive=TRUE))
  z<-expand.grid(cell_id="eta1_balanced",n_calib=c(0L,20L),rep=1L,
    family=c("m_conditional_mean","predictive_second_moment"),stringsAsFactors=FALSE)
  z$parameter<-paste0(z$family,"[1]");z$ess_bulk<-200;z$ess_tail<-150;z$rhat<-1.02
  utils::write.csv(data.frame(cell_id="eta1_balanced",rep=1L,n_calib=c(0L,20L),n_iter=c(60000L,40000L)),
    file.path(root,"combined/mcmc_diagnostics.csv"),row.names=FALSE)
  selected<-eivGP:::mixedgp_study1_trace_selection(z,root)
  expect_equal(sum(selected$selection=="Continued scientific targets"),4)
  expect_setequal(unique(selected$selection),c("Typical mean","Weakest uncertainty","Continued scientific targets"))
})

test_that("cleanup removes only verified redundant checkpoints and keeps the final fit", {
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  dir.create(file.path(root,"config"));dir.create(file.path(root,"cells/c/results/study1_publication"),recursive=TRUE)
  dir.create(file.path(root,"checkpoints/id"),recursive=TRUE)
  cfg<-list(study="study1",cells=list(list(id="c",n_rep=1L,calibration_grid=0L)))
  saveRDS(cfg,file.path(root,"config/resolved_config.rds"))
  tab<-data.frame(parameter=paste0(rep(c("m_conditional_mean","m_conditional_variance","predictive_second_moment"),each=5),"[",rep(1:5,3),"]"),
    rep=1L,n_calib=0L,diagnostic_role="target",ess_bulk=101,ess_tail=101,n_chains=4L,draw_window_complete=TRUE)
  ctx<-list(mcmc_parameter_diagnostics=tab,STUDY1_HETEROGENEITY_ETA=0,
    mc_results=data.frame(rep=1L,n_calib=0L,method="EIV-GP",RMSE=1))
  saveRDS(list(context=ctx),file.path(root,"cells/c/report_inputs.rds"))
  d<-list(x=matrix(1:3,ncol=1),y=c(1,2,3),c_ord=1:3,u_obs=rep(NA_real_,3))
  fit<-list(data=d,checkpoint=list(version=2L,iteration=30L,burn=10L,
    states=rep(list(list(state=list(u=0),rng=list(value=1))),4)),
    mcmc=list(samples_by_chain=list(u=rep(list(matrix(0,20,3)),4),logtheta=rep(list(matrix(0,20,3)),4))))
  p<-file.path(root,"cells/c/results/study1_publication/fit_rep001_cal000.rds")
  eivGP:::mixedgp_save_artifact(list(fit=fit),p);hash<-unname(tools::md5sum(p))
  model<-list(seed=301000L,n_iter=30L,X=d$x,y=d$y,C=matrix(d$c_ord,ncol=1),U_obs=matrix(d$u_obs,ncol=1))
  saveRDS(list(model=model),file.path(root,"checkpoints/id/identity.rds"))
  saveRDS(list(state=0),file.path(root,"checkpoints/id/chain-01.rds"))
  audit<-clean_study1_checkpoints(root)
  expect_equal(nrow(audit),2);expect_true(all(audit$removed))
  expect_length(list.files(file.path(root,"checkpoints"),recursive=TRUE),0)
  expect_identical(unname(tools::md5sum(p)),hash)
  expect_equal(nrow(clean_study1_checkpoints(root)),0L)
  unlink(file.path(root,"checkpoints"),recursive=TRUE)
  expect_equal(nrow(clean_study1_checkpoints(root)),0L)
})
