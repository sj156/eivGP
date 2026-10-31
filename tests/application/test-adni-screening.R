# Author: Xinye Chen. Synthetic contract checks; no participant records in tests.
repo<-Sys.getenv('EIVGP_TEST_REPO');stopifnot(nzchar(repo))
source(file.path(repo,'tools/adni.R'))
list2env(as.list(adni_load_api(repo)),envir=environment())
a<-adni_screen_api(repo)
fixture<-function(root,q=4L,plan='Plan1-6-1'){
 set.seed(819L);n<-40L
 d<-data.frame(RID=seq_len(n),Y=rnorm(n),X1=60+runif(n)*20,X2=12+runif(n)*8,U1=rnorm(n),U2=rnorm(n),calibration=0L)
 C<-list(rep(1:4,10),rep(c(1,3,2,4),10),sample(rep(1:4,10)),sample(rep(1:4,10)))
 for(j in seq_len(q))d[[paste0('C',j)]]<-C[[j]]
 frozen<-list()
 for(f in 1:5)for(cell in c('O0','O20','O50')){
  te<-d[seq(f,n,by=5),];tr<-d[-seq(f,n,by=5),];nc<-switch(cell,O0=0L,O20=6L,O50=16L);if(nc)tr$calibration[seq_len(nc)]<-1L
  path<-file.path(root,'inputs','MC',sprintf('fold%02d',f),cell);dir.create(file.path(path,'private'),recursive=TRUE)
  write.csv(tr,file.path(path,'private/train.csv'),row.names=FALSE);write.csv(te,file.path(path,'private/test.csv'),row.names=FALSE)
  p<-list(author='Xinye Chen',plan_id=plan,q=q,d=2L,p=2L,fold=f,O=cell,N=n,training_N=nrow(tr),test_N=nrow(te),observed_training_U=nc,
   U_order=c('U1','U2'),C_order=c('MEMORY3','MEMORY8','DIVATT3','MEMORY2')[seq_len(q)],X_order=c('age','education'),C_levels=1:4,m_vec=rep(4L,q),
   ident=if(cell=='O0')'lower_triangular'else'none',standardize_U=cell!='O0',test_U_available_RID=head(te$RID,switch(cell,O0=0L,O20=2L,O50=4L)),
   seeds=as.list(setNames(101:109,c('production','posterior','prediction','mixed_prediction','UC_GP','LVGP','EzGP','EzGP_draws','prospective_U'))))
  a$atomic_json(p,file.path(path,'PROTOCOL.json'))
  for(p in list.files(path,recursive=TRUE,full.names=TRUE))frozen[[normalizePath(p)]]<-digest::digest(file=p,algo='sha256')
 }
 m<-list(schema='adni-screening-v2',author='Xinye Chen',trial_id='synthetic-v2-check',smoke=TRUE,output_root=file.path(root,'workflow'),
  stages=c(24L,32L),burn=8L,folds=1L,cells='O20',cases=list(list(id=plan,input=file.path(root,'inputs'),output=file.path(root,'results'))),frozen_inputs=frozen)
 a$atomic_json(m,file.path(root,'SCREENING.json'));file.path(root,'SCREENING.json')
}
testthat::test_that('screening inputs mask hidden truth and preserve canonical identities',{
 root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE));path<-fixture(root)
 b<-adni_screen_manifest(path,a);z<-b$inputs[[1]]$cells[['1/O20']];s<-adni_screen_settings(b$manifest,1000L,2L)
 x<-a$adni_plan_fit_args(z,s);testthat::expect_identical(x$sampler,'HMC');testthat::expect_identical(x$n_chains,4L)
 testthat::expect_equal(ncol(x$C),4L);testthat::expect_equal(ncol(x$X),2L);testthat::expect_true(all(is.na(x$U_obs[z$train$calibration==0,])))
 changed<-z;changed$train[changed$train$calibration==0,z$u]<-1e8;changed$test$Y<-1e8;changed$test[,z$u]<--1e8
 testthat::expect_identical(x,a$adni_plan_fit_args(changed,s))
 testthat::expect_identical(a$adni_plan_fit_args(b$inputs[[1]]$cells[['1/O50']],s)$sampler,'ESS')
 m<-b$manifest;m$smoke<-FALSE;a$atomic_json(m,file.path(root,'bad.json'))
 testthat::expect_error(adni_screen_manifest(file.path(root,'bad.json'),a),'Production screening')
 tr<-file.path(root,'inputs/MC/fold01/O20/private/train.csv');cat('\n',file=tr,append=TRUE)
 testthat::expect_error(adni_screen_manifest(path,a),'Frozen input changed')
})
testthat::test_that('groups separate EIV, comparators and strict saved-only reporting',{
 root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE));path<-fixture(root)
 m<-jsonlite::fromJSON(path,simplifyVector=FALSE)
 e<-adni_screen_tasks(m,'eivGP','S24');c<-adni_screen_tasks(m,'competitors','BASELINE')
 testthat::expect_equal(nrow(e),1L);testthat::expect_true(all(e$kind=='eiv'))
 testthat::expect_equal(nrow(c),4L);testthat::expect_equal(sum(c$kind=='baseline_group'),1L);testthat::expect_equal(sum(c$kind=='cc'),3L)
 for(group in c('eivGP','competitors','report'))testthat::expect_error(adni_screen_run(repo,list(manifest=path,cores=2L,plan=FALSE,reuse_only=TRUE),group),'Saved-only prerequisites')
 testthat::expect_false(dir.exists(file.path(root,'results/fits')))
 testthat::expect_false(dir.exists(file.path(root,'workflow/.adni-application.lock')))
 testthat::expect_error(adni_screen_args(c(paste0('--manifest=',path),'--cores=2.5')),'positive integer')
 testthat::expect_error(adni_screen_args('--chains=2'),'Unknown')
})
testthat::test_that('public scheduler lends one total pool rather than full-chain reservations',{
 testthat::expect_false(exists('mixedgp_task_budget',envir=asNamespace('eivGP'),inherits=FALSE))
 testthat::expect_equal(a$queue$mixedgp_stage_workers(5L,'fit',10L,4L),1L)
 testthat::expect_equal(a$queue$mixedgp_stage_workers(6L,'fit',10L,4L),1L)
 testthat::expect_true(is.infinite(a$queue$mixedgp_competitor_controls()$LVGP$max_elapsed_seconds));testthat::expect_true(a$queue$mixedgp_competitor_controls()$LVGP$parallel)
 root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE));close<-a$mixedgp_budget_enter(2L);on.exit(close(),add=TRUE)
 tasks<-data.frame(plan='Plan1-6-1',fold=1:3,cell='O20',method='EIV-GP',kind='eiv',stage='S1',id=paste0('job',1:3))
 status<-adni_screen_phase(tasks,'fit',a$mixedgp_budget_cores(),a,function(t,n){stopifnot(n==a$mixedgp_budget_cores());a$mixedgp_budget_yield(function()getFromNamespace('mixedgp_stage_lapply','eivGP')(1:3,function(j)j*j));NULL},root)
 testthat::expect_true(all(status$state=='DONE'))
})
testthat::test_that('real four-chain q4 smoke continues same states and exports valid scores',{
 root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE));path<-fixture(root)
 adni_screen_run(repo,list(manifest=path,cores=4L,plan=FALSE,reuse_only=FALSE),'eivGP')
 t<-data.frame(fold=1L,cell='O20');fit<-a$adni_plan_read_fit(a$adni_plan_fit_path(file.path(root,'results'),t))
 testthat::expect_identical(fit$checkpoint$iteration,32L);testthat::expect_identical(fit$checkpoint$burn,8L)
 testthat::expect_equal(length(fit$checkpoint$states),4L);testthat::expect_equal(fit$data$q,4L);testthat::expect_identical(fit$control$sampler,'HMC')
 for(stage in c('S24','S32')){
  p<-file.path(root,'results',stage,'analysis/fold01/O20');m<-read.csv(file.path(p,'prediction_metrics.csv'))
  testthat::expect_setequal(m$protocol,c('OI','MI','CC-known'));testthat::expect_true(all(is.finite(m$RMSE_SD)))
  testthat::expect_true(file.info(file.path(p,'training_U_imputation_std.png'))$size>1000)
  testthat::expect_true(file.exists(file.path(p,'PROSPECTIVE_U_metrics.csv')))
 }
 testthat::expect_true(file.exists(file.path(root,'workflow/reports/eivGP/progress/stage_stability.csv')))
 testthat::expect_true(file.exists(file.path(root,'workflow/reports/eivGP/progress/index.html')))
 before<-tools::md5sum(a$adni_plan_fit_path(file.path(root,'results'),t))
 adni_screen_run(repo,list(manifest=path,cores=1L,plan=FALSE,reuse_only=TRUE),'eivGP')
 testthat::expect_identical(tools::md5sum(a$adni_plan_fit_path(file.path(root,'results'),t)),before)
})
testthat::test_that('CC competitors train only on calibration and resume scoring without refitting',{
 root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE));path<-fixture(root)
 e<-adni_screen_api(repo);m<-jsonlite::fromJSON(path,simplifyVector=FALSE);inp<-e$adni_plan_inputs(file.path(root,'inputs'))
 s<-adni_screen_settings(m,1000L,2L);s$profile<-'paper';captured<-NULL
 e$fit_mixedgp_competitor<-function(method,X,y,C,new_X,new_C,m_vec,n_draw,seed,...){
  captured<<-list(method=method,X=X,y=y,C=C,new_X=new_X,m_vec=m_vec,controls=list(...))
  list(predictive_mean=rep(0,nrow(new_X)),draws=matrix(0,n_draw,nrow(new_X)))
 }
 t<-data.frame(fold=1L,cell='O20',method='LVGP');out<-file.path(root,'comparison')
 e$baseline_cell(t,inp,out,s,TRUE)
 z<-inp$cells[['1/O20']];dir<-file.path(out,'CC/O20/baselines/fold01/LVGP');scores<-read.csv(file.path(dir,'metrics.csv'))
 testthat::expect_equal(nrow(captured$X),6L);testthat::expect_equal(ncol(captured$X),4L);testthat::expect_equal(ncol(captured$C),4L)
 testthat::expect_false(captured$controls$parallel);testthat::expect_true(is.infinite(captured$controls$max_elapsed_seconds))
 expected<-sqrt(mean((mean(z$train$Y[z$train$calibration==1])-z$test$Y[z$available])^2))/sd(z$train$Y)
 testthat::expect_equal(scores$RMSE_SD,expected);testthat::expect_equal(scores$n,sum(z$available));testthat::expect_equal(scores$training_N,6L)
 model<-tools::md5sum(file.path(dir,'fit.rds'));unlink(file.path(dir,'metrics.csv'))
 e$fit_mixedgp_competitor<-function(...)stop('Unexpected refitting')
 e$baseline_cell(t,inp,out,s,TRUE)
 testthat::expect_identical(tools::md5sum(file.path(dir,'fit.rds')),model)
 testthat::expect_equal(read.csv(file.path(dir,'metrics.csv'))$RMSE_SD,expected)
})
