# Author: Xinye Chen. Synthetic configuration, privacy and observation-design contracts.
repo<-Sys.getenv('EIVGP_TEST_REPO');stopifnot(nzchar(repo))
source(file.path(repo,'paper/applications/adni.R'))
list2env(as.list(adni_load_api(repo)),envir=environment())
source(file.path(repo,'data/adni/prepare_folds.R'))

testthat::test_that('application files correspond to Ocean and implementation is installed',{
 app<-file.path(repo,'paper/applications')
 ocean<-list.files(app,'^ocean([.-].*)?[.](R|Rmd|html)$')
 adni<-list.files(app,'^adni([.-].*)?[.](R|Rmd|html)$')
 testthat::expect_setequal(adni,sub('^ocean','adni',ocean))
 testthat::expect_equal(sum(grepl('[.]R$',adni)),5L)
 testthat::expect_true(all(c('adni_workflow_api','run_adni_application')%in%getNamespaceExports('eivGP')))
 api<-eivGP::adni_workflow_api(repo)
 for(n in c('adni_method_run','adni_plan_fit_task','adni_analyze_cell','adni_report_prepare','adni_saved_functionals'))
  testthat::expect_identical(body(api[[n]]),body(getFromNamespace(n,'eivGP')))
})

testthat::test_that('ADNI pilot matches screening and paper accepts explicit budgets',{
 p<-adni_settings('pilot',cores=5L)
 testthat::expect_identical(p$chains,4L);testthat::expect_identical(p$stages,c(1000L,1500L))
 testthat::expect_identical(p$iter,1500L);testthat::expect_identical(p$burn,500L)
 testthat::expect_identical(p$cells,c('O50','O20','O0'))
 testthat::expect_identical(p$draws_per_chain,250L);testthat::expect_identical(p$baseline_profile,'paper')
 testthat::expect_identical(p$sampler_o0,'HMC');testthat::expect_identical(p$sampler_o20,'HMC');testthat::expect_identical(p$sampler_o50,'ESS')
 paper<-adni_settings('paper',cores=12L)
 testthat::expect_identical(paper$chains,4L);testthat::expect_identical(paper$iter,50000L);testthat::expect_identical(paper$burn,5000L)
 custom<-adni_settings('paper',iter=20000L,burn=2000L,stages='10000,20000',chains=3L)
 testthat::expect_identical(custom$stages,c(10000L,20000L));testthat::expect_gt(custom$first_checkpoint,custom$burn)
 testthat::expect_error(adni_settings('paper',iter=1000L,burn=2000L),'exceed burn')
 testthat::expect_error(adni_settings('pilot',first_checkpoint=500L),'first-checkpoint')
 testthat::expect_error(adni_cli_params(c('--cores=2','--cores=3')),'Repeated')
 testthat::expect_error(adni_cli_params('--unknown=1'),'Unknown')
})

testthat::test_that('arbitrary C dimensions prepare portable fivefold inputs without outcome selection',{
 root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
 set.seed(306L);n<-100L
 d<-data.frame(RID=seq_len(n),Y=rnorm(n),X1=60+runif(n)*20,X2=10+runif(n)*12,U1=rnorm(n),U2=rnorm(n))
 for(j in 1:9)d[[paste0('C',j)]]<-sample(rep(1:4,length.out=n))
 for(q in c(3L,9L)) {
  meta<-list(dataset_id=paste0('synthetic',q),C_order=paste0('native_item',seq_len(q)),m_vec=rep(4L,q),
    X_order=c('age','education'),U_order=c('FDG','HIP'),fold_seed=9L,calibration_seed=10L,test_seed=11L,model_seed=12L,stratification_columns='C1')
  csv<-file.path(root,paste0('q',q,'.csv'));json<-file.path(root,paste0('q',q,'.json'));out<-file.path(root,paste0('q',q))
  write.csv(d[,c('RID','Y','X1','X2','U1','U2',paste0('C',seq_len(q)))],csv,row.names=FALSE)
  jsonlite::write_json(meta,json,auto_unbox=TRUE)
  mapping<-adni_prepare_folds(csv,json,out)
  api<-adni_screen_api(repo);inputs<-api$adni_plan_inputs(out)
  testthat::expect_equal(inputs$N,n);testthat::expect_equal(length(inputs$cells),15L)
  testthat::expect_true(all(inputs$audit$q==q));testthat::expect_true(all(table(mapping$fold)==20L))
  for(f in 1:5){x<-inputs$cells[[paste(f,'O20',sep='/')]];y<-inputs$cells[[paste(f,'O50',sep='/')]]
   testthat::expect_true(all(x$train$calibration<=y$train$calibration));testthat::expect_true(all(!x$available|y$available))
   counts<-table(x$train$C1);selected<-table(factor(x$train$C1[x$train$calibration==1],levels=names(counts)))
   testthat::expect_lte(max(abs(as.numeric(selected)-sum(selected)*as.numeric(counts)/sum(counts))),1)
  }
  changed<-d;changed$Y<-1e8;changed$U1<--1e8;changed$U2<-1e8
  testthat::expect_identical(adni_balanced_folds(d,paste0('C',seq_len(q)),c('X1','X2'),9L),
                           adni_balanced_folds(changed,paste0('C',seq_len(q)),c('X1','X2'),9L))
  testthat::expect_error(adni_prepare_folds(csv,json,out),'never overwrite')
  s<-adni_screen_settings(list(profile='paper',chains=4L,burn=2000L,folds=1:5,cells=c('O0','O20','O50')),20000L,4L)
  args<-api$adni_plan_fit_args(inputs$cells[['1/O20']],s)
  testthat::expect_length(args$ordinal_levels,q);testthat::expect_true(all(is.na(args$U_obs[inputs$cells[['1/O20']]$train$calibration==0,,drop=FALSE])))
 }
})

testthat::test_that('nested strata handle singleton groups and prohibit source-directory output',{
 d<-data.frame(RID=1:10,C1=c(1,rep(2,3),rep(3,6)))
 masks<-adni_nested_masks(d,'C1',991L)
 testthat::expect_true(all(masks$O20%in%masks$O50));testthat::expect_true(all(masks$O50%in%d$RID))
 testthat::expect_error(adni_output_check(file.path(repo,'paper/results'),repo),'ignored results')
 testthat::expect_silent(adni_output_check(file.path(repo,'results/adni'),repo))
 testthat::expect_false(file.exists(file.path(repo,'data/adni/inputs/MC/fold01/O0/private/train.csv')))
})

testthat::test_that('general entrypoint freezes custom settings and resumes with a different core count',{
 root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
 set.seed(630L);n<-40L
 d<-data.frame(RID=1:n,Y=rnorm(n),X1=60+runif(n)*20,X2=10+runif(n)*12,U1=rnorm(n),U2=rnorm(n),
   C1=sample(rep(1:4,10)),C2=sample(rep(1:4,10)),C3=sample(rep(1:4,10)))
 meta<-list(dataset_id='synthetic_general',C_order=c('A','B','C'),m_vec=rep(4L,3),X_order=c('age','education'),
   U_order=c('FDG','HIP'),fold_seed=9L,calibration_seed=10L,test_seed=11L,model_seed=12L)
 write.csv(d,file.path(root,'cohort.csv'),row.names=FALSE);jsonlite::write_json(meta,file.path(root,'meta.json'),auto_unbox=TRUE)
 input<-file.path(root,'inputs');out<-file.path(root,'results');adni_prepare_folds(file.path(root,'cohort.csv'),file.path(root,'meta.json'),input)
 e<-eivGP::adni_workflow_api(repo)
 # The runner is exercised; expensive HTML functionals are tested separately from saved files.
 e$adni_render_report<-function(...)invisible(NULL)
 args<-list(profile='quick',cores=2L,data_root=input,output=out,iter=32L,burn=8L,stages=c(24L,32L),
   folds=1L,cells='O20',chains=4L,draws_per_chain=5L)
 e$adni_method_run(repo,args,'eivGP')
 manifest<-jsonlite::fromJSON(file.path(out,'application_manifest.json'),simplifyVector=FALSE)
 testthat::expect_identical(unlist(manifest$stages),c(24L,32L));testthat::expect_identical(manifest$chains,4L)
 fit<-file.path(out,'fits/fold01/O20/fit.rds');before<-tools::md5sum(fit)
 args$cores<-1L;testthat::expect_error(e$adni_method_run(repo,args,'eivGP'),NA)
 testthat::expect_identical(tools::md5sum(fit),before)
 args$burn<-9L;testthat::expect_error(e$adni_method_run(repo,args,'eivGP'),'settings/inputs differ')
 saved<-e$adni_saved_source(out)
 testthat::expect_identical(saved$stage,'S32');testthat::expect_identical(saved$fit_profile,'quick')
 testthat::expect_error(e$adni_report_prepare(repo,input,out,adni_settings('pilot'),report_dir=file.path(repo,'paper/invalid')),'ignored results')
})
