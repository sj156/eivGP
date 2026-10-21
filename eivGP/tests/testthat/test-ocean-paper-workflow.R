test_that("cores-only CLI requests the full workflow and no arguments remain a plan", {
  repo<-normalizePath(testthat::test_path("../../.."))
  env<-new.env(parent=globalenv());sys.source(file.path(repo,"paper/applications/ocean.R"),env)
  expect_identical(env$ocean_cli_params('--cores=16'),list(cores=16L,action='fit'))
  expect_null(env$ocean_cli_params(character())$action)
  expect_identical(env$ocean_cli_params(c('--cores=16','--action=report'))$action,'report')
  p<-ocean_stage_plan(ocean_settings(cores=16))
  expect_identical(p$stage,c('fit','predict','functionals','paper_assets','html'))
  expect_equal(ocean_stage_workers(16,'fit',15,4),4)
  expect_equal(ocean_stage_workers(16,'analyze',15,4),8)
  expect_equal(ocean_stage_workers(1,'analyze',15,4),1)
})

test_that("comparator failure does not discard successful tasks or block postprocessing", {
  repo<-normalizePath(testthat::test_path("../../.."));data<-file.path(repo,'data/ocean_silicate/n200')
  root<-tempfile();on.exit(unlink(root,recursive=TRUE),add=TRUE)
  env<-new.env(parent=environment(ocean_execute));execute<-ocean_execute;environment(execute)<-env
  env$ocean_identity<-function(...)list(profile='pilot',code='fixture')
  env$ocean_run_task<-function(task,repository_dir,output_dir,settings,inputs,identity) {
    if(task$module=='original_baselines')stop('LVGP timeout')
    for(p in ocean_task_files(task,settings$profile,settings$action)) {
      if(grepl('[.]rds$',p))ocean_atomic_rds(list(chains=1:4),ocean_path(output_dir,p)) else
        ocean_atomic_csv(data.frame(rhat=1,ess_bulk=500,ess_tail=500),ocean_path(output_dir,p))
    }
  }
  cfg<-ocean_settings(profile='pilot',action='fit',folds=1L,cells='O20',modules=c('eiv','original_baselines'),cores=1L)
  expect_warning(execute(repo,data,root,cfg),'Comparator tasks failed')
  expect_true(ocean_verify_receipt(root,'fit_eiv_fold01_O20',env$ocean_identity()))
  status<-read.csv(file.path(ocean_job_dir(root,cfg$run_id),'task_status.csv'))
  expect_true(any(status$state=='FAILED'));expect_true(any(status$state=='DONE'))
  cfg$action<-'report';expect_silent(execute(repo,data,root,cfg))
})

test_that("failed native baselines are explicit and never replaced with curves", {
  root<-tempfile();on.exit(unlink(root,recursive=TRUE),add=TRUE);ocean_storage_init(root)
  p<-ocean_path(root,'full5_class4','fold01','baseline_status.csv')
  ocean_atomic_csv(data.frame(method=c('UC-GP','LVGP','EzGP'),status=c('success','unavailable_or_failed','success')),p)
  expect_equal(ocean_available_baselines(root,1),c('UC-GP','EzGP'))
  ocean_atomic_csv(data.frame(method=c('UC-GP','LVGP','EzGP'),status=rep('unavailable_or_failed',3)),p)
  ocean_atomic_rds(list(fits=list()),ocean_path(root,'full5_class4','fold01','baseline_package_outputs.rds'))
  grid<-data.frame(G2salinity=34,log_pressure=6,C=1,design='same_x')
  empty<-ocean_baseline_mean_grid(root,grid,1)
  expect_equal(nrow(empty),0);expect_true(all(c(names(grid),'posterior_mean','fold','model')%in%names(empty)))
})

test_that("paper exports retain units, fold counts and manuscript inventory", {
  registry<-ocean_paper_registry()
  expect_equal(sum(registry$kind=='figure'),8);expect_equal(sum(registry$kind=='table'),6)
  expect_false(anyDuplicated(registry$file)>0);expect_false(anyDuplicated(registry$label)>0)
  z<-data.frame(method='EIV-GP',cell='O20',prediction_mode='mixed',n_folds=5)
  for(m in c('RMSE','MAE','CRPS','NLPD','Coverage95','Width95','IntervalScore95')) {
    z[[paste0(m,'_mean')]]<-2;z[[paste0(m,'_sd')]]<-.2
  }
  z$Coverage95_mean<-.95;z$Coverage95_sd<-.02
  tab<-ocean_paper_prediction_table(z,10,TRUE)
  expect_equal(tab$Folds,5);expect_equal(tab$RMSE,'2.00 (0.20)')
  expect_equal(tab[['Coverage (%)']],'95.00 (2.00)')
  expect_equal(tab[['CRPS/sY']],'0.200 (0.020)')
  expect_equal(tab$NLPD_Z,sprintf('%.3f (0.200)',2-log(10)))
  p<-ocean_cohort_plot(data.frame(U=c(20,70,120,170),Y=1:4,G2salinity=34:37,log_pressure=6:9,C=1:4))
  expect_s3_class(p,'patchwork')
})

test_that("partial comparator summaries disclose fold counts without weakening imputation checks", {
  d<-data.frame(fold=c(1L,2L,4L,5L),method='LVGP',RMSE=c(1,2,3,4))
  expect_error(ocean_diag_fold_summary(d,'method','RMSE'),'one score per fold')
  z<-ocean_diag_fold_summary(d,'method','RMSE',require_all=FALSE)
  expect_equal(z$n_folds,4);expect_equal(z$mean,2.5)
  expect_error(ocean_diag_fold_summary(rbind(d,d[1,]),'method','RMSE',require_all=FALSE),'one score per fold')
  expect_equal(ocean_stage_workers(16,'analyze',0),1)
})

test_that("plot export uses a fork-safe raster device even with Quartz as the default", {
  skip_if_not(capabilities('cairo'))
  root<-tempfile();on.exit(unlink(root,recursive=TRUE),add=TRUE)
  old<-options(bitmapType='quartz');on.exit(options(old),add=TRUE)
  p<-ocean_cohort_plot(data.frame(U=c(20,70,120,170),Y=1:4,G2salinity=34:37,log_pressure=6:9,C=1:4))
  release<-mixedgp_budget_enter(2L);on.exit(release(),add=TRUE)
  mixedgp_save_plot_jobs(list(cohort=p),root,function(n)c(8,3))
  expect_true(all(file.exists(file.path(root,c('cohort.png','cohort.pdf')))))
  expect_true(all(file.info(file.path(root,c('cohort.png','cohort.pdf')))$size>1000))
})

test_that("appendix panels export their manuscript cross-reference and preserve raw columns", {
  root<-tempfile();on.exit(unlink(root,recursive=TRUE),add=TRUE);dir.create(root)
  d<-data.frame(Method='EIV-GP',Calibration='20%',Test='Class only',Folds=5,
    RMSE='1 (0.1)',check.names=FALSE)
  for(n in c('RMSE/sY','MAE/sY','CRPS/sY','NLPD_Z','Coverage (%)','Width/sY','Interval score/sY'))d[[n]]<-'0.2 (0.01)'
  ocean_export_prediction_panels(d,root,'test')
  tex<-readLines(file.path(root,'tab_ocean_prediction_full.tex'))
  expect_true(any(grepl('\\label{tab:ocean-si-prediction-full}',tex,fixed=TRUE)))
  expect_true(file.exists(file.path(root,'tab_ocean_prediction_full.html')))
})

test_that("HTML receives completed analysis instead of rerunning computation", {
  repo<-normalizePath(testthat::test_path('../../..'));env<-new.env(parent=globalenv())
  sys.source(file.path(repo,'paper/applications/ocean.R'),env)
  seen<-NULL;calls<-0L
  execute<-function(s,import) {calls<<-calls+1L;list(stage=s$action,token=42L)}
  render<-function(analysis) {seen<<-analysis;'ocean.html'}
  out<-env$ocean_cli_workflow(ocean_settings(profile='pilot',action='fit'),execute,render)
  expect_equal(calls,2);expect_identical(seen$out,NULL)
  expect_identical(seen$stage,'analyze');expect_identical(seen$token,42L)
  expect_identical(out$html,'ocean.html')
})
