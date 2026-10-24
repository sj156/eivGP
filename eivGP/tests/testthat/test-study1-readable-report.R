test_that('Study I reports scientific ESS strictly above 100 without using R-hat as a gate', {
  z<-data.frame(cell_id='eta0_balanced',n_calib=0,rep=1:4,family='m_conditional_mean',
    ess_bulk=c(101,100,50,NA),ess_tail=c(150,200,200,200),rhat=c(1.2,1,1,1),mcse_sd_ratio=.1)
  x<-mixedgp_study1_target_ess_summary(z,100)
  expect_equal(x$below_ess,2);expect_equal(x$unavailable,1);expect_equal(x$affected_fits,3)
  expect_equal(x$max_rhat,1.2)
})

test_that('Study I saved reporting works without fitting, frozen inputs or retained fits when traces are disabled', {
  root<-tempfile();out<-tempfile();dir.create(root);dir.create(file.path(root,'config'))
  on.exit(unlink(c(root,out),recursive=TRUE))
  saveRDS(list(study='study1',data_root='/missing'),file.path(root,'config/resolved_config.rds'))
  tab<-file.path(root,'cells/eta0_balanced/tables');dir.create(tab,recursive=TRUE)
  d<-data.frame(parameter=c('m_conditional_mean[1]','m_conditional_variance[1]'),rep=1,n_calib=0,
    ess_bulk=c(101,100),ess_tail=150,rhat=1.2,mcse_sd_ratio=.1)
  write.csv(d,file.path(tab,'study1_mcmc_target_diagnostics.csv'),row.names=FALSE)
  x<-report_study1_eivgp(root,out,mixing_plots=FALSE)
  expect_equal(nrow(x$review),1);expect_identical(x$targets$ess_acceptable,c(TRUE,FALSE))
  expect_true(file.exists(file.path(out,'report.rds')));expect_length(x$figures$paths,0)
  expect_false(dir.exists('/missing'));expect_equal(x$run_dir,normalizePath(root))
})

test_that('Study I report refuses ambiguous saved runs', {
  root<-tempfile();on.exit(unlink(root,recursive=TRUE))
  for(n in c('a','b')) {
    p<-file.path(root,n,'config');dir.create(p,recursive=TRUE)
    saveRDS(list(study='study1'),file.path(p,'resolved_config.rds'))
  }
  expect_error(mixedgp_study1_report_run_dir(root),'exact Study I')
})

test_that('summary graphics are serial even under a larger shared core budget', {
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  write.csv(data.frame(parameter='a',rhat=1,ess_bulk=200,ess_tail=200),
    file.path(root,'study1_mcmc_parameter_diagnostics.csv'),row.names=FALSE)
  f<-save_mixing_plots;e<-new.env(parent=environment(f));seen<-NULL
  e$mixedgp_stage_lapply<-function(X,FUN,...,workers) {seen<<-workers;lapply(X,FUN,...)}
  environment(f)<-e
  paths<-f(root,formats='pdf');expect_identical(seen,1L)
  expect_true(all(file.info(paths)$size>0))
})


test_that('the notebook does not discover its previous HTML output as an input resource', {
  skip_if_not_installed('rmarkdown')
  notebook<-test_path('../../../paper/numerical/study1.Rmd')
  skip_if_not(file.exists(notebook))
  resources<-rmarkdown::find_external_resources(notebook)
  expect_false('study1.html' %in% basename(resources$path))
})
