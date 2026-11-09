# Portable planner/manifest checks without model fitting or study result inputs.
repo<-Sys.getenv('EIVGP_TEST_REPO');stopifnot(nzchar(repo))
e<-new.env(parent=globalenv());sys.source(file.path(repo,'tools/study2-followup-plan.R'),e)
testthat::test_that('follow-up planning discovers one run and exports portable pilot inputs', {
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  run<-file.path(root,'results/study2-eivGP/core/run-a')
  dir.create(file.path(run,'combined'),recursive=TRUE)
  fits<-file.path(run,'cells/primary_q4_calibration');dir.create(fits,recursive=TRUE)
  data<-file.path(root,'results/data/study2/primary_q4_calibration');dir.create(data,recursive=TRUE)
  d<-expand.grid(rep=1:50,n_calib=c(0L,20L))
  d$cell_id<-'primary_q4_calibration'
  for(k in c('mcmc_pass','target_functional_pass','invariant_measurement_pass','raw_coordinate_pass'))d[[k]]<-FALSE
  write.csv(d,file.path(run,'combined/mcmc_diagnostics.csv'),row.names=FALSE)
  for(i in seq_len(nrow(d)))saveRDS(list(fixture=TRUE),file.path(fits,
    sprintf('recovery_fit_primary_rep%03d_cal%03d_fixture.rds',d$rep[i],d$n_calib[i])))
  for(r in 1:50)saveRDS(list(fixture=TRUE),file.path(data,
    sprintf('study2_primary_q4_n100_ntest100_m4_rep%03d.rds',r)))
  p<-e$study2_followup_plan(root)
  testthat::expect_equal(nrow(p$jobs),10L)
  testthat::expect_identical(p$source_run,normalizePath(run))
  testthat::expect_false(p$execution_enabled)
  testthat::expect_true(all(file.exists(p$jobs$data)))
  testthat::expect_equal(nrow(e$study2_followup_plan(root,scope='all')$jobs),100L)
  # Ambiguous run discovery requires a deliberate source selection.
  second<-file.path(dirname(run),'run-b/combined');dir.create(second,recursive=TRUE)
  file.copy(file.path(run,'combined/mcmc_diagnostics.csv'),second)
  testthat::expect_error(e$study2_followup_plan(root),'Supply --source')
  output<-file.path(root,'pilot-inputs')
  args<-c(file.path(repo,'tools/study2-followup-plan.R'),paste0('--repo=',root),
    paste0('--source=',run),paste0('--data-root=',dirname(data)),paste0('--output=',output))
  log<-file.path(root,'plan.log')
  status<-system2(file.path(R.home('bin'),'Rscript'),shQuote(args),stdout=log,stderr=log)
  testthat::expect_equal(status,0L)
  csv<-file.path(output,'plan-pilot.csv');rds<-file.path(output,'followup-plan-pilot.rds')
  testthat::expect_true(all(file.exists(c(csv,rds))))
  manifest<-read.csv(csv,stringsAsFactors=FALSE)
  testthat::expect_equal(nrow(manifest),10L)
  testthat::expect_true(all(c('rep','calibration','fit','data')%in%names(manifest)))
  before<-tools::md5sum(c(csv,rds))
  status<-suppressWarnings(system2(file.path(R.home('bin'),'Rscript'),shQuote(args),stdout=log,stderr=log))
  testthat::expect_gt(status,0L)
  testthat::expect_identical(tools::md5sum(names(before)),before)
})
