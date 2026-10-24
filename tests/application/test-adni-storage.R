# Author: Xinye Chen. Flat task storage, frozen inputs and legacy reader contracts.
repo<-Sys.getenv('EIVGP_TEST_REPO');stopifnot(nzchar(repo))
source(file.path(repo,'tools/adni.R'))
list2env(as.list(adni_load_api(repo)),envir=environment())

testthat::test_that('new ADNI roots use Ocean task and shared-directory conventions',{
 out<-tempfile();on.exit(unlink(out,recursive=TRUE))
 adni_storage_init(out,c(1000L,1500L))
 testthat::expect_true(adni_storage_flat(out))
 pairs<-c('fits/fold01/O20/fit.rds'='eiv_fold01_O20/fit.rds',
  'S1500/analysis/fold01/O20/prediction_metrics.csv'='eiv_fold01_O20/prediction_metrics.csv',
  'S1000/analysis/fold01/O20/prediction_metrics.csv'='eiv_fold01_O20/stages/S1000/prediction_metrics.csv',
  'baselines/fold01/LVGP/fit.rds'='baseline_fold01/LVGP/fit.rds',
  'CC/O20/baselines/fold01/UC-GP/fit.rds'='cc_fold01_O20_UCGP/fit.rds',
  'workflow/logs/run/fit_status.csv'='logs/run/fit_status.csv',
  'report-assets/paper/artifact_manifest.csv'='reports/paper/artifact_manifest.csv',
  'report-assets/functionals/fold1_O20.rds'='reports/diagnostics/functionals/fold1_O20.rds',
  'receipts/done.rds'='_meta/tasks/done.rds',
  'application_manifest.json'='_meta/application_manifest.json',
  'workflow_identity.rds'='_meta/workflow_identity.rds')
 testthat::expect_identical(adni_path(out,names(pairs)),unname(file.path(out,pairs)))
 testthat::expect_identical(adni_path(out,unname(pairs)),unname(file.path(out,pairs)))
 testthat::expect_identical(adni_workflow_root(out),out)
 testthat::expect_identical(adni_default_output('adni','paper',repo),file.path(repo,'results/adni'))
 testthat::expect_identical(adni_default_output('3C','paper',repo),file.path(repo,'results/3C'))
 testthat::expect_identical(adni_default_output('adni','pilot',repo),file.path(repo,'results/adni_pilot'))
})

testthat::test_that('unmarked legacy roots retain their paths and payloads',{
 out<-tempfile();dir.create(out);on.exit(unlink(out,recursive=TRUE))
 path<-file.path(out,'fits/fold01/O20/fit.rds');dir.create(dirname(path),recursive=TRUE)
 saveRDS(list(original=TRUE),path);before<-tools::md5sum(path)
 adni_storage_init(out,1500L)
 testthat::expect_false(adni_storage_flat(out))
 testthat::expect_identical(adni_plan_fit_path(out,list(fold=1L,cell='O20')),path)
 testthat::expect_identical(tools::md5sum(path),before)
 testthat::expect_false(dir.exists(file.path(out,'_meta')))
 testthat::expect_identical(adni_workflow_root(out),file.path(out,'workflow'))
})

testthat::test_that('frozen snapshots are byte-preserving and reject changed inputs',{
 root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
 input<-file.path(root,'inputs');dir.create(file.path(input,'MC/fold01/O20/private'),recursive=TRUE)
 csv<-file.path(input,'MC/fold01/O20/private/train.csv');writeLines(c('RID,Y','1,2'),csv)
 jsonlite::write_json(list(dataset_id='synthetic'),file.path(input,'DATASET.json'))
 out<-file.path(root,'results');adni_storage_init(out,1500L)
 dest<-adni_snapshot_inputs(input,out,csv);copy<-file.path(dest,'fold01_O20_train.csv')
 testthat::expect_identical(readBin(csv,'raw',n=file.info(csv)$size),readBin(copy,'raw',n=file.info(copy)$size))
 before<-tools::md5sum(copy)
 testthat::expect_identical(adni_snapshot_inputs(input,out,csv),dest)
 testthat::expect_true(file.exists(file.path(dest,'DATASET.json')))
 testthat::expect_false(dir.exists(file.path(dest,'MC')))
 writeLines(c('RID,Y','1,999'),csv)
 testthat::expect_error(adni_snapshot_inputs(input,out,csv),'Frozen input snapshot differs')
 testthat::expect_identical(tools::md5sum(copy),before)
})

testthat::test_that('flat and legacy saved reports resolve their own manifests',{
 root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
 for(flat in c(FALSE,TRUE)) {
  out<-file.path(root,if(flat)'flat'else'legacy');dir.create(out)
  if(flat)adni_storage_init(out,c(1000L,1500L))
  m<-list(cases=list(list(input='frozen-input')),stages=c(1000L,1500L),
    profile='pilot',package_version='synthetic',frozen_inputs=list())
  path<-adni_path(out,'application_manifest.json');dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  jsonlite::write_json(m,path,auto_unbox=TRUE)
  z<-adni_saved_source(out)
  testthat::expect_identical(z$stage,'S1500')
  testthat::expect_identical(z$input,'frozen-input')
  testthat::expect_identical(z$fitting_manifest,path)
 }
})


testthat::test_that('flat input filenames and dataset roots are portable',{
 testthat::expect_identical(adni_default_input('adni',repo),file.path(repo,'data/adni'))
 testthat::expect_identical(adni_default_input('3C',repo),file.path(repo,'data/adni_3C'))
 root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
 legacy<-file.path(root,'MC/fold01/O20/private/train.csv')
 testthat::expect_identical(adni_input_basename(legacy),'fold01_O20_train.csv')
 testthat::expect_identical(adni_input_basename(file.path(root,'MC/fold01/O20/PROTOCOL.json')),'fold01_O20_PROTOCOL.json')
 testthat::expect_identical(adni_input_basename('fold01_O20_train.csv'),'fold01_O20_train.csv')
 flat<-file.path(root,c('fold01_O20_PROTOCOL.json','fold01_O20_train.csv','fold01_O20_test.csv'))
 file.create(flat[1])
 testthat::expect_identical(adni_input_files(root,1L,'O20'),flat)
 # Partially supplied flat cells must fail validation instead of falling back to old files.
 testthat::expect_false(all(file.exists(adni_input_files(root,1L,'O20'))))
})
