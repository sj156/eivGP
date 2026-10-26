repo<-Sys.getenv('EIVGP_TEST_REPO');stopifnot(nzchar(repo))
source(file.path(repo,'tools/ocean-tools.R'),local=TRUE)

test_that('empty CC summaries export a readable empty pooled table', {
  root<-tempfile();on.exit(unlink(root,recursive=TRUE),add=TRUE)
  inputs<-ocean_prepare_inputs(file.path(repo,'data/ocean_silicate/n200'))
  legacy<-list(fold_metrics=data.frame(fold=integer(),variant=character(),
    cell=character(),prediction_mode=character()))
  z<-ocean_cc_prediction_tables(root,inputs,legacy,profile='quick',expected_folds=1L)
  expect_equal(nrow(z$pooled_cc),0L)
  saved<-read.csv(file.path(z$dest,'cc_pooled_rmse.csv'))
  expect_true(all(c('cell','method','pooled_RMSE')%in%names(saved)))
})

test_that('verified CC results remain reportable without EIV predictions', {
  root<-tempfile();on.exit(unlink(root,recursive=TRUE),add=TRUE)
  ocean_storage_init(root)
  inputs<-ocean_prepare_inputs(file.path(repo,'data/ocean_silicate/n200'))
  d<-inputs$folds[['1/O20']]
  tasks<-ocean_task_catalogue('quick')
  tasks<-tasks[tasks$module=='complete_case_baselines' & tasks$cell=='O20',]
  for(i in seq_len(nrow(tasks))) {
    task<-tasks[i,,drop=FALSE]
    files<-ocean_task_files(task,'quick','fit')
    for(p in files) {
      dest<-ocean_path(root,p)
      if(grepl('[.]rds$',p))mixedgp_atomic_rds(list(),dest)
      else mixedgp_atomic_csv(data.frame(value=1),dest)
    }
    dir<-dirname(ocean_path(root,files[1L]))
    design<-list(protocol='native-CC-observed-u-subset-v1',profile='quick',
      training_ids=d$train$cohort_id[d$observed_train],
      test_ids=d$test$cohort_id[d$observed_test_mixed],controls=ocean_cc_controls('quick'))
    mixedgp_atomic_rds(design,file.path(dir,'design.rds'))
    metrics<-data.frame(fold=1L,cell='O20',method=paste0(task$method,'(CC)'),
      n_train=sum(d$observed_train),n_test=sum(d$observed_test_mixed),
      RMSE=i,CRPS=i/2,Coverage95=.95,Width95=2*i)
    mixedgp_atomic_csv(metrics,file.path(dir,'metrics.csv'))
    mixedgp_atomic_csv(data.frame(method=task$method,optimization_status='converged'),
      file.path(dir,'package_status.csv'))
    mixedgp_atomic_rds(list(files=files,hashes=ocean_hashes(ocean_path(root,files))),
      ocean_path(root,'workflow_tasks',paste0(task$task_id,'.done.rds')))
  }
  legacy<-list(fold_metrics=data.frame(fold=integer(),variant=character(),
    cell=character(),prediction_mode=character()))
  z<-ocean_cc_prediction_tables(root,inputs,legacy,cc_output_dir=root,
    profile='quick',expected_folds=1L)
  expect_equal(nrow(z$cc),3L)
  expect_equal(nrow(z$pooled_cc),3L)
  expect_setequal(z$pooled_cc$method,paste0(tasks$method,'(CC)'))
  expect_false(any(z$pooled_cc$method=='EIV-GP'))
  expect_equal(z$pooled_cc$N_total,rep(sum(d$observed_test_mixed),3L))
})
