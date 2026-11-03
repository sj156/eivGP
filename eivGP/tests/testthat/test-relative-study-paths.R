test_that('new relative study roots survive the fitting working-directory change', {
  root<-tempfile('study-paths-');dir.create(root)
  old<-getwd();on.exit({setwd(old);unlink(root,recursive=TRUE)},add=TRUE)
  setwd(root)
  for(study in c('study1','study2')) {
    cfg<-study_config(study,'quick',output_root=file.path('results',study),
      data_root=file.path('results','data',study),cores=1L)
    expect_true(startsWith(cfg$data_root,paste0(normalizePath(root),'/')))
    expect_true(startsWith(cfg$output_root,paste0(normalizePath(root),'/')))
    cfg$cells<-cfg$cells[1L]
    run_numerical_study(cfg,'data',competitors=character())
    cell<-cfg$cells[[1L]]
    files<-list.files(mixedgp_data_cell_directory(cfg,cell),'[.]rds$',full.names=TRUE)
    expect_true(length(files)>0L)
    hashes<-tools::md5sum(files)
    # Exactly the directory switch performed by the real fitting driver.
    setwd(cfg$code_dir)
    controls<-get(paste0('mixedgp_cell_controls_',study))(cfg,cell,
      file.path(mixedgp_run_directory(cfg),'cells',cell$id))
    data_dir<-controls[[paste0(toupper(study),'_DATA_DIR')]]
    expect_identical(data_dir,mixedgp_data_cell_directory(cfg,cell))
    expect_true(all(file.exists(file.path(data_dir,basename(files)))))
    expect_identical(tools::md5sum(file.path(data_dir,basename(files))),hashes)
    # Loading frozen data for fitting requires no MCMC or regenerated inputs.
    expect_identical(readRDS(files[1L])$study,study)
    setwd(root)
    run_numerical_study(cfg,'data',competitors=character())
    expect_identical(tools::md5sum(files),hashes)
  }
})
