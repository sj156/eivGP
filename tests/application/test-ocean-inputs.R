repo<-Sys.getenv('EIVGP_TEST_REPO');stopifnot(nzchar(repo))
source(file.path(repo,'tools/ocean-tools.R'),local=TRUE)
testthat::test_that("frozen Ocean cohorts retain matched observation protocols", {
  root<-normalizePath(testthat::test_path("../.."),mustWork=TRUE)
  for(spec in list(c("ocean_silicate/n200",200,160,40,32,80,8,20),
                  c("ocean_silicate/n500",500,400,100,80,200,20,50))) {
    dir<-file.path(root,"data",spec[1]);testthat::skip_if_not(dir.exists(dir))
    z<-ocean_prepare_inputs(dir);v<-as.integer(spec[-1])
    testthat::expect_equal(c(z$design$N,z$design$n_train,z$design$n_test),v[1:3])
    for(f in 1:5) {
      d20<-z$folds[[paste(f,"O20",sep="/")]];d50<-z$folds[[paste(f,"O50",sep="/")]]
      testthat::expect_equal(c(sum(d20$observed_train),sum(d50$observed_train),
        sum(d20$observed_test_mixed),sum(d50$observed_test_mixed)),v[4:7])
      testthat::expect_true(all(d20$observed_train<=d50$observed_train))
      testthat::expect_true(all(d20$observed_test_mixed<=d50$observed_test_mixed))
      if(z$design$N==500L)testthat::expect_equal(sort(as.integer(table(d50$test$C[d50$observed_test_mixed]))),c(12L,12L,13L,13L))
    }
  }
})

testthat::test_that("Ocean validation rejects a changed held-out response or frozen mask", {
  root<-normalizePath(testthat::test_path("../.."),mustWork=TRUE)
  source<-file.path(root,"data/ocean_silicate/n500");testthat::skip_if_not(dir.exists(source))
  target<-tempfile("ocean-input-audit-");dir.create(target)
  on.exit(unlink(target,recursive=TRUE),add=TRUE)
  files<-list.files(source,"[.]csv$",full.names=TRUE);file.copy(files,target)
  file.copy(file.path(source,"freeze_manifest.json"),target)
  testthat::expect_silent(ocean_prepare_inputs(target))
  p<-file.path(target,"fold01_test.csv");d<-read.csv(p);d$Y[1]<-d$Y[1]+1
  write.csv(d,p,row.names=FALSE)
  testthat::expect_error(ocean_prepare_inputs(target))
  file.copy(file.path(source,"fold01_test.csv"),p,overwrite=TRUE)
  p<-file.path(target,"fold01_test_U_masks.csv");d<-read.csv(p);d$observed_U_O20[1]<-!d$observed_U_O20[1]
  write.csv(d,p,row.names=FALSE)
  testthat::expect_error(ocean_prepare_inputs(target))
})

testthat::test_that("the 1000-profile source is an unpartitioned cohort", {
  root<-normalizePath(testthat::test_path("../.."),mustWork=TRUE)
  path<-file.path(root,"data/ocean_silicate/n1000")
  testthat::skip_if_not(dir.exists(path))
  cohort<-ocean_read_inputs(path)
  testthat::expect_equal(nrow(cohort),1000L)
  testthat::expect_equal(as.integer(table(cohort$C)),rep(250L,4))
  testthat::expect_false("fold" %in% names(cohort))
  testthat::expect_length(list.files(path,"^fold"),0L)
})

testthat::test_that("Ocean inputs require the declared approximate Latin protocol", {
  source<-file.path(repo,'data/ocean_silicate/n500')
  target<-tempfile('ocean-calibration-contract-');dir.create(target)
  on.exit(unlink(target,recursive=TRUE),add=TRUE)
  files<-list.files(source,full.names=TRUE);file.copy(files,target)
  manifest<-file.path(target,'freeze_manifest.json')
  metadata<-jsonlite::read_json(manifest)
  metadata$calibration_design$method<-'unsupported selector'
  jsonlite::write_json(metadata,manifest,auto_unbox=TRUE)
  testthat::expect_error(ocean_prepare_inputs(target),'class-stratified approximate Latin')
  unlink(manifest)
  testthat::expect_error(ocean_prepare_inputs(target),'Missing freeze_manifest.json')
})
