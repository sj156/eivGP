testthat::test_that("R fold preparation is deterministic and independent of responses", {
  root<-normalizePath(testthat::test_path("../../.."),mustWork=TRUE)
  script<-file.path(root,"data/ocean_silicate/prepare_folds.R")
  testthat::skip_if_not(file.exists(script))
  entry<-new.env(parent=globalenv());sys.source(script,entry)
  for(n in c(500L,1000L)) {
    path<-file.path(root,"data/ocean_silicate",paste0("n",n),"cohort_manifest.csv")
    cohort<-read.csv(path,stringsAsFactors=FALSE);cohort$fold<-NULL
    rule<-if(n==500L)"joint_quintile_pressure" else "oxygen_interleaved"
    sample_seed<-if(n==500L)465000L else 144069L
    seed<-sample_seed+if(n==500L)100L else 120L
    a<-entry$ocean_assign_folds(cohort,rule,seed)
    no_y<-cohort[,tolower(names(cohort))!="y",drop=FALSE]
    b<-entry$ocean_assign_folds(no_y,rule,seed)
    changed<-cohort
    for(name in names(changed)[tolower(names(changed))=="y"])changed[[name]]<-rev(seq_len(n))
    c<-entry$ocean_assign_folds(changed,rule,seed)
    testthat::expect_identical(a$cohort$fold,b$cohort$fold)
    testthat::expect_identical(a$cohort$fold,c$cohort$fold)
    testthat::expect_true(all(table(a$cohort$fold,a$cohort$C)==n/20L))
    testthat::expect_lte(a$design$final_design_score,a$design$initial_design_score)
    held<-character()
    for(f in 1:5) {
      tr<-a$cohort[a$cohort$fold!=f,,drop=FALSE]
      train<-entry$ocean_fold_calibration(tr,sample_seed+200L+f)
      tr_no_y<-tr[,tolower(names(tr))!="y",drop=FALSE]
      other<-entry$ocean_fold_calibration(tr_no_y,sample_seed+200L+f)
      test<-a$cohort[a$cohort$fold==f,,drop=FALSE]
      masks<-entry$ocean_fold_test_masks(test,f)
      testthat::expect_identical(train$calib_o20,other$calib_o20)
      testthat::expect_identical(train$calib_o50,other$calib_o50)
      testthat::expect_length(intersect(train$profile,test$profile),0L)
      testthat::expect_equal(as.integer(tapply(train$calib_o20,train$C,sum)),rep(n/25L,4))
      testthat::expect_equal(as.integer(tapply(train$calib_o50,train$C,sum)),rep(n/10L,4))
      testthat::expect_true(all(train$calib_o20<=train$calib_o50))
      testthat::expect_true(all(tapply(train$calib_o20,train$cell,sum)>=n/250L))
      testthat::expect_true(all(tapply(train$calib_o50,train$cell,sum)>=3L*n/250L))
      testthat::expect_equal(sum(masks$observed_U_O20),n/25L)
      testthat::expect_equal(sum(masks$observed_U_O50),n/10L)
      testthat::expect_true(all(masks$observed_U_O20<=masks$observed_U_O50))
      testthat::expect_identical(masks,ocean_test_masks(test,f))
      held<-c(held,test$cohort_id)
    }
    testthat::expect_setequal(held,cohort$cohort_id)
    testthat::expect_false(anyDuplicated(held)>0L)
  }
})
