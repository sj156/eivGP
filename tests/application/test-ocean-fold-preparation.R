repo<-Sys.getenv('EIVGP_TEST_REPO');stopifnot(nzchar(repo))
source(file.path(repo,'tools/ocean-tools.R'),local=TRUE)
testthat::test_that("R fold preparation is deterministic and independent of responses", {
  root<-normalizePath(testthat::test_path("../.."),mustWork=TRUE)
  script<-file.path(root,"tools/ocean-tools.R")
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
      # Replacing hidden/reference information must leave the C/X-based mask unchanged.
      recorded<-tr[,c("cohort_id","C","G2salinity","G2pressure")]
      protected<-entry$ocean_fold_calibration(recorded,sample_seed+200L+f)
      testthat::expect_identical(train$calib_o20,protected$calib_o20)
      testthat::expect_identical(train$calib_o50,protected$calib_o50)
      audit<-attr(train,"calibration_design")$audit
      testthat::expect_equal(nrow(audit),16L)
      testthat::expect_true(all(audit$occupied_bins<=audit$quantile_bins))
      testthat::expect_identical(audit$exact_latin_margin,audit$squared_occupancy_loss==0)
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

testthat::test_that("nested X coverage is reproducible, audited, and invariant to hidden truth", {
  d<-read.csv(file.path(repo,"data/ocean_silicate/n500/fold01_train.csv"))
  a<-ocean_fold_calibration(d,465201L)
  changed<-d;changed$U<-rev(changed$U);changed$Y<-rnorm(nrow(d));changed$cell<-"unused"
  for(n in intersect(c("u_raw_value","u_raw_std","u_log_std","uq"),names(changed)))changed[[n]]<-NULL
  b<-ocean_fold_calibration(changed,465201L)
  testthat::expect_identical(a$calib_o20,b$calib_o20)
  testthat::expect_identical(a$calib_o50,b$calib_o50)
  shuffled<-d[rev(seq_len(nrow(d))),]
  z<-ocean_fold_calibration(shuffled,465201L);back<-match(d$cohort_id,z$cohort_id)
  testthat::expect_identical(a$calib_o20,z$calib_o20[back])
  testthat::expect_identical(a$calib_o50,z$calib_o50[back])
  audit<-attr(a,"calibration_design")$audit
  for(i in seq_len(nrow(audit))) {
    r<-audit[i,];rows<-which(a$C==r$C);k<-r$n_selected
    x<-cbind(salinity=a$G2salinity[rows],log_pressure=log(a$G2pressure[rows]))
    bins<-ocean_calibration_bins(x,k)[,r$axis==c("salinity","log_pressure")]
    occupied<-tabulate(bins[a[[paste0("calib_",tolower(r$calibration))]][rows]==1L],nbins=k)
    testthat::expect_equal(sum(occupied>0L),r$occupied_bins)
    testthat::expect_equal(max(occupied),r$max_bin_count)
    testthat::expect_identical(all(occupied==1L),r$exact_latin_margin)
  }
})

testthat::test_that("ties disclose approximate coverage instead of fabricating Latin strata", {
  d<-data.frame(cohort_id=seq_len(80),C=rep(1:4,each=20),G2salinity=35,G2pressure=700)
  a<-ocean_fold_calibration(d,123L);audit<-attr(a,"calibration_design")$audit
  testthat::expect_false(any(audit$exact_latin_margin))
  testthat::expect_true(all(audit$occupied_bins==1L))
  testthat::expect_true(all(audit$distinct_X_values==1L))
  testthat::expect_equal(sum(a$calib_o20),16L)
  testthat::expect_equal(sum(a$calib_o50),40L)
  testthat::expect_true(all(a$calib_o20<=a$calib_o50))
})

testthat::test_that("calibration-only preparation preserves frozen inputs and prevents overwrite", {
  out<-tempfile("ocean-X-calibration-");on.exit(unlink(out,recursive=TRUE),add=TRUE)
  source<-file.path(repo,"data/ocean_silicate/n500")
  before<-tools::md5sum(list.files(source,full.names=TRUE))
  testthat::expect_message(prepare_ocean_folds("ocean500",out,repo,calibration_only=TRUE),"Strict LHS is not guaranteed")
  testthat::expect_identical(before,tools::md5sum(names(before)))
  unchanged<-c("cohort_manifest.csv",sprintf("fold%02d_test.csv",1:5),sprintf("fold%02d_test_U_masks.csv",1:5))
  testthat::expect_identical(unname(tools::md5sum(file.path(source,unchanged))),unname(tools::md5sum(file.path(out,unchanged))))
  for(f in 1:5) {
    old<-read.csv(file.path(source,sprintf("fold%02d_train.csv",f)))
    new<-read.csv(file.path(out,sprintf("fold%02d_train.csv",f)))
    cols<-setdiff(names(old),c("calib_o20","calib_o50"))
    testthat::expect_equal(old[,cols],new[,cols],tolerance=1e-12)
  }
  testthat::expect_equal(nrow(read.csv(file.path(out,"calibration_X_coverage.csv"))),80L)
  testthat::expect_error(prepare_ocean_folds("ocean500",out,repo,calibration_only=TRUE),"already exists")
  meta<-jsonlite::read_json(file.path(out,"freeze_manifest.json"),simplifyVector=TRUE)
  testthat::expect_false(meta$calibration_uses_U)
  testthat::expect_false(meta$calibration_uses_Y)
  testthat::expect_false(meta$calibration_uses_cell)
  testthat::expect_equal(meta$calibration_design$seed_by_fold,465200L+1:5)
})
