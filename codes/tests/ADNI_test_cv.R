f <- grep("^--file=", commandArgs(FALSE), value=TRUE)
code <- dirname(dirname(normalizePath(sub("^--file=", "", f[1]))))
source(file.path(code,"ADNI_cv.R"))
source(file.path(code,"ADNI_case_study_helpers.R"))
dat <- data.frame(RID=1:495,R=rep(0:1,c(362,133)),y_centiloid=0,diagnosis="synthetic")
set.seed(43); old <- .Random.seed
assignment <- adni_make_folds(dat)
stopifnot(identical(.Random.seed,old),all(table(assignment$fold)==99L),
 all(table(assignment$fold[assignment$R==1]) %in% 26:27),
 identical(assignment,adni_make_folds(dat[495:1,])), !anyDuplicated(assignment$RID))
# An outcome change cannot change the assignment.
dat$y_centiloid <- seq_len(495)
stopifnot(identical(assignment,adni_make_folds(dat)))
run <- function() {
 root <- tempfile("adni-fivefold-report-");dir.create(root)
 on.exit(unlink(root,recursive=TRUE),add=TRUE)
 adni_write_csv(assignment,file.path(root,"appendix","fold_assignments.csv"))
 dir.create(file.path(root,"main"))
 for(fold in 1:5) {
  test <- dat[dat$RID %in% assignment$RID[assignment$fold==fold],]
  rows <- do.call(rbind,lapply(ADNI_METHODS,function(method) {
   draws <- matrix(rep(test$y_centiloid+fold,each=8),8,nrow(test))
   adni_prediction_rows(draws,test,method,"no_test_CSF",5L,fold,TRUE)
  }))
  status <- data.frame(method=ADNI_METHODS,optimization_status="converged",fold=fold)
  adni_write_csv(rows,file.path(root,paste0("fold_",fold),"case-study","predictions.csv"))
  adni_write_csv(status,file.path(root,paste0("fold_",fold),"case-study","method_status.csv"))
 }
 main <- adni_case_report(root)
 stopifnot(all(main$complete_comparison),all(main$n[main$group=="overall"]==495L),
   all(abs(main$RMSE[main$group=="overall"]-sqrt(mean((1:5)^2)))<1e-12))
 # Wrong fold assignment must fail even when participant counts match.
 file <- file.path(root,"appendix","fold_assignments.csv")
 bad <- assignment;bad$fold <- bad$fold %% 5L+1L;adni_write_csv(bad,file)
 stopifnot(inherits(try(adni_case_report(root),silent=TRUE),"try-error"))
 adni_write_csv(assignment,file)
 # Losing one fold of one comparator leaves an explicitly incomplete table.
 file <- file.path(root,"fold_5/case-study/predictions.csv")
 rows <- read.csv(file);adni_write_csv(rows[rows$method!="EzGP",],file)
 main <- adni_case_report(root)
 stopifnot(all(!main$complete_comparison))
}
run()
cat("PASS: deterministic balanced folds; row-order and outcome independence; pooled 495-person report; incomplete/mismatched data rejected.\n")
