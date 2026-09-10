f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
project <- dirname(dirname(normalizePath(sub("^--file=", "", f[1]))))
source(file.path(project, "adni_case_study_helpers.R"))
set.seed(91)
a <- matrix(rnorm(28), 7, 4); y <- c(-1, 0, .5, 1)
s <- adni_draw_summary(a, y)
for (j in 1:4) stopifnot(abs(s$CRPS[j] - (mean(abs(a[,j]-y[j])) - mean(abs(outer(a[,j],a[,j],"-")))/2)) < 1e-12)
stopifnot(all(adni_draw_summary(matrix(2, 8, 3), c(1,2,3))$CRPS == c(1,0,1)))
stopifnot(inherits(try(adni_draw_summary(matrix(NA_real_,2,2),1:2),silent=TRUE),"try-error"))
# Pooled RMSE is computed from participant losses, never the mean of fold RMSEs.
z <- data.frame(R = c(0,0,1), squared_error=c(1,1,100), absolute_error=c(1,1,10),
  CRPS=c(1,1,10), NLPD=NA_real_, covered95=TRUE, Width95=1,IntervalScore95=1,covered50=TRUE,covered80=TRUE)
m <- adni_metrics(z)
stopifnot(abs(m$RMSE[m$group=="overall"]-sqrt(34))<1e-12, m$n[m$group=="R0"]==2)
cat("PASS: CRPS against brute-force identity; invalid draws fail; participant-pooled and subgroup scores.\n")

d <- matrix(0, 8, 2); attr(d, "conditional_means") <- matrix(c(1,2), 1, 2)
attr(d, "conditional_vars") <- matrix(c(4,9), 1, 2)
a <- adni_draw_summary(d, c(3,4))
stopifnot(max(abs(a$NLPD + dnorm(c(3,4), c(1,2), c(2,3), log=TRUE))) < 1e-12,
          identical(a$pred_mean, c(1,2)))
cat("PASS: normal-mixture NLPD and conditional means.\n")
root <- tempfile("adni-report-test-"); dir.create(root)
for (d in c("main", "appendix", "fold_1/case-study")) dir.create(file.path(root,d), recursive=TRUE)
test <- data.frame(RID=1:6, R=c(0,0,0,1,1,1), diagnosis="test", y_centiloid=1:6)
rows <- do.call(rbind,lapply(ADNI_METHODS, function(m) adni_prediction_rows(
  matrix(rep(1:6, each=8),8,6),test,m,"no_test_CSF",2L,1L,convergence=TRUE)))
status <- data.frame(method=ADNI_METHODS,status="success",optimization_status="converged",
                     elapsed_seconds=1,message="",warnings="",repeat_id=2,fold=1)
adni_write_csv(rows,file.path(root,"fold_1/case-study/predictions.csv"))
adni_write_csv(status,file.path(root,"fold_1/case-study/method_status.csv"))
main <- adni_case_report(root,1L,TRUE)
stopifnot(nrow(main)==8L,all(main$RMSE==0),all(!main$complete_comparison))
# A missing competitor must not silently disappear from a complete-looking comparison.
adni_write_csv(rows[rows$method!="EzGP",],file.path(root,"fold_1/case-study/predictions.csv"))
status$status[status$method=="EzGP"] <- "failed"
adni_write_csv(status,file.path(root,"fold_1/case-study/method_status.csv"))
main <- adni_case_report(root,1L,TRUE)
stopifnot("status" %in% names(main))
# Same counts but different participants must fail, not yield a matched table.
rows$RID[rows$method=="UC-GP"] <- 101:106
adni_write_csv(rows,file.path(root,"fold_1/case-study/predictions.csv"))
stopifnot(inherits(try(adni_case_report(root,1L,TRUE),silent=TRUE),"try-error"))
unlink(root,recursive=TRUE)
cat("PASS: complete common-fold reporting, explicit missing competitor, mismatched-participant rejection.\n")
