#!/usr/bin/env Rscript
# Temporary generated inputs only: no Dropbox, archived results or manuscript.
script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(normalizePath(script)),"supplement.R"))
numerical_package_check()
root<-tempfile("eivgp-portable-");dir.create(root)
cfg<-settings(list(profile="quick",reps="1",calib="20",methods="FULL","data-root"=file.path(root,"data"),output=file.path(root,"output")))
run_supplement(cfg,"data",file.path(NUMERICAL_REPO_ROOT,"paper/numerical/supplement.R"))
p<-file.path(cfg$data_root,"primary_q4_calibration/study2_primary_q4_n100_ntest100_m4_rep001.rds")
a<-readRDS(p);b<-F("make_study2_synthetic_dataset")(rep_id=1L,scenario="primary",n=100L,n_test=100L,q=4L,m=4L,calib_grid=c(0L,20L,50L,80L),data_seed_base=1000000L)
stopifnot(identical(a,b),identical(a$seeds$data,1010001L),all(a$calibration_sets[["20"]]%in%a$calibration_sets[["50"]]))
before<-file.info(p)$mtime
run_supplement(cfg,"data",file.path(NUMERICAL_REPO_ROOT,"paper/numerical/supplement.R"))
stopifnot(identical(before,file.info(p)$mtime))
# Positive and negative controls on configuration validation, without fitting.
bad<-a;bad$design$n<-99L;saveRDS(bad,p)
e<-tryCatch(run_supplement(cfg,"data"),error=identity)
stopifnot(inherits(e,"error"),grepl("design mismatch",conditionMessage(e)))
unlink(root,recursive=TRUE)
cat("PASS: checkout-only generation, fixed seeds, nested calibration, reuse and design validation.\n")
