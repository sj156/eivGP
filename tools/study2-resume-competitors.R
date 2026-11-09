#!/usr/bin/env Rscript
# Defaults to a plan. Run with --action=run --cores=12; rerun to resume.
# Originals are copied and verified once, then remain untouched.
script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
if(is.na(script))stop("Run with Rscript")
env<-new.env(parent=globalenv())
sys.source(file.path(dirname(normalizePath(script)),"study2-common.R"),envir=env)
env$study2_resume_competitors_cli(script)
