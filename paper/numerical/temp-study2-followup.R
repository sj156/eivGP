# Temporary manuscript-completion entrypoint. Delete after follow-up work.
# Fresh reproduction remains in study2.R. No inference is duplicated here.
# Default action: plan. Use --action=run --profile=paper --cores=4 to execute.
# Study II additionally requires --full-fit-root or --full-fit-map.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
if(is.na(script)) stop("Run this temporary entrypoint with Rscript")
env <- new.env(parent=globalenv())
sys.source(file.path(dirname(normalizePath(script)), "common.R"), envir=env)
args <- commandArgs(trailingOnly=TRUE)
if(any(grepl("^--component=",args))) stop("This entrypoint selects follow-up mode; omit --component")
env$run_study_cli("study2", c("--component=follow-up", args))
