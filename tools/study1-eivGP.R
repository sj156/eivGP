# Study I paper reproduction: eivGP. Default is plan, never fitting.
# Complete paper workload; calibrated comparisons are included automatically.
script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
if(is.na(script))stop("Run with Rscript; notebooks use study1-common.R")
env <- new.env(parent=globalenv())
sys.source(file.path(dirname(normalizePath(script)),"study1-common.R"),envir=env)
env$study1_lane_cli("eivGP",script)
