# Study II paper reproduction: competitors. Default is plan, never fitting.
# Complete paper workload; calibrated comparisons are included automatically.
# Run/resume: Rscript tools/study2-competitors.R --action=run --cores=12
# Preserve a historical source: use study2-resume-competitors.R instead.
# Setup and output/data-root options: paper/study2-competitors-README.md
script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
if(is.na(script))stop("Run with Rscript; notebooks use study2-common.R")
env <- new.env(parent=globalenv())
sys.source(file.path(dirname(normalizePath(script)),"study2-common.R"),envir=env)
env$study2_lane_cli("competitors",script)
