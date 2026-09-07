## Reporting only. Paths may be local, copied from another machine, or synced.
args <- commandArgs(TRUE)
if(length(args)!=3L) stop("Usage: Rscript experiments/combine_results.R MCMC_RUN COMPETITOR_REPORT OUTPUT_DIR")
file_arg <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))
repo <- normalizePath(file.path(dirname(file_arg),".."))
source(file.path(repo,"codes","simulation_helpers.R"))
source(file.path(repo,"codes","combine_experiment_results.R"))
mixedgp_combine_saved_runs(args[1],args[2],args[3])
cat("Combined reporting saved to",args[3],"; no models were fitted.\n")
