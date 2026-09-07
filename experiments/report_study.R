## Usage: Rscript --vanilla experiments/report_study.R RUN_DIRECTORY [report|summarize|plot]
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Supply the saved run directory; no run is selected automatically.")
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), ".."), mustWork = TRUE)
source(file.path(repo, "codes", "00_diagnostics.R"))
source(file.path(repo, "codes", "simulation_helpers.R"))
source(file.path(repo, "codes", "experiment_reporting.R"))
status <- mixedgp_report_run(args[1L], if (length(args) > 1L) args[2L] else "report",
                            code_dir = file.path(repo, "codes"))
print(status)
if (any(status$status != "success")) quit(status = 1L)
