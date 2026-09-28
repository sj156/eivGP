#!/usr/bin/env Rscript
# Command-line access to the same functions used by the three reader notebooks.
args <- commandArgs(TRUE)
usage <- "Usage: Rscript codes/applications/run_real_data.R adni|ocean|data [--cores N] [--plan]"
if (!length(args) || !args[1] %in% c("adni", "ocean", "data")) stop(usage)
application <- args[1]; rest <- args[-1]; plan_only <- FALSE
cores <- Sys.getenv("EIVGP_CORES", Sys.getenv("EIVGP_N_CORES", "4"))
while (length(rest)) {
  if (rest[1] == "--plan") { plan_only <- TRUE; rest <- rest[-1] }
  else if (rest[1] == "--cores" && length(rest) >= 2L) {
    cores <- rest[2]; rest <- rest[-c(1,2)]
  } else if (startsWith(rest[1], "--cores=")) {
    cores <- sub("^--cores=", "", rest[1]); rest <- rest[-1]
  } else if (grepl("^--(repeats|validations)", rest[1])) {
    stop("ADNI now uses one fixed five-fold analysis. Omit --repeats/--validations.")
  } else stop(usage)
}
if (!grepl("^[1-9][0-9]*$", cores) || nchar(cores) > 6L) stop("--cores must be a positive integer.")
Sys.setenv(EIVGP_CORES = cores, MIXEDGP_CORES = cores,
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (!file.exists(f)) f <- gsub("~+~", " ", f, fixed = TRUE)
root <- dirname(normalizePath(f, mustWork = TRUE))
source(file.path(root, "real_data_paths.R"))
paths <- application_paths(root)
flag <- function(name) {
  x <- tolower(Sys.getenv(name, "0"))
  if (!x %in% c("0", "1", "true", "false", "yes", "no")) stop("Invalid ", name)
  x %in% c("1", "true", "yes")
}
smoke <- flag("EIVGP_SMOKE_TEST")
folds <- if (smoke && !flag("ADNI_SMOKE_ALL_FOLDS")) 1L else 1:5
source(file.path(root, "ADNI_parallel.R"))
max_workers <- suppressWarnings(as.numeric(Sys.getenv("ADNI_MAX_FOLD_WORKERS", "5")))
plan <- adni_core_plan(as.integer(cores), length(folds), 4L, max_workers)
if (application == "adni") {
  adni_print_core_plan(plan)
  cat("Design: one fixed 5-fold CV | train/test: 396/99 | folds:", paste(folds, collapse = ","), "\n")
  cat("Data:", paths$adni_data, "\nOutputs:", file.path(paths$adni_output, if (smoke) "smoke" else "five_fold"), "\n")
  if (any(nzchar(Sys.getenv(c("ADNI_REPEATS", "ADNI_REPEAT_ID", "ADNI_VALIDATIONS", "ADNI_VALIDATION_ID")))))
    warning("Legacy repeat/validation environment settings are ignored; this is fixed five-fold CV.")
}
if (plan_only) {
  if (application != "adni") stop("--plan is available for ADNI.")
  cat("Plan only: no data loaded and no fitting started.\n")
  quit(status = 0L)
}
source(file.path(paths$repo, "codes", "reporting", "notebook_helpers.R"))
source(file.path(root, "notebook_workflow.R"))
if (application == "adni") {
  result <- adni_run_analysis(paths$repo, if (smoke) "quick" else "paper",
    as.integer(cores), paths$adni_data, paths$adni_output)
  print(result$metrics)
} else if (application == "ocean") {
  source(file.path(root, "ocean_silicate", "workflow.R"))
  profile <- if (tolower(Sys.getenv("OCEAN_QUICK", "true")) %in% c("1", "true")) "quick" else "paper"
  out <- file.path(paths$ocean_output, profile)
  ocean_read_inputs(paths$ocean_data)
  ocean_run_analysis(paths$repo, paths$ocean_data, out, profile)
  source(file.path(root, "ocean_silicate", "reporting.R"))
  print(if (profile == "paper") ocean_read_paper_results(out)$metrics else ocean_read_results(out, profile)$metrics)
} else {
  source(file.path(root, "check_prepared_data.R"))
  check_prepared_data(root)
}
