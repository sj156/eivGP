#!/usr/bin/env Rscript
# Execute the R chunks in the selected notebook without RStudio or Pandoc.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args %in% c("adni", "ocean", "data"))
  stop("Usage: Rscript codes/real-data/run_application.R adni|ocean|data")
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = TRUE))
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
if (!requireNamespace("knitr", quietly = TRUE)) stop("Run codes/real-data/setup.R first.")
run <- function(application) {
  env <- new.env(parent = globalenv())
  code <- tempfile(fileext = ".R")
  on.exit({
    unlink(code)
    if (exists("LOCK_OWNED", env, inherits = FALSE) && isTRUE(env$LOCK_OWNED))
      unlink(env$RUN_LOCK, recursive = TRUE)
  }, add = TRUE)
  if (application == "adni") {
    project <- file.path(root, "ADNI-toledo")
    Sys.setenv(ADNI_PROJECT_DIR = project)
    notebook <- file.path(project, "ADNI_Toledo_EIVGP.Rmd")
  } else if (application == "data") {
    env$DATA_WORKFLOW_DIR <- root
    notebook <- file.path(root, "Raw_Data_Processing.Rmd")
  } else {
    env$OCEAN_REALDATA_DIR <- root
    notebook <- file.path(root, "Ocean_Validation.Rmd")
  }
  invisible(knitr::purl(notebook, output = code, quiet = TRUE))
  sys.source(code, envir = env)
}
run(args[1L])
