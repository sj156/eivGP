## Portable launcher for a publication-mode numerical study.
##
## Usage:
##   EIVGP_ARTIFACT_ROOT=/path/to/reproduction \
##     Rscript codes/cli/run_publication_study.R study1
##
## Generate frozen Study I data only:
##   EIVGP_ARTIFACT_ROOT=/path/to/reproduction \
##     Rscript codes/cli/run_publication_study.R study1-data
##
## Optional environment variables:
##   EIVGP_RUN_MODE       dry_run, data, smoke, or publication (default)
##   EIVGP_WORKERS        number of replication workers (default: 8)
##   EIVGP_ARTIFACT_ROOT  root directory for data, results, and HTML report
##   MIXEDGP_R_LIBRARY    optional R package library used for dependencies

arguments <- commandArgs(trailingOnly = TRUE)
study_target <- if (length(arguments) == 0L) "study1" else arguments[[1L]]
study_target <- match.arg(
  study_target, c("study1-data", "study1", "study2-data", "study2")
)
data_only <- grepl("-data$", study_target)
study <- sub("-data$", "", study_target)

file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_argument) != 1L) {
  stop("Run this launcher with Rscript.")
}
script_dir <- dirname(normalizePath(
  sub("^--file=", "", file_argument[[1L]]),
  winslash = "/", mustWork = TRUE
))
repository_dir <- normalizePath(
  file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE
)

dependency_library <- Sys.getenv("MIXEDGP_R_LIBRARY", unset = "")
if (nzchar(dependency_library)) {
  dependency_library <- normalizePath(
    dependency_library, winslash = "/", mustWork = TRUE
  )
  .libPaths(unique(c(dependency_library, .libPaths())))
}

run_mode <- Sys.getenv("EIVGP_RUN_MODE", unset = "publication")
run_mode <- match.arg(
  run_mode,
  c("dry_run", "smoke", "publication")
)
workers <- suppressWarnings(as.integer(Sys.getenv("EIVGP_WORKERS", unset = "8")))
if (length(workers) != 1L || is.na(workers) || workers < 1L) {
  stop("EIVGP_WORKERS must be one positive integer.")
}
artifact_root <- Sys.getenv(
  "EIVGP_ARTIFACT_ROOT", unset = file.path(repository_dir, "reproduction")
)
artifact_root <- normalizePath(
  artifact_root, winslash = "/", mustWork = FALSE
)

## R CMD INSTALL validates Imports but intentionally does not download them.
## Check all package-load and report-rendering requirements before opening the
## numerical document, so users get one actionable setup message.
required_runtime <- c(
  "posterior", "TruncatedNormal", "rmarkdown", "knitr",
  "ggplot2", "dplyr", "tidyr", "patchwork"
)
missing_runtime <- required_runtime[
  !vapply(required_runtime, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_runtime)) {
  stop(
    "Required eivGP reproduction package(s) are missing: ",
    paste(missing_runtime, collapse = ", "), ". Run:\n  Rscript --vanilla ",
    file.path(repository_dir, "codes", "cli", "install_eivgp_dependencies.R"),
    "\nthen reinstall eivGP with:\n  R CMD INSTALL eivGP",
    call. = FALSE
  )
}
study_document <- file.path(repository_dir, "replication", "01_numerical_experiments.Rmd")
report_dir <- file.path(artifact_root, "reports")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
profile <- if (run_mode == "smoke") "quick" else "paper"
action <- if (run_mode == "dry_run") "plan" else if (data_only) "data" else "run"
core_budget <- suppressWarnings(as.numeric(Sys.getenv("MIXEDGP_CORE_BUDGET", as.character(workers))))
rmarkdown::render(study_document, output_dir = report_dir,
  output_file = paste0(study, "_", profile, "_", action, ".html"),
  params = list(profile = profile, study = study, action = action, cores = core_budget,
                artifact_root = artifact_root), envir = new.env(parent = globalenv()))
