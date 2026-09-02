## Portable launcher for a publication-mode numerical study.
##
## Usage:
##   EIVGP_ARTIFACT_ROOT=/path/to/reproduction \
##     Rscript experiments/run_publication_study.R study1
##
## Generate frozen Study I data only:
##   EIVGP_ARTIFACT_ROOT=/path/to/reproduction \
##     Rscript experiments/run_publication_study.R study1-data
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
  file.path(script_dir, ".."), winslash = "/", mustWork = TRUE
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
study_document <- file.path(repository_dir, "experiments", if (data_only) {
  paste0(study, "_synthetic_data.Rmd")
} else {
  paste0(study, "_numerical_experiment.Rmd")
})
report_dir <- file.path(artifact_root, "reports")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Install rmarkdown before running a numerical-study document.")
}
rmarkdown::render(
  input = study_document,
  output_dir = report_dir,
  params = if (data_only) {
    list(
      run_mode = run_mode,
      workers = workers,
      data_root = file.path(
        artifact_root, "data",
        if (identical(run_mode, "publication")) "synthetic" else "smoke",
        study
      )
    )
  } else {
    list(
      run_mode = run_mode,
      workers = workers,
      data_root = file.path(
        artifact_root, "data",
        if (identical(run_mode, "publication")) "synthetic" else "smoke",
        study
      ),
      output_root = file.path(artifact_root, "results", study)
    )
  },
  envir = new.env(parent = globalenv())
)
