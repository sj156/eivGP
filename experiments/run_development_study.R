args <- commandArgs(trailingOnly = TRUE)
study <- if (length(args)) match.arg(args[1L], c("study1", "study2")) else "study1"
action <- if (length(args) > 1L) match.arg(args[2L], c("plan", "run")) else "plan"
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(file_arg) == 1L)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), ".."), mustWork = TRUE)
core_budget <- suppressWarnings(as.numeric(Sys.getenv("MIXEDGP_CORE_BUDGET", unset = "")))
if (length(core_budget) != 1L || !is.finite(core_budget) || core_budget < 1 || core_budget != floor(core_budget)) {
  stop("Set MIXEDGP_CORE_BUDGET to your CPU budget, for example 12 or 16.")
}
## Set before loading R's numerical dependencies. For BLAS implementations
## initialized at process startup, also set these variables in the shell.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
if (!requireNamespace("rmarkdown", quietly = TRUE)) stop("Install rmarkdown to render the development document.")
report_dir <- file.path(repo, "reproduction", "development", "reports")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
rmarkdown::render(
  file.path(repo, "experiments", "development_numerical_experiment.Rmd"),
  output_file = paste0(study, "_development_", action, ".html"),
  output_dir = report_dir,
  params = list(study = study, execute = action == "run", repository_dir = repo, core_budget = core_budget),
  envir = new.env(parent = globalenv())
)
