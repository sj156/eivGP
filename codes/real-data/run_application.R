#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
usage <- "Usage: Rscript codes/real-data/run_application.R adni|ocean|data [--cores N] [--repeats 2,3] [--plan]"
if (!length(args) || !args[1L] %in% c("adni", "ocean", "data")) stop(usage)
application <- args[1L]; rest <- args[-1L]; plan_only <- FALSE
cores <- Sys.getenv("EIVGP_CORES", "")
repeats_text <- Sys.getenv("ADNI_REPEATS", Sys.getenv("ADNI_REPEAT_ID", "2"))
while (length(rest)) {
  if (rest[1L] == "--plan") { plan_only <- TRUE; rest <- rest[-1L] }
  else if (rest[1L] %in% c("--cores", "--repeats") && length(rest) >= 2L) {
    if (rest[1L] == "--cores") cores <- rest[2L] else repeats_text <- rest[2L]
    rest <- rest[-c(1L, 2L)]
  } else if (startsWith(rest[1L], "--cores=")) {
    cores <- sub("^--cores=", "", rest[1L]); rest <- rest[-1L]
  } else if (startsWith(rest[1L], "--repeats=")) {
    repeats_text <- sub("^--repeats=", "", rest[1L]); rest <- rest[-1L]
  } else stop(usage)
}
parse_count <- function(raw, name) {
  if (!grepl("^[1-9][0-9]*$", raw) || nchar(raw) > 6L) stop("Invalid ", name)
  as.integer(raw)
}
if (nzchar(cores)) {
  budget <- parse_count(cores, "--cores")
  Sys.setenv(EIVGP_CORES = budget, MIXEDGP_CORES = budget)
} else budget <- parse_count(Sys.getenv("EIVGP_N_CORES", "4"), "EIVGP_N_CORES")
if (!grepl("^[23](,[23])*$", repeats_text)) stop("--repeats must be 2, 3, or 2,3.")
repeats <- as.integer(strsplit(repeats_text, ",", fixed = TRUE)[[1]])
if (anyDuplicated(repeats)) stop("Duplicate repeats are not allowed.")
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- sub("^--file=", "", file_arg[1L])
if (!file.exists(script_path)) script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
root <- dirname(normalizePath(script_path, mustWork = TRUE))
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
source(file.path(root, "ADNI-toledo", "adni_parallel.R"))
flag <- function(name) {
  x <- tolower(Sys.getenv(name, "0"))
  if (!x %in% c("0", "1", "true", "false", "yes", "no")) stop("Invalid ", name)
  x %in% c("1", "true", "yes")
}
smoke <- flag("EIVGP_SMOKE_TEST")
folds <- if (smoke && !flag("ADNI_SMOKE_ALL_FOLDS")) 1L else 1:3
plan <- adni_core_plan(budget, length(folds) * length(repeats), 4L,
  parse_count(Sys.getenv("ADNI_MAX_FOLD_WORKERS", "6"), "ADNI_MAX_FOLD_WORKERS"))
if (application == "adni") {
  adni_print_core_plan(plan)
  cat("Repeats:", paste(repeats, collapse = ","), "| folds per repeat:", paste(folds, collapse = ","), "\n")
}
if (plan_only) {
  if (application != "adni") stop("--plan is available for ADNI.")
  cat("Plan only: no data loaded and no fitting started.\n")
  quit(status = 0L)
}
if (!requireNamespace("knitr", quietly = TRUE)) stop("Run codes/real-data/setup.R first.")
run <- function() {
  environments <- list(); code <- tempfile(fileext = ".R")
  on.exit({
    unlink(code)
    for (env in environments) if (exists("LOCK_OWNED", env, inherits = FALSE) && isTRUE(env$LOCK_OWNED))
      unlink(env$RUN_LOCK, recursive = TRUE)
  }, add = TRUE)
  if (application %in% c("adni", "data")) {
    source(file.path(root, "adni_paths.R"), local = TRUE)
    Sys.setenv(ADNI_DATA_DIR = resolve_adni_data_dir(file.path(root, "ADNI-toledo")))
  }
  if (application != "adni") {
    env <- new.env(parent = globalenv())
    env$DATA_WORKFLOW_DIR <- root; env$OCEAN_REALDATA_DIR <- root
    notebook <- file.path(root, if (application == "data") "Raw_Data_Processing.Rmd" else "Ocean_Validation.Rmd")
    invisible(knitr::purl(notebook, output = code, quiet = TRUE))
    sys.source(code, envir = env)
    return(invisible(NULL))
  }
  project <- file.path(root, "ADNI-toledo")
  Sys.setenv(ADNI_PROJECT_DIR = project)
  invisible(knitr::purl(file.path(project, "ADNI_Toledo_EIVGP.Rmd"), output = code, quiet = TRUE))
  # Prepare each repeat in the coordinator: acquire locks, validate inputs, and write shared metadata/EDA once.
  for (r in repeats) {
    env <- new.env(parent = globalenv()); environments[[as.character(r)]] <- env
    env$ADNI_PREPARE_ONLY <- TRUE; env$ADNI_RUN_REPEAT_ID <- r; env$ADNI_RUN_CORE_PLAN <- plan
    sys.source(code, envir = env)
  }
  jobs <- do.call(rbind, lapply(repeats, function(r) data.frame(repeat_id = r, fold = folds,
    output_root = environments[[as.character(r)]]$OUTPUT_ROOT)))
  status_path <- file.path(dirname(jobs$output_root[1]),
    paste0(if (smoke) "smoke_" else "", "scheduler-repeat", paste(repeats, collapse = "-"), "-status.csv"))
  cat("Starting", nrow(jobs), "fold jobs. Follow repeatN/fold_N/worker.log or CURRENT_STATUS.txt.\n")
  values <- adni_run_jobs(jobs, function(repeat_id, fold) {
    environments[[as.character(repeat_id)]]$run_one_fold(fold)
  }, plan, status_path)
  # Combined reporting is coordinator-only, after every requested fold has returned.
  for (r in repeats) {
    env <- environments[[as.character(r)]]
    combined <- do.call(rbind, values[jobs$repeat_id == r])
    env$adni_write_csv(combined, file.path(env$OUTPUT_ROOT, "all_fold_metrics.csv"))
    env$adni_case_report(env$OUTPUT_ROOT, folds, smoke)
    writeLines(capture.output(sessionInfo()), file.path(env$OUTPUT_ROOT, "sessionInfo.txt"))
  }
  cat("All requested fold jobs finished. See repeat main/ and appendix/ reports; inspect method and convergence status.\n")
}
run()
