#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
usage <- paste(
  "Usage: Rscript codes/real-data/run_application.R adni|ocean|data",
  "[--cores N] [--validations 1-20] [--plan]"
)
if (!length(args) || !args[1L] %in% c("adni", "ocean", "data")) stop(usage)
application <- args[1L]
rest <- args[-1L]
plan_only <- FALSE
cores <- Sys.getenv("EIVGP_CORES", "")
validation_spec <- Sys.getenv("ADNI_VALIDATIONS", Sys.getenv("ADNI_VALIDATION_ID", "1"))
validation_option_seen <- FALSE
if (application == "adni" && nzchar(Sys.getenv("ADNI_REPEATS", Sys.getenv("ADNI_REPEAT_ID", "")))) {
  stop("The repeat interface has been replaced by --validations; use validation IDs 1--20.")
}

while (length(rest)) {
  if (rest[1L] == "--plan") {
    plan_only <- TRUE
    rest <- rest[-1L]
  } else if (rest[1L] %in% c("--cores", "--validations") && length(rest) >= 2L) {
    if (rest[1L] == "--cores") cores <- rest[2L] else {
      validation_spec <- rest[2L]; validation_option_seen <- TRUE
    }
    rest <- rest[-c(1L, 2L)]
  } else if (startsWith(rest[1L], "--cores=")) {
    cores <- sub("^--cores=", "", rest[1L]); rest <- rest[-1L]
  } else if (startsWith(rest[1L], "--validations=")) {
    validation_spec <- sub("^--validations=", "", rest[1L]); validation_option_seen <- TRUE
    rest <- rest[-1L]
  } else stop(usage)
}
if (application != "adni" && validation_option_seen) {
  stop("--validations applies only to the ADNI application.")
}

parse_count <- function(raw, name) {
  if (!grepl("^[1-9][0-9]*$", raw) || nchar(raw) > 6L) stop("Invalid ", name)
  as.integer(raw)
}
parse_validation_ids <- function(raw) {
  if (!nzchar(raw)) stop("--validations cannot be empty.")
  pieces <- strsplit(raw, ",", fixed = TRUE)[[1L]]
  ids <- unlist(lapply(pieces, function(piece) {
    if (grepl("^[1-9][0-9]*$", piece)) return(as.integer(piece))
    if (!grepl("^[1-9][0-9]*-[1-9][0-9]*$", piece)) {
      stop("Invalid --validations selection: ", raw)
    }
    endpoints <- as.integer(strsplit(piece, "-", fixed = TRUE)[[1L]])
    if (endpoints[1L] > endpoints[2L]) stop("Descending validation ranges are not allowed.")
    seq.int(endpoints[1L], endpoints[2L])
  }), use.names = FALSE)
  if (!length(ids) || anyNA(ids) || any(ids < 1L | ids > 20L)) {
    stop("--validations must select IDs from 1 to 20.")
  }
  if (anyDuplicated(ids)) stop("Duplicate validation IDs are not allowed.")
  as.integer(ids)
}

if (nzchar(cores)) {
  budget <- parse_count(cores, "--cores")
  Sys.setenv(EIVGP_CORES = budget, MIXEDGP_CORES = budget)
} else budget <- parse_count(Sys.getenv("EIVGP_N_CORES", "4"), "EIVGP_N_CORES")
validation_ids <- if (application == "adni") parse_validation_ids(validation_spec) else 1L
if (application == "adni") {
  Sys.setenv(ADNI_VALIDATIONS = paste(validation_ids, collapse = ","))
  Sys.unsetenv(c("ADNI_REPEAT_ID", "ADNI_REPEATS"))
}

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
plan <- adni_core_plan(
  budget, length(validation_ids), 4L,
  parse_count(Sys.getenv("ADNI_MAX_FOLD_WORKERS", "9"), "ADNI_MAX_FOLD_WORKERS")
)
if (application == "adni") {
  adni_print_core_plan(plan)
  cat("Validations:", paste(validation_ids, collapse = ","), "| one fit unit per validation\n")
}
if (plan_only) {
  if (application != "adni") stop("--plan is available for ADNI.")
  cat("Plan only: no data loaded and no fitting started.\n")
  quit(status = 0L)
}
if (!requireNamespace("knitr", quietly = TRUE)) stop("Run codes/real-data/setup.R first.")

run <- function() {
  environments <- list()
  code <- tempfile(fileext = ".R")
  on.exit({
    unlink(code)
    for (env in environments) {
      if (exists("LOCK_OWNED", env, inherits = FALSE) && isTRUE(env$LOCK_OWNED)) {
        unlink(env$RUN_LOCK, recursive = TRUE)
      }
    }
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
  for (validation_id in validation_ids) {
    key <- as.character(validation_id)
    env <- new.env(parent = globalenv())
    environments[[key]] <- env
    env$ADNI_PREPARE_ONLY <- TRUE
    env$ADNI_RUN_VALIDATION_ID <- validation_id
    env$ADNI_RUN_CORE_PLAN <- plan
    sys.source(code, envir = env)
  }

  jobs <- do.call(rbind, lapply(validation_ids, function(validation_id) data.frame(
    validation_id = validation_id, fold = 1L,
    output_root = environments[[as.character(validation_id)]]$OUTPUT_ROOT
  )))
  status_path <- file.path(dirname(jobs$output_root[1L]), paste0(
    if (smoke) "smoke_" else "", "scheduler-validation-",
    paste(validation_ids, collapse = "-"), "-status.csv"
  ))
  cat("Starting", nrow(jobs), "validation jobs. Follow validation_XX/fold_1/worker.log or CURRENT_STATUS.txt.\n")
  values <- adni_run_jobs(jobs, function(validation_id, fold) {
    environments[[as.character(validation_id)]]$run_one_fold(fold)
  }, plan, status_path)

  for (validation_id in validation_ids) {
    env <- environments[[as.character(validation_id)]]
    combined <- do.call(rbind, values[jobs$validation_id == validation_id])
    env$adni_write_csv(combined, file.path(env$OUTPUT_ROOT, "validation_metrics.csv"))
    env$adni_case_report(env$OUTPUT_ROOT, validation_id = validation_id, smoke = smoke)
    writeLines(capture.output(sessionInfo()), file.path(env$OUTPUT_ROOT, "sessionInfo.txt"))
  }
  cat("All requested validation jobs finished. Inspect method status and convergence before combined reporting.\n")
}

run()
