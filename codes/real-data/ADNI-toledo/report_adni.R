#!/usr/bin/env Rscript
# Rebuild and combine completed validation reports without model fitting.

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
project <- dirname(normalizePath(sub("^--file=", "", file_arg[1L])))
source(file.path(project, "adni_case_study_helpers.R"))
base <- Sys.getenv("ADNI_OUTPUT_DIR", file.path(project, "outputs"))

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

args <- commandArgs(TRUE)
spec <- Sys.getenv("ADNI_VALIDATIONS", Sys.getenv("ADNI_VALIDATION_ID", "1-20"))
if (length(args)) {
  if (length(args) == 2L && args[1L] == "--validations") spec <- args[2L]
  else if (length(args) == 1L && startsWith(args[1L], "--validations=")) {
    spec <- sub("^--validations=", "", args[1L])
  } else stop("Usage: Rscript report_adni.R [--validations 1-20]")
}
ids <- parse_validation_ids(spec)

tables <- lapply(ids, function(id) {
  root <- file.path(base, sprintf("validation_%02d", id))
  if (!dir.exists(root)) stop("Collect every requested output first. Missing: ", root)
  table <- adni_case_report(root, validation_id = id, smoke = FALSE)
  if (!"method" %in% names(table) || !all(table$complete_comparison)) {
    stop("Validation ", id, " is incomplete; no combined table produced.")
  }
  table$validation_id <- id
  table
})
x <- do.call(rbind, tables)
combined_dir <- file.path(base, "combined")
dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
adni_write_csv(x, file.path(combined_dir, "all_validation_main_metrics.csv"))

metric_names <- c("RMSE", "CRPS", "Coverage95", "Width95", "IntervalScore95")
summary <- do.call(rbind, lapply(split(x, interaction(x$method, x$group)), function(z) {
  out <- data.frame(
    method = z$method[1L], group = z$group[1L], validations = nrow(z),
    mean_test_n = mean(z$n), stringsAsFactors = FALSE
  )
  for (metric in metric_names) {
    out[[paste0(metric, "_mean")]] <- mean(z[[metric]])
    out[[paste0(metric, "_SD")]] <- stats::sd(z[[metric]])
  }
  out$any_diagnostic_flag <- any(z$diagnostic_flag)
  out
}))
adni_write_csv(summary, file.path(combined_dir, "table1_prediction_20validations.csv"))
writeLines(c(
  "# ADNI prediction comparison across validations", "",
  "Values are the arithmetic mean and descriptive SD of each validation-level metric.", "",
  as.character(knitr::kable(summary, format = "pipe", digits = 3))
), file.path(combined_dir, "table1_prediction_20validations.md"))
cat("Wrote combined/table1_prediction_20validations.csv and .md\n")
