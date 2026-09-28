#!/usr/bin/env Rscript
# Rebuild the pooled five-fold report without fitting.
if (length(commandArgs(TRUE))) stop("Usage: Rscript codes/applications/ADNI_report.R (no repeat/validation arguments)")
f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (!file.exists(f)) f <- gsub("~+~", " ", f, fixed = TRUE)
project <- dirname(normalizePath(f, mustWork = TRUE))
source(file.path(project, "real_data_paths.R"))
source(file.path(project, "ADNI_case_study_helpers.R"))
root <- file.path(application_paths(project)$adni_output, "five_fold")
if (dir.exists(file.path(root, ".run-lock"))) stop("Analysis is locked/running. Wait for it to finish before rebuilding reports.")
table <- adni_case_report(root, 1:5, FALSE)
cat("Wrote ", file.path(root, "main", "table1_prediction.csv"), "\n", sep = "")
if (!"method" %in% names(table) || !all(table$complete_comparison))
  stop("Incomplete method comparison. Inspect appendix/completeness.csv and method_status.csv.")
