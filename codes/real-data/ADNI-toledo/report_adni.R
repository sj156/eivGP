#!/usr/bin/env Rscript
# Rebuild reports without any model fitting. Run after collecting the selected repeats.
f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
project <- dirname(normalizePath(sub("^--file=", "", f[1])))
source(file.path(project, "adni_case_study_helpers.R"))
base <- Sys.getenv("ADNI_OUTPUT_DIR", file.path(project, "outputs"))
args <- commandArgs(TRUE)
repeat_spec <- Sys.getenv("ADNI_REPEATS", "2,3")
if (length(args)) {
  if (length(args) == 2L && args[1L] == "--repeats") repeat_spec <- args[2L]
  else if (length(args) == 1L && startsWith(args[1L], "--repeats=")) repeat_spec <- sub("^--repeats=", "", args[1L])
  else stop("Usage: Rscript report_adni.R [--repeats 1,2,3]")
}
if (!grepl("^[123](,[123])*$", repeat_spec)) stop("Select repeats from 1,2,3.")
repeats <- as.integer(strsplit(repeat_spec, ",", fixed = TRUE)[[1]])
if (anyDuplicated(repeats)) stop("Duplicate repeats are not allowed.")
tables <- list()
for (r in repeats) {
  root <- file.path(base, paste0("repeat", r))
  if (!dir.exists(root)) stop("Collect all requested repeats first. Missing: ", root)
  tb <- adni_case_report(root, 1:3, FALSE)
  if (!"method" %in% names(tb)) stop("Repeat ", r, " has no complete common-fold comparison. See its appendix.")
  if (!all(tb$complete_comparison)) stop("Repeat ", r, " is incomplete; no combined table produced.")
  tb$repeat_id <- r; tables[[as.character(r)]] <- tb
}
x <- do.call(rbind, tables)
dir.create(file.path(base, "combined"), showWarnings = FALSE)
adni_write_csv(x, file.path(base, "combined", "all_repeat_main_metrics.csv"))
summary <- do.call(rbind, lapply(split(x, interaction(x$method, x$group)), function(z) {
  data.frame(method = z$method[1], group = z$group[1], unique_participants = z$n[1],
    repeats = nrow(z), RMSE = sqrt(mean(z$RMSE^2)), CRPS = mean(z$CRPS),
    Coverage95 = mean(z$Coverage95), Width95 = mean(z$Width95),
    IntervalScore95 = mean(z$IntervalScore95), any_diagnostic_flag = any(z$diagnostic_flag))
}))
adni_write_csv(summary, file.path(base, "combined", "table1_prediction.csv"))
writeLines(c("# Complete repeated-CV prediction comparison", "",
  "Equal-repeat averages of participant losses; RMSE is sqrt(mean squared error).",
  "The same participants recur in each repeat. Repeats are not independent studies; no replication SE or significance claim is implied.",
  "Inspect repeat-specific tables, diagnostics, and paired loss differences before interpreting apparent advantages.", "",
  as.character(knitr::kable(summary, format = "pipe", digits = 3))),
  file.path(base, "combined", "table1_prediction.md"))
cat("Wrote combined/table1_prediction.csv and .md\n")
