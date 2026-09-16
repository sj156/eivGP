#!/usr/bin/env Rscript

## Aggregate completed 10-redraw EIV-GP metrics for one calibration size.
##
## Usage:
##   Rscript 02_aggregate_replicates.R N_CALIB [REPLICATE_OUT_ROOT]
##
## REPLICATE_OUT_ROOT defaults to ./outputs/replicate10 and should contain
## rep01/c10/predictive_metrics.csv, ..., rep10/c50/predictive_metrics.csv.

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
bundle_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE))
} else {
  normalizePath(getwd())
}
args <- commandArgs(trailingOnly = TRUE)
n_calib <- if (length(args) >= 1L) as.integer(args[[1L]]) else 10L
if (!n_calib %in% c(10L, 25L, 50L)) stop("N_CALIB must be 10, 25, or 50.")
root <- if (length(args) >= 2L) {
  normalizePath(args[[2L]], mustWork = TRUE)
} else {
  normalizePath(file.path(bundle_dir, "outputs", "replicate10"), mustWork = TRUE)
}

rows <- lapply(seq_len(10L), function(i) {
  p <- file.path(root, sprintf("rep%02d", i), paste0("c", n_calib), "predictive_metrics.csv")
  if (file.exists(p)) utils::read.csv(p, stringsAsFactors = FALSE) else NULL
})
rows <- Filter(Negate(is.null), rows)
if (!length(rows)) stop("No completed replicates under ", root, " for |O|=", n_calib)
all <- do.call(rbind, rows)
summary <- data.frame(
  n_calib = n_calib,
  n_complete = nrow(all),
  RMSE_mean = mean(all$RMSE), RMSE_sd = stats::sd(all$RMSE),
  MAE_mean = mean(all$MAE), MAE_sd = stats::sd(all$MAE),
  CRPS_mean = mean(all$CRPS), CRPS_sd = stats::sd(all$CRPS),
  Coverage95_mean = mean(all$Coverage95), Coverage95_sd = stats::sd(all$Coverage95),
  Width95_mean = mean(all$Width95), Width95_sd = stats::sd(all$Width95),
  IntervalScore95_mean = mean(all$IntervalScore95),
  IntervalScore95_sd = stats::sd(all$IntervalScore95)
)
agg_dir <- file.path(root, "aggregate", paste0("c", n_calib))
dir.create(agg_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(all, file.path(agg_dir, "prediction_all_reps.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(agg_dir, "prediction_mean_sd.csv"), row.names = FALSE)
print(summary)
