#!/usr/bin/env Rscript

## UC-GP and LVGP on the frozen 100/300 ocean split, through the eivGP 0.3.1
## published-competitor adapters.  EzGP is also attempted; the author predictor
## fails on this split (see 05_run_ezgp_wrapper.R).
##
## Usage:
##   Rscript 04_run_published_baselines.R [DATA_DIR] [OUT_DIR]

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
bundle_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE))
} else {
  normalizePath(getwd())
}
source(file.path(bundle_dir, "ocean_v031_paths.R"))
args <- commandArgs(trailingOnly = TRUE)
data_dir <- ocean_v031_data_dir(
  bundle_dir,
  explicit = if (length(args) >= 1L) args[[1L]] else NULL
)
out_dir <- if (length(args) >= 2L) args[[2L]] else file.path(bundle_dir, "outputs", "published_baselines")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, mustWork = TRUE)

if (!requireNamespace("eivGP", quietly = TRUE)) {
  stop("Install eivGP >= 0.3.1 before running this script.")
}
library(eivGP)

dat <- utils::read.csv(file.path(data_dir, "study1_data.csv"), stringsAsFactors = FALSE)
train <- dat[dat$split == "train", , drop = FALSE]
test <- dat[dat$split == "test", , drop = FALSE]
if (nrow(train) != 100L || nrow(test) != 300L) {
  stop("study1_data.csv must be the frozen 100/300 split.")
}
X <- as.matrix(train[, c("x_temperature", "x_salinity"), drop = FALSE])
X_new <- as.matrix(test[, c("x_temperature", "x_salinity"), drop = FALSE])
C <- matrix(as.integer(train$c), ncol = 1L, dimnames = list(NULL, "phosphate_class"))
C_new <- matrix(as.integer(test$c), ncol = 1L, dimnames = list(NULL, "phosphate_class"))

status_path <- file.path(out_dir, "status.txt")
writeLines(c(
  "state: RUNNING",
  paste("updated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  "current: fitting UC-GP, LVGP, EzGP"
), status_path)

crps_one <- function(x, y) {
  x <- sort(as.numeric(x)); n <- length(x)
  pair_mean <- 2 * sum((2 * seq_len(n) - n - 1) * x) / n^2
  mean(abs(x - y)) - 0.5 * pair_mean
}
score_draws <- function(draws, y) {
  q <- apply(draws, 2L, stats::quantile, probs = c(0.025, 0.975), names = FALSE)
  width <- q[2L, ] - q[1L, ]
  data.frame(
    RMSE = sqrt(mean((colMeans(draws) - y)^2)),
    MAE = mean(abs(colMeans(draws) - y)),
    Coverage95 = mean(y >= q[1L, ] & y <= q[2L, ]),
    Width95 = mean(width),
    CRPS = mean(vapply(seq_len(ncol(draws)), function(j) crps_one(draws[, j], y[j]), numeric(1))),
    IntervalScore95 = mean(width + 40 * pmax(q[1L, ] - y, 0) + 40 * pmax(y - q[2L, ], 0))
  )
}

methods <- c("UC-GP", "LVGP", "EzGP")
fits <- list(); metric_rows <- list(); status_rows <- list()
status_columns <- c("method", "status", "elapsed_seconds", "implementation",
                    "fitted_noise_variance", "message")
bind_status_rows <- function(rows) {
  do.call(rbind, lapply(rows, function(x) {
    missing <- setdiff(status_columns, names(x))
    for (name in missing) x[[name]] <- NA
    x[, status_columns, drop = FALSE]
  }))
}
for (j in seq_along(methods)) {
  method <- methods[[j]]
  writeLines(c(
    "state: RUNNING",
    paste("updated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste("current:", method)
  ), status_path)
  started <- proc.time()[["elapsed"]]
  ans <- tryCatch(
    fit_mixedgp_competitor(
      method = method, X = X, y = train$y, C = C,
      new_X = X_new, new_C = C_new, m_vec = 6L,
      n_draw = 600L, seed = 2026090800L + j
    ),
    error = function(e) e
  )
  elapsed <- proc.time()[["elapsed"]] - started
  if (inherits(ans, "error")) {
    status_rows[[method]] <- data.frame(
      method = method, status = "failed", elapsed_seconds = elapsed,
      message = conditionMessage(ans), stringsAsFactors = FALSE
    )
    utils::write.csv(bind_status_rows(status_rows),
                     file.path(out_dir, "baseline_fit_status.partial.csv"),
                     row.names = FALSE)
    next
  }
  fits[[method]] <- ans
  metric_rows[[method]] <- cbind(
    data.frame(method = method, target = "Y_new | X_new, C_new", stringsAsFactors = FALSE),
    score_draws(ans$draws, test$y)
  )
  status_rows[[method]] <- data.frame(
    method = method, status = "success", elapsed_seconds = elapsed,
    implementation = ans$implementation,
    fitted_noise_variance = ans$fitted_noise_variance,
    stringsAsFactors = FALSE
  )
  utils::write.csv(do.call(rbind, metric_rows),
                   file.path(out_dir, "baseline_predictive_metrics.partial.csv"),
                   row.names = FALSE)
  utils::write.csv(bind_status_rows(status_rows),
                   file.path(out_dir, "baseline_fit_status.partial.csv"),
                   row.names = FALSE)
  saveRDS(fits, file.path(out_dir, "baseline_fits.partial.rds"))
}
metrics <- if (length(metric_rows)) do.call(rbind, metric_rows) else data.frame()
status <- bind_status_rows(status_rows)
utils::write.csv(metrics, file.path(out_dir, "baseline_predictive_metrics.csv"), row.names = FALSE)
utils::write.csv(status, file.path(out_dir, "baseline_fit_status.csv"), row.names = FALSE)
saveRDS(fits, file.path(out_dir, "baseline_fits.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
writeLines(c(
  "state: DONE",
  paste("updated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  paste("successful methods:", paste(status$method[status$status == "success"], collapse = ", "))
), status_path)
print(status)
print(metrics)
