#!/usr/bin/env Rscript

## Class-adjusted linear baseline on the frozen 100/300 split:
##   y ~ temperature + salinity + factor(c)
##
## Usage:
##   Rscript 06_run_lm_class.R [DATA_DIR] [OUT_DIR]

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
out_dir <- if (length(args) >= 2L) args[[2L]] else file.path(bundle_dir, "outputs", "lm_class")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, mustWork = TRUE)

dat <- utils::read.csv(file.path(data_dir, "study1_data.csv"), stringsAsFactors = FALSE)
train <- dat[dat$split == "train", , drop = FALSE]
test <- dat[dat$split == "test", , drop = FALSE]
if (nrow(train) != 100L || nrow(test) != 300L) {
  stop("study1_data.csv must be the frozen 100/300 split.")
}
train$c <- factor(train$c, levels = seq_len(6L))
test$c <- factor(test$c, levels = seq_len(6L))

fit <- stats::lm(y ~ x_temperature + x_salinity + c, data = train)
pred <- stats::predict(fit, newdata = test, se.fit = TRUE)
sigma <- stats::sigma(fit)
pred_sd <- sqrt(pred$se.fit^2 + sigma^2)
n_draw <- 600L
set.seed(20260805L)
draws <- matrix(stats::rnorm(n_draw * nrow(test), mean = pred$fit, sd = pred_sd),
                nrow = n_draw, byrow = TRUE)

crps_one <- function(x, y) {
  x <- sort(as.numeric(x)); n <- length(x)
  pair_mean <- 2 * sum((2 * seq_len(n) - n - 1) * x) / n^2
  mean(abs(x - y)) - 0.5 * pair_mean
}
q <- apply(draws, 2L, stats::quantile, probs = c(0.025, 0.975), names = FALSE)
width <- q[2L, ] - q[1L, ]
metrics <- data.frame(
  method = "LM-Class",
  formula = "y ~ x_temperature + x_salinity + factor(c)",
  RMSE = sqrt(mean((pred$fit - test$y)^2)),
  MAE = mean(abs(pred$fit - test$y)),
  CRPS = mean(vapply(seq_len(nrow(test)), function(j) crps_one(draws[, j], test$y[j]), numeric(1))),
  Coverage95 = mean(test$y >= q[1L, ] & test$y <= q[2L, ]),
  Width95 = mean(width),
  IntervalScore95 = mean(width + 40 * pmax(q[1L, ] - test$y, 0) + 40 * pmax(test$y - q[2L, ], 0)),
  stringsAsFactors = FALSE
)
utils::write.csv(metrics, file.path(out_dir, "predictive_metrics.csv"), row.names = FALSE)
saveRDS(fit, file.path(out_dir, "lm_class_fit.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(metrics)
