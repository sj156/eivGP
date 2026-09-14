#!/usr/bin/env Rscript

## Ocean EzGP audit with an external covariance-consistency wrapper.
##
## The author package is never edited.  EzGP::EzGP_fit() is called exactly as
## before, including the training-only nugget selection used by the existing
## eivGP adapter.  The strict EzGP::EzGP_predict() result is retained as an
## audit.  For the repaired comparison only, the external predictor uses the
## same fitted covariance parameters but treats a new response at a repeated
## input as a new noisy observation: the nugget is kept on the training
## covariance diagonal and added once to the test predictive variance.
##
## Usage:
##   Rscript 05_run_ezgp_wrapper.R [DATA_DIR] [OUT_DIR]

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
out_dir <- if (length(args) >= 2L) args[[2L]] else file.path(bundle_dir, "outputs", "ezgp_wrapper")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, mustWork = TRUE)

if (!requireNamespace("eivGP", quietly = TRUE)) {
  stop("Install eivGP >= 0.3.1 before running this script.")
}
if (!requireNamespace("EzGP", quietly = TRUE)) stop("EzGP is not installed.")
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
train_x <- eivGP:::make_mixedgp_numeric_matrix(X, C, m_vec = 6L)
test_x <- eivGP:::make_mixedgp_numeric_matrix(X_new, C_new, m_vec = 6L)

status_path <- file.path(out_dir, "status.txt")
writeLines(c(
  "state: RUNNING",
  paste("updated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  "current: selecting EzGP nugget and fitting author package"
), status_path)

seed <- 2026090803L
p <- ncol(X)
q <- ncol(C)
m_vec <- 6L
tau_fractions <- c(1e-6, 0.0025, 0.01, 0.04, 0.16)
cv_folds <- 3L
maxeval <- 100L
n_draw <- 600L

fit_path <- file.path(out_dir, "ezgp_fit_author.rds")
fit_summary_path <- file.path(out_dir, "fit_summary.csv")
tau_grid_path <- file.path(out_dir, "tau_cv_grid.csv")
if (file.exists(fit_path) && file.exists(fit_summary_path) && file.exists(tau_grid_path)) {
  fit <- readRDS(fit_path)
  fit_summary_existing <- utils::read.csv(fit_summary_path, stringsAsFactors = FALSE)
  tau_grid_existing <- utils::read.csv(tau_grid_path, stringsAsFactors = FALSE)
  tau_selection <- list(
    tau = as.numeric(fit_summary_existing$fitted_tau[[1L]]),
    tau_grid = tau_grid_existing$tau,
    cv_loss = tau_grid_existing$cv_loss,
    cv_score = as.character(fit_summary_existing$cv_score[[1L]])
  )
  fit_elapsed <- as.numeric(fit_summary_existing$fit_elapsed_seconds[[1L]])
} else {
  started <- proc.time()[["elapsed"]]
  tau_selection <- eivGP:::select_ezgp_tau_cv(
    train_x = train_x,
    y_train = train$y,
    p = p,
    q = q,
    m_vec = m_vec,
    seed = seed,
    tau_fractions = tau_fractions,
    cv_folds = cv_folds,
    maxeval = maxeval,
    cv_score = "nlpd"
  )
  fit <- EzGP::EzGP_fit(
    X = train_x, Y = train$y, p = p, q = q, m = m_vec,
    tau = tau_selection$tau, maxeval = maxeval
  )
  fit_elapsed <- proc.time()[["elapsed"]] - started
}
if (length(fit$param) == 0L || any(!is.finite(fit$param))) {
  stop("EzGP returned non-finite fitted parameters.")
}
saveRDS(fit, file.path(out_dir, "ezgp_fit_author.rds"))
utils::write.csv(
  data.frame(
    package = "EzGP", version = as.character(utils::packageVersion("EzGP")),
    fitted_tau = tau_selection$tau, cv_score = tau_selection$cv_score,
    cv_folds = cv_folds, maxeval = maxeval, fit_elapsed_seconds = fit_elapsed,
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "fit_summary.csv"), row.names = FALSE
)
utils::write.csv(
  data.frame(tau = tau_selection$tau_grid, cv_loss = tau_selection$cv_loss,
             stringsAsFactors = FALSE),
  file.path(out_dir, "tau_cv_grid.csv"), row.names = FALSE
)

strict <- EzGP::EzGP_predict(test_x, fit, MSE_on = 1L)
strict_mse <- as.numeric(strict$MSE)
strict_audit <- data.frame(
  test_source_row = test$source_row,
  strict_MSE = strict_mse,
  strict_negative = strict_mse < 0,
  stringsAsFactors = FALSE
)
utils::write.csv(strict_audit, file.path(out_dir, "strict_author_prediction_audit.csv"),
                 row.names = FALSE)

ezgp_predict_covariance_consistent <- function(X_new, model) {
  if (!inherits(model, "EzGP model")) stop("model must be an EzGP model.")
  X_train <- model$data$X
  Y <- model$data$Y
  p <- model$data$p
  q <- model$data$q
  m <- model$data$m
  n <- nrow(X_train)
  tau <- as.numeric(model$data$tau)
  parv <- model$param
  K <- EzGP:::cov_m(X_train, p, q, m, n, parv, tau)
  Tm <- try(chol(K), silent = TRUE)
  if (inherits(Tm, "try-error")) stop("EzGP training covariance is not Cholesky-factorable.")
  invT <- backsolve(Tm, diag(n))
  invc <- invT %*% t(invT)
  m1 <- matrix(rep(1, n), ncol = 1L)
  y <- as.matrix(Y)
  mu <- as.numeric((1 / sum(invc)) * (t(m1) %*% invc %*% y))
  npar <- 1L + q + p + p * sum(m)
  psum <- function(x1, x2, par2) sum(-par2 * (x1 - x2)^2)
  covx <- function(w1, w2) {
    par1 <- parv[seq_len(q + 1L)]
    par2 <- parv[(q + 2L):(q + 1L + p)]
    par3 <- parv[(q + 2L + p):npar]
    x1 <- w1[seq_len(p)]
    z1 <- w1[(p + 1L):(p + q)]
    x2 <- w2[seq_len(p)]
    z2 <- w2[(p + 1L):(p + q)]
    out <- par1[1L] * exp(psum(x1, x2, par2))
    for (i in seq_len(q)) {
      if (z1[i] == z2[i]) {
        l <- z1[i]
        idx <- (sum(m[seq_len(i)]) - m[i] + (l - 1L)) * p + seq_len(p)
        out <- out + par1[i + 1L] * exp(psum(x1, x2, par3[idx]))
      }
    }
    out
  }
  one <- function(wn) {
    gamma <- numeric(n)
    for (i in seq_len(n)) {
      gamma[i] <- covx(wn, X_train[i, ])
    }
    gamma <- as.matrix(gamma)
    yhat <- mu + t(gamma) %*% invc %*% (y - mu * m1)
    latent_mse <- sum(parv[seq_len(q + 1L)]) - t(gamma) %*% invc %*% gamma
    list(
      Y_hat = as.numeric(yhat),
      latent_MSE = as.numeric(latent_mse),
      predictive_MSE = as.numeric(latent_mse) + tau
    )
  }
  vals <- lapply(seq_len(nrow(X_new)), function(i) one(X_new[i, ]))
  list(
    Y_hat = vapply(vals, `[[`, numeric(1), "Y_hat"),
    latent_MSE = vapply(vals, `[[`, numeric(1), "latent_MSE"),
    predictive_MSE = vapply(vals, `[[`, numeric(1), "predictive_MSE")
  )
}

corrected <- ezgp_predict_covariance_consistent(test_x, fit)
if (any(!is.finite(corrected$predictive_MSE)) || any(corrected$predictive_MSE < 0)) {
  stop("The covariance-consistency wrapper returned invalid predictive variances.")
}
draws <- eivGP:::sample_independent_predictive_marginals(
  corrected$Y_hat, corrected$predictive_MSE, n_draw, seed + 10000L
)
saveRDS(draws, file.path(out_dir, "predictive_draws_covariance_wrapper.rds"))
utils::write.csv(
  data.frame(
    test_source_row = test$source_row,
    corrected_latent_MSE = corrected$latent_MSE,
    corrected_predictive_MSE = corrected$predictive_MSE,
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "covariance_wrapper_prediction_components.csv"),
  row.names = FALSE
)

crps_one <- function(x, y) {
  x <- sort(as.numeric(x)); n <- length(x)
  pair_mean <- 2 * sum((2 * seq_len(n) - n - 1) * x) / n^2
  mean(abs(x - y)) - 0.5 * pair_mean
}
q <- apply(draws, 2L, stats::quantile, probs = c(0.025, 0.975), names = FALSE)
width <- q[2L, ] - q[1L, ]
metrics <- data.frame(
  method = "EzGP-covariance-wrapper",
  target = "Y_new | X_new, C_new",
  RMSE = sqrt(mean((colMeans(draws) - test$y)^2)),
  MAE = mean(abs(colMeans(draws) - test$y)),
  Coverage95 = mean(test$y >= q[1L, ] & test$y <= q[2L, ]),
  Width95 = mean(width),
  CRPS = mean(vapply(seq_len(ncol(draws)), function(j) crps_one(draws[, j], test$y[j]), numeric(1))),
  IntervalScore95 = mean(width + 40 * pmax(q[1L, ] - test$y, 0) +
                           40 * pmax(test$y - q[2L, ], 0)),
  stringsAsFactors = FALSE
)
utils::write.csv(metrics, file.path(out_dir, "predictive_metrics.csv"), row.names = FALSE)
utils::write.csv(
  data.frame(
    method = c("EzGP-author-strict", "EzGP-covariance-wrapper"),
    status = c(
      if (any(strict_mse < 0)) "failed_negative_variance" else "success",
      "success"
    ),
    min_reported_variance = c(min(strict_mse), min(corrected$predictive_MSE)),
    n_negative_variance = c(sum(strict_mse < 0), sum(corrected$predictive_MSE < 0)),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "fit_status.csv"), row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
writeLines(c(
  "state: DONE",
  paste("updated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  paste("strict_negative_variances:", sum(strict_mse < 0)),
  paste("corrected_negative_variances:", sum(corrected$predictive_MSE < 0)),
  paste("output_dir:", out_dir)
), status_path)
print(tau_selection)
print(metrics)
print(utils::read.csv(file.path(out_dir, "fit_status.csv")))
