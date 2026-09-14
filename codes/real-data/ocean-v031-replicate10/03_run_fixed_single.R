#!/usr/bin/env Rscript

## Run the frozen Ocean single case with eivGP 0.3.1.  The full mode uses the
## same fixed 20,000 / 5,000 schedule planned for the replicate analysis and
## writes continuation checkpoints every 2,000 transitions.
##
## Usage:
##   Rscript 03_run_fixed_single.R \
##     OUT_ROOT [full|smoke] [n_cores] [n_calib]
##   Rscript 03_run_fixed_single.R \
##     DATA_DIR OUT_ROOT [full|smoke] [n_cores] [n_calib]
##
## DATA_DIR defaults to ../ocean_data.

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
bundle_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE))
} else {
  normalizePath(getwd())
}
source(file.path(bundle_dir, "ocean_v031_paths.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L && file.exists(file.path(args[[1L]], "study1_data.csv"))) {
  data_dir <- normalizePath(args[[1L]], mustWork = TRUE)
  args <- args[-1L]
} else {
  data_dir <- ocean_v031_data_dir(bundle_dir)
}
if (!length(args)) {
  stop(paste(
    "Usage: Rscript 03_run_fixed_single.R",
    "[DATA_DIR] OUT_ROOT [full|smoke] [n_cores] [n_calib]"
  ))
}
out_root <- args[[1L]]
mode <- if (length(args) >= 2L) args[[2L]] else "full"
n_cores <- if (length(args) >= 3L) as.integer(args[[3L]]) else 2L
n_calib <- if (length(args) >= 4L) as.integer(args[[4L]]) else 10L
if (!mode %in% c("full", "smoke")) stop("MODE must be full or smoke.")
if (length(n_cores) != 1L || is.na(n_cores) || n_cores < 1L) {
  stop("n_cores must be a positive integer.")
}
if (length(n_calib) != 1L || is.na(n_calib) ||
    !n_calib %in% c(0L, 10L, 25L, 50L)) {
  stop("n_calib must be one of 0, 10, 25, or 50.")
}
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
out_root <- normalizePath(out_root, mustWork = TRUE)

if (!requireNamespace("eivGP", quietly = TRUE)) {
  stop("The eivGP package is not installed in the active R library.")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("jsonlite is required to read meta.json.")
}
library(eivGP)
pkg_version <- utils::packageVersion("eivGP")
if (pkg_version < "0.3.1") {
  stop("This runner requires eivGP >= 0.3.1; found ", pkg_version)
}

study_path <- file.path(data_dir, "study1_data.csv")
meta_path <- file.path(data_dir, "meta.json")
if (!file.exists(study_path) || !file.exists(meta_path)) {
  stop("The frozen study data or metadata file is missing.")
}
dat <- utils::read.csv(study_path, stringsAsFactors = FALSE)
meta <- jsonlite::fromJSON(meta_path, simplifyVector = TRUE)
train <- dat[dat$split == "train", , drop = FALSE]
test <- dat[dat$split == "test", , drop = FALSE]
if (nrow(train) != 100L || nrow(test) != 300L || as.integer(meta$m) != 6L) {
  stop("The frozen 100/300, six-class design was not recovered.")
}
calib_name <- if (n_calib == 0L) NA_character_ else paste0("calib_", n_calib)
if (n_calib == 0L) {
  calib_idx <- integer()
} else {
  if (!calib_name %in% names(train)) {
    stop("study1_data.csv has no ", calib_name, " column.")
  }
  calib_idx <- which(train[[calib_name]] == 1L)
  if (length(calib_idx) != n_calib) {
    stop("Expected ", n_calib, " fixed-design calibration rows.")
  }
}
X_train <- as.matrix(train[, c("x_temperature", "x_salinity"), drop = FALSE])
X_test <- as.matrix(test[, c("x_temperature", "x_salinity"), drop = FALSE])
C_train <- matrix(as.integer(train$c), ncol = 1L,
                  dimnames = list(NULL, "phosphate_class"))
C_test <- matrix(as.integer(test$c), ncol = 1L,
                 dimnames = list(NULL, "phosphate_class"))
U_obs <- matrix(NA_real_, nrow(train), 1L,
                dimnames = list(NULL, "log_phosphate_std"))
U_obs[calib_idx, 1L] <- train$u_log_std[calib_idx]

if (mode == "full") {
  n_iter <- if (n_calib == 0L) 30000L else 20000L
  burn <- if (n_calib == 0L) 7500L else 5000L
  n_chains <- 4L
  n_pred_draw <- 600L
} else {
  n_iter <- 300L
  burn <- 80L
  n_chains <- 2L
  n_pred_draw <- 80L
}
checkpoint_every <- 2000L
parallel_fit <- n_cores > 1L
mcmc_seed <- 2026090901L
out_dir <- file.path(out_root, mode)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
status_path <- file.path(out_dir, "status.txt")
complete_path <- file.path(out_dir, "complete.flag")
gate_pass_path <- file.path(out_dir, "gate_pass.flag")
gate_fail_path <- file.path(out_dir, "gate_fail.flag")
gate_decision_path <- file.path(out_dir, "single_gate_decision.csv")

## Rank-normalized tail ESS is not defined for some binary functionals. Allow
## that specific limitation only when every affected functional is a
## dictionary-occupancy indicator and all of its other diagnostics pass.
evaluate_single_gate <- function(diagnostic) {
  if (identical(as.character(diagnostic$status), "screen_passed")) {
    return(list(pass = TRUE, rule = "package_screen_passed",
                detail = "all package target diagnostics passed"))
  }
  tab <- diagnostic$table
  required <- c(
    "parameter", "n_chains", "min_draws_per_chain", "rhat", "ess_bulk",
    "ess_tail", "mcse_sd_ratio", "target_pass", "status"
  )
  if (!is.data.frame(tab) || !all(required %in% names(tab))) {
    return(list(pass = FALSE, rule = "package_status",
                detail = "diagnostic table is unavailable or incomplete"))
  }
  bad <- tab[!tab$target_pass, , drop = FALSE]
  allowed_discrete_tail_na <-
    nrow(bad) > 0L &&
    all(grepl("^dictionary_occupancy\\[", bad$parameter)) &&
    all(bad$n_chains >= 4L) &&
    all(bad$min_draws_per_chain >= 400L) &&
    all(is.finite(bad$rhat) & bad$rhat <= 1.01) &&
    all(is.finite(bad$ess_bulk) & bad$ess_bulk >= 400) &&
    all(is.na(bad$ess_tail)) &&
    all(is.finite(bad$mcse_sd_ratio) & bad$mcse_sd_ratio <= 0.05) &&
    all(bad$status == "insufficient_diagnostics")
  if (isTRUE(allowed_discrete_tail_na)) {
    return(list(
      pass = TRUE,
      rule = "discrete_tail_ess_exception",
      detail = paste0(
        nrow(bad),
        " binary dictionary-occupancy functionals had undefined tail ESS; ",
        "their R-hat, bulk ESS, and MCSE thresholds passed"
      )
    ))
  }
  list(pass = FALSE, rule = "package_status",
       detail = paste0("package diagnostic status: ", diagnostic$status))
}

publish_gate_decision <- function(diagnostic) {
  decision <- evaluate_single_gate(diagnostic)
  utils::write.csv(
    data.frame(
      package_status = as.character(diagnostic$status),
      gate_pass = decision$pass,
      gate_rule = decision$rule,
      detail = decision$detail,
      stringsAsFactors = FALSE
    ),
    gate_decision_path,
    row.names = FALSE
  )
  if (isTRUE(decision$pass)) {
    writeLines(decision$rule, gate_pass_path)
    if (file.exists(gate_fail_path)) unlink(gate_fail_path)
  } else {
    writeLines(as.character(diagnostic$status), gate_fail_path)
    if (file.exists(gate_pass_path)) unlink(gate_pass_path)
  }
  decision
}

write_status <- function(state, detail, iteration = NULL,
                         next_iteration = NULL) {
  lines <- c(
    paste0("state: ", state),
    paste0("updated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste0("mode: ", mode),
    paste0("n_calib: ", n_calib),
    paste0("target_n_iter: ", n_iter),
    paste0("burn: ", burn),
    paste0("n_chains: ", n_chains),
    paste0("requested_n_cores: ", n_cores),
    paste0("detail: ", detail)
  )
  if (!is.null(iteration)) lines <- c(lines, paste0("checkpoint_iteration: ", iteration))
  if (!is.null(next_iteration)) lines <- c(lines, paste0("next_iteration: ", next_iteration))
  writeLines(lines, status_path)
}

save_fit_atomic <- function(object, path) {
  tmp <- tempfile(".fit-v031-", tmpdir = out_dir, fileext = ".rds")
  saveRDS(object, tmp)
  if (!file.rename(tmp, path)) {
    if (file.exists(path)) unlink(path)
    if (!file.rename(tmp, path)) stop("Could not publish checkpoint: ", path)
  }
  invisible(path)
}

if (file.exists(complete_path)) {
  diagnostic_path <- file.path(out_dir, "target_diagnostics_v031.rds")
  if (file.exists(diagnostic_path)) {
    completed_diagnostic <- readRDS(diagnostic_path)
    completed_decision <- publish_gate_decision(completed_diagnostic)
    write_status(
      if (isTRUE(completed_decision$pass)) "DONE_GATE_PASSED" else "DONE_GATE_FAILED",
      paste0("existing single fit audited under ", completed_decision$rule),
      n_iter
    )
  }
  message("Already complete; gate decision refreshed: ", out_dir)
  quit(save = "no", status = 0L)
}

fit_path <- file.path(out_dir, "fit_v031.rds")
if (file.exists(fit_path)) {
  fit <- readRDS(fit_path)
  if (is.null(fit$checkpoint$iteration) ||
      !identical(as.integer(fit$checkpoint$burn), burn) ||
      length(fit$checkpoint$states) != n_chains ||
      as.integer(fit$checkpoint$iteration) > n_iter) {
    stop("The saved checkpoint is incompatible with this run.")
  }
  write_status("RESUMING", "loaded the saved current-package checkpoint",
               as.integer(fit$checkpoint$iteration))
} else {
  first_target <- min(n_iter, burn + 1L)
  write_status("MCMC", "starting a fresh eivGP 0.3.1 single fit",
               0L, first_target)
  fit <- fit_eivgp(
    X = X_train, y = train$y, C = C_train, U_obs = U_obs,
    engine = "univariate", m_vec = 6L, kernel = "se",
    standardize_U = FALSE,
    n_iter = first_target, burn = burn, n_chains = n_chains,
    parallel = parallel_fit,
    n_cores = if (parallel_fit) n_cores else NULL,
    seed = mcmc_seed, verbose = TRUE
  )
  save_fit_atomic(fit, fit_path)
  capture.output(summary(fit), file = file.path(out_dir, "summary_v031.txt"))
}

while (as.integer(fit$checkpoint$iteration) < n_iter) {
  current_iteration <- as.integer(fit$checkpoint$iteration)
  next_iteration <- min(n_iter, current_iteration + checkpoint_every)
  write_status("MCMC", "continuing the current-package checkpoint",
               current_iteration, next_iteration)
  fit <- continue_eivgp(
    fit, n_iter = next_iteration - current_iteration,
    parallel = parallel_fit,
    n_cores = if (parallel_fit) n_cores else NULL,
    verbose = TRUE
  )
  save_fit_atomic(fit, fit_path)
  capture.output(summary(fit), file = file.path(out_dir, "summary_v031.txt"))
  write_status("MCMC_CHECKPOINT", "checkpoint saved",
               as.integer(fit$checkpoint$iteration))
}

write_status("DIAGNOSTICS", "running the package target-level convergence audit",
             as.integer(fit$checkpoint$iteration))
if (length(calib_idx) >= 5L) {
  panel_idx <- calib_idx[unique(round(seq(
    1, length(calib_idx), length.out = 5L
  )))]
} else {
  ordered_u <- order(train$u_log_std)
  panel_idx <- ordered_u[unique(round(seq(
    1, length(ordered_u), length.out = 5L
  )))]
}
panel_U <- matrix(
  train$u_log_std[panel_idx], ncol = 1L,
  dimnames = list(NULL, "log_phosphate_std")
)
diagnostic <- tryCatch(
  diagnose_eivgp(
    fit,
    X = X_train[panel_idx, , drop = FALSE],
    C = C_train[panel_idx, , drop = FALSE],
    U = panel_U,
    n_latent = 64L,
    seed = 2026090902L
  ),
  error = function(e) e
)
if (inherits(diagnostic, "error")) {
  diagnostic_status <- data.frame(
    status = "diagnostic_error",
    recommendation = NA_character_,
    error = conditionMessage(diagnostic),
    stringsAsFactors = FALSE
  )
  gate_pass <- FALSE
} else {
  saveRDS(diagnostic, file.path(out_dir, "target_diagnostics_v031.rds"))
  write_diagnostics_eivgp(
    diagnostic,
    file.path(out_dir, "diagnostic_exports"),
    prefix = "ocean_single",
    overwrite = TRUE
  )
  diagnostic_status <- data.frame(
    status = diagnostic$status,
    recommendation = diagnostic$recommendation,
    error = NA_character_,
    stringsAsFactors = FALSE
  )
  gate_decision <- publish_gate_decision(diagnostic)
  gate_pass <- isTRUE(gate_decision$pass)
}
utils::write.csv(
  diagnostic_status,
  file.path(out_dir, "target_diagnostics_status.csv"),
  row.names = FALSE
)

write_status("PREDICTING", "drawing held-out response predictions",
             as.integer(fit$checkpoint$iteration))
n_saved <- nrow(fit$mcmc$samples_u)
draw_ids <- seq_len(n_saved)
if (length(draw_ids) > n_pred_draw) {
  set.seed(mcmc_seed + 2000L + n_calib)
  draw_ids <- sort(sample(draw_ids, n_pred_draw))
}
saveRDS(draw_ids, file.path(out_dir, "prediction_draw_ids.rds"))
pred_draws <- predict_eivgp(
  fit, new_X = X_test, new_C = C_test,
  target = "response", draw_ids = draw_ids,
  n_per_draw = 1L, joint = FALSE, seed = mcmc_seed + 3000L + n_calib
)
if (!is.matrix(pred_draws) || ncol(pred_draws) != nrow(test) ||
    any(!is.finite(pred_draws))) {
  stop("predict_eivgp returned invalid predictive draws.")
}
saveRDS(pred_draws, file.path(out_dir, "predictive_draws_v031.rds"))

crps_one <- function(x, y) {
  x <- sort(as.numeric(x))
  n <- length(x)
  pair_mean <- 2 * sum((2 * seq_len(n) - n - 1) * x) / n^2
  mean(abs(x - y)) - 0.5 * pair_mean
}
pred_q <- apply(pred_draws, 2L, stats::quantile,
                probs = c(0.025, 0.975), names = FALSE)
pred_mean <- colMeans(pred_draws)
pred_width <- pred_q[2L, ] - pred_q[1L, ]
pred_covered <- test$y >= pred_q[1L, ] & test$y <= pred_q[2L, ]
metrics <- data.frame(
  method = "EIV-GP-v0.3.1-fixed-single",
  n_calib = n_calib,
  RMSE = sqrt(mean((pred_mean - test$y)^2)),
  MAE = mean(abs(pred_mean - test$y)),
  CRPS = mean(vapply(seq_len(ncol(pred_draws)), function(j) {
    crps_one(pred_draws[, j], test$y[j])
  }, numeric(1L))),
  Coverage95 = mean(pred_covered),
  Width95 = mean(pred_width),
  IntervalScore95 = mean(
    pred_width + 40 * pmax(pred_q[1L, ] - test$y, 0) +
      40 * pmax(test$y - pred_q[2L, ], 0)
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(metrics, file.path(out_dir, "predictive_metrics.csv"), row.names = FALSE)

write_status("IMPUTING", "drawing masked training phosphate values",
             as.integer(fit$checkpoint$iteration))
miss <- which(!is.finite(U_obs[, 1L]))
u_array <- impute_eivgp(fit, rows = miss, draw_ids = draw_ids, scale = "model")
u_mat <- matrix(u_array[, , 1L], nrow = length(draw_ids), ncol = length(miss))
u_q <- apply(u_mat, 2L, stats::quantile,
             probs = c(0.025, 0.975), names = FALSE)
u_mean <- colMeans(u_mat)
u_truth <- train$u_log_std[miss]
u_summary <- data.frame(
  train_row_id = miss,
  source_row = train$source_row[miss],
  c = train$c[miss],
  true_u_log_std = u_truth,
  post_mean_u_log_std = u_mean,
  post_lo_u_log_std = u_q[1L, ],
  post_hi_u_log_std = u_q[2L, ],
  covered95 = u_truth >= u_q[1L, ] & u_truth <= u_q[2L, ],
  stringsAsFactors = FALSE
)
utils::write.csv(u_summary, file.path(out_dir, "u_imputation_summary.csv"),
                 row.names = FALSE)
u_metrics <- do.call(rbind, lapply(c(0L, seq_len(6L)), function(cc) {
  use <- if (cc == 0L) rep(TRUE, length(miss)) else train$c[miss] == cc
  data.frame(
    c = cc,
    n_hidden = sum(use),
    RMSE_u_log_std = sqrt(mean((u_mean[use] - u_truth[use])^2)),
    MAE_u_log_std = mean(abs(u_mean[use] - u_truth[use])),
    Coverage95_u_log_std = mean(u_summary$covered95[use]),
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(u_metrics, file.path(out_dir, "u_metrics_by_class.csv"),
                 row.names = FALSE)
if (!is.null(fit$diagnostics$summary)) {
  utils::write.csv(fit$diagnostics$summary,
                   file.path(out_dir, "mcmc_diagnostics.csv"), row.names = FALSE)
}

config <- data.frame(
  package_version = as.character(pkg_version),
  package_commit = Sys.getenv("EIVGP_PACKAGE_COMMIT", unset = "unknown"),
  R_version = R.version.string,
  mode = mode,
  n_iter = n_iter,
  burn = burn,
  n_chains = n_chains,
  n_cores = n_cores,
  checkpoint_every = checkpoint_every,
  mcmc_seed = mcmc_seed,
  calibration = if (n_calib == 0L) "none" else calib_name,
  gate_status = diagnostic_status$status,
  stringsAsFactors = FALSE
)
utils::write.csv(config, file.path(out_dir, "run_config.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
writeLines(capture.output(print(metrics)), file.path(out_dir, "metrics_print.txt"))
writeLines(c(
  paste0("completed=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  paste0("gate_status=", diagnostic_status$status),
  paste0("package_version=", pkg_version),
  paste0("package_commit=", Sys.getenv("EIVGP_PACKAGE_COMMIT", unset = "unknown"))
), complete_path)
if (isTRUE(gate_pass)) {
  if (!exists("gate_decision")) gate_decision <- publish_gate_decision(diagnostic)
  write_status("DONE_GATE_PASSED", "fixed-design single fit complete",
               as.integer(fit$checkpoint$iteration))
} else {
  if (!inherits(diagnostic, "error") && !exists("gate_decision")) {
    gate_decision <- publish_gate_decision(diagnostic)
  } else if (inherits(diagnostic, "error")) {
    writeLines(diagnostic_status$status, gate_fail_path)
  }
  write_status("DONE_GATE_FAILED", "fixed-design single fit complete; diagnostic gate did not pass",
               as.integer(fit$checkpoint$iteration))
}
print(metrics)
print(diagnostic_status)
