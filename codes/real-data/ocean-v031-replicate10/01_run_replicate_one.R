#!/usr/bin/env Rscript

## Run one calibration size for one Ocean replicate with the current public
## eivGP 0.3.1 API.  The script is deliberately independent of the absolute
## paths recorded in replicate_manifest.csv, because those paths belong to the
## original computer and are not valid inside Google Colab.
##
## Usage:
##   Rscript 01_run_replicate_one.R \
##     DATA_DIR REDRAW_DIR OUT_ROOT REP_ID N_CALIB MODE TARGET_DIAGNOSTICS N_CORES
##
## DATA_DIR defaults to ../ocean_data (codes/real-data/ocean_data).
## REDRAW_DIR is this bundle's calib_designs/ folder.
## MODE is "smoke" or "full".  TARGET_DIAGNOSTICS is 0 or 1.
## Set EIVGP_CHECKPOINT_EVERY (default 2000) to control continuation chunks.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5L) {
  stop(
    "Usage: Rscript 01_run_replicate_one.R ",
    "DATA_DIR REDRAW_DIR OUT_ROOT REP_ID N_CALIB ",
    "[smoke|full] [target_diagnostics 0|1] [n_cores]"
  )
}

data_dir <- normalizePath(args[[1L]], mustWork = TRUE)
redraw_dir <- normalizePath(args[[2L]], mustWork = TRUE)
out_root <- normalizePath(args[[3L]], mustWork = FALSE)
rep_id <- as.integer(args[[4L]])
n_calib <- as.integer(args[[5L]])
mode <- if (length(args) >= 6L) args[[6L]] else "full"
target_diagnostics <- if (length(args) >= 7L) as.logical(as.integer(args[[7L]])) else FALSE
n_cores_arg <- if (length(args) >= 8L) as.integer(args[[8L]]) else NA_integer_

if (!rep_id %in% seq_len(10L)) stop("REP_ID must be an integer from 1 to 10.")
if (!n_calib %in% c(10L, 25L, 50L)) stop("N_CALIB must be 10, 25, or 50.")
if (!mode %in% c("smoke", "full")) stop("MODE must be smoke or full.")
if (length(target_diagnostics) != 1L || is.na(target_diagnostics)) {
  stop("TARGET_DIAGNOSTICS must be 0 or 1.")
}

if (!requireNamespace("eivGP", quietly = TRUE)) {
  stop("The current eivGP package is not installed in the active R library.")
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
manifest_path <- file.path(redraw_dir, "replicate_manifest.csv")
allocation_candidates <- c(
  file.path(redraw_dir, "class6_proportional_calibration_allocation.csv"),
  file.path(data_dir, "class6_proportional_calibration_allocation.csv")
)
allocation_path <- allocation_candidates[file.exists(allocation_candidates)][1L]
if (is.na(allocation_path) || !nzchar(allocation_path)) {
  stop("Missing class6_proportional_calibration_allocation.csv")
}
for (path in c(study_path, meta_path, manifest_path, allocation_path)) {
  if (!file.exists(path)) stop("Missing required input: ", path)
}

dat <- utils::read.csv(study_path, stringsAsFactors = FALSE)
meta <- jsonlite::fromJSON(meta_path, simplifyVector = TRUE)
manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
allocation <- utils::read.csv(allocation_path, stringsAsFactors = FALSE)

required_data <- c(
  "source_row", "split", "x_temperature", "x_salinity", "c", "y", "u_log_std"
)
missing_data <- setdiff(required_data, names(dat))
if (length(missing_data)) stop("study1_data.csv is missing: ", paste(missing_data, collapse = ", "))
if (nrow(dat) != 400L || sum(dat$split == "train") != 100L || sum(dat$split == "test") != 300L) {
  stop("study1_data.csv must contain the frozen 100/300 train/test cohort.")
}
if (!all(as.integer(dat$c) %in% seq_len(as.integer(meta$m)))) {
  stop("Ordinal labels in study1_data.csv do not match meta$m.")
}
observed_train_counts <- as.integer(table(factor(dat$c[dat$split == "train"], levels = seq_len(as.integer(meta$m)))))
if (!identical(as.integer(allocation$train_n), observed_train_counts)) {
  stop("Calibration allocation train_n does not match study1_data.csv.")
}
if (sum(allocation$n_calib_10) != 10L || sum(allocation$n_calib_25) != 25L ||
    sum(allocation$n_calib_50) != 50L) {
  stop("Calibration allocation totals must be 10, 25, and 50.")
}

manifest_row <- manifest[manifest$rep_id == rep_id, , drop = FALSE]
if (nrow(manifest_row) != 1L) stop("Could not resolve exactly one manifest row for rep_id=", rep_id)
calib_seed <- as.integer(manifest_row$calib_seed[[1L]])
mcmc_seed <- as.integer(manifest_row$mcmc_seed[[1L]])
design_dir <- file.path(redraw_dir, sprintf("rep%02d_seed%d", rep_id, calib_seed))
flags_path <- file.path(design_dir, "train_calib_flags.csv")
if (!file.exists(flags_path)) stop("Missing replicate calibration flags: ", flags_path)
flags <- utils::read.csv(flags_path, stringsAsFactors = FALSE)

train <- dat[dat$split == "train", , drop = FALSE]
test <- dat[dat$split == "test", , drop = FALSE]
flag_train <- flags[match(train$source_row, flags$source_row), , drop = FALSE]
if (anyNA(flag_train$source_row)) stop("Calibration flags do not cover all training source rows.")
flag_name <- paste0("calib_", n_calib)
if (!flag_name %in% names(flag_train)) stop("Missing calibration flag column: ", flag_name)
calib_idx <- which(flag_train[[flag_name]] == 1L)
if (length(calib_idx) != n_calib) stop("Expected ", n_calib, " calibration rows; found ", length(calib_idx))

X_train <- as.matrix(train[, c("x_temperature", "x_salinity"), drop = FALSE])
X_test <- as.matrix(test[, c("x_temperature", "x_salinity"), drop = FALSE])
C_train <- matrix(as.integer(train$c), ncol = 1L,
                  dimnames = list(NULL, "phosphate_class"))
C_test <- matrix(as.integer(test$c), ncol = 1L,
                 dimnames = list(NULL, "phosphate_class"))
U_obs <- matrix(NA_real_, nrow(train), 1L,
                dimnames = list(NULL, "log_phosphate_std"))
U_obs[calib_idx, 1L] <- train$u_log_std[calib_idx]

if (is.na(n_cores_arg)) {
  detected <- parallel::detectCores(logical = TRUE)
  n_cores <- max(1L, min(4L, ifelse(is.na(detected), 1L, detected)))
} else {
  n_cores <- max(1L, n_cores_arg)
}

if (mode == "full") {
  n_iter <- 20000L
  burn <- 5000L
  n_chains <- 4L
  n_pred_draw <- 600L
} else {
  n_iter <- 300L
  burn <- 80L
  n_chains <- 2L
  n_pred_draw <- 80L
}
checkpoint_every <- as.integer(Sys.getenv("EIVGP_CHECKPOINT_EVERY", unset = "2000"))
if (length(checkpoint_every) != 1L || is.na(checkpoint_every) || checkpoint_every < 1L) {
  stop("EIVGP_CHECKPOINT_EVERY must be a positive integer.")
}
parallel_fit <- n_cores > 1L

rep_out <- file.path(out_root, sprintf("rep%02d", rep_id), paste0("c", n_calib))
dir.create(rep_out, recursive = TRUE, showWarnings = FALSE)
status_path <- file.path(rep_out, "status.txt")
complete_path <- file.path(rep_out, "complete.flag")
write_status <- function(state, detail, iteration = NULL, next_iteration = NULL,
                         actual_backend = NULL, actual_cores = NULL) {
  lines <- c(
    paste0("state: ", state),
    paste0("updated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste0("rep_id: ", rep_id),
    paste0("n_calib: ", n_calib),
    paste0("target_n_iter: ", n_iter),
    paste0("burn: ", burn),
    paste0("n_chains: ", n_chains),
    paste0("checkpoint_every: ", checkpoint_every),
    paste0("parallel_requested: ", parallel_fit),
    paste0("requested_n_cores: ", n_cores),
    paste0("detail: ", detail)
  )
  if (!is.null(iteration)) lines <- c(lines, paste0("checkpoint_iteration: ", iteration))
  if (!is.null(next_iteration)) lines <- c(lines, paste0("next_checkpoint_iteration: ", next_iteration))
  if (!is.null(actual_backend)) lines <- c(lines, paste0("actual_parallel_backend: ", actual_backend))
  if (!is.null(actual_cores)) lines <- c(lines, paste0("actual_parallel_cores: ", actual_cores))
  writeLines(lines, status_path)
}

## The requested n_cores value is a worker budget.  The package records the
## backend and effective worker count in the fit diagnostics; expose those
## values in status.txt/config so a Colab run can be audited directly.
fit_runtime_info <- function(object) {
  summary <- tryCatch(object$diagnostics$summary, error = function(e) NULL)
  if (!is.data.frame(summary) || nrow(summary) < 1L) {
    return(list(backend = NA_character_, cores = NA_integer_))
  }
  backend <- if ("parallel_backend" %in% names(summary)) {
    as.character(summary$parallel_backend[[1L]])
  } else {
    NA_character_
  }
  cores <- if ("parallel_cores" %in% names(summary)) {
    suppressWarnings(as.integer(summary$parallel_cores[[1L]]))
  } else {
    NA_integer_
  }
  list(backend = backend, cores = cores)
}
if (file.exists(complete_path)) {
  message("Already complete: ", rep_out)
  quit(save = "no", status = 0L)
}

write_status("RUNNING", "preparing current-package v0.3.1 fit")
fit_path <- file.path(rep_out, "fit_v031.rds")
fit <- NULL
save_fit_atomic <- function(object, path) {
  tmp <- tempfile(".fit-v031-", tmpdir = rep_out, fileext = ".rds")
  saveRDS(object, tmp)
  if (!file.rename(tmp, path)) {
    if (file.exists(path)) unlink(path)
    if (!file.rename(tmp, path)) stop("Could not publish checkpoint: ", path)
  }
  invisible(path)
}
if (file.exists(fit_path)) {
  fit <- readRDS(fit_path)
  if (is.null(fit$checkpoint) || is.null(fit$checkpoint$iteration)) {
    stop("Existing fit_v031.rds has no v0.3.1 continuation checkpoint. Use a fresh full output directory.")
  }
  if (!identical(as.integer(fit$checkpoint$burn), as.integer(burn)) ||
      length(fit$checkpoint$states) != n_chains ||
      as.integer(fit$checkpoint$iteration) > n_iter) {
    stop("Existing fit_v031.rds is incompatible with this run's full/smoke settings. Use a separate output directory.")
  }
  runtime_info <- fit_runtime_info(fit)
  write_status("RESUMING", "loaded fit_v031.rds; checking the next MCMC checkpoint",
               iteration = as.integer(fit$checkpoint$iteration),
               actual_backend = runtime_info$backend,
               actual_cores = runtime_info$cores)
} else {
  ## Burn-in must be completed before the first usable retained draw.  The
  ## first checkpoint is therefore burn + 1, followed by 2000-transition
  ## checkpoints.  Subsequent calls use the package's public continuation API.
  first_target <- min(n_iter, burn + 1L)
  write_status("MCMC", paste0("starting fit through iteration ", first_target),
               iteration = 0L, next_iteration = first_target)
  fit <- fit_eivgp(
    X = X_train,
    y = train$y,
    C = C_train,
    U_obs = U_obs,
    engine = "univariate",
    m_vec = 6L,
    kernel = "se",
    standardize_U = FALSE,
    n_iter = first_target,
    burn = burn,
    n_chains = n_chains,
    parallel = parallel_fit,
    n_cores = if (parallel_fit) n_cores else NULL,
    ## Keep the frozen replicate MCMC seed unchanged across C10/C25/C50,
    ## matching the existing replicate design; prediction and draw-selection
    ## seeds are offset by n_calib below.
    seed = mcmc_seed,
    verbose = TRUE
  )
  save_fit_atomic(fit, fit_path)
  capture.output(summary(fit), file = file.path(rep_out, "summary_v031.txt"))
  runtime_info <- fit_runtime_info(fit)
  write_status("MCMC_CHECKPOINT", "initial burn-in checkpoint saved",
               iteration = as.integer(fit$checkpoint$iteration),
               actual_backend = runtime_info$backend,
               actual_cores = runtime_info$cores)
}

while (as.integer(fit$checkpoint$iteration) < n_iter) {
  current_iteration <- as.integer(fit$checkpoint$iteration)
  next_iteration <- min(n_iter, current_iteration + checkpoint_every)
  write_status(
    "MCMC",
    paste0("continuing from iteration ", current_iteration, " to ", next_iteration),
    iteration = current_iteration,
    next_iteration = next_iteration
  )
  fit <- continue_eivgp(
    fit,
    n_iter = next_iteration - current_iteration,
    parallel = parallel_fit,
    n_cores = if (parallel_fit) n_cores else NULL,
    verbose = TRUE
  )
  save_fit_atomic(fit, fit_path)
  capture.output(summary(fit), file = file.path(rep_out, "summary_v031.txt"))
  runtime_info <- fit_runtime_info(fit)
  write_status(
    "MCMC_CHECKPOINT",
    "checkpoint saved; the next chunk can resume from fit_v031.rds",
    iteration = as.integer(fit$checkpoint$iteration),
    actual_backend = runtime_info$backend,
    actual_cores = runtime_info$cores
  )
}
runtime_info <- fit_runtime_info(fit)
write_status("FIT_DONE", "full MCMC target reached; running prediction and imputation",
             iteration = as.integer(fit$checkpoint$iteration),
             actual_backend = runtime_info$backend,
             actual_cores = runtime_info$cores)

## Use a fixed subset of retained posterior draws for scoring.  The draw IDs
## are saved so a rerun can be audited exactly.
n_saved <- nrow(fit$mcmc$samples_u)
draw_ids <- seq_len(n_saved)
if (length(draw_ids) > n_pred_draw) {
  set.seed(mcmc_seed + 2000L + n_calib)
  draw_ids <- sort(sample(draw_ids, n_pred_draw))
}
saveRDS(draw_ids, file.path(rep_out, "prediction_draw_ids.rds"))

write_status("PREDICTING", "calling predict_eivgp(target='response')",
             iteration = as.integer(fit$checkpoint$iteration),
             actual_backend = runtime_info$backend,
             actual_cores = runtime_info$cores)
pred_draws <- predict_eivgp(
  fit,
  new_X = X_test,
  new_C = C_test,
  target = "response",
  draw_ids = draw_ids,
  n_per_draw = 1L,
  joint = FALSE,
  seed = mcmc_seed + 3000L + n_calib
)
if (!is.matrix(pred_draws) || ncol(pred_draws) != nrow(test) || any(!is.finite(pred_draws))) {
  stop("predict_eivgp returned an invalid predictive-draw matrix.")
}
saveRDS(pred_draws, file.path(rep_out, "predictive_draws_v031.rds"))

crps_one <- function(x, y) {
  x <- sort(as.numeric(x))
  n <- length(x)
  pair_mean <- 2 * sum((2 * seq_len(n) - n - 1) * x) / n^2
  mean(abs(x - y)) - 0.5 * pair_mean
}
q <- apply(pred_draws, 2L, stats::quantile, probs = c(0.025, 0.975), names = FALSE)
pred_mean <- colMeans(pred_draws)
width <- q[2L, ] - q[1L, ]
covered <- test$y >= q[1L, ] & test$y <= q[2L, ]
interval_score <- width + 40 * pmax(q[1L, ] - test$y, 0) + 40 * pmax(test$y - q[2L, ], 0)
metrics <- data.frame(
  rep_id = rep_id,
  n_calib = n_calib,
  method = "EIV-GP-v0.3.1",
  RMSE = sqrt(mean((pred_mean - test$y)^2)),
  MAE = mean(abs(pred_mean - test$y)),
  Coverage95 = mean(covered),
  Width95 = mean(width),
  CRPS = mean(vapply(seq_len(ncol(pred_draws)), function(j) crps_one(pred_draws[, j], test$y[j]), numeric(1))),
  IntervalScore95 = mean(interval_score),
  stringsAsFactors = FALSE
)
utils::write.csv(metrics, file.path(rep_out, "predictive_metrics.csv"), row.names = FALSE)

write_status("IMPUTING", "calling impute_eivgp(scale='model') for hidden training U",
             iteration = as.integer(fit$checkpoint$iteration),
             actual_backend = runtime_info$backend,
             actual_cores = runtime_info$cores)
miss <- which(!is.finite(U_obs[, 1L]))
u_draws <- impute_eivgp(fit, rows = miss, draw_ids = draw_ids, scale = "model")
if (length(dim(u_draws)) != 3L || dim(u_draws)[2L] != length(miss)) {
  stop("impute_eivgp returned an unexpected array shape.")
}
u_mat <- matrix(u_draws[, , 1L], nrow = length(draw_ids), ncol = length(miss))
u_truth <- train$u_log_std[miss]
u_mean <- colMeans(u_mat)
u_q <- apply(u_mat, 2L, stats::quantile, probs = c(0.025, 0.975), names = FALSE)
u_summary <- data.frame(
  rep_id = rep_id,
  n_calib = n_calib,
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
utils::write.csv(u_summary, file.path(rep_out, "u_imputation_summary.csv"), row.names = FALSE)
u_metrics <- do.call(rbind, lapply(c(0L, seq_len(6L)), function(cc) {
  use <- if (cc == 0L) rep(TRUE, length(miss)) else train$c[miss] == cc
  data.frame(
    rep_id = rep_id,
    n_calib = n_calib,
    c = cc,
    n_hidden = sum(use),
    RMSE_u_log_std = sqrt(mean((u_mean[use] - u_truth[use])^2)),
    MAE_u_log_std = mean(abs(u_mean[use] - u_truth[use])),
    Coverage95_u_log_std = mean(u_summary$covered95[use]),
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(u_metrics, file.path(rep_out, "u_metrics_by_class.csv"), row.names = FALSE)

if (!is.null(fit$diagnostics$summary)) {
  utils::write.csv(fit$diagnostics$summary,
                   file.path(rep_out, "mcmc_diagnostics.csv"), row.names = FALSE)
}

## The target-level audit is optional because it is an additional calculation
## for every retained draw.  It uses only the current package diagnostics API.
if (isTRUE(target_diagnostics)) {
  write_status("DIAGNOSTICS", "calling diagnose_eivgp and write_diagnostics_eivgp",
               iteration = as.integer(fit$checkpoint$iteration),
               actual_backend = runtime_info$backend,
               actual_cores = runtime_info$cores)
  panel_idx <- calib_idx[unique(round(seq(1, length(calib_idx), length.out = min(5L, length(calib_idx)))))]
  diag_dir <- file.path(rep_out, "diagnostic_exports")
  diag_result <- tryCatch({
    d <- diagnose_eivgp(
      fit,
      X = X_train[panel_idx, , drop = FALSE],
      C = C_train[panel_idx, , drop = FALSE],
      U = U_obs[panel_idx, , drop = FALSE],
      n_latent = 64L,
      seed = mcmc_seed + 4000L + n_calib
    )
    saveRDS(d, file.path(rep_out, "target_diagnostics_v031.rds"))
    write_diagnostics_eivgp(d, diag_dir, prefix = "ocean", overwrite = TRUE)
    data.frame(status = d$status, recommendation = d$recommendation,
               error = NA_character_, stringsAsFactors = FALSE)
  }, error = function(e) {
    data.frame(status = NA_character_, recommendation = NA_character_,
               error = conditionMessage(e), stringsAsFactors = FALSE)
  })
  utils::write.csv(diag_result, file.path(rep_out, "target_diagnostics_status.csv"), row.names = FALSE)
}

config <- data.frame(
  rep_id = rep_id,
  n_calib = n_calib,
  calib_seed = calib_seed,
  mcmc_seed = mcmc_seed,
  mode = mode,
  package_version = as.character(pkg_version),
  package_commit = Sys.getenv("EIVGP_PACKAGE_COMMIT", unset = "unknown"),
  R_version = R.version.string,
  platform = R.version$platform,
  engine = "univariate",
  kernel = "se",
  standardize_U = FALSE,
  n_iter = n_iter,
  burn = burn,
  n_chains = n_chains,
  checkpoint_every = checkpoint_every,
  parallel = parallel_fit,
  n_cores = if (parallel_fit) n_cores else 1L,
  actual_parallel_backend = runtime_info$backend,
  actual_parallel_cores = runtime_info$cores,
  target_diagnostics = target_diagnostics,
  data_dir = data_dir,
  design_dir = design_dir,
  stringsAsFactors = FALSE
)
utils::write.csv(config, file.path(rep_out, "run_config.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(rep_out, "sessionInfo.txt"))
writeLines(capture.output(print(metrics)), file.path(rep_out, "metrics_print.txt"))
writeLines(c(
  paste0("rep_id=", rep_id), paste0("n_calib=", n_calib),
  paste0("calib_seed=", calib_seed), paste0("mcmc_seed=", mcmc_seed),
  paste0("package_version=", pkg_version),
  paste0("package_commit=", Sys.getenv("EIVGP_PACKAGE_COMMIT", unset = "unknown")),
  paste0("completed=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))
), complete_path)
write_status("DONE", "fit, prediction, latent-input imputation, and requested diagnostics saved",
             iteration = as.integer(fit$checkpoint$iteration),
             actual_backend = runtime_info$backend,
             actual_cores = runtime_info$cores)
print(metrics)
