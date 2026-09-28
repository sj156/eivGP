# ADNI implementation extracted from the former analysis notebook.
# Called by adni_run_analysis(); reader-facing ADNI work is deferred.
knitr::opts_chunk$set(echo = TRUE, message = TRUE, warning = TRUE)

## Environment variables override the notebook defaults without editing it:
## EIVGP_CORES=12 or EIVGP_SMOKE_TEST=1. Five folds and four chains are fixed.
env_int <- function(name, default, minimum = 0L) {
  raw <- Sys.getenv(name, unset = as.character(default))
  if (!grepl("^[0-9]+$", raw)) stop(name, " must be an integer.")
  x <- suppressWarnings(as.integer(raw))
  if (length(x) != 1L || is.na(x) || x < minimum) stop(name, " is invalid.")
  x
}
env_flag <- function(name, default = FALSE) {
  x <- tolower(Sys.getenv(name, unset = if (default) "true" else "false"))
  if (!x %in% c("true", "false", "1", "0", "yes", "no")) stop(name, " is invalid.")
  x %in% c("true", "1", "yes")
}

PREPARE_ONLY <- isTRUE(get0("ADNI_PREPARE_ONLY", inherits = FALSE))
CV_FOLDS <- 5L
TOTAL_CORES <- env_int("EIVGP_CORES", env_int("EIVGP_N_CORES", 4L, 1L), 1L)
MAX_FOLD_WORKERS <- env_int("ADNI_MAX_FOLD_WORKERS", 5L, 1L)
SMOKE_TEST <- env_flag("EIVGP_SMOKE_TEST", FALSE)


if (SMOKE_TEST) {
  N_CHAINS <- 4L; BURN <- 2L; INITIAL_ITER <- 8L
  EXTENSION_ITER <- 4L; CHECKPOINT_EVERY <- 3L
  FOLDS_TO_RUN <- if (env_flag("ADNI_SMOKE_ALL_FOLDS", FALSE)) 1:5 else 1L
  PREDICTION_DRAWS <- 8L
} else {
  N_CHAINS <- 4L; BURN <- 4000L; INITIAL_ITER <- 20000L
  EXTENSION_ITER <- 10000L; CHECKPOINT_EVERY <- 2000L
  FOLDS_TO_RUN <- 1:5
  PREDICTION_DRAWS <- 1000L
}

## Frozen sampler choice. A matched ADNI pilot found that block sizes 16 and
## 32 reduced wall time per transition but worsened missing-U mixing enough
## that they did not improve convergence per equal waiting time. Keep the
## package reference schedule for the confirmatory run.
SAMPLER_CONTROL <- list(
  u_block_size = 8L,
  cross_every = 10L,
  ess_max_steps = 10000L,
  dictionary_mode = "conditional",
  max_dictionary_size = 125L,
  gp_block_schur = TRUE,
  threshold_update = "gibbs"
)

R_HAT_LIMIT <- 1.01
BULK_ESS_LIMIT <- 400
TAIL_ESS_LIMIT <- 400
MCSE_SD_RATIO_LIMIT <- 0.10

input_path <- tryCatch(knitr::current_input(dir = TRUE), error = function(e) "")
if (length(input_path) != 1L || is.na(input_path)) input_path <- ""
PROJECT_DIR <- Sys.getenv("ADNI_PROJECT_DIR", if (nzchar(input_path)) normalizePath(file.path(dirname(input_path), "..", "..", "codes", "applications")) else normalizePath("codes/applications"))
source(file.path(PROJECT_DIR, "real_data_paths.R"), local = TRUE)
DATA_DIR <- resolve_adni_data_dir(PROJECT_DIR)
source(file.path(PROJECT_DIR, "ADNI_cv.R"), local = TRUE)
source(file.path(PROJECT_DIR, "ADNI_parallel.R"), local = TRUE)
CORE_PLAN <- if (PREPARE_ONLY) get("ADNI_RUN_CORE_PLAN", inherits = FALSE) else
  adni_core_plan(TOTAL_CORES, length(FOLDS_TO_RUN), N_CHAINS, MAX_FOLD_WORKERS)
N_CORES <- CORE_PLAN$chain_workers
adni_print_core_plan(CORE_PLAN)
DATA_FILE <- file.path(DATA_DIR, "toledo_adni_cohort_n495.csv")
if (exists("VALIDATION_ID")) {
  VALIDATION_PLAN_FILE <- file.path(
    DATA_DIR,
    "validation-splits",
    "validation_split_plan.csv"
  )
  validation_plan <- read.csv(VALIDATION_PLAN_FILE, stringsAsFactors = FALSE)
  VALIDATION_META <- validation_plan[validation_plan$validation_id == VALIDATION_ID, , drop = FALSE]
  if (nrow(VALIDATION_META) != 1L) stop("Validation ID is absent or duplicated in validation_split_plan.csv")
  FOLD_FILE <- file.path(DATA_DIR, "validation-splits", sprintf("validation_%02d.csv", VALIDATION_ID))
  OUTPUT_NAME <- sprintf("validation_%02d", VALIDATION_ID)
  if (SMOKE_TEST) OUTPUT_NAME <- paste0("smoke_", OUTPUT_NAME)
  OUTPUT_BASE <- Sys.getenv("ADNI_OUTPUT_DIR", application_paths(PROJECT_DIR)$adni_output)
  OUTPUT_ROOT <- file.path(OUTPUT_BASE, OUTPUT_NAME)
} else {
  OUTPUT_BASE <- application_paths(PROJECT_DIR)$adni_output
  OUTPUT_ROOT <- file.path(OUTPUT_BASE, if (SMOKE_TEST) "smoke" else "five_fold")
  FOLD_FILE <- file.path(OUTPUT_ROOT, "appendix", "fold_assignments.csv")
}

dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)
RUN_LOCK <- file.path(OUTPUT_ROOT, ".run-lock")
LOCK_OWNED <- FALSE
if (!dir.create(RUN_LOCK, showWarnings = FALSE)) stop("Run already locked: ", RUN_LOCK,
  ". If a job was killed, verify it has stopped before removing this directory.")
LOCK_OWNED <- TRUE
writeLines(paste(Sys.info()[["nodename"]], Sys.getpid()), file.path(RUN_LOCK, "owner"))


## Dependency installation is a separate, explicit setup step.
EIVGP_TESTED_COMMIT <- "a3aa0f697f97d4bdf027507ac3c9b952644affb9"
if (!requireNamespace("eivGP", quietly = TRUE)) stop("Run setup_real_data.R before analysis.")
pkg <- utils::packageDescription("eivGP")
if (!identical(pkg$Version, "0.3.1")) {
  stop("eivGP version 0.3.1 is required. Run setup_real_data.R.")
}
suppressPackageStartupMessages(library(eivGP))

source(file.path(PROJECT_DIR, "ADNI_case_study_helpers.R"), local = TRUE)
for (d in c("main", "appendix")) dir.create(file.path(OUTPUT_ROOT, d), showWarnings = FALSE)
case_packages <- adni_case_preflight(SMOKE_TEST)
adni_write_csv(as.data.frame(CORE_PLAN), file.path(OUTPUT_ROOT, "appendix", "resource_plan.csv"))
package_files <- list.files(system.file("R", package = "eivGP"), full.names = TRUE)
adni_write_csv(data.frame(file = basename(package_files), md5 = unname(tools::md5sum(package_files))),
  file.path(OUTPUT_ROOT, "appendix", "eivGP_code_fingerprints.csv"))
adni_write_csv(case_packages, file.path(OUTPUT_ROOT, "appendix", "package_versions.csv"))
adni_write_csv(data.frame(schema = ADNI_CASE_SCHEMA, smoke = SMOKE_TEST,
  prediction_draws = PREDICTION_DRAWS, run_mode = "five_fold", cv_folds = CV_FOLDS, split_seed = ADNI_CV_SEED,
  primary_target = "Y | test X,C (no test CSF)",
  competitors = paste(ADNI_METHODS, collapse = ";")),
  file.path(OUTPUT_ROOT, "appendix", "case_config.csv"))

cat("Project:", PROJECT_DIR, "\n")
cat("Data:", DATA_DIR, "\n")
cat("Five-fold CV | folds:", paste(FOLDS_TO_RUN, collapse = ","), "\n")
cat("Mode:", if (SMOKE_TEST) "SMOKE" else "PRODUCTION", "\n")
cat("eivGP:", as.character(utils::packageVersion("eivGP")), "\n")

dat <- adni_read_cohort(DATA_DIR)
folds <- adni_make_folds(dat)
if (file.exists(FOLD_FILE)) {
  previous <- read.csv(FOLD_FILE)
  if (!isTRUE(all.equal(previous, folds, check.attributes = FALSE)))
    stop("Saved fold assignment differs. Use a separate output directory.")
} else adni_write_csv(folds, FOLD_FILE)
adni_write_csv(as.data.frame(table(folds$fold, folds$R)),
               file.path(OUTPUT_ROOT, "appendix", "fold_CSF_counts.csv"))

X_COLUMNS <- c("x_age_years", "x_female", "x_apoe4_dose")
C_COLUMNS <- c("c_ab42_ab40_3", "c_gfap_3")
Y_COLUMN <- "y_centiloid"
U_COLUMN <- "u_csf_abeta42"

make_C <- function(z) {
  out <- data.frame(
    c_ab42_ab40_3 = ordered(z$c_ab42_ab40_3, levels = 1:3),
    c_gfap_3 = ordered(z$c_gfap_3, levels = 1:3),
    check.names = FALSE
  )
  names(out) <- C_COLUMNS
  out
}

adni_eda(dat, OUTPUT_ROOT, SMOKE_TEST)


atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(paste0(basename(path), "-"), tmpdir = dirname(path))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(object, tmp)
  if (!file.rename(tmp, path)) stop("Could not publish ", path)
}

atomic_write_lines <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(paste0(basename(path), "-"), tmpdir = dirname(path))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  writeLines(lines, tmp)
  if (!file.rename(tmp, path)) stop("Could not publish ", path)
}

stamp <- function(run_dir, message) {
  line <- sprintf("- [%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), message)
  cat(line, "\n", file = file.path(run_dir, "PROGRESS.md"), append = TRUE)
  cat(line, "\n")
  flush.console()
}

set_status <- function(run_dir, fold, stage, detail) {
  atomic_write_lines(c(
    paste0("updated_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("cv_folds=", CV_FOLDS), paste0("fit_unit=", fold),
    paste0("stage=", stage), paste0("detail=", detail)
  ), file.path(run_dir, "CURRENT_STATUS.txt"))
  stamp(run_dir, paste("Five-fold CV", "|", stage, "|", detail))
}

light_diagnostics <- function(fit, fold_dir, stage) {
  common <- c("parameter", "rhat", "ess_bulk", "ess_tail", "mcse_sd_ratio")
  parts <- list(
    hyper = subset(fit$diagnostics$rhat_hyper,
                   !grepl("^(J\\[|dictionary_occupancy)", parameter))[, common],
    loading = fit$diagnostics$rhat_A[, common],
    threshold = fit$diagnostics$rhat_tau[, common],
    missing_U = fit$diagnostics$rhat_U[, common]
  )
  tab <- do.call(rbind, lapply(names(parts), function(g) {
    z <- parts[[g]]; z$group <- g; z
  }))
  rownames(tab) <- NULL
  finite <- with(tab, is.finite(rhat) & is.finite(ess_bulk) &
                   is.finite(ess_tail) & is.finite(mcse_sd_ratio))
  passed <- all(finite) && all(tab$rhat <= R_HAT_LIMIT) &&
    all(tab$ess_bulk >= BULK_ESS_LIMIT) && all(tab$ess_tail >= TAIL_ESS_LIMIT) &&
    all(tab$mcse_sd_ratio <= MCSE_SD_RATIO_LIMIT)
  summary <- data.frame(
    stage = stage,
    n_iter = fit$diagnostics$summary$n_iter,
    saved_per_chain = fit$diagnostics$summary$saved_per_chain,
    max_rhat = if (all(finite)) max(tab$rhat) else NA_real_,
    min_bulk_ess = if (all(finite)) min(tab$ess_bulk) else NA_real_,
    min_tail_ess = if (all(finite)) min(tab$ess_tail) else NA_real_,
    max_mcse_sd_ratio = if (all(finite)) max(tab$mcse_sd_ratio) else NA_real_,
    light_gate_passed = passed,
    stringsAsFactors = FALSE
  )
  write.csv(tab, file.path(fold_dir, paste0("diagnostics_", stage, "_parameters.csv")),
            row.names = FALSE)
  write.csv(summary, file.path(fold_dir, paste0("diagnostics_", stage, "_summary.csv")),
            row.names = FALSE)
  list(table = tab, summary = summary, passed = passed)
}

trace_series <- function(fit) {
  out <- list(V = fit$mcmc$samples_by_chain$V, r = fit$mcmc$samples_by_chain$r)
  H <- fit$mcmc$samples_by_chain$logtheta
  for (j in seq_len(ncol(H[[1L]]))) {
    nm <- colnames(H[[1L]])[j]
    out[[nm]] <- lapply(H, function(x) x[, j])
  }
  A <- fit$mcmc$samples_by_chain$A
  for (j in seq_len(dim(A[[1L]])[2L]))
    out[[paste0("A[", j, ",1]")]] <- lapply(A, function(x) x[, j, 1L])
  T <- fit$mcmc$samples_by_chain$tau
  for (j in seq_len(ncol(T[[1L]])))
    out[[colnames(T[[1L]])[j]]] <- lapply(T, function(x) x[, j])
  out
}

draw_trace_panel <- function(series, path, title) {
  png(path, width = 2400, height = 1800, res = 180)
  on.exit(dev.off(), add = TRUE)
  old <- par(mfrow = c(4, 4), mar = c(2.4, 2.7, 2.1, 0.7), oma = c(0, 0, 2, 0))
  on.exit(par(old), add = TRUE)
  cols <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
  for (nm in names(series)) {
    chains <- series[[nm]]
    n <- length(chains[[1L]])
    ids <- unique(as.integer(round(seq(1, n, length.out = min(n, 2000L)))))
    ylim <- range(unlist(lapply(chains, function(x) x[ids])), finite = TRUE)
    plot(ids, chains[[1L]][ids], type = "l", col = cols[1L], xlab = "retained draw",
         ylab = "", main = nm, ylim = ylim)
    if (length(chains) > 1L) for (cc in 2:length(chains))
      lines(ids, chains[[cc]][ids], col = cols[cc])
  }
  mtext(title, outer = TRUE, cex = 1.2)
}

draw_convergence_plots <- function(fit, diagnostic, fold_dir, stage) {
  draw_trace_panel(trace_series(fit),
                   file.path(fold_dir, paste0("trace_key_", stage, ".png")),
                   paste("Five-fold CV", stage, "key traces"))

  udiag <- fit$diagnostics$rhat_U
  worst <- head(udiag$global_index[order(udiag$rhat, decreasing = TRUE)], 6L)
  U <- fit$mcmc$samples_by_chain$U
  useries <- setNames(lapply(worst, function(j) lapply(U, function(x) x[, j, 1L])),
                      paste0("U[", worst, "]"))
  draw_trace_panel(useries,
                   file.path(fold_dir, paste0("trace_worst_U_", stage, ".png")),
                   paste("Worst missing-U traces by raw R-hat:", stage))

  tab <- diagnostic$table
  png(file.path(fold_dir, paste0("rhat_overview_", stage, ".png")),
      width = 1800, height = 1100, res = 180)
  cols <- c(hyper = "#0072B2", loading = "#D55E00",
            threshold = "#009E73", missing_U = "#777777")
  plot(seq_len(nrow(tab)), tab$rhat, pch = 16, col = cols[tab$group],
       xlab = "raw parameter (grouped)", ylab = "rank-normalized split R-hat",
       main = paste("Five-fold CV", stage),
       ylim = range(c(1, tab$rhat, R_HAT_LIMIT), finite = TRUE))
  abline(h = R_HAT_LIMIT, col = "red", lty = 2, lwd = 2)
  legend("topright", legend = names(cols), col = cols, pch = 16, bty = "n")
  dev.off()
}

balanced_draw_ids <- function(fit, total = 500L) {
  info <- fit$mcmc$mcmc_draw_info
  chains <- sort(unique(info$chain))
  target <- max(1L, floor(total / length(chains)))
  sort(unique(unlist(lapply(chains, function(cc) {
    a <- which(info$chain == cc)
    a[unique(as.integer(round(seq(1, length(a), length.out = min(target, length(a))))))]
  }), use.names = FALSE)))
}

score_predictions <- function(z) {
  split_z <- c(list(overall = z), split(z, paste0("R", z$R)))
  do.call(rbind, lapply(names(split_z), function(g) {
    x <- split_z[[g]]
    data.frame(
      group = g, n = nrow(x),
      RMSE = sqrt(mean((x$y_true - x$pred_mean)^2)),
      MAE = mean(abs(x$y_true - x$pred_mean)),
      Coverage95 = mean(x$y_true >= x$pred_lo95 & x$y_true <= x$pred_hi95),
      Width95 = mean(x$pred_hi95 - x$pred_lo95),
      stringsAsFactors = FALSE
    )
  }))
}

run_one_fold <- function(fold) {
  fold_dir <- file.path(OUTPUT_ROOT, paste0("fold_", fold))
  dir.create(fold_dir, recursive = TRUE, showWarnings = FALSE)
  unlink(file.path(fold_dir, "CASE_STUDY_COMPLETED.txt"))
  set_status(fold_dir, fold, "prepare", "Preparing frozen train/test data.")

  test_ids <- folds$RID[folds$fold == fold]
  train <- dat[!(dat$RID %in% test_ids), , drop = FALSE]
  test <- dat[dat$RID %in% test_ids, , drop = FALSE]
  train <- train[order(train$RID), , drop = FALSE]
  test <- test[order(test$RID), , drop = FALSE]
  stopifnot(nrow(train) == 396L, nrow(test) == 99L)

  fit_seed <- 202660000L + fold
  prediction_seed <- 202670000L + fold
  case_seed <- 202680000L + fold

  X_train <- as.matrix(train[, X_COLUMNS, drop = FALSE])
  C_train <- make_C(train)
  U_train <- matrix(train[[U_COLUMN]], ncol = 1L,
                    dimnames = list(NULL, U_COLUMN))
  U_train[train$R == 0L, 1L] <- NA_real_

  config <- data.frame(
    run_mode = "five_fold", cv_folds = CV_FOLDS, fold = fold,
    assignment_source = "fixed_R_stratified_5fold", assignment_seed = ADNI_CV_SEED,
    data_md5 = unname(tools::md5sum(DATA_FILE)),
    folds_md5 = unname(tools::md5sum(FOLD_FILE)),
    n_train = nrow(train), n_test = nrow(test),
    n_train_observed_U = sum(train$R == 1L),
    n_train_missing_U = sum(train$R == 0L),
    engine = "multivariate", latent_dim = 1L, ident = "none",
    kernel = "se", standardize_U = TRUE,
    n_chains = N_CHAINS, burn = BURN, initial_iter = INITIAL_ITER,
    extension_iter = EXTENSION_ITER, checkpoint_every = CHECKPOINT_EVERY,
    u_block_size = SAMPLER_CONTROL$u_block_size,
    cross_every = SAMPLER_CONTROL$cross_every,
    gp_block_schur = SAMPLER_CONTROL$gp_block_schur,
    dictionary_mode = SAMPLER_CONTROL$dictionary_mode,
    max_dictionary_size = SAMPLER_CONTROL$max_dictionary_size,
    ess_max_steps = SAMPLER_CONTROL$ess_max_steps,
    threshold_update = SAMPLER_CONTROL$threshold_update,
    fit_seed = fit_seed, prediction_seed = prediction_seed, case_seed = case_seed,
    eivgp_version = as.character(utils::packageVersion("eivGP")),
    tested_commit = EIVGP_TESTED_COMMIT,
    installed_commit = if (is.null(pkg$RemoteSha)) "unrecorded" else pkg$RemoteSha,
    stringsAsFactors = FALSE
  )
  config_path <- file.path(fold_dir, "run_config.csv")
  if (file.exists(file.path(fold_dir, "fit_latest.rds"))) {
    if (!file.exists(config_path)) stop("Checkpoint has no run configuration.")
    old_config <- read.csv(config_path, stringsAsFactors = FALSE)
    if (!isTRUE(all.equal(old_config, config, check.attributes = FALSE)))
      stop("Checkpoint configuration mismatch; use a separate output directory.")
  }
  write.csv(config, config_path, row.names = FALSE)

  fit_path <- file.path(fold_dir, "fit_latest.rds")
  if (file.exists(fit_path)) {
    fit <- readRDS(fit_path)
    stopifnot(
      identical(fit$sampler_version, "0.3.1"),
      identical(fit$control$u_block_size, SAMPLER_CONTROL$u_block_size),
      identical(fit$control$cross_every, SAMPLER_CONTROL$cross_every),
      fit$diagnostics$summary$burn == BURN
    )
    stamp(fold_dir, paste("Loaded checkpoint at iteration", fit$diagnostics$summary$n_iter))
  } else {
    first_target <- min(INITIAL_ITER, BURN + CHECKPOINT_EVERY)
    set_status(fold_dir, fold, "fit", paste("Initial fit to iteration", first_target))
    started <- Sys.time()
    fit <- fit_eivgp(
      X = X_train, y = train[[Y_COLUMN]], C = C_train, U_obs = U_train,
      engine = "multivariate", latent_dim = 1L, ident = "none",
      standardize_U = TRUE, store_scores = FALSE, kernel = "se",
      parallel = TRUE, n_cores = N_CORES, n_chains = N_CHAINS,
      n_iter = first_target, burn = BURN, thin = 1L,
      seed = fit_seed,
      verbose = FALSE, sampler_control = SAMPLER_CONTROL
    )
    atomic_save_rds(fit, fit_path)
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    stamp(fold_dir, sprintf("Saved checkpoint at %d; segment %.1f min.",
                            first_target, elapsed / 60))
  }

  extend_to <- function(fit, target) {
    while (fit$diagnostics$summary$n_iter < target) {
      current <- fit$diagnostics$summary$n_iter
      add <- min(CHECKPOINT_EVERY, target - current)
      set_status(fold_dir, fold, "continue",
                 sprintf("Iteration %d -> %d; checkpoint will overwrite fit_latest.rds.",
                         current, current + add))
      started <- Sys.time()
      fit <- continue_eivgp(fit, n_iter = add, parallel = TRUE,
                            n_cores = N_CORES, verbose = FALSE)
      atomic_save_rds(fit, fit_path)
      elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
      remaining <- target - fit$diagnostics$summary$n_iter
      eta <- if (remaining > 0L) elapsed / add * remaining / 60 else 0
      stamp(fold_dir, sprintf(
        "Saved iteration %d; last %d took %.1f min; ETA to %d about %.1f min.",
        fit$diagnostics$summary$n_iter, add, elapsed / 60, target, eta
      ))
    }
    fit
  }

  fit <- extend_to(fit, INITIAL_ITER)
  current_iter <- fit$diagnostics$summary$n_iter
  final_stage <- paste0("iter", current_iter)
  diagnostic <- light_diagnostics(fit, fold_dir, final_stage)
  draw_convergence_plots(fit, diagnostic, fold_dir, final_stage)
  if (!diagnostic$passed && current_iter < INITIAL_ITER + EXTENSION_ITER) {
    fit <- extend_to(fit, INITIAL_ITER + EXTENSION_ITER)
    final_stage <- paste0("iter", fit$diagnostics$summary$n_iter)
    diagnostic <- light_diagnostics(fit, fold_dir, final_stage)
    draw_convergence_plots(fit, diagnostic, fold_dir, final_stage)
  }
  stamp(fold_dir, paste(final_stage, "convergence gate:", diagnostic$passed))

  set_status(fold_dir, fold, "prediction", paste("Predicting from", final_stage, "fit."))
  X_test <- as.matrix(test[, X_COLUMNS, drop = FALSE])
  C_test <- make_C(test)
  U_test <- matrix(test[[U_COLUMN]], ncol = 1L,
                   dimnames = list(NULL, U_COLUMN))
  U_test[test$R == 0L, 1L] <- NA_real_
  draw_ids <- balanced_draw_ids(fit, PREDICTION_DRAWS)
  response_draws <- predict_eivgp(
    fit, new_X = X_test, new_C = C_test, new_U = U_test,
    new_U_scale = "raw", target = "response", draw_ids = draw_ids,
    seed = prediction_seed
  )
  pred <- data.frame(
    run_mode = "five_fold",
    cv_folds = CV_FOLDS,
    fold = fold, final_stage = final_stage,
    convergence_passed = diagnostic$passed,
    RID = test$RID, diagnosis = test$diagnosis, R = test$R,
    y_true = test[[Y_COLUMN]], pred_mean = colMeans(response_draws),
    pred_sd = apply(response_draws, 2L, stats::sd),
    pred_lo95 = apply(response_draws, 2L, stats::quantile, probs = 0.025),
    pred_hi95 = apply(response_draws, 2L, stats::quantile, probs = 0.975),
    n_prediction_draws = nrow(response_draws), stringsAsFactors = FALSE
  )
  write.csv(pred, file.path(fold_dir, "eivgp_predictions.csv"), row.names = FALSE)
  metrics <- score_predictions(pred)
  metrics$convergence_passed <- diagnostic$passed
  metrics$run_mode <- "five_fold"
  metrics$cv_folds <- CV_FOLDS
  metrics$fold <- fold; metrics$final_stage <- final_stage
  write.csv(metrics, file.path(fold_dir, "eivgp_metrics.csv"), row.names = FALSE)
  atomic_write_lines(c(
    paste0("completed_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("cv_folds=", CV_FOLDS),
    paste0("fit_unit=", fold),
    paste0("convergence_passed=", diagnostic$passed),
    paste0("final_stage=", final_stage)
  ), file.path(fold_dir, "EIV_COMPLETED.txt"))
  set_status(fold_dir, fold, "EIV_complete", paste("EIV prediction completed at", final_stage))
  cat("\n### Fold", fold, "\n\n")
  print(knitr::kable(metrics, digits = 3))
  set_status(fold_dir, fold, "comparison", "Running shared prediction and CSF inference tasks.")
  adni_extended_fold(fit, train, test, draw_ids, fold_dir, CV_FOLDS, fold,
                     SMOKE_TEST, diagnostic$passed, response_draws,
                     seed_base = case_seed)
  atomic_write_lines(c(paste0("completed_at=", Sys.time()),
    paste0("schema=", ADNI_CASE_SCHEMA), paste0("convergence_passed=", diagnostic$passed)),
    file.path(fold_dir, "CASE_STUDY_COMPLETED.txt"))
  set_status(fold_dir, fold, "complete", "EIV fit, comparison, and additional inference finished; inspect method status.")
  rm(fit, response_draws)
  invisible(gc())
  metrics
}

if (!PREPARE_ONLY) {
cat("Running fit units; progress and errors are written to fold_N/worker.log.\n")
all_metrics <- adni_run_folds(FOLDS_TO_RUN, run_one_fold, CORE_PLAN, OUTPUT_ROOT)
combined_metrics <- do.call(rbind, all_metrics)
write.csv(combined_metrics, file.path(OUTPUT_ROOT, "all_fold_metrics.csv"), row.names = FALSE)
}

if (!PREPARE_ONLY) {
  main_table <- adni_case_report(OUTPUT_ROOT, folds_expected = FOLDS_TO_RUN,
                                 smoke = SMOKE_TEST)

}

writeLines(capture.output(sessionInfo()), file.path(OUTPUT_ROOT, "sessionInfo.txt"))
sessionInfo()
if (!PREPARE_ONLY) unlink(RUN_LOCK, recursive = TRUE)
