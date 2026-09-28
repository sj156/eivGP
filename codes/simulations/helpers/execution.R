# Simulation execution: loaded by codes/simulation_helpers.R.

mixedgp_cell_controls_study1 <- function(config, cell, cell_output) {
  mechanism_calib <- if (50L %in% cell$calibration_grid) {
    50L
  } else {
    max(cell$calibration_grid)
  }
  list(
    STUDY1_CONFIG = if (config$mode == "smoke") "quick" else "thorough",
    STUDY1_QUICK = config$mode == "smoke",
    STUDY1_OUT_PREFIX = cell_output,
    STUDY1_DATA_DIR = mixedgp_data_cell_directory(config, cell),
    STUDY1_SCENARIO = cell$scenario,
    STUDY1_HETEROGENEITY_ETA = cell$heterogeneity_eta,
    STUDY1_THRESHOLD_DESIGN = cell$threshold_design,
    STUDY1_CALIB_GRID = cell$calibration_grid,
    STUDY1_MC_N_TRAIN = cell$n,
    STUDY1_MC_N_TEST = cell$n_test,
    STUDY1_MC_N_REP = cell$n_rep,
    STUDY1_MC_M = cell$m,
    STUDY1_MC_N_ITER = config$mcmc$n_iter,
    STUDY1_MC_BURN = config$mcmc$burn,
    STUDY1_MC_N_CHAINS = config$mcmc$n_chains,
    STUDY1_MC_N_PRED_DRAW = config$evaluation$n_pred_draw,
    STUDY1_MC_N_M_EVAL = config$evaluation$n_m_eval,
    STUDY1_MC_N_M_DRAW = config$evaluation$n_m_draw,
    STUDY1_MC_N_M_LATENT = config$evaluation$n_m_latent,
    STUDY1_MEAS_N_ITER = config$measurement_mcmc$n_iter,
    STUDY1_MEAS_BURN = config$measurement_mcmc$burn,
    STUDY1_USE_CACHE = config$use_cache,
    STUDY1_REUSE_LOCKED_EIV = FALSE,
    STUDY1_STRICT_COMPETITORS = config$strict_competitors,
    STUDY1_PUBLISHED_COMPETITORS = config$published_methods,
    STUDY1_RUN_ABLATIONS = cell$run_ablations,
    STUDY1_EVALUATE_F = cell$evaluate_f,
    STUDY1_EVALUATE_U = cell$evaluate_u,
    STUDY1_PARALLEL_LEVEL = config$parallel$level,
    STUDY1_DATASET_WORKERS = config$parallel$workers,
    STUDY1_CHAIN_WORKERS = if (is.null(config$parallel$chain_workers)) 1L else config$parallel$chain_workers,
    STUDY1_REQUIRE_MCMC_GATE = config$mcmc$require_gate,
    STUDY1_MAX_RHAT = config$mcmc$rhat_limit,
    STUDY1_MIN_ESS = config$mcmc$ess_limit,
    STUDY1_MECHANISM_CALIB = mechanism_calib
  )
}

mixedgp_cell_controls_study2 <- function(config, cell, cell_output) {
  primary_grid <- if (cell$scenario == "primary") {
    cell$calibration_grid
  } else {
    sort(unique(c(0L, cell$calibration_grid)))
  }
  contrast_calib <- max(cell$calibration_grid)
  list(
    STUDY2_CONFIG = if (config$mode == "smoke") "quick" else "thorough",
    STUDY2_OUT_PREFIX = cell_output,
    STUDY2_DATA_DIR = mixedgp_data_cell_directory(config, cell),
    STUDY2_SCENARIOS = cell$scenario,
    STUDY2_PRIMARY_CALIB_GRID = primary_grid,
    STUDY2_CONTRAST_CALIB = contrast_calib,
    STUDY2_Q = cell$q,
    STUDY2_M = cell$m,
    STUDY2_MC_N_TRAIN = cell$n,
    STUDY2_MC_N_TEST = cell$n_test,
    STUDY2_MC_N_REP = cell$n_rep,
    STUDY2_MC_N_ITER = config$mcmc$n_iter,
    STUDY2_MC_BURN = config$mcmc$burn,
    STUDY2_MC_THIN = config$mcmc$thin,
    STUDY2_MC_N_CHAINS = config$mcmc$n_chains,
    STUDY2_MC_N_PRED_DRAW = config$evaluation$n_pred_draw,
    STUDY2_MC_N_M_EVAL = config$evaluation$n_m_eval,
    STUDY2_MC_N_M_DRAW = config$evaluation$n_m_draw,
    STUDY2_MC_N_M_LATENT = config$evaluation$n_m_latent,
    STUDY2_MC_N_M_TRUTH = config$evaluation$n_m_truth,
    STUDY2_MC_N_ORACLE_POOL = config$evaluation$n_oracle_pool,
    STUDY2_MEAS_N_ITER = config$measurement_mcmc$n_iter,
    STUDY2_MEAS_BURN = config$measurement_mcmc$burn,
    STUDY2_MEAS_THIN = config$measurement_mcmc$thin,
    STUDY2_MEAS_N_CHAINS = config$measurement_mcmc$n_chains,
    STUDY2_ABLATION_GP_N_STARTS = config$ablation_gp$n_starts,
    STUDY2_ABLATION_GP_MAXIT = config$ablation_gp$maxit,
    STUDY2_USE_CACHE = config$use_cache,
    STUDY2_MC_RESUME = TRUE,
    STUDY2_SAVE_REP_FITS = FALSE,
    STUDY2_STRICT_COMPETITORS = config$strict_competitors,
    STUDY2_PUBLISHED_COMPETITORS = config$published_methods,
    STUDY2_RUN_ABLATIONS = cell$run_ablations,
    STUDY2_ABLATION_SCENARIOS = if (cell$run_ablations) cell$scenario else character(0),
    STUDY2_EVALUATE_F = cell$evaluate_f,
    STUDY2_EVALUATE_U = cell$evaluate_u,
    STUDY2_PARALLEL_LEVEL = config$parallel$level,
    STUDY2_DATASET_WORKERS = config$parallel$workers,
    STUDY2_CHAIN_WORKERS = if (is.null(config$parallel$chain_workers)) 1L else config$parallel$chain_workers,
    STUDY2_PREDICTIVE_LATENT_SAMPLER = config$predictive_latent_sampler,
    STUDY2_ENFORCE_MCMC_GATE = config$mcmc$require_gate,
    STUDY2_MCMC_RHAT_LIMIT = config$mcmc$rhat_limit,
    STUDY2_MCMC_RAW_ESS_LIMIT = config$mcmc$raw_ess_limit,
    STUDY2_MCMC_TARGET_BULK_ESS_LIMIT = config$mcmc$target_bulk_ess_limit,
    STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT = config$mcmc$target_tail_ess_limit,
    ## Smoke mode is explicitly non-reportable: permissive measurement-chain
    ## thresholds let PI-GP/CC-GP execute so every code path is tested. The
    ## publication constructor retains the prespecified strict diagnostics.
    STUDY2_MEAS_RHAT_LIMIT = if (config$mode == "smoke") {
      3
    } else {
      config$mcmc$rhat_limit
    },
    STUDY2_MEAS_BULK_ESS_LIMIT = if (config$mode == "smoke") {
      1L
    } else {
      config$mcmc$target_bulk_ess_limit
    },
    STUDY2_MEAS_TAIL_ESS_LIMIT = if (config$mode == "smoke") {
      1L
    } else {
      config$mcmc$target_tail_ess_limit
    },
    STUDY2_SAVE_PDF = config$mode %in% c("publication", "development"),
    STUDY2_SAVE_PNG = FALSE
  )
}

mixedgp_get_cell_output <- function(envir, name) {
  get0(name, envir = envir, inherits = FALSE, ifnotfound = data.frame())
}

## Catch errors inside workers before mclapply converts them to try-error
## strings. Keep successful replications and expose the original failure.
mixedgp_save_replication <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(basename(path), ".tmp-"), dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) stop("Could not commit replication: ", path)
  invisible(path)
}

mixedgp_run_replications <- function(X, FUN, n_cores, seeds, status_path,
                                    study, parallel_map, mc.preschedule = FALSE) {
  ## Each worker owns one small status file. No shared CSV writes from forks.
  task_dir <- paste0(status_path, ".tasks")
  dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
  task_paths <- file.path(task_dir, sprintf("task-%05d.csv", seq_along(X)))
  record <- function(i, state, message = "", elapsed = NA_real_) {
    mixedgp_atomic_write_csv(data.frame(
      task = paste(X[[i]], collapse = ","), status = state,
      message = message, elapsed_seconds = elapsed,
      updated_at = format(Sys.time(), tz = "UTC", usetz = TRUE)), task_paths[i])
  }
  for (i in seq_along(X)) if (!file.exists(task_paths[i])) record(i, "pending")
  worker <- function(i) {
    record(i, "running")
    started <- proc.time()[["elapsed"]]
    result <- tryCatch(list(ok = TRUE, value = FUN(X[[i]]), message = "", call = "",
                  elapsed = proc.time()[["elapsed"]] - started),
      error = function(e) list(ok = FALSE, value = NULL,
        message = conditionMessage(e),
        call = paste(deparse(conditionCall(e)), collapse = " "),
        elapsed = proc.time()[["elapsed"]] - started))
    record(i, if (result$ok) "success" else "failed", result$message, result$elapsed)
    result
  }
  results <- parallel_map(as.list(seq_along(X)), worker, n_cores = n_cores,
    seeds = seeds, mc.preschedule = mc.preschedule)
  results <- lapply(results, function(x) {
    if (is.list(x) && is.logical(x$ok) && length(x$ok) == 1L && !is.na(x$ok)) return(x)
    list(ok = FALSE, value = NULL, message = if (inherits(x, "try-error"))
      as.character(x) else "Worker returned no valid result (possibly terminated).",
      call = "", elapsed = NA_real_)
  })
  ok <- vapply(results, function(x) isTRUE(x$ok), logical(1))
  status <- data.frame(task = vapply(X, function(x) paste(x, collapse = ","), ""),
    status = ifelse(ok, "success", "failed"),
    message = vapply(results, `[[`, "", "message"),
    call = vapply(results, `[[`, "", "call"),
    elapsed_seconds = vapply(results, `[[`, 0, "elapsed"))
  mixedgp_atomic_write_csv(status, status_path)
  if (any(!ok)) {
    detail <- paste0(study, ": ", sum(!ok), "/", length(ok),
      " replications failed. First error: ", status$message[which(!ok)[1L]],
      ". See ", status_path)
    if (!any(ok)) stop(detail, call. = FALSE)
    warning(detail, ". Aggregating successful replications only; report failures alongside estimates.",
            call. = FALSE)
  }
  lapply(results[ok], `[[`, "value")
}

mixedgp_run_study1_cell <- function(config, cell, engine, run_dir) {
  old_wd <- setwd(config$code_dir)
  on.exit(setwd(old_wd), add = TRUE)
  cell_output <- file.path(run_dir, "cells", cell$id)
  dir.create(cell_output, recursive = TRUE, showWarnings = FALSE)
  run_env <- new.env(parent = engine)
  list2env(
    mixedgp_cell_controls_study1(config, cell, cell_output),
    envir = run_env
  )
  sys.source(
    file.path(config$code_dir, "simulations/02_study1_monte_carlo.R"),
    envir = run_env,
    chdir = FALSE
  )
  outputs <- list(
    predictive_metrics = mixedgp_get_cell_output(run_env, "mc_results"),
    mean_recovery = mixedgp_get_cell_output(run_env, "mean_recovery"),
    latent_imputation = mixedgp_get_cell_output(run_env, "latent_imputation"),
    surface_recovery = mixedgp_get_cell_output(run_env, "surface_recovery"),
    ablation_predictive_metrics = mixedgp_get_cell_output(run_env, "ablation_results"),
    ablation_surface_recovery = mixedgp_get_cell_output(
      run_env, "ablation_surface_recovery"
    ),
    mcmc_diagnostics = mixedgp_get_cell_output(run_env, "mcmc_diagnostics"),
    competitor_status = mixedgp_get_cell_output(run_env, "competitor_status"),
    ablation_status = mixedgp_get_cell_output(run_env, "ablation_status"),
    sampler_control_manifest = mixedgp_get_cell_output(
      run_env, "sampler_control_manifest"
    )
  )
  list(
    cell = cell,
    outputs = outputs,
    design_tag = get0("STUDY1_DESIGN_TAG", envir = run_env, inherits = FALSE),
    output_directory = cell_output
  )
}

mixedgp_run_study2_cell <- function(config, cell, engine, run_dir) {
  old_wd <- setwd(config$code_dir)
  on.exit(setwd(old_wd), add = TRUE)
  cell_output <- file.path(run_dir, "cells", cell$id)
  dir.create(cell_output, recursive = TRUE, showWarnings = FALSE)
  run_env <- new.env(parent = engine)
  list2env(
    mixedgp_cell_controls_study2(config, cell, cell_output),
    envir = run_env
  )
  sys.source(
    file.path(config$code_dir, "simulations/02_study2_monte_carlo.R"),
    envir = run_env,
    chdir = FALSE
  )
  outputs <- get0(
    "raw_outputs", envir = run_env, inherits = FALSE, ifnotfound = list()
  )
  list(
    cell = cell,
    outputs = outputs,
    design_tag = get0("STUDY2_DESIGN_TAG", envir = run_env, inherits = FALSE),
    output_directory = cell_output
  )
}

