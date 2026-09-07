## Fitting/evaluation only. Reporting lives in report_study2_results.R.
sys.source("setup_study2_experiment.R", envir = environment(), chdir = TRUE)

run_one_study2_replication <- function(rep_id, scenario) {
  calib_grid <- scenario_calib_grid(scenario)
  scenario_id <- match(scenario, allowed_scenarios)
  ## The three Gaussian scenarios use common random numbers. For a fixed
  ## replication, they share X, U, Gaussian score innovations, and response
  ## innovations. The logistic robustness scenario shares X and U; its score
  ## errors necessarily use a different random-number transformation.
  data_seed_base <- 1000000L + 10000L * rep_id
  fit_seed_base <- 1000000L * scenario_id + 10000L * rep_id

  data_path <- file.path(
    STUDY2_DATA_DIR,
    mixedgp_dataset_filename(
      "study2", rep_id, scenario, n_train, n_test, m_vec[1L], q = length(m_vec)
    )
  )
  if (!file.exists(data_path)) {
    stop(
      "Missing frozen Study II dataset: ", data_path,
      ". Run 12_generate_synthetic_datasets.R first."
    )
  }
  frozen <- load_mixedgp_synthetic_dataset_strict(
    data_path,
    expected = list(
      study = "study2", scenario = scenario, rep_id = rep_id,
      n = n_train, n_test = n_test, q = length(m_vec), m = m_vec[1L],
      calib_grid = calib_grid
    )
  )
  dat <- frozen$data
  train <- dat$train
  test <- dat$test

  pattern_info <- classify_study2_pattern_frequency(train$C, test$C)
  pattern_counts <- as.data.frame(table(pattern_info$pattern_stratum))
  names(pattern_counts) <- c("evaluation_stratum", "n_test")
  pattern_counts$rep <- rep_id
  pattern_counts$scenario <- scenario

  calib_sets <- frozen$calibration_sets[as.character(calib_grid)]
  if (any(vapply(calib_sets, is.null, logical(1)))) {
    stop("Frozen Study II dataset does not contain every required calibration split.")
  }

  oracle_pool <- make_oracle_pool_2d(
    true_params = dat$true_params,
    n_pool = n_oracle_pool,
    seed = fit_seed_base + 3L
  )
  oracle_draws <- sample_oracle_test_y_2d(
    X_test = test$X,
    C_test = test$C,
    true_params = dat$true_params,
    sigma_eps = dat$sigma_eps,
    n_draw = n_pred_draw,
    oracle_pool = oracle_pool,
    seed = fit_seed_base + 4L
  )

  competitor_result <- mixedgp_cached_competitors(
    X_train = train$X,
    y_train = train$y,
    C_train = train$C,
    X_test = test$X,
    C_test = test$C,
    n_draw = n_pred_draw,
    seed = fit_seed_base + 100L,
    m_vec = m_vec,
    methods = STUDY2_PUBLISHED_COMPETITORS,
    strict = STUDY2_STRICT_COMPETITORS,
    controls = study2_competitor_controls
  )
  competitor_status <- competitor_result$status
  if (nrow(competitor_status) > 0L) {
    competitor_status$rep <- rep_id
    competitor_status$scenario <- scenario
  } else {
    competitor_status$rep <- integer(0)
    competitor_status$scenario <- character(0)
  }

  set.seed(data_seed_base + 150L)
  mean_eval_idx <- sort(sample(
    seq_len(nrow(test$X)), min(n_m_eval, nrow(test$X))
  ))
  target_gate_local_idx <- mixedgp_study2_diagnostic_rows(
    test$X[mean_eval_idx, , drop = FALSE], test$C[mean_eval_idx, , drop = FALSE],
    STUDY2_MCMC_TARGET_N_POINTS
  )
  target_gate_test_rows <- mean_eval_idx[target_gate_local_idx]
  mean_truth <- oracle_m0_2d(
    X = test$X[mean_eval_idx, , drop = FALSE],
    C = test$C[mean_eval_idx, , drop = FALSE],
    true_params = dat$true_params,
    oracle_pool = oracle_pool,
    n_latent = n_m_truth,
    seed = data_seed_base + 151L
  )
  mean_truth_rejection <- attr(mean_truth, "rejection_telemetry")
  mean_truth_diagnostics <- attr(mean_truth, "truth_diagnostics")
  mean_truth_rejection$rep <- rep_id
  mean_truth_rejection$scenario <- scenario
  mean_truth_diagnostics$rep <- rep_id
  mean_truth_diagnostics$scenario <- scenario
  mean_recovery <- list()

  metrics <- list(
    Oracle = summarize_predictive_samples_by_pattern(
      draw_mat = oracle_draws,
      y_true = test$y,
      pattern_stratum = pattern_info$pattern_stratum,
      method = "Oracle",
      rep_id = rep_id,
      n_calib = NA_integer_,
      scenario = scenario
    )
  )

  for (method in names(competitor_result$draws)) {
    metrics[[method]] <- summarize_predictive_samples_by_pattern(
      draw_mat = competitor_result$draws[[method]],
      y_true = test$y,
      pattern_stratum = pattern_info$pattern_stratum,
      method = method,
      rep_id = rep_id,
      n_calib = NA_integer_,
      scenario = scenario
    )
    mean_recovery[[method]] <- summarize_mean_recovery_2d(
      matrix(
        competitor_result$latent_means[[method]][mean_eval_idx],
        nrow = 1L
      ),
      m_true = mean_truth,
      method = method,
      rep_id = rep_id,
      n_calib = NA_integer_,
      scenario = scenario,
      valid_function_draws = FALSE
    )
  }

  diagnostics <- list()
  target_diagnostics <- list()
  imputation <- list()
  imputation_status <- list()
  ablation_imputation <- list()
  surface <- list()
  ablation_metrics <- list()
  ablation_surface <- list()
  measurement_diagnostics <- list()
  measurement_parameter_diagnostics <- list()
  ablation_optimizer_attempts <- list()
  ablation_status <- list()
  sampler_controls <- list()
  sampler_control_manifest <- list()

  for (n_calib in calib_grid) {
    message(
      "Study II scenario=", scenario,
      ", replication=", rep_id,
      ", |O|=", n_calib
    )

    fit <- fit_eivgp_ordprobit_fb(
      X_raw = train$X,
      y_raw = train$y,
      C_ord = train$C,
      U_obs = train$U,
      calib_idx = calib_sets[[as.character(n_calib)]],
      U_true_eval = train$U,
      d = d_latent,
      m_vec = m_vec,
      ident = ident_method,
      n_iter = mc_n_iter,
      burn = mc_burn,
      thin = mc_thin,
      n_chains = mc_n_chains,
      preset = mc_preset,
      store_scores = FALSE,
      seed = fit_seed_base + 1000L + n_calib,
      parallel_chains = parallel_chains,
      n_cores = chain_cores,
      verbose = FALSE
    )

    draw_ids <- seq_len(dim(fit$mcmc$samples_U)[1])
    if (length(draw_ids) > n_pred_draw) {
      draw_ids <- draw_ids[
        unique(round(seq(1, length(draw_ids), length.out = n_pred_draw)))
      ]
    }
    predictive_seed <- fit_seed_base + 2000L + n_calib
    mean_seed <- fit_seed_base + 2500L + n_calib

    eiv_draws <- sample_eiv_test_y_ordprobit_fb(
      X_test_raw = test$X,
      C_test = test$C,
      fit_obj = fit,
      draw_ids = draw_ids,
      n_per_draw = 1L,
      latent_sampler = predictive_latent_sampler,
      n_new_latent_gibbs = diagnostic_n_new_latent_gibbs,
      rejection_max_batches = rejection_max_batches,
      seed = predictive_seed
    )

    metrics[[paste0("EIV_", n_calib)]] <-
      summarize_predictive_samples_by_pattern(
        draw_mat = eiv_draws,
        y_true = test$y,
        pattern_stratum = pattern_info$pattern_stratum,
        method = "EIV-GP",
        rep_id = rep_id,
        n_calib = n_calib,
        scenario = scenario
      )

    m_draw_ids <- draw_ids
    if (length(m_draw_ids) > n_m_draw) {
      m_draw_ids <- m_draw_ids[
        unique(round(seq(1, length(m_draw_ids), length.out = n_m_draw)))
      ]
    }
    m_draws <- sample_eiv_m_given_xc_fb(
      X_test_raw = test$X[mean_eval_idx, , drop = FALSE],
      C_test = test$C[mean_eval_idx, , drop = FALSE],
      fit_obj = fit,
      draw_ids = m_draw_ids,
      n_latent = n_m_latent,
      include_process_uncertainty = TRUE,
      joint = FALSE,
      return_components = TRUE,
      latent_sampler = predictive_latent_sampler,
      n_new_latent_gibbs = diagnostic_n_new_latent_gibbs,
      rejection_max_batches = rejection_max_batches,
      seed = mean_seed
    )
    control_key <- as.character(n_calib)
    sampler_controls[[control_key]] <- study2_sampler_control_settings(
      fit = fit,
      prediction_seed = predictive_seed,
      mean_seed = mean_seed,
      prediction_sampler_used = attr(eiv_draws, "latent_sampler"),
      mean_sampler_used = attr(m_draws, "latent_sampler")
    )
    sampler_control_manifest[[control_key]] <-
      study2_sampler_control_rows(
        fit = fit,
        rep_id = rep_id,
        scenario = scenario,
        n_calib = n_calib,
        prediction_seed = predictive_seed,
        mean_seed = mean_seed,
        prediction_sampler_used = attr(eiv_draws, "latent_sampler"),
        mean_sampler_used = attr(m_draws, "latent_sampler")
      )
    mean_recovery[[paste0("EIV_", n_calib)]] <-
      summarize_mean_recovery_2d(
        m_draws,
        m_true = mean_truth,
        method = "EIV-GP",
        rep_id = rep_id,
        n_calib = n_calib,
        scenario = scenario,
        valid_function_draws = TRUE
      )

    target_series <- mixedgp_study2_target_series(
      fit, test$X[target_gate_test_rows, , drop = FALSE],
      test$C[target_gate_test_rows, , drop = FALSE],
      U = if (n_calib > 0L && isTRUE(STUDY2_EVALUATE_F)) {
        test$U[target_gate_test_rows, , drop = FALSE]
      } else NULL,
      max_draws_per_chain = STUDY2_MCMC_TARGET_DRAWS_PER_CHAIN,
      n_latent = STUDY2_MCMC_TARGET_N_LATENT,
      seed = fit_seed_base + 2900L
    )
    invariant_series <- mixedgp_study2_invariant_series(fit)
    raw_series <- mixedgp_study2_raw_series(fit)
    summarize_gate <- function(series, ess_min, tail_min) {
      mixedgp_summarize_diagnostic_series(
        series, STUDY2_MCMC_RHAT_LIMIT, ess_min, tail_min
      )
    }
    functional_diag <- summarize_gate(target_series, STUDY2_MCMC_TARGET_BULK_ESS_LIMIT,
                                       STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT)
    invariant_diag <- summarize_gate(invariant_series, STUDY2_MCMC_TARGET_BULK_ESS_LIMIT,
                                     STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT)
    raw_diag <- summarize_gate(raw_series, STUDY2_MCMC_RAW_ESS_LIMIT,
                               STUDY2_MCMC_RAW_ESS_LIMIT)
    functional_diag$target <- sub("\\[.*$", "", functional_diag$parameter)
    invariant_diag$target <- "ordinal_measurement_invariant"
    raw_diag$target <- "raw_coordinate"
    functional_diag$used_for_gate <- TRUE
    invariant_diag$used_for_gate <- TRUE
    raw_diag$used_for_gate <- n_calib > 0L
    target_diag <- bind_rows(functional_diag, invariant_diag, raw_diag)
    target_diag$rep <- rep_id
    target_diag$scenario <- scenario
    target_diag$n_calib <- n_calib
    target_diagnostics[[as.character(n_calib)]] <- target_diag

    raw_coordinate_pass <- mixedgp_diagnostic_table_pass(
      raw_diag, STUDY2_MCMC_RHAT_LIMIT, STUDY2_MCMC_RAW_ESS_LIMIT,
      STUDY2_MCMC_RAW_ESS_LIMIT, expected_parameters = names(raw_series)
    )
    target_functional_pass <- mixedgp_diagnostic_table_pass(
      functional_diag, STUDY2_MCMC_RHAT_LIMIT, STUDY2_MCMC_TARGET_BULK_ESS_LIMIT,
      STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT, expected_parameters = names(target_series)
    )
    invariant_measurement_pass <- mixedgp_diagnostic_table_pass(
      invariant_diag, STUDY2_MCMC_RHAT_LIMIT, STUDY2_MCMC_TARGET_BULK_ESS_LIMIT,
      STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT, expected_parameters = names(invariant_series)
    )
    diag_row <- extract_study2_diagnostics(fit, rep_id, scenario, n_calib)
    if (nrow(diag_row) != 1L) stop("Fit did not return one diagnostic summary.")
    predictive_diag <- functional_diag[
      functional_diag$target == "predictive_second_moment", , drop = FALSE
    ]
    mean_diag <- functional_diag[
      functional_diag$target == "m_conditional_mean", , drop = FALSE
    ]
    diag_row$mcmc_gate_rule <- if (n_calib == 0L) {
      "targets_and_measurement_invariants_without_calibration"
    } else {
      "targets_measurement_invariants_and_free_coordinates_with_calibration"
    }
    diag_row$mcmc_rhat_limit <- STUDY2_MCMC_RHAT_LIMIT
    diag_row$mcmc_raw_ess_limit <- STUDY2_MCMC_RAW_ESS_LIMIT
    diag_row$mcmc_target_bulk_ess_limit <- STUDY2_MCMC_TARGET_BULK_ESS_LIMIT
    diag_row$mcmc_target_tail_ess_limit <- STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT
    diag_row$mcmc_target_n_points <- length(target_gate_test_rows)
    diag_row$mcmc_target_draws_per_chain <- STUDY2_MCMC_TARGET_DRAWS_PER_CHAIN
    diag_row$mcmc_target_integration_draws <- STUDY2_MCMC_TARGET_N_LATENT
    diag_row$raw_coordinate_pass <- raw_coordinate_pass
    diag_row$target_functional_pass <- target_functional_pass
    diag_row$invariant_measurement_pass <- invariant_measurement_pass
    diag_row$predictive_max_rhat <- study2_finite_max(predictive_diag$rhat)
    diag_row$predictive_min_bulk_ess <- study2_finite_min(predictive_diag$ess_bulk)
    diag_row$predictive_min_tail_ess <- study2_finite_min(predictive_diag$ess_tail)
    diag_row$m_max_rhat <- study2_finite_max(mean_diag$rhat)
    diag_row$m_min_bulk_ess <- study2_finite_min(mean_diag$ess_bulk)
    diag_row$m_min_tail_ess <- study2_finite_min(mean_diag$ess_tail)
    diag_row$invariant_max_rhat <- study2_finite_max(invariant_diag$rhat)
    diag_row$invariant_min_bulk_ess <- study2_finite_min(invariant_diag$ess_bulk)
    diag_row$invariant_min_tail_ess <- study2_finite_min(invariant_diag$ess_tail)
    diag_row$mcmc_pass <- target_functional_pass && invariant_measurement_pass &&
      (n_calib == 0L || raw_coordinate_pass)
    diag_row$diagnostic_warning <- !isTRUE(diag_row$mcmc_pass)
    diag_row$diagnostic_advice <- if (isTRUE(diag_row$mcmc_pass)) "" else mixedgp_simulation_diagnostic_advice()
    diagnostics[[as.character(n_calib)]] <- diag_row
    if (!isTRUE(diag_row$mcmc_pass)) {
      saveRDS(list(fit = fit, diagnostics = diag_row),
        file.path(RES_DIR, sprintf("flagged_fit_%s_rep%03d_cal%03d.rds",
                                  scenario, rep_id, n_calib)))
      warning("Study II MCMC diagnostics flagged replication ", rep_id,
        ", calibration ", n_calib, ". ", mixedgp_simulation_diagnostic_advice(),
        call. = FALSE)
    }

    if (isTRUE(STUDY2_EVALUATE_U)) {
      imputation[[paste0("EIV_training_", n_calib)]] <-
        study2_latent_imputation_metrics(
          fit = fit,
          U_true = train$U,
          rep_id = rep_id,
          n_calib = n_calib,
          scenario = scenario
        )
      imputation_status[[paste0("EIV_training_", n_calib)]] <-
        study2_latent_imputation_status(
          fit = fit,
          method = "EIV-GP",
          target = "training_missing_U",
          rep_id = rep_id,
          n_calib = n_calib,
          scenario = scenario,
          n_units = length(fit$data$miss_idx)
        )
      imputation[[paste0("EIV_prospective_", n_calib)]] <-
        study2_eiv_prospective_imputation_metrics(
          fit = fit,
          C_new = test$C[mean_eval_idx, , drop = FALSE],
          U_true = test$U[mean_eval_idx, , drop = FALSE],
          rep_id = rep_id,
          n_calib = n_calib,
          scenario = scenario,
          max_draw = n_m_draw,
          latent_sampler = predictive_latent_sampler,
          n_new_latent_gibbs = diagnostic_n_new_latent_gibbs,
          rejection_max_batches = rejection_max_batches,
          seed = fit_seed_base + 2700L + n_calib
        )
      imputation_status[[paste0("EIV_prospective_", n_calib)]] <-
        study2_latent_imputation_status(
          fit = fit,
          method = "EIV-GP",
          target = "prospective_U_given_C",
          rep_id = rep_id,
          n_calib = n_calib,
          scenario = scenario,
          n_units = length(mean_eval_idx)
        )
    }

    if (isTRUE(STUDY2_EVALUATE_F)) {
      surface[[as.character(n_calib)]] <-
        study2_eiv_surface_recovery_metrics(
          fit = fit,
          scenario = scenario,
          rep_id = rep_id,
          n_calib = n_calib,
          grid_n = if (STUDY2_CONFIG == "quick") 15L else 31L,
          max_draw = if (STUDY2_CONFIG == "quick") 40L else 200L
        )
    }

    if (isTRUE(STUDY2_SAVE_REP_FITS)) {
      saveRDS(
        fit,
        file.path(
          FIT_DIR,
          sprintf("fit_%s_rep%03d_calib%03d.rds", scenario, rep_id, n_calib)
        )
      )
    }

    rm(fit, eiv_draws, m_draws)
    invisible(gc())
  }

  ##########################################################
  ## Appendix-only ablations (kept out of the main table)
  ##########################################################

  if (isTRUE(STUDY2_RUN_ABLATIONS) &&
      scenario %in% STUDY2_ABLATION_SCENARIOS) {
    if (isTRUE(STUDY2_EVALUATE_F)) {
      start_full_u <- proc.time()[3]
      full_u_gp <- tryCatch(
        fit_study2_latent_gp(
          train$X, train$y, train$U,
          gp_n_starts = STUDY2_ABLATION_GP_N_STARTS,
          gp_seed = fit_seed_base + 8000L,
          gp_maxit = STUDY2_ABLATION_GP_MAXIT
        ),
        error = function(e) e
      )
      elapsed_full_u <- proc.time()[3] - start_full_u

      if (inherits(full_u_gp, "error")) {
        ablation_optimizer_attempts[["Full-U GP"]] <-
          study2_optimizer_attempt_rows(
            full_u_gp$optimizer_attempts,
            "Full-U GP",
            NA_integer_
          )
        ablation_status[["Full-U GP"]] <- study2_ablation_status(
          "Full-U GP", NA_integer_, "failed",
          conditionMessage(full_u_gp), elapsed_full_u
        )
      } else {
        ablation_optimizer_attempts[["Full-U GP"]] <-
          study2_optimizer_attempt_rows(
            full_u_gp$optimizer_attempts,
            "Full-U GP",
            NA_integer_,
            full_u_gp$selected_start
          )
        ablation_status[["Full-U GP"]] <- study2_ablation_status(
          "Full-U GP", NA_integer_, "success", "", elapsed_full_u
        )
        ablation_surface[["Full-U GP"]] <-
          study2_latent_gp_surface_metrics(
            gp_fit = full_u_gp,
            scenario = scenario,
            method = "Full-U GP",
            rep_id = rep_id,
            n_calib = NA_integer_,
            grid_n = if (STUDY2_CONFIG == "quick") 15L else 31L
          )
      }
    }

    for (n_calib in calib_grid) {
      key <- as.character(n_calib)
      calib_idx <- calib_sets[[key]]

      start_measurement <- proc.time()[3]
      measurement_fit <- tryCatch(
        fit_ordinalprobit_measurement_fb(
          C_ord = train$C,
          U_obs = train$U,
          calib_idx = calib_idx,
          d = d_latent,
          m_vec = m_vec,
          ident = ident_method,
          n_iter = measurement_n_iter,
          burn = measurement_burn,
          thin = measurement_thin,
          n_chains = measurement_n_chains,
          seed = fit_seed_base + 5000L + n_calib,
          parallel_chains = parallel_chains,
          n_cores = chain_cores,
          rhat_limit = STUDY2_MEAS_RHAT_LIMIT,
          bulk_ess_limit = STUDY2_MEAS_BULK_ESS_LIMIT,
          tail_ess_limit = STUDY2_MEAS_TAIL_ESS_LIMIT,
          verbose = FALSE
        ),
        error = function(e) e
      )
      elapsed_measurement <- proc.time()[3] - start_measurement

      if (inherits(measurement_fit, "error")) {
        msg <- conditionMessage(measurement_fit)
        measurement_status_stub <- list(data = list(
          U_obs = train$U,
          calib_idx = calib_idx,
          miss_idx = setdiff(seq_len(nrow(train$U)), calib_idx),
          d = d_latent
        ))
        for (task_spec in list(
          list(
            target = "training_missing_U",
            n_units = length(measurement_status_stub$data$miss_idx)
          ),
          list(
            target = "prospective_U_given_C",
            n_units = length(mean_eval_idx)
          )
        )) {
          task_status <- study2_latent_imputation_status(
            fit = measurement_status_stub,
            method = "Ordinal model (no Y)",
            target = task_spec$target,
            rep_id = rep_id,
            n_calib = n_calib,
            scenario = scenario,
            n_units = task_spec$n_units
          )
          if (isTRUE(task_status$task_eligible)) {
            task_status$status <- "failed"
            task_status$task_eligible <- FALSE
            task_status$reason <- msg
          }
          imputation_status[[paste0(
            "measurement_", task_spec$target, "_", key
          )]] <- task_status
        }
        ablation_status[[paste0("measurement_", key)]] <-
          study2_ablation_status(
            "Response-free measurement model", n_calib, "failed",
            msg, elapsed_measurement
          )
        ablation_status[[paste0("PI_", key)]] <-
          study2_ablation_status(
            "PI-GP", n_calib, "failed", msg, NA_real_
          )
        ablation_status[[paste0("CC_", key)]] <-
          study2_ablation_status(
            "CC-GP", n_calib,
            if (n_calib == 0L) "not_applicable" else "failed",
            if (n_calib == 0L) "No complete cases." else msg,
            NA_real_
          )
        next
      }

      measurement_diag <- measurement_fit$diagnostics
      measurement_diag$rep <- rep_id
      measurement_diag$scenario <- scenario
      measurement_diag$n_calib <- n_calib
      measurement_diagnostics[[key]] <- measurement_diag

      key_parameter_diag <- measurement_fit$diagnostic_parameters$key
      latent_parameter_diag <- measurement_fit$diagnostic_parameters$latent_rhat
      if (nrow(latent_parameter_diag) > 0L) {
        latent_parameter_diag$block <- "latent"
      }
      parameter_diag <- bind_rows(
        key_parameter_diag[, c(
          "parameter", "block", "rhat", "ess_bulk", "ess_tail"
        )],
        latent_parameter_diag[, intersect(
          c("parameter", "block", "rhat", "ess_bulk", "ess_tail"),
          names(latent_parameter_diag)
        )]
      )
      parameter_diag$rep <- rep_id
      parameter_diag$scenario <- scenario
      parameter_diag$n_calib <- n_calib
      measurement_parameter_diagnostics[[key]] <- parameter_diag

      measurement_training_status <- study2_latent_imputation_status(
        fit = measurement_fit,
        method = "Ordinal model (no Y)",
        target = "training_missing_U",
        rep_id = rep_id,
        n_calib = n_calib,
        scenario = scenario,
        n_units = length(measurement_fit$data$miss_idx)
      )
      measurement_prospective_status <- study2_latent_imputation_status(
        fit = measurement_fit,
        method = "Ordinal model (no Y)",
        target = "prospective_U_given_C",
        rep_id = rep_id,
        n_calib = n_calib,
        scenario = scenario,
        n_units = length(mean_eval_idx)
      )

      measurement_warning <- !isTRUE(measurement_diag$convergence_pass)
      measurement_advice <- if (measurement_warning) mixedgp_simulation_diagnostic_advice() else ""
      if (measurement_warning) {
        saveRDS(list(fit = measurement_fit, diagnostics = measurement_diag),
          file.path(RES_DIR, sprintf("flagged_measurement_%s_rep%03d_cal%03d.rds",
                                    scenario, rep_id, n_calib)))
        warning("Response-free measurement model has diagnostic warnings. ",
                measurement_advice, call. = FALSE)
      }
      measurement_training_status$diagnostic_warning <- measurement_warning
      measurement_training_status$diagnostic_advice <- measurement_advice
      measurement_prospective_status$diagnostic_warning <- measurement_warning
      measurement_prospective_status$diagnostic_advice <- measurement_advice
      imputation_status[[paste0("measurement_training_", key)]] <- measurement_training_status
      imputation_status[[paste0("measurement_prospective_", key)]] <- measurement_prospective_status
      ablation_status[[paste0("measurement_", key)]] <-
        study2_ablation_status("Response-free measurement model", n_calib,
          if (measurement_warning) "completed_with_diagnostic_warnings" else "success",
          measurement_advice, elapsed_measurement)

      ablation_imputation[[paste0("measurement_training_", key)]] <-
        study2_measurement_training_imputation_metrics(
          measurement_fit = measurement_fit,
          U_true = train$U,
          rep_id = rep_id,
          n_calib = n_calib,
          scenario = scenario
        )
      ablation_imputation[[paste0("measurement_prospective_", key)]] <-
        study2_measurement_prospective_imputation_metrics(
          measurement_fit = measurement_fit,
          C_new = test$C[mean_eval_idx, , drop = FALSE],
          U_true = test$U[mean_eval_idx, , drop = FALSE],
          rep_id = rep_id,
          n_calib = n_calib,
          scenario = scenario,
          max_draw = n_m_draw,
          latent_sampler = predictive_latent_sampler,
          n_gibbs = diagnostic_n_new_latent_gibbs,
          rejection_max_batches = rejection_max_batches,
          seed = fit_seed_base + 5600L + n_calib
        )

      start_pi <- proc.time()[3]
      pi_fit <- tryCatch(
        fit_study2_pi_gp(
          train$X,
          train$y,
          measurement_fit,
          gp_n_starts = STUDY2_ABLATION_GP_N_STARTS,
          gp_seed = fit_seed_base + 6100L + n_calib,
          gp_maxit = STUDY2_ABLATION_GP_MAXIT
        ),
        error = function(e) e
      )
      pi_result <- if (inherits(pi_fit, "error")) {
        pi_fit
      } else {
        ablation_optimizer_attempts[[paste0("PI_", key)]] <-
          study2_optimizer_attempt_rows(
            pi_fit$gp$optimizer_attempts,
            "PI-GP",
            n_calib,
            pi_fit$gp$selected_start
          )
        tryCatch(
          list(
            fit = pi_fit,
            draws = sample_study2_pi_gp(
              fit = pi_fit,
              X_test = test$X,
              C_test = test$C,
              n_draw = n_pred_draw,
              n_measurement_draw = min(200L, n_pred_draw),
              latent_sampler = predictive_latent_sampler,
              n_gibbs = diagnostic_n_new_latent_gibbs,
              rejection_max_batches = rejection_max_batches,
              seed = fit_seed_base + 6000L + n_calib
            )
          ),
          error = function(e) e
        )
      }
      elapsed_pi <- proc.time()[3] - start_pi

      if (inherits(pi_result, "error")) {
        if (is.null(ablation_optimizer_attempts[[paste0("PI_", key)]])) {
          ablation_optimizer_attempts[[paste0("PI_", key)]] <-
            study2_optimizer_attempt_rows(
              pi_result$optimizer_attempts,
              "PI-GP",
              n_calib
            )
        }
        ablation_status[[paste0("PI_", key)]] <-
          study2_ablation_status(
            "PI-GP", n_calib, "failed", conditionMessage(pi_result),
            elapsed_pi
          )
      } else {
        ablation_status[[paste0("PI_", key)]] <-
          study2_ablation_status(
            "PI-GP", n_calib,
            if (measurement_warning) "completed_with_diagnostic_warnings" else "success",
            measurement_advice, elapsed_pi
          )
        ablation_metrics[[paste0("PI_", key)]] <-
          summarize_predictive_samples_by_pattern(
            draw_mat = pi_result$draws,
            y_true = test$y,
            pattern_stratum = pattern_info$pattern_stratum,
            method = "PI-GP",
            rep_id = rep_id,
            n_calib = n_calib,
            scenario = scenario
          )
        if (isTRUE(STUDY2_EVALUATE_F) && n_calib > 0L) {
          ablation_surface[[paste0("PI_", key)]] <-
            study2_latent_gp_surface_metrics(
              gp_fit = pi_result$fit$gp,
              scenario = scenario,
              method = "PI-GP",
              rep_id = rep_id,
              n_calib = n_calib,
              grid_n = if (STUDY2_CONFIG == "quick") 15L else 31L
            )
        }
      }

      if (n_calib == 0L) {
        ablation_status[[paste0("CC_", key)]] <-
          study2_ablation_status(
            "CC-GP", n_calib, "not_applicable", "No complete cases.",
            0
          )
      } else {
        start_cc <- proc.time()[3]
        cc_fit <- tryCatch(
          fit_study2_cc_gp(
            X = train$X,
            y = train$y,
            U = train$U,
            calib_idx = calib_idx,
            measurement_fit = measurement_fit,
            gp_n_starts = STUDY2_ABLATION_GP_N_STARTS,
            gp_seed = fit_seed_base + 7100L + n_calib,
            gp_maxit = STUDY2_ABLATION_GP_MAXIT
          ),
          error = function(e) e
        )
        cc_result <- if (inherits(cc_fit, "error")) {
          cc_fit
        } else {
          ablation_optimizer_attempts[[paste0("CC_", key)]] <-
            study2_optimizer_attempt_rows(
              cc_fit$gp$optimizer_attempts,
              "CC-GP",
              n_calib,
              cc_fit$gp$selected_start
            )
          tryCatch(
            list(
              fit = cc_fit,
              draws = sample_study2_cc_gp(
                fit = cc_fit,
                X_test = test$X,
                C_test = test$C,
                n_draw = n_pred_draw,
                latent_sampler = predictive_latent_sampler,
                n_gibbs = diagnostic_n_new_latent_gibbs,
                rejection_max_batches = rejection_max_batches,
                seed = fit_seed_base + 7000L + n_calib
              )
            ),
            error = function(e) e
          )
        }
        elapsed_cc <- proc.time()[3] - start_cc

        if (inherits(cc_result, "error")) {
          if (is.null(ablation_optimizer_attempts[[paste0("CC_", key)]])) {
            ablation_optimizer_attempts[[paste0("CC_", key)]] <-
              study2_optimizer_attempt_rows(
                cc_result$optimizer_attempts,
                "CC-GP",
                n_calib
              )
          }
          ablation_status[[paste0("CC_", key)]] <-
            study2_ablation_status(
              "CC-GP", n_calib, "failed", conditionMessage(cc_result),
              elapsed_cc
            )
        } else {
          ablation_status[[paste0("CC_", key)]] <-
            study2_ablation_status(
              "CC-GP", n_calib,
              if (measurement_warning) "completed_with_diagnostic_warnings" else "success",
              measurement_advice, elapsed_cc
            )
          ablation_metrics[[paste0("CC_", key)]] <-
            summarize_predictive_samples_by_pattern(
              draw_mat = cc_result$draws,
              y_true = test$y,
              pattern_stratum = pattern_info$pattern_stratum,
              method = "CC-GP",
              rep_id = rep_id,
              n_calib = n_calib,
              scenario = scenario
            )
          if (isTRUE(STUDY2_EVALUATE_F)) {
            ablation_surface[[paste0("CC_", key)]] <-
              study2_latent_gp_surface_metrics(
                gp_fit = cc_result$fit$gp,
                scenario = scenario,
                method = "CC-GP",
                rep_id = rep_id,
                n_calib = n_calib,
                grid_n = if (STUDY2_CONFIG == "quick") 15L else 31L
              )
          }
        }
      }

      rm(measurement_fit, pi_result, pi_fit)
      if (exists("cc_result")) rm(cc_result)
      if (exists("cc_fit")) rm(cc_fit)
      invisible(gc())
    }
  }

  ablation_status_df <- bind_rows(ablation_status)
  if (nrow(ablation_status_df) > 0L) {
    ablation_status_df$rep <- rep_id
    ablation_status_df$scenario <- scenario
  } else {
    ablation_status_df$rep <- integer(0)
    ablation_status_df$scenario <- character(0)
  }
  ablation_optimizer_attempts_df <- bind_rows(ablation_optimizer_attempts)
  if (nrow(ablation_optimizer_attempts_df) > 0L) {
    ablation_optimizer_attempts_df$rep <- rep_id
    ablation_optimizer_attempts_df$scenario <- scenario
  } else {
    ablation_optimizer_attempts_df$rep <- integer(0)
    ablation_optimizer_attempts_df$scenario <- character(0)
  }

  list(
    metrics = bind_rows(metrics),
    mean_recovery = bind_rows(mean_recovery),
    mean_truth_rejection = mean_truth_rejection,
    mean_truth_diagnostics = mean_truth_diagnostics,
    diagnostics = bind_rows(diagnostics),
    target_diagnostics = bind_rows(target_diagnostics),
    imputation = bind_rows(c(imputation, ablation_imputation)),
    imputation_status = bind_rows(imputation_status),
    surface = bind_rows(surface),
    ablation_metrics = bind_rows(ablation_metrics),
    ablation_surface = bind_rows(ablation_surface),
    measurement_diagnostics = bind_rows(measurement_diagnostics),
    measurement_parameter_diagnostics = bind_rows(
      measurement_parameter_diagnostics
    ),
    ablation_optimizer_attempts = ablation_optimizer_attempts_df,
    ablation_status = ablation_status_df,
    competitor_status = competitor_status,
    pattern_counts = pattern_counts,
    sampler_control_manifest = bind_rows(sampler_control_manifest),
    metadata = list(
      rep = rep_id,
      scenario = scenario,
      n_train = n_train,
      n_test = n_test,
      data_seed_base = data_seed_base,
      fit_seed_base = fit_seed_base,
      common_random_number_group = rep_id,
      calib_grid = calib_grid,
      true_params = dat$true_params,
      frozen_data_file = basename(data_path),
      frozen_data_md5 = attr(frozen, "manifest_md5"),
      frozen_manifest = attr(frozen, "manifest_path"),
      cache_tag = CACHE_TAG,
      sampler = sampler_controls,
      downstream_monte_carlo = list(
        predictive_draws_requested = n_pred_draw,
        mean_evaluation_points = n_m_eval,
        mean_posterior_draws_requested = n_m_draw,
        mean_latent_integration_draws = n_m_latent,
        mean_truth_latent_draws = n_m_truth,
        prospective_latent_sampler = predictive_latent_sampler,
        diagnostic_gibbs_sweeps = diagnostic_n_new_latent_gibbs,
        rejection_max_batches = rejection_max_batches
      )
    )
  )
}

############################################################
## Resumable execution
############################################################

run_grid <- expand.grid(
  scenario = STUDY2_SCENARIOS,
  rep = seq_len(n_rep),
  stringsAsFactors = FALSE
)

rep_files <- file.path(
  REP_DIR,
  sprintf(
    "study2_%s_rep%03d_%s.rds",
    run_grid$scenario,
    run_grid$rep,
    CACHE_TAG
  )
)

run_or_load_study2_replication <- function(ii) {
  scenario <- run_grid$scenario[ii]
  rep_id <- run_grid$rep[ii]

  if (
    isTRUE(STUDY2_USE_CACHE) &&
      isTRUE(STUDY2_MC_RESUME) &&
      file.exists(rep_files[ii])
  ) {
    cached <- readRDS(rep_files[ii])
    data_path <- file.path(
      STUDY2_DATA_DIR,
      mixedgp_dataset_filename(
        "study2", rep_id, scenario, n_train, n_test, m_vec[1L],
        q = length(m_vec)
      )
    )
    cache_valid <- is.list(cached$metadata) &&
      identical(cached$metadata$rep, as.integer(rep_id)) &&
      identical(cached$metadata$scenario, scenario) &&
      identical(cached$metadata$cache_tag, CACHE_TAG) &&
      identical(cached$metadata$frozen_data_file, basename(data_path)) &&
      file.exists(data_path) &&
      identical(
        tolower(as.character(cached$metadata$frozen_data_md5)),
        tolower(unname(tools::md5sum(data_path)))
      )
    if (isTRUE(cache_valid)) return(cached)
    message("Ignoring incompatible Study II cache: ", rep_files[ii])
  }
  out <- run_one_study2_replication(rep_id, scenario)
  saveRDS(out, rep_files[ii])
  out
}
replication_cores <- if (STUDY2_PARALLEL_LEVEL == "hybrid") {
  min(nrow(run_grid), mixedgp_as_integer_strict(STUDY2_DATASET_WORKERS, "dataset workers", 1L, 1L))
} else if (identical(STUDY2_PARALLEL_LEVEL, "replications")) {
  min(nrow(run_grid), mixedgp_resolve_cores())
} else {
  1L
}
rep_objects <- mixedgp_run_replications(
  as.list(seq_len(nrow(run_grid))),
  run_or_load_study2_replication,
  n_cores = replication_cores,
  seeds = 920000L + seq_len(nrow(run_grid)),
  status_path = file.path(TAB_DIR, "study2_replication_status.csv"),
  study = "Study II",
  parallel_map = mixedgp_parallel_lapply,
  mc.preschedule = FALSE
)

mc_results <- bind_rows(lapply(rep_objects, `[[`, "metrics"))
mc_mean_recovery <- bind_rows(lapply(rep_objects, `[[`, "mean_recovery"))
mc_mean_truth_rejection <- bind_rows(
  lapply(rep_objects, `[[`, "mean_truth_rejection")
)
mc_mean_truth_diagnostics <- bind_rows(
  lapply(rep_objects, `[[`, "mean_truth_diagnostics")
)
mc_diagnostics <- bind_rows(lapply(rep_objects, `[[`, "diagnostics"))
mc_target_diagnostics <- bind_rows(
  lapply(rep_objects, `[[`, "target_diagnostics")
)
mc_imputation <- bind_rows(lapply(rep_objects, `[[`, "imputation"))
mc_imputation_status <- bind_rows(
  lapply(rep_objects, `[[`, "imputation_status")
)
mc_surface <- bind_rows(lapply(rep_objects, `[[`, "surface"))
mc_ablation_results <- bind_rows(lapply(rep_objects, `[[`, "ablation_metrics"))
mc_ablation_surface <- bind_rows(lapply(rep_objects, `[[`, "ablation_surface"))
measurement_diagnostics <- bind_rows(
  lapply(rep_objects, `[[`, "measurement_diagnostics")
)
measurement_parameter_diagnostics <- bind_rows(
  lapply(rep_objects, `[[`, "measurement_parameter_diagnostics")
)
ablation_optimizer_attempts <- bind_rows(
  lapply(rep_objects, `[[`, "ablation_optimizer_attempts")
)
ablation_status <- bind_rows(lapply(rep_objects, `[[`, "ablation_status"))
competitor_status <- bind_rows(lapply(rep_objects, `[[`, "competitor_status"))
pattern_counts <- bind_rows(lapply(rep_objects, `[[`, "pattern_counts"))
sampler_control_manifest <- bind_rows(
  lapply(rep_objects, `[[`, "sampler_control_manifest")
)
replication_metadata <- lapply(rep_objects, `[[`, "metadata")

raw_outputs <- list(
  predictive_metrics = mc_results,
  mean_recovery = mc_mean_recovery,
  mean_truth_rejection = mc_mean_truth_rejection,
  mean_truth_diagnostics = mc_mean_truth_diagnostics,
  mcmc_diagnostics = mc_diagnostics,
  mcmc_target_diagnostics = mc_target_diagnostics,
  latent_imputation = mc_imputation,
  latent_imputation_status = mc_imputation_status,
  surface_recovery = mc_surface,
  ablation_predictive_metrics = mc_ablation_results,
  ablation_surface_recovery = mc_ablation_surface,
  measurement_diagnostics = measurement_diagnostics,
  measurement_parameter_diagnostics = measurement_parameter_diagnostics,
  ablation_optimizer_attempts = ablation_optimizer_attempts,
  ablation_status = ablation_status,
  competitor_status = competitor_status,
  pattern_counts = pattern_counts,
  sampler_control_manifest = sampler_control_manifest,
  replication_metadata = replication_metadata,
  run_grid = run_grid,
  cache_tag = CACHE_TAG,
  cache_spec = CACHE_SPEC
)
saveRDS(raw_outputs, file.path(RES_DIR, paste0("study2_results_", CACHE_TAG, ".rds")))


## Persist reporting inputs before any table or figure generation.
report_names <- ls(envir = environment(), all.names = TRUE)
report_names <- setdiff(report_names, c("rep_objects", "report_names"))
report_context <- mget(report_names, envir = environment(), inherits = FALSE)
report_context <- Filter(function(x) !is.function(x) && !is.environment(x), report_context)
mixedgp_atomic_save_rds(list(schema = 1L, study = "study2", context = report_context),
        file.path(STUDY2_OUT_PREFIX, "report_inputs.rds"))
