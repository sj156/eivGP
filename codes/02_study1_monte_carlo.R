## Fitting/evaluation only. Reporting lives in report_study1_results.R.
sys.source("setup_study1_experiment.R", envir = environment(), chdir = TRUE)

run_one_study1_replication <- function(rep_id, run_eiv = TRUE) {
  data_path <- file.path(
    STUDY1_DATA_DIR,
    mixedgp_dataset_filename(
      "study1", rep_id, STUDY1_SCENARIO, n_train, n_test, m
    )
  )
  if (!file.exists(data_path)) {
    stop(
      "Missing frozen Study I dataset: ", data_path,
      ". Run 12_generate_synthetic_datasets.R first."
    )
  }
  frozen <- load_mixedgp_synthetic_dataset_strict(
    data_path,
    expected = list(
      study = "study1", scenario = STUDY1_SCENARIO, rep_id = rep_id,
      n = n_train, n_test = n_test, m = m,
      calib_grid = calib_grid,
      threshold_design = STUDY1_THRESHOLD_DESIGN,
      heterogeneity_eta = STUDY1_HETEROGENEITY_ETA
    )
  )
  dat <- frozen$data
  train <- dat$train
  test <- dat$test
  calib_sets <- frozen$calibration_sets

  oracle_draws <- sample_oracle_test_y(
    x_test = test$x,
    c_test = test$c,
    tau_true = dat$tau_true,
    scenario = STUDY1_SCENARIO,
    sigma_eps = dat$sigma_eps,
    n_draw = n_pred_draw,
    heterogeneity_eta = STUDY1_HETEROGENEITY_ETA
  )
  metrics <- list(
    Oracle = summarize_predictive_samples_1d(
      oracle_draws, test$y, method = "Oracle", rep_id = rep_id,
      n_calib = NA_integer_, scenario = STUDY1_SCENARIO
    )
  )
  set.seed(240000L + rep_id)
  mean_eval_idx <- sort(sample(
    seq_len(nrow(test)), min(settings$n_m_eval, nrow(test))
  ))
  m_true <- m0_1d(
    test$x[mean_eval_idx], test$c[mean_eval_idx],
    tau = dat$tau_true, scenario = STUDY1_SCENARIO
  )
  mean_metrics <- list()
  ablation_metrics <- list()
  latent_imputation_metrics <- list()
  surface_recovery_metrics <- list()
  ablation_surface_metrics <- list()
  ablation_status <- list()
  measurement_fits <- list()
  mcmc_diagnostics <- list()
  mcmc_parameter_diagnostics <- list()
  sampler_controls <- list()
  sampler_control_manifest <- list()

  competitor_result <- mixedgp_cached_competitors(
    X_train = matrix(train$x, ncol = 1L),
    y_train = train$y,
    C_train = matrix(train$c, ncol = 1L),
    X_test = matrix(test$x, ncol = 1L),
    C_test = matrix(test$c, ncol = 1L),
    m_vec = m,
    n_draw = n_pred_draw,
    seed = 250000L + 1000L * rep_id,
    methods = STUDY1_PUBLISHED_COMPETITORS,
    strict = STUDY1_STRICT_COMPETITORS,
    controls = competitor_controls
  )
  for (method in names(competitor_result$draws)) {
    metrics[[method]] <- summarize_predictive_samples_1d(
      competitor_result$draws[[method]], test$y,
      method = method, rep_id = rep_id,
      n_calib = NA_integer_, scenario = STUDY1_SCENARIO
    )
    mean_metrics[[method]] <- summarize_mean_recovery_1d(
      matrix(
        competitor_result$latent_means[[method]][mean_eval_idx],
        nrow = 1L
      ),
      m_true = m_true,
      method = method,
      rep_id = rep_id,
      n_calib = NA_integer_,
      scenario = STUDY1_SCENARIO,
      valid_function_draws = FALSE
    )
  }

  if (isTRUE(STUDY1_RUN_ABLATIONS)) {
    if (isTRUE(STUDY1_EVALUATE_F)) {
      start_full_u <- proc.time()[3]
      full_u_fit <- tryCatch(
        fit_study1_full_u_gp(train$x, train$y, train$u),
        error = function(e) e
      )
      elapsed_full_u <- proc.time()[3] - start_full_u
      if (inherits(full_u_fit, "error")) {
        ablation_status[["Full-U GP"]] <- data.frame(
          method = "Full-U GP", n_calib = NA_integer_, status = "failed",
          message = conditionMessage(full_u_fit),
          elapsed_seconds = elapsed_full_u
        )
      } else {
        ablation_status[["Full-U GP"]] <- data.frame(
          method = "Full-U GP", n_calib = NA_integer_, status = "success",
          message = "", elapsed_seconds = elapsed_full_u
        )
        ablation_surface_metrics[["Full-U GP"]] <-
          study1_latent_gp_surface_metrics(
            gp_fit = full_u_fit,
            scenario = STUDY1_SCENARIO,
            method = "Full-U GP",
            rep_id = rep_id,
            n_calib = NA_integer_,
            tau_true = dat$tau_true,
            heterogeneity_eta = STUDY1_HETEROGENEITY_ETA,
            grid_n_x = if (STUDY1_QUICK) 15L else 31L,
            grid_n_u = if (STUDY1_QUICK) 21L else 51L
          )
      }
    }
    for (n_calib in calib_grid) {
      key <- as.character(n_calib)
      start_measurement <- proc.time()[3]
      measurement_fit <- tryCatch(
        fit_threshold_measurement_response_free(
          c_ord = train$c,
          u_obs = train$u,
          calib_idx = calib_sets[[key]],
          m = m,
          n_iter = measurement_n_iter,
          burn = measurement_burn,
          thin = 1L,
          seed = 270000L + 1000L * rep_id + n_calib
        ),
        error = function(e) e
      )
      elapsed_measurement <- proc.time()[3] - start_measurement
      if (inherits(measurement_fit, "error")) {
        msg <- conditionMessage(measurement_fit)
        ablation_status[[paste0("measurement_", key)]] <- data.frame(
          method = "Response-free threshold model", n_calib = n_calib,
          status = "failed", message = msg,
          elapsed_seconds = elapsed_measurement
        )
        ablation_status[[paste0("PI_", key)]] <- data.frame(
          method = "PI-GP", n_calib = n_calib,
          status = "failed", message = msg, elapsed_seconds = NA_real_
        )
        ablation_status[[paste0("CC_", key)]] <- data.frame(
          method = "CC-GP", n_calib = n_calib,
          status = if (n_calib == 0L) "not_applicable" else "failed",
          message = if (n_calib == 0L) "No complete cases." else msg,
          elapsed_seconds = NA_real_
        )
        next
      }
      ablation_status[[paste0("measurement_", key)]] <- data.frame(
        method = "Response-free threshold model", n_calib = n_calib,
        status = "success", message = "",
        elapsed_seconds = elapsed_measurement
      )
      measurement_fits[[key]] <- measurement_fit

      start_pi <- proc.time()[3]
      pi_result <- tryCatch({
        fit <- fit_study1_pi_gp(train$x, train$y, measurement_fit)
        draws <- sample_study1_pi_gp(
          fit, test$x, test$c, n_pred_draw,
          seed = 280000L + 1000L * rep_id + n_calib
        )
        list(fit = fit, draws = draws)
      }, error = function(e) e)
      elapsed_pi <- proc.time()[3] - start_pi
      if (inherits(pi_result, "error")) {
        ablation_status[[paste0("PI_", key)]] <- data.frame(
          method = "PI-GP", n_calib = n_calib,
          status = "failed", message = conditionMessage(pi_result),
          elapsed_seconds = elapsed_pi
        )
      } else {
        ablation_status[[paste0("PI_", key)]] <- data.frame(
          method = "PI-GP", n_calib = n_calib,
          status = "success", message = "", elapsed_seconds = elapsed_pi
        )
        ablation_metrics[[paste0("PI_", key)]] <- summarize_predictive_samples_1d(
          pi_result$draws, test$y, method = "PI-GP", rep_id = rep_id,
          n_calib = n_calib, scenario = STUDY1_SCENARIO
        )
        if (isTRUE(STUDY1_EVALUATE_F) && n_calib > 0L) {
          ablation_surface_metrics[[paste0("PI_", key)]] <-
            study1_latent_gp_surface_metrics(
              gp_fit = pi_result$fit,
              scenario = STUDY1_SCENARIO,
              method = "PI-GP",
              rep_id = rep_id,
              n_calib = n_calib,
              tau_true = dat$tau_true,
              heterogeneity_eta = STUDY1_HETEROGENEITY_ETA,
              grid_n_x = if (STUDY1_QUICK) 15L else 31L,
              grid_n_u = if (STUDY1_QUICK) 21L else 51L
            )
        }
      }

      if (n_calib == 0L) {
        ablation_status[[paste0("CC_", key)]] <- data.frame(
          method = "CC-GP", n_calib = n_calib,
          status = "not_applicable", message = "No complete cases.",
          elapsed_seconds = 0
        )
      } else {
        start_cc <- proc.time()[3]
        cc_result <- tryCatch({
          fit <- fit_study1_cc_gp(
            train$x, train$y, train$u, calib_sets[[key]], measurement_fit
          )
          draws <- sample_study1_cc_gp(
            fit, test$x, test$c, n_pred_draw,
            seed = 290000L + 1000L * rep_id + n_calib
          )
          list(fit = fit, draws = draws)
        }, error = function(e) e)
        elapsed_cc <- proc.time()[3] - start_cc
        if (inherits(cc_result, "error")) {
          ablation_status[[paste0("CC_", key)]] <- data.frame(
            method = "CC-GP", n_calib = n_calib,
            status = "failed", message = conditionMessage(cc_result),
            elapsed_seconds = elapsed_cc
          )
        } else {
          ablation_status[[paste0("CC_", key)]] <- data.frame(
            method = "CC-GP", n_calib = n_calib,
            status = "success", message = "", elapsed_seconds = elapsed_cc
          )
          ablation_metrics[[paste0("CC_", key)]] <- summarize_predictive_samples_1d(
            cc_result$draws, test$y, method = "CC-GP", rep_id = rep_id,
            n_calib = n_calib, scenario = STUDY1_SCENARIO
          )
          if (isTRUE(STUDY1_EVALUATE_F)) {
            ablation_surface_metrics[[paste0("CC_", key)]] <-
              study1_latent_gp_surface_metrics(
                gp_fit = cc_result$fit,
                scenario = STUDY1_SCENARIO,
                method = "CC-GP",
                rep_id = rep_id,
                n_calib = n_calib,
                tau_true = dat$tau_true,
                heterogeneity_eta = STUDY1_HETEROGENEITY_ETA,
                grid_n_x = if (STUDY1_QUICK) 15L else 31L,
                grid_n_u = if (STUDY1_QUICK) 21L else 51L
              )
          }
        }
      }
    }
  }

  if (isTRUE(run_eiv)) {
    for (n_calib in calib_grid) {
      message("Study I replication ", rep_id, ": EIV-GP |O|=", n_calib)
      fit_eiv <- fit_eivgp_1d(
        x_raw = train$x,
        y_raw = train$y,
        c_ord = train$c,
        u_true = train$u,
        calib_idx = calib_sets[[as.character(n_calib)]],
        m = m,
        tau_true = dat$tau_true,
        n_iter = settings$n_iter,
        burn = settings$burn,
        thin = 1L,
        n_chains = settings$n_chains,
        preset = settings$preset,
        seed = 300000L + 1000L * rep_id + n_calib,
        parallel_chains = parallel_chains,
        n_cores = chain_cores,
        verbose = FALSE
      )
      control_key <- as.character(n_calib)
      sampler_controls[[control_key]] <- list(
        control = fit_eiv$control,
        initialization_rule = paste(
          "threshold-compatible latent initialization with independent",
          "seeded perturbations; see fit_eivgp_1d()"
        ),
        chain_seeds = fit_eiv$mcmc$chain_stats$seed,
        covariance_jitter = fit_eiv$diagnostics$summary$covariance_jitter,
        forms_explicit_covariance_inverse =
          fit_eiv$diagnostics$summary$forms_explicit_covariance_inverse
      )
      sampler_control_manifest[[control_key]] <-
        study1_sampler_control_rows(fit_eiv, rep_id, n_calib)
      diag_row <- fit_eiv$diagnostics$summary
      parameter_diag <- mixedgp_study1_raw_diagnostics(
        fit_eiv, rhat_limit = STUDY1_MAX_RHAT,
        bulk_ess_limit = STUDY1_MIN_ESS, tail_ess_limit = STUDY1_MIN_ESS
      )
      raw_gate_pass <- mixedgp_diagnostic_table_pass(
        parameter_diag, rhat_limit = STUDY1_MAX_RHAT,
        bulk_ess_limit = STUDY1_MIN_ESS, tail_ess_limit = STUDY1_MIN_ESS,
        expected_parameters = names(mixedgp_study1_raw_series(fit_eiv))
      )
      panel <- mixedgp_study2_diagnostic_rows(
        matrix(test$x, ncol = 1L), matrix(test$c, ncol = 1L),
        max_points = STUDY1_DIAGNOSTIC_N_POINTS
      )
      target_series <- mixedgp_study1_target_series(
        fit_eiv, X = matrix(test$x[panel], ncol = 1L), C = test$c[panel],
        U = if (isTRUE(STUDY1_EVALUATE_F) && STUDY1_HETEROGENEITY_ETA == 1) {
          test$u[panel]
        } else NULL,
        max_draws_per_chain = STUDY1_DIAGNOSTIC_MAX_DRAWS,
        n_latent = STUDY1_DIAGNOSTIC_N_LATENT,
        seed = 380000L + 1000L * rep_id + n_calib
      )
      target_diag <- mixedgp_summarize_diagnostic_series(
        target_series, rhat_limit = STUDY1_MAX_RHAT,
        bulk_ess_limit = STUDY1_MIN_ESS, tail_ess_limit = STUDY1_MIN_ESS
      )
      target_gate_pass <- mixedgp_diagnostic_table_pass(
        target_diag, rhat_limit = STUDY1_MAX_RHAT,
        bulk_ess_limit = STUDY1_MIN_ESS, tail_ess_limit = STUDY1_MIN_ESS,
        expected_parameters = names(target_series)
      )
      parameter_diag$diagnostic_role <- "raw"
      target_diag$diagnostic_role <- "target"
      parameter_diag <- bind_rows(parameter_diag, target_diag)
      gate_pass <- raw_gate_pass && target_gate_pass
      max_rhat <- if (nrow(parameter_diag)) max(parameter_diag$rhat) else NA_real_
      min_ess <- if (nrow(parameter_diag)) {
        min(parameter_diag$ess_bulk, parameter_diag$ess_tail)
      } else NA_real_
      parameter_diag$rep <- rep_id
      parameter_diag$n_calib <- n_calib
      mcmc_parameter_diagnostics[[as.character(n_calib)]] <- parameter_diag
      diag_row$rep <- rep_id
      diag_row$n_calib <- n_calib
      diag_row$gate_max_rhat <- max_rhat
      diag_row$gate_min_ess <- min_ess
      diag_row$gate_pass <- gate_pass
      diag_row$raw_gate_pass <- raw_gate_pass
      diag_row$target_gate_pass <- target_gate_pass
      mcmc_diagnostics[[as.character(n_calib)]] <- diag_row
      diag_row$diagnostic_warning <- !isTRUE(gate_pass)
      diag_row$diagnostic_advice <- if (isTRUE(gate_pass)) "" else mixedgp_simulation_diagnostic_advice()
      mcmc_diagnostics[[as.character(n_calib)]] <- diag_row
      if (!isTRUE(gate_pass)) {
        failure_file <- tempfile(
          sprintf("mcmc_failure_rep%03d_cal%03d_", rep_id, n_calib),
          tmpdir = RES_DIR, fileext = ".rds"
        )
        saveRDS(list(fit = fit_eiv, diagnostics = parameter_diag,
                     diagnostic_panel_rows = panel, cache_spec = STUDY1_CACHE_SPEC),
                failure_file)
        warning(
          "Study I MCMC diagnostics flagged replication ", rep_id,
          ", |O|=", n_calib, ": max R-hat=", signif(max_rhat, 4),
          ", min ESS=", signif(min_ess, 5),
          ". Fit retained and diagnostics saved to ", failure_file, ". ",
          mixedgp_simulation_diagnostic_advice(), call. = FALSE
        )
      }
      draw_ids <- seq_len(nrow(fit_eiv$mcmc$samples_u))
      if (length(draw_ids) > n_pred_draw) {
        set.seed(350000L + 1000L * rep_id + n_calib)
        draw_ids <- sample(draw_ids, n_pred_draw)
      }
      eiv_draws <- sample_eiv_test_y(
        x_test_raw = test$x,
        c_test = test$c,
        fit_obj = fit_eiv,
        draw_ids = draw_ids,
        n_per_draw = 1L
      )
      metrics[[paste0("EIV_", n_calib)]] <- summarize_predictive_samples_1d(
        eiv_draws, test$y, method = "EIV-GP", rep_id = rep_id,
        n_calib = n_calib, scenario = STUDY1_SCENARIO
      )
      m_draw_ids <- draw_ids
      if (length(m_draw_ids) > settings$n_m_draw) {
        set.seed(360000L + 1000L * rep_id + n_calib)
        m_draw_ids <- sample(m_draw_ids, settings$n_m_draw)
      }
      m_draws <- sample_eiv_m_given_xc_1d(
        x_star_raw = test$x[mean_eval_idx],
        c_star = test$c[mean_eval_idx],
        fit_obj = fit_eiv,
        draw_ids = m_draw_ids,
        n_latent = settings$n_m_latent,
        include_process_uncertainty = TRUE,
        joint = FALSE,
        return_components = TRUE,
        seed = 370000L + 1000L * rep_id + n_calib
      )
      mean_metrics[[paste0("EIV_", n_calib)]] <-
        summarize_mean_recovery_1d(
          m_draws,
          m_true = m_true,
          method = "EIV-GP",
          rep_id = rep_id,
          n_calib = n_calib,
          scenario = STUDY1_SCENARIO,
          valid_function_draws = TRUE
        )
      if (isTRUE(STUDY1_EVALUATE_U) &&
          !is.null(measurement_fits[[as.character(n_calib)]])) {
        latent_imputation_metrics[[as.character(n_calib)]] <-
          study1_latent_imputation_metrics(
            eiv_fit = fit_eiv,
            measurement_fit = measurement_fits[[as.character(n_calib)]],
            U_true = train$u,
            rep_id = rep_id,
            n_calib = n_calib,
            scenario = STUDY1_SCENARIO
          )
      }
      if (isTRUE(STUDY1_EVALUATE_F) && n_calib > 0L) {
        surface_recovery_metrics[[as.character(n_calib)]] <-
          study1_eiv_surface_recovery_metrics(
            fit = fit_eiv,
            scenario = STUDY1_SCENARIO,
            rep_id = rep_id,
            n_calib = n_calib,
            tau_true = dat$tau_true,
            heterogeneity_eta = STUDY1_HETEROGENEITY_ETA,
            grid_n_x = if (STUDY1_QUICK) 15L else 31L,
            grid_n_u = if (STUDY1_QUICK) 21L else 51L,
            max_draw = if (STUDY1_QUICK) 40L else 200L,
            seed = 380000L + 1000L * rep_id + n_calib
          )
      }
    }
  }

  status <- competitor_result$status
  status$rep <- rep_id
  status$scenario <- STUDY1_SCENARIO
  ablation_status_df <- bind_rows(ablation_status)
  if (nrow(ablation_status_df) > 0L) {
    ablation_status_df$rep <- rep_id
    ablation_status_df$scenario <- STUDY1_SCENARIO
  }
  list(
    metrics = bind_rows(metrics),
    competitor_status = status,
    ablation_metrics = bind_rows(ablation_metrics),
    ablation_status = ablation_status_df,
    mcmc_diagnostics = bind_rows(mcmc_diagnostics),
    mcmc_parameter_diagnostics = bind_rows(mcmc_parameter_diagnostics),
    sampler_control_manifest = bind_rows(sampler_control_manifest),
    mean_recovery = bind_rows(mean_metrics),
    latent_imputation = bind_rows(latent_imputation_metrics),
    surface_recovery = bind_rows(surface_recovery_metrics),
    ablation_surface_recovery = bind_rows(ablation_surface_metrics),
    metadata = list(
      rep = rep_id,
      design_tag = STUDY1_DESIGN_TAG,
      data_seed = 100000L + rep_id,
      calibration_seed = 200000L + rep_id,
      frozen_data_file = basename(data_path),
      frozen_data_md5 = attr(frozen, "manifest_md5"),
      frozen_manifest = attr(frozen, "manifest_path"),
      sampler = sampler_controls
    )
  )
}

run_or_load_study1_replication <- function(rr) {
  cache_mode <- if (isTRUE(STUDY1_REUSE_LOCKED_EIV)) "archival" else "fresh"
  rep_file <- file.path(
    REP_DIR,
    sprintf("study1_rep%03d_%s_%s.rds", rr, STUDY1_DESIGN_TAG, cache_mode)
  )
  if (isTRUE(STUDY1_USE_CACHE) && file.exists(rep_file)) {
    cached <- readRDS(rep_file)
    data_path <- file.path(
      STUDY1_DATA_DIR,
      mixedgp_dataset_filename(
        "study1", rr, STUDY1_SCENARIO, n_train, n_test, m
      )
    )
    cache_valid <- is.list(cached$metadata) &&
      identical(cached$metadata$rep, as.integer(rr)) &&
      identical(cached$metadata$design_tag, STUDY1_DESIGN_TAG) &&
      identical(cached$metadata$frozen_data_file, basename(data_path)) &&
      file.exists(data_path) &&
      identical(
        tolower(as.character(cached$metadata$frozen_data_md5)),
        tolower(unname(tools::md5sum(data_path)))
      )
    if (isTRUE(cache_valid)) return(cached)
    message("Ignoring incompatible Study I cache: ", rep_file)
  }
  message("Study I replication ", rr, " of ", n_rep)
  out <- run_one_study1_replication(
    rr,
    run_eiv = !isTRUE(STUDY1_REUSE_LOCKED_EIV)
  )
  saveRDS(out, rep_file)
  out
}
rep_objects <- mixedgp_run_replications(
  as.list(seq_len(n_rep)),
  run_or_load_study1_replication,
  n_cores = replication_cores,
  seeds = 910000L + seq_len(n_rep),
  status_path = file.path(TAB_DIR, "study1_replication_status.csv"),
  study = "Study I",
  parallel_map = mixedgp_parallel_lapply,
  mc.preschedule = FALSE
)

new_results <- bind_rows(lapply(rep_objects, `[[`, "metrics"))
competitor_status <- bind_rows(
  lapply(rep_objects, `[[`, "competitor_status")
)
ablation_results <- bind_rows(lapply(rep_objects, `[[`, "ablation_metrics"))
ablation_status <- bind_rows(lapply(rep_objects, `[[`, "ablation_status"))
mcmc_diagnostics <- bind_rows(lapply(rep_objects, `[[`, "mcmc_diagnostics"))
mcmc_parameter_diagnostics <- bind_rows(
  lapply(rep_objects, `[[`, "mcmc_parameter_diagnostics")
)
sampler_control_manifest <- bind_rows(
  lapply(rep_objects, `[[`, "sampler_control_manifest")
)
mean_recovery <- bind_rows(lapply(rep_objects, `[[`, "mean_recovery"))
latent_imputation <- bind_rows(lapply(rep_objects, `[[`, "latent_imputation"))
surface_recovery <- bind_rows(lapply(rep_objects, `[[`, "surface_recovery"))
ablation_surface_recovery <- bind_rows(
  lapply(rep_objects, `[[`, "ablation_surface_recovery")
)
if (nrow(mcmc_diagnostics) > 0L) {
  write.csv(
    mcmc_diagnostics,
    file.path(TAB_DIR, "study1_mcmc_diagnostics.csv"),
    row.names = FALSE
  )
}
if (nrow(mcmc_parameter_diagnostics) > 0L) {
  write.csv(mcmc_parameter_diagnostics,
            file.path(TAB_DIR, "study1_mcmc_parameter_diagnostics.csv"),
            row.names = FALSE)
}
if (nrow(sampler_control_manifest) > 0L) {
  write.csv(
    sampler_control_manifest,
    file.path(TAB_DIR, "study1_sampler_control_manifest.csv"),
    row.names = FALSE
  )
}

if (isTRUE(STUDY1_REUSE_LOCKED_EIV)) {
  locked_file <- file.path(TAB_DIR, "study1_mc_raw_results.csv")
  if (!file.exists(locked_file)) {
    stop("Locked July 27 Study I results are missing: ", locked_file)
  }
  locked <- read.csv(locked_file, stringsAsFactors = FALSE)
  locked_eiv <- locked |>
    filter(method == "EIV-GP", rep %in% seq_len(n_rep))
  expected_rows <- n_rep * length(calib_grid)
  if (nrow(locked_eiv) != expected_rows ||
      !setequal(unique(locked_eiv$n_calib), calib_grid)) {
    stop("Locked EIV-GP rows do not match the archived 50-by-5 Study I design.")
  }
  ## Oracle values in the archived file are duplicated over calibration sizes.
  ## Keep one copy per replication for the revised table and curves.
  locked_oracle <- locked |>
    filter(method == "Oracle", rep %in% seq_len(n_rep)) |>
    group_by(rep, scenario, method) |>
    slice(1L) |>
    ungroup() |>
    mutate(n_calib = NA_integer_)
  mc_results <- bind_rows(
    locked_eiv,
    locked_oracle,
    new_results |> filter(!method %in% c("EIV-GP", "Oracle"))
  )
} else {
  mc_results <- new_results
}


## Persist reporting inputs before any table or figure generation.
report_names <- ls(envir = environment(), all.names = TRUE)
report_names <- setdiff(report_names, c("rep_objects", "report_names"))
report_context <- mget(report_names, envir = environment(), inherits = FALSE)
report_context <- Filter(function(x) !is.function(x) && !is.environment(x), report_context)
mixedgp_atomic_save_rds(list(schema = 1L, study = "study1", context = report_context),
        file.path(STUDY1_OUT_PREFIX, "report_inputs.rds"))
