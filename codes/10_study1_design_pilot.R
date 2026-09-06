############################################################
## Study I publication-design pilot
##
## Paired experiment varying only residual within-level heterogeneity. This is
## a design and computation pilot, not publication evidence.
############################################################

source("00_project_setup.R")
source("00_study1_functions.R")
source("03_study2_published_competitors.R")

pilot_out_dir <- file.path("..", "pilot-results", "study1")
dir.create(pilot_out_dir, recursive = TRUE, showWarnings = FALSE)

pilot_n_rep <- if (exists("PILOT_N_REP")) as.integer(PILOT_N_REP) else 2L
pilot_n <- if (exists("PILOT_STUDY1_N")) as.integer(PILOT_STUDY1_N) else 100L
pilot_n_test <- if (exists("PILOT_N_TEST")) as.integer(PILOT_N_TEST) else 250L
pilot_eta <- if (exists("PILOT_STUDY1_ETA")) {
  as.numeric(PILOT_STUDY1_ETA)
} else {
  c(0, 0.5, 1)
}
pilot_calib_fraction <- c(0, 0.10, 0.25)
pilot_calib <- unique(as.integer(round(pilot_n * pilot_calib_fraction)))
pilot_m <- 6L
pilot_n_draw <- if (exists("PILOT_N_DRAW")) as.integer(PILOT_N_DRAW) else 250L
pilot_n_iter <- if (exists("PILOT_N_ITER")) as.integer(PILOT_N_ITER) else 900L
pilot_burn <- if (exists("PILOT_BURN")) as.integer(PILOT_BURN) else 350L
pilot_n_chains <- if (exists("PILOT_N_CHAINS")) {
  as.integer(PILOT_N_CHAINS)
} else {
  2L
}
pilot_methods <- if (exists("PILOT_METHODS")) {
  as.character(PILOT_METHODS)
} else {
  c("UC-GP", "LVGP", "EzGP")
}
pilot_competitor_controls <- list(
  `UC-GP` = list(n_starts = 8L),
  LVGP = list(
    n_starts = 8L, max_retries = 2L,
    max_iter_ini = 100L, max_iter_lat = 20L,
    max_elapsed_seconds = 600, parallel = FALSE
  ),
  EzGP = list(cv_folds = 2L, maxeval = 60L, cv_score = "nlpd")
)

mixedgp_competitor_preflight(pilot_methods, strict = TRUE)

all_metrics <- list()
all_diagnostics <- list()
all_status <- list()
result_id <- 0L
diagnostic_id <- 0L
status_id <- 0L

for (rep_id in seq_len(pilot_n_rep)) {
  calibration_sets <- NULL
  for (eta in pilot_eta) {
    scenario_label <- paste0("eta_", format(eta, trim = TRUE))
    dat <- simulate_1d_data(
      n = pilot_n,
      n_test = pilot_n_test,
      m = pilot_m,
      scenario = "heterogeneity_continuum",
      heterogeneity_eta = eta,
      threshold_design = "balanced",
      min_class_count = 0L,
      seed = 100000L + rep_id
    )
    train <- dat$train
    test <- dat$test
    if (is.null(calibration_sets)) {
      calibration_sets <- make_nested_calibration_sets(
        n = pilot_n,
        calib_grid = pilot_calib,
        seed = 200000L + rep_id
      )
    }

    set.seed(300000L + 1000L * rep_id + round(100 * eta))
    oracle_draws <- sample_oracle_test_y(
      x_test = test$x,
      c_test = test$c,
      tau_true = dat$tau_true,
      scenario = "heterogeneity_continuum",
      heterogeneity_eta = eta,
      sigma_eps = dat$sigma_eps,
      n_draw = pilot_n_draw
    )
    result_id <- result_id + 1L
    all_metrics[[result_id]] <- summarize_predictive_samples_1d(
      oracle_draws, test$y, "Oracle", rep_id, NA_integer_, scenario_label
    )

    competitor <- run_published_mixedgp_competitors(
      X_train = matrix(train$x, ncol = 1L),
      y_train = train$y,
      C_train = matrix(train$c, ncol = 1L),
      X_test = matrix(test$x, ncol = 1L),
      C_test = matrix(test$c, ncol = 1L),
      m_vec = pilot_m,
      n_draw = pilot_n_draw,
      seed = 400000L + 10000L * rep_id + round(100 * eta),
      methods = pilot_methods,
      strict = FALSE,
      controls = pilot_competitor_controls
    )
    status_id <- status_id + 1L
    competitor$status$rep <- rep_id
    competitor$status$eta <- eta
    competitor$status$n_calib <- NA_integer_
    all_status[[status_id]] <- competitor$status
    for (method in names(competitor$draws)) {
      result_id <- result_id + 1L
      all_metrics[[result_id]] <- summarize_predictive_samples_1d(
        competitor$draws[[method]], test$y, method, rep_id,
        NA_integer_, scenario_label
      )
    }

    for (n_calib in pilot_calib) {
      calib_idx <- calibration_sets[[as.character(n_calib)]]
      u_obs <- rep(NA_real_, pilot_n)
      u_obs[calib_idx] <- train$u[calib_idx]
      fit_start <- proc.time()[3]
      fit <- tryCatch(
        fit_eivgp_1d(
          x_raw = train$x,
          y_raw = train$y,
          c_ord = train$c,
          u_obs = u_obs,
          calib_idx = calib_idx,
          m = pilot_m,
          n_iter = pilot_n_iter,
          burn = pilot_burn,
          thin = 2L,
          n_chains = pilot_n_chains,
          preset = "fast",
          seed = 500000L + 10000L * rep_id + 100L * round(10 * eta) + n_calib,
          parallel_chains = FALSE,
          verbose = FALSE
        ),
        error = function(e) e
      )
      elapsed <- proc.time()[3] - fit_start
      status_id <- status_id + 1L
      if (inherits(fit, "error")) {
        all_status[[status_id]] <- data.frame(
          method = "EIV-GP", status = "failed",
          optimization_status = "mcmc_failed", optimization_attempts = "",
          message = conditionMessage(fit), warnings = "",
          elapsed_seconds = elapsed, implementation = "pilot",
          fitted_noise_variance = NA_real_, rep = rep_id, eta = eta,
          n_calib = n_calib
        )
        next
      }
      all_status[[status_id]] <- data.frame(
        method = "EIV-GP", status = "success",
        optimization_status = "mcmc_completed", optimization_attempts = "",
        message = "", warnings = "",
        elapsed_seconds = elapsed, implementation = "exact collapsed sampler",
        fitted_noise_variance = NA_real_, rep = rep_id, eta = eta,
        n_calib = n_calib
      )

      diag <- as.data.frame(fit$diagnostics$summary)
      diag$rep <- rep_id
      diag$eta <- eta
      diag$n_calib <- n_calib
      diagnostic_id <- diagnostic_id + 1L
      all_diagnostics[[diagnostic_id]] <- diag

      draw_ids <- unique(round(seq(
        1, nrow(fit$mcmc$samples_u),
        length.out = min(pilot_n_draw, nrow(fit$mcmc$samples_u))
      )))
      set.seed(600000L + 10000L * rep_id + 100L * round(10 * eta) + n_calib)
      eiv_draws <- sample_eiv_test_y(
        x_test_raw = test$x,
        c_test = test$c,
        fit_obj = fit,
        draw_ids = draw_ids,
        n_per_draw = 1L
      )
      result_id <- result_id + 1L
      all_metrics[[result_id]] <- summarize_predictive_samples_1d(
        eiv_draws, test$y, "EIV-GP", rep_id, n_calib, scenario_label
      )
    }
  }
}

metrics <- do.call(rbind, all_metrics)
diagnostics <- if (length(all_diagnostics)) {
  do.call(rbind, all_diagnostics)
} else {
  data.frame()
}
status <- do.call(rbind, all_status)
metrics$eta <- as.numeric(sub("eta_", "", metrics$scenario, fixed = TRUE))

fixed <- metrics[metrics$method %in% c("UC-GP", "LVGP", "EzGP"), ]
eiv <- metrics[metrics$method == "EIV-GP", ]
paired <- merge(
  eiv[, c("rep", "eta", "n_calib", "CRPS", "IntervalScore95")],
  fixed[, c("rep", "eta", "method", "CRPS", "IntervalScore95")],
  by = c("rep", "eta"), suffixes = c("_eiv", "_competitor")
)
paired$CRPS_advantage <- paired$CRPS_competitor - paired$CRPS_eiv
paired$IntervalScore_advantage <-
  paired$IntervalScore95_competitor - paired$IntervalScore95_eiv

metric_summary <- aggregate(
  cbind(RMSE, CRPS, Coverage95, Width95, IntervalScore95) ~
    eta + method + n_calib,
  data = transform(metrics, n_calib = ifelse(is.na(n_calib), -1L, n_calib)),
  FUN = mean
)
advantage_summary <- aggregate(
  cbind(CRPS_advantage, IntervalScore_advantage) ~
    eta + n_calib + method,
  data = paired,
  FUN = mean
)

write.csv(metrics, file.path(pilot_out_dir, "study1_pilot_metrics.csv"), row.names = FALSE)
write.csv(metric_summary, file.path(pilot_out_dir, "study1_pilot_summary.csv"), row.names = FALSE)
write.csv(paired, file.path(pilot_out_dir, "study1_pilot_paired_advantages.csv"), row.names = FALSE)
write.csv(advantage_summary, file.path(pilot_out_dir, "study1_pilot_advantage_summary.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(pilot_out_dir, "study1_pilot_diagnostics.csv"), row.names = FALSE)
write.csv(status, file.path(pilot_out_dir, "study1_pilot_status.csv"), row.names = FALSE)
saveRDS(
  list(
    design = list(
      n_rep = pilot_n_rep, n = pilot_n, n_test = pilot_n_test,
      eta = pilot_eta, calibration = pilot_calib, m = pilot_m,
      n_iter = pilot_n_iter, burn = pilot_burn, n_chains = pilot_n_chains,
      n_draw = pilot_n_draw
    ),
    metrics = metrics, paired = paired, diagnostics = diagnostics,
    status = status, session = sessionInfo()
  ),
  file.path(pilot_out_dir, "study1_pilot_bundle.rds")
)

print(metric_summary, row.names = FALSE)
print(advantage_summary, row.names = FALSE)
cat("Study I design pilot completed. Results:", pilot_out_dir, "\n")
