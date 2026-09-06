############################################################
## Study II publication-design pilot
##
## Paired proxy-dimension experiment with q = 2, 4, 6 measurements of a
## shared d = 2 state, plus the q = 4 calibration curve. This is a pilot, not
## publication evidence.
############################################################

source("00_project_setup.R")
source("00_study2_functions.R")
source("03_study2_published_competitors.R")

pilot_out_dir <- file.path("..", "pilot-results", "study2")
dir.create(pilot_out_dir, recursive = TRUE, showWarnings = FALSE)

pilot_n_rep <- if (exists("PILOT_N_REP")) as.integer(PILOT_N_REP) else 2L
pilot_n <- if (exists("PILOT_STUDY2_N")) as.integer(PILOT_STUDY2_N) else 100L
pilot_n_test <- if (exists("PILOT_N_TEST")) as.integer(PILOT_N_TEST) else 250L
pilot_n_draw <- if (exists("PILOT_N_DRAW")) as.integer(PILOT_N_DRAW) else 250L
pilot_n_iter <- if (exists("PILOT_N_ITER")) as.integer(PILOT_N_ITER) else 900L
pilot_burn <- if (exists("PILOT_BURN")) as.integer(PILOT_BURN) else 350L
pilot_n_chains <- if (exists("PILOT_N_CHAINS")) {
  as.integer(PILOT_N_CHAINS)
} else {
  2L
}
pilot_m <- 4L
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

## The q = 4 rows form the calibration curve; q = 2 and 6 are evaluated at
## the common 10% calibration fraction for the proxy-dimension contrast.
pilot_design <- unique(data.frame(
  q = c(2L, 4L, 4L, 4L, 6L),
  calibration_fraction = c(0.10, 0, 0.10, 0.25, 0.10)
))
if (exists("PILOT_STUDY2_DESIGN")) {
  pilot_design <- as.data.frame(PILOT_STUDY2_DESIGN)
  if (!all(c("q", "calibration_fraction") %in% names(pilot_design))) {
    stop("PILOT_STUDY2_DESIGN needs q and calibration_fraction columns.")
  }
}
pilot_design$n_calib <- as.integer(round(pilot_n * pilot_design$calibration_fraction))

mixedgp_competitor_preflight(pilot_methods, strict = TRUE)

all_metrics <- list()
all_diagnostics <- list()
all_status <- list()
all_pattern_counts <- list()
result_id <- diagnostic_id <- status_id <- pattern_id <- 0L

for (rep_id in seq_len(pilot_n_rep)) {
  for (q in sort(unique(pilot_design$q))) {
    dat <- simulate_study2_data(
      n = pilot_n,
      n_test = pilot_n_test,
      scenario = "primary",
      q = q,
      m = pilot_m,
      seed = 1000000L + 10000L * rep_id
    )
    train <- dat$train
    test <- dat$test
    q_design <- pilot_design[pilot_design$q == q, , drop = FALSE]
    calibration_sets <- make_nested_calibration_sets_2d(
      C = train$C,
      calib_grid = q_design$n_calib,
      seed = 2000000L + rep_id,
      scheme = "random"
    )

    pattern_info <- classify_study2_pattern_frequency(train$C, test$C)
    pattern_id <- pattern_id + 1L
    all_pattern_counts[[pattern_id]] <- data.frame(
      rep = rep_id,
      q = q,
      evaluation_stratum = c("common", "rare", "unobserved"),
      n_test = as.integer(table(factor(
        pattern_info$pattern_stratum,
        levels = c("common", "rare", "unobserved")
      )))
    )

    oracle_pool <- make_oracle_pool_2d(
      dat$true_params,
      n_pool = if (q == 6L) 120000L else 80000L,
      seed = 3000000L + 10000L * rep_id + q
    )
    oracle_draws <- sample_oracle_test_y_2d(
      X_test = test$X,
      C_test = test$C,
      true_params = dat$true_params,
      sigma_eps = dat$sigma_eps,
      n_draw = pilot_n_draw,
      oracle_pool = oracle_pool,
      seed = 3100000L + 10000L * rep_id + q
    )
    result_id <- result_id + 1L
    oracle_metric <- summarize_predictive_samples_by_pattern(
      oracle_draws, test$y, pattern_info$pattern_stratum,
      "Oracle", rep_id, NA_integer_, "primary"
    )
    oracle_metric$q <- q
    all_metrics[[result_id]] <- oracle_metric

    competitor <- run_published_mixedgp_competitors(
      X_train = train$X,
      y_train = train$y,
      C_train = train$C,
      X_test = test$X,
      C_test = test$C,
      m_vec = rep(pilot_m, q),
      n_draw = pilot_n_draw,
      seed = 4000000L + 10000L * rep_id + q,
      methods = pilot_methods,
      strict = FALSE,
      controls = pilot_competitor_controls
    )
    status_id <- status_id + 1L
    competitor$status$rep <- rep_id
    competitor$status$q <- q
    competitor$status$n_calib <- NA_integer_
    all_status[[status_id]] <- competitor$status
    for (method in names(competitor$draws)) {
      result_id <- result_id + 1L
      competitor_metric <- summarize_predictive_samples_by_pattern(
        competitor$draws[[method]], test$y, pattern_info$pattern_stratum,
        method, rep_id, NA_integer_, "primary"
      )
      competitor_metric$q <- q
      all_metrics[[result_id]] <- competitor_metric
    }

    for (row_id in seq_len(nrow(q_design))) {
      n_calib <- q_design$n_calib[row_id]
      calib_idx <- calibration_sets[[as.character(n_calib)]]
      U_obs <- matrix(NA_real_, pilot_n, 2L)
      U_obs[calib_idx, ] <- train$U[calib_idx, , drop = FALSE]
      fit_start <- proc.time()[3]
      fit <- tryCatch(
        fit_eivgp_ordprobit_fb(
          X_raw = train$X,
          y_raw = train$y,
          C_ord = train$C,
          U_obs = U_obs,
          calib_idx = calib_idx,
          d = 2L,
          m_vec = rep(pilot_m, q),
          ident = "lower_triangular",
          n_iter = pilot_n_iter,
          burn = pilot_burn,
          thin = 2L,
          n_chains = pilot_n_chains,
          preset = "fast",
          sampler_strategy = "interwoven",
          store_scores = FALSE,
          seed = 5000000L + 10000L * rep_id + 100L * q + n_calib,
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
          fitted_noise_variance = NA_real_, rep = rep_id, q = q,
          n_calib = n_calib
        )
        next
      }
      all_status[[status_id]] <- data.frame(
        method = "EIV-GP", status = "success",
        optimization_status = "mcmc_completed", optimization_attempts = "",
        message = "", warnings = "",
        elapsed_seconds = elapsed,
        implementation = "exact interwoven collapsed sampler",
        fitted_noise_variance = NA_real_, rep = rep_id, q = q,
        n_calib = n_calib
      )

      diag <- as.data.frame(fit$diagnostics$summary)
      diag$rep <- rep_id
      diag$q <- q
      diag$n_calib <- n_calib
      diagnostic_id <- diagnostic_id + 1L
      all_diagnostics[[diagnostic_id]] <- diag

      n_saved <- dim(fit$mcmc$samples_U)[1L]
      draw_ids <- unique(round(seq(
        1, n_saved, length.out = min(pilot_n_draw, n_saved)
      )))
      set.seed(6000000L + 10000L * rep_id + 100L * q + n_calib)
      eiv_draws <- sample_eiv_test_y_ordprobit_fb(
        X_test_raw = test$X,
        C_test = test$C,
        fit_obj = fit,
        draw_ids = draw_ids,
        n_per_draw = 1L,
        latent_sampler = "minimax_tilting"
      )
      result_id <- result_id + 1L
      eiv_metric <- summarize_predictive_samples_by_pattern(
        eiv_draws, test$y, pattern_info$pattern_stratum,
        "EIV-GP", rep_id, n_calib, "primary"
      )
      eiv_metric$q <- q
      all_metrics[[result_id]] <- eiv_metric
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
pattern_counts <- do.call(rbind, all_pattern_counts)

overall <- metrics[metrics$evaluation_stratum == "overall", ]
fixed <- overall[overall$method %in% c("UC-GP", "LVGP", "EzGP"), ]
eiv <- overall[overall$method == "EIV-GP", ]
paired <- merge(
  eiv[, c("rep", "q", "n_calib", "CRPS", "IntervalScore95")],
  fixed[, c("rep", "q", "method", "CRPS", "IntervalScore95")],
  by = c("rep", "q"), suffixes = c("_eiv", "_competitor")
)
paired$CRPS_advantage <- paired$CRPS_competitor - paired$CRPS_eiv
paired$IntervalScore_advantage <-
  paired$IntervalScore95_competitor - paired$IntervalScore95_eiv

metric_summary <- aggregate(
  cbind(RMSE, CRPS, Coverage95, Width95, IntervalScore95) ~
    q + method + n_calib,
  data = transform(overall, n_calib = ifelse(is.na(n_calib), -1L, n_calib)),
  FUN = mean
)
advantage_summary <- aggregate(
  cbind(CRPS_advantage, IntervalScore_advantage) ~
    q + n_calib + method,
  data = paired,
  FUN = mean
)

write.csv(metrics, file.path(pilot_out_dir, "study2_pilot_metrics.csv"), row.names = FALSE)
write.csv(metric_summary, file.path(pilot_out_dir, "study2_pilot_summary.csv"), row.names = FALSE)
write.csv(paired, file.path(pilot_out_dir, "study2_pilot_paired_advantages.csv"), row.names = FALSE)
write.csv(advantage_summary, file.path(pilot_out_dir, "study2_pilot_advantage_summary.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(pilot_out_dir, "study2_pilot_diagnostics.csv"), row.names = FALSE)
write.csv(status, file.path(pilot_out_dir, "study2_pilot_status.csv"), row.names = FALSE)
write.csv(pattern_counts, file.path(pilot_out_dir, "study2_pilot_pattern_counts.csv"), row.names = FALSE)
saveRDS(
  list(
    design = list(
      n_rep = pilot_n_rep, n = pilot_n, n_test = pilot_n_test,
      grid = pilot_design, m = pilot_m, n_iter = pilot_n_iter,
      burn = pilot_burn, n_chains = pilot_n_chains, n_draw = pilot_n_draw,
      calibration_scheme = "random"
    ),
    metrics = metrics, paired = paired, diagnostics = diagnostics,
    status = status, pattern_counts = pattern_counts,
    session = sessionInfo()
  ),
  file.path(pilot_out_dir, "study2_pilot_bundle.rds")
)

print(metric_summary, row.names = FALSE)
print(advantage_summary, row.names = FALSE)
cat("Study II design pilot completed. Results:", pilot_out_dir, "\n")
