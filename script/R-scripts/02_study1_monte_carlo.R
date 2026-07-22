############################################################
## 02_study1_monte_carlo.R
##
## Monte Carlo study for revised Study I.
############################################################

if (!exists("fit_eivgp_1d")) {
  source("00_study1_functions.R")
}

if (!exists("STUDY1_QUICK")) STUDY1_QUICK <- FALSE
if (!exists("STUDY1_USE_CACHE")) STUDY1_USE_CACHE <- TRUE
if (!exists("STUDY1_OUT_PREFIX")) STUDY1_OUT_PREFIX <- ".."

FIG_DIR <- file.path(STUDY1_OUT_PREFIX, "figures")
TAB_DIR <- file.path(STUDY1_OUT_PREFIX, "tables")
RES_DIR <- file.path(STUDY1_OUT_PREFIX, "results", "study1_1d")

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(20260705)

n_train <- 100L
n_test <- if (STUDY1_QUICK) 150L else 500L
m <- 6L
calib_grid <- c(0L, 5L, 10L, 20L, 50L)

n_rep <- if (STUDY1_QUICK) 2L else 50L

mc_n_iter <- if (STUDY1_QUICK) 600L else 5000L
mc_burn <- if (STUDY1_QUICK) 200L else 1000L
mc_n_chains <- if (STUDY1_QUICK) 1L else 12L
mc_preset <- if (STUDY1_QUICK) "fast" else "balanced"

n_pred_draw <- if (STUDY1_QUICK) 150L else 600L

parallel_chains <- (
  .Platform$OS.type != "windows" &&
    parallel::detectCores(logical = TRUE) > 2L
)

method_cols <- c(
  "Oracle" = "black",
  "EIV-GP" = "firebrick",
  "GP-LearnedEmb" = "purple4",
  "GP-CondMean" = "steelblue",
  "GP-Gaussian" = "darkgreen"
)

############################################################
## One replication
############################################################

run_one_1d_replication <- function(rep_id,
                                   scenario = "active",
                                   calib_grid = c(0, 5, 10, 20, 50),
                                   n = 100,
                                   n_test = 500,
                                   m = 6,
                                   eiv_n_iter = 6000,
                                   eiv_burn = 1000,
                                   eiv_n_chains = 12L,
                                   eiv_preset = "balanced",
                                   n_pred_draw = 600,
                                   parallel_chains = FALSE) {
  dat <- simulate_1d_data(
    n = n,
    n_test = n_test,
    m = m,
    scenario = scenario,
    seed = 100000L + rep_id
  )
  
  train <- dat$train
  test <- dat$test
  
  calib_sets <- make_nested_calibration_sets(
    n = n,
    calib_grid = calib_grid,
    seed = 200000L + rep_id
  )
  
  baselines <- fit_embedding_baselines(
    x_raw = train$x,
    y_raw = train$y,
    c_ord = train$c,
    m = m,
    n_starts_learned = if (STUDY1_QUICK) 2L else 6L
  )
  
  baseline_draws <- predict_embedding_baseline_samples(
    baselines = baselines,
    x_star_raw = test$x,
    c_star = test$c,
    m = m,
    n_draw = n_pred_draw
  )
  
  oracle_draws <- sample_oracle_test_y(
    x_test = test$x,
    c_test = test$c,
    tau_true = dat$tau_true,
    scenario = scenario,
    sigma_eps = dat$sigma_eps,
    n_draw = n_pred_draw
  )
  
  out_metrics <- list()
  
  for (n_calib in calib_grid) {
    out_metrics[[paste0("oracle_", n_calib)]] <-
      summarize_predictive_samples(
        oracle_draws,
        test$y,
        method = "Oracle",
        rep_id = rep_id,
        n_calib = n_calib,
        scenario = scenario
      )
    
    for (nm in names(baseline_draws)) {
      out_metrics[[paste0(nm, "_", n_calib)]] <-
        summarize_predictive_samples(
          baseline_draws[[nm]],
          test$y,
          method = nm,
          rep_id = rep_id,
          n_calib = n_calib,
          scenario = scenario
        )
    }
  }
  
  for (n_calib in calib_grid) {
    cat("Replication", rep_id, ": fitting EIV-GP with |O| =", n_calib, "\n")
    
    fit_eiv <- fit_eivgp_1d(
      x_raw = train$x,
      y_raw = train$y,
      c_ord = train$c,
      u_true = train$u,
      calib_idx = calib_sets[[as.character(n_calib)]],
      m = m,
      tau_true = dat$tau_true,
      n_iter = eiv_n_iter,
      burn = eiv_burn,
      thin = 1L,
      n_chains = eiv_n_chains,
      preset = eiv_preset,
      seed = 300000L + 1000L * rep_id + n_calib,
      parallel_chains = parallel_chains,
      verbose = FALSE
    )
    
    draw_ids <- seq_len(nrow(fit_eiv$mcmc$samples_u))
    
    if (length(draw_ids) > n_pred_draw) {
      draw_ids <- sample(draw_ids, n_pred_draw)
    }
    
    eiv_draws <- sample_eiv_test_y(
      x_test_raw = test$x,
      c_test = test$c,
      fit_obj = fit_eiv,
      draw_ids = draw_ids,
      n_per_draw = 1L
    )
    
    out_metrics[[paste0("EIV_", n_calib)]] <-
      summarize_predictive_samples(
        eiv_draws,
        test$y,
        method = "EIV-GP",
        rep_id = rep_id,
        n_calib = n_calib,
        scenario = scenario
      )
  }
  
  do.call(rbind, out_metrics)
}

############################################################
## Run Monte Carlo
############################################################

mc_file <- file.path(
  RES_DIR,
  paste0("study1_mc_results_", ifelse(STUDY1_QUICK, "quick", "paper"), ".rds")
)

if (STUDY1_USE_CACHE && file.exists(mc_file)) {
  mc_results <- readRDS(mc_file)
} else {
  mc_results_list <- vector("list", n_rep)
  
  for (rr in seq_len(n_rep)) {
    cat("\n========== Monte Carlo replication", rr, "of", n_rep, "==========\n")
    
    mc_results_list[[rr]] <- run_one_1d_replication(
      rep_id = rr,
      scenario = "active",
      calib_grid = calib_grid,
      n = n_train,
      n_test = n_test,
      m = m,
      eiv_n_iter = mc_n_iter,
      eiv_burn = mc_burn,
      eiv_n_chains = mc_n_chains,
      eiv_preset = mc_preset,
      n_pred_draw = n_pred_draw,
      parallel_chains = parallel_chains
    )
  }
  
  mc_results <- dplyr::bind_rows(mc_results_list)
  
  saveRDS(mc_results, mc_file)
}

write.csv(
  mc_results,
  file.path(TAB_DIR, "study1_mc_raw_results.csv"),
  row.names = FALSE
)

############################################################
## Summaries, plot, LaTeX table
############################################################

method_levels <- c("Oracle", "EIV-GP", "GP-LearnedEmb", "GP-CondMean", "GP-Gaussian")

mc_results$method <- factor(mc_results$method, levels = method_levels)

mc_long <- mc_results |>
  tidyr::pivot_longer(
    cols = c(RMSE, MAE, Coverage95, Width95, CRPS, IntervalScore95),
    names_to = "metric",
    values_to = "value"
  )

metric_levels <- c("RMSE", "MAE", "CRPS", "Coverage95", "Width95", "IntervalScore95")
mc_long$metric <- factor(mc_long$metric, levels = metric_levels)

mc_summary <- mc_long |>
  dplyr::group_by(scenario, n_calib, method, metric) |>
  dplyr::summarise(
    mean = mean(value, na.rm = TRUE),
    se = sd(value, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

p_mc <- ggplot(
  mc_summary,
  aes(x = n_calib, y = mean, color = method, group = method)
) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
    width = 1.2,
    alpha = 0.7
  ) +
  geom_hline(
    data = data.frame(metric = factor("Coverage95", levels = metric_levels), yint = 0.95),
    aes(yintercept = yint),
    inherit.aes = FALSE,
    linetype = "dashed",
    color = "gray35"
  ) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  scale_color_manual(values = method_cols, name = NULL) +
  labs(
    x = "Number of calibrated latent observations",
    y = "Monte Carlo mean",
    title = "Study I: test-set performance across calibration sizes"
  ) +
  theme(legend.position = "bottom")

ggsave(
  file.path(FIG_DIR, "fig5_study1_mc_metrics.pdf"),
  p_mc,
  width = 12,
  height = 7.5
)

mc_summary_table <- mc_results |>
  dplyr::group_by(n_calib, method) |>
  dplyr::summarise(
    RMSE_mean = mean(RMSE, na.rm = TRUE),
    RMSE_se = sd(RMSE, na.rm = TRUE) / sqrt(dplyr::n()),
    CRPS_mean = mean(CRPS, na.rm = TRUE),
    CRPS_se = sd(CRPS, na.rm = TRUE) / sqrt(dplyr::n()),
    Coverage_mean = mean(Coverage95, na.rm = TRUE),
    Coverage_se = sd(Coverage95, na.rm = TRUE) / sqrt(dplyr::n()),
    Width_mean = mean(Width95, na.rm = TRUE),
    Width_se = sd(Width95, na.rm = TRUE) / sqrt(dplyr::n()),
    IntervalScore_mean = mean(IntervalScore95, na.rm = TRUE),
    IntervalScore_se = sd(IntervalScore95, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    RMSE = sprintf("%.3f (%.3f)", RMSE_mean, RMSE_se),
    CRPS = sprintf("%.3f (%.3f)", CRPS_mean, CRPS_se),
    Coverage95 = sprintf("%.3f (%.3f)", Coverage_mean, Coverage_se),
    Width95 = sprintf("%.3f (%.3f)", Width_mean, Width_se),
    IntervalScore95 = sprintf("%.3f (%.3f)", IntervalScore_mean, IntervalScore_se)
  ) |>
  dplyr::select(
    n_calib,
    Method = method,
    RMSE,
    CRPS,
    Coverage95,
    Width95,
    IntervalScore95
  )

latex_mc_table <- knitr::kable(
  mc_summary_table,
  format = "latex",
  booktabs = TRUE,
  align = "llccccc",
  escape = FALSE,
  col.names = c(
    "$|\\mathcal O|$",
    "Method",
    "RMSE",
    "CRPS",
    "Coverage95",
    "Width95",
    "IntervalScore95"
  )
)

writeLines(
  latex_mc_table,
  con = file.path(TAB_DIR, "study1_mc_summary.tex")
)

cat("\nMonte Carlo Study I figure written to:\n")
cat(normalizePath(file.path(FIG_DIR, "fig5_study1_mc_metrics.pdf")), "\n")
cat("\nMonte Carlo Study I table written to:\n")
cat(normalizePath(file.path(TAB_DIR, "study1_mc_summary.tex")), "\n")