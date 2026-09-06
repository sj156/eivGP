############################################################
## 02_ocean_replicates.R
##
## Same role as 02_study1_monte_carlo.R, but the train/test split
## is frozen. Each "replication" redraws a nested class-proportional
## calibration set on the same 100/300 ocean split.
############################################################

if (!exists("OCEAN_REALDATA_DIR")) {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script_file <- if (length(file_arg) > 0L) {
    normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = TRUE)
  } else {
    normalizePath(sys.frame(1)$ofile, mustWork = TRUE)
  }
  OCEAN_REALDATA_DIR <- dirname(script_file)
}
if (!exists("fit_eivgp_1d")) {
  source(file.path(OCEAN_REALDATA_DIR, "..", "00_study1_functions.R"))
}
if (!exists("load_ocean_prepared")) {
  source(file.path(OCEAN_REALDATA_DIR, "ocean_data_helpers.R"))
}

if (!exists("OCEAN_QUICK")) OCEAN_QUICK <- FALSE
if (!exists("OCEAN_USE_CACHE")) OCEAN_USE_CACHE <- TRUE
if (!exists("OCEAN_OUT_PREFIX")) {
  OCEAN_OUT_PREFIX <- file.path(OCEAN_REALDATA_DIR, "..", "..")
}
if (!exists("OCEAN_KERNEL")) OCEAN_KERNEL <- "se"
if (!exists("OCEAN_MATERN_NU")) OCEAN_MATERN_NU <- 2.5
if (!exists("OCEAN_CACHE_VERSION")) OCEAN_CACHE_VERSION <- "v2_same-target"
require_study1_reporting_packages(
  c("ggplot2", "dplyr", "tidyr", "knitr"),
  "ocean calibration-redraw analysis"
)

FIG_DIR <- file.path(OCEAN_OUT_PREFIX, "figures", "ocean")
TAB_DIR <- file.path(OCEAN_OUT_PREFIX, "tables", "ocean")
RES_DIR <- file.path(OCEAN_OUT_PREFIX, "results", "ocean")

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)

ocean <- load_ocean_prepared()
alloc <- read_ocean_class_allocation()
m <- ocean$m
calib_grid <- c(10L, 25L, 50L)

n_rep <- if (OCEAN_QUICK) 2L else 10L
mc_n_iter <- if (OCEAN_QUICK) 300L else 20000L
mc_burn <- if (OCEAN_QUICK) 80L else 5000L
mc_n_chains <- if (OCEAN_QUICK) 1L else 4L
mc_preset <- if (OCEAN_QUICK) "fast" else "thorough"
n_pred_draw <- if (OCEAN_QUICK) 80L else 600L
base_seed <- 20260813L

parallel_chains <- (
  .Platform$OS.type != "windows" &&
    parallel::detectCores(logical = TRUE) > 2L
)

method_cols <- c(
  "EIV-GP" = "firebrick",
  "GP-LearnedEmb" = "purple4",
  "GP-CondMean" = "steelblue",
  "GP-Gaussian" = "darkgreen",
  "GP-XOnly" = "gray40",
  "Complete-case GP" = "dodgerblue4",
  "Full-U GP oracle" = "black"
)

############################################################
## One redraw of nested calibration on the frozen split
############################################################

run_one_ocean_replication <- function(rep_id,
                                      ocean,
                                      allocation,
                                      calib_grid = c(10L, 25L, 50L),
                                      eiv_n_iter = 20000L,
                                      eiv_burn = 5000L,
                                      eiv_n_chains = 4L,
                                      eiv_preset = "thorough",
                                      n_pred_draw = 600L,
                                      parallel_chains = FALSE,
                                      base_seed = 20260813L,
                                      kernel = "se",
                                      matern_nu = 2.5) {
  calib_seed <- as.integer(base_seed * 100L + rep_id)
  mcmc_seed <- as.integer((base_seed + 1L) * 100L + rep_id)

  calib_sets <- make_ocean_nested_calibration(
    c_train = ocean$c_train,
    allocation = allocation,
    calib_grid = calib_grid,
    seed = calib_seed
  )

  ## These methods do not depend on the calibration redraw. Fix their random
  ## seeds across replications so the reported redraw variation is not Monte
  ## Carlo or optimizer noise.
  set.seed(base_seed + 10L)
  baselines <- fit_embedding_baselines(
    x_raw = ocean$X_train,
    y_raw = ocean$y_train,
    c_ord = ocean$c_train,
    m = ocean$m,
    n_starts_learned = if (exists("OCEAN_QUICK") && OCEAN_QUICK) 2L else 8L,
    kernel = kernel,
    matern_nu = matern_nu
  )
  set.seed(base_seed + 11L)
  baseline_draws <- predict_embedding_baseline_samples(
    baselines = baselines,
    x_star_raw = ocean$X_test,
    c_star = ocean$c_test,
    m = ocean$m,
    n_draw = n_pred_draw
  )

  X_train_std <- sweep(
    sweep(ocean$X_train, 2L, baselines$x_center, "-"),
    2L, baselines$x_scale, "/"
  )
  X_test_std <- sweep(
    sweep(ocean$X_test, 2L, baselines$x_center, "-"),
    2L, baselines$x_scale, "/"
  )
  y_train_std <- as.numeric(
    (ocean$y_train - baselines$y_center) / baselines$y_scale
  )

  set.seed(base_seed + 20L)
  x_only <- gp_mle_fit_1d(
    X_train_std, y_train_std, kernel = kernel, matern_nu = matern_nu
  )
  x_only_draws <- baselines$y_center + baselines$y_scale *
    sample_gp_mle_predictive(x_only, Xstar = X_test_std, n_draw = n_pred_draw)

  full_u <- fit_ocean_latent_gp(
    ocean$X_train,
    ocean$y_train,
    ocean$u_train,
    kernel = kernel,
    matern_nu = matern_nu
  )
  full_u_draws <- sample_ocean_latent_gp_y_given_xc(
    full_u,
    ocean$X_test,
    ocean$c_test,
    tau = ocean$tau_reference,
    n_draw = n_pred_draw,
    seed = base_seed + 30L
  )

  out_metrics <- list()
  for (n_calib in calib_grid) {
    for (nm in names(baseline_draws)) {
      out_metrics[[paste0(nm, "_", n_calib)]] <-
        summarize_predictive_samples_1d(
          baseline_draws[[nm]],
          ocean$y_test,
          method = nm,
          rep_id = rep_id,
          n_calib = n_calib,
          scenario = "ocean"
        )
    }
    out_metrics[[paste0("xonly_", n_calib)]] <-
      summarize_predictive_samples_1d(
        x_only_draws, ocean$y_test, method = "GP-XOnly",
        rep_id = rep_id, n_calib = n_calib, scenario = "ocean"
      )
    out_metrics[[paste0("fullu_", n_calib)]] <-
      summarize_predictive_samples_1d(
        full_u_draws, ocean$y_test, method = "Full-U GP oracle",
        rep_id = rep_id, n_calib = n_calib, scenario = "ocean"
      )
  }

  u_rows <- list()
  for (n_calib in calib_grid) {
    cat("Ocean redraw", rep_id, ": fitting EIV-GP with |O| =", n_calib, "\n")
    calib_idx <- calib_sets[[as.character(n_calib)]]
    u_obs <- rep(NA_real_, length(ocean$u_train))
    u_obs[calib_idx] <- ocean$u_train[calib_idx]
    fit_eiv <- fit_eivgp_1d(
      x_raw = ocean$X_train,
      y_raw = ocean$y_train,
      c_ord = ocean$c_train,
      calib_idx = calib_idx,
      m = ocean$m,
      tau_true = ocean$tau_reference,
      n_iter = eiv_n_iter,
      burn = eiv_burn,
      thin = 1L,
      n_chains = eiv_n_chains,
      preset = eiv_preset,
      seed = mcmc_seed + n_calib,
      parallel_chains = parallel_chains,
      verbose = FALSE,
      kernel = kernel,
      matern_nu = matern_nu,
      u_obs = u_obs
    )

    draw_ids <- seq_len(nrow(fit_eiv$mcmc$samples_u))
    if (length(draw_ids) > n_pred_draw) {
      set.seed(mcmc_seed + 2000L + n_calib)
      draw_ids <- sort(sample(draw_ids, n_pred_draw))
    }
    eiv_draws <- sample_eiv_test_y(
      x_test_raw = ocean$X_test,
      c_test = ocean$c_test,
      fit_obj = fit_eiv,
      draw_ids = draw_ids,
      n_per_draw = 1L
    )
    out_metrics[[paste0("EIV_", n_calib)]] <-
      summarize_predictive_samples_1d(
        eiv_draws, ocean$y_test, method = "EIV-GP",
        rep_id = rep_id, n_calib = n_calib, scenario = "ocean"
      )

    cc_fit <- fit_ocean_latent_gp(
      ocean$X_train,
      ocean$y_train,
      ocean$u_train,
      train_idx = calib_idx,
      kernel = kernel,
      matern_nu = matern_nu
    )
    cc_draws <- sample_ocean_latent_gp_y_given_xc(
      cc_fit,
      ocean$X_test,
      ocean$c_test,
      tau = estimate_tau_from_ordinal_codes(ocean$c_train, ocean$m),
      n_draw = n_pred_draw,
      seed = mcmc_seed + 3000L + n_calib
    )
    out_metrics[[paste0("CC_", n_calib)]] <-
      summarize_predictive_samples_1d(
        cc_draws, ocean$y_test, method = "Complete-case GP",
        rep_id = rep_id, n_calib = n_calib, scenario = "ocean"
      )

    u_sum <- summarize_ocean_u_imputation(
      fit_eiv, ocean$u_train, ocean$c_train, ocean$inverse_u
    )
    u_rows[[as.character(n_calib)]] <- data.frame(
      rep_id = rep_id,
      n_calib = n_calib,
      RMSE_u = sqrt(mean((u_sum$post_mean - u_sum$true_u)^2)),
      Coverage95_u = mean(u_sum$covered95)
    )
  }

  pred <- do.call(rbind, out_metrics)
  u_tab <- do.call(rbind, u_rows)
  list(prediction = pred, imputation = u_tab)
}

############################################################
## Run redraws
############################################################

mc_file <- file.path(
  RES_DIR,
  paste0(
    "ocean_rep_results_",
    OCEAN_CACHE_VERSION, "_",
    ifelse(OCEAN_QUICK, "quick", "paper"), "_",
    normalize_gp_kernel_1d(OCEAN_KERNEL, OCEAN_MATERN_NU)$name,
    ".rds"
  )
)

if (OCEAN_USE_CACHE && file.exists(mc_file)) {
  mc_pack <- readRDS(mc_file)
} else {
  pred_list <- vector("list", n_rep)
  u_list <- vector("list", n_rep)
  for (rr in seq_len(n_rep)) {
    cat("\n========== Ocean redraw", rr, "of", n_rep, "==========\n")
    one <- run_one_ocean_replication(
      rep_id = rr,
      ocean = ocean,
      allocation = alloc,
      calib_grid = calib_grid,
      eiv_n_iter = mc_n_iter,
      eiv_burn = mc_burn,
      eiv_n_chains = mc_n_chains,
      eiv_preset = mc_preset,
      n_pred_draw = n_pred_draw,
      parallel_chains = parallel_chains,
      base_seed = base_seed,
      kernel = OCEAN_KERNEL,
      matern_nu = OCEAN_MATERN_NU
    )
    pred_list[[rr]] <- one$prediction
    u_list[[rr]] <- one$imputation
  }
  mc_pack <- list(
    prediction = dplyr::bind_rows(pred_list),
    imputation = dplyr::bind_rows(u_list)
  )
  saveRDS(mc_pack, mc_file)
}

mc_results <- mc_pack$prediction
write.csv(
  mc_results,
  file.path(TAB_DIR, "ocean_rep_raw_prediction.csv"),
  row.names = FALSE
)
write.csv(
  mc_pack$imputation,
  file.path(TAB_DIR, "ocean_rep_raw_imputation.csv"),
  row.names = FALSE
)

############################################################
## Summaries, plot, LaTeX table
############################################################

method_levels <- c(
  "Full-U GP oracle", "EIV-GP", "Complete-case GP",
  "GP-LearnedEmb", "GP-CondMean", "GP-Gaussian", "GP-XOnly"
)
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
  dplyr::group_by(n_calib, method, metric) |>
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
    y = "Mean over calibration redraws",
    title = "Ocean: test-set performance across nested |O|"
  ) +
  theme(legend.position = "bottom")

ggsave(
  file.path(FIG_DIR, "fig5_ocean_rep_metrics.pdf"),
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

writeLines(
  knitr::kable(
    mc_summary_table,
    format = "latex",
    booktabs = TRUE,
    align = "llccccc",
    escape = FALSE,
    col.names = c(
      "$|\\mathcal O|$", "Method", "RMSE", "CRPS",
      "Coverage95", "Width95", "IntervalScore95"
    )
  ),
  con = file.path(TAB_DIR, "ocean_rep_summary.tex")
)

write.csv(mc_summary_table, file.path(TAB_DIR, "ocean_rep_summary.csv"), row.names = FALSE)

cat("\nOcean redraw figure written to:\n")
cat(normalizePath(file.path(FIG_DIR, "fig5_ocean_rep_metrics.pdf")), "\n")
cat("\nOcean redraw table written to:\n")
cat(normalizePath(file.path(TAB_DIR, "ocean_rep_summary.tex")), "\n")
