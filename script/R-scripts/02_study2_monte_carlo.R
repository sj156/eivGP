############################################################
## 02_study2_monte_carlo.R
##
## Monte Carlo study for revised Study II.
##
## Uses exact fully Bayesian ordinal-probit EIV-GP:
##
##   fit_eivgp_ordprobit_fb()
##
## This script is resumable: each replication is cached
## separately under results/study2_ordprobit_exact/mc_replications/.
############################################################

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

if (!exists("fit_eivgp_ordprobit_fb")) {
  source("00_study2_functions.R")
}

needed_pkgs <- c("ggplot2", "patchwork", "dplyr", "tidyr", "knitr")
missing_pkgs <- needed_pkgs[
  !vapply(needed_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0L) {
  stop("Please install required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(knitr)
})

if (!exists("STUDY2_CONFIG")) STUDY2_CONFIG <- "quick"
if (!exists("STUDY2_USE_CACHE")) STUDY2_USE_CACHE <- TRUE
if (!exists("STUDY2_OUT_PREFIX")) STUDY2_OUT_PREFIX <- ".."
if (!exists("STUDY2_MC_RESUME")) STUDY2_MC_RESUME <- TRUE
if (!exists("STUDY2_SAVE_REP_FITS")) STUDY2_SAVE_REP_FITS <- FALSE
if (!exists("STUDY2_SAVE_PDF")) STUDY2_SAVE_PDF <- TRUE
if (!exists("STUDY2_SAVE_PNG")) STUDY2_SAVE_PNG <- TRUE
if (!exists("STUDY2_MC_N_REP")) STUDY2_MC_N_REP <- NULL
if (!exists("STUDY2_MC_CALIB_GRID")) {
  STUDY2_MC_CALIB_GRID <- c(0L, 10L, 25L, 50L, 80L)
}

settings <- study2_config_settings(STUDY2_CONFIG)

############################################################
## Directories
############################################################

FIG_DIR <- file.path(STUDY2_OUT_PREFIX, "figures")
TAB_DIR <- file.path(STUDY2_OUT_PREFIX, "tables")
RES_DIR <- file.path(STUDY2_OUT_PREFIX, "results", "study2_ordprobit_exact")
MC_FIG_DIR <- file.path(FIG_DIR, "study2_mc_extra")
MC_REP_DIR <- file.path(RES_DIR, "mc_replications")
MC_FIT_DIR <- file.path(RES_DIR, "mc_fits")

for (dd in c(FIG_DIR, TAB_DIR, RES_DIR, MC_FIG_DIR, MC_REP_DIR, MC_FIT_DIR)) {
  dir.create(dd, showWarnings = FALSE, recursive = TRUE)
}

############################################################
## Study settings
############################################################

set.seed(20260710)

n_train <- 120L
n_test <- settings$n_test
m_vec <- rep(4L, 4L)
d_latent <- 2L
ident_method <- "lower_triangular"

calib_grid <- as.integer(STUDY2_MC_CALIB_GRID)

n_rep <- if (is.null(STUDY2_MC_N_REP)) {
  settings$n_rep
} else {
  as.integer(STUDY2_MC_N_REP)
}

mc_n_iter <- settings$mc_n_iter
mc_burn <- settings$mc_burn
mc_thin <- settings$mc_thin
mc_n_chains <- settings$mc_n_chains
mc_preset <- settings$preset

n_pred_draw <- settings$n_pred_draw
n_new_latent_gibbs <- settings$n_new_latent_gibbs
n_oracle_pool <- settings$n_oracle_pool
n_starts_learned <- settings$n_starts_learned

parallel_chains <- (
  .Platform$OS.type != "windows" &&
    parallel::detectCores(logical = TRUE) > 2L
)

CACHE_TAG <- paste0(
  STUDY2_CONFIG,
  "_", STUDY2_DESIGN_TAG,
  "_exact_fb",
  "_nrep", n_rep,
  "_ntrain", n_train,
  "_ntest", n_test,
  "_iter", mc_n_iter,
  "_burn", mc_burn,
  "_chains", mc_n_chains
)

method_levels <- c(
  "Oracle",
  "EIV-GP",
  "GP-LearnedEmb",
  "GP-CondMean",
  "GP-Gaussian"
)

embedding_methods <- c("GP-LearnedEmb", "GP-CondMean", "GP-Gaussian")

method_cols <- c(
  "Oracle" = "black",
  "EIV-GP" = "firebrick",
  "GP-LearnedEmb" = "purple4",
  "GP-CondMean" = "steelblue",
  "GP-Gaussian" = "darkgreen"
)

metric_levels <- c(
  "RMSE",
  "MAE",
  "CRPS",
  "Coverage95",
  "Width95",
  "IntervalScore95"
)

############################################################
## Helpers
############################################################

save_plot <- function(path_no_ext,
                      plot,
                      width,
                      height,
                      dpi = 320) {
  dir.create(dirname(path_no_ext), showWarnings = FALSE, recursive = TRUE)
  
  if (isTRUE(STUDY2_SAVE_PDF)) {
    ggsave(
      filename = paste0(path_no_ext, ".pdf"),
      plot = plot,
      width = width,
      height = height,
      units = "in",
      limitsize = FALSE
    )
  }
  
  if (isTRUE(STUDY2_SAVE_PNG)) {
    ggsave(
      filename = paste0(path_no_ext, ".png"),
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      bg = "white",
      limitsize = FALSE
    )
  }
  
  invisible(plot)
}

extract_fit_diagnostics <- function(fit, rep_id, n_calib) {
  out <- fit$diagnostics$summary
  out$rep <- rep_id
  out$n_calib <- n_calib
  out |>
    dplyr::select(rep, n_calib, dplyr::everything())
}

extract_chain_stats <- function(fit, rep_id, n_calib) {
  cs <- fit$mcmc$chain_stats
  
  if (is.null(cs) || nrow(cs) == 0L) {
    return(data.frame())
  }
  
  cs |>
    dplyr::mutate(
      rep = rep_id,
      n_calib = n_calib,
      u_ess_accept_rate = u_ess_accept_total / pmax(u_ess_total, 1),
      global_u_accept_rate = global_u_accept_total / pmax(global_u_total, 1),
      u_eval_per_update = u_ess_eval_total / pmax(u_ess_total, 1),
      global_u_eval_per_update = global_u_eval_total / pmax(global_u_total, 1),
      theta_eval_per_update = theta_eval_total / pmax(theta_update_total, 1)
    ) |>
    dplyr::select(rep, n_calib, dplyr::everything())
}

extract_imputation_metrics_mc <- function(fit, rep_id, n_calib) {
  if (is.null(fit$data$U_true_eval)) return(data.frame())
  
  U_sum <- posterior_U_summary_eivgp(
    fit_obj = fit,
    original_scale = FALSE,
    true_U = fit$data$U_true_eval
  )
  
  out_coord <- summarize_U_imputation_metrics(U_sum)
  if (nrow(out_coord) == 0L) return(data.frame())
  
  out_overall <- U_sum |>
    dplyr::filter(!calibrated) |>
    dplyr::summarise(
      coord = "overall",
      n = dplyr::n(),
      bias = mean(error, na.rm = TRUE),
      rmse = sqrt(mean(error^2, na.rm = TRUE)),
      mae = mean(abs_error, na.rm = TRUE),
      coverage95 = mean(covered95, na.rm = TRUE),
      mean_width95 = mean(q975 - q025, na.rm = TRUE),
      .groups = "drop"
    )
  
  dplyr::bind_rows(out_coord, out_overall) |>
    dplyr::mutate(
      rep = rep_id,
      n_calib = n_calib
    ) |>
    dplyr::select(rep, n_calib, dplyr::everything())
}

############################################################
## One replication
############################################################

run_one_study2_replication <- function(rep_id,
                                       calib_grid,
                                       n = 120,
                                       n_test = 400,
                                       m_vec = rep(4L, 4L),
                                       d_latent = 2L,
                                       ident = "lower_triangular",
                                       eiv_n_iter = 3000L,
                                       eiv_burn = 1000L,
                                       eiv_thin = 2L,
                                       eiv_n_chains = 4L,
                                       eiv_preset = "balanced",
                                       n_pred_draw = 500L,
                                       n_new_latent_gibbs = 20L,
                                       n_oracle_pool = 150000L,
                                       n_starts_learned = 4L,
                                       parallel_chains = FALSE,
                                       save_rep_fits = FALSE,
                                       fit_dir = NULL) {
  dat <- simulate_study2_data(
    n = n,
    n_test = n_test,
    seed = 100000L + rep_id
  )
  
  train <- dat$train
  test <- dat$test
  
  calib_sets <- make_cell_stratified_calibration_sets_2d(
    C = train$C,
    calib_grid = calib_grid,
    anchor_cols = c(1L, 2L),
    seed = 200000L + rep_id
  )
  
  oracle_pool <- make_oracle_pool_2d(
    true_params = dat$true_params,
    n_pool = n_oracle_pool,
    seed = 250000L + rep_id
  )
  
  ##########################################################
  ## Baselines and oracle
  ##########################################################
  
  baselines <- fit_embedding_baselines_ord(
    X_raw = train$X,
    y_raw = train$y,
    C_ord = train$C,
    m_vec = m_vec,
    n_starts_learned = n_starts_learned
  )
  
  baseline_draws <- predict_embedding_baseline_samples_ord(
    baselines = baselines,
    X_star_raw = test$X,
    C_star = test$C,
    n_draw = n_pred_draw
  )
  
  oracle_draws <- sample_oracle_test_y_2d(
    X_test = test$X,
    C_test = test$C,
    true_params = dat$true_params,
    sigma_eps = dat$sigma_eps,
    n_draw = n_pred_draw,
    oracle_pool = oracle_pool
  )
  
  metrics_list <- list()
  diagnostics_list <- list()
  chain_stats_list <- list()
  imputation_metrics_list <- list()
  
  for (n_calib in calib_grid) {
    metrics_list[[paste0("oracle_", n_calib)]] <-
      summarize_predictive_samples(
        oracle_draws,
        test$y,
        method = "Oracle",
        rep_id = rep_id,
        n_calib = n_calib,
        scenario = "study2"
      )
    
    for (nm in names(baseline_draws)) {
      metrics_list[[paste0(nm, "_", n_calib)]] <-
        summarize_predictive_samples(
          baseline_draws[[nm]],
          test$y,
          method = nm,
          rep_id = rep_id,
          n_calib = n_calib,
          scenario = "study2"
        )
    }
  }
  
  ##########################################################
  ## EIV-GP fits across calibration sizes
  ##########################################################
  
  for (n_calib in calib_grid) {
    cat("Replication", rep_id, ": fitting exact EIV-GP with |O| =", n_calib, "\n")
    
    fit_eiv <- fit_eivgp_ordprobit_fb(
      X_raw = train$X,
      y_raw = train$y,
      C_ord = train$C,
      U_obs = train$U,
      calib_idx = calib_sets[[as.character(n_calib)]],
      U_true_eval = train$U,
      d = d_latent,
      m_vec = m_vec,
      ident = ident,
      n_iter = eiv_n_iter,
      burn = eiv_burn,
      thin = eiv_thin,
      n_chains = eiv_n_chains,
      preset = eiv_preset,
      seed = 300000L + 1000L * rep_id + n_calib,
      parallel_chains = parallel_chains,
      verbose = FALSE
    )
    
    diagnostics_list[[as.character(n_calib)]] <-
      extract_fit_diagnostics(fit_eiv, rep_id = rep_id, n_calib = n_calib)
    
    chain_stats_list[[as.character(n_calib)]] <-
      extract_chain_stats(fit_eiv, rep_id = rep_id, n_calib = n_calib)
    
    imputation_metrics_list[[as.character(n_calib)]] <-
      extract_imputation_metrics_mc(
        fit = fit_eiv,
        rep_id = rep_id,
        n_calib = n_calib
      )
    
    if (isTRUE(save_rep_fits)) {
      if (is.null(fit_dir)) fit_dir <- "."
      dir.create(fit_dir, showWarnings = FALSE, recursive = TRUE)
      
      saveRDS(
        fit_eiv,
        file.path(
          fit_dir,
          sprintf("fit_eiv_rep%03d_calib%03d.rds", rep_id, n_calib)
        )
      )
    }
    
    draw_ids <- seq_len(dim(fit_eiv$mcmc$samples_U)[1])
    if (length(draw_ids) > n_pred_draw) {
      draw_ids <- sample(draw_ids, n_pred_draw)
    }
    
    eiv_draws <- sample_eiv_test_y_ordprobit_fb(
      X_test_raw = test$X,
      C_test = test$C,
      fit_obj = fit_eiv,
      draw_ids = draw_ids,
      n_per_draw = 1L,
      n_new_latent_gibbs = n_new_latent_gibbs
    )
    
    metrics_list[[paste0("EIV_", n_calib)]] <-
      summarize_predictive_samples(
        eiv_draws,
        test$y,
        method = "EIV-GP",
        rep_id = rep_id,
        n_calib = n_calib,
        scenario = "study2"
      )
    
    rm(fit_eiv, eiv_draws)
    invisible(gc())
  }
  
  list(
    metrics = dplyr::bind_rows(metrics_list),
    diagnostics = dplyr::bind_rows(diagnostics_list),
    chain_stats = dplyr::bind_rows(chain_stats_list),
    imputation_metrics = dplyr::bind_rows(imputation_metrics_list),
    metadata = list(
      rep_id = rep_id,
      n = n,
      n_test = n_test,
      calib_grid = calib_grid,
      m_vec = m_vec,
      d_latent = d_latent,
      ident = ident,
      eiv_n_iter = eiv_n_iter,
      eiv_burn = eiv_burn,
      eiv_thin = eiv_thin,
      eiv_n_chains = eiv_n_chains,
      eiv_preset = eiv_preset,
      n_pred_draw = n_pred_draw
    )
  )
}

############################################################
## Run Monte Carlo with per-replication cache
############################################################

mc_file <- file.path(
  RES_DIR,
  paste0("study2_mc_results_", CACHE_TAG, ".rds")
)

mc_diag_file <- file.path(
  RES_DIR,
  paste0("study2_mc_diagnostics_", CACHE_TAG, ".rds")
)

rep_files <- file.path(
  MC_REP_DIR,
  sprintf("study2_mc_rep_%s_rep%03d.rds", CACHE_TAG, seq_len(n_rep))
)

if (STUDY2_USE_CACHE && file.exists(mc_file) && file.exists(mc_diag_file)) {
  mc_results <- readRDS(mc_file)
  mc_diag_obj <- readRDS(mc_diag_file)
  
  mc_diagnostics <- mc_diag_obj$diagnostics
  mc_chain_stats <- mc_diag_obj$chain_stats
  mc_imputation_metrics <- mc_diag_obj$imputation_metrics
} else {
  rep_objs <- vector("list", n_rep)
  
  for (rr in seq_len(n_rep)) {
    cat("\n========== Study II Monte Carlo replication", rr, "of", n_rep, "==========\n")
    
    if (
      isTRUE(STUDY2_MC_RESUME) &&
      STUDY2_USE_CACHE &&
      file.exists(rep_files[rr])
    ) {
      cat("Using cached replication file:\n", rep_files[rr], "\n")
      rep_objs[[rr]] <- readRDS(rep_files[rr])
    } else {
      rep_objs[[rr]] <- run_one_study2_replication(
        rep_id = rr,
        calib_grid = calib_grid,
        n = n_train,
        n_test = n_test,
        m_vec = m_vec,
        d_latent = d_latent,
        ident = ident_method,
        eiv_n_iter = mc_n_iter,
        eiv_burn = mc_burn,
        eiv_thin = mc_thin,
        eiv_n_chains = mc_n_chains,
        eiv_preset = mc_preset,
        n_pred_draw = n_pred_draw,
        n_new_latent_gibbs = n_new_latent_gibbs,
        n_oracle_pool = n_oracle_pool,
        n_starts_learned = n_starts_learned,
        parallel_chains = parallel_chains,
        save_rep_fits = STUDY2_SAVE_REP_FITS,
        fit_dir = MC_FIT_DIR
      )
      
      saveRDS(rep_objs[[rr]], rep_files[rr])
    }
  }
  
  mc_results <- dplyr::bind_rows(lapply(rep_objs, function(z) z$metrics))
  mc_diagnostics <- dplyr::bind_rows(lapply(rep_objs, function(z) z$diagnostics %||% data.frame()))
  mc_chain_stats <- dplyr::bind_rows(lapply(rep_objs, function(z) z$chain_stats %||% data.frame()))
  mc_imputation_metrics <- dplyr::bind_rows(lapply(rep_objs, function(z) z$imputation_metrics %||% data.frame()))
  
  saveRDS(mc_results, mc_file)
  
  saveRDS(
    list(
      diagnostics = mc_diagnostics,
      chain_stats = mc_chain_stats,
      imputation_metrics = mc_imputation_metrics
    ),
    mc_diag_file
  )
}

############################################################
## Save raw outputs
############################################################

write.csv(
  mc_results,
  file.path(TAB_DIR, paste0("study2_mc_raw_results_", CACHE_TAG, ".csv")),
  row.names = FALSE
)

if (nrow(mc_diagnostics) > 0L) {
  write.csv(
    mc_diagnostics,
    file.path(TAB_DIR, paste0("study2_mc_fit_diagnostics_", CACHE_TAG, ".csv")),
    row.names = FALSE
  )
}

if (nrow(mc_chain_stats) > 0L) {
  write.csv(
    mc_chain_stats,
    file.path(TAB_DIR, paste0("study2_mc_chain_update_stats_", CACHE_TAG, ".csv")),
    row.names = FALSE
  )
}

if (nrow(mc_imputation_metrics) > 0L) {
  write.csv(
    mc_imputation_metrics,
    file.path(TAB_DIR, paste0("study2_mc_imputation_metrics_", CACHE_TAG, ".csv")),
    row.names = FALSE
  )
}

############################################################
## Predictive-performance summaries
############################################################

mc_results$method <- factor(mc_results$method, levels = method_levels)

mc_long <- mc_results |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(metric_levels),
    names_to = "metric",
    values_to = "value"
  )

mc_long$metric <- factor(mc_long$metric, levels = metric_levels)

mc_summary <- mc_long |>
  dplyr::group_by(scenario, n_calib, method, metric) |>
  dplyr::summarise(
    mean = mean(value, na.rm = TRUE),
    se = safe_se(value),
    n_rep_eff = sum(is.finite(value)),
    .groups = "drop"
  )

write.csv(
  mc_summary,
  file.path(TAB_DIR, paste0("study2_mc_metric_summary_long_", CACHE_TAG, ".csv")),
  row.names = FALSE
)

p_mc <- ggplot(
  mc_summary,
  aes(x = n_calib, y = mean, color = method, group = method)
) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
    width = 1.5,
    alpha = 0.7
  ) +
  geom_hline(
    data = data.frame(
      metric = factor("Coverage95", levels = metric_levels),
      yint = 0.95
    ),
    aes(yintercept = yint),
    inherit.aes = FALSE,
    linetype = "dashed",
    color = "gray35"
  ) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  scale_color_manual(values = method_cols, name = NULL, drop = FALSE) +
  labs(
    x = "Number of calibrated latent observations",
    y = "Monte Carlo mean",
    title = paste0(
      "Study II: test-set performance across calibration sizes (",
      STUDY2_CONFIG,
      ")"
    )
  ) +
  theme(legend.position = "bottom")

save_plot(
  file.path(FIG_DIR, paste0("fig5_study2_mc_metrics_", CACHE_TAG)),
  p_mc,
  width = 12,
  height = 7.5
)

############################################################
## Full LaTeX summary table
############################################################

mc_summary_table <- mc_results |>
  dplyr::group_by(n_calib, method) |>
  dplyr::summarise(
    RMSE_mean = mean(RMSE, na.rm = TRUE),
    RMSE_se = safe_se(RMSE),
    CRPS_mean = mean(CRPS, na.rm = TRUE),
    CRPS_se = safe_se(CRPS),
    Coverage_mean = mean(Coverage95, na.rm = TRUE),
    Coverage_se = safe_se(Coverage95),
    Width_mean = mean(Width95, na.rm = TRUE),
    Width_se = safe_se(Width95),
    IntervalScore_mean = mean(IntervalScore95, na.rm = TRUE),
    IntervalScore_se = safe_se(IntervalScore95),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    RMSE = format_mean_se(RMSE_mean, RMSE_se),
    CRPS = format_mean_se(CRPS_mean, CRPS_se),
    Coverage95 = format_mean_se(Coverage_mean, Coverage_se),
    Width95 = format_mean_se(Width_mean, Width_se),
    IntervalScore95 = format_mean_se(IntervalScore_mean, IntervalScore_se)
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
  con = file.path(TAB_DIR, paste0("study2_mc_summary_", CACHE_TAG, ".tex"))
)

write.csv(
  mc_summary_table,
  file.path(TAB_DIR, paste0("study2_mc_summary_", CACHE_TAG, ".csv")),
  row.names = FALSE
)

############################################################
## Compact paper-style table
############################################################

summary_wide_numeric <- mc_results |>
  dplyr::group_by(n_calib, method) |>
  dplyr::summarise(
    RMSE_mean = mean(RMSE, na.rm = TRUE),
    RMSE_se = safe_se(RMSE),
    CRPS_mean = mean(CRPS, na.rm = TRUE),
    CRPS_se = safe_se(CRPS),
    Coverage_mean = mean(Coverage95, na.rm = TRUE),
    Coverage_se = safe_se(Coverage95),
    Score_mean = mean(IntervalScore95, na.rm = TRUE),
    Score_se = safe_se(IntervalScore95),
    .groups = "drop"
  )

oracle_row <- summary_wide_numeric |>
  dplyr::filter(method == "Oracle", n_calib == min(calib_grid)) |>
  dplyr::mutate(n_calib_print = "--", Method = "Oracle")

embed_ref <- summary_wide_numeric |>
  dplyr::filter(method %in% embedding_methods, n_calib == min(calib_grid))

best_embedding_row <- data.frame(
  n_calib = min(calib_grid),
  method = "Best embedding",
  RMSE_mean = min(embed_ref$RMSE_mean, na.rm = TRUE),
  RMSE_se = embed_ref$RMSE_se[which.min(embed_ref$RMSE_mean)],
  CRPS_mean = min(embed_ref$CRPS_mean, na.rm = TRUE),
  CRPS_se = embed_ref$CRPS_se[which.min(embed_ref$CRPS_mean)],
  Coverage_mean = embed_ref$Coverage_mean[
    which.min(abs(embed_ref$Coverage_mean - 0.95))
  ],
  Coverage_se = embed_ref$Coverage_se[
    which.min(abs(embed_ref$Coverage_mean - 0.95))
  ],
  Score_mean = min(embed_ref$Score_mean, na.rm = TRUE),
  Score_se = embed_ref$Score_se[which.min(embed_ref$Score_mean)],
  n_calib_print = "--",
  Method = "Best embedding"
)

eiv_rows <- summary_wide_numeric |>
  dplyr::filter(method == "EIV-GP") |>
  dplyr::mutate(
    n_calib_print = as.character(n_calib),
    Method = "EIV-GP"
  )

compact_table <- dplyr::bind_rows(
  oracle_row,
  best_embedding_row,
  eiv_rows
) |>
  dplyr::mutate(
    RMSE = format_mean_se(RMSE_mean, RMSE_se),
    CRPS = format_mean_se(CRPS_mean, CRPS_se),
    Coverage95 = format_mean_se(Coverage_mean, Coverage_se),
    Score95 = format_mean_se(Score_mean, Score_se)
  ) |>
  dplyr::select(
    `$|\\mathcal O|$` = n_calib_print,
    Method,
    RMSE,
    CRPS,
    Coverage95,
    Score95
  )

latex_compact_table <- knitr::kable(
  compact_table,
  format = "latex",
  booktabs = TRUE,
  align = "llcccc",
  escape = FALSE
)

writeLines(
  latex_compact_table,
  con = file.path(TAB_DIR, paste0("study2_mc_compact_summary_", CACHE_TAG, ".tex"))
)

writeLines(
  latex_compact_table,
  con = file.path(TAB_DIR, "study2_mc_compact_summary_exact.tex")
)

write.csv(
  compact_table,
  file.path(TAB_DIR, paste0("study2_mc_compact_summary_", CACHE_TAG, ".csv")),
  row.names = FALSE
)

############################################################
## Latent-imputation Monte Carlo summaries
############################################################

if (nrow(mc_imputation_metrics) > 0L) {
  mc_imputation_summary <- mc_imputation_metrics |>
    dplyr::group_by(n_calib, coord) |>
    dplyr::summarise(
      rmse_mean = mean(rmse, na.rm = TRUE),
      rmse_se = safe_se(rmse),
      mae_mean = mean(mae, na.rm = TRUE),
      mae_se = safe_se(mae),
      coverage95_mean = mean(coverage95, na.rm = TRUE),
      coverage95_se = safe_se(coverage95),
      width95_mean = mean(mean_width95, na.rm = TRUE),
      width95_se = safe_se(mean_width95),
      .groups = "drop"
    )
  
  write.csv(
    mc_imputation_summary,
    file.path(TAB_DIR, paste0("study2_mc_imputation_summary_", CACHE_TAG, ".csv")),
    row.names = FALSE
  )
  
  writeLines(
    knitr::kable(
      mc_imputation_summary,
      format = "latex",
      booktabs = TRUE,
      digits = 3,
      escape = TRUE
    ),
    con = file.path(TAB_DIR, "study2_mc_imputation_summary_exact.tex")
  )
}

############################################################
## Diagnostic summaries
############################################################

if (nrow(mc_diagnostics) > 0L) {
  diag_cols <- c(
    "max_rhat_hyper",
    "max_rhat_A",
    "max_rhat_tau",
    "median_rhat_missing_U",
    "max_rhat_missing_U",
    "min_ess_key",
    "mean_u_ess_accept",
    "mean_global_u_accept",
    "time_seconds"
  )
  
  diag_cols <- intersect(diag_cols, names(mc_diagnostics))
  
  mc_diag_long <- mc_diagnostics |>
    dplyr::select(dplyr::all_of(c("rep", "n_calib", diag_cols))) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(diag_cols),
      names_to = "diagnostic",
      values_to = "value"
    )
  
  mc_diag_summary <- mc_diag_long |>
    dplyr::group_by(n_calib, diagnostic) |>
    dplyr::summarise(
      mean = mean(value, na.rm = TRUE),
      se = safe_se(value),
      median = median(value, na.rm = TRUE),
      q90 = quantile(value, 0.90, na.rm = TRUE),
      .groups = "drop"
    )
  
  write.csv(
    mc_diag_summary,
    file.path(TAB_DIR, paste0("study2_mc_diagnostic_summary_", CACHE_TAG, ".csv")),
    row.names = FALSE
  )
}

############################################################
## Console summary
############################################################

cat("\nMonte Carlo Study II metric figure written to:\n")
cat(normalizePath(file.path(FIG_DIR, paste0("fig5_study2_mc_metrics_", CACHE_TAG, ".pdf"))), "\n")

cat("\nCompact Monte Carlo Study II table written to:\n")
cat(normalizePath(file.path(TAB_DIR, "study2_mc_compact_summary_exact.tex")), "\n")

cat("\nMonte Carlo Study II raw and diagnostic CSVs written to:\n")
cat(normalizePath(TAB_DIR), "\n")