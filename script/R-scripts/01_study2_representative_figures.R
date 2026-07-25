############################################################
## 01_study2_representative_figures.R
##
## Representative-data figures and diagnostics for revised
## Study II.
##
## Uses exact fully Bayesian ordinal-probit EIV-GP.
############################################################

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
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(knitr)
})

if (!exists("STUDY2_CONFIG")) STUDY2_CONFIG <- "quick"
if (!exists("STUDY2_USE_CACHE")) STUDY2_USE_CACHE <- TRUE
if (!exists("STUDY2_OUT_PREFIX")) STUDY2_OUT_PREFIX <- ".."
if (!exists("STUDY2_SAVE_PDF")) STUDY2_SAVE_PDF <- TRUE
if (!exists("STUDY2_SAVE_PNG")) STUDY2_SAVE_PNG <- TRUE
if (!exists("STUDY2_RESEARCH_PLOTS")) STUDY2_RESEARCH_PLOTS <- TRUE

settings <- study2_config_settings(STUDY2_CONFIG)

############################################################
## Directories
############################################################

FIG_DIR <- file.path(STUDY2_OUT_PREFIX, "figures")
TAB_DIR <- file.path(STUDY2_OUT_PREFIX, "tables")
RES_DIR <- file.path(STUDY2_OUT_PREFIX, "results", "study2_ordprobit_exact")
EXP_FIG_DIR <- file.path(FIG_DIR, "study2_research_extra")
MCMC_FIG_DIR <- file.path(EXP_FIG_DIR, "mcmc")

for (dd in c(FIG_DIR, TAB_DIR, RES_DIR, EXP_FIG_DIR, MCMC_FIG_DIR)) {
  dir.create(dd, showWarnings = FALSE, recursive = TRUE)
}

############################################################
## Settings
############################################################

set.seed(20260710)

n_train <- 120L
n_test <- settings$n_test
m_vec <- rep(4L, 4L)
d_latent <- 2L
calib_grid <- c(0L, 10L, 25L, 50L, 80L)

rep_n_iter <- settings$rep_n_iter
rep_burn <- settings$rep_burn
rep_thin <- settings$rep_thin
rep_n_chains <- settings$rep_n_chains
rep_preset <- settings$preset

n_pred_draw <- settings$n_pred_draw
n_density_draw <- settings$n_density_draw
n_new_latent_gibbs <- settings$n_new_latent_gibbs
n_oracle_pool <- settings$n_oracle_pool
n_starts_learned <- settings$n_starts_learned

n_surface_draw <- if (STUDY2_CONFIG == "quick") 40L else 200L

parallel_chains <- (
  .Platform$OS.type != "windows" &&
    parallel::detectCores(logical = TRUE) > 2L
)

method_levels <- c(
  "Oracle",
  "EIV-GP",
  "GP-LearnedEmb",
  "GP-CondMean",
  "GP-Gaussian"
)

method_cols <- c(
  "Oracle" = "black",
  "EIV-GP" = "firebrick",
  "GP-LearnedEmb" = "purple4",
  "GP-CondMean" = "steelblue",
  "GP-Gaussian" = "darkgreen"
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

save_csv <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(x, path, row.names = FALSE)
  invisible(x)
}

calib_label <- function(x, levels = rep_fit_calibs) {
  factor(paste0("|O| = ", x), levels = paste0("|O| = ", levels))
}

make_long_draws <- function(draw_mat, method, grid) {
  draw_mat <- as.matrix(draw_mat)
  stopifnot(ncol(draw_mat) == nrow(grid))
  
  cbind(
    grid[rep(seq_len(nrow(grid)), each = nrow(draw_mat)), , drop = FALSE],
    data.frame(
      y = as.vector(draw_mat),
      method = method
    )
  )
}

make_trace_df <- function(fit, n_calib) {
  samples_by_chain <- fit$mcmc$samples_by_chain
  p <- fit$data$p
  d <- fit$data$d
  
  dplyr::bind_rows(
    lapply(seq_along(samples_by_chain$logtheta), function(cc) {
      lt <- as.matrix(samples_by_chain$logtheta[[cc]])
      
      out <- data.frame(
        chain = factor(cc),
        draw = seq_len(nrow(lt)),
        sigma_epsilon = sqrt(samples_by_chain$sigma2[[cc]]),
        rho = exp(lt[, 1]),
        n_calib = n_calib,
        n_calib_label = calib_label(n_calib)
      )
      
      for (j in seq_len(p)) {
        out[[paste0("theta_x", j)]] <- exp(lt[, 1 + j])
      }
      
      for (k in seq_len(d)) {
        out[[paste0("theta_u", k)]] <- exp(lt[, 1 + p + k])
      }
      
      out
    })
  )
}

############################################################
## Calibration sizes to fit
############################################################

if (exists("STUDY2_REP_FIT_CALIBS")) {
  rep_fit_calibs <- sort(unique(as.integer(STUDY2_REP_FIT_CALIBS)))
} else {
  rep_fit_calibs <- if (isTRUE(STUDY2_RESEARCH_PLOTS)) {
    calib_grid
  } else {
    c(0L, 25L, 80L)
  }
}

if (exists("STUDY2_MAIN_CALIB")) {
  main_calib <- as.integer(STUDY2_MAIN_CALIB)
} else {
  main_calib <- if (25L %in% rep_fit_calibs) 25L else rep_fit_calibs[1]
}

rep_fit_calibs <- sort(unique(c(rep_fit_calibs, main_calib)))
calib_grid_for_sets <- sort(unique(c(calib_grid, rep_fit_calibs)))

############################################################
## Cache files
############################################################

rep_data_file <- file.path(
  RES_DIR,
  paste0("representative_data_", STUDY2_DESIGN_TAG, "_", STUDY2_CONFIG, ".rds")
)

rep_fits_file <- file.path(
  RES_DIR,
  paste0("representative_eiv_fits_", STUDY2_DESIGN_TAG, "_", STUDY2_CONFIG, ".rds")
)

rep_base_file <- file.path(
  RES_DIR,
  paste0("representative_baselines_", STUDY2_DESIGN_TAG, "_", STUDY2_CONFIG, ".rds")
)

############################################################
## Simulate representative data
############################################################

if (STUDY2_USE_CACHE && file.exists(rep_data_file)) {
  rep_dat <- readRDS(rep_data_file)
} else {
  rep_dat <- simulate_study2_data(
    n = n_train,
    n_test = n_test,
    seed = 20260710
  )
  saveRDS(rep_dat, rep_data_file)
}

train <- rep_dat$train
test <- rep_dat$test

calib_sets <- make_cell_stratified_calibration_sets_2d(
  C = train$C,
  calib_grid = calib_grid_for_sets,
  anchor_cols = c(1L, 2L),
  seed = 20260711
)

############################################################
## Fit EIV-GP for requested calibration sizes
############################################################

if (STUDY2_USE_CACHE && file.exists(rep_fits_file)) {
  rep_fits <- readRDS(rep_fits_file)
  if (!is.list(rep_fits)) rep_fits <- list()
} else {
  rep_fits <- list()
}

missing_fit_calibs <- setdiff(as.character(rep_fit_calibs), names(rep_fits))

if (length(missing_fit_calibs) > 0L) {
  for (kk_chr in missing_fit_calibs) {
    kk <- as.integer(kk_chr)
    
    cat("Fitting representative exact EIV-GP with |O| =", kk, "\n")
    
    rep_fits[[kk_chr]] <- fit_eivgp_ordprobit_fb(
      X_raw = train$X,
      y_raw = train$y,
      C_ord = train$C,
      U_obs = train$U,
      calib_idx = calib_sets[[as.character(kk)]],
      U_true_eval = train$U,
      d = d_latent,
      m_vec = m_vec,
      ident = "lower_triangular",
      n_iter = rep_n_iter,
      burn = rep_burn,
      thin = rep_thin,
      n_chains = rep_n_chains,
      preset = rep_preset,
      seed = 500000L + kk,
      parallel_chains = parallel_chains,
      verbose = TRUE
    )
  }
  
  saveRDS(rep_fits, rep_fits_file)
}

rep_fits <- rep_fits[as.character(rep_fit_calibs)]
main_fit <- rep_fits[[as.character(main_calib)]]

############################################################
## Fit deterministic embedding baselines
############################################################

if (STUDY2_USE_CACHE && file.exists(rep_base_file)) {
  rep_baselines <- readRDS(rep_base_file)
} else {
  rep_baselines <- fit_embedding_baselines_ord(
    X_raw = train$X,
    y_raw = train$y,
    C_ord = train$C,
    m_vec = m_vec,
    n_starts_learned = n_starts_learned
  )
  saveRDS(rep_baselines, rep_base_file)
}

############################################################
## Oracle pool
############################################################

oracle_pool <- make_oracle_pool_2d(
  true_params = rep_dat$true_params,
  n_pool = n_oracle_pool,
  seed = 20260712
)

############################################################
## Representative data figure
############################################################

df_train <- data.frame(
  id = seq_len(n_train),
  y = train$y,
  severity = rowSums(train$C),
  pattern = make_pattern_label(train$C),
  x1 = train$X[, 1],
  x2 = train$X[, 2],
  u1 = train$U[, 1],
  u2 = train$U[, 2],
  f_true = train$f
)

p_latent <- ggplot(df_train, aes(x = u1, y = u2, color = y)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_viridis_c(name = "y") +
  labs(
    x = expression(u[1]),
    y = expression(u[2]),
    title = "Training data in true latent coordinates"
  )

p_x <- ggplot(df_train, aes(x = x1, y = x2, color = severity)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_viridis_c(name = "Ordinal severity") +
  labs(
    x = expression(x[1]),
    y = expression(x[2]),
    title = "Observed quantitative design"
  )

p_y_sev <- ggplot(df_train, aes(x = severity, y = y)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.75) +
  labs(
    x = "Observed ordinal severity score",
    y = "y",
    title = "Response versus observed severity"
  )

p_data <- (p_latent + p_x) / p_y_sev +
  patchwork::plot_layout(heights = c(1, 0.8))

save_plot(
  file.path(FIG_DIR, "fig1_study2_data_representative"),
  p_data,
  width = 10.5,
  height = 8
)

############################################################
## Latent imputation figure
############################################################

extract_imputation_df <- function(fit, n_calib) {
  samples_U <- fit$mcmc$samples_U
  miss_idx <- fit$data$miss_idx
  U_true_eval <- fit$data$U_true_eval
  
  if (length(miss_idx) == 0L || is.null(U_true_eval)) return(data.frame())
  
  d <- fit$data$d
  out_list <- vector("list", d)
  
  for (coord in seq_len(d)) {
    smat <- samples_U[, miss_idx, coord, drop = FALSE]
    smat <- matrix(smat, nrow = dim(samples_U)[1])
    
    post_mean <- colMeans(smat, na.rm = TRUE)
    post_lo <- apply(smat, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
    post_hi <- apply(smat, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
    true_u <- U_true_eval[miss_idx, coord]
    
    out_list[[coord]] <- data.frame(
      n_calib = n_calib,
      n_calib_label = calib_label(n_calib),
      id = miss_idx,
      coord = paste0("u", coord),
      true_u = true_u,
      post_mean = post_mean,
      post_lo = post_lo,
      post_hi = post_hi,
      width = post_hi - post_lo,
      error = post_mean - true_u,
      abs_error = abs(post_mean - true_u),
      covered = true_u >= post_lo & true_u <= post_hi,
      severity = rowSums(fit$data$C_ord[miss_idx, , drop = FALSE])
    )
  }
  
  dplyr::bind_rows(out_list)
}

df_imp_all <- dplyr::bind_rows(
  lapply(names(rep_fits), function(nm) {
    extract_imputation_df(rep_fits[[nm]], as.integer(nm))
  })
)

if (nrow(df_imp_all) > 0L) {
  df_imp_all$n_calib_label <- calib_label(df_imp_all$n_calib)
  
  p_imp <- ggplot(df_imp_all, aes(x = true_u, y = post_mean, color = severity)) +
    geom_errorbar(aes(ymin = post_lo, ymax = post_hi), width = 0, alpha = 0.16) +
    geom_point(size = 1.6, alpha = 0.85) +
    geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 0.6) +
    facet_grid(coord ~ n_calib_label) +
    scale_color_viridis_c(name = "Severity") +
    labs(
      x = "True latent coordinate",
      y = "Posterior mean and 95% interval",
      title = "Study II: latent imputation across calibration sizes"
    ) +
    theme(legend.position = "bottom")
  
  save_plot(
    file.path(FIG_DIR, "fig2_study2_latent_imputation_by_calibration"),
    p_imp,
    width = max(11, 2.8 * length(rep_fit_calibs)),
    height = 6.5
  )
  
  imp_metrics <- df_imp_all |>
    dplyr::group_by(n_calib, n_calib_label, coord) |>
    dplyr::summarise(
      n = dplyr::n(),
      bias = mean(error),
      rmse = sqrt(mean(error^2)),
      mae = mean(abs_error),
      coverage95 = mean(covered),
      mean_width95 = mean(width),
      .groups = "drop"
    )
  
  save_csv(
    imp_metrics,
    file.path(RES_DIR, "representative_imputation_metrics_by_coordinate.csv")
  )
  
  writeLines(
    knitr::kable(imp_metrics, format = "latex", booktabs = TRUE, digits = 3),
    con = file.path(TAB_DIR, "study2_imputation_metrics_by_coordinate.tex")
  )
}

############################################################
## Latent response-surface recovery figure
############################################################

surface_calibs_main <- intersect(c(0L, 25L, 80L), as.integer(names(rep_fits)))

df_surface_all <- dplyr::bind_rows(
  lapply(as.character(surface_calibs_main), function(nm) {
    fit <- rep_fits[[nm]]
    
    out <- make_latent_surface_recovery_2d(
      fit = fit,
      truth_fun = f0_2d,
      x_ref_raw = c(0, 0),
      u_lim = c(-2.2, 2.2),
      grid_size = if (STUDY2_CONFIG == "quick") 35L else 55L,
      max_draw = n_surface_draw
    )
    
    out$n_calib <- as.integer(nm)
    out$n_calib_label <- calib_label(as.integer(nm))
    out
  })
)

df_truth_surface <- df_surface_all |>
  dplyr::filter(n_calib == surface_calibs_main[1]) |>
  dplyr::select(u1, u2, truth) |>
  dplyr::distinct() |>
  dplyr::mutate(
    panel = "Truth",
    value = truth
  )

df_eiv_surface <- df_surface_all |>
  dplyr::mutate(
    panel = paste0("EIV-GP, |O| = ", n_calib),
    value = mean
  ) |>
  dplyr::select(u1, u2, panel, value)

df_surface_plot_main <- dplyr::bind_rows(
  df_truth_surface[, c("u1", "u2", "panel", "value")],
  df_eiv_surface
)

panel_levels <- c(
  "Truth",
  paste0("EIV-GP, |O| = ", surface_calibs_main)
)

df_surface_plot_main$panel <- factor(df_surface_plot_main$panel, levels = panel_levels)

p_surface_main <- ggplot(
  df_surface_plot_main,
  aes(x = u1, y = u2)  # ← Remove fill from global aes()
) +
  geom_raster(aes(fill = value), interpolate = TRUE) +  # ← Add fill here only
  geom_contour(
    aes(z = value),
    color = "white",
    linewidth = 0.20,
    alpha = 0.45
  ) +
  facet_wrap(~panel, nrow = 1) +
  scale_fill_gradient2(
    low = "navy",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    name = expression(f(x,u))
  ) +
  labs(
    x = expression(u[1]),
    y = expression(u[2]),
    title = "Latent response-surface recovery at x = (0,0)"
  ) +
  theme(legend.position = "bottom")

save_plot(
  file.path(FIG_DIR, "fig4_study2_surface_recovery_selected"),
  p_surface_main,
  width = max(10, 2.8 * length(panel_levels)),
  height = 4.8
)

surface_recovery_summary <- df_surface_all |>
  dplyr::group_by(n_calib) |>
  dplyr::summarise(
    surface_rmse = sqrt(mean(error^2, na.rm = TRUE)),
    surface_mae = mean(abs(error), na.rm = TRUE),
    surface_coverage95 = mean(covered95, na.rm = TRUE),
    surface_mean_width95 = mean(width95, na.rm = TRUE),
    .groups = "drop"
  )

save_csv(
  surface_recovery_summary,
  file.path(RES_DIR, "representative_surface_recovery_summary.csv")
)

writeLines(
  knitr::kable(surface_recovery_summary, format = "latex", booktabs = TRUE, digits = 3),
  con = file.path(TAB_DIR, "study2_surface_recovery_summary.tex")
)

############################################################
## Predictive-density figure
############################################################

X_selected <- rbind(
  c(-0.5, 0),
  c(0, 0),
  c(0.5, 0)
)

C_selected <- rbind(
  rep(1L, length(m_vec)),
  pmax(1L, ceiling(m_vec / 2)),
  m_vec
)

selected_grid <- expand.grid(
  x_id = seq_len(nrow(X_selected)),
  c_id = seq_len(nrow(C_selected))
)

X_star_selected <- X_selected[selected_grid$x_id, , drop = FALSE]
C_star_selected <- C_selected[selected_grid$c_id, , drop = FALSE]

selected_grid$x_label <- paste0(
  "x*=(",
  X_selected[selected_grid$x_id, 1],
  ",",
  X_selected[selected_grid$x_id, 2],
  ")"
)

selected_grid$c_label <- paste0(
  "c*=",
  make_pattern_label(C_star_selected)
)

draw_ids_density <- seq_len(dim(main_fit$mcmc$samples_U)[1])
if (length(draw_ids_density) > n_density_draw) {
  draw_ids_density <- sample(draw_ids_density, n_density_draw)
}

n_density_eff <- length(draw_ids_density)

draws_eiv_selected <- sample_eiv_test_y_ordprobit_fb(
  X_test_raw = X_star_selected,
  C_test = C_star_selected,
  fit_obj = main_fit,
  draw_ids = draw_ids_density,
  n_per_draw = 1L,
  n_new_latent_gibbs = n_new_latent_gibbs
)

draws_baseline_selected <- predict_embedding_baseline_samples_ord(
  baselines = rep_baselines,
  X_star_raw = X_star_selected,
  C_star = C_star_selected,
  n_draw = n_density_eff
)

draws_oracle_selected <- sample_oracle_test_y_2d(
  X_test = X_star_selected,
  C_test = C_star_selected,
  true_params = rep_dat$true_params,
  sigma_eps = rep_dat$sigma_eps,
  n_draw = n_density_eff,
  oracle_pool = oracle_pool
)

density_parts <- list(
  make_long_draws(draws_oracle_selected, "Oracle", selected_grid),
  make_long_draws(draws_eiv_selected, "EIV-GP", selected_grid),
  make_long_draws(draws_baseline_selected[["GP-LearnedEmb"]], "GP-LearnedEmb", selected_grid),
  make_long_draws(draws_baseline_selected[["GP-CondMean"]], "GP-CondMean", selected_grid),
  make_long_draws(draws_baseline_selected[["GP-Gaussian"]], "GP-Gaussian", selected_grid)
)

df_density <- dplyr::bind_rows(density_parts)
df_density$method <- factor(df_density$method, levels = method_levels)
df_density$x_label <- factor(df_density$x_label, levels = unique(selected_grid$x_label))
df_density$c_label <- factor(df_density$c_label, levels = unique(selected_grid$c_label))

p_density <- ggplot(df_density, aes(x = y, color = method)) +
  geom_density(linewidth = 0.85) +
  facet_grid(c_label ~ x_label, scales = "free_y") +
  scale_color_manual(values = method_cols, name = NULL, drop = FALSE) +
  labs(
    x = expression(y^"*"),
    y = "Density",
    title = "Study II: predictive distributions at selected mixed inputs"
  ) +
  theme(legend.position = "bottom")

save_plot(
  file.path(FIG_DIR, "fig3_study2_predictive_densities_selected"),
  p_density,
  width = 11,
  height = 7.5
)

############################################################
## MCMC traces and summary table
############################################################

df_trace <- make_trace_df(main_fit, main_calib)

trace_param_cols <- setdiff(
  names(df_trace),
  c("chain", "draw", "n_calib", "n_calib_label")
)

df_trace_long <- df_trace |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(trace_param_cols),
    names_to = "parameter",
    values_to = "value"
  )

p_trace <- ggplot(df_trace_long, aes(x = draw, y = value, color = chain)) +
  geom_line(linewidth = 0.35, alpha = 0.75) +
  facet_wrap(~parameter, scales = "free_y", ncol = 3) +
  labs(
    x = "Saved draw within chain",
    y = NULL,
    color = "Chain",
    title = paste0("Study II: exact MCMC trace plots, |O| = ", main_calib)
  ) +
  theme(legend.position = "bottom")

save_plot(
  file.path(FIG_DIR, "fig5_study2_mcmc_traces"),
  p_trace,
  width = 11,
  height = 6.5
)

mcmc_summary <- main_fit$diagnostics$summary

writeLines(
  knitr::kable(mcmc_summary, format = "latex", booktabs = TRUE, digits = 3),
  con = file.path(TAB_DIR, "study2_mcmc_summary.tex")
)

save_csv(
  mcmc_summary,
  file.path(RES_DIR, "representative_mcmc_summary_main_calibration.csv")
)

cat("\nRepresentative Study II figures written to:\n")
cat(normalizePath(FIG_DIR), "\n")
cat("\nRepresentative Study II tables written to:\n")
cat(normalizePath(TAB_DIR), "\n")
cat("\nRepresentative Study II results written to:\n")
cat(normalizePath(RES_DIR), "\n")