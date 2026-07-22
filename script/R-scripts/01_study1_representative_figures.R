############################################################
## 01_study1_representative_figures.R
##
## Representative-data figures and diagnostics for revised Study I.
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

rep_n_iter <- if (STUDY1_QUICK) 900L else 6000L
rep_burn <- if (STUDY1_QUICK) 300L else 2000L
rep_n_chains <- if (STUDY1_QUICK) 2L else 12L
rep_preset <- if (STUDY1_QUICK) "fast" else "thorough"

n_pred_draw <- if (STUDY1_QUICK) 150L else 600L
n_density_draw <- if (STUDY1_QUICK) 300L else 1200L

parallel_chains <- (
  .Platform$OS.type != "windows" &&
    parallel::detectCores(logical = TRUE) > 2L
)

class_cols <- setNames(
  c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"),
  as.character(seq_len(m))
)

method_cols <- c(
  "Oracle" = "black",
  "EIV-GP" = "firebrick",
  "GP-LearnedEmb" = "purple4",
  "GP-CondMean" = "steelblue",
  "GP-Gaussian" = "darkgreen"
)

############################################################
## Simulate representative data
############################################################

rep_data_file <- file.path(RES_DIR, "representative_data.rds")
rep_fits_file <- file.path(RES_DIR, "representative_eiv_fits.rds")
rep_base_file <- file.path(RES_DIR, "representative_baselines.rds")

if (STUDY1_USE_CACHE && file.exists(rep_data_file)) {
  rep_dat <- readRDS(rep_data_file)
} else {
  rep_dat <- simulate_1d_data(
    n = n_train,
    n_test = n_test,
    m = m,
    scenario = "active",
    sigma_eps = 0.1,
    seed = 20260705
  )
  
  saveRDS(rep_dat, rep_data_file)
}

train <- rep_dat$train
test <- rep_dat$test
tau_true <- rep_dat$tau_true

calib_sets <- make_nested_calibration_sets(
  n = n_train,
  calib_grid = calib_grid,
  seed = 20260706
)

############################################################
## Fit EIV-GP for representative calibration sizes
############################################################

rep_fit_calibs <- c(0L, 10L, 50L)

if (STUDY1_USE_CACHE && file.exists(rep_fits_file)) {
  rep_fits <- readRDS(rep_fits_file)
} else {
  rep_fits <- list()
  
  for (kk in rep_fit_calibs) {
    cat("Fitting representative EIV-GP with |O| =", kk, "\n")
    
    rep_fits[[as.character(kk)]] <- fit_eivgp_1d(
      x_raw = train$x,
      y_raw = train$y,
      c_ord = train$c,
      u_true = train$u,
      calib_idx = calib_sets[[as.character(kk)]],
      m = m,
      tau_true = tau_true,
      n_iter = rep_n_iter,
      burn = rep_burn,
      thin = 1L,
      n_chains = rep_n_chains,
      preset = rep_preset,
      seed = 500000L + kk,
      parallel_chains = parallel_chains,
      verbose = TRUE
    )
  }
  
  saveRDS(rep_fits, rep_fits_file)
}

main_fit <- rep_fits[["10"]]

############################################################
## Fit deterministic embedding baselines
############################################################

if (STUDY1_USE_CACHE && file.exists(rep_base_file)) {
  rep_baselines <- readRDS(rep_base_file)
} else {
  rep_baselines <- fit_embedding_baselines(
    x_raw = train$x,
    y_raw = train$y,
    c_ord = train$c,
    m = m,
    n_starts_learned = if (STUDY1_QUICK) 3L else 8L
  )
  
  saveRDS(rep_baselines, rep_base_file)
}

############################################################
## Figure 1: representative data
############################################################

df_train <- data.frame(
  x = train$x,
  u = train$u,
  c = factor(train$c),
  y = train$y
)

p_latent <- ggplot(df_train, aes(x = u, y = y, color = c)) +
  geom_point(size = 2, alpha = 0.85) +
  geom_vline(xintercept = tau_true, color = "black", linewidth = 0.5) +
  scale_color_manual(values = class_cols, name = "Class") +
  labs(
    x = "Latent u",
    y = "y",
    title = "Training data in latent coordinate"
  ) +
  theme(legend.position = "right")

p_x <- ggplot(df_train, aes(x = x, y = y, color = c)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_manual(values = class_cols, name = "Class") +
  labs(
    x = "Observed quantitative input x",
    y = "y",
    title = "Response variation over observed x"
  ) +
  theme(legend.position = "right")

p_data_active <- p_latent + p_x + patchwork::plot_layout(ncol = 2)

ggsave(
  file.path(FIG_DIR, "fig1_study1_data_active_x.pdf"),
  p_data_active,
  width = 10,
  height = 4.5
)

############################################################
## Figure 2: latent imputation across calibration sizes
############################################################

extract_imputation_df <- function(fit, n_calib) {
  samples_u <- fit$mcmc$samples_u
  miss_idx <- fit$data$miss_idx
  
  post_mean <- colMeans(samples_u)[miss_idx]
  post_lo <- apply(samples_u[, miss_idx, drop = FALSE], 2, quantile, probs = 0.025)
  post_hi <- apply(samples_u[, miss_idx, drop = FALSE], 2, quantile, probs = 0.975)
  
  true_u <- fit$data$u_true[miss_idx]
  
  data.frame(
    n_calib = n_calib,
    id = miss_idx,
    true_u = true_u,
    post_mean = post_mean,
    post_lo = post_lo,
    post_hi = post_hi,
    c = factor(fit$data$c_ord[miss_idx]),
    covered = true_u >= post_lo & true_u <= post_hi
  )
}

df_imp_all <- dplyr::bind_rows(
  lapply(names(rep_fits), function(nm) {
    extract_imputation_df(rep_fits[[nm]], as.integer(nm))
  })
)

df_imp_all$n_calib_label <- factor(
  paste0("|O| = ", df_imp_all$n_calib),
  levels = paste0("|O| = ", rep_fit_calibs)
)

p_imp_calib <- ggplot(df_imp_all, aes(x = true_u, y = post_mean, color = c)) +
  geom_errorbar(aes(ymin = post_lo, ymax = post_hi), width = 0, alpha = 0.18) +
  geom_point(size = 1.9, alpha = 0.85) +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 0.7) +
  facet_wrap(~n_calib_label, nrow = 1) +
  scale_color_manual(values = class_cols, name = "Class") +
  labs(
    x = "True latent u",
    y = "Posterior mean and 95% interval",
    title = "Latent imputation across calibration sizes"
  ) +
  theme(legend.position = "bottom")

ggsave(
  file.path(FIG_DIR, "fig2_study1_latent_imputation_by_calibration.pdf"),
  p_imp_calib,
  width = 11,
  height = 4.5
)

############################################################
## Figure 3: latent response-surface slices
############################################################

make_eiv_function_slice_df <- function(fit,
                                       x_slices = c(-1, 0, 1),
                                       u_grid = seq(-2.4, 2.4, length.out = 160),
                                       draw_ids = NULL,
                                       scenario = "active",
                                       label = "EIV-GP") {
  if (is.null(draw_ids)) {
    draw_ids <- seq_len(nrow(fit$mcmc$samples_u))
  }
  
  if (length(draw_ids) > n_pred_draw) {
    draw_ids <- sample(draw_ids, n_pred_draw)
  }
  
  grid <- expand.grid(
    x_raw = x_slices,
    u = u_grid
  )
  
  x_star <- as.numeric((grid$x_raw - fit$data$x_center) / fit$data$x_scale)
  u_star <- grid$u
  
  f_samps <- matrix(NA_real_, nrow = length(draw_ids), ncol = nrow(grid))
  
  for (ii in seq_along(draw_ids)) {
    s <- draw_ids[ii]
    
    pred <- gp_predict_draw(
      x_train = fit$data$x,
      u_train = fit$mcmc$samples_u[s, ],
      y_train = fit$data$y,
      x_star = x_star,
      u_star = u_star,
      logtheta = fit$mcmc$samples_logtheta[s, ],
      sigma2_eps = fit$mcmc$samples_sigma2[s],
      noisy = FALSE
    )
    
    f_std <- pred$mean + sqrt(pred$var) * rnorm(nrow(grid))
    f_samps[ii, ] <- fit$data$y_center + fit$data$y_scale * f_std
  }
  
  out <- grid
  out$mean <- colMeans(f_samps)
  out$lo <- apply(f_samps, 2, quantile, probs = 0.025)
  out$hi <- apply(f_samps, 2, quantile, probs = 0.975)
  out$truth <- f0_1d(out$x_raw, out$u, scenario = scenario)
  out$method <- label
  out$x_slice <- factor(paste0("x = ", out$x_raw), levels = paste0("x = ", x_slices))
  
  out
}

make_cc_function_slice_df <- function(fit,
                                      x_slices = c(-1, 0, 1),
                                      u_grid = seq(-2.4, 2.4, length.out = 160),
                                      scenario = "active",
                                      label = "Complete-case GP") {
  calib_idx <- fit$data$calib_idx
  
  if (length(calib_idx) < 3) {
    return(data.frame())
  }
  
  fit_cc <- gp_mle_fit(
    X = cbind(fit$data$x[calib_idx], fit$data$u_true[calib_idx]),
    y = fit$data$y[calib_idx]
  )
  
  grid <- expand.grid(
    x_raw = x_slices,
    u = u_grid
  )
  
  x_star <- as.numeric((grid$x_raw - fit$data$x_center) / fit$data$x_scale)
  
  pred <- gp_mle_predict(
    fit_cc,
    Xstar = cbind(x_star, grid$u),
    noisy = FALSE
  )
  
  out <- grid
  out$mean <- fit$data$y_center + fit$data$y_scale * pred$mean
  out$lo <- fit$data$y_center + fit$data$y_scale * (pred$mean - 1.96 * sqrt(pred$var))
  out$hi <- fit$data$y_center + fit$data$y_scale * (pred$mean + 1.96 * sqrt(pred$var))
  out$truth <- f0_1d(out$x_raw, out$u, scenario = scenario)
  out$method <- label
  out$x_slice <- factor(paste0("x = ", out$x_raw), levels = paste0("x = ", x_slices))
  
  out
}

df_fun_eiv <- make_eiv_function_slice_df(main_fit, scenario = rep_dat$scenario)
df_fun_cc <- make_cc_function_slice_df(main_fit, scenario = rep_dat$scenario)

df_fun <- dplyr::bind_rows(df_fun_cc, df_fun_eiv)
df_fun$method <- factor(df_fun$method, levels = c("Complete-case GP", "EIV-GP"))

p_fun_slices <- ggplot(df_fun, aes(x = u)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "skyblue", alpha = 0.45) +
  geom_line(aes(y = mean), color = "blue", linewidth = 0.85, linetype = "dashed") +
  geom_line(aes(y = truth), color = "orange", linewidth = 1.05) +
  facet_grid(method ~ x_slice) +
  labs(
    x = "Latent u",
    y = expression(f(x, u)),
    title = "Recovery of latent response-surface slices"
  )

ggsave(
  file.path(FIG_DIR, "fig3_study1_function_slices.pdf"),
  p_fun_slices,
  width = 11,
  height = 6.5
)

############################################################
## Figure 4: selected mixed-input predictive densities
############################################################

selected_grid <- expand.grid(
  x_star = c(-1, 0, 1),
  c_star = c(1, 3, 6)
)

draw_ids_density <- seq_len(nrow(main_fit$mcmc$samples_u))

if (length(draw_ids_density) > n_density_draw) {
  draw_ids_density <- sample(draw_ids_density, n_density_draw)
}

draws_eiv_selected <- sample_eiv_test_y(
  x_test_raw = selected_grid$x_star,
  c_test = selected_grid$c_star,
  fit_obj = main_fit,
  draw_ids = draw_ids_density,
  n_per_draw = 1L
)

draws_baseline_selected <- predict_embedding_baseline_samples(
  baselines = rep_baselines,
  x_star_raw = selected_grid$x_star,
  c_star = selected_grid$c_star,
  m = m,
  n_draw = n_density_draw
)

draws_oracle_selected <- sample_oracle_test_y(
  x_test = selected_grid$x_star,
  c_test = selected_grid$c_star,
  tau_true = tau_true,
  scenario = rep_dat$scenario,
  sigma_eps = rep_dat$sigma_eps,
  n_draw = n_density_draw
)

make_long_draws <- function(draw_mat, method, grid) {
  do.call(
    rbind,
    lapply(seq_len(ncol(draw_mat)), function(j) {
      data.frame(
        x_star = grid$x_star[j],
        c_star = grid$c_star[j],
        y = draw_mat[, j],
        method = method
      )
    })
  )
}

df_density <- dplyr::bind_rows(
  make_long_draws(draws_oracle_selected, "Oracle", selected_grid),
  make_long_draws(draws_eiv_selected, "EIV-GP", selected_grid),
  make_long_draws(draws_baseline_selected$`GP-Gaussian`, "GP-Gaussian", selected_grid),
  make_long_draws(draws_baseline_selected$`GP-CondMean`, "GP-CondMean", selected_grid),
  make_long_draws(draws_baseline_selected$`GP-LearnedEmb`, "GP-LearnedEmb", selected_grid)
)

df_density$method <- factor(
  df_density$method,
  levels = c("Oracle", "EIV-GP", "GP-LearnedEmb", "GP-CondMean", "GP-Gaussian")
)

df_density$x_label <- factor(
  paste0("x* = ", df_density$x_star),
  levels = paste0("x* = ", c(-1, 0, 1))
)

df_density$c_label <- factor(
  paste0("c* = ", df_density$c_star),
  levels = paste0("c* = ", c(1, 3, 6))
)

p_density <- ggplot(df_density, aes(x = y, color = method)) +
  geom_density(linewidth = 0.85) +
  facet_grid(c_label ~ x_label, scales = "free_y") +
  scale_color_manual(values = method_cols, name = NULL) +
  labs(
    x = expression(y^"*"),
    y = "Density",
    title = "Predictive distributions at selected mixed inputs"
  ) +
  theme(legend.position = "bottom")

ggsave(
  file.path(FIG_DIR, "fig4_study1_predictive_densities_selected.pdf"),
  p_density,
  width = 11,
  height = 7.5
)

############################################################
## Figure 6: MCMC trace plots for representative fit
############################################################

samples_by_chain <- main_fit$mcmc$samples_by_chain

df_trace <- dplyr::bind_rows(
  lapply(seq_along(samples_by_chain$u), function(cc) {
    data.frame(
      chain = factor(cc),
      draw = seq_along(samples_by_chain$sigma2[[cc]]),
      sigma_epsilon = sqrt(samples_by_chain$sigma2[[cc]]),
      rho = exp(samples_by_chain$logtheta[[cc]][, 1]),
      theta_x = exp(samples_by_chain$logtheta[[cc]][, 2]),
      theta_u = exp(samples_by_chain$logtheta[[cc]][, 3])
    )
  })
)

df_trace_long <- df_trace |>
  tidyr::pivot_longer(
    cols = c(sigma_epsilon, rho, theta_x, theta_u),
    names_to = "parameter",
    values_to = "value"
  )

p_trace <- ggplot(df_trace_long, aes(x = draw, y = value, color = chain)) +
  geom_line(linewidth = 0.35, alpha = 0.75) +
  facet_wrap(~parameter, scales = "free_y", ncol = 2) +
  labs(
    x = "Saved draw within chain",
    y = NULL,
    color = "Chain",
    title = "MCMC trace plots for representative EIV-GP fit"
  ) +
  theme(legend.position = "bottom")

ggsave(
  file.path(FIG_DIR, "fig6_study1_mcmc_traces.pdf"),
  p_trace,
  width = 10,
  height = 6
)

############################################################
## Tables: MCMC summary and class counts
############################################################

mcmc_summary <- main_fit$diagnostics$summary

mcmc_table <- mcmc_summary |>
  dplyr::mutate(
    max_rhat_hyper = sprintf("%.3f", max_rhat_hyper),
    max_rhat_tau = sprintf("%.3f", max_rhat_tau),
    median_rhat_missing_u = sprintf("%.3f", median_rhat_missing_u),
    max_rhat_missing_u = sprintf("%.3f", max_rhat_missing_u),
    min_ess_key = sprintf("%.1f", min_ess_key),
    time_seconds = sprintf("%.1f", time_seconds)
  )

latex_mcmc_table <- knitr::kable(
  mcmc_table,
  format = "latex",
  booktabs = TRUE,
  escape = TRUE
)

writeLines(
  latex_mcmc_table,
  con = file.path(TAB_DIR, "study1_mcmc_summary.tex")
)

class_counts <- data.frame(
  class = seq_len(m),
  count = as.integer(tabulate(train$c, nbins = m)),
  calibrated_10 = as.integer(tabulate(
    train$c[calib_sets[["10"]]],
    nbins = m
  ))
)

latex_counts_table <- knitr::kable(
  class_counts,
  format = "latex",
  booktabs = TRUE,
  align = "ccc",
  col.names = c("Ordinal level", "Training count", "Calibrated count, |O|=10")
)

writeLines(
  latex_counts_table,
  con = file.path(TAB_DIR, "study1_class_counts.tex")
)

cat("\nRepresentative Study I figures written to:\n")
cat(normalizePath(FIG_DIR), "\n")
cat("\nRepresentative Study I tables written to:\n")
cat(normalizePath(TAB_DIR), "\n")