############################################################
## 01_ocean_representative_figures.R
##
## Representative-data figures for the ocean N/P role-swap case.
## Same role as 01_study1_representative_figures.R, but x is
## (temperature, salinity) and u is standardized log(phosphate).
## Do not reuse Study I slice / oracle helpers: those assume 1d x
## and a known f0.
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
  c("ggplot2", "patchwork", "dplyr", "tidyr", "knitr"),
  "ocean representative analysis"
)

FIG_DIR <- file.path(OCEAN_OUT_PREFIX, "figures", "ocean")
TAB_DIR <- file.path(OCEAN_OUT_PREFIX, "tables", "ocean")
RES_DIR <- file.path(OCEAN_OUT_PREFIX, "results", "ocean")

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)

ocean <- load_ocean_prepared()
train <- ocean$train
test <- ocean$test
m <- ocean$m
tau_ref <- ocean$tau_reference

calib_grid <- c(10L, 25L, 50L)
rep_fit_calibs <- if (OCEAN_QUICK) 10L else c(10L, 25L, 50L)

rep_n_iter <- if (OCEAN_QUICK) 400L else 20000L
rep_burn <- if (OCEAN_QUICK) 100L else 5000L
rep_n_chains <- if (OCEAN_QUICK) 1L else 4L
rep_preset <- if (OCEAN_QUICK) "fast" else "thorough"
n_pred_draw <- if (OCEAN_QUICK) 80L else 600L

parallel_chains <- (
  .Platform$OS.type != "windows" &&
    parallel::detectCores(logical = TRUE) > 2L
)

class_cols <- setNames(
  c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"),
  as.character(seq_len(m))
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

alloc <- read_ocean_class_allocation()
calib_sets <- make_ocean_nested_calibration(
  c_train = ocean$c_train,
  allocation = alloc,
  calib_grid = calib_grid,
  seed = 20260805L
)

############################################################
## Fit EIV-GP for representative calibration sizes
############################################################

cache_tag <- paste(
  OCEAN_CACHE_VERSION,
  ifelse(OCEAN_QUICK, "quick", "paper"),
  normalize_gp_kernel_1d(OCEAN_KERNEL, OCEAN_MATERN_NU)$name,
  sep = "_"
)
rep_fits_file <- file.path(RES_DIR, paste0("representative_eiv_fits_", cache_tag, ".rds"))
rep_base_file <- file.path(RES_DIR, paste0("representative_baselines_", cache_tag, ".rds"))

if (OCEAN_USE_CACHE && file.exists(rep_fits_file)) {
  rep_fits <- readRDS(rep_fits_file)
} else {
  rep_fits <- list()
  for (kk in rep_fit_calibs) {
    cat("Fitting ocean EIV-GP with |O| =", kk, "\n")
    u_obs_kk <- rep(NA_real_, length(ocean$u_train))
    u_obs_kk[calib_sets[[as.character(kk)]]] <-
      ocean$u_train[calib_sets[[as.character(kk)]]]
    rep_fits[[as.character(kk)]] <- fit_eivgp_1d(
      x_raw = ocean$X_train,
      y_raw = ocean$y_train,
      c_ord = ocean$c_train,
      calib_idx = calib_sets[[as.character(kk)]],
      m = m,
      tau_true = tau_ref,
      n_iter = rep_n_iter,
      burn = rep_burn,
      thin = 1L,
      n_chains = rep_n_chains,
      preset = rep_preset,
      seed = 500000L + kk,
      parallel_chains = parallel_chains,
      verbose = TRUE,
      kernel = OCEAN_KERNEL,
      matern_nu = OCEAN_MATERN_NU,
      u_obs = u_obs_kk
    )
  }
  saveRDS(rep_fits, rep_fits_file)
}

main_key <- as.character(max(as.integer(names(rep_fits))))
main_fit <- rep_fits[[main_key]]

############################################################
## Embedding baselines (same names as Study I)
############################################################

if (OCEAN_USE_CACHE && file.exists(rep_base_file)) {
  rep_baselines <- readRDS(rep_base_file)
} else {
  rep_baselines <- fit_embedding_baselines(
    x_raw = ocean$X_train,
    y_raw = ocean$y_train,
    c_ord = ocean$c_train,
    m = m,
    n_starts_learned = if (OCEAN_QUICK) 3L else 8L,
    kernel = OCEAN_KERNEL,
    matern_nu = OCEAN_MATERN_NU
  )
  saveRDS(rep_baselines, rep_base_file)
}

############################################################
## Figure 1: representative data (latent u and temperature)
############################################################

df_train <- data.frame(
  u = ocean$u_train,
  x_temperature = train$x_temperature,
  x_salinity = train$x_salinity,
  c = factor(ocean$c_train),
  y = ocean$y_train
)

p_latent <- ggplot(df_train, aes(x = u, y = y, color = c)) +
  geom_vline(xintercept = tau_ref, color = "black", linewidth = 0.5) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_manual(values = class_cols, name = "Class") +
  labs(
    x = "Latent u = standardized log(phosphate)",
    y = "y = nitrate+nitrite (uM)",
    title = "Training data in latent coordinate"
  ) +
  theme(legend.position = "right")

p_x <- ggplot(df_train, aes(x = x_temperature, y = y, color = c)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_manual(values = class_cols, name = "Class") +
  labs(
    x = "Temperature (°C)",
    y = "y = nitrate+nitrite (uM)",
    title = "Response variation over temperature"
  ) +
  theme(legend.position = "right")

ggsave(
  file.path(FIG_DIR, "fig1_ocean_data_u_and_temperature.pdf"),
  p_latent + p_x + patchwork::plot_layout(ncol = 2),
  width = 10,
  height = 4.5
)

############################################################
## Figure 2: latent imputation across calibration sizes
############################################################

df_imp_all <- dplyr::bind_rows(
  lapply(names(rep_fits), function(nm) {
    ans <- summarize_ocean_u_imputation(
      fit = rep_fits[[nm]],
      u_true = ocean$u_train,
      c_ord = ocean$c_train,
      inverse_u = ocean$inverse_u
    )
    ans$n_calib <- as.integer(nm)
    ans
  })
)

df_imp_all$n_calib_label <- factor(
  paste0("|O| = ", df_imp_all$n_calib),
  levels = paste0("|O| = ", sort(unique(df_imp_all$n_calib)))
)
df_imp_all$c <- factor(df_imp_all$c)

p_imp <- ggplot(df_imp_all, aes(x = true_u, y = post_mean, color = c)) +
  geom_errorbar(aes(ymin = post_lo, ymax = post_hi), width = 0, alpha = 0.18) +
  geom_point(size = 1.9, alpha = 0.85) +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 0.7) +
  facet_wrap(~n_calib_label, nrow = 1) +
  scale_color_manual(values = class_cols, name = "Class") +
  labs(
    x = "True standardized log(phosphate)",
    y = "Posterior mean and 95% interval",
    title = "Hidden-phosphate imputation across calibration sizes"
  ) +
  theme(legend.position = "bottom")

ggsave(
  file.path(FIG_DIR, "fig2_ocean_latent_imputation_by_calibration.pdf"),
  p_imp,
  width = 11,
  height = 4.5
)

############################################################
## Figure 3: test-set predictive means (no synthetic f0)
############################################################

set.seed(20260830)
ids <- seq_len(nrow(main_fit$mcmc$samples_u))
if (length(ids) > n_pred_draw) ids <- sort(sample(ids, n_pred_draw))

eiv_draws <- sample_eiv_test_y(
  x_test_raw = ocean$X_test,
  c_test = ocean$c_test,
  fit_obj = main_fit,
  draw_ids = ids,
  n_per_draw = 1L
)
base_draws <- predict_embedding_baseline_samples(
  baselines = rep_baselines,
  x_star_raw = ocean$X_test,
  c_star = ocean$c_test,
  m = m,
  n_draw = n_pred_draw
)

X_train_std <- sweep(
  sweep(ocean$X_train, 2L, rep_baselines$x_center, "-"),
  2L, rep_baselines$x_scale, "/"
)
X_test_std <- sweep(
  sweep(ocean$X_test, 2L, rep_baselines$x_center, "-"),
  2L, rep_baselines$x_scale, "/"
)
y_train_std <- as.numeric(
  (ocean$y_train - rep_baselines$y_center) / rep_baselines$y_scale
)

set.seed(20260831)
x_only_fit <- gp_mle_fit_1d(
  X_train_std,
  y_train_std,
  kernel = OCEAN_KERNEL,
  matern_nu = OCEAN_MATERN_NU
)
x_only_draws <- rep_baselines$y_center + rep_baselines$y_scale *
  sample_gp_mle_predictive(x_only_fit, X_test_std, n_draw = n_pred_draw)

full_u_fit <- fit_ocean_latent_gp(
  ocean$X_train,
  ocean$y_train,
  ocean$u_train,
  kernel = OCEAN_KERNEL,
  matern_nu = OCEAN_MATERN_NU
)
full_u_draws <- sample_ocean_latent_gp_y_given_xc(
  full_u_fit,
  ocean$X_test,
  ocean$c_test,
  tau = tau_ref,
  n_draw = n_pred_draw,
  seed = 20260832L
)

cc_fit <- fit_ocean_latent_gp(
  ocean$X_train,
  ocean$y_train,
  ocean$u_train,
  train_idx = main_fit$data$calib_idx,
  kernel = OCEAN_KERNEL,
  matern_nu = OCEAN_MATERN_NU
)
cc_draws <- sample_ocean_latent_gp_y_given_xc(
  cc_fit,
  ocean$X_test,
  ocean$c_test,
  tau = estimate_tau_from_ordinal_codes(ocean$c_train, m),
  n_draw = n_pred_draw,
  seed = 20260833L
)

pred_df <- data.frame(
  y_true = ocean$y_test,
  c = factor(ocean$c_test),
  eiv = colMeans(eiv_draws),
  learned = colMeans(base_draws[["GP-LearnedEmb"]]),
  condmean = colMeans(base_draws[["GP-CondMean"]])
)

p_pred <- ggplot(pred_df, aes(x = y_true, y = eiv, color = c)) +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 0.7) +
  geom_point(size = 1.7, alpha = 0.8) +
  scale_color_manual(values = class_cols, name = "Class") +
  labs(
    x = "True test y (uM)",
    y = paste0("EIV-GP posterior mean, |O| = ", main_key),
    title = "Test-set prediction (ocean, multivariate x)"
  ) +
  theme(legend.position = "bottom")

ggsave(
  file.path(FIG_DIR, "fig3_ocean_test_prediction.pdf"),
  p_pred,
  width = 6.5,
  height = 5.2
)

############################################################
## Figure 4: MCMC traces (p = 2: theta_x1, theta_x2, theta_u)
############################################################

samples_by_chain <- main_fit$mcmc$samples_by_chain
theta_names <- paste0("logtheta", seq_len(ncol(samples_by_chain$logtheta[[1]])))

df_trace <- dplyr::bind_rows(
  lapply(seq_along(samples_by_chain$u), function(cc) {
    lt <- samples_by_chain$logtheta[[cc]]
    data.frame(
      chain = factor(cc),
      draw = seq_len(nrow(lt)),
      sigma_epsilon = sqrt(samples_by_chain$sigma2[[cc]]),
      rho = exp(lt[, 1]),
      theta_x1 = exp(lt[, 2]),
      theta_x2 = exp(lt[, 3]),
      theta_u = exp(lt[, ncol(lt)])
    )
  })
)

df_trace_long <- df_trace |>
  tidyr::pivot_longer(
    cols = c(sigma_epsilon, rho, theta_x1, theta_x2, theta_u),
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
    title = paste0("MCMC traces for ocean EIV-GP, |O| = ", main_key)
  ) +
  theme(legend.position = "bottom")

ggsave(
  file.path(FIG_DIR, "fig4_ocean_mcmc_traces.pdf"),
  p_trace,
  width = 10,
  height = 7
)

############################################################
## Small prediction table for this representative split
############################################################

metrics <- rbind(
  summarize_predictive_samples_1d(
    eiv_draws, ocean$y_test, method = paste0("EIV-GP |O|=", main_key),
    rep_id = 0L, n_calib = as.integer(main_key), scenario = "ocean"
  ),
  summarize_predictive_samples_1d(
    base_draws[["GP-LearnedEmb"]], ocean$y_test, method = "GP-LearnedEmb",
    rep_id = 0L, n_calib = NA_integer_, scenario = "ocean"
  ),
  summarize_predictive_samples_1d(
    base_draws[["GP-CondMean"]], ocean$y_test, method = "GP-CondMean",
    rep_id = 0L, n_calib = NA_integer_, scenario = "ocean"
  ),
  summarize_predictive_samples_1d(
    base_draws[["GP-Gaussian"]], ocean$y_test, method = "GP-Gaussian",
    rep_id = 0L, n_calib = NA_integer_, scenario = "ocean"
  ),
  summarize_predictive_samples_1d(
    x_only_draws, ocean$y_test, method = "GP-XOnly",
    rep_id = 0L, n_calib = NA_integer_, scenario = "ocean"
  ),
  summarize_predictive_samples_1d(
    cc_draws, ocean$y_test, method = "Complete-case GP",
    rep_id = 0L, n_calib = as.integer(main_key), scenario = "ocean"
  ),
  summarize_predictive_samples_1d(
    full_u_draws, ocean$y_test, method = "Full-U GP oracle",
    rep_id = 0L, n_calib = NA_integer_, scenario = "ocean"
  )
)
write.csv(metrics, file.path(TAB_DIR, "ocean_representative_prediction.csv"), row.names = FALSE)

cat("\nOcean representative figures written to:\n")
cat(normalizePath(FIG_DIR), "\n")
