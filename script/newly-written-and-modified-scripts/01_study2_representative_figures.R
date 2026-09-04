############################################################
## 01_study2_representative_figures.R
##
## Representative-data analysis for the frozen Study II design.
##
## This script produces only figures and tables used by the manuscript.  It
## deliberately excludes the home-built embedding methods used in early pilot
## runs.  Published competitors enter only through audited adapters in
## competitors.R, with failures recorded explicitly.
############################################################

if (!exists("fit_eivgp_ordprobit_fb")) {
  source("model_multivariate.R")
}
if (!exists("run_study2_published_competitors")) {
  source("competitors.R")
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
})

if (!exists("STUDY2_CONFIG")) STUDY2_CONFIG <- "quick"
if (!exists("STUDY2_USE_CACHE")) STUDY2_USE_CACHE <- TRUE
if (!exists("STUDY2_OUT_PREFIX")) STUDY2_OUT_PREFIX <- ".."
if (!exists("STUDY2_SAVE_PDF")) STUDY2_SAVE_PDF <- TRUE
if (!exists("STUDY2_SAVE_PNG")) STUDY2_SAVE_PNG <- TRUE
if (!exists("STUDY2_STRICT_COMPETITORS")) STUDY2_STRICT_COMPETITORS <- FALSE
if (!exists("STUDY2_REP_SCENARIO")) STUDY2_REP_SCENARIO <- "primary"
if (!exists("STUDY2_REP_FIT_CALIBS")) {
  STUDY2_REP_FIT_CALIBS <- c(0L, 10L, 25L, 50L, 80L)
}
if (!exists("STUDY2_MAIN_CALIB")) STUDY2_MAIN_CALIB <- 25L
if (!exists("STUDY2_PUBLISHED_COMPETITORS")) {
  STUDY2_PUBLISHED_COMPETITORS <- c("UC-GP", "LVGP", "EzGP")
}

settings <- study2_config_settings(STUDY2_CONFIG)
study2_scenario_spec(STUDY2_REP_SCENARIO)

n_train <- 120L
n_test <- settings$n_test
m_vec <- rep(4L, 4L)
d_latent <- 2L
rep_fit_calibs <- sort(unique(as.integer(STUDY2_REP_FIT_CALIBS)))
if (any(rep_fit_calibs < 0L | rep_fit_calibs > n_train)) {
  stop("Every representative calibration size must lie between 0 and 120.")
}
if (!(STUDY2_MAIN_CALIB %in% rep_fit_calibs)) {
  stop("STUDY2_MAIN_CALIB must be included in STUDY2_REP_FIT_CALIBS.")
}

rep_n_iter <- if (exists("STUDY2_REP_N_ITER")) {
  as.integer(STUDY2_REP_N_ITER)
} else settings$rep_n_iter
rep_burn <- if (exists("STUDY2_REP_BURN")) {
  as.integer(STUDY2_REP_BURN)
} else settings$rep_burn
rep_thin <- if (exists("STUDY2_REP_THIN")) {
  as.integer(STUDY2_REP_THIN)
} else settings$rep_thin
rep_n_chains <- if (exists("STUDY2_REP_N_CHAINS")) {
  as.integer(STUDY2_REP_N_CHAINS)
} else settings$rep_n_chains

n_pred_draw <- settings$n_pred_draw
predictive_latent_sampler <- if (exists("STUDY2_PREDICTIVE_LATENT_SAMPLER")) {
  as.character(STUDY2_PREDICTIVE_LATENT_SAMPLER)
} else settings$predictive_latent_sampler
if (length(predictive_latent_sampler) != 1L ||
    is.na(predictive_latent_sampler) ||
    !predictive_latent_sampler %in%
      c("minimax_tilting", "rejection", "gibbs")) {
  stop(
    "STUDY2_PREDICTIVE_LATENT_SAMPLER must be minimax_tilting, rejection, ",
    "or gibbs."
  )
}
diagnostic_n_new_latent_gibbs <- if (exists("STUDY2_DIAGNOSTIC_GIBBS_SWEEPS")) {
  as.integer(STUDY2_DIAGNOSTIC_GIBBS_SWEEPS)
} else settings$diagnostic_n_new_latent_gibbs
if (length(diagnostic_n_new_latent_gibbs) != 1L ||
    is.na(diagnostic_n_new_latent_gibbs) ||
    diagnostic_n_new_latent_gibbs < 1L) {
  stop("STUDY2_DIAGNOSTIC_GIBBS_SWEEPS must be a positive integer.")
}
rejection_max_batches <- settings$rejection_max_batches
n_oracle_pool <- settings$n_oracle_pool
study2_parallel_level <- if (exists("STUDY2_PARALLEL_LEVEL")) {
  match.arg(STUDY2_PARALLEL_LEVEL, c("chains", "replications", "none"))
} else {
  "chains"
}
parallel_chains <- mixedgp_parallel_chains_enabled(study2_parallel_level)

FIG_DIR <- file.path(STUDY2_OUT_PREFIX, "figures")
TAB_DIR <- file.path(STUDY2_OUT_PREFIX, "tables")
RES_DIR <- file.path(STUDY2_OUT_PREFIX, "results", "study2_manuscript_v5")
for (dd in c(FIG_DIR, TAB_DIR, RES_DIR)) {
  dir.create(dd, showWarnings = FALSE, recursive = TRUE)
}

calib_tag <- paste(rep_fit_calibs, collapse = "-")
CACHE_TAG <- paste0(
  "representative_", STUDY2_DESIGN_TAG,
  "_", STUDY2_REP_SCENARIO,
  "_", STUDY2_CONFIG,
  "_calib", calib_tag,
  "_iter", rep_n_iter,
  "_chains", rep_n_chains
)

save_plot <- function(path_no_ext, plot, width, height, dpi = 320) {
  if (isTRUE(STUDY2_SAVE_PDF)) {
    ggsave(
      paste0(path_no_ext, ".pdf"), plot,
      width = width, height = height, units = "in", limitsize = FALSE
    )
  }
  if (isTRUE(STUDY2_SAVE_PNG)) {
    ggsave(
      paste0(path_no_ext, ".png"), plot,
      width = width, height = height, units = "in", dpi = dpi,
      limitsize = FALSE, bg = "white"
    )
  }
  invisible(plot)
}

even_draw_ids <- function(n_available, n_keep) {
  if (n_available <= n_keep) return(seq_len(n_available))
  unique(round(seq(1, n_available, length.out = n_keep)))
}

############################################################
## Frozen representative data and calibration sets
############################################################

data_file <- file.path(RES_DIR, paste0(CACHE_TAG, "_data.rds"))
fit_file <- file.path(RES_DIR, paste0(CACHE_TAG, "_eiv_fits.rds"))

if (isTRUE(STUDY2_USE_CACHE) && file.exists(data_file)) {
  rep_dat <- readRDS(data_file)
} else {
  rep_dat <- simulate_study2_data(
    n = n_train,
    n_test = n_test,
    scenario = STUDY2_REP_SCENARIO,
    seed = 20260710L
  )
  saveRDS(rep_dat, data_file)
}

train <- rep_dat$train
test <- rep_dat$test
calib_sets <- make_stratified_calibration_sets_2d(
  C = train$C,
  calib_grid = rep_fit_calibs,
  seed = 20260711L
)

if (isTRUE(STUDY2_USE_CACHE) && file.exists(fit_file)) {
  rep_fits <- readRDS(fit_file)
} else {
  rep_fits <- list()
}

for (n_calib in rep_fit_calibs) {
  key <- as.character(n_calib)
  if (!is.null(rep_fits[[key]])) next

  message("Representative Study II fit: |O|=", n_calib)
  rep_fits[[key]] <- fit_eivgp_ordprobit_fb(
    X_raw = train$X,
    y_raw = train$y,
    C_ord = train$C,
    U_obs = train$U,
    calib_idx = calib_sets[[key]],
    U_true_eval = train$U,
    d = d_latent,
    m_vec = m_vec,
    ident = "lower_triangular",
    n_iter = rep_n_iter,
    burn = rep_burn,
    thin = rep_thin,
    n_chains = rep_n_chains,
    preset = settings$preset,
    sampler_strategy = "interwoven",
    store_scores = FALSE,
    seed = 500000L + n_calib,
    parallel_chains = parallel_chains,
    verbose = TRUE
  )
  saveRDS(rep_fits, fit_file)
}
rep_fits <- rep_fits[as.character(rep_fit_calibs)]
main_fit <- rep_fits[[as.character(STUDY2_MAIN_CALIB)]]

############################################################
## Figure 1: representative design
############################################################

df_train <- data.frame(
  x1 = train$X[, 1],
  x2 = train$X[, 2],
  u1 = train$U[, 1],
  u2 = train$U[, 2],
  y = train$y,
  severity = rowSums(train$C),
  pattern = pattern_key(train$C)
)

p_latent <- ggplot(df_train, aes(u1, u2, color = y)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_viridis_c(name = "Response") +
  labs(
    x = expression(u[1]), y = expression(u[2]),
    title = "Latent coordinates (visualization only)"
  )

p_observed <- ggplot(df_train, aes(x1, x2, color = severity)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_viridis_c(name = "Ordinal\nseverity") +
  labs(
    x = expression(x[1]), y = expression(x[2]),
    title = "Observed quantitative inputs"
  )

p_pattern <- ggplot(df_train, aes(factor(severity), y)) +
  geom_boxplot(outlier.alpha = 0.45) +
  labs(
    x = "Sum of four ordinal levels", y = "Response",
    title = "Response versus observed ordinal summary"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(
  file.path(FIG_DIR, "fig1_study2_data_representative"),
  (p_latent + p_observed) / p_pattern,
  width = 10.5,
  height = 8
)

############################################################
## Figure 2 and table: latent-state imputation
############################################################

imputation_rows <- list()
imputation_metrics <- list()

for (n_calib in rep_fit_calibs[rep_fit_calibs > 0L]) {
  fit <- rep_fits[[as.character(n_calib)]]
  miss_idx <- fit$data$miss_idx
  if (length(miss_idx) == 0L) next

  samples <- fit$mcmc$samples_U[, miss_idx, , drop = FALSE]
  for (j in seq_len(d_latent)) {
    draw_j <- samples[, , j, drop = FALSE]
    dim(draw_j) <- c(dim(samples)[1], length(miss_idx))
    post_mean <- colMeans(draw_j)
    lower <- apply(draw_j, 2, quantile, probs = 0.025, names = FALSE)
    upper <- apply(draw_j, 2, quantile, probs = 0.975, names = FALSE)
    truth <- train$U[miss_idx, j]
    error <- post_mean - truth

    imputation_rows[[paste(n_calib, j)]] <- data.frame(
      n_calib = n_calib,
      coordinate = paste0("u", j),
      id = miss_idx,
      truth = truth,
      posterior_mean = post_mean,
      lower = lower,
      upper = upper
    )
    imputation_metrics[[paste(n_calib, j)]] <- data.frame(
      n_calib = n_calib,
      coordinate = paste0("u", j),
      n_missing = length(miss_idx),
      Bias = mean(error),
      RMSE = sqrt(mean(error^2)),
      MAE = mean(abs(error)),
      Coverage95 = mean(truth >= lower & truth <= upper),
      Width95 = mean(upper - lower)
    )
  }
}

imputation_df <- bind_rows(imputation_rows) |>
  mutate(
    calibration = factor(
      paste0("|O| = ", n_calib),
      levels = paste0("|O| = ", rep_fit_calibs[rep_fit_calibs > 0L])
    ),
    coordinate = factor(coordinate, levels = c("u1", "u2"))
  )
imputation_metrics_df <- bind_rows(imputation_metrics)

if (nrow(imputation_df) > 0L) {
  p_imputation <- ggplot(
    imputation_df,
    aes(x = truth, y = posterior_mean, ymin = lower, ymax = upper)
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray35") +
    geom_linerange(alpha = 0.18, linewidth = 0.3) +
    geom_point(size = 0.9, alpha = 0.65, color = "firebrick") +
    facet_grid(coordinate ~ calibration, scales = "free") +
    labs(
      x = "True latent coordinate",
      y = "Posterior mean (95% interval)",
      title = "Study II latent-state imputation on the calibrated scale"
    ) +
    theme_bw(base_size = 10)

  save_plot(
    file.path(FIG_DIR, "fig2_study2_latent_imputation_by_calibration"),
    p_imputation,
    width = 12,
    height = 6.2
  )
}

write.csv(
  imputation_metrics_df,
  file.path(TAB_DIR, "study2_imputation_metrics_by_coordinate.csv"),
  row.names = FALSE
)
if (nrow(imputation_metrics_df) > 0L) {
  tex <- knitr::kable(
    imputation_metrics_df,
    format = "latex",
    booktabs = TRUE,
    digits = 3,
    col.names = c(
      "$|\\mathcal O|$", "Coordinate", "$n_{\\rm mis}$", "Bias",
      "RMSE", "MAE", "Coverage", "Width"
    ),
    escape = FALSE
  )
  writeLines(tex, file.path(TAB_DIR, "study2_imputation_metrics_by_coordinate.tex"))
}

############################################################
## Figure 3: honest predictive-density comparison
############################################################

pattern_info <- classify_study2_pattern_frequency(train$C, test$C)
candidate_strata <- c("common", "rare", "unobserved")
selected_ids <- unlist(lapply(candidate_strata, function(stratum) {
  ids <- which(as.character(pattern_info$pattern_stratum) == stratum)
  if (length(ids) == 0L) return(integer(0))
  severity <- rowSums(test$C[ids, , drop = FALSE])
  ids[which.min(abs(severity - stats::median(severity)))]
}))
if (length(selected_ids) == 0L) selected_ids <- 1L

oracle_pool <- make_oracle_pool_2d(
  true_params = rep_dat$true_params,
  n_pool = n_oracle_pool,
  seed = 20260712L
)
oracle_draws <- sample_oracle_test_y_2d(
  X_test = test$X[selected_ids, , drop = FALSE],
  C_test = test$C[selected_ids, , drop = FALSE],
  true_params = rep_dat$true_params,
  sigma_eps = rep_dat$sigma_eps,
  n_draw = n_pred_draw,
  oracle_pool = oracle_pool,
  seed = 20260713L
)

draw_ids <- even_draw_ids(dim(main_fit$mcmc$samples_U)[1], n_pred_draw)
eiv_draws <- sample_eiv_test_y_ordprobit_fb(
  X_test_raw = test$X[selected_ids, , drop = FALSE],
  C_test = test$C[selected_ids, , drop = FALSE],
  fit_obj = main_fit,
  draw_ids = draw_ids,
  n_per_draw = 1L,
  latent_sampler = predictive_latent_sampler,
  n_new_latent_gibbs = diagnostic_n_new_latent_gibbs,
  rejection_max_batches = rejection_max_batches
)

competitor_result <- run_study2_published_competitors(
  X_train = train$X,
  y_train = train$y,
  C_train = train$C,
  X_test = test$X[selected_ids, , drop = FALSE],
  C_test = test$C[selected_ids, , drop = FALSE],
  n_draw = nrow(eiv_draws),
  seed = 20260720L,
  m_vec = m_vec,
  methods = STUDY2_PUBLISHED_COMPETITORS,
  strict = STUDY2_STRICT_COMPETITORS,
  controls = list(
    `UC-GP` = list(n_starts = if (STUDY2_CONFIG == "quick") 2L else 8L),
    LVGP = list(
      n_starts = if (STUDY2_CONFIG == "quick") 2L else 8L,
      max_iter_ini = if (STUDY2_CONFIG == "quick") 30L else 100L,
      max_iter_lat = if (STUDY2_CONFIG == "quick") 8L else 20L
    ),
    EzGP = list(maxeval = if (STUDY2_CONFIG == "quick") 30L else 100L)
  )
)
write.csv(
  competitor_result$status,
  file.path(TAB_DIR, "study2_representative_competitor_status.csv"),
  row.names = FALSE
)

predictive_draws <- c(
  list(Oracle = oracle_draws, `EIV-GP` = eiv_draws),
  competitor_result$draws
)
predictive_long <- bind_rows(lapply(names(predictive_draws), function(method) {
  draws <- predictive_draws[[method]]
  bind_rows(lapply(seq_along(selected_ids), function(j) {
    data.frame(
      method = method,
      test_id = selected_ids[j],
      draw = draws[, j]
    )
  }))
})) |>
  mutate(
    stratum = as.character(pattern_info$pattern_stratum[test_id]),
    panel = paste0(
      "Test input ", test_id, " (", stratum, "; C = ",
      vapply(test_id, function(i) paste(test$C[i, ], collapse = ","), character(1)),
      ")"
    )
  )

observed_test <- data.frame(
  test_id = selected_ids,
  y = test$y[selected_ids]
) |>
  mutate(
    stratum = as.character(pattern_info$pattern_stratum[test_id]),
    panel = paste0(
      "Test input ", test_id, " (", stratum, "; C = ",
      vapply(test_id, function(i) paste(test$C[i, ], collapse = ","), character(1)),
      ")"
    )
  )

method_cols <- c(
  "Oracle" = "black", "EIV-GP" = "firebrick",
  "UC-GP" = "steelblue4", "LVGP" = "darkorange3",
  "EzGP" = "purple4"
)

p_density <- ggplot(predictive_long, aes(draw, color = method, fill = method)) +
  geom_density(alpha = 0.08, linewidth = 0.8) +
  geom_vline(
    data = observed_test,
    aes(xintercept = y),
    inherit.aes = FALSE,
    linetype = "dotted",
    color = "gray30"
  ) +
  facet_wrap(~panel, scales = "free", ncol = 1) +
  scale_color_manual(values = method_cols, drop = FALSE, name = NULL) +
  scale_fill_manual(values = method_cols, drop = FALSE, name = NULL) +
  labs(
    x = expression(y^"*"), y = "Predictive density",
    title = paste0(
      "Study II predictive laws at |O| = ", STUDY2_MAIN_CALIB,
      " (dotted line: realized response)"
    )
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

save_plot(
  file.path(FIG_DIR, "fig3_study2_predictive_densities_selected"),
  p_density,
  width = 9.5,
  height = 3.2 + 2.2 * length(selected_ids)
)

############################################################
## Observed-input mean m(x,c): common, rare, and new patterns
############################################################

mean_x1 <- seq(-1, 1, length.out = if (STUDY2_CONFIG == "quick") 31L else 61L)
mean_pattern_rows <- selected_ids
mean_X <- do.call(rbind, lapply(mean_pattern_rows, function(i) {
  cbind(x1 = mean_x1, x2 = 0)
}))
colnames(mean_X) <- colnames(train$X)
mean_C <- do.call(rbind, lapply(mean_pattern_rows, function(i) {
  test$C[rep(i, length(mean_x1)), , drop = FALSE]
}))
mean_truth <- oracle_m0_2d(
  X = mean_X,
  C = mean_C,
  true_params = rep_dat$true_params,
  oracle_pool = oracle_pool,
  n_latent = if (STUDY2_CONFIG == "quick") 1000L else 2000L,
  seed = 20260730L
)
write.csv(
  transform(
    attr(mean_truth, "rejection_telemetry"),
    figure = "fig_study2_mean_function_recovery"
  ),
  file.path(TAB_DIR, "study2_representative_mean_truth_rejection.csv"),
  row.names = FALSE
)
write.csv(
  transform(
    attr(mean_truth, "truth_diagnostics"),
    figure = "fig_study2_mean_function_recovery"
  ),
  file.path(TAB_DIR, "study2_representative_mean_truth_diagnostics.csv"),
  row.names = FALSE
)
mean_calibs <- unique(c(
  min(rep_fit_calibs), STUDY2_MAIN_CALIB, max(rep_fit_calibs)
))
mean_draw_cap <- if (STUDY2_CONFIG == "quick") 50L else 120L
mean_latent_draws <- if (STUDY2_CONFIG == "quick") 64L else 256L

mean_frames <- lapply(mean_calibs, function(n_calib) {
  fit <- rep_fits[[as.character(n_calib)]]
  ids <- even_draw_ids(dim(fit$mcmc$samples_U)[1L], mean_draw_cap)
  draws <- sample_eiv_m_given_xc_fb(
    X_test_raw = mean_X,
    C_test = mean_C,
    fit_obj = fit,
    draw_ids = ids,
    n_latent = mean_latent_draws,
    include_process_uncertainty = TRUE,
    joint = FALSE,
    return_components = TRUE,
    latent_sampler = predictive_latent_sampler,
    n_new_latent_gibbs = diagnostic_n_new_latent_gibbs,
    rejection_max_batches = rejection_max_batches,
    seed = 20260800L + n_calib
  )
  qs <- apply(draws, 2L, quantile, probs = c(0.025, 0.975), names = FALSE)
  data.frame(
    x1 = mean_X[, 1L],
    estimate = colMeans(attr(draws, "conditional_means")),
    lo = qs[1L, ],
    hi = qs[2L, ],
    truth = mean_truth,
    calibration = paste0("|O| = ", n_calib),
    pattern = rep(
      paste0(
        as.character(pattern_info$pattern_stratum[mean_pattern_rows]),
        ": C=",
        vapply(
          mean_pattern_rows,
          function(i) paste(test$C[i, ], collapse = ","),
          character(1)
        )
      ),
      each = length(mean_x1)
    )
  )
})
mean_df <- bind_rows(mean_frames)
mean_df$calibration <- factor(
  mean_df$calibration, levels = paste0("|O| = ", mean_calibs)
)

p_mean <- ggplot(mean_df, aes(x1)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "skyblue", alpha = 0.4) +
  geom_line(aes(y = estimate), color = "blue", linewidth = 0.75) +
  geom_line(aes(y = truth), color = "orange", linewidth = 0.9) +
  facet_grid(pattern ~ calibration, scales = "free_y") +
  labs(
    x = expression(x[1]~(x[2] == 0)),
    y = expression(m(x, c)),
    title = "Study II observed-input mean recovery",
    subtitle = "Orange: truth; blue: posterior mean and pointwise 95% band"
  ) +
  theme_bw(base_size = 9)

save_plot(
  file.path(FIG_DIR, "fig_study2_mean_function_recovery"),
  p_mean,
  width = 10.5,
  height = 3.0 + 2.3 * length(mean_pattern_rows)
)

############################################################
## Figure 4: latent response-surface recovery
############################################################

surface_frame <- function(fit, n_calib, grid, max_draw = 120L) {
  U_star <- as.matrix(grid[, c("u1", "u2")])
  X_star_raw <- matrix(0, nrow(grid), fit$data$p)
  X_star <- sweep(
    sweep(X_star_raw, 2, fit$data$X_center, "-"),
    2, fit$data$X_scale, "/"
  )
  ids <- even_draw_ids(dim(fit$mcmc$samples_U)[1], max_draw)
  draws <- matrix(NA_real_, nrow = length(ids), ncol = nrow(grid))
  kernel_spec <- kernel_spec_from_fit(fit)

  for (ii in seq_along(ids)) {
    s <- ids[ii]
    pred <- gp_predict_draw_general(
      X_train = fit$data$X,
      U_train = fit$mcmc$samples_U[s, , ],
      y_train = fit$data$y,
      X_star = X_star,
      U_star = U_star,
      logtheta = fit$mcmc$samples_logtheta[s, ],
      sigma2_eps = fit$mcmc$samples_sigma2[s],
      noisy = FALSE,
      kernel = kernel_spec$name,
      matern_nu = kernel_spec$matern_nu
    )
    draws[ii, ] <- fit$data$y_center + fit$data$y_scale * pred$mean
  }

  transform(
    grid,
    value = colMeans(draws),
    surface = paste0("EIV-GP, |O| = ", n_calib)
  )
}

grid_n <- if (STUDY2_CONFIG == "quick") 25L else 41L
surface_grid <- expand.grid(
  u1 = seq(-2.5, 2.5, length.out = grid_n),
  u2 = seq(-2.5, 2.5, length.out = grid_n)
)
true_surface <- transform(
  surface_grid,
  value = f0_2d(
    matrix(0, nrow(surface_grid), 2L),
    as.matrix(surface_grid),
    scenario = STUDY2_REP_SCENARIO
  ),
  surface = "Truth"
)

selected_surface_calibs <- intersect(c(10L, 25L, 80L), rep_fit_calibs)
if (length(selected_surface_calibs) == 0L) {
  selected_surface_calibs <- rep_fit_calibs[rep_fit_calibs > 0L][1]
}
surface_df <- bind_rows(
  list(true_surface),
  lapply(selected_surface_calibs, function(n_calib) {
    surface_frame(
      fit = rep_fits[[as.character(n_calib)]],
      n_calib = n_calib,
      grid = surface_grid,
      max_draw = if (STUDY2_CONFIG == "quick") 40L else 120L
    )
  })
)
surface_df$surface <- factor(
  surface_df$surface,
  levels = c("Truth", paste0("EIV-GP, |O| = ", selected_surface_calibs))
)

p_surface <- ggplot(surface_df, aes(u1, u2, fill = value)) +
  geom_raster() +
  geom_contour(aes(z = value), color = "white", alpha = 0.55, bins = 8) +
  facet_wrap(~surface, nrow = 1) +
  scale_fill_viridis_c(name = expression(f(x == 0, u))) +
  coord_equal() +
  labs(
    x = expression(u[1]), y = expression(u[2]),
    title = "Study II latent response-surface recovery at x = (0, 0)"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

save_plot(
  file.path(FIG_DIR, "fig4_study2_surface_recovery_selected"),
  p_surface,
  width = 3.2 + 3.0 * nlevels(surface_df$surface),
  height = 4.1
)

############################################################
## Reproducibility records
############################################################

diagnostics <- bind_rows(lapply(rep_fit_calibs, function(n_calib) {
  out <- as.data.frame(rep_fits[[as.character(n_calib)]]$diagnostics$summary)
  out$n_calib <- n_calib
  out
}))
write.csv(
  diagnostics,
  file.path(TAB_DIR, "study2_representative_mcmc_diagnostics.csv"),
  row.names = FALSE
)

manifest <- data.frame(
  design_tag = STUDY2_DESIGN_TAG,
  scenario = STUDY2_REP_SCENARIO,
  n_train = n_train,
  n_test = n_test,
  calibration_sizes = paste(rep_fit_calibs, collapse = ";"),
  main_calibration_size = STUDY2_MAIN_CALIB,
  n_iter = rep_n_iter,
  burn = rep_burn,
  thin = rep_thin,
  n_chains = rep_n_chains,
  sampler_strategy = "interwoven",
  prospective_latent_sampler = switch(
    predictive_latent_sampler,
    minimax_tilting = paste(
      "minimax-tilted accept-reject",
      "(exact on successful completion)"
    ),
    rejection = "exact prior rejection",
    gibbs = "finite Gibbs diagnostic"
  ),
  diagnostic_gibbs_sweeps = diagnostic_n_new_latent_gibbs,
  rejection_max_batches = rejection_max_batches,
  primary_target = "m(x,c)=E[f(x,U)|C=c]",
  mean_calibration_sizes = paste(mean_calibs, collapse = ";"),
  mean_posterior_draws = mean_draw_cap,
  mean_latent_integration_draws = mean_latent_draws,
  cache_tag = CACHE_TAG
)
write.csv(
  manifest,
  file.path(RES_DIR, paste0(CACHE_TAG, "_manifest.csv")),
  row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(RES_DIR, paste0(CACHE_TAG, "_sessionInfo.txt"))
)

cat("\nRepresentative Study II outputs completed under:\n")
cat(STUDY2_DESIGN_TAG, "\n")
