############################################################
## 01_study1_representative_figures.R
##
## Representative-data figures and diagnostics for revised Study I.
##
## Revised research version:
##   - keeps the original publication-style figures/tables;
##   - by default also creates many exploratory plots under
##       figures/study1_research_extra/
##   - fits EIV-GP for all calibration sizes in calib_grid by default
##     when STUDY1_RESEARCH_PLOTS = TRUE.
############################################################

if (!exists("fit_eivgp_1d")) {
  source("00_study1_functions.R")
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
})

if (!exists("STUDY1_QUICK")) STUDY1_QUICK <- FALSE
if (!exists("STUDY1_USE_CACHE")) STUDY1_USE_CACHE <- TRUE
if (!exists("STUDY1_OUT_PREFIX")) STUDY1_OUT_PREFIX <- ".."

## New switches.
## Set STUDY1_RESEARCH_PLOTS <- FALSE before sourcing this file if you only
## want the original small publication figure set.
if (!exists("STUDY1_RESEARCH_PLOTS")) STUDY1_RESEARCH_PLOTS <- TRUE
if (!exists("STUDY1_SAVE_PDF")) STUDY1_SAVE_PDF <- TRUE
if (!exists("STUDY1_SAVE_PNG")) STUDY1_SAVE_PNG <- TRUE
if (!exists("STUDY1_MAKE_TEST_PREDICTIONS")) {
  STUDY1_MAKE_TEST_PREDICTIONS <- STUDY1_RESEARCH_PLOTS
}
if (!exists("STUDY1_MAKE_SURFACE_HEATMAPS")) {
  STUDY1_MAKE_SURFACE_HEATMAPS <- STUDY1_RESEARCH_PLOTS
}

FIG_DIR <- file.path(STUDY1_OUT_PREFIX, "figures")
TAB_DIR <- file.path(STUDY1_OUT_PREFIX, "tables")
RES_DIR <- file.path(STUDY1_OUT_PREFIX, "results", "study1_1d")

## Extra research-output subdirectories.
EXP_FIG_DIR <- file.path(FIG_DIR, "study1_research_extra")
DATA_FIG_DIR <- file.path(EXP_FIG_DIR, "01_data")
IMP_FIG_DIR <- file.path(EXP_FIG_DIR, "02_latent_imputation")
FUN_FIG_DIR <- file.path(EXP_FIG_DIR, "03_function_surfaces")
PRED_FIG_DIR <- file.path(EXP_FIG_DIR, "04_predictive_distributions")
TEST_PRED_FIG_DIR <- file.path(EXP_FIG_DIR, "05_test_prediction")
MCMC_FIG_DIR <- file.path(EXP_FIG_DIR, "06_mcmc")

for (dd in c(
  FIG_DIR, TAB_DIR, RES_DIR, EXP_FIG_DIR, DATA_FIG_DIR, IMP_FIG_DIR,
  FUN_FIG_DIR, PRED_FIG_DIR, TEST_PRED_FIG_DIR, MCMC_FIG_DIR
)) {
  dir.create(dd, showWarnings = FALSE, recursive = TRUE)
}

set.seed(20260705)

n_train <- 100L
n_test <- if (STUDY1_QUICK) 150L else 500L
m <- 6L
calib_grid <- c(0L, 10L, 50L)

rep_n_iter <- if (STUDY1_QUICK) 900L else 6000L
rep_burn <- if (STUDY1_QUICK) 300L else 2000L
rep_n_chains <- if (STUDY1_QUICK) 2L else 12L
rep_preset <- if (STUDY1_QUICK) "fast" else "thorough"

n_pred_draw <- if (STUDY1_QUICK) 150L else 600L
n_density_draw <- if (STUDY1_QUICK) 300L else 1200L

## Extra research draws. These are intentionally smaller than the main
## publication predictive draw counts for expensive all-grid/all-calibration
## visualizations.
n_research_pred_draw <- if (STUDY1_QUICK) 80L else 250L
n_surface_draw <- if (STUDY1_QUICK) 50L else 200L
n_test_pred_draw <- if (STUDY1_QUICK) 150L else 600L
n_density_draw_research <- if (STUDY1_QUICK) 250L else 800L

parallel_chains <- (
  .Platform$OS.type != "windows" &&
    parallel::detectCores(logical = TRUE) > 2L
)

class_cols <- setNames(
  c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"),
  as.character(seq_len(m))
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
  "GP-Gaussian" = "darkgreen",
  "Complete-case GP" = "skyblue4"
)

############################################################
## Helper functions
############################################################

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

save_plot <- function(path_no_ext,
                      plot,
                      width,
                      height,
                      dpi = 320) {
  dir.create(dirname(path_no_ext), showWarnings = FALSE, recursive = TRUE)
  
  if (isTRUE(STUDY1_SAVE_PDF)) {
    ggsave(
      filename = paste0(path_no_ext, ".pdf"),
      plot = plot,
      width = width,
      height = height,
      units = "in",
      limitsize = FALSE
    )
  }
  
  if (isTRUE(STUDY1_SAVE_PNG)) {
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

safe_save_plot <- function(path_no_ext,
                           plot,
                           width,
                           height,
                           dpi = 320) {
  tryCatch(
    save_plot(path_no_ext, plot, width, height, dpi),
    error = function(e) {
      warning(
        "Could not save plot ", basename(path_no_ext), ": ",
        conditionMessage(e),
        call. = FALSE
      )
      invisible(NULL)
    }
  )
}

save_csv <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(x, path, row.names = FALSE)
  invisible(x)
}

select_ids_by_quantile <- function(ids, values, n = 12L) {
  if (length(ids) == 0L) return(integer(0))
  keep <- is.finite(values)
  ids <- ids[keep]
  values <- values[keep]
  if (length(ids) == 0L) return(integer(0))
  
  ord_ids <- ids[order(values)]
  ord_ids[unique(round(seq(1, length(ord_ids), length.out = min(n, length(ord_ids)))))]
}

raw_x_from_fit <- function(fit, x = fit$data$x) {
  center <- fit$data$x_center %||% 0
  scale <- fit$data$x_scale %||% 1
  as.numeric(center + scale * x)
}

raw_y_from_fit <- function(fit, y = fit$data$y) {
  center <- fit$data$y_center %||% 0
  scale <- fit$data$y_scale %||% 1
  as.numeric(center + scale * y)
}

make_long_draws <- function(draw_mat, method, grid) {
  draw_mat <- as.matrix(draw_mat)
  stopifnot(ncol(draw_mat) == nrow(grid))
  
  data.frame(
    x_star = rep(grid$x_star, each = nrow(draw_mat)),
    c_star = rep(grid$c_star, each = nrow(draw_mat)),
    y = as.vector(draw_mat),
    method = method
  )
}

summarise_draw_matrix <- function(draw_mat, method, grid) {
  draw_mat <- as.matrix(draw_mat)
  
  ## Expected orientation: rows = posterior/predictive draws,
  ##                       columns = prediction locations.
  ## If the matrix appears transposed, fix it automatically.
  if (ncol(draw_mat) != nrow(grid) && nrow(draw_mat) == nrow(grid)) {
    draw_mat <- t(draw_mat)
  }
  
  if (ncol(draw_mat) != nrow(grid)) {
    stop(
      "draw_mat has incompatible dimensions: nrow(draw_mat) = ",
      nrow(draw_mat), ", ncol(draw_mat) = ", ncol(draw_mat),
      ", but nrow(grid) = ", nrow(grid), ". ",
      "Expected columns of draw_mat to correspond to rows of grid."
    )
  }
  
  q_probs <- c(
    0.025,
    0.05,
    0.10,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95,
    0.975
  )
  
  q_names <- c(
    "lo95",
    "lo90",
    "lo80",
    "lo50",
    "med",
    "hi50",
    "hi80",
    "hi90",
    "hi95"
  )
  
  qfun <- function(z) {
    z <- z[is.finite(z)]
    if (length(z) == 0L) {
      return(rep(NA_real_, length(q_probs)))
    }
    as.numeric(stats::quantile(
      z,
      probs = q_probs,
      na.rm = TRUE,
      names = FALSE
    ))
  }
  
  qs <- vapply(
    seq_len(ncol(draw_mat)),
    function(j) qfun(draw_mat[, j]),
    numeric(length(q_probs))
  )
  
  rownames(qs) <- q_names
  
  out <- data.frame(
    grid,
    method = method,
    mean = colMeans(draw_mat, na.rm = TRUE),
    sd = apply(draw_mat, 2, stats::sd, na.rm = TRUE),
    lo95 = qs["lo95", ],
    lo90 = qs["lo90", ],
    lo80 = qs["lo80", ],
    lo50 = qs["lo50", ],
    med = qs["med", ],
    hi50 = qs["hi50", ],
    hi80 = qs["hi80", ],
    hi90 = qs["hi90", ],
    hi95 = qs["hi95", ]
  )
  
  out$width50 <- out$hi50 - out$lo50
  out$width80 <- out$hi80 - out$lo80
  out$width90 <- out$hi90 - out$lo90
  out$width95 <- out$hi95 - out$lo95
  
  out
}

############################################################
## Calibration sizes to fit
############################################################

if (exists("STUDY1_REP_FIT_CALIBS")) {
  rep_fit_calibs <- sort(unique(as.integer(STUDY1_REP_FIT_CALIBS)))
} else {
  rep_fit_calibs <- if (STUDY1_RESEARCH_PLOTS) calib_grid else c(0L, 10L, 50L)
}

if (length(rep_fit_calibs) == 0L) {
  rep_fit_calibs <- 10L
}

if (exists("STUDY1_MAIN_CALIB")) {
  main_calib <- as.integer(STUDY1_MAIN_CALIB)
} else {
  main_calib <- if (10L %in% rep_fit_calibs) {
    10L
  } else {
    rep_fit_calibs[ceiling(length(rep_fit_calibs) / 2)]
  }
}

rep_fit_calibs <- sort(unique(c(rep_fit_calibs, main_calib)))

## Ensure calibration sets exist for all requested fits.
calib_grid_for_sets <- sort(unique(c(calib_grid, rep_fit_calibs)))

calib_label <- function(x, levels = rep_fit_calibs) {
  factor(paste0("|O| = ", x), levels = paste0("|O| = ", levels))
}

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
  calib_grid = calib_grid_for_sets,
  seed = 20260706
)

############################################################
## Fit EIV-GP for representative calibration sizes
############################################################

if (STUDY1_USE_CACHE && file.exists(rep_fits_file)) {
  rep_fits <- readRDS(rep_fits_file)
  if (!is.list(rep_fits)) rep_fits <- list()
} else {
  rep_fits <- list()
}

missing_fit_calibs <- setdiff(as.character(rep_fit_calibs), names(rep_fits))

if (length(missing_fit_calibs) > 0L) {
  for (kk_chr in missing_fit_calibs) {
    kk <- as.integer(kk_chr)
    
    cat("Fitting representative EIV-GP with |O| =", kk, "\n")
    
    rep_fits[[kk_chr]] <- fit_eivgp_1d(
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

rep_fits <- rep_fits[as.character(rep_fit_calibs)]
main_fit <- rep_fits[[as.character(main_calib)]]

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
  id = seq_len(n_train),
  x = train$x,
  u = train$u,
  c = factor(train$c, levels = seq_len(m)),
  y = train$y
)

df_train$f_true <- f0_1d(df_train$x, df_train$u, scenario = rep_dat$scenario)
df_train$resid_true <- df_train$y - df_train$f_true

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

save_plot(
  file.path(FIG_DIR, "fig1_study1_data_active_x"),
  p_data_active,
  width = 10,
  height = 4.5
)

############################################################
## Extra data plots
############################################################

if (STUDY1_RESEARCH_PLOTS) {
  p_xu <- ggplot(df_train, aes(x = x, y = u, color = c)) +
    geom_point(size = 2, alpha = 0.85) +
    geom_hline(yintercept = tau_true, color = "grey30", linewidth = 0.35) +
    scale_color_manual(values = class_cols, name = "Class") +
    labs(
      x = "Observed x",
      y = "Latent u",
      title = "Latent coordinate versus observed quantitative input"
    ) +
    theme(legend.position = "bottom")
  
  p_y_truth <- ggplot(df_train, aes(x = f_true, y = y, color = c)) +
    geom_point(size = 2, alpha = 0.85) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(values = class_cols, name = "Class") +
    labs(
      x = "Noise-free truth f0(x, u)",
      y = "Observed y",
      title = "Observed response versus noise-free truth"
    ) +
    theme(legend.position = "bottom")
  
  p_noise <- ggplot(df_train, aes(x = c, y = resid_true, fill = c)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
    geom_boxplot(alpha = 0.75, outlier.alpha = 0.6) +
    scale_fill_manual(values = class_cols, name = "Class") +
    labs(
      x = "Class",
      y = "y - f0(x, u)",
      title = "Realized noise by ordinal class"
    ) +
    theme(legend.position = "none")
  
  p_y_box <- ggplot(df_train, aes(x = c, y = y, fill = c)) +
    geom_violin(alpha = 0.45, color = NA) +
    geom_boxplot(width = 0.18, alpha = 0.8, outlier.alpha = 0.5) +
    scale_fill_manual(values = class_cols, name = "Class") +
    labs(
      x = "Class",
      y = "y",
      title = "Response distribution by ordinal class"
    ) +
    theme(legend.position = "none")
  
  p_data_scatter_more <- (p_xu + p_y_truth) / (p_noise + p_y_box)
  
  safe_save_plot(
    file.path(DATA_FIG_DIR, "explore_data_scatter_noise_panels"),
    p_data_scatter_more,
    width = 12,
    height = 9
  )
  
  df_train_long <- df_train |>
    tidyr::pivot_longer(
      cols = c(x, u, y, f_true, resid_true),
      names_to = "variable",
      values_to = "value"
    )
  
  p_marginals <- ggplot(df_train_long, aes(x = value, fill = c, color = c)) +
    geom_density(alpha = 0.16, linewidth = 0.7) +
    facet_wrap(~variable, scales = "free", ncol = 3) +
    scale_fill_manual(values = class_cols, name = "Class") +
    scale_color_manual(values = class_cols, name = "Class") +
    labs(
      x = NULL,
      y = "Density",
      title = "Marginal distributions of representative training data"
    ) +
    theme(legend.position = "bottom")
  
  safe_save_plot(
    file.path(DATA_FIG_DIR, "explore_data_marginal_densities_by_class"),
    p_marginals,
    width = 12,
    height = 7
  )
  
  df_calib_all <- dplyr::bind_rows(
    lapply(rep_fit_calibs, function(kk) {
      tmp <- df_train
      tmp$n_calib <- kk
      tmp$n_calib_label <- calib_label(kk)
      tmp$calibrated <- tmp$id %in% calib_sets[[as.character(kk)]]
      tmp
    })
  )
  
  p_calib_locations <- ggplot() +
    geom_point(
      data = df_calib_all[!df_calib_all$calibrated, ],
      aes(x = x, y = u),
      color = "grey82",
      size = 1.2,
      alpha = 0.75
    ) +
    geom_point(
      data = df_calib_all[df_calib_all$calibrated, ],
      aes(x = x, y = u, color = c),
      size = 2.2,
      alpha = 0.95
    ) +
    geom_hline(yintercept = tau_true, color = "grey35", linewidth = 0.25) +
    facet_wrap(~n_calib_label, nrow = 1) +
    scale_color_manual(values = class_cols, name = "Class") +
    labs(
      x = "Observed x",
      y = "Latent u",
      title = "Calibration-set locations in the x-u plane",
      subtitle = "Grey points are uncalibrated; colored points are calibrated"
    ) +
    theme(legend.position = "bottom")
  
  safe_save_plot(
    file.path(DATA_FIG_DIR, "explore_calibration_locations_x_u"),
    p_calib_locations,
    width = 15,
    height = 4.2
  )
  
  df_calib_counts <- df_calib_all |>
    dplyr::group_by(n_calib, n_calib_label, c) |>
    dplyr::summarise(
      training_count = dplyr::n(),
      calibrated_count = sum(calibrated),
      .groups = "drop"
    )
  
  p_calib_counts <- ggplot(df_calib_counts, aes(x = c, y = calibrated_count, fill = c)) +
    geom_col(alpha = 0.88) +
    facet_wrap(~n_calib_label, nrow = 1) +
    scale_fill_manual(values = class_cols, name = "Class") +
    labs(
      x = "Class",
      y = "Number calibrated",
      title = "Class composition of nested calibration sets"
    ) +
    theme(legend.position = "none")
  
  safe_save_plot(
    file.path(DATA_FIG_DIR, "explore_calibration_class_counts"),
    p_calib_counts,
    width = 13,
    height = 3.8
  )
  
  if (!is.null(test) && all(c("x", "c") %in% names(test))) {
    df_test <- data.frame(
      x = test$x,
      u = if ("u" %in% names(test)) test$u else rep(NA_real_, length(test$x)),
      c = factor(test$c, levels = seq_len(m)),
      y = if ("y" %in% names(test)) test$y else rep(NA_real_, length(test$x)),
      set = "Test"
    )
    
    df_train_tt <- data.frame(
      x = df_train$x,
      u = df_train$u,
      c = df_train$c,
      y = df_train$y,
      set = "Train"
    )
    
    df_tt <- dplyr::bind_rows(df_train_tt, df_test)
    
    tt_cols <- c(
      "x",
      if (any(!is.na(df_tt$u))) "u",
      if (any(!is.na(df_tt$y))) "y"
    )
    
    df_tt_long <- df_tt |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(tt_cols),
        names_to = "variable",
        values_to = "value"
      )
    
    p_train_test_density <- ggplot(df_tt_long, aes(x = value, color = set, fill = set)) +
      geom_density(alpha = 0.15, linewidth = 0.8, na.rm = TRUE) +
      facet_wrap(~variable, scales = "free", nrow = 1) +
      labs(
        x = NULL,
        y = "Density",
        color = NULL,
        fill = NULL,
        title = "Train-test marginal distribution comparison"
      ) +
      theme(legend.position = "bottom")
    
    safe_save_plot(
      file.path(DATA_FIG_DIR, "explore_train_test_marginal_comparison"),
      p_train_test_density,
      width = 11,
      height = 3.8
    )
    
    df_tt_counts <- df_tt |>
      dplyr::group_by(set, c) |>
      dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
      dplyr::group_by(set) |>
      dplyr::mutate(prop = n / sum(n)) |>
      dplyr::ungroup()
    
    p_train_test_class <- ggplot(df_tt_counts, aes(x = c, y = prop, fill = set)) +
      geom_col(position = "dodge", alpha = 0.85) +
      labs(
        x = "Class",
        y = "Proportion",
        fill = NULL,
        title = "Train-test class proportion comparison"
      ) +
      theme(legend.position = "bottom")
    
    safe_save_plot(
      file.path(DATA_FIG_DIR, "explore_train_test_class_proportions"),
      p_train_test_class,
      width = 7.5,
      height = 4.2
    )
  }
  
  x_lim <- range(
    c(
      df_train$x,
      -2,
      2,
      if (!is.null(test) && "x" %in% names(test)) test$x else numeric(0)
    ),
    finite = TRUE
  )
  u_lim <- range(
    c(
      df_train$u,
      tau_true,
      -2.6,
      2.6,
      if (!is.null(test) && "u" %in% names(test)) test$u else numeric(0)
    ),
    finite = TRUE
  )
  
  truth_grid <- expand.grid(
    x_raw = seq(x_lim[1], x_lim[2], length.out = 160),
    u = seq(u_lim[1], u_lim[2], length.out = 160)
  )
  truth_grid$f <- f0_1d(truth_grid$x_raw, truth_grid$u, scenario = rep_dat$scenario)
  
  p_truth_surface <- ggplot(truth_grid, aes(x = x_raw, y = u, fill = f)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = f), color = "white", alpha = 0.45, linewidth = 0.25) +
    geom_hline(yintercept = tau_true, color = "black", linewidth = 0.35) +
    geom_point(
      data = df_train,
      aes(x = x, y = u),
      inherit.aes = FALSE,
      shape = 21,
      fill = "white",
      color = "black",
      size = 0.8,
      alpha = 0.55
    ) +
    scale_fill_gradient2(
      low = "navy",
      mid = "white",
      high = "firebrick",
      midpoint = median(truth_grid$f, na.rm = TRUE),
      name = "f0"
    ) +
    labs(
      x = "Observed x",
      y = "Latent u",
      title = "True latent response surface",
      subtitle = "Black horizontal lines are true ordinal cutpoints; white points are training inputs"
    )
  
  safe_save_plot(
    file.path(DATA_FIG_DIR, "explore_true_response_surface_with_training_points"),
    p_truth_surface,
    width = 8,
    height = 6.2
  )
}

############################################################
## Figure 2: latent imputation across calibration sizes
############################################################

extract_imputation_df <- function(fit, n_calib) {
  samples_u <- fit$mcmc$samples_u
  miss_idx <- fit$data$miss_idx
  
  if (length(miss_idx) == 0L) {
    return(data.frame())
  }
  
  post_mean <- colMeans(samples_u, na.rm = TRUE)[miss_idx]
  post_lo <- apply(
    samples_u[, miss_idx, drop = FALSE],
    2,
    stats::quantile,
    probs = 0.025,
    na.rm = TRUE
  )
  post_hi <- apply(
    samples_u[, miss_idx, drop = FALSE],
    2,
    stats::quantile,
    probs = 0.975,
    na.rm = TRUE
  )
  post_sd <- apply(
    samples_u[, miss_idx, drop = FALSE],
    2,
    stats::sd,
    na.rm = TRUE
  )
  
  true_u <- fit$data$u_true[miss_idx]
  x_raw <- raw_x_from_fit(fit)[miss_idx]
  y_raw <- raw_y_from_fit(fit)[miss_idx]
  
  data.frame(
    n_calib = n_calib,
    n_calib_label = calib_label(n_calib),
    id = miss_idx,
    x = x_raw,
    y = y_raw,
    true_u = true_u,
    post_mean = post_mean,
    post_lo = post_lo,
    post_hi = post_hi,
    post_sd = post_sd,
    width = post_hi - post_lo,
    error = post_mean - true_u,
    abs_error = abs(post_mean - true_u),
    c = factor(fit$data$c_ord[miss_idx], levels = seq_len(m)),
    covered = true_u >= post_lo & true_u <= post_hi
  )
}

df_imp_all <- dplyr::bind_rows(
  lapply(names(rep_fits), function(nm) {
    extract_imputation_df(rep_fits[[nm]], as.integer(nm))
  })
)

df_imp_all$n_calib_label <- calib_label(df_imp_all$n_calib)

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

save_plot(
  file.path(FIG_DIR, "fig2_study1_latent_imputation_by_calibration"),
  p_imp_calib,
  width = max(11, 2.8 * length(rep_fit_calibs)),
  height = 4.5
)

############################################################
## Extra latent-imputation plots
############################################################

if (STUDY1_RESEARCH_PLOTS && nrow(df_imp_all) > 0L) {
  df_imp_metrics_class <- df_imp_all |>
    dplyr::group_by(n_calib, n_calib_label, c) |>
    dplyr::summarise(
      n = dplyr::n(),
      bias = mean(error),
      rmse = sqrt(mean(error^2)),
      mae = mean(abs_error),
      coverage = mean(covered),
      mean_width = mean(width),
      median_width = median(width),
      mean_sd = mean(post_sd),
      .groups = "drop"
    )
  
  df_imp_metrics_overall <- df_imp_all |>
    dplyr::group_by(n_calib, n_calib_label) |>
    dplyr::summarise(
      n = dplyr::n(),
      bias = mean(error),
      rmse = sqrt(mean(error^2)),
      mae = mean(abs_error),
      coverage = mean(covered),
      mean_width = mean(width),
      median_width = median(width),
      mean_sd = mean(post_sd),
      .groups = "drop"
    )
  
  save_csv(df_imp_metrics_class, file.path(RES_DIR, "representative_imputation_metrics_by_class.csv"))
  save_csv(df_imp_metrics_overall, file.path(RES_DIR, "representative_imputation_metrics_overall.csv"))
  
  df_imp_metrics_long <- df_imp_metrics_overall |>
    tidyr::pivot_longer(
      cols = c(bias, rmse, mae, coverage, mean_width, mean_sd),
      names_to = "metric",
      values_to = "value"
    )
  
  p_imp_metrics <- ggplot(df_imp_metrics_long, aes(x = n_calib, y = value)) +
    geom_line(linewidth = 0.75, color = "grey25") +
    geom_point(size = 2.3, color = "firebrick") +
    facet_wrap(~metric, scales = "free_y", ncol = 3) +
    labs(
      x = "Calibration size |O|",
      y = NULL,
      title = "Overall latent-imputation metrics versus calibration size"
    )
  
  safe_save_plot(
    file.path(IMP_FIG_DIR, "explore_imputation_metrics_overall"),
    p_imp_metrics,
    width = 10,
    height = 6
  )
  
  p_imp_error_box <- ggplot(df_imp_all, aes(x = c, y = error, fill = c)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
    geom_boxplot(alpha = 0.82, outlier.alpha = 0.45) +
    facet_wrap(~n_calib_label, nrow = 1) +
    scale_fill_manual(values = class_cols, name = "Class") +
    labs(
      x = "Class",
      y = "Posterior mean error",
      title = "Latent-imputation error by class and calibration size"
    ) +
    theme(legend.position = "none")
  
  safe_save_plot(
    file.path(IMP_FIG_DIR, "explore_imputation_error_boxplots_by_class"),
    p_imp_error_box,
    width = max(12, 2.8 * length(rep_fit_calibs)),
    height = 4.2
  )
  
  p_imp_width_box <- ggplot(df_imp_all, aes(x = c, y = width, fill = c)) +
    geom_boxplot(alpha = 0.82, outlier.alpha = 0.45) +
    facet_wrap(~n_calib_label, nrow = 1) +
    scale_fill_manual(values = class_cols, name = "Class") +
    labs(
      x = "Class",
      y = "95% posterior interval width",
      title = "Latent-imputation uncertainty by class and calibration size"
    ) +
    theme(legend.position = "none")
  
  safe_save_plot(
    file.path(IMP_FIG_DIR, "explore_imputation_interval_width_boxplots_by_class"),
    p_imp_width_box,
    width = max(12, 2.8 * length(rep_fit_calibs)),
    height = 4.2
  )
  
  p_imp_coverage <- ggplot(df_imp_metrics_class, aes(x = c, y = coverage, fill = c)) +
    geom_col(alpha = 0.88) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
    facet_wrap(~n_calib_label, nrow = 1) +
    scale_fill_manual(values = class_cols, name = "Class") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "Class",
      y = "Empirical 95% coverage",
      title = "Latent-imputation interval coverage by class"
    ) +
    theme(legend.position = "none")
  
  safe_save_plot(
    file.path(IMP_FIG_DIR, "explore_imputation_coverage_by_class"),
    p_imp_coverage,
    width = max(12, 2.8 * length(rep_fit_calibs)),
    height = 4.2
  )
  
  p_imp_error_vs_true <- ggplot(df_imp_all, aes(x = true_u, y = error, color = c)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey35") +
    geom_point(alpha = 0.75, size = 1.5) +
    geom_smooth(se = FALSE, method = "loess", formula = y ~ x, linewidth = 0.7) +
    facet_wrap(~n_calib_label, nrow = 1) +
    scale_color_manual(values = class_cols, name = "Class") +
    labs(
      x = "True latent u",
      y = "Posterior mean error",
      title = "Latent-imputation error as a function of true u"
    ) +
    theme(legend.position = "bottom")
  
  safe_save_plot(
    file.path(IMP_FIG_DIR, "explore_imputation_error_vs_true_u"),
    p_imp_error_vs_true,
    width = max(12, 2.8 * length(rep_fit_calibs)),
    height = 4.4
  )
  
  p_imp_width_vs_true <- ggplot(df_imp_all, aes(x = true_u, y = width, color = c)) +
    geom_point(alpha = 0.75, size = 1.5) +
    geom_smooth(se = FALSE, method = "loess", formula = y ~ x, linewidth = 0.7) +
    facet_wrap(~n_calib_label, nrow = 1) +
    scale_color_manual(values = class_cols, name = "Class") +
    labs(
      x = "True latent u",
      y = "95% posterior interval width",
      title = "Latent-imputation uncertainty as a function of true u"
    ) +
    theme(legend.position = "bottom")
  
  safe_save_plot(
    file.path(IMP_FIG_DIR, "explore_imputation_width_vs_true_u"),
    p_imp_width_vs_true,
    width = max(12, 2.8 * length(rep_fit_calibs)),
    height = 4.4
  )
  
  df_imp_cat <- df_imp_all |>
    dplyr::group_by(n_calib, n_calib_label) |>
    dplyr::arrange(true_u, .by_group = TRUE) |>
    dplyr::mutate(order = dplyr::row_number()) |>
    dplyr::ungroup()
  
  p_imp_caterpillar <- ggplot(df_imp_cat, aes(y = order)) +
    geom_segment(
      aes(x = post_lo, xend = post_hi, yend = order, color = covered),
      alpha = 0.45,
      linewidth = 0.55
    ) +
    geom_point(aes(x = post_mean, color = covered), size = 0.9) +
    geom_point(aes(x = true_u), shape = 4, size = 0.7, color = "black") +
    facet_wrap(~n_calib_label, scales = "free_y") +
    scale_color_manual(values = c("TRUE" = "darkgreen", "FALSE" = "red3")) +
    labs(
      x = "Latent u",
      y = "Missing point ordered by true u",
      color = "Covered",
      title = "Caterpillar plot of latent-u posterior intervals"
    ) +
    theme(legend.position = "bottom")
  
  safe_save_plot(
    file.path(IMP_FIG_DIR, "explore_imputation_caterpillar_intervals"),
    p_imp_caterpillar,
    width = 14,
    height = 8.5
  )
  
  make_u_density_df <- function(fit, ids, n_calib = NULL, max_draw = 1000L) {
    if (length(ids) == 0L) return(data.frame())
    
    draws <- fit$mcmc$samples_u
    ids <- ids[ids >= 1L & ids <= ncol(draws)]
    if (length(ids) == 0L) return(data.frame())
    
    draw_ids <- seq_len(nrow(draws))
    if (length(draw_ids) > max_draw) {
      draw_ids <- sample(draw_ids, max_draw)
    }
    
    dplyr::bind_rows(
      lapply(ids, function(id) {
        data.frame(
          draw = draw_ids,
          id = id,
          id_label = paste0("id=", id, ", c=", fit$data$c_ord[id]),
          u_draw = draws[draw_ids, id],
          true_u = fit$data$u_true[id],
          c = factor(fit$data$c_ord[id], levels = seq_len(m)),
          n_calib = n_calib %||% NA_integer_
        )
      })
    )
  }
  
  selected_imp_ids <- select_ids_by_quantile(
    main_fit$data$miss_idx,
    main_fit$data$u_true[main_fit$data$miss_idx],
    n = 12L
  )
  
  df_u_density_main <- make_u_density_df(
    main_fit,
    selected_imp_ids,
    n_calib = main_calib,
    max_draw = min(1000L, nrow(main_fit$mcmc$samples_u))
  )
  
  if (nrow(df_u_density_main) > 0L) {
    p_u_density_main <- ggplot(df_u_density_main, aes(x = u_draw, fill = c, color = c)) +
      geom_density(alpha = 0.22, linewidth = 0.75) +
      geom_vline(
        data = unique(df_u_density_main[, c("id_label", "true_u")]),
        aes(xintercept = true_u),
        inherit.aes = FALSE,
        linetype = "dashed",
        color = "black"
      ) +
      facet_wrap(~id_label, scales = "free", ncol = 4) +
      scale_fill_manual(values = class_cols, name = "Class") +
      scale_color_manual(values = class_cols, name = "Class") +
      labs(
        x = "Posterior draw of latent u",
        y = "Density",
        title = paste0("Selected latent-u posterior densities for |O| = ", main_calib)
      ) +
      theme(legend.position = "bottom")
    
    safe_save_plot(
      file.path(IMP_FIG_DIR, "explore_selected_u_posterior_densities_main_fit"),
      p_u_density_main,
      width = 13,
      height = 8
    )
  }
  
  selected_imp_ids_small <- select_ids_by_quantile(
    main_fit$data$miss_idx,
    main_fit$data$u_true[main_fit$data$miss_idx],
    n = 6L
  )
  
  df_u_density_all <- dplyr::bind_rows(
    lapply(names(rep_fits), function(nm) {
      make_u_density_df(
        rep_fits[[nm]],
        selected_imp_ids_small,
        n_calib = as.integer(nm),
        max_draw = min(600L, nrow(rep_fits[[nm]]$mcmc$samples_u))
      )
    })
  )
  
  if (nrow(df_u_density_all) > 0L) {
    df_u_density_all$n_calib_label <- calib_label(df_u_density_all$n_calib)
    df_u_density_all <- df_u_density_all |>
      dplyr::group_by(id, n_calib) |>
      dplyr::filter(dplyr::n_distinct(round(u_draw, 10)) > 1L) |>
      dplyr::ungroup()
    
    if (nrow(df_u_density_all) > 0L) {
      p_u_density_all <- ggplot(df_u_density_all, aes(x = u_draw, color = n_calib_label)) +
        geom_density(linewidth = 0.8) +
        geom_vline(
          data = unique(df_u_density_all[, c("id_label", "true_u")]),
          aes(xintercept = true_u),
          inherit.aes = FALSE,
          linetype = "dashed",
          color = "black"
        ) +
        facet_wrap(~id_label, scales = "free", ncol = 3) +
        labs(
          x = "Posterior draw of latent u",
          y = "Density",
          color = "Calibration",
          title = "Selected latent-u posterior densities across calibration sizes"
        ) +
        theme(legend.position = "bottom")
      
      safe_save_plot(
        file.path(IMP_FIG_DIR, "explore_selected_u_posterior_densities_all_calibrations"),
        p_u_density_all,
        width = 12,
        height = 7.5
      )
    }
  }
}

############################################################
## Extra cutpoint/tau posterior plots, if available
############################################################

get_tau_matrix <- function(fit) {
  candidate_names <- c(
    "samples_tau",
    "samples_taus",
    "samples_cutpoints",
    "samples_alpha"
  )
  
  obj <- NULL
  for (nm in candidate_names) {
    if (!is.null(fit$mcmc[[nm]])) {
      obj <- fit$mcmc[[nm]]
      break
    }
  }
  
  if (is.null(obj)) return(NULL)
  
  if (is.null(dim(obj))) {
    if (length(tau_true) > 1L && length(obj) %% length(tau_true) == 0L) {
      mat <- matrix(as.numeric(obj), ncol = length(tau_true), byrow = TRUE)
    } else {
      mat <- matrix(as.numeric(obj), ncol = 1L)
    }
  } else {
    mat <- as.matrix(obj)
  }
  
  if (ncol(mat) != length(tau_true) && nrow(mat) == length(tau_true)) {
    mat <- t(mat)
  }
  
  mat
}

extract_tau_df <- function(fit, n_calib) {
  mat <- get_tau_matrix(fit)
  if (is.null(mat)) return(data.frame())
  
  dplyr::bind_rows(
    lapply(seq_len(ncol(mat)), function(jj) {
      data.frame(
        n_calib = n_calib,
        n_calib_label = calib_label(n_calib),
        draw = seq_len(nrow(mat)),
        tau_index = factor(jj, levels = seq_len(length(tau_true))),
        tau = mat[, jj]
      )
    })
  )
}

df_tau_all <- dplyr::bind_rows(
  lapply(names(rep_fits), function(nm) {
    extract_tau_df(rep_fits[[nm]], as.integer(nm))
  })
)

if (STUDY1_RESEARCH_PLOTS && nrow(df_tau_all) > 0L) {
  df_tau_truth <- data.frame(
    tau_index = factor(seq_along(tau_true), levels = seq_len(length(tau_true))),
    tau_true = tau_true
  )
  
  p_tau_density <- ggplot(df_tau_all, aes(x = tau, color = n_calib_label, fill = n_calib_label)) +
    geom_density(alpha = 0.10, linewidth = 0.8) +
    geom_vline(
      data = df_tau_truth,
      aes(xintercept = tau_true),
      inherit.aes = FALSE,
      linetype = "dashed",
      color = "black"
    ) +
    facet_wrap(~tau_index, scales = "free", ncol = min(3, length(tau_true))) +
    labs(
      x = "Cutpoint",
      y = "Density",
      color = "Calibration",
      fill = "Calibration",
      title = "Posterior cutpoint distributions across calibration sizes"
    ) +
    theme(legend.position = "bottom")
  
  safe_save_plot(
    file.path(IMP_FIG_DIR, "explore_tau_posterior_densities"),
    p_tau_density,
    width = 12,
    height = 6.5
  )
  
  df_tau_summary <- df_tau_all |>
    dplyr::group_by(n_calib, n_calib_label, tau_index) |>
    dplyr::summarise(
      mean = mean(tau),
      median = median(tau),
      lo = stats::quantile(tau, 0.025),
      hi = stats::quantile(tau, 0.975),
      .groups = "drop"
    )
  
  p_tau_interval <- ggplot(df_tau_summary, aes(x = n_calib, y = median, color = tau_index)) +
    geom_line(aes(group = tau_index), linewidth = 0.6) +
    geom_pointrange(aes(ymin = lo, ymax = hi), linewidth = 0.55) +
    geom_hline(
      data = df_tau_truth,
      aes(yintercept = tau_true, color = tau_index),
      linetype = "dashed",
      linewidth = 0.45
    ) +
    facet_wrap(~tau_index, scales = "free_y") +
    labs(
      x = "Calibration size |O|",
      y = "Posterior median and 95% interval",
      color = "Cutpoint",
      title = "Cutpoint posterior intervals versus calibration size"
    ) +
    theme(legend.position = "bottom")
  
  safe_save_plot(
    file.path(IMP_FIG_DIR, "explore_tau_posterior_intervals_by_calibration"),
    p_tau_interval,
    width = 11,
    height = 6
  )
}

############################################################
## Figure 3: latent response-surface slices
############################################################

make_eiv_function_slice_df <- function(fit,
                                       x_slices = c(-1, 0, 1),
                                       u_grid = seq(-2.4, 2.4, length.out = 160),
                                       draw_ids = NULL,
                                       scenario = "active",
                                       label = "EIV-GP",
                                       max_draw = n_pred_draw,
                                       include_process_uncertainty = TRUE) {
  if (is.null(draw_ids)) {
    draw_ids <- seq_len(nrow(fit$mcmc$samples_u))
  }
  
  if (!is.null(max_draw) && length(draw_ids) > max_draw) {
    draw_ids <- sample(draw_ids, max_draw)
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
    
    if (include_process_uncertainty) {
      f_std <- pred$mean + sqrt(pmax(pred$var, 0)) * rnorm(nrow(grid))
    } else {
      f_std <- pred$mean
    }
    
    f_samps[ii, ] <- fit$data$y_center + fit$data$y_scale * f_std
  }
  
  out <- grid
  out$mean <- colMeans(f_samps, na.rm = TRUE)
  out$lo <- apply(f_samps, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
  out$hi <- apply(f_samps, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
  out$width <- out$hi - out$lo
  out$truth <- f0_1d(out$x_raw, out$u, scenario = scenario)
  out$error <- out$mean - out$truth
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
  out$lo <- fit$data$y_center + fit$data$y_scale * (pred$mean - 1.96 * sqrt(pmax(pred$var, 0)))
  out$hi <- fit$data$y_center + fit$data$y_scale * (pred$mean + 1.96 * sqrt(pmax(pred$var, 0)))
  out$width <- out$hi - out$lo
  out$truth <- f0_1d(out$x_raw, out$u, scenario = scenario)
  out$error <- out$mean - out$truth
  out$method <- label
  out$x_slice <- factor(paste0("x = ", out$x_raw), levels = paste0("x = ", x_slices))
  
  out
}

df_fun_eiv <- make_eiv_function_slice_df(
  main_fit,
  scenario = rep_dat$scenario,
  max_draw = n_pred_draw
)
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

save_plot(
  file.path(FIG_DIR, "fig3_study1_function_slices"),
  p_fun_slices,
  width = 11,
  height = 6.5
)

############################################################
## Extra function-slice and surface plots
############################################################

if (STUDY1_RESEARCH_PLOTS) {
  x_slices_more <- seq(-1.5, 1.5, by = 0.5)
  u_grid_more <- seq(-2.6, 2.6, length.out = 180)
  calib_tag <- paste(rep_fit_calibs, collapse = "-")
  
  fun_slices_file <- file.path(
    RES_DIR,
    paste0(
      "representative_function_slices_all_calibs_",
      "calibs_", calib_tag,
      "_ndraw", n_research_pred_draw,
      ".rds"
    )
  )
  
  if (STUDY1_USE_CACHE && file.exists(fun_slices_file)) {
    df_fun_eiv_all <- readRDS(fun_slices_file)
  } else {
    df_fun_eiv_all <- dplyr::bind_rows(
      lapply(names(rep_fits), function(nm) {
        out <- make_eiv_function_slice_df(
          rep_fits[[nm]],
          x_slices = x_slices_more,
          u_grid = u_grid_more,
          scenario = rep_dat$scenario,
          label = "EIV-GP",
          max_draw = n_research_pred_draw,
          include_process_uncertainty = TRUE
        )
        out$n_calib <- as.integer(nm)
        out$n_calib_label <- calib_label(as.integer(nm))
        out
      })
    )
    
    saveRDS(df_fun_eiv_all, fun_slices_file)
  }
  
  df_fun_eiv_all$n_calib_label <- calib_label(df_fun_eiv_all$n_calib)
  
  p_fun_all_calib <- ggplot(df_fun_eiv_all, aes(x = u)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70", alpha = 0.38) +
    geom_line(aes(y = mean), color = "firebrick", linewidth = 0.75) +
    geom_line(aes(y = truth), color = "black", linewidth = 0.75, linetype = "dashed") +
    facet_grid(n_calib_label ~ x_slice) +
    labs(
      x = "Latent u",
      y = expression(f(x, u)),
      title = "EIV-GP response-surface slices across all calibration sizes",
      subtitle = "Red: posterior mean; grey: 95% interval; dashed black: truth"
    )
  
  safe_save_plot(
    file.path(FUN_FIG_DIR, "explore_function_slices_all_calibrations"),
    p_fun_all_calib,
    width = 16,
    height = max(7, 1.8 * length(rep_fit_calibs))
  )
  
  for (nm in names(rep_fits)) {
    df_one <- df_fun_eiv_all[df_fun_eiv_all$n_calib == as.integer(nm), ]
    
    p_one <- ggplot(df_one, aes(x = u)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70", alpha = 0.40) +
      geom_line(aes(y = mean), color = "firebrick", linewidth = 0.8) +
      geom_line(aes(y = truth), color = "black", linewidth = 0.8, linetype = "dashed") +
      facet_wrap(~x_slice, nrow = 1) +
      labs(
        x = "Latent u",
        y = expression(f(x, u)),
        title = paste0("EIV-GP response-surface slices, |O| = ", nm),
        subtitle = "Red: posterior mean; grey: 95% interval; dashed black: truth"
      )
    
    safe_save_plot(
      file.path(FUN_FIG_DIR, paste0("explore_function_slices_eiv_calib_", nm)),
      p_one,
      width = 16,
      height = 4
    )
  }
  
  make_eiv_function_surface_df <- function(fit,
                                           x_grid = seq(-1.75, 1.75, length.out = 60),
                                           u_grid = seq(-2.6, 2.6, length.out = 75),
                                           draw_ids = NULL,
                                           scenario = "active",
                                           max_draw = n_surface_draw,
                                           include_process_uncertainty = FALSE) {
    if (is.null(draw_ids)) {
      draw_ids <- seq_len(nrow(fit$mcmc$samples_u))
    }
    
    if (!is.null(max_draw) && length(draw_ids) > max_draw) {
      draw_ids <- sample(draw_ids, max_draw)
    }
    
    grid <- expand.grid(
      x_raw = x_grid,
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
      
      if (include_process_uncertainty) {
        f_std <- pred$mean + sqrt(pmax(pred$var, 0)) * rnorm(nrow(grid))
      } else {
        f_std <- pred$mean
      }
      
      f_samps[ii, ] <- fit$data$y_center + fit$data$y_scale * f_std
    }
    
    out <- grid
    out$mean <- colMeans(f_samps, na.rm = TRUE)
    out$lo <- apply(f_samps, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
    out$hi <- apply(f_samps, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
    out$width <- out$hi - out$lo
    out$truth <- f0_1d(out$x_raw, out$u, scenario = scenario)
    out$error <- out$mean - out$truth
    
    out
  }
  
  if (STUDY1_MAKE_SURFACE_HEATMAPS) {
    surface_file <- file.path(
      RES_DIR,
      paste0(
        "representative_eiv_surface_summaries_",
        "calibs_", calib_tag,
        "_ndraw", n_surface_draw,
        ".rds"
      )
    )
    
    if (STUDY1_USE_CACHE && file.exists(surface_file)) {
      df_surface_all <- readRDS(surface_file)
    } else {
      surface_x_grid <- seq(-1.75, 1.75, length.out = 60)
      surface_u_grid <- seq(-2.6, 2.6, length.out = 75)
      
      df_surface_all <- dplyr::bind_rows(
        lapply(names(rep_fits), function(nm) {
          out <- make_eiv_function_surface_df(
            rep_fits[[nm]],
            x_grid = surface_x_grid,
            u_grid = surface_u_grid,
            scenario = rep_dat$scenario,
            max_draw = n_surface_draw,
            include_process_uncertainty = FALSE
          )
          out$n_calib <- as.integer(nm)
          out$n_calib_label <- calib_label(as.integer(nm))
          out
        })
      )
      
      saveRDS(df_surface_all, surface_file)
    }
    
    df_surface_all$n_calib_label <- calib_label(df_surface_all$n_calib)
    
    p_surface_mean <- ggplot(df_surface_all, aes(x = x_raw, y = u, fill = mean)) +
      geom_raster(interpolate = TRUE) +
      geom_hline(yintercept = tau_true, color = "black", linewidth = 0.25) +
      facet_wrap(~n_calib_label, nrow = 1) +
      scale_fill_gradient2(
        low = "navy",
        mid = "white",
        high = "firebrick",
        midpoint = median(df_surface_all$truth, na.rm = TRUE),
        name = "Posterior mean"
      ) +
      labs(
        x = "Observed x",
        y = "Latent u",
        title = "EIV-GP posterior mean response surface across calibration sizes"
      )
    
    safe_save_plot(
      file.path(FUN_FIG_DIR, "explore_eiv_surface_posterior_mean_all_calibrations"),
      p_surface_mean,
      width = max(13, 3.1 * length(rep_fit_calibs)),
      height = 4.8
    )
    
    p_surface_error <- ggplot(df_surface_all, aes(x = x_raw, y = u, fill = error)) +
      geom_raster(interpolate = TRUE) +
      geom_hline(yintercept = tau_true, color = "black", linewidth = 0.25) +
      facet_wrap(~n_calib_label, nrow = 1) +
      scale_fill_gradient2(
        low = "navy",
        mid = "white",
        high = "firebrick",
        midpoint = 0,
        name = "Mean - truth"
      ) +
      labs(
        x = "Observed x",
        y = "Latent u",
        title = "EIV-GP response-surface error across calibration sizes"
      )
    
    safe_save_plot(
      file.path(FUN_FIG_DIR, "explore_eiv_surface_error_all_calibrations"),
      p_surface_error,
      width = max(13, 3.1 * length(rep_fit_calibs)),
      height = 4.8
    )
    
    p_surface_width <- ggplot(df_surface_all, aes(x = x_raw, y = u, fill = width)) +
      geom_raster(interpolate = TRUE) +
      geom_hline(yintercept = tau_true, color = "black", linewidth = 0.25) +
      facet_wrap(~n_calib_label, nrow = 1) +
      scale_fill_gradient(low = "white", high = "firebrick", name = "95% width") +
      labs(
        x = "Observed x",
        y = "Latent u",
        title = "EIV-GP response-surface posterior uncertainty"
      )
    
    safe_save_plot(
      file.path(FUN_FIG_DIR, "explore_eiv_surface_uncertainty_all_calibrations"),
      p_surface_width,
      width = max(13, 3.1 * length(rep_fit_calibs)),
      height = 4.8
    )
    
    df_surface_main <- df_surface_all[df_surface_all$n_calib == main_calib, ]
    
    p_surface_main_overlay <- ggplot(df_surface_main, aes(x = x_raw, y = u, fill = mean)) +
      geom_raster(interpolate = TRUE) +
      geom_contour(aes(z = mean), color = "white", linewidth = 0.25, alpha = 0.55) +
      geom_hline(yintercept = tau_true, color = "black", linewidth = 0.35) +
      geom_point(
        data = df_train,
        aes(x = x, y = u),
        inherit.aes = FALSE,
        shape = 21,
        fill = "white",
        color = "black",
        size = 0.85,
        alpha = 0.6
      ) +
      scale_fill_gradient2(
        low = "navy",
        mid = "white",
        high = "firebrick",
        midpoint = median(df_surface_main$truth, na.rm = TRUE),
        name = "Posterior mean"
      ) +
      labs(
        x = "Observed x",
        y = "Latent u",
        title = paste0("EIV-GP posterior mean surface with training inputs, |O| = ", main_calib)
      )
    
    safe_save_plot(
      file.path(FUN_FIG_DIR, "explore_eiv_surface_main_fit_with_training_points"),
      p_surface_main_overlay,
      width = 8,
      height = 6.2
    )
  }
}

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

n_density_eff <- length(draw_ids_density)

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
  n_draw = n_density_eff
)

draws_oracle_selected <- sample_oracle_test_y(
  x_test = selected_grid$x_star,
  c_test = selected_grid$c_star,
  tau_true = tau_true,
  scenario = rep_dat$scenario,
  sigma_eps = rep_dat$sigma_eps,
  n_draw = n_density_eff
)

baseline_method_order <- c("GP-LearnedEmb", "GP-CondMean", "GP-Gaussian")
baseline_methods_selected <- intersect(baseline_method_order, names(draws_baseline_selected))

density_parts <- list(
  make_long_draws(draws_oracle_selected, "Oracle", selected_grid),
  make_long_draws(draws_eiv_selected, "EIV-GP", selected_grid)
)
density_parts <- c(
  density_parts,
  lapply(baseline_methods_selected, function(mm) {
    make_long_draws(draws_baseline_selected[[mm]], mm, selected_grid)
  })
)

df_density <- dplyr::bind_rows(density_parts)

df_density$method <- factor(df_density$method, levels = method_levels)

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
  scale_color_manual(values = method_cols, name = NULL, drop = FALSE) +
  labs(
    x = expression(y^"*"),
    y = "Density",
    title = "Predictive distributions at selected mixed inputs"
  ) +
  theme(legend.position = "bottom")

save_plot(
  file.path(FIG_DIR, "fig4_study1_predictive_densities_selected"),
  p_density,
  width = 11,
  height = 7.5
)

############################################################
## Extra predictive-distribution plots over a larger grid
############################################################

if (STUDY1_RESEARCH_PLOTS) {
  dense_x_vals <- seq(-1.5, 1.5, by = 0.5)
  dense_c_vals <- seq_len(m)
  
  predictive_grid <- expand.grid(
    x_star = dense_x_vals,
    c_star = dense_c_vals
  )
  
  pred_grid_file <- file.path(
    RES_DIR,
    paste0(
      "representative_predictive_grid_draws_",
      "calib", main_calib,
      "_ndraw", n_density_draw_research,
      ".rds"
    )
  )
  
  if (STUDY1_USE_CACHE && file.exists(pred_grid_file)) {
    pred_grid_obj <- readRDS(pred_grid_file)
    predictive_grid <- pred_grid_obj$grid
    pred_grid_draws <- pred_grid_obj$draws
  } else {
    draw_ids_density_grid <- seq_len(nrow(main_fit$mcmc$samples_u))
    if (length(draw_ids_density_grid) > n_density_draw_research) {
      draw_ids_density_grid <- sample(draw_ids_density_grid, n_density_draw_research)
    }
    
    n_density_grid_eff <- length(draw_ids_density_grid)
    
    draws_eiv_grid <- sample_eiv_test_y(
      x_test_raw = predictive_grid$x_star,
      c_test = predictive_grid$c_star,
      fit_obj = main_fit,
      draw_ids = draw_ids_density_grid,
      n_per_draw = 1L
    )
    
    draws_baseline_grid <- predict_embedding_baseline_samples(
      baselines = rep_baselines,
      x_star_raw = predictive_grid$x_star,
      c_star = predictive_grid$c_star,
      m = m,
      n_draw = n_density_grid_eff
    )
    
    draws_oracle_grid <- sample_oracle_test_y(
      x_test = predictive_grid$x_star,
      c_test = predictive_grid$c_star,
      tau_true = tau_true,
      scenario = rep_dat$scenario,
      sigma_eps = rep_dat$sigma_eps,
      n_draw = n_density_grid_eff
    )
    
    pred_grid_draws <- c(
      stats::setNames(list(draws_oracle_grid), "Oracle"),
      stats::setNames(list(draws_eiv_grid), "EIV-GP"),
      draws_baseline_grid
    )
    
    pred_grid_obj <- list(
      grid = predictive_grid,
      draws = pred_grid_draws,
      draw_ids_eiv = draw_ids_density_grid
    )
    
    saveRDS(pred_grid_obj, pred_grid_file)
  }
  
  grid_methods <- intersect(method_levels, names(pred_grid_draws))
  
  df_density_all <- dplyr::bind_rows(
    lapply(grid_methods, function(mm) {
      make_long_draws(pred_grid_draws[[mm]], mm, predictive_grid)
    })
  )
  
  df_density_all$method <- factor(df_density_all$method, levels = method_levels)
  df_density_all$x_label <- factor(
    paste0("x* = ", df_density_all$x_star),
    levels = paste0("x* = ", dense_x_vals)
  )
  df_density_all$c_label <- factor(
    paste0("c* = ", df_density_all$c_star),
    levels = paste0("c* = ", dense_c_vals)
  )
  
  p_density_all <- ggplot(df_density_all, aes(x = y, color = method)) +
    geom_density(linewidth = 0.72) +
    facet_grid(c_label ~ x_label, scales = "free_y") +
    scale_color_manual(values = method_cols, name = NULL, drop = FALSE) +
    labs(
      x = expression(y^"*"),
      y = "Density",
      title = "Predictive distributions over a larger x-by-class grid"
    ) +
    theme(legend.position = "bottom")
  
  safe_save_plot(
    file.path(PRED_FIG_DIR, "explore_predictive_densities_full_grid"),
    p_density_all,
    width = 18,
    height = 12.5
  )
  
  df_pred_grid_summary <- dplyr::bind_rows(
    lapply(grid_methods, function(mm) {
      summarise_draw_matrix(pred_grid_draws[[mm]], mm, predictive_grid)
    })
  )
  
  df_pred_grid_summary$method <- factor(df_pred_grid_summary$method, levels = method_levels)
  df_pred_grid_summary$x_label <- factor(
    paste0("x* = ", df_pred_grid_summary$x_star),
    levels = paste0("x* = ", dense_x_vals)
  )
  df_pred_grid_summary$c_label <- factor(
    paste0("c* = ", df_pred_grid_summary$c_star),
    levels = paste0("c* = ", dense_c_vals)
  )
  
  save_csv(
    df_pred_grid_summary,
    file.path(RES_DIR, "representative_predictive_grid_summary.csv")
  )
  
  p_pred_mean_grid <- ggplot(
    df_pred_grid_summary,
    aes(x = x_star, y = mean, color = method, linetype = method)
  ) +
    geom_line(linewidth = 0.85) +
    geom_point(size = 1.5) +
    facet_wrap(~c_label, ncol = 3, scales = "free_y") +
    scale_color_manual(values = method_cols, name = NULL, drop = FALSE) +
    labs(
      x = expression(x^"*"),
      y = "Predictive mean",
      title = "Predictive mean over x for each ordinal class"
    ) +
    theme(legend.position = "bottom")
  
  safe_save_plot(
    file.path(PRED_FIG_DIR, "explore_predictive_mean_lines_by_class"),
    p_pred_mean_grid,
    width = 12,
    height = 7.5
  )
  
  p_pred_interval_grid <- ggplot(
    df_pred_grid_summary,
    aes(x = x_label, y = mean, color = method)
  ) +
    geom_pointrange(
      aes(ymin = lo95, ymax = hi95),
      position = position_dodge(width = 0.65),
      alpha = 0.78,
      linewidth = 0.5
    ) +
    facet_wrap(~c_label, ncol = 3, scales = "free_y") +
    scale_color_manual(values = method_cols, name = NULL, drop = FALSE) +
    labs(
      x = expression(x^"*"),
      y = "Predictive mean and 95% interval",
      title = "Predictive intervals over x-by-class grid"
    ) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  safe_save_plot(
    file.path(PRED_FIG_DIR, "explore_predictive_intervals_grid"),
    p_pred_interval_grid,
    width = 13,
    height = 8
  )
  
  p_pred_sd_heat <- ggplot(df_pred_grid_summary, aes(x = x_label, y = c_label, fill = sd)) +
    geom_tile(color = "white") +
    facet_wrap(~method) +
    scale_fill_gradient(low = "white", high = "firebrick", name = "SD") +
    labs(
      x = expression(x^"*"),
      y = expression(c^"*"),
      title = "Predictive standard deviation heatmaps"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  safe_save_plot(
    file.path(PRED_FIG_DIR, "explore_predictive_sd_heatmaps"),
    p_pred_sd_heat,
    width = 13,
    height = 7
  )
  
  p_pred_width_heat <- ggplot(df_pred_grid_summary, aes(x = x_label, y = c_label, fill = width95)) +
    geom_tile(color = "white") +
    facet_wrap(~method) +
    scale_fill_gradient(low = "white", high = "firebrick", name = "95% width") +
    labs(
      x = expression(x^"*"),
      y = expression(c^"*"),
      title = "Predictive 95% interval width heatmaps"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  safe_save_plot(
    file.path(PRED_FIG_DIR, "explore_predictive_width95_heatmaps"),
    p_pred_width_heat,
    width = 13,
    height = 7
  )
  
  if ("Oracle" %in% as.character(df_pred_grid_summary$method)) {
    df_oracle_ref <- df_pred_grid_summary |>
      dplyr::filter(method == "Oracle") |>
      dplyr::select(
        x_star,
        c_star,
        oracle_mean = mean,
        oracle_sd = sd,
        oracle_width95 = width95
      )
    
    df_pred_compare <- df_pred_grid_summary |>
      dplyr::filter(method != "Oracle") |>
      dplyr::left_join(df_oracle_ref, by = c("x_star", "c_star"))
    
    df_pred_compare$mean_bias <- df_pred_compare$mean - df_pred_compare$oracle_mean
    df_pred_compare$sd_ratio <- df_pred_compare$sd / df_pred_compare$oracle_sd
    df_pred_compare$width95_ratio <- df_pred_compare$width95 / df_pred_compare$oracle_width95
    
    save_csv(
      df_pred_compare,
      file.path(RES_DIR, "representative_predictive_grid_comparison_to_oracle.csv")
    )
    
    p_pred_bias_heat <- ggplot(df_pred_compare, aes(x = x_label, y = c_label, fill = mean_bias)) +
      geom_tile(color = "white") +
      facet_wrap(~method) +
      scale_fill_gradient2(
        low = "navy",
        mid = "white",
        high = "firebrick",
        midpoint = 0,
        name = "Mean - oracle"
      ) +
      labs(
        x = expression(x^"*"),
        y = expression(c^"*"),
        title = "Predictive mean bias relative to oracle"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    safe_save_plot(
      file.path(PRED_FIG_DIR, "explore_predictive_mean_bias_vs_oracle_heatmaps"),
      p_pred_bias_heat,
      width = 13,
      height = 6.5
    )
    
    p_pred_sd_ratio_heat <- ggplot(df_pred_compare, aes(x = x_label, y = c_label, fill = sd_ratio)) +
      geom_tile(color = "white") +
      facet_wrap(~method) +
      scale_fill_gradient2(
        low = "navy",
        mid = "white",
        high = "firebrick",
        midpoint = 1,
        name = "SD ratio"
      ) +
      labs(
        x = expression(x^"*"),
        y = expression(c^"*"),
        title = "Predictive SD ratio relative to oracle"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    safe_save_plot(
      file.path(PRED_FIG_DIR, "explore_predictive_sd_ratio_vs_oracle_heatmaps"),
      p_pred_sd_ratio_heat,
      width = 13,
      height = 6.5
    )
  }
}

############################################################
## Extra test-set predictive diagnostics
############################################################

if (
  isTRUE(STUDY1_MAKE_TEST_PREDICTIONS) &&
  !is.null(test) &&
  all(c("x", "c", "y") %in% names(test))
) {
  test_grid <- data.frame(
    id = seq_along(test$x),
    x_star = test$x,
    c_star = test$c
  )
  
  test_pred_file <- file.path(
    RES_DIR,
    paste0(
      "representative_test_predictive_draws_",
      "calib", main_calib,
      "_ndraw", n_test_pred_draw,
      ".rds"
    )
  )
  
  use_test_cache <- FALSE
  if (STUDY1_USE_CACHE && file.exists(test_pred_file)) {
    pred_test_obj <- readRDS(test_pred_file)
    if (!is.null(pred_test_obj$grid) && nrow(pred_test_obj$grid) == nrow(test_grid)) {
      use_test_cache <- TRUE
    }
  }
  
  if (!use_test_cache) {
    draw_ids_test <- seq_len(nrow(main_fit$mcmc$samples_u))
    if (length(draw_ids_test) > n_test_pred_draw) {
      draw_ids_test <- sample(draw_ids_test, n_test_pred_draw)
    }
    
    n_test_eff <- length(draw_ids_test)
    
    draws_eiv_test <- sample_eiv_test_y(
      x_test_raw = test_grid$x_star,
      c_test = test_grid$c_star,
      fit_obj = main_fit,
      draw_ids = draw_ids_test,
      n_per_draw = 1L
    )
    
    draws_baseline_test <- predict_embedding_baseline_samples(
      baselines = rep_baselines,
      x_star_raw = test_grid$x_star,
      c_star = test_grid$c_star,
      m = m,
      n_draw = n_test_eff
    )
    
    draws_oracle_test <- sample_oracle_test_y(
      x_test = test_grid$x_star,
      c_test = test_grid$c_star,
      tau_true = tau_true,
      scenario = rep_dat$scenario,
      sigma_eps = rep_dat$sigma_eps,
      n_draw = n_test_eff
    )
    
    pred_test_draws <- c(
      stats::setNames(list(draws_oracle_test), "Oracle"),
      stats::setNames(list(draws_eiv_test), "EIV-GP"),
      draws_baseline_test
    )
    
    pred_test_obj <- list(
      grid = test_grid,
      draws = pred_test_draws,
      draw_ids_eiv = draw_ids_test
    )
    
    saveRDS(pred_test_obj, test_pred_file)
  }
  
  pred_test_draws <- pred_test_obj$draws
  test_methods <- intersect(method_levels, names(pred_test_draws))
  
  df_test_summary <- dplyr::bind_rows(
    lapply(test_methods, function(mm) {
      summarise_draw_matrix(pred_test_draws[[mm]], mm, test_grid)
    })
  )
  
  df_test_summary$method <- factor(df_test_summary$method, levels = method_levels)
  df_test_summary$y_true <- test$y[df_test_summary$id]
  df_test_summary$x <- df_test_summary$x_star
  df_test_summary$c <- factor(df_test_summary$c_star, levels = seq_len(m))
  df_test_summary$resid <- df_test_summary$mean - df_test_summary$y_true
  df_test_summary$abs_error <- abs(df_test_summary$resid)
  df_test_summary$covered95 <- (
    df_test_summary$y_true >= df_test_summary$lo95 &
      df_test_summary$y_true <= df_test_summary$hi95
  )
  df_test_summary$u_true <- if ("u" %in% names(test)) {
    test$u[df_test_summary$id]
  } else {
    NA_real_
  }
  df_test_summary$f_true <- if ("u" %in% names(test)) {
    f0_1d(
      test$x[df_test_summary$id],
      test$u[df_test_summary$id],
      scenario = rep_dat$scenario
    )
  } else {
    NA_real_
  }
  
  save_csv(
    df_test_summary,
    file.path(RES_DIR, "representative_test_predictive_summary.csv")
  )
  
  df_test_metrics <- df_test_summary |>
    dplyr::group_by(method) |>
    dplyr::summarise(
      n = dplyr::n(),
      bias = base::mean(.data$resid, na.rm = TRUE),
      rmse = sqrt(base::mean(.data$resid^2, na.rm = TRUE)),
      mae = base::mean(.data$abs_error, na.rm = TRUE),
      coverage95 = base::mean(.data$covered95, na.rm = TRUE),
      mean_width95 = base::mean(.data$width95, na.rm = TRUE),
      median_width95 = stats::median(.data$width95, na.rm = TRUE),
      avg_pred_sd = base::mean(.data$sd, na.rm = TRUE),
      rmse_latent_truth = if (all(is.na(.data$f_true))) {
        NA_real_
      } else {
        sqrt(base::mean((.data$mean - .data$f_true)^2, na.rm = TRUE))
      },
      .groups = "drop"
    )
  
  save_csv(
    df_test_metrics,
    file.path(RES_DIR, "representative_test_predictive_metrics.csv")
  )
  
  nominal_levels <- c(0.50, 0.80, 0.90, 0.95)
  
  df_nominal <- dplyr::bind_rows(
    lapply(test_methods, function(mm) {
      mat <- as.matrix(pred_test_draws[[mm]])
      
      dplyr::bind_rows(
        lapply(nominal_levels, function(pp) {
          alpha <- 1 - pp
          lo <- apply(mat, 2, stats::quantile, probs = alpha / 2, na.rm = TRUE)
          hi <- apply(mat, 2, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE)
          
          data.frame(
            method = mm,
            nominal = pp,
            coverage = mean(test$y >= lo & test$y <= hi),
            mean_width = mean(hi - lo)
          )
        })
      )
    })
  )
  
  df_nominal$method <- factor(df_nominal$method, levels = method_levels)
  
  save_csv(
    df_nominal,
    file.path(RES_DIR, "representative_test_nominal_coverage.csv")
  )
  
  df_test_metrics_long <- df_test_metrics |>
    tidyr::pivot_longer(
      cols = c(
        bias,
        rmse,
        mae,
        coverage95,
        mean_width95,
        avg_pred_sd,
        rmse_latent_truth
      ),
      names_to = "metric",
      values_to = "value"
    )
  
  p_test_metrics <- ggplot(
    df_test_metrics_long[!is.na(df_test_metrics_long$value), ],
    aes(x = method, y = value, fill = method)
  ) +
    geom_col(show.legend = FALSE, alpha = 0.88) +
    facet_wrap(~metric, scales = "free_y", ncol = 3) +
    coord_flip() +
    scale_fill_manual(values = method_cols, drop = FALSE) +
    labs(
      x = NULL,
      y = NULL,
      title = "Test-set predictive metrics"
    )
  
  safe_save_plot(
    file.path(TEST_PRED_FIG_DIR, "explore_test_predictive_metrics"),
    p_test_metrics,
    width = 11,
    height = 7
  )
  
  p_test_pred_obs <- ggplot(df_test_summary, aes(x = y_true, y = mean, color = c)) +
    geom_point(alpha = 0.55, size = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    facet_wrap(~method, scales = "free") +
    scale_color_manual(values = class_cols, name = "Class") +
    labs(
      x = "Observed test y",
      y = "Predictive mean",
      title = "Predictive mean versus observed test response"
    ) +
    theme(legend.position = "bottom")
  
  safe_save_plot(
    file.path(TEST_PRED_FIG_DIR, "explore_test_predicted_vs_observed"),
    p_test_pred_obs,
    width = 12,
    height = 7
  )
  
  p_test_resid_x <- ggplot(df_test_summary, aes(x = x, y = resid, color = c)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey35") +
    geom_point(alpha = 0.55, size = 1.1) +
    geom_smooth(se = FALSE, method = "loess", formula = y ~ x, linewidth = 0.6) +
    facet_wrap(~method, scales = "free_y") +
    scale_color_manual(values = class_cols, name = "Class") +
    labs(
      x = "Observed x",
      y = "Predictive mean - observed y",
      title = "Test residuals versus observed x"
    ) +
    theme(legend.position = "bottom")
  
  safe_save_plot(
    file.path(TEST_PRED_FIG_DIR, "explore_test_residuals_vs_x"),
    p_test_resid_x,
    width = 12,
    height = 7
  )
  
  p_test_width <- ggplot(df_test_summary, aes(x = c, y = width95, fill = c)) +
    geom_boxplot(alpha = 0.82, outlier.alpha = 0.35) +
    facet_wrap(~method, scales = "free_y") +
    scale_fill_manual(values = class_cols, name = "Class") +
    labs(
      x = "Class",
      y = "95% predictive interval width",
      title = "Test-set predictive interval width by class"
    ) +
    theme(legend.position = "none")
  
  safe_save_plot(
    file.path(TEST_PRED_FIG_DIR, "explore_test_predictive_width_by_class"),
    p_test_width,
    width = 12,
    height = 7
  )
  
  p_nominal <- ggplot(df_nominal, aes(x = nominal, y = coverage, color = method)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.2) +
    scale_color_manual(values = method_cols, name = NULL, drop = FALSE) +
    coord_equal(xlim = c(0.45, 1), ylim = c(0.45, 1)) +
    labs(
      x = "Nominal predictive interval level",
      y = "Empirical coverage",
      title = "Test-set predictive coverage calibration"
    ) +
    theme(legend.position = "bottom")
  
  safe_save_plot(
    file.path(TEST_PRED_FIG_DIR, "explore_test_nominal_coverage"),
    p_nominal,
    width = 7,
    height = 6.2
  )
  
  p_test_abs_error_class <- ggplot(df_test_summary, aes(x = method, y = abs_error, fill = method)) +
    geom_boxplot(alpha = 0.82, outlier.alpha = 0.30) +
    facet_wrap(~c, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = method_cols, name = NULL, drop = FALSE) +
    coord_flip() +
    labs(
      x = NULL,
      y = "Absolute predictive-mean error",
      title = "Test absolute error by class"
    ) +
    theme(legend.position = "none")
  
  safe_save_plot(
    file.path(TEST_PRED_FIG_DIR, "explore_test_abs_error_by_class"),
    p_test_abs_error_class,
    width = 12,
    height = 7
  )
  
  p_test_resid_density <- ggplot(df_test_summary, aes(x = resid, color = method)) +
    geom_density(linewidth = 0.85) +
    facet_wrap(~c, scales = "free_y", ncol = 3) +
    scale_color_manual(values = method_cols, name = NULL, drop = FALSE) +
    labs(
      x = "Predictive mean - observed y",
      y = "Density",
      title = "Test residual-density comparison by class"
    ) +
    theme(legend.position = "bottom")
  
  safe_save_plot(
    file.path(TEST_PRED_FIG_DIR, "explore_test_residual_density_by_class"),
    p_test_resid_density,
    width = 12,
    height = 7
  )
}

############################################################
## Figure 6: MCMC trace plots for representative fit
############################################################

make_trace_df <- function(fit, n_calib) {
  samples_by_chain <- fit$mcmc$samples_by_chain
  
  if (
    is.null(samples_by_chain) ||
    is.null(samples_by_chain$sigma2) ||
    is.null(samples_by_chain$logtheta)
  ) {
    return(data.frame())
  }
  
  dplyr::bind_rows(
    lapply(seq_along(samples_by_chain$sigma2), function(cc) {
      logtheta_mat <- as.matrix(samples_by_chain$logtheta[[cc]])
      n_draw <- length(samples_by_chain$sigma2[[cc]])
      
      data.frame(
        chain = factor(cc),
        draw = seq_len(n_draw),
        sigma_epsilon = sqrt(samples_by_chain$sigma2[[cc]]),
        rho = if (ncol(logtheta_mat) >= 1L) exp(logtheta_mat[, 1]) else NA_real_,
        theta_x = if (ncol(logtheta_mat) >= 2L) exp(logtheta_mat[, 2]) else NA_real_,
        theta_u = if (ncol(logtheta_mat) >= 3L) exp(logtheta_mat[, 3]) else NA_real_,
        n_calib = n_calib,
        n_calib_label = calib_label(n_calib)
      )
    })
  )
}

df_trace <- make_trace_df(main_fit, main_calib)

trace_param_cols <- intersect(
  c("sigma_epsilon", "rho", "theta_x", "theta_u"),
  names(df_trace)
)

df_trace_long <- df_trace |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(trace_param_cols),
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
    title = paste0("MCMC trace plots for representative EIV-GP fit, |O| = ", main_calib)
  ) +
  theme(legend.position = "bottom")

save_plot(
  file.path(FIG_DIR, "fig6_study1_mcmc_traces"),
  p_trace,
  width = 10,
  height = 6
)

############################################################
## Extra MCMC diagnostics
############################################################

if (STUDY1_RESEARCH_PLOTS) {
  df_trace_all <- dplyr::bind_rows(
    lapply(names(rep_fits), function(nm) {
      make_trace_df(rep_fits[[nm]], as.integer(nm))
    })
  )
  
  if (nrow(df_trace_all) > 0L) {
    df_trace_all_long <- df_trace_all |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(trace_param_cols),
        names_to = "parameter",
        values_to = "value"
      )
    
    p_trace_all <- ggplot(df_trace_all_long, aes(x = draw, y = value, color = chain)) +
      geom_line(linewidth = 0.25, alpha = 0.65) +
      facet_grid(parameter ~ n_calib_label, scales = "free_y") +
      labs(
        x = "Saved draw within chain",
        y = NULL,
        color = "Chain",
        title = "MCMC trace plots across calibration sizes"
      ) +
      theme(legend.position = "bottom")
    
    safe_save_plot(
      file.path(MCMC_FIG_DIR, "explore_mcmc_traces_all_calibrations"),
      p_trace_all,
      width = max(13, 3.1 * length(rep_fit_calibs)),
      height = 9
    )
  }
  
  make_hyper_draw_df <- function(fit, n_calib) {
    logtheta_mat <- as.matrix(fit$mcmc$samples_logtheta)
    
    data.frame(
      draw = seq_len(nrow(logtheta_mat)),
      n_calib = n_calib,
      n_calib_label = calib_label(n_calib),
      sigma_epsilon = sqrt(fit$mcmc$samples_sigma2),
      rho = if (ncol(logtheta_mat) >= 1L) exp(logtheta_mat[, 1]) else NA_real_,
      theta_x = if (ncol(logtheta_mat) >= 2L) exp(logtheta_mat[, 2]) else NA_real_,
      theta_u = if (ncol(logtheta_mat) >= 3L) exp(logtheta_mat[, 3]) else NA_real_
    )
  }
  
  df_hyper_all <- dplyr::bind_rows(
    lapply(names(rep_fits), function(nm) {
      make_hyper_draw_df(rep_fits[[nm]], as.integer(nm))
    })
  )
  
  df_hyper_long <- df_hyper_all |>
    tidyr::pivot_longer(
      cols = c(sigma_epsilon, rho, theta_x, theta_u),
      names_to = "parameter",
      values_to = "value"
    )
  
  p_hyper_density <- ggplot(
    df_hyper_long,
    aes(x = value, color = n_calib_label, fill = n_calib_label)
  ) +
    geom_density(alpha = 0.10, linewidth = 0.8, na.rm = TRUE) +
    facet_wrap(~parameter, scales = "free", ncol = 2) +
    labs(
      x = NULL,
      y = "Density",
      color = "Calibration",
      fill = "Calibration",
      title = "Posterior hyperparameter densities across calibration sizes"
    ) +
    theme(legend.position = "bottom")
  
  safe_save_plot(
    file.path(MCMC_FIG_DIR, "explore_mcmc_hyperparameter_posterior_densities"),
    p_hyper_density,
    width = 11,
    height = 7
  )
  
  df_hyper_summary <- df_hyper_long |>
    dplyr::group_by(n_calib, n_calib_label, parameter) |>
    dplyr::summarise(
      mean = base::mean(value, na.rm = TRUE),
      median = stats::median(value, na.rm = TRUE),
      lo = stats::quantile(value, 0.025, na.rm = TRUE),
      hi = stats::quantile(value, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  save_csv(
    df_hyper_summary,
    file.path(RES_DIR, "representative_hyperparameter_posterior_summary.csv")
  )
  
  p_hyper_intervals <- ggplot(
    df_hyper_summary,
    aes(x = n_calib, y = median, ymin = lo, ymax = hi)
  ) +
    geom_line(color = "grey30", linewidth = 0.65) +
    geom_pointrange(color = "firebrick", linewidth = 0.55) +
    facet_wrap(~parameter, scales = "free_y", ncol = 2) +
    labs(
      x = "Calibration size |O|",
      y = "Posterior median and 95% interval",
      title = "Hyperparameter posterior intervals versus calibration size"
    )
  
  safe_save_plot(
    file.path(MCMC_FIG_DIR, "explore_mcmc_hyperparameter_intervals_by_calibration"),
    p_hyper_intervals,
    width = 10,
    height = 7
  )
  
  make_acf_df <- function(df_long, max_lag = 60L) {
    pieces <- split(df_long, list(df_long$parameter, df_long$chain), drop = TRUE)
    
    dplyr::bind_rows(
      lapply(names(pieces), function(nm) {
        dd <- pieces[[nm]]
        aa <- stats::acf(
          dd$value,
          lag.max = max_lag,
          plot = FALSE,
          na.action = stats::na.pass
        )
        
        data.frame(
          parameter = unique(dd$parameter)[1],
          chain = unique(dd$chain)[1],
          lag = seq_along(as.numeric(aa$acf)) - 1L,
          acf = as.numeric(aa$acf)
        )
      })
    )
  }
  
  df_acf_main <- make_acf_df(df_trace_long, max_lag = 60L)
  
  p_acf_main <- ggplot(df_acf_main, aes(x = lag, y = acf, color = chain)) +
    geom_hline(yintercept = 0, color = "grey40") +
    geom_line(linewidth = 0.5, alpha = 0.85) +
    facet_wrap(~parameter, scales = "free_y", ncol = 2) +
    labs(
      x = "Lag",
      y = "Autocorrelation",
      color = "Chain",
      title = paste0("MCMC autocorrelation diagnostics, |O| = ", main_calib)
    ) +
    theme(legend.position = "bottom")
  
  safe_save_plot(
    file.path(MCMC_FIG_DIR, "explore_mcmc_acf_main_fit"),
    p_acf_main,
    width = 10,
    height = 6.5
  )
  
  df_running <- df_trace_long |>
    dplyr::group_by(chain, parameter) |>
    dplyr::mutate(running_mean = cumsum(value) / seq_along(value)) |>
    dplyr::ungroup()
  
  p_running <- ggplot(df_running, aes(x = draw, y = running_mean, color = chain)) +
    geom_line(linewidth = 0.45, alpha = 0.8) +
    facet_wrap(~parameter, scales = "free_y", ncol = 2) +
    labs(
      x = "Saved draw within chain",
      y = "Running mean",
      color = "Chain",
      title = paste0("MCMC running means, |O| = ", main_calib)
    ) +
    theme(legend.position = "bottom")
  
  safe_save_plot(
    file.path(MCMC_FIG_DIR, "explore_mcmc_running_means_main_fit"),
    p_running,
    width = 10,
    height = 6.5
  )
  
  make_u_trace_df <- function(fit, ids) {
    samples_by_chain <- fit$mcmc$samples_by_chain
    if (is.null(samples_by_chain) || is.null(samples_by_chain$u)) {
      return(data.frame())
    }
    
    ids <- ids[ids >= 1L & ids <= ncol(fit$mcmc$samples_u)]
    if (length(ids) == 0L) return(data.frame())
    
    dplyr::bind_rows(
      lapply(seq_along(samples_by_chain$u), function(cc) {
        u_mat <- as.matrix(samples_by_chain$u[[cc]])
        
        dplyr::bind_rows(
          lapply(ids, function(id) {
            data.frame(
              chain = factor(cc),
              draw = seq_len(nrow(u_mat)),
              id = id,
              id_label = paste0("id=", id, ", c=", fit$data$c_ord[id]),
              u = u_mat[, id],
              true_u = fit$data$u_true[id]
            )
          })
        )
      })
    )
  }
  
  selected_trace_ids <- select_ids_by_quantile(
    main_fit$data$miss_idx,
    main_fit$data$u_true[main_fit$data$miss_idx],
    n = 6L
  )
  
  df_u_trace <- make_u_trace_df(main_fit, selected_trace_ids)
  
  if (nrow(df_u_trace) > 0L) {
    p_u_trace <- ggplot(df_u_trace, aes(x = draw, y = u, color = chain)) +
      geom_line(linewidth = 0.35, alpha = 0.78) +
      geom_hline(
        data = unique(df_u_trace[, c("id_label", "true_u")]),
        aes(yintercept = true_u),
        inherit.aes = FALSE,
        linetype = "dashed",
        color = "black"
      ) +
      facet_wrap(~id_label, scales = "free_y", ncol = 3) +
      labs(
        x = "Saved draw within chain",
        y = "Latent u",
        color = "Chain",
        title = paste0("Selected latent-u MCMC traces, |O| = ", main_calib)
      ) +
      theme(legend.position = "bottom")
    
    safe_save_plot(
      file.path(MCMC_FIG_DIR, "explore_mcmc_selected_u_traces_main_fit"),
      p_u_trace,
      width = 12,
      height = 7
    )
  }
  
  make_tau_trace_df <- function(fit) {
    samples_by_chain <- fit$mcmc$samples_by_chain
    if (is.null(samples_by_chain)) return(data.frame())
    
    tau_list <- NULL
    for (nm in c("tau", "taus", "cutpoints", "alpha")) {
      if (!is.null(samples_by_chain[[nm]])) {
        tau_list <- samples_by_chain[[nm]]
        break
      }
    }
    
    if (is.null(tau_list)) return(data.frame())
    
    dplyr::bind_rows(
      lapply(seq_along(tau_list), function(cc) {
        obj <- tau_list[[cc]]
        
        if (is.null(dim(obj))) {
          if (length(tau_true) > 1L && length(obj) %% length(tau_true) == 0L) {
            mat <- matrix(as.numeric(obj), ncol = length(tau_true), byrow = TRUE)
          } else {
            mat <- matrix(as.numeric(obj), ncol = 1L)
          }
        } else {
          mat <- as.matrix(obj)
        }
        
        if (ncol(mat) != length(tau_true) && nrow(mat) == length(tau_true)) {
          mat <- t(mat)
        }
        
        dplyr::bind_rows(
          lapply(seq_len(ncol(mat)), function(jj) {
            data.frame(
              chain = factor(cc),
              draw = seq_len(nrow(mat)),
              tau_index = factor(jj, levels = seq_len(length(tau_true))),
              tau = mat[, jj],
              tau_true = tau_true[jj]
            )
          })
        )
      })
    )
  }
  
  df_tau_trace <- make_tau_trace_df(main_fit)
  
  if (nrow(df_tau_trace) > 0L) {
    p_tau_trace <- ggplot(df_tau_trace, aes(x = draw, y = tau, color = chain)) +
      geom_line(linewidth = 0.35, alpha = 0.78) +
      geom_hline(
        data = unique(df_tau_trace[, c("tau_index", "tau_true")]),
        aes(yintercept = tau_true),
        inherit.aes = FALSE,
        linetype = "dashed",
        color = "black"
      ) +
      facet_wrap(~tau_index, scales = "free_y") +
      labs(
        x = "Saved draw within chain",
        y = "Cutpoint",
        color = "Chain",
        title = paste0("Cutpoint MCMC traces, |O| = ", main_calib)
      ) +
      theme(legend.position = "bottom")
    
    safe_save_plot(
      file.path(MCMC_FIG_DIR, "explore_mcmc_tau_traces_main_fit"),
      p_tau_trace,
      width = 11,
      height = 6.5
    )
  }
}

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

mcmc_summary_all <- dplyr::bind_rows(
  lapply(names(rep_fits), function(nm) {
    out <- rep_fits[[nm]]$diagnostics$summary
    out$n_calib <- as.integer(nm)
    out
  })
)

if (nrow(mcmc_summary_all) > 0L) {
  mcmc_summary_all <- mcmc_summary_all |>
    dplyr::select(n_calib, dplyr::everything())
  
  save_csv(
    mcmc_summary_all,
    file.path(RES_DIR, "representative_mcmc_summary_all_calibrations.csv")
  )
  
  latex_mcmc_all_table <- knitr::kable(
    mcmc_summary_all,
    format = "latex",
    booktabs = TRUE,
    escape = TRUE
  )
  
  writeLines(
    latex_mcmc_all_table,
    con = file.path(TAB_DIR, "study1_mcmc_summary_all_calibrations.tex")
  )
}

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

class_counts_all <- dplyr::bind_rows(
  lapply(rep_fit_calibs, function(kk) {
    data.frame(
      n_calib = kk,
      class = seq_len(m),
      training_count = as.integer(tabulate(train$c, nbins = m)),
      calibrated_count = as.integer(tabulate(
        train$c[calib_sets[[as.character(kk)]]],
        nbins = m
      ))
    )
  })
)

save_csv(
  class_counts_all,
  file.path(RES_DIR, "representative_class_counts_all_calibrations.csv")
)

if (exists("df_imp_metrics_overall") && nrow(df_imp_metrics_overall) > 0L) {
  imp_table <- df_imp_metrics_overall
  num_cols <- vapply(imp_table, is.numeric, logical(1))
  imp_table[num_cols] <- lapply(imp_table[num_cols], function(z) sprintf("%.3f", z))
  
  latex_imp_table <- knitr::kable(
    imp_table,
    format = "latex",
    booktabs = TRUE,
    escape = TRUE
  )
  
  writeLines(
    latex_imp_table,
    con = file.path(TAB_DIR, "study1_imputation_metrics_overall.tex")
  )
}

if (exists("df_test_metrics") && nrow(df_test_metrics) > 0L) {
  test_metric_table <- df_test_metrics
  num_cols <- vapply(test_metric_table, is.numeric, logical(1))
  test_metric_table[num_cols] <- lapply(test_metric_table[num_cols], function(z) {
    sprintf("%.3f", z)
  })
  
  latex_test_metric_table <- knitr::kable(
    test_metric_table,
    format = "latex",
    booktabs = TRUE,
    escape = TRUE
  )
  
  writeLines(
    latex_test_metric_table,
    con = file.path(TAB_DIR, "study1_test_predictive_metrics.tex")
  )
}

cat("\nRepresentative Study I figures written to:\n")
cat(normalizePath(FIG_DIR), "\n")

if (STUDY1_RESEARCH_PLOTS) {
  cat("\nExploratory Study I research figures written to:\n")
  cat(normalizePath(EXP_FIG_DIR), "\n")
}

cat("\nRepresentative Study I tables written to:\n")
cat(normalizePath(TAB_DIR), "\n")

cat("\nRepresentative Study I results and plot-data summaries written to:\n")
cat(normalizePath(RES_DIR), "\n")