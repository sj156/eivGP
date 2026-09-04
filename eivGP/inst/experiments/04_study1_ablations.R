############################################################
## 04_study1_ablations.R
##
## Response-free measurement fit, plug-in GP (PI-GP), and complete-case GP
## (CC-GP) for Study I. These are scientific ablations, not literature
## competitors. The deterministic-threshold Gibbs sampler below targets
## p(tau, U_missing | C, U_calibrated) and never receives y.
############################################################

if (!exists("rtruncnorm_vec") || !exists("gp_mle_fit_1d")) {
  source("model_univariate.R")
}

initial_thresholds_response_free <- function(c_ord,
                                             u_obs,
                                             calib_idx,
                                             m,
                                             tau_bounds = c(-8, 8)) {
  counts <- tabulate(c_ord, nbins = m)
  target <- stats::qnorm(
    pmin(pmax(cumsum(counts)[seq_len(m - 1L)] / sum(counts), 1e-4), 1 - 1e-4)
  )
  calibrated <- rep(FALSE, length(c_ord))
  calibrated[calib_idx] <- TRUE
  tau <- numeric(m - 1L)

  for (j in seq_len(m - 1L)) {
    lower_data <- u_obs[calibrated & c_ord <= j]
    upper_data <- u_obs[calibrated & c_ord > j]
    lower <- max(c(tau_bounds[1], lower_data), na.rm = TRUE)
    upper <- min(c(tau_bounds[2], upper_data), na.rm = TRUE)
    if (j > 1L) lower <- max(lower, tau[j - 1L] + 1e-6)
    if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
      stop("Calibrated latent values are incompatible with the ordinal order.")
    }
    tau[j] <- min(max(target[j], lower + 1e-5), upper - 1e-5)
  }
  tau
}

fit_threshold_measurement_response_free <- function(c_ord,
                                                    u_obs,
                                                    calib_idx,
                                                    m,
                                                    n_iter = 2500L,
                                                    burn = 750L,
                                                    thin = 1L,
                                                    seed = 1L,
                                                    tau_bounds = c(-8, 8)) {
  c_ord <- as.integer(c_ord)
  u_obs <- as.numeric(u_obs)
  calib_idx <- sort(unique(as.integer(calib_idx)))
  n <- length(c_ord)
  if (length(u_obs) != n || any(c_ord < 1L | c_ord > m)) {
    stop("Invalid Study I response-free measurement inputs.")
  }
  if (length(calib_idx) > 0L && any(!is.finite(u_obs[calib_idx]))) {
    stop("Calibrated u values must be finite.")
  }
  miss_idx <- setdiff(seq_len(n), calib_idx)
  tau <- initial_thresholds_response_free(
    c_ord, u_obs, calib_idx, m, tau_bounds
  )
  u <- numeric(n)
  if (length(calib_idx) > 0L) u[calib_idx] <- u_obs[calib_idx]
  if (length(miss_idx) > 0L) {
    lower <- c(-Inf, tau)[c_ord[miss_idx]]
    upper <- c(tau, Inf)[c_ord[miss_idx]]
    set.seed(seed)
    u[miss_idx] <- rtruncnorm_vec(0, 1, lower, upper)
  }

  keep_iter <- seq.int(burn + 1L, n_iter, by = thin)
  samples_u <- matrix(NA_real_, nrow = length(keep_iter), ncol = n)
  samples_tau <- matrix(NA_real_, nrow = length(keep_iter), ncol = m - 1L)
  set.seed(seed)
  keep_pos <- 0L
  for (iter in seq_len(n_iter)) {
    if (length(miss_idx) > 0L) {
      lower <- c(-Inf, tau)[c_ord[miss_idx]]
      upper <- c(tau, Inf)[c_ord[miss_idx]]
      u[miss_idx] <- rtruncnorm_vec(0, 1, lower, upper)
    }

    for (j in sample(seq_len(m - 1L))) {
      lower_neighbor <- if (j == 1L) tau_bounds[1] else tau[j - 1L]
      upper_neighbor <- if (j == m - 1L) tau_bounds[2] else tau[j + 1L]
      lower <- max(c(lower_neighbor, u[c_ord <= j]), na.rm = TRUE)
      upper <- min(c(upper_neighbor, u[c_ord > j]), na.rm = TRUE)
      if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
        stop("Threshold Gibbs update encountered an empty feasible interval.")
      }
      tau[j] <- stats::runif(1L, lower, upper)
    }

    if (iter > burn && ((iter - burn - 1L) %% thin == 0L)) {
      keep_pos <- keep_pos + 1L
      samples_u[keep_pos, ] <- u
      samples_tau[keep_pos, ] <- tau
    }
  }
  colnames(samples_u) <- paste0("u", seq_len(n))
  colnames(samples_tau) <- paste0("tau", seq_len(m - 1L))
  list(
    samples_u = samples_u,
    samples_tau = samples_tau,
    data = list(
      c_ord = c_ord, u_obs = u_obs, calib_idx = calib_idx,
      miss_idx = miss_idx, m = m, tau_bounds = tau_bounds
    ),
    control = list(n_iter = n_iter, burn = burn, thin = thin, seed = seed)
  )
}

truncated_standard_normal_mean <- function(lower, upper) {
  denom <- stats::pnorm(upper) - stats::pnorm(lower)
  numer <- stats::dnorm(lower) - stats::dnorm(upper)
  numer / pmax(denom, .Machine$double.eps)
}

posterior_mean_u_given_c_threshold <- function(measurement_fit, c_star) {
  c_star <- as.integer(c_star)
  tau <- measurement_fit$samples_tau
  m <- measurement_fit$data$m
  if (any(c_star < 1L | c_star > m)) stop("c_star is outside 1:m.")
  out <- numeric(length(c_star))
  for (j in seq_along(c_star)) {
    lower <- if (c_star[j] == 1L) -Inf else tau[, c_star[j] - 1L]
    upper <- if (c_star[j] == m) Inf else tau[, c_star[j]]
    out[j] <- mean(truncated_standard_normal_mean(lower, upper))
  }
  out
}

sample_u_given_c_threshold <- function(measurement_fit,
                                       c_star,
                                       n_draw,
                                       seed) {
  c_star <- as.integer(c_star)
  tau <- measurement_fit$samples_tau
  m <- measurement_fit$data$m
  if (any(c_star < 1L | c_star > m)) stop("c_star is outside 1:m.")
  set.seed(seed)
  draw_id <- sample(seq_len(nrow(tau)), n_draw, replace = TRUE)
  out <- matrix(NA_real_, nrow = n_draw, ncol = length(c_star))
  for (j in seq_along(c_star)) {
    lower <- if (c_star[j] == 1L) -Inf else tau[draw_id, c_star[j] - 1L]
    upper <- if (c_star[j] == m) Inf else tau[draw_id, c_star[j]]
    out[, j] <- rtruncnorm_vec(0, 1, lower, upper)
  }
  out
}

standardize_study1_gp_data <- function(x_raw, y_raw, row_id = seq_along(y_raw)) {
  x_raw <- as.matrix(x_raw)
  if (is.null(colnames(x_raw))) {
    colnames(x_raw) <- paste0("x", seq_len(ncol(x_raw)))
  }
  row_id <- as.integer(row_id)
  if (length(row_id) < 2L || any(row_id < 1L | row_id > nrow(x_raw))) {
    stop("row_id must select at least two valid training rows.")
  }
  x_center <- colMeans(x_raw[row_id, , drop = FALSE])
  x_scale <- apply(x_raw[row_id, , drop = FALSE], 2L, stats::sd)
  if (any(!is.finite(x_scale) | x_scale <= 0)) stop("Degenerate x scale.")
  y_center <- mean(y_raw[row_id])
  y_scale <- stats::sd(y_raw[row_id])
  if (!is.finite(y_scale) || y_scale <= 0) stop("Degenerate y scale.")
  list(
    x = sweep(sweep(x_raw, 2L, x_center, "-"), 2L, x_scale, "/"),
    y = (as.numeric(y_raw) - y_center) / y_scale,
    x_center = x_center, x_scale = x_scale,
    y_center = y_center, y_scale = y_scale
  )
}

fit_study1_pi_gp <- function(x_raw,
                             y_raw,
                             measurement_fit,
                             kernel = "se",
                             matern_nu = 2.5) {
  scaled <- standardize_study1_gp_data(x_raw, y_raw)
  u_hat <- colMeans(measurement_fit$samples_u)
  fit <- gp_mle_fit_1d(
    cbind(scaled$x, u = u_hat), scaled$y,
    kernel = kernel, matern_nu = matern_nu
  )
  list(gp = fit, measurement = measurement_fit, scale = scaled, u_hat = u_hat)
}

sample_study1_pi_gp <- function(fit,
                                x_test,
                                c_test,
                                n_draw,
                                seed) {
  x_test <- as.matrix(x_test)
  x_std <- sweep(
    sweep(x_test, 2L, fit$scale$x_center, "-"),
    2L, fit$scale$x_scale, "/"
  )
  colnames(x_std) <- colnames(fit$scale$x)
  u_hat <- posterior_mean_u_given_c_threshold(fit$measurement, c_test)
  set.seed(seed)
  draws_std <- sample_gp_mle_predictive_1d(
    fit$gp, cbind(x_std, u = u_hat), n_draw = n_draw
  )
  out <- fit$scale$y_center + fit$scale$y_scale * draws_std
  attr(out, "conditional_means") <-
    fit$scale$y_center + fit$scale$y_scale *
      attr(draws_std, "conditional_means")
  attr(out, "conditional_vars") <-
    fit$scale$y_scale^2 * attr(draws_std, "conditional_vars")
  attr(out, "mixture_components") <- "plug-in GP Gaussian predictor"
  out
}

fit_study1_cc_gp <- function(x_raw,
                             y_raw,
                             u_obs,
                             calib_idx,
                             measurement_fit,
                             kernel = "se",
                             matern_nu = 2.5) {
  calib_idx <- as.integer(calib_idx)
  if (length(calib_idx) < 3L) stop("CC-GP requires at least three complete cases.")
  scaled <- standardize_study1_gp_data(x_raw, y_raw, row_id = calib_idx)
  fit <- gp_mle_fit_1d(
    cbind(
      scaled$x[calib_idx, , drop = FALSE],
      u = u_obs[calib_idx]
    ),
    scaled$y[calib_idx],
    kernel = kernel, matern_nu = matern_nu
  )
  list(
    gp = fit, measurement = measurement_fit, scale = scaled,
    calib_idx = calib_idx
  )
}

sample_study1_cc_gp <- function(fit,
                                x_test,
                                c_test,
                                n_draw,
                                seed,
                                test_chunk_size = 50L) {
  x_test <- as.matrix(x_test)
  n_test <- nrow(x_test)
  x_std <- sweep(
    sweep(x_test, 2L, fit$scale$x_center, "-"),
    2L, fit$scale$x_scale, "/"
  )
  colnames(x_std) <- colnames(fit$scale$x)
  u_draws <- sample_u_given_c_threshold(
    fit$measurement, c_test, n_draw = n_draw, seed = seed
  )
  out <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  conditional_means <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  conditional_vars <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  chunks <- split(seq_len(n_test), ceiling(seq_len(n_test) / test_chunk_size))
  set.seed(seed + 1L)
  for (ids in chunks) {
    x_expanded <- x_std[rep(ids, each = n_draw), , drop = FALSE]
    u_expanded <- as.vector(u_draws[, ids, drop = FALSE])
    draws_std <- sample_gp_mle_predictive_1d(
      fit$gp,
      cbind(x_expanded, u = u_expanded),
      n_draw = 1L
    )
    out[, ids] <- matrix(
      as.numeric(draws_std), nrow = n_draw, ncol = length(ids)
    )
    conditional_means[, ids] <- matrix(
      as.numeric(attr(draws_std, "conditional_means")),
      nrow = n_draw, ncol = length(ids)
    )
    conditional_vars[, ids] <- matrix(
      as.numeric(attr(draws_std, "conditional_vars")),
      nrow = n_draw, ncol = length(ids)
    )
  }
  out <- fit$scale$y_center + fit$scale$y_scale * out
  attr(out, "conditional_means") <-
    fit$scale$y_center + fit$scale$y_scale * conditional_means
  attr(out, "conditional_vars") <-
    fit$scale$y_scale^2 * conditional_vars
  attr(out, "mixture_components") <-
    "measurement and complete-case GP Monte Carlo Gaussian components"
  out
}

############################################################
## Repeated-simulation evaluators for the latent targets
############################################################

study1_latent_imputation_status_rows <- function(rep_id,
                                                  n_calib,
                                                  scenario,
                                                  methods,
                                                  n_missing,
                                                  identification_status,
                                                  score_status,
                                                  reason) {
  data.frame(
    rep = as.integer(rep_id),
    scenario = as.character(scenario),
    n_calib = as.integer(n_calib),
    target = "training_latent_u",
    method = as.character(methods),
    n_missing = as.integer(n_missing),
    identification_status = as.character(identification_status),
    score_status = as.character(score_status),
    reason = as.character(reason),
    Bias = NA_real_,
    RMSE = NA_real_,
    MAE = NA_real_,
    Coverage95 = NA_real_,
    Width95 = NA_real_,
    stringsAsFactors = FALSE
  )
}

summarize_study1_training_u_draws <- function(draws,
                                               truth,
                                               method,
                                               rep_id,
                                               n_calib,
                                               scenario) {
  draws <- as.matrix(draws)
  truth <- as.numeric(truth)
  if (nrow(draws) < 2L || ncol(draws) != length(truth) ||
      length(truth) < 1L || any(!is.finite(draws)) ||
      any(!is.finite(truth))) {
    stop("Latent-U draws and truth are incompatible or nonfinite.")
  }

  post_mean <- colMeans(draws)
  lower <- apply(draws, 2L, stats::quantile, probs = 0.025, names = FALSE)
  upper <- apply(draws, 2L, stats::quantile, probs = 0.975, names = FALSE)
  error <- post_mean - truth

  data.frame(
    rep = as.integer(rep_id),
    scenario = as.character(scenario),
    n_calib = as.integer(n_calib),
    target = "training_latent_u",
    method = as.character(method),
    n_missing = length(truth),
    identification_status = "anchored_raw_scale",
    score_status = "scored",
    reason = "",
    Bias = mean(error),
    RMSE = sqrt(mean(error^2)),
    MAE = mean(abs(error)),
    Coverage95 = mean(truth >= lower & truth <= upper),
    Width95 = mean(upper - lower),
    stringsAsFactors = FALSE
  )
}

study1_latent_imputation_metrics <- function(eiv_fit,
                                              measurement_fit,
                                              U_true,
                                              rep_id,
                                              n_calib,
                                              scenario) {
  methods <- c("EIV-GP", "Response-free threshold model")
  U_true <- as.numeric(U_true)
  n_calib <- as.integer(n_calib)
  if (length(n_calib) != 1L || is.na(n_calib) || n_calib < 0L ||
      length(U_true) < 1L || any(!is.finite(U_true))) {
    stop("Invalid Study I latent-imputation inputs.")
  }

  ## Without calibration, the likelihood identifies an ordinally compatible
  ## working coordinate, not the physical scale of U.  Reporting raw-scale
  ## RMSE or coverage in that case would silently use simulation truth to align
  ## the fitted coordinate, so the case is retained but deliberately unscored.
  if (n_calib == 0L) {
    return(study1_latent_imputation_status_rows(
      rep_id = rep_id,
      n_calib = n_calib,
      scenario = scenario,
      methods = methods,
      n_missing = length(U_true),
      identification_status = "unanchored_working_scale",
      score_status = "not_scored",
      reason = paste(
        "No calibration measurements: the physical location and scale of U",
        "are not anchored."
      )
    ))
  }

  if (is.null(eiv_fit$mcmc$samples_u) ||
      is.null(measurement_fit$samples_u) ||
      is.null(eiv_fit$data$miss_idx) ||
      is.null(measurement_fit$data$miss_idx)) {
    stop("Both fitted objects must retain posterior training-U draws.")
  }
  if (!isTRUE(eiv_fit$data$latent_scale_anchored)) {
    return(study1_latent_imputation_status_rows(
      rep_id = rep_id,
      n_calib = n_calib,
      scenario = scenario,
      methods = methods,
      n_missing = length(eiv_fit$data$miss_idx),
      identification_status = "insufficient_calibration_rank",
      score_status = "not_scored",
      reason = paste(
        "Calibration values do not have the affine rank required to anchor",
        "the physical scale of U."
      )
    ))
  }

  miss_idx <- as.integer(eiv_fit$data$miss_idx)
  measurement_miss_idx <- as.integer(measurement_fit$data$miss_idx)
  if (!identical(miss_idx, measurement_miss_idx)) {
    stop("EIV-GP and response-free fits must use the same calibration rows.")
  }
  if (ncol(eiv_fit$mcmc$samples_u) != length(U_true) ||
      ncol(measurement_fit$samples_u) != length(U_true)) {
    stop("Posterior training-U draws do not match the supplied truth.")
  }
  if (length(miss_idx) == 0L) {
    return(study1_latent_imputation_status_rows(
      rep_id = rep_id,
      n_calib = n_calib,
      scenario = scenario,
      methods = methods,
      n_missing = 0L,
      identification_status = "anchored_raw_scale",
      score_status = "not_scored",
      reason = "Every training U is observed; there are no latent values to impute."
    ))
  }

  rbind(
    summarize_study1_training_u_draws(
      draws = eiv_fit$mcmc$samples_u[, miss_idx, drop = FALSE],
      truth = U_true[miss_idx],
      method = methods[1L],
      rep_id = rep_id,
      n_calib = n_calib,
      scenario = scenario
    ),
    summarize_study1_training_u_draws(
      draws = measurement_fit$samples_u[, miss_idx, drop = FALSE],
      truth = U_true[miss_idx],
      method = methods[2L],
      rep_id = rep_id,
      n_calib = n_calib,
      scenario = scenario
    )
  )
}

study1_surface_grid <- function(grid_n_x = 31L,
                                grid_n_u = 51L,
                                x_lim = c(-2, 2),
                                u_lim = c(-2.5, 2.5)) {
  grid_n_x <- as.integer(grid_n_x)
  grid_n_u <- as.integer(grid_n_u)
  if (anyNA(c(grid_n_x, grid_n_u)) || grid_n_x < 2L || grid_n_u < 2L ||
      length(x_lim) != 2L || length(u_lim) != 2L ||
      any(!is.finite(c(x_lim, u_lim))) || x_lim[1L] >= x_lim[2L] ||
      u_lim[1L] >= u_lim[2L]) {
    stop("Invalid Study I surface-grid specification.")
  }
  expand.grid(
    x = seq(x_lim[1L], x_lim[2L], length.out = grid_n_x),
    u = seq(u_lim[1L], u_lim[2L], length.out = grid_n_u),
    KEEP.OUT.ATTRS = FALSE
  )
}

study1_surface_truth <- function(grid,
                                 scenario,
                                 tau_true,
                                 heterogeneity_eta = 1) {
  f0_1d(
    x = grid$x,
    u = grid$u,
    scenario = scenario,
    c_ord = make_class(grid$u, tau_true),
    tau = tau_true,
    heterogeneity_eta = heterogeneity_eta
  )
}

study1_surface_metric_row <- function(post_mean,
                                      lower,
                                      upper,
                                      truth,
                                      method,
                                      rep_id,
                                      n_calib,
                                      scenario,
                                      grid_n_x,
                                      grid_n_u) {
  if (!all(lengths(list(post_mean, lower, upper, truth)) == length(truth)) ||
      length(truth) < 1L ||
      any(!is.finite(c(post_mean, lower, upper, truth)))) {
    stop("Invalid values supplied for Study I surface scoring.")
  }
  data.frame(
    rep = as.integer(rep_id),
    scenario = as.character(scenario),
    n_calib = as.integer(n_calib),
    target = "latent_surface_f(x,u)",
    method = as.character(method),
    grid_n_x = as.integer(grid_n_x),
    grid_n_u = as.integer(grid_n_u),
    grid_size = length(truth),
    ISE = mean((post_mean - truth)^2),
    Bias = mean(post_mean - truth),
    Coverage95 = mean(truth >= lower & truth <= upper),
    Width95 = mean(upper - lower),
    stringsAsFactors = FALSE
  )
}

study1_eiv_surface_recovery_metrics <- function(
    fit,
    scenario,
    rep_id,
    n_calib,
    tau_true,
    heterogeneity_eta = 1,
    grid_n_x = 31L,
    grid_n_u = 51L,
    max_draw = 200L,
    x_lim = c(-2, 2),
    u_lim = c(-2.5, 2.5),
    seed = NULL) {
  n_calib <- as.integer(n_calib)
  if (length(n_calib) != 1L || is.na(n_calib) || n_calib < 0L) {
    stop("n_calib must be one nonnegative integer.")
  }
  if (n_calib == 0L) return(data.frame())
  if (!isTRUE(fit$data$latent_scale_anchored)) {
    stop("The physical U scale must be anchored before scoring f(x,u).")
  }
  if (!is.null(seed)) set.seed(as.integer(seed))

  grid <- study1_surface_grid(grid_n_x, grid_n_u, x_lim, u_lim)
  n_saved <- nrow(fit$mcmc$samples_u)
  max_draw <- as.integer(max_draw)
  if (is.na(max_draw) || max_draw < 2L || n_saved < 2L) {
    stop("At least two saved posterior draws are required for surface scoring.")
  }
  draw_ids <- if (n_saved <= max_draw) {
    seq_len(n_saved)
  } else {
    unique(round(seq(1, n_saved, length.out = max_draw)))
  }

  conditional_means <- sample_eiv_f_given_xu_1d(
    x_star_raw = matrix(grid$x, ncol = 1L),
    u_star = grid$u,
    fit_obj = fit,
    draw_ids = draw_ids,
    include_gp_uncertainty = FALSE,
    u_input_scale = "raw",
    joint = FALSE
  )
  f_draws <- sample_eiv_f_given_xu_1d(
    x_star_raw = matrix(grid$x, ncol = 1L),
    u_star = grid$u,
    fit_obj = fit,
    draw_ids = draw_ids,
    include_gp_uncertainty = TRUE,
    u_input_scale = "raw",
    joint = FALSE
  )
  truth <- study1_surface_truth(
    grid, scenario, tau_true, heterogeneity_eta
  )

  study1_surface_metric_row(
    post_mean = colMeans(conditional_means),
    lower = apply(f_draws, 2L, stats::quantile, probs = 0.025, names = FALSE),
    upper = apply(f_draws, 2L, stats::quantile, probs = 0.975, names = FALSE),
    truth = truth,
    method = "EIV-GP",
    rep_id = rep_id,
    n_calib = n_calib,
    scenario = scenario,
    grid_n_x = grid_n_x,
    grid_n_u = grid_n_u
  )
}

fit_study1_full_u_gp <- function(x_raw,
                                 y_raw,
                                 u_true,
                                 kernel = "se",
                                 matern_nu = 2.5) {
  u_true <- as.numeric(u_true)
  if (length(u_true) != length(y_raw) || any(!is.finite(u_true))) {
    stop("u_true must contain one finite value per training response.")
  }
  scaled <- standardize_study1_gp_data(x_raw, y_raw)
  fit <- gp_mle_fit_1d(
    cbind(scaled$x, u = u_true),
    scaled$y,
    kernel = kernel,
    matern_nu = matern_nu
  )
  list(
    gp = fit,
    scale = scaled,
    u_train = u_true,
    train_idx = seq_along(y_raw)
  )
}

study1_latent_gp_surface_metrics <- function(
    gp_fit,
    scenario,
    method,
    rep_id,
    n_calib = NA_integer_,
    tau_true,
    heterogeneity_eta = 1,
    grid_n_x = 31L,
    grid_n_u = 51L,
    x_lim = c(-2, 2),
    u_lim = c(-2.5, 2.5)) {
  if (!is.na(n_calib) && n_calib <= 0L) return(data.frame())
  if (is.null(gp_fit$gp) || is.null(gp_fit$scale)) {
    stop("gp_fit must be a Study I PI-GP, CC-GP, or Full-U GP fit.")
  }

  grid <- study1_surface_grid(grid_n_x, grid_n_u, x_lim, u_lim)
  x_std <- sweep(
    sweep(
      matrix(grid$x, ncol = 1L),
      2L,
      gp_fit$scale$x_center,
      "-"
    ),
    2L,
    gp_fit$scale$x_scale,
    "/"
  )
  colnames(x_std) <- colnames(gp_fit$scale$x)
  pred <- gp_mle_predict_1d(
    gp_fit$gp,
    Xstar = cbind(x_std, u = grid$u),
    noisy = FALSE
  )
  post_mean <- gp_fit$scale$y_center + gp_fit$scale$y_scale * pred$mean
  post_sd <- gp_fit$scale$y_scale * sqrt(pred$var)
  truth <- study1_surface_truth(
    grid, scenario, tau_true, heterogeneity_eta
  )
  cutoff <- stats::qnorm(0.975)

  study1_surface_metric_row(
    post_mean = post_mean,
    lower = post_mean - cutoff * post_sd,
    upper = post_mean + cutoff * post_sd,
    truth = truth,
    method = method,
    rep_id = rep_id,
    n_calib = n_calib,
    scenario = scenario,
    grid_n_x = grid_n_x,
    grid_n_u = grid_n_u
  )
}
