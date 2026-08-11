# Generated from _main.Rmd: do not edit by hand

#' Extract imputation summary for missing latent variables
#'
#' For an EIV-GP fit object and a given calibration size, computes
#' posterior summaries (mean, 2.5%, 97.5% quantiles) of the latent
#' variables `u` for observations that were *not* part of the
#' calibration set.  The true (simulated) `u` values are also included.
#'
#' @param fit A fitted EIV-GP model object, as returned by
#'   `fit_eivgp_1d`.
#' @param n_calib Integer, the calibration set size used for this fit
#'   (stored in the output as a metadata column).
#'
#' @return A data frame with one row per missing observation and
#'   columns: `n_calib`, `id` (original index), `true_u`, `post_mean`,
#'   `post_lo`, `post_hi`, `c` (ordinal class factor), `covered`
#'   (whether the 95% posterior interval covers `true_u`).
#' @export
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

#' Construct EIV-GP function-slice data for plotting
#'
#' Given an EIV-GP fit, computes posterior predictive summaries of
#' the latent regression function `f(x, u)` along several fixed-`x`
#' slices over a dense grid of `u`.  Uses a subset of posterior draws
#' (controlled by the global variable `n_pred_draw` when `draw_ids`
#' is not supplied).
#'
#' @param fit EIV-GP fit object.
#' @param x_slices Numeric vector of `x` values where slices are taken
#'   (default `c(-1, 0, 1)` on the original scale).
#' @param u_grid Numeric vector of `u` values along which the
#'   function is evaluated.
#' @param draw_ids Optional integer vector of posterior draw indices
#'   to use.  If `NULL`, all saved draws are considered.
#' @param scenario Character, `"active"` or `"inactive"`, used to
#'   compute the true function for comparison.
#' @param label Character string used as the method label in the output.
#'
#' @return A data frame with columns `x_raw`, `u`, `mean` (posterior
#'   mean of `f`), `lo` (2.5% quantile), `hi` (97.5% quantile),
#'   `truth` (true `f` value), `method`, and `x_slice` (factor).
#'
#' @note This function expects the global variable `n_pred_draw` to be
#'   defined when `draw_ids` is not provided.
#' @export
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

#' Construct complete-case GP function-slice data
#'
#' For an EIV-GP fit that used calibration observations, this fits a
#' standard GP *only* to the calibrated data (using the true `u`
#' values) and computes predictive summaries for the latent function
#' along the same `x`-slices.  If fewer than 3 calibrated observations
#' are available, an empty data frame is returned.
#'
#' @param fit EIV-GP fit object (must contain `data$calib_idx` with
#'   at least 3 observations for a meaningful GP fit).
#' @param x_slices Numeric vector of `x` values (default `c(-1,0,1)`).
#' @param u_grid Numeric vector of `u` points for evaluation.
#' @param scenario Character, `"active"` or `"inactive"`.
#' @param label Character label for the method (default `"Complete-case GP"`).
#'
#' @return A data frame with columns `x_raw`, `u`, `mean`, `lo`, `hi`,
#'   `truth`, `method`, and `x_slice`, or an empty data frame if
#'   calibration data are insufficient.
#' @export
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

#' Reshape predictive draws to long format
#'
#' Converts a matrix of predictive draws (rows = draws, columns =
#' test points) into a long-format data frame suitable for plotting
#' with `ggplot2`.  Each row corresponds to one draw for one test
#' point, annotated with the corresponding covariates from `grid`.
#'
#' @param draw_mat Numeric matrix of draws (`n_draw` rows by
#'   `nrow(grid)` columns).
#' @param method Character string identifying the method (stored in
#'   the output column `method`).
#' @param grid A data frame with at least columns `x_star` and
#'   `c_star`, having one row per test point.
#'
#' @return A data frame with columns `x_star`, `c_star`, `y` (draw),
#'   and `method`.
#' @export
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
