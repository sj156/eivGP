# Generated from _main.Rmd: do not edit by hand

#' Null-coalescing operator
#'
#' Returns `x` unless it is `NULL`, in which case `y` is returned.
#'
#' @param x The preferred value.
#' @param y The fallback value.
#'
#' @return `x` if not `NULL`, otherwise `y`.
#' @export
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Save a ggplot object as PDF and/or PNG
#'
#' Saves a ggplot to one or both formats depending on the global
#' logical flags `STUDY1_SAVE_PDF` and `STUDY1_SAVE_PNG`.  The
#' output directory is created if necessary.
#'
#' @param path_no_ext File path without extension.
#' @param plot A ggplot object.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @param dpi Resolution (dots per inch) for PNG output.
#' @param save_pdf,save_plot BOOL VALUEs to save document
#'
#' @return The plot, invisibly.
#' @export
save_plot <- function(path_no_ext,
                      plot,
                      width,
                      height,
                      dpi = 320,
                      save_pdf = TRUE,
                      save_png = TRUE) {
  dir.create(dirname(path_no_ext), showWarnings = FALSE, recursive = TRUE)
  
  if (isTRUE(save_pdf)) {
    ggsave(
      filename = paste0(path_no_ext, ".pdf"),
      plot = plot, width = width, height = height,
      units = "in", limitsize = FALSE
    )
  }
  
  if (isTRUE(save_png)) {
    ggsave(
      filename = paste0(path_no_ext, ".png"),
      plot = plot, width = width, height = height,
      units = "in", dpi = dpi, bg = "white", limitsize = FALSE
    )
  }
  
  invisible(plot)
}

#' Safe plot saving with error handling
#'
#' Calls `save_plot` and catches any errors, issuing a warning instead
#' of stopping execution.
#'
#' @param path_no_ext File path without extension.
#' @param plot A ggplot object.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @param dpi Resolution for PNG output.
#' @param save_pdf,save_plot BOOL VALUEs to save document
#'
#' @return The plot invisibly on success; `NULL` invisibly on failure.
#' @export
safe_save_plot <- function(path_no_ext,
                           plot,
                           width,
                           height,
                           dpi = 320,
                           save_pdf = TRUE,
                           save_png = TRUE) {
  tryCatch(
    save_plot(path_no_ext, plot, width, height, dpi, save_pdf, save_png),
    error = function(e) {
      warning("Could not save plot ", basename(path_no_ext), ": ",
              conditionMessage(e), call. = FALSE)
      invisible(NULL)
    }
  )
}

#' Write a data frame to CSV, creating directories as needed
#'
#' @param x A data frame or matrix.
#' @param path File path for the CSV output.
#'
#' @return The input `x`, invisibly.
#' @export
save_csv <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(x, path, row.names = FALSE)
  invisible(x)
}

#' Select a subset of IDs by quantiles of a numeric vector
#'
#' Returns approximately `n` IDs whose associated values are evenly
#' spaced in sorted order, useful for choosing representative
#' observations.
#'
#' @param ids Integer vector of identifiers.
#' @param values Numeric vector of the same length as `ids`, used for
#'   sorting.
#' @param n Integer, desired number of IDs to return (default 12).
#'
#' @return A subset of `ids` (possibly smaller if fewer are available).
#' @export
select_ids_by_quantile <- function(ids, values, n = 12L) {
  if (length(ids) == 0L) return(integer(0))
  keep <- is.finite(values)
  ids <- ids[keep]
  values <- values[keep]
  if (length(ids) == 0L) return(integer(0))
  
  ord_ids <- ids[order(values)]
  ord_ids[unique(round(seq(1, length(ord_ids), length.out = min(n, length(ord_ids)))))]
}

#' Convert scaled x back to original scale using fit centering
#'
#' Given an EIV-GP fit object and (optionally) its internal scaled
#' `x` values, returns the original-scale `x` values.
#'
#' @param fit A fitted EIV-GP object containing `x_center` and
#'   `x_scale` elements.
#' @param x Numeric vector of scaled x values (defaults to
#'   `fit$data$x`).
#'
#' @return A numeric vector on the original measurement scale.
#' @export
raw_x_from_fit <- function(fit, x = fit$data$x) {
  center <- fit$data$x_center %||% 0
  scale <- fit$data$x_scale %||% 1
  as.numeric(center + scale * x)
}

#' Convert scaled y back to original scale using fit centering
#'
#' Similar to `raw_x_from_fit`, but for the response variable.
#'
#' @param fit A fitted EIV-GP object.
#' @param y Numeric vector of scaled y values (defaults to
#'   `fit$data$y`).
#'
#' @return A numeric vector on the original measurement scale.
#' @export
raw_y_from_fit <- function(fit, y = fit$data$y) {
  center <- fit$data$y_center %||% 0
  scale <- fit$data$y_scale %||% 1
  as.numeric(center + scale * y)
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
  draw_mat <- as.matrix(draw_mat)
  stopifnot(ncol(draw_mat) == nrow(grid))
  
  data.frame(
    x_star = rep(grid$x_star, each = nrow(draw_mat)),
    c_star = rep(grid$c_star, each = nrow(draw_mat)),
    y = as.vector(draw_mat),
    method = method
  )
}

#' Summarise a matrix of predictive draws with quantiles and moments
#'
#' For each column (test point) of a draw matrix, computes the mean,
#' standard deviation, and several quantiles (2.5%, 5%, 10%, 25%, 50%,
#' 75%, 90%, 95%, 97.5%).  Also adds interval widths for 50%, 80%,
#' 90%, and 95% intervals.  If the matrix appears transposed (rows
#' match the grid size), it is automatically transposed.
#'
#' @param draw_mat Numeric matrix of draws, expected to have columns
#'   corresponding to test points.
#' @param method Character string identifying the method.
#' @param grid A data frame with one row per test point; its columns
#'   are included in the output.
#'
#' @return A data frame with the same number of rows as `grid`,
#'   augmented with summary statistics.
#' @export
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

#' Extract imputation summary data frame from fitted EIV-GP model
#'
#' @param fit A fitted model object containing MCMC samples and data.
#' @param n_calib Integer, the calibration sample size.
#' @param n_calib_label Character, a label for the calibration size. If NULL, defaults to as.character(n_calib).
#' @param m Integer, number of latent classes. Defaults to config$m or 6L.
#' @param config A list of configuration parameters. Used to provide default values for m and other settings.
#'        If not supplied, the function attempts to read from getOption("study1.config").
#' @return A data.frame with columns: n_calib, n_calib_label, id, x, y, true_u, post_mean, post_lo, post_hi,
#'         post_sd, width, error, abs_error, c (factor), covered.
#' @export
#'
#' @examples
#' \dontrun{
#'   fit <- fit_eivgp_1d(train, test, ...)
#'   df <- extract_imputation_df(fit, n_calib = 10, m = 6)
#' }
extract_imputation_df <- function(
    fit,
    n_calib,
    n_calib_label = NULL,
    m = NULL,
    config = getOption("study1.config", default = list())
) {
  # Helper operator if not defined
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  # Get m from config or default
  m <- m %||% config$m %||% 6L
  
  # Default label
  if (is.null(n_calib_label)) {
    n_calib_label <- as.character(n_calib)
  }
  
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
    n_calib_label = n_calib_label,
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

#' Extract posterior density draws for latent variables
#'
#' This function extracts a subset of MCMC posterior draws for the latent
#' variables corresponding to a specified set of observation indices.
#' It returns a tidy data frame suitable for plotting density distributions
#' or comparing posterior draws to the true latent values.
#'
#' @param fit A fitted model object. It must contain:
#'   \itemize{
#'     \item \code{fit$mcmc$samples_u}: a matrix of posterior samples
#'           (rows = MCMC draws, columns = observations).
#'     \item \code{fit$data$c_ord}: integer vector of class assignments
#'           for all observations.
#'     \item \code{fit$data$u_true}: true latent values for all observations.
#'   }
#' @param ids Integer vector. Observation indices (1‑based) for which to
#'   extract draws. Indices out of range are silently ignored.
#' @param n_calib Integer (optional). The calibration sample size. If not
#'   provided, the column will contain \code{NA_integer_}.
#' @param max_draw Integer. Maximum number of MCMC draws to sample.
#'   If the total number of draws exceeds this, a random subset is taken.
#'   Defaults to \code{1000L}.
#' @param m Integer. Number of latent classes. Used to set factor levels
#'   for the class assignment. Defaults to \code{config$m} or \code{6L}
#'   if not supplied
#' @return A \code{data.frame} with one row per draw per selected observation,
#'   containing:
#'   \item{draw}{the MCMC iteration index (after thinning / subsetting)}
#'   \item{id}{the observation index}
#'   \item{id_label}{a descriptive label: \code{"id=<id>, c=<class>"}}
#'   \item{u_draw}{the posterior draw value for the latent variable}
#'   \item{true_u}{the true latent value}
#'   \item{c}{factor of class assignment, with levels \code{1:m}}
#'   \item{n_calib}{the calibration sample size (or \code{NA})}
#' @export
make_u_density_df <- function(
    fit,
    ids,
    n_calib = NA_integer_,
    max_draw = 1000L,
    m = 6L
) {
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
        n_calib = n_calib
      )
    })
  )
}


#' Retrieve the tau matrix from a fitted model
#'
#' This internal helper function searches for the MCMC samples corresponding
#' to the tau parameter (or its equivalent) in the fitted model object.
#'
#' @param fit A fitted model object containing an \code{mcmc} element.
#' @return A numeric matrix of MCMC samples, or \code{NULL} if none of the
#'   candidate names are found.
#' @keywords internal
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
  obj
}

#' Extract tau posterior samples into a tidy data frame
#'
#' This function extracts the posterior MCMC samples for the tau (or
#' cutpoint) parameters and returns them in a long-format data frame,
#' suitable for plotting or further analysis.
#'
#' @param fit A fitted model object from which the tau matrix can be
#'   extracted via \code{get_tau_matrix(fit)}.
#' @param n_calib Integer. The calibration sample size. Default is
#'   \code{NA_integer_}.
#' @param n_calib_label Character. A human-readable label for the calibration
#'   size. If \code{NULL} (default), it is set to \code{as.character(n_calib)}.
#' @param n_tau Integer. The number of tau components (used to set factor
#'   levels). If \code{NULL} (default), it is derived from
#'   \code{ncol(get_tau_matrix(fit))}.
#' @return A \code{data.frame} with columns:
#'   \item{n_calib}{calibration sample size}
#'   \item{n_calib_label}{label for the calibration size}
#'   \item{draw}{MCMC iteration index}
#'   \item{tau_index}{factor indicating the tau component (1 to \code{n_tau})}
#'   \item{tau}{the posterior sample value}
#' @export
#'
#' @examples
#' \dontrun{
#'   df <- extract_tau_df(fit, n_calib = 50, n_tau = 6)
#' }
extract_tau_df <- function(
    fit,
    n_calib = NA_integer_,
    n_calib_label = NULL,
    n_tau = NULL
) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  mat <- get_tau_matrix(fit)
  if (is.null(mat)) return(data.frame())
  
  n_tau <- n_tau %||% ncol(mat)
  n_calib_label <- n_calib_label %||% as.character(n_calib)
  
  dplyr::bind_rows(
    lapply(seq_len(ncol(mat)), function(jj) {
      data.frame(
        n_calib = n_calib,
        n_calib_label = n_calib_label,
        draw = seq_len(nrow(mat)),
        tau_index = factor(jj, levels = seq_len(n_tau)),
        tau = mat[, jj]
      )
    })
  )
}


#' Extract slice data for EIV-GP function surface
#'
#' This function computes posterior mean and credible intervals for the
#' EIV-GP predictive function over a grid of latent `u` values at fixed
#' predictor `x` slices. It can optionally include process uncertainty
#' via Monte Carlo sampling.
#'
#' @param fit A fitted EIV-GP model object. Must contain:
#'   \itemize{
#'     \item \code{fit$mcmc$samples_u}: latent variable posterior draws.
#'     \item \code{fit$mcmc$samples_logtheta}: log(theta) posterior draws.
#'     \item \code{fit$mcmc$samples_sigma2}: sigma^2 posterior draws.
#'     \item \code{fit$data$x, $y, $u_true}: training data.
#'     \item \code{fit$data$x_center, $x_scale, $y_center, $y_scale}: scaling parameters.
#'   }
#' @param x_slices Numeric vector. Values of the observed predictor `x`
#'   at which to evaluate the function. Default is \code{c(-1, 0, 1)}.
#' @param u_grid Numeric vector. Grid of latent `u` values to evaluate.
#'   Default is \code{seq(-2.4, 2.4, length.out = 160)}.
#' @param draw_ids Integer vector (optional). Indices of MCMC draws to use.
#'   If \code{NULL}, all draws are used. Default is \code{NULL}.
#' @param scenario Character. Scenario name passed to the true function
#'   \code{f0_fn}. Default is \code{"active"}.
#' @param label Character. Method label for the output data frame.
#'   Default is \code{"EIV-GP"}.
#' @param max_draw Integer or \code{NULL}. Maximum number of draws to retain.
#'   If the number of selected draws exceeds this value, a random subset is taken.
#'   Default is \code{1000L}.
#' @param include_process_uncertainty Logical. If \code{TRUE}, adds draw‑level
#'   Gaussian noise from the predictive variance; otherwise uses only the mean.
#'   Default is \code{TRUE}.
#' @param f0_fn Function. True data‑generating function with signature
#'   \code{function(x, u, scenario)}. Default is \code{f0_1d} (assumed to
#'   be available in the calling environment).
#' @return A \code{data.frame} with columns:
#'   \item{x_raw}{original predictor value}
#'   \item{u}{latent grid value}
#'   \item{mean}{posterior mean of the function}
#'   \item{lo}{2.5% posterior quantile}
#'   \item{hi}{97.5% posterior quantile}
#'   \item{width}{= \code{hi - lo}}
#'   \item{truth}{true function value from \code{f0_fn}}
#'   \item{error}{= \code{mean - truth}}
#'   \item{method}{method label}
#'   \item{x_slice}{factor of \code{x_raw} with labels \code{"x = ..."}}
#' @export
#'
#' @examples
#' \dontrun{
#'   df <- make_eiv_function_slice_df(fit, x_slices = c(-1, 0, 1),
#'                                    max_draw = 500L, f0_fn = f0_1d)
#' }
make_eiv_function_slice_df <- function(
    fit,
    x_slices = c(-1, 0, 1),
    u_grid = seq(-2.4, 2.4, length.out = 160),
    draw_ids = NULL,
    scenario = "active",
    label = "EIV-GP",
    max_draw = 1000L,
    include_process_uncertainty = TRUE,
    f0_fn = f0_1d
) {
  # If draw_ids not given, use all available draws
  if (is.null(draw_ids)) {
    draw_ids <- seq_len(nrow(fit$mcmc$samples_u))
  }
  
  # Subsample draws if needed
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
  out$truth <- f0_fn(out$x_raw, out$u, scenario = scenario)
  out$error <- out$mean - out$truth
  out$method <- label
  out$x_slice <- factor(paste0("x = ", out$x_raw), levels = paste0("x = ", x_slices))
  
  out
}


#' Extract slice data for Complete‑case GP function surface
#'
#' This function fits a standard GP to the complete‑case data (calibration
#' subset) and evaluates the predictive mean and 95% pointwise intervals
#' over a grid of latent `u` values at fixed `x` slices.
#'
#' @param fit A fitted model object. Must contain:
#'   \itemize{
#'     \item \code{fit$data$calib_idx}: indices of complete‑case observations.
#'     \item \code{fit$data$x, $u_true, $y}: full data vectors.
#'     \item \code{fit$data$x_center, $x_scale, $y_center, $y_scale}: scaling parameters.
#'   }
#' @param x_slices Numeric vector. Values of `x` at which to evaluate.
#'   Default is \code{c(-1, 0, 1)}.
#' @param u_grid Numeric vector. Grid of latent `u` values.
#'   Default is \code{seq(-2.4, 2.4, length.out = 160)}.
#' @param scenario Character. Scenario name passed to \code{f0_fn}.
#'   Default is \code{"active"}.
#' @param label Character. Method label for the output data frame.
#'   Default is \code{"Complete-case GP"}.
#' @param f0_fn Function. True data‑generating function with signature
#'   \code{function(x, u, scenario)}. Default is \code{f0_1d} (assumed to
#'   be available in the calling environment).
#' @return A \code{data.frame} with the same structure as
#'   \code{make_eiv_function_slice_df}. If fewer than 3 calibration
#'   observations are available, an empty data frame is returned.
#' @export
#'
#' @examples
#' \dontrun{
#'   df <- make_cc_function_slice_df(fit, x_slices = c(-1, 0, 1), f0_fn = f0_1d)
#' }
make_cc_function_slice_df <- function(
    fit,
    x_slices = c(-1, 0, 1),
    u_grid = seq(-2.4, 2.4, length.out = 160),
    scenario = "active",
    label = "Complete-case GP",
    f0_fn = f0_1d
) {
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
  out$truth <- f0_fn(out$x_raw, out$u, scenario = scenario)
  out$error <- out$mean - out$truth
  out$method <- label
  out$x_slice <- factor(paste0("x = ", out$x_raw), levels = paste0("x = ", x_slices))
  
  out
}

#' Extract EIV-GP function surface over a 2D grid
#'
#' This function computes the posterior mean and pointwise credible intervals
#' of the EIV-GP predictive function over a dense grid of (x, u) values.
#' It can optionally include process uncertainty via Monte Carlo sampling.
#'
#' @param fit A fitted EIV-GP model object. Must contain:
#'   \itemize{
#'     \item \code{fit$mcmc$samples_u}: latent variable posterior draws.
#'     \item \code{fit$mcmc$samples_logtheta}: log(theta) posterior draws.
#'     \item \code{fit$mcmc$samples_sigma2}: sigma^2 posterior draws.
#'     \item \code{fit$data$x, $y}: training data vectors.
#'     \item \code{fit$data$x_center, $x_scale, $y_center, $y_scale}: scaling parameters.
#'   }
#' @param x_grid Numeric vector. Grid of observed predictor `x` values.
#'   Default is \code{seq(-1.75, 1.75, length.out = 60)}.
#' @param u_grid Numeric vector. Grid of latent `u` values.
#'   Default is \code{seq(-2.6, 2.6, length.out = 75)}.
#' @param draw_ids Integer vector (optional). Indices of MCMC draws to use.
#'   If \code{NULL}, all available draws are used. Default is \code{NULL}.
#' @param scenario Character. Scenario name passed to the true function
#'   \code{f0_fn}. Default is \code{"active"}.
#' @param max_draw Integer or \code{NULL}. Maximum number of draws to retain.
#'   If the number of selected draws exceeds this value, a random subset is taken.
#'   Default is \code{500L}.
#' @param include_process_uncertainty Logical. If \code{TRUE}, adds draw‑level
#'   Gaussian noise from the predictive variance; otherwise uses only the mean.
#'   Default is \code{FALSE}.
#' @param f0_fn Function. True data‑generating function with signature
#'   \code{function(x, u, scenario)}. Default is \code{f0_1d} (assumed to
#'   be available in the calling environment).
#' @return A \code{data.frame} with columns:
#'   \item{x_raw}{original predictor value}
#'   \item{u}{latent grid value}
#'   \item{mean}{posterior mean of the function}
#'   \item{lo}{2.5% posterior quantile}
#'   \item{hi}{97.5% posterior quantile}
#'   \item{width}{= \code{hi - lo}}
#'   \item{truth}{true function value from \code{f0_fn}}
#'   \item{error}{= \code{mean - truth}}
#' @export
#'
#' @examples
#' \dontrun{
#'   df <- make_eiv_function_surface_df(fit, max_draw = 300L, f0_fn = f0_1d)
#' }
make_eiv_function_surface_df <- function(
    fit,
    x_grid = seq(-1.75, 1.75, length.out = 60),
    u_grid = seq(-2.6, 2.6, length.out = 75),
    draw_ids = NULL,
    scenario = "active",
    max_draw = 500L,
    include_process_uncertainty = FALSE,
    f0_fn = f0_1d
) {
  # Determine which draws to use
  if (is.null(draw_ids)) {
    draw_ids <- seq_len(nrow(fit$mcmc$samples_u))
  }
  
  # Subsample draws if requested
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
  out$truth <- f0_fn(out$x_raw, out$u, scenario = scenario)
  out$error <- out$mean - out$truth
  
  out
}


#' Extract MCMC trace data for hyperparameters (by chain)
#'
#' This function retrieves the posterior draws of hyperparameters (sigma_epsilon,
#' rho, theta_x, theta_u) separately for each MCMC chain, suitable for trace
#' plots or convergence diagnostics.
#'
#' @param fit A fitted model object. Must contain:
#'   \itemize{
#'     \item \code{fit$mcmc$samples_by_chain$sigma2}: list of sigma^2 samples per chain.
#'     \item \code{fit$mcmc$samples_by_chain$logtheta}: list of logtheta matrices per chain.
#'   }
#' @param n_calib Integer. Calibration sample size.
#' @param n_calib_label Character. Label for the calibration size.
#'   If \code{NULL}, defaults to \code{as.character(n_calib)}.
#' @return A \code{data.frame} with columns:
#'   \item{chain}{factor indicating chain number}
#'   \item{draw}{iteration index within chain}
#'   \item{sigma_epsilon}{posterior draw of sigma (sqrt(sigma^2))}
#'   \item{rho}{posterior draw of rho (exp(logtheta[1])), or \code{NA}}
#'   \item{theta_x}{posterior draw of theta_x (exp(logtheta[2])), or \code{NA}}
#'   \item{theta_u}{posterior draw of theta_u (exp(logtheta[3])), or \code{NA}}
#'   \item{n_calib}{calibration size}
#'   \item{n_calib_label}{label for calibration size}
#' @export
#'
#' @examples
#' \dontrun{
#'   df <- make_trace_df(fit, n_calib = 50)
#' }
make_trace_df <- function(
    fit,
    n_calib,
    n_calib_label = NULL
) {
  if (is.null(n_calib_label)) n_calib_label <- as.character(n_calib)
  
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
        n_calib_label = n_calib_label
      )
    })
  )
}


#' Extract posterior draws of hyperparameters (pooled across chains)
#'
#' This function stacks the hyperparameter posterior samples from all chains
#' into a single data frame, discarding chain information.
#'
#' @param fit A fitted model object. Must contain:
#'   \itemize{
#'     \item \code{fit$mcmc$samples_logtheta}: matrix of logtheta samples.
#'     \item \code{fit$mcmc$samples_sigma2}: vector of sigma^2 samples.
#'   }
#' @param n_calib Integer. Calibration sample size.
#' @param n_calib_label Character. Label for the calibration size.
#'   If \code{NULL}, defaults to \code{as.character(n_calib)}.
#' @return A \code{data.frame} with columns:
#'   \item{draw}{iteration index (pooled)}
#'   \item{n_calib}{calibration size}
#'   \item{n_calib_label}{label for calibration size}
#'   \item{sigma_epsilon}{sqrt(sigma^2) draw}
#'   \item{rho}{exp(logtheta[1]) draw, or \code{NA}}
#'   \item{theta_x}{exp(logtheta[2]) draw, or \code{NA}}
#'   \item{theta_u}{exp(logtheta[3]) draw, or \code{NA}}
#' @export
#'
#' @examples
#' \dontrun{
#'   df <- make_hyper_draw_df(fit, n_calib = 50)
#' }
make_hyper_draw_df <- function(
    fit,
    n_calib,
    n_calib_label = NULL
) {
  if (is.null(n_calib_label)) n_calib_label <- as.character(n_calib)
  
  logtheta_mat <- as.matrix(fit$mcmc$samples_logtheta)
  
  data.frame(
    draw = seq_len(nrow(logtheta_mat)),
    n_calib = n_calib,
    n_calib_label = n_calib_label,
    sigma_epsilon = sqrt(fit$mcmc$samples_sigma2),
    rho = if (ncol(logtheta_mat) >= 1L) exp(logtheta_mat[, 1]) else NA_real_,
    theta_x = if (ncol(logtheta_mat) >= 2L) exp(logtheta_mat[, 2]) else NA_real_,
    theta_u = if (ncol(logtheta_mat) >= 3L) exp(logtheta_mat[, 3]) else NA_real_
  )
}


#' Compute autocorrelation function (ACF) from long-format trace data
#'
#' This function takes a long-format data frame of MCMC traces and computes
#' the autocorrelation for each parameter–chain combination.
#'
#' @param df_long A data frame containing columns:
#'   \itemize{
#'     \item \code{parameter}: parameter name.
#'     \item \code{chain}: chain identifier.
#'     \item \code{value}: numeric value of the parameter at each draw.
#'   }
#' @param max_lag Integer. Maximum lag to compute. Default is \code{60L}.
#' @return A \code{data.frame} with columns:
#'   \item{parameter}{parameter name}
#'   \item{chain}{chain identifier}
#'   \item{lag}{lag index (0 = no lag)}
#'   \item{acf}{autocorrelation value}
#' @export
#'
#' @examples
#' \dontrun{
#'   df_acf <- make_acf_df(df_long, max_lag = 50)
#' }
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


#' Extract MCMC trace data for selected latent variables (by chain)
#'
#' This function retrieves the posterior draws of the latent variables
#' \code{u} for a specified set of observation indices, separately per chain.
#'
#' @param fit A fitted model object. Must contain:
#'   \itemize{
#'     \item \code{fit$mcmc$samples_by_chain$u}: list of u matrices per chain.
#'     \item \code{fit$data$c_ord}: class assignment for each observation.
#'     \item \code{fit$data$u_true}: true latent values.
#'   }
#' @param ids Integer vector. Observation indices (1‑based) to extract.
#'   Out-of-range indices are silently ignored.
#' @return A \code{data.frame} with columns:
#'   \item{chain}{factor indicating chain number}
#'   \item{draw}{iteration index within chain}
#'   \item{id}{observation index}
#'   \item{id_label}{label e.g., \code{"id=5, c=2"}}
#'   \item{u}{posterior draw of the latent variable}
#'   \item{true_u}{true latent value}
#' @export
#'
#' @examples
#' \dontrun{
#'   df <- make_u_trace_df(fit, ids = c(1, 5, 10))
#' }
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


#' Extract MCMC trace data for tau parameters (by chain)
#'
#' This function retrieves the posterior draws of the tau (cutpoint) parameters
#' separately per chain. Optionally, true tau values can be supplied for comparison.
#'
#' @param fit A fitted model object. Must contain:
#'   \itemize{
#'     \item \code{fit$mcmc$samples_by_chain} with one of:
#'       \code{tau}, \code{taus}, \code{cutpoints}, or \code{alpha}.
#'   }
#' @param tau_true_vals Numeric vector (optional). True values of the tau
#'   parameters, to be included in the output for comparison.
#'   If \code{NULL} (default), the \code{tau_true} column will be \code{NA}.
#' @return A \code{data.frame} with columns:
#'   \item{chain}{factor indicating chain number}
#'   \item{draw}{iteration index within chain}
#'   \item{tau_index}{factor indicating the tau component (1, 2, ...)}
#'   \item{tau}{posterior draw of the tau parameter}
#'   \item{tau_true}{true tau value (or \code{NA} if not supplied)}
#' @export
#'
#' @examples
#' \dontrun{
#'   df <- make_tau_trace_df(fit, tau_true_vals = c(-1, 0, 1))
#' }
make_tau_trace_df <- function(fit, tau_true_vals = NULL) {
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
  
  # If tau_true_vals is NULL, fill with NA
  if (is.null(tau_true_vals)) {
    tau_true_vals <- NA_real_
  }
  
  dplyr::bind_rows(
    lapply(seq_along(tau_list), function(cc) {
      obj <- tau_list[[cc]]
      
      if (is.null(dim(obj))) {
        # Convert vector to matrix
        if (length(obj) %% length(tau_true_vals) == 0L && length(tau_true_vals) > 1L) {
          mat <- matrix(as.numeric(obj), ncol = length(tau_true_vals), byrow = TRUE)
        } else {
          mat <- matrix(as.numeric(obj), ncol = 1L)
        }
      } else {
        mat <- as.matrix(obj)
      }
      
      # Transpose if needed to align rows = draws, cols = components
      if (ncol(mat) != length(tau_true_vals) && nrow(mat) == length(tau_true_vals)) {
        mat <- t(mat)
      }
      
      dplyr::bind_rows(
        lapply(seq_len(ncol(mat)), function(jj) {
          data.frame(
            chain = factor(cc),
            draw = seq_len(nrow(mat)),
            tau_index = factor(jj, levels = seq_len(ncol(mat))),
            tau = mat[, jj],
            tau_true = if (length(tau_true_vals) >= jj) tau_true_vals[jj] else NA_real_
          )
        })
      )
    })
  )
}
