############################################################
## Real-data interface for the multivariate ordinal EIV-GP
##
## This file deliberately does not install packages or download data.
## Package builds load the stable API before this convenience layer. When the
## file is sourced for development, load the same implementation explicitly.
############################################################

if (!exists("fit_eivgp", mode = "function")) {
  this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  code_dir <- if (is.null(this_file)) "." else dirname(normalizePath(this_file))
  if (!exists("mixedgp_parallel_lapply", mode = "function")) {
    source(file.path(code_dir, "00_parallel_utils.R"))
  }
  source(file.path(code_dir, "00_study2_functions.R"))
  source(file.path(code_dir, "00_public_api.R"))
}

mixedgp_expand_real_data_calibration <- function(X, U_obs, calib_idx, d) {
  n <- nrow(as.matrix(X))
  if (is.null(U_obs)) {
    if (is.null(d)) {
      stop("Supply d when U_obs is unavailable for every training unit.")
    }
    return(list(U_obs = NULL, d = d))
  }
  U_matrix <- as.matrix(U_obs)
  if (!is.numeric(U_matrix)) stop("U_obs must be numeric.")
  if (is.null(d)) d <- ncol(U_matrix)
  if (ncol(U_matrix) != d) stop("U_obs must have d columns.")
  if (is.null(calib_idx)) {
    if (nrow(U_matrix) != n) {
      stop("Compact U_obs requires calib_idx with one index per row.")
    }
    return(list(U_obs = U_matrix, d = d))
  }
  calib_idx <- as.integer(calib_idx)
  if (length(calib_idx) < 1L || anyNA(calib_idx) ||
      any(calib_idx < 1L | calib_idx > n) || anyDuplicated(calib_idx)) {
    stop("calib_idx must contain unique row indices between 1 and nrow(X).")
  }
  if (nrow(U_matrix) == length(calib_idx)) {
    if (any(!is.finite(U_matrix))) {
      stop("Every value in compact U_obs must be finite.")
    }
    U_full <- matrix(NA_real_, nrow = n, ncol = d)
    U_full[calib_idx, ] <- U_matrix
    colnames(U_full) <- colnames(U_matrix)
    return(list(U_obs = U_full, d = d))
  }
  if (nrow(U_matrix) != n) {
    stop("U_obs must have nrow(X) rows or length(calib_idx) compact rows.")
  }
  complete <- which(stats::complete.cases(U_matrix))
  partial <- apply(U_matrix, 1L, function(z) any(is.finite(z)) &&
    !all(is.finite(z)))
  if (any(partial) || !identical(sort(complete), sort(calib_idx))) {
    stop("For full-size U_obs, calib_idx must identify exactly its complete rows.")
  }
  list(U_obs = U_matrix, d = d)
}

#' Fit the multivariate EIV-GP to an observational dataset
#'
#' Convenience interface for a real-data analysis with quantitative inputs,
#' ordinal proxies, and optional calibrated latent inputs. New analyses may
#' also call [fit_eivgp()] directly with `engine = "multivariate"`.
#'
#' @param X Numeric training design matrix or data frame.
#' @param y Numeric response vector.
#' @param C Ordinal proxy matrix or data frame.
#' @param U_obs Optional calibrated latent inputs. Complete rows identify
#'   calibrated units; other rows must be entirely missing.
#' @param calib_idx Optional row indices for a compact calibrated `U_obs`.
#' @param d Latent dimension. Required when it cannot be inferred from
#'   `U_obs`.
#' @param ident Optional loading-identification rule.
#' @param ordinal_levels Optional low-to-high labels for the columns of `C`.
#' @param kernel Either `"se"` or `"matern"`.
#' @param matern_nu Matérn smoothness.
#' @param standardize_U Whether calibrated latent inputs are standardized;
#'   `NULL` uses the stable [fit_eivgp()] default.
#' @param store_scores Whether latent ordinal scores are retained.
#' @param parallel Whether independent MCMC chains may run in parallel.
#' @param n_cores Optional number of parallel chain workers.
#' @param ... Additional arguments passed to [fit_eivgp()].
#'
#' @return An object of class `eivgp_fit`.
#' @export
fit_eivgp_real_data <- function(X,
                                y,
                                C,
                                U_obs = NULL,
                                calib_idx = NULL,
                                d = NULL,
                                ident = NULL,
                                ordinal_levels = NULL,
                                kernel = c("se", "matern"),
                                matern_nu = 2.5,
                                standardize_U = NULL,
                                store_scores = FALSE,
                                parallel = TRUE,
                                n_cores = NULL,
                                ...) {
  calibration <- mixedgp_expand_real_data_calibration(X, U_obs, calib_idx, d)
  kernel <- match.arg(kernel)
  mixedgp_do_call(
    fit_eivgp,
    list(
      X = X, y = y, C = C, U_obs = calibration$U_obs,
      latent_dim = calibration$d, engine = "multivariate",
      ordinal_levels = ordinal_levels, kernel = kernel,
      matern_nu = matern_nu, ident = ident,
      standardize_U = standardize_U, store_scores = store_scores,
      parallel = parallel, n_cores = n_cores
    ),
    list(...)
  )
}

#' Predict responses for a real-data EIV-GP fit
#'
#' Draws from the posterior predictive distribution of a future response
#' conditional on new quantitative inputs and ordinal proxies.
#'
#' @param fit An `eivgp_fit` object from the multivariate engine.
#' @param X_new New quantitative inputs.
#' @param C_new New ordinal proxies.
#' @param draw_ids Optional posterior draw indices.
#' @param n_per_draw Number of predictive replicates per posterior draw.
#' @param U_new_obs Optional exact new latent inputs.
#' @param U_new_scale Scale of `U_new_obs`; see [predict_eivgp()].
#' @param joint Whether rows form one joint prediction target.
#' @param seed Optional random seed.
#' @param latent_sampler Method used to simulate latent inputs given proxies.
#' @param n_new_latent_gibbs Gibbs sweeps for the diagnostic Gibbs sampler.
#' @param rejection_batch_size Optional rejection-sampler batch size.
#' @param rejection_max_batches Maximum rejection-sampler batches.
#' @param ... Additional arguments passed to [predict_eivgp()].
#'
#' @return A posterior predictive draw matrix.
#' @export
predict_eivgp_y_given_xc <- function(fit,
                                     X_new,
                                     C_new,
                                     draw_ids = NULL,
                                     n_per_draw = 1L,
                                     U_new_obs = NULL,
                                     U_new_scale = c("auto", "raw", "model"),
                                     joint = FALSE,
                                     seed = NULL,
                                     latent_sampler = c(
                                       "minimax_tilting", "rejection", "gibbs"
                                     ),
                                     n_new_latent_gibbs = 100L,
                                     rejection_batch_size = NULL,
                                     rejection_max_batches = 1000L,
                                     ...) {
  common <- list(
    object = fit, new_X = X_new, new_C = C_new, new_U = U_new_obs,
    target = "response", draw_ids = draw_ids, n_per_draw = n_per_draw,
    joint = joint, seed = seed
  )
  if (!is.null(U_new_obs)) common$new_U_scale <- match.arg(U_new_scale)
  if (!missing(latent_sampler)) common$latent_sampler <- match.arg(latent_sampler)
  if (!missing(n_new_latent_gibbs)) {
    common$n_new_latent_gibbs <- n_new_latent_gibbs
  }
  if (!missing(rejection_batch_size)) {
    common$rejection_batch_size <- rejection_batch_size
  }
  if (!missing(rejection_max_batches)) {
    common$rejection_max_batches <- rejection_max_batches
  }
  mixedgp_do_call(predict_eivgp, common, list(...))
}

#' Infer the latent response surface for a real-data EIV-GP fit
#'
#' @param fit An `eivgp_fit` object from the multivariate engine.
#' @param X_new New quantitative inputs.
#' @param U_new Exact latent inputs on the requested fitted scale.
#' @param U_new_scale Scale of `U_new`; see [predict_eivgp()].
#' @param draw_ids Optional posterior draw indices.
#' @param n_per_draw Number of replicates per posterior draw.
#' @param joint Whether rows form one joint prediction target.
#' @param seed Optional random seed.
#' @param include_process_uncertainty Whether to include GP process
#'   uncertainty.
#' @param ... Additional arguments passed to [predict_eivgp()].
#'
#' @return A posterior draw matrix for the latent response surface.
#' @export
predict_eivgp_f_given_xu <- function(fit,
                                     X_new,
                                     U_new,
                                     U_new_scale = c("auto", "raw", "model"),
                                     draw_ids = NULL,
                                     n_per_draw = 1L,
                                     joint = FALSE,
                                     seed = NULL,
                                     include_process_uncertainty = TRUE,
                                     ...) {
  mixedgp_do_call(
    predict_eivgp,
    list(
      object = fit, new_X = X_new, new_U = U_new,
      new_U_scale = match.arg(U_new_scale), target = "surface",
      draw_ids = draw_ids, n_per_draw = n_per_draw, joint = joint,
      seed = seed,
      include_process_uncertainty = include_process_uncertainty
    ),
    list(...)
  )
}

#' Infer the observed-input mean for a real-data EIV-GP fit
#'
#' @param fit An `eivgp_fit` object from the multivariate engine.
#' @param X_new New quantitative inputs.
#' @param C_new New ordinal proxies.
#' @param draw_ids Optional posterior draw indices.
#' @param n_per_draw Number of replicates per posterior draw.
#' @param n_latent Number of latent draws used for integration.
#' @param joint Whether rows form one joint prediction target.
#' @param seed Optional random seed.
#' @param include_process_uncertainty Whether to include GP process
#'   uncertainty.
#' @param latent_sampler Method used to simulate latent inputs given proxies.
#' @param n_new_latent_gibbs Gibbs sweeps for the diagnostic Gibbs sampler.
#' @param rejection_batch_size Optional rejection-sampler batch size.
#' @param rejection_max_batches Maximum rejection-sampler batches.
#' @param ... Additional arguments passed to [predict_eivgp()].
#'
#' @return A posterior draw matrix for `m(x,c)`.
#' @export
predict_eivgp_m_given_xc <- function(fit,
                                     X_new,
                                     C_new,
                                     draw_ids = NULL,
                                     n_per_draw = 1L,
                                     n_latent = 256L,
                                     joint = FALSE,
                                     seed = NULL,
                                     include_process_uncertainty = TRUE,
                                     latent_sampler = c(
                                       "minimax_tilting", "rejection", "gibbs"
                                     ),
                                     n_new_latent_gibbs = 100L,
                                     rejection_batch_size = NULL,
                                     rejection_max_batches = 1000L,
                                     ...) {
  common <- list(
    object = fit, new_X = X_new, new_C = C_new, target = "mean",
    draw_ids = draw_ids, n_per_draw = n_per_draw, n_latent = n_latent,
    joint = joint, seed = seed,
    include_process_uncertainty = include_process_uncertainty
  )
  if (!missing(latent_sampler)) common$latent_sampler <- match.arg(latent_sampler)
  if (!missing(n_new_latent_gibbs)) {
    common$n_new_latent_gibbs <- n_new_latent_gibbs
  }
  if (!missing(rejection_batch_size)) {
    common$rejection_batch_size <- rejection_batch_size
  }
  if (!missing(rejection_max_batches)) {
    common$rejection_max_batches <- rejection_max_batches
  }
  mixedgp_do_call(predict_eivgp, common, list(...))
}

## Expected data layout:
##   X: numeric n x p matrix/data frame (complete)
##   y: numeric response vector (complete)
##   C: n x q ordered-factor or numeric ordinal matrix/data frame (complete)
##   U_obs: numeric n x d object; complete rows are calibrated units and
##          entirely missing rows are uncalibrated units
##
## For character or unordered-factor proxies, provide ordinal_levels as a list
## in low-to-high order, one vector per column of C.
##
## Example call:
## fit <- fit_eivgp_real_data(
##   X = X,
##   y = y,
##   C = C,
##   U_obs = U_obs,
##   kernel = "matern",
##   matern_nu = 2.5,
##   n_iter = 6000,
##   burn = 2000,
##   thin = 2,
##   n_chains = 4,
##   preset = "thorough",
##   parallel_chains = TRUE
## )
