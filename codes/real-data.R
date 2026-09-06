############################################################
## Real-data interface for the multivariate ordinal EIV-GP
##
## This file deliberately does not install packages or download data.
## Source 00_study2_functions.R first, or source this file from the
## revision/codes directory and it will load the implementation.
############################################################

if (!exists("fit_eivgp_ordprobit_fb")) {
  this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  code_dir <- if (is.null(this_file)) "." else dirname(normalizePath(this_file))
  if (!exists("mixedgp_parallel_lapply", mode = "function")) {
    source(file.path(code_dir, "00_parallel_utils.R"))
  }
  source(file.path(code_dir, "00_study2_functions.R"))
}

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
                                standardize_U = TRUE,
                                store_scores = FALSE,
                                ...) {
  if (is.null(d)) {
    if (is.null(U_obs)) {
      stop("Supply d when U_obs is unavailable for every training unit.")
    }
    d <- ncol(as.matrix(U_obs))
  }

  if (is.null(ident)) {
    n <- nrow(as.matrix(X))
    U_full <- matrix(NA_real_, n, d)
    inferred_calib <- integer(0)
    if (!is.null(U_obs)) {
      U_matrix <- as.matrix(U_obs)
      if (nrow(U_matrix) == n) {
        inferred_calib <- which(stats::complete.cases(U_matrix))
        U_full[inferred_calib, ] <- U_matrix[inferred_calib, , drop = FALSE]
      } else if (!is.null(calib_idx) && nrow(U_matrix) == length(calib_idx)) {
        inferred_calib <- as.integer(calib_idx)
        U_full[inferred_calib, ] <- U_matrix
      }
    }
    anchor_status <- mixedgp_latent_anchor_status(
      U_full, inferred_calib, d = d
    )
    ident <- if (isTRUE(anchor_status$anchored)) "none" else "lower_triangular"
  }

  fit <- fit_eivgp_ordprobit_fb(
    X_raw = X,
    y_raw = y,
    C_ord = C,
    U_obs = U_obs,
    calib_idx = calib_idx,
    d = d,
    ident = ident,
    kernel = kernel,
    matern_nu = matern_nu,
    standardize_U = standardize_U,
    ordinal_levels = ordinal_levels,
    store_scores = store_scores,
    ...
  )
  fit$interface <- list(
    engine = "multivariate",
    ordinal_level_maps = fit$data$C_level_maps
  )
  class(fit) <- unique(c("eivgp_fit", class(fit)))
  fit
}

predict_eivgp_y_given_xc <- function(fit,
                                     X_new,
                                     C_new,
                                     draw_ids = NULL,
                                     n_per_draw = 1L,
                                     U_new_obs = NULL,
                                     joint = FALSE,
                                     seed = NULL,
                                     latent_sampler = c(
                                       "minimax_tilting", "rejection", "gibbs"
                                     ),
                                     n_new_latent_gibbs = 100L,
                                     rejection_batch_size = NULL,
                                     rejection_max_batches = 1000L) {
  latent_sampler <- match.arg(latent_sampler)
  sample_eiv_test_y_ordprobit_fb(
    X_test_raw = X_new,
    C_test = C_new,
    fit_obj = fit,
    draw_ids = draw_ids,
    n_per_draw = n_per_draw,
    latent_sampler = latent_sampler,
    n_new_latent_gibbs = n_new_latent_gibbs,
    rejection_batch_size = rejection_batch_size,
    rejection_max_batches = rejection_max_batches,
    U_test_obs = U_new_obs,
    predictive_target = "response",
    joint = joint,
    seed = seed
  )
}

predict_eivgp_f_given_xu <- function(fit,
                                     X_new,
                                     U_new,
                                     draw_ids = NULL,
                                     n_per_draw = 1L,
                                     joint = FALSE,
                                     seed = NULL,
                                     include_process_uncertainty = TRUE) {
  sample_eiv_f_given_xu_fb(
    X_test_raw = X_new,
    U_test_raw = U_new,
    fit_obj = fit,
    draw_ids = draw_ids,
    n_per_draw = n_per_draw,
    include_process_uncertainty = include_process_uncertainty,
    joint = joint,
    seed = seed
  )
}

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
                                     rejection_max_batches = 1000L) {
  latent_sampler <- match.arg(latent_sampler)
  sample_eiv_m_given_xc_fb(
    X_test_raw = X_new,
    C_test = C_new,
    fit_obj = fit,
    draw_ids = draw_ids,
    n_per_draw = n_per_draw,
    n_latent = n_latent,
    include_process_uncertainty = include_process_uncertainty,
    joint = joint,
    latent_sampler = latent_sampler,
    n_new_latent_gibbs = n_new_latent_gibbs,
    rejection_batch_size = rejection_batch_size,
    rejection_max_batches = rejection_max_batches,
    seed = seed
  )
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
