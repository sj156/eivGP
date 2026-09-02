############################################################
## Package-ready public fitting interface
############################################################

mixedgp_encode_ordinal <- function(C, m_vec = NULL, ordinal_levels = NULL) {
  C_df <- if (is.data.frame(C)) C else as.data.frame(C, check.names = FALSE)
  prepared <- prepare_ordinal_matrix(
    C_df,
    m_vec = m_vec,
    level_maps = ordinal_levels,
    name = "C"
  )
  names(C_df) <- prepared$column_names
  list(
    C = prepared$C,
    raw = C_df,
    level_maps = prepared$level_maps,
    m_vec = prepared$m_vec
  )
}

mixedgp_validate_named_dots <- function(dots, argument = "...") {
  if (!is.list(dots)) stop(argument, " must be supplied as a list internally.")
  if (length(dots) > 0L &&
      (is.null(names(dots)) || anyNA(names(dots)) ||
       any(!nzchar(names(dots))))) {
    stop("Every argument passed through ", argument, " must be named.")
  }
  if (length(dots) > 0L && anyDuplicated(names(dots))) {
    stop("Arguments passed through ", argument, " must have unique names.")
  }
  invisible(dots)
}

mixedgp_do_call <- function(FUN, common, dots) {
  mixedgp_validate_named_dots(dots)
  duplicate <- intersect(names(common), names(dots))
  if (length(duplicate) > 0L) {
    stop("Do not duplicate public arguments through ...: ",
         paste(duplicate, collapse = ", "))
  }
  do.call(FUN, c(common, dots))
}

mixedgp_resolve_latent_scale <- function(object,
                                         scale = c("auto", "raw", "model"),
                                         argument = "scale") {
  scale <- match.arg(scale)
  anchored <- isTRUE(object$data$latent_scale_anchored)
  if (scale == "raw" && !anchored) {
    stop(
      argument,
      "='raw' requires calibration data that anchor the latent coordinates; ",
      "use 'model' for the identified working coordinates."
    )
  }
  internal <- if (scale == "auto") {
    if (anchored) "raw" else "model"
  } else {
    scale
  }
  list(
    internal = internal,
    label = if (!anchored && internal == "model") "working" else internal,
    anchored = anchored
  )
}

#' Fit the EIV-GP for ordinal mixed inputs
#'
#' This is the stable package-facing interface for both the univariate and
#' multivariate ordinal-input samplers. Independent chains use forked R
#' processes on macOS/Linux and automatically fall back to serial execution.
#'
#' @param X Numeric training design matrix.
#' @param y Numeric response vector.
#' @param C Ordinal inputs as factors or positive integer codes.
#' @param U_obs Optional calibrated latent inputs; unavailable rows are `NA`.
#' @param latent_dim Latent dimension. Inferred from informative `U_obs` when
#'   possible. It must be supplied explicitly for the multivariate engine when
#'   no calibrated latent input is available.
#' @param engine The scientific measurement model: `"univariate"` uses
#'   deterministic thresholding of one latent input, whereas `"multivariate"`
#'   uses the noisy ordinal-probit model. It must be selected explicitly.
#' @param m_vec Number of levels for each ordinal input.
#' @param ordinal_levels Optional ordered level labels, one vector per ordinal
#'   input. Unordered factors require this argument.
#' @param kernel `"se"` or `"matern"`.
#' @param matern_nu Smoothness of the Matérn kernel.
#' @param ident Multivariate loading identification rule. By default `"none"`
#'   is used only when the calibration design has full affine rank `d + 1`;
#'   otherwise `"lower_triangular"` is used.
#' @param standardize_U Whether calibrated latent inputs are standardized before
#'   fitting. The default is `TRUE` when calibration is available.
#' @param store_scores Whether to retain ordinal-probit score draws.
#' @param parallel Whether to run independent MCMC chains in parallel.
#' @param n_cores Number of forked workers; `NULL` uses package settings.
#' @param ... Additional sampler arguments such as `n_iter`, `burn`, and
#'   `n_chains`.
#' @return An EIV-GP fitted object.
#' @export
fit_eivgp <- function(X,
                      y,
                      C,
                      U_obs = NULL,
                      latent_dim = NULL,
                      engine = NULL,
                      m_vec = NULL,
                      ordinal_levels = NULL,
                      kernel = c("se", "matern"),
                      matern_nu = 2.5,
                      ident = NULL,
                      standardize_U = NULL,
                      store_scores = FALSE,
                      parallel = TRUE,
                      n_cores = NULL,
                      ...) {
  fit_call <- match.call()
  caller_rng <- mixedgp_rng_state()
  on.exit(mixedgp_restore_rng_state(caller_rng), add = TRUE)
  parallel <- mixedgp_validate_flag(parallel, "parallel")
  store_scores <- mixedgp_validate_flag(store_scores, "store_scores")
  standardize_U <- mixedgp_validate_flag(
    standardize_U, "standardize_U", allow_null = TRUE
  )
  if (is.null(engine)) {
    stop(
      "Select engine explicitly: 'univariate' for deterministic thresholding ",
      "or 'multivariate' for the noisy ordinal-probit measurement model."
    )
  }
  engine <- match.arg(engine, c("univariate", "multivariate"))
  kernel <- match.arg(kernel)
  X <- as_numeric_matrix_strict(X, "X")
  y_mat <- as_numeric_matrix_strict(y, "y")
  if (ncol(y_mat) != 1L) stop("y must be a univariate response.")
  y <- as.numeric(y_mat[, 1L])
  encoded <- mixedgp_encode_ordinal(C, m_vec, ordinal_levels)
  C_mat <- encoded$C
  if (nrow(X) != length(y) || nrow(C_mat) != length(y)) {
    stop("X, y, and C must have the same number of observations.")
  }
  m_vec <- encoded$m_vec
  U_input_names <- if (is.null(U_obs) || is.null(dim(U_obs))) {
    NULL
  } else {
    colnames(U_obs)
  }
  U_checked <- if (is.null(U_obs)) {
    NULL
  } else {
    as_numeric_matrix_strict(
      U_obs, "U_obs", nrow_expected = length(y), allow_na = TRUE
    )
  }

  if (engine == "multivariate" && is.null(latent_dim)) {
    if (is.null(U_checked)) {
      stop(
        "latent_dim must be supplied explicitly for the multivariate engine ",
        "when no calibrated U_obs is available."
      )
    }
    finite_probe <- is.finite(U_checked)
    if (length(finite_probe) == 0L || !any(finite_probe)) {
      stop(
        "latent_dim must be supplied explicitly for the multivariate engine ",
        "when U_obs contains no calibration values."
      )
    }
    latent_dim <- ncol(U_checked)
  }
  dots <- list(...)
  mixedgp_validate_named_dots(dots)
  if ("seed" %in% names(dots)) {
    dots$seed <- mixedgp_as_integer_strict(
      dots$seed, "seed", min_value = 0L, length_expected = 1L
    )
  }

  if (engine == "univariate") {
    if (ncol(C_mat) != 1L || length(m_vec) != 1L) {
      stop("The univariate engine requires exactly one ordinal input.")
    }
    u_obs <- if (is.null(U_checked)) {
      rep(NA_real_, length(y))
    } else {
      if (ncol(U_checked) != 1L) {
        stop("For the univariate engine, U_obs must have one column.")
      }
      as.numeric(U_checked[, 1L])
    }
    if (length(u_obs) != length(y)) {
      stop("For the univariate engine, U_obs must have one value per row.")
    }
    if (!is.null(ident) || isTRUE(store_scores)) {
      stop("ident and store_scores apply only to the multivariate engine.")
    }
    calibrated <- which(is.finite(u_obs))
    if (is.null(standardize_U)) standardize_U <- length(calibrated) > 0L
    u_center <- 0
    u_scale <- 1
    if (isTRUE(standardize_U) && length(calibrated) > 0L) {
      u_center <- mean(u_obs[calibrated])
      u_scale <- stats::sd(u_obs[calibrated])
      if (!is.finite(u_scale) || u_scale <= 0) u_scale <- 1
    }
    u_obs_raw <- u_obs
    if (length(calibrated) > 0L) {
      u_obs[calibrated] <- (u_obs[calibrated] - u_center) / u_scale
    }
    fit <- mixedgp_do_call(
      fit_eivgp_1d,
      list(
        x_raw = X, y_raw = y, c_ord = C_mat[, 1L],
        u_obs = u_obs, calib_idx = which(is.finite(u_obs)), m = m_vec,
        kernel = kernel, matern_nu = matern_nu,
        parallel_chains = isTRUE(parallel), n_cores = n_cores
      ),
      dots
    )
    fit$data$u_obs_raw <- u_obs_raw
    fit$data$U_center <- u_center
    fit$data$U_scale <- u_scale
    fit$data$standardize_U <- isTRUE(standardize_U)
    fit$data$U_names <- U_input_names
    if (is.null(fit$data$U_names) || length(fit$data$U_names) != 1L ||
        !nzchar(fit$data$U_names)) {
      fit$data$U_names <- "u1"
    }
  } else {
    latent_dim <- mixedgp_as_integer_strict(
      latent_dim, "latent_dim", min_value = 1L, length_expected = 1L
    )
    U_full <- if (is.null(U_checked)) {
      matrix(NA_real_, length(y), latent_dim)
    } else {
      U_checked
    }
    if (!identical(dim(U_full), c(length(y), latent_dim))) {
      stop("U_obs must be an n by latent_dim matrix.")
    }
    partly_observed <- apply(U_full, 1L, function(z) any(is.finite(z)) &&
      !all(is.finite(z)))
    if (any(partly_observed)) {
      stop("Each calibrated row of U_obs must observe every latent coordinate.")
    }
    calib_idx <- which(apply(U_full, 1L, function(z) all(is.finite(z))))
    anchor_status <- mixedgp_latent_anchor_status(
      U_full, calib_idx, d = latent_dim
    )
    if (is.null(ident)) {
      ident <- if (isTRUE(anchor_status$anchored)) {
        "none"
      } else {
        "lower_triangular"
      }
    }
    ident <- match.arg(ident, c("lower_triangular", "none"))
    if (ident == "none" && !isTRUE(anchor_status$anchored)) {
      stop(
        "ident='none' requires a calibration design with full affine rank ",
        "d + 1; use ident='lower_triangular' for identified working ",
        "coordinates."
      )
    }
    if (length(calib_idx) == 0L) {
      warning(
        "No calibration data anchor the latent coordinates; using identified ",
        "working coordinates, and raw/physical-scale latent outputs remain ",
        "unavailable.",
        call. = FALSE
      )
    }
    if (is.null(standardize_U)) standardize_U <- length(calib_idx) > 0L
    fit <- mixedgp_do_call(
      fit_eivgp_ordprobit_fb,
      list(
        X_raw = X, y_raw = y, C_ord = encoded$raw, U_obs = U_full,
        calib_idx = calib_idx,
        d = latent_dim, m_vec = m_vec,
        ordinal_levels = encoded$level_maps,
        kernel = kernel, matern_nu = matern_nu,
        ident = ident, standardize_U = isTRUE(standardize_U),
        store_scores = isTRUE(store_scores),
        parallel_chains = isTRUE(parallel), n_cores = n_cores
      ),
      dots
    )
  }
  fit$data$C_names <- colnames(C_mat)
  fit$data$C_level_maps <- encoded$level_maps
  fit$call <- fit_call
  fit$interface <- list(
    engine = engine,
    parallel = isTRUE(parallel),
    n_cores = fit$diagnostics$summary$parallel_cores,
    ordinal_level_maps = encoded$level_maps
  )
  class(fit) <- unique(c("eivgp_fit", class(fit)))
  fit
}

mixedgp_first_value <- function(x, default = NA) {
  if (is.null(x) || length(x) == 0L) return(default)
  x[[1L]]
}

mixedgp_diagnostic_value <- function(diagnostics, name, default = NA) {
  if (is.null(diagnostics) || !name %in% names(diagnostics) ||
      length(diagnostics[[name]]) == 0L) {
    return(default)
  }
  diagnostics[[name]][[1L]]
}

mixedgp_fit_engine <- function(object) {
  engine <- object$interface$engine
  if (is.null(engine)) {
    engine <- if (!is.null(object$data$d)) "multivariate" else "univariate"
  }
  engine <- as.character(engine[[1L]])
  if (!engine %in% c("univariate", "multivariate")) {
    stop("Malformed eivgp_fit object: unknown engine.")
  }
  engine
}

mixedgp_fit_dimensions <- function(object, engine) {
  data <- object$data
  X <- if (engine == "univariate") {
    if (!is.null(data$x_raw)) data$x_raw else data$x
  } else {
    if (!is.null(data$X_raw)) data$X_raw else data$X
  }
  n <- if (!is.null(X)) nrow(as.matrix(X)) else length(data$y)
  p <- if (!is.null(X)) ncol(as.matrix(X)) else if (engine == "univariate") {
    mixedgp_first_value(data$p_x, NA_integer_)
  } else {
    mixedgp_first_value(data$p, NA_integer_)
  }
  q <- if (engine == "univariate") {
    1L
  } else if (!is.null(data$q)) {
    data$q
  } else {
    ncol(as.matrix(data$C_ord))
  }
  d <- if (engine == "univariate") {
    1L
  } else if (!is.null(data$d)) {
    data$d
  } else {
    dim(object$mcmc$samples_U)[3L]
  }
  c(
    n = as.integer(mixedgp_first_value(n, NA_integer_)),
    p = as.integer(mixedgp_first_value(p, NA_integer_)),
    q = as.integer(mixedgp_first_value(q, NA_integer_)),
    d = as.integer(mixedgp_first_value(d, NA_integer_))
  )
}

mixedgp_kernel_label <- function(kernel) {
  name <- tolower(as.character(mixedgp_first_value(kernel$name, "unknown")))
  if (name == "se") return("squared-exponential")
  if (name == "matern") {
    nu <- mixedgp_first_value(kernel$matern_nu, NA_real_)
    return(if (is.finite(nu)) paste0("Matern (nu = ", format(nu), ")") else
      "Matern")
  }
  name
}

#' Summarize a fitted EIV-GP model
#'
#' The summary reports the fitted measurement-model engine, problem dimensions,
#' calibration rank and latent-scale interpretation, kernel, MCMC configuration,
#' and the convergence summaries already computed by the sampler.
#'
#' @param object An object returned by [fit_eivgp()].
#' @param ... Reserved for future summary controls.
#' @return An object of class `summary_eivgp_fit`.
#' @export
summary.eivgp_fit <- function(object, ...) {
  if (!inherits(object, "eivgp_fit")) {
    stop("object must be returned by fit_eivgp().")
  }
  engine <- mixedgp_fit_engine(object)
  dimensions <- mixedgp_fit_dimensions(object, engine)
  calib_idx <- object$data$calib_idx
  if (is.null(calib_idx)) calib_idx <- integer(0)
  n_calibrated <- length(calib_idx)
  anchored <- isTRUE(object$data$latent_scale_anchored)
  required_rank <- mixedgp_first_value(
    object$data$latent_anchor_required_rank,
    dimensions[["d"]] + 1L
  )
  affine_rank <- mixedgp_first_value(
    object$data$latent_anchor_rank,
    if (n_calibrated == 0L) 0L else NA_integer_
  )
  calibration_status <- if (n_calibrated == 0L) {
    "none"
  } else if (anchored) {
    "affine_anchored"
  } else {
    "rank_deficient"
  }
  interpretation <- if (n_calibrated == 0L) {
    paste(
      "No calibration data were observed. The observed-input mean m(x,c) and",
      "prediction given (x,c) remain model-defined; latent coordinates,",
      "f(x,u), and latent-input imputations are identified only on the model's",
      "working scale, so raw/physical-scale latent claims are unavailable."
    )
  } else if (!anchored) {
    paste0(
      "Calibration data inform U but have affine rank ", affine_rank, " of ",
      required_rank, ". Latent outputs remain on the identified working scale; ",
      "raw/physical-scale latent claims are unavailable."
    )
  } else {
    paste0(
      "The calibration design has full affine rank ", affine_rank, " of ",
      required_rank, "; raw/physical latent coordinates are anchored."
    )
  }

  diagnostics <- object$diagnostics$summary
  if (is.null(diagnostics)) diagnostics <- data.frame()
  if (is.data.frame(diagnostics) && nrow(diagnostics) > 1L) {
    diagnostics <- diagnostics[1L, , drop = FALSE]
  }
  convergence_names <- names(diagnostics)[
    grepl("rhat|ess", names(diagnostics), ignore.case = TRUE) &
      !grepl("accept|eval|total", names(diagnostics), ignore.case = TRUE)
  ]
  convergence <- lapply(
    convergence_names,
    function(name) mixedgp_diagnostic_value(diagnostics, name)
  )
  names(convergence) <- convergence_names
  mcmc_names <- c(
    "sampler_strategy", "n_chains", "parallel_backend", "parallel_cores",
    "n_iter", "burn", "thin", "saved_per_chain", "total_saved_draws",
    "time_seconds"
  )
  mcmc <- lapply(
    mcmc_names,
    function(name) mixedgp_diagnostic_value(diagnostics, name)
  )
  names(mcmc) <- mcmc_names

  structure(
    list(
      call = object$call,
      engine = engine,
      dimensions = dimensions,
      calibration = list(
        n = as.integer(n_calibrated),
        status = calibration_status,
        affine_rank = as.integer(affine_rank),
        required_rank = as.integer(required_rank),
        anchored = anchored,
        output_scale = if (anchored) "raw" else "working",
        interpretation = interpretation
      ),
      identification = mixedgp_first_value(object$data$ident, NA_character_),
      kernel = list(
        name = mixedgp_first_value(object$kernel$name, "unknown"),
        matern_nu = mixedgp_first_value(object$kernel$matern_nu, NA_real_),
        label = mixedgp_kernel_label(object$kernel)
      ),
      mcmc = mcmc,
      convergence = convergence,
      diagnostic_tables = object$diagnostics[
        setdiff(names(object$diagnostics), "summary")
      ]
    ),
    class = c("summary_eivgp_fit", "list")
  )
}

#' Print a fitted EIV-GP model
#'
#' @param x An object returned by [fit_eivgp()].
#' @param ... Reserved for future print controls.
#' @return `x`, invisibly.
#' @export
print.eivgp_fit <- function(x, ...) {
  fit_summary <- summary.eivgp_fit(x)
  dims <- fit_summary$dimensions
  calibration <- fit_summary$calibration
  cat("<eivgp_fit>\n")
  cat(" Engine: ", fit_summary$engine, "\n", sep = "")
  cat(
    " Data: ", dims[["n"]], " rows; ", dims[["p"]], " numeric, ",
    dims[["q"]], " ordinal, and ", dims[["d"]], " latent dimension(s)\n",
    sep = ""
  )
  cat(
    " Calibration: ", calibration$n, " rows; affine rank ",
    calibration$affine_rank, "/", calibration$required_rank, "; ",
    calibration$output_scale, " latent scale\n",
    sep = ""
  )
  cat(" Kernel: ", fit_summary$kernel$label, "\n", sep = "")
  saved <- fit_summary$mcmc$total_saved_draws
  if (length(saved) == 1L && is.finite(saved)) {
    cat(" Posterior draws: ", format(saved, scientific = FALSE), "\n", sep = "")
  }
  invisible(x)
}

#' @rdname summary.eivgp_fit
#' @param x A fitted-model summary returned by `summary()`.
#' @return `x`, invisibly, for the print method.
#' @export
print.summary_eivgp_fit <- function(x, ...) {
  dims <- x$dimensions
  calibration <- x$calibration
  cat("Ordinal mixed-input EIV-GP fit\n")
  cat(" Engine: ", x$engine, "\n", sep = "")
  cat(
    " Training design: n = ", dims[["n"]], ", p = ", dims[["p"]],
    ", q = ", dims[["q"]], ", d = ", dims[["d"]], "\n",
    sep = ""
  )
  cat(
    " Calibration: ", calibration$n, " rows; affine rank ",
    calibration$affine_rank, "/", calibration$required_rank,
    "; status = ", calibration$status, "\n",
    sep = ""
  )
  if (!is.na(x$identification)) {
    cat(" Identification: ", x$identification, "\n", sep = "")
  }
  cat(" Kernel: ", x$kernel$label, "\n", sep = "")

  mcmc <- x$mcmc
  if (is.finite(mcmc$n_chains) && is.finite(mcmc$n_iter)) {
    cat(
      " MCMC: ", mcmc$n_chains, " chain(s) x ", mcmc$n_iter,
      " iterations; burn = ", mcmc$burn, ", thin = ", mcmc$thin,
      "; saved = ", mcmc$total_saved_draws, "\n",
      sep = ""
    )
  }
  convergence <- unlist(x$convergence, use.names = TRUE)
  rhat <- convergence[
    grepl("rhat", names(convergence), ignore.case = TRUE) &
      is.finite(convergence)
  ]
  ess <- convergence[
    grepl("ess", names(convergence), ignore.case = TRUE) &
      is.finite(convergence)
  ]
  if (length(rhat) > 0L || length(ess) > 0L) {
    pieces <- character(0)
    if (length(rhat) > 0L) {
      pieces <- c(pieces, paste0("largest reported R-hat = ", format(max(rhat))))
    }
    if (length(ess) > 0L) {
      pieces <- c(pieces, paste0("smallest reported ESS = ", format(min(ess))))
    }
    cat(" Convergence: ", paste(pieces, collapse = "; "), "\n", sep = "")
  }
  cat(" Latent-scale interpretation: ", calibration$interpretation, "\n", sep = "")
  invisible(x)
}

#' Predict from a fitted EIV-GP model
#'
#' This S3 method delegates to [predict_eivgp()] without changing its posterior
#' predictive semantics.
#'
#' @param object An object returned by [fit_eivgp()].
#' @param ... Arguments passed to [predict_eivgp()], including `new_X`, `new_C`,
#'   `new_U`, and `target`.
#' @return Posterior predictive or latent-function draws from [predict_eivgp()].
#' @export
predict.eivgp_fit <- function(object, ...) {
  predict_eivgp(object, ...)
}

#' Draw from an EIV-GP predictive distribution
#'
#' @param object An object returned by [fit_eivgp()].
#' @param new_X Numeric prediction design matrix.
#' @param new_C Ordinal inputs for prediction of `Y` given `(X,C)` or
#'   inference on the observed-input mean `m(x,c)`. It may be omitted for
#'   response prediction only when every row of `new_U` is fully observed.
#' @param new_U Optional exact continuous latent inputs. For response prediction,
#'   complete rows condition on the supplied value and missing rows integrate
#'   over the ordinal measurement model. It is required for `target="surface"`
#'   and is not used for `target="mean"`.
#' @param new_U_scale Scale of supplied `new_U`. `"auto"` means the physical
#'   raw scale when calibration anchors the latent coordinates and the model's
#'   identified working scale otherwise. Requesting `"raw"` without a
#'   full-affine-rank calibration design is an error.
#' @param target `"response"` for a future noisy response, `"mean"` for
#'   `m(x,c)`, or `"surface"` for `f(x,u)`.
#' @param draw_ids Saved posterior draws to use; `NULL` uses every draw.
#' @param n_per_draw Number of predictive draws per posterior draw.
#' @param joint Whether each returned row is a coherent joint draw over all
#'   evaluation points. The default is faster and gives correct pointwise
#'   marginals, but `joint=TRUE` is required for surface realizations or
#'   simultaneous functionals.
#' @param include_process_uncertainty Whether `mean` and `surface` draws include
#'   conditional GP process uncertainty. When `FALSE`, conditional means are
#'   Rao--Blackwellized within each posterior state.
#' @param n_latent Number of draws from `U` given `C` used to approximate the
#'   integrals defining `m(x,c)` within each posterior state.
#' @param latent_sampler Prospective latent-input sampler for the multivariate
#'   engine. The default is a minimax-tilted accept--reject sampler that is exact
#'   on successful completion. Exact prior rejection is an alternative, and
#'   finite Gibbs is an explicitly requested diagnostic approximation.
#' @param n_new_latent_gibbs Number of Gibbs sweeps when
#'   `latent_sampler="gibbs"`.
#' @param rejection_batch_size,rejection_max_batches Controls for exact
#'   prospective latent-input rejection sampling.
#' @param seed Optional predictive seed.
#' @param ... Additional prediction controls.
#' @return A numeric matrix with one row per retained posterior draw times
#'   `n_per_draw` and one column per evaluation row. Attributes record the
#'   estimand, whether columns form coherent joint draws, and (when applicable)
#'   the latent-input scale and anchoring status. Row order follows `draw_ids`,
#'   with the `n_per_draw` replicates for each draw contiguous.
#' @export
predict_eivgp <- function(object,
                          new_X,
                          new_C = NULL,
                          new_U = NULL,
                          new_U_scale = c("auto", "raw", "model"),
                          target = c("response", "mean", "surface"),
                          draw_ids = NULL,
                          n_per_draw = 1L,
                          joint = FALSE,
                          include_process_uncertainty = TRUE,
                          n_latent = 256L,
                          latent_sampler = c(
                            "minimax_tilting", "rejection", "gibbs"
                          ),
                          n_new_latent_gibbs = 100L,
                          rejection_batch_size = NULL,
                          rejection_max_batches = 1000L,
                          seed = NULL,
                          ...) {
  if (!inherits(object, "eivgp_fit")) {
    stop("object must be returned by fit_eivgp().")
  }
  supplied <- list(
    new_U_scale = !missing(new_U_scale),
    include_process_uncertainty = !missing(include_process_uncertainty),
    n_latent = !missing(n_latent),
    latent_sampler = !missing(latent_sampler),
    n_new_latent_gibbs = !missing(n_new_latent_gibbs),
    rejection_batch_size = !missing(rejection_batch_size),
    rejection_max_batches = !missing(rejection_max_batches)
  )
  dots <- list(...)
  mixedgp_validate_named_dots(dots)
  target <- match.arg(target)
  engine <- mixedgp_fit_engine(object)
  latent_sampler <- match.arg(latent_sampler)
  if (target == "mean" && !is.null(new_U)) {
    stop("new_U is not used for target='mean'; omit it.")
  }
  if (target == "surface" && !is.null(new_C)) {
    stop("new_C is not used for target='surface'; omit it.")
  }
  if (is.null(new_U) && isTRUE(supplied$new_U_scale)) {
    stop("new_U_scale is relevant only when new_U is supplied.")
  }
  if (target == "response" &&
      isTRUE(supplied$include_process_uncertainty)) {
    stop(
      "include_process_uncertainty applies only to target='mean' or ",
      "target='surface'."
    )
  }
  if (target != "mean" && isTRUE(supplied$n_latent)) {
    stop("n_latent applies only to target='mean'.")
  }
  sampler_controls <- c(
    "latent_sampler", "n_new_latent_gibbs", "rejection_batch_size",
    "rejection_max_batches"
  )
  if ((engine == "univariate" || target == "surface") &&
      any(vapply(supplied[sampler_controls], isTRUE, logical(1L)))) {
    stop(
      "Prospective multivariate latent-sampler controls are unavailable for ",
      if (engine == "univariate") "the univariate engine." else
        "target='surface'."
    )
  }
  if (engine == "multivariate" && target != "surface") {
    if (isTRUE(supplied$n_new_latent_gibbs) && latent_sampler != "gibbs") {
      stop("n_new_latent_gibbs is relevant only for latent_sampler='gibbs'.")
    }
    if ((isTRUE(supplied$rejection_batch_size) ||
         isTRUE(supplied$rejection_max_batches)) &&
        latent_sampler != "rejection") {
      stop(
        "rejection_batch_size and rejection_max_batches are relevant only ",
        "for latent_sampler='rejection'."
      )
    }
  }
  joint <- mixedgp_validate_flag(joint, "joint")
  include_process_uncertainty <- mixedgp_validate_flag(
    include_process_uncertainty, "include_process_uncertainty"
  )
  n_per_draw <- mixedgp_as_integer_strict(
    n_per_draw, "n_per_draw", min_value = 1L, length_expected = 1L
  )
  if (!is.null(draw_ids)) {
    draw_ids <- mixedgp_as_integer_strict(
      draw_ids, "draw_ids", min_value = 1L
    )
  }
  n_latent <- mixedgp_as_integer_strict(
    n_latent, "n_latent", min_value = 1L, length_expected = 1L
  )
  n_new_latent_gibbs <- mixedgp_as_integer_strict(
    n_new_latent_gibbs, "n_new_latent_gibbs",
    min_value = 1L, length_expected = 1L
  )
  rejection_max_batches <- mixedgp_as_integer_strict(
    rejection_max_batches, "rejection_max_batches",
    min_value = 1L, length_expected = 1L
  )
  if (!is.null(rejection_batch_size)) {
    rejection_batch_size <- mixedgp_as_integer_strict(
      rejection_batch_size, "rejection_batch_size",
      min_value = 1L, length_expected = 1L
    )
  }
  if (!is.null(seed)) {
    seed <- mixedgp_as_integer_strict(
      seed, "seed", min_value = 0L, length_expected = 1L
    )
    caller_rng <- mixedgp_rng_state()
    on.exit(mixedgp_restore_rng_state(caller_rng), add = TRUE)
  }
  if (!is.null(new_U)) {
    if (engine == "multivariate" && is.null(dim(new_U)) &&
        length(new_U) == object$data$d) {
      new_U <- matrix(new_U, nrow = 1L)
    }
    new_U <- as_numeric_matrix_strict(
      new_U, "new_U",
      ncol_expected = if (engine == "univariate") 1L else object$data$d,
      allow_na = TRUE
    )
    if (engine == "univariate") new_U <- as.numeric(new_U[, 1L])
  }
  if (!is.null(seed)) set.seed(seed)
  latent_scale <- if (!is.null(new_U)) {
    mixedgp_resolve_latent_scale(object, new_U_scale, "new_U_scale")
  } else {
    list(internal = "model", label = "working", anchored = FALSE)
  }

  if (engine == "univariate") {
    if (is.null(draw_ids)) draw_ids <- seq_len(nrow(object$mcmc$samples_u))
    if (target == "response") {
      if (is.null(new_C)) {
        u_complete <- !is.null(new_U) && all(is.finite(as.numeric(new_U)))
        if (!u_complete) {
          stop("new_C is required whenever any new latent input is missing.")
        }
        encoded_new <- rep(1L, length(as.numeric(new_U)))
      } else {
        encoded_new <- prepare_ordinal_matrix(
          new_C,
          m_vec = object$data$m,
          level_maps = object$interface$ordinal_level_maps,
          expected_names = object$data$C_names,
          name = "new_C"
        )$C[, 1L]
      }
      ans <- mixedgp_do_call(
        sample_eiv_test_y,
        list(
          x_test_raw = new_X, c_test = encoded_new, fit_obj = object,
          draw_ids = draw_ids, n_per_draw = n_per_draw,
          u_test_obs = new_U, u_input_scale = latent_scale$internal,
          joint = isTRUE(joint)
        ),
        dots
      )
      if (!is.null(new_U)) {
        attr(ans, "latent_input_scale") <- latent_scale$label
        attr(ans, "latent_scale_anchored") <- latent_scale$anchored
      }
      return(ans)
    }
    if (target == "mean") {
      if (is.null(new_C)) stop("new_C is required for inference on m(x,c).")
      encoded_new <- prepare_ordinal_matrix(
        new_C,
        m_vec = object$data$m,
        level_maps = object$interface$ordinal_level_maps,
        expected_names = object$data$C_names,
        name = "new_C"
      )$C[, 1L]
      return(mixedgp_do_call(
        sample_eiv_m_given_xc_1d,
        list(
          x_star_raw = new_X, c_star = encoded_new, fit_obj = object,
          draw_ids = draw_ids, n_per_draw = n_per_draw,
          n_latent = n_latent,
          include_process_uncertainty = isTRUE(include_process_uncertainty),
          joint = isTRUE(joint), seed = seed
        ),
        dots
      ))
    }
    if (is.null(new_U)) stop("new_U is required for surface inference.")
    ans <- mixedgp_do_call(
      sample_eiv_f_given_xu_1d,
      list(
        x_star_raw = new_X, u_star = as.numeric(new_U), fit_obj = object,
        draw_ids = draw_ids, n_per_draw = n_per_draw,
        include_gp_uncertainty = isTRUE(include_process_uncertainty),
        u_input_scale = latent_scale$internal,
        joint = isTRUE(joint)
      ),
      dots
    )
    attr(ans, "latent_input_scale") <- latent_scale$label
    attr(ans, "latent_scale_anchored") <- latent_scale$anchored
    return(ans)
  }

  if (target == "response") {
    if (is.null(new_C)) {
      U_check <- if (is.null(new_U)) NULL else as.matrix(new_U)
      u_complete <- !is.null(U_check) && ncol(U_check) == object$data$d &&
        all(is.finite(U_check))
      if (!u_complete) {
        stop("new_C is required whenever any new latent input is missing.")
      }
      new_C <- as.data.frame(
        lapply(
          object$interface$ordinal_level_maps,
          function(levels_j) rep(levels_j[1L], nrow(U_check))
        ),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      names(new_C) <- object$data$C_names
    }
    ans <- mixedgp_do_call(
      sample_eiv_test_y_ordprobit_fb,
      list(
        X_test_raw = new_X, C_test = new_C, fit_obj = object,
        draw_ids = draw_ids, n_per_draw = n_per_draw,
        latent_sampler = latent_sampler,
        n_new_latent_gibbs = n_new_latent_gibbs,
        rejection_batch_size = rejection_batch_size,
        rejection_max_batches = rejection_max_batches,
        U_test_obs = new_U, U_input_scale = latent_scale$internal,
        joint = isTRUE(joint), seed = seed
      ),
      dots
    )
    if (!is.null(new_U)) {
      attr(ans, "latent_input_scale") <- latent_scale$label
      attr(ans, "latent_scale_anchored") <- latent_scale$anchored
    }
    return(ans)
  }
  if (target == "mean") {
    if (is.null(new_C)) stop("new_C is required for inference on m(x,c).")
    return(mixedgp_do_call(
      sample_eiv_m_given_xc_fb,
      list(
        X_test_raw = new_X, C_test = new_C, fit_obj = object,
        draw_ids = draw_ids, n_per_draw = n_per_draw,
        n_latent = n_latent,
        include_process_uncertainty = isTRUE(include_process_uncertainty),
        joint = isTRUE(joint), latent_sampler = latent_sampler,
        n_new_latent_gibbs = n_new_latent_gibbs,
        rejection_batch_size = rejection_batch_size,
        rejection_max_batches = rejection_max_batches,
        seed = seed
      ),
      dots
    ))
  }
  if (is.null(new_U)) stop("new_U is required for surface inference.")
  ans <- mixedgp_do_call(
    sample_eiv_f_given_xu_fb,
    list(
      X_test_raw = new_X, U_test_raw = new_U, fit_obj = object,
      draw_ids = draw_ids, n_per_draw = n_per_draw,
      include_process_uncertainty = isTRUE(include_process_uncertainty),
      U_input_scale = latent_scale$internal,
      joint = isTRUE(joint), seed = seed
    ),
    dots
  )
  attr(ans, "latent_input_scale") <- latent_scale$label
  attr(ans, "latent_scale_anchored") <- latent_scale$anchored
  ans
}

#' Posterior inference for latent continuous inputs
#'
#' Returns posterior draws of the latent input for either training rows or
#' prospective ordinal observations. The two cases are mutually exclusive:
#' omit `new_C` to recover training-row draws, or supply `new_C` for new units.
#'
#' @param object An object returned by [fit_eivgp()].
#' @param rows Optional training-row indices. Ignored only when `new_C` is
#'   supplied, in which case it must be `NULL`.
#' @param new_C Optional prospective ordinal inputs.
#' @param draw_ids Saved posterior draws to use; `NULL` uses every draw.
#' @param n_per_draw Number of prospective latent draws per posterior draw.
#'   Training draws already contain one latent state per posterior draw and
#'   therefore require `n_per_draw=1`.
#' @param scale Latent-coordinate scale. `"auto"` returns the physical raw
#'   scale when calibration anchors it and the identified working/model scale
#'   otherwise. `"raw"` is unavailable without a full-affine-rank calibration
#'   design.
#' @param seed Optional simulation seed for prospective units.
#' @param latent_sampler Prospective sampler for the multivariate ordinal-probit
#'   engine. The default is a minimax-tilted accept--reject sampler that is exact
#'   on successful completion. Exact prior rejection is an alternative, and
#'   finite Gibbs is a diagnostic approximation.
#' @param n_new_latent_gibbs Number of sweeps for the optional Gibbs diagnostic.
#' @param rejection_batch_size,rejection_max_batches Exact-rejection controls.
#' @param ... Controls for prospective multivariate latent simulation.
#' @return A numeric array with dimensions posterior draw by unit by latent
#'   coordinate. Training rows retain their requested row labels; prospective
#'   units follow the row order of `new_C`. Attributes `source`, `scale`, and
#'   `latent_scale_anchored` describe the returned coordinates.
#' @export
impute_eivgp <- function(object,
                         rows = NULL,
                         new_C = NULL,
                         draw_ids = NULL,
                         n_per_draw = 1L,
                         scale = c("auto", "raw", "model"),
                         latent_sampler = c(
                           "minimax_tilting", "rejection", "gibbs"
                         ),
                         n_new_latent_gibbs = 100L,
                         rejection_batch_size = NULL,
                         rejection_max_batches = 1000L,
                         seed = NULL,
                         ...) {
  if (!inherits(object, "eivgp_fit")) {
    stop("object must be returned by fit_eivgp().")
  }
  supplied <- list(
    latent_sampler = !missing(latent_sampler),
    n_new_latent_gibbs = !missing(n_new_latent_gibbs),
    rejection_batch_size = !missing(rejection_batch_size),
    rejection_max_batches = !missing(rejection_max_batches),
    seed = !missing(seed)
  )
  dots <- list(...)
  mixedgp_validate_named_dots(dots)
  engine <- mixedgp_fit_engine(object)
  prospective <- !is.null(new_C)
  latent_sampler <- match.arg(latent_sampler)
  sampler_controls <- c(
    "latent_sampler", "n_new_latent_gibbs", "rejection_batch_size",
    "rejection_max_batches"
  )
  if (!prospective &&
      (isTRUE(supplied$seed) || length(dots) > 0L ||
       any(vapply(supplied[sampler_controls], isTRUE, logical(1L))))) {
    stop(
      "seed, ..., and prospective latent-sampler controls are not used when ",
      "recovering stored training latent draws."
    )
  }
  if (prospective && engine == "univariate" &&
      any(vapply(supplied[sampler_controls], isTRUE, logical(1L)))) {
    stop(
      "Multivariate latent-sampler controls are unavailable for the ",
      "univariate engine."
    )
  }
  if (prospective && engine == "multivariate") {
    if (isTRUE(supplied$n_new_latent_gibbs) && latent_sampler != "gibbs") {
      stop("n_new_latent_gibbs is relevant only for latent_sampler='gibbs'.")
    }
    if ((isTRUE(supplied$rejection_batch_size) ||
         isTRUE(supplied$rejection_max_batches)) &&
        latent_sampler != "rejection") {
      stop(
        "rejection_batch_size and rejection_max_batches are relevant only ",
        "for latent_sampler='rejection'."
      )
    }
  }
  n_per_draw <- mixedgp_as_integer_strict(
    n_per_draw, "n_per_draw", min_value = 1L, length_expected = 1L
  )
  if (!is.null(draw_ids)) {
    draw_ids <- mixedgp_as_integer_strict(
      draw_ids, "draw_ids", min_value = 1L
    )
  }
  n_new_latent_gibbs <- mixedgp_as_integer_strict(
    n_new_latent_gibbs, "n_new_latent_gibbs",
    min_value = 1L, length_expected = 1L
  )
  rejection_max_batches <- mixedgp_as_integer_strict(
    rejection_max_batches, "rejection_max_batches",
    min_value = 1L, length_expected = 1L
  )
  if (!is.null(rejection_batch_size)) {
    rejection_batch_size <- mixedgp_as_integer_strict(
      rejection_batch_size, "rejection_batch_size",
      min_value = 1L, length_expected = 1L
    )
  }
  if (!is.null(seed)) {
    seed <- mixedgp_as_integer_strict(
      seed, "seed", min_value = 0L, length_expected = 1L
    )
    caller_rng <- mixedgp_rng_state()
    on.exit(mixedgp_restore_rng_state(caller_rng), add = TRUE)
  }
  scale_info <- mixedgp_resolve_latent_scale(object, scale, "scale")
  scale <- scale_info$internal
  if (prospective && !is.null(rows)) {
    stop("rows and new_C are mutually exclusive.")
  }
  if (!prospective && !is.null(rows)) {
    rows <- mixedgp_as_integer_strict(rows, "rows", min_value = 1L)
  }

  if (prospective) {
    if (engine == "univariate") {
      encoded <- prepare_ordinal_matrix(
        new_C,
        m_vec = object$data$m,
        level_maps = object$interface$ordinal_level_maps,
        expected_names = object$data$C_names,
        name = "new_C"
      )$C
      ans <- mixedgp_do_call(
        sample_eiv_u_given_c_1d,
        list(
          c_new = encoded[, 1L], fit_obj = object, draw_ids = draw_ids,
          n_per_draw = n_per_draw, scale = scale, seed = seed
        ),
        dots
      )
      attr(ans, "scale") <- scale_info$label
      attr(ans, "latent_scale_anchored") <- scale_info$anchored
      return(ans)
    }
    ans <- mixedgp_do_call(
      sample_eiv_u_given_c_ordprobit,
      list(
        C_new = new_C, fit_obj = object, draw_ids = draw_ids,
        n_per_draw = n_per_draw, scale = scale,
        latent_sampler = latent_sampler,
        n_new_latent_gibbs = n_new_latent_gibbs,
        rejection_batch_size = rejection_batch_size,
        rejection_max_batches = rejection_max_batches,
        seed = seed
      ),
      dots
    )
    attr(ans, "scale") <- scale_info$label
    attr(ans, "latent_scale_anchored") <- scale_info$anchored
    return(ans)
  }

  if (n_per_draw != 1L) {
    stop("Training latent draws require n_per_draw = 1.")
  }
  if (engine == "univariate") {
    n_saved <- nrow(object$mcmc$samples_u)
    n_train <- ncol(object$mcmc$samples_u)
    if (is.null(draw_ids)) draw_ids <- seq_len(n_saved)
    if (is.null(rows)) rows <- seq_len(n_train)
    draw_ids <- mixedgp_as_integer_strict(
      draw_ids, "draw_ids", min_value = 1L
    )
    rows <- mixedgp_as_integer_strict(rows, "rows", min_value = 1L)
    if (length(draw_ids) < 1L || anyNA(draw_ids) ||
        any(draw_ids < 1L | draw_ids > n_saved) ||
        length(rows) < 1L || anyNA(rows) || any(rows < 1L | rows > n_train)) {
      stop("draw_ids or rows contains invalid indices.")
    }
    values <- object$mcmc$samples_u[draw_ids, rows, drop = FALSE]
    out <- array(
      values,
      dim = c(length(draw_ids), length(rows), 1L),
      dimnames = list(NULL, as.character(rows), object$data$U_names)
    )
    if (scale == "raw") {
      center <- if (is.null(object$data$U_center)) 0 else object$data$U_center
      scale_value <- if (is.null(object$data$U_scale)) 1 else object$data$U_scale
      out[, , 1L] <- center + scale_value * out[, , 1L]
    }
  } else {
    all_draws <- posterior_u_draws(object, rows = rows, scale = scale)
    n_saved <- dim(all_draws)[1L]
    if (is.null(draw_ids)) draw_ids <- seq_len(n_saved)
    draw_ids <- mixedgp_as_integer_strict(
      draw_ids, "draw_ids", min_value = 1L
    )
    if (length(draw_ids) < 1L || anyNA(draw_ids) ||
        any(draw_ids < 1L | draw_ids > n_saved)) {
      stop("draw_ids contains invalid posterior-draw indices.")
    }
    out <- all_draws[draw_ids, , , drop = FALSE]
  }
  attr(out, "source") <- "training"
  attr(out, "scale") <- scale_info$label
  attr(out, "latent_scale_anchored") <- scale_info$anchored
  out
}

#' Fit the unrestricted-correlation GP competitor
#' @inheritParams fit_lvgp
#' @export
fit_ucgp <- function(X,
                     y,
                     C,
                     new_X,
                     new_C,
                     m_vec = NULL,
                     n_draw = 1000L,
                     seed = 1L,
                     ...) {
  n_draw <- mixedgp_as_integer_strict(
    n_draw, "n_draw", min_value = 1L, length_expected = 1L
  )
  seed <- mixedgp_as_integer_strict(
    seed, "seed", min_value = 0L, length_expected = 1L
  )
  caller_rng <- mixedgp_rng_state()
  on.exit(mixedgp_restore_rng_state(caller_rng), add = TRUE)
  mixedgp_do_call(
    mixedgp_adapter_ucgp,
    list(
      X_train = X, y_train = y, C_train = C,
      X_test = new_X, C_test = new_C, m_vec = m_vec,
      n_draw = n_draw, seed = seed
    ),
    list(...)
  )
}

#' Fit the published latent-variable GP competitor
#'
#' @param X Numeric training design matrix.
#' @param y Numeric response vector.
#' @param C Integer-coded ordinal training inputs.
#' @param new_X Numeric test design matrix.
#' @param new_C Integer-coded ordinal test inputs.
#' @param m_vec Number of levels for each ordinal input.
#' @param n_draw Number of noisy-response predictive draws.
#' @param seed Reproducible fitting seed.
#' @param ... Method-specific controls passed to the audited wrapper.
#' @return A list containing the public-package fit and predictive draws.
#' @export
fit_lvgp <- function(X,
                     y,
                     C,
                     new_X,
                     new_C,
                     m_vec = NULL,
                     n_draw = 1000L,
                     seed = 1L,
                     ...) {
  n_draw <- mixedgp_as_integer_strict(
    n_draw, "n_draw", min_value = 1L, length_expected = 1L
  )
  seed <- mixedgp_as_integer_strict(
    seed, "seed", min_value = 0L, length_expected = 1L
  )
  caller_rng <- mixedgp_rng_state()
  on.exit(mixedgp_restore_rng_state(caller_rng), add = TRUE)
  mixedgp_do_call(
    mixedgp_adapter_lvgp,
    list(
      X_train = X, y_train = y, C_train = C,
      X_test = new_X, C_test = new_C, m_vec = m_vec,
      n_draw = n_draw, seed = seed
    ),
    list(...)
  )
}

#' Fit the published easy-to-interpret GP competitor
#' @inheritParams fit_lvgp
#' @export
fit_ezgp <- function(X,
                     y,
                     C,
                     new_X,
                     new_C,
                     m_vec = NULL,
                     n_draw = 1000L,
                     seed = 1L,
                     ...) {
  n_draw <- mixedgp_as_integer_strict(
    n_draw, "n_draw", min_value = 1L, length_expected = 1L
  )
  seed <- mixedgp_as_integer_strict(
    seed, "seed", min_value = 0L, length_expected = 1L
  )
  caller_rng <- mixedgp_rng_state()
  on.exit(mixedgp_restore_rng_state(caller_rng), add = TRUE)
  mixedgp_do_call(
    mixedgp_adapter_ezgp,
    list(
      X_train = X, y_train = y, C_train = C,
      X_test = new_X, C_test = new_C, m_vec = m_vec,
      n_draw = n_draw, seed = seed
    ),
    list(...)
  )
}

#' Fit one audited published mixed-input GP competitor
#' @param method One of `"UC-GP"`, `"LVGP"`, or `"EzGP"`.
#' @param ... Arguments passed to the method-specific single-call wrapper.
#' @export
fit_mixedgp_competitor <- function(method, ...) {
  method <- match.arg(method, c("UC-GP", "LVGP", "EzGP"))
  dots <- list(...)
  mixedgp_validate_named_dots(dots)
  FUN <- switch(method, `UC-GP` = fit_ucgp, LVGP = fit_lvgp,
                EzGP = fit_ezgp)
  do.call(FUN, dots)
}
