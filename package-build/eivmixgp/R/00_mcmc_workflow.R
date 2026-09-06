## Public, explicit diagnostic-and-extension workflow. No automatic stopping
## or convergence claim: finite-chain diagnostics can miss common-mode failures.

#' Extend the existing posterior chains
#'
#' @param object An unmodified fit returned by [fit_eivgp()] in version 0.2.0
#'   or later, or by a previous call to this function.
#' @param n_iter Additional transitions per chain, not additional saved draws.
#'   Every post-warm-up transition is retained; there is no thinning.
#' @param parallel Whether to run chains in parallel. NULL keeps the fit setting.
#' @param n_cores Worker count. NULL keeps the fit setting.
#' @param verbose Whether to print sampler progress.
#' @return An updated `eivgp_fit`, containing both old and new retained draws.
#' @details The terminal state, random-number stream, original thinning phase,
#'   update schedule, and frozen post-warm-up tuning are preserved. No new
#'   warm-up is discarded and no model or sampler options can be changed here.
#'   Checkpoints survive `saveRDS()`/`readRDS()`; use the same package version
#'   and numerical environment for reproducible continuation. Older fits lacking
#'   checkpoints must be refitted. Diagnostics in the returned fit are recomputed
#'   over all retained draws; call [diagnose_eivgp()] again for target diagnostics.
#'   Running longer is not a remedy for unidentified coordinates or chains
#'   trapped in different modes. This is continuation of a completed call, not
#'   periodic crash recovery during an interrupted fit.
#' @export
continue_eivgp <- function(object, n_iter, parallel = NULL, n_cores = NULL,
                           verbose = FALSE) {
  if (!inherits(object, "eivgp_fit")) stop("object must be an eivgp_fit.")
  caller_rng <- mixedgp_rng_state()
  on.exit(mixedgp_restore_rng_state(caller_rng), add = TRUE)
  n_iter <- mixedgp_as_integer_strict(n_iter, "n_iter", 1L, 1L)
  verbose <- mixedgp_validate_flag(verbose, "verbose")
  cp <- object$checkpoint
  if (is.null(cp) || !identical(cp$version, 1L) || is.null(cp$arguments)) {
    stop("This fit has no supported continuation checkpoint; refit with version 0.2.0 or later.")
  }
  if (!identical(cp$control, object$control)) stop("The fit controls were modified; cannot continue.")
  if (!identical(cp$signature, mixedgp_checkpoint_hash(object))) {
    stop("The fit data, draws, or checkpoint were modified; cannot safely continue.")
  }
  if (cp$thin != 1L) stop("Thinned historical fits cannot recover discarded draws; refit without thinning.")
  total <- mixedgp_as_integer_strict(as.double(cp$iteration) + n_iter,
                                   "total iterations", 2L, 1L)
  expected <- floor((cp$iteration - cp$burn) / cp$thin)
  if (floor((total - cp$burn) / cp$thin) <= expected) {
    stop("The extension must reach at least one new retained draw per chain.")
  }
  by_chain <- object$mcmc$samples_by_chain
  n_chains <- length(cp$states)
  if (!n_chains || length(by_chain$logtheta) != n_chains ||
      any(vapply(by_chain$logtheta, nrow, integer(1)) != expected) ||
      nrow(object$mcmc$chain_stats) != n_chains) {
    stop("The retained chain history does not match its checkpoint.")
  }
  cp$old_chains <- lapply(seq_len(n_chains), function(i) {
    out <- lapply(by_chain, function(chains) chains[[i]])
    names(out) <- paste0("samples_", names(out))
    out$stats <- object$mcmc$chain_stats[i, , drop = FALSE]
    if (!is.null(object$mcmc$chain_initial)) out$initial_state <- object$mcmc$chain_initial[[i]]
    out
  })
  if (is.null(parallel)) parallel <- object$interface$parallel
  parallel <- mixedgp_validate_flag(parallel, "parallel")
  if (is.null(n_cores)) n_cores <- object$interface$n_cores
  if (!is.null(n_cores)) n_cores <- mixedgp_as_integer_strict(n_cores, "n_cores", 1L, 1L)
  arguments <- cp$arguments
  arguments$n_iter <- total
  arguments$parallel_chains <- parallel
  arguments$n_cores <- n_cores
  arguments$verbose <- verbose
  engine <- mixedgp_fit_engine(object)
  sampler <- if (engine == "univariate") fit_eivgp_1d else fit_eivgp_ordprobit_fb
  updated <- do.call(sampler, c(arguments, list(.resume = cp)))
  updated$checkpoint$arguments <- arguments
  updated$data <- object$data
  updated$call <- object$call
  updated$model_specification <- object$model_specification
  updated$interface <- object$interface
  updated$interface$parallel <- parallel
  updated$interface$n_cores <- updated$diagnostics$summary$parallel_cores
  prior_time <- object$diagnostics$summary$time_seconds
  updated$diagnostics$summary$time_seconds <- prior_time + updated$diagnostics$summary$time_seconds
  updated$continuation_history <- rbind(object$continuation_history, data.frame(
    from_iteration = cp$iteration, to_iteration = total,
    added_draws_per_chain = floor((total - cp$burn) / cp$thin) - expected
  ))
  class(updated) <- class(object)
  updated$checkpoint$signature <- mixedgp_checkpoint_hash(updated)
  updated
}

mixedgp_diagnostic_status <- function(tab, rhat_limit, bulk_ess_limit,
                                      tail_ess_limit, mcse_sd_ratio_limit) {
  valid <- is.finite(tab$rhat) & is.finite(tab$ess_bulk) & is.finite(tab$ess_tail) &
    is.finite(tab$mcse_sd_ratio) & tab$n_chains >= 4L & tab$draw_window_complete
  status <- rep("insufficient_diagnostics", nrow(tab))
  status[valid] <- "screen_passed"
  precision <- tab$ess_bulk < bulk_ess_limit | tab$ess_tail < tail_ess_limit |
    tab$mcse_sd_ratio > mcse_sd_ratio_limit
  status[valid & precision] <- "insufficient_precision"
  status[valid & tab$rhat > rhat_limit] <- "poor_exploration"
  status
}

#' Diagnose posterior exploration and Monte Carlo precision
#'
#' @param object An `eivgp_fit`.
#' @param X,C Diagnostic panel in the same units/ordinal labels as fitting data.
#'   Supply both or neither. By default at most five training-input locations
#'   are chosen deterministically; these are locations for evaluating m(x,c),
#'   not a claim to cover the response surface.
#' @param U Optional exact latent inputs at the panel locations, on the original
#'   calibration scale, to additionally monitor f(x,u). Requires the package's
#'   calibrated-scale reporting gate to be satisfied.
#' @param n_latent Integration points per location and posterior state.
#' @param seed Seed for common integration randomness, not for the fitted chains.
#' @param rhat_limit Maximum rank-normalized split R-hat (default 1.01).
#' @param bulk_ess_limit,tail_ess_limit Minimum bulk/tail ESS (default 400 each).
#' @param mcse_sd_ratio_limit Maximum MCSE of each monitored mean relative to
#'   its posterior standard deviation (default 0.05). This is a screening
#'   tolerance, not a decision-specific absolute accuracy guarantee.
#' @param additional_series Optional uniquely named list of scientific
#'   functionals. Each entry must be a list of one numeric vector per chain,
#'   aligned with every retained draw in that chain. For example supply a
#'   scientifically important contrast. No thinning or truncation is allowed.
#' @return An `eivgp_diagnostics` with per-functional `table`, overall `status`,
#'   `recommendation`, diagnostic `panel`, and `settings`.
#' @details All retained draws are checked, using at least four chains. Reports
#'   include free parameters, missing training U coordinates and their squares,
#'   multivariate measurement invariants, and panel-specific conditional means,
#'   conditional GP variances, and predictive second moments. These quantities
#'   avoid adding independent predictive noise that can disguise stuck chains.
#'   Mean MCSE is computed with `posterior::mcse_mean`, not bulk ESS.
#'   Common finite integration randomness makes target signatures deterministic
#'   functions of chain states, but integration error is not included in MCSE:
#'   repeat with larger `n_latent` and another integration seed for sensitivity.
#'   Quantile MCSE, simultaneous bands, and decision losses require separate
#'   checks. Passing this screen proves neither convergence, identification,
#'   nor model adequacy. Uncalibrated latent-coordinate failures may reflect
#'   symmetries; the raw and invariant panels are reported separately, without
#'   silently dismissing an unsuccessful raw-coordinate check.
#' @export
diagnose_eivgp <- function(object, X = NULL, C = NULL, U = NULL,
                           n_latent = 64L, seed = 481517L,
                           rhat_limit = 1.01, bulk_ess_limit = 400,
                           tail_ess_limit = 400, mcse_sd_ratio_limit = 0.05,
                           additional_series = list()) {
  if (!inherits(object, "eivgp_fit")) stop("object must be an eivgp_fit.")
  for (name in c("rhat_limit", "bulk_ess_limit", "tail_ess_limit", "mcse_sd_ratio_limit")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0) {
      stop(name, " must be positive and finite.")
    }
  }
  if (rhat_limit < 1) stop("rhat_limit must be at least 1.")
  n_latent <- mixedgp_as_integer_strict(n_latent, "n_latent", 2L, 1L)
  seed <- mixedgp_as_integer_strict(seed, "seed", 0L, 1L)
  engine <- mixedgp_fit_engine(object)
  if (xor(is.null(X), is.null(C))) stop("Supply both X and C, or neither.")
  if (is.null(X)) {
    X <- if (engine == "univariate") object$data$x_raw else object$data$X_raw
    C <- if (engine == "univariate") matrix(object$data$c_ord, ncol = 1L) else object$data$C_ord
    rows <- mixedgp_study2_diagnostic_rows(X, C)
    X <- X[rows, , drop = FALSE]
    C <- C[rows, , drop = FALSE]
    if (!is.null(U)) stop("Supply explicit X and C when supplying U.")
  } else {
    X <- as_numeric_matrix_strict(X, "X")
    encoded <- mixedgp_encode_ordinal(C,
      if (engine == "univariate") object$data$m else object$data$m_vec,
      object$data$C_level_maps)
    C <- encoded$C
  }
  expected_p <- if (engine == "univariate") ncol(object$data$x) else object$data$p
  expected_q <- if (engine == "univariate") 1L else object$data$q
  if (ncol(X) != expected_p || nrow(X) != nrow(C) || ncol(C) != expected_q ||
      nrow(X) < 1L || any(!is.finite(X))) stop("Invalid diagnostic panel dimensions or values.")
  if (!is.null(U)) {
    if (!isTRUE(object$data$latent_scale_anchored)) {
      stop("Physical-scale surface diagnostics require adequate calibration anchoring.")
    }
    U <- as_numeric_matrix_strict(U, "U", nrow_expected = nrow(X))
    if (ncol(U) != if (engine == "univariate") 1L else object$data$d) {
      stop("U has the wrong latent dimension.")
    }
  }
  raw <- if (engine == "univariate") mixedgp_study1_raw_series(object) else mixedgp_study2_raw_series(object)
  latent_names <- names(raw)[grepl("^[uU]\\[", names(raw))]
  for (name in latent_names) raw[[paste0(name, "^2")]] <- lapply(raw[[name]], function(z) z^2)
  invariant <- if (engine == "multivariate") mixedgp_study2_invariant_series(object) else list()
  targets <- if (engine == "univariate") {
    mixedgp_study1_target_series(object, X, C, U, n_latent = n_latent, seed = seed)
  } else {
    mixedgp_study2_target_series(object, X, C, U, n_latent = n_latent, seed = seed)
  }
  series <- c(raw, invariant, targets)
  mixedgp_validate_named_dots(additional_series, "additional_series")
  if (length(intersect(names(series), names(additional_series)))) {
    stop("additional_series names must not duplicate built-in functionals.")
  }
  lens <- vapply(object$mcmc$samples_by_chain$logtheta, nrow, integer(1))
  for (name in names(additional_series)) {
    z <- additional_series[[name]]
    if (!is.list(z) || length(z) != length(lens) ||
        !identical(as.integer(lengths(z)), as.integer(lens)) ||
        !all(vapply(z, function(v) is.numeric(v) && is.null(dim(v)), logical(1)))) {
      stop("Each additional functional must match all retained draws in every chain.")
    }
  }
  series <- c(series, additional_series)
  tab <- mixedgp_summarize_diagnostic_series(series, rhat_limit, bulk_ess_limit,
    tail_ess_limit, mcse_sd_ratio_limit = mcse_sd_ratio_limit)
  tab$group <- rep(c("raw_and_training_imputation", "measurement_invariant",
                     "scientific_panel", "user_functional"),
                   c(length(raw), length(invariant), length(targets), length(additional_series)))
  tab$status <- mixedgp_diagnostic_status(tab, rhat_limit, bulk_ess_limit,
                                        tail_ess_limit, mcse_sd_ratio_limit)
  priority <- c("poor_exploration", "insufficient_diagnostics", "insufficient_precision", "screen_passed")
  status <- priority[priority %in% tab$status][1L]
  recommendation <- switch(status,
    poor_exploration = paste("Inspect trace plots, chain-specific targets, calibration and latent symmetries.",
      "Do not assume extending alone will resolve poor exploration."),
    insufficient_diagnostics = paste("Use at least four chains and more retained draws; inspect nonfinite",
      "or constant functionals. A diagnostic pass cannot be issued."),
    insufficient_precision = paste("Chain agreement meets the R-hat screen but precision is insufficient.",
      "Consider continue_eivgp(), then recompute this report on all draws."),
    screen_passed = paste("The monitored functionals pass the specified screen. Check integration sensitivity,",
      "decision-specific MC error and model adequacy before scientific reporting."))
  structure(list(status = status, recommendation = recommendation, table = tab,
    panel = list(X = X, C = C, U = U),
    settings = list(rhat_limit = rhat_limit, bulk_ess_limit = bulk_ess_limit,
      tail_ess_limit = tail_ess_limit, mcse_sd_ratio_limit = mcse_sd_ratio_limit,
      n_latent = n_latent, seed = seed, full_retained_window = TRUE,
      physical_scale_anchored = isTRUE(object$data$latent_scale_anchored))),
    class = "eivgp_diagnostics")
}

#' @export
print.eivgp_diagnostics <- function(x, ...) {
  cat("EIV-GP diagnostic screen:", x$status, "\n")
  print(table(x$table$group, x$table$status))
  cat(x$recommendation, "\n")
  invisible(x)
}
