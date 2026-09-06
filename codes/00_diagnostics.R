## Release diagnostics must retain nonfinite entries: dropping one changes the
## proposition being checked. Only known fixed model coordinates are excluded
## when constructing the series, never because their diagnostics are inconvenient.
mixedgp_diagnostic_table_pass <- function(x, rhat_limit = 1.01,
                                           bulk_ess_limit = 400,
                                           tail_ess_limit = 400,
                                           expected_parameters = NULL,
                                           min_chains = 4L,
                                           mcse_limit = NULL,
                                           mcse_sd_ratio_limit = NULL) {
  required <- c("parameter", "rhat", "ess_bulk", "ess_tail", "n_chains")
  if (!is.data.frame(x) || nrow(x) == 0L ||
      !all(required %in% names(x)) || anyNA(x$parameter) ||
      anyDuplicated(x$parameter)) return(FALSE)
  if (!is.null(expected_parameters) &&
      !setequal(x$parameter, expected_parameters)) return(FALSE)
  if ("draw_window_complete" %in% names(x) &&
      !all(!is.na(x$draw_window_complete) & x$draw_window_complete)) return(FALSE)
  for (item in list(c("mcse", "mcse_limit"),
                    c("mcse_sd_ratio", "mcse_sd_ratio_limit"))) {
    limit <- get(item[2L])
    if (!is.null(limit)) {
      if (!is.numeric(limit) || length(limit) != 1L ||
          !is.finite(limit) || limit <= 0) stop(item[2L], " must be positive and finite.")
      if (!item[1L] %in% names(x) || any(!is.finite(x[[item[1L]]])) ||
          any(x[[item[1L]]] > limit)) return(FALSE)
    }
  }
  values <- as.matrix(x[, c("rhat", "ess_bulk", "ess_tail", "n_chains")])
  all(is.finite(values)) && all(x$n_chains >= min_chains) &&
    all(x$rhat <= rhat_limit) && all(x$ess_bulk >= bulk_ess_limit) &&
    all(x$ess_tail >= tail_ess_limit)
}

mixedgp_summarize_diagnostic_series <- function(series, rhat_limit = 1.01,
                                                bulk_ess_limit = 400,
                                                tail_ess_limit = 400,
                                                mcse_limit = NULL,
                                                mcse_sd_ratio_limit = NULL) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Release diagnostics require the posterior package.")
  }
  if (length(series) == 0L) return(data.frame())
  do.call(rbind, lapply(names(series), function(parameter) {
    chains <- lapply(series[[parameter]], as.numeric)
    lens <- lengths(chains)
    values <- unlist(chains, use.names = FALSE)
    valid <- length(chains) >= 2L && length(unique(lens)) == 1L &&
      all(lens >= 4L) && all(is.finite(values))
    mat <- if (valid) do.call(cbind, chains) else NULL
    rhat <- if (valid) as.numeric(posterior::rhat(mat)) else NA_real_
    eb <- if (valid) as.numeric(posterior::ess_bulk(mat)) else NA_real_
    et <- if (valid) as.numeric(posterior::ess_tail(mat)) else NA_real_
    em <- if (valid) as.numeric(posterior::ess_mean(mat)) else NA_real_
    mcse <- if (valid) as.numeric(posterior::mcse_mean(mat)) else NA_real_
    row <- data.frame(
      parameter = parameter, n_chains = length(chains),
      min_draws_per_chain = if (length(lens)) min(lens) else 0L,
      max_draws_per_chain = if (length(lens)) max(lens) else 0L,
      mean = if (length(values)) mean(values) else NA_real_,
      posterior_sd = if (length(values) > 1L) stats::sd(values) else NA_real_,
      rhat = rhat, ess_bulk = eb, ess_tail = et, ess_mean = em,
      draw_window_complete = !identical(attr(series, "draw_window_complete"), FALSE) &&
        !identical(attr(series[[parameter]], "draw_window_complete"), FALSE),
      stringsAsFactors = FALSE
    )
    row$ess <- min(eb, et)
    # Rank-normalized bulk ESS is a mixing screen, not the ESS of a raw mean.
    # Each series is the actual functional being averaged (including second
    # moments and conditional variances), so use its mean-specific MCSE.
    row$mcse <- mcse
    row$mcse_sd_ratio <- if (is.finite(row$posterior_sd) && row$posterior_sd > 0)
      mcse / row$posterior_sd else NA_real_
    row$target_pass <- mixedgp_diagnostic_table_pass(
      row, rhat_limit, bulk_ess_limit, tail_ess_limit,
      mcse_limit = mcse_limit, mcse_sd_ratio_limit = mcse_sd_ratio_limit
    )
    row
  }))
}

## Inf means the entire reporting window. An explicitly truncated window is
## useful for pilot costs, but its diagnostics cannot certify all retained draws.
mixedgp_diagnostic_draw_indices <- function(n_saved, max_draws_per_chain = Inf) {
  if (!is.numeric(max_draws_per_chain) || length(max_draws_per_chain) != 1L ||
      is.na(max_draws_per_chain) || max_draws_per_chain < 1 ||
      (is.finite(max_draws_per_chain) &&
       max_draws_per_chain != floor(max_draws_per_chain))) {
    stop("max_draws_per_chain must be a positive integer or Inf (all draws).")
  }
  tail(seq_len(n_saved), min(n_saved, max_draws_per_chain))
}

mixedgp_mark_diagnostic_window <- function(series, complete) {
  attr(series, "draw_window_complete") <- isTRUE(complete)
  # Per-parameter provenance survives ordinary outer-list c() and [ subsetting.
  for (name in names(series)) attr(series[[name]], "draw_window_complete") <- isTRUE(complete)
  series
}

mixedgp_study2_invariant_series <- function(fit) {
  chains <- fit$mcmc$samples_by_chain
  q <- fit$data$q
  d <- fit$data$d
  series <- list()
  score_sd <- lapply(chains$A, function(a) {
    out <- matrix(NA_real_, dim(a)[1L], q)
    for (j in seq_len(q)) {
      aj <- matrix(a[, j, , drop = FALSE], ncol = d)
      out[, j] <- sqrt(1 + rowSums(aj^2))
    }
    out
  })
  if (q > 1L) for (pair in combn(seq_len(q), 2L, simplify = FALSE)) {
    j <- pair[1L]
    k <- pair[2L]
    name <- paste0("ordinal_score_cor[", j, ",", k, "]")
    series[[name]] <- lapply(seq_along(chains$A), function(id) {
      a <- chains$A[[id]]
      aj <- matrix(a[, j, , drop = FALSE], ncol = d)
      ak <- matrix(a[, k, , drop = FALSE], ncol = d)
      rowSums(aj * ak) / (score_sd[[id]][, j] * score_sd[[id]][, k])
    })
  }
  for (j in seq_len(q)) for (r in seq_len(fit$data$m_vec[j] - 1L)) {
    name <- paste0("tau[", j, ",", r, "]")
    column <- match(name, colnames(chains$tau[[1L]]))
    if (is.na(column)) stop("Missing diagnostic cutpoint: ", name)
    series[[paste0("standardized_", name)]] <- lapply(
      seq_along(chains$tau), function(id) {
        chains$tau[[id]][, column] / score_sd[[id]][, j]
      }
    )
  }
  series
}

mixedgp_study2_raw_series <- function(fit) {
  chains <- fit$mcmc$samples_by_chain
  series <- list(sigma_epsilon = lapply(chains$sigma2, sqrt))
  theta_names <- c("rho", paste0("theta_x", seq_len(fit$data$p)),
                   paste0("theta_u", seq_len(fit$data$d)))
  for (k in seq_along(theta_names)) {
    series[[theta_names[k]]] <- lapply(chains$logtheta, function(z) exp(z[, k]))
  }
  for (j in seq_len(fit$data$q)) for (k in seq_len(fit$data$d)) {
    fixed <- identical(fit$data$ident, "lower_triangular") &&
      j <= fit$data$d && k > j
    if (!fixed) series[[paste0("A[", j, ",", k, "]")]] <-
      lapply(chains$A, function(a) a[, j, k])
  }
  for (k in seq_len(ncol(chains$tau[[1L]]))) {
    series[[colnames(chains$tau[[1L]])[k]]] <-
      lapply(chains$tau, function(z) z[, k])
  }
  for (i in fit$data$miss_idx) for (k in seq_len(fit$data$d)) {
    series[[paste0("U[", i, ",", k, "]")]] <-
      lapply(chains$U, function(z) z[, i, k])
  }
  series
}

mixedgp_study1_raw_series <- function(fit) {
  chains <- fit$mcmc$samples_by_chain
  p <- ncol(fit$data$x)
  names_theta <- c("rho", if (p == 1L) "theta_x" else paste0("theta_x", seq_len(p)),
                   "theta_u")
  series <- list(sigma_epsilon = lapply(chains$sigma2, sqrt))
  for (k in seq_along(names_theta)) {
    series[[names_theta[k]]] <- lapply(chains$logtheta, function(z) exp(z[, k]))
  }
  for (k in seq_len(ncol(chains$tau[[1L]]))) {
    series[[paste0("tau", k)]] <- lapply(chains$tau, function(z) z[, k])
  }
  for (i in fit$data$miss_idx) {
    series[[paste0("u[", i, "]")]] <- lapply(chains$u, function(z) z[, i])
  }
  series
}

mixedgp_study1_raw_diagnostics <- function(fit, rhat_limit = 1.01,
                                            bulk_ess_limit = 400,
                                            tail_ess_limit = 400) {
  mixedgp_summarize_diagnostic_series(
    mixedgp_study1_raw_series(fit), rhat_limit, bulk_ess_limit, tail_ess_limit
  )
}

mixedgp_study1_target_series <- function(fit, X, C, U = NULL,
                                         max_draws_per_chain = Inf,
                                         n_latent = 64L, seed = 481517L) {
  caller_rng <- mixedgp_rng_state()
  on.exit(mixedgp_restore_rng_state(caller_rng), add = TRUE)
  X <- as.matrix(X)
  C <- as.integer(C)
  G <- nrow(X)
  M <- as.integer(n_latent)
  if (G < 1L || length(C) != G || anyNA(C) ||
      any(C < 1L | C > fit$data$m) || M < 2L) stop("Invalid Study I diagnostic panel.")
  X_std <- sweep(sweep(X, 2L, fit$data$x_center, "-"),
                 2L, fit$data$x_scale, "/")
  colnames(X_std) <- colnames(fit$data$x)
  calibrated <- length(fit$data$calib_idx) > 0L
  U_std <- if (!is.null(U)) {
    uc <- if (is.null(fit$data$U_center)) 0 else fit$data$U_center
    us <- if (is.null(fit$data$U_scale)) 1 else fit$data$U_scale
    (as.numeric(U) - uc) / us
  } else NULL
  if (!is.null(U_std) && (length(U_std) != G || any(!is.finite(U_std)))) {
    stop("U does not match the Study I diagnostic panel.")
  }
  nms <- c(paste0("m_conditional_mean[", seq_len(G), "]"),
           paste0("m_conditional_variance[", seq_len(G), "]"),
           paste0("predictive_second_moment[", seq_len(G), "]"))
  if (calibrated) nms <- c(nms,
    paste0("prospective_U_mean[", seq_len(G), "]"),
    paste0("prospective_U_second_moment[", seq_len(G), "]"))
  if (calibrated && !is.null(U_std)) nms <- c(nms,
    paste0("f_conditional_mean[", seq_len(G), "]"),
    paste0("f_conditional_variance[", seq_len(G), "]"))
  chains <- fit$mcmc$samples_by_chain
  series <- setNames(lapply(nms, function(z) vector("list", length(chains$u))), nms)
  attr(series, "draw_window_complete") <- TRUE
  kernel <- kernel_spec_from_fit_1d(fit)
  repeated <- rep(seq_len(G), each = M)
  for (id in seq_along(chains$u)) {
    n_saved <- nrow(chains$u[[id]])
    keep <- mixedgp_diagnostic_draw_indices(n_saved, max_draws_per_chain)
    if (length(keep) < n_saved) attr(series, "draw_window_complete") <- FALSE
    values <- matrix(NA_real_, length(keep), length(nms), dimnames = list(NULL, nms))
    for (ii in seq_along(keep)) {
      s <- keep[ii]
      tau <- chains$tau[[id]][s, ]
      # Common integration randomness avoids injecting fresh outcome/GP noise.
      # This finite-M approximation still requires an integration sensitivity check.
      set.seed(seed)
      u_mc <- matrix(rtruncnorm_vec(0, 1,
        c(-Inf, tau)[C[repeated]], c(tau, Inf)[C[repeated]]), M, G)
      args_gp <- list(x_train = fit$data$x, u_train = chains$u[[id]][s, ],
        y_train = fit$data$y, logtheta = chains$logtheta[[id]][s, ],
        sigma2_eps = chains$sigma2[[id]][s], kernel = kernel$name,
        matern_nu = kernel$matern_nu)
      mp <- do.call(gp_integrated_mean_state_1d,
                    c(args_gp, list(x_star = X_std, u_mc = u_mc)))
      values[ii, paste0("m_conditional_mean[", seq_len(G), "]")] <-
        fit$data$y_center + fit$data$y_scale * mp$mean
      values[ii, paste0("m_conditional_variance[", seq_len(G), "]")] <-
        fit$data$y_scale^2 * mp$var
      yp <- do.call(gp_predict_draw, c(args_gp, list(
        x_star = X_std[repeated, , drop = FALSE], u_star = as.numeric(u_mc),
        noisy = TRUE)))
      mu <- fit$data$y_center + fit$data$y_scale * yp$mean
      values[ii, paste0("predictive_second_moment[", seq_len(G), "]")] <-
        colMeans(matrix(fit$data$y_scale^2 * yp$var + mu^2, M, G))
      if (calibrated) {
        uc <- if (is.null(fit$data$U_center)) 0 else fit$data$U_center
        us <- if (is.null(fit$data$U_scale)) 1 else fit$data$U_scale
        u_physical <- uc + us * u_mc
        values[ii, paste0("prospective_U_mean[", seq_len(G), "]")] <- colMeans(u_physical)
        values[ii, paste0("prospective_U_second_moment[", seq_len(G), "]")] <- colMeans(u_physical^2)
      }
      if (calibrated && !is.null(U_std)) {
        fp <- do.call(gp_predict_draw, c(args_gp, list(x_star = X_std, u_star = U_std)))
        values[ii, paste0("f_conditional_mean[", seq_len(G), "]")] <-
          fit$data$y_center + fit$data$y_scale * fp$mean
        values[ii, paste0("f_conditional_variance[", seq_len(G), "]")] <-
          fit$data$y_scale^2 * fp$var
      }
    }
    for (name in nms) series[[name]][[id]] <- values[, name]
  }
  mixedgp_mark_diagnostic_window(series, attr(series, "draw_window_complete"))
}

mixedgp_study2_diagnostic_rows <- function(X, C, max_points = 5L) {
  X <- as.matrix(X)
  C <- as.matrix(C)
  keys <- apply(C, 1L, paste, collapse = ":")
  counts <- sort(table(keys), decreasing = TRUE)
  patterns <- names(counts)[unique(round(seq(
    1L, length(counts), length.out = min(3L, length(counts))
  )))]
  rows <- vapply(patterns, function(key) which(keys == key)[1L], integer(1))
  for (j in seq_len(ncol(X))) rows <- c(rows, which.min(X[, j]), which.max(X[, j]))
  rows <- unique(c(rows, seq_len(nrow(X))))
  head(rows, min(as.integer(max_points), nrow(X)))
}

## Fixed common integration randomness makes each signature a deterministic
## function of the retained state. Fresh observation or GP noise must not make
## a poorly mixing parameter chain appear to have a large effective sample size.
mixedgp_study2_target_series <- function(fit, X, C, U = NULL,
                                          max_draws_per_chain = Inf,
                                          n_latent = 64L, seed = 481516L) {
  caller_rng <- mixedgp_rng_state()
  on.exit(mixedgp_restore_rng_state(caller_rng), add = TRUE)
  X <- as.matrix(X)
  C <- as.matrix(C)
  G <- nrow(X)
  d <- fit$data$d
  M <- as.integer(n_latent)
  if (G < 1L || M < 2L || nrow(C) != G) stop("Invalid target diagnostic panel.")
  X_std <- sweep(sweep(X, 2L, fit$data$X_center, "-"),
                 2L, fit$data$X_scale, "/")
  repeated <- rep(seq_len(G), each = M)
  calibrated <- length(fit$data$calib_idx) > 0L
  U_std <- if (!is.null(U)) {
    uc <- if (is.null(fit$data$U_center)) rep(0, d) else fit$data$U_center
    us <- if (is.null(fit$data$U_scale)) rep(1, d) else fit$data$U_scale
    sweep(sweep(as.matrix(U), 2L, uc, "-"), 2L, us, "/")
  } else NULL
  if (!is.null(U_std) && !identical(dim(U_std), c(G, d))) {
    stop("U does not match the target diagnostic panel.")
  }
  by_chain <- fit$mcmc$samples_by_chain
  names_target <- c(paste0("m_conditional_mean[", seq_len(G), "]"),
                    paste0("m_conditional_variance[", seq_len(G), "]"),
                    paste0("predictive_second_moment[", seq_len(G), "]"))
  if (calibrated && !is.null(U_std)) names_target <- c(
    names_target, paste0("f_conditional_mean[", seq_len(G), "]"),
    paste0("f_conditional_variance[", seq_len(G), "]")
  )
  if (calibrated) for (k in seq_len(d)) names_target <- c(
    names_target, paste0("prospective_U_mean[", seq_len(G), ",", k, "]"),
    paste0("prospective_U_second_moment[", seq_len(G), ",", k, "]")
  )
  series <- setNames(lapply(names_target, function(z) vector("list", length(by_chain$U))),
                     names_target)
  attr(series, "draw_window_complete") <- TRUE
  kernel <- kernel_spec_from_fit(fit)
  for (id in seq_along(by_chain$U)) {
    n_saved <- dim(by_chain$U[[id]])[1L]
    keep <- mixedgp_diagnostic_draw_indices(n_saved, max_draws_per_chain)
    if (length(keep) < n_saved) attr(series, "draw_window_complete") <- FALSE
    values <- matrix(NA_real_, length(keep), length(names_target),
                     dimnames = list(NULL, names_target))
    for (ii in seq_along(keep)) {
      s <- keep[ii]
      tau <- unflatten_tau(by_chain$tau[[id]][s, ], fit$data$m_vec)
      A <- matrix(by_chain$A[[id]][s, , ], fit$data$q, d)
      set.seed(seed)
      u_mc <- sample_u_given_c_ordprobit_minimax(C[repeated, , drop = FALSE],
                                                A, tau, fit$data$m_vec)
      U_train <- matrix(by_chain$U[[id]][s, , ], nrow(fit$data$X), d)
      pred <- gp_predict_draw_general(
        X_train = fit$data$X, U_train = U_train, y_train = fit$data$y,
        X_star = X_std[repeated, , drop = FALSE], U_star = u_mc,
        logtheta = by_chain$logtheta[[id]][s, ],
        sigma2_eps = by_chain$sigma2[[id]][s], noisy = TRUE,
        kernel = kernel$name, matern_nu = kernel$matern_nu, return_cov = FALSE
      )
      mu <- fit$data$y_center + fit$data$y_scale * pred$mean
      variance <- fit$data$y_scale^2 * pred$var
      mp <- gp_integrated_mean_state_general(
        X_train = fit$data$X, U_train = U_train, y_train = fit$data$y,
        X_star = X_std, U_mc = array(u_mc, c(M, G, d)),
        logtheta = by_chain$logtheta[[id]][s, ],
        sigma2_eps = by_chain$sigma2[[id]][s], kernel = kernel$name,
        matern_nu = kernel$matern_nu
      )
      values[ii, paste0("m_conditional_mean[", seq_len(G), "]")] <-
        fit$data$y_center + fit$data$y_scale * mp$mean
      values[ii, paste0("m_conditional_variance[", seq_len(G), "]")] <-
        fit$data$y_scale^2 * mp$var
      values[ii, paste0("predictive_second_moment[", seq_len(G), "]")] <-
        colMeans(matrix(variance + mu^2, M, G))
      if (calibrated && !is.null(U_std)) {
        fp <- gp_predict_draw_general(
          X_train = fit$data$X, U_train = U_train, y_train = fit$data$y,
          X_star = X_std, U_star = U_std,
          logtheta = by_chain$logtheta[[id]][s, ],
          sigma2_eps = by_chain$sigma2[[id]][s], noisy = FALSE,
          kernel = kernel$name, matern_nu = kernel$matern_nu, return_cov = FALSE
        )
        values[ii, paste0("f_conditional_mean[", seq_len(G), "]")] <-
          fit$data$y_center + fit$data$y_scale * fp$mean
        values[ii, paste0("f_conditional_variance[", seq_len(G), "]")] <-
          fit$data$y_scale^2 * fp$var
      }
      if (calibrated) for (k in seq_len(d)) {
        uc <- if (is.null(fit$data$U_center)) 0 else fit$data$U_center[k]
        us <- if (is.null(fit$data$U_scale)) 1 else fit$data$U_scale[k]
        u_physical <- uc + us * u_mc[, k]
        values[ii, paste0("prospective_U_mean[", seq_len(G), ",", k, "]")] <-
          colMeans(matrix(u_physical, M, G))
        values[ii, paste0("prospective_U_second_moment[", seq_len(G), ",", k, "]")] <-
          colMeans(matrix(u_physical^2, M, G))
      }
    }
    for (name in names_target) series[[name]][[id]] <- values[, name]
  }
  mixedgp_mark_diagnostic_window(series, attr(series, "draw_window_complete"))
}
