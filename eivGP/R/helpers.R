

# ---- core/00_study1_functions.R ----
############################################################
## 00_study1_functions.R
##
## EIV-GP sampler: univariate latent u / ordinal thresholds,
## covariates x of any dimension. A numeric vector x is treated
## as n x 1 and recovers the original Study I model.
############################################################


safe_max <- function(x) {
  if (length(x) == 0 || all(!is.finite(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}


safe_min <- function(x) {
  if (length(x) == 0 || all(!is.finite(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}


safe_median <- function(x) {
  if (length(x) == 0 || all(!is.finite(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}


normalize_gp_kernel_1d <- function(kernel = c("se", "matern"), matern_nu = 2.5) {
  if (length(kernel) < 1L || is.na(kernel[1])) stop("kernel must be specified.")
  key <- tolower(gsub("[^a-z0-9]", "", as.character(kernel[1])))

  if (key %in% c("se", "rbf", "gaussian", "squaredexponential")) {
    return(list(name = "se", matern_nu = NA_real_))
  }
  if (!key %in% c("matern", "maternkernel")) {
    stop("kernel must be 'se' or 'matern'.")
  }

  matern_nu <- as.numeric(matern_nu)
  if (length(matern_nu) != 1L || !is.finite(matern_nu) || matern_nu <= 0) {
    stop("matern_nu must be one positive finite number.")
  }

  list(name = "matern", matern_nu = matern_nu)
}


kernel_spec_from_fit_1d <- function(fit_obj) {
  if (is.null(fit_obj$kernel)) return(normalize_gp_kernel_1d("se"))
  normalize_gp_kernel_1d(fit_obj$kernel$name, fit_obj$kernel$matern_nu)
}


kernel_from_weighted_sqdist_1d <- function(D2, kernel, matern_nu = 2.5) {
  spec <- normalize_gp_kernel_1d(kernel, matern_nu)
  D2 <- pmax(as.matrix(D2), 0)

  if (spec$name == "se") return(exp(-D2))

  r <- sqrt(D2)
  nu <- spec$matern_nu
  if (isTRUE(all.equal(nu, 0.5))) return(exp(-r))
  if (isTRUE(all.equal(nu, 1.5))) {
    z <- sqrt(3) * r
    return((1 + z) * exp(-z))
  }
  if (isTRUE(all.equal(nu, 2.5))) {
    z <- sqrt(5) * r
    return((1 + z + z^2 / 3) * exp(-z))
  }

  z <- sqrt(2 * nu) * r
  out <- matrix(1, nrow(D2), ncol(D2))
  positive <- z > 0
  out[positive] <-
    2^(1 - nu) / gamma(nu) *
    z[positive]^nu * besselK(z[positive], nu = nu)
  out[!is.finite(out)] <- 0
  out
}


safe_chol_1d <- function(A, jitter = 0) {
  n <- nrow(A)

  ans <- try(chol(A), silent = TRUE)
  if (!inherits(ans, "try-error")) {
    attr(ans, "jitter") <- 0
    return(ans)
  }

  if (!is.finite(jitter) || jitter <= 0) {
    stop("Cholesky decomposition failed without numerical regularization.")
  }
  for (k in 0:8) {
    jj <- jitter * 10^k
    ans <- try(chol(A + jj * diag(n)), silent = TRUE)
    if (!inherits(ans, "try-error")) {
      attr(ans, "jitter") <- jj
      return(ans)
    }
  }
  stop("Cholesky decomposition failed.")
}


weighted_distance_multix <- function(Dlist, theta) {
  if (length(Dlist) != length(theta)) stop("Distance and theta dimensions differ")
  out <- matrix(0, nrow(Dlist[[1]]), ncol(Dlist[[1]]))
  for (j in seq_along(Dlist)) out <- out + theta[j] * Dlist[[j]]
  out
}


prepare_gp_block_1d <- function(y, u, Dx_list, logtheta, theta_spec,
                                block_idx, kernel = "se", matern_nu = 2.5,
                                gp_cache = NULL) {
  n <- length(y)
  b <- length(block_idx)
  ## Dense calculations are cheaper for very small data or large blocks.
  if (is.null(gp_cache) || !isTRUE(gp_cache$use_blocks) || n < 40L ||
      b < 1L || b >= n || b > ceiling(sqrt(n))) return(NULL)
  if (anyNA(block_idx) || anyDuplicated(block_idx) ||
      any(block_idx < 1L | block_idx > n)) stop("Invalid GP block indices.")
  if (!is.list(Dx_list)) Dx_list <- list(as.matrix(Dx_list))
  if (is.null(theta_spec)) theta_spec <- make_theta_spec_multix(length(Dx_list))
  current <- gp_state_1d(
    y, u, Dx_list, logtheta, 1, theta_spec, kernel, matern_nu,
    gp_cache = gp_cache, require_factor = FALSE
  )
  rest <- setdiff(seq_len(n), block_idx)
  chol_rest <- try(chol(current$A[rest, rest, drop = FALSE]), silent = TRUE)
  if (inherits(chol_rest, "try-error")) {
    gp_cache$block_fallbacks <- gp_cache$block_fallbacks + 1L
    return(NULL)
  }
  gp_cache$block_setups <- gp_cache$block_setups + 1L
  alpha_rest <- solve_chol(chol_rest, y[rest])
  Dx <- weighted_distance_multix(Dx_list, exp(logtheta[theta_spec$x_index]))
  list(y = y, u = u, Dx_list = Dx_list, logtheta = logtheta,
       theta_spec = theta_spec, kernel = kernel, matern_nu = matern_nu,
       block_idx = block_idx, rest = rest, n = n, b = b,
       rho2 = exp(logtheta[1L])^2,
       theta_u = exp(logtheta[theta_spec$u_index]),
       Dx_cross = Dx[block_idx, rest, drop = FALSE],
       Dx_block = Dx[block_idx, block_idx, drop = FALSE],
       A = current$A, chol_rest = chol_rest, alpha_rest = alpha_rest,
       quad_rest = sum(y[rest] * alpha_rest),
       logdet_rest = 2 * sum(log(diag(chol_rest))))
}


gp_predict_draw <- function(x_train, u_train, y_train,
                            x_star, u_star,
                            logtheta, sigma2_eps,
                            noisy = FALSE,
                            kernel = "se",
                            matern_nu = 2.5,
                            return_cov = FALSE) {
  x_train <- as_study1_numeric_matrix(x_train, "x_train")
  n <- length(y_train)
  if (nrow(x_train) != n) stop("x_train and y_train have different row counts.")
  p_x <- ncol(x_train)
  x_star <- as_study1_numeric_matrix(
    x_star,
    "x_star",
    ncol_expected = p_x,
    expected_names = colnames(x_train)
  )
  if (length(u_train) != n || any(!is.finite(u_train))) {
    stop("u_train must contain one finite value per training row.")
  }
  if (length(u_star) != nrow(x_star) || any(!is.finite(u_star))) {
    stop("u_star must contain one finite value per prediction row.")
  }
  if (length(logtheta) != p_x + 2L || any(!is.finite(logtheta))) {
    stop("logtheta must contain log(rho), one x parameter per column, and theta_u.")
  }
  if (length(sigma2_eps) != 1L || !is.finite(sigma2_eps) || sigma2_eps <= 0) {
    stop("sigma2_eps must be positive and finite.")
  }
  kernel_spec <- normalize_gp_kernel_1d(kernel, matern_nu)
  theta_spec <- make_theta_spec_multix(p_x)

  rho <- exp(logtheta[1])
  theta_x <- exp(logtheta[theta_spec$x_index])
  theta_u <- exp(logtheta[theta_spec$u_index])

  Dx_train <- lapply(seq_len(p_x), function(j) {
    pairwise_sqdist(x_train[, j, drop = FALSE])
  })
  Du_train <- pairwise_sqdist(matrix(u_train, ncol = 1))

  D2_train <- weighted_distance_multix(Dx_train, theta_x) + theta_u * Du_train
  R <- kernel_from_weighted_sqdist_1d(
    D2_train, kernel_spec$name, kernel_spec$matern_nu
  )
  K <- rho^2 * sigma2_eps * R
  C <- K + sigma2_eps * diag(n)
  
  U <- safe_chol_1d(C)
  alpha <- solve_chol(U, y_train)
  
  Dxs <- lapply(seq_len(p_x), function(j) {
    pairwise_sqdist(x_star[, j, drop = FALSE], x_train[, j, drop = FALSE])
  })
  Dus <- pairwise_sqdist(
    matrix(u_star, ncol = 1),
    matrix(u_train, ncol = 1)
  )
  
  D2_star <- weighted_distance_multix(Dxs, theta_x) + theta_u * Dus
  R_star <- kernel_from_weighted_sqdist_1d(
    D2_star, kernel_spec$name, kernel_spec$matern_nu
  )
  K_star <- rho^2 * sigma2_eps * R_star
  
  mu <- as.numeric(K_star %*% alpha)
  
  v <- forwardsolve(t(U), t(K_star))
  var_lat <- rho^2 * sigma2_eps - colSums(v^2)
  var_tolerance <- 100 * (n + nrow(x_star)) * .Machine$double.eps *
    max(1, rho^2 * sigma2_eps)
  if (any(var_lat < -var_tolerance)) {
    stop("GP conditional variance is materially negative.")
  }
  var_lat <- pmax(var_lat, 0)
  
  if (!isTRUE(return_cov)) {
    if (noisy) var_lat <- var_lat + sigma2_eps
    return(list(mean = mu, var = var_lat))
  }

  Dxx <- lapply(seq_len(p_x), function(j) {
    pairwise_sqdist(x_star[, j, drop = FALSE])
  })
  Duu <- pairwise_sqdist(matrix(u_star, ncol = 1L))
  D2_ss <- weighted_distance_multix(Dxx, theta_x) + theta_u * Duu
  R_ss <- kernel_from_weighted_sqdist_1d(
    D2_ss, kernel_spec$name, kernel_spec$matern_nu
  )
  cov_lat <- rho^2 * sigma2_eps * R_ss - crossprod(v)
  cov_lat <- 0.5 * (cov_lat + t(cov_lat))
  if (any(diag(cov_lat) < -var_tolerance)) {
    stop("GP conditional covariance has a materially negative diagonal.")
  }
  diag(cov_lat) <- pmax(diag(cov_lat), 0)
  if (noisy) cov_lat <- cov_lat + sigma2_eps * diag(nrow(x_star))

  list(mean = mu, var = diag(cov_lat), cov = cov_lat)
}


gp_mle_predict_1d <- function(fit, Xstar, noisy = FALSE) {
  X <- fit$X
  y <- fit$y

  n <- nrow(X)
  d <- ncol(X)
  Xstar <- as_study1_numeric_matrix(
    Xstar,
    "Xstar",
    ncol_expected = d,
    expected_names = colnames(X)
  )
  N <- nrow(Xstar)
  kernel_spec <- kernel_spec_from_fit_1d(fit)
  
  par <- fit$par
  
  sigma2 <- exp(par[1])
  rho <- exp(par[2])
  theta <- exp(par[-c(1, 2)])
  
  Rexp <- matrix(0, n, n)
  
  for (j in seq_len(d)) {
    Rexp <- Rexp + theta[j] * pairwise_sqdist(X[, j, drop = FALSE])
  }
  
  R <- kernel_from_weighted_sqdist_1d(
    Rexp, kernel_spec$name, kernel_spec$matern_nu
  )
  K <- rho^2 * sigma2 * R
  C <- K + sigma2 * diag(n)
  
  U <- safe_chol_1d(C)
  alpha <- solve_chol(U, y)
  
  Rstar_exp <- matrix(0, N, n)
  
  for (j in seq_len(d)) {
    Rstar_exp <- Rstar_exp +
      theta[j] * pairwise_sqdist(
        Xstar[, j, drop = FALSE],
        X[, j, drop = FALSE]
      )
  }
  
  Rstar <- kernel_from_weighted_sqdist_1d(
    Rstar_exp, kernel_spec$name, kernel_spec$matern_nu
  )
  Kstar <- rho^2 * sigma2 * Rstar
  
  mu <- as.numeric(Kstar %*% alpha)
  
  v <- forwardsolve(t(U), t(Kstar))
  var_lat <- rho^2 * sigma2 - colSums(v^2)
  var_tolerance <- 100 * (n + N) * .Machine$double.eps *
    max(1, rho^2 * sigma2)
  if (any(var_lat < -var_tolerance)) {
    stop("GP conditional variance is materially negative.")
  }
  var_lat <- pmax(var_lat, 0)
  
  if (noisy) {
    var_lat <- var_lat + sigma2
  }
  
  list(mean = mu, var = var_lat)
}


sample_gp_mle_predictive_1d <- function(fit, Xstar, n_draw = 1000) {
  pred <- gp_mle_predict_1d(fit, Xstar = Xstar, noisy = TRUE)
  
  n_test <- length(pred$mean)
  
  out <- matrix(
    rnorm(
      n_draw * n_test,
      mean = rep(pred$mean, each = n_draw),
      sd = rep(sqrt(pred$var), each = n_draw)
    ),
    nrow = n_draw,
    ncol = n_test
  )
  attr(out, "conditional_means") <- matrix(pred$mean, nrow = 1L)
  attr(out, "conditional_vars") <- matrix(pred$var, nrow = 1L)
  attr(out, "mixture_components") <- "fitted GP Gaussian predictor"
  out
}


conditional_mean_scores <- function(c_ord, m) {
  counts <- tabulate(c_ord, nbins = m)
  probs <- cumsum(counts) / sum(counts)
  probs <- pmin(pmax(probs[seq_len(m - 1)], 1e-4), 1 - 1e-4)
  
  tau_hat <- qnorm(probs)
  
  lower <- c(-Inf, tau_hat)
  upper <- c(tau_hat, Inf)
  
  denom <- pnorm(upper) - pnorm(lower)
  numer <- dnorm(lower) - dnorm(upper)
  
  numer / pmax(denom, .Machine$double.eps)
}


make_monotone_scores_1d <- function(a, m) {
  if (m == 1) {
    return(0)
  }
  
  if (m == 2) {
    return(c(0, 1))
  }
  
  logits <- c(a, 0)
  e <- exp(logits - max(logits))
  inc <- e / sum(e)
  
  c(0, cumsum(inc))
}


gp_mle_predict_learned_embedding <- function(fit, x_star, c_star,
                                             noisy = TRUE) {
  X_train <- cbind(fit$x, fit$z_scores[fit$c_ord])
  X_star <- cbind(x_star, fit$z_scores[c_star])
  predictor_names <- colnames(fit$x)
  if (is.null(predictor_names)) {
    predictor_names <- paste0("x", seq_len(ncol(fit$x)))
  }
  embedding_names <- c(predictor_names, ".ordinal_score")
  colnames(X_train) <- embedding_names
  colnames(X_star) <- embedding_names
  
  fake_fit <- list(
    par = fit$par,
    X = X_train,
    y = fit$y,
    kernel = fit$kernel
  )
  
  gp_mle_predict_1d(fake_fit, Xstar = X_star, noisy = noisy)
}


sample_gp_learned_embedding_predictive <- function(fit, x_star, c_star,
                                                   n_draw = 1000) {
  pred <- gp_mle_predict_learned_embedding(
    fit,
    x_star = x_star,
    c_star = c_star,
    noisy = TRUE
  )
  
  n_test <- length(pred$mean)
  
  matrix(
    rnorm(
      n_draw * n_test,
      mean = rep(pred$mean, each = n_draw),
      sd = rep(sqrt(pred$var), each = n_draw)
    ),
    nrow = n_draw,
    ncol = n_test
  )
}


predict_embedding_baseline_samples <- function(baselines,
                                               x_star_raw,
                                               c_star,
                                               m,
                                               n_draw = 1000) {
  x_star_raw <- as_study1_numeric_matrix(
    x_star_raw,
    "x_star_raw",
    ncol_expected = length(baselines$x_center),
    expected_names = baselines$predictor_names
  )
  x_star <- sweep(
    sweep(x_star_raw, 2, baselines$x_center, "-"),
    2, baselines$x_scale, "/"
  )
  
  z_gauss_star <- qnorm(c_star / (m + 1))
  X_gauss_star <- cbind(x_star, z_gauss_star)
  colnames(X_gauss_star) <- colnames(baselines$fit_gauss$X)
  
  draws_gauss_std <- sample_gp_mle_predictive_1d(
    baselines$fit_gauss,
    Xstar = X_gauss_star,
    n_draw = n_draw
  )
  
  draws_gauss <- baselines$y_center + baselines$y_scale * draws_gauss_std
  
  z_cm_star <- baselines$z_cm_scores[c_star]
  X_cm_star <- cbind(x_star, z_cm_star)
  colnames(X_cm_star) <- colnames(baselines$fit_cm$X)
  
  draws_cm_std <- sample_gp_mle_predictive_1d(
    baselines$fit_cm,
    Xstar = X_cm_star,
    n_draw = n_draw
  )
  
  draws_cm <- baselines$y_center + baselines$y_scale * draws_cm_std
  
  draws_learned_std <- sample_gp_learned_embedding_predictive(
    baselines$fit_learned,
    x_star = x_star,
    c_star = c_star,
    n_draw = n_draw
  )
  
  draws_learned <- baselines$y_center + baselines$y_scale * draws_learned_std
  
  list(
    `GP-Gaussian` = draws_gauss,
    `GP-CondMean` = draws_cm,
    `GP-LearnedEmb` = draws_learned
  )
}


interval_score <- function(lo, hi, y, alpha = 0.05) {
  width <- hi - lo
  
  width +
    2 / alpha * (lo - y) * (y < lo) +
    2 / alpha * (y - hi) * (y > hi)
}


summarize_predictive_samples_1d <- function(draw_mat, y_true,
                                         method,
                                         rep_id,
                                         n_calib,
                                         scenario) {
  component_means <- attr(draw_mat, "conditional_means", exact = TRUE)
  component_vars <- attr(draw_mat, "conditional_vars", exact = TRUE)
  draw_mat <- as.matrix(draw_mat)
  y_true <- as.numeric(y_true)
  if (ncol(draw_mat) != length(y_true)) {
    stop("draw_mat must have one column for every element of y_true.")
  }
  nlpd <- normal_mixture_nlpd(y_true, component_means, component_vars)
  pred_mean <- colMeans(draw_mat)
  
  qs <- apply(
    draw_mat,
    2,
    quantile,
    probs = c(0.025, 0.975),
    names = FALSE
  )
  
  lo <- qs[1, ]
  hi <- qs[2, ]
  
  data.frame(
    rep = rep_id,
    scenario = scenario,
    n_calib = n_calib,
    method = method,
    RMSE = sqrt(mean((pred_mean - y_true)^2)),
    MAE = mean(abs(pred_mean - y_true)),
    Coverage95 = mean(y_true >= lo & y_true <= hi),
    Width95 = mean(hi - lo),
    CRPS = mean(crps_sample_matrix(draw_mat, y_true)),
    NLPD = nlpd$value,
    NLPD_reason = nlpd$reason,
    IntervalScore95 = mean(interval_score(lo, hi, y_true)),
    stringsAsFactors = FALSE
  )
}


validate_control_ordprobit <- function(control, sampler_strategy, n_mis, n_iter) {
  integer_control <- function(name, minimum) {
    value <- control[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value != floor(value) || value < minimum) {
      stop(name, " must be one integer >= ", minimum, ".")
    }
  }
  for (name in c("score_update_every", "A_update_every", "tau_update_every",
                 "theta_update_every", "adapt_every", "adapt_window",
                 "max_ess_try")) integer_control(name, 1L)
  for (name in c("n_u_blocks_per_iter", "joint_local_blocks_per_iter",
                 "global_u_every", "joint_theta_every", "joint_global_every",
                 "joint_collapsed_every", "marginal_measurement_every",
                 "joint_measurement_every", "loading_transport_every")) {
    integer_control(name, 0L)
  }
  integer_control("u_block_size", if (n_mis > 0L) 1L else 0L)
  for (name in c("s_A", "tau_bound", "theta_width_min", "theta_width_max")) {
    value <- control[[name]]
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value <= 0) stop(name, " must be positive and finite.")
  }
  if (!is.logical(control$adapt_theta_width) ||
      length(control$adapt_theta_width) != 1L || is.na(control$adapt_theta_width)) {
    stop("adapt_theta_width must be one non-missing logical value.")
  }
  if (control$theta_width_min > control$theta_width_max) {
    stop("theta_width_min must not exceed theta_width_max.")
  }
  scheduled <- function(value) value > 0L && value <= n_iter
  for (name in c("score_update_every", "A_update_every", "tau_update_every")) {
    if (!scheduled(control[[name]])) stop(name, " must not exceed n_iter.")
  }
  if (sampler_strategy == "interwoven") {
    theta_active <- any(vapply(control[c("joint_theta_every", "joint_global_every",
                                        "joint_collapsed_every")], scheduled, logical(1)))
    u_active <- control$joint_local_blocks_per_iter > 0L ||
      any(vapply(control[c("joint_global_every", "joint_collapsed_every")],
                 scheduled, logical(1)))
  } else {
    theta_active <- scheduled(control$theta_update_every)
    u_active <- control$n_u_blocks_per_iter > 0L || scheduled(control$global_u_every)
    if (control$joint_measurement_every > 0L || control$loading_transport_every > 0L) {
      stop("joint_measurement_every and loading_transport_every require interwoven sampling.")
    }
  }
  if (!theta_active) stop("The schedule freezes all GP hyperparameter moves.")
  if (n_mis > 0L && !u_active) stop("The schedule freezes all missing-U moves.")
  invisible(TRUE)
}


############################################################
## General utilities
############################################################

safe_max <- function(x) {
  if (length(x) == 0 || all(!is.finite(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}


safe_min <- function(x) {
  if (length(x) == 0 || all(!is.finite(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}


safe_median <- function(x) {
  if (length(x) == 0 || all(!is.finite(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}


prepare_ordinal_matrix <- function(C,
                                   m_vec = NULL,
                                   level_maps = NULL,
                                   expected_names = NULL,
                                   name = "C_ord") {
  input_names <- if (is.null(dim(C))) NULL else colnames(C)
  if (!is.null(input_names) &&
      (anyNA(input_names) || any(!nzchar(input_names)) ||
       anyDuplicated(input_names))) {
    stop(
      name,
      " column names must be complete and unique when any names are supplied."
    )
  }
  has_explicit_names <- !is.null(input_names)

  if (is.null(dim(C))) {
    C_df <- data.frame(C, check.names = FALSE)
  } else if (is.data.frame(C)) {
    C_df <- C
  } else {
    C_df <- as.data.frame(C, stringsAsFactors = FALSE, check.names = FALSE)
  }

  q <- ncol(C_df)
  if (q < 1L) stop(name, " must contain at least one ordinal column.")

  if (!is.null(expected_names)) {
    if (length(expected_names) != q) {
      stop(name, " has the wrong number of columns.")
    }
    current_names <- names(C_df)
    if (has_explicit_names &&
        setequal(current_names, expected_names)) {
      C_df <- C_df[, expected_names, drop = FALSE]
    } else if (has_explicit_names && !identical(current_names, expected_names)) {
      stop(name, " column names do not match the training ordinal proxies.")
    } else if (!has_explicit_names) {
      names(C_df) <- expected_names
    }
  }

  if (is.null(names(C_df)) || any(!nzchar(names(C_df)))) {
    names(C_df) <- paste0("c", seq_len(q))
  }

  if (!is.null(m_vec)) {
    m_vec <- mixedgp_as_integer_strict(
      m_vec, "m_vec", min_value = 2L, length_expected = q
    )
  }

  if (!is.null(level_maps)) {
    if (!is.list(level_maps) || length(level_maps) != q) {
      stop("level_maps must be a list with one element per ordinal proxy.")
    }
  } else {
    level_maps <- vector("list", q)
  }

  encoded <- matrix(NA_integer_, nrow(C_df), q)
  colnames(encoded) <- names(C_df)
  names(level_maps) <- names(C_df)

  for (j in seq_len(q)) {
    value <- C_df[[j]]
    if (anyNA(value)) stop(name, " contains missing ordinal values.")

    levels_j <- level_maps[[j]]
    if (is.null(levels_j)) {
      if (is.ordered(value)) {
        levels_j <- levels(value)
      } else if (is.numeric(value) || is.integer(value) || is.logical(value)) {
        if (!is.null(m_vec)) {
          levels_j <- seq_len(m_vec[j])
        } else {
          integer_codes <- is.numeric(value) &&
            all(is.finite(value)) &&
            all(value == floor(value)) &&
            min(value) == 1 &&
            max(value) <= 50
          levels_j <- if (integer_codes) {
            seq_len(max(value))
          } else {
            sort(unique(value))
          }
        }
      } else {
        stop(
          name, " column '", names(C_df)[j],
          "' must be an ordered factor, numeric code, or have an explicit level map."
        )
      }
    }

    if (length(levels_j) < 2L || anyDuplicated(as.character(levels_j))) {
      stop("Each ordinal level map must contain at least two distinct levels.")
    }
    if (!is.null(m_vec) && length(levels_j) != m_vec[j]) {
      stop("The level map for '", names(C_df)[j], "' does not match m_vec.")
    }

    code_j <- match(as.character(value), as.character(levels_j))
    if (anyNA(code_j)) {
      unknown <- unique(as.character(value[is.na(code_j)]))
      stop(
        name, " column '", names(C_df)[j], "' contains unknown level(s): ",
        paste(unknown, collapse = ", ")
      )
    }

    encoded[, j] <- code_j
    level_maps[[j]] <- levels_j
  }

  inferred_m <- vapply(level_maps, length, integer(1))
  if (is.null(m_vec)) m_vec <- inferred_m

  list(
    C = encoded,
    m_vec = as.integer(m_vec),
    level_maps = level_maps,
    column_names = colnames(encoded)
  )
}


normalize_gp_kernel <- function(kernel = c("se", "matern"), matern_nu = 2.5) {
  if (length(kernel) < 1L || is.na(kernel[1])) stop("kernel must be specified.")
  key <- tolower(gsub("[^a-z0-9]", "", as.character(kernel[1])))

  if (key %in% c("se", "rbf", "gaussian", "squaredexponential")) {
    return(list(name = "se", matern_nu = NA_real_))
  }
  if (!key %in% c("matern", "maternkernel")) {
    stop("kernel must be 'se' or 'matern'.")
  }

  matern_nu <- as.numeric(matern_nu)
  if (length(matern_nu) != 1L || !is.finite(matern_nu) || matern_nu <= 0) {
    stop("matern_nu must be one positive finite number.")
  }

  list(name = "matern", matern_nu = matern_nu)
}


kernel_spec_from_fit <- function(fit_obj) {
  if (is.null(fit_obj$kernel)) return(normalize_gp_kernel("se"))
  normalize_gp_kernel(fit_obj$kernel$name, fit_obj$kernel$matern_nu)
}


kernel_from_weighted_sqdist <- function(D2, kernel, matern_nu = 2.5) {
  spec <- normalize_gp_kernel(kernel, matern_nu)
  D2 <- pmax(as.matrix(D2), 0)

  if (spec$name == "se") return(exp(-D2))

  r <- sqrt(D2)
  nu <- spec$matern_nu
  if (isTRUE(all.equal(nu, 0.5))) return(exp(-r))
  if (isTRUE(all.equal(nu, 1.5))) {
    z <- sqrt(3) * r
    return((1 + z) * exp(-z))
  }
  if (isTRUE(all.equal(nu, 2.5))) {
    z <- sqrt(5) * r
    return((1 + z + z^2 / 3) * exp(-z))
  }

  z <- sqrt(2 * nu) * r
  out <- matrix(1, nrow(D2), ncol(D2))
  positive <- z > 0
  out[positive] <-
    2^(1 - nu) / gamma(nu) *
    z[positive]^nu * besselK(z[positive], nu = nu)
  out[!is.finite(out)] <- 0
  out
}


safe_chol <- function(A, jitter = 0) {
  n <- nrow(A)

  ans <- try(chol(A), silent = TRUE)
  if (!inherits(ans, "try-error")) {
    attr(ans, "jitter") <- 0
    return(ans)
  }
  
  if (!is.finite(jitter) || jitter <= 0) {
    stop("Cholesky decomposition failed without numerical regularization.")
  }
  for (k in 0:8) {
    jj <- jitter * 10^k
    ans <- try(chol(A + jj * diag(n)), silent = TRUE)
    if (!inherits(ans, "try-error")) {
      attr(ans, "jitter") <- jj
      return(ans)
    }
  }
  
  stop("Cholesky decomposition failed.")
}


############################################################
## General ordinal-probit measurement model
############################################################

empirical_normal_scores <- function(C, m_vec = NULL) {
  C <- as.matrix(C)
  q <- ncol(C)
  
  if (is.null(m_vec)) {
    m_vec <- apply(C, 2, max, na.rm = TRUE)
  }
  
  scores <- vector("list", q)
  
  for (j in seq_len(q)) {
    counts <- tabulate(C[, j], nbins = m_vec[j])
    n <- sum(counts)
    cum_counts <- cumsum(counts)
    
    mid_probs <- (cum_counts - 0.5 * counts) / n
    mid_probs <- pmin(pmax(mid_probs, 1e-4), 1 - 1e-4)
    
    scores[[j]] <- qnorm(mid_probs)
  }
  
  scores
}


sample_scores_ord <- function(C, U, A, tau) {
  C <- as.matrix(C)
  U <- as.matrix(U)
  A <- as.matrix(A)
  
  n <- nrow(C)
  q <- ncol(C)
  
  S <- matrix(NA_real_, n, q)
  
  for (j in seq_len(q)) {
    mu <- as.numeric(U %*% A[j, ])
    
    tau_j <- tau[[j]]
    
    lower <- c(-Inf, tau_j)[C[, j]]
    upper <- c(tau_j, Inf)[C[, j]]
    
    S[, j] <- rtruncnorm_vec(mu, 1, lower, upper)
  }
  
  S
}


gp_predict_draw_general <- function(X_train,
                                    U_train,
                                    y_train,
                                    X_star,
                                    U_star,
                                    logtheta,
                                    sigma2_eps,
                                    noisy = FALSE,
                                    kernel = "se",
                                    matern_nu = 2.5,
                                    return_cov = FALSE) {
  X_train <- as.matrix(X_train)
  U_train <- as.matrix(U_train)
  X_star <- as.matrix(X_star)
  U_star <- as.matrix(U_star)
  
  n <- nrow(X_train)
  N <- nrow(X_star)
  p <- ncol(X_train)
  d <- ncol(U_train)

  if (nrow(U_train) != n || length(y_train) != n) {
    stop("X_train, U_train, and y_train have incompatible row counts.")
  }
  if (nrow(U_star) != N || ncol(X_star) != p || ncol(U_star) != d) {
    stop("X_star and U_star have incompatible dimensions.")
  }
  if (any(!is.finite(c(X_train, U_train, y_train, X_star, U_star)))) {
    stop("GP prediction inputs must be finite.")
  }

  kernel_spec <- normalize_gp_kernel(kernel, matern_nu)
  
  rho <- exp(logtheta[1])
  theta_x <- exp(logtheta[1L + seq_len(p)])
  theta_u <- exp(logtheta[1L + p + seq_len(d)])
  
  D2 <- weighted_sqdist_general(X_train, U_train, theta_x, theta_u)
  R <- kernel_from_weighted_sqdist(
    D2,
    kernel = kernel_spec$name,
    matern_nu = kernel_spec$matern_nu
  )
  
  K <- rho^2 * sigma2_eps * R
  C <- K + sigma2_eps * diag(n)
  
  Uchol <- safe_chol(C)
  alpha <- solve_chol(Uchol, y_train)
  
  D2_star <- weighted_sqdist_general(
    X_star, U_star, theta_x, theta_u,
    X_b = X_train,
    U_b = U_train
  )
  Rstar <- kernel_from_weighted_sqdist(
    D2_star,
    kernel = kernel_spec$name,
    matern_nu = kernel_spec$matern_nu
  )
  Kstar <- rho^2 * sigma2_eps * Rstar
  
  mu <- as.numeric(Kstar %*% alpha)
  
  v <- forwardsolve(t(Uchol), t(Kstar))
  
  var_lat <- rho^2 * sigma2_eps - colSums(v^2)
  var_tolerance <- 100 * (n + N) * .Machine$double.eps *
    max(1, rho^2 * sigma2_eps)
  if (any(var_lat < -var_tolerance)) {
    stop("GP conditional variance is materially negative.")
  }
  var_lat <- pmax(var_lat, 0)

  if (!isTRUE(return_cov)) {
    if (noisy) var_lat <- var_lat + sigma2_eps
    return(list(mean = mu, var = var_lat))
  }

  D2_ss <- weighted_sqdist_general(X_star, U_star, theta_x, theta_u)
  Rss <- kernel_from_weighted_sqdist(
    D2_ss,
    kernel = kernel_spec$name,
    matern_nu = kernel_spec$matern_nu
  )
  cov_lat <- rho^2 * sigma2_eps * Rss - crossprod(v)
  cov_lat <- 0.5 * (cov_lat + t(cov_lat))
  if (any(diag(cov_lat) < -var_tolerance)) {
    stop("GP conditional covariance has a materially negative diagonal.")
  }
  diag(cov_lat) <- pmax(diag(cov_lat), 0)
  if (noisy) cov_lat <- cov_lat + sigma2_eps * diag(N)

  list(mean = mu, var = diag(cov_lat), cov = cov_lat)
}


loading_transport_score_loglik <- function(S, U, A, miss_idx) {
  n <- nrow(U)
  q <- nrow(A)
  calibrated <- setdiff(seq_len(n), miss_idx)
  answer <- 0
  if (length(miss_idx)) {
    covariance_chol <- safe_chol(diag(q) + tcrossprod(A))
    standardized <- forwardsolve(t(covariance_chol), t(S[miss_idx, , drop = FALSE]))
    answer <- -0.5 * (length(miss_idx) *
      (q * log(2 * base::pi) + 2 * sum(log(diag(covariance_chol)))) +
      sum(standardized^2))
  }
  if (length(calibrated)) {
    residual <- S[calibrated, , drop = FALSE] -
      U[calibrated, , drop = FALSE] %*% t(A)
    answer <- answer - 0.5 *
      (length(calibrated) * q * log(2 * base::pi) + sum(residual^2))
  }
  answer
}


draw_study2_score_errors <- function(n, true_params) {
  q <- true_params$q

  if (true_params$score_error == "gaussian") {
    return(matrix(rnorm(n * q), nrow = n, ncol = q))
  }

  if (true_params$score_error == "logistic") {
    ## scale sqrt(3)/pi gives unit variance before multiplication by lambda.
    return(matrix(
      rlogis(n * q, scale = sqrt(3) / base::pi),
      nrow = n,
      ncol = q
    ))
  }

  stop("Unknown Study II score-error distribution: ", true_params$score_error)
}


ordinal_from_scores_matrix_tau <- function(S, tau_mat) {
  S <- as.matrix(S)
  n <- nrow(S)
  q <- ncol(S)
  
  C <- matrix(NA_integer_, n, q)
  
  for (j in seq_len(q)) {
    C[, j] <- as.integer(cut(
      S[, j],
      breaks = c(-Inf, tau_mat[j, ], Inf),
      labels = FALSE
    ))
  }
  
  C
}


validate_oracle_pool_2d <- function(oracle_pool, true_params) {
  expected_provenance <- study2_oracle_pool_provenance(true_params)
  if (!is.list(oracle_pool) ||
      !identical(oracle_pool$provenance, expected_provenance) ||
      !is.matrix(oracle_pool$U) ||
      ncol(oracle_pool$U) != true_params$d ||
      !is.matrix(oracle_pool$C) ||
      nrow(oracle_pool$C) != nrow(oracle_pool$U) ||
      ncol(oracle_pool$C) != true_params$q ||
      anyNA(oracle_pool$U) || any(!is.finite(oracle_pool$U)) ||
      anyNA(oracle_pool$C) ||
      any(oracle_pool$C < 1L | oracle_pool$C > true_params$m)) {
    stop("oracle_pool provenance or dimensions do not match the oracle DGP.")
  }
  expected_key <- pattern_key(oracle_pool$C)
  if (!identical(as.character(oracle_pool$key), expected_key) ||
      !is.list(oracle_pool$split_idx)) {
    stop("oracle_pool category keys are inconsistent with its stored C matrix.")
  }
  expected_split <- split(seq_len(nrow(oracle_pool$C)), expected_key)
  if (!setequal(names(oracle_pool$split_idx), names(expected_split)) ||
      any(!vapply(names(expected_split), function(key) {
        identical(
          as.integer(oracle_pool$split_idx[[key]]),
          as.integer(expected_split[[key]])
        )
      }, logical(1L)))) {
    stop("oracle_pool split indices are inconsistent with its category keys.")
  }
  invisible(TRUE)
}


############################################################
## Legacy deterministic embeddings (not publication competitors)
##
## Retained only to reproduce early pilot analyses.  The publication drivers
## 01_study2_representative_figures.R and 02_study2_monte_carlo.R must not call
## these functions or report them under literature-method labels.
############################################################

conditional_mean_scores_ord <- function(c_ord, m) {
  counts <- tabulate(c_ord, nbins = m)
  probs <- cumsum(counts) / sum(counts)
  probs <- pmin(pmax(probs[seq_len(m - 1)], 1e-4), 1 - 1e-4)
  
  tau_hat <- qnorm(probs)
  
  lower <- c(-Inf, tau_hat)
  upper <- c(tau_hat, Inf)
  
  denom <- pnorm(upper) - pnorm(lower)
  numer <- dnorm(lower) - dnorm(upper)
  
  numer / pmax(denom, .Machine$double.eps)
}


gp_mle_predict <- function(fit, Xstar, noisy = FALSE) {
  X <- fit$X
  y <- fit$y
  Xstar <- as.matrix(Xstar)
  
  n <- nrow(X)
  N <- nrow(Xstar)
  d <- ncol(X)
  
  par <- fit$par
  
  sigma2 <- exp(par[1])
  rho <- exp(par[2])
  theta <- exp(par[-c(1, 2)])
  
  Rexp <- matrix(0, n, n)
  
  for (j in seq_len(d)) {
    Rexp <- Rexp + theta[j] * pairwise_sqdist(X[, j, drop = FALSE])
  }
  
  R <- exp(-Rexp)
  K <- rho^2 * sigma2 * R
  Cmat <- K + sigma2 * diag(n)
  
  Uchol <- safe_chol(Cmat)
  alpha <- solve_chol(Uchol, y)
  
  Rstar_exp <- matrix(0, N, n)
  
  for (j in seq_len(d)) {
    Rstar_exp <- Rstar_exp +
      theta[j] * pairwise_sqdist(
        Xstar[, j, drop = FALSE],
        X[, j, drop = FALSE]
      )
  }
  
  Rstar <- exp(-Rstar_exp)
  Kstar <- rho^2 * sigma2 * Rstar
  
  mu <- as.numeric(Kstar %*% alpha)
  
  v <- forwardsolve(t(Uchol), t(Kstar))
  
  var_lat <- rho^2 * sigma2 - colSums(v^2)
  var_tolerance <- 100 * (n + N) * .Machine$double.eps *
    max(1, rho^2 * sigma2)
  if (any(var_lat < -var_tolerance)) {
    stop("GP conditional variance is materially negative.")
  }
  var_lat <- pmax(var_lat, 0)
  
  if (noisy) var_lat <- var_lat + sigma2
  
  list(mean = mu, var = var_lat)
}


sample_gp_mle_predictive <- function(fit, Xstar, n_draw = 1000) {
  pred <- gp_mle_predict(fit, Xstar = Xstar, noisy = TRUE)
  
  n_test <- length(pred$mean)
  
  out <- matrix(
    rnorm(
      n_draw * n_test,
      mean = rep(pred$mean, each = n_draw),
      sd = rep(sqrt(pred$var), each = n_draw)
    ),
    nrow = n_draw,
    ncol = n_test
  )
  attr(out, "conditional_means") <- matrix(pred$mean, nrow = 1L)
  attr(out, "conditional_vars") <- matrix(pred$var, nrow = 1L)
  attr(out, "mixture_components") <- "fitted GP Gaussian predictor"
  out
}


make_monotone_scores <- function(a, m) {
  if (m == 1) return(0)
  if (m == 2) return(c(0, 1))
  
  logits <- c(a, 0)
  e <- exp(logits - max(logits))
  inc <- e / sum(e)
  
  c(0, cumsum(inc))
}


gp_mle_predict_learned_embedding_ord <- function(fit, X_star, C_star,
                                                 noisy = TRUE) {
  X_star <- as.matrix(X_star)
  C_star <- as.matrix(C_star)
  
  q <- ncol(C_star)
  Z_star <- matrix(NA_real_, nrow(C_star), q)
  
  for (j in seq_len(q)) {
    Z_star[, j] <- fit$z_scores_list[[j]][C_star[, j]]
  }
  
  fake_fit <- list(
    par = fit$par,
    X = fit$X_aug,
    y = fit$y
  )
  
  gp_mle_predict(fake_fit, Xstar = cbind(X_star, Z_star), noisy = noisy)
}


sample_gp_learned_embedding_predictive_ord <- function(fit,
                                                       X_star,
                                                       C_star,
                                                       n_draw = 1000) {
  pred <- gp_mle_predict_learned_embedding_ord(
    fit,
    X_star = X_star,
    C_star = C_star,
    noisy = TRUE
  )
  
  n_test <- length(pred$mean)
  
  matrix(
    rnorm(
      n_draw * n_test,
      mean = rep(pred$mean, each = n_draw),
      sd = rep(sqrt(pred$var), each = n_draw)
    ),
    nrow = n_draw,
    ncol = n_test
  )
}


predict_embedding_baseline_samples_ord <- function(baselines,
                                                   X_star_raw,
                                                   C_star,
                                                   n_draw = 1000) {
  X_star_raw <- as.matrix(X_star_raw)
  C_star <- as.matrix(C_star)
  
  m_vec <- baselines$m_vec
  q <- ncol(C_star)
  
  X_star <- sweep(
    sweep(X_star_raw, 2, baselines$X_center, "-"),
    2,
    baselines$X_scale,
    "/"
  )
  
  Z_gauss_star <- matrix(NA_real_, nrow(C_star), q)
  
  for (j in seq_len(q)) {
    Z_gauss_star[, j] <- qnorm(C_star[, j] / (m_vec[j] + 1))
  }
  
  draws_gauss_std <- sample_gp_mle_predictive(
    baselines$fit_gauss,
    Xstar = cbind(X_star, Z_gauss_star),
    n_draw = n_draw
  )
  
  draws_gauss <- baselines$y_center + baselines$y_scale * draws_gauss_std
  
  Z_cm_star <- matrix(NA_real_, nrow(C_star), q)
  
  for (j in seq_len(q)) {
    Z_cm_star[, j] <- baselines$z_cm_scores[[j]][C_star[, j]]
  }
  
  draws_cm_std <- sample_gp_mle_predictive(
    baselines$fit_cm,
    Xstar = cbind(X_star, Z_cm_star),
    n_draw = n_draw
  )
  
  draws_cm <- baselines$y_center + baselines$y_scale * draws_cm_std
  
  draws_learned_std <- sample_gp_learned_embedding_predictive_ord(
    baselines$fit_learned,
    X_star = X_star,
    C_star = C_star,
    n_draw = n_draw
  )
  
  draws_learned <- baselines$y_center + baselines$y_scale * draws_learned_std
  
  list(
    `GP-Gaussian` = draws_gauss,
    `GP-CondMean` = draws_cm,
    `GP-LearnedEmb` = draws_learned
  )
}


interval_score <- function(lo, hi, y, alpha = 0.05) {
  width <- hi - lo
  
  width +
    2 / alpha * (lo - y) * (y < lo) +
    2 / alpha * (y - hi) * (y > hi)
}


subset_predictive_samples <- function(draw_mat, columns) {
  component_means <- attr(draw_mat, "conditional_means", exact = TRUE)
  component_vars <- attr(draw_mat, "conditional_vars", exact = TRUE)
  out <- as.matrix(draw_mat)[, columns, drop = FALSE]
  if (!is.null(component_means)) {
    attr(out, "conditional_means") <-
      as.matrix(component_means)[, columns, drop = FALSE]
  }
  if (!is.null(component_vars)) {
    attr(out, "conditional_vars") <-
      as.matrix(component_vars)[, columns, drop = FALSE]
  }
  out
}


summarize_predictive_samples <- function(draw_mat,
                                         y_true,
                                         method,
                                         rep_id,
                                         n_calib,
                                         scenario = "study2",
                                         evaluation_stratum = "overall") {
  component_means <- attr(draw_mat, "conditional_means", exact = TRUE)
  component_vars <- attr(draw_mat, "conditional_vars", exact = TRUE)
  draw_mat <- as.matrix(draw_mat)
  y_true <- as.numeric(y_true)

  if (ncol(draw_mat) != length(y_true)) {
    stop("draw_mat must have one column for every element of y_true.")
  }
  nlpd <- normal_mixture_nlpd(y_true, component_means, component_vars)

  pred_mean <- colMeans(draw_mat)
  
  qs <- apply(
    draw_mat,
    2,
    quantile,
    probs = c(0.025, 0.975),
    names = FALSE
  )
  
  lo <- qs[1, ]
  hi <- qs[2, ]
  
  data.frame(
    rep = rep_id,
    scenario = scenario,
    evaluation_stratum = evaluation_stratum,
    n_test_eval = length(y_true),
    n_calib = n_calib,
    method = method,
    RMSE = sqrt(mean((pred_mean - y_true)^2)),
    MAE = mean(abs(pred_mean - y_true)),
    Coverage95 = mean(y_true >= lo & y_true <= hi),
    Width95 = mean(hi - lo),
    CRPS = mean(crps_sample_matrix(draw_mat, y_true)),
    NLPD = nlpd$value,
    NLPD_reason = nlpd$reason,
    IntervalScore95 = mean(interval_score(lo, hi, y_true)),
    stringsAsFactors = FALSE
  )
}


summarize_predictive_samples_by_pattern <- function(draw_mat,
                                                    y_true,
                                                    pattern_stratum,
                                                    method,
                                                    rep_id,
                                                    n_calib,
                                                    scenario) {
  pattern_stratum <- factor(
    as.character(pattern_stratum),
    levels = c("common", "rare", "unobserved")
  )

  if (length(pattern_stratum) != length(y_true)) {
    stop("pattern_stratum and y_true must have the same length.")
  }

  out <- list(
    overall = summarize_predictive_samples(
      draw_mat = draw_mat,
      y_true = y_true,
      method = method,
      rep_id = rep_id,
      n_calib = n_calib,
      scenario = scenario,
      evaluation_stratum = "overall"
    )
  )

  for (stratum in levels(pattern_stratum)) {
    idx <- which(pattern_stratum == stratum)
    if (length(idx) == 0L) next

    out[[stratum]] <- summarize_predictive_samples(
      draw_mat = subset_predictive_samples(draw_mat, idx),
      y_true = y_true[idx],
      method = method,
      rep_id = rep_id,
      n_calib = n_calib,
      scenario = scenario,
      evaluation_stratum = stratum
    )
  }

  do.call(rbind, out)
}


study2_latent_imputation_metrics <- function(fit,
                                             U_true,
                                             rep_id,
                                             n_calib,
                                             scenario,
                                             method = "EIV-GP") {
  miss_idx <- fit$data$miss_idx
  status <- study2_latent_imputation_status(
    fit = fit,
    method = method,
    target = "training_missing_U",
    rep_id = rep_id,
    n_calib = n_calib,
    scenario = scenario,
    n_units = length(miss_idx)
  )
  if (!isTRUE(status$task_eligible)) return(data.frame())

  samples <- fit$mcmc$samples_U[, miss_idx, , drop = FALSE]
  U_center <- fit$data$U_center
  U_scale <- fit$data$U_scale
  if (!is.null(U_center) && !is.null(U_scale)) {
    if (length(U_center) != dim(samples)[3L] ||
        length(U_scale) != dim(samples)[3L]) {
      stop("Stored latent standardization parameters have the wrong length.")
    }
    for (j in seq_len(dim(samples)[3L])) {
      samples[, , j] <- U_center[j] + U_scale[j] * samples[, , j]
    }
  }
  truth <- as.matrix(U_true)[miss_idx, , drop = FALSE]
  summarize_study2_latent_draws(
    draws = samples,
    U_true = truth,
    method = method,
    target = "training_missing_U",
    rep_id = rep_id,
    n_calib = n_calib,
    scenario = scenario
  )
}


study2_eiv_prospective_imputation_metrics <- function(
    fit,
    C_new,
    U_true,
    rep_id,
    n_calib,
    scenario,
    max_draw = 200L,
    latent_sampler = c("minimax_tilting", "rejection", "gibbs"),
    n_new_latent_gibbs = 100L,
    rejection_batch_size = NULL,
    rejection_max_batches = 1000L,
    seed = NULL) {
  latent_sampler <- match.arg(latent_sampler)
  U_true <- as.matrix(U_true)
  if (nrow(as.matrix(C_new)) != nrow(U_true)) {
    stop("C_new and U_true must contain the same prospective units.")
  }
  status <- study2_latent_imputation_status(
    fit = fit,
    method = "EIV-GP",
    target = "prospective_U_given_C",
    rep_id = rep_id,
    n_calib = n_calib,
    scenario = scenario,
    n_units = nrow(U_true)
  )
  if (!isTRUE(status$task_eligible)) return(data.frame())

  n_available <- dim(fit$mcmc$samples_U)[1L]
  max_draw <- as.integer(max_draw)
  if (length(max_draw) != 1L || is.na(max_draw) || max_draw < 1L) {
    stop("max_draw must be one positive integer.")
  }
  draw_ids <- if (n_available <= max_draw) {
    seq_len(n_available)
  } else {
    unique(round(seq(1, n_available, length.out = max_draw)))
  }
  draws <- sample_eiv_u_given_c_ordprobit(
    C_new = C_new,
    fit_obj = fit,
    draw_ids = draw_ids,
    n_per_draw = 1L,
    scale = "raw",
    latent_sampler = latent_sampler,
    n_new_latent_gibbs = n_new_latent_gibbs,
    rejection_batch_size = rejection_batch_size,
    rejection_max_batches = rejection_max_batches,
    seed = seed
  )
  summarize_study2_latent_draws(
    draws = draws,
    U_true = U_true,
    method = "EIV-GP",
    target = "prospective_U_given_C",
    rep_id = rep_id,
    n_calib = n_calib,
    scenario = scenario
  )
}


study2_eiv_surface_recovery_metrics <- function(
    fit,
    scenario,
    rep_id,
    n_calib,
    grid_n = 31L,
    max_draw = 200L,
    u_lim = c(-2.5, 2.5),
    seed = NULL) {
  if (n_calib == 0L) return(data.frame())
  if (!is.null(seed)) set.seed(as.integer(seed))

  grid <- expand.grid(
    u1 = seq(u_lim[1], u_lim[2], length.out = grid_n),
    u2 = seq(u_lim[1], u_lim[2], length.out = grid_n)
  )
  U_star <- as.matrix(grid)
  X_star_raw <- matrix(0, nrow(grid), fit$data$p)
  X_star <- sweep(
    sweep(X_star_raw, 2, fit$data$X_center, "-"),
    2,
    fit$data$X_scale,
    "/"
  )

  draw_ids <- seq_len(dim(fit$mcmc$samples_U)[1])
  if (length(draw_ids) > max_draw) {
    draw_ids <- draw_ids[
      unique(round(seq(1, length(draw_ids), length.out = max_draw)))
    ]
  }

  f_draws <- matrix(NA_real_, length(draw_ids), nrow(grid))
  conditional_means <- matrix(NA_real_, length(draw_ids), nrow(grid))
  kernel_spec <- kernel_spec_from_fit(fit)
  for (ii in seq_along(draw_ids)) {
    s <- draw_ids[ii]
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

    f_draws[ii, ] <-
      fit$data$y_center + fit$data$y_scale *
      (as.numeric(pred$mean) + sqrt(pred$var) * rnorm(nrow(grid)))
    conditional_means[ii, ] <-
      fit$data$y_center + fit$data$y_scale * as.numeric(pred$mean)
  }

  truth <- f0_2d(X_star_raw, U_star, scenario = scenario)
  post_mean <- colMeans(conditional_means)
  lo <- apply(f_draws, 2, quantile, probs = 0.025, names = FALSE)
  hi <- apply(f_draws, 2, quantile, probs = 0.975, names = FALSE)

  data.frame(
    rep = rep_id,
    scenario = scenario,
    n_calib = n_calib,
    method = "EIV-GP",
    grid_n = nrow(grid),
    ISE = mean((post_mean - truth)^2),
    Bias = mean(post_mean - truth),
    Coverage95 = mean(truth >= lo & truth <= hi),
    Width95 = mean(hi - lo),
    stringsAsFactors = FALSE
  )
}


validate_mixedgp_inputs <- function(X, C, m_vec = NULL) {
  X <- as.matrix(X)
  C <- as.matrix(C)
  storage.mode(X) <- "double"
  storage.mode(C) <- "integer"

  if (nrow(X) != nrow(C) || ncol(X) < 1L || ncol(C) < 1L) {
    stop("X and C must have the same number of rows and at least one column.")
  }
  if (any(!is.finite(X)) || any(!is.finite(C))) {
    stop("Competitor inputs must be finite.")
  }

  if (is.null(m_vec)) m_vec <- apply(C, 2L, max)
  m_vec <- as.integer(m_vec)
  if (length(m_vec) != ncol(C) || any(m_vec < 2L)) {
    stop("m_vec must give at least two ordered levels for every column of C.")
  }
  for (j in seq_len(ncol(C))) {
    if (any(C[, j] < 1L | C[, j] > m_vec[j])) {
      stop("Column ", j, " of C contains a value outside 1:m_vec[j].")
    }
  }

  list(X = X, C = C, m_vec = m_vec)
}


predictive_variance_diagonal <- function(x,
                                         n_test,
                                         context = "package prediction") {
  if (is.null(x)) stop("The package prediction did not return uncertainty.")
  if (is.matrix(x)) {
    if (!identical(dim(x), c(as.integer(n_test), as.integer(n_test)))) {
      stop("Unexpected predictive covariance dimension.")
    }
    x <- diag(x)
  }
  x <- as.numeric(x)
  if (length(x) != n_test || any(!is.finite(x))) {
    stop("Predictive variances must be finite with length n_test.")
  }
  tolerance <- 100 * (as.integer(n_test) + 1L) * .Machine$double.eps *
    max(1, max(abs(x)))
  if (any(x < -tolerance)) {
    stop(context, " returned a materially negative predictive variance.")
  }
  pmax(x, 0)
}


sample_independent_predictive_marginals <- function(mean, variance, n_draw, seed) {
  mean <- as.numeric(mean)
  variance <- as.numeric(variance)
  if (length(mean) != length(variance) || any(!is.finite(mean)) ||
      any(!is.finite(variance)) || any(variance < 0)) {
    stop("Invalid predictive means or variances.")
  }
  set.seed(seed)
  matrix(
    stats::rnorm(
      n_draw * length(mean),
      mean = rep(mean, each = n_draw),
      sd = rep(sqrt(variance), each = n_draw)
    ),
    nrow = n_draw,
    ncol = length(mean)
  )
}


attach_predictive_normal_components <- function(draws, mean, variance) {
  draws <- as.matrix(draws)
  mean <- as.numeric(mean)
  variance <- as.numeric(variance)
  if (length(mean) != ncol(draws) || length(variance) != ncol(draws) ||
      any(!is.finite(mean)) || any(!is.finite(variance)) ||
      any(variance < 0)) {
    stop(
      "A Gaussian predictive distribution requires one finite mean and one ",
      "nonnegative variance per test case."
    )
  }
  ## Each package adapter supplies one fitted Gaussian predictive marginal per
  ## test case.  A single row is therefore the exact mixture representation;
  ## it must not be repeated n_draw times when computing log density.
  attr(draws, "conditional_means") <- matrix(mean, nrow = 1L)
  attr(draws, "conditional_vars") <- matrix(variance, nrow = 1L)
  attr(draws, "mixture_components") <- "public-package Gaussian predictor"
  draws
}


validate_mixedgp_competitor_result <- function(result, n_draw, n_test) {
  if (!is.list(result) || is.null(result$draws)) {
    stop("A competitor adapter must return a list containing 'draws'.")
  }
  component_means <- attr(result$draws, "conditional_means", exact = TRUE)
  component_vars <- attr(result$draws, "conditional_vars", exact = TRUE)
  result$draws <- as.matrix(result$draws)
  if (!is.null(component_means)) {
    attr(result$draws, "conditional_means") <- component_means
  }
  if (!is.null(component_vars)) {
    attr(result$draws, "conditional_vars") <- component_vars
  }
  expected <- c(as.integer(n_draw), as.integer(n_test))
  if (!identical(dim(result$draws), expected)) {
    stop(
      "Competitor predictive draws must have dimension n_draw by n_test; got ",
      paste(dim(result$draws), collapse = " by "), "."
    )
  }
  if (any(!is.finite(result$draws))) {
    stop("Competitor adapter returned non-finite predictive draws.")
  }

  ## Adapters in this file return exact package point predictors.  The
  ## fallbacks preserve compatibility with older/custom adapter results that
  ## contained only noisy-response draws; they are not used by the audited
  ## publication adapters above.
  if (is.null(result$latent_mean)) {
    result$latent_mean <- result$predictive_mean
  }
  if (is.null(result$latent_mean)) {
    result$latent_mean <- colMeans(result$draws)
    result$mean_source <- "monte_carlo_draw_average_compatibility_fallback"
  }
  result$latent_mean <- as.numeric(result$latent_mean)
  if (length(result$latent_mean) != n_test ||
      any(!is.finite(result$latent_mean))) {
    stop("Competitor latent mean must be finite with length n_test.")
  }
  if (is.null(result$predictive_mean)) {
    result$predictive_mean <- result$latent_mean
  }
  result$predictive_mean <- as.numeric(result$predictive_mean)
  if (length(result$predictive_mean) != n_test ||
      any(!is.finite(result$predictive_mean))) {
    stop("Competitor predictive mean must be finite with length n_test.")
  }
  if (is.null(result$predictive_variance) && !is.null(component_vars) &&
      nrow(as.matrix(component_vars)) == 1L) {
    result$predictive_variance <- as.numeric(component_vars)
  }
  if (!is.null(result$predictive_variance)) {
    result$predictive_variance <- as.numeric(result$predictive_variance)
    if (length(result$predictive_variance) != n_test ||
        any(!is.finite(result$predictive_variance)) ||
        any(result$predictive_variance < 0)) {
      stop("Competitor predictive variance must be nonnegative with length n_test.")
    }
    if (is.null(component_means) || is.null(component_vars)) {
      result$draws <- attach_predictive_normal_components(
        result$draws,
        result$predictive_mean,
        result$predictive_variance
      )
    }
  }
  result$mean_source <- result$mean_source %||% "public_package_point_predictor"
  result
}


validate_study2_competitor_result <- validate_mixedgp_competitor_result


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
#' @details Numerical and application workflows use a shared core pool. Within
#'   that pool, prediction stream schema 1 assigns seeds to fixed blocks of eight
#'   selected states, independent of core count. Joint prediction rows are never
#'   split. The matrix records `prediction_stream_schema = 1L`. These Monte Carlo
#'   realizations can differ from calls outside the pool and older serial runs;
#'   fitted MCMC draws are unchanged. Workflow cache identities remain enforced.
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
  if (!is.null(mixedgp_budget()) && !isTRUE(getOption("eivGP.prediction_block",FALSE))) {
    call <- match.call(expand.dots=TRUE)
    args <- lapply(as.list(call)[-1L],eval,envir=parent.frame())
    return(mixedgp_predict_blocks(args))
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


# ---- core/00_fit_storage.R ----
# Lossless repository artifact storage. Keep one chain history; rebuild pooled
# views on read. This is not thinning or a change to the sampler/checkpoint.
mixedgp_storage_join <- function(parts) {
  if (is.null(parts[[1L]])) return(NULL)
  dims <- dim(parts[[1L]])
  if (length(dims) == 3L) {
    out <- array(NA_real_, c(sum(vapply(parts, function(x) dim(x)[1L], integer(1))), dims[-1L]))
    pos <- 0L
    for (x in parts) { n <- dim(x)[1L]; out[pos + seq_len(n), , ] <- x; pos <- pos + n }
    return(out)
  }
  if (is.matrix(parts[[1L]])) return(do.call(rbind, parts))
  unlist(parts, use.names = FALSE)
}

mixedgp_pack_artifact <- function(x) {
  if (inherits(x, "mixedgp_compact_fit")) return(x)
  if (!is.list(x) || is.data.frame(x)) return(x)
  if (is.list(x$mcmc) && length(x$mcmc$samples_by_chain)) {
    original_names <- names(x$mcmc)
    removed <- character()
    for (field in names(x$mcmc$samples_by_chain)) {
      name <- paste0("samples_", field)
      parts <- x$mcmc$samples_by_chain[[field]]
      # Preserve unusual/legacy views unless reconstruction is exactly equal.
      if (name %in% original_names && length(parts) &&
          identical(x$mcmc[[name]], mixedgp_storage_join(parts))) {
        x$mcmc[name] <- NULL
        removed <- c(removed, field)
      }
    }
    return(structure(list(version = 1L, fit = x, pooled_fields = removed,
      mcmc_names = original_names), class = "mixedgp_compact_fit"))
  }
  # Traverse ordinary containers, preserving names/classes/NULL entries.
  for (i in seq_along(x)) x[i] <- list(mixedgp_pack_artifact(x[[i]]))
  x
}

mixedgp_unpack_artifact <- function(x) {
  if (inherits(x, "mixedgp_compact_fit")) {
    if (!identical(x$version, 1L)) stop("Unsupported compact-fit version.")
    fit <- x$fit
    for (field in x$pooled_fields) {
      parts <- fit$mcmc$samples_by_chain[[field]]
      if (!length(parts)) stop("Compact fit is missing chain history: ", field)
      fit$mcmc[paste0("samples_", field)] <- list(mixedgp_storage_join(parts))
    }
    fit$mcmc <- fit$mcmc[x$mcmc_names]
    return(fit)
  }
  if (!is.list(x) || is.data.frame(x)) return(x)
  for (i in seq_along(x)) x[i] <- list(mixedgp_unpack_artifact(x[[i]]))
  x
}

mixedgp_read_artifact <- function(path) mixedgp_unpack_artifact(readRDS(path))

mixedgp_save_artifact <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(".fit-storage-", tmpdir = dirname(path))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(mixedgp_pack_artifact(object), tmp, version = 3L)
  # Never remove a valid prior checkpoint before a successful replacement.
  if (!file.rename(tmp, path)) stop("Cannot publish artifact; prior file retained: ", path)
  mixedgp_retire_saved_checkpoints(object, path)
  invisible(path)
}


# ---- simulations/00_synthetic_data.R ----
############################################################
## Frozen synthetic datasets for reproducible experiments
############################################################

MIXEDGP_DATA_SCHEMA_VERSION <- "1.1.0"
MIXEDGP_STUDY1_GENERATOR_TAG <- "study1-synthetic-heterogeneity-v2"
MIXEDGP_STUDY2_GENERATOR_TAG <- "study2-synthetic-shared-state-v2"

mixedgp_object_fingerprint <- function(object) {
  path <- tempfile("mixedgp-object-", fileext = ".rds")
  on.exit(if (file.exists(path)) unlink(path), add = TRUE)
  saveRDS(object, path, version = 3L, compress = FALSE)
  unname(tools::md5sum(path))
}

mixedgp_expected_generator_tag <- function(study) {
  study <- match.arg(study, c("study1", "study2"))
  if (study == "study1") return(MIXEDGP_STUDY1_GENERATOR_TAG)
  MIXEDGP_STUDY2_GENERATOR_TAG
}

mixedgp_synthetic_design_fingerprint <- function(study,
                                                  scenario,
                                                  design,
                                                  generator_tag) {
  payload <- list(
    schema_version = MIXEDGP_DATA_SCHEMA_VERSION,
    study = as.character(study),
    scenario = as.character(scenario),
    generator_tag = as.character(generator_tag),
    design = design
  )
  path <- tempfile("mixedgp-design-", fileext = ".rds")
  on.exit(if (file.exists(path)) unlink(path), add = TRUE)
  saveRDS(payload, path, version = 3L, compress = FALSE)
  unname(tools::md5sum(path))
}

mixedgp_dataset_filename <- function(study,
                                     rep_id,
                                     scenario,
                                     n,
                                     n_test,
                                     m,
                                     q = NULL) {
  study <- match.arg(study, c("study1", "study2"))
  scenario <- gsub("[^[:alnum:]_]+", "-", as.character(scenario))
  q_tag <- if (is.null(q)) "" else paste0("_q", as.integer(q))
  sprintf(
    "%s_%s%s_n%d_ntest%d_m%d_rep%03d.rds",
    study, scenario, q_tag, as.integer(n), as.integer(n_test),
    as.integer(m), as.integer(rep_id)
  )
}

make_study1_synthetic_dataset <- function(
    rep_id,
    n = 100L,
    n_test = 500L,
    m = 6L,
    scenario = "active",
    sigma_eps = 0.1,
    threshold_design = "imbalanced",
    min_class_count = NULL,
    heterogeneity_eta = 1,
    calib_grid = c(0L, 5L, 10L, 20L, 50L),
    data_seed_base = 100000L,
    calibration_seed_base = 200000L) {
  rep_id <- mixedgp_as_integer_strict(
    rep_id, "rep_id", min_value = 1L, length_expected = 1L
  )
  n <- mixedgp_as_integer_strict(n, "n", min_value = 1L, length_expected = 1L)
  n_test <- mixedgp_as_integer_strict(
    n_test, "n_test", min_value = 1L, length_expected = 1L
  )
  m <- mixedgp_as_integer_strict(m, "m", min_value = 2L, length_expected = 1L)
  calib_grid <- mixedgp_as_integer_strict(
    calib_grid, "calib_grid", min_value = 0L
  )
  if (any(calib_grid > n)) stop("calib_grid cannot exceed n.")
  data_seed_base <- mixedgp_as_integer_strict(
    data_seed_base, "data_seed_base", min_value = 0L, length_expected = 1L
  )
  calibration_seed_base <- mixedgp_as_integer_strict(
    calibration_seed_base, "calibration_seed_base",
    min_value = 0L, length_expected = 1L
  )
  data_seed <- mixedgp_as_integer_strict(
    as.numeric(data_seed_base) + rep_id, "data seed", min_value = 0L,
    length_expected = 1L
  )
  calibration_seed <- mixedgp_as_integer_strict(
    as.numeric(calibration_seed_base) + rep_id, "calibration seed",
    min_value = 0L, length_expected = 1L
  )
  dat <- simulate_1d_data(
    n = n,
    n_test = n_test,
    m = m,
    scenario = scenario,
    sigma_eps = sigma_eps,
    seed = data_seed,
    threshold_design = threshold_design,
    min_class_count = min_class_count,
    heterogeneity_eta = heterogeneity_eta
  )
  calibration_sets <- make_nested_calibration_sets(
    n = n,
    calib_grid = as.integer(calib_grid),
    seed = calibration_seed
  )
  design <- list(
    n = as.integer(n), n_test = as.integer(n_test), m = as.integer(m),
    sigma_eps = sigma_eps, threshold_design = threshold_design,
    min_class_count = dat$min_class_count,
    heterogeneity_eta = heterogeneity_eta,
    calib_grid = as.integer(calib_grid),
    generator_tag = mixedgp_expected_generator_tag("study1")
  )
  design$design_fingerprint <- mixedgp_synthetic_design_fingerprint(
    "study1", scenario, design, design$generator_tag
  )
  structure(
    list(
      schema_version = MIXEDGP_DATA_SCHEMA_VERSION,
      study = "study1",
      rep_id = rep_id,
      scenario = scenario,
      seeds = list(data = data_seed, calibration = calibration_seed),
      design = design,
      data = dat,
      calibration_sets = calibration_sets
    ),
    class = c("mixedgp_synthetic_dataset", "list")
  )
}

make_study2_synthetic_dataset <- function(
    rep_id,
    scenario = "primary",
    n = 120L,
    n_test = 500L,
    q = 4L,
    m = 4L,
    sigma_eps = NULL,
    calib_grid = c(0L, 10L, 25L, 50L, 80L),
    data_seed_base = 1000000L,
    calibration_seed_base = NULL) {
  rep_id <- mixedgp_as_integer_strict(
    rep_id, "rep_id", min_value = 1L, length_expected = 1L
  )
  n <- mixedgp_as_integer_strict(n, "n", min_value = 1L, length_expected = 1L)
  n_test <- mixedgp_as_integer_strict(
    n_test, "n_test", min_value = 1L, length_expected = 1L
  )
  q <- mixedgp_as_integer_strict(q, "q", min_value = 2L, length_expected = 1L)
  m <- mixedgp_as_integer_strict(m, "m", min_value = 2L, length_expected = 1L)
  calib_grid <- mixedgp_as_integer_strict(
    calib_grid, "calib_grid", min_value = 0L
  )
  if (any(calib_grid > n)) stop("calib_grid cannot exceed n.")
  data_seed_base <- mixedgp_as_integer_strict(
    data_seed_base, "data_seed_base", min_value = 0L, length_expected = 1L
  )
  data_seed_base_rep <- mixedgp_as_integer_strict(
    as.numeric(data_seed_base) + 10000 * rep_id,
    "derived Study II data seed base", min_value = 0L, length_expected = 1L
  )
  if (data_seed_base_rep > .Machine$integer.max - 8L) {
    stop("Derived Study II seed base leaves no room for component seeds.")
  }
  data_seed <- data_seed_base_rep + 1L
  if (is.null(calibration_seed_base)) {
    calibration_seed <- data_seed_base_rep + 2L
  } else {
    calibration_seed_base <- mixedgp_as_integer_strict(
      calibration_seed_base, "calibration_seed_base",
      min_value = 0L, length_expected = 1L
    )
    calibration_seed <- mixedgp_as_integer_strict(
      as.numeric(calibration_seed_base) + rep_id, "calibration seed",
      min_value = 0L, length_expected = 1L
    )
  }
  dat <- simulate_study2_data(
    n = n,
    n_test = n_test,
    seed = data_seed,
    scenario = scenario,
    sigma_eps = sigma_eps,
    q = q,
    m = m
  )
  calibration_sets <- make_stratified_calibration_sets_2d(
    C = dat$train$C,
    calib_grid = as.integer(calib_grid),
    seed = calibration_seed
  )
  generator_tag <- mixedgp_expected_generator_tag("study2")
  design <- list(
    n = as.integer(n), n_test = as.integer(n_test), q = as.integer(q),
    m = as.integer(m), d = as.integer(dat$true_params$d),
    sigma_eps = dat$sigma_eps,
    calib_grid = as.integer(calib_grid), generator_tag = generator_tag
  )
  design$design_fingerprint <- mixedgp_synthetic_design_fingerprint(
    "study2", scenario, design, generator_tag
  )
  structure(
    list(
      schema_version = MIXEDGP_DATA_SCHEMA_VERSION,
      study = "study2",
      rep_id = rep_id,
      scenario = scenario,
      seeds = list(data = data_seed, calibration = calibration_seed),
      design = design,
      data = dat,
      calibration_sets = calibration_sets
    ),
    class = c("mixedgp_synthetic_dataset", "list")
  )
}

validate_mixedgp_synthetic_dataset <- function(x) {
  if (!inherits(x, "mixedgp_synthetic_dataset")) {
    stop("Object is not a mixedgp_synthetic_dataset.")
  }
  required <- c(
    "schema_version", "study", "rep_id", "scenario", "seeds", "design",
    "data", "calibration_sets"
  )
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop("Synthetic dataset is missing: ", paste(missing, collapse = ", "))
  }
  if (!identical(x$schema_version, MIXEDGP_DATA_SCHEMA_VERSION)) {
    stop(
      "Unsupported synthetic-data schema ", x$schema_version,
      "; expected ", MIXEDGP_DATA_SCHEMA_VERSION, "."
    )
  }
  if (length(x$study) != 1L || !x$study %in% c("study1", "study2") ||
      length(x$scenario) != 1L || is.na(x$scenario) || !nzchar(x$scenario)) {
    stop("Synthetic dataset has invalid study, scenario, or replication metadata.")
  }
  rep_id <- mixedgp_as_integer_strict(
    x$rep_id, "rep_id", min_value = 1L, length_expected = 1L
  )
  if (!is.list(x$seeds) ||
      !identical(names(x$seeds), c("data", "calibration"))) {
    stop("Synthetic dataset seeds must contain data and calibration entries.")
  }
  invisible(mixedgp_as_integer_strict(
    unlist(x$seeds, use.names = FALSE), "synthetic seeds",
    min_value = 0L, length_expected = 2L
  ))
  design_required <- c(
    "n", "n_test", "m", "calib_grid", "generator_tag",
    "design_fingerprint"
  )
  if (x$study == "study1") {
    design_required <- c(
      design_required, "sigma_eps", "threshold_design", "min_class_count",
      "heterogeneity_eta"
    )
  } else {
    design_required <- c(design_required, "q", "d", "sigma_eps")
  }
  missing_design <- setdiff(design_required, names(x$design))
  if (length(missing_design) > 0L) {
    stop("Synthetic design is missing: ", paste(missing_design, collapse = ", "))
  }
  generator_tag <- x$design$generator_tag
  design_fingerprint <- x$design$design_fingerprint
  design_without_fingerprint <- x$design[
    setdiff(names(x$design), "design_fingerprint")
  ]
  expected_fingerprint <- mixedgp_synthetic_design_fingerprint(
    x$study, x$scenario, design_without_fingerprint, generator_tag
  )
  if (length(generator_tag) != 1L || is.na(generator_tag) ||
      !nzchar(generator_tag) || length(design_fingerprint) != 1L ||
      is.na(design_fingerprint) ||
      !identical(
        as.character(generator_tag), mixedgp_expected_generator_tag(x$study)
      ) ||
      !identical(as.character(design_fingerprint), expected_fingerprint)) {
    stop("Synthetic generator tag or design fingerprint is invalid.")
  }
  n <- mixedgp_as_integer_strict(
    x$design$n, "design$n", min_value = 1L, length_expected = 1L
  )
  n_test <- mixedgp_as_integer_strict(
    x$design$n_test, "design$n_test", min_value = 1L, length_expected = 1L
  )
  m <- mixedgp_as_integer_strict(
    x$design$m, "design$m", min_value = 2L, length_expected = 1L
  )
  calib_grid <- mixedgp_as_integer_strict(
    x$design$calib_grid, "design$calib_grid", min_value = 0L
  )
  if (any(calib_grid > n)) {
    stop("Synthetic dataset calibration sizes cannot exceed n.")
  }
  if (!is.numeric(x$design$sigma_eps) || length(x$design$sigma_eps) != 1L ||
      !is.finite(x$design$sigma_eps) || x$design$sigma_eps < 0) {
    stop("design$sigma_eps must be one nonnegative finite number.")
  }
  if (!is.list(x$calibration_sets) ||
      !identical(names(x$calibration_sets), as.character(calib_grid)) ||
      length(x$calibration_sets) != length(calib_grid)) {
    stop("Calibration sets do not match the declared calibration grid.")
  }
  valid_calibration <- vapply(seq_along(calib_grid), function(ii) {
    idx <- x$calibration_sets[[ii]]
    if (calib_grid[ii] == 0L) return(length(idx) == 0L)
    valid_idx <- try(
      mixedgp_as_integer_strict(idx, "calibration row indices", min_value = 1L),
      silent = TRUE
    )
    !inherits(valid_idx, "try-error") &&
      identical(length(valid_idx), calib_grid[ii]) &&
      length(unique(valid_idx)) == length(valid_idx) && all(valid_idx <= n)
  }, logical(1L))
  if (!all(valid_calibration)) {
    stop("A calibration set has invalid size, duplicates, or row indices.")
  }
  grid_order <- order(calib_grid)
  if (length(grid_order) > 1L) {
    nested <- vapply(seq_len(length(grid_order) - 1L), function(ii) {
      all(x$calibration_sets[[grid_order[ii]]] %in%
            x$calibration_sets[[grid_order[ii + 1L]]])
    }, logical(1L))
    if (!all(nested)) stop("Calibration sets must be nested by size.")
  }

  if (x$study == "study1") {
    required_values <- c("x", "u", "c", "y", "f")
    numeric_train <- is.data.frame(x$data$train) &&
      all(required_values %in% names(x$data$train)) &&
      all(vapply(x$data$train[required_values], is.numeric, logical(1L)))
    numeric_test <- is.data.frame(x$data$test) &&
      all(required_values %in% names(x$data$test)) &&
      all(vapply(x$data$test[required_values], is.numeric, logical(1L)))
    if (!numeric_train || !numeric_test) {
      stop("Study I train/test truth columns must be numeric.")
    }
    if (!is.numeric(x$data$tau_true) || length(x$data$tau_true) != m - 1L ||
        any(!is.finite(x$data$tau_true)) ||
        is.unsorted(x$data$tau_true, strictly = TRUE) ||
        !identical(as.character(x$data$scenario), as.character(x$scenario)) ||
        !identical(as.integer(x$data$m), m) ||
        !isTRUE(all.equal(
          as.numeric(x$data$sigma_eps), as.numeric(x$design$sigma_eps),
          tolerance = 0
        ))) {
      stop("Study I DGP metadata or cutpoints are inconsistent with the design.")
    }
    if (!is.data.frame(x$data$train) || !is.data.frame(x$data$test) ||
        nrow(x$data$train) != n || nrow(x$data$test) != n_test ||
        anyNA(x$data$train[, required_values]) ||
        anyNA(x$data$test[, required_values]) ||
        any(!is.finite(as.matrix(x$data$train[, required_values]))) ||
        any(!is.finite(as.matrix(x$data$test[, required_values]))) ||
        any(x$data$train$c != floor(x$data$train$c)) ||
        any(x$data$test$c != floor(x$data$test$c)) ||
        any(x$data$train$c < 1L | x$data$train$c > m) ||
        any(x$data$test$c < 1L | x$data$test$c > m)) {
      stop("Study I data dimensions or values do not match the design.")
    }
    expected_train_f <- f0_1d(
      x$data$train$x, x$data$train$u, scenario = x$scenario,
      c_ord = x$data$train$c, tau = x$data$tau_true,
      heterogeneity_eta = x$design$heterogeneity_eta
    )
    expected_test_f <- f0_1d(
      x$data$test$x, x$data$test$u, scenario = x$scenario,
      c_ord = x$data$test$c, tau = x$data$tau_true,
      heterogeneity_eta = x$design$heterogeneity_eta
    )
    if (!isTRUE(all.equal(x$data$train$f, expected_train_f, tolerance = 0)) ||
        !isTRUE(all.equal(x$data$test$f, expected_test_f, tolerance = 0))) {
      stop("Study I stored response-surface truth is inconsistent with the DGP.")
    }
  } else {
    q <- mixedgp_as_integer_strict(
      x$design$q, "design$q", min_value = 2L, length_expected = 1L
    )
    d <- mixedgp_as_integer_strict(
      x$design$d, "design$d", min_value = 1L, length_expected = 1L
    )
    train_U <- as.matrix(x$data$train$U)
    test_U <- as.matrix(x$data$test$U)
    train_f <- x$data$train$f
    test_f <- x$data$test$f
    true_params <- x$data$true_params
    numeric_objects <- list(
      x$data$train$X, x$data$test$X, train_U, test_U,
      x$data$train$C, x$data$test$C, x$data$train$y,
      x$data$test$y, train_f, test_f
    )
    if (any(!vapply(
      numeric_objects,
      function(value) is.numeric(value) || is.integer(value),
      logical(1L)
    ))) {
      stop("Study II data and truth arrays must be numeric.")
    }
    true_required <- c(
      "scenario", "score_error", "response_interactions", "lambda", "A",
      "Omega", "tau", "q", "m", "d", "sigma_eps"
    )
    if (!is.list(true_params) ||
        length(setdiff(true_required, names(true_params))) > 0L ||
        any(!vapply(
          true_params[c("lambda", "A", "Omega", "tau", "q", "m", "d",
                        "sigma_eps")],
          function(value) is.numeric(value) || is.integer(value),
          logical(1L)
        ))) {
      stop("Study II true_params is incomplete or nonnumeric.")
    }
    if (length(q) != 1L || is.na(q) || q < 2L ||
        length(d) != 1L || is.na(d) || d < 1L ||
        !is.list(x$data$train) || !is.list(x$data$test) ||
        nrow(as.matrix(x$data$train$X)) != n ||
        nrow(as.matrix(x$data$test$X)) != n_test ||
        !identical(dim(as.matrix(x$data$train$C)), c(n, q)) ||
        !identical(dim(as.matrix(x$data$test$C)), c(n_test, q)) ||
        !identical(dim(train_U), c(n, d)) ||
        !identical(dim(test_U), c(n_test, d)) ||
        length(x$data$train$y) != n || length(x$data$test$y) != n_test ||
        length(train_f) != n || length(test_f) != n_test ||
        anyNA(x$data$train$X) || anyNA(x$data$test$X) ||
        anyNA(train_U) || anyNA(test_U) ||
        anyNA(x$data$train$C) || anyNA(x$data$test$C) ||
        anyNA(x$data$train$y) || anyNA(x$data$test$y) ||
        anyNA(train_f) || anyNA(test_f) ||
        any(x$data$train$C != floor(x$data$train$C)) ||
        any(x$data$test$C != floor(x$data$test$C)) ||
        any(!is.finite(c(
          x$data$train$X, x$data$test$X, train_U, test_U,
          x$data$train$y, x$data$test$y, train_f, test_f
        ))) ||
        any(x$data$train$C < 1L | x$data$train$C > m) ||
        any(x$data$test$C < 1L | x$data$test$C > m) ||
        !is.list(true_params) || !identical(as.integer(true_params$q), q) ||
        !identical(as.integer(true_params$m), m) ||
        !identical(as.integer(true_params$d), d) ||
        !identical(dim(as.matrix(true_params$A)), c(q, d)) ||
        !identical(dim(as.matrix(true_params$tau)), c(q, m - 1L)) ||
        !identical(dim(as.matrix(true_params$Omega)), c(q, q)) ||
        any(!is.finite(c(
          true_params$A, true_params$tau, true_params$Omega,
          true_params$lambda, true_params$sigma_eps
        ))) ||
        !identical(as.character(true_params$scenario), as.character(x$scenario))) {
      stop("Study II data dimensions or values do not match the design.")
    }
    expected_params <- make_study2_true_params(
      scenario = x$scenario, q = q, m = m
    )
    if (!isTRUE(all.equal(true_params, expected_params, tolerance = 0)) ||
        !identical(as.character(x$data$scenario), as.character(x$scenario)) ||
        !isTRUE(all.equal(
          as.numeric(x$data$sigma_eps), as.numeric(x$design$sigma_eps),
          tolerance = 0
        ))) {
      stop("Study II true parameters or DGP metadata do not match the design.")
    }
    expected_train_f <- f0_2d(
      x$data$train$X, train_U, scenario = x$scenario
    )
    expected_test_f <- f0_2d(x$data$test$X, test_U, scenario = x$scenario)
    if (!isTRUE(all.equal(train_f, expected_train_f, tolerance = 0)) ||
        !isTRUE(all.equal(test_f, expected_test_f, tolerance = 0))) {
      stop("Study II stored response-surface truth is inconsistent with the DGP.")
    }
  }
  invisible(TRUE)
}

store_mixedgp_synthetic_dataset <- function(x,
                                            directory,
                                            overwrite = FALSE,
                                            compress = "xz") {
  overwrite <- mixedgp_validate_flag(overwrite, "overwrite")
  validate_mixedgp_synthetic_dataset(x)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  design <- x$design
  path <- file.path(
    directory,
    mixedgp_dataset_filename(
      study = x$study,
      rep_id = x$rep_id,
      scenario = x$scenario,
      n = design$n,
      n_test = design$n_test,
      m = design$m,
      q = design$q
    )
  )
  if (file.exists(path) && !isTRUE(overwrite)) {
    existing <- readRDS(path)
    validate_mixedgp_synthetic_dataset(existing)
    if (!identical(existing, x)) {
      stop(
        "A different synthetic dataset already exists at ", path, ". ",
        "The filename is not a sufficient design fingerprint; choose a new ",
        "directory or set overwrite = TRUE after verifying the requested design."
      )
    }
  } else {
    temporary <- tempfile(
      pattern = paste0(basename(path), "."),
      tmpdir = directory
    )
    on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
    saveRDS(x, temporary, version = 3L, compress = compress)
    if (!file.rename(temporary, path)) {
      stop("Could not atomically move the synthetic dataset to ", path, ".")
    }
  }
  info <- file.info(path)
  data.frame(
    schema_version = x$schema_version,
    study = x$study,
    scenario = x$scenario,
    rep = x$rep_id,
    n = design$n,
    n_test = design$n_test,
    q = if (is.null(design$q)) NA_integer_ else design$q,
    d = if (is.null(design$d)) NA_integer_ else design$d,
    m = design$m,
    data_seed = x$seeds$data,
    calibration_seed = x$seeds$calibration,
    file = basename(path),
    bytes = as.numeric(info$size),
    md5 = unname(tools::md5sum(path)),
    generator_tag = as.character(design$generator_tag),
    design_fingerprint = as.character(design$design_fingerprint),
    stringsAsFactors = FALSE
  )
}

load_mixedgp_synthetic_dataset <- function(path, expected_md5 = NULL) {
  if (!file.exists(path)) stop("Synthetic dataset does not exist: ", path)
  if (!is.null(expected_md5)) {
    observed <- unname(tools::md5sum(path))
    if (!identical(tolower(observed), tolower(as.character(expected_md5)))) {
      stop("Checksum mismatch for synthetic dataset: ", path)
    }
  }
  out <- readRDS(path)
  validate_mixedgp_synthetic_dataset(out)
  out
}

load_mixedgp_synthetic_dataset_strict <- function(path,
                                                   expected,
                                                   manifest_path = file.path(
                                                     dirname(path),
                                                     "manifest.rds"
                                                   )) {
  if (!is.list(expected) || is.null(names(expected))) {
    stop("expected must be a named list of frozen-design fields.")
  }
  if (!file.exists(manifest_path)) {
    stop("Frozen-data manifest does not exist: ", manifest_path)
  }
  manifest_object <- readRDS(manifest_path)
  manifest <- if (is.data.frame(manifest_object)) {
    manifest_object
  } else if (is.list(manifest_object) &&
             is.data.frame(manifest_object$manifest)) {
    manifest_object$manifest
  } else {
    stop("Frozen-data manifest has an unsupported structure: ", manifest_path)
  }
  required_columns <- c(
    "schema_version", "study", "scenario", "rep", "n", "n_test", "m",
    "file", "md5", "generator_tag", "design_fingerprint"
  )
  missing_columns <- setdiff(required_columns, names(manifest))
  if (length(missing_columns) > 0L) {
    stop(
      "Frozen-data manifest is missing: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  row_id <- which(as.character(manifest$file) == basename(path))
  if (length(row_id) != 1L) {
    stop(
      "Frozen-data manifest must contain exactly one row for ", basename(path),
      "; found ", length(row_id), "."
    )
  }
  manifest_row <- manifest[row_id, , drop = FALSE]
  out <- load_mixedgp_synthetic_dataset(
    path,
    expected_md5 = as.character(manifest_row$md5[[1L]])
  )

  observed <- c(
    list(
      study = out$study,
      scenario = out$scenario,
      rep_id = out$rep_id
    ),
    out$design
  )
  for (field in names(expected)) {
    if (!field %in% names(observed)) {
      stop("Unknown expected frozen-design field: ", field)
    }
    expected_value <- expected[[field]]
    observed_value <- observed[[field]]
    equal <- if (is.numeric(expected_value) && is.numeric(observed_value)) {
      isTRUE(all.equal(
        as.numeric(observed_value), as.numeric(expected_value),
        tolerance = 0, check.attributes = FALSE
      ))
    } else {
      identical(as.character(observed_value), as.character(expected_value))
    }
    if (!equal) {
      stop(
        "Frozen dataset ", basename(path), " has ", field, "=",
        paste(observed_value, collapse = ","), "; expected ",
        paste(expected_value, collapse = ","), "."
      )
    }
  }

  manifest_checks <- list(
    schema_version = out$schema_version,
    study = out$study,
    scenario = out$scenario,
    rep = out$rep_id,
    n = out$design$n,
    n_test = out$design$n_test,
    m = out$design$m,
    generator_tag = out$design$generator_tag,
    design_fingerprint = out$design$design_fingerprint
  )
  if ("q" %in% names(manifest) && !is.null(out$design$q)) {
    manifest_checks$q <- out$design$q
  }
  if ("d" %in% names(manifest) && !is.null(out$design$d)) {
    manifest_checks$d <- out$design$d
  }
  for (field in names(manifest_checks)) {
    left <- manifest_row[[field]][[1L]]
    right <- manifest_checks[[field]]
    if (is.numeric(left) && is.numeric(right)) {
      equal <- isTRUE(all.equal(
        as.numeric(left), as.numeric(right), tolerance = 0,
        check.attributes = FALSE
      ))
    } else {
      equal <- identical(as.character(left), as.character(right))
    }
    if (!equal) {
      stop(
        "Manifest/artifact mismatch for ", basename(path), " field ", field,
        "."
      )
    }
  }
  attr(out, "manifest_path") <- normalizePath(manifest_path)
  attr(out, "manifest_md5") <- as.character(manifest_row$md5[[1L]])
  out
}

mixedgp_validate_generation_dots <- function(dots, forbidden = "rep_id") {
  if (length(dots) > 0L &&
      (is.null(names(dots)) || anyNA(names(dots)) ||
       any(!nzchar(names(dots))) || anyDuplicated(names(dots)))) {
    stop("Every synthetic-generator argument passed through ... must be named uniquely.")
  }
  duplicate <- intersect(names(dots), forbidden)
  if (length(duplicate) > 0L) {
    stop(
      "Do not pass generator-controlled arguments through ...: ",
      paste(duplicate, collapse = ", "), "."
    )
  }
  invisible(dots)
}

mixedgp_stop_on_worker_errors <- function(rows, context) {
  failed <- which(vapply(rows, inherits, logical(1L), what = "try-error"))
  if (length(failed) == 0L) return(invisible(rows))
  details <- vapply(failed, function(ii) {
    condition <- attr(rows[[ii]], "condition")
    message <- if (inherits(condition, "condition")) {
      conditionMessage(condition)
    } else {
      trimws(as.character(rows[[ii]])[1L])
    }
    paste0("task ", ii, ": ", message)
  }, character(1L))
  stop(
    context, " failed in ", length(failed), " worker(s); no manifest was written. ",
    paste(details, collapse = " | ")
  )
}

mixedgp_bind_manifest_rows <- function(rows) {
  if (length(rows) == 0L) return(data.frame())
  if (!all(vapply(rows, is.data.frame, logical(1L)))) {
    stop("Every successful synthetic-data worker must return a manifest row.")
  }
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(row) {
    missing <- setdiff(columns, names(row))
    for (name in missing) row[[name]] <- NA
    row[, columns, drop = FALSE]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

mixedgp_merge_manifests <- function(existing, incoming) {
  manifests <- Filter(function(x) !is.null(x) && nrow(x) > 0L,
                      list(existing, incoming))
  if (length(manifests) == 0L) return(data.frame())
  if (!all(vapply(manifests, is.data.frame, logical(1L)))) {
    stop("Synthetic-data manifests must be data frames.")
  }
  combined <- mixedgp_bind_manifest_rows(manifests)
  if (!"file" %in% names(combined) || anyNA(combined$file) ||
      any(!nzchar(as.character(combined$file)))) {
    stop("Synthetic-data manifest rows require nonempty file names.")
  }
  groups <- split(seq_len(nrow(combined)), as.character(combined$file))
  keep <- integer(length(groups))
  group_id <- 0L
  for (indices in groups) {
    group_id <- group_id + 1L
    reference <- combined[indices[1L], , drop = FALSE]
    rownames(reference) <- NULL
    consistent <- vapply(indices, function(ii) {
      candidate <- combined[ii, , drop = FALSE]
      rownames(candidate) <- NULL
      identical(as.list(candidate), as.list(reference))
    }, logical(1L))
    if (!all(consistent)) {
      stop(
        "Manifest metadata disagree for existing file ",
        as.character(reference$file[[1L]]), "."
      )
    }
    keep[group_id] <- indices[1L]
  }
  out <- combined[keep, , drop = FALSE]
  ordering_fields <- intersect(c("study", "scenario", "rep", "file"), names(out))
  if (length(ordering_fields) > 0L) {
    ordering <- do.call(order, unname(out[ordering_fields]))
    out <- out[ordering, , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

mixedgp_atomic_replace <- function(temporary, target) {
  if (file.rename(temporary, target)) return(invisible(target))
  backup <- NULL
  if (file.exists(target)) {
    backup <- tempfile(paste0(basename(target), ".backup-"), dirname(target))
    if (!file.rename(target, backup)) {
      stop("Could not prepare atomic replacement of ", target, ".")
    }
  }
  installed <- file.rename(temporary, target)
  if (!installed) {
    if (!is.null(backup) && file.exists(backup)) file.rename(backup, target)
    stop("Could not atomically replace ", target, ".")
  }
  if (!is.null(backup) && file.exists(backup)) unlink(backup)
  invisible(target)
}

mixedgp_write_manifest_atomic <- function(manifest, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  csv_target <- file.path(directory, "manifest.csv")
  rds_target <- file.path(directory, "manifest.rds")
  csv_temporary <- tempfile("manifest-", tmpdir = directory, fileext = ".csv")
  rds_temporary <- tempfile("manifest-", tmpdir = directory, fileext = ".rds")
  on.exit({
    if (file.exists(csv_temporary)) unlink(csv_temporary)
    if (file.exists(rds_temporary)) unlink(rds_temporary)
  }, add = TRUE)
  utils::write.csv(manifest, csv_temporary, row.names = FALSE)
  saveRDS(manifest, rds_temporary, version = 3L)
  mixedgp_atomic_replace(rds_temporary, rds_target)
  mixedgp_atomic_replace(csv_temporary, csv_target)
  invisible(manifest)
}

mixedgp_existing_manifest <- function(directory) {
  rds_path <- file.path(directory, "manifest.rds")
  csv_path <- file.path(directory, "manifest.csv")
  if (!file.exists(rds_path)) {
    if (file.exists(csv_path)) {
      stop("manifest.csv exists without the canonical manifest.rds in ", directory)
    }
    return(NULL)
  }
  manifest <- readRDS(rds_path)
  if (!is.data.frame(manifest)) {
    stop("Existing per-study manifest.rds must contain a data frame.")
  }
  manifest
}

generate_study1_synthetic_datasets <- function(
    n_rep = 50L,
    directory = file.path("..", "data-synthetic", "study1"),
    n_cores = NULL,
    overwrite = FALSE,
    ...) {
  n_rep <- mixedgp_as_integer_strict(
    n_rep, "n_rep", min_value = 1L, length_expected = 1L
  )
  overwrite <- mixedgp_validate_flag(overwrite, "overwrite")
  dots <- list(...)
  mixedgp_validate_generation_dots(dots, forbidden = "rep_id")
  rep_ids <- seq_len(n_rep)
  rows <- mixedgp_parallel_mapply(
    function(rep_id) {
      artifact <- do.call(
        make_study1_synthetic_dataset, c(list(rep_id = rep_id), dots)
      )
      store_mixedgp_synthetic_dataset(artifact, directory, overwrite)
    },
    rep_id = rep_ids,
    n_cores = n_cores,
    seeds = 710000L + rep_ids,
    SIMPLIFY = FALSE,
    mc.preschedule = FALSE
  )
  mixedgp_stop_on_worker_errors(rows, "Study I synthetic-data generation")
  incoming <- mixedgp_bind_manifest_rows(rows)
  manifest <- mixedgp_merge_manifests(
    mixedgp_existing_manifest(directory), incoming
  )
  mixedgp_write_manifest_atomic(manifest, directory)
  manifest
}

generate_study2_synthetic_datasets <- function(
    n_rep = 50L,
    scenarios = c(
      "primary", "latent_additive_control", "high_uncertainty",
      "logistic_misspec"
    ),
    directory = file.path("..", "data-synthetic", "study2"),
    n_cores = NULL,
    overwrite = FALSE,
    ...) {
  n_rep <- mixedgp_as_integer_strict(
    n_rep, "n_rep", min_value = 1L, length_expected = 1L
  )
  overwrite <- mixedgp_validate_flag(overwrite, "overwrite")
  if (!is.character(scenarios) || length(scenarios) < 1L ||
      anyNA(scenarios) || any(!nzchar(scenarios)) || anyDuplicated(scenarios)) {
    stop("scenarios must contain one or more unique, nonempty names.")
  }
  dots <- list(...)
  mixedgp_validate_generation_dots(dots, forbidden = c("rep_id", "scenario"))
  grid <- expand.grid(
    rep_id = seq_len(n_rep),
    scenario = scenarios,
    stringsAsFactors = FALSE
  )
  rows <- mixedgp_parallel_mapply(
    function(rep_id, scenario) {
      artifact <- do.call(
        make_study2_synthetic_dataset,
        c(list(rep_id = rep_id, scenario = scenario), dots)
      )
      store_mixedgp_synthetic_dataset(artifact, directory, overwrite)
    },
    rep_id = grid$rep_id,
    scenario = grid$scenario,
    n_cores = n_cores,
    seeds = 720000L + seq_len(nrow(grid)),
    SIMPLIFY = FALSE,
    mc.preschedule = FALSE
  )
  mixedgp_stop_on_worker_errors(rows, "Study II synthetic-data generation")
  incoming <- mixedgp_bind_manifest_rows(rows)
  manifest <- mixedgp_merge_manifests(
    mixedgp_existing_manifest(directory), incoming
  )
  mixedgp_write_manifest_atomic(manifest, directory)
  manifest
}


# ---- simulations/competitor_cache.R ----
## Repository-only competitor protocol/cache. No eivGP MCMC or installation.
mixedgp_competitor_protocol <- function(study) {
  study <- match.arg(study, c("study1", "study2"))
  list(`UC-GP` = list(n_starts = 8L, rescue_starts = c(32L,64L)),
    LVGP = list(n_starts = 8L, max_retries = 3L, max_iter_ini = 100L,
      max_iter_lat = 20L, rescue_iter_ini = 300L, rescue_iter_lat = 100L,
      max_elapsed_seconds = if (study == "study1") 1800 else 3600, parallel = FALSE),
    EzGP = list(tau_fractions = c(1e-6, .0025, .01, .04, .16), cv_folds = 3L, maxeval = 100L, rescue_maxeval = c(300L,1000L)))
}

mixedgp_release_competitor_lock <- function(lock, token) {
  if (!file.exists(file.path(lock,"owner.rds"))) return(invisible(NULL))
  owner <- suppressWarnings(tryCatch(readRDS(file.path(lock,"owner.rds")),error=function(e)NULL))
  if (identical(owner$token,token)) unlink(lock,recursive=TRUE)
  invisible(NULL)
}

mixedgp_competitor_hash <- function(x) {
  f <- tempfile(); on.exit(unlink(f), add = TRUE)
  saveRDS(x, f, version = 3L, compress = FALSE)
  unname(tools::md5sum(f))
}

mixedgp_cached_competitors <- function(X_train, y_train, C_train, X_test, C_test,
    n_draw, seed, m_vec = NULL, methods = c("UC-GP", "LVGP", "EzGP"),
    strict = FALSE, controls = list(), cache_root = Sys.getenv("EIVGP_COMPETITOR_CACHE", file.path(tempdir(),"eivGP-competitor-cache")),
    allow_fit = FALSE, retry_failed = FALSE) {
  checked <- validate_mixedgp_inputs(X_train, C_train, m_vec)
  test <- validate_mixedgp_inputs(X_test, C_test, checked$m_vec)
  n_draw <- mixedgp_as_integer_strict(n_draw, "n_draw", 1L, 1L)
  registry <- mixedgp_competitor_registry()
  if (anyDuplicated(methods) || any(!methods %in% registry$method)) stop("Invalid competitor selection.")
  out <- list(draws=list(), means=list(), latent_means=list(), predictive_means=list(),
              predictive_variances=list(), fits=list(), status=NULL)
  statuses <- list()
  for (method in methods) {
    j <- match(method, registry$method)
    package <- registry$package[j]
    available <- requireNamespace(package, quietly = TRUE)
    identity <- list(schema = "competitor-cache-v1", method = method,
      X = unname(checked$X), y = as.numeric(y_train), C = unname(checked$C),
      X_test = unname(test$X), C_test = unname(test$C), m_vec = checked$m_vec,
      seed = as.integer(seed + j), controls = controls[[method]],
      package_version = if (available) as.character(utils::packageVersion(package)) else "missing",
      R_version = paste(R.version$major, R.version$minor, R.version$platform),
      adapter_md5 = eivgp_module_identity())
    key <- mixedgp_competitor_hash(identity)
    path <- file.path(cache_root, method, paste0(key, ".rds"))
    cached <- if (file.exists(path)) readRDS(path) else NULL
    if (!is.null(cached) && !identical(cached$identity, identity)) stop("Competitor cache identity mismatch: ", path)
    hit <- !is.null(cached)
    need_fit <- is.null(cached) || (retry_failed && !all(cached$result$status$status == "success"))
    if (need_fit && allow_fit && available) {
      message("Competitor ", method, ": fitting ", substr(key,1L,12L))
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      lock <- paste0(path, ".lock")
      if (!dir.create(lock, showWarnings = FALSE)) {
        status <- mixedgp_competitor_status(method,status="unavailable_or_failed",
          optimization_status="cache_locked",message=paste("Competitor cache is locked:",lock,
            "No lock was removed; other methods remain eligible."))
        status$cache_key <- key; status$cache_hit <- FALSE; status$cache_file <- path
        statuses[[method]] <- status
        if (strict) stop(status$message)
        next
      }
      # Release our own lock on interrupts as well as ordinary errors.
      token <- basename(tempfile("owner-"))
      do.call(on.exit,list(substitute(mixedgp_release_competitor_lock(P,T),list(P=lock,T=token)),add=TRUE),envir=environment())
      saveRDS(list(token=token,pid=Sys.getpid(),host=Sys.info()[["nodename"]],created=as.character(Sys.time())),
        file.path(lock,"owner.rds"))
      latest <- if(file.exists(path)) readRDS(path) else NULL
      if (!is.null(latest) && !identical(latest$identity,identity)) stop("Competitor cache identity mismatch: ",path)
      result <- if(!is.null(latest) && all(latest$result$status$status=="success")) latest$result else tryCatch({
        ## The adapter returns a Gaussian predictor. Retain the model and exact
        ## moments; draws can subsequently be generated at any evaluation budget.
        run_published_mixedgp_competitors(checked$X, as.numeric(y_train), checked$C,
          test$X, test$C, n_draw = 1L, seed = as.integer(seed + j - 1L),
          m_vec = checked$m_vec, methods = method, strict = FALSE, controls = controls)
      }, error = function(e) list(status = mixedgp_competitor_status(method,
        status="unavailable_or_failed", optimization_status="error", message=conditionMessage(e))))
      tryCatch({
        if (file.exists(path)) {
          backup <- tempfile(paste0(key,"-previous-"), tmpdir=dirname(path), fileext=".rds")
          if (!file.copy(path, backup)) stop("Cannot preserve previous attempt.")
        }
        cached <- list(identity=identity, result=result, created=as.character(Sys.time()))
        tmp <- tempfile(tmpdir=dirname(path)); saveRDS(cached,tmp,version=3L)
        if (!file.rename(tmp,path)) stop("Cannot save competitor cache: ",path)
      }, finally = mixedgp_release_competitor_lock(lock, token))
      hit <- FALSE
      message("Competitor ", method, ": saved ", paste(result$status$status,collapse=","))
    }
    if (is.null(cached)) {
      status <- mixedgp_competitor_status(method, status="unavailable_or_failed",
        optimization_status=if (available) "missing_cache" else "missing_package",
        message=if (available) "Use the run action to include competitors; the fit stage did not fit them."
          else paste("Install experiment dependencies; missing package", package))
    } else {
      result <- cached$result; status <- result$status
      if (all(status$status == "success")) {
        mu <- result$predictive_means[[method]]
        variance <- result$predictive_variances[[method]]
        draw <- sample_independent_predictive_marginals(mu, variance, n_draw, as.integer(seed+j+10000L))
        out$draws[[method]] <- attach_predictive_normal_components(draw, mu, variance)
        out$means[[method]] <- out$latent_means[[method]] <- result$latent_means[[method]]
        out$predictive_means[[method]] <- mu
        out$predictive_variances[[method]] <- variance
        out$fits[[method]] <- result$fits[[method]]
      }
    }
    status$cache_key <- key; status$cache_hit <- hit
    status$cache_file <- if (file.exists(path)) normalizePath(path) else path
    statuses[[method]] <- status
    if (strict && !all(status$status == "success")) stop(status$message[1L])
  }
  out$status <- if (length(statuses)) do.call(rbind,statuses) else
    mixedgp_competitor_status("UC-GP",status="success")[FALSE,,drop=FALSE]
  out
}


# ---- simulations/study1/04_study1_ablations.R ----
############################################################
## 04_study1_ablations.R
##
## Response-free measurement fit, plug-in GP (PI-GP), and complete-case GP
## (CC-GP) for Study I. These are scientific ablations, not literature
## competitors. The deterministic-threshold Gibbs sampler below targets
## p(tau, U_missing | C, U_calibrated) and never receives y.
############################################################



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
                                                    tau_bounds = c(-8, 8),
                                                    .checkpoint_path = NULL) {
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
  start <- 0L; elapsed_before <- 0; clock_start <- proc.time()[["elapsed"]]
  enabled <- is.finite(getOption("eivGP.checkpoint_seconds",Inf)) ||
    is.finite(getOption("eivGP.checkpoint_iterations",Inf))
  if(enabled && !is.null(.checkpoint_path) && file.exists(.checkpoint_path)) {
    cp<-readRDS(.checkpoint_path)
    if(!identical(cp$schema,1L) || cp$iteration<0L || cp$iteration>n_iter ||
       cp$saved!=sum(keep_iter<=cp$iteration))stop("Invalid Study I measurement checkpoint")
    restored<-mixedgp_checkpoint_restore_draws(cp$draws,list(samples_u=samples_u,samples_tau=samples_tau))
    samples_u<-restored$samples_u;samples_tau<-restored$samples_tau
    start<-cp$iteration;keep_pos<-cp$saved;u<-cp$u;tau<-cp$tau
    elapsed_before<-cp$elapsed_seconds;mixedgp_restore_rng_state(cp$rng)
  }
  last_save<-last_progress<-proc.time()[["elapsed"]];last_iteration<-start
  checkpoint<-function(iter,force=FALSE) {
    if(is.null(.checkpoint_path))return(invisible(NULL))
    now<-proc.time()[["elapsed"]]
    saving<-mixedgp_checkpoint_due(now,last_save,iter-last_iteration,force)
    if(saving) {
      mixedgp_checkpoint_write(list(schema=1L,iteration=iter,saved=keep_pos,u=u,tau=tau,
        draws=mixedgp_checkpoint_prefix(list(samples_u=samples_u,samples_tau=samples_tau),keep_pos),
        rng=mixedgp_rng_state(),elapsed_seconds=elapsed_before+now-clock_start),.checkpoint_path)
      last_save<<-proc.time()[["elapsed"]];last_iteration<<-iter
    }
    if(force || saving || now-last_progress>=60) {
      elapsed<-elapsed_before+proc.time()[["elapsed"]]-clock_start
      remaining<-if(iter>0)elapsed/iter*(n_iter-iter)else NA_real_
      mixedgp_timing_record(sub("[.]rds$","-timing.csv",.checkpoint_path),data.frame(
        updated_at=format(Sys.time(),tz="UTC",usetz=TRUE),iteration=iter,target_iteration=n_iter,
        checkpoint_iteration=last_iteration,phase=if(iter==n_iter)"complete"else if(iter<=burn)"warmup"else"sampling",
        elapsed_seconds=elapsed,remaining_seconds=remaining,
        estimated_finish_utc=if(is.finite(remaining))format(Sys.time()+remaining,tz="UTC",usetz=TRUE)else NA_character_))
      last_progress<<-proc.time()[["elapsed"]]
    }
  }
  if(start<n_iter)for (iter in seq.int(start+1L,n_iter)) {
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
    checkpoint(iter)
  }
  checkpoint(n_iter,TRUE)
  colnames(samples_u) <- paste0("u", seq_len(n))
  colnames(samples_tau) <- paste0("tau", seq_len(m - 1L))
  list(
    samples_u = samples_u,
    samples_tau = samples_tau,
    data = list(
      c_ord = c_ord, u_obs = u_obs, calib_idx = calib_idx,
      miss_idx = miss_idx, m = m, tau_bounds = tau_bounds
    ),
    control = list(n_iter = n_iter, burn = burn, thin = thin, seed = seed),
    checkpoint = list(iteration=n_iter,u=u,tau=tau,rng=mixedgp_rng_state())
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
  methods <- if(is.null(measurement_fit)) "EIV-GP" else c("EIV-GP", "Response-free threshold model")
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
      (!is.null(measurement_fit) && is.null(measurement_fit$samples_u)) ||
      is.null(eiv_fit$data$miss_idx) ||
      (!is.null(measurement_fit) && is.null(measurement_fit$data$miss_idx))) {
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
  if (!is.null(measurement_fit) && !identical(miss_idx, measurement_miss_idx)) {
    stop("EIV-GP and response-free fits must use the same calibration rows.")
  }
  if (ncol(eiv_fit$mcmc$samples_u) != length(U_true) ||
      (!is.null(measurement_fit) && ncol(measurement_fit$samples_u) != length(U_true))) {
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
    if (!is.null(measurement_fit)) summarize_study1_training_u_draws(
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


# ---- simulations/study2/04_study2_ablations.R ----
############################################################
## 04_study2_ablations.R
##
## Appendix-only ablations for Study II.
##
## PI-GP and CC-GP are diagnostic decompositions of EIV-GP, not published
## competitors.  Their ordinal-probit measurement model is fitted without y.
## Full-U GP is an infeasible fitted benchmark that observes every training U.
############################################################


fit_ordinalprobit_measurement_fb <- function(
    C_ord,
    U_obs = NULL,
    calib_idx = integer(0),
    d = 2L,
    m_vec = NULL,
    ident = c("lower_triangular", "none"),
    n_iter = 2000L,
    burn = 600L,
    thin = 2L,
    n_chains = 4L,
    seed = 1L,
    s_A = 2.5,
    tau_bound = 6,
    parallel_chains = TRUE,
    n_cores = NULL,
    rhat_limit = 1.05,
    bulk_ess_limit = 50,
    tail_ess_limit = 50,
    verbose = FALSE) {
  ident <- match.arg(ident)
  C_ord <- as.matrix(C_ord)
  n <- nrow(C_ord)
  q <- ncol(C_ord)

  if (anyNA(C_ord)) stop("C_ord must be complete.")
  if (is.null(m_vec)) m_vec <- apply(C_ord, 2, max)
  m_vec <- as.integer(m_vec)
  if (length(m_vec) != q) stop("m_vec has the wrong length.")
  if (any(C_ord < 1L) || any(C_ord > matrix(rep(m_vec, each = n), n, q))) {
    stop("C_ord contains levels outside 1:m_j.")
  }
  if (ident == "lower_triangular" && q < d) {
    stop("lower_triangular identification requires q >= d.")
  }

  calib_idx <- sort(unique(as.integer(calib_idx)))
  if (any(calib_idx < 1L | calib_idx > n)) stop("Invalid calib_idx.")
  miss_idx <- setdiff(seq_len(n), calib_idx)
  U_obs_full <- matrix(NA_real_, n, d)

  if (length(calib_idx) > 0L) {
    if (is.null(U_obs)) stop("U_obs is required when calib_idx is nonempty.")
    U_obs <- as.matrix(U_obs)
    if (nrow(U_obs) == n && ncol(U_obs) == d) {
      U_obs_full[calib_idx, ] <- U_obs[calib_idx, , drop = FALSE]
    } else if (nrow(U_obs) == length(calib_idx) && ncol(U_obs) == d) {
      U_obs_full[calib_idx, ] <- U_obs
    } else {
      stop("U_obs must have n rows or length(calib_idx) rows and d columns.")
    }
  }
  anchor_status <- mixedgp_latent_anchor_status(
    U_obs_full,
    calib_idx,
    d = d
  )

  n_iter <- as.integer(n_iter)
  burn <- as.integer(burn)
  thin <- as.integer(thin)
  n_chains <- as.integer(n_chains)
  n_save <- floor((n_iter - burn) / thin)
  if (n_save < 1L || n_chains < 1L) {
    stop("The MCMC settings imply no saved measurement-model draws.")
  }
  rhat_limit <- as.numeric(rhat_limit)
  bulk_ess_limit <- as.numeric(bulk_ess_limit)
  tail_ess_limit <- as.numeric(tail_ess_limit)
  if (length(rhat_limit) != 1L || !is.finite(rhat_limit) || rhat_limit <= 1 ||
      length(bulk_ess_limit) != 1L || !is.finite(bulk_ess_limit) ||
      bulk_ess_limit <= 0 ||
      length(tail_ess_limit) != 1L || !is.finite(tail_ess_limit) ||
      tail_ess_limit <= 0) {
    stop("Measurement-model diagnostic limits are invalid.")
  }

  tau_names <- tau_names_from_mvec(m_vec)

  run_chain <- function(chain_id) {
    set.seed(seed + 10000L * chain_id)
    tau <- initialize_tau_ord(C_ord, m_vec, tau_bound)
    U <- init_U_from_ordinal(
      C_ord, d, m_vec,
      U_obs_full = U_obs_full,
      calib_idx = calib_idx
    )
    U <- U + matrix(rnorm(n * d, 0, 0.15), n, d)
    if (length(calib_idx) > 0L) U[calib_idx, ] <- U_obs_full[calib_idx, ]
    A <- initialize_A_ord(C_ord, U, tau, m_vec, ident)
    S <- sample_scores_ord(C_ord, U, A, tau)

    samples_U <- array(NA_real_, c(n_save, n, d))
    samples_A <- array(NA_real_, c(n_save, q, d))
    samples_tau <- matrix(NA_real_, n_save, length(tau_names))
    colnames(samples_tau) <- tau_names
    save_id <- 0L

    for (iter in seq_len(n_iter)) {
      S <- sample_scores_ord(C_ord, U, A, tau)
      tau <- update_tau_ord(tau, S, C_ord, m_vec, tau_bound)
      A <- update_A_ord(S, U, s_A = s_A, ident = ident)

      if (length(miss_idx) > 0L) {
        ref <- latent_reference_params(S, A)
        U[miss_idx, ] <- rmvnorm_rows_common(
          ref$mean[miss_idx, , drop = FALSE],
          ref$V
        )
      }
      if (length(calib_idx) > 0L) U[calib_idx, ] <- U_obs_full[calib_idx, ]

      if (iter > burn && (iter - burn) %% thin == 0L) {
        save_id <- save_id + 1L
        samples_U[save_id, , ] <- U
        samples_A[save_id, , ] <- A
        samples_tau[save_id, ] <- flatten_tau(tau)
      }
    }

    list(U = samples_U, A = samples_A, tau = samples_tau)
  }

  if (verbose) {
    message("Response-free ordinal-probit fit with ", n_chains, " chain(s).")
  }
  chain_seeds <- seed + 10000L * seq_len(n_chains)
  backend <- mixedgp_parallel_backend(
    if (isTRUE(parallel_chains)) n_cores else 1L
  )
  elapsed <- system.time({
    chains <- mixedgp_parallel_lapply(
      as.list(seq_len(n_chains)),
      run_chain,
      n_cores = if (isTRUE(parallel_chains)) {
        min(n_chains, backend$cores)
      } else {
        1L
      },
      seeds = chain_seeds,
      mc.preschedule = FALSE
    )
  })

  by_chain <- list(
    U = lapply(chains, `[[`, "U"),
    A = lapply(chains, `[[`, "A"),
    tau = lapply(chains, `[[`, "tau")
  )
  samples_U <- combine_chain_arrays(by_chain$U)
  samples_A <- combine_chain_arrays(by_chain$A)
  samples_tau <- do.call(rbind, by_chain$tau)

  key_chains <- list()
  key_blocks <- character(0)
  for (j in seq_len(q)) {
    for (k in seq_len(d)) {
      fixed_zero <- ident == "lower_triangular" && j <= d && k > j
      if (!fixed_zero) {
        parameter <- paste0("A[", j, ",", k, "]")
        key_chains[[parameter]] <- lapply(by_chain$A, function(z) z[, j, k])
        key_blocks[[parameter]] <- "loading"
      }
    }
  }
  for (k in seq_len(ncol(samples_tau))) {
    parameter <- colnames(samples_tau)[k]
    key_chains[[parameter]] <- lapply(by_chain$tau, function(z) z[, k])
    key_blocks[[parameter]] <- "cutpoint"
  }
  key_diagnostics <- data.frame(
    parameter = names(key_chains),
    block = unname(key_blocks[names(key_chains)]),
    rhat = vapply(key_chains, rank_rhat, numeric(1)),
    ess_bulk = vapply(key_chains, bulk_ess, numeric(1)),
    ess_tail = vapply(key_chains, tail_ess, numeric(1)),
    stringsAsFactors = FALSE
  )

  if (length(miss_idx) > 0L) {
    latent_diagnostics <- do.call(rbind, lapply(seq_len(d), function(k) {
      do.call(rbind, lapply(miss_idx, function(ii) {
        chains_ii <- lapply(by_chain$U, function(z) z[, ii, k])
        data.frame(
          parameter = paste0("U[", ii, ",", k, "]"),
          global_index = ii,
          coordinate = k,
          rhat = rank_rhat(chains_ii),
          ess_bulk = bulk_ess(chains_ii),
          ess_tail = tail_ess(chains_ii),
          stringsAsFactors = FALSE
        )
      }))
    }))
  } else {
    latent_diagnostics <- data.frame(
      parameter = character(0), global_index = integer(0),
      coordinate = integer(0), rhat = numeric(0),
      ess_bulk = numeric(0), ess_tail = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  diagnostic_backend <- if (requireNamespace("posterior", quietly = TRUE)) {
    "posterior_rank_normalized"
  } else {
    "fallback_without_tail_ess"
  }
  all_rhat <- c(key_diagnostics$rhat, latent_diagnostics$rhat)
  all_bulk <- c(key_diagnostics$ess_bulk, latent_diagnostics$ess_bulk)
  all_tail <- c(key_diagnostics$ess_tail, latent_diagnostics$ess_tail)
  rhat_pass <- length(all_rhat) > 0L && all(is.finite(all_rhat)) &&
    all(all_rhat <= rhat_limit)
  bulk_pass <- length(all_bulk) > 0L && all(is.finite(all_bulk)) &&
    all(all_bulk >= bulk_ess_limit)
  tail_pass <- length(all_tail) > 0L && all(is.finite(all_tail)) &&
    all(all_tail >= tail_ess_limit)
  convergence_pass <- identical(diagnostic_backend, "posterior_rank_normalized") &&
    n_chains >= 2L && rhat_pass && bulk_pass && tail_pass

  rhat_A <- key_diagnostics[key_diagnostics$block == "loading", , drop = FALSE]
  rhat_tau <- key_diagnostics[key_diagnostics$block == "cutpoint", , drop = FALSE]
  min_bulk <- safe_min(key_diagnostics$ess_bulk)
  min_tail <- safe_min(key_diagnostics$ess_tail)
  min_bulk_latent <- safe_min(latent_diagnostics$ess_bulk)
  min_tail_latent <- safe_min(latent_diagnostics$ess_tail)
  min_bulk_all <- safe_min(all_bulk)
  min_tail_all <- safe_min(all_tail)
  min_ess <- safe_min(pmin(
    all_bulk, all_tail, na.rm = TRUE
  ))

  list(
    data = list(
      C_ord = C_ord,
      U_obs = U_obs_full,
      calib_idx = calib_idx,
      miss_idx = miss_idx,
      latent_scale_anchored = isTRUE(anchor_status$anchored),
      latent_anchor_rank = anchor_status$affine_rank,
      latent_anchor_required_rank = anchor_status$required_rank,
      m_vec = m_vec,
      q = q,
      d = d,
      ident = ident
    ),
    mcmc = list(
      samples_U = samples_U,
      samples_A = samples_A,
      samples_tau = samples_tau,
      samples_by_chain = by_chain
    ),
    diagnostic_parameters = list(
      key = key_diagnostics,
      latent_rhat = latent_diagnostics
    ),
    diagnostics = data.frame(
      n_chains = n_chains,
      n_iter = n_iter,
      burn = burn,
      thin = thin,
      total_saved_draws = dim(samples_U)[1],
      max_rhat_A = safe_max(rhat_A$rhat),
      max_rhat_tau = safe_max(rhat_tau$rhat),
      max_rhat_missing_U = safe_max(latent_diagnostics$rhat),
      min_bulk_ess_key = min_bulk,
      min_tail_ess_key = min_tail,
      min_bulk_ess_missing_U = min_bulk_latent,
      min_tail_ess_missing_U = min_tail_latent,
      min_bulk_ess_all = min_bulk_all,
      min_tail_ess_all = min_tail_all,
      min_ess_key = min_ess,
      diagnostic_backend = diagnostic_backend,
      rhat_limit = rhat_limit,
      bulk_ess_limit = bulk_ess_limit,
      tail_ess_limit = tail_ess_limit,
      rhat_pass = rhat_pass,
      bulk_ess_pass = bulk_pass,
      tail_ess_pass = tail_pass,
      convergence_pass = convergence_pass,
      time_seconds = as.numeric(elapsed["elapsed"])
    )
  )
}

aggregate_rejection_pattern_telemetry <- function(telemetry_list) {
  telemetry_list <- Filter(
    function(x) is.data.frame(x) && nrow(x) > 0L,
    telemetry_list
  )
  if (length(telemetry_list) == 0L) return(NULL)
  raw <- do.call(rbind, telemetry_list)
  rownames(raw) <- NULL
  patterns <- unique(raw$pattern)
  out <- do.call(rbind, lapply(patterns, function(key) {
    z <- raw[raw$pattern == key, , drop = FALSE]
    total_candidates <- sum(z$total_candidates)
    proposal_hits <- sum(z$proposal_hits)
    retained_draws <- sum(z$retained_draws)
    data.frame(
      pattern = key,
      requested_draws = sum(z$requested_draws),
      retained_draws = retained_draws,
      proposal_hits = proposal_hits,
      shortfall = sum(z$shortfall),
      complete = all(z$complete),
      total_candidates = total_candidates,
      batches_used = sum(z$batches_used),
      empirical_pattern_probability = proposal_hits / total_candidates,
      retained_fraction_of_candidates = retained_draws / total_candidates,
      posterior_states = nrow(z),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

measurement_posterior_mean_u_test <- function(measurement_fit,
                                               C_test,
                                               max_draw = 200L,
                                               n_gibbs = 100L,
                                               seed = NULL,
                                               latent_sampler = c(
                                                 "minimax_tilting", "rejection", "gibbs"
                                               ),
                                               rejection_batch_size = NULL,
                                               rejection_max_batches = 1000L) {
  latent_sampler <- match.arg(latent_sampler)
  if (!is.null(seed)) set.seed(seed)
  C_test <- as.matrix(C_test)
  samples_A <- measurement_fit$mcmc$samples_A
  samples_tau <- measurement_fit$mcmc$samples_tau
  m_vec <- measurement_fit$data$m_vec
  n_available <- dim(samples_A)[1]
  ids <- if (n_available <= max_draw) {
    seq_len(n_available)
  } else {
    unique(round(seq(1, n_available, length.out = max_draw)))
  }

  U_sum <- matrix(0, nrow(C_test), measurement_fit$data$d)
  candidate_draws <- 0
  rejection_telemetry <- vector("list", length(ids))
  for (draw_position in seq_along(ids)) {
    s <- ids[[draw_position]]
    U_draw <- sample_u_given_c_ordprobit_dispatch(
      C_new = C_test,
      A = samples_A[s, , ],
      tau = unflatten_tau(samples_tau[s, ], m_vec),
      m_vec = m_vec,
      latent_sampler = latent_sampler,
      n_gibbs = n_gibbs,
      rejection_batch_size = rejection_batch_size,
      rejection_max_batches = rejection_max_batches
    )
    if (latent_sampler == "rejection") {
      candidate_draws <- candidate_draws + attr(U_draw, "candidate_draws")
      rejection_telemetry[[draw_position]] <-
        attr(U_draw, "pattern_telemetry")
    }
    U_sum <- U_sum + U_draw
  }
  out <- U_sum / length(ids)
  attr(out, "latent_sampler") <- switch(
    latent_sampler,
    minimax_tilting = "exact_minimax_tilting",
    rejection = "exact_rejection",
    gibbs = "finite_gibbs_diagnostic"
  )
  if (latent_sampler == "rejection") {
    requested_draws <- nrow(C_test) * length(ids)
    attr(out, "candidate_draws") <- candidate_draws
    attr(out, "overall_acceptance") <- requested_draws / candidate_draws
    pattern_telemetry <- aggregate_rejection_pattern_telemetry(
      rejection_telemetry
    )
    attr(out, "pattern_telemetry") <- pattern_telemetry
    attr(out, "rejection_telemetry") <- pattern_telemetry
  } else {
    attr(out, "gibbs_sweeps") <- as.integer(n_gibbs)
  }
  out
}

sample_measurement_u_given_c <- function(
    measurement_fit,
    C_new,
    draw_ids = NULL,
    max_draw = 200L,
    n_per_draw = 1L,
    scale = c("auto", "raw", "model"),
    latent_sampler = c("minimax_tilting", "rejection", "gibbs"),
    n_gibbs = 100L,
    rejection_batch_size = NULL,
    rejection_max_batches = 1000L,
    seed = NULL) {
  scale <- match.arg(scale)
  latent_sampler <- match.arg(latent_sampler)
  if (!is.null(seed)) set.seed(as.integer(seed))

  anchor <- study2_latent_anchor_from_fit(measurement_fit)
  if (scale == "raw" && !isTRUE(anchor$anchored)) {
    stop("scale='raw' requires a full-affine-rank calibration design.")
  }
  if (scale == "auto") scale <- if (isTRUE(anchor$anchored)) "raw" else "model"

  C_new <- as.matrix(C_new)
  m_vec <- measurement_fit$data$m_vec
  q <- measurement_fit$data$q
  d <- measurement_fit$data$d
  if (nrow(C_new) < 1L || ncol(C_new) != q || anyNA(C_new) ||
      any(C_new < 1L) ||
      any(C_new > matrix(rep(m_vec, each = nrow(C_new)), nrow(C_new), q))) {
    stop("C_new is incompatible with the fitted ordinal measurement model.")
  }

  samples_A <- measurement_fit$mcmc$samples_A
  samples_tau <- measurement_fit$mcmc$samples_tau
  n_available <- dim(samples_A)[1L]
  if (is.null(draw_ids)) {
    max_draw <- as.integer(max_draw)
    if (length(max_draw) != 1L || is.na(max_draw) || max_draw < 1L) {
      stop("max_draw must be one positive integer.")
    }
    draw_ids <- if (n_available <= max_draw) {
      seq_len(n_available)
    } else {
      unique(round(seq(1, n_available, length.out = max_draw)))
    }
  }
  draw_ids <- as.integer(draw_ids)
  if (length(draw_ids) < 1L || anyNA(draw_ids) ||
      any(draw_ids < 1L | draw_ids > n_available)) {
    stop("draw_ids contains invalid measurement-posterior indices.")
  }
  n_per_draw <- as.integer(n_per_draw)
  if (length(n_per_draw) != 1L || is.na(n_per_draw) || n_per_draw < 1L) {
    stop("n_per_draw must be one positive integer.")
  }

  out <- array(
    NA_real_,
    dim = c(length(draw_ids) * n_per_draw, nrow(C_new), d),
    dimnames = list(NULL, NULL, paste0("u", seq_len(d)))
  )
  acceptance <- numeric(dim(out)[1L])
  row_id <- 0L
  for (s in draw_ids) {
    A_s <- matrix(samples_A[s, , ], nrow = q, ncol = d)
    tau_s <- unflatten_tau(samples_tau[s, ], m_vec)
    for (rr in seq_len(n_per_draw)) {
      row_id <- row_id + 1L
      U_new <- sample_u_given_c_ordprobit_dispatch(
        C_new = C_new,
        A = A_s,
        tau = tau_s,
        m_vec = m_vec,
        latent_sampler = latent_sampler,
        n_gibbs = n_gibbs,
        rejection_batch_size = rejection_batch_size,
        rejection_max_batches = rejection_max_batches
      )
      out[row_id, , ] <- U_new
      acceptance[row_id] <- if (latent_sampler == "rejection") {
        attr(U_new, "overall_acceptance")
      } else {
        NA_real_
      }
    }
  }
  attr(out, "source") <- "prospective_response_free_measurement_model"
  attr(out, "latent_sampler") <- latent_sampler
  attr(out, "rejection_acceptance") <- acceptance
  attr(out, "scale") <- if (!isTRUE(anchor$anchored) && scale == "model") {
    "working"
  } else {
    scale
  }
  attr(out, "latent_scale_anchored") <- isTRUE(anchor$anchored)
  out
}

study2_measurement_training_imputation_metrics <- function(
    measurement_fit,
    U_true,
    rep_id,
    n_calib,
    scenario,
    method = "Ordinal model (no Y)") {
  miss_idx <- measurement_fit$data$miss_idx
  status <- study2_latent_imputation_status(
    fit = measurement_fit,
    method = method,
    target = "training_missing_U",
    rep_id = rep_id,
    n_calib = n_calib,
    scenario = scenario,
    n_units = length(miss_idx)
  )
  if (!isTRUE(status$task_eligible)) return(data.frame())
  summarize_study2_latent_draws(
    draws = measurement_fit$mcmc$samples_U[, miss_idx, , drop = FALSE],
    U_true = as.matrix(U_true)[miss_idx, , drop = FALSE],
    method = method,
    target = "training_missing_U",
    rep_id = rep_id,
    n_calib = n_calib,
    scenario = scenario
  )
}

study2_measurement_prospective_imputation_metrics <- function(
    measurement_fit,
    C_new,
    U_true,
    rep_id,
    n_calib,
    scenario,
    max_draw = 200L,
    latent_sampler = c("minimax_tilting", "rejection", "gibbs"),
    n_gibbs = 100L,
    rejection_batch_size = NULL,
    rejection_max_batches = 1000L,
    seed = NULL,
    method = "Ordinal model (no Y)") {
  latent_sampler <- match.arg(latent_sampler)
  U_true <- as.matrix(U_true)
  if (nrow(as.matrix(C_new)) != nrow(U_true)) {
    stop("C_new and U_true must contain the same prospective units.")
  }
  status <- study2_latent_imputation_status(
    fit = measurement_fit,
    method = method,
    target = "prospective_U_given_C",
    rep_id = rep_id,
    n_calib = n_calib,
    scenario = scenario,
    n_units = nrow(U_true)
  )
  if (!isTRUE(status$task_eligible)) return(data.frame())
  draws <- sample_measurement_u_given_c(
    measurement_fit = measurement_fit,
    C_new = C_new,
    max_draw = max_draw,
    n_per_draw = 1L,
    scale = "raw",
    latent_sampler = latent_sampler,
    n_gibbs = n_gibbs,
    rejection_batch_size = rejection_batch_size,
    rejection_max_batches = rejection_max_batches,
    seed = seed
  )
  summarize_study2_latent_draws(
    draws = draws,
    U_true = U_true,
    method = method,
    target = "prospective_U_given_C",
    rep_id = rep_id,
    n_calib = n_calib,
    scenario = scenario
  )
}

fit_study2_latent_gp <- function(X,
                                 y,
                                 U,
                                 train_idx = seq_len(nrow(X)),
                                 gp_n_starts = 5L,
                                 gp_seed = 1L,
                                 gp_maxit = 500L) {
  X <- as.matrix(X)
  U <- as.matrix(U)
  y <- as.numeric(y)
  train_idx <- as.integer(train_idx)
  if (length(train_idx) < 3L) stop("At least three training observations are required.")

  X_center <- colMeans(X)
  X_scale <- apply(X, 2, sd)
  X_scale[!is.finite(X_scale) | X_scale <= 0] <- 1
  X_std <- sweep(sweep(X, 2, X_center, "-"), 2, X_scale, "/")
  y_center <- mean(y[train_idx])
  y_scale <- sd(y[train_idx])
  if (!is.finite(y_scale) || y_scale <= 0) y_scale <- 1

  fit <- gp_mle_fit(
    X = cbind(X_std[train_idx, , drop = FALSE], U[train_idx, , drop = FALSE]),
    y = (y[train_idx] - y_center) / y_scale,
    n_starts = gp_n_starts,
    seed = gp_seed,
    maxit = gp_maxit
  )

  list(
    fit = fit,
    X_center = X_center,
    X_scale = X_scale,
    y_center = y_center,
    y_scale = y_scale,
    train_idx = train_idx,
    optimizer_attempts = fit$optimizer_attempts,
    selected_start = fit$selected_start,
    optimizer_control = fit$optimizer_control
  )
}

prepare_gp_mle_prediction <- function(fit) {
  X <- fit$X
  y <- fit$y
  n <- nrow(X)
  d <- ncol(X)
  sigma2 <- exp(fit$par[1])
  rho <- exp(fit$par[2])
  theta <- exp(fit$par[-c(1, 2)])

  Rexp <- matrix(0, n, n)
  for (j in seq_len(d)) {
    Rexp <- Rexp + theta[j] * pairwise_sqdist(X[, j, drop = FALSE])
  }
  K <- rho^2 * sigma2 * exp(-Rexp)
  chol_train <- safe_chol(K + sigma2 * diag(n))

  list(
    X = X,
    alpha = solve_chol(chol_train, y),
    chol_train = chol_train,
    sigma2 = sigma2,
    rho = rho,
    theta = theta
  )
}

gp_mle_predict_prepared <- function(prepared, Xstar, noisy = FALSE) {
  Xstar <- as.matrix(Xstar)
  Rstar_exp <- matrix(0, nrow(Xstar), nrow(prepared$X))
  for (j in seq_along(prepared$theta)) {
    Rstar_exp <- Rstar_exp + prepared$theta[j] * pairwise_sqdist(
      Xstar[, j, drop = FALSE],
      prepared$X[, j, drop = FALSE]
    )
  }
  Kstar <- prepared$rho^2 * prepared$sigma2 * exp(-Rstar_exp)
  mean <- as.numeric(Kstar %*% prepared$alpha)
  v <- forwardsolve(t(prepared$chol_train), t(Kstar))
  var <- prepared$rho^2 * prepared$sigma2 - colSums(v^2)
  var_tolerance <- 100 *
    (nrow(prepared$X) + nrow(Xstar)) * .Machine$double.eps *
    max(1, prepared$rho^2 * prepared$sigma2)
  if (any(var < -var_tolerance)) {
    stop("Appendix GP conditional variance is materially negative.")
  }
  var <- pmax(var, 0)
  if (isTRUE(noisy)) var <- var + prepared$sigma2
  list(mean = mean, var = var)
}

sample_study2_latent_gp <- function(gp_fit,
                                    X_test,
                                    U_test,
                                    n_draw = 500L,
                                    noisy = TRUE,
                                    prepared = NULL) {
  X_test <- as.matrix(X_test)
  U_test <- as.matrix(U_test)
  X_std <- sweep(
    sweep(X_test, 2, gp_fit$X_center, "-"),
    2, gp_fit$X_scale, "/"
  )
  if (is.null(prepared)) prepared <- prepare_gp_mle_prediction(gp_fit$fit)
  pred <- gp_mle_predict_prepared(
    prepared,
    Xstar = cbind(X_std, U_test),
    noisy = noisy
  )
  draws_std <- matrix(
    rnorm(
      n_draw * nrow(X_test),
      mean = rep(pred$mean, each = n_draw),
      sd = rep(sqrt(pred$var), each = n_draw)
    ),
    nrow = n_draw,
    ncol = nrow(X_test)
  )
  out <- gp_fit$y_center + gp_fit$y_scale * draws_std
  attr(out, "conditional_means") <- matrix(
    gp_fit$y_center + gp_fit$y_scale * pred$mean,
    nrow = 1L
  )
  attr(out, "conditional_vars") <- matrix(
    gp_fit$y_scale^2 * pred$var,
    nrow = 1L
  )
  attr(out, "mixture_components") <- "fitted latent-input GP Gaussian predictor"
  out
}

fit_study2_pi_gp <- function(X,
                             y,
                             measurement_fit,
                             gp_n_starts = 5L,
                             gp_seed = 1L,
                             gp_maxit = 500L) {
  U_hat <- apply(measurement_fit$mcmc$samples_U, c(2, 3), mean)
  list(
    gp = fit_study2_latent_gp(
      X, y, U_hat,
      gp_n_starts = gp_n_starts,
      gp_seed = gp_seed,
      gp_maxit = gp_maxit
    ),
    measurement = measurement_fit,
    U_hat_train = U_hat
  )
}

sample_study2_pi_gp <- function(fit,
                                X_test,
                                C_test,
                                n_draw,
                                n_measurement_draw = 200L,
                                n_gibbs = 100L,
                                seed = NULL,
                                latent_sampler = c(
                                  "minimax_tilting", "rejection", "gibbs"
                                ),
                                rejection_batch_size = NULL,
                                rejection_max_batches = 1000L) {
  latent_sampler <- match.arg(latent_sampler)
  if (!is.null(seed)) set.seed(seed)
  U_hat_test <- measurement_posterior_mean_u_test(
    fit$measurement,
    C_test,
    max_draw = n_measurement_draw,
    n_gibbs = n_gibbs,
    latent_sampler = latent_sampler,
    rejection_batch_size = rejection_batch_size,
    rejection_max_batches = rejection_max_batches
  )
  out <- sample_study2_latent_gp(
    fit$gp, X_test, U_hat_test, n_draw, noisy = TRUE
  )
  attr(out, "latent_sampler") <- attr(U_hat_test, "latent_sampler")
  attr(out, "candidate_draws") <- attr(U_hat_test, "candidate_draws")
  attr(out, "overall_acceptance") <- attr(U_hat_test, "overall_acceptance")
  attr(out, "pattern_telemetry") <- attr(U_hat_test, "pattern_telemetry")
  attr(out, "rejection_telemetry") <- attr(U_hat_test, "rejection_telemetry")
  attr(out, "gibbs_sweeps") <- attr(U_hat_test, "gibbs_sweeps")
  out
}

fit_study2_cc_gp <- function(X,
                             y,
                             U,
                             calib_idx,
                             measurement_fit,
                             gp_n_starts = 5L,
                             gp_seed = 1L,
                             gp_maxit = 500L) {
  if (length(calib_idx) < 3L) {
    stop("CC-GP requires at least three calibrated observations.")
  }
  list(
    gp = fit_study2_latent_gp(
      X, y, U,
      train_idx = calib_idx,
      gp_n_starts = gp_n_starts,
      gp_seed = gp_seed,
      gp_maxit = gp_maxit
    ),
    measurement = measurement_fit
  )
}

sample_study2_cc_gp <- function(fit,
                                X_test,
                                C_test,
                                n_draw,
                                n_gibbs = 100L,
                                seed = NULL,
                                latent_sampler = c(
                                  "minimax_tilting", "rejection", "gibbs"
                                ),
                                rejection_batch_size = NULL,
                                rejection_max_batches = 1000L) {
  latent_sampler <- match.arg(latent_sampler)
  if (!is.null(seed)) set.seed(seed)
  samples_A <- fit$measurement$mcmc$samples_A
  samples_tau <- fit$measurement$mcmc$samples_tau
  m_vec <- fit$measurement$data$m_vec
  ids <- sample(seq_len(dim(samples_A)[1]), n_draw, replace = TRUE)
  out <- matrix(NA_real_, n_draw, nrow(X_test))
  conditional_means <- matrix(NA_real_, n_draw, nrow(X_test))
  conditional_vars <- matrix(NA_real_, n_draw, nrow(X_test))
  prepared <- prepare_gp_mle_prediction(fit$gp$fit)
  candidate_draws <- 0
  rejection_telemetry <- vector("list", n_draw)

  for (ii in seq_len(n_draw)) {
    s <- ids[ii]
    U_test <- sample_u_given_c_ordprobit_dispatch(
      C_new = C_test,
      A = samples_A[s, , ],
      tau = unflatten_tau(samples_tau[s, ], m_vec),
      m_vec = m_vec,
      latent_sampler = latent_sampler,
      n_gibbs = n_gibbs,
      rejection_batch_size = rejection_batch_size,
      rejection_max_batches = rejection_max_batches
    )
    if (latent_sampler == "rejection") {
      candidate_draws <- candidate_draws + attr(U_test, "candidate_draws")
      rejection_telemetry[[ii]] <- attr(U_test, "pattern_telemetry")
    }
    draw_ii <- sample_study2_latent_gp(
      fit$gp, X_test, U_test, n_draw = 1L, noisy = TRUE,
      prepared = prepared
    )
    out[ii, ] <- draw_ii
    conditional_means[ii, ] <- attr(draw_ii, "conditional_means")[1L, ]
    conditional_vars[ii, ] <- attr(draw_ii, "conditional_vars")[1L, ]
  }
  attr(out, "latent_sampler") <- switch(
    latent_sampler,
    minimax_tilting = "exact_minimax_tilting",
    rejection = "exact_rejection",
    gibbs = "finite_gibbs_diagnostic"
  )
  if (latent_sampler == "rejection") {
    requested_draws <- n_draw * nrow(X_test)
    attr(out, "candidate_draws") <- candidate_draws
    attr(out, "overall_acceptance") <- requested_draws / candidate_draws
    pattern_telemetry <- aggregate_rejection_pattern_telemetry(
      rejection_telemetry
    )
    attr(out, "pattern_telemetry") <- pattern_telemetry
    attr(out, "rejection_telemetry") <- pattern_telemetry
  } else {
    attr(out, "gibbs_sweeps") <- as.integer(n_gibbs)
  }
  attr(out, "conditional_means") <- conditional_means
  attr(out, "conditional_vars") <- conditional_vars
  attr(out, "mixture_components") <-
    "measurement and complete-case GP Monte Carlo Gaussian components"
  out
}

study2_latent_gp_surface_metrics <- function(gp_fit,
                                             scenario,
                                             method,
                                             rep_id,
                                             n_calib = NA_integer_,
                                             grid_n = 31L,
                                             u_lim = c(-2.5, 2.5)) {
  grid <- expand.grid(
    u1 = seq(u_lim[1], u_lim[2], length.out = grid_n),
    u2 = seq(u_lim[1], u_lim[2], length.out = grid_n)
  )
  X_star <- matrix(0, nrow(grid), 2L)
  U_star <- as.matrix(grid)
  X_std <- sweep(
    sweep(X_star, 2, gp_fit$X_center, "-"),
    2, gp_fit$X_scale, "/"
  )
  pred <- gp_mle_predict(
    gp_fit$fit,
    Xstar = cbind(X_std, U_star),
    noisy = FALSE
  )
  post_mean <- gp_fit$y_center + gp_fit$y_scale * pred$mean
  post_sd <- gp_fit$y_scale * sqrt(pred$var)
  truth <- f0_2d(X_star, U_star, scenario = scenario)
  lower <- post_mean - qnorm(0.975) * post_sd
  upper <- post_mean + qnorm(0.975) * post_sd

  data.frame(
    rep = rep_id,
    scenario = scenario,
    n_calib = n_calib,
    method = method,
    grid_n = nrow(grid),
    ISE = mean((post_mean - truth)^2),
    Bias = mean(post_mean - truth),
    Coverage95 = mean(truth >= lower & truth <= upper),
    Width95 = mean(upper - lower)
  )
}


# ---- simulation_helpers.R ----
MIXEDGP_SIMULATION_SCHEMA <- "2.1.1"
MIXEDGP_PUBLISHED_METHODS <- c("UC-GP", "LVGP", "EzGP")


mixedgp_simulation_code_dir <- function(code_dir = NULL) {
  if (is.null(code_dir)) code_dir <- system.file(package = "eivGP")
  normalizePath(code_dir, mustWork = TRUE)
}
mixedgp_simulation_modules <- function() character()
mixedgp_simulation_engine <- function(code_dir = NULL) {
  mixedgp_simulation_code_dir(code_dir)
  ns <- asNamespace("eivGP")
  env <- new.env(parent = ns)
  for (n in ls(ns, all.names=TRUE)) {
    value <- get(n, ns)
    if(is.function(value) && identical(environment(value),ns)) environment(value) <- env
    assign(n,value,env)
  }
  env

}


# ---- simulations/helpers/design.R ----
# Simulation design: package implementation.

mixedgp_validate_scalar_integer <- function(x, name, minimum = 0L) {
  if (length(x) != 1L || is.na(x) || !is.numeric(x) || x != as.integer(x) ||
      x < minimum) {
    stop(name, " must be one integer at least ", minimum, ".")
  }
  as.integer(x)
}

mixedgp_validate_boolean <- function(x, name) {
  if (length(x) != 1L || !is.logical(x) || is.na(x)) {
    stop(name, " must be TRUE or FALSE.")
  }
  isTRUE(x)
}

mixedgp_validate_calibration_grid <- function(x, n, name = "calibration_grid") {
  if (!is.numeric(x) || length(x) < 1L || anyNA(x) ||
      any(x != as.integer(x)) || any(x < 0L | x > n) || anyDuplicated(x)) {
    stop(name, " must contain unique integer sizes between zero and n.")
  }
  sort(as.integer(x))
}

mixedgp_available_workers <- function(cap = 8L) {
  cap <- mixedgp_validate_scalar_integer(cap, "cap", minimum = 1L)
  detected <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (!is.finite(detected) || detected < 1L) detected <- 1L
  as.integer(min(cap, detected))
}

mixedgp_simulation_mode <- function(mode) {
  match.arg(mode, c("dry_run", "data", "smoke", "development", "publication"))
}

mixedgp_development_profile <- function(config) {
  integer_env <- function(name, default, minimum = 1L) {
    value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
    mixedgp_validate_scalar_integer(value, name, minimum)
  }
  selected <- vapply(config$cells, `[[`, character(1L), "id")
  requested <- Sys.getenv("MIXEDGP_DEV_CELLS", unset = "")
  ids <- vapply(config$cells, `[[`, character(1L), "id")
  if (nzchar(requested)) selected <- trimws(strsplit(requested, ",", fixed = TRUE)[[1L]])
  if (identical(selected, "all")) selected <- ids
  if (anyDuplicated(selected) || any(!selected %in% ids)) {
    stop("MIXEDGP_DEV_CELLS must name unique cells from: ", paste(ids, collapse = ", "))
  }
  config$cells <- config$cells[match(selected, ids)]
  reps <- integer_env("MIXEDGP_DEV_REPS", 3L)
  burn <- integer_env("MIXEDGP_DEV_BURN", 500L, 0L)
  draws <- integer_env("MIXEDGP_DEV_DRAWS", 1250L)
  config$cells <- lapply(config$cells, function(cell) {
    cell$n_rep <- reps
    cell
  })
  config$run_id <- paste0(config$study, "-development")
  config$stages <- c("data", "fit", "aggregate")
  config$strict_competitors <- FALSE
  config$fail_closed <- FALSE
  config$mcmc$n_iter <- burn + draws
  config$mcmc$burn <- burn
  config$mcmc$thin <- 1L
  config$mcmc$n_chains <- 4L
  config$mcmc$require_gate <- FALSE
  config$mcmc$automatic_continuation <- FALSE
  config$measurement_mcmc$n_iter <- burn + draws
  config$measurement_mcmc$burn <- burn
  if (config$study == "study2") config$measurement_mcmc$thin <- 1L
  config$evaluation$n_pred_draw <- 100L
  config$evaluation$n_m_eval <- 30L
  config$evaluation$n_m_draw <- 50L
  config$evaluation$n_m_latent <- 64L
  if (config$study == "study2") {
    config$evaluation$n_m_truth <- 500L
    config$evaluation$n_oracle_pool <- 10000L
  }
  config$development_note <- paste(
    "DEVELOPMENT RESULTS: reduced replication and Monte Carlo budgets.",
    "For draft figures and workflow inspection; review diagnostics and rerun",
    "the frozen publication design before reporting scientific conclusions."
  )
  config
}

mixedgp_apply_core_budget <- function(config, core_budget = NULL) {
  if (is.null(core_budget)) {
    value <- Sys.getenv("MIXEDGP_CORE_BUDGET", unset = "")
    if (!nzchar(value)) return(config)
    core_budget <- suppressWarnings(as.numeric(value))
  }
  allocation <- eivgp_run_settings(
    core_budget, max(vapply(config$cells, `[[`, integer(1), "n_rep")),
    chains = config$mcmc$n_chains,
    sampling_iterations = config$mcmc$n_iter - config$mcmc$burn,
    warmup = config$mcmc$burn)
  config$parallel <- allocation[c("level", "workers", "chain_workers", "core_budget", "active_chain_limit")]
  config
}

mixedgp_number_slug <- function(x) {
  gsub("[^0-9A-Za-z]+", "p", format(x, scientific = FALSE, trim = TRUE))
}

mixedgp_config_fingerprint <- function(config) {
  payload <- config
  payload$created_at <- NULL
  payload$stages <- NULL
  payload$source_identity <- eivgp_module_identity()
  package_names <- c(
    ggplot2 = "ggplot2", dplyr = "dplyr", tidyr = "tidyr",
    knitr = "knitr", posterior = "posterior"
  )
  if (identical(config$study, "study2")) {
    package_names <- c(package_names, TruncatedNormal = "TruncatedNormal")
  }
  payload$runtime_package_identity <- vapply(
    package_names,
    function(package) {
      if (requireNamespace(package, quietly = TRUE)) {
        as.character(utils::packageVersion(package))
      } else {
        "not-installed"
      }
    },
    character(1L)
  )
  path <- tempfile("mixedgp-config-", fileext = ".rds")
  on.exit(if (file.exists(path)) unlink(path), add = TRUE)
  saveRDS(payload, path, version = 3L, compress = FALSE)
  unname(tools::md5sum(path))
}

mixedgp_estimand_method_matrix <- function(study = c("study1", "study2")) {
  study <- match.arg(study)
  rows <- list(
    c("surface", "m(x,c)", "EIV-GP", "proposed", "yes", "yes", "all calibration sizes", ""),
    c("surface", "m(x,c)", "UC-GP", "published competitor", "yes", "point only", "all calibration sizes", "does not define physical u"),
    c("surface", "m(x,c)", "LVGP", "published competitor", "yes", "point only", "all calibration sizes", "does not define physical u"),
    c("surface", "m(x,c)", "EzGP", "published competitor", "yes", "point only", "all calibration sizes", "does not define physical u"),
    c("surface", "f(x,u)", "EIV-GP", "proposed", "yes", "yes", "positive, anchoring calibration", "physical u scale is not identified at k=0"),
    c("surface", "f(x,u)", "PI-GP", "scientific ablation", "yes", "optimizer approximation", "positive, anchoring calibration", "plug-in ignores latent-input uncertainty"),
    c("surface", "f(x,u)", "CC-GP", "scientific ablation", "yes", "optimizer approximation", "at least three complete cases", "uses only calibrated cases"),
    c("surface", "f(x,u)", "Full-U GP", "infeasible benchmark", "yes", "optimizer approximation", "all training u observed", "benchmark, not a deployable competitor"),
    c("prediction", "Y*|x*,c*", "EIV-GP", "proposed", "yes", "yes", "all calibration sizes", ""),
    c("prediction", "Y*|x*,c*", "UC-GP", "published competitor", "yes", "predictive", "all calibration sizes", ""),
    c("prediction", "Y*|x*,c*", "LVGP", "published competitor", "yes", "predictive", "all calibration sizes", ""),
    c("prediction", "Y*|x*,c*", "EzGP", "published competitor", "yes", "predictive", "all calibration sizes", ""),
    c("prediction", "Y*|x*,c*", "PI-GP", "scientific ablation", "yes", "predictive", "cells with ablations", "imputes a latent input and then plugs it into a GP"),
    c("prediction", "Y*|x*,c*", "CC-GP", "scientific ablation", "yes", "predictive", "positive calibration in cells with ablations", "fits the response GP only to complete cases"),
    c("prediction", "Y*|x*,c*", "DGM oracle", "truth benchmark", "yes", "truth", "simulation only", "never ranked as a fitted method"),
    c("latent state", "training U|X,C,Y", "EIV-GP", "proposed", "yes", "yes", "positive, anchoring calibration", "physical coordinates are not identified at k=0"),
    c("latent state", "training U|C", "RF-OM", "response-free ablation", "yes", "yes", "positive, anchoring calibration", "does not use y")
  )
  if (study == "study2") {
    rows <- c(rows, list(
      c("latent state", "prospective U*|X*,C*,Y-data", "EIV-GP", "proposed", "yes", "yes", "positive, anchoring calibration", "y* is not conditioned on"),
      c("latent state", "prospective U*|C*", "RF-OM", "response-free ablation", "yes", "yes", "positive, anchoring calibration", "does not use response data")
    ))
  }
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(out) <- c(
    "task", "estimand", "method", "role", "eligible",
    "uncertainty_supported", "condition", "reason_or_limit"
  )
  out$study <- study
  out[, c("study", setdiff(names(out), "study")), drop = FALSE]
}

mixedgp_study1_cells <- function(mode) {
  smoke <- identical(mode, "smoke")
  cells <- list()
  for (threshold_design in "balanced") {
    for (eta in c(0, 1)) {
      cells[[length(cells) + 1L]] <- list(
        id = paste0("eta", eta, "_", threshold_design),
        role = "primary_heterogeneity",
        scenario = "heterogeneity_continuum",
        heterogeneity_eta = eta, threshold_design = threshold_design,
        min_class_count = if (threshold_design == "imbalanced") 3L else 0L,
        n = 100L, n_test = if (smoke) 80L else 100L,
        n_rep = if (smoke) 1L else 50L, m = 6L,
        calibration_grid = if (smoke) c(0L, 5L) else c(0L, 20L, 50L),
        evaluate_f = identical(eta, 1), evaluate_u = TRUE,
        run_ablations = TRUE)
    }
  }
  cells
}

mixedgp_study2_cells <- function(mode) {
  smoke <- identical(mode, "smoke")
  make_cell <- function(id, role, scenario, q, grid,
                        run_ablations, evaluate_f = FALSE, evaluate_u = TRUE) {
    list(id = id, role = role, scenario = scenario, q = as.integer(q),
         d = 2L, m = 4L, n = 100L, n_test = if (smoke) 60L else 100L,
         n_rep = if (smoke) 1L else 50L, calibration_grid = as.integer(grid),
         run_ablations = run_ablations, evaluate_f = evaluate_f,
         evaluate_u = evaluate_u)
  }
  list(
    make_cell("primary_q2", "correct_measurement_specification", "primary", 2L,
              if (smoke) 6L else 50L, TRUE),
    make_cell("primary_q4_calibration", "primary_calibration_curve", "primary", 4L,
              if (smoke) c(0L, 6L) else c(0L, 20L, 50L, 80L), TRUE,
              evaluate_f = TRUE),
    make_cell("logistic_q4", "measurement_misspecification", "logistic_misspec", 4L,
              if (smoke) 6L else 50L, FALSE, evaluate_u = FALSE)
  )
}

study1_simulation_config <- function(
    mode = "dry_run",
    code_dir = NULL,
    workers = mixedgp_available_workers(),
    output_root = NULL,
    data_root = NULL, core_budget = NULL) {
  mode <- mixedgp_simulation_mode(mode)
  code_dir <- mixedgp_simulation_code_dir(code_dir)
  revision_dir <- dirname(code_dir)
  if (mode == "development") {
    if (is.null(output_root)) output_root <- file.path(revision_dir, "reproduction", "development", "results", "study1")
    if (is.null(data_root)) data_root <- file.path(revision_dir, "reproduction", "development", "data", "study1")
  }
  if (is.null(output_root)) {
    output_root <- file.path(revision_dir, "simulation-runs", "study1")
  }
  if (is.null(data_root)) {
    data_root <- file.path(
      revision_dir, "data-synthetic", "publication-v2", "study1"
    )
  }
  smoke <- identical(mode, "smoke")
  config <- list(
    schema_version = MIXEDGP_SIMULATION_SCHEMA,
    study = "study1",
    mode = mode,
    run_id = paste0("study1-", mode, "-v2"),
    code_dir = normalizePath(code_dir, winslash = "/"),
    output_root = normalizePath(output_root, winslash = "/", mustWork = FALSE),
    data_root = normalizePath(data_root, winslash = "/", mustWork = FALSE),
    stages = switch(
      mode,
      dry_run = character(0), data = "data",
      smoke = c("data", "fit", "aggregate"),
      publication = c("data", "fit", "aggregate")
    ),
    cells = mixedgp_study1_cells(mode),
    published_methods = character(),
    strict_competitors = identical(mode, "publication"),
    fail_closed = identical(mode, "publication"),
    use_cache = TRUE,
    parallel = list(level = "replications", workers = as.integer(workers)),
    mcmc = if (smoke) {
      list(n_iter = 120L, burn = 40L, thin = 1L, n_chains = 1L,
           rhat_limit = 1.05, ess_limit = 20L, require_gate = FALSE)
    } else {
      list(n_iter = if (mode == "publication") 20000L else 5000L,
           burn = if (mode == "publication") 5000L else 1000L,
           thin = 1L, n_chains = 4L,
           rhat_limit = 1.01, ess_limit = 400L,
           require_gate = identical(mode, "publication"))
    },
    measurement_mcmc = if (smoke) {
      list(n_iter = 120L, burn = 40L)
    } else {
      list(n_iter = 4000L, burn = 1000L)
    },
    evaluation = if (smoke) {
      list(n_pred_draw = 40L, n_m_eval = 12L, n_m_draw = 20L,
           n_m_latent = 24L)
    } else {
      list(n_pred_draw = 500L, n_m_eval = 200L, n_m_draw = 250L,
           n_m_latent = 512L)
    },
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  class(config) <- c("mixedgp_simulation_config", "list")
  if (mode == "development") config <- mixedgp_development_profile(config)
  config <- mixedgp_apply_core_budget(config, core_budget)
  validate_simulation_config(config)
}

study2_simulation_config <- function(
    mode = "dry_run",
    code_dir = NULL,
    workers = mixedgp_available_workers(),
    output_root = NULL,
    data_root = NULL, core_budget = NULL) {
  mode <- mixedgp_simulation_mode(mode)
  code_dir <- mixedgp_simulation_code_dir(code_dir)
  revision_dir <- dirname(code_dir)
  if (mode == "development") {
    if (is.null(output_root)) output_root <- file.path(revision_dir, "reproduction", "development", "results", "study2")
    if (is.null(data_root)) data_root <- file.path(revision_dir, "reproduction", "development", "data", "study2")
  }
  if (is.null(output_root)) {
    output_root <- file.path(revision_dir, "simulation-runs", "study2")
  }
  if (is.null(data_root)) {
    data_root <- file.path(
      revision_dir, "data-synthetic", "publication-v2", "study2"
    )
  }
  smoke <- identical(mode, "smoke")
  config <- list(
    schema_version = MIXEDGP_SIMULATION_SCHEMA,
    study = "study2",
    mode = mode,
    run_id = paste0("study2-", mode, "-v2"),
    code_dir = normalizePath(code_dir, winslash = "/"),
    output_root = normalizePath(output_root, winslash = "/", mustWork = FALSE),
    data_root = normalizePath(data_root, winslash = "/", mustWork = FALSE),
    stages = switch(
      mode,
      dry_run = character(0), data = "data",
      smoke = c("data", "fit", "aggregate"),
      publication = c("data", "fit", "aggregate")
    ),
    cells = mixedgp_study2_cells(mode),
    published_methods = character(),
    strict_competitors = identical(mode, "publication"),
    fail_closed = identical(mode, "publication"),
    use_cache = TRUE,
    parallel = list(level = "replications", workers = as.integer(workers)),
    predictive_latent_sampler = "minimax_tilting",
    mcmc = if (smoke) {
      list(n_iter = 120L, burn = 40L, thin = 1L, n_chains = 1L,
           rhat_limit = 1.05, raw_ess_limit = 20L,
           target_bulk_ess_limit = 10L, target_tail_ess_limit = 10L,
           require_gate = FALSE)
    } else {
      list(n_iter = if (mode == "publication") 20000L else 4000L,
           burn = if (mode == "publication") 5000L else 1000L,
           thin = 1L, n_chains = 4L,
           rhat_limit = 1.01, raw_ess_limit = 400L,
           target_bulk_ess_limit = 400L, target_tail_ess_limit = 400L,
           ess_min = if(mode == "publication")100 else NULL,
           require_gate = identical(mode, "publication"))
    },
    measurement_mcmc = if (smoke) {
      list(n_iter = 120L, burn = 40L, thin = 1L, n_chains = 2L)
    } else {
      list(n_iter = 4000L, burn = 1000L, thin = 1L, n_chains = 4L)
    },
    evaluation = if (smoke) {
      list(n_pred_draw = 30L, n_m_eval = 8L, n_m_draw = 12L,
           n_m_latent = 16L, n_m_truth = 100L, n_oracle_pool = 5000L)
    } else {
      list(n_pred_draw = 500L, n_m_eval = 150L, n_m_draw = 200L,
           n_m_latent = 256L, n_m_truth = 5000L,
           n_oracle_pool = 250000L)
    },
    ablation_gp = list(
      n_starts = if (smoke) 3L else 8L,
      maxit = if (smoke) 200L else 500L
    ),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  class(config) <- c("mixedgp_simulation_config", "list")
  if (mode == "development") config <- mixedgp_development_profile(config)
  config <- mixedgp_apply_core_budget(config, core_budget)
  validate_simulation_config(config)
}

validate_simulation_config <- function(config) {
  if (!is.list(config) || !identical(
    config$schema_version, MIXEDGP_SIMULATION_SCHEMA
  )) {
    stop("Unsupported or malformed simulation configuration.")
  }
  config$study <- match.arg(config$study, c("study1", "study2"))
  config$mode <- mixedgp_simulation_mode(config$mode)
  allowed_stages <- c("data", "fit", "aggregate")
  if (anyDuplicated(config$stages) || any(!config$stages %in% allowed_stages)) {
    stop("stages must be a unique subset of data, fit, and aggregate.")
  }
  if ("aggregate" %in% config$stages && !"fit" %in% config$stages) {
    stop("The aggregate stage requires fit in the same invocation.")
  }
  if ("fit" %in% config$stages && !"data" %in% config$stages) {
    ## Fitting frozen data without regenerating it is valid, but the artifacts
    ## are verified before dispatch. This branch is intentionally allowed.
    invisible(NULL)
  }
  if (!is.list(config$cells) || length(config$cells) < 1L) {
    stop("config$cells must be a nonempty list.")
  }
  ids <- vapply(config$cells, `[[`, character(1L), "id")
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Simulation cell ids must be nonempty and unique.")
  }
  required <- c("id", "role", "scenario", "n", "n_test", "n_rep", "m",
                "calibration_grid", "run_ablations", "evaluate_f", "evaluate_u")
  for (cell in config$cells) {
    study_required <- if (config$study == "study1") {
      c("heterogeneity_eta", "threshold_design", "min_class_count")
    } else {
      c("q", "d")
    }
    required_cell <- c(required, study_required)
    missing <- setdiff(required_cell, names(cell))
    if (length(missing) > 0L) {
      stop("Cell ", cell$id, " is missing: ", paste(missing, collapse = ", "))
    }
    n <- mixedgp_validate_scalar_integer(cell$n, paste0(cell$id, "$n"), 1L)
    mixedgp_validate_scalar_integer(cell$n_test, paste0(cell$id, "$n_test"), 1L)
    mixedgp_validate_scalar_integer(cell$n_rep, paste0(cell$id, "$n_rep"), 1L)
    mixedgp_validate_scalar_integer(cell$m, paste0(cell$id, "$m"), 2L)
    calibration_grid <- mixedgp_validate_calibration_grid(
      cell$calibration_grid, n, paste0(cell$id, "$calibration_grid")
    )
    if (!identical(as.integer(cell$calibration_grid), calibration_grid)) {
      stop(cell$id, "$calibration_grid must be sorted increasingly.")
    }
    for (flag in c("run_ablations", "evaluate_f", "evaluate_u")) {
      mixedgp_validate_boolean(cell[[flag]], paste0(cell$id, "$", flag))
    }
    if (config$study == "study1") {
      eta <- as.numeric(cell$heterogeneity_eta)
      if (length(eta) != 1L || !is.finite(eta) || eta < 0) {
        stop(cell$id, "$heterogeneity_eta must be finite and nonnegative.")
      }
      if (!cell$threshold_design %in% c("balanced", "imbalanced")) {
        stop(cell$id, "$threshold_design must be balanced or imbalanced.")
      }
      mixedgp_validate_scalar_integer(
        cell$min_class_count, paste0(cell$id, "$min_class_count"), 0L
      )
      if (cell$threshold_design == "imbalanced" && cell$m != 6L) {
        stop("The imbalanced Study I design requires m=6 in ", cell$id, ".")
      }
    } else {
      mixedgp_validate_scalar_integer(cell$q, paste0(cell$id, "$q"), 2L)
      if (cell$q > 6L) stop("Study II publication cells require q <= 6.")
      mixedgp_validate_scalar_integer(cell$d, paste0(cell$id, "$d"), 1L)
      if (cell$d != 2L) stop("The Study II DGM currently requires d=2.")
    }
  }
  if (length(config$published_methods) && !identical(config$published_methods, MIXEDGP_PUBLISHED_METHODS)) {
    stop(
      "The publication competitor set is frozen as: ",
      paste(MIXEDGP_PUBLISHED_METHODS, collapse = ", "), "."
    )
  }
  config$parallel$workers <- mixedgp_validate_scalar_integer(
    config$parallel$workers, "parallel$workers", 1L
  )
  if (!config$parallel$level %in% c("replications", "hybrid")) {
    stop("parallel$level must be replications or hybrid.")
  }
  if (config$parallel$level == "hybrid") {
    budget <- mixedgp_validate_scalar_integer(config$parallel$core_budget, "core_budget", 1L)
    chain_workers <- mixedgp_validate_scalar_integer(config$parallel$chain_workers, "chain_workers", 1L)
    if (config$parallel$workers * chain_workers > budget) {
      stop("Dataset workers times chain workers exceeds core_budget.")
    }
  }
  for (flag in c("strict_competitors", "fail_closed", "use_cache")) {
    mixedgp_validate_boolean(config[[flag]], flag)
  }
  for (field in c("n_iter", "burn", "thin", "n_chains")) {
    minimum <- if (field == "burn") 0L else 1L
    mixedgp_validate_scalar_integer(
      config$mcmc[[field]], paste0("mcmc$", field), minimum
    )
  }
  if (config$mcmc$burn >= config$mcmc$n_iter) {
    stop("mcmc$burn must be smaller than mcmc$n_iter.")
  }
  mixedgp_validate_boolean(config$mcmc$require_gate, "mcmc$require_gate")
  if (!is.numeric(config$mcmc$rhat_limit) ||
      length(config$mcmc$rhat_limit) != 1L ||
      !is.finite(config$mcmc$rhat_limit) || config$mcmc$rhat_limit <= 1) {
    stop("mcmc$rhat_limit must be one finite number greater than one.")
  }
  measurement_fields <- intersect(
    c("n_iter", "burn", "thin", "n_chains"),
    names(config$measurement_mcmc)
  )
  for (field in measurement_fields) {
    minimum <- if (field == "burn") 0L else 1L
    mixedgp_validate_scalar_integer(
      config$measurement_mcmc[[field]],
      paste0("measurement_mcmc$", field), minimum
    )
  }
  if (config$measurement_mcmc$burn >= config$measurement_mcmc$n_iter) {
    stop("measurement_mcmc$burn must be smaller than measurement_mcmc$n_iter.")
  }
  evaluation_fields <- if (config$study == "study1") {
    c("n_pred_draw", "n_m_eval", "n_m_draw", "n_m_latent")
  } else {
    c(
      "n_pred_draw", "n_m_eval", "n_m_draw", "n_m_latent",
      "n_m_truth", "n_oracle_pool"
    )
  }
  missing_evaluation <- setdiff(evaluation_fields, names(config$evaluation))
  if (length(missing_evaluation) > 0L) {
    stop("evaluation is missing: ", paste(missing_evaluation, collapse = ", "))
  }
  for (field in evaluation_fields) {
    mixedgp_validate_scalar_integer(
      config$evaluation[[field]], paste0("evaluation$", field), 1L
    )
  }
  if (config$study == "study1") {
    mixedgp_validate_scalar_integer(config$mcmc$ess_limit, "mcmc$ess_limit", 1L)
  } else {
    for (field in c(
      "raw_ess_limit", "target_bulk_ess_limit", "target_tail_ess_limit"
    )) {
      mixedgp_validate_scalar_integer(
        config$mcmc[[field]], paste0("mcmc$", field), 1L
      )
    }
    if (!identical(config$predictive_latent_sampler, "minimax_tilting")) {
      stop("Publication Study II uses the exact minimax_tilting latent sampler.")
    }
    mixedgp_validate_scalar_integer(
      config$ablation_gp$n_starts, "ablation_gp$n_starts", 1L
    )
    mixedgp_validate_scalar_integer(
      config$ablation_gp$maxit, "ablation_gp$maxit", 1L
    )
  }
  ## Diagnostics assess reliability; they do not terminate either run mode.
  if (config$mode %in% c("publication", "development")) {
    config$strict_competitors <- FALSE
    config$fail_closed <- FALSE
    config$mcmc$require_gate <- FALSE
  }
  config
}



# ---- simulations/helpers/data.R ----
# Simulation data: package implementation.

mixedgp_simulation_diagnostic_advice <- function() {
  paste("Inspect diagnostic tables and chain traces before interpretation.",
        "If mixing is stable but ESS is low or MCSE is high, consider explicit",
        "continuation with continue_eivgp() for a compatible public fit, then",
        "recompute diagnostics. High R-hat or separated chains require checking",
        "identification, initialization and sampler behavior; more iterations",
        "alone may not help. Missing outputs/packages require repair or rerunning",
        "the affected method, not longer MCMC. Continuation is never automatic.")
}

mixedgp_bind_rows_base <- function(rows) {
  rows <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, rows)
  if (length(rows) == 0L) return(data.frame())
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  aligned <- lapply(rows, function(x) {
    missing <- setdiff(columns, names(x))
    for (name in missing) x[[name]] <- NA
    x[, columns, drop = FALSE]
  })
  out <- do.call(rbind, aligned)
  rownames(out) <- NULL
  out
}

mixedgp_task_plan <- function(config) {
  config <- validate_simulation_config(config)
  rows <- list()
  add <- function(cell, task_type, method, n_calib = NA_integer_) {
    n <- cell$n_rep * if (all(is.na(n_calib))) 1L else length(n_calib)
    rows[[length(rows) + 1L]] <<- data.frame(
      study = config$study,
      cell_id = cell$id,
      design_role = cell$role,
      task_type = task_type,
      method = method,
      n_calib = if (all(is.na(n_calib))) NA_integer_ else rep(n_calib, each = cell$n_rep),
      rep = if (all(is.na(n_calib))) seq_len(cell$n_rep) else rep(seq_len(cell$n_rep), times = length(n_calib)),
      planned_units = rep(1L, n),
      stringsAsFactors = FALSE
    )
  }
  for (cell in config$cells) {
    add(cell, "frozen_data", "DGM")
    add(cell, "prediction_truth", "DGM oracle")
    for (method in config$published_methods) {
      add(cell, "published_fit", method)
    }
    for (n_calib in cell$calibration_grid) {
      add(cell, "eiv_fit", "EIV-GP", n_calib)
    }
    if (isTRUE(cell$run_ablations)) {
      for (method in c("RF-OM", "PI-GP", "CC-GP")) {
        for (n_calib in cell$calibration_grid) {
          add(cell, "ablation_fit", method, n_calib)
        }
      }
      if (isTRUE(cell$evaluate_f)) add(cell, "benchmark_fit", "Full-U GP")
    }
  }
  mixedgp_bind_rows_base(rows)
}

mixedgp_cell_summary <- function(config) {
  mixedgp_bind_rows_base(lapply(config$cells, function(cell) {
    data.frame(
      cell_id = cell$id,
      role = cell$role,
      scenario = cell$scenario,
      eta = if (is.null(cell$heterogeneity_eta)) NA_real_ else cell$heterogeneity_eta,
      q = if (is.null(cell$q)) NA_integer_ else cell$q,
      n = cell$n,
      n_test = cell$n_test,
      n_rep = cell$n_rep,
      calibration_grid = paste(cell$calibration_grid, collapse = ";"),
      evaluate_f = cell$evaluate_f,
      evaluate_u = cell$evaluate_u,
      run_ablations = cell$run_ablations,
      stringsAsFactors = FALSE
    )
  }))
}

mixedgp_runtime_preflight <- function(config) {
  packages <- c("ggplot2", "dplyr", "tidyr", "knitr", "posterior")
  if (config$study == "study2") packages <- c(packages, "TruncatedNormal")
  available <- vapply(packages, requireNamespace, logical(1L), quietly = TRUE)
  version <- vapply(seq_along(packages), function(ii) {
    if (available[[ii]]) {
      as.character(utils::packageVersion(packages[[ii]]))
    } else {
      NA_character_
    }
  }, character(1L))
  data.frame(
    package = packages,
    role = c(
      "figures", "summaries and diagnostics", "summary reshaping",
      "LaTeX tables", "rank-normalized MCMC diagnostics",
      if (config$study == "study2") "exact minimax-tilting sampler"
    ),
    available = available,
    version = version,
    stringsAsFactors = FALSE
  )
}

mixedgp_run_directory <- function(config) {
  fingerprint <- mixedgp_config_fingerprint(config)
  file.path(config$output_root, paste0(config$run_id, "-", substr(fingerprint, 1L, 12L)))
}

mixedgp_data_cell_directory <- function(config, cell) {
  file.path(config$data_root, cell$id)
}

mixedgp_code_hashes <- function(config) {
  ids <- eivgp_module_identity()
  data.frame(file = names(ids), exists=TRUE, md5 = unname(ids), stringsAsFactors = FALSE)
}

mixedgp_atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(object, temporary, version = 3L)
  if (file.exists(path) && !file.remove(path)) {
    stop("Could not replace existing file: ", path)
  }
  if (!file.rename(temporary, path)) stop("Could not install file: ", path)
  invisible(path)
}

mixedgp_atomic_write_csv <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  utils::write.csv(object, temporary, row.names = FALSE, na = "")
  if (file.exists(path) && !file.remove(path)) {
    stop("Could not replace existing file: ", path)
  }
  if (!file.rename(temporary, path)) stop("Could not install file: ", path)
  invisible(path)
}

mixedgp_write_run_metadata <- function(config,
                                       run_dir,
                                       preflight,
                                       runtime_preflight,
                                       task_plan) {
  metadata_dir <- file.path(run_dir, "config")
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
  mixedgp_atomic_save_rds(config, file.path(metadata_dir, "resolved_config.rds"))
  config_text <- capture.output(dput(config))
  writeLines(config_text, file.path(metadata_dir, "resolved_config.R"))
  mixedgp_atomic_write_csv(
    mixedgp_cell_summary(config), file.path(metadata_dir, "design_cells.csv")
  )
  mixedgp_atomic_write_csv(
    mixedgp_estimand_method_matrix(config$study),
    file.path(metadata_dir, "estimand_method_matrix.csv")
  )
  mixedgp_atomic_write_csv(task_plan, file.path(metadata_dir, "task_plan.csv"))
  mixedgp_atomic_write_csv(preflight, file.path(metadata_dir, "competitor_preflight.csv"))
  mixedgp_atomic_write_csv(
    runtime_preflight, file.path(metadata_dir, "runtime_preflight.csv")
  )
  mixedgp_atomic_write_csv(
    mixedgp_code_hashes(config), file.path(metadata_dir, "code_hashes.csv")
  )
  capture.output(
    sessionInfo(), file = file.path(metadata_dir, "sessionInfo.txt")
  )
  invisible(metadata_dir)
}

mixedgp_existing_data_status <- function(config) {
  mixedgp_bind_rows_base(lapply(config$cells, function(cell) {
    directory <- mixedgp_data_cell_directory(config, cell)
    manifest_path <- file.path(directory, "manifest.rds")
    manifest <- if (file.exists(manifest_path)) {
      tryCatch(readRDS(manifest_path), error = function(e) NULL)
    } else {
      NULL
    }
    compatible <- if (is.data.frame(manifest)) {
      keep <- manifest$scenario == cell$scenario &
        manifest$n == cell$n & manifest$n_test == cell$n_test &
        manifest$m == cell$m & manifest$rep %in% seq_len(cell$n_rep)
      if (config$study == "study2") keep <- keep & manifest$q == cell$q
      manifest[keep, , drop = FALSE]
    } else {
      data.frame()
    }
    complete <- nrow(compatible) == cell$n_rep &&
      setequal(compatible$rep, seq_len(cell$n_rep))
    data.frame(
      cell_id = cell$id,
      directory = directory,
      manifest_exists = file.exists(manifest_path),
      frozen_artifacts = nrow(compatible),
      expected_artifacts = cell$n_rep,
      complete_count = complete,
      stringsAsFactors = FALSE
    )
  }))
}

mixedgp_print_dry_run <- function(config, preflight, runtime_preflight, task_plan) {
  cat(
    "Publication simulation dry run\n",
    "  study: ", config$study, "\n",
    "  schema: ", config$schema_version, "\n",
    "  run directory: ", mixedgp_run_directory(config), "\n",
    "  parallelism: ", config$parallel$workers,
    " dataset worker(s); chain workers per dataset: ",
    if (is.null(config$parallel$chain_workers)) 1L else config$parallel$chain_workers, "\n",
    sep = ""
  )
  print(mixedgp_cell_summary(config), row.names = FALSE)
  cat("\nPlanned fit/data units by type and method:\n")
  task_counts <- aggregate(
    planned_units ~ task_type + method, task_plan, sum
  )
  print(task_counts, row.names = FALSE)
  cat("\nPublished-package preflight:\n")
  print(preflight, row.names = FALSE)
  cat("\nRuntime-package preflight:\n")
  print(runtime_preflight, row.names = FALSE)
  cat("\nFrozen-data status (informational during dry run):\n")
  print(mixedgp_existing_data_status(config), row.names = FALSE)
  invisible(config)
}

mixedgp_generate_cell_data <- function(config, cell, engine) {
  directory <- mixedgp_data_cell_directory(config, cell)
  generator <- if (config$study == "study1") {
    get("generate_study1_synthetic_datasets", envir = engine)
  } else {
    get("generate_study2_synthetic_datasets", envir = engine)
  }
  common <- list(
    n_rep = cell$n_rep,
    directory = directory,
    n_cores = config$parallel$workers,
    overwrite = FALSE,
    n = cell$n,
    n_test = cell$n_test,
    m = cell$m,
    calib_grid = cell$calibration_grid
  )
  args <- if (config$study == "study1") {
    c(common, list(
      scenario = cell$scenario,
      threshold_design = cell$threshold_design,
      min_class_count = cell$min_class_count,
      heterogeneity_eta = cell$heterogeneity_eta,
      data_seed_base = 100000L,
      calibration_seed_base = 200000L
    ))
  } else {
    c(common, list(
      scenarios = cell$scenario,
      q = cell$q,
      data_seed_base = 1000000L
    ))
  }
  do.call(generator, args)
}

mixedgp_verify_cell_data <- function(config, cell, engine = NULL) {
  directory <- mixedgp_data_cell_directory(config, cell)
  manifest_path <- file.path(directory, "manifest.rds")
  if (!file.exists(manifest_path)) {
    stop("Frozen-data manifest is missing for cell ", cell$id, ": ", manifest_path)
  }
  manifest <- readRDS(manifest_path)
  if (!is.data.frame(manifest)) stop("Invalid manifest: ", manifest_path)
  matching <- manifest$scenario == cell$scenario &
    manifest$n == cell$n & manifest$n_test == cell$n_test &
    manifest$m == cell$m & manifest$rep %in% seq_len(cell$n_rep)
  if (config$study == "study2") {
    matching <- matching & manifest$q == cell$q
  }
  selected <- manifest[matching, , drop = FALSE]
  if (nrow(selected) != cell$n_rep ||
      !setequal(selected$rep, seq_len(cell$n_rep))) {
    stop(
      "Frozen-data manifest for ", cell$id, " contains ", nrow(selected),
      " compatible artifacts; expected ", cell$n_rep, "."
    )
  }
  paths <- file.path(directory, selected$file)
  if (any(!file.exists(paths))) stop("A frozen artifact is missing for ", cell$id, ".")
  observed_md5 <- unname(tools::md5sum(paths))
  if (!identical(tolower(observed_md5), tolower(as.character(selected$md5)))) {
    stop("Frozen-data checksum mismatch for cell ", cell$id, ".")
  }
  if (!is.null(engine)) {
    loader <- get(
      "load_mixedgp_synthetic_dataset_strict",
      envir = engine,
      inherits = FALSE
    )
    for (ii in seq_len(nrow(selected))) {
      expected <- list(
        study = config$study,
        scenario = cell$scenario,
        rep_id = as.integer(selected$rep[[ii]]),
        n = cell$n,
        n_test = cell$n_test,
        m = cell$m,
        calib_grid = cell$calibration_grid
      )
      if (config$study == "study1") {
        expected$threshold_design <- cell$threshold_design
        expected$heterogeneity_eta <- cell$heterogeneity_eta
      } else {
        expected$q <- cell$q
      }
      loader(
        paths[[ii]], expected = expected,
        manifest_path = manifest_path
      )
    }
  }
  selected
}

mixedgp_read_cell_replication <- function(config, cell, rep_id) {
  directory <- mixedgp_data_cell_directory(config, cell)
  manifest <- readRDS(file.path(directory, "manifest.rds"))
  matching <- manifest$rep == rep_id & manifest$scenario == cell$scenario &
    manifest$n == cell$n & manifest$n_test == cell$n_test &
    manifest$m == cell$m
  if (config$study == "study2") matching <- matching & manifest$q == cell$q
  row <- manifest[matching, , drop = FALSE]
  if (nrow(row) != 1L) {
    stop(
      "Expected one full-design frozen artifact for ", cell$id,
      " replication ", rep_id, "; found ", nrow(row), "."
    )
  }
  path <- file.path(directory, row$file[[1L]])
  observed_md5 <- unname(tools::md5sum(path))
  if (!identical(tolower(observed_md5), tolower(as.character(row$md5[[1L]])))) {
    stop("Frozen-data checksum mismatch for ", cell$id, " replication ", rep_id, ".")
  }
  artifact <- readRDS(path)
  if (!identical(as.integer(artifact$design$calib_grid),
                 as.integer(cell$calibration_grid))) {
    stop("Frozen calibration grid does not match cell ", cell$id, ".")
  }
  if (config$study == "study1" &&
      (!identical(as.character(artifact$design$threshold_design),
                  as.character(cell$threshold_design)) ||
       !isTRUE(all.equal(
         as.numeric(artifact$design$heterogeneity_eta),
         as.numeric(cell$heterogeneity_eta), tolerance = 0
       )))) {
    stop("Frozen Study I mechanism does not match cell ", cell$id, ".")
  }
  artifact
}

mixedgp_validate_common_random_numbers <- function(config, engine = NULL) {
  checks <- list()
  same_numeric <- function(x, y) {
    isTRUE(all.equal(
      x, y, tolerance = sqrt(.Machine$double.eps),
      check.attributes = FALSE
    ))
  }
  add_check <- function(reference, comparison, component, pass, n_rep) {
    checks[[length(checks) + 1L]] <<- data.frame(
      reference_cell = reference,
      comparison_cell = comparison,
      component = component,
      replications_checked = n_rep,
      pass = isTRUE(pass),
      stringsAsFactors = FALSE
    )
    if (!isTRUE(pass)) {
      stop(
        "Common-random-number validation failed for ", comparison,
        " component ", component, "."
      )
    }
  }
  by_id <- setNames(config$cells, vapply(config$cells, `[[`, character(1L), "id"))
  if (config$study == "study1") {
    for (design in c("balanced", "imbalanced")) {
    primary_ids <- intersect(
      paste0("eta", c(0, 1), "_", design), names(by_id)
    )
    if (length(primary_ids) >= 2L) {
      reference <- by_id[[primary_ids[[1L]]]]
      for (id in primary_ids[-1L]) {
        comparison <- by_id[[id]]
        n_check <- min(reference$n_rep, comparison$n_rep)
        component_pass <- c(inputs = TRUE, noise = TRUE, calibration = TRUE,
                            observable_mean = TRUE)
        for (rr in seq_len(n_check)) {
          a <- mixedgp_read_cell_replication(config, reference, rr)
          b <- mixedgp_read_cell_replication(config, comparison, rr)
          component_pass[["inputs"]] <- component_pass[["inputs"]] &&
            identical(a$data$train$x, b$data$train$x) &&
            identical(a$data$train$u, b$data$train$u) &&
            identical(a$data$train$c, b$data$train$c) &&
            identical(a$data$test$x, b$data$test$x) &&
            identical(a$data$test$u, b$data$test$u) &&
            identical(a$data$test$c, b$data$test$c)
          component_pass[["noise"]] <- component_pass[["noise"]] &&
            same_numeric(a$data$train$y - a$data$train$f,
                         b$data$train$y - b$data$train$f) &&
            same_numeric(a$data$test$y - a$data$test$f,
                         b$data$test$y - b$data$test$f)
          component_pass[["calibration"]] <-
            component_pass[["calibration"]] &&
            identical(a$calibration_sets, b$calibration_sets)
          if (is.null(engine) || !exists("m0_1d", envir = engine, inherits = FALSE)) {
            stop("Study I paired-design validation requires m0_1d in engine.")
          }
          ma <- get("m0_1d", envir = engine, inherits = FALSE)(
            a$data$test$x, a$data$test$c, a$data$tau_true,
            scenario = a$scenario
          )
          mb <- get("m0_1d", envir = engine, inherits = FALSE)(
            b$data$test$x, b$data$test$c, b$data$tau_true,
            scenario = b$scenario
          )
          component_pass[["observable_mean"]] <-
          component_pass[["observable_mean"]] && same_numeric(ma, mb)
        }
        for (component in names(component_pass)) {
          add_check(reference$id, id, component, component_pass[[component]], n_check)
        }
      }
    }
    }
  } else {
    reference_id <- if ("primary_q4_calibration" %in% names(by_id)) {
      "primary_q4_calibration"
    } else {
      names(by_id)[[1L]]
    }
    reference <- by_id[[reference_id]]
    for (id in setdiff(names(by_id), reference_id)) {
      comparison <- by_id[[id]]
      n_check <- min(reference$n_rep, comparison$n_rep)
      same_n_test <- identical(reference$n_test, comparison$n_test)
      component_pass <- c(
        train_X_U = TRUE,
        train_response_noise = TRUE,
        test_X_U_same_size = if (same_n_test) TRUE else NA,
        test_response_noise_same_size = if (same_n_test) TRUE else NA
      )
      if (comparison$scenario == "primary") {
        component_pass <- c(
          component_pass,
          nested_train_C = TRUE,
          nested_test_C_same_size = if (same_n_test) TRUE else NA
        )
      } else if (comparison$scenario == "latent_additive_control" &&
                 identical(reference$q, comparison$q)) {
        component_pass <- c(
          component_pass,
          identical_measurement_C = TRUE,
          identical_common_calibration_sets = TRUE
        )
      }
      for (rr in seq_len(n_check)) {
        a <- mixedgp_read_cell_replication(config, reference, rr)
        b <- mixedgp_read_cell_replication(config, comparison, rr)
        component_pass[["train_X_U"]] <- component_pass[["train_X_U"]] &&
          identical(a$data$train$X, b$data$train$X) &&
          identical(a$data$train$U, b$data$train$U)
        component_pass[["train_response_noise"]] <-
          component_pass[["train_response_noise"]] &&
          same_numeric(a$data$train$y - a$data$train$f,
                       b$data$train$y - b$data$train$f)
        if (same_n_test) {
          component_pass[["test_X_U_same_size"]] <-
            component_pass[["test_X_U_same_size"]] &&
            identical(a$data$test$X, b$data$test$X) &&
            identical(a$data$test$U, b$data$test$U)
          component_pass[["test_response_noise_same_size"]] <-
            component_pass[["test_response_noise_same_size"]] &&
            same_numeric(a$data$test$y - a$data$test$f,
                         b$data$test$y - b$data$test$f)
        }
        if (comparison$scenario == "primary") {
          q_shared <- min(ncol(a$data$train$C), ncol(b$data$train$C))
          component_pass[["nested_train_C"]] <-
            component_pass[["nested_train_C"]] &&
            identical(a$data$train$C[, seq_len(q_shared), drop = FALSE],
                      b$data$train$C[, seq_len(q_shared), drop = FALSE])
          if (same_n_test) {
            component_pass[["nested_test_C_same_size"]] <-
              component_pass[["nested_test_C_same_size"]] &&
              identical(a$data$test$C[, seq_len(q_shared), drop = FALSE],
                        b$data$test$C[, seq_len(q_shared), drop = FALSE])
          }
        } else if (comparison$scenario == "latent_additive_control" &&
                   identical(reference$q, comparison$q)) {
          component_pass[["identical_measurement_C"]] <-
            component_pass[["identical_measurement_C"]] &&
            identical(a$data$train$C, b$data$train$C) &&
            (!same_n_test || identical(a$data$test$C, b$data$test$C))
          common_calibration <- intersect(
            names(a$calibration_sets), names(b$calibration_sets)
          )
          component_pass[["identical_common_calibration_sets"]] <-
            component_pass[["identical_common_calibration_sets"]] &&
            length(common_calibration) > 0L &&
            all(vapply(common_calibration, function(k) {
              identical(a$calibration_sets[[k]], b$calibration_sets[[k]])
            }, logical(1L)))
        }
      }
      for (component in names(component_pass)[!is.na(component_pass)]) {
        add_check(reference$id, id, component, component_pass[[component]], n_check)
      }
    }
  }
  mixedgp_bind_rows_base(checks)
}



# ---- simulations/helpers/execution.R ----
# Simulation execution: package implementation.

mixedgp_cell_controls_study1 <- function(config, cell, cell_output) {
  mechanism_calib <- if (50L %in% cell$calibration_grid) {
    50L
  } else {
    max(cell$calibration_grid)
  }
  list(
    STUDY1_CONFIG = if (config$mode == "smoke") "quick" else "thorough",
    STUDY1_QUICK = config$mode == "smoke",
    STUDY1_OUT_PREFIX = cell_output,
    STUDY1_DATA_DIR = mixedgp_data_cell_directory(config, cell),
    STUDY1_SCENARIO = cell$scenario,
    STUDY1_HETEROGENEITY_ETA = cell$heterogeneity_eta,
    STUDY1_THRESHOLD_DESIGN = cell$threshold_design,
    STUDY1_CALIB_GRID = cell$calibration_grid,
    STUDY1_MC_N_TRAIN = cell$n,
    STUDY1_MC_N_TEST = cell$n_test,
    STUDY1_MC_N_REP = cell$n_rep,
    STUDY1_MC_M = cell$m,
    STUDY1_MC_N_ITER = config$mcmc$n_iter,
    STUDY1_MC_BURN = config$mcmc$burn,
    STUDY1_MC_N_CHAINS = config$mcmc$n_chains,
    STUDY1_MC_N_PRED_DRAW = config$evaluation$n_pred_draw,
    STUDY1_MC_N_M_EVAL = config$evaluation$n_m_eval,
    STUDY1_MC_N_M_DRAW = config$evaluation$n_m_draw,
    STUDY1_MC_N_M_LATENT = config$evaluation$n_m_latent,
    STUDY1_MEAS_N_ITER = config$measurement_mcmc$n_iter,
    STUDY1_MEAS_BURN = config$measurement_mcmc$burn,
    STUDY1_USE_CACHE = config$use_cache,
    STUDY1_REUSE_LOCKED_EIV = FALSE,
    STUDY1_STRICT_COMPETITORS = config$strict_competitors,
    STUDY1_PUBLISHED_COMPETITORS = config$published_methods,
    STUDY1_RUN_ABLATIONS = cell$run_ablations,
    STUDY1_EVALUATE_F = cell$evaluate_f,
    STUDY1_EVALUATE_U = cell$evaluate_u,
    STUDY1_PARALLEL_LEVEL = config$parallel$level,
    STUDY1_DATASET_WORKERS = config$parallel$workers,
    STUDY1_CHAIN_WORKERS = if (is.null(config$parallel$chain_workers)) 1L else config$parallel$chain_workers,
    STUDY1_REQUIRE_MCMC_GATE = config$mcmc$require_gate,
    STUDY1_MAX_RHAT = config$mcmc$rhat_limit,
    STUDY1_MIN_ESS = config$mcmc$ess_limit,
    STUDY1_MECHANISM_CALIB = mechanism_calib
  )
}

mixedgp_cell_controls_study2 <- function(config, cell, cell_output) {
  primary_grid <- if (cell$scenario == "primary") {
    cell$calibration_grid
  } else {
    sort(unique(c(0L, cell$calibration_grid)))
  }
  contrast_calib <- max(cell$calibration_grid)
  list(
    STUDY2_CONFIG = if (config$mode == "smoke") "quick" else "thorough",
    STUDY2_OUT_PREFIX = cell_output,
    STUDY2_DATA_DIR = mixedgp_data_cell_directory(config, cell),
    STUDY2_SCENARIOS = cell$scenario,
    STUDY2_PRIMARY_CALIB_GRID = primary_grid,
    STUDY2_CONTRAST_CALIB = contrast_calib,
    STUDY2_Q = cell$q,
    STUDY2_M = cell$m,
    STUDY2_MC_N_TRAIN = cell$n,
    STUDY2_MC_N_TEST = cell$n_test,
    STUDY2_MC_N_REP = cell$n_rep,
    STUDY2_MC_N_ITER = config$mcmc$n_iter,
    STUDY2_MC_BURN = config$mcmc$burn,
    STUDY2_MC_THIN = config$mcmc$thin,
    STUDY2_MC_N_CHAINS = config$mcmc$n_chains,
    STUDY2_MC_N_PRED_DRAW = config$evaluation$n_pred_draw,
    STUDY2_MC_N_M_EVAL = config$evaluation$n_m_eval,
    STUDY2_MC_N_M_DRAW = config$evaluation$n_m_draw,
    STUDY2_MC_N_M_LATENT = config$evaluation$n_m_latent,
    STUDY2_MC_N_M_TRUTH = config$evaluation$n_m_truth,
    STUDY2_MC_N_ORACLE_POOL = config$evaluation$n_oracle_pool,
    STUDY2_MEAS_N_ITER = config$measurement_mcmc$n_iter,
    STUDY2_MEAS_BURN = config$measurement_mcmc$burn,
    STUDY2_MEAS_THIN = config$measurement_mcmc$thin,
    STUDY2_MEAS_N_CHAINS = config$measurement_mcmc$n_chains,
    STUDY2_ABLATION_GP_N_STARTS = config$ablation_gp$n_starts,
    STUDY2_ABLATION_GP_MAXIT = config$ablation_gp$maxit,
    STUDY2_USE_CACHE = config$use_cache,
    STUDY2_MC_RESUME = TRUE,
    STUDY2_SAVE_REP_FITS = FALSE,
    STUDY2_STRICT_COMPETITORS = config$strict_competitors,
    STUDY2_PUBLISHED_COMPETITORS = config$published_methods,
    STUDY2_RUN_ABLATIONS = cell$run_ablations,
    STUDY2_ABLATION_SCENARIOS = if (cell$run_ablations) cell$scenario else character(0),
    STUDY2_EVALUATE_F = cell$evaluate_f,
    STUDY2_EVALUATE_U = cell$evaluate_u,
    STUDY2_PARALLEL_LEVEL = config$parallel$level,
    STUDY2_DATASET_WORKERS = config$parallel$workers,
    STUDY2_CHAIN_WORKERS = if (is.null(config$parallel$chain_workers)) 1L else config$parallel$chain_workers,
    STUDY2_PREDICTIVE_LATENT_SAMPLER = config$predictive_latent_sampler,
    STUDY2_ENFORCE_MCMC_GATE = config$mcmc$require_gate,
    STUDY2_MCMC_RHAT_LIMIT = config$mcmc$rhat_limit,
    STUDY2_ESS_MIN = config$mcmc$ess_min,
    STUDY2_MCMC_RAW_ESS_LIMIT = config$mcmc$raw_ess_limit,
    STUDY2_MCMC_TARGET_BULK_ESS_LIMIT = config$mcmc$target_bulk_ess_limit,
    STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT = config$mcmc$target_tail_ess_limit,
    ## Smoke mode is explicitly non-reportable: permissive measurement-chain
    ## thresholds let PI-GP/CC-GP execute so every code path is tested. The
    ## publication constructor retains the prespecified strict diagnostics.
    STUDY2_MEAS_RHAT_LIMIT = if (config$mode == "smoke") {
      3
    } else {
      config$mcmc$rhat_limit
    },
    STUDY2_MEAS_BULK_ESS_LIMIT = if (config$mode == "smoke") {
      1L
    } else {
      config$mcmc$target_bulk_ess_limit
    },
    STUDY2_MEAS_TAIL_ESS_LIMIT = if (config$mode == "smoke") {
      1L
    } else {
      config$mcmc$target_tail_ess_limit
    },
    STUDY2_SAVE_PDF = config$mode %in% c("publication", "development"),
    STUDY2_SAVE_PNG = FALSE
  )
}

mixedgp_get_cell_output <- function(envir, name) {
  get0(name, envir = envir, inherits = FALSE, ifnotfound = data.frame())
}

## Catch errors inside workers before mclapply converts them to try-error
## strings. Keep successful replications and expose the original failure.
mixedgp_save_replication <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(basename(path), ".tmp-"), dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) stop("Could not commit replication: ", path)
  invisible(path)
}

mixedgp_run_replications <- function(X, FUN, n_cores, seeds, status_path,
                                    study, parallel_map, mc.preschedule = FALSE) {
  ## Each worker owns one small status file. No shared CSV writes from forks.
  task_dir <- paste0(status_path, ".tasks")
  dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
  task_paths <- file.path(task_dir, sprintf("task-%05d.csv", seq_along(X)))
  record <- function(i, state, message = "", elapsed = NA_real_) {
    mixedgp_atomic_write_csv(data.frame(
      task = paste(X[[i]], collapse = ","), status = state,
      message = message, elapsed_seconds = elapsed,
      updated_at = format(Sys.time(), tz = "UTC", usetz = TRUE)), task_paths[i])
  }
  for (i in seq_along(X)) if (!file.exists(task_paths[i])) record(i, "pending")
  timing <- mixedgp_batch_timing(paste0(status_path,".timing"),
    vapply(X,function(x)paste(x,collapse=","),""), workers=n_cores)
  worker <- function(i) {
    record(i, "running")
    timing(i,"running")
    started <- proc.time()[["elapsed"]]
    result <- tryCatch(list(ok = TRUE, value = FUN(X[[i]]), message = "", call = "",
                  elapsed = proc.time()[["elapsed"]] - started),
      error = function(e) list(ok = FALSE, value = NULL,
        message = conditionMessage(e),
        call = paste(deparse(conditionCall(e)), collapse = " "),
        elapsed = proc.time()[["elapsed"]] - started))
    record(i, if (result$ok) "success" else "failed", result$message, result$elapsed)
    timing(i,if(result$ok && isTRUE(attr(result$value,"eivgp_cache_hit"))) "cached" else
      if(result$ok) "complete" else "failed",result$elapsed)
    if(result$ok) attr(result$value,"eivgp_cache_hit") <- NULL
    result
  }
  results <- parallel_map(as.list(seq_along(X)), worker, n_cores = n_cores,
    seeds = seeds, mc.preschedule = mc.preschedule)
  results <- lapply(results, function(x) {
    if (is.list(x) && is.logical(x$ok) && length(x$ok) == 1L && !is.na(x$ok)) return(x)
    list(ok = FALSE, value = NULL, message = if (inherits(x, "try-error"))
      as.character(x) else "Worker returned no valid result (possibly terminated).",
      call = "", elapsed = NA_real_)
  })
  ok <- vapply(results, function(x) isTRUE(x$ok), logical(1))
  status <- data.frame(task = vapply(X, function(x) paste(x, collapse = ","), ""),
    status = ifelse(ok, "success", "failed"),
    message = vapply(results, `[[`, "", "message"),
    call = vapply(results, `[[`, "", "call"),
    elapsed_seconds = vapply(results, `[[`, 0, "elapsed"))
  mixedgp_atomic_write_csv(status, status_path)
  if (any(!ok)) {
    detail <- paste0(study, ": ", sum(!ok), "/", length(ok),
      " replications failed. First error: ", status$message[which(!ok)[1L]],
      ". See ", status_path)
    if (!any(ok)) stop(detail, call. = FALSE)
    warning(detail, ". Aggregating successful replications only; report failures alongside estimates.",
            call. = FALSE)
  }
  lapply(results[ok], `[[`, "value")
}

mixedgp_run_study1_cell <- function(config, cell, engine, run_dir) {
  old_wd <- setwd(config$code_dir)
  on.exit(setwd(old_wd), add = TRUE)
  cell_output <- file.path(run_dir, "cells", cell$id)
  dir.create(cell_output, recursive = TRUE, showWarnings = FALSE)
  run_env <- new.env(parent = engine)
  list2env(
    mixedgp_cell_controls_study1(config, cell, cell_output),
    envir = run_env
  )
  eval(body(eivgp_study1_driver), envir = run_env)
  outputs <- list(
    predictive_metrics = mixedgp_get_cell_output(run_env, "mc_results"),
    mean_recovery = mixedgp_get_cell_output(run_env, "mean_recovery"),
    latent_imputation = mixedgp_get_cell_output(run_env, "latent_imputation"),
    surface_recovery = mixedgp_get_cell_output(run_env, "surface_recovery"),
    ablation_predictive_metrics = mixedgp_get_cell_output(run_env, "ablation_results"),
    ablation_surface_recovery = mixedgp_get_cell_output(
      run_env, "ablation_surface_recovery"
    ),
    mcmc_diagnostics = mixedgp_get_cell_output(run_env, "mcmc_diagnostics"),
    competitor_status = mixedgp_get_cell_output(run_env, "competitor_status"),
    ablation_status = mixedgp_get_cell_output(run_env, "ablation_status"),
    sampler_control_manifest = mixedgp_get_cell_output(
      run_env, "sampler_control_manifest"
    )
  )
  list(
    cell = cell,
    outputs = outputs,
    design_tag = get0("STUDY1_DESIGN_TAG", envir = run_env, inherits = FALSE),
    output_directory = cell_output
  )
}

mixedgp_run_study2_cell <- function(config, cell, engine, run_dir) {
  old_wd <- setwd(config$code_dir)
  on.exit(setwd(old_wd), add = TRUE)
  cell_output <- file.path(run_dir, "cells", cell$id)
  dir.create(cell_output, recursive = TRUE, showWarnings = FALSE)
  run_env <- new.env(parent = engine)
  list2env(
    mixedgp_cell_controls_study2(config, cell, cell_output),
    envir = run_env
  )
  eval(body(eivgp_study2_driver), envir = run_env)
  outputs <- get0(
    "raw_outputs", envir = run_env, inherits = FALSE, ifnotfound = list()
  )
  list(
    cell = cell,
    outputs = outputs,
    design_tag = get0("STUDY2_DESIGN_TAG", envir = run_env, inherits = FALSE),
    output_directory = cell_output
  )
}



# ---- simulations/helpers/gates.R ----
# Simulation gates: package implementation.

mixedgp_competitor_gate <- function(result, config) {
  if (!length(config$published_methods)) return(data.frame())
  status <- result$outputs$competitor_status
  expected <- result$cell$n_rep * length(config$published_methods)
  if (!is.data.frame(status) || nrow(status) == 0L) {
    failure_rate <- 1
    n_success <- 0L
    n_attempted <- 0L
  } else {
    n_attempted <- nrow(status)
    n_success <- sum(status$status == "success", na.rm = TRUE)
    failure_rate <- 1 - n_success / expected
  }
  data.frame(
    cell_id = result$cell$id,
    gate = "published_competitors",
    expected = expected,
    attempted = n_attempted,
    success = n_success,
    failure_rate = failure_rate,
    pass = n_attempted == expected && if (isTRUE(config$strict_competitors)) {
      failure_rate == 0
    } else {
      failure_rate <= 0.05
    },
    stringsAsFactors = FALSE
  )
}

mixedgp_mcmc_gate <- function(result, config) {
  diagnostics <- result$outputs$mcmc_diagnostics
  pass_column <- if (config$study == "study1") "gate_pass" else "mcmc_pass"
  expected <- result$cell$n_rep * length(result$cell$calibration_grid)
  values <- if (is.data.frame(diagnostics) && pass_column %in% names(diagnostics)) {
    as.logical(diagnostics[[pass_column]])
  } else {
    logical(0)
  }
  if (identical(config$study, "study2") && length(values) > 0L) {
    components <- c("target_functional_pass", "invariant_measurement_pass",
                    "raw_coordinate_pass", "n_calib")
    if (!all(components %in% names(diagnostics))) {
      values[] <- FALSE
    } else {
      values <- (values %in% TRUE) &
        (diagnostics$target_functional_pass %in% TRUE) &
        (diagnostics$invariant_measurement_pass %in% TRUE) &
        ((diagnostics$n_calib %in% 0L) | (diagnostics$raw_coordinate_pass %in% TRUE))
    }
  }
  data.frame(
    cell_id = result$cell$id,
    gate = "eiv_mcmc",
    expected = expected,
    attempted = length(values),
    success = sum(values %in% TRUE),
    failure_rate = if (expected > 0L) 1 - sum(values %in% TRUE) / expected else 0,
    pass = length(values) == expected && all(values %in% TRUE),
    stringsAsFactors = FALSE
  )
}

mixedgp_expected_method_keys <- function(cell, method_calibrations) {
  rows <- lapply(names(method_calibrations), function(method) {
    calibration <- method_calibrations[[method]]
    expand.grid(
      rep = seq_len(cell$n_rep),
      method = method,
      n_calib = calibration,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  })
  mixedgp_bind_rows_base(rows)
}

mixedgp_encode_completeness_keys <- function(x, columns) {
  if (!is.data.frame(x) || nrow(x) == 0L ||
      any(!columns %in% names(x))) return(character(0))
  do.call(paste, c(lapply(x[columns], function(z) {
    z <- as.character(z)
    z[is.na(z)] <- "<NA>"
    z
  }), sep = "\r"))
}

mixedgp_completeness_gate_row <- function(cell_id,
                                          gate,
                                          expected,
                                          actual,
                                          columns) {
  expected_key <- mixedgp_encode_completeness_keys(expected, columns)
  actual_key <- mixedgp_encode_completeness_keys(actual, columns)
  expected_unique <- unique(expected_key)
  actual_unique <- unique(actual_key)
  n_success <- length(intersect(expected_unique, actual_unique))
  pass <- length(actual_key) == length(actual_unique) &&
    setequal(expected_unique, actual_unique)
  data.frame(
    cell_id = cell_id,
    gate = gate,
    expected = length(expected_unique),
    attempted = length(actual_key),
    success = n_success,
    failure_rate = if (length(expected_unique) > 0L) {
      1 - n_success / length(expected_unique)
    } else if (length(actual_key) == 0L) {
      0
    } else {
      1
    },
    pass = pass,
    stringsAsFactors = FALSE
  )
}

mixedgp_normalize_rf_method <- function(x) {
  x <- as.character(x)
  x[x %in% c(
    "Response-free threshold model",
    "Response-free measurement model",
    "Ordinal model (no Y)"
  )] <- "RF-OM"
  x
}

mixedgp_output_completeness_gates <- function(result, config) {
  cell <- result$cell
  calibs <- as.integer(cell$calibration_grid)
  positive_calibs <- calibs[calibs > 0L]
  overall <- function(x) {
    if (is.data.frame(x) && "evaluation_stratum" %in% names(x)) {
      x <- x[as.character(x$evaluation_stratum) == "overall", , drop = FALSE]
    }
    x
  }
  filter_expected_methods <- function(x, methods) {
    if (!is.data.frame(x) || nrow(x) == 0L || !"method" %in% names(x)) {
      return(data.frame())
    }
    x$method <- as.character(x$method)
    x[x$method %in% methods, , drop = FALSE]
  }

  central_spec <- c(
    list(`EIV-GP` = calibs),
    setNames(rep(list(NA_integer_), length(config$published_methods)),
             config$published_methods)
  )
  predictive_spec <- c(central_spec, list(Oracle = NA_integer_))
  predictive_expected <- mixedgp_expected_method_keys(cell, predictive_spec)
  predictive_actual <- filter_expected_methods(
    overall(result$outputs$predictive_metrics), names(predictive_spec)
  )
  mean_expected <- mixedgp_expected_method_keys(cell, central_spec)
  mean_actual <- filter_expected_methods(
    result$outputs$mean_recovery, names(central_spec)
  )
  rows <- list(
    mixedgp_completeness_gate_row(
      cell$id, "prediction_output_completeness",
      predictive_expected, predictive_actual, c("rep", "method", "n_calib")
    ),
    mixedgp_completeness_gate_row(
      cell$id, "observable_mean_output_completeness",
      mean_expected, mean_actual, c("rep", "method", "n_calib")
    )
  )

  if (isTRUE(cell$run_ablations)) {
    ablation_prediction_spec <- list(`PI-GP` = calibs, `CC-GP` = positive_calibs)
    ablation_prediction_expected <- mixedgp_expected_method_keys(
      cell, ablation_prediction_spec
    )
    ablation_prediction_actual <- filter_expected_methods(
      overall(result$outputs$ablation_predictive_metrics),
      names(ablation_prediction_spec)
    )
    rows[[length(rows) + 1L]] <- mixedgp_completeness_gate_row(
      cell$id, "ablation_prediction_output_completeness",
      ablation_prediction_expected, ablation_prediction_actual,
      c("rep", "method", "n_calib")
    )
  }

  if (isTRUE(cell$evaluate_f)) {
    eiv_surface_expected <- mixedgp_expected_method_keys(
      cell, list(`EIV-GP` = positive_calibs)
    )
    eiv_surface_actual <- filter_expected_methods(
      result$outputs$surface_recovery, "EIV-GP"
    )
    ablation_surface_spec <- list(
      `PI-GP` = positive_calibs,
      `CC-GP` = positive_calibs,
      `Full-U GP` = NA_integer_
    )
    ablation_surface_expected <- mixedgp_expected_method_keys(
      cell, ablation_surface_spec
    )
    ablation_surface_actual <- filter_expected_methods(
      result$outputs$ablation_surface_recovery,
      names(ablation_surface_spec)
    )
    rows[[length(rows) + 1L]] <- mixedgp_completeness_gate_row(
      cell$id, "latent_surface_eiv_output_completeness",
      eiv_surface_expected, eiv_surface_actual,
      c("rep", "method", "n_calib")
    )
    rows[[length(rows) + 1L]] <- mixedgp_completeness_gate_row(
      cell$id, "latent_surface_benchmark_output_completeness",
      ablation_surface_expected, ablation_surface_actual,
      c("rep", "method", "n_calib")
    )
  }

  if (isTRUE(cell$evaluate_u)) {
    if (config$study == "study1") {
      latent_expected <- mixedgp_expected_method_keys(
        cell, list(`EIV-GP` = calibs, `RF-OM` = calibs)
      )
      latent_actual <- result$outputs$latent_imputation
      if (is.data.frame(latent_actual) && nrow(latent_actual) > 0L) {
        latent_actual$method <- mixedgp_normalize_rf_method(latent_actual$method)
        latent_actual <- latent_actual[
          latent_actual$method %in% c("EIV-GP", "RF-OM"), , drop = FALSE
        ]
      }
      rows[[length(rows) + 1L]] <- mixedgp_completeness_gate_row(
        cell$id, "training_latent_u_output_completeness",
        latent_expected, latent_actual, c("rep", "method", "n_calib")
      )
    } else {
      targets <- c("training_missing_U", "prospective_U_given_C")
      latent_expected <- expand.grid(
        rep = seq_len(cell$n_rep),
        method = c("EIV-GP", "RF-OM"),
        n_calib = calibs,
        target = targets,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
      latent_actual <- result$outputs$latent_imputation_status
      if (is.data.frame(latent_actual) && nrow(latent_actual) > 0L) {
        latent_actual$method <- mixedgp_normalize_rf_method(latent_actual$method)
        latent_actual <- latent_actual[
          latent_actual$method %in% c("EIV-GP", "RF-OM") &
            latent_actual$target %in% targets, , drop = FALSE
        ]
      }
      rows[[length(rows) + 1L]] <- mixedgp_completeness_gate_row(
        cell$id, "latent_u_task_status_completeness",
        latent_expected, latent_actual,
        c("rep", "method", "n_calib", "target")
      )
    }
  }
  mixedgp_bind_rows_base(rows)
}

mixedgp_ablation_gate <- function(result) {
  status <- result$outputs$ablation_status
  if (!isTRUE(result$cell$run_ablations)) {
    return(data.frame(
      cell_id = result$cell$id, gate = "scientific_ablations",
      expected = 0L, attempted = 0L, success = 0L,
      failure_rate = 0, pass = TRUE, stringsAsFactors = FALSE
    ))
  }
  calibs <- as.integer(result$cell$calibration_grid)
  expected_spec <- list(`RF-OM` = calibs, `PI-GP` = calibs, `CC-GP` = calibs)
  if (isTRUE(result$cell$evaluate_f)) {
    expected_spec[["Full-U GP"]] <- NA_integer_
  }
  expected <- mixedgp_expected_method_keys(result$cell, expected_spec)
  if (!is.data.frame(status) || nrow(status) == 0L) {
    status <- data.frame()
  } else {
    status$method <- mixedgp_normalize_rf_method(status$method)
    status <- status[status$method %in% names(expected_spec), , drop = FALSE]
  }
  completeness <- mixedgp_completeness_gate_row(
    result$cell$id, "scientific_ablation_status_completeness",
    expected, status, c("rep", "method", "n_calib")
  )
  acceptable <- is.data.frame(status) && nrow(status) > 0L &&
    all(status$status %in% c("success", "not_applicable"))
  completeness$gate <- "scientific_ablations"
  expected_key <- mixedgp_encode_completeness_keys(
    expected, c("rep", "method", "n_calib")
  )
  acceptable_status <- if (is.data.frame(status) && nrow(status) > 0L) {
    status[status$status %in% c("success", "not_applicable"), , drop = FALSE]
  } else {
    data.frame()
  }
  acceptable_key <- mixedgp_encode_completeness_keys(
    acceptable_status, c("rep", "method", "n_calib")
  )
  completeness$success <- length(intersect(
    unique(expected_key), unique(acceptable_key)
  ))
  completeness$failure_rate <- if (completeness$expected > 0L) {
    1 - completeness$success / completeness$expected
  } else {
    0
  }
  completeness$pass <- isTRUE(completeness$pass) && acceptable
  completeness
}

mixedgp_result_gates <- function(result, config) {
  mixedgp_bind_rows_base(list(
    mixedgp_competitor_gate(result, config),
    mixedgp_mcmc_gate(result, config),
    mixedgp_ablation_gate(result),
    mixedgp_output_completeness_gates(result, config)
  ))
}



# ---- simulations/helpers/summaries.R ----
# Simulation summaries: package implementation.

mixedgp_tag_output <- function(x, cell, study) {
  if (!is.data.frame(x) || nrow(x) == 0L) return(x)
  design <- data.frame(
    study = rep(study, nrow(x)),
    cell_id = rep(cell$id, nrow(x)),
    design_role = rep(cell$role, nrow(x)),
    design_eta = rep(
      if (is.null(cell$heterogeneity_eta)) NA_real_ else cell$heterogeneity_eta,
      nrow(x)
    ),
    design_q = rep(if (is.null(cell$q)) NA_integer_ else cell$q, nrow(x)),
    stringsAsFactors = FALSE
  )
  duplicate <- intersect(names(design), names(x))
  if (length(duplicate) > 0L) design[duplicate] <- NULL
  cbind(design, x)
}

mixedgp_group_summary <- function(x, value, groups) {
  if (!is.data.frame(x) || nrow(x) == 0L ||
      !value %in% names(x)) return(data.frame())
  groups <- intersect(groups, names(x))
  if (length(groups) == 0L) stop("At least one grouping column is required.")
  group_key <- do.call(
    paste,
    c(lapply(x[groups], function(z) {
      z <- as.character(z)
      z[is.na(z)] <- "<NA>"
      z
    }), sep = "\r")
  )
  index <- split(seq_len(nrow(x)), group_key, drop = TRUE)
  rows <- lapply(index, function(ii) {
    values <- as.numeric(x[[value]][ii])
    values <- values[is.finite(values)]
    n_eff <- length(values)
    estimate <- if (n_eff == 0L) NA_real_ else mean(values)
    se <- if (n_eff < 2L) NA_real_ else stats::sd(values) / sqrt(n_eff)
    cbind(
      x[ii[[1L]], groups, drop = FALSE],
      data.frame(
        mean = estimate,
        se = se,
        ci_lower = if (is.finite(se)) estimate - 1.96 * se else NA_real_,
        ci_upper = if (is.finite(se)) estimate + 1.96 * se else NA_real_,
        n_pairs = n_eff,
        stringsAsFactors = FALSE
      )
    )
  })
  mixedgp_bind_rows_base(rows)
}

mixedgp_paired_eiv_published <- function(x, task, metrics) {
  required <- c("cell_id", "rep", "scenario", "n_calib", "method")
  if (!is.data.frame(x) || nrow(x) == 0L ||
      any(!required %in% names(x))) return(data.frame())
  metrics <- intersect(metrics, names(x))
  if (length(metrics) == 0L) return(data.frame())
  x$method <- as.character(x$method)
  eiv <- x[x$method == "EIV-GP" & !is.na(x$n_calib), , drop = FALSE]
  competitor <- x[
    x$method %in% MIXEDGP_PUBLISHED_METHODS & is.na(x$n_calib),
    , drop = FALSE
  ]
  if (nrow(eiv) == 0L || nrow(competitor) == 0L) return(data.frame())
  keys <- intersect(
    c(
      "study", "cell_id", "design_role", "design_eta", "design_q",
      "rep", "scenario", "evaluation_stratum"
    ),
    intersect(names(eiv), names(competitor))
  )
  eiv_keep <- unique(c(keys, "n_calib", metrics))
  competitor_keep <- unique(c(keys, "method", metrics))
  eiv <- eiv[, eiv_keep, drop = FALSE]
  competitor <- competitor[, competitor_keep, drop = FALSE]
  names(competitor)[names(competitor) == "method"] <- "competitor"
  names(eiv)[match(metrics, names(eiv))] <- paste0("EIV_", metrics)
  names(competitor)[match(metrics, names(competitor))] <-
    paste0("Competitor_", metrics)
  paired <- merge(eiv, competitor, by = keys, all = FALSE, sort = FALSE)
  if (nrow(paired) == 0L) return(data.frame())
  if (!"evaluation_stratum" %in% names(paired)) {
    paired$evaluation_stratum <- "overall"
  }
  rows <- lapply(metrics, function(metric) {
    eiv_value <- as.numeric(paired[[paste0("EIV_", metric)]])
    competitor_value <- as.numeric(
      paired[[paste0("Competitor_", metric)]]
    )
    id_columns <- setdiff(
      names(paired),
      c(paste0("EIV_", metrics), paste0("Competitor_", metrics))
    )
    cbind(
      paired[, id_columns, drop = FALSE],
      data.frame(
        task = task,
        metric = metric,
        EIV_value = eiv_value,
        competitor_value = competitor_value,
        EIV_advantage = competitor_value - eiv_value,
        advantage_direction = "positive favors EIV-GP",
        stringsAsFactors = FALSE
      )
    )
  })
  mixedgp_bind_rows_base(rows)
}

mixedgp_paired_eiv_same_calibration <- function(eiv_x,
                                                 competitor_x,
                                                 task,
                                                 methods,
                                                 metrics) {
  required <- c("cell_id", "rep", "scenario", "n_calib", "method")
  if (!is.data.frame(eiv_x) || !is.data.frame(competitor_x) ||
      nrow(eiv_x) == 0L || nrow(competitor_x) == 0L ||
      any(!required %in% names(eiv_x)) ||
      any(!required %in% names(competitor_x))) return(data.frame())
  metrics <- intersect(metrics, intersect(names(eiv_x), names(competitor_x)))
  if (length(metrics) == 0L) return(data.frame())
  eiv_x$method <- as.character(eiv_x$method)
  competitor_x$method <- as.character(competitor_x$method)
  eiv <- eiv_x[eiv_x$method == "EIV-GP" & !is.na(eiv_x$n_calib), , drop = FALSE]
  competitor <- competitor_x[
    competitor_x$method %in% methods & !is.na(competitor_x$n_calib),
    , drop = FALSE
  ]
  if (nrow(eiv) == 0L || nrow(competitor) == 0L) return(data.frame())
  keys <- intersect(
    c(
      "study", "cell_id", "design_role", "design_eta", "design_q",
      "rep", "scenario", "n_calib", "evaluation_stratum"
    ),
    intersect(names(eiv), names(competitor))
  )
  eiv <- eiv[, unique(c(keys, metrics)), drop = FALSE]
  competitor <- competitor[, unique(c(keys, "method", metrics)), drop = FALSE]
  names(competitor)[names(competitor) == "method"] <- "competitor"
  names(eiv)[match(metrics, names(eiv))] <- paste0("EIV_", metrics)
  names(competitor)[match(metrics, names(competitor))] <-
    paste0("Competitor_", metrics)
  paired <- merge(eiv, competitor, by = keys, all = FALSE, sort = FALSE)
  if (nrow(paired) == 0L) return(data.frame())
  if (!"evaluation_stratum" %in% names(paired)) {
    paired$evaluation_stratum <- "overall"
  }
  rows <- lapply(metrics, function(metric) {
    eiv_value <- as.numeric(paired[[paste0("EIV_", metric)]])
    competitor_value <- as.numeric(
      paired[[paste0("Competitor_", metric)]]
    )
    id_columns <- setdiff(
      names(paired),
      c(paste0("EIV_", metrics), paste0("Competitor_", metrics))
    )
    cbind(
      paired[, id_columns, drop = FALSE],
      data.frame(
        task = task,
        metric = metric,
        EIV_value = eiv_value,
        competitor_value = competitor_value,
        EIV_advantage = competitor_value - eiv_value,
        advantage_direction = "positive favors EIV-GP",
        stringsAsFactors = FALSE
      )
    )
  })
  mixedgp_bind_rows_base(rows)
}

mixedgp_paired_latent_u <- function(x) {
  required <- c(
    "cell_id", "rep", "scenario", "n_calib", "method", "RMSE", "MAE"
  )
  if (!is.data.frame(x) || nrow(x) == 0L ||
      any(!required %in% names(x))) return(data.frame())
  x$method <- as.character(x$method)
  response_free_names <- c(
    "Response-free threshold model", "Ordinal model (no Y)"
  )
  if ("score_status" %in% names(x)) {
    x <- x[is.na(x$score_status) | x$score_status == "scored", , drop = FALSE]
  }
  eiv <- x[x$method == "EIV-GP", , drop = FALSE]
  response_free <- x[x$method %in% response_free_names, , drop = FALSE]
  if (nrow(eiv) == 0L || nrow(response_free) == 0L) return(data.frame())
  keys <- intersect(
    c(
      "study", "cell_id", "design_role", "design_eta", "design_q",
      "rep", "scenario", "n_calib", "target", "coordinate"
    ),
    intersect(names(eiv), names(response_free))
  )
  metrics <- c("RMSE", "MAE")
  eiv <- eiv[, unique(c(keys, metrics)), drop = FALSE]
  response_free <- response_free[, unique(c(keys, "method", metrics)), drop = FALSE]
  names(response_free)[names(response_free) == "method"] <- "competitor"
  names(eiv)[match(metrics, names(eiv))] <- paste0("EIV_", metrics)
  names(response_free)[match(metrics, names(response_free))] <-
    paste0("Competitor_", metrics)
  paired <- merge(eiv, response_free, by = keys, all = FALSE, sort = FALSE)
  if (nrow(paired) == 0L) return(data.frame())
  rows <- lapply(metrics, function(metric) {
    eiv_value <- as.numeric(paired[[paste0("EIV_", metric)]])
    competitor_value <- as.numeric(
      paired[[paste0("Competitor_", metric)]]
    )
    id_columns <- setdiff(
      names(paired),
      c(paste0("EIV_", metrics), paste0("Competitor_", metrics))
    )
    cbind(
      paired[, id_columns, drop = FALSE],
      data.frame(
        task = "latent U inference",
        metric = metric,
        EIV_value = eiv_value,
        competitor_value = competitor_value,
        EIV_advantage = competitor_value - eiv_value,
        advantage_direction = "positive favors EIV-GP",
        stringsAsFactors = FALSE
      )
    )
  })
  mixedgp_bind_rows_base(rows)
}

mixedgp_paired_latent_surface <- function(eiv, ablations) {
  required <- c("cell_id", "rep", "scenario", "n_calib", "method", "ISE")
  if (!is.data.frame(eiv) || !is.data.frame(ablations) ||
      nrow(eiv) == 0L || nrow(ablations) == 0L ||
      any(!required %in% names(eiv)) || any(!required %in% names(ablations))) {
    return(data.frame())
  }
  eiv$method <- as.character(eiv$method)
  ablations$method <- as.character(ablations$method)
  eiv <- eiv[eiv$method == "EIV-GP" & !is.na(eiv$n_calib), , drop = FALSE]
  ablations <- ablations[
    ablations$method %in% c("PI-GP", "CC-GP", "Full-U GP"),
    , drop = FALSE
  ]
  if (nrow(eiv) == 0L || nrow(ablations) == 0L) return(data.frame())
  keys <- intersect(
    c(
      "study", "cell_id", "design_role", "design_eta", "design_q",
      "rep", "scenario"
    ),
    intersect(names(eiv), names(ablations))
  )
  eiv <- eiv[, unique(c(keys, "n_calib", "ISE")), drop = FALSE]
  names(eiv)[names(eiv) == "ISE"] <- "EIV_ISE"
  same_k <- ablations[!is.na(ablations$n_calib), , drop = FALSE]
  benchmark <- ablations[is.na(ablations$n_calib), , drop = FALSE]
  paired_same <- data.frame()
  if (nrow(same_k) > 0L) {
    paired_same <- merge(
      eiv,
      same_k[, unique(c(keys, "n_calib", "method", "ISE")), drop = FALSE],
      by = c(keys, "n_calib"), all = FALSE, sort = FALSE
    )
  }
  paired_benchmark <- data.frame()
  if (nrow(benchmark) > 0L) {
    paired_benchmark <- merge(
      eiv,
      benchmark[, unique(c(keys, "method", "ISE")), drop = FALSE],
      by = keys, all = FALSE, sort = FALSE
    )
  }
  paired <- mixedgp_bind_rows_base(list(paired_same, paired_benchmark))
  if (nrow(paired) == 0L) return(data.frame())
  names(paired)[names(paired) == "method"] <- "competitor"
  names(paired)[names(paired) == "ISE"] <- "competitor_value"
  paired$task <- "latent surface f(x,u)"
  paired$metric <- "ISE"
  paired$EIV_value <- paired$EIV_ISE
  paired$EIV_advantage <- paired$competitor_value - paired$EIV_value
  paired$advantage_direction <- "positive favors EIV-GP"
  paired$EIV_ISE <- NULL
  paired
}

mixedgp_method_comparisons <- function(combined) {
  predictive <- combined$predictive_metrics
  ablation_predictive <- combined$ablation_predictive_metrics
  mean_recovery <- combined$mean_recovery
  latent_u <- combined$latent_imputation
  surface <- combined$surface_recovery
  ablation_surface <- combined$ablation_surface_recovery
  mixedgp_bind_rows_base(list(
    mixedgp_paired_eiv_published(
      predictive,
      task = "response prediction Y*|x*,c*",
      metrics = c("RMSE", "MAE", "CRPS", "NLPD", "IntervalScore95")
    ),
    mixedgp_paired_eiv_same_calibration(
      predictive,
      ablation_predictive,
      task = "response prediction Y*|x*,c*",
      methods = c("PI-GP", "CC-GP"),
      metrics = c("RMSE", "MAE", "CRPS", "NLPD", "IntervalScore95")
    ),
    mixedgp_paired_eiv_published(
      mean_recovery,
      task = "observable mean m(x,c)",
      metrics = c("RMSE", "MAE")
    ),
    mixedgp_paired_latent_u(latent_u),
    mixedgp_paired_latent_surface(surface, ablation_surface)
  ))
}

mixedgp_study1_mechanism_contrasts <- function(comparisons, threshold_design = NULL) {
  if (is.null(threshold_design)) {
    return(mixedgp_bind_rows_base(lapply(c("balanced", "imbalanced"), function(design) {
      mixedgp_study1_mechanism_contrasts(comparisons, design)
    })))
  }
  if (!is.data.frame(comparisons) || nrow(comparisons) == 0L ||
      !"cell_id" %in% names(comparisons)) return(data.frame())
  control <- comparisons[
    comparisons$cell_id == paste0("eta0_", threshold_design) &
      comparisons$task %in% c(
        "response prediction Y*|x*,c*", "observable mean m(x,c)"
      ), , drop = FALSE
  ]
  active <- comparisons[
    comparisons$cell_id == paste0("eta1_", threshold_design) &
      comparisons$task %in% c(
        "response prediction Y*|x*,c*", "observable mean m(x,c)"
      ), , drop = FALSE
  ]
  keys <- intersect(
    c(
      "rep", "n_calib", "task", "competitor", "metric",
      "evaluation_stratum"
    ),
    intersect(names(control), names(active))
  )
  if (nrow(control) == 0L || nrow(active) == 0L) return(data.frame())
  control <- control[, c(keys, "EIV_advantage"), drop = FALSE]
  active <- active[, c(keys, "EIV_advantage"), drop = FALSE]
  names(control)[names(control) == "EIV_advantage"] <- "advantage_eta0"
  names(active)[names(active) == "EIV_advantage"] <- "advantage_eta1"
  out <- merge(control, active, by = keys, all = FALSE, sort = FALSE)
  out$contrast <- paste0("heterogeneity (", threshold_design, "): eta=1 minus eta=0")
  out$advantage_change <- out$advantage_eta1 - out$advantage_eta0
  out$contrast_direction <- "positive means heterogeneity strengthens EIV advantage"
  out
}

mixedgp_study2_design_contrasts <- function(comparisons) {
  if (!is.data.frame(comparisons) || nrow(comparisons) == 0L ||
      !"cell_id" %in% names(comparisons)) return(data.frame())
  x <- comparisons[
    comparisons$task %in% c(
      "response prediction Y*|x*,c*", "observable mean m(x,c)"
    ), , drop = FALSE
  ]
  if ("evaluation_stratum" %in% names(x)) {
    x <- x[x$evaluation_stratum == "overall", , drop = FALSE]
  }
  keys <- intersect(
    c("rep", "n_calib", "task", "competitor", "metric"), names(x)
  )
  extract_advantage <- function(cell_id, label) {
    out <- x[x$cell_id == cell_id, c(keys, "EIV_advantage"), drop = FALSE]
    names(out)[names(out) == "EIV_advantage"] <- label
    out
  }
  primary <- extract_advantage("primary_q4_calibration", "advantage_primary")
  additive <- extract_advantage("additive_q4", "advantage_additive")
  uncertain <- extract_advantage(
    "high_uncertainty_q4", "advantage_high_uncertainty"
  )
  controls <- list()
  if (nrow(primary) > 0L && nrow(additive) > 0L) {
    z <- merge(primary, additive, by = keys, all = FALSE, sort = FALSE)
    z$contrast <- "interactions: primary minus additive"
    z$advantage_change <- z$advantage_primary - z$advantage_additive
    controls[[length(controls) + 1L]] <- z
  }
  if (nrow(primary) > 0L && nrow(uncertain) > 0L) {
    z <- merge(uncertain, primary, by = keys, all = FALSE, sort = FALSE)
    z$contrast <- "latent uncertainty: high minus primary"
    z$advantage_change <-
      z$advantage_high_uncertainty - z$advantage_primary
    controls[[length(controls) + 1L]] <- z
  }
  primary_ids <- c("primary_q2", "primary_q3", "primary_q4_calibration")
  q_rows <- x[x$cell_id %in% primary_ids, , drop = FALSE]
  q_rows$design_q <- as.integer(q_rows$design_q)
  q_keys <- intersect(
    c("rep", "n_calib", "task", "competitor", "metric"), names(q_rows)
  )
  for (pair in list(c(2L, 3L), c(3L, 4L), c(2L, 4L))) {
    low <- q_rows[q_rows$design_q == pair[[1L]],
                  c(q_keys, "EIV_advantage"), drop = FALSE]
    high <- q_rows[q_rows$design_q == pair[[2L]],
                   c(q_keys, "EIV_advantage"), drop = FALSE]
    names(low)[names(low) == "EIV_advantage"] <- "advantage_low_q"
    names(high)[names(high) == "EIV_advantage"] <- "advantage_high_q"
    if (nrow(low) > 0L && nrow(high) > 0L) {
      z <- merge(low, high, by = q_keys, all = FALSE, sort = FALSE)
      z$contrast <- paste0("proxy dimension: q=", pair[[2L]],
                           " minus q=", pair[[1L]])
      z$advantage_change <- z$advantage_high_q - z$advantage_low_q
      controls[[length(controls) + 1L]] <- z
    }
  }
  out <- mixedgp_bind_rows_base(controls)
  if (nrow(out) > 0L) {
    out$contrast_direction <-
      "positive means the named design change strengthens EIV advantage"
  }
  out
}

mixedgp_write_publication_comparisons <- function(combined, config, combined_dir) {
  comparisons <- mixedgp_method_comparisons(combined)
  if (nrow(comparisons) > 0L) {
    mixedgp_atomic_write_csv(
      comparisons, file.path(combined_dir, "method_comparisons_paired_raw.csv")
    )
    comparison_groups <- c(
      "study", "cell_id", "design_role", "design_eta", "design_q",
      "task", "n_calib", "competitor", "metric", "evaluation_stratum",
      "target", "coordinate", "advantage_direction"
    )
    comparison_summary <- mixedgp_group_summary(
      comparisons, "EIV_advantage", comparison_groups
    )
    mixedgp_atomic_write_csv(
      comparison_summary,
      file.path(combined_dir, "method_comparisons_paired_summary.csv")
    )
  }
  contrasts <- if (config$study == "study1") {
    mixedgp_study1_mechanism_contrasts(comparisons)
  } else {
    mixedgp_study2_design_contrasts(comparisons)
  }
  if (nrow(contrasts) > 0L) {
    mixedgp_atomic_write_csv(
      contrasts, file.path(combined_dir, "design_contrasts_paired_raw.csv")
    )
    contrast_summary <- mixedgp_group_summary(
      contrasts, "advantage_change",
      c(
        "task", "n_calib", "competitor", "metric",
        "evaluation_stratum", "contrast", "contrast_direction"
      )
    )
    mixedgp_atomic_write_csv(
      contrast_summary,
      file.path(combined_dir, "design_contrasts_paired_summary.csv")
    )
  }
  invisible(list(comparisons = comparisons, contrasts = contrasts))
}

mixedgp_aggregate_results <- function(results, config, run_dir) {
  combined_dir <- file.path(run_dir, "combined")
  dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
  non_tabular <- lapply(results, function(result) {
    result$outputs[!vapply(result$outputs, is.data.frame, logical(1L))]
  })
  non_tabular <- Filter(function(x) length(x) > 0L, non_tabular)
  if (length(non_tabular) > 0L) {
    mixedgp_atomic_save_rds(
      non_tabular, file.path(combined_dir, "cell_non_tabular_outputs.rds")
    )
  }
  output_names <- unique(unlist(lapply(results, function(z) names(z$outputs))))
  combined <- setNames(vector("list", length(output_names)), output_names)
  for (name in output_names) {
    rows <- lapply(results, function(result) {
      mixedgp_tag_output(result$outputs[[name]], result$cell, config$study)
    })
    combined[[name]] <- mixedgp_bind_rows_base(rows)
    if (is.data.frame(combined[[name]]) && ncol(combined[[name]]) > 0L) {
      mixedgp_atomic_write_csv(
        combined[[name]], file.path(combined_dir, paste0(name, ".csv"))
      )
    }
  }
  mixedgp_atomic_save_rds(combined, file.path(combined_dir, "all_raw_outputs.rds"))
  mixedgp_write_publication_comparisons(combined, config, combined_dir)
  output_files <- setdiff(
    list.files(combined_dir, full.names = TRUE),
    file.path(combined_dir, "result_manifest.csv")
  )
  manifest <- data.frame(
    file = basename(output_files),
    bytes = file.info(output_files)$size,
    md5 = unname(tools::md5sum(output_files)),
    stringsAsFactors = FALSE
  )
  mixedgp_atomic_write_csv(manifest, file.path(combined_dir, "result_manifest.csv"))
  combined
}



# ---- simulations/helpers/workflow.R ----
# Simulation workflow: package implementation.

mixedgp_run_simulation <- function(config) {
  config <- validate_simulation_config(config)
  engine <- mixedgp_simulation_engine(config$code_dir)
  preflight_fun <- get("mixedgp_competitor_preflight", envir = engine)
  preflight <- preflight_fun(config$published_methods, strict = FALSE)
  runtime_preflight <- mixedgp_runtime_preflight(config)
  if (any(!preflight$available) || any(!runtime_preflight$available)) {
    message("Install missing experiment dependencies explicitly from the repository root:\n",
      "  install.packages(c('kergp', 'LVGP', 'EzGP'))\n",
      "Use the same Rscript and MIXEDGP_R_LIBRARY setting for setup and execution.\n",
      "No packages are installed automatically. Existing missing-method caches are not repaired by installation.")
  }
  task_plan <- mixedgp_task_plan(config)
  if (identical(config$mode, "dry_run")) {
    mixedgp_print_dry_run(config, preflight, runtime_preflight, task_plan)
    return(invisible(list(
      config = config, preflight = preflight,
      runtime_preflight = runtime_preflight, task_plan = task_plan,
      data_status = mixedgp_existing_data_status(config)
    )))
  }

  run_dir <- mixedgp_run_directory(config)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  mixedgp_write_run_metadata(
    config, run_dir, preflight, runtime_preflight, task_plan
  )
  if (identical(config$mode, "development")) {
    writeLines(config$development_note, file.path(run_dir, "DEVELOPMENT_RESULTS.txt"))
    message(config$development_note)
  }
  if ("fit" %in% config$stages && any(!preflight$available)) {
    unavailable <- preflight$method[!preflight$available]
    warning(
      "Competitor packages unavailable: ", paste(unavailable, collapse = ", "),
      ". Their missing results will be recorded; other methods may proceed.",
      call. = FALSE
    )
  }
  if ("fit" %in% config$stages && any(!runtime_preflight$available)) {
    unavailable <- runtime_preflight$package[!runtime_preflight$available]
    stop(
      "Run stopped before fitting because runtime packages are unavailable: ",
      paste(unavailable, collapse = ", "), "."
    )
  }

  old_cores <- getOption("mixedgp.cores", NULL)
  options(mixedgp.cores = config$parallel$workers)
  on.exit(options(mixedgp.cores = old_cores), add = TRUE)
  thread_vars <- c(
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"
  )
  old_threads <- Sys.getenv(thread_vars, unset = NA_character_)
  do.call(Sys.setenv, as.list(setNames(rep("1", length(thread_vars)), thread_vars)))
  on.exit({
    present <- !is.na(old_threads)
    if (any(present)) {
      do.call(Sys.setenv, as.list(old_threads[present]))
    }
    if (any(!present)) Sys.unsetenv(thread_vars[!present])
  }, add = TRUE)

  generation_manifests <- list()
  if ("data" %in% config$stages) {
    for (cell in config$cells) {
      message("Freezing synthetic data: ", cell$id)
      generation_manifests[[cell$id]] <- mixedgp_generate_cell_data(
        config, cell, engine
      )
    }
  }
  paired_validation <- data.frame()
  selected_manifests <- list()
  if (any(c("data", "fit") %in% config$stages)) {
    for (cell in config$cells) {
      selected <- mixedgp_verify_cell_data(config, cell, engine)
      selected$cell_id <- cell$id
      selected_manifests[[cell$id]] <- selected
    }
    paired_validation <- mixedgp_validate_common_random_numbers(config, engine)
    if (is.data.frame(paired_validation) && ncol(paired_validation) > 0L) {
      mixedgp_atomic_write_csv(
        paired_validation,
        file.path(run_dir, "config", "paired_design_validation.csv")
      )
    }
  }
  data_manifest <- mixedgp_bind_rows_base(selected_manifests)
  if (nrow(data_manifest) > 0L) {
    mixedgp_atomic_write_csv(
      data_manifest, file.path(run_dir, "config", "input_manifest.csv")
    )
  }

  if (!"fit" %in% config$stages) {
    return(invisible(list(
      config = config, run_dir = run_dir, data_manifest = data_manifest,
      preflight = preflight, runtime_preflight = runtime_preflight,
      paired_validation = paired_validation
    )))
  }
  results <- vector("list", length(config$cells))
  names(results) <- vapply(config$cells, `[[`, character(1L), "id")
  gate_rows <- list()
  status_rows <- list()
  for (ii in seq_along(config$cells)) {
    cell <- config$cells[[ii]]
    message(
      "Running ", config$study, " cell ", ii, "/", length(config$cells),
      ": ", cell$id
    )
    started <- Sys.time()
    answer <- tryCatch(
      if (config$study == "study1") {
        mixedgp_run_study1_cell(config, cell, engine, run_dir)
      } else {
        mixedgp_run_study2_cell(config, cell, engine, run_dir)
      },
      error = function(e) e
    )
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (inherits(answer, "error")) {
      status_rows[[cell$id]] <- data.frame(
        cell_id = cell$id, status = "failed",
        message = conditionMessage(answer), elapsed_seconds = elapsed,
        stringsAsFactors = FALSE
      )
      status <- mixedgp_bind_rows_base(status_rows)
      mixedgp_atomic_write_csv(
        status, file.path(run_dir, "config", "cell_status.csv")
      )
      stop("Simulation cell ", cell$id, " failed: ", conditionMessage(answer))
    }
    results[[cell$id]] <- answer
    gate_rows[[cell$id]] <- mixedgp_result_gates(answer, config)
    failed_gates <- gate_rows[[cell$id]]$gate[!gate_rows[[cell$id]]$pass]
    gate_failed <- length(failed_gates) > 0L
    gate_rows[[cell$id]]$advice <- ifelse(gate_rows[[cell$id]]$pass, "",
      mixedgp_simulation_diagnostic_advice())
    if (gate_failed) warning("Cell ", cell$id, " has flagged results: ",
      paste(failed_gates, collapse = "; "), ". ",
      mixedgp_simulation_diagnostic_advice(), call. = FALSE)
    cell_status <- if (gate_failed && isTRUE(config$fail_closed)) {
      "gate_failed"
    } else if (gate_failed) {
      "completed_with_gate_warnings"
    } else {
      "success"
    }
    status_rows[[cell$id]] <- data.frame(
      cell_id = cell$id,
      status = cell_status,
      message = if (gate_failed) paste(failed_gates, collapse = "; ") else "",
      elapsed_seconds = elapsed, stringsAsFactors = FALSE
    )
    mixedgp_atomic_write_csv(
      mixedgp_bind_rows_base(status_rows),
      file.path(run_dir, "config", "cell_status.csv")
    )
    mixedgp_atomic_write_csv(
      mixedgp_bind_rows_base(gate_rows),
      file.path(run_dir, "config", "diagnostic_gates.csv")
    )
    if (isTRUE(config$fail_closed) && gate_failed) {
      stop(
        "Publication gates failed for ", cell$id, ": ",
        paste(failed_gates, collapse = ", "), "."
      )
    }
  }

  combined <- if ("aggregate" %in% config$stages) {
    mixedgp_aggregate_results(results, config, run_dir)
  } else {
    list()
  }
  out <- structure(
    list(
      config = config,
      run_dir = run_dir,
      preflight = preflight,
      runtime_preflight = runtime_preflight,
      data_manifest = data_manifest,
      paired_validation = paired_validation,
      gates = mixedgp_bind_rows_base(gate_rows),
      results = results,
      combined = combined
    ),
    class = c("mixedgp_simulation_run", "list")
  )
  mixedgp_atomic_save_rds(
    out[c(
      "config", "run_dir", "preflight", "runtime_preflight", "data_manifest",
      "paired_validation", "gates"
    )],
    file.path(run_dir, "run_summary.rds")
  )
  message("Completed simulation run: ", normalizePath(run_dir))
  if ("aggregate" %in% config$stages) {
    report_env <- new.env(parent = environment())
    report_env$mixedgp_report_run <- mixedgp_report_run
    report_status <- tryCatch(report_env$mixedgp_report_run(run_dir, "report", code_dir = config$code_dir),
      error = function(e) warning("Fitting completed, but reporting failed: ", conditionMessage(e),
        ". Use run_numerical_study with action=report to retry without fitting.", call. = FALSE))
  }
  if(identical(config$study,"study1") && exists("report_status",inherits=FALSE) &&
     is.data.frame(report_status) && all(report_status$status=="success") &&
     all(vapply(results,function(z)"scientific_ess_pass" %in% names(z$outputs$mcmc_diagnostics) &&
       all(z$outputs$mcmc_diagnostics$scientific_ess_pass %in% TRUE),logical(1))))
    clean_study1_checkpoints(run_dir)
  invisible(out)
}

run_study1_simulation <- function(config = study1_simulation_config()) {
  if (!identical(config$study, "study1")) stop("Expected a Study I config.")
  mixedgp_run_simulation(config)
}

run_study2_simulation <- function(config = study2_simulation_config()) {
  if (!identical(config$study, "study2")) stop("Expected a Study II config.")
  mixedgp_run_simulation(config)
}

## Stable content identity shared by independent fitting/reporting layers.
mixedgp_dataset_identity <- function(data) {
  path <- tempfile(); on.exit(unlink(path), add=TRUE)
  saveRDS(data,path,version=3L,compress=FALSE)
  unname(tools::md5sum(path))
}


# ---- simulations/study1/setup_study1_experiment.R ----
eivgp_study1_setup <- function() {
if (!exists("STUDY1_RUN_ABLATIONS")) STUDY1_RUN_ABLATIONS <- TRUE
if (!exists("STUDY1_EVALUATE_F")) STUDY1_EVALUATE_F <- TRUE
if (!exists("STUDY1_EVALUATE_U")) STUDY1_EVALUATE_U <- TRUE
needed_pkgs <- c("ggplot2", "dplyr", "tidyr", "knitr")
missing_pkgs <- needed_pkgs[
  !vapply(needed_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0L) {
  stop("Please install required packages: ", paste(missing_pkgs, collapse = ", "))
}
suppressPackageStartupMessages({
})

if (!exists("STUDY1_CONFIG")) {
  STUDY1_CONFIG <- if (exists("STUDY1_QUICK") && isTRUE(STUDY1_QUICK)) {
    "quick"
  } else {
    "thorough"
  }
}
if (!STUDY1_CONFIG %in% c("quick", "balanced", "thorough")) {
  stop("STUDY1_CONFIG must be 'quick', 'balanced', or 'thorough'.")
}
STUDY1_QUICK <- identical(STUDY1_CONFIG, "quick")
if (!exists("STUDY1_USE_CACHE")) STUDY1_USE_CACHE <- TRUE
if (!exists("STUDY1_OUT_PREFIX")) STUDY1_OUT_PREFIX <- ".."
if (!exists("STUDY1_REUSE_LOCKED_EIV")) {
  STUDY1_REUSE_LOCKED_EIV <- FALSE
}
if (!exists("STUDY1_REQUIRE_MCMC_GATE")) {
  STUDY1_REQUIRE_MCMC_GATE <- !identical(STUDY1_CONFIG, "quick")
}
if (!exists("STUDY1_MAX_RHAT")) STUDY1_MAX_RHAT <- 1.01
if (!exists("STUDY1_MIN_ESS")) STUDY1_MIN_ESS <- 400
if (!exists("STUDY1_DIAGNOSTIC_MAX_DRAWS")) STUDY1_DIAGNOSTIC_MAX_DRAWS <- Inf
if (!exists("STUDY1_DIAGNOSTIC_N_LATENT")) STUDY1_DIAGNOSTIC_N_LATENT <- 64L
if (!exists("STUDY1_DIAGNOSTIC_N_POINTS")) STUDY1_DIAGNOSTIC_N_POINTS <- 5L
if (!exists("STUDY1_STRICT_COMPETITORS")) {
  STUDY1_STRICT_COMPETITORS <- identical(STUDY1_CONFIG, "thorough")
}
STUDY1_PUBLISHED_COMPETITORS <- character() # Fitted only by the standalone runner.
if (!exists("STUDY1_MECHANISM_CALIB")) STUDY1_MECHANISM_CALIB <- 20L
if (!exists("STUDY1_DATA_DIR")) {
  STUDY1_DATA_DIR <- file.path("..", "data-synthetic", "study1")
}

settings <- switch(
  STUDY1_CONFIG,
  quick = list(
    n_test = 150L, n_rep = 2L, n_iter = 600L, burn = 200L,
    n_chains = 1L, preset = "fast", n_pred_draw = 150L,
    n_m_eval = 40L, n_m_draw = 60L, n_m_latent = 64L
  ),
  balanced = list(
    n_test = 300L, n_rep = 10L, n_iter = 2500L, burn = 750L,
    n_chains = 4L, preset = "balanced", n_pred_draw = 300L,
    n_m_eval = 80L, n_m_draw = 150L, n_m_latent = 128L
  ),
  thorough = list(
    n_test = 500L, n_rep = 50L, n_iter = 5000L, burn = 1000L,
    n_chains = 12L, preset = "balanced", n_pred_draw = 600L,
    n_m_eval = 120L, n_m_draw = 200L, n_m_latent = 256L
  )
)

## Publication masters may override the run size without changing the
## scientific engine below. Keeping these controls explicit also prevents a
## smoke run from being mistaken for a reportable Monte Carlo experiment.
study1_setting_override <- function(name, fallback, minimum = 1L) {
  value <- get0(name, inherits = TRUE, ifnotfound = fallback)
  value <- as.integer(value)
  if (length(value) != 1L || is.na(value) || value < minimum) {
    stop(name, " must be one integer at least ", minimum, ".")
  }
  value
}
settings$n_test <- study1_setting_override(
  "STUDY1_MC_N_TEST", settings$n_test
)
settings$n_rep <- study1_setting_override(
  "STUDY1_MC_N_REP", settings$n_rep
)
settings$n_iter <- study1_setting_override(
  "STUDY1_MC_N_ITER", settings$n_iter, minimum = 2L
)
settings$burn <- study1_setting_override(
  "STUDY1_MC_BURN", settings$burn, minimum = 0L
)
settings$n_chains <- study1_setting_override(
  "STUDY1_MC_N_CHAINS", settings$n_chains
)
settings$n_pred_draw <- study1_setting_override(
  "STUDY1_MC_N_PRED_DRAW", settings$n_pred_draw
)
settings$n_m_eval <- study1_setting_override(
  "STUDY1_MC_N_M_EVAL", settings$n_m_eval
)
settings$n_m_draw <- study1_setting_override(
  "STUDY1_MC_N_M_DRAW", settings$n_m_draw
)
settings$n_m_latent <- study1_setting_override(
  "STUDY1_MC_N_M_LATENT", settings$n_m_latent
)
if (settings$burn >= settings$n_iter) {
  stop("STUDY1_MC_BURN must be smaller than STUDY1_MC_N_ITER.")
}

if (!exists("STUDY1_SCENARIO")) {
  STUDY1_SCENARIO <- "heterogeneity_continuum"
}
if (!STUDY1_SCENARIO %in% c(
  "active", "inactive", "category_sufficient", "heterogeneity_continuum"
)) {
  stop("STUDY1_SCENARIO is not supported by simulate_1d_data().")
}
if (!exists("STUDY1_HETEROGENEITY_ETA")) STUDY1_HETEROGENEITY_ETA <- 1
STUDY1_HETEROGENEITY_ETA <- as.numeric(STUDY1_HETEROGENEITY_ETA)
if (length(STUDY1_HETEROGENEITY_ETA) != 1L ||
    !is.finite(STUDY1_HETEROGENEITY_ETA) ||
    STUDY1_HETEROGENEITY_ETA < 0) {
  stop("STUDY1_HETEROGENEITY_ETA must be one finite nonnegative number.")
}
if (!exists("STUDY1_THRESHOLD_DESIGN")) STUDY1_THRESHOLD_DESIGN <- "balanced"
STUDY1_THRESHOLD_DESIGN <- match.arg(
  STUDY1_THRESHOLD_DESIGN, c("balanced", "imbalanced")
)

FIG_DIR <- file.path(STUDY1_OUT_PREFIX, "figures")
TAB_DIR <- file.path(STUDY1_OUT_PREFIX, "tables")
RES_DIR <- file.path(STUDY1_OUT_PREFIX, "results", "study1_publication")
REP_DIR <- file.path(RES_DIR, "competitor_replications")
for (dd in c(FIG_DIR, TAB_DIR, RES_DIR, REP_DIR)) {
  dir.create(dd, showWarnings = FALSE, recursive = TRUE)
}

set.seed(20260705)
n_train <- study1_setting_override("STUDY1_MC_N_TRAIN", 100L)
n_test <- settings$n_test
n_rep <- settings$n_rep
m <- study1_setting_override("STUDY1_MC_M", 6L, minimum = 2L)
calib_grid <- as.integer(get0(
  "STUDY1_CALIB_GRID", inherits = TRUE,
  ifnotfound = c(0L, 20L, 50L)
))
if (length(calib_grid) < 1L || anyNA(calib_grid) ||
    any(calib_grid < 0L | calib_grid > n_train) || anyDuplicated(calib_grid)) {
  stop("STUDY1_CALIB_GRID must contain unique sizes between zero and n_train.")
}
calib_grid <- sort(calib_grid)
n_pred_draw <- settings$n_pred_draw
if (!exists("STUDY1_PARALLEL_LEVEL")) STUDY1_PARALLEL_LEVEL <- "chains"
if (!STUDY1_PARALLEL_LEVEL %in% c("chains", "replications", "none", "hybrid")) {
  stop("STUDY1_PARALLEL_LEVEL must be chains, replications, or none.")
}
parallel_chains <- STUDY1_PARALLEL_LEVEL %in% c("chains", "hybrid")
chain_cores <- if (STUDY1_PARALLEL_LEVEL == "hybrid") {
  mixedgp_as_integer_strict(STUDY1_CHAIN_WORKERS, "chain workers", 1L, 1L)
} else NULL
replication_cores <- if (STUDY1_PARALLEL_LEVEL == "hybrid") {
  min(n_rep, mixedgp_as_integer_strict(STUDY1_DATASET_WORKERS, "dataset workers", 1L, 1L))
} else if (identical(STUDY1_PARALLEL_LEVEL, "replications")) {
  min(n_rep, mixedgp_resolve_cores())
} else {
  1L
}

STUDY1_DESIGN_TAG <- paste0(
  "study1_", STUDY1_SCENARIO,
  "_eta", gsub("[^0-9]+", "p", format(STUDY1_HETEROGENEITY_ETA)),
  "_", STUDY1_THRESHOLD_DESIGN,
  "_n", n_train, "_m", m,
  "_audited_sampler_v7_nlpd_all_tasks_",
  STUDY1_CONFIG,
  "_meval", settings$n_m_eval,
  "_mdraw", settings$n_m_draw,
  "_mint", settings$n_m_latent,
  "_", paste(STUDY1_PUBLISHED_COMPETITORS, collapse = "-")
)

method_levels <- c(
  "Oracle", "EIV-GP", "UC-GP", "LVGP", "EzGP"
)
method_cols <- c(
  "Oracle" = "black",
  "EIV-GP" = "firebrick",
  "UC-GP" = "steelblue4",
  "LVGP" = "darkorange3",
  "EzGP" = "purple4"
)
metric_levels <- c(
  "RMSE", "MAE", "CRPS", "NLPD", "Coverage95", "Width95",
  "IntervalScore95"
)

safe_se <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

format_mean_se <- function(mean, se, digits = 3L) {
  if (!is.finite(mean)) return("--")
  if (!is.finite(se)) return(sprintf(paste0("%.", digits, "f"), mean))
  paste0(
    sprintf(paste0("%.", digits, "f"), mean),
    " (", sprintf(paste0("%.", digits, "f"), se), ")"
  )
}

study1_write_csv_safe <- function(x, path) {
  if (!is.data.frame(x)) stop("Study I CSV output must be a data frame.")
  if (ncol(x) == 0L) x <- data.frame(note = character(0))
  write.csv(x, path, row.names = FALSE)
  invisible(x)
}

study1_sampler_control_rows <- function(fit, rep_id, n_calib) {
  encode <- function(value) {
    if (is.numeric(value)) {
      return(paste(format(value, digits = 17L, scientific = FALSE), collapse = ";"))
    }
    paste(as.character(value), collapse = ";")
  }
  settings <- c(
    fit$control,
    list(
      initialization_rule = paste(
        "threshold-compatible latent initialization with independent",
        "seeded perturbations; see fit_eivgp_1d()"
      ),
      chain_seeds = fit$mcmc$chain_stats$seed,
      covariance_jitter = fit$diagnostics$summary$covariance_jitter,
      forms_explicit_covariance_inverse =
        fit$diagnostics$summary$forms_explicit_covariance_inverse
    )
  )
  data.frame(
    rep = as.integer(rep_id),
    n_calib = as.integer(n_calib),
    setting = names(settings),
    value = vapply(settings, encode, character(1L)),
    stringsAsFactors = FALSE
  )
}

competitor_controls <- NULL

measurement_n_iter <- if (STUDY1_QUICK) 500L else if (
  STUDY1_CONFIG == "balanced"
) 1500L else 3000L
measurement_burn <- if (STUDY1_QUICK) 150L else if (
  STUDY1_CONFIG == "balanced"
) 500L else 1000L
measurement_n_iter <- study1_setting_override(
  "STUDY1_MEAS_N_ITER", measurement_n_iter, minimum = 2L
)
measurement_burn <- study1_setting_override(
  "STUDY1_MEAS_BURN", measurement_burn, minimum = 0L
)
if (measurement_burn >= measurement_n_iter) {
  stop("STUDY1_MEAS_BURN must be smaller than STUDY1_MEAS_N_ITER.")
}
STUDY1_DESIGN_TAG <- paste0(
  STUDY1_DESIGN_TAG,
  "_nrep", n_rep,
  "_ntest", n_test,
  "_iter", settings$n_iter,
  "_burn", settings$burn,
  "_chains", settings$n_chains,
  "_pred", n_pred_draw,
  "_cal", paste(calib_grid, collapse = "-"),
  "_meas", measurement_n_iter, "b", measurement_burn,
  "_f", as.integer(isTRUE(STUDY1_EVALUATE_F)),
  "_u", as.integer(isTRUE(STUDY1_EVALUATE_U))
)
STUDY1_AUTO_CONTINUE <- get0("STUDY1_AUTO_CONTINUE", inherits=TRUE, ifnotfound=!STUDY1_QUICK)
STUDY1_TARGET_ESS <- get0("STUDY1_TARGET_ESS", inherits=TRUE, ifnotfound=100)
STUDY1_CONTINUE_CHUNK <- get0("STUDY1_CONTINUE_CHUNK", inherits=TRUE, ifnotfound=5000L)
if (isTRUE(STUDY1_AUTO_CONTINUE) && settings$n_chains<4L)
  stop("Automatic scientific-target ESS continuation requires four chains.")
STUDY1_CACHE_SPEC <- list(
  schema = "s1v10_scientific_ess_continuation",
  continuation=list(enabled=STUDY1_AUTO_CONTINUE,ess_min=STUDY1_TARGET_ESS,chunk_iter=STUDY1_CONTINUE_CHUNK),
  design_label = STUDY1_DESIGN_TAG,
  study1_config = STUDY1_CONFIG,
  scenario = STUDY1_SCENARIO,
  heterogeneity_eta = STUDY1_HETEROGENEITY_ETA,
  threshold_design = STUDY1_THRESHOLD_DESIGN,
  n_train = n_train,
  n_test = n_test,
  n_rep = n_rep,
  m = m,
  calibration_grid = calib_grid,
  mcmc = list(
    n_iter = settings$n_iter, burn = settings$burn,
    n_chains = settings$n_chains, preset = settings$preset,
    rhat_limit = STUDY1_MAX_RHAT, ess_limit = STUDY1_MIN_ESS
  ),
  diagnostic_panel = list(
    max_draws_per_chain = STUDY1_DIAGNOSTIC_MAX_DRAWS,
    n_latent = STUDY1_DIAGNOSTIC_N_LATENT,
    n_points = STUDY1_DIAGNOSTIC_N_POINTS
  ),
  prediction = list(n_draw = n_pred_draw),
  mean_recovery = list(
    n_eval = settings$n_m_eval,
    n_draw = settings$n_m_draw,
    n_latent = settings$n_m_latent
  ),
  measurement = list(n_iter = measurement_n_iter, burn = measurement_burn),
  task_flags = list(
    run_ablations = isTRUE(STUDY1_RUN_ABLATIONS),
    evaluate_f = isTRUE(STUDY1_EVALUATE_F),
    evaluate_u = isTRUE(STUDY1_EVALUATE_U)
  ),
  published_competitors = STUDY1_PUBLISHED_COMPETITORS,
  competitor_controls = competitor_controls
)
STUDY1_CACHE_SPEC$fingerprint <-
  mixedgp_object_fingerprint(STUDY1_CACHE_SPEC)
STUDY1_DESIGN_TAG <- paste0(
  "s1v7-", gsub("[^0-9A-Za-z]+", "-", STUDY1_SCENARIO),
  "-eta", gsub("[^0-9A-Za-z]+", "p", STUDY1_HETEROGENEITY_ETA),
  "-", substr(STUDY1_CACHE_SPEC$fingerprint, 1L, 16L)
)
saveRDS(
  STUDY1_CACHE_SPEC,
  file.path(RES_DIR, paste0("cache_spec_", STUDY1_DESIGN_TAG, ".rds")),
  version = 3L
)
writeLines(
  capture.output(dput(STUDY1_CACHE_SPEC)),
  file.path(RES_DIR, paste0("cache_spec_", STUDY1_DESIGN_TAG, ".R"))
)

preflight <- mixedgp_competitor_preflight(
  STUDY1_PUBLISHED_COMPETITORS,
  strict = STUDY1_STRICT_COMPETITORS
)
write.csv(
  preflight,
  file.path(TAB_DIR, "study1_competitor_preflight.csv"),
  row.names = FALSE
)

}


# ---- simulations/study1/02_study1_monte_carlo.R ----
eivgp_study1_driver <- function() {
## Fitting/evaluation only. Reporting lives in report_study1_results.R.
eval(body(eivgp_study1_setup), envir = environment())

run_one_study1_replication <- function(rep_id, run_eiv = TRUE) {
  data_path <- file.path(
    STUDY1_DATA_DIR,
    mixedgp_dataset_filename(
      "study1", rep_id, STUDY1_SCENARIO, n_train, n_test, m
    )
  )
  if (!file.exists(data_path)) {
    stop(
      "Missing frozen Study I dataset: ", data_path,
      ". Run 12_generate_synthetic_datasets.R first."
    )
  }
  frozen <- load_mixedgp_synthetic_dataset_strict(
    data_path,
    expected = list(
      study = "study1", scenario = STUDY1_SCENARIO, rep_id = rep_id,
      n = n_train, n_test = n_test, m = m,
      calib_grid = calib_grid,
      threshold_design = STUDY1_THRESHOLD_DESIGN,
      heterogeneity_eta = STUDY1_HETEROGENEITY_ETA
    )
  )
  dat <- frozen$data
  train <- dat$train
  test <- dat$test
  calib_sets <- frozen$calibration_sets

  oracle_draws <- sample_oracle_test_y(
    x_test = test$x,
    c_test = test$c,
    tau_true = dat$tau_true,
    scenario = STUDY1_SCENARIO,
    sigma_eps = dat$sigma_eps,
    n_draw = n_pred_draw,
    heterogeneity_eta = STUDY1_HETEROGENEITY_ETA
  )
  metrics <- list(
    Oracle = summarize_predictive_samples_1d(
      oracle_draws, test$y, method = "Oracle", rep_id = rep_id,
      n_calib = NA_integer_, scenario = STUDY1_SCENARIO
    )
  )
  set.seed(240000L + rep_id)
  mean_eval_idx <- sort(sample(
    seq_len(nrow(test)), min(settings$n_m_eval, nrow(test))
  ))
  m_true <- m0_1d(
    test$x[mean_eval_idx], test$c[mean_eval_idx],
    tau = dat$tau_true, scenario = STUDY1_SCENARIO
  )
  mean_metrics <- list()
  ablation_metrics <- list()
  latent_imputation_metrics <- list()
  surface_recovery_metrics <- list()
  ablation_surface_metrics <- list()
  ablation_status <- list()
  measurement_fits <- list()
  mcmc_diagnostics <- list()
  mcmc_parameter_diagnostics <- list()
  sampler_controls <- list()
  sampler_control_manifest <- list()


  if (isTRUE(STUDY1_RUN_ABLATIONS)) {
    if (isTRUE(STUDY1_EVALUATE_F)) {
      start_full_u <- proc.time()[3]
      full_u_fit <- tryCatch(
        fit_study1_full_u_gp(train$x, train$y, train$u),
        error = function(e) e
      )
      elapsed_full_u <- proc.time()[3] - start_full_u
      if (inherits(full_u_fit, "error")) {
        ablation_status[["Full-U GP"]] <- data.frame(
          method = "Full-U GP", n_calib = NA_integer_, status = "failed",
          message = conditionMessage(full_u_fit),
          elapsed_seconds = elapsed_full_u
        )
      } else {
        ablation_status[["Full-U GP"]] <- data.frame(
          method = "Full-U GP", n_calib = NA_integer_, status = "success",
          message = "", elapsed_seconds = elapsed_full_u
        )
        ablation_surface_metrics[["Full-U GP"]] <-
          study1_latent_gp_surface_metrics(
            gp_fit = full_u_fit,
            scenario = STUDY1_SCENARIO,
            method = "Full-U GP",
            rep_id = rep_id,
            n_calib = NA_integer_,
            tau_true = dat$tau_true,
            heterogeneity_eta = STUDY1_HETEROGENEITY_ETA,
            grid_n_x = if (STUDY1_QUICK) 15L else 31L,
            grid_n_u = if (STUDY1_QUICK) 21L else 51L
          )
      }
    }
    for (n_calib in calib_grid) {
      key <- as.character(n_calib)
      start_measurement <- proc.time()[3]
      measurement_fit <- tryCatch(
        fit_threshold_measurement_response_free(
          c_ord = train$c,
          u_obs = train$u,
          calib_idx = calib_sets[[key]],
          m = m,
          n_iter = measurement_n_iter,
          burn = measurement_burn,
          thin = 1L,
          seed = 270000L + 1000L * rep_id + n_calib
        ),
        error = function(e) e
      )
      elapsed_measurement <- proc.time()[3] - start_measurement
      if (inherits(measurement_fit, "error")) {
        msg <- conditionMessage(measurement_fit)
        ablation_status[[paste0("measurement_", key)]] <- data.frame(
          method = "Response-free threshold model", n_calib = n_calib,
          status = "failed", message = msg,
          elapsed_seconds = elapsed_measurement
        )
        ablation_status[[paste0("PI_", key)]] <- data.frame(
          method = "PI-GP", n_calib = n_calib,
          status = "failed", message = msg, elapsed_seconds = NA_real_
        )
        ablation_status[[paste0("CC_", key)]] <- data.frame(
          method = "CC-GP", n_calib = n_calib,
          status = if (n_calib == 0L) "not_applicable" else "failed",
          message = if (n_calib == 0L) "No complete cases." else msg,
          elapsed_seconds = NA_real_
        )
        next
      }
      ablation_status[[paste0("measurement_", key)]] <- data.frame(
        method = "Response-free threshold model", n_calib = n_calib,
        status = "success", message = "",
        elapsed_seconds = elapsed_measurement
      )
      measurement_fits[[key]] <- measurement_fit

      start_pi <- proc.time()[3]
      pi_result <- tryCatch({
        fit <- fit_study1_pi_gp(train$x, train$y, measurement_fit)
        draws <- sample_study1_pi_gp(
          fit, test$x, test$c, n_pred_draw,
          seed = 280000L + 1000L * rep_id + n_calib
        )
        list(fit = fit, draws = draws)
      }, error = function(e) e)
      elapsed_pi <- proc.time()[3] - start_pi
      if (inherits(pi_result, "error")) {
        ablation_status[[paste0("PI_", key)]] <- data.frame(
          method = "PI-GP", n_calib = n_calib,
          status = "failed", message = conditionMessage(pi_result),
          elapsed_seconds = elapsed_pi
        )
      } else {
        ablation_status[[paste0("PI_", key)]] <- data.frame(
          method = "PI-GP", n_calib = n_calib,
          status = "success", message = "", elapsed_seconds = elapsed_pi
        )
        ablation_metrics[[paste0("PI_", key)]] <- summarize_predictive_samples_1d(
          pi_result$draws, test$y, method = "PI-GP", rep_id = rep_id,
          n_calib = n_calib, scenario = STUDY1_SCENARIO
        )
        if (isTRUE(STUDY1_EVALUATE_F) && n_calib > 0L) {
          ablation_surface_metrics[[paste0("PI_", key)]] <-
            study1_latent_gp_surface_metrics(
              gp_fit = pi_result$fit,
              scenario = STUDY1_SCENARIO,
              method = "PI-GP",
              rep_id = rep_id,
              n_calib = n_calib,
              tau_true = dat$tau_true,
              heterogeneity_eta = STUDY1_HETEROGENEITY_ETA,
              grid_n_x = if (STUDY1_QUICK) 15L else 31L,
              grid_n_u = if (STUDY1_QUICK) 21L else 51L
            )
        }
      }

      if (n_calib == 0L) {
        ablation_status[[paste0("CC_", key)]] <- data.frame(
          method = "CC-GP", n_calib = n_calib,
          status = "not_applicable", message = "No complete cases.",
          elapsed_seconds = 0
        )
      } else {
        start_cc <- proc.time()[3]
        cc_result <- tryCatch({
          fit <- fit_study1_cc_gp(
            train$x, train$y, train$u, calib_sets[[key]], measurement_fit
          )
          draws <- sample_study1_cc_gp(
            fit, test$x, test$c, n_pred_draw,
            seed = 290000L + 1000L * rep_id + n_calib
          )
          list(fit = fit, draws = draws)
        }, error = function(e) e)
        elapsed_cc <- proc.time()[3] - start_cc
        if (inherits(cc_result, "error")) {
          ablation_status[[paste0("CC_", key)]] <- data.frame(
            method = "CC-GP", n_calib = n_calib,
            status = "failed", message = conditionMessage(cc_result),
            elapsed_seconds = elapsed_cc
          )
        } else {
          ablation_status[[paste0("CC_", key)]] <- data.frame(
            method = "CC-GP", n_calib = n_calib,
            status = "success", message = "", elapsed_seconds = elapsed_cc
          )
          ablation_metrics[[paste0("CC_", key)]] <- summarize_predictive_samples_1d(
            cc_result$draws, test$y, method = "CC-GP", rep_id = rep_id,
            n_calib = n_calib, scenario = STUDY1_SCENARIO
          )
          if (isTRUE(STUDY1_EVALUATE_F)) {
            ablation_surface_metrics[[paste0("CC_", key)]] <-
              study1_latent_gp_surface_metrics(
                gp_fit = cc_result$fit,
                scenario = STUDY1_SCENARIO,
                method = "CC-GP",
                rep_id = rep_id,
                n_calib = n_calib,
                tau_true = dat$tau_true,
                heterogeneity_eta = STUDY1_HETEROGENEITY_ETA,
                grid_n_x = if (STUDY1_QUICK) 15L else 31L,
                grid_n_u = if (STUDY1_QUICK) 21L else 51L
              )
          }
        }
      }
    }
  }

  if (isTRUE(run_eiv)) {
    for (n_calib in calib_grid) {
      message("Study I replication ", rep_id, ": EIV-GP |O|=", n_calib)
      continuation_file <- file.path(RES_DIR,sprintf("continuation_%s_rep%03d_cal%03d.rds",STUDY1_DESIGN_TAG,rep_id,n_calib))
      final_file <- file.path(RES_DIR,sprintf("fit_rep%03d_cal%03d.rds",rep_id,n_calib))
      final_saved <- if(file.exists(final_file))read_fit(final_file)else NULL
      fit_eiv <- if(file.exists(continuation_file))read_fit(continuation_file)$fit else
        if(!is.null(final_saved) && identical(final_saved$cache_spec$fingerprint,STUDY1_CACHE_SPEC$fingerprint))final_saved$fit else fit_eivgp_1d(
        x_raw = train$x,
        y_raw = train$y,
        c_ord = train$c,
        u_true = train$u,
        calib_idx = calib_sets[[as.character(n_calib)]],
        m = m,
        tau_true = dat$tau_true,
        n_iter = settings$n_iter,
        burn = settings$burn,
        thin = 1L,
        n_chains = settings$n_chains,
        preset = settings$preset,
        seed = 300000L + 1000L * rep_id + n_calib,
        parallel_chains = parallel_chains,
        n_cores = chain_cores,
        verbose = FALSE
      )
      observed_u<-rep(NA_real_,nrow(train));observed_u[calib_sets[[as.character(n_calib)]]]<-train$u[calib_sets[[as.character(n_calib)]]]
      if(!identical(as.numeric(fit_eiv$data$x_raw),as.numeric(train$x)) ||
         !identical(as.numeric(fit_eiv$data$y_raw),as.numeric(train$y)) ||
         !identical(as.integer(fit_eiv$data$c_ord),as.integer(train$c)) ||
         !identical(as.numeric(fit_eiv$data$u_obs),observed_u))
        stop("Saved Study I fit does not match this frozen dataset and calibration set.")
      control_key <- as.character(n_calib)
      sampler_controls[[control_key]] <- list(
        control = fit_eiv$control,
        initialization_rule = paste(
          "threshold-compatible latent initialization with independent",
          "seeded perturbations; see fit_eivgp_1d()"
        ),
        chain_seeds = fit_eiv$mcmc$chain_stats$seed,
        covariance_jitter = fit_eiv$diagnostics$summary$covariance_jitter,
        forms_explicit_covariance_inverse =
          fit_eiv$diagnostics$summary$forms_explicit_covariance_inverse
      )
      sampler_control_manifest[[control_key]] <-
        study1_sampler_control_rows(fit_eiv, rep_id, n_calib)
      panel <- mixedgp_study2_diagnostic_rows(
        matrix(test$x, ncol = 1L), matrix(test$c, ncol = 1L),
        max_points = STUDY1_DIAGNOSTIC_N_POINTS
      )
      if(isTRUE(STUDY1_AUTO_CONTINUE)) {
        if(file.exists(continuation_file)) {
          pending<-read_fit(continuation_file)
          if(!identical(pending$fit$data,fit_eiv$data) || !identical(pending$fit$control,fit_eiv$control))
            stop("Saved continuation is incompatible with this Study I fit.")
          fit_eiv<-pending$fit
        }
        extended <- continue_eivgp_until_ess(fit_eiv,
          X=matrix(test$x[panel],ncol=1L),C=test$c[panel],
          U=if(isTRUE(STUDY1_EVALUATE_F) && STUDY1_HETEROGENEITY_ETA==1)test$u[panel] else NULL,
          ess_min=STUDY1_TARGET_ESS,chunk_iter=STUDY1_CONTINUE_CHUNK,
          n_latent=STUDY1_DIAGNOSTIC_N_LATENT,seed=380000L+1000L*rep_id+n_calib,
          parallel=parallel_chains,n_cores=chain_cores,checkpoint_file=continuation_file)
        fit_eiv<-extended$fit;target_series<-extended$target_series;target_diag<-extended$diagnostics
      } else {
      target_series <- mixedgp_study1_target_series(
        fit_eiv, X = matrix(test$x[panel], ncol = 1L), C = test$c[panel],
        U = if (isTRUE(STUDY1_EVALUATE_F) && STUDY1_HETEROGENEITY_ETA == 1) {
          test$u[panel]
        } else NULL,
        max_draws_per_chain = STUDY1_DIAGNOSTIC_MAX_DRAWS,
        n_latent = STUDY1_DIAGNOSTIC_N_LATENT,
        seed = 380000L + 1000L * rep_id + n_calib
      )
      target_diag <- mixedgp_summarize_diagnostic_series(
        target_series, rhat_limit = STUDY1_MAX_RHAT,
        bulk_ess_limit = STUDY1_MIN_ESS, tail_ess_limit = STUDY1_MIN_ESS
      )
      }
      diag_row <- fit_eiv$diagnostics$summary
      parameter_diag <- mixedgp_study1_raw_diagnostics(
        fit_eiv, rhat_limit = STUDY1_MAX_RHAT,
        bulk_ess_limit = STUDY1_MIN_ESS, tail_ess_limit = STUDY1_MIN_ESS
      )
      raw_gate_pass <- mixedgp_diagnostic_table_pass(
        parameter_diag, rhat_limit = STUDY1_MAX_RHAT,
        bulk_ess_limit = STUDY1_MIN_ESS, tail_ess_limit = STUDY1_MIN_ESS,
        expected_parameters = names(mixedgp_study1_raw_series(fit_eiv))
      )
      target_gate_pass <- mixedgp_diagnostic_table_pass(
        target_diag, rhat_limit = STUDY1_MAX_RHAT,
        bulk_ess_limit = STUDY1_MIN_ESS, tail_ess_limit = STUDY1_MIN_ESS,
        expected_parameters = names(target_series)
      )
      # Retain small scientific-target traces for every fit, including good fits.
      # Full-window diagnostics above are never recomputed on display-thinned traces.
      trace_n <- length(target_series[[1L]][[1L]])
      trace_keep <- unique(c(seq.int(1L, trace_n, by=max(1L, ceiling(trace_n/750L))), trace_n))
      trace_series <- lapply(target_series, function(chains) lapply(chains, function(v) v[trace_keep]))
      trace_series <- mixedgp_mark_diagnostic_window(trace_series, length(trace_keep)==trace_n)
      mixedgp_save_artifact(list(schema_version=1L, series=trace_series,
        iteration=settings$burn+trace_keep, retained_index=trace_keep,
        panel_rows=panel, test_panel=test[panel,], diagnostics=target_diag,
        source="retained target functionals from this fit"),
        file.path(RES_DIR, sprintf("target_traces_rep%03d_cal%03d.rds",rep_id,n_calib)))
      parameter_diag$diagnostic_role <- "raw"
      target_diag$diagnostic_role <- "target"
      parameter_diag <- bind_rows(parameter_diag, target_diag)
      scientific_ess_pass <- mixedgp_scientific_ess_pass(target_diag,STUDY1_TARGET_ESS)
      gate_pass <- scientific_ess_pass
      max_rhat <- if (nrow(parameter_diag)) max(parameter_diag$rhat) else NA_real_
      min_ess <- if (nrow(parameter_diag)) {
        min(parameter_diag$ess_bulk, parameter_diag$ess_tail)
      } else NA_real_
      parameter_diag$rep <- rep_id
      parameter_diag$n_calib <- n_calib
      mcmc_parameter_diagnostics[[as.character(n_calib)]] <- parameter_diag
      diag_row$rep <- rep_id
      diag_row$n_calib <- n_calib
      diag_row$gate_max_rhat <- max_rhat
      diag_row$gate_min_ess <- min_ess
      diag_row$gate_pass <- gate_pass
      diag_row$raw_gate_pass <- raw_gate_pass
      diag_row$target_gate_pass <- target_gate_pass
      diag_row$scientific_ess_pass <- scientific_ess_pass
      mcmc_diagnostics[[as.character(n_calib)]] <- diag_row
      diag_row$diagnostic_warning <- !isTRUE(gate_pass)
      diag_row$diagnostic_advice <- if (isTRUE(gate_pass)) "" else mixedgp_simulation_diagnostic_advice()
      mcmc_diagnostics[[as.character(n_calib)]] <- diag_row
      # Keep a final fit with terminal states for every completed inference.
      mixedgp_checkpoint_write(mixedgp_pack_artifact(list(fit=fit_eiv,diagnostics=parameter_diag,
        diagnostic_panel_rows=panel,cache_spec=STUDY1_CACHE_SPEC)),
        file.path(RES_DIR,sprintf("fit_rep%03d_cal%03d.rds",rep_id,n_calib)))
      if(file.exists(continuation_file))unlink(continuation_file)
      if (!isTRUE(gate_pass)) {
        failure_file <- file.path(RES_DIR,
          sprintf("mcmc_failure_rep%03d_cal%03d.rds", rep_id, n_calib))
        mixedgp_save_artifact(list(fit = fit_eiv, diagnostics = parameter_diag,
                     diagnostic_panel_rows = panel, cache_spec = STUDY1_CACHE_SPEC),
                failure_file)
        warning(
          "Study I MCMC diagnostics flagged replication ", rep_id,
          ", |O|=", n_calib, ": max R-hat=", signif(max_rhat, 4),
          ", min ESS=", signif(min_ess, 5),
          ". Fit retained and diagnostics saved to ", failure_file, ". ",
          mixedgp_simulation_diagnostic_advice(), call. = FALSE
        )
      }
      scored <- mixedgp_score_study1_fit(fit_eiv, frozen, STUDY1_CACHE_SPEC, rep_id, n_calib,
        measurement_fits[[as.character(n_calib)]])
      metrics[[paste0("EIV_",n_calib)]] <- scored$metrics
      mean_metrics[[paste0("EIV_",n_calib)]] <- scored$mean_recovery
      latent_imputation_metrics[[as.character(n_calib)]] <- scored$latent_imputation
      surface_recovery_metrics[[as.character(n_calib)]] <- scored$surface_recovery

    }
  }

  status <- data.frame(method=character(),status=character(),rep=integer(),scenario=character())
  ablation_status_df <- bind_rows(ablation_status)
  if (nrow(ablation_status_df) > 0L) {
    ablation_status_df$rep <- rep_id
    ablation_status_df$scenario <- STUDY1_SCENARIO
  }
  list(
    dataset_identity = data.frame(rep=rep_id, dataset_md5=mixedgp_dataset_identity(frozen$data)),
    metrics = bind_rows(metrics),
    competitor_status = status,
    ablation_metrics = bind_rows(ablation_metrics),
    ablation_status = ablation_status_df,
    mcmc_diagnostics = bind_rows(mcmc_diagnostics),
    mcmc_parameter_diagnostics = bind_rows(mcmc_parameter_diagnostics),
    sampler_control_manifest = bind_rows(sampler_control_manifest),
    mean_recovery = bind_rows(mean_metrics),
    latent_imputation = bind_rows(latent_imputation_metrics),
    surface_recovery = bind_rows(surface_recovery_metrics),
    ablation_surface_recovery = bind_rows(ablation_surface_metrics),
    metadata = list(
      rep = rep_id,
      design_tag = STUDY1_DESIGN_TAG,
      data_seed = 100000L + rep_id,
      calibration_seed = 200000L + rep_id,
      frozen_data_file = basename(data_path),
      frozen_data_md5 = attr(frozen, "manifest_md5"),
      frozen_manifest = attr(frozen, "manifest_path"),
      sampler = sampler_controls
    )
  )
}

run_or_load_study1_replication <- function(rr) {
  cache_mode <- if (isTRUE(STUDY1_REUSE_LOCKED_EIV)) "archival" else "fresh"
  rep_file <- file.path(
    REP_DIR,
    sprintf("study1_rep%03d_%s_%s.rds", rr, STUDY1_DESIGN_TAG, cache_mode)
  )
  if (isTRUE(STUDY1_USE_CACHE) && file.exists(rep_file)) {
    cached <- tryCatch(readRDS(rep_file), error = function(e) NULL)
    data_path <- file.path(
      STUDY1_DATA_DIR,
      mixedgp_dataset_filename(
        "study1", rr, STUDY1_SCENARIO, n_train, n_test, m
      )
    )
    cache_valid <- is.list(cached$metadata) &&
      identical(cached$metadata$rep, as.integer(rr)) &&
      identical(cached$metadata$design_tag, STUDY1_DESIGN_TAG) &&
      identical(cached$metadata$frozen_data_file, basename(data_path)) &&
      file.exists(data_path) &&
      identical(
        tolower(as.character(cached$metadata$frozen_data_md5)),
        tolower(unname(tools::md5sum(data_path)))
      )
    if(!is.null(STUDY2_ESS_MIN))cache_valid<-isTRUE(cache_valid) &&
      identical(cached$metadata$scientific_ess_min,STUDY2_ESS_MIN) &&
      is.data.frame(cached$diagnostics) && nrow(cached$diagnostics)==length(scenario_calib_grid(scenario)) && "scientific_ess_pass" %in% names(cached$diagnostics) &&
      all(cached$diagnostics$scientific_ess_pass %in% TRUE)
    if (isTRUE(cache_valid)) { attr(cached,"eivgp_cache_hit") <- TRUE; return(cached) }
    message("Ignoring incompatible Study I cache: ", rep_file)
  }
  message("Study I replication ", rr, " of ", n_rep)
  out <- run_one_study1_replication(
    rr,
    run_eiv = !isTRUE(STUDY1_REUSE_LOCKED_EIV)
  )
  mixedgp_save_replication(out, rep_file)
  out
}
rep_objects <- mixedgp_run_replications(
  as.list(seq_len(n_rep)),
  run_or_load_study1_replication,
  n_cores = replication_cores,
  seeds = 910000L + seq_len(n_rep),
  status_path = file.path(TAB_DIR, "study1_replication_status.csv"),
  study = "Study I",
  parallel_map = mixedgp_parallel_lapply,
  mc.preschedule = FALSE
)

dataset_identity <- bind_rows(lapply(rep_objects, `[[`, "dataset_identity"))
write.csv(dataset_identity, file.path(STUDY1_OUT_PREFIX, "dataset_identity.csv"), row.names=FALSE)

new_results <- bind_rows(lapply(rep_objects, `[[`, "metrics"))
competitor_status <- bind_rows(
  lapply(rep_objects, `[[`, "competitor_status")
)
ablation_results <- bind_rows(lapply(rep_objects, `[[`, "ablation_metrics"))
ablation_status <- bind_rows(lapply(rep_objects, `[[`, "ablation_status"))
mcmc_diagnostics <- bind_rows(lapply(rep_objects, `[[`, "mcmc_diagnostics"))
mcmc_parameter_diagnostics <- bind_rows(
  lapply(rep_objects, `[[`, "mcmc_parameter_diagnostics")
)
sampler_control_manifest <- bind_rows(
  lapply(rep_objects, `[[`, "sampler_control_manifest")
)
mean_recovery <- bind_rows(lapply(rep_objects, `[[`, "mean_recovery"))
latent_imputation <- bind_rows(lapply(rep_objects, `[[`, "latent_imputation"))
surface_recovery <- bind_rows(lapply(rep_objects, `[[`, "surface_recovery"))
ablation_surface_recovery <- bind_rows(
  lapply(rep_objects, `[[`, "ablation_surface_recovery")
)
mixedgp_write_diagnostic_tables(mcmc_diagnostics, mcmc_parameter_diagnostics,
  TAB_DIR, "study1", overwrite = TRUE)
if (nrow(sampler_control_manifest) > 0L) {
  write.csv(
    sampler_control_manifest,
    file.path(TAB_DIR, "study1_sampler_control_manifest.csv"),
    row.names = FALSE
  )
}

if (isTRUE(STUDY1_REUSE_LOCKED_EIV)) {
  locked_file <- file.path(TAB_DIR, "study1_mc_raw_results.csv")
  if (!file.exists(locked_file)) {
    stop("Locked July 27 Study I results are missing: ", locked_file)
  }
  locked <- read.csv(locked_file, stringsAsFactors = FALSE)
  locked_eiv <- locked |>
    filter(method == "EIV-GP", rep %in% seq_len(n_rep))
  expected_rows <- n_rep * length(calib_grid)
  if (nrow(locked_eiv) != expected_rows ||
      !setequal(unique(locked_eiv$n_calib), calib_grid)) {
    stop("Locked EIV-GP rows do not match the archived 50-by-5 Study I design.")
  }
  ## Oracle values in the archived file are duplicated over calibration sizes.
  ## Keep one copy per replication for the revised table and curves.
  locked_oracle <- locked |>
    filter(method == "Oracle", rep %in% seq_len(n_rep)) |>
    group_by(rep, scenario, method) |>
    slice(1L) |>
    ungroup() |>
    mutate(n_calib = NA_integer_)
  mc_results <- bind_rows(
    locked_eiv,
    locked_oracle,
    new_results |> filter(!method %in% c("EIV-GP", "Oracle"))
  )
} else {
  mc_results <- new_results
}


## Persist reporting inputs before any table or figure generation.
report_names <- c(
  "ablation_results", "ablation_status", "ablation_surface_recovery", "calib_grid",
  "competitor_status", "latent_imputation", "m", "mc_results",
  "mcmc_diagnostics", "mcmc_parameter_diagnostics", "mean_recovery", "method_cols",
  "method_levels", "metric_levels", "n_rep", "n_test",
  "n_train", "sampler_control_manifest", "settings", "STUDY1_DESIGN_TAG",
  "STUDY1_HETEROGENEITY_ETA", "STUDY1_MAX_RHAT", "STUDY1_MECHANISM_CALIB", "STUDY1_MIN_ESS",
  "STUDY1_PUBLISHED_COMPETITORS", "STUDY1_REQUIRE_MCMC_GATE", "STUDY1_REUSE_LOCKED_EIV", "STUDY1_RUN_ABLATIONS",
  "STUDY1_SCENARIO", "STUDY1_THRESHOLD_DESIGN", "surface_recovery")
report_names <- intersect(report_names, ls(envir = environment(), all.names = TRUE))
report_context <- mget(report_names, envir = environment(), inherits = FALSE)
report_context <- Filter(function(x) !is.function(x) && !is.environment(x), report_context)
mixedgp_atomic_save_rds(list(schema = 1L, study = "study1", context = report_context),
        file.path(STUDY1_OUT_PREFIX, "report_inputs.rds"))

}


# ---- simulations/study2/setup_study2_experiment.R ----
eivgp_study2_setup <- function() {
if (!exists("STUDY2_RUN_ABLATIONS")) STUDY2_RUN_ABLATIONS <- TRUE
needed_pkgs <- c("ggplot2", "dplyr", "tidyr", "knitr", "posterior", "TruncatedNormal")
if ((exists("STUDY2_CONFIG") && !identical(STUDY2_CONFIG, "quick")) ||
    isTRUE(STUDY2_RUN_ABLATIONS)) {
  needed_pkgs <- c(needed_pkgs, "posterior")
}
missing_pkgs <- needed_pkgs[
  !vapply(needed_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0L) {
  stop("Please install required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
})

if (!exists("STUDY2_CONFIG")) STUDY2_CONFIG <- "quick"
if (!exists("STUDY2_USE_CACHE")) STUDY2_USE_CACHE <- TRUE
if (!exists("STUDY2_MC_RESUME")) STUDY2_MC_RESUME <- TRUE
if (!exists("STUDY2_REP_IDS")) STUDY2_REP_IDS <- NULL
if (!exists("STUDY2_CHECKPOINT_FITS")) STUDY2_CHECKPOINT_FITS <- TRUE
if (!exists("STUDY2_RECOVERY_MODEL")) STUDY2_RECOVERY_MODEL <- NULL

if (!exists("STUDY2_ORACLE_MEAN_METHOD")) STUDY2_ORACLE_MEAN_METHOD <- "quadrature"
STUDY2_ORACLE_MEAN_METHOD <- match.arg(STUDY2_ORACLE_MEAN_METHOD,c("quadrature","rejection"))

if (!exists("STUDY2_OUT_PREFIX")) STUDY2_OUT_PREFIX <- ".."
if (!exists("STUDY2_SAVE_REP_FITS")) STUDY2_SAVE_REP_FITS <- FALSE
if (!exists("STUDY2_STRICT_COMPETITORS")) STUDY2_STRICT_COMPETITORS <- FALSE
if (!exists("STUDY2_ENFORCE_MCMC_GATE")) {
  STUDY2_ENFORCE_MCMC_GATE <- !identical(STUDY2_CONFIG, "quick")
}
if (!exists("STUDY2_SAVE_PDF")) STUDY2_SAVE_PDF <- TRUE
if (!exists("STUDY2_SAVE_PNG")) STUDY2_SAVE_PNG <- TRUE
if (!exists("STUDY2_ABLATION_SCENARIOS")) {
  STUDY2_ABLATION_SCENARIOS <- "primary"
}

if (!exists("STUDY2_SCENARIOS")) {
  STUDY2_SCENARIOS <- c(
    "primary",
    "latent_additive_control",
    "high_uncertainty",
    "logistic_misspec"
  )
}

STUDY2_PUBLISHED_COMPETITORS <- character() # Fitted only by the standalone runner.

if (!exists("STUDY2_PRIMARY_CALIB_GRID")) {
  STUDY2_PRIMARY_CALIB_GRID <- c(0L, 10L, 25L, 50L, 80L)
}
if (!exists("STUDY2_CONTRAST_CALIB")) STUDY2_CONTRAST_CALIB <- 25L
if (!exists("STUDY2_MC_N_REP")) STUDY2_MC_N_REP <- NULL
if (!exists("STUDY2_DATA_DIR")) {
  STUDY2_DATA_DIR <- file.path("..", "data-synthetic", "study2")
}

settings <- study2_config_settings(STUDY2_CONFIG)
study2_competitor_controls <- NULL

allowed_scenarios <- c(
  "primary",
  "latent_additive_control",
  "high_uncertainty",
  "logistic_misspec"
)
unknown_scenarios <- setdiff(STUDY2_SCENARIOS, allowed_scenarios)
if (length(STUDY2_SCENARIOS) == 0L || length(unknown_scenarios) > 0L) {
  stop(
    "STUDY2_SCENARIOS must be a nonempty subset of: ",
    paste(allowed_scenarios, collapse = ", "),
    if (length(unknown_scenarios) > 0L) {
      paste0(". Unknown: ", paste(unknown_scenarios, collapse = ", "))
    } else {
      ""
    }
  )
}
unknown_ablation_scenarios <- setdiff(
  STUDY2_ABLATION_SCENARIOS,
  allowed_scenarios
)
if (length(unknown_ablation_scenarios) > 0L) {
  stop(
    "Unknown STUDY2_ABLATION_SCENARIOS: ",
    paste(unknown_ablation_scenarios, collapse = ", ")
  )
}

############################################################
## Frozen settings and output locations
############################################################

n_train <- 120L
n_train <- as.integer(get0(
  "STUDY2_MC_N_TRAIN", inherits = TRUE, ifnotfound = n_train
))
n_test <- as.integer(get0(
  "STUDY2_MC_N_TEST", inherits = TRUE, ifnotfound = settings$n_test
))
n_rep <- if (is.null(STUDY2_MC_N_REP)) settings$n_rep else as.integer(STUDY2_MC_N_REP)
if (anyNA(c(n_train, n_test, n_rep)) ||
    any(c(n_train, n_test, n_rep) < 1L)) {
  stop("Study II n_train, n_test, and n_rep must be positive integers.")
}
study2_q <- as.integer(get0("STUDY2_Q", inherits = TRUE, ifnotfound = 4L))
study2_m <- as.integer(get0("STUDY2_M", inherits = TRUE, ifnotfound = 4L))
if (length(study2_q) != 1L || is.na(study2_q) ||
    study2_q < 2L || study2_q > 6L ||
    length(study2_m) != 1L || is.na(study2_m) || study2_m < 2L) {
  stop("STUDY2_Q must be in 2:6 and STUDY2_M must be at least two.")
}
m_vec <- rep(study2_m, study2_q)
d_latent <- 2L
ident_method <- "lower_triangular"
STUDY2_DESIGN_TAG <- paste(
  "study2-publication-v14-alltasks-rfimp",
  paste0("q", study2_q, "-d", d_latent),
  "A-fixed-nested",
  paste0("balanced-", study2_m, "-level"),
  "sigma0.12",
  "exact-minimax-interwoven",
  "common-random-numbers",
  sep = "_"
)

mc_n_iter <- if (exists("STUDY2_MC_N_ITER")) {
  as.integer(STUDY2_MC_N_ITER)
} else settings$mc_n_iter
mc_burn <- if (exists("STUDY2_MC_BURN")) {
  as.integer(STUDY2_MC_BURN)
} else settings$mc_burn
mc_thin <- if (exists("STUDY2_MC_THIN")) {
  as.integer(STUDY2_MC_THIN)
} else settings$mc_thin
mc_n_chains <- if (exists("STUDY2_MC_N_CHAINS")) {
  as.integer(STUDY2_MC_N_CHAINS)
} else settings$mc_n_chains
mc_preset <- settings$preset
measurement_n_iter <- if (exists("STUDY2_MEAS_N_ITER")) {
  as.integer(STUDY2_MEAS_N_ITER)
} else mc_n_iter
measurement_burn <- if (exists("STUDY2_MEAS_BURN")) {
  as.integer(STUDY2_MEAS_BURN)
} else mc_burn
measurement_thin <- if (exists("STUDY2_MEAS_THIN")) {
  as.integer(STUDY2_MEAS_THIN)
} else mc_thin
measurement_n_chains <- if (exists("STUDY2_MEAS_N_CHAINS")) {
  as.integer(STUDY2_MEAS_N_CHAINS)
} else max(2L, mc_n_chains)
n_pred_draw <- as.integer(get0(
  "STUDY2_MC_N_PRED_DRAW", inherits = TRUE,
  ifnotfound = settings$n_pred_draw
))
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
n_oracle_pool <- as.integer(get0(
  "STUDY2_MC_N_ORACLE_POOL", inherits = TRUE,
  ifnotfound = settings$n_oracle_pool
))
n_m_eval <- if (STUDY2_CONFIG == "quick") 30L else if (
  STUDY2_CONFIG == "balanced"
) 60L else 80L
n_m_draw <- if (STUDY2_CONFIG == "quick") 50L else if (
  STUDY2_CONFIG == "balanced"
) 120L else 200L
n_m_latent <- if (STUDY2_CONFIG == "quick") 64L else if (
  STUDY2_CONFIG == "balanced"
) 128L else 256L
n_m_truth <- if (STUDY2_CONFIG == "quick") 1000L else 2000L
n_m_eval <- as.integer(get0(
  "STUDY2_MC_N_M_EVAL", inherits = TRUE, ifnotfound = n_m_eval
))
n_m_draw <- as.integer(get0(
  "STUDY2_MC_N_M_DRAW", inherits = TRUE, ifnotfound = n_m_draw
))
n_m_latent <- as.integer(get0(
  "STUDY2_MC_N_M_LATENT", inherits = TRUE, ifnotfound = n_m_latent
))
n_m_truth <- as.integer(get0(
  "STUDY2_MC_N_M_TRUTH", inherits = TRUE, ifnotfound = n_m_truth
))
if (anyNA(c(
  n_pred_draw, n_m_eval, n_m_draw, n_m_latent, n_m_truth, n_oracle_pool
)) || any(c(
  n_pred_draw, n_m_eval, n_m_draw, n_m_latent, n_m_truth, n_oracle_pool
) < 1L)) {
  stop("Study II prediction and mean-integration sizes must be positive integers.")
}
if (!exists("STUDY2_MCMC_RHAT_LIMIT")) STUDY2_MCMC_RHAT_LIMIT <- 1.01
if (!exists("STUDY2_ESS_MIN")) STUDY2_ESS_MIN <- NULL
if (!exists("STUDY2_MCMC_RAW_ESS_LIMIT")) STUDY2_MCMC_RAW_ESS_LIMIT <- 400
if (!exists("STUDY2_MCMC_TARGET_BULK_ESS_LIMIT")) {
  STUDY2_MCMC_TARGET_BULK_ESS_LIMIT <- 400
}
if (!exists("STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT")) {
  STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT <- 400
}
if (!exists("STUDY2_MCMC_TARGET_N_POINTS")) {
  STUDY2_MCMC_TARGET_N_POINTS <- 5L
}
if (!exists("STUDY2_MCMC_TARGET_DRAWS_PER_CHAIN")) {
  STUDY2_MCMC_TARGET_DRAWS_PER_CHAIN <- Inf
}
if (!exists("STUDY2_MCMC_TARGET_N_LATENT")) {
  STUDY2_MCMC_TARGET_N_LATENT <- 64L
}
if (!exists("STUDY2_MEAS_RHAT_LIMIT")) {
  STUDY2_MEAS_RHAT_LIMIT <- STUDY2_MCMC_RHAT_LIMIT
}
if (!exists("STUDY2_MEAS_BULK_ESS_LIMIT")) {
  STUDY2_MEAS_BULK_ESS_LIMIT <- STUDY2_MCMC_TARGET_BULK_ESS_LIMIT
}
if (!exists("STUDY2_MEAS_TAIL_ESS_LIMIT")) {
  STUDY2_MEAS_TAIL_ESS_LIMIT <- STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT
}
if (!exists("STUDY2_ABLATION_GP_N_STARTS")) {
  STUDY2_ABLATION_GP_N_STARTS <- if (STUDY2_CONFIG == "quick") 2L else 8L
}
if (!exists("STUDY2_ABLATION_GP_MAXIT")) {
  STUDY2_ABLATION_GP_MAXIT <- 500L
}
if (!is.finite(STUDY2_MCMC_RHAT_LIMIT) ||
    STUDY2_MCMC_RHAT_LIMIT <= 1 ||
    !is.finite(STUDY2_MCMC_RAW_ESS_LIMIT) ||
    STUDY2_MCMC_RAW_ESS_LIMIT < 1 ||
    !is.finite(STUDY2_MCMC_TARGET_BULK_ESS_LIMIT) ||
    STUDY2_MCMC_TARGET_BULK_ESS_LIMIT < 1 ||
    !is.finite(STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT) ||
    STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT < 1 ||
    length(STUDY2_MCMC_TARGET_N_POINTS) != 1L ||
    is.na(STUDY2_MCMC_TARGET_N_POINTS) ||
    STUDY2_MCMC_TARGET_N_POINTS < 1L ||
    length(STUDY2_MCMC_TARGET_DRAWS_PER_CHAIN) != 1L ||
    is.na(STUDY2_MCMC_TARGET_DRAWS_PER_CHAIN) ||
    STUDY2_MCMC_TARGET_DRAWS_PER_CHAIN < 4L ||
    length(STUDY2_MCMC_TARGET_N_LATENT) != 1L ||
    !is.finite(STUDY2_MCMC_TARGET_N_LATENT) ||
    STUDY2_MCMC_TARGET_N_LATENT < 2L ||
    !is.finite(STUDY2_MEAS_RHAT_LIMIT) || STUDY2_MEAS_RHAT_LIMIT <= 1 ||
    !is.finite(STUDY2_MEAS_BULK_ESS_LIMIT) ||
    STUDY2_MEAS_BULK_ESS_LIMIT < 1 ||
    !is.finite(STUDY2_MEAS_TAIL_ESS_LIMIT) ||
    STUDY2_MEAS_TAIL_ESS_LIMIT < 1 ||
    length(STUDY2_ABLATION_GP_N_STARTS) != 1L ||
    is.na(STUDY2_ABLATION_GP_N_STARTS) ||
    STUDY2_ABLATION_GP_N_STARTS < 1L ||
    length(STUDY2_ABLATION_GP_MAXIT) != 1L ||
    is.na(STUDY2_ABLATION_GP_MAXIT) ||
    STUDY2_ABLATION_GP_MAXIT < 0L ||
    length(measurement_n_iter) != 1L || is.na(measurement_n_iter) ||
    measurement_n_iter < 1L ||
    length(measurement_burn) != 1L || is.na(measurement_burn) ||
    measurement_burn < 0L || measurement_burn >= measurement_n_iter ||
    length(measurement_thin) != 1L || is.na(measurement_thin) ||
    measurement_thin < 1L ||
    length(measurement_n_chains) != 1L || is.na(measurement_n_chains) ||
    (isTRUE(STUDY2_RUN_ABLATIONS) && measurement_n_chains < 2L)) {
  stop("Invalid Study II publication MCMC-gate settings.")
}
measurement_n_iter <- as.integer(measurement_n_iter)
measurement_burn <- as.integer(measurement_burn)
measurement_thin <- as.integer(measurement_thin)
measurement_n_chains <- as.integer(measurement_n_chains)
STUDY2_MCMC_TARGET_N_POINTS <- as.integer(STUDY2_MCMC_TARGET_N_POINTS)
STUDY2_ABLATION_GP_N_STARTS <- as.integer(STUDY2_ABLATION_GP_N_STARTS)
STUDY2_ABLATION_GP_MAXIT <- as.integer(STUDY2_ABLATION_GP_MAXIT)

if (!exists("STUDY2_PARALLEL_LEVEL")) STUDY2_PARALLEL_LEVEL <- "chains"
if (!STUDY2_PARALLEL_LEVEL %in% c("chains", "replications", "none", "hybrid")) {
  stop("STUDY2_PARALLEL_LEVEL must be chains, replications, or none.")
}
parallel_chains <- STUDY2_PARALLEL_LEVEL %in% c("chains", "hybrid")
chain_cores <- if (STUDY2_PARALLEL_LEVEL == "hybrid") {
  mixedgp_as_integer_strict(STUDY2_CHAIN_WORKERS, "chain workers", 1L, 1L)
} else NULL

scenario_calib_grid <- function(scenario) {
  if (identical(scenario, "primary")) {
    return(as.integer(STUDY2_PRIMARY_CALIB_GRID))
  }
  as.integer(STUDY2_CONTRAST_CALIB)
}

scenario_label <- function(scenario) {
  c(
    primary = "Primary interactive",
    latent_additive_control = "Latent-additive control",
    high_uncertainty = "High latent uncertainty",
    logistic_misspec = "Logistic-error robustness"
  )[[scenario]]
}

FIG_DIR <- file.path(STUDY2_OUT_PREFIX, "figures")
TAB_DIR <- file.path(STUDY2_OUT_PREFIX, "tables")
RES_DIR <- file.path(STUDY2_OUT_PREFIX, "results", "study2_manuscript_v5")
REP_DIR <- file.path(RES_DIR, "mc_replications")
FIT_DIR <- file.path(RES_DIR, "mc_fits")

for (dd in c(FIG_DIR, TAB_DIR, RES_DIR, REP_DIR, FIT_DIR)) {
  dir.create(dd, showWarnings = FALSE, recursive = TRUE)
}

scenario_code <- c(
  primary = "P",
  latent_additive_control = "A",
  high_uncertainty = "H",
  logistic_misspec = "L"
)
scenario_tag <- paste(unname(scenario_code[STUDY2_SCENARIOS]), collapse = "")
CACHE_SPEC <- list(
  schema = "s2v18_nonfatal_reference_evaluation",
  oracle_mean_method = STUDY2_ORACLE_MEAN_METHOD,
  recovery_model = STUDY2_RECOVERY_MODEL,
  design_tag = STUDY2_DESIGN_TAG,
  study2_config = STUDY2_CONFIG,
  scenario_code = scenario_tag,
  scenarios = STUDY2_SCENARIOS,
  ablation_scenarios = STUDY2_ABLATION_SCENARIOS,
  n_train = n_train,
  n_test = n_test,
  n_rep = n_rep,
  q = length(m_vec),
  m = m_vec[1L],
  primary_calibration_grid = STUDY2_PRIMARY_CALIB_GRID,
  contrast_calibration = STUDY2_CONTRAST_CALIB,
  mcmc = list(
    n_iter = mc_n_iter, burn = mc_burn, thin = mc_thin,
    n_chains = mc_n_chains,
    rhat_limit = STUDY2_MCMC_RHAT_LIMIT,
    raw_ess_limit = STUDY2_MCMC_RAW_ESS_LIMIT,
    target_bulk_ess_limit = STUDY2_MCMC_TARGET_BULK_ESS_LIMIT,
    target_tail_ess_limit = STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT,
    target_n_points = STUDY2_MCMC_TARGET_N_POINTS,
    target_draws_per_chain = STUDY2_MCMC_TARGET_DRAWS_PER_CHAIN,
    target_integration_draws = STUDY2_MCMC_TARGET_N_LATENT
  ),
  prediction = list(
    n_draw = n_pred_draw,
    latent_sampler = predictive_latent_sampler,
    oracle_pool = n_oracle_pool
  ),
  mean_recovery = list(
    n_eval = n_m_eval, n_draw = n_m_draw,
    n_latent = n_m_latent, n_truth = n_m_truth
  ),
  measurement = list(
    n_iter = measurement_n_iter, burn = measurement_burn,
    thin = measurement_thin, n_chains = measurement_n_chains,
    rhat_limit = STUDY2_MEAS_RHAT_LIMIT,
    bulk_ess_limit = STUDY2_MEAS_BULK_ESS_LIMIT,
    tail_ess_limit = STUDY2_MEAS_TAIL_ESS_LIMIT
  ),
  ablation_gp = list(
    n_starts = STUDY2_ABLATION_GP_N_STARTS,
    maxit = STUDY2_ABLATION_GP_MAXIT
  ),
  task_flags = list(
    run_ablations = isTRUE(STUDY2_RUN_ABLATIONS),
    evaluate_f = isTRUE(STUDY2_EVALUATE_F),
    evaluate_u = isTRUE(STUDY2_EVALUATE_U)
  ),
  published_competitors = STUDY2_PUBLISHED_COMPETITORS,
  competitor_controls = study2_competitor_controls
)
CACHE_SPEC$fingerprint <- mixedgp_object_fingerprint(CACHE_SPEC)
CACHE_TAG <- paste0(
  "s2v18-q", length(m_vec), "-", substr(CACHE_SPEC$fingerprint, 1L, 16L)
)
saveRDS(
  CACHE_SPEC,
  file.path(RES_DIR, paste0("cache_spec_", CACHE_TAG, ".rds")),
  version = 3L
)
writeLines(
  capture.output(dput(CACHE_SPEC)),
  file.path(RES_DIR, paste0("cache_spec_", CACHE_TAG, ".R"))
)

method_levels <- c("Oracle", "EIV-GP", "UC-GP", "LVGP", "EzGP")
method_cols <- c(
  "Oracle" = "black",
  "EIV-GP" = "firebrick",
  "UC-GP" = "steelblue4",
  "LVGP" = "darkorange3",
  "EzGP" = "purple4"
)
metric_levels <- c(
  "RMSE", "MAE", "CRPS", "NLPD", "Coverage95", "Width95",
  "IntervalScore95"
)

safe_se <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  sd(x) / sqrt(length(x))
}

format_mean_se <- function(mean, se, digits = 3L) {
  if (!is.finite(mean)) return("--")
  if (!is.finite(se)) return(sprintf(paste0("%.", digits, "f"), mean))
  paste0(
    sprintf(paste0("%.", digits, "f"), mean),
    " (", sprintf(paste0("%.", digits, "f"), se), ")"
  )
}

save_plot <- function(path_no_ext, plot, width, height, dpi = 320) {
  if (isTRUE(STUDY2_SAVE_PDF)) {
    ggsave(paste0(path_no_ext, ".pdf"), plot, width = width, height = height)
  }
  if (isTRUE(STUDY2_SAVE_PNG)) {
    ggsave(
      paste0(path_no_ext, ".png"), plot,
      width = width, height = height, dpi = dpi, bg = "white"
    )
  }
  invisible(plot)
}

write_csv_safe <- function(x, path) {
  if (ncol(x) == 0L) x <- data.frame(note = character(0))
  write.csv(x, path, row.names = FALSE)
  invisible(x)
}

extract_study2_diagnostics <- function(fit, rep_id, scenario, n_calib) {
  out <- fit$diagnostics$summary
  if (is.null(out)) return(data.frame())
  out <- as.data.frame(out)
  out$rep <- rep_id
  out$scenario <- scenario
  out$n_calib <- n_calib
  out
}

study2_sampler_control_settings <- function(fit,
                                             prediction_seed,
                                             mean_seed,
                                             prediction_sampler_used,
                                             mean_sampler_used) {
  diag <- fit$diagnostics$summary
  c(
    fit$control,
    list(
      sampler_strategy = fit$sampler_strategy,
      identification = fit$data$ident,
      kernel = fit$kernel$name,
      matern_nu = fit$kernel$matern_nu,
      n_iter = diag$n_iter,
      burn = diag$burn,
      thin = diag$thin,
      n_chains = diag$n_chains,
      saved_per_chain = diag$saved_per_chain,
      parallel_backend = diag$parallel_backend,
      parallel_cores = diag$parallel_cores,
      initialization_rule = paste(
        "threshold-compatible ordinal-score latent initialization with",
        "independent seeded perturbations; see",
        "fit_eivgp_ordprobit_fb()"
      ),
      chain_seeds = fit$mcmc$chain_stats$seed,
      covariance_jitter = diag$covariance_jitter,
      forms_explicit_covariance_inverse =
        diag$forms_explicit_covariance_inverse,
      prospective_latent_sampler_requested = predictive_latent_sampler,
      predictive_latent_sampler_used = prediction_sampler_used,
      mean_latent_sampler_used = mean_sampler_used,
      diagnostic_gibbs_sweeps = diagnostic_n_new_latent_gibbs,
      rejection_max_batches = rejection_max_batches,
      predictive_seed = prediction_seed,
      mean_seed = mean_seed,
      predictive_draws_requested = n_pred_draw,
      mean_evaluation_points = n_m_eval,
      mean_posterior_draws_requested = n_m_draw,
      mean_latent_integration_draws = n_m_latent,
      mean_truth_latent_draws = n_m_truth
    )
  )
}

study2_sampler_control_rows <- function(fit,
                                        rep_id,
                                        scenario,
                                        n_calib,
                                        prediction_seed,
                                        mean_seed,
                                        prediction_sampler_used,
                                        mean_sampler_used) {
  encode <- function(value) {
    if (is.numeric(value)) {
      return(paste(
        format(value, digits = 17L, scientific = FALSE),
        collapse = ";"
      ))
    }
    paste(as.character(value), collapse = ";")
  }
  settings <- study2_sampler_control_settings(
    fit = fit,
    prediction_seed = prediction_seed,
    mean_seed = mean_seed,
    prediction_sampler_used = prediction_sampler_used,
    mean_sampler_used = mean_sampler_used
  )
  data.frame(
    rep = as.integer(rep_id),
    scenario = as.character(scenario),
    n_calib = as.integer(n_calib),
    setting = names(settings),
    value = vapply(settings, encode, character(1L)),
    stringsAsFactors = FALSE
  )
}

study2_finite_max <- function(x) {
  if (length(x) == 0L || any(!is.finite(x))) NA_real_ else max(x)
}

study2_finite_min <- function(x) {
  if (length(x) == 0L || any(!is.finite(x))) NA_real_ else min(x)
}

study2_ablation_status <- function(method,
                                   n_calib,
                                   status,
                                   message = "",
                                   elapsed_seconds = NA_real_) {
  data.frame(
    method = method,
    n_calib = n_calib,
    status = status,
    message = message,
    elapsed_seconds = elapsed_seconds,
    stringsAsFactors = FALSE
  )
}

study2_optimizer_attempt_rows <- function(attempts,
                                          method,
                                          n_calib,
                                          selected_start = NA_integer_) {
  if (is.null(attempts) || !is.data.frame(attempts) || nrow(attempts) == 0L) {
    return(data.frame())
  }
  out <- attempts
  out$method <- method
  out$n_calib <- n_calib
  out$selected <- !is.na(selected_start) & out$start == selected_start
  leading <- c("method", "n_calib", "start", "selected")
  out[, c(leading, setdiff(names(out), leading)), drop = FALSE]
}

############################################################
## One scenario-replication pair
############################################################

}


# ---- simulations/study2/02_study2_monte_carlo.R ----
eivgp_study2_driver <- function() {
## Fitting/evaluation only. Reporting lives in report_study2_results.R.
eval(body(eivgp_study2_setup), envir = environment())

run_one_study2_replication <- function(rep_id, scenario) {
  calib_grid <- scenario_calib_grid(scenario)
  scenario_id <- match(scenario, allowed_scenarios)
  ## The three Gaussian scenarios use common random numbers. For a fixed
  ## replication, they share X, U, Gaussian score innovations, and response
  ## innovations. The logistic robustness scenario shares X and U; its score
  ## errors necessarily use a different random-number transformation.
  data_seed_base <- 1000000L + 10000L * rep_id
  fit_seed_base <- 1000000L * scenario_id + 10000L * rep_id

  data_path <- file.path(
    STUDY2_DATA_DIR,
    mixedgp_dataset_filename(
      "study2", rep_id, scenario, n_train, n_test, m_vec[1L], q = length(m_vec)
    )
  )
  if (!file.exists(data_path)) {
    stop(
      "Missing frozen Study II dataset: ", data_path,
      ". Run 12_generate_synthetic_datasets.R first."
    )
  }
  frozen <- load_mixedgp_synthetic_dataset_strict(
    data_path,
    expected = list(
      study = "study2", scenario = scenario, rep_id = rep_id,
      n = n_train, n_test = n_test, q = length(m_vec), m = m_vec[1L],
      calib_grid = calib_grid
    )
  )
  dat <- frozen$data
  train <- dat$train
  test <- dat$test

  pattern_info <- classify_study2_pattern_frequency(train$C, test$C)
  pattern_counts <- as.data.frame(table(pattern_info$pattern_stratum))
  names(pattern_counts) <- c("evaluation_stratum", "n_test")
  pattern_counts$rep <- rep_id
  pattern_counts$scenario <- scenario

  calib_sets <- frozen$calibration_sets[as.character(calib_grid)]
  if (any(vapply(calib_sets, is.null, logical(1)))) {
    stop("Frozen Study II dataset does not contain every required calibration split.")
  }

  reference_status <- list()
  pool_task <- mixedgp_reference_task("oracle_pool", make_oracle_pool_2d(
    true_params=dat$true_params, n_pool=n_oracle_pool, seed=fit_seed_base+3L))
  oracle_pool <- pool_task$value
  reference_status$pool <- pool_task$status
  prediction_task <- mixedgp_reference_task("oracle_prediction", sample_oracle_test_y_2d(
    X_test=test$X, C_test=test$C, true_params=dat$true_params,
    sigma_eps=dat$sigma_eps, n_draw=n_pred_draw, oracle_pool=oracle_pool, seed=fit_seed_base+4L))
  oracle_draws <- prediction_task$value
  reference_status$prediction <- prediction_task$status

  competitor_status <- data.frame(method=character(),status=character(),rep=integer(),scenario=character())

  set.seed(data_seed_base + 150L)
  mean_eval_idx <- sort(sample(
    seq_len(nrow(test$X)), min(n_m_eval, nrow(test$X))
  ))
  target_gate_local_idx <- mixedgp_study2_diagnostic_rows(
    test$X[mean_eval_idx, , drop = FALSE], test$C[mean_eval_idx, , drop = FALSE],
    STUDY2_MCMC_TARGET_N_POINTS
  )
  target_gate_test_rows <- mean_eval_idx[target_gate_local_idx]
  mean_task <- mixedgp_reference_task("oracle_mean", {
    if (STUDY2_ORACLE_MEAN_METHOD == "quadrature") {
      oracle_m0_quadrature_2d(test$X[mean_eval_idx,,drop=FALSE],
        test$C[mean_eval_idx,,drop=FALSE], dat$true_params)
    } else {
      oracle_m0_2d(X=test$X[mean_eval_idx,,drop=FALSE], C=test$C[mean_eval_idx,,drop=FALSE],
        true_params=dat$true_params, oracle_pool=oracle_pool, n_latent=n_m_truth,
        seed=data_seed_base+151L)
    }
  })
  mean_truth <- mean_task$value
  mean_task$status$method <- STUDY2_ORACLE_MEAN_METHOD
  reference_status$mean <- mean_task$status
  mean_truth_rejection <- attr(mean_truth, "rejection_telemetry")
  mean_truth_diagnostics <- attr(mean_truth, "truth_diagnostics")
  tag_truth <- function(x) {
    if (is.null(x) || !nrow(x)) return(data.frame())
    x$rep <- rep_id; x$scenario <- scenario; x
  }
  mean_truth_rejection <- tag_truth(mean_truth_rejection)
  mean_truth_diagnostics <- tag_truth(mean_truth_diagnostics)
  reference_status <- tag_truth(bind_rows(reference_status))
  write.csv(reference_status, file.path(TAB_DIR, sprintf("reference_status_%s_rep%03d.csv",scenario,rep_id)),row.names=FALSE)
  mean_recovery <- list()

  metrics <- if (is.null(oracle_draws)) list() else list(
    Oracle = summarize_predictive_samples_by_pattern(
      draw_mat = oracle_draws,
      y_true = test$y,
      pattern_stratum = pattern_info$pattern_stratum,
      method = "Oracle",
      rep_id = rep_id,
      n_calib = NA_integer_,
      scenario = scenario
    )
  )


  diagnostics <- list()
  target_diagnostics <- list()
  imputation <- list()
  imputation_status <- list()
  ablation_imputation <- list()
  surface <- list()
  ablation_metrics <- list()
  ablation_surface <- list()
  measurement_diagnostics <- list()
  measurement_parameter_diagnostics <- list()
  ablation_optimizer_attempts <- list()
  ablation_status <- list()
  sampler_controls <- list()
  sampler_control_manifest <- list()

  for (n_calib in calib_grid) {
    message(
      "Study II scenario=", scenario,
      ", replication=", rep_id,
      ", |O|=", n_calib
    )

    fit_file <- file.path(FIT_DIR, sprintf("recovery_fit_%s_rep%03d_cal%03d_%s.rds",scenario,rep_id,n_calib,CACHE_TAG))
    fit_identity <- list(cache_spec=CACHE_SPEC, data_md5=unname(tools::md5sum(data_path)),
      rep=rep_id, scenario=scenario, calibration=n_calib)
    saved_fit <- if (file.exists(fit_file)) mixedgp_read_artifact(fit_file) else NULL
    if (!is.null(saved_fit) && !identical(saved_fit$identity,fit_identity)) stop("Saved-fit identity mismatch: ",fit_file)
    fit <- if (!is.null(saved_fit)) saved_fit$fit else fit_eivgp_ordprobit_fb(
      X_raw = train$X,
      y_raw = train$y,
      C_ord = train$C,
      U_obs = train$U,
      calib_idx = calib_sets[[as.character(n_calib)]],
      U_true_eval = train$U,
      d = d_latent,
      m_vec = m_vec,
      ident = ident_method,
      priors = if(is.null(STUDY2_RECOVERY_MODEL))list()else STUDY2_RECOVERY_MODEL$priors,
      sampler_control = if(is.null(STUDY2_RECOVERY_MODEL))list()else STUDY2_RECOVERY_MODEL$control,
      kernel = if(is.null(STUDY2_RECOVERY_MODEL))"se"else STUDY2_RECOVERY_MODEL$kernel$name,
      matern_nu = if(is.null(STUDY2_RECOVERY_MODEL))2.5 else STUDY2_RECOVERY_MODEL$kernel$matern_nu,
      n_iter = mc_n_iter,
      burn = mc_burn,
      thin = mc_thin,
      n_chains = mc_n_chains,
      preset = mc_preset,
      store_scores = FALSE,
      seed = fit_seed_base + 1000L + n_calib,
      parallel_chains = parallel_chains,
      n_cores = chain_cores,
      verbose = FALSE
    )

    # Save immediately: downstream diagnostics/evaluation must not lose MCMC work.
    if (isTRUE(STUDY2_CHECKPOINT_FITS) && is.null(saved_fit))
      mixedgp_save_artifact(list(identity=fit_identity,fit=fit),fit_file)
    rm(saved_fit)
    if(!is.null(STUDY2_ESS_MIN)) {
      panel<-list(X=test$X[target_gate_test_rows,,drop=FALSE],
        C=test$C[target_gate_test_rows,,drop=FALSE],
        U=if(n_calib>0L && isTRUE(STUDY2_EVALUATE_F))test$U[target_gate_test_rows,,drop=FALSE]else NULL,
        n_latent=STUDY2_MCMC_TARGET_N_LATENT,seed=fit_seed_base+2900L)
      extension_file<-paste0(fit_file,".ess-state.rds")
      if(file.exists(extension_file)) {
        resumed<-mixedgp_read_artifact(extension_file)
        if(!identical(resumed$settings$X,panel$X)||!identical(resumed$settings$C,panel$C)||
           !identical(resumed$settings$U,panel$U))stop("Scientific continuation panel mismatch.")
        fit<-resumed$fit
      }
      fit<-mixedgp_finish_ess(fit,STUDY2_ESS_MIN,panel,parallel_chains,chain_cores,extension_file)
      mixedgp_save_artifact(list(identity=fit_identity,fit=fit),fit_file)
      mixedgp_save_scientific_trace(fit,file.path(FIT_DIR,sprintf("scientific-trace-%s-%03d-%03d",scenario,rep_id,n_calib)),
        sprintf("%s | q=%d | replication %d | calibration %d",scenario,length(m_vec),rep_id,n_calib))
      unlink(extension_file)
    }
    draw_ids <- seq_len(dim(fit$mcmc$samples_U)[1])
    if (length(draw_ids) > n_pred_draw) {
      draw_ids <- draw_ids[
        unique(round(seq(1, length(draw_ids), length.out = n_pred_draw)))
      ]
    }
    predictive_seed <- fit_seed_base + 2000L + n_calib
    mean_seed <- fit_seed_base + 2500L + n_calib

    eiv_draws <- sample_eiv_test_y_ordprobit_fb(
      X_test_raw = test$X,
      C_test = test$C,
      fit_obj = fit,
      draw_ids = draw_ids,
      n_per_draw = 1L,
      latent_sampler = predictive_latent_sampler,
      n_new_latent_gibbs = diagnostic_n_new_latent_gibbs,
      rejection_max_batches = rejection_max_batches,
      seed = predictive_seed
    )

    metrics[[paste0("EIV_", n_calib)]] <-
      summarize_predictive_samples_by_pattern(
        draw_mat = eiv_draws,
        y_true = test$y,
        pattern_stratum = pattern_info$pattern_stratum,
        method = "EIV-GP",
        rep_id = rep_id,
        n_calib = n_calib,
        scenario = scenario
      )

    mixedgp_atomic_save_rds(list(identity=fit_identity,metrics=metrics[[paste0("EIV_",n_calib)]]),
      file.path(RES_DIR,sprintf("predictive_checkpoint_%s_rep%03d_cal%03d_%s.rds",scenario,rep_id,n_calib,CACHE_TAG)))

    m_draw_ids <- draw_ids
    if (length(m_draw_ids) > n_m_draw) {
      m_draw_ids <- m_draw_ids[
        unique(round(seq(1, length(m_draw_ids), length.out = n_m_draw)))
      ]
    }
    m_draws <- sample_eiv_m_given_xc_fb(
      X_test_raw = test$X[mean_eval_idx, , drop = FALSE],
      C_test = test$C[mean_eval_idx, , drop = FALSE],
      fit_obj = fit,
      draw_ids = m_draw_ids,
      n_latent = n_m_latent,
      include_process_uncertainty = TRUE,
      joint = FALSE,
      return_components = TRUE,
      latent_sampler = predictive_latent_sampler,
      n_new_latent_gibbs = diagnostic_n_new_latent_gibbs,
      rejection_max_batches = rejection_max_batches,
      seed = mean_seed
    )
    control_key <- as.character(n_calib)
    sampler_controls[[control_key]] <- study2_sampler_control_settings(
      fit = fit,
      prediction_seed = predictive_seed,
      mean_seed = mean_seed,
      prediction_sampler_used = attr(eiv_draws, "latent_sampler"),
      mean_sampler_used = attr(m_draws, "latent_sampler")
    )
    sampler_control_manifest[[control_key]] <-
      study2_sampler_control_rows(
        fit = fit,
        rep_id = rep_id,
        scenario = scenario,
        n_calib = n_calib,
        prediction_seed = predictive_seed,
        mean_seed = mean_seed,
        prediction_sampler_used = attr(eiv_draws, "latent_sampler"),
        mean_sampler_used = attr(m_draws, "latent_sampler")
      )
    if (!is.null(mean_truth)) mean_recovery[[paste0("EIV_", n_calib)]] <-
      summarize_mean_recovery_2d(
        m_draws,
        m_true = mean_truth,
        method = "EIV-GP",
        rep_id = rep_id,
        n_calib = n_calib,
        scenario = scenario,
        valid_function_draws = TRUE
      )

    target_series <- if(!is.null(fit$scientific_ess))fit$scientific_ess$target_series else mixedgp_study2_target_series(
      fit, test$X[target_gate_test_rows, , drop = FALSE],
      test$C[target_gate_test_rows, , drop = FALSE],
      U = if (n_calib > 0L && isTRUE(STUDY2_EVALUATE_F)) {
        test$U[target_gate_test_rows, , drop = FALSE]
      } else NULL,
      max_draws_per_chain = STUDY2_MCMC_TARGET_DRAWS_PER_CHAIN,
      n_latent = STUDY2_MCMC_TARGET_N_LATENT,
      seed = fit_seed_base + 2900L
    )
    invariant_series <- mixedgp_study2_invariant_series(fit)
    raw_series <- mixedgp_study2_raw_series(fit)
    summarize_gate <- function(series, ess_min, tail_min) {
      mixedgp_summarize_diagnostic_series(
        series, STUDY2_MCMC_RHAT_LIMIT, ess_min, tail_min
      )
    }
    functional_diag <- summarize_gate(target_series, STUDY2_MCMC_TARGET_BULK_ESS_LIMIT,
                                       STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT)
    invariant_diag <- summarize_gate(invariant_series, STUDY2_MCMC_TARGET_BULK_ESS_LIMIT,
                                     STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT)
    raw_diag <- summarize_gate(raw_series, STUDY2_MCMC_RAW_ESS_LIMIT,
                               STUDY2_MCMC_RAW_ESS_LIMIT)
    functional_diag$target <- sub("\\[.*$", "", functional_diag$parameter)
    invariant_diag$target <- "ordinal_measurement_invariant"
    raw_diag$target <- "raw_coordinate"
    functional_diag$used_for_gate <- TRUE
    invariant_diag$used_for_gate <- is.null(STUDY2_ESS_MIN)
    raw_diag$used_for_gate <- is.null(STUDY2_ESS_MIN) && n_calib > 0L
    target_diag <- bind_rows(functional_diag, invariant_diag, raw_diag)
    target_diag$rep <- rep_id
    target_diag$scenario <- scenario
    target_diag$n_calib <- n_calib
    target_diagnostics[[as.character(n_calib)]] <- target_diag

    raw_coordinate_pass <- mixedgp_diagnostic_table_pass(
      raw_diag, STUDY2_MCMC_RHAT_LIMIT, STUDY2_MCMC_RAW_ESS_LIMIT,
      STUDY2_MCMC_RAW_ESS_LIMIT, expected_parameters = names(raw_series)
    )
    target_functional_pass <- mixedgp_diagnostic_table_pass(
      functional_diag, STUDY2_MCMC_RHAT_LIMIT, STUDY2_MCMC_TARGET_BULK_ESS_LIMIT,
      STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT, expected_parameters = names(target_series)
    )
    invariant_measurement_pass <- mixedgp_diagnostic_table_pass(
      invariant_diag, STUDY2_MCMC_RHAT_LIMIT, STUDY2_MCMC_TARGET_BULK_ESS_LIMIT,
      STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT, expected_parameters = names(invariant_series)
    )
    diag_row <- extract_study2_diagnostics(fit, rep_id, scenario, n_calib)
    if (nrow(diag_row) != 1L) stop("Fit did not return one diagnostic summary.")
    predictive_diag <- functional_diag[
      functional_diag$target == "predictive_second_moment", , drop = FALSE
    ]
    mean_diag <- functional_diag[
      functional_diag$target == "m_conditional_mean", , drop = FALSE
    ]
    diag_row$mcmc_gate_rule <- if (n_calib == 0L) {
      "targets_and_measurement_invariants_without_calibration"
    } else {
      "targets_measurement_invariants_and_free_coordinates_with_calibration"
    }
    diag_row$mcmc_rhat_limit <- STUDY2_MCMC_RHAT_LIMIT
    diag_row$mcmc_raw_ess_limit <- STUDY2_MCMC_RAW_ESS_LIMIT
    diag_row$mcmc_target_bulk_ess_limit <- STUDY2_MCMC_TARGET_BULK_ESS_LIMIT
    diag_row$mcmc_target_tail_ess_limit <- STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT
    diag_row$mcmc_target_n_points <- length(target_gate_test_rows)
    diag_row$mcmc_target_draws_per_chain <- STUDY2_MCMC_TARGET_DRAWS_PER_CHAIN
    diag_row$mcmc_target_integration_draws <- STUDY2_MCMC_TARGET_N_LATENT
    diag_row$raw_coordinate_pass <- raw_coordinate_pass
    diag_row$target_functional_pass <- target_functional_pass
    diag_row$invariant_measurement_pass <- invariant_measurement_pass
    diag_row$predictive_max_rhat <- study2_finite_max(predictive_diag$rhat)
    diag_row$predictive_min_bulk_ess <- study2_finite_min(predictive_diag$ess_bulk)
    diag_row$predictive_min_tail_ess <- study2_finite_min(predictive_diag$ess_tail)
    diag_row$m_max_rhat <- study2_finite_max(mean_diag$rhat)
    diag_row$m_min_bulk_ess <- study2_finite_min(mean_diag$ess_bulk)
    diag_row$m_min_tail_ess <- study2_finite_min(mean_diag$ess_tail)
    diag_row$invariant_max_rhat <- study2_finite_max(invariant_diag$rhat)
    diag_row$invariant_min_bulk_ess <- study2_finite_min(invariant_diag$ess_bulk)
    diag_row$invariant_min_tail_ess <- study2_finite_min(invariant_diag$ess_tail)
    diag_row$mcmc_pass <- target_functional_pass && invariant_measurement_pass &&
      (n_calib == 0L || raw_coordinate_pass)
    if(!is.null(STUDY2_ESS_MIN)) {
      diag_row$scientific_ess_pass<-mixedgp_scientific_ess_pass(functional_diag,STUDY2_ESS_MIN)
      diag_row$mcmc_pass<-diag_row$scientific_ess_pass
      diag_row$mcmc_gate_rule<-paste0("scientific_bulk_and_tail_ESS_strictly_greater_than_",STUDY2_ESS_MIN)
    }
    diag_row$diagnostic_warning <- !isTRUE(diag_row$mcmc_pass)
    diag_row$diagnostic_advice <- if (isTRUE(diag_row$mcmc_pass)) "" else mixedgp_simulation_diagnostic_advice()
    diagnostics[[as.character(n_calib)]] <- diag_row
    if (!isTRUE(diag_row$mcmc_pass)) {
      # One canonical fit; warnings contain diagnostics and a relative reference.
      if (!file.exists(fit_file))
        mixedgp_save_artifact(list(identity = fit_identity, fit = fit), fit_file)
      mixedgp_atomic_save_rds(list(diagnostics = diag_row,
        fit_file = file.path("mc_fits", basename(fit_file))),
        file.path(RES_DIR, sprintf("flagged_fit_%s_rep%03d_cal%03d.rds",
                                  scenario, rep_id, n_calib)))
      warning("Study II MCMC diagnostics flagged replication ", rep_id,
        ", calibration ", n_calib, ". ", mixedgp_simulation_diagnostic_advice(),
        call. = FALSE)
    }

    if (isTRUE(STUDY2_EVALUATE_U)) {
      imputation[[paste0("EIV_training_", n_calib)]] <-
        study2_latent_imputation_metrics(
          fit = fit,
          U_true = train$U,
          rep_id = rep_id,
          n_calib = n_calib,
          scenario = scenario
        )
      imputation_status[[paste0("EIV_training_", n_calib)]] <-
        study2_latent_imputation_status(
          fit = fit,
          method = "EIV-GP",
          target = "training_missing_U",
          rep_id = rep_id,
          n_calib = n_calib,
          scenario = scenario,
          n_units = length(fit$data$miss_idx)
        )
      imputation[[paste0("EIV_prospective_", n_calib)]] <-
        study2_eiv_prospective_imputation_metrics(
          fit = fit,
          C_new = test$C[mean_eval_idx, , drop = FALSE],
          U_true = test$U[mean_eval_idx, , drop = FALSE],
          rep_id = rep_id,
          n_calib = n_calib,
          scenario = scenario,
          max_draw = n_m_draw,
          latent_sampler = predictive_latent_sampler,
          n_new_latent_gibbs = diagnostic_n_new_latent_gibbs,
          rejection_max_batches = rejection_max_batches,
          seed = fit_seed_base + 2700L + n_calib
        )
      imputation_status[[paste0("EIV_prospective_", n_calib)]] <-
        study2_latent_imputation_status(
          fit = fit,
          method = "EIV-GP",
          target = "prospective_U_given_C",
          rep_id = rep_id,
          n_calib = n_calib,
          scenario = scenario,
          n_units = length(mean_eval_idx)
        )
    }

    if (isTRUE(STUDY2_EVALUATE_F)) {
      surface[[as.character(n_calib)]] <-
        study2_eiv_surface_recovery_metrics(
          fit = fit,
          scenario = scenario,
          rep_id = rep_id,
          n_calib = n_calib,
          grid_n = if (STUDY2_CONFIG == "quick") 15L else 31L,
          max_draw = if (STUDY2_CONFIG == "quick") 40L else 200L
        )
    }

    if (isTRUE(STUDY2_SAVE_REP_FITS)) {
      if (!file.exists(fit_file))
        mixedgp_save_artifact(list(identity = fit_identity, fit = fit), fit_file)
    }

    rm(fit, eiv_draws, m_draws)
    invisible(gc())
  }

  ##########################################################
  ## Appendix-only ablations (kept out of the main table)
  ##########################################################

  if (isTRUE(STUDY2_RUN_ABLATIONS) &&
      scenario %in% STUDY2_ABLATION_SCENARIOS) {
    if (isTRUE(STUDY2_EVALUATE_F)) {
      start_full_u <- proc.time()[3]
      full_u_gp <- tryCatch(
        fit_study2_latent_gp(
          train$X, train$y, train$U,
          gp_n_starts = STUDY2_ABLATION_GP_N_STARTS,
          gp_seed = fit_seed_base + 8000L,
          gp_maxit = STUDY2_ABLATION_GP_MAXIT
        ),
        error = function(e) e
      )
      elapsed_full_u <- proc.time()[3] - start_full_u

      if (inherits(full_u_gp, "error")) {
        ablation_optimizer_attempts[["Full-U GP"]] <-
          study2_optimizer_attempt_rows(
            full_u_gp$optimizer_attempts,
            "Full-U GP",
            NA_integer_
          )
        ablation_status[["Full-U GP"]] <- study2_ablation_status(
          "Full-U GP", NA_integer_, "failed",
          conditionMessage(full_u_gp), elapsed_full_u
        )
      } else {
        ablation_optimizer_attempts[["Full-U GP"]] <-
          study2_optimizer_attempt_rows(
            full_u_gp$optimizer_attempts,
            "Full-U GP",
            NA_integer_,
            full_u_gp$selected_start
          )
        ablation_status[["Full-U GP"]] <- study2_ablation_status(
          "Full-U GP", NA_integer_, "success", "", elapsed_full_u
        )
        ablation_surface[["Full-U GP"]] <-
          study2_latent_gp_surface_metrics(
            gp_fit = full_u_gp,
            scenario = scenario,
            method = "Full-U GP",
            rep_id = rep_id,
            n_calib = NA_integer_,
            grid_n = if (STUDY2_CONFIG == "quick") 15L else 31L
          )
      }
    }

    for (n_calib in calib_grid) {
      key <- as.character(n_calib)
      calib_idx <- calib_sets[[key]]

      start_measurement <- proc.time()[3]
      measurement_fit <- tryCatch(
        fit_ordinalprobit_measurement_fb(
          C_ord = train$C,
          U_obs = train$U,
          calib_idx = calib_idx,
          d = d_latent,
          m_vec = m_vec,
          ident = ident_method,
          n_iter = measurement_n_iter,
          burn = measurement_burn,
          thin = measurement_thin,
          n_chains = measurement_n_chains,
          seed = fit_seed_base + 5000L + n_calib,
          parallel_chains = parallel_chains,
          n_cores = chain_cores,
          rhat_limit = STUDY2_MEAS_RHAT_LIMIT,
          bulk_ess_limit = STUDY2_MEAS_BULK_ESS_LIMIT,
          tail_ess_limit = STUDY2_MEAS_TAIL_ESS_LIMIT,
          verbose = FALSE
        ),
        error = function(e) e
      )
      elapsed_measurement <- proc.time()[3] - start_measurement

      if (inherits(measurement_fit, "error")) {
        msg <- conditionMessage(measurement_fit)
        measurement_status_stub <- list(data = list(
          U_obs = train$U,
          calib_idx = calib_idx,
          miss_idx = setdiff(seq_len(nrow(train$U)), calib_idx),
          d = d_latent
        ))
        for (task_spec in list(
          list(
            target = "training_missing_U",
            n_units = length(measurement_status_stub$data$miss_idx)
          ),
          list(
            target = "prospective_U_given_C",
            n_units = length(mean_eval_idx)
          )
        )) {
          task_status <- study2_latent_imputation_status(
            fit = measurement_status_stub,
            method = "Ordinal model (no Y)",
            target = task_spec$target,
            rep_id = rep_id,
            n_calib = n_calib,
            scenario = scenario,
            n_units = task_spec$n_units
          )
          if (isTRUE(task_status$task_eligible)) {
            task_status$status <- "failed"
            task_status$task_eligible <- FALSE
            task_status$reason <- msg
          }
          imputation_status[[paste0(
            "measurement_", task_spec$target, "_", key
          )]] <- task_status
        }
        ablation_status[[paste0("measurement_", key)]] <-
          study2_ablation_status(
            "Response-free measurement model", n_calib, "failed",
            msg, elapsed_measurement
          )
        ablation_status[[paste0("PI_", key)]] <-
          study2_ablation_status(
            "PI-GP", n_calib, "failed", msg, NA_real_
          )
        ablation_status[[paste0("CC_", key)]] <-
          study2_ablation_status(
            "CC-GP", n_calib,
            if (n_calib == 0L) "not_applicable" else "failed",
            if (n_calib == 0L) "No complete cases." else msg,
            NA_real_
          )
        next
      }

      measurement_diag <- measurement_fit$diagnostics
      measurement_diag$rep <- rep_id
      measurement_diag$scenario <- scenario
      measurement_diag$n_calib <- n_calib
      measurement_diagnostics[[key]] <- measurement_diag

      key_parameter_diag <- measurement_fit$diagnostic_parameters$key
      latent_parameter_diag <- measurement_fit$diagnostic_parameters$latent_rhat
      if (nrow(latent_parameter_diag) > 0L) {
        latent_parameter_diag$block <- "latent"
      }
      parameter_diag <- bind_rows(
        key_parameter_diag[, c(
          "parameter", "block", "rhat", "ess_bulk", "ess_tail"
        )],
        latent_parameter_diag[, intersect(
          c("parameter", "block", "rhat", "ess_bulk", "ess_tail"),
          names(latent_parameter_diag)
        )]
      )
      parameter_diag$rep <- rep_id
      parameter_diag$scenario <- scenario
      parameter_diag$n_calib <- n_calib
      measurement_parameter_diagnostics[[key]] <- parameter_diag

      measurement_training_status <- study2_latent_imputation_status(
        fit = measurement_fit,
        method = "Ordinal model (no Y)",
        target = "training_missing_U",
        rep_id = rep_id,
        n_calib = n_calib,
        scenario = scenario,
        n_units = length(measurement_fit$data$miss_idx)
      )
      measurement_prospective_status <- study2_latent_imputation_status(
        fit = measurement_fit,
        method = "Ordinal model (no Y)",
        target = "prospective_U_given_C",
        rep_id = rep_id,
        n_calib = n_calib,
        scenario = scenario,
        n_units = length(mean_eval_idx)
      )

      measurement_warning <- !isTRUE(measurement_diag$convergence_pass)
      measurement_advice <- if (measurement_warning) mixedgp_simulation_diagnostic_advice() else ""
      if (measurement_warning) {
        mixedgp_save_artifact(list(fit = measurement_fit, diagnostics = measurement_diag),
          file.path(RES_DIR, sprintf("flagged_measurement_%s_rep%03d_cal%03d.rds",
                                    scenario, rep_id, n_calib)))
        warning("Response-free measurement model has diagnostic warnings. ",
                measurement_advice, call. = FALSE)
      }
      measurement_training_status$diagnostic_warning <- measurement_warning
      measurement_training_status$diagnostic_advice <- measurement_advice
      measurement_prospective_status$diagnostic_warning <- measurement_warning
      measurement_prospective_status$diagnostic_advice <- measurement_advice
      imputation_status[[paste0("measurement_training_", key)]] <- measurement_training_status
      imputation_status[[paste0("measurement_prospective_", key)]] <- measurement_prospective_status
      ablation_status[[paste0("measurement_", key)]] <-
        study2_ablation_status("Response-free measurement model", n_calib,
          if (measurement_warning) "completed_with_diagnostic_warnings" else "success",
          measurement_advice, elapsed_measurement)

      ablation_imputation[[paste0("measurement_training_", key)]] <-
        study2_measurement_training_imputation_metrics(
          measurement_fit = measurement_fit,
          U_true = train$U,
          rep_id = rep_id,
          n_calib = n_calib,
          scenario = scenario
        )
      ablation_imputation[[paste0("measurement_prospective_", key)]] <-
        study2_measurement_prospective_imputation_metrics(
          measurement_fit = measurement_fit,
          C_new = test$C[mean_eval_idx, , drop = FALSE],
          U_true = test$U[mean_eval_idx, , drop = FALSE],
          rep_id = rep_id,
          n_calib = n_calib,
          scenario = scenario,
          max_draw = n_m_draw,
          latent_sampler = predictive_latent_sampler,
          n_gibbs = diagnostic_n_new_latent_gibbs,
          rejection_max_batches = rejection_max_batches,
          seed = fit_seed_base + 5600L + n_calib
        )

      start_pi <- proc.time()[3]
      pi_fit <- tryCatch(
        fit_study2_pi_gp(
          train$X,
          train$y,
          measurement_fit,
          gp_n_starts = STUDY2_ABLATION_GP_N_STARTS,
          gp_seed = fit_seed_base + 6100L + n_calib,
          gp_maxit = STUDY2_ABLATION_GP_MAXIT
        ),
        error = function(e) e
      )
      pi_result <- if (inherits(pi_fit, "error")) {
        pi_fit
      } else {
        ablation_optimizer_attempts[[paste0("PI_", key)]] <-
          study2_optimizer_attempt_rows(
            pi_fit$gp$optimizer_attempts,
            "PI-GP",
            n_calib,
            pi_fit$gp$selected_start
          )
        tryCatch(
          list(
            fit = pi_fit,
            draws = sample_study2_pi_gp(
              fit = pi_fit,
              X_test = test$X,
              C_test = test$C,
              n_draw = n_pred_draw,
              n_measurement_draw = min(200L, n_pred_draw),
              latent_sampler = predictive_latent_sampler,
              n_gibbs = diagnostic_n_new_latent_gibbs,
              rejection_max_batches = rejection_max_batches,
              seed = fit_seed_base + 6000L + n_calib
            )
          ),
          error = function(e) e
        )
      }
      elapsed_pi <- proc.time()[3] - start_pi

      if (inherits(pi_result, "error")) {
        if (is.null(ablation_optimizer_attempts[[paste0("PI_", key)]])) {
          ablation_optimizer_attempts[[paste0("PI_", key)]] <-
            study2_optimizer_attempt_rows(
              pi_result$optimizer_attempts,
              "PI-GP",
              n_calib
            )
        }
        ablation_status[[paste0("PI_", key)]] <-
          study2_ablation_status(
            "PI-GP", n_calib, "failed", conditionMessage(pi_result),
            elapsed_pi
          )
      } else {
        ablation_status[[paste0("PI_", key)]] <-
          study2_ablation_status(
            "PI-GP", n_calib,
            if (measurement_warning) "completed_with_diagnostic_warnings" else "success",
            measurement_advice, elapsed_pi
          )
        ablation_metrics[[paste0("PI_", key)]] <-
          summarize_predictive_samples_by_pattern(
            draw_mat = pi_result$draws,
            y_true = test$y,
            pattern_stratum = pattern_info$pattern_stratum,
            method = "PI-GP",
            rep_id = rep_id,
            n_calib = n_calib,
            scenario = scenario
          )
        if (isTRUE(STUDY2_EVALUATE_F) && n_calib > 0L) {
          ablation_surface[[paste0("PI_", key)]] <-
            study2_latent_gp_surface_metrics(
              gp_fit = pi_result$fit$gp,
              scenario = scenario,
              method = "PI-GP",
              rep_id = rep_id,
              n_calib = n_calib,
              grid_n = if (STUDY2_CONFIG == "quick") 15L else 31L
            )
        }
      }

      if (n_calib == 0L) {
        ablation_status[[paste0("CC_", key)]] <-
          study2_ablation_status(
            "CC-GP", n_calib, "not_applicable", "No complete cases.",
            0
          )
      } else {
        start_cc <- proc.time()[3]
        cc_fit <- tryCatch(
          fit_study2_cc_gp(
            X = train$X,
            y = train$y,
            U = train$U,
            calib_idx = calib_idx,
            measurement_fit = measurement_fit,
            gp_n_starts = STUDY2_ABLATION_GP_N_STARTS,
            gp_seed = fit_seed_base + 7100L + n_calib,
            gp_maxit = STUDY2_ABLATION_GP_MAXIT
          ),
          error = function(e) e
        )
        cc_result <- if (inherits(cc_fit, "error")) {
          cc_fit
        } else {
          ablation_optimizer_attempts[[paste0("CC_", key)]] <-
            study2_optimizer_attempt_rows(
              cc_fit$gp$optimizer_attempts,
              "CC-GP",
              n_calib,
              cc_fit$gp$selected_start
            )
          tryCatch(
            list(
              fit = cc_fit,
              draws = sample_study2_cc_gp(
                fit = cc_fit,
                X_test = test$X,
                C_test = test$C,
                n_draw = n_pred_draw,
                latent_sampler = predictive_latent_sampler,
                n_gibbs = diagnostic_n_new_latent_gibbs,
                rejection_max_batches = rejection_max_batches,
                seed = fit_seed_base + 7000L + n_calib
              )
            ),
            error = function(e) e
          )
        }
        elapsed_cc <- proc.time()[3] - start_cc

        if (inherits(cc_result, "error")) {
          if (is.null(ablation_optimizer_attempts[[paste0("CC_", key)]])) {
            ablation_optimizer_attempts[[paste0("CC_", key)]] <-
              study2_optimizer_attempt_rows(
                cc_result$optimizer_attempts,
                "CC-GP",
                n_calib
              )
          }
          ablation_status[[paste0("CC_", key)]] <-
            study2_ablation_status(
              "CC-GP", n_calib, "failed", conditionMessage(cc_result),
              elapsed_cc
            )
        } else {
          ablation_status[[paste0("CC_", key)]] <-
            study2_ablation_status(
              "CC-GP", n_calib,
              if (measurement_warning) "completed_with_diagnostic_warnings" else "success",
              measurement_advice, elapsed_cc
            )
          ablation_metrics[[paste0("CC_", key)]] <-
            summarize_predictive_samples_by_pattern(
              draw_mat = cc_result$draws,
              y_true = test$y,
              pattern_stratum = pattern_info$pattern_stratum,
              method = "CC-GP",
              rep_id = rep_id,
              n_calib = n_calib,
              scenario = scenario
            )
          if (isTRUE(STUDY2_EVALUATE_F)) {
            ablation_surface[[paste0("CC_", key)]] <-
              study2_latent_gp_surface_metrics(
                gp_fit = cc_result$fit$gp,
                scenario = scenario,
                method = "CC-GP",
                rep_id = rep_id,
                n_calib = n_calib,
                grid_n = if (STUDY2_CONFIG == "quick") 15L else 31L
              )
          }
        }
      }

      rm(measurement_fit, pi_result, pi_fit)
      if (exists("cc_result")) rm(cc_result)
      if (exists("cc_fit")) rm(cc_fit)
      invisible(gc())
    }
  }

  ablation_status_df <- bind_rows(ablation_status)
  if (nrow(ablation_status_df) > 0L) {
    ablation_status_df$rep <- rep_id
    ablation_status_df$scenario <- scenario
  } else {
    ablation_status_df$rep <- integer(0)
    ablation_status_df$scenario <- character(0)
  }
  ablation_optimizer_attempts_df <- bind_rows(ablation_optimizer_attempts)
  if (nrow(ablation_optimizer_attempts_df) > 0L) {
    ablation_optimizer_attempts_df$rep <- rep_id
    ablation_optimizer_attempts_df$scenario <- scenario
  } else {
    ablation_optimizer_attempts_df$rep <- integer(0)
    ablation_optimizer_attempts_df$scenario <- character(0)
  }

  list(
    dataset_identity = data.frame(rep=rep_id, dataset_md5=mixedgp_dataset_identity(frozen$data)),
    metrics = bind_rows(metrics),
    mean_recovery = bind_rows(mean_recovery),
    reference_status = reference_status,
    mean_truth_rejection = mean_truth_rejection,
    mean_truth_diagnostics = mean_truth_diagnostics,
    diagnostics = bind_rows(diagnostics),
    target_diagnostics = bind_rows(target_diagnostics),
    imputation = bind_rows(c(imputation, ablation_imputation)),
    imputation_status = bind_rows(imputation_status),
    surface = bind_rows(surface),
    ablation_metrics = bind_rows(ablation_metrics),
    ablation_surface = bind_rows(ablation_surface),
    measurement_diagnostics = bind_rows(measurement_diagnostics),
    measurement_parameter_diagnostics = bind_rows(
      measurement_parameter_diagnostics
    ),
    ablation_optimizer_attempts = ablation_optimizer_attempts_df,
    ablation_status = ablation_status_df,
    competitor_status = competitor_status,
    pattern_counts = pattern_counts,
    sampler_control_manifest = bind_rows(sampler_control_manifest),
    metadata = list(
      rep = rep_id,
      scenario = scenario,
      n_train = n_train,
      n_test = n_test,
      data_seed_base = data_seed_base,
      fit_seed_base = fit_seed_base,
      common_random_number_group = rep_id,
      calib_grid = calib_grid,
      true_params = dat$true_params,
      frozen_data_file = basename(data_path),
      frozen_data_md5 = attr(frozen, "manifest_md5"),
      frozen_manifest = attr(frozen, "manifest_path"),
      cache_tag = CACHE_TAG,
      scientific_ess_min = STUDY2_ESS_MIN,
      sampler = sampler_controls,
      downstream_monte_carlo = list(
        predictive_draws_requested = n_pred_draw,
        mean_evaluation_points = n_m_eval,
        mean_posterior_draws_requested = n_m_draw,
        mean_latent_integration_draws = n_m_latent,
        mean_truth_latent_draws = n_m_truth,
        mean_truth_method = STUDY2_ORACLE_MEAN_METHOD,
        mean_truth_quadrature_tolerance = if(STUDY2_ORACLE_MEAN_METHOD=="quadrature")1e-5 else NA_real_,
        prospective_latent_sampler = predictive_latent_sampler,
        diagnostic_gibbs_sweeps = diagnostic_n_new_latent_gibbs,
        rejection_max_batches = rejection_max_batches
      )
    )
  )
}

############################################################
## Resumable execution
############################################################

run_grid <- expand.grid(
  scenario = STUDY2_SCENARIOS,
  rep = seq_len(n_rep),
  stringsAsFactors = FALSE
)

if (!is.null(STUDY2_REP_IDS)) {
  ids <- mixedgp_as_integer_strict(STUDY2_REP_IDS,"STUDY2_REP_IDS",1L)
  if (anyDuplicated(ids) || any(ids > n_rep)) stop("Invalid original replication IDs.")
  run_grid <- run_grid[run_grid$rep %in% ids,,drop=FALSE]
}

rep_files <- file.path(
  REP_DIR,
  sprintf(
    "study2_%s_rep%03d_%s.rds",
    run_grid$scenario,
    run_grid$rep,
    CACHE_TAG
  )
)

run_or_load_study2_replication <- function(ii) {
  scenario <- run_grid$scenario[ii]
  rep_id <- run_grid$rep[ii]

  if (
    isTRUE(STUDY2_USE_CACHE) &&
      isTRUE(STUDY2_MC_RESUME) &&
      file.exists(rep_files[ii])
  ) {
    cached <- tryCatch(readRDS(rep_files[ii]), error = function(e) NULL)
    data_path <- file.path(
      STUDY2_DATA_DIR,
      mixedgp_dataset_filename(
        "study2", rep_id, scenario, n_train, n_test, m_vec[1L],
        q = length(m_vec)
      )
    )
    cache_valid <- is.list(cached$metadata) &&
      identical(cached$metadata$rep, as.integer(rep_id)) &&
      identical(cached$metadata$scenario, scenario) &&
      identical(cached$metadata$cache_tag, CACHE_TAG) &&
      identical(cached$metadata$frozen_data_file, basename(data_path)) &&
      file.exists(data_path) &&
      identical(
        tolower(as.character(cached$metadata$frozen_data_md5)),
        tolower(unname(tools::md5sum(data_path)))
      )
    if (isTRUE(cache_valid)) { attr(cached,"eivgp_cache_hit") <- TRUE; return(cached) }
    message("Ignoring incompatible Study II cache: ", rep_files[ii])
  }
  out <- run_one_study2_replication(rep_id, scenario)
  mixedgp_save_replication(out, rep_files[ii])
  out
}
replication_cores <- if (STUDY2_PARALLEL_LEVEL == "hybrid") {
  min(nrow(run_grid), mixedgp_as_integer_strict(STUDY2_DATASET_WORKERS, "dataset workers", 1L, 1L))
} else if (identical(STUDY2_PARALLEL_LEVEL, "replications")) {
  min(nrow(run_grid), mixedgp_resolve_cores())
} else {
  1L
}
rep_objects <- mixedgp_run_replications(
  as.list(seq_len(nrow(run_grid))),
  run_or_load_study2_replication,
  n_cores = replication_cores,
  seeds = 920000L + seq_len(nrow(run_grid)),
  status_path = file.path(TAB_DIR, "study2_replication_status.csv"),
  study = "Study II",
  parallel_map = mixedgp_parallel_lapply,
  mc.preschedule = FALSE
)

dataset_identity <- bind_rows(lapply(rep_objects, `[[`, "dataset_identity"))
write.csv(dataset_identity, file.path(STUDY2_OUT_PREFIX, "dataset_identity.csv"), row.names=FALSE)

mc_results <- bind_rows(lapply(rep_objects, `[[`, "metrics"))
mc_mean_recovery <- bind_rows(lapply(rep_objects, `[[`, "mean_recovery"))
mc_reference_status <- bind_rows(lapply(rep_objects, `[[`, "reference_status"))
mc_mean_truth_rejection <- bind_rows(
  lapply(rep_objects, `[[`, "mean_truth_rejection")
)
mc_mean_truth_diagnostics <- bind_rows(
  lapply(rep_objects, `[[`, "mean_truth_diagnostics")
)
mc_diagnostics <- bind_rows(lapply(rep_objects, `[[`, "diagnostics"))
mc_target_diagnostics <- bind_rows(
  lapply(rep_objects, `[[`, "target_diagnostics")
)
mc_imputation <- bind_rows(lapply(rep_objects, `[[`, "imputation"))
mc_imputation_status <- bind_rows(
  lapply(rep_objects, `[[`, "imputation_status")
)
mc_surface <- bind_rows(lapply(rep_objects, `[[`, "surface"))
mc_ablation_results <- bind_rows(lapply(rep_objects, `[[`, "ablation_metrics"))
mc_ablation_surface <- bind_rows(lapply(rep_objects, `[[`, "ablation_surface"))
measurement_diagnostics <- bind_rows(
  lapply(rep_objects, `[[`, "measurement_diagnostics")
)
measurement_parameter_diagnostics <- bind_rows(
  lapply(rep_objects, `[[`, "measurement_parameter_diagnostics")
)
ablation_optimizer_attempts <- bind_rows(
  lapply(rep_objects, `[[`, "ablation_optimizer_attempts")
)
ablation_status <- bind_rows(lapply(rep_objects, `[[`, "ablation_status"))
competitor_status <- bind_rows(lapply(rep_objects, `[[`, "competitor_status"))
pattern_counts <- bind_rows(lapply(rep_objects, `[[`, "pattern_counts"))
sampler_control_manifest <- bind_rows(
  lapply(rep_objects, `[[`, "sampler_control_manifest")
)
replication_metadata <- lapply(rep_objects, `[[`, "metadata")

raw_outputs <- list(
  predictive_metrics = mc_results,
  reference_status = mc_reference_status,
  mean_recovery = mc_mean_recovery,
  mean_truth_rejection = mc_mean_truth_rejection,
  mean_truth_diagnostics = mc_mean_truth_diagnostics,
  mcmc_diagnostics = mc_diagnostics,
  mcmc_target_diagnostics = mc_target_diagnostics,
  latent_imputation = mc_imputation,
  latent_imputation_status = mc_imputation_status,
  surface_recovery = mc_surface,
  ablation_predictive_metrics = mc_ablation_results,
  ablation_surface_recovery = mc_ablation_surface,
  measurement_diagnostics = measurement_diagnostics,
  measurement_parameter_diagnostics = measurement_parameter_diagnostics,
  ablation_optimizer_attempts = ablation_optimizer_attempts,
  ablation_status = ablation_status,
  competitor_status = competitor_status,
  pattern_counts = pattern_counts,
  sampler_control_manifest = sampler_control_manifest,
  replication_metadata = replication_metadata,
  run_grid = run_grid,
  cache_tag = CACHE_TAG,
  cache_spec = CACHE_SPEC
)
saveRDS(raw_outputs, file.path(RES_DIR, paste0("study2_results_", CACHE_TAG, ".rds")))


## Persist reporting inputs before any table or figure generation.
report_names <- c(
  "ablation_optimizer_attempts", "ablation_status", "CACHE_SPEC", "CACHE_TAG",
  "competitor_status", "d_latent", "m_vec", "mc_ablation_results",
  "mc_ablation_surface", "mc_diagnostics", "mc_imputation", "mc_imputation_status",
  "mc_mean_recovery", "mc_mean_truth_diagnostics", "mc_mean_truth_rejection", "mc_results",
  "mc_surface", "mc_target_diagnostics", "mcmc_gate_rule", "mcmc_pass",
  "mcmc_raw_ess_limit", "mcmc_rhat_limit", "mcmc_target_bulk_ess_limit", "mcmc_target_tail_ess_limit",
  "measurement_burn", "measurement_diagnostics", "measurement_n_chains", "measurement_n_iter",
  "measurement_parameter_diagnostics", "measurement_thin", "method_cols", "method_levels",
  "metric_levels", "n_m_draw", "n_m_eval", "n_m_latent",
  "n_m_truth", "n_rep_eff", "n_test", "n_train",
  "pattern_counts", "predictive_latent_sampler", "raw_outputs", "rejection_max_batches",
  "sampler_control_manifest", "scenario_label", "STUDY2_ABLATION_GP_MAXIT", "STUDY2_ABLATION_GP_N_STARTS",
  "STUDY2_ESS_MIN", "STUDY2_ABLATION_SCENARIOS", "STUDY2_CONTRAST_CALIB", "STUDY2_DESIGN_TAG", "STUDY2_MCMC_RAW_ESS_LIMIT",
  "STUDY2_MCMC_RHAT_LIMIT", "STUDY2_MCMC_TARGET_BULK_ESS_LIMIT", "STUDY2_MCMC_TARGET_DRAWS_PER_CHAIN", "STUDY2_MCMC_TARGET_N_LATENT",
  "STUDY2_MCMC_TARGET_N_POINTS", "STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT", "STUDY2_MEAS_BULK_ESS_LIMIT", "STUDY2_MEAS_RHAT_LIMIT",
  "STUDY2_MEAS_TAIL_ESS_LIMIT", "STUDY2_PUBLISHED_COMPETITORS", "STUDY2_RUN_ABLATIONS", "STUDY2_SCENARIOS")
report_names <- intersect(report_names, ls(envir = environment(), all.names = TRUE))
report_context <- mget(report_names, envir = environment(), inherits = FALSE)
report_context <- Filter(function(x) !is.function(x) && !is.environment(x), report_context)
mixedgp_atomic_save_rds(list(schema = 1L, study = "study2", context = report_context),
        file.path(STUDY2_OUT_PREFIX, "report_inputs.rds"))
mixedgp_write_diagnostic_tables(mc_diagnostics, mc_target_diagnostics,
  TAB_DIR, "study2", overwrite = TRUE)

}


# ---- simulations/competitor_workflow.R ----
# Shared competitor workflow for notebooks and the command-line runner.
# Results are local; manuscript export is a separate caller decision.
mixedgp_run_competitor_workflow <- function(cfg, output, workers = 1L,
    methods = c("UC-GP", "LVGP", "EzGP"), action = c("run", "retry", "export")) {
  action <- match.arg(action)
  cfg <- validate_simulation_config(cfg)
  workers <- mixedgp_validate_scalar_integer(workers, "workers", 1L)
  if (anyDuplicated(methods) || !length(methods) || any(!methods %in% MIXEDGP_PUBLISHED_METHODS)) stop("Invalid competitor methods.")
  study <- cfg$study
  engine <- mixedgp_simulation_engine(cfg$code_dir)
  controls <- engine$mixedgp_competitor_protocol(study)
if(!requireNamespace("knitr",quietly=TRUE)) stop("Install experiment dependencies first (knitr is missing).")
if(action %in% c("run","retry")) {
  preflight <- engine$mixedgp_competitor_preflight(methods,strict=FALSE)
  if(any(!preflight$available)) stop("Install missing comparison packages first: install.packages(c('kergp', 'LVGP', 'EzGP'))")
}
## Validate the complete selected frozen collection before spending time fitting.
for(cell in cfg$cells) mixedgp_verify_cell_data(cfg,cell,engine)
dir.create(output,recursive=TRUE,showWarnings=FALSE)
statuses <- metrics <- predictions <- list()
for(cell in cfg$cells) {
  cat("Competitors:",cell$id,"with",cell$n_rep,"datasets\n")
  work <- function(rep_id) {
    frozen <- mixedgp_read_cell_replication(cfg,cell,rep_id)
    dat <- frozen$data; train <- dat$train; test <- dat$test
    seed <- if(study=="study1") 250000L+1000L*rep_id else
      1000000L*match(cell$scenario,c("primary","latent_additive_control","high_uncertainty","logistic_misspec"))+10000L*rep_id+100L
    X <- if(study=="study1") matrix(train$x,ncol=1L) else train$X
    C <- if(study=="study1") matrix(train$c,ncol=1L) else train$C
    Xt <- if(study=="study1") matrix(test$x,ncol=1L) else test$X
    Ct <- if(study=="study1") matrix(test$c,ncol=1L) else test$C
    cat(cell$id,"replication",rep_id,"starting\n")
    ans <- engine$mixedgp_cached_competitors(X,train$y,C,Xt,Ct,n_draw=500L,seed=seed,
      m_vec=rep(cell$m,ncol(C)),methods=methods,controls=controls,
      allow_fit=action %in% c("run","retry"),retry_failed=action=="retry",
      cache_root=file.path(output,"cache"))
    s <- ans$status; s$cell_id<-cell$id; s$rep<-rep_id
    s$dataset_md5 <- engine$mixedgp_competitor_hash(frozen$data)
    m <- p <- list()
    for(method in names(ans$draws)) {
      m[[method]] <- engine$summarize_predictive_samples_1d(ans$draws[[method]],test$y,
        method=method,rep_id=rep_id,n_calib=NA_integer_,scenario=cell$scenario)
      m[[method]]$cell_id <- cell$id
      p[[method]] <- data.frame(cell_id=cell$id,rep=rep_id,method=method,test_id=seq_along(test$y),
        y=test$y,mean=ans$predictive_means[[method]],variance=ans$predictive_variances[[method]])
    }
    cat(cell$id,"replication",rep_id,"saved; successful methods:",length(m),"\n")
    list(status=s,metrics=mixedgp_bind_rows_base(m),predictions=mixedgp_bind_rows_base(p))
  }
  safe_work <- function(rep_id) tryCatch(work(rep_id), error=function(e)
    structure(conditionMessage(e), class="try-error"))
  answers <- engine$mixedgp_parallel_lapply(as.list(seq_len(cell$n_rep)),safe_work,
    n_cores=workers,seeds=700000L+seq_len(cell$n_rep),mc.preschedule=FALSE)
  for(i in seq_along(answers)) {
    a <- answers[[i]]; id<-paste(cell$id,i)
    if(inherits(a,"try-error") || is.null(a)) {
      statuses[[id]]<-data.frame(cell_id=cell$id,rep=i,method="dataset",status="worker_error",
        message=paste(a,collapse=" ")); next
    }
    statuses[[id]]<-a$status; metrics[[id]]<-a$metrics; predictions[[id]]<-a$predictions
  }
  ## Write incremental tables after each cell. The per-method cache is saved sooner.
  for(item in c("statuses","metrics","predictions"))
    mixedgp_atomic_write_csv(mixedgp_bind_rows_base(get(item)),file.path(output,paste0(item,".csv")))
}
status <- mixedgp_bind_rows_base(statuses); scores <- mixedgp_bind_rows_base(metrics)
summary <- list()
if(nrow(scores)) for(cell_id in unique(scores$cell_id)) for(method in unique(scores$method)) {
  d <- scores[scores$cell_id==cell_id & scores$method==method,,drop=FALSE]
  if(!nrow(d)) next
  for(metric in intersect(c("RMSE","MAE","CRPS","NLPD","Coverage95","Width95","IntervalScore95"),names(d))) {
    v<-d[[metric]]; v<-v[is.finite(v)]
    summary[[paste(cell_id,method,metric)]]<-data.frame(cell_id=cell_id,method=method,metric=metric,
      mean=if(length(v)) mean(v) else NA_real_,mcse=if(length(v)>1L) sd(v)/sqrt(length(v)) else NA_real_,n_success=length(v))
  }
}
summary<-mixedgp_bind_rows_base(summary)
mixedgp_atomic_write_csv(summary,file.path(output,"predictive_summary.csv"))
saveRDS(list(config=cfg,controls=controls,methods=methods,cache=engine$.mixedgp_competitor_cache_root,
  session=sessionInfo()),file.path(output,"provenance.rds"))
writeLines(if(nrow(summary)) knitr::kable(summary,format="latex",booktabs=TRUE,digits=3) else
  "% No successful competitor predictions. Inspect statuses.csv.",
  file.path(output,"predictive_summary.tex"))
if(any(status$status!="success")) warning("Incomplete comparisons: inspect statuses.csv; successful caches are retained.",call.=FALSE)
invisible(list(status=status,summary=summary,output=output))
}


# ---- reporting/combine_experiment_results.R ----
## Reporting only: combine portable numeric outputs, never load/fix model fits.



# ---- reporting/notebook_helpers.R ----
# Keep HTML focused on analysis while preserving complete execution logs.
notebook_logged <- function(path, expr) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  connection <- file(path, open = "wt")
  output_depth <- sink.number()
  sink(connection)
  on.exit({
    while (sink.number() > output_depth) sink()
    close(connection)
  }, add = TRUE)
  warnings <- character()
  value <- withCallingHandlers(eval(substitute(expr), envir = parent.frame()),
    message = function(m) { cat(conditionMessage(m)); invokeRestart("muffleMessage") },
    warning = function(w) {
      warnings <<- unique(c(warnings, conditionMessage(w)))
      cat("WARNING:", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    })
  attr(value, "notebook_warnings") <- warnings
  value
}

notebook_with_env <- function(values, expr) {
  old <- Sys.getenv(names(values), unset = NA_character_)
  on.exit({
    if (any(!is.na(old))) do.call(Sys.setenv, as.list(old[!is.na(old)]))
    if (any(is.na(old))) Sys.unsetenv(names(old)[is.na(old)])
  }, add = TRUE)
  do.call(Sys.setenv, as.list(values))
  eval(substitute(expr), envir = parent.frame())
}

notebook_cache_identity <- function(directory, identity) {
  path <- file.path(directory, "notebook_identity.rds")
  if (file.exists(path)) {
    if (!identical(readRDS(path), identity)) stop("Inputs or settings changed. Use a new output directory: ", directory)
  } else {
    if (dir.exists(directory) && length(list.files(directory, all.files = TRUE, no.. = TRUE)))
      stop("Existing outputs lack notebook provenance. Choose a fresh output directory: ", directory)
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    saveRDS(identity, path)
  }
  invisible(identity)
}


# ---- applications/ADNI_cv.R ----
# Fixed design: no outcome or fitted performance is used to choose the split.
ADNI_CV_FOLDS <- 5L
ADNI_CV_SEED <- 20260914L
adni_make_folds <- function(dat) {
  stopifnot(nrow(dat) == 495L, !anyNA(dat$RID), !anyDuplicated(dat$RID),
            !anyNA(dat$R), all(dat$R %in% 0:1), sum(dat$R == 1L) == 133L)
  old_kind <- RNGkind()
  old_seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({ do.call(RNGkind, as.list(old_kind))
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(ADNI_CV_SEED)
  z <- dat[order(dat$RID), c("RID", "R")]
  z$fold <- NA_integer_; offset <- 0L
  # Continuing the cycle across strata gives 99 people in every fold and
  # either 26 or 27 observed-CSF participants in each held-out group.
  for (r in 0:1) {
    idx <- which(z$R == r)
    idx <- idx[sample.int(length(idx))]
    z$fold[idx] <- (offset + seq_along(idx) - 1L) %% ADNI_CV_FOLDS + 1L
    offset <- offset + length(idx)
  }
  stopifnot(all(table(z$fold) == 99L), all(table(z$fold[z$R == 1L]) %in% 26:27))
  rownames(z) <- NULL
  z
}
adni_read_cohort <- function(data_dir) {
  path <- file.path(data_dir, "toledo_adni_cohort_n495.csv")
  if (unname(tools::md5sum(path)) != "d10b980253304620efea9f2a013166a9")
    stop("The cleaned cohort differs from the frozen ADNI cohort: ", path)
  dat <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("RID", "R", "y_centiloid", "u_csf_abeta42", "x_age_years",
                "x_female", "x_apoe4_dose", "c_ab42_ab40_3", "c_gfap_3", "diagnosis")
  stopifnot(all(required %in% names(dat)))
  dat
}


# ---- applications/ADNI_parallel.R ----
# Fold/chain scheduling. Statistical settings and seed formulas are unchanged.
adni_core_plan <- function(total_cores, n_folds = 5L, n_chains = 4L,
                           max_fold_workers = 5L, os_type = .Platform$OS.type) {
  valid <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 && x == floor(x)
  if (!all(vapply(list(total_cores, n_folds, n_chains, max_fold_workers), valid, logical(1))))
    stop("Core, fold, and chain counts must be positive integers.")
  candidates <- seq_len(min(n_folds, max_fold_workers))
  chains <- pmin(n_chains, floor(total_cores / candidates))
  active <- candidates * chains
  # Maximize concurrent chain workers; ties prefer fewer in-memory folds.
  best <- which.max(active)
  f <- candidates[best]; c <- chains[best]
  if (os_type == "windows") { f <- 1L; c <- 1L }
  list(total_cores = as.integer(total_cores), fold_workers = as.integer(f),
       chain_workers = as.integer(c), active_chain_limit = as.integer(f * c),
       n_folds = as.integer(n_folds), n_chains = as.integer(n_chains))
}

adni_print_core_plan <- function(plan) {
  cat("Total core budget:", plan$total_cores, "| concurrent folds:", plan$fold_workers,
      "| chain workers per fold:", plan$chain_workers,
      "| maximum concurrent chains:", plan$active_chain_limit, "\n")
}

adni_run_jobs <- function(jobs, FUN, plan, status_path) {
  release_budget <- mixedgp_budget_enter(plan$total_cores)
  on.exit(release_budget(),add=TRUE)
  worker <- function(i) {
    fold <- jobs$fold[i]
    dir <- file.path(jobs$output_root[i], paste0("fold_", fold))
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    lock <- file.path(dir, ".fold-lock")
    if (!dir.create(lock, showWarnings = FALSE))
      return(list(ok = FALSE, fold = fold, message = paste("Fold is locked:", lock)))
    on.exit(unlink(lock, recursive = TRUE), add = TRUE)
    writeLines(paste(Sys.info()[["nodename"]], Sys.getpid()), file.path(lock, "owner"))
    log <- file(file.path(dir, "worker.log"), open = "at")
    sink(log); sink(log, type = "message")
    on.exit({ sink(type = "message"); sink(); close(log) }, add = TRUE)
    cat("\nWorker start:", as.character(Sys.time()), "PID", Sys.getpid(), "fold", fold, "\n")
    unlink(file.path(dir, "WORKER_FAILURE.txt"))
    tryCatch({
      value <- FUN(fold)
      cat("Worker finished:", as.character(Sys.time()), "\n")
      list(ok = TRUE, fold = fold, value = value)
    }, error = function(e) {
      text <- paste("Fold", fold, "failed:", conditionMessage(e))
      writeLines(text, file.path(dir, "WORKER_FAILURE.txt"))
      cat(text, "\n")
      list(ok = FALSE, fold = fold, message = text)
    })
  }
  results <- mixedgp_parallel_lapply(seq_len(nrow(jobs)),worker,
    n_cores=plan$fold_workers,.preserve_rng=TRUE)
  good <- vapply(results, function(x) is.list(x) && isTRUE(x$ok), logical(1))
  status <- data.frame(fold = jobs$fold, success = good,
    message = vapply(seq_along(results), function(i) {
      x <- results[[i]]
      if (good[i]) "finished" else if (is.list(x) && !is.null(x$message)) x$message else
        "Worker exited without a result; inspect worker.log and any stale fold lock."
    }, character(1)))
  write.csv(status, status_path, row.names = FALSE)
  if (any(!good)) stop("One or more folds failed. Completed fold checkpoints are retained. See ",
                       status_path, ". No combined report generated.")
  lapply(results, `[[`, "value")
}

# Rscript and direct Rmd rendering use the same fold scheduler.
adni_run_folds <- function(folds, FUN, plan, output_root) {
  jobs <- data.frame(fold = folds, output_root = output_root)
  values <- adni_run_jobs(jobs, function(fold) FUN(fold), plan,
                         file.path(output_root, "fold_worker_status.csv"))
  setNames(values, as.character(folds))
}


# ---- applications/ADNI_case_study_helpers.R ----
# Reporting and comparison helpers. Inference remains in eivGP 0.3.1.
ADNI_CASE_SCHEMA <- "case-study-v5-five-fold"
ADNI_METHODS <- c("EIV-GP", "UC-GP", "LVGP", "EzGP")

adni_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  write.csv(x, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) stop("Could not write ", path)
}

adni_case_preflight <- function(smoke = FALSE) {
  pkgs <- c("ggplot2", "patchwork", "kergp", "LVGP", "EzGP")
  tab <- data.frame(package = pkgs,
    installed = vapply(pkgs, requireNamespace, logical(1), quietly = TRUE))
  if (any(!tab$installed)) stop("Install case-study dependencies via setup_real_data.R: ",
                               paste(tab$package[!tab$installed], collapse = ", "))
  tab$version <- vapply(pkgs, function(p) as.character(utils::packageVersion(p)), character(1))
  tab
}

adni_eda <- function(dat, root, smoke) {
  z <- dat
  z$CSF <- factor(ifelse(z$R == 1, "Observed CSF", "Missing CSF"),
                  levels = c("Observed CSF", "Missing CSF"))
  z$pattern <- interaction(z$c_ab42_ab40_3, z$c_gfap_3, sep = "/", lex.order = TRUE)
  p1 <- ggplot2::ggplot(z, ggplot2::aes(x = CSF, fill = CSF)) +
    ggplot2::geom_bar(width = .65) + ggplot2::labs(x = NULL, y = "Participants", title = "A  CSF availability") +
    ggplot2::theme(legend.position = "none")
  p2 <- ggplot2::ggplot(z[z$R == 1, ], ggplot2::aes(x = u_csf_abeta42, y = y_centiloid)) +
    ggplot2::geom_point(alpha = .55, size = 1.3, color = "#0072B2") +
    ggplot2::labs(x = "Observed CSF A-beta42 (source units)", y = "Amyloid PET (Centiloid)",
                  title = "B  Observed-CSF subset")
  p3 <- ggplot2::ggplot(z, ggplot2::aes(x = pattern, y = y_centiloid)) +
    ggplot2::geom_boxplot(outlier.size = .7, fill = "#B8D8E8") +
    ggplot2::labs(x = "Ordinal A-beta42/40 / GFAP pattern", y = "Amyloid PET (Centiloid)",
                  title = "C  Variation within proxy patterns")
  fig <- patchwork::wrap_plots(p1, p2, p3, nrow = 1) +
    patchwork::plot_annotation(title = if (smoke) "SMOKE RUN - descriptive cohort figure" else
      "ADNI: incomplete CSF measurements and ordinal plasma proxies",
      caption = "Exploratory associations, not causal effects. Proxy classes are constructed from continuous assays.")
  ggplot2::ggsave(file.path(root, "main", "fig1_exploration.pdf"), fig, width = 12, height = 3.8)
  ggplot2::ggsave(file.path(root, "main", "fig1_exploration.png"), fig, width = 12, height = 3.8, dpi = 150)
  counts <- as.data.frame(table(z$pattern, z$CSF))
  names(counts) <- c("pattern", "CSF_availability", "n")
  adni_write_csv(counts, file.path(root, "appendix", "pattern_counts.csv"))
  desc <- do.call(rbind, lapply(split(z, z$CSF), function(a) data.frame(
    group = as.character(a$CSF[1]), n = nrow(a), mean_age = mean(a$x_age_years),
    female_fraction = mean(a$x_female), mean_apoe4_dose = mean(a$x_apoe4_dose),
    mean_centiloid = mean(a$y_centiloid), sd_centiloid = sd(a$y_centiloid))))
  adni_write_csv(desc, file.path(root, "appendix", "cohort_by_CSF_availability.csv"))
}

# Empirical-distribution CRPS uses every draw, in O(S log S) per participant.
adni_draw_summary <- function(draws, truth, point_mean = NULL) {
  stopifnot(is.matrix(draws), ncol(draws) == length(truth), nrow(draws) > 1L,
            all(is.finite(draws)), all(is.finite(truth)))
  s <- nrow(draws)
  cm <- attr(draws, "conditional_means", exact = TRUE)
  cv <- attr(draws, "conditional_vars", exact = TRUE)
  if (is.null(point_mean) && is.matrix(cm) && ncol(cm) == length(truth) && all(is.finite(cm)))
    point_mean <- colMeans(cm)
  mu <- if (is.null(point_mean)) colMeans(draws) else as.numeric(point_mean)
  nlpd <- rep(NA_real_, length(truth))
  if (is.matrix(cm) && is.matrix(cv) && identical(dim(cm), dim(cv)) &&
      ncol(cm) == length(truth) && all(is.finite(cm)) && all(is.finite(cv)) && all(cv > 0)) {
    nlpd <- vapply(seq_along(truth), function(j) {
      ll <- dnorm(truth[j], cm[, j], sqrt(cv[, j]), log = TRUE)
      mx <- max(ll); -(mx + log(mean(exp(ll - mx))))
    }, numeric(1))
  }
  stopifnot(length(mu) == length(truth), all(is.finite(mu)))
  lo <- apply(draws, 2L, stats::quantile, probs = .025)
  hi <- apply(draws, 2L, stats::quantile, probs = .975)
  crps <- vapply(seq_along(truth), function(j) {
    a <- sort(draws[, j])
    mean(abs(a - truth[j])) - sum((2 * seq_len(s) - s - 1) * a) / s^2
  }, numeric(1))
  score <- hi - lo + 40 * pmax(lo - truth, 0) + 40 * pmax(truth - hi, 0)
  data.frame(y_true = truth, pred_mean = mu, pred_sd = apply(draws, 2L, stats::sd),
    pred_lo95 = lo, pred_hi95 = hi, squared_error = (truth - mu)^2,
    absolute_error = abs(truth - mu), CRPS = crps, NLPD = nlpd, IntervalScore95 = score,
    covered95 = truth >= lo & truth <= hi, Width95 = hi - lo,
    covered50 = truth >= apply(draws, 2L, quantile, .25) & truth <= apply(draws, 2L, quantile, .75),
    covered80 = truth >= apply(draws, 2L, quantile, .10) & truth <= apply(draws, 2L, quantile, .90))
}

adni_prediction_rows <- function(draws, test, method, scenario, cv_folds, fold,
                                 convergence = NA, point_mean = NULL) {
  z <- adni_draw_summary(draws, test$y_centiloid, point_mean)
  data.frame(cv_folds = cv_folds, fold = fold,
    RID = test$RID, R = test$R,
    diagnosis = test$diagnosis, method = method, scenario = scenario,
    convergence_passed = convergence, n_draws = nrow(draws), z)
}

adni_metrics <- function(z) {
  if (!nrow(z)) return(data.frame())
  do.call(rbind, lapply(c("overall", "R0", "R1"), function(g) {
    a <- if (g == "overall") z else z[z$R == as.integer(sub("R", "", g)), ]
    data.frame(group = g, n = nrow(a), RMSE = sqrt(mean(a$squared_error)),
      MAE = mean(a$absolute_error), CRPS = mean(a$CRPS),
      NLPD = if (all(is.finite(a$NLPD))) mean(a$NLPD) else NA_real_,
      Coverage95 = mean(a$covered95), Width95 = mean(a$Width95),
      IntervalScore95 = mean(a$IntervalScore95), Coverage50 = mean(a$covered50),
      Coverage80 = mean(a$covered80))
  }))
}

adni_competitors <- function(train, test, n_draw, seed, smoke, path) {
  xc <- c("x_age_years", "x_female", "x_apoe4_dose")
  cc <- c("c_ab42_ab40_3", "c_gfap_3")
  # All transformations use training rows only. Category codes are not scaled.
  center <- colMeans(train[, xc]); scale <- vapply(train[, xc], sd, numeric(1))
  scale[!is.finite(scale) | scale == 0] <- 1
  X <- sweep(sweep(as.matrix(train[, xc]), 2L, center), 2L, scale, "/")
  NX <- sweep(sweep(as.matrix(test[, xc]), 2L, center), 2L, scale, "/")
  yc <- mean(train$y_centiloid); ys <- sd(train$y_centiloid)
  stopifnot(is.finite(ys), ys > 0)
  key <- list(schema = ADNI_CASE_SCHEMA, train = train[, c("RID", xc, cc, "y_centiloid")],
    test = test[, c("RID", xc, cc)], n_draw = n_draw, seed = seed, smoke = smoke,
    eiv_code_md5 = unname(tools::md5sum(list.files(system.file("R", package = "eivGP"), full.names = TRUE))),
    packages = vapply(c("eivGP", "kergp", "LVGP", "EzGP"), function(p) as.character(packageVersion(p)), character(1)))
  if (file.exists(path)) {
    old <- readRDS(path)
    if (identical(old$key, key)) return(old$value)
    stop("Comparator cache configuration changed; remove or relocate: ", path)
  }
  out <- list(draws = list(), means = list(), status = list())
  controls <- if (smoke) list(
    `UC-GP` = list(n_starts = 1L),
    LVGP = list(n_starts = 1L, max_retries = 1L, max_iter_ini = 5L, max_iter_lat = 2L,
                rescue_iter_ini = 5L, rescue_iter_lat = 2L, max_elapsed_seconds = 20),
    EzGP = list(maxeval = 5L, tau_fractions = .01, cv_folds = 2L)
  ) else list()
  for (j in seq_along(ADNI_METHODS[-1])) {
    method <- ADNI_METHODS[-1][j]; started <- proc.time()[3]
    warnings <- character()
    cat("ADNI comparator:", method, "\n"); flush.console()
    ans <- tryCatch(withCallingHandlers(do.call(eivGP::fit_mixedgp_competitor,
      c(list(method = method, X = X, y = (train$y_centiloid - yc) / ys,
        C = as.matrix(train[, cc]), new_X = NX, new_C = as.matrix(test[, cc]),
        m_vec = c(3L, 3L), n_draw = as.integer(n_draw), seed = as.integer(seed + j)), controls[[method]])),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) e)
    good <- !inherits(ans, "error")
    if (good) {
      good <- is.matrix(ans$draws) && identical(dim(ans$draws), c(as.integer(n_draw), as.integer(nrow(test)))) &&
        all(is.finite(ans$draws)) && length(ans$predictive_mean) == nrow(test) && all(is.finite(ans$predictive_mean))
      if (!good) ans <- simpleError("Invalid public-package prediction dimensions or values.")
    }
    if (good) {
      d <- yc + ys * ans$draws
      cm <- attr(ans$draws, "conditional_means", exact = TRUE)
      cv <- attr(ans$draws, "conditional_vars", exact = TRUE)
      attr(d, "conditional_means") <- if (is.null(cm)) NULL else yc + ys * cm
      attr(d, "conditional_vars") <- if (is.null(cv)) NULL else ys^2 * cv
      out$draws[[method]] <- d
      out$means[[method]] <- yc + ys * ans$predictive_mean
    }
    out$status[[method]] <- data.frame(method = method, status = if (good) "success" else "failed",
      optimization_status = if (good && !is.null(ans$optimization_status)) ans$optimization_status else NA_character_,
      elapsed_seconds = unname(proc.time()[3] - started),
      message = if (good) "" else conditionMessage(ans), warnings = paste(unique(warnings), collapse = " | "))
  }
  out$status <- do.call(rbind, out$status)
  save_fit(list(key = key, value = out), path)
  out
}

adni_extended_fold <- function(fit, train, test, draw_ids, fold_dir, cv_folds,
                                fold, smoke, diagnostic_passed, available_draws,
                                seed_base = NULL) {
  out_dir <- file.path(fold_dir, "case-study")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  seed <- if (is.null(seed_base)) as.integer(202680000 + cv_folds) else
    as.integer(seed_base)
  cc <- c("c_ab42_ab40_3", "c_gfap_3")
  xc <- c("x_age_years", "x_female", "x_apoe4_dose")
  Cnew <- adni_make_C(test)
  # Primary target: same test information X,C for every method.
  no_csf <- eivGP::predict_eivgp(fit, new_X = as.matrix(test[, xc]), new_C = Cnew,
    target = "response", draw_ids = draw_ids, seed = seed)
  # Check prediction sensitivity to the retained chain subset; these are not independent CV replications.
  chains <- fit$mcmc$mcmc_draw_info$chain[draw_ids]
  chain_scores <- do.call(rbind, lapply(sort(unique(chains)), function(k) {
    dd <- no_csf[chains == k, , drop = FALSE]
    if (nrow(dd) < 2L) return(NULL)
    cm <- attr(no_csf, "conditional_means", exact = TRUE)
    cv <- attr(no_csf, "conditional_vars", exact = TRUE)
    if (is.matrix(cm) && nrow(cm) == length(chains)) attr(dd, "conditional_means") <- cm[chains == k, , drop = FALSE]
    if (is.matrix(cv) && nrow(cv) == length(chains)) attr(dd, "conditional_vars") <- cv[chains == k, , drop = FALSE]
    zz <- adni_prediction_rows(dd, test, "EIV-GP", "no_test_CSF", cv_folds, fold,
                               diagnostic_passed)
    mm <- adni_metrics(zz); mm$chain <- k; mm
  }))
  if (nrow(chain_scores)) adni_write_csv(chain_scores, file.path(out_dir, "prediction_scores_by_chain.csv"))
  rows <- list(adni_prediction_rows(no_csf, test, "EIV-GP", "no_test_CSF", cv_folds,
    fold, diagnostic_passed),
    adni_prediction_rows(available_draws, test, "EIV-GP", "available_CSF",
      cv_folds, fold, diagnostic_passed))
  comparison <- adni_competitors(train, test, nrow(no_csf), seed + 1000L, smoke,
                                file.path(out_dir, "competitors.rds"))
  for (method in names(comparison$draws)) for (scenario in c("no_test_CSF", "available_CSF")) {
    rows[[length(rows) + 1L]] <- adni_prediction_rows(comparison$draws[[method]], test, method,
      scenario, cv_folds, fold, point_mean = comparison$means[[method]])
  }
  adni_write_csv(do.call(rbind, rows), file.path(out_dir, "predictions.csv"))
  status <- rbind(data.frame(method = "EIV-GP", status = "success",
    optimization_status = if (diagnostic_passed) "mcmc_gate_passed" else "mcmc_gate_failed",
    elapsed_seconds = NA_real_, message = "Sampling time in PROGRESS.md", warnings = ""), comparison$status)
  status$cv_folds <- cv_folds; status$fold <- fold
  adni_write_csv(status, file.path(out_dir, "method_status.csv"))
  # Save only the bounded predictive draw subset, not another full MCMC fit.
  save_fit(list(no_test_CSF = no_csf, available_CSF = available_draws,
    draw_ids = draw_ids, RID = test$RID), file.path(out_dir, "eiv_predictive_draws.rds"))

  observed <- which(test$R == 1L)
  if (length(observed)) {
    # Prospective U|C; neither held-out U nor held-out Y is passed to imputation.
    ud <- eivGP::impute_eivgp(fit, new_C = Cnew[observed, , drop = FALSE],
      draw_ids = draw_ids, scale = "raw", seed = seed + 2000L)
    ud <- matrix(ud[, , 1L], nrow = length(draw_ids), ncol = length(observed))
    scores <- adni_draw_summary(ud, test$u_csf_abeta42[observed])
    scores$probability_nonpositive_CSF <- colMeans(ud <= 0)
    urows <- data.frame(cv_folds = cv_folds, fold = fold,
      RID = test$RID[observed],
      method = "EIV-GP: prospective U|C", convergence_passed = diagnostic_passed, scores)
    # Simple response-free benchmark illustrates that imputation is not exclusive to EIV-GP.
    tr <- train[train$R == 1L, ]
    for (v in cc) { tr[[v]] <- factor(tr[[v]], levels = 1:3) }
    te <- test[observed, ]
    for (v in cc) te[[v]] <- factor(te[[v]], levels = 1:3)
    ub <- tryCatch({
      model <- lm(u_csf_abeta42 ~ c_ab42_ab40_3 + c_gfap_3, data = tr)
      pr <- predict(model, te, se.fit = TRUE)
      set.seed(seed + 2001L)
      v <- matrix(rt(length(draw_ids) * nrow(te), df = model$df.residual), nrow = length(draw_ids))
      v <- sweep(sweep(v, 2, sqrt(pr$se.fit^2 + pr$residual.scale^2), "*"), 2, pr$fit, "+")
      result <- adni_draw_summary(v, test$u_csf_abeta42[observed], pr$fit)
      result$probability_nonpositive_CSF <- colMeans(v <= 0)
      result
    }, error = function(e) e)
    if (!inherits(ub, "error")) urows <- rbind(urows, data.frame(cv_folds = cv_folds,
      fold = fold, RID = test$RID[observed],
      method = "LM-CSF: response-free", convergence_passed = NA, ub))
    writeLines(if (inherits(ub, "error")) conditionMessage(ub) else "success",
               file.path(out_dir, "CSF_benchmark_status.txt"))
    adni_write_csv(urows, file.path(out_dir, "hidden_CSF_validation.csv"))
  }
  # Draws for missing training U are conditional on observed training Y; no truth score exists.
  missing <- which(train$R == 0L)
  ud <- eivGP::impute_eivgp(fit, rows = missing, draw_ids = draw_ids, scale = "raw")
  umat <- matrix(ud[, , 1L], nrow = length(draw_ids), ncol = length(missing))
  adni_write_csv(data.frame(RID = train$RID[missing],
    mean = colMeans(umat), probability_nonpositive_CSF = colMeans(umat <= 0), lo95 = apply(umat, 2, quantile, .025),
    hi95 = apply(umat, 2, quantile, .975)), file.path(out_dir, "training_missing_CSF_posterior.csv"))

  # Fixed, prespecified fold 1 for a physical-coordinate surface, never best-fold selection.
  if (fold == 1L) {
    ug <- seq(quantile(train$u_csf_abeta42[train$R == 1], .1),
              quantile(train$u_csf_abeta42[train$R == 1], .9), length.out = if (smoke) 8L else 40L)
    reference <- c(x_age_years = median(train$x_age_years), x_female = 0, x_apoe4_dose = 1)
    xx <- matrix(rep(reference, each = length(ug)), ncol = 3,
                 dimnames = list(NULL, names(reference)))
    fs <- eivGP::predict_eivgp(fit, new_X = xx, new_U = matrix(ug, ncol = 1),
      new_U_scale = "raw", target = "surface", draw_ids = draw_ids, seed = seed + 3000L)
    adni_write_csv(data.frame(U = ug, mean = colMeans(fs), lo95 = apply(fs, 2, quantile, .025),
      hi95 = apply(fs, 2, quantile, .975), age = unname(reference[1]), female = 0, apoe4_dose = 1,
      cv_folds = cv_folds, fold = fold,
      convergence_passed = diagnostic_passed),
      file.path(out_dir, "CSF_PET_surface.csv"))
  }
  # Include the dictionary states deliberately omitted from the lightweight stopping gate.
  adni_write_csv(fit$diagnostics$rhat_hyper, file.path(out_dir, "all_hyper_diagnostics.csv"))
  j <- fit$mcmc$samples_by_chain$J
  if (length(j)) {
    occupancy <- do.call(rbind, lapply(seq_along(j), function(k) do.call(rbind,
      lapply(seq_len(ncol(j[[k]])), function(d) {
        tb <- table(j[[k]][, d]); data.frame(chain = k, coordinate = d,
          index = as.integer(names(tb)), fraction = as.numeric(tb) / sum(tb))
      }))))
    adni_write_csv(occupancy, file.path(out_dir, "dictionary_occupancy.csv"))
  }
  invisible(NULL)
}



# Public package workflows. Numerical settings remain explicit and results external.
eivgp_module_identity <- function() {
  p <- system.file("module-manifest.csv", package = "eivGP")
  if (!nzchar(p)) stop("Package module manifest is missing; reinstall eivGP.")
  z <- utils::read.csv(p, stringsAsFactors = FALSE)
  setNames(z$md5, z$module)
}

study_config <- function(study = c("study1", "study2"), profile = c("quick", "development", "paper"),
                         output_root, data_root, cores = 4L) {
  study <- match.arg(study); profile <- match.arg(profile)
  mode <- c(quick="smoke",development="development",paper="publication")[[profile]]
  constructor <- get(paste0(study,"_simulation_config"))
  constructor(mode, code_dir=system.file(package="eivGP"), output_root=output_root,
              data_root=data_root, core_budget=cores)
}

run_numerical_study <- function(config, action=c("plan","data","fit","run","report"),
                                competitors=c("UC-GP","LVGP","EzGP"), run_dir=NULL) {
  release_budget <- mixedgp_budget_enter(if(is.null(config$parallel$core_budget)) config$parallel$workers else config$parallel$core_budget)
  on.exit(release_budget(),add=TRUE)
  action<-match.arg(action)
  if(action=="plan") {
    config$mode<-"dry_run";config$stages<-character()
    return(mixedgp_run_simulation(config))
  }
  if(action=="report") {
    if(is.null(run_dir)) run_dir<-mixedgp_run_directory(config)
    status<-mixedgp_report_run(run_dir,"report",system.file(package="eivGP"))
    if(any(status$status!="success")) stop("Reporting failed; inspect status outputs.")
    return(invisible(list(run_dir=run_dir,status=status)))
  }
  old_checkpoint <- options(eivGP.checkpoint_root=if(identical(config$study,"study1"))
    file.path(mixedgp_run_directory(config),"checkpoints") else getOption("eivGP.checkpoint_root",
    file.path(mixedgp_run_directory(config),"checkpoints")))
  on.exit(options(old_checkpoint),add=TRUE)
  config$stages<-switch(action,data="data",fit=c("data","fit"),run=c("data","fit","aggregate"))
  result<-mixedgp_run_simulation(config)
  if(action=="run" && length(competitors)) {
    comp_dir<-file.path(result$run_dir,"competitors")
    comp<-mixedgp_run_competitor_workflow(config,comp_dir,workers=config$parallel$core_budget,methods=competitors)
    result$comparison<-mixedgp_combine_saved_runs(result$run_dir,comp_dir,file.path(result$run_dir,"reporting","comparison"))
  }
  invisible(result)
}

read_study_results <- function(path, competitor_dir = NULL) {
  # A selected merged export is already complete; never append original rows.
  read <- function(name) {p<-file.path(path,paste0(name,".csv"));if(file.exists(p)) utils::read.csv(p,check.names=FALSE) else data.frame()}
  competitors <- if(is.null(competitor_dir)) read("competitor_metrics") else
    utils::read.csv(file.path(competitor_dir,"metrics.csv"),check.names=FALSE)
  list(prediction=read("predictive_metrics"),competitors=competitors,
       ablations=read("ablation_predictive_metrics"),diagnostics=read("mcmc_diagnostics"),
       mean=read("mean_recovery"),surface=read("surface_recovery"),imputation=read("latent_imputation"))
}


read_fit <- function(path) mixedgp_read_artifact(path)
save_fit <- function(object,path) mixedgp_save_artifact(object,path)

eivgp_adni_driver <- function(data_dir,output_dir) {
# ADNI implementation extracted from the former analysis notebook.
# Called by adni_run_analysis(); reader-facing ADNI work is deferred.


## Environment variables override the notebook defaults without editing it:
## EIVGP_CORES=12 or EIVGP_SMOKE_TEST=1. Five folds and four chains are fixed.
env_int <- function(name, default, minimum = 0L) {
  raw <- Sys.getenv(name, unset = as.character(default))
  if (!grepl("^[0-9]+$", raw)) stop(name, " must be an integer.")
  x <- suppressWarnings(as.integer(raw))
  if (length(x) != 1L || is.na(x) || x < minimum) stop(name, " is invalid.")
  x
}
env_flag <- function(name, default = FALSE) {
  x <- tolower(Sys.getenv(name, unset = if (default) "true" else "false"))
  if (!x %in% c("true", "false", "1", "0", "yes", "no")) stop(name, " is invalid.")
  x %in% c("true", "1", "yes")
}

PREPARE_ONLY <- isTRUE(get0("ADNI_PREPARE_ONLY", inherits = FALSE))
CV_FOLDS <- 5L
TOTAL_CORES <- env_int("EIVGP_CORES", env_int("EIVGP_N_CORES", 4L, 1L), 1L)
release_budget <- mixedgp_budget_enter(TOTAL_CORES)
on.exit(release_budget(),add=TRUE)
MAX_FOLD_WORKERS <- env_int("ADNI_MAX_FOLD_WORKERS", 5L, 1L)
SMOKE_TEST <- env_flag("EIVGP_SMOKE_TEST", FALSE)


if (SMOKE_TEST) {
  N_CHAINS <- 4L; BURN <- 2L; INITIAL_ITER <- 8L
  EXTENSION_ITER <- 4L; CHECKPOINT_EVERY <- 3L
  FOLDS_TO_RUN <- if (env_flag("ADNI_SMOKE_ALL_FOLDS", FALSE)) 1:5 else 1L
  PREDICTION_DRAWS <- 8L
} else {
  N_CHAINS <- 4L; BURN <- 4000L; INITIAL_ITER <- 20000L
  EXTENSION_ITER <- 10000L; CHECKPOINT_EVERY <- 2000L
  FOLDS_TO_RUN <- 1:5
  PREDICTION_DRAWS <- 1000L
}

## Frozen sampler choice. A matched ADNI pilot found that block sizes 16 and
## 32 reduced wall time per transition but worsened missing-U mixing enough
## that they did not improve convergence per equal waiting time. Keep the
## package reference schedule for the confirmatory run.
SAMPLER_CONTROL <- list(
  u_block_size = 8L,
  cross_every = 10L,
  ess_max_steps = 10000L,
  dictionary_mode = "conditional",
  max_dictionary_size = 125L,
  gp_block_schur = TRUE,
  threshold_update = "gibbs"
)

R_HAT_LIMIT <- 1.01
BULK_ESS_LIMIT <- 400
TAIL_ESS_LIMIT <- 400
MCSE_SD_RATIO_LIMIT <- 0.10

PROJECT_DIR <- system.file(package="eivGP")
DATA_DIR <- data_dir
CORE_PLAN <- if (PREPARE_ONLY) get("ADNI_RUN_CORE_PLAN", inherits = FALSE) else
  adni_core_plan(TOTAL_CORES, length(FOLDS_TO_RUN), N_CHAINS, MAX_FOLD_WORKERS)
N_CORES <- CORE_PLAN$chain_workers
adni_print_core_plan(CORE_PLAN)
DATA_FILE <- file.path(DATA_DIR, "toledo_adni_cohort_n495.csv")
OUTPUT_ROOT <- file.path(output_dir, if (SMOKE_TEST) "smoke" else "five_fold")
FOLD_FILE <- file.path(OUTPUT_ROOT,"appendix","fold_assignments.csv")

dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)
RUN_LOCK <- file.path(OUTPUT_ROOT, ".run-lock")
LOCK_OWNED <- FALSE
if (!dir.create(RUN_LOCK, showWarnings = FALSE)) stop("Run already locked: ", RUN_LOCK,
  ". If a job was killed, verify it has stopped before removing this directory.")
LOCK_OWNED <- TRUE
on.exit(if (LOCK_OWNED && dir.exists(RUN_LOCK)) unlink(RUN_LOCK,recursive=TRUE),add=TRUE)
writeLines(paste(Sys.info()[["nodename"]], Sys.getpid()), file.path(RUN_LOCK, "owner"))


## Dependency installation is a separate, explicit setup step.
EIVGP_TESTED_COMMIT <- paste(eivgp_module_identity(),collapse=";")
if (!requireNamespace("eivGP", quietly = TRUE)) stop("Run setup_real_data.R before analysis.")
pkg <- utils::packageDescription("eivGP")




for (d in c("main", "appendix")) dir.create(file.path(OUTPUT_ROOT, d), showWarnings = FALSE)
case_packages <- adni_case_preflight(SMOKE_TEST)
adni_write_csv(as.data.frame(CORE_PLAN), file.path(OUTPUT_ROOT, "appendix", "resource_plan.csv"))
package_files <- list.files(system.file("R", package = "eivGP"), full.names = TRUE)
adni_write_csv(data.frame(file = basename(package_files), md5 = unname(tools::md5sum(package_files))),
  file.path(OUTPUT_ROOT, "appendix", "eivGP_code_fingerprints.csv"))
adni_write_csv(case_packages, file.path(OUTPUT_ROOT, "appendix", "package_versions.csv"))
adni_write_csv(data.frame(schema = ADNI_CASE_SCHEMA, smoke = SMOKE_TEST,
  prediction_draws = PREDICTION_DRAWS, run_mode = "five_fold", cv_folds = CV_FOLDS, split_seed = ADNI_CV_SEED,
  primary_target = "Y | test X,C (no test CSF)",
  competitors = paste(ADNI_METHODS, collapse = ";")),
  file.path(OUTPUT_ROOT, "appendix", "case_config.csv"))

cat("Project:", PROJECT_DIR, "\n")
cat("Data:", DATA_DIR, "\n")
cat("Five-fold CV | folds:", paste(FOLDS_TO_RUN, collapse = ","), "\n")
cat("Mode:", if (SMOKE_TEST) "SMOKE" else "PRODUCTION", "\n")
cat("eivGP:", as.character(utils::packageVersion("eivGP")), "\n")

dat <- adni_read_cohort(DATA_DIR)
folds <- adni_make_folds(dat)
if (file.exists(FOLD_FILE)) {
  previous <- read.csv(FOLD_FILE)
  if (!isTRUE(all.equal(previous, folds, check.attributes = FALSE)))
    stop("Saved fold assignment differs. Use a separate output directory.")
} else adni_write_csv(folds, FOLD_FILE)
adni_write_csv(as.data.frame(table(folds$fold, folds$R)),
               file.path(OUTPUT_ROOT, "appendix", "fold_CSF_counts.csv"))

X_COLUMNS <- c("x_age_years", "x_female", "x_apoe4_dose")
C_COLUMNS <- c("c_ab42_ab40_3", "c_gfap_3")
Y_COLUMN <- "y_centiloid"
U_COLUMN <- "u_csf_abeta42"


adni_eda(dat, OUTPUT_ROOT, SMOKE_TEST)



atomic_write_lines <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(paste0(basename(path), "-"), tmpdir = dirname(path))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  writeLines(lines, tmp)
  if (!file.rename(tmp, path)) stop("Could not publish ", path)
}

stamp <- function(run_dir, message) {
  line <- sprintf("- [%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), message)
  cat(line, "\n", file = file.path(run_dir, "PROGRESS.md"), append = TRUE)
  cat(line, "\n")
  flush.console()
}

set_status <- function(run_dir, fold, stage, detail) {
  atomic_write_lines(c(
    paste0("updated_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("cv_folds=", CV_FOLDS), paste0("fit_unit=", fold),
    paste0("stage=", stage), paste0("detail=", detail)
  ), file.path(run_dir, "CURRENT_STATUS.txt"))
  stamp(run_dir, paste("Five-fold CV", "|", stage, "|", detail))
}

light_diagnostics <- function(fit, fold_dir, stage) {
  common <- c("parameter", "rhat", "ess_bulk", "ess_tail", "mcse_sd_ratio")
  parts <- list(
    hyper = subset(fit$diagnostics$rhat_hyper,
                   !grepl("^(J\\[|dictionary_occupancy)", parameter))[, common],
    loading = fit$diagnostics$rhat_A[, common],
    threshold = fit$diagnostics$rhat_tau[, common],
    missing_U = fit$diagnostics$rhat_U[, common]
  )
  tab <- do.call(rbind, lapply(names(parts), function(g) {
    z <- parts[[g]]; z$group <- g; z
  }))
  rownames(tab) <- NULL
  finite <- with(tab, is.finite(rhat) & is.finite(ess_bulk) &
                   is.finite(ess_tail) & is.finite(mcse_sd_ratio))
  passed <- all(finite) && all(tab$rhat <= R_HAT_LIMIT) &&
    all(tab$ess_bulk >= BULK_ESS_LIMIT) && all(tab$ess_tail >= TAIL_ESS_LIMIT) &&
    all(tab$mcse_sd_ratio <= MCSE_SD_RATIO_LIMIT)
  summary <- data.frame(
    stage = stage,
    n_iter = fit$diagnostics$summary$n_iter,
    saved_per_chain = fit$diagnostics$summary$saved_per_chain,
    max_rhat = if (all(finite)) max(tab$rhat) else NA_real_,
    min_bulk_ess = if (all(finite)) min(tab$ess_bulk) else NA_real_,
    min_tail_ess = if (all(finite)) min(tab$ess_tail) else NA_real_,
    max_mcse_sd_ratio = if (all(finite)) max(tab$mcse_sd_ratio) else NA_real_,
    light_gate_passed = passed,
    stringsAsFactors = FALSE
  )
  write.csv(tab, file.path(fold_dir, paste0("diagnostics_", stage, "_parameters.csv")),
            row.names = FALSE)
  write.csv(summary, file.path(fold_dir, paste0("diagnostics_", stage, "_summary.csv")),
            row.names = FALSE)
  list(table = tab, summary = summary, passed = passed)
}

trace_series <- function(fit) {
  out <- list(V = fit$mcmc$samples_by_chain$V, r = fit$mcmc$samples_by_chain$r)
  H <- fit$mcmc$samples_by_chain$logtheta
  for (j in seq_len(ncol(H[[1L]]))) {
    nm <- colnames(H[[1L]])[j]
    out[[nm]] <- lapply(H, function(x) x[, j])
  }
  A <- fit$mcmc$samples_by_chain$A
  for (j in seq_len(dim(A[[1L]])[2L]))
    out[[paste0("A[", j, ",1]")]] <- lapply(A, function(x) x[, j, 1L])
  T <- fit$mcmc$samples_by_chain$tau
  for (j in seq_len(ncol(T[[1L]])))
    out[[colnames(T[[1L]])[j]]] <- lapply(T, function(x) x[, j])
  out
}

draw_trace_panel <- function(series, path, title) {
  png(path, width = 2400, height = 1800, res = 180)
  on.exit(dev.off(), add = TRUE)
  old <- par(mfrow = c(4, 4), mar = c(2.4, 2.7, 2.1, 0.7), oma = c(0, 0, 2, 0))
  on.exit(par(old), add = TRUE)
  cols <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
  for (nm in names(series)) {
    chains <- series[[nm]]
    n <- length(chains[[1L]])
    ids <- unique(as.integer(round(seq(1, n, length.out = min(n, 2000L)))))
    ylim <- range(unlist(lapply(chains, function(x) x[ids])), finite = TRUE)
    plot(ids, chains[[1L]][ids], type = "l", col = cols[1L], xlab = "retained draw",
         ylab = "", main = nm, ylim = ylim)
    if (length(chains) > 1L) for (cc in 2:length(chains))
      lines(ids, chains[[cc]][ids], col = cols[cc])
  }
  mtext(title, outer = TRUE, cex = 1.2)
}

draw_convergence_plots <- function(fit, diagnostic, fold_dir, stage) {
  draw_trace_panel(trace_series(fit),
                   file.path(fold_dir, paste0("trace_key_", stage, ".png")),
                   paste("Five-fold CV", stage, "key traces"))

  udiag <- fit$diagnostics$rhat_U
  worst <- head(udiag$global_index[order(udiag$rhat, decreasing = TRUE)], 6L)
  U <- fit$mcmc$samples_by_chain$U
  useries <- setNames(lapply(worst, function(j) lapply(U, function(x) x[, j, 1L])),
                      paste0("U[", worst, "]"))
  draw_trace_panel(useries,
                   file.path(fold_dir, paste0("trace_worst_U_", stage, ".png")),
                   paste("Worst missing-U traces by raw R-hat:", stage))

  tab <- diagnostic$table
  png(file.path(fold_dir, paste0("rhat_overview_", stage, ".png")),
      width = 1800, height = 1100, res = 180)
  cols <- c(hyper = "#0072B2", loading = "#D55E00",
            threshold = "#009E73", missing_U = "#777777")
  plot(seq_len(nrow(tab)), tab$rhat, pch = 16, col = cols[tab$group],
       xlab = "raw parameter (grouped)", ylab = "rank-normalized split R-hat",
       main = paste("Five-fold CV", stage),
       ylim = range(c(1, tab$rhat, R_HAT_LIMIT), finite = TRUE))
  abline(h = R_HAT_LIMIT, col = "red", lty = 2, lwd = 2)
  legend("topright", legend = names(cols), col = cols, pch = 16, bty = "n")
  dev.off()
}

balanced_draw_ids <- function(fit, total = 500L) {
  info <- fit$mcmc$mcmc_draw_info
  chains <- sort(unique(info$chain))
  target <- max(1L, floor(total / length(chains)))
  sort(unique(unlist(lapply(chains, function(cc) {
    a <- which(info$chain == cc)
    a[unique(as.integer(round(seq(1, length(a), length.out = min(target, length(a))))))]
  }), use.names = FALSE)))
}

score_predictions <- function(z) {
  split_z <- c(list(overall = z), split(z, paste0("R", z$R)))
  do.call(rbind, lapply(names(split_z), function(g) {
    x <- split_z[[g]]
    data.frame(
      group = g, n = nrow(x),
      RMSE = sqrt(mean((x$y_true - x$pred_mean)^2)),
      MAE = mean(abs(x$y_true - x$pred_mean)),
      Coverage95 = mean(x$y_true >= x$pred_lo95 & x$y_true <= x$pred_hi95),
      Width95 = mean(x$pred_hi95 - x$pred_lo95),
      stringsAsFactors = FALSE
    )
  }))
}

run_one_fold <- function(fold) {
  fold_dir <- file.path(OUTPUT_ROOT, paste0("fold_", fold))
  dir.create(fold_dir, recursive = TRUE, showWarnings = FALSE)
  unlink(file.path(fold_dir, "CASE_STUDY_COMPLETED.txt"))
  set_status(fold_dir, fold, "prepare", "Preparing frozen train/test data.")

  test_ids <- folds$RID[folds$fold == fold]
  train <- dat[!(dat$RID %in% test_ids), , drop = FALSE]
  test <- dat[dat$RID %in% test_ids, , drop = FALSE]
  train <- train[order(train$RID), , drop = FALSE]
  test <- test[order(test$RID), , drop = FALSE]
  stopifnot(nrow(train) == 396L, nrow(test) == 99L)

  fit_seed <- 202660000L + fold
  prediction_seed <- 202670000L + fold
  case_seed <- 202680000L + fold

  X_train <- as.matrix(train[, X_COLUMNS, drop = FALSE])
  C_train <- adni_make_C(train)
  U_train <- matrix(train[[U_COLUMN]], ncol = 1L,
                    dimnames = list(NULL, U_COLUMN))
  U_train[train$R == 0L, 1L] <- NA_real_

  config <- data.frame(
    run_mode = "five_fold", cv_folds = CV_FOLDS, fold = fold,
    assignment_source = "fixed_R_stratified_5fold", assignment_seed = ADNI_CV_SEED,
    data_md5 = unname(tools::md5sum(DATA_FILE)),
    folds_md5 = unname(tools::md5sum(FOLD_FILE)),
    n_train = nrow(train), n_test = nrow(test),
    n_train_observed_U = sum(train$R == 1L),
    n_train_missing_U = sum(train$R == 0L),
    engine = "multivariate", latent_dim = 1L, ident = "none",
    kernel = "se", standardize_U = TRUE,
    n_chains = N_CHAINS, burn = BURN, initial_iter = INITIAL_ITER,
    extension_iter = EXTENSION_ITER, checkpoint_every = CHECKPOINT_EVERY,
    u_block_size = SAMPLER_CONTROL$u_block_size,
    cross_every = SAMPLER_CONTROL$cross_every,
    gp_block_schur = SAMPLER_CONTROL$gp_block_schur,
    dictionary_mode = SAMPLER_CONTROL$dictionary_mode,
    max_dictionary_size = SAMPLER_CONTROL$max_dictionary_size,
    ess_max_steps = SAMPLER_CONTROL$ess_max_steps,
    threshold_update = SAMPLER_CONTROL$threshold_update,
    fit_seed = fit_seed, prediction_seed = prediction_seed, case_seed = case_seed,
    eivgp_version = as.character(utils::packageVersion("eivGP")),
    tested_commit = EIVGP_TESTED_COMMIT,
    installed_commit = if (is.null(pkg$RemoteSha)) "unrecorded" else pkg$RemoteSha,
    stringsAsFactors = FALSE
  )
  config_path <- file.path(fold_dir, "run_config.csv")
  if (file.exists(file.path(fold_dir, "fit_latest.rds"))) {
    if (!file.exists(config_path)) stop("Checkpoint has no run configuration.")
    old_config <- read.csv(config_path, stringsAsFactors = FALSE)
    if (!isTRUE(all.equal(old_config, config, check.attributes = FALSE)))
      stop("Checkpoint configuration mismatch; use a separate output directory.")
  }
  write.csv(config, config_path, row.names = FALSE)

  fit_path <- file.path(fold_dir, "fit_latest.rds")
  if (file.exists(fit_path)) {
    fit <- mixedgp_read_artifact(fit_path)
    stopifnot(
      identical(fit$sampler_version, "0.3.1"),
      identical(fit$control$u_block_size, SAMPLER_CONTROL$u_block_size),
      identical(fit$control$cross_every, SAMPLER_CONTROL$cross_every),
      fit$diagnostics$summary$burn == BURN
    )
    stamp(fold_dir, paste("Loaded checkpoint at iteration", fit$diagnostics$summary$n_iter))
  } else {
    first_target <- min(INITIAL_ITER, BURN + CHECKPOINT_EVERY)
    set_status(fold_dir, fold, "fit", paste("Initial fit to iteration", first_target))
    started <- Sys.time()
    fit <- fit_eivgp(
      X = X_train, y = train[[Y_COLUMN]], C = C_train, U_obs = U_train,
      engine = "multivariate", latent_dim = 1L, ident = "none",
      standardize_U = TRUE, store_scores = FALSE, kernel = "se",
      parallel = TRUE, n_cores = N_CORES, n_chains = N_CHAINS,
      n_iter = first_target, burn = BURN, thin = 1L,
      seed = fit_seed,
      verbose = FALSE, sampler_control = SAMPLER_CONTROL
    )
    save_fit(fit, fit_path)
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    stamp(fold_dir, sprintf("Saved checkpoint at %d; segment %.1f min.",
                            first_target, elapsed / 60))
  }

  extend_to <- function(fit, target) {
    while (fit$diagnostics$summary$n_iter < target) {
      current <- fit$diagnostics$summary$n_iter
      add <- min(CHECKPOINT_EVERY, target - current)
      set_status(fold_dir, fold, "continue",
                 sprintf("Iteration %d -> %d; checkpoint will overwrite fit_latest.rds.",
                         current, current + add))
      started <- Sys.time()
      fit <- continue_eivgp(fit, n_iter = add, parallel = TRUE,
                            n_cores = N_CORES, verbose = FALSE)
      save_fit(fit, fit_path)
      elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
      remaining <- target - fit$diagnostics$summary$n_iter
      eta <- if (remaining > 0L) elapsed / add * remaining / 60 else 0
      stamp(fold_dir, sprintf(
        "Saved iteration %d; last %d took %.1f min; ETA to %d about %.1f min.",
        fit$diagnostics$summary$n_iter, add, elapsed / 60, target, eta
      ))
    }
    fit
  }

  fit <- extend_to(fit, INITIAL_ITER)
  current_iter <- fit$diagnostics$summary$n_iter
  final_stage <- paste0("iter", current_iter)
  diagnostic <- light_diagnostics(fit, fold_dir, final_stage)
  draw_convergence_plots(fit, diagnostic, fold_dir, final_stage)
  if (!diagnostic$passed && current_iter < INITIAL_ITER + EXTENSION_ITER) {
    fit <- extend_to(fit, INITIAL_ITER + EXTENSION_ITER)
    final_stage <- paste0("iter", fit$diagnostics$summary$n_iter)
    diagnostic <- light_diagnostics(fit, fold_dir, final_stage)
    draw_convergence_plots(fit, diagnostic, fold_dir, final_stage)
  }
  stamp(fold_dir, paste(final_stage, "convergence gate:", diagnostic$passed))

  set_status(fold_dir, fold, "prediction", paste("Predicting from", final_stage, "fit."))
  X_test <- as.matrix(test[, X_COLUMNS, drop = FALSE])
  C_test <- adni_make_C(test)
  U_test <- matrix(test[[U_COLUMN]], ncol = 1L,
                   dimnames = list(NULL, U_COLUMN))
  U_test[test$R == 0L, 1L] <- NA_real_
  draw_ids <- balanced_draw_ids(fit, PREDICTION_DRAWS)
  response_draws <- predict_eivgp(
    fit, new_X = X_test, new_C = C_test, new_U = U_test,
    new_U_scale = "raw", target = "response", draw_ids = draw_ids,
    seed = prediction_seed
  )
  pred <- data.frame(
    run_mode = "five_fold",
    cv_folds = CV_FOLDS,
    fold = fold, final_stage = final_stage,
    convergence_passed = diagnostic$passed,
    RID = test$RID, diagnosis = test$diagnosis, R = test$R,
    y_true = test[[Y_COLUMN]], pred_mean = colMeans(response_draws),
    pred_sd = apply(response_draws, 2L, stats::sd),
    pred_lo95 = apply(response_draws, 2L, stats::quantile, probs = 0.025),
    pred_hi95 = apply(response_draws, 2L, stats::quantile, probs = 0.975),
    n_prediction_draws = nrow(response_draws), stringsAsFactors = FALSE
  )
  write.csv(pred, file.path(fold_dir, "eivgp_predictions.csv"), row.names = FALSE)
  metrics <- score_predictions(pred)
  metrics$convergence_passed <- diagnostic$passed
  metrics$run_mode <- "five_fold"
  metrics$cv_folds <- CV_FOLDS
  metrics$fold <- fold; metrics$final_stage <- final_stage
  write.csv(metrics, file.path(fold_dir, "eivgp_metrics.csv"), row.names = FALSE)
  atomic_write_lines(c(
    paste0("completed_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("cv_folds=", CV_FOLDS),
    paste0("fit_unit=", fold),
    paste0("convergence_passed=", diagnostic$passed),
    paste0("final_stage=", final_stage)
  ), file.path(fold_dir, "EIV_COMPLETED.txt"))
  set_status(fold_dir, fold, "EIV_complete", paste("EIV prediction completed at", final_stage))
  cat("\n### Fold", fold, "\n\n")
  print(knitr::kable(metrics, digits = 3))
  set_status(fold_dir, fold, "comparison", "Running shared prediction and CSF inference tasks.")
  adni_extended_fold(fit, train, test, draw_ids, fold_dir, CV_FOLDS, fold,
                     SMOKE_TEST, diagnostic$passed, response_draws,
                     seed_base = case_seed)
  atomic_write_lines(c(paste0("completed_at=", Sys.time()),
    paste0("schema=", ADNI_CASE_SCHEMA), paste0("convergence_passed=", diagnostic$passed)),
    file.path(fold_dir, "CASE_STUDY_COMPLETED.txt"))
  set_status(fold_dir, fold, "complete", "EIV fit, comparison, and additional inference finished; inspect method status.")
  rm(fit, response_draws)
  invisible(gc())
  metrics
}

if (!PREPARE_ONLY) {
cat("Running fit units; progress and errors are written to fold_N/worker.log.\n")
all_metrics <- adni_run_folds(FOLDS_TO_RUN, run_one_fold, CORE_PLAN, OUTPUT_ROOT)
combined_metrics <- do.call(rbind, all_metrics)
write.csv(combined_metrics, file.path(OUTPUT_ROOT, "all_fold_metrics.csv"), row.names = FALSE)
}

if (!PREPARE_ONLY) {
  main_table <- adni_case_report(OUTPUT_ROOT, folds_expected = FOLDS_TO_RUN,
                                 smoke = SMOKE_TEST)

}

writeLines(capture.output(sessionInfo()), file.path(OUTPUT_ROOT, "sessionInfo.txt"))
sessionInfo()
if (!PREPARE_ONLY) unlink(RUN_LOCK, recursive = TRUE)

}

run_adni_analysis <- function(data_dir,output_dir,action=c("plan","run","report"),
                              profile=c("paper","quick"),cores=4L) {
  action<-match.arg(action);profile<-match.arg(profile)
  folds<-if(profile=="paper")1:5 else 1L
  if(action=="plan") return(data.frame(profile=profile,fold=folds,action="planned",
    data_available=file.exists(file.path(data_dir,"toledo_adni_cohort_n495.csv"))))
  target<-file.path(output_dir,if(profile=="paper")"five_fold" else "smoke")
  if(action=="report") return(adni_case_report(target,folds_expected=folds,smoke=profile=="quick"))
  notebook_with_env(c(EIVGP_SMOKE_TEST=if(profile=="quick")"1" else "0",
    EIVGP_CORES=as.character(cores)),eivgp_adni_driver(data_dir,output_dir))
}

adni_make_C <- function(z) {
  out <- data.frame(
    c_ab42_ab40_3 = ordered(z$c_ab42_ab40_3, levels = 1:3),
    c_gfap_3 = ordered(z$c_gfap_3, levels = 1:3),
    check.names = FALSE
  )
  names(out) <- c("c_ab42_ab40_3", "c_gfap_3")
  out
}


# ============================================================================
# Ocean: settings, frozen inputs and task execution
# ============================================================================

# 01: Configuration, fixed inputs, task selection, execution and provenance.

# 1. Numerical settings are fixed; selectors only choose work, never alter a seed.
ocean_controls <- function() list(
  `UC-GP` = list(n_starts = 8L),
  LVGP = list(n_starts = 8L, max_retries = 2L, max_iter_ini = 100L,
    max_iter_lat = 20L, rescue_iter_ini = 250L, rescue_iter_lat = 80L,
    max_elapsed_seconds = Inf, parallel = TRUE),
  EzGP = list(cv_folds = 2L, maxeval = 100L, cv_score = "nlpd"))

ocean_settings <- function(profile = "paper", action = "plan", folds = 1:5,
  cells = c("O0", "O20", "O50"),
  modules = c("eiv", "original_baselines", "continuation", "complete_case_baselines",
    "prediction", "diagnostics", "surface_comparison"),
  task_ids = character(), cores = 4L, run_id = "local", continuation_iter = NULL, chains = NULL, iter = NULL, burn = NULL, ess_min = NULL) {
  ess_min <- mixedgp_parse_ess_min(ess_min)
  profile <- match.arg(profile, c("paper", "quick", "pilot"))
  action <- match.arg(action, c("plan", "fit", "run", "analyze", "continue", "report"))
  if(action=="run")action<-"fit"
  if(is.null(iter))iter<-if(profile=="paper")50000L else 80L
  if(is.null(burn))burn<-if(profile=="paper")5000L else 30L
  iter<-mixedgp_as_integer_strict(iter,"iter",1L,1L)
  burn<-mixedgp_as_integer_strict(burn,"burn",0L,1L)
  if(burn>=iter)stop("burn must be smaller than iter.")
  if(is.null(chains))chains<-if(profile=="quick")1L else 4L
  chains<-mixedgp_as_integer_strict(chains,"chains",1L,1L)
  if(!is.null(ess_min) && chains<4L)stop("Scientific ESS requires at least four chains.")
  if(profile=="quick" && chains!=1L)stop("quick uses one chain; choose pilot or paper for multiple chains.")
  if(profile=="paper" && action%in%c("plan","fit","continue") &&
     !identical(c(chains,iter,burn),c(4L,50000L,5000L)))
    stop("paper fixes chains=4, iter=50000 and burn=5000; use the profile settings.")
  if(!is.null(continuation_iter))
    continuation_iter <- ocean_continuation_policy(profile,continuation_iter)$n_iter
  folds <- as.integer(unlist(folds)); cells <- as.character(unlist(cells))
  modules <- as.character(unlist(modules))
  task_ids <- as.character(unlist(task_ids))
  task_ids <- task_ids[nzchar(task_ids)]
  stopifnot(length(folds) > 0L, !anyNA(folds), all(folds %in% 1:5),
    length(cells) > 0L, all(cells %in% c("O0", "O20", "O50")),
    all(modules %in% c("eiv", "original_baselines", "continuation", "complete_case_baselines",
      "prediction", "diagnostics", "surface_comparison")),
    length(cores) == 1L,
    is.finite(cores), cores >= 1L, cores == as.integer(cores),
    length(run_id) == 1L, grepl("^[A-Za-z0-9_-]+$", run_id))
  if (profile == "quick") {
    folds <- intersect(folds, 1L); cells <- intersect(cells, c("O0", "O20", "O50"))
    modules <- intersect(modules, c("eiv", "complete_case_baselines", "prediction"))
    if (!length(folds) || !length(cells)) stop("Quick supports fold 1 only.")
  }
  list(profile = profile, action = action, folds = unique(folds), cells = unique(cells),
    modules = unique(modules), task_ids = unique(task_ids), cores = as.integer(cores),
    run_id = run_id, continuation_iter = continuation_iter, ess_min=ess_min,
    iter = iter, burn = burn,
    chains = chains, u_block = 8L)
}

# 2. A complete, stable task catalogue, irrespective of the device's subset.
ocean_task_catalogue <- function(profile = "paper") {
  rows <- list(); add <- function(id, module, fold, cell, method = "EIV-GP", seed = NA_integer_) {
    rows[[length(rows) + 1L]] <<- data.frame(task_id = id, module, fold, cell,
      method, fit_seed = seed, stringsAsFactors = FALSE)
  }
  for (f in if (profile != "quick") 1:5 else 1L) {
    for (cell in if (profile != "quick") c("O50", "O20", "O0") else c("O50", "O20", "O0"))
      add(sprintf("eiv_fold%02d_%s", f, cell), "eiv", f, cell,
        seed = 2026092680L + 100L*f + c(O0=0L, O20=20L, O50=50L)[[cell]])
    if (profile != "quick") add(sprintf("original_fold%02d", f), "original_baselines", f,
      "baseline", "UC-GP+LVGP+EzGP", seed = 2026092660L + 100L*f)
  }
  for(f in if(profile!="quick") 1:5 else 1L) for(cell in c("O20","O50"))
    for(method in c("UC-GP","LVGP","EzGP"))
      add(sprintf("cc_fold%02d_%s_%s",f,cell,gsub("-","",method)),
        "complete_case_baselines",f,cell,method,
        seed=2026101400L+100L*f+as.integer(sub("O","",cell))+match(method,c("UC-GP","LVGP","EzGP")))
  if (profile != "quick") {
    # Every EIV cell is a candidate; current initial diagnostics decide execution.
    for (f in 1:5) for (cell in c("O0","O20","O50"))
      add(sprintf("continue_fold%02d_%s", f, cell), "continuation", f, cell)
  }
  do.call(rbind, rows)
}

ocean_select_tasks <- function(settings) {
  all <- ocean_task_catalogue(settings$profile)
  chosen <- all[all$module %in% settings$modules & all$fold %in% settings$folds &
    (all$cell %in% settings$cells | all$cell == "baseline"), , drop = FALSE]
  if (length(settings$task_ids)) {
    if (length(setdiff(settings$task_ids, chosen$task_id)))
      stop("task_ids not included by the other selectors: ", paste(setdiff(settings$task_ids, chosen$task_id), collapse=", "))
    chosen <- chosen[chosen$task_id %in% settings$task_ids, , drop=FALSE]
  }
  rownames(chosen) <- NULL
  chosen
}

# Adaptive continuation: the rule is fixed before fitting, the cells are not.
# Finite Rhat and bulk ESS determine this limited screen. Undefined
# diagnostics (e.g. constant dictionary occupancy) and tail ESS are also
# recorded for inspection.
ocean_continuation_policy <- function(profile="paper",n_iter=NULL) {
  if(is.null(n_iter))n_iter<-if(profile=="pilot")80L else 20000L
  if(!is.numeric(n_iter) || length(n_iter)!=1L || !is.finite(n_iter) ||
    n_iter<1 || n_iter>.Machine$integer.max || n_iter!=floor(n_iter))
    stop("continuation_iter must be a positive integer")
  list(max_rhat=1.01,min_ess_bulk=400,n_iter=as.integer(n_iter),max_rounds=1L,
    screen="finite Rhat and bulk ESS")
}

# Extension lengths belong to individual continuation tasks, not initial fits.
# Keep historical checkpoint paths stable; record the actual budget separately.
ocean_continuation_iterations <- function(output_dir,profile="paper",fold,cell,
  requested=NULL,register=FALSE) {
  path<-ocean_path(output_dir,"continuation_plus20k",sprintf("fold%02d",fold),cell,
    "continuation_settings.rds")
  saved<-if(file.exists(path))readRDS(path)$n_iter else NULL
  if(is.null(saved) && file.exists(ocean_path(dirname(path),"fit_plus20k.rds")))
    saved<-ocean_continuation_policy(profile)$n_iter # Legacy fixed-budget outputs.
  if(!is.null(saved))saved<-ocean_continuation_policy(profile,saved)$n_iter
  if(register && !is.null(saved) && !is.null(requested) && saved!=requested)
    stop("Continuation length differs from the saved task for fold ",fold," / ",cell,
      ": saved ",saved,", requested ",requested,
      ". Reuse its saved length; a second extension is not supported.")
  n<-if(!is.null(saved))saved else ocean_continuation_policy(profile,requested)$n_iter
  if(register && !file.exists(path))ocean_atomic_rds(list(n_iter=n),path)
  n
}

ocean_continuation_lengths <- function(output_dir,profile=ocean_saved_profile(output_dir)) {
  paths<-list.files(if(ocean_storage_flat(output_dir))output_dir else ocean_path(output_dir,"continuation_plus20k"),
    "^continuation_settings[.]rds$",recursive=TRUE,full.names=TRUE)
  if(!length(paths))return(ocean_continuation_policy(profile)$n_iter)
  sort(unique(vapply(paths,function(path)
    ocean_continuation_policy(profile,readRDS(path)$n_iter)$n_iter,integer(1))))
}

ocean_read_screen <- function(path) {
  p <- ocean_continuation_policy()
  if(!file.exists(path)) return(data.frame(ready=FALSE,max_rhat=NA_real_,min_ess_bulk=NA_real_,
    n_nonfinite_rhat=NA_integer_,n_nonfinite_ess_bulk=NA_integer_,n_nonfinite_ess_tail=NA_integer_,
    screen_pass=NA,diagnostic_md5=NA_character_))
  z <- read.csv(path);stopifnot(all(c("rhat","ess_bulk","ess_tail") %in% names(z)))
  r <- z$rhat[is.finite(z$rhat)];b <- z$ess_bulk[is.finite(z$ess_bulk)]
  mx <- if(length(r)) max(r) else NA_real_;mn <- if(length(b)) min(b) else NA_real_
  data.frame(ready=TRUE,max_rhat=mx,min_ess_bulk=mn,
    n_nonfinite_rhat=sum(!is.finite(z$rhat)),n_nonfinite_ess_bulk=sum(!is.finite(z$ess_bulk)),
    n_nonfinite_ess_tail=sum(!is.finite(z$ess_tail)),
    screen_pass=length(r)>0L && length(b)>0L && mx<=p$max_rhat && mn>=p$min_ess_bulk,
    diagnostic_md5=unname(tools::md5sum(path)))
}

ocean_continuation_decisions <- function(output_dir,profile="paper",write=FALSE,continuation_iter=NULL) {
  rows<-list()
  for(f in if(profile!="quick")1:5 else 1L) for(cell in if(profile!="quick")c("O0","O20","O50") else c("O0","O20","O50")) {
    path<-ocean_path(if(profile!="quick")"full5_class4" else "smoke_class4",sprintf("fold%02d",f),cell,"convergence_parameters.csv")
    path<-ocean_result_relative(output_dir,path)
    screen<-ocean_read_screen(ocean_path(output_dir,path))
    needs<-if(profile=="quick")FALSE else if(screen$ready)!screen$screen_pass else NA
    rows[[length(rows)+1L]]<-cbind(data.frame(fold=f,cell,initial_diagnostics=path),screen,
      needs_continuation=needs,decision=if(is.na(needs))"awaiting_initial" else if(needs)paste0("append_",ocean_continuation_iterations(output_dir,profile,f,cell,continuation_iter)) else "keep_initial")
  }
  d<-do.call(rbind,rows)
  if(write)ocean_atomic_csv(d,ocean_path(output_dir,"continuation_decisions.csv"))
  d
}

ocean_final_eiv_sources <- function(output_dir,profile="paper") {
  decisions<-ocean_continuation_decisions(output_dir,profile)
  rows<-lapply(seq_len(nrow(decisions)),function(i) {
    d<-decisions[i,];continued<-file.exists(ocean_path(output_dir,"continuation_plus20k",
      sprintf("fold%02d",d$fold),d$cell,"fit_plus20k.rds"))
    source<-ocean_path(if(continued)"continuation_plus20k" else if(profile!="quick")"full5_class4" else "smoke_class4",
      sprintf("fold%02d",d$fold),d$cell)
    source<-ocean_result_relative(output_dir,source)
    final<-ocean_read_screen(ocean_path(output_dir,source,"convergence_parameters.csv"))
    data.frame(fold=d$fold,cell=d$cell,source,continued,
      fit_file=if(continued)"fit_plus20k.rds" else "fit.rds",
      available=d$ready && file.exists(ocean_path(output_dir,source,"predictive_metrics.csv")),
      initial_screen_pass=d$screen_pass,final_screen_pass=final$screen_pass)
  })
  do.call(rbind,rows)
}

# Execution and provenance
# Load definitions only. Rmd and 01_setup.R both enter through ocean_execute().
ocean_atomic_rds <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  tmp <- tempfile(".atomic-", tmpdir=dirname(path)); saveRDS(x,tmp)
  if (!file.rename(tmp,path)) stop("Cannot publish ",path)
}
ocean_atomic_csv <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  tmp <- tempfile(".atomic-", tmpdir=dirname(path)); write.csv(x,tmp,row.names=FALSE,na="")
  if (!file.rename(tmp,path)) stop("Cannot publish ",path)
}
ocean_atomic_lines <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  tmp <- tempfile(".atomic-", tmpdir=dirname(path)); writeLines(x,tmp)
  if (!file.rename(tmp,path)) stop("Cannot publish ",path)
}
ocean_hashes <- function(paths) unname(tools::md5sum(paths))
ocean_resolve_path <- function(path,repository_dir) {
  path <- path.expand(as.character(path))
  if(!length(path)) return(character())
  absolute <- grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)",path)
  ifelse(absolute,path,ocean_path(repository_dir,path))
}

# Frozen design and test masks. This never draws a new cohort or calibration set.
ocean_test_masks <- function(test,fold) {
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(202609260L+fold)
  k20<-k50<-rep(FALSE,nrow(test))
  for(c in 1:4) {
    ids<-sample(which(test$C==c));n<-length(ids)
    q20<-round(.2*n);q50<-floor(.5*n)+as.integer(n%%2L==1L && (c+fold)%%2L==0L)
    k20[head(ids,q20)]<-TRUE;k50[head(ids,q50)]<-TRUE
  }
  data.frame(cohort_id=test$cohort_id,C=test$C,observed_U_O20=k20,
    observed_U_O50=k50,mask_seed=202609260L+fold)
}
ocean_read_inputs <- function(data_dir) {
  cohort<-read.csv(ocean_path(data_dir,"cohort_manifest.csv"))
  stopifnot(nrow(cohort)%in%c(200L,500L,1000L),!anyDuplicated(cohort$profile),!anyDuplicated(cohort$cohort_id))
  cohort$C<-findInterval(cohort$U,c(50,100,150))+1L
  cohort$oxygen_class<-factor(cohort$C,levels=1:4,labels=c("<50","50-100","100-150","150-180"))
  cohort
}
ocean_prepare_inputs <- function(data_dir) {
  cohort<-ocean_read_inputs(data_dir);N<-nrow(cohort);ntrain<-as.integer(.8*N);ntest<-N-ntrain
  stopifnot(all(table(factor(cohort$C,levels=1:4))==N/4))
  folds<-list();held<-character();audit<-list()
  for(f in 1:5) {
    tr<-read.csv(ocean_path(data_dir,sprintf("fold%02d_train.csv",f)))
    te<-read.csv(ocean_path(data_dir,sprintf("fold%02d_test.csv",f)))
    tr$C<-findInterval(tr$U,c(50,100,150))+1L;te$C<-findInterval(te$U,c(50,100,150))+1L
    tr$calib_o20<-ocean_diag_flag(tr$calib_o20);tr$calib_o50<-ocean_diag_flag(tr$calib_o50)
    stopifnot(nrow(tr)==ntrain,nrow(te)==ntest,!length(intersect(tr$profile,te$profile)),
      setequal(c(tr$cohort_id,te$cohort_id),cohort$cohort_id),!anyDuplicated(c(tr$cohort_id,te$cohort_id)),
      all(table(factor(tr$C,levels=1:4))==ntrain/4),all(table(factor(te$C,levels=1:4))==ntest/4),
      all(tr$calib_o20<=tr$calib_o50),
      all(table(factor(tr$C[tr$calib_o20],levels=1:4))==.2*ntrain/4),
      all(table(factor(tr$C[tr$calib_o50],levels=1:4))==.5*ntrain/4))
    columns<-c("cohort_id","profile","G2salinity","G2pressure","log_pressure","U","Y")
    for(d in list(tr,te)) {
      ref<-cohort[match(d$cohort_id,cohort$cohort_id),columns]
      stopifnot(isTRUE(all.equal(d[,columns],ref,check.attributes=FALSE)),
        all(is.finite(as.matrix(d[,columns[3:7]]))),all(d$U>0 & d$U<180),
        all(d$G2pressure>=300 & d$G2pressure<=1500),max(abs(d$log_pressure-log(d$G2pressure)))<1e-10)
    }
    masks<-ocean_test_masks(te,f);saved<-ocean_path(data_dir,sprintf("fold%02d_test_U_masks.csv",f))
    stopifnot(file.exists(saved),isTRUE(all.equal(read.csv(saved),masks,check.attributes=FALSE)),
      sum(masks$observed_U_O20)==.2*ntest,sum(masks$observed_U_O50)==.5*ntest)
    for(cell in c("O20","O50")) {
      obs<-tr[[paste0("calib_",tolower(cell))]];known<-masks[[paste0("observed_U_",cell)]]
      center<-mean(tr$U[obs]);scale<-sd(tr$U[obs]);stopifnot(is.finite(scale),scale>0)
      folds[[paste(f,cell,sep="/")]]<-list(train=tr,test=te,observed_train=obs,
        observed_test_mixed=known,center=center,scale=scale,masks=masks)
      audit[[length(audit)+1L]]<-data.frame(fold=f,cell,n_train=ntrain,n_test=ntest,
        n_train_exact=sum(obs),n_test_exact=sum(known),oxygen_center=center,oxygen_sd=scale)
    }
    held<-c(held,te$cohort_id)
  }
  stopifnot(length(held)==N,!anyDuplicated(held))
  list(folds=folds,cohort=cohort,response_sd=sd(cohort$Y),audit=do.call(rbind,audit),
    design=list(N=N,n_train=ntrain,n_test=ntest,n_common_hidden=ntrain/2,
      train_exact=c(O0=0L,O20=.2*ntrain,O50=.5*ntrain),test_exact=c(O0=0L,O20=.2*ntest,O50=.5*ntest)))
}

# Portable identity excludes paths, run_id, selectors, device cores and the
# requested extension length (checked in each continuation task receipt).
ocean_identity <- function(repository_dir,data_dir,settings) {
  packages<-c("eivGP","kergp","LVGP","EzGP")
  versions<-setNames(vapply(packages,function(p)as.character(utils::packageVersion(p)),""),packages)
  paths<-sort(list.files(data_dir,"[.]csv$",full.names=TRUE),method="radix")
  functions<-c("fit_eivgp","continue_eivgp","predict_eivgp","impute_eivgp",
    "run_published_mixedgp_competitors","mixedgp_adapter_ucgp","mixedgp_adapter_lvgp","mixedgp_adapter_ezgp")
  list(schema="ocean-package-CC-v2",profile=settings$profile,
    upstream_commit="677791de6637298ead40d87a534ccd761722d1d5",
    data=setNames(ocean_hashes(paths),basename(paths)),code=eivgp_module_identity(),packages=versions,
    package_functions=setNames(vapply(functions,function(n)digest::digest(
      paste(deparse(get(n,asNamespace("eivGP"))),collapse="\n"),algo="sha256"),""),functions),
    numerics=list(iter=settings$iter,burn=settings$burn,chains=settings$chains,u_block=8L,ess_min=settings$ess_min,
      controls=ocean_cc_controls(settings$profile),n_draw=400L,
      continuation_policy=ocean_continuation_policy(settings$profile)),
    RNG=c("Mersenne-Twister","Inversion","Rejection"))
}


# Explain incompatible saved results before any task or checkpoint is changed.
ocean_assert_identity <- function(saved,current,output_dir) {
  if(identical(saved,current))return(invisible(TRUE))
  fields<-union(names(saved),names(current))
  changed<-fields[!vapply(fields,function(n)identical(saved[[n]],current[[n]]),logical(1))]
  details<-character()
  for(key in c("chains","iter","burn")) {
    old<-saved$numerics[[key]];new<-current$numerics[[key]]
    if(!identical(old,new)) {
      show<-function(x)if(is.null(x))"missing" else paste(x,collapse=",")
      details<-c(details,paste0("  ",key,": saved=",show(old),"; requested=",show(new)))
    }
  }
  # Suggest a fresh root without creating, renaming or deleting any directory.
  fresh<-paste0(output_dir,"_new");i<-1L
  while(file.exists(fresh)) {i<-i+1L;fresh<-paste0(output_dir,"_new",i)}
  stop(paste(c("Existing Ocean results are incompatible with this run.",
    paste0("Output directory: ",output_dir),
    paste0("Changed identity fields: ",paste(changed,collapse=", ")),details,
    "Existing fits and metadata have been preserved; no fitting started.",
    "To start a separate run, repeat your command with:",
    paste0("  --output=",shQuote(fresh)),
    "To resume the existing run, use its original data, package/code and initial sampling settings.",
    "Do not delete workflow_identity.rds to bypass this check."),collapse="\n"),call.=FALSE)
}

ocean_task_files <- function(task,profile="paper",stage="analyze") {
  if(stage %in% c("fit","continue") && task$module %in% c("eiv","continuation")) {
    sub<-ocean_path(if(task$module=="continuation")"continuation_plus20k" else
      if(profile=="quick")"smoke_class4" else "full5_class4",sprintf("fold%02d",task$fold),task$cell)
    return(ocean_path(sub,c(if(task$module=="continuation")"fit_plus20k.rds" else "fit.rds",
      "convergence_parameters.csv","convergence_summary.csv",
      if(task$module=="continuation")"continuation_settings.rds")))
  }
  initial <- if(profile!="quick") "full5_class4" else "smoke_class4"
  if(task$module=="original_baselines") return(ocean_path(initial,sprintf("fold%02d",task$fold),
    c("baseline_metrics.csv","baseline_package_outputs.rds","baseline_status.csv","baseline_predictions.csv")))
  if(task$module=="complete_case_baselines") return(ocean_path("complete_case_baselines/fits",task$task_id,
    c("package_output.rds","metrics.csv","predictions.csv","package_status.csv","design.rds")))
  sub <- ocean_path(if(task$module=="continuation") "continuation_plus20k" else initial,
    sprintf("fold%02d",task$fold),task$cell)
  ocean_path(sub,c(if(task$module=="continuation") "fit_plus20k.rds" else "fit.rds",
    "metrics.csv","predictive_metrics.csv","predictions.csv","mixed_predictions.csv",
    "prediction_distributions.rds","U_posterior_intervals.csv","convergence_parameters.csv",
    "convergence_summary.csv",if(task$module=="continuation")"continuation_settings.rds",
    if(task$cell=="O0") "imputations_model_scale.csv" else "imputations.csv"))
}
ocean_verify_receipt <- function(output_dir,id,identity=NULL) {
  path <- ocean_path(output_dir,"workflow_tasks",paste0(id,".done.rds"))
  if(!file.exists(path)) return(FALSE)
  receipt <- readRDS(path)
  if(!is.null(identity) && !identical(receipt$identity,identity)) stop("Task identity mismatch: ",id)
  files <- ocean_path(output_dir,receipt$files)
  if(!all(file.exists(files)) || !identical(ocean_hashes(files),receipt$hashes))
    stop("Incomplete or changed completed task: ",id)
  TRUE
}

# A single task-local environment prevents initial/continuation helpers colliding.
ocean_run_task <- function(task,repository_dir,output_dir,settings,inputs,identity) {
  env <- new.env(parent=environment())
  env$root <- output_dir; env$smoke <- settings$profile=="quick"
  env$stage <- settings$action
  env$out <- ocean_path(output_dir,if(env$smoke) "smoke_class4" else "full5_class4")
  env$ess_min <- settings$ess_min
  env$iter <- settings$iter; env$burn <- settings$burn; env$chains <- settings$chains
  env$cores <- if(env$smoke) 1L else settings$cores; env$u_block <- 8L
  env$baseline_controls <- ocean_cc_controls(settings$profile)
  env$continuation_iter <- ocean_continuation_iterations(output_dir,settings$profile,
    task$fold,task$cell,settings$continuation_iter)
  env$atomic_csv <- ocean_atomic_csv; env$atomic_lines <- ocean_atomic_lines; env$atomic_rds <- ocean_atomic_rds
  env$set_task <- function(...) invisible(NULL)
  env$write_progress <- function(...) {cat(paste(...,collapse=" "),"\n");flush.console()}
  if(task$module=="complete_case_baselines") {
    ocean_run_complete_case(task,output_dir,settings,inputs)
  } else if(task$module %in% c("eiv","original_baselines")) {
    ocean_initial_functions(env)
    d <- env$read_fold(task$fold); fold_out <- ocean_path(env$out,sprintf("fold%02d",task$fold))
    dir.create(fold_out,recursive=TRUE,showWarnings=FALSE)
    if(task$module=="eiv") env$run_eiv(task$fold,task$cell,d$train,d$test,fold_out)
    else env$run_baselines(task$fold,d$train,d$test,fold_out)
  } else if(task$module=="continuation") {
    dependency <- sprintf("fit_eiv_fold%02d_%s",task$fold,task$cell)
    if(!ocean_verify_receipt(output_dir,dependency,identity)) stop("Missing initial task: ",dependency)
    env$out <- ocean_path(output_dir,"continuation_plus20k")
    env$original_path <- function(fold,cell) ocean_path(output_dir,"full5_class4",sprintf("fold%02d",fold),cell,"fit.rds")
    ocean_continuation_functions(env); env$run_one(task$fold,task$cell)
  } else stop("Unknown application task module: ",task$module)
  invisible(NULL)
}

# Preview the first wave of EIV fits; idle resources require pending work.
# This changes execution concurrency only, never the requested chain count.
# The single CLI budget bounds active workers; outer limits bound resident fit objects.
ocean_stage_workers <- function(cores,stage,n_tasks,chains=4L) {
  per_fit <- if(stage %in% c("fit","continue")) chains else 2L
  as.integer(max(1L,min(n_tasks,ceiling(cores/per_fit))))
}
ocean_stage_plan <- function(settings) {
  data.frame(stage=c("fit","predict","functionals","paper_assets","html"),
    depends_on=c("frozen inputs","saved fits","predictions","scores and diagnostics","paper assets"),
    unit=c("independent chains and native comparator tasks","fixed seeded posterior-draw blocks",
      "fold/calibration and posterior-draw blocks","independent tables and figure files","one document"),
    core_budget=rep(settings$cores,5L),
    concurrency=c("shared pool; resident fits bounded","shared pool; resident fits bounded",
      "shared pool; resident fits bounded","shared pool","serial assembly"),
    stringsAsFactors=FALSE)
}

ocean_core_plan <- function(cores,n_fits,chains=4L) {
  cores<-mixedgp_as_integer_strict(cores,"cores",1L,1L)
  n_fits<-mixedgp_as_integer_strict(n_fits,"n_fits",0L,1L)
  chains<-mixedgp_as_integer_strict(chains,"chains",1L,1L)
  n<-as.integer(min(n_fits,ceiling(cores/chains)))
  allocation<-if(n)as.integer(pmin(chains,pmax(0,cores-(seq_len(n)-1)*chains))) else integer()
  list(total_cores=cores,chains_per_fit=chains,concurrent_fits=n,
    task_cores=allocation,active_cores=sum(allocation),pending_fits=n_fits-n)
}

# Run independent fits against one shared CPU budget. Each fork reserves only
# its assigned chain slots; the coordinator alone publishes aggregate status.
ocean_schedule <- function(indices, caps, cores, FUN, on_start, on_done) {
  stopifnot(length(indices)==length(caps),all(caps>=1L),cores>=1L)
  if(!length(indices))return(invisible(NULL))
  if(.Platform$OS.type!="unix" || cores==1L) {
    for(k in seq_along(indices)) {
      on_start(indices[k],1L,Sys.getpid(),1L)
      on_done(indices[k],FUN(indices[k],1L),0L)
    }
    return(invisible(NULL))
  }
  running<-list();next_task<-1L;active<-0L
  on.exit({
    if(length(running)) {
      for(x in running)try(tools::pskill(x$job$pid),silent=TRUE)
      try(parallel::mccollect(lapply(running,`[[`,"job"),wait=FALSE),silent=TRUE)
    }
  },add=TRUE)
  while(next_task<=length(indices) || length(running)) {
    while(next_task<=length(indices) && active<cores) {
      k<-next_task;i<-indices[k];allocated<-ocean_core_plan(cores-active,1L,caps[k])$task_cores[1L]
      child<-parallel::mcparallel(FUN(i,allocated),mc.set.seed=FALSE,silent=TRUE)
      running[[as.character(child$pid)]]<-list(job=child,index=i,cores=allocated)
      active<-active+allocated
      on_start(i,allocated,child$pid,active)
      next_task<-next_task+1L
    }
    completed<-parallel::mccollect(lapply(running,`[[`,"job"),wait=FALSE)
    if(is.null(completed)) {Sys.sleep(.05);next}
    for(pid in names(completed)) {
      item<-running[[pid]];active<-active-item$cores;running[[pid]]<-NULL
      on_done(item$index,completed[[pid]],active)
    }
  }
  invisible(NULL)
}

# Keep existing execution logs in place; new runs need no jobs/ wrapper.
ocean_job_dir <- function(output_dir,run_id) {
  if(length(run_id)!=1L || !grepl("^[A-Za-z0-9_-]+$",run_id) ||
     run_id %in% c("data","full5_class4","smoke_class4","prediction","diagnostics",
                   "continuation_plus20k","figures","receipts"))
    stop("Choose a distinct execution run_id, such as paper_fit.")
  if(ocean_storage_flat(output_dir))return(ocean_path(output_dir,"logs",run_id))
  legacy <- ocean_path(output_dir,"jobs",run_id)
  if(dir.exists(legacy))legacy else ocean_path(output_dir,run_id)
}

# Unified entry; report never invokes a fit. All cached computations are validated.
ocean_execute <- function(repository_dir,data_dir,output_dir,settings,
  import_dirs=character()) {
  if(settings$action=="fit")ocean_storage_init(output_dir)
  release_budget <- mixedgp_budget_enter(settings$cores)
  on.exit(release_budget(),add=TRUE)
  settings$cores <- mixedgp_budget_cores()
  RNGkind("Mersenne-Twister","Inversion","Rejection")
  import_dirs <- as.character(unlist(import_dirs)); import_dirs <- import_dirs[nzchar(import_dirs)]
  inputs <- ocean_prepare_inputs(data_dir); selected <- ocean_select_tasks(settings)
  stage<-settings$action
  if(stage=="fit")selected<-selected[selected$module!="continuation",,drop=FALSE]
  if(stage=="continue")selected<-selected[selected$module=="continuation",,drop=FALSE]
  if(stage=="analyze")selected<-selected[selected$module%in%c("eiv","continuation"),,drop=FALSE]
  execute<-stage%in%c("fit","analyze","continue")
  if(settings$action=="plan") return(list(tasks=selected,inputs=inputs,settings=settings,stages=ocean_stage_plan(settings)))
  if(execute || length(import_dirs)) {
    suppressPackageStartupMessages(library(eivGP))
    identity <- ocean_identity(repository_dir,data_dir,settings)
    dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
    lock <- ocean_path(output_dir,".workflow.lock")
    if(!dir.create(lock,showWarnings=FALSE)) stop("Output directory is locked: ",lock,
      ". Confirm the previous process has stopped before removing a stale lock.")
    on.exit(unlink(lock,recursive=TRUE),add=TRUE)
    idfile <- ocean_path(output_dir,"workflow_identity.rds")
    if(file.exists(idfile)) {
      ocean_assert_identity(readRDS(idfile),identity,output_dir)
    } else {
      if(length(setdiff(list.files(output_dir,all.files=TRUE,no..=TRUE),c(".workflow.lock",".application.lock","jobs","_meta","logs",basename(ocean_job_dir(output_dir,settings$run_id))))))
        stop("Existing outputs lack distributed provenance. Read them with action=report, or use a fresh directory.")
      ocean_atomic_rds(identity,idfile)
    }
    dir.create(ocean_path(output_dir,"data"),showWarnings=FALSE)
    csv <- list.files(data_dir,"[.]csv$",full.names=TRUE)
    for(f in csv) {
      to <- ocean_path(output_dir,"data",basename(f))
      if(file.exists(to)) stopifnot(identical(ocean_hashes(f),ocean_hashes(to))) else stopifnot(file.copy(f,to))
    }
    for(f in 1:5) {
      path <- ocean_path(output_dir,if(settings$profile!="quick") "full5_class4" else "smoke_class4",sprintf("fold%02d",f),"test_U_masks.csv")
      ocean_atomic_csv(inputs$folds[[paste(f,"O20",sep="/")]]$masks,path)
    }
    ocean_atomic_rds(settings,ocean_path(output_dir,"execution_settings.rds"))
    ocean_atomic_rds(settings,ocean_path(output_dir,paste0("execution_settings_",stage,".rds")))
    ocean_combine_results(import_dirs,output_dir,identity)
    job <- ocean_job_dir(output_dir,settings$run_id)
    status <- if(execute) selected else selected[FALSE,,drop=FALSE]
    status$state <- rep("PENDING",nrow(status)); status$detail <- rep("",nrow(status))
    status$elapsed_seconds <- rep(NA_real_,nrow(status))
    status$allocated_cores<-rep(NA_integer_,nrow(status))
    status$worker_pid<-rep(NA_integer_,nrow(status))
    publish <- function(state,current="none") {
      ocean_atomic_csv(status,ocean_path(job,"task_status.csv"))
      ocean_atomic_lines(c(paste("state:",state),paste("updated:",format(Sys.time(),"%Y-%m-%d %H:%M:%S %z")),
        paste("pid:",Sys.getpid()),paste("current:",current),
        sprintf("executed done: %d/%d; failed: %d; conditional skips: %d",sum(status$state %in% c("DONE","REUSED")),nrow(status)-sum(status$state=="SKIPPED_CONVERGED"),sum(status$state=="FAILED"),sum(status$state=="SKIPPED_CONVERGED")),
        paste("Total core budget:",settings$cores,"; shared across chains and analysis. See execution_schedule.csv and task logs.")),
        ocean_path(job,"PROGRESS.txt"))
    }
    ocean_atomic_csv(inputs$audit,ocean_path(job,"design_audit.csv"))
    capture.output(utils::sessionInfo(),file=ocean_path(job,"sessionInfo.txt"))
    pending<-integer();ids<-character(nrow(selected))
    for(i in if(execute) seq_len(nrow(selected)) else integer()) {
      task <- selected[i,,drop=FALSE]; id <- task$task_id
      if(stage%in%c("fit","continue") && task$module%in%c("eiv","continuation"))id<-paste0("fit_",id)
      ids[i]<-id
      rownames(task) <- NULL # Assignment position is not part of task identity.
      if(task$module=="continuation") {
        dependency<-sprintf("fit_eiv_fold%02d_%s",task$fold,task$cell)
        if(!ocean_verify_receipt(output_dir,dependency,identity)) {
          status$state[i]<-"FAILED";status$detail[i]<-paste("Missing verified initial fit:",dependency)
          publish("RUNNING",id);next
        }
        checkpoint<-ocean_path(output_dir,"continuation_plus20k",sprintf("fold%02d",task$fold),task$cell,"fit_plus20k.rds")
        if(stage=="analyze" && !file.exists(checkpoint)) {
          status$state[i]<-"SKIPPED_NOT_REQUESTED";publish("RUNNING",id);next
        }
        if(stage=="continue") {
          decisions<-ocean_continuation_decisions(output_dir,settings$profile,write=TRUE,
            continuation_iter=settings$continuation_iter)
          d<-decisions[decisions$fold==task$fold & decisions$cell==task$cell,]
          if(nrow(d)==1L && isTRUE(d$ready) && !isTRUE(d$needs_continuation)) {
            status$state[i]<-"SKIPPED_CONVERGED";publish("RUNNING",id);next
          }
          ocean_continuation_iterations(output_dir,settings$profile,task$fold,task$cell,
            requested=settings$continuation_iter,register=TRUE)
        }
      }
      if(ocean_verify_receipt(output_dir,id,identity)) {status$state[i]<-"REUSED";publish("RUNNING",id);next}
      cfg <- ocean_path(output_dir,"workflow_tasks",paste0(id,".config.rds"))
      if(file.exists(cfg) && !identical(readRDS(cfg),list(identity=identity,task=task))) stop("Partial-task identity mismatch: ",id)
      ocean_atomic_rds(list(identity=identity,task=task),cfg)
      pending<-c(pending,i)
    }
    events<-data.frame(task_id=character(),event=character(),time=numeric(),
      worker_pid=integer(),allocated_cores=integer(),active_cores=integer())
    ocean_atomic_csv(events,ocean_path(job,"execution_schedule.csv"))
    event<-function(i,type,active) {
      events<<-rbind(events,data.frame(task_id=ids[i],event=type,time=as.numeric(Sys.time()),
        worker_pid=status$worker_pid[i],allocated_cores=status$allocated_cores[i],active_cores=active))
      ocean_atomic_csv(events,ocean_path(job,"execution_schedule.csv"))
    }
    work<-function(i,allocated) {
      task<-selected[i,,drop=FALSE];rownames(task)<-NULL;id<-ids[i]
      local_settings<-settings;local_settings$cores<-allocated
      log_path<-ocean_path(job,"tasks",paste0(id,".log"))
      dir.create(dirname(log_path),recursive=TRUE,showWarnings=FALSE)
      log<-file(log_path,open="at");old_message<-getConnection(sink.number(type="message"))
      sink(log);sink(log,type="message")
      on.exit({sink(old_message,type="message");sink();close(log)},add=TRUE)
      cat("Task:",id,"PID:",Sys.getpid(),"allocated cores:",allocated,"stage:",stage,"\n")
      t0<-proc.time()[3]
      outcome<-tryCatch({
        ocean_run_task(task,repository_dir,output_dir,local_settings,inputs,identity)
        files<-ocean_result_relative(output_dir,ocean_task_files(task,settings$profile,stage))
        stopifnot(all(file.exists(ocean_path(output_dir,files))))
        ocean_atomic_rds(list(identity=identity,task=task,files=files,
          hashes=ocean_hashes(ocean_path(output_dir,files)),finished=as.character(Sys.time())),
          ocean_path(output_dir,"workflow_tasks",paste0(id,".done.rds")))
        NULL
      },error=function(e)e)
      detail<-if(inherits(outcome,"error"))conditionMessage(outcome) else ""
      if(nzchar(detail))cat("FAILED:",detail,"\n")
      list(state=if(inherits(outcome,"error"))"FAILED" else "DONE",detail=detail,
        elapsed_seconds=as.numeric(proc.time()[3]-t0))
    }
    start_task<-function(i,allocated,pid,active) {
      status$state[i]<<-"RUNNING";status$allocated_cores[i]<<-allocated;status$worker_pid[i]<<-pid
      event(i,"START",active);publish("RUNNING",ids[i])
    }
    finish_task<-function(i,result,active) {
      if(!is.list(result) || is.null(result$state))result<-list(state="FAILED",
        detail="Worker exited without a result; inspect its task log and memory limits.",elapsed_seconds=NA_real_)
      status$state[i]<<-result$state;status$detail[i]<<-result$detail
      status$elapsed_seconds[i]<<-result$elapsed_seconds
      event(i,"FINISH",active);publish("RUNNING",ids[i])
    }
    caps<-ifelse(selected$module[pending]%in%c("eiv","continuation"),settings$chains,1L)
    # Bound simultaneously resident fits, but allow each stage to borrow all
    # free slots. A waiting fit coordinator releases its slot to chain workers.
    outer <- ocean_stage_workers(settings$cores,stage,length(pending),settings$chains)
    mixedgp_budget_yield(function() ocean_schedule(pending,rep(1L,length(pending)),outer,
      function(i,allocated) mixedgp_budget_worker(function() work(i,settings$cores)),
      start_task,finish_task))
    publish(if(any(status$state=="FAILED")) "SELECTED_COMPLETE_WITH_FAILURES" else "SELECTED_COMPLETE")
    ocean_atomic_lines(capture.output(print(status,row.names=FALSE)),ocean_path(job,"SUMMARY.txt"))
    failed <- status$state=="FAILED"
    # Native comparator failures are reported, never allowed to discard valid EIV results.
    if(any(failed & !status$module %in% c("original_baselines","complete_case_baselines")))
      stop("Some EIV tasks failed; inspect ",ocean_path(job,"task_status.csv"))
    if(any(failed))warning("Comparator tasks failed; continuing analysis of available results. See ",
      ocean_path(job,"task_status.csv"),call.=FALSE)
  }
  decisions<-ocean_continuation_decisions(output_dir,settings$profile,write=dir.exists(output_dir),
    continuation_iter=settings$continuation_iter)
  result <- if(!stage%in%c("fit","continue") && "prediction" %in% settings$modules) ocean_prediction_report(output_dir,inputs,
    settings$profile) else NULL
  list(tasks=selected,inputs=inputs,settings=settings,prediction=result,output_dir=output_dir,continuation_decisions=decisions)
}

# Read numerical metadata without refitting or changing saved results.
ocean_saved_profile <- function(output_dir) {
  p<-ocean_path(output_dir,"workflow_identity.rds")
  if(file.exists(p))readRDS(p)$profile else "paper"
}
ocean_continuation_label <- function(output_dir)
  paste(paste0("+",format(ocean_continuation_lengths(output_dir),big.mark=",",scientific=FALSE)),collapse=" / ")


# ============================================================================
# Ocean: core-model and native-comparator calls
# ============================================================================

# 02: Fit -> predict -> impute, using package functions; sourcing never fits.
# Entry points used by 01_setup.R:
#   run_eiv()       : eivGP::fit_eivgp -> predict_eivgp -> impute_eivgp
#   run_baselines(): eivGP:::run_published_mixedgp_competitors, m_vec=4
#   run_one()      : eivGP::continue_eivgp, then predict and impute again
# Native calls inside that pinned eivGP wrapper:
#   UC-GP: kergp::gp / stats::predict (kGauss * q1Symm covariance)
#   LVGP : LVGP::LVGP_fit / LVGP::LVGP_predict
#   EzGP : EzGP::EzGP_fit / EzGP::EzGP_predict (training-only nugget CV)
# RMD-BEGIN/END markers let the Rmd display these exact calls without refitting.
# Numerical arguments, task seeds, masks and optimizer controls are frozen.


# A. Initial EIV and original baseline group
ocean_initial_functions <- function(envir) {
  eval(quote({
# Initial-fit functions are defined in a task-local environment.
# root/out/smoke/iter/burn/chains/cores/u_block are supplied by 01_setup.R.
score <- function(pred, y, fold, method, cell) {
  pred <- as.numeric(pred)
  stopifnot(length(pred) == length(y), all(is.finite(pred)))
  data.frame(fold = fold, cell = cell, method = method, n_test = length(y),
             RMSE = sqrt(mean((pred - y)^2)), MAE = mean(abs(pred - y)))
}

probability_score <- function(draws,point,y,fold,method,cell) {
  qs<-apply(as.matrix(draws),2L,quantile,c(.025,.975))
  nlpd<-eivGP:::normal_mixture_nlpd(y,attr(draws,"conditional_means",exact=TRUE),attr(draws,"conditional_vars",exact=TRUE))
  cbind(score(point,y,fold,method,cell),CRPS=mean(eivGP:::crps_sample_matrix(draws,y)),
    NLPD=nlpd$value,Coverage95=mean(y>=qs[1,]&y<=qs[2,]),Width95=mean(qs[2,]-qs[1,]),
    IntervalScore95=mean(eivGP:::interval_score(qs[1,],qs[2,],y)),NLPD_reason=nlpd$reason)
}

read_fold <- function(fold) {
  train <- read.csv(ocean_path(root, "data", sprintf("fold%02d_train.csv", fold)))
  test <- read.csv(ocean_path(root, "data", sprintf("fold%02d_test.csv", fold)))
  train$C <- as.integer(findInterval(train$U, c(50, 100, 150)) + 1L)
  test$C <- as.integer(findInterval(test$U, c(50, 100, 150)) + 1L)
  design<-ocean_prepare_inputs(ocean_path(root,"data"))
  stopifnot(identical(train$cohort_id,design$folds[[paste(fold,"O20",sep="/")]]$train$cohort_id))
  list(train = train, test = test)
}
run_baselines <- function(fold, train, test, fold_out) {
  metrics_path <- ocean_path(fold_out, "baseline_metrics.csv")
  fit_path <- ocean_path(fold_out,"baseline_package_outputs.rds")
  methods <- c("UC-GP","LVGP","EzGP")
  fit <- if(file.exists(fit_path))readRDS(fit_path) else NULL
  done <- if(is.null(fit))character() else fit$status$method[fit$status$status=="success"]
  pending <- setdiff(methods,done)
  if(!length(pending) && file.exists(metrics_path) &&
     setequal(read.csv(metrics_path)$method,methods))return(invisible("reused output"))
  # Retry only missing methods; preserve their original position-dependent seeds.
  for(method in pending) {
    replacement <- run_published_mixedgp_competitors(
      X_train=as.matrix(train[,c("G2salinity","log_pressure")]),y_train=train$Y,
      C_train=matrix(train$C,ncol=1L,dimnames=list(NULL,"oxygen_class")),
      X_test=as.matrix(test[,c("G2salinity","log_pressure")]),
      C_test=matrix(test$C,ncol=1L,dimnames=list(NULL,"oxygen_class")),
      m_vec=4L,n_draw=400L,seed=2026092660L+fold*100L+match(method,methods)-1L,
      methods=method,strict=FALSE,controls=baseline_controls)
    if(is.null(fit))fit<-replacement else {
      for(n in setdiff(names(replacement),"status"))fit[[n]][[method]]<-replacement[[n]][[method]]
      fit$status<-rbind(fit$status[fit$status$method!=method,,drop=FALSE],replacement$status)
    }
    # Atomic replacement also avoids writing through a read-only source symlink.
    atomic_rds(fit,fit_path)
  }
  fit$status<-fit$status[match(methods,fit$status$method),,drop=FALSE]
  atomic_rds(fit,fit_path)
  atomic_csv(fit$status, ocean_path(fold_out, "baseline_status.csv"))
  metrics <- list()
  predictions <- list()
  for (method in c("UC-GP", "LVGP", "EzGP")) {
    if (!method %in% names(fit$means)) next
    point <- as.numeric(fit$means[[method]])
    metrics[[method]] <- probability_score(fit$draws[[method]],point,test$Y,fold,method,"baseline")
    predictions[[method]] <- data.frame(cohort_id = test$cohort_id,
      method = method, observed_Y = test$Y, predicted_Y = point)
  }
  if (length(metrics)) atomic_csv(do.call(rbind, metrics), metrics_path)
  if (length(predictions)) atomic_csv(do.call(rbind, predictions),
    ocean_path(fold_out, "baseline_predictions.csv"))
  failed <- fit$status$method[fit$status$status != "success"]
  if (length(failed)) stop("Baseline failures: ", paste(failed, collapse = ","))
  invisible("all three baselines fitted")
}
run_eiv <- function(fold, cell, train, test, fold_out) {
  cell_out <- ocean_path(fold_out, cell)
  dir.create(cell_out, showWarnings = FALSE)
  metric_path <- ocean_path(cell_out, "metrics.csv")
  if (file.exists(metric_path) && !identical(get0("stage",inherits=TRUE),"fit")) return(invisible("reused output"))
  # 1. Assemble X/C and reveal U only on the fixed training-calibration mask.
  x_train <- as.matrix(train[, c("G2salinity", "log_pressure")])
  x_test <- as.matrix(test[, c("G2salinity", "log_pressure")])
  c_train <- matrix(train$C, ncol = 1L, dimnames = list(NULL, "oxygen_class"))
  c_test <- matrix(test$C, ncol = 1L, dimnames = list(NULL, "oxygen_class"))
  u_obs <- matrix(NA_real_, nrow(train), 1L, dimnames = list(NULL, "oxygen"))
  mask <- if (cell == "O0") rep(FALSE, nrow(train)) else
    train[[paste0("calib_", tolower(cell))]] == 1L
  u_obs[mask, 1L] <- train$U[mask]
  seed_offset <- c(O0 = 0L, O20 = 20L, O50 = 50L)[[cell]]
  fit_path <- ocean_path(cell_out, "fit.rds")
  if(identical(get0("stage",inherits=TRUE),"analyze") && !file.exists(fit_path))
    stop("Missing fit; run ocean.R --action=fit first: ",fit_path)
  fit <- if (file.exists(fit_path)) readRDS(fit_path) else {
    # 2. Fit: missing U is NA; exact calibration entries are in u_obs.
    # RMD-BEGIN eiv-fit
    model <- eivGP::fit_eivgp(
      X = x_train, y = train$Y, C = c_train, U_obs = u_obs,
      engine = "univariate", m_vec = 4L, standardize_U = TRUE,
      kernel = "se", n_iter = iter, burn = burn, thin = 1L,
      n_chains = chains, parallel = !smoke, n_cores = cores,
      sampler_control = list(u_block_size = u_block),
      seed = 2026092680L + fold * 100L + seed_offset,
      verbose = TRUE
    )
    # RMD-END eiv-fit
    saveRDS(model, fit_path)
    capture.output(summary(model), file = ocean_path(cell_out, "fit_summary.txt"))
    model
  }
  if(identical(get0("stage",inherits=TRUE),"fit") && length(fit$checkpoint$states)>=4L) {
    rows<-mixedgp_study2_diagnostic_rows(x_test,c_test,5L)
    panel<-list(X=x_test[rows,,drop=FALSE],C=as.integer(c_test[rows,1L]),
      U=if(any(mask))test$U[rows]else NULL,n_latent=64L,seed=202610240L+fold*100L+seed_offset)
    fit<-mixedgp_application_ess(fit,get0("ess_min",inherits=TRUE),panel,
      ocean_path(cell_out,"scientific-ess-state.rds"),cores)
    atomic_rds(fit,fit_path)
    mixedgp_save_scientific_trace(fit,ocean_path(cell_out,"scientific-trace-targets"),
      sprintf("Ocean | fold %d | %s",fold,cell))
    unlink(ocean_path(cell_out,"scientific-ess-state.rds"))
  }
  diagnostic_parts <- lapply(
    fit$diagnostics[c("rhat_hyper", "rhat_tau", "rhat_u")],
    function(z) z[, c("parameter", "rhat", "ess_bulk", "ess_tail", "target_pass"),
                  drop = FALSE])
  diagnostics <- do.call(rbind, diagnostic_parts)
  atomic_csv(diagnostics, ocean_path(cell_out, "convergence_parameters.csv"))
  finite_rhat <- diagnostics$rhat[is.finite(diagnostics$rhat)]
  finite_ess <- diagnostics$ess_bulk[is.finite(diagnostics$ess_bulk)]
  atomic_csv(data.frame(iterations=fit$checkpoint$iteration,
                        max_finite_rhat = if (length(finite_rhat)) max(finite_rhat) else NA_real_,
                        min_finite_ess_bulk = if (length(finite_ess)) min(finite_ess) else NA_real_,
                        limited_screen_pass = length(finite_rhat) > 0L && length(finite_ess) > 0L &&
                          max(finite_rhat) <= 1.01 && min(finite_ess) >= 400),
             ocean_path(cell_out, "convergence_summary.csv"))
  if(identical(get0("stage",inherits=TRUE),"fit"))return(invisible("fit saved; analysis deferred"))
  n_draw <- nrow(fit$mcmc$samples_tau)
  stopifnot(is.finite(n_draw), n_draw > 0L)
  set.seed(2026092690L + fold * 100L + seed_offset)
  draw_ids <- sort(sample(seq_len(n_draw), min(if (smoke) 10L else 1000L, n_draw)))
  # 3. Predict withheld Y from (X,C), integrating over unobserved U.
  # RMD-BEGIN eiv-predict
  response_draws <- eivGP::predict_eivgp(fit, new_X = x_test, new_C = c_test,
    target = "response", draw_ids = draw_ids, n_per_draw = 1L, joint = FALSE,
    seed = 2026092700L + fold * 100L + seed_offset)
  # RMD-END eiv-predict
  point <- colMeans(as.matrix(response_draws))
  atomic_csv(data.frame(cohort_id = test$cohort_id, observed_Y = test$Y,
    predicted_Y = point), ocean_path(cell_out, "predictions.csv"))
  metrics <- score(point, test$Y, fold, "EIV-GP", cell)

  # Freeze nested test masks once per fold using the prespecified seed, without Y.
  masks<-ocean_test_masks(test,fold);known20<-masks$observed_U_O20;known50<-masks$observed_U_O50
  known<-switch(cell,O0=rep(FALSE,nrow(test)),O20=known20,O50=known50)
  atomic_csv(data.frame(cohort_id=test$cohort_id,C=test$C,observed_U_O20=known20,observed_U_O50=known50,mask_seed=202609260L+fold),
    ocean_path(fold_out,"test_U_masks.csv"))
  mixed<-response_draws
  if (any(known)) {
    xu<-eivGP::predict_eivgp(fit,new_X=x_test[known,,drop=FALSE],new_U=matrix(test$U[known],ncol=1L),new_U_scale="raw",
      target="response",draw_ids=draw_ids,n_per_draw=1L,joint=FALSE,seed=2026093000L+fold*100L+match(cell,c("O0","O20","O50")))
    mixed[,known]<-xu
    for (a in c("conditional_means","conditional_vars")) {
      comp<-attr(response_draws,a,exact=TRUE);comp[,known]<-attr(xu,a,exact=TRUE);attr(mixed,a)<-comp
    }
  }
  stopifnot(identical(mixed[,!known,drop=FALSE],response_draws[,!known,drop=FALSE]))
  saveRDS(list(xc=response_draws,mixed=mixed,known=known,draw_ids=draw_ids),
    ocean_path(cell_out,"prediction_distributions.rds"))
  atomic_csv(data.frame(cohort_id=test$cohort_id,observed_Y=test$Y,predicted_Y=colMeans(mixed),test_U_observed=known),
    ocean_path(cell_out,"mixed_predictions.csv"))
  atomic_csv(rbind(probability_score(response_draws,point,test$Y,fold,"EIV-GP xc",cell),
    probability_score(mixed,colMeans(mixed),test$Y,fold,"EIV-GP mixed",cell)),ocean_path(cell_out,"predictive_metrics.csv"))

  # 4. Impute hidden training U; O0 uses model coordinates, O20/O50 raw units.
  hidden <- which(!is.finite(u_obs[, 1L]))
  if (cell == "O0") {
    ## Without calibration, raw U is not affinely identified by the package.
    imp <- eivGP::impute_eivgp(fit, rows = hidden, draw_ids = draw_ids, scale = "model")
    imp_point <- colMeans(matrix(imp[, , 1L], nrow = length(draw_ids)))
    atomic_csv(data.frame(cohort_id = train$cohort_id[hidden], C = train$C[hidden],
      true_U_raw = train$U[hidden], EIV_U_model_mean = imp_point),
      ocean_path(cell_out, "imputations_model_scale.csv"))
  } else {
    # RMD-BEGIN eiv-impute
    imp <- eivGP::impute_eivgp(fit, rows = hidden, draw_ids = draw_ids, scale = "raw")
    # RMD-END eiv-impute
    imp_point <- colMeans(matrix(imp[, , 1L], nrow = length(draw_ids)))
    calib_means <- tapply(train$U[mask], train$C[mask], mean)
    imp_base <- as.numeric(calib_means[as.character(train$C[hidden])])
    atomic_csv(data.frame(cohort_id = train$cohort_id[hidden], C = train$C[hidden],
      true_U = train$U[hidden], EIV_U_mean = imp_point,
      class_calibration_mean = imp_base), ocean_path(cell_out, "imputations.csv"))
    metrics <- rbind(metrics,
      score(imp_point, train$U[hidden], fold, "EIV U imputation", cell),
      score(imp_base, train$U[hidden], fold, "calibration class mean", cell))
  }
  iq<-apply(matrix(imp[,,1L],nrow=length(draw_ids)),2L,quantile,c(.025,.975))
  atomic_csv(data.frame(cohort_id=train$cohort_id[hidden],C=train$C[hidden],posterior_mean=imp_point,q025=iq[1,],q975=iq[2,],scale=if(cell=="O0")"model" else "raw"),ocean_path(cell_out,"U_posterior_intervals.csv"))
  atomic_csv(metrics, metric_path)
  tau <- as.matrix(fit$mcmc$samples_tau)
  atomic_csv(data.frame(boundary = seq_len(ncol(tau)),
    q025_model_scale = apply(tau, 2L, quantile, .025),
    median_model_scale = apply(tau, 2L, median),
    q975_model_scale = apply(tau, 2L, quantile, .975)),
    ocean_path(cell_out, "threshold_posterior_model_scale.csv"))
  invisible(paste0("Y RMSE ", round(metrics$RMSE[1L], 3)))
}


  }), envir=envir)
  invisible(envir)
}

# B. Diagnostic-selected +20k continuation
ocean_continuation_functions <- function(envir) {
  eval(quote({
# Diagnostic-selected +20k continuation; independent of unrelated folds.
# Selection uses current initial diagnostics in 01_setup.R; prediction seeds and saved RNG states are retained.
diagnostic_row <- function(fit) {
  tables <- fit$diagnostics[c("rhat_hyper", "rhat_tau", "rhat_u")]
  rhats <- unlist(lapply(tables, function(z) z$rhat))
  bulk <- unlist(lapply(tables, function(z) z$ess_bulk))
  max_rhat <- max(rhats[is.finite(rhats)])
  min_ess <- min(bulk[is.finite(bulk)])
  data.frame(max_rhat = max_rhat, min_ess_bulk = min_ess,
    pass_screen = max_rhat <= 1.01 && min_ess >= 400)
}
score <- function(pred, y, fold, cell, method) {
  pred <- as.numeric(pred)
  stopifnot(length(pred) == length(y), all(is.finite(pred)))
  data.frame(fold = fold, cell = cell, method = method, n_test = length(y),
    RMSE = sqrt(mean((pred - y)^2)), MAE = mean(abs(pred - y)))
}


probability_score<-function(draws,point,y,fold,cell,method) {
  q<-apply(draws,2L,quantile,c(.025,.975))
  nlpd<-eivGP:::normal_mixture_nlpd(y,attr(draws,"conditional_means",exact=TRUE),attr(draws,"conditional_vars",exact=TRUE))
  cbind(score(point,y,fold,cell,method),CRPS=mean(eivGP:::crps_sample_matrix(draws,y)),NLPD=nlpd$value,
    Coverage95=mean(y>=q[1,]&y<=q[2,]),Width95=mean(q[2,]-q[1,]),
    IntervalScore95=mean(eivGP:::interval_score(q[1,],q[2,],y)),NLPD_reason=nlpd$reason)
}

run_one <- function(fold, cell) {
  task_out <- ocean_path(out, sprintf("fold%02d", fold), cell)
  dir.create(task_out, recursive = TRUE, showWarnings = FALSE)
  metrics_path <- ocean_path(task_out, "metrics.csv")
  fit_path <- ocean_path(task_out, "fit_plus20k.rds")
  if (file.exists(metrics_path)) {
    fit <- readRDS(fit_path)
    set_task(fold, cell, "DONE", "reused complete continuation",
             diag = diagnostic_row(fit))
    return(invisible(NULL))
  }
  set_task(fold, cell, "RUNNING", "continuing saved chains")
  write_progress("RUNNING", sprintf("fold%02d %s", fold, cell), "continuing original chains")
  t0 <- Sys.time()
  if(identical(get0("stage",inherits=TRUE),"analyze") && !file.exists(fit_path))
    stop("Missing continuation; run ocean.R --action=continue first.")
  fit <- if (file.exists(fit_path)) readRDS(fit_path) else {
    original <- readRDS(original_path(fold, cell))
    # Continue the saved chains/RNG states, then repeat prediction and imputation.
    # RMD-BEGIN eiv-continue
    continued <- eivGP::continue_eivgp(original, n_iter = continuation_iter,
      parallel = cores > 1L, n_cores = cores, verbose = TRUE)
    # RMD-END eiv-continue
    atomic_rds(continued, fit_path)
    rm(original); gc(FALSE)
    continued
  }
  diag <- diagnostic_row(fit)
  diagnostic_parts<-lapply(fit$diagnostics[c("rhat_hyper","rhat_tau","rhat_u")],function(z)z[,c("parameter","rhat","ess_bulk","ess_tail","target_pass"),drop=FALSE])
  atomic_csv(do.call(rbind,diagnostic_parts),ocean_path(task_out,"convergence_parameters.csv"))
  atomic_csv(data.frame(max_finite_rhat=diag$max_rhat,min_finite_ess_bulk=diag$min_ess_bulk,limited_screen_pass=diag$pass_screen),ocean_path(task_out,"convergence_summary.csv"))
  atomic_csv(data.frame(fold = fold, cell = cell, stage = paste0("plus",continuation_iter),
    max_rhat = diag$max_rhat, min_ess_bulk = diag$min_ess_bulk,
    pass_screen = diag$pass_screen), ocean_path(task_out, "diagnostics.csv"))
  capture.output(summary(fit), file = ocean_path(task_out, "fit_summary.txt"))
  if(identical(get0("stage",inherits=TRUE),"continue"))return(invisible("continued fit saved; analysis deferred"))

  train <- read.csv(ocean_path(root, "data", sprintf("fold%02d_train.csv", fold)))
  test <- read.csv(ocean_path(root, "data", sprintf("fold%02d_test.csv", fold)))
  train$C <- as.integer(findInterval(train$U, c(50, 100, 150)) + 1L)
  test$C <- as.integer(findInterval(test$U, c(50, 100, 150)) + 1L)
  stopifnot((nrow(train)+nrow(test))%in%c(200L,500L,1000L),
    !length(intersect(train$profile,test$profile)),all(table(test$C)==nrow(test)/4))
  x_test <- as.matrix(test[, c("G2salinity", "log_pressure")])
  c_test <- matrix(as.integer(test$C), ncol = 1L,
                   dimnames = list(NULL, "oxygen_class"))
  n_draw <- nrow(fit$mcmc$samples_tau)
  set.seed(2026092800L + fold * 100L + match(cell, c("O0", "O20", "O50")))
  draw_ids <- sort(sample(seq_len(n_draw), min(2000L, n_draw)))
  response_draws <- eivGP::predict_eivgp(fit, new_X = x_test, new_C = c_test,
    target = "response", draw_ids = draw_ids, n_per_draw = 1L,
    joint = FALSE, seed = 2026092900L + fold * 100L + match(cell, c("O0", "O20", "O50")))
  point <- colMeans(as.matrix(response_draws))
  metrics <- score(point, test$Y, fold, cell, "EIV-GP")
  atomic_csv(data.frame(cohort_id = test$cohort_id, observed_Y = test$Y,
    predicted_Y = point), ocean_path(task_out, "predictions.csv"))

  mask <- if (cell == "O0") rep(FALSE, nrow(train)) else
    train[[paste0("calib_", tolower(cell))]] == 1L

  # Use the fixed test observation masks for both initial and continued fits.
  saved_masks<-read.csv(ocean_path(root,"full5_class4",sprintf("fold%02d",fold),"test_U_masks.csv"))
  stopifnot(identical(as.character(saved_masks$cohort_id),as.character(test$cohort_id)))
  known<-if(cell=="O0")rep(FALSE,nrow(test)) else as.logical(saved_masks[[paste0("observed_U_",cell)]])
  stopifnot(sum(known)==switch(cell,O0=0,O20=.2,O50=.5)*nrow(test))
  mixed<-response_draws
  if(any(known)) {
    xu<-eivGP::predict_eivgp(fit,new_X=x_test[known,,drop=FALSE],new_U=matrix(test$U[known],ncol=1L),new_U_scale="raw",
      target="response",draw_ids=draw_ids,n_per_draw=1L,joint=FALSE,seed=2026093000L+fold*100L+match(cell,c("O0","O20","O50")))
    mixed[,known]<-xu
    for(a in c("conditional_means","conditional_vars")) {
      comp<-attr(response_draws,a,exact=TRUE);comp[,known]<-attr(xu,a,exact=TRUE);attr(mixed,a)<-comp
    }
  }
  stopifnot(identical(mixed[,!known,drop=FALSE],response_draws[,!known,drop=FALSE]))
  atomic_rds(list(xc=response_draws,mixed=mixed,known=known,draw_ids=draw_ids),ocean_path(task_out,"prediction_distributions.rds"))
  atomic_csv(data.frame(cohort_id=test$cohort_id,observed_Y=test$Y,predicted_Y=colMeans(mixed),test_U_observed=known),ocean_path(task_out,"mixed_predictions.csv"))
  atomic_csv(rbind(probability_score(response_draws,point,test$Y,fold,cell,"EIV-GP xc"),
    probability_score(mixed,colMeans(mixed),test$Y,fold,cell,"EIV-GP mixed")),ocean_path(task_out,"predictive_metrics.csv"))
  hidden <- which(!mask)
  imp <- eivGP::impute_eivgp(fit, rows = hidden, draw_ids = draw_ids,
                      scale = if (cell == "O0") "model" else "raw")
  imp_point <- colMeans(matrix(imp[, , 1L], nrow = length(draw_ids)))
  if (cell == "O0") {
    atomic_csv(data.frame(cohort_id = train$cohort_id[hidden], C = train$C[hidden],
      true_U_raw = train$U[hidden], EIV_U_model_mean = imp_point),
      ocean_path(task_out, "imputations_model_scale.csv"))
  } else {
    calib_means <- tapply(train$U[mask], train$C[mask], mean)
    imp_base <- as.numeric(calib_means[as.character(train$C[hidden])])
    atomic_csv(data.frame(cohort_id = train$cohort_id[hidden], C = train$C[hidden],
      true_U = train$U[hidden], EIV_U_mean = imp_point,
      class_calibration_mean = imp_base),
      ocean_path(task_out, "imputations.csv"))
    metrics <- rbind(metrics,
      score(imp_point, train$U[hidden], fold, cell, "EIV U imputation"),
      score(imp_base, train$U[hidden], fold, cell, "calibration class mean"))
  }
  iq<-apply(matrix(imp[,,1L],nrow=length(draw_ids)),2L,quantile,c(.025,.975))
  atomic_csv(data.frame(cohort_id=train$cohort_id[hidden],C=train$C[hidden],posterior_mean=imp_point,q025=iq[1,],q975=iq[2,],
    scale=if(cell=="O0")"model" else "raw"),ocean_path(task_out,"U_posterior_intervals.csv"))
  atomic_csv(metrics, metrics_path)
  minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  set_task(fold, cell, "DONE", paste0("Y RMSE ", round(metrics$RMSE[1], 3)),
           minutes = minutes, diag = diag)
  write_progress("RUNNING", sprintf("fold%02d %s", fold, cell), "continuation and scoring complete")
}


  }), envir=envir)
  invisible(envir)
}

# E. Native complete-case comparators. No hidden U enters native fit/predict.
# Pilot budgets are deliberately separate from paper budgets and recorded.
ocean_cc_controls <- function(profile="paper") {
  if(profile=="paper") return(ocean_controls())
  list(`UC-GP`=list(n_starts=1L),
    LVGP=list(n_starts=1L,max_retries=1L,max_iter_ini=20L,max_iter_lat=5L,
      rescue_iter_ini=30L,rescue_iter_lat=10L,max_elapsed_seconds=120,parallel=FALSE),
    EzGP=list(cv_folds=2L,maxeval=15L,cv_score="nlpd"))
}
ocean_run_complete_case <- function(task,output_dir,settings,inputs) {
  d <- inputs$folds[[paste(task$fold,task$cell,sep="/")]]
  train <- d$train[d$observed_train,,drop=FALSE]
  test <- d$test[d$observed_test_mixed,,drop=FALSE]
  stopifnot(nrow(train)==sum(d$observed_train),
    nrow(test)==sum(d$observed_test_mixed),!anyNA(train$U),!anyNA(test$U),
    !length(intersect(train$cohort_id,test$cohort_id)))
  # Standardize within the complete training subset; neither test U nor hidden
  # training U defines transformations, optimization or model selection.
  cols <- c("G2salinity","log_pressure","U")
  center <- colMeans(train[,cols]); spread <- vapply(train[,cols],sd,0)
  stopifnot(all(is.finite(spread)),all(spread>0))
  transform <- function(z) sweep(sweep(as.matrix(z[,cols]),2L,center,"-"),2L,spread,"/")
  dir <- ocean_path(output_dir,"complete_case_baselines/fits",task$task_id)
  dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  design <- list(protocol="native-CC-observed-u-subset-v1",profile=settings$profile,
    fold=task$fold,cell=task$cell,method=task$method,seed=task$fit_seed,
    training_ids=train$cohort_id,test_ids=test$cohort_id,center=center,scale=spread,
    controls=ocean_cc_controls(settings$profile),n_draw=400L,
    omitted_test_ids=d$test$cohort_id[!d$observed_test_mixed])
  ocean_atomic_rds(design,ocean_path(dir,"design.rds"))
  ans <- eivGP:::run_published_mixedgp_competitors(
    X_train=transform(train),y_train=train$Y,C_train=matrix(train$C,ncol=1L),
    X_test=transform(test),C_test=matrix(test$C,ncol=1L),m_vec=4L,n_draw=400L,
    seed=task$fit_seed,methods=task$method,strict=FALSE,controls=design$controls)
  ocean_atomic_rds(ans,ocean_path(dir,"package_output.rds"))
  ocean_atomic_csv(ans$status,ocean_path(dir,"package_status.csv"))
  if(!task$method %in% names(ans$means) || any(ans$status$status!="success"))
    stop("CC native fitting failed: ",paste(ans$status$message,collapse="; "))
  # The original probability scoring function accepts any matched test size.
  env <- new.env(parent=environment());ocean_initial_functions(env)
  point <- as.numeric(ans$means[[task$method]])
  score <- env$probability_score(ans$draws[[task$method]],point,test$Y,
    task$fold,paste0(task$method,"(CC)"),task$cell)
  score$profile <- settings$profile;score$n_train <- nrow(train)
  ocean_atomic_csv(score,ocean_path(dir,"metrics.csv"))
  ocean_atomic_csv(data.frame(cohort_id=test$cohort_id,observed_Y=test$Y,
    predicted_Y=point),ocean_path(dir,"predictions.csv"))
  invisible(design)
}


# ============================================================================
# Ocean: public workflow API
# ============================================================================

# Ocean package API; scientific core fitting functions are defined in the unchanged core modules.
ocean_workflow_api <- function() {
  ns<-asNamespace("eivGP");n<-grep("^ocean_",ls(ns,all.names=TRUE),value=TRUE)
  setNames(lapply(n,function(x)get(x,ns,inherits=FALSE)),n)
}
run_ocean_analysis <- function(data_dir,output_dir,action="plan",profile="paper",cores=4L,
  repository_dir=NULL,settings=NULL,import_dirs=character(),cc_output_dir="",diagnostics_dir="",
  representative_fold=1L,...) {
  if(is.null(repository_dir))repository_dir<-getwd()
  if(is.null(settings))settings<-ocean_settings(profile=profile,action=action,cores=cores,...)
  ocean_analysis_impl(repository_dir,data_dir,output_dir,settings,import_dirs,
    cc_output_dir,diagnostics_dir,representative_fold)
}
ocean_analysis_impl <- function(repository_dir,data_dir,output_dir,settings,import_dirs=character(),
  cc_output_dir="",diagnostics_dir="",representative_fold=1L) {
  release_budget <- mixedgp_budget_enter(settings$cores)
  on.exit(release_budget(),add=TRUE)
  settings$cores <- mixedgp_budget_cores()
  a<-ocean_execute(repository_dir,data_dir,output_dir,settings,import_dirs)
  if(settings$action%in%c("plan","fit","continue"))return(a)
  if("prediction" %in% settings$modules && !is.null(a$prediction)) {
    a$cc_review<-ocean_cc_prediction_tables(output_dir,a$inputs,a$prediction,
      cc_output_dir=if(nzchar(cc_output_dir))cc_output_dir else output_dir,
      profile=settings$profile,expected_folds=settings$folds)
    a$cc_review$SD_Y<-a$inputs$response_sd
  }
  if(settings$profile!="quick" && any(c("diagnostics","surface_comparison") %in% settings$modules)) {
    available<-a$prediction
    if(is.null(available))available<-ocean_prediction_report(output_dir,a$inputs,settings$profile)
    ready<-all(available$availability$available) && all(vapply(1:5,function(f)
      file.exists(ocean_path(output_dir,"full5_class4",sprintf("fold%02d",f),"baseline_package_outputs.rds")),logical(1)))
    if(ready) {
      d<-ocean_run_diagnostics(output_dir,data_dir,
        report_dir=if(nzchar(diagnostics_dir))diagnostics_dir else ocean_path(output_dir,"diagnostics"),
        representative_fold=representative_fold,
        build_prediction=FALSE,prediction_report=available,
        build_surfaces="surface_comparison" %in% settings$modules)
      a$diagnostic_run<-d;a$diagnostics<-d$diagnostics;a$surfaces<-d$surfaces
    }
  }
  a$paper_assets <- ocean_export_application(output_dir,a,settings,
    data_dir=data_dir,diagnostics_dir=if(nzchar(diagnostics_dir))diagnostics_dir else ocean_path(output_dir,"diagnostics"))
  a
}

# ---- numerical paper workflow: isolated helpers avoid application name collisions ----
eivgp_numerical_api <- function(repo_root) {
  NUMERICAL_REPO_ROOT <- normalizePath(repo_root,mustWork=TRUE)
  numerical_default_data <- function(study)ocean_path(path.expand("~/eivGP-results"),"data",study)
  numerical_default_output <- function(study)file.path(path.expand("~/eivGP-results"),study)
  numerical_package_check <- function()invisible(TRUE) # checkout CLI validates locked environment
F <- function(name) getFromNamespace(name, "eivGP")
atomic_rds <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp<-tempfile(tmpdir=dirname(path));on.exit(unlink(tmp));saveRDS(F("mixedgp_pack_artifact")(x),tmp)
  if(!file.rename(tmp,path)) stop("Cannot publish checkpoint: ",path)
  mixedgp_retire_saved_checkpoints(x,path)
}
retire_conditional_fit <- function(path, states_path) {
  saved<-readRDS(states_path)
  # Preserve every conditional chain for later target diagnostics and continuation.
  # The selected prediction states are an additional integration cache, not a replacement.
  saved
}
hash_object <- function(x) {
  p<-tempfile();on.exit(unlink(p));saveRDS(x,p,version=2);unname(tools::md5sum(p))
}
parse_cli <- function(args) {
  out<-list(); for(a in args) {
    if(!grepl("^--[^=]+=.+$",a)) stop("Use --name=value: ",a)
    k<-sub("^--([^=]+)=.*$","\\1",a);out[[k]]<-sub("^--[^=]+=","",a)
  };out
}
integers <- function(x) {
  z<-unlist(lapply(strsplit(x,",",fixed=TRUE)[[1]],function(a) {
    b<-as.integer(strsplit(a,":",fixed=TRUE)[[1]]);if(length(b)==2)seq.int(b[1],b[2]) else b
  }));if(anyNA(z)||any(z<1)||anyDuplicated(z))stop("Invalid integer selection: ",x);z
}
settings <- function(a) {
  allowed<-c("action","profile","data-root","output","reps","calib","methods","outer","inner-draws",
    "latent","gp-response-latent","cc-response-latent","iter","warmup","inner-iter","inner-warmup","cores","seed","full-fit-map")
  if(length(setdiff(names(a),allowed)))stop("Unknown option: ",paste(setdiff(names(a),allowed),collapse=","))
  get<-function(k,default)if(is.null(a[[k]]))default else a[[k]]
  quick<-identical(get("profile","paper"),"quick")
  pilot<-identical(get("profile","paper"),"pilot")
  if(!get("profile","paper")%in%c("paper","quick","pilot"))stop("profile must be paper, quick or pilot")
  cfg<-list(full_fit_map=get("full-fit-map",""),profile=get("profile","paper"),data_root=get("data-root",Sys.getenv("EIVGP_DATA_ROOT",numerical_default_data("study2"))),
    output=get("output",Sys.getenv("EIVGP_MECHANISM_ROOT",numerical_default_output("study2/mechanism"))),reps=integers(get("reps",if(quick)"1"else if(pilot)"1,25"else"1:50")),
    calib=integers(get("calib","20,50")),methods=strsplit(get("methods",if(pilot)"FULL,CUT,UC-GP-CC,LVGP-CC,EzGP-CC"else"FULL,CUT,UC-GP,LVGP,EzGP,UC-GP-CC,LVGP-CC,EzGP-CC"),",",fixed=TRUE)[[1]])
  defaults<-if(quick)c(outer=4,inner_draws=1,latent=4,iter=40,warmup=20,inner_iter=40,inner_warmup=20,cores=1,seed=77101) else if(pilot)
    c(outer=4,inner_draws=4,latent=16,iter=6000,warmup=2000,inner_iter=1500,inner_warmup=500,cores=4,seed=77101) else
    c(outer=64,inner_draws=4,latent=64,iter=20000,warmup=5000,inner_iter=4000,inner_warmup=1000,cores=4,seed=77101)
  defaults<-c(defaults,gp_response_latent=1L,cc_response_latent=4L)
  for(k in names(defaults))cfg[[k]]<-as.integer(get(gsub("_","-",k),defaults[[k]]))
  if(is.null(a[["cc-response-latent"]]))cfg$cc_response_latent<-min(4L,cfg$latent)
  if(any(!is.finite(unlist(cfg[names(defaults)])))||any(unlist(cfg[names(defaults)])<1))stop("Numeric controls must be positive")
  if(cfg$iter<=cfg$warmup||cfg$inner_iter<=cfg$inner_warmup||cfg$outer%%4!=0)stop("iterations must exceed warmup; outer must be divisible by four")
  if(cfg$gp_response_latent>cfg$latent||cfg$cc_response_latent>cfg$latent)stop("Response integration draws cannot exceed latent draws")
  if(cfg$latent<2||cfg$inner_draws>cfg$inner_iter-cfg$inner_warmup)stop("Need at least two latent draws and enough retained conditional draws")
  if(cfg$outer/4>cfg$iter-cfg$warmup)stop("Too few measurement draws for outer selection")
  if(!all(cfg$calib%in%c(20,50)))stop("Study II mechanism component uses calibration 20 and 50 only")
  if(!all(cfg$methods%in%c("CUT","FULL","UC-GP","LVGP","EzGP","UC-GP-CC","LVGP-CC","EzGP-CC")))stop("Unknown method")
  if(!nzchar(cfg$output)||!nzchar(cfg$data_root))stop("Supply --data-root and --output")
  if(grepl("Dropbox|Overleaf",cfg$output,ignore.case=TRUE))stop("Keep output outside the synced manuscript tree")
  cfg
}
# This API intentionally does not accept y, X, or evaluation truth.
measurement_cut <- function(C,U_obs,cfg,seed) {
  disk <- mixedgp_sampler_disk(list(component="cut-measurement",C=C,U_obs=U_obs,
    iter=cfg$iter,warmup=cfg$warmup,outer=cfg$outer,seed=seed))
  if(!is.null(disk))on.exit(unlink(file.path(disk,"LOCK"),recursive=TRUE),add=TRUE)
  checkpoint_enabled <- is.finite(getOption("eivGP.checkpoint_seconds",Inf)) ||
    is.finite(getOption("eivGP.checkpoint_iterations",Inf))
  run<-function(chain) {
    set.seed(seed+10000L*chain)
    ctx<-list(C=C,U_obs=U_obs,n=nrow(C),p=2L,q=ncol(C),d=ncol(U_obs),m_vec=rep(4L,ncol(C)),
      observed=which(complete.cases(U_obs)),missing=which(!complete.cases(U_obs)),measurement="probit",ident="lower_triangular")
    ctx$priors<-F("mixedgp_v030_priors")(list(),ctx$p,ctx$d,ctx$m_vec)
    st<-F("mixedgp_v030_initial_state")(ctx)
    logmeas<-function(s)F("mixedgp_v030_measurement_loglik")(s,F("mixedgp_v030_decode")(s,ctx),ctx)
    K<-cfg$iter-cfg$warmup;states<-vector("list",cfg$outer/4)
    selected<-floor(((seq_len(cfg$outer/4)-.5)/(cfg$outer/4))*K)+1L
    trace<-matrix(NA_real_,K,length(ctx$missing)*ctx$d+ctx$q*ctx$d+sum(ctx$m_vec))
    path <- if(!is.null(disk))file.path(disk,sprintf("chain-%02d.rds",chain))else NULL
    start <- 0L; elapsed_before <- 0; started <- proc.time()[["elapsed"]]
    if(checkpoint_enabled&&!is.null(path)&&file.exists(path)) {
      cp<-readRDS(path); start<-cp$iteration;st<-cp$state;states<-cp$states
      trace<-mixedgp_checkpoint_restore_draws(list(trace=cp$trace),list(trace=trace))$trace
      elapsed_before<-cp$elapsed_seconds;F("mixedgp_restore_rng_state")(cp$rng)
    }
    last_save<-last_progress<-proc.time()[["elapsed"]];last_iteration<-start
    checkpoint<-function(it,force=FALSE) {
      if(is.null(path))return(invisible(NULL))
      now<-proc.time()[["elapsed"]];elapsed<-elapsed_before+now-started
      saving<-mixedgp_checkpoint_due(now,last_save,it-last_iteration,force)
      if(saving) {
        mixedgp_checkpoint_write(list(iteration=it,state=st,states=states,trace=trace[seq_len(max(0L,it-cfg$warmup)),,drop=FALSE],
          rng=F("mixedgp_rng_state")(),elapsed_seconds=elapsed),path)
        last_save<<-proc.time()[["elapsed"]];last_iteration<<-it
      }
      if(force || saving || now-last_progress>=60) {
        elapsed<-elapsed_before+proc.time()[["elapsed"]]-started
        remaining<-if(it>0)elapsed/it*(cfg$iter-it)else NA_real_
        mixedgp_timing_record(sub("[.]rds$","-timing.csv",path),data.frame(
          updated_at=format(Sys.time(),tz="UTC",usetz=TRUE),chain=chain,iteration=it,
          target_iteration=cfg$iter,checkpoint_iteration=last_iteration,
          phase=if(it==cfg$iter)"complete" else if(it<=cfg$warmup)"warmup" else "sampling",
          elapsed_seconds=elapsed,remaining_seconds=remaining,
          estimated_finish_utc=if(is.finite(remaining))format(Sys.time()+remaining,tz="UTC",usetz=TRUE)else NA_character_))
        last_progress<<-proc.time()[["elapsed"]]
      }
    }
    if(start<cfg$iter)for(it in seq.int(start+1L,cfg$iter)) {
      for(rows in split(ctx$missing,ceiling(seq_along(ctx$missing)/8))) {
        old<-as.vector(st$U[rows,,drop=FALSE]);base<-st
        proposal<-F("mixedgp_v030_ess")(old,function(v){s<-base;s$U[rows,]<-v;logmeas(s)})
        st$U[rows,]<-proposal$value
      }
      for(j in seq_len(ctx$q)) {
        active<-seq_len(min(j,ctx$d));np<-length(st$e_pi[[j]]);base<-st
        update<-function(v){s<-base;s$e_pi[[j]]<-v[seq_len(np)];s$e_A[j,active]<-v[np+seq_along(active)];s}
        v<-F("mixedgp_v030_ess")(c(st$e_pi[[j]],st$e_A[j,active]),function(v)logmeas(update(v)))$value
        st<-update(v)
      }
      if(it>cfg$warmup) {
        k<-it-cfg$warmup;dec<-F("mixedgp_v030_decode")(st,ctx)
        trace[k,]<-c(st$U[ctx$missing,,drop=FALSE],dec$A,unlist(dec$pi))
        if(k%in%selected)states[[match(k,selected)]]<-list(U=st$U,A=dec$A,tau=lapply(dec$tau,function(t)t[-c(1,length(t))]),chain=chain,iteration=it)
      }
      checkpoint(it)
    }
    checkpoint(cfg$iter,force=TRUE)
    list(states=states,trace=trace)
  }
  z<-mixedgp_parallel_lapply(1:4,run,n_cores=min(4L,cfg$cores),.preserve_rng=TRUE)
  if(any(vapply(z,inherits,logical(1),"try-error")))stop("Measurement chain failed")
  mis<-which(!complete.cases(U_obs))
  labels<-c(unlist(lapply(1:2,function(k)sprintf("U[%d,%d]",mis,k))),
    unlist(lapply(1:2,function(k)sprintf("A[%d,%d]",seq_len(ncol(C)),k))),
    unlist(lapply(seq_len(ncol(C)),function(j)sprintf("pi[%d,%d]",j,1:4))))
  arr<-array(NA_real_,c(nrow(z[[1]]$trace),4,ncol(z[[1]]$trace)),
    dimnames=list(NULL,paste0("chain",1:4),labels))
  for(j in 1:4)arr[,j,]<-z[[j]]$trace
  # Constant constrained loading entries have no diagnostic interpretation.
  keep<-apply(arr,3,function(x)stats::sd(as.vector(x))>0)
  diagnostics<-posterior::summarise_draws(posterior::as_draws_array(arr[,,keep,drop=FALSE]))
  # Interleave chains, so first/second-half integration checks include every chain.
  states<-unlist(lapply(seq_len(cfg$outer/4),function(i)lapply(z,function(c)c$states[[i]])),recursive=FALSE)
  list(states=states,diagnostics=diagnostics,trace=arr,
    target="p(U,psi | C,U_obs); same category-probability and loading priors as FULL")
}
diagnostic_pass <- function(z) {
  nrow(z)>0 && all(is.finite(z$rhat) & z$rhat<=1.01 &
    is.finite(z$ess_bulk) & z$ess_bulk>=400 & is.finite(z$ess_tail) & z$ess_tail>=400)
}
record_gp_diagnostics <- function(fit,path,conditional=FALSE) {
  z<-fit$diagnostics$ess_key
  if(conditional)z<-z[grepl("^(V$|r$|theta_x|J\\[|dictionary_occupancy)",z$parameter),,drop=FALSE]
  else z<-z[is.finite(z$posterior_sd)&z$posterior_sd>0,,drop=FALSE]
  write.csv(z,path,row.names=FALSE)
  finite<-is.finite(z$rhat)&is.finite(z$ess_bulk)&is.finite(z$ess_tail)
  write.csv(data.frame(n_parameter=nrow(z),n_nonfinite=sum(!finite),
    n_finite_failure=sum(finite & (z$rhat>1.01|z$ess_bulk<400|z$ess_tail<400)),
    worst_rhat=if(any(is.finite(z$rhat)))max(z$rhat,na.rm=TRUE)else NA_real_,
    min_bulk=if(any(is.finite(z$ess_bulk)))min(z$ess_bulk,na.rm=TRUE)else NA_real_,
    min_tail=if(any(is.finite(z$ess_tail)))min(z$ess_tail,na.rm=TRUE)else NA_real_),
    sub("[.]csv$","_summary.csv",path),row.names=FALSE)
  diagnostic_pass(z)
}
fit_full <- function(train,U_obs,cfg,seed,conditional=FALSE) {
  eivGP::fit_eivgp(train$X,train$y,train$C,U_obs=U_obs,latent_dim=2L,engine="multivariate",
    m_vec=rep(4L,4),ident="lower_triangular",standardize_U=FALSE,n_chains=4L,
    n_iter=if(conditional)cfg$inner_iter else cfg$iter,burn=if(conditional)cfg$inner_warmup else cfg$warmup,
    seed=seed,parallel=cfg$cores>1,n_cores=cfg$cores,verbose=FALSE)
}
read_matched_full <- function(path,train,U_obs) {
  fit<-F("mixedgp_read_artifact")(path)
  if(is.null(fit$mcmc)&&!is.null(fit$fit))fit<-fit$fit
  equal<-function(a,b)isTRUE(all.equal(unname(a),unname(b),check.attributes=FALSE,tolerance=1e-12))
  expected<-F("mixedgp_v030_priors")(list(),ncol(train$X),2L,rep(4L,4))
  if(!identical(fit$sampler_version,"0.3.1")||!equal(fit$priors,expected)||
     !equal(fit$data$X_raw,train$X)||!equal(fit$data$y_raw,train$y)||
     !equal(fit$data$C_ord,train$C)||!equal(fit$data$U_obs,U_obs)||
     !identical(fit$data$ident,"lower_triangular")||
     !all(fit$data$U_center==0)||!all(fit$data$U_scale==1)||
     !identical(fit$kernel$name,"se"))stop("Saved FULL fit does not match data, calibration, priors or model: ",path)
  fit
}
select_gp <- function(fit,number,measurement=NULL) {
  mm<-fit$mcmc;info<-mm$mcmc_draw_info
  ids<-unlist(lapply(split(seq_len(nrow(info)),info$chain),function(i)i[unique(round(seq(1,length(i),length.out=number)))]))
  lapply(ids,function(i) {
    st<-if(is.null(measurement))list(U=mm$samples_U[i,,],A=mm$samples_A[i,,],
      tau=F("unflatten_tau")(mm$samples_tau[i,],rep(4L,4)))else measurement
    c(st,list(logtheta=mm$samples_logtheta[i,],sigma2=mm$samples_sigma2[i],gp_chain=info$chain[i]))
  })
}
# Same integration points for all CC methods. Only C_test enters this function.
prospective_u <- function(states,C_test,L,seed) {
  set.seed(seed)
  lapply(states,function(s) {
    a<-F("sample_u_given_c_ordprobit_dispatch")(C_test[rep(seq_len(nrow(C_test)),each=L),,drop=FALSE],
      s$A,s$tau,rep(4L,4),latent_sampler="minimax_tilting")
    array(as.numeric(a),c(L,nrow(C_test),2L))
  })
}
# Components of the posterior distribution of m and f, plus predictive y.
# The covariance of the integrated GP is retained: mean(cov), not mean(diag(cov)).
gp_evaluate <- function(states,train,test,cfg,seed) {
  xy<-F("mixedgp_v030_standardize_xy")(train$X,train$y)
  Xt<-sweep(sweep(test$X,2,xy$X_center,"-"),2,xy$X_scale,"/")
  us<-prospective_u(states,test$C,cfg$latent,seed)
  N<-nrow(test$X);S<-length(states)
  response_draws<-if(is.null(cfg$gp_response_latent))1L else cfg$gp_response_latent
  if(response_draws<1L||response_draws>cfg$latent)stop("Invalid GP response integration budget")
  out<-lapply(c("y_mean","y_var","m_mean","m_var","f_mean","f_var"),function(x)matrix(NA_real_,N,S))
  names(out)<-c("y_mean","y_var","m_mean","m_var","f_mean","f_var")
  out$y_mean<-out$y_var<-matrix(NA_real_,N,S*response_draws)
  parts <- mixedgp_stage_lapply(seq_len(S),function(s) {
    ans <- list()
    st<-states[[s]]
    predict<-function(x,u,cov=FALSE,noisy=FALSE)F("gp_predict_draw_general")(xy$X,st$U,xy$y,x,u,
      st$logtheta,st$sigma2,noisy=noisy,return_cov=cov)
    f<-predict(Xt,test$U);ans$f_mean<-xy$y_center+xy$y_scale*f$mean;ans$f_var<-xy$y_scale^2*f$var
    # Equal weight per posterior state and prospective input draw. This budget
    # is distinct from integration of m and does not alter fitted states.
    ans$y_mean <- ans$y_var <- matrix(NA_real_,N,response_draws)
    ans$m_mean <- ans$m_var <- numeric(N)
    for(rr in seq_len(response_draws)) {
      yp<-predict(Xt,matrix(us[[s]][rr,,],N,2),noisy=TRUE)
      col<-(s-1L)*response_draws+rr
      ans$y_mean[,rr]<-xy$y_center+xy$y_scale*yp$mean;ans$y_var[,rr]<-xy$y_scale^2*yp$var
    }
    for(i in seq_len(N)) {
      mp<-predict(Xt[rep(i,cfg$latent),,drop=FALSE],matrix(us[[s]][,i,],cfg$latent,2),cov=TRUE)
      ans$m_mean[i]<-xy$y_center+xy$y_scale*mean(mp$mean)
      ans$m_var[i]<-xy$y_scale^2*max(0,mean(mp$cov))
    }
    ans
  })
  for(name in names(out)) out[[name]] <- do.call(cbind,lapply(parts,`[[`,name))
  out
}
categorical_evaluate <- function(method,train,test,cfg,seed,path) {
  if(file.exists(path))return(readRDS(path))
  controls<-if(cfg$profile=="quick")switch(method,`UC-GP`=list(n_starts=2L),
    LVGP=list(n_starts=1L,max_retries=1L,max_iter_ini=10L,max_iter_lat=5L),
    EzGP=list(tau_fractions=c(.01,.04),cv_folds=2L,maxeval=10L))else if(!is.null(cfg$controls))cfg$controls[[method]]else list()
  model<-stage_time(dirname(path),"optimizer",do.call(eivGP::fit_mixedgp_competitor,c(list(method=method,X=train$X,y=train$y,C=train$C,
    new_X=test$X,new_C=test$C,m_vec=if(is.null(cfg$m_vec))rep(4L,4)else cfg$m_vec,n_draw=1L,seed=seed),controls)))
  atomic_rds(model,sub("predictions.rds$","fit.rds",path))
  out<-list(y_mean=matrix(model$predictive_mean,ncol=1),y_var=matrix(model$predictive_variance,ncol=1),
    m_point=as.numeric(model$predictive_mean),optimization_status=model$optimization_status)
  atomic_rds(out,path);out
}
cc_predict <- function(method,model,X,C,m_vec=rep(4L,4)) {
  if(method=="UC-GP") {
    nd<-F("make_mixedgp_competitor_data")(X,C,m_vec)
    p<-stats::predict(model$fit,newdata=nd,type="UK",seCompute=TRUE,covCompute=FALSE,lightReturn=TRUE)
    mu<-as.numeric(p$mean);v<-as.numeric(p$sd)^2
  } else {
    p<-if(method=="LVGP")LVGP::LVGP_predict(cbind(X,C),model$fit,MSE_on=1) else EzGP::EzGP_predict(cbind(X,C),model$fit,MSE_on=1)
    mu<-as.numeric(p$Y_hat);v<-F("predictive_variance_diagonal")(p$MSE,nrow(X),context=method)
  }
  list(mean=mu,var=v+model$fitted_noise_variance)
}
cc_evaluate <- function(method,train,test,cal,us,cfg,seed,path) {
  base<-sub("-CC$","",method)
  d<-ncol(train$U);m_vec<-if(is.null(cfg$m_vec))rep(4L,4)else cfg$m_vec
  if(file.exists(path))model<-readRDS(path)else {
    controls<-if(cfg$profile=="quick")switch(base,`UC-GP`=list(n_starts=2L),LVGP=list(n_starts=1L,max_retries=1L,max_iter_ini=10L,max_iter_lat=5L),EzGP=list(tau_fractions=c(.01,.04),cv_folds=2L,maxeval=10L))else if(!is.null(cfg$controls))cfg$controls[[base]]else list()
    model<-stage_time(dirname(path),"optimizer",do.call(eivGP::fit_mixedgp_competitor,c(list(method=base,X=cbind(train$X[cal,,drop=FALSE],train$U[cal,,drop=FALSE]),
      y=train$y[cal],C=train$C[cal,,drop=FALSE],new_X=cbind(test$X[rep(1L,2L),,drop=FALSE],matrix(rep(us[[1]][1,1,],each=2L),2,d)),
      new_C=test$C[rep(1L,2L),,drop=FALSE],m_vec=if(is.null(cfg$m_vec))rep(4L,4)else cfg$m_vec,n_draw=1L,seed=seed),controls)))
    atomic_rds(model,path)
  }
  N<-nrow(test$X);S<-length(us)*cfg$latent;mu<-v<-matrix(NA_real_,N,S)
  parts<-mixedgp_stage_lapply(seq_along(us),function(s) {
    a<-b<-matrix(NA_real_,N,cfg$latent)
    for(i in seq_len(N)) {
      p<-cc_predict(base,model,cbind(test$X[rep(i,cfg$latent),,drop=FALSE],matrix(us[[s]][,i,],cfg$latent,d)),test$C[rep(i,cfg$latent),,drop=FALSE],m_vec=m_vec)
      a[i,]<-p$mean;b[i,]<-p$var
    }
    list(mean=a,var=b)
  })
  mu<-do.call(cbind,lapply(parts,`[[`,"mean"));v<-do.call(cbind,lapply(parts,`[[`,"var"))
  response_draws<-if(is.null(cfg$cc_response_latent))min(cfg$latent,4L)else cfg$cc_response_latent
  if(response_draws<1L||response_draws>cfg$latent)stop("Invalid CC response integration budget")
  ids<-unlist(lapply(seq_along(us),function(s)(s-1)*cfg$latent+unique(round(seq(1,cfg$latent,length.out=response_draws)))))
  list(y_mean=mu[,ids,drop=FALSE],y_var=v[,ids,drop=FALSE],m_point=rowMeans(mu),optimization_status=model$optimization_status)
}
normal_abs <- function(a,s) {s<-pmax(s,1e-12);2*s*dnorm(a/s)+a*(2*pnorm(a/s)-1)}
mixture_score <- function(y,mu,v) {
  if(any(!is.finite(c(y,mu,v)))||any(v< -1e-10))stop("Invalid predictive components")
  v<-pmax(v,1e-14);sd<-sqrt(v);lp<-dnorm(y,mu,sd,log=TRUE);mx<-max(lp)
  q<-function(p)uniroot(function(x)mean(pnorm(x,mu,sd))-p,c(min(mu-12*sd),max(mu+12*sd)),tol=1e-8)$root
  lo<-q(.025);hi<-q(.975)
  # Block pairwise calculation bounds memory for CC mixtures (64*64 components).
  pair<-0;for(ids in split(seq_along(mu),ceiling(seq_along(mu)/128)))
    pair<-pair+sum(normal_abs(outer(mu[ids],mu,"-"),sqrt(outer(v[ids],v,"+"))))
  c(sq_error=(mean(mu)-y)^2,crps=mean(normal_abs(y-mu,sd))-.5*pair/length(mu)^2,
    nlpd=-(mx+log(mean(exp(lp-mx)))),coverage=as.numeric(y>=lo&&y<=hi),width=hi-lo)
}
score_result <- function(pred,dat,cal,truth,method) {
  train<-dat$train;test<-dat$test;key<-function(C)apply(C,1,paste,collapse=":")
  strata<-ifelse(key(test$C)%in%key(train$C),"seen","unseen")
  cc_support<-key(test$C)%in%key(train$C[cal,,drop=FALSE])
  rows<-list()
  add<-function(target,values,group=strata,component=rep(NA_integer_,length(group))) {
    z<-as.data.frame(values);z$target<-target;z$id<-seq_len(nrow(z));z$stratum<-group;z$coordinate<-component
    z$seen_in_complete_cases<-if(nrow(z)==length(cc_support))cc_support else NA
    rows[[length(rows)+1L]]<<-z
  }
  add("y",do.call(rbind,mixedgp_stage_lapply(seq_along(test$y),function(i)mixture_score(test$y[i],pred$y_mean[i,],pred$y_var[i,]))))
  mp<-if(!is.null(pred$m_point))pred$m_point else rowMeans(pred$m_mean)
  add("m",data.frame(sq_error=(mp-truth$m)^2,crps=NA_real_,nlpd=NA_real_,coverage=NA_real_,width=NA_real_))
  if(!is.null(pred$f_mean)) {
    add("f",do.call(rbind,mixedgp_stage_lapply(seq_along(truth$f),function(i)mixture_score(truth$f[i],pred$f_mean[i,],pred$f_var[i,]))))
    # m intervals use GP covariance under integration, with finite-L approximation.
    rows[[2]]$coverage<-vapply(seq_along(truth$m),function(i)mixture_score(truth$m[i],pred$m_mean[i,],pred$m_var[i,])["coverage"],numeric(1))
    rows[[2]]$width<-vapply(seq_along(truth$m),function(i)mixture_score(truth$m[i],pred$m_mean[i,],pred$m_var[i,])["width"],numeric(1))
  }
  if(!is.null(pred$U)) {
    miss<-setdiff(seq_len(nrow(train$X)),cal)
    for(d in 1:2) {
      u<-pred$U[,miss,d,drop=FALSE];dim(u)<-c(dim(pred$U)[1],length(miss))
      lo<-apply(u,2,quantile,.025);hi<-apply(u,2,quantile,.975);t<-train$U[miss,d]
      add("U",data.frame(sq_error=(colMeans(u)-t)^2,crps=NA_real_,nlpd=NA_real_,coverage=as.numeric(t>=lo&t<=hi),width=hi-lo),
        rep("missing_training",length(miss)),rep(d,length(miss)))
      rows[[length(rows)]]$id<-miss
    }
  }
  out<-do.call(rbind,rows);out$method<-method;out
}
aggregate_scores <- function(jobs,root) {
  rows<-list();status<-list()
  for(j in seq_len(nrow(jobs))) {
    p<-file.path(root,jobs$job[j],"scores.rds");err<-file.path(root,jobs$job[j],"failure.rds")
    dg<-file.path(root,jobs$job[j],"diagnostics.rds")
    passed<-if(file.exists(dg))readRDS(dg)$pass else NA
    state<-if(file.exists(p))"complete"else if(file.exists(err))"failed"else"pending"
    status[[j]]<-data.frame(jobs[j,],status=state,diagnostics_pass=passed,message=if(state=="failed")readRDS(err)$message else "")
    if(state=="complete") {z<-readRDS(p);z$rep<-jobs$rep[j];z$calib<-jobs$calib[j];z$diagnostics_pass<-passed;rows[[length(rows)+1]]<-z}
  }
  write.csv(do.call(rbind,status),file.path(root,"status.csv"),row.names=FALSE)
  if(file.exists(file.path(root,"job_registry.rds")))
    F("report_numerical_convergence")(setNames(root,basename(root)),paste0(root,"-convergence-report"))
  if(!length(rows))return(invisible(NULL))
  points<-do.call(rbind,rows);points$coordinate[is.na(points$coordinate)]<-0L;write.csv(points,file.path(root,"point_scores.csv"),row.names=FALSE)
  both<-rbind(points,transform(points,stratum="overall"))
  groups<-split(both,interaction(both$rep,both$calib,both$method,both$target,both$stratum,both$coordinate,drop=TRUE))
  reps<-do.call(rbind,lapply(groups,function(z){
    out<-z[1,c("rep","calib","method","target","stratum","coordinate")];out$n<-nrow(z);out$diagnostics_pass<-all(z$diagnostics_pass)
    for(k in c("sq_error","crps","nlpd","coverage","width"))out[[k]]<-if(all(is.na(z[[k]])))NA_real_ else mean(z[[k]],na.rm=TRUE)
    out$rmse<-sqrt(out$sq_error);out
  }))
  write.csv(reps,file.path(root,"replicate_scores.csv"),row.names=FALSE)
  groups<-split(reps,interaction(reps$calib,reps$method,reps$target,reps$stratum,reps$coordinate,drop=TRUE))
  sums<-do.call(rbind,lapply(groups,function(z){
    out<-z[1,c("calib","method","target","stratum","coordinate")];out$n_rep<-nrow(z);out$n_diagnostics_pass<-sum(z$diagnostics_pass,na.rm=TRUE);out$n_points<-sum(z$n)
    for(k in c("rmse","crps","nlpd","coverage","width")) {
      v<-z[[k]];out[[k]]<-if(all(is.na(v)))NA_real_ else mean(v,na.rm=TRUE)
      out[[paste0(k,"_se")]]<-if(sum(is.finite(v))>1)sd(v,na.rm=TRUE)/sqrt(sum(is.finite(v)))else NA_real_
    };out
  }))
  write.csv(sums,file.path(root,"summary.csv"),row.names=FALSE)
  # Pair on replication, calibration, target and stratum; never pool test rows as independent replications.
  pairs<-list();methods<-unique(reps$method)
  if(length(methods)>1)for(a in seq_len(length(methods)-1))for(b in (a+1):length(methods)) {
    z<-merge(reps[reps$method==methods[a],],reps[reps$method==methods[b],],by=c("rep","calib","target","stratum","coordinate"))
    if(nrow(z)) {z$method_a<-methods[a];z$method_b<-methods[b];for(k in c("rmse","crps","nlpd","coverage","width"))z[[paste0(k,"_difference")]]<-z[[paste0(k,".x")]]-z[[paste0(k,".y")]];pairs[[length(pairs)+1]]<-z}
  }
  if(length(pairs)) {
    paired<-do.call(rbind,pairs)
    write.csv(paired,file.path(root,"paired_differences.csv"),row.names=FALSE)
    summarize_pairs<-function(z,keys) {
      gs<-split(z,do.call(interaction,c(z[keys],list(drop=TRUE))))
      do.call(rbind,lapply(gs,function(d){
        out<-d[1,keys,drop=FALSE];out$n_rep<-nrow(d)
        for(k in paste0(c("rmse","crps","nlpd","coverage","width"),"_difference")) {
          v<-d[[k]];out[[k]]<-if(any(is.finite(v)))mean(v,na.rm=TRUE)else NA_real_
          out[[paste0(k,"_mcse")]]<-if(sum(is.finite(v))>1)sd(v,na.rm=TRUE)/sqrt(sum(is.finite(v)))else NA_real_
        };out
      }))
    }
    keys<-c("calib","method_a","method_b","target","stratum","coordinate")
    write.csv(summarize_pairs(paired,keys),file.path(root,"paired_summary.csv"),row.names=FALSE)
    si<-merge(paired[paired$stratum=="seen",],paired[paired$stratum=="unseen",],
      by=c("rep","calib","method_a","method_b","target","coordinate"),suffixes=c("_seen","_unseen"))
    if(nrow(si)) {
      for(k in paste0(c("rmse","crps","nlpd","coverage","width"),"_difference"))si[[k]]<-si[[paste0(k,"_unseen")]]-si[[paste0(k,"_seen")]]
      si$stratum<-"unseen_minus_seen"
      write.csv(si,file.path(root,"paired_interactions.csv"),row.names=FALSE)
      write.csv(summarize_pairs(si,keys),file.path(root,"interaction_summary.csv"),row.names=FALSE)
    }
  }
  F("mixedgp_cut_trace_report")(root)
  mechanism_report(root)
  invisible(sums)
}
# Stage-level wall time includes parallel work; it is not summed CPU time.
stage_time <- function(folder,stage,expr) {
  dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  start<-proc.time()[["elapsed"]];ok<-FALSE
  on.exit({
    row<-data.frame(stage=stage,seconds=proc.time()[["elapsed"]]-start,success=ok,
      finished_utc=format(Sys.time(),tz="UTC",usetz=TRUE))
    path<-file.path(folder,paste0("stage-",stage,".csv"))
    if(file.exists(path))row<-rbind(read.csv(path),row)
    write.csv(row,path,row.names=FALSE)
  })
  value<-tryCatch(force(expr),error=function(e){if(is.null(e$stage))e$stage<-stage;stop(e)});ok<-TRUE;value
}
mechanism_report <- function(root) {
  escape<-function(x)gsub(">","&gt;",gsub("<","&lt;",gsub("&","&amp;",as.character(x),fixed=TRUE),fixed=TRUE),fixed=TRUE)
  table<-function(d)paste0("<table><tr>",paste0("<th>",escape(names(d)),"</th>",collapse=""),"</tr>",
    paste(apply(d,1,function(r)paste0("<tr>",paste0("<td>",escape(r),"</td>",collapse=""),"</tr>")),collapse=""),"</table>")
  status<-read.csv(file.path(root,"status.csv"));protocol<-readRDS(file.path(root,"protocol.rds"))
  paths<-list.files(root,pattern="^stage-.*[.]csv$",recursive=TRUE,full.names=TRUE)
  timing<-if(length(paths))do.call(rbind,lapply(paths,function(p){z<-read.csv(p);z$job<-basename(dirname(p));z}))else data.frame()
  if(nrow(timing))write.csv(timing,file.path(root,"stage_timing.csv"),row.names=FALSE)
  body<-c("<!doctype html><meta charset='utf-8'><title>Numerical comparison report</title>",
    "<style>body{font:16px system-ui;max-width:1300px;margin:35px auto;padding:20px}table{border-collapse:collapse;font-size:13px;display:block;overflow:auto}th,td{border:1px solid #ddd;padding:6px}th{background:#eef}p{max-width:950px}</style>",
    "<h1>Paired numerical comparisons</h1>",
    paste0("<p>Profile: <strong>",escape(protocol$settings$profile),"</strong>. No fits are excluded for diagnostic warnings. This saved-results report does not fit models.</p>"),
    "<p><strong>Caution:</strong> A pilot measures execution, mixing and cost; a small replication count and limited outer states cannot establish comparative performance or nested Monte Carlo accuracy. Replication MCSEs are not MCMC errors. Latent coordinates, when available, are summarized separately (coordinate 0 denotes non-coordinate targets). Complete-case m intervals are unavailable because the adapters do not supply the required joint covariance.</p>",
    "<h2>Settings</h2>",paste0("<pre>",escape(paste(capture.output(str(protocol$settings)),collapse="\n")),"</pre>"),
    "<h2>Completion and diagnostic flags</h2>",table(status))
  for(f in c("summary.csv","paired_summary.csv","interaction_summary.csv"))if(file.exists(file.path(root,f))) {
    d<-read.csv(file.path(root,f));if("stratum"%in%names(d)&&f=="summary.csv")d<-d[d$stratum=="overall",]
    body<-c(body,paste0("<h2>",escape(f),"</h2>"),table(d))
  }
  if(nrow(timing))body<-c(body,"<h2>Stage wall times (seconds)</h2>",table(timing))
  body<-c(body,"<p>Inspect measurement_diagnostics.csv and full/conditional diagnostic CSVs alongside flags. integration_stability.csv compares two subsets of prediction components; it is an integration sensitivity screen, not an R-hat or a proof of stationarity. Full chains, failure records and checkpoints are retained.</p>")
  figures<-list.files(file.path(root,"cut-traces"),"[.]png$",full.names=TRUE)
  if(length(figures))body<-c(body,"<h2>CUT scientific traces (descriptive; no ESS requirement)</h2>",
    "<p>Conditional GP traces cover every retained iteration with the outer measurement state fixed. Outer traces summarize the retained integration states, with finite inner Monte Carlo error; they are not the full measurement-chain history. No fit is continued or excluded on these diagnostics.</p>",
    vapply(figures,function(p)paste0("<img style='width:100%' src='cut-traces/",basename(p),"'>"),character(1)))
  writeLines(body,file.path(root,"report.html"))
  invisible(file.path(root,"report.html"))
}
run_mechanism <- function(cfg,action="plan",script_path=NULL) {
  release_budget <- mixedgp_budget_enter(cfg$cores)
  on.exit(release_budget(),add=TRUE)
  cfg$cores <- mixedgp_budget_cores()
  if(!action%in%c("plan","data","run","report"))stop("action must be plan, data, run or report")
  numerical_package_check()
  jobs<-expand.grid(rep=cfg$reps,calib=cfg$calib,method=cfg$methods,stringsAsFactors=FALSE)
  jobs$job<-sprintf("rep%03d-cal%02d-%s",jobs$rep,jobs$calib,jobs$method)
  dir.create(cfg$output,recursive=TRUE,showWarnings=FALSE)
  files<-setNames(file.path(cfg$data_root,"primary_q4_calibration",sprintf("study2_primary_q4_n100_ntest100_m4_rep%03d.rds",cfg$reps)),cfg$reps)
  if(action=="plan") {
    print(jobs);cat(nrow(jobs),"jobs; missing datasets will be generated on run. No fits started.\n")
    return(invisible(jobs))
  }
  if(action %in% c("data","run")) {
    for(r in cfg$reps)if(!file.exists(files[as.character(r)])) {
      artifact<-F("make_study2_synthetic_dataset")(rep_id=r,scenario="primary",n=100L,n_test=100L,
        q=4L,m=4L,calib_grid=c(0L,20L,50L,80L),data_seed_base=1000000L)
      F("store_mixedgp_synthetic_dataset")(artifact,dirname(files[as.character(r)]),overwrite=FALSE)
    }
  }
  if(any(!file.exists(files)))stop("Missing datasets; run with --action=data first")
  for(r in cfg$reps) {
    artifact<-readRDS(files[as.character(r)])
    expected<-list(n=100L,n_test=100L,q=4L,m=4L,d=2L,calib_grid=c(0L,20L,50L,80L))
    if(!identical(artifact$study,"study2")||!identical(artifact$scenario,"primary")||
       !identical(artifact$rep_id,as.integer(r))||
       !isTRUE(all.equal(artifact$design[names(expected)],expected,check.attributes=FALSE)))stop("Frozen design mismatch: ",files[as.character(r)])
  }
  if(action=="data")return(invisible(files))
  pkgs<-c("eivGP","TruncatedNormal","posterior","kergp","LVGP","EzGP")
  versions<-setNames(vapply(pkgs,function(p)if(requireNamespace(p,quietly=TRUE))as.character(packageVersion(p))else"missing",character(1)),pkgs)
  # A root is one scientific configuration; job selection/cores can change on resume.
  scientific<-cfg[setdiff(names(cfg),c("reps","methods","calib","output","data_root","cores"))]
  identity<-list(settings=scientific,versions=versions,R=as.character(getRversion()),
    code=c(mechanism=if(is.null(script_path))hash_object(body(run_mechanism)) else unname(tools::md5sum(script_path)),
      common=unname(tools::md5sum(file.path(NUMERICAL_REPO_ROOT,"paper/numerical/common.R")))) ,
    package_manifest=F("eivgp_module_identity")(),
    dependency_lock=unname(tools::md5sum(file.path(NUMERICAL_REPO_ROOT,"renv.lock"))),
    full_map_hash=if(nzchar(cfg$full_fit_map))unname(tools::md5sum(cfg$full_fit_map))else NULL)
  manifest<-file.path(cfg$output,"protocol.rds")
  if(file.exists(manifest)&&!mixedgp_followup_identity_equal(readRDS(manifest),identity))stop("Configuration/code changed: use a new output directory")
  if(!file.exists(manifest))atomic_rds(identity,manifest)
  capture.output(sessionInfo(),file=file.path(cfg$output,"sessionInfo.txt"))
  revision<-tryCatch(system2("git",c("-C",shQuote(NUMERICAL_REPO_ROOT),"rev-parse","HEAD"),stdout=TRUE,stderr=FALSE),error=function(e)"unavailable")
  writeLines(revision,file.path(cfg$output,"git_revision.txt"))
  write.csv(jobs,file.path(cfg$output,"requested_jobs.csv"),row.names=FALSE)
  for(r in cfg$reps) {
    fp<-file.path(cfg$output,sprintf("data-rep%03d.rds",r));record<-list(file=basename(files[as.character(r)]),md5=unname(tools::md5sum(files[as.character(r)])))
    if(file.exists(fp)&&!identical(readRDS(fp),record))stop("Frozen dataset changed for replication ",r)
    if(!file.exists(fp))atomic_rds(record,fp)
  }
  if(action=="plan") {print(jobs);cat(nrow(jobs),"jobs; completed scores will be skipped. Profile:",cfg$profile,"\n");return(invisible(jobs))}
  registry<-file.path(cfg$output,"job_registry.rds")
  alljobs<-if(file.exists(registry))unique(rbind(readRDS(registry),jobs))else jobs
  if(action=="report")return(aggregate_scores(alljobs,cfg$output))
  runlock<-file.path(cfg$output,"RUN.lock")
  if(!dir.create(runlock,showWarnings=FALSE))stop("Run already locked: ",runlock)
  on.exit(unlink(runlock,recursive=TRUE),add=TRUE)
  writeLines(as.character(Sys.getpid()),file.path(runlock,"pid"))
  atomic_rds(alljobs,registry)
  allocation <- mixedgp_followup_allocation(cfg$cores,length(cfg$reps))
  timing <- mixedgp_batch_timing(file.path(cfg$output,"timing"),jobs$job,groups=jobs$method,workers=allocation$workers)
  run_rep <- function(indices) {
  cfg$cores <- mixedgp_budget_cores(allocation$chains)
  failed<-0L
  for(j in indices) {
    job<-jobs[j,];folder<-file.path(cfg$output,job$job);scorefile<-file.path(folder,"scores.rds")
    if(file.exists(scorefile)) { timing(j,"cached"); next }
    dir.create(folder,recursive=TRUE,showWarnings=FALSE)
    # Atomic lock prevents two workers writing the same job. Remove only after checking no process owns it.
    lock<-paste0(folder,".lock");if(!dir.create(lock,showWarnings=FALSE))stop("Job locked: ",lock)
    writeLines(as.character(Sys.getpid()),file.path(lock,"pid"))
    message("Running ",job$job)
    timing(j,"running"); job_started <- proc.time()[["elapsed"]]
    tryCatch({
      artifact<-readRDS(files[as.character(job$rep)]);dat<-artifact$data
      stopifnot(artifact$study=="study2",artifact$scenario=="primary",artifact$design$q==4,artifact$design$d==2)
      cal<-artifact$calibration_sets[[as.character(job$calib)]];stopifnot(length(cal)==job$calib)
      if(qr(cbind(1,dat$train$U[cal,,drop=FALSE]))$rank<3L)stop("Calibration fails affine-rank anchoring")
      Uobs<-matrix(NA_real_,nrow(dat$train$U),2);Uobs[cal,]<-dat$train$U[cal,]
      base_seed<-cfg$seed+job$rep*1000L+job$calib
      shared<-file.path(cfg$output,sprintf("shared-rep%03d-cal%02d",job$rep,job$calib));dir.create(shared,showWarnings=FALSE)
      mp<-file.path(shared,"measurement.rds")
      if(job$method=="CUT"||grepl("-CC$",job$method)) {
        if(!file.exists(mp))atomic_rds(stage_time(shared,"measurement",measurement_cut(dat$train$C,Uobs,cfg,base_seed)),mp)
        meas<-readRDS(mp)
        write.csv(meas$diagnostics,file.path(shared,"measurement_diagnostics.csv"),row.names=FALSE)
      }
      tp<-file.path(shared,"truth.rds")
      if(!file.exists(tp))atomic_rds(list(m=F("oracle_m0_quadrature_2d")(dat$test$X,dat$test$C,dat$true_params),
        f=F("f0_2d")(dat$test$X,dat$test$U)),tp)
      truth<-readRDS(tp)
      if(job$method=="CUT") {
        passes<-diagnostic_pass(meas$diagnostics);states<-list()
        for(b in seq_along(meas$states)) {
          cp<-file.path(folder,sprintf("conditional-%03d.rds",b))
          sp<-sub("[.]rds$","-states.rds",cp)
          if(!file.exists(sp)) {
            if(!file.exists(cp))atomic_rds(stage_time(folder,sprintf("conditional-%03d",b),fit_full(dat$train,meas$states[[b]]$U,cfg,base_seed+200000L+b,TRUE)),cp)
            fit<-F("mixedgp_read_artifact")(cp)
            passed<-record_gp_diagnostics(fit,sub("rds$","diagnostics.csv",cp),TRUE)
            atomic_rds(list(pass=passed,states=select_gp(fit,cfg$inner_draws,meas$states[[b]])),sp)
          }
          if(job$rep==min(cfg$reps) && b==1L) {
            fit<-F("mixedgp_read_artifact")(cp)
            rows<-F("mixedgp_study2_diagnostic_rows")(dat$test$X,dat$test$C,5L)
            series<-F("mixedgp_cut_conditional_series")(fit,meas$states[[b]],
              dat$test$X[rows,,drop=FALSE],dat$test$C[rows,,drop=FALSE],dat$test$U[rows,,drop=FALSE],
              n_latent=cfg$latent,seed=base_seed+2900L)
            diagnostic_fit<-fit
            diagnostic_fit$scientific_ess<-list(target_series=series,diagnostics=F("mixedgp_summarize_diagnostic_series")(series),
              scope="conditional_gp",success=NA)
            F("mixedgp_save_scientific_trace")(diagnostic_fit,file.path(folder,"scientific-trace-cut-conditional"),
              sprintf("CUT conditional GP | rep %d | calibration %d | outer state 1",job$rep,job$calib))
          }
          # Preserve full chain histories as well as selected integration states.
          saved<-retire_conditional_fit(cp,sp)
          passes<-c(passes,saved$pass)
          states<-c(states,saved$states)
        }
        pred<-stage_time(folder,"prediction",gp_evaluate(states,dat$train,dat$test,cfg,base_seed+300000L))
        F("mixedgp_cut_outer_trace")(pred,states,dat$test,cfg,folder,
          sprintf("CUT outer integration | rep %d | calibration %d",job$rep,job$calib))
        # Keep every measurement draw for imputation, not just integration states.
        mis<-which(!complete.cases(Uobs));arr<-meas$trace[,,seq_len(length(mis)*2),drop=FALSE]
        pred$U<-array(NA_real_,c(dim(arr)[1]*4,nrow(Uobs),2));pred$U[,mis,]<-matrix(arr,ncol=length(mis)*2)
      } else if(job$method=="FULL") {
        fp<-file.path(folder,"fit.rds")
        if(nzchar(cfg$full_fit_map)) {
          map<-read.csv(cfg$full_fit_map,stringsAsFactors=FALSE)
          if(!all(c("rep","calib","path")%in%names(map)))stop("FULL map requires rep,calib,path columns")
          selected<-map[map$rep==job$rep & map$calib==job$calib,,drop=FALSE]
          if(nrow(selected)!=1||!file.exists(selected$path))stop("FULL map must supply exactly one existing fit for this job")
          provenance<-list(path=normalizePath(selected$path),md5=unname(tools::md5sum(selected$path)))
          sourcefile<-file.path(folder,"fit_source.rds")
          if(file.exists(sourcefile)&&!identical(readRDS(sourcefile),provenance))stop("Saved FULL source changed")
          fit<-read_matched_full(selected$path,dat$train,Uobs);atomic_rds(provenance,sourcefile)
        } else {
          if(!file.exists(fp))atomic_rds(fit_full(dat$train,Uobs,cfg,base_seed+400000L),fp)
          fit<-F("mixedgp_read_artifact")(fp)
        }
        passes<-record_gp_diagnostics(fit,file.path(folder,"full_diagnostics.csv"));states<-select_gp(fit,cfg$outer*cfg$inner_draws)
        pred<-stage_time(folder,"prediction",gp_evaluate(states,dat$train,dat$test,cfg,base_seed+300000L));pred$U<-fit$mcmc$samples_U
      } else if(job$method%in%c("UC-GP","LVGP","EzGP")) {
        cp<-file.path(cfg$output,sprintf("categorical-rep%03d",job$rep),job$method,"predictions.rds")
        pred<-categorical_evaluate(job$method,dat$train,dat$test,cfg,cfg$seed+job$rep*1000L+600000L,cp)
        passes<-is.null(pred$optimization_status)||!grepl("nonconverged",pred$optimization_status)
      } else {
        up<-file.path(shared,"prospective.rds")
        if(!file.exists(up))atomic_rds(prospective_u(meas$states,dat$test$C,cfg$latent,base_seed+300000L),up)
        passes<-diagnostic_pass(meas$diagnostics)
        pred<-stage_time(folder,"cc-fit-predict",cc_evaluate(job$method,dat$train,dat$test,cal,readRDS(up),cfg,base_seed+500000L,file.path(folder,"fit.rds")))
        if(!is.null(pred$optimization_status)&&grepl("nonconverged",pred$optimization_status))passes<-c(passes,FALSE)
      }
      atomic_rds(list(pass=all(passes),checks=passes,profile=cfg$profile,ess_required=FALSE,
        diagnostic_scope=if(job$method=="CUT")"Measurement and conditional screens are descriptive, not final CUT precision certification"else"descriptive",
        optimizer_status=pred$optimization_status,
        measurement_pass=if(job$method=="CUT"||grepl("-CC$",job$method))diagnostic_pass(meas$diagnostics)else NA,
        note="MCMC screening only; also inspect target stability and nested integration error before publication"),file.path(folder,"diagnostics.rds"))
      if(!is.null(pred$m_mean)) {
        cols<-seq_len(ncol(pred$m_mean));half<-length(cols)%/%2
        stability<-data.frame(test_id=seq_len(nrow(pred$m_mean)),
          m_first=rowMeans(pred$m_mean[,cols<=half,drop=FALSE]),m_second=rowMeans(pred$m_mean[,cols>half,drop=FALSE]))
        stability$difference<-stability$m_second-stability$m_first
        write.csv(stability,file.path(folder,"integration_stability.csv"),row.names=FALSE)
      }
      scores<-stage_time(folder,"scoring",score_result(pred,dat,cal,truth,job$method))
      # Imputation summaries are in scores; full draws remain in the source fit
      # or shared measurement trace, so do not duplicate them per prediction file.
      pred$U<-NULL
      atomic_rds(pred,file.path(folder,"predictions.rds"))
      atomic_rds(scores,scorefile);unlink(file.path(folder,"failure.rds"))
      timing(j,"complete",proc.time()[["elapsed"]]-job_started)
    },error=function(e){timing(j,"failed",proc.time()[["elapsed"]]-job_started);failed<<-failed+1L;atomic_rds(list(message=conditionMessage(e),stage=if(is.null(e$stage))"unknown"else e$stage,time=Sys.time()),file.path(folder,"failure.rds"));message("FAILED: ",conditionMessage(e))},finally=unlink(lock,recursive=TRUE))
  }
  failed
  }
  failed <- sum(mixedgp_followup_apply(split(seq_len(nrow(jobs)),jobs$rep),run_rep,allocation))
  aggregate_scores(alljobs,cfg$output)
  if(failed)stop(failed," jobs failed; see status.csv. Rerun the identical command to retry.")
  invisible(jobs)
}
  environment()
}

# Study I: published predictors with continuous (X,U) on calibrated records.
# This is a response-free prediction adaptation, not the matched CUT ablation.
mixedgp_study1_cc_measurement <- function(C,U_obs,cfg,seed) {
  stopifnot(length(C)==length(U_obs),all(is.finite(U_obs[!is.na(U_obs)])))
  cal<-which(!is.na(U_obs));miss<-which(is.na(U_obs))
  disk<-mixedgp_sampler_disk(list(component="study1-measurement",C=C,U_obs=U_obs,
    iter=cfg$iter,warmup=cfg$warmup,outer=cfg$outer,seed=seed))
  if(!is.null(disk))on.exit(unlink(file.path(disk,"LOCK"),recursive=TRUE),add=TRUE)
  run<-function(ch)fit_threshold_measurement_response_free(C,U_obs,cal,6L,
    n_iter=cfg$iter,burn=cfg$warmup,seed=seed+10000L*ch,
    .checkpoint_path=if(is.null(disk))NULL else file.path(disk,sprintf("chain-%02d.rds",ch)))
  chains<-mixedgp_parallel_lapply(1:4,run,n_cores=min(cfg$cores,4L),.preserve_rng=TRUE)
  if(any(vapply(chains,inherits,logical(1),"try-error")))stop("Study I measurement chain failed")
  trace<-array(NA_real_,c(cfg$iter-cfg$warmup,4L,5L+length(miss)),
    dimnames=list(NULL,paste0("chain",1:4),c(paste0("tau",1:5),paste0("u",miss))))
  for(ch in 1:4)trace[,ch,]<-cbind(chains[[ch]]$samples_tau,chains[[ch]]$samples_u[,miss,drop=FALSE])
  diagnostics<-posterior::summarise_draws(posterior::as_draws_array(trace))
  ids<-floor(((seq_len(cfg$outer/4)-.5)/(cfg$outer/4))*(cfg$iter-cfg$warmup))+1L
  states<-unlist(lapply(ids,function(i)lapply(1:4,function(ch)
    list(tau=chains[[ch]]$samples_tau[i,],chain=ch,iteration=cfg$warmup+i))),recursive=FALSE)
  result<-list(chains=chains,states=states,diagnostics=diagnostics,
    target="Response-free deterministic threshold model; U ~ N(0,1), ordered uniform thresholds on [-8,8]. Same measurement adapter as Study I PI/CC references; not a matched FULL-CUT ablation.")
  if(!is.null(disk)&&length(list.files(disk,"^chain-[0-9]+[.]rds$")))
    attr(result,"eivgp_study1_restart_storage")<-mixedgp_checkpoint_reference(disk)
  result
}

mixedgp_study1_cc_workflow <- function(args,repo_root,root,data_root,action) {
  release_budget <- mixedgp_budget_enter(as.integer(if(is.null(args$cores))4L else args$cores))
  on.exit(release_budget(),add=TRUE)
  api<-eivgp_numerical_api(repo_root)
  get<-function(k,d)if(is.null(args[[k]]))d else args[[k]]
  profile<-get("profile","paper");pilot<-profile=="pilot"
  a<-args[setdiff(names(args),c("component","competitors","checkpoint-seconds","checkpoint-iterations","checkpoint-after-seconds","eta"))]
  a$output<-root;a[["data-root"]]<-data_root
  if(is.null(a$methods))a$methods<-if(pilot)"UC-GP,LVGP,EzGP,UC-GP-CC,LVGP-CC,EzGP-CC"else"UC-GP-CC,LVGP-CC,EzGP-CC"
  if(is.null(a$reps)&&pilot)a$reps<-"1"
  if(is.null(a$calib)&&pilot)a$calib<-"50"
  cfg<-api$settings(a);cfg$m_vec<-6L;cfg$controls<-mixedgp_competitor_protocol("study1")
  if(any(!cfg$methods%in%c("UC-GP","LVGP","EzGP","UC-GP-CC","LVGP-CC","EzGP-CC")))
    stop("Study I complete-case component supports published predictors only")
  if(any(names(args)%in%c("full-fit-map","inner-iter","inner-warmup","inner-draws")))
    stop("Study I complete-case component has no conditional response MCMC or FULL map")
  eta<-as.integer(strsplit(get("eta","0,1"),",",fixed=TRUE)[[1]])
  if(anyNA(eta)||!length(eta)||anyDuplicated(eta)||any(!eta%in%0:1))stop("eta must select 0, 1 or 0,1")
  jobs<-expand.grid(rep=cfg$reps,calib=cfg$calib,method=cfg$methods,eta=eta,stringsAsFactors=FALSE)
  jobs$job<-sprintf("rep%03d-cal%02d-%s",jobs$rep,jobs$calib,jobs$method)
  if(action=="plan") {print(jobs);cat(nrow(jobs),"Study I jobs; no fits started. Paper default adds 600 CC fits; pilot pairs original/CC versions on replication 1, both regimes, calibration 50.\n");return(invisible(jobs))}
  if(!action%in%c("data","run","report"))stop("Use action=run for complete-case fitting, or component=core for fit-only")
  dir.create(root,recursive=TRUE,showWarnings=FALSE)
  if(action=="run") {
    lock<-file.path(root,"RUN.lock");if(!dir.create(lock,showWarnings=FALSE))stop("Study I complete-case workflow locked: ",lock)
    writeLines(as.character(Sys.getpid()),file.path(lock,"pid"));on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  }
  result<-list();failures<-0L
  for(e in eta) {
    cell<-sprintf("eta%d_balanced",e);out<-file.path(root,cell);dir.create(out,recursive=TRUE,showWarnings=FALSE)
    selected<-jobs[jobs$eta==e,setdiff(names(jobs),"eta"),drop=FALSE]
    files<-setNames(file.path(data_root,"complete-case",cell,sprintf("rep%03d.rds",cfg$reps)),cfg$reps)
    for(r in cfg$reps) {
      expected<-make_study1_synthetic_dataset(rep_id=r,n=100L,n_test=100L,m=6L,
        scenario="heterogeneity_continuum",threshold_design="balanced",heterogeneity_eta=e,calib_grid=c(0L,20L,50L))
      p<-files[as.character(r)]
      if(!file.exists(p)&&action%in%c("data","run"))api$atomic_rds(expected,p)
      if(!file.exists(p)||!identical(readRDS(p),expected))stop("Missing or altered frozen Study I data: ",p)
    }
    if(action=="data")next
    identity<-list(settings=cfg[setdiff(names(cfg),c("reps","calib","methods","cores","output","data_root"))],
      eta=e,package_manifest=eivgp_module_identity(),R=as.character(getRversion()),
      versions=setNames(vapply(c("kergp","LVGP","EzGP","posterior"),function(p)as.character(packageVersion(p)),""),c("kergp","LVGP","EzGP","posterior")))
    protocol<-file.path(out,"protocol.rds")
    if(file.exists(protocol)&&!mixedgp_followup_identity_equal(readRDS(protocol),identity))stop("Configuration/code changed; use a new Study I output directory")
    if(!file.exists(protocol))api$atomic_rds(identity,protocol)
    registry<-file.path(out,"job_registry.rds");alljobs<-if(file.exists(registry))unique(rbind(readRDS(registry),selected))else selected
    if(action=="run")api$atomic_rds(alljobs,registry)
    allocation <- mixedgp_followup_allocation(cfg$cores,length(cfg$reps))
    run_rep <- function(indices) {
    cfg$cores <- mixedgp_budget_cores(allocation$chains)
    failures <- 0L
    for(j in indices) {
      job<-selected[j,];folder<-file.path(out,job$job);dir.create(folder,recursive=TRUE,showWarnings=FALSE)
      scorefile<-file.path(folder,"scores.rds");if(file.exists(scorefile))next
      message("Running Study I ",cell," ",job$job)
      tryCatch({
        artifact<-readRDS(files[as.character(job$rep)]);d<-artifact$data
        convert<-function(d)list(X=matrix(d$x,ncol=1L),U=matrix(d$u,ncol=1L),C=matrix(d$c,ncol=1L),y=d$y)
        dat<-list(train=convert(d$train),test=convert(d$test));cal<-artifact$calibration_sets[[as.character(job$calib)]]
        stopifnot(length(cal)==job$calib,qr(cbind(1,dat$train$U[cal,,drop=FALSE]))$rank==2L)
        seed<-cfg$seed+job$rep*1000L+job$calib
        if(grepl("-CC$",job$method)) {
          # No hidden training U or any Y enters the measurement likelihood.
          Uobs<-rep(NA_real_,100L);Uobs[cal]<-dat$train$U[cal,1]
          shared<-file.path(out,sprintf("shared-rep%03d-cal%02d",job$rep,job$calib));dir.create(shared,showWarnings=FALSE)
          mp<-file.path(shared,"measurement.rds")
          if(!file.exists(mp))api$atomic_rds(api$stage_time(shared,"measurement",mixedgp_study1_cc_measurement(d$train$c,Uobs,cfg,seed)),mp)
          meas<-readRDS(mp);write.csv(meas$diagnostics,file.path(shared,"measurement_diagnostics.csv"),row.names=FALSE)
          up<-file.path(shared,"prospective.rds")
          if(!file.exists(up)) {
            set.seed(seed+300000L)
            us<-lapply(meas$states,function(st){z<-matrix(NA_real_,cfg$latent,100L)
              for(i in 1:100)z[,i]<-rtruncnorm_vec(rep(0,cfg$latent),rep(1,cfg$latent),
                rep(c(-Inf,st$tau)[d$test$c[i]],cfg$latent),rep(c(st$tau,Inf)[d$test$c[i]],cfg$latent))
              array(z,c(cfg$latent,100L,1L))})
            api$atomic_rds(us,up)
          }
          pred<-api$stage_time(folder,"cc-fit-predict",api$cc_evaluate(job$method,dat$train,dat$test,cal,readRDS(up),cfg,seed+500000L,file.path(folder,"fit.rds")))
          passes<-api$diagnostic_pass(meas$diagnostics)
        } else {
          cp<-file.path(out,sprintf("original-rep%03d",job$rep),job$method,"predictions.rds")
          pred<-api$stage_time(folder,"original-fit-predict",api$categorical_evaluate(job$method,dat$train,dat$test,cfg,250000L+1000L*job$rep,cp))
          passes<-TRUE
        }
        passes<-passes&&(is.null(pred$optimization_status)||!grepl("nonconverged",pred$optimization_status))
        truth<-list(m=m0_1d(d$test$x,d$test$c,tau=d$tau_true,scenario="heterogeneity_continuum"))
        scores<-api$stage_time(folder,"scoring",api$score_result(pred,dat,cal,truth,job$method))
        api$atomic_rds(pred,file.path(folder,"predictions.rds"));api$atomic_rds(scores,scorefile)
        api$atomic_rds(list(pass=passes,optimizer_status=pred$optimization_status,
          measurement_pass=if(grepl("-CC$",job$method))api$diagnostic_pass(meas$diagnostics)else NA),file.path(folder,"diagnostics.rds"));unlink(file.path(folder,"failure.rds"))
      },error=function(err){failures<<-failures+1L;api$atomic_rds(list(message=conditionMessage(err),stage=if(is.null(err$stage))"unknown"else err$stage,time=Sys.time()),file.path(folder,"failure.rds"));message("FAILED: ",conditionMessage(err))})
    }
    failures
    }
    if(action=="run") failures <- failures + sum(mixedgp_followup_apply(
      split(seq_len(nrow(selected)),selected$rep),run_rep,allocation))
    result[[cell]]<-api$aggregate_scores(alljobs,out)
  }
  if(failures)stop(failures," Study I jobs failed; inspect per-regime status.csv; identical reruns retry failures")
  invisible(result)
}

# Permanent orchestration for both paper studies. Historical readers are separate.
run_paper_study <- function(study=c("study1","study2"), args=list(), repo_root=getwd()) {
  old_clock<-options(eivGP.checkpoint_script_start=getOption("eivGP.checkpoint_script_start",proc.time()[["elapsed"]]))
  on.exit(options(old_clock),add=TRUE)
  release_budget <- mixedgp_budget_enter(as.integer(if(is.null(args$cores)) 4L else args$cores))
  on.exit(release_budget(),add=TRUE)
  study<-match.arg(study)
  # Follow-up mode adds missing comparisons; it never routes through core fitting.
  if(identical(args$component,"follow-up")) {
    args$component <- "publication"
    if(is.null(args$stage)) args$stage <- "full"
    return(mixedgp_publication_workflow(study,args,repo_root))
  }
  if(identical(args$component,"publication")) return(mixedgp_publication_workflow(study,args,repo_root))
  api<-eivgp_numerical_api(repo_root)
  get<-function(k,d)if(is.null(args[[k]]))d else args[[k]]
  mixing <- get("mixing-plots", TRUE)
  if(length(mixing)!=1L || is.na(mixing) || !tolower(as.character(mixing)) %in% c("true","false"))
    stop("mixing-plots must be true or false")
  mixing <- identical(tolower(as.character(mixing)),"true")
  args[["mixing-plots"]] <- NULL
  component<-get("component","all")
  if(!component%in%c("all","core","core-competitors","mechanism","complete-case"))stop("component must be all, core, core-competitors, mechanism or complete-case")
  if(study=="study1" && component=="mechanism")stop("Use complete-case for Study I; matched CUT belongs to Study II")
  if(study=="study2" && component=="complete-case")stop("Use mechanism with methods=UC-GP-CC,LVGP-CC,EzGP-CC for Study II")
  action<-get("action","plan");profile<-get("profile","paper")
  if(!action%in%c("plan","data","fit","run","report"))stop("Unknown action")
  if(action=="fit" && study=="study2" && component!="core")stop("Use action=run for the complete Study II workflow, or component=core for fit-only")
  if(action=="fit" && study=="study1" && component!="core")stop("Use action=run for Study I with complete-case additions, or component=core for fit-only")
  if(profile=="pilot" && !((study=="study2"&&component=="mechanism")||(study=="study1"&&component=="complete-case")))
    stop("Pilot requires Study II mechanism or Study I complete-case")
  root<-get("output",Sys.getenv("EIVGP_RESULTS_ROOT",api$numerical_default_output(study)))
  data<-get("data-root",Sys.getenv("EIVGP_DATA_ROOT",api$numerical_default_data(study)))
  if(!nzchar(root)||!nzchar(data)||grepl("Dropbox|Overleaf",root,ignore.case=TRUE))stop("Supply a non-synced local output root and data root")
  if(action!="plan" && grepl("Dropbox|Overleaf",data,ignore.case=TRUE))stop("Copy frozen inputs to a local data-root before running")
  root <- path.expand(root); data <- path.expand(data)
  checkpoint_seconds <- as.numeric(get("checkpoint-seconds","10800"))
  checkpoint_iterations <- as.numeric(get("checkpoint-iterations","Inf"))
  if(length(checkpoint_seconds)!=1L || is.na(checkpoint_seconds) || checkpoint_seconds<=0 ||
     length(checkpoint_iterations)!=1L || is.na(checkpoint_iterations) || checkpoint_iterations<1 ||
     checkpoint_iterations!=floor(checkpoint_iterations)) stop("Checkpoint intervals must be positive or Inf (finite iterations must be integral)")
  checkpoint_after <- as.numeric(get("checkpoint-after-seconds",
    if(any(c("checkpoint-seconds","checkpoint-iterations","checkpoint-after-seconds") %in% names(args)))"0"else"21600"))
  if(length(checkpoint_after)!=1L || !is.finite(checkpoint_after) || checkpoint_after<0)
    stop("checkpoint-after-seconds must be finite and nonnegative")
  old_checkpoint <- options(eivGP.checkpoint_root=file.path(root,"checkpoints"),
    eivGP.checkpoint_seconds=checkpoint_seconds,eivGP.checkpoint_iterations=checkpoint_iterations,
    eivGP.checkpoint_after_seconds=checkpoint_after)
  on.exit(options(old_checkpoint),add=TRUE)
  if(action %in% c("fit","run")) {
    if(is.infinite(checkpoint_seconds)&&is.infinite(checkpoint_iterations))
      message("Outputs: ",root,"; completed fits only; interrupted fits restart from the beginning.")
    else message("Outputs: ",root,"; checkpoints every ",checkpoint_seconds,
      " seconds or ",checkpoint_iterations," iterations, after ",checkpoint_after,
      " seconds of script runtime; completed fits always saved immediately.")
  }
  common<-c("component","action","profile","output","data-root","cores","competitors",
    "checkpoint-seconds","checkpoint-iterations","checkpoint-after-seconds")
  targeted<-c("reps","calib","methods","outer","inner-draws","latent","iter","warmup","inner-iter","inner-warmup","seed","full-fit-map","eta","gp-response-latent","cc-response-latent")
  if(length(setdiff(names(args),c(common,targeted))))stop("Unknown options: ",paste(setdiff(names(args),c(common,targeted)),collapse=", "))
  if(!component%in%c("mechanism","complete-case") && any(names(args)%in%targeted))stop("Replication/method/sampler selections require component=mechanism; the core grid uses its fixed profile")
  if(component%in%c("mechanism","complete-case") && "competitors"%in%names(args))stop("Use methods for the mechanism component")
  result<-list()
  # Independent original competitors use only frozen data, never joint EIV fits.
  if(component=="core-competitors") {
    if(action=="fit")stop("Use action=run for original competitors")
    config<-study_config(study,profile,output_root=file.path(root,"core"),data_root=data,cores=as.integer(get("cores","4")))
    methods<-strsplit(get("competitors","UC-GP,LVGP,EzGP"),",",fixed=TRUE)[[1]]
    if(!length(methods)||anyDuplicated(methods)||any(!methods%in%c("UC-GP","LVGP","EzGP")))stop("Invalid original competitors")
    if(action=="plan")return(run_numerical_study(config,"plan",competitors=character()))
    if(action%in%c("data","run")) {
      prepared<-run_numerical_study(config,"data",competitors=character())
      if(action=="data")return(prepared)
    }
    return(mixedgp_run_competitor_workflow(config,file.path(root,"original-competitors"),
      workers=as.integer(get("cores","4")),methods=methods,action=if(action=="report")"export"else"run"))
  }
  if(component%in%c("all","core")) {
    config<-study_config(study,profile,output_root=file.path(root,"core"),data_root=data,cores=as.integer(get("cores","4")))
    competitors<-strsplit(get("competitors",if(profile=="quick")"none"else"UC-GP,LVGP,EzGP"),",",fixed=TRUE)[[1]]
    if(identical(competitors,"none"))competitors<-character()
    result$core<-run_numerical_study(config,action,competitors=competitors)
  }
  if(study=="study2" && component%in%c("all","mechanism")) {
    a<-args[setdiff(names(args),c("component","competitors","checkpoint-seconds","checkpoint-iterations","checkpoint-after-seconds"))]
    a$output<-file.path(root,"mechanism");a[["data-root"]]<-data
    result$mechanism<-api$run_mechanism(api$settings(a),action)
  }
  if(study=="study1" && component%in%c("all","complete-case")) {
    if(action=="fit")stop("Use action=run for the complete Study I workflow, or component=core for fit-only")
    cc_args<-args
    if(profile=="paper")cc_args<-modifyList(list(iter="12000",warmup="4000",outer="16",
      latent="128",`cc-response-latent`="64"),args)
    result$complete_case<-mixedgp_study1_cc_workflow(cc_args,repo_root,file.path(root,"complete-case"),data,action)
  }
  if(mixing && action %in% c("report","run")) {
    directories <- c(if(component %in% c("all","core"))file.path(root,"core"),
      if(study=="study1" && component %in% c("all","complete-case"))file.path(root,"complete-case"),
      if(study=="study2" && component %in% c("all","mechanism"))file.path(root,"mechanism"))
    attr(result,"mixing_plots") <- unlist(lapply(directories,function(p) {
      if(!dir.exists(p)) return(character())
      source <- p
      if(identical(p,file.path(root,"core")) && !is.null(result$core$run_dir)) {
        cells <- file.path(result$core$run_dir,"cells")
        if(dir.exists(cells)) source <- cells
      }
      save_mixing_plots(source, output_dir=file.path(p,"mixing-plots"))
    }),use.names=FALSE)
  }
  invisible(result)
}

# Restart files are committed by rename on the same filesystem. Never truncate
# a valid checkpoint in place. Timing is operational metadata, not MCMC evidence.
mixedgp_checkpoint_write <- function(object, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  tmp <- tempfile(".checkpoint-", tmpdir=dirname(path))
  on.exit(unlink(tmp), add=TRUE)
  saveRDS(object, tmp, compress="gzip", version=3L)
  if (!file.rename(tmp, path)) stop("Cannot publish checkpoint; prior file retained: ", path)
  invisible(path)
}

mixedgp_timing_record <- function(path, row) {
  # Each chain/task owns its own history; no concurrent append to a shared file.
  history <- sub("[.]csv$", "-history.csv", path)
  old <- if(file.exists(history)) utils::read.csv(history, stringsAsFactors=FALSE) else NULL
  mixedgp_atomic_write_csv(rbind(old, row), history)
  mixedgp_atomic_write_csv(row, path)
}

mixedgp_sampler_disk <- function(identity) {
  root <- getOption("eivGP.checkpoint_root", NULL)
  if (is.null(root)) return(NULL)
  identity <- list(schema=1L, model=identity, modules=mixedgp_followup_checkpoint_modules(),
    R=R.version.string, platform=R.version$platform, rng_kind=RNGkind())
  tmp <- tempfile(); on.exit(unlink(tmp), add=TRUE)
  saveRDS(identity, tmp, compress=FALSE, version=3L)
  directory <- file.path(root, unname(tools::md5sum(tmp)))
  dir.create(directory, recursive=TRUE, showWarnings=FALSE)
  lock <- file.path(directory, "LOCK")
  if (!dir.create(lock, showWarnings=FALSE)) stop(
    "Checkpoint is locked: ", lock,
    ". Check owner.csv; remove this lock directory only after its process has stopped.")
  ok <- FALSE
  on.exit(if(!ok) unlink(lock,recursive=TRUE), add=TRUE)
  utils::write.csv(data.frame(pid=Sys.getpid(),host=Sys.info()[["nodename"]],
    started_utc=format(Sys.time(),tz="UTC",usetz=TRUE)), file.path(lock,"owner.csv"),row.names=FALSE)
  manifest <- file.path(directory,"identity.rds")
  if(file.exists(manifest)) {
    if(!identical(readRDS(manifest),identity)) stop("Checkpoint identity mismatch: ",directory)
  } else mixedgp_checkpoint_write(identity,manifest)
  ok <- TRUE
  directory
}

# A batch forecast uses observed completed-task wall times. Separate groups
# prevent fast comparator fits from setting the forecast for expensive CUT jobs.
mixedgp_batch_timing <- function(directory, ids, groups=rep("replication",length(ids)), workers=1L) {
  dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  paths <- file.path(directory, paste0(make.names(ids,unique=TRUE),".rds"))
  update <- function(i,status,seconds=NA_real_) {
    previous <- if(file.exists(paths[i]))readRDS(paths[i])else NULL
    duration <- if(status %in% c("cached","running") && !is.null(previous))previous$duration_seconds else if(status=="cached")NA_real_ else seconds
    attempts <- if(is.null(previous))0 else previous$total_seconds
    if(status %in% c("complete","failed") && is.finite(seconds))attempts <- attempts+seconds
    row <- list(id=ids[i],group=groups[i],status=status,duration_seconds=duration,
      total_seconds=attempts,updated_at=format(Sys.time(),tz="UTC",usetz=TRUE))
    mixedgp_checkpoint_write(row,paths[i])
    rows <- lapply(paths,function(p)if(file.exists(p))readRDS(p)else NULL)
    pending <- vapply(rows,function(x)is.null(x)||!x$status%in%c("complete","cached"),logical(1))
    estimate <- numeric(length(ids)); observed <- integer(length(ids))
    for(g in unique(groups)) {
      finished <- Filter(function(x)!is.null(x)&&x$group==g&&x$status%in%c("complete","cached")&&
        is.finite(x$duration_seconds)&&x$duration_seconds>0,rows)
      rate <- if(length(finished))mean(vapply(finished,`[[`,0,"duration_seconds"))else NA_real_
      estimate[groups==g & pending] <- rate
      observed[groups==g] <- length(finished)
    }
    remaining <- if(anyNA(estimate))NA_real_ else sum(estimate)/max(1L,workers)
    # Each worker owns a forecast snapshot: concurrent workers never overwrite history.
    mixedgp_timing_record(file.path(directory,paste0("forecast-",Sys.getpid(),".csv")),data.frame(
      updated_at=row$updated_at, completed=sum(!pending),total=length(ids),
      observed_completed_tasks=sum(vapply(rows,function(x)!is.null(x)&&is.finite(x$duration_seconds)&&
        x$status%in%c("complete","cached"),logical(1))),
      recorded_task_seconds=sum(vapply(rows,function(x)if(is.null(x))0 else x$total_seconds,0)),
      assumed_workers=workers,remaining_seconds=remaining,
      estimated_finish_utc=if(is.finite(remaining))format(Sys.time()+remaining,tz="UTC",usetz=TRUE)else NA_character_))
    invisible(row)
  }
  update
}

# Publication completion workflow: orchestration only; inference remains in the
# authoritative study components. Planning never fits or reads chain archives.
mixedgp_publication_workflow <- function(study, args, repo_root) {
  get <- function(k,d) if(is.null(args[[k]])) d else args[[k]]
  allowed <- c("component","stage","action","profile","output","data-root","reps","calib","methods","cores","full-fit-root",
    "full-fit-map","outer","latent","inner-draws","iter","warmup","inner-iter",
    "inner-warmup","gp-response-latent","cc-response-latent","mixing-plots",
    "checkpoint-seconds","checkpoint-iterations","checkpoint-after-seconds")
  if(length(setdiff(names(args),allowed))) stop("Unknown publication options: ",paste(setdiff(names(args),allowed),collapse=", "))
  if(!identical(get("profile","paper"),"paper"))stop("Publication follow-ups use profile=paper; stage selects the pilot size")
  stage <- match.arg(get("stage","pilot"),c("pilot","integration","full"))
  action <- match.arg(get("action","plan"),c("plan","data","run","report"))
  parse_reps <- function(x) {
    tokens<-strsplit(x,",",fixed=TRUE)[[1]]
    if(!length(tokens)||any(!grepl("^[0-9]+(:[0-9]+)?$",tokens)))stop("Invalid reps")
    z<-unlist(lapply(tokens,function(t){v<-as.integer(strsplit(t,":",fixed=TRUE)[[1]]);if(length(v)==2L)seq.int(v[1],v[2])else v}))
    if(anyNA(z)||any(z<1L|z>50L)||anyDuplicated(z))stop("reps must be distinct indices in 1:50")
    sort(z)
  }
  reps<-parse_reps(get("reps",if(stage=="full")"1:50" else "1,25"))
  root<-path.expand(get("output",Sys.getenv("EIVGP_RESULTS_ROOT",file.path(path.expand("~/eivGP-results"),study))))
  data<-path.expand(get("data-root",Sys.getenv("EIVGP_DATA_ROOT",file.path(path.expand("~/eivGP-results"),"data",study))))
  if(!nzchar(root)||!nzchar(data)||any(grepl("Dropbox|Overleaf",c(root,data),ignore.case=TRUE)))stop("Supply non-synced output and data-root")
  cores<-as.integer(get("cores","4"));if(length(cores)!=1L||is.na(cores)||cores<1L)stop("cores must be positive")
  defaults<-if(study=="study1")list(iter="12000",warmup="4000",outer="16",latent="128",`cc-response-latent`="64")else
    list(iter="100000",warmup="20000",outer="32",latent="128",`inner-iter`="4000",`inner-warmup`="1000",`inner-draws`="4",`gp-response-latent`="16",`cc-response-latent`="64")
  for(k in intersect(names(defaults),names(args)))defaults[[k]]<-args[[k]]
  for(k in names(defaults)) {
    v<-suppressWarnings(as.numeric(defaults[[k]]))
    if(length(v)!=1L||!is.finite(v)||v<1||v!=floor(v))stop("Invalid ",k)
  }
  if(as.integer(defaults$iter)<=as.integer(defaults$warmup)||as.integer(defaults$outer)%%4L!=0L)stop("Invalid iterations or outer count")
  if(as.integer(defaults[["cc-response-latent"]])>as.integer(defaults$latent))stop("Response budget exceeds latent budget")
  if(study=="study1"&&any(c("inner-iter","inner-warmup","inner-draws","full-fit-map","full-fit-root","gp-response-latent")%in%names(args)))stop("Study I does not use conditional response chains or a FULL map")
  methods<-if(study=="study1")if(stage=="full")"UC-GP-CC,LVGP-CC,EzGP-CC"else"UC-GP,LVGP,EzGP,UC-GP-CC,LVGP-CC,EzGP-CC" else "FULL,CUT,UC-GP-CC,LVGP-CC,EzGP-CC"
  available <- strsplit(methods,",",fixed=TRUE)[[1]]
  selected <- strsplit(get("methods",methods),",",fixed=TRUE)[[1]]
  if(!length(selected)||any(!nzchar(selected))||anyDuplicated(selected)||any(!selected%in%available))
    stop("methods must be a distinct subset of: ",paste(available,collapse=", "))
  methods <- paste(selected,collapse=",")
  calib <- parse_reps(get("calib","20,50"))
  if(any(!calib%in%c(20L,50L)))stop("Publication calib must be 20, 50, or 20,50")
  base<-c(list(component=if(study=="study1")"complete-case"else"mechanism",action=action,profile="paper",reps=paste(reps,collapse=","),calib=paste(calib,collapse=","),methods=methods,cores=as.character(cores),`data-root`=data),defaults)
  for(k in intersect(c("full-fit-map","mixing-plots","checkpoint-seconds","checkpoint-iterations","checkpoint-after-seconds"),names(args)))base[[k]]<-args[[k]]
  variants<-list(reference=list())
  if(stage=="integration") {
    variants$outer_doubled<-list(outer=as.character(2L*as.integer(defaults$outer)))
    variants$latent_doubled<-list(latent=as.character(2L*as.integer(defaults$latent)))
    variants$response_doubled<-list(`cc-response-latent`=as.character(min(as.integer(defaults$latent),2L*as.integer(defaults[["cc-response-latent"]]))))
    if(study=="study2") {
      variants$response_doubled[["gp-response-latent"]]<-as.character(min(as.integer(defaults$latent),2L*as.integer(defaults[["gp-response-latent"]])))
      variants$gp_states_doubled<-list(`inner-draws`=as.character(2L*as.integer(defaults[["inner-draws"]])))
    }
  }
  # Separate roots prevent cross-machine locks and incompatible-cache reuse.
  batch<-paste0("reps-",paste(sprintf("%02d",reps),collapse="-"))
  plans<-lapply(names(variants),function(v){a<-modifyList(base,variants[[v]]);a$output<-file.path(root,"publication",stage,batch,v);a})
  names(plans)<-names(variants)
  jobs<-do.call(rbind,lapply(names(plans),function(v){a<-plans[[v]];d<-expand.grid(rep=reps,calib=calib,method=strsplit(methods,",",fixed=TRUE)[[1]],eta=if(study=="study1")c(0L,1L)else NA_integer_,stringsAsFactors=FALSE);d$study<-study;d$stage<-stage;d$variant<-v;d$output<-a$output;d}))
  message(nrow(jobs)," method evaluations; settings are provisional, not a convergence claim.")
  if(action=="plan") {
    print(jobs[,c("study","stage","variant","rep","calib","method","eta")],row.names=FALSE)
    return(invisible(list(jobs=jobs,commands=plans)))
  }
  # FULL must rescore matching archived fits, never silently launch replacements.
  if(study=="study2"&&"FULL"%in%selected&&action%in%c("run","report")) {
    mapfile<-get("full-fit-map","")
    fitroot<-get("full-fit-root","")
    if(nzchar(mapfile)&&nzchar(fitroot))stop("Choose full-fit-map or full-fit-root, not both")
    if(nzchar(fitroot)) {
      fitroot<-path.expand(fitroot)
      if(!dir.exists(fitroot))stop("FULL archive directory does not exist: ",fitroot)
      candidates<-list.files(fitroot,recursive=TRUE,full.names=TRUE)
      map<-expand.grid(rep=reps,calib=calib)
      map$path<-vapply(seq_len(nrow(map)),function(i) {
        name<-sprintf("recovery_fit_primary_rep%03d_cal%03d_s2v18-q4-50fcae689bed9194.rds",map$rep[i],map$calib[i])
        found<-candidates[basename(candidates)==name]
        if(length(found)!=1L)stop("Need exactly one saved FULL archive named ",name,"; found ",length(found))
        normalizePath(found,mustWork=TRUE)
      },character(1))
      mapdir<-file.path(root,"publication",stage,batch)
      dir.create(mapdir,recursive=TRUE,showWarnings=FALSE)
      mapfile<-file.path(mapdir,"full-fit-map.csv")
      write.csv(map,mapfile,row.names=FALSE)
      for(v in names(plans))plans[[v]][["full-fit-map"]]<-mapfile
    }
    if(!nzchar(mapfile)||!file.exists(mapfile))stop("Supply full-fit-map for matched saved FULL rescoring")
    map<-read.csv(mapfile,stringsAsFactors=FALSE)
    if(!all(c("rep","calib","path")%in%names(map)))stop("FULL map requires rep,calib,path")
    for(rr in reps)for(cc in calib) {
      z<-map[map$rep==rr&map$calib==cc,,drop=FALSE]
      if(nrow(z)!=1L||!file.exists(z$path))stop("Missing unique FULL archive for rep ",rr," calibration ",cc)
    }
  }
  out<-file.path(root,"publication",stage,batch);dir.create(out,recursive=TRUE,showWarnings=FALSE)
  write.csv(jobs,file.path(out,"planned_jobs.csv"),row.names=FALSE)
  results<-lapply(plans,function(a)run_paper_study(study,a,repo_root))
  invisible(list(jobs=jobs,results=results))
}


# Independent replications own their caches; chains share each worker's budget.
mixedgp_followup_allocation <- function(cores, replications) {
  stopifnot(length(cores)==1L,is.finite(cores),cores>=1L,cores==floor(cores),
    length(replications)==1L,replications>=1L)
  cores <- mixedgp_resolve_cores(cores)
  workers <- if(.Platform$OS.type=="unix") min(replications,max(1L,cores%/%4L)) else 1L
  list(workers=as.integer(workers),chains=as.integer(min(4L,cores%/%workers)),cores=cores)
}

mixedgp_followup_apply <- function(groups, fun, allocation) {
  message("CPU budget ",allocation$cores,": ",allocation$workers,
    " independent replication workers, up to ",allocation$chains," chains per worker.")
  run <- function(x) {
    # Explicit per-job seeds/checkpoint RNG determine the scientific random streams.
    fun(x)
  }
  release_budget <- mixedgp_budget_enter(allocation$cores)
  on.exit(release_budget(),add=TRUE)
  ans <- mixedgp_parallel_lapply(groups,run,n_cores=allocation$workers,
    mc.preschedule=FALSE,.preserve_rng=TRUE)
  if(any(vapply(ans,function(x)inherits(x,"try-error")||is.null(x),logical(1))))
    stop("A replication worker failed or was interrupted; retained checkpoints can be resumed.")
  unlist(ans,use.names=FALSE)
}

mixedgp_followup_identity_equal <- function(saved,current,compat=NULL) {
  if(identical(saved,current))return(TRUE)
  # Whitelist a verified orchestration-only edit, never arbitrary inference changes.
  if(is.null(compat)) {
    path <- system.file("followup-scheduler-compatibility.csv",package="eivGP")
    if(!nzchar(path))return(FALSE)
    compat <- utils::read.csv(path,stringsAsFactors=FALSE)
  }
  normalize <- function(x) {
    for(i in seq_len(nrow(compat))) {
      row <- compat[i,]
      if(row$field=="helpers.R" && identical(unname(x$package_manifest["helpers.R"]),row$new))
        x$package_manifest["helpers.R"] <- row$old
      if(row$field=="mechanism" && identical(unname(x$code["mechanism"]),row$new))
        x$code["mechanism"] <- row$old
    }
    x
  }
  identical(normalize(saved),normalize(current))
}


# Preserve checkpoint directory keys across the audited scheduling-only edit.
mixedgp_followup_checkpoint_modules <- function() {
  modules <- eivgp_module_identity()
  path <- system.file("followup-scheduler-compatibility.csv",package="eivGP")
  if(nzchar(path)) {
    compat <- utils::read.csv(path,stringsAsFactors=FALSE)
    row <- compat[compat$field=="helpers.R",,drop=FALSE]
    if(nrow(row)==1L && identical(unname(modules["helpers.R"]),row$new))
      modules["helpers.R"] <- row$old
  }
  modules
}


# Prediction stream schema 1: fixed blocks of eight selected posterior states,
# independent of worker count. Joint prediction rows stay together within a
# state; never partition a multivariate Gaussian draw across prediction rows.
mixedgp_predict_blocks <- function(args,FUN=predict_eivgp,object_name="object") {
  object <- args[[object_name]]
  ids <- args$draw_ids
  if(is.null(ids)) ids <- seq_len(nrow(object$mcmc$samples_tau))
  ids <- mixedgp_as_integer_strict(ids,"draw_ids",1L)
  if(any(ids>nrow(object$mcmc$samples_tau)))stop("draw_ids must index saved posterior draws.")
  old_rng <- mixedgp_rng_state()
  if(!is.null(args$seed)) {
    seed <- mixedgp_as_integer_strict(args$seed,"seed",0L,1L)
    on.exit(mixedgp_restore_rng_state(old_rng),add=TRUE)
  } else seed <- sample.int(.Machine$integer.max,1L)
  groups <- split(seq_along(ids),ceiling(seq_along(ids)/8L))
  # Task seeds are tied to block position, not PID or assigned CPU.
  seeds <- as.integer((as.double(seed)+104729*(seq_along(groups)-1)) %% .Machine$integer.max)
  parts <- mixedgp_parallel_lapply(seq_along(groups),function(i) {
    old <- options(eivGP.prediction_block=TRUE)
    on.exit(options(old),add=TRUE)
    a <- args;a$draw_ids <- ids[groups[[i]]]
    if("seed" %in% names(formals(FUN)))a$seed <- seeds[i]
    do.call(FUN,a)
  },n_cores=mixedgp_budget_cores(),seeds=seeds,mc.preschedule=TRUE)
  ans <- do.call(rbind,parts)
  attrs <- setdiff(names(attributes(parts[[1L]])),c("dim","dimnames"))
  for(name in attrs) {
    values <- lapply(parts,attr,which=name,exact=TRUE)
    if(name %in% c("conditional_means","conditional_vars")) value <- do.call(rbind,values)
    else if(name %in% c("rejection_acceptance","rejection_candidate_draws")) value <- unlist(values,use.names=FALSE)
    else {
      if(!all(vapply(values,identical,logical(1),values[[1L]])))
        stop("Prediction block metadata disagree: ",name)
      value <- values[[1L]]
    }
    attr(ans,name) <- value
  }
  attr(ans,"prediction_stream_schema") <- 1L
  ans
}

# Ocean storage schema 1. Model identities and RNG streams are independent of
# this physical layout. Unmarked roots retain their historical layout.
ocean_storage_flat <- function(root) {
  marker <- file.path(root,"_meta","storage-layout.rds")
  if(!file.exists(marker))return(FALSE)
  schema <- readRDS(marker)
  if(!identical(schema$schema,"ocean-flat-v1"))stop("Unknown ocean storage schema: ",marker)
  TRUE
}
ocean_flat_relative <- function(path) {
  path <- gsub("\\\\","/",path)
  path <- sub("^_meta/initial(?=/|$)","full5_class4",path,perl=TRUE)
  path <- sub("^_meta/continuation(?=/|$)","continuation_plus20k",path,perl=TRUE)
  path <- sub("^_meta/cc(?=/|$)","complete_case_baselines/fits",path,perl=TRUE)
  path <- sub("^baseline_(fold[0-9]+)/(?=O[0-9]+(?:/|$))","full5_class4/\\1/",path,perl=TRUE)
  path <- sub("^(full5_class4|smoke_class4)/(fold[0-9]+)/(O[0-9]+)(?=/|$)","eiv_\\2_\\3",path,perl=TRUE)
  path <- sub("^(full5_class4|smoke_class4)/(fold[0-9]+)(?=/|$)","baseline_\\2",path,perl=TRUE)
  path <- sub("^(full5_class4|smoke_class4)(?=/|$)","_meta/initial",path,perl=TRUE)
  path <- sub("^continuation_plus20k/(fold[0-9]+)/(O[0-9]+)(?=/|$)","eiv_\\1_\\2_continued",path,perl=TRUE)
  path <- sub("^continuation_plus20k(?=/|$)","_meta/continuation",path,perl=TRUE)
  path <- sub("^complete_case_baselines/fits/(cc_[^/]+)(?=/|$)","\\1",path,perl=TRUE)
  path <- sub("^complete_case_baselines/fits(?=/|$)","_meta/cc",path,perl=TRUE)
  path <- sub("^workflow_tasks(?=/|$)","_meta/tasks",path,perl=TRUE)
  path <- sub("^jobs/([^/]+)/tasks/","logs/\\1/",path)
  path <- sub("^jobs(?=/|$)","logs",path,perl=TRUE)
  path <- sub("^logs/([^/]+)/tasks/","logs/\\1/",path)
  path <- sub("^(prediction_report|diagnostics|surface_report)(?=/|$)","reports/\\1",path,perl=TRUE)
  meta <- grepl("^(workflow_identity[.]rds|execution_settings[^/]*[.]rds|continuation_decisions[.]csv)$",path)
  path[meta] <- paste0("_meta/",path[meta])
  path
}
ocean_result_relative <- function(root,path) {
  if(ocean_storage_flat(root))ocean_flat_relative(path)else path
}
# Assemble first, then resolve against the nearest explicitly marked root.
# This also handles paths built incrementally by task-local closures.
ocean_path <- function(...,fsep=.Platform$file.sep) {
  paths <- file.path(...,fsep=fsep)
  vapply(paths,function(path) {
    candidate <- path
    repeat {
      if(dir.exists(candidate) && file.exists(file.path(candidate,"_meta","storage-layout.rds"))) {
        if(!ocean_storage_flat(candidate))return(path)
        if(identical(path,candidate))return(path)
        rel <- substring(path,nchar(candidate)+2L)
        return(file.path(candidate,ocean_flat_relative(rel)))
      }
      parent <- dirname(candidate)
      if(identical(parent,candidate)||candidate==".")return(path)
      candidate <- parent
    }
  },character(1),USE.NAMES=FALSE)
}
ocean_storage_init <- function(root) {
  if(ocean_storage_flat(root))return(invisible(root))
  if(dir.exists(file.path(root,".workflow.lock")))stop("Output directory is locked: ",root)
  existing <- if(dir.exists(root))list.files(root,all.files=TRUE,no..=TRUE)else character()
  # Never opt an existing, unmarked result tree into a new schema implicitly.
  if(length(existing))return(invisible(root))
  dir.create(file.path(root,"_meta"),recursive=TRUE,showWarnings=FALSE)
  saveRDS(list(schema="ocean-flat-v1"),file.path(root,"_meta","storage-layout.rds"))
  invisible(root)
}

# Dry-run by default. Copy into a new sibling, verify bytes, then publish it.
# The source remains usable; never rewrite model provenance to make it current.
ocean_migrate_storage <- function(source,destination,dry_run=TRUE) {
  source <- normalizePath(source,mustWork=TRUE)
  destination <- path.expand(destination)
  if(!dir.exists(dirname(destination)))stop("Destination parent must exist.")
  destination <- file.path(normalizePath(dirname(destination)),basename(destination))
  if(file.exists(destination)||dir.exists(destination))stop("Destination already exists: ",destination)
  if(startsWith(paste0(destination,"/"),paste0(source,"/")))stop("Destination must be outside source.")
  lock <- file.path(source,".workflow.lock")
  if(dir.exists(lock) || dir.exists(file.path(source,".application.lock")))
    stop("Source is locked. Wait for its job to finish; do not migrate an active run.")
  if(ocean_storage_flat(source))stop("Source already uses the flat layout.")
  files <- character(); pending <- ""
  while(length(pending)) {
    directory <- pending[1L];pending <- pending[-1L]
    entries <- list.files(file.path(source,directory),all.files=TRUE,no..=TRUE)
    relative <- if(nzchar(directory))file.path(directory,entries)else entries
    full <- file.path(source,relative)
    if(any(nzchar(Sys.readlink(full))))stop("Migration does not follow symbolic links.")
    dirs <- dir.exists(full)
    files <- c(files,relative[!dirs]);pending <- c(pending,relative[dirs])
  }
  files <- sort(files)
  mapped <- ocean_flat_relative(files)
  if(anyDuplicated(mapped))stop("Destination collision in migration plan.")
  plan <- data.frame(source=files,destination=mapped,bytes=file.info(file.path(source,files))$size)
  if(isTRUE(dry_run))return(plan)
  if(!dir.create(lock,showWarnings=FALSE))stop("Source became locked.")
  on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  current_files <- sort(list.files(source,recursive=TRUE,all.files=TRUE,no..=TRUE))
  current_files <- current_files[!dir.exists(file.path(source,current_files))]
  if(!identical(files,current_files))stop("Source inventory changed before migration acquired its lock.")
  stage <- tempfile(".ocean-migration-",tmpdir=dirname(destination))
  dir.create(stage);on.exit(if(dir.exists(stage))unlink(stage,recursive=TRUE),add=TRUE)
  before <- unname(tools::md5sum(file.path(source,files)))
  if(anyNA(before))stop("Cannot hash all source files.")
  for(i in seq_along(files)) {
    target <- file.path(stage,mapped[i]);dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
    if(!file.copy(file.path(source,files[i]),target,copy.mode=TRUE,copy.date=TRUE))stop("Copy failed: ",files[i])
  }
  if(!identical(before,unname(tools::md5sum(file.path(stage,mapped)))) ||
     !identical(before,unname(tools::md5sum(file.path(source,files)))))stop("Source changed or copy verification failed.")
  # Receipts record actual relative locations; scientific identities/hashes stay intact.
  receipts <- mapped[grepl("^_meta/tasks/.*[.]done[.]rds$",mapped)]
  for(p in receipts) {
    receipt <- readRDS(file.path(stage,p))
    if(any(grepl("(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$))",receipt$files)))stop("Unsafe receipt path: ",p)
    receipt$files <- ocean_flat_relative(receipt$files)
    if(!identical(unname(tools::md5sum(file.path(stage,receipt$files))),receipt$hashes))stop("Invalid source receipt: ",p)
    saveRDS(receipt,file.path(stage,p))
  }
  dir.create(file.path(stage,"_meta"),showWarnings=FALSE)
  saveRDS(list(schema="ocean-flat-v1"),file.path(stage,"_meta","storage-layout.rds"))
  plan$source_md5 <- before
  plan$destination_md5 <- unname(tools::md5sum(file.path(stage,mapped)))
  write.csv(plan,file.path(stage,"_meta","migration-manifest.csv"),row.names=FALSE)
  saveRDS(list(source=source,source_identity=readRDS(file.path(source,"workflow_identity.rds")),
    migrated_at=as.character(Sys.time()),note="Original scientific provenance retained; code compatibility is still enforced."),
    file.path(stage,"_meta","migration.rds"))
  if(file.exists(destination)||!file.rename(stage,destination))stop("Cannot publish migration destination.")
  invisible(plan)
}

# Task receipts may be imported between legacy and flat roots in either direction.
ocean_legacy_task_relative <- function(path,profile="paper") {
  initial <- if(isTRUE(profile=="quick"))"smoke_class4"else"full5_class4"
  path <- sub("^eiv_(fold[0-9]+)_(O[0-9]+)_continued(?=/|$)","continuation_plus20k/\\1/\\2",path,perl=TRUE)
  path <- sub("^eiv_(fold[0-9]+)_(O[0-9]+)(?=/|$)",paste0(initial,"/\\1/\\2"),path,perl=TRUE)
  path <- sub("^baseline_(fold[0-9]+)(?=/|$)",paste0(initial,"/\\1"),path,perl=TRUE)
  sub("^(cc_[^/]+)(?=/|$)","complete_case_baselines/fits/\\1",path,perl=TRUE)
}

# Keep only rows already generated; restore into freshly allocated arrays on
# restart. Handles legacy checkpoints containing the full allocated arrays too.
mixedgp_checkpoint_prefix <- function(draws, n) {
  lapply(draws,function(x) {
    if(is.null(x))return(NULL)
    dims<-dim(x)
    if(is.null(dims))return(x[seq_len(n)])
    do.call(`[`,c(list(x,seq_len(n)),rep(list(TRUE),length(dims)-1L),list(drop=FALSE)))
  })
}
mixedgp_checkpoint_restore_draws <- function(saved, allocated) {
  if(!identical(names(saved),names(allocated)))stop("Checkpoint draw fields differ")
  for(name in names(allocated)) {
    x<-saved[[name]];target<-allocated[[name]]
    if(is.null(x)&&is.null(target))next
    dims<-dim(target);n<-if(is.null(dim(x)))length(x)else dim(x)[1L]
    limit<-if(is.null(dims))length(target)else dims[1L]
    if(n>limit || !identical(dim(x)[-1L],dims[-1L]))stop("Checkpoint draw dimensions differ")
    if(n) {
      if(is.null(dims))target[seq_len(n)]<-x else
        target<-do.call(`[<-`,c(list(target,seq_len(n)),rep(list(TRUE),length(dims)-1L),list(value=x)))
    }
    allocated[name]<-list(target)
  }
  allocated
}
mixedgp_checkpoint_reference <- function(directory) {
  if(is.null(directory))return(NULL)
  list(directory=normalizePath(directory,mustWork=TRUE),
       identity_md5=unname(tools::md5sum(file.path(directory,"identity.rds"))))
}

# Retirement is allowed only when an exact read-back of a published full fit
# succeeds. Selected prediction states, scores and incomplete fits never qualify.
mixedgp_retire_saved_checkpoints <- function(object, path) {
  if(!is.null(attr(object,"eivgp_study1_restart_storage",exact=TRUE)))
    return(mixedgp_retire_study1_measurement(object,path))
  refs<-function(x) {
    if(!is.list(x)||is.data.frame(x))return(list())
    ref<-attr(x,"eivgp_restart_storage",exact=TRUE)
    if(!is.null(ref) && is.list(x$checkpoint) && length(x$mcmc$samples_by_chain)) {
      expected<-x$checkpoint$iteration-x$checkpoint$burn
      parts<-unlist(x$mcmc$samples_by_chain,recursive=FALSE)
      complete<-is.finite(expected)&&expected>0 && all(vapply(parts,function(a)
        is.null(a)||(if(is.null(dim(a)))length(a)else dim(a)[1L])==expected,logical(1)))
      if(complete)return(list(c(ref,list(checkpoint=x$checkpoint))))
      return(list())
    }
    if(is.list(x$mcmc))return(list())
    unlist(lapply(x,refs),recursive=FALSE)
  }
  candidates<-refs(object)
  root<-getOption("eivGP.checkpoint_root",NULL)
  if(!length(candidates)||is.null(root))return(invisible(NULL))
  # Do not risk successful fitting because cleanup encounters an I/O problem.
  tryCatch({
    if(!identical(readRDS(path),mixedgp_pack_artifact(object)))
      stop("Published artifact did not match read-back; restart files retained")
    root<-normalizePath(root,mustWork=TRUE)
    final<-normalizePath(path,mustWork=TRUE)
    for(ref in candidates) local({
      dir<-normalizePath(ref$directory,mustWork=FALSE)
      if(!startsWith(dir,paste0(root,.Platform$file.sep)) ||
         startsWith(final,paste0(dir,.Platform$file.sep)))return(NULL)
      manifest<-file.path(dir,"identity.rds")
      if(!file.exists(manifest)||!identical(unname(tools::md5sum(manifest)),ref$identity_md5))return(NULL)
      lock<-file.path(dir,"LOCK")
      if(!dir.create(lock,showWarnings=FALSE))return(NULL)
      on.exit(unlink(lock,recursive=TRUE),add=TRUE)
      verification<-mixedgp_verify_checkpoint_pair(dir,final)
      files<-verification$files$path
      verification$status<-"verified_before_deletion"
      mixedgp_checkpoint_write(verification,file.path(dir,"retired.rds"))
      if(!identical(unname(tools::md5sum(files)),verification$files$md5)||
         !identical(unname(tools::md5sum(final)),verification$fit_md5))
        stop("Files changed before retirement; restart files retained")
      if(length(files) && unlink(files)!=0L)stop("Some redundant restart files could not be removed")
      verification$status<-"removed"
      mixedgp_checkpoint_write(verification,file.path(dir,"retired.rds"))
      # Preserve the small identity and timing files for provenance and audit.

    })
  },error=function(e)warning("Checkpoint cleanup: ",conditionMessage(e),call.=FALSE))
  invisible(NULL)
}

# Exact comparison shared by automatic fit retirement and explicit legacy cleanup.
# Continuation-segment checkpoints are conservatively retained; their raw arrays
# cover only a suffix of the full fit and need a different history comparison.
mixedgp_verify_checkpoint_pair <- function(directory, fit_path) {
  need <- function(ok, why) if(!isTRUE(ok))stop(why,call.=FALSE)
  plain <- function(p) !isTRUE(nzchar(Sys.readlink(p)))
  need(plain(directory)&&plain(fit_path),"symlink_preserved")
  identity_path<-file.path(directory,"identity.rds")
  need(file.exists(identity_path)&&plain(identity_path),"missing_identity")
  identity_md5<-unname(tools::md5sum(identity_path))
  id<-readRDS(identity_path);m<-id$model
  need(is.list(m)&&is.null(m$component),"measurement_or_unknown_checkpoint_preserved")
  need(is.null(m$resume),"continuation_checkpoint_preserved")
  fit_md5<-unname(tools::md5sum(fit_path))
  fit<-readRDS(fit_path)
  for(j in seq_len(4L)) {
    if(inherits(fit,"mixedgp_compact_fit")) {
      need(identical(fit$version,1L),"unsupported_compact_fit")
      fit<-fit[["fit"]]
    } else if(is.list(fit)&&is.null(fit[["mcmc"]])&&is.list(fit[["fit"]])) {
      fit<-fit[["fit"]]
    } else break
  }
  need(is.list(fit),"not_a_full_fit")
  cp<-fit$checkpoint;by<-fit$mcmc$samples_by_chain
  need(is.list(cp)&&is.list(by)&&length(by)>0L,"not_a_full_fit")
  need(identical(cp$version,2L)&&identical(cp$sampler_version,"0.3.1"),"unsupported_fit_checkpoint")
  need(cp$iteration==m$n_iter&&cp$burn==m$burn&&cp$thin==1L,"budget_mismatch")
  nc<-as.integer(m$n_chains);n<-as.integer(m$n_iter-m$burn)
  need(length(nc)==1L&&nc>0L&&n>0L&&length(cp$states)==nc,"chain_count_mismatch")
  num_equal<-function(a,b)identical(dim(a),dim(b))&&identical(as.numeric(a),as.numeric(b))
  univariate<-identical(m$measurement,"threshold")
  d<-fit$data
  need(num_equal(m$X,if(univariate)as.matrix(d$x)else d$X)&&
    num_equal(m$y,if(univariate)as.numeric(d$y)else d$y)&&
    num_equal(m$C,if(univariate)matrix(d$c_ord,ncol=1L)else d$C_ord)&&
    num_equal(m$U_obs,if(univariate)matrix(d$u_obs,ncol=1L)else d$U_obs),"input_mismatch")
  need(identical(as.integer(m$m_vec),as.integer(if(univariate)d$m else d$m_vec))&&
    identical(m$kernel,normalize_gp_kernel(fit$kernel$name,fit$kernel$matern_nu))&&identical(m$priors,fit$priors)&&
    identical(m$control,fit$control)&&identical(m$control,cp$control)&&
    (univariate||identical(m$ident,d$ident)),"model_mismatch")
  need(identical(as.integer(m$seed),as.integer(cp$arguments$seed)),"seed_mismatch")
  need(all(vapply(by,function(parts)length(parts)==nc&&all(vapply(parts,function(a)
    is.null(a)||(if(is.null(dim(a)))length(a)else dim(a)[1L])==n,logical(1))),logical(1))),"incomplete_draw_history")
  paths<-file.path(directory,sprintf("chain-%02d.rds",seq_len(nc)))
  found<-list.files(directory,"^chain-[0-9]+[.]rds$",full.names=TRUE)
  need(setequal(paths,found)&&all(file.exists(paths))&&all(vapply(paths,plain,logical(1))),"incomplete_or_extra_chain_files")
  info<-file.info(paths);hashes<-unname(tools::md5sum(paths));elapsed<-numeric(nc)
  need(!anyNA(c(hashes,fit_md5,identity_md5)),"unreadable_file")
  for(i in seq_len(nc)) {
    raw<-readRDS(paths[i]);terminal<-cp$states[[i]]
    need(identical(raw$schema,1L)&&raw$iteration==m$n_iter&&raw$saved==n,"unfinished_chain")
    need(identical(raw$state,terminal$state)&&identical(raw$rng,terminal$rng)&&
      identical(raw$cache,terminal$cache)&&identical(raw$initial,fit$mcmc$chain_initial[[i]]),"restart_state_mismatch")
    for(field in names(raw$draws)) {
      name<-sub("^samples_","",field);if(univariate&&name=="U")name<-"u"
      a<-raw$draws[[field]];b<-if(name%in%names(by))by[[name]][[i]]else NULL
      if(univariate&&field=="samples_U"&&!is.null(a))a<-matrix(a,nrow=n)
      need(identical(dim(a),dim(b))&&identical(as.vector(a),as.vector(b)),paste0("draw_mismatch:",field,"/chain",i))
    }
    expected<-paste0("samples_",ifelse(names(by)=="u","U",names(by)))
    need(all(expected%in%names(raw$draws)),"missing_draw_fields")
    stats<-fit$mcmc$chain_stats
    need(nrow(stats)==nc&&stats$chain[i]==i&&stats$saved[i]==n,"chain_stats_mismatch")
    for(k in names(raw$counts))need(k%in%names(stats)&&
      identical(as.numeric(raw$counts[[k]]),as.numeric(stats[[k]][i])),paste0("counter_mismatch:",k))
    elapsed[i]<-raw$elapsed_seconds
  }
  # Recheck published fit and identities; timestamps alone are not integrity checks.
  need(identical(unname(tools::md5sum(fit_path)),fit_md5)&&
    identical(unname(tools::md5sum(identity_path)),identity_md5)&&
    identical(file.info(paths)[,c("size","mtime")],info[,c("size","mtime")]),"files_changed_during_verification")
  list(schema=1L,verified_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
    directory=normalizePath(directory),fit_path=normalizePath(fit_path),fit_md5=fit_md5,
    identity_md5=identity_md5,files=data.frame(path=paths,md5=hashes,bytes=info$size,
      elapsed_seconds=elapsed),checks="all inputs, model, budgets, draws, terminal state, RNG, cache, initial state and counters")
}

#' Verify and prune redundant conditional GP checkpoints
#'
#' Audit legacy Study II CUT restart files against published full conditional fits.
#' Filenames/seeds find candidates only; every chain and model input is compared.
#' Measurement checkpoints, continuation segments, missing/mismatched/locked fits,
#' final fits, identity files and timing logs are retained. No fitting is performed.
#' @param output_dir Study output directory containing checkpoints/ and mechanism/.
#' @param dry_run Logical; TRUE audits without modifying the output tree.
#' @param audit_dir Local directory outside output_dir for manifests and receipts.
#'   Required; use a fresh directory for each invocation.
#' @param progress_every Number of directories between progress messages.
#' @return A data frame with a status, reason and verified/removed bytes per directory.
#' @export
prune_study_checkpoints <- function(output_dir,dry_run=TRUE,audit_dir,progress_every=25L) {
  if(!is.logical(dry_run)||length(dry_run)!=1L||is.na(dry_run))stop("dry_run must be TRUE or FALSE")
  root<-normalizePath(path.expand(output_dir),mustWork=TRUE)
  checkpoint_root<-file.path(root,"checkpoints");mechanism<-file.path(root,"mechanism")
  if(!dir.exists(checkpoint_root)||!dir.exists(mechanism))stop("Expected checkpoints/ and mechanism/ under output_dir")
  audit<-normalizePath(path.expand(audit_dir),mustWork=FALSE)
  if(audit==root||startsWith(audit,paste0(root,.Platform$file.sep)))stop("Keep the audit directory outside the output tree")
  if(dir.exists(audit)&&length(list.files(audit,all.files=TRUE,no..=TRUE)))stop("Use a fresh audit directory")
  dir.create(audit,recursive=TRUE,showWarnings=FALSE);audit<-normalizePath(audit,mustWork=TRUE)
  dir.create(file.path(audit,"receipts"))
  progress_every<-mixedgp_as_integer_strict(progress_every,"progress_every",1L,1L)
  runlock<-file.path(mechanism,"RUN.lock")
  if(file.exists(runlock))stop("Run lock exists; no files changed: ",runlock)
  if(!dry_run) {
    if(!dir.create(runlock,showWarnings=FALSE))stop("Cannot acquire run lock")
    on.exit(unlink(runlock,recursive=TRUE),add=TRUE)
    writeLines(as.character(Sys.getpid()),file.path(runlock,"pid"))
  }
  protocol<-readRDS(file.path(mechanism,"protocol.rds"))
  base_seed<-protocol$settings$seed
  if(length(base_seed)!=1L||!is.finite(base_seed))stop("Missing protocol seed")
  paths<-list.files(mechanism,"^conditional-[0-9]+[.]rds$",recursive=TRUE,full.names=TRUE)
  jobs<-basename(dirname(paths));valid<-grepl("^rep[0-9]+-cal[0-9]+-CUT$",jobs)
  paths<-paths[valid];jobs<-jobs[valid]
  rep<-as.integer(sub("^rep([0-9]+)-.*$","\\1",jobs))
  cal<-as.integer(sub("^.*-cal([0-9]+)-CUT$","\\1",jobs))
  inner<-as.integer(sub("^conditional-([0-9]+)[.]rds$","\\1",basename(paths)))
  seeds<-base_seed+rep*1000L+cal+200000L+inner
  dirs<-sort(list.dirs(checkpoint_root,recursive=FALSE,full.names=TRUE))
  dirs<-dirs[grepl("^[0-9a-f]{32}$",basename(dirs))]
  results<-list()
  process<-function(directory) {
    row<-data.frame(directory=directory,status="preserved",reason="",fit_path="",
      verified_bytes=0,removed_bytes=0,chains_removed=0,stringsAsFactors=FALSE)
    manifest<-file.path(directory,"identity.rds");lock<-file.path(directory,"LOCK")
    if(nzchar(Sys.readlink(directory))||!file.exists(manifest)) {row$reason<-"symlink_or_missing_identity";return(row)}
    if(file.exists(lock)){row$reason<-"locked";return(row)}
    if(!dry_run) {
      if(!dir.create(lock,showWarnings=FALSE)){row$reason<-"locked";return(row)}
      on.exit(unlink(lock,recursive=TRUE),add=TRUE)
      writeLines(as.character(Sys.getpid()),file.path(lock,"pid"))
    }
    tryCatch({
      m<-readRDS(manifest)$model
      if(!is.null(m$component)){row$reason<-"measurement_checkpoint";return(row)}
      if(!is.null(m$resume)){row$reason<-"continuation_checkpoint";return(row)}
      if(!length(list.files(directory,"^chain-[0-9]+[.]rds$"))) {row$reason<-"no_chain_files";return(row)}
      candidates<-paths[which(seeds==m$seed)]
      if(!length(candidates)){row$reason<-"no_conditional_fit_candidate";return(row)}
      successes<-list();errors<-character()
      for(p in candidates) {
        canonical<-normalizePath(p,mustWork=TRUE)
        if(!startsWith(canonical,paste0(normalizePath(mechanism),.Platform$file.sep)))next
        v<-tryCatch(mixedgp_verify_checkpoint_pair(directory,p),error=function(e){errors<<-c(errors,conditionMessage(e));NULL})
        if(!is.null(v))successes[[length(successes)+1L]]<-v
      }
      if(length(successes)!=1L) {row$reason<-if(length(successes)>1L)"ambiguous_verified_match"else paste(unique(errors),collapse="; ");return(row)}
      v<-successes[[1L]];row$fit_path<-v$fit_path;row$verified_bytes<-sum(v$files$bytes)
      receipt<-file.path(audit,"receipts",paste0(basename(directory),".rds"))
      v$status<-if(dry_run)"verified_dry_run"else"verified_before_deletion"
      mixedgp_checkpoint_write(v,receipt)
      row$status<-"verified";row$reason<-"exact_match"
      if(!dry_run) {
        # Revalidate bytes immediately before unlinking, while holding both locks.
        if(!identical(unname(tools::md5sum(v$fit_path)),v$fit_md5)||
           !identical(unname(tools::md5sum(manifest)),v$identity_md5)||
           !identical(unname(tools::md5sum(v$files$path)),v$files$md5))stop("files_changed_before_deletion")
        result<-unlink(v$files$path) # Exact verified files only; never recursive.
        removed<-!file.exists(v$files$path)
        row$removed_bytes<-sum(v$files$bytes[removed]);row$chains_removed<-sum(removed)
        row$status<-if(result==0L&&all(removed))"removed"else"partial_removal"
        v$status<-row$status;v$removed<-removed
        mixedgp_checkpoint_write(v,receipt)
      }
      row
    },error=function(e){if(row$removed_bytes==0)row$status<-"preserved";row$reason<-conditionMessage(e);row})
  }
  for(i in seq_along(dirs)) {
    results[[i]]<-process(dirs[i])
    if(i%%progress_every==0L||i==length(dirs)) {
      table<-do.call(rbind,results)
      mixedgp_atomic_write_csv(table,file.path(audit,"manifest.csv"))
      message(i,"/",length(dirs)," directories; ",sum(table$status%in%c("verified","removed")),
        " verified fits; ",round(sum(table$removed_bytes)/1e9,3)," GB removed")
    }
  }
  if(!length(results))return(data.frame())
  do.call(rbind,results)
}

# Delayed checkpoint policy shared by GP and measurement samplers. With a delay,
# use script-wide clock boundaries so new fits do not reset the save schedule.
mixedgp_checkpoint_due <- function(now, last_save, iterations_since, force=FALSE) {
  seconds<-getOption("eivGP.checkpoint_seconds",Inf)
  iterations<-getOption("eivGP.checkpoint_iterations",Inf)
  if(!is.finite(seconds)&&!is.finite(iterations))return(FALSE)
  delay<-getOption("eivGP.checkpoint_after_seconds",0)
  start<-getOption("eivGP.checkpoint_script_start",0)
  if(now-start<delay)return(FALSE)
  time_due<-if(delay>0 && is.finite(seconds))
    floor((now-start)/seconds)>floor((last_save-start)/seconds) else now-last_save>=seconds
  isTRUE(force)||time_due||iterations_since>=iterations
}

# Retire Study I measurement restart copies only after a full artifact read-back.
mixedgp_retire_study1_measurement <- function(object,path) {
  tryCatch({
    ref<-attr(object,"eivgp_study1_restart_storage",exact=TRUE)
    root<-getOption("eivGP.checkpoint_root",NULL)
    if(is.null(root)||!dir.exists(ref$directory))return(invisible(NULL))
    directory<-normalizePath(ref$directory);final<-normalizePath(path,mustWork=TRUE)
    if(!startsWith(directory,paste0(normalizePath(root),.Platform$file.sep)) ||
       startsWith(final,paste0(directory,.Platform$file.sep)))return(invisible(NULL))
    lock<-file.path(directory,"LOCK")
    if(!dir.create(lock,showWarnings=FALSE))return(invisible(NULL))
    on.exit(unlink(lock,recursive=TRUE),add=TRUE)
    final_hash<-unname(tools::md5sum(final))
    if(!identical(readRDS(final),mixedgp_pack_artifact(object)))stop("Full measurement fit read-back mismatch")
    identity<-file.path(directory,"identity.rds")
    if(!identical(unname(tools::md5sum(identity)),ref$identity_md5))stop("Measurement checkpoint identity changed")
    model<-readRDS(identity)$model
    if(!identical(model$component,"study1-measurement"))stop("Unexpected measurement identity")
    files<-file.path(directory,sprintf("chain-%02d.rds",seq_along(object$chains)))
    found<-list.files(directory,"^chain-[0-9]+[.]rds$",full.names=TRUE)
    if(!setequal(files,found))return(invisible(NULL))
    hashes<-unname(tools::md5sum(files));sizes<-file.info(files)$size
    for(i in seq_along(files)) {
      raw<-readRDS(files[i]);fit<-object$chains[[i]];state<-fit$checkpoint
      same_draws<-function(a,b)identical(dim(a),dim(b))&&identical(as.vector(a),as.vector(b))
      if(!identical(raw$schema,1L)||raw$iteration!=model$iter||raw$saved!=model$iter-model$warmup||
         !identical(fit$data$c_ord,as.integer(model$C))||!identical(fit$data$u_obs,as.numeric(model$U_obs))||
         fit$control$seed!=model$seed+10000L*i||fit$control$n_iter!=model$iter||fit$control$burn!=model$warmup||
         !identical(raw$u,state$u)||!identical(raw$tau,state$tau)||!identical(raw$rng,state$rng)||
         !same_draws(raw$draws$samples_u,fit$samples_u)||!same_draws(raw$draws$samples_tau,fit$samples_tau))
        stop("Measurement checkpoint differs from full fit")
    }
    receipt<-list(schema=1L,artifact=final,artifact_md5=final_hash,
      files=data.frame(path=files,md5=hashes,bytes=sizes),status="verified_before_deletion",
      verified_utc=format(Sys.time(),tz="UTC",usetz=TRUE))
    mixedgp_checkpoint_write(receipt,file.path(directory,"retired.rds"))
    if(!identical(unname(tools::md5sum(final)),final_hash)||
       !identical(unname(tools::md5sum(identity)),ref$identity_md5)||
       !identical(unname(tools::md5sum(files)),hashes))stop("Measurement files changed before retirement")
    if(unlink(files)!=0L || any(file.exists(files)))stop("Some measurement checkpoints could not be removed")
    receipt$status<-"removed";mixedgp_checkpoint_write(receipt,file.path(directory,"retired.rds"))
  },error=function(e)warning("Measurement checkpoint cleanup: ",conditionMessage(e),call.=FALSE))
  invisible(NULL)
}


# ============================================================================
# ADNI prepared-fold application: adni-data.R
# ============================================================================

# Author: Xinye Chen. ADNI measurement, validation and scoring adapters.
# Models, prediction, optimizers and CPU scheduling remain in the unmodified installed upstream package.
adni_ezgp_mixed_matrix <- function(X, C, m_vec) {
    X <- as.matrix(X)
    C <- as.matrix(C)
    if (!is.numeric(X) || any(!is.finite(X)))
        stop("X must be a finite numeric matrix.")
    if (nrow(X) != nrow(C) || ncol(C) != length(m_vec))
        stop("X/C dimensions do not match.")
    if (any(!is.finite(C)) || any(C != round(C)))
        stop("C must contain finite integer codes.")
    for (j in seq_along(m_vec)) {
        if (any(C[, j] < 1L | C[, j] > m_vec[j]))
            stop("C contains an out-of-range category.")
    }
    out <- cbind(X, C)
    storage.mode(out) <- "double"
    colnames(out) <- c(paste0("x", seq_len(ncol(X))), paste0("c", seq_len(ncol(C))))
    out
}

adni_ezgp_cross_collisions <- function(train_x, test_x, digits = 5L) {
    train_key <- apply(round(as.matrix(train_x), digits), 1L, paste, collapse = "\r")
    test_key <- apply(round(as.matrix(test_x), digits), 1L, paste, collapse = "\r")
    sum(test_key %in% train_key)
}

adni_ezgp_predict_covariance_consistent <- function(X_new, model) {
    if (!inherits(model, "EzGP model"))
        stop("model must be an EzGP model.")
    if (!requireNamespace("EzGP", quietly = TRUE))
        stop("Package 'EzGP' is not installed.")
    X_new <- as.matrix(X_new)
    X_train <- as.matrix(model$data$X)
    Y <- as.numeric(model$data$Y)
    p <- as.integer(model$data$p)
    q <- as.integer(model$data$q)
    m <- as.integer(model$data$m)
    n <- nrow(X_train)
    tau <- as.numeric(model$data$tau)
    parv <- as.numeric(model$param)
    if (ncol(X_new) != ncol(X_train))
        stop("X_new has the wrong number of columns.")
    if (length(parv) == 0L || any(!is.finite(parv)) || !is.finite(tau) || tau < 0) {
        stop("EzGP fit has invalid parameters or nugget.")
    }
    K <- EzGP:::cov_m(X_train, p, q, m, n, parv, tau)
    Tm <- try(chol(K), silent = TRUE)
    if (inherits(Tm, "try-error"))
        stop("EzGP training covariance is not Cholesky-factorable.")
    invT <- backsolve(Tm, diag(n))
    invc <- invT %*% t(invT)
    one <- matrix(1, nrow = n, ncol = 1L)
    y <- matrix(Y, ncol = 1L)
    mu <- as.numeric((t(one) %*% invc %*% y)/sum(invc))
    npar <- 1L + q + p + p * sum(m)
    psum <- function(x1, x2, par2) sum(-par2 * (x1 - x2)^2)
    covx <- function(w1, w2) {
        par1 <- parv[seq_len(q + 1L)]
        par2 <- parv[(q + 2L):(q + 1L + p)]
        par3 <- parv[(q + 2L + p):npar]
        x1 <- w1[seq_len(p)]
        z1 <- w1[p + seq_len(q)]
        x2 <- w2[seq_len(p)]
        z2 <- w2[p + seq_len(q)]
        value <- par1[1L] * exp(psum(x1, x2, par2))
        for (j in seq_len(q)) {
            if (z1[j] == z2[j]) {
                level <- as.integer(z1[j])
                offset <- (sum(m[seq_len(j)]) - m[j] + level - 1L) * p
                idx <- offset + seq_len(p)
                value <- value + par1[j + 1L] * exp(psum(x1, x2, par3[idx]))
            }
        }
        value
    }
    means <- latent_var <- numeric(nrow(X_new))
    prior_var <- sum(parv[seq_len(q + 1L)])
    centered_y <- y - mu * one
    for (j in seq_len(nrow(X_new))) {
        gamma <- matrix(vapply(seq_len(n), function(i) covx(X_new[j, ], X_train[i, ]), numeric(1)), ncol = 1L)
        means[j] <- as.numeric(mu + t(gamma) %*% invc %*% centered_y)
        latent_var[j] <- as.numeric(prior_var - t(gamma) %*% invc %*% gamma)
    }
    predictive_var <- latent_var + tau
    tolerance <- 100 * (nrow(X_new) + 1L) * .Machine$double.eps * max(1, max(abs(predictive_var)))
    if (any(!is.finite(means)) || any(!is.finite(predictive_var)) || any(predictive_var < -tolerance)) {
        stop("Covariance-consistent EzGP prediction returned invalid predictive variances.")
    }
    predictive_var <- pmax(predictive_var, 0)
    list(Y_hat = means, latent_MSE = latent_var, predictive_MSE = predictive_var, tau = tau)
}

adni_ezgp_select_tau_cv <- function(train_x, y_train, p, q, m_vec, seed, tau_fractions, cv_folds = 3L, maxeval = 100L,
    progress = NULL, initial_audit = NULL) {
    n <- nrow(train_x)
    cv_folds <- max(2L, min(as.integer(cv_folds), n))
    set.seed(as.integer(seed))
    fold_id <- sample(rep(seq_len(cv_folds), length.out = n))
    response_var <- max(stats::var(y_train), .Machine$double.eps)
    tau_fractions <- unique(as.numeric(tau_fractions))
    tau_grid <- tau_fractions * response_var
    audit <- list()
    index <- 0L
    for (tt in seq_along(tau_grid)) {
        for (fold in seq_len(cv_folds)) {
            index <- index + 1L
            train_id <- fold_id != fold
            validation_id <- !train_id
            prior <- if (is.data.frame(initial_audit)) {
                initial_audit[initial_audit$tau_index == tt & initial_audit$inner_fold == fold, , drop = FALSE]
            }
            else data.frame()
            if (nrow(prior) == 1L) {
                audit[[index]] <- prior
                if (is.function(progress))
                  progress(do.call(rbind, audit), index, length(tau_grid) * cv_folds)
                next
            }
            started <- proc.time()[["elapsed"]]
            warnings <- character()
            ans <- tryCatch(withCallingHandlers({
                fit <- EzGP::EzGP_fit(X = train_x[train_id, , drop = FALSE], Y = as.numeric(y_train[train_id]),
                  p = p, q = q, m = m_vec, tau = tau_grid[tt], maxeval = as.integer(maxeval))
                if (length(fit$param) == 0L || any(!is.finite(fit$param))) {
                  stop("EzGP returned non-finite fitted parameters.")
                }
                pred <- adni_ezgp_predict_covariance_consistent(train_x[validation_id, , drop = FALSE], fit)
                error <- pred$Y_hat - y_train[validation_id]
                if (any(pred$predictive_MSE <= 0))
                  stop("EzGP CV returned zero predictive variance.")
                list(nlpd = mean(0.5 * log(2 * pi * pred$predictive_MSE) + 0.5 * error^2/pred$predictive_MSE),
                  mse = mean(error^2), message = "")
            }, warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            }), error = function(e) list(nlpd = NA_real_, mse = NA_real_, message = conditionMessage(e)))
            audit[[index]] <- data.frame(tau_index = tt, tau_fraction = tau_fractions[tt], tau = tau_grid[tt],
                inner_fold = fold, n_train = sum(train_id), n_validation = sum(validation_id), rounded_cross_collisions = adni_ezgp_cross_collisions(train_x[train_id,
                  , drop = FALSE], train_x[validation_id, , drop = FALSE]), nlpd = ans$nlpd, mse = ans$mse, elapsed_seconds = unname(proc.time()[["elapsed"]] -
                  started), message = ans$message, warnings = paste(unique(warnings), collapse = " | "), stringsAsFactors = FALSE)
            if (is.function(progress))
                progress(do.call(rbind, audit), index, length(tau_grid) * cv_folds)
        }
    }
    audit <- do.call(rbind, audit)
    candidate <- do.call(rbind, lapply(seq_along(tau_grid), function(tt) {
        z <- audit[audit$tau_index == tt, , drop = FALSE]
        data.frame(tau_index = tt, tau_fraction = tau_fractions[tt], tau = tau_grid[tt], cv_nlpd = if (all(is.finite(z$nlpd)))
            mean(z$nlpd)
        else Inf, cv_mse = if (all(is.finite(z$mse)))
            mean(z$mse)
        else Inf, successful_folds = sum(is.finite(z$nlpd)), stringsAsFactors = FALSE)
    }))
    if (!any(is.finite(candidate$cv_nlpd))) {
        stop("Covariance-consistent EzGP training-only CV failed for every nugget value.")
    }
    best <- which.min(candidate$cv_nlpd)
    list(tau = candidate$tau[best], selected_index = candidate$tau_index[best], cv_score = "nlpd", fold_id = fold_id,
        audit = audit, candidates = candidate)
}

adni_ezgp_draws <- function(mean, variance, n_draw, seed) {
    if (any(!is.finite(mean)) || any(!is.finite(variance)) || any(variance < 0)) {
        stop("Invalid predictive mean or variance.")
    }
    set.seed(as.integer(seed))
    matrix(stats::rnorm(n_draw * length(mean), rep(mean, each = n_draw), rep(sqrt(variance), each = n_draw)), nrow = n_draw,
        ncol = length(mean))
}

adni_plan_require <- function() {
    for (p in c("jsonlite", "digest", "kergp", "LVGP", "EzGP")) if (!requireNamespace(p, quietly = TRUE))
        stop("Required ADNI dependency: ", p)
}

adni_plan_input_root <- function(path) {
    path <- path.expand(path)
    if (dir.exists(file.path(path, "MC")))
        path <- file.path(path, "MC")
    normalizePath(path, mustWork = TRUE)
}

adni_plan_inputs <- function(data_dir) {
    if (!requireNamespace("jsonlite", quietly = TRUE))
        stop("Install jsonlite to read ADNI protocols.")
    root <- adni_plan_input_root(data_dir)
    all <- list()
    audit <- list()
    files <- character()
    held <- list()
    reference <- NULL
    for (f in 1:5) {
        by_cell <- list()
        for (cell in c("O0", "O20", "O50")) {
            dir <- file.path(root, sprintf("fold%02d", f), cell)
            paths <- file.path(dir, c("PROTOCOL.json", "private/train.csv", "private/test.csv"))
            if (!all(file.exists(paths)))
                stop("Missing private fold inputs: ", dir)
            p <- jsonlite::fromJSON(paths[1], simplifyVector = TRUE)
            tr <- read.csv(paths[2], stringsAsFactors = FALSE)
            te <- read.csv(paths[3], stringsAsFactors = FALSE)
            q <- as.integer(p$q)
            d <- as.integer(p$d)
            xc <- grep("^X[0-9]+$", names(tr), value = TRUE)
            xc <- xc[order(as.integer(sub("X", "", xc)))]
            cc <- paste0("C", seq_len(q))
            uc <- paste0("U", seq_len(d))
            needed <- c("RID", "Y", xc, cc, uc, "calibration")
            if (!length(xc) || !all(needed %in% names(tr)) || !all(needed %in% names(te)))
                stop("Invalid fold columns: ", dir)
            if(length(p$X_order)!=length(xc)||p$p!=length(xc))stop('Continuous-X labels/dimensions differ: ',dir)
            if(any(vapply(tr[,xc,drop=FALSE],function(v)length(unique(v))<=2L,TRUE)))
                stop('X must be continuous; fix binary/categorical population attributes before fitting: ',dir)
            if(any(te$calibration!=0L))stop('Test calibration must be zero: ',dir)
            if (anyNA(tr[, needed]) || anyNA(te[, needed]) || any(!is.finite(as.matrix(tr[, setdiff(needed, "RID")]))) ||
                any(!is.finite(as.matrix(te[, setdiff(needed, "RID")]))))
                stop("Non-finite private inputs: ", dir)
            if (anyDuplicated(tr$RID) || anyDuplicated(te$RID) || length(intersect(tr$RID, te$RID)))
                stop("Training/test ID overlap: ", dir)
            if (nrow(tr) != p$training_N || nrow(te) != p$test_N || p$fold != f || p$O != cell || nrow(tr) + nrow(te) !=
                p$N)
                stop("Protocol row counts do not match: ", dir)
            if (!all(tr$calibration %in% 0:1) || sum(tr$calibration) != p$observed_training_U)
                stop("Calibration mask differs: ", dir)
            if (!is.null(p$train_md5) && !identical(tolower(p$train_md5), unname(tools::md5sum(paths[2]))))
                stop("Training checksum differs: ", dir)
            if (!is.null(p$test_md5) && !identical(tolower(p$test_md5), unname(tools::md5sum(paths[3]))))
                stop("Test checksum differs: ", dir)
            if (length(p$U_order) != d || length(p$C_order) != q || length(p$m_vec) != q)
                stop("Protocol coordinate labels differ: ", dir)
            for (j in seq_len(q)) if (!all(tr[[cc[j]]] %in% seq_len(p$m_vec[j])) || !all(te[[cc[j]]] %in% seq_len(p$m_vec[j])) ||
                !all(seq_len(p$m_vec[j]) %in% tr[[cc[j]]]))
                stop("Invalid native ordinal levels: ", dir)
            observed <- tr$calibration == 1L
            available <- as.character(p$test_U_available_RID)
            if (cell == "O0") {
                if (any(observed) || length(available) || p$ident != "lower_triangular" || isTRUE(p$standardize_U))
                  stop("O0 must use unanchored model coordinates.")
            }
            else {
                if (!any(observed) || !isTRUE(p$standardize_U) || p$ident != "none")
                  stop("Calibrated U must use the stated raw-coordinate anchoring.")
                if (!all(available %in% as.character(te$RID)) || anyDuplicated(available))
                  stop("Invalid fixed test-U mask: ", dir)
                if (qr(cbind(1, as.matrix(tr[observed, uc, drop = FALSE])))$rank < d + 1L)
                  stop("Calibration vectors have insufficient affine rank.")
            }
            seed_names <- c("production", "posterior", "prediction", "mixed_prediction", "UC_GP", "LVGP", "EzGP",
                "EzGP_draws")
            if (!all(seed_names %in% names(p$seeds)))
                stop("Missing fixed protocol seeds: ", dir)
            for (n in seed_names) mixedgp_as_integer_strict(p$seeds[[n]], paste("seed", n), 0L, 1L)
            by_cell[[cell]] <- list(train = tr, test = te, protocol = p, x = xc, c = cc, u = uc, available = te$RID %in%
                available)
            files <- c(files, paths)
            audit[[length(audit) + 1L]] <- data.frame(fold = f, cell = cell, n_train = nrow(tr), n_test = nrow(te),
                n_calibration = sum(observed), n_test_U = sum(te$RID %in% available), q = q, d = d)
        }
        a <- by_cell$O0
        for (cell in c("O20", "O50")) {
            b <- by_cell[[cell]]
            if (!identical(a$train[, setdiff(names(a$train), "calibration")], b$train[, setdiff(names(b$train),
                "calibration")]) || !identical(a$test, b$test))
                stop("Cells must share the same frozen fold records.")
        }
        if (any(by_cell$O20$train$calibration > by_cell$O50$train$calibration) || any(by_cell$O20$available & !by_cell$O50$available))
            stop("Calibration/test-U masks must be nested.")
        ids <- sort(as.character(c(a$train$RID, a$test$RID)))
        if (is.null(reference))
            reference <- ids
        else if (!identical(reference, ids))
            stop("Fold populations differ.")
        held[[f]] <- as.character(a$test$RID)
        for (cell in names(by_cell)) all[[paste(f, cell, sep = "/")]] <- by_cell[[cell]]
    }
    if (anyDuplicated(unlist(held)) || !identical(sort(unlist(held)), reference))
        stop("Each person must be held out exactly once.")
    list(root = root, cells = all, audit = do.call(rbind, audit), files = files, N = length(reference))
}

adni_plan_fit_args <- function(z, s) {
    p <- z$protocol
    u <- as.matrix(z$train[, z$u, drop = FALSE])
    colnames(u) <- p$U_order
    u[z$train$calibration == 0L, ] <- NA_real_
    cmat <- as.matrix(z$train[, z$c, drop = FALSE])
    colnames(cmat) <- p$C_order
    # The input contract encodes each native item as1:m; JSON may simplify nested arrays.
    levels <- lapply(as.integer(p$m_vec), function(m) as.numeric(seq_len(m)))
    list(X = as.matrix(z$train[, z$x, drop = FALSE]), y = z$train$Y, C = cmat, U_obs = u, latent_dim = p$d, engine = "multivariate",
        m_vec = as.integer(p$m_vec), ordinal_levels = levels, ident = p$ident, standardize_U = isTRUE(p$standardize_U),
        kernel = "se", store_scores = TRUE, thin = 1L,
        sampler = if (!is.null(s$samplers)) s$samplers[[p$O]] else if (p$O == "O50") "ESS" else "HMC",
        parallel = s$cores > 1L, n_cores = s$cores,
        n_chains = s$chains, seed = p$seeds$production, verbose = TRUE)
}

adni_plan_fit_path <- function(out, task) file.path(out, "fits", sprintf("fold%02d", task$fold), task$cell, "fit.rds")

adni_plan_read_fit <- function(path) mixedgp_read_artifact(path)

adni_plan_save_fit <- function(fit, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(".adni-fit-", tmpdir = dirname(path))
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(mixedgp_pack_artifact(fit), tmp)
    restored <- mixedgp_read_artifact(tmp)
    if (!identical(restored, fit))
        stop("ADNI checkpoint read-back differs; previous checkpoint retained.")
    rm(restored)
    gc(FALSE)
    if (!file.rename(tmp, path))
        stop("Could not publish ADNI checkpoint; previous retained: ", path)
    invisible(path)
}

adni_plan_fit_task <- function(task, inputs, out, s, extend = FALSE) {
    z <- inputs$cells[[paste(task$fold, task$cell, sep = "/")]]
    path <- adni_plan_fit_path(out, task)
    old <- options(eivGP.checkpoint_root = NULL)
    on.exit(options(old), add = TRUE)
    fit <- if (file.exists(path))
        adni_plan_read_fit(path)
    else NULL
    target <- s$iter
    request <- file.path(dirname(path), "continuation_request.rds")
    if (extend) {
        if (is.null(fit))
            stop("Missing initial fit: ", path)
        tab <- adni_plan_diagnostics(fit)
        if (!file.exists(request)) {
            if (isTRUE(tab$summary$core_pass))
                return(invisible("SKIPPED_CONVERGED"))
            ocean_atomic_rds(list(from = fit$checkpoint$iteration, target = fit$checkpoint$iteration + s$continuation_iter,
                added = s$continuation_iter), request)
        }
        req <- readRDS(request)
        if (req$added != s$continuation_iter)
            stop("Saved continuation length differs; use its original --continuation-iter.")
        target <- req$target
    }
    else if (file.exists(request))
        target <- readRDS(request)$target
    if (is.null(fit)) {
        first <- min(target, s$first_checkpoint)
        fit <- do.call(fit_eivgp, c(adni_plan_fit_args(z, s), list(n_iter = first, burn = s$burn)))
        adni_plan_save_fit(fit, path)
    }
    while (fit$checkpoint$iteration < target) {
        current <- fit$checkpoint$iteration
        next_target <- min(target, (floor(current/s$checkpoint_step) + 1L) * s$checkpoint_step)
        fit <- continue_eivgp(fit, n_iter = next_target - current, parallel = s$cores > 1L, n_cores = s$cores,
            verbose = TRUE)
        adni_plan_save_fit(fit, path)
    }
    if (fit$checkpoint$burn != s$burn || length(fit$checkpoint$states) != s$chains)
        stop("Saved fit settings differ.")
    if(isTRUE(s$final_stage) && length(fit$checkpoint$states)>=4L) {
      X<-as.matrix(z$test[,z$x,drop=FALSE]);C<-as.matrix(z$test[,z$c,drop=FALSE])
      rows<-mixedgp_study2_diagnostic_rows(X,C,5L)
      panel<-list(X=X[rows,,drop=FALSE],C=C[rows,,drop=FALSE],
        U=if(task$cell!="O0")as.matrix(z$test[rows,z$u,drop=FALSE])else NULL,
        n_latent=64L,seed=z$protocol$seeds$prediction)
      fit<-mixedgp_application_ess(fit,s$ess_min,panel,file.path(dirname(path),"scientific-ess-state.rds"),s$cores)
      adni_plan_save_fit(fit,path)
      mixedgp_save_scientific_trace(fit,file.path(dirname(path),"scientific-trace-targets"),
        sprintf("ADNI | fold %d | %s",task$fold,task$cell))
      unlink(file.path(dirname(path),"scientific-ess-state.rds"))
    }
    dg <- adni_plan_diagnostics(fit)
    ocean_atomic_csv(dg$parameters, file.path(dirname(path), "diagnostic_parameters.csv"))
    ocean_atomic_csv(dg$summary, file.path(dirname(path), "diagnostic_summary.csv"))
    invisible("DONE")
}

adni_plan_diagnostics <- function(fit) {
    parts <- fit$diagnostics[grep("^rhat_", names(fit$diagnostics))]
    parts <- Filter(function(z) is.data.frame(z) && nrow(z)>0L && all(c("parameter", "rhat", "ess_bulk", "ess_tail") %in% names(z)),
        parts)
    d <- do.call(rbind, lapply(names(parts), function(n) {
        z <- parts[[n]][, c("parameter", "rhat", "ess_bulk", "ess_tail")]
        z$group <- n
        z
    }))
    core <- d[!grepl("^dictionary_occupancy", d$parameter), , drop = FALSE]
    pass <- function(z) nrow(z) > 0L && all(is.finite(z$rhat) & is.finite(z$ess_bulk) & is.finite(z$ess_tail) &
        z$rhat <= 1.01 & z$ess_bulk >= 400 & z$ess_tail >= 400)
    mx <- function(x) if (any(is.finite(x)))
        max(x[is.finite(x)])
    else NA_real_
    mn <- function(x) if (any(is.finite(x)))
        min(x[is.finite(x)])
    else NA_real_
    list(parameters = d, summary = data.frame(iterations = fit$checkpoint$iteration, burn = fit$checkpoint$burn,
        chains = length(fit$checkpoint$states), max_rhat = mx(core$rhat), min_bulk_ESS = mn(core$ess_bulk), min_tail_ESS = mn(core$ess_tail),
        core_nonfinite = sum(!is.finite(core$rhat) | !is.finite(core$ess_bulk) | !is.finite(core$ess_tail)), core_pass = pass(core),
        strict_all_pass = pass(d)))
}

adni_plan_score <- function(draws, point, truth, training_sd) {
    draws <- as.matrix(draws)
    point <- as.numeric(point)
    if (ncol(draws) != length(truth) || any(!is.finite(draws)) || any(!is.finite(point)) || !is.finite(training_sd) ||
        training_sd <= 0)
        stop("Invalid held-out predictive scores.")
    qq <- apply(draws, 2, quantile, c(0.025, 0.975))
    r <- sqrt(mean((point - truth)^2))
    a <- mean(abs(point - truth))
    crps <- mean(vapply(seq_along(truth), function(j) {
        v <- sort(draws[, j])
        n <- length(v)
        mean(abs(v - truth[j])) - sum((2 * seq_len(n) - n - 1) * v)/n^2
    }, 0))
    width <- mean(qq[2, ] - qq[1, ])
    interval <- mean(qq[2, ] - qq[1, ] + 40 * pmax(qq[1, ] - truth, 0) + 40 * pmax(truth - qq[2, ], 0))
    data.frame(n = length(truth), RMSE = r, RMSE_SD = r/training_sd, MAE = a, MAE_SD = a/training_sd, CRPS = crps,
        CRPS_SD = crps/training_sd, Coverage95 = mean(truth >= qq[1, ] & truth <= qq[2, ]), Width95 = width, Width95_SD = width/training_sd,
        IntervalScore95 = interval, IntervalScore95_SD = interval/training_sd)
}

adni_plan_baseline <- function(task, inputs, out, s) {
    dir <- file.path(out, "baselines", sprintf("fold%02d", task$fold), task$method)
    path <- file.path(dir, "fit.rds")
    if (all(file.exists(file.path(dir, c("fit.rds", "prediction_draws.rds")))))
        return(invisible("DONE"))
    z <- inputs$cells[[paste(task$fold, "O50", sep = "/")]]
    tr <- z$train
    te <- z$test
    p <- z$protocol
    mu <- colMeans(tr[, z$x, drop = FALSE])
    ss <- vapply(tr[, z$x, drop = FALSE], sd, 0)
    if (any(ss <= 0))
        stop("Constant ADNI training X.")
    sx <- function(d) sweep(sweep(as.matrix(d[, z$x, drop = FALSE]), 2, mu, "-"), 2, ss, "/")
    xc <- sx(tr)
    xs <- sx(te)
    yc <- mean(tr$Y)
    ys <- sd(tr$Y)
    yt <- (tr$Y - yc)/ys
    ct <- as.matrix(tr[, z$c, drop = FALSE])
    cs <- as.matrix(te[, z$c, drop = FALSE])
    mvec <- as.integer(p$m_vec)
    if (task$method != "EzGP") {
        controls <- if (s$profile == "paper") ocean_controls()[[task$method]] else
            if (task$method == "LVGP") list(n_starts=1L,max_retries=1L,max_iter_ini=10L,max_iter_lat=3L,
                rescue_iter_ini=10L,rescue_iter_lat=3L,max_elapsed_seconds=120,parallel=TRUE) else
                list(n_starts=1L,rescue_starts=integer())
        ans <- do.call(fit_mixedgp_competitor, c(list(method = task$method, X = xc, y = yt, C = ct, new_X = xs,
            new_C = cs, m_vec = mvec, n_draw = 4L * s$draws_per_chain, seed = p$seeds[[if (task$method == "UC-GP") "UC_GP" else "LVGP"]]),
            controls))
        point <- yc + ys * as.numeric(ans$predictive_mean)
        draws <- yc + ys * as.matrix(ans$draws)
        model <- ans
    }
    else {
        tx <- adni_ezgp_mixed_matrix(xc, ct, mvec)
        vx <- adni_ezgp_mixed_matrix(xs, cs, mvec)
        sel <- adni_ezgp_select_tau_cv(tx, yt, ncol(xc), length(z$c), mvec, p$seeds$EzGP, if (s$profile == "paper")
            c(1e-06, 0.0025, 0.01, 0.04, 0.16)
        else c(0.0025, 0.04), cv_folds = if (s$profile == "paper")
            3L
        else 2L, maxeval = if (s$profile == "paper")
            100L
        else 8L)
        model <- EzGP::EzGP_fit(tx, yt, p = ncol(xc), q = length(z$c), m = mvec, tau = sel$tau, maxeval = if (s$profile ==
            "paper")
            100L
        else 8L)
        pp <- adni_ezgp_predict_covariance_consistent(vx, model)
        point <- yc + ys * pp$Y_hat
        draws <- yc + ys * adni_ezgp_draws(pp$Y_hat, pp$predictive_MSE, 4L * s$draws_per_chain, p$seeds$EzGP_draws)
        ocean_atomic_csv(sel$audit, file.path(dir, "tau_cv_audit.csv"))
        ocean_atomic_csv(sel$candidates, file.path(dir, "tau_cv_candidates.csv"))
    }
    ocean_atomic_rds(list(model = model, center = mu, scale = ss, y_center = yc, y_scale = ys), path)
    mm <- cbind(data.frame(fold = task$fold, cell = "OI", method = task$method, protocol = "OI"), adni_plan_score(draws,
        point, te$Y, ys))
    ocean_atomic_csv(mm, file.path(dir, "metrics.csv"))
    ocean_atomic_rds(list(draws = draws, point = point), file.path(dir, "prediction_draws.rds"))
    invisible("DONE")
}

adni_plan_analyze_cell <- function(task, inputs, out, s) {
    path <- adni_plan_fit_path(out, task)
    if (!file.exists(path))
        stop("Missing saved ADNI fit: ", path)
    fit <- adni_plan_read_fit(path)
    z <- inputs$cells[[paste(task$fold, task$cell, sep = "/")]]
    p <- z$protocol
    tr <- z$train
    te <- z$test
    dir <- file.path(out, "analysis", sprintf("fold%02d", task$fold), task$cell)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    n <- nrow(fit$mcmc$samples_tau)
    cn <- n%/%s$chains
    set.seed(p$seeds$posterior)
    ids <- sort(unlist(lapply(0:(s$chains - 1L), function(i) i * cn + sample.int(cn, min(s$draws_per_chain, cn)))))
    ct <- as.matrix(te[, z$c, drop = FALSE])
    colnames(ct) <- p$C_order
    yd <- predict_eivgp(fit, new_X = as.matrix(te[, z$x, drop = FALSE]), new_C = ct, target = "response", draw_ids = ids,
        n_per_draw = 1L, joint = FALSE, seed = p$seeds$prediction)
    conditional <- attr(yd, "conditional_means", exact = TRUE)
    point <- if (is.null(conditional))
        colMeans(as.matrix(yd))
    else colMeans(conditional)
    score <- function(draws, point, protocol) cbind(data.frame(fold = task$fold, cell = task$cell, method = "EIV-GP",
        protocol = protocol), adni_plan_score(draws, point, te$Y, sd(tr$Y)))
    metrics <- score(yd, point, "OI")
    predictions <- data.frame(row = seq_len(nrow(te)), observed_Y = te$Y, predicted_Y = point)
    ocean_atomic_csv(predictions, file.path(dir, "Y_predictions.csv"))
    if (any(z$available)) {
        if (!isTRUE(fit$data$latent_scale_anchored))
            stop("Unanchored fit cannot receive raw test U.")
        sel <- which(z$available)
        u <- as.matrix(te[sel, z$u, drop = FALSE])
        colnames(u) <- p$U_order
        exact <- predict_eivgp(fit, new_X = as.matrix(te[sel, z$x, drop = FALSE]), new_C = ct[sel, , drop = FALSE],
            new_U = u, new_U_scale = "raw", target = "response", draw_ids = ids, n_per_draw = 1L, joint = FALSE,
            seed = p$seeds$mixed_prediction)
        mixed <- as.matrix(yd)
        mixed[, sel] <- as.matrix(exact)
        cm <- attr(exact, "conditional_means", exact = TRUE)
        mp <- point
        mp[sel] <- if (is.null(cm))
            colMeans(as.matrix(exact))
        else colMeans(cm)
        metrics <- rbind(metrics, score(mixed, mp, "MI"))
        ocean_atomic_csv(data.frame(row = seq_len(nrow(te)), observed_Y = te$Y, predicted_Y = mp, U_available = z$available),
            file.path(dir, "mixed_Y_predictions.csv"))
    }
    ocean_atomic_csv(metrics, file.path(dir, "prediction_metrics.csv"))
    dg <- adni_plan_diagnostics(fit)
    ocean_atomic_csv(cbind(fold = task$fold, cell = task$cell, dg$summary), file.path(dir, "convergence.csv"))
    ocean_atomic_csv(dg$parameters, file.path(dir, "diagnostic_parameters.csv"))
    scale <- if (any(tr$calibration == 1L))
        "raw"
    else "model"
    imp <- impute_eivgp(fit, rows = seq_len(nrow(tr)), draw_ids = ids, scale = scale)
    ocean_atomic_rds(list(draws = imp, scale = scale, draw_ids = ids, coordinates = p$U_order), file.path(dir,
        "U_draws.rds"))
    if (scale == "raw") {
        plots <- tables <- list()
        hidden <- which(tr$calibration == 0)
        cal <- which(tr$calibration == 1)
        common <- inputs$cells[[paste(task$fold, "O50", sep = "/")]]$train$calibration == 0
        for (j in seq_along(z$u)) {
            ud <- matrix(imp[, , j], nrow = length(ids))
            up <- colMeans(ud)
            qq <- apply(ud, 2, quantile, c(0.025, 0.975))
            truth <- tr[[z$u[j]]]
            sd0 <- sd(truth)
            tables[[j]] <- cbind(data.frame(fold = task$fold, cell = task$cell, coordinate = p$U_order[j]), adni_plan_score(ud[,
                hidden, drop = FALSE], up[hidden], truth[hidden], sd0), calibration_mean_RMSE = sqrt(mean((mean(truth[cal]) -
                truth[hidden])^2)))
            plots[[j]] <- data.frame(fold = task$fold, cell = task$cell, coordinate = p$U_order[j], row = seq_len(nrow(tr)),
                calibration = tr$calibration, common_hidden = common, true_U = truth, posterior_mean = up, q025 = qq[1,
                  ], q975 = qq[2, ])
        }
        ocean_atomic_csv(do.call(rbind, tables), file.path(dir, "U_metrics.csv"))
        ocean_atomic_csv(do.call(rbind, plots), file.path(dir, "U_predictions.csv"))
    }
    if (task$fold == s$folds[1]) {
        series <- mixedgp_v030_raw_series(fit)
        parameters <- head(grep("^(V$|r$|theta_x|tau|A\\[)", names(series), value = TRUE), 6L)
        if (length(parameters) && length(fit$checkpoint$states) >= 2L)
            ggplot2::ggsave(file.path(dir, "mcmc_trace.png"), plot_mcmc(fit, parameters), width = 10, height = 7,
                dpi = 130)
        else if (length(fit$checkpoint$states) < 2L)
            writeLines("Single-chain workflow check: multi-chain trace diagnostics are unavailable; convergence is not established.",
                file.path(dir, "MCMC_TRACE_UNAVAILABLE.txt"))
    }
    invisible("DONE")
}


# ============================================================================
# ADNI prepared-fold application: adni-evaluation.R
# ============================================================================

# Author: Xinye Chen. Preserve approved held-out and latent-recovery protocols as package code.
adni_save_plot <- function(filename, plot, width, height, dpi = 130, ...) {
    input <- paste0(filename, ".plot.rds")
    saveRDS(plot, input)
    code <- "lib<-Sys.getenv('EIVGP_R_LIBRARY','');if(nzchar(lib)).libPaths(c(lib,.libPaths()));eivGP::adni_workflow_api()$adni_render_plot(commandArgs(TRUE))"
    args <- c('-e', shQuote(code), shQuote(input), shQuote(filename), width, height, dpi)
    rc <- system2(file.path(R.home("bin"), "Rscript"), args, stdout = paste0(filename, ".render.log"), stderr = paste0(filename,
        ".render.log"))
    stopifnot(rc == 0L, file.exists(filename), file.info(filename)$size > 1000)
    unlink(input)
    invisible(filename)
}

adni_replace_plot_device <- function(x) {
    if (is.call(x) && identical(x, quote(ggplot2::ggsave)))
        return(as.name("adni_save_plot"))
    if (is.call(x))
        for (i in seq_along(x)) x[[i]] <- adni_replace_plot_device(x[[i]])
    x
}

adni_reference_points <- function(calibration, targets, c_names, u_name) {
    global <- rep(mean(calibration[[u_name]]), nrow(targets))
    c1 <- tapply(calibration[[u_name]], calibration$C1, mean)
    one <- as.numeric(c1[as.character(targets$C1)])
    one[!is.finite(one)] <- global[!is.finite(one)]
    key <- function(d) apply(d[, c_names, drop = FALSE], 1, paste, collapse = "/")
    joint <- tapply(calibration[[u_name]], key(calibration), mean)
    jp <- as.numeric(joint[key(targets)])
    fallback <- !is.finite(jp)
    jp[fallback] <- one[fallback]
    list(global = global, C1 = one, joint = jp, fallback = sum(fallback))
}

adni_analyze_cell <- function(task, inputs, out, s, stage) {
    z <- inputs$cells[[paste(task$fold, task$cell, sep = "/")]]
    tr <- z$train
    te <- z$test
    p <- z$protocol
    stageout <- file.path(out, stage)
    dir <- file.path(stageout, "analysis", sprintf("fold%02d", task$fold), task$cell)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    cache <- new.env(parent = emptyenv())
    env <- new.env(parent = environment(api$adni_plan_analyze_cell))
    env$adni_plan_fit_path <- function(ignored, task) api$adni_plan_fit_path(out, task)
    env$adni_plan_read_fit <- function(path) {
        cache$fit <- api$adni_plan_read_fit(path)
        stopifnot(cache$fit$checkpoint$iteration == s$iter)
        cache$fit
    }
    env$predict_eivgp <- function(...) {
        a <- list(...)
        r <- eivGP::predict_eivgp(...)
        if (is.null(a$new_U)) {
            stopifnot(is.null(a$new_y))
            cache$ids <- a$draw_ids
            cache$pure <- r
            atomic_rds(list(draws = r, draw_ids = a$draw_ids), file.path(dir, "Y_draws.rds"))
        }
        else {
            cache$exact <- r
            atomic_rds(list(draws = r, draw_ids = a$draw_ids), file.path(dir, "known_U_Y_draws.rds"))
        }
        r
    }
    env$adni_save_plot <- adni_save_plot
    fun <- api$adni_plan_analyze_cell
    body(fun) <- adni_replace_plot_device(body(fun))
    environment(fun) <- env
    fun(task, inputs, stageout, s)
    fit <- cache$fit
    stopifnot(fit$data$p == length(z$x), fit$data$q == p$q, fit$data$d == length(z$u), all(is.na(fit$data$U_obs_raw[tr$calibration ==
        0, , drop = FALSE])))
    met <- read.csv(file.path(dir, "prediction_metrics.csv"))
    met$stage <- stage
    met$iterations <- s$iter
    met$sampler <- fit$control$sampler
    if (exists("exact", cache, inherits = FALSE)) {
        ids <- which(z$available)
        cm <- attr(cache$exact, "conditional_means", exact = TRUE)
        stopifnot(!is.null(cm))
        extra <- cbind(data.frame(fold = task$fold, cell = task$cell, method = "EIV-GP", protocol = "CC-known"),
            api$adni_plan_score(cache$exact, colMeans(cm), te$Y[ids], sd(tr$Y)))
        extra$stage <- stage
        extra$iterations <- s$iter
        extra$sampler <- fit$control$sampler
        met <- rbind(met, extra)
    }
    atomic_csv(met, file.path(dir, "prediction_metrics.csv"))
    dg <- read.csv(file.path(dir, "convergence.csv"))
    dg$stage <- stage
    dg$sampler <- fit$control$sampler
    if (!is.null(fit$diagnostics$hmc))
        atomic_csv(fit$diagnostics$hmc, file.path(dir, "HMC_DIAGNOSTICS.csv"))
    atomic_csv(dg, file.path(dir, "convergence.csv"))
    cm <- attr(cache$pure, "conditional_means", exact = TRUE)
    per <- nrow(cm)%/%s$chains
    target_diagnostics <- do.call(rbind, lapply(seq_len(ncol(cm)), function(j) {
        a <- matrix(cm[, j], per, s$chains)
        data.frame(RID = te$RID[j], predicted_Y = mean(a), selected_Rhat = posterior::rhat(a), selected_ESS = posterior::ess_bulk(a),
            selected_MCSE = posterior::mcse_mean(a))
    }))
    atomic_csv(target_diagnostics, file.path(dir, "PREDICTION_FUNCTION_DIAGNOSTICS.csv"))
    if (task$cell != "O0") {
        cal <- tr[tr$calibration == 1, , drop = FALSE]
        hidden <- which(tr$calibration == 0)
        ui <- read.csv(file.path(dir, "U_predictions.csv"))
        um <- read.csv(file.path(dir, "U_metrics.csv"))
        for (j in seq_along(z$u)) {
            truth <- tr[[z$u[j]]]
            ref <- adni_reference_points(cal, tr[hidden, , drop = FALSE], z$c, z$u[j])
            row <- which(um$coordinate == p$U_order[j])
            um$C1_mean_RMSE[row] <- sqrt(mean((ref$C1 - truth[hidden])^2))
            um$joint_C_mean_RMSE[row] <- sqrt(mean((ref$joint - truth[hidden])^2))
            um$joint_C_fallback_N[row] <- ref$fallback
            ui$C1[ui$coordinate == p$U_order[j]] <- tr$C1
            hit <- ui$coordinate == p$U_order[j]
            ui$true_std[hit] <- (ui$true_U[hit] - mean(truth))/sd(truth)
            ui$mean_std[hit] <- (ui$posterior_mean[hit] - mean(truth))/sd(truth)
            ui$lo_std[hit] <- (ui$q025[hit] - mean(truth))/sd(truth)
            ui$hi_std[hit] <- (ui$q975[hit] - mean(truth))/sd(truth)
        }
        um$stage <- stage
        atomic_csv(um, file.path(dir, "U_metrics.csv"))
        atomic_csv(ui, file.path(dir, "U_predictions.csv"))
        ct <- as.matrix(te[, z$c, drop = FALSE])
        colnames(ct) <- p$C_order
        prospective <- impute_eivgp(fit, new_C = ct, draw_ids = cache$ids, scale = "raw", seed = p$seeds$prospective_U)
        stopifnot(identical(dim(prospective), c(length(cache$ids), nrow(te), length(z$u))))
        pm <- list()
        pp <- list()
        for (j in seq_along(z$u)) {
            ud <- matrix(prospective[, , j], nrow = length(cache$ids))
            point <- colMeans(ud)
            truth <- te[[z$u[j]]]
            ref <- adni_reference_points(cal, te, z$c, z$u[j])
            qq <- apply(ud, 2, quantile, c(0.025, 0.975))
            pm[[j]] <- cbind(data.frame(stage = stage, fold = task$fold, cell = task$cell, coordinate = p$U_order[j]),
                api$adni_plan_score(ud, point, truth, sd(tr[[z$u[j]]])), calibration_mean_RMSE = sqrt(mean((ref$global -
                  truth)^2)), C1_mean_RMSE = sqrt(mean((ref$C1 - truth)^2)), joint_C_mean_RMSE = sqrt(mean((ref$joint -
                  truth)^2)), joint_C_fallback_N = ref$fallback)
            pp[[j]] <- data.frame(RID = te$RID, coordinate = p$U_order[j], true_U = truth, posterior_mean = point,
                q025 = qq[1, ], q975 = qq[2, ])
        }
        atomic_csv(do.call(rbind, pm), file.path(dir, "PROSPECTIVE_U_metrics.csv"))
        atomic_csv(do.call(rbind, pp), file.path(dir, "PROSPECTIVE_U_predictions.csv"))
        ui$role <- factor(ifelse(ui$calibration == 1, "Observed calibration", "Masked training U"), levels = c("Masked training U",
            "Observed calibration"))
        plot <- ggplot(ui, aes(true_std, mean_std, colour = factor(C1))) + geom_abline(slope = 1, intercept = 0,
            linetype = 2, colour = "grey50") + geom_linerange(data = ui[ui$calibration == 0, ], aes(ymin = lo_std,
            ymax = hi_std), alpha = 0.15, linewidth = 0.3) + geom_point(aes(shape = role), size = 1.2, alpha = 0.65) +
            scale_shape_manual(values = c(16, 1)) + facet_wrap(~coordinate) + theme_bw(base_size = 11) + theme(legend.position = "bottom") +
            labs(x = "Reference U (training SD units; evaluation only)", y = "Posterior mean (training SD units)",
                colour = "MEMORY3", shape = NULL, title = paste(p$plan_id, "fold", task$fold, task$cell, stage),
                subtitle = paste(fit$control$sampler, "; CORE pass:", dg$core_pass))
        adni_save_plot(file.path(dir, "training_U_imputation_std.png"), plot, width = 10, height = 5.3, dpi = 140,
            bg = "white")
    }
    else writeLines("O0 is unanchored: no physical-U RMSE or truth-based alignment.", file.path(dir, "O0_SCALE_LIMITATION.txt"))
    atomic_json(list(state = "PASS", stage = stage, iterations = s$iter, primary_no_test_U_Y = TRUE, hidden_training_U_masked = TRUE,
        calibration_excluded_from_U_scoring = TRUE, prospective_U_uses_C_not_test_Y = TRUE, sampler = fit$control$sampler),
        file.path(dir, "EVALUATION_VALIDATION.json"))
    rm(fit, cache)
    gc(FALSE)
    invisible("DONE")
}

adni_baseline_cell <- function(task, inputs, out, s, cc = FALSE) {
    source <- inputs$cells[[paste(task$fold, if (cc)
        task$cell
    else "O50", sep = "/")]]
    tr <- source$train
    te <- source$test
    root <- if (cc)
        file.path(out, "CC", task$cell)
    else out
    inp <- inputs
    if (cc) {
        z <- source
        z$train <- tr[tr$calibration == 1, , drop = FALSE]
        z$test <- te
        for (j in seq_along(z$u)) {
            name <- paste0("X", length(z$x) + j)
            z$train[[name]] <- z$train[[z$u[j]]]
            z$test[[name]] <- z$test[[z$u[j]]]
        }
        z$x <- paste0("X", seq_len(length(z$x) + length(z$u)))
        for (n in c("UC_GP", "LVGP", "EzGP", "EzGP_draws")) z$protocol$seeds[[n]] <- as.integer((as.double(z$protocol$seeds[[n]]) +
            if (task$cell == "O20") 17041 else 35023)%%(2^30 - 1) + 1)
        inp$cells[[paste(task$fold, "O50", sep = "/")]] <- z
    }
    api$adni_plan_baseline(task, inp, root, s)
    dir <- file.path(root, "baselines", sprintf("fold%02d", task$fold), task$method)
    saved <- readRDS(file.path(dir, "prediction_draws.rds"))
    selected <- if (cc)
        which(source$available)
    else seq_len(nrow(te))
    mm <- cbind(data.frame(stage = "BASELINE", fold = task$fold, cell = if (cc)
        task$cell
    else "OI", method = paste0(task$method, if (cc)
        "-CC"
    else ""), protocol = if (cc)
        "CC-known"
    else "OI"), api$adni_plan_score(saved$draws[, selected, drop = FALSE], saved$point[selected], te$Y[selected],
        sd(tr$Y)))
    mm$training_N <- if (cc)
        sum(tr$calibration == 1)
    else nrow(tr)
    atomic_csv(mm, file.path(dir, "metrics.csv"))
    atomic_csv(data.frame(RID = te$RID, observed_Y = te$Y, predicted_Y = saved$point, scored = seq_len(nrow(te)) %in%
        selected), file.path(dir, "Y_predictions.csv"))
    invisible("DONE")
}


# ============================================================================
# ADNI prepared-fold application: adni-reporting.R
# ============================================================================

# Author: Xinye Chen. Saved-only ADNI report data, metrics and figures.
adni_table_html <- function(d,caption='',note='') {
  if(!is.data.frame(d)||!nrow(d))return('<p>No saved result available for this quantity.</p>')
  for(j in seq_along(d))if(is.numeric(d[[j]]))d[[j]]<-round(d[[j]],4)
  # Reuse Ocean's padding, rules, header shading, escaping and bold-mean markup.
  eivGP::ocean_workflow_api()$ocean_table_html(d,caption,note)
}
adni_bind <- function(parts) {
  parts<-Filter(function(x)is.data.frame(x)&&nrow(x)>0L,parts)
  if(!length(parts))return(data.frame())
  cols<-unique(unlist(lapply(parts,names)))
  do.call(rbind,lapply(parts,function(x){for(n in setdiff(cols,names(x)))x[[n]]<-NA;x[,cols,drop=FALSE]}))
}
adni_csv <- function(path,required=TRUE) {
  if(!file.exists(path)){if(required)stop('Missing saved report prerequisite: ',path);return(data.frame())}
  read.csv(path,stringsAsFactors=FALSE,check.names=FALSE)
}
adni_summary <- function(d,groups,metrics) {
  if(!nrow(d))return(data.frame())
  split_key<-interaction(d[,groups,drop=FALSE],drop=TRUE,lex.order=TRUE)
  adni_bind(lapply(split(d,split_key),function(x){
    out<-x[1,groups,drop=FALSE];out$folds<-length(unique(x$fold))
    for(k in metrics){v<-as.numeric(x[[k]]);out[[paste0(k,'_mean')]]<-mean(v);out[[paste0(k,'_SD')]]<-if(length(v)>1)sd(v)else NA_real_}
    out
  }))
}
adni_C_dictionary <- function(items,levels,labels=NULL) {
  meanings<-c(MEMORY1='Remembering shopping items without a list',
    MEMORY2='Remembering recent events',MEMORY3='Recalling conversations a few days later',
    MEMORY4='Remembering where objects were placed',MEMORY5='Repeating stories or questions',
    MEMORY6='Remembering the current date or day',MEMORY7='Remembering already having told someone something',
    MEMORY8='Remembering appointments, meetings or engagements',
    DIVATT3='Concentrating without distraction from the environment')
  if(length(labels))for(n in names(labels))meanings[n]<-as.character(labels[[n]])
  text<-unname(meanings[items]);text[is.na(text)]<-'Description not supplied in dataset metadata'
  data.frame(Model_column=paste0('C',seq_along(items)),ECog_field=items,
    Item_meaning=text,Native_grades=vapply(levels,function(n)paste(seq_len(n),collapse='/'),''),check.names=FALSE)
}
adni_saved_source <- function(output_dir,source=NULL) {
  if(!is.null(source)&&nzchar(source)) {
    m<-jsonlite::fromJSON(source,simplifyVector=FALSE)
    if(!identical(m$schema,'adni-saved-report-v1'))stop('Unsupported saved-report source')
    return(m)
  }
  p<-file.path(output_dir,'application_manifest.json')
  if(!file.exists(p))stop('Saved report requires application_manifest.json or --report-source=FILE')
  m<-jsonlite::fromJSON(p,simplifyVector=FALSE)
  list(schema='adni-saved-report-v1',output=output_dir,input=m$cases[[1]]$input,
    stage=paste0('S',tail(unlist(m$stages),1)),fit_package_version=m$package_version,fit_profile=m$profile,
    fitting_manifest=p,input_sha256=m$frozen_inputs)
}
adni_report_prepare <- function(repo,data_dir,output_dir,settings,report_source=NULL,report_dir=output_dir,group='report') {
  for(p in c('jsonlite','digest','ggplot2','patchwork','posterior'))if(!requireNamespace(p,quietly=TRUE))stop('Report dependency: ',p)
  api<-adni_load_api(repo);api$application$adni_output_check(report_dir,repo)
  m<-adni_saved_source(output_dir,report_source)
  fitting<-if(!is.null(m$fitting_manifest))jsonlite::fromJSON(m$fitting_manifest,simplifyVector=FALSE)else list()
  for(n in setdiff(c('folds','cells'),attr(settings,'explicit_selection')))if(length(fitting[[n]]))
    settings[[n]]<-if(n=='folds')as.integer(unlist(fitting[[n]]))else unlist(fitting[[n]],use.names=FALSE)
  for(n in c('draws_per_chain','baseline_profile'))if(!is.null(fitting[[n]]))settings[[n]]<-fitting[[n]]
  if(length(fitting$samplers))for(cell in names(fitting$samplers))
    settings[[paste0('sampler_',tolower(cell))]]<-fitting$samplers[[cell]]
  if(!is.null(m$fit_profile))settings$profile<-m$fit_profile
  if(!is.null(m$input))data_dir<-m$input
  if(!is.null(m$output))output_dir<-m$output
  input<-api$adni_plan_inputs(data_dir)
  if(length(m$input_sha256))for(p in names(m$input_sha256))
    if(digest::digest(file=p,algo='sha256')!=m$input_sha256[[p]])stop('Frozen report input changed: ',p)
  # Historical fitting identities are read, not rewritten to match reporting code.
  identity<-file.path(output_dir,'workflow_identity.rds')
  if(file.exists(identity)) {
    id<-readRDS(identity)
    if(!is.null(id$inputs)&&!identical(unname(tools::md5sum(input$files)),id$inputs))stop('Saved fit input identity differs')
  }
  if(is.null(m$stage)) {
    stages<-list.dirs(output_dir,recursive=FALSE,full.names=FALSE)
    stages<-stages[grepl('^S[0-9]+K?$',stages)]
    amount<-as.numeric(gsub('K','',sub('^S','',stages)))*ifelse(grepl('K$',stages),1000,1)
    if(!length(stages))stop('No saved EIV stages')
    m$stage<-stages[which.max(amount)]
  }
  pred<-conv<-up<-um<-pp<-pm<-parameters<-statuses<-paths<-list()
  final<-m$stage
  selected<-expand.grid(fold=settings$folds,cell=settings$cells,stringsAsFactors=FALSE)
  want_eiv<-group!='competitors';want_baseline<-group!='eivGP'
  for(i in seq_len(nrow(selected))) {
    f<-selected$fold[i];cell<-selected$cell[i];z<-input$cells[[paste(f,cell,sep='/')]]
    dir<-file.path(output_dir,final,'analysis',sprintf('fold%02d',f),cell)
    if(want_eiv) {
      d<-adni_csv(file.path(dir,'prediction_metrics.csv'));d<-d[d$protocol%in%c('OI','CC-known'),,drop=FALSE]
      pred[[length(pred)+1L]]<-d
      cv<-adni_csv(file.path(dir,'convergence.csv'));conv[[length(conv)+1L]]<-cv
      par<-adni_csv(file.path(dir,'diagnostic_parameters.csv'));par$fold<-f;par$cell<-cell
      parameters[[length(parameters)+1L]]<-par
      paths[[paste(f,cell,sep='/')]]<-list(fold=f,cell=cell,dir=dir,
        fit=file.path(output_dir,'fits',sprintf('fold%02d',f),cell,'fit.rds'))
      if(!file.exists(paths[[paste(f,cell,sep='/')]]$fit))stop('Missing saved EIV fit; report never fits')
      if(cell!='O0') {
        u<-adni_csv(file.path(dir,'U_predictions.csv'));u$fold<-f;u$cell<-cell
        up[[length(up)+1L]]<-u;um[[length(um)+1L]]<-adni_csv(file.path(dir,'U_metrics.csv'))
        p<-adni_csv(file.path(dir,'PROSPECTIVE_U_predictions.csv'));p$fold<-f;p$cell<-cell
        pp[[length(pp)+1L]]<-p;pm[[length(pm)+1L]]<-adni_csv(file.path(dir,'PROSPECTIVE_U_metrics.csv'))
      }
    }
  }
  if(want_baseline)for(f in settings$folds)for(cell in c('OI',intersect(settings$cells,c('O20','O50'))))for(method in c('UC-GP','LVGP','EzGP')) {
    dir<-if(cell=='OI')file.path(output_dir,'baselines',sprintf('fold%02d',f),method)else
      file.path(output_dir,'CC',cell,'baselines',sprintf('fold%02d',f),method)
    version<-m$fit_package_version
    if(length(m$baseline_overrides))for(b in m$baseline_overrides)
      if(b$fold==f&&b$cell==cell&&b$method==method){dir<-b$directory;version<-b$package_version}
    valid<-all(file.exists(file.path(dir,c('metrics.csv','fit.rds','prediction_draws.rds'))))
    state<-'MISSING';optimization<-'not_available'
    if(valid) {
      d<-adni_csv(file.path(dir,'metrics.csv'));d$method<-paste0(method,if(cell=='OI')''else'-CC');d$cell<-cell
      d$stage<-'BASELINE';pred[[length(pred)+1L]]<-d
      a<-readRDS(file.path(dir,'fit.rds'));optimization<-a$model$optimization_status
      if(is.null(optimization))optimization<-'NOT_EXPOSED'
      state<-if(grepl('nonconverged',optimization))'OPTIMIZER_FLAGGED'else'AVAILABLE'
      rm(a)
    }
    statuses[[length(statuses)+1L]]<-data.frame(fold=f,cell=cell,method=paste0(method,if(cell=='OI')''else'-CC'),
      state=state,optimization_status=optimization,fit_package_version=if(is.null(version))'recorded in saved provenance'else version,
      source_directory=dir)
  }
  p<-input$cells[[1]]$protocol
  z<-input$cells[['1/O50']];cohort<-rbind(z$train,z$test);cohort<-cohort[order(cohort$RID),]
  metadata_path<-file.path(data_dir,'DATASET.json')
  meta<-if(file.exists(metadata_path))jsonlite::fromJSON(metadata_path,simplifyVector=FALSE)else list()
  meta<-modifyList(meta,if(is.null(m$metadata))list()else m$metadata)
  if(is.null(meta))meta<-list()
  if(is.null(meta$title))meta$title<-'ADNI memory: partially observed brain imaging'
  if(is.null(meta$population))meta$population<-if(is.null(p$population))'Population specified in private preparation records'else p$population
  if(is.null(meta$calibration_design))meta$calibration_design<-if(is.null(p$calibration_design))
    'C1-proportional stratified sampling, fixed seeded within-stratum order; O20 nested in O50'else
    paste(p$calibration_design,'; strata:',paste(unlist(p$stratification_columns),collapse=', '))
  analysis<-list(repo=repo,api=api,source=m,input=input,output_dir=output_dir,report_dir=report_dir,
    settings=settings,group=group,stage=final,metadata=meta,protocol=p,cohort=cohort,
    prediction=adni_bind(pred),convergence=adni_bind(conv),parameters=adni_bind(parameters),
    U_predictions=adni_bind(up),U_metrics=adni_bind(um),prospective_predictions=adni_bind(pp),
    prospective_metrics=adni_bind(pm),method_status=adni_bind(statuses),paths=paths)
  if(nrow(analysis$convergence)) {
    nums<-unique(analysis$convergence[,c('iterations','burn','chains')])
    nums<-mixedgp_application_numerics(nums)
    analysis$numerics<-nums
    analysis$settings$chains<-as.integer(nums$chains)
    analysis$settings$iter<-as.integer(nums$iterations)
    analysis$settings$burn<-as.integer(nums$burn)
  }else analysis$numerics<-data.frame(iterations=NA,burn=NA,chains=NA)
  if(!is.null(m$fitting_manifest)) {
    stages<-paste0('S',unlist(fitting$stages))
    analysis$settings$stages<-as.integer(unlist(fitting$stages))
  }else {
    stages<-paste0('S',unlist(m$evaluation_iterations))
    if(length(m$evaluation_iterations))analysis$settings$stages<-as.integer(unlist(m$evaluation_iterations))
  }
  # Historical folders may use explicit stage labels rather than S<iteration>.
  if(length(m$stage_directories))stages<-unlist(m$stage_directories)
  stability<-list()
  for(st in stages)for(i in seq_len(nrow(selected))) {
    path<-file.path(output_dir,st,'analysis',sprintf('fold%02d',selected$fold[i]),selected$cell[i],'prediction_metrics.csv')
    d<-adni_csv(path,FALSE)
    if(nrow(d)){d<-d[d$protocol=='OI',];d$stage<-st;stability[[length(stability)+1L]]<-d}
  }
  analysis$stage_scores<-adni_bind(stability)
  adni_report_exports(analysis)
  analysis$plots<-adni_report_plots(analysis)
  if(want_eiv) {
    analysis$functionals<-adni_saved_functionals(analysis)
  }
  analysis$artifact_manifest<-adni_paper_assets(analysis)
  analysis
}
adni_report_exports <- function(a) {
  dest<-file.path(a$report_dir,'report-assets');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  for(n in c('prediction','convergence','parameters','U_metrics','prospective_metrics','method_status','stage_scores'))
    if(nrow(a[[n]]))write.csv(a[[n]],file.path(dest,paste0(n,'.csv')),row.names=FALSE)
  a$api$atomic_json(list(author='Xinye Chen',saved_only=TRUE,new_MCMC=FALSE,new_optimization=FALSE,
    fitting_source=a$source,stage=a$stage,saved_fitting_numerics=a$numerics,report_settings=a$settings,
    input_md5=as.list(tools::md5sum(a$input$files)),analysis_modules=a$api$eivgp_module_identity(),
    reporting_source=as.list(tools::md5sum(c(file.path(a$repo,'paper/applications/adni.Rmd'),file.path(a$repo,'eivGP/R/helpers.R'))))),
    file.path(dest,'REPORT_PROVENANCE.json'))
}
adni_metric_table <- function(a,protocol='OI',full=FALSE) {
  d<-a$prediction[a$prediction$protocol==protocol,,drop=FALSE]
  metrics<-c('RMSE','RMSE_SD','MAE_SD','CRPS_SD','Coverage95','Width95_SD','IntervalScore95_SD')
  tab<-adni_summary(d,c('cell','method'),metrics)
  if(!nrow(tab))return('<p>No available saved scores for this protocol.</p>')
  tab$CORE<-vapply(seq_len(nrow(tab)),function(i){
    if(tab$method[i]!='EIV-GP')return('N/A')
    q<-a$convergence[a$convergence$cell==tab$cell[i],];paste0(sum(q$core_pass),'/',nrow(q))},'')
  method<-tab$method
  for(i in seq_len(nrow(tab))) {
    flagged<-if(method[i]=='EIV-GP')tab$CORE[i]!=paste0(tab$folds[i],'/',tab$folds[i])else
      any(a$method_status$cell==tab$cell[i]&a$method_status$method==method[i]&a$method_status$state=='OPTIMIZER_FLAGGED')
    if(flagged)method[i]<-paste0(method[i],'*')
  }
  result<-data.frame(Method=method,Calibration=tab$cell,
    Test=if(protocol=='OI')'(X,C)'else ifelse(tab$method=='EIV-GP','(X,U)','(X,U,C)'),
    Folds=tab$folds,check.names=FALSE)
  shown<-if(full)metrics else c('RMSE','RMSE_SD','CRPS_SD','Coverage95','Width95_SD')
  labels<-c(RMSE='RMSE',RMSE_SD='RMSE/sY',MAE_SD='MAE/sY',CRPS_SD='CRPS/sY',
    Coverage95='Coverage (%)',Width95_SD='Width/sY',IntervalScore95_SD='Interval score/sY')
  best_mask<-matrix(FALSE,nrow(result),length(shown),dimnames=list(NULL,unname(labels[shown])))
  for(k in shown) {
    scale<-if(k=='Coverage95')100 else 1
    values<-sprintf(if(k=='Coverage95')'%.2f (%.2f)'else'%.4f (%.4f)',
      scale*tab[[paste0(k,'_mean')]],scale*tab[[paste0(k,'_SD')]])
    if(k%in%c('RMSE','RMSE_SD','MAE_SD','CRPS_SD','IntervalScore95_SD')) {
      strata<-if(protocol=='OI')rep('shared',nrow(tab))else tab$cell
      for(g in unique(strata)) {
        ix<-which(strata==g);best<-min(tab[[paste0(k,'_mean')]][ix])
        hit<-ix[tab[[paste0(k,'_mean')]][ix]==best]
        # Unconverged methods remain eligible; the asterisk records their status.
        best_mask[hit,labels[[k]]]<-TRUE
      }
    }
    result[[labels[[k]]]]<-values
  }
  result$CORE<-tab$CORE
  attr(result,'best')<-best_mask
  adni_table_html(result,if(protocol=='OI')'Ordinal-input prediction: fold mean (sample SD)'else
    'Complete-case prediction: fold mean (sample SD)')
}
adni_u_table <- function(d) {
  if(!nrow(d))return(data.frame())
  tab<-adni_summary(d,c('cell','coordinate'),c('RMSE','RMSE_SD','MAE_SD','CRPS_SD','Coverage95','Width95_SD','IntervalScore95_SD'))
  tab$gain_vs_C1_percent<-vapply(seq_len(nrow(tab)),function(i){
    v<-d[d$cell==tab$cell[i]&d$coordinate==tab$coordinate[i],]
    100*(1-sqrt(mean(v$RMSE^2))/sqrt(mean(v$C1_mean_RMSE^2)))},0)
  out<-tab[,c('cell','coordinate','folds')]
  for(k in c('RMSE','RMSE_SD','MAE_SD','CRPS_SD','Coverage95','Width95_SD','IntervalScore95_SD'))
    out[[k]]<-sprintf('%.4f (%.4f)',tab[[paste0(k,'_mean')]],tab[[paste0(k,'_SD')]])
  out$gain_vs_C1_percent<-sprintf('%+.2f',tab$gain_vs_C1_percent);out
}
adni_target_table <- function(a,all_folds=FALSE) {
  if(is.null(a$functionals))return(data.frame())
  d<-a$functionals$precision
  if(!all_folds)d<-d[d$fold==a$settings$representative_fold,,drop=FALSE]
  z<-a$input$cells[[1]];target<-rep('m(x,c)',nrow(d))
  for(j in seq_along(z$u)) {
    name<-paste0('Mean hidden ',a$protocol$U_order[j]);at<-d$target==name
    label<-if(grepl('FDG',a$protocol$U_order[j]))'FDG'else if(grepl('HIP',a$protocol$U_order[j]))'HIP1000/ICV'else a$protocol$U_order[j]
    target[at]<-paste0(label,' mean (',ifelse(d$cell[at]=='O0','model','raw'),')')
  }
  d$Target<-target
  d<-d[order(match(d$target,unique(d$target)),match(d$cell,c('O0','O20','O50')),d$fold),,drop=FALSE]
  out<-data.frame(Target=d$Target,Calibration=d$cell,`Mean (SD)`=sprintf('%.4f (%.4f)',d$mean,d$posterior_SD),
    `R-hat`=sprintf('%.3f',d$Rhat),`Bulk ESS`=sprintf('%.1f',d$bulk_ESS),`Tail ESS`=sprintf('%.1f',d$tail_ESS),
    MCSE=ifelse(abs(d$MCSE)<1e-4,formatC(d$MCSE,format='e',digits=2),sprintf('%.4f',d$MCSE)),check.names=FALSE)
  if(all_folds)out<-data.frame(Fold=d$fold,out,check.names=FALSE)
  out
}
adni_integration_table <- function(a,all_folds=FALSE) {
  if(is.null(a$functionals))return(data.frame())
  d<-a$functionals$integration
  if(!all_folds)d<-d[d$fold==a$settings$representative_fold,,drop=FALSE]
  data.frame(Fold=d$fold,Calibration=d$cell,Nodes=d$n_latent,`Check nodes`=d$n_latent_check,
    `Mean at nodes`=sprintf('%.5f',d$mean),`Mean at check nodes`=sprintf('%.5f',d$mean_check),
    Difference=sprintf('%+.5f',d$difference),check.names=FALSE)
}
adni_report_plots <- function(a) {
  library(ggplot2);library(patchwork)
  z<-a$input$cells[['1/O50']];d<-a$cohort;p<-a$protocol;plots<-list()
  long<-adni_bind(lapply(seq_along(z$u),function(j)data.frame(U=d[[z$u[j]]],Y=d$Y,C1=factor(d$C1),coordinate=p$U_order[j])))
  plots$cohort<-ggplot(long,aes(U,Y,colour=C1))+geom_point(alpha=.55,size=1.3)+
    facet_wrap(~coordinate,scales='free_x')+labs(x='Observed imaging reference',y='ADNI memory score',colour=p$C_order[1])+
    theme_minimal(base_size=11)+theme(legend.position='bottom')
  xlong<-adni_bind(lapply(seq_along(z$x),function(j)data.frame(X=d[[z$x[j]]],Y=d$Y,C1=factor(d$C1),coordinate=p$X_order[j])))
  plots$covariates<-ggplot(xlong,aes(X,Y,colour=C1))+geom_point(alpha=.55,size=1.3)+
    facet_wrap(~coordinate,scales='free_x')+labs(x='Observed covariate',y='ADNI memory score',colour=p$C_order[1])+
    theme_minimal(base_size=11)+theme(legend.position='bottom')
  shown<-seq_len(min(2L,length(z$c)))
  cu<-adni_bind(lapply(shown,function(k)adni_bind(lapply(seq_along(z$u),function(j)
    data.frame(grade=factor(d[[z$c[k]]],levels=seq_len(p$m_vec[k])),U=d[[z$u[j]]],
      item=paste0('C',k,': ',p$C_order[k]),coordinate=p$U_order[j])))))
  plots$C_U<-ggplot(cu,aes(grade,U))+geom_boxplot(width=.55,outlier.shape=NA,fill='grey90',colour='grey35')+
    geom_point(aes(colour=grade),position=position_jitter(width=.12,height=0,seed=6203L),alpha=.45,size=.9)+
    facet_wrap(~item+coordinate,ncol=length(z$u),scales='free_y',labeller=label_wrap_gen(36))+
    labs(x='Native ECog grade (higher = greater functional difficulty)',y='Measured imaging reference (raw units)')+
    theme_minimal(base_size=11)+theme(legend.position='none')
  counts<-adni_bind(lapply(seq_along(z$c),function(j){t<-as.data.frame(table(factor(d[[z$c[j]]],levels=seq_len(p$m_vec[j]))));
    data.frame(item=p$C_order[j],grade=t$Var1,N=t$Freq)}))
  plots$C_counts<-ggplot(counts,aes(grade,N,fill=grade))+geom_col()+facet_wrap(~item,ncol=3,scales='free_y')+
    labs(x='Native ECog grade',y='Participants')+theme_bw()+theme(legend.position='none')
  u<-a$U_predictions
  if(nrow(u)) {
    u$role<-ifelse(u$calibration==1,'Observed calibration','Masked training U')
    u$C1<-factor(u$C1)
    draw_u<-function(x,all=FALSE){
      ggplot(x,aes(true_std,mean_std,colour=C1))+geom_abline(slope=1,intercept=0,linetype=2,colour='grey50')+
        geom_linerange(data=x[x$calibration==0,],aes(ymin=lo_std,ymax=hi_std),alpha=.15)+
        geom_point(aes(shape=role),size=1,alpha=.6)+scale_shape_manual(values=c('Masked training U'=16,'Observed calibration'=1))+
        facet_grid(if(all)fold~coordinate else cell~coordinate,scales='free')+
        labs(x='Reference U (training SD units; evaluation only)',y='Posterior mean (training SD units)',colour=p$C_order[1],shape=NULL)+theme_bw()
    }
    plots$imputation<-draw_u(u[u$fold==a$settings$representative_fold,])
    for(cell in unique(u$cell))plots[[paste0('imputation_',cell)]]<-draw_u(u[u$cell==cell,],TRUE)
    hidden<-u[u$calibration==0,]
    density<-adni_bind(list(data.frame(hidden[,c('cell','coordinate')],value=hidden$true_std,distribution='Reference'),
      data.frame(hidden[,c('cell','coordinate')],value=hidden$mean_std,distribution='Posterior means')))
    plots$input_distribution<-ggplot(density,aes(value,colour=distribution,fill=distribution))+
      geom_density(alpha=.12)+facet_grid(cell~coordinate,scales='free')+
      labs(x='Training SD units',y='Density',colour=NULL,fill=NULL)+theme_bw()
    class_rows<-adni_bind(lapply(split(hidden,interaction(hidden$fold,hidden$cell,hidden$coordinate,hidden$C1,drop=TRUE)),function(x)
      data.frame(fold=x$fold[1],cell=x$cell[1],coordinate=x$coordinate[1],C1=as.character(x$C1[1]),N=nrow(x),
        RMSE_SD=sqrt(mean((x$mean_std-x$true_std)^2)),Coverage95=mean(x$true_U>=x$q025&x$true_U<=x$q975))))
    plots$class_table<-adni_summary(class_rows,c('cell','coordinate','C1'),c('RMSE_SD','Coverage95'))
    common<-hidden[as.character(hidden$common_hidden)%in%c('TRUE','True','1'),]
    common_rows<-adni_bind(lapply(split(common,interaction(common$fold,common$cell,common$coordinate,drop=TRUE)),function(x)
      data.frame(fold=x$fold[1],cell=x$cell[1],coordinate=x$coordinate[1],N=nrow(x),RMSE_SD=sqrt(mean((x$mean_std-x$true_std)^2)))))
    plots$common_hidden<-adni_summary(common_rows,c('cell','coordinate'),'RMSE_SD')
    plots$support<-adni_bind(lapply(split(hidden,interaction(hidden$cell,hidden$coordinate,drop=TRUE)),function(x)
      data.frame(cell=x$cell[1],coordinate=x$coordinate[1],fold_rows=nrow(x),lower_nonpositive=sum(x$q025<=0),
        mean_nonpositive=sum(x$posterior_mean<=0),interval_outside_C_grade='Not applicable: C is not an imaging bin')))
  }
  if(nrow(u))plots<-modifyList(plots,adni_input_distributions(a))
  plots
}

adni_input_distributions <- function(a) {
  library(ggplot2);parts<-list()
  for(key in names(a$paths)) {
    path<-a$paths[[key]];if(path$cell=='O0'||path$fold!=a$settings$representative_fold)next
    z<-a$input$cells[[key]]
    hidden<-which(a$input$cells[[paste(path$fold,'O50',sep='/')]]$train$calibration==0L)
    saved<-readRDS(file.path(path$dir,'U_draws.rds'))
    if(!identical(saved$scale,'raw'))stop('Raw saved U draws required for physical distribution checks')
    ids<-unique(round(seq(1,dim(saved$draws)[1],length.out=min(128,dim(saved$draws)[1]))))
    for(j in seq_along(z$u)) {
      ref<-z$train[[z$u[j]]];mu<-mean(ref);ss<-sd(ref)
      draw<-matrix(saved$draws[ids,hidden,j],nrow=length(ids))
      tag<-data.frame(cell=path$cell,coordinate=z$protocol$U_order[j])
      parts[[length(parts)+1L]]<-data.frame(tag,C1=rep(as.character(z$train$C1[hidden]),each=length(ids)),
        value=as.vector((draw-mu)/ss),distribution='Posterior mixture')
      parts[[length(parts)+1L]]<-data.frame(tag,C1=as.character(z$train$C1[hidden]),value=(ref[hidden]-mu)/ss,distribution='Reference')
      means<-colMeans(matrix(saved$draws[,hidden,j],ncol=length(hidden)))
      parts[[length(parts)+1L]]<-data.frame(tag,C1=as.character(z$train$C1[hidden]),value=(means-mu)/ss,distribution='Posterior means')
    }
  }
  d<-adni_bind(parts);if(!nrow(d))return(list())
  d$C1<-factor(d$C1);reference<-d[d$distribution=='Reference',];post<-d[d$distribution!='Reference',]
  draw_distribution<-function(class=FALSE) {
    g<-ggplot()+geom_histogram(data=reference,aes(value,after_stat(density)),bins=15,
      fill='grey80',colour='white')+geom_rug(data=reference,aes(value),alpha=.3)+
      geom_density(data=post,aes(value,colour=distribution,linetype=distribution),linewidth=.7)+
      scale_linetype_manual(values=c('Posterior mixture'='solid','Posterior means'='dashed'))+
      theme_minimal(base_size=11)+theme(legend.position='bottom')+
      labs(x='Common-hidden U (training SD units)',y='Density',colour=NULL,linetype=NULL)
    if(class)g+facet_grid(cell+coordinate~C1,scales='free')else g+facet_grid(cell~coordinate,scales='free')
  }
  list(input_distribution=draw_distribution(),class_distribution=draw_distribution(TRUE))
}

adni_paper_assets <- function(a) {
  dest<-file.path(a$report_dir,'report-assets','paper');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  plots<-list(fig_adni_cohort=a$plots$cohort,fig_adni_C_support=a$plots$C_counts,
    fig_adni_Y_X=a$plots$covariates,fig_adni_C_U=a$plots$C_U,
    fig_adni_target_traces=a$functionals$trace_plot,fig_adni_mean_functions=a$functionals$mean_plot,
    fig_adni_response_surfaces=a$functionals$surface_plot,fig_adni_U_imputation=a$plots$imputation,
    fig_adni_U_distributions=a$plots$input_distribution,fig_adni_C1_U_distributions=a$plots$class_distribution)
  manifest<-list()
  for(name in names(plots)) {
    file<-paste0(name,'.pdf');available<-!is.null(plots[[name]])
    if(available)ggplot2::ggsave(file.path(dest,file),plots[[name]],width=12,height=if(grepl('support|traces|C1',name))9 else 6,
      device=grDevices::pdf,limitsize=FALSE)
    manifest[[length(manifest)+1L]]<-data.frame(label=name,file=file,status=if(available)'available'else'not_applicable',
      N=a$input$N,profile=a$settings$profile,MD5=if(available)unname(tools::md5sum(file.path(dest,file)))else'')
  }
  oi<-a$prediction[a$prediction$protocol=='OI',,drop=FALSE]
  cc<-a$prediction[a$prediction$protocol=='CC-known',,drop=FALSE]
  metrics<-c('RMSE','RMSE_SD','MAE_SD','CRPS_SD','Coverage95','Width95_SD','IntervalScore95_SD')
  tables<-list(tab_adni_design=a$input$audit,
    tab_adni_C_definitions=adni_C_dictionary(a$protocol$C_order,a$protocol$m_vec,a$metadata$C_labels),
    tab_adni_prediction_OI=adni_summary(oi,c('cell','method'),metrics),
    tab_adni_prediction_CC=adni_summary(cc,c('cell','method'),metrics),tab_adni_training_U=adni_u_table(a$U_metrics),
    tab_adni_prospective_U=adni_u_table(a$prospective_metrics),tab_adni_convergence=a$convergence)
  for(name in names(tables)) {
    d<-tables[[name]];available<-nrow(d)>0L;file<-paste0(name,'.csv')
    if(available) {
      write.csv(d,file.path(dest,file),row.names=FALSE)
      writeLines(adni_table_html(d,name),file.path(dest,paste0(name,'.html')))
    }
    manifest[[length(manifest)+1L]]<-data.frame(label=name,file=file,status=if(available)'available'else'unavailable',
      N=a$input$N,profile=a$settings$profile,MD5=if(available)unname(tools::md5sum(file.path(dest,file)))else'')
  }
  out<-adni_bind(manifest);write.csv(out,file.path(dest,'artifact_manifest.csv'),row.names=FALSE);out
}


# ============================================================================
# ADNI prepared-fold application: adni-functionals.R
# ============================================================================

# Author: Xinye Chen. Read-only scientific targets and conditional response slices.
adni_trace_diagnostic <- function(values,chain) {
  groups<-split(values,chain);n<-min(lengths(groups));x<-do.call(cbind,lapply(groups,head,n))
  data.frame(states_per_chain=n,chains=ncol(x),mean=mean(values),posterior_SD=sd(values),
    Rhat=posterior::rhat(x),bulk_ESS=posterior::ess_bulk(x),tail_ESS=posterior::ess_tail(x),MCSE=posterior::mcse_mean(x))
}
adni_baseline_mean <- function(api,dir,X,C,mvec,method) {
  a<-readRDS(file.path(dir,'fit.rds'));x<-sweep(sweep(X,2,a$center,'-'),2,a$scale,'/')
  if(method=='UC-GP') {
    nd<-getFromNamespace('make_mixedgp_competitor_data','eivGP')(x,C,mvec)
    mean<-as.numeric(stats::predict(a$model$fit,newdata=nd,type='UK',seCompute=FALSE,lightReturn=TRUE)$mean)
  }else if(method=='LVGP')mean<-as.numeric(LVGP::LVGP_predict(cbind(x,C),a$model$fit,MSE_on=0)$Y_hat)
  else mean<-api$adni_ezgp_predict_covariance_consistent(api$adni_ezgp_mixed_matrix(x,C,mvec),a$model)$Y_hat
  a$y_center+a$y_scale*mean
}
adni_saved_functionals <- function(a) {
  library(ggplot2)
  dest<-file.path(a$report_dir,'report-assets','functionals');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  traces<-raw_traces<-curves<-surfaces<-precision<-integration<-list()
  code<-file.path(a$repo,'eivGP/R/helpers.R')
  summarize_draws<-function(x)data.frame(mean=colMeans(x),lo=apply(x,2,quantile,.025),hi=apply(x,2,quantile,.975))
  for(key in names(a$paths)) {
    path<-a$paths[[key]];f<-path$fold;cell<-path$cell;z<-a$input$cells[[key]];p<-z$protocol
    cache<-file.path(dest,paste0('fold',f,'_',cell,'.rds'))
    identity<-list(fit=unname(tools::md5sum(path$fit)),inputs=unname(tools::md5sum(a$input$files)),
      code=unname(tools::md5sum(code)),states=a$settings$functional_draws_per_chain,n_latent=a$settings$n_latent,
      analysis_modules=a$api$eivgp_module_identity())
    if(file.exists(cache)) {
      saved<-readRDS(cache)
      if(!identical(saved$identity,identity))stop('Derived functional cache differs; choose a fresh report directory: ',cache)
      result<-saved$result
    }else {
      cat('Read-only functionals fold',f,cell,'\n');flush.console()
      fit<-a$api$adni_plan_read_fit(path$fit)
      stopifnot(if(is.null(fit$scientific_ess$settings$ess_min))fit$checkpoint$iteration==a$numerics$iterations else fit$checkpoint$iteration>=a$numerics$iterations,fit$checkpoint$burn==a$numerics$burn)
      info<-fit$mcmc$mcmc_draw_info
      ids<-unlist(lapply(split(seq_len(nrow(info)),info$chain),function(v)
        v[unique(round(seq(1,length(v),length.out=min(length(v),a$settings$functional_draws_per_chain))))]),use.names=FALSE)
      chains<-info$chain[ids]
      xref<-vapply(z$x,function(n)median(z$train[[n]]),0)
      patterns<-apply(z$train[,z$c,drop=FALSE],1,paste,collapse='/');mode<-sort(names(table(patterns))[table(patterns)==max(table(patterns))])[1]
      cref<-as.numeric(strsplit(mode,'/',fixed=TRUE)[[1]])
      age<-seq(quantile(z$train[[z$x[1]]],.1),quantile(z$train[[z$x[1]]],.9),length.out=7)
      grid<-expand.grid(x=age,C1=seq_len(p$m_vec[1]))
      xx<-matrix(rep(xref,each=nrow(grid)),ncol=length(xref));xx[,1]<-grid$x
      cc<-matrix(rep(cref,each=nrow(grid)),ncol=length(cref));cc[,1]<-grid$C1;colnames(cc)<-p$C_order
      seed<-as.integer((as.double(p$seeds$posterior)+63001)%%(2^30-1)+1)
      md<-eivGP::predict_eivgp(fit,new_X=xx,new_C=cc,target='mean',draw_ids=ids,n_latent=a$settings$n_latent,
        include_process_uncertainty=TRUE,return_components=TRUE,seed=seed)
      cm<-attr(md,'conditional_means',exact=TRUE);stopifnot(is.matrix(cm))
      cv<-a$convergence[a$convergence$fold==f&a$convergence$cell==cell,]
      method_label<-if(nrow(cv)&&all(cv$core_pass))'EIV-GP'else'EIV-GP*'
      curve<-cbind(data.frame(fold=f,cell=cell,method=method_label,grid),summarize_draws(md))
      target_grade<-min(2L,p$m_vec[1]);target<-which(grid$C1==target_grade)[4L]
      mt<-data.frame(fold=f,cell=cell,target='m(X reference, C profile)',chain=chains,
        iteration=fit$checkpoint$burn+info$draw_within_chain[ids],value=cm[,target])
      pr<-cbind(data.frame(fold=f,cell=cell,target='m(X reference, C profile)',units='ADNI_MEM'),adni_trace_diagnostic(cm[,target],chains))
      # Repeat only the fixed target with twice as many integration nodes.
      md2<-eivGP::predict_eivgp(fit,new_X=xx[target,,drop=FALSE],new_C=cc[target,,drop=FALSE],target='mean',
        draw_ids=ids,n_latent=2L*a$settings$n_latent,include_process_uncertainty=FALSE,seed=seed)
      check<-data.frame(fold=f,cell=cell,n_latent=a$settings$n_latent,n_latent_check=2L*a$settings$n_latent,
        mean=mean(cm[,target]),mean_check=mean(md2),difference=mean(md2)-mean(cm[,target]))
      refcell<-a$input$cells[[paste(f,'O50',sep='/')]]
      hidden<-which(refcell$train$calibration==0L&refcell$train$C1==target_grade)
      if(!length(hidden))hidden<-which(refcell$train$calibration==0L)
      for(j in seq_along(z$u)) {
        v<-rowMeans(matrix(fit$mcmc$samples_U[,hidden,j],nrow=nrow(info)))
        anchored<-isTRUE(fit$data$latent_scale_anchored)
        if(anchored)v<-fit$data$U_center[j]+fit$data$U_scale[j]*v
        name<-paste0('Mean hidden ',p$U_order[j])
        mt<-rbind(mt,data.frame(fold=f,cell=cell,target=name,chain=info$chain,
          iteration=fit$checkpoint$burn+info$draw_within_chain,value=v))
        pr<-rbind(pr,cbind(data.frame(fold=f,cell=cell,target=name,
          units=if(anchored)'physical U'else'unanchored model U'),adni_trace_diagnostic(v,info$chain)))
      }
      surface<-data.frame()
      if(cell!='O0') {
        cal<-z$train[z$train$calibration==1,,drop=FALSE]
        ug<-lapply(z$u,function(n)seq(quantile(cal[[n]],.1),quantile(cal[[n]],.9),length.out=11))
        if(length(z$u)>=2L)uv<-expand.grid(U1=ug[[1]],U2=ug[[2]])else uv<-data.frame(U1=ug[[1]],U2=0)
        uu<-matrix(rep(vapply(z$u,function(n)median(cal[[n]]),0),each=nrow(uv)),ncol=length(z$u))
        uu[,1]<-uv$U1;if(length(z$u)>=2L)uu[,2]<-uv$U2;colnames(uu)<-p$U_order
        x<-matrix(rep(xref,each=nrow(uv)),ncol=length(xref))
        fd<-eivGP::predict_eivgp(fit,new_X=x,new_U=uu,new_U_scale='raw',target='surface',draw_ids=ids,
          include_process_uncertainty=TRUE,seed=seed+3000L)
        surface<-cbind(data.frame(fold=f,cell=cell,uv),summarize_draws(fd))
      }
      raw<-adni_bind(lapply(c('V','sigma2','r'),function(n){
        v<-fit$mcmc[[paste0('samples_',n)]]
        if(is.null(v))return(data.frame())
        data.frame(fold=f,cell=cell,target=n,chain=info$chain,iteration=fit$checkpoint$burn+info$draw_within_chain,value=as.numeric(v))
      }))
      result<-list(curves=curve,surfaces=surface,traces=mt,raw_traces=raw,precision=pr,integration=check,
        reference=list(X=xref,C=cref,target_C1=target_grade,common_hidden_N=length(hidden)))
      saveRDS(list(identity=identity,result=result),cache)
      rm(fit,md,md2);gc(FALSE)
    }
    curves[[key]]<-result$curves;surfaces[[key]]<-result$surfaces;traces[[key]]<-result$traces
    raw_traces[[key]]<-result$raw_traces;precision[[key]]<-result$precision;integration[[key]]<-result$integration
  }
  # Native observable-mean curves on the same fixed training-defined slices.
  native<-list()
  for(f in a$settings$folds) {
    z<-a$input$cells[[paste(f,'O50',sep='/')]];p<-z$protocol
    grid<-curves[[paste(f,a$settings$cells[1],sep='/')]][,c('x','C1')]
    xref<-vapply(z$x,function(n)median(z$train[[n]]),0)
    patterns<-apply(z$train[,z$c,drop=FALSE],1,paste,collapse='/');mode<-sort(names(table(patterns))[table(patterns)==max(table(patterns))])[1]
    cref<-as.numeric(strsplit(mode,'/',fixed=TRUE)[[1]])
    xx<-matrix(rep(xref,each=nrow(grid)),ncol=length(xref));xx[,1]<-grid$x
    cc<-matrix(rep(cref,each=nrow(grid)),ncol=length(cref));cc[,1]<-grid$C1
    for(method in c('UC-GP','LVGP','EzGP')) {
      status<-a$method_status[a$method_status$fold==f&a$method_status$cell=='OI'&a$method_status$method==method,]
      if(!nrow(status)||status$state=='MISSING')next
      value<-adni_baseline_mean(a$api,status$source_directory,xx,cc,p$m_vec,method)
      native[[length(native)+1L]]<-data.frame(fold=f,cell='native OI',
        method=paste0(method,if(status$state=='OPTIMIZER_FLAGGED')'*'else''),grid,mean=value,lo=NA_real_,hi=NA_real_)
    }
  }
  out<-list(curves=adni_bind(c(curves,native)),surfaces=adni_bind(surfaces),traces=adni_bind(traces),
    raw_traces=adni_bind(raw_traces),precision=adni_bind(precision),integration=adni_bind(integration))
  for(n in names(out))if(nrow(out[[n]]))write.csv(out[[n]],file.path(dest,paste0(n,'.csv')),row.names=FALSE)
  reference_fold<-a$settings$representative_fold
  plot_curve<-function(x)ggplot(x,aes(x,mean,colour=interaction(method,cell),group=interaction(method,cell)))+
    geom_ribbon(aes(ymin=lo,ymax=hi,fill=interaction(method,cell)),alpha=.06,colour=NA,na.rm=TRUE)+
    geom_line(linewidth=.65)+facet_wrap(~C1,ncol=2)+theme_bw()+
    labs(x=a$protocol$X_order[1],y='m(x,c): ADNI_MEM',colour='Method / calibration',fill='Method / calibration')
  out$mean_plot<-plot_curve(out$curves[out$curves$fold==reference_fold,])
  out$mean_all_plot<-plot_curve(out$curves)+facet_grid(fold~C1)
  if(nrow(out$surfaces)) {
    out$surface_plot<-ggplot(out$surfaces[out$surfaces$fold==reference_fold,],aes(U1,U2,fill=mean))+
      geom_raster()+geom_contour(aes(z=mean),colour='white',alpha=.6)+facet_wrap(~cell)+
      scale_fill_viridis_c()+theme_bw()+labs(x=a$protocol$U_order[1],y=a$protocol$U_order[min(2,a$protocol$d)],fill='f(x,u)')
    out$surface_all_plot<-ggplot(out$surfaces,aes(U1,U2,fill=mean))+geom_raster()+
      facet_grid(fold~cell)+scale_fill_viridis_c()+theme_bw()+labs(x=a$protocol$U_order[1],y=a$protocol$U_order[min(2,a$protocol$d)],fill='f(x,u)')
  }
  traceplot<-function(d)ggplot(d[d$fold==reference_fold,],aes(iteration,value,colour=factor(chain)))+
    geom_line(alpha=.6,linewidth=.25)+facet_grid(target~cell,scales='free_y')+theme_bw()+
    labs(x='Saved iteration',y=NULL,colour='Chain')
  out$trace_plot<-traceplot(out$traces);out$raw_trace_plot<-traceplot(out$raw_traces)
  out
}


# ============================================================================
# ADNI prepared-fold application: adni-render.R
# ============================================================================

# Author: Xinye Chen. Fresh-process plotting avoids macOS font initialization after fork.
adni_render_plot <- function(args=commandArgs(TRUE)) {
  stopifnot(length(args)==5L)
  plot<-readRDS(args[1])
  ggplot2::ggsave(args[2],plot,width=as.numeric(args[3]),height=as.numeric(args[4]),
    dpi=as.numeric(args[5]),bg='white')
  stopifnot(file.exists(args[2]),file.info(args[2])$size>1000)
  invisible(args[2])
}


# ============================================================================
# ADNI prepared-fold application: adni.R
# ============================================================================

# Author: Xinye Chen. Application settings and routing; scientific kernels are eivGP APIs.
adni_resolve_path <- function(path,repo) {
  path<-path.expand(path)
  if(!grepl('^(/|[A-Za-z]:[/\\\\])',path))path<-file.path(repo,path)
  base<-path;suffix<-character()
  while(!file.exists(base)&&dirname(base)!=base) {
    suffix<-c(basename(base),suffix);base<-dirname(base)
  }
  resolved<-normalizePath(base,mustWork=TRUE)
  if(length(suffix))file.path(resolved,paste(suffix,collapse='/'))else resolved
}
adni_cli_params <- function(args) {
  allowed<-c('cores','profile','dataset','data_root','data_dir','output','output_dir','config','manifest',
    'chains','iter','burn','stages','folds','cells','draws_per_chain','first_checkpoint','checkpoint_step',
    'baseline_profile','sampler_o0','sampler_o20','sampler_o50','representative_fold','appendix',
    'plan','reuse_only','report_source','report_dir','functional_draws_per_chain','n_latent','run_id','ess_min')
  p<-list()
  for(a in args) {
    if(!grepl('^--[^=]+=.+$',a))stop('Use --name=value: ',a)
    k<-gsub('-','_',sub('^--([^=]+)=.*$','\\1',a));v<-sub('^[^=]+=','',a)
    if(!k%in%allowed)stop('Unknown ADNI option: ',k)
    if(k%in%names(p))stop('Repeated option: ',k)
    if(k%in%c('plan','reuse_only','appendix')) {
      if(!v%in%c('true','false'))stop(k,' must be true or false')
      v<-v=='true'
    }
    p[[k]]<-v
  }
  p
}
adni_integer <- function(value,name,minimum=1L) {
  n<-suppressWarnings(as.numeric(value))
  if(length(n)!=1L||!is.finite(n)||n<minimum||n!=floor(n)||n>.Machine$integer.max)
    stop(name,' must be an integer >= ',minimum)
  as.integer(n)
}
adni_integer_list <- function(value,name) {
  if(is.character(value)&&length(value)==1L) {
    if(grepl('^[0-9]+:[0-9]+$',value)) {
      z<-as.integer(strsplit(value,':',fixed=TRUE)[[1]]);value<-seq.int(z[1],z[2])
    }else value<-strsplit(value,',',fixed=TRUE)[[1]]
  }
  unname(vapply(unlist(value),adni_integer,1L,name=name))
}
adni_settings <- function(profile='pilot',cores=4L,...) {
  overrides<-list(...);profile<-match.arg(profile,c('pilot','paper','quick'))
  defaults<-switch(profile,paper=list(chains=4L,iter=50000L,burn=5000L,stages=50000L),
    pilot=list(chains=4L,iter=1500L,burn=500L,stages=c(1000L,1500L)),
    quick=list(chains=1L,iter=80L,burn=30L,stages=80L))
  s<-modifyList(c(defaults,list(profile=profile,cores=cores,folds=if(profile=='quick')1L else 1:5,
    cells=c('O50','O20','O0'),draws_per_chain=if(profile=='quick')10L else 250L,
    checkpoint_step=1000L,baseline_profile=if(profile=='quick')'quick'else'paper',
    sampler_o0='HMC',sampler_o20='HMC',sampler_o50='ESS',representative_fold=1L,
    appendix=TRUE,functional_draws_per_chain=64L,n_latent=256L,run_id='local')),overrides)
  s$ess_min<-mixedgp_parse_ess_min(s$ess_min)
  if(!is.null(s$ess_min) && as.numeric(s$chains)<4L)stop('Scientific ESS requires at least four chains.')
  for(n in c('cores','chains','iter','draws_per_chain','checkpoint_step','representative_fold',
             'functional_draws_per_chain','n_latent'))s[[n]]<-adni_integer(s[[n]],n)
  s$burn<-adni_integer(s$burn,'burn',0L);s$folds<-adni_integer_list(s$folds,'folds')
  if(any(!s$folds%in%1:5)||anyDuplicated(s$folds))stop('folds must be distinct members of1:5')
  s$cells<-if(length(s$cells)==1L)strsplit(s$cells,',',fixed=TRUE)[[1]]else unlist(s$cells)
  if(!length(s$cells)||any(!s$cells%in%c('O0','O20','O50'))||anyDuplicated(s$cells))stop('Invalid cells')
  if(!is.null(overrides$stages)) {
    s$stages<-adni_integer_list(overrides$stages,'stages')
    if(!is.null(overrides$iter)&&tail(s$stages,1)!=s$iter)stop('iter must equal the last stage')
    s$iter<-tail(s$stages,1)
  }else if(!is.null(overrides$iter))s$stages<-s$iter
  if(any(diff(s$stages)<=0)||min(s$stages)<=s$burn)stop('Stages must increase and exceed burn')
  s$first_checkpoint<-adni_integer(if(is.null(s$first_checkpoint))min(s$stages[1],s$burn+250L)else s$first_checkpoint,'first_checkpoint')
  if(s$first_checkpoint<=s$burn||s$first_checkpoint>min(s$stages))stop('first-checkpoint must exceed burn and not exceed first stage')
  if(!s$baseline_profile%in%c('paper','quick'))stop('baseline-profile must be paper or quick')
  for(n in c('sampler_o0','sampler_o20','sampler_o50'))if(!s[[n]]%in%c('HMC','ESS'))stop('Unknown sampler')
  attr(s,'explicit_selection')<-intersect(names(overrides),c('folds','cells'))
  s
}
adni_load_api <- function(repo) adni_workflow_api(repo)
adni_method_run <- function(repo,params,group) {
  if(!is.null(params$config)) {
    config<-jsonlite::fromJSON(adni_resolve_path(params$config,repo),simplifyVector=TRUE)
    params<-modifyList(config,params)
  }
  get<-function(n,default)if(is.null(params[[n]]))default else params[[n]]
  if(!is.null(params$manifest)) {
    if(!is.null(params$ess_min)) {
      declared<-jsonlite::fromJSON(adni_resolve_path(params$manifest,repo))$ess_min
      if(!identical(mixedgp_parse_ess_min(params$ess_min),mixedgp_parse_ess_min(declared)))
        stop('--ess-min must match the supplied manifest; set ess_min in that manifest.')
    }
    env<-adni_load_api(repo)
    return(env$adni_screen_run(repo,list(manifest=adni_resolve_path(params$manifest,repo),
      cores=adni_integer(get('cores',4L),'cores'),plan=isTRUE(params$plan),reuse_only=isTRUE(params$reuse_only)),group))
  }
  fields<-c('chains','iter','burn','stages','folds','cells','draws_per_chain','first_checkpoint','checkpoint_step',
    'baseline_profile','sampler_o0','sampler_o20','sampler_o50','representative_fold','appendix',
    'functional_draws_per_chain','n_latent','run_id','ess_min')
  s<-do.call(adni_settings,c(list(profile=get('profile','pilot'),cores=get('cores',4L)),params[intersect(names(params),fields)]))
  dataset<-get('dataset','adni')
  if(!grepl('^[A-Za-z0-9][A-Za-z0-9_-]*$',dataset))stop('Invalid dataset name')
  default_input<-if(dataset=='adni')file.path('data','adni','inputs')else file.path('data','adni',dataset,'inputs')
  input<-adni_resolve_path(get('data_root',get('data_dir',default_input)),repo)
  out<-adni_resolve_path(get('output',get('output_dir',file.path('results',paste(dataset,s$profile,sep='_')))),repo)
  if(isTRUE(params$plan)&&!dir.exists(file.path(input,'MC'))) {
    cat('Dataset:',dataset,'; data:',input,'\n');print(s[c('profile','chains','iter','burn','cores')])
    cat('Supply private inputs before fitting. No computation started.\n');return(invisible(s))
  }
  report_source<-if(!is.null(params$report_source))adni_resolve_path(params$report_source,repo)else NULL
  saved_only<-group=='report'||isTRUE(params$reuse_only)
  if(saved_only&&!isTRUE(params$plan))return(adni_render_report(repo,input,out,s,report_source,
    report_dir=adni_resolve_path(get('report_dir',out),repo),group=group))
  api<-adni_load_api(repo);inp<-api$adni_plan_inputs(input);p<-inp$cells[[1]]$protocol
  api$application$adni_output_check(out,repo)
  files<-setNames(lapply(inp$files,function(x)digest::digest(file=x,algo='sha256')),inp$files)
  m<-list(schema='adni-application-v1',author='Xinye Chen',profile=s$profile,
    package_version=as.character(packageVersion('eivGP')),
    baseline_profile=s$baseline_profile,trial_id=dataset,smoke=s$profile=='quick',
    output_root=file.path(out,'workflow'),cores=s$cores,chains=s$chains,stages=s$stages,burn=s$burn,
    first_checkpoint=s$first_checkpoint,checkpoint_step=s$checkpoint_step,draws_per_chain=s$draws_per_chain,ess_min=s$ess_min,
    samplers=list(O0=s$sampler_o0,O20=s$sampler_o20,O50=s$sampler_o50),folds=s$folds,cells=s$cells,
    cases=list(list(id=p$plan_id,input=input,output=out,C=p$C_order,N=inp$N,q=p$q)),frozen_inputs=files)
  if(isTRUE(params$plan)) {
    e<-adni_load_api(repo)
    print(e$adni_screen_tasks(m,group),row.names=FALSE);print(s[c('profile','chains','iter','burn','cores')]);return(invisible(m))
  }
  dir.create(out,recursive=TRUE,showWarnings=FALSE);manifest<-file.path(out,'application_manifest.json')
  if(file.exists(manifest)) {
    prior<-jsonlite::fromJSON(manifest,simplifyVector=FALSE);a<-m;b<-prior;a$cores<-b$cores<-NULL
    canonical<-function(x)jsonlite::fromJSON(jsonlite::toJSON(x,auto_unbox=TRUE,null='null'),simplifyVector=FALSE)
    if(!identical(canonical(a),canonical(b)))
      stop('Saved ADNI settings/inputs differ; select a new output directory. ',
        paste(all.equal(canonical(a),canonical(b)),collapse='; '))
  }else api$atomic_json(m,manifest)
  e<-adni_load_api(repo)
  result<-e$adni_screen_run(repo,list(manifest=manifest,cores=s$cores,plan=FALSE,reuse_only=FALSE),group)
  adni_render_report(repo,input,out,s,report_dir=file.path(out,'reports',group),group=group)
  invisible(result)
}
adni_render_report <- function(repo,data_dir,output_dir,settings,report_source=NULL,
                              report_dir=output_dir,group='report',rmd=NULL,filename='adni.html') {
  if(!requireNamespace('rmarkdown',quietly=TRUE)||!rmarkdown::pandoc_available())
    stop('HTML requires rmarkdown and Pandoc on PATH or RSTUDIO_PANDOC')
  env<-adni_load_api(repo)
  analysis<-env$adni_report_prepare(repo,data_dir,output_dir,settings,report_source,report_dir,group)
  render_env<-new.env(parent=globalenv());render_env$.adni_precomputed_analysis<-analysis
  render_env$.adni_repository_dir<-repo
  if(is.null(rmd))rmd<-file.path(repo,'paper/applications/adni.Rmd')
  rmarkdown::render(rmd,output_file=filename,output_dir=report_dir,
    params=list(profile=settings$profile,action='report',cores=settings$cores,data_dir=data_dir,
      output_dir=output_dir,report_source=if(is.null(report_source))''else report_source,
      representative_fold=settings$representative_fold,appendix=settings$appendix),envir=render_env,quiet=TRUE)
  cat('Saved ADNI HTML: ',file.path(report_dir,filename),'\n',sep='')
  invisible(list(html=file.path(report_dir,filename),analysis=analysis))
}


# ============================================================================
# ADNI prepared-fold application: adni-workflow.R
# ============================================================================

# Author: Xinye Chen. ADNI screening orchestration on the unchanged public package.
adni_output_check <- function(path,repo) {
  path<-normalizePath(path.expand(path),mustWork=FALSE);repo<-normalizePath(repo)
  if(startsWith(paste0(path,'/'),paste0(repo,'/'))&&
     !startsWith(paste0(path,'/'),paste0(repo,'/results/')))
    stop('Choose ignored results/ or an external output directory; never write to versioned source/data.')
  invisible(path)
}
adni_screen_api <- function(repo) adni_workflow_api(repo)

adni_screen_args <- function(args) {
  p<-list(cores=5L,plan=FALSE,reuse_only=FALSE)
  for(a in args){
    if(!grepl('^--[^=]+=.+$',a))stop('Use --name=value: ',a)
    k<-gsub('-','_',sub('^--([^=]+)=.*$','\\1',a));v<-sub('^[^=]+=','',a)
    if(!k%in%c('manifest','cores','plan','reuse_only'))stop('Unknown screening option: ',k)
    if(k%in%c('plan','reuse_only')){if(!v%in%c('true','false'))stop(k,' must be true or false');v<-v=='true'}
    if(k=='cores'){n<-suppressWarnings(as.numeric(v));if(!is.finite(n)||n<1||n!=floor(n)||n>.Machine$integer.max)stop('cores must be a positive integer');v<-as.integer(n)}
    if(k%in%names(p)&&!k%in%c('cores','plan','reuse_only'))stop('Repeated option: ',k)
    p[[k]]<-v
  }
  if(is.null(p$manifest))stop('Supply --manifest=/external/SCREENING.json')
  p
}
adni_screen_settings <- function(m,iter,cores) {
  smoke<-isTRUE(m$smoke)
  get<-function(n,default)if(is.null(m[[n]]))default else m[[n]]
  list(profile=get('baseline_profile',if(smoke)'pilot'else'paper'),run_profile=get('profile',if(smoke)'quick'else'pilot'),
    action='fit',folds=as.integer(unlist(m$folds)),cells=unlist(m$cells),
    cores=as.integer(cores),chains=as.integer(get('chains',4L)),iter=as.integer(iter),burn=as.integer(m$burn),
    first_checkpoint=as.integer(min(iter,get('first_checkpoint',if(smoke)iter else m$burn+250L))),
    checkpoint_step=as.integer(get('checkpoint_step',1000L)),
    ess_min=if(iter==max(unlist(m$stages)))mixedgp_parse_ess_min(m$ess_min)else NULL,
    final_stage=iter==max(unlist(m$stages)),
    draws_per_chain=as.integer(get('draws_per_chain',if(smoke)5L else 250L)),protocol=m$trial_id,
    samplers=get('samplers',list(O0='HMC',O20='HMC',O50='ESS')))
}
adni_screen_manifest <- function(path,api) {
  m<-jsonlite::fromJSON(normalizePath(path),simplifyVector=FALSE)
  if(!m$schema%in%c('adni-screening-v2','adni-application-v1')||!length(m$cases))stop('Unsupported ADNI manifest')
  if(m$schema=='adni-screening-v2'&&!isTRUE(m$smoke)&&
    (m$burn!=500L||!identical(as.integer(unlist(m$folds)),1:5)))stop('Production screening requires all five folds and burn500')
  m$ess_min<-mixedgp_parse_ess_min(m$ess_min)
  if(!is.null(m$ess_min) && !is.null(m$chains) && m$chains<4L)stop('Scientific ESS requires four chains.')
  stages<-unlist(m$stages)
  if(!length(stages)||any(!is.finite(stages))||any(stages!=as.integer(stages))||any(diff(stages)<=0)||
    length(m$burn)!=1L||!is.finite(m$burn)||m$burn<0||m$burn!=as.integer(m$burn)||min(stages)<=m$burn)stop('Invalid stage budgets')
  if(!length(m$folds)||any(!unlist(m$folds)%in%1:5)||anyDuplicated(unlist(m$folds)))stop('Invalid folds')
  if(!is.null(m$chains)&&(length(m$chains)!=1L||m$chains<1||m$chains!=as.integer(m$chains)))stop('Invalid chain count')
  if(!all(unlist(m$cells)%in%c('O0','O20','O50')))stop('Invalid calibration cells')
  if(anyDuplicated(vapply(m$cases,`[[`,'','id')))stop('Repeated Plan identity')
  for(k in names(m$frozen_inputs))if(!identical(digest::digest(file=k,algo='sha256'),m$frozen_inputs[[k]]))stop('Frozen input changed: ',k)
  inputs<-list()
  for(z in m$cases){
    if(!grepl('^[A-Za-z0-9][A-Za-z0-9_-]*$',z$id))stop('Invalid dataset identity')
    adni_output_check(z$output,api$repository_dir)
    inp<-api$adni_plan_inputs(z$input)
    for(cell in inp$cells)if(cell$protocol$plan_id!=z$id||any(cell$test$calibration!=0L))stop('Dataset identity or test mask differs')
    inputs[[z$id]]<-inp
  }
  if(length(inputs)>1L)for(i in 2:length(inputs))for(key in names(inputs[[1]]$cells)){
    a<-inputs[[1]]$cells[[key]];b<-inputs[[i]]$cells[[key]]
    if(!identical(a$train$RID,b$train$RID)||!identical(a$test$RID,b$test$RID)||
       !identical(a$train[,c('Y','X1','X2','U1','U2','calibration')],b$train[,c('Y','X1','X2','U1','U2','calibration')])||
       !identical(a$available,b$available))stop('Paired people, anchors or masks differ')
  }
  list(manifest=m,inputs=inputs)
}
adni_screen_tasks <- function(m,group,stage='') {
  rows<-list();add<-function(plan,fold,cell,method,kind)rows[[length(rows)+1L]]<<-data.frame(plan=plan,fold=as.integer(fold),cell=cell,method=method,kind=kind,stage=stage)
  for(z in m$cases)for(f in unlist(m$folds)){
    if(group=='eivGP')for(cell in unlist(m$cells))add(z$id,f,cell,'EIV-GP','eiv')
    if(group=='competitors'){
      add(z$id,f,'OI','UC-GP+LVGP+EzGP','baseline_group')
      for(cell in intersect(unlist(m$cells),c('O20','O50')))for(method in c('UC-GP','LVGP','EzGP'))add(z$id,f,cell,method,'cc')
    }
  }
  if(!length(rows))return(data.frame())
  x<-do.call(rbind,rows);x$id<-paste(x$plan,x$stage,x$kind,x$fold,x$cell,gsub('[+]','_',x$method),sep='__');x
}
adni_screen_receipt <- function(out,task,phase)file.path(out,'receipts',paste0(phase,'_',task$id,'.rds'))
adni_screen_verified <- function(path,identity) {
  if(!file.exists(path))return(FALSE)
  r<-readRDS(path);if(!identical(r$identity,identity))stop('Receipt identity differs: ',path)
  if(length(r$files)&&(!all(file.exists(r$files))||!identical(unname(tools::md5sum(r$files)),r$md5)))stop('Completed output changed: ',path)
  TRUE
}
adni_screen_phase <- function(tasks,phase,cores,api,FUN,status_dir,chains=4L) {
  dir.create(status_dir,recursive=TRUE,showWarnings=FALSE)
  status<-tasks;status$state<-'PENDING';status$detail<-'';status$worker_pid<-NA_integer_;status$seconds<-NA_real_
  outer<-api$queue$ocean_stage_workers(cores,if(phase=='fit')'fit'else'analyze',nrow(tasks),chains)
  publish<-function()api$atomic_csv(status,file.path(status_dir,paste0(phase,'_status.csv')))
  publish()
  # Match the public Ocean executor: resident jobs share the full active-worker pool.
  api$mixedgp_budget_yield(function()api$queue$ocean_schedule(seq_len(nrow(tasks)),rep(1L,nrow(tasks)),outer,
    function(i,allocated)api$mixedgp_budget_worker(function(){
      logpath<-file.path(status_dir,paste0(phase,'_',make.names(tasks$id[i]),'.log'));con<-file(logpath,'at')
      old_message<-getConnection(sink.number(type='message'));sink(con);sink(con,type='message')
      on.exit({sink(old_message,type='message');sink();close(con)},add=TRUE);t<-proc.time()[3]
      err<-tryCatch({FUN(tasks[i,,drop=FALSE],cores);NULL},error=function(e)conditionMessage(e))
      if(!is.null(err))cat('FAILED:',err,'\n')
      list(state=if(is.null(err))'DONE'else'FAILED',detail=if(is.null(err))''else err,seconds=as.numeric(proc.time()[3]-t))
    }),
    function(i,n,pid,active){status$state[i]<<-'RUNNING';status$worker_pid[i]<<-pid;publish()},
    function(i,r,active){if(!is.list(r)||is.null(r$state))r<-list(state='FAILED',detail='Worker exited without a result',seconds=NA_real_)
      status$state[i]<<-r$state;status$detail[i]<<-r$detail;status$seconds[i]<<-r$seconds;publish()}))
  status
}
adni_screen_collect <- function(m,api) {
  pred<-conv<-u<-pros<-list()
  for(z in m$cases){
    for(mark in list.files(z$output,'^EVALUATION_VALIDATION[.]json$',recursive=TRUE,full.names=TRUE)){
      dir<-dirname(mark)
      for(pair in list(c('prediction_metrics.csv','pred'),c('convergence.csv','conv'),c('U_metrics.csv','u'),c('PROSPECTIVE_U_metrics.csv','pros'))){
        p<-file.path(dir,pair[1]);if(!file.exists(p))next
        val<-cbind(plan=z$id,read.csv(p));dest<-get(pair[2]);dest[[length(dest)+1L]]<-val;assign(pair[2],dest)
      }
    }
    files<-list.files(z$output,'^metrics[.]csv$',recursive=TRUE,full.names=TRUE)
    for(p in files){x<-read.csv(p);if('stage'%in%names(x)&&all(x$stage=='BASELINE'))pred[[length(pred)+1L]]<-cbind(plan=z$id,x)}
  }
  bind<-function(x)if(length(x))api$queue$ocean_bind(x)else data.frame()
  tables<-list(prediction=bind(pred),convergence=bind(conv),training_U=bind(u),prospective_U=bind(pros))
  dest<-file.path(m$output_root,'reports');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  for(k in names(tables))if(nrow(tables[[k]]))api$atomic_csv(tables[[k]],file.path(dest,paste0(k,'_fold_metrics.csv')))
  p<-tables$prediction;summary<-list()
  if(nrow(p))for(ix in split(seq_len(nrow(p)),interaction(p$plan,p$stage,p$cell,p$method,p$protocol,drop=TRUE))){
    x<-p[ix,,drop=FALSE];v<-x$RMSE_SD
    summary[[length(summary)+1L]]<-data.frame(x[1,c('plan','stage','cell','method','protocol')],completed_folds=length(unique(x$fold)),
       RMSE_SD_mean=mean(v),RMSE_SD_SD=if(length(v)>1)sd(v)else NA_real_)
  }
  five<-bind(summary);if(nrow(five))api$atomic_csv(five,file.path(dest,'fivefold_prediction.csv'))
  stability<-list()
  for(kind in names(tables)[names(tables)%in%c('prediction','training_U','prospective_U')]){
    tab<-tables[[kind]];if(!nrow(tab)||!'stage'%in%names(tab))next
    stage_names<-paste0('S',unlist(m$stages));if(length(stage_names)<2L)next
    x<-tab[tab$stage==stage_names[1],,drop=FALSE];y<-tab[tab$stage==tail(stage_names,1),,drop=FALSE]
    keys<-if(kind=='prediction')c('plan','fold','cell','method','protocol')else c('plan','fold','cell','coordinate')
    if(nrow(x)&&nrow(y)){
      z<-merge(x[,c(keys,'RMSE_SD'),drop=FALSE],y[,c(keys,'RMSE_SD'),drop=FALSE],by=keys,suffixes=c('_first','_last'))
      z$first_stage<-stage_names[1];z$last_stage<-tail(stage_names,1)
      z$target<-kind;z$delta_last_minus_first<-z$RMSE_SD_last-z$RMSE_SD_first;stability[[kind]]<-z
    }
  }
  stable<-bind(stability);if(nrow(stable))api$atomic_csv(stable,file.path(dest,'stage_stability.csv'))
  lines<-c('# ADNI screening: saved results',paste('Author: Xinye Chen. Trial',m$trial_id),
    'Partial fold counts and CORE failures remain explicit. Budget completion does not establish convergence.',
    if(nrow(five))capture.output(print(five,row.names=FALSE))else'No completed scores yet.')
  writeLines(lines,file.path(dest,'REPORT.md'))
  html<-c('<!doctype html><meta charset="utf-8"><title>ADNI screening</title>',
    '<h1>ADNI screening: saved results</h1><p>Author: Xinye Chen. Partial folds and CORE failures remain explicit; completing a budget does not establish convergence.</p>',
    api$queue$ocean_table_html(five,'Response prediction: fold means and sample SD'),
    api$queue$ocean_table_html(tables$convergence,'Full-chain CORE diagnostics'),
    api$queue$ocean_table_html(stable,'Paired stage stability; descriptive screening only'))
  writeLines(html,file.path(dest,'index.html'));invisible(tables)
}
adni_screen_run <- function(repo,p,group) {
  api<-adni_screen_api(repo);bundle<-adni_screen_manifest(p$manifest,api);m<-bundle$manifest;inputs<-bundle$inputs
  cases<-setNames(m$cases,vapply(m$cases,`[[`,'','id'));api$B<-dirname(normalizePath(p$manifest))
  if(p$plan){
    for(g in if(group=='report')c('eivGP','competitors')else group){x<-adni_screen_tasks(m,g);print(x,row.names=FALSE)}
    cat('Public shared CPU pool:',p$cores,'; fit resident limit:',api$queue$ocean_stage_workers(p$cores,'fit',100L,4L),
      '; full-chain reservations/backfill disabled.\n');return(invisible(bundle))
  }
  root<-m$output_root
  adni_output_check(root,repo)
  dir.create(root,recursive=TRUE,showWarnings=FALSE);lock<-file.path(root,'.adni-application.lock')
  if(!dir.create(lock,showWarnings=FALSE))stop('ADNI workflow already active: ',lock)
  writeLines(as.character(Sys.getpid()),file.path(lock,'pid'));on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  code<-file.path(repo,'eivGP/R/helpers.R')
  base_id<-list(schema='adni-upstream-screening-v2',trial=m$trial_id,modules=api$eivgp_module_identity(),
    application=unname(tools::md5sum(code)),burn=m$burn,stages=unlist(m$stages),smoke=isTRUE(m$smoke),ess_min=m$ess_min,
    chains=if(is.null(m$chains))4L else as.integer(m$chains),
    samplers=if(is.null(m$samplers))list(O0='HMC',O20='HMC',O50='ESS')else m$samplers,
    baseline_profile=if(is.null(m$baseline_profile))if(isTRUE(m$smoke))'pilot'else'paper'else m$baseline_profile,
    draws_per_chain=if(is.null(m$draws_per_chain))if(isTRUE(m$smoke))5L else 250L else m$draws_per_chain)
  identities<-list()
  for(z in m$cases){
    id<-c(base_id,list(plan=z$id,inputs=unname(tools::md5sum(inputs[[z$id]]$files))))
    identities[[z$id]]<-id;path<-file.path(z$output,'workflow_identity.rds')
    if(file.exists(path)){if(!identical(readRDS(path),id))stop('Saved workflow identity differs; use a new output root')}else if(group!='report'&&!p$reuse_only)api$atomic_rds(id,path)
  }
  required<-list()
  for(g in if(group=='report')c('eivGP','competitors')else group){
    for(st in if(g=='eivGP')paste0('S',unlist(m$stages))else'BASELINE'){
      tasks<-adni_screen_tasks(m,g,st)
      for(i in seq_len(nrow(tasks))){t<-tasks[i,];out<-cases[[t$plan]]$output
        required[[length(required)+1L]]<-list(path=adni_screen_receipt(out,t,if(g=='eivGP')'evaluate'else'fit'),identity=identities[[t$plan]])}
    }
  }
  if(group=='report'||p$reuse_only){
    if(!all(vapply(required,function(r)adni_screen_verified(r$path,r$identity),logical(1))))stop('Saved-only prerequisites incomplete; no fitting or continuation started')
    return(adni_screen_collect(m,api))
  }
  api$adni_plan_require();job<-file.path(root,'logs',paste0(group,'_',format(Sys.time(),'%Y%m%d_%H%M%S')))
  old<-options(eivGP.core_events=file.path(job,'core-events'));on.exit(options(old),add=TRUE)
  close<-api$mixedgp_budget_enter(p$cores);on.exit(close(),add=TRUE,after=FALSE);cores<-api$mixedgp_budget_cores()
  writeLines(paste('Effective shared worker budget',cores,'; public resident-fit scheduling.'),file.path(root,paste0(group,'_budget.txt')))
  stages<-if(group=='eivGP')unlist(m$stages)else 0L
  for(iter in stages){
    stage<-if(group=='eivGP')paste0('S',iter)else'BASELINE';tasks<-adni_screen_tasks(m,group,stage)
    work<-function(t,n,evaluate=FALSE){
      out<-cases[[t$plan]]$output;inp<-inputs[[t$plan]];id<-identities[[t$plan]];s<-adni_screen_settings(m,if(iter)iter else max(unlist(m$stages)),n)
      phase<-if(evaluate)'evaluate'else'fit';receipt<-adni_screen_receipt(out,t,phase)
      if(adni_screen_verified(receipt,id))return(invisible(NULL))
      if(t$kind=='eiv'){
        path<-api$adni_plan_fit_path(out,t)
        if(evaluate){
          if(!file.exists(path))stop('Stage checkpoint unavailable for evaluation')
          actual<-api$adni_plan_read_fit(path)$checkpoint$iteration
          if(is.null(s$ess_min) && actual!=iter || !is.null(s$ess_min) && actual<iter)stop('Stage checkpoint unavailable for evaluation')
          s$iter<-actual
          api$analyze_cell(t,inp,out,s,stage);folder<-file.path(out,stage,'analysis',sprintf('fold%02d',t$fold),t$cell)
        }else{
          if(file.exists(path)){
            fit<-api$adni_plan_read_fit(path)
            if(fit$checkpoint$iteration>iter && is.null(s$ess_min))stop('Cannot refit an earlier stage after continuation; restore its completed evaluation receipt')
            if(fit$checkpoint$burn!=s$burn||fit$control$sampler!=s$samplers[[t$cell]])stop('Checkpoint sampler/burn differs')
          }
          api$adni_plan_fit_task(t,inp,out,s,extend=FALSE);folder<-NULL
        }
      }else{
        methods<-if(t$kind=='baseline_group')c('UC-GP','LVGP','EzGP')else t$method
        failures<-character()
        for(method in methods){one<-t;one$method<-method;tryCatch(api$baseline_cell(one,inp,out,s,cc=t$kind=='cc'),error=function(e)failures<<-c(failures,paste(method,conditionMessage(e))))}
        if(length(failures))stop(paste(failures,collapse='; '))
        folder<-if(t$kind=='cc')file.path(out,'CC',t$cell,'baselines',sprintf('fold%02d',t$fold),t$method)else file.path(out,'baselines',sprintf('fold%02d',t$fold))
      }
      files<-if(is.null(folder))character()else sort(list.files(folder,recursive=TRUE,full.names=TRUE))
      api$atomic_rds(list(identity=id,task=t,files=files,md5=unname(tools::md5sum(files)),finished=as.character(Sys.time())),receipt)
    }
    nchains<-if(is.null(m$chains))4L else as.integer(m$chains)
    status<-adni_screen_phase(tasks,'fit',cores,api,function(t,n)work(t,n),file.path(job,stage),nchains)
    if(group=='eivGP'&&any(status$state=='FAILED'))stop('EIV stage incomplete; checkpoints retained. Inspect ',job)
    if(group=='eivGP'){
      status<-adni_screen_phase(tasks,'evaluate',cores,api,function(t,n)work(t,n,TRUE),file.path(job,stage),nchains)
      if(any(status$state=='FAILED'))stop('EIV evaluation incomplete; continuation not started. Inspect ',job)
    }else if(any(status$state=='FAILED'))warning('Comparator failures retained; saved results remain available.')
    adni_screen_collect(m,api)
  }
  for(z in m$cases)mixedgp_application_trace_report(z$output,file.path(z$output,'reports','scientific-traces'))
  invisible(adni_screen_collect(m,api))
}

# Author: Xinye Chen. Prepared-fold ADNI package API; legacy run_adni_analysis remains available.
adni_workflow_api <- function(repository_dir=NULL) {
  ns<-asNamespace('eivGP')
  if(!is.null(repository_dir)) {
    repository_dir<-normalizePath(repository_dir)
    manifest<-utils::read.csv(system.file('module-manifest.csv',package='eivGP'))
    native<-grepl('[.]cpp$|^Makevars',manifest$module)
    paths<-file.path(repository_dir,'eivGP',ifelse(native,'src','R'),manifest$module)
    if(!identical(unname(tools::md5sum(paths)),manifest$md5))
      stop('Installed eivGP differs from this checkout; rebuild and reinstall it first.')
  }
  e<-new.env(parent=ns)
  names<-setdiff(grep('^adni_',ls(ns,all.names=TRUE),value=TRUE),'adni_workflow_api')
  for(n in names) {
    value<-get(n,ns,inherits=FALSE)
    if(is.function(value))environment(value)<-e
    e[[n]]<-value
  }
  e$repository_dir<-repository_dir;e$api<-e;e$application<-e
  e$queue<-ocean_workflow_api()
  for(n in c('mixedgp_budget_yield','mixedgp_budget_worker','mixedgp_budget_enter',
             'mixedgp_budget_cores','eivgp_module_identity'))e[[n]]<-get(n,ns)
  e$atomic_json<-function(value,path) {
    dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
    p<-tempfile(tmpdir=dirname(path));on.exit(unlink(p))
    jsonlite::write_json(value,p,auto_unbox=TRUE,pretty=TRUE,null='null',na='null')
    stopifnot(file.rename(p,path))
  }
  e$atomic_csv<-ocean_atomic_csv;e$atomic_rds<-ocean_atomic_rds
  e$analyze_cell<-e$adni_analyze_cell;e$baseline_cell<-e$adni_baseline_cell
  e$adni_load_api<-e$adni_screen_api<-function(repo)e
  e
}
run_adni_application <- function(repository_dir,params=list(),group=c('report','eivGP','competitors')) {
  group<-match.arg(group)
  api<-adni_workflow_api(repository_dir)
  api$adni_method_run(repository_dir,params,group)
}

# Score one existing fit; this function never fits any model.
mixedgp_score_study1_fit <- function(fit_eiv, frozen, cache_spec, rep_id, n_calib,
                                    measurement_fit=NULL) {
  # Use the original fixed eight-state prediction streams even in serial replay.
  release_budget<-mixedgp_budget_enter(1L);on.exit(release_budget(),add=TRUE)
  dat<-frozen$data;train<-dat$train;test<-dat$test
  settings<-list(n_m_eval=cache_spec$mean_recovery$n_eval,
    n_m_draw=cache_spec$mean_recovery$n_draw,n_m_latent=cache_spec$mean_recovery$n_latent)
  n_pred_draw<-cache_spec$prediction$n_draw
  STUDY1_SCENARIO<-cache_spec$scenario
  STUDY1_HETEROGENEITY_ETA<-cache_spec$heterogeneity_eta
  STUDY1_EVALUATE_F<-cache_spec$task_flags$evaluate_f
  STUDY1_EVALUATE_U<-cache_spec$task_flags$evaluate_u
  STUDY1_QUICK<-identical(cache_spec$study1_config,"quick")
  set.seed(240000L+rep_id)
  mean_eval_idx<-sort(sample(seq_len(nrow(test)),min(settings$n_m_eval,nrow(test))))
  m_true<-m0_1d(test$x[mean_eval_idx],test$c[mean_eval_idx],tau=dat$tau_true,scenario=STUDY1_SCENARIO)
  metrics<-mean_metrics<-latent_imputation_metrics<-surface_recovery_metrics<-list()
  measurement_fits<-setNames(list(measurement_fit),as.character(n_calib))
      draw_ids <- seq_len(nrow(fit_eiv$mcmc$samples_u))
      if (length(draw_ids) > n_pred_draw) {
        set.seed(350000L + 1000L * rep_id + n_calib)
        draw_ids <- sample(draw_ids, n_pred_draw)
      }
      eiv_draws <- sample_eiv_test_y(
        x_test_raw = test$x,
        c_test = test$c,
        fit_obj = fit_eiv,
        draw_ids = draw_ids,
        n_per_draw = 1L
      )
      metrics[[paste0("EIV_", n_calib)]] <- summarize_predictive_samples_1d(
        eiv_draws, test$y, method = "EIV-GP", rep_id = rep_id,
        n_calib = n_calib, scenario = STUDY1_SCENARIO
      )
      m_draw_ids <- draw_ids
      if (length(m_draw_ids) > settings$n_m_draw) {
        set.seed(360000L + 1000L * rep_id + n_calib)
        m_draw_ids <- sample(m_draw_ids, settings$n_m_draw)
      }
      m_draws <- sample_eiv_m_given_xc_1d(
        x_star_raw = test$x[mean_eval_idx],
        c_star = test$c[mean_eval_idx],
        fit_obj = fit_eiv,
        draw_ids = m_draw_ids,
        n_latent = settings$n_m_latent,
        include_process_uncertainty = TRUE,
        joint = FALSE,
        return_components = TRUE,
        seed = 370000L + 1000L * rep_id + n_calib
      )
      mean_metrics[[paste0("EIV_", n_calib)]] <-
        summarize_mean_recovery_1d(
          m_draws,
          m_true = m_true,
          method = "EIV-GP",
          rep_id = rep_id,
          n_calib = n_calib,
          scenario = STUDY1_SCENARIO,
          valid_function_draws = TRUE
        )
      if (isTRUE(STUDY1_EVALUATE_U)) {
        latent_imputation_metrics[[as.character(n_calib)]] <-
          study1_latent_imputation_metrics(
            eiv_fit = fit_eiv,
            measurement_fit = measurement_fits[[as.character(n_calib)]],
            U_true = train$u,
            rep_id = rep_id,
            n_calib = n_calib,
            scenario = STUDY1_SCENARIO
          )
      }
      if (isTRUE(STUDY1_EVALUATE_F) && n_calib > 0L) {
        surface_recovery_metrics[[as.character(n_calib)]] <-
          study1_eiv_surface_recovery_metrics(
            fit = fit_eiv,
            scenario = STUDY1_SCENARIO,
            rep_id = rep_id,
            n_calib = n_calib,
            tau_true = dat$tau_true,
            heterogeneity_eta = STUDY1_HETEROGENEITY_ETA,
            grid_n_x = if (STUDY1_QUICK) 15L else 31L,
            grid_n_u = if (STUDY1_QUICK) 21L else 51L,
            max_draw = if (STUDY1_QUICK) 40L else 200L,
            seed = 380000L + 1000L * rep_id + n_calib
          )
      }
  list(metrics=bind_rows(metrics),mean_recovery=bind_rows(mean_metrics),
    latent_imputation=bind_rows(latent_imputation_metrics),surface_recovery=bind_rows(surface_recovery_metrics))
}
# Refresh a single saved replication after exact continuation; never fit baselines.
mixedgp_refresh_study1_continuation <- function(run_dir,cell_id,rep_id,n_calib,extended,data_root,ess_min=100) {
  if(!isTRUE(extended$success) || !mixedgp_scientific_ess_pass(extended$diagnostics,ess_min))
    stop("Continuation is not scientifically ESS-complete.")
  cell_dir<-file.path(run_dir,"cells",cell_id)
  res<-file.path(cell_dir,"results","study1_publication")
  fit_path<-file.path(res,sprintf("mcmc_failure_rep%03d_cal%03d.rds",rep_id,n_calib))
  if(!file.exists(fit_path))fit_path<-file.path(res,sprintf("fit_rep%03d_cal%03d.rds",rep_id,n_calib))
  saved<-read_fit(fit_path);old<-saved$fit;fit<-extended$fit
  if(!identical(old$data,fit$data) || !identical(old$control,fit$control) || !identical(old$priors,fit$priors))
    stop("Continuation changed the data, controls or priors.")
  nold<-old$checkpoint$iteration-old$checkpoint$burn
  for(name in names(old$mcmc$samples_by_chain))for(i in seq_along(old$mcmc$samples_by_chain[[name]])) {
    a<-old$mcmc$samples_by_chain[[name]][[i]];b<-fit$mcmc$samples_by_chain[[name]][[i]]
    if(is.null(a))next
    prefix<-if(is.null(dim(b)))b[seq_len(nold)] else b[seq_len(nold),,drop=FALSE]
    if(!identical(a,prefix))stop("Continuation did not preserve the original draws: ",name)
  }
  spec<-saved$cache_spec
  data_path<-file.path(data_root,cell_id,mixedgp_dataset_filename("study1",rep_id,
    spec$scenario,spec$n_train,spec$n_test,spec$m))
  manifest<-utils::read.csv(file.path(run_dir,"config","input_manifest.csv"),stringsAsFactors=FALSE)
  matchrow<-manifest[manifest$cell_id==cell_id & basename(manifest$file)==basename(data_path),,drop=FALSE]
  if(nrow(matchrow)!=1L || !identical(tolower(as.character(matchrow$md5)),tolower(unname(tools::md5sum(data_path)))))
    stop("Frozen data identity mismatch.")
  frozen<-readRDS(data_path)
  panel<-saved$diagnostic_panel_rows;t<-frozen$data$test
  expected_U<-if(spec$heterogeneity_eta==1 && isTRUE(spec$task_flags$evaluate_f))t$u[panel]else NULL
  expected_targets<-saved$diagnostics$parameter[saved$diagnostics$diagnostic_role=="target"]
  if(!setequal(extended$diagnostics$parameter,expected_targets) ||
     anyDuplicated(extended$diagnostics$parameter) ||
     !identical(extended$settings$U,expected_U) ||
     !identical(extended$settings$X,matrix(t$x[panel],ncol=1L)) ||
     !identical(extended$settings$C,t$c[panel]) ||
     extended$settings$n_latent!=spec$diagnostic_panel$n_latent ||
     extended$settings$seed!=380000L+1000L*rep_id+n_calib)
    stop("Continuation target panel or integration settings changed.")
  scored<-mixedgp_score_study1_fit(fit,frozen,spec,rep_id,n_calib)
  raw<-mixedgp_study1_raw_diagnostics(fit,rhat_limit=spec$mcmc$rhat_limit,
    bulk_ess_limit=spec$mcmc$ess_limit,tail_ess_limit=spec$mcmc$ess_limit)
  targets<-extended$diagnostics
  raw$diagnostic_role<-"raw";targets$diagnostic_role<-"target"
  tab<-bind_rows(raw,targets);tab$rep<-rep_id;tab$n_calib<-n_calib
  row<-fit$diagnostics$summary;row$rep<-rep_id;row$n_calib<-n_calib
  row$gate_max_rhat<-max(tab$rhat);row$gate_min_ess<-min(tab$ess_bulk,tab$ess_tail)
  row$raw_gate_pass<-mixedgp_diagnostic_table_pass(raw)
  row$target_gate_pass<-mixedgp_diagnostic_table_pass(targets)
  row$scientific_ess_pass<-TRUE;row$gate_pass<-TRUE
  row$diagnostic_warning<-FALSE;row$diagnostic_advice<-""
  replace<-function(x,y,method=FALSE) {
    if(!is.data.frame(y)||!nrow(y))return(x)
    hit<-x$rep==rep_id & x$n_calib==n_calib;hit[is.na(hit)]<-FALSE
    if(method)hit<-hit & x$method=="EIV-GP"
    out<-bind_rows(x[!hit,,drop=FALSE],y)
    out[order(out$rep,out$n_calib,na.last=TRUE),,drop=FALSE]
  }
  rep_path<-list.files(file.path(res,"competitor_replications"),
    pattern=sprintf("^study1_rep%03d_.*[.]rds$",rep_id),full.names=TRUE)
  if(length(rep_path)!=1L)stop("Expected one saved replication bundle.")
  bundle<-read_fit(rep_path)
  for(name in names(scored))bundle[[name]]<-replace(bundle[[name]],scored[[name]],TRUE)
  bundle$mcmc_diagnostics<-replace(bundle$mcmc_diagnostics,row)
  bundle$mcmc_parameter_diagnostics<-replace(bundle$mcmc_parameter_diagnostics,tab)
  provenance<-list(updated_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
    source_fit_md5=unname(tools::md5sum(fit_path)),initial_iteration=old$checkpoint$iteration,
    final_iteration=fit$checkpoint$iteration,ess_min=ess_min,history=extended$history,
    modules=eivgp_module_identity(),data_md5=unname(tools::md5sum(data_path)))
  bundle$metadata$continuations[[as.character(n_calib)]]<-provenance
  report_path<-file.path(cell_dir,"report_inputs.rds");report<-readRDS(report_path)
  mapping<-c(metrics="mc_results",mean_recovery="mean_recovery",latent_imputation="latent_imputation",surface_recovery="surface_recovery")
  for(name in names(mapping))report$context[[mapping[[name]]]]<-replace(report$context[[mapping[[name]]]],scored[[name]],TRUE)
  report$context$mcmc_diagnostics<-replace(report$context$mcmc_diagnostics,row)
  report$context$mcmc_parameter_diagnostics<-replace(report$context$mcmc_parameter_diagnostics,tab)
  # Persist scores/report inputs before publishing the replacement fit. If interrupted,
  # the same completed extension can safely replay this deterministic update.
  mixedgp_save_replication(bundle,rep_path)
  mixedgp_atomic_save_rds(report,report_path)
  saved$fit<-fit;saved$diagnostics<-tab;saved$continuation<-provenance
  mixedgp_checkpoint_write(mixedgp_pack_artifact(saved),fit_path)
  n<-length(extended$target_series[[1]][[1]]);keep<-unique(c(seq.int(1L,n,by=max(1L,ceiling(n/750L))),n))
  traces<-lapply(extended$target_series,function(chains)lapply(chains,function(x)x[keep]))
  traces<-mixedgp_mark_diagnostic_window(traces,length(keep)==n)
  mixedgp_save_artifact(list(schema_version=1L,series=traces,iteration=fit$checkpoint$burn+keep,
    retained_index=keep,panel_rows=panel,test_panel=t[panel,],diagnostics=targets,
    source="full retained history after scientific-target ESS continuation"),
    file.path(res,sprintf("target_traces_rep%03d_cal%03d.rds",rep_id,n_calib)))
  invisible(provenance)
}

#' Continue only Study I fits below the scientific-target ESS requirement
#' @param path Saved Study I run directory.
#' @param data_root Frozen Study I data root.
#' @param continuation_dir Directory for resumable extension states.
#' @param ess_min Strict bulk and tail ESS lower bound.
#' @param chunk_iter Additional transitions per chain between checks.
#' @param n_cores Number of chain workers; fits are processed sequentially.
#' @param cleanup Remove files under this run's checkpoints after all required
#'   fits pass and refreshed result exports have been verified. Final fit states remain.
#' @return Updated run path and an audit of continued fits.
#' @export
continue_study1_results <- function(path,data_root,continuation_dir=file.path(path,"continuation"),
    ess_min=100,chunk_iter=5000L,n_cores=1L,cleanup=TRUE) {
  path<-mixedgp_study1_report_run_dir(path)
  release_budget<-mixedgp_budget_enter(n_cores);on.exit(release_budget(),add=TRUE)
  checkpoint_options<-options(eivGP.checkpoint_root=file.path(path,"checkpoints"),
    eivGP.checkpoint_seconds=getOption("eivGP.checkpoint_seconds",120))
  on.exit(options(checkpoint_options),add=TRUE)
  cfg<-readRDS(file.path(path,"config","resolved_config.rds"))
  dir.create(continuation_dir,recursive=TRUE,showWarnings=FALSE)
  audit<-list()
  for(cell in cfg$cells) {
    report<-readRDS(file.path(path,"cells",cell$id,"report_inputs.rds"))
    all<-report$context$mcmc_parameter_diagnostics
    targets<-all[all$diagnostic_role=="target",,drop=FALSE]
    keys<-unique(targets[,c("rep","n_calib")])
    if(nrow(keys)!=cell$n_rep*length(cell$calibration_grid))stop("Missing target diagnostics in ",cell$id)
    for(j in seq_len(nrow(keys))) {
      k<-keys[j,];tab<-targets[targets$rep==k$rep & targets$n_calib==k$n_calib,,drop=FALSE]
      if(mixedgp_scientific_ess_pass(tab,ess_min))next
      res<-file.path(path,"cells",cell$id,"results","study1_publication")
      p<-file.path(res,sprintf("mcmc_failure_rep%03d_cal%03d.rds",k$rep,k$n_calib))
      if(!file.exists(p))p<-file.path(res,sprintf("fit_rep%03d_cal%03d.rds",k$rep,k$n_calib))
      saved<-read_fit(p);spec<-saved$cache_spec
      dest<-file.path(continuation_dir,cell$id,sprintf("rep%03d-cal%03d.rds",k$rep,k$n_calib))
      extended<-if(file.exists(dest))read_fit(dest)else NULL
      fit<-if(is.null(extended))saved$fit else extended$fit
      data_path<-file.path(data_root,cell$id,mixedgp_dataset_filename("study1",k$rep,spec$scenario,spec$n_train,spec$n_test,spec$m))
      t<-readRDS(data_path)$data$test;panel<-saved$diagnostic_panel_rows
      if(is.null(extended)||!isTRUE(extended$success)||!mixedgp_scientific_ess_pass(extended$diagnostics,ess_min))
        extended<-continue_eivgp_until_ess(fit,matrix(t$x[panel],ncol=1L),t$c[panel],
          if(spec$heterogeneity_eta==1 && isTRUE(spec$task_flags$evaluate_f))t$u[panel]else NULL,
          ess_min=ess_min,chunk_iter=chunk_iter,n_latent=spec$diagnostic_panel$n_latent,
          seed=380000L+1000L*k$rep+k$n_calib,parallel=n_cores>1,n_cores=n_cores,checkpoint_file=dest)
      audit[[paste(cell$id,k$rep,k$n_calib,sep="/")]]<-mixedgp_refresh_study1_continuation(path,cell$id,k$rep,k$n_calib,extended,data_root,ess_min)
    }
  }
  # Rebuild all summaries from the updated saved rows; no model is refitted.
  results<-list()
  mapping<-c(predictive_metrics="mc_results",mean_recovery="mean_recovery",latent_imputation="latent_imputation",
    surface_recovery="surface_recovery",ablation_predictive_metrics="ablation_results",ablation_surface_recovery="ablation_surface_recovery",
    mcmc_diagnostics="mcmc_diagnostics",competitor_status="competitor_status",ablation_status="ablation_status",sampler_control_manifest="sampler_control_manifest")
  for(cell in cfg$cells) {
    rp<-file.path(path,"cells",cell$id,"report_inputs.rds");report<-readRDS(rp);ctx<-report$context
    target<-ctx$mcmc_parameter_diagnostics;target<-target[target$diagnostic_role=="target",,drop=FALSE]
    d<-ctx$mcmc_diagnostics
    d$scientific_ess_pass<-vapply(seq_len(nrow(d)),function(i)
      mixedgp_scientific_ess_pass(target[target$rep==d$rep[i]&target$n_calib==d$n_calib[i],,drop=FALSE],ess_min),logical(1))
    if(!all(d$scientific_ess_pass)||nrow(d)!=cell$n_rep*length(cell$calibration_grid))stop("Incomplete scientific-target results.")
    d$gate_pass<-d$scientific_ess_pass;d$diagnostic_warning<-!d$scientific_ess_pass;d$diagnostic_advice<-""
    ctx$mcmc_diagnostics<-d;report$context<-ctx;mixedgp_atomic_save_rds(report,rp)
    tables<-file.path(path,"cells",cell$id,"tables")
    mixedgp_write_diagnostic_tables(d,ctx$mcmc_parameter_diagnostics,tables,"study1",overwrite=TRUE)
    results[[cell$id]]<-list(cell=cell,outputs=setNames(lapply(mapping,function(n)ctx[[n]]),names(mapping)))
  }
  mixedgp_aggregate_results(results,cfg,path)
  gates<-mixedgp_bind_rows_base(lapply(results,mixedgp_result_gates,config=cfg))
  mixedgp_atomic_write_csv(gates,file.path(path,"config","diagnostic_gates.csv"))
  summary_path<-file.path(path,"run_summary.rds");summary<-readRDS(summary_path);summary$gates<-gates
  summary$continuation_audit<-c(summary$continuation_audit,audit);mixedgp_atomic_save_rds(summary,summary_path)
  status_path<-file.path(path,"config","cell_status.csv");status<-utils::read.csv(status_path,stringsAsFactors=FALSE)
  for(cell in cfg$cells)if(all(gates$pass[gates$cell_id==cell$id])) {
    status$status[status$cell_id==cell$id]<-"success";status$message[status$cell_id==cell$id]<-""
  }
  mixedgp_atomic_write_csv(status,status_path)
  reporting_status<-mixedgp_report_run(path,"report",code_dir=system.file(package="eivGP"))
  if(any(reporting_status$status!="success"))stop("Reporting incomplete; checkpoints retained.")
  if(cleanup)clean_study1_checkpoints(path,ess_min=ess_min)
  invisible(list(run_dir=path,audit=audit))
}

#' Remove redundant Study I checkpoint files after successful inference
#' @param path Saved Study I run directory.
#' @param ess_min Strict scientific-target ESS lower bound.
#' @return File removal audit, stored outside checkpoints.
#' @details Requires complete replication scores and full-window target diagnostics
#'   for every expected fit, and checks final continuation states when present.
#'   Active locks stop cleanup. Only the run's checkpoints subtree is removed;
#'   final fit artifacts and continuation states outside that subtree are retained.
#' @export
clean_study1_checkpoints <- function(path,ess_min=100) {
  cfg<-readRDS(file.path(path,"config","resolved_config.rds"))
  if(!identical(cfg$study,"study1"))stop("Expected a Study I run.")
  for(cell in cfg$cells) {
    ctx<-readRDS(file.path(path,"cells",cell$id,"report_inputs.rds"))$context
    d<-ctx$mcmc_parameter_diagnostics;d<-d[d$diagnostic_role=="target",,drop=FALSE]
    for(rep in seq_len(cell$n_rep))for(cal in cell$calibration_grid) {
      tab<-d[d$rep==rep & d$n_calib==cal,,drop=FALSE]
      fam<-c("m_conditional_mean","m_conditional_variance","predictive_second_moment")
      if(cal>0L)fam<-c(fam,"prospective_U_mean","prospective_U_second_moment")
      if(cal>0L && ctx$STUDY1_HETEROGENEITY_ETA==1)fam<-c(fam,"f_conditional_mean","f_conditional_variance")
      expected<-unlist(lapply(fam,function(f)paste0(f,"[",seq_len(5L),"]")),use.names=FALSE)
      if(!setequal(tab$parameter,expected) || anyDuplicated(tab$parameter) ||
         !mixedgp_scientific_ess_pass(tab,ess_min))stop("Checkpoint cleanup requires complete successful scientific-target diagnostics for every fit.")
      rows<-ctx$mc_results[ctx$mc_results$rep==rep & ctx$mc_results$method=="EIV-GP" & !is.na(ctx$mc_results$n_calib) & ctx$mc_results$n_calib==cal,,drop=FALSE]
      if(nrow(rows)!=1L || !is.finite(rows$RMSE))stop("Missing final prediction scores; checkpoints retained.")
    }
  }
  root<-file.path(path,"checkpoints")
  files<-list.files(root,recursive=TRUE,full.names=TRUE,all.files=TRUE,no..=TRUE)
  if(any(basename(list.dirs(root,recursive=TRUE))=="LOCK"))stop("Active checkpoint lock; cleanup refused.")
  if(length(files)) {
    identities<-files[basename(files)=="identity.rds"]
    if(!length(identities))stop("Unidentified checkpoint files; cleanup refused.")
    checked<-new.env(parent=emptyenv())
    for(identity_path in identities) {
      model<-readRDS(identity_path)$model
      seed<-model$seed
      if(length(seed)!=1L || !is.finite(seed))stop("Unidentified checkpoint model; cleanup refused.")
      rep<-as.integer((seed-300000L)%/%1000L);cal<-as.integer((seed-300000L)%%1000L)
      verified<-FALSE
      for(cell in cfg$cells) {
        if(!rep %in% seq_len(cell$n_rep)||!cal %in% cell$calibration_grid)next
        folder<-file.path(path,"cells",cell$id,"results","study1_publication")
        p<-file.path(folder,sprintf("fit_rep%03d_cal%03d.rds",rep,cal))
        if(!file.exists(p))p<-file.path(folder,sprintf("mcmc_failure_rep%03d_cal%03d.rds",rep,cal))
        if(!file.exists(p))next
        if(!exists(p,checked,inherits=FALSE)) {
          fit<-read_fit(p)$fit;cp<-fit$checkpoint
          valid<-identical(cp$version,2L) && !is.null(cp$states) && length(cp$states)>=4L &&
            length(fit$mcmc$samples_by_chain$u)==length(cp$states) &&
            length(fit$mcmc$samples_by_chain$logtheta)==length(cp$states) &&
            all(vapply(cp$states,function(st)!is.null(st$state)&&!is.null(st$rng),logical(1))) &&
            all(vapply(fit$mcmc$samples_by_chain$u,nrow,integer(1))==cp$iteration-cp$burn) &&
            all(vapply(fit$mcmc$samples_by_chain$logtheta,nrow,integer(1))==cp$iteration-cp$burn)
          assign(p,list(data=fit$data,iteration=cp$iteration,valid=valid),checked)
        }
        saved<-get(p,checked,inherits=FALSE);d<-saved$data
        if(isTRUE(saved$valid) && saved$iteration>=model$n_iter &&
           identical(d$x,model$X) && identical(d$y,model$y) &&
           identical(as.numeric(d$u_obs),as.numeric(model$U_obs)) &&
           identical(as.integer(d$c_ord),as.integer(model$C)))verified<-TRUE
      }
      if(!verified)stop("No matching final fit with complete terminal states; checkpoints retained: ",identity_path)
    }
  }
  audit<-data.frame(file=files,bytes=file.info(files)$size,md5=unname(tools::md5sum(files)),removed=rep(FALSE,length(files)))
  audit_path<-file.path(path,"config","checkpoint_cleanup.csv")
  mixedgp_atomic_write_csv(audit,audit_path)
  if(length(files)) {
    for(p in files)if(unlink(p)!=0L || file.exists(p))stop("Could not remove checkpoint: ",p)
    audit$removed<-TRUE;mixedgp_atomic_write_csv(audit,audit_path)
  }
  invisible(audit)
}


mixedgp_parse_ess_min <- function(x=NULL) {
  if(is.null(x)||identical(x,"NULL"))return(NULL)
  if(is.character(x))x<-suppressWarnings(as.numeric(x))
  mixedgp_validate_ess_request(x)
  x
}

mixedgp_application_ess <- function(fit,ess_min,panel,checkpoint_file,cores=1L) {
  ess_min<-mixedgp_parse_ess_min(ess_min)
  if(!is.null(ess_min)) {
    if(file.exists(checkpoint_file)) {
      old<-mixedgp_read_artifact(checkpoint_file)
      if(!identical(old$settings$X,panel$X)||!identical(old$settings$C,panel$C)||!identical(old$settings$U,panel$U)||
         !identical(old$fit$data,fit$data))stop("Saved scientific continuation identity differs.")
      if(old$fit$checkpoint$iteration>=fit$checkpoint$iteration)fit<-old$fit
    }
    return(mixedgp_finish_ess(fit,ess_min,panel,cores>1L,cores,checkpoint_file))
  }
  series<-mixedgp_continuation_target_series(fit,panel$X,panel$C,panel$U,panel$n_latent,panel$seed,cores>1L,cores)
  fit$scientific_ess<-list(target_series=series,diagnostics=mixedgp_summarize_diagnostic_series(series),
    settings=c(panel,list(ess_min=NULL)),success=NA,history=data.frame())
  fit
}

mixedgp_application_trace_report <- function(root,dest) {
  index<-list.files(root,'^scientific-trace-.*[.]csv$',recursive=TRUE,full.names=TRUE)
  if(!length(index))return(character())
  mixedgp_save_scientific_trace_plots(index,dest)
}


# CUT functionals must hold the outer measurement state fixed; the incidental
# A/tau draws from the all-inputs-observed conditional GP fit are not CUT draws.
mixedgp_cut_conditional_series <- function(fit,measurement,X,C,U=NULL,n_latent=64L,seed=1L) {
  display<-fit
  by<-display$mcmc$samples_by_chain
  tau<-flatten_tau(measurement$tau)
  for(ch in seq_along(by$A)) {
    n<-dim(by$A[[ch]])[1L]
    by$A[[ch]]<-array(rep(as.vector(measurement$A),each=n),dim(by$A[[ch]]))
    by$tau[[ch]]<-matrix(rep(tau,each=n),nrow=n)
  }
  display$mcmc$samples_by_chain<-by
  mixedgp_study2_target_series(display,X,C,U,max_draws_per_chain=Inf,n_latent=n_latent,seed=seed)
}

mixedgp_cut_outer_trace <- function(pred,states,test,cfg,folder,label) {
  rows<-mixedgp_study2_diagnostic_rows(test$X,test$C,5L)
  chain<-vapply(states,function(s)as.integer(s$chain),integer(1))
  iteration<-vapply(states,function(s)as.integer(s$iteration),integer(1))
  iterations<-sort(unique(iteration))
  groups<-split(seq_along(states),paste(chain,iteration,sep="/"))
  series<-list()
  for(i in seq_along(rows)) {
    r<-rows[i]
    for(family in c("m_conditional_mean","f_conditional_mean","predictive_second_moment")) {
      values<-if(family=="m_conditional_mean")pred$m_mean[r,]else if(family=="f_conditional_mean")pred$f_mean[r,]else {
        v<-pred$y_var[r,]+pred$y_mean[r,]^2
        colMeans(matrix(v,nrow=cfg$gp_response_latent))
      }
      stopifnot(length(values)==length(states))
      series[[paste0(family,"[",i,"]")]]<-lapply(1:4,function(ch)vapply(iterations,function(it) {
        ids<-groups[[paste(ch,it,sep="/")]]
        if(!length(ids))stop("Incomplete CUT outer chain history.")
        mean(values[ids])
      },numeric(1)))
    }
  }
  series<-mixedgp_mark_diagnostic_window(series,TRUE)
  tab<-mixedgp_summarize_diagnostic_series(series)
  # This iteration vector is the complete retained outer integration sequence,
  # not every measurement transition. Keep that distinction in the plot caption.
  fake<-list(checkpoint=list(burn=0L,iteration=length(iterations)),
    scientific_ess=list(target_series=series,diagnostics=tab,scope="outer_integration",success=NA))
  prefix<-file.path(folder,"scientific-trace-cut-outer")
  mixedgp_save_scientific_trace(fake,prefix,label)
  obj<-readRDS(paste0(prefix,".rds"));obj$iteration<-iterations
  mixedgp_atomic_save_rds(obj,paste0(prefix,".rds"))
  utils::write.csv(tab,file.path(folder,"cut_outer_target_diagnostics.csv"),row.names=FALSE)
}

mixedgp_cut_trace_report <- function(root,dest=file.path(root,"cut-traces"),formats=c("png","pdf")) {
  paths<-character()
  for(kind in c("conditional","outer")) {
    index<-list.files(root,paste0("^scientific-trace-cut-",kind,"[.]csv$"),recursive=TRUE,full.names=TRUE)
    if(!length(index))next
    made<-mixedgp_save_scientific_trace_plots(index,file.path(dest,kind),formats)
    for(p in made){target<-file.path(dest,paste0(kind,"_",basename(p)));file.copy(p,target,overwrite=TRUE);paths<-c(paths,target)}
  }
  paths
}


mixedgp_application_numerics <- function(nums) {
  if(nrow(unique(nums[,c("burn","chains"),drop=FALSE]))!=1L)
    stop("Selected report stage mixes warmup or chain counts.")
  if(nrow(nums)==1L)return(nums)
  out<-nums[1,,drop=FALSE]
  out$iterations<-min(nums$iterations);out$iterations_max<-max(nums$iterations)
  out
}
