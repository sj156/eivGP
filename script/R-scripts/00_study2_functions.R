############################################################
## 00_study2_functions.R
##
## Exact fully Bayesian ordinal-probit EIV-GP code
## plus revised Study II synthetic data utilities.
##
## Revised Study II design:
##   n_train = 120
##   U ~ N_2(0, I)
##   S = A U + zeta, zeta ~ N_4(0, I)
##   four four-level ordinal proxies
##   unknown A and cutpoints in fitted EIV-GP
##   exact GP-marginalized posterior sampler
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))

load_or_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

load_or_install("ggplot2")
load_or_install("patchwork")
load_or_install("dplyr")
load_or_install("tidyr")
load_or_install("knitr")

theme_set(theme_bw(base_size = 12))

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

STUDY2_DESIGN_TAG <- "exact_v3_n120_Iomega_balanced"

############################################################
## Configurations
############################################################

study2_config_settings <- function(config = c("quick", "balanced", "thorough")) {
  config <- match.arg(config)
  
  if (config == "quick") {
    return(list(
      config = config,
      n_test = 150L,
      n_rep = 2L,
      rep_n_iter = 1000L,
      rep_burn = 400L,
      rep_thin = 2L,
      rep_n_chains = 2L,
      mc_n_iter = 900L,
      mc_burn = 400L,
      mc_thin = 2L,
      mc_n_chains = 2L,
      n_pred_draw = 120L,
      n_density_draw = 200L,
      n_new_latent_gibbs = 12L,
      n_oracle_pool = 50000L,
      n_starts_learned = 2L,
      preset = "fast"
    ))
  }
  
  if (config == "balanced") {
    return(list(
      config = config,
      n_test = 400L,
      n_rep = 30L,
      rep_n_iter = 6000L,
      rep_burn = 2000L,
      rep_thin = 1L,
      rep_n_chains = 8L,
      mc_n_iter = 6000L,
      mc_burn = 2000L,
      mc_thin = 1L,
      mc_n_chains = 8L,
      n_pred_draw = 700L,
      n_density_draw = 1000L,
      n_new_latent_gibbs = 20L,
      n_oracle_pool = 150000L,
      n_starts_learned = 4L,
      preset = "balanced"
    ))
  }
  
  list(
    config = config,
    n_test = 500L,
    n_rep = 100L,
    rep_n_iter = 10000L,
    rep_burn = 3000L,
    rep_thin = 1L,
    rep_n_chains = 12L,
    mc_n_iter = 10000L,
    mc_burn = 3000L,
    mc_thin = 1L,
    mc_n_chains = 12L,
    n_pred_draw = 800L,
    n_density_draw = 1200L,
    n_new_latent_gibbs = 30L,
    n_oracle_pool = 300000L,
    n_starts_learned = 6L,
    preset = "thorough"
  )
}

make_default_control_ordprobit <- function(n,
                                           n_mis,
                                           preset = c("fast", "balanced", "thorough"),
                                           d = 2L) {
  preset <- match.arg(preset)
  
  if (preset == "fast") {
    theta_update_every <- 6L
    n_u_blocks_per_iter <- 2L
    u_block_size <- min(n_mis, max(8L, ceiling(1.25 * sqrt(max(n_mis, 1L)))))
    global_u_every <- 40L
  }
  
  if (preset == "balanced") {
    theta_update_every <- 3L
    n_u_blocks_per_iter <- 3L
    u_block_size <- min(n_mis, max(10L, ceiling(1.75 * sqrt(max(n_mis, 1L)))))
    global_u_every <- 20L
  }
  
  if (preset == "thorough") {
    theta_update_every <- 2L
    n_u_blocks_per_iter <- 4L
    u_block_size <- min(n_mis, max(12L, ceiling(2.25 * sqrt(max(n_mis, 1L)))))
    global_u_every <- 10L
  }
  
  list(
    preset = preset,
    
    theta_update_every = theta_update_every,
    n_u_blocks_per_iter = n_u_blocks_per_iter,
    u_block_size = u_block_size,
    global_u_every = global_u_every,
    
    score_update_every = 1L,
    A_update_every = 1L,
    tau_update_every = 1L,
    
    theta_slice_width_init = rep(0.8, 1L + d + 10L),
    adapt_theta_width = TRUE,
    adapt_every = 100L,
    adapt_window = 500L,
    theta_width_min = 0.15,
    theta_width_max = 2.50,
    
    tau_bound = 6,
    s_A = 2.5,
    
    max_ess_try = 300L
  )
}

############################################################
## General utilities
############################################################

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

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

safe_se <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(0)
  stats::sd(x) / sqrt(length(x))
}

format_mean_se <- function(mean, se, digits = 3) {
  paste0(
    sprintf(paste0("%.", digits, "f"), mean),
    " (",
    sprintf(paste0("%.", digits, "f"), se),
    ")"
  )
}

pairwise_sqdist <- function(a, b = NULL) {
  a <- as.matrix(a)
  if (is.null(b)) b <- a
  b <- as.matrix(b)
  aa <- rowSums(a^2)
  bb <- rowSums(b^2)
  out <- outer(aa, bb, "+") - 2 * tcrossprod(a, b)
  pmax(out, 0)
}

safe_chol <- function(A, jitter = 1e-8) {
  n <- nrow(A)
  for (k in 0:8) {
    jj <- jitter * 10^k
    ans <- try(chol(A + jj * diag(n)), silent = TRUE)
    if (!inherits(ans, "try-error")) return(ans)
  }
  stop("Cholesky decomposition failed.")
}

solve_chol <- function(U, b) {
  backsolve(U, forwardsolve(t(U), b))
}

rmvnorm_chol <- function(n, mean, Sigma) {
  mean <- as.numeric(mean)
  d <- length(mean)
  U <- safe_chol(Sigma)
  Z <- matrix(rnorm(n * d), n, d)
  sweep(Z %*% U, 2, mean, "+")
}

rmvnorm_rows_common <- function(mean_mat, Sigma) {
  mean_mat <- as.matrix(mean_mat)
  n <- nrow(mean_mat)
  d <- ncol(mean_mat)
  U <- safe_chol(Sigma)
  Z <- matrix(rnorm(n * d), n, d)
  mean_mat + Z %*% U
}

rtruncnorm_vec <- function(mean, sd, lower, upper) {
  n <- max(length(mean), length(lower), length(upper))
  mean <- rep(mean, length.out = n)
  sd <- rep(sd, length.out = n)
  lower <- rep(lower, length.out = n)
  upper <- rep(upper, length.out = n)
  
  pl <- pnorm((lower - mean) / sd)
  pu <- pnorm((upper - mean) / sd)
  
  pl <- pmax(pl, 0)
  pu <- pmin(pu, 1)
  
  width <- pmax(pu - pl, .Machine$double.eps)
  uu <- pl + runif(n) * width
  uu <- pmin(pmax(uu, .Machine$double.eps), 1 - .Machine$double.eps)
  
  mean + sd * qnorm(uu)
}

rtruncnorm_one <- function(mean, sd, lower = -Inf, upper = Inf) {
  rtruncnorm_vec(mean, sd, lower, upper)[1]
}

maximin_lhs_nd <- function(n, d, lower = -1, upper = 1) {
  X <- matrix(NA_real_, n, d)
  for (j in seq_len(d)) {
    z <- (seq_len(n) - runif(n)) / n
    X[, j] <- sample(z)
  }
  lower + (upper - lower) * X
}

pattern_key <- function(C) {
  C <- as.matrix(C)
  apply(C, 1, paste, collapse = "_")
}

make_pattern_label <- function(C) {
  apply(as.matrix(C), 1, function(z) paste0("(", paste(z, collapse = ","), ")"))
}

combine_chain_arrays <- function(arr_list) {
  dims <- dim(arr_list[[1]])
  total <- sum(vapply(arr_list, function(a) dim(a)[1], integer(1)))
  out <- array(NA_real_, dim = c(total, dims[-1]))
  pos <- 0L
  
  for (a in arr_list) {
    ns <- dim(a)[1]
    idx <- pos + seq_len(ns)
    
    if (length(dims) == 3L) {
      out[idx, , ] <- a
    } else if (length(dims) == 4L) {
      out[idx, , , ] <- a
    } else if (length(dims) == 2L) {
      out[idx, ] <- a
    } else {
      out[idx] <- a
    }
    
    pos <- pos + ns
  }
  
  out
}

############################################################
## Tau helpers
############################################################

tau_names_from_mvec <- function(m_vec) {
  unlist(
    lapply(seq_along(m_vec), function(j) {
      paste0("tau[", j, ",", seq_len(m_vec[j] - 1L), "]")
    }),
    use.names = FALSE
  )
}

flatten_tau <- function(tau_list) {
  unlist(tau_list, use.names = FALSE)
}

unflatten_tau <- function(tau_vec, m_vec) {
  tau_list <- vector("list", length(m_vec))
  pos <- 0L
  
  for (j in seq_along(m_vec)) {
    len <- m_vec[j] - 1L
    if (len > 0L) {
      tau_list[[j]] <- tau_vec[pos + seq_len(len)]
      pos <- pos + len
    } else {
      tau_list[[j]] <- numeric(0)
    }
  }
  
  tau_list
}

############################################################
## Slice sampler
############################################################

bounded_slice_update <- function(x0,
                                 logf,
                                 w = 1,
                                 lower = -Inf,
                                 upper = Inf,
                                 max_steps_out = 50,
                                 max_iter = 200) {
  if (upper <= lower) return(list(x = x0, n_eval = 0L))
  
  eps <- 1e-12
  if (is.finite(lower)) x0 <- max(x0, lower + eps)
  if (is.finite(upper)) x0 <- min(x0, upper - eps)
  
  f0 <- logf(x0)
  n_eval <- 1L
  
  if (!is.finite(f0)) return(list(x = x0, n_eval = n_eval))
  
  logy <- f0 + log(runif(1))
  L <- x0 - runif(1) * w
  R <- L + w
  
  if (is.finite(lower)) L <- max(L, lower)
  if (is.finite(upper)) R <- min(R, upper)
  
  J <- floor(runif(1) * max_steps_out)
  K <- max_steps_out - 1L - J
  
  while (J > 0 && (!is.finite(lower) || L > lower)) {
    fL <- logf(L)
    n_eval <- n_eval + 1L
    if (!is.finite(fL) || fL <= logy) break
    L <- L - w
    if (is.finite(lower)) L <- max(L, lower)
    J <- J - 1L
  }
  
  while (K > 0 && (!is.finite(upper) || R < upper)) {
    fR <- logf(R)
    n_eval <- n_eval + 1L
    if (!is.finite(fR) || fR <= logy) break
    R <- R + w
    if (is.finite(upper)) R <- min(R, upper)
    K <- K - 1L
  }
  
  for (iter in seq_len(max_iter)) {
    x1 <- runif(1, L, R)
    f1 <- logf(x1)
    n_eval <- n_eval + 1L
    
    if (is.finite(f1) && f1 >= logy) {
      return(list(x = x1, n_eval = n_eval))
    }
    
    if (x1 < x0) {
      L <- x1
    } else {
      R <- x1
    }
  }
  
  list(x = x0, n_eval = n_eval)
}

############################################################
## Ordinal-probit measurement model
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

score_ordinal_matrix <- function(C, level_scores) {
  C <- as.matrix(C)
  n <- nrow(C)
  q <- ncol(C)
  Z <- matrix(NA_real_, n, q)
  
  for (j in seq_len(q)) {
    Z[, j] <- level_scores[[j]][C[, j]]
  }
  
  Z
}

initialize_tau_ord <- function(C, m_vec, tau_bound = 6) {
  C <- as.matrix(C)
  q <- ncol(C)
  tau <- vector("list", q)
  
  for (j in seq_len(q)) {
    m_j <- m_vec[j]
    
    if (m_j <= 1L) {
      tau[[j]] <- numeric(0)
      next
    }
    
    counts <- tabulate(C[, j], nbins = m_j)
    counts_s <- counts + 0.5
    probs <- cumsum(counts_s)[seq_len(m_j - 1L)] / sum(counts_s)
    probs <- pmin(pmax(probs, 0.02), 0.98)
    tau_j <- qnorm(probs)
    tau_j <- pmin(pmax(tau_j, -tau_bound + 1e-3), tau_bound - 1e-3)
    tau[[j]] <- tau_j
  }
  
  tau
}

init_U_from_ordinal <- function(C,
                                d,
                                m_vec,
                                U_obs_full = NULL,
                                calib_idx = integer(0)) {
  C <- as.matrix(C)
  n <- nrow(C)
  
  level_scores <- empirical_normal_scores(C, m_vec = m_vec)
  Z <- score_ordinal_matrix(C, level_scores)
  Zs <- scale(Z)
  Zs[!is.finite(Zs)] <- 0
  
  pc <- prcomp(Zs, center = FALSE, scale. = FALSE)
  
  U_work <- matrix(0, n, d)
  k <- min(d, ncol(pc$x))
  
  if (k > 0L) {
    U_work[, seq_len(k)] <- pc$x[, seq_len(k), drop = FALSE]
  }
  
  if (d > k) {
    U_work[, (k + 1L):d] <- matrix(rnorm(n * (d - k)), n, d - k)
  }
  
  U_work <- scale(U_work)
  U_work[!is.finite(U_work)] <- 0
  
  calib_idx <- sort(as.integer(calib_idx))
  
  if (length(calib_idx) >= d + 1L && !is.null(U_obs_full)) {
    X_cal <- cbind(1, U_work[calib_idx, , drop = FALSE])
    Y_cal <- U_obs_full[calib_idx, , drop = FALSE]
    ridge <- diag(c(0, rep(1e-4, d)), d + 1L)
    Beta <- solve(crossprod(X_cal) + ridge, crossprod(X_cal, Y_cal))
    U_work <- cbind(1, U_work) %*% Beta
  }
  
  if (length(calib_idx) > 0L && !is.null(U_obs_full)) {
    U_work[calib_idx, ] <- U_obs_full[calib_idx, ]
  }
  
  U_work
}

initialize_A_ord <- function(C,
                             U,
                             tau,
                             m_vec,
                             ident = c("lower_triangular", "none")) {
  ident <- match.arg(ident)
  
  C <- as.matrix(C)
  U <- as.matrix(U)
  n <- nrow(C)
  q <- ncol(C)
  d <- ncol(U)
  
  A <- matrix(0, q, d)
  Z <- matrix(NA_real_, n, q)
  
  for (j in seq_len(q)) {
    tau_j <- tau[[j]]
    lower <- c(-Inf, tau_j)[C[, j]]
    upper <- c(tau_j, Inf)[C[, j]]
    denom <- pnorm(upper) - pnorm(lower)
    numer <- dnorm(lower) - dnorm(upper)
    Z[, j] <- numer / pmax(denom, .Machine$double.eps)
  }
  
  for (j in seq_len(q)) {
    fit <- try(lm.fit(x = cbind(1, U), y = Z[, j]), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      A[j, ] <- fit$coefficients[-1]
    } else {
      A[j, ] <- rnorm(d, 0, 0.2)
    }
  }
  
  if (ident == "lower_triangular") {
    if (q < d) {
      stop("lower_triangular identification requires q >= d.")
    }
    
    for (j in seq_len(d)) {
      if (j < d) {
        A[j, (j + 1L):d] <- 0
      }
      A[j, j] <- max(abs(A[j, j]), 0.2)
    }
  }
  
  A[!is.finite(A)] <- 0
  A
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

update_tau_ord <- function(tau, S, C, m_vec, tau_bound = 6) {
  tau_new <- tau
  q <- ncol(C)
  
  for (j in seq_len(q)) {
    m_j <- m_vec[j]
    if (m_j <= 1L) next
    
    for (r in seq_len(m_j - 1L)) {
      tau_j <- tau_new[[j]]
      
      L <- max(
        c(
          -tau_bound,
          if (r > 1L) tau_j[r - 1L] else -Inf,
          S[C[, j] <= r, j]
        ),
        na.rm = TRUE
      )
      
      U <- min(
        c(
          tau_bound,
          if (r < m_j - 1L) tau_j[r + 1L] else Inf,
          S[C[, j] > r, j]
        ),
        na.rm = TRUE
      )
      
      if (is.finite(L) && is.finite(U) && L < U) {
        tau_j[r] <- runif(1, L, U)
      }
      
      tau_new[[j]] <- tau_j
    }
  }
  
  tau_new
}

rmvnorm_single_truncated_last_positive <- function(mean, Sigma) {
  mean <- as.numeric(mean)
  r <- length(mean)
  
  if (r == 1L) {
    return(rtruncnorm_one(mean[1], sqrt(Sigma[1, 1]), lower = 0, upper = Inf))
  }
  
  idx_minus <- seq_len(r - 1L)
  idx_last <- r
  
  mu_l <- mean[idx_last]
  var_l <- Sigma[idx_last, idx_last]
  beta_l <- rtruncnorm_one(mu_l, sqrt(var_l), lower = 0, upper = Inf)
  
  mu_m <- mean[idx_minus]
  Sigma_mm <- Sigma[idx_minus, idx_minus, drop = FALSE]
  Sigma_ml <- Sigma[idx_minus, idx_last, drop = FALSE]
  
  cond_mean <- mu_m + as.numeric(Sigma_ml) / var_l * (beta_l - mu_l)
  cond_cov <- Sigma_mm - Sigma_ml %*% t(Sigma_ml) / var_l
  
  beta_m <- as.numeric(rmvnorm_chol(1, cond_mean, cond_cov))
  c(beta_m, beta_l)
}

update_A_ord <- function(S,
                         U,
                         s_A = 2.5,
                         ident = c("lower_triangular", "none")) {
  ident <- match.arg(ident)
  
  S <- as.matrix(S)
  U <- as.matrix(U)
  
  q <- ncol(S)
  d <- ncol(U)
  
  if (ident == "lower_triangular" && q < d) {
    stop("lower_triangular identification requires q >= d.")
  }
  
  A <- matrix(0, q, d)
  prior_prec <- 1 / s_A^2
  
  for (j in seq_len(q)) {
    if (ident == "lower_triangular" && j <= d) {
      active <- seq_len(j)
      Xj <- U[, active, drop = FALSE]
      Vj <- solve(crossprod(Xj) + prior_prec * diag(length(active)))
      mj <- as.numeric(Vj %*% crossprod(Xj, S[, j]))
      beta <- rmvnorm_single_truncated_last_positive(mj, Vj)
      
      A[j, active] <- beta
      
      if (j < d) {
        A[j, (j + 1L):d] <- 0
      }
    } else {
      Xj <- U
      Vj <- solve(crossprod(Xj) + prior_prec * diag(d))
      mj <- as.numeric(Vj %*% crossprod(Xj, S[, j]))
      A[j, ] <- rmvnorm_chol(1, mj, Vj)
    }
  }
  
  A
}

latent_reference_params <- function(S, A) {
  S <- as.matrix(S)
  A <- as.matrix(A)
  
  d <- ncol(A)
  V <- solve(diag(d) + crossprod(A))
  M <- S %*% A %*% V
  
  list(mean = M, V = V)
}

############################################################
## GP likelihood and prediction
############################################################

a_eps0 <- 2
b_eps0 <- 0.05

make_gp_prior <- function(p, d) {
  npar <- 1L + p + d
  list(
    mean = c(log(3), rep(log(0.5), p + d)),
    sd = c(1.5, rep(1.5, p + d)),
    lower = c(log(0.05), rep(log(1e-4), p + d)),
    upper = c(log(100), rep(log(100), p + d))
  )
}

log_prior_logtheta_gp <- function(logtheta, gp_prior) {
  if (length(logtheta) != length(gp_prior$mean)) return(-Inf)
  
  if (any(logtheta < gp_prior$lower) || any(logtheta > gp_prior$upper)) {
    return(-Inf)
  }
  
  sum(dnorm(logtheta, mean = gp_prior$mean, sd = gp_prior$sd, log = TRUE))
}

gp_corr_general <- function(X, U, logtheta) {
  X <- as.matrix(X)
  U <- as.matrix(U)
  
  p <- ncol(X)
  d <- ncol(U)
  
  expected_len <- 1L + p + d
  if (length(logtheta) != expected_len) {
    stop("logtheta has wrong length.")
  }
  
  rho <- exp(logtheta[1])
  theta_x <- exp(logtheta[1L + seq_len(p)])
  theta_u <- exp(logtheta[1L + p + seq_len(d)])
  
  n <- nrow(X)
  Rexp <- matrix(0, n, n)
  
  for (j in seq_len(p)) {
    Rexp <- Rexp + theta_x[j] * pairwise_sqdist(X[, j, drop = FALSE])
  }
  
  for (k in seq_len(d)) {
    Rexp <- Rexp + theta_u[k] * pairwise_sqdist(U[, k, drop = FALSE])
  }
  
  list(
    R = exp(-Rexp),
    rho = rho,
    theta_x = theta_x,
    theta_u = theta_u
  )
}

gp_state_general <- function(y, X, U, logtheta, sigma2_eps) {
  n <- length(y)
  
  cc <- gp_corr_general(X, U, logtheta)
  A_mat <- diag(n) + cc$rho^2 * cc$R
  Uchol <- safe_chol(A_mat)
  
  Ainv_y <- solve_chol(Uchol, y)
  logdetA <- 2 * sum(log(diag(Uchol)))
  quad <- sum(y * Ainv_y)
  
  loglik <- -0.5 * (
    n * log(2 * pi * sigma2_eps) +
      logdetA +
      quad / sigma2_eps
  )
  
  list(
    loglik = loglik,
    R = cc$R,
    A = A_mat,
    cholA = Uchol,
    Ainv_y = Ainv_y,
    logdetA = logdetA,
    quad = quad
  )
}

theta_logpost_integrated_general <- function(y, X, U, logtheta, gp_prior) {
  lp <- log_prior_logtheta_gp(logtheta, gp_prior)
  if (!is.finite(lp)) return(-Inf)
  
  n <- length(y)
  cc <- gp_corr_general(X, U, logtheta)
  A_mat <- diag(n) + cc$rho^2 * cc$R
  
  Uchol <- try(safe_chol(A_mat), silent = TRUE)
  if (inherits(Uchol, "try-error")) return(-Inf)
  
  Ainv_y <- solve_chol(Uchol, y)
  logdetA <- 2 * sum(log(diag(Uchol)))
  quad <- sum(y * Ainv_y)
  
  lp -
    0.5 * logdetA -
    (a_eps0 + n / 2) * log(b_eps0 + 0.5 * quad)
}

sample_sigma2_eps_general <- function(y, X, U, logtheta) {
  n <- length(y)
  
  cc <- gp_corr_general(X, U, logtheta)
  A_mat <- diag(n) + cc$rho^2 * cc$R
  Uchol <- safe_chol(A_mat)
  Ainv_y <- solve_chol(Uchol, y)
  quad <- sum(y * Ainv_y)
  
  shape <- a_eps0 + n / 2
  rate <- b_eps0 + 0.5 * quad
  
  1 / rgamma(1, shape = shape, rate = rate)
}

update_logtheta_slice_general <- function(y,
                                          X,
                                          U,
                                          logtheta,
                                          theta_slice_width,
                                          gp_prior) {
  n_eval_total <- 0L
  
  for (j in seq_along(logtheta)) {
    logf_j <- function(val) {
      lt <- logtheta
      lt[j] <- val
      theta_logpost_integrated_general(y, X, U, lt, gp_prior)
    }
    
    ans <- bounded_slice_update(
      x0 = logtheta[j],
      logf = logf_j,
      w = theta_slice_width[j],
      lower = gp_prior$lower[j],
      upper = gp_prior$upper[j],
      max_steps_out = 30L,
      max_iter = 100L
    )
    
    logtheta[j] <- ans$x
    n_eval_total <- n_eval_total + ans$n_eval
  }
  
  list(logtheta = logtheta, n_eval = n_eval_total)
}

gp_predict_draw_general <- function(X_train,
                                    U_train,
                                    y_train,
                                    X_star,
                                    U_star,
                                    logtheta,
                                    sigma2_eps,
                                    noisy = FALSE) {
  X_train <- as.matrix(X_train)
  U_train <- as.matrix(U_train)
  X_star <- as.matrix(X_star)
  U_star <- as.matrix(U_star)
  
  n <- nrow(X_train)
  N <- nrow(X_star)
  p <- ncol(X_train)
  d <- ncol(U_train)
  
  rho <- exp(logtheta[1])
  theta_x <- exp(logtheta[1L + seq_len(p)])
  theta_u <- exp(logtheta[1L + p + seq_len(d)])
  
  Rexp <- matrix(0, n, n)
  for (j in seq_len(p)) {
    Rexp <- Rexp + theta_x[j] * pairwise_sqdist(X_train[, j, drop = FALSE])
  }
  for (k in seq_len(d)) {
    Rexp <- Rexp + theta_u[k] * pairwise_sqdist(U_train[, k, drop = FALSE])
  }
  
  R <- exp(-Rexp)
  K <- rho^2 * sigma2_eps * R
  C <- K + sigma2_eps * diag(n)
  
  Uchol <- safe_chol(C)
  alpha <- solve_chol(Uchol, y_train)
  
  Rstar_exp <- matrix(0, N, n)
  for (j in seq_len(p)) {
    Rstar_exp <- Rstar_exp +
      theta_x[j] * pairwise_sqdist(
        X_star[, j, drop = FALSE],
        X_train[, j, drop = FALSE]
      )
  }
  for (k in seq_len(d)) {
    Rstar_exp <- Rstar_exp +
      theta_u[k] * pairwise_sqdist(
        U_star[, k, drop = FALSE],
        U_train[, k, drop = FALSE]
      )
  }
  
  Rstar <- exp(-Rstar_exp)
  Kstar <- rho^2 * sigma2_eps * Rstar
  
  mu <- as.numeric(Kstar %*% alpha)
  
  v <- forwardsolve(t(Uchol), t(Kstar))
  var_lat <- rho^2 * sigma2_eps - colSums(v^2)
  var_lat <- pmax(var_lat, 1e-10)
  
  if (noisy) var_lat <- var_lat + sigma2_eps
  
  list(mean = mu, var = var_lat)
}

############################################################
## Latent U update by elliptical slice
############################################################

update_U_ess_block_general <- function(y,
                                       X,
                                       U_curr,
                                       S,
                                       A,
                                       logtheta,
                                       sigma2_eps,
                                       block_idx,
                                       max_try = 300L) {
  block_idx <- as.integer(block_idx)
  
  if (length(block_idx) == 0L) {
    return(list(U = U_curr, n_eval = 0L, accepted = TRUE))
  }
  
  ref <- latent_reference_params(S, A)
  Mref <- ref$mean
  Vref <- ref$V
  
  U_block <- U_curr[block_idx, , drop = FALSE]
  M_block <- Mref[block_idx, , drop = FALSE]
  
  loglik_fun <- function(U_block_prop) {
    U_prop <- U_curr
    U_prop[block_idx, ] <- U_block_prop
    gp_state_general(y, X, U_prop, logtheta, sigma2_eps)$loglik
  }
  
  loglik_cur <- loglik_fun(U_block)
  
  if (!is.finite(loglik_cur)) {
    return(list(U = U_curr, n_eval = 1L, accepted = FALSE))
  }
  
  Nu <- rmvnorm_rows_common(
    mean_mat = matrix(0, nrow = length(block_idx), ncol = ncol(U_curr)),
    Sigma = Vref
  )
  
  Zcur <- U_block - M_block
  logy <- loglik_cur + log(runif(1))
  
  angle <- runif(1, 0, 2 * pi)
  angle_min <- angle - 2 * pi
  angle_max <- angle
  
  n_eval <- 1L
  
  for (try_id in seq_len(max_try)) {
    U_prop_block <- M_block + Zcur * cos(angle) + Nu * sin(angle)
    loglik_prop <- loglik_fun(U_prop_block)
    n_eval <- n_eval + 1L
    
    if (is.finite(loglik_prop) && loglik_prop >= logy) {
      U_new <- U_curr
      U_new[block_idx, ] <- U_prop_block
      return(list(U = U_new, n_eval = n_eval, accepted = TRUE))
    }
    
    if (angle < 0) {
      angle_min <- angle
    } else {
      angle_max <- angle
    }
    
    angle <- runif(1, angle_min, angle_max)
  }
  
  list(U = U_curr, n_eval = n_eval, accepted = FALSE)
}

############################################################
## MCMC diagnostics
############################################################

split_rhat <- function(chain_list) {
  chain_list <- lapply(chain_list, function(x) as.numeric(x[is.finite(x)]))
  lens <- vapply(chain_list, length, integer(1))
  n0 <- min(lens)
  
  if (length(chain_list) < 2 || n0 < 20) return(NA_real_)
  
  chain_list <- lapply(chain_list, function(x) tail(x, n0))
  n_half <- floor(n0 / 2)
  
  if (n_half < 10) return(NA_real_)
  
  split_mat <- do.call(
    rbind,
    lapply(chain_list, function(x) {
      rbind(
        x[seq_len(n_half)],
        x[(n0 - n_half + 1):n0]
      )
    })
  )
  
  n_split <- ncol(split_mat)
  chain_means <- rowMeans(split_mat)
  chain_vars <- apply(split_mat, 1, var)
  
  W <- mean(chain_vars)
  B <- n_split * var(chain_means)
  
  if (!is.finite(W) || W <= 0) return(NA_real_)
  
  var_hat <- ((n_split - 1) / n_split) * W + B / n_split
  sqrt(var_hat / W)
}

ess_ips <- function(x, max_lag = NULL) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  
  if (n < 5) return(NA_real_)
  
  sx <- sd(x)
  if (!is.finite(sx)) return(NA_real_)
  if (sx == 0) return(n)
  
  if (is.null(max_lag)) {
    max_lag <- min(n - 1L, 1000L)
  } else {
    max_lag <- min(n - 1L, as.integer(max_lag))
  }
  
  if (max_lag < 1L) return(n)
  
  ac <- tryCatch(
    as.numeric(stats::acf(
      x,
      lag.max = max_lag,
      plot = FALSE,
      demean = TRUE
    )$acf),
    error = function(e) NA_real_
  )
  
  if (length(ac) <= 1 || all(!is.finite(ac))) return(NA_real_)
  
  ac <- ac[-1]
  ac <- ac[is.finite(ac)]
  if (length(ac) == 0) return(n)
  
  mm <- floor(length(ac) / 2)
  
  if (mm >= 1) {
    pair_sums <- ac[2 * seq_len(mm) - 1] + ac[2 * seq_len(mm)]
    first_nonpos <- which(pair_sums <= 0)[1]
    
    if (is.na(first_nonpos)) {
      use_lag <- 2 * mm
    } else {
      use_lag <- 2 * (first_nonpos - 1)
    }
  } else {
    use_lag <- ifelse(ac[1] > 0, 1, 0)
  }
  
  if (use_lag > 0) {
    tau_int <- 1 + 2 * sum(ac[seq_len(use_lag)])
  } else {
    tau_int <- 1
  }
  
  tau_int <- max(tau_int, 1)
  ess <- n / tau_int
  min(max(ess, 1), n)
}

classify_parameter_block <- function(parameter) {
  out <- rep("other", length(parameter))
  out[grepl("^sigma|^rho|^theta", parameter)] <- "GP/noise"
  out[grepl("^A\\[|^A[0-9]", parameter)] <- "loading A"
  out[grepl("^tau\\[|^tau", parameter)] <- "cutpoints"
  out[grepl("^U\\[", parameter)] <- "latent U"
  out
}

############################################################
## Exact fully Bayesian ordinal-probit EIV-GP sampler
############################################################

fit_eivgp_ordprobit_fb <- function(X_raw,
                                   y_raw,
                                   C_ord,
                                   U_obs = NULL,
                                   calib_idx = integer(0),
                                   U_true_eval = NULL,
                                   d = 2L,
                                   m_vec = NULL,
                                   ident = c("lower_triangular", "none"),
                                   n_iter = 3000L,
                                   burn = 1000L,
                                   thin = 2L,
                                   n_chains = 4L,
                                   preset = "balanced",
                                   seed = 1L,
                                   parallel_chains = FALSE,
                                   verbose = FALSE) {
  ident <- match.arg(ident)
  set.seed(seed)
  
  X_raw <- as.matrix(X_raw)
  C_ord <- as.matrix(C_ord)
  y_raw <- as.numeric(y_raw)
  
  if (anyNA(X_raw) || anyNA(y_raw) || anyNA(C_ord)) {
    stop("This implementation currently requires complete X, y, and C.")
  }
  
  n <- length(y_raw)
  p <- ncol(X_raw)
  q <- ncol(C_ord)
  
  if (is.null(m_vec)) {
    m_vec <- apply(C_ord, 2, max)
  }
  
  m_vec <- as.integer(m_vec)
  
  if (length(m_vec) != q) {
    stop("m_vec must have length equal to ncol(C_ord).")
  }
  
  for (j in seq_len(q)) {
    if (any(C_ord[, j] < 1) || any(C_ord[, j] > m_vec[j])) {
      stop("C_ord contains levels outside 1:m_j.")
    }
  }
  
  if (ident == "lower_triangular" && q < d) {
    stop("lower_triangular identification requires q >= d.")
  }
  
  calib_idx <- sort(as.integer(calib_idx))
  U_obs_full <- matrix(NA_real_, n, d)
  
  if (!is.null(U_obs)) {
    U_obs <- as.matrix(U_obs)
    
    if (nrow(U_obs) == n) {
      if (ncol(U_obs) != d) stop("U_obs has wrong number of columns.")
      
      if (length(calib_idx) == 0L) {
        calib_idx <- which(stats::complete.cases(U_obs))
      }
      
      U_obs_full[calib_idx, ] <- U_obs[calib_idx, , drop = FALSE]
    } else if (nrow(U_obs) == length(calib_idx)) {
      if (ncol(U_obs) != d) stop("U_obs has wrong number of columns.")
      U_obs_full[calib_idx, ] <- U_obs
    } else {
      stop("U_obs must have either n rows or length(calib_idx) rows.")
    }
  } else {
    if (length(calib_idx) > 0L) {
      stop("calib_idx was supplied but U_obs is NULL.")
    }
  }
  
  miss_idx <- setdiff(seq_len(n), calib_idx)
  
  X_center <- colMeans(X_raw)
  X_scale <- apply(X_raw, 2, sd)
  X_scale[!is.finite(X_scale) | X_scale <= 0] <- 1
  X <- sweep(sweep(X_raw, 2, X_center, "-"), 2, X_scale, "/")
  
  y_center <- mean(y_raw)
  y_scale <- sd(y_raw)
  if (!is.finite(y_scale) || y_scale <= 0) y_scale <- 1
  y <- as.numeric((y_raw - y_center) / y_scale)
  
  control <- make_default_control_ordprobit(
    n = n,
    n_mis = length(miss_idx),
    preset = preset,
    d = d
  )
  
  gp_prior <- make_gp_prior(p = p, d = d)
  n_logtheta <- 1L + p + d
  control$theta_slice_width_init <- rep(0.8, n_logtheta)
  
  n_save <- floor((n_iter - burn) / thin)
  if (n_save <= 0L) stop("n_iter, burn, and thin imply no saved draws.")
  
  tau_flat_names <- tau_names_from_mvec(m_vec)
  tau_total <- length(tau_flat_names)
  
  initialize_chain_state <- function(chain_seed) {
    set.seed(chain_seed)
    
    tau0 <- initialize_tau_ord(
      C = C_ord,
      m_vec = m_vec,
      tau_bound = control$tau_bound
    )
    
    U0 <- init_U_from_ordinal(
      C = C_ord,
      d = d,
      m_vec = m_vec,
      U_obs_full = U_obs_full,
      calib_idx = calib_idx
    )
    
    U0 <- U0 + matrix(rnorm(n * d, 0, 0.15), n, d)
    
    if (length(calib_idx) > 0L) {
      U0[calib_idx, ] <- U_obs_full[calib_idx, ]
    }
    
    A0 <- initialize_A_ord(
      C = C_ord,
      U = U0,
      tau = tau0,
      m_vec = m_vec,
      ident = ident
    )
    
    S0 <- sample_scores_ord(C_ord, U0, A0, tau0)
    
    sigma2_eps0 <- exp(log(0.05) + rnorm(1, 0, 0.5))
    sigma2_eps0 <- min(max(sigma2_eps0, 1e-4), 2)
    
    logtheta0 <- gp_prior$mean + rnorm(n_logtheta, 0, 0.5)
    logtheta0 <- pmin(pmax(logtheta0, gp_prior$lower), gp_prior$upper)
    
    list(
      U_curr = U0,
      S_curr = S0,
      A_curr = A0,
      tau_curr = tau0,
      sigma2_eps = sigma2_eps0,
      logtheta = logtheta0,
      theta_slice_width = control$theta_slice_width_init
    )
  }
  
  run_one_chain <- function(chain_id, chain_seed) {
    set.seed(chain_seed)
    
    state <- initialize_chain_state(chain_seed)
    
    U_curr <- state$U_curr
    S_curr <- state$S_curr
    A_curr <- state$A_curr
    tau_curr <- state$tau_curr
    sigma2_eps <- state$sigma2_eps
    logtheta <- state$logtheta
    theta_slice_width <- state$theta_slice_width
    
    samples_U_chain <- array(NA_real_, dim = c(n_save, n, d))
    samples_S_chain <- array(NA_real_, dim = c(n_save, n, q))
    samples_A_chain <- array(NA_real_, dim = c(n_save, q, d))
    samples_tau_chain <- matrix(NA_real_, n_save, tau_total)
    colnames(samples_tau_chain) <- tau_flat_names
    samples_logtheta_chain <- matrix(NA_real_, n_save, n_logtheta)
    samples_sigma2_chain <- numeric(n_save)
    
    logtheta_trace_all <- matrix(NA_real_, n_iter, n_logtheta)
    
    u_ess_eval_total <- 0L
    u_ess_accept_total <- 0L
    u_ess_total <- 0L
    
    global_u_eval_total <- 0L
    global_u_accept_total <- 0L
    global_u_total <- 0L
    
    theta_eval_total <- 0L
    theta_update_total <- 0L
    
    save_id <- 0L
    
    for (iter in seq_len(n_iter)) {
      if (iter %% control$score_update_every == 0L) {
        S_curr <- sample_scores_ord(C_ord, U_curr, A_curr, tau_curr)
      }
      
      if (iter %% control$tau_update_every == 0L) {
        tau_curr <- update_tau_ord(
          tau = tau_curr,
          S = S_curr,
          C = C_ord,
          m_vec = m_vec,
          tau_bound = control$tau_bound
        )
      }
      
      if (iter %% control$A_update_every == 0L) {
        A_curr <- update_A_ord(
          S = S_curr,
          U = U_curr,
          s_A = control$s_A,
          ident = ident
        )
      }
      
      if (length(miss_idx) > 0L) {
        for (bb in seq_len(control$n_u_blocks_per_iter)) {
          block_idx <- sample(
            miss_idx,
            min(control$u_block_size, length(miss_idx))
          )
          
          uu <- update_U_ess_block_general(
            y = y,
            X = X,
            U_curr = U_curr,
            S = S_curr,
            A = A_curr,
            logtheta = logtheta,
            sigma2_eps = sigma2_eps,
            block_idx = block_idx,
            max_try = control$max_ess_try
          )
          
          U_curr <- uu$U
          u_ess_eval_total <- u_ess_eval_total + uu$n_eval
          u_ess_accept_total <- u_ess_accept_total + as.integer(uu$accepted)
          u_ess_total <- u_ess_total + 1L
        }
        
        if (control$global_u_every > 0L &&
            iter %% control$global_u_every == 0L) {
          gu <- update_U_ess_block_general(
            y = y,
            X = X,
            U_curr = U_curr,
            S = S_curr,
            A = A_curr,
            logtheta = logtheta,
            sigma2_eps = sigma2_eps,
            block_idx = miss_idx,
            max_try = control$max_ess_try
          )
          
          U_curr <- gu$U
          global_u_eval_total <- global_u_eval_total + gu$n_eval
          global_u_accept_total <- global_u_accept_total + as.integer(gu$accepted)
          global_u_total <- global_u_total + 1L
        }
      }
      
      if (length(calib_idx) > 0L) {
        U_curr[calib_idx, ] <- U_obs_full[calib_idx, ]
      }
      
      if (iter %% control$theta_update_every == 0L) {
        th <- update_logtheta_slice_general(
          y = y,
          X = X,
          U = U_curr,
          logtheta = logtheta,
          theta_slice_width = theta_slice_width,
          gp_prior = gp_prior
        )
        
        logtheta <- th$logtheta
        theta_eval_total <- theta_eval_total + th$n_eval
        theta_update_total <- theta_update_total + 1L
      }
      
      sigma2_eps <- sample_sigma2_eps_general(y, X, U_curr, logtheta)
      logtheta_trace_all[iter, ] <- logtheta
      
      if (control$adapt_theta_width &&
          iter <= burn &&
          iter >= 200 &&
          iter %% control$adapt_every == 0L) {
        lo <- max(1, iter - control$adapt_window + 1)
        recent <- logtheta_trace_all[lo:iter, , drop = FALSE]
        recent <- recent[complete.cases(recent), , drop = FALSE]
        
        if (nrow(recent) >= 50) {
          sds <- apply(recent, 2, sd)
          new_width <- 2 * sds
          new_width <- pmin(
            pmax(new_width, control$theta_width_min),
            control$theta_width_max
          )
          
          if (all(is.finite(new_width))) {
            theta_slice_width <- new_width
          }
        }
      }
      
      if (iter > burn && ((iter - burn) %% thin == 0L)) {
        save_id <- save_id + 1L
        
        samples_U_chain[save_id, , ] <- U_curr
        samples_S_chain[save_id, , ] <- S_curr
        samples_A_chain[save_id, , ] <- A_curr
        samples_tau_chain[save_id, ] <- flatten_tau(tau_curr)
        samples_logtheta_chain[save_id, ] <- logtheta
        samples_sigma2_chain[save_id] <- sigma2_eps
      }
    }
    
    if (save_id < n_save) {
      samples_U_chain <- samples_U_chain[seq_len(save_id), , , drop = FALSE]
      samples_S_chain <- samples_S_chain[seq_len(save_id), , , drop = FALSE]
      samples_A_chain <- samples_A_chain[seq_len(save_id), , , drop = FALSE]
      samples_tau_chain <- samples_tau_chain[seq_len(save_id), , drop = FALSE]
      samples_logtheta_chain <- samples_logtheta_chain[seq_len(save_id), , drop = FALSE]
      samples_sigma2_chain <- samples_sigma2_chain[seq_len(save_id)]
    }
    
    stats <- data.frame(
      chain = chain_id,
      seed = chain_seed,
      saved = save_id,
      u_ess_eval_total = u_ess_eval_total,
      u_ess_accept_total = u_ess_accept_total,
      u_ess_total = u_ess_total,
      global_u_eval_total = global_u_eval_total,
      global_u_accept_total = global_u_accept_total,
      global_u_total = global_u_total,
      theta_eval_total = theta_eval_total,
      theta_update_total = theta_update_total
    )
    
    list(
      chain_id = chain_id,
      seed = chain_seed,
      samples_U = samples_U_chain,
      samples_S = samples_S_chain,
      samples_A = samples_A_chain,
      samples_tau = samples_tau_chain,
      samples_logtheta = samples_logtheta_chain,
      samples_sigma2 = samples_sigma2_chain,
      theta_slice_width_final = theta_slice_width,
      stats = stats
    )
  }
  
  chain_seeds <- seed + 10000L * seq_len(n_chains)
  
  use_mclapply <- (
    parallel_chains &&
      .Platform$OS.type != "windows" &&
      n_chains > 1L
  )
  
  mc_cores <- min(
    n_chains,
    max(1L, parallel::detectCores(logical = TRUE) - 2L)
  )
  
  if (verbose) {
    cat("Running exact ordinal-probit EIV-GP with", n_chains, "chain(s).\n")
  }
  
  mcmc_time <- system.time({
    if (use_mclapply) {
      chains <- parallel::mclapply(
        seq_len(n_chains),
        function(cc) {
          run_one_chain(chain_id = cc, chain_seed = chain_seeds[cc])
        },
        mc.cores = mc_cores,
        mc.set.seed = TRUE,
        mc.preschedule = FALSE
      )
    } else {
      chains <- lapply(
        seq_len(n_chains),
        function(cc) {
          run_one_chain(chain_id = cc, chain_seed = chain_seeds[cc])
        }
      )
    }
  })
  
  samples_by_chain <- list(
    U = lapply(chains, function(z) z$samples_U),
    S = lapply(chains, function(z) z$samples_S),
    A = lapply(chains, function(z) z$samples_A),
    tau = lapply(chains, function(z) z$samples_tau),
    logtheta = lapply(chains, function(z) z$samples_logtheta),
    sigma2 = lapply(chains, function(z) z$samples_sigma2)
  )
  
  samples_U <- combine_chain_arrays(samples_by_chain$U)
  samples_S <- combine_chain_arrays(samples_by_chain$S)
  samples_A <- combine_chain_arrays(samples_by_chain$A)
  samples_tau <- do.call(rbind, samples_by_chain$tau)
  samples_logtheta <- do.call(rbind, samples_by_chain$logtheta)
  samples_sigma2 <- unlist(samples_by_chain$sigma2)
  
  mcmc_draw_info <- data.frame(
    chain = rep(
      seq_len(n_chains),
      times = vapply(samples_by_chain$logtheta, nrow, integer(1))
    ),
    draw_within_chain = unlist(
      lapply(samples_by_chain$logtheta, function(mat) seq_len(nrow(mat)))
    )
  )
  
  chain_stats <- do.call(rbind, lapply(chains, function(z) z$stats))
  
  hyper_names <- c(
    "sigma_epsilon",
    "rho",
    paste0("theta_x", seq_len(p)),
    paste0("theta_u", seq_len(d))
  )
  
  rhat_hyper <- data.frame(
    parameter = hyper_names,
    rhat = c(
      split_rhat(lapply(samples_by_chain$sigma2, function(v) sqrt(v))),
      sapply(seq_len(n_logtheta), function(k) {
        split_rhat(lapply(samples_by_chain$logtheta, function(mat) exp(mat[, k])))
      })
    )
  )
  
  rhat_A <- do.call(
    rbind,
    lapply(seq_len(q), function(j) {
      do.call(
        rbind,
        lapply(seq_len(d), function(k) {
          data.frame(
            parameter = paste0("A[", j, ",", k, "]"),
            rhat = split_rhat(lapply(samples_by_chain$A, function(arr) arr[, j, k]))
          )
        })
      )
    })
  )
  
  rhat_tau <- data.frame(
    parameter = colnames(samples_tau),
    rhat = sapply(seq_len(ncol(samples_tau)), function(k) {
      split_rhat(lapply(samples_by_chain$tau, function(mat) mat[, k]))
    })
  )
  
  if (length(miss_idx) > 0L) {
    rhat_U <- do.call(
      rbind,
      lapply(miss_idx, function(ii) {
        do.call(
          rbind,
          lapply(seq_len(d), function(k) {
            data.frame(
              parameter = paste0("U[", ii, ",", k, "]"),
              global_index = ii,
              coord = k,
              rhat = split_rhat(lapply(samples_by_chain$U, function(arr) arr[, ii, k]))
            )
          })
        )
      })
    )
  } else {
    rhat_U <- data.frame(
      parameter = character(0),
      global_index = integer(0),
      coord = integer(0),
      rhat = numeric(0)
    )
  }
  
  hyper_mat <- cbind(
    sigma_epsilon = sqrt(samples_sigma2),
    rho = exp(samples_logtheta[, 1])
  )
  
  if (p > 0L) {
    tmp_x <- exp(samples_logtheta[, 1L + seq_len(p), drop = FALSE])
    colnames(tmp_x) <- paste0("theta_x", seq_len(p))
    hyper_mat <- cbind(hyper_mat, tmp_x)
  }
  
  if (d > 0L) {
    tmp_u <- exp(samples_logtheta[, 1L + p + seq_len(d), drop = FALSE])
    colnames(tmp_u) <- paste0("theta_u", seq_len(d))
    hyper_mat <- cbind(hyper_mat, tmp_u)
  }
  
  n_samp_A <- dim(samples_A)[1]
  A_mat <- matrix(NA_real_, nrow = n_samp_A, ncol = q * d)
  A_names <- character(q * d)
  col_id <- 0L
  
  for (j in seq_len(q)) {
    for (k in seq_len(d)) {
      col_id <- col_id + 1L
      A_mat[, col_id] <- samples_A[, j, k]
      A_names[col_id] <- paste0("A", j, k)
    }
  }
  
  colnames(A_mat) <- A_names
  
  ess_input <- cbind(hyper_mat, A_mat, samples_tau)
  ess_key <- data.frame(
    parameter = colnames(ess_input),
    ess = apply(ess_input, 2, ess_ips)
  )
  
  diagnostics_summary <- data.frame(
    n_chains = n_chains,
    n_iter = n_iter,
    burn = burn,
    thin = thin,
    saved_per_chain = mean(vapply(samples_by_chain$logtheta, nrow, integer(1))),
    total_saved_draws = nrow(samples_logtheta),
    max_rhat_hyper = safe_max(rhat_hyper$rhat),
    max_rhat_A = safe_max(rhat_A$rhat),
    max_rhat_tau = safe_max(rhat_tau$rhat),
    median_rhat_missing_U = safe_median(rhat_U$rhat),
    max_rhat_missing_U = safe_max(rhat_U$rhat),
    min_ess_key = safe_min(ess_key$ess),
    mean_u_ess_accept = sum(chain_stats$u_ess_accept_total) /
      max(sum(chain_stats$u_ess_total), 1),
    mean_global_u_accept = sum(chain_stats$global_u_accept_total) /
      max(sum(chain_stats$global_u_total), 1),
    time_seconds = as.numeric(mcmc_time["elapsed"])
  )
  
  list(
    data = list(
      X_raw = X_raw,
      X = X,
      X_center = X_center,
      X_scale = X_scale,
      y_raw = y_raw,
      y = y,
      y_center = y_center,
      y_scale = y_scale,
      C_ord = C_ord,
      U_obs = U_obs_full,
      U_true_eval = U_true_eval,
      calib_idx = calib_idx,
      miss_idx = miss_idx,
      m_vec = m_vec,
      q = q,
      d = d,
      p = p,
      ident = ident
    ),
    gp_prior = gp_prior,
    control = control,
    mcmc = list(
      samples_U = samples_U,
      samples_S = samples_S,
      samples_A = samples_A,
      samples_tau = samples_tau,
      samples_logtheta = samples_logtheta,
      samples_sigma2 = samples_sigma2,
      samples_by_chain = samples_by_chain,
      mcmc_draw_info = mcmc_draw_info,
      chain_stats = chain_stats
    ),
    diagnostics = list(
      rhat_hyper = rhat_hyper,
      rhat_A = rhat_A,
      rhat_tau = rhat_tau,
      rhat_U = rhat_U,
      ess_key = ess_key,
      summary = diagnostics_summary
    )
  )
}

############################################################
## Prediction for ordinal-probit EIV-GP
############################################################

sample_u_given_c_ordprobit <- function(C_new,
                                       A,
                                       tau,
                                       m_vec,
                                       n_gibbs = 20L,
                                       U_init = NULL) {
  C_new <- as.matrix(C_new)
  A <- as.matrix(A)
  
  N <- nrow(C_new)
  q <- ncol(C_new)
  d <- ncol(A)
  
  if (length(m_vec) != q) {
    stop("m_vec has wrong length.")
  }
  
  if (is.null(U_init)) {
    U <- matrix(rnorm(N * d), N, d)
  } else {
    U <- as.matrix(U_init)
  }
  
  S <- matrix(0, N, q)
  V <- solve(diag(d) + crossprod(A))
  
  for (gg in seq_len(n_gibbs)) {
    for (j in seq_len(q)) {
      mu <- as.numeric(U %*% A[j, ])
      tau_j <- tau[[j]]
      lower <- c(-Inf, tau_j)[C_new[, j]]
      upper <- c(tau_j, Inf)[C_new[, j]]
      S[, j] <- rtruncnorm_vec(mu, 1, lower, upper)
    }
    
    M <- S %*% A %*% V
    U <- rmvnorm_rows_common(M, V)
  }
  
  U
}

sample_eiv_test_y_ordprobit_fb <- function(X_test_raw,
                                           C_test,
                                           fit_obj,
                                           draw_ids,
                                           n_per_draw = 1L,
                                           n_new_latent_gibbs = 20L) {
  X_test_raw <- as.matrix(X_test_raw)
  C_test <- as.matrix(C_test)
  
  samples_U <- fit_obj$mcmc$samples_U
  samples_A <- fit_obj$mcmc$samples_A
  samples_tau <- fit_obj$mcmc$samples_tau
  samples_logtheta <- fit_obj$mcmc$samples_logtheta
  samples_sigma2 <- fit_obj$mcmc$samples_sigma2
  
  X_train <- fit_obj$data$X
  y_train <- fit_obj$data$y
  y_center <- fit_obj$data$y_center
  y_scale <- fit_obj$data$y_scale
  X_center <- fit_obj$data$X_center
  X_scale <- fit_obj$data$X_scale
  
  m_vec <- fit_obj$data$m_vec
  d <- fit_obj$data$d
  q <- fit_obj$data$q
  n_train <- nrow(X_train)
  
  X_test <- sweep(sweep(X_test_raw, 2, X_center, "-"), 2, X_scale, "/")
  n_test <- nrow(C_test)
  
  n_draw <- length(draw_ids) * n_per_draw
  out <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  
  row_id <- 0L
  
  for (s in draw_ids) {
    A_s <- matrix(samples_A[s, , ], nrow = q, ncol = d)
    tau_s <- unflatten_tau(samples_tau[s, ], m_vec)
    U_train_s <- matrix(samples_U[s, , ], nrow = n_train, ncol = d)
    
    for (rr in seq_len(n_per_draw)) {
      row_id <- row_id + 1L
      
      U_star <- sample_u_given_c_ordprobit(
        C_new = C_test,
        A = A_s,
        tau = tau_s,
        m_vec = m_vec,
        n_gibbs = n_new_latent_gibbs
      )
      
      pred <- gp_predict_draw_general(
        X_train = X_train,
        U_train = U_train_s,
        y_train = y_train,
        X_star = X_test,
        U_star = U_star,
        logtheta = samples_logtheta[s, ],
        sigma2_eps = samples_sigma2[s],
        noisy = TRUE
      )
      
      y_std <- pred$mean + sqrt(pred$var) * rnorm(n_test)
      out[row_id, ] <- y_center + y_scale * y_std
    }
  }
  
  out
}

############################################################
## Revised Study II data-generating mechanism
############################################################

make_balanced_tau_from_A <- function(A,
                                     probs = c(0.25, 0.50, 0.75),
                                     residual_sd = 1) {
  A <- as.matrix(A)
  score_sd <- sqrt(residual_sd^2 + rowSums(A^2))
  outer(score_sd, qnorm(probs))
}

make_study2_true_params <- function() {
  A <- rbind(
    c(1.70, 0.00),
    c(0.20, 1.70),
    c(1.20, 0.70),
    c(0.70, 1.20)
  )
  
  Omega <- diag(4)
  
  tau <- make_balanced_tau_from_A(
    A = A,
    probs = c(0.25, 0.50, 0.75),
    residual_sd = 1
  )
  
  list(
    A = A,
    Omega = Omega,
    tau = tau,
    q = 4L,
    m = 4L,
    d = 2L,
    sigma_eps = 0.12
  )
}

f0_2d <- function(X, U) {
  X <- as.matrix(X)
  U <- as.matrix(U)
  
  x1 <- X[, 1]
  x2 <- X[, 2]
  u1 <- U[, 1]
  u2 <- U[, 2]
  
  0.55 * sin(pi * x1 / 2) +
    0.35 * x2 +
    0.85 * tanh(u1) +
    0.70 * tanh(u2) +
    0.20 * tanh(0.80 * u1 * u2) +
    0.15 * x1 * tanh(u1) -
    0.12 * x2 * tanh(u2)
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

simulate_study2_data <- function(n = 120,
                                 n_test = 400,
                                 seed = NULL,
                                 sigma_eps = NULL,
                                 min_count_per_level = 6L) {
  if (!is.null(seed)) set.seed(seed)
  
  pars <- make_study2_true_params()
  A <- pars$A
  Omega <- pars$Omega
  tau <- pars$tau
  
  if (is.null(sigma_eps)) sigma_eps <- pars$sigma_eps
  
  repeat {
    X <- maximin_lhs_nd(n, d = 2, lower = -1, upper = 1)
    U <- matrix(rnorm(n * 2), n, 2)
    
    Zeta <- rmvnorm_chol(n, rep(0, 4), Omega)
    S <- U %*% t(A) + Zeta
    C <- ordinal_from_scores_matrix_tau(S, tau)
    
    ok <- all(apply(C, 2, function(z) {
      all(tabulate(z, nbins = 4) >= min_count_per_level)
    }))
    
    if (ok) break
  }
  
  f <- f0_2d(X, U)
  y <- f + rnorm(n, 0, sigma_eps)
  
  X_test <- matrix(runif(n_test * 2, -1, 1), n_test, 2)
  U_test <- matrix(rnorm(n_test * 2), n_test, 2)
  
  Zeta_test <- rmvnorm_chol(n_test, rep(0, 4), Omega)
  S_test <- U_test %*% t(A) + Zeta_test
  C_test <- ordinal_from_scores_matrix_tau(S_test, tau)
  
  f_test <- f0_2d(X_test, U_test)
  y_test <- f_test + rnorm(n_test, 0, sigma_eps)
  
  list(
    train = list(
      X = X,
      U = U,
      C = C,
      y = y,
      f = f
    ),
    test = list(
      X = X_test,
      U = U_test,
      C = C_test,
      y = y_test,
      f = f_test
    ),
    true_params = pars,
    sigma_eps = sigma_eps
  )
}

############################################################
## Cell-stratified calibration design
############################################################

make_cell_stratified_calibration_sets_2d <- function(C,
                                                     calib_grid = c(0, 10, 25, 50, 80),
                                                     anchor_cols = c(1L, 2L),
                                                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  C <- as.matrix(C)
  n <- nrow(C)
  anchor_cols <- as.integer(anchor_cols)
  
  if (length(anchor_cols) != 2L) {
    stop("anchor_cols must contain exactly two ordinal proxy columns.")
  }
  
  if (any(anchor_cols < 1L) || any(anchor_cols > ncol(C))) {
    stop("anchor_cols contains invalid column indices.")
  }
  
  cell_key <- apply(C[, anchor_cols, drop = FALSE], 1, paste, collapse = "_")
  cell_names <- sort(unique(cell_key))
  cell_names <- sample(cell_names, length(cell_names))
  
  idx_by_cell <- lapply(cell_names, function(cc) {
    sample(which(cell_key == cc))
  })
  names(idx_by_cell) <- cell_names
  
  ordering <- integer(0)
  
  repeat {
    active <- which(vapply(idx_by_cell, length, integer(1)) > 0L)
    if (length(active) == 0L) break
    
    active <- sample(active, length(active))
    
    for (aa in active) {
      ordering <- c(ordering, idx_by_cell[[aa]][1])
      idx_by_cell[[aa]] <- idx_by_cell[[aa]][-1]
    }
  }
  
  if (length(ordering) != n) {
    missing_idx <- setdiff(seq_len(n), ordering)
    ordering <- c(ordering, sample(missing_idx))
  }
  
  out <- lapply(calib_grid, function(k) {
    k <- min(as.integer(k), n)
    if (k <= 0L) integer(0) else sort(ordering[seq_len(k)])
  })
  
  names(out) <- as.character(calib_grid)
  out
}

make_stratified_calibration_sets_2d <- function(C,
                                                calib_grid = c(0, 10, 25, 50, 80),
                                                seed = NULL,
                                                anchor_cols = c(1L, 2L)) {
  make_cell_stratified_calibration_sets_2d(
    C = C,
    calib_grid = calib_grid,
    anchor_cols = anchor_cols,
    seed = seed
  )
}

############################################################
## Oracle prediction
############################################################

make_oracle_pool_2d <- function(true_params,
                                n_pool = 200000L,
                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  A <- true_params$A
  Omega <- true_params$Omega
  tau <- true_params$tau
  
  U <- matrix(rnorm(n_pool * 2), n_pool, 2)
  Zeta <- rmvnorm_chol(n_pool, rep(0, 4), Omega)
  S <- U %*% t(A) + Zeta
  C <- ordinal_from_scores_matrix_tau(S, tau)
  
  key <- pattern_key(C)
  split_idx <- split(seq_len(n_pool), key)
  
  list(
    U = U,
    C = C,
    key = key,
    split_idx = split_idx
  )
}

sample_oracle_u_rejection <- function(c_star,
                                      true_params,
                                      n_draw,
                                      max_batches = 100L,
                                      batch_size = 20000L) {
  A <- true_params$A
  Omega <- true_params$Omega
  tau <- true_params$tau
  
  U_keep <- matrix(NA_real_, 0, 2)
  key_star <- paste(c_star, collapse = "_")
  
  for (bb in seq_len(max_batches)) {
    U <- matrix(rnorm(batch_size * 2), batch_size, 2)
    Zeta <- rmvnorm_chol(batch_size, rep(0, 4), Omega)
    S <- U %*% t(A) + Zeta
    C <- ordinal_from_scores_matrix_tau(S, tau)
    idx <- which(pattern_key(C) == key_star)
    
    if (length(idx) > 0) {
      U_keep <- rbind(U_keep, U[idx, , drop = FALSE])
    }
    
    if (nrow(U_keep) >= n_draw) break
  }
  
  if (nrow(U_keep) == 0) {
    U_keep <- matrix(rnorm(n_draw * 2), n_draw, 2)
  }
  
  U_keep[sample(seq_len(nrow(U_keep)), n_draw, replace = TRUE), , drop = FALSE]
}

sample_oracle_test_y_2d <- function(X_test,
                                    C_test,
                                    true_params,
                                    sigma_eps,
                                    n_draw = 1000L,
                                    oracle_pool = NULL,
                                    n_pool = 200000L,
                                    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  X_test <- as.matrix(X_test)
  C_test <- as.matrix(C_test)
  
  if (is.null(oracle_pool)) {
    oracle_pool <- make_oracle_pool_2d(true_params, n_pool = n_pool)
  }
  
  n_test <- nrow(C_test)
  out <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  
  keys_test <- pattern_key(C_test)
  
  for (j in seq_len(n_test)) {
    key_j <- keys_test[j]
    idx_pool <- oracle_pool$split_idx[[key_j]]
    
    if (!is.null(idx_pool) && length(idx_pool) >= 5) {
      idx <- sample(idx_pool, n_draw, replace = TRUE)
      U_draw <- oracle_pool$U[idx, , drop = FALSE]
    } else {
      U_draw <- sample_oracle_u_rejection(
        c_star = C_test[j, ],
        true_params = true_params,
        n_draw = n_draw
      )
    }
    
    X_rep <- matrix(rep(X_test[j, ], each = n_draw), n_draw, 2)
    f_draw <- f0_2d(X_rep, U_draw)
    out[, j] <- f_draw + rnorm(n_draw, 0, sigma_eps)
  }
  
  out
}

############################################################
## Deterministic embedding GP baselines
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

gp_mle_fit <- function(X, y) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  
  n <- nrow(X)
  d <- ncol(X)
  
  Dlist <- lapply(seq_len(d), function(j) {
    pairwise_sqdist(X[, j, drop = FALSE])
  })
  
  nll <- function(par) {
    log_sigma2 <- par[1]
    log_rho <- par[2]
    log_theta <- par[-c(1, 2)]
    
    sigma2 <- exp(log_sigma2)
    rho <- exp(log_rho)
    theta <- exp(log_theta)
    
    Rexp <- matrix(0, n, n)
    for (j in seq_len(d)) {
      Rexp <- Rexp + theta[j] * Dlist[[j]]
    }
    
    R <- exp(-Rexp)
    A_mat <- rho^2 * R + diag(n)
    
    Uchol <- try(safe_chol(A_mat), silent = TRUE)
    if (inherits(Uchol, "try-error")) return(1e20)
    
    Ainv_y <- solve_chol(Uchol, y)
    logdetA <- 2 * sum(log(diag(Uchol)))
    quad <- sum(y * Ainv_y)
    
    0.5 * (n * log(2 * pi * sigma2) + logdetA + quad / sigma2)
  }
  
  init <- c(log(0.05), log(3), rep(log(0.5), d))
  lower <- c(log(1e-5), log(0.05), rep(log(1e-4), d))
  upper <- c(log(5), log(100), rep(log(100), d))
  
  opt <- optim(
    par = init,
    fn = nll,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(maxit = 500)
  )
  
  list(
    par = opt$par,
    value = opt$value,
    convergence = opt$convergence,
    X = X,
    y = y
  )
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
  var_lat <- pmax(var_lat, 1e-10)
  
  if (noisy) var_lat <- var_lat + sigma2
  
  list(mean = mu, var = var_lat)
}

sample_gp_mle_predictive <- function(fit, Xstar, n_draw = 1000) {
  pred <- gp_mle_predict(fit, Xstar = Xstar, noisy = TRUE)
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

make_monotone_scores <- function(a, m) {
  if (m == 1) return(0)
  if (m == 2) return(c(0, 1))
  
  logits <- c(a, 0)
  e <- exp(logits - max(logits))
  inc <- e / sum(e)
  
  c(0, cumsum(inc))
}

gp_mle_fit_learned_embedding_ord <- function(X,
                                             C_ord,
                                             y,
                                             m_vec,
                                             n_starts = 4L) {
  X <- as.matrix(X)
  C_ord <- as.matrix(C_ord)
  y <- as.numeric(y)
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(C_ord)
  
  n_feat <- p + q
  n_embed <- sum(pmax(m_vec - 2L, 0L))
  
  nll <- function(par) {
    log_sigma2 <- par[1]
    log_rho <- par[2]
    log_theta <- par[3:(2 + n_feat)]
    
    sigma2 <- exp(log_sigma2)
    rho <- exp(log_rho)
    theta <- exp(log_theta)
    
    embed_par <- if (n_embed > 0L) {
      par[(3 + n_feat):length(par)]
    } else {
      numeric(0)
    }
    
    Z <- matrix(NA_real_, n, q)
    pos <- 0L
    
    for (j in seq_len(q)) {
      mj <- m_vec[j]
      
      if (mj > 2L) {
        a_j <- embed_par[pos + seq_len(mj - 2L)]
        pos <- pos + (mj - 2L)
      } else {
        a_j <- numeric(0)
      }
      
      z_scores <- make_monotone_scores(a_j, mj)
      Z[, j] <- z_scores[C_ord[, j]]
    }
    
    X_aug <- cbind(X, Z)
    
    Rexp <- matrix(0, n, n)
    for (j in seq_len(n_feat)) {
      Rexp <- Rexp + theta[j] * pairwise_sqdist(X_aug[, j, drop = FALSE])
    }
    
    R <- exp(-Rexp)
    A_mat <- rho^2 * R + diag(n)
    
    Uchol <- try(safe_chol(A_mat), silent = TRUE)
    if (inherits(Uchol, "try-error")) return(1e20)
    
    Ainv_y <- solve_chol(Uchol, y)
    logdetA <- 2 * sum(log(diag(Uchol)))
    quad <- sum(y * Ainv_y)
    
    0.5 * (n * log(2 * pi * sigma2) + logdetA + quad / sigma2)
  }
  
  lower <- c(
    log(1e-5),
    log(0.05),
    rep(log(1e-4), n_feat),
    rep(-6, n_embed)
  )
  
  upper <- c(
    log(5),
    log(100),
    rep(log(100), n_feat),
    rep(6, n_embed)
  )
  
  make_init <- function() {
    c(
      log(0.05),
      log(3),
      rep(log(0.5), n_feat),
      rnorm(n_embed, 0, 0.5)
    )
  }
  
  starts <- replicate(n_starts, make_init(), simplify = FALSE)
  
  opts <- lapply(starts, function(init) {
    optim(
      par = init,
      fn = nll,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = 500)
    )
  })
  
  best_id <- which.min(vapply(opts, function(z) z$value, numeric(1)))
  opt <- opts[[best_id]]
  
  embed_par <- if (n_embed > 0L) {
    opt$par[(3 + n_feat):length(opt$par)]
  } else {
    numeric(0)
  }
  
  z_scores_list <- vector("list", q)
  pos <- 0L
  
  for (j in seq_len(q)) {
    mj <- m_vec[j]
    
    if (mj > 2L) {
      a_j <- embed_par[pos + seq_len(mj - 2L)]
      pos <- pos + (mj - 2L)
    } else {
      a_j <- numeric(0)
    }
    
    z_scores_list[[j]] <- make_monotone_scores(a_j, mj)
  }
  
  Z <- matrix(NA_real_, n, q)
  for (j in seq_len(q)) {
    Z[, j] <- z_scores_list[[j]][C_ord[, j]]
  }
  
  list(
    par = opt$par[1:(2 + n_feat)],
    embed_par = embed_par,
    z_scores_list = z_scores_list,
    value = opt$value,
    convergence = opt$convergence,
    X_aug = cbind(X, Z),
    X = X,
    C_ord = C_ord,
    y = y,
    m_vec = m_vec
  )
}

gp_mle_predict_learned_embedding_ord <- function(fit,
                                                 X_star,
                                                 C_star,
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

fit_embedding_baselines_ord <- function(X_raw,
                                        y_raw,
                                        C_ord,
                                        m_vec,
                                        n_starts_learned = 4L) {
  X_raw <- as.matrix(X_raw)
  C_ord <- as.matrix(C_ord)
  
  X_center <- colMeans(X_raw)
  X_scale <- apply(X_raw, 2, sd)
  X_scale[!is.finite(X_scale) | X_scale <= 0] <- 1
  
  X <- sweep(sweep(X_raw, 2, X_center, "-"), 2, X_scale, "/")
  
  y_center <- mean(y_raw)
  y_scale <- sd(y_raw)
  if (!is.finite(y_scale) || y_scale <= 0) y_scale <- 1
  y <- as.numeric((y_raw - y_center) / y_scale)
  
  q <- ncol(C_ord)
  
  Z_gauss <- matrix(NA_real_, nrow(C_ord), q)
  for (j in seq_len(q)) {
    Z_gauss[, j] <- qnorm(C_ord[, j] / (m_vec[j] + 1))
  }
  
  fit_gauss <- gp_mle_fit(
    X = cbind(X, Z_gauss),
    y = y
  )
  
  z_cm_scores <- lapply(seq_len(q), function(j) {
    conditional_mean_scores_ord(C_ord[, j], m_vec[j])
  })
  
  Z_cm <- matrix(NA_real_, nrow(C_ord), q)
  for (j in seq_len(q)) {
    Z_cm[, j] <- z_cm_scores[[j]][C_ord[, j]]
  }
  
  fit_cm <- gp_mle_fit(
    X = cbind(X, Z_cm),
    y = y
  )
  
  fit_learned <- gp_mle_fit_learned_embedding_ord(
    X = X,
    C_ord = C_ord,
    y = y,
    m_vec = m_vec,
    n_starts = n_starts_learned
  )
  
  list(
    X_center = X_center,
    X_scale = X_scale,
    y_center = y_center,
    y_scale = y_scale,
    fit_gauss = fit_gauss,
    fit_cm = fit_cm,
    fit_learned = fit_learned,
    z_cm_scores = z_cm_scores,
    m_vec = m_vec
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

############################################################
## Predictive scoring
############################################################

crps_sample_one <- function(draws, y) {
  draws <- sort(as.numeric(draws))
  S <- length(draws)
  
  term1 <- mean(abs(draws - y))
  weights <- 2 * seq_len(S) - S - 1
  mean_abs_pair <- 2 * sum(weights * draws) / S^2
  
  term1 - 0.5 * mean_abs_pair
}

crps_sample_matrix <- function(draw_mat, y) {
  vapply(seq_along(y), function(j) {
    crps_sample_one(draw_mat[, j], y[j])
  }, numeric(1))
}

interval_score <- function(lo, hi, y, alpha = 0.05) {
  width <- hi - lo
  width +
    2 / alpha * (lo - y) * (y < lo) +
    2 / alpha * (y - hi) * (y > hi)
}

summarize_predictive_samples <- function(draw_mat,
                                         y_true,
                                         method,
                                         rep_id,
                                         n_calib,
                                         scenario = "study2") {
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
    IntervalScore95 = mean(interval_score(lo, hi, y_true)),
    stringsAsFactors = FALSE
  )
}

############################################################
## Real-data wrappers and latent summaries
############################################################

standardize_U_obs_for_fit <- function(U_obs,
                                      calib_idx = NULL,
                                      standardize = TRUE) {
  if (is.null(U_obs)) {
    return(list(
      U_std = NULL,
      calib_idx = integer(0),
      center = NULL,
      scale = NULL
    ))
  }
  
  U_obs <- as.matrix(U_obs)
  n <- nrow(U_obs)
  d <- ncol(U_obs)
  
  if (is.null(calib_idx) || length(calib_idx) == 0L) {
    calib_idx <- which(stats::complete.cases(U_obs))
  } else {
    calib_idx <- sort(as.integer(calib_idx))
  }
  
  if (length(calib_idx) == 0L) {
    return(list(
      U_std = U_obs,
      calib_idx = integer(0),
      center = rep(0, d),
      scale = rep(1, d)
    ))
  }
  
  center <- rep(0, d)
  scale <- rep(1, d)
  
  if (isTRUE(standardize)) {
    center <- colMeans(U_obs[calib_idx, , drop = FALSE], na.rm = TRUE)
    scale <- apply(U_obs[calib_idx, , drop = FALSE], 2, sd, na.rm = TRUE)
    scale[!is.finite(scale) | scale <= 0] <- 1
  }
  
  U_std <- U_obs
  for (j in seq_len(d)) {
    U_std[, j] <- (U_obs[, j] - center[j]) / scale[j]
  }
  
  list(
    U_std = U_std,
    calib_idx = calib_idx,
    center = center,
    scale = scale
  )
}

fit_eivgp_ordprobit_exact_realdata <- function(X_raw,
                                               y_raw,
                                               C_ord,
                                               U_obs = NULL,
                                               calib_idx = NULL,
                                               d = NULL,
                                               m_vec = NULL,
                                               ident = c("lower_triangular", "none"),
                                               standardize_U = TRUE,
                                               n_iter = 6000L,
                                               burn = 2000L,
                                               thin = 1L,
                                               n_chains = 8L,
                                               preset = "balanced",
                                               seed = 1L,
                                               parallel_chains = FALSE,
                                               verbose = TRUE) {
  ident <- match.arg(ident)
  
  X_raw <- as.matrix(X_raw)
  C_ord <- as.matrix(C_ord)
  
  if (is.null(d)) {
    if (!is.null(U_obs)) {
      d <- ncol(as.matrix(U_obs))
    } else {
      stop("Please supply d when U_obs is NULL.")
    }
  }
  
  U_std_obj <- standardize_U_obs_for_fit(
    U_obs = U_obs,
    calib_idx = calib_idx,
    standardize = standardize_U
  )
  
  fit <- fit_eivgp_ordprobit_fb(
    X_raw = X_raw,
    y_raw = y_raw,
    C_ord = C_ord,
    U_obs = U_std_obj$U_std,
    calib_idx = U_std_obj$calib_idx,
    U_true_eval = NULL,
    d = d,
    m_vec = m_vec,
    ident = ident,
    n_iter = n_iter,
    burn = burn,
    thin = thin,
    n_chains = n_chains,
    preset = preset,
    seed = seed,
    parallel_chains = parallel_chains,
    verbose = verbose
  )
  
  fit$u_standardization <- list(
    standardized = isTRUE(standardize_U),
    center = U_std_obj$center,
    scale = U_std_obj$scale
  )
  fit$data$U_obs_raw <- if (is.null(U_obs)) NULL else as.matrix(U_obs)
  
  fit
}

predict_eivgp_ordprobit_exact <- function(fit_obj,
                                          X_new_raw,
                                          C_new,
                                          n_draw = 1000L,
                                          n_per_draw = 1L,
                                          n_new_latent_gibbs = 20L,
                                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  n_saved <- dim(fit_obj$mcmc$samples_U)[1]
  draw_ids <- seq_len(n_saved)
  
  if (length(draw_ids) > n_draw) {
    draw_ids <- sample(draw_ids, n_draw)
  }
  
  sample_eiv_test_y_ordprobit_fb(
    X_test_raw = X_new_raw,
    C_test = C_new,
    fit_obj = fit_obj,
    draw_ids = draw_ids,
    n_per_draw = n_per_draw,
    n_new_latent_gibbs = n_new_latent_gibbs
  )
}

posterior_U_summary_eivgp <- function(fit_obj,
                                      probs = c(0.025, 0.50, 0.975),
                                      original_scale = FALSE,
                                      true_U = NULL) {
  samples_U <- fit_obj$mcmc$samples_U
  
  n_draw <- dim(samples_U)[1]
  n <- dim(samples_U)[2]
  d <- dim(samples_U)[3]
  
  arr <- samples_U
  
  if (isTRUE(original_scale) &&
      !is.null(fit_obj$u_standardization) &&
      !is.null(fit_obj$u_standardization$center)) {
    cc <- fit_obj$u_standardization$center
    ss <- fit_obj$u_standardization$scale
    
    for (k in seq_len(d)) {
      arr[, , k] <- cc[k] + ss[k] * arr[, , k]
    }
  }
  
  out <- dplyr::bind_rows(
    lapply(seq_len(d), function(k) {
      mat <- arr[, , k, drop = FALSE]
      mat <- matrix(mat, nrow = n_draw, ncol = n)
      
      qs <- apply(
        mat,
        2,
        stats::quantile,
        probs = probs,
        na.rm = TRUE,
        names = FALSE
      )
      
      tmp <- data.frame(
        id = seq_len(n),
        coord = paste0("u", k),
        mean = colMeans(mat, na.rm = TRUE),
        sd = apply(mat, 2, stats::sd, na.rm = TRUE),
        q025 = qs[1, ],
        q500 = qs[2, ],
        q975 = qs[3, ],
        calibrated = seq_len(n) %in% fit_obj$data$calib_idx
      )
      
      if (!is.null(true_U)) {
        true_U <- as.matrix(true_U)
        tmp$true_u <- true_U[, k]
        tmp$error <- tmp$mean - tmp$true_u
        tmp$abs_error <- abs(tmp$error)
        tmp$covered95 <- tmp$true_u >= tmp$q025 & tmp$true_u <= tmp$q975
      }
      
      tmp
    })
  )
  
  out
}

summarize_U_imputation_metrics <- function(U_summary) {
  if (!all(c("true_u", "error", "abs_error", "covered95") %in% names(U_summary))) {
    stop("U_summary must include true latent values.")
  }
  
  U_summary |>
    dplyr::filter(!calibrated) |>
    dplyr::group_by(coord) |>
    dplyr::summarise(
      n = dplyr::n(),
      bias = mean(error, na.rm = TRUE),
      rmse = sqrt(mean(error^2, na.rm = TRUE)),
      mae = mean(abs_error, na.rm = TRUE),
      coverage95 = mean(covered95, na.rm = TRUE),
      mean_width95 = mean(q975 - q025, na.rm = TRUE),
      .groups = "drop"
    )
}

raw_y_from_fit <- function(fit, y = fit$data$y) {
  fit$data$y_center + fit$data$y_scale * y
}

make_latent_surface_recovery_2d <- function(fit,
                                            truth_fun = f0_2d,
                                            x_ref_raw = NULL,
                                            u_lim = c(-2.2, 2.2),
                                            grid_size = 45L,
                                            max_draw = 200L) {
  if (fit$data$d != 2L) {
    stop("make_latent_surface_recovery_2d currently requires d = 2.")
  }
  
  if (is.null(x_ref_raw)) {
    x_ref_raw <- rep(0, fit$data$p)
  }
  
  u_grid <- seq(u_lim[1], u_lim[2], length.out = grid_size)
  grid <- expand.grid(u1 = u_grid, u2 = u_grid)
  
  U_star <- as.matrix(grid[, c("u1", "u2")])
  
  X_star_raw <- matrix(
    rep(x_ref_raw, each = nrow(grid)),
    nrow = nrow(grid),
    ncol = fit$data$p
  )
  
  X_star <- sweep(
    sweep(X_star_raw, 2, fit$data$X_center, "-"),
    2,
    fit$data$X_scale,
    "/"
  )
  
  draw_ids <- seq_len(dim(fit$mcmc$samples_U)[1])
  
  if (length(draw_ids) > max_draw) {
    draw_ids <- sample(draw_ids, max_draw)
  }
  
  f_draws <- matrix(NA_real_, nrow = length(draw_ids), ncol = nrow(grid))
  
  for (ii in seq_along(draw_ids)) {
    s <- draw_ids[ii]
    
    U_train_s <- as.matrix(fit$mcmc$samples_U[s, , ])
    
    pred <- gp_predict_draw_general(
      X_train = fit$data$X,
      U_train = U_train_s,
      y_train = fit$data$y,
      X_star = X_star,
      U_star = U_star,
      logtheta = fit$mcmc$samples_logtheta[s, ],
      sigma2_eps = fit$mcmc$samples_sigma2[s],
      noisy = FALSE
    )
    
    f_draws[ii, ] <- raw_y_from_fit(fit, pred$mean)
  }
  
  out <- grid
  out$mean <- colMeans(f_draws, na.rm = TRUE)
  out$lo <- apply(f_draws, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
  out$hi <- apply(f_draws, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
  out$width95 <- out$hi - out$lo
  out$truth <- truth_fun(X_star_raw, U_star)
  out$error <- out$mean - out$truth
  out$covered95 <- out$truth >= out$lo & out$truth <= out$hi
  
  attr(out, "surface_rmse") <- sqrt(mean(out$error^2, na.rm = TRUE))
  attr(out, "surface_mae") <- mean(abs(out$error), na.rm = TRUE)
  attr(out, "surface_coverage95") <- mean(out$covered95, na.rm = TRUE)
  attr(out, "surface_mean_width95") <- mean(out$width95, na.rm = TRUE)
  
  out
}

summarize_latent_surface_recovery_2d <- function(surface_df) {
  data.frame(
    surface_rmse = sqrt(mean(surface_df$error^2, na.rm = TRUE)),
    surface_mae = mean(abs(surface_df$error), na.rm = TRUE),
    surface_coverage95 = mean(surface_df$covered95, na.rm = TRUE),
    surface_mean_width95 = mean(surface_df$width95, na.rm = TRUE)
  )
}