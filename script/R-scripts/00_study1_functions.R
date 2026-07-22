############################################################
## 00_study1_functions.R
##
## Functions for revised Study I:
## one-dimensional deterministic-threshold EIV-GP example.
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))

load_or_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
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

maximin_lhs_1d <- function(n, lower = -2, upper = 2) {
  z <- (seq_len(n) - runif(n)) / n
  z <- sample(z)
  lower + (upper - lower) * z
}

make_class <- function(u, tau) {
  as.integer(cut(u, breaks = c(-Inf, tau, Inf), labels = FALSE))
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

############################################################
## Slice samplers
############################################################

bounded_slice_update <- function(x0, logf, w = 1,
                                 lower = -Inf, upper = Inf,
                                 max_steps_out = 50,
                                 max_iter = 200) {
  if (upper <= lower) {
    return(list(x = x0, n_eval = 0))
  }
  
  eps <- 1e-12
  if (is.finite(lower)) x0 <- max(x0, lower + eps)
  if (is.finite(upper)) x0 <- min(x0, upper - eps)
  
  f0 <- logf(x0)
  n_eval <- 1L
  
  if (!is.finite(f0)) {
    return(list(x = x0, n_eval = n_eval))
  }
  
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

full_interval_slice_update <- function(x0, logf, lower, upper,
                                       max_iter = 200) {
  if (!is.finite(lower) || !is.finite(upper)) {
    stop("full_interval_slice_update requires finite lower and upper bounds.")
  }
  
  if (upper <= lower) {
    return(list(x = x0, n_eval = 0))
  }
  
  eps <- 1e-12
  lower_eps <- lower + eps
  upper_eps <- upper - eps
  
  if (upper_eps <= lower_eps) {
    return(list(x = x0, n_eval = 0))
  }
  
  x0 <- min(max(x0, lower_eps), upper_eps)
  
  f0 <- logf(x0)
  n_eval <- 1L
  
  if (!is.finite(f0)) {
    return(list(x = x0, n_eval = n_eval))
  }
  
  logy <- f0 + log(runif(1))
  L <- lower_eps
  R <- upper_eps
  
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
## Data-generating mechanism
############################################################

f0_1d <- function(x, u, scenario = c("active", "inactive")) {
  scenario <- match.arg(scenario)
  
  if (scenario == "inactive") {
    return(cos(0.9 * pi * u))
  }
  
  0.5 * sin(pi * x / 2) +
    cos(0.9 * pi * u) +
    0.35 * (x / 2) * sin(0.9 * pi * u)
}

simulate_1d_data <- function(n = 100,
                             n_test = 500,
                             m = 6,
                             scenario = c("active", "inactive"),
                             sigma_eps = 0.1,
                             seed = NULL) {
  scenario <- match.arg(scenario)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  tau_true <- c(-5 / 3, -10 / 9, 0, 10 / 9, 5 / 3)
  
  repeat {
    x <- maximin_lhs_1d(n, lower = -2, upper = 2)
    u <- rnorm(n)
    c_ord <- make_class(u, tau_true)
    
    if (all(tabulate(c_ord, nbins = m) >= 3)) break
  }
  
  f <- f0_1d(x, u, scenario = scenario)
  y <- f + rnorm(n, mean = 0, sd = sigma_eps)
  
  x_test <- runif(n_test, -2, 2)
  u_test <- rnorm(n_test)
  c_test <- make_class(u_test, tau_true)
  f_test <- f0_1d(x_test, u_test, scenario = scenario)
  y_test <- f_test + rnorm(n_test, mean = 0, sd = sigma_eps)
  
  list(
    train = data.frame(
      x = x,
      u = u,
      c = c_ord,
      y = y,
      f = f
    ),
    test = data.frame(
      x = x_test,
      u = u_test,
      c = c_test,
      y = y_test,
      f = f_test
    ),
    tau_true = tau_true,
    sigma_eps = sigma_eps,
    scenario = scenario,
    m = m
  )
}

make_nested_calibration_sets <- function(n, calib_grid, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  ord <- sample(seq_len(n))
  
  out <- lapply(calib_grid, function(k) {
    if (k == 0) integer(0) else sort(ord[seq_len(k)])
  })
  
  names(out) <- as.character(calib_grid)
  out
}

############################################################
## GP likelihood and prediction
############################################################

a_eps0 <- 2
b_eps0 <- 0.05

logtheta_prior_mean <- c(log(3), log(0.5), log(0.5))
logtheta_prior_sd <- c(1.5, 1.5, 1.5)

theta_log_bounds <- rbind(
  log_rho     = c(log(0.05), log(100)),
  log_theta_x = c(log(1e-4), log(100)),
  log_theta_u = c(log(1e-4), log(100))
)

gp_state_1d <- function(y, u, Dx, logtheta, sigma2_eps) {
  n <- length(y)
  
  rho <- exp(logtheta[1])
  theta_x <- exp(logtheta[2])
  theta_u <- exp(logtheta[3])
  
  Du <- pairwise_sqdist(matrix(u, ncol = 1))
  R <- exp(-theta_x * Dx - theta_u * Du)
  
  A <- rho^2 * R + diag(n)
  
  U <- safe_chol(A)
  Ainv_y <- solve_chol(U, y)
  
  logdetA <- 2 * sum(log(diag(U)))
  quad <- sum(y * Ainv_y)
  
  loglik <- -0.5 * (
    n * log(2 * pi * sigma2_eps) +
      logdetA +
      quad / sigma2_eps
  )
  
  list(
    loglik = loglik,
    R = R,
    A = A,
    cholA = U,
    Ainv_y = Ainv_y,
    logdetA = logdetA,
    quad = quad
  )
}

gp_predict_draw <- function(x_train, u_train, y_train,
                            x_star, u_star,
                            logtheta, sigma2_eps,
                            noisy = FALSE) {
  n <- length(y_train)
  
  rho <- exp(logtheta[1])
  theta_x <- exp(logtheta[2])
  theta_u <- exp(logtheta[3])
  
  Dx_train <- pairwise_sqdist(matrix(x_train, ncol = 1))
  Du_train <- pairwise_sqdist(matrix(u_train, ncol = 1))
  
  R <- exp(-theta_x * Dx_train - theta_u * Du_train)
  K <- rho^2 * sigma2_eps * R
  C <- K + sigma2_eps * diag(n)
  
  U <- safe_chol(C)
  alpha <- solve_chol(U, y_train)
  
  Dxs <- pairwise_sqdist(
    matrix(x_star, ncol = 1),
    matrix(x_train, ncol = 1)
  )
  Dus <- pairwise_sqdist(
    matrix(u_star, ncol = 1),
    matrix(u_train, ncol = 1)
  )
  
  R_star <- exp(-theta_x * Dxs - theta_u * Dus)
  K_star <- rho^2 * sigma2_eps * R_star
  
  mu <- as.numeric(K_star %*% alpha)
  
  v <- forwardsolve(t(U), t(K_star))
  var_lat <- rho^2 * sigma2_eps - colSums(v^2)
  var_lat <- pmax(var_lat, 1e-10)
  
  if (noisy) {
    var_lat <- var_lat + sigma2_eps
  }
  
  list(mean = mu, var = var_lat)
}

############################################################
## EIV-GP sampler
############################################################

make_default_control <- function(n, n_mis,
                                 preset = c("fast", "balanced", "thorough")) {
  preset <- match.arg(preset)
  
  if (preset == "fast") {
    local_frac <- 0.03
    theta_update_every <- 10L
    block_ess_every <- 2L
    n_blocks_per_iter <- 1L
    global_ess_every <- 50L
    full_local_every <- 100L
  }
  
  if (preset == "balanced") {
    local_frac <- 0.06
    theta_update_every <- 5L
    block_ess_every <- 1L
    n_blocks_per_iter <- 1L
    global_ess_every <- 25L
    full_local_every <- 50L
  }
  
  if (preset == "thorough") {
    local_frac <- 0.10
    theta_update_every <- 3L
    block_ess_every <- 1L
    n_blocks_per_iter <- 2L
    global_ess_every <- 10L
    full_local_every <- 25L
  }
  
  local_per_iter <- min(
    n_mis,
    max(3L, ceiling(local_frac * max(n_mis, 1L)))
  )
  
  ess_block_size <- min(
    n_mis,
    max(5L, ceiling(sqrt(max(n_mis, 1L))))
  )
  
  list(
    preset = preset,
    local_per_iter = local_per_iter,
    full_local_every = full_local_every,
    ess_block_size = ess_block_size,
    block_ess_every = block_ess_every,
    n_blocks_per_iter = n_blocks_per_iter,
    global_ess_every = global_ess_every,
    theta_update_every = theta_update_every,
    theta_slice_width_init = c(1.0, 1.0, 1.0),
    adapt_theta_width = TRUE,
    adapt_every = 100L,
    adapt_window = 500L,
    theta_width_min = 0.20,
    theta_width_max = 2.50
  )
}

check_constraints_1d <- function(u, c_ord, tau) {
  lower <- c(-Inf, tau)[c_ord]
  upper <- c(tau, Inf)[c_ord]
  
  all(u > lower & u <= upper)
}

log_prior_logtheta <- function(logtheta) {
  if (any(logtheta < theta_log_bounds[, 1]) ||
      any(logtheta > theta_log_bounds[, 2])) {
    return(-Inf)
  }
  
  sum(dnorm(
    logtheta,
    mean = logtheta_prior_mean,
    sd = logtheta_prior_sd,
    log = TRUE
  ))
}

gp_loglik_with_constraints_1d <- function(y, u, Dx, c_ord, tau,
                                          logtheta, sigma2_eps) {
  if (!check_constraints_1d(u, c_ord, tau)) {
    return(-Inf)
  }
  
  gp_state_1d(y, u, Dx, logtheta, sigma2_eps)$loglik
}

theta_logpost_1d <- function(y, u, Dx, logtheta, sigma2_eps) {
  lp <- log_prior_logtheta(logtheta)
  
  if (!is.finite(lp)) {
    return(-Inf)
  }
  
  gp_state_1d(y, u, Dx, logtheta, sigma2_eps)$loglik + lp
}

initialize_tau_1d <- function(c_ord, u_obs, calib_idx, m, tau_bound = 8) {
  n <- length(c_ord)
  
  counts <- tabulate(c_ord, nbins = m)
  probs <- cumsum(counts)[1:(m - 1)] / n
  probs <- pmin(pmax(probs, 0.03), 0.97)
  
  tau <- qnorm(probs)
  tau <- pmin(pmax(tau, -tau_bound + 1e-3), tau_bound - 1e-3)
  
  if (length(calib_idx) > 0) {
    obs_c <- c_ord[calib_idx]
    obs_u <- u_obs[calib_idx]
    
    eps_ord <- 1e-4
    
    for (j in seq_len(m - 1)) {
      lower <- max(c(-tau_bound, obs_u[obs_c <= j]), na.rm = TRUE)
      upper <- min(c( tau_bound, obs_u[obs_c >  j]), na.rm = TRUE)
      
      if (lower >= upper) {
        stop("Calibrated u values are incompatible with ordinal labels.")
      }
      
      tau[j] <- min(max(tau[j], lower + eps_ord), upper - eps_ord)
    }
  }
  
  tau
}

update_tau_1d <- function(tau, u, c_ord, m, tau_bound = 8) {
  tau_new <- tau
  
  for (j in seq_len(m - 1)) {
    L <- max(
      c(
        -tau_bound,
        if (j > 1) tau_new[j - 1] else -Inf,
        u[c_ord <= j]
      ),
      na.rm = TRUE
    )
    
    U <- min(
      c(
        tau_bound,
        if (j < m - 1) tau_new[j + 1] else Inf,
        u[c_ord > j]
      ),
      na.rm = TRUE
    )
    
    if (is.finite(L) && is.finite(U) && L < U) {
      tau_new[j] <- runif(1, L, U)
    }
  }
  
  tau_new
}

sample_sigma2_eps_1d <- function(y, u, Dx, logtheta) {
  n <- length(y)
  
  rho <- exp(logtheta[1])
  theta_x <- exp(logtheta[2])
  theta_u <- exp(logtheta[3])
  
  Du <- pairwise_sqdist(matrix(u, ncol = 1))
  R <- exp(-theta_x * Dx - theta_u * Du)
  
  A <- rho^2 * R + diag(n)
  
  U <- safe_chol(A)
  Ainv_y <- solve_chol(U, y)
  quad <- sum(y * Ainv_y)
  
  shape <- a_eps0 + n / 2
  rate <- b_eps0 + 0.5 * quad
  
  1 / rgamma(1, shape = shape, rate = rate)
}

update_u_ess_block_1d <- function(y, u, Dx, c_ord, tau, logtheta, sigma2_eps,
                                  block_idx, max_try = 300L) {
  if (length(block_idx) == 0) {
    return(list(u = u, n_eval = 0L, accepted = TRUE))
  }
  
  u_block <- u[block_idx]
  
  loglik_fun <- function(u_block_prop) {
    u_prop <- u
    u_prop[block_idx] <- u_block_prop
    
    gp_loglik_with_constraints_1d(
      y = y,
      u = u_prop,
      Dx = Dx,
      c_ord = c_ord,
      tau = tau,
      logtheta = logtheta,
      sigma2_eps = sigma2_eps
    )
  }
  
  loglik_cur <- loglik_fun(u_block)
  
  if (!is.finite(loglik_cur)) {
    stop("Current latent state violates constraints or has invalid likelihood.")
  }
  
  nu <- rnorm(length(u_block))
  
  logy <- loglik_cur + log(runif(1))
  
  angle <- runif(1, 0, 2 * pi)
  angle_min <- angle - 2 * pi
  angle_max <- angle
  
  n_eval <- 1L
  
  for (try_id in seq_len(max_try)) {
    u_block_prop <- u_block * cos(angle) + nu * sin(angle)
    
    loglik_prop <- loglik_fun(u_block_prop)
    n_eval <- n_eval + 1L
    
    if (is.finite(loglik_prop) && loglik_prop >= logy) {
      u_new <- u
      u_new[block_idx] <- u_block_prop
      return(list(u = u_new, n_eval = n_eval, accepted = TRUE))
    }
    
    if (angle < 0) {
      angle_min <- angle
    } else {
      angle_max <- angle
    }
    
    angle <- runif(1, angle_min, angle_max)
  }
  
  list(u = u, n_eval = n_eval, accepted = FALSE)
}

update_u_local_z_slice_1d <- function(y, u, Dx, c_ord, tau,
                                      logtheta, sigma2_eps,
                                      update_idx) {
  n_eval_total <- 0L
  eps_z <- 1e-12
  
  update_idx <- sample(update_idx)
  
  for (i in update_idx) {
    lower_u <- c(-Inf, tau)[c_ord[i]]
    upper_u <- c(tau, Inf)[c_ord[i]]
    
    z_lower <- pnorm(lower_u)
    z_upper <- pnorm(upper_u)
    
    z_lower <- max(z_lower, eps_z)
    z_upper <- min(z_upper, 1 - eps_z)
    
    if (z_upper <= z_lower) next
    
    z0 <- pnorm(u[i])
    z0 <- min(max(z0, z_lower + eps_z), z_upper - eps_z)
    
    logf_z <- function(z) {
      if (z <= z_lower || z >= z_upper) {
        return(-Inf)
      }
      
      u_prop_i <- qnorm(z)
      
      if (!(u_prop_i > lower_u && u_prop_i <= upper_u)) {
        return(-Inf)
      }
      
      u_prop <- u
      u_prop[i] <- u_prop_i
      
      gp_state_1d(y, u_prop, Dx, logtheta, sigma2_eps)$loglik
    }
    
    ans <- full_interval_slice_update(
      x0 = z0,
      logf = logf_z,
      lower = z_lower,
      upper = z_upper,
      max_iter = 200L
    )
    
    n_eval_total <- n_eval_total + ans$n_eval
    u[i] <- qnorm(ans$x)
  }
  
  list(u = u, n_eval = n_eval_total)
}

update_logtheta_slice_1d <- function(y, u, Dx, logtheta, sigma2_eps,
                                     theta_slice_width) {
  n_eval_total <- 0L
  
  for (j in seq_along(logtheta)) {
    logf_j <- function(val) {
      lt <- logtheta
      lt[j] <- val
      
      theta_logpost_1d(y, u, Dx, lt, sigma2_eps)
    }
    
    ans <- bounded_slice_update(
      x0 = logtheta[j],
      logf = logf_j,
      w = theta_slice_width[j],
      lower = theta_log_bounds[j, 1],
      upper = theta_log_bounds[j, 2],
      max_steps_out = 30L,
      max_iter = 100L
    )
    
    logtheta[j] <- ans$x
    n_eval_total <- n_eval_total + ans$n_eval
  }
  
  list(logtheta = logtheta, n_eval = n_eval_total)
}

############################################################
## MCMC diagnostics
############################################################

split_rhat <- function(chain_list) {
  chain_list <- lapply(chain_list, function(x) as.numeric(x[is.finite(x)]))
  
  lens <- vapply(chain_list, length, integer(1))
  n0 <- min(lens)
  
  if (length(chain_list) < 2 || n0 < 20) {
    return(NA_real_)
  }
  
  chain_list <- lapply(chain_list, function(x) tail(x, n0))
  
  n_half <- floor(n0 / 2)
  
  if (n_half < 10) {
    return(NA_real_)
  }
  
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
  
  if (!is.finite(W) || W <= 0) {
    return(NA_real_)
  }
  
  var_hat <- ((n_split - 1) / n_split) * W + B / n_split
  
  sqrt(var_hat / W)
}

ess_ips <- function(x, max_lag = NULL) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  n <- length(x)
  
  if (n < 5) {
    return(NA_real_)
  }
  
  sx <- sd(x)
  
  if (!is.finite(sx)) {
    return(NA_real_)
  }
  
  if (sx == 0) {
    return(n)
  }
  
  if (is.null(max_lag)) {
    max_lag <- min(n - 1L, 1000L)
  } else {
    max_lag <- min(n - 1L, as.integer(max_lag))
  }
  
  if (max_lag < 1L) {
    return(n)
  }
  
  ac <- tryCatch(
    as.numeric(stats::acf(
      x,
      lag.max = max_lag,
      plot = FALSE,
      demean = TRUE
    )$acf),
    error = function(e) NA_real_
  )
  
  if (length(ac) <= 1 || all(!is.finite(ac))) {
    return(NA_real_)
  }
  
  ac <- ac[-1]
  ac <- ac[is.finite(ac)]
  
  if (length(ac) == 0) {
    return(n)
  }
  
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

############################################################
## Main EIV-GP fitting function
############################################################

fit_eivgp_1d <- function(x_raw,
                         y_raw,
                         c_ord,
                         u_true = NULL,
                         calib_idx = integer(0),
                         m = 6L,
                         tau_true = NULL,
                         n_iter = 6000L,
                         burn = 2000L,
                         thin = 1L,
                         n_chains = 4L,
                         preset = "balanced",
                         seed = 1L,
                         parallel_chains = FALSE,
                         verbose = FALSE) {
  set.seed(seed)
  
  n <- length(y_raw)
  
  calib_idx <- sort(as.integer(calib_idx))
  miss_idx <- setdiff(seq_len(n), calib_idx)
  
  if (length(calib_idx) > 0 && is.null(u_true)) {
    stop("u_true must be provided when calibration observations are used.")
  }
  
  x_center <- mean(x_raw)
  x_scale <- sd(x_raw)
  x <- as.numeric((x_raw - x_center) / x_scale)
  
  y_center <- mean(y_raw)
  y_scale <- sd(y_raw)
  y <- as.numeric((y_raw - y_center) / y_scale)
  
  u_obs <- rep(NA_real_, n)
  
  if (length(calib_idx) > 0) {
    u_obs[calib_idx] <- u_true[calib_idx]
  }
  
  tau_bound <- 8
  
  control <- make_default_control(
    n = n,
    n_mis = length(miss_idx),
    preset = preset
  )
  
  Dx <- pairwise_sqdist(matrix(x, ncol = 1))
  n_save <- floor((n_iter - burn) / thin)
  
  initialize_chain_state <- function(chain_seed) {
    set.seed(chain_seed)
    
    tau0 <- initialize_tau_1d(
      c_ord = c_ord,
      u_obs = u_obs,
      calib_idx = calib_idx,
      m = m,
      tau_bound = tau_bound
    )
    
    u0 <- u_obs
    
    if (length(miss_idx) > 0) {
      lower_all0 <- c(-Inf, tau0)[c_ord]
      upper_all0 <- c(tau0, Inf)[c_ord]
      
      u0[miss_idx] <- rtruncnorm_vec(
        mean = rep(0, length(miss_idx)),
        sd = 1,
        lower = lower_all0[miss_idx],
        upper = upper_all0[miss_idx]
      )
    }
    
    if (length(miss_idx) > 0) {
      for (rr in seq_len(3)) {
        tau0 <- update_tau_1d(tau0, u0, c_ord, m, tau_bound)
        
        lower_all0 <- c(-Inf, tau0)[c_ord]
        upper_all0 <- c(tau0, Inf)[c_ord]
        
        u0[miss_idx] <- rtruncnorm_vec(
          mean = rep(0, length(miss_idx)),
          sd = 1,
          lower = lower_all0[miss_idx],
          upper = upper_all0[miss_idx]
        )
      }
    }
    
    sigma2_eps0 <- exp(log(0.05) + rnorm(1, 0, 0.5))
    sigma2_eps0 <- min(max(sigma2_eps0, 1e-4), 2)
    
    Du_init <- pairwise_sqdist(matrix(u0, ncol = 1))
    
    med_dx <- median(Dx[upper.tri(Dx)], na.rm = TRUE)
    med_du <- median(Du_init[upper.tri(Du_init)], na.rm = TRUE)
    
    theta_x_init <- ifelse(
      is.finite(med_dx) && med_dx > 0,
      1 / med_dx,
      0.5
    )
    
    theta_u_init <- ifelse(
      is.finite(med_du) && med_du > 0,
      1 / med_du,
      0.5
    )
    
    rho_init <- sqrt(max(var(y) / sigma2_eps0 - 1, 1))
    
    logtheta0 <- c(log(rho_init), log(theta_x_init), log(theta_u_init))
    logtheta0 <- logtheta0 + rnorm(3, mean = 0, sd = c(0.4, 0.7, 0.7))
    logtheta0 <- pmin(pmax(logtheta0, theta_log_bounds[, 1]), theta_log_bounds[, 2])
    
    list(
      u_curr = u0,
      tau = tau0,
      sigma2_eps = sigma2_eps0,
      logtheta = logtheta0,
      theta_slice_width = control$theta_slice_width_init
    )
  }
  
  run_one_chain <- function(chain_id, chain_seed) {
    set.seed(chain_seed)
    
    state <- initialize_chain_state(chain_seed)
    
    u_curr <- state$u_curr
    tau <- state$tau
    sigma2_eps <- state$sigma2_eps
    logtheta <- state$logtheta
    theta_slice_width <- state$theta_slice_width
    
    samples_u_chain <- matrix(NA_real_, n_save, n)
    samples_tau_chain <- matrix(NA_real_, n_save, m - 1)
    samples_logtheta_chain <- matrix(NA_real_, n_save, 3)
    samples_sigma2_chain <- numeric(n_save)
    
    logtheta_trace_all <- matrix(NA_real_, n_iter, 3)
    
    block_ess_eval_total <- 0L
    block_ess_accept_total <- 0L
    block_ess_total <- 0L
    
    global_ess_eval_total <- 0L
    global_ess_accept_total <- 0L
    global_ess_total <- 0L
    
    local_eval_total <- 0L
    theta_eval_total <- 0L
    theta_update_total <- 0L
    
    save_id <- 0L
    
    for (iter in seq_len(n_iter)) {
      sigma2_eps <- sample_sigma2_eps_1d(y, u_curr, Dx, logtheta)
      
      tau <- update_tau_1d(tau, u_curr, c_ord, m, tau_bound)
      
      if (length(miss_idx) > 0 &&
          iter %% control$block_ess_every == 0) {
        for (bb in seq_len(control$n_blocks_per_iter)) {
          block_idx <- sample(
            miss_idx,
            min(control$ess_block_size, length(miss_idx))
          )
          
          ess <- update_u_ess_block_1d(
            y = y,
            u = u_curr,
            Dx = Dx,
            c_ord = c_ord,
            tau = tau,
            logtheta = logtheta,
            sigma2_eps = sigma2_eps,
            block_idx = block_idx
          )
          
          u_curr <- ess$u
          
          block_ess_eval_total <- block_ess_eval_total + ess$n_eval
          block_ess_accept_total <- block_ess_accept_total + as.integer(ess$accepted)
          block_ess_total <- block_ess_total + 1L
        }
      }
      
      if (length(miss_idx) > 0 &&
          control$global_ess_every > 0 &&
          iter %% control$global_ess_every == 0) {
        gess <- update_u_ess_block_1d(
          y = y,
          u = u_curr,
          Dx = Dx,
          c_ord = c_ord,
          tau = tau,
          logtheta = logtheta,
          sigma2_eps = sigma2_eps,
          block_idx = miss_idx
        )
        
        u_curr <- gess$u
        
        global_ess_eval_total <- global_ess_eval_total + gess$n_eval
        global_ess_accept_total <- global_ess_accept_total + as.integer(gess$accepted)
        global_ess_total <- global_ess_total + 1L
      }
      
      if (length(miss_idx) > 0) {
        if (iter %% control$full_local_every == 0) {
          local_idx <- miss_idx
        } else {
          local_idx <- sample(
            miss_idx,
            min(control$local_per_iter, length(miss_idx))
          )
        }
        
        loc <- update_u_local_z_slice_1d(
          y = y,
          u = u_curr,
          Dx = Dx,
          c_ord = c_ord,
          tau = tau,
          logtheta = logtheta,
          sigma2_eps = sigma2_eps,
          update_idx = local_idx
        )
        
        u_curr <- loc$u
        local_eval_total <- local_eval_total + loc$n_eval
      }
      
      if (iter %% control$theta_update_every == 0) {
        th <- update_logtheta_slice_1d(
          y = y,
          u = u_curr,
          Dx = Dx,
          logtheta = logtheta,
          sigma2_eps = sigma2_eps,
          theta_slice_width = theta_slice_width
        )
        
        logtheta <- th$logtheta
        
        theta_eval_total <- theta_eval_total + th$n_eval
        theta_update_total <- theta_update_total + 1L
      }
      
      logtheta_trace_all[iter, ] <- logtheta
      
      if (control$adapt_theta_width &&
          iter <= burn &&
          iter >= 200 &&
          iter %% control$adapt_every == 0) {
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
      
      if (iter > burn && ((iter - burn) %% thin == 0)) {
        save_id <- save_id + 1L
        
        samples_u_chain[save_id, ] <- u_curr
        samples_tau_chain[save_id, ] <- tau
        samples_logtheta_chain[save_id, ] <- logtheta
        samples_sigma2_chain[save_id] <- sigma2_eps
      }
    }
    
    samples_u_chain <- samples_u_chain[seq_len(save_id), , drop = FALSE]
    samples_tau_chain <- samples_tau_chain[seq_len(save_id), , drop = FALSE]
    samples_logtheta_chain <- samples_logtheta_chain[seq_len(save_id), , drop = FALSE]
    samples_sigma2_chain <- samples_sigma2_chain[seq_len(save_id)]
    
    stats <- data.frame(
      chain = chain_id,
      seed = chain_seed,
      saved = save_id,
      block_ess_eval_total = block_ess_eval_total,
      block_ess_accept_total = block_ess_accept_total,
      block_ess_total = block_ess_total,
      global_ess_eval_total = global_ess_eval_total,
      global_ess_accept_total = global_ess_accept_total,
      global_ess_total = global_ess_total,
      local_eval_total = local_eval_total,
      theta_eval_total = theta_eval_total,
      theta_update_total = theta_update_total
    )
    
    list(
      chain_id = chain_id,
      seed = chain_seed,
      samples_u = samples_u_chain,
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
    cat("Running EIV-GP with", n_chains, "chain(s).\n")
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
    u = lapply(chains, function(z) z$samples_u),
    tau = lapply(chains, function(z) z$samples_tau),
    logtheta = lapply(chains, function(z) z$samples_logtheta),
    sigma2 = lapply(chains, function(z) z$samples_sigma2)
  )
  
  samples_u <- do.call(rbind, samples_by_chain$u)
  samples_tau <- do.call(rbind, samples_by_chain$tau)
  samples_logtheta <- do.call(rbind, samples_by_chain$logtheta)
  samples_sigma2 <- unlist(samples_by_chain$sigma2)
  
  mcmc_draw_info <- data.frame(
    chain = rep(
      seq_len(n_chains),
      times = vapply(samples_by_chain$u, nrow, integer(1))
    ),
    draw_within_chain = unlist(
      lapply(samples_by_chain$u, function(mat) seq_len(nrow(mat)))
    )
  )
  
  chain_stats <- do.call(rbind, lapply(chains, function(z) z$stats))
  
  rhat_hyper <- data.frame(
    parameter = c("sigma_epsilon", "rho", "theta_x", "theta_u"),
    rhat = c(
      split_rhat(lapply(samples_by_chain$sigma2, function(v) sqrt(v))),
      split_rhat(lapply(samples_by_chain$logtheta, function(mat) exp(mat[, 1]))),
      split_rhat(lapply(samples_by_chain$logtheta, function(mat) exp(mat[, 2]))),
      split_rhat(lapply(samples_by_chain$logtheta, function(mat) exp(mat[, 3])))
    )
  )
  
  rhat_tau <- data.frame(
    parameter = paste0("tau", seq_len(m - 1)),
    rhat = sapply(seq_len(m - 1), function(j) {
      split_rhat(
        lapply(samples_by_chain$tau, function(mat) mat[, j])
      )
    })
  )
  
  if (length(miss_idx) > 0) {
    rhat_u <- data.frame(
      parameter = paste0("u[", miss_idx, "]"),
      global_index = miss_idx,
      rhat = sapply(seq_along(miss_idx), function(k) {
        jj <- miss_idx[k]
        
        split_rhat(
          lapply(samples_by_chain$u, function(mat) mat[, jj])
        )
      })
    )
  } else {
    rhat_u <- data.frame(
      parameter = character(0),
      global_index = integer(0),
      rhat = numeric(0)
    )
  }
  
  hyper_mat <- cbind(
    sigma_epsilon = sqrt(samples_sigma2),
    rho = exp(samples_logtheta[, 1]),
    theta_x = exp(samples_logtheta[, 2]),
    theta_u = exp(samples_logtheta[, 3])
  )
  
  tau_mat <- samples_tau
  colnames(tau_mat) <- paste0("tau", seq_len(ncol(tau_mat)))
  
  ess_key <- data.frame(
    parameter = colnames(cbind(hyper_mat, tau_mat)),
    ess = apply(cbind(hyper_mat, tau_mat), 2, ess_ips)
  )
  
  diagnostics_summary <- data.frame(
    n_chains = n_chains,
    n_iter = n_iter,
    burn = burn,
    saved_per_chain = mean(vapply(samples_by_chain$u, nrow, integer(1))),
    total_saved_draws = nrow(samples_u),
    max_rhat_hyper = safe_max(rhat_hyper$rhat),
    max_rhat_tau = safe_max(rhat_tau$rhat),
    median_rhat_missing_u = safe_median(rhat_u$rhat),
    max_rhat_missing_u = safe_max(rhat_u$rhat),
    min_ess_key = safe_min(ess_key$ess),
    time_seconds = as.numeric(mcmc_time["elapsed"])
  )
  
  list(
    data = list(
      x_raw = x_raw,
      x = x,
      x_center = x_center,
      x_scale = x_scale,
      y_raw = y_raw,
      y = y,
      y_center = y_center,
      y_scale = y_scale,
      c_ord = c_ord,
      u_true = u_true,
      u_obs = u_obs,
      calib_idx = calib_idx,
      miss_idx = miss_idx,
      tau_true = tau_true,
      m = m
    ),
    control = control,
    mcmc = list(
      samples_u = samples_u,
      samples_tau = samples_tau,
      samples_logtheta = samples_logtheta,
      samples_sigma2 = samples_sigma2,
      samples_by_chain = samples_by_chain,
      mcmc_draw_info = mcmc_draw_info,
      chain_stats = chain_stats
    ),
    diagnostics = list(
      rhat_hyper = rhat_hyper,
      rhat_tau = rhat_tau,
      rhat_u = rhat_u,
      ess_key = ess_key,
      summary = diagnostics_summary
    )
  )
}

############################################################
## Baseline GP models
############################################################

gp_mle_fit <- function(X, y) {
  X <- as.matrix(X)
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
    A <- rho^2 * R + diag(n)
    
    U <- try(safe_chol(A), silent = TRUE)
    
    if (inherits(U, "try-error")) {
      return(1e20)
    }
    
    Ainv_y <- solve_chol(U, y)
    logdetA <- 2 * sum(log(diag(U)))
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
  C <- K + sigma2 * diag(n)
  
  U <- safe_chol(C)
  alpha <- solve_chol(U, y)
  
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
  
  v <- forwardsolve(t(U), t(Kstar))
  var_lat <- rho^2 * sigma2 - colSums(v^2)
  var_lat <- pmax(var_lat, 1e-10)
  
  if (noisy) {
    var_lat <- var_lat + sigma2
  }
  
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

make_monotone_scores <- function(a, m) {
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

gp_mle_fit_learned_embedding <- function(x, c_ord, y, m,
                                         n_starts = 8L) {
  n <- length(y)
  
  Dx <- pairwise_sqdist(matrix(x, ncol = 1))
  
  nll <- function(par) {
    log_sigma2 <- par[1]
    log_rho <- par[2]
    log_theta_x <- par[3]
    log_theta_z <- par[4]
    
    a <- par[-(1:4)]
    
    sigma2 <- exp(log_sigma2)
    rho <- exp(log_rho)
    theta_x <- exp(log_theta_x)
    theta_z <- exp(log_theta_z)
    
    z_scores <- make_monotone_scores(a, m)
    z <- z_scores[c_ord]
    
    Dz <- pairwise_sqdist(matrix(z, ncol = 1))
    
    R <- exp(-theta_x * Dx - theta_z * Dz)
    A <- rho^2 * R + diag(n)
    
    U <- try(safe_chol(A), silent = TRUE)
    
    if (inherits(U, "try-error")) {
      return(1e20)
    }
    
    Ainv_y <- solve_chol(U, y)
    logdetA <- 2 * sum(log(diag(U)))
    quad <- sum(y * Ainv_y)
    
    0.5 * (n * log(2 * pi * sigma2) + logdetA + quad / sigma2)
  }
  
  lower <- c(
    log(1e-5),
    log(0.05),
    log(1e-4),
    log(1e-4),
    rep(-6, max(m - 2, 0))
  )
  
  upper <- c(
    log(5),
    log(100),
    log(100),
    log(100),
    rep(6, max(m - 2, 0))
  )
  
  make_init <- function() {
    c(
      log(0.05),
      log(3),
      log(0.5),
      log(0.5),
      rnorm(max(m - 2, 0), 0, 0.5)
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
  
  z_scores <- make_monotone_scores(opt$par[-(1:4)], m)
  
  list(
    par = opt$par[1:4],
    a = opt$par[-(1:4)],
    z_scores = z_scores,
    value = opt$value,
    convergence = opt$convergence,
    x = x,
    c_ord = c_ord,
    y = y,
    m = m
  )
}

gp_mle_predict_learned_embedding <- function(fit, x_star, c_star,
                                             noisy = TRUE) {
  X_train <- cbind(fit$x, fit$z_scores[fit$c_ord])
  X_star <- cbind(x_star, fit$z_scores[c_star])
  
  fake_fit <- list(
    par = fit$par,
    X = X_train,
    y = fit$y
  )
  
  gp_mle_predict(fake_fit, Xstar = X_star, noisy = noisy)
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

fit_embedding_baselines <- function(x_raw, y_raw, c_ord, m,
                                    n_starts_learned = 8L) {
  x_center <- mean(x_raw)
  x_scale <- sd(x_raw)
  x <- as.numeric((x_raw - x_center) / x_scale)
  
  y_center <- mean(y_raw)
  y_scale <- sd(y_raw)
  y <- as.numeric((y_raw - y_center) / y_scale)
  
  z_gauss <- qnorm(c_ord / (m + 1))
  
  fit_gauss <- gp_mle_fit(
    X = cbind(x, z_gauss),
    y = y
  )
  
  z_cm_scores <- conditional_mean_scores(c_ord, m)
  z_cm <- z_cm_scores[c_ord]
  
  fit_cm <- gp_mle_fit(
    X = cbind(x, z_cm),
    y = y
  )
  
  fit_learned <- gp_mle_fit_learned_embedding(
    x = x,
    c_ord = c_ord,
    y = y,
    m = m,
    n_starts = n_starts_learned
  )
  
  list(
    x_center = x_center,
    x_scale = x_scale,
    y_center = y_center,
    y_scale = y_scale,
    fit_gauss = fit_gauss,
    fit_cm = fit_cm,
    fit_learned = fit_learned,
    z_cm_scores = z_cm_scores
  )
}

predict_embedding_baseline_samples <- function(baselines,
                                               x_star_raw,
                                               c_star,
                                               m,
                                               n_draw = 1000) {
  x_star <- as.numeric((x_star_raw - baselines$x_center) / baselines$x_scale)
  
  z_gauss_star <- qnorm(c_star / (m + 1))
  
  draws_gauss_std <- sample_gp_mle_predictive(
    baselines$fit_gauss,
    Xstar = cbind(x_star, z_gauss_star),
    n_draw = n_draw
  )
  
  draws_gauss <- baselines$y_center + baselines$y_scale * draws_gauss_std
  
  z_cm_star <- baselines$z_cm_scores[c_star]
  
  draws_cm_std <- sample_gp_mle_predictive(
    baselines$fit_cm,
    Xstar = cbind(x_star, z_cm_star),
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

summarize_predictive_samples <- function(draw_mat, y_true,
                                         method,
                                         rep_id,
                                         n_calib,
                                         scenario) {
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

sample_oracle_y_star <- function(x_star, c_star, tau_true,
                                 scenario = "active",
                                 sigma_eps = 0.1,
                                 n_draw = 1000) {
  lower <- c(-Inf, tau_true)[c_star]
  upper <- c(tau_true, Inf)[c_star]
  
  u_star <- rtruncnorm_vec(
    mean = rep(0, n_draw),
    sd = 1,
    lower = rep(lower, n_draw),
    upper = rep(upper, n_draw)
  )
  
  f_star <- f0_1d(x_star, u_star, scenario = scenario)
  
  f_star + rnorm(n_draw, 0, sigma_eps)
}

sample_oracle_test_y <- function(x_test,
                                 c_test,
                                 tau_true,
                                 scenario,
                                 sigma_eps,
                                 n_draw = 1000) {
  n_test <- length(c_test)
  out <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  
  for (j in seq_len(n_test)) {
    out[, j] <- sample_oracle_y_star(
      x_star = x_test[j],
      c_star = c_test[j],
      tau_true = tau_true,
      scenario = scenario,
      sigma_eps = sigma_eps,
      n_draw = n_draw
    )
  }
  
  out
}

sample_eiv_test_y <- function(x_test_raw,
                              c_test,
                              fit_obj,
                              draw_ids,
                              n_per_draw = 1L) {
  samples_u <- fit_obj$mcmc$samples_u
  samples_tau <- fit_obj$mcmc$samples_tau
  samples_logtheta <- fit_obj$mcmc$samples_logtheta
  samples_sigma2 <- fit_obj$mcmc$samples_sigma2
  
  x_train <- fit_obj$data$x
  y_train <- fit_obj$data$y
  y_center <- fit_obj$data$y_center
  y_scale <- fit_obj$data$y_scale
  x_center <- fit_obj$data$x_center
  x_scale <- fit_obj$data$x_scale
  
  x_test <- as.numeric((x_test_raw - x_center) / x_scale)
  
  n_test <- length(c_test)
  n_draw <- length(draw_ids) * n_per_draw
  
  out <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  
  row_id <- 0L
  
  for (s in draw_ids) {
    for (rr in seq_len(n_per_draw)) {
      row_id <- row_id + 1L
      
      tau_s <- samples_tau[s, ]
      
      lower <- c(-Inf, tau_s)[c_test]
      upper <- c(tau_s, Inf)[c_test]
      
      u_star <- rtruncnorm_vec(
        mean = rep(0, n_test),
        sd = 1,
        lower = lower,
        upper = upper
      )
      
      pred <- gp_predict_draw(
        x_train = x_train,
        u_train = samples_u[s, ],
        y_train = y_train,
        x_star = x_test,
        u_star = u_star,
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