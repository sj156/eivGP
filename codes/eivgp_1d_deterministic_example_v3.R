############################################################
## replicate_eivgp_ordinal_robust.R
##
## Robust revised ordinal EIV-GP numerical example.
##
## Design goals:
##   1. No regression u | x.
##   2. Working latent prior u ~ N(0,1), after standardizing u_obs.
##   3. Exact GP MCMC.
##   4. Low-tuning sampler:
##        - blocked elliptical slice updates for u_mis;
##        - full-interval slice updates in z = Phi(u);
##        - occasional global elliptical slice updates;
##        - slice updates for GP hyperparameters on log scale;
##        - automatic burn-in adaptation of hyperparameter slice widths.
##   5. ggplot2 figures.
############################################################

rm(list = ls())
set.seed(20260701)

options(repos = c(CRAN = "https://cloud.r-project.org"))

load_or_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

load_or_install("ggplot2")
load_or_install("patchwork")
load_or_install("matrixStats")

theme_set(theme_bw(base_size = 12))

############################################################
## Multicore settings
############################################################

## Important on Mac: if BLAS/LAPACK is already multithreaded, running many
## forked chains can oversubscribe the CPU. These settings ask each R process
## to use one linear-algebra thread, while parallelism is across chains.
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

## Number of parallel MCMC chains.
## On a 16-core MacBook Pro, 4 is a safe default.
## You can try 6 or 8 if memory and thermals are fine.
n_chains <- 12L

## Leave a couple of cores free for the system.
mc_cores <- min(n_chains, max(1L, parallel::detectCores(logical = TRUE) - 2L))

## mclapply uses fork and works on macOS/Linux. It is not available on Windows.
use_mclapply <- (.Platform$OS.type != "windows" && mc_cores > 1L)

############################################################
## 0. User controls
############################################################

out_dir <- "eivgp_ordinal_output_robust"
fig_dir <- file.path(out_dir, "figures")
dir.create(out_dir, showWarnings = FALSE)
dir.create(fig_dir, showWarnings = FALSE)

## For non-expert users:
##   "fast"      : quicker, slightly less thorough MCMC;
##   "balanced"  : recommended default;
##   "thorough"  : slower, better exploration.
mcmc_preset <- "balanced"
mcmc_preset <- "thorough"

## Data size
n <- 100
n_calib <- 10
m <- 6

## True data-generating parameters
sigma_u_true <- 0.9 * pi
sigma_eps_true <- 0.1
tau_true_raw <- c(-1.5 * pi, -pi, 0, pi, 1.5 * pi)

## MCMC length
n_iter <- 6000
burn <- 2000
thin <- 1L 

## Number of posterior draws used in plotting/prediction
max_plot_draws <- 600
n_pred_truth <- 5000
n_pred_baseline <- 5000

############################################################
## 1. Default sampler controls
############################################################

make_default_control <- function(n, n_mis,
                                 preset = c("fast", "balanced", "thorough")) {
  preset <- match.arg(preset)
  
  if (preset == "fast") {
    local_frac <- 0.03
    theta_update_every <- 10
    block_ess_every <- 2
    n_blocks_per_iter <- 1
    global_ess_every <- 50
    full_local_every <- 100
  }
  
  if (preset == "balanced") {
    local_frac <- 0.06
    theta_update_every <- 5
    block_ess_every <- 1
    n_blocks_per_iter <- 1
    global_ess_every <- 25
    full_local_every <- 50
  }
  
  if (preset == "thorough") {
    local_frac <- 0.10
    theta_update_every <- 3
    block_ess_every <- 1
    n_blocks_per_iter <- 2
    global_ess_every <- 10
    full_local_every <- 25
  }
  
  local_per_iter <- min(
    n_mis,
    max(3, ceiling(local_frac * n_mis))
  )
  
  ess_block_size <- min(
    n_mis,
    max(5, ceiling(sqrt(n_mis)))
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
    adapt_every = 100,
    adapt_window = 500,
    theta_width_min = 0.20,
    theta_width_max = 2.50
  )
}

############################################################
## 2. Priors and numerical bounds
############################################################

## Since y is standardized, this is a weakly informative noise prior.
a_eps0 <- 2
b_eps0 <- 0.05

## Priors on log GP hyperparameters.
## Kernel: rho^2 sigma_eps^2 exp{-theta_x dx^2 - theta_u du^2}.
## Because x and u are standardized, theta_x and theta_u around 0.1--2
## are usually reasonable a priori.
logtheta_prior_mean <- c(log(3), log(0.5), log(0.5))
logtheta_prior_sd <- c(1.5, 1.5, 1.5)

## Broad numerical bounds for log hyperparameters.
theta_log_bounds <- rbind(
  log_rho     = c(log(0.05), log(100)),
  log_theta_x = c(log(1e-4), log(100)),
  log_theta_u = c(log(1e-4), log(100))
)

############################################################
## 3. Utility functions
############################################################

maximin_lhs_1d <- function(n, lower = -2, upper = 2) {
  z <- (seq_len(n) - runif(n)) / n
  z <- sample(z)
  lower + (upper - lower) * z
}

make_class <- function(u, tau) {
  as.integer(cut(u, breaks = c(-Inf, tau, Inf), labels = FALSE))
}

logdiffexp <- function(logx, logy) {
  logx + log1p(-exp(logy - logx))
}

rtruncnorm_vec <- function(mean, sd, lower, upper) {
  mean <- as.numeric(mean)
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  n <- length(mean)
  
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
  for (k in 0:7) {
    jj <- jitter * 10^k
    ans <- try(chol(A + jj * diag(n)), silent = TRUE)
    if (!inherits(ans, "try-error")) return(ans)
  }
  stop("Cholesky failed even after jitter.")
}

solve_chol <- function(U, b) {
  backsolve(U, forwardsolve(t(U), b))
}

############################################################
## 4. Slice samplers
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
  n_eval <- 1
  
  if (!is.finite(f0)) {
    return(list(x = x0, n_eval = n_eval))
  }
  
  logy <- f0 + log(runif(1))
  
  L <- x0 - runif(1) * w
  R <- L + w
  
  if (is.finite(lower)) L <- max(L, lower)
  if (is.finite(upper)) R <- min(R, upper)
  
  J <- floor(runif(1) * max_steps_out)
  K <- max_steps_out - 1 - J
  
  while (J > 0 && (!is.finite(lower) || L > lower)) {
    fL <- logf(L)
    n_eval <- n_eval + 1
    if (!is.finite(fL) || fL <= logy) break
    L <- L - w
    if (is.finite(lower)) L <- max(L, lower)
    J <- J - 1
  }
  
  while (K > 0 && (!is.finite(upper) || R < upper)) {
    fR <- logf(R)
    n_eval <- n_eval + 1
    if (!is.finite(fR) || fR <= logy) break
    R <- R + w
    if (is.finite(upper)) R <- min(R, upper)
    K <- K - 1
  }
  
  for (iter in seq_len(max_iter)) {
    x1 <- runif(1, L, R)
    f1 <- logf(x1)
    n_eval <- n_eval + 1
    
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
  ## Slice sampler with initial bracket equal to the full bounded support.
  ## This is slower than a tuned local bracket, but removes width tuning.
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
  n_eval <- 1
  
  if (!is.finite(f0)) {
    return(list(x = x0, n_eval = n_eval))
  }
  
  logy <- f0 + log(runif(1))
  
  L <- lower_eps
  R <- upper_eps
  
  for (iter in seq_len(max_iter)) {
    x1 <- runif(1, L, R)
    f1 <- logf(x1)
    n_eval <- n_eval + 1
    
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
## 5. Simulate data
############################################################

repeat {
  x_raw <- maximin_lhs_1d(n, lower = -2, upper = 2)
  u_true_raw <- rnorm(n, mean = 0, sd = sigma_u_true)
  c_ord <- make_class(u_true_raw, tau_true_raw)
  
  if (all(tabulate(c_ord, nbins = m) >= 3)) {
    break
  }
}

eps <- rnorm(n, mean = 0, sd = sigma_eps_true)
y_raw <- cos(u_true_raw) + eps

calib_idx <- sort(sample(seq_len(n), n_calib))
miss_idx <- setdiff(seq_len(n), calib_idx)

u_obs_raw <- rep(NA_real_, n)
u_obs_raw[calib_idx] <- u_true_raw[calib_idx]

## Standardize x for robust default length-scale priors.
x_center <- mean(x_raw)
x_scale <- sd(x_raw)
x <- as.numeric((x_raw - x_center) / x_scale)

## Standardize y.
y_center <- mean(y_raw)
y_scale <- sd(y_raw)
y <- as.numeric((y_raw - y_center) / y_scale)

## Standardize u using calibrated u observations.
## The model's working latent prior is N(0,1) on this scale.
u_center <- mean(u_obs_raw[calib_idx])
u_scale <- sd(u_obs_raw[calib_idx])

u_true <- as.numeric((u_true_raw - u_center) / u_scale)
u_obs <- rep(NA_real_, n)
u_obs[calib_idx] <- as.numeric((u_obs_raw[calib_idx] - u_center) / u_scale)

tau_true <- as.numeric((tau_true_raw - u_center) / u_scale)
tau_bound <- max(8, max(abs(u_obs[calib_idx]), na.rm = TRUE) + 2)

control <- make_default_control(
  n = n,
  n_mis = length(miss_idx),
  preset = mcmc_preset
)

cat("Class counts:\n")
print(tabulate(c_ord, nbins = m))
cat("Calibration indices:\n")
print(calib_idx)
cat("x standardization center/scale:\n")
print(c(center = x_center, scale = x_scale))
cat("y standardization center/scale:\n")
print(c(center = y_center, scale = y_scale))
cat("u standardization center/scale:\n")
print(c(center = u_center, scale = u_scale))
cat("Sampler preset:\n")
print(control)

############################################################
## 6. Initialize EIV-GP state
############################################################

initialize_tau <- function(c_ord, u_obs, calib_idx, m, tau_bound = 8) {
  n <- length(c_ord)
  counts <- tabulate(c_ord, nbins = m)
  probs <- cumsum(counts)[1:(m - 1)] / n
  probs <- pmin(pmax(probs, 0.03), 0.97)
  
  tau <- qnorm(probs)
  tau <- pmin(pmax(tau, -tau_bound + 1e-3), tau_bound - 1e-3)
  
  obs_c <- c_ord[calib_idx]
  obs_u <- u_obs[calib_idx]
  
  eps_ord <- 1e-4
  
  for (j in seq_len(m - 1)) {
    lower <- max(c(-tau_bound, obs_u[obs_c <= j]), na.rm = TRUE)
    upper <- min(c( tau_bound, obs_u[obs_c >  j]), na.rm = TRUE)
    
    if (lower >= upper) {
      stop("Calibrated u values are incompatible with ordinal classes.")
    }
    
    tau[j] <- min(max(tau[j], lower + eps_ord), upper - eps_ord)
  }
  
  tau
}

tau <- initialize_tau(c_ord, u_obs, calib_idx, m, tau_bound)

u_curr <- u_obs

lower_all <- c(-Inf, tau)[c_ord]
upper_all <- c(tau, Inf)[c_ord]

u_curr[miss_idx] <- rtruncnorm_vec(
  mean = rep(0, length(miss_idx)),
  sd = 1,
  lower = lower_all[miss_idx],
  upper = upper_all[miss_idx]
)

## Initialize GP parameters using standardized data.
sigma2_eps <- 0.05

## Median-distance heuristic for initial inverse length-scales.
Dx_init <- pairwise_sqdist(matrix(x, ncol = 1))
du_init <- pairwise_sqdist(matrix(u_curr, ncol = 1))

med_dx <- median(Dx_init[upper.tri(Dx_init)], na.rm = TRUE)
med_du <- median(du_init[upper.tri(du_init)], na.rm = TRUE)

theta_x_init <- ifelse(is.finite(med_dx) && med_dx > 0, 1 / med_dx, 0.5)
theta_u_init <- ifelse(is.finite(med_du) && med_du > 0, 1 / med_du, 0.5)
rho_init <- sqrt(max(var(y) / sigma2_eps - 1, 1))

logtheta <- c(log(rho_init), log(theta_x_init), log(theta_u_init))
logtheta <- pmin(pmax(logtheta, theta_log_bounds[, 1]), theta_log_bounds[, 2])

Dx <- pairwise_sqdist(matrix(x, ncol = 1))

############################################################
## 7. GP likelihood and posterior components
############################################################

gp_state <- function(y, u, Dx, logtheta, sigma2_eps) {
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

check_constraints <- function(u, c_ord, tau) {
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

gp_loglik_with_constraints <- function(y, u, Dx, c_ord, tau,
                                       logtheta, sigma2_eps) {
  if (!check_constraints(u, c_ord, tau)) return(-Inf)
  gp_state(y, u, Dx, logtheta, sigma2_eps)$loglik
}

theta_logpost <- function(y, u, Dx, logtheta, sigma2_eps) {
  lp <- log_prior_logtheta(logtheta)
  if (!is.finite(lp)) return(-Inf)
  
  gp_state(y, u, Dx, logtheta, sigma2_eps)$loglik + lp
}

############################################################
## 8. MCMC update functions
############################################################

update_tau <- function(tau, u, c_ord, m, tau_bound = 8) {
  tau_new <- tau
  
  for (j in seq_len(m - 1)) {
    ## Cumulative constraints:
    ## max_{c <= j} u_i < tau_j < min_{c > j} u_i.
    L <- max(c(-tau_bound, u[c_ord <= j]), na.rm = TRUE)
    U <- min(c( tau_bound, u[c_ord >  j]), na.rm = TRUE)
    
    if (is.finite(L) && is.finite(U) && L < U) {
      tau_new[j] <- runif(1, L, U)
    }
  }
  
  tau_new
}

sample_sigma2_eps <- function(y, u, Dx, logtheta) {
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

update_u_ess_block <- function(y, u, Dx, c_ord, tau, logtheta, sigma2_eps,
                               block_idx, max_try = 300) {
  if (length(block_idx) == 0) {
    return(list(u = u, n_eval = 0, accepted = TRUE))
  }
  
  u_block <- u[block_idx]
  
  loglik_fun <- function(u_block_prop) {
    u_prop <- u
    u_prop[block_idx] <- u_block_prop
    
    gp_loglik_with_constraints(
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
  
  n_eval <- 1
  
  for (try_id in seq_len(max_try)) {
    u_block_prop <- u_block * cos(angle) + nu * sin(angle)
    loglik_prop <- loglik_fun(u_block_prop)
    n_eval <- n_eval + 1
    
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

update_u_local_z_slice <- function(y, u, Dx, c_ord, tau,
                                   logtheta, sigma2_eps,
                                   update_idx) {
  ## Local update in z_i = Phi(u_i).
  ## Under u_i ~ N(0,1), z_i is uniform on (0,1).
  ## Conditional on class c_i, z_i is restricted to
  ## [Phi(tau_{c_i-1}), Phi(tau_{c_i})].
  ##
  ## Hence the prior and Jacobian cancel, and the log target in z_i
  ## is just the GP log likelihood on the allowed interval.
  n_eval_total <- 0
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
      if (z <= z_lower || z >= z_upper) return(-Inf)
      
      u_prop_i <- qnorm(z)
      
      if (!(u_prop_i > lower_u && u_prop_i <= upper_u)) return(-Inf)
      
      u_prop <- u
      u_prop[i] <- u_prop_i
      
      gp_state(y, u_prop, Dx, logtheta, sigma2_eps)$loglik
    }
    
    ans <- full_interval_slice_update(
      x0 = z0,
      logf = logf_z,
      lower = z_lower,
      upper = z_upper,
      max_iter = 200
    )
    
    n_eval_total <- n_eval_total + ans$n_eval
    u[i] <- qnorm(ans$x)
  }
  
  list(u = u, n_eval = n_eval_total)
}

update_logtheta_slice <- function(y, u, Dx, logtheta, sigma2_eps,
                                  theta_slice_width) {
  n_eval_total <- 0
  
  for (j in seq_along(logtheta)) {
    logf_j <- function(val) {
      lt <- logtheta
      lt[j] <- val
      theta_logpost(y, u, Dx, lt, sigma2_eps)
    }
    
    ans <- bounded_slice_update(
      x0 = logtheta[j],
      logf = logf_j,
      w = theta_slice_width[j],
      lower = theta_log_bounds[j, 1],
      upper = theta_log_bounds[j, 2],
      max_steps_out = 30,
      max_iter = 100
    )
    
    logtheta[j] <- ans$x
    n_eval_total <- n_eval_total + ans$n_eval
  }
  
  list(logtheta = logtheta, n_eval = n_eval_total)
}

############################################################
## 9. Run MCMC
############################################################
############################################################
## 9. Run MCMC using parallel independent chains
############################################################

n_save <- floor((n_iter - burn) / thin)

initialize_chain_state <- function(chain_id, chain_seed) {
  set.seed(chain_seed)
  
  ## Initialize cut points.
  tau0 <- initialize_tau(c_ord, u_obs, calib_idx, m, tau_bound)
  
  ## Initialize latent variables from the truncated working prior.
  u0 <- u_obs
  
  lower_all0 <- c(-Inf, tau0)[c_ord]
  upper_all0 <- c(tau0, Inf)[c_ord]
  
  u0[miss_idx] <- rtruncnorm_vec(
    mean = rep(0, length(miss_idx)),
    sd = 1,
    lower = lower_all0[miss_idx],
    upper = upper_all0[miss_idx]
  )
  
  ## A few cheap prior-threshold refreshes create more dispersed starts.
  for (rr in seq_len(3)) {
    tau0 <- update_tau(tau0, u0, c_ord, m, tau_bound)
    
    lower_all0 <- c(-Inf, tau0)[c_ord]
    upper_all0 <- c(tau0, Inf)[c_ord]
    
    u0[miss_idx] <- rtruncnorm_vec(
      mean = rep(0, length(miss_idx)),
      sd = 1,
      lower = lower_all0[miss_idx],
      upper = upper_all0[miss_idx]
    )
  }
  
  ## Initialize sigma_epsilon^2 on standardized y scale.
  sigma2_eps0 <- exp(log(0.05) + rnorm(1, 0, 0.5))
  sigma2_eps0 <- min(max(sigma2_eps0, 1e-4), 2)
  
  ## Median-distance initialization for length scales.
  Dx_init <- pairwise_sqdist(matrix(x, ncol = 1))
  Du_init <- pairwise_sqdist(matrix(u0, ncol = 1))
  
  med_dx <- median(Dx_init[upper.tri(Dx_init)], na.rm = TRUE)
  med_du <- median(Du_init[upper.tri(Du_init)], na.rm = TRUE)
  
  theta_x_init <- ifelse(is.finite(med_dx) && med_dx > 0, 1 / med_dx, 0.5)
  theta_u_init <- ifelse(is.finite(med_du) && med_du > 0, 1 / med_du, 0.5)
  rho_init <- sqrt(max(var(y) / sigma2_eps0 - 1, 1))
  
  logtheta0 <- c(log(rho_init), log(theta_x_init), log(theta_u_init))
  
  ## Overdisperse chain starts mildly.
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
  
  state <- initialize_chain_state(chain_id, chain_seed)
  
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
  
  block_ess_eval_total <- 0
  block_ess_accept_total <- 0
  block_ess_total <- 0
  
  global_ess_eval_total <- 0
  global_ess_accept_total <- 0
  global_ess_total <- 0
  
  local_eval_total <- 0
  theta_eval_total <- 0
  theta_update_total <- 0
  
  save_id <- 0
  
  for (iter in seq_len(n_iter)) {
    
    ## 1. Update sigma_epsilon^2.
    sigma2_eps <- sample_sigma2_eps(y, u_curr, Dx, logtheta)
    
    ## 2. Update cut points.
    tau <- update_tau(tau, u_curr, c_ord, m, tau_bound)
    
    ## 3. Blocked elliptical slice updates for missing latent variables.
    if (iter %% control$block_ess_every == 0) {
      for (bb in seq_len(control$n_blocks_per_iter)) {
        block_idx <- sample(
          miss_idx,
          min(control$ess_block_size, length(miss_idx))
        )
        
        ess <- update_u_ess_block(
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
        block_ess_total <- block_ess_total + 1
      }
    }
    
    ## 4. Occasional global ESS move.
    if (control$global_ess_every > 0 &&
        iter %% control$global_ess_every == 0) {
      
      gess <- update_u_ess_block(
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
      global_ess_total <- global_ess_total + 1
    }
    
    ## 5. Local full-interval z-slice updates.
    if (iter %% control$full_local_every == 0) {
      local_idx <- miss_idx
    } else {
      local_idx <- sample(
        miss_idx,
        min(control$local_per_iter, length(miss_idx))
      )
    }
    
    loc <- update_u_local_z_slice(
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
    
    ## 6. Slice update GP hyperparameters on log scale.
    if (iter %% control$theta_update_every == 0) {
      th <- update_logtheta_slice(
        y = y,
        u = u_curr,
        Dx = Dx,
        logtheta = logtheta,
        sigma2_eps = sigma2_eps,
        theta_slice_width = theta_slice_width
      )
      
      logtheta <- th$logtheta
      theta_eval_total <- theta_eval_total + th$n_eval
      theta_update_total <- theta_update_total + 1
    }
    
    logtheta_trace_all[iter, ] <- logtheta
    
    ## 7. Burn-in-only adaptation of hyperparameter slice widths.
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
    
    ## 8. Save posterior draws.
    if (iter > burn && ((iter - burn) %% thin == 0)) {
      save_id <- save_id + 1
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

chain_seeds <- 20260701L + 10000L * seq_len(n_chains)

cat("\nStarting parallel MCMC...\n")
cat("Number of chains:", n_chains, "\n")
cat("Number of mclapply cores:", ifelse(use_mclapply, mc_cores, 1L), "\n")

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

cat("\nParallel MCMC completed.\n")
print(mcmc_time)

## Store chain-specific samples.
samples_by_chain <- list(
  u = lapply(chains, function(z) z$samples_u),
  tau = lapply(chains, function(z) z$samples_tau),
  logtheta = lapply(chains, function(z) z$samples_logtheta),
  sigma2 = lapply(chains, function(z) z$samples_sigma2)
)

## Pooled samples for posterior summaries and plots.
samples_u <- do.call(rbind, samples_by_chain$u)
samples_tau <- do.call(rbind, samples_by_chain$tau)
samples_logtheta <- do.call(rbind, samples_by_chain$logtheta)
samples_sigma2 <- unlist(samples_by_chain$sigma2)

## Draw metadata, useful for diagnostics.
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
print(chain_stats)

## Aggregate counters used later by the diagnostic code.
block_ess_eval_total <- sum(chain_stats$block_ess_eval_total)
block_ess_accept_total <- sum(chain_stats$block_ess_accept_total)
block_ess_total <- sum(chain_stats$block_ess_total)

global_ess_eval_total <- sum(chain_stats$global_ess_eval_total)
global_ess_accept_total <- sum(chain_stats$global_ess_accept_total)
global_ess_total <- sum(chain_stats$global_ess_total)

local_eval_total <- sum(chain_stats$local_eval_total)

theta_eval_total <- sum(chain_stats$theta_eval_total)
theta_update_total <- sum(chain_stats$theta_update_total)

theta_slice_width <- Reduce(
  "+",
  lapply(chains, function(z) z$theta_slice_width_final)
) / n_chains

cat("\nCombined posterior draws:", nrow(samples_u), "\n")
cat("Average saved draws per chain:", mean(vapply(samples_by_chain$u, nrow, integer(1))), "\n")

cat("\nAggregated sampler behavior:\n")
cat("Block ESS accepted fraction:",
    block_ess_accept_total / max(block_ess_total, 1), "\n")
cat("Average block ESS likelihood evaluations:",
    block_ess_eval_total / max(block_ess_total, 1), "\n")
cat("Global ESS accepted fraction:",
    global_ess_accept_total / max(global_ess_total, 1), "\n")
cat("Average global ESS likelihood evaluations:",
    global_ess_eval_total / max(global_ess_total, 1), "\n")
cat("Average local-slice likelihood evaluations per chain iteration:",
    local_eval_total / max(n_iter * n_chains, 1), "\n")
cat("Average theta-slice likelihood evaluations per theta update:",
    theta_eval_total / max(theta_update_total, 1), "\n")
cat("Average final theta slice widths:\n")
print(theta_slice_width)

############################################################
## 10. Posterior summaries for latent u
############################################################

u_draws_raw <- u_center + u_scale * samples_u

u_post_mean_raw <- colMeans(u_draws_raw)
u_post_lo_raw <- apply(u_draws_raw, 2, quantile, probs = 0.025)
u_post_hi_raw <- apply(u_draws_raw, 2, quantile, probs = 0.975)

coverage_miss <- u_true_raw[miss_idx] >= u_post_lo_raw[miss_idx] &
  u_true_raw[miss_idx] <= u_post_hi_raw[miss_idx]

############################################################
## 11. GP prediction functions
############################################################

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
  
  Dxs <- pairwise_sqdist(matrix(x_star, ncol = 1), matrix(x_train, ncol = 1))
  Dus <- pairwise_sqdist(matrix(u_star, ncol = 1), matrix(u_train, ncol = 1))
  
  R_star <- exp(-theta_x * Dxs - theta_u * Dus)
  K_star <- rho^2 * sigma2_eps * R_star
  
  mu <- as.numeric(K_star %*% alpha)
  
  v <- forwardsolve(t(U), t(K_star))
  var_lat <- rho^2 * sigma2_eps - colSums(v^2)
  var_lat <- pmax(var_lat, 1e-10)
  
  if (noisy) var_lat <- var_lat + sigma2_eps
  
  list(mean = mu, var = var_lat)
}

############################################################
## 12. Baseline GP models
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
    if (inherits(U, "try-error")) return(1e20)
    
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
    upper = upper
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
  
  if (noisy) var_lat <- var_lat + sigma2
  
  list(mean = mu, var = var_lat)
}

fit_gp_calib <- gp_mle_fit(
  X = cbind(x[calib_idx], u_obs[calib_idx]),
  y = y[calib_idx]
)

z_unif <- c_ord / (m + 1)
z_gauss <- qnorm(c_ord / (m + 1))

fit_gp_unif <- gp_mle_fit(
  X = cbind(x, z_unif),
  y = y
)

fit_gp_gauss <- gp_mle_fit(
  X = cbind(x, z_gauss),
  y = y
)

cat("\nBaseline GP convergence codes:\n")
print(c(
  calib = fit_gp_calib$convergence,
  unif = fit_gp_unif$convergence,
  gauss = fit_gp_gauss$convergence
))

############################################################
## 13. Latent function inference at x_raw = 0
############################################################

x_star_raw <- 0
x_star <- as.numeric((x_star_raw - x_center) / x_scale)

u_grid_raw <- seq(-2 * pi, 2 * pi, length.out = 160)
u_grid <- as.numeric((u_grid_raw - u_center) / u_scale)
x_grid <- rep(x_star, length(u_grid))

plot_draw_ids <- seq_len(nrow(samples_u))
if (length(plot_draw_ids) > max_plot_draws) {
  plot_draw_ids <- sample(plot_draw_ids, max_plot_draws)
}

post_cores <- if (use_mclapply) {
  min(mc_cores, length(plot_draw_ids))
} else {
  1L
}

f_eiv_list <- if (use_mclapply && post_cores > 1L) {
  parallel::mclapply(
    seq_along(plot_draw_ids),
    function(ii) {
      s <- plot_draw_ids[ii]
      
      pred <- gp_predict_draw(
        x_train = x,
        u_train = samples_u[s, ],
        y_train = y,
        x_star = x_grid,
        u_star = u_grid,
        logtheta = samples_logtheta[s, ],
        sigma2_eps = samples_sigma2[s],
        noisy = FALSE
      )
      
      f_std <- pred$mean + sqrt(pred$var) * rnorm(length(u_grid))
      y_center + y_scale * f_std
    },
    mc.cores = post_cores,
    mc.set.seed = TRUE
  )
} else {
  lapply(
    seq_along(plot_draw_ids),
    function(ii) {
      s <- plot_draw_ids[ii]
      
      pred <- gp_predict_draw(
        x_train = x,
        u_train = samples_u[s, ],
        y_train = y,
        x_star = x_grid,
        u_star = u_grid,
        logtheta = samples_logtheta[s, ],
        sigma2_eps = samples_sigma2[s],
        noisy = FALSE
      )
      
      f_std <- pred$mean + sqrt(pred$var) * rnorm(length(u_grid))
      y_center + y_scale * f_std
    }
  )
}

f_eiv_samps <- do.call(rbind, f_eiv_list)

f_eiv_mean <- colMeans(f_eiv_samps)
f_eiv_lo <- apply(f_eiv_samps, 2, quantile, probs = 0.025)
f_eiv_hi <- apply(f_eiv_samps, 2, quantile, probs = 0.975)


pred_calib <- gp_mle_predict(
  fit_gp_calib,
  Xstar = cbind(x_grid, u_grid),
  noisy = FALSE
)

f_gp_mean <- y_center + y_scale * pred_calib$mean
f_gp_lo <- y_center + y_scale * (pred_calib$mean - 1.96 * sqrt(pred_calib$var))
f_gp_hi <- y_center + y_scale * (pred_calib$mean + 1.96 * sqrt(pred_calib$var))

############################################################
## 14. Mixed-input prediction at x_raw* = 0, c* = 1,...,6
############################################################

sample_eiv_y_star <- function(c_star, draw_ids) {
  out <- numeric(length(draw_ids))
  
  for (ii in seq_along(draw_ids)) {
    s <- draw_ids[ii]
    tau_s <- samples_tau[s, ]
    
    lower <- c(-Inf, tau_s)[c_star]
    upper <- c(tau_s, Inf)[c_star]
    
    u_star <- rtruncnorm_vec(
      mean = 0,
      sd = 1,
      lower = lower,
      upper = upper
    )
    
    pred <- gp_predict_draw(
      x_train = x,
      u_train = samples_u[s, ],
      y_train = y,
      x_star = x_star,
      u_star = u_star,
      logtheta = samples_logtheta[s, ],
      sigma2_eps = samples_sigma2[s],
      noisy = TRUE
    )
    
    y_std <- pred$mean + sqrt(pred$var) * rnorm(1)
    out[ii] <- y_center + y_scale * y_std
  }
  
  out
}

sample_truth_y_star <- function(c_star, n_draw = n_pred_truth) {
  lower <- c(-Inf, tau_true_raw)[c_star]
  upper <- c(tau_true_raw, Inf)[c_star]
  
  uu <- rtruncnorm_vec(
    mean = rep(0, n_draw),
    sd = sigma_u_true,
    lower = rep(lower, n_draw),
    upper = rep(upper, n_draw)
  )
  
  cos(uu) + rnorm(n_draw, 0, sigma_eps_true)
}

sample_embed_y_star <- function(fit, z_star, n_draw = n_pred_baseline) {
  pred <- gp_mle_predict(
    fit,
    Xstar = cbind(x_star, z_star),
    noisy = TRUE
  )
  
  y_std <- rnorm(n_draw, pred$mean, sqrt(pred$var))
  y_center + y_scale * y_std
}

pred_cores <- if (use_mclapply) {
  min(mc_cores, m)
} else {
  1L
}

pred_list <- if (use_mclapply && pred_cores > 1L) {
  parallel::mclapply(
    seq_len(m),
    function(cc) {
      list(
        truth = sample_truth_y_star(cc, n_pred_truth),
        eiv = sample_eiv_y_star(cc, plot_draw_ids),
        gp_unif = sample_embed_y_star(
          fit_gp_unif,
          cc / (m + 1),
          n_pred_baseline
        ),
        gp_gauss = sample_embed_y_star(
          fit_gp_gauss,
          qnorm(cc / (m + 1)),
          n_pred_baseline
        )
      )
    },
    mc.cores = pred_cores,
    mc.set.seed = TRUE
  )
} else {
  lapply(
    seq_len(m),
    function(cc) {
      list(
        truth = sample_truth_y_star(cc, n_pred_truth),
        eiv = sample_eiv_y_star(cc, plot_draw_ids),
        gp_unif = sample_embed_y_star(
          fit_gp_unif,
          cc / (m + 1),
          n_pred_baseline
        ),
        gp_gauss = sample_embed_y_star(
          fit_gp_gauss,
          qnorm(cc / (m + 1)),
          n_pred_baseline
        )
      )
    }
  )
}

############################################################
## 15. Save results
############################################################

res <- list(
  data = list(
    x_raw = x_raw,
    x = x,
    x_center = x_center,
    x_scale = x_scale,
    y_raw = y_raw,
    y = y,
    y_center = y_center,
    y_scale = y_scale,
    u_true_raw = u_true_raw,
    u_true = u_true,
    u_obs_raw = u_obs_raw,
    u_obs = u_obs,
    u_center = u_center,
    u_scale = u_scale,
    c_ord = c_ord,
    calib_idx = calib_idx,
    miss_idx = miss_idx,
    tau_true_raw = tau_true_raw,
    tau_true = tau_true
  ),
  control = control,
  mcmc = list(
    samples_u = samples_u,
    samples_tau = samples_tau,
    samples_logtheta = samples_logtheta,
    samples_sigma2 = samples_sigma2,
    theta_slice_width_final = theta_slice_width
  ),
  summaries = list(
    u_post_mean_raw = u_post_mean_raw,
    u_post_lo_raw = u_post_lo_raw,
    u_post_hi_raw = u_post_hi_raw,
    coverage_miss = coverage_miss
  ),
  baselines = list(
    fit_gp_calib = fit_gp_calib,
    fit_gp_unif = fit_gp_unif,
    fit_gp_gauss = fit_gp_gauss
  ),
  predictions = list(
    u_grid_raw = u_grid_raw,
    f_eiv_mean = f_eiv_mean,
    f_eiv_lo = f_eiv_lo,
    f_eiv_hi = f_eiv_hi,
    f_gp_mean = f_gp_mean,
    f_gp_lo = f_gp_lo,
    f_gp_hi = f_gp_hi,
    pred_list = pred_list
  )
)

saveRDS(res, file.path(out_dir, "eivgp_ordinal_results_robust.rds"))

############################################################
## 16. ggplot2 figures
############################################################

class_cols <- setNames(
  c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"),
  as.character(seq_len(m))
)

############################################################
## Figure 1: ground truth in latent space
############################################################

df_data <- data.frame(
  u = u_true_raw,
  y = y_raw,
  class = factor(c_ord)
)

df_curve <- data.frame(
  u = u_grid_raw,
  y = cos(u_grid_raw)
)

p_ground <- ggplot(df_data, aes(x = u, y = y, color = class)) +
  geom_point(size = 2, alpha = 0.85) +
  geom_line(
    data = df_curve,
    aes(x = u, y = y),
    inherit.aes = FALSE,
    color = "orange",
    linewidth = 1.2
  ) +
  geom_vline(
    xintercept = tau_true_raw,
    color = "black",
    linewidth = 0.6
  ) +
  scale_color_manual(values = class_cols, name = "Class") +
  labs(
    x = "Latent u",
    y = "y",
    title = "Synthetic data in latent space"
  ) +
  theme(legend.position = "right")

ggsave(
  file.path(fig_dir, "fig1_ground_truth_latent_space.pdf"),
  p_ground,
  width = 7,
  height = 4.5
)

############################################################
## Figure 2: observed data
############################################################

df_obs <- data.frame(
  x = x_raw,
  c = factor(c_ord),
  y = y_raw
)

df_cal <- data.frame(
  x = x_raw[calib_idx],
  u_obs = u_obs_raw[calib_idx],
  y = y_raw[calib_idx]
)

p_obs1 <- ggplot(df_obs, aes(x = c, y = x, color = y)) +
  geom_jitter(width = 0.12, height = 0, size = 2.2, alpha = 0.9) +
  scale_color_gradient2(
    low = "navy",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    name = "y"
  ) +
  labs(
    x = "Ordinal input c",
    y = "x",
    title = "Observed mixed input"
  )

p_obs2 <- ggplot(df_cal, aes(x = u_obs, y = x, color = y)) +
  geom_point(size = 2.6, alpha = 0.9) +
  scale_color_gradient2(
    low = "navy",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    name = "y"
  ) +
  labs(
    x = expression(u[obs]),
    y = "x",
    title = "Calibrated latent observations"
  )

p_obs <- p_obs1 + p_obs2 + plot_layout(ncol = 2)

ggsave(
  file.path(fig_dir, "fig2_observed_data.pdf"),
  p_obs,
  width = 10,
  height = 4.5
)

############################################################
## Figure 3: latent imputation
############################################################

df_imp <- data.frame(
  true_u = u_true_raw[miss_idx],
  mean_u = u_post_mean_raw[miss_idx],
  lo_u = u_post_lo_raw[miss_idx],
  hi_u = u_post_hi_raw[miss_idx],
  class = factor(c_ord[miss_idx]),
  covered = factor(
    coverage_miss,
    levels = c(TRUE, FALSE),
    labels = c("Covered", "Not covered")
  )
)

p_imp_scatter <- ggplot(df_imp, aes(x = true_u, y = mean_u, color = class)) +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 0.8) +
  scale_color_manual(values = class_cols, name = "Class") +
  labs(
    x = "True u",
    y = "Posterior mean of imputed u",
    title = "Latent imputation"
  )

df_imp_ord <- df_imp[order(df_imp$true_u), ]
df_imp_ord$id <- seq_len(nrow(df_imp_ord))

p_imp_ci <- ggplot(df_imp_ord, aes(y = id)) +
  geom_segment(
    aes(x = lo_u, xend = hi_u, yend = id, color = covered),
    linewidth = 0.7
  ) +
  geom_point(aes(x = true_u), size = 1.6) +
  scale_color_manual(
    values = c("Covered" = "steelblue", "Not covered" = "firebrick"),
    name = NULL
  ) +
  labs(
    x = "u",
    y = "Missing observations, sorted by true u",
    title = "95% credible intervals for u"
  )

p_imp <- p_imp_scatter + p_imp_ci + plot_layout(ncol = 2)

ggsave(
  file.path(fig_dir, "fig3_latent_imputation.pdf"),
  p_imp,
  width = 11,
  height = 4.8
)

############################################################
## Figure 4: latent function inference
############################################################

df_fun <- rbind(
  data.frame(
    method = "Complete-case GP",
    u = u_grid_raw,
    truth = cos(u_grid_raw),
    mean = f_gp_mean,
    lo = f_gp_lo,
    hi = f_gp_hi
  ),
  data.frame(
    method = "EIV-GP",
    u = u_grid_raw,
    truth = cos(u_grid_raw),
    mean = f_eiv_mean,
    lo = f_eiv_lo,
    hi = f_eiv_hi
  )
)

df_fun$method <- factor(df_fun$method, levels = c("Complete-case GP", "EIV-GP"))

p_fun <- ggplot(df_fun, aes(x = u)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "skyblue", alpha = 0.45) +
  geom_line(aes(y = mean), color = "blue", linewidth = 0.9, linetype = "dashed") +
  geom_line(aes(y = truth), color = "orange", linewidth = 1.1) +
  facet_wrap(~method, ncol = 2) +
  labs(
    x = "u",
    y = expression(f(0, u)),
    title = "Posterior inference for the latent regression function"
  )

ggsave(
  file.path(fig_dir, "fig4_function_inference.pdf"),
  p_fun,
  width = 10,
  height = 4.5
)

############################################################
## Figure 5: mixed-input predictive distributions
############################################################

method_labels <- c(
  truth = "Truth",
  eiv = "EIV-GP",
  gp_unif = "GP-Unif",
  gp_gauss = "GP-Gaussian"
)

df_pred <- do.call(
  rbind,
  lapply(seq_len(m), function(cc) {
    dd <- pred_list[[cc]]
    
    do.call(
      rbind,
      lapply(names(dd), function(nm) {
        data.frame(
          class = factor(
            paste0("c* = ", cc),
            levels = paste0("c* = ", seq_len(m))
          ),
          method = method_labels[[nm]],
          y = as.numeric(dd[[nm]])
        )
      })
    )
  })
)

df_pred$method <- factor(
  df_pred$method,
  levels = c("Truth", "EIV-GP", "GP-Unif", "GP-Gaussian")
)

method_cols <- c(
  "Truth" = "black",
  "EIV-GP" = "firebrick",
  "GP-Unif" = "steelblue",
  "GP-Gaussian" = "darkgreen"
)

p_pred_density <- ggplot(df_pred, aes(x = y, color = method)) +
  geom_density(linewidth = 0.9) +
  facet_wrap(~class, ncol = 3, scales = "free_y") +
  scale_color_manual(values = method_cols, name = NULL) +
  labs(
    x = expression(y^"*"),
    y = "Density",
    title = expression(paste("Mixed-input prediction at ", x^"*", " = 0"))
  ) +
  theme(legend.position = "bottom")

ggsave(
  file.path(fig_dir, "fig5_mixed_input_prediction_density.pdf"),
  p_pred_density,
  width = 12,
  height = 7.5
)

p_pred_box <- ggplot(df_pred, aes(x = method, y = y, fill = method)) +
  geom_boxplot(outlier.size = 0.25, alpha = 0.75) +
  facet_wrap(~class, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = method_cols, name = NULL) +
  labs(
    x = NULL,
    y = expression(y^"*"),
    title = expression(paste("Mixed-input predictive distributions at ", x^"*", " = 0"))
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

ggsave(
  file.path(fig_dir, "fig5_mixed_input_prediction_boxplot.pdf"),
  p_pred_box,
  width = 12,
  height = 7.5
)

############################################################
## Figure 6: MCMC traces
############################################################

df_trace <- data.frame(
  draw = seq_len(nrow(samples_u)),
  sigma_epsilon = sqrt(samples_sigma2),
  rho = exp(samples_logtheta[, 1]),
  theta_x = exp(samples_logtheta[, 2]),
  theta_u = exp(samples_logtheta[, 3])
)

df_trace_long <- do.call(
  rbind,
  lapply(names(df_trace)[-1], function(nm) {
    data.frame(
      draw = df_trace$draw,
      parameter = nm,
      value = df_trace[[nm]]
    )
  })
)

p_trace <- ggplot(df_trace_long, aes(x = draw, y = value)) +
  geom_line(linewidth = 0.35) +
  facet_wrap(~parameter, scales = "free_y", ncol = 2) +
  labs(
    x = "Saved draw",
    y = NULL,
    title = "MCMC trace plots"
  )

ggsave(
  file.path(fig_dir, "fig6_mcmc_traces_hyperparameters.pdf"),
  p_trace,
  width = 9,
  height = 6.5
)

df_tau_trace <- do.call(
  rbind,
  lapply(seq_len(m - 1), function(j) {
    data.frame(
      draw = seq_len(nrow(samples_tau)),
      cutpoint = paste0("tau", j),
      value = u_center + u_scale * samples_tau[, j]
    )
  })
)

p_tau_trace <- ggplot(df_tau_trace, aes(x = draw, y = value, color = cutpoint)) +
  geom_line(linewidth = 0.45, alpha = 0.85) +
  geom_hline(
    data = data.frame(
      cutpoint = paste0("tau", seq_len(m - 1)),
      true_value = tau_true_raw
    ),
    aes(yintercept = true_value, color = cutpoint),
    linetype = "dashed"
  ) +
  labs(
    x = "Saved draw",
    y = "Cut point on original u scale",
    color = NULL,
    title = "Cut-point trace plots"
  ) +
  theme(legend.position = "bottom")

ggsave(
  file.path(fig_dir, "fig6_mcmc_traces_cutpoints.pdf"),
  p_tau_trace,
  width = 9,
  height = 4.8
)

############################################################
## 17. Console summary
############################################################

cat("\nPosterior summary:\n")
cat("Mean sigma_epsilon:", mean(sqrt(samples_sigma2)), "\n")
cat("Mean rho:", mean(exp(samples_logtheta[, 1])), "\n")
cat("Mean theta_x:", mean(exp(samples_logtheta[, 2])), "\n")
cat("Mean theta_u:", mean(exp(samples_logtheta[, 3])), "\n")
cat("Latent u 95% coverage for missing cases:",
    mean(coverage_miss), "\n")

cat("\nOutput written to:\n")
cat(normalizePath(out_dir), "\n")
cat("Figures written to:\n")
cat(normalizePath(fig_dir), "\n")

############################################################
## 18. MCMC diagnostics and automatic advice
############################################################

diag_dir <- file.path(out_dir, "diagnostics")
dir.create(diag_dir, showWarnings = FALSE)

############################################################
## Multi-chain split R-hat diagnostics
############################################################

split_rhat <- function(chain_list) {
  ## chain_list is a list of numeric vectors, one per chain.
  chain_list <- lapply(chain_list, function(x) as.numeric(x[is.finite(x)]))
  
  lens <- vapply(chain_list, length, integer(1))
  n0 <- min(lens)
  
  if (length(chain_list) < 2 || n0 < 20) {
    return(NA_real_)
  }
  
  ## Use equal-length tails.
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
  
  m_split <- nrow(split_mat)
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
      lapply(samples_by_chain$tau, function(mat) {
        u_center + u_scale * mat[, j]
      })
    )
  })
)

rhat_u_mis <- data.frame(
  parameter = paste0("u[", miss_idx, "]"),
  global_index = miss_idx,
  rhat = sapply(seq_along(miss_idx), function(k) {
    jj <- miss_idx[k]
    split_rhat(
      lapply(samples_by_chain$u, function(mat) {
        u_center + u_scale * mat[, jj]
      })
    )
  })
)

rhat_key <- rbind(
  transform(rhat_hyper, group = "kernel/noise"),
  transform(rhat_tau, group = "cutpoint")
)

write.csv(
  rhat_hyper,
  file.path(diag_dir, "mcmc_rhat_hyperparameters.csv"),
  row.names = FALSE
)

write.csv(
  rhat_tau,
  file.path(diag_dir, "mcmc_rhat_cutpoints.csv"),
  row.names = FALSE
)

write.csv(
  rhat_u_mis,
  file.path(diag_dir, "mcmc_rhat_missing_u.csv"),
  row.names = FALSE
)

p_rhat_key <- ggplot(rhat_key, aes(x = reorder(parameter, rhat), y = rhat, fill = group)) +
  geom_col(alpha = 0.85) +
  geom_hline(yintercept = 1.01, linetype = "dotted", color = "darkgreen") +
  geom_hline(yintercept = 1.05, linetype = "dashed", color = "firebrick") +
  coord_flip() +
  scale_fill_manual(
    values = c("kernel/noise" = "steelblue", "cutpoint" = "orange"),
    name = NULL
  ) +
  labs(
    x = NULL,
    y = expression("split-" * hat(R)),
    title = expression("Multi-chain split-" * hat(R) * " for key parameters"),
    subtitle = expression("Dotted green: " * hat(R) * "= 1.01" * "; dashed red: " * hat(R) * "= 1.05")
  )

ggsave(
  file.path(diag_dir, "diag_rhat_key_parameters.pdf"),
  p_rhat_key,
  width = 8,
  height = 5.5
)

p_rhat_u <- ggplot(rhat_u_mis, aes(x = rhat)) +
  geom_histogram(bins = 30, fill = "gray70", color = "white") +
  geom_vline(xintercept = 1.01, linetype = "dotted", color = "darkgreen") +
  geom_vline(xintercept = 1.05, linetype = "dashed", color = "firebrick") +
  labs(
    x = expression("split-" * hat(R)),
    y = "Number of missing latent variables",
    title = expression("Split-" * hat(R) * " distribution for missing latent variables")
  )

ggsave(
  file.path(diag_dir, "diag_rhat_missing_u_histogram.pdf"),
  p_rhat_u,
  width = 7,
  height = 4.5
)

cat("\nMulti-chain split-Rhat summary:\n")
cat("Max Rhat among kernel/noise parameters:",
    max(rhat_hyper$rhat, na.rm = TRUE), "\n")
cat("Max Rhat among cut points:",
    max(rhat_tau$rhat, na.rm = TRUE), "\n")
cat("Median Rhat among missing u:",
    median(rhat_u_mis$rhat, na.rm = TRUE), "\n")
cat("Max Rhat among missing u:",
    max(rhat_u_mis$rhat, na.rm = TRUE), "\n")

############################################################
## 18.1 Diagnostic helper functions
############################################################

to_numeric_matrix <- function(x, name = "matrix") {
  if (is.null(dim(x))) {
    x <- matrix(x, ncol = 1)
  } else if (is.data.frame(x)) {
    x <- as.matrix(x)
  } else {
    x <- as.matrix(x)
  }
  
  dn <- dimnames(x)
  
  x_num <- suppressWarnings(matrix(
    as.numeric(x),
    nrow = nrow(x),
    ncol = ncol(x)
  ))
  
  dimnames(x_num) <- dn
  
  if (any(!is.finite(x_num))) {
    warning(sprintf(
      "%s contains %d non-finite values; diagnostics will ignore them.",
      name,
      sum(!is.finite(x_num))
    ))
  }
  
  x_num
}

safe_var <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  
  if (length(x) < 2) return(NA_real_)
  
  out <- tryCatch(
    stats::var(x),
    error = function(e) NA_real_
  )
  
  if (!is.finite(out)) return(NA_real_)
  out
}

safe_sd <- function(x) {
  v <- safe_var(x)
  if (!is.finite(v) || v < 0) return(NA_real_)
  sqrt(v)
}

ess_ips <- function(x, max_lag = NULL) {
  ## Initial-positive-sequence style ESS estimate.
  ## Important: max_lag is capped by default for speed.
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  n <- length(x)
  
  if (n < 5) return(NA_real_)
  
  sx <- safe_sd(x)
  if (!is.finite(sx)) return(NA_real_)
  if (sx == 0) return(n)
  
  if (is.null(max_lag)) {
    ## Capped default. This is much faster than using n/2 for long chains.
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
  
  m <- floor(length(ac) / 2)
  
  if (m >= 1) {
    pair_sums <- ac[2 * seq_len(m) - 1] + ac[2 * seq_len(m)]
    first_nonpos <- which(pair_sums <= 0)[1]
    
    if (is.na(first_nonpos)) {
      use_lag <- 2 * m
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
  ess <- min(max(ess, 1), n)
  
  ess
}

lag1_acf <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  
  if (length(x) < 5) return(NA_real_)
  
  sx <- safe_sd(x)
  
  if (!is.finite(sx)) return(NA_real_)
  if (sx == 0) return(0)
  
  ac <- tryCatch(
    as.numeric(stats::acf(
      x,
      lag.max = 1,
      plot = FALSE,
      demean = TRUE
    )$acf),
    error = function(e) NA_real_
  )
  
  if (length(ac) < 2 || !is.finite(ac[2])) return(NA_real_)
  
  ac[2]
}

geweke_z <- function(x, first = 0.1, last = 0.5, max_lag = 500L) {
  ## Simple Geweke-style z-score using ESS-adjusted variance of segment means.
  ## Uses capped-lag ESS for speed.
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  n <- length(x)
  
  if (n < 50) return(NA_real_)
  
  n_first <- max(10, floor(first * n))
  n_last <- max(10, floor(last * n))
  
  if (n_first + n_last >= n) return(NA_real_)
  
  x_first <- x[seq_len(n_first)]
  x_last <- x[(n - n_last + 1):n]
  
  ess_first <- ess_ips(x_first, max_lag = min(max_lag, n_first - 1L))
  ess_last <- ess_ips(x_last, max_lag = min(max_lag, n_last - 1L))
  
  if (!is.finite(ess_first) || !is.finite(ess_last)) return(NA_real_)
  if (ess_first <= 1 || ess_last <= 1) return(NA_real_)
  
  v_first <- safe_var(x_first)
  v_last <- safe_var(x_last)
  
  if (!is.finite(v_first) || !is.finite(v_last)) return(NA_real_)
  
  var_first_mean <- v_first / ess_first
  var_last_mean <- v_last / ess_last
  
  denom <- sqrt(var_first_mean + var_last_mean)
  
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  
  (mean(x_first) - mean(x_last)) / denom
}

empty_diag_row <- function(parameter,
                           n_total = NA_integer_,
                           n_finite = NA_integer_,
                           error_message = NA_character_) {
  data.frame(
    parameter = parameter,
    n = n_total,
    n_finite = n_finite,
    mean = NA_real_,
    sd = NA_real_,
    q025 = NA_real_,
    q500 = NA_real_,
    q975 = NA_real_,
    ess = NA_real_,
    rel_ess = NA_real_,
    mcse_mean = NA_real_,
    mcse_over_sd = NA_real_,
    acf1 = NA_real_,
    geweke_z = NA_real_,
    diagnostic_error = error_message,
    stringsAsFactors = FALSE
  )
}

diag_one_series <- function(x, parameter) {
  out <- tryCatch({
    
    x_raw <- suppressWarnings(as.numeric(x))
    n_total <- length(x_raw)
    
    x_fin <- x_raw[is.finite(x_raw)]
    n_finite <- length(x_fin)
    
    if (n_finite == 0) {
      return(empty_diag_row(
        parameter = parameter,
        n_total = n_total,
        n_finite = 0,
        error_message = "No finite samples"
      ))
    }
    
    sx <- safe_sd(x_fin)
    
    qs <- tryCatch(
      as.numeric(stats::quantile(
        x_fin,
        probs = c(0.025, 0.5, 0.975),
        na.rm = TRUE,
        names = FALSE
      )),
      error = function(e) c(NA_real_, NA_real_, NA_real_)
    )
    
    ess <- ess_ips(x_fin)
    
    mcse <- if (is.finite(ess) && ess > 0 && is.finite(sx)) {
      sx / sqrt(ess)
    } else {
      NA_real_
    }
    
    data.frame(
      parameter = parameter,
      n = n_total,
      n_finite = n_finite,
      mean = mean(x_fin),
      sd = sx,
      q025 = qs[1],
      q500 = qs[2],
      q975 = qs[3],
      ess = ess,
      rel_ess = ifelse(is.finite(ess), ess / n_finite, NA_real_),
      mcse_mean = mcse,
      mcse_over_sd = ifelse(is.finite(sx) && sx > 0, mcse / sx, NA_real_),
      acf1 = lag1_acf(x_fin),
      geweke_z = geweke_z(x_fin),
      diagnostic_error = NA_character_,
      stringsAsFactors = FALSE
    )
    
  }, error = function(e) {
    empty_diag_row(
      parameter = parameter,
      n_total = length(x),
      n_finite = sum(is.finite(suppressWarnings(as.numeric(x)))),
      error_message = conditionMessage(e)
    )
  })
  
  out
}

make_diag_table <- function(mat, parameter_names = NULL, matrix_name = "matrix") {
  mat <- to_numeric_matrix(mat, name = matrix_name)
  
  if (ncol(mat) == 0) {
    return(data.frame())
  }
  
  if (is.null(parameter_names)) {
    parameter_names <- colnames(mat)
  }
  
  if (is.null(parameter_names)) {
    parameter_names <- paste0("param", seq_len(ncol(mat)))
  }
  
  stopifnot(length(parameter_names) == ncol(mat))
  
  out <- do.call(
    rbind,
    lapply(seq_len(ncol(mat)), function(j) {
      diag_one_series(mat[, j], parameter_names[j])
    })
  )
  
  rownames(out) <- NULL
  out
}

make_diag_table_missing_u_fast <- function(mat,
                                           global_index,
                                           parameter_names = NULL,
                                           matrix_name = "u_mis_mat_raw",
                                           ess_max_lag = 500L,
                                           acf1_max_lag = 1L,
                                           compute_geweke = FALSE,
                                           cores = 1L) {
  mat <- to_numeric_matrix(mat, name = matrix_name)
  
  if (ncol(mat) == 0) {
    return(data.frame())
  }
  
  n_total <- nrow(mat)
  p <- ncol(mat)
  
  if (is.null(parameter_names)) {
    parameter_names <- colnames(mat)
  }
  
  if (is.null(parameter_names)) {
    parameter_names <- paste0("u[", global_index, "]")
  }
  
  stopifnot(length(parameter_names) == p)
  
  ## Convert non-finite values to NA for vectorized summaries.
  mat_na <- mat
  mat_na[!is.finite(mat_na)] <- NA_real_
  
  n_finite <- colSums(is.finite(mat))
  
  means <- matrixStats::colMeans2(mat_na, na.rm = TRUE)
  sds <- matrixStats::colSds(mat_na, na.rm = TRUE)
  
  qs <- matrixStats::colQuantiles(
    mat_na,
    probs = c(0.025, 0.5, 0.975),
    na.rm = TRUE,
    drop = FALSE
  )
  
  ## Column-wise ESS and lag-1 ACF.
  idx <- seq_len(p)
  
  use_parallel <- (
    exists("use_mclapply", inherits = TRUE) &&
      isTRUE(use_mclapply) &&
      cores > 1L &&
      .Platform$OS.type != "windows"
  )
  
  ess_fun <- function(j) {
    ess_ips(mat[, j], max_lag = ess_max_lag)
  }
  
  acf1_fun <- function(j) {
    lag1_acf(mat[, j])
  }
  
  if (use_parallel) {
    ess_vals <- unlist(parallel::mclapply(
      idx,
      ess_fun,
      mc.cores = cores,
      mc.preschedule = TRUE
    ))
    
    acf1_vals <- unlist(parallel::mclapply(
      idx,
      acf1_fun,
      mc.cores = cores,
      mc.preschedule = TRUE
    ))
  } else {
    ess_vals <- vapply(idx, ess_fun, numeric(1))
    acf1_vals <- vapply(idx, acf1_fun, numeric(1))
  }
  
  if (compute_geweke) {
    geweke_fun <- function(j) {
      geweke_z(mat[, j], max_lag = min(ess_max_lag, 500L))
    }
    
    if (use_parallel) {
      geweke_vals <- unlist(parallel::mclapply(
        idx,
        geweke_fun,
        mc.cores = cores,
        mc.preschedule = TRUE
      ))
    } else {
      geweke_vals <- vapply(idx, geweke_fun, numeric(1))
    }
  } else {
    geweke_vals <- rep(NA_real_, p)
  }
  
  mcse_mean <- sds / sqrt(pmax(ess_vals, 1))
  mcse_over_sd <- mcse_mean / sds
  mcse_over_sd[!is.finite(mcse_over_sd)] <- NA_real_
  
  out <- data.frame(
    parameter = parameter_names,
    n = n_total,
    n_finite = n_finite,
    mean = means,
    sd = sds,
    q025 = qs[, 1],
    q500 = qs[, 2],
    q975 = qs[, 3],
    ess = ess_vals,
    rel_ess = ess_vals / pmax(n_finite, 1),
    mcse_mean = mcse_mean,
    mcse_over_sd = mcse_over_sd,
    acf1 = acf1_vals,
    geweke_z = geweke_vals,
    diagnostic_error = NA_character_,
    global_index = global_index,
    stringsAsFactors = FALSE
  )
  
  out
}

make_acf_df <- function(mat, parameter_names = NULL, max_lag = 50,
                        matrix_name = "matrix") {
  mat <- to_numeric_matrix(mat, name = matrix_name)
  
  if (ncol(mat) == 0) {
    return(data.frame())
  }
  
  if (is.null(parameter_names)) {
    parameter_names <- colnames(mat)
  }
  
  if (is.null(parameter_names)) {
    parameter_names <- paste0("param", seq_len(ncol(mat)))
  }
  
  out <- do.call(
    rbind,
    lapply(seq_len(ncol(mat)), function(j) {
      x <- suppressWarnings(as.numeric(mat[, j]))
      x <- x[is.finite(x)]
      
      if (length(x) < 5 || !is.finite(safe_sd(x)) || safe_sd(x) == 0) {
        ac <- c(1, rep(0, max_lag))
      } else {
        ac <- tryCatch(
          as.numeric(stats::acf(
            x,
            lag.max = max_lag,
            plot = FALSE,
            demean = TRUE
          )$acf),
          error = function(e) c(1, rep(NA_real_, max_lag))
        )
      }
      
      data.frame(
        parameter = parameter_names[j],
        lag = seq_along(ac) - 1,
        acf = ac,
        stringsAsFactors = FALSE
      )
    })
  )
  
  rownames(out) <- NULL
  out
}

make_running_mean_df <- function(mat, parameter_names = NULL,
                                 matrix_name = "matrix") {
  mat <- to_numeric_matrix(mat, name = matrix_name)
  
  if (ncol(mat) == 0) {
    return(data.frame())
  }
  
  if (is.null(parameter_names)) {
    parameter_names <- colnames(mat)
  }
  
  if (is.null(parameter_names)) {
    parameter_names <- paste0("param", seq_len(ncol(mat)))
  }
  
  out <- do.call(
    rbind,
    lapply(seq_len(ncol(mat)), function(j) {
      x <- suppressWarnings(as.numeric(mat[, j]))
      finite <- is.finite(x)
      
      cs <- cumsum(ifelse(finite, x, 0))
      cn <- cumsum(finite)
      
      running_mean <- cs / pmax(cn, 1)
      running_mean[cn == 0] <- NA_real_
      
      data.frame(
        parameter = parameter_names[j],
        draw = seq_along(x),
        running_mean = running_mean,
        stringsAsFactors = FALSE
      )
    })
  )
  
  rownames(out) <- NULL
  out
}

############################################################
## 18.2 Build diagnostic matrices
############################################################

hyper_mat <- cbind(
  sigma_epsilon = sqrt(samples_sigma2),
  rho = exp(samples_logtheta[, 1]),
  theta_x = exp(samples_logtheta[, 2]),
  theta_u = exp(samples_logtheta[, 3])
)

hyper_mat <- to_numeric_matrix(hyper_mat, "hyper_mat")

tau_mat_raw <- u_center + u_scale * samples_tau
tau_mat_raw <- to_numeric_matrix(tau_mat_raw, "tau_mat_raw")
colnames(tau_mat_raw) <- paste0("tau", seq_len(ncol(tau_mat_raw)))

u_mis_mat_raw <- u_center + u_scale * samples_u[, miss_idx, drop = FALSE]
u_mis_mat_raw <- to_numeric_matrix(u_mis_mat_raw, "u_mis_mat_raw")
colnames(u_mis_mat_raw) <- paste0("u[", miss_idx, "]")

cat("\nDiagnostic matrix dimensions:\n")
cat("hyper_mat:     ", dim(hyper_mat), "\n")
cat("tau_mat_raw:   ", dim(tau_mat_raw), "\n")
cat("u_mis_mat_raw: ", dim(u_mis_mat_raw), "\n")

cat("\nNon-finite counts:\n")
cat("hyper_mat:     ", sum(!is.finite(hyper_mat)), "\n")
cat("tau_mat_raw:   ", sum(!is.finite(tau_mat_raw)), "\n")
cat("u_mis_mat_raw: ", sum(!is.finite(u_mis_mat_raw)), "\n")

diag_hyper <- make_diag_table(
  hyper_mat,
  matrix_name = "hyper_mat"
)

diag_tau <- make_diag_table(
  tau_mat_raw,
  matrix_name = "tau_mat_raw"
)

diag_cores <- if (
  exists("use_mclapply") &&
  isTRUE(use_mclapply) &&
  exists("mc_cores")
) {
  min(mc_cores, 8L)
} else {
  1L
}

cat("\nComputing fast diagnostics for missing latent variables...\n")
cat("Using ESS max lag = 500 and cores =", diag_cores, "\n")

diag_u_mis <- make_diag_table_missing_u_fast(
  mat = u_mis_mat_raw,
  global_index = miss_idx,
  matrix_name = "u_mis_mat_raw",
  ess_max_lag = 500L,
  compute_geweke = FALSE,
  cores = diag_cores
)

diag_hyper$group <- "kernel/noise"
diag_tau$group <- "cutpoint"
diag_u_mis$group <- "missing latent u"

diag_key <- rbind(
  diag_hyper,
  diag_tau
)
############################################################
## 18.3 Diagnostic plots
############################################################

## ESS bar plot for key parameters
p_ess_key <- ggplot(
  diag_key,
  aes(x = reorder(parameter, ess), y = ess, fill = group)
) +
  geom_col(alpha = 0.85) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "firebrick") +
  geom_hline(yintercept = 400, linetype = "dotted", color = "darkgreen") +
  coord_flip() +
  scale_fill_manual(
    values = c("kernel/noise" = "steelblue", "cutpoint" = "orange"),
    name = NULL
  ) +
  labs(
    x = NULL,
    y = "Effective sample size",
    title = "MCMC effective sample sizes for key parameters",
    subtitle = "Dashed red: ESS = 100; dotted green: ESS = 400"
  ) +
  theme(legend.position = "bottom")

ggsave(
  file.path(diag_dir, "diag_ess_key_parameters.pdf"),
  p_ess_key,
  width = 8,
  height = 5.5
)

## ACF plot for key parameters
acf_key_df <- make_acf_df(
  cbind(hyper_mat, tau_mat_raw),
  max_lag = min(60, floor(nrow(samples_u) / 4))
)

p_acf_key <- ggplot(acf_key_df, aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0, color = "gray40") +
  geom_col(width = 0.75, fill = "steelblue", alpha = 0.8) +
  facet_wrap(~parameter, scales = "free_y", ncol = 3) +
  labs(
    x = "Lag",
    y = "Autocorrelation",
    title = "Autocorrelation functions for key MCMC parameters"
  )

ggsave(
  file.path(diag_dir, "diag_acf_key_parameters.pdf"),
  p_acf_key,
  width = 10,
  height = 7
)

## Running means for key parameters
running_key_df <- make_running_mean_df(
  cbind(hyper_mat, tau_mat_raw)
)

p_running_key <- ggplot(running_key_df, aes(x = draw, y = running_mean)) +
  geom_line(linewidth = 0.45, color = "steelblue") +
  facet_wrap(~parameter, scales = "free_y", ncol = 3) +
  labs(
    x = "Saved draw",
    y = "Running mean",
    title = "Running posterior means for key parameters"
  )

ggsave(
  file.path(diag_dir, "diag_running_means_key_parameters.pdf"),
  p_running_key,
  width = 10,
  height = 7
)

## ESS distribution for missing latent variables
p_u_ess_hist <- ggplot(diag_u_mis, aes(x = ess)) +
  geom_histogram(bins = 30, fill = "gray70", color = "white") +
  geom_vline(xintercept = 100, linetype = "dashed", color = "firebrick") +
  geom_vline(xintercept = 400, linetype = "dotted", color = "darkgreen") +
  labs(
    x = "Effective sample size",
    y = "Number of missing latent variables",
    title = "ESS distribution for missing latent variables",
    subtitle = "Dashed red: ESS = 100; dotted green: ESS = 400"
  )

ggsave(
  file.path(diag_dir, "diag_ess_missing_u_histogram.pdf"),
  p_u_ess_hist,
  width = 7,
  height = 4.5
)

## Trace plots for worst-mixing missing u's
n_worst_u <- min(6, nrow(diag_u_mis))
worst_u <- diag_u_mis[order(diag_u_mis$ess), ][seq_len(n_worst_u), ]

df_worst_u_trace <- do.call(
  rbind,
  lapply(seq_len(nrow(worst_u)), function(k) {
    idx <- worst_u$global_index[k]
    col_id <- which(miss_idx == idx)
    
    data.frame(
      draw = seq_len(nrow(samples_u)),
      parameter = paste0("u[", idx, "]"),
      value = u_mis_mat_raw[, col_id],
      true_value = u_true_raw[idx]
    )
  })
)

df_worst_u_hline <- unique(df_worst_u_trace[, c("parameter", "true_value")])

p_worst_u_trace <- ggplot(df_worst_u_trace, aes(x = draw, y = value)) +
  geom_line(linewidth = 0.35, color = "steelblue") +
  geom_hline(
    data = df_worst_u_hline,
    aes(yintercept = true_value),
    color = "firebrick",
    linetype = "dashed"
  ) +
  facet_wrap(~parameter, scales = "free_y", ncol = 2) +
  labs(
    x = "Saved draw",
    y = "Latent u on original scale",
    title = "Trace plots for the worst-mixing missing latent variables",
    subtitle = "Dashed red lines show true latent values"
  )

ggsave(
  file.path(diag_dir, "diag_worst_missing_u_traces.pdf"),
  p_worst_u_trace,
  width = 9,
  height = 7
)

## ACF plots for worst-mixing missing u's
worst_u_mat <- u_mis_mat_raw[, match(worst_u$global_index, miss_idx), drop = FALSE]
colnames(worst_u_mat) <- paste0("u[", worst_u$global_index, "]")

acf_worst_u_df <- make_acf_df(
  worst_u_mat,
  max_lag = min(60, floor(nrow(samples_u) / 4))
)

p_worst_u_acf <- ggplot(acf_worst_u_df, aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0, color = "gray40") +
  geom_col(width = 0.75, fill = "steelblue", alpha = 0.8) +
  facet_wrap(~parameter, scales = "free_y", ncol = 2) +
  labs(
    x = "Lag",
    y = "Autocorrelation",
    title = "ACF plots for the worst-mixing missing latent variables"
  )

ggsave(
  file.path(diag_dir, "diag_worst_missing_u_acf.pdf"),
  p_worst_u_acf,
  width = 9,
  height = 7
)

############################################################
## 18.4 Automatic mixing assessment and advice
############################################################

make_mcmc_advice <- function(diag_hyper, diag_tau, diag_u_mis,
                             n_iter, burn, thin, control,
                             block_ess_eval_total,
                             block_ess_total,
                             global_ess_eval_total,
                             global_ess_total,
                             theta_update_total,
                             theta_eval_total) {
  lines <- character(0)
  
  add <- function(...) {
    lines <<- c(lines, sprintf(...))
  }
  
  n_saved <- unique(diag_hyper$n)[1]
  
  key <- rbind(diag_hyper, diag_tau)
  
  min_key_ess <- min(key$ess, na.rm = TRUE)
  min_hyper_ess <- min(diag_hyper$ess, na.rm = TRUE)
  min_tau_ess <- min(diag_tau$ess, na.rm = TRUE)
  
  median_u_ess <- median(diag_u_mis$ess, na.rm = TRUE)
  min_u_ess <- min(diag_u_mis$ess, na.rm = TRUE)
  frac_u_ess_lt_100 <- mean(diag_u_mis$ess < 100, na.rm = TRUE)
  frac_u_ess_lt_200 <- mean(diag_u_mis$ess < 200, na.rm = TRUE)
  
  max_key_acf1 <- max(abs(key$acf1), na.rm = TRUE)
  max_hyper_acf1 <- max(abs(diag_hyper$acf1), na.rm = TRUE)
  max_tau_acf1 <- max(abs(diag_tau$acf1), na.rm = TRUE)
  
  key_geweke_flag <- key$parameter[is.finite(key$geweke_z) & abs(key$geweke_z) > 2.5]
  hyper_geweke_flag <- diag_hyper$parameter[
    is.finite(diag_hyper$geweke_z) & abs(diag_hyper$geweke_z) > 2.5
  ]
  tau_geweke_flag <- diag_tau$parameter[
    is.finite(diag_tau$geweke_z) & abs(diag_tau$geweke_z) > 2.5
  ]
  
  avg_block_eval <- block_ess_eval_total / max(block_ess_total, 1)
  avg_global_eval <- global_ess_eval_total / max(global_ess_total, 1)
  avg_theta_eval <- theta_eval_total / max(theta_update_total, 1)
  
  add("MCMC diagnostic report")
  add("======================")
  add("")
  add("This report is based on a single MCMC chain.")
  add("Single-chain diagnostics can reveal poor mixing, but they cannot prove convergence.")
  add("For final reported results, it is advisable to run at least 3 dispersed chains and compare posterior summaries.")
  add("")
  add("Run settings")
  add("------------")
  add("Preset: %s", control$preset)
  add("Total iterations: %d", n_iter)
  add("Burn-in: %d", burn)
  add("Thinning interval: %d", thin)
  add("Saved posterior draws: %d", n_saved)
  add("Theta update frequency: every %d iteration(s)", control$theta_update_every)
  add("Block ESS block size: %d", control$ess_block_size)
  add("Local z-slice updates per iteration: %d", control$local_per_iter)
  add("")
  add("Mixing summary")
  add("--------------")
  add("Minimum ESS among kernel/noise parameters: %.1f", min_hyper_ess)
  add("Minimum ESS among cut points: %.1f", min_tau_ess)
  add("Minimum ESS among all key parameters: %.1f", min_key_ess)
  add("Median ESS among missing latent variables: %.1f", median_u_ess)
  add("Minimum ESS among missing latent variables: %.1f", min_u_ess)
  add("Fraction of missing latent variables with ESS < 100: %.2f", frac_u_ess_lt_100)
  add("Fraction of missing latent variables with ESS < 200: %.2f", frac_u_ess_lt_200)
  add("Maximum |lag-1 ACF| among key parameters: %.3f", max_key_acf1)
  add("Maximum |lag-1 ACF| among kernel/noise parameters: %.3f", max_hyper_acf1)
  add("Maximum |lag-1 ACF| among cut points: %.3f", max_tau_acf1)
  add("")
  add("Computational behavior")
  add("----------------------")
  add("Average likelihood evaluations per block ESS update: %.1f", avg_block_eval)
  if (global_ess_total > 0) {
    add("Average likelihood evaluations per global ESS update: %.1f", avg_global_eval)
  } else {
    add("No global ESS updates were performed.")
  }
  add("Average likelihood evaluations per theta-slice update: %.1f", avg_theta_eval)
  add("")
  
  add("Geweke-style checks")
  add("-------------------")
  if (length(key_geweke_flag) == 0) {
    add("No key parameters have |Geweke z| > 2.5.")
  } else {
    add("Parameters with |Geweke z| > 2.5: %s",
        paste(key_geweke_flag, collapse = ", "))
  }
  add("")
  
  add("Automatic advice")
  add("----------------")
  
  good_key <- is.finite(min_key_ess) && min_key_ess >= 400
  good_u <- is.finite(median_u_ess) &&
    median_u_ess >= 200 &&
    frac_u_ess_lt_100 <= 0.10
  good_acf <- is.finite(max_key_acf1) && max_key_acf1 < 0.85
  good_geweke <- length(key_geweke_flag) == 0
  
  if (good_key && good_u && good_acf && good_geweke) {
    add("Overall assessment: mixing looks acceptable for this numerical example.")
    add("The current preset '%s' appears adequate.", control$preset)
    add("For publication-quality results, still consider running multiple chains.")
  } else {
    add("Overall assessment: some signs of slow mixing remain.")
    
    ## Suggested longer run length based on key ESS
    target_key_ess <- 400
    target_u_median_ess <- 200
    
    factor_key <- if (is.finite(min_key_ess) && min_key_ess > 0) {
      target_key_ess / min_key_ess
    } else {
      2
    }
    
    factor_u <- if (is.finite(median_u_ess) && median_u_ess > 0) {
      target_u_median_ess / median_u_ess
    } else {
      2
    }
    
    factor_needed <- max(1, factor_key, factor_u)
    suggested_n_iter <- ceiling(n_iter * min(5, 1.25 * factor_needed))
    suggested_burn <- ceiling(suggested_n_iter / 3)
    
    if (factor_needed > 1.2) {
      add("If computation permits, increase n_iter from %d to about %d and burn to about %d.",
          n_iter, suggested_n_iter, suggested_burn)
    }
    
    if (min_hyper_ess < 200 || max_hyper_acf1 > 0.90 ||
        length(hyper_geweke_flag) > 0) {
      add("Kernel/noise hyperparameters mix slowly.")
      add("Recommended actions:")
      add("  - update GP hyperparameters more often; for example, set theta_update_every to %d.",
          max(1, floor(control$theta_update_every / 2)))
      add("  - run a longer chain, because length-scale parameters often have high posterior dependence.")
      add("  - check the trace and ACF plots in the diagnostics folder.")
    }
    
    if (min_tau_ess < 200 || max_tau_acf1 > 0.90 ||
        length(tau_geweke_flag) > 0) {
      add("Some cut points mix slowly.")
      add("Recommended actions:")
      add("  - use mcmc_preset = 'thorough', or reduce full_local_every to update all latent variables more frequently.")
      add("  - run a longer chain, because cut points are strongly coupled with the latent u_i's.")
    }
    
    if (frac_u_ess_lt_100 > 0.10 || median_u_ess < 200) {
      add("Some missing latent variables mix slowly.")
      add("Recommended actions:")
      add("  - use mcmc_preset = 'thorough'.")
      add("  - increase local_per_iter or update all missing u_i's more frequently.")
      add("  - inspect 'diag_worst_missing_u_traces.pdf' and 'diag_worst_missing_u_acf.pdf'.")
    }
    
    if (avg_block_eval > 40) {
      add("Blocked elliptical slice updates require many likelihood evaluations.")
      add("This usually means the ordinal constraints make large prior-ellipse moves difficult.")
      add("Recommended actions:")
      add("  - reduce ess_block_size slightly.")
      add("  - rely more on local z-slice updates by increasing local_per_iter.")
      add("  - keep occasional global ESS moves, but do not make them too frequent.")
    }
    
    if (avg_block_eval < 8 && frac_u_ess_lt_100 > 0.10) {
      add("Blocked ESS updates are inexpensive but latent ESS is low.")
      add("Recommended action: increase ess_block_size or use the 'thorough' preset.")
    }
    
    if (length(key_geweke_flag) > 0) {
      add("Some key parameters differ between early and late parts of the chain.")
      add("Recommended actions:")
      add("  - increase burn-in.")
      add("  - initialize another chain from an overdispersed starting point.")
      add("  - compare posterior summaries across chains.")
    }
  }
  
  add("")
  add("Interpretation guide")
  add("--------------------")
  add("ESS > 400 is usually comfortable for posterior means and rough intervals.")
  add("ESS between 100 and 400 is often usable but should be interpreted with caution.")
  add("ESS < 100 suggests that more MCMC effort is needed.")
  add("High lag-1 autocorrelation, say above 0.9, indicates slow local movement.")
  add("Large |Geweke z| values are a warning sign, not a formal proof of non-convergence.")
  
  lines
}

diagnostic_report <- make_mcmc_advice(
  diag_hyper = diag_hyper,
  diag_tau = diag_tau,
  diag_u_mis = diag_u_mis,
  n_iter = n_iter,
  burn = burn,
  thin = thin,
  control = control,
  block_ess_eval_total = block_ess_eval_total,
  block_ess_total = block_ess_total,
  global_ess_eval_total = global_ess_eval_total,
  global_ess_total = global_ess_total,
  theta_update_total = theta_update_total,
  theta_eval_total = theta_eval_total
)

writeLines(
  diagnostic_report,
  con = file.path(diag_dir, "mcmc_diagnostic_report.txt")
)

cat("\n")
cat(paste(diagnostic_report, collapse = "\n"))
cat("\n")

############################################################
## 18.5 Add diagnostics to saved RDS object
############################################################

res$diagnostics <- list(
  hyperparameters = diag_hyper,
  cutpoints = diag_tau,
  missing_u = diag_u_mis,
  report = diagnostic_report
)

saveRDS(res, file.path(out_dir, "eivgp_ordinal_results_robust.rds"))

cat("\nDiagnostics written to:\n")
cat(normalizePath(diag_dir), "\n")

