############################################################
## 00_study1_functions.R
##
## EIV-GP sampler: univariate latent u / ordinal thresholds,
## covariates x of any dimension. A numeric vector x is treated
## as n x 1 and recovers the original Study I model.
############################################################

if (!exists("mixedgp_parallel_lapply", mode = "function")) {
  this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  sibling_dir <- if (is.null(this_file)) character(0) else {
    dirname(normalizePath(this_file))
  }
  parallel_utility <- c(
    file.path(sibling_dir, "00_parallel_utils.R"),
    "00_parallel_utils.R",
    file.path("codes", "00_parallel_utils.R")
  )
  parallel_utility <- parallel_utility[file.exists(parallel_utility)][1L]
  if (is.na(parallel_utility)) {
    stop("Source 00_parallel_utils.R before 00_study1_functions.R.")
  }
  sys.source(parallel_utility, envir = environment())
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

as_study1_numeric_matrix <- function(x,
                                     name,
                                     nrow_expected = NULL,
                                     ncol_expected = NULL,
                                     expected_names = NULL) {
  input_names <- if (is.null(dim(x))) NULL else colnames(x)

  if (is.null(dim(x))) {
    if (is.null(ncol_expected) || ncol_expected == 1L) {
      x <- matrix(x, ncol = 1L)
    } else if (length(x) == ncol_expected) {
      x <- matrix(x, nrow = 1L)
    } else {
      stop(
        name, " is a vector, but a matrix with ", ncol_expected,
        " columns is required."
      )
    }
  }

  if (is.data.frame(x)) {
    numeric_columns <- vapply(
      x,
      function(z) is.numeric(z) || is.integer(z) || is.logical(z),
      logical(1)
    )
    if (any(!numeric_columns)) {
      stop(name, " must contain only numeric columns; encode factors first.")
    }
  } else if (!(is.numeric(x) || is.integer(x) || is.logical(x))) {
    stop(name, " must be numeric.")
  }

  out <- as.matrix(x)
  storage.mode(out) <- "double"

  if (ncol(out) < 1L) stop(name, " must contain at least one column.")
  mixedgp_validate_column_names(out, name)
  if (!is.null(nrow_expected) && nrow(out) != nrow_expected) {
    stop(name, " must have ", nrow_expected, " rows.")
  }
  if (!is.null(ncol_expected) && ncol(out) != ncol_expected) {
    stop(name, " must have ", ncol_expected, " columns.")
  }
  if (any(!is.finite(out))) stop(name, " contains missing or non-finite values.")

  if (!is.null(expected_names)) {
    if (length(expected_names) != ncol(out)) {
      stop(name, " has the wrong number of columns.")
    }
    if (!is.null(input_names) && all(nzchar(input_names))) {
      if (!setequal(input_names, expected_names)) {
        stop(name, " column names do not match the training predictors.")
      }
      out <- out[, match(expected_names, colnames(out)), drop = FALSE]
    }
    colnames(out) <- expected_names
  }

  out
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

  if (any(!is.finite(mean)) || any(!is.finite(sd)) || any(sd <= 0) ||
      anyNA(lower) || anyNA(upper) || any(lower >= upper)) {
    stop("Invalid parameters in rtruncnorm_vec().")
  }

  a <- (lower - mean) / sd
  b <- (upper - mean) / sd
  z <- numeric(n)

  ## Draw uniformly between two probabilities without subtracting nearly
  ## equal CDF values.  Positive-tail intervals use survival probabilities;
  ## negative-tail intervals use lower-tail probabilities.  qnorm(log.p=TRUE)
  ## then remains accurate even when the ordinary probabilities underflow.
  log_uniform_between <- function(log_lo, log_hi, size) {
    ratio <- exp(log_lo - log_hi)
    log_hi + log(ratio + runif(size) * (1 - ratio))
  }

  right <- a >= 0
  if (any(right)) {
    log_lo <- pnorm(b[right], lower.tail = FALSE, log.p = TRUE)
    log_hi <- pnorm(a[right], lower.tail = FALSE, log.p = TRUE)
    log_prob <- log_uniform_between(log_lo, log_hi, sum(right))
    z[right] <- qnorm(log_prob, lower.tail = FALSE, log.p = TRUE)
  }

  left <- !right & b <= 0
  if (any(left)) {
    log_lo <- pnorm(a[left], log.p = TRUE)
    log_hi <- pnorm(b[left], log.p = TRUE)
    log_prob <- log_uniform_between(log_lo, log_hi, sum(left))
    z[left] <- qnorm(log_prob, log.p = TRUE)
  }

  middle <- !right & !left
  if (any(middle)) {
    p_lo <- pnorm(a[middle])
    p_hi <- pnorm(b[middle])
    prob <- p_lo + runif(sum(middle)) * (p_hi - p_lo)
    z[middle] <- qnorm(prob)
  }

  out <- mean + sd * z
  tolerance <- 64 * .Machine$double.eps * pmax(1, abs(lower), abs(upper))
  violates <- !is.finite(out) |
    out < lower - tolerance | out > upper + tolerance
  if (any(violates)) {
    stop("Numerical failure in tail-stable truncated-normal sampling.")
  }
  pmin(pmax(out, lower), upper)
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

solve_chol <- function(U, b) {
  backsolve(U, forwardsolve(t(U), b))
}

rmvnorm_psd <- function(n, mean, Sigma, tolerance = NULL) {
  mean <- as.numeric(mean)
  Sigma <- as.matrix(Sigma)
  if (!identical(dim(Sigma), c(length(mean), length(mean)))) {
    stop("Sigma has incompatible dimensions in rmvnorm_psd().")
  }
  Sigma <- 0.5 * (Sigma + t(Sigma))
  eig <- eigen(Sigma, symmetric = TRUE)
  scale <- max(1, max(abs(eig$values)))
  if (is.null(tolerance)) {
    tolerance <- 100 * length(mean) * .Machine$double.eps * scale
  }
  if (min(eig$values) < -tolerance) {
    stop("Predictive covariance is not positive semidefinite within tolerance.")
  }
  values <- pmax(eig$values, 0)
  factor <- diag(sqrt(values), nrow = length(values)) %*% t(eig$vectors)
  Z <- matrix(rnorm(as.integer(n) * length(mean)), as.integer(n), length(mean))
  out <- sweep(Z %*% factor, 2, mean, "+")
  attr(out, "psd_tolerance") <- tolerance
  attr(out, "clipped_eigenvalues") <- sum(eig$values < 0)
  out
}

bounded_slice_update_1d <- function(x0, logf, w = 1,
                                 lower = -Inf, upper = Inf,
                                 max_steps_out = 50,
                                 max_iter = 200,
                                 fail_on_limit = FALSE) {
  if (length(w) != 1L || !is.finite(w) || w <= 0 ||
      max_steps_out < 1L || max_iter < 1L) {
    stop("Slice width and search budgets must be positive.")
  }
  if (upper <= lower) {
    if (isTRUE(fail_on_limit)) stop("Slice interval is empty.")
    return(list(x = x0, n_eval = 0))
  }
  
  if (length(x0) != 1L || !is.finite(x0) || x0 < lower || x0 > upper) {
    stop("Current slice state is outside its declared interval.")
  }
  
  f0 <- logf(x0)
  n_eval <- 1L
  
  if (!is.finite(f0)) {
    if (isTRUE(fail_on_limit)) stop("Current slice state has non-finite density.")
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
  
  if (isTRUE(fail_on_limit)) {
    stop("Slice shrinkage exceeded max_iter.")
  }
  list(x = x0, n_eval = n_eval)
}

full_interval_slice_update <- function(x0, logf, lower, upper,
                                       max_iter = 200,
                                       fail_on_limit = FALSE) {
  if (!is.finite(lower) || !is.finite(upper)) {
    stop("full_interval_slice_update requires finite lower and upper bounds.")
  }
  
  if (upper <= lower) {
    if (isTRUE(fail_on_limit)) stop("Full slice interval is empty.")
    return(list(x = x0, n_eval = 0))
  }
  
  if (length(x0) != 1L || !is.finite(x0) || x0 < lower || x0 > upper) {
    stop("Current full-interval slice state is outside its declared interval.")
  }
  
  f0 <- logf(x0)
  n_eval <- 1L
  
  if (!is.finite(f0)) {
    if (isTRUE(fail_on_limit)) {
      stop("Current full-interval slice state has non-finite density.")
    }
    return(list(x = x0, n_eval = n_eval))
  }
  
  logy <- f0 + log(runif(1))
  L <- lower
  R <- upper
  
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
  
  if (isTRUE(fail_on_limit)) {
    stop("Full-interval slice shrinkage exceeded max_iter.")
  }
  list(x = x0, n_eval = n_eval)
}

.study1_category_moment_cache <- new.env(parent = emptyenv())

study1_category_trig_moments <- function(tau,
                                         frequency = 0.9 * base::pi) {
  tau <- as.numeric(tau)
  if (length(tau) < 1L || any(!is.finite(tau)) ||
      is.unsorted(tau, strictly = TRUE)) {
    stop("tau must be a finite, strictly increasing vector.")
  }

  cache_key <- paste(
    formatC(c(frequency, tau), digits = 16L, format = "fg"),
    collapse = "|"
  )
  if (exists(cache_key, envir = .study1_category_moment_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .study1_category_moment_cache, inherits = FALSE))
  }

  bounds <- c(-Inf, tau, Inf)
  m <- length(bounds) - 1L
  probability <- diff(pnorm(bounds))
  if (any(probability <= 0)) stop("tau creates an empty ordinal category.")

  integrate_moment <- function(fun, lower, upper) {
    stats::integrate(
      function(z) fun(frequency * z) * dnorm(z),
      lower = lower,
      upper = upper,
      rel.tol = 1e-11,
      subdivisions = 200L,
      stop.on.error = TRUE
    )$value
  }
  cosine <- sine <- numeric(m)
  for (r in seq_len(m)) {
    cosine[r] <- integrate_moment(cos, bounds[r], bounds[r + 1L]) /
      probability[r]
    sine[r] <- integrate_moment(sin, bounds[r], bounds[r + 1L]) /
      probability[r]
  }

  out <- list(cosine = cosine, sine = sine, probability = probability)
  assign(cache_key, out, envir = .study1_category_moment_cache)
  out
}

f0_1d <- function(x,
                  u,
                  scenario = c(
                    "active", "inactive", "category_sufficient",
                    "heterogeneity_continuum"
                  ),
                  c_ord = NULL,
                  tau = c(-5 / 3, -10 / 9, 0, 10 / 9, 5 / 3),
                  heterogeneity_eta = 1) {
  scenario <- match.arg(scenario)
  
  if (scenario == "inactive") {
    return(cos(0.9 * base::pi * u))
  }

  x_main <- 0.5 * sin(base::pi * x / 2)
  if (scenario %in% c("category_sufficient", "heterogeneity_continuum")) {
    if (is.null(c_ord)) c_ord <- make_class(u, tau)
    c_ord <- as.integer(c_ord)
    moments <- study1_category_trig_moments(tau)
    if (length(c_ord) != max(length(x), length(u)) ||
        anyNA(c_ord) || any(c_ord < 1L | c_ord > length(moments$cosine))) {
      stop("c_ord must contain one valid ordinal level per response input.")
    }
    conditional_mean <- moments$cosine[c_ord] +
      0.35 * (x / 2) * moments$sine[c_ord]
    if (scenario == "category_sufficient") return(x_main + conditional_mean)

    heterogeneity_eta <- as.numeric(heterogeneity_eta)
    if (length(heterogeneity_eta) != 1L || !is.finite(heterogeneity_eta) ||
        heterogeneity_eta < 0) {
      stop("heterogeneity_eta must be one finite nonnegative number.")
    }
    continuous_effect <- cos(0.9 * base::pi * u) +
      0.35 * (x / 2) * sin(0.9 * base::pi * u)
    return(
      x_main + conditional_mean +
        heterogeneity_eta * (continuous_effect - conditional_mean)
    )
  }
  
  x_main + cos(0.9 * base::pi * u) +
    0.35 * (x / 2) * sin(0.9 * base::pi * u)
}

m0_1d <- function(x,
                  c_ord,
                  tau = c(-5 / 3, -10 / 9, 0, 10 / 9, 5 / 3),
                  scenario = c(
                    "active", "inactive", "category_sufficient",
                    "heterogeneity_continuum"
                  )) {
  scenario <- match.arg(scenario)
  x <- as.numeric(x)
  c_ord <- as.integer(c_ord)
  n_out <- max(length(x), length(c_ord))
  x <- rep(x, length.out = n_out)
  c_ord <- rep(c_ord, length.out = n_out)
  moments <- study1_category_trig_moments(tau)
  if (anyNA(c_ord) || any(c_ord < 1L | c_ord > length(moments$cosine))) {
    stop("c_ord contains invalid ordinal levels in m0_1d().")
  }
  conditional_continuous <- moments$cosine[c_ord] +
    0.35 * (x / 2) * moments$sine[c_ord]
  if (scenario == "inactive") return(moments$cosine[c_ord])
  0.5 * sin(base::pi * x / 2) + conditional_continuous
}

simulate_1d_data <- function(n = 100,
                             n_test = 500,
                             m = 6,
                             scenario = c(
                               "active", "inactive", "category_sufficient",
                               "heterogeneity_continuum"
                             ),
                             sigma_eps = 0.1,
                             seed = NULL,
                             threshold_design = c("imbalanced", "balanced"),
                             min_class_count = NULL,
                             heterogeneity_eta = 1) {
  scenario <- match.arg(scenario)
  threshold_design <- match.arg(threshold_design)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (threshold_design == "imbalanced") {
    if (m != 6L) stop("The imbalanced Study I design is defined for m = 6.")
    tau_true <- c(-5 / 3, -10 / 9, 0, 10 / 9, 5 / 3)
    if (is.null(min_class_count)) min_class_count <- 3L
  } else {
    tau_true <- qnorm(seq_len(m - 1L) / m)
    if (is.null(min_class_count)) min_class_count <- 0L
  }
  min_class_count <- as.integer(min_class_count)
  if (length(min_class_count) != 1L || is.na(min_class_count) ||
      min_class_count < 0L) {
    stop("min_class_count must be one nonnegative integer.")
  }
  
  repeat {
    x <- maximin_lhs_1d(n, lower = -2, upper = 2)
    u <- rnorm(n)
    c_ord <- make_class(u, tau_true)
    
    if (min_class_count == 0L ||
        all(tabulate(c_ord, nbins = m) >= min_class_count)) break
  }
  
  f <- f0_1d(
    x, u, scenario = scenario, c_ord = c_ord, tau = tau_true,
    heterogeneity_eta = heterogeneity_eta
  )
  y <- f + rnorm(n, mean = 0, sd = sigma_eps)
  
  x_test <- runif(n_test, -2, 2)
  u_test <- rnorm(n_test)
  c_test <- make_class(u_test, tau_true)
  f_test <- f0_1d(
    x_test, u_test, scenario = scenario, c_ord = c_test, tau = tau_true,
    heterogeneity_eta = heterogeneity_eta
  )
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
    m = m,
    threshold_design = threshold_design,
    min_class_count = min_class_count,
    heterogeneity_eta = heterogeneity_eta
  )
}

make_nested_calibration_sets <- function(n, calib_grid, seed = NULL) {
  n <- as.integer(n)
  calib_grid <- as.integer(calib_grid)
  if (length(n) != 1L || is.na(n) || n < 1L ||
      length(calib_grid) < 1L || anyNA(calib_grid) ||
      any(calib_grid < 0L | calib_grid > n)) {
    stop("Calibration sizes must be integers between zero and n.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  ord <- sample(seq_len(n))
  
  out <- lapply(calib_grid, function(k) {
    if (k == 0L) integer(0) else sort(ord[seq_len(k)])
  })
  
  names(out) <- as.character(calib_grid)
  out
}

make_theta_spec_multix <- function(p_x) {
  if (!is.finite(p_x) || p_x < 1L) stop("p_x must be positive")
  parameter <- c("log_rho", paste0("log_theta_x", seq_len(p_x)), "log_theta_u")
  bounds <- rbind(
    c(log(0.05), log(100)),
    matrix(rep(c(log(1e-4), log(100)), p_x + 1L), ncol = 2, byrow = TRUE)
  )
  rownames(bounds) <- parameter
  list(
    p_x = p_x,
    parameter = parameter,
    prior_mean = c(log(3), rep(log(0.5), p_x + 1L)),
    prior_sd = rep(1.5, p_x + 2L),
    bounds = bounds,
    x_index = 1L + seq_len(p_x),
    u_index = p_x + 2L
  )
}

weighted_distance_multix <- function(Dlist, theta) {
  if (length(Dlist) != length(theta)) stop("Distance and theta dimensions differ")
  out <- matrix(0, nrow(Dlist[[1]]), ncol(Dlist[[1]]))
  for (j in seq_along(Dlist)) out <- out + theta[j] * Dlist[[j]]
  out
}

new_gp_cache_1d <- function(use_blocks = TRUE) {
  cache <- new.env(parent = emptyenv())
  cache$entry <- NULL
  cache$use_blocks <- isTRUE(use_blocks)
  cache$full_factorizations <- 0L
  cache$cache_hits <- 0L
  cache$block_setups <- 0L
  cache$block_evaluations <- 0L
  cache$block_fallbacks <- 0L
  class(cache) <- "eivgp_cache_1d"
  cache
}

gp_cache_matches_1d <- function(cache, y, u, Dx_list, logtheta,
                                theta_spec, kernel, matern_nu) {
  if (is.null(cache)) return(FALSE)
  if (!is.environment(cache) || !inherits(cache, "eivgp_cache_1d")) {
    stop("gp_cache must be created by new_gp_cache_1d().")
  }
  entry <- cache$entry
  !is.null(entry) && identical(entry$y, y) && identical(entry$u, u) &&
    identical(entry$Dx_list, Dx_list) && identical(entry$logtheta, logtheta) &&
    identical(entry$theta_spec, theta_spec) && identical(entry$kernel, kernel) &&
    identical(entry$matern_nu, matern_nu)
}

gp_cache_store_1d <- function(cache, state, y, u, Dx_list, logtheta,
                              theta_spec, kernel, matern_nu) {
  if (!is.null(cache)) {
    cache$entry <- list(state = state, y = y, u = u, Dx_list = Dx_list,
                        logtheta = logtheta, theta_spec = theta_spec,
                        kernel = kernel, matern_nu = matern_nu)
  }
  invisible(state)
}

gp_update_loglik_1d <- function(state, n, collapse_noise = FALSE,
                                a_eps0 = 2, b_eps0 = 0.05) {
  if (!collapse_noise) return(state$loglik)
  ## A = I + rho^2 R and sigma2 ~ IG(a_eps0, b_eps0). Integrating
  ## sigma2 removes neither the latent-input prior nor the theta prior.
  ## The determinant and quadratic form are shared by dense/Schur caches.
  -0.5 * (n * log(2 * base::pi) + state$logdetA) +
    a_eps0 * log(b_eps0) - lgamma(a_eps0) + lgamma(a_eps0 + n / 2) -
    (a_eps0 + n / 2) * log(b_eps0 + 0.5 * state$quad)
}

gp_state_1d <- function(y, u, Dx_list, logtheta, sigma2_eps,
                        theta_spec = NULL,
                        kernel = "se",
                        matern_nu = 2.5,
                        gp_cache = NULL,
                        require_factor = TRUE) {
  if (!is.list(Dx_list)) Dx_list <- list(as.matrix(Dx_list))
  if (is.null(theta_spec)) {
    theta_spec <- make_theta_spec_multix(length(Dx_list))
  }
  n <- length(y)
  if (gp_cache_matches_1d(gp_cache, y, u, Dx_list, logtheta,
                          theta_spec, kernel, matern_nu) &&
      (!isTRUE(require_factor) || !is.null(gp_cache$entry$state$cholA))) {
    state <- gp_cache$entry$state
    gp_cache$cache_hits <- gp_cache$cache_hits + 1L
    state$loglik <- -0.5 * (n * log(2 * base::pi * sigma2_eps) +
                             state$logdetA + state$quad / sigma2_eps)
    return(state)
  }

  rho <- exp(logtheta[1])
  theta_x <- exp(logtheta[theta_spec$x_index])
  theta_u <- exp(logtheta[theta_spec$u_index])

  Du <- pairwise_sqdist(matrix(u, ncol = 1))
  D2 <- weighted_distance_multix(Dx_list, theta_x) + theta_u * Du
  R <- kernel_from_weighted_sqdist_1d(D2, kernel, matern_nu)
  
  A <- rho^2 * R + diag(n)
  
  U <- safe_chol_1d(A)
  if (!is.null(gp_cache)) {
    gp_cache$full_factorizations <- gp_cache$full_factorizations + 1L
  }
  Ainv_y <- solve_chol(U, y)
  
  logdetA <- 2 * sum(log(diag(U)))
  quad <- sum(y * Ainv_y)
  
  loglik <- -0.5 * (
    n * log(2 * base::pi * sigma2_eps) +
      logdetA +
      quad / sigma2_eps
  )
  
  state <- list(
    loglik = loglik,
    R = R,
    A = A,
    cholA = U,
    Ainv_y = Ainv_y,
    logdetA = logdetA,
    quad = quad
  )
  gp_cache_store_1d(gp_cache, state, y, u, Dx_list, logtheta,
                    theta_spec, kernel, matern_nu)
  state
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

gp_block_state_1d <- function(prepared, u_block, sigma2_eps, gp_cache) {
  a <- prepared
  if (length(u_block) != a$b || any(!is.finite(u_block))) {
    stop("A GP block proposal has incompatible or nonfinite latent values.")
  }
  u_prop <- a$u
  u_prop[a$block_idx] <- u_block
  if (gp_cache_matches_1d(gp_cache, a$y, u_prop, a$Dx_list, a$logtheta,
                          a$theta_spec, a$kernel, a$matern_nu)) {
    return(gp_state_1d(
      a$y, u_prop, a$Dx_list, a$logtheta, sigma2_eps, a$theta_spec,
      a$kernel, a$matern_nu, gp_cache = gp_cache, require_factor = FALSE
    ))
  }
  gp_cache$block_evaluations <- gp_cache$block_evaluations + 1L
  ## Only the proposed rows/columns change. Factor the fixed complement once,
  ## then evaluate the exact Gaussian conditional through its Schur complement.
  candidate <- tryCatch({
    cross <- a$rho2 * kernel_from_weighted_sqdist_1d(
      a$Dx_cross + a$theta_u * pairwise_sqdist(matrix(u_block, ncol = 1),
                                                matrix(a$u[a$rest], ncol = 1)),
      a$kernel, a$matern_nu
    )
    block <- diag(a$b) + a$rho2 * kernel_from_weighted_sqdist_1d(
      a$Dx_block + a$theta_u * pairwise_sqdist(matrix(u_block, ncol = 1)),
      a$kernel, a$matern_nu
    )
    V <- forwardsolve(t(a$chol_rest), t(cross))
    schur <- block - crossprod(V)
    schur <- 0.5 * (schur + t(schur))
    H <- chol(schur)
    residual <- a$y[a$block_idx] - as.numeric(cross %*% a$alpha_rest)
    quad <- a$quad_rest + sum(residual * solve_chol(H, residual))
    logdet <- a$logdet_rest + 2 * sum(log(diag(H)))
    if (!is.finite(quad) || quad < 0 || !is.finite(logdet)) {
      stop("Nonfinite Schur-complement likelihood calculation.")
    }
    A <- a$A
    A[a$block_idx, a$rest] <- cross
    A[a$rest, a$block_idx] <- t(cross)
    A[a$block_idx, a$block_idx] <- block
    list(A = A, cholA = NULL, R = NULL, Ainv_y = NULL,
         quad = quad, logdetA = logdet,
         loglik = -0.5 * (a$n * log(2 * base::pi * sigma2_eps) +
                           logdet + quad / sigma2_eps))
  }, error = function(e) NULL)
  if (is.null(candidate)) {
    ## Retry the same covariance with the dense factorization. This does not
    ## regularize the matrix or reject a proposal because the shortcut failed.
    gp_cache$block_fallbacks <- gp_cache$block_fallbacks + 1L
    return(gp_state_1d(
      a$y, u_prop, a$Dx_list, a$logtheta, sigma2_eps, a$theta_spec,
      a$kernel, a$matern_nu, gp_cache = gp_cache, require_factor = TRUE
    ))
  }
  gp_cache_store_1d(gp_cache, candidate, a$y, u_prop, a$Dx_list, a$logtheta,
                    a$theta_spec, a$kernel, a$matern_nu)
  candidate
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

make_default_control <- function(n, n_mis,
                                 preset = c("fast", "balanced", "thorough"),
                                 p_x = 1L) {
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
    theta_slice_width_init = rep(1.0, p_x + 2L),
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

log_prior_logtheta <- function(logtheta, theta_spec = NULL,
                               mean = NULL, sd = NULL, bounds = NULL) {
  if (is.null(theta_spec)) {
    p_x <- length(logtheta) - 2L
    if (!is.finite(p_x) || p_x < 1L) {
      stop("logtheta must have length p_x + 2")
    }
    theta_spec <- make_theta_spec_multix(p_x)
  }
  if (!is.null(mean)) theta_spec$prior_mean <- mean
  if (!is.null(sd)) theta_spec$prior_sd <- sd
  if (!is.null(bounds)) theta_spec$bounds <- bounds
  if (any(logtheta < theta_spec$bounds[, 1]) ||
      any(logtheta > theta_spec$bounds[, 2])) {
    return(-Inf)
  }
  
  sum(dnorm(
    logtheta,
    mean = theta_spec$prior_mean,
    sd = theta_spec$prior_sd,
    log = TRUE
  ))
}

gp_loglik_with_constraints_1d <- function(y, u, Dx_list, c_ord, tau,
                                          logtheta, sigma2_eps,
                                          theta_spec = NULL,
                                          kernel = "se",
                                          matern_nu = 2.5,
                                          gp_cache = NULL,
                                          collapse_noise = FALSE) {
  if (!is.list(Dx_list)) Dx_list <- list(as.matrix(Dx_list))
  if (is.null(theta_spec)) {
    theta_spec <- make_theta_spec_multix(length(Dx_list))
  }
  if (!check_constraints_1d(u, c_ord, tau)) {
    return(-Inf)
  }
  
  state <- gp_state_1d(
    y, u, Dx_list, logtheta, if (collapse_noise) 1 else sigma2_eps, theta_spec,
    kernel = kernel, matern_nu = matern_nu,
    gp_cache = gp_cache, require_factor = FALSE
  )
  gp_update_loglik_1d(state, length(y), collapse_noise)
}

theta_logpost_1d <- function(y, u, Dx_list, logtheta, sigma2_eps,
                             theta_spec = NULL,
                             prior_mean = NULL, prior_sd = NULL, bounds = NULL,
                             kernel = "se", matern_nu = 2.5,
                             gp_cache = NULL,
                             collapse_noise = FALSE) {
  if (!is.list(Dx_list)) Dx_list <- list(as.matrix(Dx_list))
  if (is.null(theta_spec)) {
    theta_spec <- make_theta_spec_multix(length(Dx_list))
  }
  lp <- log_prior_logtheta(
    logtheta, theta_spec,
    mean = prior_mean, sd = prior_sd, bounds = bounds
  )
  
  if (!is.finite(lp)) {
    return(-Inf)
  }
  
  state <- gp_state_1d(
    y, u, Dx_list, logtheta, if (collapse_noise) 1 else sigma2_eps, theta_spec,
    kernel = kernel, matern_nu = matern_nu,
    gp_cache = gp_cache, require_factor = FALSE
  )
  gp_update_loglik_1d(state, length(y), collapse_noise) + lp
}

initialize_tau_1d <- function(c_ord, u_obs, calib_idx, m, tau_bound = 8) {
  n <- length(c_ord)
  
  counts <- tabulate(c_ord, nbins = m)
  probs <- cumsum(counts)[1:(m - 1)] / n
  probs <- pmin(pmax(probs, 0.03), 0.97)
  
  tau <- qnorm(probs)
  obs_c <- c_ord[calib_idx]
  obs_u <- u_obs[calib_idx]
  for (j in seq_len(m - 1)) {
    lower <- max(c(-tau_bound, if (j > 1L) tau[j - 1L],
                   obs_u[obs_c <= j]))
    upper <- min(c(tau_bound, obs_u[obs_c > j]))
    if (lower >= upper) {
      stop("Calibrated u values are incompatible with ordinal labels or threshold bounds.")
    }
    ## An absolute inward margin can exceed a valid calibration gap. A
    ## midpoint also leaves room for later thresholds sharing this gap.
    if (!is.finite(tau[j]) || tau[j] <= lower || tau[j] >= upper) {
      tau[j] <- lower + (upper - lower) / 2
    }
    if (tau[j] <= lower || tau[j] >= upper) {
      stop("No representable interior value exists in a threshold initialization interval.")
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
    
    if (!is.finite(L) || !is.finite(U) || L >= U) {
      stop("Threshold Gibbs update encountered an empty conditional interval.")
    }
    tau_new[j] <- runif(1, L, U)
  }
  
  tau_new
}

assert_threshold_state_1d <- function(u,
                                      tau,
                                      c_ord,
                                      m,
                                      calib_idx = integer(0),
                                      u_obs = NULL,
                                      logtheta = NULL,
                                      sigma2_eps = NULL,
                                      tolerance = 1e-10) {
  if (length(u) != length(c_ord) || any(!is.finite(u)) ||
      length(tau) != m - 1L || any(!is.finite(tau)) ||
      is.unsorted(tau, strictly = TRUE)) {
    stop("The deterministic-threshold sampler state is invalid.")
  }
  lower <- c(-Inf, tau)[c_ord]
  upper <- c(tau, Inf)[c_ord]
  if (any(u <= lower | u > upper)) {
    stop("A latent input left the interval implied by its ordinal category.")
  }
  if (length(calib_idx) > 0L) {
    if (is.null(u_obs) || anyNA(u_obs[calib_idx]) ||
        any(u[calib_idx] != u_obs[calib_idx])) {
      stop("A calibrated latent input moved away from its observed value.")
    }
  }
  if (!is.null(logtheta) && any(!is.finite(logtheta))) {
    stop("A GP hyperparameter is nonfinite.")
  }
  if (!is.null(sigma2_eps) &&
      (length(sigma2_eps) != 1L || !is.finite(sigma2_eps) || sigma2_eps <= 0)) {
    stop("The observation-noise variance is not positive and finite.")
  }
  invisible(TRUE)
}

sample_sigma2_eps_1d <- function(y, u, Dx_list, logtheta, theta_spec = NULL,
                                 a_eps0 = 2, b_eps0 = 0.05,
                                 kernel = "se", matern_nu = 2.5,
                                 gp_cache = NULL) {
  if (!is.list(Dx_list)) Dx_list <- list(as.matrix(Dx_list))
  if (is.null(theta_spec)) {
    theta_spec <- make_theta_spec_multix(length(Dx_list))
  }
  n <- length(y)
  
  ## The factor, determinant and quadratic form do not depend on sigma2.
  quad <- gp_state_1d(
    y, u, Dx_list, logtheta, sigma2_eps = 1, theta_spec = theta_spec,
    kernel = kernel, matern_nu = matern_nu,
    gp_cache = gp_cache, require_factor = FALSE
  )$quad
  
  shape <- a_eps0 + n / 2
  rate <- b_eps0 + 0.5 * quad
  
  1 / rgamma(1, shape = shape, rate = rate)
}

update_u_ess_block_1d <- function(y, u, Dx_list, c_ord, tau, logtheta,
                                      sigma2_eps, theta_spec,
                                      block_idx, max_try = 300L,
                                      kernel = "se", matern_nu = 2.5,
                                      gp_cache = NULL,
                                      collapse_noise = FALSE) {
  if (length(block_idx) == 0) {
    return(list(u = u, n_eval = 0L, accepted = TRUE))
  }
  
  u_block <- u[block_idx]
  block_state <- prepare_gp_block_1d(
    y, u, Dx_list, logtheta, theta_spec, block_idx,
    kernel, matern_nu, gp_cache
  )
  
  loglik_fun <- function(u_block_prop) {
    u_prop <- u
    u_prop[block_idx] <- u_block_prop
    if (!is.null(block_state)) {
      if (!check_constraints_1d(u_prop, c_ord, tau)) return(-Inf)
      candidate <- gp_block_state_1d(
        block_state, u_block_prop, if (collapse_noise) 1 else sigma2_eps, gp_cache
      )
      return(gp_update_loglik_1d(candidate, length(y), collapse_noise))
    }
    
    gp_loglik_with_constraints_1d(
      y = y,
      u = u_prop,
      Dx_list = Dx_list,
      c_ord = c_ord,
      tau = tau,
      logtheta = logtheta,
      sigma2_eps = sigma2_eps,
      theta_spec = theta_spec,
      kernel = kernel,
      matern_nu = matern_nu,
      gp_cache = gp_cache,
      collapse_noise = collapse_noise
    )
  }
  
  loglik_cur <- loglik_fun(u_block)
  
  if (!is.finite(loglik_cur)) {
    stop("Current latent state violates constraints or has invalid likelihood.")
  }
  
  nu <- rnorm(length(u_block))
  
  logy <- loglik_cur + log(runif(1))
  
  angle <- runif(1, 0, 2 * base::pi)
  angle_min <- angle - 2 * base::pi
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
  
  stop("Elliptical-slice update exceeded max_try.")
}

update_u_local_z_slice_1d <- function(y, u, Dx_list, c_ord, tau,
                                          logtheta, sigma2_eps, theta_spec,
                                          update_idx,
                                          kernel = "se", matern_nu = 2.5,
                                          gp_cache = NULL,
                                          collapse_noise = FALSE) {
  n_eval_total <- 0L
  update_idx <- update_idx[sample.int(length(update_idx))]
  
  for (i in update_idx) {
    lower_u <- c(-Inf, tau)[c_ord[i]]
    upper_u <- c(tau, Inf)[c_ord[i]]
    block_state <- prepare_gp_block_1d(
      y, u, Dx_list, logtheta, theta_spec, i,
      kernel, matern_nu, gp_cache
    )
    local_gp_loglik <- function(value) {
      if (!is.null(block_state)) {
        candidate <- gp_block_state_1d(
          block_state, value, if (collapse_noise) 1 else sigma2_eps, gp_cache
        )
        return(gp_update_loglik_1d(candidate, length(y), collapse_noise))
      }
      u_prop <- u
      u_prop[i] <- value
      candidate <- gp_state_1d(
        y, u_prop, Dx_list, logtheta, if (collapse_noise) 1 else sigma2_eps, theta_spec,
        kernel = kernel, matern_nu = matern_nu,
        gp_cache = gp_cache, require_factor = FALSE
      )
      gp_update_loglik_1d(candidate, length(y), collapse_noise)
    }
    
    z_lower <- pnorm(lower_u)
    z_upper <- pnorm(upper_u)
    
    ## Choose the coordinate system using only the fixed category bounds,
    ## not the current latent value. Infinite/tail/narrow intervals use an
    ## ordinary u-scale slice including the normal prior. This retains the
    ## whole interval instead of clipping its CDF tails or projecting u.
    probability_scale <- is.finite(lower_u) && is.finite(upper_u) &&
      z_lower > 1e-8 && z_upper < 1 - 1e-8 && z_upper - z_lower > 1e-8
    if (!probability_scale) {
      logf_u <- function(value) {
        if (!is.finite(value) || value <= lower_u || value > upper_u) return(-Inf)
        local_gp_loglik(value) + dnorm(value, log = TRUE)
      }
      width <- if (is.finite(lower_u) && is.finite(upper_u)) {
        min(1, (upper_u - lower_u) / 2)
      } else 1
      ans <- bounded_slice_update_1d(
        x0 = u[i], logf = logf_u, w = width,
        lower = lower_u, upper = upper_u, max_steps_out = 30L,
        max_iter = 200L, fail_on_limit = TRUE
      )
      n_eval_total <- n_eval_total + ans$n_eval
      u[i] <- ans$x
      next
    }
    z0 <- pnorm(u[i])
    
    logf_z <- function(z) {
      if (z <= z_lower || z > z_upper) {
        return(-Inf)
      }
      
      u_prop_i <- qnorm(z)
      
      if (!(u_prop_i > lower_u && u_prop_i <= upper_u)) {
        return(-Inf)
      }
      
      local_gp_loglik(u_prop_i)
    }
    
    ans <- full_interval_slice_update(
      x0 = z0,
      logf = logf_z,
      lower = z_lower,
      upper = z_upper,
      max_iter = 200L,
      fail_on_limit = TRUE
    )
    
    n_eval_total <- n_eval_total + ans$n_eval
    u[i] <- qnorm(ans$x)
  }
  
  list(u = u, n_eval = n_eval_total)
}

update_logtheta_slice_1d <- function(y, u, Dx_list, logtheta, sigma2_eps,
                                         theta_slice_width, theta_spec,
                                         kernel = "se", matern_nu = 2.5,
                                         gp_cache = NULL,
                                         collapse_noise = FALSE) {
  n_eval_total <- 0L
  
  for (j in seq_along(logtheta)) {
    logf_j <- function(val) {
      lt <- logtheta
      lt[j] <- val
      
      theta_logpost_1d(
        y, u, Dx_list, lt, sigma2_eps, theta_spec,
        kernel = kernel, matern_nu = matern_nu, gp_cache = gp_cache,
        collapse_noise = collapse_noise
      )
    }
    
    ans <- bounded_slice_update_1d(
      x0 = logtheta[j],
      logf = logf_j,
      w = theta_slice_width[j],
      lower = theta_spec$bounds[j, 1],
      upper = theta_spec$bounds[j, 2],
      max_steps_out = 30L,
      max_iter = 100L,
      fail_on_limit = TRUE
    )
    
    logtheta[j] <- ans$x
    n_eval_total <- n_eval_total + ans$n_eval
  }
  
  list(logtheta = logtheta, n_eval = n_eval_total)
}

split_rhat_1d <- function(chain_list) {
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

aligned_chain_matrix <- function(chain_list) {
  chain_list <- lapply(chain_list, as.numeric)
  lens <- vapply(chain_list, length, integer(1))
  n0 <- min(lens)
  if (length(chain_list) < 2L || n0 < 4L) return(NULL)
  do.call(cbind, lapply(chain_list, function(x) tail(x, n0)))
}

rank_rhat_1d <- function(chain_list) {
  mat <- aligned_chain_matrix(chain_list)
  if (is.null(mat)) return(NA_real_)
  if (requireNamespace("posterior", quietly = TRUE)) {
    return(as.numeric(posterior::rhat(mat)))
  }
  split_rhat_1d(chain_list)
}

bulk_ess_1d <- function(chain_list) {
  mat <- aligned_chain_matrix(chain_list)
  if (is.null(mat)) return(NA_real_)
  if (requireNamespace("posterior", quietly = TRUE)) {
    return(as.numeric(posterior::ess_bulk(mat)))
  }
  sum(vapply(chain_list, ess_ips_1d, numeric(1)))
}

tail_ess <- function(chain_list) {
  mat <- aligned_chain_matrix(chain_list)
  if (is.null(mat)) return(NA_real_)
  if (requireNamespace("posterior", quietly = TRUE)) {
    return(as.numeric(posterior::ess_tail(mat)))
  }
  NA_real_
}

ess_ips_1d <- function(x, max_lag = NULL) {
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

fit_eivgp_1d <- function(...) {
  if (!exists("fit_eivgp_1d_v030", mode = "function")) {
    stop("Load the current eivGP package or source codes/load_mixedgp.R before fitting.")
  }
  fit_eivgp_1d_v030(...)
}

gp_mle_fit_1d <- function(X, y, kernel = c("se", "matern"), matern_nu = 2.5,
                       n_starts = 5L, seed = 1L) {
  X <- as_study1_numeric_matrix(X, "X", nrow_expected = length(y))
  if (any(!is.finite(y)) || length(y) < 2L) {
    stop("y must contain at least two finite observations.")
  }
  y <- as.numeric(y)
  n <- nrow(X)
  d <- ncol(X)
  kernel_spec <- normalize_gp_kernel_1d(kernel, matern_nu)
  
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
    
    R <- kernel_from_weighted_sqdist_1d(
      Rexp, kernel_spec$name, kernel_spec$matern_nu
    )
    A <- rho^2 * R + diag(n)
    
    U <- try(safe_chol_1d(A), silent = TRUE)
    
    if (inherits(U, "try-error")) {
      return(1e20)
    }
    
    Ainv_y <- solve_chol(U, y)
    logdetA <- 2 * sum(log(diag(U)))
    quad <- sum(y * Ainv_y)
    
    0.5 * (n * log(2 * base::pi * sigma2) + logdetA + quad / sigma2)
  }
  
  lower <- c(log(1e-5), log(0.05), rep(log(1e-4), d))
  upper <- c(log(5), log(100), rep(log(100), d))
  n_starts <- as.integer(n_starts)
  if (length(n_starts) != 1L || is.na(n_starts) || n_starts < 1L) {
    stop("n_starts must be one positive integer.")
  }
  set.seed(as.integer(seed))
  initial <- matrix(NA_real_, nrow = n_starts, ncol = length(lower))
  initial[1L, ] <- c(log(0.05), log(3), rep(log(0.5), d))
  if (n_starts > 1L) {
    for (ss in 2:n_starts) {
      initial[ss, ] <- stats::runif(length(lower), lower, upper)
    }
  }

  attempts <- lapply(seq_len(n_starts), function(ss) {
    tryCatch(
      optim(
        par = initial[ss, ], fn = nll, method = "L-BFGS-B",
        lower = lower, upper = upper, control = list(maxit = 500)
      ),
      error = function(e) e
    )
  })
  report <- data.frame(
    start = seq_len(n_starts), convergence = NA_integer_, objective = Inf,
    message = "", stringsAsFactors = FALSE
  )
  for (ss in seq_along(attempts)) {
    if (inherits(attempts[[ss]], "error")) {
      report$message[ss] <- conditionMessage(attempts[[ss]])
    } else {
      report$convergence[ss] <- attempts[[ss]]$convergence
      report$objective[ss] <- attempts[[ss]]$value
      report$message[ss] <- if (is.null(attempts[[ss]]$message)) {
        ""
      } else {
        attempts[[ss]]$message
      }
    }
  }
  valid <- which(
    report$convergence == 0L & is.finite(report$objective)
  )
  if (length(valid) == 0L) {
    stop("Internal GP MLE did not converge from any prespecified start.")
  }
  best <- valid[which.min(report$objective[valid])]
  opt <- attempts[[best]]
  
  list(
    par = opt$par,
    value = opt$value,
    convergence = opt$convergence,
    optimizer_attempts = report,
    selected_start = best,
    X = X,
    y = y,
    kernel = kernel_spec
  )
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

gp_mle_fit_learned_embedding <- function(x, c_ord, y, m,
                                         n_starts = 8L,
                                         kernel = c("se", "matern"),
                                         matern_nu = 2.5) {
  n <- length(y)
  x <- as_study1_numeric_matrix(x, "x", nrow_expected = n)
  p_x <- ncol(x)
  kernel_spec <- normalize_gp_kernel_1d(kernel, matern_nu)
  Dx_list <- lapply(seq_len(p_x), function(j) {
    pairwise_sqdist(x[, j, drop = FALSE])
  })
  idx_x <- 2L + seq_len(p_x)
  idx_z <- p_x + 3L
  idx_a <- if (m > 2L) seq.int(p_x + 4L, p_x + m + 1L) else integer(0)
  
  nll <- function(par) {
    log_sigma2 <- par[1]
    log_rho <- par[2]
    log_theta_x <- par[idx_x]
    log_theta_z <- par[idx_z]

    a <- par[idx_a]
    
    sigma2 <- exp(log_sigma2)
    rho <- exp(log_rho)
    theta_x <- exp(log_theta_x)
    theta_z <- exp(log_theta_z)
    
    z_scores <- make_monotone_scores_1d(a, m)
    z <- z_scores[c_ord]
    
    Dz <- pairwise_sqdist(matrix(z, ncol = 1))
    
    D2 <- weighted_distance_multix(Dx_list, theta_x) + theta_z * Dz
    R <- kernel_from_weighted_sqdist_1d(
      D2, kernel_spec$name, kernel_spec$matern_nu
    )
    A <- rho^2 * R + diag(n)
    
    U <- try(safe_chol_1d(A), silent = TRUE)
    
    if (inherits(U, "try-error")) {
      return(1e20)
    }
    
    Ainv_y <- solve_chol(U, y)
    logdetA <- 2 * sum(log(diag(U)))
    quad <- sum(y * Ainv_y)
    
    0.5 * (n * log(2 * base::pi * sigma2) + logdetA + quad / sigma2)
  }
  
  lower <- c(
    log(1e-5),
    log(0.05),
    rep(log(1e-4), p_x),
    log(1e-4),
    rep(-6, max(m - 2, 0))
  )
  
  upper <- c(
    log(5),
    log(100),
    rep(log(100), p_x),
    log(100),
    rep(6, max(m - 2, 0))
  )
  
  make_init <- function() {
    c(
      log(0.05),
      log(3),
      rep(log(0.5), p_x),
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
  
  z_scores <- make_monotone_scores_1d(opt$par[idx_a], m)
  
  list(
    par = opt$par[seq_len(idx_z)],
    a = opt$par[idx_a],
    z_scores = z_scores,
    value = opt$value,
    convergence = opt$convergence,
    x = x,
    c_ord = c_ord,
    y = y,
    m = m,
    kernel = kernel_spec
  )
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

fit_embedding_baselines <- function(x_raw, y_raw, c_ord, m,
                                    n_starts_learned = 8L,
                                    kernel = c("se", "matern"),
                                    matern_nu = 2.5) {
  x_raw <- as_study1_numeric_matrix(
    x_raw, "x_raw", nrow_expected = length(y_raw)
  )
  kernel_spec <- normalize_gp_kernel_1d(kernel, matern_nu)
  x_center <- colMeans(x_raw)
  x_scale <- apply(x_raw, 2, sd)
  if (any(!is.finite(x_scale) | x_scale <= 0)) {
    stop("Every x_raw column must have positive finite variation")
  }
  x <- sweep(sweep(x_raw, 2, x_center, "-"), 2, x_scale, "/")
  
  y_center <- mean(y_raw)
  y_scale <- sd(y_raw)
  if (!is.finite(y_scale) || y_scale <= 0) {
    stop("y_raw must have positive finite variation.")
  }
  y <- as.numeric((y_raw - y_center) / y_scale)
  
  z_gauss <- qnorm(c_ord / (m + 1))
  
  fit_gauss <- gp_mle_fit_1d(
    X = cbind(x, z_gauss),
    y = y,
    kernel = kernel_spec$name,
    matern_nu = kernel_spec$matern_nu
  )
  
  z_cm_scores <- conditional_mean_scores(c_ord, m)
  z_cm <- z_cm_scores[c_ord]
  
  fit_cm <- gp_mle_fit_1d(
    X = cbind(x, z_cm),
    y = y,
    kernel = kernel_spec$name,
    matern_nu = kernel_spec$matern_nu
  )
  
  fit_learned <- gp_mle_fit_learned_embedding(
    x = x,
    c_ord = c_ord,
    y = y,
    m = m,
    n_starts = n_starts_learned,
    kernel = kernel_spec$name,
    matern_nu = kernel_spec$matern_nu
  )
  
  list(
    x_center = x_center,
    x_scale = x_scale,
    y_center = y_center,
    y_scale = y_scale,
    predictor_names = colnames(x_raw),
    kernel = kernel_spec,
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

normal_mixture_nlpd <- function(y_true,
                                component_means,
                                component_vars) {
  y_true <- as.numeric(y_true)
  if (is.null(component_means) || is.null(component_vars)) {
    missing <- c(
      if (is.null(component_means)) "conditional_means" else NULL,
      if (is.null(component_vars)) "conditional_vars" else NULL
    )
    return(list(
      value = NA_real_,
      pointwise = rep(NA_real_, length(y_true)),
      reason = paste0(
        "predictive normal-mixture components unavailable: ",
        paste(missing, collapse = " and ")
      )
    ))
  }

  component_means <- as.matrix(component_means)
  component_vars <- as.matrix(component_vars)
  if (ncol(component_means) != length(y_true) ||
      ncol(component_vars) != length(y_true)) {
    stop("Predictive mixture components must have one column per test case.")
  }
  n_component <- max(nrow(component_means), nrow(component_vars))
  if (nrow(component_means) == 1L && n_component > 1L) {
    component_means <- component_means[
      rep(1L, n_component), , drop = FALSE
    ]
  }
  if (nrow(component_vars) == 1L && n_component > 1L) {
    component_vars <- component_vars[
      rep(1L, n_component), , drop = FALSE
    ]
  }
  if (nrow(component_means) != n_component ||
      nrow(component_vars) != n_component || n_component < 1L) {
    stop("Predictive component means and variances have incompatible rows.")
  }
  if (any(!is.finite(component_means)) ||
      any(!is.finite(component_vars)) || any(component_vars < 0)) {
    stop("Predictive normal-mixture components must have finite, nonnegative variances.")
  }
  if (any(component_vars == 0)) {
    return(list(
      value = NA_real_,
      pointwise = rep(NA_real_, length(y_true)),
      reason = paste0(
        "continuous NLPD undefined: predictive mixture contains ",
        "zero-variance components"
      )
    ))
  }

  pointwise <- vapply(seq_along(y_true), function(j) {
    log_component <- stats::dnorm(
      y_true[j],
      mean = component_means[, j],
      sd = sqrt(component_vars[, j]),
      log = TRUE
    )
    largest <- max(log_component)
    if (is.infinite(largest) && largest < 0) return(Inf)
    -(largest + log(mean(exp(log_component - largest))))
  }, numeric(1L))
  list(value = mean(pointwise), pointwise = pointwise, reason = "")
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

summarize_mean_recovery_1d <- function(draw_mat,
                                       m_true,
                                       method,
                                       rep_id,
                                       n_calib,
                                       scenario,
                                       valid_function_draws = FALSE) {
  draw_mat <- as.matrix(draw_mat)
  m_true <- as.numeric(m_true)
  if (ncol(draw_mat) != length(m_true) || any(!is.finite(m_true))) {
    stop("draw_mat and m_true are incompatible in summarize_mean_recovery_1d().")
  }
  conditional_means <- attr(draw_mat, "conditional_means")
  post_mean <- if (!is.null(conditional_means)) {
    colMeans(conditional_means)
  } else {
    colMeans(draw_mat)
  }
  lo <- hi <- rep(NA_real_, length(m_true))
  if (isTRUE(valid_function_draws)) {
    qs <- apply(
      draw_mat, 2L, stats::quantile,
      probs = c(0.025, 0.975), names = FALSE
    )
    lo <- qs[1L, ]
    hi <- qs[2L, ]
  }
  error <- post_mean - m_true
  data.frame(
    rep = rep_id,
    scenario = scenario,
    n_calib = n_calib,
    method = method,
    n_eval = length(m_true),
    RMSE = sqrt(mean(error^2)),
    MAE = mean(abs(error)),
    Bias = mean(error),
    Coverage95 = if (isTRUE(valid_function_draws)) {
      mean(m_true >= lo & m_true <= hi)
    } else {
      NA_real_
    },
    Width95 = if (isTRUE(valid_function_draws)) mean(hi - lo) else NA_real_,
    uncertainty_target = if (isTRUE(valid_function_draws)) "m(x,c)" else "none",
    stringsAsFactors = FALSE
  )
}

sample_oracle_y_star <- function(x_star, c_star, tau_true,
                                 scenario = "active",
                                 sigma_eps = 0.1,
                                 n_draw = 1000,
                                 heterogeneity_eta = 1) {
  lower <- c(-Inf, tau_true)[c_star]
  upper <- c(tau_true, Inf)[c_star]
  
  u_star <- rtruncnorm_vec(
    mean = rep(0, n_draw),
    sd = 1,
    lower = rep(lower, n_draw),
    upper = rep(upper, n_draw)
  )
  
  f_star <- f0_1d(
    x_star,
    u_star,
    scenario = scenario,
    c_ord = rep(c_star, n_draw),
    tau = tau_true,
    heterogeneity_eta = heterogeneity_eta
  )
  
  out <- f_star + rnorm(n_draw, 0, sigma_eps)
  ## Conditional on each Monte Carlo draw of U*, the DGM predictive law is
  ## exactly Gaussian.  Retain those components so that log predictive density
  ## can be evaluated without applying a KDE to the simulated Y* draws.
  attr(out, "conditional_means") <- matrix(f_star, ncol = 1L)
  attr(out, "conditional_vars") <- matrix(
    rep(sigma_eps^2, n_draw), ncol = 1L
  )
  out
}

sample_oracle_test_y <- function(x_test,
                                 c_test,
                                 tau_true,
                                 scenario,
                                 sigma_eps,
                                 n_draw = 1000,
                                 heterogeneity_eta = 1) {
  n_test <- length(c_test)
  out <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  conditional_means <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  conditional_vars <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  
  for (j in seq_len(n_test)) {
    draw_j <- sample_oracle_y_star(
      x_star = x_test[j],
      c_star = c_test[j],
      tau_true = tau_true,
      scenario = scenario,
      sigma_eps = sigma_eps,
      n_draw = n_draw,
      heterogeneity_eta = heterogeneity_eta
    )
    out[, j] <- draw_j
    conditional_means[, j] <- attr(draw_j, "conditional_means")[, 1L]
    conditional_vars[, j] <- attr(draw_j, "conditional_vars")[, 1L]
  }

  attr(out, "conditional_means") <- conditional_means
  attr(out, "conditional_vars") <- conditional_vars
  attr(out, "mixture_components") <- "DGM conditional on Monte Carlo U*"
  out
}

sample_eiv_test_y <- function(x_test_raw,
                              c_test,
                              fit_obj,
                              draw_ids,
                              n_per_draw = 1L,
                              u_test_obs = NULL,
                              u_input_scale = c("raw", "model"),
                              joint = FALSE) {
  u_input_scale <- match.arg(u_input_scale)
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
  predictor_names <- fit_obj$data$predictor_names
  if (is.null(predictor_names)) predictor_names <- colnames(x_train)
  kernel_spec <- kernel_spec_from_fit_1d(fit_obj)

  x_test_raw <- as_study1_numeric_matrix(
    x_test_raw,
    "x_test_raw",
    ncol_expected = ncol(x_train),
    expected_names = predictor_names
  )
  x_test <- sweep(sweep(x_test_raw, 2, x_center, "-"), 2, x_scale, "/")

  n_test <- nrow(x_test)
  if (!(is.numeric(c_test) || is.integer(c_test)) ||
      length(c_test) != n_test || anyNA(c_test) ||
      any(c_test != as.integer(c_test)) ||
      any(c_test < 1L | c_test > fit_obj$data$m)) {
    stop("c_test must contain one ordinal code 1, ..., m per prediction row.")
  }
  c_test <- as.integer(c_test)
  observed_u <- integer(0)
  u_test_model <- rep(NA_real_, n_test)
  if (!is.null(u_test_obs)) {
    u_test_obs <- as.numeric(u_test_obs)
    if (length(u_test_obs) != n_test || any(is.infinite(u_test_obs))) {
      stop("u_test_obs must contain one finite value or NA per prediction row.")
    }
    observed_u <- which(is.finite(u_test_obs))
    u_center <- if (is.null(fit_obj$data$U_center)) 0 else fit_obj$data$U_center
    u_scale <- if (is.null(fit_obj$data$U_scale)) 1 else fit_obj$data$U_scale
    u_test_model[observed_u] <- if (u_input_scale == "raw") {
      (u_test_obs[observed_u] - u_center) / u_scale
    } else {
      u_test_obs[observed_u]
    }
  }
  draw_ids <- as.integer(draw_ids)
  if (length(draw_ids) < 1L || anyNA(draw_ids) ||
      any(draw_ids < 1L | draw_ids > nrow(samples_u))) {
    stop("draw_ids must index saved posterior draws.")
  }
  n_per_draw <- as.integer(n_per_draw)
  if (length(n_per_draw) != 1L || is.na(n_per_draw) || n_per_draw < 1L) {
    stop("n_per_draw must be one positive integer.")
  }
  n_draw <- length(draw_ids) * n_per_draw
  
  out <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  conditional_means <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  conditional_vars <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  
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
      if (length(observed_u) > 0L) u_star[observed_u] <- u_test_model[observed_u]
      
      pred <- gp_predict_draw(
        x_train = x_train,
        u_train = samples_u[s, ],
        y_train = y_train,
        x_star = x_test,
        u_star = u_star,
        logtheta = samples_logtheta[s, ],
        sigma2_eps = samples_sigma2[s],
        noisy = TRUE,
        kernel = kernel_spec$name,
        matern_nu = kernel_spec$matern_nu,
        return_cov = isTRUE(joint)
      )

      y_std <- if (isTRUE(joint)) {
        as.numeric(rmvnorm_psd(1L, pred$mean, pred$cov))
      } else {
        pred$mean + sqrt(pred$var) * rnorm(n_test)
      }
      out[row_id, ] <- y_center + y_scale * y_std
      conditional_means[row_id, ] <- y_center + y_scale * pred$mean
      conditional_vars[row_id, ] <- y_scale^2 * pred$var
    }
  }

  attr(out, "joint") <- isTRUE(joint)
  attr(out, "latent_input_scale") <- u_input_scale
  attr(out, "conditional_means") <- conditional_means
  attr(out, "conditional_vars") <- conditional_vars
  attr(out, "mixture_components") <-
    "posterior and latent-input Monte Carlo Gaussian components"
  out
}

sample_eiv_f_given_xu_1d <- function(x_star_raw,
                                      u_star,
                                      fit_obj,
                                      draw_ids,
                                      n_per_draw = 1L,
                                      include_gp_uncertainty = TRUE,
                                      u_input_scale = c("raw", "model"),
                                      joint = FALSE) {
  u_input_scale <- match.arg(u_input_scale)
  samples_u <- fit_obj$mcmc$samples_u
  samples_logtheta <- fit_obj$mcmc$samples_logtheta
  samples_sigma2 <- fit_obj$mcmc$samples_sigma2
  x_train <- fit_obj$data$x
  predictor_names <- fit_obj$data$predictor_names
  if (is.null(predictor_names)) predictor_names <- colnames(x_train)
  kernel_spec <- kernel_spec_from_fit_1d(fit_obj)

  x_star_raw <- as_study1_numeric_matrix(
    x_star_raw,
    "x_star_raw",
    ncol_expected = ncol(x_train),
    expected_names = predictor_names
  )
  u_star <- as.numeric(u_star)
  if (length(u_star) != nrow(x_star_raw) || any(!is.finite(u_star))) {
    stop("u_star must contain one finite value per prediction row.")
  }
  u_center <- if (is.null(fit_obj$data$U_center)) 0 else fit_obj$data$U_center
  u_scale <- if (is.null(fit_obj$data$U_scale)) 1 else fit_obj$data$U_scale
  if (u_input_scale == "raw") {
    u_star <- (u_star - u_center) / u_scale
  }
  x_star <- sweep(
    sweep(x_star_raw, 2L, fit_obj$data$x_center, "-"),
    2L, fit_obj$data$x_scale, "/"
  )
  draw_ids <- as.integer(draw_ids)
  if (length(draw_ids) < 1L || anyNA(draw_ids) ||
      any(draw_ids < 1L | draw_ids > nrow(samples_u))) {
    stop("draw_ids must index saved posterior draws.")
  }
  n_per_draw <- as.integer(n_per_draw)
  if (length(n_per_draw) != 1L || is.na(n_per_draw) || n_per_draw < 1L) {
    stop("n_per_draw must be one positive integer.")
  }

  out <- matrix(
    NA_real_,
    nrow = length(draw_ids) * n_per_draw,
    ncol = nrow(x_star)
  )
  row_id <- 0L
  for (s in draw_ids) {
    pred <- gp_predict_draw(
      x_train = x_train,
      u_train = samples_u[s, ],
      y_train = fit_obj$data$y,
      x_star = x_star,
      u_star = u_star,
      logtheta = samples_logtheta[s, ],
      sigma2_eps = samples_sigma2[s],
      noisy = FALSE,
      kernel = kernel_spec$name,
      matern_nu = kernel_spec$matern_nu,
      return_cov = isTRUE(joint) && isTRUE(include_gp_uncertainty)
    )
    for (rr in seq_len(n_per_draw)) {
      row_id <- row_id + 1L
      f_std <- if (isTRUE(include_gp_uncertainty)) {
        if (isTRUE(joint)) {
          as.numeric(rmvnorm_psd(1L, pred$mean, pred$cov))
        } else {
          pred$mean + sqrt(pred$var) * rnorm(nrow(x_star))
        }
      } else {
        pred$mean
      }
      out[row_id, ] <- fit_obj$data$y_center + fit_obj$data$y_scale * f_std
    }
  }
  attr(out, "joint") <- isTRUE(joint)
  attr(out, "include_process_uncertainty") <- isTRUE(include_gp_uncertainty)
  attr(out, "latent_input_scale") <- u_input_scale
  out
}

gp_integrated_mean_state_1d <- function(x_train,
                                        u_train,
                                        y_train,
                                        x_star,
                                        u_mc,
                                        logtheta,
                                        sigma2_eps,
                                        kernel = "se",
                                        matern_nu = 2.5,
                                        return_cov = FALSE) {
  x_train <- as_study1_numeric_matrix(x_train, "x_train")
  x_star <- as_study1_numeric_matrix(
    x_star,
    "x_star",
    ncol_expected = ncol(x_train),
    expected_names = colnames(x_train)
  )
  u_mc <- as.matrix(u_mc)
  n <- nrow(x_train)
  N <- nrow(x_star)
  M <- nrow(u_mc)
  p_x <- ncol(x_train)
  if (length(u_train) != n || length(y_train) != n || ncol(u_mc) != N ||
      M < 1L || any(!is.finite(c(x_train, u_train, y_train, x_star, u_mc)))) {
    stop("Incompatible or nonfinite inputs in gp_integrated_mean_state_1d().")
  }

  kernel_spec <- normalize_gp_kernel_1d(kernel, matern_nu)
  theta_spec <- make_theta_spec_multix(p_x)
  rho <- exp(logtheta[1L])
  theta_x <- exp(logtheta[theta_spec$x_index])
  theta_u <- exp(logtheta[theta_spec$u_index])

  Dx_train <- lapply(seq_len(p_x), function(j) {
    pairwise_sqdist(x_train[, j, drop = FALSE])
  })
  D2_train <- weighted_distance_multix(Dx_train, theta_x) +
    theta_u * pairwise_sqdist(matrix(u_train, ncol = 1L))
  R_train <- kernel_from_weighted_sqdist_1d(
    D2_train, kernel_spec$name, kernel_spec$matern_nu
  )
  C_train <- sigma2_eps * (rho^2 * R_train + diag(n))
  chol_train <- safe_chol_1d(C_train)
  alpha <- solve_chol(chol_train, y_train)

  K_bar <- matrix(NA_real_, N, n)
  prior_diag <- numeric(N)
  for (i in seq_len(N)) {
    x_i <- x_star[rep(i, M), , drop = FALSE]
    u_i <- u_mc[, i]
    Dxi <- lapply(seq_len(p_x), function(j) {
      pairwise_sqdist(x_i[, j, drop = FALSE], x_train[, j, drop = FALSE])
    })
    D2_cross <- weighted_distance_multix(Dxi, theta_x) +
      theta_u * pairwise_sqdist(
        matrix(u_i, ncol = 1L), matrix(u_train, ncol = 1L)
      )
    K_cross <- rho^2 * sigma2_eps * kernel_from_weighted_sqdist_1d(
      D2_cross, kernel_spec$name, kernel_spec$matern_nu
    )
    K_bar[i, ] <- colMeans(K_cross)

    D2_ii <- theta_u * pairwise_sqdist(matrix(u_i, ncol = 1L))
    prior_diag[i] <- mean(
      rho^2 * sigma2_eps * kernel_from_weighted_sqdist_1d(
        D2_ii, kernel_spec$name, kernel_spec$matern_nu
      )
    )
  }

  conditional_mean <- as.numeric(K_bar %*% alpha)
  v <- forwardsolve(t(chol_train), t(K_bar))
  var_tolerance <- 100 * (n + N + M) * .Machine$double.eps *
    max(1, prior_diag)
  if (!isTRUE(return_cov)) {
    conditional_var <- prior_diag - colSums(v^2)
    if (any(conditional_var < -var_tolerance)) {
      stop("Integrated GP conditional variance is materially negative.")
    }
    conditional_var <- pmax(conditional_var, 0)
    return(list(mean = conditional_mean, var = conditional_var))
  }

  prior_cov <- matrix(NA_real_, N, N)
  diag(prior_cov) <- prior_diag
  if (N > 1L) {
    for (i in seq_len(N - 1L)) {
      u_i <- u_mc[, i]
      for (j in (i + 1L):N) {
        u_j <- u_mc[, j]
        D2_ij <- matrix(
          sum(theta_x * (x_star[i, ] - x_star[j, ])^2),
          M,
          M
        ) + theta_u * pairwise_sqdist(
          matrix(u_i, ncol = 1L), matrix(u_j, ncol = 1L)
        )
        value <- mean(
          rho^2 * sigma2_eps * kernel_from_weighted_sqdist_1d(
            D2_ij, kernel_spec$name, kernel_spec$matern_nu
          )
        )
        prior_cov[i, j] <- value
        prior_cov[j, i] <- value
      }
    }
  }
  conditional_cov <- prior_cov - crossprod(v)
  conditional_cov <- 0.5 * (conditional_cov + t(conditional_cov))
  if (any(diag(conditional_cov) < -var_tolerance)) {
    stop("Integrated GP conditional covariance has a materially negative diagonal.")
  }
  diag(conditional_cov) <- pmax(diag(conditional_cov), 0)
  list(
    mean = conditional_mean,
    var = diag(conditional_cov),
    cov = conditional_cov
  )
}

sample_eiv_m_given_xc_1d <- function(x_star_raw,
                                      c_star,
                                      fit_obj,
                                      draw_ids = NULL,
                                      n_per_draw = 1L,
                                      n_latent = 256L,
                                      include_process_uncertainty = TRUE,
                                      joint = FALSE,
                                      return_components = FALSE,
                                      seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  x_train <- fit_obj$data$x
  predictor_names <- fit_obj$data$predictor_names
  if (is.null(predictor_names)) predictor_names <- colnames(x_train)
  x_star_raw <- as_study1_numeric_matrix(
    x_star_raw,
    "x_star_raw",
    ncol_expected = ncol(x_train),
    expected_names = predictor_names
  )
  N <- nrow(x_star_raw)
  if (!(is.numeric(c_star) || is.integer(c_star)) || length(c_star) != N ||
      anyNA(c_star) || any(c_star != as.integer(c_star)) ||
      any(c_star < 1L | c_star > fit_obj$data$m)) {
    stop("c_star must contain one ordinal code 1, ..., m per evaluation row.")
  }
  c_star <- as.integer(c_star)
  x_star <- sweep(
    sweep(x_star_raw, 2L, fit_obj$data$x_center, "-"),
    2L, fit_obj$data$x_scale, "/"
  )

  n_saved <- nrow(fit_obj$mcmc$samples_u)
  if (is.null(draw_ids)) draw_ids <- seq_len(n_saved)
  draw_ids <- as.integer(draw_ids)
  if (length(draw_ids) < 1L || anyNA(draw_ids) ||
      any(draw_ids < 1L | draw_ids > n_saved)) {
    stop("draw_ids must index saved posterior draws.")
  }
  n_per_draw <- as.integer(n_per_draw)
  n_latent <- as.integer(n_latent)
  if (length(n_per_draw) != 1L || is.na(n_per_draw) || n_per_draw < 1L ||
      length(n_latent) != 1L || is.na(n_latent) || n_latent < 1L) {
    stop("n_per_draw and n_latent must be positive integers.")
  }

  kernel_spec <- kernel_spec_from_fit_1d(fit_obj)
  out <- matrix(NA_real_, length(draw_ids) * n_per_draw, N)
  conditional_means <- matrix(NA_real_, length(draw_ids), N)
  conditional_vars <- matrix(NA_real_, length(draw_ids), N)
  row_id <- 0L
  for (ii in seq_along(draw_ids)) {
    s <- draw_ids[ii]
    tau_s <- fit_obj$mcmc$samples_tau[s, ]
    lower <- c(-Inf, tau_s)[c_star]
    upper <- c(tau_s, Inf)[c_star]
    u_mc <- matrix(
      rtruncnorm_vec(
        mean = 0,
        sd = 1,
        lower = rep(lower, each = n_latent),
        upper = rep(upper, each = n_latent)
      ),
      nrow = n_latent,
      ncol = N
    )
    pred <- gp_integrated_mean_state_1d(
      x_train = x_train,
      u_train = fit_obj$mcmc$samples_u[s, ],
      y_train = fit_obj$data$y,
      x_star = x_star,
      u_mc = u_mc,
      logtheta = fit_obj$mcmc$samples_logtheta[s, ],
      sigma2_eps = fit_obj$mcmc$samples_sigma2[s],
      kernel = kernel_spec$name,
      matern_nu = kernel_spec$matern_nu,
      return_cov = isTRUE(joint) && isTRUE(include_process_uncertainty)
    )
    conditional_means[ii, ] <-
      fit_obj$data$y_center + fit_obj$data$y_scale * pred$mean
    conditional_vars[ii, ] <- fit_obj$data$y_scale^2 * pred$var
    for (rr in seq_len(n_per_draw)) {
      row_id <- row_id + 1L
      value_std <- if (!isTRUE(include_process_uncertainty)) {
        pred$mean
      } else if (isTRUE(joint)) {
        as.numeric(rmvnorm_psd(1L, pred$mean, pred$cov))
      } else {
        pred$mean + sqrt(pred$var) * rnorm(N)
      }
      out[row_id, ] <- fit_obj$data$y_center + fit_obj$data$y_scale * value_std
    }
  }
  attr(out, "target") <- "m(x,c)"
  attr(out, "n_latent") <- n_latent
  attr(out, "joint") <- isTRUE(joint)
  attr(out, "include_process_uncertainty") <- isTRUE(include_process_uncertainty)
  if (isTRUE(return_components)) {
    attr(out, "conditional_means") <- conditional_means
    attr(out, "conditional_vars") <- conditional_vars
  }
  out
}

sample_eiv_u_given_c_1d <- function(c_new,
                                     fit_obj,
                                     draw_ids = NULL,
                                     n_per_draw = 1L,
                                     scale = c("auto", "raw", "model"),
                                     seed = NULL) {
  scale <- match.arg(scale)
  anchored <- isTRUE(fit_obj$data$latent_scale_anchored)
  if (scale == "raw" && !anchored) {
    stop("scale='raw' requires calibration data that anchor latent U.")
  }
  if (scale == "auto") scale <- if (anchored) "raw" else "model"
  if (!is.null(seed)) set.seed(as.integer(seed))
  c_new <- as.integer(c_new)
  if (length(c_new) < 1L || anyNA(c_new) ||
      any(c_new < 1L | c_new > fit_obj$data$m)) {
    stop("c_new contains invalid ordinal codes.")
  }
  n_saved <- nrow(fit_obj$mcmc$samples_tau)
  if (is.null(draw_ids)) draw_ids <- seq_len(n_saved)
  draw_ids <- as.integer(draw_ids)
  if (length(draw_ids) < 1L || anyNA(draw_ids) ||
      any(draw_ids < 1L | draw_ids > n_saved)) {
    stop("draw_ids must index saved posterior draws.")
  }
  n_per_draw <- as.integer(n_per_draw)
  if (length(n_per_draw) != 1L || is.na(n_per_draw) || n_per_draw < 1L) {
    stop("n_per_draw must be a positive integer.")
  }
  out <- array(
    NA_real_,
    dim = c(length(draw_ids) * n_per_draw, length(c_new), 1L),
    dimnames = list(NULL, NULL, fit_obj$data$U_names)
  )
  row_id <- 0L
  for (s in draw_ids) {
    tau_s <- fit_obj$mcmc$samples_tau[s, ]
    lower <- c(-Inf, tau_s)[c_new]
    upper <- c(tau_s, Inf)[c_new]
    for (rr in seq_len(n_per_draw)) {
      row_id <- row_id + 1L
      out[row_id, , 1L] <- rtruncnorm_vec(0, 1, lower, upper)
    }
  }
  if (scale == "raw") {
    center <- if (is.null(fit_obj$data$U_center)) 0 else fit_obj$data$U_center
    scale_value <- if (is.null(fit_obj$data$U_scale)) 1 else fit_obj$data$U_scale
    out[, , 1L] <- center + scale_value * out[, , 1L]
  }
  attr(out, "source") <- "prospective"
  attr(out, "latent_sampler") <- "exact_truncated_normal"
  attr(out, "scale") <- if (!anchored && scale == "model") "working" else scale
  attr(out, "latent_scale_anchored") <- anchored
  out
}
