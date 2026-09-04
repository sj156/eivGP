############################################################
## Shared numerical, diagnostic, and scoring utilities
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
  pmax(outer(aa, bb, "+") - 2 * tcrossprod(a, b), 0)
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

aligned_chain_matrix <- function(chain_list) {
  chain_list <- lapply(chain_list, as.numeric)
  lens <- vapply(chain_list, length, integer(1))
  n0 <- min(lens)
  if (length(chain_list) < 2L || n0 < 4L) return(NULL)
  do.call(cbind, lapply(chain_list, function(x) tail(x, n0)))
}

tail_ess <- function(chain_list) {
  mat <- aligned_chain_matrix(chain_list)
  if (is.null(mat)) return(NA_real_)
  if (requireNamespace("posterior", quietly = TRUE)) {
    return(as.numeric(posterior::ess_tail(mat)))
  }
  NA_real_
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
    stop(
      "Predictive normal-mixture components must have finite, ",
      "nonnegative variances."
    )
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
