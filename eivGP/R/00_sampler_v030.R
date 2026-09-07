## eivGP 0.3.0: exact dense-GP, variance-collapsed, score-marginal sampler.
## This module is data-generic. Study designs and run budgets belong to callers.

mixedgp_v030_logadd <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log1p(exp(-abs(a - b)))
  out[is.infinite(m) & m < 0] <- -Inf
  out
}

mixedgp_v030_priors <- function(priors = list(), p, d, m_vec) {
  defaults <- list(
    variance_shape = 3, variance_rate = 2, signal_shape = c(32, 8),
    log_theta_x_mean = rep(log(0.5), p), log_theta_x_sd = rep(1.5, p),
    u_dictionary = rep(list(exp(log(0.5) + 1.5 *
      stats::qnorm((seq_len(5L) - 0.5) / 5))), d),
    u_weights = NULL, category_alpha = lapply(m_vec, function(m) rep(2, m)),
    loading_sd = 2.5
  )
  if (!is.list(priors) || (length(priors) && (is.null(names(priors)) ||
      any(!nzchar(names(priors))) || anyDuplicated(names(priors))))) {
    stop("priors must be a uniquely named list.")
  }
  bad <- setdiff(names(priors), names(defaults))
  if (length(bad)) stop("Unknown prior fields: ", paste(bad, collapse = ", "))
  out <- defaults
  out[names(priors)] <- priors
  positive <- function(x, n, name) {
    if (!is.numeric(x) || length(x) != n || any(!is.finite(x)) || any(x <= 0))
      stop(name, " must contain ", n, " positive finite value(s).")
    as.numeric(x)
  }
  out$variance_shape <- positive(out$variance_shape, 1L, "variance_shape")
  out$variance_rate <- positive(out$variance_rate, 1L, "variance_rate")
  out$signal_shape <- positive(out$signal_shape, 2L, "signal_shape")
  out$loading_sd <- positive(out$loading_sd, 1L, "loading_sd")
  if (length(out$log_theta_x_mean) == 1L) out$log_theta_x_mean <- rep(out$log_theta_x_mean, p)
  if (length(out$log_theta_x_sd) == 1L) out$log_theta_x_sd <- rep(out$log_theta_x_sd, p)
  if (!is.numeric(out$log_theta_x_mean) || length(out$log_theta_x_mean) != p ||
      any(!is.finite(out$log_theta_x_mean))) stop("log_theta_x_mean must be finite and have length p or 1.")
  out$log_theta_x_sd <- positive(out$log_theta_x_sd, p, "log_theta_x_sd")
  grid <- out$u_dictionary
  if (is.numeric(grid) && is.null(dim(grid))) grid <- rep(list(grid), d)
  if (is.matrix(grid)) {
    if (!is.numeric(grid) || ncol(grid) != d || !nrow(grid) ||
        any(!is.finite(grid)) || any(grid <= 0)) stop("u_dictionary matrix must be positive, with d columns.")
    if (anyDuplicated(as.data.frame(grid))) stop("u_dictionary vectors must be unique.")
    out$dictionary_type <- "vector"
    w <- out$u_weights
    if (is.null(w)) w <- rep(1, nrow(grid))
    w <- positive(w, nrow(grid), "u_weights")
    out$u_weights <- w / sum(w)
  } else {
    if (!is.list(grid) || length(grid) != d) stop("u_dictionary must be a vector, a list of d grids, or an M by d matrix.")
    for (k in seq_len(d)) {
      grid[[k]] <- positive(grid[[k]], length(grid[[k]]), "u_dictionary")
      if (!length(grid[[k]]) || anyDuplicated(grid[[k]])) stop("Coordinate grids must be nonempty and unique.")
    }
    out$dictionary_type <- "product"
    w <- out$u_weights
    if (is.null(w)) w <- lapply(grid, function(g) rep(1, length(g)))
    if (d == 1L && is.numeric(w)) w <- list(w)
    if (!is.list(w) || length(w) != d) stop("u_weights must be a list with d entries.")
    out$u_weights <- lapply(seq_len(d), function(k) {
      a <- positive(w[[k]], length(grid[[k]]), "u_weights")
      a / sum(a)
    })
  }
  out$u_dictionary <- grid
  if (length(m_vec) == 1L && is.numeric(out$category_alpha))
    out$category_alpha <- list(out$category_alpha)
  if (!is.list(out$category_alpha) || length(out$category_alpha) != length(m_vec))
    stop("category_alpha must supply one positive vector per ordinal item.")
  out$category_alpha <- lapply(seq_along(m_vec), function(j)
    positive(out$category_alpha[[j]], m_vec[j], "category_alpha"))
  out
}

mixedgp_v030_control <- function(control = list(), n_mis) {
  out <- list(u_block_size = 8L, cross_every = 10L, ess_max_steps = 10000L,
              dictionary_mode = "conditional", max_dictionary_size = 125L,
              gp_block_schur = TRUE)
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      any(!nzchar(names(control))) || anyDuplicated(names(control)))))
    stop("sampler_control must be a uniquely named list.")
  bad <- setdiff(names(control), names(out))
  if (length(bad)) stop("Unknown sampler_control fields: ", paste(bad, collapse = ", "))
  out <- utils::modifyList(out, control)
  for (nm in c("u_block_size", "ess_max_steps", "max_dictionary_size"))
    out[[nm]] <- mixedgp_as_integer_strict(out[[nm]], nm, 1L, 1L)
  out$cross_every <- mixedgp_as_integer_strict(out$cross_every, "cross_every", 0L, 1L)
  out$dictionary_mode <- match.arg(out$dictionary_mode, c("conditional", "marginal"))
  out$gp_block_schur <- mixedgp_validate_flag(out$gp_block_schur, "gp_block_schur")
  out
}

mixedgp_v030_beta_forward <- function(e, shape) {
  lower <- e <= 0
  prob <- stats::pnorm(abs(e), lower.tail = FALSE, log.p = TRUE)
  ans <- numeric(length(e))
  ans[lower] <- stats::qbeta(prob[lower], shape[1L], shape[2L], log.p = TRUE)
  ans[!lower] <- stats::qbeta(prob[!lower], shape[1L], shape[2L],
                             lower.tail = FALSE, log.p = TRUE)
  if (any(!is.finite(ans)) || any(ans <= 0 | ans >= 1))
    stop("Beta transformation exceeded floating-point resolution; no boundary clipping is applied.")
  ans
}

mixedgp_v030_stick_forward <- function(e, alpha) {
  if (length(e) != length(alpha) - 1L || any(!is.finite(e)) ||
      any(!is.finite(alpha) | alpha <= 0)) stop("Invalid Gaussian stick coordinates.")
  v <- vapply(seq_along(e), function(k)
    mixedgp_v030_beta_forward(e[k], c(alpha[k], sum(alpha[-seq_len(k)]))), numeric(1))
  log_rem <- c(0, cumsum(log1p(-v)))
  pi <- exp(c(log(v) + log_rem[seq_along(v)], tail(log_rem, 1L)))
  if (any(!is.finite(pi)) || any(pi <= 0)) stop("Category probabilities underflowed.")
  pi / sum(pi)
}

mixedgp_v030_stick_inverse <- function(pi, alpha) {
  if (length(pi) != length(alpha) || any(!is.finite(pi) | pi <= 0) ||
      abs(sum(pi) - 1) > 1e-10 || any(!is.finite(alpha) | alpha <= 0))
    stop("Invalid category probabilities.")
  vapply(seq_len(length(pi) - 1L), function(k) {
    v <- pi[k] / sum(pi[k:length(pi)])
    a <- alpha[k]; b <- sum(alpha[-seq_len(k)])
    lp <- stats::pbeta(v, a, b, log.p = TRUE)
    if (lp < log(0.5)) stats::qnorm(lp, log.p = TRUE) else
      stats::qnorm(stats::pbeta(v, a, b, lower.tail = FALSE, log.p = TRUE),
                   lower.tail = FALSE, log.p = TRUE)
  }, numeric(1))
}

mixedgp_v030_cutpoints <- function(pi, scale = 1) {
  b <- vapply(seq_len(length(pi) - 1L), function(k) {
    lo <- sum(pi[seq_len(k)]); hi <- sum(pi[-seq_len(k)])
    if (lo <= hi) stats::qnorm(log(lo), log.p = TRUE) else
      stats::qnorm(log(hi), lower.tail = FALSE, log.p = TRUE)
  }, numeric(1))
  c(-Inf, scale * b, Inf)
}

mixedgp_v030_interval_forward <- function(z, lower, upper) {
  n <- max(length(z), length(lower), length(upper))
  z <- rep_len(z, n); lower <- rep_len(lower, n); upper <- rep_len(upper, n)
  if (any(!is.finite(z)) || anyNA(lower) || anyNA(upper) || any(lower >= upper))
    stop("Invalid interval transform arguments.")
  mass <- log_normal_interval_prob(lower, upper)
  lp <- mixedgp_v030_logadd(stats::pnorm(lower, log.p = TRUE),
                           mass + stats::pnorm(z, log.p = TRUE))
  lq <- mixedgp_v030_logadd(stats::pnorm(upper, lower.tail = FALSE, log.p = TRUE),
                           mass + stats::pnorm(z, lower.tail = FALSE, log.p = TRUE))
  ans <- ifelse(lp <= lq, stats::qnorm(lp, log.p = TRUE),
                stats::qnorm(lq, lower.tail = FALSE, log.p = TRUE))
  if (any(!is.finite(ans)) || any(ans <= lower | ans >= upper))
    stop("Interval transformation exceeded floating-point resolution.")
  ans
}

mixedgp_v030_interval_inverse <- function(u, lower, upper) {
  n <- max(length(u), length(lower), length(upper))
  u <- rep_len(u, n); lower <- rep_len(lower, n); upper <- rep_len(upper, n)
  if (any(!is.finite(u)) || anyNA(lower) || anyNA(upper) ||
      any(u <= lower | u >= upper)) stop("Inputs must lie strictly inside their intervals.")
  mass <- log_normal_interval_prob(lower, upper)
  lp <- log_normal_interval_prob(lower, u) - mass
  lq <- log_normal_interval_prob(u, upper) - mass
  ans <- ifelse(lp <= lq, stats::qnorm(lp, log.p = TRUE),
                stats::qnorm(lq, lower.tail = FALSE, log.p = TRUE))
  if (any(!is.finite(ans))) stop("Inverse interval transformation exceeded floating-point resolution.")
  ans
}

mixedgp_v030_collapsed <- function(y, R, r, a = 3, b = 2) {
  n <- length(y)
  if (!is.matrix(R) || !identical(dim(R), c(n, n)) ||
      any(!is.finite(R)) || any(!is.finite(y)) ||
      length(r) != 1L || !is.finite(r) || r <= 0 || r >= 1 ||
      !is.finite(a) || !is.finite(b) || a <= 0 || b <= 0)
    stop("Invalid collapsed GP arguments.")
  B <- r * R
  diag(B) <- diag(B) + 1 - r
  ch <- tryCatch(chol(B), error = function(e) NULL)
  if (is.null(ch)) stop("Dense GP factorization failed; the covariance was not modified.")
  v <- forwardsolve(t(ch), y)
  Q <- sum(v * v)
  ld <- 2 * sum(log(diag(ch)))
  list(loglik = -0.5 * ld - (a + n / 2) * log(b + Q / 2),
       Q = Q, logdet = ld, chol = ch, B = B)
}

mixedgp_v030_ess <- function(v, loglik, max_steps = 10000L) {
  if (!length(v)) return(list(value = v, loglik = loglik(v), evaluations = 1L, angle = 0))
  if (any(!is.finite(v))) stop("ESS initial coordinates must be finite.")
  current <- loglik(v)
  if (length(current) != 1L || !is.finite(current))
    stop("ESS initial state has nonfinite log likelihood.")
  nu <- stats::rnorm(length(v))
  height <- current + log(stats::runif(1L))
  theta <- stats::runif(1L, 0, 2 * pi)
  lo <- theta - 2 * pi; hi <- theta
  for (step in seq_len(max_steps)) {
    proposal <- v * cos(theta) + nu * sin(theta)
    lp <- loglik(proposal)
    if (length(lp) != 1L || is.na(lp) || lp == Inf)
      stop("ESS likelihood evaluation is undefined.")
    if (lp >= height) return(list(value = proposal, loglik = lp,
                                  evaluations = step + 1L,
                                  angle = atan2(sin(theta), cos(theta))))
    if (theta < 0) lo <- theta else hi <- theta
    theta <- stats::runif(1L, lo, hi)
  }
  stop("ESS search budget exhausted; no successful transition was produced.")
}

mixedgp_v030_decode <- function(state, ctx) {
  prior <- ctx$priors
  pi <- lapply(seq_len(ctx$q), function(j)
    mixedgp_v030_stick_forward(state$e_pi[[j]], prior$category_alpha[[j]]))
  A <- matrix(0, ctx$q, ctx$d)
  if (ctx$measurement == "probit") {
    for (j in seq_len(ctx$q)) {
      active <- if (ctx$ident == "lower_triangular") seq_len(min(j, ctx$d)) else seq_len(ctx$d)
      A[j, active] <- prior$loading_sd * state$e_A[j, active]
      if (ctx$ident == "lower_triangular" && j <= ctx$d) {
        ## Half-normal quantile via log survival: P(|Z|>a)=2 P(Z>a).
        lq <- stats::pnorm(state$e_A[j, j], lower.tail = FALSE, log.p = TRUE) - log(2)
        A[j, j] <- prior$loading_sd * stats::qnorm(lq, lower.tail = FALSE, log.p = TRUE)
      }
    }
  }
  tau <- lapply(seq_len(ctx$q), function(j)
    mixedgp_v030_cutpoints(pi[[j]], if (ctx$measurement == "threshold") 1 else sqrt(1 + sum(A[j, ]^2))))
  theta_u <- if (prior$dictionary_type == "product")
    vapply(seq_len(ctx$d), function(k) prior$u_dictionary[[k]][state$J[k]], numeric(1)) else
    as.numeric(prior$u_dictionary[state$J, ])
  r <- mixedgp_v030_beta_forward(state$e_r, prior$signal_shape)
  list(pi = pi, tau = tau, A = A, r = r, theta_u = theta_u)
}

mixedgp_v030_corr <- function(state, dec, ctx) {
  weights <- c(exp(state$logtheta_x), dec$theta_u)
  if (any(!is.finite(weights))) stop("Kernel coefficients exceeded floating-point resolution.")
  Z <- cbind(ctx$X, state$U)
  D <- matrix(0, ctx$n, ctx$n)
  for (j in seq_along(weights)) D <- D + weights[j] * outer(Z[, j], Z[, j], "-")^2
  kernel_from_weighted_sqdist(D, ctx$kernel, ctx$matern_nu)
}

mixedgp_v030_measurement_loglik <- function(state, dec, ctx, rows = seq_len(ctx$n),
                                           items = seq_len(ctx$q)) {
  if (ctx$measurement == "threshold") {
    c <- ctx$C[rows, 1L]; u <- state$U[rows, 1L]; tau <- dec$tau[[1L]]
    return(if (all(u > tau[c] & u <= tau[c + 1L])) 0 else -Inf)
  }
  ans <- 0
  for (j in items) {
    mu <- as.numeric(state$U[rows, , drop = FALSE] %*% dec$A[j, ])
    c <- ctx$C[rows, j]
    ans <- ans + sum(log_normal_interval_prob(dec$tau[[j]][c] - mu,
                                              dec$tau[[j]][c + 1L] - mu))
  }
  ans
}

mixedgp_v030_dictionary <- function(priors, limit = Inf) {
  if (priors$dictionary_type == "vector") {
    if (nrow(priors$u_dictionary) > limit) stop("Dictionary exceeds max_dictionary_size for marginal updates.")
    return(list(indices = matrix(seq_len(nrow(priors$u_dictionary)), ncol = 1L),
                logweights = log(priors$u_weights)))
  }
  sizes <- lengths(priors$u_dictionary)
  if (sum(log(sizes)) > log(limit) + 1e-12)
    stop("Product dictionary exceeds max_dictionary_size; use conditional coordinate updates.")
  ind <- as.matrix(do.call(expand.grid, lapply(sizes, seq_len)))
  lw <- numeric(nrow(ind))
  for (k in seq_along(sizes)) lw <- lw + log(priors$u_weights[[k]][ind[, k]])
  list(indices = ind, logweights = lw)
}

mixedgp_v030_gp <- function(state, ctx, prepared = NULL) {
  dec <- mixedgp_v030_decode(state, ctx)
  key <- list(U = state$U, logtheta_x = state$logtheta_x, J = state$J, r = dec$r)
  if (!is.null(ctx$cache$key) && identical(ctx$cache$key, key)) return(ctx$cache$value)
  R <- mixedgp_v030_corr(state, dec, ctx)
  value <- NULL
  if (!is.null(prepared) && !identical(prepared$key, list(
      U = state$U[prepared$J, , drop = FALSE], logtheta_x = state$logtheta_x,
      J = state$J, r = dec$r))) prepared <- NULL
  if (!is.null(prepared)) {
    value <- tryCatch({
      B <- dec$r * R; diag(B) <- diag(B) + 1 - dec$r
      I <- prepared$I; J <- prepared$J
      cross <- B[I, J, drop = FALSE]
      z <- forwardsolve(t(prepared$chol), t(cross))
      S <- B[I, I, drop = FALSE] - crossprod(z)
      ch <- chol(S)
      e <- ctx$y[I] - as.numeric(cross %*% prepared$solve_y)
      v <- forwardsolve(t(ch), e)
      Q <- prepared$Q + sum(v * v)
      ld <- prepared$logdet + 2 * sum(log(diag(ch)))
      ll <- -0.5 * ld - (ctx$priors$variance_shape + ctx$n / 2) *
        log(ctx$priors$variance_rate + Q / 2)
      if (any(!is.finite(c(Q, ld, ll))) || Q < 0) stop("Invalid Schur calculation.")
      list(loglik = ll, Q = Q, logdet = ld, B = B, chol = NULL)
    }, error = function(e) NULL)
    if (is.null(value)) {
      ctx$counts$gp_block_fallbacks <- ctx$counts$gp_block_fallbacks + 1L
    } else ctx$counts$gp_block_evaluations <- ctx$counts$gp_block_evaluations + 1L
  }
  if (is.null(value)) {
    ctx$counts$gp_full_factorizations <- ctx$counts$gp_full_factorizations + 1L
    value <- mixedgp_v030_collapsed(ctx$y, R, dec$r,
                                    ctx$priors$variance_shape, ctx$priors$variance_rate)
  }
  ctx$cache$key <- key; ctx$cache$value <- value
  value
}

mixedgp_v030_prepare_block <- function(state, ctx, rows) {
  if (!ctx$control$gp_block_schur || ctx$control$dictionary_mode == "marginal" ||
      !length(rows) || length(rows) == ctx$n) return(NULL)
  base <- mixedgp_v030_gp(state, ctx)
  J <- setdiff(seq_len(ctx$n), rows)
  ctx$counts$gp_block_setups <- ctx$counts$gp_block_setups + 1L
  tryCatch({
    ch <- chol(base$B[J, J, drop = FALSE])
    v <- forwardsolve(t(ch), ctx$y[J])
    list(I = rows, J = J, chol = ch, solve_y = backsolve(ch, v),
         Q = sum(v * v), logdet = 2 * sum(log(diag(ch))),
         key = list(U = state$U[J, , drop = FALSE], logtheta_x = state$logtheta_x,
                    J = state$J, r = mixedgp_v030_decode(state, ctx)$r))
  }, error = function(e) NULL)
}

mixedgp_v030_loggp <- function(state, ctx, prepared = NULL) {
  if (ctx$control$dictionary_mode == "conditional")
    return(mixedgp_v030_gp(state, ctx, prepared)$loglik)
  dict <- ctx$dictionary
  lp <- vapply(seq_len(nrow(dict$indices)), function(h) {
    st <- state; st$J <- as.integer(dict$indices[h, ])
    mixedgp_v030_gp(st, ctx)$loglik + dict$logweights[h]
  }, numeric(1))
  max(lp) + log(sum(exp(lp - max(lp))))
}

mixedgp_v030_update_dictionary <- function(state, ctx) {
  prior <- ctx$priors
  draw_one <- function(indices, weights) {
    lp <- vapply(seq_len(nrow(indices)), function(h) {
      st <- state; st$J <- as.integer(indices[h, ])
      mixedgp_v030_gp(st, ctx)$loglik + log(weights[h])
    }, numeric(1))
    if (any(!is.finite(lp))) stop("Nonfinite dictionary likelihood.")
    ctx$counts$dictionary_evaluations <- ctx$counts$dictionary_evaluations + length(lp)
    probs <- exp(lp - max(lp))
    sample.int(length(lp), 1L, prob = probs)
  }
  before <- state$J
  if (prior$dictionary_type == "vector") {
    state$J <- as.integer(draw_one(matrix(seq_len(nrow(prior$u_dictionary)), ncol = 1L),
                                    prior$u_weights))
  } else {
    for (k in seq_len(ctx$d)) {
      sz <- length(prior$u_dictionary[[k]])
      choices <- matrix(rep(state$J, each = sz), nrow = sz)
      choices[, k] <- seq_len(sz)
      state$J[k] <- draw_one(choices, prior$u_weights[[k]])
    }
  }
  ctx$counts$dictionary_switches <- ctx$counts$dictionary_switches + sum(state$J != before)
  state
}

mixedgp_v030_initial_state <- function(ctx) {
  prior <- ctx$priors
  state <- list(U = matrix(stats::rnorm(ctx$n * ctx$d), ctx$n, ctx$d),
                e_pi = lapply(ctx$m_vec, function(m) stats::rnorm(m - 1L, sd = 0.5)),
                e_A = matrix(stats::rnorm(ctx$q * ctx$d, sd = 0.5), ctx$q, ctx$d),
                logtheta_x = stats::rnorm(ctx$p, prior$log_theta_x_mean,
                                          prior$log_theta_x_sd / 2),
                e_r = stats::rnorm(1L, sd = 0.5))
  state$J <- if (prior$dictionary_type == "product")
    vapply(seq_len(ctx$d), function(k) sample.int(length(prior$u_dictionary[[k]]),
                                                 1L, prob = prior$u_weights[[k]]), integer(1)) else
    sample.int(nrow(prior$u_dictionary), 1L, prob = prior$u_weights)
  if (length(ctx$observed)) state$U[ctx$observed, ] <- ctx$U_obs[ctx$observed, , drop = FALSE]
  if (ctx$measurement == "threshold") {
    ## Calibration-compatible increasing cutpoints, including unobserved levels.
    ## Work between consecutive distinct calibrated classes and interpolate all
    ## intervening boundaries; this avoids arbitrary raw-cutpoint bounds.
    obs <- ctx$observed
    if (length(obs)) {
      cls <- sort(unique(ctx$C[obs, 1L]))
      mins <- vapply(cls, function(c) min(ctx$U_obs[obs[ctx$C[obs, 1L] == c], 1L]), numeric(1))
      maxs <- vapply(cls, function(c) max(ctx$U_obs[obs[ctx$C[obs, 1L] == c], 1L]), numeric(1))
      if (length(cls) > 1L && any(head(maxs, -1L) >= tail(mins, -1L)))
        stop("Calibration values are incompatible with deterministic ordinal thresholds.")
      tau <- numeric(ctx$m_vec[1L] - 1L)
      left <- seq_len(cls[1L] - 1L)
      if (length(left)) tau[left] <- mins[1L] - rev(seq_along(left))
      if (length(cls) > 1L) for (j in seq_len(length(cls) - 1L)) {
        ids <- seq.int(cls[j], cls[j + 1L] - 1L)
        tau[ids] <- maxs[j] + (mins[j + 1L] - maxs[j]) *
          seq_along(ids) / (length(ids) + 1)
      }
      right <- if (tail(cls, 1L) < ctx$m_vec[1L])
        seq.int(tail(cls, 1L), ctx$m_vec[1L] - 1L) else integer(0)
      if (length(right)) tau[right] <- tail(maxs, 1L) + seq_along(right)
      edges <- c(-Inf, tau, Inf)
      pi <- exp(log_normal_interval_prob(head(edges, -1L), tail(edges, -1L)))
      if (any(pi <= 0) || any(!is.finite(pi))) stop("Calibrated cutpoints exceed numerical prior support.")
      pi <- pi / sum(pi)
      state$e_pi[[1L]] <- mixedgp_v030_stick_inverse(pi, prior$category_alpha[[1L]])
    }
    dec <- mixedgp_v030_decode(state, ctx)
    mis <- ctx$missing
    if (length(mis)) {
      c <- ctx$C[mis, 1L]
      state$U[mis, 1L] <- mixedgp_v030_interval_forward(stats::rnorm(length(mis)),
                                                       dec$tau[[1L]][c], dec$tau[[1L]][c + 1L])
    }
  }
  state
}

mixedgp_v030_record_ess <- function(ctx, result, block) {
  nm <- paste0(block, "_ess_total")
  en <- paste0(block, "_ess_eval_total")
  ctx$counts[[nm]] <- ctx$counts[[nm]] + 1L
  ctx$counts[[en]] <- ctx$counts[[en]] + result$evaluations
  ctx$counts$abs_angle_sum <- ctx$counts$abs_angle_sum + abs(result$angle)
  ctx$counts$angle_count <- ctx$counts$angle_count + 1L
}

mixedgp_v030_sweep <- function(state, ctx, iteration) {
  mis <- ctx$missing
  blocks <- if (length(mis)) split(mis[sample.int(length(mis))],
                ceiling(seq_along(mis) / ctx$control$u_block_size)) else list()
  for (rows in blocks) {
    dec <- mixedgp_v030_decode(state, ctx)
    prep <- mixedgp_v030_prepare_block(state, ctx, rows)
    if (ctx$measurement == "threshold") {
      c <- ctx$C[rows, 1L]
      lo <- dec$tau[[1L]][c]; hi <- dec$tau[[1L]][c + 1L]
      z <- mixedgp_v030_interval_inverse(state$U[rows, 1L], lo, hi)
      propose <- function(v) {
        st <- state
        st$U[rows, 1L] <- mixedgp_v030_interval_forward(v, lo, hi)
        st
      }
      residual <- function(v) mixedgp_v030_loggp(propose(v), ctx, prep)
    } else {
      z <- as.numeric(state$U[rows, , drop = FALSE])
      propose <- function(v) {
        st <- state; st$U[rows, ] <- matrix(v, length(rows), ctx$d); st
      }
      residual <- function(v) {
        st <- propose(v)
        mixedgp_v030_loggp(st, ctx, prep) +
          mixedgp_v030_measurement_loglik(st, dec, ctx, rows = rows)
      }
    }
    update <- mixedgp_v030_ess(z, residual, ctx$control$ess_max_steps)
    state <- propose(update$value)
    mixedgp_v030_record_ess(ctx, update, "u")
  }

  for (j in seq_len(ctx$q)) {
    active <- if (ctx$ident == "lower_triangular") seq_len(min(j, ctx$d)) else seq_len(ctx$d)
    z <- state$e_pi[[j]]
    if (ctx$measurement == "probit") z <- c(state$e_A[j, active], z)
    propose <- function(v) {
      st <- state
      if (ctx$measurement == "probit") {
        st$e_A[j, active] <- v[seq_along(active)]
        v <- v[-seq_along(active)]
      }
      st$e_pi[[j]] <- v
      st
    }
    residual <- function(v) {
      st <- propose(v)
      mixedgp_v030_measurement_loglik(st, mixedgp_v030_decode(st, ctx), ctx, items = j)
    }
    update <- mixedgp_v030_ess(z, residual, ctx$control$ess_max_steps)
    state <- propose(update$value)
    mixedgp_v030_record_ess(ctx, update, "measurement")
  }

  if (ctx$control$cross_every > 0L && iteration %% ctx$control$cross_every == 0L &&
      length(mis)) {
    if (ctx$measurement == "threshold") {
      dec <- mixedgp_v030_decode(state, ctx)
      c <- ctx$C[mis, 1L]
      zmis <- mixedgp_v030_interval_inverse(state$U[mis, 1L],
                                             dec$tau[[1L]][c], dec$tau[[1L]][c + 1L])
      propose <- function(v) {
        st <- state; st$e_pi[[1L]] <- v
        dd <- mixedgp_v030_decode(st, ctx)
        st$U[mis, 1L] <- mixedgp_v030_interval_forward(zmis, dd$tau[[1L]][c],
                                                      dd$tau[[1L]][c + 1L])
        st
      }
      residual <- function(v) {
        st <- propose(v); dd <- mixedgp_v030_decode(st, ctx)
        ll <- mixedgp_v030_measurement_loglik(st, dd, ctx, rows = ctx$observed)
        if (!is.finite(ll)) return(-Inf)
        mixedgp_v030_loggp(st, ctx) + sum(log(dd$pi[[1L]][c]))
      }
      z <- state$e_pi[[1L]]
    } else {
      ## Cycle measurement rows deterministically; the selected subject block is
      ## random independently of its values. Every coordinate of a subject moves.
      j <- 1L + ((iteration %/% ctx$control$cross_every - 1L) %% ctx$q)
      rows <- blocks[[1L]]
      active <- if (ctx$ident == "lower_triangular") seq_len(min(j, ctx$d)) else seq_len(ctx$d)
      na <- length(active); np <- length(state$e_pi[[j]])
      z <- c(state$e_A[j, active], state$e_pi[[j]], as.numeric(state$U[rows, , drop = FALSE]))
      propose <- function(v) {
        st <- state
        st$e_A[j, active] <- v[seq_len(na)]
        st$e_pi[[j]] <- v[na + seq_len(np)]
        st$U[rows, ] <- matrix(v[-seq_len(na + np)], length(rows), ctx$d)
        st
      }
      residual <- function(v) {
        st <- propose(v); dd <- mixedgp_v030_decode(st, ctx)
        ## Full ordinal product counts the affected-factor union once.
        mixedgp_v030_loggp(st, ctx) + mixedgp_v030_measurement_loglik(st, dd, ctx)
      }
    }
    update <- mixedgp_v030_ess(z, residual, ctx$control$ess_max_steps)
    state <- propose(update$value)
    mixedgp_v030_record_ess(ctx, update, "cross")
  }

  prior <- ctx$priors
  z <- c((state$logtheta_x - prior$log_theta_x_mean) / prior$log_theta_x_sd, state$e_r)
  propose <- function(v) {
    st <- state
    st$logtheta_x <- prior$log_theta_x_mean + prior$log_theta_x_sd * v[seq_len(ctx$p)]
    st$e_r <- v[ctx$p + 1L]
    st
  }
  update <- mixedgp_v030_ess(z, function(v) mixedgp_v030_loggp(propose(v), ctx),
                             ctx$control$ess_max_steps)
  state <- propose(update$value)
  mixedgp_v030_record_ess(ctx, update, "theta")
  ## In marginal mode ALL preceding blocks are J-marginal; restore J before
  ## conditional recovery/retention. Product conditional mode uses coordinate Gibbs.
  if (ctx$control$dictionary_mode == "marginal") {
    dict <- ctx$dictionary
    lp <- vapply(seq_len(nrow(dict$indices)), function(h) {
      st <- state; st$J <- as.integer(dict$indices[h, ])
      mixedgp_v030_gp(st, ctx)$loglik + dict$logweights[h]
    }, numeric(1))
    h <- sample.int(length(lp), 1L, prob = exp(lp - max(lp)))
    old <- state$J
    state$J <- as.integer(dict$indices[h, ])
    ctx$counts$dictionary_evaluations <- ctx$counts$dictionary_evaluations + length(lp)
    ctx$counts$dictionary_switches <- ctx$counts$dictionary_switches + sum(state$J != old)
  } else state <- mixedgp_v030_update_dictionary(state, ctx)
  state
}

mixedgp_v030_fit <- function(X, y, C, U_obs, m_vec,
                             measurement = c("threshold", "probit"), ident = "none",
                             kernel = "se", matern_nu = 2.5,
                             n_iter = 1750L, burn = 500L, n_chains = 4L, seed = 1L,
                             parallel_chains = TRUE, n_cores = NULL,
                             priors = list(), sampler_control = list(),
                             store_scores = FALSE, .resume = NULL) {
  measurement <- match.arg(measurement)
  ident <- match.arg(ident, c("none", "lower_triangular"))
  spec <- normalize_gp_kernel(kernel, matern_nu)
  X <- as_numeric_matrix_strict(X, "X")
  U_obs <- as_numeric_matrix_strict(U_obs, "U_obs", nrow_expected = nrow(X), allow_na = TRUE)
  n <- nrow(X); d <- ncol(U_obs); p <- ncol(X)
  C <- as.matrix(C)
  if (!is.numeric(y) || length(y) != n || any(!is.finite(y)) ||
      !is.numeric(C) || nrow(C) != n || anyNA(C) ||
      length(m_vec) != ncol(C) || any(C != floor(C))) stop("Invalid sampler data.")
  for (j in seq_along(m_vec)) if (m_vec[j] < 2 || any(C[, j] < 1 | C[, j] > m_vec[j]))
    stop("Ordinal codes must lie within the declared categories.")
  count <- rowSums(is.finite(U_obs))
  if (any(count != 0 & count != d)) stop("Partially observed U rows are not supported.")
  observed <- which(count == d); missing_rows <- setdiff(seq_len(n), observed)
  if (measurement == "threshold" && (d != 1L || ncol(C) != 1L))
    stop("Deterministic threshold measurement requires one latent input and one ordinal item.")
  if (ident == "lower_triangular" && ncol(C) < d) stop("Triangular loadings require q >= d.")
  for (nm in c("n_iter", "n_chains", "seed", "burn")) {
    value <- mixedgp_as_integer_strict(get(nm), nm, if (nm %in% c("burn", "seed")) 0L else 1L, 1L)
    assign(nm, value)
  }
  if (n_iter <= burn) stop("n_iter must exceed burn.")
  priors <- mixedgp_v030_priors(priors, p, d, m_vec)
  control <- mixedgp_v030_control(sampler_control, length(missing_rows))
  dictionary <- if (control$dictionary_mode == "marginal")
    mixedgp_v030_dictionary(priors, control$max_dictionary_size) else NULL
  start <- if (is.null(.resume)) 0L else .resume$iteration
  if (!is.null(.resume) && (!identical(.resume$version, 2L) ||
      !identical(.resume$sampler_version, "0.3.0") || length(.resume$states) != n_chains ||
      .resume$burn != burn || start >= n_iter)) stop("Incompatible 0.3.0 continuation state.")
  seeds <- as.integer((as.double(seed) + 10000 * seq_len(n_chains)) %% .Machine$integer.max)
  backend <- mixedgp_parallel_backend(if (parallel_chains) n_cores else 1L)
  workers <- min(n_chains, backend$cores)
  use_fork <- parallel_chains && backend$backend == "fork" && workers > 1L
  if (!use_fork) workers <- 1L
  run_chain <- function(cc) {
    ctx <- new.env(parent = environment())
    ctx$X <- X; ctx$y <- y; ctx$C <- C; ctx$U_obs <- U_obs
    ctx$n <- n; ctx$p <- p; ctx$q <- ncol(C); ctx$d <- d; ctx$m_vec <- m_vec
    ctx$observed <- observed; ctx$missing <- missing_rows
    ctx$measurement <- measurement; ctx$ident <- ident
    ctx$kernel <- spec$name; ctx$matern_nu <- spec$matern_nu
    ctx$priors <- priors; ctx$control <- control; ctx$dictionary <- dictionary
    ctx$cache <- new.env(parent = emptyenv())
    ctx$counts <- new.env(parent = emptyenv())
    names <- c("u_ess_total", "u_ess_eval_total",
      "measurement_ess_total", "measurement_ess_eval_total", "cross_ess_total",
      "cross_ess_eval_total", "theta_ess_total", "theta_ess_eval_total",
      "dictionary_evaluations", "dictionary_switches", "gp_full_factorizations",
      "gp_block_setups", "gp_block_evaluations", "gp_block_fallbacks",
      "abs_angle_sum", "angle_count")
    for (name in names) ctx$counts[[name]] <- 0
    if (is.null(.resume)) {
      set.seed(seeds[cc])
      state <- mixedgp_v030_initial_state(ctx)
    } else {
      cp <- .resume$states[[cc]]
      mixedgp_restore_rng_state(cp$rng)
      state <- cp$state
      if (!is.null(cp$cache)) {
        ctx$cache$key <- cp$cache$key
        ctx$cache$value <- cp$cache$value
      }
    }
    initial <- state
    dec <- mixedgp_v030_decode(state, ctx)
    if (!is.finite(mixedgp_v030_measurement_loglik(state, dec, ctx)))
      stop("Initial measurement likelihood is nonfinite.")
    if (is.null(.resume)) mixedgp_v030_gp(state, ctx)
    keep <- n_iter - max(burn, start)
    draws <- list(
      samples_U = array(NA_real_, c(keep, n, d)),
      samples_A = if (measurement == "probit") array(NA_real_, c(keep, ncol(C), d)) else NULL,
      samples_S = if (store_scores && measurement == "probit") array(NA_real_, c(keep, n, ncol(C))) else NULL,
      samples_tau = matrix(NA_real_, keep, sum(m_vec - 1L)),
      samples_logtheta = matrix(NA_real_, keep, 1L + p + d),
      samples_sigma2 = numeric(keep), samples_V = numeric(keep), samples_r = numeric(keep),
      samples_J = matrix(NA_integer_, keep, if (priors$dictionary_type == "product") d else 1L),
      samples_pi = matrix(NA_real_, keep, sum(m_vec))
    )
    colnames(draws$samples_tau) <- if (measurement == "threshold") paste0("tau", seq_len(m_vec - 1L)) else
      tau_names_from_mvec(m_vec)
    colnames(draws$samples_logtheta) <- c("log_rho", paste0("log_theta_x", seq_len(p)),
                                         paste0("log_theta_u", seq_len(d)))
    saved <- 0L
    for (iteration in seq.int(start + 1L, n_iter)) {
      state <- mixedgp_v030_sweep(state, ctx, iteration)
      if (iteration <= burn) next
      saved <- saved + 1L
      dec <- mixedgp_v030_decode(state, ctx)
      if (length(observed) && !isTRUE(all.equal(unname(state$U[observed, , drop = FALSE]),
            unname(U_obs[observed, , drop = FALSE]), tolerance = 0)))
        stop("Calibrated inputs changed.")
      if (!is.finite(mixedgp_v030_measurement_loglik(state, dec, ctx)))
        stop("Retained state has invalid ordinal likelihood.")
      gp <- mixedgp_v030_gp(state, ctx)
      V <- 1 / stats::rgamma(1L, shape = priors$variance_shape + n / 2,
                             rate = priors$variance_rate + gp$Q / 2)
      if (!is.finite(V) || V <= 0) stop("Conditional variance recovery failed.")
      draws$samples_U[saved, , ] <- state$U
      if (measurement == "probit") draws$samples_A[saved, , ] <- dec$A
      if (store_scores && measurement == "probit") {
        S <- matrix(NA_real_, n, ncol(C))
        for (j in seq_len(ncol(C))) {
          mu <- as.numeric(state$U %*% dec$A[j, ])
          S[, j] <- mu + mixedgp_v030_interval_forward(stats::rnorm(n),
            dec$tau[[j]][C[, j]] - mu, dec$tau[[j]][C[, j] + 1L] - mu)
        }
        draws$samples_S[saved, , ] <- S
      }
      draws$samples_tau[saved, ] <- unlist(lapply(dec$tau, function(t) t[-c(1L, length(t))]), use.names = FALSE)
      draws$samples_logtheta[saved, ] <- c(0.5 * (log(dec$r) - log1p(-dec$r)),
                                           state$logtheta_x, log(dec$theta_u))
      draws$samples_sigma2[saved] <- V * (1 - dec$r)
      draws$samples_V[saved] <- V; draws$samples_r[saved] <- dec$r
      draws$samples_J[saved, ] <- state$J
      draws$samples_pi[saved, ] <- unlist(dec$pi, use.names = FALSE)
    }
    counters <- as.list(ctx$counts)
    counters$theta_eval_total <- counters$theta_ess_eval_total
    counters$mean_abs_angle <- counters$abs_angle_sum / max(1, counters$angle_count)
    stats <- as.data.frame(c(list(chain = cc, seed = seeds[cc], saved = saved), counters))
    c(draws, list(chain_id = cc, seed = seeds[cc], initial_state = initial,
      checkpoint = list(state = state, rng = mixedgp_rng_state(),
                        cache = list(key = ctx$cache$key, value = ctx$cache$value)), stats = stats))
  }
  timing <- system.time({
    chains <- mixedgp_parallel_lapply(as.list(seq_len(n_chains)), run_chain,
      n_cores = workers, seeds = seeds, mc.preschedule = FALSE)
  })
  list(chains = chains, priors = priors, control = control,
       iteration = n_iter, burn = burn, sampler_version = "0.3.0",
       time_seconds = unname(timing["elapsed"]),
       parallel_backend = if (use_fork) "fork" else "serial", parallel_cores = workers)
}
