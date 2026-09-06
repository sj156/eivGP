## Run from repository root:
## Rscript revision/codes/tests/test_study2_joint_measurement.R
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(script_arg) == 1L)
codes_dir <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
source(file.path(codes_dir, "00_study2_functions.R"))
expect_error <- function(expr, pattern) {
  result <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(result, "error"), grepl(pattern, conditionMessage(result)))
}
row_logprior <- function(row, j, d, ident, s_A, B) {
  active <- measurement_row_active(j, d, ident)
  sum(dnorm(row$A[active], 0, s_A, log = TRUE)) +
    ifelse(ident == "lower_triangular" && j <= d, log(2), 0) +
    lgamma(length(row$tau) + 1L) - length(row$tau) * log(2 * B)
}

## Independent change-of-variables identity, numerical Jacobian, and support.
set.seed(260906)
identity_error <- inverse_error <- jacobian_error <- 0
for (ident in c("none", "lower_triangular")) {
  for (j in 1:3) for (K in 0:3) {
    d <- 2L; s_A <- 1.3; B <- 3
    k_A <- length(measurement_row_active(j, d, ident))
    z <- rnorm(k_A + K)
    row <- measurement_row_from_gaussian(z, j, d, K, ident, s_A, B)
    inverse <- measurement_row_to_gaussian(row$A, row$tau, j, ident, s_A, B)
    inverse_error <- max(inverse_error, abs(inverse - z))
    identity_error <- max(identity_error, abs(row_logprior(row, j, d, ident, s_A, B) +
      row$log_jacobian - sum(dnorm(z, log = TRUE))))
    physical <- function(value) {
      out <- measurement_row_from_gaussian(value, j, d, K, ident, s_A, B)
      c(out$A[measurement_row_active(j, d, ident)], out$tau)
    }
    J <- matrix(vapply(seq_along(z), function(k) {
      plus <- minus <- z; plus[k] <- plus[k] + 1e-5; minus[k] <- minus[k] - 1e-5
      (physical(plus) - physical(minus)) / 2e-5
    }, numeric(length(z))), nrow = length(z))
    jacobian_error <- max(jacobian_error,
      abs(as.numeric(determinant(J, logarithm = TRUE)$modulus) - row$log_jacobian))
  }
}
stopifnot(identity_error < 1e-10, inverse_error < 1e-10, jacobian_error < 1e-7)
cat("Prior/Jacobian, inverse, finite-difference errors:",
    identity_error, inverse_error, jacobian_error, "\n")
for (A_diag in c(1e-12, 1e-6, .2, 12)) {
  tau <- c(-3 + 1e-8, .25, .25 + 1e-8, 3 - 1e-8)
  z <- measurement_row_to_gaussian(A_diag, tau, 1L, "lower_triangular", 1.3, 3)
  row <- measurement_row_from_gaussian(z, 1L, 1L, 4L, "lower_triangular", 1.3, 3)
  stopifnot(all(is.finite(z)), abs(row$A / A_diag - 1) < 1e-8,
            max(abs(row$tau - tau)) < 1e-12)
}
stopifnot(is.null(measurement_row_from_gaussian(c(0, 1000), 1L, 1L, 1L,
                                               "lower_triangular")))
cat("Tiny positive diagonal and near-adjacent/boundary cut-point round trips: PASS\n")

## Independent two-dimensional quadrature conditional reference. This is a
## stationary-target check for the cheap measurement kernel, not a GP mixing claim.
U <- matrix(c(-1.4, -.9, -.5, -.2, .2, .5, .9, 1.4), ncol = 1L)
C <- matrix(c(1L, 1L, 2L, 1L, 2L, 1L, 2L, 2L), ncol = 1L)
s_A <- 1.1; B <- 2
grid <- expand.grid(A = seq(0, 6, length.out = 401),
                    tau = seq(-B, B, length.out = 401))
logw <- dnorm(grid$A, 0, s_A, log = TRUE)
for (i in seq_len(nrow(U))) {
  logw <- logw + pnorm(grid$tau - grid$A * U[i, 1L],
                       lower.tail = C[i, 1L] == 1L, log.p = TRUE)
}
w <- exp(logw - max(logw))
w[grid$A %in% range(grid$A)] <- w[grid$A %in% range(grid$A)] / 2
w[grid$tau %in% range(grid$tau)] <- w[grid$tau %in% range(grid$tau)] / 2
w <- w / sum(w)
reference <- colSums(as.matrix(grid) * w)
A <- matrix(.7, 1L, 1L); tau <- list(.2)
draws <- matrix(NA_real_, 10000L, 2L, dimnames = list(NULL, c("A", "tau")))
evaluations <- 0L
for (iteration in seq_len(nrow(draws))) {
  result <- update_measurement_joint_ess(C, U, A, tau, 2L, s_A, B,
                                         "lower_triangular")
  A <- result$A; tau <- result$tau
  draws[iteration, ] <- c(A, tau[[1L]])
  evaluations <- evaluations + result$n_eval
}
draws <- draws[-seq_len(2000L), , drop = FALSE]
batch_means <- rowsum(draws, rep(seq_len(80L), each = 100L)) / 100
mcse <- apply(batch_means, 2L, sd) / sqrt(nrow(batch_means))
observed <- colMeans(draws)
print(data.frame(target = names(reference), reference, observed, mcse), row.names = FALSE)
stopifnot(all(abs(observed - reference) < pmax(.045, 7 * mcse)))
cat("Joint row conditional quadrature reference: PASS; evaluations", evaluations, "\n")
expect_error(update_measurement_joint_ess(C, U, A, tau, 2L, max_try = 0L),
             "exceeded max_try")

## Independent transport change-of-variables identity for no, singleton,
## multiple, and all missing rows. The difference must not depend on A.
set.seed(119)
n <- 6L; d <- 2L; q <- 3L
S <- matrix(rnorm(n * q), n, q)
U0 <- matrix(rnorm(n * d), n, d)
A0 <- rbind(c(.9, 0), c(-.2, 1.3), c(.6, -.8))
transport_error <- 0
for (missing in list(integer(0), 6L, c(2L, 4L, 6L), seq_len(n))) {
  ref0 <- loading_transport_reference(S, A0)
  W <- (U0[missing, , drop = FALSE] - ref0$mean[missing, , drop = FALSE]) %*%
    t(ref0$precision_chol)
  differences <- numeric(20L)
  for (i in seq_along(differences)) {
    Ap <- rbind(c(exp(rnorm(1)), 0), c(rnorm(1), exp(rnorm(1))), rnorm(2))
    ref <- loading_transport_reference(S, Ap)
    Up <- U0
    if (length(missing)) Up[missing, ] <- ref$mean[missing, , drop = FALSE] +
      t(backsolve(ref$precision_chol, t(W)))
    old_density <- sum(dnorm(Up[missing, , drop = FALSE], log = TRUE)) +
      sum(dnorm(S - Up %*% t(Ap), log = TRUE))
    logJ <- -length(missing) * sum(log(diag(ref$precision_chol)))
    residual <- loading_transport_score_loglik(S, Up, Ap, missing)
    differences[i] <- old_density + logJ - residual
    stopifnot(all(Up[setdiff(seq_len(n), missing), , drop = FALSE] ==
                    U0[setdiff(seq_len(n), missing), , drop = FALSE]))
  }
  transport_error <- max(transport_error, abs(differences - sum(dnorm(W, log = TRUE))))
}
stopifnot(transport_error < 1e-10)
cat("Transport density/Jacobian cancellation error:", transport_error, "\n")

## A separate one-dimensional quadrature reference in transported coordinates.
X <- matrix(c(-.8, -.1, .4, 1), ncol = 1L)
y <- c(-.7, .2, .4, .9)
S <- matrix(c(-1.2, -.4, .5, 1.4), ncol = 1L)
U <- matrix(c(-.8, .1, .3, 1.1), ncol = 1L)
A <- matrix(.9, 1L, 1L)
lt <- log(c(1.4, .7, .8)); missing <- 2:4
ref <- loading_transport_reference(S, A)
W_fixed <- (U[missing, , drop = FALSE] - ref$mean[missing, , drop = FALSE]) %*%
  t(ref$precision_chol)
a_grid <- seq(.0001, 8, length.out = 1601L)
target <- vapply(a_grid, function(a) {
  ap <- matrix(a, 1L, 1L)
  rp <- loading_transport_reference(S, ap)
  up <- U
  up[missing, ] <- rp$mean[missing, , drop = FALSE] +
    t(backsolve(rp$precision_chol, t(W_fixed)))
  dnorm(a, 0, 1.3, log = TRUE) + loading_transport_score_loglik(S, up, ap, missing) +
    gp_loglik_integrated_general(y, X, up, lt)
}, numeric(1))
weights <- exp(target - max(target)); weights[c(1L, length(weights))] <-
  weights[c(1L, length(weights))] / 2
reference <- sum(a_grid * weights) / sum(weights)
evaluator <- make_gp_evaluator_general(y, X, calib_idx = 1L, U_obs = U)
trace_A <- numeric(6000L)
max_whiten_error <- 0
for (iteration in seq_along(trace_A)) {
  move <- update_loading_transport_ess(y, X, U, S, A, lt, missing,
    s_A = 1.3, ident = "lower_triangular", gp_evaluator = evaluator)
  U <- move$U; A <- move$A; trace_A[iteration] <- A[1L, 1L]
  ref <- loading_transport_reference(S, A)
  W_current <- (U[missing, , drop = FALSE] - ref$mean[missing, , drop = FALSE]) %*%
    t(ref$precision_chol)
  max_whiten_error <- max(max_whiten_error, abs(W_current - W_fixed))
}
kept <- trace_A[-seq_len(1000L)]
batch <- rowsum(matrix(kept, ncol = 1L), rep(seq_len(50L), each = 100L)) / 100
mcse <- sd(batch) / sqrt(nrow(batch))
stopifnot(abs(mean(kept) - reference) < max(.045, 7 * mcse),
          max_whiten_error < 1e-10)
before <- evaluator$counts()
state <- evaluator$evaluate(U, lt)
after <- evaluator$counts()
stopifnot(identical(state$U, U), after["gp_cache_hits"] == before["gp_cache_hits"] + 1L,
          after["gp_dense_evaluations"] == before["gp_dense_evaluations"],
          abs(gp_loglik_from_moments_general(state) -
                gp_loglik_integrated_general(y, X, U, lt)) < 1e-10)
set.seed(262)
cached_sigma <- sample_sigma2_eps_general(y, X, U, lt, gp_evaluator = evaluator)
set.seed(262)
dense_sigma <- sample_sigma2_eps_general(y, X, U, lt)
stopifnot(identical(cached_sigma, dense_sigma))
cat("Transport conditional reference, observed, MCSE:", reference, mean(kept), mcse,
    "; whiten error", max_whiten_error, "\n")
cat("Accepted transport GP cache and following exact sigma draw: PASS\n")
expect_error(update_loading_transport_ess(y, X, U, S, A, lt, missing, max_try = 0L),
             "exceeded max_try")

## Instrument every score-conditioned update: scores must be regenerated after
## a marginal row update, and calibrated U must never move at any transition.
run_small <- function(missing, use_schur) {
  env <- new.env(parent = globalenv())
  source(file.path(codes_dir, "00_parallel_utils.R"), local = env)
  source(file.path(codes_dir, "00_study2_functions.R"), local = env)
  n <- 8L
  X <- matrix(seq(-1, 1, length.out = n), ncol = 1L)
  U <- cbind(seq(-1.5, 1.5, length.out = n), sin(seq_len(n)))
  C <- cbind(rep(1:2, each = 4L), rep(c(1L, 2L, 3L, 2L), 2L))
  y <- X[, 1L] + U[, 1L]^2 + U[, 2L] / 3
  calibrated <- setdiff(seq_len(n), missing)
  needs_scores <- FALSE
  checks <- 0L
  original_joint <- env$update_measurement_joint_ess
  env$update_measurement_joint_ess <- function(...) {
    out <- original_joint(...); needs_scores <<- TRUE; out
  }
  original_scores <- env$sample_scores_ord
  env$sample_scores_ord <- function(...) {
    out <- original_scores(...); needs_scores <<- FALSE; out
  }
  original_transport <- env$update_loading_transport_ess
  needs_sigma <- FALSE
  env$update_loading_transport_ess <- function(...) {
    stopifnot(!needs_scores)
    out <- original_transport(...)
    stopifnot(all(out$U[calibrated, , drop = FALSE] == U[calibrated, , drop = FALSE]))
    checks <<- checks + 1L
    needs_sigma <<- TRUE
    out
  }
  original_sigma <- env$sample_sigma2_eps_general
  env$sample_sigma2_eps_general <- function(...) {
    out <- original_sigma(...); needs_sigma <<- FALSE; out
  }
  original_assert <- env$assert_ordprobit_state
  env$assert_ordprobit_state <- function(...) {
    stopifnot(!needs_sigma); original_assert(...)
  }
  env$fit_eivgp_ordprobit_fb(X, y, C, U_obs = U, calib_idx = calibrated, d = 2L,
    m_vec = c(2L, 3L), n_iter = 48L, burn = 16L, thin = 1L, n_chains = 1L,
    parallel_chains = FALSE, preset = "fast", seed = 971L,
    control_overrides = list(joint_measurement_every = 2L,
      loading_transport_every = 3L, gp_use_block_schur = use_schur))
}
for (missing in list(integer(0), 8L, c(2L, 4L, 8L), seq_len(8L))) {
  dense <- run_small(missing, FALSE)
  schur <- run_small(missing, TRUE)
  stopifnot(max(abs(dense$mcmc$samples_U - schur$mcmc$samples_U)) < 1e-8,
            max(abs(dense$mcmc$samples_logtheta - schur$mcmc$samples_logtheta)) < 1e-8,
            max(abs(dense$mcmc$samples_A - schur$mcmc$samples_A)) < 1e-8,
            all(is.finite(dense$mcmc$samples_sigma2)),
            dense$mcmc$chain_stats$joint_measurement_update_total == 24L,
            dense$mcmc$chain_stats$loading_transport_update_total == 16L)
}
cat("Combined-kernel zero/one/multiple missing, regeneration, dense/Schur: PASS\n")

control <- make_default_control_ordprobit(8L, 4L, "fast", 2L)
stopifnot(control$joint_measurement_every == 0L, control$loading_transport_every == 0L)
for (name in c("score_update_every", "A_update_every", "tau_update_every",
               "theta_update_every", "max_ess_try")) {
  bad <- control; bad[[name]] <- 0L
  expect_error(validate_control_ordprobit(bad, "interwoven", 4L, 100L), name)
}
for (name in c("score_update_every", "A_update_every", "tau_update_every")) {
  bad <- control; bad[[name]] <- 101L
  expect_error(validate_control_ordprobit(bad, "interwoven", 4L, 100L), name)
}
for (value in list(-1, NA_real_, .5, "2")) {
  bad <- control; bad$joint_measurement_every <- value
  expect_error(validate_control_ordprobit(bad, "interwoven", 4L, 100L),
               "joint_measurement_every")
}
bad <- control
bad$joint_theta_every <- bad$joint_global_every <- bad$joint_collapsed_every <- 0L
expect_error(validate_control_ordprobit(bad, "interwoven", 4L, 100L), "freezes all GP")
bad <- control
bad$joint_local_blocks_per_iter <- bad$joint_global_every <- bad$joint_collapsed_every <- 0L
bad$loading_transport_every <- 1L
expect_error(validate_control_ordprobit(bad, "interwoven", 4L, 100L), "freezes all missing-U")
bad <- control
bad$joint_theta_every <- bad$joint_global_every <- bad$joint_collapsed_every <- 200L
expect_error(validate_control_ordprobit(bad, "interwoven", 4L, 100L), "freezes all GP")
bad <- control; bad$joint_measurement_every <- 2L
expect_error(validate_control_ordprobit(bad, "legacy", 4L, 100L), "require interwoven")
for (overrides in list(list(1L), list(joint_measurement_every = 1L,
                                    joint_measurement_every = 2L))) {
  expect_error(fit_eivgp_ordprobit_fb(matrix(1:4, ncol = 1L), c(1, 4, 2, 3),
    matrix(c(1L, 1L, 2L, 2L), ncol = 1L), d = 1L, n_iter = 10L, burn = 2L,
    n_chains = 1L, parallel_chains = FALSE, control_overrides = overrides),
    "unique, nonempty names")
}
cat("Schedule safety and opt-in defaults: PASS\n")
