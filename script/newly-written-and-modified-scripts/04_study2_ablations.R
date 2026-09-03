############################################################
## 04_study2_ablations.R
##
## Appendix-only ablations for Study II.
##
## PI-GP and CC-GP are diagnostic decompositions of EIV-GP, not published
## competitors.  Their ordinal-probit measurement model is fitted without y.
## Full-U GP is an infeasible fitted benchmark that observes every training U.
############################################################

if (!exists("sample_scores_ord")) source("00_study2_functions.R")

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
