## General-purpose adapters for the 0.3.1 collapsed sampler.
## Existing prediction helpers use an algebraically equivalent covariance:
## sigma2 = V * (1-r), rho^2 = r/(1-r).

mixedgp_v030_validate_arguments <- function(n_iter, burn, thin, n_chains,
                                           seed, parallel_chains, n_cores,
                                           verbose, preset, dots) {
  mixedgp_validate_named_dots(dots)
  if (length(dots)) {
    stop("Unsupported 0.3.1 sampler argument(s): ", paste(names(dots), collapse = ", "),
         ". Configure the new model with priors and sampler_control; historical sampler controls are incompatible.")
  }
  n_iter <- mixedgp_as_integer_strict(n_iter, "n_iter", 1L, 1L)
  burn <- mixedgp_as_integer_strict(burn, "burn", 0L, 1L)
  thin <- mixedgp_as_integer_strict(thin, "thin", 1L, 1L)
  n_chains <- mixedgp_as_integer_strict(n_chains, "n_chains", 1L, 1L)
  seed <- mixedgp_as_integer_strict(seed, "seed", 0L, 1L)
  if (n_iter <= burn) stop("n_iter must exceed burn.")
  if (thin != 1L) stop("Every post-warm-up draw is retained; thin must be 1.")
  parallel_chains <- mixedgp_validate_flag(parallel_chains, "parallel_chains")
  verbose <- mixedgp_validate_flag(verbose, "verbose")
  if (!is.null(n_cores)) n_cores <- mixedgp_as_integer_strict(n_cores, "n_cores", 1L, 1L)
  if (!is.null(preset)) {
    if (!is.character(preset) || length(preset) != 1L || is.na(preset) ||
        !preset %in% c("fast", "balanced", "robust", "thorough")) {
      stop("preset must be NULL or a recognized historical budget label.")
    }
    warning("preset is a deprecated budget label; 0.3.1 uses explicit n_iter, burn, priors, and sampler_control.",
            call. = FALSE)
  }
  list(n_iter = n_iter, burn = burn, n_chains = n_chains, seed = seed,
       parallel_chains = parallel_chains, n_cores = n_cores, verbose = verbose)
}

mixedgp_v030_standardize_xy <- function(X, y) {
  y <- as_numeric_matrix_strict(y, "y_raw")
  if (ncol(y) != 1L || nrow(y) < 2L) {
    stop("y_raw must contain at least two finite univariate observations.")
  }
  y <- as.numeric(y[, 1L])
  X <- as_numeric_matrix_strict(X, "X_raw", nrow_expected = length(y))
  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
  X_center <- colMeans(X)
  X_scale <- apply(X, 2L, stats::sd)
  y_center <- mean(y)
  y_scale <- stats::sd(y)
  if (any(!is.finite(X_scale) | X_scale <= 0)) {
    stop("Every numeric predictor must have positive finite variation.")
  }
  if (!is.finite(y_scale) || y_scale <= 0) stop("y_raw must have positive finite variation.")
  list(X_raw = X, X = sweep(sweep(X, 2L, X_center, "-"), 2L, X_scale, "/"),
       X_center = X_center, X_scale = X_scale, y_raw = y,
       y = (y - y_center) / y_scale, y_center = y_center, y_scale = y_scale)
}

mixedgp_v030_calibration <- function(U_obs, n, d, calib_idx = NULL,
                                     standardize_U = FALSE) {
  standardize_U <- mixedgp_validate_flag(standardize_U, "standardize_U")
  explicit_indices <- !is.null(calib_idx)
  if (explicit_indices) {
    calib_idx <- if (length(calib_idx)) mixedgp_as_integer_strict(calib_idx, "calib_idx", 1L) else integer(0)
    if (anyDuplicated(calib_idx) || any(calib_idx > n)) stop("Invalid calib_idx.")
    calib_idx <- sort(calib_idx)
  }
  full <- matrix(NA_real_, n, d, dimnames = list(NULL, paste0("u", seq_len(d))))
  if (!is.null(U_obs)) {
    U_obs <- as_numeric_matrix_strict(U_obs, "U_obs", ncol_expected = d, allow_na = TRUE)
    if (!is.null(colnames(U_obs))) colnames(full) <- colnames(U_obs)
    if (nrow(U_obs) == n) {
      count <- rowSums(is.finite(U_obs))
      if (any(count > 0L & count < d)) stop("Each U_obs row must contain every latent coordinate or only NA.")
      observed <- which(count == d)
      if (explicit_indices && !all(calib_idx %in% observed)) {
        stop("Every row selected by calib_idx must be completely observed in U_obs.")
      }
      if (!explicit_indices) calib_idx <- observed
      ## Low-level callers may supply a complete calibration matrix and select
      ## its observed rows explicitly. Values outside that selection are hidden.
      full[calib_idx, ] <- U_obs[calib_idx, , drop = FALSE]
    } else if (explicit_indices && nrow(U_obs) == length(calib_idx) && all(is.finite(U_obs))) {
      full[calib_idx, ] <- U_obs
    } else stop("U_obs must have n rows, or one complete row per calib_idx.")
  } else if (explicit_indices && length(calib_idx)) {
    stop("calib_idx requires observed latent values.")
  }
  if (is.null(calib_idx)) calib_idx <- integer(0)
  anchor <- mixedgp_latent_anchor_status(full, calib_idx, d)
  center <- rep(0, d)
  scale <- rep(1, d)
  if (standardize_U && length(calib_idx)) {
    center <- colMeans(full[calib_idx, , drop = FALSE])
    scale <- apply(full[calib_idx, , drop = FALSE], 2L, stats::sd)
    scale[!is.finite(scale) | scale <= 0] <- 1
  }
  model <- full
  if (length(calib_idx)) {
    model[calib_idx, ] <- sweep(sweep(full[calib_idx, , drop = FALSE],
                                    2L, center, "-"), 2L, scale, "/")
  }
  list(U_obs = model, U_obs_raw = full, U_center = center, U_scale = scale,
       U_names = colnames(full), standardize_U = standardize_U,
       calib_idx = calib_idx, miss_idx = setdiff(seq_len(n), calib_idx),
       latent_scale_anchored = anchor$anchored,
       latent_anchor_rank = anchor$affine_rank,
       latent_anchor_required_rank = anchor$required_rank)
}

mixedgp_v030_join_draws <- function(parts) {
  if (is.null(parts[[1L]])) return(NULL)
  if (length(dim(parts[[1L]])) == 3L) return(combine_chain_arrays(parts))
  if (is.matrix(parts[[1L]])) return(do.call(rbind, parts))
  unlist(parts, use.names = FALSE)
}

mixedgp_v030_append_chains <- function(chains, resume) {
  if (is.null(resume)) return(chains)
  if (length(resume$old_chains) != length(chains)) stop("Checkpoint chain count changed.")
  for (i in seq_along(chains)) {
    old <- resume$old_chains[[i]]
    for (name in names(old)[startsWith(names(old), "samples_")]) {
      if (is.null(old[[name]])) next
      if (is.null(chains[[i]][[name]])) stop("Checkpoint draw fields changed: ", name)
      chains[[i]][[name]] <- mixedgp_v030_join_draws(list(old[[name]], chains[[i]][[name]]))
    }
    if (!is.null(old$initial_state)) chains[[i]]$initial_state <- old$initial_state
    counters <- setdiff(names(old$stats), c("chain", "seed", "mean_abs_angle"))
    for (name in counters) chains[[i]]$stats[[name]] <-
      old$stats[[name]] + chains[[i]]$stats[[name]]
    if (all(c("abs_angle_sum", "angle_count") %in% names(chains[[i]]$stats))) {
      chains[[i]]$stats$mean_abs_angle <- chains[[i]]$stats$abs_angle_sum /
        max(1, chains[[i]]$stats$angle_count)
    } else if ("mean_abs_angle" %in% names(old$stats)) {
      before <- resume$iteration
      added <- resume$target_iteration - before
      chains[[i]]$stats$mean_abs_angle <- (old$stats$mean_abs_angle * before +
        chains[[i]]$stats$mean_abs_angle * added) / (before + added)
    }
  }
  chains
}

mixedgp_v030_raw_series <- function(fit) {
  engine <- mixedgp_fit_engine(fit)
  raw <- if (engine == "univariate") mixedgp_study1_raw_series(fit) else mixedgp_study2_raw_series(fit)
  ## These compatibility quantities are derived from V, r and dictionary states.
  raw <- raw[!grepl("^(rho$|sigma_epsilon$|theta_u)", names(raw))]
  chains <- fit$mcmc$samples_by_chain
  raw <- c(list(V = chains$V, r = chains$r), raw)
  dictionary <- fit$priors$u_dictionary
  for (k in seq_len(ncol(chains$J[[1L]]))) {
    size <- if (is.matrix(dictionary)) nrow(dictionary) else length(dictionary[[k]])
    if (size > 1L) raw[[paste0("J[", k, "]")]] <- lapply(chains$J, function(z) as.numeric(z[, k]))
    ## A fixed singleton is known from the model; any unvisited nontrivial
    ## dictionary state remains present and cannot be silently screened out.
    if (size > 1L) for (j in seq_len(size)) {
      raw[[paste0("dictionary_occupancy[", k, ",", j, "]")]] <-
        lapply(chains$J, function(z) as.numeric(z[, k] == j))
    }
  }
  if (length(chains$pi) && !is.null(chains$pi[[1L]])) {
    for (j in seq_len(ncol(chains$pi[[1L]]))) raw[[paste0("category_probability[", j, "]")]] <-
      lapply(chains$pi, function(z) z[, j])
  }
  raw
}

mixedgp_v030_checkpoint_hash <- function(object) {
  cp <- object$checkpoint
  cp$signature <- NULL
  path <- tempfile("eivgp-v030-checkpoint-", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(list(data = object$data, control = object$control, priors = object$priors,
               kernel = object$kernel, mcmc = object$mcmc, checkpoint = cp,
               interface = object$interface, model_specification = object$model_specification),
          path, compress = FALSE, version = 2L)
  unname(tools::md5sum(path))
}

mixedgp_v030_fit_object <- function(result, data, engine, kernel, matern_nu,
                                    resume = NULL) {
  if (!identical(result$sampler_version, "0.3.1")) stop("Unexpected sampler version.")
  if (!is.null(resume)) resume$target_iteration <- result$iteration
  chains <- result$chains
  m_vec <- if (engine == "univariate") data$m else data$m_vec
  tau_names <- if (engine == "univariate") paste0("tau", seq_len(m_vec - 1L)) else
    unlist(lapply(seq_along(m_vec), function(j) paste0("tau[", j, ",", seq_len(m_vec[j] - 1L), "]")),
           use.names = FALSE)
  for (i in seq_along(chains)) colnames(chains[[i]]$samples_tau) <- tau_names
  if (engine == "univariate") for (i in seq_along(chains)) {
    dims <- dim(chains[[i]]$samples_U)
    chains[[i]]$samples_u <- matrix(chains[[i]]$samples_U, dims[1L], dims[2L])
  }
  chains <- mixedgp_v030_append_chains(chains, resume)
  fields <- c(if (engine == "univariate") "u" else c("U", "S", "A"),
              "tau", "logtheta", "sigma2", "V", "r", "J", "pi")
  by_chain <- setNames(lapply(fields, function(name) {
    lapply(chains, function(chain) chain[[paste0("samples_", name)]])
  }), fields)
  mcmc <- setNames(lapply(by_chain, mixedgp_v030_join_draws), paste0("samples_", fields))
  mcmc$samples_by_chain <- by_chain
  mcmc$chain_initial <- lapply(chains, function(chain) chain$initial_state)
  lengths <- vapply(by_chain$logtheta, nrow, integer(1L))
  mcmc$mcmc_draw_info <- data.frame(
    chain = rep(seq_along(chains), times = lengths),
    draw_within_chain = unlist(lapply(lengths, seq_len), use.names = FALSE))
  mcmc$chain_stats <- do.call(rbind, lapply(chains, function(chain) chain$stats))
  fit <- list(data = data, kernel = list(name = kernel, matern_nu = matern_nu),
              priors = result$priors, gp_prior = result$priors, control = result$control,
              sampler_strategy = "collapsed_dictionary_ess",
              sampler_version = "0.3.1", mcmc = mcmc,
              checkpoint = list(version = 2L, sampler_version = "0.3.1",
                iteration = result$iteration, burn = result$burn, thin = 1L,
                control = result$control,
                states = lapply(chains, function(chain) chain$checkpoint)))
  raw <- mixedgp_v030_raw_series(fit)
  table <- suppressWarnings(mixedgp_summarize_diagnostic_series(raw))
  get_rows <- function(pattern) table[grepl(pattern, table$parameter), , drop = FALSE]
  hyper <- get_rows("^(V$|r$|theta_x|J\\[|dictionary_occupancy)")
  tau <- get_rows("^tau")
  latent <- get_rows("^[uU]\\[")
  latent$global_index <- as.integer(sub("^[uU]\\[([0-9]+).*$", "\\1", latent$parameter))
  if (engine == "multivariate") {
    latent$coord <- as.integer(sub("^U\\[[0-9]+,([0-9]+)\\]$", "\\1", latent$parameter))
  }
  finite_summary <- function(x, f) if (length(x) && all(is.finite(x))) f(x) else NA_real_
  diagnostic_summary <- data.frame(
    sampler_version = "0.3.1", sampler_strategy = "collapsed_dictionary_ess",
    kernel = kernel, matern_nu = matern_nu, n_chains = length(chains),
    parallel_backend = result$parallel_backend, parallel_cores = result$parallel_cores,
    n_iter = result$iteration, burn = result$burn, thin = 1L,
    saved_per_chain = lengths[1L], total_saved_draws = sum(lengths),
    max_rhat_hyper = finite_summary(hyper$rhat, max),
    max_rhat_tau = finite_summary(tau$rhat, max),
    min_ess_key = finite_summary(table$ess, min),
    covariance_jitter = 0, forms_explicit_covariance_inverse = FALSE,
    time_seconds = as.numeric(result$time_seconds))
  diagnostic_summary[[if (engine == "univariate") "max_rhat_missing_u" else "max_rhat_missing_U"]] <-
    finite_summary(latent$rhat, max)
  diagnostic_summary[[if (engine == "univariate") "median_rhat_missing_u" else "median_rhat_missing_U"]] <-
    finite_summary(latent$rhat, stats::median)
  stats <- mcmc$chain_stats
  for (name in intersect(c("gp_full_factorizations", "gp_block_setups", "gp_block_evaluations",
                           "dictionary_evaluations", "dictionary_switches"), names(stats))) {
    diagnostic_summary[[name]] <- sum(stats[[name]])
  }
  fit$diagnostics <- list(rhat_hyper = hyper, rhat_tau = tau, ess_key = table,
                          summary = diagnostic_summary)
  fit$diagnostics[[if (engine == "univariate") "rhat_u" else "rhat_U"]] <- latent
  if (engine == "multivariate") {
    fit$diagnostics$rhat_A <- get_rows("^A\\[")
    fit$diagnostics$summary$max_rhat_A <- finite_summary(fit$diagnostics$rhat_A$rhat, max)
  }
  fit
}

fit_eivgp_1d_v030 <- function(x_raw, y_raw, c_ord, u_true = NULL,
                              calib_idx = integer(0), m = 6L, tau_true = NULL,
                              n_iter = 1750L, burn = 500L, thin = 1L,
                              n_chains = 4L, preset = NULL, seed = 1L,
                              parallel_chains = TRUE, n_cores = NULL,
                              verbose = FALSE, kernel = c("se", "matern"),
                              matern_nu = 2.5, u_obs = NULL, priors = list(),
                              sampler_control = list(), .resume = NULL, ...) {
  args <- mixedgp_v030_validate_arguments(n_iter, burn, thin, n_chains, seed,
    parallel_chains, n_cores, verbose, preset, list(...))
  xy <- mixedgp_v030_standardize_xy(x_raw, y_raw)
  n <- length(xy$y)
  m <- mixedgp_as_integer_strict(m, "m", 2L, 1L)
  C <- prepare_ordinal_matrix(c_ord, m_vec = m, name = "c_ord")$C
  if (nrow(C) != n || ncol(C) != 1L) stop("c_ord must have one ordinal code per observation.")
  kernel <- match.arg(kernel)
  if (is.null(u_obs) && length(calib_idx)) {
    if (is.null(u_true) || length(u_true) != n) stop("Supply calibration values through u_obs.")
    u_obs <- rep(NA_real_, n)
    u_obs[calib_idx] <- u_true[calib_idx]
  }
  calibration <- mixedgp_v030_calibration(u_obs, n, 1L,
    if (length(calib_idx)) calib_idx else NULL)
  args$verbose <- NULL
  if (verbose) cat("Running collapsed dictionary EIV-GP (threshold measurement).\n")
  result <- do.call(mixedgp_v030_fit, c(list(X = xy$X, y = xy$y, C = C,
    U_obs = calibration$U_obs, m_vec = m, measurement = "threshold", ident = "none",
    kernel = kernel, matern_nu = matern_nu, priors = priors,
    sampler_control = sampler_control, store_scores = FALSE, .resume = .resume), args))
  data <- list(x_raw = xy$X_raw, x = xy$X, x_center = xy$X_center, x_scale = xy$X_scale,
    y_raw = xy$y_raw, y = xy$y, y_center = xy$y_center, y_scale = xy$y_scale,
    c_ord = as.integer(C[, 1L]), u_obs = as.numeric(calibration$U_obs),
    u_true = u_true, tau_true = tau_true, m = m, p_x = ncol(xy$X),
    predictor_names = colnames(xy$X), theta_spec = list(
      p_x = ncol(xy$X), x_index = 1L + seq_len(ncol(xy$X)),
      u_index = 2L + ncol(xy$X)))
  data <- c(data, calibration[setdiff(names(calibration), "U_obs")])
  mixedgp_v030_fit_object(result, data, "univariate", kernel, matern_nu, .resume)
}

fit_eivgp_ordprobit_v030 <- function(X_raw, y_raw, C_ord, U_obs = NULL,
                                     calib_idx = NULL, U_true_eval = NULL,
                                     d = 2L, m_vec = NULL,
                                     ident = c("lower_triangular", "none"),
                                     n_iter = 1750L, burn = 500L, thin = 1L,
                                     n_chains = 4L, preset = NULL, seed = 1L,
                                     parallel_chains = TRUE, n_cores = NULL,
                                     verbose = FALSE, kernel = c("se", "matern"),
                                     matern_nu = 2.5, standardize_U = FALSE,
                                     ordinal_levels = NULL, store_scores = FALSE,
                                     priors = list(), sampler_control = list(),
                                     .resume = NULL, ...) {
  args <- mixedgp_v030_validate_arguments(n_iter, burn, thin, n_chains, seed,
    parallel_chains, n_cores, verbose, preset, list(...))
  xy <- mixedgp_v030_standardize_xy(X_raw, y_raw)
  n <- length(xy$y)
  if (missing(d) && !is.null(U_obs)) d <- ncol(as.matrix(U_obs))
  d <- mixedgp_as_integer_strict(d, "d", 1L, 1L)
  ordinal <- prepare_ordinal_matrix(C_ord, m_vec = m_vec, level_maps = ordinal_levels, name = "C_ord")
  if (nrow(ordinal$C) != n) stop("C_ord must have one row per observation.")
  ident <- match.arg(ident)
  if (ident == "lower_triangular" && ncol(ordinal$C) < d) stop("lower_triangular requires q >= d.")
  kernel <- match.arg(kernel)
  store_scores <- mixedgp_validate_flag(store_scores, "store_scores")
  calibration <- mixedgp_v030_calibration(U_obs, n, d, calib_idx, standardize_U)
  args$verbose <- NULL
  if (verbose) cat("Running collapsed dictionary EIV-GP (ordinal-probit measurement).\n")
  result <- do.call(mixedgp_v030_fit, c(list(X = xy$X, y = xy$y, C = ordinal$C,
    U_obs = calibration$U_obs, m_vec = ordinal$m_vec, measurement = "probit", ident = ident,
    kernel = kernel, matern_nu = matern_nu, priors = priors,
    sampler_control = sampler_control, store_scores = store_scores, .resume = .resume), args))
  truth <- if (is.null(U_true_eval)) NULL else as_numeric_matrix_strict(
    U_true_eval, "U_true_eval", nrow_expected = n, ncol_expected = d)
  truth_model <- if (is.null(truth)) NULL else sweep(sweep(truth, 2L,
    calibration$U_center, "-"), 2L, calibration$U_scale, "/")
  data <- c(xy, calibration, list(C_ord = ordinal$C, C_level_maps = ordinal$level_maps,
    C_names = colnames(ordinal$C), U_true_eval = truth_model, U_true_eval_raw = truth,
    m_vec = ordinal$m_vec, q = ncol(ordinal$C), d = d, p = ncol(xy$X), ident = ident))
  mixedgp_v030_fit_object(result, data, "multivariate", kernel, matern_nu, .resume)
}
