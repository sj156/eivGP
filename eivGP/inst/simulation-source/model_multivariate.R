############################################################
## model_multivariate.R
##
## General-purpose fully Bayesian ordinal-probit EIV-GP code,
## plus Study II synthetic data utilities.
############################################################

if (!exists("mixedgp_parallel_lapply", mode = "function")) {
  this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  sibling_dir <- if (is.null(this_file)) character(0) else {
    dirname(normalizePath(this_file))
  }
  parallel_utility <- c(
    file.path(sibling_dir, "core_parallel.R"),
    "core_parallel.R",
    file.path("codes", "core_parallel.R")
  )
  parallel_utility <- parallel_utility[file.exists(parallel_utility)][1L]
  if (is.na(parallel_utility)) {
    stop("Source core_parallel.R before model_multivariate.R.")
  }
  sys.source(parallel_utility, envir = environment())
}
if (!exists("pairwise_sqdist", mode = "function")) {
  this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  sibling_dir <- if (is.null(this_file)) character(0) else {
    dirname(normalizePath(this_file))
  }
  numeric_utility <- c(
    file.path(sibling_dir, "core_numerics.R"),
    "core_numerics.R",
    file.path("codes", "core_numerics.R")
  )
  numeric_utility <- numeric_utility[file.exists(numeric_utility)][1L]
  if (is.na(numeric_utility)) {
    stop("Source core_numerics.R before model_multivariate.R.")
  }
  sys.source(numeric_utility, envir = environment())
}

require_study2_reporting_packages <- function(packages, context) {
  missing <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing) > 0L) {
    stop(
      "The ", context, " requires optional package(s): ",
      paste(missing, collapse = ", ")
    )
  }
  invisible(TRUE)
}

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
      
      rep_n_iter = 800L,
      rep_burn = 300L,
      rep_thin = 2L,
      rep_n_chains = 2L,
      
      mc_n_iter = 700L,
      mc_burn = 300L,
      mc_thin = 2L,
      mc_n_chains = 1L,
      
      n_pred_draw = 120L,
      n_density_draw = 200L,
      predictive_latent_sampler = "minimax_tilting",
      diagnostic_n_new_latent_gibbs = 100L,
      rejection_max_batches = 1000L,
      n_oracle_pool = 50000L,
      n_starts_learned = 2L,
      preset = "fast"
    ))
  }
  
  if (config == "balanced") {
    return(list(
      config = config,
      n_test = 400L,
      n_rep = 10L,
      
      rep_n_iter = 6000L,
      rep_burn = 2000L,
      rep_thin = 1L,
      rep_n_chains = 12L,
      
      mc_n_iter = 2500L,
      mc_burn = 800L,
      mc_thin = 1L,
      mc_n_chains = 12L,
      
      n_pred_draw = 400L,
      n_density_draw = 800L,
      predictive_latent_sampler = "minimax_tilting",
      diagnostic_n_new_latent_gibbs = 100L,
      rejection_max_batches = 1000L,
      n_oracle_pool = 150000L,
      n_starts_learned = 4L,
      preset = "balanced"
    ))
  }
  
  list(
    config = config,
    n_test = 500L,
    n_rep = 50L,
    
    rep_n_iter = 6000L,
    rep_burn = 2000L,
    rep_thin = 4L,
    rep_n_chains = 12L,
    
    mc_n_iter = 5000L,
    mc_burn = 1500L,
    mc_thin = 5L,
    mc_n_chains = 12L,
    
    n_pred_draw = 600L,
    n_density_draw = 1200L,
    predictive_latent_sampler = "minimax_tilting",
    diagnostic_n_new_latent_gibbs = 100L,
    rejection_max_batches = 1000L,
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
    theta_update_every <- 8L
    n_u_blocks_per_iter <- 1L
    u_block_size <- min(n_mis, max(8L, ceiling(sqrt(max(n_mis, 1L)))))
    global_u_every <- 50L
    joint_local_blocks_per_iter <- 1L
    joint_global_every <- 10L
    joint_collapsed_every <- 10L
    marginal_measurement_every <- 5L
  }
  
  if (preset == "balanced") {
    theta_update_every <- 4L
    n_u_blocks_per_iter <- 2L
    u_block_size <- min(n_mis, max(10L, ceiling(1.5 * sqrt(max(n_mis, 1L)))))
    global_u_every <- 25L
    joint_local_blocks_per_iter <- 2L
    joint_global_every <- 5L
    joint_collapsed_every <- 5L
    marginal_measurement_every <- 2L
  }
  
  if (preset == "thorough") {
    theta_update_every <- 2L
    n_u_blocks_per_iter <- 3L
    u_block_size <- min(n_mis, max(12L, ceiling(2.0 * sqrt(max(n_mis, 1L)))))
    global_u_every <- 15L
    joint_local_blocks_per_iter <- 3L
    joint_global_every <- 3L
    joint_collapsed_every <- 3L
    marginal_measurement_every <- 1L
  }
  
  list(
    preset = preset,
    
    ## Main MCMC knobs.
    theta_update_every = theta_update_every,
    n_u_blocks_per_iter = n_u_blocks_per_iter,
    u_block_size = u_block_size,
    global_u_every = global_u_every,

    ## Exact partially collapsed/interwoven transition frequencies.
    joint_local_blocks_per_iter = joint_local_blocks_per_iter,
    joint_theta_every = 1L,
    joint_global_every = joint_global_every,
    joint_collapsed_every = joint_collapsed_every,
    marginal_measurement_every = marginal_measurement_every,
    
    ## Ordinal-probit update frequencies.
    score_update_every = 1L,
    A_update_every = 1L,
    tau_update_every = 1L,
    
    ## Slice sampler knobs for GP hyperparameters.
    theta_slice_width_init = rep(0.8, 1L + d + 10L), # trimmed inside fitter
    adapt_theta_width = TRUE,
    adapt_every = 100L,
    adapt_window = 500L,
    theta_width_min = 0.15,
    theta_width_max = 2.50,
    
    ## Measurement-model prior and constraints.
    tau_bound = 6,
    s_A = 2.5,
    
    ## Elliptical-slice safeguard.
    max_ess_try = 5000L
  )
}

############################################################
## General utilities
############################################################

as_numeric_matrix_strict <- function(x,
                                     name,
                                     nrow_expected = NULL,
                                     ncol_expected = NULL,
                                     allow_na = FALSE) {
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)

  if (is.data.frame(x)) {
    is_numeric_column <- vapply(
      x,
      function(z) is.numeric(z) || is.integer(z) || is.logical(z),
      logical(1)
    )
    if (any(!is_numeric_column)) {
      stop(
        name,
        " must contain only numeric columns; encode factors before fitting."
      )
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
  if (isTRUE(allow_na)) {
    if (any(is.infinite(out))) stop(name, " contains infinite values.")
  } else if (any(!is.finite(out))) {
    stop(name, " contains missing or non-finite values.")
  }

  out
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

rtruncnorm_one <- function(mean, sd, lower = -Inf, upper = Inf) {
  rtruncnorm_vec(mean, sd, lower, upper)[1]
}

maximin_lhs_nd <- function(n,
                           d,
                           lower = -1,
                           upper = 1,
                           n_starts = 50L) {
  n <- as.integer(n)
  d <- as.integer(d)
  n_starts <- as.integer(n_starts)
  if (n < 1L || d < 1L || n_starts < 1L) {
    stop("n, d, and n_starts must be positive integers.")
  }

  random_lhs <- function() {
    X <- matrix(NA_real_, n, d)
    for (j in seq_len(d)) {
      z <- (seq_len(n) - runif(n)) / n
      X[, j] <- sample(z)
    }
    X
  }

  best <- NULL
  best_min_distance <- -Inf
  for (start in seq_len(n_starts)) {
    candidate <- random_lhs()
    min_distance <- if (n == 1L) {
      Inf
    } else {
      min(stats::dist(candidate))
    }
    if (min_distance > best_min_distance) {
      best <- candidate
      best_min_distance <- min_distance
    }
  }

  lower + (upper - lower) * best
}

pattern_key <- function(C) {
  C <- as.matrix(C)
  apply(C, 1, paste, collapse = "_")
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
## Tau helpers for proxy-specific level counts
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

bounded_slice_update <- function(x0, logf, w = 1,
                                 lower = -Inf, upper = Inf,
                                 max_steps_out = 50,
                                 max_iter = 200,
                                 fail_on_limit = FALSE) {
  if (upper <= lower) {
    if (isTRUE(fail_on_limit)) stop("Slice interval is empty.")
    return(list(x = x0, n_eval = 0L))
  }
  
  eps <- 1e-12
  if (is.finite(lower)) x0 <- max(x0, lower + eps)
  if (is.finite(upper)) x0 <- min(x0, upper - eps)
  
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
    
    ## Add small smoothing to handle empty levels at initialization.
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
  q <- ncol(C)
  
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

assert_ordprobit_state <- function(C,
                                   U,
                                   S,
                                   A,
                                   tau,
                                   m_vec,
                                   ident,
                                   calib_idx = integer(0),
                                   U_obs = NULL,
                                   logtheta = NULL,
                                   sigma2_eps = NULL,
                                   tolerance = 1e-10) {
  C <- as.matrix(C)
  U <- as.matrix(U)
  S <- as.matrix(S)
  A <- as.matrix(A)
  n <- nrow(C)
  q <- ncol(C)
  d <- ncol(U)
  if (!identical(dim(S), c(n, q)) || nrow(U) != n ||
      !identical(dim(A), c(q, d)) || length(tau) != q ||
      any(!is.finite(c(U, S, A)))) {
    stop("Ordinal-probit sampler state has incompatible or nonfinite values.")
  }
  for (j in seq_len(q)) {
    tau_j <- as.numeric(tau[[j]])
    if (length(tau_j) != m_vec[j] - 1L || any(!is.finite(tau_j)) ||
        is.unsorted(tau_j, strictly = TRUE)) {
      stop("Ordinal-probit cut points are not finite and strictly ordered.")
    }
    lower <- c(-Inf, tau_j)[C[, j]]
    upper <- c(tau_j, Inf)[C[, j]]
    if (any(S[, j] < lower - tolerance | S[, j] > upper + tolerance)) {
      stop("Ordinal-probit scores left their observed-category rectangles.")
    }
  }
  if (ident == "lower_triangular") {
    for (j in seq_len(d)) {
      if (A[j, j] <= 0) stop("A constrained diagonal loading is not positive.")
      if (j < d && any(abs(A[j, (j + 1L):d]) > tolerance)) {
        stop("A violates the lower-triangular identification constraint.")
      }
    }
  }
  if (length(calib_idx) > 0L) {
    if (is.null(U_obs) ||
        max(abs(U[calib_idx, , drop = FALSE] -
                  U_obs[calib_idx, , drop = FALSE])) > tolerance) {
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
      
      if (!is.finite(L) || !is.finite(U) || L >= U) {
        stop("Cut-point Gibbs update encountered an empty conditional interval.")
      }
      tau_j[r] <- runif(1, L, U)
      
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
## General GP likelihood and prediction
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
  
  sum(dnorm(
    logtheta,
    mean = gp_prior$mean,
    sd = gp_prior$sd,
    log = TRUE
  ))
}

weighted_sqdist_general <- function(X_a,
                                    U_a,
                                    theta_x,
                                    theta_u,
                                    X_b = NULL,
                                    U_b = NULL) {
  X_a <- as.matrix(X_a)
  U_a <- as.matrix(U_a)
  if (is.null(X_b)) X_b <- X_a
  if (is.null(U_b)) U_b <- U_a
  X_b <- as.matrix(X_b)
  U_b <- as.matrix(U_b)

  if (nrow(X_a) != nrow(U_a) || nrow(X_b) != nrow(U_b)) {
    stop("X and U must have the same number of rows within each data set.")
  }
  if (ncol(X_a) != ncol(X_b) || ncol(U_a) != ncol(U_b)) {
    stop("Training and prediction inputs have incompatible dimensions.")
  }
  if (length(theta_x) != ncol(X_a) || length(theta_u) != ncol(U_a)) {
    stop("Kernel parameter dimensions do not match X and U.")
  }

  D2 <- matrix(0, nrow(X_a), nrow(X_b))
  for (j in seq_len(ncol(X_a))) {
    D2 <- D2 + theta_x[j] * pairwise_sqdist(
      X_a[, j, drop = FALSE],
      X_b[, j, drop = FALSE]
    )
  }
  for (k in seq_len(ncol(U_a))) {
    D2 <- D2 + theta_u[k] * pairwise_sqdist(
      U_a[, k, drop = FALSE],
      U_b[, k, drop = FALSE]
    )
  }

  pmax(D2, 0)
}

gp_corr_general <- function(X,
                            U,
                            logtheta,
                            kernel = "se",
                            matern_nu = 2.5) {
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
  
  D2 <- weighted_sqdist_general(X, U, theta_x, theta_u)
  kernel_spec <- normalize_gp_kernel(kernel, matern_nu)
  
  list(
    R = kernel_from_weighted_sqdist(
      D2,
      kernel = kernel_spec$name,
      matern_nu = kernel_spec$matern_nu
    ),
    rho = rho,
    theta_x = theta_x,
    theta_u = theta_u,
    kernel = kernel_spec
  )
}

gp_state_general <- function(y,
                             X,
                             U,
                             logtheta,
                             sigma2_eps,
                             kernel = "se",
                             matern_nu = 2.5) {
  n <- length(y)
  
  cc <- gp_corr_general(X, U, logtheta, kernel, matern_nu)
  A_mat <- diag(n) + cc$rho^2 * cc$R
  
  Uchol <- safe_chol(A_mat)
  Ainv_y <- solve_chol(Uchol, y)
  
  logdetA <- 2 * sum(log(diag(Uchol)))
  quad <- sum(y * Ainv_y)
  
  loglik <- -0.5 * (
    n * log(2 * base::pi * sigma2_eps) +
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

theta_logpost_integrated_general <- function(y,
                                             X,
                                             U,
                                             logtheta,
                                             gp_prior,
                                             kernel = "se",
                                             matern_nu = 2.5) {
  lp <- log_prior_logtheta_gp(logtheta, gp_prior)
  if (!is.finite(lp)) return(-Inf)

  lp + gp_loglik_integrated_general(
    y, X, U, logtheta,
    kernel = kernel,
    matern_nu = matern_nu
  )
}

gp_loglik_integrated_general <- function(y,
                                         X,
                                         U,
                                         logtheta,
                                         kernel = "se",
                                         matern_nu = 2.5) {
  n <- length(y)
  
  cc <- gp_corr_general(X, U, logtheta, kernel, matern_nu)
  A_mat <- diag(n) + cc$rho^2 * cc$R
  
  ## A_mat is positive definite under the stated covariance model. A failure
  ## here is therefore a numerical failure, not an out-of-support proposal;
  ## abort rather than silently changing the transition kernel.
  Uchol <- safe_chol(A_mat)
  
  Ainv_y <- solve_chol(Uchol, y)
  
  logdetA <- 2 * sum(log(diag(Uchol)))
  quad <- sum(y * Ainv_y)
  
  -0.5 * logdetA -
    (a_eps0 + n / 2) * log(b_eps0 + 0.5 * quad)
}

logspace_sub <- function(log_x, log_y) {
  log_x <- as.numeric(log_x)
  log_y <- as.numeric(log_y)
  if (length(log_x) != length(log_y)) {
    stop("log_x and log_y must have the same length.")
  }

  out <- rep(-Inf, length(log_x))
  only_x <- is.finite(log_x) & !is.finite(log_y)
  out[only_x] <- log_x[only_x]
  both <- is.finite(log_x) & is.finite(log_y) & log_y < log_x
  out[both] <- log_x[both] + log1p(-exp(log_y[both] - log_x[both]))
  out
}

log_normal_interval_prob <- function(lower, upper) {
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  if (length(lower) != length(upper)) {
    stop("lower and upper must have the same length.")
  }

  out <- rep(-Inf, length(lower))
  full <- is.infinite(lower) & lower < 0 & is.infinite(upper) & upper > 0
  out[full] <- 0
  valid <- is.finite(lower) | is.finite(upper)
  valid <- valid & lower < upper
  if (!any(valid)) return(out)

  right <- valid & lower >= 0
  if (any(right)) {
    out[right] <- logspace_sub(
      pnorm(lower[right], lower.tail = FALSE, log.p = TRUE),
      pnorm(upper[right], lower.tail = FALSE, log.p = TRUE)
    )
  }

  left <- valid & lower < 0
  if (any(left)) {
    out[left] <- logspace_sub(
      pnorm(upper[left], log.p = TRUE),
      pnorm(lower[left], log.p = TRUE)
    )
  }

  out
}

ordinal_loglik_marginal <- function(C, U, A, tau) {
  C <- as.matrix(C)
  U <- as.matrix(U)
  A <- as.matrix(A)

  if (nrow(C) != nrow(U) || ncol(C) != nrow(A) || ncol(U) != ncol(A)) {
    stop("Incompatible dimensions in ordinal_loglik_marginal().")
  }

  ans <- 0
  for (j in seq_len(ncol(C))) {
    mu <- as.numeric(U %*% A[j, ])
    tau_j <- tau[[j]]
    lower <- c(-Inf, tau_j)[C[, j]] - mu
    upper <- c(tau_j, Inf)[C[, j]] - mu
    log_prob <- log_normal_interval_prob(lower, upper)
    if (any(!is.finite(log_prob))) return(-Inf)
    ans <- ans + sum(log_prob)
  }

  ans
}

sample_sigma2_eps_general <- function(y,
                                      X,
                                      U,
                                      logtheta,
                                      kernel = "se",
                                      matern_nu = 2.5) {
  n <- length(y)
  
  cc <- gp_corr_general(X, U, logtheta, kernel, matern_nu)
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
                                          gp_prior,
                                          kernel = "se",
                                          matern_nu = 2.5) {
  n_eval_total <- 0L
  
  for (j in seq_along(logtheta)) {
    logf_j <- function(val) {
      lt <- logtheta
      lt[j] <- val
      theta_logpost_integrated_general(
        y, X, U, lt, gp_prior,
        kernel = kernel,
        matern_nu = matern_nu
      )
    }
    
    ans <- bounded_slice_update(
      x0 = logtheta[j],
      logf = logf_j,
      w = theta_slice_width[j],
      lower = gp_prior$lower[j],
      upper = gp_prior$upper[j],
      max_steps_out = 30L,
      max_iter = 100L,
      fail_on_limit = TRUE
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

############################################################
## Latent U update
############################################################

update_U_ess_block_general <- function(y,
                                       X,
                                       U_curr,
                                       S,
                                       A,
                                       logtheta,
                                       sigma2_eps,
                                       block_idx,
                                       max_try = 200L,
                                       kernel = "se",
                                       matern_nu = 2.5) {
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
    gp_state_general(
      y, X, U_prop, logtheta, sigma2_eps,
      kernel = kernel,
      matern_nu = matern_nu
    )$loglik
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
  
  angle <- runif(1, 0, 2 * base::pi)
  angle_min <- angle - 2 * base::pi
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

update_U_theta_ess_integrated_general <- function(y,
                                                  X,
                                                  C,
                                                  U_curr,
                                                  logtheta_curr,
                                                  gp_prior,
                                                  block_idx,
                                                  reference = c("score", "prior"),
                                                  S = NULL,
                                                  A = NULL,
                                                  tau = NULL,
                                                  update_theta = TRUE,
                                                  max_try = 5000L,
                                                  kernel = "se",
                                                  matern_nu = 2.5) {
  reference <- match.arg(reference)
  block_idx <- as.integer(block_idx)
  d <- ncol(U_curr)

  if (reference == "score") {
    if (is.null(S) || is.null(A)) {
      stop("S and A are required for a score-reference update.")
    }
    ref <- latent_reference_params(S, A)
    M_block <- ref$mean[block_idx, , drop = FALSE]
    V_block <- ref$V
  } else {
    if (is.null(A) || is.null(tau)) {
      stop("A and tau are required for a prior-reference update.")
    }
    M_block <- matrix(0, nrow = length(block_idx), ncol = d)
    V_block <- diag(d)
  }

  U_block <- U_curr[block_idx, , drop = FALSE]
  theta_mean <- if (isTRUE(update_theta)) gp_prior$mean else numeric(0)
  theta_sd <- if (isTRUE(update_theta)) gp_prior$sd else numeric(0)

  residual_log_density <- function(U_block_prop, logtheta_prop) {
    if (any(logtheta_prop < gp_prior$lower) ||
        any(logtheta_prop > gp_prior$upper)) {
      return(-Inf)
    }

    U_prop <- U_curr
    if (length(block_idx) > 0L) {
      U_prop[block_idx, ] <- U_block_prop
    }

    out <- gp_loglik_integrated_general(
      y, X, U_prop, logtheta_prop,
      kernel = kernel,
      matern_nu = matern_nu
    )
    if (!is.finite(out)) return(-Inf)

    if (reference == "prior") {
      out <- out + ordinal_loglik_marginal(C, U_prop, A, tau)
    }

    out
  }

  loglik_cur <- residual_log_density(U_block, logtheta_curr)
  if (!is.finite(loglik_cur)) {
    stop("Current state has a non-finite density in joint ESS update.")
  }

  Nu_U <- if (length(block_idx) > 0L) {
    rmvnorm_rows_common(
      mean_mat = matrix(0, nrow = length(block_idx), ncol = d),
      Sigma = V_block
    )
  } else {
    matrix(numeric(0), nrow = 0L, ncol = d)
  }
  Nu_theta <- if (isTRUE(update_theta)) {
    rnorm(length(logtheta_curr), 0, theta_sd)
  } else {
    numeric(0)
  }

  Z_U <- U_block - M_block
  Z_theta <- if (isTRUE(update_theta)) logtheta_curr - theta_mean else numeric(0)
  logy <- loglik_cur + log(runif(1))

  angle <- runif(1, 0, 2 * base::pi)
  angle_min <- angle - 2 * base::pi
  angle_max <- angle
  n_eval <- 1L

  for (try_id in seq_len(max_try)) {
    U_prop_block <- M_block + Z_U * cos(angle) + Nu_U * sin(angle)
    logtheta_prop <- if (isTRUE(update_theta)) {
      theta_mean + Z_theta * cos(angle) + Nu_theta * sin(angle)
    } else {
      logtheta_curr
    }

    loglik_prop <- residual_log_density(U_prop_block, logtheta_prop)
    n_eval <- n_eval + 1L

    if (is.finite(loglik_prop) && loglik_prop >= logy) {
      U_new <- U_curr
      if (length(block_idx) > 0L) {
        U_new[block_idx, ] <- U_prop_block
      }

      return(list(
        U = U_new,
        logtheta = logtheta_prop,
        n_eval = n_eval,
        accepted = TRUE,
        angle = angle,
        reference = reference
      ))
    }

    if (angle < 0) {
      angle_min <- angle
    } else {
      angle_max <- angle
    }
    angle <- runif(1, angle_min, angle_max)
  }

  stop("Joint elliptical-slice update exceeded max_try.")
}

update_measurement_marginal_slice <- function(C,
                                              U,
                                              A,
                                              tau,
                                              m_vec,
                                              s_A = 2.5,
                                              tau_bound = 6,
                                              ident = c("lower_triangular", "none"),
                                              A_width = 0.5,
                                              tau_width = 0.5) {
  ident <- match.arg(ident)
  q <- nrow(A)
  d <- ncol(A)
  n_eval <- 0L

  active_A <- lapply(seq_len(q), function(j) {
    if (ident == "lower_triangular" && j <= d) seq_len(j) else seq_len(d)
  })

  for (j in seq_len(q)) {
    for (k in active_A[[j]]) {
      lower <- if (ident == "lower_triangular" && j == k && j <= d) 0 else -Inf
      logf <- function(value) {
        A_prop <- A
        A_prop[j, k] <- value
        ordinal_loglik_marginal(C, U, A_prop, tau) -
          0.5 * sum(A_prop^2) / s_A^2
      }
      ans <- bounded_slice_update(
        x0 = A[j, k],
        logf = logf,
        w = A_width,
        lower = lower,
        upper = Inf,
        max_steps_out = 50L,
        max_iter = 500L,
        fail_on_limit = TRUE
      )
      A[j, k] <- ans$x
      n_eval <- n_eval + ans$n_eval
    }
  }

  for (j in seq_len(q)) {
    if (m_vec[j] <= 1L) next
    for (r in seq_len(m_vec[j] - 1L)) {
      lower <- max(-tau_bound, if (r > 1L) tau[[j]][r - 1L] else -Inf)
      upper <- min(tau_bound, if (r < m_vec[j] - 1L) tau[[j]][r + 1L] else Inf)
      logf <- function(value) {
        tau_prop <- tau
        tau_prop[[j]][r] <- value
        ordinal_loglik_marginal(C, U, A, tau_prop)
      }
      ans <- bounded_slice_update(
        x0 = tau[[j]][r],
        logf = logf,
        w = tau_width,
        lower = lower,
        upper = upper,
        max_steps_out = 50L,
        max_iter = 500L,
        fail_on_limit = TRUE
      )
      tau[[j]][r] <- ans$x
      n_eval <- n_eval + ans$n_eval
    }
  }

  list(A = A, tau = tau, n_eval = n_eval)
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

rank_rhat <- function(chain_list) {
  mat <- aligned_chain_matrix(chain_list)
  if (is.null(mat)) return(NA_real_)
  if (requireNamespace("posterior", quietly = TRUE)) {
    return(as.numeric(posterior::rhat(mat)))
  }
  split_rhat(chain_list)
}

bulk_ess <- function(chain_list) {
  mat <- aligned_chain_matrix(chain_list)
  if (is.null(mat)) return(NA_real_)
  if (requireNamespace("posterior", quietly = TRUE)) {
    return(as.numeric(posterior::ess_bulk(mat)))
  }
  sum(vapply(chain_list, ess_ips, numeric(1)))
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

############################################################
## Fully Bayesian ordinal-probit EIV-GP sampler
############################################################

fit_eivgp_ordprobit_fb <- function(X_raw,
                                   y_raw,
                                   C_ord,
                                   U_obs = NULL,
                                   calib_idx = NULL,
                                   U_true_eval = NULL,
                                   d = 2L,
                                   m_vec = NULL,
                                   ident = c("lower_triangular", "none"),
                                   n_iter = 3000L,
                                   burn = 1000L,
                                   thin = 2L,
                                   n_chains = 4L,
                                   preset = "balanced",
                                   sampler_strategy = c("interwoven", "legacy"),
                                   control_overrides = list(),
                                   seed = 1L,
                                   parallel_chains = TRUE,
                                   n_cores = NULL,
                                   verbose = FALSE,
                                   progress_every = 100L,
                                   progress_label = "EIV-GP Study II",
                                   progress_file = NULL,
                                   adaptive_control = NULL,
                                   kernel = c("se", "matern"),
                                   matern_nu = 2.5,
                                   standardize_U = FALSE,
                                   ordinal_levels = NULL,
                                   store_scores = TRUE) {
  fit_call <- match.call()
  ident <- match.arg(ident)
  sampler_strategy <- match.arg(sampler_strategy)
  kernel_spec <- normalize_gp_kernel(kernel, matern_nu)

  if (missing(d) && !is.null(U_obs)) {
    d <- ncol(as.matrix(U_obs))
  }
  d <- as.integer(d)
  if (length(d) != 1L || is.na(d) || d < 1L) {
    stop("d must be one positive integer.")
  }
  standardize_U <- isTRUE(standardize_U)
  store_scores <- isTRUE(store_scores)

  n_iter <- as.integer(n_iter)
  burn <- as.integer(burn)
  thin <- as.integer(thin)
  n_chains <- as.integer(n_chains)
  seed <- as.integer(seed)
  if (anyNA(c(n_iter, burn, thin, n_chains, seed)) ||
      n_iter < 2L || burn < 0L || burn >= n_iter || thin < 1L ||
      n_chains < 1L) {
    stop(
      "Require n_iter >= 2, 0 <= burn < n_iter, thin >= 1, ",
      "n_chains >= 1, and a finite integer seed."
    )
  }
  adaptive_mcmc <- if (is.null(adaptive_control)) {
    NULL
  } else {
    if (!is.list(adaptive_control)) {
      stop("adaptive_control must be NULL or a named list.")
    }
    do.call(mixedgp_adaptive_mcmc_control, adaptive_control)
  }
  max_retained_capacity <- floor((n_iter - burn) / thin)
  if (!is.null(adaptive_mcmc) &&
      adaptive_mcmc$max_draws > max_retained_capacity) {
    stop(
      "adaptive_control$max_draws exceeds the retained-draw capacity implied ",
      "by n_iter, burn, and thin."
    )
  }
  
  set.seed(seed)

  y_matrix <- as_numeric_matrix_strict(y_raw, "y_raw")
  if (ncol(y_matrix) != 1L) stop("y_raw must be univariate.")
  y_raw <- as.numeric(y_matrix[, 1])
  n <- length(y_raw)
  if (n < 2L) stop("At least two training observations are required.")

  X_raw <- as_numeric_matrix_strict(X_raw, "X_raw", nrow_expected = n)
  if (is.null(colnames(X_raw))) colnames(X_raw) <- paste0("x", seq_len(ncol(X_raw)))

  ordinal_data <- prepare_ordinal_matrix(
    C_ord,
    m_vec = m_vec,
    level_maps = ordinal_levels,
    name = "C_ord"
  )
  C_ord <- ordinal_data$C
  m_vec <- ordinal_data$m_vec
  ordinal_levels <- ordinal_data$level_maps
  if (nrow(C_ord) != n) stop("C_ord must have the same number of rows as y_raw.")

  p <- ncol(X_raw)
  q <- ncol(C_ord)
  
  if (ident == "lower_triangular" && q < d) {
    stop("lower_triangular identification requires q >= d.")
  }
  
  infer_calib_idx <- is.null(calib_idx)
  if (infer_calib_idx) calib_idx <- integer(0)
  calib_idx <- sort(as.integer(calib_idx))
  if (anyNA(calib_idx) || any(calib_idx < 1L | calib_idx > n)) {
    stop("calib_idx must contain valid training-row indices.")
  }
  if (anyDuplicated(calib_idx)) {
    stop("calib_idx must not contain duplicates.")
  }
  
  U_obs_full_raw <- matrix(NA_real_, n, d)
  colnames(U_obs_full_raw) <- paste0("u", seq_len(d))
  
  if (!is.null(U_obs)) {
    U_obs <- as_numeric_matrix_strict(
      U_obs,
      "U_obs",
      ncol_expected = d,
      allow_na = TRUE
    )
    if (!is.null(colnames(U_obs))) colnames(U_obs_full_raw) <- colnames(U_obs)
    
    if (nrow(U_obs) == n) {
      if (infer_calib_idx) {
        observed_per_row <- rowSums(!is.na(U_obs))
        if (any(observed_per_row > 0L & observed_per_row < d)) {
          stop(
            "Partially observed rows of U_obs are not yet supported; each row ",
            "must contain either all d coordinates or none."
          )
        }
        calib_idx <- which(stats::complete.cases(U_obs))
      }

      if (length(calib_idx) > 0L && anyNA(U_obs[calib_idx, , drop = FALSE])) {
        stop("Every row selected by calib_idx must contain all d latent coordinates.")
      }
      U_obs_full_raw[calib_idx, ] <- U_obs[calib_idx, , drop = FALSE]
    } else if (nrow(U_obs) == length(calib_idx)) {
      if (anyNA(U_obs)) stop("Calibration-only U_obs must be complete.")
      U_obs_full_raw[calib_idx, ] <- U_obs
    } else {
      stop("U_obs must have either n rows or length(calib_idx) rows.")
    }
  } else {
    if (length(calib_idx) > 0L) {
      stop("calib_idx was supplied but U_obs is NULL.")
    }
  }

  anchor_status <- mixedgp_latent_anchor_status(
    U_obs_full_raw, calib_idx, d = d
  )
  if (length(calib_idx) > 0L && !isTRUE(anchor_status$anchored)) {
    warning(
      "The calibration data inform the latent state but do not have affine ",
      "rank d + 1; raw-scale latent outputs remain unavailable.",
      call. = FALSE
    )
  }

  if (ident == "none" && length(calib_idx) == 0L && d > 1L) {
    warning(
      "With multivariate completely latent U, ident = 'none' leaves rotational ",
      "nonidentifiability. Use a loading constraint unless only predictions ",
      "in the unidentified coordinate system are needed."
    )
  }

  U_center <- rep(0, d)
  U_scale <- rep(1, d)
  if (standardize_U && length(calib_idx) > 0L) {
    U_center <- colMeans(U_obs_full_raw[calib_idx, , drop = FALSE])
    U_scale <- apply(U_obs_full_raw[calib_idx, , drop = FALSE], 2, sd)
    U_scale[!is.finite(U_scale) | U_scale <= 0] <- 1
  }

  U_obs_full <- U_obs_full_raw
  if (length(calib_idx) > 0L) {
    U_obs_full[calib_idx, ] <- sweep(
      sweep(U_obs_full_raw[calib_idx, , drop = FALSE], 2, U_center, "-"),
      2,
      U_scale,
      "/"
    )
  }

  U_true_eval_raw <- U_true_eval
  if (!is.null(U_true_eval)) {
    U_true_eval_raw <- as_numeric_matrix_strict(
      U_true_eval,
      "U_true_eval",
      nrow_expected = n,
      ncol_expected = d
    )
    U_true_eval <- sweep(
      sweep(U_true_eval_raw, 2, U_center, "-"),
      2,
      U_scale,
      "/"
    )
  }
  
  miss_idx <- setdiff(seq_len(n), calib_idx)
  
  X_center <- colMeans(X_raw)
  X_scale <- apply(X_raw, 2, sd)
  X_scale[!is.finite(X_scale) | X_scale <= 0] <- 1
  
  X <- sweep(sweep(X_raw, 2, X_center, "-"), 2, X_scale, "/")
  
  y_center <- mean(y_raw)
  y_scale <- sd(y_raw)
  if (!is.finite(y_scale) || y_scale <= 0) {
    stop("y_raw must have positive finite sample variance.")
  }
  y <- as.numeric((y_raw - y_center) / y_scale)
  
  control <- make_default_control_ordprobit(
    n = n,
    n_mis = length(miss_idx),
    preset = preset,
    d = d
  )
  if (length(control_overrides) > 0L) {
    unknown_control <- setdiff(names(control_overrides), names(control))
    if (length(unknown_control) > 0L) {
      stop(
        "Unknown control override(s): ",
        paste(unknown_control, collapse = ", ")
      )
    }
    control[names(control_overrides)] <- control_overrides
  }
  
  gp_prior <- make_gp_prior(p = p, d = d)
  
  n_logtheta <- 1L + p + d
  control$theta_slice_width_init <- rep(0.8, n_logtheta)
  
  n_save <- floor((n_iter - burn) / thin)
  
  if (n_save <= 0L) {
    stop("n_iter, burn, and thin imply no saved draws.")
  }
  
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
    assert_ordprobit_state(
      C = C_ord, U = U0, S = S0, A = A0, tau = tau0,
      m_vec = m_vec, ident = ident, calib_idx = calib_idx,
      U_obs = U_obs_full, logtheta = logtheta0,
      sigma2_eps = sigma2_eps0
    )
    
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
    chain_start_time <- proc.time()[["elapsed"]]
    format_duration <- function(seconds) {
      seconds <- max(0, as.numeric(seconds))
      sprintf(
        "%02d:%02d:%02d",
        floor(seconds / 3600),
        floor((seconds %% 3600) / 60),
        floor(seconds %% 60)
      )
    }
    set.seed(chain_seed)
    
    state <- initialize_chain_state(chain_seed)
    initial_state <- state
    
    U_curr <- state$U_curr
    S_curr <- state$S_curr
    A_curr <- state$A_curr
    tau_curr <- state$tau_curr
    sigma2_eps <- state$sigma2_eps
    logtheta <- state$logtheta
    theta_slice_width <- state$theta_slice_width
    
    samples_U_chain <- array(NA_real_, dim = c(n_save, n, d))
    samples_S_chain <- if (store_scores) {
      array(NA_real_, dim = c(n_save, n, q))
    } else {
      NULL
    }
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
    joint_score_eval_total <- 0L
    joint_score_update_total <- 0L
    joint_global_eval_total <- 0L
    joint_global_update_total <- 0L
    joint_theta_eval_total <- 0L
    joint_theta_update_total <- 0L
    joint_collapsed_eval_total <- 0L
    joint_collapsed_update_total <- 0L
    marginal_measurement_eval_total <- 0L
    marginal_measurement_update_total <- 0L
    
    save_id <- 0L
    adaptive_schedule <- list()
    adaptive_extensions <- 0L
    adaptive_termination <- if (is.null(adaptive_mcmc)) {
      "fixed_iterations"
    } else {
      NA_character_
    }
    planned_draws <- if (is.null(adaptive_mcmc)) {
      n_save
    } else {
      adaptive_mcmc$initial_draws
    }
    iter_limit <- min(n_iter, burn + planned_draws * thin)
    iter <- 0L
    
    while (iter < iter_limit) {
      iter <- iter + 1L
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
      
      if (sampler_strategy == "legacy") {
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
              max_try = control$max_ess_try,
              kernel = kernel_spec$name,
              matern_nu = kernel_spec$matern_nu
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
              max_try = control$max_ess_try,
              kernel = kernel_spec$name,
              matern_nu = kernel_spec$matern_nu
            )
            
            U_curr <- gu$U
            
            global_u_eval_total <- global_u_eval_total + gu$n_eval
            global_u_accept_total <- global_u_accept_total + as.integer(gu$accepted)
            global_u_total <- global_u_total + 1L
          }
        }
        
        if (iter %% control$theta_update_every == 0L) {
          th <- update_logtheta_slice_general(
            y = y,
            X = X,
            U = U_curr,
            logtheta = logtheta,
            theta_slice_width = theta_slice_width,
            gp_prior = gp_prior,
            kernel = kernel_spec$name,
            matern_nu = kernel_spec$matern_nu
          )
          
          logtheta <- th$logtheta
          
          theta_eval_total <- theta_eval_total + th$n_eval
          theta_update_total <- theta_update_total + 1L
        }
      } else {
        if (length(miss_idx) > 0L) {
          for (bb in seq_len(control$joint_local_blocks_per_iter)) {
            block_idx <- sample(
              miss_idx,
              min(control$u_block_size, length(miss_idx))
            )
            ju <- update_U_theta_ess_integrated_general(
              y = y,
              X = X,
              C = C_ord,
              U_curr = U_curr,
              logtheta_curr = logtheta,
              gp_prior = gp_prior,
              block_idx = block_idx,
              reference = "score",
              S = S_curr,
              A = A_curr,
              update_theta = FALSE,
              max_try = control$max_ess_try,
              kernel = kernel_spec$name,
              matern_nu = kernel_spec$matern_nu
            )
            U_curr <- ju$U
            joint_score_eval_total <- joint_score_eval_total + ju$n_eval
            joint_score_update_total <- joint_score_update_total + 1L
          }
        }

        if (control$joint_theta_every > 0L &&
            iter %% control$joint_theta_every == 0L) {
          jt <- update_U_theta_ess_integrated_general(
            y = y,
            X = X,
            C = C_ord,
            U_curr = U_curr,
            logtheta_curr = logtheta,
            gp_prior = gp_prior,
            block_idx = integer(0),
            reference = "score",
            S = S_curr,
            A = A_curr,
            update_theta = TRUE,
            max_try = control$max_ess_try,
            kernel = kernel_spec$name,
            matern_nu = kernel_spec$matern_nu
          )
          logtheta <- jt$logtheta
          joint_theta_eval_total <- joint_theta_eval_total + jt$n_eval
          joint_theta_update_total <- joint_theta_update_total + 1L
        }

        if (control$joint_global_every > 0L &&
            iter %% control$joint_global_every == 0L) {
          jg <- update_U_theta_ess_integrated_general(
            y = y,
            X = X,
            C = C_ord,
            U_curr = U_curr,
            logtheta_curr = logtheta,
            gp_prior = gp_prior,
            block_idx = miss_idx,
            reference = "score",
            S = S_curr,
            A = A_curr,
            update_theta = TRUE,
            max_try = control$max_ess_try,
            kernel = kernel_spec$name,
            matern_nu = kernel_spec$matern_nu
          )
          U_curr <- jg$U
          logtheta <- jg$logtheta
          joint_global_eval_total <- joint_global_eval_total + jg$n_eval
          joint_global_update_total <- joint_global_update_total + 1L
        }

        if (control$joint_collapsed_every > 0L &&
            iter %% control$joint_collapsed_every == 0L) {
          jc <- update_U_theta_ess_integrated_general(
            y = y,
            X = X,
            C = C_ord,
            U_curr = U_curr,
            logtheta_curr = logtheta,
            gp_prior = gp_prior,
            block_idx = miss_idx,
            reference = "prior",
            A = A_curr,
            tau = tau_curr,
            max_try = control$max_ess_try,
            kernel = kernel_spec$name,
            matern_nu = kernel_spec$matern_nu
          )
          U_curr <- jc$U
          logtheta <- jc$logtheta
          joint_collapsed_eval_total <- joint_collapsed_eval_total + jc$n_eval
          joint_collapsed_update_total <- joint_collapsed_update_total + 1L
          S_curr <- sample_scores_ord(C_ord, U_curr, A_curr, tau_curr)
        }

        if (control$marginal_measurement_every > 0L &&
            iter %% control$marginal_measurement_every == 0L) {
          mm <- update_measurement_marginal_slice(
            C = C_ord,
            U = U_curr,
            A = A_curr,
            tau = tau_curr,
            m_vec = m_vec,
            s_A = control$s_A,
            tau_bound = control$tau_bound,
            ident = ident
          )
          A_curr <- mm$A
          tau_curr <- mm$tau
          marginal_measurement_eval_total <-
            marginal_measurement_eval_total + mm$n_eval
          marginal_measurement_update_total <-
            marginal_measurement_update_total + 1L
          S_curr <- sample_scores_ord(C_ord, U_curr, A_curr, tau_curr)
        }

      }
      
      if (length(calib_idx) > 0L) {
        U_curr[calib_idx, ] <- U_obs_full[calib_idx, ]
      }
      
      sigma2_eps <- sample_sigma2_eps_general(
        y, X, U_curr, logtheta,
        kernel = kernel_spec$name,
        matern_nu = kernel_spec$matern_nu
      )
      
      logtheta_trace_all[iter, ] <- logtheta
      
      if (sampler_strategy == "legacy" && control$adapt_theta_width &&
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
        assert_ordprobit_state(
          C = C_ord, U = U_curr, S = S_curr, A = A_curr,
          tau = tau_curr, m_vec = m_vec, ident = ident,
          calib_idx = calib_idx, U_obs = U_obs_full,
          logtheta = logtheta, sigma2_eps = sigma2_eps
        )
        
        samples_U_chain[save_id, , ] <- U_curr
        if (store_scores) samples_S_chain[save_id, , ] <- S_curr
        samples_A_chain[save_id, , ] <- A_curr
        samples_tau_chain[save_id, ] <- flatten_tau(tau_curr)
        samples_logtheta_chain[save_id, ] <- logtheta
        samples_sigma2_chain[save_id] <- sigma2_eps
      }

      if (!is.null(adaptive_mcmc) && save_id == planned_draws) {
        a_trace <- matrix(
          samples_A_chain[seq_len(save_id), , , drop = FALSE],
          nrow = save_id
        )
        key_trace <- cbind(
          sigma_epsilon = sqrt(samples_sigma2_chain[seq_len(save_id)]),
          exp(samples_logtheta_chain[seq_len(save_id), , drop = FALSE]),
          a_trace,
          samples_tau_chain[seq_len(save_id), , drop = FALSE]
        )
        key_names <- c(
          "sigma_epsilon", "rho",
          paste0("theta_x", seq_len(p)),
          paste0("theta_u", seq_len(d)),
          paste0(
            "A[", rep(seq_len(q), times = d), ",",
            rep(seq_len(d), each = q), "]"
          ),
          tau_flat_names
        )
        colnames(key_trace) <- key_names
        key_ess <- apply(key_trace, 2L, ess_ips)
        finite_key <- is.finite(key_ess)
        min_chain_ess <- if (any(finite_key)) {
          min(key_ess[finite_key])
        } else {
          NA_real_
        }
        limiting_parameter <- if (any(finite_key)) {
          names(which.min(replace(key_ess, !finite_key, Inf)))[1L]
        } else {
          NA_character_
        }
        next_step <- mixedgp_adaptive_next_draws(
          current_draws = save_id,
          observed_ess = min_chain_ess,
          control = adaptive_mcmc,
          n_chains = n_chains
        )
        target_met <- is.finite(min_chain_ess) &&
          min_chain_ess >= next_step$target_per_chain
        at_cap <- save_id >= adaptive_mcmc$max_draws
        termination <- if (target_met) {
          "ess_target_met"
        } else if (at_cap || next_step$next_draws <= save_id) {
          "max_draws_reached"
        } else {
          "continued"
        }
        adaptive_schedule[[length(adaptive_schedule) + 1L]] <- data.frame(
          chain = chain_id,
          check = length(adaptive_schedule) + 1L,
          retained_draws = save_id,
          min_chain_ess = min_chain_ess,
          target_ess_per_chain = next_step$target_per_chain,
          limiting_parameter = limiting_parameter,
          next_retained_draws = next_step$next_draws,
          elapsed_seconds = proc.time()[["elapsed"]] - chain_start_time,
          decision = termination,
          stringsAsFactors = FALSE
        )
        if (termination == "continued") {
          adaptive_extensions <- adaptive_extensions + 1L
          planned_draws <- next_step$next_draws
          iter_limit <- min(n_iter, burn + planned_draws * thin)
        } else {
          adaptive_termination <- termination
          iter_limit <- iter
        }
      }

      if ((verbose || !is.null(progress_file)) &&
          progress_every > 0L &&
          (iter == 1L || iter %% progress_every == 0L || iter == iter_limit)) {
        elapsed <- proc.time()[["elapsed"]] - chain_start_time
        remaining_this_chain <- elapsed / iter * (iter_limit - iter)
        remaining_later_chains <- if (use_mclapply) {
          0
        } else {
          elapsed / iter * iter_limit * max(0L, n_chains - chain_id)
        }
        pct <- 100 * iter / iter_limit
        bar_width <- 24L
        n_done <- min(bar_width, floor(bar_width * iter / iter_limit))
        progress_bar <- paste0(
          "[",
          paste(rep("=", n_done), collapse = ""),
          paste(rep(" ", bar_width - n_done), collapse = ""),
          "]"
        )
        progress_message <- sprintf(
          "%s %s chain %d/%d %s %5.1f%% iter %d/%d elapsed %s ETA %s",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          progress_label,
          chain_id,
          n_chains,
          progress_bar,
          pct,
          iter,
          iter_limit,
          format_duration(elapsed),
          format_duration(remaining_this_chain + remaining_later_chains)
        )
        show_console_progress <- verbose && (!use_mclapply || chain_id == 1L)
        if (show_console_progress) cat("\r", progress_message, sep = "")
        if (!is.null(progress_file)) {
          cat(progress_message, "\n", file = progress_file, append = TRUE)
        }
        if (show_console_progress && iter == iter_limit) cat("\n")
        if (show_console_progress) flush.console()
      }
    }
    
    if (save_id < n_save) {
      samples_U_chain <- samples_U_chain[seq_len(save_id), , , drop = FALSE]
      if (store_scores) {
        samples_S_chain <- samples_S_chain[seq_len(save_id), , , drop = FALSE]
      }
      samples_A_chain <- samples_A_chain[seq_len(save_id), , , drop = FALSE]
      samples_tau_chain <- samples_tau_chain[seq_len(save_id), , drop = FALSE]
      samples_logtheta_chain <- samples_logtheta_chain[seq_len(save_id), , drop = FALSE]
      samples_sigma2_chain <- samples_sigma2_chain[seq_len(save_id)]
    }
    
    if (is.null(adaptive_mcmc)) {
      adaptive_schedule_df <- data.frame()
    } else {
      adaptive_schedule_df <- do.call(rbind, adaptive_schedule)
      if (is.na(adaptive_termination)) {
        adaptive_termination <- if (save_id >= adaptive_mcmc$max_draws) {
          "max_draws_reached"
        } else {
          "iteration_capacity_reached"
        }
      }
    }

    stats <- data.frame(
      chain = chain_id,
      seed = chain_seed,
      iterations_run = iter,
      saved = save_id,
      adaptive_extensions = adaptive_extensions,
      adaptive_termination = adaptive_termination,
      u_ess_eval_total = u_ess_eval_total,
      u_ess_accept_total = u_ess_accept_total,
      u_ess_total = u_ess_total,
      global_u_eval_total = global_u_eval_total,
      global_u_accept_total = global_u_accept_total,
      global_u_total = global_u_total,
      theta_eval_total = theta_eval_total,
      theta_update_total = theta_update_total,
      joint_score_eval_total = joint_score_eval_total,
      joint_score_update_total = joint_score_update_total,
      joint_global_eval_total = joint_global_eval_total,
      joint_global_update_total = joint_global_update_total,
      joint_theta_eval_total = joint_theta_eval_total,
      joint_theta_update_total = joint_theta_update_total,
      joint_collapsed_eval_total = joint_collapsed_eval_total,
      joint_collapsed_update_total = joint_collapsed_update_total,
      marginal_measurement_eval_total = marginal_measurement_eval_total,
      marginal_measurement_update_total = marginal_measurement_update_total
    )
    
    list(
      chain_id = chain_id,
      seed = chain_seed,
      sampler_strategy = sampler_strategy,
      samples_U = samples_U_chain,
      samples_S = samples_S_chain,
      samples_A = samples_A_chain,
      samples_tau = samples_tau_chain,
      samples_logtheta = samples_logtheta_chain,
      samples_sigma2 = samples_sigma2_chain,
      initial_state = initial_state,
      final_state = list(
        U_curr = U_curr,
        S_curr = S_curr,
        A_curr = A_curr,
        tau_curr = tau_curr,
        sigma2_eps = sigma2_eps,
        logtheta = logtheta,
        theta_slice_width = theta_slice_width,
        ## Publication simulations source this module into an isolated
        ## environment whose parent intentionally excludes .GlobalEnv.
        ## The RNG state itself nevertheless belongs to .GlobalEnv.
        rng_state = mixedgp_rng_state()$value
      ),
      adaptive_schedule = adaptive_schedule_df,
      theta_slice_width_final = theta_slice_width,
      stats = stats
    )
  }
  
  chain_seeds <- seed + 10000L * seq_len(n_chains)
  
  parallel_backend <- mixedgp_parallel_backend(
    if (isTRUE(parallel_chains)) n_cores else 1L
  )
  use_mclapply <- parallel_backend$backend == "fork" && n_chains > 1L
  mc_cores <- min(n_chains, parallel_backend$cores)
  
  if (verbose) {
    cat("Running fully Bayesian ordinal-probit EIV-GP with", n_chains, "chain(s).\n")
  }
  
  mcmc_time <- system.time({
    chains <- mixedgp_parallel_lapply(
        as.list(seq_len(n_chains)),
        function(cc) {
          run_one_chain(chain_id = cc, chain_seed = chain_seeds[cc])
        },
        n_cores = if (use_mclapply) mc_cores else 1L,
        seeds = chain_seeds,
        mc.preschedule = FALSE
    )
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
  samples_S <- if (store_scores) {
    combine_chain_arrays(samples_by_chain$S)
  } else {
    NULL
  }
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
  adaptive_schedule <- if (is.null(adaptive_mcmc)) {
    data.frame()
  } else {
    do.call(rbind, lapply(chains, function(z) z$adaptive_schedule))
  }
  
  hyper_names <- c(
    "sigma_epsilon",
    "rho",
    paste0("theta_x", seq_len(p)),
    paste0("theta_u", seq_len(d))
  )
  
  rhat_hyper <- data.frame(
    parameter = hyper_names,
    rhat = c(
      rank_rhat(lapply(samples_by_chain$sigma2, function(v) sqrt(v))),
      sapply(seq_len(n_logtheta), function(k) {
        rank_rhat(lapply(samples_by_chain$logtheta, function(mat) exp(mat[, k])))
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
            rhat = rank_rhat(lapply(samples_by_chain$A, function(arr) arr[, j, k]))
          )
        })
      )
    })
  )
  
  rhat_tau <- data.frame(
    parameter = colnames(samples_tau),
    rhat = sapply(seq_len(ncol(samples_tau)), function(k) {
      rank_rhat(lapply(samples_by_chain$tau, function(mat) mat[, k]))
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
              rhat = rank_rhat(lapply(samples_by_chain$U, function(arr) arr[, ii, k]))
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
  
  key_chains <- list(
    sigma_epsilon = lapply(samples_by_chain$sigma2, sqrt),
    rho = lapply(samples_by_chain$logtheta, function(mat) exp(mat[, 1]))
  )

  for (j in seq_len(p)) {
    key_chains[[paste0("theta_x", j)]] <- lapply(
      samples_by_chain$logtheta,
      function(mat) exp(mat[, 1L + j])
    )
  }
  for (k in seq_len(d)) {
    key_chains[[paste0("theta_u", k)]] <- lapply(
      samples_by_chain$logtheta,
      function(mat) exp(mat[, 1L + p + k])
    )
  }
  for (j in seq_len(q)) {
    for (k in seq_len(d)) {
      fixed_zero <- ident == "lower_triangular" && j <= d && k > j
      if (!fixed_zero) {
        key_chains[[paste0("A", j, k)]] <- lapply(
          samples_by_chain$A,
          function(arr) arr[, j, k]
        )
      }
    }
  }
  for (k in seq_len(tau_total)) {
    key_chains[[tau_flat_names[k]]] <- lapply(
      samples_by_chain$tau,
      function(mat) mat[, k]
    )
  }

  ess_key <- data.frame(
    parameter = names(key_chains),
    rhat = vapply(key_chains, rank_rhat, numeric(1)),
    ess_bulk = vapply(key_chains, bulk_ess, numeric(1)),
    ess_tail = vapply(key_chains, tail_ess, numeric(1)),
    stringsAsFactors = FALSE
  )
  ess_key$ess <- pmin(ess_key$ess_bulk, ess_key$ess_tail, na.rm = TRUE)
  
  diagnostics_summary <- data.frame(
    sampler_strategy = sampler_strategy,
    kernel = kernel_spec$name,
    matern_nu = kernel_spec$matern_nu,
    n_chains = n_chains,
    parallel_backend = if (use_mclapply) "fork" else "serial",
    parallel_cores = if (use_mclapply) mc_cores else 1L,
    n_iter = max(chain_stats$iterations_run),
    n_iter_capacity = n_iter,
    burn = burn,
    thin = thin,
    adaptive = !is.null(adaptive_mcmc),
    adaptive_initial_draws = if (is.null(adaptive_mcmc)) NA_integer_ else adaptive_mcmc$initial_draws,
    adaptive_max_draws = if (is.null(adaptive_mcmc)) NA_integer_ else adaptive_mcmc$max_draws,
    adaptive_target_ess = if (is.null(adaptive_mcmc)) NA_real_ else adaptive_mcmc$target_ess,
    adaptive_extensions = sum(chain_stats$adaptive_extensions),
    adaptive_termination = paste(unique(chain_stats$adaptive_termination), collapse = ";"),
    min_saved_per_chain = min(chain_stats$saved),
    max_saved_per_chain = max(chain_stats$saved),
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
    total_gp_evaluations = sum(
      chain_stats$u_ess_eval_total,
      chain_stats$global_u_eval_total,
      chain_stats$theta_eval_total,
      chain_stats$joint_score_eval_total,
      chain_stats$joint_global_eval_total,
      chain_stats$joint_theta_eval_total,
      chain_stats$joint_collapsed_eval_total
    ),
    covariance_jitter = 0,
    forms_explicit_covariance_inverse = FALSE,
    time_seconds = as.numeric(mcmc_time["elapsed"])
  )
  
  list(
    call = fit_call,
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
      C_level_maps = ordinal_levels,
      C_names = colnames(C_ord),
      U_obs = U_obs_full,
      U_obs_raw = U_obs_full_raw,
      U_names = colnames(U_obs_full_raw),
      U_center = U_center,
      U_scale = U_scale,
      standardize_U = standardize_U,
      U_true_eval = U_true_eval,
      U_true_eval_raw = U_true_eval_raw,
      calib_idx = calib_idx,
      miss_idx = miss_idx,
      latent_scale_anchored = isTRUE(anchor_status$anchored),
      latent_anchor_rank = anchor_status$affine_rank,
      latent_anchor_required_rank = anchor_status$required_rank,
      m_vec = m_vec,
      q = q,
      d = d,
      p = p,
      ident = ident
    ),
    gp_prior = gp_prior,
    kernel = kernel_spec,
    control = control,
    sampler_strategy = sampler_strategy,
    mcmc = list(
      samples_U = samples_U,
      samples_S = samples_S,
      samples_A = samples_A,
      samples_tau = samples_tau,
      samples_logtheta = samples_logtheta,
      samples_sigma2 = samples_sigma2,
      samples_by_chain = samples_by_chain,
      chain_initial = lapply(chains, function(z) z$initial_state),
      chain_final = lapply(chains, function(z) z$final_state),
      mcmc_draw_info = mcmc_draw_info,
      chain_stats = chain_stats,
      adaptive_control = adaptive_mcmc,
      adaptive_schedule = adaptive_schedule
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
## Prediction for general ordinal-probit EIV-GP
############################################################

sample_u_given_c_ordprobit <- function(C_new,
                                       A,
                                       tau,
                                       m_vec,
                                       n_gibbs = 100L,
                                       U_init = NULL) {
  C_new <- as.matrix(C_new)
  A <- as.matrix(A)
  
  N <- nrow(C_new)
  q <- ncol(C_new)
  d <- ncol(A)
  
  if (nrow(A) != q || length(m_vec) != q || length(tau) != q) {
    stop("m_vec has wrong length.")
  }
  if (any(C_new < 1L) ||
      any(C_new > matrix(rep(m_vec, each = N), N, q))) {
    stop("C_new contains invalid ordinal codes.")
  }
  n_gibbs <- as.integer(n_gibbs)
  if (length(n_gibbs) != 1L || is.na(n_gibbs) || n_gibbs < 1L) {
    stop("n_gibbs must be a positive integer.")
  }
  
  if (is.null(U_init)) {
    U <- matrix(rnorm(N * d), N, d)
  } else {
    U <- as_numeric_matrix_strict(
      U_init,
      "U_init",
      nrow_expected = N,
      ncol_expected = d
    )
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

sample_u_given_c_ordprobit_rejection <- function(
    C_new,
    A,
    tau,
    m_vec,
    batch_size = NULL,
    max_batches = 1000L) {
  C_new <- as.matrix(C_new)
  storage.mode(C_new) <- "integer"
  A <- as.matrix(A)

  N <- nrow(C_new)
  q <- ncol(C_new)
  d <- ncol(A)
  if (N < 1L || nrow(A) != q || length(m_vec) != q || length(tau) != q) {
    stop("Incompatible inputs in sample_u_given_c_ordprobit_rejection().")
  }
  if (anyNA(C_new) || any(C_new < 1L) ||
      any(C_new > matrix(rep(m_vec, each = N), N, q))) {
    stop("C_new contains invalid ordinal codes.")
  }
  if (is.null(batch_size)) batch_size <- max(5000L, 20L * N)
  batch_size <- as.integer(batch_size)
  max_batches <- as.integer(max_batches)
  if (is.na(batch_size) || batch_size < 1L ||
      is.na(max_batches) || max_batches < 1L) {
    stop("batch_size and max_batches must be positive integers.")
  }

  target_keys <- pattern_key(C_new)
  unique_keys <- unique(target_keys)
  required <- table(factor(target_keys, levels = unique_keys))
  accepted <- setNames(vector("list", length(unique_keys)), unique_keys)
  accepted_n <- setNames(integer(length(unique_keys)), unique_keys)
  proposal_hits <- setNames(numeric(length(unique_keys)), unique_keys)
  total_candidates <- 0
  batches_used <- 0L

  for (batch_id in seq_len(max_batches)) {
    batches_used <- batch_id
    U_batch <- matrix(rnorm(batch_size * d), batch_size, d)
    S_batch <- U_batch %*% t(A) + matrix(rnorm(batch_size * q), batch_size, q)
    C_batch <- matrix(NA_integer_, batch_size, q)
    for (j in seq_len(q)) {
      C_batch[, j] <- findInterval(
        S_batch[, j], c(-Inf, tau[[j]], Inf), rightmost.closed = TRUE
      )
    }
    candidate_keys <- pattern_key(C_batch)
    total_candidates <- total_candidates + batch_size

    for (key in unique_keys) {
      hit <- which(candidate_keys == key)
      proposal_hits[[key]] <- proposal_hits[[key]] + length(hit)
      remaining <- required[[key]] - accepted_n[[key]]
      if (remaining <= 0L) next
      if (length(hit) == 0L) next
      take <- hit[seq_len(min(remaining, length(hit)))]
      accepted[[key]][[length(accepted[[key]]) + 1L]] <-
        U_batch[take, , drop = FALSE]
      accepted_n[[key]] <- accepted_n[[key]] + length(take)
    }

    if (all(accepted_n >= required)) break
  }

  telemetry <- data.frame(
    pattern = unique_keys,
    requested_draws = as.integer(required),
    retained_draws = as.integer(accepted_n),
    proposal_hits = as.numeric(proposal_hits),
    shortfall = pmax(as.integer(required) - as.integer(accepted_n), 0L),
    complete = as.integer(accepted_n) >= as.integer(required),
    total_candidates = rep.int(total_candidates, length(unique_keys)),
    batches_used = rep.int(batches_used, length(unique_keys)),
    empirical_pattern_probability =
      as.numeric(proposal_hits) / total_candidates,
    retained_fraction_of_candidates =
      as.numeric(accepted_n) / total_candidates,
    stringsAsFactors = FALSE
  )

  if (any(accepted_n < required)) {
    missing <- unique_keys[accepted_n < required]
    message <- paste0(
      "Exact rejection sampling did not fill ordinal pattern(s): ",
      paste(missing, collapse = ", "), ". Use the exact default ",
      "latent_sampler = 'minimax_tilting', or increase ",
      "rejection_max_batches. Finite Gibbs is diagnostic only."
    )
    stop(structure(
      list(
        message = message,
        call = NULL,
        pattern_telemetry = telemetry,
        candidate_draws = total_candidates,
        retained_draws = sum(accepted_n),
        overall_acceptance = sum(accepted_n) / total_candidates
      ),
      class = c("mixedgp_rejection_sampling_error", "error", "condition")
    ))
  }

  out <- matrix(NA_real_, N, d)
  for (key in unique_keys) {
    rows <- which(target_keys == key)
    draws <- do.call(rbind, accepted[[key]])
    out[rows, ] <- draws[seq_along(rows), , drop = FALSE]
  }
  attr(out, "latent_sampler") <- "exact_rejection"
  attr(out, "candidate_draws") <- total_candidates
  attr(out, "overall_acceptance") <- N / total_candidates
  attr(out, "pattern_telemetry") <- telemetry
  attr(out, "rejection_telemetry") <- telemetry
  out
}

sample_u_given_c_ordprobit_minimax <- function(C_new,
                                                A,
                                                tau,
                                                m_vec) {
  if (!requireNamespace("TruncatedNormal", quietly = TRUE)) {
    stop(
      "latent_sampler='minimax_tilting' requires the CRAN package ",
      "TruncatedNormal.",
      call. = FALSE
    )
  }
  C_new <- as.matrix(C_new)
  storage.mode(C_new) <- "integer"
  A <- as.matrix(A)
  N <- nrow(C_new)
  q <- ncol(C_new)
  d <- ncol(A)
  m_vec <- as.integer(m_vec)
  if (N < 1L || nrow(A) != q || length(m_vec) != q || length(tau) != q ||
      anyNA(C_new) || any(!is.finite(A)) ||
      any(C_new < 1L) ||
      any(C_new > matrix(rep(m_vec, each = N), N, q))) {
    stop("Incompatible inputs in sample_u_given_c_ordprobit_minimax().")
  }
  for (j in seq_len(q)) {
    if (length(tau[[j]]) != m_vec[j] - 1L ||
        any(!is.finite(tau[[j]])) ||
        is.unsorted(tau[[j]], strictly = TRUE)) {
      stop("tau contains an invalid cutpoint vector.")
    }
  }

  ## If S = A U + epsilon with U ~ N(0,I_d) and epsilon ~ N(0,I_q),
  ## then S ~ N(0, A A' + I_q). Use the minimax-tilted accept--reject sampler
  ## for S in its ordinal rectangle; on successful completion this draw and the
  ## following Gaussian U | S update are exact.
  Sigma_S <- tcrossprod(A) + diag(q)
  precision_U_given_S <- diag(d) + crossprod(A)
  chol_precision <- safe_chol(precision_U_given_S)
  V <- solve_chol(chol_precision, diag(d))
  B <- V %*% t(A)
  V <- 0.5 * (V + t(V))
  chol_V <- safe_chol(V)

  target_keys <- pattern_key(C_new)
  unique_keys <- unique(target_keys)
  out <- matrix(NA_real_, N, d)
  telemetry <- vector("list", length(unique_keys))
  for (ii in seq_along(unique_keys)) {
    key <- unique_keys[ii]
    rows <- which(target_keys == key)
    c_key <- C_new[rows[1L], ]
    lower <- upper <- numeric(q)
    for (j in seq_len(q)) {
      lower[j] <- c(-Inf, tau[[j]])[c_key[j]]
      upper[j] <- c(tau[[j]], Inf)[c_key[j]]
    }
    S_draw <- TruncatedNormal::rtmvnorm(
      n = length(rows),
      mu = rep(0, q),
      sigma = Sigma_S,
      lb = lower,
      ub = upper
    )
    if (length(S_draw) != length(rows) * q || any(!is.finite(S_draw))) {
      stop("Minimax tilting returned malformed or non-finite score draws.")
    }
    S_draw <- matrix(S_draw, nrow = length(rows), ncol = q)
    lower_mat <- matrix(lower, nrow = length(rows), ncol = q, byrow = TRUE)
    upper_mat <- matrix(upper, nrow = length(rows), ncol = q, byrow = TRUE)
    if (any(S_draw < lower_mat | S_draw > upper_mat)) {
      stop("Minimax tilting returned a score outside its ordinal rectangle.")
    }
    U_mean <- S_draw %*% t(B)
    U_draw <- U_mean + matrix(rnorm(length(rows) * d), length(rows), d) %*%
      chol_V
    if (any(!is.finite(U_draw))) {
      stop("The exact Gaussian U | S update returned a non-finite draw.")
    }
    out[rows, ] <- U_draw
    telemetry[[ii]] <- data.frame(
      pattern = key,
      n_draw = length(rows),
      stringsAsFactors = FALSE
    )
  }
  attr(out, "latent_sampler") <- "exact_minimax_tilting"
  attr(out, "pattern_telemetry") <- do.call(rbind, telemetry)
  out
}

sample_u_given_c_ordprobit_dispatch <- function(
    C_new,
    A,
    tau,
    m_vec,
    latent_sampler = c("minimax_tilting", "rejection", "gibbs"),
    n_gibbs = 100L,
    rejection_batch_size = NULL,
    rejection_max_batches = 1000L) {
  latent_sampler <- match.arg(latent_sampler)
  if (latent_sampler == "minimax_tilting") {
    return(sample_u_given_c_ordprobit_minimax(C_new, A, tau, m_vec))
  }
  if (latent_sampler == "rejection") {
    return(sample_u_given_c_ordprobit_rejection(
      C_new = C_new,
      A = A,
      tau = tau,
      m_vec = m_vec,
      batch_size = rejection_batch_size,
      max_batches = rejection_max_batches
    ))
  }
  out <- sample_u_given_c_ordprobit(
    C_new = C_new,
    A = A,
    tau = tau,
    m_vec = m_vec,
    n_gibbs = n_gibbs
  )
  attr(out, "latent_sampler") <- "finite_gibbs_diagnostic"
  attr(out, "gibbs_sweeps") <- as.integer(n_gibbs)
  out
}

sample_eiv_test_y_ordprobit_fb <- function(X_test_raw,
                                           C_test,
                                           fit_obj,
                                           draw_ids = NULL,
                                           n_per_draw = 1L,
                                           n_new_latent_gibbs = 100L,
                                           latent_sampler = c(
                                             "minimax_tilting", "rejection", "gibbs"
                                           ),
                                           rejection_batch_size = NULL,
                                           rejection_max_batches = 1000L,
                                           U_test_obs = NULL,
                                           U_input_scale = c("raw", "model"),
                                           predictive_target = c("response", "latent"),
                                           joint = FALSE,
                                           seed = NULL) {
  predictive_target <- match.arg(predictive_target)
  latent_sampler <- match.arg(latent_sampler)
  U_input_scale <- match.arg(U_input_scale)
  if (!is.null(seed)) set.seed(seed)
  
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
  p <- fit_obj$data$p
  n_train <- nrow(X_train)

  fit_level_maps <- fit_obj$data$C_level_maps
  if (is.null(fit_level_maps)) {
    fit_level_maps <- lapply(m_vec, seq_len)
  }
  fit_C_names <- fit_obj$data$C_names
  if (is.null(fit_C_names)) fit_C_names <- colnames(fit_obj$data$C_ord)
  fit_U_center <- fit_obj$data$U_center
  if (is.null(fit_U_center)) fit_U_center <- rep(0, d)
  fit_U_scale <- fit_obj$data$U_scale
  if (is.null(fit_U_scale)) fit_U_scale <- rep(1, d)
  kernel_spec <- kernel_spec_from_fit(fit_obj)

  ordinal_test <- prepare_ordinal_matrix(
    C_test,
    m_vec = m_vec,
    level_maps = fit_level_maps,
    expected_names = fit_C_names,
    name = "C_test"
  )
  C_test <- ordinal_test$C
  n_test <- nrow(C_test)

  X_test_names <- if (is.null(dim(X_test_raw))) names(X_test_raw) else colnames(X_test_raw)
  X_test_raw <- as_numeric_matrix_strict(
    X_test_raw,
    "X_test_raw",
    nrow_expected = n_test,
    ncol_expected = p
  )
  train_x_names <- colnames(fit_obj$data$X_raw)
  if (!is.null(X_test_names) && all(nzchar(X_test_names)) &&
      !is.null(train_x_names)) {
    if (!setequal(X_test_names, train_x_names)) {
      stop("X_test_raw column names do not match the training X columns.")
    }
    X_test_raw <- X_test_raw[, train_x_names, drop = FALSE]
  }

  U_test_complete <- integer(0)
  U_test_standardized <- matrix(NA_real_, n_test, d)
  if (!is.null(U_test_obs)) {
    U_test_names <- if (is.null(dim(U_test_obs))) names(U_test_obs) else colnames(U_test_obs)
    U_test_obs <- as_numeric_matrix_strict(
      U_test_obs,
      "U_test_obs",
      nrow_expected = n_test,
      ncol_expected = d,
      allow_na = TRUE
    )
    if (!is.null(U_test_names) && all(nzchar(U_test_names)) &&
        !is.null(fit_obj$data$U_names)) {
      if (!setequal(U_test_names, fit_obj$data$U_names)) {
        stop("U_test_obs column names do not match the training U columns.")
      }
      U_test_obs <- U_test_obs[, fit_obj$data$U_names, drop = FALSE]
    }
    observed_per_row <- rowSums(!is.na(U_test_obs))
    if (any(observed_per_row > 0L & observed_per_row < d)) {
      stop(
        "Partially observed rows of U_test_obs are not supported; provide all ",
        "d coordinates for a row or leave the entire row missing."
      )
    }
    U_test_complete <- which(observed_per_row == d)
    if (length(U_test_complete) > 0L) {
      U_test_standardized[U_test_complete, ] <- if (U_input_scale == "raw") {
        sweep(
          sweep(
            U_test_obs[U_test_complete, , drop = FALSE],
            2,
            fit_U_center,
            "-"
          ),
          2,
          fit_U_scale,
          "/"
        )
      } else {
        U_test_obs[U_test_complete, , drop = FALSE]
      }
    }
  }

  if (is.null(draw_ids)) draw_ids <- seq_len(nrow(samples_logtheta))
  draw_ids <- as.integer(draw_ids)
  if (length(draw_ids) < 1L || anyNA(draw_ids) ||
      any(draw_ids < 1L | draw_ids > nrow(samples_logtheta))) {
    stop("draw_ids contains invalid posterior-draw indices.")
  }
  n_per_draw <- as.integer(n_per_draw)
  if (length(n_per_draw) != 1L || is.na(n_per_draw) || n_per_draw < 1L) {
    stop("n_per_draw must be a positive integer.")
  }
  
  X_test <- sweep(sweep(X_test_raw, 2, X_center, "-"), 2, X_scale, "/")
  n_draw <- length(draw_ids) * n_per_draw
  
  out <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  conditional_means <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  conditional_vars <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  latent_acceptance <- rep(NA_real_, n_draw)
  latent_candidates <- integer(n_draw)
  
  row_id <- 0L
  
  for (s in draw_ids) {
    A_s <- matrix(samples_A[s, , ], nrow = q, ncol = d)
    tau_s <- unflatten_tau(samples_tau[s, ], m_vec)
    
    U_train_s <- matrix(samples_U[s, , ], nrow = n_train, ncol = d)
    
    for (rr in seq_len(n_per_draw)) {
      row_id <- row_id + 1L
      
      missing_test <- setdiff(seq_len(n_test), U_test_complete)
      U_star <- U_test_standardized
      if (length(missing_test) > 0L) {
        U_missing <- sample_u_given_c_ordprobit_dispatch(
          C_new = C_test[missing_test, , drop = FALSE],
          A = A_s,
          tau = tau_s,
          m_vec = m_vec,
          latent_sampler = latent_sampler,
          n_gibbs = n_new_latent_gibbs,
          rejection_batch_size = rejection_batch_size,
          rejection_max_batches = rejection_max_batches
        )
        U_star[missing_test, ] <- U_missing
        if (latent_sampler == "rejection") {
          latent_acceptance[row_id] <- attr(U_missing, "overall_acceptance")
          latent_candidates[row_id] <- attr(U_missing, "candidate_draws")
        }
      }
      if (length(U_test_complete) > 0L) {
        U_star[U_test_complete, ] <-
          U_test_standardized[U_test_complete, , drop = FALSE]
      }
      
      pred <- gp_predict_draw_general(
        X_train = X_train,
        U_train = U_train_s,
        y_train = y_train,
        X_star = X_test,
        U_star = U_star,
        logtheta = samples_logtheta[s, ],
        sigma2_eps = samples_sigma2[s],
        noisy = predictive_target == "response",
        kernel = kernel_spec$name,
        matern_nu = kernel_spec$matern_nu,
        return_cov = isTRUE(joint)
      )
      
      if (isTRUE(joint)) {
        y_std <- as.numeric(rmvnorm_psd(1L, pred$mean, pred$cov))
      } else {
        y_std <- pred$mean + sqrt(pred$var) * rnorm(n_test)
      }
      out[row_id, ] <- y_center + y_scale * y_std
      conditional_means[row_id, ] <- y_center + y_scale * pred$mean
      conditional_vars[row_id, ] <- y_scale^2 * pred$var
    }
  }
  
  attr(out, "predictive_target") <- predictive_target
  attr(out, "joint") <- isTRUE(joint)
  attr(out, "latent_sampler") <- latent_sampler
  attr(out, "latent_input_scale") <- U_input_scale
  attr(out, "rejection_acceptance") <- latent_acceptance
  attr(out, "rejection_candidate_draws") <- latent_candidates
  attr(out, "conditional_means") <- conditional_means
  attr(out, "conditional_vars") <- conditional_vars
  attr(out, "mixture_components") <-
    "posterior and latent-input Monte Carlo Gaussian components"
  out
}

posterior_u_draws <- function(fit_obj,
                              rows = NULL,
                              scale = c("auto", "raw", "model")) {
  scale <- match.arg(scale)
  anchored <- isTRUE(fit_obj$data$latent_scale_anchored)
  if (scale == "raw" && !anchored) {
    stop("scale='raw' requires calibration data that anchor latent U.")
  }
  if (scale == "auto") scale <- if (anchored) "raw" else "model"
  draws <- fit_obj$mcmc$samples_U
  n <- dim(draws)[2]
  d <- dim(draws)[3]

  if (is.null(rows)) rows <- seq_len(n)
  rows <- as.integer(rows)
  if (length(rows) < 1L || anyNA(rows) || any(rows < 1L | rows > n)) {
    stop("rows contains invalid training-row indices.")
  }
  draws <- draws[, rows, , drop = FALSE]

  if (scale == "raw") {
    U_center <- fit_obj$data$U_center
    if (is.null(U_center)) U_center <- rep(0, d)
    U_scale <- fit_obj$data$U_scale
    if (is.null(U_scale)) U_scale <- rep(1, d)
    for (k in seq_len(d)) {
      draws[, , k] <- U_center[k] + U_scale[k] * draws[, , k]
    }
  }

  draw_dimnames <- dimnames(draws)
  if (is.null(draw_dimnames)) draw_dimnames <- vector("list", 3L)
  draw_dimnames[[2]] <- as.character(rows)
  draw_dimnames[[3]] <- fit_obj$data$U_names
  dimnames(draws) <- draw_dimnames
  attr(draws, "scale") <- if (!anchored && scale == "model") "working" else scale
  attr(draws, "latent_scale_anchored") <- anchored
  draws
}

sample_eiv_f_given_xu_fb <- function(X_test_raw,
                                     U_test_raw,
                                     fit_obj,
                                     draw_ids = NULL,
                                     n_per_draw = 1L,
                                     include_process_uncertainty = TRUE,
                                     U_input_scale = c("raw", "model"),
                                     joint = FALSE,
                                     seed = NULL) {
  U_input_scale <- match.arg(U_input_scale)
  if (!is.null(seed)) set.seed(seed)

  p <- fit_obj$data$p
  d <- fit_obj$data$d
  X_input_names <- if (is.null(dim(X_test_raw))) names(X_test_raw) else colnames(X_test_raw)
  U_input_names <- if (is.null(dim(U_test_raw))) names(U_test_raw) else colnames(U_test_raw)

  X_test_raw <- as_numeric_matrix_strict(
    X_test_raw,
    "X_test_raw",
    ncol_expected = p
  )
  U_test_raw <- as_numeric_matrix_strict(
    U_test_raw,
    "U_test_raw",
    nrow_expected = nrow(X_test_raw),
    ncol_expected = d
  )

  train_x_names <- colnames(fit_obj$data$X_raw)
  if (!is.null(X_input_names) && all(nzchar(X_input_names)) &&
      !is.null(train_x_names)) {
    if (!setequal(X_input_names, train_x_names)) {
      stop("X_test_raw column names do not match the training X columns.")
    }
    X_test_raw <- X_test_raw[, train_x_names, drop = FALSE]
  }
  train_u_names <- fit_obj$data$U_names
  if (!is.null(U_input_names) && all(nzchar(U_input_names)) &&
      !is.null(train_u_names)) {
    if (!setequal(U_input_names, train_u_names)) {
      stop("U_test_raw column names do not match the training U columns.")
    }
    U_test_raw <- U_test_raw[, train_u_names, drop = FALSE]
  }

  U_center <- fit_obj$data$U_center
  if (is.null(U_center)) U_center <- rep(0, d)
  U_scale <- fit_obj$data$U_scale
  if (is.null(U_scale)) U_scale <- rep(1, d)
  X_test <- sweep(
    sweep(X_test_raw, 2, fit_obj$data$X_center, "-"),
    2,
    fit_obj$data$X_scale,
    "/"
  )
  U_test <- if (U_input_scale == "raw") {
    sweep(sweep(U_test_raw, 2, U_center, "-"), 2, U_scale, "/")
  } else {
    U_test_raw
  }

  samples_logtheta <- fit_obj$mcmc$samples_logtheta
  if (is.null(draw_ids)) draw_ids <- seq_len(nrow(samples_logtheta))
  draw_ids <- as.integer(draw_ids)
  if (length(draw_ids) < 1L || anyNA(draw_ids) ||
      any(draw_ids < 1L | draw_ids > nrow(samples_logtheta))) {
    stop("draw_ids contains invalid posterior-draw indices.")
  }
  n_per_draw <- as.integer(n_per_draw)
  if (length(n_per_draw) != 1L || is.na(n_per_draw) || n_per_draw < 1L) {
    stop("n_per_draw must be a positive integer.")
  }

  kernel_spec <- kernel_spec_from_fit(fit_obj)
  n_train <- nrow(fit_obj$data$X)
  out <- matrix(NA_real_, length(draw_ids) * n_per_draw, nrow(X_test))

  row_id <- 0L
  for (ii in seq_along(draw_ids)) {
    s <- draw_ids[ii]
    U_train_s <- matrix(
      fit_obj$mcmc$samples_U[s, , ],
      nrow = n_train,
      ncol = d
    )
    pred <- gp_predict_draw_general(
      X_train = fit_obj$data$X,
      U_train = U_train_s,
      y_train = fit_obj$data$y,
      X_star = X_test,
      U_star = U_test,
      logtheta = samples_logtheta[s, ],
      sigma2_eps = fit_obj$mcmc$samples_sigma2[s],
      noisy = FALSE,
      kernel = kernel_spec$name,
      matern_nu = kernel_spec$matern_nu,
      return_cov = isTRUE(joint) && isTRUE(include_process_uncertainty)
    )

    for (rr in seq_len(n_per_draw)) {
      row_id <- row_id + 1L
      if (!isTRUE(include_process_uncertainty)) {
        f_std <- pred$mean
      } else if (isTRUE(joint)) {
        f_std <- as.numeric(rmvnorm_psd(1L, pred$mean, pred$cov))
      } else {
        f_std <- pred$mean + sqrt(pred$var) * rnorm(nrow(X_test))
      }
      out[row_id, ] <- fit_obj$data$y_center + fit_obj$data$y_scale * f_std
    }
  }

  attr(out, "joint") <- isTRUE(joint)
  attr(out, "include_process_uncertainty") <- isTRUE(include_process_uncertainty)
  attr(out, "latent_input_scale") <- U_input_scale
  out
}

gp_integrated_mean_state_general <- function(X_train,
                                             U_train,
                                             y_train,
                                             X_star,
                                             U_mc,
                                             logtheta,
                                             sigma2_eps,
                                             kernel = "se",
                                             matern_nu = 2.5,
                                             return_cov = FALSE) {
  X_train <- as.matrix(X_train)
  U_train <- as.matrix(U_train)
  X_star <- as.matrix(X_star)
  U_mc <- as.array(U_mc)
  n <- nrow(X_train)
  N <- nrow(X_star)
  p <- ncol(X_train)
  d <- ncol(U_train)
  if (length(dim(U_mc)) != 3L || dim(U_mc)[2L] != N || dim(U_mc)[3L] != d) {
    stop("U_mc must be an M by N by d array.")
  }
  M <- dim(U_mc)[1L]
  if (nrow(U_train) != n || length(y_train) != n || ncol(X_star) != p ||
      M < 1L || any(!is.finite(c(X_train, U_train, y_train, X_star, U_mc)))) {
    stop("Incompatible or nonfinite inputs in gp_integrated_mean_state_general().")
  }

  kernel_spec <- normalize_gp_kernel(kernel, matern_nu)
  rho <- exp(logtheta[1L])
  theta_x <- exp(logtheta[1L + seq_len(p)])
  theta_u <- exp(logtheta[1L + p + seq_len(d)])
  D2_train <- weighted_sqdist_general(X_train, U_train, theta_x, theta_u)
  R_train <- kernel_from_weighted_sqdist(
    D2_train, kernel_spec$name, kernel_spec$matern_nu
  )
  C_train <- sigma2_eps * (rho^2 * R_train + diag(n))
  chol_train <- safe_chol(C_train)
  alpha <- solve_chol(chol_train, y_train)

  K_bar <- matrix(NA_real_, N, n)
  prior_diag <- numeric(N)
  for (i in seq_len(N)) {
    U_i <- matrix(U_mc[, i, ], nrow = M, ncol = d)
    X_i <- X_star[rep(i, M), , drop = FALSE]
    D2_cross <- weighted_sqdist_general(
      X_i, U_i, theta_x, theta_u,
      X_b = X_train, U_b = U_train
    )
    K_cross <- rho^2 * sigma2_eps * kernel_from_weighted_sqdist(
      D2_cross, kernel_spec$name, kernel_spec$matern_nu
    )
    K_bar[i, ] <- colMeans(K_cross)

    D2_ii <- weighted_sqdist_general(X_i, U_i, theta_x, theta_u)
    prior_diag[i] <- mean(
      rho^2 * sigma2_eps * kernel_from_weighted_sqdist(
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
      U_i <- matrix(U_mc[, i, ], nrow = M, ncol = d)
      X_i <- X_star[rep(i, M), , drop = FALSE]
      for (j in (i + 1L):N) {
        U_j <- matrix(U_mc[, j, ], nrow = M, ncol = d)
        X_j <- X_star[rep(j, M), , drop = FALSE]
        D2_ij <- weighted_sqdist_general(
          X_i, U_i, theta_x, theta_u,
          X_b = X_j, U_b = U_j
        )
        value <- mean(
          rho^2 * sigma2_eps * kernel_from_weighted_sqdist(
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

sample_eiv_m_given_xc_fb <- function(X_test_raw,
                                      C_test,
                                      fit_obj,
                                      draw_ids = NULL,
                                      n_per_draw = 1L,
                                      n_latent = 256L,
                                      include_process_uncertainty = TRUE,
                                      joint = FALSE,
                                      return_components = FALSE,
                                      latent_sampler = c(
                                        "minimax_tilting", "rejection", "gibbs"
                                      ),
                                      n_new_latent_gibbs = 100L,
                                      rejection_batch_size = NULL,
                                      rejection_max_batches = 1000L,
                                      seed = NULL) {
  latent_sampler <- match.arg(latent_sampler)
  if (!is.null(seed)) set.seed(as.integer(seed))
  p <- fit_obj$data$p
  d <- fit_obj$data$d
  q <- fit_obj$data$q
  m_vec <- fit_obj$data$m_vec
  X_input_names <- if (is.null(dim(X_test_raw))) names(X_test_raw) else colnames(X_test_raw)
  X_test_raw <- as_numeric_matrix_strict(
    X_test_raw, "X_test_raw", ncol_expected = p
  )
  train_x_names <- colnames(fit_obj$data$X_raw)
  if (!is.null(X_input_names) && all(nzchar(X_input_names)) &&
      !is.null(train_x_names)) {
    if (!setequal(X_input_names, train_x_names)) {
      stop("X_test_raw column names do not match the training X columns.")
    }
    X_test_raw <- X_test_raw[, train_x_names, drop = FALSE]
  }
  level_maps <- fit_obj$data$C_level_maps
  if (is.null(level_maps)) level_maps <- lapply(m_vec, seq_len)
  C_names <- fit_obj$data$C_names
  if (is.null(C_names)) C_names <- colnames(fit_obj$data$C_ord)
  C_test <- prepare_ordinal_matrix(
    C_test,
    m_vec = m_vec,
    level_maps = level_maps,
    expected_names = C_names,
    name = "C_test"
  )$C
  N <- nrow(C_test)
  if (nrow(X_test_raw) != N || ncol(C_test) != q) {
    stop("X_test_raw and C_test have incompatible dimensions.")
  }
  X_test <- sweep(
    sweep(X_test_raw, 2L, fit_obj$data$X_center, "-"),
    2L, fit_obj$data$X_scale, "/"
  )

  n_saved <- nrow(fit_obj$mcmc$samples_logtheta)
  if (is.null(draw_ids)) draw_ids <- seq_len(n_saved)
  draw_ids <- as.integer(draw_ids)
  if (length(draw_ids) < 1L || anyNA(draw_ids) ||
      any(draw_ids < 1L | draw_ids > n_saved)) {
    stop("draw_ids contains invalid posterior-draw indices.")
  }
  n_per_draw <- as.integer(n_per_draw)
  n_latent <- as.integer(n_latent)
  if (length(n_per_draw) != 1L || is.na(n_per_draw) || n_per_draw < 1L ||
      length(n_latent) != 1L || is.na(n_latent) || n_latent < 1L) {
    stop("n_per_draw and n_latent must be positive integers.")
  }

  kernel_spec <- kernel_spec_from_fit(fit_obj)
  n_train <- nrow(fit_obj$data$X)
  out <- matrix(NA_real_, length(draw_ids) * n_per_draw, N)
  conditional_means <- matrix(NA_real_, length(draw_ids), N)
  conditional_vars <- matrix(NA_real_, length(draw_ids), N)
  acceptance <- numeric(length(draw_ids))
  row_id <- 0L
  C_repeated <- C_test[rep(seq_len(N), each = n_latent), , drop = FALSE]
  for (ii in seq_along(draw_ids)) {
    s <- draw_ids[ii]
    A_s <- matrix(fit_obj$mcmc$samples_A[s, , ], nrow = q, ncol = d)
    tau_s <- unflatten_tau(fit_obj$mcmc$samples_tau[s, ], m_vec)
    U_flat <- sample_u_given_c_ordprobit_dispatch(
      C_new = C_repeated,
      A = A_s,
      tau = tau_s,
      m_vec = m_vec,
      latent_sampler = latent_sampler,
      n_gibbs = n_new_latent_gibbs,
      rejection_batch_size = rejection_batch_size,
      rejection_max_batches = rejection_max_batches
    )
    acceptance[ii] <- if (latent_sampler == "rejection") {
      attr(U_flat, "overall_acceptance")
    } else {
      NA_real_
    }
    U_mc <- array(NA_real_, dim = c(n_latent, N, d))
    for (k in seq_len(d)) {
      U_mc[, , k] <- matrix(U_flat[, k], nrow = n_latent, ncol = N)
    }
    U_train_s <- matrix(
      fit_obj$mcmc$samples_U[s, , ], nrow = n_train, ncol = d
    )
    pred <- gp_integrated_mean_state_general(
      X_train = fit_obj$data$X,
      U_train = U_train_s,
      y_train = fit_obj$data$y,
      X_star = X_test,
      U_mc = U_mc,
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
  attr(out, "latent_sampler") <- latent_sampler
  attr(out, "rejection_acceptance") <- acceptance
  if (isTRUE(return_components)) {
    attr(out, "conditional_means") <- conditional_means
    attr(out, "conditional_vars") <- conditional_vars
  }
  out
}

sample_eiv_u_given_c_ordprobit <- function(C_new,
                                            fit_obj,
                                            draw_ids = NULL,
                                            n_per_draw = 1L,
                                            scale = c("auto", "raw", "model"),
                                            latent_sampler = c(
                                              "minimax_tilting", "rejection", "gibbs"
                                            ),
                                            n_new_latent_gibbs = 100L,
                                            rejection_batch_size = NULL,
                                            rejection_max_batches = 1000L,
                                            seed = NULL) {
  scale <- match.arg(scale)
  anchored <- isTRUE(fit_obj$data$latent_scale_anchored)
  if (scale == "raw" && !anchored) {
    stop("scale='raw' requires calibration data that anchor latent U.")
  }
  if (scale == "auto") scale <- if (anchored) "raw" else "model"
  latent_sampler <- match.arg(latent_sampler)
  if (!is.null(seed)) set.seed(as.integer(seed))
  m_vec <- fit_obj$data$m_vec
  q <- fit_obj$data$q
  d <- fit_obj$data$d
  level_maps <- fit_obj$data$C_level_maps
  if (is.null(level_maps)) level_maps <- lapply(m_vec, seq_len)
  C_names <- fit_obj$data$C_names
  if (is.null(C_names)) C_names <- colnames(fit_obj$data$C_ord)
  C_new <- prepare_ordinal_matrix(
    C_new,
    m_vec = m_vec,
    level_maps = level_maps,
    expected_names = C_names,
    name = "C_new"
  )$C
  if (ncol(C_new) != q) stop("C_new has the wrong number of ordinal inputs.")
  n_saved <- nrow(fit_obj$mcmc$samples_logtheta)
  if (is.null(draw_ids)) draw_ids <- seq_len(n_saved)
  draw_ids <- as.integer(draw_ids)
  if (length(draw_ids) < 1L || anyNA(draw_ids) ||
      any(draw_ids < 1L | draw_ids > n_saved)) {
    stop("draw_ids contains invalid posterior-draw indices.")
  }
  n_per_draw <- as.integer(n_per_draw)
  if (length(n_per_draw) != 1L || is.na(n_per_draw) || n_per_draw < 1L) {
    stop("n_per_draw must be a positive integer.")
  }

  N <- nrow(C_new)
  out <- array(
    NA_real_,
    dim = c(length(draw_ids) * n_per_draw, N, d),
    dimnames = list(NULL, NULL, fit_obj$data$U_names)
  )
  acceptance <- numeric(length(draw_ids) * n_per_draw)
  row_id <- 0L
  for (s in draw_ids) {
    A_s <- matrix(fit_obj$mcmc$samples_A[s, , ], nrow = q, ncol = d)
    tau_s <- unflatten_tau(fit_obj$mcmc$samples_tau[s, ], m_vec)
    for (rr in seq_len(n_per_draw)) {
      row_id <- row_id + 1L
      U_new <- sample_u_given_c_ordprobit_dispatch(
        C_new = C_new,
        A = A_s,
        tau = tau_s,
        m_vec = m_vec,
        latent_sampler = latent_sampler,
        n_gibbs = n_new_latent_gibbs,
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
  if (scale == "raw") {
    center <- fit_obj$data$U_center
    if (is.null(center)) center <- rep(0, d)
    scale_value <- fit_obj$data$U_scale
    if (is.null(scale_value)) scale_value <- rep(1, d)
    for (k in seq_len(d)) {
      out[, , k] <- center[k] + scale_value[k] * out[, , k]
    }
  }
  attr(out, "source") <- "prospective"
  attr(out, "latent_sampler") <- latent_sampler
  attr(out, "rejection_acceptance") <- acceptance
  attr(out, "scale") <- if (!anchored && scale == "model") "working" else scale
  attr(out, "latent_scale_anchored") <- anchored
  out
}

############################################################
## Study II synthetic data-generating mechanism
############################################################

STUDY2_DESIGN_TAG <- paste(
  "study2-manuscript-v10-tailstable",
  "q4-d2",
  "A-fixed",
  "balanced-4-level",
  "sigma0.12",
  "exact-minimax-interwoven",
  "common-random-numbers",
  sep = "_"
)

study2_scenario_spec <- function(
    scenario = c(
      "primary",
      "latent_additive_control",
      "high_uncertainty",
      "logistic_misspec"
    )) {
  scenario <- match.arg(scenario)

  switch(
    scenario,
    primary = list(
      scenario = scenario,
      response_interactions = TRUE,
      score_error = "gaussian",
      lambda = 1.0
    ),
    latent_additive_control = list(
      scenario = scenario,
      response_interactions = FALSE,
      score_error = "gaussian",
      lambda = 1.0
    ),
    high_uncertainty = list(
      scenario = scenario,
      response_interactions = TRUE,
      score_error = "gaussian",
      lambda = 1.5
    ),
    logistic_misspec = list(
      scenario = scenario,
      response_interactions = TRUE,
      score_error = "logistic",
      lambda = 1.0
    )
  )
}

make_study2_balanced_tau <- function(A, lambda, score_error, m = 4L) {
  m <- as.integer(m)
  if (length(m) != 1L || is.na(m) || m < 2L) {
    stop("m must be one integer at least two.")
  }
  probs <- seq_len(m - 1L) / m

  if (score_error == "gaussian") {
    score_sd <- sqrt(rowSums(A^2) + lambda^2)
    return(outer(score_sd, qnorm(probs)))
  }

  ## For the robustness setting, each score is the convolution of a centered
  ## Gaussian loading term and a variance-one logistic error. Compute fixed
  ## marginal quantiles numerically so that changing the error distribution is
  ## not confounded with category imbalance.
  logistic_scale <- sqrt(3) / base::pi

  tau <- t(vapply(rowSums(A^2), function(var_normal) {
    cdf_score <- function(t) {
      if (var_normal <= .Machine$double.eps) {
        return(plogis(t, scale = lambda * logistic_scale))
      }

      integrate(
        function(z) {
          plogis(
            (t - sqrt(var_normal) * z) / (lambda * logistic_scale)
          ) * dnorm(z)
        },
        lower = -9,
        upper = 9,
        rel.tol = 1e-9
      )$value
    }

    vapply(probs, function(prob) {
      if (prob == 0.5) return(0)
      uniroot(
        function(t) cdf_score(t) - prob,
        interval = c(-10, 10),
        tol = 1e-9
      )$root
    }, numeric(1))
  }, numeric(length(probs))))

  tau
}

make_study2_loading_matrix <- function(q = 4L) {
  q <- as.integer(q)
  if (length(q) != 1L || is.na(q) || q < 2L || q > 6L) {
    stop("The publication proxy-dimension design requires 2 <= q <= 6.")
  }

  ## The first two rows identify the two latent axes. Later rows are mixed or
  ## redundant measurements of the same state. Designs with q = 2, 4, and 6
  ## are therefore nested rather than unrelated DGMs.
  A_master <- rbind(
    c(1.70, 0.00),
    c(0.20, 1.70),
    c(1.20, 0.70),
    c(0.70, 1.20),
    c(1.45, -0.45),
    c(-0.45, 1.45)
  )
  A_master[seq_len(q), , drop = FALSE]
}

make_study2_true_params <- function(scenario = "primary", q = 4L, m = 4L) {
  spec <- study2_scenario_spec(scenario)
  q <- as.integer(q)
  m <- as.integer(m)
  A <- make_study2_loading_matrix(q)
  Omega <- diag(q)
  tau <- make_study2_balanced_tau(
    A = A,
    lambda = spec$lambda,
    score_error = spec$score_error,
    m = m
  )

  c(
    spec,
    list(
      A = A,
      Omega = Omega,
      tau = tau,
      q = q,
      m = m,
      d = 2L,
      sigma_eps = 0.12
    )
  )
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

f0_2d <- function(X, U, scenario = "primary") {
  X <- as.matrix(X)
  U <- as.matrix(U)
  spec <- study2_scenario_spec(scenario)
  
  x1 <- X[, 1]
  x2 <- X[, 2]
  u1 <- U[, 1]
  u2 <- U[, 2]
  
  f_main <- 0.65 * tanh(1.10 * u1) +
    0.50 * tanh(0.90 * u2) +
    0.30 * x1 -
    0.25 * x2 +
    0.10 * x1 * x2

  if (!isTRUE(spec$response_interactions)) {
    return(f_main)
  }

  f_main +
    0.30 * tanh(0.60 * u1 * u2) +
    0.25 * x1 * tanh(u1) -
    0.20 * x2 * tanh(u2)
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
                                 scenario = "primary",
                                 sigma_eps = NULL,
                                 q = 4L,
                                 m = 4L) {

  n <- as.integer(n)
  n_test <- as.integer(n_test)
  if (n < 1L || n_test < 1L) {
    stop("n and n_test must be positive integers.")
  }
  
  seed <- if (is.null(seed)) sample.int(.Machine$integer.max, 1L) else as.integer(seed)
  pars <- make_study2_true_params(scenario = scenario, q = q, m = m)
  A <- pars$A
  tau <- pars$tau
  lambda <- pars$lambda
  d <- pars$d

  if (is.null(sigma_eps)) sigma_eps <- pars$sigma_eps
  
  ## Generate exactly once.  The experiment must retain naturally sparse
  ## realized patterns rather than condition on a favorable category count.
  set.seed(seed + 1L)
  X <- maximin_lhs_nd(n, d = 2L, lower = -1, upper = 1)
  set.seed(seed + 2L)
  U <- matrix(rnorm(n * d), n, d)
  set.seed(seed + 3L)
  eps_y <- rnorm(n, 0, sigma_eps)
  set.seed(seed + 4L)
  X_test <- matrix(runif(n_test * 2L, -1, 1), n_test, 2L)
  set.seed(seed + 5L)
  U_test <- matrix(rnorm(n_test * d), n_test, d)
  set.seed(seed + 6L)
  eps_y_test <- rnorm(n_test, 0, sigma_eps)

  set.seed(seed + 7L)
  Zeta <- draw_study2_score_errors(n, pars)
  S <- U %*% t(A) + lambda * Zeta
  C <- ordinal_from_scores_matrix_tau(S, tau)
  
  f <- f0_2d(X, U, scenario = scenario)
  y <- f + eps_y

  set.seed(seed + 8L)
  Zeta_test <- draw_study2_score_errors(n_test, pars)
  S_test <- U_test %*% t(A) + lambda * Zeta_test
  C_test <- ordinal_from_scores_matrix_tau(S_test, tau)
  
  f_test <- f0_2d(X_test, U_test, scenario = scenario)
  y_test <- f_test + eps_y_test
  
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
    sigma_eps = sigma_eps,
    scenario = scenario
  )
}

classify_study2_pattern_frequency <- function(C_train,
                                              C_test,
                                              common_min = 5L) {
  train_key <- pattern_key(C_train)
  test_key <- pattern_key(C_test)
  train_counts <- table(train_key)
  n_train_pattern <- as.integer(train_counts[test_key])
  n_train_pattern[is.na(n_train_pattern)] <- 0L

  stratum <- ifelse(
    n_train_pattern == 0L,
    "unobserved",
    ifelse(n_train_pattern < common_min, "rare", "common")
  )

  data.frame(
    pattern = test_key,
    n_train_pattern = n_train_pattern,
    pattern_stratum = factor(
      stratum,
      levels = c("common", "rare", "unobserved")
    ),
    stringsAsFactors = FALSE
  )
}

make_stratified_calibration_sets_2d <- function(C,
                                                calib_grid = c(0, 10, 25, 50, 80),
                                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  C <- as.matrix(C)
  n <- nrow(C)
  calib_grid <- as.integer(calib_grid)
  if (ncol(C) < 2L) {
    stop("Study II calibration stratification requires at least two proxies.")
  }
  if (any(calib_grid < 0L | calib_grid > n)) {
    stop("Calibration sizes must lie between 0 and n.")
  }

  ## The first two proxies primarily anchor u1 and u2.  Randomize units within
  ## each observed (c1,c2) cell, randomize the order of the observed cells, and
  ## cycle through cells.  Prefixes of the resulting ordering give nested
  ## calibration sets without consulting y or the latent truth.
  cell <- interaction(C[, 1], C[, 2], drop = TRUE, lex.order = TRUE)
  idx_by_stratum <- split(seq_len(n), cell)
  idx_by_stratum <- idx_by_stratum[sample(seq_along(idx_by_stratum))]
  idx_by_stratum <- lapply(idx_by_stratum, function(idx) {
    if (length(idx) <= 1L) idx else sample(idx, size = length(idx), replace = FALSE)
  })
  
  ordering <- integer(0)
  
  repeat {
    added <- FALSE
    
    for (ss in seq_along(idx_by_stratum)) {
      if (length(idx_by_stratum[[ss]]) > 0) {
        ordering <- c(ordering, idx_by_stratum[[ss]][1])
        idx_by_stratum[[ss]] <- idx_by_stratum[[ss]][-1]
        added <- TRUE
      }
    }
    
    if (!added) break
  }
  
  out <- lapply(calib_grid, function(k) {
    if (k == 0) integer(0) else sort(ordering[seq_len(k)])
  })
  stopifnot(
    length(ordering) == n,
    length(unique(ordering)) == n,
    all(ordering %in% seq_len(n)),
    all(vapply(out, function(idx) length(idx) == length(unique(idx)), logical(1)))
  )
  
  names(out) <- as.character(calib_grid)
  out
}

make_nested_calibration_sets_2d <- function(C,
                                            calib_grid,
                                            seed = NULL,
                                            scheme = c("random", "anchor_stratified")) {
  scheme <- match.arg(scheme)
  C <- as.matrix(C)
  n <- nrow(C)
  calib_grid <- as.integer(calib_grid)
  if (anyNA(calib_grid) || any(calib_grid < 0L | calib_grid > n)) {
    stop("Calibration sizes must lie between zero and the training size.")
  }
  if (scheme == "anchor_stratified") {
    return(make_stratified_calibration_sets_2d(C, calib_grid, seed))
  }

  if (!is.null(seed)) set.seed(seed)
  ordering <- sample.int(n)
  out <- lapply(calib_grid, function(k) {
    if (k == 0L) integer(0) else sort(ordering[seq_len(k)])
  })
  names(out) <- as.character(calib_grid)
  out
}

############################################################
## Oracle prediction for Study II simulation
############################################################

study2_oracle_pool_provenance <- function(true_params) {
  list(
    schema_version = "study2-oracle-pool-1.0.0",
    scenario = as.character(true_params$scenario),
    score_error = as.character(true_params$score_error),
    lambda = as.numeric(true_params$lambda),
    A = unname(as.matrix(true_params$A)),
    tau = unname(as.matrix(true_params$tau)),
    Omega = unname(as.matrix(true_params$Omega)),
    d = as.integer(true_params$d),
    q = as.integer(true_params$q),
    m = as.integer(true_params$m)
  )
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

make_oracle_pool_2d <- function(true_params,
                                n_pool = 200000L,
                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  A <- true_params$A
  tau <- true_params$tau
  lambda <- true_params$lambda
  
  d <- true_params$d
  U <- matrix(rnorm(n_pool * d), n_pool, d)
  Zeta <- draw_study2_score_errors(n_pool, true_params)
  
  S <- U %*% t(A) + lambda * Zeta
  C <- ordinal_from_scores_matrix_tau(S, tau)
  
  key <- pattern_key(C)
  split_idx <- split(seq_len(n_pool), key)
  
  list(
    U = U,
    C = C,
    key = key,
    split_idx = split_idx,
    provenance = study2_oracle_pool_provenance(true_params)
  )
}

sample_oracle_u_rejection <- function(c_star,
                                      true_params,
                                      n_draw,
                                      max_batches = 100L,
                                      batch_size = 20000L) {
  A <- true_params$A
  tau <- true_params$tau
  lambda <- true_params$lambda
  
  d <- true_params$d
  U_keep <- matrix(NA_real_, 0, d)
  key_star <- paste(c_star, collapse = "_")
  
  for (bb in seq_len(max_batches)) {
    U <- matrix(rnorm(batch_size * d), batch_size, d)
    Zeta <- draw_study2_score_errors(batch_size, true_params)
    S <- U %*% t(A) + lambda * Zeta
    C <- ordinal_from_scores_matrix_tau(S, tau)
    
    idx <- which(pattern_key(C) == key_star)
    
    if (length(idx) > 0) {
      U_keep <- rbind(U_keep, U[idx, , drop = FALSE])
    }
    
    if (nrow(U_keep) >= n_draw) break
  }
  
  if (nrow(U_keep) < n_draw) {
    stop(
      "Oracle rejection sampling produced only ", nrow(U_keep), " of ",
      n_draw, " required draws for ordinal pattern ",
      key_star,
      ". Increase max_batches or batch_size; do not bootstrap a short pool."
    )
  }
  
  U_keep[seq_len(n_draw), , drop = FALSE]
}

sample_oracle_u_patterns_shared_2d <- function(
    C,
    true_params,
    n_per_pattern,
    oracle_pool = NULL,
    n_pool = 200000L,
    batch_size = 200000L,
    max_candidates = 20000000L,
    seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  C <- as.matrix(C)
  if (nrow(C) < 1L || ncol(C) != true_params$q || anyNA(C) ||
      any(!is.finite(C)) || any(C != floor(C)) ||
      any(C < 1L | C > true_params$m)) {
    stop("C is incompatible with the Study II oracle model.")
  }

  keys <- unique(pattern_key(C))
  if (length(n_per_pattern) == 1L) {
    required <- setNames(rep(as.integer(n_per_pattern), length(keys)), keys)
  } else {
    if (is.null(names(n_per_pattern)) || !all(keys %in% names(n_per_pattern))) {
      stop("A vector n_per_pattern must be named for every requested pattern.")
    }
    required <- as.integer(n_per_pattern[keys])
    names(required) <- keys
  }
  if (anyNA(required) || any(required < 1L)) {
    stop("Every requested oracle pattern must require at least one draw.")
  }

  n_pool <- as.integer(n_pool)
  batch_size <- as.integer(batch_size)
  max_candidates <- as.numeric(max_candidates)
  if (is.na(n_pool) || n_pool < 1L || is.na(batch_size) || batch_size < 1L ||
      !is.finite(max_candidates) || max_candidates < batch_size) {
    stop("Oracle rejection budgets are invalid.")
  }
  if (is.null(oracle_pool)) {
    oracle_pool <- make_oracle_pool_2d(true_params, n_pool = n_pool)
  }
  validate_oracle_pool_2d(oracle_pool, true_params)

  accepted <- setNames(vector("list", length(keys)), keys)
  accepted_n <- setNames(integer(length(keys)), keys)
  proposal_hits <- setNames(integer(length(keys)), keys)
  initial_candidates <- nrow(oracle_pool$U)
  initial_pool_hits <- setNames(integer(length(keys)), keys)
  for (key in keys) {
    idx <- oracle_pool$split_idx[[key]]
    if (is.null(idx)) idx <- integer(0)
    initial_pool_hits[[key]] <- length(idx)
    proposal_hits[[key]] <- length(idx)
    if (length(idx) > 0L) {
      take <- idx[seq_len(min(length(idx), required[[key]]))]
      accepted[[key]] <- list(oracle_pool$U[take, , drop = FALSE])
      accepted_n[[key]] <- length(take)
    }
  }

  additional_candidates <- 0
  while (any(accepted_n < required)) {
    if (additional_candidates + batch_size > max_candidates) {
      missing <- keys[accepted_n < required]
      detail <- paste0(
        missing, " (", accepted_n[missing], "/", required[missing], ")"
      )
      stop(
        "Shared exact oracle rejection reached max_candidates before filling ",
        "pattern(s): ", paste(detail, collapse = ", "), "."
      )
    }
    U_batch <- matrix(
      rnorm(batch_size * true_params$d), batch_size, true_params$d
    )
    Zeta_batch <- draw_study2_score_errors(batch_size, true_params)
    S_batch <- U_batch %*% t(true_params$A) +
      true_params$lambda * Zeta_batch
    C_batch <- ordinal_from_scores_matrix_tau(S_batch, true_params$tau)
    batch_split <- split(seq_len(batch_size), pattern_key(C_batch))
    additional_candidates <- additional_candidates + batch_size

    for (key in keys) {
      idx <- batch_split[[key]]
      if (is.null(idx) || length(idx) == 0L) next
      proposal_hits[[key]] <- proposal_hits[[key]] + length(idx)
      need <- required[[key]] - accepted_n[[key]]
      if (need <= 0L) next
      take <- idx[seq_len(min(need, length(idx)))]
      accepted[[key]][[length(accepted[[key]]) + 1L]] <-
        U_batch[take, , drop = FALSE]
      accepted_n[[key]] <- accepted_n[[key]] + length(take)
    }
  }

  out <- setNames(lapply(keys, function(key) {
    do.call(rbind, accepted[[key]])[
      seq_len(required[[key]]), , drop = FALSE
    ]
  }), keys)
  attr(out, "rejection_telemetry") <- data.frame(
    pattern = keys,
    required = as.integer(required[keys]),
    retained = as.integer(accepted_n[keys]),
    initial_pool_hits = as.integer(initial_pool_hits[keys]),
    initial_candidates = initial_candidates,
    additional_candidates = additional_candidates,
    total_candidates = initial_candidates + additional_candidates,
    empirical_acceptance = as.integer(proposal_hits[keys]) /
      (initial_candidates + additional_candidates),
    stringsAsFactors = FALSE
  )
  out
}

oracle_m0_2d <- function(X,
                         C,
                         true_params,
                         oracle_pool = NULL,
                         n_latent = 5000L,
                         batch_size = 200000L,
                         max_candidates = 20000000L,
                         seed = NULL) {
  X <- as.matrix(X)
  C <- as.matrix(C)
  n_latent <- as.integer(n_latent)
  if (nrow(X) != nrow(C) || ncol(C) != true_params$q ||
      n_latent < 1L || any(!is.finite(X)) || anyNA(C)) {
    stop("Incompatible inputs in oracle_m0_2d().")
  }
  keys <- pattern_key(C)
  unique_keys <- unique(keys)
  pattern_draws <- sample_oracle_u_patterns_shared_2d(
    C = C,
    true_params = true_params,
    n_per_pattern = n_latent,
    oracle_pool = oracle_pool,
    n_pool = max(200000L, 20L * n_latent),
    batch_size = batch_size,
    max_candidates = max_candidates,
    seed = seed
  )

  out <- numeric(nrow(X))
  truth_mcse <- numeric(nrow(X))
  nested_half_difference <- numeric(nrow(X))
  for (key in unique_keys) {
    rows <- which(keys == key)
    U_draw <- pattern_draws[[key]]
    for (i in rows) {
      X_rep <- X[rep(i, n_latent), , drop = FALSE]
      f_value <- f0_2d(X_rep, U_draw, scenario = true_params$scenario)
      out[i] <- mean(f_value)
      truth_mcse[i] <- if (n_latent > 1L) {
        stats::sd(f_value) / sqrt(n_latent)
      } else {
        NA_real_
      }
      n_half <- floor(n_latent / 2L)
      nested_half_difference[i] <- if (n_half > 0L) {
        abs(mean(f_value[seq_len(n_half)]) - out[i])
      } else {
        NA_real_
      }
    }
  }
  attr(out, "n_latent_truth") <- n_latent
  attr(out, "rejection_telemetry") <-
    attr(pattern_draws, "rejection_telemetry")
  attr(out, "truth_diagnostics") <- data.frame(
    evaluation_row = seq_len(nrow(X)),
    pattern = keys,
    n_latent = n_latent,
    mcse = truth_mcse,
    nested_half_difference = nested_half_difference,
    stringsAsFactors = FALSE
  )
  out
}

summarize_mean_recovery_2d <- function(draw_mat,
                                       m_true,
                                       method,
                                       rep_id,
                                       n_calib,
                                       scenario,
                                       valid_function_draws = FALSE) {
  draw_mat <- as.matrix(draw_mat)
  m_true <- as.numeric(m_true)
  if (ncol(draw_mat) != length(m_true) || any(!is.finite(m_true))) {
    stop("draw_mat and m_true are incompatible in summarize_mean_recovery_2d().")
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

sample_oracle_test_y_2d <- function(X_test,
                                    C_test,
                                    true_params,
                                    sigma_eps,
                                    n_draw = 1000L,
                                    oracle_pool = NULL,
                                    n_pool = 200000L,
                                    batch_size = 200000L,
                                    max_candidates = 20000000L,
                                    seed = NULL) {
  X_test <- as.matrix(X_test)
  C_test <- as.matrix(C_test)
  n_draw <- as.integer(n_draw)
  if (nrow(X_test) != nrow(C_test) || ncol(C_test) != true_params$q ||
      n_draw < 1L || length(sigma_eps) != 1L || !is.finite(sigma_eps) ||
      sigma_eps < 0 || any(!is.finite(X_test)) || anyNA(C_test)) {
    stop("Incompatible inputs in sample_oracle_test_y_2d().")
  }

  n_test <- nrow(C_test)
  out <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  conditional_means <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  conditional_vars <- matrix(
    sigma_eps^2, nrow = n_draw, ncol = n_test
  )
  keys_test <- pattern_key(C_test)
  pattern_counts <- table(factor(keys_test, levels = unique(keys_test)))
  required <- n_draw * as.integer(pattern_counts)
  names(required) <- names(pattern_counts)
  pattern_draws <- sample_oracle_u_patterns_shared_2d(
    C = C_test,
    true_params = true_params,
    n_per_pattern = required,
    oracle_pool = oracle_pool,
    n_pool = n_pool,
    batch_size = batch_size,
    max_candidates = max_candidates,
    seed = seed
  )
  cursor <- setNames(integer(length(required)), names(required))

  for (j in seq_len(n_test)) {
    key_j <- keys_test[j]
    draw_idx <- cursor[[key_j]] + seq_len(n_draw)
    U_draw <- pattern_draws[[key_j]][draw_idx, , drop = FALSE]
    cursor[[key_j]] <- cursor[[key_j]] + n_draw

    X_rep <- matrix(
      rep(X_test[j, ], each = n_draw),
      nrow = n_draw,
      ncol = ncol(X_test)
    )
    f_draw <- f0_2d(
      X_rep,
      U_draw,
      scenario = true_params$scenario
    )
    
    out[, j] <- f_draw + rnorm(n_draw, 0, sigma_eps)
    conditional_means[, j] <- f_draw
  }

  attr(out, "latent_sampler") <- "exact_shared_rejection"
  attr(out, "rejection_telemetry") <-
    attr(pattern_draws, "rejection_telemetry")
  attr(out, "conditional_means") <- conditional_means
  attr(out, "conditional_vars") <- conditional_vars
  attr(out, "mixture_components") <- "DGM conditional on Monte Carlo U*"
  out
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

## Internal small-data GP engine used by the appendix Full-U, PI, and CC
## benchmarks as well as the archived deterministic-embedding pilots.  Every
## fit uses reproducible multistart optimization and fails closed unless at
## least one L-BFGS-B run reports convergence.
gp_mle_fit <- function(X,
                       y,
                       n_starts = 5L,
                       seed = 1L,
                       maxit = 500L) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  y <- as.numeric(y)

  n <- nrow(X)
  d <- ncol(X)
  if (n < 2L || d < 1L || length(y) != n || any(!is.finite(X)) ||
      any(!is.finite(y))) {
    stop("X and y must contain compatible finite GP training data.")
  }
  n_starts <- as.integer(n_starts)
  seed <- as.integer(seed)
  maxit <- as.integer(maxit)
  if (length(n_starts) != 1L || is.na(n_starts) || n_starts < 1L ||
      length(seed) != 1L || is.na(seed) ||
      length(maxit) != 1L || is.na(maxit) || maxit < 0L) {
    stop("n_starts, seed, and maxit must be finite nonnegative integers, with n_starts >= 1.")
  }
  
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
    
    0.5 * (n * log(2 * base::pi * sigma2) + logdetA + quad / sigma2)
  }
  
  lower <- c(log(1e-5), log(0.05), rep(log(1e-4), d))
  upper <- c(log(5), log(100), rep(log(100), d))

  caller_rng <- mixedgp_rng_state()
  on.exit(mixedgp_restore_rng_state(caller_rng), add = TRUE)
  set.seed(seed)
  initial <- matrix(NA_real_, nrow = n_starts, ncol = length(lower))
  initial[1L, ] <- c(log(0.05), log(3), rep(log(0.5), d))
  if (n_starts > 1L) {
    for (ss in 2:n_starts) {
      initial[ss, ] <- stats::runif(length(lower), lower, upper)
    }
  }

  attempts <- lapply(seq_len(n_starts), function(ss) {
    tryCatch(
      stats::optim(
        par = initial[ss, ],
        fn = nll,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        control = list(maxit = maxit)
      ),
      error = function(e) e
    )
  })
  parameter_names <- c(
    "log_sigma2", "log_rho", paste0("log_theta", seq_len(d))
  )
  initial_report <- as.data.frame(initial, stringsAsFactors = FALSE)
  names(initial_report) <- paste0("initial_", parameter_names)
  final_report <- as.data.frame(
    matrix(NA_real_, nrow = n_starts, ncol = length(parameter_names)),
    stringsAsFactors = FALSE
  )
  names(final_report) <- paste0("final_", parameter_names)
  report <- cbind(
    data.frame(
      start = seq_len(n_starts),
      n_starts = rep.int(n_starts, n_starts),
      optimizer_seed = rep.int(seed, n_starts),
      maxit = rep.int(maxit, n_starts),
      stringsAsFactors = FALSE
    ),
    initial_report,
    final_report,
    data.frame(
      convergence = rep.int(NA_integer_, n_starts),
      objective = rep.int(Inf, n_starts),
      function_evaluations = rep.int(NA_integer_, n_starts),
      gradient_evaluations = rep.int(NA_integer_, n_starts),
      message = rep.int("", n_starts),
      stringsAsFactors = FALSE
    )
  )
  for (ss in seq_along(attempts)) {
    attempt <- attempts[[ss]]
    if (inherits(attempt, "error")) {
      report$message[ss] <- conditionMessage(attempt)
    } else {
      report$convergence[ss] <- as.integer(attempt$convergence)
      report$objective[ss] <- as.numeric(attempt$value)
      if (!is.null(attempt$counts)) {
        function_count <- attempt$counts[["function"]]
        gradient_count <- attempt$counts[["gradient"]]
        if (length(function_count) == 1L) {
          report$function_evaluations[ss] <- as.integer(function_count)
        }
        if (length(gradient_count) == 1L) {
          report$gradient_evaluations[ss] <- as.integer(gradient_count)
        }
      }
      if (length(attempt$par) == length(parameter_names)) {
        report[ss, paste0("final_", parameter_names)] <- attempt$par
      }
      attempt_message <- attempt$message
      report$message[ss] <- if (
        is.null(attempt_message) || length(attempt_message) == 0L
      ) {
        ""
      } else {
        as.character(attempt_message[[1L]])
      }
    }
  }
  valid <- which(
    report$convergence == 0L & is.finite(report$objective) &
      report$objective < 1e19
  )
  if (length(valid) == 0L) {
    failure <- structure(
      list(
        message = paste0(
          "Internal Study II GP MLE did not converge from any of the ",
          n_starts, " prespecified starts."
        ),
        call = sys.call(),
        optimizer_attempts = report,
        optimizer_control = list(
          n_starts = n_starts, seed = seed, maxit = maxit
        )
      ),
      class = c("mixedgp_gp_mle_error", "error", "condition")
    )
    stop(failure)
  }
  best <- valid[which.min(report$objective[valid])]
  opt <- attempts[[best]]

  list(
    par = opt$par,
    value = opt$value,
    convergence = opt$convergence,
    optimizer_attempts = report,
    selected_start = best,
    optimizer_control = list(
      n_starts = n_starts, seed = seed, maxit = maxit
    ),
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

gp_mle_fit_learned_embedding_ord <- function(X,
                                             C_ord,
                                             y,
                                             m_vec,
                                             n_starts = 4L) {
  X <- as.matrix(X)
  C_ord <- as.matrix(C_ord)
  
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
    
    0.5 * (n * log(2 * base::pi * sigma2) + logdetA + quad / sigma2)
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

study2_latent_anchor_from_fit <- function(fit) {
  if (!is.list(fit) || !is.list(fit$data)) {
    stop("fit must contain a data component.")
  }
  rank <- fit$data$latent_anchor_rank
  required_rank <- fit$data$latent_anchor_required_rank
  anchored <- fit$data$latent_scale_anchored
  if (!is.null(rank) && !is.null(required_rank) && !is.null(anchored)) {
    return(list(
      anchored = isTRUE(anchored),
      affine_rank = as.integer(rank),
      required_rank = as.integer(required_rank),
      n_calibrated = length(fit$data$calib_idx)
    ))
  }
  if (is.null(fit$data$U_obs) || is.null(fit$data$calib_idx) ||
      is.null(fit$data$d)) {
    stop("fit does not contain enough information to assess latent anchoring.")
  }
  mixedgp_latent_anchor_status(
    fit$data$U_obs,
    fit$data$calib_idx,
    d = fit$data$d
  )
}

study2_latent_imputation_status <- function(fit,
                                            method,
                                            target,
                                            rep_id,
                                            n_calib,
                                            scenario,
                                            n_units) {
  target <- match.arg(
    target,
    c("training_missing_U", "prospective_U_given_C")
  )
  n_units <- as.integer(n_units)
  if (length(n_units) != 1L || is.na(n_units) || n_units < 0L) {
    stop("n_units must be one nonnegative integer.")
  }
  anchor <- study2_latent_anchor_from_fit(fit)
  if (!isTRUE(anchor$anchored)) {
    status <- "not_identified"
    reason <- paste0(
      "Raw-coordinate latent error is not identified: calibration affine ",
      "rank is ", anchor$affine_rank, " but ", anchor$required_rank,
      " is required."
    )
  } else if (n_units == 0L) {
    status <- "not_applicable"
    reason <- if (target == "training_missing_U") {
      "Every training latent input is calibrated."
    } else {
      "No prospective units were supplied."
    }
  } else {
    status <- "eligible"
    reason <- ""
  }
  data.frame(
    rep = as.integer(rep_id),
    scenario = as.character(scenario),
    n_calib = as.integer(n_calib),
    method = as.character(method),
    target = target,
    status = status,
    task_eligible = identical(status, "eligible"),
    reason = reason,
    n_units = n_units,
    latent_scale = if (isTRUE(anchor$anchored)) "physical" else "working",
    calibration_affine_rank = as.integer(anchor$affine_rank),
    calibration_required_rank = as.integer(anchor$required_rank),
    stringsAsFactors = FALSE
  )
}

summarize_study2_latent_draws <- function(draws,
                                          U_true,
                                          method,
                                          target,
                                          rep_id,
                                          n_calib,
                                          scenario) {
  target <- match.arg(
    target,
    c("training_missing_U", "prospective_U_given_C")
  )
  U_true <- as.matrix(U_true)
  draw_dim <- dim(draws)
  if (length(draw_dim) != 3L || draw_dim[1L] < 1L ||
      draw_dim[2L] != nrow(U_true) || draw_dim[3L] != ncol(U_true)) {
    stop("draws must have dimensions draw by unit by latent coordinate.")
  }
  if (nrow(U_true) < 1L || ncol(U_true) < 1L ||
      any(!is.finite(U_true)) || any(!is.finite(draws))) {
    stop("Latent draws and truth must be nonempty and finite.")
  }

  out <- lapply(seq_len(ncol(U_true)), function(j) {
    draw_j <- draws[, , j, drop = FALSE]
    dim(draw_j) <- c(draw_dim[1L], draw_dim[2L])
    post_mean <- colMeans(draw_j)
    lo <- apply(draw_j, 2L, stats::quantile, probs = 0.025, names = FALSE)
    hi <- apply(draw_j, 2L, stats::quantile, probs = 0.975, names = FALSE)
    err <- post_mean - U_true[, j]

    data.frame(
      rep = as.integer(rep_id),
      scenario = as.character(scenario),
      n_calib = as.integer(n_calib),
      method = as.character(method),
      target = target,
      coordinate = paste0("u", j),
      n_units = nrow(U_true),
      n_missing = if (target == "training_missing_U") {
        nrow(U_true)
      } else {
        NA_integer_
      },
      Bias = mean(err),
      RMSE = sqrt(mean(err^2)),
      MAE = mean(abs(err)),
      Coverage95 = mean(U_true[, j] >= lo & U_true[, j] <= hi),
      Width95 = mean(hi - lo),
      stringsAsFactors = FALSE
    )
  })
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

############################################################
## MCMC review helpers for ordinal-probit EIV-GP
############################################################

classify_parameter_block <- function(parameter) {
  out <- rep("other", length(parameter))
  
  out[grepl("^sigma|^rho|^theta", parameter)] <- "GP/noise"
  out[grepl("^A\\[|^A[0-9]", parameter)] <- "loading A"
  out[grepl("^tau\\[|^tau", parameter)] <- "cutpoints"
  out[grepl("^U\\[", parameter)] <- "latent U"
  
  out
}

summarise_rhat_block <- function(rhat_df) {
  rhat_df |>
    dplyr::filter(is.finite(rhat)) |>
    dplyr::group_by(block) |>
    dplyr::summarise(
      n = dplyr::n(),
      median = median(rhat),
      q90 = quantile(rhat, 0.90),
      q95 = quantile(rhat, 0.95),
      max = max(rhat),
      prop_gt_1.01 = mean(rhat > 1.01),
      prop_gt_1.05 = mean(rhat > 1.05),
      prop_gt_1.10 = mean(rhat > 1.10),
      .groups = "drop"
    )
}

summarise_ess_block <- function(ess_df) {
  ess_df |>
    dplyr::filter(is.finite(ess)) |>
    dplyr::group_by(block) |>
    dplyr::summarise(
      n = dplyr::n(),
      median = median(ess),
      q10 = quantile(ess, 0.10),
      q25 = quantile(ess, 0.25),
      min = min(ess),
      .groups = "drop"
    )
}

make_mcmc_review_ordprobit <- function(fit,
                                       label = "fit",
                                       out_dir = NULL,
                                       save_plots = TRUE,
                                       save_csvs = TRUE) {
  stopifnot(!is.null(fit$diagnostics))
  required <- c("dplyr", "tidyr")
  if (isTRUE(save_plots)) required <- c(required, "ggplot2")
  require_study2_reporting_packages(required, "MCMC review")
  
  rhat_parts <- list()
  
  if (!is.null(fit$diagnostics$rhat_hyper)) {
    rhat_parts$hyper <- fit$diagnostics$rhat_hyper |>
      dplyr::mutate(block = "GP/noise")
  }
  
  if (!is.null(fit$diagnostics$rhat_A)) {
    rhat_parts$A <- fit$diagnostics$rhat_A |>
      dplyr::mutate(block = "loading A")
  }
  
  if (!is.null(fit$diagnostics$rhat_tau)) {
    rhat_parts$tau <- fit$diagnostics$rhat_tau |>
      dplyr::mutate(block = "cutpoints")
  }
  
  if (!is.null(fit$diagnostics$rhat_U)) {
    rhat_parts$U <- fit$diagnostics$rhat_U |>
      dplyr::mutate(block = "latent U")
  }
  
  rhat_df <- dplyr::bind_rows(rhat_parts)
  
  if (!"block" %in% names(rhat_df)) {
    rhat_df$block <- classify_parameter_block(rhat_df$parameter)
  }
  
  rhat_summary <- summarise_rhat_block(rhat_df)
  
  ess_df <- fit$diagnostics$ess_key
  
  if (!is.null(ess_df) && nrow(ess_df) > 0L) {
    ess_df$block <- classify_parameter_block(ess_df$parameter)
    ess_summary <- summarise_ess_block(ess_df)
  } else {
    ess_summary <- data.frame()
  }
  
  chain_stats <- fit$mcmc$chain_stats
  
  if (!is.null(chain_stats) && nrow(chain_stats) > 0L) {
    chain_stats <- chain_stats |>
      dplyr::mutate(
        u_ess_accept_rate = u_ess_accept_total / pmax(u_ess_total, 1),
        global_u_accept_rate = global_u_accept_total / pmax(global_u_total, 1),
        u_eval_per_update = u_ess_eval_total / pmax(u_ess_total, 1),
        global_u_eval_per_update = global_u_eval_total / pmax(global_u_total, 1),
        theta_eval_per_update = theta_eval_total / pmax(theta_update_total, 1)
      )
  }
  
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    
    if (save_csvs) {
      utils::write.csv(
        rhat_df,
        file.path(out_dir, paste0(label, "_rhat_all.csv")),
        row.names = FALSE
      )
      
      utils::write.csv(
        rhat_summary,
        file.path(out_dir, paste0(label, "_rhat_summary_by_block.csv")),
        row.names = FALSE
      )
      
      if (!is.null(ess_df) && nrow(ess_df) > 0L) {
        utils::write.csv(
          ess_df,
          file.path(out_dir, paste0(label, "_ess_all.csv")),
          row.names = FALSE
        )
        
        utils::write.csv(
          ess_summary,
          file.path(out_dir, paste0(label, "_ess_summary_by_block.csv")),
          row.names = FALSE
        )
      }
      
      if (!is.null(chain_stats) && nrow(chain_stats) > 0L) {
        utils::write.csv(
          chain_stats,
          file.path(out_dir, paste0(label, "_chain_update_stats.csv")),
          row.names = FALSE
        )
      }
    }
  }
  
  plots <- list()
  
  if (nrow(rhat_df) > 0L) {
    if (nrow(rhat_df) > 0L) {
      rhat_plot_df <- rhat_df |>
        dplyr::filter(is.finite(rhat))
      
      rhat_ref <- data.frame(
        xintercept = c(1.01, 1.05, 1.10),
        ref = factor(
          c("1.01", "1.05", "1.10"),
          levels = c("1.01", "1.05", "1.10")
        )
      )
      
      plots$rhat_hist <- ggplot2::ggplot(
        rhat_plot_df,
        ggplot2::aes(x = rhat)
      ) +
        ggplot2::geom_histogram(
          bins = 40,
          fill = "steelblue",
          color = "white",
          alpha = 0.85
        ) +
        ggplot2::geom_vline(
          data = rhat_ref,
          ggplot2::aes(
            xintercept = xintercept,
            color = ref,
            linetype = ref
          ),
          inherit.aes = FALSE,
          linewidth = 0.6
        ) +
        ggplot2::scale_color_manual(
          values = c(
            "1.01" = "grey30",
            "1.05" = "orange3",
            "1.10" = "red3"
          ),
          name = "Reference"
        ) +
        ggplot2::scale_linetype_manual(
          values = c(
            "1.01" = "dotted",
            "1.05" = "dashed",
            "1.10" = "solid"
          ),
          name = "Reference"
        ) +
        ggplot2::facet_wrap(~block, scales = "free_y") +
        ggplot2::labs(
          x = expression(hat(R)),
          y = "Count",
          title = paste0("MCMC review: R-hat distribution by block, ", label),
          subtitle = "Vertical lines: 1.01, 1.05, 1.10"
        ) +
        ggplot2::theme(legend.position = "bottom")
    }
  }
  
  if (!is.null(ess_df) && nrow(ess_df) > 0L) {
    plots$ess_box <- ggplot2::ggplot(
      ess_df |> dplyr::filter(is.finite(ess), ess > 0),
      ggplot2::aes(x = block, y = ess, fill = block)
    ) +
      ggplot2::geom_boxplot(alpha = 0.80, outlier.alpha = 0.45) +
      ggplot2::scale_y_log10() +
      ggplot2::labs(
        x = NULL,
        y = "ESS, log scale",
        title = paste0("MCMC review: ESS distribution by block, ", label)
      ) +
      ggplot2::theme(legend.position = "none")
  }
  
  if (!is.null(chain_stats) && nrow(chain_stats) > 0L) {
    
    chain_stat_cols <- c(
      "u_ess_accept_rate",
      "global_u_accept_rate",
      "u_eval_per_update",
      "global_u_eval_per_update",
      "theta_eval_per_update"
    )
    
    chain_stat_cols <- intersect(chain_stat_cols, names(chain_stats))
    
    chain_long <- chain_stats |>
      dplyr::select(
        dplyr::all_of(c("chain", chain_stat_cols))
      ) |>
      tidyr::pivot_longer(
        cols = -chain,
        names_to = "statistic",
        values_to = "value"
      )
    
    plots$chain_update_stats <- ggplot2::ggplot(
      chain_long, ggplot2::aes(x = factor(chain), y = value)
    ) +
      ggplot2::geom_col(fill = "steelblue", alpha = 0.85) +
      ggplot2::facet_wrap(~statistic, scales = "free_y", ncol = 3) +
      ggplot2::labs(
        x = "Chain",
        y = NULL,
        title = paste0("MCMC review: chain update statistics, ", label)
      )
  }
  
  if (!is.null(out_dir) && save_plots) {
    if (!is.null(plots$rhat_hist)) {
      ggplot2::ggsave(
        file.path(out_dir, paste0(label, "_rhat_hist_by_block.pdf")),
        plots$rhat_hist,
        width = 10,
        height = 6
      )
    }
    
    if (!is.null(plots$ess_box)) {
      ggplot2::ggsave(
        file.path(out_dir, paste0(label, "_ess_box_by_block.pdf")),
        plots$ess_box,
        width = 8,
        height = 5
      )
    }
    
    if (!is.null(plots$chain_update_stats)) {
      ggplot2::ggsave(
        file.path(out_dir, paste0(label, "_chain_update_stats.pdf")),
        plots$chain_update_stats,
        width = 10,
        height = 6
      )
    }
  }
  
  list(
    rhat_df = rhat_df,
    rhat_summary = rhat_summary,
    ess_df = ess_df,
    ess_summary = ess_summary,
    chain_stats = chain_stats,
    plots = plots
  )
}
