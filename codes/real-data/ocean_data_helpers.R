############################################################
## Ocean data helpers (not part of 00_study1_functions.R).
## 01 / 02 source 00_study1_functions.R first, then this file.
############################################################

ocean_current_script_dir <- function(fallback = getwd()) {
  frames <- sys.frames()
  for (ii in rev(seq_along(frames))) {
    ofile <- frames[[ii]]$ofile
    if (!is.null(ofile) && nzchar(ofile)) {
      return(dirname(normalizePath(ofile, mustWork = FALSE)))
    }
  }
  normalizePath(fallback, mustWork = FALSE)
}

if (!exists("OCEAN_REALDATA_DIR")) {
  OCEAN_REALDATA_DIR <- ocean_current_script_dir()
}

default_ocean_profile_dir <- function() {
  candidates <- c(
    Sys.getenv("OCEAN_DATA_DIR", unset = ""),
    if (exists("OCEAN_DATA_DIR")) OCEAN_DATA_DIR else "",
    file.path(OCEAN_REALDATA_DIR, "data", "prepared"),
    file.path(
      OCEAN_REALDATA_DIR, "..", "..",
      "code", "New-version", "ocean_new",
      "bcodmo_np_role_swap_class6_split_2026-08-05",
      "data", "derived", "positive_primary_np_swap_class6"
    ),
    file.path(
      "code", "New-version", "ocean_new",
      "bcodmo_np_role_swap_class6_split_2026-08-05",
      "data", "derived", "positive_primary_np_swap_class6"
    )
  )
  candidates <- candidates[nzchar(candidates)]
  for (p in candidates) {
    if (file.exists(file.path(p, "study1_data.csv"))) {
      return(normalizePath(p, mustWork = TRUE))
    }
  }
  stop(
    "Could not find ocean study1_data.csv. Set OCEAN_DATA_DIR to the ",
    "positive_primary_np_swap_class6 folder."
  )
}

load_ocean_prepared <- function(data_dir = NULL) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to read ocean meta.json")
  }
  if (is.null(data_dir)) data_dir <- default_ocean_profile_dir()
  data_dir <- normalizePath(data_dir, mustWork = TRUE)
  dat <- read.csv(file.path(data_dir, "study1_data.csv"), stringsAsFactors = FALSE)
  meta <- jsonlite::fromJSON(file.path(data_dir, "meta.json"), simplifyVector = TRUE)
  x_cols <- as.character(meta$x_columns)
  required <- c("split", "y", "u_log_std", "c", x_cols)
  missing_columns <- setdiff(required, names(dat))
  if (length(missing_columns) > 0L) {
    stop("Ocean data are missing column(s): ", paste(missing_columns, collapse = ", "))
  }
  if (length(x_cols) < 1L || anyDuplicated(x_cols)) {
    stop("meta.json must specify at least one distinct x column.")
  }
  numeric_columns <- c("y", "u_log_std", "c", x_cols)
  if (any(!vapply(dat[numeric_columns], is.numeric, logical(1)))) {
    stop("Ocean y, u, c, and x columns must be numeric.")
  }
  if (any(!is.finite(as.matrix(dat[numeric_columns])))) {
    stop("Ocean y, u, c, and x columns must be complete and finite.")
  }
  train <- dat[dat$split == "train", , drop = FALSE]
  test <- dat[dat$split == "test", , drop = FALSE]
  if (nrow(train) < 2L || nrow(test) < 1L) {
    stop("Ocean data must contain nonempty train and test splits.")
  }
  m <- as.integer(meta$m)
  if (length(m) != 1L || is.na(m) || m < 2L ||
      any(dat$c != as.integer(dat$c)) || any(dat$c < 1L | dat$c > m)) {
    stop("Ocean ordinal codes must be integers 1, ..., meta$m.")
  }
  tau_reference <- as.numeric(meta$log_tau_reference)
  if (length(tau_reference) != m - 1L ||
      any(!is.finite(tau_reference)) ||
      is.unsorted(tau_reference, strictly = TRUE)) {
    stop("meta$log_tau_reference must contain m - 1 increasing thresholds.")
  }
  list(
    data_dir = data_dir,
    meta = meta,
    train = train,
    test = test,
    x_cols = x_cols,
    X_train = as.matrix(train[, x_cols, drop = FALSE]),
    X_test = as.matrix(test[, x_cols, drop = FALSE]),
    u_train = train$u_log_std,
    u_test = test$u_log_std,
    y_train = train$y,
    y_test = test$y,
    c_train = train$c,
    c_test = test$c,
    m = m,
    tau_reference = tau_reference,
    log_center = as.numeric(meta$log_center),
    log_scale = as.numeric(meta$log_scale),
    inverse_u = function(z) exp(meta$log_center + meta$log_scale * z)
  )
}

read_ocean_class_allocation <- function(path = NULL) {
  if (is.null(path)) {
    local <- file.path(
      OCEAN_REALDATA_DIR, "data", "class6_proportional_calibration_allocation.csv"
    )
    parent <- file.path(
      default_ocean_profile_dir(), "..", "..", "calibration_designs",
      "class6_proportional_calibration_allocation.csv"
    )
    path <- if (file.exists(local)) local else parent
  }
  path <- normalizePath(path, mustWork = TRUE)
  out <- read.csv(path, stringsAsFactors = FALSE)
  if (!"c" %in% names(out) || anyDuplicated(out$c)) {
    stop("Calibration allocation must contain one distinct row per class c.")
  }
  out
}

## Nested within-class proportional seats: 10 ⊂ 25 ⊂ 50.
make_ocean_nested_calibration <- function(c_train,
                                          allocation,
                                          calib_grid = c(10L, 25L, 50L),
                                          seed = 20260805L) {
  c_train <- as.integer(c_train)
  calib_grid <- sort(unique(as.integer(calib_grid)))
  if (length(c_train) < 1L || anyNA(c_train) ||
      length(calib_grid) < 1L || anyNA(calib_grid) || any(calib_grid < 1L)) {
    stop("c_train and calib_grid must contain valid positive integer values.")
  }
  m <- nrow(allocation)
  if (!setequal(allocation$c, seq_len(m))) {
    stop("Allocation rows must represent ordinal classes 1, ..., m.")
  }
  allocation <- allocation[match(seq_len(m), allocation$c), , drop = FALSE]
  if ("train_n" %in% names(allocation)) {
    observed_counts <- tabulate(c_train, nbins = m)
    if (!identical(as.integer(allocation$train_n), as.integer(observed_counts))) {
      stop("Allocation train_n does not match the observed training class counts.")
    }
  }
  allocation_columns <- paste0("n_calib_", calib_grid)
  if (any(!allocation_columns %in% names(allocation))) {
    stop(
      "Allocation is missing column(s): ",
      paste(setdiff(allocation_columns, names(allocation)), collapse = ", ")
    )
  }
  allocation_counts <- as.matrix(allocation[, allocation_columns, drop = FALSE])
  if (anyNA(allocation_counts) || any(allocation_counts < 0) ||
      any(allocation_counts != floor(allocation_counts))) {
    stop("Calibration allocations must be nonnegative integers.")
  }
  if (ncol(allocation_counts) > 1L &&
      any(apply(allocation_counts, 1L, function(z) any(diff(z) < 0)))) {
    stop("Within-class allocations must be nondecreasing over calib_grid.")
  }
  out <- setNames(vector("list", length(calib_grid)), as.character(calib_grid))
  for (nm in names(out)) out[[nm]] <- integer(0)

  for (cc in seq_len(m)) {
    class_idx <- which(c_train == allocation$c[cc])
    set.seed(as.integer(seed) + cc)
    ordering <- sample(class_idx, length(class_idx), replace = FALSE)
    for (k in calib_grid) {
      col <- paste0("n_calib_", k)
      if (!col %in% names(allocation)) {
        stop("Allocation is missing column ", col)
      }
      n_take <- as.integer(allocation[[col]][cc])
      if (n_take > length(ordering)) {
        stop("Class ", cc, " has fewer observations than requested calibrations.")
      }
      if (n_take > 0L) {
        out[[as.character(k)]] <- c(
          out[[as.character(k)]],
          ordering[seq_len(n_take)]
        )
      }
    }
  }
  for (nm in names(out)) {
    out[[nm]] <- sort(as.integer(out[[nm]]))
    expected <- sum(allocation[[paste0("n_calib_", nm)]])
    if (length(out[[nm]]) != expected) {
      stop("Calibration |O|=", nm, " has ", length(out[[nm]]),
           " seats, expected ", expected)
    }
  }
  out
}

summarize_ocean_u_imputation <- function(fit, u_true, c_ord, inverse_u = NULL) {
  miss <- fit$data$miss_idx
  samples <- fit$mcmc$samples_u[, miss, drop = FALSE]
  truth <- u_true[miss]
  post_mean <- colMeans(samples)
  qs <- apply(samples, 2L, quantile, probs = c(0.025, 0.975))
  ans <- data.frame(
    train_row_id = miss,
    c = c_ord[miss],
    true_u = truth,
    post_mean = post_mean,
    post_lo = qs[1L, ],
    post_hi = qs[2L, ],
    covered95 = truth >= qs[1L, ] & truth <= qs[2L, ]
  )
  if (!is.null(inverse_u)) {
    raw_samples <- inverse_u(samples)
    raw_truth <- inverse_u(truth)
    raw_mean <- colMeans(raw_samples)
    raw_qs <- apply(raw_samples, 2L, quantile, probs = c(0.025, 0.975))
    ans$true_u_raw <- raw_truth
    ans$post_mean_raw <- raw_mean
    ans$post_lo_raw <- raw_qs[1L, ]
    ans$post_hi_raw <- raw_qs[2L, ]
    ans$covered95_raw <- raw_truth >= raw_qs[1L, ] & raw_truth <= raw_qs[2L, ]
  }
  ans
}

prepare_gp_mle_prediction_1d <- function(fit) {
  X <- fit$X
  y <- fit$y
  n <- nrow(X)
  d <- ncol(X)
  sigma2 <- exp(fit$par[1])
  rho <- exp(fit$par[2])
  theta <- exp(fit$par[-c(1, 2)])
  kernel_spec <- kernel_spec_from_fit_1d(fit)

  D2 <- matrix(0, n, n)
  for (j in seq_len(d)) {
    D2 <- D2 + theta[j] * pairwise_sqdist(X[, j, drop = FALSE])
  }
  R <- kernel_from_weighted_sqdist_1d(
    D2, kernel_spec$name, kernel_spec$matern_nu
  )
  K <- rho^2 * sigma2 * R
  chol_train <- safe_chol_1d(K + sigma2 * diag(n))

  list(
    X = X,
    alpha = solve_chol(chol_train, y),
    chol_train = chol_train,
    sigma2 = sigma2,
    rho = rho,
    theta = theta,
    kernel = kernel_spec
  )
}

estimate_tau_from_ordinal_codes <- function(c_ord, m, tail_clip = 1e-4) {
  c_ord <- as.integer(c_ord)
  m <- as.integer(m)
  if (length(m) != 1L || is.na(m) || m < 2L ||
      length(c_ord) < 1L || anyNA(c_ord) || any(c_ord < 1L | c_ord > m)) {
    stop("c_ord must contain valid codes 1, ..., m.")
  }
  probs <- cumsum(tabulate(c_ord, nbins = m))[seq_len(m - 1L)] / length(c_ord)
  probs <- pmin(pmax(probs, tail_clip), 1 - tail_clip)
  tau <- qnorm(probs)
  if (is.unsorted(tau, strictly = TRUE)) {
    stop("All ordinal levels must occur to estimate distinct cutpoints.")
  }
  tau
}

gp_mle_predict_prepared_1d <- function(prepared, Xstar, noisy = FALSE) {
  Xstar <- as_study1_numeric_matrix(
    Xstar,
    "Xstar",
    ncol_expected = ncol(prepared$X),
    expected_names = colnames(prepared$X)
  )
  D2_star <- matrix(0, nrow(Xstar), nrow(prepared$X))
  for (j in seq_along(prepared$theta)) {
    D2_star <- D2_star + prepared$theta[j] * pairwise_sqdist(
      Xstar[, j, drop = FALSE],
      prepared$X[, j, drop = FALSE]
    )
  }
  Rstar <- kernel_from_weighted_sqdist_1d(
    D2_star, prepared$kernel$name, prepared$kernel$matern_nu
  )
  Kstar <- prepared$rho^2 * prepared$sigma2 * Rstar
  mean <- as.numeric(Kstar %*% prepared$alpha)
  v <- forwardsolve(t(prepared$chol_train), t(Kstar))
  prior_var <- prepared$rho^2 * prepared$sigma2
  variance_reduction <- colSums(v^2)
  var <- prior_var - variance_reduction
  tolerance <- 100 * (nrow(prepared$X) + nrow(Xstar) + 1L) *
    .Machine$double.eps *
    max(1, abs(prior_var), max(abs(variance_reduction)))
  if (any(!is.finite(var)) || any(var < -tolerance)) {
    stop(
      "Ocean GP prediction produced a materially negative or non-finite ",
      "latent predictive variance."
    )
  }
  var <- pmax(var, 0)
  if (isTRUE(noisy)) var <- var + prepared$sigma2
  list(mean = mean, var = var)
}

fit_ocean_latent_gp <- function(X,
                                y,
                                u,
                                train_idx = seq_len(nrow(X)),
                                kernel = "se",
                                matern_nu = 2.5) {
  X <- as_study1_numeric_matrix(X, "X", nrow_expected = length(y))
  u <- as.numeric(u)
  y <- as.numeric(y)
  train_idx <- sort(unique(as.integer(train_idx)))
  if (length(u) != nrow(X) || any(!is.finite(u[train_idx])) ||
      any(!is.finite(y)) || length(train_idx) < 3L ||
      any(train_idx < 1L | train_idx > nrow(X))) {
    stop("Latent-input GP requires at least three valid indexed observations.")
  }

  X_center <- colMeans(X)
  X_scale <- apply(X, 2L, sd)
  if (any(!is.finite(X_scale) | X_scale <= 0)) {
    stop("Every ocean x column must have positive finite variation.")
  }
  X_std <- sweep(sweep(X, 2L, X_center, "-"), 2L, X_scale, "/")
  y_center <- mean(y[train_idx])
  y_scale <- sd(y[train_idx])
  if (!is.finite(y_scale) || y_scale <= 0) {
    stop("Complete-case y values must have positive finite variation.")
  }
  X_aug <- cbind(X_std[train_idx, , drop = FALSE], u = u[train_idx])

  list(
    fit = gp_mle_fit_1d(
      X_aug,
      (y[train_idx] - y_center) / y_scale,
      kernel = kernel,
      matern_nu = matern_nu
    ),
    X_center = X_center,
    X_scale = X_scale,
    y_center = y_center,
    y_scale = y_scale,
    train_idx = train_idx
  )
}

sample_ocean_latent_gp_y_given_xc <- function(gp_fit,
                                               X_test,
                                               c_test,
                                               tau,
                                               n_draw = 500L,
                                               seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  predictor_names <- colnames(gp_fit$fit$X)[seq_along(gp_fit$X_center)]
  X_test <- as_study1_numeric_matrix(
    X_test,
    "X_test",
    ncol_expected = length(gp_fit$X_center),
    expected_names = predictor_names
  )
  c_test <- as.integer(c_test)
  m <- length(tau) + 1L
  if (length(c_test) != nrow(X_test) || anyNA(c_test) ||
      any(c_test < 1L | c_test > m)) {
    stop("c_test must contain one valid ordinal code per test row.")
  }
  n_draw <- as.integer(n_draw)
  if (length(n_draw) != 1L || is.na(n_draw) || n_draw < 1L) {
    stop("n_draw must be one positive integer.")
  }

  X_std <- sweep(
    sweep(X_test, 2L, gp_fit$X_center, "-"),
    2L, gp_fit$X_scale, "/"
  )
  lower <- c(-Inf, tau)[c_test]
  upper <- c(tau, Inf)[c_test]
  prepared <- prepare_gp_mle_prediction_1d(gp_fit$fit)
  out <- matrix(NA_real_, n_draw, nrow(X_test))

  for (ii in seq_len(n_draw)) {
    u_star <- rtruncnorm_vec(
      mean = rep(0, nrow(X_test)),
      sd = 1,
      lower = lower,
      upper = upper
    )
    X_aug <- cbind(X_std, u = u_star)
    colnames(X_aug) <- colnames(gp_fit$fit$X)
    pred <- gp_mle_predict_prepared_1d(prepared, X_aug, noisy = TRUE)
    y_std <- pred$mean + sqrt(pred$var) * rnorm(nrow(X_test))
    out[ii, ] <- gp_fit$y_center + gp_fit$y_scale * y_std
  }
  out
}
