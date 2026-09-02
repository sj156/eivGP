############################################################
## 03_study2_published_competitors.R
##
## Common, auditable competitor layer for Studies I and II.
## The historical filename is retained so older entry points still work.
##
## Publication rules:
##   * use only the public implementation attached to the cited method;
##   * never relabel a home-built approximation as a published method;
##   * record unavailable and failed fits explicitly;
##   * return predictive draws for the noisy-response target Y* | X*, C*;
##   * retain the package's fitted latent/predictive mean separately, so that
##     recovery of m(x,c) is not scored by averaging Monte Carlo-noisy Y draws.
############################################################

`%||%` <- function(x, y) if (is.null(x)) y else x

mixedgp_competitor_registry <- function() {
  data.frame(
    method = c("UC-GP", "LVGP", "EzGP"),
    reference = c(
      "Qian et al. (2008); Zhou et al. (2011)",
      "Zhang et al. (2020)",
      "Xiao et al. (2021)"
    ),
    package = c("kergp", "LVGP", "EzGP"),
    implementation = c(
      "kergp::q1Symm",
      "LVGP::LVGP_fit and LVGP::LVGP_predict",
      "EzGP::EzGP_fit and EzGP::EzGP_predict"
    ),
    stringsAsFactors = FALSE
  )
}

## Backward-compatible name used by existing Study II scripts.
study2_competitor_registry <- mixedgp_competitor_registry

mixedgp_competitor_preflight <- function(
    methods = mixedgp_competitor_registry()$method,
    strict = FALSE) {
  registry <- mixedgp_competitor_registry()
  unknown <- setdiff(methods, registry$method)
  if (length(unknown) > 0L) {
    stop("Unknown published competitor(s): ", paste(unknown, collapse = ", "))
  }

  out <- registry[match(methods, registry$method), , drop = FALSE]
  out$available <- vapply(
    out$package,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
  out$version <- vapply(seq_len(nrow(out)), function(i) {
    if (!out$available[i]) return(NA_character_)
    as.character(utils::packageVersion(out$package[i]))
  }, character(1))

  if (isTRUE(strict) && any(!out$available)) {
    stop(
      "Publication preflight failed. Install the missing public package(s): ",
      paste(out$package[!out$available], collapse = ", "),
      ". No substitute implementation will be used."
    )
  }
  out
}

mixedgp_competitor_status <- function(method,
                                      status,
                                      optimization_status = "",
                                      optimization_attempts = "",
                                      message = "",
                                      warnings = "",
                                      elapsed = NA_real_,
                                      implementation = "",
                                      fitted_noise_variance = NA_real_) {
  data.frame(
    method = method,
    status = status,
    optimization_status = optimization_status,
    optimization_attempts = optimization_attempts,
    message = message,
    warnings = warnings,
    elapsed_seconds = elapsed,
    implementation = implementation,
    fitted_noise_variance = fitted_noise_variance,
    stringsAsFactors = FALSE
  )
}

study2_competitor_status <- mixedgp_competitor_status

format_mixedgp_optimizer_attempts <- function(attempt) {
  if (is.null(attempt) || nrow(attempt) == 0L) return("")
  fmt_num <- function(x, digits = 6L) {
    ifelse(
      is.finite(x),
      trimws(formatC(x, digits = digits, format = "fg")),
      "NA"
    )
  }
  paste(
    paste0(
      "retry=", attempt$retry,
      ",seed=", attempt$seed,
      ",limits=", attempt$max_iter_ini, "/", attempt$max_iter_lat,
      ",elapsed=", fmt_num(attempt$elapsed_seconds, 4L),
      ",status=", attempt$status,
      ",flag=", ifelse(is.na(attempt$selected_flag), "NA", attempt$selected_flag),
      ",objective=", fmt_num(attempt$reported_objective)
    ),
    collapse = " | "
  )
}

validate_mixedgp_inputs <- function(X, C, m_vec = NULL) {
  X <- as.matrix(X)
  C <- as.matrix(C)
  storage.mode(X) <- "double"
  storage.mode(C) <- "integer"

  if (nrow(X) != nrow(C) || ncol(X) < 1L || ncol(C) < 1L) {
    stop("X and C must have the same number of rows and at least one column.")
  }
  if (any(!is.finite(X)) || any(!is.finite(C))) {
    stop("Competitor inputs must be finite.")
  }

  if (is.null(m_vec)) m_vec <- apply(C, 2L, max)
  m_vec <- as.integer(m_vec)
  if (length(m_vec) != ncol(C) || any(m_vec < 2L)) {
    stop("m_vec must give at least two ordered levels for every column of C.")
  }
  for (j in seq_len(ncol(C))) {
    if (any(C[, j] < 1L | C[, j] > m_vec[j])) {
      stop("Column ", j, " of C contains a value outside 1:m_vec[j].")
    }
  }

  list(X = X, C = C, m_vec = m_vec)
}

make_mixedgp_competitor_data <- function(X, C, m_vec, y = NULL) {
  checked <- validate_mixedgp_inputs(X, C, m_vec)
  X <- checked$X
  C <- checked$C
  m_vec <- checked$m_vec

  out <- as.data.frame(X)
  names(out) <- paste0("x", seq_len(ncol(X)))
  for (j in seq_len(ncol(C))) {
    out[[paste0("c", j)]] <- factor(C[, j], levels = seq_len(m_vec[j]))
  }
  if (!is.null(y)) {
    if (length(y) != nrow(out) || any(!is.finite(y))) {
      stop("y must be finite and have one value per training row.")
    }
    out$y <- as.numeric(y)
  }
  out
}

make_mixedgp_numeric_matrix <- function(X, C, m_vec) {
  checked <- validate_mixedgp_inputs(X, C, m_vec)
  out <- cbind(checked$X, checked$C)
  storage.mode(out) <- "double"
  colnames(out) <- c(
    paste0("x", seq_len(ncol(checked$X))),
    paste0("c", seq_len(ncol(checked$C)))
  )
  out
}

predictive_variance_diagonal <- function(x,
                                         n_test,
                                         context = "package prediction") {
  if (is.null(x)) stop("The package prediction did not return uncertainty.")
  if (is.matrix(x)) {
    if (!identical(dim(x), c(as.integer(n_test), as.integer(n_test)))) {
      stop("Unexpected predictive covariance dimension.")
    }
    x <- diag(x)
  }
  x <- as.numeric(x)
  if (length(x) != n_test || any(!is.finite(x))) {
    stop("Predictive variances must be finite with length n_test.")
  }
  tolerance <- 100 * (as.integer(n_test) + 1L) * .Machine$double.eps *
    max(1, max(abs(x)))
  if (any(x < -tolerance)) {
    stop(context, " returned a materially negative predictive variance.")
  }
  pmax(x, 0)
}

sample_independent_predictive_marginals <- function(mean, variance, n_draw, seed) {
  mean <- as.numeric(mean)
  variance <- as.numeric(variance)
  if (length(mean) != length(variance) || any(!is.finite(mean)) ||
      any(!is.finite(variance)) || any(variance < 0)) {
    stop("Invalid predictive means or variances.")
  }
  set.seed(seed)
  matrix(
    stats::rnorm(
      n_draw * length(mean),
      mean = rep(mean, each = n_draw),
      sd = rep(sqrt(variance), each = n_draw)
    ),
    nrow = n_draw,
    ncol = length(mean)
  )
}

attach_predictive_normal_components <- function(draws, mean, variance) {
  draws <- as.matrix(draws)
  mean <- as.numeric(mean)
  variance <- as.numeric(variance)
  if (length(mean) != ncol(draws) || length(variance) != ncol(draws) ||
      any(!is.finite(mean)) || any(!is.finite(variance)) ||
      any(variance < 0)) {
    stop(
      "A Gaussian predictive distribution requires one finite mean and one ",
      "nonnegative variance per test case."
    )
  }
  ## Each package adapter supplies one fitted Gaussian predictive marginal per
  ## test case.  A single row is therefore the exact mixture representation;
  ## it must not be repeated n_draw times when computing log density.
  attr(draws, "conditional_means") <- matrix(mean, nrow = 1L)
  attr(draws, "conditional_vars") <- matrix(variance, nrow = 1L)
  attr(draws, "mixture_components") <- "public-package Gaussian predictor"
  draws
}

mixedgp_adapter_ucgp <- function(X_train,
                                 y_train,
                                 C_train,
                                 X_test,
                                 C_test,
                                 m_vec,
                                 n_draw,
                                 seed,
                                 n_starts = 8L) {
  if (!requireNamespace("kergp", quietly = TRUE)) {
    stop("Package 'kergp' is not installed.")
  }
  checked <- validate_mixedgp_inputs(X_train, C_train, m_vec)
  p <- ncol(checked$X)
  q <- ncol(checked$C)
  train <- make_mixedgp_competitor_data(
    checked$X, checked$C, checked$m_vec, y_train
  )
  test <- make_mixedgp_competitor_data(X_test, C_test, checked$m_vec)

  kernel_env <- new.env(parent = parent.frame())
  kernel_env$kx <- kergp::kGauss(d = p)
  if (!identical(kergp::inputNames(kernel_env$kx), paste0("x", seq_len(p)))) {
    stop("Unexpected quantitative-input names returned by kergp::kGauss().")
  }
  for (j in seq_len(q)) {
    kernel_env[[paste0("kc", j)]] <- kergp::q1Symm(
      factor = factor(seq_len(checked$m_vec[j]),
                      levels = seq_len(checked$m_vec[j])),
      input = paste0("c", j),
      cov = "corr"
    )
  }
  kernel_formula <- stats::as.formula(
    paste0(
      "~ ",
      paste(c("kx()", paste0("kc", seq_len(q), "()")), collapse = " * ")
    ),
    env = kernel_env
  )
  covariance <- kergp::covComp(formula = kernel_formula, where = kernel_env)

  set.seed(seed)
  multistart_fit <- kergp::gp(
    y ~ 1,
    data = train,
    cov = covariance,
    noise = TRUE,
    multistart = as.integer(n_starts),
    trace = 0
  )
  report <- multistart_fit$MLE$report
  finite_rows <- apply(report$par, 1L, function(z) all(is.finite(z)))
  eligible <- which(report$convergence & is.finite(report$logLik) & finite_rows)
  if (length(eligible) == 0L) {
    stop("No kergp likelihood start converged to a finite fit.")
  }
  selected <- eligible[which.max(report$logLik[eligible])]

  ## kergp's multistart routine selects the largest likelihood even when that
  ## start has a nonconvergence flag. Reconstruct the public-package GP object
  ## at the best *converged* reported solution instead of accepting that flag
  ## or discarding converged alternatives.
  selected_par <- as.numeric(report$par[selected, ])
  covariance_par <- methods::getMethod("coef", "covComp")(covariance)
  n_covariance_par <- length(covariance_par)
  covariance_selected <- kergp::`coef<-`(
    covariance,
    value = selected_par[seq_len(n_covariance_par)]
  )
  noise_var <- selected_par[n_covariance_par + 1L]
  fit <- kergp::gp(
    y ~ 1,
    data = train,
    cov = covariance_selected,
    estim = FALSE,
    varNoise = noise_var
  )
  fit$MLE <- multistart_fit$MLE
  fit$logLik <- as.numeric(report$logLik[selected])
  fit$selected_multistart_index <- selected
  fit$selected_multistart_report <- data.frame(
    start = seq_len(nrow(report$par)),
    logLik = as.numeric(report$logLik),
    convergence = as.logical(report$convergence)
  )
  pred <- stats::predict(
    fit,
    newdata = test,
    type = "UK",
    seCompute = TRUE,
    covCompute = FALSE,
    lightReturn = TRUE
  )
  noise_var <- as.numeric(fit$varNoise)
  if (length(noise_var) != 1L || !is.finite(noise_var) || noise_var < 0) {
    stop("kergp returned an invalid fitted noise variance.")
  }
  pred_var_y <- predictive_variance_diagonal(
    as.numeric(pred$sd)^2 + noise_var,
    nrow(test),
    context = "kergp"
  )
  latent_mean <- as.numeric(pred$mean)

  predictive_draws <- sample_independent_predictive_marginals(
    latent_mean, pred_var_y, n_draw, seed + 10000L
  )
  predictive_draws <- attach_predictive_normal_components(
    predictive_draws, latent_mean, pred_var_y
  )

  list(
    draws = predictive_draws,
    ## The fitted response nugget is zero mean, so the latent-function and
    ## noisy-response predictive means coincide.  Keep both names explicit:
    ## the former is used for m(x,c), whereas draws above target Y*.
    latent_mean = latent_mean,
    predictive_mean = latent_mean,
    predictive_variance = pred_var_y,
    fit = fit,
    fitted_noise_variance = noise_var,
    implementation = paste0(
      "kergp ", utils::packageVersion("kergp"),
      "::q1Symm; best of ", sum(report$convergence),
      " converged fits among ", n_starts, " likelihood starts"
    )
  )
}

mixedgp_adapter_lvgp <- function(X_train,
                                 y_train,
                                 C_train,
                                 X_test,
                                 C_test,
                                 m_vec,
                                 n_draw,
                                 seed,
                                 n_starts = 8L,
                                 max_retries = 3L,
                                 dim_z = 2L,
                                 max_iter_ini = 100L,
                                 max_iter_lat = 20L,
                                 rescue_iter_ini = 300L,
                                 rescue_iter_lat = 100L,
                                 max_elapsed_seconds = Inf,
                                 parallel = FALSE) {
  if (!requireNamespace("LVGP", quietly = TRUE)) {
    stop("Package 'LVGP' is not installed.")
  }
  checked <- validate_mixedgp_inputs(X_train, C_train, m_vec)
  p <- ncol(checked$X)
  q <- ncol(checked$C)
  train_x <- make_mixedgp_numeric_matrix(
    checked$X, checked$C, checked$m_vec
  )
  test_x <- make_mixedgp_numeric_matrix(X_test, C_test, checked$m_vec)

  ## LVGP_fit() selects the smallest objective without requiring its optim
  ## convergence code to be zero. Because the returned fit object is tied to
  ## that selected solution, do not reconstruct an undocumented alternative by
  ## hand. Use a prespecified sequence of complete author-package multistart
  ## fits, prefer a converged return, and—only if every rescue still reports a
  ## nonzero code—retain the best finite package-selected fit with an explicit
  ## optimizer warning. This follows the public package more faithfully than
  ## discarding every fit that the authors' prediction routine accepts.
  n_starts <- as.integer(n_starts)
  max_retries <- as.integer(max_retries)
  if (n_starts < 1L || max_retries < 1L) {
    stop("n_starts and max_retries must be positive integers.")
  }
  if (length(max_elapsed_seconds) != 1L ||
      is.na(max_elapsed_seconds) || max_elapsed_seconds <= 0) {
    stop("max_elapsed_seconds must be a positive number or Inf.")
  }
  candidate_fits <- vector("list", max_retries)
  attempt <- data.frame(
    retry = seq_len(max_retries),
    seed = as.integer(seed) + 1000L * (seq_len(max_retries) - 1L),
    max_iter_ini = c(
      as.integer(max_iter_ini),
      rep(as.integer(rescue_iter_ini), max_retries - 1L)
    ),
    max_iter_lat = c(
      as.integer(max_iter_lat),
      rep(as.integer(rescue_iter_lat), max_retries - 1L)
    ),
    elapsed_seconds = NA_real_, selected_flag = NA_integer_, converged = FALSE,
    status = "not_run",
    reported_objective = NA_real_, selection_objective = Inf, message = "",
    stringsAsFactors = FALSE
  )
  for (ss in seq_len(max_retries)) {
    attempt_start <- proc.time()[["elapsed"]]
    candidate <- tryCatch(
      local({
        if (is.finite(max_elapsed_seconds)) {
          setTimeLimit(
            elapsed = as.numeric(max_elapsed_seconds), transient = TRUE
          )
          on.exit(
            setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE),
            add = TRUE
          )
        }
        LVGP::LVGP_fit(
          X = train_x,
          Y = as.numeric(y_train),
          ind_qual = p + seq_len(q),
          dim_z = as.integer(dim_z),
          n_opt = n_starts,
          max_iter_ini = attempt$max_iter_ini[ss],
          max_iter_lat = attempt$max_iter_lat[ss],
          seed = as.integer(attempt$seed[ss]),
          progress = FALSE,
          parallel = isTRUE(parallel),
          noise = TRUE
        )
      }),
      error = function(e) e
    )
    attempt$elapsed_seconds[ss] <-
      proc.time()[["elapsed"]] - attempt_start
    if (inherits(candidate, "error")) {
      attempt$message[ss] <- conditionMessage(candidate)
      attempt$status[ss] <- if (grepl(
        "elapsed time limit|reached elapsed", attempt$message[ss],
        ignore.case = TRUE
      )) "timeout" else "error"
      next
    }
    best_try <- which.min(unlist(candidate$optim_hist$obj_sol_best))
    best_flag <- as.numeric(
      candidate$optim_hist$flag_sol_unique[[best_try]][1L]
    )
    finite_package_fit <- is.finite(candidate$fit_detail$min_n_log_l) &&
      is.finite(candidate$fit_detail$raw_min_eig)
    valid <- identical(best_flag, 0) && finite_package_fit
    attempt$selected_flag[ss] <- best_flag
    attempt$converged[ss] <- valid
    attempt$status[ss] <- if (valid) "converged" else "nonconverged"
    attempt$reported_objective[ss] <- candidate$fit_detail$min_n_log_l
    attempt$selection_objective[ss] <- if (finite_package_fit) {
      candidate$fit_detail$min_n_log_l
    } else {
      Inf
    }
    if (!valid) {
      attempt$message[ss] <- paste0(
        "Package-selected optim convergence code ", best_flag,
        "; extended-limit rescue attempted when available."
      )
    }
    if (finite_package_fit) candidate_fits[[ss]] <- candidate
    if (valid) break
  }
  finite_attempt <- is.finite(attempt$selection_objective)
  if (!any(finite_attempt)) {
    attempted <- attempt[!is.na(attempt$selected_flag), , drop = FALSE]
    detail <- if (nrow(attempted) == 0L) {
      paste(unique(attempt$message[nzchar(attempt$message)]), collapse = " | ")
    } else {
      paste0(
        "seed=", attempted$seed,
        ", limits=", attempted$max_iter_ini, "/", attempted$max_iter_lat,
        ", flag=", attempted$selected_flag,
        collapse = "; "
      )
    }
    failure <- simpleError(paste0(
      "LVGP did not return a converged finite-likelihood fit after the ",
      "prespecified default and extended-limit attempts. ", detail
    ))
    failure$optimizer_attempts <- attempt
    stop(failure)
  }
  converged_attempt <- which(attempt$converged)
  best_retry <- if (length(converged_attempt) > 0L) {
    converged_attempt[which.min(attempt$selection_objective[converged_attempt])]
  } else {
    which.min(attempt$selection_objective)
  }
  fit <- candidate_fits[[best_retry]]
  accepted_package_nonconvergence <- !isTRUE(attempt$converged[best_retry])
  if (accepted_package_nonconvergence) {
    warning(
      "LVGP returned a finite author-package-selected solution after the ",
      "prespecified attempt sequence, but optim reported code ",
      attempt$selected_flag[best_retry],
      ". The package-selected fit is retained and flagged; no internal ",
      "replacement or hand-edited latent embedding is used.",
      call. = FALSE
    )
  }
  pred <- LVGP::LVGP_predict(test_x, fit, MSE_on = 1)
  latent_var <- predictive_variance_diagonal(
    pred$MSE, nrow(test_x), context = "LVGP"
  )

  ## LVGP represents the fitted nugget on the internally standardized response
  ## correlation scale. Convert it back to the original response scale before
  ## adding it to the latent-mean MSE returned by LVGP_predict().
  response_range <- as.numeric(fit$data$Y_max - fit$data$Y_min)
  noise_var <- as.numeric(
    fit$fit_detail$sigma2 * fit$fit_detail$nug_opt * response_range^2
  )
  if (length(noise_var) != 1L || !is.finite(noise_var) || noise_var < 0) {
    stop("LVGP returned an invalid fitted noise variance.")
  }
  pred_var_y <- predictive_variance_diagonal(
    latent_var + noise_var,
    nrow(test_x),
    context = "LVGP noisy prediction"
  )
  latent_mean <- as.numeric(pred$Y_hat)

  predictive_draws <- sample_independent_predictive_marginals(
    latent_mean, pred_var_y, n_draw, seed + 10000L
  )
  predictive_draws <- attach_predictive_normal_components(
    predictive_draws, latent_mean, pred_var_y
  )

  list(
    draws = predictive_draws,
    latent_mean = latent_mean,
    predictive_mean = latent_mean,
    predictive_variance = pred_var_y,
    fit = fit,
    optimizer_attempts = attempt,
    optimization_status = if (accepted_package_nonconvergence) {
      if (max_retries > 1L) {
        "package_selected_nonconverged_after_rescue"
      } else {
        "package_selected_nonconverged"
      }
    } else if (best_retry == 1L) {
      "default_converged"
    } else {
      "rescued_converged"
    },
    fitted_noise_variance = noise_var,
    implementation = paste0(
      "LVGP ", utils::packageVersion("LVGP"),
      "::LVGP_fit(noise=TRUE, dim_z=", dim_z,
      ", n_opt=", n_starts, "); author-package-selected finite fit after ",
      "a prespecified sequence of default and extended-limit attempts; ",
      "nonzero optim codes are retained only after rescue exhaustion and ",
      "are explicitly flagged"
    )
  )
}

fit_ezgp_at_tau <- function(train_x,
                            y_train,
                            p,
                            q,
                            m_vec,
                            tau,
                            maxeval) {
  EzGP::EzGP_fit(
    X = train_x,
    Y = as.numeric(y_train),
    p = p,
    q = q,
    m = m_vec,
    tau = tau,
    maxeval = as.integer(maxeval)
  )
}

select_ezgp_tau_cv <- function(train_x,
                              y_train,
                              p,
                              q,
                              m_vec,
                              seed,
                              tau_fractions,
                              cv_folds,
                              maxeval,
                              cv_score = c("nlpd", "mse")) {
  cv_score <- match.arg(cv_score)
  n <- nrow(train_x)
  cv_folds <- max(2L, min(as.integer(cv_folds), n))
  set.seed(seed)
  fold_id <- sample(rep(seq_len(cv_folds), length.out = n))
  response_var <- max(stats::var(y_train), .Machine$double.eps)
  tau_grid <- unique(as.numeric(tau_fractions) * response_var)
  cv_loss <- rep(Inf, length(tau_grid))

  for (tt in seq_along(tau_grid)) {
    fold_loss <- rep(NA_real_, cv_folds)
    for (fold in seq_len(cv_folds)) {
      train_id <- fold_id != fold
      validation_id <- !train_id
      ans <- tryCatch({
        fit <- fit_ezgp_at_tau(
          train_x[train_id, , drop = FALSE],
          y_train[train_id], p, q, m_vec, tau_grid[tt], maxeval
        )
        pred <- EzGP::EzGP_predict(
          train_x[validation_id, , drop = FALSE], fit,
          MSE_on = as.integer(cv_score == "nlpd")
        )
        error <- as.numeric(pred$Y_hat) - y_train[validation_id]
        if (cv_score == "mse") {
          mean(error^2)
        } else {
          latent_var <- predictive_variance_diagonal(
            pred$MSE, sum(validation_id), context = "EzGP cross-validation"
          )
          predictive_var <- predictive_variance_diagonal(
            latent_var + tau_grid[tt],
            sum(validation_id),
            context = "EzGP cross-validation noisy prediction"
          )
          if (any(predictive_var <= 0)) {
            stop("EzGP cross-validation returned zero predictive variance.")
          }
          mean(
            0.5 * log(2 * base::pi * predictive_var) +
              0.5 * error^2 / predictive_var
          )
        }
      }, error = function(e) NA_real_)
      fold_loss[fold] <- ans
    }
    if (all(is.finite(fold_loss))) cv_loss[tt] <- mean(fold_loss)
  }
  if (!any(is.finite(cv_loss))) {
    stop("EzGP training-only cross-validation failed for every nugget value.")
  }
  best <- which.min(cv_loss)
  list(
    tau = tau_grid[best], tau_grid = tau_grid, cv_loss = cv_loss,
    cv_score = cv_score
  )
}

mixedgp_adapter_ezgp <- function(X_train,
                                 y_train,
                                 C_train,
                                 X_test,
                                 C_test,
                                 m_vec,
                                 n_draw,
                                 seed,
                                 tau_fractions = c(1e-6, 0.0025, 0.01, 0.04, 0.16),
                                 cv_folds = 3L,
                                 maxeval = 100L,
                                 cv_score = c("nlpd", "mse")) {
  cv_score <- match.arg(cv_score)
  if (!requireNamespace("EzGP", quietly = TRUE)) {
    stop("Package 'EzGP' is not installed.")
  }
  checked <- validate_mixedgp_inputs(X_train, C_train, m_vec)
  p <- ncol(checked$X)
  q <- ncol(checked$C)
  train_x <- make_mixedgp_numeric_matrix(
    checked$X, checked$C, checked$m_vec
  )
  test_x <- make_mixedgp_numeric_matrix(X_test, C_test, checked$m_vec)

  tau_selection <- select_ezgp_tau_cv(
    train_x = train_x,
    y_train = as.numeric(y_train),
    p = p,
    q = q,
    m_vec = checked$m_vec,
    seed = seed,
    tau_fractions = tau_fractions,
    cv_folds = cv_folds,
    maxeval = maxeval,
    cv_score = cv_score
  )
  fit <- fit_ezgp_at_tau(
    train_x, y_train, p, q, checked$m_vec,
    tau_selection$tau, maxeval
  )
  if (length(fit$param) == 0L || any(!is.finite(fit$param))) {
    stop("EzGP returned non-finite fitted parameters.")
  }
  pred <- EzGP::EzGP_predict(test_x, fit, MSE_on = 1)
  latent_var <- predictive_variance_diagonal(
    pred$MSE, nrow(test_x), context = "EzGP"
  )
  pred_var_y <- predictive_variance_diagonal(
    latent_var + tau_selection$tau,
    nrow(test_x),
    context = "EzGP noisy prediction"
  )
  latent_mean <- as.numeric(pred$Y_hat)

  predictive_draws <- sample_independent_predictive_marginals(
    latent_mean, pred_var_y, n_draw, seed + 10000L
  )
  predictive_draws <- attach_predictive_normal_components(
    predictive_draws, latent_mean, pred_var_y
  )

  list(
    draws = predictive_draws,
    latent_mean = latent_mean,
    predictive_mean = latent_mean,
    predictive_variance = pred_var_y,
    fit = fit,
    fitted_noise_variance = tau_selection$tau,
    tuning = tau_selection,
    implementation = paste0(
      "EzGP ", utils::packageVersion("EzGP"),
      "::EzGP_fit; nugget selected by ", cv_folds,
      "-fold training-only ", toupper(cv_score),
      " CV over prespecified variance fractions"
    )
  )
}

validate_mixedgp_competitor_result <- function(result, n_draw, n_test) {
  if (!is.list(result) || is.null(result$draws)) {
    stop("A competitor adapter must return a list containing 'draws'.")
  }
  component_means <- attr(result$draws, "conditional_means", exact = TRUE)
  component_vars <- attr(result$draws, "conditional_vars", exact = TRUE)
  result$draws <- as.matrix(result$draws)
  if (!is.null(component_means)) {
    attr(result$draws, "conditional_means") <- component_means
  }
  if (!is.null(component_vars)) {
    attr(result$draws, "conditional_vars") <- component_vars
  }
  expected <- c(as.integer(n_draw), as.integer(n_test))
  if (!identical(dim(result$draws), expected)) {
    stop(
      "Competitor predictive draws must have dimension n_draw by n_test; got ",
      paste(dim(result$draws), collapse = " by "), "."
    )
  }
  if (any(!is.finite(result$draws))) {
    stop("Competitor adapter returned non-finite predictive draws.")
  }

  ## Adapters in this file return exact package point predictors.  The
  ## fallbacks preserve compatibility with older/custom adapter results that
  ## contained only noisy-response draws; they are not used by the audited
  ## publication adapters above.
  if (is.null(result$latent_mean)) {
    result$latent_mean <- result$predictive_mean
  }
  if (is.null(result$latent_mean)) {
    result$latent_mean <- colMeans(result$draws)
    result$mean_source <- "monte_carlo_draw_average_compatibility_fallback"
  }
  result$latent_mean <- as.numeric(result$latent_mean)
  if (length(result$latent_mean) != n_test ||
      any(!is.finite(result$latent_mean))) {
    stop("Competitor latent mean must be finite with length n_test.")
  }
  if (is.null(result$predictive_mean)) {
    result$predictive_mean <- result$latent_mean
  }
  result$predictive_mean <- as.numeric(result$predictive_mean)
  if (length(result$predictive_mean) != n_test ||
      any(!is.finite(result$predictive_mean))) {
    stop("Competitor predictive mean must be finite with length n_test.")
  }
  if (is.null(result$predictive_variance) && !is.null(component_vars) &&
      nrow(as.matrix(component_vars)) == 1L) {
    result$predictive_variance <- as.numeric(component_vars)
  }
  if (!is.null(result$predictive_variance)) {
    result$predictive_variance <- as.numeric(result$predictive_variance)
    if (length(result$predictive_variance) != n_test ||
        any(!is.finite(result$predictive_variance)) ||
        any(result$predictive_variance < 0)) {
      stop("Competitor predictive variance must be nonnegative with length n_test.")
    }
    if (is.null(component_means) || is.null(component_vars)) {
      result$draws <- attach_predictive_normal_components(
        result$draws,
        result$predictive_mean,
        result$predictive_variance
      )
    }
  }
  result$mean_source <- result$mean_source %||% "public_package_point_predictor"
  result
}

validate_study2_competitor_result <- validate_mixedgp_competitor_result

run_published_mixedgp_competitors <- function(
    X_train,
    y_train,
    C_train,
    X_test,
    C_test,
    n_draw,
    seed,
    m_vec = NULL,
    methods = mixedgp_competitor_registry()$method,
    strict = FALSE,
    controls = list()) {
  adapters <- list(
    `UC-GP` = mixedgp_adapter_ucgp,
    LVGP = mixedgp_adapter_lvgp,
    EzGP = mixedgp_adapter_ezgp
  )
  unknown <- setdiff(methods, names(adapters))
  if (length(unknown) > 0L) {
    stop("Unknown published competitor(s): ", paste(unknown, collapse = ", "))
  }
  checked <- validate_mixedgp_inputs(X_train, C_train, m_vec)
  validate_mixedgp_inputs(X_test, C_test, checked$m_vec)

  if (length(methods) == 0L) {
    return(list(
      draws = list(), means = list(), latent_means = list(),
      predictive_means = list(), predictive_variances = list(), fits = list(),
      status = mixedgp_competitor_status(
        method = character(0), status = character(0),
        optimization_status = character(0),
        optimization_attempts = character(0), message = character(0),
        warnings = character(0), elapsed = numeric(0),
        implementation = character(0), fitted_noise_variance = numeric(0)
      )
    ))
  }

  draws <- list()
  latent_means <- list()
  predictive_means <- list()
  predictive_variances <- list()
  fits <- list()
  status <- list()
  for (j in seq_along(methods)) {
    method <- methods[j]
    start <- proc.time()[3]
    common_args <- list(
      X_train = checked$X,
      y_train = as.numeric(y_train),
      C_train = checked$C,
      X_test = as.matrix(X_test),
      C_test = as.matrix(C_test),
      m_vec = checked$m_vec,
      n_draw = as.integer(n_draw),
      seed = as.integer(seed + j)
    )
    method_args <- controls[[method]] %||% list()
    duplicated_args <- intersect(names(common_args), names(method_args))
    if (length(duplicated_args) > 0L) {
      stop(
        "Do not override common competitor arguments in controls: ",
        paste(duplicated_args, collapse = ", ")
      )
    }
    captured_warnings <- character(0)
    ans <- tryCatch(
      withCallingHandlers(
        do.call(adapters[[method]], c(common_args, method_args)),
        warning = function(w) {
          captured_warnings <<- c(captured_warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )
    elapsed <- proc.time()[3] - start
    warning_text <- paste(
      unique(trimws(gsub("[[:space:]]+", " ", captured_warnings))),
      collapse = " | "
    )

    if (inherits(ans, "error")) {
      failure_attempts <- ans$optimizer_attempts %||% NULL
      failure_status <- if (!is.null(failure_attempts) &&
          any(failure_attempts$status == "timeout")) {
        "timeout"
      } else if (!is.null(failure_attempts) &&
                 any(failure_attempts$status == "nonconverged")) {
        "nonconverged"
      } else {
        "failed"
      }
      status[[method]] <- mixedgp_competitor_status(
        method = method,
        status = "unavailable_or_failed",
        optimization_status = failure_status,
        optimization_attempts = format_mixedgp_optimizer_attempts(
          failure_attempts
        ),
        message = conditionMessage(ans),
        warnings = warning_text,
        elapsed = elapsed
      )
      if (isTRUE(strict)) stop(ans)
      next
    }

    ans <- validate_mixedgp_competitor_result(
      ans, n_draw = n_draw, n_test = nrow(X_test)
    )
    draws[[method]] <- ans$draws
    latent_means[[method]] <- ans$latent_mean
    predictive_means[[method]] <- ans$predictive_mean
    predictive_variances[[method]] <- ans$predictive_variance %||% NA_real_
    fits[[method]] <- ans$fit %||% NULL
    status[[method]] <- mixedgp_competitor_status(
      method = method,
      status = "success",
      optimization_status = ans$optimization_status %||% "converged",
      optimization_attempts = format_mixedgp_optimizer_attempts(
        ans$optimizer_attempts %||% NULL
      ),
      warnings = warning_text,
      elapsed = elapsed,
      implementation = ans$implementation %||% "public package",
      fitted_noise_variance = ans$fitted_noise_variance %||% NA_real_
    )
  }

  list(
    draws = draws,
    ## `means` is a concise alias for callers interested in m(x,c).  The two
    ## explicit names make the estimand transparent and remain identical for
    ## the zero-mean response-noise models fitted by these packages.
    means = latent_means,
    latent_means = latent_means,
    predictive_means = predictive_means,
    predictive_variances = predictive_variances,
    fits = fits,
    status = do.call(rbind, status)
  )
}

run_study2_published_competitors <- run_published_mixedgp_competitors
run_study1_published_competitors <- run_published_mixedgp_competitors

## Backward-compatible adapter name used in older validation code.
study2_adapter_ucgp_qian <- mixedgp_adapter_ucgp
