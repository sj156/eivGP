############################################################
## 08_published_competitor_validation.R
##
## Exact small-data interface validation for the public-package competitor
## wrappers used in Studies I and II. The validation is defined as a function
## so publication runners can execute it in an isolated source environment.
############################################################

run_published_competitor_validation <- function(require_all = FALSE,
                                                output_file = NULL) {
  if (!is.logical(require_all) || length(require_all) != 1L ||
      is.na(require_all)) {
    stop("require_all must be TRUE or FALSE.")
  }
  if (!is.null(output_file) &&
      (!is.character(output_file) || length(output_file) != 1L ||
       is.na(output_file) || !nzchar(output_file))) {
    stop("output_file must be NULL or one nonempty path.")
  }

  caller_rng <- mixedgp_rng_state()
  on.exit(mixedgp_restore_rng_state(caller_rng), add = TRUE)

  set.seed(20260830L)
  n_train <- 32L
  n_test <- 7L
  m_vec <- rep(4L, 4L)

  X_train <- cbind(
    x1 = stats::rnorm(n_train),
    x2 = stats::runif(n_train, -1, 1)
  )
  C_train <- cbind(
    c1 = rep(1:4, length.out = n_train),
    c2 = rep(1:4, each = 2L, length.out = n_train),
    c3 = rep(c(1L, 3L, 2L, 4L), each = 4L, length.out = n_train),
    c4 = rep(4:1, each = 8L, length.out = n_train)
  )
  y_factor_effect <- matrix(
    c(
      -0.50, 0.10, 0.45, 0.80,
      -0.20, 0.25, 0.40, -0.10,
      -0.15, 0.05, 0.20, 0.35,
      0.30, 0.10, -0.10, -0.25
    ),
    nrow = 4L,
    byrow = TRUE
  )
  y_train <-
    sin(X_train[, 1L]) +
    0.4 * X_train[, 2L] +
    rowSums(vapply(seq_len(ncol(C_train)), function(j) {
      y_factor_effect[j, C_train[, j]]
    }, numeric(n_train))) +
    stats::rnorm(n_train, sd = 0.08)
  X_test <- cbind(
    x1 = seq(-0.8, 0.8, length.out = n_test),
    x2 = seq(0.7, -0.7, length.out = n_test)
  )
  C_test <- cbind(
    c1 = rep(1:4, length.out = n_test),
    c2 = rep(4:1, length.out = n_test),
    c3 = rep(c(1L, 3L, 2L, 4L), length.out = n_test),
    c4 = rep(c(2L, 4L, 1L, 3L), length.out = n_test)
  )

  checked <- validate_mixedgp_inputs(X_train, C_train, m_vec)
  stopifnot(
    identical(dim(checked$X), c(n_train, 2L)),
    identical(dim(checked$C), c(n_train, 4L)),
    identical(checked$m_vec, m_vec)
  )
  factor_data <- make_mixedgp_competitor_data(X_test, C_test, m_vec)
  stopifnot(
    identical(names(factor_data), c("x1", "x2", "c1", "c2", "c3", "c4")),
    identical(levels(factor_data$c1), as.character(seq_len(m_vec[1L]))),
    identical(levels(factor_data$c2), as.character(seq_len(m_vec[2L])))
  )
  numeric_data <- make_mixedgp_numeric_matrix(X_test, C_test, m_vec)
  stopifnot(
    identical(dim(numeric_data), c(n_test, 6L)),
    all(is.finite(numeric_data))
  )

  wrapper_formals <- names(formals(run_published_mixedgp_competitors))
  stopifnot(!any(grepl("(^|_)U($|_)", wrapper_formals)))
  stopifnot(all(vapply(
    c("fit_ucgp", "fit_lvgp", "fit_ezgp", "fit_mixedgp_competitor"),
    function(name) is.function(get(name, mode = "function")),
    logical(1L)
  )))

  draw_check <- sample_independent_predictive_marginals(
    mean = seq_len(n_test), variance = rep(0.25, n_test),
    n_draw = 11L, seed = 1L
  )
  stopifnot(identical(dim(draw_check), c(11L, n_test)))

  bad_levels <- C_test
  bad_levels[1L, 1L] <- m_vec[1L] + 1L
  bad_level_error <- try(
    validate_mixedgp_inputs(X_test, bad_levels, m_vec),
    silent = TRUE
  )
  stopifnot(inherits(bad_level_error, "try-error"))

  preflight <- mixedgp_competitor_preflight(strict = FALSE)
  fit_check <- run_published_mixedgp_competitors(
    X_train = X_train,
    y_train = y_train,
    C_train = C_train,
    X_test = X_test,
    C_test = C_test,
    n_draw = 11L,
    seed = 20260831L,
    m_vec = m_vec,
    strict = FALSE,
    controls = list(
      `UC-GP` = list(n_starts = 8L),
      LVGP = list(n_starts = 8L, max_iter_ini = 100L, max_iter_lat = 20L),
      EzGP = list(tau_fractions = 0.01, cv_folds = 2L, maxeval = 30L)
    )
  )

  stopifnot(
    identical(sort(fit_check$status$method), sort(preflight$method)),
    all(fit_check$status$status %in% c("success", "unavailable_or_failed"))
  )
  for (method in names(fit_check$draws)) {
    component_means <- attr(
      fit_check$draws[[method]], "conditional_means", exact = TRUE
    )
    component_vars <- attr(
      fit_check$draws[[method]], "conditional_vars", exact = TRUE
    )
    stopifnot(
      identical(dim(fit_check$draws[[method]]), c(11L, n_test)),
      all(is.finite(fit_check$draws[[method]])),
      identical(dim(component_means), c(1L, n_test)),
      identical(dim(component_vars), c(1L, n_test)),
      all(component_vars > 0),
      isTRUE(all.equal(
        as.numeric(component_vars),
        as.numeric(fit_check$predictive_variances[[method]])
      ))
    )
  }

  preflight_for_merge <- preflight
  names(preflight_for_merge)[
    names(preflight_for_merge) == "implementation"
  ] <- "registered_implementation"
  status_for_merge <- fit_check$status
  names(status_for_merge)[
    names(status_for_merge) == "implementation"
  ] <- "executed_implementation"
  validation <- merge(
    preflight_for_merge, status_for_merge,
    by = "method", all.x = TRUE, sort = FALSE
  )
  validation <- validation[match(preflight$method, validation$method), ]

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(validation, output_file, row.names = FALSE, na = "")
  }
  if (isTRUE(require_all) && any(validation$status != "success")) {
    failed <- validation$method[validation$status != "success"]
    stop(
      "Strict competitor validation failed for: ",
      paste(failed, collapse = ", "), "."
    )
  }
  validation
}

if (sys.nframe() == 0L) {
  validator_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  validator_dir <- if (length(validator_arg) > 0L) {
    dirname(normalizePath(sub("^--file=", "", validator_arg[[1L]])))
  } else {
    normalizePath(getwd())
  }
  configured_library <- Sys.getenv("MIXEDGP_R_LIBRARY", unset = "")
  if (nzchar(configured_library)) {
    .libPaths(unique(c(normalizePath(configured_library), .libPaths())))
  }
  for (module in c(
    "00_parallel_utils.R", "00_study1_functions.R", "00_study2_functions.R",
    "03_study2_published_competitors.R", "00_public_api.R"
  )) {
    sys.source(file.path(validator_dir, module), envir = .GlobalEnv)
  }
  strict <- tolower(Sys.getenv(
    "MIXEDGP_REQUIRE_ALL_COMPETITORS", unset = "false"
  )) %in% c("1", "true", "yes", "y")
  print(run_published_competitor_validation(require_all = strict))
}
