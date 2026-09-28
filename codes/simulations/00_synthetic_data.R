############################################################
## Frozen synthetic datasets for reproducible experiments
############################################################

MIXEDGP_DATA_SCHEMA_VERSION <- "1.1.0"
MIXEDGP_STUDY1_GENERATOR_TAG <- "study1-synthetic-heterogeneity-v2"
MIXEDGP_STUDY2_GENERATOR_TAG <- "study2-synthetic-shared-state-v2"

mixedgp_object_fingerprint <- function(object) {
  path <- tempfile("mixedgp-object-", fileext = ".rds")
  on.exit(if (file.exists(path)) unlink(path), add = TRUE)
  saveRDS(object, path, version = 3L, compress = FALSE)
  unname(tools::md5sum(path))
}

mixedgp_expected_generator_tag <- function(study) {
  study <- match.arg(study, c("study1", "study2"))
  if (study == "study1") return(MIXEDGP_STUDY1_GENERATOR_TAG)
  MIXEDGP_STUDY2_GENERATOR_TAG
}

mixedgp_synthetic_design_fingerprint <- function(study,
                                                  scenario,
                                                  design,
                                                  generator_tag) {
  payload <- list(
    schema_version = MIXEDGP_DATA_SCHEMA_VERSION,
    study = as.character(study),
    scenario = as.character(scenario),
    generator_tag = as.character(generator_tag),
    design = design
  )
  path <- tempfile("mixedgp-design-", fileext = ".rds")
  on.exit(if (file.exists(path)) unlink(path), add = TRUE)
  saveRDS(payload, path, version = 3L, compress = FALSE)
  unname(tools::md5sum(path))
}

mixedgp_dataset_filename <- function(study,
                                     rep_id,
                                     scenario,
                                     n,
                                     n_test,
                                     m,
                                     q = NULL) {
  study <- match.arg(study, c("study1", "study2"))
  scenario <- gsub("[^[:alnum:]_]+", "-", as.character(scenario))
  q_tag <- if (is.null(q)) "" else paste0("_q", as.integer(q))
  sprintf(
    "%s_%s%s_n%d_ntest%d_m%d_rep%03d.rds",
    study, scenario, q_tag, as.integer(n), as.integer(n_test),
    as.integer(m), as.integer(rep_id)
  )
}

make_study1_synthetic_dataset <- function(
    rep_id,
    n = 100L,
    n_test = 500L,
    m = 6L,
    scenario = "active",
    sigma_eps = 0.1,
    threshold_design = "imbalanced",
    min_class_count = NULL,
    heterogeneity_eta = 1,
    calib_grid = c(0L, 5L, 10L, 20L, 50L),
    data_seed_base = 100000L,
    calibration_seed_base = 200000L) {
  rep_id <- mixedgp_as_integer_strict(
    rep_id, "rep_id", min_value = 1L, length_expected = 1L
  )
  n <- mixedgp_as_integer_strict(n, "n", min_value = 1L, length_expected = 1L)
  n_test <- mixedgp_as_integer_strict(
    n_test, "n_test", min_value = 1L, length_expected = 1L
  )
  m <- mixedgp_as_integer_strict(m, "m", min_value = 2L, length_expected = 1L)
  calib_grid <- mixedgp_as_integer_strict(
    calib_grid, "calib_grid", min_value = 0L
  )
  if (any(calib_grid > n)) stop("calib_grid cannot exceed n.")
  data_seed_base <- mixedgp_as_integer_strict(
    data_seed_base, "data_seed_base", min_value = 0L, length_expected = 1L
  )
  calibration_seed_base <- mixedgp_as_integer_strict(
    calibration_seed_base, "calibration_seed_base",
    min_value = 0L, length_expected = 1L
  )
  data_seed <- mixedgp_as_integer_strict(
    as.numeric(data_seed_base) + rep_id, "data seed", min_value = 0L,
    length_expected = 1L
  )
  calibration_seed <- mixedgp_as_integer_strict(
    as.numeric(calibration_seed_base) + rep_id, "calibration seed",
    min_value = 0L, length_expected = 1L
  )
  dat <- simulate_1d_data(
    n = n,
    n_test = n_test,
    m = m,
    scenario = scenario,
    sigma_eps = sigma_eps,
    seed = data_seed,
    threshold_design = threshold_design,
    min_class_count = min_class_count,
    heterogeneity_eta = heterogeneity_eta
  )
  calibration_sets <- make_nested_calibration_sets(
    n = n,
    calib_grid = as.integer(calib_grid),
    seed = calibration_seed
  )
  design <- list(
    n = as.integer(n), n_test = as.integer(n_test), m = as.integer(m),
    sigma_eps = sigma_eps, threshold_design = threshold_design,
    min_class_count = dat$min_class_count,
    heterogeneity_eta = heterogeneity_eta,
    calib_grid = as.integer(calib_grid),
    generator_tag = mixedgp_expected_generator_tag("study1")
  )
  design$design_fingerprint <- mixedgp_synthetic_design_fingerprint(
    "study1", scenario, design, design$generator_tag
  )
  structure(
    list(
      schema_version = MIXEDGP_DATA_SCHEMA_VERSION,
      study = "study1",
      rep_id = rep_id,
      scenario = scenario,
      seeds = list(data = data_seed, calibration = calibration_seed),
      design = design,
      data = dat,
      calibration_sets = calibration_sets
    ),
    class = c("mixedgp_synthetic_dataset", "list")
  )
}

make_study2_synthetic_dataset <- function(
    rep_id,
    scenario = "primary",
    n = 120L,
    n_test = 500L,
    q = 4L,
    m = 4L,
    sigma_eps = NULL,
    calib_grid = c(0L, 10L, 25L, 50L, 80L),
    data_seed_base = 1000000L,
    calibration_seed_base = NULL) {
  rep_id <- mixedgp_as_integer_strict(
    rep_id, "rep_id", min_value = 1L, length_expected = 1L
  )
  n <- mixedgp_as_integer_strict(n, "n", min_value = 1L, length_expected = 1L)
  n_test <- mixedgp_as_integer_strict(
    n_test, "n_test", min_value = 1L, length_expected = 1L
  )
  q <- mixedgp_as_integer_strict(q, "q", min_value = 2L, length_expected = 1L)
  m <- mixedgp_as_integer_strict(m, "m", min_value = 2L, length_expected = 1L)
  calib_grid <- mixedgp_as_integer_strict(
    calib_grid, "calib_grid", min_value = 0L
  )
  if (any(calib_grid > n)) stop("calib_grid cannot exceed n.")
  data_seed_base <- mixedgp_as_integer_strict(
    data_seed_base, "data_seed_base", min_value = 0L, length_expected = 1L
  )
  data_seed_base_rep <- mixedgp_as_integer_strict(
    as.numeric(data_seed_base) + 10000 * rep_id,
    "derived Study II data seed base", min_value = 0L, length_expected = 1L
  )
  if (data_seed_base_rep > .Machine$integer.max - 8L) {
    stop("Derived Study II seed base leaves no room for component seeds.")
  }
  data_seed <- data_seed_base_rep + 1L
  if (is.null(calibration_seed_base)) {
    calibration_seed <- data_seed_base_rep + 2L
  } else {
    calibration_seed_base <- mixedgp_as_integer_strict(
      calibration_seed_base, "calibration_seed_base",
      min_value = 0L, length_expected = 1L
    )
    calibration_seed <- mixedgp_as_integer_strict(
      as.numeric(calibration_seed_base) + rep_id, "calibration seed",
      min_value = 0L, length_expected = 1L
    )
  }
  dat <- simulate_study2_data(
    n = n,
    n_test = n_test,
    seed = data_seed,
    scenario = scenario,
    sigma_eps = sigma_eps,
    q = q,
    m = m
  )
  calibration_sets <- make_stratified_calibration_sets_2d(
    C = dat$train$C,
    calib_grid = as.integer(calib_grid),
    seed = calibration_seed
  )
  generator_tag <- mixedgp_expected_generator_tag("study2")
  design <- list(
    n = as.integer(n), n_test = as.integer(n_test), q = as.integer(q),
    m = as.integer(m), d = as.integer(dat$true_params$d),
    sigma_eps = dat$sigma_eps,
    calib_grid = as.integer(calib_grid), generator_tag = generator_tag
  )
  design$design_fingerprint <- mixedgp_synthetic_design_fingerprint(
    "study2", scenario, design, generator_tag
  )
  structure(
    list(
      schema_version = MIXEDGP_DATA_SCHEMA_VERSION,
      study = "study2",
      rep_id = rep_id,
      scenario = scenario,
      seeds = list(data = data_seed, calibration = calibration_seed),
      design = design,
      data = dat,
      calibration_sets = calibration_sets
    ),
    class = c("mixedgp_synthetic_dataset", "list")
  )
}

validate_mixedgp_synthetic_dataset <- function(x) {
  if (!inherits(x, "mixedgp_synthetic_dataset")) {
    stop("Object is not a mixedgp_synthetic_dataset.")
  }
  required <- c(
    "schema_version", "study", "rep_id", "scenario", "seeds", "design",
    "data", "calibration_sets"
  )
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop("Synthetic dataset is missing: ", paste(missing, collapse = ", "))
  }
  if (!identical(x$schema_version, MIXEDGP_DATA_SCHEMA_VERSION)) {
    stop(
      "Unsupported synthetic-data schema ", x$schema_version,
      "; expected ", MIXEDGP_DATA_SCHEMA_VERSION, "."
    )
  }
  if (length(x$study) != 1L || !x$study %in% c("study1", "study2") ||
      length(x$scenario) != 1L || is.na(x$scenario) || !nzchar(x$scenario)) {
    stop("Synthetic dataset has invalid study, scenario, or replication metadata.")
  }
  rep_id <- mixedgp_as_integer_strict(
    x$rep_id, "rep_id", min_value = 1L, length_expected = 1L
  )
  if (!is.list(x$seeds) ||
      !identical(names(x$seeds), c("data", "calibration"))) {
    stop("Synthetic dataset seeds must contain data and calibration entries.")
  }
  invisible(mixedgp_as_integer_strict(
    unlist(x$seeds, use.names = FALSE), "synthetic seeds",
    min_value = 0L, length_expected = 2L
  ))
  design_required <- c(
    "n", "n_test", "m", "calib_grid", "generator_tag",
    "design_fingerprint"
  )
  if (x$study == "study1") {
    design_required <- c(
      design_required, "sigma_eps", "threshold_design", "min_class_count",
      "heterogeneity_eta"
    )
  } else {
    design_required <- c(design_required, "q", "d", "sigma_eps")
  }
  missing_design <- setdiff(design_required, names(x$design))
  if (length(missing_design) > 0L) {
    stop("Synthetic design is missing: ", paste(missing_design, collapse = ", "))
  }
  generator_tag <- x$design$generator_tag
  design_fingerprint <- x$design$design_fingerprint
  design_without_fingerprint <- x$design[
    setdiff(names(x$design), "design_fingerprint")
  ]
  expected_fingerprint <- mixedgp_synthetic_design_fingerprint(
    x$study, x$scenario, design_without_fingerprint, generator_tag
  )
  if (length(generator_tag) != 1L || is.na(generator_tag) ||
      !nzchar(generator_tag) || length(design_fingerprint) != 1L ||
      is.na(design_fingerprint) ||
      !identical(
        as.character(generator_tag), mixedgp_expected_generator_tag(x$study)
      ) ||
      !identical(as.character(design_fingerprint), expected_fingerprint)) {
    stop("Synthetic generator tag or design fingerprint is invalid.")
  }
  n <- mixedgp_as_integer_strict(
    x$design$n, "design$n", min_value = 1L, length_expected = 1L
  )
  n_test <- mixedgp_as_integer_strict(
    x$design$n_test, "design$n_test", min_value = 1L, length_expected = 1L
  )
  m <- mixedgp_as_integer_strict(
    x$design$m, "design$m", min_value = 2L, length_expected = 1L
  )
  calib_grid <- mixedgp_as_integer_strict(
    x$design$calib_grid, "design$calib_grid", min_value = 0L
  )
  if (any(calib_grid > n)) {
    stop("Synthetic dataset calibration sizes cannot exceed n.")
  }
  if (!is.numeric(x$design$sigma_eps) || length(x$design$sigma_eps) != 1L ||
      !is.finite(x$design$sigma_eps) || x$design$sigma_eps < 0) {
    stop("design$sigma_eps must be one nonnegative finite number.")
  }
  if (!is.list(x$calibration_sets) ||
      !identical(names(x$calibration_sets), as.character(calib_grid)) ||
      length(x$calibration_sets) != length(calib_grid)) {
    stop("Calibration sets do not match the declared calibration grid.")
  }
  valid_calibration <- vapply(seq_along(calib_grid), function(ii) {
    idx <- x$calibration_sets[[ii]]
    if (calib_grid[ii] == 0L) return(length(idx) == 0L)
    valid_idx <- try(
      mixedgp_as_integer_strict(idx, "calibration row indices", min_value = 1L),
      silent = TRUE
    )
    !inherits(valid_idx, "try-error") &&
      identical(length(valid_idx), calib_grid[ii]) &&
      length(unique(valid_idx)) == length(valid_idx) && all(valid_idx <= n)
  }, logical(1L))
  if (!all(valid_calibration)) {
    stop("A calibration set has invalid size, duplicates, or row indices.")
  }
  grid_order <- order(calib_grid)
  if (length(grid_order) > 1L) {
    nested <- vapply(seq_len(length(grid_order) - 1L), function(ii) {
      all(x$calibration_sets[[grid_order[ii]]] %in%
            x$calibration_sets[[grid_order[ii + 1L]]])
    }, logical(1L))
    if (!all(nested)) stop("Calibration sets must be nested by size.")
  }

  if (x$study == "study1") {
    required_values <- c("x", "u", "c", "y", "f")
    numeric_train <- is.data.frame(x$data$train) &&
      all(required_values %in% names(x$data$train)) &&
      all(vapply(x$data$train[required_values], is.numeric, logical(1L)))
    numeric_test <- is.data.frame(x$data$test) &&
      all(required_values %in% names(x$data$test)) &&
      all(vapply(x$data$test[required_values], is.numeric, logical(1L)))
    if (!numeric_train || !numeric_test) {
      stop("Study I train/test truth columns must be numeric.")
    }
    if (!is.numeric(x$data$tau_true) || length(x$data$tau_true) != m - 1L ||
        any(!is.finite(x$data$tau_true)) ||
        is.unsorted(x$data$tau_true, strictly = TRUE) ||
        !identical(as.character(x$data$scenario), as.character(x$scenario)) ||
        !identical(as.integer(x$data$m), m) ||
        !isTRUE(all.equal(
          as.numeric(x$data$sigma_eps), as.numeric(x$design$sigma_eps),
          tolerance = 0
        ))) {
      stop("Study I DGP metadata or cutpoints are inconsistent with the design.")
    }
    if (!is.data.frame(x$data$train) || !is.data.frame(x$data$test) ||
        nrow(x$data$train) != n || nrow(x$data$test) != n_test ||
        anyNA(x$data$train[, required_values]) ||
        anyNA(x$data$test[, required_values]) ||
        any(!is.finite(as.matrix(x$data$train[, required_values]))) ||
        any(!is.finite(as.matrix(x$data$test[, required_values]))) ||
        any(x$data$train$c != floor(x$data$train$c)) ||
        any(x$data$test$c != floor(x$data$test$c)) ||
        any(x$data$train$c < 1L | x$data$train$c > m) ||
        any(x$data$test$c < 1L | x$data$test$c > m)) {
      stop("Study I data dimensions or values do not match the design.")
    }
    expected_train_f <- f0_1d(
      x$data$train$x, x$data$train$u, scenario = x$scenario,
      c_ord = x$data$train$c, tau = x$data$tau_true,
      heterogeneity_eta = x$design$heterogeneity_eta
    )
    expected_test_f <- f0_1d(
      x$data$test$x, x$data$test$u, scenario = x$scenario,
      c_ord = x$data$test$c, tau = x$data$tau_true,
      heterogeneity_eta = x$design$heterogeneity_eta
    )
    if (!isTRUE(all.equal(x$data$train$f, expected_train_f, tolerance = 0)) ||
        !isTRUE(all.equal(x$data$test$f, expected_test_f, tolerance = 0))) {
      stop("Study I stored response-surface truth is inconsistent with the DGP.")
    }
  } else {
    q <- mixedgp_as_integer_strict(
      x$design$q, "design$q", min_value = 2L, length_expected = 1L
    )
    d <- mixedgp_as_integer_strict(
      x$design$d, "design$d", min_value = 1L, length_expected = 1L
    )
    train_U <- as.matrix(x$data$train$U)
    test_U <- as.matrix(x$data$test$U)
    train_f <- x$data$train$f
    test_f <- x$data$test$f
    true_params <- x$data$true_params
    numeric_objects <- list(
      x$data$train$X, x$data$test$X, train_U, test_U,
      x$data$train$C, x$data$test$C, x$data$train$y,
      x$data$test$y, train_f, test_f
    )
    if (any(!vapply(
      numeric_objects,
      function(value) is.numeric(value) || is.integer(value),
      logical(1L)
    ))) {
      stop("Study II data and truth arrays must be numeric.")
    }
    true_required <- c(
      "scenario", "score_error", "response_interactions", "lambda", "A",
      "Omega", "tau", "q", "m", "d", "sigma_eps"
    )
    if (!is.list(true_params) ||
        length(setdiff(true_required, names(true_params))) > 0L ||
        any(!vapply(
          true_params[c("lambda", "A", "Omega", "tau", "q", "m", "d",
                        "sigma_eps")],
          function(value) is.numeric(value) || is.integer(value),
          logical(1L)
        ))) {
      stop("Study II true_params is incomplete or nonnumeric.")
    }
    if (length(q) != 1L || is.na(q) || q < 2L ||
        length(d) != 1L || is.na(d) || d < 1L ||
        !is.list(x$data$train) || !is.list(x$data$test) ||
        nrow(as.matrix(x$data$train$X)) != n ||
        nrow(as.matrix(x$data$test$X)) != n_test ||
        !identical(dim(as.matrix(x$data$train$C)), c(n, q)) ||
        !identical(dim(as.matrix(x$data$test$C)), c(n_test, q)) ||
        !identical(dim(train_U), c(n, d)) ||
        !identical(dim(test_U), c(n_test, d)) ||
        length(x$data$train$y) != n || length(x$data$test$y) != n_test ||
        length(train_f) != n || length(test_f) != n_test ||
        anyNA(x$data$train$X) || anyNA(x$data$test$X) ||
        anyNA(train_U) || anyNA(test_U) ||
        anyNA(x$data$train$C) || anyNA(x$data$test$C) ||
        anyNA(x$data$train$y) || anyNA(x$data$test$y) ||
        anyNA(train_f) || anyNA(test_f) ||
        any(x$data$train$C != floor(x$data$train$C)) ||
        any(x$data$test$C != floor(x$data$test$C)) ||
        any(!is.finite(c(
          x$data$train$X, x$data$test$X, train_U, test_U,
          x$data$train$y, x$data$test$y, train_f, test_f
        ))) ||
        any(x$data$train$C < 1L | x$data$train$C > m) ||
        any(x$data$test$C < 1L | x$data$test$C > m) ||
        !is.list(true_params) || !identical(as.integer(true_params$q), q) ||
        !identical(as.integer(true_params$m), m) ||
        !identical(as.integer(true_params$d), d) ||
        !identical(dim(as.matrix(true_params$A)), c(q, d)) ||
        !identical(dim(as.matrix(true_params$tau)), c(q, m - 1L)) ||
        !identical(dim(as.matrix(true_params$Omega)), c(q, q)) ||
        any(!is.finite(c(
          true_params$A, true_params$tau, true_params$Omega,
          true_params$lambda, true_params$sigma_eps
        ))) ||
        !identical(as.character(true_params$scenario), as.character(x$scenario))) {
      stop("Study II data dimensions or values do not match the design.")
    }
    expected_params <- make_study2_true_params(
      scenario = x$scenario, q = q, m = m
    )
    if (!isTRUE(all.equal(true_params, expected_params, tolerance = 0)) ||
        !identical(as.character(x$data$scenario), as.character(x$scenario)) ||
        !isTRUE(all.equal(
          as.numeric(x$data$sigma_eps), as.numeric(x$design$sigma_eps),
          tolerance = 0
        ))) {
      stop("Study II true parameters or DGP metadata do not match the design.")
    }
    expected_train_f <- f0_2d(
      x$data$train$X, train_U, scenario = x$scenario
    )
    expected_test_f <- f0_2d(x$data$test$X, test_U, scenario = x$scenario)
    if (!isTRUE(all.equal(train_f, expected_train_f, tolerance = 0)) ||
        !isTRUE(all.equal(test_f, expected_test_f, tolerance = 0))) {
      stop("Study II stored response-surface truth is inconsistent with the DGP.")
    }
  }
  invisible(TRUE)
}

store_mixedgp_synthetic_dataset <- function(x,
                                            directory,
                                            overwrite = FALSE,
                                            compress = "xz") {
  overwrite <- mixedgp_validate_flag(overwrite, "overwrite")
  validate_mixedgp_synthetic_dataset(x)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  design <- x$design
  path <- file.path(
    directory,
    mixedgp_dataset_filename(
      study = x$study,
      rep_id = x$rep_id,
      scenario = x$scenario,
      n = design$n,
      n_test = design$n_test,
      m = design$m,
      q = design$q
    )
  )
  if (file.exists(path) && !isTRUE(overwrite)) {
    existing <- readRDS(path)
    validate_mixedgp_synthetic_dataset(existing)
    if (!identical(existing, x)) {
      stop(
        "A different synthetic dataset already exists at ", path, ". ",
        "The filename is not a sufficient design fingerprint; choose a new ",
        "directory or set overwrite = TRUE after verifying the requested design."
      )
    }
  } else {
    temporary <- tempfile(
      pattern = paste0(basename(path), "."),
      tmpdir = directory
    )
    on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
    saveRDS(x, temporary, version = 3L, compress = compress)
    if (!file.rename(temporary, path)) {
      stop("Could not atomically move the synthetic dataset to ", path, ".")
    }
  }
  info <- file.info(path)
  data.frame(
    schema_version = x$schema_version,
    study = x$study,
    scenario = x$scenario,
    rep = x$rep_id,
    n = design$n,
    n_test = design$n_test,
    q = if (is.null(design$q)) NA_integer_ else design$q,
    d = if (is.null(design$d)) NA_integer_ else design$d,
    m = design$m,
    data_seed = x$seeds$data,
    calibration_seed = x$seeds$calibration,
    file = basename(path),
    bytes = as.numeric(info$size),
    md5 = unname(tools::md5sum(path)),
    generator_tag = as.character(design$generator_tag),
    design_fingerprint = as.character(design$design_fingerprint),
    stringsAsFactors = FALSE
  )
}

load_mixedgp_synthetic_dataset <- function(path, expected_md5 = NULL) {
  if (!file.exists(path)) stop("Synthetic dataset does not exist: ", path)
  if (!is.null(expected_md5)) {
    observed <- unname(tools::md5sum(path))
    if (!identical(tolower(observed), tolower(as.character(expected_md5)))) {
      stop("Checksum mismatch for synthetic dataset: ", path)
    }
  }
  out <- readRDS(path)
  validate_mixedgp_synthetic_dataset(out)
  out
}

load_mixedgp_synthetic_dataset_strict <- function(path,
                                                   expected,
                                                   manifest_path = file.path(
                                                     dirname(path),
                                                     "manifest.rds"
                                                   )) {
  if (!is.list(expected) || is.null(names(expected))) {
    stop("expected must be a named list of frozen-design fields.")
  }
  if (!file.exists(manifest_path)) {
    stop("Frozen-data manifest does not exist: ", manifest_path)
  }
  manifest_object <- readRDS(manifest_path)
  manifest <- if (is.data.frame(manifest_object)) {
    manifest_object
  } else if (is.list(manifest_object) &&
             is.data.frame(manifest_object$manifest)) {
    manifest_object$manifest
  } else {
    stop("Frozen-data manifest has an unsupported structure: ", manifest_path)
  }
  required_columns <- c(
    "schema_version", "study", "scenario", "rep", "n", "n_test", "m",
    "file", "md5", "generator_tag", "design_fingerprint"
  )
  missing_columns <- setdiff(required_columns, names(manifest))
  if (length(missing_columns) > 0L) {
    stop(
      "Frozen-data manifest is missing: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  row_id <- which(as.character(manifest$file) == basename(path))
  if (length(row_id) != 1L) {
    stop(
      "Frozen-data manifest must contain exactly one row for ", basename(path),
      "; found ", length(row_id), "."
    )
  }
  manifest_row <- manifest[row_id, , drop = FALSE]
  out <- load_mixedgp_synthetic_dataset(
    path,
    expected_md5 = as.character(manifest_row$md5[[1L]])
  )

  observed <- c(
    list(
      study = out$study,
      scenario = out$scenario,
      rep_id = out$rep_id
    ),
    out$design
  )
  for (field in names(expected)) {
    if (!field %in% names(observed)) {
      stop("Unknown expected frozen-design field: ", field)
    }
    expected_value <- expected[[field]]
    observed_value <- observed[[field]]
    equal <- if (is.numeric(expected_value) && is.numeric(observed_value)) {
      isTRUE(all.equal(
        as.numeric(observed_value), as.numeric(expected_value),
        tolerance = 0, check.attributes = FALSE
      ))
    } else {
      identical(as.character(observed_value), as.character(expected_value))
    }
    if (!equal) {
      stop(
        "Frozen dataset ", basename(path), " has ", field, "=",
        paste(observed_value, collapse = ","), "; expected ",
        paste(expected_value, collapse = ","), "."
      )
    }
  }

  manifest_checks <- list(
    schema_version = out$schema_version,
    study = out$study,
    scenario = out$scenario,
    rep = out$rep_id,
    n = out$design$n,
    n_test = out$design$n_test,
    m = out$design$m,
    generator_tag = out$design$generator_tag,
    design_fingerprint = out$design$design_fingerprint
  )
  if ("q" %in% names(manifest) && !is.null(out$design$q)) {
    manifest_checks$q <- out$design$q
  }
  if ("d" %in% names(manifest) && !is.null(out$design$d)) {
    manifest_checks$d <- out$design$d
  }
  for (field in names(manifest_checks)) {
    left <- manifest_row[[field]][[1L]]
    right <- manifest_checks[[field]]
    if (is.numeric(left) && is.numeric(right)) {
      equal <- isTRUE(all.equal(
        as.numeric(left), as.numeric(right), tolerance = 0,
        check.attributes = FALSE
      ))
    } else {
      equal <- identical(as.character(left), as.character(right))
    }
    if (!equal) {
      stop(
        "Manifest/artifact mismatch for ", basename(path), " field ", field,
        "."
      )
    }
  }
  attr(out, "manifest_path") <- normalizePath(manifest_path)
  attr(out, "manifest_md5") <- as.character(manifest_row$md5[[1L]])
  out
}

mixedgp_validate_generation_dots <- function(dots, forbidden = "rep_id") {
  if (length(dots) > 0L &&
      (is.null(names(dots)) || anyNA(names(dots)) ||
       any(!nzchar(names(dots))) || anyDuplicated(names(dots)))) {
    stop("Every synthetic-generator argument passed through ... must be named uniquely.")
  }
  duplicate <- intersect(names(dots), forbidden)
  if (length(duplicate) > 0L) {
    stop(
      "Do not pass generator-controlled arguments through ...: ",
      paste(duplicate, collapse = ", "), "."
    )
  }
  invisible(dots)
}

mixedgp_stop_on_worker_errors <- function(rows, context) {
  failed <- which(vapply(rows, inherits, logical(1L), what = "try-error"))
  if (length(failed) == 0L) return(invisible(rows))
  details <- vapply(failed, function(ii) {
    condition <- attr(rows[[ii]], "condition")
    message <- if (inherits(condition, "condition")) {
      conditionMessage(condition)
    } else {
      trimws(as.character(rows[[ii]])[1L])
    }
    paste0("task ", ii, ": ", message)
  }, character(1L))
  stop(
    context, " failed in ", length(failed), " worker(s); no manifest was written. ",
    paste(details, collapse = " | ")
  )
}

mixedgp_bind_manifest_rows <- function(rows) {
  if (length(rows) == 0L) return(data.frame())
  if (!all(vapply(rows, is.data.frame, logical(1L)))) {
    stop("Every successful synthetic-data worker must return a manifest row.")
  }
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(row) {
    missing <- setdiff(columns, names(row))
    for (name in missing) row[[name]] <- NA
    row[, columns, drop = FALSE]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

mixedgp_merge_manifests <- function(existing, incoming) {
  manifests <- Filter(function(x) !is.null(x) && nrow(x) > 0L,
                      list(existing, incoming))
  if (length(manifests) == 0L) return(data.frame())
  if (!all(vapply(manifests, is.data.frame, logical(1L)))) {
    stop("Synthetic-data manifests must be data frames.")
  }
  combined <- mixedgp_bind_manifest_rows(manifests)
  if (!"file" %in% names(combined) || anyNA(combined$file) ||
      any(!nzchar(as.character(combined$file)))) {
    stop("Synthetic-data manifest rows require nonempty file names.")
  }
  groups <- split(seq_len(nrow(combined)), as.character(combined$file))
  keep <- integer(length(groups))
  group_id <- 0L
  for (indices in groups) {
    group_id <- group_id + 1L
    reference <- combined[indices[1L], , drop = FALSE]
    rownames(reference) <- NULL
    consistent <- vapply(indices, function(ii) {
      candidate <- combined[ii, , drop = FALSE]
      rownames(candidate) <- NULL
      identical(as.list(candidate), as.list(reference))
    }, logical(1L))
    if (!all(consistent)) {
      stop(
        "Manifest metadata disagree for existing file ",
        as.character(reference$file[[1L]]), "."
      )
    }
    keep[group_id] <- indices[1L]
  }
  out <- combined[keep, , drop = FALSE]
  ordering_fields <- intersect(c("study", "scenario", "rep", "file"), names(out))
  if (length(ordering_fields) > 0L) {
    ordering <- do.call(order, unname(out[ordering_fields]))
    out <- out[ordering, , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

mixedgp_atomic_replace <- function(temporary, target) {
  if (file.rename(temporary, target)) return(invisible(target))
  backup <- NULL
  if (file.exists(target)) {
    backup <- tempfile(paste0(basename(target), ".backup-"), dirname(target))
    if (!file.rename(target, backup)) {
      stop("Could not prepare atomic replacement of ", target, ".")
    }
  }
  installed <- file.rename(temporary, target)
  if (!installed) {
    if (!is.null(backup) && file.exists(backup)) file.rename(backup, target)
    stop("Could not atomically replace ", target, ".")
  }
  if (!is.null(backup) && file.exists(backup)) unlink(backup)
  invisible(target)
}

mixedgp_write_manifest_atomic <- function(manifest, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  csv_target <- file.path(directory, "manifest.csv")
  rds_target <- file.path(directory, "manifest.rds")
  csv_temporary <- tempfile("manifest-", tmpdir = directory, fileext = ".csv")
  rds_temporary <- tempfile("manifest-", tmpdir = directory, fileext = ".rds")
  on.exit({
    if (file.exists(csv_temporary)) unlink(csv_temporary)
    if (file.exists(rds_temporary)) unlink(rds_temporary)
  }, add = TRUE)
  utils::write.csv(manifest, csv_temporary, row.names = FALSE)
  saveRDS(manifest, rds_temporary, version = 3L)
  mixedgp_atomic_replace(rds_temporary, rds_target)
  mixedgp_atomic_replace(csv_temporary, csv_target)
  invisible(manifest)
}

mixedgp_existing_manifest <- function(directory) {
  rds_path <- file.path(directory, "manifest.rds")
  csv_path <- file.path(directory, "manifest.csv")
  if (!file.exists(rds_path)) {
    if (file.exists(csv_path)) {
      stop("manifest.csv exists without the canonical manifest.rds in ", directory)
    }
    return(NULL)
  }
  manifest <- readRDS(rds_path)
  if (!is.data.frame(manifest)) {
    stop("Existing per-study manifest.rds must contain a data frame.")
  }
  manifest
}

generate_study1_synthetic_datasets <- function(
    n_rep = 50L,
    directory = file.path("..", "data-synthetic", "study1"),
    n_cores = NULL,
    overwrite = FALSE,
    ...) {
  n_rep <- mixedgp_as_integer_strict(
    n_rep, "n_rep", min_value = 1L, length_expected = 1L
  )
  overwrite <- mixedgp_validate_flag(overwrite, "overwrite")
  dots <- list(...)
  mixedgp_validate_generation_dots(dots, forbidden = "rep_id")
  rep_ids <- seq_len(n_rep)
  rows <- mixedgp_parallel_mapply(
    function(rep_id) {
      artifact <- do.call(
        make_study1_synthetic_dataset, c(list(rep_id = rep_id), dots)
      )
      store_mixedgp_synthetic_dataset(artifact, directory, overwrite)
    },
    rep_id = rep_ids,
    n_cores = n_cores,
    seeds = 710000L + rep_ids,
    SIMPLIFY = FALSE,
    mc.preschedule = FALSE
  )
  mixedgp_stop_on_worker_errors(rows, "Study I synthetic-data generation")
  incoming <- mixedgp_bind_manifest_rows(rows)
  manifest <- mixedgp_merge_manifests(
    mixedgp_existing_manifest(directory), incoming
  )
  mixedgp_write_manifest_atomic(manifest, directory)
  manifest
}

generate_study2_synthetic_datasets <- function(
    n_rep = 50L,
    scenarios = c(
      "primary", "latent_additive_control", "high_uncertainty",
      "logistic_misspec"
    ),
    directory = file.path("..", "data-synthetic", "study2"),
    n_cores = NULL,
    overwrite = FALSE,
    ...) {
  n_rep <- mixedgp_as_integer_strict(
    n_rep, "n_rep", min_value = 1L, length_expected = 1L
  )
  overwrite <- mixedgp_validate_flag(overwrite, "overwrite")
  if (!is.character(scenarios) || length(scenarios) < 1L ||
      anyNA(scenarios) || any(!nzchar(scenarios)) || anyDuplicated(scenarios)) {
    stop("scenarios must contain one or more unique, nonempty names.")
  }
  dots <- list(...)
  mixedgp_validate_generation_dots(dots, forbidden = c("rep_id", "scenario"))
  grid <- expand.grid(
    rep_id = seq_len(n_rep),
    scenario = scenarios,
    stringsAsFactors = FALSE
  )
  rows <- mixedgp_parallel_mapply(
    function(rep_id, scenario) {
      artifact <- do.call(
        make_study2_synthetic_dataset,
        c(list(rep_id = rep_id, scenario = scenario), dots)
      )
      store_mixedgp_synthetic_dataset(artifact, directory, overwrite)
    },
    rep_id = grid$rep_id,
    scenario = grid$scenario,
    n_cores = n_cores,
    seeds = 720000L + seq_len(nrow(grid)),
    SIMPLIFY = FALSE,
    mc.preschedule = FALSE
  )
  mixedgp_stop_on_worker_errors(rows, "Study II synthetic-data generation")
  incoming <- mixedgp_bind_manifest_rows(rows)
  manifest <- mixedgp_merge_manifests(
    mixedgp_existing_manifest(directory), incoming
  )
  mixedgp_write_manifest_atomic(manifest, directory)
  manifest
}
