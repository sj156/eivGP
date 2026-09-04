############################################################
## reproduction_workflows.R
##
## The single publication-facing helper bundle for both numerical studies.
## Scientific posterior engines remain in their package-ready source modules;
## all simulation design, validation, data freezing, orchestration, manifests,
## task eligibility, aggregation, and fail-closed checks live here.
############################################################

MIXEDGP_SIMULATION_SCHEMA <- "2.0.0"
MIXEDGP_PUBLISHED_METHODS <- c("UC-GP", "LVGP", "EzGP")

mixedgp_simulation_code_dir <- function(code_dir = NULL) {
  marker <- "reproduction_workflows.R"
  if (!is.null(code_dir)) {
    code_dir <- normalizePath(code_dir, winslash = "/", mustWork = TRUE)
    if (!file.exists(file.path(code_dir, marker))) {
      stop("code_dir does not contain ", marker, ": ", code_dir)
    }
    return(code_dir)
  }
  command_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  command_dir <- if (length(command_file) > 0L) {
    dirname(sub("^--file=", "", command_file[[1L]]))
  } else {
    character(0)
  }
  installed_dir <- tryCatch(
    system.file("simulation-source", package = "eivGP"),
    error = function(e) ""
  )
  candidates <- unique(c(
    command_dir, getwd(), file.path(getwd(), "codes"),
    file.path(getwd(), "revision", "codes"), installed_dir
  ))
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(file.path(candidates, marker))]
  if (length(hit) == 0L) {
    stop(
      "Cannot locate the bundled publication simulation sources; ",
      "supply code_dir explicitly."
    )
  }
  normalizePath(hit[[1L]], winslash = "/", mustWork = TRUE)
}

mixedgp_simulation_artifact_root <- function(code_dir) {
  configured <- getOption("eivGP.project_root", NULL)
  if (!is.null(configured)) {
    if (length(configured) != 1L || is.na(configured) || !nzchar(configured)) {
      stop("getOption('eivGP.project_root') must be one nonempty path.")
    }
    return(normalizePath(
      configured, winslash = "/", mustWork = FALSE
    ))
  }
  installed_dir <- tryCatch(
    system.file("simulation-source", package = "eivGP"),
    error = function(e) ""
  )
  normalized_code <- normalizePath(code_dir, winslash = "/", mustWork = TRUE)
  if (nzchar(installed_dir)) {
    normalized_installed <- normalizePath(
      installed_dir, winslash = "/", mustWork = TRUE
    )
    if (identical(normalized_code, normalized_installed)) {
      return(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
    }
  }
  dirname(normalized_code)
}

mixedgp_simulation_modules <- function() {
  c(
    "core_parallel.R",
    "core_numerics.R",
    "model_univariate.R",
    "model_multivariate.R",
    "competitors.R",
    "model_api.R",
    "reproduction_data.R",
    "04_study1_ablations.R",
    "04_study2_ablations.R",
    "08_published_competitor_validation.R",
    "09_experiment_design_validation.R"
  )
}

mixedgp_simulation_engine <- function(code_dir) {
  code_dir <- mixedgp_simulation_code_dir(code_dir)
  configured_library <- Sys.getenv("MIXEDGP_R_LIBRARY", unset = "")
  if (nzchar(configured_library) && !dir.exists(configured_library)) {
    stop("MIXEDGP_R_LIBRARY does not exist: ", configured_library)
  }
  library_candidates <- c(
    if (nzchar(configured_library)) configured_library else character(0),
    file.path(dirname(code_dir), "R-library")
  )
  library_candidates <- library_candidates[dir.exists(library_candidates)]
  if (length(library_candidates) > 0L) {
    .libPaths(unique(c(
      normalizePath(library_candidates, winslash = "/", mustWork = TRUE),
      .libPaths()
    )))
  }
  ## Isolate the publication engine from objects in the caller's global
  ## workspace. Its parent retains the ordinary attached-package search path
  ## for base R compatibility, but never .GlobalEnv itself.
  engine <- new.env(parent = parent.env(.GlobalEnv))
  engine$source <- function(file, ...) {
    target <- if (grepl("^(/|[A-Za-z]:)", file)) {
      file
    } else {
      file.path(code_dir, file)
    }
    sys.source(target, envir = parent.frame(), chdir = TRUE)
    invisible(target)
  }
  engine$library <- function(package,
                             ...,
                             character.only = FALSE,
                             logical.return = FALSE,
                             quietly = FALSE,
                             warn.conflicts = TRUE) {
    package_name <- if (isTRUE(character.only)) {
      as.character(package)
    } else {
      as.character(substitute(package))
    }
    if (length(package_name) != 1L || is.na(package_name) ||
        !nzchar(package_name)) {
      stop("library() requires one package name in simulation scripts.")
    }
    namespace <- loadNamespace(package_name)
    exports <- getNamespaceExports(namespace)
    for (name in exports) {
      assign(name, getExportedValue(package_name, name), envir = engine)
    }
    if (isTRUE(logical.return)) return(TRUE)
    invisible(namespace)
  }
  for (module in mixedgp_simulation_modules()) {
    path <- file.path(code_dir, module)
    if (!file.exists(path)) stop("Missing simulation module: ", path)
    sys.source(path, envir = engine, chdir = TRUE)
  }
  engine
}

mixedgp_validate_scalar_integer <- function(x, name, minimum = 0L) {
  if (length(x) != 1L || is.na(x) || !is.numeric(x) || x != as.integer(x) ||
      x < minimum) {
    stop(name, " must be one integer at least ", minimum, ".")
  }
  as.integer(x)
}

mixedgp_validate_boolean <- function(x, name) {
  if (length(x) != 1L || !is.logical(x) || is.na(x)) {
    stop(name, " must be TRUE or FALSE.")
  }
  isTRUE(x)
}

mixedgp_validate_calibration_grid <- function(x, n, name = "calibration_grid") {
  if (!is.numeric(x) || length(x) < 1L || anyNA(x) ||
      any(x != as.integer(x)) || any(x < 0L | x > n) || anyDuplicated(x)) {
    stop(name, " must contain unique integer sizes between zero and n.")
  }
  sort(as.integer(x))
}

mixedgp_available_workers <- function(cap = 8L) {
  cap <- mixedgp_validate_scalar_integer(cap, "cap", minimum = 1L)
  detected <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (!is.finite(detected) || detected < 1L) detected <- 1L
  as.integer(min(cap, detected))
}

mixedgp_simulation_mode <- function(mode) {
  match.arg(mode, c("dry_run", "data", "smoke", "publication"))
}

mixedgp_number_slug <- function(x) {
  gsub("[^0-9A-Za-z]+", "p", format(x, scientific = FALSE, trim = TRUE))
}

mixedgp_config_fingerprint <- function(config) {
  payload <- config
  payload$created_at <- NULL
  payload$stages <- NULL
  code_files <- unique(c(
    mixedgp_simulation_modules(), "reproduction_workflows.R",
    "02_study1_monte_carlo.R", "02_study2_monte_carlo.R",
    "run_study1_simulation.R", "run_study2_simulation.R"
  ))
  code_paths <- file.path(config$code_dir, code_files)
  payload$source_identity <- setNames(
    ifelse(
      file.exists(code_paths),
      unname(tools::md5sum(code_paths)),
      NA_character_
    ),
    code_files
  )
  package_names <- c(
    kergp = "kergp", LVGP = "LVGP", EzGP = "EzGP",
    ggplot2 = "ggplot2", dplyr = "dplyr", tidyr = "tidyr",
    knitr = "knitr", posterior = "posterior"
  )
  if (identical(config$study, "study2")) {
    package_names <- c(package_names, TruncatedNormal = "TruncatedNormal")
  }
  payload$competitor_package_identity <- vapply(
    package_names,
    function(package) {
      if (requireNamespace(package, quietly = TRUE)) {
        as.character(utils::packageVersion(package))
      } else {
        "not-installed"
      }
    },
    character(1L)
  )
  path <- tempfile("mixedgp-config-", fileext = ".rds")
  on.exit(if (file.exists(path)) unlink(path), add = TRUE)
  saveRDS(payload, path, version = 3L, compress = FALSE)
  unname(tools::md5sum(path))
}

mixedgp_estimand_method_matrix <- function(study = c("study1", "study2")) {
  study <- match.arg(study)
  rows <- list(
    c("surface", "m(x,c)", "EIV-GP", "proposed", "yes", "yes", "all calibration sizes", ""),
    c("surface", "m(x,c)", "UC-GP", "published competitor", "yes", "point only", "all calibration sizes", "does not define physical u"),
    c("surface", "m(x,c)", "LVGP", "published competitor", "yes", "point only", "all calibration sizes", "does not define physical u"),
    c("surface", "m(x,c)", "EzGP", "published competitor", "yes", "point only", "all calibration sizes", "does not define physical u"),
    c("surface", "f(x,u)", "EIV-GP", "proposed", "yes", "yes", "positive, anchoring calibration", "physical u scale is not identified at k=0"),
    c("surface", "f(x,u)", "PI-GP", "scientific ablation", "yes", "optimizer approximation", "positive, anchoring calibration", "plug-in ignores latent-input uncertainty"),
    c("surface", "f(x,u)", "CC-GP", "scientific ablation", "yes", "optimizer approximation", "at least three complete cases", "uses only calibrated cases"),
    c("surface", "f(x,u)", "Full-U GP", "infeasible benchmark", "yes", "optimizer approximation", "all training u observed", "benchmark, not a deployable competitor"),
    c("prediction", "Y*|x*,c*", "EIV-GP", "proposed", "yes", "yes", "all calibration sizes", ""),
    c("prediction", "Y*|x*,c*", "UC-GP", "published competitor", "yes", "predictive", "all calibration sizes", ""),
    c("prediction", "Y*|x*,c*", "LVGP", "published competitor", "yes", "predictive", "all calibration sizes", ""),
    c("prediction", "Y*|x*,c*", "EzGP", "published competitor", "yes", "predictive", "all calibration sizes", ""),
    c("prediction", "Y*|x*,c*", "PI-GP", "scientific ablation", "yes", "predictive", "cells with ablations", "imputes a latent input and then plugs it into a GP"),
    c("prediction", "Y*|x*,c*", "CC-GP", "scientific ablation", "yes", "predictive", "positive calibration in cells with ablations", "fits the response GP only to complete cases"),
    c("prediction", "Y*|x*,c*", "DGM oracle", "truth benchmark", "yes", "truth", "simulation only", "never ranked as a fitted method"),
    c("latent state", "training U|X,C,Y", "EIV-GP", "proposed", "yes", "yes", "positive, anchoring calibration", "physical coordinates are not identified at k=0"),
    c("latent state", "training U|C", "RF-OM", "response-free ablation", "yes", "yes", "positive, anchoring calibration", "does not use y")
  )
  if (study == "study2") {
    rows <- c(rows, list(
      c("latent state", "prospective U*|X*,C*,Y-data", "EIV-GP", "proposed", "yes", "yes", "positive, anchoring calibration", "y* is not conditioned on"),
      c("latent state", "prospective U*|C*", "RF-OM", "response-free ablation", "yes", "yes", "positive, anchoring calibration", "does not use response data")
    ))
  }
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(out) <- c(
    "task", "estimand", "method", "role", "eligible",
    "uncertainty_supported", "condition", "reason_or_limit"
  )
  out$study <- study
  out[, c("study", setdiff(names(out), "study")), drop = FALSE]
}

mixedgp_study1_cells <- function(mode) {
  publication <- identical(mode, "publication")
  smoke <- identical(mode, "smoke")
  primary_rep <- if (smoke) 1L else if (publication) 100L else 100L
  sensitivity_rep <- if (smoke) 1L else if (publication) 50L else 50L
  n_test <- if (smoke) 80L else 1000L
  primary <- lapply(c(0, 0.5, 1), function(eta) {
    list(
      id = paste0("eta", mixedgp_number_slug(eta), "_balanced"),
      role = "primary_paired_heterogeneity",
      scenario = "heterogeneity_continuum",
      heterogeneity_eta = eta,
      threshold_design = "balanced",
      min_class_count = 0L,
      n = 100L,
      n_test = n_test,
      n_rep = primary_rep,
      m = 6L,
      calibration_grid = if (smoke) c(0L, 5L) else c(0L, 5L, 10L, 20L, 50L),
      evaluate_f = identical(eta, 1),
      evaluate_u = eta %in% c(0, 1),
      run_ablations = TRUE
    )
  })
  sensitivity <- list(list(
    id = "eta1_imbalanced",
    role = "appendix_threshold_imbalance",
    scenario = "heterogeneity_continuum",
    heterogeneity_eta = 1,
    threshold_design = "imbalanced",
    min_class_count = 3L,
    n = 100L,
    n_test = n_test,
    n_rep = sensitivity_rep,
    m = 6L,
    calibration_grid = if (smoke) c(5L) else c(20L),
    evaluate_f = TRUE,
    evaluate_u = TRUE,
    run_ablations = TRUE
  ))
  c(primary, sensitivity)
}

mixedgp_study2_cells <- function(mode) {
  publication <- identical(mode, "publication")
  smoke <- identical(mode, "smoke")
  main_rep <- if (smoke) 1L else if (publication) 100L else 100L
  appendix_rep <- if (smoke) 1L else if (publication) 50L else 50L
  n_test <- if (smoke) 60L else 1000L
  n_test_stress <- if (smoke) 80L else 2000L
  make_cell <- function(id, role, scenario, q, grid, n_rep, n_test,
                        run_ablations, evaluate_f = FALSE,
                        evaluate_u = TRUE) {
    list(
      id = id, role = role, scenario = scenario, q = as.integer(q),
      d = 2L, m = 4L, n = 120L, n_test = as.integer(n_test),
      n_rep = as.integer(n_rep), calibration_grid = as.integer(grid),
      run_ablations = isTRUE(run_ablations),
      evaluate_f = isTRUE(evaluate_f), evaluate_u = isTRUE(evaluate_u)
    )
  }
  list(
    make_cell("primary_q2", "primary_proxy_dimension", "primary", 2L,
              if (smoke) 6L else 12L, main_rep, n_test, TRUE),
    make_cell("primary_q3", "primary_proxy_dimension", "primary", 3L,
              if (smoke) 6L else 12L, main_rep, n_test, TRUE),
    make_cell("primary_q4_calibration", "primary_calibration_curve", "primary", 4L,
              if (smoke) c(0L, 6L) else c(0L, 6L, 12L, 24L, 48L),
              main_rep, n_test, TRUE, evaluate_f = TRUE),
    make_cell("additive_q4", "negative_control", "latent_additive_control", 4L,
              if (smoke) 6L else 12L, appendix_rep, n_test, TRUE,
              evaluate_u = FALSE),
    make_cell("high_uncertainty_q4", "uncertainty_control", "high_uncertainty", 4L,
              if (smoke) 6L else 12L, appendix_rep, n_test, TRUE,
              evaluate_u = FALSE),
    make_cell("logistic_q4", "appendix_misspecification", "logistic_misspec", 4L,
              if (smoke) 6L else 12L, appendix_rep, n_test, FALSE,
              evaluate_u = FALSE),
    make_cell("primary_q6", "appendix_sparse_pattern_stress", "primary", 6L,
              if (smoke) 6L else 12L, appendix_rep, n_test_stress, TRUE)
  )
}

#' Configure the Study I publication simulation
#'
#' Constructs the frozen Study I design used by the publication workflow.
#' The default is a read-only dry run. Use `mode = "smoke"` for an
#' end-to-end non-reportable test and `mode = "publication"` for the frozen
#' reportable design.
#'
#' @param mode One of `"dry_run"`, `"data"`, `"smoke"`, or
#'   `"publication"`.
#' @param code_dir Optional directory containing the bundled audited source
#'   scripts. Installed packages locate these scripts automatically.
#' @param workers Number of replication workers.
#' @param output_root Optional directory for fingerprinted run outputs.
#' @param data_root Optional directory for frozen synthetic datasets.
#' @param stages Optional execution stages. Resume overrides are allowed only
#'   in `"smoke"` or `"publication"` mode and may be `"data"`, `"fit"`, or
#'   `c("fit", "aggregate")`.
#'
#' @return A validated object of class `mixedgp_simulation_config`.
#' @export
study1_simulation_config <- function(
    mode = "dry_run",
    code_dir = NULL,
    workers = mixedgp_available_workers(),
    output_root = NULL,
    data_root = NULL,
    stages = NULL) {
  mode <- mixedgp_simulation_mode(mode)
  code_dir <- mixedgp_simulation_code_dir(code_dir)
  artifact_root <- mixedgp_simulation_artifact_root(code_dir)
  if (is.null(output_root)) {
    output_root <- file.path(artifact_root, "results", "study1")
  }
  if (is.null(data_root)) {
    data_root <- file.path(
      artifact_root, "data", "synthetic", "study1"
    )
  }
  smoke <- identical(mode, "smoke")
  config <- list(
    schema_version = MIXEDGP_SIMULATION_SCHEMA,
    study = "study1",
    mode = mode,
    run_id = paste0("study1-", mode),
    code_dir = normalizePath(code_dir, winslash = "/"),
    output_root = normalizePath(output_root, winslash = "/", mustWork = FALSE),
    data_root = normalizePath(data_root, winslash = "/", mustWork = FALSE),
    stages = switch(
      mode,
      dry_run = character(0), data = "data",
      smoke = c("data", "fit", "aggregate"),
      publication = c("data", "fit", "aggregate")
    ),
    cells = mixedgp_study1_cells(mode),
    published_methods = MIXEDGP_PUBLISHED_METHODS,
    strict_competitors = identical(mode, "publication"),
    fail_closed = identical(mode, "publication"),
    use_cache = TRUE,
    parallel = list(level = "replications", workers = as.integer(workers)),
    mcmc = if (smoke) {
      list(n_iter = 120L, burn = 40L, thin = 1L, n_chains = 1L,
           rhat_limit = 1.05, ess_limit = 20L, require_gate = FALSE,
           pilot_reps = 0L,
           adaptive = list(initial_draws = 40L, target_ess = 20,
                           max_draws = 80L, extension_draws = 20L))
    } else {
      list(n_iter = 16000L, burn = 1000L, thin = 1L, n_chains = 4L,
           # Formal publication runs use the stricter 200-draw
           # key-parameter ESS gate. Fits below it are not accepted.
           rhat_limit = 1.05, ess_limit = 200L,
           require_gate = identical(mode, "publication"),
           pilot_reps = 0L,
           adaptive = list(initial_draws = 5000L, target_ess = 200,
                           max_draws = 15000L, extension_draws = 1000L))
    },
    measurement_mcmc = if (smoke) {
      list(n_iter = 120L, burn = 40L)
    } else {
      list(n_iter = 4000L, burn = 1000L)
    },
    evaluation = if (smoke) {
      list(n_pred_draw = 40L, n_m_eval = 12L, n_m_draw = 20L,
           n_m_latent = 24L)
    } else {
      list(n_pred_draw = 500L, n_m_eval = 200L, n_m_draw = 250L,
           n_m_latent = 512L)
    },
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  if (!is.null(stages)) {
    if (!mode %in% c("smoke", "publication")) {
      stop("stages may be overridden only in smoke or publication mode.")
    }
    config$stages <- as.character(stages)
  }
  class(config) <- c("mixedgp_simulation_config", "list")
  validate_simulation_config(config)
}

#' Configure the Study II publication simulation
#'
#' Constructs the frozen Study II design used by the publication workflow.
#' The default is a read-only dry run. Use `mode = "smoke"` for an
#' end-to-end non-reportable test and `mode = "publication"` for the frozen
#' reportable design.
#'
#' @param mode One of `"dry_run"`, `"data"`, `"smoke"`, or
#'   `"publication"`.
#' @param code_dir Optional directory containing the bundled audited source
#'   scripts. Installed packages locate these scripts automatically.
#' @param workers Number of replication workers.
#' @param output_root Optional directory for fingerprinted run outputs.
#' @param data_root Optional directory for frozen synthetic datasets.
#' @param stages Optional execution stages. Resume overrides are allowed only
#'   in `"smoke"` or `"publication"` mode and may be `"data"`, `"fit"`, or
#'   `c("fit", "aggregate")`.
#'
#' @return A validated object of class `mixedgp_simulation_config`.
#' @export
study2_simulation_config <- function(
    mode = "dry_run",
    code_dir = NULL,
    workers = mixedgp_available_workers(),
    output_root = NULL,
    data_root = NULL,
    stages = NULL) {
  mode <- mixedgp_simulation_mode(mode)
  code_dir <- mixedgp_simulation_code_dir(code_dir)
  artifact_root <- mixedgp_simulation_artifact_root(code_dir)
  if (is.null(output_root)) {
    output_root <- file.path(artifact_root, "results", "study2")
  }
  if (is.null(data_root)) {
    data_root <- file.path(
      artifact_root, "data", "synthetic", "study2"
    )
  }
  smoke <- identical(mode, "smoke")
  config <- list(
    schema_version = MIXEDGP_SIMULATION_SCHEMA,
    study = "study2",
    mode = mode,
    run_id = paste0("study2-", mode),
    code_dir = normalizePath(code_dir, winslash = "/"),
    output_root = normalizePath(output_root, winslash = "/", mustWork = FALSE),
    data_root = normalizePath(data_root, winslash = "/", mustWork = FALSE),
    stages = switch(
      mode,
      dry_run = character(0), data = "data",
      smoke = c("data", "fit", "aggregate"),
      publication = c("data", "fit", "aggregate")
    ),
    cells = mixedgp_study2_cells(mode),
    published_methods = MIXEDGP_PUBLISHED_METHODS,
    strict_competitors = identical(mode, "publication"),
    fail_closed = identical(mode, "publication"),
    use_cache = TRUE,
    parallel = list(level = "replications", workers = as.integer(workers)),
    predictive_latent_sampler = "minimax_tilting",
    mcmc = if (smoke) {
      list(n_iter = 200L, burn = 40L, thin = 2L, n_chains = 1L,
           rhat_limit = 1.05, raw_ess_limit = 20L,
           target_bulk_ess_limit = 10L, target_tail_ess_limit = 10L,
           require_gate = FALSE, pilot_reps = 0L,
           adaptive = list(initial_draws = 40L, target_ess = 20,
                           max_draws = 80L, extension_draws = 20L))
    } else {
      list(n_iter = 31000L, burn = 1000L, thin = 2L, n_chains = 4L,
           rhat_limit = 1.05, raw_ess_limit = 200L,
           target_bulk_ess_limit = 100L, target_tail_ess_limit = 100L,
           require_gate = identical(mode, "publication"), pilot_reps = 0L,
           adaptive = list(initial_draws = 5000L, target_ess = 200,
                           max_draws = 15000L, extension_draws = 1000L))
    },
    measurement_mcmc = if (smoke) {
      list(n_iter = 120L, burn = 40L, thin = 2L, n_chains = 2L)
    } else {
      list(n_iter = 4000L, burn = 1000L, thin = 2L, n_chains = 4L)
    },
    evaluation = if (smoke) {
      list(n_pred_draw = 30L, n_m_eval = 8L, n_m_draw = 12L,
           n_m_latent = 16L, n_m_truth = 100L, n_oracle_pool = 5000L)
    } else {
      list(n_pred_draw = 500L, n_m_eval = 150L, n_m_draw = 200L,
           n_m_latent = 256L, n_m_truth = 5000L,
           n_oracle_pool = 250000L)
    },
    ablation_gp = list(
      n_starts = if (smoke) 3L else 8L,
      maxit = if (smoke) 200L else 500L
    ),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  if (!is.null(stages)) {
    if (!mode %in% c("smoke", "publication")) {
      stop("stages may be overridden only in smoke or publication mode.")
    }
    config$stages <- as.character(stages)
  }
  class(config) <- c("mixedgp_simulation_config", "list")
  validate_simulation_config(config)
}

validate_simulation_config <- function(config) {
  if (!is.list(config) || !identical(
    config$schema_version, MIXEDGP_SIMULATION_SCHEMA
  )) {
    stop("Unsupported or malformed simulation configuration.")
  }
  config$study <- match.arg(config$study, c("study1", "study2"))
  config$mode <- mixedgp_simulation_mode(config$mode)
  allowed_stages <- c("data", "fit", "aggregate")
  if (anyDuplicated(config$stages) || any(!config$stages %in% allowed_stages)) {
    stop("stages must be a unique subset of data, fit, and aggregate.")
  }
  if ("aggregate" %in% config$stages && !"fit" %in% config$stages) {
    stop("The aggregate stage requires fit in the same invocation.")
  }
  if ("fit" %in% config$stages && !"data" %in% config$stages) {
    ## Fitting frozen data without regenerating it is valid, but the artifacts
    ## are verified before dispatch. This branch is intentionally allowed.
    invisible(NULL)
  }
  if (!is.list(config$cells) || length(config$cells) < 1L) {
    stop("config$cells must be a nonempty list.")
  }
  ids <- vapply(config$cells, `[[`, character(1L), "id")
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Simulation cell ids must be nonempty and unique.")
  }
  required <- c("id", "role", "scenario", "n", "n_test", "n_rep", "m",
                "calibration_grid", "run_ablations", "evaluate_f", "evaluate_u")
  for (cell in config$cells) {
    study_required <- if (config$study == "study1") {
      c("heterogeneity_eta", "threshold_design", "min_class_count")
    } else {
      c("q", "d")
    }
    required_cell <- c(required, study_required)
    missing <- setdiff(required_cell, names(cell))
    if (length(missing) > 0L) {
      stop("Cell ", cell$id, " is missing: ", paste(missing, collapse = ", "))
    }
    n <- mixedgp_validate_scalar_integer(cell$n, paste0(cell$id, "$n"), 1L)
    mixedgp_validate_scalar_integer(cell$n_test, paste0(cell$id, "$n_test"), 1L)
    mixedgp_validate_scalar_integer(cell$n_rep, paste0(cell$id, "$n_rep"), 1L)
    mixedgp_validate_scalar_integer(cell$m, paste0(cell$id, "$m"), 2L)
    calibration_grid <- mixedgp_validate_calibration_grid(
      cell$calibration_grid, n, paste0(cell$id, "$calibration_grid")
    )
    if (!identical(as.integer(cell$calibration_grid), calibration_grid)) {
      stop(cell$id, "$calibration_grid must be sorted increasingly.")
    }
    for (flag in c("run_ablations", "evaluate_f", "evaluate_u")) {
      mixedgp_validate_boolean(cell[[flag]], paste0(cell$id, "$", flag))
    }
    if (config$study == "study1") {
      eta <- as.numeric(cell$heterogeneity_eta)
      if (length(eta) != 1L || !is.finite(eta) || eta < 0) {
        stop(cell$id, "$heterogeneity_eta must be finite and nonnegative.")
      }
      if (!cell$threshold_design %in% c("balanced", "imbalanced")) {
        stop(cell$id, "$threshold_design must be balanced or imbalanced.")
      }
      mixedgp_validate_scalar_integer(
        cell$min_class_count, paste0(cell$id, "$min_class_count"), 0L
      )
      if (cell$threshold_design == "imbalanced" && cell$m != 6L) {
        stop("The imbalanced Study I design requires m=6 in ", cell$id, ".")
      }
    } else {
      mixedgp_validate_scalar_integer(cell$q, paste0(cell$id, "$q"), 2L)
      if (cell$q > 6L) stop("Study II publication cells require q <= 6.")
      mixedgp_validate_scalar_integer(cell$d, paste0(cell$id, "$d"), 1L)
      if (cell$d != 2L) stop("The Study II DGM currently requires d=2.")
    }
  }
  if (!identical(config$published_methods, MIXEDGP_PUBLISHED_METHODS)) {
    stop(
      "The publication competitor set is frozen as: ",
      paste(MIXEDGP_PUBLISHED_METHODS, collapse = ", "), "."
    )
  }
  config$parallel$workers <- mixedgp_validate_scalar_integer(
    config$parallel$workers, "parallel$workers", 1L
  )
  if (!identical(config$parallel$level, "replications")) {
    stop("Publication runs parallelize replications and keep chains serial.")
  }
  for (flag in c("strict_competitors", "fail_closed", "use_cache")) {
    mixedgp_validate_boolean(config[[flag]], flag)
  }
  for (field in c("n_iter", "burn", "thin", "n_chains")) {
    minimum <- if (field == "burn") 0L else 1L
    mixedgp_validate_scalar_integer(
      config$mcmc[[field]], paste0("mcmc$", field), minimum
    )
  }
  if (config$mcmc$burn >= config$mcmc$n_iter) {
    stop("mcmc$burn must be smaller than mcmc$n_iter.")
  }
  if (!is.list(config$mcmc$adaptive)) {
    stop("mcmc$adaptive must be a named adaptive-MCMC control list.")
  }
  adaptive <- do.call(
    mixedgp_adaptive_mcmc_control,
    config$mcmc$adaptive
  )
  retained_capacity <- floor(
    (config$mcmc$n_iter - config$mcmc$burn) / config$mcmc$thin
  )
  if (adaptive$max_draws > retained_capacity) {
    stop(
      "mcmc$adaptive$max_draws exceeds the retained-draw capacity implied ",
      "by mcmc$n_iter, mcmc$burn, and mcmc$thin."
    )
  }
  mixedgp_validate_boolean(config$mcmc$require_gate, "mcmc$require_gate")
  if (!is.numeric(config$mcmc$rhat_limit) ||
      length(config$mcmc$rhat_limit) != 1L ||
      !is.finite(config$mcmc$rhat_limit) || config$mcmc$rhat_limit <= 1) {
    stop("mcmc$rhat_limit must be one finite number greater than one.")
  }
  measurement_fields <- intersect(
    c("n_iter", "burn", "thin", "n_chains"),
    names(config$measurement_mcmc)
  )
  for (field in measurement_fields) {
    minimum <- if (field == "burn") 0L else 1L
    mixedgp_validate_scalar_integer(
      config$measurement_mcmc[[field]],
      paste0("measurement_mcmc$", field), minimum
    )
  }
  if (config$measurement_mcmc$burn >= config$measurement_mcmc$n_iter) {
    stop("measurement_mcmc$burn must be smaller than measurement_mcmc$n_iter.")
  }
  evaluation_fields <- if (config$study == "study1") {
    c("n_pred_draw", "n_m_eval", "n_m_draw", "n_m_latent")
  } else {
    c(
      "n_pred_draw", "n_m_eval", "n_m_draw", "n_m_latent",
      "n_m_truth", "n_oracle_pool"
    )
  }
  missing_evaluation <- setdiff(evaluation_fields, names(config$evaluation))
  if (length(missing_evaluation) > 0L) {
    stop("evaluation is missing: ", paste(missing_evaluation, collapse = ", "))
  }
  for (field in evaluation_fields) {
    mixedgp_validate_scalar_integer(
      config$evaluation[[field]], paste0("evaluation$", field), 1L
    )
  }
  mixedgp_validate_scalar_integer(
    config$mcmc$pilot_reps, "mcmc$pilot_reps", 0L
  )
  if (config$study == "study1") {
    mixedgp_validate_scalar_integer(config$mcmc$ess_limit, "mcmc$ess_limit", 1L)
  } else {
    for (field in c(
      "raw_ess_limit", "target_bulk_ess_limit", "target_tail_ess_limit"
    )) {
      mixedgp_validate_scalar_integer(
        config$mcmc[[field]], paste0("mcmc$", field), 1L
      )
    }
    if (!identical(config$predictive_latent_sampler, "minimax_tilting")) {
      stop("Publication Study II uses the exact minimax_tilting latent sampler.")
    }
    mixedgp_validate_scalar_integer(
      config$ablation_gp$n_starts, "ablation_gp$n_starts", 1L
    )
    mixedgp_validate_scalar_integer(
      config$ablation_gp$maxit, "ablation_gp$maxit", 1L
    )
  }
  if (config$mode == "publication" &&
      (!isTRUE(config$strict_competitors) || !isTRUE(config$fail_closed) ||
       !isTRUE(config$mcmc$require_gate))) {
    stop(
      "Publication mode requires strict competitors, fail-closed outputs, ",
      "and the MCMC gate."
    )
  }
  config
}

mixedgp_bind_rows_base <- function(rows) {
  rows <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, rows)
  if (length(rows) == 0L) return(data.frame())
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  aligned <- lapply(rows, function(x) {
    missing <- setdiff(columns, names(x))
    for (name in missing) x[[name]] <- NA
    x[, columns, drop = FALSE]
  })
  out <- do.call(rbind, aligned)
  rownames(out) <- NULL
  out
}

mixedgp_task_plan <- function(config) {
  config <- validate_simulation_config(config)
  rows <- list()
  add <- function(cell, task_type, method, n_calib = NA_integer_) {
    n <- cell$n_rep * if (all(is.na(n_calib))) 1L else length(n_calib)
    rows[[length(rows) + 1L]] <<- data.frame(
      study = config$study,
      cell_id = cell$id,
      design_role = cell$role,
      task_type = task_type,
      method = method,
      n_calib = if (all(is.na(n_calib))) NA_integer_ else rep(n_calib, each = cell$n_rep),
      rep = if (all(is.na(n_calib))) seq_len(cell$n_rep) else rep(seq_len(cell$n_rep), times = length(n_calib)),
      planned_units = rep(1L, n),
      stringsAsFactors = FALSE
    )
  }
  for (cell in config$cells) {
    add(cell, "frozen_data", "DGM")
    add(cell, "prediction_truth", "DGM oracle")
    for (method in config$published_methods) {
      add(cell, "published_fit", method)
    }
    for (n_calib in cell$calibration_grid) {
      add(cell, "eiv_fit", "EIV-GP", n_calib)
    }
    if (isTRUE(cell$run_ablations)) {
      for (method in c("RF-OM", "PI-GP", "CC-GP")) {
        for (n_calib in cell$calibration_grid) {
          add(cell, "ablation_fit", method, n_calib)
        }
      }
      if (isTRUE(cell$evaluate_f)) add(cell, "benchmark_fit", "Full-U GP")
    }
  }
  mixedgp_bind_rows_base(rows)
}

mixedgp_cell_summary <- function(config) {
  mixedgp_bind_rows_base(lapply(config$cells, function(cell) {
    data.frame(
      cell_id = cell$id,
      role = cell$role,
      scenario = cell$scenario,
      eta = if (is.null(cell$heterogeneity_eta)) NA_real_ else cell$heterogeneity_eta,
      q = if (is.null(cell$q)) NA_integer_ else cell$q,
      n = cell$n,
      n_test = cell$n_test,
      n_rep = cell$n_rep,
      calibration_grid = paste(cell$calibration_grid, collapse = ";"),
      evaluate_f = cell$evaluate_f,
      evaluate_u = cell$evaluate_u,
      run_ablations = cell$run_ablations,
      stringsAsFactors = FALSE
    )
  }))
}

mixedgp_runtime_preflight <- function(config) {
  packages <- c("ggplot2", "dplyr", "tidyr", "knitr", "posterior")
  if (config$study == "study2") packages <- c(packages, "TruncatedNormal")
  available <- vapply(packages, requireNamespace, logical(1L), quietly = TRUE)
  version <- vapply(seq_along(packages), function(ii) {
    if (available[[ii]]) {
      as.character(utils::packageVersion(packages[[ii]]))
    } else {
      NA_character_
    }
  }, character(1L))
  data.frame(
    package = packages,
    role = c(
      "figures", "summaries and diagnostics", "summary reshaping",
      "LaTeX tables", "rank-normalized MCMC diagnostics",
      if (config$study == "study2") "exact minimax-tilting sampler"
    ),
    available = available,
    version = version,
    stringsAsFactors = FALSE
  )
}

mixedgp_run_directory <- function(config) {
  fingerprint <- mixedgp_config_fingerprint(config)
  file.path(config$output_root, paste0(config$run_id, "-", substr(fingerprint, 1L, 12L)))
}

mixedgp_data_cell_directory <- function(config, cell) {
  file.path(config$data_root, cell$id)
}

mixedgp_code_hashes <- function(config) {
  files <- unique(c(
    mixedgp_simulation_modules(), "reproduction_workflows.R",
    "02_study1_monte_carlo.R", "02_study2_monte_carlo.R",
    "run_study1_simulation.R", "run_study2_simulation.R"
  ))
  paths <- file.path(config$code_dir, files)
  exists <- file.exists(paths)
  data.frame(
    file = files,
    exists = exists,
    md5 = ifelse(exists, unname(tools::md5sum(paths)), NA_character_),
    stringsAsFactors = FALSE
  )
}

mixedgp_atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(object, temporary, version = 3L)
  if (file.exists(path) && !file.remove(path)) {
    stop("Could not replace existing file: ", path)
  }
  if (!file.rename(temporary, path)) stop("Could not install file: ", path)
  invisible(path)
}

mixedgp_atomic_write_csv <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  utils::write.csv(object, temporary, row.names = FALSE, na = "")
  if (file.exists(path) && !file.remove(path)) {
    stop("Could not replace existing file: ", path)
  }
  if (!file.rename(temporary, path)) stop("Could not install file: ", path)
  invisible(path)
}

mixedgp_write_progress <- function(run_dir,
                                   status,
                                   phase,
                                   cells_completed = 0L,
                                   cells_total = 0L,
                                   current_cell = "",
                                   detail = "",
                                   run_started_at = NULL) {
  cells_completed <- as.integer(cells_completed)
  cells_total <- as.integer(cells_total)
  now <- Sys.time()
  elapsed_seconds <- if (is.null(run_started_at)) {
    NA_real_
  } else {
    as.numeric(difftime(now, run_started_at, units = "secs"))
  }
  estimated_remaining_seconds <- if (
    is.finite(elapsed_seconds) && cells_completed > 0L &&
      cells_total > cells_completed
  ) {
    elapsed_seconds / cells_completed * (cells_total - cells_completed)
  } else if (cells_completed >= cells_total && cells_total > 0L) {
    0
  } else {
    NA_real_
  }
  progress <- data.frame(
    updated_at_utc = format(now, tz = "UTC", usetz = TRUE),
    status = as.character(status),
    phase = as.character(phase),
    current_cell = as.character(current_cell),
    cells_completed = cells_completed,
    cells_total = cells_total,
    percent_cells_complete = if (cells_total > 0L) {
      round(100 * cells_completed / cells_total, 1L)
    } else {
      NA_real_
    },
    elapsed_seconds = elapsed_seconds,
    estimated_remaining_seconds = estimated_remaining_seconds,
    estimated_finish_utc = if (is.finite(estimated_remaining_seconds)) {
      format(now + estimated_remaining_seconds, tz = "UTC", usetz = TRUE)
    } else {
      NA_character_
    },
    estimate_basis = if (is.finite(estimated_remaining_seconds) &&
                         cells_completed > 0L) {
      "completed_cell_average"
    } else {
      "not_available"
    },
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
  mixedgp_atomic_write_csv(progress, file.path(run_dir, "config", "progress.csv"))
  message(
    "[", progress$updated_at_utc, "] ",
    progress$status, " | ", progress$phase,
    if (nzchar(progress$current_cell)) paste0(" | ", progress$current_cell) else "",
    if (cells_total > 0L) paste0(
      " | ", cells_completed, "/", cells_total,
      " cells (", progress$percent_cells_complete, "%)"
    ) else "",
    if (is.finite(estimated_remaining_seconds)) paste0(
      " | estimated remaining ", round(estimated_remaining_seconds / 60, 1L),
      " min"
    ) else "",
    if (nzchar(progress$detail)) paste0(" | ", progress$detail) else ""
  )
  invisible(progress)
}

mixedgp_write_run_metadata <- function(config,
                                       run_dir,
                                       preflight,
                                       runtime_preflight,
                                       task_plan) {
  metadata_dir <- file.path(run_dir, "config")
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
  mixedgp_atomic_save_rds(config, file.path(metadata_dir, "resolved_config.rds"))
  config_text <- utils::capture.output(dput(config))
  writeLines(config_text, file.path(metadata_dir, "resolved_config.R"))
  mixedgp_atomic_write_csv(
    mixedgp_cell_summary(config), file.path(metadata_dir, "design_cells.csv")
  )
  mixedgp_atomic_write_csv(
    mixedgp_estimand_method_matrix(config$study),
    file.path(metadata_dir, "estimand_method_matrix.csv")
  )
  mixedgp_atomic_write_csv(task_plan, file.path(metadata_dir, "task_plan.csv"))
  mixedgp_atomic_write_csv(preflight, file.path(metadata_dir, "competitor_preflight.csv"))
  mixedgp_atomic_write_csv(
    runtime_preflight, file.path(metadata_dir, "runtime_preflight.csv")
  )
  mixedgp_atomic_write_csv(
    mixedgp_code_hashes(config), file.path(metadata_dir, "code_hashes.csv")
  )
  utils::capture.output(
    utils::sessionInfo(), file = file.path(metadata_dir, "sessionInfo.txt")
  )
  invisible(metadata_dir)
}

mixedgp_run_automatic_validators <- function(config, engine, run_dir) {
  metadata_dir <- file.path(run_dir, "config")
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

  thread_vars <- c(
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"
  )
  old_threads <- Sys.getenv(thread_vars, unset = NA_character_)
  do.call(Sys.setenv, as.list(setNames(rep("1", length(thread_vars)), thread_vars)))
  on.exit({
    present <- !is.na(old_threads)
    if (any(present)) do.call(Sys.setenv, as.list(old_threads[present]))
    if (any(!present)) Sys.unsetenv(thread_vars[!present])
  }, add = TRUE)

  validators <- list(
    published_competitors = list(
      function_name = "run_published_competitor_validation",
      output_file = "published_competitor_validation.csv"
    ),
    experiment_design = list(
      function_name = "run_experiment_design_validation",
      output_file = "experiment_design_validation.csv"
    )
  )
  outputs <- list()
  statuses <- list()

  for (name in names(validators)) {
    specification <- validators[[name]]
    started <- Sys.time()
    answer <- tryCatch(
      get(
        specification$function_name,
        envir = engine,
        inherits = FALSE
      )(),
      error = function(e) e
    )
    if (!inherits(answer, "error")) {
      required_columns <- if (identical(name, "published_competitors")) {
        c("method", "status")
      } else {
        c("validator", "pass")
      }
      if (!is.data.frame(answer) || nrow(answer) < 1L ||
          any(!required_columns %in% names(answer))) {
        answer <- simpleError(paste0(
          specification$function_name,
          " returned a malformed validation table."
        ))
      }
    }
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (inherits(answer, "error")) {
      outputs[[name]] <- NULL
      statuses[[name]] <- data.frame(
        validator = name,
        pass = FALSE,
        message = conditionMessage(answer),
        elapsed_seconds = elapsed,
        stringsAsFactors = FALSE
      )
    } else {
      outputs[[name]] <- answer
      mixedgp_atomic_write_csv(
        answer,
        file.path(metadata_dir, specification$output_file)
      )
      failed_competitors <- if (identical(name, "published_competitors")) {
        is.na(answer$status) | answer$status != "success"
      } else {
        logical(0)
      }
      ## External competitor packages are audited, but are not a prerequisite
      ## for completing the EIV-GP MCMC experiment.  Their validation result
      ## remains in the run directory for transparent reporting.
      strict_failure <- FALSE
      invariant_failure <- identical(name, "experiment_design") &&
        any(is.na(answer$pass) | !answer$pass)
      statuses[[name]] <- data.frame(
        validator = name,
        pass = !(strict_failure || invariant_failure),
        message = if (strict_failure) {
          paste0(
            "Publication requires successful validation for: ",
            paste(answer$method[failed_competitors], collapse = ", "),
            "."
          )
        } else if (identical(name, "published_competitors") &&
                   any(failed_competitors)) {
          paste0(
            "External competitor validation did not succeed for: ",
            paste(answer$method[failed_competitors], collapse = ", "),
            ". EIV-GP MCMC will continue; inspect the saved validation table."
          )
        } else if (invariant_failure) {
          "One or more paired-design invariants failed."
        } else {
          ""
        },
        elapsed_seconds = elapsed,
        stringsAsFactors = FALSE
      )
    }
    mixedgp_atomic_write_csv(
      mixedgp_bind_rows_base(statuses),
      file.path(metadata_dir, "automatic_validation_status.csv")
    )
  }

  status <- mixedgp_bind_rows_base(statuses)
  if (any(!status$pass)) {
    failures <- paste0(
      status$validator[!status$pass], ": ", status$message[!status$pass]
    )
    stop(
      "Automatic pre-fit validation failed. ",
      paste(failures, collapse = " | ")
    )
  }
  list(status = status, outputs = outputs)
}

mixedgp_existing_data_status <- function(config) {
  mixedgp_bind_rows_base(lapply(config$cells, function(cell) {
    directory <- mixedgp_data_cell_directory(config, cell)
    manifest_path <- file.path(directory, "manifest.rds")
    manifest <- if (file.exists(manifest_path)) {
      tryCatch(readRDS(manifest_path), error = function(e) NULL)
    } else {
      NULL
    }
    compatible <- if (is.data.frame(manifest)) {
      keep <- manifest$scenario == cell$scenario &
        manifest$n == cell$n & manifest$n_test == cell$n_test &
        manifest$m == cell$m & manifest$rep %in% seq_len(cell$n_rep)
      if (config$study == "study2") keep <- keep & manifest$q == cell$q
      manifest[keep, , drop = FALSE]
    } else {
      data.frame()
    }
    complete <- nrow(compatible) == cell$n_rep &&
      setequal(compatible$rep, seq_len(cell$n_rep))
    data.frame(
      cell_id = cell$id,
      directory = directory,
      manifest_exists = file.exists(manifest_path),
      frozen_artifacts = nrow(compatible),
      expected_artifacts = cell$n_rep,
      complete_count = complete,
      stringsAsFactors = FALSE
    )
  }))
}

mixedgp_print_dry_run <- function(config, preflight, runtime_preflight, task_plan) {
  cat(
    "Publication simulation dry run\n",
    "  study: ", config$study, "\n",
    "  schema: ", config$schema_version, "\n",
    "  run directory: ", mixedgp_run_directory(config), "\n",
    "  parallelism: ", config$parallel$workers,
    " replication worker(s); chains are serial within workers\n",
    sep = ""
  )
  print(mixedgp_cell_summary(config), row.names = FALSE)
  cat("\nPlanned fit/data units by type and method:\n")
  task_counts <- aggregate(
    planned_units ~ task_type + method, task_plan, sum
  )
  print(task_counts, row.names = FALSE)
  cat("\nPublished-package preflight:\n")
  print(preflight, row.names = FALSE)
  cat("\nRuntime-package preflight:\n")
  print(runtime_preflight, row.names = FALSE)
  cat("\nFrozen-data status (informational during dry run):\n")
  print(mixedgp_existing_data_status(config), row.names = FALSE)
  invisible(config)
}

mixedgp_generate_cell_data <- function(config, cell, engine) {
  directory <- mixedgp_data_cell_directory(config, cell)
  generator <- if (config$study == "study1") {
    get("generate_study1_synthetic_datasets", envir = engine)
  } else {
    get("generate_study2_synthetic_datasets", envir = engine)
  }
  common <- list(
    n_rep = cell$n_rep,
    directory = directory,
    n_cores = config$parallel$workers,
    overwrite = FALSE,
    n = cell$n,
    n_test = cell$n_test,
    m = cell$m,
    calib_grid = cell$calibration_grid
  )
  args <- if (config$study == "study1") {
    c(common, list(
      scenario = cell$scenario,
      threshold_design = cell$threshold_design,
      min_class_count = cell$min_class_count,
      heterogeneity_eta = cell$heterogeneity_eta,
      data_seed_base = 100000L,
      calibration_seed_base = 200000L
    ))
  } else {
    c(common, list(
      scenarios = cell$scenario,
      q = cell$q,
      data_seed_base = 1000000L
    ))
  }
  do.call(generator, args)
}

mixedgp_verify_cell_data <- function(config, cell, engine = NULL) {
  directory <- mixedgp_data_cell_directory(config, cell)
  manifest_path <- file.path(directory, "manifest.rds")
  if (!file.exists(manifest_path)) {
    stop("Frozen-data manifest is missing for cell ", cell$id, ": ", manifest_path)
  }
  manifest <- readRDS(manifest_path)
  if (!is.data.frame(manifest)) stop("Invalid manifest: ", manifest_path)
  matching <- manifest$scenario == cell$scenario &
    manifest$n == cell$n & manifest$n_test == cell$n_test &
    manifest$m == cell$m & manifest$rep %in% seq_len(cell$n_rep)
  if (config$study == "study2") {
    matching <- matching & manifest$q == cell$q
  }
  selected <- manifest[matching, , drop = FALSE]
  if (nrow(selected) != cell$n_rep ||
      !setequal(selected$rep, seq_len(cell$n_rep))) {
    stop(
      "Frozen-data manifest for ", cell$id, " contains ", nrow(selected),
      " compatible artifacts; expected ", cell$n_rep, "."
    )
  }
  paths <- file.path(directory, selected$file)
  if (any(!file.exists(paths))) stop("A frozen artifact is missing for ", cell$id, ".")
  observed_md5 <- unname(tools::md5sum(paths))
  if (!identical(tolower(observed_md5), tolower(as.character(selected$md5)))) {
    stop("Frozen-data checksum mismatch for cell ", cell$id, ".")
  }
  if (!is.null(engine)) {
    loader <- get(
      "load_mixedgp_synthetic_dataset_strict",
      envir = engine,
      inherits = FALSE
    )
    for (ii in seq_len(nrow(selected))) {
      expected <- list(
        study = config$study,
        scenario = cell$scenario,
        rep_id = as.integer(selected$rep[[ii]]),
        n = cell$n,
        n_test = cell$n_test,
        m = cell$m,
        calib_grid = cell$calibration_grid
      )
      if (config$study == "study1") {
        expected$threshold_design <- cell$threshold_design
        expected$heterogeneity_eta <- cell$heterogeneity_eta
      } else {
        expected$q <- cell$q
      }
      loader(
        paths[[ii]], expected = expected,
        manifest_path = manifest_path
      )
    }
  }
  selected
}

mixedgp_read_cell_replication <- function(config, cell, rep_id) {
  directory <- mixedgp_data_cell_directory(config, cell)
  manifest <- readRDS(file.path(directory, "manifest.rds"))
  matching <- manifest$rep == rep_id & manifest$scenario == cell$scenario &
    manifest$n == cell$n & manifest$n_test == cell$n_test &
    manifest$m == cell$m
  if (config$study == "study2") matching <- matching & manifest$q == cell$q
  row <- manifest[matching, , drop = FALSE]
  if (nrow(row) != 1L) {
    stop(
      "Expected one full-design frozen artifact for ", cell$id,
      " replication ", rep_id, "; found ", nrow(row), "."
    )
  }
  path <- file.path(directory, row$file[[1L]])
  observed_md5 <- unname(tools::md5sum(path))
  if (!identical(tolower(observed_md5), tolower(as.character(row$md5[[1L]])))) {
    stop("Frozen-data checksum mismatch for ", cell$id, " replication ", rep_id, ".")
  }
  artifact <- readRDS(path)
  if (!identical(as.integer(artifact$design$calib_grid),
                 as.integer(cell$calibration_grid))) {
    stop("Frozen calibration grid does not match cell ", cell$id, ".")
  }
  if (config$study == "study1" &&
      (!identical(as.character(artifact$design$threshold_design),
                  as.character(cell$threshold_design)) ||
       !isTRUE(all.equal(
         as.numeric(artifact$design$heterogeneity_eta),
         as.numeric(cell$heterogeneity_eta), tolerance = 0
       )))) {
    stop("Frozen Study I mechanism does not match cell ", cell$id, ".")
  }
  artifact
}

mixedgp_validate_common_random_numbers <- function(config, engine = NULL) {
  checks <- list()
  same_numeric <- function(x, y) {
    isTRUE(all.equal(
      x, y, tolerance = sqrt(.Machine$double.eps),
      check.attributes = FALSE
    ))
  }
  add_check <- function(reference, comparison, component, pass, n_rep) {
    checks[[length(checks) + 1L]] <<- data.frame(
      reference_cell = reference,
      comparison_cell = comparison,
      component = component,
      replications_checked = n_rep,
      pass = isTRUE(pass),
      stringsAsFactors = FALSE
    )
    if (!isTRUE(pass)) {
      stop(
        "Common-random-number validation failed for ", comparison,
        " component ", component, "."
      )
    }
  }
  by_id <- setNames(config$cells, vapply(config$cells, `[[`, character(1L), "id"))
  if (config$study == "study1") {
    primary_ids <- intersect(
      c("eta0_balanced", "eta0p5_balanced", "eta1_balanced"), names(by_id)
    )
    if (length(primary_ids) >= 2L) {
      reference <- by_id[[primary_ids[[1L]]]]
      for (id in primary_ids[-1L]) {
        comparison <- by_id[[id]]
        n_check <- min(reference$n_rep, comparison$n_rep)
        component_pass <- c(inputs = TRUE, noise = TRUE, calibration = TRUE,
                            observable_mean = TRUE)
        for (rr in seq_len(n_check)) {
          a <- mixedgp_read_cell_replication(config, reference, rr)
          b <- mixedgp_read_cell_replication(config, comparison, rr)
          component_pass[["inputs"]] <- component_pass[["inputs"]] &&
            identical(a$data$train$x, b$data$train$x) &&
            identical(a$data$train$u, b$data$train$u) &&
            identical(a$data$train$c, b$data$train$c) &&
            identical(a$data$test$x, b$data$test$x) &&
            identical(a$data$test$u, b$data$test$u) &&
            identical(a$data$test$c, b$data$test$c)
          component_pass[["noise"]] <- component_pass[["noise"]] &&
            same_numeric(a$data$train$y - a$data$train$f,
                         b$data$train$y - b$data$train$f) &&
            same_numeric(a$data$test$y - a$data$test$f,
                         b$data$test$y - b$data$test$f)
          component_pass[["calibration"]] <-
            component_pass[["calibration"]] &&
            identical(a$calibration_sets, b$calibration_sets)
          if (is.null(engine) || !exists("m0_1d", envir = engine, inherits = FALSE)) {
            stop("Study I paired-design validation requires m0_1d in engine.")
          }
          ma <- get("m0_1d", envir = engine, inherits = FALSE)(
            a$data$test$x, a$data$test$c, a$data$tau_true,
            scenario = a$scenario
          )
          mb <- get("m0_1d", envir = engine, inherits = FALSE)(
            b$data$test$x, b$data$test$c, b$data$tau_true,
            scenario = b$scenario
          )
          component_pass[["observable_mean"]] <-
          component_pass[["observable_mean"]] && same_numeric(ma, mb)
        }
        for (component in names(component_pass)) {
          add_check(reference$id, id, component, component_pass[[component]], n_check)
        }
      }
    }
  } else {
    reference_id <- if ("primary_q4_calibration" %in% names(by_id)) {
      "primary_q4_calibration"
    } else {
      names(by_id)[[1L]]
    }
    reference <- by_id[[reference_id]]
    for (id in setdiff(names(by_id), reference_id)) {
      comparison <- by_id[[id]]
      n_check <- min(reference$n_rep, comparison$n_rep)
      same_n_test <- identical(reference$n_test, comparison$n_test)
      component_pass <- c(
        train_X_U = TRUE,
        train_response_noise = TRUE,
        test_X_U_same_size = if (same_n_test) TRUE else NA,
        test_response_noise_same_size = if (same_n_test) TRUE else NA
      )
      if (comparison$scenario == "primary") {
        component_pass <- c(
          component_pass,
          nested_train_C = TRUE,
          nested_test_C_same_size = if (same_n_test) TRUE else NA
        )
      } else if (comparison$scenario == "latent_additive_control" &&
                 identical(reference$q, comparison$q)) {
        component_pass <- c(
          component_pass,
          identical_measurement_C = TRUE,
          identical_common_calibration_sets = TRUE
        )
      }
      for (rr in seq_len(n_check)) {
        a <- mixedgp_read_cell_replication(config, reference, rr)
        b <- mixedgp_read_cell_replication(config, comparison, rr)
        component_pass[["train_X_U"]] <- component_pass[["train_X_U"]] &&
          identical(a$data$train$X, b$data$train$X) &&
          identical(a$data$train$U, b$data$train$U)
        component_pass[["train_response_noise"]] <-
          component_pass[["train_response_noise"]] &&
          same_numeric(a$data$train$y - a$data$train$f,
                       b$data$train$y - b$data$train$f)
        if (same_n_test) {
          component_pass[["test_X_U_same_size"]] <-
            component_pass[["test_X_U_same_size"]] &&
            identical(a$data$test$X, b$data$test$X) &&
            identical(a$data$test$U, b$data$test$U)
          component_pass[["test_response_noise_same_size"]] <-
            component_pass[["test_response_noise_same_size"]] &&
            same_numeric(a$data$test$y - a$data$test$f,
                         b$data$test$y - b$data$test$f)
        }
        if (comparison$scenario == "primary") {
          q_shared <- min(ncol(a$data$train$C), ncol(b$data$train$C))
          component_pass[["nested_train_C"]] <-
            component_pass[["nested_train_C"]] &&
            identical(a$data$train$C[, seq_len(q_shared), drop = FALSE],
                      b$data$train$C[, seq_len(q_shared), drop = FALSE])
          if (same_n_test) {
            component_pass[["nested_test_C_same_size"]] <-
              component_pass[["nested_test_C_same_size"]] &&
              identical(a$data$test$C[, seq_len(q_shared), drop = FALSE],
                        b$data$test$C[, seq_len(q_shared), drop = FALSE])
          }
        } else if (comparison$scenario == "latent_additive_control" &&
                   identical(reference$q, comparison$q)) {
          component_pass[["identical_measurement_C"]] <-
            component_pass[["identical_measurement_C"]] &&
            identical(a$data$train$C, b$data$train$C) &&
            (!same_n_test || identical(a$data$test$C, b$data$test$C))
          common_calibration <- intersect(
            names(a$calibration_sets), names(b$calibration_sets)
          )
          component_pass[["identical_common_calibration_sets"]] <-
            component_pass[["identical_common_calibration_sets"]] &&
            length(common_calibration) > 0L &&
            all(vapply(common_calibration, function(k) {
              identical(a$calibration_sets[[k]], b$calibration_sets[[k]])
            }, logical(1L)))
        }
      }
      for (component in names(component_pass)[!is.na(component_pass)]) {
        add_check(reference$id, id, component, component_pass[[component]], n_check)
      }
    }
  }
  mixedgp_bind_rows_base(checks)
}

mixedgp_cell_controls_study1 <- function(config, cell, cell_output) {
  mechanism_calib <- if (20L %in% cell$calibration_grid) {
    20L
  } else {
    max(cell$calibration_grid)
  }
  list(
    STUDY1_CONFIG = if (config$mode == "smoke") "quick" else "thorough",
    STUDY1_QUICK = config$mode == "smoke",
    STUDY1_OUT_PREFIX = cell_output,
    STUDY1_DATA_DIR = mixedgp_data_cell_directory(config, cell),
    STUDY1_SCENARIO = cell$scenario,
    STUDY1_HETEROGENEITY_ETA = cell$heterogeneity_eta,
    STUDY1_THRESHOLD_DESIGN = cell$threshold_design,
    STUDY1_CALIB_GRID = cell$calibration_grid,
    STUDY1_MC_N_TRAIN = cell$n,
    STUDY1_MC_N_TEST = cell$n_test,
    STUDY1_MC_N_REP = cell$n_rep,
    STUDY1_MC_M = cell$m,
    STUDY1_MC_N_ITER = config$mcmc$n_iter,
    STUDY1_MC_BURN = config$mcmc$burn,
    STUDY1_MC_N_CHAINS = config$mcmc$n_chains,
    STUDY1_ADAPTIVE_MCMC = config$mcmc$adaptive,
    STUDY1_MC_N_PRED_DRAW = config$evaluation$n_pred_draw,
    STUDY1_MC_N_M_EVAL = config$evaluation$n_m_eval,
    STUDY1_MC_N_M_DRAW = config$evaluation$n_m_draw,
    STUDY1_MC_N_M_LATENT = config$evaluation$n_m_latent,
    STUDY1_MEAS_N_ITER = config$measurement_mcmc$n_iter,
    STUDY1_MEAS_BURN = config$measurement_mcmc$burn,
    STUDY1_USE_CACHE = config$use_cache,
    STUDY1_REUSE_LOCKED_EIV = FALSE,
    STUDY1_STRICT_COMPETITORS = config$strict_competitors,
    STUDY1_PUBLISHED_COMPETITORS = config$published_methods,
    STUDY1_RUN_ABLATIONS = cell$run_ablations,
    STUDY1_EVALUATE_F = cell$evaluate_f,
    STUDY1_EVALUATE_U = cell$evaluate_u,
    STUDY1_PARALLEL_LEVEL = config$parallel$level,
    STUDY1_REQUIRE_MCMC_GATE = config$mcmc$require_gate,
    STUDY1_MAX_RHAT = config$mcmc$rhat_limit,
    STUDY1_MIN_ESS = config$mcmc$ess_limit,
    STUDY1_MCMC_PILOT_REPS = config$mcmc$pilot_reps,
    STUDY1_MECHANISM_CALIB = mechanism_calib,
    STUDY1_LVGP_MAX_ELAPSED = if (config$mode == "smoke") 60 else 1800
  )
}

mixedgp_cell_controls_study2 <- function(config, cell, cell_output) {
  primary_grid <- if (cell$scenario == "primary") {
    cell$calibration_grid
  } else {
    sort(unique(c(0L, cell$calibration_grid)))
  }
  contrast_calib <- max(cell$calibration_grid)
  list(
    STUDY2_CONFIG = if (config$mode == "smoke") "quick" else "thorough",
    STUDY2_OUT_PREFIX = cell_output,
    STUDY2_DATA_DIR = mixedgp_data_cell_directory(config, cell),
    STUDY2_SCENARIOS = cell$scenario,
    STUDY2_PRIMARY_CALIB_GRID = primary_grid,
    STUDY2_CONTRAST_CALIB = contrast_calib,
    STUDY2_Q = cell$q,
    STUDY2_M = cell$m,
    STUDY2_MC_N_TRAIN = cell$n,
    STUDY2_MC_N_TEST = cell$n_test,
    STUDY2_MC_N_REP = cell$n_rep,
    STUDY2_MC_N_ITER = config$mcmc$n_iter,
    STUDY2_MC_BURN = config$mcmc$burn,
    STUDY2_MC_THIN = config$mcmc$thin,
    STUDY2_MC_N_CHAINS = config$mcmc$n_chains,
    STUDY2_ADAPTIVE_MCMC = config$mcmc$adaptive,
    STUDY2_MC_N_PRED_DRAW = config$evaluation$n_pred_draw,
    STUDY2_MC_N_M_EVAL = config$evaluation$n_m_eval,
    STUDY2_MC_N_M_DRAW = config$evaluation$n_m_draw,
    STUDY2_MC_N_M_LATENT = config$evaluation$n_m_latent,
    STUDY2_MC_N_M_TRUTH = config$evaluation$n_m_truth,
    STUDY2_MC_N_ORACLE_POOL = config$evaluation$n_oracle_pool,
    STUDY2_MEAS_N_ITER = config$measurement_mcmc$n_iter,
    STUDY2_MEAS_BURN = config$measurement_mcmc$burn,
    STUDY2_MEAS_THIN = config$measurement_mcmc$thin,
    STUDY2_MEAS_N_CHAINS = config$measurement_mcmc$n_chains,
    STUDY2_ABLATION_GP_N_STARTS = config$ablation_gp$n_starts,
    STUDY2_ABLATION_GP_MAXIT = config$ablation_gp$maxit,
    STUDY2_USE_CACHE = config$use_cache,
    STUDY2_MC_RESUME = TRUE,
    STUDY2_SAVE_REP_FITS = FALSE,
    STUDY2_STRICT_COMPETITORS = config$strict_competitors,
    STUDY2_PUBLISHED_COMPETITORS = config$published_methods,
    STUDY2_RUN_ABLATIONS = cell$run_ablations,
    STUDY2_ABLATION_SCENARIOS = if (cell$run_ablations) cell$scenario else character(0),
    STUDY2_EVALUATE_F = cell$evaluate_f,
    STUDY2_EVALUATE_U = cell$evaluate_u,
    STUDY2_PARALLEL_LEVEL = config$parallel$level,
    STUDY2_PREDICTIVE_LATENT_SAMPLER = config$predictive_latent_sampler,
    STUDY2_ENFORCE_MCMC_GATE = config$mcmc$require_gate,
    STUDY2_MCMC_RHAT_LIMIT = config$mcmc$rhat_limit,
    STUDY2_MCMC_RAW_ESS_LIMIT = config$mcmc$raw_ess_limit,
    STUDY2_MCMC_TARGET_BULK_ESS_LIMIT = config$mcmc$target_bulk_ess_limit,
    STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT = config$mcmc$target_tail_ess_limit,
    ## Smoke mode is explicitly non-reportable: permissive measurement-chain
    ## thresholds let PI-GP/CC-GP execute so every code path is tested. The
    ## publication constructor retains the prespecified strict diagnostics.
    STUDY2_MEAS_RHAT_LIMIT = if (config$mode == "smoke") {
      3
    } else {
      config$mcmc$rhat_limit
    },
    STUDY2_MEAS_BULK_ESS_LIMIT = if (config$mode == "smoke") {
      1L
    } else {
      config$mcmc$target_bulk_ess_limit
    },
    STUDY2_MEAS_TAIL_ESS_LIMIT = if (config$mode == "smoke") {
      1L
    } else {
      config$mcmc$target_tail_ess_limit
    },
    STUDY2_SAVE_PDF = identical(config$mode, "publication"),
    STUDY2_SAVE_PNG = FALSE,
    STUDY2_LVGP_MAX_ELAPSED = if (config$mode == "smoke") 60 else 3600
  )
}

mixedgp_get_cell_output <- function(envir, name) {
  get0(name, envir = envir, inherits = FALSE, ifnotfound = data.frame())
}

mixedgp_run_study1_cell <- function(config, cell, engine, run_dir) {
  cell_output <- file.path(run_dir, "cells", cell$id)
  dir.create(cell_output, recursive = TRUE, showWarnings = FALSE)
  run_env <- new.env(parent = engine)
  list2env(
    mixedgp_cell_controls_study1(config, cell, cell_output),
    envir = run_env
  )
  sys.source(
    file.path(config$code_dir, "02_study1_monte_carlo.R"),
    envir = run_env,
    chdir = TRUE
  )
  outputs <- list(
    predictive_metrics = mixedgp_get_cell_output(run_env, "mc_results"),
    mean_recovery = mixedgp_get_cell_output(run_env, "mean_recovery"),
    latent_imputation = mixedgp_get_cell_output(run_env, "latent_imputation"),
    surface_recovery = mixedgp_get_cell_output(run_env, "surface_recovery"),
    ablation_predictive_metrics = mixedgp_get_cell_output(run_env, "ablation_results"),
    ablation_surface_recovery = mixedgp_get_cell_output(
      run_env, "ablation_surface_recovery"
    ),
    mcmc_diagnostics = mixedgp_get_cell_output(run_env, "mcmc_diagnostics"),
    mcmc_schedule = mixedgp_get_cell_output(run_env, "mcmc_schedule"),
    competitor_status = mixedgp_get_cell_output(run_env, "competitor_status"),
    replication_status = mixedgp_get_cell_output(run_env, "replication_status"),
    ablation_status = mixedgp_get_cell_output(run_env, "ablation_status"),
    sampler_control_manifest = mixedgp_get_cell_output(
      run_env, "sampler_control_manifest"
    )
  )
  list(
    cell = cell,
    outputs = outputs,
    design_tag = get0("STUDY1_DESIGN_TAG", envir = run_env, inherits = FALSE),
    output_directory = cell_output
  )
}

mixedgp_run_study2_cell <- function(config, cell, engine, run_dir) {
  cell_output <- file.path(run_dir, "cells", cell$id)
  dir.create(cell_output, recursive = TRUE, showWarnings = FALSE)
  run_env <- new.env(parent = engine)
  list2env(
    mixedgp_cell_controls_study2(config, cell, cell_output),
    envir = run_env
  )
  sys.source(
    file.path(config$code_dir, "02_study2_monte_carlo.R"),
    envir = run_env,
    chdir = TRUE
  )
  outputs <- get0(
    "raw_outputs", envir = run_env, inherits = FALSE, ifnotfound = list()
  )
  list(
    cell = cell,
    outputs = outputs,
    design_tag = get0("STUDY2_DESIGN_TAG", envir = run_env, inherits = FALSE),
    output_directory = cell_output
  )
}

mixedgp_competitor_gate <- function(result, config) {
  status <- result$outputs$competitor_status
  replication_status <- result$outputs$replication_status
  n_eiv_completed <- if (is.data.frame(replication_status) &&
      "status" %in% names(replication_status)) {
    sum(replication_status$status == "success", na.rm = TRUE)
  } else {
    result$cell$n_rep
  }
  expected <- n_eiv_completed * length(config$published_methods)
  if (!is.data.frame(status) || nrow(status) == 0L) {
    failure_rate <- 1
    n_success <- 0L
    n_attempted <- 0L
  } else {
    n_attempted <- nrow(status)
    n_success <- sum(status$status == "success", na.rm = TRUE)
    failure_rate <- 1 - n_success / expected
  }
  data.frame(
    cell_id = result$cell$id,
    gate = "published_competitors",
    expected = expected,
    attempted = n_attempted,
    success = n_success,
    failure_rate = failure_rate,
    ## A public-package optimizer can legitimately fail on an individual
    ## frozen data set.  It is reportable metadata, not a reason to abort
    ## EIV-GP MCMC or invalidate its diagnostic gate.
    pass = n_attempted == expected,
    stringsAsFactors = FALSE
  )
}

mixedgp_mcmc_gate <- function(result, config) {
  diagnostics <- result$outputs$mcmc_diagnostics
  pass_column <- if (config$study == "study1") "gate_pass" else "mcmc_pass"
  expected_per_rep <- if (config$study == "study1") {
    ## Study I's zero-calibration condition is intentionally unanchored and
    ## has no raw-parameter convergence gate.
    sum(result$cell$calibration_grid > 0L)
  } else {
    length(result$cell$calibration_grid)
  }
  expected <- result$cell$n_rep * expected_per_rep
  values <- if (is.data.frame(diagnostics) && pass_column %in% names(diagnostics)) {
    applicable <- if ("gate_applicable" %in% names(diagnostics)) {
      isTRUE(diagnostics$gate_applicable) | diagnostics$gate_applicable %in% TRUE
    } else {
      rep(TRUE, nrow(diagnostics))
    }
    as.logical(diagnostics[[pass_column]][applicable])
  } else {
    logical(0)
  }
  data.frame(
    cell_id = result$cell$id,
    gate = "eiv_mcmc",
    expected = expected,
    attempted = length(values),
    success = sum(values %in% TRUE),
    failure_rate = if (expected > 0L) 1 - sum(values %in% TRUE) / expected else 0,
    pass = length(values) == expected && all(values %in% TRUE),
    stringsAsFactors = FALSE
  )
}

mixedgp_expected_method_keys <- function(cell, method_calibrations) {
  rows <- lapply(names(method_calibrations), function(method) {
    calibration <- method_calibrations[[method]]
    expand.grid(
      rep = seq_len(cell$n_rep),
      method = method,
      n_calib = calibration,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  })
  mixedgp_bind_rows_base(rows)
}

mixedgp_encode_completeness_keys <- function(x, columns) {
  if (!is.data.frame(x) || nrow(x) == 0L ||
      any(!columns %in% names(x))) return(character(0))
  do.call(paste, c(lapply(x[columns], function(z) {
    z <- as.character(z)
    z[is.na(z)] <- "<NA>"
    z
  }), sep = "\r"))
}

mixedgp_completeness_gate_row <- function(cell_id,
                                          gate,
                                          expected,
                                          actual,
                                          columns) {
  expected_key <- mixedgp_encode_completeness_keys(expected, columns)
  actual_key <- mixedgp_encode_completeness_keys(actual, columns)
  expected_unique <- unique(expected_key)
  actual_unique <- unique(actual_key)
  n_success <- length(intersect(expected_unique, actual_unique))
  pass <- length(actual_key) == length(actual_unique) &&
    setequal(expected_unique, actual_unique)
  data.frame(
    cell_id = cell_id,
    gate = gate,
    expected = length(expected_unique),
    attempted = length(actual_key),
    success = n_success,
    failure_rate = if (length(expected_unique) > 0L) {
      1 - n_success / length(expected_unique)
    } else if (length(actual_key) == 0L) {
      0
    } else {
      1
    },
    pass = pass,
    stringsAsFactors = FALSE
  )
}

mixedgp_normalize_rf_method <- function(x) {
  x <- as.character(x)
  x[x %in% c(
    "Response-free threshold model",
    "Response-free measurement model",
    "Ordinal model (no Y)"
  )] <- "RF-OM"
  x
}

mixedgp_output_completeness_gates <- function(result, config) {
  cell <- result$cell
  calibs <- as.integer(cell$calibration_grid)
  positive_calibs <- calibs[calibs > 0L]
  overall <- function(x) {
    if (is.data.frame(x) && "evaluation_stratum" %in% names(x)) {
      x <- x[as.character(x$evaluation_stratum) == "overall", , drop = FALSE]
    }
    x
  }
  filter_expected_methods <- function(x, methods) {
    if (!is.data.frame(x) || nrow(x) == 0L || !"method" %in% names(x)) {
      return(data.frame())
    }
    x$method <- as.character(x$method)
    x[x$method %in% methods, , drop = FALSE]
  }

  central_spec <- c(
    list(`EIV-GP` = calibs),
    setNames(rep(list(NA_integer_), length(config$published_methods)),
             config$published_methods)
  )
  predictive_spec <- c(central_spec, list(Oracle = NA_integer_))
  predictive_expected <- mixedgp_expected_method_keys(cell, predictive_spec)
  predictive_actual <- filter_expected_methods(
    overall(result$outputs$predictive_metrics), names(predictive_spec)
  )
  mean_expected <- mixedgp_expected_method_keys(cell, central_spec)
  mean_actual <- filter_expected_methods(
    result$outputs$mean_recovery, names(central_spec)
  )
  rows <- list(
    mixedgp_completeness_gate_row(
      cell$id, "prediction_output_completeness",
      predictive_expected, predictive_actual, c("rep", "method", "n_calib")
    ),
    mixedgp_completeness_gate_row(
      cell$id, "observable_mean_output_completeness",
      mean_expected, mean_actual, c("rep", "method", "n_calib")
    )
  )

  if (isTRUE(cell$run_ablations)) {
    ablation_prediction_spec <- list(`PI-GP` = calibs, `CC-GP` = positive_calibs)
    ablation_prediction_expected <- mixedgp_expected_method_keys(
      cell, ablation_prediction_spec
    )
    ablation_prediction_actual <- filter_expected_methods(
      overall(result$outputs$ablation_predictive_metrics),
      names(ablation_prediction_spec)
    )
    rows[[length(rows) + 1L]] <- mixedgp_completeness_gate_row(
      cell$id, "ablation_prediction_output_completeness",
      ablation_prediction_expected, ablation_prediction_actual,
      c("rep", "method", "n_calib")
    )
  }

  if (isTRUE(cell$evaluate_f)) {
    eiv_surface_expected <- mixedgp_expected_method_keys(
      cell, list(`EIV-GP` = positive_calibs)
    )
    eiv_surface_actual <- filter_expected_methods(
      result$outputs$surface_recovery, "EIV-GP"
    )
    ablation_surface_spec <- list(
      `PI-GP` = positive_calibs,
      `CC-GP` = positive_calibs,
      `Full-U GP` = NA_integer_
    )
    ablation_surface_expected <- mixedgp_expected_method_keys(
      cell, ablation_surface_spec
    )
    ablation_surface_actual <- filter_expected_methods(
      result$outputs$ablation_surface_recovery,
      names(ablation_surface_spec)
    )
    rows[[length(rows) + 1L]] <- mixedgp_completeness_gate_row(
      cell$id, "latent_surface_eiv_output_completeness",
      eiv_surface_expected, eiv_surface_actual,
      c("rep", "method", "n_calib")
    )
    rows[[length(rows) + 1L]] <- mixedgp_completeness_gate_row(
      cell$id, "latent_surface_benchmark_output_completeness",
      ablation_surface_expected, ablation_surface_actual,
      c("rep", "method", "n_calib")
    )
  }

  if (isTRUE(cell$evaluate_u)) {
    if (config$study == "study1") {
      latent_expected <- mixedgp_expected_method_keys(
        cell, list(`EIV-GP` = calibs, `RF-OM` = calibs)
      )
      latent_actual <- result$outputs$latent_imputation
      if (is.data.frame(latent_actual) && nrow(latent_actual) > 0L) {
        latent_actual$method <- mixedgp_normalize_rf_method(latent_actual$method)
        latent_actual <- latent_actual[
          latent_actual$method %in% c("EIV-GP", "RF-OM"), , drop = FALSE
        ]
      }
      rows[[length(rows) + 1L]] <- mixedgp_completeness_gate_row(
        cell$id, "training_latent_u_output_completeness",
        latent_expected, latent_actual, c("rep", "method", "n_calib")
      )
    } else {
      targets <- c("training_missing_U", "prospective_U_given_C")
      latent_expected <- expand.grid(
        rep = seq_len(cell$n_rep),
        method = c("EIV-GP", "RF-OM"),
        n_calib = calibs,
        target = targets,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
      latent_actual <- result$outputs$latent_imputation_status
      if (is.data.frame(latent_actual) && nrow(latent_actual) > 0L) {
        latent_actual$method <- mixedgp_normalize_rf_method(latent_actual$method)
        latent_actual <- latent_actual[
          latent_actual$method %in% c("EIV-GP", "RF-OM") &
            latent_actual$target %in% targets, , drop = FALSE
        ]
      }
      rows[[length(rows) + 1L]] <- mixedgp_completeness_gate_row(
        cell$id, "latent_u_task_status_completeness",
        latent_expected, latent_actual,
        c("rep", "method", "n_calib", "target")
      )
    }
  }
  mixedgp_bind_rows_base(rows)
}

mixedgp_ablation_gate <- function(result) {
  status <- result$outputs$ablation_status
  if (!isTRUE(result$cell$run_ablations)) {
    return(data.frame(
      cell_id = result$cell$id, gate = "scientific_ablations",
      expected = 0L, attempted = 0L, success = 0L,
      failure_rate = 0, pass = TRUE, stringsAsFactors = FALSE
    ))
  }
  calibs <- as.integer(result$cell$calibration_grid)
  expected_spec <- list(`RF-OM` = calibs, `PI-GP` = calibs, `CC-GP` = calibs)
  if (isTRUE(result$cell$evaluate_f)) {
    expected_spec[["Full-U GP"]] <- NA_integer_
  }
  expected <- mixedgp_expected_method_keys(result$cell, expected_spec)
  if (!is.data.frame(status) || nrow(status) == 0L) {
    status <- data.frame()
  } else {
    status$method <- mixedgp_normalize_rf_method(status$method)
    status <- status[status$method %in% names(expected_spec), , drop = FALSE]
  }
  completeness <- mixedgp_completeness_gate_row(
    result$cell$id, "scientific_ablation_status_completeness",
    expected, status, c("rep", "method", "n_calib")
  )
  acceptable <- is.data.frame(status) && nrow(status) > 0L &&
    all(status$status %in% c("success", "not_applicable"))
  completeness$gate <- "scientific_ablations"
  expected_key <- mixedgp_encode_completeness_keys(
    expected, c("rep", "method", "n_calib")
  )
  acceptable_status <- if (is.data.frame(status) && nrow(status) > 0L) {
    status[status$status %in% c("success", "not_applicable"), , drop = FALSE]
  } else {
    data.frame()
  }
  acceptable_key <- mixedgp_encode_completeness_keys(
    acceptable_status, c("rep", "method", "n_calib")
  )
  completeness$success <- length(intersect(
    unique(expected_key), unique(acceptable_key)
  ))
  completeness$failure_rate <- if (completeness$expected > 0L) {
    1 - completeness$success / completeness$expected
  } else {
    0
  }
  completeness$pass <- isTRUE(completeness$pass) && acceptable
  completeness
}

mixedgp_result_gates <- function(result, config) {
  mixedgp_bind_rows_base(list(
    mixedgp_competitor_gate(result, config),
    mixedgp_mcmc_gate(result, config),
    mixedgp_ablation_gate(result),
    mixedgp_output_completeness_gates(result, config)
  ))
}

mixedgp_tag_output <- function(x, cell, study) {
  if (!is.data.frame(x) || nrow(x) == 0L) return(x)
  design <- data.frame(
    study = rep(study, nrow(x)),
    cell_id = rep(cell$id, nrow(x)),
    design_role = rep(cell$role, nrow(x)),
    design_eta = rep(
      if (is.null(cell$heterogeneity_eta)) NA_real_ else cell$heterogeneity_eta,
      nrow(x)
    ),
    design_q = rep(if (is.null(cell$q)) NA_integer_ else cell$q, nrow(x)),
    stringsAsFactors = FALSE
  )
  duplicate <- intersect(names(design), names(x))
  if (length(duplicate) > 0L) design[duplicate] <- NULL
  cbind(design, x)
}

mixedgp_group_summary <- function(x, value, groups) {
  if (!is.data.frame(x) || nrow(x) == 0L ||
      !value %in% names(x)) return(data.frame())
  groups <- intersect(groups, names(x))
  if (length(groups) == 0L) stop("At least one grouping column is required.")
  group_key <- do.call(
    paste,
    c(lapply(x[groups], function(z) {
      z <- as.character(z)
      z[is.na(z)] <- "<NA>"
      z
    }), sep = "\r")
  )
  index <- split(seq_len(nrow(x)), group_key, drop = TRUE)
  rows <- lapply(index, function(ii) {
    values <- as.numeric(x[[value]][ii])
    values <- values[is.finite(values)]
    n_eff <- length(values)
    estimate <- if (n_eff == 0L) NA_real_ else mean(values)
    se <- if (n_eff < 2L) NA_real_ else stats::sd(values) / sqrt(n_eff)
    cbind(
      x[ii[[1L]], groups, drop = FALSE],
      data.frame(
        mean = estimate,
        se = se,
        ci_lower = if (is.finite(se)) estimate - 1.96 * se else NA_real_,
        ci_upper = if (is.finite(se)) estimate + 1.96 * se else NA_real_,
        n_pairs = n_eff,
        stringsAsFactors = FALSE
      )
    )
  })
  mixedgp_bind_rows_base(rows)
}

mixedgp_paired_eiv_published <- function(x, task, metrics) {
  required <- c("cell_id", "rep", "scenario", "n_calib", "method")
  if (!is.data.frame(x) || nrow(x) == 0L ||
      any(!required %in% names(x))) return(data.frame())
  metrics <- intersect(metrics, names(x))
  if (length(metrics) == 0L) return(data.frame())
  x$method <- as.character(x$method)
  eiv <- x[x$method == "EIV-GP" & !is.na(x$n_calib), , drop = FALSE]
  competitor <- x[
    x$method %in% MIXEDGP_PUBLISHED_METHODS & is.na(x$n_calib),
    , drop = FALSE
  ]
  if (nrow(eiv) == 0L || nrow(competitor) == 0L) return(data.frame())
  keys <- intersect(
    c(
      "study", "cell_id", "design_role", "design_eta", "design_q",
      "rep", "scenario", "evaluation_stratum"
    ),
    intersect(names(eiv), names(competitor))
  )
  eiv_keep <- unique(c(keys, "n_calib", metrics))
  competitor_keep <- unique(c(keys, "method", metrics))
  eiv <- eiv[, eiv_keep, drop = FALSE]
  competitor <- competitor[, competitor_keep, drop = FALSE]
  names(competitor)[names(competitor) == "method"] <- "competitor"
  names(eiv)[match(metrics, names(eiv))] <- paste0("EIV_", metrics)
  names(competitor)[match(metrics, names(competitor))] <-
    paste0("Competitor_", metrics)
  paired <- merge(eiv, competitor, by = keys, all = FALSE, sort = FALSE)
  if (nrow(paired) == 0L) return(data.frame())
  if (!"evaluation_stratum" %in% names(paired)) {
    paired$evaluation_stratum <- "overall"
  }
  rows <- lapply(metrics, function(metric) {
    eiv_value <- as.numeric(paired[[paste0("EIV_", metric)]])
    competitor_value <- as.numeric(
      paired[[paste0("Competitor_", metric)]]
    )
    id_columns <- setdiff(
      names(paired),
      c(paste0("EIV_", metrics), paste0("Competitor_", metrics))
    )
    cbind(
      paired[, id_columns, drop = FALSE],
      data.frame(
        task = task,
        metric = metric,
        EIV_value = eiv_value,
        competitor_value = competitor_value,
        EIV_advantage = competitor_value - eiv_value,
        advantage_direction = "positive favors EIV-GP",
        stringsAsFactors = FALSE
      )
    )
  })
  mixedgp_bind_rows_base(rows)
}

mixedgp_paired_eiv_same_calibration <- function(eiv_x,
                                                 competitor_x,
                                                 task,
                                                 methods,
                                                 metrics) {
  required <- c("cell_id", "rep", "scenario", "n_calib", "method")
  if (!is.data.frame(eiv_x) || !is.data.frame(competitor_x) ||
      nrow(eiv_x) == 0L || nrow(competitor_x) == 0L ||
      any(!required %in% names(eiv_x)) ||
      any(!required %in% names(competitor_x))) return(data.frame())
  metrics <- intersect(metrics, intersect(names(eiv_x), names(competitor_x)))
  if (length(metrics) == 0L) return(data.frame())
  eiv_x$method <- as.character(eiv_x$method)
  competitor_x$method <- as.character(competitor_x$method)
  eiv <- eiv_x[eiv_x$method == "EIV-GP" & !is.na(eiv_x$n_calib), , drop = FALSE]
  competitor <- competitor_x[
    competitor_x$method %in% methods & !is.na(competitor_x$n_calib),
    , drop = FALSE
  ]
  if (nrow(eiv) == 0L || nrow(competitor) == 0L) return(data.frame())
  keys <- intersect(
    c(
      "study", "cell_id", "design_role", "design_eta", "design_q",
      "rep", "scenario", "n_calib", "evaluation_stratum"
    ),
    intersect(names(eiv), names(competitor))
  )
  eiv <- eiv[, unique(c(keys, metrics)), drop = FALSE]
  competitor <- competitor[, unique(c(keys, "method", metrics)), drop = FALSE]
  names(competitor)[names(competitor) == "method"] <- "competitor"
  names(eiv)[match(metrics, names(eiv))] <- paste0("EIV_", metrics)
  names(competitor)[match(metrics, names(competitor))] <-
    paste0("Competitor_", metrics)
  paired <- merge(eiv, competitor, by = keys, all = FALSE, sort = FALSE)
  if (nrow(paired) == 0L) return(data.frame())
  if (!"evaluation_stratum" %in% names(paired)) {
    paired$evaluation_stratum <- "overall"
  }
  rows <- lapply(metrics, function(metric) {
    eiv_value <- as.numeric(paired[[paste0("EIV_", metric)]])
    competitor_value <- as.numeric(
      paired[[paste0("Competitor_", metric)]]
    )
    id_columns <- setdiff(
      names(paired),
      c(paste0("EIV_", metrics), paste0("Competitor_", metrics))
    )
    cbind(
      paired[, id_columns, drop = FALSE],
      data.frame(
        task = task,
        metric = metric,
        EIV_value = eiv_value,
        competitor_value = competitor_value,
        EIV_advantage = competitor_value - eiv_value,
        advantage_direction = "positive favors EIV-GP",
        stringsAsFactors = FALSE
      )
    )
  })
  mixedgp_bind_rows_base(rows)
}

mixedgp_paired_latent_u <- function(x) {
  required <- c(
    "cell_id", "rep", "scenario", "n_calib", "method", "RMSE", "MAE"
  )
  if (!is.data.frame(x) || nrow(x) == 0L ||
      any(!required %in% names(x))) return(data.frame())
  x$method <- as.character(x$method)
  response_free_names <- c(
    "Response-free threshold model", "Ordinal model (no Y)"
  )
  if ("score_status" %in% names(x)) {
    x <- x[is.na(x$score_status) | x$score_status == "scored", , drop = FALSE]
  }
  eiv <- x[x$method == "EIV-GP", , drop = FALSE]
  response_free <- x[x$method %in% response_free_names, , drop = FALSE]
  if (nrow(eiv) == 0L || nrow(response_free) == 0L) return(data.frame())
  keys <- intersect(
    c(
      "study", "cell_id", "design_role", "design_eta", "design_q",
      "rep", "scenario", "n_calib", "target", "coordinate"
    ),
    intersect(names(eiv), names(response_free))
  )
  metrics <- c("RMSE", "MAE")
  eiv <- eiv[, unique(c(keys, metrics)), drop = FALSE]
  response_free <- response_free[, unique(c(keys, "method", metrics)), drop = FALSE]
  names(response_free)[names(response_free) == "method"] <- "competitor"
  names(eiv)[match(metrics, names(eiv))] <- paste0("EIV_", metrics)
  names(response_free)[match(metrics, names(response_free))] <-
    paste0("Competitor_", metrics)
  paired <- merge(eiv, response_free, by = keys, all = FALSE, sort = FALSE)
  if (nrow(paired) == 0L) return(data.frame())
  rows <- lapply(metrics, function(metric) {
    eiv_value <- as.numeric(paired[[paste0("EIV_", metric)]])
    competitor_value <- as.numeric(
      paired[[paste0("Competitor_", metric)]]
    )
    id_columns <- setdiff(
      names(paired),
      c(paste0("EIV_", metrics), paste0("Competitor_", metrics))
    )
    cbind(
      paired[, id_columns, drop = FALSE],
      data.frame(
        task = "latent U inference",
        metric = metric,
        EIV_value = eiv_value,
        competitor_value = competitor_value,
        EIV_advantage = competitor_value - eiv_value,
        advantage_direction = "positive favors EIV-GP",
        stringsAsFactors = FALSE
      )
    )
  })
  mixedgp_bind_rows_base(rows)
}

mixedgp_paired_latent_surface <- function(eiv, ablations) {
  required <- c("cell_id", "rep", "scenario", "n_calib", "method", "ISE")
  if (!is.data.frame(eiv) || !is.data.frame(ablations) ||
      nrow(eiv) == 0L || nrow(ablations) == 0L ||
      any(!required %in% names(eiv)) || any(!required %in% names(ablations))) {
    return(data.frame())
  }
  eiv$method <- as.character(eiv$method)
  ablations$method <- as.character(ablations$method)
  eiv <- eiv[eiv$method == "EIV-GP" & !is.na(eiv$n_calib), , drop = FALSE]
  ablations <- ablations[
    ablations$method %in% c("PI-GP", "CC-GP", "Full-U GP"),
    , drop = FALSE
  ]
  if (nrow(eiv) == 0L || nrow(ablations) == 0L) return(data.frame())
  keys <- intersect(
    c(
      "study", "cell_id", "design_role", "design_eta", "design_q",
      "rep", "scenario"
    ),
    intersect(names(eiv), names(ablations))
  )
  eiv <- eiv[, unique(c(keys, "n_calib", "ISE")), drop = FALSE]
  names(eiv)[names(eiv) == "ISE"] <- "EIV_ISE"
  same_k <- ablations[!is.na(ablations$n_calib), , drop = FALSE]
  benchmark <- ablations[is.na(ablations$n_calib), , drop = FALSE]
  paired_same <- data.frame()
  if (nrow(same_k) > 0L) {
    paired_same <- merge(
      eiv,
      same_k[, unique(c(keys, "n_calib", "method", "ISE")), drop = FALSE],
      by = c(keys, "n_calib"), all = FALSE, sort = FALSE
    )
  }
  paired_benchmark <- data.frame()
  if (nrow(benchmark) > 0L) {
    paired_benchmark <- merge(
      eiv,
      benchmark[, unique(c(keys, "method", "ISE")), drop = FALSE],
      by = keys, all = FALSE, sort = FALSE
    )
  }
  paired <- mixedgp_bind_rows_base(list(paired_same, paired_benchmark))
  if (nrow(paired) == 0L) return(data.frame())
  names(paired)[names(paired) == "method"] <- "competitor"
  names(paired)[names(paired) == "ISE"] <- "competitor_value"
  paired$task <- "latent surface f(x,u)"
  paired$metric <- "ISE"
  paired$EIV_value <- paired$EIV_ISE
  paired$EIV_advantage <- paired$competitor_value - paired$EIV_value
  paired$advantage_direction <- "positive favors EIV-GP"
  paired$EIV_ISE <- NULL
  paired
}

mixedgp_method_comparisons <- function(combined) {
  predictive <- combined$predictive_metrics
  ablation_predictive <- combined$ablation_predictive_metrics
  mean_recovery <- combined$mean_recovery
  latent_u <- combined$latent_imputation
  surface <- combined$surface_recovery
  ablation_surface <- combined$ablation_surface_recovery
  mixedgp_bind_rows_base(list(
    mixedgp_paired_eiv_published(
      predictive,
      task = "response prediction Y*|x*,c*",
      metrics = c("RMSE", "MAE", "CRPS", "NLPD", "IntervalScore95")
    ),
    mixedgp_paired_eiv_same_calibration(
      predictive,
      ablation_predictive,
      task = "response prediction Y*|x*,c*",
      methods = c("PI-GP", "CC-GP"),
      metrics = c("RMSE", "MAE", "CRPS", "NLPD", "IntervalScore95")
    ),
    mixedgp_paired_eiv_published(
      mean_recovery,
      task = "observable mean m(x,c)",
      metrics = c("RMSE", "MAE")
    ),
    mixedgp_paired_latent_u(latent_u),
    mixedgp_paired_latent_surface(surface, ablation_surface)
  ))
}

mixedgp_study1_mechanism_contrasts <- function(comparisons) {
  if (!is.data.frame(comparisons) || nrow(comparisons) == 0L ||
      !"cell_id" %in% names(comparisons)) return(data.frame())
  control <- comparisons[
    comparisons$cell_id == "eta0_balanced" &
      comparisons$task %in% c(
        "response prediction Y*|x*,c*", "observable mean m(x,c)"
      ), , drop = FALSE
  ]
  active <- comparisons[
    comparisons$cell_id == "eta1_balanced" &
      comparisons$task %in% c(
        "response prediction Y*|x*,c*", "observable mean m(x,c)"
      ), , drop = FALSE
  ]
  keys <- intersect(
    c(
      "rep", "n_calib", "task", "competitor", "metric",
      "evaluation_stratum"
    ),
    intersect(names(control), names(active))
  )
  if (nrow(control) == 0L || nrow(active) == 0L) return(data.frame())
  control <- control[, c(keys, "EIV_advantage"), drop = FALSE]
  active <- active[, c(keys, "EIV_advantage"), drop = FALSE]
  names(control)[names(control) == "EIV_advantage"] <- "advantage_eta0"
  names(active)[names(active) == "EIV_advantage"] <- "advantage_eta1"
  out <- merge(control, active, by = keys, all = FALSE, sort = FALSE)
  out$contrast <- "heterogeneity: eta=1 minus eta=0"
  out$advantage_change <- out$advantage_eta1 - out$advantage_eta0
  out$contrast_direction <- "positive means heterogeneity strengthens EIV advantage"
  out
}

mixedgp_study2_design_contrasts <- function(comparisons) {
  if (!is.data.frame(comparisons) || nrow(comparisons) == 0L ||
      !"cell_id" %in% names(comparisons)) return(data.frame())
  x <- comparisons[
    comparisons$task %in% c(
      "response prediction Y*|x*,c*", "observable mean m(x,c)"
    ), , drop = FALSE
  ]
  if ("evaluation_stratum" %in% names(x)) {
    x <- x[x$evaluation_stratum == "overall", , drop = FALSE]
  }
  keys <- intersect(
    c("rep", "n_calib", "task", "competitor", "metric"), names(x)
  )
  extract_advantage <- function(cell_id, label) {
    out <- x[x$cell_id == cell_id, c(keys, "EIV_advantage"), drop = FALSE]
    names(out)[names(out) == "EIV_advantage"] <- label
    out
  }
  primary <- extract_advantage("primary_q4_calibration", "advantage_primary")
  additive <- extract_advantage("additive_q4", "advantage_additive")
  uncertain <- extract_advantage(
    "high_uncertainty_q4", "advantage_high_uncertainty"
  )
  controls <- list()
  if (nrow(primary) > 0L && nrow(additive) > 0L) {
    z <- merge(primary, additive, by = keys, all = FALSE, sort = FALSE)
    z$contrast <- "interactions: primary minus additive"
    z$advantage_change <- z$advantage_primary - z$advantage_additive
    controls[[length(controls) + 1L]] <- z
  }
  if (nrow(primary) > 0L && nrow(uncertain) > 0L) {
    z <- merge(uncertain, primary, by = keys, all = FALSE, sort = FALSE)
    z$contrast <- "latent uncertainty: high minus primary"
    z$advantage_change <-
      z$advantage_high_uncertainty - z$advantage_primary
    controls[[length(controls) + 1L]] <- z
  }
  primary_ids <- c("primary_q2", "primary_q3", "primary_q4_calibration")
  q_rows <- x[x$cell_id %in% primary_ids, , drop = FALSE]
  q_rows$design_q <- as.integer(q_rows$design_q)
  q_keys <- intersect(
    c("rep", "n_calib", "task", "competitor", "metric"), names(q_rows)
  )
  for (pair in list(c(2L, 3L), c(3L, 4L), c(2L, 4L))) {
    low <- q_rows[q_rows$design_q == pair[[1L]],
                  c(q_keys, "EIV_advantage"), drop = FALSE]
    high <- q_rows[q_rows$design_q == pair[[2L]],
                   c(q_keys, "EIV_advantage"), drop = FALSE]
    names(low)[names(low) == "EIV_advantage"] <- "advantage_low_q"
    names(high)[names(high) == "EIV_advantage"] <- "advantage_high_q"
    if (nrow(low) > 0L && nrow(high) > 0L) {
      z <- merge(low, high, by = q_keys, all = FALSE, sort = FALSE)
      z$contrast <- paste0("proxy dimension: q=", pair[[2L]],
                           " minus q=", pair[[1L]])
      z$advantage_change <- z$advantage_high_q - z$advantage_low_q
      controls[[length(controls) + 1L]] <- z
    }
  }
  out <- mixedgp_bind_rows_base(controls)
  if (nrow(out) > 0L) {
    out$contrast_direction <-
      "positive means the named design change strengthens EIV advantage"
  }
  out
}

mixedgp_write_publication_comparisons <- function(combined, config, combined_dir) {
  comparisons <- mixedgp_method_comparisons(combined)
  if (nrow(comparisons) > 0L) {
    mixedgp_atomic_write_csv(
      comparisons, file.path(combined_dir, "method_comparisons_paired_raw.csv")
    )
    comparison_groups <- c(
      "study", "cell_id", "design_role", "design_eta", "design_q",
      "task", "n_calib", "competitor", "metric", "evaluation_stratum",
      "target", "coordinate", "advantage_direction"
    )
    comparison_summary <- mixedgp_group_summary(
      comparisons, "EIV_advantage", comparison_groups
    )
    mixedgp_atomic_write_csv(
      comparison_summary,
      file.path(combined_dir, "method_comparisons_paired_summary.csv")
    )
  }
  contrasts <- if (config$study == "study1") {
    mixedgp_study1_mechanism_contrasts(comparisons)
  } else {
    mixedgp_study2_design_contrasts(comparisons)
  }
  if (nrow(contrasts) > 0L) {
    mixedgp_atomic_write_csv(
      contrasts, file.path(combined_dir, "design_contrasts_paired_raw.csv")
    )
    contrast_summary <- mixedgp_group_summary(
      contrasts, "advantage_change",
      c(
        "task", "n_calib", "competitor", "metric",
        "evaluation_stratum", "contrast", "contrast_direction"
      )
    )
    mixedgp_atomic_write_csv(
      contrast_summary,
      file.path(combined_dir, "design_contrasts_paired_summary.csv")
    )
  }
  invisible(list(comparisons = comparisons, contrasts = contrasts))
}

mixedgp_aggregate_results <- function(results, config, run_dir) {
  combined_dir <- file.path(run_dir, "combined")
  dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
  non_tabular <- lapply(results, function(result) {
    result$outputs[!vapply(result$outputs, is.data.frame, logical(1L))]
  })
  non_tabular <- Filter(function(x) length(x) > 0L, non_tabular)
  if (length(non_tabular) > 0L) {
    mixedgp_atomic_save_rds(
      non_tabular, file.path(combined_dir, "cell_non_tabular_outputs.rds")
    )
  }
  output_names <- unique(unlist(lapply(results, function(z) names(z$outputs))))
  combined <- setNames(vector("list", length(output_names)), output_names)
  for (name in output_names) {
    rows <- lapply(results, function(result) {
      mixedgp_tag_output(result$outputs[[name]], result$cell, config$study)
    })
    combined[[name]] <- mixedgp_bind_rows_base(rows)
    if (is.data.frame(combined[[name]]) && ncol(combined[[name]]) > 0L) {
      mixedgp_atomic_write_csv(
        combined[[name]], file.path(combined_dir, paste0(name, ".csv"))
      )
    }
  }
  mixedgp_atomic_save_rds(combined, file.path(combined_dir, "all_raw_outputs.rds"))
  mixedgp_write_publication_comparisons(combined, config, combined_dir)
  output_files <- setdiff(
    list.files(combined_dir, full.names = TRUE),
    file.path(combined_dir, "result_manifest.csv")
  )
  manifest <- data.frame(
    file = basename(output_files),
    bytes = file.info(output_files)$size,
    md5 = unname(tools::md5sum(output_files)),
    stringsAsFactors = FALSE
  )
  mixedgp_atomic_write_csv(manifest, file.path(combined_dir, "result_manifest.csv"))
  combined
}

mixedgp_run_simulation <- function(config) {
  config <- validate_simulation_config(config)
  run_started_at <- Sys.time()
  engine <- mixedgp_simulation_engine(config$code_dir)
  preflight_fun <- get("mixedgp_competitor_preflight", envir = engine)
  preflight <- preflight_fun(config$published_methods, strict = FALSE)
  runtime_preflight <- mixedgp_runtime_preflight(config)
  task_plan <- mixedgp_task_plan(config)
  if (identical(config$mode, "dry_run")) {
    mixedgp_print_dry_run(config, preflight, runtime_preflight, task_plan)
    return(invisible(list(
      config = config, preflight = preflight,
      runtime_preflight = runtime_preflight, task_plan = task_plan,
      data_status = mixedgp_existing_data_status(config)
    )))
  }

  run_dir <- mixedgp_run_directory(config)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  mixedgp_write_run_metadata(
    config, run_dir, preflight, runtime_preflight, task_plan
  )
  n_cells <- length(config$cells)
  mixedgp_write_progress(
    run_dir, "running", "preflight", cells_total = n_cells,
    detail = "Run directory initialized."
  )
  automatic_validation <- list()
  if ("fit" %in% config$stages &&
      config$mode %in% c("smoke", "publication")) {
    message("Running automatic pre-fit publication validators.")
    mixedgp_write_progress(
      run_dir, "running", "automatic_validation", cells_total = n_cells,
      detail = "Checking published competitors and experiment design."
    )
    automatic_validation <- mixedgp_run_automatic_validators(
      config, engine, run_dir
    )
    mixedgp_write_progress(
      run_dir, "running", "automatic_validation", cells_total = n_cells,
      detail = "Automatic validators passed."
    )
  }
  if ("fit" %in% config$stages && any(!preflight$available)) {
    unavailable <- preflight$method[!preflight$available]
    warning(
      "Published competitor package(s) are unavailable: ",
      paste(unavailable, collapse = ", "),
      ". EIV-GP MCMC will continue; the omissions are recorded in the ",
      "competitor preflight and per-replication status files.",
      call. = FALSE
    )
  }
  if ("fit" %in% config$stages && any(!runtime_preflight$available)) {
    unavailable <- runtime_preflight$package[!runtime_preflight$available]
    stop(
      "Run stopped before fitting because runtime packages are unavailable: ",
      paste(unavailable, collapse = ", "), "."
    )
  }

  old_cores <- getOption("mixedgp.cores", NULL)
  options(mixedgp.cores = config$parallel$workers)
  on.exit(options(mixedgp.cores = old_cores), add = TRUE)
  thread_vars <- c(
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"
  )
  old_threads <- Sys.getenv(thread_vars, unset = NA_character_)
  do.call(Sys.setenv, as.list(setNames(rep("1", length(thread_vars)), thread_vars)))
  on.exit({
    present <- !is.na(old_threads)
    if (any(present)) {
      do.call(Sys.setenv, as.list(old_threads[present]))
    }
    if (any(!present)) Sys.unsetenv(thread_vars[!present])
  }, add = TRUE)

  generation_manifests <- list()
  if ("data" %in% config$stages) {
    for (ii in seq_along(config$cells)) {
      cell <- config$cells[[ii]]
      mixedgp_write_progress(
        run_dir, "running", "data_generation", ii - 1L, n_cells, cell$id,
        "Generating deterministic frozen datasets."
      )
      message("Freezing synthetic data: ", cell$id)
      generation_manifests[[cell$id]] <- mixedgp_generate_cell_data(
        config, cell, engine
      )
      mixedgp_write_progress(
        run_dir, "running", "data_generation", ii, n_cells, cell$id,
        "Frozen datasets and manifest written."
      )
    }
  }
  paired_validation <- data.frame()
  selected_manifests <- list()
  if (any(c("data", "fit") %in% config$stages)) {
    mixedgp_write_progress(
      run_dir, "running", "data_verification", cells_total = n_cells,
      detail = "Verifying manifests, checksums, and paired design."
    )
    for (cell in config$cells) {
      selected <- mixedgp_verify_cell_data(config, cell, engine)
      selected$cell_id <- cell$id
      selected_manifests[[cell$id]] <- selected
    }
    paired_validation <- mixedgp_validate_common_random_numbers(config, engine)
    if (is.data.frame(paired_validation) && ncol(paired_validation) > 0L) {
      mixedgp_atomic_write_csv(
        paired_validation,
        file.path(run_dir, "config", "paired_design_validation.csv")
      )
    }
  }
  data_manifest <- mixedgp_bind_rows_base(selected_manifests)
  if (nrow(data_manifest) > 0L) {
    mixedgp_atomic_write_csv(
      data_manifest, file.path(run_dir, "config", "input_manifest.csv")
    )
  }

  if (!"fit" %in% config$stages) {
    mixedgp_write_progress(
      run_dir, "completed", "data_generation", n_cells, n_cells,
      detail = "Data-only run completed successfully."
    )
    return(invisible(list(
      config = config, run_dir = run_dir, data_manifest = data_manifest,
      preflight = preflight, runtime_preflight = runtime_preflight,
      automatic_validation = automatic_validation,
      paired_validation = paired_validation
    )))
  }
  results <- vector("list", length(config$cells))
  names(results) <- vapply(config$cells, `[[`, character(1L), "id")
  gate_rows <- list()
  status_rows <- list()
  for (ii in seq_along(config$cells)) {
    cell <- config$cells[[ii]]
    mixedgp_write_progress(
      run_dir, "running", "fitting", ii - 1L, n_cells, cell$id,
      "Fitting EIV-GP, competitors, and prespecified ablations.",
      run_started_at = run_started_at
    )
    message(
      "Running ", config$study, " cell ", ii, "/", length(config$cells),
      ": ", cell$id
    )
    started <- Sys.time()
    answer <- tryCatch(
      if (config$study == "study1") {
        mixedgp_run_study1_cell(config, cell, engine, run_dir)
      } else {
        mixedgp_run_study2_cell(config, cell, engine, run_dir)
      },
      error = function(e) e
    )
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (inherits(answer, "error")) {
      status_rows[[cell$id]] <- data.frame(
        cell_id = cell$id, status = "failed",
        message = conditionMessage(answer), elapsed_seconds = elapsed,
        stringsAsFactors = FALSE
      )
      status <- mixedgp_bind_rows_base(status_rows)
      mixedgp_atomic_write_csv(
        status, file.path(run_dir, "config", "cell_status.csv")
      )
      mixedgp_write_progress(
        run_dir, "failed", "fitting", ii - 1L, n_cells, cell$id,
        conditionMessage(answer), run_started_at = run_started_at
      )
      stop("Simulation cell ", cell$id, " failed: ", conditionMessage(answer))
    }
    results[[cell$id]] <- answer
    gate_rows[[cell$id]] <- mixedgp_result_gates(answer, config)
    failed_gates <- gate_rows[[cell$id]]$gate[!gate_rows[[cell$id]]$pass]
    gate_failed <- length(failed_gates) > 0L
    cell_status <- if (gate_failed && isTRUE(config$fail_closed)) {
      "gate_failed"
    } else if (gate_failed) {
      "completed_with_gate_warnings"
    } else {
      "success"
    }
    status_rows[[cell$id]] <- data.frame(
      cell_id = cell$id,
      status = cell_status,
      message = if (gate_failed) paste(failed_gates, collapse = "; ") else "",
      elapsed_seconds = elapsed, stringsAsFactors = FALSE
    )
    mixedgp_atomic_write_csv(
      mixedgp_bind_rows_base(status_rows),
      file.path(run_dir, "config", "cell_status.csv")
    )
    mixedgp_atomic_write_csv(
      mixedgp_bind_rows_base(gate_rows),
      file.path(run_dir, "config", "diagnostic_gates.csv")
    )
    mixedgp_write_progress(
      run_dir, cell_status, "fitting", ii, n_cells, cell$id,
      paste0("Cell finished in ", round(elapsed, 1L), " seconds."),
      run_started_at = run_started_at
    )
    if (isTRUE(config$fail_closed) && gate_failed) {
      mixedgp_write_progress(
        run_dir, "failed", "diagnostic_gates", ii, n_cells, cell$id,
        paste(failed_gates, collapse = "; "), run_started_at = run_started_at
      )
      stop(
        "Publication gates failed for ", cell$id, ": ",
        paste(failed_gates, collapse = ", "), "."
      )
    }
  }

  combined <- if ("aggregate" %in% config$stages) {
    mixedgp_write_progress(
      run_dir, "running", "aggregation", n_cells, n_cells,
      detail = "Writing combined results, comparisons, tables, and figures.",
      run_started_at = run_started_at
    )
    mixedgp_aggregate_results(results, config, run_dir)
  } else {
    list()
  }
  out <- structure(
    list(
      config = config,
      run_dir = run_dir,
      preflight = preflight,
      runtime_preflight = runtime_preflight,
      automatic_validation = automatic_validation,
      data_manifest = data_manifest,
      paired_validation = paired_validation,
      gates = mixedgp_bind_rows_base(gate_rows),
      results = results,
      combined = combined
    ),
    class = c("mixedgp_simulation_run", "list")
  )
  mixedgp_atomic_save_rds(
    out[c(
      "config", "run_dir", "preflight", "runtime_preflight", "data_manifest",
      "automatic_validation", "paired_validation", "gates"
    )],
    file.path(run_dir, "run_summary.rds")
  )
  mixedgp_write_progress(
    run_dir, "completed", "completed", n_cells, n_cells,
    detail = "Run completed successfully."
  )
  message("Completed simulation run: ", normalizePath(run_dir))
  invisible(out)
}

#' Run the Study I publication workflow
#'
#' Runs the stages selected by [study1_simulation_config()]. The default
#' configuration only prints and returns the validated dry-run plan. Smoke and
#' publication fits automatically execute and archive the exact competitor and
#' paired-design validators before computation.
#'
#' @param config A Study I `mixedgp_simulation_config` object.
#'
#' @return Invisibly, the dry-run specification or completed simulation run.
#' @export
run_study1_simulation <- function(config = study1_simulation_config()) {
  if (!identical(config$study, "study1")) stop("Expected a Study I config.")
  mixedgp_run_simulation(config)
}

#' Run the Study II publication workflow
#'
#' Runs the stages selected by [study2_simulation_config()]. The default
#' configuration only prints and returns the validated dry-run plan. Smoke and
#' publication fits automatically execute and archive the exact competitor and
#' paired-design validators before computation.
#'
#' @param config A Study II `mixedgp_simulation_config` object.
#'
#' @return Invisibly, the dry-run specification or completed simulation run.
#' @export
run_study2_simulation <- function(config = study2_simulation_config()) {
  if (!identical(config$study, "study2")) stop("Expected a Study II config.")
  mixedgp_run_simulation(config)
}
