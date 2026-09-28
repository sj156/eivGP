# Simulation design: loaded by codes/simulation_helpers.R.

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
  match.arg(mode, c("dry_run", "data", "smoke", "development", "publication"))
}

mixedgp_development_profile <- function(config) {
  integer_env <- function(name, default, minimum = 1L) {
    value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
    mixedgp_validate_scalar_integer(value, name, minimum)
  }
  selected <- vapply(config$cells, `[[`, character(1L), "id")
  requested <- Sys.getenv("MIXEDGP_DEV_CELLS", unset = "")
  ids <- vapply(config$cells, `[[`, character(1L), "id")
  if (nzchar(requested)) selected <- trimws(strsplit(requested, ",", fixed = TRUE)[[1L]])
  if (identical(selected, "all")) selected <- ids
  if (anyDuplicated(selected) || any(!selected %in% ids)) {
    stop("MIXEDGP_DEV_CELLS must name unique cells from: ", paste(ids, collapse = ", "))
  }
  config$cells <- config$cells[match(selected, ids)]
  reps <- integer_env("MIXEDGP_DEV_REPS", 3L)
  burn <- integer_env("MIXEDGP_DEV_BURN", 500L, 0L)
  draws <- integer_env("MIXEDGP_DEV_DRAWS", 1250L)
  config$cells <- lapply(config$cells, function(cell) {
    cell$n_rep <- reps
    cell
  })
  config$run_id <- paste0(config$study, "-development")
  config$stages <- c("data", "fit", "aggregate")
  config$strict_competitors <- FALSE
  config$fail_closed <- FALSE
  config$mcmc$n_iter <- burn + draws
  config$mcmc$burn <- burn
  config$mcmc$thin <- 1L
  config$mcmc$n_chains <- 4L
  config$mcmc$require_gate <- FALSE
  config$mcmc$automatic_continuation <- FALSE
  config$measurement_mcmc$n_iter <- burn + draws
  config$measurement_mcmc$burn <- burn
  if (config$study == "study2") config$measurement_mcmc$thin <- 1L
  config$evaluation$n_pred_draw <- 100L
  config$evaluation$n_m_eval <- 30L
  config$evaluation$n_m_draw <- 50L
  config$evaluation$n_m_latent <- 64L
  if (config$study == "study2") {
    config$evaluation$n_m_truth <- 500L
    config$evaluation$n_oracle_pool <- 10000L
  }
  config$development_note <- paste(
    "DEVELOPMENT RESULTS: reduced replication and Monte Carlo budgets.",
    "For draft figures and workflow inspection; review diagnostics and rerun",
    "the frozen publication design before reporting scientific conclusions."
  )
  config
}

mixedgp_apply_core_budget <- function(config, core_budget = NULL) {
  if (is.null(core_budget)) {
    value <- Sys.getenv("MIXEDGP_CORE_BUDGET", unset = "")
    if (!nzchar(value)) return(config)
    core_budget <- suppressWarnings(as.numeric(value))
  }
  allocation_env <- new.env(parent = baseenv())
  sys.source(file.path(config$code_dir, "core/00_parallel_utils.R"), envir = allocation_env)
  allocation <- allocation_env$eivgp_run_settings(
    core_budget, max(vapply(config$cells, `[[`, integer(1), "n_rep")),
    chains = config$mcmc$n_chains,
    sampling_iterations = config$mcmc$n_iter - config$mcmc$burn,
    warmup = config$mcmc$burn)
  config$parallel <- allocation[c("level", "workers", "chain_workers", "core_budget", "active_chain_limit")]
  config
}

mixedgp_number_slug <- function(x) {
  gsub("[^0-9A-Za-z]+", "p", format(x, scientific = FALSE, trim = TRUE))
}

mixedgp_config_fingerprint <- function(config) {
  payload <- config
  payload$created_at <- NULL
  payload$stages <- NULL
  code_files <- unique(c(
    mixedgp_simulation_modules(), "simulation_helpers.R",
    "simulations/02_study1_monte_carlo.R", "simulations/02_study2_monte_carlo.R",
    "simulations/setup_study1_experiment.R", "simulations/setup_study2_experiment.R",
    "simulations/run_study1_simulation.R", "simulations/run_study2_simulation.R"
  ))
  code_files <- setdiff(code_files, c("core/03_study2_published_competitors.R", "simulations/competitor_cache.R"))
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
    ggplot2 = "ggplot2", dplyr = "dplyr", tidyr = "tidyr",
    knitr = "knitr", posterior = "posterior"
  )
  if (identical(config$study, "study2")) {
    package_names <- c(package_names, TruncatedNormal = "TruncatedNormal")
  }
  payload$runtime_package_identity <- vapply(
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
  smoke <- identical(mode, "smoke")
  cells <- list()
  for (threshold_design in "balanced") {
    for (eta in c(0, 1)) {
      cells[[length(cells) + 1L]] <- list(
        id = paste0("eta", eta, "_", threshold_design),
        role = "primary_heterogeneity",
        scenario = "heterogeneity_continuum",
        heterogeneity_eta = eta, threshold_design = threshold_design,
        min_class_count = if (threshold_design == "imbalanced") 3L else 0L,
        n = 100L, n_test = if (smoke) 80L else 100L,
        n_rep = if (smoke) 1L else 50L, m = 6L,
        calibration_grid = if (smoke) c(0L, 5L) else c(0L, 20L, 50L),
        evaluate_f = identical(eta, 1), evaluate_u = TRUE,
        run_ablations = TRUE)
    }
  }
  cells
}

mixedgp_study2_cells <- function(mode) {
  smoke <- identical(mode, "smoke")
  make_cell <- function(id, role, scenario, q, grid,
                        run_ablations, evaluate_f = FALSE, evaluate_u = TRUE) {
    list(id = id, role = role, scenario = scenario, q = as.integer(q),
         d = 2L, m = 4L, n = 100L, n_test = if (smoke) 60L else 100L,
         n_rep = if (smoke) 1L else 50L, calibration_grid = as.integer(grid),
         run_ablations = run_ablations, evaluate_f = evaluate_f,
         evaluate_u = evaluate_u)
  }
  list(
    make_cell("primary_q2", "correct_measurement_specification", "primary", 2L,
              if (smoke) 6L else 50L, TRUE),
    make_cell("primary_q4_calibration", "primary_calibration_curve", "primary", 4L,
              if (smoke) c(0L, 6L) else c(0L, 20L, 50L, 80L), TRUE,
              evaluate_f = TRUE),
    make_cell("logistic_q4", "measurement_misspecification", "logistic_misspec", 4L,
              if (smoke) 6L else 50L, FALSE, evaluate_u = FALSE)
  )
}

study1_simulation_config <- function(
    mode = "dry_run",
    code_dir = NULL,
    workers = mixedgp_available_workers(),
    output_root = NULL,
    data_root = NULL, core_budget = NULL) {
  mode <- mixedgp_simulation_mode(mode)
  code_dir <- mixedgp_simulation_code_dir(code_dir)
  revision_dir <- dirname(code_dir)
  if (mode == "development") {
    if (is.null(output_root)) output_root <- file.path(revision_dir, "reproduction", "development", "results", "study1")
    if (is.null(data_root)) data_root <- file.path(revision_dir, "reproduction", "development", "data", "study1")
  }
  if (is.null(output_root)) {
    output_root <- file.path(revision_dir, "simulation-runs", "study1")
  }
  if (is.null(data_root)) {
    data_root <- file.path(
      revision_dir, "data-synthetic", "publication-v2", "study1"
    )
  }
  smoke <- identical(mode, "smoke")
  config <- list(
    schema_version = MIXEDGP_SIMULATION_SCHEMA,
    study = "study1",
    mode = mode,
    run_id = paste0("study1-", mode, "-v2"),
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
    published_methods = character(),
    strict_competitors = identical(mode, "publication"),
    fail_closed = identical(mode, "publication"),
    use_cache = TRUE,
    parallel = list(level = "replications", workers = as.integer(workers)),
    mcmc = if (smoke) {
      list(n_iter = 120L, burn = 40L, thin = 1L, n_chains = 1L,
           rhat_limit = 1.05, ess_limit = 20L, require_gate = FALSE)
    } else {
      list(n_iter = if (mode == "publication") 20000L else 5000L,
           burn = if (mode == "publication") 5000L else 1000L,
           thin = 1L, n_chains = 4L,
           rhat_limit = 1.01, ess_limit = 400L,
           require_gate = identical(mode, "publication"))
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
  class(config) <- c("mixedgp_simulation_config", "list")
  if (mode == "development") config <- mixedgp_development_profile(config)
  config <- mixedgp_apply_core_budget(config, core_budget)
  validate_simulation_config(config)
}

study2_simulation_config <- function(
    mode = "dry_run",
    code_dir = NULL,
    workers = mixedgp_available_workers(),
    output_root = NULL,
    data_root = NULL, core_budget = NULL) {
  mode <- mixedgp_simulation_mode(mode)
  code_dir <- mixedgp_simulation_code_dir(code_dir)
  revision_dir <- dirname(code_dir)
  if (mode == "development") {
    if (is.null(output_root)) output_root <- file.path(revision_dir, "reproduction", "development", "results", "study2")
    if (is.null(data_root)) data_root <- file.path(revision_dir, "reproduction", "development", "data", "study2")
  }
  if (is.null(output_root)) {
    output_root <- file.path(revision_dir, "simulation-runs", "study2")
  }
  if (is.null(data_root)) {
    data_root <- file.path(
      revision_dir, "data-synthetic", "publication-v2", "study2"
    )
  }
  smoke <- identical(mode, "smoke")
  config <- list(
    schema_version = MIXEDGP_SIMULATION_SCHEMA,
    study = "study2",
    mode = mode,
    run_id = paste0("study2-", mode, "-v2"),
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
    published_methods = character(),
    strict_competitors = identical(mode, "publication"),
    fail_closed = identical(mode, "publication"),
    use_cache = TRUE,
    parallel = list(level = "replications", workers = as.integer(workers)),
    predictive_latent_sampler = "minimax_tilting",
    mcmc = if (smoke) {
      list(n_iter = 120L, burn = 40L, thin = 1L, n_chains = 1L,
           rhat_limit = 1.05, raw_ess_limit = 20L,
           target_bulk_ess_limit = 10L, target_tail_ess_limit = 10L,
           require_gate = FALSE)
    } else {
      list(n_iter = if (mode == "publication") 20000L else 4000L,
           burn = if (mode == "publication") 5000L else 1000L,
           thin = 1L, n_chains = 4L,
           rhat_limit = 1.01, raw_ess_limit = 400L,
           target_bulk_ess_limit = 400L, target_tail_ess_limit = 400L,
           require_gate = identical(mode, "publication"))
    },
    measurement_mcmc = if (smoke) {
      list(n_iter = 120L, burn = 40L, thin = 1L, n_chains = 2L)
    } else {
      list(n_iter = 4000L, burn = 1000L, thin = 1L, n_chains = 4L)
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
  class(config) <- c("mixedgp_simulation_config", "list")
  if (mode == "development") config <- mixedgp_development_profile(config)
  config <- mixedgp_apply_core_budget(config, core_budget)
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
  if (length(config$published_methods) && !identical(config$published_methods, MIXEDGP_PUBLISHED_METHODS)) {
    stop(
      "The publication competitor set is frozen as: ",
      paste(MIXEDGP_PUBLISHED_METHODS, collapse = ", "), "."
    )
  }
  config$parallel$workers <- mixedgp_validate_scalar_integer(
    config$parallel$workers, "parallel$workers", 1L
  )
  if (!config$parallel$level %in% c("replications", "hybrid")) {
    stop("parallel$level must be replications or hybrid.")
  }
  if (config$parallel$level == "hybrid") {
    budget <- mixedgp_validate_scalar_integer(config$parallel$core_budget, "core_budget", 1L)
    chain_workers <- mixedgp_validate_scalar_integer(config$parallel$chain_workers, "chain_workers", 1L)
    if (config$parallel$workers * chain_workers > budget) {
      stop("Dataset workers times chain workers exceeds core_budget.")
    }
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
  ## Diagnostics assess reliability; they do not terminate either run mode.
  if (config$mode %in% c("publication", "development")) {
    config$strict_competitors <- FALSE
    config$fail_closed <- FALSE
    config$mcmc$require_gate <- FALSE
  }
  config
}

