############################################################
## Unified Study I / Study II experiment runner
############################################################

mixedgp_find_code_dir <- function(code_dir = NULL) {
  marker <- "01_study1_representative_figures.R"
  if (!is.null(code_dir)) {
    code_dir <- normalizePath(code_dir, winslash = "/", mustWork = TRUE)
    if (!file.exists(file.path(code_dir, marker))) {
      stop("code_dir does not contain the mixed-input GP scripts: ", code_dir)
    }
    return(code_dir)
  }

  installed_dir <- system.file("experiments", package = "eivGP")
  in_package <- isNamespace(environment(mixedgp_find_code_dir))
  if (in_package) {
    if (!nzchar(installed_dir) ||
        !file.exists(file.path(installed_dir, marker))) {
      stop("The installed eivGP experiment scripts are unavailable.")
    }
    return(normalizePath(installed_dir, winslash = "/", mustWork = TRUE))
  }

  candidates <- c(
    getwd(),
    file.path(getwd(), "script", "newly-written-and-modified-scripts"),
    file.path(getwd(), "codes"),
    file.path(getwd(), "revision", "codes"), installed_dir
  )
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(file.path(candidates, marker))]
  if (length(hit) == 0L) {
    stop("Cannot locate the mixed-input GP experiment sources; ",
         "supply code_dir explicitly.")
  }
  normalizePath(hit[1L], winslash = "/", mustWork = TRUE)
}

mixedgp_experiment_option_names <- function(study) {
  study <- match.arg(study, c("study1", "study2"))
  if (study == "study1") {
    return(c(
      "STUDY1_CONFIG", "STUDY1_QUICK", "STUDY1_USE_CACHE",
      "STUDY1_REUSE_LOCKED_EIV", "STUDY1_STRICT_COMPETITORS",
      "STUDY1_PUBLISHED_COMPETITORS", "STUDY1_RUN_ABLATIONS",
      "STUDY1_RUN_TARGETED_CONTROLS", "STUDY1_PARALLEL_LEVEL",
      "STUDY1_OUT_PREFIX", "STUDY1_DATA_DIR", "STUDY1_LVGP_MAX_ELAPSED",
      "STUDY1_MAX_RHAT", "STUDY1_MECHANISM_CALIB", "STUDY1_MIN_ESS",
      "STUDY1_REQUIRE_MCMC_GATE", "STUDY1_CONTROL_CALIB"
    ))
  }
  c(
    "STUDY2_CONFIG", "STUDY2_USE_CACHE", "STUDY2_MC_RESUME",
    "STUDY2_OUT_PREFIX", "STUDY2_DATA_DIR", "STUDY2_RUN_REPRESENTATIVE",
    "STUDY2_RUN_MONTE_CARLO", "STUDY2_STRICT_COMPETITORS",
    "STUDY2_RUN_ABLATIONS", "STUDY2_PARALLEL_LEVEL",
    "STUDY2_PUBLISHED_COMPETITORS", "STUDY2_SCENARIOS",
    "STUDY2_PRIMARY_CALIB_GRID", "STUDY2_CONTRAST_CALIB",
    "STUDY2_ABLATION_SCENARIOS", "STUDY2_ABLATION_GP_MAXIT",
    "STUDY2_ABLATION_GP_N_STARTS", "STUDY2_DIAGNOSTIC_GIBBS_SWEEPS",
    "STUDY2_ENFORCE_MCMC_GATE", "STUDY2_LVGP_MAX_ELAPSED",
    "STUDY2_MCMC_RAW_ESS_LIMIT", "STUDY2_MCMC_RHAT_LIMIT",
    "STUDY2_MCMC_TARGET_BULK_ESS_LIMIT", "STUDY2_MCMC_TARGET_N_POINTS",
    "STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT", "STUDY2_MC_BURN",
    "STUDY2_MC_N_CHAINS", "STUDY2_MC_N_ITER", "STUDY2_MC_N_REP",
    "STUDY2_MC_THIN", "STUDY2_MEAS_BULK_ESS_LIMIT", "STUDY2_MEAS_BURN",
    "STUDY2_MEAS_N_CHAINS", "STUDY2_MEAS_N_ITER",
    "STUDY2_MEAS_RHAT_LIMIT", "STUDY2_MEAS_TAIL_ESS_LIMIT",
    "STUDY2_MEAS_THIN", "STUDY2_PREDICTIVE_LATENT_SAMPLER",
    "STUDY2_REP_BURN", "STUDY2_REP_FIT_CALIBS", "STUDY2_REP_N_CHAINS",
    "STUDY2_REP_N_ITER", "STUDY2_REP_SCENARIO", "STUDY2_REP_THIN",
    "STUDY2_SAVE_PDF", "STUDY2_SAVE_PNG", "STUDY2_SAVE_REP_FITS",
    "STUDY2_MAIN_CALIB"
  )
}

mixedgp_validate_study_options <- function(study_options, study) {
  if (length(study_options) == 0L) return(study_options)
  flag_names <- if (study == "study1") {
    c(
      "STUDY1_QUICK", "STUDY1_USE_CACHE", "STUDY1_REUSE_LOCKED_EIV",
      "STUDY1_STRICT_COMPETITORS", "STUDY1_RUN_ABLATIONS",
      "STUDY1_RUN_TARGETED_CONTROLS", "STUDY1_REQUIRE_MCMC_GATE"
    )
  } else {
    c(
      "STUDY2_USE_CACHE", "STUDY2_MC_RESUME", "STUDY2_RUN_REPRESENTATIVE",
      "STUDY2_RUN_MONTE_CARLO", "STUDY2_STRICT_COMPETITORS",
      "STUDY2_RUN_ABLATIONS", "STUDY2_ENFORCE_MCMC_GATE",
      "STUDY2_SAVE_PDF", "STUDY2_SAVE_PNG", "STUDY2_SAVE_REP_FITS"
    )
  }
  for (name in intersect(names(study_options), flag_names)) {
    study_options[[name]] <- mixedgp_validate_flag(study_options[[name]], name)
  }

  positive_integer_names <- if (study == "study1") {
    character()
  } else {
    c(
      "STUDY2_ABLATION_GP_MAXIT", "STUDY2_ABLATION_GP_N_STARTS",
      "STUDY2_DIAGNOSTIC_GIBBS_SWEEPS", "STUDY2_MCMC_TARGET_N_POINTS",
      "STUDY2_MC_BURN", "STUDY2_MC_N_CHAINS", "STUDY2_MC_N_ITER",
      "STUDY2_MC_N_REP", "STUDY2_MC_THIN", "STUDY2_MEAS_BURN",
      "STUDY2_MEAS_N_CHAINS", "STUDY2_MEAS_N_ITER", "STUDY2_MEAS_THIN",
      "STUDY2_REP_BURN", "STUDY2_REP_N_CHAINS", "STUDY2_REP_N_ITER",
      "STUDY2_REP_THIN"
    )
  }
  for (name in intersect(names(study_options), positive_integer_names)) {
    study_options[[name]] <- mixedgp_as_integer_strict(
      study_options[[name]], name, min_value = 1L, length_expected = 1L
    )
  }
  calibration_names <- if (study == "study1") {
    c("STUDY1_MECHANISM_CALIB", "STUDY1_CONTROL_CALIB")
  } else {
    c(
      "STUDY2_PRIMARY_CALIB_GRID", "STUDY2_CONTRAST_CALIB",
      "STUDY2_REP_FIT_CALIBS", "STUDY2_MAIN_CALIB"
    )
  }
  for (name in intersect(names(study_options), calibration_names)) {
    study_options[[name]] <- mixedgp_as_integer_strict(
      study_options[[name]], name, min_value = 0L
    )
  }

  config_name <- paste0(toupper(study), "_CONFIG")
  if (config_name %in% names(study_options)) {
    study_options[[config_name]] <- match.arg(
      study_options[[config_name]], c("quick", "balanced", "thorough")
    )
  }
  parallel_name <- paste0(toupper(study), "_PARALLEL_LEVEL")
  if (parallel_name %in% names(study_options)) {
    study_options[[parallel_name]] <- match.arg(
      study_options[[parallel_name]], c("chains", "replications", "none")
    )
  }
  competitor_name <- paste0(toupper(study), "_PUBLISHED_COMPETITORS")
  if (competitor_name %in% names(study_options)) {
    values <- study_options[[competitor_name]]
    if (!is.character(values) || anyNA(values) || anyDuplicated(values) ||
        any(!values %in% c("UC-GP", "LVGP", "EzGP"))) {
      stop(competitor_name, " contains an unsupported or duplicate method.")
    }
  }
  if (study == "study2") {
    allowed_scenarios <- c(
      "primary", "latent_additive_control", "high_uncertainty",
      "logistic_misspec"
    )
    for (name in intersect(
      names(study_options), c("STUDY2_SCENARIOS", "STUDY2_ABLATION_SCENARIOS")
    )) {
      values <- study_options[[name]]
      if (!is.character(values) || length(values) < 1L || anyNA(values) ||
          anyDuplicated(values) || any(!values %in% allowed_scenarios)) {
        stop(name, " contains an unsupported or duplicate scenario.")
      }
    }
    if ("STUDY2_PREDICTIVE_LATENT_SAMPLER" %in% names(study_options)) {
      study_options$STUDY2_PREDICTIVE_LATENT_SAMPLER <- match.arg(
        study_options$STUDY2_PREDICTIVE_LATENT_SAMPLER,
        c("minimax_tilting", "rejection", "gibbs")
      )
    }
  }
  study_options
}

mixedgp_experiment_execution_plan <- function(parallel_level, n_cores = NULL) {
  parallel_level <- match.arg(
    parallel_level, c("chains", "replications", "none")
  )
  available <- mixedgp_parallel_backend(n_cores)
  serial <- list(backend = "serial", cores = 1L)
  active <- list(backend = available$backend, cores = available$cores)
  list(
    requested_level = parallel_level,
    requested_cores = available$requested_cores,
    chains = if (parallel_level == "chains") active else serial,
    replications = if (parallel_level == "replications") active else serial,
    data = if (parallel_level == "replications") active else serial,
    os_type = available$os_type
  )
}

mixedgp_experiment_spec <- function(
    study = c("study1", "study2"),
    config = c("quick", "balanced", "thorough"),
    stages = c("data", "representative", "monte_carlo"),
    output_dir = NULL,
    data_dir = NULL,
    competitors = c("UC-GP", "LVGP", "EzGP"),
    parallel_level = c("chains", "replications", "none"),
    n_cores = NULL,
    use_cache = FALSE,
    run_ablations = TRUE,
    strict_competitors = NULL,
    code_dir = NULL,
    study_options = list()) {
  study <- match.arg(study)
  config <- match.arg(config)
  parallel_level <- match.arg(parallel_level)
  use_cache <- mixedgp_validate_flag(use_cache, "use_cache")
  run_ablations <- mixedgp_validate_flag(run_ablations, "run_ablations")
  if (is.null(strict_competitors)) {
    strict_competitors <- identical(config, "thorough")
  } else {
    strict_competitors <- mixedgp_validate_flag(
      strict_competitors, "strict_competitors"
    )
  }
  allowed_stages <- c("data", "representative", "monte_carlo", "controls")
  stages <- unique(as.character(stages))
  unknown_stages <- setdiff(stages, allowed_stages)
  if (length(unknown_stages) > 0L) {
    stop("Unknown stages: ", paste(unknown_stages, collapse = ", "))
  }
  if (study == "study2" && "controls" %in% stages) {
    stop("The separate controls stage applies only to Study I.")
  }
  competitors <- unique(as.character(competitors))
  unknown_competitors <- setdiff(competitors, c("UC-GP", "LVGP", "EzGP"))
  if (length(unknown_competitors) > 0L) {
    stop("Unknown competitors: ", paste(unknown_competitors, collapse = ", "))
  }
  if (!is.list(study_options) || is.null(names(study_options)) &&
      length(study_options) > 0L) {
    stop("study_options must be a named list.")
  }
  if (length(study_options) > 0L &&
      (anyNA(names(study_options)) || any(!nzchar(names(study_options))) ||
       anyDuplicated(names(study_options)))) {
    stop("study_options names must be complete and unique.")
  }
  prefix <- if (study == "study1") "STUDY1_" else "STUDY2_"
  option_names <- names(study_options)
  if (length(option_names) == 0L) option_names <- character()
  bad_options <- option_names[!startsWith(option_names, prefix)]
  if (length(bad_options) > 0L) {
    stop(
      "Study-specific option names must start with ", prefix,
      ". Invalid: ", paste(bad_options, collapse = ", ")
    )
  }
  unknown_options <- setdiff(option_names, mixedgp_experiment_option_names(study))
  if (length(unknown_options) > 0L) {
    stop(
      "Unknown study_options: ", paste(unknown_options, collapse = ", "),
      ". Check the option spelling."
    )
  }
  study_options <- mixedgp_validate_study_options(study_options, study)
  code_dir <- mixedgp_find_code_dir(code_dir)
  required_scripts <- if (study == "study1") {
    c(
      "01_study1_representative_figures.R", "02_study1_monte_carlo.R",
      "02_study1_targeted_controls.R", "04_study1_ablations.R"
    )
  } else {
    c(
      "01_study2_representative_figures.R", "02_study2_monte_carlo.R",
      "04_study2_ablations.R"
    )
  }
  missing_scripts <- required_scripts[
    !file.exists(file.path(code_dir, required_scripts))
  ]
  if (length(missing_scripts) > 0L) {
    stop("Missing study scripts: ", paste(missing_scripts, collapse = ", "))
  }
  if (is.null(output_dir)) output_dir <- getwd()
  output_dir <- normalizePath(
    output_dir, winslash = "/", mustWork = FALSE
  )
  if (is.null(data_dir)) {
    data_dir <- file.path(output_dir, "data-synthetic", study)
  }
  data_dir <- normalizePath(data_dir, winslash = "/", mustWork = FALSE)
  if (!is.null(n_cores)) n_cores <- mixedgp_resolve_cores(n_cores)
  execution <- mixedgp_experiment_execution_plan(parallel_level, n_cores)

  structure(
    list(
      schema_version = "1.0.0",
      study = study,
      config = config,
      stages = stages,
      output_dir = output_dir,
      data_dir = data_dir,
      competitors = competitors,
      parallel_level = parallel_level,
      n_cores = n_cores,
      execution = execution,
      use_cache = use_cache,
      run_ablations = run_ablations,
      strict_competitors = strict_competitors,
      code_dir = code_dir,
      study_options = study_options
    ),
    class = c("mixedgp_experiment_spec", "list")
  )
}

#' @export
#' @noRd
print.mixedgp_experiment_spec <- function(x, ...) {
  cat(
    "Mixed-input GP experiment\n",
    "  study: ", x$study, "\n",
    "  configuration: ", x$config, "\n",
    "  stages: ", paste(x$stages, collapse = ", "), "\n",
    "  parallel level: ", x$parallel_level, "\n",
    "  chain backend/cores: ", x$execution$chains$backend, "/",
    x$execution$chains$cores, "\n",
    "  replication backend/cores: ", x$execution$replications$backend, "/",
    x$execution$replications$cores, "\n",
    "  output: ", x$output_dir, "\n",
    sep = ""
  )
  invisible(x)
}

mixedgp_experiment_environment <- function(spec) {
  run_env <- new.env(parent = parent.env(.GlobalEnv))
  owner_env <- environment(run_mixedgp_experiment)
  if (!identical(owner_env, .GlobalEnv)) {
    owned_names <- ls(owner_env, all.names = TRUE)
    for (name in owned_names) {
      assign(name, get(name, envir = owner_env, inherits = FALSE), run_env)
    }
  }
  run_env$exists <- function(
      x, where = -1, envir = NULL,
      frame, mode = "any", inherits = FALSE) {
    if (!missing(frame)) {
      envir <- sys.frame(frame)
    } else if (is.null(envir)) {
      envir <- if (identical(where, -1)) parent.frame() else as.environment(where)
    }
    base::exists(x, envir = envir, mode = mode, inherits = inherits)
  }
  code_dir <- spec$code_dir
  run_env$source <- function(file, ...) {
    target <- if (grepl("^(/|[A-Za-z]:)", file)) {
      file
    } else {
      file.path(code_dir, file)
    }
    sys.source(target, envir = run_env, chdir = TRUE)
    invisible(target)
  }
  run_env$library <- function(package,
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
      stop("library() requires one package name in experiment scripts.")
    }
    namespace <- loadNamespace(package_name)
    exports <- getNamespaceExports(namespace)
    for (name in exports) {
      assign(name, getExportedValue(package_name, name), envir = run_env)
    }
    if (isTRUE(logical.return)) return(TRUE)
    invisible(namespace)
  }
  run_env
}

mixedgp_assign_experiment_controls <- function(spec, run_env) {
  common <- if (spec$study == "study1") {
    list(
      STUDY1_CONFIG = spec$config,
      STUDY1_QUICK = identical(spec$config, "quick"),
      STUDY1_USE_CACHE = spec$use_cache,
      STUDY1_REUSE_LOCKED_EIV = FALSE,
      STUDY1_STRICT_COMPETITORS = spec$strict_competitors,
      STUDY1_PUBLISHED_COMPETITORS = spec$competitors,
      STUDY1_RUN_ABLATIONS = spec$run_ablations,
      STUDY1_RUN_TARGETED_CONTROLS = "controls" %in% spec$stages,
      STUDY1_PARALLEL_LEVEL = spec$parallel_level,
      STUDY1_OUT_PREFIX = spec$output_dir,
      STUDY1_DATA_DIR = spec$data_dir
    )
  } else {
    list(
      STUDY2_CONFIG = spec$config,
      STUDY2_USE_CACHE = spec$use_cache,
      STUDY2_MC_RESUME = TRUE,
      STUDY2_OUT_PREFIX = spec$output_dir,
      STUDY2_DATA_DIR = spec$data_dir,
      STUDY2_RUN_REPRESENTATIVE = "representative" %in% spec$stages,
      STUDY2_RUN_MONTE_CARLO = "monte_carlo" %in% spec$stages,
      STUDY2_STRICT_COMPETITORS = spec$strict_competitors,
      STUDY2_RUN_ABLATIONS = spec$run_ablations,
      STUDY2_PARALLEL_LEVEL = spec$parallel_level,
      STUDY2_PUBLISHED_COMPETITORS = spec$competitors,
      STUDY2_SCENARIOS = c(
        "primary", "latent_additive_control", "high_uncertainty",
        "logistic_misspec"
      ),
      STUDY2_PRIMARY_CALIB_GRID = c(0L, 10L, 25L, 50L, 80L),
      STUDY2_CONTRAST_CALIB = 25L,
      STUDY2_ABLATION_SCENARIOS = "primary"
    )
  }
  duplicate <- intersect(names(common), names(spec$study_options))
  if (length(duplicate) > 0L) common[duplicate] <- NULL
  list2env(c(common, spec$study_options), envir = run_env)
  invisible(run_env)
}

mixedgp_experiment_data_cores <- function(spec) {
  as.integer(spec$execution$data$cores)
}

mixedgp_generate_experiment_data <- function(spec) {
  data_cores <- mixedgp_experiment_data_cores(spec)
  settings <- if (spec$study == "study1") {
    switch(
      spec$config,
      quick = list(n_rep = 2L, n_test = 150L),
      balanced = list(n_rep = 10L, n_test = 300L),
      thorough = list(n_rep = 50L, n_test = 500L)
    )
  } else {
    study2_config_settings(spec$config)[c("n_rep", "n_test")]
  }
  if (spec$study == "study1") {
    return(generate_study1_synthetic_datasets(
      n_rep = settings$n_rep,
      directory = spec$data_dir,
      n_cores = data_cores,
      overwrite = FALSE,
      n = 100L,
      n_test = settings$n_test,
      m = 6L,
      scenario = "active",
      threshold_design = "imbalanced",
      calib_grid = c(0L, 5L, 10L, 20L, 50L)
    ))
  }
  scenarios <- spec$study_options$STUDY2_SCENARIOS
  if (is.null(scenarios)) {
    scenarios <- c(
      "primary", "latent_additive_control", "high_uncertainty",
      "logistic_misspec"
    )
  }
  n_rep <- spec$study_options$STUDY2_MC_N_REP
  if (is.null(n_rep)) n_rep <- settings$n_rep
  primary_grid <- spec$study_options$STUDY2_PRIMARY_CALIB_GRID
  if (is.null(primary_grid)) primary_grid <- c(0L, 10L, 25L, 50L, 80L)
  contrast_calib <- spec$study_options$STUDY2_CONTRAST_CALIB
  if (is.null(contrast_calib)) contrast_calib <- 25L
  generate_study2_synthetic_datasets(
    n_rep = as.integer(n_rep),
    scenarios = scenarios,
    directory = spec$data_dir,
    n_cores = data_cores,
    overwrite = FALSE,
    n = 120L,
    n_test = settings$n_test,
    q = 4L,
    m = 4L,
    calib_grid = sort(unique(as.integer(c(primary_grid, contrast_calib))))
  )
}

#' Run either numerical study through one interface
#'
#' The orchestration is common across the two studies. Study-specific data
#' generators, EIV-GP engines, estimands, and reporting code remain separate.
#'
#' @param study `"study1"` or `"study2"`.
#' @param config `"quick"`, `"balanced"`, or `"thorough"`.
#' @param stages Any of `"data"`, `"representative"`, `"monte_carlo"`, and
#'   (Study I only) `"controls"`.
#' @param output_dir Root directory for figures, tables, results, and data.
#' @param data_dir Optional frozen-data directory.
#' @param competitors Published competitors to fit.
#' @param parallel_level Parallelize `"chains"`, `"replications"`, or neither.
#' @param n_cores Number of workers.
#' @param use_cache Whether validated fitted-result caches may be reused.
#' @param run_ablations Whether to run the appendix ablations.
#' @param strict_competitors Whether missing or failed competitor packages stop.
#'   `NULL` (the default) enables strictness for `config="thorough"` and
#'   disables it otherwise.
#' @param code_dir Location of the manuscript experiment scripts.
#' @param study_options Named study-prefixed overrides, such as
#'   `list(STUDY2_MC_N_REP = 2L)`.
#' @param dry_run Return the validated specification without doing work.
#' @return A `mixedgp_experiment_run`, or a specification for a dry run.
#' @export
run_mixedgp_experiment <- function(
    study = c("study1", "study2"),
    config = c("quick", "balanced", "thorough"),
    stages = c("data", "representative", "monte_carlo"),
    output_dir = NULL,
    data_dir = NULL,
    competitors = c("UC-GP", "LVGP", "EzGP"),
    parallel_level = c("chains", "replications", "none"),
    n_cores = NULL,
    use_cache = FALSE,
    run_ablations = TRUE,
    strict_competitors = NULL,
    code_dir = NULL,
    study_options = list(),
    dry_run = FALSE) {
  dry_run <- mixedgp_validate_flag(dry_run, "dry_run")
  spec <- mixedgp_experiment_spec(
    study = study, config = config, stages = stages,
    output_dir = output_dir, data_dir = data_dir,
    competitors = competitors, parallel_level = parallel_level,
    n_cores = n_cores, use_cache = use_cache,
    run_ablations = run_ablations,
    strict_competitors = strict_competitors, code_dir = code_dir,
    study_options = study_options
  )
  if (isTRUE(dry_run)) return(spec)

  caller_rng <- mixedgp_rng_state()
  on.exit(mixedgp_restore_rng_state(caller_rng), add = TRUE)
  search_before <- search()
  on.exit({
    newly_attached <- setdiff(search(), search_before)
    for (name in newly_attached) {
      try(detach(name, character.only = TRUE, unload = FALSE), silent = TRUE)
    }
  }, add = TRUE)
  old_cores <- getOption("mixedgp.cores", NULL)
  options(mixedgp.cores = spec$execution$requested_cores)
  on.exit(options(mixedgp.cores = old_cores), add = TRUE)
  dir.create(spec$output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(spec$data_dir, recursive = TRUE, showWarnings = FALSE)

  run_env <- mixedgp_experiment_environment(spec)
  mixedgp_assign_experiment_controls(spec, run_env)
  manifest <- NULL
  if ("data" %in% spec$stages) {
    manifest <- mixedgp_generate_experiment_data(spec)
  }

  fitting_stages <- intersect(
    spec$stages, c("representative", "monte_carlo", "controls")
  )
  preflight <- if (length(fitting_stages) > 0L) {
    mixedgp_competitor_preflight(
      spec$competitors, strict = spec$strict_competitors
    )
  } else {
    data.frame()
  }
  if (spec$study == "study1" && length(fitting_stages) > 0L) {
    run_env$source("04_study1_ablations.R")
    if ("representative" %in% spec$stages) {
      run_env$source("01_study1_representative_figures.R")
    }
    if ("monte_carlo" %in% spec$stages) {
      run_env$source("02_study1_monte_carlo.R")
    }
    if ("controls" %in% spec$stages) {
      run_env$source("02_study1_targeted_controls.R")
    }
  } else if (spec$study == "study2" && length(fitting_stages) > 0L) {
    run_env$source("04_study2_ablations.R")
    if ("representative" %in% spec$stages) {
      run_env$source("01_study2_representative_figures.R")
    }
    if ("monte_carlo" %in% spec$stages) {
      run_env$source("02_study2_monte_carlo.R")
    }
  }

  result <- structure(
    list(
      spec = spec,
      execution = spec$execution,
      manifest = manifest,
      competitor_preflight = preflight,
      design_tag = get0(
        if (spec$study == "study1") "STUDY1_DESIGN_TAG" else "STUDY2_DESIGN_TAG",
        envir = run_env, inherits = FALSE
      ),
      cache_tag = get0("CACHE_TAG", envir = run_env, inherits = FALSE),
      results = get0(
        if (spec$study == "study1") "mc_results" else "raw_outputs",
        envir = run_env, inherits = FALSE
      )
    ),
    class = c("mixedgp_experiment_run", "list")
  )
  result
}
