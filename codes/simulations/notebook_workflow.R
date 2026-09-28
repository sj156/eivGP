# Reader-facing profiles share the canonical experiment constructors.
replication_config <- function(study, profile, core_budget, repository_dir,
                               artifact_root = file.path(repository_dir, "reproduction")) {
  study <- match.arg(study, c("study1", "study2"))
  profile <- match.arg(profile, c("quick", "development", "paper"))
  mode <- c(quick = "smoke", development = "development", paper = "publication")[[profile]]
  if (profile == "development" && any(nzchar(Sys.getenv(c(
      "MIXEDGP_DEV_CELLS", "MIXEDGP_DEV_REPS", "MIXEDGP_DEV_BURN", "MIXEDGP_DEV_DRAWS"))))) {
    stop("Unset MIXEDGP_DEV_* overrides to use the documented development profile.")
  }
  data_root <- switch(profile,
    paper = file.path(artifact_root, "data", "synthetic", study),
    quick = file.path(artifact_root, "data", "smoke", study),
    development = file.path(artifact_root, "development", "data", study))
  output_root <- if (profile == "paper") file.path(artifact_root, "results", study) else
    file.path(artifact_root, profile, "results", study)
  config <- get(paste0(study, "_simulation_config"))(
    mode = mode, code_dir = file.path(repository_dir, "codes"),
    data_root = data_root, output_root = output_root, core_budget = core_budget)
  if (profile == "paper") {
    stopifnot(config$mcmc$n_chains == 4L, config$mcmc$n_iter == 20000L,
      config$mcmc$burn == 5000L, all(vapply(config$cells, `[[`, 1L, "n_rep") == 50L))
  }
  config
}

replication_data_figures <- function(config, artifact_root) {
  figure_env <- new.env(parent = environment())
  figure_env$DATA_ROOT <- config$data_root
  figure_env$ARTIFACT_ROOT <- artifact_root
  figure_env$params <- list(figure_root = file.path(artifact_root, "figures", config$mode))
  sys.source(file.path(config$code_dir, "reporting", paste0(config$study, "_data_figures.R")),
             envir = figure_env)
  invisible(NULL)
}

replication_config_table <- function(config) {
  cells <- mixedgp_cell_summary(config)
  cells$chains <- config$mcmc$n_chains
  cells$iterations_per_chain <- config$mcmc$n_iter
  cells$warmup_per_chain <- config$mcmc$burn
  cells
}
