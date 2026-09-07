## Stand-alone reporting. Never sources a Monte Carlo driver or invokes a fit.
mixedgp_report_run <- function(run_dir, action = c("report", "summarize", "plot"),
                               code_dir, output_dir = file.path(run_dir, "reporting")) {
  action <- match.arg(action)
  run_dir <- normalizePath(run_dir, mustWork = TRUE)
  config <- readRDS(file.path(run_dir, "config", "resolved_config.rds"))
  code_dir <- normalizePath(code_dir, mustWork = TRUE)
  engine <- mixedgp_simulation_engine(code_dir)
  rows <- list()
  collected <- list()
  for (cell in config$cells) {
    cell_dir <- file.path(run_dir, "cells", cell$id)
    checkpoint <- file.path(cell_dir, "report_inputs.rds")
    legacy <- list.files(file.path(cell_dir, "results"),
      pattern = "^study2_results_.*[.]rds$", recursive = TRUE, full.names = TRUE)
    if (!file.exists(checkpoint) && !(config$study == "study2" && length(legacy) == 1L)) {
      rows[[cell$id]] <- data.frame(cell_id = cell$id, status = "missing_inputs",
        message = "No unique saved reporting input; no fitting was attempted.")
      next
    }
    target <- file.path(output_dir, cell$id)
    dir.create(target, recursive = TRUE, showWarnings = FALSE)
    outcome <- tryCatch({
      env <- new.env(parent = engine)
      controls <- get(paste0("mixedgp_cell_controls_", config$study))(config, cell, target)
      list2env(controls, env)
      ## Setup supplies labels, palettes and helpers, but its bookkeeping writes
      ## must not replace the original fit provenance during report replay.
      env$saveRDS <- env$writeLines <- env$write.csv <- function(...) invisible(NULL)
      env$capture.output <- function(..., file = NULL) {
        if (is.null(file)) utils::capture.output(...) else invisible(NULL)
      }
      sys.source(file.path(code_dir, paste0("setup_", config$study, "_experiment.R")),
                 envir = env, chdir = TRUE)
      if (file.exists(checkpoint)) {
        saved <- readRDS(checkpoint)
        if (!identical(saved$schema, 1L) || !identical(saved$study, config$study))
          stop("Incompatible reporting checkpoint.")
        list2env(saved$context, env)
      } else {
        ## Recover pre-separation Study II bundles without recomputing fits.
        raw <- readRDS(legacy)
        mapping <- c(mc_results = "predictive_metrics", mc_mean_recovery = "mean_recovery",
          mc_mean_truth_rejection = "mean_truth_rejection", mc_mean_truth_diagnostics = "mean_truth_diagnostics",
          mc_diagnostics = "mcmc_diagnostics", mc_target_diagnostics = "mcmc_target_diagnostics",
          mc_imputation = "latent_imputation", mc_imputation_status = "latent_imputation_status",
          mc_surface = "surface_recovery", mc_ablation_results = "ablation_predictive_metrics",
          mc_ablation_surface = "ablation_surface_recovery")
        for (name in names(mapping)) assign(name, raw[[mapping[[name]]]], env)
        for (name in c("measurement_diagnostics", "measurement_parameter_diagnostics",
          "ablation_optimizer_attempts", "ablation_status", "competitor_status",
          "pattern_counts", "sampler_control_manifest", "replication_metadata", "run_grid"))
          assign(name, raw[[name]], env)
        env$raw_outputs <- raw
        env$CACHE_TAG <- raw$cache_tag
        env$CACHE_SPEC <- raw$cache_spec
      }
      outputs <- if (config$study == "study2") env$raw_outputs else {
        mapping <- c(predictive_metrics="mc_results", mean_recovery="mean_recovery",
          latent_imputation="latent_imputation", surface_recovery="surface_recovery",
          ablation_predictive_metrics="ablation_results", ablation_surface_recovery="ablation_surface_recovery",
          mcmc_diagnostics="mcmc_diagnostics", competitor_status="competitor_status",
          ablation_status="ablation_status", sampler_control_manifest="sampler_control_manifest")
        setNames(lapply(mapping, function(name) get(name, env)), names(mapping))
      }
      collected[[cell$id]] <- list(cell=cell, outputs=outputs)
      ## All generated files go to the reporting tree, never the fitted cache.
      env$FIG_DIR <- file.path(target, "figures")
      env$TAB_DIR <- file.path(target, "tables")
      env$RES_DIR <- file.path(target, "metadata")
      assign(paste0(toupper(config$study), "_OUT_PREFIX"), target, env)
      for (dir in c(env$FIG_DIR, env$TAB_DIR, env$RES_DIR)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      env$write.csv <- if (action == "plot") function(...) invisible(NULL) else utils::write.csv
      if (action != "plot" && config$study == "study1")
        mixedgp_write_diagnostic_tables(env$mcmc_diagnostics, env$mcmc_parameter_diagnostics,
          env$TAB_DIR, "study1", overwrite = TRUE)
      if (action != "plot" && config$study == "study2")
        mixedgp_write_diagnostic_tables(env$mc_diagnostics, env$mc_target_diagnostics,
          env$TAB_DIR, "study2", overwrite = TRUE)
      env$writeLines <- if (action == "plot") function(...) invisible(NULL) else base::writeLines
      env$capture.output <- function(..., file = NULL) {
        if (action == "plot" && !is.null(file)) return(invisible(NULL))
        utils::capture.output(..., file = file)
      }
      plot_status <- list()
      env$ggsave <- function(filename, plot, ...) {
        if (action == "summarize") return(invisible(NULL))
        if (is.data.frame(plot$data) && nrow(plot$data) == 0L) {
          plot_status[[filename]] <<- data.frame(file = basename(filename),
            status = "skipped_empty", message = "No observations for this plot in the saved setting.")
          message("Skipping empty plot: ", basename(filename))
          return(invisible(NULL))
        }
        ggplot2::ggsave(filename, plot, ...)
        plot_status[[filename]] <<- data.frame(file = basename(filename), status = "written", message = "")
      }
      sys.source(file.path(code_dir, paste0("report_", config$study, "_results.R")), envir = env)
      if (length(plot_status)) utils::write.csv(do.call(rbind, plot_status),
        file.path(target, "plot_status.csv"), row.names = FALSE)
      utils::write.csv(data.frame(source_run = run_dir, source_cell = cell$id,
        source_input = if (file.exists(checkpoint)) checkpoint else legacy,
        source_md5 = unname(tools::md5sum(if (file.exists(checkpoint)) checkpoint else legacy)),
        action = action), file.path(target, paste0(action, "_provenance.csv")), row.names = FALSE)
      report_sources <- file.path(code_dir, c("experiment_reporting.R",
        paste0("setup_", config$study, "_experiment.R"), paste0("report_", config$study, "_results.R")))
      utils::write.csv(data.frame(file = basename(report_sources),
        md5 = unname(tools::md5sum(report_sources))),
        file.path(target, paste0(action, "_code_hashes.csv")), row.names = FALSE)
      NULL
    }, error = identity)
    rows[[cell$id]] <- data.frame(cell_id = cell$id,
      status = if (inherits(outcome, "error")) "report_failed" else "success",
      message = if (inherits(outcome, "error")) conditionMessage(outcome) else "")
  }
  if (action != "plot" && length(collected)) {
    combined_error <- tryCatch({
      mixedgp_aggregate_results(collected, config, output_dir)
      NULL
    }, error=identity)
    if (inherits(combined_error,"error")) rows[["combined"]] <- data.frame(
      cell_id="combined",status="report_failed",message=conditionMessage(combined_error))
  }
  status <- do.call(rbind, rows)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(status, file.path(output_dir, paste0(action, "_status.csv")), row.names = FALSE)
  if (any(status$status != "success")) warning("Reporting incomplete; inspect ",
    file.path(output_dir, paste0(action, "_status.csv")), ". Saved fits were not changed.", call. = FALSE)
  invisible(status)
}
