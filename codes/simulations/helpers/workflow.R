# Simulation workflow: loaded by codes/simulation_helpers.R.

mixedgp_run_simulation <- function(config) {
  config <- validate_simulation_config(config)
  engine <- mixedgp_simulation_engine(config$code_dir)
  preflight_fun <- get("mixedgp_competitor_preflight", envir = engine)
  preflight <- preflight_fun(config$published_methods, strict = FALSE)
  runtime_preflight <- mixedgp_runtime_preflight(config)
  if (any(!preflight$available) || any(!runtime_preflight$available)) {
    message("Install missing experiment dependencies explicitly from the repository root:\n",
      "  Rscript --vanilla codes/cli/install_eivgp_dependencies.R\n",
      "Use the same Rscript and MIXEDGP_R_LIBRARY setting for setup and execution.\n",
      "No packages are installed automatically. Existing missing-method caches are not repaired by installation.")
  }
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
  if (identical(config$mode, "development")) {
    writeLines(config$development_note, file.path(run_dir, "DEVELOPMENT_RESULTS.txt"))
    message(config$development_note)
  }
  if ("fit" %in% config$stages && any(!preflight$available)) {
    unavailable <- preflight$method[!preflight$available]
    warning(
      "Competitor packages unavailable: ", paste(unavailable, collapse = ", "),
      ". Their missing results will be recorded; other methods may proceed.",
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
    for (cell in config$cells) {
      message("Freezing synthetic data: ", cell$id)
      generation_manifests[[cell$id]] <- mixedgp_generate_cell_data(
        config, cell, engine
      )
    }
  }
  paired_validation <- data.frame()
  selected_manifests <- list()
  if (any(c("data", "fit") %in% config$stages)) {
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
    return(invisible(list(
      config = config, run_dir = run_dir, data_manifest = data_manifest,
      preflight = preflight, runtime_preflight = runtime_preflight,
      paired_validation = paired_validation
    )))
  }
  results <- vector("list", length(config$cells))
  names(results) <- vapply(config$cells, `[[`, character(1L), "id")
  gate_rows <- list()
  status_rows <- list()
  for (ii in seq_along(config$cells)) {
    cell <- config$cells[[ii]]
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
      stop("Simulation cell ", cell$id, " failed: ", conditionMessage(answer))
    }
    results[[cell$id]] <- answer
    gate_rows[[cell$id]] <- mixedgp_result_gates(answer, config)
    failed_gates <- gate_rows[[cell$id]]$gate[!gate_rows[[cell$id]]$pass]
    gate_failed <- length(failed_gates) > 0L
    gate_rows[[cell$id]]$advice <- ifelse(gate_rows[[cell$id]]$pass, "",
      mixedgp_simulation_diagnostic_advice())
    if (gate_failed) warning("Cell ", cell$id, " has flagged results: ",
      paste(failed_gates, collapse = "; "), ". ",
      mixedgp_simulation_diagnostic_advice(), call. = FALSE)
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
    if (isTRUE(config$fail_closed) && gate_failed) {
      stop(
        "Publication gates failed for ", cell$id, ": ",
        paste(failed_gates, collapse = ", "), "."
      )
    }
  }

  combined <- if ("aggregate" %in% config$stages) {
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
      "paired_validation", "gates"
    )],
    file.path(run_dir, "run_summary.rds")
  )
  message("Completed simulation run: ", normalizePath(run_dir))
  if ("aggregate" %in% config$stages) {
    report_env <- new.env(parent = environment())
    sys.source(file.path(config$code_dir, "reporting/experiment_reporting.R"), report_env)
    tryCatch(report_env$mixedgp_report_run(run_dir, "report", code_dir = config$code_dir),
      error = function(e) warning("Fitting completed, but reporting failed: ", conditionMessage(e),
        ". Use codes/cli/report_study.R to retry without fitting.", call. = FALSE))
  }
  invisible(out)
}

run_study1_simulation <- function(config = study1_simulation_config()) {
  if (!identical(config$study, "study1")) stop("Expected a Study I config.")
  mixedgp_run_simulation(config)
}

run_study2_simulation <- function(config = study2_simulation_config()) {
  if (!identical(config$study, "study2")) stop("Expected a Study II config.")
  mixedgp_run_simulation(config)
}

## Stable content identity shared by independent fitting/reporting layers.
mixedgp_dataset_identity <- function(data) {
  path <- tempfile(); on.exit(unlink(path), add=TRUE)
  saveRDS(data,path,version=3L,compress=FALSE)
  unname(tools::md5sum(path))
}
