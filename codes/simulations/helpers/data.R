# Simulation data: loaded by codes/simulation_helpers.R.

mixedgp_simulation_diagnostic_advice <- function() {
  paste("Inspect diagnostic tables and chain traces before interpretation.",
        "If mixing is stable but ESS is low or MCSE is high, consider explicit",
        "continuation with continue_eivgp() for a compatible public fit, then",
        "recompute diagnostics. High R-hat or separated chains require checking",
        "identification, initialization and sampler behavior; more iterations",
        "alone may not help. Missing outputs/packages require repair or rerunning",
        "the affected method, not longer MCMC. Continuation is never automatic.")
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
    mixedgp_simulation_modules(), "simulation_helpers.R",
    paste0("simulations/helpers/", c("design", "data", "execution", "gates", "summaries", "workflow"), ".R"),
    "simulations/02_study1_monte_carlo.R", "simulations/02_study2_monte_carlo.R",
    "simulations/setup_study1_experiment.R", "simulations/setup_study2_experiment.R",
    "reporting/report_study1_results.R", "reporting/report_study2_results.R", "reporting/experiment_reporting.R",
    "simulations/run_study1_simulation.R", "simulations/run_study2_simulation.R"
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

mixedgp_write_run_metadata <- function(config,
                                       run_dir,
                                       preflight,
                                       runtime_preflight,
                                       task_plan) {
  metadata_dir <- file.path(run_dir, "config")
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
  mixedgp_atomic_save_rds(config, file.path(metadata_dir, "resolved_config.rds"))
  config_text <- capture.output(dput(config))
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
  capture.output(
    sessionInfo(), file = file.path(metadata_dir, "sessionInfo.txt")
  )
  invisible(metadata_dir)
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
    " dataset worker(s); chain workers per dataset: ",
    if (is.null(config$parallel$chain_workers)) 1L else config$parallel$chain_workers, "\n",
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
    for (design in c("balanced", "imbalanced")) {
    primary_ids <- intersect(
      paste0("eta", c(0, 1), "_", design), names(by_id)
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

