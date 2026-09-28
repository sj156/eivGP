# Simulation gates: loaded by codes/simulation_helpers.R.

mixedgp_competitor_gate <- function(result, config) {
  if (!length(config$published_methods)) return(data.frame())
  status <- result$outputs$competitor_status
  expected <- result$cell$n_rep * length(config$published_methods)
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
    pass = n_attempted == expected && if (isTRUE(config$strict_competitors)) {
      failure_rate == 0
    } else {
      failure_rate <= 0.05
    },
    stringsAsFactors = FALSE
  )
}

mixedgp_mcmc_gate <- function(result, config) {
  diagnostics <- result$outputs$mcmc_diagnostics
  pass_column <- if (config$study == "study1") "gate_pass" else "mcmc_pass"
  expected <- result$cell$n_rep * length(result$cell$calibration_grid)
  values <- if (is.data.frame(diagnostics) && pass_column %in% names(diagnostics)) {
    as.logical(diagnostics[[pass_column]])
  } else {
    logical(0)
  }
  if (identical(config$study, "study2") && length(values) > 0L) {
    components <- c("target_functional_pass", "invariant_measurement_pass",
                    "raw_coordinate_pass", "n_calib")
    if (!all(components %in% names(diagnostics))) {
      values[] <- FALSE
    } else {
      values <- (values %in% TRUE) &
        (diagnostics$target_functional_pass %in% TRUE) &
        (diagnostics$invariant_measurement_pass %in% TRUE) &
        ((diagnostics$n_calib %in% 0L) | (diagnostics$raw_coordinate_pass %in% TRUE))
    }
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

