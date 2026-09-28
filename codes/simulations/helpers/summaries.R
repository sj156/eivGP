# Simulation summaries: loaded by codes/simulation_helpers.R.

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

mixedgp_study1_mechanism_contrasts <- function(comparisons, threshold_design = NULL) {
  if (is.null(threshold_design)) {
    return(mixedgp_bind_rows_base(lapply(c("balanced", "imbalanced"), function(design) {
      mixedgp_study1_mechanism_contrasts(comparisons, design)
    })))
  }
  if (!is.data.frame(comparisons) || nrow(comparisons) == 0L ||
      !"cell_id" %in% names(comparisons)) return(data.frame())
  control <- comparisons[
    comparisons$cell_id == paste0("eta0_", threshold_design) &
      comparisons$task %in% c(
        "response prediction Y*|x*,c*", "observable mean m(x,c)"
      ), , drop = FALSE
  ]
  active <- comparisons[
    comparisons$cell_id == paste0("eta1_", threshold_design) &
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
  out$contrast <- paste0("heterogeneity (", threshold_design, "): eta=1 minus eta=0")
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

