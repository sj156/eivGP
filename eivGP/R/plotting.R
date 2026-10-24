
# ---- reporting/report_study1_results.R ----
eivgp_study1_report <- function() {
## Reporting only: evaluated in a restored reporting context, never in the fitting driver.

mc_results$method <- factor(mc_results$method, levels = method_levels)
write.csv(
  mc_results,
  file.path(TAB_DIR, "study1_mc_raw_results_revised.csv"),
  row.names = FALSE
)
write.csv(
  competitor_status,
  file.path(TAB_DIR, "study1_competitor_status.csv"),
  row.names = FALSE
)
write.csv(
  ablation_results,
  file.path(TAB_DIR, "study1_ablation_predictive_raw.csv"),
  row.names = FALSE
)
write.csv(
  ablation_status,
  file.path(TAB_DIR, "study1_ablation_status.csv"),
  row.names = FALSE
)
write.csv(
  mean_recovery,
  file.path(TAB_DIR, "study1_mean_recovery_raw.csv"),
  row.names = FALSE
)
study1_write_csv_safe(
  latent_imputation,
  file.path(TAB_DIR, "study1_latent_imputation_raw.csv")
)
study1_write_csv_safe(
  surface_recovery,
  file.path(TAB_DIR, "study1_surface_recovery_raw.csv")
)
study1_write_csv_safe(
  ablation_surface_recovery,
  file.path(TAB_DIR, "study1_ablation_surface_recovery_raw.csv")
)

failure_summary <- competitor_status |>
  group_by(method) |>
  summarise(
    n_attempted = n(),
    n_success = sum(status == "success"),
    failure_rate = mean(status != "success"),
    .groups = "drop"
  )
write.csv(
  failure_summary,
  file.path(TAB_DIR, "study1_competitor_failure_summary.csv"),
  row.names = FALSE
)
if (any(failure_summary$failure_rate > 0.05)) {
  warning(
    "A Study I competitor failed or was unavailable in more than 5% of fits. ",
    "Do not report its performance without the operational warning.",
    call. = FALSE
  )
}

if (nrow(latent_imputation) > 0L) {
  latent_imputation_summary <- latent_imputation |>
    filter(score_status == "scored") |>
    pivot_longer(
      cols = all_of(c("Bias", "RMSE", "MAE", "Coverage95", "Width95")),
      names_to = "metric", values_to = "value"
    ) |>
    group_by(scenario, n_calib, method, metric) |>
    summarise(
      mean = mean(value, na.rm = TRUE),
      se = safe_se(value),
      n_rep_eff = sum(is.finite(value)),
      .groups = "drop"
    )
  study1_write_csv_safe(
    latent_imputation_summary,
    file.path(TAB_DIR, "study1_latent_imputation_summary.csv")
  )
}

surface_all <- bind_rows(surface_recovery, ablation_surface_recovery)
if (nrow(surface_all) > 0L) {
  surface_summary <- surface_all |>
    pivot_longer(
      cols = all_of(c("ISE", "Bias", "Coverage95", "Width95")),
      names_to = "metric", values_to = "value"
    ) |>
    group_by(scenario, n_calib, method, metric) |>
    summarise(
      mean = mean(value, na.rm = TRUE),
      se = safe_se(value),
      n_rep_eff = sum(is.finite(value)),
      .groups = "drop"
    )
  study1_write_csv_safe(
    surface_summary,
    file.path(TAB_DIR, "study1_surface_recovery_summary.csv")
  )
}

############################################################
## Summaries, paired differences, figure, and LaTeX tables
############################################################

mean_point_metrics <- c("RMSE", "MAE", "Bias")
mean_curve <- bind_rows(
  mean_recovery |> filter(!is.na(n_calib)),
  mean_recovery |>
    filter(is.na(n_calib)) |>
    select(-n_calib) |>
    tidyr::crossing(n_calib = calib_grid)
)
mean_summary <- mean_curve |>
  pivot_longer(
    cols = all_of(c(mean_point_metrics, "Coverage95", "Width95")),
    names_to = "metric",
    values_to = "value"
  ) |>
  group_by(scenario, n_calib, method, metric) |>
  summarise(
    mean = if (all(!is.finite(value))) NA_real_ else mean(value, na.rm = TRUE),
    se = safe_se(value),
    n_success = sum(is.finite(value)),
    .groups = "drop"
  )
write.csv(
  mean_summary,
  file.path(TAB_DIR, "study1_mean_recovery_summary.csv"),
  row.names = FALSE
)

p_mean <- mean_summary |>
  filter(metric %in% mean_point_metrics) |>
  ggplot(aes(x = n_calib, y = mean, color = method, group = method)) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
    width = 1.2, alpha = 0.7
  ) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_color_manual(values = method_cols, name = NULL, drop = FALSE) +
  labs(
    x = "Number of calibrated latent observations",
    y = "Monte Carlo mean",
    title = "Study I: recovery of the observed-input mean m(x,c)"
  ) +
  theme(legend.position = "bottom")
ggsave(
  file.path(FIG_DIR, "fig6_study1_mean_recovery.pdf"),
  p_mean, width = 11, height = 4.2
)

curve_results <- bind_rows(
  mc_results |> filter(!is.na(n_calib)),
  mc_results |>
    filter(is.na(n_calib)) |>
    select(-n_calib) |>
    tidyr::crossing(n_calib = calib_grid)
)
curve_long <- curve_results |>
  pivot_longer(
    cols = all_of(metric_levels),
    names_to = "metric",
    values_to = "value"
  ) |>
  mutate(metric = factor(metric, levels = metric_levels))
curve_summary <- curve_long |>
  group_by(scenario, n_calib, method, metric) |>
  summarise(mean = mean(value, na.rm = TRUE), se = safe_se(value), .groups = "drop")

p_mc <- ggplot(
  curve_summary,
  aes(x = n_calib, y = mean, color = method, group = method)
) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
    width = 1.2, alpha = 0.7
  ) +
  geom_hline(
    data = data.frame(
      metric = factor("Coverage95", levels = metric_levels), yint = 0.95
    ),
    aes(yintercept = yint), inherit.aes = FALSE,
    linetype = "dashed", color = "gray35"
  ) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  scale_color_manual(values = method_cols, name = NULL, drop = FALSE) +
  labs(
    x = "Number of calibrated latent observations",
    y = "Monte Carlo mean",
    title = "Study I: prediction of Y* given x* and c*"
  ) +
  theme(legend.position = "bottom")
ggsave(
  file.path(FIG_DIR, "fig5_study1_mc_metrics_revised.pdf"),
  p_mc, width = 12, height = 7.5
)

summary_long <- mc_results |>
  pivot_longer(
    cols = all_of(metric_levels),
    names_to = "metric",
    values_to = "value"
  ) |>
  group_by(n_calib, method, metric) |>
  summarise(
    mean = mean(value, na.rm = TRUE),
    se = safe_se(value),
    n_success = sum(is.finite(value)),
    .groups = "drop"
  )

summary_table <- summary_long |>
  mutate(value = mapply(format_mean_se, mean, se, USE.NAMES = FALSE)) |>
  select(n_calib, method, metric, value) |>
  pivot_wider(names_from = metric, values_from = value) |>
  arrange(is.na(n_calib), n_calib, method) |>
  mutate(
    calibration = ifelse(is.na(n_calib), "--", as.character(n_calib)),
    Method = as.character(method)
  ) |>
  select(
    `$|\\mathcal O|$` = calibration, Method,
    RMSE, CRPS, NLPD, Coverage95, Width95, IntervalScore95
  )
writeLines(
  knitr::kable(
    summary_table, format = "latex", booktabs = TRUE,
    align = "llcccccc", escape = FALSE
  ),
  file.path(TAB_DIR, "study1_mc_summary_revised.tex")
)

main_table <- summary_table |>
  filter(
    Method != "EIV-GP" |
      `$|\\mathcal O|$` %in% c("0", "10", "50")
  )
writeLines(
  knitr::kable(
    main_table, format = "latex", booktabs = TRUE,
    align = "llcccccc", escape = FALSE
  ),
  file.path(TAB_DIR, "study1_main_summary.tex")
)

fixed_competitors <- mc_results |>
  filter(method %in% STUDY1_PUBLISHED_COMPETITORS) |>
  select(rep, competitor = method, CRPS_comp = CRPS,
         IntervalScore95_comp = IntervalScore95)
paired_raw <- mc_results |>
  filter(method == "EIV-GP") |>
  select(rep, n_calib, CRPS_eiv = CRPS,
         IntervalScore95_eiv = IntervalScore95) |>
  inner_join(fixed_competitors, by = "rep") |>
  mutate(
    CRPS_difference = CRPS_eiv - CRPS_comp,
    IntervalScore95_difference =
      IntervalScore95_eiv - IntervalScore95_comp
  )
write.csv(
  paired_raw,
  file.path(TAB_DIR, "study1_paired_differences_raw.csv"),
  row.names = FALSE
)
paired_summary <- paired_raw |>
  group_by(n_calib, competitor) |>
  summarise(
    CRPS_difference = format_mean_se(
      mean(CRPS_difference), safe_se(CRPS_difference)
    ),
    IntervalScore95_difference = format_mean_se(
      mean(IntervalScore95_difference),
      safe_se(IntervalScore95_difference)
    ),
    n_pairs = n(),
    .groups = "drop"
  )
writeLines(
  knitr::kable(
    paired_summary, format = "latex", booktabs = TRUE, escape = FALSE,
    col.names = c(
      "$|\\mathcal O|$", "Competitor", "$\\Delta$ CRPS",
      "$\\Delta$ interval score", "Pairs"
    )
  ),
  file.path(TAB_DIR, "study1_paired_differences.tex")
)

if (nrow(ablation_results) > 0L) {
  ablation_summary <- ablation_results |>
    pivot_longer(
      cols = all_of(metric_levels),
      names_to = "metric", values_to = "value"
    ) |>
    group_by(n_calib, method, metric) |>
    summarise(
      mean = mean(value, na.rm = TRUE), se = safe_se(value),
      .groups = "drop"
    ) |>
    mutate(value = mapply(format_mean_se, mean, se, USE.NAMES = FALSE)) |>
    select(n_calib, method, metric, value) |>
    pivot_wider(names_from = metric, values_from = value) |>
    arrange(n_calib, method) |>
    select(
      `$|\\mathcal O|$` = n_calib,
      Method = method,
      RMSE, CRPS, Coverage95, Width95, IntervalScore95
    )
  writeLines(
    knitr::kable(
      ablation_summary, format = "latex", booktabs = TRUE,
      escape = FALSE
    ),
    file.path(TAB_DIR, "study1_ablation_summary.tex")
  )

  mechanism_table <- bind_rows(
    mc_results |>
      filter(method == "EIV-GP", n_calib == STUDY1_MECHANISM_CALIB),
    ablation_results |>
      filter(
        method %in% c("PI-GP", "CC-GP"),
        n_calib == STUDY1_MECHANISM_CALIB
      )
  ) |>
    pivot_longer(
      cols = all_of(metric_levels), names_to = "metric", values_to = "value"
    ) |>
    group_by(method, metric) |>
    summarise(
      mean = mean(value, na.rm = TRUE),
      se = safe_se(value),
      .groups = "drop"
    ) |>
    mutate(value = mapply(format_mean_se, mean, se, USE.NAMES = FALSE)) |>
    select(method, metric, value) |>
    pivot_wider(names_from = metric, values_from = value) |>
    mutate(
      Method = factor(method, levels = c("EIV-GP", "PI-GP", "CC-GP"))
    ) |>
    arrange(Method) |>
    select(Method, RMSE, CRPS, Coverage95, Width95, IntervalScore95)
  writeLines(
    knitr::kable(
      mechanism_table,
      format = "latex",
      booktabs = TRUE,
      escape = FALSE
    ),
    file.path(TAB_DIR, "study1_mechanism_summary.tex")
  )
}

design_manifest <- data.frame(
  design_tag = STUDY1_DESIGN_TAG,
  scenario = STUDY1_SCENARIO,
  heterogeneity_eta = STUDY1_HETEROGENEITY_ETA,
  threshold_design = STUDY1_THRESHOLD_DESIGN,
  n_train = n_train,
  n_test = n_test,
  n_rep = n_rep,
  m = m,
  sigma_eps = 0.10,
  calibration_grid = paste(calib_grid, collapse = ";"),
  mechanism_calibration = STUDY1_MECHANISM_CALIB,
  primary_target = "m(x,c)=E[f(x,U)|C=c]",
  latent_surface_target = "f(x,u)",
  predictive_target = "Y_star_given_X_star_C_star",
  latent_state_target = "U_given_X_C_Y_and_calibration",
  mean_evaluation_points = settings$n_m_eval,
  mean_posterior_draws = settings$n_m_draw,
  mean_latent_integration_draws = settings$n_m_latent,
  published_competitors = paste(STUDY1_PUBLISHED_COMPETITORS, collapse = ";"),
  appendix_ablations = if (isTRUE(STUDY1_RUN_ABLATIONS)) "PI-GP;CC-GP" else "",
  locked_eiv_results = isTRUE(STUDY1_REUSE_LOCKED_EIV),
  mcmc_gate_required = isTRUE(STUDY1_REQUIRE_MCMC_GATE),
  mcmc_max_rhat = STUDY1_MAX_RHAT,
  mcmc_min_ess = STUDY1_MIN_ESS,
  covariance_jitter = 0,
  forms_explicit_covariance_inverse = FALSE,
  sampler_control_manifest = "study1_sampler_control_manifest.csv",
  locked_eiv_file = if (isTRUE(STUDY1_REUSE_LOCKED_EIV)) {
    "tables/study1_mc_raw_results.csv"
  } else {
    ""
  },
  stringsAsFactors = FALSE
)
write.csv(
  design_manifest,
  file.path(TAB_DIR, "study1_design_manifest.csv"),
  row.names = FALSE
)
capture.output(
  sessionInfo(),
  file = file.path(RES_DIR, paste0("sessionInfo_", STUDY1_DESIGN_TAG, ".txt"))
)

message("Study I reporting outputs written under ", normalizePath(TAB_DIR))

}


# ---- reporting/report_study2_results.R ----
eivgp_study2_report <- function() {
## Reporting only: evaluated in a restored reporting context, never in the fitting driver.

write.csv(mc_results, file.path(TAB_DIR, paste0("study2_predictive_raw_", CACHE_TAG, ".csv")), row.names = FALSE)
write.csv(mc_mean_recovery, file.path(TAB_DIR, paste0("study2_mean_recovery_raw_", CACHE_TAG, ".csv")), row.names = FALSE)
write.csv(mc_mean_truth_rejection, file.path(TAB_DIR, paste0("study2_mean_truth_rejection_", CACHE_TAG, ".csv")), row.names = FALSE)
write.csv(mc_mean_truth_diagnostics, file.path(TAB_DIR, paste0("study2_mean_truth_diagnostics_", CACHE_TAG, ".csv")), row.names = FALSE)
write.csv(mc_diagnostics, file.path(TAB_DIR, paste0("study2_mcmc_diagnostics_", CACHE_TAG, ".csv")), row.names = FALSE)
write.csv(
  mc_target_diagnostics,
  file.path(
    TAB_DIR,
    paste0("study2_mcmc_target_diagnostics_", CACHE_TAG, ".csv")
  ),
  row.names = FALSE
)
write.csv(mc_imputation, file.path(TAB_DIR, paste0("study2_imputation_raw_", CACHE_TAG, ".csv")), row.names = FALSE)
write.csv(
  mc_imputation_status,
  file.path(
    TAB_DIR,
    paste0("study2_imputation_status_", CACHE_TAG, ".csv")
  ),
  row.names = FALSE
)
write.csv(mc_surface, file.path(TAB_DIR, paste0("study2_surface_raw_", CACHE_TAG, ".csv")), row.names = FALSE)
write_csv_safe(mc_ablation_results, file.path(TAB_DIR, paste0("study2_ablation_predictive_raw_", CACHE_TAG, ".csv")))
write_csv_safe(mc_ablation_surface, file.path(TAB_DIR, paste0("study2_ablation_surface_raw_", CACHE_TAG, ".csv")))
write_csv_safe(measurement_diagnostics, file.path(TAB_DIR, paste0("study2_measurement_diagnostics_", CACHE_TAG, ".csv")))
write_csv_safe(
  measurement_parameter_diagnostics,
  file.path(
    TAB_DIR,
    paste0("study2_measurement_parameter_diagnostics_", CACHE_TAG, ".csv")
  )
)
write_csv_safe(
  ablation_optimizer_attempts,
  file.path(
    TAB_DIR,
    paste0("study2_ablation_optimizer_attempts_", CACHE_TAG, ".csv")
  )
)
write_csv_safe(ablation_status, file.path(TAB_DIR, paste0("study2_ablation_status_", CACHE_TAG, ".csv")))
write.csv(competitor_status, file.path(TAB_DIR, paste0("study2_competitor_status_", CACHE_TAG, ".csv")), row.names = FALSE)
write.csv(pattern_counts, file.path(TAB_DIR, paste0("study2_pattern_counts_", CACHE_TAG, ".csv")), row.names = FALSE)
write_csv_safe(
  sampler_control_manifest,
  file.path(
    TAB_DIR,
    paste0("study2_sampler_control_manifest_", CACHE_TAG, ".csv")
  )
)

if (nrow(competitor_status) > 0L) {
  competitor_failure_summary <- competitor_status |>
    group_by(method) |>
    summarise(
      n_attempted = n(),
      n_success = sum(status == "success"),
      failure_rate = mean(status != "success"),
      .groups = "drop"
    )
  write.csv(
    competitor_failure_summary,
    file.path(
      TAB_DIR,
      paste0("study2_competitor_failure_summary_", CACHE_TAG, ".csv")
    ),
    row.names = FALSE
  )
  if (any(competitor_failure_summary$failure_rate > 0.05)) {
    warning(
      "At least one published Study II competitor failed or was unavailable ",
      "in more than 5% of attempted fits; inspect the saved failure summary ",
      "before reporting performance.",
      call. = FALSE
    )
  }
}

############################################################
## Design manifest
############################################################

design_manifest <- bind_rows(lapply(STUDY2_SCENARIOS, function(scenario) {
  pars <- make_study2_true_params(
    scenario = scenario, q = length(m_vec), m = m_vec[1L]
  )
  data.frame(
    design_tag = STUDY2_DESIGN_TAG,
    scenario = scenario,
    q = length(m_vec),
    d = d_latent,
    m = m_vec[1L],
    n_train = n_train,
    n_test = n_test,
    sigma_eps = pars$sigma_eps,
    lambda = pars$lambda,
    score_error = pars$score_error,
    response_interactions = pars$response_interactions,
    A = paste(formatC(as.vector(t(pars$A)), digits = 3, format = "f"), collapse = ";"),
    Omega = paste0("I", length(m_vec)),
    common_random_numbers = if (scenario == "logistic_misspec") {
      "shared X and U with primary"
    } else {
      "shared X, U, score innovations, and response innovations"
    },
    sampler_strategy = "collapsed_dictionary_ess",
    prospective_latent_sampler = switch(
      predictive_latent_sampler,
      minimax_tilting = paste(
        "minimax-tilted accept-reject",
        "(exact on successful completion)"
      ),
      rejection = "exact prior rejection",
      gibbs = "finite Gibbs diagnostic"
    ),
    diagnostic_gibbs_sweeps = diagnostic_n_new_latent_gibbs,
    rejection_max_batches = rejection_max_batches,
    primary_target = "m(x,c)=E[f(x,U)|C=c]",
    latent_surface_target = "f(x,u)",
    predictive_target = "Y_star_given_X_star_C_star",
    latent_state_target = "U_given_X_C_Y_and_calibration",
    mean_evaluation_points = n_m_eval,
    mean_posterior_draws = n_m_draw,
    mean_latent_integration_draws = n_m_latent,
    mean_truth_latent_draws = n_m_truth,
    covariance_jitter = 0,
    forms_explicit_covariance_inverse = FALSE,
    sampler_control_manifest = paste0(
      "study2_sampler_control_manifest_", CACHE_TAG, ".csv"
    ),
    mcmc_gate_rule = paste0(
      "All calibrations: Rao-Blackwellized response moments and every ordinal ",
      "score correlation/standardized cutpoint; calibrated fits additionally ",
      "require free-coordinate, f(x,u), and prospective-U diagnostics."
    ),
    mcmc_rhat_limit = STUDY2_MCMC_RHAT_LIMIT,
    mcmc_raw_ess_limit = STUDY2_MCMC_RAW_ESS_LIMIT,
    mcmc_target_bulk_ess_limit = STUDY2_MCMC_TARGET_BULK_ESS_LIMIT,
    mcmc_target_tail_ess_limit = STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT,
    mcmc_target_n_points = STUDY2_MCMC_TARGET_N_POINTS,
    mcmc_target_draws_per_chain = STUDY2_MCMC_TARGET_DRAWS_PER_CHAIN,
    mcmc_target_integration_draws = STUDY2_MCMC_TARGET_N_LATENT,
    measurement_n_iter = measurement_n_iter,
    measurement_burn = measurement_burn,
    measurement_thin = measurement_thin,
    measurement_n_chains = measurement_n_chains,
    measurement_rhat_limit = STUDY2_MEAS_RHAT_LIMIT,
    measurement_bulk_ess_limit = STUDY2_MEAS_BULK_ESS_LIMIT,
    measurement_tail_ess_limit = STUDY2_MEAS_TAIL_ESS_LIMIT,
    ablation_gp_n_starts = STUDY2_ABLATION_GP_N_STARTS,
    ablation_gp_maxit = STUDY2_ABLATION_GP_MAXIT,
    tau = paste(formatC(as.vector(t(pars$tau)), digits = 6, format = "f"), collapse = ";"),
    calibration_grid = paste(scenario_calib_grid(scenario), collapse = ";"),
    published_competitors = paste(STUDY2_PUBLISHED_COMPETITORS, collapse = ";"),
    appendix_ablations = if (
      isTRUE(STUDY2_RUN_ABLATIONS) && scenario %in% STUDY2_ABLATION_SCENARIOS
    ) {
      "PI-GP;CC-GP;Full-U GP"
    } else {
      ""
    },
    stringsAsFactors = FALSE
  )
}))
write.csv(design_manifest, file.path(TAB_DIR, paste0("study2_design_manifest_", CACHE_TAG, ".csv")), row.names = FALSE)

capture.output(
  sessionInfo(),
  file = file.path(RES_DIR, paste0("sessionInfo_", CACHE_TAG, ".txt"))
)

############################################################
## Publication summaries
############################################################

mean_calibration_grid <- bind_rows(lapply(STUDY2_SCENARIOS, function(sc) {
  data.frame(scenario = sc, n_calib = scenario_calib_grid(sc))
}))
mean_curve <- bind_rows(
  mc_mean_recovery |> filter(!is.na(n_calib)),
  mc_mean_recovery |>
    filter(is.na(n_calib)) |>
    select(-n_calib) |>
    inner_join(mean_calibration_grid, by = "scenario")
)
mean_summary <- mean_curve |>
  pivot_longer(
    cols = all_of(c("RMSE", "MAE", "Bias", "Coverage95", "Width95")),
    names_to = "metric",
    values_to = "value"
  ) |>
  group_by(scenario, n_calib, method, metric) |>
  summarise(
    mean = if (all(!is.finite(value))) NA_real_ else mean(value, na.rm = TRUE),
    se = safe_se(value),
    n_success = sum(is.finite(value)),
    .groups = "drop"
  ) |>
  mutate(scenario_label = vapply(scenario, scenario_label, character(1)))
write.csv(
  mean_summary,
  file.path(TAB_DIR, paste0("study2_mean_recovery_summary_", CACHE_TAG, ".csv")),
  row.names = FALSE
)

p_mean <- mean_summary |>
  filter(metric %in% c("RMSE", "MAE", "Bias")) |>
  ggplot(aes(n_calib, mean, color = method, group = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  geom_errorbar(
    aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
    width = 1.2, alpha = 0.7
  ) +
  facet_grid(scenario_label ~ metric, scales = "free_y") +
  scale_color_manual(values = method_cols, name = NULL, drop = FALSE) +
  labs(
    x = "Number of calibrated latent observations",
    y = "Monte Carlo mean",
    title = "Study II: recovery of the observed-input mean m(x,c)"
  ) +
  theme(legend.position = "bottom")
save_plot(
  file.path(FIG_DIR, paste0("study2_mean_recovery_", CACHE_TAG)),
  p_mean,
  width = 11,
  height = 9
)

mc_overall <- mc_results |>
  filter(evaluation_stratum == "overall")

metric_summary <- mc_overall |>
  pivot_longer(
    cols = all_of(metric_levels),
    names_to = "metric",
    values_to = "value"
  ) |>
  group_by(scenario, n_calib, method, metric) |>
  summarise(
    mean = mean(value, na.rm = TRUE),
    se = safe_se(value),
    n_rep_eff = sum(is.finite(value)),
    .groups = "drop"
  )
write.csv(metric_summary, file.path(TAB_DIR, paste0("study2_metric_summary_", CACHE_TAG, ".csv")), row.names = FALSE)

pattern_summary <- mc_results |>
  filter(evaluation_stratum != "overall") |>
  pivot_longer(
    cols = all_of(metric_levels),
    names_to = "metric",
    values_to = "value"
  ) |>
  group_by(scenario, evaluation_stratum, n_calib, method, metric) |>
  summarise(
    mean = mean(value, na.rm = TRUE),
    se = safe_se(value),
    n_rep_eff = sum(is.finite(value)),
    .groups = "drop"
  )
write.csv(pattern_summary, file.path(TAB_DIR, paste0("study2_pattern_summary_", CACHE_TAG, ".csv")), row.names = FALSE)

if (nrow(mc_ablation_results) > 0L) {
  ablation_metric_summary <- mc_ablation_results |>
    filter(evaluation_stratum == "overall") |>
    pivot_longer(
      cols = all_of(metric_levels),
      names_to = "metric",
      values_to = "value"
    ) |>
    group_by(scenario, n_calib, method, metric) |>
    summarise(
      mean = mean(value, na.rm = TRUE),
      se = safe_se(value),
      n_rep_eff = sum(is.finite(value)),
      .groups = "drop"
    )
  write.csv(
    ablation_metric_summary,
    file.path(
      TAB_DIR,
      paste0("study2_ablation_metric_summary_", CACHE_TAG, ".csv")
    ),
    row.names = FALSE
  )
} else {
  ablation_metric_summary <- data.frame()
}

if (nrow(mc_imputation) > 0L) {
  imputation_summary <- mc_imputation |>
    pivot_longer(
      cols = all_of(c("Bias", "RMSE", "MAE", "Coverage95", "Width95")),
      names_to = "metric",
      values_to = "value"
    ) |>
    group_by(scenario, target, n_calib, method, coordinate, metric) |>
    summarise(
      mean = mean(value, na.rm = TRUE),
      se = safe_se(value),
      n_rep_eff = sum(is.finite(value)),
      .groups = "drop"
    )
  write.csv(
    imputation_summary,
    file.path(
      TAB_DIR,
      paste0("study2_imputation_summary_", CACHE_TAG, ".csv")
    ),
    row.names = FALSE
  )
} else {
  imputation_summary <- data.frame()
}

surface_all <- bind_rows(mc_surface, mc_ablation_surface)
if (nrow(surface_all) > 0L) {
  surface_summary <- surface_all |>
    pivot_longer(
      cols = c(ISE, Bias, Coverage95, Width95),
      names_to = "metric",
      values_to = "value"
    ) |>
    group_by(scenario, n_calib, method, metric) |>
    summarise(
      mean = mean(value, na.rm = TRUE),
      se = safe_se(value),
      n_rep_eff = sum(is.finite(value)),
      .groups = "drop"
    )
  write.csv(
    surface_summary,
    file.path(
      TAB_DIR,
      paste0("study2_surface_recovery_summary_", CACHE_TAG, ".csv")
    ),
    row.names = FALSE
  )
} else {
  surface_summary <- data.frame()
}

if (nrow(mc_ablation_results) > 0L) {
  ablation_table <- mc_ablation_results |>
    filter(scenario == "primary", evaluation_stratum == "overall") |>
    group_by(n_calib, method) |>
    summarise(
      RMSE = format_mean_se(mean(RMSE), safe_se(RMSE)),
      CRPS = format_mean_se(mean(CRPS), safe_se(CRPS)),
      Coverage95 = format_mean_se(mean(Coverage95), safe_se(Coverage95)),
      Width95 = format_mean_se(mean(Width95), safe_se(Width95)),
      IntervalScore95 = format_mean_se(
        mean(IntervalScore95), safe_se(IntervalScore95)
      ),
      .groups = "drop"
    ) |>
    arrange(factor(method, levels = c("PI-GP", "CC-GP")), n_calib) |>
    select(Calibration = n_calib, Method = method, everything())

  writeLines(
    knitr::kable(
      ablation_table,
      format = "latex",
      booktabs = TRUE,
      escape = FALSE,
      col.names = c(
        "$|\\mathcal O|$", "Ablation", "RMSE", "CRPS",
        "Coverage95", "Width95", "IntervalScore95"
      )
    ),
    file.path(
      TAB_DIR,
      paste0("study2_ablation_summary_", CACHE_TAG, ".tex")
    )
  )
}

if (nrow(surface_all) > 0L) {
  surface_table <- surface_all |>
    filter(scenario == "primary") |>
    group_by(n_calib, method) |>
    summarise(
      ISE = format_mean_se(mean(ISE), safe_se(ISE)),
      Bias = format_mean_se(mean(Bias), safe_se(Bias)),
      Coverage95 = format_mean_se(
        mean(Coverage95), safe_se(Coverage95)
      ),
      Width95 = format_mean_se(mean(Width95), safe_se(Width95)),
      .groups = "drop"
    ) |>
    mutate(Calibration = ifelse(is.na(n_calib), "--", n_calib)) |>
    arrange(
      factor(method, levels = c("EIV-GP", "PI-GP", "CC-GP", "Full-U GP")),
      n_calib
    ) |>
    select(Calibration, Method = method, ISE, Bias, Coverage95, Width95)

  writeLines(
    knitr::kable(
      surface_table,
      format = "latex",
      booktabs = TRUE,
      escape = FALSE,
      col.names = c(
        "$|\\mathcal O|$", "Method", "ISE", "Bias",
        "Coverage95", "Width95"
      )
    ),
    file.path(
      TAB_DIR,
      paste0("study2_surface_recovery_mc_summary_", CACHE_TAG, ".tex")
    )
  )
}

paired_differences <- mc_overall |>
  select(rep, scenario, n_calib, method, CRPS, IntervalScore95) |>
  filter(method == "EIV-GP") |>
  rename(EIV_CRPS = CRPS, EIV_IntervalScore95 = IntervalScore95) |>
  inner_join(
    mc_overall |>
      filter(!method %in% c("EIV-GP", "Oracle")) |>
      select(rep, scenario, competitor = method, Comp_CRPS = CRPS,
             Comp_IntervalScore95 = IntervalScore95),
    by = c("rep", "scenario")
  ) |>
  mutate(
    CRPS_difference = EIV_CRPS - Comp_CRPS,
    IntervalScore95_difference =
      EIV_IntervalScore95 - Comp_IntervalScore95
  )
write.csv(paired_differences, file.path(TAB_DIR, paste0("study2_paired_differences_", CACHE_TAG, ".csv")), row.names = FALSE)

paired_summary <- paired_differences |>
  group_by(scenario, n_calib, competitor) |>
  summarise(
    CRPS_difference = mean(CRPS_difference),
    CRPS_difference_se = safe_se(CRPS_difference),
    IntervalScore95_difference = mean(IntervalScore95_difference),
    IntervalScore95_difference_se = safe_se(IntervalScore95_difference),
    .groups = "drop"
  )
write.csv(paired_summary, file.path(TAB_DIR, paste0("study2_paired_summary_", CACHE_TAG, ".csv")), row.names = FALSE)

############################################################
## Main primary-setting table: every method is its own row
############################################################

primary_summary <- mc_overall |>
  filter(scenario == "primary") |>
  group_by(n_calib, method) |>
  summarise(
    RMSE_mean = mean(RMSE), RMSE_se = safe_se(RMSE),
    CRPS_mean = mean(CRPS), CRPS_se = safe_se(CRPS),
    NLPD_mean = mean(NLPD), NLPD_se = safe_se(NLPD),
    Coverage_mean = mean(Coverage95), Coverage_se = safe_se(Coverage95),
    Width_mean = mean(Width95), Width_se = safe_se(Width95),
    Score_mean = mean(IntervalScore95), Score_se = safe_se(IntervalScore95),
    .groups = "drop"
  ) |>
  mutate(
    Calibration = ifelse(is.na(n_calib), "--", as.character(n_calib)),
    RMSE = mapply(format_mean_se, RMSE_mean, RMSE_se),
    CRPS = mapply(format_mean_se, CRPS_mean, CRPS_se),
    NLPD = mapply(format_mean_se, NLPD_mean, NLPD_se),
    Coverage95 = mapply(format_mean_se, Coverage_mean, Coverage_se),
    Width95 = mapply(format_mean_se, Width_mean, Width_se),
    IntervalScore95 = mapply(format_mean_se, Score_mean, Score_se)
  ) |>
  arrange(factor(method, levels = method_levels), n_calib) |>
  select(
    Calibration, Method = method, RMSE, CRPS, NLPD,
    Coverage95, Width95, IntervalScore95
  )

write.csv(primary_summary, file.path(TAB_DIR, paste0("study2_main_summary_", CACHE_TAG, ".csv")), row.names = FALSE)
writeLines(
  knitr::kable(
    primary_summary,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    align = "llcccccc",
    col.names = c(
      "$|\\mathcal O|$", "Method", "RMSE", "CRPS", "NLPD",
      "Coverage95", "Width95", "IntervalScore95"
    )
  ),
  file.path(TAB_DIR, paste0("study2_main_summary_", CACHE_TAG, ".tex"))
)

############################################################
## Prespecified design contrasts and sparse-pattern summaries
############################################################

contrast_eiv <- mc_overall |>
  filter(
    scenario %in% c(
      "primary", "latent_additive_control", "high_uncertainty"
    ),
    method == "EIV-GP",
    n_calib == STUDY2_CONTRAST_CALIB
  ) |>
  select(
    rep, scenario, EIV_CRPS = CRPS,
    EIV_IntervalScore95 = IntervalScore95
  )
contrast_competitors <- mc_overall |>
  filter(
    scenario %in% c(
      "primary", "latent_additive_control", "high_uncertainty"
    ),
    method %in% STUDY2_PUBLISHED_COMPETITORS,
    is.na(n_calib)
  ) |>
  select(
    rep, scenario, competitor = method, Comp_CRPS = CRPS,
    Comp_IntervalScore95 = IntervalScore95
  )

scenario_advantages <- contrast_eiv |>
  inner_join(contrast_competitors, by = c("rep", "scenario")) |>
  mutate(
    CRPS_advantage = Comp_CRPS - EIV_CRPS,
    IntervalScore95_advantage =
      Comp_IntervalScore95 - EIV_IntervalScore95
  )
write.csv(
  scenario_advantages,
  file.path(
    TAB_DIR,
    paste0("study2_scenario_advantages_", CACHE_TAG, ".csv")
  ),
  row.names = FALSE
)

advantage_wide <- scenario_advantages |>
  select(
    rep, competitor, scenario, CRPS_advantage,
    IntervalScore95_advantage
  ) |>
  pivot_wider(
    names_from = scenario,
    values_from = c(CRPS_advantage, IntervalScore95_advantage)
  )

required_advantage_columns <- c(
  "CRPS_advantage_primary",
  "CRPS_advantage_latent_additive_control",
  "CRPS_advantage_high_uncertainty",
  "IntervalScore95_advantage_primary",
  "IntervalScore95_advantage_latent_additive_control",
  "IntervalScore95_advantage_high_uncertainty"
)
if (nrow(advantage_wide) > 0L &&
    all(required_advantage_columns %in% names(advantage_wide))) {
  advantage_contrasts <- bind_rows(
    advantage_wide |>
      transmute(
        rep, competitor,
        contrast = "Interactions: primary minus additive",
        CRPS_advantage_change =
          CRPS_advantage_primary - CRPS_advantage_latent_additive_control,
        IntervalScore95_advantage_change =
          IntervalScore95_advantage_primary -
            IntervalScore95_advantage_latent_additive_control
      ),
    advantage_wide |>
      transmute(
        rep, competitor,
        contrast = "Uncertainty: high minus primary",
        CRPS_advantage_change =
          CRPS_advantage_high_uncertainty - CRPS_advantage_primary,
        IntervalScore95_advantage_change =
          IntervalScore95_advantage_high_uncertainty -
            IntervalScore95_advantage_primary
      )
  )
  write.csv(
    advantage_contrasts,
    file.path(
      TAB_DIR,
      paste0("study2_advantage_contrasts_raw_", CACHE_TAG, ".csv")
    ),
    row.names = FALSE
  )

  main_contrast_table <- advantage_contrasts |>
    group_by(contrast, competitor) |>
    summarise(
      `Change in CRPS advantage` = format_mean_se(
        mean(CRPS_advantage_change), safe_se(CRPS_advantage_change)
      ),
      `Change in interval-score advantage` = format_mean_se(
        mean(IntervalScore95_advantage_change),
        safe_se(IntervalScore95_advantage_change)
      ),
      Pairs = sum(is.finite(CRPS_advantage_change)),
      .groups = "drop"
    ) |>
    mutate(
      Contrast = factor(
        contrast,
        levels = c(
          "Interactions: primary minus additive",
          "Uncertainty: high minus primary"
        )
      ),
      Competitor = factor(competitor, levels = STUDY2_PUBLISHED_COMPETITORS)
    ) |>
    arrange(Contrast, Competitor) |>
    select(
      Contrast, Competitor, `Change in CRPS advantage`,
      `Change in interval-score advantage`, Pairs
    )
  write.csv(
    main_contrast_table,
    file.path(
      TAB_DIR,
      paste0("study2_design_advantage_contrasts_", CACHE_TAG, ".csv")
    ),
    row.names = FALSE
  )
  writeLines(
    knitr::kable(
      main_contrast_table,
      format = "latex",
      booktabs = TRUE,
      escape = FALSE
    ),
    file.path(TAB_DIR, "study2_design_advantage_contrasts.tex")
  )
}

logistic_rows <- mc_overall |>
  filter(
    scenario == "logistic_misspec",
    (method == "EIV-GP" & n_calib == STUDY2_CONTRAST_CALIB) |
      (method != "EIV-GP" & is.na(n_calib))
  )
if (nrow(logistic_rows) > 0L) {
  logistic_table <- logistic_rows |>
    group_by(method) |>
    summarise(
      n_rep_eff = sum(is.finite(CRPS)),
      RMSE = format_mean_se(mean(RMSE), safe_se(RMSE)),
      CRPS = format_mean_se(mean(CRPS), safe_se(CRPS)),
      Coverage95 = format_mean_se(
        mean(Coverage95), safe_se(Coverage95)
      ),
      Width95 = format_mean_se(mean(Width95), safe_se(Width95)),
      IntervalScore95 = format_mean_se(
        mean(IntervalScore95), safe_se(IntervalScore95)
      ),
      .groups = "drop"
    ) |>
    mutate(Method = factor(method, levels = method_levels)) |>
    arrange(Method) |>
    select(
      Method, RMSE, CRPS, Coverage95, Width95, IntervalScore95,
      `Successful replications` = n_rep_eff
    )
  write.csv(
    logistic_table,
    file.path(
      TAB_DIR,
      paste0("study2_logistic_robustness_summary_", CACHE_TAG, ".csv")
    ),
    row.names = FALSE
  )
  writeLines(
    knitr::kable(
      logistic_table,
      format = "latex",
      booktabs = TRUE,
      escape = FALSE
    ),
    file.path(TAB_DIR, "study2_logistic_robustness_summary.tex")
  )
}

pattern_rows <- mc_results |>
  filter(
    scenario == "primary",
    evaluation_stratum %in% c("common", "rare", "unobserved"),
    (method == "EIV-GP" & n_calib == STUDY2_CONTRAST_CALIB) |
      (method != "EIV-GP" & is.na(n_calib))
  )
if (nrow(pattern_rows) > 0L) {
  pattern_sizes <- pattern_counts |>
    filter(scenario == "primary") |>
    group_by(evaluation_stratum) |>
    summarise(mean_n_test = mean(n_test), .groups = "drop")
  pattern_table <- pattern_rows |>
    group_by(evaluation_stratum, method) |>
    summarise(
      n_rep_eff = sum(is.finite(CRPS)),
      CRPS = format_mean_se(mean(CRPS), safe_se(CRPS)),
      Coverage95 = format_mean_se(
        mean(Coverage95), safe_se(Coverage95)
      ),
      Width95 = format_mean_se(mean(Width95), safe_se(Width95)),
      .groups = "drop"
    ) |>
    left_join(pattern_sizes, by = "evaluation_stratum") |>
    mutate(
      Pattern = factor(
        evaluation_stratum,
        levels = c("common", "rare", "unobserved")
      ),
      Method = factor(method, levels = method_levels),
      `Mean test-set size` = sprintf("%.1f", mean_n_test)
    ) |>
    arrange(Pattern, Method) |>
    select(
      Pattern, Method, `Mean test-set size`, CRPS, Coverage95, Width95,
      `Successful replications` = n_rep_eff
    )
  write.csv(
    pattern_table,
    file.path(
      TAB_DIR,
      paste0("study2_pattern_strata_summary_", CACHE_TAG, ".csv")
    ),
    row.names = FALSE
  )
  writeLines(
    knitr::kable(
      pattern_table,
      format = "latex",
      booktabs = TRUE,
      escape = FALSE
    ),
    file.path(TAB_DIR, "study2_pattern_strata_summary.tex")
  )

  pattern_advantages <- pattern_rows |>
    filter(method == "EIV-GP") |>
    select(
      rep, evaluation_stratum, EIV_CRPS = CRPS,
      EIV_IntervalScore95 = IntervalScore95
    ) |>
    inner_join(
      pattern_rows |>
        filter(method %in% STUDY2_PUBLISHED_COMPETITORS) |>
        select(
          rep, evaluation_stratum, competitor = method,
          Comp_CRPS = CRPS,
          Comp_IntervalScore95 = IntervalScore95
        ),
      by = c("rep", "evaluation_stratum")
    ) |>
    mutate(
      CRPS_advantage = Comp_CRPS - EIV_CRPS,
      IntervalScore95_advantage =
        Comp_IntervalScore95 - EIV_IntervalScore95
    )
  write.csv(
    pattern_advantages,
    file.path(
      TAB_DIR,
      paste0("study2_pattern_advantages_raw_", CACHE_TAG, ".csv")
    ),
    row.names = FALSE
  )
  pattern_advantage_table <- pattern_advantages |>
    group_by(evaluation_stratum, competitor) |>
    summarise(
      `EIV CRPS advantage` = format_mean_se(
        mean(CRPS_advantage), safe_se(CRPS_advantage)
      ),
      `EIV interval-score advantage` = format_mean_se(
        mean(IntervalScore95_advantage),
        safe_se(IntervalScore95_advantage)
      ),
      Pairs = sum(is.finite(CRPS_advantage)),
      .groups = "drop"
    ) |>
    mutate(
      Pattern = factor(
        evaluation_stratum,
        levels = c("common", "rare", "unobserved")
      ),
      Competitor = factor(competitor, levels = STUDY2_PUBLISHED_COMPETITORS)
    ) |>
    arrange(Pattern, Competitor) |>
    select(
      Pattern, Competitor, `EIV CRPS advantage`,
      `EIV interval-score advantage`, Pairs
    )
  writeLines(
    knitr::kable(
      pattern_advantage_table,
      format = "latex",
      booktabs = TRUE,
      escape = FALSE
    ),
    file.path(TAB_DIR, "study2_pattern_advantage_summary.tex")
  )
}

############################################################
## Figures: expand fixed competitors only for plotting
############################################################

plot_summary <- metric_summary |>
  filter(scenario == "primary")

fixed_rows <- plot_summary |> filter(is.na(n_calib))
varying_rows <- plot_summary |> filter(!is.na(n_calib))
if (nrow(fixed_rows) > 0L) {
  fixed_rows <- merge(
    fixed_rows |> select(-n_calib),
    data.frame(n_calib = scenario_calib_grid("primary")),
    by = NULL
  )
}
plot_summary <- bind_rows(varying_rows, fixed_rows)
plot_summary$metric <- factor(plot_summary$metric, levels = metric_levels)
plot_summary$method <- factor(plot_summary$method, levels = method_levels)

p_primary <- ggplot(
  plot_summary,
  aes(x = n_calib, y = mean, color = method, group = method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  geom_errorbar(
    aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
    width = 1.2,
    alpha = 0.65
  ) +
  geom_hline(
    data = data.frame(
      metric = factor("Coverage95", levels = metric_levels),
      yint = 0.95
    ),
    aes(yintercept = yint),
    inherit.aes = FALSE,
    linetype = "dashed",
    color = "gray35"
  ) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  scale_color_manual(values = method_cols, drop = FALSE, name = NULL) +
  labs(
    x = "Number of calibrated latent observations",
    y = "Monte Carlo mean",
    title = "Study II: primary multivariate interactive setting"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

save_plot(
  file.path(FIG_DIR, paste0("fig5_study2_mc_metrics_", CACHE_TAG)),
  p_primary,
  width = 12,
  height = 7.5
)

contrast_crps <- metric_summary |>
  filter(
    scenario %in% c("primary", "logistic_misspec"),
    metric == "CRPS",
    (method == "EIV-GP" & n_calib == STUDY2_CONTRAST_CALIB) |
      method != "EIV-GP"
  ) |>
  mutate(
    scenario = factor(
      scenario,
      levels = c("primary", "logistic_misspec"),
      labels = vapply(
        c("primary", "logistic_misspec"),
        scenario_label,
        character(1)
      )
    ),
    method = factor(method, levels = method_levels)
  )

p_contrast <- ggplot(
  contrast_crps,
  aes(x = method, y = mean, color = method)
) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), width = 0.15) +
  facet_wrap(~scenario, scales = "free_y") +
  scale_color_manual(values = method_cols, drop = FALSE, guide = "none") +
  labs(x = NULL, y = "CRPS", title = "Study II: prespecified design contrasts") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

save_plot(
  file.path(FIG_DIR, paste0("fig_study2_scenario_crps_", CACHE_TAG)),
  p_contrast,
  width = 11,
  height = 6.5
)

############################################################
## Publication gates
############################################################

if (nrow(mc_diagnostics) > 0L) {
  mcmc_gate <- mc_diagnostics |>
    group_by(mcmc_gate_rule) |>
    summarise(
      n_fits = n(),
      n_pass = sum(mcmc_pass),
      pass_rate = mean(mcmc_pass),
      rhat_limit = first(mcmc_rhat_limit),
      raw_ess_limit = first(mcmc_raw_ess_limit),
      target_bulk_ess_limit = first(mcmc_target_bulk_ess_limit),
      target_tail_ess_limit = first(mcmc_target_tail_ess_limit),
      .groups = "drop"
    )
  write.csv(mcmc_gate, file.path(TAB_DIR, paste0("study2_mcmc_gate_", CACHE_TAG, ".csv")), row.names = FALSE)

  if (!all(mc_diagnostics$mcmc_pass)) {
    warning("MCMC results include diagnostic warnings; finite estimates are retained. ",
            mixedgp_simulation_diagnostic_advice(), call. = FALSE)
  }
}

cat("\nStudy II reporting completed for saved cache tag:\n", CACHE_TAG, "\n")
cat("Primary table:\n", file.path(TAB_DIR, paste0("study2_main_summary_", CACHE_TAG, ".tex")), "\n")

}


# ---- reporting/experiment_reporting.R ----
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
      local({
        previous_dir <- setwd(code_dir)
        on.exit(setwd(previous_dir), add = TRUE)
        eval(body(get(paste0("eivgp_", config$study, "_setup"))), envir = env)
      })
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
      eval(body(get(paste0("eivgp_", config$study, "_report"))), envir = env)
      if (length(plot_status)) utils::write.csv(do.call(rbind, plot_status),
        file.path(target, "plot_status.csv"), row.names = FALSE)
      utils::write.csv(data.frame(source_run = run_dir, source_cell = cell$id,
        source_input = if (file.exists(checkpoint)) checkpoint else legacy,
        source_md5 = unname(tools::md5sum(if (file.exists(checkpoint)) checkpoint else legacy)),
        action = action), file.path(target, paste0(action, "_provenance.csv")), row.names = FALSE)
      ids <- eivgp_module_identity()
      utils::write.csv(data.frame(file = names(ids), md5 = unname(ids)),
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


plot_comparison <- function(summary, metric="RMSE") {
  if(!requireNamespace("ggplot2",quietly=TRUE)) stop("Install ggplot2.")
  z<-summary[summary$metric==metric,,drop=FALSE]
  ggplot2::ggplot(z,ggplot2::aes(x=factor(n_calib),y=mean,color=method,group=method))+
    ggplot2::geom_point()+ggplot2::geom_line()+ggplot2::facet_wrap(~cell_id)+
    ggplot2::scale_x_discrete(labels=function(x) ifelse(is.na(x),"Not used",x))+
    ggplot2::theme_minimal()+ggplot2::labs(x="Calibration size",y=metric,color="Method",
      caption="Not used: the method does not use calibration measurements.")
}
report_table <- function(x, caption=NULL, digits=3L) {
  if(!requireNamespace("knitr",quietly=TRUE)) stop("Install knitr.")
  knitr::kable(x,caption=caption,digits=digits)
}

summarize_diagnostics <- function(diagnostics) {
  d<-diagnostics
  if(!nrow(d))return(data.frame())
  group<-intersect(c("cell_id","n_calib"),names(d))
  keys<-do.call(paste,c(d[group],sep="|"))
  rows<-lapply(unique(keys),function(key){
    z<-d[keys==key,,drop=FALSE];out<-z[1,group,drop=FALSE];out$fits<-nrow(z)
    for(n in intersect(c("target_functional_pass","target_gate_pass","mcmc_pass","gate_pass"),names(z)))
      out[[paste0(n,"_count")]]<-sum(z[[n]] %in% c(TRUE,"TRUE"))
    for(n in intersect(c("predictive_min_bulk_ess","predictive_min_tail_ess","predictive_max_rhat"),names(z))) {
      v<-as.numeric(z[[n]]);v<-v[is.finite(v)]
      out[[n]]<-if(!length(v))NA_real_ else if(grepl("max",n))max(v) else min(v)
    }
    out
  })
  do.call(rbind,rows)
}

mixedgp_combine_predictions <- function(own, competitors, status, identity, cells) {
  metrics <- c("RMSE","MAE","CRPS","NLPD","Coverage95","Width95","IntervalScore95")
  key <- function(x) paste(x$cell_id,x$rep,sep="/")
  if(anyDuplicated(key(identity))) stop("Duplicate MCMC dataset identities.")
  expected <- mixedgp_bind_rows_base(lapply(cells,function(c)
    data.frame(cell_id=c$id,rep=seq_len(c$n_rep))))
  for(name in c("own","competitors","status")) {
    x <- get(name)
    if(nrow(x)) x <- x[key(x) %in% key(expected),,drop=FALSE]
    assign(name,x)
  }
  if(nrow(own) && any(!key(own) %in% key(identity)))
    stop("MCMC metrics lack dataset identity; do not infer identity from replication number alone.")
  if(nrow(competitors)) {
    required <- c("cell_id","rep","method","dataset_md5","status")
    if(!all(required %in% names(status))) stop("Competitor statuses lack dataset checksums.")
    sk <- paste(key(status),status$method)
    if(anyDuplicated(sk)) stop("Duplicate competitor statuses.")
    at <- match(paste(key(competitors),competitors$method),sk)
    if(anyNA(at) || any(status$status[at]!="success")) stop("Competitor metrics lack successful status.")
    ik <- match(key(competitors),key(identity))
    shared <- !is.na(ik)
    if(anyNA(status$dataset_md5[at]) || any(!nzchar(status$dataset_md5[at]))) stop("Missing competitor checksum.")
    if(any(status$dataset_md5[at[shared]] != identity$dataset_md5[ik[shared]]))
      stop("Dataset checksum mismatch: results cannot be paired.")
  }
  normalize <- function(x,source) {
    if(!nrow(x)) return(data.frame())
    if(!"evaluation_stratum" %in% names(x)) x$evaluation_stratum <- "overall"
    x <- x[x$evaluation_stratum=="overall",,drop=FALSE]
    x$source <- source
    if(anyDuplicated(paste(key(x),x$method,x$n_calib))) stop("Duplicate method/dataset/calibration metrics.")
    x
  }
  own <- normalize(own,"mcmc_experiment")
  competitors <- normalize(competitors,"competitor_experiment")
  combined <- mixedgp_bind_rows_base(list(own,competitors))
  summary <- paired <- list()
  if(nrow(combined)) {
    groups <- split(combined,paste(combined$cell_id,combined$method,combined$n_calib,sep="/"))
    for(g in groups) for(metric in intersect(metrics,names(g))) {
      v <- g[[metric]]; v<-v[is.finite(v)]
      c <- cells[[match(g$cell_id[1],vapply(cells,`[[`,"","id"))]]
      summary[[length(summary)+1L]] <- data.frame(cell_id=c$id,method=g$method[1],
        n_calib=g$n_calib[1],metric=metric,mean=if(length(v)) mean(v) else NA_real_,
        mcse=if(length(v)>1L) sd(v)/sqrt(length(v)) else NA_real_,
        n_success=length(v),n_expected=c$n_rep)
    }
  }
  if(nrow(own) && nrow(competitors)) {
    main <- own[own$method=="EIV-GP",,drop=FALSE]
    for(method in unique(competitors$method)) {
      comp <- competitors[competitors$method==method,,drop=FALSE]
      pairs <- merge(main,comp,by=c("cell_id","rep"),suffixes=c("_eiv","_competitor"))
      if(!nrow(pairs)) next
      for(g in split(pairs,paste(pairs$cell_id,pairs$n_calib_eiv))) for(metric in metrics) {
        a<-g[[paste0(metric,"_eiv")]]; b<-g[[paste0(metric,"_competitor")]]
        if(is.null(a)||is.null(b)) next
        delta<-a-b; delta<-delta[is.finite(a)&is.finite(b)]
        paired[[length(paired)+1L]]<-data.frame(cell_id=g$cell_id[1],n_calib=g$n_calib_eiv[1],
          competitor=method,metric=metric,mean_eiv_minus_competitor=if(length(delta)) mean(delta) else NA_real_,
          mcse=if(length(delta)>1L) sd(delta)/sqrt(length(delta)) else NA_real_,n_pairs=length(delta))
      }
    }
  }
  availability <- expected
  availability$mcmc_dataset_available <- key(expected) %in% key(identity)
  availability$competitor_metrics_available <- key(expected) %in% key(competitors)
  list(metrics=combined,summary=mixedgp_bind_rows_base(summary),
       paired=mixedgp_bind_rows_base(paired),competitor_status=status,availability=availability)
}

mixedgp_combine_saved_runs <- function(run_dir, competitor_dir, output_dir) {
  config <- readRDS(file.path(run_dir,"config","resolved_config.rds"))
  provenance <- readRDS(file.path(competitor_dir,"provenance.rds"))
  if(!identical(config$study,provenance$config$study)) stop("Study mismatch.")
  read_csv <- function(p) {x<-read.csv(p,check.names=FALSE); x}
  own <- ids <- diagnostics <- list()
  for(cell in config$cells) {
    p <- file.path(run_dir,"cells",cell$id)
    if(!file.exists(file.path(p,"report_inputs.rds"))) next
    if(!file.exists(file.path(p,"dataset_identity.csv")))
      stop("This older run lacks dataset_identity.csv: ",p,". Keep it; identity must be audited before merging.")
    saved <- readRDS(file.path(p,"report_inputs.rds"))
    if(!identical(saved$study,config$study)) stop("Reporting checkpoint study mismatch.")
    x <- saved$context$mc_results
    if(nrow(x)) {x$cell_id<-cell$id; own[[cell$id]]<-x}
    x <- read_csv(file.path(p,"dataset_identity.csv")); x$cell_id<-cell$id; ids[[cell$id]]<-x
    x <- saved$context[[if(config$study=="study1") "mcmc_diagnostics" else "mc_diagnostics"]]
    if(is.data.frame(x)&&nrow(x)) {x$cell_id<-cell$id; diagnostics[[cell$id]]<-x}
  }
  statuses <- read_csv(file.path(competitor_dir,"statuses.csv"))
  methods <- provenance$methods
  if(is.null(methods)) methods <- c("UC-GP","LVGP","EzGP")
  present <- if(nrow(statuses)) paste(statuses$cell_id,statuses$rep,statuses$method) else character()
  absent <- list()
  for(cell in config$cells) for(method in methods) for(rep in seq_len(cell$n_rep)) {
    if(!paste(cell$id,rep,method) %in% present)
      absent[[length(absent)+1L]] <- data.frame(cell_id=cell$id,rep=rep,method=method,
        status="missing",message="No saved competitor status",dataset_md5=NA_character_)
  }
  statuses <- mixedgp_bind_rows_base(c(list(statuses),absent))
  ans <- mixedgp_combine_predictions(mixedgp_bind_rows_base(own),
    read_csv(file.path(competitor_dir,"metrics.csv")),statuses,
    mixedgp_bind_rows_base(ids),config$cells)
  ans$mcmc_diagnostics <- mixedgp_bind_rows_base(diagnostics)
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  for(name in names(ans)) mixedgp_atomic_write_csv(ans[[name]],file.path(output_dir,paste0(name,".csv")))
  writeLines(if(nrow(ans$summary)) knitr::kable(ans$summary,format="latex",booktabs=TRUE,digits=3) else
    "% No metrics available.",file.path(output_dir,"predictive_summary.tex"))
  saveRDS(list(mcmc_config=config,competitor_provenance=provenance,
    mcmc_run=normalizePath(run_dir),competitor_report=normalizePath(competitor_dir),
    created=Sys.time(),session=sessionInfo()),file.path(output_dir,"provenance.rds"))
  invisible(ans)
}


summarize_comparison <- function(results, metrics=c("RMSE","CRPS","Coverage95","Width95")) {
  frames<-if(is.data.frame(results))list(results) else results[c("prediction","competitors","ablations")]
  rows<-list();idx<-0L
  for(d in frames) {
    if(!is.data.frame(d)||!nrow(d)) next
    if("evaluation_stratum" %in% names(d)) d<-d[is.na(d$evaluation_stratum)|d$evaluation_stratum=="overall",,drop=FALSE]
    for(n in c("cell_id","method","n_calib")) if(!n %in% names(d)) d[[n]]<-NA
    key<-do.call(paste,c(lapply(d[c("cell_id","method","n_calib")],as.character),sep="|"))
    for(k in unique(key)) {
      z<-d[key==k,,drop=FALSE]
      if("rep" %in% names(z)&&anyDuplicated(z$rep)) stop("Duplicate replication IDs in comparison group: ",k)
      for(m in intersect(metrics,names(z))) {
        v<-as.numeric(z[[m]]);v<-v[is.finite(v)];idx<-idx+1L
        rows[[idx]]<-data.frame(z[1,c("cell_id","method","n_calib"),drop=FALSE],metric=m,n=length(v),
          mean=if(length(v))mean(v) else NA_real_,mcse=if(length(v)>1)stats::sd(v)/sqrt(length(v)) else NA_real_)
      }
    }
  }
  if(length(rows))do.call(rbind,rows) else data.frame()
}


# ============================================================================
# ============================================================================

# 03: Combine compatible completed tasks; summarize and plot prediction scores.

# Capture console output, messages and warnings in a per-run log.
# Keep HTML focused on analysis while preserving complete execution logs.


# Import only completed, compatible task payloads; validate before any copying.
# Independent machine roots are preferable to concurrent writers on shared storage.

# Read saved scores without invoking a fit.
mixedgp_bind <- function(parts) {
  if(!length(parts)) return(data.frame())
  cols <- unique(unlist(lapply(parts,names)))
  do.call(rbind,lapply(parts,function(x) {for(n in setdiff(cols,names(x))) x[[n]]<-NA; x[,cols,drop=FALSE]}))
}
# Shared table formatting for CSV/Markdown/LaTeX/HTML. Formatting does not
# change the fold scores or their aggregation; paper coverage is a percentage.
mixedgp_table_align <- function(d) {
  a <- rep("r",ncol(d)); k <- if("Method" %in% names(d)) 3L else if("Class" %in% names(d)) 2L else 1L
  a[seq_len(min(k,ncol(d)))] <- "l"; a
}
# Rank unrounded means within comparable groups. Incomplete groups are not ranked.
mixedgp_html_escape <- function(x) {
  x <- gsub("&","&amp;",as.character(x),fixed=TRUE)
  x <- gsub("<","&lt;",x,fixed=TRUE)
  gsub(">","&gt;",x,fixed=TRUE)
}
# Formatting attributes never enter CSV scores; only the mean is emphasized.
mixedgp_table_markup <- function(d, format=c("html","pipe","latex")) {
  format <- match.arg(format)
  mask <- attr(d,"best")
  out <- d
  for(j in seq_len(ncol(d))) {
    out[[j]] <- if(format=="html") mixedgp_html_escape(d[[j]]) else as.character(d[[j]])
    if(!is.null(mask) && names(d)[j] %in% colnames(mask)) {
      rows <- which(mask[,names(d)[j]])
      for(i in rows) {
        value <- out[[j]][i]
        mean <- sub(" \\(.*$","",value)
        suffix <- substring(value,nchar(mean)+1L)
        marked <- switch(format,html=paste0("<strong>",mean,"</strong>"),
          pipe=paste0("**",mean,"**"),latex=paste0("\\textbf{",mean,"}"))
        out[[j]][i] <- paste0(marked,suffix)
      }
    }
  }
  out
}
# Preserve semantic highlighting when presentation columns are omitted.
mixedgp_table_html <- function(d,caption,note="") {
  table <- as.character(knitr::kable(mixedgp_table_markup(d,"html"),format="html",caption=caption,
    escape=FALSE,row.names=FALSE,align=mixedgp_table_align(d),
    table.attr='class="ocean-summary table" style="width:100%;border-collapse:collapse"'))
  table <- sub("<caption>",'<caption style="text-align:left;padding:8px 0;color:#555">',table,fixed=TRUE)
  table <- sub("<thead>",'<thead style="background:#f2f4f6">',table,fixed=TRUE)
  table <- gsub('style="text-align:',
    'style="padding:8px 10px;white-space:nowrap;border-bottom:1px solid #ddd;text-align:',
    table,fixed=TRUE)
  paste0('<div class="ocean-table-scroll" style="max-width:100%;overflow-x:auto">',table,
    '</div><p class="note">',mixedgp_html_escape(note),"</p>")
}
# Training and test information define each comparison panel, as in the paper.
# Group training/test information in panel headings, keeping the paper LaTeX
# body compact. CSV retains explicit columns for filtering and checking.
# One row per fold combines fit availability, continuation and final screening.

# Response-level evidence reconstructed from saved predictions; never fits.


# Current OI/CC/MI predictions.
# Labels are presentation parameters; never infer fitting inputs from suffixes.
# OI / CC / MI are data-use protocols, with training/test scope defined together.
# Test proportions are calculated from saved masks, not calibration labels.
# Dynamic eligibility: the prescribed folds and computation profiles must match.
# Numeric minima and ties are determined from unrounded losses, never row IDs.


# ============================================================================
# ============================================================================

# 04: U imputation, category means, continuous functions, and prediction diagnostics.
# Loading defines functions only. Direct Rscript execution uses a parameter RDS.

# Final paper-profile diagnostics. Source this file after 01_setup.R; sourcing
# only defines functions. It never fits models, changes masks, or installs packages.

# 1. Frozen inputs and final-fit provenance ------------------------------------





# 2. Training U imputation: five metrics on the common-hidden rows ------------


# Display all training rows in each panel. This is deliberately separate
# from the common-hidden scoring set below: calibration values are NOT imputations.


# 3. Y-blind display domains for m(x,c) and f(x,u) ------------------------------



# 4. Package-only functional posterior calculations, cached separately per fit -

# 5. Original held-out Y tasks; reuse saved scores without new predictions -----

# 6. Figures: one fixed fold in the main text, all folds in the appendix --------

# 7. Rmd entry point: postprocess and export; no fitting ------------------------

# 8. Compact HTML/Rmd presentation: the four diagnostic blocks ----------


# The original three baselines predict a category response directly. Evaluate
# their native means on the fixed X grids; no R/W slice is relabeled as m.


# Input-based display selection. Use the fitted grid and its local
# salinity pool; show at most 3 nearest-salinity observations per log-pressure
# bin, within the curve range. Y, residuals and held-out scores are never used.


# Nested training-calibration roles for descriptive response-function displays.

# Descriptive data analysis: measured nearby responses provide data context.

# Unified postprocessing entry used by both Rmd and this script's CLI.
# No fitting routine is called. Optional copied caches retain their full identity

# Direct execution: Rscript .../04_diagnostics.R diagnostics_settings.rds
# Source this file from Rmd/R without executing the CLI below.
# Compact saved-result checks, separate from fitting and function caches.


# Class-specific marginal imputed-U distributions, a supplement to individual
# interval calibration. Equal weight per common-hidden training profile.

# ---- Public mixing visualization API (saved results only) ----

# Normalize without dropping nonfinite values or changing diagnostic flags.
mixedgp_plot_diagnostics <- function(object) {
  z <- if(inherits(object,"eivgp_diagnostics")) object$table else object
  if(!is.data.frame(z) || !nrow(z)) stop("Supply a nonempty diagnostic data frame or eivgp_diagnostics.")
  if(!"parameter" %in% names(z) && "variable" %in% names(z)) z$parameter <- z$variable
  if(!all(c("parameter","rhat","ess_bulk","ess_tail") %in% names(z)))
    stop("Diagnostics need parameter (or variable), rhat, ess_bulk, and ess_tail columns; run-level gate summaries are insufficient.")
  for(n in c("rhat","ess_bulk","ess_tail")) if(!is.numeric(z[[n]])) stop(n," must be numeric.")
  if(!"mcse_sd_ratio" %in% names(z)) z$mcse_sd_ratio <- NA_real_
  if(!is.numeric(z$mcse_sd_ratio)) stop("mcse_sd_ratio must be numeric.")
  if(!"diagnostic_role" %in% names(z)) z$diagnostic_role <- if("group" %in% names(z)) as.character(z$group) else "Monitored quantities"
  z$diagnostic_role <- as.character(z$diagnostic_role)
  z$diagnostic_role[is.na(z$diagnostic_role)] <- "Unspecified"
  labels <- c(raw="Raw coordinates",target="Scientific targets",invariant="Measurement invariants",
    raw_and_training_imputation="Raw / training imputation",scientific_panel="Scientific targets",
    measurement_invariant="Measurement invariants",user_functional="User functionals")
  known <- z$diagnostic_role %in% names(labels)
  z$diagnostic_role[known] <- unname(labels[z$diagnostic_role[known]])
  if(!"fit_id" %in% names(z)) {
    keys <- intersect(c("scenario","rep","n_calib"),names(z))
    z$fit_id <- if(length(keys)) do.call(paste,c(z[keys],sep=" / ")) else "Fit"
  }
  if(!"panel" %in% names(z)) {
    keys <- intersect(c("scenario","n_calib"),names(z))
    z$panel <- if(length(keys)) do.call(paste,c(lapply(keys,function(k)if(k=="n_calib")paste("calibration",z[[k]]) else z[[k]]),sep=" / ")) else "Diagnostics"
  }
  for(n in c("fit_id","panel")) z[[n]][is.na(z[[n]])] <- "Unspecified"
  z
}

#' Visualize saved MCMC diagnostics with ggplot2
#' @export
plot_mixing <- function(object, type=c("overview","precision","coverage"),
                        rhat_limit=NULL, ess_limit=NULL, mcse_limit=NULL) {
  type <- match.arg(type); z <- mixedgp_plot_diagnostics(object)
  settings <- if(inherits(object,"eivgp_diagnostics")) object$settings else list()
  limit <- function(x,setting,default) {
    if(is.null(x)) x <- if(is.null(setting)) default else setting
    if(length(x)!=1L || !is.numeric(x) || !is.finite(x) || x<=0) stop("Plot limits must be positive finite scalars.")
    x
  }
  rhat_limit <- limit(rhat_limit,settings$rhat_limit,1.01)
  ess_limit <- limit(ess_limit,settings$bulk_ess_limit,400)
  mcse_limit <- limit(mcse_limit,settings$mcse_sd_ratio_limit,.05)
  if(rhat_limit<=1) stop("rhat_limit must exceed 1.")
  extreme <- function(x,fun) if(any(is.finite(x))) fun(x[is.finite(x)]) else NA_real_
  groups <- interaction(z$fit_id,z$panel,z$diagnostic_role,drop=TRUE,lex.order=TRUE)
  s <- do.call(rbind,lapply(split(z,groups),function(q) data.frame(
    fit_id=q$fit_id[1],panel=q$panel[1],role=q$diagnostic_role[1],
    max_rhat=extreme(q$rhat,max),min_bulk=extreme(q$ess_bulk,min),
    max_mcse_sd=extreme(q$mcse_sd_ratio,max),n=nrow(q),
    unavailable=sum(!is.finite(q$rhat)|!is.finite(q$ess_bulk)|!is.finite(q$ess_tail)),
    mcse_unavailable=sum(!is.finite(q$mcse_sd_ratio)))))
  if(type=="overview") {
    p <- ggplot2::ggplot(s,ggplot2::aes(x=.data$max_rhat,y=.data$min_bulk,color=.data$role))+
      ggplot2::geom_vline(xintercept=rhat_limit,linetype=2,color="grey45")+
      ggplot2::geom_hline(yintercept=ess_limit,linetype=2,color="grey45")+
      ggplot2::geom_point(alpha=.75,size=2,na.rm=TRUE)+ggplot2::scale_y_log10()+
      ggplot2::labs(title="Chain agreement and effective sample size",x="Largest finite R-hat",y="Smallest finite bulk ESS (log scale)")
  } else if(type=="precision") {
    p <- ggplot2::ggplot(s,ggplot2::aes(x=.data$max_rhat,y=.data$max_mcse_sd,color=.data$role))+
      ggplot2::geom_vline(xintercept=rhat_limit,linetype=2,color="grey45")+
      ggplot2::geom_hline(yintercept=mcse_limit,linetype=2,color="grey45")+
      ggplot2::geom_point(alpha=.75,size=2,na.rm=TRUE)+
      ggplot2::labs(title="Monte Carlo precision of monitored means",x="Largest finite R-hat",y="Largest finite MCSE / posterior SD")
  } else {
    s <- rbind(transform(s,metric="R-hat / bulk or tail ESS",count=s$unavailable),
               transform(s,metric="Mean MCSE / posterior SD",count=s$mcse_unavailable))
    p <- ggplot2::ggplot(s,ggplot2::aes(x=.data$role,y=.data$count,color=.data$metric))+
      ggplot2::geom_point(position=ggplot2::position_dodge(width=.5),alpha=.6,size=2)+
      ggplot2::labs(title="Unavailable diagnostic values",x=NULL,y="Unavailable quantities per fit",color=NULL)+
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=20,hjust=1))
  }
  p <- p+ggplot2::facet_wrap(~panel,scales="free_x")+ggplot2::theme_minimal(base_size=11)+
    ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),legend.position="bottom")+
    ggplot2::labs(color=NULL,subtitle=paste(nrow(z),"saved diagnostic rows;",length(unique(z$fit_id)),"fits"),
      caption="Extrema use finite values; missing diagnostics are not passes. Reference lines are screens, not convergence guarantees.\nMCSE excludes integration error. See the coverage plot for unavailable values.")
  if(type=="precision" && !any(is.finite(s$max_mcse_sd))) p <- p+
    ggplot2::annotate("text",x=rhat_limit,y=mcse_limit,label="MCSE/SD was not saved; precision cannot be assessed",vjust=-1,size=3)
  p
}

#' Plot retained chains without refitting or computing predictions
#' @export
plot_mcmc <- function(object, parameters=NULL,
                      type=c("trace","rank","acf","running_mean"),
                      max_points=2000L, lag_max=100L, bins=20L) {
  type <- match.arg(type)
  for(v in list(max_points,lag_max,bins)) if(length(v)!=1L || !is.numeric(v) || !is.finite(v) || v<2 || v!=floor(v))
    stop("max_points, lag_max, and bins must be integers >= 2.")
  if(inherits(object,"eivgp_fit") || (is.list(object) && !is.null(object$mcmc) && !is.null(object$mcmc$samples_by_chain))) {
    series <- if(isTRUE(object$sampler_version %in% c("0.3.0","0.3.1"))) mixedgp_v030_raw_series(object) else
      if(mixedgp_fit_engine(object)=="univariate") mixedgp_study1_raw_series(object) else mixedgp_study2_raw_series(object)
  } else series <- object
  if(!is.list(series) || !length(series) || is.null(names(series)) || any(!nzchar(names(series))) || anyDuplicated(names(series)))
    stop("Supply a fit or a uniquely named list of functionals, each containing one numeric vector per chain.")
  if(is.null(parameters)) parameters <- head(names(series),6L)
  if(!is.character(parameters) || !length(parameters) || anyNA(parameters) || anyDuplicated(parameters) || any(!parameters %in% names(series)))
    stop("parameters must name distinct available functionals.")
  rows <- list(); max_stride <- 1L; nc <- NULL
  for(name in parameters) {
    chains <- series[[name]]
    if(!is.list(chains) || length(chains)<2L || !all(vapply(chains,function(v)is.numeric(v)&&is.null(dim(v))&&length(v)>=4L&&all(is.finite(v)),logical(1))))
      stop("Each functional needs >=2 chains of >=4 finite numeric draws; no draws are silently dropped.")
    if(is.null(nc)) nc <- length(chains)
    if(length(chains)!=nc) stop("Functionals must have the same chain count.")
    if(type=="rank") {
      ranks <- rank(unlist(chains,use.names=FALSE),ties.method="average")
      ranks <- split(ranks,rep(seq_along(chains),lengths(chains)))
    }
    for(ch in seq_along(chains)) {
      v <- chains[[ch]]; n <- length(v)
      if(type %in% c("trace","running_mean")) {
        stride <- max(1L,ceiling(n/max_points)); max_stride <- max(max_stride,stride)
        ix <- unique(c(seq.int(1L,n,by=stride),n));y <- if(type=="running_mean")cumsum(v)/seq_len(n) else v
        xx <- ix; yy <- y[ix]
      } else if(type=="rank") {
        xx <- seq_len(bins); yy <- hist(ranks[[ch]],breaks=seq(0,sum(lengths(chains)),length.out=bins+1L),plot=FALSE)$counts/n
      } else {
        xx <- 0:min(lag_max,n-1L)
        yy <- if(stats::sd(v)==0)rep(NA_real_,length(xx)) else as.numeric(stats::acf(v,lag.max=max(xx),plot=FALSE)$acf)
      }
      rows[[length(rows)+1L]] <- data.frame(parameter=name,chain=factor(ch),x=xx,value=yy)
    }
  }
  z <- do.call(rbind,rows)
  p <- ggplot2::ggplot(z,ggplot2::aes(x=.data$x,y=.data$value,color=.data$chain,group=.data$chain))+
    ggplot2::geom_line(alpha=.75,linewidth=.4,na.rm=TRUE)+
    ggplot2::scale_color_manual(values=if(nc<=8L)c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9","#F0E442","#000000")[seq_len(nc)] else grDevices::hcl.colors(nc,"Dark 3"))+ggplot2::facet_wrap(~parameter,scales="free_y")+
    ggplot2::theme_minimal(base_size=11)+ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),legend.position="bottom")+
    ggplot2::labs(title=switch(type,trace="Retained-chain traces",rank="Between-chain pooled-rank balance",acf="Within-chain autocorrelation",running_mean="Cumulative chain means"),
      x=switch(type,rank="Pooled-rank bin",acf="Lag in retained draws","Retained draw index (after warmup)"),
      y=switch(type,rank="Fraction of chain",acf="Autocorrelation",running_mean="Running mean","Value"),color="Chain",
      caption=paste0("Display stride up to ",max_stride,"; calculations use all retained draws. No prediction, refitting, or integration is performed.\nRanks use average ties; discrete/constant ranks require caution. ACF is undefined for constant chains."))
  if(type=="rank") p <- p+ggplot2::geom_hline(yintercept=1/bins,linetype=2,color="grey50")
  if(type=="acf" && anyNA(z$value)) p <- p+ggplot2::labs(subtitle="Constant-chain ACFs are unavailable and omitted; inspect the trace plot.")
  p
}

#' Save a small default mixing report from diagnostic CSV files
#' @export
save_mixing_plots <- function(path, output_dir=file.path(path,"mixing-plots"),
                              formats=c("pdf","png")) {
  if(length(path)!=1L || !dir.exists(path)) stop("path must be an existing results directory.")
  if(!length(formats) || any(!formats %in% c("pdf","png"))) stop("formats must contain pdf and/or png.")
  trace_index <- list.files(path,pattern="^scientific-trace-.*[.]csv$",recursive=TRUE,full.names=TRUE)
  if(length(trace_index) && all(grepl("^scientific-trace-cut-",basename(trace_index))))
    return(mixedgp_cut_trace_report(path,output_dir,formats))
  if(length(trace_index))return(mixedgp_save_scientific_trace_plots(trace_index,output_dir,formats))
  files <- list.files(path,pattern="[.]csv$",recursive=TRUE,full.names=TRUE)
  files <- files[!grepl("/\\._",files)]
  names <- basename(files)
  preferred <- grepl("(^|_)mcmc_(parameter|target)_diagnostics[.]csv$|^measurement_diagnostics[.]csv$|^conditional-[0-9]+[.]diagnostics[.]csv$|^full_diagnostics[.]csv$",names)
  fallback <- grepl("(^|_)mcmc_diagnostic_details[.]csv$",names) & !dirname(files) %in% dirname(files[preferred])
  files <- files[preferred|fallback]
  if(!length(files)) {
    warning("No parameter-level diagnostic CSVs found; aggregate gate summaries cannot replace them. No plots written.",call.=FALSE)
    return(invisible(character()))
  }
  tables <- lapply(files,function(f) {
    z <- utils::read.csv(f,stringsAsFactors=FALSE)
    if(!nrow(z)) return(NULL)
    z <- mixedgp_plot_diagnostics(z)
    relative <- substring(f,nchar(sub("/+$","",path))+2L)
    role <- if(grepl("_target_",basename(f)))"Scientific targets" else if(grepl("_parameter_",basename(f)))"Raw coordinates" else if(grepl("measurement",basename(f)))"Measurement" else "Saved quantities"
    if(all(z$diagnostic_role=="Monitored quantities")) z$diagnostic_role <- role
    leaf <- dirname(relative)
    # Parameter and target CSVs describe the same fits; conditional files do not.
    id <- if(grepl("^conditional-",basename(f))) relative else leaf
    z$fit_id <- paste(id,z$fit_id,sep=" / ")
    z$panel <- if(grepl("^(shared-)?rep[0-9]+-cal",basename(leaf)))
      paste(basename(dirname(leaf)),sub(".*cal([0-9]+).*","cal \\1",basename(leaf)),role) else if(leaf==".") z$panel else paste(leaf,z$panel,sep=" / ")
    z[,c("parameter","rhat","ess_bulk","ess_tail","mcse_sd_ratio","diagnostic_role","fit_id","panel")]
  })
  z <- do.call(rbind,tables)
  if(is.null(z)||!nrow(z)) {warning("Diagnostic files are empty; no plots written.",call.=FALSE);return(invisible(character()))}
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  jobs <- expand.grid(ext=unique(formats),type=c("overview","precision","coverage"),stringsAsFactors=FALSE)
  paths <- unlist(mixedgp_stage_lapply(seq_len(nrow(jobs)),function(i) {
    dest <- file.path(output_dir,paste0("mixing_",jobs$type[i],".",jobs$ext[i]))
    ggplot2::ggsave(dest,plot_mixing(z,jobs$type[i]),width=11,height=7,dpi=160,bg="white")
    dest
  },workers=1L),use.names=FALSE)
  utils::write.csv(data.frame(path=normalizePath(files),md5=unname(tools::md5sum(files))),
    file.path(output_dir,"diagnostic_sources.csv"),row.names=FALSE)
  invisible(paths)
}

# Each figure/format owns its graphics device and destination.
mixedgp_save_plot_jobs <- function(plots,dest,dimensions) {
  jobs <- expand.grid(ext=c("png","pdf"),name=names(plots),stringsAsFactors=FALSE)
  if(!nrow(jobs))return(invisible(character()))
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  paths <- mixedgp_stage_lapply(seq_len(nrow(jobs)),function(i) {
    name <- jobs$name[i];size <- dimensions(name)
    path <- file.path(dest,paste0(name,".",jobs$ext[i]))
    # Quartz is unsafe in forked macOS workers. Use Cairo for PNG and base PDF.
    device <- if(jobs$ext[i]=="png")function(filename,width,height,bg,...) {
      grDevices::png(filename,width=width,height=height,units="in",res=180,bg=bg,
        type=if(capabilities("cairo"))"cairo" else getOption("bitmapType"))
    } else grDevices::pdf
    ggplot2::ggsave(path,plots[[name]],device=device,width=size[1],height=size[2],dpi=180,bg="white")
    path
  },workers=if(Sys.info()[["sysname"]]=="Darwin" && !capabilities("cairo"))1L else mixedgp_budget_cores())
  invisible(unlist(paths,use.names=FALSE))
}

# Stable manuscript labels; every artifact is tied to the selected cohort.

# Split the appendix's point and distributional metrics into readable panels.

# Scientific-target diagnostics are postprocessing, never a convergence gate.




# Manuscript-ready LaTeX is derived from the same unrounded-score exports as HTML.
# Raw CSVs retain all conditions, fold counts and missing-method information.

# Read-only audit of registered calibrated numerical jobs. Output is always separate.
# Named roots sharing a design label may be paired only after the caller verifies
# identical datasets, test targets, calibration masks and scoring protocols.
report_numerical_convergence <- function(roots, output) {
  if(is.null(names(roots)) || any(!nzchar(names(roots)))) stop("Name each root with its verified study/design identifier")
  roots<-vapply(roots,normalizePath,character(1),mustWork=TRUE)
  canonical<-function(p){p<-path.expand(p);if(file.exists(p))normalizePath(p,mustWork=TRUE)else file.path(canonical(dirname(p)),basename(p))}
  output<-canonical(output)
  if(any(output==roots | startsWith(output,paste0(roots,.Platform$file.sep))))
    stop("Report output must be outside all experiment roots")
  read_errors<-list()
  read<-function(p)if(file.exists(p))tryCatch(readRDS(p),error=function(e){
    read_errors[[length(read_errors)+1L]]<<-data.frame(path=p,message=conditionMessage(e));NULL
  })else NULL
  diag_state<-function(p) {
    if(!file.exists(p))return("unknown")
    z<-tryCatch(read.csv(p),error=function(e){read_errors[[length(read_errors)+1L]]<<-data.frame(path=p,message=conditionMessage(e));NULL})
    if(is.null(z))return("unknown")
    if(!all(c("rhat","ess_bulk","ess_tail")%in%names(z)))return("unknown")
    if(nrow(z)>0 && all(is.finite(z$rhat)&z$rhat<=1.01&is.finite(z$ess_bulk)&z$ess_bulk>=400&is.finite(z$ess_tail)&z$ess_tail>=400))"passed"else"failed"
  }
  rows<-scores<-list()
  for(i in seq_along(roots)) {
    root<-roots[i];jobs<-read(file.path(root,"job_registry.rds"))
    if(is.null(jobs))stop("Missing job_registry.rds: ",root)
    for(j in seq_len(nrow(jobs))) {
      job<-jobs[j,];folder<-file.path(root,job$job)
      # Atomic score publication is the completion marker; partial runs remain pending.
      z<-read(file.path(folder,"scores.rds"));err<-read(file.path(folder,"failure.rds"))
      if(!is.null(z))err<-NULL
      pred<-read(file.path(folder,"predictions.rds"));dg<-read(file.path(folder,"diagnostics.rds"))
      cc<-grepl("-CC$",job$method);native<-job$method%in%c("UC-GP","LVGP","EzGP")||cc
      stage<-if(is.null(err$stage))"unknown"else as.character(err$stage)
      opt<-if(!native)"not_applicable"else"unknown"
      os<-pred$optimization_status
      if(is.null(os))os<-dg$optimizer_status
      if(native&&!is.null(os)&&length(os)==1&&!is.na(os)) {
        if(os%in%c("converged","default_converged","rescued_converged"))opt<-"passed"
        else if(grepl("nonconverged",os,fixed=TRUE))opt<-"nonconverged"
      }
      if(native&&stage=="optimizer")opt<-"error"
      measurement<-if(cc||job$method=="CUT")diag_state(file.path(root,sprintf("shared-rep%03d-cal%02d",job$rep,job$calib),"measurement_diagnostics.csv"))else"not_applicable"
      if(stage=="measurement")measurement<-"error"
      response<-"not_applicable"
      if(job$method=="FULL")response<-diag_state(file.path(folder,"full_diagnostics.csv"))
      if(job$method=="CUT") {
        # Combined gate covers every selected conditional response fit. A positive
        # gate plus independently verified measurement screen establishes both.
        response<-if(isTRUE(dg$pass)&&measurement=="passed")"passed"else"unknown"
        checks<-dg$checks
        if(length(checks)>1)response<-if(all(!is.na(checks[-1])&checks[-1]))"passed"else"failed"
      }
      state<-if(!is.null(z))"complete"else if(!is.null(err))"failed"else if(any(file.exists(file.path(folder,c("scores.rds","failure.rds")))))"artifact_error"else"pending"
      passed<-state=="complete"&&all(c(opt,measurement,response)%in%c("passed","not_applicable"))
      row<-data.frame(design=names(roots)[i],rep=job$rep,calib=job$calib,method=job$method,
        status=state,optimizer=opt,measurement=measurement,response=response,
        converged=passed,error_stage=if(is.null(err))""else stage,
        message=if(is.null(err))""else err$message,source=folder,stringsAsFactors=FALSE)
      rows[[length(rows)+1]]<-row
      if(!is.null(z)) {
        z$coordinate[is.na(z$coordinate)]<-0L
        z<-rbind(z,transform(z,stratum="overall"))
        gs<-split(z,interaction(z$target,z$stratum,z$coordinate,drop=TRUE))
        for(d in gs) {
          r<-row[,c("design","rep","calib","method","converged")]
          r$target<-d$target[1];r$stratum<-d$stratum[1];r$coordinate<-d$coordinate[1]
          for(k in c("sq_error","crps","nlpd","coverage","width"))r[[k]]<-if(all(is.finite(d[[k]])))mean(d[[k]])else NA_real_
          r$rmse<-sqrt(r$sq_error);scores[[length(scores)+1]]<-r
        }
      }
    }
  }
  status<-do.call(rbind,rows)
  if(is.null(status))stop("No registered jobs")
  if(anyDuplicated(status[c("design","rep","calib","method")]))stop("Duplicate method/replication keys across sources")
  rates<-do.call(rbind,lapply(split(status,interaction(status$design,status$calib,status$method,drop=TRUE)),function(d){
    o<-d[1,c("design","calib","method")];o$n_planned<-nrow(d);o$n_resolved<-sum(d$status%in%c("complete","failed"))
    o$n_pending<-sum(d$status=="pending");o$n_hard_failure<-sum(d$status=="failed");o$n_complete<-sum(d$status=="complete")
    o$n_artifact_error<-sum(d$status=="artifact_error")
    o$n_converged<-sum(d$converged)
    o$hard_failure_rate<-if(o$n_resolved)o$n_hard_failure/o$n_resolved else NA_real_
    for(k in c("optimizer","measurement","response")) {
      known<-d[[k]]%in%c("passed","failed","nonconverged","error")
      fail<-d[[k]]%in%c("failed","nonconverged","error")
      o[[paste0(k,"_n_nonconverged")]]<-sum(d[[k]]=="nonconverged")
      o[[paste0(k,"_n_error")]]<-sum(d[[k]]=="error")
      o[[paste0(k,"_n_assessed")]]<-sum(known);o[[paste0(k,"_n_failed")]]<-sum(fail)
      o[[paste0(k,"_n_unknown")]]<-sum(d[[k]]=="unknown")
      o[[paste0(k,"_failure_rate")]]<-if(any(known))sum(fail)/sum(known)else NA_real_
    };o
  }))
  reps<-if(length(scores))do.call(rbind,scores)else data.frame()
  pairs<-list();summaries<-list();metrics<-c("rmse","crps","nlpd","coverage","width")
  if(nrow(reps))for(des in unique(reps$design)) {
    rr<-reps[reps$design==des & reps$converged,,drop=FALSE];methods<-sort(unique(rr$method))
    if(length(methods)>1)for(a in seq_len(length(methods)-1))for(b in (a+1):length(methods)) {
      p<-merge(rr[rr$method==methods[a],],rr[rr$method==methods[b],],by=c("design","rep","calib","target","stratum","coordinate"))
      for(k in metrics) {
        d<-p[is.finite(p[[paste0(k,".x")]])&is.finite(p[[paste0(k,".y")]]),,drop=FALSE]
        if(!nrow(d))next
        d$metric<-k;d$value_a<-d[[paste0(k,".x")]];d$value_b<-d[[paste0(k,".y")]];d$difference<-d$value_a-d$value_b
        d<-d[,c("design","rep","calib","target","stratum","coordinate","method.x","method.y","metric","value_a","value_b","difference")]
        names(d)[7:8]<-c("method_a","method_b");pairs[[length(pairs)+1]]<-d
        for(g in split(d,interaction(d$calib,d$target,d$stratum,d$coordinate,drop=TRUE))) {
          s<-g[1,c("design","calib","target","stratum","coordinate","method_a","method_b","metric")]
          s$n_matched<-nrow(g);s$mean_a<-mean(g$value_a);s$mean_b<-mean(g$value_b);s$difference<-mean(g$difference)
          s$mcse<-if(nrow(g)>1)sd(g$difference)/sqrt(nrow(g))else NA_real_
          summaries[[length(summaries)+1]]<-s
        }
      }
    }
  }
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  tables<-list(read_errors=if(length(read_errors))do.call(rbind,read_errors)else data.frame(path=character(),message=character()),status=status,failure_rates=rates,replicate_scores=reps,
    matched_converged_pairs=if(length(pairs))do.call(rbind,pairs)else data.frame(),
    matched_converged_summary=if(length(summaries))do.call(rbind,summaries)else data.frame())
  for(n in names(tables))write.csv(tables[[n]],file.path(output,paste0(n,".csv")),row.names=FALSE)
  writeLines(c("Diagnostic screening is not proof of convergence.",
    "Failure rates use assessed-stage denominators; unknown, not-applicable and pending remain separate.",
    "Hard-failure rate uses resolved jobs. Rates for shared measurement fits are method-job exposure rates, not independent fits.",
    "Matched results are pairwise, metric-specific complete cases passing every applicable diagnostic; they can be selection-biased.",
    "MCSE describes replication variability, not within-fit error. Unfiltered replicate scores are retained.",
    "Equal design labels assert verified matching data, masks, targets and protocols across roots.",
    "This snapshot reads saved artifacts only; source folders were not modified."),file.path(output,"README.txt"))
  invisible(tables)
}

# Study I saved-result reporting. Numerical inference is never invoked here.
mixedgp_study1_report_run_dir <- function(path) {
  path <- normalizePath(path,mustWork=TRUE)
  if(file.exists(file.path(path,'config/resolved_config.rds')))return(path)
  candidates <- unique(c(list.dirs(path,recursive=FALSE,full.names=TRUE),
    if(dir.exists(file.path(path,'core')))list.dirs(file.path(path,'core'),recursive=FALSE,full.names=TRUE)))
  candidates <- candidates[file.exists(file.path(candidates,'config/resolved_config.rds'))]
  if(length(candidates)!=1L)stop('Supply one exact Study I run directory; found ',length(candidates),' runs under ',path)
  candidates
}

mixedgp_study1_score_summary <- function(x,metrics) {
  if(!nrow(x))return(data.frame())
  keys<-intersect(c('cell_id','n_calib','method'),names(x))
  x[keys]<-lapply(x[keys],function(v)ifelse(is.na(v),'Not applicable',as.character(v)))
  groups<-interaction(x[keys],drop=TRUE,lex.order=TRUE)
  do.call(rbind,lapply(split(x,groups),function(z) {
    row<-z[1,keys,drop=FALSE]
    for(m in intersect(metrics,names(z))) {
      v<-z[[m]][is.finite(z[[m]])];row[[m]]<-if(length(v))mean(v)else NA_real_
      row[[paste0(m,'_MCSE')]]<-if(length(v)>1)stats::sd(v)/sqrt(length(v))else NA_real_
      row[[paste0(m,'_n')]]<-length(v)
    };row
  }))
}

mixedgp_study1_target_ess_summary <- function(z,ess_min=100) {
  stopifnot(is.finite(ess_min),ess_min>0)
  do.call(rbind,lapply(split(z,interaction(z$cell_id,z$n_calib,z$family,drop=TRUE)),function(x) {
    finite_summary<-function(v,fun)if(any(is.finite(v)))fun(v[is.finite(v)])else NA_real_
    unavailable<-!is.finite(x$ess_bulk)|!is.finite(x$ess_tail)
    low<-!unavailable & (x$ess_bulk<=ess_min|x$ess_tail<=ess_min)
    data.frame(cell_id=x$cell_id[1],n_calib=x$n_calib[1],family=x$family[1],n_quantities=nrow(x),
      median_bulk_ess=finite_summary(x$ess_bulk,median),min_bulk_ess=finite_summary(x$ess_bulk,min),
      min_tail_ess=finite_summary(x$ess_tail,min),max_rhat=finite_summary(x$rhat,max),
      max_mcse_sd_ratio=finite_summary(x$mcse_sd_ratio,max),below_ess= sum(low),unavailable=sum(unavailable),
      affected_fits=length(unique(x$rep[low|unavailable])))
  }))
}

mixedgp_study1_trace_selection <- function(z,run_dir) {
  rows<-list()
  for(cell in unique(z$cell_id))for(ca in sort(unique(z$n_calib[z$cell_id==cell]))) {
    a<-z[z$cell_id==cell & z$n_calib==ca,]
    folder<-file.path(run_dir,'cells',cell,'results/study1_publication')
    a$trace_file<-file.path(folder,sprintf('target_traces_rep%03d_cal%03d.rds',a$rep,ca))
    a$fit_file<-file.path(folder,sprintf('mcmc_failure_rep%03d_cal%03d.rds',a$rep,ca))
    current<-file.path(folder,sprintf('fit_rep%03d_cal%03d.rds',a$rep,ca))
    a$fit_file[!file.exists(a$fit_file)]<-current[!file.exists(a$fit_file)]
    a$trace_available<-file.exists(a$trace_file)|file.exists(a$fit_file)
    m<-a[a$family=='m_conditional_mean',]
    mins<-aggregate(ess_bulk~rep,m,min);mid<-stats::median(mins$ess_bulk)
    mins<-mins[order(abs(log(mins$ess_bulk/mid)),mins$rep),]
    eligible<-mins$rep[mins$rep %in% m$rep[m$trace_available]]
    b<-if(length(eligible))m[m$rep==eligible[1],]else m[m$rep==mins$rep[1],]
    b<-b[which.min(b$ess_bulk),];b$selection<-'Typical mean';rows[[length(rows)+1L]]<-b
    for(kind in c('Weakest mean','Weakest uncertainty')) {
      families<-if(kind=='Weakest mean')c('m_conditional_mean','f_conditional_mean')else
        c('m_conditional_variance','f_conditional_variance','predictive_second_moment')
      b<-a[a$family %in% families,];b<-b[order(pmin(b$ess_bulk,b$ess_tail),-b$rhat),][1,,drop=FALSE]
      b$selection<-kind;rows[[length(rows)+1L]]<-b
    }
  }
  selected<-do.call(rbind,rows)
  sampling_file<-file.path(run_dir,'combined/mcmc_diagnostics.csv')
  if(file.exists(sampling_file)) {
    sampling<-utils::read.csv(sampling_file)
    continued<-sampling[is.finite(sampling$n_iter)&sampling$n_iter>20000L,,drop=FALSE]
    extra<-list()
    if(nrow(continued))for(i in seq_len(nrow(continued))) {
      key<-continued[i,];a<-z[z$cell_id==key$cell_id & z$rep==key$rep & z$n_calib==key$n_calib,]
      folder<-file.path(run_dir,'cells',key$cell_id,'results/study1_publication')
      for(family in c('m_conditional_mean','predictive_second_moment')) {
        b<-a[a$family==family,];if(!nrow(b))next
        b<-b[which.min(pmin(b$ess_bulk,b$ess_tail)),,drop=FALSE]
        b$trace_file<-file.path(folder,sprintf('target_traces_rep%03d_cal%03d.rds',key$rep,key$n_calib))
        b$fit_file<-file.path(folder,sprintf('mcmc_failure_rep%03d_cal%03d.rds',key$rep,key$n_calib))
        if(!file.exists(b$fit_file))b$fit_file<-file.path(folder,sprintf('fit_rep%03d_cal%03d.rds',key$rep,key$n_calib))
        b$trace_available<-file.exists(b$trace_file)|file.exists(b$fit_file)
        b$selection<-'Continued scientific targets';extra[[length(extra)+1L]]<-b
      }
    }
    if(length(extra))selected<-bind_rows(selected[selected$selection!='Weakest mean',],bind_rows(extra))
  }
  selected
}

mixedgp_study1_trace_data <- function(row,run_dir,data_root,cache_dir,max_points=751L) {
  cache<-file.path(cache_dir,sprintf('%s_rep%03d_cal%03d.rds',row$cell_id,row$rep,row$n_calib))
  # Small direct target traces are saved for all new fits, including diagnostically strong fits.
  if(file.exists(row$trace_file))return(mixedgp_read_artifact(row$trace_file))
  if(!file.exists(row$fit_file))stop('Retained chains unavailable for ',row$cell_id,' rep ',row$rep,' calibration ',row$n_calib)
  source_md5<-unname(tools::md5sum(row$fit_file))
  data_file<-file.path(data_root,row$cell_id,sprintf('study1_heterogeneity_continuum_n100_ntest100_m6_rep%03d.rds',row$rep))
  if(!file.exists(data_file))stop('Frozen dataset missing: ',data_file)
  manifest<-utils::read.csv(file.path(run_dir,'config/input_manifest.csv'))
  expected<-manifest[manifest$cell_id==row$cell_id & manifest$rep==row$rep,]
  data_md5<-unname(tools::md5sum(data_file))
  if(nrow(expected)!=1L || !identical(data_md5,expected$md5))stop('Frozen input checksum mismatch: ',data_file)
  identity<-list(schema=1L,fit=source_md5,data=data_md5,max_points=max_points,
    mcmc=unname(eivgp_module_identity()['mcmc.R']))
  if(file.exists(cache)) {
    old<-readRDS(cache);if(identical(old$identity,identity))return(old)
  }
  message('Reading retained target chains: ',row$cell_id,' rep ',row$rep,' calibration ',row$n_calib)
  saved<-mixedgp_read_artifact(row$fit_file);fit<-saved$fit
  test<-readRDS(data_file)$data$test;panel<-saved$diagnostic_panel_rows
  n<-nrow(fit$mcmc$samples_by_chain$u[[1]])
  keep<-unique(c(seq.int(1L,n,by=max(1L,ceiling(n/(max_points-1L)))),n))
  fit$mcmc$samples_by_chain<-lapply(fit$mcmc$samples_by_chain,function(chains)lapply(chains,function(v) {
    if(is.null(dim(v)))v[keep]else v[keep,,drop=FALSE]
  }))
  eta<-saved$cache_spec$heterogeneity_eta
  series<-mixedgp_study1_target_series(fit,X=matrix(test$x[panel],ncol=1),C=test$c[panel],
    U=if(eta==1)test$u[panel]else NULL,n_latent=saved$cache_spec$diagnostic_panel$n_latent,
    seed=380000L+1000L*row$rep+row$n_calib)
  # These thinned series are for display only; ESS comes from the full-window CSV.
  series<-mixedgp_mark_diagnostic_window(series,FALSE)
  burn<-saved$cache_spec$mcmc$burn;thin<-fit$diagnostics$summary$thin
  if(is.null(thin))thin<-1L
  result<-list(identity=identity,series=series,iteration=burn+keep*thin,
    retained_index=keep,panel_rows=panel,test_panel=test[panel,],source=row$fit_file)
  saveRDS(result,cache);result
}

mixedgp_study1_save_trace_figures <- function(selection,run_dir,data_root,output_dir,ess_min) {
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  cache_dir<-file.path(output_dir,'trace-data');dir.create(cache_dir,showWarnings=FALSE)
  all_rows<-list();status<-list();chain_summary<-list()
  for(i in seq_len(nrow(selection))) {
    row<-selection[i,]
    tryCatch({
      trace<-mixedgp_study1_trace_data(row,run_dir,data_root,cache_dir)
      s<-trace$series[[row$parameter]]
      if(is.null(s)||length(s)!=4L||any(!is.finite(unlist(s))))stop('Missing or invalid target chain series')
      label<-paste0(row$cell_id,'; |O|=',row$n_calib,'; rep ',row$rep,'\n',
        row$parameter,'\nBulk ESS ',round(row$ess_bulk),'; tail ESS ',round(row$ess_tail),
        '; R-hat ',sprintf('%.3f',row$rhat))
      all_rows[[i]]<-do.call(rbind,lapply(seq_along(s),function(ch)data.frame(
        selection=row$selection,panel=label,iteration=trace$iteration,value=s[[ch]],chain=factor(ch))))
      chain_summary[[i]]<-do.call(rbind,lapply(seq_along(s),function(ch) {
        v<-s[[ch]];half<-floor(length(v)/2)
        data.frame(cell_id=row$cell_id,n_calib=row$n_calib,rep=row$rep,parameter=row$parameter,
          selection=row$selection,chain=ch,display_mean=mean(v),first_half_mean=mean(v[seq_len(half)]),
          second_half_mean=mean(v[(half+1):length(v)]),display_sd=stats::sd(v))
      }))
      status[[i]]<-data.frame(row[,c('cell_id','n_calib','rep','parameter','selection')],status='written',message='')
    },error=function(e)status[[i]]<<-data.frame(row[,c('cell_id','n_calib','rep','parameter','selection')],
      status='unavailable',message=conditionMessage(e)))
  }
  status<-do.call(rbind,status);utils::write.csv(status,file.path(output_dir,'trace_status.csv'),row.names=FALSE)
  if(length(chain_summary))utils::write.csv(do.call(rbind,chain_summary),file.path(output_dir,'display_chain_summaries.csv'),row.names=FALSE)
  z<-do.call(rbind,all_rows);paths<-character()
  if(!is.null(z))for(kind in unique(selection$selection)) {
    d<-z[z$selection==kind,];if(!nrow(d))next
    d$panel<-factor(d$panel,levels=unique(d$panel))
    p<-ggplot2::ggplot(d,ggplot2::aes(iteration,value,color=chain))+
      ggplot2::geom_line(linewidth=.25,alpha=.65)+ggplot2::facet_wrap(~panel,ncol=2,scales=if(kind=='Continued scientific targets')'free'else'free_y')+
      ggplot2::scale_color_manual(values=c('#0072B2','#D55E00','#009E73','#CC79A7'))+
      ggplot2::theme_minimal(base_size=11)+ggplot2::theme(legend.position='bottom',panel.grid.minor=ggplot2::element_blank())+
      ggplot2::labs(title=paste('Study I:',tolower(kind),'target traces'),x='MCMC iteration (after warmup)',y='Target value',color='Chain',
        subtitle='Four retained chains; target values reconstructed without refitting',
        caption=paste0('Traces display evenly spaced retained states across the full window. ESS and R-hat use all saved draws.\nESS criterion: bulk and tail ESS > ',ess_min,'. Panel y-scales differ. Selection is documented in trace_selection.csv.'))
    if(kind=='Continued scientific targets')p<-p+ggplot2::geom_vline(xintercept=20000,linetype=2,color='gray45')
    base<-file.path(output_dir,paste0('mixing_',gsub(' ','_',tolower(kind))))
    # Never fork graphics: on macOS even PDF font initialization may enter Objective-C.
    for(ext in c('pdf','png')) {
      dest<-paste0(base,'.',ext)
      device<-if(ext=='pdf')grDevices::pdf else function(filename,width,height,bg,...) {
        grDevices::png(filename,width=width,height=height,units='in',res=150,bg=bg,
          type=if(capabilities('cairo'))'cairo'else getOption('bitmapType'))
      }
      ggplot2::ggsave(dest,p,device=device,width=11,height=3*ceiling(length(unique(d$panel))/2)+1.6,bg='white')
      paths<-c(paths,dest)
    }
  }
  if(length(paths) && all(status$status=='written')) {
    known<-as.vector(outer(c('mixing_typical_mean','mixing_weakest_mean','mixing_weakest_uncertainty','mixing_continued_scientific_targets'),c('.pdf','.png'),paste0))
    stale<-setdiff(file.path(output_dir,known),paths)
    unlink(stale[file.exists(stale)])
  }
  list(paths=paths,status=status)
}

#' Report saved Study I EIV-GP results and scientific-target mixing
#' @export
report_study1_eivgp <- function(path,output_dir=NULL,data_root=NULL,ess_min=100,mixing_plots=TRUE) {
  run_dir<-mixedgp_study1_report_run_dir(path)
  config<-readRDS(file.path(run_dir,'config/resolved_config.rds'))
  if(!identical(config$study,'study1'))stop('Expected a Study I run')
  if(is.null(output_dir))output_dir<-file.path(dirname(run_dir),'study1-report')
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE);output_dir<-normalizePath(output_dir)
  if(identical(output_dir,run_dir)||startsWith(output_dir,paste0(run_dir,'/')))stop('Use a report directory outside the saved run')
  if(is.null(data_root))data_root<-config$data_root
  read<-function(name) {
    f<-file.path(run_dir,name);if(!file.exists(f))return(data.frame())
    utils::read.csv(f,stringsAsFactors=FALSE)
  }
  files<-list.files(file.path(run_dir,'cells'),pattern='^study1_mcmc_target_diagnostics[.]csv$',recursive=TRUE,full.names=TRUE)
  if(!length(files))stop('No saved scientific-target diagnostics')
  targets<-do.call(rbind,lapply(files,function(f) {
    z<-utils::read.csv(f);z$cell_id<-basename(dirname(dirname(f)));z$family<-sub('\\[.*','',z$parameter);z
  }))
  targets$ess_acceptable<-is.finite(targets$ess_bulk)&is.finite(targets$ess_tail)&targets$ess_bulk>ess_min&targets$ess_tail>ess_min
  ess_summary<-mixedgp_study1_target_ess_summary(targets,ess_min)
  review<-targets[!targets$ess_acceptable,];review<-review[order(pmin(review$ess_bulk,review$ess_tail)),]
  summaries<-list(
    prediction=mixedgp_study1_score_summary(read('combined/predictive_metrics.csv'),c('RMSE','CRPS','NLPD','Coverage95','Width95','IntervalScore95')),
    observed_mean=mixedgp_study1_score_summary(read('combined/mean_recovery.csv'),c('RMSE','Bias','Coverage95','Width95')),
    latent_surface=mixedgp_study1_score_summary(read('combined/surface_recovery.csv'),c('ISE','Bias','Coverage95','Width95')),
    latent_imputation=mixedgp_study1_score_summary(read('combined/latent_imputation.csv'),c('RMSE','Bias','Coverage95','Width95')),
    predictive_controls=mixedgp_study1_score_summary(read('combined/ablation_predictive_metrics.csv'),c('RMSE','CRPS','Coverage95','Width95')),
    surface_controls=mixedgp_study1_score_summary(read('combined/ablation_surface_recovery.csv'),c('ISE','Coverage95','Width95')))
  selection<-mixedgp_study1_trace_selection(targets,run_dir)
  utils::write.csv(selection,file.path(output_dir,'trace_selection.csv'),row.names=FALSE)
  figures<-if(isTRUE(mixing_plots))mixedgp_study1_save_trace_figures(selection,run_dir,data_root,file.path(output_dir,'mixing-plots'),ess_min)else list(paths=character(),status=data.frame())
  source_files<-c(files,list.files(file.path(run_dir,'combined'),pattern='[.]csv$',full.names=TRUE))
  provenance<-data.frame(path=normalizePath(source_files),md5=unname(tools::md5sum(source_files)))
  inspection<-NULL;review_file<-file.path(output_dir,'trace_review.rds')
  if(file.exists(review_file)) { reviewed<-readRDS(review_file);if(identical(reviewed$provenance,provenance))inspection<-reviewed$text }
  result<-list(run_dir=run_dir,output_dir=output_dir,config=config,ess_min=ess_min,
    completion=read('config/cell_status.csv'),completeness=read('config/diagnostic_gates.csv'),
    sampling=read('combined/mcmc_diagnostics.csv'),
    summaries=summaries,target_summary=ess_summary,targets=targets,review=review,
    selection=selection,figures=figures,provenance=provenance,
    ablation_status=read('combined/ablation_status.csv'),created_at=format(Sys.time(),tz='UTC'),
    reporting_module_identity=eivgp_module_identity(),inspection=inspection)
  utils::write.csv(ess_summary,file.path(output_dir,'target_ess_summary.csv'),row.names=FALSE)
  utils::write.csv(review,file.path(output_dir,'targets_needing_review.csv'),row.names=FALSE)
  utils::write.csv(provenance,file.path(output_dir,'source_manifest.csv'),row.names=FALSE)
  for(name in names(summaries))utils::write.csv(summaries[[name]],file.path(output_dir,paste0(name,'_summary.csv')),row.names=FALSE)
  saveRDS(result,file.path(output_dir,'report.rds'))
  result
}


mixedgp_save_scientific_trace <- function(fit, prefix, label) {
  z<-fit$scientific_ess;tab<-z$diagnostics
  groups<-list(means=grep("^(m|f)_conditional_mean",tab$parameter),
    uncertainty=grep("variance|second_moment",tab$parameter))
  selected<-vapply(groups,function(rows) {
    score<-pmin(tab$ess_bulk[rows],tab$ess_tail[rows]);
    if(!any(is.finite(score)))rows[1L]else rows[which.min(replace(score,!is.finite(score),Inf))]
  },integer(1))
  tab<-tab[selected,,drop=FALSE];tab$figure<-names(groups);tab$label<-label
  tab$scope<-if(is.null(z$scope))"full_history"else z$scope
  mixedgp_atomic_save_rds(list(table=tab,series=z$target_series[tab$parameter],
    iteration=seq.int(fit$checkpoint$burn+1L,fit$checkpoint$iteration)),paste0(prefix,".rds"))
  utils::write.csv(tab,paste0(prefix,".csv"),row.names=FALSE)
}

mixedgp_save_scientific_trace_plots <- function(index,output_dir,formats=c("png","pdf")) {
  tab<-do.call(rbind,lapply(index,function(p){z<-utils::read.csv(p);z$file<-sub("[.]csv$",".rds",p);z}))
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  paths<-character()
  for(kind in unique(tab$figure)) {
    z<-tab[tab$figure==kind,,drop=FALSE]
    z<-head(z[order(pmin(z$ess_bulk,z$ess_tail)),,drop=FALSE],4L)
    dat<-do.call(rbind,lapply(seq_len(nrow(z)),function(i) {
      row<-z[i,];a<-readRDS(row$file);keep<-unique(as.integer(round(seq(1,length(a$iteration),length.out=min(1500,length(a$iteration))))))
      title<-sprintf("%s\n%s | bulk %.1f | tail %.1f | R-hat %.3f",row$label,row$parameter,row$ess_bulk,row$ess_tail,row$rhat)
      do.call(rbind,lapply(seq_along(a$series[[row$parameter]]),function(ch)
        data.frame(iteration=a$iteration[keep],value=a$series[[row$parameter]][[ch]][keep],chain=factor(ch),panel=title)))
    }))
    outer<-"scope" %in% names(z) && any(z$scope=="outer_integration")
    conditional<-"scope" %in% names(z) && any(z$scope=="conditional_gp")
    caption<-if(outer)"All retained outer integration states; ESS includes finite inner Monte Carlo error. No ESS stopping requirement." else if(conditional)
      "Full retained conditional-GP history, with the outer measurement state fixed. ESS is descriptive; no stopping requirement." else
      "Display spans the full retained history; ESS uses every retained draw. Panel y-scales differ."
    p<-ggplot2::ggplot(dat,ggplot2::aes(x=iteration,y=value,color=chain))+
      ggplot2::geom_line(linewidth=.25,alpha=.7)+ggplot2::facet_wrap(~panel,ncol=2,scales="free")+
      ggplot2::theme_minimal(base_size=10)+ggplot2::theme(legend.position="bottom")+
      ggplot2::labs(title=paste("Scientific-target traces:",kind),x="Iteration (warmup excluded)",y="Target value",
        caption=caption)
    for(ext in formats) {
      dest<-file.path(output_dir,paste0("scientific_traces_",kind,".",ext))
      device<-if(ext=="png")function(filename,width,height,bg,...) {
        grDevices::png(filename,width=width,height=height,units="in",res=150,bg=bg,
          type=if(capabilities("cairo"))"cairo"else getOption("bitmapType"))
      }else grDevices::pdf
      ggplot2::ggsave(dest,p,device=device,width=13,height=7.5,dpi=150);paths<-c(paths,dest)
    }
  }
  invisible(paths)
}
