
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


# ---- applications/ocean_silicate/03_reporting.R ----
# 03: Combine compatible completed tasks; summarize and plot prediction scores.

# Capture console output, messages and warnings in a per-run log.
# Keep HTML focused on analysis while preserving complete execution logs.
ocean_logged <- function(path, expr) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  connection <- file(path, open = "wt")
  output_depth <- sink.number()
  sink(connection)
  on.exit({
    while (sink.number() > output_depth) sink()
    close(connection)
  }, add = TRUE)
  warnings <- character()
  value <- withCallingHandlers(eval(substitute(expr), envir = parent.frame()),
    message = function(m) { cat(conditionMessage(m)); invokeRestart("muffleMessage") },
    warning = function(w) {
      warnings <<- unique(c(warnings, conditionMessage(w)))
      cat("WARNING:", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    })
  attr(value, "notebook_warnings") <- warnings
  value
}


# Import only completed, compatible task payloads; validate before any copying.
# Independent machine roots are preferable to concurrent writers on shared storage.
ocean_combine_results <- function(source_dirs, output_dir, identity) {
  sources <- as.character(unlist(source_dirs)); sources <- sources[nzchar(sources)]
  plan <- list(); receipts <- list()
  for(source in sources) {
    source <- normalizePath(source,mustWork=TRUE)
    if(identical(source,normalizePath(output_dir,mustWork=TRUE))) stop("Cannot import a directory into itself.")
    path <- file.path(source,"workflow_identity.rds")
    if(!file.exists(path) || !identical(readRDS(path),identity)) stop("Incompatible imported run: ",source)
    done <- list.files(file.path(source,"workflow_tasks"),"[.]done[.]rds$",full.names=TRUE)
    for(receipt_path in done) {
      receipt <- readRDS(receipt_path); id <- receipt$task$task_id
      ocean_verify_receipt(source,id,identity)
      for(rel in receipt$files) {
        from <- file.path(source,rel); to <- file.path(output_dir,rel)
        key <- rel
        hash <- ocean_hashes(from)
        if(file.exists(to) && !identical(ocean_hashes(to),hash)) stop("Conflicting existing payload: ",rel)
        if(!is.null(plan[[key]]) && !identical(plan[[key]]$hash,hash)) stop("Conflicting imported payloads: ",rel)
        plan[[key]] <- list(from=from,to=to,hash=hash)
      }
      if(!is.null(receipts[[id]]) && !identical(receipts[[id]]$hashes,receipt$hashes))
        stop("Conflicting completed task: ",id)
      receipts[[id]] <- receipt
    }
  }
  for(item in plan) {
    dir.create(dirname(item$to),recursive=TRUE,showWarnings=FALSE)
    if(!file.exists(item$to)) stopifnot(file.copy(item$from,item$to))
    stopifnot(identical(ocean_hashes(item$to),item$hash))
  }
  for(id in names(receipts)) ocean_atomic_rds(receipts[[id]],
    file.path(output_dir,"workflow_tasks",paste0(id,".done.rds")))
  invisible(names(receipts))
}

# Read scores only; no fitting or test-score-based choice of placeholder.
ocean_bind <- function(parts) {
  if(!length(parts)) return(data.frame())
  cols <- unique(unlist(lapply(parts,names)))
  do.call(rbind,lapply(parts,function(x) {for(n in setdiff(cols,names(x))) x[[n]]<-NA; x[,cols,drop=FALSE]}))
}
ocean_prediction_report <- function(output_dir,inputs,profile="paper",augmented_output_dir="",report_dir=file.path(output_dir,"prediction_report")) {
  rows <- list(); inventory <- list(); original <- if(profile=="paper") "full5_class4" else "smoke_class4"
  final_sources <- ocean_final_eiv_sources(output_dir,profile)
  add <- function(path,f,cell,variant,b=NA_real_) {
    if(!file.exists(path)) return(invisible(NULL))
    z <- read.csv(path,stringsAsFactors=FALSE)
    stopifnot(nrow(z)>0L,all(z$fold==f),all(z$n_test==40L))
    if(variant=="EIV") {
      z$prediction_mode <- ifelse(z$method=="EIV-GP xc","class_only","mixed")
      if(cell=="O0") z <- z[z$prediction_mode=="class_only",,drop=FALSE]
      z$method <- "EIV-GP"
    } else if(variant=="Original baseline") z$prediction_mode <- "class_only"
    z$variant <- variant; z$placeholder <- b
    z$source <- gsub("\\\\","/",path)
    rows[[length(rows)+1L]] <<- z
  }
  for(f in if(profile=="paper") 1:5 else 1L) {
    initial <- file.path(output_dir,original,sprintf("fold%02d",f))
    add(file.path(initial,"baseline_metrics.csv"),f,"baseline","Original baseline")
    for(cell in if(profile=="paper") c("O0","O20","O50") else c("O0","O50")) {
      final<-final_sources[final_sources$fold==f & final_sources$cell==cell,]
      path<-file.path(output_dir,final$source)
      available<-isTRUE(final$available)
      inventory[[length(inventory)+1L]]<-data.frame(fold=f,cell,
        stage=if(final$continued)"plus20k" else "initial",available,
        initial_screen_pass=final$initial_screen_pass,final_screen_pass=final$final_screen_pass)
      # A cell that needs continuation never falls back to its initial fit.
      if(available)add(file.path(path,"predictive_metrics.csv"),f,cell,"EIV")
    }
  }
  aug <- if(nzchar(augmented_output_dir)) augmented_output_dir else file.path(output_dir,"augmented_baselines")
  if(profile=="paper") {
    tasks <- ocean_task_catalogue("paper"); tasks <- tasks[tasks$module=="augmented_baselines",]
    for(i in seq_len(nrow(tasks))) {
      task <- tasks[i,]; add(file.path(aug,"fits",task$task_id,"metrics.csv"),
        task$fold,task$cell,"Expanded baseline",task$placeholder)
    }
  }
  raw <- ocean_bind(rows); summaries <- list()
  metrics <- c("RMSE","MAE","CRPS","NLPD","Coverage95","Width95","IntervalScore95")
  if(nrow(raw)) {
    key <- paste(raw$variant,raw$method,raw$cell,raw$prediction_mode,raw$placeholder,sep="/")
    if(anyDuplicated(paste(key,raw$fold))) stop("Duplicate fold/task scores in report.")
    for(group in split(seq_len(nrow(raw)),key)) {
      z <- raw[group,,drop=FALSE]
      ans <- z[1,c("variant","method","cell","prediction_mode","placeholder"),drop=FALSE]
      ans$n_folds <- nrow(z); ans$complete_5fold <- nrow(z)==5L && setequal(z$fold,1:5)
      for(metric in metrics) {
        stopifnot(all(is.finite(z[[metric]])))
        ans[[paste0(metric,"_mean")]] <- mean(z[[metric]])
        ans[[paste0(metric,"_sd")]] <- if(nrow(z)>1L) sd(z[[metric]]) else NA_real_
      }
      summaries[[length(summaries)+1L]] <- ans
    }
  }
  summary_raw <- ocean_bind(summaries); scaled <- raw
  if(nrow(scaled)) {
    for(m in c("RMSE","MAE","CRPS","Width95","IntervalScore95")) scaled[[m]] <- scaled[[m]]/inputs$response_sd
    scaled$NLPD <- scaled$NLPD-log(inputs$response_sd)
  }
  # Presentation is restricted to the prespecified primary encoding b=0.
  # Keep all saved scores in raw audit exports; do not select b from test scores.
  plots <- list()
  if(nrow(summary_raw) && requireNamespace("ggplot2",quietly=TRUE)) {
    d <- summary_raw[(summary_raw$variant=="Expanded baseline" &
      !is.na(summary_raw$placeholder) & summary_raw$placeholder==0 |
      summary_raw$variant=="EIV") & summary_raw$cell %in% c("O20","O50") &
      summary_raw$complete_5fold,,drop=FALSE]
    if(nrow(d)) {
      d$method <- factor(d$method,levels=c("UC-GP","LVGP","EzGP","EIV-GP"))
      d$calibration <- factor(d$cell,levels=c("O20","O50"),labels=c("20% exact training oxygen","50% exact training oxygen"))
      d$test_input <- factor(d$prediction_mode,levels=c("class_only","mixed"),
        labels=c("Class-only test inputs","Matched mixed test inputs"))
      plots$primary_comparison <- ggplot2::ggplot(d,ggplot2::aes(RMSE_mean,method,color=method)) +
        ggplot2::geom_linerange(ggplot2::aes(xmin=RMSE_mean-RMSE_sd,xmax=RMSE_mean+RMSE_sd),
          orientation="y",linewidth=.65) + ggplot2::geom_point(size=2.8) +
        ggplot2::facet_grid(test_input~calibration) +
        ggplot2::scale_color_manual(values=c("UC-GP"="#7b3294",LVGP="#008837",EzGP="#a6761d","EIV-GP"="#2365a3"),guide="none") +
        ggplot2::labs(x="Silicate prediction RMSE (micromol/kg)",y=NULL,
          title="Prediction with matched available oxygen measurements",
          caption="Expanded baselines use the prespecified b=0 encoding. Points: five-fold means; bars: sample SD.\nMixed test inputs reveal exact oxygen for 20%/50% of test profiles, matching training calibration.") +
        ggplot2::theme_minimal(base_size=11) +
        ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),panel.grid.major.y=ggplot2::element_blank(),
          plot.caption=ggplot2::element_text(hjust=0))
    }
  }
  report <- list(fold_metrics=raw,summary=summary_raw,availability=ocean_bind(inventory),
    SD_Y=inputs$response_sd,plots=plots,fold_metrics_SD_Y=scaled)
  if(!is.null(report_dir)) ocean_export_prediction_report(report,report_dir)
  report

}
# Shared table formatting for CSV/Markdown/LaTeX/HTML. Formatting does not
# change the fold scores or their aggregation; paper coverage is a percentage.
ocean_mean_sd <- function(mean,sd,digits=2L) {
  ifelse(is.na(sd),sprintf(paste0("%.",digits,"f (SD unavailable)"),mean),
    sprintf(paste0("%.",digits,"f (%.",digits,"f)"),mean,sd))
}
ocean_table_align <- function(d) {
  a <- rep("r",ncol(d)); k <- if("Method" %in% names(d)) 3L else if("Class" %in% names(d)) 2L else 1L
  a[seq_len(min(k,ncol(d)))] <- "l"; a
}
ocean_tex_note <- function(note) gsub("%","\\%",note,fixed=TRUE)
# Rank unrounded means within comparable groups. Incomplete groups are not ranked.
ocean_best_mask <- function(values, groups, eligible = rep(TRUE, nrow(values))) {
  best <- matrix(FALSE, nrow(values), ncol(values), dimnames=list(NULL,names(values)))
  for(rows in split(seq_len(nrow(values)), groups)) {
    if(length(rows)<2L || !isTRUE(all(eligible[rows]))) next
    for(j in seq_len(ncol(values))) {
      v <- values[rows,j]
      if(all(is.finite(v)))
        best[rows,j] <- abs(v-min(v)) <= 1e-12*max(1,abs(min(v)))
    }
  }
  best
}
ocean_html_escape <- function(x) {
  x <- gsub("&","&amp;",as.character(x),fixed=TRUE)
  x <- gsub("<","&lt;",x,fixed=TRUE)
  gsub(">","&gt;",x,fixed=TRUE)
}
# Formatting attributes never enter CSV scores; only the mean is emphasized.
ocean_table_markup <- function(d, format=c("html","pipe","latex")) {
  format <- match.arg(format)
  mask <- attr(d,"best")
  out <- d
  for(j in seq_len(ncol(d))) {
    out[[j]] <- if(format=="html") ocean_html_escape(d[[j]]) else as.character(d[[j]])
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
ocean_table_html <- function(d,caption,note="") {
  table <- as.character(knitr::kable(ocean_table_markup(d,"html"),format="html",caption=caption,
    escape=FALSE,row.names=FALSE,align=ocean_table_align(d),
    table.attr='class="ocean-summary table" style="width:100%;border-collapse:collapse"'))
  table <- sub("<caption>",'<caption style="text-align:left;padding:8px 0;color:#555">',table,fixed=TRUE)
  table <- sub("<thead>",'<thead style="background:#f2f4f6">',table,fixed=TRUE)
  table <- gsub('style="text-align:',
    'style="padding:8px 10px;white-space:nowrap;border-bottom:1px solid #ddd;text-align:',
    table,fixed=TRUE)
  paste0('<div class="ocean-table-scroll" style="max-width:100%;overflow-x:auto">',table,
    '</div><p class="note">',ocean_html_escape(note),"</p>")
}
# Training and test information define each comparison panel, as in the paper.
ocean_prediction_panel_html <- function(d, caption, note="") {
  if(!nrow(d)) return("")
  columns <- setdiff(names(d),c("Training exact U","Test oxygen"))
  styled <- ocean_table_markup(d,"html")
  training_groups <- unique(d[["Training exact U"]])
  out <- c("<table class='ocean-prediction table' style='width:100%;border-collapse:collapse'>",
    paste0("<caption style='text-align:left'>",ocean_html_escape(caption),"</caption>"),
    "<thead><tr>",paste0("<th style='text-align:",ifelse(columns=="Method","left","right"),
      ";padding:7px 10px'>",ocean_html_escape(columns),"</th>"),"</tr></thead><tbody>")
  for(j in seq_along(training_groups)) {
    training_rows <- which(d[["Training exact U"]]==training_groups[j])
    tests <- unique(d[["Test oxygen"]][training_rows])
    title <- paste0(LETTERS[j],". ",training_groups[j]," exact training oxygen")
    if(length(tests)==1L) title <- paste0(title,"; test: ",tests)
    out <- c(out,paste0("<tr class='comparison-panel'><th scope='rowgroup' colspan='",length(columns),
      "' style='text-align:left;background:#f2f4f6;border-top:2px solid #bbb;padding:9px 10px'>",
      ocean_html_escape(title),"</th></tr>"))
    for(test in tests) {
      rows <- training_rows[d[["Test oxygen"]][training_rows]==test]
      if(length(tests)>1L) out <- c(out,paste0(
        "<tr class='comparison-condition'><th scope='rowgroup' colspan='",length(columns),
        "' style='text-align:left;font-weight:500;padding:7px 10px;border-bottom:1px solid #ddd'>",
        "Test oxygen: ",ocean_html_escape(test),"</th></tr>"))
      for(i in rows) out <- c(out,"<tr>",paste0("<td style='text-align:",
        ifelse(columns=="Method","left","right"),";padding:6px 10px'>",
        unlist(styled[i,columns,drop=FALSE],use.names=FALSE),"</td>"),"</tr>")
    }
  }
  paste(c(out,"</tbody></table>",paste0("<p class='note'>",ocean_html_escape(note),"</p>")),collapse="\n")
}
ocean_export_paper_table <- function(d,name,caption,dest,note="") {
  if(!nrow(d)) return(invisible(NULL))
  ocean_atomic_csv(d,file.path(dest,paste0(name,".csv")))
  ocean_atomic_lines(c(as.character(knitr::kable(ocean_table_markup(d,"pipe"),format="pipe",row.names=FALSE,align=ocean_table_align(d))),"",note),file.path(dest,paste0(name,".md")))
  tex_data <- ocean_table_markup(d,"latex")
  tex_data[] <- lapply(tex_data,ocean_tex_note)
  names(tex_data) <- ocean_tex_note(names(tex_data))
  tex <- as.character(knitr::kable(tex_data,format="latex",booktabs=TRUE,caption=caption,label=name,
    escape=FALSE,row.names=FALSE,align=ocean_table_align(d)))
  tex <- paste0(strsplit(tex,"\\end{table}",fixed=TRUE)[[1]][1],
    "\\par\\smallskip\n\\begin{minipage}{\\linewidth}\\footnotesize\n",ocean_tex_note(note),
    "\n\\end{minipage}\n\\end{table}")
  ocean_atomic_lines(tex,file.path(dest,paste0(name,".tex")))
  html <- ocean_table_html(d,caption,note)
  ocean_atomic_lines(c("<!doctype html><meta charset='utf-8'>",
    "<style>body{font:15px/1.5 sans-serif;margin:24px}table{border-collapse:collapse}td,th{padding:7px 10px;border-bottom:1px solid #ddd}caption{text-align:left;font-weight:600}.note{color:#555}</style>",html),
    file.path(dest,paste0(name,".html")))
  invisible(html)
}
# Group training/test information in panel headings, keeping the paper LaTeX
# body compact. CSV retains explicit columns for filtering and checking.
ocean_prediction_panel_tex <- function(d,columns,headers,label,caption,note) {
  styled <- ocean_table_markup(d,"latex")
  tex <- c("\\begin{table}[t]","\\centering\\small",paste0("\\caption{",caption,"}"),
    paste0("\\label{tab:",label,"}"),"\\setlength{\\tabcolsep}{4pt}","\\renewcommand{\\arraystretch}{1.1}",
    paste0("\\begin{tabular}{l",paste(rep("r",length(columns)),collapse=""),"}"),
    "\\toprule",paste0("Method & ",paste(headers,collapse=" & "),"\\\\"),"\\midrule")
  training_groups <- unique(d[["Training exact U"]])
  for(j in seq_along(training_groups)) {
    training_rows <- which(d[["Training exact U"]]==training_groups[j])
    tests <- unique(d[["Test oxygen"]][training_rows])
    if(j>1) tex <- c(tex,"\\addlinespace[4pt]")
    title <- paste0(LETTERS[j],". ",training_groups[j]," exact training oxygen")
    if(length(tests)==1L) title <- paste0(title,"; test: ",tests)
    tex <- c(tex,paste0("\\multicolumn{",length(columns)+1L,"}{l}{\\textit{",ocean_tex_note(title),"}}\\\\[2pt]"))
    for(test in tests) {
      rows <- training_rows[d[["Test oxygen"]][training_rows]==test]
      if(length(tests)>1L) tex <- c(tex,paste0("\\multicolumn{",length(columns)+1L,
        "}{l}{Test oxygen: ",ocean_tex_note(test),"}\\\\[2pt]"))
      for(i in rows) tex <- c(tex,paste0(sub("*","$^{*}$",d$Method[i],fixed=TRUE)," & ",
        paste(unlist(styled[i,columns,drop=FALSE],use.names=FALSE),collapse=" & "),"\\\\"))
    }
  }
  c(tex,"\\bottomrule","\\end{tabular}","\\par\\smallskip","\\begin{minipage}{\\linewidth}\\footnotesize",
    ocean_tex_note(note),"\\end{minipage}","\\end{table}")
}
ocean_paper_prediction_tables <- function(report) {
  z <- report$summary
  if(!nrow(z)) return(list(main=data.frame(),full=data.frame(),numeric=data.frame()))
  z <- z[z$variant!="Expanded baseline" | (!is.na(z$placeholder) & z$placeholder==0),,drop=FALSE]
  z$training_cell <- ifelse(z$variant=="Original baseline","O0",z$cell)
  z <- z[order(match(z$training_cell,c("O0","O20","O50")),match(z$prediction_mode,c("class_only","mixed")),
    match(z$method,c("UC-GP","LVGP","EzGP","EIV-GP"))),,drop=FALSE]
  base <- data.frame("Training exact U"=paste0(sub("O","",z$training_cell),"%"),
    "Test oxygen"=ifelse(z$prediction_mode=="class_only","Class only",paste0(sub("O","",z$training_cell),"% exact")),
    Method=ifelse(z$variant=="Expanded baseline",paste0(z$method,"*"),z$method),check.names=FALSE)
  if(any(!z$complete_5fold)) base$Method <- ifelse(z$complete_5fold,base$Method,
    paste0(base$Method," [n=",z$n_folds," folds]"))
  format_metric <- function(metric,multiplier=1,digits=2,offset=0)
    ocean_mean_sd(multiplier*z[[paste0(metric,"_mean")]]-offset,multiplier*z[[paste0(metric,"_sd")]],digits)
  sy <- report$SD_Y
  main <- base; main$RMSE <- format_metric("RMSE")
  main[["CRPS / sY"]] <- format_metric("CRPS",1/sy,3)
  main[["Coverage (%)"]] <- format_metric("Coverage95",100,1)
  main[["Width / sY"]] <- format_metric("Width95",1/sy,3)
  full <- base; full$RMSE <- format_metric("RMSE")
  full[["RMSE / sY"]] <- format_metric("RMSE",1/sy,3)
  full[["MAE / sY"]] <- format_metric("MAE",1/sy,3)
  full[["CRPS / sY"]] <- format_metric("CRPS",1/sy,3)
  full[["NLPD standardized"]] <- format_metric("NLPD",1,3,log(sy))
  full[["Coverage (%)"]] <- format_metric("Coverage95",100,1)
  full[["Width / sY"]] <- format_metric("Width95",1/sy,3)
  full[["Interval score / sY"]] <- format_metric("IntervalScore95",1/sy,3)
  groups <- paste(base[["Training exact U"]],base[["Test oxygen"]],sep="/")
  main_values <- data.frame(RMSE=z$RMSE_mean,"CRPS / sY"=z$CRPS_mean/sy,check.names=FALSE)
  attr(main,"best") <- ocean_best_mask(main_values,groups,z$complete_5fold)
  full_values <- data.frame(RMSE=z$RMSE_mean,"RMSE / sY"=z$RMSE_mean/sy,
    "MAE / sY"=z$MAE_mean/sy,"CRPS / sY"=z$CRPS_mean/sy,
    "NLPD standardized"=z$NLPD_mean-log(sy),
    "Interval score / sY"=z$IntervalScore95_mean/sy,check.names=FALSE)
  attr(full,"best") <- ocean_best_mask(full_values,groups,z$complete_5fold)
  list(main=main,full=full,numeric=z)
}
ocean_write_prediction_tables <- function(tables,dest,sy) {
  if(!nrow(tables$main)) return(invisible(NULL))
  note <- paste0("Entries are fold means (sample SDs). sY = ",sprintf("%.4f",sy),
    " micromol/kg, the frozen cohort response SD. RMSE is in micromol/kg; coverage is percent. ",
    "* marks baselines with the same recorded oxygen information as EIV-GP: retain (x,c), add categorical observation status R and continuous oxygen record W, and use b=0 for hidden W. Native kernels and fitting procedures are used. ",
    "At 0%, they use their original (x,c) inputs. Test inputs and held-out profiles are matched within each block. ",
    "Bold marks the smallest unrounded mean error or predictive score within a complete five-fold block; SDs are not bold. Coverage and width are descriptive.")
  ocean_export_paper_table(tables$main,"paper_y_prediction_main","Silicate prediction with matched oxygen information",dest,note)
  ocean_export_paper_table(tables$full,"paper_y_prediction_full","Additional silicate prediction summaries: b=0",dest,note)
  for(name in c("main","full")) {
    content <- ocean_prediction_panel_html(tables[[name]],"Silicate prediction with matched oxygen information",note)
    ocean_atomic_lines(c("<!doctype html><meta charset='utf-8'>",
      "<style>body{font:15px/1.5 sans-serif;margin:24px}.note{color:#555}</style>",content),
      file.path(dest,paste0("paper_y_prediction_",name,".html")))
  }
  ocean_atomic_lines(ocean_prediction_panel_tex(tables$main,
    c("RMSE","CRPS / sY","Coverage (%)","Width / sY"),
    c("RMSE","CRPS/$s_Y$","Cov. (\\%)","Width/$s_Y$"),"paper_y_prediction_main",
    "Silicate prediction under matched oxygen information",note),file.path(dest,"paper_y_prediction_main.tex"))
  ocean_atomic_lines(c(ocean_prediction_panel_tex(tables$full,
    c("RMSE","RMSE / sY","MAE / sY"),c("RMSE","RMSE/$s_Y$","MAE/$s_Y$"),
    "paper_y_prediction_full_errors","Point prediction summaries: primary encoding b=0",note),"",
    ocean_prediction_panel_tex(tables$full,
    c("CRPS / sY","NLPD standardized","Coverage (%)","Width / sY","Interval score / sY"),
    c("CRPS/$s_Y$","NLPD (std.)","Cov. (\\%)","Width/$s_Y$","IS/$s_Y$"),
    "paper_y_prediction_full_distribution","Predictive distribution summaries: primary encoding b=0",note)),
    file.path(dest,"paper_y_prediction_full.tex"))
  ocean_atomic_csv(tables$numeric,file.path(dest,"primary_summary_numeric.csv"))
}
ocean_export_prediction_report <- function(report,dest) {
  ocean_atomic_csv(report$fold_metrics,file.path(dest,"fold_metrics_raw.csv"))
  ocean_atomic_csv(report$summary,file.path(dest,"summary_raw.csv"))
  ocean_atomic_csv(report$fold_metrics_SD_Y,file.path(dest,"fold_metrics_SD_Y.csv"))
  ocean_atomic_csv(report$availability,file.path(dest,"eiv_availability.csv"))
  for(name in names(report$plots)) for(ext in c("png","pdf")) ggplot2::ggsave(
    file.path(dest,paste0(name,".",ext)),report$plots[[name]],width=10,height=6,dpi=180,bg="white")
  ocean_write_prediction_tables(ocean_paper_prediction_tables(report),dest,report$SD_Y)
  invisible(report)
}
ocean_show_predictions <- function(report, show_plot = FALSE) {
  if(is.null(report) || !nrow(report$summary)) {cat("No completed prediction scores yet.\n");return(invisible(NULL))}
  tables <- ocean_paper_prediction_tables(report)
  cat("\n\n### Silicate prediction\n\n",
    "Entries are fold means (sample SDs). RMSE is in micromol/kg; sY is the frozen cohort response SD. ",
    "* marks UC-GP, LVGP and EzGP with the same recorded training/test oxygen measurements as EIV-GP. ",
    "Incomplete groups state their fold count in the method label.\n\n",sep="")
  cat(ocean_prediction_panel_html(tables$main,"Silicate prediction under matched oxygen information",
    "Bold marks the smallest RMSE and CRPS mean within each complete five-fold information block, using unrounded values. Coverage and width are descriptive."),"\n\n")
  if(show_plot) for(p in report$plots) print(p)
  invisible(report)
}

# One row per fold combines fit availability, continuation and final screening.
ocean_fit_overview <- function(analysis, output_dir) {
  settings <- analysis$settings
  status <- ocean_final_eiv_sources(output_dir, settings$profile)
  overview <- data.frame(Fold = settings$folds, check.names = FALSE)
  for(cell in settings$cells) {
    overview[[paste0(sub("O", "", cell), "% exact U")]] <- vapply(settings$folds, function(fold) {
      row <- status[status$fold == fold & status$cell == cell, , drop = FALSE]
      if(settings$action == "plan") return("Planned")
      if(!nrow(row) || !isTRUE(row$available))
        return(if(nrow(row) && isTRUE(row$continued)) "Continuation pending" else "Pending")
      if(settings$profile == "quick") return("Initial saved")
      stage <- if(isTRUE(row$continued)) "+20k" else "Initial"
      paste(stage, if(isTRUE(row$final_screen_pass)) "pass" else "review", sep = "; ")
    }, "")
  }
  cat(as.character(knitr::kable(overview, format = "pipe", row.names = FALSE)), sep = "\n")
  cat("\n\n")
  if(settings$profile == "paper" && settings$action != "plan") {
    cat("Initial / +20k identifies the saved fit used in the report. Pass refers to the prespecified finite-Rhat/bulk-ESS screen; full chain diagnostics remain in the result files.\n\n")
    available <- status$available & !is.na(status$final_screen_pass)
    if(any(available & !status$final_screen_pass, na.rm = TRUE))
      cat("Fits marked review remain outside this limited screen.\n\n")
  }
  invisible(overview)
}


# ---- applications/ocean_silicate/04_diagnostics.R ----
# 04: U imputation, category means, continuous functions, and prediction diagnostics.
# Loading defines functions only. Direct Rscript execution uses a parameter RDS.

# Final paper-profile diagnostics. Source this file after 01_setup.R; sourcing
# only defines functions. It never fits models, changes masks, or installs packages.

# 1. Frozen inputs and final-fit provenance ------------------------------------
ocean_diag_require <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

ocean_diag_csv <- function(path, columns) {
  ocean_diag_require(file.exists(path), paste("Missing input:", path))
  x <- utils::read.csv(path, stringsAsFactors = FALSE)
  ocean_diag_require(all(columns %in% names(x)), paste("Unexpected columns:", path))
  x
}

ocean_diag_flag <- function(x) {
  ocean_diag_require(!anyNA(x) && all(x %in% c(0, 1, FALSE, TRUE)),
    "Observation masks must contain only 0/1 or FALSE/TRUE.")
  as.logical(x)
}

ocean_diag_inputs <- function(data_dir) {
  columns <- c("cohort_id", "profile", "G2salinity", "G2pressure", "log_pressure", "U", "Y")
  cohort <- ocean_diag_csv(file.path(data_dir, "cohort_manifest.csv"), columns)
  ocean_diag_require(nrow(cohort) == 200L && !anyDuplicated(cohort$cohort_id) &&
    !anyDuplicated(cohort$profile), "Expected the frozen 200-profile cohort.")
  check_values <- function(d) {
    ocean_diag_require(all(is.finite(as.matrix(d[, columns[3:7]]))) &&
      all(d$U > 0 & d$U < 180) && all(d$G2pressure >= 300 & d$G2pressure <= 1500) &&
      max(abs(d$log_pressure - log(d$G2pressure))) < 1e-10,
      "Inputs must use physical units and natural-log pressure in the frozen domain.")
    d$C <- findInterval(d$U, c(50, 100, 150)) + 1L # Derive the four oxygen classes from measured oxygen.
    d
  }
  cohort <- check_values(cohort)
  ocean_diag_require(identical(as.integer(table(factor(cohort$C, levels = 1:4))), rep(50L, 4)),
    "Expected 50 reference observations in each recorded oxygen band.")
  train <- test <- audit <- vector("list", 5L)
  for (fold in 1:5) {
    tr <- check_values(ocean_diag_csv(file.path(data_dir, sprintf("fold%02d_train.csv", fold)),
      c(columns, "calib_o20", "calib_o50")))
    te <- check_values(ocean_diag_csv(file.path(data_dir, sprintf("fold%02d_test.csv", fold)), columns))
    tr$calib_o20 <- ocean_diag_flag(tr$calib_o20)
    tr$calib_o50 <- ocean_diag_flag(tr$calib_o50)
    ocean_diag_require(nrow(tr) == 160L && nrow(te) == 40L &&
      !anyDuplicated(c(tr$cohort_id, te$cohort_id)) &&
      setequal(c(tr$cohort_id, te$cohort_id), cohort$cohort_id) &&
      !length(intersect(tr$profile, te$profile)), "Train/test separation is inconsistent.")
    for (d in list(tr, te)) {
      ref <- cohort[match(d$cohort_id, cohort$cohort_id), columns]
      ocean_diag_require(isTRUE(all.equal(d[, columns], ref, check.attributes = FALSE)),
        "Fold data differ from the frozen cohort.")
    }
    ocean_diag_require(all(tr$calib_o20 <= tr$calib_o50) &&
      all(table(factor(tr$C[tr$calib_o20], levels = 1:4)) == 8L) &&
      all(table(factor(tr$C[tr$calib_o50], levels = 1:4)) == 20L) &&
      all(table(factor(tr$C, levels = 1:4)) == 40L) &&
      all(table(factor(te$C, levels = 1:4)) == 10L), "Calibration or class counts changed.")
    train[[fold]] <- tr; test[[fold]] <- te
    audit[[fold]] <- data.frame(fold = fold, n_train = nrow(tr), n_test = nrow(te),
      n_calib_O20 = sum(tr$calib_o20), n_calib_O50 = sum(tr$calib_o50),
      n_common_missing = sum(!tr$calib_o50))
  }
  ocean_diag_require(!anyDuplicated(unlist(lapply(test, `[[`, "cohort_id"))),
    "Each cohort observation must be held out exactly once across five folds.")
  list(cohort = cohort, train = train, test = test, audit = do.call(rbind, audit))
}

ocean_diag_sources <- function(output_dir) {
  selected<-ocean_final_eiv_sources(output_dir,"paper")
  ocean_diag_require(all(selected$available),"Final fits selected from current initial diagnostics are incomplete.")
  sources<-selected[,c("fold","cell","source","continued","fit_file")]
  for(i in seq_len(nrow(sources))) {
    s<-sources[i,];path<-file.path(output_dir,s$source)
    for(file in c(s$fit_file,"predictive_metrics.csv","prediction_distributions.rds",
      "U_posterior_intervals.csv","convergence_summary.csv"))
      ocean_diag_require(file.exists(file.path(path,file)),paste("Incomplete paper output:",path,file))
  }
  manifest <- file.path(output_dir, "continuation_plus20k", "final_review_20260927", "fit_sources.csv")
  if (file.exists(manifest)) {
    recorded <- ocean_diag_csv(manifest, c("fold", "cell", "source", "continued"))
    i <- match(paste(sources$fold, sources$cell), paste(recorded$fold, recorded$cell))
    ocean_diag_require(nrow(recorded) == 15L && !anyNA(i) &&
      identical(sources$source, recorded$source[i]), "Final-review source manifest disagrees with the paper run.")
  }
  sources
}

# 2. Training U imputation: five metrics on the same 80 hidden rows ------------
ocean_diag_five_metrics <- function(d) {
  error <- d$posterior_mean - d$reference_U
  data.frame(n = nrow(d), Bias = mean(error), RMSE = sqrt(mean(error^2)),
    MAE = mean(abs(error)), Coverage95 = mean(d$reference_U >= d$q025 & d$reference_U <= d$q975),
    Width95 = mean(d$q975 - d$q025))
}

ocean_diag_fold_summary <- function(d, groups, metrics) {
  # Fold SD is descriptive variability across overlapping CV training sets.
  # It is neither posterior SD nor an independent-replicate Monte Carlo SE.
  key <- do.call(paste, c(d[groups], sep = "/"))
  do.call(rbind, lapply(split(seq_len(nrow(d)), key), function(i) {
    z <- d[i, , drop = FALSE]
    ocean_diag_require(nrow(z) == 5L && setequal(z$fold, 1:5), "Expected one score per fold.")
    do.call(rbind, lapply(metrics, function(metric) {
      ocean_diag_require(all(is.finite(z[[metric]])), paste("Non-finite score:", metric))
      data.frame(z[1, groups, drop = FALSE], metric = metric, mean = mean(z[[metric]]),
        sd = stats::sd(z[[metric]]), n_folds = nrow(z), row.names = NULL)
    }))
  }))
}

# Display all 160 training rows in each panel. This is deliberately separate
# from the common-hidden scoring set below: calibration values are NOT imputations.
ocean_diag_u_display <- function(output_dir, inputs, sources) {
  parts <- scales <- counts <- list()
  for (i in seq_len(nrow(sources))) {
    s <- sources[i, ]; tr <- inputs$train[[s$fold]]
    observed <- switch(s$cell, O0 = rep(FALSE, nrow(tr)),
      O20 = tr$calib_o20, O50 = tr$calib_o50)
    z <- ocean_diag_csv(file.path(output_dir, s$source, "U_posterior_intervals.csv"),
      c("cohort_id", "C", "posterior_mean", "q025", "q975", "scale"))
    expected_scale <- if (s$cell == "O0") "model" else "raw"
    ocean_diag_require(!anyDuplicated(z$cohort_id) &&
      setequal(z$cohort_id, tr$cohort_id[!observed]) && all(z$scale == expected_scale),
      "U display rows or scales disagree with the fitted observation mask.")
    # The full TRAINING reference mean/SD are an external plotting convention.
    # Hidden U is used here only for evaluation/axes, never to fit/align O0.
    center <- mean(tr$U); scale <- stats::sd(tr$U)
    ocean_diag_require(is.finite(scale) && scale > 0, "Invalid training reference SD.")
    d <- data.frame(fold = s$fold, cell = s$cell, cohort_id = tr$cohort_id,
      C = tr$C, observed = observed, common_hidden = !tr$calib_o50,
      reference_U = tr$U, posterior_mean = NA_real_, q025 = NA_real_, q975 = NA_real_,
      posterior_scale = expected_scale, reference_center = center, reference_sd = scale,
      plot_reference = (tr$U - center) / scale,
      plot_estimate = NA_real_, plot_q025 = NA_real_, plot_q975 = NA_real_)
    j <- match(z$cohort_id, tr$cohort_id)
    ocean_diag_require(!anyNA(j) && all(z$C == tr$C[j]) &&
      all(is.finite(as.matrix(z[, c("posterior_mean", "q025", "q975")]))) &&
      all(z$q975 > z$q025), "Invalid saved U posterior display intervals.")
    d$posterior_mean[j] <- z$posterior_mean
    d$q025[j] <- z$q025; d$q975[j] <- z$q975
    # O0 is left in its ORIGINAL model coordinate. Standardizing its posterior
    # or fitting a truth-based affine alignment would conceal scale ambiguity.
    shift <- if (s$cell == "O0") 0 else center
    divisor <- if (s$cell == "O0") 1 else scale
    d$plot_estimate[j] <- (z$posterior_mean - shift) / divisor
    d$plot_q025[j] <- (z$q025 - shift) / divisor
    d$plot_q975[j] <- (z$q975 - shift) / divisor
    d$plot_estimate[observed] <- d$plot_reference[observed]
    # NA posterior summaries on observed rows prevent them being mistaken for
    # inferred values or entering imputation scores. Only display values coincide.
    parts[[i]] <- d
    scales[[i]] <- data.frame(fold = s$fold, cell = s$cell,
      reference_center = center, reference_sd = scale, posterior_scale = expected_scale,
      interpretation = if (s$cell == "O0") "Unanchored model coordinate; illustrative only" else
        "Raw posterior transformed to training-reference SD units")
    counts[[i]] <- data.frame(fold = s$fold, cell = s$cell, n_train = nrow(tr),
      n_masked = sum(!observed), n_observed = sum(observed), n_common_hidden_scored = if (s$cell == "O0") 0L else 80L)
  }
  list(points = do.call(rbind, parts), scales = do.call(rbind, scales), counts = do.call(rbind, counts))
}

ocean_diag_u <- function(output_dir, inputs, sources) {
  points <- list()
  for (i in which(sources$cell != "O0")) {
    s <- sources[i, ]; tr <- inputs$train[[s$fold]]
    z <- ocean_diag_csv(file.path(output_dir, s$source, "U_posterior_intervals.csv"),
      c("cohort_id", "C", "posterior_mean", "q025", "q975", "scale"))
    mask <- tr[[if (s$cell == "O20") "calib_o20" else "calib_o50"]]
    ocean_diag_require(!anyDuplicated(z$cohort_id) && setequal(z$cohort_id, tr$cohort_id[!mask]) &&
      all(z$scale == "raw"), "Saved U intervals do not match the training observation mask.")
    # Nested masks imply that !calib_o50 is hidden under BOTH fitted models.
    ref <- tr[!tr$calib_o50, ]; z <- z[match(ref$cohort_id, z$cohort_id), ]
    ocean_diag_require(nrow(z) == 80L && all(z$C == ref$C) &&
      all(is.finite(as.matrix(z[, c("posterior_mean", "q025", "q975")]))) &&
      all(z$q975 > z$q025), "Invalid common-missing U intervals.")
    points[[length(points) + 1L]] <- data.frame(fold = s$fold, cell = s$cell,
      cohort_id = ref$cohort_id, C = ref$C, reference_U = ref$U,
      posterior_mean = z$posterior_mean, q025 = z$q025, q975 = z$q975)
  }
  points <- do.call(rbind, points)
  scores <- list()
  for (fold in 1:5) for (cell in c("O20", "O50")) for (C in 0:4) {
    z <- points[points$fold == fold & points$cell == cell & (C == 0 | points$C == C), ]
    scores[[length(scores) + 1L]] <- data.frame(fold = fold, cell = cell,
      group = if (C == 0) "overall" else paste0("class", C), ocean_diag_five_metrics(z))
  }
  scores <- do.call(rbind, scores)
  display <- ocean_diag_u_display(output_dir, inputs, sources)
  list(points = points, fold_metrics = scores, summary = ocean_diag_fold_summary(scores,
    c("cell", "group"), c("Bias", "RMSE", "MAE", "Coverage95", "Width95")),
    display_points = display$points, display_scales = display$scales, display_counts = display$counts)
}

# 3. Y-blind display domains for m(x,c) and f(x,u) ------------------------------
ocean_diag_grids <- function(cohort) {
  # Fixed S=34.20. Within each recorded band, select the 25
  # profiles closest in salinity; use their central 80% log-pressure range.
  # Ties are resolved by frozen cohort ID, without looking at Y.
  parts <- neighbors <- list()
  for (C in 1:4) {
    d <- cohort[cohort$C == C, ]
    near <- d[order(abs(d$G2salinity - 34.20), d$cohort_id)[1:25], ]
    bounds <- stats::quantile(near$log_pressure, c(.1, .9), type = 7)
    lp <- seq(bounds[1], bounds[2], length.out = 9L)
    parts[[C]] <- data.frame(design = "candidate_A", C = C,
      G2salinity = 34.20, log_pressure = lp, pressure = exp(lp))
    neighbors[[C]] <- data.frame(design = "candidate_A", class = C, cohort_id = near$cohort_id)
  }
  m <- do.call(rbind, parts)
  # Selected fixed X: S=34.25, P=650. Reference U is used ONLY to set an
  # offline illustration domain and rug; it is never supplied as calibration.
  center <- c(34.25, log(650))
  x <- as.matrix(cohort[, c("G2salinity", "log_pressure")])
  distance <- sqrt(rowSums(sweep(sweep(x, 2, center, "-"), 2, apply(x, 2, stats::sd), "/")^2))
  near <- cohort[order(distance, cohort$cohort_id)[1:30], ]
  bounds <- stats::quantile(near$U, c(.1, .9), type = 7)
  f <- data.frame(G2salinity = center[1], log_pressure = center[2], pressure = 650,
    U = seq(bounds[1], bounds[2], length.out = 41L))
  same_x <- data.frame(design = "same_x", C = 1:4, G2salinity = center[1],
    log_pressure = center[2], pressure = 650)
  m <- rbind(m, same_x); m$grid_id <- sprintf("m_%02d", seq_len(nrow(m)))
  f$grid_id <- sprintf("f_%02d", seq_len(nrow(f)))
  list(m = m, f = f, near = near, m_neighbors = do.call(rbind, neighbors),
    support = data.frame(S = 34.25, pressure = 650, n_neighbors = 30L,
      U_low = bounds[1], U_high = bounds[2],
      C1 = sum(near$C == 1L), C2 = sum(near$C == 2L), C3 = sum(near$C == 3L), C4 = sum(near$C == 4L)))
}

ocean_diag_balanced_ids <- function(fit, n_per_chain) {
  lens <- vapply(fit$mcmc$samples_by_chain$u, nrow, integer(1))
  ocean_diag_require(length(lens) == 4L && all(lens >= n_per_chain),
    "Final diagnostics require four retained paper chains.")
  offsets <- c(0L, head(cumsum(lens), -1L))
  unlist(Map(function(n, off) off + unique(round(seq(1, n, length.out = n_per_chain))), lens, offsets),
    use.names = FALSE)
}

ocean_diag_point_summary <- function(draws, means, grid) {
  ocean_diag_require(is.matrix(draws) && is.matrix(means) && ncol(draws) == nrow(grid) &&
    ncol(means) == nrow(grid) && all(is.finite(draws)) && all(is.finite(means)),
    "Package returned invalid function draws.")
  q <- apply(draws, 2L, stats::quantile, probs = c(.025, .975), type = 7)
  # Average conditional means for the point estimate; do not add GP simulation
  # noise to this estimate. Quantiles include conditional GP uncertainty.
  data.frame(grid, posterior_mean = colMeans(means), q025 = q[1, ], q975 = q[2, ],
    Width95 = q[2, ] - q[1, ])
}

# 4. Package-only functional posterior calculations, cached separately per fit -
ocean_diag_functionals <- function(i, sources, output_dir, inputs, grids, settings,
  cache_dir, overwrite) {
  s <- sources[i, ]; path <- file.path(output_dir, s$source)
  fit_path <- file.path(path, s$fit_file)
  prediction_path <- file.path(path, "prediction_distributions.rds")
  tr <- inputs$train[[s$fold]]
  mg <- grids$m[grids$m$design == "candidate_A" | s$cell != "O0", ]
  identity <- list(version = "ocean-diagnostics-1", package = as.character(utils::packageVersion("eivGP")),
    source = s, hashes = unname(tools::md5sum(c(fit_path, prediction_path))),
    train = tr[, c("cohort_id", "G2salinity", "log_pressure", "C", "Y", "U", "calib_o20", "calib_o50")],
    m_grid = mg, f_grid = grids$f, settings = settings,
    code = list(body(ocean_diag_functionals), body(ocean_diag_balanced_ids), body(ocean_diag_point_summary)))
  cache_path <- file.path(cache_dir, sprintf("fold%02d_%s.rds", s$fold, s$cell))
  if (file.exists(cache_path) && !overwrite) {
    saved <- readRDS(cache_path)
    ocean_diag_require(identical(saved$identity, identity),
      paste("Diagnostic cache inputs changed. Choose a new report_dir or set overwrite=TRUE:", cache_path))
    return(saved$result)
  }
  message("Function diagnostics: fold ", s$fold, "/", s$cell)
  fit <- mixedgp_read_artifact(fit_path)
  mask <- switch(s$cell, O0 = rep(FALSE, nrow(tr)), O20 = tr$calib_o20, O50 = tr$calib_o50)
  xmat <- function(d) as.matrix(d[, c("G2salinity", "log_pressure"), drop = FALSE])
  ocean_diag_require(inherits(fit, "eivgp_fit") && identical(fit$sampler_version, "0.3.1") &&
    isTRUE(all.equal(unname(fit$data$x_raw), unname(xmat(tr)), tolerance = 1e-10)) &&
    isTRUE(all.equal(as.numeric(fit$data$y_raw), tr$Y, tolerance = 1e-10)) &&
    identical(as.integer(fit$data$c_ord), tr$C) &&
    setequal(fit$data$calib_idx, which(mask)) && setequal(fit$data$miss_idx, which(!mask)),
    "Fit object does not match the frozen training rows, model version, or calibration mask.")
  if (any(mask)) ocean_diag_require(
    isTRUE(all.equal(as.numeric(fit$data$u_obs_raw[mask]), tr$U[mask], tolerance = 1e-10)),
    "Exact oxygen measurements in the fitted object disagree with the calibration rows.")
  ids <- ocean_diag_balanced_ids(fit, settings$mean_states_per_chain)
  seed <- settings$seed + 100L * s$fold + match(s$cell, c("O0", "O20", "O50"))
  md <- eivGP::predict_eivgp(fit, new_X = xmat(mg), new_C = matrix(mg$C, ncol = 1),
    target = "mean", draw_ids = ids, n_latent = settings$n_latent,
    n_per_draw = settings$mean_gp_draws, include_process_uncertainty = TRUE,
    joint = FALSE, return_components = TRUE, seed = seed)
  conditional_means <- attr(md, "conditional_means", exact = TRUE)
  ms <- ocean_diag_point_summary(md, conditional_means, mg)
  ms$fold <- s$fold; ms$cell <- s$cell
  # Small numerical integration sensitivity check, not a recovery score. The
  # same 120 parameter states are evaluated with 256 and 512 latent nodes.
  check_columns <- match(1:4, mg$C)
  check_states <- unique(round(seq(1, length(ids), length.out = 120L)))
  check_grid <- mg[check_columns, ]
  alt <- eivGP::predict_eivgp(fit, new_X = xmat(check_grid),
    new_C = matrix(check_grid$C, ncol = 1), target = "mean", draw_ids = ids[check_states],
    n_latent = 2L * settings$n_latent, n_per_draw = 1L,
    include_process_uncertainty = FALSE, return_components = TRUE, seed = seed + 20000L)
  delta <- attr(alt, "conditional_means", exact = TRUE) -
    conditional_means[check_states, check_columns, drop = FALSE]
  check <- data.frame(fold = s$fold, cell = s$cell, check_grid,
    delta_mean_double_nodes = colMeans(delta), paired_state_mcse = apply(delta, 2, stats::sd) / sqrt(nrow(delta)))
  fs <- NULL; fids <- integer()
  if (s$cell != "O0") {
    fids <- readRDS(prediction_path)$draw_ids
    ocean_diag_require(length(fids) >= 1000L && !anyDuplicated(fids) &&
      all(fids >= 1 & fids <= nrow(fit$mcmc$samples_tau)), "Unexpected saved posterior state IDs.")
    fm <- eivGP::predict_eivgp(fit, new_X = xmat(grids$f), new_U = matrix(grids$f$U, ncol = 1),
      new_U_scale = "raw", target = "surface", draw_ids = fids, include_process_uncertainty = FALSE)
    fd <- eivGP::predict_eivgp(fit, new_X = xmat(grids$f), new_U = matrix(grids$f$U, ncol = 1),
      new_U_scale = "raw", target = "surface", draw_ids = fids, n_per_draw = settings$surface_gp_draws,
      include_process_uncertainty = TRUE, joint = FALSE, seed = seed + 10000L)
    fs <- ocean_diag_point_summary(fd, fm, grids$f)
    fs$fold <- s$fold; fs$cell <- s$cell
  }
  result <- list(m = ms, f = fs, integration_check = check,
    draw_ids = list(mean = ids, surface = fids), seed = seed)
  saveRDS(list(identity = identity, result = result), cache_path)
  result
}

# 5. Original held-out Y tasks; reuse saved scores without new predictions -----
ocean_diag_prediction <- function(output_dir, inputs, sources) {
  metrics <- c("RMSE", "MAE", "CRPS", "NLPD", "Coverage95", "Width95", "IntervalScore95")
  rows <- masks <- list()
  for (fold in 1:5) {
    path <- file.path(output_dir, "full5_class4", sprintf("fold%02d", fold))
    z <- ocean_diag_csv(file.path(path, "baseline_metrics.csv"), c("fold", "cell", "method", "n_test", metrics))
    ocean_diag_require(nrow(z) == 3L && setequal(z$method, c("UC-GP", "LVGP", "EzGP")) &&
      all(z$fold == fold & z$n_test == 40L), "Unexpected baseline scores.")
    z$task <- "class_only"; z$test_exact_fraction <- 0
    rows[[length(rows) + 1L]] <- z[, c("fold", "cell", "method", "task", "test_exact_fraction", "n_test", metrics)]
    mask <- ocean_diag_csv(file.path(path, "test_U_masks.csv"),
      c("cohort_id", "C", "observed_U_O20", "observed_U_O50"))
    te <- inputs$test[[fold]]
    ocean_diag_require(nrow(mask) == 40L && !anyDuplicated(mask$cohort_id) &&
      setequal(mask$cohort_id, te$cohort_id), "Unexpected held-out oxygen mask rows.")
    mask <- mask[match(te$cohort_id, mask$cohort_id), ]
    m20 <- ocean_diag_flag(mask$observed_U_O20); m50 <- ocean_diag_flag(mask$observed_U_O50)
    ocean_diag_require(all(mask$C == te$C) && all(m20 <= m50) &&
      all(table(factor(te$C[m20], levels = 1:4)) == 2L) &&
      all(table(factor(te$C[m50], levels = 1:4)) == 5L), "Held-out masks changed.")
    masks[[fold]] <- data.frame(fold = fold, mask)
    for (i in which(sources$fold == fold)) {
      s <- sources[i, ]
      z <- ocean_diag_csv(file.path(output_dir, s$source, "predictive_metrics.csv"),
        c("fold", "cell", "method", "n_test", metrics))
      ocean_diag_require(nrow(z) == 2L && setequal(z$method, c("EIV-GP xc", "EIV-GP mixed")) &&
        all(z$fold == fold & z$cell == s$cell & z$n_test == 40L), "Unexpected EIV prediction scores.")
      saved <- readRDS(file.path(output_dir, s$source, "prediction_distributions.rds"))
      known <- switch(s$cell, O0 = rep(FALSE, 40L), O20 = m20, O50 = m50)
      ocean_diag_require(identical(as.logical(saved$known), known) &&
        ncol(saved$xc) == 40L && ncol(saved$mixed) == 40L &&
        identical(saved$xc[, !known, drop = FALSE], saved$mixed[, !known, drop = FALSE]),
        "Saved mixed predictions disagree with the frozen exact-oxygen masks.")
      if (s$cell == "O0") z <- z[z$method == "EIV-GP xc", ] # O0 mixed is identical.
      z$task <- ifelse(z$method == "EIV-GP xc", "class_only", "matched_mixed")
      z$test_exact_fraction <- ifelse(z$task == "class_only", 0, if (s$cell == "O20") .2 else .5)
      z$method <- "EIV-GP"
      rows[[length(rows) + 1L]] <- z[, c("fold", "cell", "method", "task", "test_exact_fraction", "n_test", metrics)]
    }
  }
  raw <- do.call(rbind, rows)
  normalized <- raw
  s_y <- stats::sd(inputs$cohort$Y) # Scale by the SD of Y in the frozen cohort.
  for (metric in c("RMSE", "MAE", "CRPS", "Width95", "IntervalScore95"))
    normalized[[metric]] <- normalized[[metric]] / s_y
  normalized$NLPD <- normalized$NLPD - log(s_y) # Density transformation, not division.
  list(fold_metrics_raw = raw, fold_metrics_SD_Y = normalized,
    summary_raw = ocean_diag_fold_summary(raw, c("task", "method", "cell"), metrics),
    summary_SD_Y = ocean_diag_fold_summary(normalized, c("task", "method", "cell"), metrics),
    test_masks = do.call(rbind, masks), SD_Y = s_y)
}

# 6. Figures: one fixed fold in the main text, all folds in the appendix --------
ocean_diag_plots <- function(u, m, f, grids, representative_fold, baseline_same_x = NULL) {
  colors <- c(O0 = "#777777", O20 = "#2365a3", O50 = "#c76326")
  labels <- c(O0 = "0% exact training U", O20 = "20% exact training U", O50 = "50% exact training U")
  wrap <- function(x, width = 105L) paste(strwrap(x, width = width), collapse = "\n")
  style <- function(p, fill = FALSE) {
    p <- p + ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
        plot.caption = ggplot2::element_text(hjust = 0)) +
      ggplot2::scale_color_manual(values = colors, labels = labels, name = NULL)
    if (fill) p <- p + ggplot2::scale_fill_manual(values = colors, labels = labels, name = NULL)
    p
  }
  interval_note <- "Pointwise 95% credible intervals for the mean; response noise excluded."
  u_plot <- function(all_folds) {
    d <- u$display_points
    if (!all_folds) d <- d[d$fold == representative_fold, ]
    d$class <- factor(d$C, levels = 1:4)
    d$role <- factor(ifelse(d$observed, "Observed calibration (not imputed)", "Masked training U"),
      levels = c("Masked training U", "Observed calibration (not imputed)"))
    panel_labels <- c(O0 = "0% calibration\n160 masked / 0 observed\nUnanchored model coordinate",
      O20 = "20% calibration\n128 masked / 32 observed\nTraining-reference SD units",
      O50 = "50% calibration\n80 masked / 80 observed\nTraining-reference SD units")
    d$calibration <- factor(d$cell, levels = names(panel_labels), labels = unname(panel_labels))
    d$fold_label <- factor(paste("Fold", d$fold), levels = paste("Fold", 1:5))
    p <- ggplot2::ggplot(d, ggplot2::aes(plot_reference, plot_estimate, color = class)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
      ggplot2::geom_linerange(data = d[!d$observed, ],
        ggplot2::aes(ymin = plot_q025, ymax = plot_q975),
        alpha = .20, linewidth = .35, show.legend = FALSE) +
      ggplot2::geom_point(ggplot2::aes(shape = role), size = 1.5, alpha = .85, stroke = .65) +
      ggplot2::scale_color_manual(values = c("1" = "#1b9e77", "2" = "#d95f02",
        "3" = "#7570b3", "4" = "#0072B2"), name = "Oxygen class", drop = FALSE) +
      ggplot2::scale_shape_manual(values = c("Masked training U" = 16,
        "Observed calibration (not imputed)" = 1), name = NULL) +
      ggplot2::guides(shape = ggplot2::guide_legend(order = 1, override.aes = list(color = "grey25", alpha = 1)),
        color = ggplot2::guide_legend(order = 2, override.aes = list(shape = 16, alpha = 1))) +
      ggplot2::coord_cartesian(xlim = range(d$plot_reference) + c(-.1, .1),
        ylim = range(c(d$plot_q025, d$plot_q975, d$plot_estimate), na.rm = TRUE) + c(-.1, .1)) +
      ggplot2::labs(x = "Reference oxygen (training-reference SD units)", y = "U estimate (display units)",
        title = "Training oxygen: masked estimates and observed calibration inputs",
        subtitle = if (all_folds) "All five folds; 160 training profiles in every panel" else
          paste("Fixed fold", representative_fold, "| 160 training profiles in every panel"),
        caption = wrap(paste("Solid: masked training U with 95% posterior intervals. Hollow: exact calibration, not imputed.",
          "Reference is centered/scaled by each fold's full training reference U; O20/O50 posteriors use that same transformation.",
          "O0 stays on its unanchored model scale: identity-line agreement is illustrative, not physical-scale recovery.",
          "Training Y participates. Scores remain restricted to the common 80 hidden rows at O20/O50."))) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
        plot.caption = ggplot2::element_text(hjust = 0))
    if (all_folds) p + ggplot2::facet_grid(fold_label ~ calibration) else
      p + ggplot2::facet_wrap(~calibration, nrow = 1)
  }
  m_plot <- function(all_folds) {
    d <- m[m$design == "candidate_A", ]
    if (!all_folds) d <- d[d$fold == representative_fold, ]
    d$class <- factor(d$C, levels = 1:4, labels = paste("Class", 1:4))
    p <- ggplot2::ggplot(d, ggplot2::aes(pressure, posterior_mean, color = cell, fill = cell)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = q025, ymax = q975), alpha = .13, color = NA) +
      ggplot2::geom_line(linewidth = .65) +
      ggplot2::labs(x = "Pressure (dbar)", y = "Mean silicate m(x,c) (micromol/kg)",
        title = "Category-averaged response at salinity 34.20",
        subtitle = if (all_folds) "All five folds; each band belongs to its own fitted model" else
          paste("Fixed fold", representative_fold, "| candidate A"),
        caption = wrap(paste(interval_note, "Pressure ranges follow the prespecified local display rule.")))
    p <- if (all_folds) p + ggplot2::facet_wrap(~fold + class, ncol = 4, scales = "free") else
      p + ggplot2::facet_wrap(~class, ncol = 2, scales = "free")
    style(p, fill = TRUE)
  }
  fm_plot <- function(all_folds) {
    df <- f; dm <- m[m$design == "same_x", ]; db <- baseline_same_x
    if (!all_folds) {
      df <- df[df$fold == representative_fold, ]; dm <- dm[dm$fold == representative_fold, ]
      if(!is.null(db)) db <- db[db$fold == representative_fold, ]
    }
    df$model <- paste0("EIV-GP ",df$cell); dm$model <- paste0("EIV-GP ",dm$cell)
    means <- ocean_bind(list(dm,db))
    palette <- c("EIV-GP O20"="#2365a3","EIV-GP O50"="#c76326",
      "UC-GP (x,c)"="#7b3294","LVGP (x,c)"="#008837","EzGP (x,c)"="#a6761d")
    models <- names(palette)[names(palette) %in% means$model]
    offsets <- setNames(seq(-.26,.26,length.out=length(models)),models)
    means$class_position <- means$C + offsets[means$model]
    means$type <- ifelse(grepl("^EIV",means$model),"EIV posterior mean","Baseline fitted mean")
    # Match the physical response scale across panels; give the discrete panel
    # less width, and use larger points/capped intervals to reduce empty space.
    ylim <- range(c(df$q025,df$q975,dm$q025,dm$q975,means$posterior_mean),na.rm=TRUE)
    rug <- grids$near[grids$near$U >= min(df$U) & grids$near$U <= max(df$U), ]
    pf <- ggplot2::ggplot(df, ggplot2::aes(U, posterior_mean, color = model, fill = model)) +
      ggplot2::geom_vline(xintercept = c(50, 100, 150), linetype = 3, color = "grey65", linewidth = .4) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = q025, ymax = q975), alpha = .14, color = NA) +
      ggplot2::geom_line(linewidth = .9) +
      ggplot2::geom_rug(data = rug, ggplot2::aes(x = U), inherit.aes = FALSE,
        sides = "b", color = "grey45", alpha = .65) +
      ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=.025)) +
      ggplot2::coord_cartesian(xlim = range(df$U), ylim = ylim) +
      ggplot2::labs(x = "Oxygen (micromol/kg)", y = "Mean silicate (micromol/kg)", title = "(a) Continuous response f(x,u)")
    pm <- ggplot2::ggplot(means, ggplot2::aes(class_position,posterior_mean,color=model)) +
      ggplot2::geom_errorbar(data=means[grepl("^EIV",means$model),],
        ggplot2::aes(ymin=q025,ymax=q975),width=.10,linewidth=.85,show.legend=FALSE) +
      ggplot2::geom_point(ggplot2::aes(shape=type),size=3,stroke=.8) +
      ggplot2::scale_shape_manual(values=c("EIV posterior mean"=16,"Baseline fitted mean"=17),guide="none") +
      ggplot2::scale_x_continuous(breaks=1:4,limits=c(.5,4.5),expand=c(0,0)) +
      ggplot2::coord_cartesian(ylim = ylim) +
      ggplot2::labs(x = "Recorded oxygen class", y = NULL, title = "(b) Category mean m(x,c)")
    theme_fm <- function(p) p + ggplot2::theme_minimal(base_size=12) +
      ggplot2::scale_color_manual(values=palette,limits=models,breaks=models,drop=FALSE,name=NULL) +
      ggplot2::scale_fill_manual(values=palette,guide="none") +
      ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),legend.position="bottom",
        plot.title=ggplot2::element_text(size=12),legend.key.width=grid::unit(1.1,"cm")) +
      ggplot2::guides(color=ggplot2::guide_legend(nrow=1,override.aes=list(shape=16,linewidth=.8)))
    if (all_folds) { pf <- pf + ggplot2::facet_grid(fold ~ .); pm <- pm + ggplot2::facet_grid(fold ~ .) }
    patchwork::wrap_plots(theme_fm(pf)+ggplot2::guides(color="none"),
      theme_fm(pm)+ggplot2::guides(color=ggplot2::guide_legend(nrow=1,
        override.aes=list(shape=ifelse(grepl("^EIV",models),16,17),size=3))),
      widths=c(1.65,1),guides="collect") +
      patchwork::plot_annotation(title = "Salinity 34.25, pressure 650 dbar",
        subtitle = if (all_folds) "All five folds; separate fitted posteriors" else paste("Fixed fold", representative_fold),
        caption = paste("EIV bars/shading: pointwise 95% mean credible intervals. Triangles: original (x,c) baseline fitted means.",
          "\nm integrates the full fitted oxygen distribution. Dotted lines: recorded band boundaries; rug: nearby reference oxygen.")) &
      ggplot2::theme(legend.position="bottom",plot.caption=ggplot2::element_text(hjust=0,size=10))
  }
  list(u_main = u_plot(FALSE), m_main = m_plot(FALSE), fm_main = fm_plot(FALSE),
    u_all_folds = u_plot(TRUE), m_all_folds = m_plot(TRUE), fm_all_folds = fm_plot(TRUE))
}

# 7. Rmd entry point: postprocess and export; no fitting ------------------------
ocean_build_diagnostics <- function(output_dir, data_dir = file.path(output_dir, "data"),
  report_dir = file.path(output_dir, "diagnostics"), representative_fold = 1L,
  workers = 1L, overwrite = FALSE) {
  for (package in c("eivGP", "ggplot2", "patchwork"))
    ocean_diag_require(requireNamespace(package, quietly = TRUE), paste("Required package:", package))
  ocean_diag_require(utils::packageVersion("eivGP") >= "0.4.0",
    "Install eivGP 0.4.0 or later; retained fits must use the supported 0.3.1 sampler.")
  ocean_diag_require(length(representative_fold) == 1L && representative_fold %in% 1:5,
    "representative_fold must be one of 1:5; the default main display uses fold 1.")
  ocean_diag_require(length(workers) == 1L && is.finite(workers) && workers >= 1 && workers == as.integer(workers),
    "workers must be a positive integer.")
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  inputs <- ocean_diag_inputs(data_dir)
  sources <- ocean_diag_sources(output_dir)
  grids <- ocean_diag_grids(inputs$cohort)
  u <- ocean_diag_u(output_dir, inputs, sources)
  prediction <- ocean_diag_prediction(output_dir, inputs, sources)
  settings <- list(mean_states_per_chain = 60L, n_latent = 256L, mean_gp_draws = 20L,
    surface_gp_draws = 4L, seed = 202610080L)
  for (subdir in c("tables", "figures", "internal", "cache"))
    dir.create(file.path(report_dir, subdir), recursive = TRUE, showWarnings = FALSE)
  run_one <- function(i) ocean_diag_functionals(i, sources, output_dir, inputs, grids, settings,
    cache_dir = file.path(report_dir, "cache"), overwrite = overwrite)
  if (workers > 1L && .Platform$OS.type != "windows") {
    parts <- parallel::mclapply(seq_len(nrow(sources)), run_one, mc.cores = min(workers, nrow(sources)),
      mc.preschedule = FALSE, mc.set.seed = FALSE) # Explicit seed for every package call.
    ocean_diag_require(!any(vapply(parts, inherits, logical(1), "try-error")),
      "At least one function calculation failed. Inspect messages; completed cell caches can be reused.")
  } else parts <- lapply(seq_len(nrow(sources)), run_one)
  m <- do.call(rbind, lapply(parts, `[[`, "m"))
  f <- do.call(rbind, lapply(parts, `[[`, "f"))
  checks <- do.call(rbind, lapply(parts, `[[`, "integration_check"))
  ocean_diag_require(nrow(m) == 15L * 36L + 10L * 4L && nrow(f) == 10L * 41L &&
    all(m$Width95 > 0) && all(f$Width95 > 0), "Incomplete functional diagnostics.")
  # Internal width averages are descriptive precision checks, not true m/f
  # recovery metrics or new main-text tables. Equal weights over the fixed grid.
  widths <- stats::aggregate(Width95 ~ fold + cell + C, m[m$design == "candidate_A", ], mean)
  width_summary <- ocean_diag_fold_summary(widths, c("cell", "C"), "Width95")
  diag <- lapply(seq_len(nrow(sources)), function(i) {
    s <- sources[i, ]; d <- utils::read.csv(file.path(output_dir, s$source, "convergence_summary.csv"))
    data.frame(fold = s$fold, cell = s$cell, d)
  })
  baseline_paths <- file.path(output_dir,"full5_class4",sprintf("fold%02d",1:5),"baseline_package_outputs.rds")
  baseline_same_x <- if(all(file.exists(baseline_paths)))
    ocean_baseline_mean_grid(output_dir,grids$m[grids$m$design=="same_x",]) else NULL
  plots <- ocean_diag_plots(u,m,f,grids,representative_fold,baseline_same_x)
  write_csv <- function(x, subdir, name) utils::write.csv(x, file.path(report_dir, subdir, name), row.names = FALSE)
  write_csv(u$points, "tables", "u_common_missing_points.csv")
  write_csv(u$display_points, "tables", "u_display_points.csv")
  write_csv(u$display_scales, "internal", "u_display_scales.csv")
  write_csv(u$display_counts, "internal", "u_display_counts.csv")
  write_csv(u$fold_metrics, "tables", "u_fold_metrics.csv")
  write_csv(u$summary, "tables", "u_summary.csv")
  write_csv(m, "tables", "m_points.csv"); write_csv(f, "tables", "f_points.csv")
  if(!is.null(baseline_same_x)) write_csv(baseline_same_x,"tables","baseline_m_same_x.csv")
  paper_u <- ocean_paper_u_tables(u$summary)
  ocean_export_paper_table(paper_u$main,"paper_u_imputation_main","Oxygen imputation: common 80 masked training profiles",file.path(report_dir,"tables"),paper_u$note)
  ocean_export_paper_table(paper_u$class,"paper_u_imputation_class","Class-specific oxygen imputation",file.path(report_dir,"tables"),paper_u$note)
  write_csv(prediction$fold_metrics_raw, "tables", "y_prediction_fold_raw.csv")
  write_csv(prediction$fold_metrics_SD_Y, "tables", "y_prediction_fold_SD_Y.csv")
  write_csv(prediction$summary_raw, "tables", "y_prediction_summary_raw.csv")
  write_csv(prediction$summary_SD_Y, "tables", "y_prediction_summary_SD_Y.csv")
  write_csv(inputs$audit, "internal", "input_audit.csv")
  write_csv(sources, "internal", "fit_sources.csv")
  write_csv(grids$m, "internal", "m_grid.csv"); write_csv(grids$f, "internal", "f_grid.csv")
  write_csv(grids$m_neighbors, "internal", "m_display_neighbors.csv")
  write_csv(grids$near, "internal", "f_display_neighbors.csv")
  write_csv(grids$support, "internal", "f_display_support.csv")
  write_csv(checks, "internal", "m_integration_sensitivity.csv")
  write_csv(widths, "internal", "m_width_fold.csv")
  write_csv(width_summary, "internal", "m_width_summary.csv")
  write_csv(prediction$test_masks, "internal", "test_U_masks.csv")
  write_csv(do.call(rbind, diag), "internal", "mcmc_existing_summary.csv")
  dimensions <- list(u_main = c(13.2, 5.2), m_main = c(10, 7), fm_main = c(11, 4.9),
    u_all_folds = c(13.2, 14.5), m_all_folds = c(13, 13), fm_all_folds = c(11, 14))
  for (name in names(plots)) for (extension in c("png", "pdf"))
    ggplot2::ggsave(file.path(report_dir, "figures", paste0(name, ".", extension)), plots[[name]],
      width = dimensions[[name]][1], height = dimensions[[name]][2], dpi = 180, bg = "white")
  provenance <- list(settings = settings, representative_fold = representative_fold,
    package = as.character(utils::packageVersion("eivGP")), upstream_commit = "671ce98339807cd3a13347e8468f36acf097d50a",
    input_hashes = tools::md5sum(file.path(data_dir, c("cohort_manifest.csv",
      sprintf("fold%02d_train.csv", 1:5), sprintf("fold%02d_test.csv", 1:5)))),
    source_hashes = lapply(seq_len(nrow(sources)), function(i) {
      s <- sources[i, ]
      tools::md5sum(file.path(output_dir, s$source, c(s$fit_file,
        "prediction_distributions.rds", "predictive_metrics.csv", "U_posterior_intervals.csv")))
    }),
    draw_ids = lapply(parts, `[[`, "draw_ids"), SD_Y = prediction$SD_Y,
    main_targets = c("training U imputation", "m(x,c)", "f(x,u) and same-X m(x,c)", "held-out Y prediction"))
  saveRDS(provenance, file.path(report_dir, "internal", "provenance.rds"))
  capture.output(utils::sessionInfo(), file = file.path(report_dir, "internal", "sessionInfo.txt"))
  invisible(list(u = u, m = m, f = f, prediction = prediction, grids = grids, plots = plots, baseline_same_x=baseline_same_x,
    integration_check = checks, sources = sources, report_dir = report_dir, provenance = provenance))
}

# 8. Compact HTML/Rmd presentation: the four diagnostic blocks ----------
ocean_diag_display_table <- function(summary, groups, metrics, digits = 3L) {
  rows <- unique(summary[, groups, drop = FALSE])
  key <- function(d) do.call(paste, c(d[groups], sep = "/"))
  for (metric in metrics) {
    z <- summary[summary$metric == metric, ]
    i <- match(key(rows), key(z))
    rows[[metric]] <- sprintf(paste0("%.", digits, "f +/- %.", digits, "f"), z$mean[i], z$sd[i])
  }
  rownames(rows) <- NULL
  rows
}

ocean_show_diagnostics <- function(report, appendix = FALSE, show_prediction = TRUE) {
  ocean_diag_require(requireNamespace("knitr", quietly = TRUE), "knitr is required for notebook presentation.")
  table <- function(x) cat(as.character(knitr::kable(x, format = "pipe", row.names = FALSE)), sep = "\n")
  cat("\n\n### Training oxygen imputation\n\n",
    "Both calibration levels are scored on the same 80 hidden training observations per fold. ",
    "Reference oxygen is the masked source measurement. Entries are the mean +/- sample SD of five fold scores; ",
    "training sets overlap, and these SDs are descriptive. Training Y contributes to imputation.\n\n", sep = "")
  table(ocean_paper_u_tables(report$u$summary)$main)
  cat("\n\nErrors and widths use micromol/kg; coverage is percent. O0 has no anchored raw oxygen scale.\n\n")
  cat("\n\nThe figure displays all 160 training rows per calibration level (O0/O20/O50). ",
    "Solid points with intervals are masked inputs; hollow points are observed calibration values, not imputations. ",
    "Reference U is centered/scaled by the full training-reference mean/SD within each fold. ",
    "O20/O50 posteriors are transformed to those same units; O0 remains in its original unanchored model coordinate. ",
    "The O0 identity-line comparison is illustrative only and supplies no physical-scale recovery score. ",
    "Hidden references enter only this offline display, never fitting or posterior alignment. ",
    "Displayed masked-row counts are 160/128/80; the table still scores the same common 80 rows at O20/O50.\n\n", sep = "")
  print(report$plots$u_main)
  cat("\n\n### Category-averaged silicate response m(x,c)\n\n",
    "The pressure slice fixes salinity at 34.20 and varies pressure within a class-specific local display range. ",
    "The package integrates f against the fitted conditional oxygen distribution. Bands describe uncertainty ",
    "in this mean and exclude response noise. True m is unknown, so no mean-recovery error is reported.\n\n", sep = "")
  print(report$plots$m_main)
  cat("\n\n### Continuous response f(x,u) and the category average at the same X\n\n",
    "Salinity is 34.25 and pressure is 650 dbar. The curve varies oxygen continuously; the four category means ",
    "average over their full fitted conditional oxygen distributions. The displayed oxygen range is the ",
    "10th-90th percentile among 30 nearest reference profiles in standardized (salinity, log pressure). ",
    "Reference oxygen defines only this offline display domain and rug, not additional calibration. ",
    "Nearby profiles indicate the available input support but do not guarantee dense support at every joint input. ",
    "The fitted end-class distributions may extend beyond 0-180 and beyond the plotted range. ",
    "True f is unknown; this is a descriptive fitted association, without a surface-recovery score or a causal interpretation.\n\n", sep = "")
  print(report$plots$fm_main)
  if (show_prediction) {
  cat("\n\n### Held-out silicate prediction\n\n",
    "The original prediction tasks are retained. Class-only uses no exact test oxygen. Matched mixed uses ",
    "20% or 50% exact test oxygen with the corresponding training calibration level. O0 mixed duplicates ",
    "O0 class-only and is shown once. The three categorical comparators use class-only information. ",
    "All entries are five-fold mean +/- SD; no new held-out predictions are generated here.\n\n", sep = "")
  metrics <- c("RMSE", "CRPS", "NLPD", "Coverage95", "Width95", "IntervalScore95")
  for (task in c("class_only", "matched_mixed")) {
    cat("\n\n**", if (task == "class_only") "Class-only" else "Matched mixed", "**\n\n", sep = "")
    s <- report$prediction$summary_SD_Y
    table(ocean_diag_display_table(s[s$task == task, ], c("method", "cell"), metrics))
  }
  cat("\n\nRMSE, CRPS, width and interval score are divided by the frozen cohort SD(Y); ",
    "NLPD is transformed by subtracting log(SD(Y)); coverage is unchanged. ",
    "Raw-unit tables and MAE are also exported. Mixed-versus-baseline differences combine model and input-information differences.\n\n", sep = "")
  }
  if (appendix) {
    cat("\n\n### Appendix: all five fitted folds\n\n",
      "Each interval belongs to one fitted fold. These plots do not pool the folds into a posterior.\n\n", sep = "")
    print(report$plots$u_all_folds); print(report$plots$m_all_folds); print(report$plots$fm_all_folds)
    cat("\n\n**Class-specific U imputation (same five metrics)**\n\n")
    table(ocean_diag_display_table(report$u$summary[report$u$summary$group != "overall", ], c("cell", "group"),
      c("Bias", "RMSE", "MAE", "Coverage95", "Width95")))
  }
  invisible(report)
}

# Native baseline mean predictors; no custom covariance or GP posterior algebra.
ocean_native_mean <- function(method,fit,grid) {
  X <- as.matrix(grid[,c("G2salinity","log_pressure")]); C <- matrix(grid$C,ncol=1L)
  if(method=="UC-GP") {
    d <- make_mixedgp_competitor_data(X,C,4L)
    as.numeric(stats::predict(fit,newdata=d,type="UK",seCompute=FALSE,
      covCompute=FALSE,lightReturn=TRUE)$mean)
  } else {
    d <- make_mixedgp_numeric_matrix(X,C,4L)
    if(method=="LVGP") as.numeric(LVGP::LVGP_predict(d,fit,MSE_on=0)$Y_hat)
    else as.numeric(EzGP::EzGP_predict(d,fit,MSE_on=0)$Y_hat)
  }
}

# The original three baselines predict a category response directly. Evaluate
# their native means on the fixed X grids; no R/W slice is relabeled as m.
ocean_baseline_mean_grid <- function(output_dir,grid,folds=1:5) {
  out <- list()
  for(fold in folds) {
    saved <- readRDS(file.path(output_dir,"full5_class4",sprintf("fold%02d",fold),"baseline_package_outputs.rds"))
    for(method in c("UC-GP","LVGP","EzGP")) {
      mean <- ocean_native_mean(method,saved$fits[[method]],grid)
      stopifnot(length(mean)==nrow(grid),all(is.finite(mean)))
      out[[length(out)+1L]] <- data.frame(grid,posterior_mean=mean,fold,model=paste0(method," (x,c)"))
    }
  }
  ocean_bind(out)
}

# Input-based display selection. Use the fitted grid and its local
# salinity pool; show at most 3 nearest-salinity observations per log-pressure
# bin, within the curve range. Y, residuals and held-out scores are never used.
ocean_m_display_points <- function(cohort,grids,n_bins=6L,per_bin=3L) {
  out <- list()
  for(C in 1:4) {
    g <- grids$m[grids$m$design=="candidate_A" & grids$m$C==C,]
    ids <- grids$m_neighbors$cohort_id[grids$m_neighbors$class==C]
    d <- cohort[cohort$cohort_id %in% ids & cohort$log_pressure>=min(g$log_pressure) &
      cohort$log_pressure<=max(g$log_pressure),]
    d$pressure_bin <- cut(d$log_pressure,breaks=seq(min(g$log_pressure),max(g$log_pressure),length.out=n_bins+1L),
      include.lowest=TRUE,labels=FALSE)
    d$salinity_offset <- d$G2salinity-34.20
    keep <- unlist(lapply(split(seq_len(nrow(d)),d$pressure_bin),function(i)
      head(i[order(abs(d$salinity_offset[i]),d$cohort_id[i])],per_bin)),use.names=FALSE)
    out[[C]] <- d[keep,,drop=FALSE]
  }
  ocean_bind(out)
}

ocean_paper_u_tables <- function(summary) {
  format_u <- function(groups) {
    z <- summary[summary$group %in% groups,]
    keys <- unique(z[c("cell","group")]); keys <- keys[order(keys$cell,keys$group),]
    d <- data.frame(Calibration=paste0(sub("O","",keys$cell),"%"),check.names=FALSE)
    if(length(groups)>1L) d$Class <- sub("class","",keys$group)
    ranked <- as.data.frame(matrix(NA_real_,nrow(keys),3L))
    names(ranked) <- c("Bias","RMSE","MAE")
    eligible <- rep(TRUE,nrow(keys))
    for(metric in c("Bias","RMSE","MAE","Coverage95","Width95")) {
      v <- z[z$metric==metric,]; i <- match(paste(keys$cell,keys$group),paste(v$cell,v$group))
      mult <- if(metric=="Coverage95") 100 else 1
      label <- if(metric=="Coverage95") "Coverage (%)" else if(metric=="Width95") "Width" else metric
      d[[label]] <- ocean_mean_sd(mult*v$mean[i],mult*v$sd[i],if(metric=="Coverage95") 1 else 2)
      eligible <- eligible & v$n_folds[i]==5L
      if(metric %in% names(ranked))
        ranked[[metric]] <- if(metric=="Bias") abs(v$mean[i]) else v$mean[i]
      if(any(v$n_folds[i]!=5L)) d[["n folds"]] <- v$n_folds[i]
    }
    attr(d,"best") <- ocean_best_mask(ranked,keys$group,eligible)
    d
  }
  list(main=format_u("overall"),class=format_u(paste0("class",1:4)),
    note=paste("Entries are five-fold means (sample SDs), computed on the same 80 masked training profiles per fold.",
      "Bias, RMSE, MAE and width are in micromol/kg; coverage is percent for individual-U 95% posterior intervals.",
      "Bold marks the smaller absolute mean bias, RMSE and MAE across calibration levels for the same target group, using unrounded values. Coverage and width are descriptive."))
}

# Descriptive data analysis: measured nearby responses provide data context.
ocean_surface_report <- function(output_dir,data_dir,diagnostics_dir=file.path(output_dir,"diagnostics"),
  representative_fold=1L,folds=1:5) {
  inputs <- ocean_diag_inputs(data_dir); grids <- ocean_diag_grids(inputs$cohort)
  m <- read.csv(file.path(diagnostics_dir,"tables/m_points.csv"))
  f <- read.csv(file.path(diagnostics_dir,"tables/f_points.csv"))
  baseline_all <- ocean_baseline_mean_grid(output_dir,grids$m,folds)
  baseline <- baseline_all[baseline_all$design=="candidate_A",]
  baseline_same <- baseline_all[baseline_all$design=="same_x",]
  selected_m <- ocean_m_display_points(inputs$cohort,grids)
  points_m <- points_f <- list()
  for(fold in folds) {
    training <- inputs$train[[fold]]$cohort_id
    dm <- selected_m; dm$role <- ifelse(dm$cohort_id %in% training,"Training","Held-out")
    dm$fold <- fold; points_m[[length(points_m)+1L]] <- dm
    df <- grids$near[grids$near$U>=min(grids$f$U) & grids$near$U<=max(grids$f$U),]
    df$role <- ifelse(df$cohort_id %in% training,"Training","Held-out")
    df$fold <- fold; points_f[[length(points_f)+1L]] <- df
  }
  points_m <- ocean_bind(points_m); points_f <- ocean_bind(points_f)
  em <- m[m$design=="candidate_A" & m$fold %in% folds,]; em$model <- paste0("EIV-GP ",em$cell)
  lines <- ocean_bind(list(em,baseline))
  palette <- c("EIV-GP O0"="#777777","EIV-GP O20"="#2365a3","EIV-GP O50"="#c76326",
    "UC-GP (x,c)"="#7b3294","LVGP (x,c)"="#008837","EzGP (x,c)"="#a6761d")
  plot_m <- function(all=FALSE) {
    dl <- lines; de <- em; dp <- points_m
    if(!all) {dl<-dl[dl$fold==representative_fold,];de<-de[de$fold==representative_fold,];dp<-dp[dp$fold==representative_fold,]}
    dl$class_label <- paste("Class",dl$C); de$class_label <- paste("Class",de$C); dp$class_label <- paste("Class",dp$C)
    p <- ggplot2::ggplot(dl,ggplot2::aes(exp(log_pressure),posterior_mean,color=model)) +
      ggplot2::geom_ribbon(data=de,ggplot2::aes(ymin=q025,ymax=q975,fill=model),alpha=.08,color=NA) +
      ggplot2::geom_point(data=dp,ggplot2::aes(G2pressure,Y,shape=role),inherit.aes=FALSE,
        color="grey45",size=2,alpha=.7,stroke=.65) +
      ggplot2::geom_line(linewidth=.8) +
      ggplot2::scale_color_manual(values=palette,name=NULL) +
      ggplot2::scale_fill_manual(values=palette,guide="none") +
      ggplot2::scale_shape_manual(values=c(Training=16,"Held-out"=1),name="Nearby profiles") +
      ggplot2::guides(color=ggplot2::guide_legend(nrow=2,order=2),shape=ggplot2::guide_legend(order=1)) +
      ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=.05)) +
      ggplot2::labs(x="Pressure (dbar)",y="Silicate (micromol/kg)",
        title="Category response: fitted means and nearby observations",
        subtitle="Fixed salinity 34.20; separate oxygen classes",
        caption="Points: up to 3 nearest-salinity profiles per pressure bin, within the displayed curve range.\nEIV shading: pointwise 95% mean credible intervals. Baseline lines: original (x,c) fitted means.") +
      ggplot2::theme_minimal(base_size=12) + ggplot2::theme(legend.position="bottom",
        panel.grid.minor=ggplot2::element_blank(),plot.caption=ggplot2::element_text(hjust=0,size=10))
    if(all) p+ggplot2::facet_wrap(~fold+class_label,ncol=4,scales="free") else p+ggplot2::facet_wrap(~class_label,ncol=2,scales="free")
  }
  plot_f <- function(all=FALSE) {
    df <- f[f$fold %in% folds,]; dp <- points_f
    if(!all) {df<-df[df$fold==representative_fold,];dp<-dp[dp$fold==representative_fold,]}
    p <- ggplot2::ggplot(df,ggplot2::aes(U,posterior_mean,color=cell,fill=cell)) +
      ggplot2::geom_vline(xintercept=c(50,100,150),linetype=3,color="grey60") +
      ggplot2::geom_ribbon(ggplot2::aes(ymin=q025,ymax=q975),alpha=.12,color=NA) +
      ggplot2::geom_line(linewidth=.7) +
      ggplot2::geom_point(data=dp,ggplot2::aes(U,Y,shape=role),inherit.aes=FALSE,
        color="grey35",size=1.8,alpha=.7) +
      ggplot2::scale_color_manual(values=c(O20="#2365a3",O50="#c76326"),name="Training calibration") +
      ggplot2::scale_fill_manual(values=c(O20="#2365a3",O50="#c76326"),guide="none") +
      ggplot2::scale_shape_manual(values=c(Training=16,"Held-out"=1),name="Nearby observed profiles") +
      ggplot2::labs(x="Oxygen (micromol/kg)",y="Silicate (micromol/kg)",
        title="Continuous response f(x,u): salinity 34.25, pressure 650 dbar",
        caption="Points: nearby profiles with measured oxygen and silicate; their salinity/pressure vary.\nShading: pointwise 95% credible intervals for the mean, excluding response noise.") +
      ggplot2::theme_minimal(base_size=11) + ggplot2::theme(legend.position="bottom")
    if(all) p+ggplot2::facet_wrap(~fold,ncol=2) else p
  }
  plots <- list(m_comparison=plot_m(),f_with_observations=plot_f(),
    m_comparison_all_folds=plot_m(TRUE),f_with_observations_all_folds=plot_f(TRUE))
  dest <- file.path(diagnostics_dir,"surface_comparison")
  ocean_atomic_csv(baseline,file.path(dest,"baseline_m_points.csv"))
  ocean_atomic_csv(baseline_same,file.path(dest,"baseline_m_same_x.csv"))
  ocean_atomic_csv(points_m,file.path(dest,"nearby_m_observations.csv"))
  ocean_atomic_csv(selected_m,file.path(dest,"m_selection_audit.csv"))
  ocean_atomic_csv(points_f,file.path(dest,"nearby_f_observations.csv"))
  for(n in names(plots)) for(ext in c("png","pdf")) ggplot2::ggsave(file.path(dest,paste0(n,".",ext)),plots[[n]],
    width=if(grepl("fm_",n)) 11 else if(grepl("m_comparison",n)) 11 else 10,
    height=if(grepl("all_folds",n)) 13 else if(grepl("fm_",n)) 4.9 else 7,dpi=180,bg="white")
  list(plots=plots,baseline=baseline,baseline_same_x=baseline_same,points_m=points_m,points_f=points_f)
}

# Unified postprocessing entry used by both Rmd and this script's CLI.
# No fitting routine is called. Optional copied caches retain their full identity
# checks inside ocean_build_diagnostics; changed fits/settings are rejected.
ocean_run_diagnostics <- function(output_dir,data_dir,
  report_dir=file.path(output_dir,"diagnostics"),augmented_output_dir="",
  representative_fold=1L,workers=1L,cache_from="",
  build_prediction=TRUE,build_surfaces=TRUE,prediction_report=NULL) {
  for(p in c("eivGP","ggplot2","patchwork","digest"))
    ocean_diag_require(requireNamespace(p,quietly=TRUE),paste("Required package:",p))
  output_dir <- normalizePath(output_dir,mustWork=TRUE)
  data_dir <- normalizePath(data_dir,mustWork=TRUE)
  dir.create(report_dir,recursive=TRUE,showWarnings=FALSE)
  report_dir <- normalizePath(report_dir,mustWork=TRUE)
  progress <- function(stage) ocean_atomic_lines(c(paste("stage:",stage),
    paste("updated:",format(Sys.time(),"%Y-%m-%d %H:%M:%S %z")),
    paste("pid:",Sys.getpid()),"Postprocessing saved fits only; no MCMC or baseline refit."),
    file.path(report_dir,"PROGRESS.txt"))
  progress("input and source audit")
  inputs <- ocean_prepare_inputs(data_dir)
  sources <- ocean_diag_sources(output_dir)
  copied <- character()
  if(nzchar(cache_from)) {
    from <- normalizePath(cache_from,mustWork=TRUE)
    dir.create(file.path(report_dir,"cache"),showWarnings=FALSE)
    for(n in sprintf("fold%02d_%s.rds",sources$fold,sources$cell)) {
      source <- file.path(from,n); target <- file.path(report_dir,"cache",n)
      if(file.exists(source) && !file.exists(target)) {
        stopifnot(file.copy(source,target)); copied <- c(copied,n)
      }
    }
  }
  progress("U imputation and cached package m/f calculations")
  diagnostics <- ocean_build_diagnostics(output_dir,data_dir,report_dir=report_dir,
    representative_fold=representative_fold,workers=workers)
  prediction <- prediction_report; surfaces <- NULL
  if(build_prediction) {
    progress("original and O20/O50 expanded baseline prediction summaries")
    prediction <- ocean_prediction_report(output_dir,inputs,augmented_output_dir=augmented_output_dir,
      report_dir=file.path(report_dir,"prediction"))
  }
  if(!build_prediction && !is.null(prediction))
    ocean_export_prediction_report(prediction,file.path(report_dir,"prediction"))
  if(build_surfaces) {
    progress("native baseline m comparison and nearby measured-profile overlays")
    surfaces <- ocean_surface_report(output_dir,data_dir,report_dir,
      representative_fold=representative_fold)
  }
  # Paper-facing outputs use the selected b=0 encoding only. Raw audit
  # exports retain scores for every available encoding.
  main <- c("U imputation (O0/O20/O50, hollow calibration)"="figures/u_main.png",
    "m: native baselines and pressure-balanced nearby profiles"="surface_comparison/m_comparison.png",
    "f: S=34.25, P=650, measured nearby profiles"="surface_comparison/f_with_observations.png")
  main <- main[file.exists(file.path(report_dir,main))]
  ut <- ocean_paper_u_tables(diagnostics$u$summary)
  yt <- if(!is.null(prediction)) ocean_paper_prediction_tables(prediction) else NULL
  unote <- ut$note
  ynote <- if(!is.null(prediction)) paste0("Entries are fold means (sample SDs). sY = ",sprintf("%.4f",prediction$SD_Y),
    " micromol/kg. Coverage is percent. * marks baselines receiving the same available exact oxygen measurements using (x,c,R,W) inputs and b=0.") else ""
  text <- c("# Ocean diagnostics review","","Saved final fits; matched-oxygen baselines are marked with *.","",
    "## Paper table: U imputation","",ocean_table_html(ut$main,"Training oxygen imputation",unote),"",
    "[CSV](tables/paper_u_imputation_main.csv) | [LaTeX](tables/paper_u_imputation_main.tex)","")
  if(!is.null(yt) && nrow(yt$main)) text <- c(text,"## Paper table: Y prediction","",
    ocean_prediction_panel_html(yt$main,"Silicate prediction under matched oxygen information",
      paste(ynote,"Bold marks the smallest RMSE and CRPS mean within each complete five-fold information block, using unrounded values. Coverage and width are descriptive.")),"",
    "[CSV](prediction/paper_y_prediction_main.csv) | [LaTeX](prediction/paper_y_prediction_main.tex)","")
  text <- c(text,"## Main displays","",unlist(Map(function(label,path) c(paste0("### ",label),"",
    paste0("![",label,"](",path,")"),""),names(main),unname(main))),
    "## Additional tables","",
    "- [Class-specific U imputation](tables/paper_u_imputation_class.csv)",
    "- [Full Y-prediction summaries, b=0](prediction/paper_y_prediction_full.csv)",
    "- [Unrounded primary Y summary](prediction/primary_summary_numeric.csv)",
    "- [Same-X native baseline means](tables/baseline_m_same_x.csv)",
    "- [m measured-point selection audit](surface_comparison/m_selection_audit.csv)","",
    "Nearby points are measured responses at varying X. Curves and fitted posteriors are unchanged by the display-only point selection.",
    "Baseline category means are original (x,c) fits. Expanded baselines appear in matched Y-prediction comparisons.")
  ocean_atomic_lines(text,file.path(report_dir,"INDEX.md"))
  escape <- function(x) {x<-gsub("&","&amp;",x,fixed=TRUE);x<-gsub("<","&lt;",x,fixed=TRUE);gsub('"',"&quot;",x,fixed=TRUE)}
  html <- c('<!doctype html><meta charset="utf-8"><title>Ocean diagnostics review</title>',
    '<style>body{max-width:1150px;margin:32px auto;padding:0 16px;font:16px/1.5 sans-serif;color:#222}img{max-width:100%}section{margin:32px 0}table{border-collapse:collapse;width:100%;font-size:14px}td,th{padding:7px 9px;border-bottom:1px solid #ddd;text-align:right}td:nth-child(-n+3),th:nth-child(-n+3){text-align:left}caption{text-align:left;font-weight:600;padding:10px 0}.note{font-size:13px;color:#555}.table-scroll{overflow-x:auto}a{color:#2365a3}</style>',
    '<h1>Ocean diagnostics review</h1><p>Saved final fits; matched-oxygen baselines are marked with *. Tables show fold means (sample SDs).</p>',
    '<section><h2>Paper table: oxygen imputation</h2><div class="table-scroll">',
    ocean_table_html(ut$main,"Common 80 masked training profiles per fold",unote),
    '</div><p><a href="tables/paper_u_imputation_main.csv">CSV</a> \u00b7 <a href="tables/paper_u_imputation_main.tex">LaTeX</a> \u00b7 <a href="tables/paper_u_imputation_class.html">Class-specific table</a></p></section>')
  if(!is.null(yt) && nrow(yt$main)) html <- c(html,
    '<section><h2>Paper table: silicate prediction</h2><div class="table-scroll">',
    ocean_table_html(yt$main,"Original O0 comparators and matched-oxygen O20/O50 baselines",ynote),
    '</div><p><a href="prediction/paper_y_prediction_main.csv">CSV</a> \u00b7 <a href="prediction/paper_y_prediction_main.tex">LaTeX</a> \u00b7 <a href="prediction/paper_y_prediction_full.html">Full metric table</a></p></section>')
  html <- c(html,unlist(Map(function(label,path) paste0('<section><h2>',escape(label),'</h2><img src="',escape(path),'"></section>'),names(main),unname(main))))
  ocean_atomic_lines(html,file.path(report_dir,"preview.html"))
  saveRDS(list(diagnostic_plots=diagnostics$plots,
    surface_plots=if(!is.null(surfaces)) surfaces$plots else NULL,
    prediction_plots=if(!is.null(prediction)) prediction$plots else NULL,
    u_tables=ut,y_tables=yt),file.path(report_dir,"internal/presentation_plots.rds"))
  saveRDS(list(source_fit_root=output_dir,augmented_output_root=augmented_output_dir,
    report_dir=report_dir,copied_cache_files=copied,representative_fold=representative_fold,
    time=as.character(Sys.time()),no_refit=TRUE),file.path(report_dir,"internal/review_run.rds"))
  progress("COMPLETE")
  invisible(list(diagnostics=diagnostics,prediction=prediction,surfaces=surfaces,report_dir=report_dir))
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

adni_case_report <- function(root, folds_expected = 1:5, smoke = FALSE) {
  stopifnot(!anyDuplicated(folds_expected), all(folds_expected %in% 1:5))
  case_dirs <- file.path(root, paste0("fold_", folds_expected), "case-study")
  files <- file.path(case_dirs, "predictions.csv")
  if (!all(file.exists(files))) stop("Missing fold prediction files; reporting requires every requested fold.")
  z <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  status <- do.call(rbind, lapply(file.path(case_dirs, "method_status.csv"), read.csv,
                                stringsAsFactors = FALSE))
  adni_write_csv(status, file.path(root, "appendix", "method_status.csv"))
  stopifnot(!anyDuplicated(z[, c("RID", "method", "scenario")]))
  reference <- z[z$method == "EIV-GP" & z$scenario == "no_test_CSF", c("RID", "fold")]
  for (m in ADNI_METHODS[-1]) {
    a <- z[z$method == m & z$scenario == "no_test_CSF", c("RID", "fold")]
    for (f in unique(a$fold)) if (!setequal(a$RID[a$fold == f], reference$RID[reference$fold == f]))
      stop("Participant mismatch for ", m, " fold ", f)
  }
  if (!smoke) {
    assignment_file <- file.path(root, "appendix", "fold_assignments.csv")
    if (!file.exists(assignment_file)) stop("Missing saved fold assignments.")
    assignment <- read.csv(assignment_file)
    stopifnot(nrow(assignment) == 495L, !anyDuplicated(assignment$RID),
              all(table(factor(assignment$fold, levels = 1:5)) == 99L))
    expected <- assignment[assignment$fold %in% folds_expected, c("RID", "fold")]
    if (nrow(reference) != nrow(expected) || !setequal(reference$RID, expected$RID) ||
        any(reference$fold != expected$fold[match(reference$RID, expected$RID)]))
      stop("Predictions do not match the saved held-out assignments.")
  }
  coverage <- do.call(rbind, lapply(ADNI_METHODS, function(m) {
    a <- z[z$method == m & z$scenario == "no_test_CSF", ]
    complete_unit <- !smoke && nrow(a) == 495L && setequal(unique(a$fold), 1:5)
    data.frame(method = m, n_folds = length(unique(a$fold)), n_predictions = nrow(a),
      all_requested_folds = setequal(unique(a$fold), folds_expected),
      complete_cv = complete_unit,
      flagged_folds = sum(status$method == m & (status$optimization_status != "converged" &
        status$optimization_status != "mcmc_gate_passed"), na.rm = TRUE))
  }))
  adni_write_csv(coverage, file.path(root, "appendix", "completeness.csv"))
  # Main table requires all methods on the same fold set. Never compare each method's favorable subset.
  common_folds <- Reduce(intersect, lapply(ADNI_METHODS, function(m)
    unique(z$fold[z$method == m & z$scenario == "no_test_CSF"])))
  main <- list(); detailed <- list(); fold_metrics <- list()
  for (m in ADNI_METHODS) for (scenario in unique(z$scenario)) {
    a <- z[z$method == m & z$scenario == scenario, ]
    if (!nrow(a)) next
    metric <- adni_metrics(a); metric$method <- m; metric$scenario <- scenario
    metric$complete_cv <- coverage$complete_cv[coverage$method == m]
    detailed[[length(detailed) + 1L]] <- metric
    for (f in unique(a$fold)) {
      fm <- adni_metrics(a[a$fold == f, ]); fm$method <- m; fm$scenario <- scenario; fm$fold <- f
      fold_metrics[[length(fold_metrics) + 1L]] <- fm
    }
    if (scenario == "no_test_CSF") {
      matched <- a[a$fold %in% common_folds, ]
      if (nrow(matched)) {
        mm <- adni_metrics(matched)
        mm <- mm[mm$group %in% c("overall", "R0", "R1"), ]
        mm$method <- m; mm$matched_folds <- length(common_folds)
        mm$complete_comparison <- !smoke && all(coverage$complete_cv) &&
          setequal(common_folds, 1:5)
        mm$diagnostic_flag <- coverage$flagged_folds[coverage$method == m] > 0
        main[[length(main) + 1L]] <- mm
      }
    }
  }
  adni_write_csv(do.call(rbind, detailed), file.path(root, "appendix", "all_scenario_metrics.csv"))
  adni_write_csv(do.call(rbind, fold_metrics), file.path(root, "appendix", "fold_metrics.csv"))
  main_table <- if (length(main)) do.call(rbind, main) else data.frame(
    status = "No fold has valid predictions from all four methods; no comparative ranking available.")
  if ("method" %in% names(main_table)) main_table <- main_table[, c("method", "group", "n",
    "RMSE", "CRPS", "Coverage95", "Width95", "IntervalScore95", "complete_comparison", "diagnostic_flag")]
  adni_write_csv(main_table, file.path(root, "main", "table1_prediction.csv"))
  writeLines(c("# ADNI prediction comparison", "",
    if (smoke) "SMOKE OUTPUT - no inferential interpretation." else
      if (!all(coverage$complete_cv) || !setequal(common_folds, 1:5))
        "INCOMPLETE: descriptive common-fit subset only; no full-run superiority claim." else
        "All methods evaluated on the identical held-out participants.", "",
    "Primary task: PET response prediction with no test CSF. Overall, R=0, and R=1 scores are reported.",
    "EIV-GP additionally uses calibrated training CSF; competitors use training X,C,Y. This compares available workflows, not an isolated algorithmic effect.",
    "Diagnostic warnings are retained. Coverage and width are not ranked. See appendix for all failures and settings.", "",
    as.character(knitr::kable(main_table, format = "pipe", digits = 3))),
    file.path(root, "main", "table1_prediction.md"))

  # Per-participant paired contrasts. Bootstrap conditions on fitted folds, not independent refitted studies.
  paired <- list(); set.seed(202650000L)
  for (m in ADNI_METHODS[-1]) for (group in c("overall", "R0")) {
    e <- z[z$method == "EIV-GP" & z$scenario == "no_test_CSF", ]
    a <- z[z$method == m & z$scenario == "no_test_CSF", ]
    if (group == "R0") { e <- e[e$R == 0, ]; a <- a[a$R == 0, ] }
    merge_cols <- c("cv_folds", "fold", "RID")
    p <- merge(e[, c(merge_cols, "CRPS", "absolute_error")],
               a[, c(merge_cols, "CRPS", "absolute_error")], by = merge_cols, suffixes = c("_eiv", "_other"))
    if (!nrow(p)) next
    for (loss in c("CRPS", "absolute_error")) {
      difference <- p[[paste0(loss, "_other")]] - p[[paste0(loss, "_eiv")]]
      bs <- replicate(if (smoke) 50L else 1000L, mean(sample(difference, replace = TRUE)))
      paired[[length(paired) + 1L]] <- data.frame(method = m, group = group, loss = loss,
        n_pairs = length(difference), competitor_minus_EIV = mean(difference),
        lo_descriptive = quantile(bs, .025), hi_descriptive = quantile(bs, .975))
    }
  }
  if (length(paired)) adni_write_csv(do.call(rbind, paired), file.path(root, "appendix", "paired_losses.csv"))
  # Observable task shared by every method: held-out response predictions.
  zp <- z[z$scenario == "no_test_CSF", ]; zp$method <- factor(zp$method, levels = ADNI_METHODS)
  p <- ggplot2::ggplot(zp, ggplot2::aes(x = y_true, y = pred_mean)) +
    ggplot2::geom_point(alpha = .45, size = .8) + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::facet_wrap(~method, nrow = 1) + ggplot2::theme_bw() +
    ggplot2::labs(x = "Observed Centiloid", y = "Out-of-fold prediction", title = "Common response prediction task")
  ggplot2::ggsave(file.path(root, "appendix", "prediction_scatter.pdf"), p, width = 11, height = 3.3)
  uf <- file.path(case_dirs, "hidden_CSF_validation.csv"); uf <- uf[file.exists(uf)]
  if (length(uf)) {
    uz <- do.call(rbind, lapply(uf, read.csv, stringsAsFactors = FALSE))
    adni_write_csv(uz, file.path(root, "appendix", "hidden_CSF_predictions.csv"))
    um <- do.call(rbind, lapply(split(uz, uz$method), function(a) data.frame(method = a$method[1],
      n = nrow(a), RMSE = sqrt(mean(a$squared_error)), CRPS = mean(a$CRPS),
      Coverage95 = mean(a$covered95), Width95 = mean(a$Width95),
      mean_probability_nonpositive_CSF = mean(a$probability_nonpositive_CSF))))
    adni_write_csv(um, file.path(root, "appendix", "hidden_CSF_metrics.csv"))
    e <- uz[uz$method == "EIV-GP: prospective U|C", ]
    pi <- ggplot2::ggplot(e, ggplot2::aes(x = y_true, y = pred_mean)) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = pred_lo95, ymax = pred_hi95), alpha = .17) +
      ggplot2::geom_point(size = 1.2, color = "#0072B2") + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
      ggplot2::labs(x = "Held-out observed CSF A-beta42", y = "Prospective U|C mean and 95% interval",
                    title = "A  Masked-CSF validation (observed subset)") + ggplot2::theme_bw()
    sf <- file.path(root, "fold_1", "case-study", "CSF_PET_surface.csv")
    if (file.exists(sf)) {
      s <- read.csv(sf)
      ps <- ggplot2::ggplot(s, ggplot2::aes(x = U, y = mean)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lo95, ymax = hi95), fill = "#0072B2", alpha = .2) +
        ggplot2::geom_line(color = "#0072B2") + ggplot2::theme_bw() +
        ggplot2::labs(x = "CSF A-beta42 (source units)", y = "Latent PET response surface",
          title = "B  Physical-coordinate inference",
          subtitle = paste0("Prespecified fold 1; age ", round(s$age[1], 1), ", female=0, APOE4 dose=1"))
      combined <- patchwork::wrap_plots(pi, ps) + patchwork::plot_annotation(
        title = if (smoke) "SMOKE OUTPUT - not inferential results" else "EIV-GP: additional inference on the CSF scale",
        caption = paste("Pointwise posterior intervals; inspect convergence diagnostics.",
          "Surface is an association, not a causal effect. Masked-CSF validation does not verify naturally missing CSF."))
      ggplot2::ggsave(file.path(root, "main", "fig2_CSF_inference.pdf"), combined, width = 10, height = 4.5)
      ggplot2::ggsave(file.path(root, "main", "fig2_CSF_inference.png"), combined, width = 10, height = 4.5, dpi = 150)
    }
  }
  writeLines(c("# Reporting boundaries", "",
    "Do not claim superiority from smoke, incomplete, or inadequately mixed fits.",
    "Positive paired differences favor EIV-GP. Bootstrap intervals are descriptive conditional-on-fits intervals; they omit refitting and model-selection uncertainty.",
    "Each participant is tested once. Training sets overlap, so fold scores are not independent replications.",
    "The literature competitors do not estimate physical CSF posteriors or a CSF-coordinate surface. Other calibrated models can; these capabilities are not universally unique to EIV-GP.",
    "Background and task specifications are in the notebook. EDA is descriptive and must not drive post-hoc test-set tuning."),
    file.path(root, "appendix", "INTERPRETATION.md"))
  invisible(main_table)
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