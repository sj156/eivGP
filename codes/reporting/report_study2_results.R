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
