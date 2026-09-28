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
