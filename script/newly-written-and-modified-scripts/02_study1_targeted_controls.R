############################################################
## 02_study1_targeted_controls.R
##
## Two prespecified Study I controls at one calibration size:
##   1. category_sufficient: remove within-category response heterogeneity
##      while preserving E{f_active(x, U) | C};
##   2. balanced_threshold: retain the active response but use marginal
##      quantile cut points and no minimum-cell-count regeneration.
##
## The category-sufficient control uses the July 27 replication seeds, so its
## X, U, C, and response innovations are paired with the locked primary study.
############################################################

if (!exists("fit_eivgp_1d")) source("model_univariate.R")
if (!exists("run_study1_published_competitors")) {
  source("competitors.R")
}

needed_pkgs <- c("dplyr", "tidyr", "knitr")
missing_pkgs <- needed_pkgs[
  !vapply(needed_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0L) {
  stop("Please install required packages: ", paste(missing_pkgs, collapse = ", "))
}
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(knitr)
})

if (!exists("STUDY1_CONFIG")) STUDY1_CONFIG <- "quick"
if (!STUDY1_CONFIG %in% c("quick", "balanced", "thorough")) {
  stop("STUDY1_CONFIG must be 'quick', 'balanced', or 'thorough'.")
}
if (!exists("STUDY1_USE_CACHE")) STUDY1_USE_CACHE <- TRUE
if (!exists("STUDY1_OUT_PREFIX")) STUDY1_OUT_PREFIX <- ".."
if (!exists("STUDY1_STRICT_COMPETITORS")) {
  STUDY1_STRICT_COMPETITORS <- identical(STUDY1_CONFIG, "thorough")
}
if (!exists("STUDY1_PUBLISHED_COMPETITORS")) {
  STUDY1_PUBLISHED_COMPETITORS <- c("UC-GP", "LVGP", "EzGP")
}
if (!exists("STUDY1_CONTROL_CALIB")) STUDY1_CONTROL_CALIB <- 20L

settings <- switch(
  STUDY1_CONFIG,
  quick = list(
    n_test = 150L, n_rep = 2L, n_iter = 600L, burn = 200L,
    n_chains = 1L, preset = "fast", n_pred_draw = 150L
  ),
  balanced = list(
    n_test = 300L, n_rep = 10L, n_iter = 2500L, burn = 750L,
    n_chains = 4L, preset = "balanced", n_pred_draw = 300L
  ),
  thorough = list(
    n_test = 500L, n_rep = 50L, n_iter = 5000L, burn = 1000L,
    n_chains = 12L, preset = "balanced", n_pred_draw = 600L
  )
)

n_train <- 100L
m <- 6L
n_calib <- as.integer(STUDY1_CONTROL_CALIB)
if (length(n_calib) != 1L || is.na(n_calib) || n_calib <= 0L ||
    n_calib > n_train) {
  stop("STUDY1_CONTROL_CALIB must be between 1 and n_train.")
}

control_specs <- list(
  category_sufficient = list(
    response_scenario = "category_sufficient",
    threshold_design = "imbalanced",
    min_class_count = 3L,
    label = "Category-sufficient control"
  ),
  balanced_threshold = list(
    response_scenario = "active",
    threshold_design = "balanced",
    min_class_count = 0L,
    label = "Balanced-threshold sensitivity"
  )
)

FIG_DIR <- file.path(STUDY1_OUT_PREFIX, "figures")
TAB_DIR <- file.path(STUDY1_OUT_PREFIX, "tables")
RES_DIR <- file.path(STUDY1_OUT_PREFIX, "results", "study1_targeted_controls")
REP_DIR <- file.path(RES_DIR, "replications")
for (dd in c(FIG_DIR, TAB_DIR, RES_DIR, REP_DIR)) {
  dir.create(dd, showWarnings = FALSE, recursive = TRUE)
}

CONTROL_TAG <- paste0(
  "study1_targeted_controls_v1_", STUDY1_CONFIG,
  "_k", n_calib,
  "_", paste(STUDY1_PUBLISHED_COMPETITORS, collapse = "-")
)
study1_parallel_level <- if (exists("STUDY1_PARALLEL_LEVEL")) {
  match.arg(STUDY1_PARALLEL_LEVEL, c("chains", "replications", "none"))
} else {
  "chains"
}
parallel_chains <- mixedgp_parallel_chains_enabled(study1_parallel_level)

method_levels <- c("Oracle", "EIV-GP", "UC-GP", "LVGP", "EzGP")
metric_levels <- c(
  "RMSE", "MAE", "CRPS", "Coverage95", "Width95", "IntervalScore95"
)
safe_se <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}
format_mean_se <- function(mean, se, digits = 3L) {
  if (!is.finite(mean)) return("--")
  if (!is.finite(se)) return(sprintf(paste0("%.", digits, "f"), mean))
  paste0(
    sprintf(paste0("%.", digits, "f"), mean),
    " (", sprintf(paste0("%.", digits, "f"), se), ")"
  )
}

if (!exists("STUDY1_LVGP_MAX_ELAPSED")) {
  STUDY1_LVGP_MAX_ELAPSED <- if (STUDY1_CONFIG == "quick") 180 else 900
}
competitor_controls <- list(
  `UC-GP` = list(n_starts = if (STUDY1_CONFIG == "quick") 2L else 8L),
  LVGP = list(
    n_starts = if (STUDY1_CONFIG == "quick") 2L else 8L,
    max_retries = if (STUDY1_CONFIG == "quick") 1L else 3L,
    max_iter_ini = if (STUDY1_CONFIG == "quick") 30L else 100L,
    max_iter_lat = if (STUDY1_CONFIG == "quick") 8L else 20L,
    rescue_iter_ini = 300L,
    rescue_iter_lat = 100L,
    max_elapsed_seconds = STUDY1_LVGP_MAX_ELAPSED,
    parallel = FALSE
  ),
  EzGP = list(
    tau_fractions = c(1e-6, 0.0025, 0.01, 0.04, 0.16),
    cv_folds = 3L,
    maxeval = if (STUDY1_CONFIG == "quick") 30L else 100L
  )
)

run_one_control <- function(rep_id, control_name) {
  spec <- control_specs[[control_name]]
  scenario_id <- match(control_name, names(control_specs))
  data_seed <- 100000L + rep_id
  fit_seed_base <- 1000000L * scenario_id + 10000L * rep_id

  dat <- simulate_1d_data(
    n = n_train,
    n_test = settings$n_test,
    m = m,
    scenario = spec$response_scenario,
    seed = data_seed,
    threshold_design = spec$threshold_design,
    min_class_count = spec$min_class_count
  )
  train <- dat$train
  test <- dat$test
  calib_idx <- make_nested_calibration_sets(
    n_train,
    calib_grid = n_calib,
    seed = 200000L + rep_id
  )[[as.character(n_calib)]]

  set.seed(fit_seed_base + 10L)
  oracle_draws <- sample_oracle_test_y(
    x_test = test$x,
    c_test = test$c,
    tau_true = dat$tau_true,
    scenario = spec$response_scenario,
    sigma_eps = dat$sigma_eps,
    n_draw = settings$n_pred_draw
  )
  metrics <- list(
    Oracle = summarize_predictive_samples_1d(
      oracle_draws, test$y, "Oracle", rep_id, n_calib, control_name
    )
  )

  competitors <- run_study1_published_competitors(
    X_train = matrix(train$x, ncol = 1L),
    y_train = train$y,
    C_train = matrix(train$c, ncol = 1L),
    X_test = matrix(test$x, ncol = 1L),
    C_test = matrix(test$c, ncol = 1L),
    m_vec = m,
    n_draw = settings$n_pred_draw,
    seed = fit_seed_base + 100L,
    methods = STUDY1_PUBLISHED_COMPETITORS,
    strict = STUDY1_STRICT_COMPETITORS,
    controls = competitor_controls
  )
  for (method in names(competitors$draws)) {
    metrics[[method]] <- summarize_predictive_samples_1d(
      competitors$draws[[method]], test$y, method, rep_id, n_calib,
      control_name
    )
  }

  fit_eiv <- fit_eivgp_1d(
    x_raw = train$x,
    y_raw = train$y,
    c_ord = train$c,
    u_true = train$u,
    calib_idx = calib_idx,
    m = m,
    tau_true = dat$tau_true,
    n_iter = settings$n_iter,
    burn = settings$burn,
    thin = 1L,
    n_chains = settings$n_chains,
    preset = settings$preset,
    seed = fit_seed_base + 1000L,
    parallel_chains = parallel_chains,
    verbose = FALSE
  )
  draw_ids <- seq_len(nrow(fit_eiv$mcmc$samples_u))
  if (length(draw_ids) > settings$n_pred_draw) {
    set.seed(fit_seed_base + 2000L)
    draw_ids <- sample(draw_ids, settings$n_pred_draw)
  }
  eiv_draws <- sample_eiv_test_y(
    x_test_raw = test$x,
    c_test = test$c,
    fit_obj = fit_eiv,
    draw_ids = draw_ids,
    n_per_draw = 1L
  )
  metrics[["EIV-GP"]] <- summarize_predictive_samples_1d(
    eiv_draws, test$y, "EIV-GP", rep_id, n_calib, control_name
  )

  status <- competitors$status
  status$rep <- rep_id
  status$scenario <- control_name
  list(
    metrics = bind_rows(metrics),
    status = status,
    metadata = list(
      rep = rep_id,
      scenario = control_name,
      data_seed = data_seed,
      calibration_seed = 200000L + rep_id,
      fit_seed_base = fit_seed_base,
      tau_true = dat$tau_true,
      min_class_count = spec$min_class_count
    )
  )
}

run_grid <- expand.grid(
  scenario = names(control_specs),
  rep = seq_len(settings$n_rep),
  stringsAsFactors = FALSE
)
rep_objects <- vector("list", nrow(run_grid))
for (ii in seq_len(nrow(run_grid))) {
  scenario <- run_grid$scenario[ii]
  rep_id <- run_grid$rep[ii]
  rep_file <- file.path(
    REP_DIR,
    sprintf("%s_%s_rep%03d.rds", CONTROL_TAG, scenario, rep_id)
  )
  if (isTRUE(STUDY1_USE_CACHE) && file.exists(rep_file)) {
    rep_objects[[ii]] <- readRDS(rep_file)
  } else {
    message(
      "Study I targeted control ", scenario, ", replication ", rep_id,
      " of ", settings$n_rep
    )
    rep_objects[[ii]] <- run_one_control(rep_id, scenario)
    saveRDS(rep_objects[[ii]], rep_file)
  }
}

control_results <- bind_rows(lapply(rep_objects, `[[`, "metrics"))
control_status <- bind_rows(lapply(rep_objects, `[[`, "status"))
write.csv(
  control_results,
  file.path(TAB_DIR, paste0("study1_targeted_controls_raw_", CONTROL_TAG, ".csv")),
  row.names = FALSE
)
write.csv(
  control_status,
  file.path(TAB_DIR, paste0("study1_targeted_controls_status_", CONTROL_TAG, ".csv")),
  row.names = FALSE
)

control_summary <- control_results |>
  pivot_longer(
    cols = all_of(metric_levels), names_to = "metric", values_to = "value"
  ) |>
  group_by(scenario, method, metric) |>
  summarise(
    mean = mean(value, na.rm = TRUE),
    se = safe_se(value),
    .groups = "drop"
  ) |>
  mutate(value = mapply(format_mean_se, mean, se, USE.NAMES = FALSE)) |>
  select(scenario, method, metric, value) |>
  pivot_wider(names_from = metric, values_from = value) |>
  mutate(
    Setting = vapply(scenario, function(z) control_specs[[z]]$label, character(1)),
    Method = factor(method, levels = method_levels)
  ) |>
  arrange(factor(scenario, levels = names(control_specs)), Method) |>
  select(Setting, Method, RMSE, CRPS, Coverage95, Width95, IntervalScore95)
writeLines(
  knitr::kable(
    control_summary,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE
  ),
  file.path(TAB_DIR, "study1_targeted_controls.tex")
)

control_advantage <- control_results |>
  filter(method == "EIV-GP") |>
  select(rep, scenario, CRPS_eiv = CRPS,
         IntervalScore95_eiv = IntervalScore95) |>
  inner_join(
    control_results |>
      filter(method %in% STUDY1_PUBLISHED_COMPETITORS) |>
      select(rep, scenario, competitor = method,
             CRPS_comp = CRPS,
             IntervalScore95_comp = IntervalScore95),
    by = c("rep", "scenario")
  ) |>
  mutate(
    CRPS_advantage = CRPS_comp - CRPS_eiv,
    IntervalScore95_advantage =
      IntervalScore95_comp - IntervalScore95_eiv
  )
write.csv(
  control_advantage,
  file.path(TAB_DIR, paste0("study1_targeted_advantage_raw_", CONTROL_TAG, ".csv")),
  row.names = FALSE
)

active_file <- file.path(TAB_DIR, "study1_mc_raw_results_revised.csv")
if (file.exists(active_file)) {
  active_results <- read.csv(active_file, stringsAsFactors = FALSE)
  active_advantage <- active_results |>
    filter(method == "EIV-GP", n_calib == .env$n_calib) |>
    select(rep, CRPS_eiv = CRPS,
           IntervalScore95_eiv = IntervalScore95) |>
    inner_join(
      active_results |>
        filter(method %in% STUDY1_PUBLISHED_COMPETITORS) |>
        select(rep, competitor = method, CRPS_comp = CRPS,
               IntervalScore95_comp = IntervalScore95),
      by = "rep"
    ) |>
    mutate(
      active_CRPS_advantage = CRPS_comp - CRPS_eiv,
      active_IntervalScore95_advantage =
        IntervalScore95_comp - IntervalScore95_eiv
    ) |>
    select(rep, competitor, starts_with("active_"))

  heterogeneity_contrast <- control_advantage |>
    filter(scenario == "category_sufficient") |>
    inner_join(active_advantage, by = c("rep", "competitor")) |>
    mutate(
      CRPS_advantage_increase =
        active_CRPS_advantage - CRPS_advantage,
      IntervalScore95_advantage_increase =
        active_IntervalScore95_advantage - IntervalScore95_advantage
    )
  write.csv(
    heterogeneity_contrast,
    file.path(
      TAB_DIR,
      paste0("study1_heterogeneity_contrast_raw_", CONTROL_TAG, ".csv")
    ),
    row.names = FALSE
  )

  heterogeneity_table <- heterogeneity_contrast |>
    group_by(competitor) |>
    summarise(
      `Increase in EIV CRPS advantage` = format_mean_se(
        mean(CRPS_advantage_increase), safe_se(CRPS_advantage_increase)
      ),
      `Increase in EIV interval-score advantage` = format_mean_se(
        mean(IntervalScore95_advantage_increase),
        safe_se(IntervalScore95_advantage_increase)
      ),
      Pairs = sum(is.finite(CRPS_advantage_increase)),
      .groups = "drop"
    ) |>
    rename(Competitor = competitor)
  writeLines(
    knitr::kable(
      heterogeneity_table,
      format = "latex",
      booktabs = TRUE,
      escape = FALSE
    ),
    file.path(TAB_DIR, "study1_heterogeneity_advantage_contrast.tex")
  )
} else {
  warning(
    "Primary Study I revised results were not found. The targeted-control table ",
    "was written, but the paired heterogeneity contrast was not computed.",
    call. = FALSE
  )
}

design_manifest <- bind_rows(lapply(names(control_specs), function(scenario) {
  spec <- control_specs[[scenario]]
  data.frame(
    design_tag = CONTROL_TAG,
    scenario = scenario,
    response_scenario = spec$response_scenario,
    threshold_design = spec$threshold_design,
    min_class_count = spec$min_class_count,
    n_train = n_train,
    n_test = settings$n_test,
    n_rep = settings$n_rep,
    n_calib = n_calib,
    data_seed_rule = "100000 + replication",
    calibration_seed_rule = "200000 + replication",
    primary_target = "Y_star_given_X_star_C_star",
    stringsAsFactors = FALSE
  )
}))
write.csv(
  design_manifest,
  file.path(TAB_DIR, paste0("study1_targeted_control_manifest_", CONTROL_TAG, ".csv")),
  row.names = FALSE
)
capture.output(
  sessionInfo(),
  file = file.path(RES_DIR, paste0("sessionInfo_", CONTROL_TAG, ".txt"))
)

message("Study I targeted-control outputs written under ", normalizePath(STUDY1_OUT_PREFIX))
