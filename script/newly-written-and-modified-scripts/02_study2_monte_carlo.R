############################################################
## 02_study2_monte_carlo.R
##
## Publication Monte Carlo driver for Study II.
##
## Canonical design:
##   * n = 120 training observations;
##   * four ordinal proxies of a two-dimensional latent state;
##   * primary, latent-additive, high-uncertainty, and logistic
##     misspecification scenarios;
##   * primary recovery of m(x,c)=E{f(x,U)|C=c};
##   * prediction of Y* | X*, C* with no test U*;
##   * separate latent-imputation and physically anchored f(x, u) summaries;
##   * published competitors only, with explicit availability records.
############################################################

if (!exists("fit_eivgp_ordprobit_fb")) {
  source("00_study2_functions.R")
}
if (!exists("load_mixedgp_synthetic_dataset")) source("00_synthetic_data.R")
if (!exists("run_study2_published_competitors")) {
  source("03_study2_published_competitors.R")
}
if (!exists("STUDY2_RUN_ABLATIONS")) STUDY2_RUN_ABLATIONS <- TRUE
if (!exists("STUDY2_EVALUATE_F")) STUDY2_EVALUATE_F <- TRUE
if (!exists("STUDY2_EVALUATE_U")) STUDY2_EVALUATE_U <- TRUE
if (isTRUE(STUDY2_RUN_ABLATIONS) &&
    !exists("fit_ordinalprobit_measurement_fb")) {
  source("04_study2_ablations.R")
}

needed_pkgs <- c("ggplot2", "dplyr", "tidyr", "knitr")
if ((exists("STUDY2_CONFIG") && !identical(STUDY2_CONFIG, "quick")) ||
    isTRUE(STUDY2_RUN_ABLATIONS)) {
  needed_pkgs <- c(needed_pkgs, "posterior")
}
missing_pkgs <- needed_pkgs[
  !vapply(needed_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0L) {
  stop("Please install required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(knitr)
})

if (!exists("STUDY2_CONFIG")) STUDY2_CONFIG <- "quick"
if (!exists("STUDY2_USE_CACHE")) STUDY2_USE_CACHE <- TRUE
if (!exists("STUDY2_MC_RESUME")) STUDY2_MC_RESUME <- TRUE
if (!exists("STUDY2_OUT_PREFIX")) STUDY2_OUT_PREFIX <- ".."
if (!exists("STUDY2_SAVE_REP_FITS")) STUDY2_SAVE_REP_FITS <- FALSE
if (!exists("STUDY2_STRICT_COMPETITORS")) STUDY2_STRICT_COMPETITORS <- FALSE
if (!exists("STUDY2_ENFORCE_MCMC_GATE")) {
  STUDY2_ENFORCE_MCMC_GATE <- !identical(STUDY2_CONFIG, "quick")
}
if (!exists("STUDY2_SAVE_PDF")) STUDY2_SAVE_PDF <- TRUE
if (!exists("STUDY2_SAVE_PNG")) STUDY2_SAVE_PNG <- TRUE
if (!exists("STUDY2_ABLATION_SCENARIOS")) {
  STUDY2_ABLATION_SCENARIOS <- "primary"
}

if (!exists("STUDY2_SCENARIOS")) {
  STUDY2_SCENARIOS <- c(
    "primary",
    "latent_additive_control",
    "high_uncertainty",
    "logistic_misspec"
  )
}

if (!exists("STUDY2_PUBLISHED_COMPETITORS")) {
  STUDY2_PUBLISHED_COMPETITORS <- c("UC-GP", "LVGP", "EzGP")
}

if (!exists("STUDY2_PRIMARY_CALIB_GRID")) {
  STUDY2_PRIMARY_CALIB_GRID <- c(0L, 10L, 25L, 50L, 80L)
}
if (!exists("STUDY2_CONTRAST_CALIB")) STUDY2_CONTRAST_CALIB <- 25L
if (!exists("STUDY2_MC_N_REP")) STUDY2_MC_N_REP <- NULL
if (!exists("STUDY2_DATA_DIR")) {
  STUDY2_DATA_DIR <- file.path("..", "data-synthetic", "study2")
}

settings <- study2_config_settings(STUDY2_CONFIG)
if (!exists("STUDY2_LVGP_MAX_ELAPSED")) {
  STUDY2_LVGP_MAX_ELAPSED <- if (STUDY2_CONFIG == "quick") 180 else 1800
}

study2_competitor_controls <- list(
  `UC-GP` = list(n_starts = if (STUDY2_CONFIG == "quick") 2L else 8L),
  LVGP = list(
    n_starts = if (STUDY2_CONFIG == "quick") 2L else 8L,
    max_retries = if (STUDY2_CONFIG == "quick") 1L else 3L,
    max_iter_ini = if (STUDY2_CONFIG == "quick") 30L else 100L,
    max_iter_lat = if (STUDY2_CONFIG == "quick") 8L else 20L,
    rescue_iter_ini = 300L,
    rescue_iter_lat = 100L,
    max_elapsed_seconds = STUDY2_LVGP_MAX_ELAPSED,
    parallel = FALSE
  ),
  EzGP = list(
    tau_fractions = c(1e-6, 0.0025, 0.01, 0.04, 0.16),
    cv_folds = 3L,
    maxeval = if (STUDY2_CONFIG == "quick") 30L else 100L
  )
)

allowed_scenarios <- c(
  "primary",
  "latent_additive_control",
  "high_uncertainty",
  "logistic_misspec"
)
unknown_scenarios <- setdiff(STUDY2_SCENARIOS, allowed_scenarios)
if (length(STUDY2_SCENARIOS) == 0L || length(unknown_scenarios) > 0L) {
  stop(
    "STUDY2_SCENARIOS must be a nonempty subset of: ",
    paste(allowed_scenarios, collapse = ", "),
    if (length(unknown_scenarios) > 0L) {
      paste0(". Unknown: ", paste(unknown_scenarios, collapse = ", "))
    } else {
      ""
    }
  )
}
unknown_ablation_scenarios <- setdiff(
  STUDY2_ABLATION_SCENARIOS,
  allowed_scenarios
)
if (length(unknown_ablation_scenarios) > 0L) {
  stop(
    "Unknown STUDY2_ABLATION_SCENARIOS: ",
    paste(unknown_ablation_scenarios, collapse = ", ")
  )
}

############################################################
## Frozen settings and output locations
############################################################

n_train <- 120L
n_train <- as.integer(get0(
  "STUDY2_MC_N_TRAIN", inherits = TRUE, ifnotfound = n_train
))
n_test <- as.integer(get0(
  "STUDY2_MC_N_TEST", inherits = TRUE, ifnotfound = settings$n_test
))
n_rep <- if (is.null(STUDY2_MC_N_REP)) settings$n_rep else as.integer(STUDY2_MC_N_REP)
if (anyNA(c(n_train, n_test, n_rep)) ||
    any(c(n_train, n_test, n_rep) < 1L)) {
  stop("Study II n_train, n_test, and n_rep must be positive integers.")
}
study2_q <- as.integer(get0("STUDY2_Q", inherits = TRUE, ifnotfound = 4L))
study2_m <- as.integer(get0("STUDY2_M", inherits = TRUE, ifnotfound = 4L))
if (length(study2_q) != 1L || is.na(study2_q) ||
    study2_q < 2L || study2_q > 6L ||
    length(study2_m) != 1L || is.na(study2_m) || study2_m < 2L) {
  stop("STUDY2_Q must be in 2:6 and STUDY2_M must be at least two.")
}
m_vec <- rep(study2_m, study2_q)
d_latent <- 2L
ident_method <- "lower_triangular"
STUDY2_DESIGN_TAG <- paste(
  "study2-publication-v14-alltasks-rfimp",
  paste0("q", study2_q, "-d", d_latent),
  "A-fixed-nested",
  paste0("balanced-", study2_m, "-level"),
  "sigma0.12",
  "exact-minimax-interwoven",
  "common-random-numbers",
  sep = "_"
)

mc_n_iter <- if (exists("STUDY2_MC_N_ITER")) {
  as.integer(STUDY2_MC_N_ITER)
} else settings$mc_n_iter
mc_burn <- if (exists("STUDY2_MC_BURN")) {
  as.integer(STUDY2_MC_BURN)
} else settings$mc_burn
mc_thin <- if (exists("STUDY2_MC_THIN")) {
  as.integer(STUDY2_MC_THIN)
} else settings$mc_thin
mc_n_chains <- if (exists("STUDY2_MC_N_CHAINS")) {
  as.integer(STUDY2_MC_N_CHAINS)
} else settings$mc_n_chains
mc_preset <- settings$preset
measurement_n_iter <- if (exists("STUDY2_MEAS_N_ITER")) {
  as.integer(STUDY2_MEAS_N_ITER)
} else mc_n_iter
measurement_burn <- if (exists("STUDY2_MEAS_BURN")) {
  as.integer(STUDY2_MEAS_BURN)
} else mc_burn
measurement_thin <- if (exists("STUDY2_MEAS_THIN")) {
  as.integer(STUDY2_MEAS_THIN)
} else mc_thin
measurement_n_chains <- if (exists("STUDY2_MEAS_N_CHAINS")) {
  as.integer(STUDY2_MEAS_N_CHAINS)
} else max(2L, mc_n_chains)
n_pred_draw <- as.integer(get0(
  "STUDY2_MC_N_PRED_DRAW", inherits = TRUE,
  ifnotfound = settings$n_pred_draw
))
predictive_latent_sampler <- if (exists("STUDY2_PREDICTIVE_LATENT_SAMPLER")) {
  as.character(STUDY2_PREDICTIVE_LATENT_SAMPLER)
} else settings$predictive_latent_sampler
if (length(predictive_latent_sampler) != 1L ||
    is.na(predictive_latent_sampler) ||
    !predictive_latent_sampler %in%
      c("minimax_tilting", "rejection", "gibbs")) {
  stop(
    "STUDY2_PREDICTIVE_LATENT_SAMPLER must be minimax_tilting, rejection, ",
    "or gibbs."
  )
}
diagnostic_n_new_latent_gibbs <- if (exists("STUDY2_DIAGNOSTIC_GIBBS_SWEEPS")) {
  as.integer(STUDY2_DIAGNOSTIC_GIBBS_SWEEPS)
} else settings$diagnostic_n_new_latent_gibbs
if (length(diagnostic_n_new_latent_gibbs) != 1L ||
    is.na(diagnostic_n_new_latent_gibbs) ||
    diagnostic_n_new_latent_gibbs < 1L) {
  stop("STUDY2_DIAGNOSTIC_GIBBS_SWEEPS must be a positive integer.")
}
rejection_max_batches <- settings$rejection_max_batches
n_oracle_pool <- as.integer(get0(
  "STUDY2_MC_N_ORACLE_POOL", inherits = TRUE,
  ifnotfound = settings$n_oracle_pool
))
n_m_eval <- if (STUDY2_CONFIG == "quick") 30L else if (
  STUDY2_CONFIG == "balanced"
) 60L else 80L
n_m_draw <- if (STUDY2_CONFIG == "quick") 50L else if (
  STUDY2_CONFIG == "balanced"
) 120L else 200L
n_m_latent <- if (STUDY2_CONFIG == "quick") 64L else if (
  STUDY2_CONFIG == "balanced"
) 128L else 256L
n_m_truth <- if (STUDY2_CONFIG == "quick") 1000L else 2000L
n_m_eval <- as.integer(get0(
  "STUDY2_MC_N_M_EVAL", inherits = TRUE, ifnotfound = n_m_eval
))
n_m_draw <- as.integer(get0(
  "STUDY2_MC_N_M_DRAW", inherits = TRUE, ifnotfound = n_m_draw
))
n_m_latent <- as.integer(get0(
  "STUDY2_MC_N_M_LATENT", inherits = TRUE, ifnotfound = n_m_latent
))
n_m_truth <- as.integer(get0(
  "STUDY2_MC_N_M_TRUTH", inherits = TRUE, ifnotfound = n_m_truth
))
if (anyNA(c(
  n_pred_draw, n_m_eval, n_m_draw, n_m_latent, n_m_truth, n_oracle_pool
)) || any(c(
  n_pred_draw, n_m_eval, n_m_draw, n_m_latent, n_m_truth, n_oracle_pool
) < 1L)) {
  stop("Study II prediction and mean-integration sizes must be positive integers.")
}
if (!exists("STUDY2_MCMC_RHAT_LIMIT")) STUDY2_MCMC_RHAT_LIMIT <- 1.05
if (!exists("STUDY2_MCMC_RAW_ESS_LIMIT")) STUDY2_MCMC_RAW_ESS_LIMIT <- 200
if (!exists("STUDY2_MCMC_TARGET_BULK_ESS_LIMIT")) {
  STUDY2_MCMC_TARGET_BULK_ESS_LIMIT <- 50
}
if (!exists("STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT")) {
  STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT <- 50
}
if (!exists("STUDY2_MCMC_TARGET_N_POINTS")) {
  STUDY2_MCMC_TARGET_N_POINTS <- 5L
}
if (!exists("STUDY2_MEAS_RHAT_LIMIT")) {
  STUDY2_MEAS_RHAT_LIMIT <- STUDY2_MCMC_RHAT_LIMIT
}
if (!exists("STUDY2_MEAS_BULK_ESS_LIMIT")) {
  STUDY2_MEAS_BULK_ESS_LIMIT <- STUDY2_MCMC_TARGET_BULK_ESS_LIMIT
}
if (!exists("STUDY2_MEAS_TAIL_ESS_LIMIT")) {
  STUDY2_MEAS_TAIL_ESS_LIMIT <- STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT
}
if (!exists("STUDY2_ABLATION_GP_N_STARTS")) {
  STUDY2_ABLATION_GP_N_STARTS <- if (STUDY2_CONFIG == "quick") 2L else 8L
}
if (!exists("STUDY2_ABLATION_GP_MAXIT")) {
  STUDY2_ABLATION_GP_MAXIT <- 500L
}
if (!is.finite(STUDY2_MCMC_RHAT_LIMIT) ||
    STUDY2_MCMC_RHAT_LIMIT <= 1 ||
    !is.finite(STUDY2_MCMC_RAW_ESS_LIMIT) ||
    STUDY2_MCMC_RAW_ESS_LIMIT < 1 ||
    !is.finite(STUDY2_MCMC_TARGET_BULK_ESS_LIMIT) ||
    STUDY2_MCMC_TARGET_BULK_ESS_LIMIT < 1 ||
    !is.finite(STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT) ||
    STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT < 1 ||
    length(STUDY2_MCMC_TARGET_N_POINTS) != 1L ||
    is.na(STUDY2_MCMC_TARGET_N_POINTS) ||
    STUDY2_MCMC_TARGET_N_POINTS < 1L ||
    !is.finite(STUDY2_MEAS_RHAT_LIMIT) || STUDY2_MEAS_RHAT_LIMIT <= 1 ||
    !is.finite(STUDY2_MEAS_BULK_ESS_LIMIT) ||
    STUDY2_MEAS_BULK_ESS_LIMIT < 1 ||
    !is.finite(STUDY2_MEAS_TAIL_ESS_LIMIT) ||
    STUDY2_MEAS_TAIL_ESS_LIMIT < 1 ||
    length(STUDY2_ABLATION_GP_N_STARTS) != 1L ||
    is.na(STUDY2_ABLATION_GP_N_STARTS) ||
    STUDY2_ABLATION_GP_N_STARTS < 1L ||
    length(STUDY2_ABLATION_GP_MAXIT) != 1L ||
    is.na(STUDY2_ABLATION_GP_MAXIT) ||
    STUDY2_ABLATION_GP_MAXIT < 0L ||
    length(measurement_n_iter) != 1L || is.na(measurement_n_iter) ||
    measurement_n_iter < 1L ||
    length(measurement_burn) != 1L || is.na(measurement_burn) ||
    measurement_burn < 0L || measurement_burn >= measurement_n_iter ||
    length(measurement_thin) != 1L || is.na(measurement_thin) ||
    measurement_thin < 1L ||
    length(measurement_n_chains) != 1L || is.na(measurement_n_chains) ||
    (isTRUE(STUDY2_RUN_ABLATIONS) && measurement_n_chains < 2L)) {
  stop("Invalid Study II publication MCMC-gate settings.")
}
measurement_n_iter <- as.integer(measurement_n_iter)
measurement_burn <- as.integer(measurement_burn)
measurement_thin <- as.integer(measurement_thin)
measurement_n_chains <- as.integer(measurement_n_chains)
STUDY2_MCMC_TARGET_N_POINTS <- as.integer(STUDY2_MCMC_TARGET_N_POINTS)
STUDY2_ABLATION_GP_N_STARTS <- as.integer(STUDY2_ABLATION_GP_N_STARTS)
STUDY2_ABLATION_GP_MAXIT <- as.integer(STUDY2_ABLATION_GP_MAXIT)

if (!exists("STUDY2_PARALLEL_LEVEL")) STUDY2_PARALLEL_LEVEL <- "chains"
if (!STUDY2_PARALLEL_LEVEL %in% c("chains", "replications", "none")) {
  stop("STUDY2_PARALLEL_LEVEL must be chains, replications, or none.")
}
parallel_chains <- identical(STUDY2_PARALLEL_LEVEL, "chains")

scenario_calib_grid <- function(scenario) {
  if (identical(scenario, "primary")) {
    return(as.integer(STUDY2_PRIMARY_CALIB_GRID))
  }
  as.integer(STUDY2_CONTRAST_CALIB)
}

scenario_label <- function(scenario) {
  c(
    primary = "Primary interactive",
    latent_additive_control = "Latent-additive control",
    high_uncertainty = "High latent uncertainty",
    logistic_misspec = "Logistic-error robustness"
  )[[scenario]]
}

FIG_DIR <- file.path(STUDY2_OUT_PREFIX, "figures")
TAB_DIR <- file.path(STUDY2_OUT_PREFIX, "tables")
RES_DIR <- file.path(STUDY2_OUT_PREFIX, "results", "study2_manuscript_v5")
REP_DIR <- file.path(RES_DIR, "mc_replications")
FIT_DIR <- file.path(RES_DIR, "mc_fits")

for (dd in c(FIG_DIR, TAB_DIR, RES_DIR, REP_DIR, FIT_DIR)) {
  dir.create(dd, showWarnings = FALSE, recursive = TRUE)
}

scenario_code <- c(
  primary = "P",
  latent_additive_control = "A",
  high_uncertainty = "H",
  logistic_misspec = "L"
)
scenario_tag <- paste(unname(scenario_code[STUDY2_SCENARIOS]), collapse = "")
CACHE_SPEC <- list(
  schema = "s2v15_nlpd_alltasks_rfimp",
  design_tag = STUDY2_DESIGN_TAG,
  study2_config = STUDY2_CONFIG,
  scenario_code = scenario_tag,
  scenarios = STUDY2_SCENARIOS,
  ablation_scenarios = STUDY2_ABLATION_SCENARIOS,
  n_train = n_train,
  n_test = n_test,
  n_rep = n_rep,
  q = length(m_vec),
  m = m_vec[1L],
  primary_calibration_grid = STUDY2_PRIMARY_CALIB_GRID,
  contrast_calibration = STUDY2_CONTRAST_CALIB,
  mcmc = list(
    n_iter = mc_n_iter, burn = mc_burn, thin = mc_thin,
    n_chains = mc_n_chains,
    rhat_limit = STUDY2_MCMC_RHAT_LIMIT,
    raw_ess_limit = STUDY2_MCMC_RAW_ESS_LIMIT,
    target_bulk_ess_limit = STUDY2_MCMC_TARGET_BULK_ESS_LIMIT,
    target_tail_ess_limit = STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT,
    target_n_points = STUDY2_MCMC_TARGET_N_POINTS
  ),
  prediction = list(
    n_draw = n_pred_draw,
    latent_sampler = predictive_latent_sampler,
    oracle_pool = n_oracle_pool
  ),
  mean_recovery = list(
    n_eval = n_m_eval, n_draw = n_m_draw,
    n_latent = n_m_latent, n_truth = n_m_truth
  ),
  measurement = list(
    n_iter = measurement_n_iter, burn = measurement_burn,
    thin = measurement_thin, n_chains = measurement_n_chains,
    rhat_limit = STUDY2_MEAS_RHAT_LIMIT,
    bulk_ess_limit = STUDY2_MEAS_BULK_ESS_LIMIT,
    tail_ess_limit = STUDY2_MEAS_TAIL_ESS_LIMIT
  ),
  ablation_gp = list(
    n_starts = STUDY2_ABLATION_GP_N_STARTS,
    maxit = STUDY2_ABLATION_GP_MAXIT
  ),
  task_flags = list(
    run_ablations = isTRUE(STUDY2_RUN_ABLATIONS),
    evaluate_f = isTRUE(STUDY2_EVALUATE_F),
    evaluate_u = isTRUE(STUDY2_EVALUATE_U)
  ),
  published_competitors = STUDY2_PUBLISHED_COMPETITORS,
  competitor_controls = study2_competitor_controls
)
CACHE_SPEC$fingerprint <- mixedgp_object_fingerprint(CACHE_SPEC)
CACHE_TAG <- paste0(
  "s2v15-q", length(m_vec), "-", substr(CACHE_SPEC$fingerprint, 1L, 16L)
)
saveRDS(
  CACHE_SPEC,
  file.path(RES_DIR, paste0("cache_spec_", CACHE_TAG, ".rds")),
  version = 3L
)
writeLines(
  capture.output(dput(CACHE_SPEC)),
  file.path(RES_DIR, paste0("cache_spec_", CACHE_TAG, ".R"))
)

method_levels <- c("Oracle", "EIV-GP", "UC-GP", "LVGP", "EzGP")
method_cols <- c(
  "Oracle" = "black",
  "EIV-GP" = "firebrick",
  "UC-GP" = "steelblue4",
  "LVGP" = "darkorange3",
  "EzGP" = "purple4"
)
metric_levels <- c(
  "RMSE", "MAE", "CRPS", "NLPD", "Coverage95", "Width95",
  "IntervalScore95"
)

safe_se <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  sd(x) / sqrt(length(x))
}

format_mean_se <- function(mean, se, digits = 3L) {
  if (!is.finite(mean)) return("--")
  if (!is.finite(se)) return(sprintf(paste0("%.", digits, "f"), mean))
  paste0(
    sprintf(paste0("%.", digits, "f"), mean),
    " (", sprintf(paste0("%.", digits, "f"), se), ")"
  )
}

save_plot <- function(path_no_ext, plot, width, height, dpi = 320) {
  if (isTRUE(STUDY2_SAVE_PDF)) {
    ggsave(paste0(path_no_ext, ".pdf"), plot, width = width, height = height)
  }
  if (isTRUE(STUDY2_SAVE_PNG)) {
    ggsave(
      paste0(path_no_ext, ".png"), plot,
      width = width, height = height, dpi = dpi, bg = "white"
    )
  }
  invisible(plot)
}

write_csv_safe <- function(x, path) {
  if (ncol(x) == 0L) x <- data.frame(note = character(0))
  write.csv(x, path, row.names = FALSE)
  invisible(x)
}

extract_study2_diagnostics <- function(fit, rep_id, scenario, n_calib) {
  out <- fit$diagnostics$summary
  if (is.null(out)) return(data.frame())
  out <- as.data.frame(out)
  out$rep <- rep_id
  out$scenario <- scenario
  out$n_calib <- n_calib
  out
}

study2_sampler_control_settings <- function(fit,
                                             prediction_seed,
                                             mean_seed,
                                             prediction_sampler_used,
                                             mean_sampler_used) {
  diag <- fit$diagnostics$summary
  c(
    fit$control,
    list(
      sampler_strategy = fit$sampler_strategy,
      identification = fit$data$ident,
      kernel = fit$kernel$name,
      matern_nu = fit$kernel$matern_nu,
      n_iter = diag$n_iter,
      burn = diag$burn,
      thin = diag$thin,
      n_chains = diag$n_chains,
      saved_per_chain = diag$saved_per_chain,
      parallel_backend = diag$parallel_backend,
      parallel_cores = diag$parallel_cores,
      initialization_rule = paste(
        "threshold-compatible ordinal-score latent initialization with",
        "independent seeded perturbations; see",
        "fit_eivgp_ordprobit_fb()"
      ),
      chain_seeds = fit$mcmc$chain_stats$seed,
      covariance_jitter = diag$covariance_jitter,
      forms_explicit_covariance_inverse =
        diag$forms_explicit_covariance_inverse,
      prospective_latent_sampler_requested = predictive_latent_sampler,
      predictive_latent_sampler_used = prediction_sampler_used,
      mean_latent_sampler_used = mean_sampler_used,
      diagnostic_gibbs_sweeps = diagnostic_n_new_latent_gibbs,
      rejection_max_batches = rejection_max_batches,
      predictive_seed = prediction_seed,
      mean_seed = mean_seed,
      predictive_draws_requested = n_pred_draw,
      mean_evaluation_points = n_m_eval,
      mean_posterior_draws_requested = n_m_draw,
      mean_latent_integration_draws = n_m_latent,
      mean_truth_latent_draws = n_m_truth
    )
  )
}

study2_sampler_control_rows <- function(fit,
                                        rep_id,
                                        scenario,
                                        n_calib,
                                        prediction_seed,
                                        mean_seed,
                                        prediction_sampler_used,
                                        mean_sampler_used) {
  encode <- function(value) {
    if (is.numeric(value)) {
      return(paste(
        format(value, digits = 17L, scientific = FALSE),
        collapse = ";"
      ))
    }
    paste(as.character(value), collapse = ";")
  }
  settings <- study2_sampler_control_settings(
    fit = fit,
    prediction_seed = prediction_seed,
    mean_seed = mean_seed,
    prediction_sampler_used = prediction_sampler_used,
    mean_sampler_used = mean_sampler_used
  )
  data.frame(
    rep = as.integer(rep_id),
    scenario = as.character(scenario),
    n_calib = as.integer(n_calib),
    setting = names(settings),
    value = vapply(settings, encode, character(1L)),
    stringsAsFactors = FALSE
  )
}

study2_mcmc_pass <- function(diag_row,
                             rhat_limit = 1.05,
                             ess_limit = 200) {
  if (nrow(diag_row) == 0L) return(FALSE)
  rhat_cols <- intersect(
    c("max_rhat_hyper", "max_rhat_A", "max_rhat_tau", "max_rhat_missing_U"),
    names(diag_row)
  )
  rhat_ok <- length(rhat_cols) > 0L && all(
    vapply(diag_row[rhat_cols], function(z) is.finite(z) && z <= rhat_limit, logical(1))
  )
  ess_ok <- "min_ess_key" %in% names(diag_row) &&
    is.finite(diag_row$min_ess_key) && diag_row$min_ess_key >= ess_limit
  rhat_ok && ess_ok
}

study2_target_chain_diagnostics <- function(draw_matrix,
                                            draw_ids,
                                            draw_info,
                                            columns,
                                            test_rows,
                                            target,
                                            rhat_limit,
                                            bulk_ess_limit,
                                            tail_ess_limit) {
  draw_matrix <- as.matrix(draw_matrix)
  draw_ids <- as.integer(draw_ids)
  columns <- as.integer(columns)
  test_rows <- as.integer(test_rows)
  if (nrow(draw_matrix) != length(draw_ids) ||
      length(columns) != length(test_rows) ||
      anyNA(draw_ids) || any(draw_ids < 1L | draw_ids > nrow(draw_info)) ||
      anyNA(columns) || any(columns < 1L | columns > ncol(draw_matrix))) {
    stop("Incompatible target-diagnostic draws, indices, or MCMC metadata.")
  }
  if (!all(c("chain", "draw_within_chain") %in% names(draw_info))) {
    stop("mcmc_draw_info does not contain chain labels.")
  }
  selected_info <- draw_info[draw_ids, , drop = FALSE]
  chain_ids <- sort(unique(draw_info$chain))

  do.call(rbind, lapply(seq_along(columns), function(ii) {
    chains <- lapply(chain_ids, function(chain_id) {
      draw_matrix[selected_info$chain == chain_id, columns[ii]]
    })
    rhat <- rank_rhat(chains)
    ess_bulk <- bulk_ess(chains)
    ess_tail <- tail_ess(chains)
    pass <- is.finite(rhat) && rhat <= rhat_limit &&
      is.finite(ess_bulk) && ess_bulk >= bulk_ess_limit &&
      is.finite(ess_tail) && ess_tail >= tail_ess_limit
    data.frame(
      target = target,
      evaluation_position = ii,
      draw_matrix_column = columns[ii],
      test_row = test_rows[ii],
      n_chains = length(chain_ids),
      min_draws_per_chain = min(lengths(chains)),
      max_draws_per_chain = max(lengths(chains)),
      rhat = rhat,
      ess_bulk = ess_bulk,
      ess_tail = ess_tail,
      rhat_limit = rhat_limit,
      bulk_ess_limit = bulk_ess_limit,
      tail_ess_limit = tail_ess_limit,
      target_pass = pass,
      stringsAsFactors = FALSE
    )
  }))
}

study2_finite_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else max(x)
}

study2_finite_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else min(x)
}

study2_ablation_status <- function(method,
                                   n_calib,
                                   status,
                                   message = "",
                                   elapsed_seconds = NA_real_) {
  data.frame(
    method = method,
    n_calib = n_calib,
    status = status,
    message = message,
    elapsed_seconds = elapsed_seconds,
    stringsAsFactors = FALSE
  )
}

study2_optimizer_attempt_rows <- function(attempts,
                                          method,
                                          n_calib,
                                          selected_start = NA_integer_) {
  if (is.null(attempts) || !is.data.frame(attempts) || nrow(attempts) == 0L) {
    return(data.frame())
  }
  out <- attempts
  out$method <- method
  out$n_calib <- n_calib
  out$selected <- !is.na(selected_start) & out$start == selected_start
  leading <- c("method", "n_calib", "start", "selected")
  out[, c(leading, setdiff(names(out), leading)), drop = FALSE]
}

############################################################
## One scenario-replication pair
############################################################

run_one_study2_replication <- function(rep_id, scenario) {
  calib_grid <- scenario_calib_grid(scenario)
  scenario_id <- match(scenario, allowed_scenarios)
  ## The three Gaussian scenarios use common random numbers. For a fixed
  ## replication, they share X, U, Gaussian score innovations, and response
  ## innovations. The logistic robustness scenario shares X and U; its score
  ## errors necessarily use a different random-number transformation.
  data_seed_base <- 1000000L + 10000L * rep_id
  fit_seed_base <- 1000000L * scenario_id + 10000L * rep_id

  data_path <- file.path(
    STUDY2_DATA_DIR,
    mixedgp_dataset_filename(
      "study2", rep_id, scenario, n_train, n_test, m_vec[1L], q = length(m_vec)
    )
  )
  if (!file.exists(data_path)) {
    stop(
      "Missing frozen Study II dataset: ", data_path,
      ". Run 12_generate_synthetic_datasets.R first."
    )
  }
  frozen <- load_mixedgp_synthetic_dataset_strict(
    data_path,
    expected = list(
      study = "study2", scenario = scenario, rep_id = rep_id,
      n = n_train, n_test = n_test, q = length(m_vec), m = m_vec[1L],
      calib_grid = calib_grid
    )
  )
  dat <- frozen$data
  train <- dat$train
  test <- dat$test

  pattern_info <- classify_study2_pattern_frequency(train$C, test$C)
  pattern_counts <- as.data.frame(table(pattern_info$pattern_stratum))
  names(pattern_counts) <- c("evaluation_stratum", "n_test")
  pattern_counts$rep <- rep_id
  pattern_counts$scenario <- scenario

  calib_sets <- frozen$calibration_sets[as.character(calib_grid)]
  if (any(vapply(calib_sets, is.null, logical(1)))) {
    stop("Frozen Study II dataset does not contain every required calibration split.")
  }

  oracle_pool <- make_oracle_pool_2d(
    true_params = dat$true_params,
    n_pool = n_oracle_pool,
    seed = fit_seed_base + 3L
  )
  oracle_draws <- sample_oracle_test_y_2d(
    X_test = test$X,
    C_test = test$C,
    true_params = dat$true_params,
    sigma_eps = dat$sigma_eps,
    n_draw = n_pred_draw,
    oracle_pool = oracle_pool,
    seed = fit_seed_base + 4L
  )

  competitor_result <- run_study2_published_competitors(
    X_train = train$X,
    y_train = train$y,
    C_train = train$C,
    X_test = test$X,
    C_test = test$C,
    n_draw = n_pred_draw,
    seed = fit_seed_base + 100L,
    m_vec = m_vec,
    methods = STUDY2_PUBLISHED_COMPETITORS,
    ## Record failures in external competitor optimizers and continue to run
    ## EIV-GP MCMC for this frozen replication.
    strict = FALSE,
    controls = study2_competitor_controls
  )
  competitor_status <- competitor_result$status
  if (nrow(competitor_status) > 0L) {
    competitor_status$rep <- rep_id
    competitor_status$scenario <- scenario
  } else {
    competitor_status$rep <- integer(0)
    competitor_status$scenario <- character(0)
  }

  set.seed(data_seed_base + 150L)
  mean_eval_idx <- sort(sample(
    seq_len(nrow(test$X)), min(n_m_eval, nrow(test$X))
  ))
  target_gate_local_idx <- unique(round(seq(
    1L,
    length(mean_eval_idx),
    length.out = min(STUDY2_MCMC_TARGET_N_POINTS, length(mean_eval_idx))
  )))
  target_gate_test_rows <- mean_eval_idx[target_gate_local_idx]
  mean_truth <- oracle_m0_2d(
    X = test$X[mean_eval_idx, , drop = FALSE],
    C = test$C[mean_eval_idx, , drop = FALSE],
    true_params = dat$true_params,
    oracle_pool = oracle_pool,
    n_latent = n_m_truth,
    seed = data_seed_base + 151L
  )
  mean_truth_rejection <- attr(mean_truth, "rejection_telemetry")
  mean_truth_diagnostics <- attr(mean_truth, "truth_diagnostics")
  mean_truth_rejection$rep <- rep_id
  mean_truth_rejection$scenario <- scenario
  mean_truth_diagnostics$rep <- rep_id
  mean_truth_diagnostics$scenario <- scenario
  mean_recovery <- list()

  metrics <- list(
    Oracle = summarize_predictive_samples_by_pattern(
      draw_mat = oracle_draws,
      y_true = test$y,
      pattern_stratum = pattern_info$pattern_stratum,
      method = "Oracle",
      rep_id = rep_id,
      n_calib = NA_integer_,
      scenario = scenario
    )
  )

  for (method in names(competitor_result$draws)) {
    metrics[[method]] <- summarize_predictive_samples_by_pattern(
      draw_mat = competitor_result$draws[[method]],
      y_true = test$y,
      pattern_stratum = pattern_info$pattern_stratum,
      method = method,
      rep_id = rep_id,
      n_calib = NA_integer_,
      scenario = scenario
    )
    mean_recovery[[method]] <- summarize_mean_recovery_2d(
      matrix(
        competitor_result$latent_means[[method]][mean_eval_idx],
        nrow = 1L
      ),
      m_true = mean_truth,
      method = method,
      rep_id = rep_id,
      n_calib = NA_integer_,
      scenario = scenario,
      valid_function_draws = FALSE
    )
  }

  diagnostics <- list()
  target_diagnostics <- list()
  imputation <- list()
  imputation_status <- list()
  ablation_imputation <- list()
  surface <- list()
  ablation_metrics <- list()
  ablation_surface <- list()
  measurement_diagnostics <- list()
  measurement_parameter_diagnostics <- list()
  ablation_optimizer_attempts <- list()
  ablation_status <- list()
  sampler_controls <- list()
  sampler_control_manifest <- list()

  for (n_calib in calib_grid) {
    message(
      "Study II scenario=", scenario,
      ", replication=", rep_id,
      ", |O|=", n_calib
    )

    fit <- fit_eivgp_ordprobit_fb(
      X_raw = train$X,
      y_raw = train$y,
      C_ord = train$C,
      U_obs = train$U,
      calib_idx = calib_sets[[as.character(n_calib)]],
      U_true_eval = train$U,
      d = d_latent,
      m_vec = m_vec,
      ident = ident_method,
      n_iter = mc_n_iter,
      burn = mc_burn,
      thin = mc_thin,
      n_chains = mc_n_chains,
      preset = mc_preset,
      sampler_strategy = "interwoven",
      store_scores = FALSE,
      seed = fit_seed_base + 1000L + n_calib,
      parallel_chains = parallel_chains,
      verbose = FALSE
    )

    draw_ids <- seq_len(dim(fit$mcmc$samples_U)[1])
    if (length(draw_ids) > n_pred_draw) {
      draw_ids <- draw_ids[
        unique(round(seq(1, length(draw_ids), length.out = n_pred_draw)))
      ]
    }
    predictive_seed <- fit_seed_base + 2000L + n_calib
    mean_seed <- fit_seed_base + 2500L + n_calib

    eiv_draws <- sample_eiv_test_y_ordprobit_fb(
      X_test_raw = test$X,
      C_test = test$C,
      fit_obj = fit,
      draw_ids = draw_ids,
      n_per_draw = 1L,
      latent_sampler = predictive_latent_sampler,
      n_new_latent_gibbs = diagnostic_n_new_latent_gibbs,
      rejection_max_batches = rejection_max_batches,
      seed = predictive_seed
    )

    metrics[[paste0("EIV_", n_calib)]] <-
      summarize_predictive_samples_by_pattern(
        draw_mat = eiv_draws,
        y_true = test$y,
        pattern_stratum = pattern_info$pattern_stratum,
        method = "EIV-GP",
        rep_id = rep_id,
        n_calib = n_calib,
        scenario = scenario
      )

    m_draw_ids <- draw_ids
    if (length(m_draw_ids) > n_m_draw) {
      m_draw_ids <- m_draw_ids[
        unique(round(seq(1, length(m_draw_ids), length.out = n_m_draw)))
      ]
    }
    m_draws <- sample_eiv_m_given_xc_fb(
      X_test_raw = test$X[mean_eval_idx, , drop = FALSE],
      C_test = test$C[mean_eval_idx, , drop = FALSE],
      fit_obj = fit,
      draw_ids = m_draw_ids,
      n_latent = n_m_latent,
      include_process_uncertainty = TRUE,
      joint = FALSE,
      return_components = TRUE,
      latent_sampler = predictive_latent_sampler,
      n_new_latent_gibbs = diagnostic_n_new_latent_gibbs,
      rejection_max_batches = rejection_max_batches,
      seed = mean_seed
    )
    control_key <- as.character(n_calib)
    sampler_controls[[control_key]] <- study2_sampler_control_settings(
      fit = fit,
      prediction_seed = predictive_seed,
      mean_seed = mean_seed,
      prediction_sampler_used = attr(eiv_draws, "latent_sampler"),
      mean_sampler_used = attr(m_draws, "latent_sampler")
    )
    sampler_control_manifest[[control_key]] <-
      study2_sampler_control_rows(
        fit = fit,
        rep_id = rep_id,
        scenario = scenario,
        n_calib = n_calib,
        prediction_seed = predictive_seed,
        mean_seed = mean_seed,
        prediction_sampler_used = attr(eiv_draws, "latent_sampler"),
        mean_sampler_used = attr(m_draws, "latent_sampler")
      )
    mean_recovery[[paste0("EIV_", n_calib)]] <-
      summarize_mean_recovery_2d(
        m_draws,
        m_true = mean_truth,
        method = "EIV-GP",
        rep_id = rep_id,
        n_calib = n_calib,
        scenario = scenario,
        valid_function_draws = TRUE
      )

    m_conditional_means <- attr(m_draws, "conditional_means")
    if (is.null(m_conditional_means)) {
      stop("m(x,c) draws did not retain the conditional means needed by the gate.")
    }
    predictive_target_diag <- study2_target_chain_diagnostics(
      draw_matrix = eiv_draws,
      draw_ids = draw_ids,
      draw_info = fit$mcmc$mcmc_draw_info,
      columns = target_gate_test_rows,
      test_rows = target_gate_test_rows,
      target = "predictive_Y_given_X_C",
      rhat_limit = STUDY2_MCMC_RHAT_LIMIT,
      bulk_ess_limit = STUDY2_MCMC_TARGET_BULK_ESS_LIMIT,
      tail_ess_limit = STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT
    )
    mean_target_diag <- study2_target_chain_diagnostics(
      draw_matrix = m_conditional_means,
      draw_ids = m_draw_ids,
      draw_info = fit$mcmc$mcmc_draw_info,
      columns = target_gate_local_idx,
      test_rows = target_gate_test_rows,
      target = "m_conditional_mean",
      rhat_limit = STUDY2_MCMC_RHAT_LIMIT,
      bulk_ess_limit = STUDY2_MCMC_TARGET_BULK_ESS_LIMIT,
      tail_ess_limit = STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT
    )
    target_diag <- bind_rows(predictive_target_diag, mean_target_diag)
    target_diag$rep <- rep_id
    target_diag$scenario <- scenario
    target_diag$n_calib <- n_calib
    target_diag$used_for_gate <- n_calib == 0L
    target_diagnostics[[as.character(n_calib)]] <- target_diag

    diag_row <- extract_study2_diagnostics(fit, rep_id, scenario, n_calib)
    raw_coordinate_pass <- study2_mcmc_pass(
      diag_row,
      rhat_limit = STUDY2_MCMC_RHAT_LIMIT,
      ess_limit = STUDY2_MCMC_RAW_ESS_LIMIT
    )
    target_functional_pass <- nrow(target_diag) ==
      2L * length(target_gate_local_idx) && all(target_diag$target_pass)
    predictive_diag <- target_diag[
      target_diag$target == "predictive_Y_given_X_C", , drop = FALSE
    ]
    mean_diag <- target_diag[
      target_diag$target == "m_conditional_mean", , drop = FALSE
    ]
    diag_row$mcmc_gate_rule <- if (n_calib == 0L) {
      "target_functionals_without_calibration"
    } else {
      "full_raw_coordinates_with_calibration"
    }
    diag_row$mcmc_rhat_limit <- STUDY2_MCMC_RHAT_LIMIT
    diag_row$mcmc_raw_ess_limit <- STUDY2_MCMC_RAW_ESS_LIMIT
    diag_row$mcmc_target_bulk_ess_limit <-
      STUDY2_MCMC_TARGET_BULK_ESS_LIMIT
    diag_row$mcmc_target_tail_ess_limit <-
      STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT
    diag_row$mcmc_target_n_points <- length(target_gate_local_idx)
    diag_row$raw_coordinate_pass <- raw_coordinate_pass
    diag_row$target_functional_pass <- target_functional_pass
    diag_row$predictive_max_rhat <- study2_finite_max(predictive_diag$rhat)
    diag_row$predictive_min_bulk_ess <-
      study2_finite_min(predictive_diag$ess_bulk)
    diag_row$predictive_min_tail_ess <-
      study2_finite_min(predictive_diag$ess_tail)
    diag_row$m_max_rhat <- study2_finite_max(mean_diag$rhat)
    diag_row$m_min_bulk_ess <- study2_finite_min(mean_diag$ess_bulk)
    diag_row$m_min_tail_ess <- study2_finite_min(mean_diag$ess_tail)
    diag_row$mcmc_pass <- if (n_calib == 0L) {
      target_functional_pass
    } else {
      raw_coordinate_pass
    }
    diagnostics[[as.character(n_calib)]] <- diag_row

    if (isTRUE(STUDY2_EVALUATE_U)) {
      imputation[[paste0("EIV_training_", n_calib)]] <-
        study2_latent_imputation_metrics(
          fit = fit,
          U_true = train$U,
          rep_id = rep_id,
          n_calib = n_calib,
          scenario = scenario
        )
      imputation_status[[paste0("EIV_training_", n_calib)]] <-
        study2_latent_imputation_status(
          fit = fit,
          method = "EIV-GP",
          target = "training_missing_U",
          rep_id = rep_id,
          n_calib = n_calib,
          scenario = scenario,
          n_units = length(fit$data$miss_idx)
        )
      imputation[[paste0("EIV_prospective_", n_calib)]] <-
        study2_eiv_prospective_imputation_metrics(
          fit = fit,
          C_new = test$C[mean_eval_idx, , drop = FALSE],
          U_true = test$U[mean_eval_idx, , drop = FALSE],
          rep_id = rep_id,
          n_calib = n_calib,
          scenario = scenario,
          max_draw = n_m_draw,
          latent_sampler = predictive_latent_sampler,
          n_new_latent_gibbs = diagnostic_n_new_latent_gibbs,
          rejection_max_batches = rejection_max_batches,
          seed = fit_seed_base + 2700L + n_calib
        )
      imputation_status[[paste0("EIV_prospective_", n_calib)]] <-
        study2_latent_imputation_status(
          fit = fit,
          method = "EIV-GP",
          target = "prospective_U_given_C",
          rep_id = rep_id,
          n_calib = n_calib,
          scenario = scenario,
          n_units = length(mean_eval_idx)
        )
    }

    if (isTRUE(STUDY2_EVALUATE_F)) {
      surface[[as.character(n_calib)]] <-
        study2_eiv_surface_recovery_metrics(
          fit = fit,
          scenario = scenario,
          rep_id = rep_id,
          n_calib = n_calib,
          grid_n = if (STUDY2_CONFIG == "quick") 15L else 31L,
          max_draw = if (STUDY2_CONFIG == "quick") 40L else 200L
        )
    }

    if (isTRUE(STUDY2_SAVE_REP_FITS)) {
      saveRDS(
        fit,
        file.path(
          FIT_DIR,
          sprintf("fit_%s_rep%03d_calib%03d.rds", scenario, rep_id, n_calib)
        )
      )
    }

    rm(fit, eiv_draws, m_draws)
    invisible(gc())
  }

  ##########################################################
  ## Appendix-only ablations (kept out of the main table)
  ##########################################################

  if (isTRUE(STUDY2_RUN_ABLATIONS) &&
      scenario %in% STUDY2_ABLATION_SCENARIOS) {
    if (isTRUE(STUDY2_EVALUATE_F)) {
      start_full_u <- proc.time()[3]
      full_u_gp <- tryCatch(
        fit_study2_latent_gp(
          train$X, train$y, train$U,
          gp_n_starts = STUDY2_ABLATION_GP_N_STARTS,
          gp_seed = fit_seed_base + 8000L,
          gp_maxit = STUDY2_ABLATION_GP_MAXIT
        ),
        error = function(e) e
      )
      elapsed_full_u <- proc.time()[3] - start_full_u

      if (inherits(full_u_gp, "error")) {
        ablation_optimizer_attempts[["Full-U GP"]] <-
          study2_optimizer_attempt_rows(
            full_u_gp$optimizer_attempts,
            "Full-U GP",
            NA_integer_
          )
        ablation_status[["Full-U GP"]] <- study2_ablation_status(
          "Full-U GP", NA_integer_, "failed",
          conditionMessage(full_u_gp), elapsed_full_u
        )
      } else {
        ablation_optimizer_attempts[["Full-U GP"]] <-
          study2_optimizer_attempt_rows(
            full_u_gp$optimizer_attempts,
            "Full-U GP",
            NA_integer_,
            full_u_gp$selected_start
          )
        ablation_status[["Full-U GP"]] <- study2_ablation_status(
          "Full-U GP", NA_integer_, "success", "", elapsed_full_u
        )
        ablation_surface[["Full-U GP"]] <-
          study2_latent_gp_surface_metrics(
            gp_fit = full_u_gp,
            scenario = scenario,
            method = "Full-U GP",
            rep_id = rep_id,
            n_calib = NA_integer_,
            grid_n = if (STUDY2_CONFIG == "quick") 15L else 31L
          )
      }
    }

    for (n_calib in calib_grid) {
      key <- as.character(n_calib)
      calib_idx <- calib_sets[[key]]

      start_measurement <- proc.time()[3]
      measurement_fit <- tryCatch(
        fit_ordinalprobit_measurement_fb(
          C_ord = train$C,
          U_obs = train$U,
          calib_idx = calib_idx,
          d = d_latent,
          m_vec = m_vec,
          ident = ident_method,
          n_iter = measurement_n_iter,
          burn = measurement_burn,
          thin = measurement_thin,
          n_chains = measurement_n_chains,
          seed = fit_seed_base + 5000L + n_calib,
          parallel_chains = parallel_chains,
          rhat_limit = STUDY2_MEAS_RHAT_LIMIT,
          bulk_ess_limit = STUDY2_MEAS_BULK_ESS_LIMIT,
          tail_ess_limit = STUDY2_MEAS_TAIL_ESS_LIMIT,
          verbose = FALSE
        ),
        error = function(e) e
      )
      elapsed_measurement <- proc.time()[3] - start_measurement

      if (inherits(measurement_fit, "error")) {
        msg <- conditionMessage(measurement_fit)
        measurement_status_stub <- list(data = list(
          U_obs = train$U,
          calib_idx = calib_idx,
          miss_idx = setdiff(seq_len(nrow(train$U)), calib_idx),
          d = d_latent
        ))
        for (task_spec in list(
          list(
            target = "training_missing_U",
            n_units = length(measurement_status_stub$data$miss_idx)
          ),
          list(
            target = "prospective_U_given_C",
            n_units = length(mean_eval_idx)
          )
        )) {
          task_status <- study2_latent_imputation_status(
            fit = measurement_status_stub,
            method = "Ordinal model (no Y)",
            target = task_spec$target,
            rep_id = rep_id,
            n_calib = n_calib,
            scenario = scenario,
            n_units = task_spec$n_units
          )
          if (isTRUE(task_status$task_eligible)) {
            task_status$status <- "failed"
            task_status$task_eligible <- FALSE
            task_status$reason <- msg
          }
          imputation_status[[paste0(
            "measurement_", task_spec$target, "_", key
          )]] <- task_status
        }
        ablation_status[[paste0("measurement_", key)]] <-
          study2_ablation_status(
            "Response-free measurement model", n_calib, "failed",
            msg, elapsed_measurement
          )
        ablation_status[[paste0("PI_", key)]] <-
          study2_ablation_status(
            "PI-GP", n_calib, "failed", msg, NA_real_
          )
        ablation_status[[paste0("CC_", key)]] <-
          study2_ablation_status(
            "CC-GP", n_calib,
            if (n_calib == 0L) "not_applicable" else "failed",
            if (n_calib == 0L) "No complete cases." else msg,
            NA_real_
          )
        next
      }

      measurement_diag <- measurement_fit$diagnostics
      measurement_diag$rep <- rep_id
      measurement_diag$scenario <- scenario
      measurement_diag$n_calib <- n_calib
      measurement_diagnostics[[key]] <- measurement_diag

      key_parameter_diag <- measurement_fit$diagnostic_parameters$key
      latent_parameter_diag <- measurement_fit$diagnostic_parameters$latent_rhat
      if (nrow(latent_parameter_diag) > 0L) {
        latent_parameter_diag$block <- "latent"
      }
      parameter_diag <- bind_rows(
        key_parameter_diag[, c(
          "parameter", "block", "rhat", "ess_bulk", "ess_tail"
        )],
        latent_parameter_diag[, intersect(
          c("parameter", "block", "rhat", "ess_bulk", "ess_tail"),
          names(latent_parameter_diag)
        )]
      )
      parameter_diag$rep <- rep_id
      parameter_diag$scenario <- scenario
      parameter_diag$n_calib <- n_calib
      measurement_parameter_diagnostics[[key]] <- parameter_diag

      measurement_training_status <- study2_latent_imputation_status(
        fit = measurement_fit,
        method = "Ordinal model (no Y)",
        target = "training_missing_U",
        rep_id = rep_id,
        n_calib = n_calib,
        scenario = scenario,
        n_units = length(measurement_fit$data$miss_idx)
      )
      measurement_prospective_status <- study2_latent_imputation_status(
        fit = measurement_fit,
        method = "Ordinal model (no Y)",
        target = "prospective_U_given_C",
        rep_id = rep_id,
        n_calib = n_calib,
        scenario = scenario,
        n_units = length(mean_eval_idx)
      )

      if (!isTRUE(measurement_diag$convergence_pass)) {
        msg <- paste0(
          "Response-free measurement-model convergence gate failed: ",
          "backend=", measurement_diag$diagnostic_backend,
          ", max R-hat=", signif(study2_finite_max(c(
            measurement_diag$max_rhat_A,
            measurement_diag$max_rhat_tau,
            measurement_diag$max_rhat_missing_U
          )), 4),
          ", min bulk/tail ESS=",
          signif(measurement_diag$min_bulk_ess_all, 5), "/",
          signif(measurement_diag$min_tail_ess_all, 5), "."
        )
        ablation_status[[paste0("measurement_", key)]] <-
          study2_ablation_status(
            "Response-free measurement model", n_calib,
            "failed_convergence", msg, elapsed_measurement
          )
        ablation_status[[paste0("PI_", key)]] <-
          study2_ablation_status(
            "PI-GP", n_calib, "failed_convergence", msg, NA_real_
          )
        ablation_status[[paste0("CC_", key)]] <-
          study2_ablation_status(
            "CC-GP", n_calib,
            if (n_calib == 0L) "not_applicable" else "failed_convergence",
            if (n_calib == 0L) "No complete cases." else msg,
            NA_real_
          )
        for (task_status in list(
          measurement_training_status,
          measurement_prospective_status
        )) {
          if (isTRUE(task_status$task_eligible)) {
            task_status$status <- "failed_convergence"
            task_status$task_eligible <- FALSE
            task_status$reason <- msg
          }
          imputation_status[[paste0(
            "measurement_", task_status$target, "_", key
          )]] <- task_status
        }
        next
      }
      imputation_status[[paste0("measurement_training_", key)]] <-
        measurement_training_status
      imputation_status[[paste0("measurement_prospective_", key)]] <-
        measurement_prospective_status
      ablation_status[[paste0("measurement_", key)]] <-
        study2_ablation_status(
          "Response-free measurement model", n_calib, "success", "",
          elapsed_measurement
        )

      ablation_imputation[[paste0("measurement_training_", key)]] <-
        study2_measurement_training_imputation_metrics(
          measurement_fit = measurement_fit,
          U_true = train$U,
          rep_id = rep_id,
          n_calib = n_calib,
          scenario = scenario
        )
      ablation_imputation[[paste0("measurement_prospective_", key)]] <-
        study2_measurement_prospective_imputation_metrics(
          measurement_fit = measurement_fit,
          C_new = test$C[mean_eval_idx, , drop = FALSE],
          U_true = test$U[mean_eval_idx, , drop = FALSE],
          rep_id = rep_id,
          n_calib = n_calib,
          scenario = scenario,
          max_draw = n_m_draw,
          latent_sampler = predictive_latent_sampler,
          n_gibbs = diagnostic_n_new_latent_gibbs,
          rejection_max_batches = rejection_max_batches,
          seed = fit_seed_base + 5600L + n_calib
        )

      start_pi <- proc.time()[3]
      pi_fit <- tryCatch(
        fit_study2_pi_gp(
          train$X,
          train$y,
          measurement_fit,
          gp_n_starts = STUDY2_ABLATION_GP_N_STARTS,
          gp_seed = fit_seed_base + 6100L + n_calib,
          gp_maxit = STUDY2_ABLATION_GP_MAXIT
        ),
        error = function(e) e
      )
      pi_result <- if (inherits(pi_fit, "error")) {
        pi_fit
      } else {
        ablation_optimizer_attempts[[paste0("PI_", key)]] <-
          study2_optimizer_attempt_rows(
            pi_fit$gp$optimizer_attempts,
            "PI-GP",
            n_calib,
            pi_fit$gp$selected_start
          )
        tryCatch(
          list(
            fit = pi_fit,
            draws = sample_study2_pi_gp(
              fit = pi_fit,
              X_test = test$X,
              C_test = test$C,
              n_draw = n_pred_draw,
              n_measurement_draw = min(200L, n_pred_draw),
              latent_sampler = predictive_latent_sampler,
              n_gibbs = diagnostic_n_new_latent_gibbs,
              rejection_max_batches = rejection_max_batches,
              seed = fit_seed_base + 6000L + n_calib
            )
          ),
          error = function(e) e
        )
      }
      elapsed_pi <- proc.time()[3] - start_pi

      if (inherits(pi_result, "error")) {
        if (is.null(ablation_optimizer_attempts[[paste0("PI_", key)]])) {
          ablation_optimizer_attempts[[paste0("PI_", key)]] <-
            study2_optimizer_attempt_rows(
              pi_result$optimizer_attempts,
              "PI-GP",
              n_calib
            )
        }
        ablation_status[[paste0("PI_", key)]] <-
          study2_ablation_status(
            "PI-GP", n_calib, "failed", conditionMessage(pi_result),
            elapsed_pi
          )
      } else {
        ablation_status[[paste0("PI_", key)]] <-
          study2_ablation_status(
            "PI-GP", n_calib, "success", "", elapsed_pi
          )
        ablation_metrics[[paste0("PI_", key)]] <-
          summarize_predictive_samples_by_pattern(
            draw_mat = pi_result$draws,
            y_true = test$y,
            pattern_stratum = pattern_info$pattern_stratum,
            method = "PI-GP",
            rep_id = rep_id,
            n_calib = n_calib,
            scenario = scenario
          )
        if (isTRUE(STUDY2_EVALUATE_F) && n_calib > 0L) {
          ablation_surface[[paste0("PI_", key)]] <-
            study2_latent_gp_surface_metrics(
              gp_fit = pi_result$fit$gp,
              scenario = scenario,
              method = "PI-GP",
              rep_id = rep_id,
              n_calib = n_calib,
              grid_n = if (STUDY2_CONFIG == "quick") 15L else 31L
            )
        }
      }

      if (n_calib == 0L) {
        ablation_status[[paste0("CC_", key)]] <-
          study2_ablation_status(
            "CC-GP", n_calib, "not_applicable", "No complete cases.",
            0
          )
      } else {
        start_cc <- proc.time()[3]
        cc_fit <- tryCatch(
          fit_study2_cc_gp(
            X = train$X,
            y = train$y,
            U = train$U,
            calib_idx = calib_idx,
            measurement_fit = measurement_fit,
            gp_n_starts = STUDY2_ABLATION_GP_N_STARTS,
            gp_seed = fit_seed_base + 7100L + n_calib,
            gp_maxit = STUDY2_ABLATION_GP_MAXIT
          ),
          error = function(e) e
        )
        cc_result <- if (inherits(cc_fit, "error")) {
          cc_fit
        } else {
          ablation_optimizer_attempts[[paste0("CC_", key)]] <-
            study2_optimizer_attempt_rows(
              cc_fit$gp$optimizer_attempts,
              "CC-GP",
              n_calib,
              cc_fit$gp$selected_start
            )
          tryCatch(
            list(
              fit = cc_fit,
              draws = sample_study2_cc_gp(
                fit = cc_fit,
                X_test = test$X,
                C_test = test$C,
                n_draw = n_pred_draw,
                latent_sampler = predictive_latent_sampler,
                n_gibbs = diagnostic_n_new_latent_gibbs,
                rejection_max_batches = rejection_max_batches,
                seed = fit_seed_base + 7000L + n_calib
              )
            ),
            error = function(e) e
          )
        }
        elapsed_cc <- proc.time()[3] - start_cc

        if (inherits(cc_result, "error")) {
          if (is.null(ablation_optimizer_attempts[[paste0("CC_", key)]])) {
            ablation_optimizer_attempts[[paste0("CC_", key)]] <-
              study2_optimizer_attempt_rows(
                cc_result$optimizer_attempts,
                "CC-GP",
                n_calib
              )
          }
          ablation_status[[paste0("CC_", key)]] <-
            study2_ablation_status(
              "CC-GP", n_calib, "failed", conditionMessage(cc_result),
              elapsed_cc
            )
        } else {
          ablation_status[[paste0("CC_", key)]] <-
            study2_ablation_status(
              "CC-GP", n_calib, "success", "", elapsed_cc
            )
          ablation_metrics[[paste0("CC_", key)]] <-
            summarize_predictive_samples_by_pattern(
              draw_mat = cc_result$draws,
              y_true = test$y,
              pattern_stratum = pattern_info$pattern_stratum,
              method = "CC-GP",
              rep_id = rep_id,
              n_calib = n_calib,
              scenario = scenario
            )
          if (isTRUE(STUDY2_EVALUATE_F)) {
            ablation_surface[[paste0("CC_", key)]] <-
              study2_latent_gp_surface_metrics(
                gp_fit = cc_result$fit$gp,
                scenario = scenario,
                method = "CC-GP",
                rep_id = rep_id,
                n_calib = n_calib,
                grid_n = if (STUDY2_CONFIG == "quick") 15L else 31L
              )
          }
        }
      }

      rm(measurement_fit, pi_result, pi_fit)
      if (exists("cc_result")) rm(cc_result)
      if (exists("cc_fit")) rm(cc_fit)
      invisible(gc())
    }
  }

  ablation_status_df <- bind_rows(ablation_status)
  if (nrow(ablation_status_df) > 0L) {
    ablation_status_df$rep <- rep_id
    ablation_status_df$scenario <- scenario
  } else {
    ablation_status_df$rep <- integer(0)
    ablation_status_df$scenario <- character(0)
  }
  ablation_optimizer_attempts_df <- bind_rows(ablation_optimizer_attempts)
  if (nrow(ablation_optimizer_attempts_df) > 0L) {
    ablation_optimizer_attempts_df$rep <- rep_id
    ablation_optimizer_attempts_df$scenario <- scenario
  } else {
    ablation_optimizer_attempts_df$rep <- integer(0)
    ablation_optimizer_attempts_df$scenario <- character(0)
  }

  list(
    metrics = bind_rows(metrics),
    mean_recovery = bind_rows(mean_recovery),
    mean_truth_rejection = mean_truth_rejection,
    mean_truth_diagnostics = mean_truth_diagnostics,
    diagnostics = bind_rows(diagnostics),
    target_diagnostics = bind_rows(target_diagnostics),
    imputation = bind_rows(c(imputation, ablation_imputation)),
    imputation_status = bind_rows(imputation_status),
    surface = bind_rows(surface),
    ablation_metrics = bind_rows(ablation_metrics),
    ablation_surface = bind_rows(ablation_surface),
    measurement_diagnostics = bind_rows(measurement_diagnostics),
    measurement_parameter_diagnostics = bind_rows(
      measurement_parameter_diagnostics
    ),
    ablation_optimizer_attempts = ablation_optimizer_attempts_df,
    ablation_status = ablation_status_df,
    competitor_status = competitor_status,
    pattern_counts = pattern_counts,
    sampler_control_manifest = bind_rows(sampler_control_manifest),
    metadata = list(
      rep = rep_id,
      scenario = scenario,
      n_train = n_train,
      n_test = n_test,
      data_seed_base = data_seed_base,
      fit_seed_base = fit_seed_base,
      common_random_number_group = rep_id,
      calib_grid = calib_grid,
      true_params = dat$true_params,
      frozen_data_file = basename(data_path),
      frozen_data_md5 = attr(frozen, "manifest_md5"),
      frozen_manifest = attr(frozen, "manifest_path"),
      cache_tag = CACHE_TAG,
      sampler = sampler_controls,
      downstream_monte_carlo = list(
        predictive_draws_requested = n_pred_draw,
        mean_evaluation_points = n_m_eval,
        mean_posterior_draws_requested = n_m_draw,
        mean_latent_integration_draws = n_m_latent,
        mean_truth_latent_draws = n_m_truth,
        prospective_latent_sampler = predictive_latent_sampler,
        diagnostic_gibbs_sweeps = diagnostic_n_new_latent_gibbs,
        rejection_max_batches = rejection_max_batches
      )
    )
  )
}

############################################################
## Resumable execution
############################################################

run_grid <- expand.grid(
  scenario = STUDY2_SCENARIOS,
  rep = seq_len(n_rep),
  stringsAsFactors = FALSE
)

rep_files <- file.path(
  REP_DIR,
  sprintf(
    "study2_%s_rep%03d_%s.rds",
    run_grid$scenario,
    run_grid$rep,
    CACHE_TAG
  )
)

run_or_load_study2_replication <- function(ii) {
  scenario <- run_grid$scenario[ii]
  rep_id <- run_grid$rep[ii]

  if (
    isTRUE(STUDY2_USE_CACHE) &&
      isTRUE(STUDY2_MC_RESUME) &&
      file.exists(rep_files[ii])
  ) {
    cached <- readRDS(rep_files[ii])
    data_path <- file.path(
      STUDY2_DATA_DIR,
      mixedgp_dataset_filename(
        "study2", rep_id, scenario, n_train, n_test, m_vec[1L],
        q = length(m_vec)
      )
    )
    cache_valid <- is.list(cached$metadata) &&
      identical(cached$metadata$rep, as.integer(rep_id)) &&
      identical(cached$metadata$scenario, scenario) &&
      identical(cached$metadata$cache_tag, CACHE_TAG) &&
      identical(cached$metadata$frozen_data_file, basename(data_path)) &&
      file.exists(data_path) &&
      identical(
        tolower(as.character(cached$metadata$frozen_data_md5)),
        tolower(unname(tools::md5sum(data_path)))
      )
    if (isTRUE(cache_valid)) return(cached)
    message("Ignoring incompatible Study II cache: ", rep_files[ii])
  }
  out <- run_one_study2_replication(rep_id, scenario)
  saveRDS(out, rep_files[ii])
  out
}
run_study2_replication_safely <- function(ii) {
  tryCatch(
    run_or_load_study2_replication(ii),
    error = function(e) structure(
      list(
        row = as.integer(ii),
        status = "failed",
        message = conditionMessage(e)
      ),
      class = c("mixedgp_replication_failure", "list")
    )
  )
}
replication_cores <- if (identical(STUDY2_PARALLEL_LEVEL, "replications")) {
  min(nrow(run_grid), mixedgp_resolve_cores())
} else {
  1L
}
rep_objects <- mixedgp_parallel_lapply(
  as.list(seq_len(nrow(run_grid))),
  run_study2_replication_safely,
  n_cores = replication_cores,
  seeds = 920000L + seq_len(nrow(run_grid)),
  mc.preschedule = FALSE
)

failed_replications <- vapply(
  rep_objects, inherits, logical(1L), what = "mixedgp_replication_failure"
)
replication_status <- data.frame(
  run_grid,
  status = ifelse(failed_replications, "failed", "success"),
  message = vapply(seq_len(nrow(run_grid)), function(ii) {
    if (failed_replications[ii]) rep_objects[[ii]]$message else ""
  }, character(1L)),
  stringsAsFactors = FALSE
)
write.csv(
  replication_status,
  file.path(TAB_DIR, paste0("study2_replication_status_", CACHE_TAG, ".csv")),
  row.names = FALSE
)
successful_rep_objects <- rep_objects[!failed_replications]
if (length(successful_rep_objects) == 0L) {
  stop("All Study II replications failed; see study2_replication_status CSV.")
}

mc_results <- bind_rows(lapply(successful_rep_objects, `[[`, "metrics"))
mc_mean_recovery <- bind_rows(lapply(successful_rep_objects, `[[`, "mean_recovery"))
mc_mean_truth_rejection <- bind_rows(
  lapply(successful_rep_objects, `[[`, "mean_truth_rejection")
)
mc_mean_truth_diagnostics <- bind_rows(
  lapply(successful_rep_objects, `[[`, "mean_truth_diagnostics")
)
mc_diagnostics <- bind_rows(lapply(successful_rep_objects, `[[`, "diagnostics"))
mc_target_diagnostics <- bind_rows(
  lapply(successful_rep_objects, `[[`, "target_diagnostics")
)
mc_imputation <- bind_rows(lapply(successful_rep_objects, `[[`, "imputation"))
mc_imputation_status <- bind_rows(
  lapply(successful_rep_objects, `[[`, "imputation_status")
)
mc_surface <- bind_rows(lapply(successful_rep_objects, `[[`, "surface"))
mc_ablation_results <- bind_rows(lapply(successful_rep_objects, `[[`, "ablation_metrics"))
mc_ablation_surface <- bind_rows(lapply(successful_rep_objects, `[[`, "ablation_surface"))
measurement_diagnostics <- bind_rows(
  lapply(successful_rep_objects, `[[`, "measurement_diagnostics")
)
measurement_parameter_diagnostics <- bind_rows(
  lapply(successful_rep_objects, `[[`, "measurement_parameter_diagnostics")
)
ablation_optimizer_attempts <- bind_rows(
  lapply(successful_rep_objects, `[[`, "ablation_optimizer_attempts")
)
ablation_status <- bind_rows(lapply(successful_rep_objects, `[[`, "ablation_status"))
competitor_status <- bind_rows(lapply(successful_rep_objects, `[[`, "competitor_status"))
pattern_counts <- bind_rows(lapply(successful_rep_objects, `[[`, "pattern_counts"))
sampler_control_manifest <- bind_rows(
  lapply(successful_rep_objects, `[[`, "sampler_control_manifest")
)
replication_metadata <- lapply(successful_rep_objects, `[[`, "metadata")

raw_outputs <- list(
  predictive_metrics = mc_results,
  mean_recovery = mc_mean_recovery,
  mean_truth_rejection = mc_mean_truth_rejection,
  mean_truth_diagnostics = mc_mean_truth_diagnostics,
  mcmc_diagnostics = mc_diagnostics,
  mcmc_target_diagnostics = mc_target_diagnostics,
  latent_imputation = mc_imputation,
  latent_imputation_status = mc_imputation_status,
  surface_recovery = mc_surface,
  ablation_predictive_metrics = mc_ablation_results,
  ablation_surface_recovery = mc_ablation_surface,
  measurement_diagnostics = measurement_diagnostics,
  measurement_parameter_diagnostics = measurement_parameter_diagnostics,
  ablation_optimizer_attempts = ablation_optimizer_attempts,
  ablation_status = ablation_status,
  competitor_status = competitor_status,
  pattern_counts = pattern_counts,
  sampler_control_manifest = sampler_control_manifest,
  replication_metadata = replication_metadata,
  run_grid = run_grid,
  cache_tag = CACHE_TAG,
  cache_spec = CACHE_SPEC
)
saveRDS(raw_outputs, file.path(RES_DIR, paste0("study2_results_", CACHE_TAG, ".rds")))

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
    sampler_strategy = "interwoven",
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
      "n_calib=0: predictive_Y_and_m_functionals; n_calib>0: ",
      "hyperparameters_A_tau_and_missing_U"
    ),
    mcmc_rhat_limit = STUDY2_MCMC_RHAT_LIMIT,
    mcmc_raw_ess_limit = STUDY2_MCMC_RAW_ESS_LIMIT,
    mcmc_target_bulk_ess_limit = STUDY2_MCMC_TARGET_BULK_ESS_LIMIT,
    mcmc_target_tail_ess_limit = STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT,
    mcmc_target_n_points = STUDY2_MCMC_TARGET_N_POINTS,
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
    scenario %in% c("latent_additive_control", "high_uncertainty"),
    metric == "CRPS",
    (method == "EIV-GP" & n_calib == STUDY2_CONTRAST_CALIB) |
      method != "EIV-GP"
  ) |>
  mutate(
    scenario = factor(
      scenario,
      levels = c("latent_additive_control", "high_uncertainty"),
      labels = vapply(
        c("latent_additive_control", "high_uncertainty"),
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

  if (isTRUE(STUDY2_ENFORCE_MCMC_GATE) && !all(mc_diagnostics$mcmc_pass)) {
    stop(
      "The publication MCMC gate failed. Results were saved, but no run should ",
      "be reported until every failed fit is diagnosed and rerun."
    )
  }
}

cat("\nStudy II run completed under design tag:\n", STUDY2_DESIGN_TAG, "\n")
cat("Raw result bundle:\n", file.path(RES_DIR, paste0("study2_results_", CACHE_TAG, ".rds")), "\n")
cat("Primary table:\n", file.path(TAB_DIR, paste0("study2_main_summary_", CACHE_TAG, ".tex")), "\n")
