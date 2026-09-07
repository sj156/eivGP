############################################################
## 02_study2_monte_carlo.R
##
## Publication Monte Carlo driver for Study II.
##
## Canonical design:
##   * n = 100 training and 100 test observations in the current masters;
##   * two or four ordinal proxies of a two-dimensional latent state;
##   * primary and logistic misspecification scenarios;
##   * primary recovery of m(x,c)=E{f(x,U)|C=c};
##   * prediction of Y* | X*, C* with no test U*;
##   * separate latent-imputation and physically anchored f(x, u) summaries;
##   * published competitors only, with explicit availability records.
############################################################

if (!exists("mixedgp_v030_fit")) {
  source("load_mixedgp.R")
}
if (!exists("load_mixedgp_synthetic_dataset")) source("00_synthetic_data.R")
if (!exists("mixedgp_study2_target_series") || !exists("mixedgp_simulation_diagnostic_advice")) {
  source_paths <- unlist(lapply(sys.frames(), function(frame) {
    get0("ofile", envir = frame, inherits = FALSE, ifnotfound = character(0))
  }))
  command_paths <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  helper_candidates <- unique(c(
    file.path(dirname(c(rev(source_paths), command_paths)), "simulation_helpers.R"),
    file.path(getwd(), "simulation_helpers.R")
  ))
  helper_hits <- helper_candidates[file.exists(helper_candidates)]
  if (length(helper_hits) == 0L) stop("Cannot locate simulation_helpers.R beside the experiment driver.")
  source(helper_hits[1L], local = TRUE)
}
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

needed_pkgs <- c("ggplot2", "dplyr", "tidyr", "knitr", "posterior", "TruncatedNormal")
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
if (!exists("mixedgp_cached_competitors")) source("load_mixedgp.R")
study2_competitor_controls <- mixedgp_competitor_protocol("study2")

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
if (!exists("STUDY2_MCMC_RHAT_LIMIT")) STUDY2_MCMC_RHAT_LIMIT <- 1.01
if (!exists("STUDY2_MCMC_RAW_ESS_LIMIT")) STUDY2_MCMC_RAW_ESS_LIMIT <- 400
if (!exists("STUDY2_MCMC_TARGET_BULK_ESS_LIMIT")) {
  STUDY2_MCMC_TARGET_BULK_ESS_LIMIT <- 400
}
if (!exists("STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT")) {
  STUDY2_MCMC_TARGET_TAIL_ESS_LIMIT <- 400
}
if (!exists("STUDY2_MCMC_TARGET_N_POINTS")) {
  STUDY2_MCMC_TARGET_N_POINTS <- 5L
}
if (!exists("STUDY2_MCMC_TARGET_DRAWS_PER_CHAIN")) {
  STUDY2_MCMC_TARGET_DRAWS_PER_CHAIN <- Inf
}
if (!exists("STUDY2_MCMC_TARGET_N_LATENT")) {
  STUDY2_MCMC_TARGET_N_LATENT <- 64L
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
    length(STUDY2_MCMC_TARGET_DRAWS_PER_CHAIN) != 1L ||
    is.na(STUDY2_MCMC_TARGET_DRAWS_PER_CHAIN) ||
    STUDY2_MCMC_TARGET_DRAWS_PER_CHAIN < 4L ||
    length(STUDY2_MCMC_TARGET_N_LATENT) != 1L ||
    !is.finite(STUDY2_MCMC_TARGET_N_LATENT) ||
    STUDY2_MCMC_TARGET_N_LATENT < 2L ||
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
if (!STUDY2_PARALLEL_LEVEL %in% c("chains", "replications", "none", "hybrid")) {
  stop("STUDY2_PARALLEL_LEVEL must be chains, replications, or none.")
}
parallel_chains <- STUDY2_PARALLEL_LEVEL %in% c("chains", "hybrid")
chain_cores <- if (STUDY2_PARALLEL_LEVEL == "hybrid") {
  mixedgp_as_integer_strict(STUDY2_CHAIN_WORKERS, "chain workers", 1L, 1L)
} else NULL

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
  schema = "s2v17_nonfatal_diagnostic_targets",
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
    target_n_points = STUDY2_MCMC_TARGET_N_POINTS,
    target_draws_per_chain = STUDY2_MCMC_TARGET_DRAWS_PER_CHAIN,
    target_integration_draws = STUDY2_MCMC_TARGET_N_LATENT
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
  "s2v16-q", length(m_vec), "-", substr(CACHE_SPEC$fingerprint, 1L, 16L)
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

study2_finite_max <- function(x) {
  if (length(x) == 0L || any(!is.finite(x))) NA_real_ else max(x)
}

study2_finite_min <- function(x) {
  if (length(x) == 0L || any(!is.finite(x))) NA_real_ else min(x)
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
