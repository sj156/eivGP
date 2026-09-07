############################################################
## 02_study1_monte_carlo.R
##
## Publication comparison for Study I.
##
## The publication path refits every method, including EIV-GP, with the current
## audited sampler. The historical July 27 rows remain available only through
## the explicit archival switch STUDY1_REUSE_LOCKED_EIV=TRUE.
############################################################

if (!exists("mixedgp_v030_fit")) source("load_mixedgp.R")
if (!exists("load_mixedgp_synthetic_dataset")) source("00_synthetic_data.R")
if (!exists("mixedgp_study1_raw_diagnostics") || !exists("mixedgp_simulation_diagnostic_advice")) {
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
if (!exists("run_study1_published_competitors")) {
  source("03_study2_published_competitors.R")
}
if (!exists("STUDY1_RUN_ABLATIONS")) STUDY1_RUN_ABLATIONS <- TRUE
if (!exists("STUDY1_EVALUATE_F")) STUDY1_EVALUATE_F <- TRUE
if (!exists("STUDY1_EVALUATE_U")) STUDY1_EVALUATE_U <- TRUE
if (isTRUE(STUDY1_RUN_ABLATIONS) &&
    !exists("fit_threshold_measurement_response_free")) {
  source("04_study1_ablations.R")
}

needed_pkgs <- c("ggplot2", "dplyr", "tidyr", "knitr")
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

if (!exists("STUDY1_CONFIG")) {
  STUDY1_CONFIG <- if (exists("STUDY1_QUICK") && isTRUE(STUDY1_QUICK)) {
    "quick"
  } else {
    "thorough"
  }
}
if (!STUDY1_CONFIG %in% c("quick", "balanced", "thorough")) {
  stop("STUDY1_CONFIG must be 'quick', 'balanced', or 'thorough'.")
}
STUDY1_QUICK <- identical(STUDY1_CONFIG, "quick")
if (!exists("STUDY1_USE_CACHE")) STUDY1_USE_CACHE <- TRUE
if (!exists("STUDY1_OUT_PREFIX")) STUDY1_OUT_PREFIX <- ".."
if (!exists("STUDY1_REUSE_LOCKED_EIV")) {
  STUDY1_REUSE_LOCKED_EIV <- FALSE
}
if (!exists("STUDY1_REQUIRE_MCMC_GATE")) {
  STUDY1_REQUIRE_MCMC_GATE <- !identical(STUDY1_CONFIG, "quick")
}
if (!exists("STUDY1_MAX_RHAT")) STUDY1_MAX_RHAT <- 1.01
if (!exists("STUDY1_MIN_ESS")) STUDY1_MIN_ESS <- 400
if (!exists("STUDY1_DIAGNOSTIC_MAX_DRAWS")) STUDY1_DIAGNOSTIC_MAX_DRAWS <- Inf
if (!exists("STUDY1_DIAGNOSTIC_N_LATENT")) STUDY1_DIAGNOSTIC_N_LATENT <- 64L
if (!exists("STUDY1_DIAGNOSTIC_N_POINTS")) STUDY1_DIAGNOSTIC_N_POINTS <- 5L
if (!exists("STUDY1_STRICT_COMPETITORS")) {
  STUDY1_STRICT_COMPETITORS <- identical(STUDY1_CONFIG, "thorough")
}
if (!exists("STUDY1_PUBLISHED_COMPETITORS")) {
  STUDY1_PUBLISHED_COMPETITORS <- c("UC-GP", "LVGP", "EzGP")
}
if (!exists("STUDY1_MECHANISM_CALIB")) STUDY1_MECHANISM_CALIB <- 20L
if (!exists("STUDY1_DATA_DIR")) {
  STUDY1_DATA_DIR <- file.path("..", "data-synthetic", "study1")
}

settings <- switch(
  STUDY1_CONFIG,
  quick = list(
    n_test = 150L, n_rep = 2L, n_iter = 600L, burn = 200L,
    n_chains = 1L, preset = "fast", n_pred_draw = 150L,
    n_m_eval = 40L, n_m_draw = 60L, n_m_latent = 64L
  ),
  balanced = list(
    n_test = 300L, n_rep = 10L, n_iter = 2500L, burn = 750L,
    n_chains = 4L, preset = "balanced", n_pred_draw = 300L,
    n_m_eval = 80L, n_m_draw = 150L, n_m_latent = 128L
  ),
  thorough = list(
    n_test = 500L, n_rep = 50L, n_iter = 5000L, burn = 1000L,
    n_chains = 12L, preset = "balanced", n_pred_draw = 600L,
    n_m_eval = 120L, n_m_draw = 200L, n_m_latent = 256L
  )
)

## Publication masters may override the run size without changing the
## scientific engine below. Keeping these controls explicit also prevents a
## smoke run from being mistaken for a reportable Monte Carlo experiment.
study1_setting_override <- function(name, fallback, minimum = 1L) {
  value <- get0(name, inherits = TRUE, ifnotfound = fallback)
  value <- as.integer(value)
  if (length(value) != 1L || is.na(value) || value < minimum) {
    stop(name, " must be one integer at least ", minimum, ".")
  }
  value
}
settings$n_test <- study1_setting_override(
  "STUDY1_MC_N_TEST", settings$n_test
)
settings$n_rep <- study1_setting_override(
  "STUDY1_MC_N_REP", settings$n_rep
)
settings$n_iter <- study1_setting_override(
  "STUDY1_MC_N_ITER", settings$n_iter, minimum = 2L
)
settings$burn <- study1_setting_override(
  "STUDY1_MC_BURN", settings$burn, minimum = 0L
)
settings$n_chains <- study1_setting_override(
  "STUDY1_MC_N_CHAINS", settings$n_chains
)
settings$n_pred_draw <- study1_setting_override(
  "STUDY1_MC_N_PRED_DRAW", settings$n_pred_draw
)
settings$n_m_eval <- study1_setting_override(
  "STUDY1_MC_N_M_EVAL", settings$n_m_eval
)
settings$n_m_draw <- study1_setting_override(
  "STUDY1_MC_N_M_DRAW", settings$n_m_draw
)
settings$n_m_latent <- study1_setting_override(
  "STUDY1_MC_N_M_LATENT", settings$n_m_latent
)
if (settings$burn >= settings$n_iter) {
  stop("STUDY1_MC_BURN must be smaller than STUDY1_MC_N_ITER.")
}

if (!exists("STUDY1_SCENARIO")) {
  STUDY1_SCENARIO <- "heterogeneity_continuum"
}
if (!STUDY1_SCENARIO %in% c(
  "active", "inactive", "category_sufficient", "heterogeneity_continuum"
)) {
  stop("STUDY1_SCENARIO is not supported by simulate_1d_data().")
}
if (!exists("STUDY1_HETEROGENEITY_ETA")) STUDY1_HETEROGENEITY_ETA <- 1
STUDY1_HETEROGENEITY_ETA <- as.numeric(STUDY1_HETEROGENEITY_ETA)
if (length(STUDY1_HETEROGENEITY_ETA) != 1L ||
    !is.finite(STUDY1_HETEROGENEITY_ETA) ||
    STUDY1_HETEROGENEITY_ETA < 0) {
  stop("STUDY1_HETEROGENEITY_ETA must be one finite nonnegative number.")
}
if (!exists("STUDY1_THRESHOLD_DESIGN")) STUDY1_THRESHOLD_DESIGN <- "balanced"
STUDY1_THRESHOLD_DESIGN <- match.arg(
  STUDY1_THRESHOLD_DESIGN, c("balanced", "imbalanced")
)

FIG_DIR <- file.path(STUDY1_OUT_PREFIX, "figures")
TAB_DIR <- file.path(STUDY1_OUT_PREFIX, "tables")
RES_DIR <- file.path(STUDY1_OUT_PREFIX, "results", "study1_publication")
REP_DIR <- file.path(RES_DIR, "competitor_replications")
for (dd in c(FIG_DIR, TAB_DIR, RES_DIR, REP_DIR)) {
  dir.create(dd, showWarnings = FALSE, recursive = TRUE)
}

set.seed(20260705)
n_train <- study1_setting_override("STUDY1_MC_N_TRAIN", 100L)
n_test <- settings$n_test
n_rep <- settings$n_rep
m <- study1_setting_override("STUDY1_MC_M", 6L, minimum = 2L)
calib_grid <- as.integer(get0(
  "STUDY1_CALIB_GRID", inherits = TRUE,
  ifnotfound = c(0L, 5L, 10L, 20L, 50L)
))
if (length(calib_grid) < 1L || anyNA(calib_grid) ||
    any(calib_grid < 0L | calib_grid > n_train) || anyDuplicated(calib_grid)) {
  stop("STUDY1_CALIB_GRID must contain unique sizes between zero and n_train.")
}
calib_grid <- sort(calib_grid)
n_pred_draw <- settings$n_pred_draw
if (!exists("STUDY1_PARALLEL_LEVEL")) STUDY1_PARALLEL_LEVEL <- "chains"
if (!STUDY1_PARALLEL_LEVEL %in% c("chains", "replications", "none", "hybrid")) {
  stop("STUDY1_PARALLEL_LEVEL must be chains, replications, or none.")
}
parallel_chains <- STUDY1_PARALLEL_LEVEL %in% c("chains", "hybrid")
chain_cores <- if (STUDY1_PARALLEL_LEVEL == "hybrid") {
  mixedgp_as_integer_strict(STUDY1_CHAIN_WORKERS, "chain workers", 1L, 1L)
} else NULL
replication_cores <- if (STUDY1_PARALLEL_LEVEL == "hybrid") {
  min(n_rep, mixedgp_as_integer_strict(STUDY1_DATASET_WORKERS, "dataset workers", 1L, 1L))
} else if (identical(STUDY1_PARALLEL_LEVEL, "replications")) {
  min(n_rep, mixedgp_resolve_cores())
} else {
  1L
}

STUDY1_DESIGN_TAG <- paste0(
  "study1_", STUDY1_SCENARIO,
  "_eta", gsub("[^0-9]+", "p", format(STUDY1_HETEROGENEITY_ETA)),
  "_", STUDY1_THRESHOLD_DESIGN,
  "_n", n_train, "_m", m,
  "_audited_sampler_v7_nlpd_all_tasks_",
  STUDY1_CONFIG,
  "_meval", settings$n_m_eval,
  "_mdraw", settings$n_m_draw,
  "_mint", settings$n_m_latent,
  "_", paste(STUDY1_PUBLISHED_COMPETITORS, collapse = "-")
)

method_levels <- c(
  "Oracle", "EIV-GP", "UC-GP", "LVGP", "EzGP"
)
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

study1_write_csv_safe <- function(x, path) {
  if (!is.data.frame(x)) stop("Study I CSV output must be a data frame.")
  if (ncol(x) == 0L) x <- data.frame(note = character(0))
  write.csv(x, path, row.names = FALSE)
  invisible(x)
}

study1_sampler_control_rows <- function(fit, rep_id, n_calib) {
  encode <- function(value) {
    if (is.numeric(value)) {
      return(paste(format(value, digits = 17L, scientific = FALSE), collapse = ";"))
    }
    paste(as.character(value), collapse = ";")
  }
  settings <- c(
    fit$control,
    list(
      initialization_rule = paste(
        "threshold-compatible latent initialization with independent",
        "seeded perturbations; see fit_eivgp_1d()"
      ),
      chain_seeds = fit$mcmc$chain_stats$seed,
      covariance_jitter = fit$diagnostics$summary$covariance_jitter,
      forms_explicit_covariance_inverse =
        fit$diagnostics$summary$forms_explicit_covariance_inverse
    )
  )
  data.frame(
    rep = as.integer(rep_id),
    n_calib = as.integer(n_calib),
    setting = names(settings),
    value = vapply(settings, encode, character(1L)),
    stringsAsFactors = FALSE
  )
}

if (!exists("STUDY1_LVGP_MAX_ELAPSED")) {
  STUDY1_LVGP_MAX_ELAPSED <- if (STUDY1_QUICK) 180 else 900
}
competitor_controls <- list(
  `UC-GP` = list(n_starts = if (STUDY1_QUICK) 2L else 8L),
  LVGP = list(
    n_starts = if (STUDY1_QUICK) 2L else 8L,
    max_retries = if (STUDY1_QUICK) 1L else 3L,
    max_iter_ini = if (STUDY1_QUICK) 30L else 100L,
    max_iter_lat = if (STUDY1_QUICK) 8L else 20L,
    rescue_iter_ini = 300L,
    rescue_iter_lat = 100L,
    max_elapsed_seconds = STUDY1_LVGP_MAX_ELAPSED,
    parallel = FALSE
  ),
  EzGP = list(
    tau_fractions = c(1e-6, 0.0025, 0.01, 0.04, 0.16),
    cv_folds = 3L,
    maxeval = if (STUDY1_QUICK) 30L else 100L
  )
)

measurement_n_iter <- if (STUDY1_QUICK) 500L else if (
  STUDY1_CONFIG == "balanced"
) 1500L else 3000L
measurement_burn <- if (STUDY1_QUICK) 150L else if (
  STUDY1_CONFIG == "balanced"
) 500L else 1000L
measurement_n_iter <- study1_setting_override(
  "STUDY1_MEAS_N_ITER", measurement_n_iter, minimum = 2L
)
measurement_burn <- study1_setting_override(
  "STUDY1_MEAS_BURN", measurement_burn, minimum = 0L
)
if (measurement_burn >= measurement_n_iter) {
  stop("STUDY1_MEAS_BURN must be smaller than STUDY1_MEAS_N_ITER.")
}
STUDY1_DESIGN_TAG <- paste0(
  STUDY1_DESIGN_TAG,
  "_nrep", n_rep,
  "_ntest", n_test,
  "_iter", settings$n_iter,
  "_burn", settings$burn,
  "_chains", settings$n_chains,
  "_pred", n_pred_draw,
  "_cal", paste(calib_grid, collapse = "-"),
  "_meas", measurement_n_iter, "b", measurement_burn,
  "_f", as.integer(isTRUE(STUDY1_EVALUATE_F)),
  "_u", as.integer(isTRUE(STUDY1_EVALUATE_U))
)
STUDY1_CACHE_SPEC <- list(
  schema = "s1v9_nonfatal_target_diagnostics",
  design_label = STUDY1_DESIGN_TAG,
  study1_config = STUDY1_CONFIG,
  scenario = STUDY1_SCENARIO,
  heterogeneity_eta = STUDY1_HETEROGENEITY_ETA,
  threshold_design = STUDY1_THRESHOLD_DESIGN,
  n_train = n_train,
  n_test = n_test,
  n_rep = n_rep,
  m = m,
  calibration_grid = calib_grid,
  mcmc = list(
    n_iter = settings$n_iter, burn = settings$burn,
    n_chains = settings$n_chains, preset = settings$preset,
    rhat_limit = STUDY1_MAX_RHAT, ess_limit = STUDY1_MIN_ESS
  ),
  diagnostic_panel = list(
    max_draws_per_chain = STUDY1_DIAGNOSTIC_MAX_DRAWS,
    n_latent = STUDY1_DIAGNOSTIC_N_LATENT,
    n_points = STUDY1_DIAGNOSTIC_N_POINTS
  ),
  prediction = list(n_draw = n_pred_draw),
  mean_recovery = list(
    n_eval = settings$n_m_eval,
    n_draw = settings$n_m_draw,
    n_latent = settings$n_m_latent
  ),
  measurement = list(n_iter = measurement_n_iter, burn = measurement_burn),
  task_flags = list(
    run_ablations = isTRUE(STUDY1_RUN_ABLATIONS),
    evaluate_f = isTRUE(STUDY1_EVALUATE_F),
    evaluate_u = isTRUE(STUDY1_EVALUATE_U)
  ),
  published_competitors = STUDY1_PUBLISHED_COMPETITORS,
  competitor_controls = competitor_controls
)
STUDY1_CACHE_SPEC$fingerprint <-
  mixedgp_object_fingerprint(STUDY1_CACHE_SPEC)
STUDY1_DESIGN_TAG <- paste0(
  "s1v7-", gsub("[^0-9A-Za-z]+", "-", STUDY1_SCENARIO),
  "-eta", gsub("[^0-9A-Za-z]+", "p", STUDY1_HETEROGENEITY_ETA),
  "-", substr(STUDY1_CACHE_SPEC$fingerprint, 1L, 16L)
)
saveRDS(
  STUDY1_CACHE_SPEC,
  file.path(RES_DIR, paste0("cache_spec_", STUDY1_DESIGN_TAG, ".rds")),
  version = 3L
)
writeLines(
  capture.output(dput(STUDY1_CACHE_SPEC)),
  file.path(RES_DIR, paste0("cache_spec_", STUDY1_DESIGN_TAG, ".R"))
)

preflight <- mixedgp_competitor_preflight(
  STUDY1_PUBLISHED_COMPETITORS,
  strict = STUDY1_STRICT_COMPETITORS
)
write.csv(
  preflight,
  file.path(TAB_DIR, "study1_competitor_preflight.csv"),
  row.names = FALSE
)
