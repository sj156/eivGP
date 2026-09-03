############################################################
## 02_study1_monte_carlo.R
##
## Publication comparison for Study I.
##
## The publication path refits every method, including EIV-GP, with the current
## audited sampler. The historical July 27 rows remain available only through
## the explicit archival switch STUDY1_REUSE_LOCKED_EIV=TRUE.
############################################################

if (!exists("fit_eivgp_1d")) source("00_study1_functions.R")
if (!exists("load_mixedgp_synthetic_dataset")) source("00_synthetic_data.R")
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
if (!exists("STUDY1_MAX_RHAT")) STUDY1_MAX_RHAT <- 1.05
if (!exists("STUDY1_MIN_ESS")) STUDY1_MIN_ESS <- 100
if (!exists("STUDY1_MCMC_PILOT_REPS")) STUDY1_MCMC_PILOT_REPS <- 0L
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
    n_chains = 12L, preset = "thorough", n_pred_draw = 600L,
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
MCMC_PROGRESS_DIR <- file.path(RES_DIR, "mcmc_progress")
for (dd in c(FIG_DIR, TAB_DIR, RES_DIR, REP_DIR, MCMC_PROGRESS_DIR)) {
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
if (!STUDY1_PARALLEL_LEVEL %in% c("chains", "replications", "none")) {
  stop("STUDY1_PARALLEL_LEVEL must be chains, replications, or none.")
}
parallel_chains <- identical(STUDY1_PARALLEL_LEVEL, "chains")
replication_cores <- if (identical(STUDY1_PARALLEL_LEVEL, "replications")) {
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
  schema = "s1v7_nlpd_alltasks",
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
  strict = FALSE
)
write.csv(
  preflight,
  file.path(TAB_DIR, "study1_competitor_preflight.csv"),
  row.names = FALSE
)

run_one_study1_replication <- function(rep_id, run_eiv = TRUE) {
  data_path <- file.path(
    STUDY1_DATA_DIR,
    mixedgp_dataset_filename(
      "study1", rep_id, STUDY1_SCENARIO, n_train, n_test, m
    )
  )
  if (!file.exists(data_path)) {
    stop(
      "Missing frozen Study I dataset: ", data_path,
      ". Run 12_generate_synthetic_datasets.R first."
    )
  }
  frozen <- load_mixedgp_synthetic_dataset_strict(
    data_path,
    expected = list(
      study = "study1", scenario = STUDY1_SCENARIO, rep_id = rep_id,
      n = n_train, n_test = n_test, m = m,
      calib_grid = calib_grid,
      threshold_design = STUDY1_THRESHOLD_DESIGN,
      heterogeneity_eta = STUDY1_HETEROGENEITY_ETA
    )
  )
  dat <- frozen$data
  train <- dat$train
  test <- dat$test
  calib_sets <- frozen$calibration_sets

  oracle_draws <- sample_oracle_test_y(
    x_test = test$x,
    c_test = test$c,
    tau_true = dat$tau_true,
    scenario = STUDY1_SCENARIO,
    sigma_eps = dat$sigma_eps,
    n_draw = n_pred_draw,
    heterogeneity_eta = STUDY1_HETEROGENEITY_ETA
  )
  metrics <- list(
    Oracle = summarize_predictive_samples_1d(
      oracle_draws, test$y, method = "Oracle", rep_id = rep_id,
      n_calib = NA_integer_, scenario = STUDY1_SCENARIO
    )
  )
  set.seed(240000L + rep_id)
  mean_eval_idx <- sort(sample(
    seq_len(nrow(test)), min(settings$n_m_eval, nrow(test))
  ))
  m_true <- m0_1d(
    test$x[mean_eval_idx], test$c[mean_eval_idx],
    tau = dat$tau_true, scenario = STUDY1_SCENARIO
  )
  mean_metrics <- list()
  ablation_metrics <- list()
  latent_imputation_metrics <- list()
  surface_recovery_metrics <- list()
  ablation_surface_metrics <- list()
  ablation_status <- list()
  measurement_fits <- list()
  mcmc_diagnostics <- list()
  sampler_controls <- list()
  sampler_control_manifest <- list()

  competitor_result <- run_study1_published_competitors(
    X_train = matrix(train$x, ncol = 1L),
    y_train = train$y,
    C_train = matrix(train$c, ncol = 1L),
    X_test = matrix(test$x, ncol = 1L),
    C_test = matrix(test$c, ncol = 1L),
    m_vec = m,
    n_draw = n_pred_draw,
    seed = 250000L + 1000L * rep_id,
    methods = STUDY1_PUBLISHED_COMPETITORS,
    ## Published packages are external numerical optimizers.  Their individual
    ## failures are recorded in competitor_status, but must not prevent the
    ## EIV-GP MCMC from being run for this frozen replication.
    strict = FALSE,
    controls = competitor_controls
  )
  for (method in names(competitor_result$draws)) {
    metrics[[method]] <- summarize_predictive_samples_1d(
      competitor_result$draws[[method]], test$y,
      method = method, rep_id = rep_id,
      n_calib = NA_integer_, scenario = STUDY1_SCENARIO
    )
    mean_metrics[[method]] <- summarize_mean_recovery_1d(
      matrix(
        competitor_result$latent_means[[method]][mean_eval_idx],
        nrow = 1L
      ),
      m_true = m_true,
      method = method,
      rep_id = rep_id,
      n_calib = NA_integer_,
      scenario = STUDY1_SCENARIO,
      valid_function_draws = FALSE
    )
  }

  if (isTRUE(STUDY1_RUN_ABLATIONS)) {
    if (isTRUE(STUDY1_EVALUATE_F)) {
      start_full_u <- proc.time()[3]
      full_u_fit <- tryCatch(
        fit_study1_full_u_gp(train$x, train$y, train$u),
        error = function(e) e
      )
      elapsed_full_u <- proc.time()[3] - start_full_u
      if (inherits(full_u_fit, "error")) {
        ablation_status[["Full-U GP"]] <- data.frame(
          method = "Full-U GP", n_calib = NA_integer_, status = "failed",
          message = conditionMessage(full_u_fit),
          elapsed_seconds = elapsed_full_u
        )
      } else {
        ablation_status[["Full-U GP"]] <- data.frame(
          method = "Full-U GP", n_calib = NA_integer_, status = "success",
          message = "", elapsed_seconds = elapsed_full_u
        )
        ablation_surface_metrics[["Full-U GP"]] <-
          study1_latent_gp_surface_metrics(
            gp_fit = full_u_fit,
            scenario = STUDY1_SCENARIO,
            method = "Full-U GP",
            rep_id = rep_id,
            n_calib = NA_integer_,
            tau_true = dat$tau_true,
            heterogeneity_eta = STUDY1_HETEROGENEITY_ETA,
            grid_n_x = if (STUDY1_QUICK) 15L else 31L,
            grid_n_u = if (STUDY1_QUICK) 21L else 51L
          )
      }
    }
    for (n_calib in calib_grid) {
      key <- as.character(n_calib)
      start_measurement <- proc.time()[3]
      measurement_fit <- tryCatch(
        fit_threshold_measurement_response_free(
          c_ord = train$c,
          u_obs = train$u,
          calib_idx = calib_sets[[key]],
          m = m,
          n_iter = measurement_n_iter,
          burn = measurement_burn,
          thin = 1L,
          seed = 270000L + 1000L * rep_id + n_calib
        ),
        error = function(e) e
      )
      elapsed_measurement <- proc.time()[3] - start_measurement
      if (inherits(measurement_fit, "error")) {
        msg <- conditionMessage(measurement_fit)
        ablation_status[[paste0("measurement_", key)]] <- data.frame(
          method = "Response-free threshold model", n_calib = n_calib,
          status = "failed", message = msg,
          elapsed_seconds = elapsed_measurement
        )
        ablation_status[[paste0("PI_", key)]] <- data.frame(
          method = "PI-GP", n_calib = n_calib,
          status = "failed", message = msg, elapsed_seconds = NA_real_
        )
        ablation_status[[paste0("CC_", key)]] <- data.frame(
          method = "CC-GP", n_calib = n_calib,
          status = if (n_calib == 0L) "not_applicable" else "failed",
          message = if (n_calib == 0L) "No complete cases." else msg,
          elapsed_seconds = NA_real_
        )
        next
      }
      ablation_status[[paste0("measurement_", key)]] <- data.frame(
        method = "Response-free threshold model", n_calib = n_calib,
        status = "success", message = "",
        elapsed_seconds = elapsed_measurement
      )
      measurement_fits[[key]] <- measurement_fit

      start_pi <- proc.time()[3]
      pi_result <- tryCatch({
        fit <- fit_study1_pi_gp(train$x, train$y, measurement_fit)
        draws <- sample_study1_pi_gp(
          fit, test$x, test$c, n_pred_draw,
          seed = 280000L + 1000L * rep_id + n_calib
        )
        list(fit = fit, draws = draws)
      }, error = function(e) e)
      elapsed_pi <- proc.time()[3] - start_pi
      if (inherits(pi_result, "error")) {
        ablation_status[[paste0("PI_", key)]] <- data.frame(
          method = "PI-GP", n_calib = n_calib,
          status = "failed", message = conditionMessage(pi_result),
          elapsed_seconds = elapsed_pi
        )
      } else {
        ablation_status[[paste0("PI_", key)]] <- data.frame(
          method = "PI-GP", n_calib = n_calib,
          status = "success", message = "", elapsed_seconds = elapsed_pi
        )
        ablation_metrics[[paste0("PI_", key)]] <- summarize_predictive_samples_1d(
          pi_result$draws, test$y, method = "PI-GP", rep_id = rep_id,
          n_calib = n_calib, scenario = STUDY1_SCENARIO
        )
        if (isTRUE(STUDY1_EVALUATE_F) && n_calib > 0L) {
          ablation_surface_metrics[[paste0("PI_", key)]] <-
            study1_latent_gp_surface_metrics(
              gp_fit = pi_result$fit,
              scenario = STUDY1_SCENARIO,
              method = "PI-GP",
              rep_id = rep_id,
              n_calib = n_calib,
              tau_true = dat$tau_true,
              heterogeneity_eta = STUDY1_HETEROGENEITY_ETA,
              grid_n_x = if (STUDY1_QUICK) 15L else 31L,
              grid_n_u = if (STUDY1_QUICK) 21L else 51L
            )
        }
      }

      if (n_calib == 0L) {
        ablation_status[[paste0("CC_", key)]] <- data.frame(
          method = "CC-GP", n_calib = n_calib,
          status = "not_applicable", message = "No complete cases.",
          elapsed_seconds = 0
        )
      } else {
        start_cc <- proc.time()[3]
        cc_result <- tryCatch({
          fit <- fit_study1_cc_gp(
            train$x, train$y, train$u, calib_sets[[key]], measurement_fit
          )
          draws <- sample_study1_cc_gp(
            fit, test$x, test$c, n_pred_draw,
            seed = 290000L + 1000L * rep_id + n_calib
          )
          list(fit = fit, draws = draws)
        }, error = function(e) e)
        elapsed_cc <- proc.time()[3] - start_cc
        if (inherits(cc_result, "error")) {
          ablation_status[[paste0("CC_", key)]] <- data.frame(
            method = "CC-GP", n_calib = n_calib,
            status = "failed", message = conditionMessage(cc_result),
            elapsed_seconds = elapsed_cc
          )
        } else {
          ablation_status[[paste0("CC_", key)]] <- data.frame(
            method = "CC-GP", n_calib = n_calib,
            status = "success", message = "", elapsed_seconds = elapsed_cc
          )
          ablation_metrics[[paste0("CC_", key)]] <- summarize_predictive_samples_1d(
            cc_result$draws, test$y, method = "CC-GP", rep_id = rep_id,
            n_calib = n_calib, scenario = STUDY1_SCENARIO
          )
          if (isTRUE(STUDY1_EVALUATE_F)) {
            ablation_surface_metrics[[paste0("CC_", key)]] <-
              study1_latent_gp_surface_metrics(
                gp_fit = cc_result$fit,
                scenario = STUDY1_SCENARIO,
                method = "CC-GP",
                rep_id = rep_id,
                n_calib = n_calib,
                tau_true = dat$tau_true,
                heterogeneity_eta = STUDY1_HETEROGENEITY_ETA,
                grid_n_x = if (STUDY1_QUICK) 15L else 31L,
                grid_n_u = if (STUDY1_QUICK) 21L else 51L
              )
          }
        }
      }
    }
  }

  if (isTRUE(run_eiv)) {
    for (n_calib in calib_grid) {
      message("Study I replication ", rep_id, ": EIV-GP |O|=", n_calib)
      fit_eiv <- fit_eivgp_1d(
        x_raw = train$x,
        y_raw = train$y,
        c_ord = train$c,
        u_true = train$u,
        calib_idx = calib_sets[[as.character(n_calib)]],
        m = m,
        tau_true = dat$tau_true,
        n_iter = settings$n_iter,
        burn = settings$burn,
        thin = 1L,
        n_chains = settings$n_chains,
        preset = settings$preset,
        seed = 300000L + 1000L * rep_id + n_calib,
        parallel_chains = parallel_chains,
        verbose = FALSE,
        progress_every = 1000L,
        progress_label = sprintf("Study I rep %03d |O|=%d", rep_id, n_calib),
        progress_file = file.path(
          MCMC_PROGRESS_DIR,
          sprintf("rep%03d_calib%03d.log", rep_id, n_calib)
        )
      )
      control_key <- as.character(n_calib)
      sampler_controls[[control_key]] <- list(
        control = fit_eiv$control,
        initialization_rule = paste(
          "threshold-compatible latent initialization with independent",
          "seeded perturbations; see fit_eivgp_1d()"
        ),
        chain_seeds = fit_eiv$mcmc$chain_stats$seed,
        covariance_jitter = fit_eiv$diagnostics$summary$covariance_jitter,
        forms_explicit_covariance_inverse =
          fit_eiv$diagnostics$summary$forms_explicit_covariance_inverse
      )
      sampler_control_manifest[[control_key]] <-
        study1_sampler_control_rows(fit_eiv, rep_id, n_calib)
      diag_row <- fit_eiv$diagnostics$summary
      rhat_values <- unlist(
        diag_row[c("max_rhat_hyper", "max_rhat_tau", "max_rhat_missing_u")],
        use.names = FALSE
      )
      rhat_values <- rhat_values[is.finite(rhat_values)]
      max_rhat <- if (length(rhat_values) > 0L) max(rhat_values) else NA_real_
      min_ess <- as.numeric(diag_row$min_ess_key)
      ## With no calibrated latent values, the raw latent scale and threshold
      ## coordinates are not identified.  Their parameter-level R-hat/ESS is
      ## therefore descriptive only; publication gating starts once the
      ## calibration set anchors the raw latent scale.
      gate_applicable <- isTRUE(fit_eiv$data$latent_scale_anchored)
      gate_pass <- if (gate_applicable) {
        is.finite(max_rhat) && max_rhat <= STUDY1_MAX_RHAT &&
          is.finite(min_ess) && min_ess >= STUDY1_MIN_ESS
      } else {
        NA
      }
      diag_row$rep <- rep_id
      diag_row$n_calib <- n_calib
      diag_row$gate_max_rhat <- max_rhat
      diag_row$gate_min_ess <- min_ess
      diag_row$gate_applicable <- gate_applicable
      diag_row$gate_pass <- gate_pass
      mcmc_diagnostics[[as.character(n_calib)]] <- diag_row
      draw_ids <- seq_len(nrow(fit_eiv$mcmc$samples_u))
      if (length(draw_ids) > n_pred_draw) {
        set.seed(350000L + 1000L * rep_id + n_calib)
        draw_ids <- sample(draw_ids, n_pred_draw)
      }
      eiv_draws <- sample_eiv_test_y(
        x_test_raw = test$x,
        c_test = test$c,
        fit_obj = fit_eiv,
        draw_ids = draw_ids,
        n_per_draw = 1L
      )
      metrics[[paste0("EIV_", n_calib)]] <- summarize_predictive_samples_1d(
        eiv_draws, test$y, method = "EIV-GP", rep_id = rep_id,
        n_calib = n_calib, scenario = STUDY1_SCENARIO
      )
      m_draw_ids <- draw_ids
      if (length(m_draw_ids) > settings$n_m_draw) {
        set.seed(360000L + 1000L * rep_id + n_calib)
        m_draw_ids <- sample(m_draw_ids, settings$n_m_draw)
      }
      m_draws <- sample_eiv_m_given_xc_1d(
        x_star_raw = test$x[mean_eval_idx],
        c_star = test$c[mean_eval_idx],
        fit_obj = fit_eiv,
        draw_ids = m_draw_ids,
        n_latent = settings$n_m_latent,
        include_process_uncertainty = TRUE,
        joint = FALSE,
        return_components = TRUE,
        seed = 370000L + 1000L * rep_id + n_calib
      )
      mean_metrics[[paste0("EIV_", n_calib)]] <-
        summarize_mean_recovery_1d(
          m_draws,
          m_true = m_true,
          method = "EIV-GP",
          rep_id = rep_id,
          n_calib = n_calib,
          scenario = STUDY1_SCENARIO,
          valid_function_draws = TRUE
        )
      if (isTRUE(STUDY1_EVALUATE_U) &&
          !is.null(measurement_fits[[as.character(n_calib)]])) {
        latent_imputation_metrics[[as.character(n_calib)]] <-
          study1_latent_imputation_metrics(
            eiv_fit = fit_eiv,
            measurement_fit = measurement_fits[[as.character(n_calib)]],
            U_true = train$u,
            rep_id = rep_id,
            n_calib = n_calib,
            scenario = STUDY1_SCENARIO
          )
      }
      if (isTRUE(STUDY1_EVALUATE_F) && n_calib > 0L) {
        surface_recovery_metrics[[as.character(n_calib)]] <-
          study1_eiv_surface_recovery_metrics(
            fit = fit_eiv,
            scenario = STUDY1_SCENARIO,
            rep_id = rep_id,
            n_calib = n_calib,
            tau_true = dat$tau_true,
            heterogeneity_eta = STUDY1_HETEROGENEITY_ETA,
            grid_n_x = if (STUDY1_QUICK) 15L else 31L,
            grid_n_u = if (STUDY1_QUICK) 21L else 51L,
            max_draw = if (STUDY1_QUICK) 40L else 200L,
            seed = 380000L + 1000L * rep_id + n_calib
          )
      }
    }
  }

  status <- competitor_result$status
  status$rep <- rep_id
  status$scenario <- STUDY1_SCENARIO
  ablation_status_df <- bind_rows(ablation_status)
  if (nrow(ablation_status_df) > 0L) {
    ablation_status_df$rep <- rep_id
    ablation_status_df$scenario <- STUDY1_SCENARIO
  }
  list(
    metrics = bind_rows(metrics),
    competitor_status = status,
    ablation_metrics = bind_rows(ablation_metrics),
    ablation_status = ablation_status_df,
    mcmc_diagnostics = bind_rows(mcmc_diagnostics),
    sampler_control_manifest = bind_rows(sampler_control_manifest),
    mean_recovery = bind_rows(mean_metrics),
    latent_imputation = bind_rows(latent_imputation_metrics),
    surface_recovery = bind_rows(surface_recovery_metrics),
    ablation_surface_recovery = bind_rows(ablation_surface_metrics),
    metadata = list(
      rep = rep_id,
      design_tag = STUDY1_DESIGN_TAG,
      data_seed = 100000L + rep_id,
      calibration_seed = 200000L + rep_id,
      frozen_data_file = basename(data_path),
      frozen_data_md5 = attr(frozen, "manifest_md5"),
      frozen_manifest = attr(frozen, "manifest_path"),
      sampler = sampler_controls
    )
  )
}

run_study1_mcmc_pilot <- function(rep_id) {
  data_path <- file.path(
    STUDY1_DATA_DIR,
    mixedgp_dataset_filename(
      "study1", rep_id, STUDY1_SCENARIO, n_train, n_test, m
    )
  )
  frozen <- load_mixedgp_synthetic_dataset_strict(
    data_path,
    expected = list(
      study = "study1", scenario = STUDY1_SCENARIO, rep_id = rep_id,
      n = n_train, n_test = n_test, m = m,
      calib_grid = calib_grid,
      threshold_design = STUDY1_THRESHOLD_DESIGN,
      heterogeneity_eta = STUDY1_HETEROGENEITY_ETA
    )
  )
  dat <- frozen$data
  anchored_grid <- calib_grid[calib_grid > 0L]
  rows <- lapply(anchored_grid, function(n_calib) {
    tryCatch({
      fit <- fit_eivgp_1d(
        x_raw = dat$train$x, y_raw = dat$train$y, c_ord = dat$train$c,
        u_true = dat$train$u,
        calib_idx = frozen$calibration_sets[[as.character(n_calib)]],
        m = m, tau_true = dat$tau_true,
        n_iter = settings$n_iter, burn = settings$burn, thin = 1L,
        n_chains = settings$n_chains, preset = settings$preset,
        seed = 300000L + 1000L * rep_id + n_calib,
        parallel_chains = parallel_chains, verbose = FALSE,
        progress_every = 1000L,
        progress_label = sprintf(
          "Study I pilot rep %03d |O|=%d", rep_id, n_calib
        ),
        progress_file = file.path(
          MCMC_PROGRESS_DIR,
          sprintf("pilot_rep%03d_calib%03d.log", rep_id, n_calib)
        )
      )
      summary <- fit$diagnostics$summary
      rhat_values <- unlist(
        summary[c("max_rhat_hyper", "max_rhat_tau", "max_rhat_missing_u")],
        use.names = FALSE
      )
      rhat_values <- rhat_values[is.finite(rhat_values)]
      max_rhat <- if (length(rhat_values)) max(rhat_values) else NA_real_
      min_ess <- as.numeric(summary$min_ess_key)
      rhat_table <- bind_rows(
        transform(fit$diagnostics$rhat_hyper, family = "hyperparameter"),
        transform(fit$diagnostics$rhat_tau, family = "threshold"),
        transform(fit$diagnostics$rhat_u, family = "latent_u")
      )
      finite_rhat <- which(is.finite(rhat_table$rhat))
      worst_rhat_parameter <- if (length(finite_rhat)) {
        rhat_table$parameter[finite_rhat[which.max(rhat_table$rhat[finite_rhat])]]
      } else {
        NA_character_
      }
      finite_ess <- which(is.finite(fit$diagnostics$ess_key$ess))
      min_ess_parameter <- if (length(finite_ess)) {
        fit$diagnostics$ess_key$parameter[
          finite_ess[which.min(fit$diagnostics$ess_key$ess[finite_ess])]
        ]
      } else {
        NA_character_
      }
      diagnostic_calibrations <- unique(c(min(anchored_grid), max(anchored_grid)))
      if (rep_id == 1L && n_calib %in% diagnostic_calibrations) {
        diagnostic_plots <- plot_eivgp_mcmc_diagnostics(
          fit, max_draws = 1000L, max_lag = 60L, max_latent = 2L
        )
        diagnostic_stub <- sprintf(
          "study1_mcmc_pilot_rep%03d_calib%03d", rep_id, n_calib
        )
        ggplot2::ggsave(
          file.path(FIG_DIR, paste0(diagnostic_stub, "_trace.pdf")),
          diagnostic_plots$trace, width = 12, height = 9
        )
        ggplot2::ggsave(
          file.path(FIG_DIR, paste0(diagnostic_stub, "_acf.pdf")),
          diagnostic_plots$autocorrelation, width = 12, height = 9
        )
        ggplot2::ggsave(
          file.path(FIG_DIR, paste0(diagnostic_stub, "_rank.pdf")),
          diagnostic_plots$rank, width = 12, height = 9
        )
      }
      data.frame(
        rep = rep_id, n_calib = n_calib, status = "completed",
        sampler_strategy = summary$sampler_strategy,
        gate_max_rhat = max_rhat, gate_min_ess = min_ess,
        worst_rhat_parameter = worst_rhat_parameter,
        min_ess_parameter = min_ess_parameter,
        gate_pass = isTRUE(fit$data$latent_scale_anchored) &&
          is.finite(max_rhat) && max_rhat <= STUDY1_MAX_RHAT &&
          is.finite(min_ess) && min_ess >= STUDY1_MIN_ESS,
        message = "", stringsAsFactors = FALSE
      )
    }, error = function(e) data.frame(
      rep = rep_id, n_calib = n_calib, status = "failed",
      sampler_strategy = settings$preset,
      gate_max_rhat = NA_real_, gate_min_ess = NA_real_, gate_pass = FALSE,
      worst_rhat_parameter = NA_character_, min_ess_parameter = NA_character_,
      message = conditionMessage(e), stringsAsFactors = FALSE
    ))
  })
  bind_rows(rows)
}

if (STUDY1_MCMC_PILOT_REPS > 0L) {
  pilot_replications <- seq_len(min(STUDY1_MCMC_PILOT_REPS, n_rep))
  message(
    "Study I MCMC publication pilot: ", length(pilot_replications),
    " frozen replication(s), anchored calibration sizes only."
  )
  pilot_rows <- mixedgp_parallel_lapply(
    as.list(pilot_replications),
    function(rr) tryCatch(
      run_study1_mcmc_pilot(rr),
      error = function(e) data.frame(
        rep = rr, n_calib = NA_integer_, status = "failed",
        sampler_strategy = settings$preset,
        gate_max_rhat = NA_real_, gate_min_ess = NA_real_, gate_pass = FALSE,
        worst_rhat_parameter = NA_character_, min_ess_parameter = NA_character_,
        message = conditionMessage(e), stringsAsFactors = FALSE
      )
    ),
    n_cores = min(replication_cores, length(pilot_replications)),
    seeds = 890000L + pilot_replications,
    mc.preschedule = FALSE
  )
  pilot_results <- bind_rows(pilot_rows)
  write.csv(
    pilot_results,
    file.path(TAB_DIR, "study1_mcmc_pilot.csv"),
    row.names = FALSE
  )
  if (nrow(pilot_results) == 0L || any(!pilot_results$gate_pass)) {
    stop(
      "Study I MCMC publication pilot failed; see study1_mcmc_pilot.csv. ",
      "The full replication grid was not started."
    )
  }
}

run_or_load_study1_replication <- function(rr) {
  cache_mode <- if (isTRUE(STUDY1_REUSE_LOCKED_EIV)) "archival" else "fresh"
  rep_file <- file.path(
    REP_DIR,
    sprintf("study1_rep%03d_%s_%s.rds", rr, STUDY1_DESIGN_TAG, cache_mode)
  )
  if (isTRUE(STUDY1_USE_CACHE) && file.exists(rep_file)) {
    cached <- readRDS(rep_file)
    data_path <- file.path(
      STUDY1_DATA_DIR,
      mixedgp_dataset_filename(
        "study1", rr, STUDY1_SCENARIO, n_train, n_test, m
      )
    )
    cache_valid <- is.list(cached$metadata) &&
      identical(cached$metadata$rep, as.integer(rr)) &&
      identical(cached$metadata$design_tag, STUDY1_DESIGN_TAG) &&
      identical(cached$metadata$frozen_data_file, basename(data_path)) &&
      file.exists(data_path) &&
      identical(
        tolower(as.character(cached$metadata$frozen_data_md5)),
        tolower(unname(tools::md5sum(data_path)))
      )
    if (isTRUE(cache_valid)) return(cached)
    message("Ignoring incompatible Study I cache: ", rep_file)
  }
  message("Study I replication ", rr, " of ", n_rep)
  out <- run_one_study1_replication(
    rr,
    run_eiv = !isTRUE(STUDY1_REUSE_LOCKED_EIV)
  )
  saveRDS(out, rep_file)
  out
}
run_study1_replication_safely <- function(rr) {
  tryCatch(
    run_or_load_study1_replication(rr),
    error = function(e) structure(
      list(
        rep = as.integer(rr),
        status = "failed",
        message = conditionMessage(e)
      ),
      class = c("mixedgp_replication_failure", "list")
    )
  )
}
rep_objects <- mixedgp_parallel_lapply(
  as.list(seq_len(n_rep)),
  run_study1_replication_safely,
  n_cores = replication_cores,
  seeds = 910000L + seq_len(n_rep),
  mc.preschedule = FALSE
)

failed_replications <- vapply(
  rep_objects, inherits, logical(1L), what = "mixedgp_replication_failure"
)
replication_status <- data.frame(
  rep = seq_len(n_rep),
  status = ifelse(failed_replications, "failed", "success"),
  message = vapply(seq_len(n_rep), function(ii) {
    if (failed_replications[ii]) rep_objects[[ii]]$message else ""
  }, character(1L)),
  stringsAsFactors = FALSE
)
write.csv(
  replication_status,
  file.path(TAB_DIR, "study1_replication_status.csv"),
  row.names = FALSE
)
successful_rep_objects <- rep_objects[!failed_replications]
if (length(successful_rep_objects) == 0L) {
  stop("All Study I replications failed; see study1_replication_status.csv.")
}

new_results <- bind_rows(lapply(successful_rep_objects, `[[`, "metrics"))
competitor_status <- bind_rows(
  lapply(successful_rep_objects, `[[`, "competitor_status")
)
ablation_results <- bind_rows(lapply(successful_rep_objects, `[[`, "ablation_metrics"))
ablation_status <- bind_rows(lapply(successful_rep_objects, `[[`, "ablation_status"))
mcmc_diagnostics <- bind_rows(lapply(successful_rep_objects, `[[`, "mcmc_diagnostics"))
sampler_control_manifest <- bind_rows(
  lapply(successful_rep_objects, `[[`, "sampler_control_manifest")
)
mean_recovery <- bind_rows(lapply(successful_rep_objects, `[[`, "mean_recovery"))
latent_imputation <- bind_rows(lapply(successful_rep_objects, `[[`, "latent_imputation"))
surface_recovery <- bind_rows(lapply(successful_rep_objects, `[[`, "surface_recovery"))
ablation_surface_recovery <- bind_rows(
  lapply(successful_rep_objects, `[[`, "ablation_surface_recovery")
)
if (nrow(mcmc_diagnostics) > 0L) {
  write.csv(
    mcmc_diagnostics,
    file.path(TAB_DIR, "study1_mcmc_diagnostics.csv"),
    row.names = FALSE
  )
}
if (nrow(sampler_control_manifest) > 0L) {
  write.csv(
    sampler_control_manifest,
    file.path(TAB_DIR, "study1_sampler_control_manifest.csv"),
    row.names = FALSE
  )
}

if (isTRUE(STUDY1_REUSE_LOCKED_EIV)) {
  locked_file <- file.path(TAB_DIR, "study1_mc_raw_results.csv")
  if (!file.exists(locked_file)) {
    stop("Locked July 27 Study I results are missing: ", locked_file)
  }
  locked <- read.csv(locked_file, stringsAsFactors = FALSE)
  locked_eiv <- locked |>
    filter(method == "EIV-GP", rep %in% seq_len(n_rep))
  expected_rows <- n_rep * length(calib_grid)
  if (nrow(locked_eiv) != expected_rows ||
      !setequal(unique(locked_eiv$n_calib), calib_grid)) {
    stop("Locked EIV-GP rows do not match the archived 50-by-5 Study I design.")
  }
  ## Oracle values in the archived file are duplicated over calibration sizes.
  ## Keep one copy per replication for the revised table and curves.
  locked_oracle <- locked |>
    filter(method == "Oracle", rep %in% seq_len(n_rep)) |>
    group_by(rep, scenario, method) |>
    slice(1L) |>
    ungroup() |>
    mutate(n_calib = NA_integer_)
  mc_results <- bind_rows(
    locked_eiv,
    locked_oracle,
    new_results |> filter(!method %in% c("EIV-GP", "Oracle"))
  )
} else {
  mc_results <- new_results
}

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

message("Study I revised outputs written under ", normalizePath(STUDY1_OUT_PREFIX))
