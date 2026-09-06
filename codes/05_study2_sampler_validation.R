############################################################
## 05_study2_sampler_validation.R
##
## Release-gate validation of the exact interwoven Study II sampler.
## This is a computational validation script, not a predictive comparison.
############################################################

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) == 1L) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}

source(file.path(script_dir, "00_parallel_utils.R"))
source(file.path(script_dir, "00_study2_functions.R"))
source(file.path(script_dir, "simulation_helpers.R"))

read_int_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.integer(default))
  as.integer(value)
}

read_num_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.numeric(default))
  as.numeric(value)
}

read_bool_env <- function(name, default = TRUE) {
  value <- tolower(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) return(isTRUE(default))
  if (value %in% c("1", "true", "yes", "y")) return(TRUE)
  if (value %in% c("0", "false", "no", "n")) return(FALSE)
  stop(name, " must be true or false.")
}

n_train <- read_int_env("EIVGP_VALIDATION_N", 24L)
n_iter <- read_int_env("EIVGP_VALIDATION_ITER", 3000L)
burn <- read_int_env("EIVGP_VALIDATION_BURN", 1000L)
thin <- read_int_env("EIVGP_VALIDATION_THIN", 1L)
n_chains <- read_int_env("EIVGP_VALIDATION_CHAINS", 4L)
seed <- read_int_env("EIVGP_VALIDATION_SEED", 20260829L)
functional_draws_per_chain <- read_int_env(
"EIVGP_VALIDATION_FUNCTIONAL_DRAWS", ceiling((n_iter - burn) / thin)
)
m_latent_draws <- read_int_env("EIVGP_VALIDATION_M_LATENT_DRAWS", 64L)
grid_transition_iter <- read_int_env("EIVGP_VALIDATION_GRID_ITER", 6000L)
grid_transition_burn <- read_int_env("EIVGP_VALIDATION_GRID_BURN", 1000L)
rhat_limit <- read_num_env("EIVGP_VALIDATION_RHAT_MAX", 1.01)
ess_limit <- read_num_env("EIVGP_VALIDATION_ESS_MIN", 400)
target_ess_limit <- read_num_env("EIVGP_VALIDATION_TARGET_ESS_MIN", 400)
max_mcse_z <- read_num_env("EIVGP_VALIDATION_MAX_MCSE_Z", 4.5)
enforce_gate <- read_bool_env("EIVGP_VALIDATION_ENFORCE", TRUE)
parallel_chains <- read_bool_env("EIVGP_VALIDATION_PARALLEL", TRUE)
strategy_env <- Sys.getenv("EIVGP_VALIDATION_STRATEGIES", unset = "")
sampler_strategies <- if (nzchar(strategy_env)) {
  unique(trimws(strsplit(strategy_env, ",", fixed = TRUE)[[1L]]))
} else {
  c("legacy", "interwoven")
}
calibration_env <- Sys.getenv("EIVGP_VALIDATION_CALIB_COUNTS", unset = "")

if (burn >= n_iter || thin < 1L || n_chains < 2L) {
  stop("Validation requires burn < n_iter and at least two chains.")
}
if (length(sampler_strategies) < 1L ||
    any(!sampler_strategies %in% c("legacy", "interwoven"))) {
  stop("EIVGP_VALIDATION_STRATEGIES must contain legacy and/or interwoven.")
}
if (functional_draws_per_chain < 20L || m_latent_draws < 8L ||
    grid_transition_burn >= grid_transition_iter) {
  stop("Functional/grid validation settings are too small or inconsistent.")
}
if (!is.finite(rhat_limit) || rhat_limit <= 1 ||
    !is.finite(ess_limit) || ess_limit < 1 ||
    !is.finite(target_ess_limit) || target_ess_limit < 1 ||
    !is.finite(max_mcse_z) || max_mcse_z <= 0) {
  stop("Invalid release-gate threshold.")
}

############################################################
## Independent scalar checks
############################################################

## Check the score-marginal ordinal likelihood against direct probability
## calculations in a numerically benign case.
set.seed(seed)
U_check <- matrix(rnorm(12), 6, 2)
A_check <- matrix(c(1.2, 0, 0.3, 1.1, 0.8, 0.4, 0.4, 0.9), 4, 2,
                  byrow = TRUE)
tau_check <- replicate(4, c(-0.7, 0, 0.7), simplify = FALSE)
C_check <- matrix(sample(1:4, 24, replace = TRUE), 6, 4)
ll_stable <- ordinal_loglik_marginal(C_check, U_check, A_check, tau_check)
ll_direct <- 0
for (j in seq_len(ncol(C_check))) {
  mu <- as.numeric(U_check %*% A_check[j, ])
  lower <- c(-Inf, tau_check[[j]])[C_check[, j]] - mu
  upper <- c(tau_check[[j]], Inf)[C_check[, j]] - mu
  ll_direct <- ll_direct + sum(log(pnorm(upper) - pnorm(lower)))
}
stopifnot(isTRUE(all.equal(ll_stable, ll_direct, tolerance = 1e-10)))

## A q=1, d=1 collapsed latent-state transition has a one-dimensional target
## that can be evaluated independently on a dense grid. This directly checks
## the prior-reference score-collapsed ESS transition used in the interweaving.
run_grid_target_check <- function() {
  X <- matrix(c(-0.7, 0.9), ncol = 1L)
  y <- c(0.85, -0.35)
  C <- matrix(c(3L, 2L), ncol = 1L)
  A <- matrix(1.25, nrow = 1L, ncol = 1L)
  tau <- list(c(-0.45, 0.55))
  logtheta <- log(c(2.2, 0.8, 1.1))
  gp_prior <- make_gp_prior(p = 1L, d = 1L)
  fixed_u2 <- -0.25

  grid <- seq(-5.5, 5.5, length.out = 22001L)
  grid_log_density <- vapply(grid, function(value) {
    U <- matrix(c(value, fixed_u2), ncol = 1L)
    dnorm(value, log = TRUE) +
      gp_loglik_integrated_general(y, X, U, logtheta) +
      ordinal_loglik_marginal(C, U, A, tau)
  }, numeric(1))
  weight <- exp(grid_log_density - max(grid_log_density))
  weight <- weight / sum(weight)
  target_mean <- sum(grid * weight)
  target_var <- sum((grid - target_mean)^2 * weight)

  set.seed(seed + 91L)
  U_state <- matrix(c(0.8, fixed_u2), ncol = 1L)
  kept <- numeric(grid_transition_iter - grid_transition_burn)
  for (iter in seq_len(grid_transition_iter)) {
    update <- update_U_theta_ess_integrated_general(
      y = y,
      X = X,
      C = C,
      U_curr = U_state,
      logtheta_curr = logtheta,
      gp_prior = gp_prior,
      block_idx = 1L,
      reference = "prior",
      A = A,
      tau = tau,
      update_theta = FALSE,
      max_try = 5000L
    )
    U_state <- update$U
    if (iter > grid_transition_burn) {
      kept[iter - grid_transition_burn] <- U_state[1L, 1L]
    }
  }

  mean_ess <- as.numeric(posterior::ess_mean(kept))
  squared_error <- (kept - target_mean)^2
  variance_ess <- as.numeric(posterior::ess_mean(squared_error))
  mean_mcse <- as.numeric(posterior::mcse_mean(kept))
  variance_mcse <- as.numeric(posterior::mcse_mean(squared_error))
  sample_mean <- mean(kept)
  sample_var <- mean(squared_error)
  mean_z <- abs(sample_mean - target_mean) / max(mean_mcse, 1e-12)
  variance_z <- abs(sample_var - target_var) / max(variance_mcse, 1e-12)
  pass <- is.finite(mean_z) && is.finite(variance_z) &&
    mean_ess >= ess_limit && variance_ess >= ess_limit &&
    mean_z <= max_mcse_z && variance_z <= max_mcse_z

  data.frame(
    check = "q1_d1_collapsed_U_grid",
    target_mean = target_mean,
    sample_mean = sample_mean,
    target_variance = target_var,
    sample_variance = sample_var,
    mean_ess = mean_ess,
    variance_ess = variance_ess,
    mean_mcse_z = mean_z,
    variance_mcse_z = variance_z,
    pass = pass,
    stringsAsFactors = FALSE
  )
}

grid_target_table <- run_grid_target_check()

############################################################
## Legacy versus interwoven posterior comparison
############################################################

dat <- simulate_study2_data(
  n = n_train,
  n_test = 20L,
  seed = seed,
  scenario = "primary"
)

calib_grid <- if (nzchar(calibration_env)) {
  unique(as.integer(trimws(strsplit(
    calibration_env, ",", fixed = TRUE
  )[[1L]])))
} else {
  unique(as.integer(c(0L, round(n_train / 4), round(n_train / 2))))
}
if (length(calib_grid) < 1L || anyNA(calib_grid) ||
    any(calib_grid < 0L | calib_grid > n_train)) {
  stop("EIVGP_VALIDATION_CALIB_COUNTS must be integers from zero to n.")
}
calib_sets <- make_stratified_calibration_sets_2d(
  dat$train$C,
  calib_grid = calib_grid,
  seed = seed + 1L
)

evaluation_rows <- mixedgp_study2_diagnostic_rows(dat$train$X, dat$train$C, 5L)
X_eval <- dat$train$X[evaluation_rows, , drop = FALSE]
C_eval <- dat$train$C[evaluation_rows, , drop = FALSE]
U_eval <- dat$train$U[evaluation_rows, , drop = FALSE]

posterior_series <- function(fit, n_calib) {
  targets <- mixedgp_study2_target_series(
    fit, X_eval, C_eval, U = if (n_calib > 0L) U_eval else NULL,
    max_draws_per_chain = functional_draws_per_chain,
    n_latent = m_latent_draws, seed = seed + 100000L
  )
  out <- c(
    mixedgp_study2_raw_series(fit),
    mixedgp_study2_invariant_series(fit),
    targets
  )
  attr(out, "draw_window_complete") <- attr(targets, "draw_window_complete")
  out
}

summarize_series <- function(series, n_calib, strategy) {
  out <- mixedgp_summarize_diagnostic_series(
    series, rhat_limit, target_ess_limit, target_ess_limit
  )
  out$n_calib <- n_calib
  out$sampler_strategy <- strategy
  out
}

summarize_chain_series <- function(series, n_calib, strategy) {
  do.call(rbind, lapply(names(series), function(parameter) {
    chains <- series[[parameter]]
    do.call(rbind, lapply(seq_along(chains), function(chain_id) {
      values <- as.numeric(chains[[chain_id]])
      data.frame(
        n_calib = n_calib,
        sampler_strategy = strategy,
        parameter = parameter,
        chain_id = chain_id,
        n_saved = length(values),
        mean = mean(values),
        sd = sd(values),
        q05 = unname(quantile(values, 0.05)),
        median = median(values),
        q95 = unname(quantile(values, 0.95)),
        ess = ess_ips(values),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

fit_summaries <- list()
signature_summaries <- list()
chain_signature_summaries <- list()
row_id <- 0L

for (k in calib_grid) {
  calib_idx <- calib_sets[[as.character(k)]]

  for (strategy in sampler_strategies) {
    fit <- fit_eivgp_ordprobit_fb(
      X_raw = dat$train$X,
      y_raw = dat$train$y,
      C_ord = dat$train$C,
      U_obs = dat$train$U,
      calib_idx = calib_idx,
      U_true_eval = dat$train$U,
      d = 2L,
      n_iter = n_iter,
      burn = burn,
      thin = thin,
      n_chains = n_chains,
      preset = "balanced",
      sampler_strategy = strategy,
      store_scores = FALSE,
      seed = seed + 100L * k + ifelse(strategy == "legacy", 1L, 2L),
      parallel_chains = parallel_chains,
      verbose = FALSE
    )

    ## Structural, calibration, and saved-state invariants.
    stopifnot(length(fit$data$calib_idx) == k)
    stopifnot(length(fit$data$miss_idx) == n_train - k)
    stopifnot(all(is.finite(fit$mcmc$samples_logtheta)))
    stopifnot(all(fit$mcmc$samples_sigma2 > 0))
    stopifnot(all(fit$mcmc$samples_A[, 1, 2] == 0))
    stopifnot(all(fit$mcmc$samples_A[, 1, 1] > 0))
    stopifnot(all(fit$mcmc$samples_A[, 2, 2] > 0))
    if (k > 0L) {
      for (ii in calib_idx) {
        stopifnot(all(fit$mcmc$samples_U[, ii, 1] == dat$train$U[ii, 1]))
        stopifnot(all(fit$mcmc$samples_U[, ii, 2] == dat$train$U[ii, 2]))
      }
    }

    tau_draws <- fit$mcmc$samples_tau
    for (j in seq_len(ncol(dat$train$C))) {
      tau_cols <- grep(paste0("^tau\\[", j, ","), colnames(tau_draws))
      stopifnot(all(apply(tau_draws[, tau_cols, drop = FALSE], 1, diff) > 0))
    }

    row_id <- row_id + 1L
    summary_row <- fit$diagnostics$summary
    summary_row$n_train <- n_train
    summary_row$n_calib <- k
    summary_row$gp_evaluations_per_chain_iteration <-
      summary_row$total_gp_evaluations / (n_chains * n_iter)
    summary_row$min_ess_per_1000_gp_evaluations <-
      1000 * summary_row$min_ess_key / summary_row$total_gp_evaluations
    summary_row$min_ess_per_second <-
      summary_row$min_ess_key / summary_row$time_seconds
    fit_summaries[[row_id]] <- summary_row

    series <- posterior_series(fit, n_calib = k)
    signature_summaries[[paste(k, strategy, sep = ":")]] <- summarize_series(
      series,
      n_calib = k,
      strategy = strategy
    )
    chain_signature_summaries[[paste(k, strategy, sep = ":")]] <-
      summarize_chain_series(series, n_calib = k, strategy = strategy)
  }
}

summary_table <- do.call(rbind, fit_summaries)
signature_table <- do.call(rbind, signature_summaries)
chain_signature_table <- do.call(rbind, chain_signature_summaries)

compare_strategies <- all(c("legacy", "interwoven") %in% sampler_strategies)
agreement_table <- if (compare_strategies) {
  do.call(rbind, lapply(calib_grid, function(k) {
    legacy <- signature_table[
      signature_table$n_calib == k &
        signature_table$sampler_strategy == "legacy",
    ]
    interwoven <- signature_table[
      signature_table$n_calib == k &
        signature_table$sampler_strategy == "interwoven",
    ]
    merged <- merge(legacy, interwoven, by = c("n_calib", "parameter"),
                    suffixes = c("_legacy", "_interwoven"))
    if (k == 0L) {
      ## Unanchored raw coordinates are descriptive, not an agreement target.
      merged <- merged[grepl(
        "^(ordinal_score_cor|standardized_tau|m_conditional_mean|predictive_second_moment)\\[",
        merged$parameter
      ), , drop = FALSE]
    }
    merged$absolute_mean_difference <- abs(
      merged$mean_legacy - merged$mean_interwoven
    )
    merged$combined_mcse <- sqrt(
      merged$mcse_legacy^2 + merged$mcse_interwoven^2
    )
    merged$mcse_z <- merged$absolute_mean_difference / merged$combined_mcse
    merged$agreement_pass <- is.finite(merged$mcse_z) &
      merged$mcse_z <= max_mcse_z
    merged
  }))
} else {
  data.frame()
}

diagnostic_gate <- do.call(rbind, lapply(seq_len(nrow(summary_table)), function(i) {
  row <- summary_table[i, , drop = FALSE]
  signatures <- signature_table[
    signature_table$n_calib == row$n_calib &
      signature_table$sampler_strategy == row$sampler_strategy, , drop = FALSE
  ]
  zero_calibration <- identical(as.integer(row$n_calib), 0L)
  invariant <- grepl("^(ordinal_score_cor|standardized_tau)\\[", signatures$parameter)
  target <- grepl(
    "^(m_conditional_mean|predictive_second_moment|f_conditional_mean|f_conditional_variance|prospective_U_mean|prospective_U_second_moment)\\[",
    signatures$parameter
  )
  raw <- !invariant & !target
  invariant_pass <- mixedgp_diagnostic_table_pass(
    signatures[invariant, , drop = FALSE], rhat_limit, target_ess_limit, target_ess_limit
  )
  target_pass <- mixedgp_diagnostic_table_pass(
    signatures[target, , drop = FALSE], rhat_limit, target_ess_limit, target_ess_limit
  )
  raw_pass <- mixedgp_diagnostic_table_pass(
    signatures[raw, , drop = FALSE], rhat_limit, ess_limit, ess_limit
  )
  selected <- signatures[target | invariant | (!zero_calibration & raw), , drop = FALSE]
  data.frame(
    n_calib = row$n_calib, sampler_strategy = row$sampler_strategy,
    selected_max_rhat = max(selected$rhat),
    selected_min_bulk_ess = min(selected$ess_bulk),
    selected_min_tail_ess = min(selected$ess_tail),
    selected_ess_limit = target_ess_limit,
    raw_coordinate_pass = raw_pass, invariant_measurement_pass = invariant_pass,
    target_functional_pass = target_pass,
    gate_basis = if (zero_calibration) {
      "response_moments_and_ordinal_measurement_invariants"
    } else {
      "all_reported_targets_measurement_invariants_and_free_coordinates"
    },
    required_for_release = row$sampler_strategy == "interwoven",
    rhat_limit = rhat_limit, ess_limit = ess_limit,
    pass = invariant_pass && target_pass && (zero_calibration || raw_pass),
    stringsAsFactors = FALSE
  )
}))

agreement_gate <- if (compare_strategies) {
  gate <- aggregate(
    cbind(max_mcse_z_observed = agreement_table$mcse_z,
          n_failed = as.integer(!agreement_table$agreement_pass)),
    by = list(n_calib = agreement_table$n_calib),
    FUN = function(x) if (all(is.finite(x))) max(x) else Inf
  )
  names(gate)[names(gate) == "n_failed"] <- "any_failed"
  gate$any_failed <- gate$any_failed > 0
  gate$mcse_z_limit <- max_mcse_z
  gate$pass <- !gate$any_failed & gate$max_mcse_z_observed <= max_mcse_z
  gate
} else {
  data.frame()
}

out_dir_env <- Sys.getenv("EIVGP_VALIDATION_OUT_DIR", unset = "")
out_dir <- if (nzchar(out_dir_env)) {
  dir.create(out_dir_env, recursive = TRUE, showWarnings = FALSE)
  normalizePath(out_dir_env, mustWork = TRUE)
} else {
  normalizePath(file.path(script_dir, "..", "tables"), mustWork = TRUE)
}
write.csv(
  summary_table,
  file.path(out_dir, "study2_sampler_validation_small.csv"),
  row.names = FALSE
)
write.csv(
  signature_table,
  file.path(out_dir, "study2_sampler_validation_signatures.csv"),
  row.names = FALSE
)
write.csv(
  chain_signature_table,
  file.path(out_dir, "study2_sampler_validation_chain_signatures.csv"),
  row.names = FALSE
)
write.csv(
  agreement_table,
  file.path(out_dir, "study2_sampler_validation_agreement_small.csv"),
  row.names = FALSE
)
write.csv(
  diagnostic_gate,
  file.path(out_dir, "study2_sampler_validation_diagnostic_gate.csv"),
  row.names = FALSE
)
write.csv(
  agreement_gate,
  file.path(out_dir, "study2_sampler_validation_agreement_gate.csv"),
  row.names = FALSE
)
write.csv(
  grid_target_table,
  file.path(out_dir, "study2_sampler_validation_grid_target.csv"),
  row.names = FALSE
)

print(summary_table)
print(diagnostic_gate)
print(agreement_gate)
print(grid_target_table)

agreement_pass <- !compare_strategies || all(agreement_gate$pass)
overall_pass <- isTRUE(grid_target_table$pass) &&
  any(diagnostic_gate$required_for_release) &&
  all(diagnostic_gate$pass[diagnostic_gate$required_for_release]) &&
  agreement_pass
if (!overall_pass) {
  message_text <- paste0(
    "Study II sampler release gate failed: grid=", grid_target_table$pass,
    ", required diagnostics=", paste(
      diagnostic_gate$pass[diagnostic_gate$required_for_release],
      collapse = "/"
    ),
    ", agreement=", if (compare_strategies) {
      paste(agreement_gate$pass, collapse = "/")
    } else {
      "not-run"
    },
    ". Inspect CSV files in ", out_dir, "."
  )
  if (isTRUE(enforce_gate)) stop(message_text)
  warning(message_text, call. = FALSE)
} else {
  cat("Study II sampler release gate passed.\n")
}
