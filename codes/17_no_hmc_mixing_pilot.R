## Reproducible, bounded comparison of exact no-HMC transition strategies.
## This is a computation pilot, never a publication simulation or certification.
## Run both studies: Rscript revision/codes/17_no_hmc_mixing_pilot.R
## Controls: EIVGP_MIX_N=24, _ITER=1200, _BURN=400, _CHAINS=4, _CORES=4,
## _CALIB=0,6, _STUDIES=1,2, _KERNELS=se, _LATENT_DRAWS=16,
## _TRANSPORT=true, _OUT_DIR=/absolute/new/directory, _SEED=20260906.
args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
code_dir <- dirname(normalizePath(sub("^--file=", "", args[1L])))
int_env <- function(name, default, minimum = 0L) {
  value <- suppressWarnings(as.numeric(Sys.getenv(paste0("EIVGP_MIX_", name), as.character(default))))
  if (length(value) != 1L || !is.finite(value) || value < minimum ||
      value != floor(value)) stop("Invalid EIVGP_MIX_", name)
  as.integer(value)
}
csv_env <- function(name, default) {
  unique(trimws(strsplit(Sys.getenv(paste0("EIVGP_MIX_", name), default), ",", fixed = TRUE)[[1L]]))
}
cfg <- list(n = int_env("N", 24L, 8L), n_iter = int_env("ITER", 1200L, 8L),
  burn = int_env("BURN", 400L), n_chains = int_env("CHAINS", 4L, 2L),
  n_cores = int_env("CORES", 4L, 1L), seed = int_env("SEED", 20260906L),
  n_latent = int_env("LATENT_DRAWS", 16L, 8L),
  studies = csv_env("STUDIES", "1,2"), kernels = csv_env("KERNELS", "se"),
  calibration = suppressWarnings(as.numeric(csv_env("CALIB", "0,6"))),
  transport = tolower(Sys.getenv("EIVGP_MIX_TRANSPORT", "true")))
if (cfg$n_iter - cfg$burn < 4L || length(cfg$calibration) < 1L ||
    length(cfg$studies) < 1L || length(cfg$kernels) < 1L ||
    any(!is.finite(cfg$calibration)) || any(cfg$calibration != floor(cfg$calibration)) ||
    any(cfg$calibration < 0L | cfg$calibration > cfg$n) ||
    any(!cfg$studies %in% c("1", "2")) || any(!cfg$kernels %in% c("se", "matern")) ||
    !cfg$transport %in% c("true", "false")) stop("Invalid mixing pilot configuration.")
cfg$transport <- cfg$transport == "true"
out <- Sys.getenv("EIVGP_MIX_OUT_DIR", "")
if (!nzchar(out)) out <- tempfile("no-hmc-mixing-", tmpdir = tempdir())
if (grepl("(^|/)(paper-ordinal|tables|figures|publication-results)(/|$)", out)) {
  stop("Use a fresh pilot directory outside publication artifacts.")
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out, mustWork = TRUE)
if (file.exists(file.path(out, "manifest.rds"))) stop("Pilot output already exists; choose a new directory.")
source_files <- file.path(code_dir, c("00_parallel_utils.R", "00_study1_functions.R",
  "00_study2_functions.R", "00_diagnostics.R", "simulation_helpers.R", "17_no_hmc_mixing_pilot.R"))
frozen_code_dir <- file.path(out, "source")
dir.create(frozen_code_dir)
if (!all(file.copy(source_files, frozen_code_dir, overwrite = FALSE))) {
  stop("Could not freeze pilot source files.")
}
saveRDS(list(config = cfg, purpose = "computation_pilot_not_publication_evidence",
  source_md5 = tools::md5sum(source_files),
  created = Sys.time(), session = sessionInfo()), file.path(out, "manifest.rds"))
cat("Pilot output:", out, "\n")

summary_rows <- list()
for (study in cfg$studies) {
  env <- new.env(parent = globalenv())
  for (module in c("00_parallel_utils.R", paste0("00_study", study, "_functions.R"),
                   "simulation_helpers.R")) source(file.path(frozen_code_dir, module), local = env)
  dat <- if (study == "1") env$simulate_1d_data(n = cfg$n, n_test = 4L,
    m = 4L, threshold_design = "balanced", min_class_count = 1L, seed = cfg$seed + 1L) else
    env$simulate_study2_data(n = cfg$n, n_test = 4L, seed = cfg$seed + 2L)
  saveRDS(dat, file.path(out, paste0("study", study, "-frozen-data.rds")))
  set.seed(cfg$seed + 301L)
  calibration_order <- sample.int(cfg$n)
  modes <- if (study == "1") c("conditional", "collapsed") else
    c("baseline", "joint_rows", if (cfg$transport) "joint_rows_transport")
  for (kernel in cfg$kernels) for (n_calib in cfg$calibration) for (mode in modes) {
    key <- paste0("study", study, "-", kernel, "-cal", n_calib, "-", mode)
    cat("Fitting", key, "\n"); flush.console()
    idx <- head(calibration_order, n_calib)
    warnings <- character()
    started <- proc.time()[3L]
    fit <- tryCatch(withCallingHandlers({
      common <- list(n_iter = cfg$n_iter, burn = cfg$burn, thin = 1L,
        n_chains = cfg$n_chains, preset = "balanced", seed = cfg$seed + 1000L * as.integer(study) + n_calib,
        parallel_chains = cfg$n_cores > 1L, n_cores = cfg$n_cores, kernel = kernel)
      if (study == "1") do.call(env$fit_eivgp_1d, c(list(x_raw = dat$train$x,
        y_raw = dat$train$y, c_ord = dat$train$c, u_true = dat$train$u,
        calib_idx = idx, m = 4L, noise_strategy = mode), common)) else
        do.call(env$fit_eivgp_ordprobit_fb, c(list(X_raw = dat$train$X,
          y_raw = dat$train$y, C_ord = dat$train$C, U_obs = dat$train$U,
          calib_idx = idx, d = 2L, m_vec = rep(4L, 4L), store_scores = FALSE,
          sampler_strategy = "interwoven", control_overrides = list(
            joint_measurement_every = if (mode == "baseline") 0L else 2L,
            loading_transport_every = if (mode == "joint_rows_transport") 5L else 0L)), common))
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = identity)
    seconds <- unname(proc.time()[3L] - started)
    if (inherits(fit, "error")) {
      saveRDS(list(error = conditionMessage(fit), warnings = warnings, seconds = seconds),
        file.path(out, paste0(key, "-FAILED.rds")))
      stop(key, ": ", conditionMessage(fit), "; completed fits preserved in ", out)
    }
    # Save before potentially costly target integration so fitting is recoverable.
    saveRDS(list(fit = fit, seconds = seconds, warnings = warnings), file.path(out, paste0(key, ".rds")))
    raw <- if (study == "1") env$mixedgp_study1_raw_series(fit) else env$mixedgp_study2_raw_series(fit)
    raw_diag <- env$mixedgp_summarize_diagnostic_series(raw)
    raw_diag$scope <- "raw"
    invariant <- if (study == "2") env$mixedgp_summarize_diagnostic_series(
      env$mixedgp_study2_invariant_series(fit)) else raw_diag[FALSE, setdiff(names(raw_diag), "scope"), drop = FALSE]
    invariant$scope <- rep("invariant", nrow(invariant))
    targets <- if (study == "1") env$mixedgp_study1_target_series(fit,
      matrix(dat$test$x[1:2], ncol = 1L), dat$test$c[1:2],
      U = if (n_calib) dat$test$u[1:2] else NULL, n_latent = cfg$n_latent, seed = cfg$seed + 700L) else
      env$mixedgp_study2_target_series(fit, dat$test$X[1:2, , drop = FALSE],
        dat$test$C[1:2, , drop = FALSE], U = if (n_calib) dat$test$U[1:2, , drop = FALSE] else NULL,
        n_latent = cfg$n_latent, seed = cfg$seed + 700L)
    target_diag <- env$mixedgp_summarize_diagnostic_series(targets)
    target_diag$scope <- "target"
    diagnostics <- rbind(raw_diag, invariant, target_diag)
    diagnostics$required <- diagnostics$scope != "raw" | study == "1" | n_calib > 0L
    diagnostics$ess_mean_per_second <- diagnostics$ess_mean / seconds
    diagnostics$ess_bulk_per_second <- diagnostics$ess_bulk / seconds
    write.csv(diagnostics, file.path(out, paste0(key, "-diagnostics.csv")), row.names = FALSE)
    required <- diagnostics[diagnostics$required, , drop = FALSE]
    pass <- env$mixedgp_diagnostic_table_pass(required)
    row <- data.frame(study = study, kernel = kernel, n = cfg$n, n_calib = n_calib,
      mode = mode, fit_seconds = seconds, max_rhat = max(required$rhat),
      min_bulk_ess = min(required$ess_bulk), min_tail_ess = min(required$ess_tail),
      min_mean_ess_per_second = min(required$ess_mean) / seconds,
      target_max_rhat = max(target_diag$rhat), target_min_bulk_ess = min(target_diag$ess_bulk),
      target_max_mcse_sd_ratio = max(target_diag$mcse_sd_ratio),
      strict_pass = pass, integration_draws = cfg$n_latent,
      reporting_window_complete = all(target_diag$draw_window_complete),
      publication_evidence = FALSE, warning_count = length(warnings))
    summary_rows[[length(summary_rows) + 1L]] <- row
    write.csv(do.call(rbind, summary_rows), file.path(out, "summary.csv"), row.names = FALSE)
    cat(key, ": seconds=", round(seconds, 2), "; max Rhat=", round(row$max_rhat, 3),
      "; min bulk ESS=", round(row$min_bulk_ess, 1), "; strict pass=", pass, "\n", sep = "")
    flush.console()
  }
}
cat("Finished. Interpret target ESS/time jointly with Rhat, MCSE and integration sensitivity.\n")
