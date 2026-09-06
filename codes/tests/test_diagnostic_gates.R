## Run from any directory: Rscript revision/codes/tests/test_diagnostic_gates.R
args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", args[1L]), mustWork = TRUE)
code_dir <- dirname(dirname(script))
project_library <- file.path(dirname(code_dir), "R-library")
if (dir.exists(project_library)) .libPaths(c(project_library, .libPaths()))
source(file.path(code_dir, "simulation_helpers.R"))
source(file.path(code_dir, "00_parallel_utils.R"))
source(file.path(code_dir, "00_study2_functions.R"))

set.seed(230905)
independent <- list(m = lapply(seq_len(4L), function(i) rnorm(2000L)))
good <- mixedgp_summarize_diagnostic_series(independent)
stopifnot(mixedgp_diagnostic_table_pass(good), good$target_pass)

## Old audit signatures had acceptable noisy-Y diagnostics while latent ordinal
## correlations mixed poorly. A release gate must include the latter as well.
old_signature <- data.frame(
  parameter = c("m_conditional_mean[1]", "predictive_second_moment[1]",
                "ordinal_score_cor[1,2]", "standardized_tau[1,1]"),
  rhat = c(1.001, 1.002, 1.18, 1.002),
  ess_bulk = c(1000, 1200, 18, 950),
  ess_tail = c(900, 1000, 28, 850), n_chains = 4L
)
stopifnot(!mixedgp_diagnostic_table_pass(old_signature))
for (field in c("rhat", "ess_bulk", "ess_tail", "n_chains")) {
  broken <- good
  broken[[field]][1L] <- NA_real_
  stopifnot(!mixedgp_diagnostic_table_pass(broken))
}
for (field in c("ess_bulk", "ess_tail")) {
  broken <- good
  broken[[field]][1L] <- 399
  stopifnot(!mixedgp_diagnostic_table_pass(broken))
}
broken <- good
broken$rhat <- 1.011
stopifnot(!mixedgp_diagnostic_table_pass(broken))
stopifnot(!mixedgp_diagnostic_table_pass(good, expected_parameters = c("m", "missing")))
stopifnot(!mixedgp_diagnostic_table_pass(good[FALSE, ]))

bad_draws <- independent
bad_draws$m[[2L]][9L] <- NA_real_
stopifnot(!mixedgp_summarize_diagnostic_series(bad_draws)$target_pass)
constant <- list(stuck = replicate(4L, rep(1, 2000L), simplify = FALSE))
stopifnot(!mixedgp_summarize_diagnostic_series(constant)$target_pass)
shifted <- list(m = lapply(seq_len(4L), function(i) rnorm(2000L, i)))
stopifnot(!mixedgp_summarize_diagnostic_series(shifted)$target_pass)
sticky <- list(m = lapply(seq_len(4L), function(i) {
  as.numeric(stats::arima.sim(list(ar = 0.995), n = 2000L))
}))
stopifnot(!mixedgp_summarize_diagnostic_series(sticky)$target_pass)

study1_fit <- list(
  data = list(x = matrix(0, 3L, 1L), miss_idx = 2:3),
  mcmc = list(samples_by_chain = list(
    sigma2 = lapply(independent$m, exp),
    logtheta = lapply(1:4, function(i) matrix(rnorm(6000L), 2000L, 3L)),
    tau = lapply(1:4, function(i) matrix(rnorm(4000L), 2000L, 2L)),
    u = lapply(1:4, function(i) cbind(0, matrix(rnorm(4000L), 2000L, 2L)))
  ))
)
study1_series <- mixedgp_study1_raw_series(study1_fit)
stopifnot(all(c("u[2]", "u[3]") %in% names(study1_series)),
          !"u[1]" %in% names(study1_series))
stopifnot(mixedgp_diagnostic_table_pass(mixedgp_study1_raw_diagnostics(study1_fit)))
study1_fit$mcmc$samples_by_chain$u[[1L]][100L, 3L] <- NA_real_
stopifnot(!mixedgp_diagnostic_table_pass(mixedgp_study1_raw_diagnostics(study1_fit)))

## Aggregation must reject stale caches carrying only the former pass flag.
fake_result <- list(cell = list(id = "test", n_rep = 1L, calibration_grid = 0L),
                    outputs = list(mcmc_diagnostics = data.frame(mcmc_pass = TRUE)))
stopifnot(!mixedgp_mcmc_gate(fake_result, list(study = "study2"))$pass)
fake_result$outputs$mcmc_diagnostics <- data.frame(
  mcmc_pass = TRUE, target_functional_pass = TRUE,
  invariant_measurement_pass = FALSE, raw_coordinate_pass = TRUE, n_calib = 0L
)
stopifnot(!mixedgp_mcmc_gate(fake_result, list(study = "study2"))$pass)
fake_result$outputs$mcmc_diagnostics$invariant_measurement_pass <- TRUE
stopifnot(mixedgp_mcmc_gate(fake_result, list(study = "study2"))$pass)
fake_result$outputs$mcmc_diagnostics$n_calib <- 2L
fake_result$outputs$mcmc_diagnostics$raw_coordinate_pass <- FALSE
stopifnot(!mixedgp_mcmc_gate(fake_result, list(study = "study2"))$pass)
fake_result$outputs$mcmc_diagnostics$raw_coordinate_pass <- TRUE
fake_result$outputs$mcmc_diagnostics$target_functional_pass <- FALSE
stopifnot(!mixedgp_mcmc_gate(fake_result, list(study = "study2"))$pass)

## Every proxy pair and every cutpoint must be represented, not selected pairs.
make_iid_fit <- function(n_draw = 2000L, q = 4L, d = 2L) {
  A <- lapply(seq_len(4L), function(id) {
    a <- array(rnorm(n_draw * q * d, 0.5, 0.15), c(n_draw, q, d))
    a[, 1L, 2L] <- 0
    a
  })
  tau <- lapply(seq_len(4L), function(id) {
    out <- matrix(rnorm(n_draw * q * 2L, 0, 0.05), n_draw, q * 2L)
    out <- sweep(out, 2L, rep(c(-0.6, 0.6), q), "+")
    colnames(out) <- unlist(lapply(seq_len(q), function(j) paste0("tau[", j, ",", 1:2, "]")))
    out
  })
  list(data = list(q = q, d = d, m_vec = rep(3L, q)),
       mcmc = list(samples_by_chain = list(A = A, tau = tau)))
}
iid_fit <- make_iid_fit()
invariants <- mixedgp_study2_invariant_series(iid_fit)
stopifnot(length(invariants) == choose(4L, 2L) + 4L * 2L)
stopifnot("ordinal_score_cor[3,4]" %in% names(invariants))
stopifnot(mixedgp_diagnostic_table_pass(mixedgp_summarize_diagnostic_series(invariants)))

## A deterministic posterior state stays constant under the diagnostic
## integration, and postprocessing restores the caller's random-number stream.
tiny <- make_iid_fit(n_draw = 8L, q = 2L, d = 2L)
tiny$data <- c(tiny$data, list(
  p = 1L, X = matrix(c(-1, 0, 1), ncol = 1L), y = c(-0.5, 0.1, 0.4),
  X_center = 0, X_scale = 1, y_center = 0, y_scale = 1,
  U_center = c(0, 0), U_scale = c(1, 1), calib_idx = 1L,
  miss_idx = 2:3, ident = "lower_triangular"
))
tiny$kernel <- list(name = "se", matern_nu = 2.5)
for (id in 1:4) {
  tiny$mcmc$samples_by_chain$A[[id]][] <- rep(c(1, 0.4, 0, 0.9), each = 8L)
  tiny$mcmc$samples_by_chain$tau[[id]][] <- rep(c(-0.6, 0.6, -0.6, 0.6), each = 8L)
}
tiny$mcmc$samples_by_chain$U <- replicate(4L,
  array(rep(c(0, -0.5, 0.7, 0, 0.4, -0.3), each = 8L), c(8L, 3L, 2L)),
  simplify = FALSE)
tiny$mcmc$samples_by_chain$logtheta <- replicate(4L,
  matrix(rep(log(c(2, 0.5, 0.5, 0.5)), each = 8L), 8L, 4L), simplify = FALSE)
tiny$mcmc$samples_by_chain$sigma2 <- replicate(4L, rep(0.1, 8L), simplify = FALSE)
set.seed(4321)
rng_before <- .Random.seed
target_series <- mixedgp_study2_target_series(
  tiny, X = matrix(c(-0.25, 0.25), ncol = 1L), C = matrix(c(1, 2, 2, 3), 2L, 2L),
  U = matrix(c(-0.1, 0.1, 0.2, -0.2), 2L, 2L),
  max_draws_per_chain = 8L, n_latent = 8L, seed = 900L
)
stopifnot(identical(rng_before, .Random.seed))
stopifnot(all(vapply(target_series, function(chains) length(unique(unlist(chains))) == 1L,
                     logical(1L))))
stopifnot(!any(mixedgp_summarize_diagnostic_series(target_series)$target_pass))
cat("Strict diagnostic gate regression tests passed.\n")
