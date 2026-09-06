############################################################
## 13_study2_prospective_latent_validation.R
##
## Focused release-gate checks for prospective U | C draws.
## The publication path must use an exact draw under U ~ N(0, I).
## It samples the truncated Gaussian score by minimax tilting and then
## samples U | S exactly. Prior rejection is an independent validation
## path; finite Gibbs is retained only as a diagnostic.
############################################################

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
codes_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
} else {
  getwd()
}

source(file.path(codes_dir, "00_project_setup.R"), chdir = TRUE)
source(file.path(codes_dir, "00_parallel_utils.R"))
source(file.path(codes_dir, "00_study2_functions.R"))
source(file.path(codes_dir, "04_study2_ablations.R"))

settings_check <- lapply(c("quick", "balanced", "thorough"), study2_config_settings)
stopifnot(
  all(vapply(settings_check, function(x) {
    identical(x$predictive_latent_sampler, "minimax_tilting")
  }, logical(1))),
  all(vapply(settings_check, function(x) {
    identical(x$diagnostic_n_new_latent_gibbs, 100L)
  }, logical(1)))
)

## Tail-stability and truncation-support checks used by the Gibbs diagnostic.
set.seed(17001L)
upper_tail <- rtruncnorm_vec(rep(10, 2000L), 1, -Inf, 0)
lower_tail <- rtruncnorm_vec(rep(-10, 2000L), 1, 0, Inf)
stopifnot(
  all(is.finite(upper_tail)), all(upper_tail <= 0),
  all(is.finite(lower_tail)), all(lower_tail >= 0)
)

## A one-dimensional case with a directly integrated reference distribution.
## Four concordant upper categories make short fresh-start Gibbs visibly biased.
A <- matrix(3, nrow = 4L, ncol = 1L)
m_vec <- rep(3L, 4L)
tau <- replicate(4L, c(-0.5, 0.5), simplify = FALSE)
n_draw <- 20000L
C_upper <- matrix(3L, nrow = n_draw, ncol = 4L)

grid <- seq(-6, 6, length.out = 120001L)
du <- grid[[2L]] - grid[[1L]]
unnormalized <- dnorm(grid) * pnorm(3 * grid - 0.5)^4
normalizer <- sum(unnormalized) * du
reference_mean <- sum(grid * unnormalized) * du / normalizer
reference_var <- sum((grid - reference_mean)^2 * unnormalized) * du / normalizer
reference_mu4 <- sum((grid - reference_mean)^4 * unnormalized) * du / normalizer

set.seed(17002L)
minimax_draw <- sample_u_given_c_ordprobit_minimax(
  C_new = C_upper,
  A = A,
  tau = tau,
  m_vec = m_vec
)
minimax_mean <- mean(minimax_draw)
minimax_var <- var(as.numeric(minimax_draw))
mean_tolerance <- 6 * sqrt(reference_var / n_draw) + 0.002
var_tolerance <- 6 * sqrt((reference_mu4 - reference_var^2) / n_draw) + 0.004
stopifnot(
  identical(attr(minimax_draw, "latent_sampler"),
            "exact_minimax_tilting"),
  abs(minimax_mean - reference_mean) <= mean_tolerance,
  abs(minimax_var - reference_var) <= var_tolerance
)

## Compare the publication sampler with independent exact prior rejection.
set.seed(17003L)
rejection_draw <- sample_u_given_c_ordprobit_rejection(
  C_new = C_upper,
  A = A,
  tau = tau,
  m_vec = m_vec,
  batch_size = 100000L,
  max_batches = 50L
)
rejection_mean <- mean(rejection_draw)
rejection_var <- var(as.numeric(rejection_draw))
rejection_telemetry <- attr(rejection_draw, "pattern_telemetry")
stopifnot(
  identical(attr(rejection_draw, "latent_sampler"), "exact_rejection"),
  is.finite(attr(rejection_draw, "overall_acceptance")),
  is.data.frame(rejection_telemetry),
  nrow(rejection_telemetry) == 1L,
  rejection_telemetry$requested_draws == n_draw,
  rejection_telemetry$retained_draws == n_draw,
  rejection_telemetry$proposal_hits >= rejection_telemetry$retained_draws,
  rejection_telemetry$total_candidates ==
    attr(rejection_draw, "candidate_draws"),
  rejection_telemetry$complete,
  abs(rejection_mean - reference_mean) <= mean_tolerance,
  abs(rejection_var - reference_var) <= var_tolerance
)

## A frozen posterior state from the development audit contains an ordinal
## pattern whose probability is about 6.21e-8. Prior rejection has about a
## 73% chance of returning no hit in five million proposals. Independent
## two-dimensional quadrature supplies the reference moments below.
A_rare <- structure(c(
  0.0333061757923723, -2.72318504957344, -4.78716379467788,
  -4.14661836467995, 0, 1.75291234472192,
  -2.77181281706168, -0.905880468957085
), dim = c(4L, 2L))
tau_rare <- list(
  c(-0.87147255687319, 0.0754277803401278, 2.73715395786259),
  c(-1.99550503041999, 0.507117989058624, 1.65433442066491),
  c(-4.55704694724826, -2.47351849275641, 2.99549722775762),
  c(-3.75172907520982, -1.11104566026423, -0.629214632472905)
)
n_rare <- 20000L
C_rare <- matrix(rep(c(4L, 3L, 4L, 2L), each = n_rare),
                 nrow = n_rare)
rare_mean_ref <- c(-0.185109899540878, -0.362959439102707)
rare_cov_ref <- matrix(c(
  0.0307388638448463, -0.0272495852840934,
  -0.0272495852840934, 0.119875778131897
), 2L, 2L, byrow = TRUE)
set.seed(17004L)
rare_draw <- sample_u_given_c_ordprobit_minimax(
  C_new = C_rare,
  A = A_rare,
  tau = tau_rare,
  m_vec = rep(4L, 4L)
)
stopifnot(
  identical(dim(rare_draw), c(n_rare, 2L)),
  all(is.finite(rare_draw)),
  identical(attr(rare_draw, "latent_sampler"), "exact_minimax_tilting")
)
rare_mean <- colMeans(rare_draw)
rare_centered <- sweep(rare_draw, 2L, rare_mean_ref, "-")
rare_mean_se <- sqrt(diag(rare_cov_ref) / n_rare)
stopifnot(all(abs(rare_mean - rare_mean_ref) <= 6 * rare_mean_se + 1e-4))
for (j in 1:2) {
  for (k in j:2) {
    cross_moment <- rare_centered[, j] * rare_centered[, k]
    cross_se <- stats::sd(cross_moment) / sqrt(n_rare)
    stopifnot(
      abs(mean(cross_moment) - rare_cov_ref[j, k]) <= 6 * cross_se + 1e-4
    )
  }
}

## Underfilled rejection proposals must fail, not return a partial/biased draw.
underfilled <- tryCatch(
  sample_u_given_c_ordprobit_rejection(
    C_new = matrix(3L, nrow = 2L, ncol = 4L),
    A = A,
    tau = tau,
    m_vec = m_vec,
    batch_size = 1L,
    max_batches = 1L
  ),
  error = function(e) e
)
stopifnot(
  inherits(underfilled, "mixedgp_rejection_sampling_error"),
  is.data.frame(underfilled$pattern_telemetry),
  any(underfilled$pattern_telemetry$shortfall > 0L),
  any(!underfilled$pattern_telemetry$complete),
  underfilled$candidate_draws == 1
)

## EIV-GP, PI-GP, and CC-GP all default to minimax-tilted accept--reject
## sampling, which is exact on successful completion.
stopifnot(
  identical(eval(formals(measurement_posterior_mean_u_test)$latent_sampler)[1L],
            "minimax_tilting"),
  identical(eval(formals(sample_study2_pi_gp)$latent_sampler)[1L],
            "minimax_tilting"),
  identical(eval(formals(sample_study2_cc_gp)$latent_sampler)[1L],
            "minimax_tilting")
)

## Exercise both explicitly labeled branches of the measurement-model helper.
samples_A <- array(NA_real_, dim = c(2L, 4L, 1L))
samples_A[1L, , 1L] <- A[, 1L]
samples_A[2L, , 1L] <- A[, 1L]
samples_tau <- rbind(flatten_tau(tau), flatten_tau(tau))
measurement_fit <- list(
  data = list(
    m_vec = m_vec,
    q = 4L,
    d = 1L,
    U_obs = matrix(c(-1, 1, NA, NA), ncol = 1L),
    calib_idx = 1:2,
    miss_idx = 3:4,
    latent_scale_anchored = TRUE,
    latent_anchor_rank = 2L,
    latent_anchor_required_rank = 2L
  ),
  mcmc = list(
    samples_A = samples_A,
    samples_tau = samples_tau,
    samples_U = array(
      c(-1, -1, 1, 1, -0.2, -0.1, 0.4, 0.5),
      dim = c(2L, 4L, 1L)
    )
  )
)
C_small <- matrix(3L, nrow = 12L, ncol = 4L)

set.seed(17005L)
pi_exact <- measurement_posterior_mean_u_test(
  measurement_fit, C_small, max_draw = 2L
)
stopifnot(
  identical(attr(pi_exact, "latent_sampler"), "exact_minimax_tilting"),
  all(is.finite(pi_exact))
)

set.seed(17006L)
pi_gibbs <- measurement_posterior_mean_u_test(
  measurement_fit,
  C_small,
  max_draw = 2L,
  n_gibbs = 2L,
  latent_sampler = "gibbs"
)
stopifnot(
  identical(attr(pi_gibbs, "latent_sampler"), "finite_gibbs_diagnostic"),
  identical(attr(pi_gibbs, "gibbs_sweeps"), 2L),
  all(is.finite(pi_gibbs))
)

## The response-free imputation comparator retains full posterior draws rather
## than only the conditional mean used by PI-GP. Raw-coordinate scoring is
## available only with a full-affine-rank calibration design.
measurement_draws <- sample_measurement_u_given_c(
  measurement_fit,
  C_small,
  max_draw = 2L,
  scale = "raw",
  latent_sampler = "minimax_tilting",
  seed = 17007L
)
stopifnot(
  identical(dim(measurement_draws), c(2L, nrow(C_small), 1L)),
  identical(attr(measurement_draws, "scale"), "raw"),
  isTRUE(attr(measurement_draws, "latent_scale_anchored")),
  all(is.finite(measurement_draws))
)
training_score <- study2_measurement_training_imputation_metrics(
  measurement_fit = measurement_fit,
  U_true = matrix(c(-1, 1, -0.15, 0.45), ncol = 1L),
  rep_id = 1L,
  n_calib = 2L,
  scenario = "validation"
)
stopifnot(
  nrow(training_score) == 1L,
  training_score$method == "Ordinal model (no Y)",
  training_score$target == "training_missing_U",
  all(is.finite(unlist(training_score[c(
    "Bias", "RMSE", "MAE", "Coverage95", "Width95"
  )])))
)

measurement_fit_unanchored <- measurement_fit
measurement_fit_unanchored$data$U_obs[,] <- NA_real_
measurement_fit_unanchored$data$calib_idx <- integer(0)
measurement_fit_unanchored$data$miss_idx <- seq_len(4L)
measurement_fit_unanchored$data$latent_scale_anchored <- FALSE
measurement_fit_unanchored$data$latent_anchor_rank <- 0L
unanchored_status <- study2_latent_imputation_status(
  fit = measurement_fit_unanchored,
  method = "Ordinal model (no Y)",
  target = "prospective_U_given_C",
  rep_id = 1L,
  n_calib = 0L,
  scenario = "validation",
  n_units = nrow(C_small)
)
unanchored_raw <- tryCatch(
  sample_measurement_u_given_c(
    measurement_fit_unanchored,
    C_small,
    max_draw = 2L,
    scale = "raw",
    seed = 17008L
  ),
  error = function(e) e
)
stopifnot(
  identical(unanchored_status$status, "not_identified"),
  !unanchored_status$task_eligible,
  inherits(unanchored_raw, "error")
)

cat(
  sprintf(
    paste0(
      "Prospective-latent validation passed: reference mean %.5f, ",
      "minimax/rejection means %.5f/%.5f; reference variance %.5f, ",
      "minimax/rejection variances %.5f/%.5f; ",
      "overall acceptance %.5g.\n"
    ),
    reference_mean,
    minimax_mean,
    rejection_mean,
    reference_var,
    minimax_var,
    rejection_var,
    attr(rejection_draw, "overall_acceptance")
  )
)
