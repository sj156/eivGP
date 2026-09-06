############################################################
## Small end-to-end validation of the real-data interface
############################################################

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  getwd()
}

source(file.path(script_dir, "00_parallel_utils.R"))
source(file.path(script_dir, "00_study2_functions.R"))
source(file.path(script_dir, "real-data.R"))

set.seed(20260830)
n <- 18L
p <- 3L
d <- 2L
q <- 4L

X <- data.frame(
  age = rnorm(n, 65, 8),
  biomarker = rnorm(n, 0, 1.5),
  exposure = runif(n, -2, 2)
)
U_model <- cbind(
  severity = rnorm(n),
  reserve = rnorm(n)
)
U_raw <- cbind(
  severity = 50 + 8 * U_model[, 1],
  reserve = -3 + 2.5 * U_model[, 2]
)
A_true <- rbind(
  c(1.5, 0.0),
  c(0.3, 1.4),
  c(1.0, 0.6),
  c(0.5, 1.1)
)
S <- U_model %*% t(A_true) + matrix(rnorm(n * q), n, q)
level_labels <- list(
  c("none", "moderate", "severe"),
  c("none", "mild", "moderate", "severe"),
  c("none", "trace", "mild", "moderate", "severe"),
  c("none", "mild", "moderate", "severe")
)
cut_points <- lapply(level_labels, function(z) {
  qnorm(seq_len(length(z) - 1L) / length(z))
})
C_code <- vapply(seq_len(q), function(j) {
  as.integer(cut(S[, j], c(-Inf, cut_points[[j]], Inf), labels = FALSE))
}, integer(n))
C <- as.data.frame(lapply(seq_len(q), function(j) {
  ordered(level_labels[[j]][C_code[, j]], levels = level_labels[[j]])
}))
names(C) <- c("proxy_a", "proxy_b", "proxy_c", "proxy_d")
y <-
  0.03 * X$age +
  0.45 * X$biomarker -
  0.2 * X$exposure +
  sin(U_model[, 1]) +
  0.5 * U_model[, 2]^2 +
  rnorm(n, 0, 0.15)

calib_idx <- c(2L, 4L, 7L, 10L, 14L, 17L)
U_obs <- matrix(NA_real_, n, d, dimnames = list(NULL, colnames(U_raw)))
U_obs[calib_idx, ] <- U_raw[calib_idx, ]

## Kernel matrices must be symmetric, have unit diagonal, and be positive
## semidefinite for both supported covariance families.
X_std <- scale(as.matrix(X))
U_std <- scale(U_raw)
logtheta <- log(c(2, rep(0.7, p + d)))
for (kernel_name in c("se", "matern")) {
  R <- gp_corr_general(
    X_std,
    U_std,
    logtheta,
    kernel = kernel_name,
    matern_nu = 2.5
  )$R
  stopifnot(max(abs(R - t(R))) < 1e-10)
  stopifnot(max(abs(diag(R) - 1)) < 1e-10)
  stopifnot(min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) > -1e-8)
}

fits <- list()
for (kernel_name in c("se", "matern")) {
  fits[[kernel_name]] <- fit_eivgp_real_data(
    X = X,
    y = y,
    C = C,
    U_obs = U_obs,
    kernel = kernel_name,
    matern_nu = 2.5,
    n_iter = 80L,
    burn = 40L,
    thin = 2L,
    n_chains = 1L,
    preset = "fast",
    sampler_strategy = "interwoven",
    control_overrides = list(
      joint_global_every = 8L,
      joint_collapsed_every = 8L,
      marginal_measurement_every = 4L
    ),
    seed = 7300L + match(kernel_name, c("se", "matern")),
    parallel_chains = FALSE,
    verbose = FALSE
  )

  fit <- fits[[kernel_name]]
  stopifnot(fit$data$p == p, fit$data$q == q, fit$data$d == d)
  stopifnot(fit$kernel$name == kernel_name)
  stopifnot(is.null(fit$mcmc$samples_S))
  stopifnot(identical(fit$data$calib_idx, calib_idx))
  stopifnot(all(is.finite(fit$mcmc$samples_logtheta)))
  stopifnot(all(fit$mcmc$samples_sigma2 > 0))

  standardized_calibration <- sweep(
    sweep(U_raw[calib_idx, , drop = FALSE], 2, fit$data$U_center, "-"),
    2,
    fit$data$U_scale,
    "/"
  )
  stopifnot(max(abs(
    fit$mcmc$samples_U[1, calib_idx, ] - standardized_calibration
  )) < 1e-10)

  test_idx <- 1:4
  X_new <- X[test_idx, c("exposure", "age", "biomarker")]
  C_new <- C[test_idx, c("proxy_d", "proxy_b", "proxy_a", "proxy_c")]
  U_new_obs <- matrix(
    NA_real_,
    length(test_idx),
    d,
    dimnames = list(NULL, rev(colnames(U_raw)))
  )
  U_new_obs[1, ] <- U_raw[test_idx[1], rev(colnames(U_raw))]

  ids <- seq_len(min(4L, nrow(fit$mcmc$samples_logtheta)))
  y_draws <- predict_eivgp_y_given_xc(
    fit,
    X_new = X_new,
    C_new = C_new,
    U_new_obs = U_new_obs,
    draw_ids = ids,
    joint = TRUE,
    seed = 8100L
  )
  f_draws <- predict_eivgp_f_given_xu(
    fit,
    X_new = X_new,
    U_new = U_raw[test_idx, , drop = FALSE],
    draw_ids = ids,
    joint = TRUE,
    seed = 8200L
  )
  raw_u_draws <- posterior_u_draws(fit, rows = calib_idx, scale = "raw")

  stopifnot(identical(dim(y_draws), c(length(ids), length(test_idx))))
  stopifnot(identical(dim(f_draws), c(length(ids), length(test_idx))))
  stopifnot(all(is.finite(y_draws)), all(is.finite(f_draws)))
  stopifnot(max(abs(raw_u_draws[1, , ] - U_raw[calib_idx, ])) < 1e-8)
}

## A partially observed latent row must fail loudly rather than be silently
## treated as either calibrated or uncalibrated.
U_bad <- U_obs
U_bad[1, 1] <- U_raw[1, 1]
partial_error <- try(
  fit_eivgp_real_data(
    X = X,
    y = y,
    C = C,
    U_obs = U_bad,
    n_iter = 4L,
    burn = 2L,
    n_chains = 1L
  ),
  silent = TRUE
)
stopifnot(inherits(partial_error, "try-error"))

cat("Real-data interface checks passed for SE and Matern kernels.\n")
