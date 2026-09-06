############################################################
## Validation of the m(x,c) linear-functional implementation
############################################################

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L])))
} else {
  getwd()
}
source(file.path(script_dir, "00_parallel_utils.R"))
source(file.path(script_dir, "00_study1_functions.R"))
source(file.path(script_dir, "00_study2_functions.R"))

if (!exists("MEAN_VALIDATION_OUT_DIR")) {
  MEAN_VALIDATION_OUT_DIR <- file.path(dirname(script_dir), "tables")
}
dir.create(MEAN_VALIDATION_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

############################################################
## 1. Empirical-Q formulas versus brute-force GP aggregation
############################################################

X_train <- cbind(x1 = c(-1, 0, 1))
u_train <- c(-0.8, 0.1, 1.2)
y_train <- c(-0.4, 0.2, 0.7)
X_star <- cbind(x1 = c(-0.3, 0.5))
u_mc <- matrix(c(-1, -0.5, 0, 0.4, 0.2, 0.7, 1.1, 1.5), 4, 2)
M <- nrow(u_mc)
W <- matrix(0, 2, 2 * M)
W[1, seq_len(M)] <- 1 / M
W[2, M + seq_len(M)] <- 1 / M
X_expanded <- X_star[rep(seq_len(nrow(X_star)), each = M), , drop = FALSE]

formula_checks <- lapply(c("se", "matern"), function(kernel_name) {
  logtheta_1d <- log(c(1.1, 0.7, 0.9))
  sigma2_eps <- 0.25
  integrated_1d <- gp_integrated_mean_state_1d(
    X_train, u_train, y_train, X_star, u_mc,
    logtheta_1d, sigma2_eps,
    kernel = kernel_name, matern_nu = 2.5, return_cov = TRUE
  )
  expanded_1d <- gp_predict_draw(
    X_train, u_train, y_train,
    X_expanded, as.vector(u_mc),
    logtheta_1d, sigma2_eps,
    noisy = FALSE, kernel = kernel_name, matern_nu = 2.5,
    return_cov = TRUE
  )

  U_train <- cbind(u1 = u_train, u2 = c(0.4, -0.2, 0.8))
  U_mc <- array(NA_real_, c(M, 2, 2))
  U_mc[, , 1] <- u_mc
  U_mc[, , 2] <- matrix(
    c(-0.7, -0.2, 0.3, 0.8, -1, 0, 0.5, 1), M, 2
  )
  U_expanded <- cbind(as.vector(U_mc[, , 1]), as.vector(U_mc[, , 2]))
  logtheta_2d <- log(c(1.1, 0.7, 0.9, 0.6))
  integrated_2d <- gp_integrated_mean_state_general(
    X_train, U_train, y_train, X_star, U_mc,
    logtheta_2d, sigma2_eps,
    kernel = kernel_name, matern_nu = 2.5, return_cov = TRUE
  )
  expanded_2d <- gp_predict_draw_general(
    X_train, U_train, y_train,
    X_expanded, U_expanded,
    logtheta_2d, sigma2_eps,
    noisy = FALSE, kernel = kernel_name, matern_nu = 2.5,
    return_cov = TRUE
  )

  data.frame(
    kernel = kernel_name,
    univariate_mean_error = max(abs(
      integrated_1d$mean - as.numeric(W %*% expanded_1d$mean)
    )),
    univariate_cov_error = max(abs(
      integrated_1d$cov - W %*% expanded_1d$cov %*% t(W)
    )),
    multivariate_mean_error = max(abs(
      integrated_2d$mean - as.numeric(W %*% expanded_2d$mean)
    )),
    multivariate_cov_error = max(abs(
      integrated_2d$cov - W %*% expanded_2d$cov %*% t(W)
    ))
  )
})
formula_checks <- do.call(rbind, formula_checks)
if (max(as.matrix(formula_checks[, -1L])) > 1e-9) {
  stop("The empirical-Q GP linear-functional formula failed brute-force validation.")
}

############################################################
## 2. Study I finite-M integration against analytic m0
############################################################

set.seed(20260902L)
tau_1d <- c(-5 / 3, -10 / 9, 0, 10 / 9, 5 / 3)
grid_1d <- expand.grid(x = c(-1, 0, 1), c = c(1L, 3L, 5L, 6L))
m_truth_1d <- m0_1d(grid_1d$x, grid_1d$c, tau_1d, "active")
M_grid_1d <- c(32L, 64L, 96L, 256L, 1024L)
M_max_1d <- max(M_grid_1d)
u_pool_1d <- lapply(seq_len(nrow(grid_1d)), function(i) {
  lower <- c(-Inf, tau_1d)[grid_1d$c[i]]
  upper <- c(tau_1d, Inf)[grid_1d$c[i]]
  rtruncnorm_vec(
    rep(0, M_max_1d), rep(1, M_max_1d),
    rep(lower, M_max_1d), rep(upper, M_max_1d)
  )
})
study1_sensitivity <- do.call(rbind, lapply(M_grid_1d, function(M_now) {
  estimate <- vapply(seq_len(nrow(grid_1d)), function(i) {
    mean(f0_1d(
      rep(grid_1d$x[i], M_now),
      u_pool_1d[[i]][seq_len(M_now)],
      scenario = "active"
    ))
  }, numeric(1))
  data.frame(
    n_latent = M_now,
    RMSE = sqrt(mean((estimate - m_truth_1d)^2)),
    max_abs_error = max(abs(estimate - m_truth_1d))
  )
}))
if (tail(study1_sensitivity$RMSE, 1L) > 0.025) {
  stop("Study I high-M integration did not approach analytic m0 closely enough.")
}

############################################################
## 3. Study II nested exact-rejection sensitivity
############################################################

set.seed(20260903L)
true_params <- make_study2_true_params("primary", q = 4L, m = 4L)
pilot_pool <- make_oracle_pool_2d(true_params, n_pool = 100000L)
common_keys <- names(sort(table(pilot_pool$key), decreasing = TRUE))[seq_len(3L)]
C_patterns <- do.call(rbind, strsplit(common_keys, "_", fixed = TRUE))
storage.mode(C_patterns) <- "integer"
X_patterns <- matrix(c(-0.5, 0, 0.5, 0, 0, 0), 3, 2)
M_grid_2d <- c(32L, 64L, 96L, 256L, 1024L)
M_max_2d <- max(M_grid_2d)
U_nested <- lapply(seq_len(nrow(C_patterns)), function(i) {
  sample_oracle_u_rejection(
    C_patterns[i, ], true_params, n_draw = M_max_2d,
    max_batches = 200L, batch_size = 50000L
  )
})
reference_2d <- vapply(seq_len(nrow(C_patterns)), function(i) {
  mean(f0_2d(
    X_patterns[rep(i, M_max_2d), , drop = FALSE],
    U_nested[[i]], scenario = "primary"
  ))
}, numeric(1))
study2_sensitivity <- do.call(rbind, lapply(M_grid_2d, function(M_now) {
  estimate <- vapply(seq_len(nrow(C_patterns)), function(i) {
    mean(f0_2d(
      X_patterns[rep(i, M_now), , drop = FALSE],
      U_nested[[i]][seq_len(M_now), , drop = FALSE],
      scenario = "primary"
    ))
  }, numeric(1))
  data.frame(
    n_latent = M_now,
    RMSE_to_1024 = sqrt(mean((estimate - reference_2d)^2)),
    max_abs_difference_to_1024 = max(abs(estimate - reference_2d))
  )
}))

write.csv(
  formula_checks,
  file.path(MEAN_VALIDATION_OUT_DIR, "mean_function_formula_validation.csv"),
  row.names = FALSE
)
write.csv(
  study1_sensitivity,
  file.path(MEAN_VALIDATION_OUT_DIR, "study1_mean_integration_sensitivity.csv"),
  row.names = FALSE
)
write.csv(
  study2_sensitivity,
  file.path(MEAN_VALIDATION_OUT_DIR, "study2_mean_integration_sensitivity.csv"),
  row.names = FALSE
)

cat("m(x,c) formula and finite-integration validation passed.\n")
print(formula_checks)
print(study1_sensitivity)
print(study2_sensitivity)
