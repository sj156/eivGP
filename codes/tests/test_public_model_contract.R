# Public-model contract tests; run with Rscript from any working directory.
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
stopifnot(length(script) == 1L)
codes <- dirname(dirname(normalizePath(script)))
for (module in c("00_parallel_utils.R", "00_study1_functions.R",
                 "00_study2_functions.R", "00_public_api.R")) {
  source(file.path(codes, module))
}
expect_error_text <- function(expr, pattern) {
  result <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(result, "error"), grepl(pattern, conditionMessage(result)))
}
stopifnot(
  identical(mixedgp_resolve_loading_rule(NULL, 0L, FALSE), "lower_triangular"),
  identical(mixedgp_resolve_loading_rule(NULL, 4L, TRUE), "none"),
  identical(mixedgp_resolve_loading_rule("none", 1L, FALSE), "none"),
  identical(mixedgp_resolve_loading_rule("lower_triangular", 1L, FALSE),
            "lower_triangular")
)
expect_error_text(mixedgp_resolve_loading_rule(NULL, 1L, FALSE),
                  "explicit ident choice")
expect_error_text(mixedgp_resolve_loading_rule("none", 0L, FALSE),
                  "supported-workflow restriction")

set.seed(912L)
n <- 12L
X <- matrix(seq(-1, 1, length.out = n), ncol = 1L)
U <- matrix(rnorm(2L * n), n, 2L)
A <- rbind(c(1, 0.4), c(-0.3, 1), c(0.7, -0.4))
S <- U %*% t(A) + matrix(rnorm(3L * n), n, 3L)
C <- matrix(as.integer(S > 0) + 1L, n, 3L)
y <- sin(X[, 1L]) + 0.3 * U[, 1L] + rnorm(n, sd = 0.15)
U_obs <- U
U_obs[-c(1L, 2L), ] <- NA_real_
expect_error_text(
  fit_eivgp(X, y, C, U_obs, engine = "multivariate", m_vec = rep(2L, 3L)),
  "explicit ident choice"
)
caller_rng <- .Random.seed
fit <- suppressWarnings(fit_eivgp(
  X, y, C, U_obs, engine = "multivariate", m_vec = rep(2L, 3L),
  ident = "none", standardize_U = FALSE, parallel = FALSE,
  n_chains = 2L, n_iter = 60L, burn = 20L, thin = 1L, seed = 14L
))
stopifnot(
  identical(.Random.seed, caller_rng),
  identical(fit$data$ident, "none"),
  !isTRUE(fit$data$latent_scale_anchored),
  identical(fit$model_specification$prospective_convention, "fixed_fit"),
  identical(fit$model_specification$loading_structure, "none"),
  !fit$model_specification$standardize_U,
  identical(summary(fit)$model_specification, fit$model_specification)
)
printed <- capture.output(print(summary(fit)))
stopifnot(any(grepl("Latent population:", printed)),
          any(grepl("not a necessary-condition theorem", printed)))
message("Public model-contract tests passed; no convergence claim from this short fit.")
