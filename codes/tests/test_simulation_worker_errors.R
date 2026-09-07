## Run from repository root in a clean R session; never loads an installed eivGP.
source("codes/simulation_helpers.R")
engine <- mixedgp_simulation_engine(normalizePath("codes"))
stopifnot(exists("mixedgp_validate_named_dots", engine, inherits = FALSE))
for (measurement in c("univariate", "multivariate")) {
  u <- seq(-1, 1, length.out = 10L)
  C <- matrix(as.integer(u > 0) + 1L, ncol = 1L)
  if (measurement == "multivariate") C <- cbind(C, rep(1:2, 5L))
  U <- matrix(u, ncol = 1L); U[3:8, ] <- NA_real_
  args <- list(X = matrix(seq(-1, 1, length.out = 10L), ncol = 1L),
    y = sin(u), C = C, U_obs = U, engine = measurement,
    n_iter = 8L, burn = 4L, n_chains = 2L, parallel = FALSE, seed = 913L)
  if (measurement == "multivariate") args$ident <- "none"
  fit <- do.call(engine$fit_eivgp, args)
  stopifnot(identical(fit$sampler_version, "0.3.0"))
}
## Helper needs the same parallel function that cell drivers inherit.
for (cores in c(1L, 2L)) {
  path <- tempfile(fileext = ".csv")
  result <- suppressWarnings(mixedgp_run_replications(as.list(1:3),
    function(i) {if (i == 2L) stop("original worker failure"); list(metrics = i)},
    n_cores = cores, seeds = 11:13, status_path = path, study = "test",
    parallel_map = engine$mixedgp_parallel_lapply))
  status <- read.csv(path)
  stopifnot(identical(vapply(result, `[[`, 0L, "metrics"), c(1L, 3L)),
    identical(status$status, c("success", "failed", "success")),
    identical(status$message[2L], "original worker failure"))
  failure <- tryCatch(mixedgp_run_replications(as.list(1:2),
    function(i) stop("original all-failed error"), n_cores = cores,
    seeds = 11:12, status_path = path, study = "test",
    parallel_map = engine$mixedgp_parallel_lapply), error = identity)
  stopifnot(inherits(failure, "error"),
    grepl("original all-failed error", conditionMessage(failure)),
    !grepl("subscript out of bounds", conditionMessage(failure)),
    all(read.csv(path)$status == "failed"))
}
message("Clean-loader fits for both engines and serial/parallel failure reporting passed.")
