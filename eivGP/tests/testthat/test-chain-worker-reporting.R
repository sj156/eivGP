test_that("fits report the effective shared chain-worker ceiling", {
  old <- options(eivGP.checkpoint_root = NULL,
                 eivGP.checkpoint_seconds = Inf,
                 eivGP.checkpoint_iterations = Inf)
  on.exit(options(old), add = TRUE)
  run <- function(cores, expr) {
    release <- mixedgp_budget_enter(cores)
    on.exit(release(), add = TRUE)
    force(expr)
  }
  u <- seq(-1, 1, length.out = 8L)
  U <- matrix(u, ncol = 1L); U[3:6, ] <- NA_real_
  for (engine in c("univariate", "multivariate")) {
    C <- matrix(as.integer(u > 0) + 1L, ncol = 1L)
    if (engine == "multivariate") C <- cbind(C, rep(1:2, 4L))
    args <- list(X = matrix(u, ncol = 1L), y = sin(u), C = C, U_obs = U,
                 engine = engine, m_vec = rep(2L, ncol(C)),
                 n_iter = 20L, burn = 8L, n_chains = 4L, seed = 736L,
                 verbose = FALSE)
    if (engine == "multivariate") args$ident <- "none"
    for (sampler in c("ESS", "HMC", "HYBRID")) {
      args$sampler <- sampler
      reference <- run(1L, do.call(fit_eivgp, c(args, list(n_cores = 1L))))
      for (budget in c(1L, 2L, 4L)) {
        fit <- run(budget, do.call(fit_eivgp, c(args, list(n_cores = 12L))))
        capacity <- if (.Platform$OS.type == "unix") {
          min(budget, mixedgp_available_cores())
        } else 1L
        expect_identical(fit$diagnostics$summary$n_chains, 4L)
        expect_identical(fit$diagnostics$summary$parallel_cores, capacity)
        expect_identical(fit$diagnostics$summary$parallel_backend,
                         if (capacity > 1L) "fork" else "serial")
        expect_identical(fit$interface$n_cores, capacity)
        expect_identical(fit$mcmc$samples_by_chain, reference$mcmc$samples_by_chain)
        expect_identical(fit$checkpoint$states, reference$checkpoint$states)
      }
      extended <- run(1L, continue_eivgp(fit, n_iter = 4L, n_cores = 12L))
      expect_identical(extended$diagnostics$summary$parallel_cores, 1L)
      expect_identical(extended$diagnostics$summary$parallel_backend, "serial")
      expect_identical(extended$interface$n_cores, 1L)
      expect_identical(extended$diagnostics$summary$n_chains, 4L)
    }
  }
  expect_null(mixedgp_budget())
})

test_that("complete fit groups retain the full budget for later stages", {
  check_pool <- function(cores) {
    release <- mixedgp_budget_enter(cores)
    on.exit(release(), add = TRUE)
    capacity <- if (.Platform$OS.type == "unix") {
      min(as.integer(cores), mixedgp_available_cores())
    } else 1L
    expect_identical(mixedgp_budget_cores(), capacity)
  }
  for (cores in 13:15) {
    allocation <- mixedgp_application_allocation(cores)
    expect_identical(allocation$chains, 4L)
    expect_identical(allocation$jobs, 3L)
    expect_identical(allocation$requested_cores, as.integer(cores))
    expect_identical(mixedgp_stage_workers(cores, "fit", 20L, 4L), 3L)
    check_pool(cores)
  }
  expect_null(mixedgp_budget())
})
