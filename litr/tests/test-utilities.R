## Reusable regression tests retained from the literate package source.
## These exercise current numerical and input contracts without old-fit fixtures.

testthat::test_that("parallel maps are reproducible", {
  f <- function(i) c(i = i, u = stats::runif(1))
  serial <- mixedgp_parallel_lapply(
    as.list(1:3), f, n_cores = 1, seeds = 101:103
  )
  parallel <- mixedgp_parallel_lapply(
    as.list(1:3), f, n_cores = 2, seeds = 101:103
  )
  testthat::expect_identical(serial, parallel)

  serial_m <- mixedgp_parallel_mapply(
    function(a, b) c(sum = a + b, u = stats::runif(1)),
    a = 1:3, b = 3:1, n_cores = 1, seeds = 201:203
  )
  parallel_m <- mixedgp_parallel_mapply(
    function(a, b) c(sum = a + b, u = stats::runif(1)),
    a = 1:3, b = 3:1, n_cores = 2, seeds = 201:203
  )
  testthat::expect_identical(serial_m, parallel_m)
})

testthat::test_that("seeded parallel maps preserve the caller RNG stream", {
  set.seed(4101)
  expected_lapply_continuation <- stats::runif(4)
  set.seed(4101)
  invisible(mixedgp_parallel_lapply(
    as.list(1:3), function(i) stats::runif(i),
    n_cores = 2, seeds = 501:503
  ))
  testthat::expect_identical(stats::runif(4), expected_lapply_continuation)

  set.seed(4102)
  expected_mapply_continuation <- stats::runif(4)
  set.seed(4102)
  invisible(mixedgp_parallel_mapply(
    function(a, b) c(a + b, stats::runif(1)),
    a = 1:3, b = 3:1, n_cores = 2, seeds = 601:603
  ))
  testthat::expect_identical(stats::runif(4), expected_mapply_continuation)
})

testthat::test_that("tail-stable truncated normals respect extreme support", {
  lower <- rep(c(8, -Inf, 40, -41), each = 64L)
  upper <- rep(c(Inf, -8, 41, -40), each = 64L)
  set.seed(4103)
  z <- rtruncnorm_vec(0, 1, lower, upper)
  testthat::expect_true(all(is.finite(z)))
  testthat::expect_true(all(z >= lower & z <= upper))
})

testthat::test_that("ordinal factor maps survive dropped prediction levels", {
  training <- data.frame(severity = ordered(
    c("low", "high", "medium"),
    levels = c("low", "medium", "high")
  ))
  encoded <- prepare_ordinal_matrix(training, name = "C")
  prediction <- data.frame(severity = ordered(
    c("high", "low"), levels = c("low", "high")
  ))
  encoded_new <- prepare_ordinal_matrix(
    prediction,
    m_vec = encoded$m_vec,
    level_maps = encoded$level_maps,
    expected_names = encoded$column_names,
    name = "new_C"
  )
  testthat::expect_identical(unname(encoded_new$C[, 1L]), c(3L, 1L))
  testthat::expect_error(
    prepare_ordinal_matrix(
      data.frame(severity = ordered("unseen")),
      m_vec = encoded$m_vec,
      level_maps = encoded$level_maps,
      expected_names = encoded$column_names,
      name = "new_C"
    ),
    "unknown level"
  )
})

testthat::test_that("multivariate latent dimension is never silently guessed", {
  X <- matrix(c(-1, 0, 1), ncol = 1L)
  y <- c(-0.5, 0, 0.5)
  C <- data.frame(
    proxy_a = ordered(c("low", "medium", "high")),
    proxy_b = ordered(c("none", "mild", "severe"))
  )
  testthat::expect_error(
    fit_eivgp(X, y, C, engine = "multivariate", U_obs = NULL),
    "latent_dim must be supplied explicitly"
  )
  testthat::expect_error(
    fit_eivgp(
      X, y, C, engine = "multivariate",
      U_obs = matrix(NA_real_, nrow = 3L, ncol = 2L)
    ),
    "U_obs contains no calibration values"
  )

  rank_deficient <- mixedgp_latent_anchor_status(
    matrix(c(0, 0, 1, 1), ncol = 2L, byrow = TRUE), 1:2, d = 2L
  )
  anchored <- mixedgp_latent_anchor_status(
    matrix(c(0, 0, 1, 0, 0, 1), ncol = 2L, byrow = TRUE), 1:3, d = 2L
  )
  testthat::expect_false(rank_deficient$anchored)
  testthat::expect_identical(rank_deficient$affine_rank, 2L)
  testthat::expect_true(anchored$anchored)
  testthat::expect_identical(anchored$affine_rank, 3L)
})

testthat::test_that("duplicate-point GP predictions remain positive semidefinite", {
  prediction <- gp_predict_draw_general(
    X_train = matrix(c(-1, 0, 1), ncol = 1L),
    U_train = matrix(c(-0.5, 0, 0.5), ncol = 1L),
    y_train = c(-0.7, 0, 0.8),
    X_star = matrix(c(0.25, 0.25), ncol = 1L),
    U_star = matrix(c(0.1, 0.1), ncol = 1L),
    logtheta = log(c(1, 1, 1)),
    sigma2_eps = 0.1,
    return_cov = TRUE
  )
  testthat::expect_equal(prediction$mean[1L], prediction$mean[2L])
  testthat::expect_true(
    min(eigen(prediction$cov, symmetric = TRUE, only.values = TRUE)$values) >=
      -1e-10
  )
  set.seed(4104)
  draws <- rmvnorm_psd(32L, prediction$mean, prediction$cov)
  testthat::expect_true(all(is.finite(draws)))
  testthat::expect_equal(draws[, 1L], draws[, 2L], tolerance = 1e-7)
})

