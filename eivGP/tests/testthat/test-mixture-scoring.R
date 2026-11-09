test_that("blocked mixture CRPS agrees with the full pairwise formula", {
  api <- eivgp_numerical_api(system.file(package = "eivGP"))
  set.seed(1907)
  for (n in c(1L, 2L, 63L, 64L, 65L, 129L, 257L)) {
    mu <- rnorm(n); v <- exp(rnorm(n)); y <- .7
    expected <- mean(api$normal_abs(y - mu, sqrt(v))) - .5 *
      mean(api$normal_abs(outer(mu, mu, "-"), sqrt(outer(v, v, "+"))))
    actual <- api$mixture_score(y, mu, v)
    expect_equal(unname(actual["crps"]), expected, tolerance = 1e-12)
    expect_true(all(is.finite(actual)))
    expect_gte(unname(actual["crps"]), 0)
    expect_equal(api$mixture_score(y, rep(mu, 2), rep(v, 2)), actual,
      tolerance = 1e-12)
  }
  expect_equal(unname(api$mixture_score(0, 0, 1)["crps"]),
    (sqrt(2) - 1) / sqrt(pi), tolerance = 1e-12)
  for (args in list(list(0, numeric(), numeric()), list(0, 1:2, 1),
      list(1:2, 0, 1), list(0, Inf, 1), list(0, 0, -1)))
    expect_error(do.call(api$mixture_score, args), "Invalid predictive components")
})

test_that("prediction scoring is serial within a replicate", {
  api <- eivgp_numerical_api(system.file(package = "eivGP"))
  api$mixedgp_stage_lapply <- function(...) stop("Unexpected scoring fork")
  dat <- list(train = list(C = matrix(1:2)),
    test = list(C = matrix(1:3), y = c(-1, 0, 1)))
  pred <- list(y_mean = matrix(0, 3, 2), y_var = matrix(1, 3, 2),
    f_mean = matrix(0, 3, 2), f_var = matrix(1, 3, 2),
    m_mean = matrix(0, 3, 2), m_var = matrix(1, 3, 2))
  truth <- list(m = c(0, .1, .2), f = c(-.1, 0, .1))
  set.seed(2); before <- .Random.seed
  scores <- api$score_result(pred, dat, 1L, truth, "test")
  expect_identical(.Random.seed, before)
  expect_equal(nrow(scores), 9L)
  expect_setequal(scores$target, c("y", "m", "f"))
  expect_equal(scores$stratum[scores$target == "y"], c("seen", "seen", "unseen"))
  expect_true(all(is.finite(scores$coverage)))
})
