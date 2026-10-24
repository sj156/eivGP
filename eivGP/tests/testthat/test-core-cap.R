test_that("requested budgets are capped without changing chain counts", {
  # Virtual hosts: no workers or experiments are launched.
  env <- new.env(parent=environment(mixedgp_resolve_cores))
  resolve <- mixedgp_resolve_cores; environment(resolve) <- env
  env$mixedgp_available_cores <- function()16L
  for(request in c(1L,2L,4L,8L,16L,17L,56L))
    expect_equal(resolve(request),min(request,16L))
  env$mixedgp_available_cores <- function()56L
  expect_equal(resolve(56L),56L)
  expect_equal(resolve(17L),17L)
  env$mixedgp_available_cores <- function()2L
  expect_equal(resolve(56L),2L)
  expect_error(resolve(0L))
  expect_error(resolve(1.5))
  env$mixedgp_available_cores <- function()16L
  env$mixedgp_resolve_cores <- resolve
  enter <- mixedgp_budget_enter; environment(enter) <- env
  expect_message(close <- enter(56L),"requested 56; using")
  expect_equal(mixedgp_budget_cores(),if(.Platform$OS.type=="unix")16L else 1L)
  close()
  expect_null(mixedgp_budget())

})
