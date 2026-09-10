testthat::test_that("reference failures preserve RNG and return an explicit status", {
  set.seed(194);before<-.Random.seed
  x<-eivGP:::mixedgp_reference_task("reference",{runif(3);stop("budget exhausted")})
  testthat::expect_null(x$value)
  testthat::expect_identical(x$status$status,"failed")
  testthat::expect_identical(.Random.seed,before)
})
testthat::test_that("quadrature agrees with an analytic independent-proxy mean", {
  p<-list(q=2L,d=2L,m=4L,Omega=diag(2),A=matrix(0,2,2),tau=matrix(rep(qnorm(c(.25,.5,.75)),each=2),2,3),
    score_error="gaussian",lambda=1,scenario="primary")
  X<-rbind(c(.5,-.2),c(-.3,.4));C<-rbind(c(1L,4L),c(3L,2L))
  for(error in c("gaussian","logistic")){
    p$score_error<-error
    m<-eivGP:::oracle_m0_quadrature_2d(X,C,p)
    testthat::expect_equal(as.numeric(m),.3*X[,1]-.25*X[,2]+.1*X[,1]*X[,2],tolerance=1e-10)
    testthat::expect_true(all(is.na(attr(m,"truth_diagnostics")$mcse)))
  }
})
testthat::test_that("optimizer rescue is bounded and stops on the first success", {
  calls<-0L
  mock<-function(seed,n_starts){calls<<-calls+1L;if(n_starts<32L)stop("no valid fit");list(fit=list())}
  x<-eivGP:::mixedgp_retry_adapter(mock,list(seed=1L),"n_starts",c(8L,32L,64L))
  testthat::expect_identical(calls,2L)
  testthat::expect_identical(x$optimization_status,"rescued_converged")
  testthat::expect_equal(x$optimizer_attempts$seed,c(1,1001))
})
