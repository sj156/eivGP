## Run from repository root: Rscript revision/codes/tests/test_study2_sampler_repairs.R
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(script_arg) == 1L)
codes_dir <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
source(file.path(codes_dir, "00_study2_functions.R"))
expect_error <- function(expr, pattern) {
  ans <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(ans, "error"), grepl(pattern, conditionMessage(ans)))
}

set.seed(20509)
n <- 24L
X <- matrix(rnorm(2*n), n, 2L)
U <- matrix(rnorm(2*n), n, 2L)
y <- rnorm(n)
lt <- log(c(3, .4, .7, .5, .8))
max_error <- 0
for (kernel in c("se", "matern")) {
  for (b in c(1L,3L,11L)) {
    idx <- seq_len(b)
    evaluator <- make_gp_evaluator_general(y, X, kernel = kernel)
    block <- evaluator$prepare_block(U, lt, idx)
    for (j in seq_len(15L)) {
      proposed <- U
      proposed[idx,] <- proposed[idx,,drop=FALSE] + matrix(rnorm(2*b,sd=.2),b,2L)
      exact <- gp_loglik_integrated_general(y,X,proposed,lt,kernel)
      cached <- gp_loglik_from_moments_general(evaluator$evaluate(proposed,lt,block))
      max_error <- max(max_error,abs(exact-cached))
      stopifnot(abs(exact-cached) < 1e-8)
    }
    expect_error(evaluator$evaluate(U,lt+.01,block), "fixed inputs changed")
  }
}
cat("Maximum dense versus Schur log-likelihood error:",max_error,"\n")

## Count actual calls to Cholesky independently of the evaluator's counters.
isolated <- new.env(parent=globalenv())
source(file.path(codes_dir, "00_parallel_utils.R"),local=isolated)
source(file.path(codes_dir, "00_study2_functions.R"),local=isolated)
chol_calls <- 0L
original_chol <- isolated$safe_chol
isolated$safe_chol <- function(...) {
  chol_calls <<- chol_calls+1L
  original_chol(...)
}
ev <- isolated$make_gp_evaluator_general(y,X,calib_idx=2:n,U_obs=U)
state <- ev$evaluate(U,lt);ev$accept(state)
block <- ev$prepare_block(U,lt,1L)
proposed <- U;proposed[1,] <- proposed[1,]+.1
state <- ev$evaluate(proposed,lt,block);ev$accept(state)
before <- ev$counts()
set.seed(911)
noise_cached <- isolated$sample_sigma2_eps_general(y,X,proposed,lt,gp_evaluator=ev)
set.seed(911)
noise_dense <- sample_sigma2_eps_general(y,X,proposed,lt)
stopifnot(abs(noise_cached-noise_dense)<1e-12,
          is.null(ev$prepare_block(proposed,lt,integer(0))))
after <- ev$counts()
stopifnot(chol_calls==3L,
  sum(after[c("gp_full_cholesky","gp_block_setup_cholesky","gp_schur_cholesky")])==chol_calls,
  after["gp_dense_evaluations"]==before["gp_dense_evaluations"],
  after["gp_block_evaluations"]==before["gp_block_evaluations"],
  after["gp_cache_hits"]==before["gp_cache_hits"]+1L)
bad <- proposed;bad[2,1] <- bad[2,1]+1
expect_error(ev$evaluate(bad,lt),"calibrated latent input changed")
expect_error(ev$validate_block(2L),"calibrated row")
cat("Independent factorization counters and noise-update reuse: PASS\n")

## Check calibrated coordinates before every transition and density evaluation.
run_small <- function(missing_idx,strategy,use_schur=TRUE) {
  env <- new.env(parent=globalenv())
  source(file.path(codes_dir, "00_parallel_utils.R"),local=env)
  source(file.path(codes_dir, "00_study2_functions.R"),local=env)
  n <- 8L
  X <- matrix(seq(-1,1,length.out=n),ncol=1L)
  U <- matrix(seq(-1.5,1.5,length.out=n),ncol=1L)
  C <- matrix(rep(1:2,each=4),ncol=1L)
  y <- X[,1]+U[,1]^2+seq(.01,.08,length.out=n)
  calib <- setdiff(seq_len(n),missing_idx)
  checks <- 0L
  install <- function(name,u_arg,block_arg=NULL) {
    original <- get(name,env)
    wrapper <- function(...) {
      caller <- parent.frame()
      matched <- match.call(original,sys.call(),expand.dots=TRUE)
      u <- if(is.null(u_arg)) get("U_curr",caller,inherits=FALSE) else {
        eval(matched[[u_arg]],caller)
      }
      checks <<- checks+1L
      stopifnot(all(u[calib,,drop=FALSE]==U[calib,,drop=FALSE]))
      if(!is.null(block_arg)) {
        block <- eval(matched[[block_arg]],caller)
        stopifnot(all(block %in% missing_idx),!anyDuplicated(block))
      }
      original(...)
    }
    assign(name,wrapper,env)
  }
  for(name in c("sample_scores_ord","update_A_ord",
                "update_logtheta_slice_general","update_measurement_marginal_slice",
                "sample_sigma2_eps_general","ordinal_loglik_marginal")) install(name,"U")
  install("update_tau_ord",NULL)
  install("update_U_ess_block_general","U_curr","block_idx")
  install("update_U_theta_ess_integrated_general","U_curr","block_idx")
  fit <- env$fit_eivgp_ordprobit_fb(X,y,C,U_obs=U,calib_idx=calib,d=1L,
    n_iter=60L,burn=20L,thin=1L,n_chains=1L,parallel_chains=FALSE,
    preset="fast",seed=712L,sampler_strategy=strategy,
    control_overrides=if(is.null(use_schur)) list() else list(gp_use_block_schur=use_schur))
  stopifnot(all(vapply(missing_idx,function(i)
    var(fit$mcmc$samples_U[,i,1L])>0,logical(1))))
  list(fit=fit,checks=checks)
}
rows <- list()
for(strategy in c("interwoven","legacy")) {
  for(missing_idx in list(integer(0),8L,c(2L,4L,8L),seq_len(8L))) {
    schur <- run_small(missing_idx,strategy,TRUE)
    dense <- run_small(missing_idx,strategy,FALSE)
    stopifnot(max(abs(schur$fit$mcmc$samples_U-dense$fit$mcmc$samples_U))<1e-8,
      max(abs(schur$fit$mcmc$samples_logtheta-dense$fit$mcmc$samples_logtheta))<1e-8,
      max(abs(schur$fit$mcmc$samples_sigma2-dense$fit$mcmc$samples_sigma2))<1e-8)
    counts <- schur$fit$diagnostics$summary
    if(length(missing_idx)==1L) stopifnot(counts$gp_block_evaluations>0L)
    rows[[length(rows)+1L]] <- data.frame(strategy=strategy,
      missing=length(missing_idx),pretransition_checks=schur$checks,
      gp_evaluations=counts$total_gp_evaluations,
      full_cholesky=counts$gp_full_cholesky,cache_hits=counts$gp_cache_hits)
  }
}
print(do.call(rbind,rows),row.names=FALSE)
cat("Zero/one/multiple missing regressions and dense/Schur path agreement: PASS\n")
auto_small <- run_small(8L,"interwoven",NULL)$fit
stopifnot(identical(auto_small$control$gp_block_schur_mode,"automatic"),
  identical(auto_small$control$gp_use_block_schur,FALSE),
  auto_small$diagnostics$summary$gp_block_evaluations==0L,
  auto_small$diagnostics$summary$gp_block_setup_cholesky==0L)
cat("Automatic small-n dense route and explicit small-n Schur override: PASS\n")

## Legacy moves must fail explicitly on exhaustion or a nonfinite current target.
X <- matrix(c(-.4,.7),ncol=1);U <- matrix(c(.3,.5),ncol=1)
A <- matrix(1.6,1,1);S <- matrix(c(-.2,.9),ncol=1)
lt <- log(c(2,.8,1.2));y <- c(-.7,.9)
expect_error(update_U_ess_block_general(y,X,U,S,A,lt,.1,2L,max_try=0L),"exceeded max_try")
expect_error(update_U_ess_block_general(c(NA,.9),X,U,S,A,lt,.1,2L),"non-finite density")

## Solver/degeneracy warnings cannot produce output labelled exact.
for(message in c("Did not find a solution to the nonlinear system in `mvrandn`!",
                 "Some variables have a degenerate distribution.",
                 "Method may fail as covariance matrix is singular!",
                 "Sample of size smaller than n returned.","Unexpected numerical warning")) {
  expect_error(checked_minimax_tilting_draw(function(){warning(message);1}),
               "cannot be certified exact")
}
low <- checked_minimax_tilting_draw(function(){
  warning("Acceptance probability smaller than 0.001");1
})
stopifnot(identical(low$draws,1),length(low$low_acceptance_warnings)==1L)
if(requireNamespace("TruncatedNormal",quietly=TRUE)) {
  expect_error(sample_u_given_c_ordprobit_minimax(matrix(2L,1L,1L),
    matrix(1,1L,1L),list(c(0,1e-11)),3L),"cannot be certified exact")
  draws <- sample_u_given_c_ordprobit_minimax(matrix(c(1L,2L),2L,1L),
    matrix(1,1L,1L),list(0),2L)
  stopifnot(identical(attr(draws,"latent_sampler"),"exact_minimax_tilting"),
    all(attr(draws,"pattern_telemetry")$low_acceptance_warning_count==0L))
}
cat("Explicit ESS failure and minimax warning classification: PASS\n")

## Independent quadrature target for a Schur-evaluated scalar U transition.
set.seed(912)
X <- matrix(c(-.4,.1,.5,.7),ncol=1L)
y <- c(-.7,.2,.3,.9)
U0 <- matrix(c(.3,-.3,.1,.5),ncol=1L)
C <- matrix(c(1L,1L,2L,2L),ncol=1L)
A <- matrix(1.6,1L,1L);tau <- list(.2)
lt <- log(c(2,.8,1.2));prior <- make_gp_prior(1L,1L)
grid <- seq(-7,7,length.out=7001L)
lp <- vapply(grid,function(u){
  up <- U0;up[4,1] <- u
  dnorm(u,log=TRUE)+gp_loglik_integrated_general(y,X,up,lt)+
    ordinal_loglik_marginal(C,up,A,tau)
},numeric(1))
w <- exp(lp-max(lp));w <- w/sum(w)
truth <- sum(grid*w)
values <- numeric(12000L)
u <- U0
ev <- make_gp_evaluator_general(y,X,calib_idx=1:3,U_obs=U0)
for(i in seq_along(values)) {
  ans <- update_U_theta_ess_integrated_general(y,X,C,u,lt,prior,4L,
    reference="prior",A=A,tau=tau,update_theta=FALSE,gp_evaluator=ev)
  u <- ans$U;values[i] <- u[4,1]
}
values <- tail(values,10000L)
batch <- colMeans(matrix(values,nrow=250L))
mcse <- sd(batch)/sqrt(length(batch))
stopifnot(abs(mean(values)-truth)<5*mcse+.005)
cat("Independent scalar posterior mean:",mean(values),"; grid:",truth,
    "; batch MCSE:",mcse,"\n")
cat("All Study II sampler repair tests passed.\n")
