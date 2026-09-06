## Run from the project root:
## Rscript revision/codes/tests/test_study1_sampler_repairs.R
## These bounded checks validate transition targets and implementation;
## they are not convergence certification for a scientific analysis.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(script_arg) == 1L)
codes_dir <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
source(file.path(codes_dir, "00_study1_functions.R"))
stop_if_not_error <- function(expr, pattern) {
  result <- tryCatch(force(expr), error = conditionMessage)
  stopifnot(is.character(result), grepl(pattern, result, fixed = TRUE))
}
mcse_batch <- function(x, batch = 100L) {
  n <- length(x) %/% batch
  sd(colMeans(matrix(x[seq_len(n * batch)], nrow = batch))) / sqrt(n)
}

## Feasible tiny gaps and shared calibration gaps initialize in support.
tight_u <- c(.2, .200000001)
tight_tau <- initialize_tau_1d(1:2, tight_u, 1:2, 2L)
stopifnot(check_constraints_1d(tight_u, 1:2, tight_tau))
shared_u <- c(.2, NA_real_, NA_real_, .200000001)
shared_tau <- initialize_tau_1d(1:4, shared_u, c(1L,4L), 4L)
stopifnot(all(diff(shared_tau) > 0), shared_tau[1] > shared_u[1],
          shared_tau[3] < shared_u[4])
stop_if_not_error(initialize_tau_1d(1:2, c(.4,.3), 1:2, 2L),
                  "incompatible")
stop_if_not_error(assert_threshold_state_1d(c(-.5,.5+1e-12),0,1:2,2L,
                    calib_idx=1:2,u_obs=c(-.5,.5)), "calibrated")
stop_if_not_error(assert_threshold_state_1d(c(0,.5),0,c(2L,2L),2L),
                  "interval")

## Slice helpers neither shift a finite current state nor erode support.
seen <- numeric(0)
invisible(bounded_slice_update_1d(.2, function(x) {
  seen <<- c(seen,x); if(x>=.2 && x<=.200000001) 0 else -Inf
}, w=1e-10, lower=.2,upper=.200000001,fail_on_limit=TRUE))
stopifnot(identical(seen[1], .2))
seen <- numeric(0)
invisible(full_interval_slice_update(.2, function(x) {
  seen <<- c(seen,x); 0
}, lower=.2,upper=.200000001,fail_on_limit=TRUE))
stopifnot(identical(seen[1], .2))

## Flat residual likelihood isolates the latent prior in each coordinate.
flat_env <- new.env(parent=globalenv())
flat_env$gp_state_1d <- function(...) list(loglik=0)
flat_local <- update_u_local_z_slice_1d
environment(flat_local) <- flat_env
N <- 6000L
target_rows <- list()
for (case in list(list(tau=0,c=1L,initial=-.5,lo=-Inf,hi=0),
                 list(tau=c(-.5,.8),c=2L,initial=.2,lo=-.5,hi=.8),
                 list(tau=7.99,c=2L,initial=8.05,lo=7.99,hi=Inf),
                 list(tau=-7.99,c=1L,initial=-8.05,lo=-Inf,hi=-7.99),
                 list(tau=c(.2,.200000001),c=2L,initial=.2000000005,
                      lo=.2,hi=.200000001))) {
  set.seed(617L + length(target_rows))
  value <- case$initial
  draws <- numeric(N)
  for (i in seq_len(N)) {
    value <- flat_local(0,value,list(),case$c,case$tau,c(0,0,0),1,NULL,1L)$u
    draws[i] <- value
  }
  stopifnot(all(draws>case$lo), all(draws<=case$hi))
  if (is.finite(case$lo) && is.finite(case$hi) && case$hi-case$lo < 1e-8) {
    expected <- (case$hi+case$lo)/2
    tolerance <- 0.06*(case$hi-case$lo)
  } else {
    mass <- if(case$lo>=0) pnorm(case$lo,lower.tail=FALSE)-pnorm(case$hi,lower.tail=FALSE) else
      pnorm(case$hi)-pnorm(case$lo)
    expected <- (dnorm(case$lo)-dnorm(case$hi))/mass
    tolerance <- 7*mcse_batch(draws)+1e-4
  }
  stopifnot(abs(mean(draws)-expected) < tolerance)
  target_rows[[length(target_rows)+1L]] <- data.frame(lower=case$lo,upper=case$hi,
      empirical_mean=mean(draws),reference=expected,tolerance=tolerance)
}
print(do.call(rbind,target_rows),row.names=FALSE)

## One actual missing GP input: numerical quadrature is an independent target.
set.seed(621)
gx <- matrix(c(-1,0,1),ncol=1)
Dx <- list(pairwise_sqdist(gx))
u <- c(-.7,.3,.9)
y <- c(-.8,.1,.7)
lt <- log(c(2,.7,.4))
spec <- make_theta_spec_multix(1L)
log_density <- function(z) vapply(z, function(value) {
  proposed <- u; proposed[1] <- value
  A <- diag(3) + 4*exp(-.7*outer(gx[,1],gx[,1],"-")^2-
                        .4*outer(proposed,proposed,"-")^2)
  -.5*as.numeric(determinant(A,logarithm=TRUE)$modulus)-
    drop(crossprod(y,solve(A,y)))/.2 + dnorm(value,log=TRUE)
},numeric(1))
normalizer <- integrate(function(z) exp(log_density(z)),-8,0,rel.tol=1e-10)$value
reference_mean <- integrate(function(z) z*exp(log_density(z)),-8,0,
                            rel.tol=1e-10)$value/normalizer
actual_draws <- numeric(N)
state <- u
for (i in seq_len(N)) {
  state <- update_u_local_z_slice_1d(y,state,Dx,c(1L,2L,2L),0,lt,.1,spec,1L)$u
  actual_draws[i] <- state[1]
}
stopifnot(abs(mean(actual_draws)-reference_mean) < 7*mcse_batch(actual_draws))
cat("Actual GP conditional mean/reference/MCSE:",mean(actual_draws),reference_mean,
    mcse_batch(actual_draws),"\n")
ess_draws <- numeric(N)
state <- u
for (i in seq_len(N)) {
  state <- update_u_ess_block_1d(y,state,Dx,c(1L,2L,2L),0,lt,.1,spec,1L)$u
  ess_draws[i] <- state[1]
}
stopifnot(abs(mean(ess_draws)-reference_mean) < 7*mcse_batch(ess_draws))
cat("Actual GP ESS mean/reference/MCSE:",mean(ess_draws),reference_mean,
    mcse_batch(ess_draws),"\n")
noise_A <- diag(3)+4*exp(-.7*Dx[[1]]-.4*outer(u,u,"-")^2)
noise_rate <- .05+drop(crossprod(y,solve(noise_A,y)))/2
noise_mean <- noise_rate/(2+3/2-1)
noise_draws <- replicate(N,sample_sigma2_eps_1d(y,u,Dx,lt,spec))
stopifnot(abs(mean(noise_draws)-noise_mean) < 7*sd(noise_draws)/sqrt(N))

## Cache hits retain the sigma dependence and invalidate every model input.
cache <- new_gp_cache_1d()
s1 <- gp_state_1d(y,u,Dx,lt,.1,spec,gp_cache=cache)
s2 <- gp_state_1d(y,u,Dx,lt,.3,spec,gp_cache=cache)
stopifnot(cache$full_factorizations==1L,cache$cache_hits==1L)
stopifnot(isTRUE(all.equal(s2$loglik,gp_state_1d(y,u,Dx,lt,.3,spec)$loglik)))
before <- cache$full_factorizations
invisible(sample_sigma2_eps_1d(y,u,Dx,lt,spec,gp_cache=cache))
stopifnot(cache$full_factorizations==before)
for (args in list(list(y=y+c(.01,0,0),u=u,Dx_list=Dx,logtheta=lt),
                 list(y=y,u=u+c(.01,0,0),Dx_list=Dx,logtheta=lt),
                 list(y=y,u=u,Dx_list=list(Dx[[1]]*1.1),logtheta=lt),
                 list(y=y,u=u,Dx_list=Dx,logtheta=lt+c(.01,0,0)))) {
  before <- cache$full_factorizations
  cached <- do.call(gp_state_1d,c(args,list(sigma2_eps=.1,gp_cache=cache)))
  direct <- do.call(gp_state_1d,c(args,list(sigma2_eps=.1)))
  stopifnot(cache$full_factorizations==before+1L,
            isTRUE(all.equal(cached$loglik,direct$loglik)))
}

## Exact small-block Schur calculations match dense likelihoods across
## kernels, block sizes and allowed hyperparameter extremes.
set.seed(622)
n <- 60L
X <- matrix(rnorm(2*n),n,2L)
U <- rnorm(n)
Y <- rnorm(n)
D <- lapply(1:2,function(j) pairwise_sqdist(X[,j,drop=FALSE]))
spec2 <- make_theta_spec_multix(2L)
maximum_error <- 0
for (kernel in c("se","matern")) for (nu in c(.5,1.5,2.5)) {
  for (theta in list(c(3,.4,.7,.8),c(.05,1e-4,100,1e-4),c(100,1e-4,1e-4,1e-4))) {
    for (idx in list(37L,c(2L,17L,44L),seq_len(7L))) {
      block_cache <- new_gp_cache_1d()
      prepared <- prepare_gp_block_1d(Y,U,D,log(theta),spec2,idx,kernel,nu,block_cache)
      stopifnot(!is.null(prepared))
      for (r in 1:4) {
        proposed <- U; proposed[idx] <- U[idx]+rnorm(length(idx),sd=.3)
        accelerated <- gp_block_state_1d(prepared,proposed[idx],.1,block_cache)
        dense <- gp_state_1d(Y,proposed,D,log(theta),.1,spec2,kernel,nu)
        err <- abs(accelerated$loglik-dense$loglik)
        maximum_error <- max(maximum_error,err)
        stopifnot(err < 2e-7,
                  max(abs(accelerated$A-dense$A)) < 1e-9)
        ## The accepted block summary serves the next variance update.
        before <- block_cache$full_factorizations
        invisible(sample_sigma2_eps_1d(Y,proposed,D,log(theta),spec2,
                                      kernel=kernel,matern_nu=nu,gp_cache=block_cache))
        stopifnot(block_cache$full_factorizations==before)
      }
    }
  }
}
cat("Maximum Schur versus dense log-likelihood error:",maximum_error,"\n")
## If the shortcut factor fails, recompute the same candidate densely.
bad_prepared <- prepared
bad_prepared$chol_rest[,] <- NA_real_
proposed <- U; proposed[idx] <- proposed[idx]+.01
fallback_count <- block_cache$block_fallbacks
fallback <- gp_block_state_1d(bad_prepared,proposed[idx],.1,block_cache)
dense <- gp_state_1d(Y,proposed,D,log(theta),.1,spec2,kernel,nu)
stopifnot(block_cache$block_fallbacks==fallback_count+1L,
          isTRUE(all.equal(fallback$loglik,dense$loglik)))

## Fit edge cases, asserting calibration before/after every latent transition.
fit_env <- new.env(parent=globalenv())
source(file.path(codes_dir, "00_study1_functions.R"),local=fit_env)
u8 <- c(-1.5,-1,-.5,-.2,.2,.5,1,1.5)
missing <- integer(0)
calibrated <- seq_len(8)
checks <- 0L
check_calibration <- function(value) {
  stopifnot(identical(as.numeric(value[calibrated]),as.numeric(u8[calibrated])))
}
for (name in c("update_u_ess_block_1d","update_u_local_z_slice_1d")) {
  transition <- get(name,envir=fit_env)
  index_name <- if(name=="update_u_ess_block_1d") "block_idx" else "update_idx"
  wrapped <- local({fn <- transition; key <- index_name; function(...) {
    args <- list(...)
    stopifnot(all(args[[key]] %in% missing))
    check_calibration(args$u)
    ans <- do.call(fn,args)
    check_calibration(ans$u)
    checks <<- checks+1L
    ans
  }})
  assign(name,wrapped,envir=fit_env)
}
fit_rows <- list()
for (indices in list(integer(0),8L,1L,c(2L,5L,8L),seq_len(8))) {
  missing <- indices; calibrated <- setdiff(seq_len(8),missing)
  obs <- u8; obs[missing] <- NA_real_
  checks <- 0L
  fit <- fit_env$fit_eivgp_1d(1:8,c(0,1,0,1,2,3,2,3),rep(1:2,each=4),
        u_obs=obs,m=2L,n_iter=120L,burn=20L,n_chains=2L,
        parallel_chains=FALSE,seed=623L)
  stopifnot(nrow(fit$mcmc$samples_u)==200L,
            all(is.finite(fit$mcmc$samples_u)),all(fit$mcmc$samples_sigma2>0),
            identical(fit$control$gp_block_schur,FALSE),
            identical(fit$control$gp_block_schur_requested,"auto"),
            sum(fit$mcmc$chain_stats$gp_block_evaluations)==0L)
  for (i in seq_len(200L)) {
    check_calibration(fit$mcmc$samples_u[i,])
    stopifnot(check_constraints_1d(fit$mcmc$samples_u[i,],rep(1:2,each=4),
                                    fit$mcmc$samples_tau[i,]))
  }
  stopifnot(sum(fit$mcmc$chain_stats$gp_cache_hits)>0)
  fit_rows[[length(fit_rows)+1L]] <- data.frame(missing=length(missing),
      saved=200L,checked_transitions=checks,cache_hits=sum(fit$mcmc$chain_stats$gp_cache_hits))
}
print(do.call(rbind,fit_rows),row.names=FALSE)

## A moderate fit exercises production Schur proposals and saved invariants.
set.seed(624)
um <- sort(rnorm(40)); cm <- as.integer(um>0)+1L
om <- um; om[c(3L,10L,28L,37L)] <- NA_real_
moderate <- fit_eivgp_1d(matrix(rnorm(80),40,2),sin(um)+rnorm(40,sd=.3),cm,
             u_obs=om,m=2L,n_iter=70L,burn=20L,n_chains=2L,
             parallel_chains=FALSE,seed=624L,kernel="matern",matern_nu=2.5,
             gp_block_schur=TRUE)
stopifnot(sum(moderate$mcmc$chain_stats$gp_block_evaluations)>0,
          all(is.finite(moderate$mcmc$samples_u)),
          identical(moderate$control$gp_block_schur_requested,"enabled"))
cached_dense <- fit_eivgp_1d(moderate$data$x_raw,moderate$data$y_raw,cm,
             u_obs=om,m=2L,n_iter=70L,burn=20L,n_chains=2L,
             parallel_chains=FALSE,seed=624L,kernel="matern",matern_nu=2.5,
             gp_block_schur=FALSE)
stopifnot(identical(cached_dense$control$gp_block_schur,FALSE),
          sum(cached_dense$mcmc$chain_stats$gp_block_evaluations)==0L,
          sum(cached_dense$mcmc$chain_stats$gp_cache_hits)>0L,
          all(is.finite(cached_dense$mcmc$samples_u)))
stop_if_not_error(fit_eivgp_1d(1:8,1:8,rep(1:2,each=4),gp_block_schur=NA),
                  "gp_block_schur")
for (bad_flag in list(logical(0),c(TRUE,FALSE),1,"auto",list(TRUE))) {
  stop_if_not_error(fit_eivgp_1d(1:8,1:8,rep(1:2,each=4),gp_block_schur=bad_flag),
                    "gp_block_schur")
}
auto_small <- fit_eivgp_1d(moderate$data$x_raw,moderate$data$y_raw,cm,
             u_obs=om,m=2L,n_iter=70L,burn=20L,n_chains=2L,
             parallel_chains=FALSE,seed=624L,kernel="matern",matern_nu=2.5)
stopifnot(identical(auto_small$control$gp_block_schur,FALSE),
          identical(auto_small$control$gp_block_schur_requested,"auto"),
          identical(auto_small$mcmc$samples_u,cached_dense$mcmc$samples_u),
          identical(auto_small$mcmc$samples_logtheta,cached_dense$mcmc$samples_logtheta))
## A tiny all-calibrated fit verifies the exact default threshold without
## spending simulation effort on a mixing comparison.
for (sample_size in c(119L,120L)) {
  known_u <- seq(-2,2,length.out=sample_size)
  auto_fit <- fit_eivgp_1d(known_u,sin(known_u),as.integer(known_u>0)+1L,
                u_obs=known_u,m=2L,n_iter=3L,burn=1L,n_chains=1L,
                parallel_chains=FALSE,seed=625L)
  stopifnot(identical(auto_fit$control$gp_block_schur,sample_size>=120L),
            identical(auto_fit$control$gp_block_schur_requested,"auto"))
}
print(moderate$mcmc$chain_stats[,c("gp_full_factorizations","gp_cache_hits",
       "gp_block_setups","gp_block_evaluations","gp_block_fallbacks")],row.names=FALSE)
cat("PASS: Study I support repairs, calibration invariants, conditional targets, cache and exact blocks.\n")
