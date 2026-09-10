source("codes/simulation_helpers.R")
e <- mixedgp_simulation_engine("codes")
e$.mixedgp_competitor_cache_root <- tempfile("competitor-test-")
e$calls <- 0L
e$fail <- FALSE
## Mock fitting, not the cache, to check reuse without expensive optimization.
e$run_published_mixedgp_competitors <- eval(quote(function(X_train,y_train,C_train,
    X_test,C_test,n_draw,seed,m_vec,methods,strict,controls) {
  calls <<- calls + 1L
  if (fail) stop("deliberate optimizer failure")
  method <- methods[1L]
  mu <- rep(mean(y_train),nrow(X_test)); variance <- rep(1,nrow(X_test))
  list(status=mixedgp_competitor_status(method,status="success"),
    predictive_means=setNames(list(mu),method),
    predictive_variances=setNames(list(variance),method),
    latent_means=setNames(list(mu),method), fits=setNames(list(list(seed=seed)),method))
}),e)
a <- list(X_train=matrix(seq_len(12)/12),y_train=seq_len(12)/12,
  C_train=matrix(rep(1:3,4)),X_test=matrix(c(.1,.9)),C_test=matrix(c(1,3)),
  n_draw=10L,seed=100L,m_vec=3L,controls=e$mixedgp_competitor_protocol("study1"))
run <- function(...) do.call(e$mixedgp_cached_competitors,modifyList(a,list(...)))
miss <- run()
stopifnot(e$calls==0L,all(miss$status$optimization_status=="missing_cache"))
first <- run(allow_fit=TRUE)
stopifnot(e$calls==3L,all(first$status$status=="success"),!any(first$status$cache_hit))
subset <- run(n_draw=4L,methods=c("EzGP","UC-GP"))
stopifnot(e$calls==3L,all(subset$status$cache_hit),
  identical(dim(subset$draws[[1L]]),c(4L,2L)),
  identical(attr(subset$draws[[1L]],"conditional_vars"),matrix(1,nrow=1,ncol=2)))
changed <- run(y_train=a$y_train+1)
stopifnot(e$calls==3L,all(changed$status$optimization_status=="missing_cache"))
e$fail <- TRUE
failed <- run(seed=200L,allow_fit=TRUE)
stopifnot(e$calls==6L,all(failed$status$status=="unavailable_or_failed"))
e$fail <- FALSE
again <- run(seed=200L,allow_fit=TRUE)
stopifnot(e$calls==6L,all(again$status$cache_hit))
retried <- run(seed=200L,allow_fit=TRUE,retry_failed=TRUE)
stopifnot(e$calls==9L,all(retried$status$status=="success"))
unchanged <- run(seed=200L,allow_fit=TRUE,retry_failed=TRUE)
stopifnot(e$calls==9L,all(unchanged$status$cache_hit),nrow(run(methods=character())$status)==0L)
message("Competitor cache: read-only misses, reuse, selection, draw budgets, invalidation and retry passed.")

# A locked method cannot erase successes from the same dataset.
miss <- run(seed=300L)
lock <- paste0(miss$status$cache_file[miss$status$method=="LVGP"],".lock")
dir.create(dirname(lock),recursive=TRUE,showWarnings=FALSE);dir.create(lock)
locked <- run(seed=300L,allow_fit=TRUE)
stopifnot(locked$status$optimization_status[locked$status$method=="LVGP"]=="cache_locked",
  all(locked$status$status[locked$status$method!="LVGP"]=="success"),dir.exists(lock),
  setequal(names(locked$draws),c("UC-GP","EzGP")))
unlink(lock,recursive=TRUE)
# Cleanup must not remove a lock subsequently owned by another coordinator.
dir.create(lock);saveRDS(list(token="new-owner"),file.path(lock,"owner.rds"))
e$mixedgp_release_competitor_lock(lock,"old-owner");stopifnot(dir.exists(lock))
e$mixedgp_release_competitor_lock(lock,"new-owner");stopifnot(!dir.exists(lock))
message("Cache-lock isolation and ownership checks passed.")
