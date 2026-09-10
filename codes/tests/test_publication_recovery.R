source("codes/simulation_helpers.R")
source("codes/publication_recovery.R")
e <- mixedgp_simulation_engine("codes")
# Reference errors preserve RNG and never become fictitious numerical truth.
set.seed(718); before <- .Random.seed
bad <- e$mixedgp_reference_task("oracle",{runif(10);stop("rare pattern exhausted")})
stopifnot(is.null(bad$value),bad$status$status=="failed",identical(.Random.seed,before))
good <- e$mixedgp_reference_task("oracle",{runif(1);42})
stopifnot(good$value==42,identical(.Random.seed,before))
# Independent analytic reference when proxies carry no information about U.
d <- e$make_study2_synthetic_dataset(1,n=30,n_test=8,q=4,calib_grid=6L)$data
p <- d$true_params;p$A[]<-0
for(error in c("gaussian","logistic")){
 p$score_error<-error
 m<-e$oracle_m0_quadrature_2d(d$test$X,d$test$C,p)
 X<-d$test$X;truth<-.30*X[,1]-.25*X[,2]+.10*X[,1]*X[,2]
 stopifnot(max(abs(m-truth))<1e-10,all(is.na(attr(m,"truth_diagnostics")$mcse)))
}
# Failure at an inadequate refinement budget is explicit, not a plausible number.
stopifnot(inherits(try(e$oracle_m0_quadrature_2d(d$test$X,d$test$C,d$true_params,
  orders=c(2L,3L),tolerance=1e-14),silent=TRUE),"try-error"))
# Coverage audit catches duplicates even when there are 50 rows, missing IDs,
# and a non-finite predictive score. A diagnostic flag alone does not exclude.
cfg<-study2_simulation_config("publication",code_dir="codes",core_budget=1L)
cfg$cells<-cfg$cells[1L]
base<-data.frame(cell_id=cfg$cells[[1]]$id,rep=1:50,n_calib=50L,method="EIV-GP",
 RMSE=.2,MAE=.1,CRPS=.1,NLPD=0,Coverage95=.95,Width95=1,IntervalScore95=1,
 diagnostic_warning=TRUE)
comp<-do.call(rbind,lapply(c("UC-GP","LVGP","EzGP"),function(method){z<-base;z$method<-method;z$n_calib<-NA_integer_;z}))
stopifnot(all(mixedgp_recovery_audit(base,comp,cfg)$complete))
b<-base;b$rep[50]<-49L;stopifnot(!mixedgp_recovery_audit(b,comp,cfg)$complete[1])
b<-base;b$CRPS[1]<-NA;stopifnot(mixedgp_recovery_audit(b,comp,cfg)$missing_ids[1]=="1")
stopifnot(nrow(mixedgp_recovery_upsert(base,base[1,,drop=FALSE],c("cell_id","rep","method","n_calib")))==50L)
# Bounded rescues keep the original successful first attempt and record failures.
calls<-0L
mock<-function(seed,n_starts){calls<<-calls+1L;if(n_starts<32L)stop("no converged start");list(fit=list(),value=seed)}
x<-e$mixedgp_retry_adapter(mock,list(seed=100L),"n_starts",c(8L,32L,64L))
stopifnot(calls==2L,x$value==1100L,nrow(x$optimizer_attempts)==2L,x$optimization_status=="rescued_converged")
calls<-0L;x<-e$mixedgp_retry_adapter(mock,list(seed=100L),"n_starts",c(32L,64L))
stopifnot(calls==1L,x$value==100L)
# Genuine tiny end-to-end driver: all three reference tasks throw before fitting.
root<-tempfile();dir.create(root)
cfg<-study2_simulation_config("development",code_dir="codes",core_budget=1L,
 output_root=root,data_root=file.path(root,"data"))
cfg$cells<-cfg$cells[1L];cell<-cfg$cells[[1]]
cell$n_rep<-1L;cell$n<-30L;cell$n_test<-8L;cell$calibration_grid<-6L
cell$run_ablations<-FALSE;cell$evaluate_f<-FALSE;cell$evaluate_u<-FALSE;cfg$cells[[1]]<-cell
cfg$mcmc$n_iter<-12L;cfg$mcmc$burn<-4L;cfg$mcmc$n_chains<-2L
cfg$evaluation<-list(n_pred_draw=4L,n_m_eval=2L,n_m_draw=4L,n_m_latent=4L,n_m_truth=16L,n_oracle_pool=30L)
mixedgp_generate_cell_data(cfg,cell,e)
for(name in c("make_oracle_pool_2d","sample_oracle_test_y_2d","oracle_m0_quadrature_2d"))
 assign(name,function(...)stop("deliberate reference failure"),e)
e$STUDY2_RECOVERY_MODEL<-list(priors=list(signal_shape=c(11,4)),control=list(),kernel=list(name="se",matern_nu=2.5))
result<-suppressWarnings(mixedgp_run_study2_cell(cfg,cell,e,root))
stopifnot(all(result$outputs$reference_status$status=="failed"),
 any(result$outputs$predictive_metrics$method=="EIV-GP"),
 !any(result$outputs$predictive_metrics$method=="Oracle"),nrow(result$outputs$mean_recovery)==0L)
files<-list.files(root,pattern="^recovery_fit_.*rds$",recursive=TRUE,full.names=TRUE)
stopifnot(length(files)==1L,identical(readRDS(files[1])$fit$priors$signal_shape,c(11,4)))
# Resume must use saved work, even if sampling is unavailable on the next call.
e$fit_eivgp_ordprobit_fb<-function(...)stop("Sampler must not run on resume")
again<-suppressWarnings(mixedgp_run_study2_cell(cfg,cell,e,root))
stopifnot(identical(result$outputs$predictive_metrics,again$outputs$predictive_metrics))
message("Recovery checks passed: analytic truth, RNG isolation, completeness, rescues, nonfatal reference failures, fit checkpoint/resume.")
