source("codes/simulation_helpers.R")
root<-tempfile();dir.create(root)
cfg<-study2_simulation_config("development",code_dir="codes",core_budget=1L,
  output_root=root,data_root=file.path(root,"data"))
cfg$cells<-cfg$cells[1L]
cfg$cells[[1]]$n_rep<-1L;cfg$cells[[1]]$n<-30L;cfg$cells[[1]]$n_test<-8L
cfg$cells[[1]]$calibration_grid<-6L
cfg$cells[[1]]$run_ablations<-FALSE;cfg$cells[[1]]$evaluate_u<-FALSE
cfg$mcmc$n_iter<-12L;cfg$mcmc$burn<-4L;cfg$mcmc$n_chains<-2L
cfg$evaluation<-list(n_pred_draw=4L,n_m_eval=2L,n_m_draw=4L,n_m_latent=4L,
  n_m_truth=16L,n_oracle_pool=30L)
e<-mixedgp_simulation_engine("codes")
e$mixedgp_cached_competitors<-function(...) stop("MCMC accessed competitor cache")
mixedgp_generate_cell_data(cfg,cfg$cells[[1]],e)
result<-mixedgp_run_study2_cell(cfg,cfg$cells[[1]],e,root)
stopifnot(nrow(result$outputs$competitor_status)==0L,
  file.exists(file.path(root,"cells",cfg$cells[[1]]$id,"dataset_identity.csv")))
message("Study II tiny independent fit passed without accessing competitors.")
