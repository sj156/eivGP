## End-to-end exports for BOTH studies, with fixture fits and temporary outputs.
## No optimizer or MCMC runs; no real frozen data or reports are changed.
source("codes/tests/test_competitor_cache.R")
root <- tempfile("competitor-export-test-"); dir.create(root)
paper <- file.path(root,"paper"); dir.create(paper)
Sys.setenv(EIVGP_COMPETITOR_CACHE=e$.mixedgp_competitor_cache_root,
  EIVGP_COMPETITOR_REPORT_ROOT=file.path(root,"reports"),
  EIVGP_OVERLEAF_ROOT=paper,EIVGP_WORKERS="1")
for (study in c("study1","study2")) {
  data_root <- file.path(root,"data",study)
  Sys.setenv(EIVGP_DATA_ROOT=data_root)
  cfg <- get(paste0(study,"_simulation_config"))("development",code_dir="codes",data_root=data_root)
  for (cell in cfg$cells) {
    mixedgp_generate_cell_data(cfg,cell,e)
    for (rep_id in seq_len(cell$n_rep)) {
      d <- mixedgp_read_cell_replication(cfg,cell,rep_id)$data
      X <- if(study=="study1") matrix(d$train$x) else d$train$X
      C <- if(study=="study1") matrix(d$train$c) else d$train$C
      Xt <- if(study=="study1") matrix(d$test$x) else d$test$X
      Ct <- if(study=="study1") matrix(d$test$c) else d$test$C
      seed <- if(study=="study1") 250000L+1000L*rep_id else
        1000000L*match(cell$scenario,c("primary","latent_additive_control","high_uncertainty","logistic_misspec"))+10000L*rep_id+100L
      e$mixedgp_cached_competitors(X,d$train$y,C,Xt,Ct,n_draw=1L,seed=seed,
        m_vec=rep(cell$m,ncol(C)),controls=e$mixedgp_competitor_protocol(study),allow_fit=TRUE)
    }
  }
  log <- file.path(root,paste0(study,".log"))
  status <- system2(file.path(R.home("bin"),"Rscript"),
    c("--vanilla","experiments/run_competitors.R",study,"export","development"),stdout=log,stderr=log)
  if(status!=0L) stop(paste(readLines(log),collapse="\n"))
  target <- file.path(paper,"tables","competitors",study,"development")
  s <- read.csv(file.path(target,"statuses.csv"))
  summary <- read.csv(file.path(target,"predictive_summary.csv"))
  stopifnot(all(s$status=="success"),all(s$cache_hit),nrow(summary)>0L,
    all(summary$n_success==3L),file.exists(file.path(target,"predictive_summary.tex")))
}
message("Both standalone study exports passed using fixture caches and temporary paper folders.")
