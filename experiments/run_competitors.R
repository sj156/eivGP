## Rscript --vanilla experiments/run_competitors.R study1 [plan|run|retry|export] [publication|development]
## Run install_eivgp_dependencies.R explicitly before this module. Never runs MCMC.
args <- commandArgs(trailingOnly=TRUE)
study <- match.arg(if(length(args)) args[1L] else "study1", c("study1","study2"))
action <- match.arg(if(length(args)>1L) args[2L] else "plan", c("plan","run","retry","export"))
mode <- match.arg(if(length(args)>2L) args[3L] else "publication", c("publication","development"))
cli <- grep("^--file=",commandArgs(FALSE),value=TRUE)
repo <- normalizePath(file.path(dirname(sub("^--file=","",cli[1L])),".."),mustWork=TRUE)
Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",VECLIB_MAXIMUM_THREADS="1")
source(file.path(repo,"codes","00_diagnostics.R"))
source(file.path(repo,"codes","simulation_helpers.R"))
engine <- mixedgp_simulation_engine(file.path(repo,"codes"))
workers <- suppressWarnings(as.integer(Sys.getenv("EIVGP_WORKERS","1")))
if(length(workers)!=1L || is.na(workers) || workers<1L) stop("EIVGP_WORKERS must be positive.")
data_root <- Sys.getenv("EIVGP_DATA_ROOT",file.path(repo,"reproduction","data","synthetic",study))
cfg <- get(paste0(study,"_simulation_config"))(mode,code_dir=file.path(repo,"codes"),data_root=data_root,workers=workers)
cfg$parallel$workers <- workers
methods <- strsplit(Sys.getenv("EIVGP_COMPETITOR_METHODS","UC-GP,LVGP,EzGP"),",",fixed=TRUE)[[1L]]
methods <- trimws(methods)
if(anyDuplicated(methods) || any(!methods %in% c("UC-GP","LVGP","EzGP"))) stop("Invalid EIVGP_COMPETITOR_METHODS.")
controls <- engine$mixedgp_competitor_protocol(study)
output <- file.path(Sys.getenv("EIVGP_COMPETITOR_REPORT_ROOT",
  file.path(repo,"reproduction","competitor-reports")),study,mode)
overleaf <- Sys.getenv("EIVGP_OVERLEAF_ROOT","")
print(list(study=study,selected_datasets=mode,cells=mixedgp_cell_summary(cfg),
  data_root=data_root,cache=engine$.mixedgp_competitor_cache_root,
  workers=workers,methods=methods,controls=controls,output=output,overleaf=overleaf))
cat("Competitor settings are shared across modes. LVGP's existing time limit is PER ATTEMPT.\n")
if(action=="plan") {
  print(mixedgp_existing_data_status(cfg)); quit(status=0L)
}
if(!nzchar(overleaf) || !dir.exists(overleaf)) stop("Set EIVGP_OVERLEAF_ROOT to the existing paper project folder.")
if(!requireNamespace("knitr",quietly=TRUE)) stop("Install experiment dependencies first (knitr is missing).")
if(action %in% c("run","retry")) {
  preflight <- engine$mixedgp_competitor_preflight(methods,strict=FALSE)
  if(any(!preflight$available)) stop("Install missing comparison packages first: Rscript --vanilla experiments/install_eivgp_dependencies.R")
}
## Validate the complete selected frozen collection before spending time fitting.
for(cell in cfg$cells) mixedgp_verify_cell_data(cfg,cell,engine)
dir.create(output,recursive=TRUE,showWarnings=FALSE)
statuses <- metrics <- predictions <- list()
for(cell in cfg$cells) {
  cat("Competitors:",cell$id,"with",cell$n_rep,"datasets\n")
  work <- function(rep_id) {
    frozen <- mixedgp_read_cell_replication(cfg,cell,rep_id)
    dat <- frozen$data; train <- dat$train; test <- dat$test
    seed <- if(study=="study1") 250000L+1000L*rep_id else
      1000000L*match(cell$scenario,c("primary","latent_additive_control","high_uncertainty","logistic_misspec"))+10000L*rep_id+100L
    X <- if(study=="study1") matrix(train$x,ncol=1L) else train$X
    C <- if(study=="study1") matrix(train$c,ncol=1L) else train$C
    Xt <- if(study=="study1") matrix(test$x,ncol=1L) else test$X
    Ct <- if(study=="study1") matrix(test$c,ncol=1L) else test$C
    cat(cell$id,"replication",rep_id,"starting\n")
    ans <- engine$mixedgp_cached_competitors(X,train$y,C,Xt,Ct,n_draw=500L,seed=seed,
      m_vec=rep(cell$m,ncol(C)),methods=methods,controls=controls,
      allow_fit=action %in% c("run","retry"),retry_failed=action=="retry")
    s <- ans$status; s$cell_id<-cell$id; s$rep<-rep_id
    s$dataset_md5 <- engine$mixedgp_competitor_hash(frozen$data)
    m <- p <- list()
    for(method in names(ans$draws)) {
      m[[method]] <- engine$summarize_predictive_samples_1d(ans$draws[[method]],test$y,
        method=method,rep_id=rep_id,n_calib=NA_integer_,scenario=cell$scenario)
      m[[method]]$cell_id <- cell$id
      p[[method]] <- data.frame(cell_id=cell$id,rep=rep_id,method=method,test_id=seq_along(test$y),
        y=test$y,mean=ans$predictive_means[[method]],variance=ans$predictive_variances[[method]])
    }
    cat(cell$id,"replication",rep_id,"saved; successful methods:",length(m),"\n")
    list(status=s,metrics=mixedgp_bind_rows_base(m),predictions=mixedgp_bind_rows_base(p))
  }
  safe_work <- function(rep_id) tryCatch(work(rep_id), error=function(e)
    structure(conditionMessage(e), class="try-error"))
  answers <- engine$mixedgp_parallel_lapply(as.list(seq_len(cell$n_rep)),safe_work,
    n_cores=workers,seeds=700000L+seq_len(cell$n_rep),mc.preschedule=FALSE)
  for(i in seq_along(answers)) {
    a <- answers[[i]]; id<-paste(cell$id,i)
    if(inherits(a,"try-error") || is.null(a)) {
      statuses[[id]]<-data.frame(cell_id=cell$id,rep=i,method="dataset",status="worker_error",
        message=paste(a,collapse=" ")); next
    }
    statuses[[id]]<-a$status; metrics[[id]]<-a$metrics; predictions[[id]]<-a$predictions
  }
  ## Write incremental tables after each cell. The per-method cache is saved sooner.
  for(item in c("statuses","metrics","predictions"))
    mixedgp_atomic_write_csv(mixedgp_bind_rows_base(get(item)),file.path(output,paste0(item,".csv")))
}
status <- mixedgp_bind_rows_base(statuses); scores <- mixedgp_bind_rows_base(metrics)
summary <- list()
if(nrow(scores)) for(cell_id in unique(scores$cell_id)) for(method in unique(scores$method)) {
  d <- scores[scores$cell_id==cell_id & scores$method==method,,drop=FALSE]
  if(!nrow(d)) next
  for(metric in intersect(c("RMSE","MAE","CRPS","NLPD","Coverage95","Width95","IntervalScore95"),names(d))) {
    v<-d[[metric]]; v<-v[is.finite(v)]
    summary[[paste(cell_id,method,metric)]]<-data.frame(cell_id=cell_id,method=method,metric=metric,
      mean=if(length(v)) mean(v) else NA_real_,mcse=if(length(v)>1L) sd(v)/sqrt(length(v)) else NA_real_,n_success=length(v))
  }
}
summary<-mixedgp_bind_rows_base(summary)
mixedgp_atomic_write_csv(summary,file.path(output,"predictive_summary.csv"))
saveRDS(list(config=cfg,controls=controls,methods=methods,cache=engine$.mixedgp_competitor_cache_root,
  session=sessionInfo()),file.path(output,"provenance.rds"))
target <- file.path(overleaf,"tables","competitors",study,mode)
dir.create(target,recursive=TRUE,showWarnings=FALSE)
for(f in list.files(output,full.names=TRUE)) if(!file.copy(f,target,overwrite=TRUE)) stop("Export failed: ",f)
writeLines(if(nrow(summary)) knitr::kable(summary,format="latex",booktabs=TRUE,digits=3) else
  "% No successful competitor predictions. See statuses.csv; do not use an older table.",
  file.path(target,"predictive_summary.tex"))
cat("Saved competitor-only tables to",target,"\nNo EIV-GP fits were run.\n")
if(any(status$status!="success")) {
  warning("Incomplete comparisons: inspect statuses.csv; successful caches are retained.",call.=FALSE)
  quit(status=1L)
}
