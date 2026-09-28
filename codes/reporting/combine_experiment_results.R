## Reporting only: combine portable numeric outputs, never load/fix model fits.
mixedgp_combine_predictions <- function(own, competitors, status, identity, cells) {
  metrics <- c("RMSE","MAE","CRPS","NLPD","Coverage95","Width95","IntervalScore95")
  key <- function(x) paste(x$cell_id,x$rep,sep="/")
  if(anyDuplicated(key(identity))) stop("Duplicate MCMC dataset identities.")
  expected <- mixedgp_bind_rows_base(lapply(cells,function(c)
    data.frame(cell_id=c$id,rep=seq_len(c$n_rep))))
  for(name in c("own","competitors","status")) {
    x <- get(name)
    if(nrow(x)) x <- x[key(x) %in% key(expected),,drop=FALSE]
    assign(name,x)
  }
  if(nrow(own) && any(!key(own) %in% key(identity)))
    stop("MCMC metrics lack dataset identity; do not infer identity from replication number alone.")
  if(nrow(competitors)) {
    required <- c("cell_id","rep","method","dataset_md5","status")
    if(!all(required %in% names(status))) stop("Competitor statuses lack dataset checksums.")
    sk <- paste(key(status),status$method)
    if(anyDuplicated(sk)) stop("Duplicate competitor statuses.")
    at <- match(paste(key(competitors),competitors$method),sk)
    if(anyNA(at) || any(status$status[at]!="success")) stop("Competitor metrics lack successful status.")
    ik <- match(key(competitors),key(identity))
    shared <- !is.na(ik)
    if(anyNA(status$dataset_md5[at]) || any(!nzchar(status$dataset_md5[at]))) stop("Missing competitor checksum.")
    if(any(status$dataset_md5[at[shared]] != identity$dataset_md5[ik[shared]]))
      stop("Dataset checksum mismatch: results cannot be paired.")
  }
  normalize <- function(x,source) {
    if(!nrow(x)) return(data.frame())
    if(!"evaluation_stratum" %in% names(x)) x$evaluation_stratum <- "overall"
    x <- x[x$evaluation_stratum=="overall",,drop=FALSE]
    x$source <- source
    if(anyDuplicated(paste(key(x),x$method,x$n_calib))) stop("Duplicate method/dataset/calibration metrics.")
    x
  }
  own <- normalize(own,"mcmc_experiment")
  competitors <- normalize(competitors,"competitor_experiment")
  combined <- mixedgp_bind_rows_base(list(own,competitors))
  summary <- paired <- list()
  if(nrow(combined)) {
    groups <- split(combined,paste(combined$cell_id,combined$method,combined$n_calib,sep="/"))
    for(g in groups) for(metric in intersect(metrics,names(g))) {
      v <- g[[metric]]; v<-v[is.finite(v)]
      c <- cells[[match(g$cell_id[1],vapply(cells,`[[`,"","id"))]]
      summary[[length(summary)+1L]] <- data.frame(cell_id=c$id,method=g$method[1],
        n_calib=g$n_calib[1],metric=metric,mean=if(length(v)) mean(v) else NA_real_,
        mcse=if(length(v)>1L) sd(v)/sqrt(length(v)) else NA_real_,
        n_success=length(v),n_expected=c$n_rep)
    }
  }
  if(nrow(own) && nrow(competitors)) {
    main <- own[own$method=="EIV-GP",,drop=FALSE]
    for(method in unique(competitors$method)) {
      comp <- competitors[competitors$method==method,,drop=FALSE]
      pairs <- merge(main,comp,by=c("cell_id","rep"),suffixes=c("_eiv","_competitor"))
      if(!nrow(pairs)) next
      for(g in split(pairs,paste(pairs$cell_id,pairs$n_calib_eiv))) for(metric in metrics) {
        a<-g[[paste0(metric,"_eiv")]]; b<-g[[paste0(metric,"_competitor")]]
        if(is.null(a)||is.null(b)) next
        delta<-a-b; delta<-delta[is.finite(a)&is.finite(b)]
        paired[[length(paired)+1L]]<-data.frame(cell_id=g$cell_id[1],n_calib=g$n_calib_eiv[1],
          competitor=method,metric=metric,mean_eiv_minus_competitor=if(length(delta)) mean(delta) else NA_real_,
          mcse=if(length(delta)>1L) sd(delta)/sqrt(length(delta)) else NA_real_,n_pairs=length(delta))
      }
    }
  }
  availability <- expected
  availability$mcmc_dataset_available <- key(expected) %in% key(identity)
  availability$competitor_metrics_available <- key(expected) %in% key(competitors)
  list(metrics=combined,summary=mixedgp_bind_rows_base(summary),
       paired=mixedgp_bind_rows_base(paired),competitor_status=status,availability=availability)
}

mixedgp_combine_saved_runs <- function(run_dir, competitor_dir, output_dir) {
  config <- readRDS(file.path(run_dir,"config","resolved_config.rds"))
  provenance <- readRDS(file.path(competitor_dir,"provenance.rds"))
  if(!identical(config$study,provenance$config$study)) stop("Study mismatch.")
  read_csv <- function(p) {x<-read.csv(p,check.names=FALSE); x}
  own <- ids <- diagnostics <- list()
  for(cell in config$cells) {
    p <- file.path(run_dir,"cells",cell$id)
    if(!file.exists(file.path(p,"report_inputs.rds"))) next
    if(!file.exists(file.path(p,"dataset_identity.csv")))
      stop("This older run lacks dataset_identity.csv: ",p,". Keep it; identity must be audited before merging.")
    saved <- readRDS(file.path(p,"report_inputs.rds"))
    if(!identical(saved$study,config$study)) stop("Reporting checkpoint study mismatch.")
    x <- saved$context$mc_results
    if(nrow(x)) {x$cell_id<-cell$id; own[[cell$id]]<-x}
    x <- read_csv(file.path(p,"dataset_identity.csv")); x$cell_id<-cell$id; ids[[cell$id]]<-x
    x <- saved$context[[if(config$study=="study1") "mcmc_diagnostics" else "mc_diagnostics"]]
    if(is.data.frame(x)&&nrow(x)) {x$cell_id<-cell$id; diagnostics[[cell$id]]<-x}
  }
  statuses <- read_csv(file.path(competitor_dir,"statuses.csv"))
  methods <- provenance$methods
  if(is.null(methods)) methods <- c("UC-GP","LVGP","EzGP")
  present <- if(nrow(statuses)) paste(statuses$cell_id,statuses$rep,statuses$method) else character()
  absent <- list()
  for(cell in config$cells) for(method in methods) for(rep in seq_len(cell$n_rep)) {
    if(!paste(cell$id,rep,method) %in% present)
      absent[[length(absent)+1L]] <- data.frame(cell_id=cell$id,rep=rep,method=method,
        status="missing",message="No saved competitor status",dataset_md5=NA_character_)
  }
  statuses <- mixedgp_bind_rows_base(c(list(statuses),absent))
  ans <- mixedgp_combine_predictions(mixedgp_bind_rows_base(own),
    read_csv(file.path(competitor_dir,"metrics.csv")),statuses,
    mixedgp_bind_rows_base(ids),config$cells)
  ans$mcmc_diagnostics <- mixedgp_bind_rows_base(diagnostics)
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  for(name in names(ans)) mixedgp_atomic_write_csv(ans[[name]],file.path(output_dir,paste0(name,".csv")))
  writeLines(if(nrow(ans$summary)) knitr::kable(ans$summary,format="latex",booktabs=TRUE,digits=3) else
    "% No metrics available.",file.path(output_dir,"predictive_summary.tex"))
  saveRDS(list(mcmc_config=config,competitor_provenance=provenance,
    mcmc_run=normalizePath(run_dir),competitor_report=normalizePath(competitor_dir),
    created=Sys.time(),session=sessionInfo()),file.path(output_dir,"provenance.rds"))
  invisible(ans)
}
