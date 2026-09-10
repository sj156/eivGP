# Repository-only recovery orchestration. No model fitting occurs when sourced.
mixedgp_recovery_source_cell <- function(path, env) sys.source(path,env,chdir=TRUE)
mixedgp_recovery_model <- function(run, cell) {
  directory <- file.path(run,"cells",cell$id,"results")
  files <- list.files(directory,pattern="^(flagged_fit|fit)_.*rep001_cal(ib)?050.*[.]rds$",
    recursive=TRUE,full.names=TRUE)
  if(length(files)!=1L)stop("Need one saved replication-1/calibration-50 fit to preserve archived priors and controls for ",cell$id,
    ". Supply the full source archive; no current prior defaults were substituted.")
  saved <- readRDS(files[1]);f <- if(is.null(saved$fit))saved else saved$fit
  if(!identical(f$sampler_version,"0.3.1") || is.null(f$priors) || is.null(f$control) ||
     is.null(f$kernel) || !identical(as.integer(f$data$q),as.integer(cell$q)) ||
     !identical(as.integer(f$data$d),2L) || !identical(f$data$ident,"lower_triangular"))
    stop("Saved model is incompatible with publication recovery: ",files[1])
  prior <- f$priors;prior$dictionary_type <- NULL
  list(priors=prior,control=f$control,kernel=f$kernel,source=files[1],
    source_md5=unname(tools::md5sum(files[1])),sampler_version=f$sampler_version)
}

mixedgp_recovery_valid <- function(d) {
  metrics <- c("RMSE","MAE","CRPS","NLPD","Coverage95","Width95","IntervalScore95")
  if (!is.data.frame(d) || !all(c("rep",metrics) %in% names(d))) return(rep(FALSE,nrow(d)))
  ok <- apply(as.matrix(d[,metrics,drop=FALSE]),1,function(x)all(is.finite(x)))
  ok & d$RMSE>=0 & d$CRPS>=0 & d$Coverage95>=0 & d$Coverage95<=1 & d$Width95>0
}

mixedgp_recovery_overall <- function(d) {
  if (!nrow(d)) return(d)
  if ("evaluation_stratum" %in% names(d)) d <- d[is.na(d$evaluation_stratum)|d$evaluation_stratum=="overall",,drop=FALSE]
  d
}

mixedgp_recovery_audit <- function(eiv, competitors, config) {
  rows <- list(); eiv <- mixedgp_recovery_overall(eiv)
  for (cell in config$cells) for (method in c("EIV-GP","UC-GP","LVGP","EzGP")) {
    for (cal in if(method=="EIV-GP")cell$calibration_grid else NA_integer_) {
      d <- if(method=="EIV-GP")eiv else competitors
      if(nrow(d)) d <- d[d$cell_id==cell$id & d$method==method,,drop=FALSE]
      if(nrow(d) && method=="EIV-GP") d <- d[!is.na(d$n_calib)&d$n_calib==cal,,drop=FALSE]
      ids <- if(nrow(d))as.integer(d$rep)else integer()
      valid <- if(nrow(d))ids[mixedgp_recovery_valid(d)]else integer()
      duplicate <- anyDuplicated(ids)>0L
      unexpected <- any(!ids %in% seq_len(cell$n_rep)) ||
        (nrow(d)>0L && any(!is.finite(d$rep) | d$rep != floor(d$rep)))
      missing <- setdiff(seq_len(cell$n_rep),valid)
      rows[[length(rows)+1L]] <- data.frame(study=config$study,cell_id=cell$id,method=method,
        n_calib=cal,expected=cell$n_rep,n_valid=length(unique(valid)),
        missing_ids=paste(missing,collapse=","),duplicate_ids=duplicate,unexpected_ids=unexpected,
        complete=!duplicate&&!unexpected&&!length(missing)&&length(ids)==cell$n_rep)
    }
  }
  do.call(rbind,rows)
}

mixedgp_recovery_upsert <- function(old, fresh, keys) {
  if(!nrow(fresh)) return(old)
  if(!nrow(old)) return(fresh)
  stopifnot(all(keys %in% names(old)),all(keys %in% names(fresh)))
  key <- function(d)do.call(paste,c(lapply(d[,keys,drop=FALSE],as.character),sep="\034"))
  mixedgp_bind_rows_base(list(old[!key(old)%in%key(fresh),,drop=FALSE],fresh))
}

mixedgp_recovery_read_csv <- function(path) {
  if(file.exists(path))read.csv(path,stringsAsFactors=FALSE)else data.frame()
}

mixedgp_recovery_save <- function(state, output) {
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  mixedgp_atomic_save_rds(state,file.path(output,"state.rds"))
  for(study in names(state)) {
    folder <- file.path(output,study);dir.create(folder,recursive=TRUE,showWarnings=FALSE)
    for(name in names(state[[study]]$raw)) {
      d <- state[[study]]$raw[[name]]
      if(is.data.frame(d)&&ncol(d))mixedgp_atomic_write_csv(d,file.path(folder,paste0(name,".csv")))
    }
    mixedgp_atomic_save_rds(state[[study]]$raw,file.path(folder,"all_raw_outputs.rds"))
    for(name in c("metrics","statuses","predictions"))
      if(ncol(state[[study]][[name]]))mixedgp_atomic_write_csv(state[[study]][[name]],file.path(folder,paste0("competitor_",name,".csv")))
  }
}

mixedgp_recovery_tables <- function(state, configs, output, certify=FALSE) {
  audits <- lapply(names(configs),function(study)mixedgp_recovery_audit(
    state[[study]]$raw$predictive_metrics,state[[study]]$metrics,configs[[study]]))
  audit <- do.call(rbind,audits)
  mixedgp_atomic_write_csv(audit,file.path(output,"completeness.csv"))
  table2 <- audit[audit$study=="study2" & (audit$method!="EIV-GP"|(!is.na(audit$n_calib)&audit$n_calib==50)),,drop=FALSE]
  mixedgp_atomic_write_csv(table2,file.path(output,"table2_completeness.csv"))
  summaries <- list()
  for(study in names(configs)) {
    d <- mixedgp_bind_rows_base(list(mixedgp_recovery_overall(state[[study]]$raw$predictive_metrics),state[[study]]$metrics))
    if(!nrow(d))next
    d <- d[d$method %in% c("EIV-GP","UC-GP","LVGP","EzGP"),,drop=FALSE]
    groups <- split(d,paste(d$cell_id,d$method,ifelse(is.na(d$n_calib),"none",d$n_calib)))
    for(z in groups)for(metric in c("RMSE","MAE","CRPS","NLPD","Coverage95","Width95","IntervalScore95")) {
      v <- z[[metric]]; v<-v[is.finite(v)]
      summaries[[length(summaries)+1L]]<-data.frame(study=study,cell_id=z$cell_id[1],method=z$method[1],
        n_calib=z$n_calib[1],metric=metric,n=length(v),mean=if(length(v))mean(v)else NA_real_,
        mcse=if(length(v)>1L)sd(v)/sqrt(length(v))else NA_real_)
    }
  }
  summary <- mixedgp_bind_rows_base(summaries)
  mixedgp_atomic_write_csv(summary,file.path(output,"predictive_summary.csv"))
  destination <- file.path(output,"table2_complete.tex")
  if(isTRUE(certify) && nrow(table2)==12L && all(table2$complete) && all(table2$expected==50L)) {
    lines<-c("\\begin{table}[H]","\\centering\\small",
      "\\caption{Study II predictive comparison at calibration size 50. All four methods use the same 50 datasets per setting. Entries are means (MCSEs); bold marks the lowest displayed RMSE and CRPS, including ties.}",
      "\\label{tab:study2_comparison}","\\begin{tabular}{lrrr}","\\toprule",
      "Method & RMSE & CRPS & Coverage (\\%)\\\\","\\midrule")
    cells<-c(primary_q2="Gaussian, $q=2$",primary_q4_calibration="Gaussian, $q=4$",logistic_q4="Logistic, $q=4$")
    for(cell in names(cells)) {
      lines<-c(lines,paste0("\\multicolumn{4}{l}{\\textit{",cells[[cell]],"}}\\\\"))
      z<-summary[summary$study=="study2"&summary$cell_id==cell&
        (summary$method!="EIV-GP"|(!is.na(summary$n_calib)&summary$n_calib==50)),,drop=FALSE]
      best<-sapply(c("RMSE","CRPS"),function(m)min(round(z$mean[z$metric==m],3)))
      for(method in c("EIV-GP","UC-GP","LVGP","EzGP")) {
        values<-sapply(c("RMSE","CRPS","Coverage95"),function(m){
          row<-z[z$method==method&z$metric==m,,drop=FALSE];stopifnot(nrow(row)==1L,row$n==50L)
          scale<-if(m=="Coverage95")100 else 1; digits<-if(m=="Coverage95")1 else 3
          value<-formatC(row$mean*scale,format="f",digits=digits)
          if(m!="Coverage95"&&round(row$mean,3)==best[[m]])value<-paste0("\\textbf{",value,"}")
          paste0(value," (",formatC(row$mcse*scale,format="f",digits=digits),")")
        })
        lines<-c(lines,paste0(gsub("-","--",method,fixed=TRUE)," & ",paste(values,collapse=" & "),"\\\\"))
      }
      lines<-c(lines,"\\addlinespace")
    }
    writeLines(c(lines,"\\bottomrule","\\end{tabular}","\\end{table}"),destination)
  } else {
    # Invalidate an earlier complete export if later checks no longer pass.
    writeLines("% INCOMPLETE: inspect table2_completeness.csv. No R=50 table is certified.",destination)
  }
  audit
}

mixedgp_recover_publication <- function(repo, archive_root, data_root, output,
    action=c("plan","run","check"), chain_workers=4L, studies=c("study1","study2")) {
  action<-match.arg(action)
  chain_workers<-mixedgp_validate_scalar_integer(chain_workers,"chain_workers",1L)
  repo<-normalizePath(repo,mustWork=TRUE)
  archive_root<-normalizePath(archive_root,mustWork=TRUE)
  if(!all(studies%in%c("study1","study2"))||anyDuplicated(studies))stop("Invalid study selection.")
  dir.create(output,recursive=TRUE,showWarnings=FALSE);output<-normalizePath(output,mustWork=TRUE)
  if(output==archive_root||startsWith(archive_root,paste0(output,"/")))stop("Recovery output must be separate from the source archives.")
  # One recovery coordinator per output. Never delete or steal another run's lock.
  lock<-file.path(output,"recovery.lock")
  if(!dir.create(lock,showWarnings=FALSE))stop("Recovery output is locked: ",lock,". Check the owner; no lock was removed.")
  on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  saveRDS(list(pid=Sys.getpid(),host=Sys.info()[["nodename"]],created=as.character(Sys.time())),file.path(lock,"owner.rds"))
  engine<-mixedgp_simulation_engine(file.path(repo,"codes"))
  configs<-sources<-list();state<-list()
  locate<-function(study){
    roots<-c(archive_root,file.path(archive_root,study),file.path(archive_root,"results",study))
    candidates<-unique(unlist(lapply(roots,function(root)Sys.glob(file.path(root,paste0(study,"-publication*"))))))
    candidates<-candidates[file.exists(file.path(candidates,"run_summary.rds"))]
    if(length(candidates)!=1L)stop("Expected exactly one ",study," publication archive under ",archive_root,"; found ",length(candidates))
    candidates
  }
  for(study in studies){
    run<-locate(study)
    if(output==normalizePath(run)||startsWith(output,paste0(normalizePath(run),"/")))stop("Output must be outside source run directories.")
    archive<-readRDS(file.path(run,"run_summary.rds"));cfg<-archive$config
    if(any(vapply(cfg$cells,function(c)c$n_rep!=50L,logical(1))))stop("Recovery requires the publication 50-dataset design.")
    cfg$code_dir<-file.path(repo,"codes");cfg$data_root<-file.path(data_root,study)
    cfg$use_cache<-TRUE;cfg$parallel$level<-"chains";cfg$parallel$workers<-1L;cfg$parallel$chain_workers<-chain_workers
    configs[[study]]<-cfg
    report<-file.path(archive_root,"competitor-reports",study,"publication")
    if(!file.exists(file.path(report,"metrics.csv")))stop("Missing original competitor report: ",report)
    sources[[study]]<-list(run=run,report=report,input_manifest=archive$data_manifest)
    state[[study]]<-list(raw=readRDS(file.path(run,"combined","all_raw_outputs.rds")),
      metrics=mixedgp_recovery_read_csv(file.path(report,"metrics.csv")),
      statuses=mixedgp_recovery_read_csv(file.path(report,"statuses.csv")),
      predictions=mixedgp_recovery_read_csv(file.path(report,"predictions.csv")))
  }
  source_files<-unlist(lapply(sources,function(z)c(file.path(z$run,"run_summary.rds"),file.path(z$run,"combined","all_raw_outputs.rds"),file.path(z$report,c("metrics.csv","statuses.csv","predictions.csv")))))
  code_files<-file.path(repo,"codes",c(mixedgp_simulation_modules(),"simulation_helpers.R","publication_recovery.R","02_study2_monte_carlo.R","setup_study2_experiment.R"))
  identity<-list(schema="publication-recovery-v1",source_md5=tools::md5sum(source_files),
    code_md5=tools::md5sum(code_files),
    R_version=R.version.string, package_versions=vapply(c("kergp","LVGP","EzGP","TruncatedNormal"),
      function(p)if(requireNamespace(p,quietly=TRUE))as.character(utils::packageVersion(p))else"missing",character(1)),data_root=normalizePath(data_root,mustWork=FALSE),studies=studies)
  provenance_file<-file.path(output,"provenance.rds")
  if(file.exists(provenance_file)) {
    previous<-readRDS(provenance_file)
    if(!identical(previous$identity,identity))stop("Recovery sources/code changed. Use a new output folder; completed results were not overwritten.")
  } else mixedgp_atomic_save_rds(list(identity=identity,configs=configs,sources=sources,
    rescue_protocol=lapply(studies,engine$mixedgp_competitor_protocol),session=sessionInfo()),provenance_file)
  if(file.exists(file.path(output,"state.rds")))state<-readRDS(file.path(output,"state.rds"))
  audit<-mixedgp_recovery_tables(state,configs,output)
  print(audit[,c("study","cell_id","method","n_calib","n_valid","missing_ids","complete")],row.names=FALSE)
  if(action=="plan") {
    message("Plan only: no fitting. Frozen inputs expected under ",data_root,
      ". Use 'run' with the same arguments; 'check' verifies outputs.")
    return(invisible(audit))
  }
  # Validate original frozen files, including agreement with the archived run's
  # manifest. Never regenerate datasets or silently substitute new random seeds.
  for(study in studies)for(cell in configs[[study]]$cells){
    cfg<-configs[[study]];mixedgp_verify_cell_data(cfg,cell,engine)
    archived<-sources[[study]]$input_manifest
    current<-readRDS(file.path(mixedgp_data_cell_directory(cfg,cell),"manifest.rds"))
    a<-archived[archived$cell_id==cell$id,,drop=FALSE]
    current<-current[match(a$file,current$file),,drop=FALSE]
    if(nrow(a)!=50L||anyNA(current$md5)||!identical(tolower(as.character(a$md5)),tolower(as.character(current$md5))))
      stop("Frozen datasets differ from the archived 50-dataset collection: ",cell$id,
        ". Point --data-root to the original Linux inputs; do not regenerate them.")
    # CSV decimals/R serialization can vary by platform; validate actual held-out
    # responses before reusing the published competitor metrics.
    preds<-state[[study]]$predictions
    for(rep_id in seq_len(cell$n_rep)){
      d<-mixedgp_read_cell_replication(cfg,cell,rep_id)
      for(method in c("UC-GP","LVGP","EzGP")){
        old<-preds[preds$cell_id==cell$id&preds$rep==rep_id&preds$method==method,,drop=FALSE]
        met<-state[[study]]$metrics
        met<-met[met$cell_id==cell$id&met$rep==rep_id&met$method==method,,drop=FALSE]
        if(nrow(met)&&any(mixedgp_recovery_valid(met))){
          if(nrow(old)!=cell$n_test||anyDuplicated(old$test_id)||!setequal(old$test_id,seq_len(cell$n_test))||
             !isTRUE(all.equal(as.numeric(old$y[order(old$test_id)]),as.numeric(d$data$test$y),tolerance=1e-12)))
            stop("Stored competitor predictions do not match frozen test responses: ",study," ",cell$id," ",rep_id," ",method)
        }
      }
    }
  }
  if(action=="check")return(invisible(mixedgp_recovery_tables(state,configs,output,certify=TRUE)))
  preflight<-engine$mixedgp_competitor_preflight(c("UC-GP","LVGP","EzGP"),strict=FALSE)
  if(any(!preflight$available))stop("Install missing experiment dependencies before recovery: ",paste(preflight$method[!preflight$available],collapse=", "))
  attempts<-mixedgp_recovery_read_csv(file.path(output,"attempts.csv"))
  record<-function(study,cell,rep,method,status,message=""){
    attempts<<-mixedgp_bind_rows_base(list(attempts,data.frame(study=study,cell_id=cell,rep=rep,method=method,status=status,message=message,time=as.character(Sys.time()))))
    mixedgp_atomic_write_csv(attempts,file.path(output,"attempts.csv"))
  }
  for(study in studies)for(cell in configs[[study]]$cells){
    cfg<-configs[[study]]
    audit<-mixedgp_recovery_audit(state[[study]]$raw$predictive_metrics,state[[study]]$metrics,cfg)
    if(any(audit$duplicate_ids|audit$unexpected_ids))stop("Duplicate or unexpected replication IDs; inspect completeness.csv before fitting.")
    missing<-unique(unlist(strsplit(audit$missing_ids[audit$cell_id==cell$id&audit$method=="EIV-GP"],",",fixed=TRUE)))
    missing<-as.integer(missing[nzchar(missing)])
    if(length(missing)&&study!="study2")stop("Unexpected missing Study I EIV-GP outputs; this recovery targets the audited Study II failures.")
    model <- if(length(missing))mixedgp_recovery_model(sources[[study]]$run,cell)else NULL
    if(!is.null(model))mixedgp_atomic_save_rds(model,file.path(output,paste0("model_",cell$id,".rds")))
    for(rep_id in missing){
      message("Recovering EIV-GP ",cell$id," dataset ",rep_id," (original calibration grid)")
      directory<-file.path(output,"fits",study,cell$id,sprintf("rep%03d",rep_id))
      env<-new.env(parent=engine);list2env(mixedgp_cell_controls_study2(cfg,cell,directory),env)
      env$STUDY2_REP_IDS<-rep_id;env$STUDY2_CHECKPOINT_FITS<-TRUE;env$STUDY2_ORACLE_MEAN_METHOD<-"quadrature"
      env$STUDY2_RECOVERY_MODEL<-model
      err<-tryCatch({mixedgp_recovery_source_cell(file.path(repo,"codes","02_study2_monte_carlo.R"),env);NULL},error=function(e)conditionMessage(e))
      fresh<-get0("raw_outputs",envir=env,inherits=FALSE)
      if(!is.null(fresh))for(name in names(fresh))if(is.data.frame(fresh[[name]])&&nrow(fresh[[name]])&&"rep"%in%names(fresh[[name]])){
        z<-mixedgp_tag_output(fresh[[name]],cell,study)
        old<-state[[study]]$raw[[name]];if(is.null(old))old<-data.frame()
        state[[study]]$raw[[name]]<-mixedgp_recovery_upsert(old,z,c("cell_id","rep"))
      }
      # Salvage predictive checkpoints even if a later scientific task failed.
      checkpoint_files<-list.files(directory,pattern="^predictive_checkpoint_.*\\.rds$",recursive=TRUE,full.names=TRUE)
      for(file in checkpoint_files){
        z<-readRDS(file);if(!identical(z$identity$cache_spec,env$CACHE_SPEC)||z$identity$rep!=rep_id)stop("Predictive checkpoint provenance mismatch: ",file)
        d<-mixedgp_tag_output(z$metrics,cell,study)
        state[[study]]$raw$predictive_metrics<-mixedgp_recovery_upsert(state[[study]]$raw$predictive_metrics,d,
          c("cell_id","rep","method","n_calib","evaluation_stratum"))
      }
      record(study,cell$id,rep_id,"EIV-GP",if(is.null(err))"completed"else"downstream_or_fit_error",if(is.null(err))""else err)
      mixedgp_recovery_save(state,output);mixedgp_recovery_tables(state,configs,output,certify=TRUE)
      rm(env);gc()
    }
    for(rep_id in seq_len(cell$n_rep))for(method in c("UC-GP","LVGP","EzGP")){
      old<-state[[study]]$metrics;old<-old[old$cell_id==cell$id&old$rep==rep_id&old$method==method,,drop=FALSE]
      if(nrow(old)==1L&&mixedgp_recovery_valid(old))next
      frozen<-mixedgp_read_cell_replication(cfg,cell,rep_id);train<-frozen$data$train;test<-frozen$data$test
      X<-if(study=="study1")matrix(train$x,ncol=1)else train$X;C<-if(study=="study1")matrix(train$c,ncol=1)else train$C
      Xt<-if(study=="study1")matrix(test$x,ncol=1)else test$X;Ct<-if(study=="study1")matrix(test$c,ncol=1)else test$C
      seed<-if(study=="study1")250000L+1000L*rep_id else 1000000L*match(cell$scenario,c("primary","latent_additive_control","high_uncertainty","logistic_misspec"))+10000L*rep_id+100L
      ans<-engine$mixedgp_cached_competitors(X,train$y,C,Xt,Ct,n_draw=500L,seed=seed,
        m_vec=rep(cell$m,ncol(C)),methods=method,controls=engine$mixedgp_competitor_protocol(study),
        cache_root=file.path(output,"competitor-cache"),allow_fit=TRUE,retry_failed=TRUE)
      status<-ans$status;status$cell_id<-cell$id;status$rep<-rep_id;status$dataset_md5<-engine$mixedgp_competitor_hash(frozen$data)
      state[[study]]$statuses<-mixedgp_recovery_upsert(state[[study]]$statuses,status,c("cell_id","rep","method"))
      if(all(status$status=="success")){
        met<-engine$summarize_predictive_samples_1d(ans$draws[[method]],test$y,method,rep_id,NA_integer_,cell$scenario)
        met$cell_id<-cell$id
        state[[study]]$metrics<-mixedgp_recovery_upsert(state[[study]]$metrics,met,c("cell_id","rep","method"))
        pred<-data.frame(cell_id=cell$id,rep=rep_id,method=method,test_id=seq_along(test$y),y=test$y,
          mean=ans$predictive_means[[method]],variance=ans$predictive_variances[[method]])
        state[[study]]$predictions<-mixedgp_recovery_upsert(state[[study]]$predictions,pred,c("cell_id","rep","method"))
      }
      record(study,cell$id,rep_id,method,status$status[1],status$message[1])
      mixedgp_recovery_save(state,output);mixedgp_recovery_tables(state,configs,output,certify=TRUE)
    }
  }
  audit<-mixedgp_recovery_tables(state,configs,output,certify=TRUE)
  message(if(all(audit$complete))"All requested predictive comparisons have 50 valid original datasets." else
    "Recovery remains incomplete. Inspect completeness.csv and attempts.csv; no missing scores were invented.")
  invisible(audit)
}
