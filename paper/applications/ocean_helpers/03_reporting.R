# 03: Combine compatible completed tasks; summarize and plot prediction scores.

# Capture console output, messages and warnings in a per-run log.
# Keep HTML focused on analysis while preserving complete execution logs.
ocean_logged <- function(path, expr) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  connection <- file(path, open = "wt")
  output_depth <- sink.number()
  sink(connection)
  on.exit({
    while (sink.number() > output_depth) sink()
    close(connection)
  }, add = TRUE)
  warnings <- character()
  value <- withCallingHandlers(eval(substitute(expr), envir = parent.frame()),
    message = function(m) { cat(conditionMessage(m)); invokeRestart("muffleMessage") },
    warning = function(w) {
      warnings <<- unique(c(warnings, conditionMessage(w)))
      cat("WARNING:", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    })
  attr(value, "notebook_warnings") <- warnings
  value
}


# Import only completed, compatible task payloads; validate before any copying.
# Independent machine roots are preferable to concurrent writers on shared storage.
ocean_combine_results <- function(source_dirs, output_dir, identity) {
  sources <- as.character(unlist(source_dirs)); sources <- sources[nzchar(sources)]
  plan <- list(); receipts <- list()
  for(source in sources) {
    source <- normalizePath(source,mustWork=TRUE)
    if(identical(source,normalizePath(output_dir,mustWork=TRUE))) stop("Cannot import a directory into itself.")
    path <- file.path(source,"workflow_identity.rds")
    if(!file.exists(path) || !identical(readRDS(path),identity)) stop("Incompatible imported run: ",source)
    done <- list.files(file.path(source,"workflow_tasks"),"[.]done[.]rds$",full.names=TRUE)
    for(receipt_path in done) {
      receipt <- readRDS(receipt_path); id <- receipt$task$task_id
      ocean_verify_receipt(source,id,identity)
      for(rel in receipt$files) {
        from <- file.path(source,rel); to <- file.path(output_dir,rel)
        key <- rel
        hash <- ocean_hashes(from)
        if(file.exists(to) && !identical(ocean_hashes(to),hash)) stop("Conflicting existing payload: ",rel)
        if(!is.null(plan[[key]]) && !identical(plan[[key]]$hash,hash)) stop("Conflicting imported payloads: ",rel)
        plan[[key]] <- list(from=from,to=to,hash=hash)
      }
      if(!is.null(receipts[[id]]) && !identical(receipts[[id]]$hashes,receipt$hashes))
        stop("Conflicting completed task: ",id)
      receipts[[id]] <- receipt
    }
  }
  for(item in plan) {
    dir.create(dirname(item$to),recursive=TRUE,showWarnings=FALSE)
    if(!file.exists(item$to)) stopifnot(file.copy(item$from,item$to))
    stopifnot(identical(ocean_hashes(item$to),item$hash))
  }
  for(id in names(receipts)) ocean_atomic_rds(receipts[[id]],
    file.path(output_dir,"workflow_tasks",paste0(id,".done.rds")))
  invisible(names(receipts))
}

# Read saved scores without invoking a fit.
ocean_bind <- function(parts) {
  if(!length(parts)) return(data.frame())
  cols <- unique(unlist(lapply(parts,names)))
  do.call(rbind,lapply(parts,function(x) {for(n in setdiff(cols,names(x))) x[[n]]<-NA; x[,cols,drop=FALSE]}))
}
ocean_prediction_report <- function(output_dir,inputs,profile="paper",report_dir=file.path(output_dir,"prediction_report")) {
  rows <- list(); inventory <- list(); original <- if(profile!="quick") "full5_class4" else "smoke_class4"
  final_sources <- ocean_final_eiv_sources(output_dir,profile)
  add <- function(path,f,cell,variant) {
    if(!file.exists(path)) return(invisible(NULL))
    z <- read.csv(path,stringsAsFactors=FALSE)
    stopifnot(nrow(z)>0L,all(z$fold==f),all(z$n_test==40L))
    if(variant=="EIV") {
      z$prediction_mode <- ifelse(z$method=="EIV-GP xc","class_only","mixed")
      if(cell=="O0") z <- z[z$prediction_mode=="class_only",,drop=FALSE]
      z$method <- "EIV-GP"
    } else if(variant=="Original baseline") z$prediction_mode <- "class_only"
    z$variant <- variant
    z$source <- gsub("\\\\","/",path)
    rows[[length(rows)+1L]] <<- z
  }
  for(f in if(profile!="quick") 1:5 else 1L) {
    initial <- file.path(output_dir,original,sprintf("fold%02d",f))
    add(file.path(initial,"baseline_metrics.csv"),f,"baseline","Original baseline")
    for(cell in if(profile!="quick") c("O0","O20","O50") else c("O0","O20","O50")) {
      final<-final_sources[final_sources$fold==f & final_sources$cell==cell,]
      path<-file.path(output_dir,final$source)
      available<-isTRUE(final$available)
      inventory[[length(inventory)+1L]]<-data.frame(fold=f,cell,
        stage=if(final$continued)ocean_continuation_label(output_dir) else "initial",available,
        initial_screen_pass=final$initial_screen_pass,final_screen_pass=final$final_screen_pass)
      # A cell that needs continuation never falls back to its initial fit.
      if(available)add(file.path(path,"predictive_metrics.csv"),f,cell,"EIV")
    }
  }
  raw <- ocean_bind(rows); summaries <- list()
  metrics <- c("RMSE","MAE","CRPS","NLPD","Coverage95","Width95","IntervalScore95")
  if(nrow(raw)) {
    key <- paste(raw$variant,raw$method,raw$cell,raw$prediction_mode,sep="/")
    if(anyDuplicated(paste(key,raw$fold))) stop("Duplicate fold/task scores in report.")
    for(group in split(seq_len(nrow(raw)),key)) {
      z <- raw[group,,drop=FALSE]
      ans <- z[1,c("variant","method","cell","prediction_mode"),drop=FALSE]
      ans$n_folds <- nrow(z); ans$complete_5fold <- nrow(z)==5L && setequal(z$fold,1:5)
      for(metric in metrics) {
        stopifnot(all(is.finite(z[[metric]])))
        ans[[paste0(metric,"_mean")]] <- mean(z[[metric]])
        ans[[paste0(metric,"_sd")]] <- if(nrow(z)>1L) sd(z[[metric]]) else NA_real_
      }
      summaries[[length(summaries)+1L]] <- ans
    }
  }
  summary_raw <- ocean_bind(summaries); scaled <- raw
  if(nrow(scaled)) {
    for(m in c("RMSE","MAE","CRPS","Width95","IntervalScore95")) scaled[[m]] <- scaled[[m]]/inputs$response_sd
    scaled$NLPD <- scaled$NLPD-log(inputs$response_sd)
  }
  plots <- list()
  report <- list(fold_metrics=raw,summary=summary_raw,availability=ocean_bind(inventory),
    SD_Y=inputs$response_sd,plots=plots,fold_metrics_SD_Y=scaled)
  if(!is.null(report_dir)) ocean_export_prediction_report(report,report_dir)
  report

}
# Shared table formatting for CSV/Markdown/LaTeX/HTML. Formatting does not
# change the fold scores or their aggregation; paper coverage is a percentage.
ocean_mean_sd <- function(mean,sd,digits=2L) {
  ifelse(is.na(sd),sprintf(paste0("%.",digits,"f (SD unavailable)"),mean),
    sprintf(paste0("%.",digits,"f (%.",digits,"f)"),mean,sd))
}
ocean_table_align <- function(d) {
  a <- rep("r",ncol(d)); k <- if("Method" %in% names(d)) 3L else if("Class" %in% names(d)) 2L else 1L
  a[seq_len(min(k,ncol(d)))] <- "l"; a
}
ocean_tex_note <- function(note) gsub("%","\\%",note,fixed=TRUE)
# Rank unrounded means within comparable groups. Incomplete groups are not ranked.
ocean_best_mask <- function(values, groups, eligible = rep(TRUE, nrow(values))) {
  best <- matrix(FALSE, nrow(values), ncol(values), dimnames=list(NULL,names(values)))
  for(rows in split(seq_len(nrow(values)), groups)) {
    if(length(rows)<2L || !isTRUE(all(eligible[rows]))) next
    for(j in seq_len(ncol(values))) {
      v <- values[rows,j]
      if(all(is.finite(v)))
        best[rows,j] <- abs(v-min(v)) <= 1e-12*max(1,abs(min(v)))
    }
  }
  best
}
ocean_html_escape <- function(x) {
  x <- gsub("&","&amp;",as.character(x),fixed=TRUE)
  x <- gsub("<","&lt;",x,fixed=TRUE)
  gsub(">","&gt;",x,fixed=TRUE)
}
# Formatting attributes never enter CSV scores; only the mean is emphasized.
ocean_table_markup <- function(d, format=c("html","pipe","latex")) {
  format <- match.arg(format)
  mask <- attr(d,"best")
  out <- d
  for(j in seq_len(ncol(d))) {
    out[[j]] <- if(format=="html") ocean_html_escape(d[[j]]) else as.character(d[[j]])
    if(!is.null(mask) && names(d)[j] %in% colnames(mask)) {
      rows <- which(mask[,names(d)[j]])
      for(i in rows) {
        value <- out[[j]][i]
        mean <- sub(" \\(.*$","",value)
        suffix <- substring(value,nchar(mean)+1L)
        marked <- switch(format,html=paste0("<strong>",mean,"</strong>"),
          pipe=paste0("**",mean,"**"),latex=paste0("\\textbf{",mean,"}"))
        out[[j]][i] <- paste0(marked,suffix)
      }
    }
  }
  out
}
# Preserve semantic highlighting when presentation columns are omitted.
ocean_table_columns <- function(d,columns) {
  mask<-attr(d,"best")
  out<-d[,columns,drop=FALSE]
  if(!is.null(mask)) attr(out,"best")<-mask[,intersect(columns,colnames(mask)),drop=FALSE]
  out
}
ocean_table_html <- function(d,caption,note="") {
  table <- as.character(knitr::kable(ocean_table_markup(d,"html"),format="html",caption=caption,
    escape=FALSE,row.names=FALSE,align=ocean_table_align(d),
    table.attr='class="ocean-summary table" style="width:100%;border-collapse:collapse"'))
  table <- sub("<caption>",'<caption style="text-align:left;padding:8px 0;color:#555">',table,fixed=TRUE)
  table <- sub("<thead>",'<thead style="background:#f2f4f6">',table,fixed=TRUE)
  table <- gsub('style="text-align:',
    'style="padding:8px 10px;white-space:nowrap;border-bottom:1px solid #ddd;text-align:',
    table,fixed=TRUE)
  paste0('<div class="ocean-table-scroll" style="max-width:100%;overflow-x:auto">',table,
    '</div><p class="note">',ocean_html_escape(note),"</p>")
}
# Training and test information define each comparison panel, as in the paper.
ocean_prediction_panel_html <- function(d, caption, note="") {
  if(!nrow(d)) return("")
  columns <- setdiff(names(d),c("Training exact U","Test oxygen"))
  styled <- ocean_table_markup(d,"html")
  training_groups <- unique(d[["Training exact U"]])
  out <- c("<table class='ocean-prediction table' style='width:100%;border-collapse:collapse'>",
    paste0("<caption style='text-align:left'>",ocean_html_escape(caption),"</caption>"),
    "<thead><tr>",paste0("<th style='text-align:",ifelse(columns=="Method","left","right"),
      ";padding:7px 10px'>",ocean_html_escape(columns),"</th>"),"</tr></thead><tbody>")
  for(j in seq_along(training_groups)) {
    training_rows <- which(d[["Training exact U"]]==training_groups[j])
    tests <- unique(d[["Test oxygen"]][training_rows])
    panel <- match(training_groups[j],c("0%","20%","50%"))
    if(is.na(panel))panel<-j
    title <- paste0(LETTERS[panel],". Training calibration: ",training_groups[j])

    out <- c(out,paste0("<tr class='comparison-panel'><th scope='rowgroup' colspan='",length(columns),
      "' style='text-align:left;background:#f2f4f6;border-top:2px solid #bbb;padding:9px 10px'>",
      ocean_html_escape(title),"</th></tr>"))
    for(test in tests) {
      rows <- training_rows[d[["Test oxygen"]][training_rows]==test]
      out <- c(out,paste0(
        "<tr class='comparison-condition'><th scope='rowgroup' colspan='",length(columns),
        "' style='text-align:left;font-weight:500;padding:7px 10px;border-bottom:1px solid #ddd'>",
        "Test inputs: ",ocean_html_escape(test),"</th></tr>"))
      for(i in rows) out <- c(out,"<tr>",paste0("<td style='text-align:",
        ifelse(columns=="Method","left","right"),";padding:6px 10px'>",
        unlist(styled[i,columns,drop=FALSE],use.names=FALSE),"</td>"),"</tr>")
    }
  }
  paste(c(out,"</tbody></table>",paste0("<p class='note'>",ocean_html_escape(note),"</p>")),collapse="\n")
}
ocean_export_paper_table <- function(d,name,caption,dest,note="") {
  if(!nrow(d)) return(invisible(NULL))
  ocean_atomic_csv(d,file.path(dest,paste0(name,".csv")))
  ocean_atomic_lines(c(as.character(knitr::kable(ocean_table_markup(d,"pipe"),format="pipe",row.names=FALSE,align=ocean_table_align(d))),"",note),file.path(dest,paste0(name,".md")))
  tex_data <- ocean_table_markup(d,"latex")
  tex_data[] <- lapply(tex_data,ocean_tex_note)
  names(tex_data) <- ocean_tex_note(names(tex_data))
  tex <- as.character(knitr::kable(tex_data,format="latex",booktabs=TRUE,caption=caption,label=name,
    escape=FALSE,row.names=FALSE,align=ocean_table_align(d)))
  tex <- paste0(strsplit(tex,"\\end{table}",fixed=TRUE)[[1]][1],
    "\\par\\smallskip\n\\begin{minipage}{\\linewidth}\\footnotesize\n",ocean_tex_note(note),
    "\n\\end{minipage}\n\\end{table}")
  ocean_atomic_lines(tex,file.path(dest,paste0(name,".tex")))
  html <- ocean_table_html(d,caption,note)
  ocean_atomic_lines(c("<!doctype html><meta charset='utf-8'>",
    "<style>body{font:15px/1.5 sans-serif;margin:24px}table{border-collapse:collapse}td,th{padding:7px 10px;border-bottom:1px solid #ddd}caption{text-align:left;font-weight:600}.note{color:#555}</style>",html),
    file.path(dest,paste0(name,".html")))
  invisible(html)
}
# Group training/test information in panel headings, keeping the paper LaTeX
# body compact. CSV retains explicit columns for filtering and checking.
ocean_export_prediction_report <- function(report,dest) {
  ocean_atomic_csv(report$fold_metrics,file.path(dest,"fold_metrics_raw.csv"))
  ocean_atomic_csv(report$summary,file.path(dest,"summary_raw.csv"))
  ocean_atomic_csv(report$fold_metrics_SD_Y,file.path(dest,"fold_metrics_SD_Y.csv"))
  ocean_atomic_csv(report$availability,file.path(dest,"eiv_availability.csv"))
  for(name in names(report$plots)) for(ext in c("png","pdf")) ggplot2::ggsave(
    file.path(dest,paste0(name,".",ext)),report$plots[[name]],width=10,height=6,dpi=180,bg="white")
  invisible(report)
}
# One row per fold combines fit availability, continuation and final screening.
ocean_fit_overview <- function(analysis, output_dir) {
  settings <- analysis$settings
  status <- ocean_final_eiv_sources(output_dir, settings$profile)
  overview <- data.frame(Fold = settings$folds, check.names = FALSE)
  for(cell in settings$cells) {
    overview[[paste0(sub("O", "", cell), "% exact U")]] <- vapply(settings$folds, function(fold) {
      row <- status[status$fold == fold & status$cell == cell, , drop = FALSE]
      if(settings$action == "plan") return("Planned")
      if(!nrow(row) || !isTRUE(row$available))
        return(if(nrow(row) && isTRUE(row$continued)) "Continuation pending" else "Pending")
      if(settings$profile == "quick") return("Initial saved")
      stage <- if(isTRUE(row$continued)) ocean_continuation_label(output_dir) else "Initial"
      paste(stage, if(isTRUE(row$final_screen_pass)) "finite screen met" else "review required", sep = "; ")
    }, "")
  }
  cat(as.character(knitr::kable(overview, format = "pipe", row.names = FALSE)), sep = "\n")
  cat("\n\n")
  if(settings$profile != "quick" && settings$action != "plan") {
    cat("Initial / continuation identifies the saved fit used in the report. Finite screen met means maximum finite Rhat <= 1.01 and minimum finite bulk ESS >= 400; undefined diagnostics and tail ESS are reviewed separately. This screen does not establish convergence of every posterior quantity.\n\n")
    available <- status$available & !is.na(status$final_screen_pass)
    if(any(available & !status$final_screen_pass, na.rm = TRUE))
      cat("Fits marked review remain outside this limited screen.\n\n")
  }
  invisible(overview)
}

# Response-level evidence reconstructed from saved predictions; never fits.
ocean_prediction_evidence <- function(output_dir, inputs, report, profile="paper",
  dest=file.path(output_dir,"prediction_report","evidence")) {
  if(is.null(report) || !nrow(report$fold_metrics)) return(NULL)
  raw <- report$fold_metrics
  sources <- ocean_final_eiv_sources(output_dir,profile)
  root <- if(profile!="quick") "full5_class4" else "smoke_class4"
  cache <- new.env(parent=emptyenv())
  load_saved <- function(path) {
    if(!exists(path,cache,inherits=FALSE)) assign(path,readRDS(path),cache)
    get(path,cache,inherits=FALSE)
  }
  rows <- list()
  for(i in seq_len(nrow(raw))) {
    r <- raw[i,]; cell <- if(r$variant=="Original baseline") "O0" else r$cell
    d <- inputs$folds[[paste(r$fold,if(cell=="O0") "O20" else cell,sep="/")]]; te <- d$test
    known <- if(cell=="O0") rep(FALSE,nrow(te)) else d$observed_test_mixed
    if(r$variant=="Original baseline") {
      saved <- load_saved(file.path(output_dir,root,sprintf("fold%02d",r$fold),"baseline_package_outputs.rds"))
      draws <- saved$draws[[r$method]]; point <- saved$means[[r$method]]
    } else if(r$variant=="EIV") {
      s <- sources[sources$fold==r$fold & sources$cell==cell,]
      path <- file.path(output_dir,s$source)
      saved <- load_saved(file.path(path,"prediction_distributions.rds"))
      stopifnot(identical(as.logical(saved$known),as.logical(known)))
      draws <- saved[[if(r$prediction_mode=="class_only") "xc" else "mixed"]]
      pred <- read.csv(file.path(path,if(r$prediction_mode=="class_only") "predictions.csv" else "mixed_predictions.csv"))
      stopifnot(identical(as.character(pred$cohort_id),as.character(te$cohort_id)))
      point <- pred$predicted_Y
    } else stop("Unsupported prediction variant: ",r$variant)
    stopifnot(ncol(draws)==nrow(te),length(point)==nrow(te),all(is.finite(point)),all(is.finite(draws)))
    se <- (point-te$Y)^2; crps <- eivGP:::crps_sample_matrix(draws,te$Y)
    stopifnot(abs(sqrt(mean(se))-r$RMSE)<1e-8,abs(mean(crps)-r$CRPS)<1e-8)
    rows[[i]] <- data.frame(fold=r$fold,cell,variant=r$variant,method=r$method,
      prediction_mode=r$prediction_mode,cohort_id=te$cohort_id,observed_Y=te$Y,
      predicted_Y=as.numeric(point),mixed_U_observed=known,
      test_U_observed=r$prediction_mode=="mixed" & known,squared_error=se,CRPS=crps)
  }
  oof <- ocean_bind(rows)
  keys <- paste(oof$cell,oof$variant,oof$method,oof$prediction_mode,sep="/")
  pooled <- lapply(split(seq_len(nrow(oof)),keys),function(ids) {
    z <- oof[ids,]; stopifnot(!anyDuplicated(z$cohort_id))
    fold_rmse <- vapply(split(z$squared_error,z$fold),function(v)sqrt(mean(v)),0)
    data.frame(z[1,c("cell","variant","method","prediction_mode")],
      n_folds=length(fold_rmse),N=nrow(z),mean_fold_RMSE=mean(fold_rmse),
      pooled_RMSE=sqrt(mean(z$squared_error)),pooled_CRPS=mean(z$CRPS))
  })
  pooled <- ocean_bind(pooled)
  subgroups <- pairs <- info <- list()
  for(ids in split(which(oof$prediction_mode=="mixed"),keys[oof$prediction_mode=="mixed"])) {
    z <- oof[ids,]; common <- z[1,c("cell","variant","method")]
    ordinal <- oof[oof$cell==z$cell[1] & oof$variant==z$variant[1] &
      oof$method==z$method[1] & oof$prediction_mode=="class_only",]
    j <- match(paste(z$fold,z$cohort_id),paste(ordinal$fold,ordinal$cohort_id))
    stopifnot(!anyNA(j)); ordinal <- ordinal[j,]
    stopifnot(max(abs(z$predicted_Y[!z$mixed_U_observed]-ordinal$predicted_Y[!z$mixed_U_observed]))<1e-8,
      max(abs(z$CRPS[!z$mixed_U_observed]-ordinal$CRPS[!z$mixed_U_observed]))<1e-8)
    for(g in c("Exact U","Ordinal U","All rows")) {
      take <- if(g=="Exact U") z$mixed_U_observed else if(g=="Ordinal U") !z$mixed_U_observed else rep(TRUE,nrow(z))
      subgroups[[length(subgroups)+1L]] <- data.frame(common,group=g,N=sum(take),
        n_folds=length(unique(z$fold[take])),RMSE=sqrt(mean(z$squared_error[take])),CRPS=mean(z$CRPS[take]))
    }
    k <- z$mixed_U_observed
    info[[length(info)+1L]] <- data.frame(common,N=sum(k),n_folds=length(unique(z$fold[k])),
      ordinal_RMSE=sqrt(mean(ordinal$squared_error[k])),exact_RMSE=sqrt(mean(z$squared_error[k])),
      delta_squared_error=mean(z$squared_error[k]-ordinal$squared_error[k]),
      delta_CRPS=mean(z$CRPS[k]-ordinal$CRPS[k]),
      mean_abs_prediction_change=mean(abs(z$predicted_Y[k]-ordinal$predicted_Y[k])))
  }
  for(ids in split(which(oof$variant=="EIV"),keys[oof$variant=="EIV"])) {
    z <- oof[ids,]
    b <- oof[oof$cell==z$cell[1] & oof$prediction_mode==z$prediction_mode[1] & oof$variant!="EIV",]
    for(method in unique(b$method)) {
      bb <- b[b$method==method,]; j<-match(paste(z$fold,z$cohort_id),paste(bb$fold,bb$cohort_id))
      take <- !is.na(j); if(!any(take)) next
      zz <- z[take,]; bb <- bb[j[take],]
      pairs[[length(pairs)+1L]] <- data.frame(cell=z$cell[1],prediction_mode=z$prediction_mode[1],
        comparator=method,N=nrow(zz),n_folds=length(unique(zz$fold)),
        delta_squared_error=mean(zz$squared_error-bb$squared_error),delta_CRPS=mean(zz$CRPS-bb$CRPS))
    }
  }
  ans <- list(oof=oof,pooled=pooled,subgroups=ocean_bind(subgroups),
    information=ocean_bind(info),paired=ocean_bind(pairs),SD_Y=report$SD_Y)
  for(n in c("oof","pooled","subgroups","information","paired"))
    if(nrow(ans[[n]])) ocean_atomic_csv(ans[[n]],file.path(dest,paste0(n,".csv")))
  ans
}

ocean_evidence_table <- function(evidence, type, eiv_only=FALSE, labels=ocean_prediction_labels()) {
  z <- evidence[[type]]
  if(!nrow(z)) return(data.frame())
  if(eiv_only) z <- z[z$variant=="EIV",]
  method_order <- if("method" %in% names(z)) z$method else z$comparator
  z <- z[order(match(z$cell,c("O0","O20","O50")),method_order),]
  out <- data.frame("Training exact U"=paste0(sub("O","",z$cell),"%"),check.names=FALSE)
  observed_fraction <- function(cell,method) {
    d<-evidence$oof[evidence$oof$cell==cell & evidence$oof$method==method & evidence$oof$prediction_mode=="mixed",]
    if(nrow(d))mean(d$test_U_observed) else 0
  }
  if("method" %in% names(z)) {
    route<-if(type=="information")rep("xc",nrow(z)) else if(type=="subgroups")rep("mixed",nrow(z)) else
      ifelse(z$prediction_mode=="class_only","xc","mixed")
    out$Method<-ocean_prediction_method(z$method,route,labels)
    if(type=="information")out$Method<-paste0(z$method,"(",labels[["xc"]],"/",labels[["mixed"]],")")
  }
  if("prediction_mode" %in% names(z))out$Test<-vapply(seq_len(nrow(z)),function(i) {
    fraction<-if(z$prediction_mode[i]=="class_only")0 else observed_fraction(z$cell[i],z$method[i])
    ocean_test_information(fraction,labels,include_label=TRUE)
  },character(1))
  if("group" %in% names(z))out[["Evaluation subset"]]<-vapply(seq_len(nrow(z)),function(i) {
    if(z$group[i]=="Exact U")return("100% (x,u)")
    if(z$group[i]=="Ordinal U")return("100% (x,c)")
    ocean_test_information(observed_fraction(z$cell[i],z$method[i]),labels,include_label=FALSE)
  },character(1))
  if("comparator" %in% names(z)) out$Comparator <- z$comparator
  out$N <- z$N; out$Folds <- z$n_folds
  score_labels <- c(mean_fold_RMSE="Mean fold RMSE",pooled_RMSE="Pooled RMSE",pooled_CRPS="CRPS / sY",
    RMSE="Pooled RMSE",CRPS="CRPS / sY",ordinal_RMSE=paste0("RMSE: ",labels[["xc"]]),exact_RMSE=paste0("RMSE: ",labels[["mixed"]]),
    delta_squared_error="Mean delta SE",delta_CRPS="Delta CRPS / sY",mean_abs_prediction_change="Mean |prediction change|")
  for(n in intersect(names(score_labels),names(z))) {
    v <- z[[n]]; if(grepl("CRPS",n)) v<-v/evidence$SD_Y
    out[[score_labels[[n]]]] <- sprintf(if(grepl("CRPS",n)) "%.4f" else "%.2f",v)
  }
  required_folds<-if(!is.null(evidence$expected_folds))length(evidence$expected_folds) else length(unique(evidence$oof$fold))
  # Highlight only losses on comparable, complete groups. Delta losses are not ranked.
  if(type %in% c("subgroups","pooled")) {
    groups<-if(type=="subgroups") paste(z$cell,z$group) else paste(z$cell,z$prediction_mode)
    complete<-vapply(groups,function(g) length(unique(z$method[groups==g]))==4L,logical(1))
    vals<-if(type=="subgroups") data.frame("Pooled RMSE"=z$RMSE,"CRPS / sY"=z$CRPS/evidence$SD_Y,check.names=FALSE) else
      data.frame("Pooled RMSE"=z$pooled_RMSE,"CRPS / sY"=z$pooled_CRPS/evidence$SD_Y,check.names=FALSE)
    attr(out,"best")<-ocean_best_mask(vals,groups,z$n_folds==required_folds & complete)
  } else if(type=="information") {
    vals<-cbind(z$ordinal_RMSE,z$exact_RMSE)
    mask<-matrix(FALSE,nrow(vals),2L,dimnames=list(NULL,c(paste0("RMSE: ",labels[["xc"]]),paste0("RMSE: ",labels[["mixed"]]))))
    for(i in seq_len(nrow(vals))) if(z$n_folds[i]==required_folds && all(is.finite(vals[i,])))
      mask[i,]<-abs(vals[i,]-min(vals[i,]))<=1e-12*max(1,abs(min(vals[i,])))
    attr(out,"best")<-mask
  }
  out
}

# Current OI/CC/MI predictions.
ocean_cc_prediction_tables <- function(output_dir,inputs,legacy,cc_output_dir="",
  profile="paper",dest=file.path(output_dir,"prediction_report","cc_revision"),
  expected_folds=NULL) {
  roots <- ocean_final_eiv_sources(output_dir,profile)
  if(is.null(expected_folds))expected_folds<-sort(unique(roots$fold))
  roots<-roots[roots$fold %in% expected_folds,,drop=FALSE]
  score_env <- new.env(parent=environment());ocean_initial_functions(score_env)
  rows <- list();add <- function(z,scope,training,test,source_profile=profile) {
    if(!nrow(z))return(invisible(NULL))
    if(!"n_train" %in% names(z))z$n_train<-vapply(seq_len(nrow(z)),function(i) {
      cell<-if(z$cell[i] %in% c("O0","baseline"))"O20" else z$cell[i]
      nrow(inputs$folds[[paste(z$fold[i],cell,sep="/")]]$train)
    },integer(1))
    z$scope<-scope;z$Training<-training;z$Test<-test;z$source_profile<-source_profile
    if(scope!="CC matched subset")z$Test<-vapply(seq_len(nrow(z)),function(i) {
      fraction<-0
      if(z$prediction_mode[i]=="mixed") {
        d<-inputs$folds[[paste(z$fold[i],z$cell[i],sep="/")]]
        fraction<-sum(d$observed_test_mixed)/z$n_test[i]
      }
      ocean_test_information(fraction,include_label=FALSE)
    },character(1))
    rows[[length(rows)+1L]]<<-z
  }
  base <- legacy$fold_metrics[legacy$fold_metrics$fold %in% expected_folds,,drop=FALSE]
  ac <- base[base$variant=="Original baseline",];ac$cell <- rep("O0",nrow(ac))
  add(ac,"AC","AC: 160 (x,c)","40 (x,c)")
  for(cell in c("O0","O20","O50")) {
    z<-base[base$variant=="EIV" & base$cell==cell & base$prediction_mode=="class_only",]
    if(nrow(z))add(z,if(cell=="O0")"AC" else "EIV all-profile",
      paste0(sub("O","",cell),"% calibration; 160 rows"),"40 (x,c)")
    if(cell!="O0") {
      z<-base[base$variant=="EIV" & base$cell==cell & base$prediction_mode=="mixed",]
      if(nrow(z))add(z,"EIV all-profile",paste0(sub("O","",cell),"% calibration; 160 rows"),
        if(cell=="O20")"20% (x,u) + 80% (x,c)" else "50% (x,u) + 50% (x,c)")
    }
  }
  cc<-details<-list()
  if(nzchar(cc_output_dir)) {
    tasks<-ocean_task_catalogue("paper");tasks<-tasks[tasks$module=="complete_case_baselines" & tasks$fold %in% expected_folds,]
    for(i in seq_len(nrow(tasks))) {
      t<-tasks[i,];dir<-file.path(cc_output_dir,"complete_case_baselines/fits",t$task_id)
      if(!file.exists(file.path(dir,"metrics.csv")))next
      if(!ocean_verify_receipt(cc_output_dir,t$task_id))stop("Unverified CC task: ",t$task_id)
      design<-readRDS(file.path(dir,"design.rds"));d<-inputs$folds[[paste(t$fold,t$cell,sep="/")]]
      stopifnot(identical(as.character(design$training_ids),as.character(d$train$cohort_id[d$observed_train])),
        identical(as.character(design$test_ids),as.character(d$test$cohort_id[d$observed_test_mixed])),
        identical(design$protocol,"native-CC-observed-u-subset-v1"))
      z<-read.csv(file.path(dir,"metrics.csv"));z$source_profile<-design$profile
      status<-read.csv(file.path(dir,"package_status.csv"))
      z$optimization_status<-status$optimization_status[match(t$method,status$method)]
      control<-design$controls[[t$method]]
      spec<-switch(t$method,`UC-GP`=paste0("starts=",control$n_starts),
        LVGP=paste0("starts=",control$n_starts,"; limits=",control$max_iter_ini,"/",control$max_iter_lat),
        EzGP=paste0("CV folds=",control$cv_folds,"; maxeval=",control$maxeval))
      details[[length(details)+1L]]<-data.frame(Fold=t$fold,Calibration=paste0(sub("O","",t$cell),"%"),
        Method=ocean_prediction_method(t$method,"cc",ocean_prediction_labels()),"Training N"=length(design$training_ids),
        "Test N"=length(design$test_ids),Controls=spec,Optimizer=z$optimization_status,check.names=FALSE)
      cc[[length(cc)+1L]]<-z
    }
  }
  cc<-ocean_bind(cc)
  if(!nrow(cc))cc<-data.frame(fold=integer(),cell=character(),method=character())
  for(cell in c("O20","O50")) {
    groups<-lapply(c("UC-GP(CC)","LVGP(CC)","EzGP(CC)"),function(m)cc$fold[cc$cell==cell & cc$method==m])
    eiv_folds<-roots$fold[roots$cell==cell & roots$available]
    matched<-Reduce(intersect,c(groups,list(eiv_folds)))
    if(!length(matched))next
    for(f in matched) {
      z<-cc[cc$cell==cell & cc$fold==f,]
      stopifnot(nrow(z)==3L,length(unique(z$source_profile))==1L)
      sp<-z$source_profile[1]
      add(z,"CC matched subset",paste0(sub("O","",cell),"% CC: ",if(cell=="O20")32 else 80," rows"),
        paste0(if(cell=="O20")8 else 20," (x,u,c)"),sp)
      s<-roots[roots$cell==cell & roots$fold==f,];dir<-file.path(output_dir,s$source)
      saved<-readRDS(file.path(dir,"prediction_distributions.rds"));k<-saved$known
      pred<-read.csv(file.path(dir,"mixed_predictions.csv"));d<-inputs$folds[[paste(f,cell,sep="/")]]
      stopifnot(identical(as.logical(k),as.logical(d$observed_test_mixed)),
        identical(as.character(pred$cohort_id),as.character(d$test$cohort_id)))
      draws<-saved$mixed[,k,drop=FALSE]
      for(a in c("conditional_means","conditional_vars"))attr(draws,a)<-attr(saved$mixed,a)[,k,drop=FALSE]
      z<-score_env$probability_score(draws,pred$predicted_Y[k],d$test$Y[k],f,"EIV-GP",cell)
      # Preserve the actual computation profile separately from display labels.
      add(z,"CC matched subset",paste0(sub("O","",cell),"% calibration; 160 rows"),
        paste0(sum(k)," (x,u)"),profile)
    }
  }
  raw<-ocean_bind(rows);summary<-list()
  if(nrow(raw))for(ids in split(seq_len(nrow(raw)),paste(raw$scope,raw$cell,raw$method,raw$Test,raw$source_profile))) {
    z<-raw[ids,];a<-z[1,c("scope","cell","method","Training","Test","source_profile")]
    a$n_folds<-length(unique(z$fold));a$N_per_fold<-unique(z$n_test)
    a$n_train<-unique(z$n_train)
    a$prediction_mode<-if("prediction_mode" %in% names(z))paste(unique(stats::na.omit(z$prediction_mode)),collapse="") else ""
    a$fold_ids<-paste(sort(unique(z$fold)),collapse=",")
    a$complete_folds<-setequal(z$fold,expected_folds) && !anyDuplicated(z$fold)
    a$N_total<-sum(z$n_test);a$pooled_RMSE<-sqrt(weighted.mean(z$RMSE^2,z$n_test))
    a$optimization_status<-paste(unique(stats::na.omit(z$optimization_status)),collapse="; ")
    for(m in c("RMSE","CRPS","Coverage95","Width95")) {
      a[[paste0(m,"_mean")]]<-mean(z[[m]]);a[[paste0(m,"_sd")]]<-if(nrow(z)>1)sd(z[[m]]) else NA_real_
    }
    summary[[length(summary)+1L]]<-a
  }
  summary<-ocean_bind(summary)
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  ocean_atomic_csv(raw,file.path(dest,"fold_metrics.csv"));ocean_atomic_csv(summary,file.path(dest,"summary.csv"))
  pooled_cc<-summary[summary$scope=="CC matched subset",c("cell","method","source_profile","n_folds","fold_ids","complete_folds","N_total","RMSE_mean","RMSE_sd","pooled_RMSE")]
  ocean_atomic_csv(pooled_cc,file.path(dest,"cc_pooled_rmse.csv"))
  list(raw=raw,summary=summary,cc=cc,dest=dest,pooled_cc=pooled_cc,computation_details=ocean_bind(details),expected_folds=expected_folds)
}
# Labels are presentation parameters; never infer fitting inputs from suffixes.
# OI / CC / MI are data-use protocols, with training/test scope defined together.
ocean_prediction_labels <- function(xc="OI",mixed="MI",exact="MI",cc="CC") {
  labels<-c(xc=xc,mixed=mixed,exact=exact,cc=cc)
  stopifnot(all(nzchar(labels)),!anyNA(labels));labels
}
# Test proportions are calculated from saved masks, not calibration labels.
ocean_test_information <- function(fraction,labels=ocean_prediction_labels(),include_label=TRUE) {
  stopifnot(length(fraction)==1L,is.finite(fraction),fraction>=0,fraction<=1)
  percent<-function(p)paste0(format(round(100*p,2L),trim=TRUE,scientific=FALSE),"%")
  if(fraction==0)return(paste0(if(include_label)paste0(labels[["xc"]],": ") else "","100% (x,c)"))
  if(fraction==1)return("100% (x,u)")
  paste0(if(include_label)paste0(labels[["mixed"]],": ") else "",percent(fraction)," (x,u) + ",percent(1-fraction)," (x,c)")
}
ocean_prediction_protocols <- function(labels=ocean_prediction_labels()) {
  descriptions<-c(xc="Use all training profiles; predict from (x,c) for every test profile.",
    cc="Fit only training profiles with observed u; evaluate only test profiles with observed u.",
    mixed="Use all training profiles; use (x,u) where test u is observed and (x,c) otherwise.")
  full_names<-c(OM="observable-mean protocol",MF="mixed m/f protocol",CC="complete-case protocol",
    OC="ordinal-input case",MC="mixed-input case",OI="ordinal-input protocol",MI="mixed-input protocol",
    `AT-OR`="all-training / ordinal-test protocol",`CC-EX`="complete-training / exact-test protocol",
    `AT-MX`="all-training / mixed-test protocol")
  keys<-c("xc","cc","mixed");titles<-unname(full_names[labels[keys]])
  titles[is.na(titles)]<-"data-use protocol"
  paste0("- **",labels[keys]," (",titles,"):** ",descriptions[keys],collapse="\n")
}
ocean_prediction_method <- function(method,route,labels) {
  base<-sub("(CC)","",method,fixed=TRUE)
  paste0(base,"(",unname(labels[route]),")")
}
# Dynamic eligibility: the prescribed folds and computation profiles must match.
# Numeric minima and ties are determined from unrounded losses, never row IDs.
ocean_prediction_ranks <- function(values,d,groups) {
  eligible<-d$complete_folds
  for(rows in split(seq_len(nrow(d)),groups))
    eligible[rows]<-eligible[rows] & length(unique(d$source_profile[rows]))==1L &
      length(unique(d$fold_ids[rows]))==1L
  ocean_best_mask(values,groups,eligible)
}
ocean_show_cc_predictions <- function(review,section=c("all","full","cc"),labels=ocean_prediction_labels()) {
  section<-match.arg(section);z<-review$summary
  display<-function(m,s,digits=2L)ifelse(is.na(s),sprintf(paste0("%.",digits,"f"),m),ocean_mean_sd(m,s,digits))
  full<-z[z$scope %in% c("AC","EIV all-profile"),,drop=FALSE]
  if(section %in% c("all","full") && nrow(full)) {
    full<-full[order(full$cell,full$Test,match(full$method,c("UC-GP","LVGP","EzGP","EIV-GP"))),,drop=FALSE]
    route<-ifelse(full$prediction_mode=="mixed","mixed","xc")
    out<-data.frame("Training exact U"=paste0(sub("O","",full$cell),"%"),
      "Test oxygen"=paste0(unname(labels[route]),": ",full$Test),
      Method=ocean_prediction_method(full$method,route,labels),"Training N"=full$n_train,
      Folds=full$n_folds,RMSE=display(full$RMSE_mean,full$RMSE_sd),
      "CRPS/sY"=display(full$CRPS_mean/review$SD_Y,full$CRPS_sd/review$SD_Y,3L),
      "Coverage (%)"=display(100*full$Coverage95_mean,100*full$Coverage95_sd,1L),
      "Width/sY"=display(full$Width95_mean/review$SD_Y,full$Width95_sd/review$SD_Y,3L),check.names=FALSE)
    vals<-data.frame(RMSE=full$RMSE_mean,"CRPS/sY"=full$CRPS_mean,check.names=FALSE)
    attr(out,"best")<-ocean_prediction_ranks(vals,full,paste(full$cell,full$Test))
    caption<-"Response prediction on the full test set"
    note<-paste("Panels indicate training calibration; test-input conditions are listed separately.",
      if(any(full$method!="EIV-GP"))"The (x,c) baselines ignore measured training u and are shown once in A as a shared reference." else "",
      "Entries are fold means (sample SD, where available). Bold marks the smallest unrounded RMSE/CRPS among methods with matching test profiles, completed prescribed folds and the same computation profile. Coverage and width are not ranked.")
    cat(ocean_prediction_panel_html(out,caption,note),"\n\n")
  }
  if(section=="full")return(invisible(NULL))
  d<-z[z$scope=="CC matched subset",,drop=FALSE]
  if(!nrow(d)) {cat("No matched complete-case results are available.\n");return(invisible(NULL))}
  d<-d[order(d$cell,match(sub("(CC)","",d$method,fixed=TRUE),c("UC-GP","LVGP","EzGP","EIV-GP"))),,drop=FALSE]
  method<-ocean_prediction_method(d$method,ifelse(d$method=="EIV-GP","exact","cc"),labels)
  flagged<-grepl("nonconverged",d$optimization_status,fixed=TRUE)

  out<-data.frame("Training exact U"=paste0(sub("O","",d$cell),"%"),
    "Test oxygen"=rep("100% (x,u); observed-u subset",nrow(d)),Method=method,"Training N"=d$n_train,
    Folds=d$n_folds,RMSE=display(d$RMSE_mean,d$RMSE_sd),
    "CRPS/sY"=display(d$CRPS_mean/review$SD_Y,d$CRPS_sd/review$SD_Y,3L),
    "Coverage (%)"=display(100*d$Coverage95_mean,100*d$Coverage95_sd,1L),
    "Width/sY"=display(d$Width95_mean/review$SD_Y,d$Width95_sd/review$SD_Y,3L),check.names=FALSE)
  attr(out,"best")<-ocean_prediction_ranks(data.frame(RMSE=d$RMSE_mean,"CRPS/sY"=d$CRPS_mean,check.names=FALSE),d,d$cell)
  note<-"Each block evaluates methods on identical observed-u test profiles. Native CC methods retain c as a categorical input. The mixed-input EIV prediction rule is evaluated only on its observed-u subset here, using f(x,u). Training N distinguishes CC fits from the full-data EIV fit. Entries are fold means (sample SD, where available); bold follows the same matching rule as the all-profile table. Pooled errors and actual fitting controls are in Supplementary."
  if(any(flagged))note<-paste(note,"Native fits with nonzero optimizer termination codes are documented in Supplementary.")
  cat(ocean_prediction_panel_html(out,"Complete-case comparators: matched observed-u test subset",note),"\n\n")
}
ocean_show_cc_pooled <- function(review,labels=ocean_prediction_labels()) {
  d<-review$pooled_cc;if(!nrow(d))return(invisible(NULL))
  d<-d[order(d$cell,match(sub("(CC)","",d$method,fixed=TRUE),c("UC-GP","LVGP","EzGP","EIV-GP"))),,drop=FALSE]
  out<-data.frame(Calibration=paste0(sub("O","",d$cell),"%"),
    Method=ocean_prediction_method(d$method,ifelse(d$method=="EIV-GP","exact","cc"),labels),
    Folds=d$n_folds,"Total test N"=d$N_total,"Pooled RMSE"=sprintf("%.2f",d$pooled_RMSE),check.names=FALSE)
  attr(out,"best")<-ocean_prediction_ranks(data.frame("Pooled RMSE"=d$pooled_RMSE,check.names=FALSE),d,d$cell)
  cat(ocean_table_html(out,"Pooled RMSE on matched observed-u test profiles",
    "Pooled RMSE is sqrt(sum(n_fold * RMSE_fold^2) / sum(n_fold)). Total test N and fold counts reflect available matched results. With one fold, pooled and fold RMSE coincide. Bold uses the same matching rule as the main tables."))
  if(nrow(review$computation_details))cat(ocean_table_html(review$computation_details,
    "Native complete-case fitting controls and optimizer status",
    "Controls are read from each saved fit design; termination codes are retained from the native packages. EIV numerical settings and posterior diagnostics are described in Section 3."))
}
