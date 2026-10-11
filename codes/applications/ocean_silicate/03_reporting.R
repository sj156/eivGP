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

# Read scores only; no fitting or test-score-based choice of placeholder.
ocean_bind <- function(parts) {
  if(!length(parts)) return(data.frame())
  cols <- unique(unlist(lapply(parts,names)))
  do.call(rbind,lapply(parts,function(x) {for(n in setdiff(cols,names(x))) x[[n]]<-NA; x[,cols,drop=FALSE]}))
}
ocean_prediction_report <- function(output_dir,inputs,profile="paper",augmented_output_dir="",report_dir=file.path(output_dir,"prediction_report")) {
  rows <- list(); inventory <- list(); original <- if(profile=="paper") "full5_class4" else "smoke_class4"
  final_sources <- ocean_final_eiv_sources(output_dir,profile)
  add <- function(path,f,cell,variant,b=NA_real_) {
    if(!file.exists(path)) return(invisible(NULL))
    z <- read.csv(path,stringsAsFactors=FALSE)
    stopifnot(nrow(z)>0L,all(z$fold==f),all(z$n_test==40L))
    if(variant=="EIV") {
      z$prediction_mode <- ifelse(z$method=="EIV-GP xc","class_only","mixed")
      if(cell=="O0") z <- z[z$prediction_mode=="class_only",,drop=FALSE]
      z$method <- "EIV-GP"
    } else if(variant=="Original baseline") z$prediction_mode <- "class_only"
    z$variant <- variant; z$placeholder <- b
    z$source <- gsub("\\\\","/",path)
    rows[[length(rows)+1L]] <<- z
  }
  for(f in if(profile=="paper") 1:5 else 1L) {
    initial <- file.path(output_dir,original,sprintf("fold%02d",f))
    add(file.path(initial,"baseline_metrics.csv"),f,"baseline","Original baseline")
    for(cell in if(profile=="paper") c("O0","O20","O50") else c("O0","O50")) {
      final<-final_sources[final_sources$fold==f & final_sources$cell==cell,]
      path<-file.path(output_dir,final$source)
      available<-isTRUE(final$available)
      inventory[[length(inventory)+1L]]<-data.frame(fold=f,cell,
        stage=if(final$continued)"plus20k" else "initial",available,
        initial_screen_pass=final$initial_screen_pass,final_screen_pass=final$final_screen_pass)
      # A cell that needs continuation never falls back to its initial fit.
      if(available)add(file.path(path,"predictive_metrics.csv"),f,cell,"EIV")
    }
  }
  aug <- if(nzchar(augmented_output_dir)) augmented_output_dir else file.path(output_dir,"augmented_baselines")
  if(profile=="paper") {
    tasks <- ocean_task_catalogue("paper"); tasks <- tasks[tasks$module=="augmented_baselines",]
    for(i in seq_len(nrow(tasks))) {
      task <- tasks[i,]; add(file.path(aug,"fits",task$task_id,"metrics.csv"),
        task$fold,task$cell,"Expanded baseline",task$placeholder)
    }
  }
  raw <- ocean_bind(rows); summaries <- list()
  metrics <- c("RMSE","MAE","CRPS","NLPD","Coverage95","Width95","IntervalScore95")
  if(nrow(raw)) {
    key <- paste(raw$variant,raw$method,raw$cell,raw$prediction_mode,raw$placeholder,sep="/")
    if(anyDuplicated(paste(key,raw$fold))) stop("Duplicate fold/task scores in report.")
    for(group in split(seq_len(nrow(raw)),key)) {
      z <- raw[group,,drop=FALSE]
      ans <- z[1,c("variant","method","cell","prediction_mode","placeholder"),drop=FALSE]
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
  # Presentation is restricted to the prespecified primary encoding b=0.
  # Keep all saved scores in raw audit exports; do not select b from test scores.
  plots <- list()
  if(nrow(summary_raw) && requireNamespace("ggplot2",quietly=TRUE)) {
    d <- summary_raw[(summary_raw$variant=="Expanded baseline" &
      !is.na(summary_raw$placeholder) & summary_raw$placeholder==0 |
      summary_raw$variant=="EIV") & summary_raw$cell %in% c("O20","O50") &
      summary_raw$complete_5fold,,drop=FALSE]
    if(nrow(d)) {
      d$method <- factor(d$method,levels=c("UC-GP","LVGP","EzGP","EIV-GP"))
      d$calibration <- factor(d$cell,levels=c("O20","O50"),labels=c("20% exact training oxygen","50% exact training oxygen"))
      d$test_input <- factor(d$prediction_mode,levels=c("class_only","mixed"),
        labels=c("Class-only test inputs","Matched mixed test inputs"))
      plots$primary_comparison <- ggplot2::ggplot(d,ggplot2::aes(RMSE_mean,method,color=method)) +
        ggplot2::geom_linerange(ggplot2::aes(xmin=RMSE_mean-RMSE_sd,xmax=RMSE_mean+RMSE_sd),
          orientation="y",linewidth=.65) + ggplot2::geom_point(size=2.8) +
        ggplot2::facet_grid(test_input~calibration) +
        ggplot2::scale_color_manual(values=c("UC-GP"="#7b3294",LVGP="#008837",EzGP="#a6761d","EIV-GP"="#2365a3"),guide="none") +
        ggplot2::labs(x="Silicate prediction RMSE (micromol/kg)",y=NULL,
          title="Prediction with matched available oxygen measurements",
          caption="Expanded baselines use the prespecified b=0 encoding. Points: five-fold means; bars: sample SD.\nMixed test inputs reveal exact oxygen for 20%/50% of test profiles, matching training calibration.") +
        ggplot2::theme_minimal(base_size=11) +
        ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),panel.grid.major.y=ggplot2::element_blank(),
          plot.caption=ggplot2::element_text(hjust=0))
    }
  }
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
    title <- paste0(LETTERS[j],". ",training_groups[j]," exact training oxygen")
    if(length(tests)==1L) title <- paste0(title,"; test: ",tests)
    out <- c(out,paste0("<tr class='comparison-panel'><th scope='rowgroup' colspan='",length(columns),
      "' style='text-align:left;background:#f2f4f6;border-top:2px solid #bbb;padding:9px 10px'>",
      ocean_html_escape(title),"</th></tr>"))
    for(test in tests) {
      rows <- training_rows[d[["Test oxygen"]][training_rows]==test]
      if(length(tests)>1L) out <- c(out,paste0(
        "<tr class='comparison-condition'><th scope='rowgroup' colspan='",length(columns),
        "' style='text-align:left;font-weight:500;padding:7px 10px;border-bottom:1px solid #ddd'>",
        "Test oxygen: ",ocean_html_escape(test),"</th></tr>"))
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
ocean_prediction_panel_tex <- function(d,columns,headers,label,caption,note) {
  styled <- ocean_table_markup(d,"latex")
  tex <- c("\\begin{table}[t]","\\centering\\small",paste0("\\caption{",caption,"}"),
    paste0("\\label{tab:",label,"}"),"\\setlength{\\tabcolsep}{4pt}","\\renewcommand{\\arraystretch}{1.1}",
    paste0("\\begin{tabular}{l",paste(rep("r",length(columns)),collapse=""),"}"),
    "\\toprule",paste0("Method & ",paste(headers,collapse=" & "),"\\\\"),"\\midrule")
  training_groups <- unique(d[["Training exact U"]])
  for(j in seq_along(training_groups)) {
    training_rows <- which(d[["Training exact U"]]==training_groups[j])
    tests <- unique(d[["Test oxygen"]][training_rows])
    if(j>1) tex <- c(tex,"\\addlinespace[4pt]")
    title <- paste0(LETTERS[j],". ",training_groups[j]," exact training oxygen")
    if(length(tests)==1L) title <- paste0(title,"; test: ",tests)
    tex <- c(tex,paste0("\\multicolumn{",length(columns)+1L,"}{l}{\\textit{",ocean_tex_note(title),"}}\\\\[2pt]"))
    for(test in tests) {
      rows <- training_rows[d[["Test oxygen"]][training_rows]==test]
      if(length(tests)>1L) tex <- c(tex,paste0("\\multicolumn{",length(columns)+1L,
        "}{l}{Test oxygen: ",ocean_tex_note(test),"}\\\\[2pt]"))
      for(i in rows) tex <- c(tex,paste0(sub("*","$^{*}$",d$Method[i],fixed=TRUE)," & ",
        paste(unlist(styled[i,columns,drop=FALSE],use.names=FALSE),collapse=" & "),"\\\\"))
    }
  }
  c(tex,"\\bottomrule","\\end{tabular}","\\par\\smallskip","\\begin{minipage}{\\linewidth}\\footnotesize",
    ocean_tex_note(note),"\\end{minipage}","\\end{table}")
}
ocean_paper_prediction_tables <- function(report) {
  z <- report$summary
  if(!nrow(z)) return(list(main=data.frame(),full=data.frame(),numeric=data.frame()))
  z <- z[z$variant!="Expanded baseline" | (!is.na(z$placeholder) & z$placeholder==0),,drop=FALSE]
  z$training_cell <- ifelse(z$variant=="Original baseline","O0",z$cell)
  z <- z[order(match(z$training_cell,c("O0","O20","O50")),match(z$prediction_mode,c("class_only","mixed")),
    match(z$method,c("UC-GP","LVGP","EzGP","EIV-GP"))),,drop=FALSE]
  base <- data.frame("Training exact U"=paste0(sub("O","",z$training_cell),"%"),
    "Test oxygen"=ifelse(z$prediction_mode=="class_only","Class only",paste0(sub("O","",z$training_cell),"% exact")),
    Method=ifelse(z$variant=="Expanded baseline",paste0(z$method,"*"),z$method),check.names=FALSE)
  if(any(!z$complete_5fold)) base$Method <- ifelse(z$complete_5fold,base$Method,
    paste0(base$Method," [n=",z$n_folds," folds]"))
  format_metric <- function(metric,multiplier=1,digits=2,offset=0)
    ocean_mean_sd(multiplier*z[[paste0(metric,"_mean")]]-offset,multiplier*z[[paste0(metric,"_sd")]],digits)
  sy <- report$SD_Y
  main <- base; main$RMSE <- format_metric("RMSE")
  main[["CRPS / sY"]] <- format_metric("CRPS",1/sy,3)
  main[["Coverage (%)"]] <- format_metric("Coverage95",100,1)
  main[["Width / sY"]] <- format_metric("Width95",1/sy,3)
  full <- base; full$RMSE <- format_metric("RMSE")
  full[["RMSE / sY"]] <- format_metric("RMSE",1/sy,3)
  full[["MAE / sY"]] <- format_metric("MAE",1/sy,3)
  full[["CRPS / sY"]] <- format_metric("CRPS",1/sy,3)
  full[["NLPD standardized"]] <- format_metric("NLPD",1,3,log(sy))
  full[["Coverage (%)"]] <- format_metric("Coverage95",100,1)
  full[["Width / sY"]] <- format_metric("Width95",1/sy,3)
  full[["Interval score / sY"]] <- format_metric("IntervalScore95",1/sy,3)
  groups <- paste(base[["Training exact U"]],base[["Test oxygen"]],sep="/")
  main_values <- data.frame(RMSE=z$RMSE_mean,"CRPS / sY"=z$CRPS_mean/sy,check.names=FALSE)
  attr(main,"best") <- ocean_best_mask(main_values,groups,z$complete_5fold)
  full_values <- data.frame(RMSE=z$RMSE_mean,"RMSE / sY"=z$RMSE_mean/sy,
    "MAE / sY"=z$MAE_mean/sy,"CRPS / sY"=z$CRPS_mean/sy,
    "NLPD standardized"=z$NLPD_mean-log(sy),
    "Interval score / sY"=z$IntervalScore95_mean/sy,check.names=FALSE)
  attr(full,"best") <- ocean_best_mask(full_values,groups,z$complete_5fold)
  list(main=main,full=full,numeric=z)
}
ocean_write_prediction_tables <- function(tables,dest,sy) {
  if(!nrow(tables$main)) return(invisible(NULL))
  note <- paste0("Entries are fold means (sample SDs). sY = ",sprintf("%.4f",sy),
    " micromol/kg, the frozen cohort response SD. RMSE is in micromol/kg; coverage is percent. ",
    "* marks baselines with the same recorded oxygen information as EIV-GP: retain (x,c), add categorical observation status R and continuous oxygen record W, and use b=0 for hidden W. Native kernels and fitting procedures are used. ",
    "At 0%, they use their original (x,c) inputs. Test inputs and held-out profiles are matched within each block. ",
    "Bold marks the smallest unrounded mean error or predictive score within a complete five-fold block; SDs are not bold. Coverage and width are descriptive.")
  ocean_export_paper_table(tables$main,"paper_y_prediction_main","Silicate prediction with matched oxygen information",dest,note)
  ocean_export_paper_table(tables$full,"paper_y_prediction_full","Additional silicate prediction summaries: b=0",dest,note)
  for(name in c("main","full")) {
    content <- ocean_prediction_panel_html(tables[[name]],"Silicate prediction with matched oxygen information",note)
    ocean_atomic_lines(c("<!doctype html><meta charset='utf-8'>",
      "<style>body{font:15px/1.5 sans-serif;margin:24px}.note{color:#555}</style>",content),
      file.path(dest,paste0("paper_y_prediction_",name,".html")))
  }
  ocean_atomic_lines(ocean_prediction_panel_tex(tables$main,
    c("RMSE","CRPS / sY","Coverage (%)","Width / sY"),
    c("RMSE","CRPS/$s_Y$","Cov. (\\%)","Width/$s_Y$"),"paper_y_prediction_main",
    "Silicate prediction under matched oxygen information",note),file.path(dest,"paper_y_prediction_main.tex"))
  ocean_atomic_lines(c(ocean_prediction_panel_tex(tables$full,
    c("RMSE","RMSE / sY","MAE / sY"),c("RMSE","RMSE/$s_Y$","MAE/$s_Y$"),
    "paper_y_prediction_full_errors","Point prediction summaries: primary encoding b=0",note),"",
    ocean_prediction_panel_tex(tables$full,
    c("CRPS / sY","NLPD standardized","Coverage (%)","Width / sY","Interval score / sY"),
    c("CRPS/$s_Y$","NLPD (std.)","Cov. (\\%)","Width/$s_Y$","IS/$s_Y$"),
    "paper_y_prediction_full_distribution","Predictive distribution summaries: primary encoding b=0",note)),
    file.path(dest,"paper_y_prediction_full.tex"))
  ocean_atomic_csv(tables$numeric,file.path(dest,"primary_summary_numeric.csv"))
}
ocean_export_prediction_report <- function(report,dest) {
  ocean_atomic_csv(report$fold_metrics,file.path(dest,"fold_metrics_raw.csv"))
  ocean_atomic_csv(report$summary,file.path(dest,"summary_raw.csv"))
  ocean_atomic_csv(report$fold_metrics_SD_Y,file.path(dest,"fold_metrics_SD_Y.csv"))
  ocean_atomic_csv(report$availability,file.path(dest,"eiv_availability.csv"))
  for(name in names(report$plots)) for(ext in c("png","pdf")) ggplot2::ggsave(
    file.path(dest,paste0(name,".",ext)),report$plots[[name]],width=10,height=6,dpi=180,bg="white")
  ocean_write_prediction_tables(ocean_paper_prediction_tables(report),dest,report$SD_Y)
  invisible(report)
}
ocean_show_predictions <- function(report, show_plot = FALSE) {
  if(is.null(report) || !nrow(report$summary)) {cat("No completed prediction scores yet.\n");return(invisible(NULL))}
  tables <- ocean_paper_prediction_tables(report)
  cat("\n\n### Silicate prediction\n\n",
    "Entries are fold means (sample SDs). RMSE is in micromol/kg; sY is the frozen cohort response SD. ",
    "* marks UC-GP, LVGP and EzGP with the same recorded training/test oxygen measurements as EIV-GP. ",
    "Incomplete groups state their fold count in the method label.\n\n",sep="")
  cat(ocean_prediction_panel_html(tables$main,"Silicate prediction under matched oxygen information",
    "Bold marks the smallest RMSE and CRPS mean within each complete five-fold information block, using unrounded values. Coverage and width are descriptive."),"\n\n")
  if(show_plot) for(p in report$plots) print(p)
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
      stage <- if(isTRUE(row$continued)) "+20k" else "Initial"
      paste(stage, if(isTRUE(row$final_screen_pass)) "pass" else "review", sep = "; ")
    }, "")
  }
  cat(as.character(knitr::kable(overview, format = "pipe", row.names = FALSE)), sep = "\n")
  cat("\n\n")
  if(settings$profile == "paper" && settings$action != "plan") {
    cat("Initial / +20k identifies the saved fit used in the report. Pass refers to the prespecified finite-Rhat/bulk-ESS screen; full chain diagnostics remain in the result files.\n\n")
    available <- status$available & !is.na(status$final_screen_pass)
    if(any(available & !status$final_screen_pass, na.rm = TRUE))
      cat("Fits marked review remain outside this limited screen.\n\n")
  }
  invisible(overview)
}
