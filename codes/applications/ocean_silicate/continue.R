#!/usr/bin/env Rscript

## Continue existing eivGP 0.3.1 chains, selected only by training-chain
## diagnostics. No test-Y-based selection, refitting, or package modification.
options(warn = 1)
root <- OCEAN_CASE_ROOT
suppressPackageStartupMessages(library(eivGP))
stopifnot(as.character(packageVersion("eivGP")) == "0.3.1")

out <- file.path(root, "continuation_plus20k")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
progress_path <- file.path(out, "PROGRESS.txt")
status_path <- file.path(out, "task_status.csv")
started <- Sys.time()

atomic_lines <- function(x, path) {
  tmp <- tempfile(".atomic-", tmpdir = dirname(path))
  writeLines(x, tmp)
  if (!file.rename(tmp, path)) stop("Cannot publish ", path)
}
atomic_csv <- function(x, path) {
  tmp <- tempfile(".atomic-", tmpdir = dirname(path), fileext = ".csv")
  write.csv(x, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) stop("Cannot publish ", path)
}
atomic_rds <- function(x, path) {
  tmp <- tempfile(".atomic-", tmpdir = dirname(path), fileext = ".rds")
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) stop("Cannot publish ", path)
}
original_path <- function(fold, cell) {
  file.path(root, "full5_class4", sprintf("fold%02d", fold), cell, "fit.rds")
}
diagnostic_row <- function(fit) {
  tables <- fit$diagnostics[c("rhat_hyper", "rhat_tau", "rhat_u")]
  rhats <- unlist(lapply(tables, function(z) z$rhat))
  bulk <- unlist(lapply(tables, function(z) z$ess_bulk))
  max_rhat <- max(rhats[is.finite(rhats)])
  min_ess <- min(bulk[is.finite(bulk)])
  data.frame(max_rhat = max_rhat, min_ess_bulk = min_ess,
    pass_screen = max_rhat <= 1.01 && min_ess >= 400)
}


if (!file.exists(file.path(root, "full5_class4", "complete.flag"))) {
  atomic_lines(c(
    "state: WAITING", paste0("updated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    "dependency: waiting for original full-five-fold experiment to finish",
    "plan: use installed continue_eivgp() to add 20000 transitions per chain to fits failing Rhat<=1.01 and bulk ESS>=400",
    "selection: chain diagnostics only; no held-out Y or RMSE used"
  ), progress_path)
  stop("Complete the initial ocean fits before continuation.")
}

tasks <- if (file.exists(status_path)) {
  read.csv(status_path, stringsAsFactors = FALSE)
} else {
  rows <- list()
  for (fold in 1:5) for (cell in c("O0", "O20", "O50")) {
    source_diag <- read.csv(file.path(dirname(original_path(fold, cell)), "convergence_summary.csv"))
    stopifnot(nrow(source_diag) == 1L)
    d <- data.frame(max_rhat = source_diag$max_finite_rhat,
      min_ess_bulk = source_diag$min_finite_ess_bulk,
      pass_screen = source_diag$limited_screen_pass)
    stopifnot(is.finite(d$max_rhat), is.finite(d$min_ess_bulk))
    rows[[length(rows) + 1L]] <- data.frame(
      fold = fold, cell = cell, original_max_rhat = d$max_rhat,
      original_min_ess_bulk = d$min_ess_bulk,
      selected = !d$pass_screen,
      state = if (d$pass_screen) "NOT_SELECTED" else "PENDING",
      updated = "", minutes = NA_real_, continued_max_rhat = NA_real_,
      continued_min_ess_bulk = NA_real_, continued_pass_screen = NA,
      detail = "", stringsAsFactors = FALSE)

  }
  x <- do.call(rbind, rows)
  atomic_csv(x, status_path)
  x
}
stopifnot(nrow(tasks) == 15L, sum(tasks$selected) == 3L)
expected_selected <- c("1 O0", "3 O20", "4 O0")
stopifnot(setequal(paste(tasks$fold[tasks$selected], tasks$cell[tasks$selected]), expected_selected))

set_task <- function(fold, cell, state, detail = "", minutes = NA_real_, diag = NULL) {
  i <- which(tasks$fold == fold & tasks$cell == cell)
  stopifnot(length(i) == 1L)
  tasks$state[i] <<- state
  tasks$updated[i] <<- format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
  tasks$detail[i] <<- detail
  if (is.finite(minutes)) tasks$minutes[i] <<- minutes
  if (!is.null(diag)) {
    tasks$continued_max_rhat[i] <<- diag$max_rhat
    tasks$continued_min_ess_bulk[i] <<- diag$min_ess_bulk
    tasks$continued_pass_screen[i] <<- diag$pass_screen
  }
  atomic_csv(tasks, status_path)
}
write_progress <- function(state, current, detail = "") {
  selected <- tasks$selected
  initial <- c(O0 = 48, O20 = 38, O50 = 30)
  estimates <- vapply(names(initial), function(cell) {
    completed <- tasks$minutes[tasks$cell == cell & tasks$state == "DONE" & is.finite(tasks$minutes)]
    if (length(completed)) median(completed) else initial[[cell]]
  }, numeric(1L))
  remaining <- sum(estimates[tasks$cell[selected & tasks$state != "DONE"]])
  atomic_lines(c(
    paste0("state: ", state),
    paste0("updated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste0("pid: ", Sys.getpid()),
    "dependency: original U<180 pressure300-1500 joint-stratified full-five-fold experiment complete",
    "selection: max finite Rhat > 1.01 or min finite bulk ESS < 400; no test Y used",
    "method: installed eivGP 0.3.1 continue_eivgp(), +20000 transitions per chain, 4 cores",
    "controls: frozen data/folds/training calibration/test masks; hard-threshold engine; u_block_size=8 retained",
    "chain length: original20000 + new20000 =40000 per chain; original burn5000 retained; no second burn-in",
    "launch: one-shot launchd RunAtLoad=true, KeepAlive=false; no recurring restarts",
    "original fit files and baseline outputs are never overwritten",
    paste0("selected_fits: ", sum(selected), "/15"),
    paste0("completed_selected: ", sum(selected & tasks$state == "DONE"), "/", sum(selected)),
    paste0("continued_passing_screen: ", sum(tasks$continued_pass_screen %in% TRUE)),
    paste0("current: ", current),
    paste0("detail: ", detail),
    paste0("elapsed_this_launch_hours: ", round(as.numeric(difftime(Sys.time(), started, units = "hours")), 2)),
    paste0("estimated_remaining_hours: ", round(remaining / 60, 2)),
    paste0("estimated_finish: ", format(Sys.time() + remaining * 60, "%Y-%m-%d %H:%M:%S %z")),
    "eta_note: task-boundary estimate; longer chains do not guarantee convergence",
    "", capture.output(print(tasks, row.names = FALSE))
  ), progress_path)
}

score <- function(pred, y, fold, cell, method) {
  pred <- as.numeric(pred)
  stopifnot(length(pred) == length(y), all(is.finite(pred)))
  data.frame(fold = fold, cell = cell, method = method, n_test = length(y),
    RMSE = sqrt(mean((pred - y)^2)), MAE = mean(abs(pred - y)))
}


probability_score<-function(draws,point,y,fold,cell,method) {
  q<-apply(draws,2L,quantile,c(.025,.975))
  nlpd<-eivGP:::normal_mixture_nlpd(y,attr(draws,"conditional_means",exact=TRUE),attr(draws,"conditional_vars",exact=TRUE))
  cbind(score(point,y,fold,cell,method),CRPS=mean(eivGP:::crps_sample_matrix(draws,y)),NLPD=nlpd$value,
    Coverage95=mean(y>=q[1,]&y<=q[2,]),Width95=mean(q[2,]-q[1,]),
    IntervalScore95=mean(eivGP:::interval_score(q[1,],q[2,],y)),NLPD_reason=nlpd$reason)
}
publish_latest_summary<-function() {
  metrics<-list();convergence<-list()
  for(f in 1:5) {
    b<-read.csv(file.path(root,"full5_class4",sprintf("fold%02d",f),"baseline_metrics.csv"))
    b$source_stage<-"original";metrics[[length(metrics)+1L]]<-b
    for(cell in c("O0","O20","O50")) {
      new<-file.path(out,sprintf("fold%02d",f),cell)
      extended<-file.exists(file.path(new,"predictive_metrics.csv"))
      p<-if(extended)new else file.path(root,"full5_class4",sprintf("fold%02d",f),cell)
      m<-read.csv(file.path(p,"predictive_metrics.csv"));m$source_stage<-if(extended)"plus20k" else "original"
      metrics[[length(metrics)+1L]]<-m
      d<-read.csv(file.path(p,"convergence_summary.csv"));d$fold<-f;d$cell<-cell;d$source_stage<-if(extended)"plus20k" else "original"
      convergence[[length(convergence)+1L]]<-d
    }
  }
  all<-do.call(rbind,metrics);stopifnot(nrow(all)==45L)
  atomic_csv(all,file.path(out,"latest_fold_predictive_metrics.csv"))
  atomic_csv(do.call(rbind,convergence),file.path(out,"latest_convergence_summary.csv"))
  label<-ifelse(all$cell=="baseline",all$method,paste(all$cell,all$method))
  ysd<-sd(read.csv(file.path(root,"data","cohort_manifest.csv"))$Y)
  sumrows<-lapply(split(all,label),function(z) {
    ans<-data.frame(method=if(z$cell[1]=="baseline")z$method[1] else paste(z$cell[1],z$method[1]),
      n_test=sum(z$n_test),common_Y_SD=ysd,pooled_RMSE=sqrt(weighted.mean(z$RMSE^2,z$n_test)),
      fold_RMSE_mean=mean(z$RMSE),fold_RMSE_sd=sd(z$RMSE))
    ans$RMSE_over_common_Y_SD<-ans$pooled_RMSE/ysd
    for(k in c("MAE","CRPS","NLPD","Coverage95","Width95","IntervalScore95")) {
      ans[[paste0("pooled_",k)]]<-weighted.mean(z[[k]],z$n_test)
      ans[[paste0("fold_",k,"_mean")]]<-mean(z[[k]])
      ans[[paste0("fold_",k,"_sd")]]<-sd(z[[k]])
    }
    ans
  })
  atomic_csv(do.call(rbind,sumrows),file.path(out,"latest_summary_metrics.csv"))
}

run_one <- function(fold, cell) {
  task_out <- file.path(out, sprintf("fold%02d", fold), cell)
  dir.create(task_out, recursive = TRUE, showWarnings = FALSE)
  metrics_path <- file.path(task_out, "metrics.csv")
  fit_path <- file.path(task_out, "fit_plus20k.rds")
  if (file.exists(metrics_path)) {
    fit <- readRDS(fit_path)
    set_task(fold, cell, "DONE", "reused complete continuation",
             diag = diagnostic_row(fit))
    return(invisible(NULL))
  }
  set_task(fold, cell, "RUNNING", "continuing four chains")
  write_progress("RUNNING", sprintf("fold%02d %s", fold, cell), "continuing original chains")
  t0 <- Sys.time()
  fit <- if (file.exists(fit_path)) readRDS(fit_path) else {
    original <- readRDS(original_path(fold, cell))
    continued <- continue_eivgp(original, n_iter = 20000L,
      parallel = TRUE, n_cores = 4L, verbose = TRUE)
    atomic_rds(continued, fit_path)
    rm(original); gc(FALSE)
    continued
  }
  diag <- diagnostic_row(fit)
  diagnostic_parts<-lapply(fit$diagnostics[c("rhat_hyper","rhat_tau","rhat_u")],function(z)z[,c("parameter","rhat","ess_bulk","ess_tail","target_pass"),drop=FALSE])
  atomic_csv(do.call(rbind,diagnostic_parts),file.path(task_out,"convergence_parameters.csv"))
  atomic_csv(data.frame(max_finite_rhat=diag$max_rhat,min_finite_ess_bulk=diag$min_ess_bulk,limited_screen_pass=diag$pass_screen),file.path(task_out,"convergence_summary.csv"))
  atomic_csv(data.frame(fold = fold, cell = cell, stage = "plus20k",
    max_rhat = diag$max_rhat, min_ess_bulk = diag$min_ess_bulk,
    pass_screen = diag$pass_screen), file.path(task_out, "diagnostics.csv"))
  capture.output(summary(fit), file = file.path(task_out, "fit_summary.txt"))

  train <- read.csv(file.path(root, "data", sprintf("fold%02d_train.csv", fold)))
  test <- read.csv(file.path(root, "data", sprintf("fold%02d_test.csv", fold)))
  train$C <- as.integer(findInterval(train$U, c(50, 100, 150)) + 1L)
  test$C <- as.integer(findInterval(test$U, c(50, 100, 150)) + 1L)
  stopifnot(nrow(train) == 160L, nrow(test) == 40L,
    !length(intersect(train$profile, test$profile)),
    all(table(test$C) == 10L))
  x_test <- as.matrix(test[, c("G2salinity", "log_pressure")])
  c_test <- matrix(as.integer(test$C), ncol = 1L,
                   dimnames = list(NULL, "oxygen_class"))
  n_draw <- nrow(fit$mcmc$samples_tau)
  set.seed(2026092800L + fold * 100L + match(cell, c("O0", "O20", "O50")))
  draw_ids <- sort(sample(seq_len(n_draw), min(2000L, n_draw)))
  response_draws <- predict_eivgp(fit, new_X = x_test, new_C = c_test,
    target = "response", draw_ids = draw_ids, n_per_draw = 1L,
    joint = FALSE, seed = 2026092900L + fold * 100L + match(cell, c("O0", "O20", "O50")))
  point <- colMeans(as.matrix(response_draws))
  metrics <- score(point, test$Y, fold, cell, "EIV-GP")
  atomic_csv(data.frame(cohort_id = test$cohort_id, observed_Y = test$Y,
    predicted_Y = point), file.path(task_out, "predictions.csv"))

  mask <- if (cell == "O0") rep(FALSE, nrow(train)) else
    train[[paste0("calib_", tolower(cell))]] == 1L

  # Reuse the exact original test observation masks, not a new random selection.
  saved_masks<-read.csv(file.path(root,"full5_class4",sprintf("fold%02d",fold),"test_U_masks.csv"))
  stopifnot(identical(as.character(saved_masks$cohort_id),as.character(test$cohort_id)))
  known<-if(cell=="O0")rep(FALSE,nrow(test)) else as.logical(saved_masks[[paste0("observed_U_",cell)]])
  stopifnot(sum(known)==switch(cell,O0=0L,O20=8L,O50=20L))
  mixed<-response_draws
  if(any(known)) {
    xu<-predict_eivgp(fit,new_X=x_test[known,,drop=FALSE],new_U=matrix(test$U[known],ncol=1L),new_U_scale="raw",
      target="response",draw_ids=draw_ids,n_per_draw=1L,joint=FALSE,seed=2026093000L+fold*100L+match(cell,c("O0","O20","O50")))
    mixed[,known]<-xu
    for(a in c("conditional_means","conditional_vars")) {
      comp<-attr(response_draws,a,exact=TRUE);comp[,known]<-attr(xu,a,exact=TRUE);attr(mixed,a)<-comp
    }
  }
  stopifnot(identical(mixed[,!known,drop=FALSE],response_draws[,!known,drop=FALSE]))
  atomic_rds(list(xc=response_draws,mixed=mixed,known=known,draw_ids=draw_ids),file.path(task_out,"prediction_distributions.rds"))
  atomic_csv(data.frame(cohort_id=test$cohort_id,observed_Y=test$Y,predicted_Y=colMeans(mixed),test_U_observed=known),file.path(task_out,"mixed_predictions.csv"))
  atomic_csv(rbind(probability_score(response_draws,point,test$Y,fold,cell,"EIV-GP xc"),
    probability_score(mixed,colMeans(mixed),test$Y,fold,cell,"EIV-GP mixed")),file.path(task_out,"predictive_metrics.csv"))
  hidden <- which(!mask)
  imp <- impute_eivgp(fit, rows = hidden, draw_ids = draw_ids,
                      scale = if (cell == "O0") "model" else "raw")
  imp_point <- colMeans(matrix(imp[, , 1L], nrow = length(draw_ids)))
  if (cell == "O0") {
    atomic_csv(data.frame(cohort_id = train$cohort_id[hidden], C = train$C[hidden],
      true_U_raw = train$U[hidden], EIV_U_model_mean = imp_point),
      file.path(task_out, "imputations_model_scale.csv"))
  } else {
    calib_means <- tapply(train$U[mask], train$C[mask], mean)
    imp_base <- as.numeric(calib_means[as.character(train$C[hidden])])
    atomic_csv(data.frame(cohort_id = train$cohort_id[hidden], C = train$C[hidden],
      true_U = train$U[hidden], EIV_U_mean = imp_point,
      class_calibration_mean = imp_base),
      file.path(task_out, "imputations.csv"))
    metrics <- rbind(metrics,
      score(imp_point, train$U[hidden], fold, cell, "EIV U imputation"),
      score(imp_base, train$U[hidden], fold, cell, "calibration class mean"))
  }
  iq<-apply(matrix(imp[,,1L],nrow=length(draw_ids)),2L,quantile,c(.025,.975))
  atomic_csv(data.frame(cohort_id=train$cohort_id[hidden],C=train$C[hidden],posterior_mean=imp_point,q025=iq[1,],q975=iq[2,],
    scale=if(cell=="O0")"model" else "raw"),file.path(task_out,"U_posterior_intervals.csv"))
  atomic_csv(metrics, metrics_path)
  minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  set_task(fold, cell, "DONE", paste0("Y RMSE ", round(metrics$RMSE[1], 3)),
           minutes = minutes, diag = diag)
  write_progress("RUNNING", sprintf("fold%02d %s", fold, cell), "continuation and scoring complete")
}

tryCatch({
  write_progress("RUNNING", "initializing", "selected from original diagnostics")
  order <- data.frame(fold=c(1L,3L,4L),cell=c("O0","O20","O0"))
  for (j in seq_len(nrow(order))) run_one(order$fold[j], order$cell[j])
  publish_latest_summary()
  final_state <- if (all(tasks$continued_pass_screen[tasks$selected] %in% TRUE))
    "DONE_SCREEN_PASSED" else "DONE_WITH_CONVERGENCE_WARNINGS"
  write_progress(final_state, "none", "all selected chains extended; interpret only with diagnostics")
  atomic_lines(format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
               file.path(out, "complete.flag"))
}, error = function(e) {
  atomic_lines(conditionMessage(e), file.path(out, "failed.flag"))
  write_progress("FAILED", "none", conditionMessage(e))
  stop(e)
})
