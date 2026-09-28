#!/usr/bin/env Rscript

## Frozen U<180, pressure300-1500 joint-quintile-pressure five-fold package experiment.
## No package function is modified. O0/O20/O50 and UC-GP/LVGP/EzGP are run.
options(warn = 1)
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", arg[[1L]]), mustWork = TRUE)
root <- dirname(dirname(script))
suppressPackageStartupMessages(library(eivGP))
stopifnot(as.character(packageVersion("eivGP")) == "0.3.1")

smoke <- identical(Sys.getenv("SI_O2_CLASS4_SMOKE"), "1")
out <- file.path(root, if (smoke) "smoke_class4" else "full5_class4")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(out, "complete.flag"))) quit(save="no", status=0L)
progress_path <- file.path(out, "PROGRESS.txt")
status_path <- file.path(out, "task_status.csv")
started <- Sys.time()
iter <- if (smoke) 80L else 20000L
burn <- if (smoke) 30L else 5000L
chains <- if (smoke) 1L else 4L
cores <- if (smoke) 1L else 4L
u_block <- 8L   # Installed-package default; avoids carrying over the untested 32 setting.
fold_order <- if (smoke) 1L else 1:5
task_order <- if (smoke) c("O50", "O0") else c("baselines", "O50", "O20", "O0")

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
tasks <- if (file.exists(status_path)) {
  read.csv(status_path, stringsAsFactors = FALSE)
} else {
  x <- expand.grid(fold = fold_order, cell = task_order,
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  x <- x[order(match(x$fold, fold_order), match(x$cell, task_order)), ]
  x$state <- "PENDING"
  x$updated <- ""
  x$minutes <- NA_real_
  x$detail <- ""
  x
}
## A previous launch may have begun before O0 was added. Extend its status
## without changing any existing task or output; interrupted RUNNING is retried.
expected <- expand.grid(fold = fold_order, cell = task_order,
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
missing <- expected[!paste(expected$fold, expected$cell) %in%
                      paste(tasks$fold, tasks$cell), , drop = FALSE]
if (nrow(missing)) {
  missing$state <- "PENDING"
  missing$updated <- ""
  missing$minutes <- NA_real_
  missing$detail <- "added O0 to ongoing experiment"
  tasks <- rbind(tasks, missing)
}
tasks <- tasks[order(match(tasks$fold, fold_order), match(tasks$cell, task_order)), ]
stopifnot(nrow(tasks) == length(fold_order) * length(task_order))

set_task <- function(fold, cell, state, detail = "", minutes = NA_real_) {
  i <- which(tasks$fold == fold & tasks$cell == cell)
  stopifnot(length(i) == 1L)
  tasks$state[i] <<- state
  tasks$updated[i] <<- format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
  tasks$detail[i] <<- detail
  if (is.finite(minutes)) tasks$minutes[i] <<- minutes
  atomic_csv(tasks, status_path)
}
write_progress <- function(state, current = "none", detail = "") {
  initial <- c(baselines = 5, O50 = 25, O20 = 37, O0 = 56)
  estimates <- vapply(task_order, function(cell) {
    v <- tasks$minutes[tasks$cell == cell & tasks$state == "DONE" &
                       is.finite(tasks$minutes)]
    if (length(v)) median(v) else initial[[cell]]
  }, numeric(1L))
  remaining <- sum(estimates[tasks$cell[tasks$state == "PENDING"]])
  if (state == "RUNNING" && current != "none") {
    part <- strsplit(current, " ")[[1L]]
    if (length(part) == 2L && part[2L] %in% names(estimates)) {
      remaining <- remaining + estimates[[part[2L]]]
    }
  }
  atomic_lines(c(
    paste0("state: ", state),
    paste0("updated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste0("pid: ", Sys.getpid()),
    "design: retained original 200-profile cohort; U<180; pressure300-1500; classes50/100/150; frozen joint_quintile_pressure folds",
    "purpose: test whether historical gates translate to package performance; only3/5folds pass all historical gates; exploratory",
    "fold order: 1,2,3,4,5; no fold reroll or cohort reselection",
    "launch: one-shot launchd Standard, KeepAlive=false; MCMC4cores; no console-lifetime dependency",
    "outputs: xc and mixed predictions, probability scores, baseline fits/distributions, EIV posterior draws and U imputation",
    "methods: installed eivGP 0.3.1 O0/O20/O50 plus UC-GP/LVGP/EzGP",
    paste0("MCMC: ", chains, " chains x ", iter, "; burn ", burn,
           "; u_block_size ", u_block, "; cores ", cores),
    paste0("completed_tasks: ", sum(tasks$state == "DONE"), "/", nrow(tasks)),
    paste0("current_task: ", current),
    paste0("detail: ", detail),
    paste0("elapsed_this_launch_hours: ", round(as.numeric(difftime(Sys.time(), started, units = "hours")), 2)),
    paste0("estimated_remaining_hours: ", round(remaining / 60, 2)),
    paste0("estimated_finish: ", format(Sys.time() + remaining * 60,
                                        "%Y-%m-%d %H:%M:%S %z")),
    "eta_note: provisional task-boundary estimate; updated from measured durations; see worker.out.log during a fit",
    "task_state_note: RUNNING can remain unchanged during a long MCMC fit; package iteration messages are in logs",
    "", capture.output(print(tasks, row.names = FALSE))
  ), progress_path)
}
score <- function(pred, y, fold, method, cell) {
  pred <- as.numeric(pred)
  stopifnot(length(pred) == length(y), all(is.finite(pred)))
  data.frame(fold = fold, cell = cell, method = method, n_test = length(y),
             RMSE = sqrt(mean((pred - y)^2)), MAE = mean(abs(pred - y)))
}

probability_score <- function(draws,point,y,fold,method,cell) {
  qs<-apply(as.matrix(draws),2L,quantile,c(.025,.975))
  nlpd<-eivGP:::normal_mixture_nlpd(y,attr(draws,"conditional_means",exact=TRUE),attr(draws,"conditional_vars",exact=TRUE))
  cbind(score(point,y,fold,method,cell),CRPS=mean(eivGP:::crps_sample_matrix(draws,y)),
    NLPD=nlpd$value,Coverage95=mean(y>=qs[1,]&y<=qs[2,]),Width95=mean(qs[2,]-qs[1,]),
    IntervalScore95=mean(eivGP:::interval_score(qs[1,],qs[2,],y)),NLPD_reason=nlpd$reason)
}

read_fold <- function(fold) {
  train <- read.csv(file.path(root, "data", sprintf("fold%02d_train.csv", fold)))
  test <- read.csv(file.path(root, "data", sprintf("fold%02d_test.csv", fold)))
  train$C <- as.integer(findInterval(train$U, c(50, 100, 150)) + 1L)
  test$C <- as.integer(findInterval(test$U, c(50, 100, 150)) + 1L)
  stopifnot(nrow(train) == 160L, nrow(test) == 40L,
            !length(intersect(train$profile, test$profile)),
            identical(as.integer(table(test$C)), rep(10L, 4L)),
            sum(train$calib_o20) == 32L, sum(train$calib_o50) == 80L,
            all(train$calib_o20 <= train$calib_o50))
  stopifnot(all(train$U<180),all(test$U<180),all(train$G2pressure>=300 & train$G2pressure<=1500),all(test$G2pressure>=300 & test$G2pressure<=1500))
  for (c in 1:4) stopifnot(sum(train$calib_o20[train$C==c])==8L,sum(train$calib_o50[train$C==c])==20L)
  list(train = train, test = test)
}
run_baselines <- function(fold, train, test, fold_out) {
  metrics_path <- file.path(fold_out, "baseline_metrics.csv")
  if (file.exists(metrics_path) && file.exists(file.path(fold_out,"baseline_package_outputs.rds"))) {
    previous<-read.csv(metrics_path)
    if(setequal(previous$method,c("UC-GP","LVGP","EzGP"))) return(invisible("reused output"))
  }
  fit <- eivGP:::run_published_mixedgp_competitors(
    X_train = as.matrix(train[, c("G2salinity", "log_pressure")]),
    y_train = train$Y,
    C_train = matrix(train$C, ncol = 1L, dimnames = list(NULL, "oxygen_class")),
    X_test = as.matrix(test[, c("G2salinity", "log_pressure")]),
    C_test = matrix(test$C, ncol = 1L, dimnames = list(NULL, "oxygen_class")),
    m_vec = 4L, n_draw = 400L, seed = 2026092660L + fold * 100L,
    methods = c("UC-GP", "LVGP", "EzGP"), strict = FALSE,
    controls = list(
      `UC-GP` = list(n_starts = 8L),
      LVGP = list(n_starts = 8L, max_retries = 2L, max_iter_ini = 100L,
        max_iter_lat = 20L, rescue_iter_ini = 250L,
        rescue_iter_lat = 80L, max_elapsed_seconds = 900,
        parallel = FALSE),
      EzGP = list(cv_folds = 2L, maxeval = 100L, cv_score = "nlpd")
    )
  )
  saveRDS(fit, file.path(fold_out, "baseline_package_outputs.rds"))
  atomic_csv(fit$status, file.path(fold_out, "baseline_status.csv"))
  metrics <- list()
  predictions <- list()
  for (method in c("UC-GP", "LVGP", "EzGP")) {
    if (!method %in% names(fit$means)) next
    point <- as.numeric(fit$means[[method]])
    metrics[[method]] <- probability_score(fit$draws[[method]],point,test$Y,fold,method,"baseline")
    predictions[[method]] <- data.frame(cohort_id = test$cohort_id,
      method = method, observed_Y = test$Y, predicted_Y = point)
  }
  if (length(metrics)) atomic_csv(do.call(rbind, metrics), metrics_path)
  if (length(predictions)) atomic_csv(do.call(rbind, predictions),
    file.path(fold_out, "baseline_predictions.csv"))
  failed <- fit$status$method[fit$status$status != "success"]
  if (length(failed)) stop("Baseline failures: ", paste(failed, collapse = ","))
  invisible("all three baselines fitted")
}
run_eiv <- function(fold, cell, train, test, fold_out) {
  cell_out <- file.path(fold_out, cell)
  dir.create(cell_out, showWarnings = FALSE)
  metric_path <- file.path(cell_out, "metrics.csv")
  if (file.exists(metric_path)) return(invisible("reused output"))
  x_train <- as.matrix(train[, c("G2salinity", "log_pressure")])
  x_test <- as.matrix(test[, c("G2salinity", "log_pressure")])
  c_train <- matrix(train$C, ncol = 1L, dimnames = list(NULL, "oxygen_class"))
  c_test <- matrix(test$C, ncol = 1L, dimnames = list(NULL, "oxygen_class"))
  u_obs <- matrix(NA_real_, nrow(train), 1L, dimnames = list(NULL, "oxygen"))
  mask <- if (cell == "O0") rep(FALSE, nrow(train)) else
    train[[paste0("calib_", tolower(cell))]] == 1L
  u_obs[mask, 1L] <- train$U[mask]
  seed_offset <- c(O0 = 0L, O20 = 20L, O50 = 50L)[[cell]]
  fit_path <- file.path(cell_out, "fit.rds")
  fit <- if (file.exists(fit_path)) readRDS(fit_path) else {
    model <- fit_eivgp(
      X = x_train, y = train$Y, C = c_train, U_obs = u_obs,
      engine = "univariate", m_vec = 4L, standardize_U = TRUE,
      kernel = "se", n_iter = iter, burn = burn, thin = 1L,
      n_chains = chains, parallel = !smoke, n_cores = cores,
      sampler_control = list(u_block_size = u_block),
      seed = 2026092680L + fold * 100L + seed_offset,
      verbose = TRUE
    )
    saveRDS(model, fit_path)
    capture.output(summary(model), file = file.path(cell_out, "fit_summary.txt"))
    model
  }
  diagnostic_parts <- lapply(
    fit$diagnostics[c("rhat_hyper", "rhat_tau", "rhat_u")],
    function(z) z[, c("parameter", "rhat", "ess_bulk", "ess_tail", "target_pass"),
                  drop = FALSE])
  diagnostics <- do.call(rbind, diagnostic_parts)
  atomic_csv(diagnostics, file.path(cell_out, "convergence_parameters.csv"))
  finite_rhat <- diagnostics$rhat[is.finite(diagnostics$rhat)]
  finite_ess <- diagnostics$ess_bulk[is.finite(diagnostics$ess_bulk)]
  atomic_csv(data.frame(max_finite_rhat = if (length(finite_rhat)) max(finite_rhat) else NA_real_,
                        min_finite_ess_bulk = if (length(finite_ess)) min(finite_ess) else NA_real_,
                        limited_screen_pass = length(finite_rhat) > 0L && length(finite_ess) > 0L &&
                          max(finite_rhat) <= 1.01 && min(finite_ess) >= 400),
             file.path(cell_out, "convergence_summary.csv"))
  n_draw <- nrow(fit$mcmc$samples_tau)
  stopifnot(is.finite(n_draw), n_draw > 0L)
  set.seed(2026092690L + fold * 100L + seed_offset)
  draw_ids <- sort(sample(seq_len(n_draw), min(if (smoke) 10L else 1000L, n_draw)))
  response_draws <- predict_eivgp(fit, new_X = x_test, new_C = c_test,
    target = "response", draw_ids = draw_ids, n_per_draw = 1L, joint = FALSE,
    seed = 2026092700L + fold * 100L + seed_offset)
  point <- colMeans(as.matrix(response_draws))
  atomic_csv(data.frame(cohort_id = test$cohort_id, observed_Y = test$Y,
    predicted_Y = point), file.path(cell_out, "predictions.csv"))
  metrics <- score(point, test$Y, fold, "EIV-GP", cell)

  # Freeze Y-blind nested test masks once per fold, with the same rule as the reference review.
  set.seed(202609260L+fold)
  known20<-known50<-rep(FALSE,nrow(test))
  for (c in 1:4) {idx<-sample(which(test$C==c));known20[idx[1:2]]<-TRUE;known50[idx[1:5]]<-TRUE}
  known<-switch(cell,O0=rep(FALSE,nrow(test)),O20=known20,O50=known50)
  atomic_csv(data.frame(cohort_id=test$cohort_id,C=test$C,observed_U_O20=known20,observed_U_O50=known50,mask_seed=202609260L+fold),
    file.path(fold_out,"test_U_masks.csv"))
  mixed<-response_draws
  if (any(known)) {
    xu<-predict_eivgp(fit,new_X=x_test[known,,drop=FALSE],new_U=matrix(test$U[known],ncol=1L),new_U_scale="raw",
      target="response",draw_ids=draw_ids,n_per_draw=1L,joint=FALSE,seed=2026093000L+fold*100L+match(cell,c("O0","O20","O50")))
    mixed[,known]<-xu
    for (a in c("conditional_means","conditional_vars")) {
      comp<-attr(response_draws,a,exact=TRUE);comp[,known]<-attr(xu,a,exact=TRUE);attr(mixed,a)<-comp
    }
  }
  stopifnot(identical(mixed[,!known,drop=FALSE],response_draws[,!known,drop=FALSE]))
  saveRDS(list(xc=response_draws,mixed=mixed,known=known,draw_ids=draw_ids),
    file.path(cell_out,"prediction_distributions.rds"))
  atomic_csv(data.frame(cohort_id=test$cohort_id,observed_Y=test$Y,predicted_Y=colMeans(mixed),test_U_observed=known),
    file.path(cell_out,"mixed_predictions.csv"))
  atomic_csv(rbind(probability_score(response_draws,point,test$Y,fold,"EIV-GP xc",cell),
    probability_score(mixed,colMeans(mixed),test$Y,fold,"EIV-GP mixed",cell)),file.path(cell_out,"predictive_metrics.csv"))

  hidden <- which(!is.finite(u_obs[, 1L]))
  if (cell == "O0") {
    ## Without calibration, raw U is not affinely identified by the package.
    imp <- impute_eivgp(fit, rows = hidden, draw_ids = draw_ids, scale = "model")
    imp_point <- colMeans(matrix(imp[, , 1L], nrow = length(draw_ids)))
    atomic_csv(data.frame(cohort_id = train$cohort_id[hidden], C = train$C[hidden],
      true_U_raw = train$U[hidden], EIV_U_model_mean = imp_point),
      file.path(cell_out, "imputations_model_scale.csv"))
  } else {
    imp <- impute_eivgp(fit, rows = hidden, draw_ids = draw_ids, scale = "raw")
    imp_point <- colMeans(matrix(imp[, , 1L], nrow = length(draw_ids)))
    calib_means <- tapply(train$U[mask], train$C[mask], mean)
    imp_base <- as.numeric(calib_means[as.character(train$C[hidden])])
    atomic_csv(data.frame(cohort_id = train$cohort_id[hidden], C = train$C[hidden],
      true_U = train$U[hidden], EIV_U_mean = imp_point,
      class_calibration_mean = imp_base), file.path(cell_out, "imputations.csv"))
    metrics <- rbind(metrics,
      score(imp_point, train$U[hidden], fold, "EIV U imputation", cell),
      score(imp_base, train$U[hidden], fold, "calibration class mean", cell))
  }
  iq<-apply(matrix(imp[,,1L],nrow=length(draw_ids)),2L,quantile,c(.025,.975))
  atomic_csv(data.frame(cohort_id=train$cohort_id[hidden],C=train$C[hidden],posterior_mean=imp_point,q025=iq[1,],q975=iq[2,],scale=if(cell=="O0")"model" else "raw"),file.path(cell_out,"U_posterior_intervals.csv"))
  atomic_csv(metrics, metric_path)
  tau <- as.matrix(fit$mcmc$samples_tau)
  atomic_csv(data.frame(boundary = seq_len(ncol(tau)),
    q025_model_scale = apply(tau, 2L, quantile, .025),
    median_model_scale = apply(tau, 2L, median),
    q975_model_scale = apply(tau, 2L, quantile, .975)),
    file.path(cell_out, "threshold_posterior_model_scale.csv"))
  invisible(paste0("Y RMSE ", round(metrics$RMSE[1L], 3)))
}

write_progress("RUNNING", "initializing", "starting/resuming exploratory five folds")
for (fold in fold_order) {
  d <- read_fold(fold)
  fold_out <- file.path(out, sprintf("fold%02d", fold))
  dir.create(fold_out, showWarnings = FALSE)
  for (cell in task_order) {
    i <- which(tasks$fold == fold & tasks$cell == cell)
    if (tasks$state[i] == "DONE") next
    set_task(fold, cell, "RUNNING", "package fit/prediction in progress")
    write_progress("RUNNING", sprintf("fold%02d %s", fold, cell), "task started")
    t0 <- Sys.time()
    outcome <- tryCatch(
      if (cell == "baselines") run_baselines(fold, d$train, d$test, fold_out)
      else run_eiv(fold, cell, d$train, d$test, fold_out),
      error = function(e) structure(conditionMessage(e), class = "task_error")
    )
    minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    if (inherits(outcome, "task_error")) {
      set_task(fold, cell, "ERROR", as.character(outcome), minutes)
    } else {
      set_task(fold, cell, "DONE", as.character(outcome), minutes)
    }
    write_progress("RUNNING", sprintf("fold%02d %s", fold, cell),
                   if (inherits(outcome, "task_error")) paste("error:", outcome) else "task completed")
  }
}
final_state <- if (all(tasks$state == "DONE")) "DONE" else "PARTIAL"
write_progress(final_state, "none", if (final_state == "DONE")
  if (smoke) "fold 5 O50/O0 smoke tests complete" else "all five folds complete"
  else "one or more tasks errored; inspect status and logs")
if(final_state=="DONE") {
  atomic_lines(format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),file.path(out,"complete.flag"))
} else {
  atomic_lines("One or more tasks errored; inspect task_status.csv",file.path(out,"failed.flag"))
}
