# 02: Fit -> predict -> impute, using package functions; sourcing never fits.
# Entry points used by 01_setup.R:
#   run_eiv()       : eivGP::fit_eivgp -> predict_eivgp -> impute_eivgp
#   run_baselines(): eivGP:::run_published_mixedgp_competitors, m_vec=4
#   run_one()      : eivGP::continue_eivgp, then predict and impute again
#   run_task()     : the same baseline wrapper, m_vec=c(4,2), expanded X/C
# Native calls inside that pinned eivGP wrapper:
#   UC-GP: kergp::gp / stats::predict (kGauss * q1Symm covariance)
#   LVGP : LVGP::LVGP_fit / LVGP::LVGP_predict
#   EzGP : EzGP::EzGP_fit / EzGP::EzGP_predict (training-only nugget CV)
# RMD-BEGIN/END markers let the Rmd display these exact calls without refitting.
# Numerical arguments, task seeds, masks and optimizer controls are frozen.


# A. Initial EIV and original baseline group
ocean_initial_functions <- function(envir) {
  eval(quote({
# Initial-fit functions are defined in a task-local environment.
# root/out/smoke/iter/burn/chains/cores/u_block are supplied by 01_setup.R.
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
  # Fit/predict all three published baselines with original continuous X and class C.
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
  # 1. Assemble X/C and reveal U only on the fixed training-calibration mask.
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
    # 2. Fit: missing U is NA; exact calibration entries are in u_obs.
    # RMD-BEGIN eiv-fit
    model <- eivGP::fit_eivgp(
      X = x_train, y = train$Y, C = c_train, U_obs = u_obs,
      engine = "univariate", m_vec = 4L, standardize_U = TRUE,
      kernel = "se", n_iter = iter, burn = burn, thin = 1L,
      n_chains = chains, parallel = !smoke, n_cores = cores,
      sampler_control = list(u_block_size = u_block),
      seed = 2026092680L + fold * 100L + seed_offset,
      verbose = TRUE
    )
    # RMD-END eiv-fit
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
  # 3. Predict withheld Y from (X,C), integrating over unobserved U.
  # RMD-BEGIN eiv-predict
  response_draws <- eivGP::predict_eivgp(fit, new_X = x_test, new_C = c_test,
    target = "response", draw_ids = draw_ids, n_per_draw = 1L, joint = FALSE,
    seed = 2026092700L + fold * 100L + seed_offset)
  # RMD-END eiv-predict
  point <- colMeans(as.matrix(response_draws))
  atomic_csv(data.frame(cohort_id = test$cohort_id, observed_Y = test$Y,
    predicted_Y = point), file.path(cell_out, "predictions.csv"))
  metrics <- score(point, test$Y, fold, "EIV-GP", cell)

  # Freeze nested test masks once per fold using the prespecified seed, without Y.
  set.seed(202609260L+fold)
  known20<-known50<-rep(FALSE,nrow(test))
  for (c in 1:4) {idx<-sample(which(test$C==c));known20[idx[1:2]]<-TRUE;known50[idx[1:5]]<-TRUE}
  known<-switch(cell,O0=rep(FALSE,nrow(test)),O20=known20,O50=known50)
  atomic_csv(data.frame(cohort_id=test$cohort_id,C=test$C,observed_U_O20=known20,observed_U_O50=known50,mask_seed=202609260L+fold),
    file.path(fold_out,"test_U_masks.csv"))
  mixed<-response_draws
  if (any(known)) {
    xu<-eivGP::predict_eivgp(fit,new_X=x_test[known,,drop=FALSE],new_U=matrix(test$U[known],ncol=1L),new_U_scale="raw",
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

  # 4. Impute hidden training U; O0 uses model coordinates, O20/O50 raw units.
  hidden <- which(!is.finite(u_obs[, 1L]))
  if (cell == "O0") {
    ## Without calibration, raw U is not affinely identified by the package.
    imp <- eivGP::impute_eivgp(fit, rows = hidden, draw_ids = draw_ids, scale = "model")
    imp_point <- colMeans(matrix(imp[, , 1L], nrow = length(draw_ids)))
    atomic_csv(data.frame(cohort_id = train$cohort_id[hidden], C = train$C[hidden],
      true_U_raw = train$U[hidden], EIV_U_model_mean = imp_point),
      file.path(cell_out, "imputations_model_scale.csv"))
  } else {
    # RMD-BEGIN eiv-impute
    imp <- eivGP::impute_eivgp(fit, rows = hidden, draw_ids = draw_ids, scale = "raw")
    # RMD-END eiv-impute
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


  }), envir=envir)
  invisible(envir)
}

# B. Diagnostic-selected +20k continuation
ocean_continuation_functions <- function(envir) {
  eval(quote({
# Diagnostic-selected +20k continuation; independent of unrelated folds.
# Selection uses current initial diagnostics in 01_setup.R; prediction seeds and saved RNG states are retained.
diagnostic_row <- function(fit) {
  tables <- fit$diagnostics[c("rhat_hyper", "rhat_tau", "rhat_u")]
  rhats <- unlist(lapply(tables, function(z) z$rhat))
  bulk <- unlist(lapply(tables, function(z) z$ess_bulk))
  max_rhat <- max(rhats[is.finite(rhats)])
  min_ess <- min(bulk[is.finite(bulk)])
  data.frame(max_rhat = max_rhat, min_ess_bulk = min_ess,
    pass_screen = max_rhat <= 1.01 && min_ess >= 400)
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
    # Continue the saved chains/RNG states, then repeat prediction and imputation.
    # RMD-BEGIN eiv-continue
    continued <- eivGP::continue_eivgp(original, n_iter = 20000L,
      parallel = cores > 1L, n_cores = cores, verbose = TRUE)
    # RMD-END eiv-continue
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
  response_draws <- eivGP::predict_eivgp(fit, new_X = x_test, new_C = c_test,
    target = "response", draw_ids = draw_ids, n_per_draw = 1L,
    joint = FALSE, seed = 2026092900L + fold * 100L + match(cell, c("O0", "O20", "O50")))
  point <- colMeans(as.matrix(response_draws))
  metrics <- score(point, test$Y, fold, cell, "EIV-GP")
  atomic_csv(data.frame(cohort_id = test$cohort_id, observed_Y = test$Y,
    predicted_Y = point), file.path(task_out, "predictions.csv"))

  mask <- if (cell == "O0") rep(FALSE, nrow(train)) else
    train[[paste0("calib_", tolower(cell))]] == 1L

  # Use the fixed test observation masks for both initial and continued fits.
  saved_masks<-read.csv(file.path(root,"full5_class4",sprintf("fold%02d",fold),"test_U_masks.csv"))
  stopifnot(identical(as.character(saved_masks$cohort_id),as.character(test$cohort_id)))
  known<-if(cell=="O0")rep(FALSE,nrow(test)) else as.logical(saved_masks[[paste0("observed_U_",cell)]])
  stopifnot(sum(known)==switch(cell,O0=0L,O20=8L,O50=20L))
  mixed<-response_draws
  if(any(known)) {
    xu<-eivGP::predict_eivgp(fit,new_X=x_test[known,,drop=FALSE],new_U=matrix(test$U[known],ncol=1L),new_U_scale="raw",
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
  imp <- eivGP::impute_eivgp(fit, rows = hidden, draw_ids = draw_ids,
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


  }), envir=envir)
  invisible(envir)
}

# C. Expanded native baseline fits
ocean_augmented_functions <- function(envir) {
  eval(quote({
# Native baseline fitting, prediction and scoring for expanded X/C inputs.
# Task-local inputs/output_root/manifest/controls are supplied by 01_setup.R.
subset_draws <- function(draws, idx) {
  out <- as.matrix(draws)[,idx,drop=FALSE]
  for (name in c("conditional_means", "conditional_vars")) {
    component <- attr(draws,name,exact=TRUE)
    stopifnot(!is.null(component))
    attr(out,name) <- as.matrix(component)[,idx,drop=FALSE]
  }
  out
}
score_draws <- function(draws, point, y) {
  stopifnot(all(is.finite(draws)), all(is.finite(point)), length(point)==length(y))
  q <- apply(draws,2L,quantile,c(.025,.975))
  density <- eivGP:::normal_mixture_nlpd(y, attr(draws,"conditional_means",exact=TRUE),
    attr(draws,"conditional_vars",exact=TRUE))
  stopifnot(is.finite(density$value))
  data.frame(n_test=length(y), RMSE=sqrt(mean((point-y)^2)), MAE=mean(abs(point-y)),
    CRPS=mean(eivGP:::crps_sample_matrix(draws,y)), NLPD=density$value,
    Coverage95=mean(y>=q[1,]&y<=q[2,]), Width95=mean(q[2,]-q[1,]),
    IntervalScore95=mean(eivGP:::interval_score(q[1,],q[2,],y)), NLPD_reason=density$reason)
}

run_task <- function(task) {
  # Keep oxygen class in C; add categorical R to C and continuous W to X.
  p <- prepare_task(task,inputs)
  dir <- file.path(output_root,"fits",task$task_id);dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  signature <- digest::digest(list(task=task[,c("task_id","fit_seed","placeholder")],
    X_train=p$X_train,C_train=p$C_train,y_train=p$y_train,X_test=p$X_test,C_test=p$C_test,
    manifest=manifest),algo="sha256")
  atomic_csv(data.frame(cohort_id=p$data$train$cohort_id,C=p$C_train[,1],
    observed=p$data$observed_train,W=p$X_train[,3]),file.path(dir,"encoded_training.csv"))
  atomic_csv(data.frame(cohort_id=rep(p$data$test$cohort_id,2), prediction_mode=rep(c("class_only","mixed"),each=40),
    C=p$C_test[,1],observed=p$C_test[,2]==2L,W=p$X_test[,3]),file.path(dir,"encoded_test.csv"))
  cached <- file.path(dir,"package_output.rds")
  if (file.exists(cached)) {
    saved <- readRDS(cached);stopifnot(identical(saved$signature,signature));ans<-saved$output
  } else {
    # Same native wrapper/controls as original baselines; only input columns differ.
    # RMD-BEGIN native-baselines
    ans <- eivGP:::run_published_mixedgp_competitors(
      X_train=p$X_train,y_train=p$y_train,C_train=p$C_train,
      X_test=p$X_test,C_test=p$C_test,m_vec=c(4L,2L),n_draw=400L,
      seed=as.integer(task$wrapper_seed),methods=task$method,strict=FALSE,controls=controls)
    # RMD-END native-baselines
    atomic_rds(list(signature=signature,output=ans,fit_seed=task$fit_seed),cached)
  }
  atomic_csv(ans$status,file.path(dir,"package_status.csv"))
  if (!task$method %in% names(ans$draws)) stop(ans$status$message[[1]])
  point <- as.numeric(ans$means[[task$method]])
  draws <- ans$draws[[task$method]]
  class <- subset_draws(draws,1:40); mixed <- subset_draws(draws,41:80)
  unknown <- !p$data$observed_test_mixed
  stopifnot(max(abs(point[1:40][unknown]-point[41:80][unknown]))<1e-8)
  # Share Monte Carlo draws between class-only and mixed test predictions.
  mixed[,unknown] <- class[,unknown,drop=FALSE]
  for (name in c("conditional_means","conditional_vars")) {
    component<-attr(mixed,name,exact=TRUE);component[,unknown]<-attr(class,name,exact=TRUE)[,unknown,drop=FALSE]
    attr(mixed,name)<-component
  }
  metric_rows<-prediction_rows<-list()
  for (mode in c("class_only","mixed")) {
    idx<-if(mode=="class_only")1:40 else 41:80
    dd<-if(mode=="class_only")class else mixed
    metadata<-data.frame(task_id=task$task_id,fold=task$fold,cell=task$cell,method=task$method,
      placeholder=task$placeholder,prediction_mode=mode,fit_seed=task$fit_seed,
      optimizer_status=ans$status$optimization_status[[1]])
    metric_rows[[mode]]<-cbind(metadata,score_draws(dd,point[idx],p$y_test))
    prediction_rows[[mode]]<-data.frame(metadata,cohort_id=p$data$test$cohort_id,
      observed_Y=p$y_test,predicted_Y=point[idx],observed_U=p$C_test[idx,2]==2L)
  }
  atomic_csv(do.call(rbind,metric_rows),file.path(dir,"metrics.csv"))
  atomic_csv(do.call(rbind,prediction_rows),file.path(dir,"predictions.csv"))
  atomic_rds(list(class_only=class,mixed=mixed,known=p$data$observed_test_mixed),file.path(dir,"prediction_distributions.rds"))
  list(optimizer_status=ans$status$optimization_status[[1]],detail=ans$status$warnings[[1]])
}

  }), envir=envir)
  invisible(envir)
}
