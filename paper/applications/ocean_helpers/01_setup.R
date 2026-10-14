# 01: Configuration, fixed inputs, task selection, execution and provenance.

# 1. Numerical settings are fixed; selectors only choose work, never alter a seed.
ocean_controls <- function() list(
  `UC-GP` = list(n_starts = 8L),
  LVGP = list(n_starts = 8L, max_retries = 2L, max_iter_ini = 100L,
    max_iter_lat = 20L, rescue_iter_ini = 250L, rescue_iter_lat = 80L,
    max_elapsed_seconds = 900, parallel = FALSE),
  EzGP = list(cv_folds = 2L, maxeval = 100L, cv_score = "nlpd"))

ocean_settings <- function(profile = "paper", action = "plan", folds = 1:5,
  cells = c("O0", "O20", "O50"),
  modules = c("eiv", "original_baselines", "continuation", "complete_case_baselines",
    "prediction", "diagnostics", "surface_comparison"),
  task_ids = character(), cores = 4L, run_id = "local") {
  profile <- match.arg(profile, c("paper", "quick", "pilot"))
  action <- match.arg(action, c("plan", "run", "report"))
  folds <- as.integer(unlist(folds)); cells <- as.character(unlist(cells))
  modules <- as.character(unlist(modules))
  task_ids <- as.character(unlist(task_ids))
  task_ids <- task_ids[nzchar(task_ids)]
  stopifnot(length(folds) > 0L, !anyNA(folds), all(folds %in% 1:5),
    length(cells) > 0L, all(cells %in% c("O0", "O20", "O50")),
    all(modules %in% c("eiv", "original_baselines", "continuation", "complete_case_baselines",
      "prediction", "diagnostics", "surface_comparison")),
    length(cores) == 1L,
    is.finite(cores), cores >= 1L, cores == as.integer(cores),
    length(run_id) == 1L, grepl("^[A-Za-z0-9_-]+$", run_id))
  if (profile == "quick") {
    folds <- intersect(folds, 1L); cells <- intersect(cells, c("O0", "O20", "O50"))
    modules <- intersect(modules, c("eiv", "complete_case_baselines", "prediction"))
    if (!length(folds) || !length(cells)) stop("Quick supports fold 1 only.")
  }
  list(profile = profile, action = action, folds = unique(folds), cells = unique(cells),
    modules = unique(modules), task_ids = unique(task_ids), cores = as.integer(cores),
    run_id = run_id, iter = if (profile == "paper") 20000L else 80L,
    burn = if (profile == "paper") 5000L else 30L,
    chains = if (profile != "quick") 4L else 1L, u_block = 8L)
}

# 2. A complete, stable task catalogue, irrespective of the device's subset.
ocean_task_catalogue <- function(profile = "paper") {
  rows <- list(); add <- function(id, module, fold, cell, method = "EIV-GP", seed = NA_integer_) {
    rows[[length(rows) + 1L]] <<- data.frame(task_id = id, module, fold, cell,
      method, fit_seed = seed, stringsAsFactors = FALSE)
  }
  for (f in if (profile != "quick") 1:5 else 1L) {
    for (cell in if (profile != "quick") c("O50", "O20", "O0") else c("O50", "O20", "O0"))
      add(sprintf("eiv_fold%02d_%s", f, cell), "eiv", f, cell,
        seed = 2026092680L + 100L*f + c(O0=0L, O20=20L, O50=50L)[[cell]])
    if (profile != "quick") add(sprintf("original_fold%02d", f), "original_baselines", f,
      "baseline", "UC-GP+LVGP+EzGP", seed = 2026092660L + 100L*f)
  }
  for(f in if(profile!="quick") 1:5 else 1L) for(cell in c("O20","O50"))
    for(method in c("UC-GP","LVGP","EzGP"))
      add(sprintf("cc_fold%02d_%s_%s",f,cell,gsub("-","",method)),
        "complete_case_baselines",f,cell,method,
        seed=2026101400L+100L*f+as.integer(sub("O","",cell))+match(method,c("UC-GP","LVGP","EzGP")))
  if (profile != "quick") {
    # Every EIV cell is a candidate; current initial diagnostics decide execution.
    for (f in 1:5) for (cell in c("O0","O20","O50"))
      add(sprintf("continue_fold%02d_%s", f, cell), "continuation", f, cell)
  }
  do.call(rbind, rows)
}

ocean_select_tasks <- function(settings) {
  all <- ocean_task_catalogue(settings$profile)
  chosen <- all[all$module %in% settings$modules & all$fold %in% settings$folds &
    (all$cell %in% settings$cells | all$cell == "baseline"), , drop = FALSE]
  if (length(settings$task_ids)) {
    if (length(setdiff(settings$task_ids, chosen$task_id)))
      stop("task_ids not included by the other selectors: ", paste(setdiff(settings$task_ids, chosen$task_id), collapse=", "))
    chosen <- chosen[chosen$task_id %in% settings$task_ids, , drop=FALSE]
  }
  rownames(chosen) <- NULL
  chosen
}

# Adaptive continuation: the rule is fixed before fitting, the cells are not.
# Finite Rhat and bulk ESS determine this limited screen. Undefined
# diagnostics (e.g. constant dictionary occupancy) and tail ESS are also
# recorded for inspection.
ocean_continuation_policy <- function(profile="paper") list(max_rhat=1.01,min_ess_bulk=400,
  n_iter=if(profile=="pilot")80L else 20000L,max_rounds=1L,screen="finite Rhat and bulk ESS")

ocean_read_screen <- function(path) {
  p <- ocean_continuation_policy()
  if(!file.exists(path)) return(data.frame(ready=FALSE,max_rhat=NA_real_,min_ess_bulk=NA_real_,
    n_nonfinite_rhat=NA_integer_,n_nonfinite_ess_bulk=NA_integer_,n_nonfinite_ess_tail=NA_integer_,
    screen_pass=NA,diagnostic_md5=NA_character_))
  z <- read.csv(path);stopifnot(all(c("rhat","ess_bulk","ess_tail") %in% names(z)))
  r <- z$rhat[is.finite(z$rhat)];b <- z$ess_bulk[is.finite(z$ess_bulk)]
  mx <- if(length(r)) max(r) else NA_real_;mn <- if(length(b)) min(b) else NA_real_
  data.frame(ready=TRUE,max_rhat=mx,min_ess_bulk=mn,
    n_nonfinite_rhat=sum(!is.finite(z$rhat)),n_nonfinite_ess_bulk=sum(!is.finite(z$ess_bulk)),
    n_nonfinite_ess_tail=sum(!is.finite(z$ess_tail)),
    screen_pass=length(r)>0L && length(b)>0L && mx<=p$max_rhat && mn>=p$min_ess_bulk,
    diagnostic_md5=unname(tools::md5sum(path)))
}

ocean_continuation_decisions <- function(output_dir,profile="paper",write=FALSE) {
  rows<-list()
  for(f in if(profile!="quick")1:5 else 1L) for(cell in if(profile!="quick")c("O0","O20","O50") else c("O0","O20","O50")) {
    path<-file.path(if(profile!="quick")"full5_class4" else "smoke_class4",sprintf("fold%02d",f),cell,"convergence_parameters.csv")
    screen<-ocean_read_screen(file.path(output_dir,path))
    needs<-if(profile=="quick")FALSE else if(screen$ready)!screen$screen_pass else NA
    rows[[length(rows)+1L]]<-cbind(data.frame(fold=f,cell,initial_diagnostics=path),screen,
      needs_continuation=needs,decision=if(is.na(needs))"awaiting_initial" else if(needs)paste0("append_",ocean_continuation_policy(profile)$n_iter) else "keep_initial")
  }
  d<-do.call(rbind,rows)
  if(write)ocean_atomic_csv(d,file.path(output_dir,"continuation_decisions.csv"))
  d
}

ocean_final_eiv_sources <- function(output_dir,profile="paper") {
  decisions<-ocean_continuation_decisions(output_dir,profile)
  rows<-lapply(seq_len(nrow(decisions)),function(i) {
    d<-decisions[i,];continued<-isTRUE(d$needs_continuation)
    source<-file.path(if(continued)"continuation_plus20k" else if(profile!="quick")"full5_class4" else "smoke_class4",
      sprintf("fold%02d",d$fold),d$cell)
    final<-ocean_read_screen(file.path(output_dir,source,"convergence_parameters.csv"))
    data.frame(fold=d$fold,cell=d$cell,source,continued,
      fit_file=if(continued)"fit_plus20k.rds" else "fit.rds",
      available=d$ready && file.exists(file.path(output_dir,source,"predictive_metrics.csv")),
      initial_screen_pass=d$screen_pass,final_screen_pass=final$screen_pass)
  })
  do.call(rbind,rows)
}

# Execution and provenance
# Load definitions only. Rmd and 01_setup.R both enter through ocean_execute().
ocean_source_modules <- function(repository_dir, envir = parent.frame()) {
  app <- file.path(repository_dir, "paper", "applications", "ocean_helpers")
  for (name in c("02_fit_models.R", "03_reporting.R", "04_diagnostics.R"))
    sys.source(file.path(app, name), envir)
  invisible(app)
}

ocean_atomic_rds <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  tmp <- tempfile(".atomic-", tmpdir=dirname(path)); saveRDS(x,tmp)
  if (!file.rename(tmp,path)) stop("Cannot publish ",path)
}
ocean_atomic_csv <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  tmp <- tempfile(".atomic-", tmpdir=dirname(path)); write.csv(x,tmp,row.names=FALSE,na="")
  if (!file.rename(tmp,path)) stop("Cannot publish ",path)
}
ocean_atomic_lines <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  tmp <- tempfile(".atomic-", tmpdir=dirname(path)); writeLines(x,tmp)
  if (!file.rename(tmp,path)) stop("Cannot publish ",path)
}
ocean_hashes <- function(paths) unname(tools::md5sum(paths))
ocean_resolve_path <- function(path,repository_dir) {
  path <- path.expand(as.character(path))
  if(!length(path)) return(character())
  absolute <- grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)",path)
  ifelse(absolute,path,file.path(repository_dir,path))
}

# Frozen design and test masks. This never draws a new cohort or calibration set.
ocean_read_inputs <- function(data_dir) {
  cohort <- read.csv(file.path(data_dir,"cohort_manifest.csv"))
  stopifnot(nrow(cohort)==200L,!anyDuplicated(cohort$profile),!anyDuplicated(cohort$cohort_id))
  cohort$oxygen_class <- factor(findInterval(cohort$U,c(50,100,150))+1L,
    levels=1:4, labels=c("<50","50–100","100–150","150–180"))
  cohort
}
ocean_prepare_inputs <- function(data_dir) {
  cohort <- ocean_read_inputs(data_dir); folds <- list(); held <- character(); audit <- list()
  for (f in 1:5) {
    tr <- read.csv(file.path(data_dir,sprintf("fold%02d_train.csv",f)))
    te <- read.csv(file.path(data_dir,sprintf("fold%02d_test.csv",f)))
    tr$C <- findInterval(tr$U,c(50,100,150))+1L; te$C <- findInterval(te$U,c(50,100,150))+1L
    tr$calib_o20 <- as.logical(tr$calib_o20); tr$calib_o50 <- as.logical(tr$calib_o50)
    stopifnot(nrow(tr)==160L,nrow(te)==40L,!length(intersect(tr$profile,te$profile)),
      setequal(c(tr$cohort_id,te$cohort_id),cohort$cohort_id),
      !anyDuplicated(c(tr$cohort_id,te$cohort_id)),
      all(table(factor(tr$C,levels=1:4))==40L),all(table(factor(te$C,levels=1:4))==10L),
      all(tr$calib_o20<=tr$calib_o50),
      all(table(factor(tr$C[tr$calib_o20],levels=1:4))==8L),
      all(table(factor(tr$C[tr$calib_o50],levels=1:4))==20L))
    columns <- c("cohort_id","profile","G2salinity","G2pressure","log_pressure","U","Y")
    for (d in list(tr,te)) {
      ref <- cohort[match(d$cohort_id,cohort$cohort_id),columns]
      stopifnot(isTRUE(all.equal(d[,columns],ref,check.attributes=FALSE)),
        all(d$U>0 & d$U<180),all(d$G2pressure>=300 & d$G2pressure<=1500),
        max(abs(d$log_pressure-log(d$G2pressure)))<1e-10)
    }
    set.seed(202609260L+f); k20 <- k50 <- rep(FALSE,40L)
    for (c in 1:4) {idx <- sample(which(te$C==c)); k20[idx[1:2]]<-TRUE;k50[idx[1:5]]<-TRUE}
    masks <- data.frame(cohort_id=te$cohort_id,C=te$C,observed_U_O20=k20,
      observed_U_O50=k50,mask_seed=202609260L+f)
    # The same seed rule is independently checked against bundled saved masks.
    saved <- file.path(data_dir,sprintf("fold%02d_test_U_masks.csv",f))
    if (file.exists(saved)) stopifnot(isTRUE(all.equal(read.csv(saved),masks,check.attributes=FALSE)))
    for (cell in c("O20","O50")) {
      obs <- tr[[paste0("calib_",tolower(cell))]]; known <- if(cell=="O20") k20 else k50
      center <- mean(tr$U[obs]); scale <- sd(tr$U[obs]); stopifnot(is.finite(scale),scale>0)
      folds[[paste(f,cell,sep="/")]] <- list(train=tr,test=te,observed_train=obs,
        observed_test_mixed=known,center=center,scale=scale,masks=masks)
      audit[[length(audit)+1L]] <- data.frame(fold=f,cell,n_train=160,n_test=40,
        n_train_exact=sum(obs),n_test_exact=sum(known),oxygen_center=center,oxygen_sd=scale)
    }
    held <- c(held,te$cohort_id)
  }
  stopifnot(length(held)==200L,!anyDuplicated(held))
  inputs <- list(folds=folds,cohort=cohort,response_sd=sd(cohort$Y),audit=do.call(rbind,audit))
  inputs
}

# Portable identity excludes absolute paths, run_id, selectors and device cores.
ocean_identity <- function(repository_dir,data_dir,settings) {
  packages <- c("eivGP","kergp","LVGP","EzGP")
  versions <- setNames(vapply(packages,function(p) as.character(utils::packageVersion(p)),""),packages)
  if (versions[["eivGP"]] != "0.4.0") stop("Install the pinned upstream eivGP 0.4.0 package; see ocean_helpers/README.md.")
  names <- c("cohort_manifest.csv",sprintf("fold%02d_train.csv",1:5),sprintf("fold%02d_test.csv",1:5))
  app <- file.path(repository_dir,"paper/applications/ocean_helpers")
  pin<-read.csv(file.path(app,"UPSTREAM_PACKAGE.csv"),stringsAsFactors=FALSE)
  installed<-read.csv(system.file("module-manifest.csv",package="eivGP"),stringsAsFactors=FALSE)
  stopifnot(setequal(pin$module,installed$module),identical(pin$md5,installed$md5[match(pin$module,installed$module)]))
  code <- c("01_setup.R","02_fit_models.R")
  functions <- c("fit_eivgp","continue_eivgp","predict_eivgp","impute_eivgp",
    "run_published_mixedgp_competitors","mixedgp_adapter_ucgp","mixedgp_adapter_lvgp","mixedgp_adapter_ezgp")
  list(schema="ocean-helpers-v040-CC-v1",profile=settings$profile,upstream_commit=unique(pin$commit),
    data=setNames(ocean_hashes(file.path(data_dir,names)),names),
    code=setNames(ocean_hashes(file.path(app,code)),code),packages=versions,
    package_functions=setNames(vapply(functions,function(n) digest::digest(
      paste(deparse(get(n,asNamespace("eivGP"))),collapse="\n"),algo="sha256"),""),functions),
    numerics=list(iter=settings$iter,burn=settings$burn,chains=settings$chains,
      u_block=8L,controls=ocean_cc_controls(settings$profile),cc_controls=ocean_cc_controls(settings$profile),n_draw=400L,continuation_policy=ocean_continuation_policy(settings$profile)),
    RNG=c("Mersenne-Twister","Inversion","Rejection"))
}
ocean_task_files <- function(task,profile="paper") {
  initial <- if(profile!="quick") "full5_class4" else "smoke_class4"
  if(task$module=="original_baselines") return(file.path(initial,sprintf("fold%02d",task$fold),
    c("baseline_metrics.csv","baseline_package_outputs.rds","baseline_status.csv","baseline_predictions.csv")))
  if(task$module=="complete_case_baselines") return(file.path("complete_case_baselines/fits",task$task_id,
    c("package_output.rds","metrics.csv","predictions.csv","package_status.csv","design.rds")))
  sub <- file.path(if(task$module=="continuation") "continuation_plus20k" else initial,
    sprintf("fold%02d",task$fold),task$cell)
  file.path(sub,c(if(task$module=="continuation") "fit_plus20k.rds" else "fit.rds",
    "metrics.csv","predictive_metrics.csv","predictions.csv","mixed_predictions.csv",
    "prediction_distributions.rds","U_posterior_intervals.csv","convergence_parameters.csv",
    "convergence_summary.csv",if(task$cell=="O0") "imputations_model_scale.csv" else "imputations.csv"))
}
ocean_verify_receipt <- function(output_dir,id,identity=NULL) {
  path <- file.path(output_dir,"workflow_tasks",paste0(id,".done.rds"))
  if(!file.exists(path)) return(FALSE)
  receipt <- readRDS(path)
  if(!is.null(identity) && !identical(receipt$identity,identity)) stop("Task identity mismatch: ",id)
  files <- file.path(output_dir,receipt$files)
  if(!all(file.exists(files)) || !identical(ocean_hashes(files),receipt$hashes))
    stop("Incomplete or changed completed task: ",id)
  TRUE
}

# A single task-local environment prevents initial/continuation helpers colliding.
ocean_run_task <- function(task,repository_dir,output_dir,settings,inputs,identity) {
  env <- new.env(parent=environment()); app <- file.path(repository_dir,"paper/applications/ocean_helpers")
  env$root <- output_dir; env$smoke <- settings$profile=="quick"
  env$out <- file.path(output_dir,if(env$smoke) "smoke_class4" else "full5_class4")
  env$iter <- settings$iter; env$burn <- settings$burn; env$chains <- settings$chains
  env$cores <- if(env$smoke) 1L else settings$cores; env$u_block <- 8L
  env$baseline_controls <- ocean_cc_controls(settings$profile)
  env$continuation_iter <- ocean_continuation_policy(settings$profile)$n_iter
  env$atomic_csv <- ocean_atomic_csv; env$atomic_lines <- ocean_atomic_lines; env$atomic_rds <- ocean_atomic_rds
  env$set_task <- function(...) invisible(NULL)
  env$write_progress <- function(...) {cat(paste(...,collapse=" "),"\n");flush.console()}
  if(task$module=="complete_case_baselines") {
    ocean_run_complete_case(task,output_dir,settings,inputs)
  } else if(task$module %in% c("eiv","original_baselines")) {
    ocean_initial_functions(env)
    d <- env$read_fold(task$fold); fold_out <- file.path(env$out,sprintf("fold%02d",task$fold))
    dir.create(fold_out,recursive=TRUE,showWarnings=FALSE)
    if(task$module=="eiv") env$run_eiv(task$fold,task$cell,d$train,d$test,fold_out)
    else env$run_baselines(task$fold,d$train,d$test,fold_out)
  } else if(task$module=="continuation") {
    dependency <- sprintf("eiv_fold%02d_%s",task$fold,task$cell)
    if(!ocean_verify_receipt(output_dir,dependency,identity)) stop("Missing initial task: ",dependency)
    env$out <- file.path(output_dir,"continuation_plus20k")
    env$original_path <- function(fold,cell) file.path(output_dir,"full5_class4",sprintf("fold%02d",fold),cell,"fit.rds")
    ocean_continuation_functions(env); env$run_one(task$fold,task$cell)
  } else stop("Unknown application task module: ",task$module)
  invisible(NULL)
}

# Unified entry; report never invokes a fit. All cached computations are validated.
ocean_execute <- function(repository_dir,data_dir,output_dir,settings,
  import_dirs=character()) {
  RNGkind("Mersenne-Twister","Inversion","Rejection")
  import_dirs <- as.character(unlist(import_dirs)); import_dirs <- import_dirs[nzchar(import_dirs)]
  inputs <- ocean_prepare_inputs(data_dir); selected <- ocean_select_tasks(settings)
  if(settings$action=="plan") return(list(tasks=selected,inputs=inputs,settings=settings))
  if(settings$action=="run" || length(import_dirs)) {
    suppressPackageStartupMessages(library(eivGP))
    identity <- ocean_identity(repository_dir,data_dir,settings)
    dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
    lock <- file.path(output_dir,".workflow.lock")
    if(!dir.create(lock,showWarnings=FALSE)) stop("Output directory is locked: ",lock,
      ". Confirm the previous process has stopped before removing a stale lock.")
    on.exit(unlink(lock,recursive=TRUE),add=TRUE)
    idfile <- file.path(output_dir,"workflow_identity.rds")
    if(file.exists(idfile)) {
      if(!identical(readRDS(idfile),identity)) stop("Changed inputs/code/package/settings; choose a new output_dir.")
    } else {
      if(length(setdiff(list.files(output_dir,all.files=TRUE,no..=TRUE),c(".workflow.lock","jobs"))))
        stop("Existing outputs lack distributed provenance. Read them with action=report, or use a fresh directory.")
      ocean_atomic_rds(identity,idfile)
    }
    dir.create(file.path(output_dir,"data"),showWarnings=FALSE)
    csv <- list.files(data_dir,"[.]csv$",full.names=TRUE)
    for(f in csv) {
      to <- file.path(output_dir,"data",basename(f))
      if(file.exists(to)) stopifnot(identical(ocean_hashes(f),ocean_hashes(to))) else stopifnot(file.copy(f,to))
    }
    for(f in 1:5) {
      path <- file.path(output_dir,if(settings$profile!="quick") "full5_class4" else "smoke_class4",sprintf("fold%02d",f),"test_U_masks.csv")
      ocean_atomic_csv(inputs$folds[[paste(f,"O20",sep="/")]]$masks,path)
    }
    ocean_atomic_rds(settings,file.path(output_dir,"execution_settings.rds"))
    ocean_combine_results(import_dirs,output_dir,identity)
    job <- file.path(output_dir,"jobs",settings$run_id)
    status <- if(settings$action=="run") selected else selected[FALSE,,drop=FALSE]
    status$state <- rep("PENDING",nrow(status)); status$detail <- rep("",nrow(status))
    status$elapsed_seconds <- rep(NA_real_,nrow(status))
    publish <- function(state,current="none") {
      ocean_atomic_csv(status,file.path(job,"task_status.csv"))
      ocean_atomic_lines(c(paste("state:",state),paste("updated:",format(Sys.time(),"%Y-%m-%d %H:%M:%S %z")),
        paste("pid:",Sys.getpid()),paste("current:",current),
        sprintf("executed done: %d/%d; failed: %d; conditional skips: %d",sum(status$state %in% c("DONE","REUSED")),nrow(status)-sum(status$state=="SKIPPED_CONVERGED"),sum(status$state=="FAILED"),sum(status$state=="SKIPPED_CONVERGED")),
        "Unselected tasks are not counted as completed. Fits are serial; EIV chains may use multiple cores."),
        file.path(job,"PROGRESS.txt"))
    }
    ocean_atomic_csv(inputs$audit,file.path(job,"design_audit.csv"))
    capture.output(utils::sessionInfo(),file=file.path(job,"sessionInfo.txt"))
    for(i in if(settings$action=="run") seq_len(nrow(selected)) else integer()) {
      task <- selected[i,,drop=FALSE]; id <- task$task_id
      rownames(task) <- NULL # Assignment position is not part of task identity.
      if(task$module=="continuation") {
        dependency<-sprintf("eiv_fold%02d_%s",task$fold,task$cell)
        if(!ocean_verify_receipt(output_dir,dependency,identity)) {
          status$state[i]<-"FAILED";status$detail[i]<-paste("Missing verified initial task:",dependency)
          publish("RUNNING",id);next
        }
        decisions<-ocean_continuation_decisions(output_dir,settings$profile,write=TRUE)
        d<-decisions[decisions$fold==task$fold & decisions$cell==task$cell,]
        if(nrow(d)==1L && isTRUE(d$ready) && !isTRUE(d$needs_continuation)) {
          status$state[i]<-"SKIPPED_CONVERGED"
          status$detail[i]<-"Initial finite-Rhat/bulk-ESS screen passed; no continuation."
          publish("RUNNING",id);next
        }
      }
      if(ocean_verify_receipt(output_dir,id,identity)) {status$state[i]<-"REUSED";publish("RUNNING",id);next}
      cfg <- file.path(output_dir,"workflow_tasks",paste0(id,".config.rds"))
      if(file.exists(cfg) && !identical(readRDS(cfg),list(identity=identity,task=task))) stop("Partial-task identity mismatch: ",id)
      ocean_atomic_rds(list(identity=identity,task=task),cfg)
      status$state[i]<-"RUNNING";publish("RUNNING",id)
      t0 <- proc.time()[3]
      outcome <- tryCatch({
        ocean_run_task(task,repository_dir,output_dir,settings,inputs,identity)
        files <- ocean_task_files(task,settings$profile)
        stopifnot(all(file.exists(file.path(output_dir,files))))
        ocean_atomic_rds(list(identity=identity,task=task,files=files,
          hashes=ocean_hashes(file.path(output_dir,files)),finished=as.character(Sys.time())),
          file.path(output_dir,"workflow_tasks",paste0(id,".done.rds")))
        NULL
      },error=function(e)e)
      status$elapsed_seconds[i]<-as.numeric(proc.time()[3]-t0)
      status$state[i]<-if(inherits(outcome,"error")) "FAILED" else "DONE"
      status$detail[i]<-if(inherits(outcome,"error")) conditionMessage(outcome) else ""
      publish("RUNNING",id)
    }
    publish(if(any(status$state=="FAILED")) "SELECTED_COMPLETE_WITH_FAILURES" else "SELECTED_COMPLETE")
    ocean_atomic_lines(capture.output(print(status,row.names=FALSE)),file.path(job,"SUMMARY.txt"))
    if(any(status$state=="FAILED")) stop("Some selected tasks failed; inspect ",file.path(job,"task_status.csv"))
  }
  decisions<-ocean_continuation_decisions(output_dir,settings$profile,write=dir.exists(output_dir))
  result <- if("prediction" %in% settings$modules) ocean_prediction_report(output_dir,inputs,
    settings$profile) else NULL
  list(tasks=selected,inputs=inputs,settings=settings,prediction=result,output_dir=output_dir,continuation_decisions=decisions)
}

# Read numerical metadata without refitting or changing saved results.
ocean_saved_profile <- function(output_dir) {
  p<-file.path(output_dir,"workflow_identity.rds")
  if(file.exists(p))readRDS(p)$profile else "paper"
}
ocean_continuation_label <- function(output_dir)
  paste0("+",format(ocean_continuation_policy(ocean_saved_profile(output_dir))$n_iter,big.mark=",",scientific=FALSE))


# Optional direct execution: Rscript .../01_setup.R workflow_settings.rds
.ocean_setup_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(sys.nframe()==0L && length(.ocean_setup_arg) &&
   basename(sub("^--file=","",.ocean_setup_arg[[1]]))=="01_setup.R") {
  entry <- normalizePath(sub("^--file=","",.ocean_setup_arg[[1]]),mustWork=TRUE)
  repository_dir <- normalizePath(file.path(dirname(entry),"../../.."),mustWork=TRUE)
  ocean_source_modules(repository_dir,globalenv())
  args <- commandArgs(trailingOnly=TRUE)
  p <- if(length(args)) readRDS(args[[1]]) else list(action="plan")
  settings <- do.call(ocean_settings,p[intersect(names(p),names(formals(ocean_settings)))])
  data_dir <- if(!is.null(p$data_dir) && nzchar(p$data_dir)) ocean_resolve_path(p$data_dir,repository_dir) else file.path(repository_dir,"data/ocean_silicate")
  output_dir <- if(!is.null(p$output_dir) && nzchar(p$output_dir)) ocean_resolve_path(p$output_dir,repository_dir) else file.path(repository_dir,"reproduction/ocean",settings$profile)
  imports <- as.character(unlist(p$import_dirs)); imports <- ocean_resolve_path(imports[nzchar(imports)],repository_dir)
  result <- ocean_execute(repository_dir,data_dir,output_dir,settings,imports)
  if(settings$action=="plan") print(result$tasks,row.names=FALSE)
  if(!is.null(result$prediction)) print(result$prediction$summary,row.names=FALSE)
  if(settings$action!="plan" && settings$profile!="quick" &&
     any(c("diagnostics","surface_comparison") %in% settings$modules)) {
    dir <- if(!is.null(p$diagnostics_dir) && nzchar(p$diagnostics_dir)) ocean_resolve_path(p$diagnostics_dir,repository_dir) else file.path(output_dir,"diagnostics")
    ocean_run_diagnostics(output_dir,data_dir,dir,
      build_prediction="prediction" %in% settings$modules,build_surfaces="surface_comparison" %in% settings$modules)
  }
}
