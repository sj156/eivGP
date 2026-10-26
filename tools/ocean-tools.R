# Ocean application functions; shared numerics and utilities remain in eivGP.
.ocean_repository <- local({
  files<-as.character(unlist(lapply(sys.frames(),function(frame) {
    c(if(is.character(frame$ofile))frame$ofile,if(is.character(frame$file))frame$file)
  }),use.names=FALSE))
  files<-files[basename(files)=='ocean-tools.R']
  path<-if(length(files))tail(files,1L)else sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
  normalizePath(file.path(dirname(normalizePath(path)),'..'))
})
ocean_workflow_api <- function(repository_dir=.ocean_repository) {
  lib<-Sys.getenv('EIVGP_R_LIBRARY',file.path(repository_dir,'.r-library'))
  if(nzchar(lib)&&dir.exists(lib)).libPaths(c(normalizePath(lib),.libPaths()))
  if(!requireNamespace('eivGP',quietly=TRUE))stop('Install eivGP before running Ocean.')
  ns<-asNamespace('eivGP');e<-new.env(parent=ns)
  evalq({
ocean_settings <- function(profile = "paper", action = "plan", folds = 1:5,
  cells = c("O0", "O20", "O50"),
  modules = c("eiv", "original_baselines", "continuation", "complete_case_baselines",
    "prediction", "diagnostics", "surface_comparison"),
  task_ids = character(), cores = 4L, run_id = "local", continuation_iter = NULL, chains = NULL, iter = NULL, burn = NULL, ess_min = NULL) {
  ess_min <- mixedgp_parse_ess_min(ess_min)
  profile <- match.arg(profile, c("paper", "quick", "pilot"))
  action <- match.arg(action, c("plan", "fit", "run", "analyze", "continue", "report"))
  if(action=="run")action<-"fit"
  if(is.null(iter))iter<-if(profile=="paper")50000L else 80L
  if(is.null(burn))burn<-if(profile=="paper")5000L else 30L
  iter<-mixedgp_as_integer_strict(iter,"iter",1L,1L)
  burn<-mixedgp_as_integer_strict(burn,"burn",0L,1L)
  if(burn>=iter)stop("burn must be smaller than iter.")
  if(is.null(chains))chains<-if(profile=="quick")1L else 4L
  chains<-mixedgp_as_integer_strict(chains,"chains",1L,1L)
  if(!is.null(ess_min) && chains<4L)stop("Scientific ESS requires at least four chains.")
  if(profile=="quick" && chains!=1L)stop("quick uses one chain; choose pilot or paper for multiple chains.")
  if(profile=="paper" && action%in%c("plan","fit","continue") &&
     !identical(c(chains,iter,burn),c(4L,50000L,5000L)))
    stop("paper fixes chains=4, iter=50000 and burn=5000; use the profile settings.")
  if(!is.null(continuation_iter))
    continuation_iter <- ocean_continuation_policy(profile,continuation_iter)$n_iter
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
    run_id = run_id, continuation_iter = continuation_iter, ess_min=ess_min,
    iter = iter, burn = burn,
    chains = chains, u_block = 8L)
}

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

ocean_continuation_policy <- function(profile="paper",n_iter=NULL) {
  if(is.null(n_iter))n_iter<-if(profile=="pilot")80L else 20000L
  if(!is.numeric(n_iter) || length(n_iter)!=1L || !is.finite(n_iter) ||
    n_iter<1 || n_iter>.Machine$integer.max || n_iter!=floor(n_iter))
    stop("continuation_iter must be a positive integer")
  list(max_rhat=1.01,min_ess_bulk=400,n_iter=as.integer(n_iter),max_rounds=1L,
    screen="finite Rhat and bulk ESS")
}

ocean_continuation_iterations <- function(output_dir,profile="paper",fold,cell,
  requested=NULL,register=FALSE) {
  path<-ocean_path(output_dir,"continuation_plus20k",sprintf("fold%02d",fold),cell,
    "continuation_settings.rds")
  saved<-if(file.exists(path))readRDS(path)$n_iter else NULL
  if(is.null(saved) && file.exists(ocean_path(dirname(path),"fit_plus20k.rds")))
    saved<-ocean_continuation_policy(profile)$n_iter # Legacy fixed-budget outputs.
  if(!is.null(saved))saved<-ocean_continuation_policy(profile,saved)$n_iter
  if(register && !is.null(saved) && !is.null(requested) && saved!=requested)
    stop("Continuation length differs from the saved task for fold ",fold," / ",cell,
      ": saved ",saved,", requested ",requested,
      ". Reuse its saved length; a second extension is not supported.")
  n<-if(!is.null(saved))saved else ocean_continuation_policy(profile,requested)$n_iter
  if(register && !file.exists(path))mixedgp_atomic_rds(list(n_iter=n),path)
  n
}

ocean_continuation_lengths <- function(output_dir,profile=ocean_saved_profile(output_dir)) {
  paths<-list.files(if(ocean_storage_flat(output_dir))output_dir else ocean_path(output_dir,"continuation_plus20k"),
    "^continuation_settings[.]rds$",recursive=TRUE,full.names=TRUE)
  if(!length(paths))return(ocean_continuation_policy(profile)$n_iter)
  sort(unique(vapply(paths,function(path)
    ocean_continuation_policy(profile,readRDS(path)$n_iter)$n_iter,integer(1))))
}

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

ocean_continuation_decisions <- function(output_dir,profile="paper",write=FALSE,continuation_iter=NULL) {
  rows<-list()
  for(f in if(profile!="quick")1:5 else 1L) for(cell in if(profile!="quick")c("O0","O20","O50") else c("O0","O20","O50")) {
    path<-ocean_path(if(profile!="quick")"full5_class4" else "smoke_class4",sprintf("fold%02d",f),cell,"convergence_parameters.csv")
    path<-ocean_result_relative(output_dir,path)
    screen<-ocean_read_screen(ocean_path(output_dir,path))
    needs<-if(profile=="quick")FALSE else if(screen$ready)!screen$screen_pass else NA
    rows[[length(rows)+1L]]<-cbind(data.frame(fold=f,cell,initial_diagnostics=path),screen,
      needs_continuation=needs,decision=if(is.na(needs))"awaiting_initial" else if(needs)paste0("append_",ocean_continuation_iterations(output_dir,profile,f,cell,continuation_iter)) else "keep_initial")
  }
  d<-do.call(rbind,rows)
  if(write)mixedgp_atomic_csv(d,ocean_path(output_dir,"continuation_decisions.csv"))
  d
}

ocean_final_eiv_sources <- function(output_dir,profile="paper") {
  decisions<-ocean_continuation_decisions(output_dir,profile)
  rows<-lapply(seq_len(nrow(decisions)),function(i) {
    d<-decisions[i,];continued<-file.exists(ocean_path(output_dir,"continuation_plus20k",
      sprintf("fold%02d",d$fold),d$cell,"fit_plus20k.rds"))
    source<-ocean_path(if(continued)"continuation_plus20k" else if(profile!="quick")"full5_class4" else "smoke_class4",
      sprintf("fold%02d",d$fold),d$cell)
    source<-ocean_result_relative(output_dir,source)
    final<-ocean_read_screen(ocean_path(output_dir,source,"convergence_parameters.csv"))
    data.frame(fold=d$fold,cell=d$cell,source,continued,
      fit_file=if(continued)"fit_plus20k.rds" else "fit.rds",
      available=d$ready && file.exists(ocean_path(output_dir,source,"predictive_metrics.csv")),
      initial_screen_pass=d$screen_pass,final_screen_pass=final$screen_pass)
  })
  do.call(rbind,rows)
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
  ifelse(absolute,path,ocean_path(repository_dir,path))
}

ocean_test_masks <- function(test,fold) {
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(202609260L+fold)
  k20<-k50<-rep(FALSE,nrow(test))
  for(c in 1:4) {
    ids<-sample(which(test$C==c));n<-length(ids)
    q20<-round(.2*n);q50<-floor(.5*n)+as.integer(n%%2L==1L && (c+fold)%%2L==0L)
    k20[head(ids,q20)]<-TRUE;k50[head(ids,q50)]<-TRUE
  }
  data.frame(cohort_id=test$cohort_id,C=test$C,observed_U_O20=k20,
    observed_U_O50=k50,mask_seed=202609260L+fold)
}

ocean_read_inputs <- function(data_dir) {
  cohort<-read.csv(ocean_path(data_dir,"cohort_manifest.csv"))
  stopifnot(nrow(cohort)%in%c(200L,500L,1000L),!anyDuplicated(cohort$profile),!anyDuplicated(cohort$cohort_id))
  cohort$C<-findInterval(cohort$U,c(50,100,150))+1L
  cohort$oxygen_class<-factor(cohort$C,levels=1:4,labels=c("<50","50-100","100-150","150-180"))
  cohort
}

ocean_prepare_inputs <- function(data_dir) {
  cohort<-ocean_read_inputs(data_dir);N<-nrow(cohort);ntrain<-as.integer(.8*N);ntest<-N-ntrain
  stopifnot(all(table(factor(cohort$C,levels=1:4))==N/4))
  folds<-list();held<-character();audit<-list()
  for(f in 1:5) {
    tr<-read.csv(ocean_path(data_dir,sprintf("fold%02d_train.csv",f)))
    te<-read.csv(ocean_path(data_dir,sprintf("fold%02d_test.csv",f)))
    tr$C<-findInterval(tr$U,c(50,100,150))+1L;te$C<-findInterval(te$U,c(50,100,150))+1L
    tr$calib_o20<-ocean_diag_flag(tr$calib_o20);tr$calib_o50<-ocean_diag_flag(tr$calib_o50)
    stopifnot(nrow(tr)==ntrain,nrow(te)==ntest,!length(intersect(tr$profile,te$profile)),
      setequal(c(tr$cohort_id,te$cohort_id),cohort$cohort_id),!anyDuplicated(c(tr$cohort_id,te$cohort_id)),
      all(table(factor(tr$C,levels=1:4))==ntrain/4),all(table(factor(te$C,levels=1:4))==ntest/4),
      all(tr$calib_o20<=tr$calib_o50),
      all(table(factor(tr$C[tr$calib_o20],levels=1:4))==.2*ntrain/4),
      all(table(factor(tr$C[tr$calib_o50],levels=1:4))==.5*ntrain/4))
    columns<-c("cohort_id","profile","G2salinity","G2pressure","log_pressure","U","Y")
    for(d in list(tr,te)) {
      ref<-cohort[match(d$cohort_id,cohort$cohort_id),columns]
      stopifnot(isTRUE(all.equal(d[,columns],ref,check.attributes=FALSE)),
        all(is.finite(as.matrix(d[,columns[3:7]]))),all(d$U>0 & d$U<180),
        all(d$G2pressure>=300 & d$G2pressure<=1500),max(abs(d$log_pressure-log(d$G2pressure)))<1e-10)
    }
    masks<-ocean_test_masks(te,f);saved<-ocean_path(data_dir,sprintf("fold%02d_test_U_masks.csv",f))
    stopifnot(file.exists(saved),isTRUE(all.equal(read.csv(saved),masks,check.attributes=FALSE)),
      sum(masks$observed_U_O20)==.2*ntest,sum(masks$observed_U_O50)==.5*ntest)
    for(cell in c("O20","O50")) {
      obs<-tr[[paste0("calib_",tolower(cell))]];known<-masks[[paste0("observed_U_",cell)]]
      center<-mean(tr$U[obs]);scale<-sd(tr$U[obs]);stopifnot(is.finite(scale),scale>0)
      folds[[paste(f,cell,sep="/")]]<-list(train=tr,test=te,observed_train=obs,
        observed_test_mixed=known,center=center,scale=scale,masks=masks)
      audit[[length(audit)+1L]]<-data.frame(fold=f,cell,n_train=ntrain,n_test=ntest,
        n_train_exact=sum(obs),n_test_exact=sum(known),oxygen_center=center,oxygen_sd=scale)
    }
    held<-c(held,te$cohort_id)
  }
  stopifnot(length(held)==N,!anyDuplicated(held))
  list(folds=folds,cohort=cohort,response_sd=sd(cohort$Y),audit=do.call(rbind,audit),
    design=list(N=N,n_train=ntrain,n_test=ntest,n_common_hidden=ntrain/2,
      train_exact=c(O0=0L,O20=.2*ntrain,O50=.5*ntrain),test_exact=c(O0=0L,O20=.2*ntest,O50=.5*ntest)))
}

ocean_identity <- function(repository_dir,data_dir,settings) {
  packages<-c("eivGP","kergp","LVGP","EzGP")
  versions<-setNames(vapply(packages,function(p)as.character(utils::packageVersion(p)),""),packages)
  paths<-sort(list.files(data_dir,"[.]csv$",full.names=TRUE),method="radix")
  functions<-c("fit_eivgp","continue_eivgp","predict_eivgp","impute_eivgp",
    "run_published_mixedgp_competitors","mixedgp_adapter_ucgp","mixedgp_adapter_lvgp","mixedgp_adapter_ezgp")
  list(schema="ocean-package-CC-v2",profile=settings$profile,
    upstream_commit="677791de6637298ead40d87a534ccd761722d1d5",
    data=setNames(ocean_hashes(paths),basename(paths)),code=eivgp_module_identity(),
    application=unname(tools::md5sum(file.path(repository_dir,"tools/ocean-tools.R"))),packages=versions,
    package_functions=setNames(vapply(functions,function(n)digest::digest(
      paste(deparse(get(n,asNamespace("eivGP"))),collapse="\n"),algo="sha256"),""),functions),
    numerics=list(iter=settings$iter,burn=settings$burn,chains=settings$chains,u_block=8L,ess_min=settings$ess_min,
      controls=ocean_cc_controls(settings$profile),n_draw=400L,
      continuation_policy=ocean_continuation_policy(settings$profile)),
    RNG=c("Mersenne-Twister","Inversion","Rejection"))
}

ocean_assert_identity <- function(saved,current,output_dir) {
  if(identical(saved,current))return(invisible(TRUE))
  fields<-union(names(saved),names(current))
  changed<-fields[!vapply(fields,function(n)identical(saved[[n]],current[[n]]),logical(1))]
  details<-character()
  for(key in c("chains","iter","burn")) {
    old<-saved$numerics[[key]];new<-current$numerics[[key]]
    if(!identical(old,new)) {
      show<-function(x)if(is.null(x))"missing" else paste(x,collapse=",")
      details<-c(details,paste0("  ",key,": saved=",show(old),"; requested=",show(new)))
    }
  }
  # Suggest a fresh root without creating, renaming or deleting any directory.
  fresh<-paste0(output_dir,"_new");i<-1L
  while(file.exists(fresh)) {i<-i+1L;fresh<-paste0(output_dir,"_new",i)}
  stop(paste(c("Existing Ocean results are incompatible with this run.",
    paste0("Output directory: ",output_dir),
    paste0("Changed identity fields: ",paste(changed,collapse=", ")),details,
    "Existing fits and metadata have been preserved; no fitting started.",
    "To start a separate run, repeat your command with:",
    paste0("  --output=",shQuote(fresh)),
    "To resume the existing run, use its original data, package/code and initial sampling settings.",
    "Do not delete workflow_identity.rds to bypass this check."),collapse="\n"),call.=FALSE)
}

ocean_task_files <- function(task,profile="paper",stage="analyze") {
  if(stage %in% c("fit","continue") && task$module %in% c("eiv","continuation")) {
    sub<-ocean_path(if(task$module=="continuation")"continuation_plus20k" else
      if(profile=="quick")"smoke_class4" else "full5_class4",sprintf("fold%02d",task$fold),task$cell)
    return(ocean_path(sub,c(if(task$module=="continuation")"fit_plus20k.rds" else "fit.rds",
      "convergence_parameters.csv","convergence_summary.csv",
      if(task$module=="continuation")"continuation_settings.rds")))
  }
  initial <- if(profile!="quick") "full5_class4" else "smoke_class4"
  if(task$module=="original_baselines") return(ocean_path(initial,sprintf("fold%02d",task$fold),
    c("baseline_metrics.csv","baseline_package_outputs.rds","baseline_status.csv","baseline_predictions.csv")))
  if(task$module=="complete_case_baselines") return(ocean_path("complete_case_baselines/fits",task$task_id,
    c("package_output.rds","metrics.csv","predictions.csv","package_status.csv","design.rds")))
  sub <- ocean_path(if(task$module=="continuation") "continuation_plus20k" else initial,
    sprintf("fold%02d",task$fold),task$cell)
  ocean_path(sub,c(if(task$module=="continuation") "fit_plus20k.rds" else "fit.rds",
    "metrics.csv","predictive_metrics.csv","predictions.csv","mixed_predictions.csv",
    "prediction_distributions.rds","U_posterior_intervals.csv","convergence_parameters.csv",
    "convergence_summary.csv",if(task$module=="continuation")"continuation_settings.rds",
    if(task$cell=="O0") "imputations_model_scale.csv" else "imputations.csv"))
}

ocean_verify_receipt <- function(output_dir,id,identity=NULL) {
  path <- ocean_path(output_dir,"workflow_tasks",paste0(id,".done.rds"))
  if(!file.exists(path)) return(FALSE)
  receipt <- readRDS(path)
  if(!is.null(identity) && !identical(receipt$identity,identity)) stop("Task identity mismatch: ",id)
  files <- ocean_path(output_dir,receipt$files)
  if(!all(file.exists(files)) || !identical(ocean_hashes(files),receipt$hashes))
    stop("Incomplete or changed completed task: ",id)
  TRUE
}

ocean_run_task <- function(task,repository_dir,output_dir,settings,inputs,identity) {
  env <- new.env(parent=environment())
  env$root <- output_dir; env$smoke <- settings$profile=="quick"
  env$stage <- settings$action
  env$out <- ocean_path(output_dir,if(env$smoke) "smoke_class4" else "full5_class4")
  env$ess_min <- settings$ess_min
  env$iter <- settings$iter; env$burn <- settings$burn; env$chains <- settings$chains
  env$cores <- if(env$smoke) 1L else settings$cores; env$u_block <- 8L
  env$baseline_controls <- ocean_cc_controls(settings$profile)
  env$continuation_iter <- ocean_continuation_iterations(output_dir,settings$profile,
    task$fold,task$cell,settings$continuation_iter)
  env$atomic_csv <- mixedgp_atomic_csv; env$atomic_lines <- ocean_atomic_lines; env$atomic_rds <- mixedgp_atomic_rds
  env$set_task <- function(...) invisible(NULL)
  env$write_progress <- function(...) {cat(paste(...,collapse=" "),"\n");flush.console()}
  if(task$module=="complete_case_baselines") {
    ocean_run_complete_case(task,output_dir,settings,inputs)
  } else if(task$module %in% c("eiv","original_baselines")) {
    ocean_initial_functions(env)
    d <- env$read_fold(task$fold); fold_out <- ocean_path(env$out,sprintf("fold%02d",task$fold))
    dir.create(fold_out,recursive=TRUE,showWarnings=FALSE)
    if(task$module=="eiv") env$run_eiv(task$fold,task$cell,d$train,d$test,fold_out)
    else env$run_baselines(task$fold,d$train,d$test,fold_out)
  } else if(task$module=="continuation") {
    dependency <- sprintf("fit_eiv_fold%02d_%s",task$fold,task$cell)
    if(!ocean_verify_receipt(output_dir,dependency,identity)) stop("Missing initial task: ",dependency)
    env$out <- ocean_path(output_dir,"continuation_plus20k")
    env$original_path <- function(fold,cell) ocean_path(output_dir,"full5_class4",sprintf("fold%02d",fold),cell,"fit.rds")
    ocean_continuation_functions(env); env$run_one(task$fold,task$cell)
  } else stop("Unknown application task module: ",task$module)
  invisible(NULL)
}

ocean_stage_plan <- function(settings) {
  data.frame(stage=c("fit","predict","functionals","paper_assets","html"),
    depends_on=c("frozen inputs","saved fits","predictions","scores and diagnostics","paper assets"),
    unit=c("independent chains and native comparator tasks","fixed seeded posterior-draw blocks",
      "fold/calibration and posterior-draw blocks","independent tables and figure files","one document"),
    core_budget=rep(settings$cores,5L),
    concurrency=c("shared pool; resident fits bounded","shared pool; resident fits bounded",
      "shared pool; resident fits bounded","shared pool","serial assembly"),
    stringsAsFactors=FALSE)
}

ocean_job_dir <- function(output_dir,run_id) {
  if(length(run_id)!=1L || !grepl("^[A-Za-z0-9_-]+$",run_id) ||
     run_id %in% c("data","full5_class4","smoke_class4","prediction","diagnostics",
                   "continuation_plus20k","figures","receipts"))
    stop("Choose a distinct execution run_id, such as paper_fit.")
  if(ocean_storage_flat(output_dir))return(ocean_path(output_dir,"logs",run_id))
  legacy <- ocean_path(output_dir,"jobs",run_id)
  if(dir.exists(legacy))legacy else ocean_path(output_dir,run_id)
}

ocean_execute <- function(repository_dir,data_dir,output_dir,settings,
  import_dirs=character()) {
  if(settings$action=="fit")ocean_storage_init(output_dir)
  release_budget <- mixedgp_budget_enter(settings$cores)
  on.exit(release_budget(),add=TRUE)
  settings$cores <- mixedgp_budget_cores()
  RNGkind("Mersenne-Twister","Inversion","Rejection")
  import_dirs <- as.character(unlist(import_dirs)); import_dirs <- import_dirs[nzchar(import_dirs)]
  inputs <- ocean_prepare_inputs(data_dir); selected <- ocean_select_tasks(settings)
  stage<-settings$action
  if(stage=="fit")selected<-selected[selected$module!="continuation",,drop=FALSE]
  if(stage=="continue")selected<-selected[selected$module=="continuation",,drop=FALSE]
  if(stage=="analyze")selected<-selected[selected$module%in%c("eiv","continuation"),,drop=FALSE]
  execute<-stage%in%c("fit","analyze","continue")
  if(settings$action=="plan") return(list(tasks=selected,inputs=inputs,settings=settings,stages=ocean_stage_plan(settings)))
  if(execute || length(import_dirs)) {
    suppressPackageStartupMessages(library(eivGP))
    identity <- ocean_identity(repository_dir,data_dir,settings)
    dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
    lock <- ocean_path(output_dir,".workflow.lock")
    if(!dir.create(lock,showWarnings=FALSE)) stop("Output directory is locked: ",lock,
      ". Confirm the previous process has stopped before removing a stale lock.")
    on.exit(unlink(lock,recursive=TRUE),add=TRUE)
    idfile <- ocean_path(output_dir,"workflow_identity.rds")
    if(file.exists(idfile)) {
      ocean_assert_identity(readRDS(idfile),identity,output_dir)
    } else {
      if(length(setdiff(list.files(output_dir,all.files=TRUE,no..=TRUE),c(".workflow.lock",".application.lock","jobs","_meta","logs",basename(ocean_job_dir(output_dir,settings$run_id))))))
        stop("Existing outputs lack distributed provenance. Read them with action=report, or use a fresh directory.")
      mixedgp_atomic_rds(identity,idfile)
    }
    dir.create(ocean_path(output_dir,"data"),showWarnings=FALSE)
    csv <- list.files(data_dir,"[.]csv$",full.names=TRUE)
    for(f in csv) {
      to <- ocean_path(output_dir,"data",basename(f))
      if(file.exists(to)) stopifnot(identical(ocean_hashes(f),ocean_hashes(to))) else stopifnot(file.copy(f,to))
    }
    for(f in 1:5) {
      path <- ocean_path(output_dir,if(settings$profile!="quick") "full5_class4" else "smoke_class4",sprintf("fold%02d",f),"test_U_masks.csv")
      mixedgp_atomic_csv(inputs$folds[[paste(f,"O20",sep="/")]]$masks,path)
    }
    mixedgp_atomic_rds(settings,ocean_path(output_dir,"execution_settings.rds"))
    mixedgp_atomic_rds(settings,ocean_path(output_dir,paste0("execution_settings_",stage,".rds")))
    ocean_combine_results(import_dirs,output_dir,identity)
    job <- ocean_job_dir(output_dir,settings$run_id)
    status <- if(execute) selected else selected[FALSE,,drop=FALSE]
    status$state <- rep("PENDING",nrow(status)); status$detail <- rep("",nrow(status))
    status$elapsed_seconds <- rep(NA_real_,nrow(status))
    status$allocated_cores<-rep(NA_integer_,nrow(status))
    status$worker_pid<-rep(NA_integer_,nrow(status))
    publish <- function(state,current="none") {
      mixedgp_atomic_csv(status,ocean_path(job,"task_status.csv"))
      ocean_atomic_lines(c(paste("state:",state),paste("updated:",format(Sys.time(),"%Y-%m-%d %H:%M:%S %z")),
        paste("pid:",Sys.getpid()),paste("current:",current),
        sprintf("executed done: %d/%d; failed: %d; conditional skips: %d",sum(status$state %in% c("DONE","REUSED")),nrow(status)-sum(status$state=="SKIPPED_CONVERGED"),sum(status$state=="FAILED"),sum(status$state=="SKIPPED_CONVERGED")),
        paste("Total core budget:",settings$cores,"; shared across chains and analysis. See execution_schedule.csv and task logs.")),
        ocean_path(job,"PROGRESS.txt"))
    }
    mixedgp_atomic_csv(inputs$audit,ocean_path(job,"design_audit.csv"))
    capture.output(utils::sessionInfo(),file=ocean_path(job,"sessionInfo.txt"))
    pending<-integer();ids<-character(nrow(selected))
    for(i in if(execute) seq_len(nrow(selected)) else integer()) {
      task <- selected[i,,drop=FALSE]; id <- task$task_id
      if(stage%in%c("fit","continue") && task$module%in%c("eiv","continuation"))id<-paste0("fit_",id)
      ids[i]<-id
      rownames(task) <- NULL # Assignment position is not part of task identity.
      if(task$module=="continuation") {
        dependency<-sprintf("fit_eiv_fold%02d_%s",task$fold,task$cell)
        if(!ocean_verify_receipt(output_dir,dependency,identity)) {
          status$state[i]<-"FAILED";status$detail[i]<-paste("Missing verified initial fit:",dependency)
          publish("RUNNING",id);next
        }
        checkpoint<-ocean_path(output_dir,"continuation_plus20k",sprintf("fold%02d",task$fold),task$cell,"fit_plus20k.rds")
        if(stage=="analyze" && !file.exists(checkpoint)) {
          status$state[i]<-"SKIPPED_NOT_REQUESTED";publish("RUNNING",id);next
        }
        if(stage=="continue") {
          decisions<-ocean_continuation_decisions(output_dir,settings$profile,write=TRUE,
            continuation_iter=settings$continuation_iter)
          d<-decisions[decisions$fold==task$fold & decisions$cell==task$cell,]
          if(nrow(d)==1L && isTRUE(d$ready) && !isTRUE(d$needs_continuation)) {
            status$state[i]<-"SKIPPED_CONVERGED";publish("RUNNING",id);next
          }
          ocean_continuation_iterations(output_dir,settings$profile,task$fold,task$cell,
            requested=settings$continuation_iter,register=TRUE)
        }
      }
      if(ocean_verify_receipt(output_dir,id,identity)) {status$state[i]<-"REUSED";publish("RUNNING",id);next}
      cfg <- ocean_path(output_dir,"workflow_tasks",paste0(id,".config.rds"))
      if(file.exists(cfg) && !identical(readRDS(cfg),list(identity=identity,task=task))) stop("Partial-task identity mismatch: ",id)
      mixedgp_atomic_rds(list(identity=identity,task=task),cfg)
      pending<-c(pending,i)
    }
    events<-data.frame(task_id=character(),event=character(),time=numeric(),
      worker_pid=integer(),allocated_cores=integer(),active_cores=integer())
    mixedgp_atomic_csv(events,ocean_path(job,"execution_schedule.csv"))
    event<-function(i,type,active) {
      events<<-rbind(events,data.frame(task_id=ids[i],event=type,time=as.numeric(Sys.time()),
        worker_pid=status$worker_pid[i],allocated_cores=status$allocated_cores[i],active_cores=active))
      mixedgp_atomic_csv(events,ocean_path(job,"execution_schedule.csv"))
    }
    work<-function(i,allocated) {
      task<-selected[i,,drop=FALSE];rownames(task)<-NULL;id<-ids[i]
      local_settings<-settings;local_settings$cores<-allocated
      log_path<-ocean_path(job,"tasks",paste0(id,".log"))
      dir.create(dirname(log_path),recursive=TRUE,showWarnings=FALSE)
      log<-file(log_path,open="at");old_message<-getConnection(sink.number(type="message"))
      sink(log);sink(log,type="message")
      on.exit({sink(old_message,type="message");sink();close(log)},add=TRUE)
      cat("Task:",id,"PID:",Sys.getpid(),"allocated cores:",allocated,"stage:",stage,"\n")
      t0<-proc.time()[3]
      outcome<-tryCatch({
        ocean_run_task(task,repository_dir,output_dir,local_settings,inputs,identity)
        files<-ocean_result_relative(output_dir,ocean_task_files(task,settings$profile,stage))
        stopifnot(all(file.exists(ocean_path(output_dir,files))))
        mixedgp_atomic_rds(list(identity=identity,task=task,files=files,
          hashes=ocean_hashes(ocean_path(output_dir,files)),finished=as.character(Sys.time())),
          ocean_path(output_dir,"workflow_tasks",paste0(id,".done.rds")))
        NULL
      },error=function(e)e)
      detail<-if(inherits(outcome,"error"))conditionMessage(outcome) else ""
      if(nzchar(detail))cat("FAILED:",detail,"\n")
      list(state=if(inherits(outcome,"error"))"FAILED" else "DONE",detail=detail,
        elapsed_seconds=as.numeric(proc.time()[3]-t0))
    }
    start_task<-function(i,allocated,pid,active) {
      status$state[i]<<-"RUNNING";status$allocated_cores[i]<<-allocated;status$worker_pid[i]<<-pid
      event(i,"START",active);publish("RUNNING",ids[i])
    }
    finish_task<-function(i,result,active) {
      if(!is.list(result) || is.null(result$state))result<-list(state="FAILED",
        detail="Worker exited without a result; inspect its task log and memory limits.",elapsed_seconds=NA_real_)
      status$state[i]<<-result$state;status$detail[i]<<-result$detail
      status$elapsed_seconds[i]<<-result$elapsed_seconds
      event(i,"FINISH",active);publish("RUNNING",ids[i])
    }
    caps<-ifelse(selected$module[pending]%in%c("eiv","continuation"),settings$chains,1L)
    # Bound simultaneously resident fits, but allow each stage to borrow all
    # free slots. A waiting fit coordinator releases its slot to chain workers.
    outer <- mixedgp_stage_workers(settings$cores,stage,length(pending),settings$chains)
    mixedgp_budget_yield(function() mixedgp_schedule(pending,rep(1L,length(pending)),outer,
      function(i,allocated) mixedgp_budget_worker(function() work(i,settings$cores)),
      start_task,finish_task))
    publish(if(any(status$state=="FAILED")) "SELECTED_COMPLETE_WITH_FAILURES" else "SELECTED_COMPLETE")
    ocean_atomic_lines(capture.output(print(status,row.names=FALSE)),ocean_path(job,"SUMMARY.txt"))
    failed <- status$state=="FAILED"
    # Native comparator failures are reported, never allowed to discard valid EIV results.
    if(any(failed & !status$module %in% c("original_baselines","complete_case_baselines")))
      stop("Some EIV tasks failed; inspect ",ocean_path(job,"task_status.csv"))
    if(any(failed))warning("Comparator tasks failed; continuing analysis of available results. See ",
      ocean_path(job,"task_status.csv"),call.=FALSE)
  }
  decisions<-ocean_continuation_decisions(output_dir,settings$profile,write=dir.exists(output_dir),
    continuation_iter=settings$continuation_iter)
  result <- if(!stage%in%c("fit","continue") && "prediction" %in% settings$modules) ocean_prediction_report(output_dir,inputs,
    settings$profile) else NULL
  list(tasks=selected,inputs=inputs,settings=settings,prediction=result,output_dir=output_dir,continuation_decisions=decisions)
}

ocean_saved_profile <- function(output_dir) {
  p<-ocean_path(output_dir,"workflow_identity.rds")
  if(file.exists(p))readRDS(p)$profile else "paper"
}

ocean_continuation_label <- function(output_dir)
  paste(paste0("+",format(ocean_continuation_lengths(output_dir),big.mark=",",scientific=FALSE)),collapse=" / ")

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
  train <- read.csv(ocean_path(root, "data", sprintf("fold%02d_train.csv", fold)))
  test <- read.csv(ocean_path(root, "data", sprintf("fold%02d_test.csv", fold)))
  train$C <- as.integer(findInterval(train$U, c(50, 100, 150)) + 1L)
  test$C <- as.integer(findInterval(test$U, c(50, 100, 150)) + 1L)
  design<-ocean_prepare_inputs(ocean_path(root,"data"))
  stopifnot(identical(train$cohort_id,design$folds[[paste(fold,"O20",sep="/")]]$train$cohort_id))
  list(train = train, test = test)
}
run_baselines <- function(fold, train, test, fold_out) {
  metrics_path <- ocean_path(fold_out, "baseline_metrics.csv")
  fit_path <- ocean_path(fold_out,"baseline_package_outputs.rds")
  methods <- c("UC-GP","LVGP","EzGP")
  fit <- if(file.exists(fit_path))readRDS(fit_path) else NULL
  done <- if(is.null(fit))character() else fit$status$method[fit$status$status=="success"]
  pending <- setdiff(methods,done)
  if(!length(pending) && file.exists(metrics_path) &&
     setequal(read.csv(metrics_path)$method,methods))return(invisible("reused output"))
  # Retry only missing methods; preserve their original position-dependent seeds.
  for(method in pending) {
    replacement <- run_published_mixedgp_competitors(
      X_train=as.matrix(train[,c("G2salinity","log_pressure")]),y_train=train$Y,
      C_train=matrix(train$C,ncol=1L,dimnames=list(NULL,"oxygen_class")),
      X_test=as.matrix(test[,c("G2salinity","log_pressure")]),
      C_test=matrix(test$C,ncol=1L,dimnames=list(NULL,"oxygen_class")),
      m_vec=4L,n_draw=400L,seed=2026092660L+fold*100L+match(method,methods)-1L,
      methods=method,strict=FALSE,controls=baseline_controls)
    if(is.null(fit))fit<-replacement else {
      for(n in setdiff(names(replacement),"status"))fit[[n]][[method]]<-replacement[[n]][[method]]
      fit$status<-rbind(fit$status[fit$status$method!=method,,drop=FALSE],replacement$status)
    }
    # Atomic replacement also avoids writing through a read-only source symlink.
    atomic_rds(fit,fit_path)
  }
  fit$status<-fit$status[match(methods,fit$status$method),,drop=FALSE]
  atomic_rds(fit,fit_path)
  atomic_csv(fit$status, ocean_path(fold_out, "baseline_status.csv"))
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
    ocean_path(fold_out, "baseline_predictions.csv"))
  failed <- fit$status$method[fit$status$status != "success"]
  if (length(failed)) stop("Baseline failures: ", paste(failed, collapse = ","))
  invisible("all three baselines fitted")
}
run_eiv <- function(fold, cell, train, test, fold_out) {
  cell_out <- ocean_path(fold_out, cell)
  dir.create(cell_out, showWarnings = FALSE)
  metric_path <- ocean_path(cell_out, "metrics.csv")
  if (file.exists(metric_path) && !identical(get0("stage",inherits=TRUE),"fit")) return(invisible("reused output"))
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
  fit_path <- ocean_path(cell_out, "fit.rds")
  if(identical(get0("stage",inherits=TRUE),"analyze") && !file.exists(fit_path))
    stop("Missing fit; run ocean.R --action=fit first: ",fit_path)
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
    capture.output(summary(model), file = ocean_path(cell_out, "fit_summary.txt"))
    model
  }
  if(identical(get0("stage",inherits=TRUE),"fit") && length(fit$checkpoint$states)>=4L) {
    rows<-mixedgp_study2_diagnostic_rows(x_test,c_test,5L)
    panel<-list(X=x_test[rows,,drop=FALSE],C=as.integer(c_test[rows,1L]),
      U=if(any(mask))test$U[rows]else NULL,n_latent=64L,seed=202610240L+fold*100L+seed_offset)
    fit<-mixedgp_application_ess(fit,get0("ess_min",inherits=TRUE),panel,
      ocean_path(cell_out,"scientific-ess-state.rds"),cores)
    atomic_rds(fit,fit_path)
    mixedgp_save_scientific_trace(fit,ocean_path(cell_out,"scientific-trace-targets"),
      sprintf("Ocean | fold %d | %s",fold,cell))
    unlink(ocean_path(cell_out,"scientific-ess-state.rds"))
  }
  diagnostic_parts <- lapply(
    fit$diagnostics[c("rhat_hyper", "rhat_tau", "rhat_u")],
    function(z) z[, c("parameter", "rhat", "ess_bulk", "ess_tail", "target_pass"),
                  drop = FALSE])
  diagnostics <- do.call(rbind, diagnostic_parts)
  atomic_csv(diagnostics, ocean_path(cell_out, "convergence_parameters.csv"))
  finite_rhat <- diagnostics$rhat[is.finite(diagnostics$rhat)]
  finite_ess <- diagnostics$ess_bulk[is.finite(diagnostics$ess_bulk)]
  atomic_csv(data.frame(iterations=fit$checkpoint$iteration,
                        max_finite_rhat = if (length(finite_rhat)) max(finite_rhat) else NA_real_,
                        min_finite_ess_bulk = if (length(finite_ess)) min(finite_ess) else NA_real_,
                        limited_screen_pass = length(finite_rhat) > 0L && length(finite_ess) > 0L &&
                          max(finite_rhat) <= 1.01 && min(finite_ess) >= 400),
             ocean_path(cell_out, "convergence_summary.csv"))
  if(identical(get0("stage",inherits=TRUE),"fit"))return(invisible("fit saved; analysis deferred"))
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
    predicted_Y = point), ocean_path(cell_out, "predictions.csv"))
  metrics <- score(point, test$Y, fold, "EIV-GP", cell)

  # Freeze nested test masks once per fold using the prespecified seed, without Y.
  masks<-ocean_test_masks(test,fold);known20<-masks$observed_U_O20;known50<-masks$observed_U_O50
  known<-switch(cell,O0=rep(FALSE,nrow(test)),O20=known20,O50=known50)
  atomic_csv(data.frame(cohort_id=test$cohort_id,C=test$C,observed_U_O20=known20,observed_U_O50=known50,mask_seed=202609260L+fold),
    ocean_path(fold_out,"test_U_masks.csv"))
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
    ocean_path(cell_out,"prediction_distributions.rds"))
  atomic_csv(data.frame(cohort_id=test$cohort_id,observed_Y=test$Y,predicted_Y=colMeans(mixed),test_U_observed=known),
    ocean_path(cell_out,"mixed_predictions.csv"))
  atomic_csv(rbind(probability_score(response_draws,point,test$Y,fold,"EIV-GP xc",cell),
    probability_score(mixed,colMeans(mixed),test$Y,fold,"EIV-GP mixed",cell)),ocean_path(cell_out,"predictive_metrics.csv"))

  # 4. Impute hidden training U; O0 uses model coordinates, O20/O50 raw units.
  hidden <- which(!is.finite(u_obs[, 1L]))
  if (cell == "O0") {
    ## Without calibration, raw U is not affinely identified by the package.
    imp <- eivGP::impute_eivgp(fit, rows = hidden, draw_ids = draw_ids, scale = "model")
    imp_point <- colMeans(matrix(imp[, , 1L], nrow = length(draw_ids)))
    atomic_csv(data.frame(cohort_id = train$cohort_id[hidden], C = train$C[hidden],
      true_U_raw = train$U[hidden], EIV_U_model_mean = imp_point),
      ocean_path(cell_out, "imputations_model_scale.csv"))
  } else {
    # RMD-BEGIN eiv-impute
    imp <- eivGP::impute_eivgp(fit, rows = hidden, draw_ids = draw_ids, scale = "raw")
    # RMD-END eiv-impute
    imp_point <- colMeans(matrix(imp[, , 1L], nrow = length(draw_ids)))
    calib_means <- tapply(train$U[mask], train$C[mask], mean)
    imp_base <- as.numeric(calib_means[as.character(train$C[hidden])])
    atomic_csv(data.frame(cohort_id = train$cohort_id[hidden], C = train$C[hidden],
      true_U = train$U[hidden], EIV_U_mean = imp_point,
      class_calibration_mean = imp_base), ocean_path(cell_out, "imputations.csv"))
    metrics <- rbind(metrics,
      score(imp_point, train$U[hidden], fold, "EIV U imputation", cell),
      score(imp_base, train$U[hidden], fold, "calibration class mean", cell))
  }
  iq<-apply(matrix(imp[,,1L],nrow=length(draw_ids)),2L,quantile,c(.025,.975))
  atomic_csv(data.frame(cohort_id=train$cohort_id[hidden],C=train$C[hidden],posterior_mean=imp_point,q025=iq[1,],q975=iq[2,],scale=if(cell=="O0")"model" else "raw"),ocean_path(cell_out,"U_posterior_intervals.csv"))
  atomic_csv(metrics, metric_path)
  tau <- as.matrix(fit$mcmc$samples_tau)
  atomic_csv(data.frame(boundary = seq_len(ncol(tau)),
    q025_model_scale = apply(tau, 2L, quantile, .025),
    median_model_scale = apply(tau, 2L, median),
    q975_model_scale = apply(tau, 2L, quantile, .975)),
    ocean_path(cell_out, "threshold_posterior_model_scale.csv"))
  invisible(paste0("Y RMSE ", round(metrics$RMSE[1L], 3)))
}


  }), envir=envir)
  invisible(envir)
}

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
  task_out <- ocean_path(out, sprintf("fold%02d", fold), cell)
  dir.create(task_out, recursive = TRUE, showWarnings = FALSE)
  metrics_path <- ocean_path(task_out, "metrics.csv")
  fit_path <- ocean_path(task_out, "fit_plus20k.rds")
  if (file.exists(metrics_path)) {
    fit <- readRDS(fit_path)
    set_task(fold, cell, "DONE", "reused complete continuation",
             diag = diagnostic_row(fit))
    return(invisible(NULL))
  }
  set_task(fold, cell, "RUNNING", "continuing saved chains")
  write_progress("RUNNING", sprintf("fold%02d %s", fold, cell), "continuing original chains")
  t0 <- Sys.time()
  if(identical(get0("stage",inherits=TRUE),"analyze") && !file.exists(fit_path))
    stop("Missing continuation; run ocean.R --action=continue first.")
  fit <- if (file.exists(fit_path)) readRDS(fit_path) else {
    original <- readRDS(original_path(fold, cell))
    # Continue the saved chains/RNG states, then repeat prediction and imputation.
    # RMD-BEGIN eiv-continue
    continued <- eivGP::continue_eivgp(original, n_iter = continuation_iter,
      parallel = cores > 1L, n_cores = cores, verbose = TRUE)
    # RMD-END eiv-continue
    atomic_rds(continued, fit_path)
    rm(original); gc(FALSE)
    continued
  }
  diag <- diagnostic_row(fit)
  diagnostic_parts<-lapply(fit$diagnostics[c("rhat_hyper","rhat_tau","rhat_u")],function(z)z[,c("parameter","rhat","ess_bulk","ess_tail","target_pass"),drop=FALSE])
  atomic_csv(do.call(rbind,diagnostic_parts),ocean_path(task_out,"convergence_parameters.csv"))
  atomic_csv(data.frame(max_finite_rhat=diag$max_rhat,min_finite_ess_bulk=diag$min_ess_bulk,limited_screen_pass=diag$pass_screen),ocean_path(task_out,"convergence_summary.csv"))
  atomic_csv(data.frame(fold = fold, cell = cell, stage = paste0("plus",continuation_iter),
    max_rhat = diag$max_rhat, min_ess_bulk = diag$min_ess_bulk,
    pass_screen = diag$pass_screen), ocean_path(task_out, "diagnostics.csv"))
  capture.output(summary(fit), file = ocean_path(task_out, "fit_summary.txt"))
  if(identical(get0("stage",inherits=TRUE),"continue"))return(invisible("continued fit saved; analysis deferred"))

  train <- read.csv(ocean_path(root, "data", sprintf("fold%02d_train.csv", fold)))
  test <- read.csv(ocean_path(root, "data", sprintf("fold%02d_test.csv", fold)))
  train$C <- as.integer(findInterval(train$U, c(50, 100, 150)) + 1L)
  test$C <- as.integer(findInterval(test$U, c(50, 100, 150)) + 1L)
  stopifnot((nrow(train)+nrow(test))%in%c(200L,500L,1000L),
    !length(intersect(train$profile,test$profile)),all(table(test$C)==nrow(test)/4))
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
    predicted_Y = point), ocean_path(task_out, "predictions.csv"))

  mask <- if (cell == "O0") rep(FALSE, nrow(train)) else
    train[[paste0("calib_", tolower(cell))]] == 1L

  # Use the fixed test observation masks for both initial and continued fits.
  saved_masks<-read.csv(ocean_path(root,"full5_class4",sprintf("fold%02d",fold),"test_U_masks.csv"))
  stopifnot(identical(as.character(saved_masks$cohort_id),as.character(test$cohort_id)))
  known<-if(cell=="O0")rep(FALSE,nrow(test)) else as.logical(saved_masks[[paste0("observed_U_",cell)]])
  stopifnot(sum(known)==switch(cell,O0=0,O20=.2,O50=.5)*nrow(test))
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
  atomic_rds(list(xc=response_draws,mixed=mixed,known=known,draw_ids=draw_ids),ocean_path(task_out,"prediction_distributions.rds"))
  atomic_csv(data.frame(cohort_id=test$cohort_id,observed_Y=test$Y,predicted_Y=colMeans(mixed),test_U_observed=known),ocean_path(task_out,"mixed_predictions.csv"))
  atomic_csv(rbind(probability_score(response_draws,point,test$Y,fold,cell,"EIV-GP xc"),
    probability_score(mixed,colMeans(mixed),test$Y,fold,cell,"EIV-GP mixed")),ocean_path(task_out,"predictive_metrics.csv"))
  hidden <- which(!mask)
  imp <- eivGP::impute_eivgp(fit, rows = hidden, draw_ids = draw_ids,
                      scale = if (cell == "O0") "model" else "raw")
  imp_point <- colMeans(matrix(imp[, , 1L], nrow = length(draw_ids)))
  if (cell == "O0") {
    atomic_csv(data.frame(cohort_id = train$cohort_id[hidden], C = train$C[hidden],
      true_U_raw = train$U[hidden], EIV_U_model_mean = imp_point),
      ocean_path(task_out, "imputations_model_scale.csv"))
  } else {
    calib_means <- tapply(train$U[mask], train$C[mask], mean)
    imp_base <- as.numeric(calib_means[as.character(train$C[hidden])])
    atomic_csv(data.frame(cohort_id = train$cohort_id[hidden], C = train$C[hidden],
      true_U = train$U[hidden], EIV_U_mean = imp_point,
      class_calibration_mean = imp_base),
      ocean_path(task_out, "imputations.csv"))
    metrics <- rbind(metrics,
      score(imp_point, train$U[hidden], fold, cell, "EIV U imputation"),
      score(imp_base, train$U[hidden], fold, cell, "calibration class mean"))
  }
  iq<-apply(matrix(imp[,,1L],nrow=length(draw_ids)),2L,quantile,c(.025,.975))
  atomic_csv(data.frame(cohort_id=train$cohort_id[hidden],C=train$C[hidden],posterior_mean=imp_point,q025=iq[1,],q975=iq[2,],
    scale=if(cell=="O0")"model" else "raw"),ocean_path(task_out,"U_posterior_intervals.csv"))
  atomic_csv(metrics, metrics_path)
  minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  set_task(fold, cell, "DONE", paste0("Y RMSE ", round(metrics$RMSE[1], 3)),
           minutes = minutes, diag = diag)
  write_progress("RUNNING", sprintf("fold%02d %s", fold, cell), "continuation and scoring complete")
}


  }), envir=envir)
  invisible(envir)
}

ocean_cc_controls <- function(profile="paper") {
  if(profile=="paper") return(mixedgp_competitor_controls())
  list(`UC-GP`=list(n_starts=1L),
    LVGP=list(n_starts=1L,max_retries=1L,max_iter_ini=20L,max_iter_lat=5L,
      rescue_iter_ini=30L,rescue_iter_lat=10L,max_elapsed_seconds=120,parallel=FALSE),
    EzGP=list(cv_folds=2L,maxeval=15L,cv_score="nlpd"))
}

ocean_run_complete_case <- function(task,output_dir,settings,inputs) {
  d <- inputs$folds[[paste(task$fold,task$cell,sep="/")]]
  train <- d$train[d$observed_train,,drop=FALSE]
  test <- d$test[d$observed_test_mixed,,drop=FALSE]
  stopifnot(nrow(train)==sum(d$observed_train),
    nrow(test)==sum(d$observed_test_mixed),!anyNA(train$U),!anyNA(test$U),
    !length(intersect(train$cohort_id,test$cohort_id)))
  # Standardize within the complete training subset; neither test U nor hidden
  # training U defines transformations, optimization or model selection.
  cols <- c("G2salinity","log_pressure","U")
  center <- colMeans(train[,cols]); spread <- vapply(train[,cols],sd,0)
  stopifnot(all(is.finite(spread)),all(spread>0))
  transform <- function(z) sweep(sweep(as.matrix(z[,cols]),2L,center,"-"),2L,spread,"/")
  dir <- ocean_path(output_dir,"complete_case_baselines/fits",task$task_id)
  dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  design <- list(protocol="native-CC-observed-u-subset-v1",profile=settings$profile,
    fold=task$fold,cell=task$cell,method=task$method,seed=task$fit_seed,
    training_ids=train$cohort_id,test_ids=test$cohort_id,center=center,scale=spread,
    controls=ocean_cc_controls(settings$profile),n_draw=400L,
    omitted_test_ids=d$test$cohort_id[!d$observed_test_mixed])
  mixedgp_atomic_rds(design,ocean_path(dir,"design.rds"))
  ans <- eivGP:::run_published_mixedgp_competitors(
    X_train=transform(train),y_train=train$Y,C_train=matrix(train$C,ncol=1L),
    X_test=transform(test),C_test=matrix(test$C,ncol=1L),m_vec=4L,n_draw=400L,
    seed=task$fit_seed,methods=task$method,strict=FALSE,controls=design$controls)
  mixedgp_atomic_rds(ans,ocean_path(dir,"package_output.rds"))
  mixedgp_atomic_csv(ans$status,ocean_path(dir,"package_status.csv"))
  if(!task$method %in% names(ans$means) || any(ans$status$status!="success"))
    stop("CC native fitting failed: ",paste(ans$status$message,collapse="; "))
  # The original probability scoring function accepts any matched test size.
  env <- new.env(parent=environment());ocean_initial_functions(env)
  point <- as.numeric(ans$means[[task$method]])
  score <- env$probability_score(ans$draws[[task$method]],point,test$Y,
    task$fold,paste0(task$method,"(CC)"),task$cell)
  score$profile <- settings$profile;score$n_train <- nrow(train)
  mixedgp_atomic_csv(score,ocean_path(dir,"metrics.csv"))
  mixedgp_atomic_csv(data.frame(cohort_id=test$cohort_id,observed_Y=test$Y,
    predicted_Y=point),ocean_path(dir,"predictions.csv"))
  invisible(design)
}

run_ocean_analysis <- function(data_dir,output_dir,action="plan",profile="paper",cores=4L,
  repository_dir=NULL,settings=NULL,import_dirs=character(),cc_output_dir="",diagnostics_dir="",
  representative_fold=1L,...) {
  if(is.null(repository_dir))repository_dir<-getwd()
  if(is.null(settings))settings<-ocean_settings(profile=profile,action=action,cores=cores,...)
  ocean_analysis_impl(repository_dir,data_dir,output_dir,settings,import_dirs,
    cc_output_dir,diagnostics_dir,representative_fold)
}

ocean_analysis_impl <- function(repository_dir,data_dir,output_dir,settings,import_dirs=character(),
  cc_output_dir="",diagnostics_dir="",representative_fold=1L) {
  release_budget <- mixedgp_budget_enter(settings$cores)
  on.exit(release_budget(),add=TRUE)
  settings$cores <- mixedgp_budget_cores()
  a<-ocean_execute(repository_dir,data_dir,output_dir,settings,import_dirs)
  if(settings$action%in%c("plan","fit","continue"))return(a)
  if("prediction" %in% settings$modules && !is.null(a$prediction)) {
    a$cc_review<-ocean_cc_prediction_tables(output_dir,a$inputs,a$prediction,
      cc_output_dir=if(nzchar(cc_output_dir))cc_output_dir else output_dir,
      profile=settings$profile,expected_folds=settings$folds)
    a$cc_review$SD_Y<-a$inputs$response_sd
  }
  if(settings$profile!="quick" && any(c("diagnostics","surface_comparison") %in% settings$modules)) {
    available<-a$prediction
    if(is.null(available))available<-ocean_prediction_report(output_dir,a$inputs,settings$profile)
    ready<-all(available$availability$available) && all(vapply(1:5,function(f)
      file.exists(ocean_path(output_dir,"full5_class4",sprintf("fold%02d",f),"baseline_package_outputs.rds")),logical(1)))
    if(ready) {
      d<-ocean_run_diagnostics(output_dir,data_dir,
        report_dir=if(nzchar(diagnostics_dir))diagnostics_dir else ocean_path(output_dir,"diagnostics"),
        representative_fold=representative_fold,
        build_prediction=FALSE,prediction_report=available,
        build_surfaces="surface_comparison" %in% settings$modules)
      a$diagnostic_run<-d;a$diagnostics<-d$diagnostics;a$surfaces<-d$surfaces
    }
  }
  a$paper_assets <- ocean_export_application(output_dir,a,settings,
    data_dir=data_dir,diagnostics_dir=if(nzchar(diagnostics_dir))diagnostics_dir else ocean_path(output_dir,"diagnostics"))
  a
}

ocean_storage_flat <- function(root) {
  marker <- file.path(root,"_meta","storage-layout.rds")
  if(!file.exists(marker))return(FALSE)
  schema <- readRDS(marker)
  if(!identical(schema$schema,"ocean-flat-v1"))stop("Unknown ocean storage schema: ",marker)
  TRUE
}

ocean_flat_relative <- function(path) {
  path <- gsub("\\\\","/",path)
  path <- sub("^_meta/initial(?=/|$)","full5_class4",path,perl=TRUE)
  path <- sub("^_meta/continuation(?=/|$)","continuation_plus20k",path,perl=TRUE)
  path <- sub("^_meta/cc(?=/|$)","complete_case_baselines/fits",path,perl=TRUE)
  path <- sub("^baseline_(fold[0-9]+)/(?=O[0-9]+(?:/|$))","full5_class4/\\1/",path,perl=TRUE)
  path <- sub("^(full5_class4|smoke_class4)/(fold[0-9]+)/(O[0-9]+)(?=/|$)","eiv_\\2_\\3",path,perl=TRUE)
  path <- sub("^(full5_class4|smoke_class4)/(fold[0-9]+)(?=/|$)","baseline_\\2",path,perl=TRUE)
  path <- sub("^(full5_class4|smoke_class4)(?=/|$)","_meta/initial",path,perl=TRUE)
  path <- sub("^continuation_plus20k/(fold[0-9]+)/(O[0-9]+)(?=/|$)","eiv_\\1_\\2_continued",path,perl=TRUE)
  path <- sub("^continuation_plus20k(?=/|$)","_meta/continuation",path,perl=TRUE)
  path <- sub("^complete_case_baselines/fits/(cc_[^/]+)(?=/|$)","\\1",path,perl=TRUE)
  path <- sub("^complete_case_baselines/fits(?=/|$)","_meta/cc",path,perl=TRUE)
  path <- sub("^workflow_tasks(?=/|$)","_meta/tasks",path,perl=TRUE)
  path <- sub("^jobs/([^/]+)/tasks/","logs/\\1/",path)
  path <- sub("^jobs(?=/|$)","logs",path,perl=TRUE)
  path <- sub("^logs/([^/]+)/tasks/","logs/\\1/",path)
  path <- sub("^(prediction_report|diagnostics|surface_report)(?=/|$)","reports/\\1",path,perl=TRUE)
  meta <- grepl("^(workflow_identity[.]rds|execution_settings[^/]*[.]rds|continuation_decisions[.]csv)$",path)
  path[meta] <- paste0("_meta/",path[meta])
  path
}

ocean_result_relative <- function(root,path) {
  if(ocean_storage_flat(root))ocean_flat_relative(path)else path
}

ocean_path <- function(...,fsep=.Platform$file.sep) {
  paths <- file.path(...,fsep=fsep)
  vapply(paths,function(path) {
    candidate <- path
    repeat {
      if(dir.exists(candidate) && file.exists(file.path(candidate,"_meta","storage-layout.rds"))) {
        if(!ocean_storage_flat(candidate))return(path)
        if(identical(path,candidate))return(path)
        rel <- substring(path,nchar(candidate)+2L)
        return(file.path(candidate,ocean_flat_relative(rel)))
      }
      parent <- dirname(candidate)
      if(identical(parent,candidate)||candidate==".")return(path)
      candidate <- parent
    }
  },character(1),USE.NAMES=FALSE)
}

ocean_storage_init <- function(root) {
  if(ocean_storage_flat(root))return(invisible(root))
  if(dir.exists(file.path(root,".workflow.lock")))stop("Output directory is locked: ",root)
  existing <- if(dir.exists(root))list.files(root,all.files=TRUE,no..=TRUE)else character()
  # Never opt an existing, unmarked result tree into a new schema implicitly.
  if(length(existing))return(invisible(root))
  dir.create(file.path(root,"_meta"),recursive=TRUE,showWarnings=FALSE)
  saveRDS(list(schema="ocean-flat-v1"),file.path(root,"_meta","storage-layout.rds"))
  invisible(root)
}

ocean_migrate_storage <- function(source,destination,dry_run=TRUE) {
  source <- normalizePath(source,mustWork=TRUE)
  destination <- path.expand(destination)
  if(!dir.exists(dirname(destination)))stop("Destination parent must exist.")
  destination <- file.path(normalizePath(dirname(destination)),basename(destination))
  if(file.exists(destination)||dir.exists(destination))stop("Destination already exists: ",destination)
  if(startsWith(paste0(destination,"/"),paste0(source,"/")))stop("Destination must be outside source.")
  lock <- file.path(source,".workflow.lock")
  if(dir.exists(lock) || dir.exists(file.path(source,".application.lock")))
    stop("Source is locked. Wait for its job to finish; do not migrate an active run.")
  if(ocean_storage_flat(source))stop("Source already uses the flat layout.")
  files <- character(); pending <- ""
  while(length(pending)) {
    directory <- pending[1L];pending <- pending[-1L]
    entries <- list.files(file.path(source,directory),all.files=TRUE,no..=TRUE)
    relative <- if(nzchar(directory))file.path(directory,entries)else entries
    full <- file.path(source,relative)
    if(any(nzchar(Sys.readlink(full))))stop("Migration does not follow symbolic links.")
    dirs <- dir.exists(full)
    files <- c(files,relative[!dirs]);pending <- c(pending,relative[dirs])
  }
  files <- sort(files)
  mapped <- ocean_flat_relative(files)
  if(anyDuplicated(mapped))stop("Destination collision in migration plan.")
  plan <- data.frame(source=files,destination=mapped,bytes=file.info(file.path(source,files))$size)
  if(isTRUE(dry_run))return(plan)
  if(!dir.create(lock,showWarnings=FALSE))stop("Source became locked.")
  on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  current_files <- sort(list.files(source,recursive=TRUE,all.files=TRUE,no..=TRUE))
  current_files <- current_files[!dir.exists(file.path(source,current_files))]
  if(!identical(files,current_files))stop("Source inventory changed before migration acquired its lock.")
  stage <- tempfile(".ocean-migration-",tmpdir=dirname(destination))
  dir.create(stage);on.exit(if(dir.exists(stage))unlink(stage,recursive=TRUE),add=TRUE)
  before <- unname(tools::md5sum(file.path(source,files)))
  if(anyNA(before))stop("Cannot hash all source files.")
  for(i in seq_along(files)) {
    target <- file.path(stage,mapped[i]);dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
    if(!file.copy(file.path(source,files[i]),target,copy.mode=TRUE,copy.date=TRUE))stop("Copy failed: ",files[i])
  }
  if(!identical(before,unname(tools::md5sum(file.path(stage,mapped)))) ||
     !identical(before,unname(tools::md5sum(file.path(source,files)))))stop("Source changed or copy verification failed.")
  # Receipts record actual relative locations; scientific identities/hashes stay intact.
  receipts <- mapped[grepl("^_meta/tasks/.*[.]done[.]rds$",mapped)]
  for(p in receipts) {
    receipt <- readRDS(file.path(stage,p))
    if(any(grepl("(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$))",receipt$files)))stop("Unsafe receipt path: ",p)
    receipt$files <- ocean_flat_relative(receipt$files)
    if(!identical(unname(tools::md5sum(file.path(stage,receipt$files))),receipt$hashes))stop("Invalid source receipt: ",p)
    saveRDS(receipt,file.path(stage,p))
  }
  dir.create(file.path(stage,"_meta"),showWarnings=FALSE)
  saveRDS(list(schema="ocean-flat-v1"),file.path(stage,"_meta","storage-layout.rds"))
  plan$source_md5 <- before
  plan$destination_md5 <- unname(tools::md5sum(file.path(stage,mapped)))
  write.csv(plan,file.path(stage,"_meta","migration-manifest.csv"),row.names=FALSE)
  saveRDS(list(source=source,source_identity=readRDS(file.path(source,"workflow_identity.rds")),
    migrated_at=as.character(Sys.time()),note="Original scientific provenance retained; code compatibility is still enforced."),
    file.path(stage,"_meta","migration.rds"))
  if(file.exists(destination)||!file.rename(stage,destination))stop("Cannot publish migration destination.")
  invisible(plan)
}

ocean_legacy_task_relative <- function(path,profile="paper") {
  initial <- if(isTRUE(profile=="quick"))"smoke_class4"else"full5_class4"
  path <- sub("^eiv_(fold[0-9]+)_(O[0-9]+)_continued(?=/|$)","continuation_plus20k/\\1/\\2",path,perl=TRUE)
  path <- sub("^eiv_(fold[0-9]+)_(O[0-9]+)(?=/|$)",paste0(initial,"/\\1/\\2"),path,perl=TRUE)
  path <- sub("^baseline_(fold[0-9]+)(?=/|$)",paste0(initial,"/\\1"),path,perl=TRUE)
  sub("^(cc_[^/]+)(?=/|$)","complete_case_baselines/fits/\\1",path,perl=TRUE)
}

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

ocean_combine_results <- function(source_dirs, output_dir, identity) {
  sources <- as.character(unlist(source_dirs)); sources <- sources[nzchar(sources)]
  plan <- list(); receipts <- list()
  for(source in sources) {
    source <- normalizePath(source,mustWork=TRUE)
    if(identical(source,normalizePath(output_dir,mustWork=TRUE))) stop("Cannot import a directory into itself.")
    path <- ocean_path(source,"workflow_identity.rds")
    if(!file.exists(path) || !identical(readRDS(path),identity)) stop("Incompatible imported run: ",source)
    done <- list.files(ocean_path(source,"workflow_tasks"),"[.]done[.]rds$",full.names=TRUE)
    for(receipt_path in done) {
      receipt <- readRDS(receipt_path); id <- sub("[.]done[.]rds$","",basename(receipt_path))
      ocean_verify_receipt(source,id,identity)
      destination_files <- ocean_result_relative(output_dir,
        ocean_legacy_task_relative(receipt$files,identity$profile))
      for(i in seq_along(receipt$files)) {
        rel <- receipt$files[i]
        from <- ocean_path(source,rel); to <- ocean_path(output_dir,destination_files[i])
        key <- destination_files[i]
        hash <- ocean_hashes(from)
        if(file.exists(to) && !identical(ocean_hashes(to),hash)) stop("Conflicting existing payload: ",rel)
        if(!is.null(plan[[key]]) && !identical(plan[[key]]$hash,hash)) stop("Conflicting imported payloads: ",rel)
        plan[[key]] <- list(from=from,to=to,hash=hash)
      }
      if(!is.null(receipts[[id]]) && !identical(receipts[[id]]$hashes,receipt$hashes))
        stop("Conflicting completed task: ",id)
      receipt$files <- destination_files
      receipts[[id]] <- receipt
    }
  }
  for(item in plan) {
    dir.create(dirname(item$to),recursive=TRUE,showWarnings=FALSE)
    if(!file.exists(item$to)) stopifnot(file.copy(item$from,item$to))
    stopifnot(identical(ocean_hashes(item$to),item$hash))
  }
  for(id in names(receipts)) mixedgp_atomic_rds(receipts[[id]],
    ocean_path(output_dir,"workflow_tasks",paste0(id,".done.rds")))
  invisible(names(receipts))
}

ocean_prediction_report <- function(output_dir,inputs,profile="paper",report_dir=ocean_path(output_dir,"prediction_report")) {
  rows <- list(); inventory <- list(); original <- if(profile!="quick") "full5_class4" else "smoke_class4"
  final_sources <- ocean_final_eiv_sources(output_dir,profile)
  add <- function(path,f,cell,variant) {
    if(!file.exists(path)) return(invisible(NULL))
    z <- read.csv(path,stringsAsFactors=FALSE)
    stopifnot(nrow(z)>0L,all(z$fold==f),all(z$n_test==nrow(inputs$folds[[paste(f,"O20",sep="/")]]$test)))
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
    initial <- ocean_path(output_dir,original,sprintf("fold%02d",f))
    add(ocean_path(initial,"baseline_metrics.csv"),f,"baseline","Original baseline")
    for(cell in if(profile!="quick") c("O0","O20","O50") else c("O0","O20","O50")) {
      final<-final_sources[final_sources$fold==f & final_sources$cell==cell,]
      path<-ocean_path(output_dir,final$source)
      available<-isTRUE(final$available)
      inventory[[length(inventory)+1L]]<-data.frame(fold=f,cell,
        stage=if(final$continued)ocean_continuation_label(output_dir) else "initial",available,
        initial_screen_pass=final$initial_screen_pass,final_screen_pass=final$final_screen_pass)
      # Use the latest explicitly saved stage; diagnostic flags remain visible.
      if(available)add(ocean_path(path,"predictive_metrics.csv"),f,cell,"EIV")
    }
  }
  raw <- mixedgp_bind(rows); summaries <- list()
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
  summary_raw <- mixedgp_bind(summaries); scaled <- raw
  if(nrow(scaled)) {
    for(m in c("RMSE","MAE","CRPS","Width95","IntervalScore95")) scaled[[m]] <- scaled[[m]]/inputs$response_sd
    scaled$NLPD <- scaled$NLPD-log(inputs$response_sd)
  }
  plots <- list()
  report <- list(fold_metrics=raw,summary=summary_raw,availability=mixedgp_bind(inventory),
    SD_Y=inputs$response_sd,plots=plots,fold_metrics_SD_Y=scaled)
  if(!is.null(report_dir)) ocean_export_prediction_report(report,report_dir)
  report

}

ocean_mean_sd <- function(mean,sd,digits=2L) {
  ifelse(is.na(sd),sprintf(paste0("%.",digits,"f (SD unavailable)"),mean),
    sprintf(paste0("%.",digits,"f (%.",digits,"f)"),mean,sd))
}

ocean_tex_note <- function(note) gsub("%","\\%",note,fixed=TRUE)

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

ocean_table_columns <- function(d,columns) {
  mask<-attr(d,"best")
  out<-d[,columns,drop=FALSE]
  if(!is.null(mask)) attr(out,"best")<-mask[,intersect(columns,colnames(mask)),drop=FALSE]
  out
}

ocean_prediction_panel_html <- function(d, caption, note="") {
  if(!nrow(d)) return("")
  columns <- setdiff(names(d),c("Training exact U","Test oxygen"))
  styled <- mixedgp_table_markup(d,"html")
  training_groups <- unique(d[["Training exact U"]])
  out <- c("<table class='ocean-prediction table' style='width:100%;border-collapse:collapse'>",
    paste0("<caption style='text-align:left'>",mixedgp_html_escape(caption),"</caption>"),
    "<thead><tr>",paste0("<th style='text-align:",ifelse(columns=="Method","left","right"),
      ";padding:7px 10px'>",mixedgp_html_escape(columns),"</th>"),"</tr></thead><tbody>")
  for(j in seq_along(training_groups)) {
    training_rows <- which(d[["Training exact U"]]==training_groups[j])
    tests <- unique(d[["Test oxygen"]][training_rows])
    panel <- match(training_groups[j],c("0%","20%","50%"))
    if(is.na(panel))panel<-j
    title <- paste0(LETTERS[panel],". Training calibration: ",training_groups[j])

    out <- c(out,paste0("<tr class='comparison-panel'><th scope='rowgroup' colspan='",length(columns),
      "' style='text-align:left;background:#f2f4f6;border-top:2px solid #bbb;padding:9px 10px'>",
      mixedgp_html_escape(title),"</th></tr>"))
    for(test in tests) {
      rows <- training_rows[d[["Test oxygen"]][training_rows]==test]
      out <- c(out,paste0(
        "<tr class='comparison-condition'><th scope='rowgroup' colspan='",length(columns),
        "' style='text-align:left;font-weight:500;padding:7px 10px;border-bottom:1px solid #ddd'>",
        "Test inputs: ",mixedgp_html_escape(test),"</th></tr>"))
      for(i in rows) out <- c(out,"<tr>",paste0("<td style='text-align:",
        ifelse(columns=="Method","left","right"),";padding:6px 10px'>",
        unlist(styled[i,columns,drop=FALSE],use.names=FALSE),"</td>"),"</tr>")
    }
  }
  paste(c(out,"</tbody></table>",paste0("<p class='note'>",mixedgp_html_escape(note),"</p>")),collapse="\n")
}

ocean_export_paper_table <- function(d,name,caption,dest,note="") {
  if(!nrow(d)) return(invisible(NULL))
  mixedgp_atomic_csv(d,ocean_path(dest,paste0(name,".csv")))
  ocean_atomic_lines(c(as.character(knitr::kable(mixedgp_table_markup(d,"pipe"),format="pipe",row.names=FALSE,align=mixedgp_table_align(d))),"",note),ocean_path(dest,paste0(name,".md")))
  tex_data <- mixedgp_table_markup(d,"latex")
  tex_data[] <- lapply(tex_data,ocean_tex_note)
  names(tex_data) <- ocean_tex_note(names(tex_data))
  tex <- as.character(knitr::kable(tex_data,format="latex",booktabs=TRUE,caption=caption,label=name,
    escape=FALSE,row.names=FALSE,align=mixedgp_table_align(d)))
  tex <- paste0(strsplit(tex,"\\end{table}",fixed=TRUE)[[1]][1],
    "\\par\\smallskip\n\\begin{minipage}{\\linewidth}\\footnotesize\n",ocean_tex_note(note),
    "\n\\end{minipage}\n\\end{table}")
  ocean_atomic_lines(tex,ocean_path(dest,paste0(name,".tex")))
  html <- mixedgp_table_html(d,caption,note)
  ocean_atomic_lines(c("<!doctype html><meta charset='utf-8'>",
    "<style>body{font:15px/1.5 sans-serif;margin:24px}table{border-collapse:collapse}td,th{padding:7px 10px;border-bottom:1px solid #ddd}caption{text-align:left;font-weight:600}.note{color:#555}</style>",html),
    ocean_path(dest,paste0(name,".html")))
  invisible(html)
}

ocean_export_prediction_report <- function(report,dest) {
  mixedgp_atomic_csv(report$fold_metrics,ocean_path(dest,"fold_metrics_raw.csv"))
  mixedgp_atomic_csv(report$summary,ocean_path(dest,"summary_raw.csv"))
  mixedgp_atomic_csv(report$fold_metrics_SD_Y,ocean_path(dest,"fold_metrics_SD_Y.csv"))
  mixedgp_atomic_csv(report$availability,ocean_path(dest,"eiv_availability.csv"))
  mixedgp_save_plot_jobs(report$plots,dest,function(name)c(10,6))
  invisible(report)
}

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

ocean_prediction_evidence <- function(output_dir, inputs, report, profile="paper",
  dest=ocean_path(output_dir,"prediction_report","evidence")) {
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
      saved <- load_saved(ocean_path(output_dir,root,sprintf("fold%02d",r$fold),"baseline_package_outputs.rds"))
      draws <- saved$draws[[r$method]]; point <- saved$means[[r$method]]
    } else if(r$variant=="EIV") {
      s <- sources[sources$fold==r$fold & sources$cell==cell,]
      path <- ocean_path(output_dir,s$source)
      saved <- load_saved(ocean_path(path,"prediction_distributions.rds"))
      stopifnot(identical(as.logical(saved$known),as.logical(known)))
      draws <- saved[[if(r$prediction_mode=="class_only") "xc" else "mixed"]]
      pred <- read.csv(ocean_path(path,if(r$prediction_mode=="class_only") "predictions.csv" else "mixed_predictions.csv"))
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
  oof <- mixedgp_bind(rows)
  keys <- paste(oof$cell,oof$variant,oof$method,oof$prediction_mode,sep="/")
  pooled <- lapply(split(seq_len(nrow(oof)),keys),function(ids) {
    z <- oof[ids,]; stopifnot(!anyDuplicated(z$cohort_id))
    fold_rmse <- vapply(split(z$squared_error,z$fold),function(v)sqrt(mean(v)),0)
    data.frame(z[1,c("cell","variant","method","prediction_mode")],
      n_folds=length(fold_rmse),N=nrow(z),mean_fold_RMSE=mean(fold_rmse),
      pooled_RMSE=sqrt(mean(z$squared_error)),pooled_CRPS=mean(z$CRPS))
  })
  pooled <- mixedgp_bind(pooled)
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
  ans <- list(oof=oof,pooled=pooled,subgroups=mixedgp_bind(subgroups),
    information=mixedgp_bind(info),paired=mixedgp_bind(pairs),SD_Y=report$SD_Y)
  for(n in c("oof","pooled","subgroups","information","paired"))
    if(nrow(ans[[n]])) mixedgp_atomic_csv(ans[[n]],ocean_path(dest,paste0(n,".csv")))
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

ocean_cc_prediction_tables <- function(output_dir,inputs,legacy,cc_output_dir="",
  profile="paper",dest=ocean_path(output_dir,"prediction_report","cc_revision"),
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
  add(ac,"AC",paste0("OI: ",inputs$design$n_train," (x,c)"),"100% (x,c)")
  for(cell in c("O0","O20","O50")) {
    z<-base[base$variant=="EIV" & base$cell==cell & base$prediction_mode=="class_only",]
    if(nrow(z))add(z,if(cell=="O0")"AC" else "EIV all-profile",
      paste0(sub("O","",cell),paste0("% calibration; ",inputs$design$n_train," rows")),"100% (x,c)")
    if(cell!="O0") {
      z<-base[base$variant=="EIV" & base$cell==cell & base$prediction_mode=="mixed",]
      if(nrow(z))add(z,"EIV all-profile",paste0(sub("O","",cell),paste0("% calibration; ",inputs$design$n_train," rows")),
        if(cell=="O20")"20% (x,u) + 80% (x,c)" else "50% (x,u) + 50% (x,c)")
    }
  }
  cc<-details<-list()
  if(nzchar(cc_output_dir)) {
    tasks<-ocean_task_catalogue("paper");tasks<-tasks[tasks$module=="complete_case_baselines" & tasks$fold %in% expected_folds,]
    for(i in seq_len(nrow(tasks))) {
      t<-tasks[i,];dir<-ocean_path(cc_output_dir,"complete_case_baselines/fits",t$task_id)
      if(!file.exists(ocean_path(dir,"metrics.csv")))next
      if(!ocean_verify_receipt(cc_output_dir,t$task_id))stop("Unverified CC task: ",t$task_id)
      design<-readRDS(ocean_path(dir,"design.rds"));d<-inputs$folds[[paste(t$fold,t$cell,sep="/")]]
      stopifnot(identical(as.character(design$training_ids),as.character(d$train$cohort_id[d$observed_train])),
        identical(as.character(design$test_ids),as.character(d$test$cohort_id[d$observed_test_mixed])),
        identical(design$protocol,"native-CC-observed-u-subset-v1"))
      z<-read.csv(ocean_path(dir,"metrics.csv"));z$source_profile<-design$profile
      status<-read.csv(ocean_path(dir,"package_status.csv"))
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
  cc<-mixedgp_bind(cc)
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
      add(z,"CC matched subset",paste0(sub("O","",cell),"% CC: ",inputs$design$train_exact[[cell]]," rows"),
        paste0(inputs$design$test_exact[[cell]]," (x,u,c)"),sp)
      s<-roots[roots$cell==cell & roots$fold==f,];dir<-ocean_path(output_dir,s$source)
      saved<-readRDS(ocean_path(dir,"prediction_distributions.rds"));k<-saved$known
      pred<-read.csv(ocean_path(dir,"mixed_predictions.csv"));d<-inputs$folds[[paste(f,cell,sep="/")]]
      stopifnot(identical(as.logical(k),as.logical(d$observed_test_mixed)),
        identical(as.character(pred$cohort_id),as.character(d$test$cohort_id)))
      draws<-saved$mixed[,k,drop=FALSE]
      for(a in c("conditional_means","conditional_vars"))attr(draws,a)<-attr(saved$mixed,a)[,k,drop=FALSE]
      z<-score_env$probability_score(draws,pred$predicted_Y[k],d$test$Y[k],f,"EIV-GP",cell)
      # Preserve the actual computation profile separately from display labels.
      add(z,"CC matched subset",paste0(sub("O","",cell),paste0("% calibration; ",inputs$design$n_train," rows")),
        paste0(sum(k)," (x,u)"),profile)
    }
  }
  raw<-mixedgp_bind(rows);summary<-list()
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
  summary<-mixedgp_bind(summary)
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  mixedgp_atomic_csv(raw,ocean_path(dest,"fold_metrics.csv"));mixedgp_atomic_csv(summary,ocean_path(dest,"summary.csv"))
  pooled_cc<-summary[summary$scope=="CC matched subset",c("cell","method","source_profile","n_folds","fold_ids","complete_folds","N_total","RMSE_mean","RMSE_sd","pooled_RMSE")]
  mixedgp_atomic_csv(pooled_cc,ocean_path(dest,"cc_pooled_rmse.csv"))
  list(raw=raw,summary=summary,cc=cc,dest=dest,pooled_cc=pooled_cc,computation_details=mixedgp_bind(details),expected_folds=expected_folds)
}

ocean_prediction_labels <- function(xc="OI",mixed="MI",exact="MI",cc="CC") {
  labels<-c(xc=xc,mixed=mixed,exact=exact,cc=cc)
  stopifnot(all(nzchar(labels)),!anyNA(labels));labels
}

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
  cat(mixedgp_table_html(out,"Pooled RMSE on matched observed-u test profiles",
    "Pooled RMSE is sqrt(sum(n_fold * RMSE_fold^2) / sum(n_fold)). Total test N and fold counts reflect available matched results. With one fold, pooled and fold RMSE coincide. Bold uses the same matching rule as the main tables."))
  if(nrow(review$computation_details))cat(mixedgp_table_html(review$computation_details,
    "Native complete-case fitting controls and optimizer status",
    "Controls are read from each saved fit design; termination codes are retained from the native packages. EIV numerical settings and posterior diagnostics are described in Section 3."))
}

ocean_diag_require <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

ocean_diag_csv <- function(path, columns) {
  ocean_diag_require(file.exists(path), paste("Missing input:", path))
  x <- utils::read.csv(path, stringsAsFactors = FALSE)
  ocean_diag_require(all(columns %in% names(x)), paste("Unexpected columns:", path))
  x
}

ocean_diag_flag <- function(x) {
  ocean_diag_require(!anyNA(x) && all(x %in% c(0, 1, FALSE, TRUE)),
    "Observation masks must contain only 0/1 or FALSE/TRUE.")
  as.logical(x)
}

ocean_diag_inputs <- function(data_dir) {
  z<-ocean_prepare_inputs(data_dir)
  train<-lapply(1:5,function(f)z$folds[[paste(f,"O50",sep="/")]]$train)
  test<-lapply(1:5,function(f)z$folds[[paste(f,"O50",sep="/")]]$test)
  list(cohort=z$cohort,train=train,test=test,design=z$design,
    audit=do.call(rbind,lapply(1:5,function(f)data.frame(fold=f,n_train=nrow(train[[f]]),
      n_test=nrow(test[[f]]),n_calib_O20=sum(train[[f]]$calib_o20),
      n_calib_O50=sum(train[[f]]$calib_o50),n_common_missing=sum(!train[[f]]$calib_o50)))))
}

ocean_diag_sources <- function(output_dir) {
  selected<-ocean_final_eiv_sources(output_dir,ocean_saved_profile(output_dir))
  ocean_diag_require(all(selected$available),"Final fits selected from current initial diagnostics are incomplete.")
  sources<-selected[,c("fold","cell","source","continued","fit_file")]
  for(i in seq_len(nrow(sources))) {
    s<-sources[i,];path<-ocean_path(output_dir,s$source)
    for(file in c(s$fit_file,"predictive_metrics.csv","prediction_distributions.rds",
      "U_posterior_intervals.csv","convergence_summary.csv"))
      ocean_diag_require(file.exists(ocean_path(path,file)),paste("Incomplete paper output:",path,file))
  }
  manifest <- ocean_path(output_dir, "continuation_plus20k", "final_review_20260927", "fit_sources.csv")
  if (file.exists(manifest)) {
    recorded <- ocean_diag_csv(manifest, c("fold", "cell", "source", "continued"))
    i <- match(paste(sources$fold, sources$cell), paste(recorded$fold, recorded$cell))
    ocean_diag_require(nrow(recorded) == 15L && !anyNA(i) &&
      identical(sources$source, ocean_result_relative(output_dir,recorded$source[i])), "Final-review source manifest disagrees with the paper run.")
  }
  sources
}

ocean_diag_five_metrics <- function(d) {
  error <- d$posterior_mean - d$reference_U
  data.frame(n = nrow(d), Bias = mean(error), RMSE = sqrt(mean(error^2)),
    MAE = mean(abs(error)), Coverage95 = mean(d$reference_U >= d$q025 & d$reference_U <= d$q975),
    Width95 = mean(d$q975 - d$q025))
}

ocean_diag_fold_summary <- function(d, groups, metrics, require_all=TRUE) {
  # Fold SD is descriptive variability across overlapping CV training sets.
  # It is neither posterior SD nor an independent-replicate Monte Carlo SE.
  key <- do.call(paste, c(d[groups], sep = "/"))
  do.call(rbind, lapply(split(seq_len(nrow(d)), key), function(i) {
    z <- d[i, , drop = FALSE]
    ocean_diag_require(!anyDuplicated(z$fold) && all(z$fold %in% 1:5) &&
      (!require_all || (nrow(z)==5L && setequal(z$fold,1:5))), "Expected one score per fold.")
    do.call(rbind, lapply(metrics, function(metric) {
      ocean_diag_require(all(is.finite(z[[metric]])), paste("Non-finite score:", metric))
      data.frame(z[1, groups, drop = FALSE], metric = metric, mean = mean(z[[metric]]),
        sd = stats::sd(z[[metric]]), n_folds = nrow(z), row.names = NULL)
    }))
  }))
}

ocean_diag_u_display <- function(output_dir, inputs, sources) {
  parts <- scales <- counts <- list()
  for (i in seq_len(nrow(sources))) {
    s <- sources[i, ]; tr <- inputs$train[[s$fold]]
    observed <- switch(s$cell, O0 = rep(FALSE, nrow(tr)),
      O20 = tr$calib_o20, O50 = tr$calib_o50)
    z <- ocean_diag_csv(ocean_path(output_dir, s$source, "U_posterior_intervals.csv"),
      c("cohort_id", "C", "posterior_mean", "q025", "q975", "scale"))
    expected_scale <- if (s$cell == "O0") "model" else "raw"
    ocean_diag_require(!anyDuplicated(z$cohort_id) &&
      setequal(z$cohort_id, tr$cohort_id[!observed]) && all(z$scale == expected_scale),
      "U display rows or scales disagree with the fitted observation mask.")
    # The full TRAINING reference mean/SD are an external plotting convention.
    # Hidden U is used here only for evaluation/axes, never to fit/align O0.
    center <- mean(tr$U); scale <- stats::sd(tr$U)
    ocean_diag_require(is.finite(scale) && scale > 0, "Invalid training reference SD.")
    d <- data.frame(fold = s$fold, cell = s$cell, cohort_id = tr$cohort_id,
      C = tr$C, observed = observed, common_hidden = !tr$calib_o50,
      reference_U = tr$U, posterior_mean = NA_real_, q025 = NA_real_, q975 = NA_real_,
      posterior_scale = expected_scale, reference_center = center, reference_sd = scale,
      plot_reference = (tr$U - center) / scale,
      plot_estimate = NA_real_, plot_q025 = NA_real_, plot_q975 = NA_real_)
    j <- match(z$cohort_id, tr$cohort_id)
    ocean_diag_require(!anyNA(j) && all(z$C == tr$C[j]) &&
      all(is.finite(as.matrix(z[, c("posterior_mean", "q025", "q975")]))) &&
      all(z$q975 > z$q025), "Invalid saved U posterior display intervals.")
    d$posterior_mean[j] <- z$posterior_mean
    d$q025[j] <- z$q025; d$q975[j] <- z$q975
    # O0 is left in its ORIGINAL model coordinate. Standardizing its posterior
    # or fitting a truth-based affine alignment would conceal scale ambiguity.
    shift <- if (s$cell == "O0") 0 else center
    divisor <- if (s$cell == "O0") 1 else scale
    d$plot_estimate[j] <- (z$posterior_mean - shift) / divisor
    d$plot_q025[j] <- (z$q025 - shift) / divisor
    d$plot_q975[j] <- (z$q975 - shift) / divisor
    d$plot_estimate[observed] <- d$plot_reference[observed]
    # NA posterior summaries on observed rows prevent them being mistaken for
    # inferred values or entering imputation scores. Only display values coincide.
    parts[[i]] <- d
    scales[[i]] <- data.frame(fold = s$fold, cell = s$cell,
      reference_center = center, reference_sd = scale, posterior_scale = expected_scale,
      interpretation = if (s$cell == "O0") "Unanchored model coordinate; illustrative only" else
        "Raw posterior transformed to training-reference SD units")
    counts[[i]] <- data.frame(fold = s$fold, cell = s$cell, n_train = nrow(tr),
      n_masked = sum(!observed), n_observed = sum(observed), n_common_hidden_scored = if (s$cell == "O0") 0L else sum(!tr$calib_o50))
  }
  list(points = do.call(rbind, parts), scales = do.call(rbind, scales), counts = do.call(rbind, counts))
}

ocean_diag_u <- function(output_dir, inputs, sources) {
  points <- list()
  for (i in which(sources$cell != "O0")) {
    s <- sources[i, ]; tr <- inputs$train[[s$fold]]
    z <- ocean_diag_csv(ocean_path(output_dir, s$source, "U_posterior_intervals.csv"),
      c("cohort_id", "C", "posterior_mean", "q025", "q975", "scale"))
    mask <- tr[[if (s$cell == "O20") "calib_o20" else "calib_o50"]]
    ocean_diag_require(!anyDuplicated(z$cohort_id) && setequal(z$cohort_id, tr$cohort_id[!mask]) &&
      all(z$scale == "raw"), "Saved U intervals do not match the training observation mask.")
    # Nested masks imply that !calib_o50 is hidden under BOTH fitted models.
    ref <- tr[!tr$calib_o50, ]; z <- z[match(ref$cohort_id, z$cohort_id), ]
    ocean_diag_require(nrow(z) == sum(!tr$calib_o50) && all(z$C == ref$C) &&
      all(is.finite(as.matrix(z[, c("posterior_mean", "q025", "q975")]))) &&
      all(z$q975 > z$q025), "Invalid common-missing U intervals.")
    points[[length(points) + 1L]] <- data.frame(fold = s$fold, cell = s$cell,
      cohort_id = ref$cohort_id, C = ref$C, reference_U = ref$U,
      posterior_mean = z$posterior_mean, q025 = z$q025, q975 = z$q975)
  }
  points <- do.call(rbind, points)
  scores <- list()
  for (fold in 1:5) for (cell in c("O20", "O50")) for (C in 0:4) {
    z <- points[points$fold == fold & points$cell == cell & (C == 0 | points$C == C), ]
    scores[[length(scores) + 1L]] <- data.frame(fold = fold, cell = cell,
      group = if (C == 0) "overall" else paste0("class", C), ocean_diag_five_metrics(z))
  }
  scores <- do.call(rbind, scores)
  display <- ocean_diag_u_display(output_dir, inputs, sources)
  list(points = points, fold_metrics = scores, summary = ocean_diag_fold_summary(scores,
    c("cell", "group"), c("Bias", "RMSE", "MAE", "Coverage95", "Width95")),
    display_points = display$points, display_scales = display$scales, display_counts = display$counts)
}

ocean_diag_grids <- function(cohort) {
  # Fixed S=34.20. Within each recorded band, select the 25
  # profiles closest in salinity; use their central 80% log-pressure range.
  # Ties are resolved by frozen cohort ID, without looking at Y.
  parts <- neighbors <- list()
  for (C in 1:4) {
    d <- cohort[cohort$C == C, ]
    near <- d[order(abs(d$G2salinity - 34.20), d$cohort_id)[1:25], ]
    bounds <- stats::quantile(near$log_pressure, c(.1, .9), type = 7)
    lp <- seq(bounds[1], bounds[2], length.out = 9L)
    parts[[C]] <- data.frame(design = "candidate_A", C = C,
      G2salinity = 34.20, log_pressure = lp, pressure = exp(lp))
    neighbors[[C]] <- data.frame(design = "candidate_A", class = C, cohort_id = near$cohort_id)
  }
  m <- do.call(rbind, parts)
  # Selected fixed X: S=34.25, P=650. Reference U is used ONLY to set an
  # offline illustration domain and rug; it is never supplied as calibration.
  center <- c(34.25, log(650))
  x <- as.matrix(cohort[, c("G2salinity", "log_pressure")])
  distance <- sqrt(rowSums(sweep(sweep(x, 2, center, "-"), 2, apply(x, 2, stats::sd), "/")^2))
  near <- cohort[order(distance, cohort$cohort_id)[1:30], ]
  bounds <- stats::quantile(near$U, c(.1, .9), type = 7)
  f <- data.frame(G2salinity = center[1], log_pressure = center[2], pressure = 650,
    U = seq(bounds[1], bounds[2], length.out = 41L))
  same_x <- data.frame(design = "same_x", C = 1:4, G2salinity = center[1],
    log_pressure = center[2], pressure = 650)
  m <- rbind(m, same_x); m$grid_id <- sprintf("m_%02d", seq_len(nrow(m)))
  f$grid_id <- sprintf("f_%02d", seq_len(nrow(f)))
  list(m = m, f = f, near = near, m_neighbors = do.call(rbind, neighbors),
    support = data.frame(S = 34.25, pressure = 650, n_neighbors = 30L,
      U_low = bounds[1], U_high = bounds[2],
      C1 = sum(near$C == 1L), C2 = sum(near$C == 2L), C3 = sum(near$C == 3L), C4 = sum(near$C == 4L)))
}

ocean_diag_balanced_ids <- function(fit, n_per_chain) {
  lens <- vapply(fit$mcmc$samples_by_chain$u, nrow, integer(1))
  ocean_diag_require(length(lens) == 4L && all(lens >= n_per_chain),
    "Final diagnostics require four retained paper chains.")
  offsets <- c(0L, head(cumsum(lens), -1L))
  unlist(Map(function(n, off) off + unique(round(seq(1, n, length.out = n_per_chain))), lens, offsets),
    use.names = FALSE)
}

ocean_diag_point_summary <- function(draws, means, grid) {
  ocean_diag_require(is.matrix(draws) && is.matrix(means) && ncol(draws) == nrow(grid) &&
    ncol(means) == nrow(grid) && all(is.finite(draws)) && all(is.finite(means)),
    "Package returned invalid function draws.")
  q <- apply(draws, 2L, stats::quantile, probs = c(.025, .975), type = 7)
  # Average conditional means for the point estimate; do not add GP simulation
  # noise to this estimate. Quantiles include conditional GP uncertainty.
  data.frame(grid, posterior_mean = colMeans(means), q025 = q[1, ], q975 = q[2, ],
    Width95 = q[2, ] - q[1, ])
}

ocean_diag_functionals <- function(i, sources, output_dir, inputs, grids, settings,
  cache_dir, overwrite) {
  s <- sources[i, ]; path <- ocean_path(output_dir, s$source)
  fit_path <- ocean_path(path, s$fit_file)
  prediction_path <- ocean_path(path, "prediction_distributions.rds")
  tr <- inputs$train[[s$fold]]
  mg <- grids$m[grids$m$design == "candidate_A" | s$cell != "O0", ]
  identity <- list(version = "ocean-diagnostics-1", package = as.character(utils::packageVersion("eivGP")),
    source = s, hashes = unname(tools::md5sum(c(fit_path, prediction_path))),
    train = tr[, c("cohort_id", "G2salinity", "log_pressure", "C", "Y", "U", "calib_o20", "calib_o50")],
    prediction_stream_schema = if(is.null(mixedgp_budget())) 0L else 1L,
    m_grid = mg, f_grid = grids$f, settings = settings,
    code = list(body(ocean_diag_functionals), body(ocean_diag_balanced_ids), body(ocean_diag_point_summary)))
  cache_path <- ocean_path(cache_dir, sprintf("fold%02d_%s.rds", s$fold, s$cell))
  if (file.exists(cache_path) && !overwrite) {
    saved <- readRDS(cache_path)
    saved$identity$source$source <- ocean_result_relative(output_dir,saved$identity$source$source)
    ocean_diag_require(identical(saved$identity, identity),
      paste("Diagnostic cache inputs changed. Choose a new report_dir or set overwrite=TRUE:", cache_path))
    return(saved$result)
  }
  message("Function diagnostics: fold ", s$fold, "/", s$cell)
  fit <- readRDS(fit_path)
  mask <- switch(s$cell, O0 = rep(FALSE, nrow(tr)), O20 = tr$calib_o20, O50 = tr$calib_o50)
  xmat <- function(d) as.matrix(d[, c("G2salinity", "log_pressure"), drop = FALSE])
  ocean_diag_require(inherits(fit, "eivgp_fit") && identical(fit$sampler_version, "0.3.1") &&
    isTRUE(all.equal(unname(fit$data$x_raw), unname(xmat(tr)), tolerance = 1e-10)) &&
    isTRUE(all.equal(as.numeric(fit$data$y_raw), tr$Y, tolerance = 1e-10)) &&
    identical(as.integer(fit$data$c_ord), tr$C) &&
    setequal(fit$data$calib_idx, which(mask)) && setequal(fit$data$miss_idx, which(!mask)),
    "Fit object does not match the frozen training rows, model version, or calibration mask.")
  if (any(mask)) ocean_diag_require(
    isTRUE(all.equal(as.numeric(fit$data$u_obs_raw[mask]), tr$U[mask], tolerance = 1e-10)),
    "Exact oxygen measurements in the fitted object disagree with the calibration rows.")
  ids <- ocean_diag_balanced_ids(fit, settings$mean_states_per_chain)
  seed <- settings$seed + 100L * s$fold + match(s$cell, c("O0", "O20", "O50"))
  md <- eivGP::predict_eivgp(fit, new_X = xmat(mg), new_C = matrix(mg$C, ncol = 1),
    target = "mean", draw_ids = ids, n_latent = settings$n_latent,
    n_per_draw = settings$mean_gp_draws, include_process_uncertainty = TRUE,
    joint = FALSE, return_components = TRUE, seed = seed)
  conditional_means <- attr(md, "conditional_means", exact = TRUE)
  ms <- ocean_diag_point_summary(md, conditional_means, mg)
  ms$fold <- s$fold; ms$cell <- s$cell
  # Small numerical integration sensitivity check, not a recovery score. The
  # same 120 parameter states are evaluated with 256 and 512 latent nodes.
  check_columns <- match(1:4, mg$C)
  check_states <- unique(round(seq(1, length(ids), length.out = 120L)))
  check_grid <- mg[check_columns, ]
  alt <- eivGP::predict_eivgp(fit, new_X = xmat(check_grid),
    new_C = matrix(check_grid$C, ncol = 1), target = "mean", draw_ids = ids[check_states],
    n_latent = 2L * settings$n_latent, n_per_draw = 1L,
    include_process_uncertainty = FALSE, return_components = TRUE, seed = seed + 20000L)
  delta <- attr(alt, "conditional_means", exact = TRUE) -
    conditional_means[check_states, check_columns, drop = FALSE]
  check <- data.frame(fold = s$fold, cell = s$cell, check_grid,
    delta_mean_double_nodes = colMeans(delta), paired_state_mcse = apply(delta, 2, stats::sd) / sqrt(nrow(delta)))
  fs <- NULL; fids <- integer()
  if (s$cell != "O0") {
    fids <- readRDS(prediction_path)$draw_ids
    ocean_diag_require(length(fids) >= (if(ocean_saved_profile(output_dir)=="pilot") 20L else 1000L) && !anyDuplicated(fids) &&
      all(fids >= 1 & fids <= nrow(fit$mcmc$samples_tau)), "Unexpected saved posterior state IDs.")
    fm <- eivGP::predict_eivgp(fit, new_X = xmat(grids$f), new_U = matrix(grids$f$U, ncol = 1),
      new_U_scale = "raw", target = "surface", draw_ids = fids, include_process_uncertainty = FALSE)
    fd <- eivGP::predict_eivgp(fit, new_X = xmat(grids$f), new_U = matrix(grids$f$U, ncol = 1),
      new_U_scale = "raw", target = "surface", draw_ids = fids, n_per_draw = settings$surface_gp_draws,
      include_process_uncertainty = TRUE, joint = FALSE, seed = seed + 10000L)
    fs <- ocean_diag_point_summary(fd, fm, grids$f)
    fs$fold <- s$fold; fs$cell <- s$cell
  }
  result <- list(m = ms, f = fs, integration_check = check,
    draw_ids = list(mean = ids, surface = fids), seed = seed)
  saveRDS(list(identity = identity, result = result), cache_path)
  result
}

ocean_diag_prediction <- function(output_dir, inputs, sources) {
  metrics <- c("RMSE", "MAE", "CRPS", "NLPD", "Coverage95", "Width95", "IntervalScore95")
  rows <- masks <- list()
  for (fold in 1:5) {
    path <- ocean_path(output_dir, "full5_class4", sprintf("fold%02d", fold))
    expected <- ocean_available_baselines(output_dir,fold)
    if(length(expected)) {
      z <- ocean_diag_csv(ocean_path(path, "baseline_metrics.csv"), c("fold", "cell", "method", "n_test", metrics))
      ocean_diag_require(!anyDuplicated(z$method) && setequal(z$method,expected) &&
        all(z$fold == fold & z$n_test == nrow(inputs$test[[fold]])), "Unexpected baseline scores.")
      z$task <- "class_only"; z$test_exact_fraction <- 0
      rows[[length(rows) + 1L]] <- z[, c("fold", "cell", "method", "task", "test_exact_fraction", "n_test", metrics)]
    }
    mask <- ocean_diag_csv(ocean_path(path, "test_U_masks.csv"),
      c("cohort_id", "C", "observed_U_O20", "observed_U_O50"))
    te <- inputs$test[[fold]]
    ocean_diag_require(nrow(mask) == nrow(te) && !anyDuplicated(mask$cohort_id) &&
      setequal(mask$cohort_id, te$cohort_id), "Unexpected held-out oxygen mask rows.")
    mask <- mask[match(te$cohort_id, mask$cohort_id), ]
    m20 <- ocean_diag_flag(mask$observed_U_O20); m50 <- ocean_diag_flag(mask$observed_U_O50)
    expected<-ocean_test_masks(te,fold)
    ocean_diag_require(all(mask$C==te$C) && identical(m20,expected$observed_U_O20) &&
      identical(m50,expected$observed_U_O50),"Held-out masks changed.")
    masks[[fold]] <- data.frame(fold = fold, mask)
    for (i in which(sources$fold == fold)) {
      s <- sources[i, ]
      z <- ocean_diag_csv(ocean_path(output_dir, s$source, "predictive_metrics.csv"),
        c("fold", "cell", "method", "n_test", metrics))
      ocean_diag_require(nrow(z) == 2L && setequal(z$method, c("EIV-GP xc", "EIV-GP mixed")) &&
        all(z$fold == fold & z$cell == s$cell & z$n_test == nrow(inputs$test[[fold]])), "Unexpected EIV prediction scores.")
      saved <- readRDS(ocean_path(output_dir, s$source, "prediction_distributions.rds"))
      known <- switch(s$cell, O0 = rep(FALSE, nrow(te)), O20 = m20, O50 = m50)
      ocean_diag_require(identical(as.logical(saved$known), known) &&
        ncol(saved$xc) == nrow(te) && ncol(saved$mixed) == nrow(te) &&
        identical(saved$xc[, !known, drop = FALSE], saved$mixed[, !known, drop = FALSE]),
        "Saved mixed predictions disagree with the frozen exact-oxygen masks.")
      if (s$cell == "O0") z <- z[z$method == "EIV-GP xc", ] # O0 mixed is identical.
      z$task <- ifelse(z$method == "EIV-GP xc", "class_only", "matched_mixed")
      z$test_exact_fraction <- ifelse(z$task == "class_only", 0, if (s$cell == "O20") .2 else .5)
      z$method <- "EIV-GP"
      rows[[length(rows) + 1L]] <- z[, c("fold", "cell", "method", "task", "test_exact_fraction", "n_test", metrics)]
    }
  }
  raw <- do.call(rbind, rows)
  normalized <- raw
  s_y <- stats::sd(inputs$cohort$Y) # Scale by the SD of Y in the frozen cohort.
  for (metric in c("RMSE", "MAE", "CRPS", "Width95", "IntervalScore95"))
    normalized[[metric]] <- normalized[[metric]] / s_y
  normalized$NLPD <- normalized$NLPD - log(s_y) # Density transformation, not division.
  list(fold_metrics_raw = raw, fold_metrics_SD_Y = normalized,
    summary_raw = ocean_diag_fold_summary(raw, c("task", "method", "cell"), metrics,require_all=FALSE),
    summary_SD_Y = ocean_diag_fold_summary(normalized, c("task", "method", "cell"), metrics,require_all=FALSE),
    test_masks = do.call(rbind, masks), SD_Y = s_y)
}

ocean_diag_plots <- function(u, m, f, grids, representative_fold, baseline_same_x = NULL) {
  colors <- c(O0 = "#777777", O20 = "#2365a3", O50 = "#c76326")
  labels <- c(O0 = "0% exact training U", O20 = "20% exact training U", O50 = "50% exact training U")
  wrap <- function(x, width = 105L) paste(strwrap(x, width = width), collapse = "\n")
  style <- function(p, fill = FALSE) {
    p <- p + ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
        plot.caption = ggplot2::element_text(hjust = 0)) +
      ggplot2::scale_color_manual(values = colors, labels = labels, name = NULL)
    if (fill) p <- p + ggplot2::scale_fill_manual(values = colors, labels = labels, name = NULL)
    p
  }
  interval_note <- "Pointwise 95% credible intervals for the mean; response noise excluded."
  u_plot <- function(all_folds) {
    d <- u$display_points
    if (!all_folds) d <- d[d$fold == representative_fold, ]
    d$class <- factor(d$C, levels = 1:4)
    d$role <- factor(ifelse(d$observed, "Observed calibration (not imputed)", "Masked training U"),
      levels = c("Masked training U", "Observed calibration (not imputed)"))
    panel_labels<-setNames(vapply(c("O0","O20","O50"),function(cell) {
      z<-d[d$fold==d$fold[1] & d$cell==cell,]
      sprintf("%s%% calibration\n%d masked / %d observed\n%s",sub("O","",cell),sum(!z$observed),sum(z$observed),
        if(cell=="O0")"Unanchored model coordinate" else "Training-reference SD units")
    },""),c("O0","O20","O50"))
    d$calibration <- factor(d$cell, levels = names(panel_labels), labels = unname(panel_labels))
    d$fold_label <- factor(paste("Fold", d$fold), levels = paste("Fold", 1:5))
    p <- ggplot2::ggplot(d, ggplot2::aes(plot_reference, plot_estimate, color = class)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
      ggplot2::geom_linerange(data = d[!d$observed, ],
        ggplot2::aes(ymin = plot_q025, ymax = plot_q975),
        alpha = .20, linewidth = .35, show.legend = FALSE) +
      ggplot2::geom_point(ggplot2::aes(shape = role), size = 1.5, alpha = .85, stroke = .65) +
      ggplot2::scale_color_manual(values = c("1" = "#1b9e77", "2" = "#d95f02",
        "3" = "#7570b3", "4" = "#0072B2"), name = "Oxygen class", drop = FALSE) +
      ggplot2::scale_shape_manual(values = c("Masked training U" = 16,
        "Observed calibration (not imputed)" = 1), name = NULL) +
      ggplot2::guides(shape = ggplot2::guide_legend(order = 1, override.aes = list(color = "grey25", alpha = 1)),
        color = ggplot2::guide_legend(order = 2, override.aes = list(shape = 16, alpha = 1))) +
      ggplot2::coord_cartesian(xlim = range(d$plot_reference) + c(-.1, .1),
        ylim = range(c(d$plot_q025, d$plot_q975, d$plot_estimate), na.rm = TRUE) + c(-.1, .1)) +
      ggplot2::labs(x = "Reference oxygen (training-reference SD units)", y = "U estimate (display units)",
        title = "Training oxygen: masked estimates and observed calibration inputs",
        subtitle = if (all_folds) sprintf("All five folds; %d training profiles in every panel",sum(d$fold==d$fold[1] & d$cell=="O0")) else
          paste("Fixed fold", representative_fold, paste0("| ",sum(d$fold==representative_fold & d$cell=="O0")," training profiles in every panel")),
        caption = wrap(paste("Solid: masked training U with 95% posterior intervals. Hollow: exact calibration, not imputed.",
          "Reference is centered/scaled by each fold's full training reference U; O20/O50 posteriors use that same transformation.",
          "O0 stays on its unanchored model scale: identity-line agreement is illustrative, not physical-scale recovery.",
          "Training Y participates. Scores use the same profiles hidden under both calibrated conditions."))) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
        plot.caption = ggplot2::element_text(hjust = 0))
    if (all_folds) p + ggplot2::facet_grid(fold_label ~ calibration) else
      p + ggplot2::facet_wrap(~calibration, nrow = 1)
  }
  m_plot <- function(all_folds) {
    d <- m[m$design == "candidate_A", ]
    if (!all_folds) d <- d[d$fold == representative_fold, ]
    d$class <- factor(d$C, levels = 1:4, labels = paste("Class", 1:4))
    p <- ggplot2::ggplot(d, ggplot2::aes(pressure, posterior_mean, color = cell, fill = cell)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = q025, ymax = q975), alpha = .13, color = NA) +
      ggplot2::geom_line(linewidth = .65) +
      ggplot2::labs(x = "Pressure (dbar)", y = "Mean silicate m(x,c) (micromol/kg)",
        title = "Observable mean m(x,c) at salinity 34.20",
        subtitle = if (all_folds) "All five folds; each band belongs to its own fitted model" else
          paste("Fixed fold", representative_fold, "| candidate A"),
        caption = wrap(paste(interval_note, "Pressure ranges follow the prespecified local display rule.")))
    p <- if (all_folds) p + ggplot2::facet_wrap(~fold + class, ncol = 4, scales = "free") else
      p + ggplot2::facet_wrap(~class, ncol = 2, scales = "free")
    style(p, fill = TRUE)
  }
  fm_plot <- function(all_folds) {
    df <- f; dm <- m[m$design == "same_x", ]; db <- baseline_same_x
    if (!all_folds) {
      df <- df[df$fold == representative_fold, ]; dm <- dm[dm$fold == representative_fold, ]
      if(!is.null(db)) db <- db[db$fold == representative_fold, ]
    }
    df$model <- paste0("EIV-GP ",df$cell); dm$model <- paste0("EIV-GP ",dm$cell)
    means <- mixedgp_bind(list(dm,db))
    palette <- c("EIV-GP O20"="#2365a3","EIV-GP O50"="#c76326",
      "UC-GP (x,c)"="#7b3294","LVGP (x,c)"="#008837","EzGP (x,c)"="#a6761d")
    models <- names(palette)[names(palette) %in% means$model]
    offsets <- setNames(seq(-.26,.26,length.out=length(models)),models)
    means$class_position <- means$C + offsets[means$model]
    means$type <- ifelse(grepl("^EIV",means$model),"EIV posterior mean","Baseline fitted mean")
    # Match the physical response scale across panels; give the discrete panel
    # less width, and use larger points/capped intervals to reduce empty space.
    ylim <- range(c(df$q025,df$q975,dm$q025,dm$q975,means$posterior_mean),na.rm=TRUE)
    rug <- grids$near[grids$near$U >= min(df$U) & grids$near$U <= max(df$U), ]
    pf <- ggplot2::ggplot(df, ggplot2::aes(U, posterior_mean, color = model, fill = model)) +
      ggplot2::geom_vline(xintercept = c(50, 100, 150), linetype = 3, color = "grey65", linewidth = .4) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = q025, ymax = q975), alpha = .14, color = NA) +
      ggplot2::geom_line(linewidth = .9) +
      ggplot2::geom_rug(data = rug, ggplot2::aes(x = U), inherit.aes = FALSE,
        sides = "b", color = "grey45", alpha = .65) +
      ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=.025)) +
      ggplot2::coord_cartesian(xlim = range(df$U), ylim = ylim) +
      ggplot2::labs(x = "Oxygen (micromol/kg)", y = "Mean silicate (micromol/kg)", title = "(a) Calibrated response f(x,u)")
    pm <- ggplot2::ggplot(means, ggplot2::aes(class_position,posterior_mean,color=model)) +
      ggplot2::geom_errorbar(data=means[grepl("^EIV",means$model),],
        ggplot2::aes(ymin=q025,ymax=q975),width=.10,linewidth=.85,show.legend=FALSE) +
      ggplot2::geom_point(ggplot2::aes(shape=type),size=3,stroke=.8) +
      ggplot2::scale_shape_manual(values=c("EIV posterior mean"=16,"Baseline fitted mean"=17),guide="none") +
      ggplot2::scale_x_continuous(breaks=1:4,limits=c(.5,4.5),expand=c(0,0)) +
      ggplot2::coord_cartesian(ylim = ylim) +
      ggplot2::labs(x = "Recorded oxygen class", y = NULL, title = "(b) Observable mean m(x,c)")
    theme_fm <- function(p) p + ggplot2::theme_minimal(base_size=12) +
      ggplot2::scale_color_manual(values=palette,limits=models,breaks=models,drop=FALSE,name=NULL) +
      ggplot2::scale_fill_manual(values=palette,guide="none") +
      ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),legend.position="bottom",
        plot.title=ggplot2::element_text(size=12),legend.key.width=grid::unit(1.1,"cm")) +
      ggplot2::guides(color=ggplot2::guide_legend(nrow=1,override.aes=list(shape=16,linewidth=.8)))
    if (all_folds) { pf <- pf + ggplot2::facet_grid(fold ~ .); pm <- pm + ggplot2::facet_grid(fold ~ .) }
    patchwork::wrap_plots(theme_fm(pf)+ggplot2::guides(color="none"),
      theme_fm(pm)+ggplot2::guides(color=ggplot2::guide_legend(nrow=1,
        override.aes=list(shape=ifelse(grepl("^EIV",models),16,17),size=3))),
      widths=c(1.65,1),guides="collect") +
      patchwork::plot_annotation(title = "Salinity 34.25, pressure 650 dbar",
        subtitle = if (all_folds) "All five folds; separate fitted posteriors" else paste("Fixed fold", representative_fold),
        caption = paste("EIV bars/shading: pointwise 95% mean credible intervals. Triangles: original (x,c) baseline fitted means.",
          "\nm integrates the full fitted oxygen distribution. Dotted lines: recorded band boundaries; rug: nearby reference oxygen.")) &
      ggplot2::theme(legend.position="bottom",plot.caption=ggplot2::element_text(hjust=0,size=10))
  }
  list(u_main = u_plot(FALSE), m_main = m_plot(FALSE), fm_main = fm_plot(FALSE),
    u_all_folds = u_plot(TRUE), m_all_folds = m_plot(TRUE), fm_all_folds = fm_plot(TRUE))
}

ocean_build_diagnostics <- function(output_dir, data_dir = ocean_path(output_dir, "data"),
  report_dir = ocean_path(output_dir, "diagnostics"), representative_fold = 1L,
  overwrite = FALSE) {
  for (package in c("eivGP", "ggplot2", "patchwork"))
    ocean_diag_require(requireNamespace(package, quietly = TRUE), paste("Required package:", package))
  ocean_diag_require(utils::packageVersion("eivGP") >= "0.4.0.9004",
    "Install this checkout of eivGP before running the Ocean workflow.")
  ocean_diag_require(length(representative_fold) == 1L && representative_fold %in% 1:5,
    "representative_fold must be one of 1:5; the default main display uses fold 1.")
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  inputs <- ocean_diag_inputs(data_dir)
  sources <- ocean_diag_sources(output_dir)
  grids <- ocean_diag_grids(inputs$cohort)
  u <- ocean_diag_u(output_dir, inputs, sources)
  prediction <- ocean_diag_prediction(output_dir, inputs, sources)
  settings <- list(mean_states_per_chain = 60L, n_latent = 256L, mean_gp_draws = 20L,
    surface_gp_draws = 4L, seed = 202610080L)
  if(ocean_saved_profile(output_dir)=="pilot") settings<-list(
    mean_states_per_chain=8L,n_latent=16L,mean_gp_draws=4L,
    surface_gp_draws=2L,seed=202610080L)
  for (subdir in c("tables", "figures", "internal", "cache"))
    dir.create(ocean_path(report_dir, subdir), recursive = TRUE, showWarnings = FALSE)
  run_one <- function(i) ocean_diag_functionals(i, sources, output_dir, inputs, grids, settings,
    cache_dir = ocean_path(report_dir, "cache"), overwrite = overwrite)
  parts <- mixedgp_stage_lapply(seq_len(nrow(sources)), run_one,
    workers=mixedgp_stage_workers(mixedgp_budget_cores(),"functionals",nrow(sources)))
  m <- do.call(rbind, lapply(parts, `[[`, "m"))
  f <- do.call(rbind, lapply(parts, `[[`, "f"))
  checks <- do.call(rbind, lapply(parts, `[[`, "integration_check"))
  ocean_diag_require(nrow(m) == 15L * 36L + 10L * 4L && nrow(f) == 10L * 41L &&
    all(m$Width95 > 0) && all(f$Width95 > 0), "Incomplete functional diagnostics.")
  # Internal width averages are descriptive precision checks, not true m/f
  # recovery metrics or new main-text tables. Equal weights over the fixed grid.
  widths <- stats::aggregate(Width95 ~ fold + cell + C, m[m$design == "candidate_A", ], mean)
  width_summary <- ocean_diag_fold_summary(widths, c("cell", "C"), "Width95")
  diag <- lapply(seq_len(nrow(sources)), function(i) {
    s <- sources[i, ]; d <- utils::read.csv(ocean_path(output_dir, s$source, "convergence_summary.csv"))
    data.frame(fold = s$fold, cell = s$cell, d)
  })
  baseline_paths <- ocean_path(output_dir,"full5_class4",sprintf("fold%02d",1:5),"baseline_package_outputs.rds")
  baseline_same_x <- if(all(file.exists(baseline_paths)))
    ocean_baseline_mean_grid(output_dir,grids$m[grids$m$design=="same_x",]) else NULL
  plots <- ocean_diag_plots(u,m,f,grids,representative_fold,baseline_same_x)
  write_csv <- function(x, subdir, name) utils::write.csv(x, ocean_path(report_dir, subdir, name), row.names = FALSE)
  write_csv(u$points, "tables", "u_common_missing_points.csv")
  write_csv(u$display_points, "tables", "u_display_points.csv")
  write_csv(u$display_scales, "internal", "u_display_scales.csv")
  write_csv(u$display_counts, "internal", "u_display_counts.csv")
  write_csv(u$fold_metrics, "tables", "u_fold_metrics.csv")
  write_csv(u$summary, "tables", "u_summary.csv")
  write_csv(m, "tables", "m_points.csv"); write_csv(f, "tables", "f_points.csv")
  if(!is.null(baseline_same_x)) write_csv(baseline_same_x,"tables","baseline_m_same_x.csv")
  paper_u <- ocean_paper_u_tables(u$summary)
  ocean_export_paper_table(paper_u$main,"paper_u_imputation_main","Oxygen imputation: common masked training profiles",ocean_path(report_dir,"tables"),paper_u$note)
  ocean_export_paper_table(paper_u$class,"paper_u_imputation_class","Class-specific oxygen imputation",ocean_path(report_dir,"tables"),paper_u$note)
  write_csv(prediction$fold_metrics_raw, "tables", "y_prediction_fold_raw.csv")
  write_csv(prediction$fold_metrics_SD_Y, "tables", "y_prediction_fold_SD_Y.csv")
  write_csv(prediction$summary_raw, "tables", "y_prediction_summary_raw.csv")
  write_csv(prediction$summary_SD_Y, "tables", "y_prediction_summary_SD_Y.csv")
  write_csv(inputs$audit, "internal", "input_audit.csv")
  write_csv(sources, "internal", "fit_sources.csv")
  write_csv(grids$m, "internal", "m_grid.csv"); write_csv(grids$f, "internal", "f_grid.csv")
  write_csv(grids$m_neighbors, "internal", "m_display_neighbors.csv")
  write_csv(grids$near, "internal", "f_display_neighbors.csv")
  write_csv(grids$support, "internal", "f_display_support.csv")
  write_csv(checks, "internal", "m_integration_sensitivity.csv")
  write_csv(widths, "internal", "m_width_fold.csv")
  write_csv(width_summary, "internal", "m_width_summary.csv")
  write_csv(prediction$test_masks, "internal", "test_U_masks.csv")
  write_csv(do.call(rbind, diag), "internal", "mcmc_existing_summary.csv")
  dimensions <- list(u_main = c(13.2, 5.2), m_main = c(10, 7), fm_main = c(11, 4.9),
    u_all_folds = c(13.2, 14.5), m_all_folds = c(13, 13), fm_all_folds = c(11, 14))
  mixedgp_save_plot_jobs(plots,ocean_path(report_dir,"figures"),function(name)dimensions[[name]])
  provenance <- list(settings = settings, representative_fold = representative_fold,
    package = as.character(utils::packageVersion("eivGP")), upstream_commit = if(file.exists(ocean_path(output_dir,"workflow_identity.rds"))) readRDS(ocean_path(output_dir,"workflow_identity.rds"))$upstream_commit else NA_character_,
    input_hashes = tools::md5sum(ocean_path(data_dir, c("cohort_manifest.csv",
      sprintf("fold%02d_train.csv", 1:5), sprintf("fold%02d_test.csv", 1:5)))),
    source_hashes = lapply(seq_len(nrow(sources)), function(i) {
      s <- sources[i, ]
      tools::md5sum(ocean_path(output_dir, s$source, c(s$fit_file,
        "prediction_distributions.rds", "predictive_metrics.csv", "U_posterior_intervals.csv")))
    }),
    draw_ids = lapply(parts, `[[`, "draw_ids"), SD_Y = prediction$SD_Y,
    main_targets = c("training U imputation", "m(x,c)", "f(x,u) and same-X m(x,c)", "held-out Y prediction"))
  saveRDS(provenance, ocean_path(report_dir, "internal", "provenance.rds"))
  capture.output(utils::sessionInfo(), file = ocean_path(report_dir, "internal", "sessionInfo.txt"))
  invisible(list(u = u, m = m, f = f, prediction = prediction, grids = grids, plots = plots, baseline_same_x=baseline_same_x,
    integration_check = checks, sources = sources, report_dir = report_dir, provenance = provenance))
}

ocean_diag_display_table <- function(summary, groups, metrics, digits = 3L) {
  rows <- unique(summary[, groups, drop = FALSE])
  key <- function(d) do.call(paste, c(d[groups], sep = "/"))
  for (metric in metrics) {
    z <- summary[summary$metric == metric, ]
    i <- match(key(rows), key(z))
    rows[[metric]] <- sprintf(paste0("%.", digits, "f +/- %.", digits, "f"), z$mean[i], z$sd[i])
  }
  rownames(rows) <- NULL
  rows
}

ocean_native_mean <- function(method,fit,grid) {
  X <- as.matrix(grid[,c("G2salinity","log_pressure")]); C <- matrix(grid$C,ncol=1L)
  if(method=="UC-GP") {
    d <- eivGP:::make_mixedgp_competitor_data(X,C,4L)
    as.numeric(stats::predict(fit,newdata=d,type="UK",seCompute=FALSE,
      covCompute=FALSE,lightReturn=TRUE)$mean)
  } else {
    d <- eivGP:::make_mixedgp_numeric_matrix(X,C,4L)
    if(method=="LVGP") as.numeric(LVGP::LVGP_predict(d,fit,MSE_on=0)$Y_hat)
    else as.numeric(mixedgp_ezgp_predict(d,fit)$Y_hat)
  }
}

ocean_available_baselines <- function(output_dir,fold) {
  path <- ocean_path(output_dir,"full5_class4",sprintf("fold%02d",fold),"baseline_status.csv")
  if(!file.exists(path))return(character())
  status <- read.csv(path,stringsAsFactors=FALSE)
  stopifnot(all(c("method","status") %in% names(status)),!anyDuplicated(status$method),
    all(status$method %in% c("UC-GP","LVGP","EzGP")))
  status$method[status$status=="success"]
}

ocean_baseline_mean_grid <- function(output_dir,grid,folds=1:5) {
  out <- list()
  for(fold in folds) {
    saved <- readRDS(ocean_path(output_dir,"full5_class4",sprintf("fold%02d",fold),"baseline_package_outputs.rds"))
    for(method in ocean_available_baselines(output_dir,fold)) {
      mean <- ocean_native_mean(method,saved$fits[[method]],grid)
      stopifnot(length(mean)==nrow(grid),all(is.finite(mean)))
      out[[length(out)+1L]] <- data.frame(grid,posterior_mean=mean,fold,model=paste0(method," (x,c)"))
    }
  }
  if(!length(out))return(data.frame(grid[FALSE,,drop=FALSE],posterior_mean=numeric(),fold=integer(),model=character()))
  mixedgp_bind(out)
}

ocean_m_display_points <- function(cohort,grids,n_bins=6L,per_bin=3L) {
  out <- list()
  for(C in 1:4) {
    g <- grids$m[grids$m$design=="candidate_A" & grids$m$C==C,]
    ids <- grids$m_neighbors$cohort_id[grids$m_neighbors$class==C]
    d <- cohort[cohort$cohort_id %in% ids & cohort$log_pressure>=min(g$log_pressure) &
      cohort$log_pressure<=max(g$log_pressure),]
    d$pressure_bin <- cut(d$log_pressure,breaks=seq(min(g$log_pressure),max(g$log_pressure),length.out=n_bins+1L),
      include.lowest=TRUE,labels=FALSE)
    d$salinity_offset <- d$G2salinity-34.20
    keep <- unlist(lapply(split(seq_len(nrow(d)),d$pressure_bin),function(i)
      head(i[order(abs(d$salinity_offset[i]),d$cohort_id[i])],per_bin)),use.names=FALSE)
    out[[C]] <- d[keep,,drop=FALSE]
  }
  mixedgp_bind(out)
}

ocean_paper_u_tables <- function(summary) {
  format_u <- function(groups) {
    z <- summary[summary$group %in% groups,]
    keys <- unique(z[c("cell","group")]); keys <- keys[order(keys$cell,keys$group),]
    d <- data.frame(Calibration=paste0(sub("O","",keys$cell),"%"),check.names=FALSE)
    if(length(groups)>1L) d$Class <- sub("class","",keys$group)
    ranked <- as.data.frame(matrix(NA_real_,nrow(keys),2L))
    names(ranked) <- c("RMSE","MAE")
    eligible <- rep(TRUE,nrow(keys))
    for(metric in c("Bias","RMSE","MAE","Coverage95","Width95")) {
      v <- z[z$metric==metric,]; i <- match(paste(keys$cell,keys$group),paste(v$cell,v$group))
      mult <- if(metric=="Coverage95") 100 else 1
      label <- if(metric=="Coverage95") "Coverage (%)" else if(metric=="Width95") "Width" else metric
      d[[label]] <- ocean_mean_sd(mult*v$mean[i],mult*v$sd[i],if(metric=="Coverage95") 1 else 2)
      eligible <- eligible & v$n_folds[i]==5L
      if(metric %in% names(ranked))
        ranked[[metric]] <- if(metric=="Bias") abs(v$mean[i]) else v$mean[i]
      if(any(v$n_folds[i]!=5L)) d[["n folds"]] <- v$n_folds[i]
    }
    attr(d,"best") <- ocean_best_mask(ranked,keys$group,eligible)
    d
  }
  list(main=format_u("overall"),class=format_u(paste0("class",1:4)),
    note=paste("Entries are five-fold means (sample SDs), computed on the same common-hidden training profiles per fold.",
      "Bias, RMSE, MAE and width are in micromol/kg; coverage is percent for individual-U 95% posterior intervals.",
      "Bold marks smaller RMSE and MAE across calibration levels for the same target group, using unrounded values. Bias, coverage and width are descriptive; opposing class biases can cancel."))
}

ocean_surface_point_roles <- function(d,train) {
  i <- match(d$cohort_id,train$cohort_id)
  training <- !is.na(i)
  d$role <- ifelse(training,"Training","Held-out")
  d$calibration_O20 <- d$calibration_O50 <- rep(FALSE,nrow(d))
  d$calibration_O20[training] <- train$calib_o20[i[training]]
  d$calibration_O50[training] <- train$calib_o50[i[training]]
  levels <- c("Training: U masked in both","Training: U observed in O20/O50",
    "Training: U added in O50","Held-out")
  status <- rep(levels[1],nrow(d))
  status[d$calibration_O50] <- levels[3]
  status[d$calibration_O20] <- levels[2]
  status[!training] <- levels[4]
  d$point_status <- factor(status,levels=levels)
  d$point_outline <- c("grey45","#2365a3","#c76326","grey45")[as.integer(d$point_status)]
  d
}

ocean_surface_report <- function(output_dir,data_dir,diagnostics_dir=ocean_path(output_dir,"diagnostics"),
  representative_fold=1L,folds=1:5) {
  inputs <- ocean_diag_inputs(data_dir); grids <- ocean_diag_grids(inputs$cohort)
  m <- read.csv(ocean_path(diagnostics_dir,"tables/m_points.csv"))
  f <- read.csv(ocean_path(diagnostics_dir,"tables/f_points.csv"))
  baseline_all <- ocean_baseline_mean_grid(output_dir,grids$m,folds)
  baseline <- baseline_all[baseline_all$design=="candidate_A",]
  baseline_same <- baseline_all[baseline_all$design=="same_x",]
  selected_m <- ocean_m_display_points(inputs$cohort,grids)
  points_m <- points_f <- list()
  for(fold in folds) {
    training <- inputs$train[[fold]]
    dm <- ocean_surface_point_roles(selected_m,training)
    dm$fold <- fold; points_m[[length(points_m)+1L]] <- dm
    df <- grids$near[grids$near$U>=min(grids$f$U) & grids$near$U<=max(grids$f$U),]
    df <- ocean_surface_point_roles(df,training)
    df$fold <- fold; points_f[[length(points_f)+1L]] <- df
  }
  points_m <- mixedgp_bind(points_m); points_f <- mixedgp_bind(points_f)
  em <- m[m$design=="candidate_A" & m$fold %in% folds,]; em$model <- paste0("EIV-GP ",em$cell)
  lines <- mixedgp_bind(list(em,baseline))
  palette <- c("EIV-GP O0"="#777777","EIV-GP O20"="#2365a3","EIV-GP O50"="#c76326",
    "UC-GP (x,c)"="#7b3294","LVGP (x,c)"="#008837","EzGP (x,c)"="#a6761d")
  point_shapes <- c("Training: U masked in both"=16,"Training: U observed in O20/O50"=21,
    "Training: U added in O50"=21,"Held-out"=1)
  point_legend <- ggplot2::guide_legend(nrow=2,order=1,override.aes=list(
    color=c("grey45","#2365a3","#c76326","grey45"),fill="grey55",alpha=1,size=2.8,stroke=1))
  plot_m <- function(all=FALSE) {
    dl <- lines; de <- em; dp <- points_m
    if(!all) {dl<-dl[dl$fold==representative_fold,];de<-de[de$fold==representative_fold,];dp<-dp[dp$fold==representative_fold,]}
    dl$class_label <- paste("Class",dl$C); de$class_label <- paste("Class",de$C); dp$class_label <- paste("Class",dp$C)
    p <- ggplot2::ggplot(dl,ggplot2::aes(exp(log_pressure),posterior_mean,color=model)) +
      ggplot2::geom_ribbon(data=de,ggplot2::aes(ymin=q025,ymax=q975,fill=model),alpha=.08,color=NA) +
      ggplot2::geom_line(linewidth=.8) +
      ggplot2::geom_point(data=dp,ggplot2::aes(G2pressure,Y,shape=point_status,color=I(point_outline)),inherit.aes=FALSE,
        fill="grey55",size=2.5,alpha=.85,stroke=1,show.legend=c(shape=TRUE,colour=FALSE)) +
      ggplot2::scale_color_manual(values=palette,name=NULL) +
      ggplot2::scale_fill_manual(values=palette,guide="none") +
      ggplot2::scale_shape_manual(values=point_shapes,limits=names(point_shapes),drop=FALSE,name="Nearby profiles") +
      ggplot2::guides(color=ggplot2::guide_legend(nrow=2,order=2),shape=point_legend) +
      ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=.05)) +
      ggplot2::labs(x="Pressure (dbar)",y="Silicate (micromol/kg)",
        title="Observable mean m(x,c): fitted means and nearby observations",
        subtitle="Fixed salinity 34.20; separate oxygen classes",
        caption="Points: up to 3 nearest-salinity profiles per pressure bin, within the displayed curve range.\nBlue rings: training calibration in O20 and O50; orange rings: training calibration added in O50.\nEIV shading: pointwise 95% mean credible intervals. Baseline lines: original (x,c) fitted means.") +
      ggplot2::theme_minimal(base_size=12) + ggplot2::theme(legend.position="bottom",
        panel.grid.minor=ggplot2::element_blank(),plot.caption=ggplot2::element_text(hjust=0,size=10))
    if(all) p+ggplot2::facet_wrap(~fold+class_label,ncol=4,scales="free") else p+ggplot2::facet_wrap(~class_label,ncol=2,scales="free")
  }
  plot_f <- function(all=FALSE) {
    df <- f[f$fold %in% folds,]; dp <- points_f
    if(!all) {df<-df[df$fold==representative_fold,];dp<-dp[dp$fold==representative_fold,]}
    p <- ggplot2::ggplot(df,ggplot2::aes(U,posterior_mean,color=cell,fill=cell)) +
      ggplot2::geom_vline(xintercept=c(50,100,150),linetype=3,color="grey60") +
      ggplot2::geom_ribbon(ggplot2::aes(ymin=q025,ymax=q975),alpha=.12,color=NA) +
      ggplot2::geom_line(linewidth=.7) +
      ggplot2::geom_point(data=dp,ggplot2::aes(U,Y,shape=point_status,color=I(point_outline)),inherit.aes=FALSE,
        fill="grey55",size=2.5,alpha=.85,stroke=1,show.legend=c(shape=TRUE,colour=FALSE)) +
      ggplot2::scale_color_manual(values=c(O20="#2365a3",O50="#c76326"),name="Training calibration") +
      ggplot2::scale_fill_manual(values=c(O20="#2365a3",O50="#c76326"),guide="none") +
      ggplot2::scale_shape_manual(values=point_shapes,limits=names(point_shapes),drop=FALSE,name="Nearby profiles") +
      ggplot2::guides(color=ggplot2::guide_legend(order=2),shape=point_legend) +
      ggplot2::labs(x="Oxygen (micromol/kg)",y="Silicate (micromol/kg)",
        title="Calibrated response f(x,u): salinity 34.25, pressure 650 dbar",
        caption="Points: nearby profiles; salinity/pressure vary. Masked training U is shown only as an offline reference.\nBlue rings: training calibration in O20 and O50; orange rings: training calibration added in O50.\nShading: pointwise 95% credible intervals for the mean, excluding response noise.") +
      ggplot2::theme_minimal(base_size=11) + ggplot2::theme(legend.position="bottom",
        plot.caption=ggplot2::element_text(hjust=0,size=10))
    if(all) p+ggplot2::facet_wrap(~fold,ncol=2) else p
  }
  plots <- list(m_comparison=plot_m(),f_with_observations=plot_f(),
    m_comparison_all_folds=plot_m(TRUE),f_with_observations_all_folds=plot_f(TRUE))
  dest <- ocean_path(diagnostics_dir,"surface_comparison")
  mixedgp_atomic_csv(baseline,ocean_path(dest,"baseline_m_points.csv"))
  mixedgp_atomic_csv(baseline_same,ocean_path(dest,"baseline_m_same_x.csv"))
  mixedgp_atomic_csv(points_m,ocean_path(dest,"nearby_m_observations.csv"))
  mixedgp_atomic_csv(selected_m,ocean_path(dest,"m_selection_audit.csv"))
  mixedgp_atomic_csv(points_f,ocean_path(dest,"nearby_f_observations.csv"))
  mixedgp_save_plot_jobs(plots,dest,function(n)c(
    if(grepl("fm_",n)||grepl("m_comparison",n))11 else 10,
    if(grepl("all_folds",n))13 else if(grepl("fm_",n))4.9 else 7))
  list(plots=plots,baseline=baseline,baseline_same_x=baseline_same,points_m=points_m,points_f=points_f)
}

ocean_run_diagnostics <- function(output_dir,data_dir,
  report_dir=ocean_path(output_dir,"diagnostics"),
  representative_fold=1L,cache_from="",
  build_prediction=TRUE,build_surfaces=TRUE,prediction_report=NULL) {
  for(p in c("eivGP","ggplot2","patchwork","digest"))
    ocean_diag_require(requireNamespace(p,quietly=TRUE),paste("Required package:",p))
  output_dir <- normalizePath(output_dir,mustWork=TRUE)
  data_dir <- normalizePath(data_dir,mustWork=TRUE)
  dir.create(report_dir,recursive=TRUE,showWarnings=FALSE)
  report_dir <- normalizePath(report_dir,mustWork=TRUE)
  progress <- function(stage) ocean_atomic_lines(c(paste("stage:",stage),
    paste("updated:",format(Sys.time(),"%Y-%m-%d %H:%M:%S %z")),
    paste("pid:",Sys.getpid()),"Postprocessing saved fits only; no MCMC or baseline refit."),
    ocean_path(report_dir,"PROGRESS.txt"))
  progress("input and source audit")
  inputs <- ocean_prepare_inputs(data_dir)
  sources <- ocean_diag_sources(output_dir)
  copied <- character()
  if(nzchar(cache_from)) {
    from <- normalizePath(cache_from,mustWork=TRUE)
    dir.create(ocean_path(report_dir,"cache"),showWarnings=FALSE)
    for(n in sprintf("fold%02d_%s.rds",sources$fold,sources$cell)) {
      source <- ocean_path(from,n); target <- ocean_path(report_dir,"cache",n)
      if(file.exists(source) && !file.exists(target)) {
        stopifnot(file.copy(source,target)); copied <- c(copied,n)
      }
    }
  }
  progress("U imputation and cached package m/f calculations")
  diagnostics <- ocean_build_diagnostics(output_dir,data_dir,report_dir=report_dir,
    representative_fold=representative_fold)
  prediction <- prediction_report; surfaces <- NULL
  if(build_prediction) {
    progress("saved response prediction summaries")
    prediction <- ocean_prediction_report(output_dir,inputs,
      report_dir=ocean_path(report_dir,"prediction"))
  }
  if(!build_prediction && !is.null(prediction))
    ocean_export_prediction_report(prediction,ocean_path(report_dir,"prediction"))
  if(build_surfaces) {
    progress("native baseline m comparison and nearby measured-profile overlays")
    surfaces <- ocean_surface_report(output_dir,data_dir,report_dir,
      representative_fold=representative_fold)
  }
  main <- c("U imputation (O0/O20/O50, hollow calibration)"="figures/u_main.png",
    "m: native baselines and pressure-balanced nearby profiles"="surface_comparison/m_comparison.png",
    "f: S=34.25, P=650, measured nearby profiles"="surface_comparison/f_with_observations.png")
  main <- main[file.exists(ocean_path(report_dir,main))]
  ut <- ocean_paper_u_tables(diagnostics$u$summary)
  ocean_atomic_lines(c("# Ocean diagnostics review","",mixedgp_table_html(ut$main,"Training oxygen imputation",ut$note)),ocean_path(report_dir,"INDEX.md"))
  ocean_atomic_lines(c('<!doctype html><meta charset="utf-8"><title>Ocean diagnostics</title>',
    mixedgp_table_html(ut$main,"Training oxygen imputation",ut$note)),ocean_path(report_dir,"preview.html"))
  saveRDS(list(diagnostic_plots=diagnostics$plots,surface_plots=if(!is.null(surfaces))surfaces$plots else NULL,
    prediction_plots=if(!is.null(prediction))prediction$plots else NULL,u_tables=ut),ocean_path(report_dir,"internal/presentation_plots.rds"))
  saveRDS(list(source_fit_root=output_dir,report_dir=report_dir,representative_fold=representative_fold,
    no_refit=TRUE,time=as.character(Sys.time())),ocean_path(report_dir,"internal/review_run.rds"))
  progress("COMPLETE")
  invisible(list(diagnostics=diagnostics,prediction=prediction,surfaces=surfaces,report_dir=report_dir))
}

ocean_refinement_review <- function(output_dir, inputs, diagnostics,
  dest=ocean_path(output_dir,"diagnostics","refinement")) {
  if(is.null(diagnostics)) return(NULL)
  sources <- ocean_final_eiv_sources(output_dir,ocean_saved_profile(output_dir))
  stopifnot(all(sources$available))
  signature <- list(schema=3L,sources=sources,
    files=lapply(seq_len(nrow(sources)),function(i) {
      p<-ocean_path(output_dir,sources$source[i])
      list(diag=unname(tools::md5sum(ocean_path(p,"convergence_parameters.csv"))),
        fit=file.info(ocean_path(p,sources$fit_file[i]))[,c("size","mtime")])
    }),code=body(ocean_refinement_review))
  cached <- ocean_path(dest,"saved_chain_review.rds")
  if(file.exists(cached)) {
    old<-readRDS(cached)
    if(identical(old$signature,signature)) return(old$result)
  }
  p <- diagnostics$u$points; support <- list()
  for(cell in c("O20","O50")) for(C in 1:4) {
    z<-p[p$cell==cell & p$C==C,]; lower<-c(0,50,100,150)[C];upper<-c(50,100,150,180)[C]
    support[[length(support)+1L]]<-data.frame(cell,C,N=nrow(z),unique_profiles=length(unique(z$cohort_id)),
      interval_outside_band=sum(z$q025<lower | z$q975>upper),
      lower_below_zero=sum(z$q025<0),upper_above_180=sum(z$q975>180),negative_mean=sum(z$posterior_mean<0))
  }
  summary <- exceptions <- list(); tables<-list()
  finite <- function(x,fun) {x<-x[is.finite(x)];if(length(x))fun(x) else NA_real_}
  for(i in seq_len(nrow(sources))) {
    s<-sources[i,];z<-read.csv(ocean_path(output_dir,s$source,"convergence_parameters.csv"))
    tables[[i]]<-z
    summary[[i]]<-data.frame(fold=s$fold,cell=s$cell,stage=if(s$continued)ocean_continuation_label(output_dir) else "Initial",
      max_finite_Rhat=finite(z$rhat,max),min_finite_bulk_ESS=finite(z$ess_bulk,min),
      min_finite_tail_ESS=finite(z$ess_tail,min),undefined_Rhat=sum(!is.finite(z$rhat)),
      undefined_bulk=sum(!is.finite(z$ess_bulk)),undefined_tail=sum(!is.finite(z$ess_tail)),
      finite_screen_met=s$final_screen_pass)
  }
  summary <- mixedgp_bind(summary)
  selected<-which.max(summary$max_finite_Rhat)
  trace_selected<-selected
  trace<-NULL; precision<-NULL; trace_parameters<-NULL
  for(i in seq_len(nrow(sources))) {
    s<-sources[i,];z<-tables[[i]]
    bad<-z[!is.finite(z$rhat) | !is.finite(z$ess_bulk) | !is.finite(z$ess_tail),]
    if(!nrow(bad) && i!=selected && i!=trace_selected) next
    cat(sprintf("Reviewing saved chains: fold %d / %s\n",s$fold,s$cell))
    fit<-readRDS(ocean_path(output_dir,s$source,s$fit_file));by<-fit$mcmc$samples_by_chain
    for(j in seq_len(nrow(bad))) {
      name<-bad$parameter[j]
      if(grepl("^dictionary_occupancy",name)) {
        state<-as.integer(sub("^dictionary_occupancy\\[1,([0-9]+)\\]$","\\1",name))
        series<-lapply(by$J,function(x)as.numeric(as.matrix(x)[,1]==state))
      } else if(name=="J[1]") series<-lapply(by$J,function(x)as.numeric(as.matrix(x)[,1])) else series<-NULL
      reason<-"Requires separate series review";range_chain<-NA_character_;quantiles<-NA_character_
      if(!is.null(series)) {
        v<-unlist(series);constant<-length(unique(v))==1L
        reason<-if(constant && grepl("occupancy",name) && v[1]==0)
          "State never visited; rare-state mixing not established" else if(constant)
          "Constant in saved chains; not evidence of global mixing" else
          "Discrete series; inspect tail thresholds and chain frequencies"
        range_chain<-paste(sprintf("%.4f",vapply(series,mean,0)),collapse=" / ")
        quantiles<-paste(quantile(v,c(.05,.95),names=FALSE),collapse=" / ")
      }
      exceptions[[length(exceptions)+1L]]<-data.frame(fold=s$fold,cell=s$cell,parameter=name,
        undefined=paste(c("Rhat","bulk ESS","tail ESS")[!is.finite(as.numeric(unlist(bad[j,c("rhat","ess_bulk","ess_tail")])) )],collapse=", "),
        chain_means=range_chain,q05_q95=quantiles,assessment=reason)
    }
    if(i==trace_selected) {
      # Selection uses saved diagnostics, never response prediction scores.
      latent<-z[grepl("^u\\[[0-9]+\\]$",z$parameter),]
      latent$input_id<-as.integer(sub("^u\\[([0-9]+)\\]$","\\1",latent$parameter))
      latent<-latent[latent$input_id %in% fit$data$miss_idx,]
      rh<-latent[is.finite(latent$rhat),];eb<-latent[is.finite(latent$ess_bulk),]
      stopifnot(nrow(rh)>0L,nrow(eb)>0L)
      worst_rhat<-rh$input_id[which.max(rh$rhat)]
      worst_bulk<-eb$input_id[which.min(eb$ess_bulk)]
      series<-list(V=lapply(by$V,as.numeric),r=lapply(by$r,as.numeric),
        log_theta_salinity=lapply(by$logtheta,function(x)as.matrix(x)[,fit$data$theta_spec$x_index[1]]),
        tau2=lapply(by$tau,function(x)as.matrix(x)[,2]))
      series[[sprintf("u[%d] (max R-hat)",worst_rhat)]]<-lapply(by$u,function(x)as.matrix(x)[,worst_rhat])
      series[[sprintf("u[%d] (min bulk ESS)",worst_bulk)]]<-lapply(by$u,function(x)as.matrix(x)[,worst_bulk])
      precision<-eivGP:::mixedgp_summarize_diagnostic_series(series)
      trace_parameters<-names(series)
      trace<-mixedgp_bind(lapply(names(series),function(n) mixedgp_bind(lapply(seq_along(series[[n]]),function(c) {
        v<-series[[n]][[c]];idx<-seq(1L,length(v),by=max(1L,as.integer(floor(length(v)/500L))))
        displayed<-v[idx]
        window<-min(15L,length(displayed));if(window%%2L==0L)window<-window-1L
        smoothed<-as.numeric(stats::filter(displayed,rep(1/window,window),sides=2))
        data.frame(parameter=n,chain=factor(c),saved_iteration=idx,value=displayed,rolling_mean=smoothed)
      }))))
    }
    rm(fit,by);invisible(gc(FALSE))
  }
  trace$parameter<-factor(trace$parameter,levels=trace_parameters)
  trace_plot<-ggplot2::ggplot(trace,ggplot2::aes(saved_iteration,value,color=chain))+
    ggplot2::geom_line(linewidth=.25,alpha=.18)+
    ggplot2::geom_line(ggplot2::aes(y=rolling_mean),linewidth=.5,na.rm=TRUE)+
    ggplot2::facet_wrap(~parameter,ncol=3,scales="free_y")+
    ggplot2::scale_color_manual(values=c("#356A8A","#B54A63","#4E8B57","#C48A2A"))+
    ggplot2::labs(x="Post-warmup retained draw",y=NULL,color="Chain",
      title=sprintf("MCMC trace diagnostics: fold %d, %s",sources$fold[trace_selected],sources$cell[trace_selected]),
      caption="Fit selected by the largest finite Rhat across all 15 final fits; latent panels use max Rhat/min bulk ESS.\nFaint: up to 500 retained draws per chain. Heavy: centered rolling means (up to 15 points). Diagnostics use full chains.\nO0 latent inputs and cutpoints are in working-model units; traces do not assess input-model adequacy.")+
    ggplot2::theme_bw(base_size=11)+ggplot2::theme(legend.position="bottom",
      plot.caption=ggplot2::element_text(hjust=0))
  result<-list(support=mixedgp_bind(support),mcmc=summary,exceptions=mixedgp_bind(exceptions),
    precision=precision,trace=trace_plot,selected=sources[selected,],
    trace_selected=sources[trace_selected,],trace_parameters=trace_parameters)
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  for(n in c("support","mcmc","exceptions","precision")) mixedgp_atomic_csv(result[[n]],ocean_path(dest,paste0(n,".csv")))
  ggplot2::ggsave(ocean_path(dest,"selected_trace.png"),trace_plot,width=11,height=7.5,dpi=160)
  ggplot2::ggsave(ocean_path(dest,"selected_trace.pdf"),trace_plot,width=11,height=7.5)
  saveRDS(list(signature=signature,result=result),cached)
  result
}

ocean_design_review <- function(inputs) {
  rows<-overlap<-list()
  for(f in 1:5) for(cell in c("O0","O20","O50")) {
    d<-inputs$folds[[paste(f,if(cell=="O0") "O20" else cell,sep="/")]]
    rows[[length(rows)+1L]]<-data.frame(fold=f,cell,training_rows=nrow(d$train),test_rows=nrow(d$test),
      exact_training_U=if(cell=="O0") 0L else sum(d$observed_train),
      exact_mixed_test_U=if(cell=="O0") 0L else sum(d$observed_test_mixed))
    if(cell=="O0") {
      cruises<-unique(d$train$G2cruise)
      overlap[[length(overlap)+1L]]<-data.frame(Fold=f,"Test profiles"=nrow(d$test),
        "Test cruises"=length(unique(d$test$G2cruise)),
        "Test profiles sharing a training cruise"=sum(d$test$G2cruise %in% cruises),check.names=FALSE)
    }
  }
  all<-mixedgp_bind(rows)
  design<-mixedgp_bind(lapply(c("O0","O20","O50"),function(cell) {
    z<-all[all$cell==cell,];stopifnot(length(unique(z$exact_training_U))==1,length(unique(z$exact_mixed_test_U))==1)
    data.frame("Training exact U"=paste0(sub("O","",cell),"%"),"Training profiles"=z$training_rows[1],
      "Exact training U"=z$exact_training_U[1],"Test N: (x,c)"=z$test_rows[1],
      "Partial-u test: (x,u) N"=z$exact_mixed_test_U[1],
      "Partial-u test: (x,c) N"=z$test_rows[1]-z$exact_mixed_test_U[1],check.names=FALSE)
  }))
  list(design=design,overlap=mixedgp_bind(overlap))
}

ocean_u_distribution_report <- function(output_dir,inputs,representative_fold=1L,
  dest=ocean_path(output_dir,"diagnostics","u_distributions")) {
  sources<-ocean_final_eiv_sources(output_dir,ocean_saved_profile(output_dir))
  sources<-sources[sources$fold==representative_fold,]
  if(nrow(sources)!=3L || !all(sources$available))return(NULL)
  tr<-inputs$folds[[paste(representative_fold,"O50",sep="/")]]$train
  ids<-which(!tr$calib_o50);stopifnot(length(ids)==sum(!tr$calib_o50),all(!tr$calib_o20[ids]))
  files<-ocean_path(output_dir,sources$source,sources$fit_file)
  signature<-list(schema=1L,sources=sources,files=file.info(files)[,c("size","mtime")],
    ids=tr$cohort_id[ids],reference=tr$U[ids],classes=tr$C[ids],draws_per_chain=1000L)
  dir.create(dest,recursive=TRUE,showWarnings=FALSE);cache<-ocean_path(dest,"posterior_inputs.rds")
  if(file.exists(cache) && identical(readRDS(cache)$signature,signature)) {
    object<-readRDS(cache)$object
  } else {
    object<-list()
    for(cell in c("O0","O20","O50")) {
      s<-sources[sources$cell==cell,];fit<-readRDS(ocean_path(output_dir,s$source,s$fit_file))
      draw<-do.call(rbind,lapply(fit$mcmc$samples_by_chain$u,function(u) {
        idx<-unique(round(seq(1,nrow(u),length.out=min(1000L,nrow(u)))))
        as.matrix(u)[idx,ids,drop=FALSE]
      }))
      if(cell!="O0")draw<-draw*as.numeric(fit$data$U_scale)+as.numeric(fit$data$U_center)
      object[[cell]]<-draw;rm(fit,draw);gc(FALSE)
    }
    mixedgp_atomic_rds(list(signature=signature,object=object),cache)
  }
  colors<-c(O0="#777777",O20="#2365a3",O50="#c76326")
  curves<-histograms<-list()
  prepare_group<-function(C=NULL) {
    keep<-if(is.null(C))rep(TRUE,length(ids)) else tr$C[ids]==C
    ref<-tr$U[ids[keep]]
    breaks<-if(is.null(C))seq(0,180,by=15) else seq(c(0,50,100,150)[C],c(50,100,150,180)[C],length.out=6L)
    h<-hist(ref,breaks=breaks,plot=FALSE,right=FALSE,include.lowest=TRUE)
    hd<-data.frame(C=if(is.null(C))0L else C,left=head(h$breaks,-1),right=tail(h$breaks,-1),
      count=h$counts,density=h$density)
    histograms[[if(is.null(C))"overall" else paste0("class",C)]]<<-hd
    pair<-unlist(lapply(object[c("O20","O50")],function(o)as.vector(o[,keep,drop=FALSE])),use.names=FALSE)
    bw<-bw.nrd0(ref);rng<-range(c(quantile(pair,c(.001,.999)),ref,breaks))+c(-3,3)*bw
    # Same physical x/y limits and posterior KDE bandwidth across O20/O50.
    dens<-lapply(object[c("O20","O50")],function(o)density(as.vector(o[,keep,drop=FALSE]),bw=bw,from=rng[1],to=rng[2],n=512L))
    ymax<-1.04*max(c(hd$density,unlist(lapply(dens,function(d)d$y))))
    list(C=C,keep=keep,ref=ref,histogram=hd,bw=bw,range=rng,ymax=ymax)
  }
  make_panel<-function(cell,g) {
    v<-as.vector(object[[cell]][,g$keep,drop=FALSE])
    if(cell=="O0") {
      bw<-bw.nrd0(v);rng<-range(quantile(v,c(.001,.999)))+c(-3,3)*bw
    } else {bw<-g$bw;rng<-g$range}
    z<-density(v,bw=bw,from=rng[1],to=rng[2],n=512L)
    d<-data.frame(fold=representative_fold,cell,C=if(is.null(g$C))0L else g$C,u=z$x,density=z$y)
    curves[[paste(cell,if(is.null(g$C))"overall" else g$C)]]<<-d
    p<-ggplot2::ggplot()
    if(cell!="O0")p<-p+ggplot2::geom_rect(data=g$histogram,
      ggplot2::aes(xmin=left,xmax=right,ymin=0,ymax=density),fill="#c5c5c5",color="white",alpha=.65)
    p<-p+ggplot2::geom_line(data=d,ggplot2::aes(u,density),color=colors[[cell]],linewidth=.85)+
      ggplot2::labs(title=paste0(if(is.null(g$C))"" else paste0("c = ",g$C," | "),sub("O","",cell),"% calibration"),
        subtitle=if(cell=="O0")"Model units: no physical-scale overlay" else paste0("Measured reference: n = ",length(g$ref)),
        x=if(cell=="O0")"u (model units)" else "Oxygen u (micromol/kg)",y="Density")+
      ggplot2::theme_bw(base_size=10)+ggplot2::theme(legend.position="none")
    if(cell!="O0")p<-p+ggplot2::geom_rug(data=data.frame(u=g$ref),ggplot2::aes(x=u),sides="b",alpha=.35)+
      ggplot2::scale_x_continuous(limits=g$range)+ggplot2::scale_y_continuous(limits=c(0,g$ymax),expand=ggplot2::expansion(mult=c(0,.02)))
    p
  }
  overall<-prepare_group()
  pooled<-patchwork::wrap_plots(lapply(c("O0","O20","O50"),function(cell)make_panel(cell,overall)),nrow=1L)+
    patchwork::plot_annotation(title=paste0("Marginal imputed-U distribution | fold ",representative_fold),
      subtitle=sprintf("Solid: equally weighted posterior mixture. Gray histogram/rug: %d withheld measured references.",length(ids)),
      caption="Same common-hidden profiles in every panel. O0 is unanchored; measured raw u cannot be overlaid on model units.\nO20/O50 share histogram bins, posterior smoothing and axes. Marginal spread includes differences between individuals.")
  plots<-lapply(1:4,function(C) {
    g<-prepare_group(C)
    patchwork::wrap_plots(lapply(c("O0","O20","O50"),function(cell)make_panel(cell,g)),nrow=1L)
  })
  p<-patchwork::wrap_plots(plots,ncol=1L)+patchwork::plot_annotation(
    title=paste0("Class-specific imputed-U distributions | fold ",representative_fold),
    subtitle="Solid: equally weighted individual posteriors. Gray histogram/rug: the same measured references within each class.",
    caption="Same common-hidden profiles in each c. O0 has no raw-U overlay.\nO20/O50 share bins, posterior KDE bandwidths and axes within class; marginal spread is not individual uncertainty alone.")
  mixedgp_atomic_csv(do.call(rbind,curves),ocean_path(dest,"density_curves.csv"))
  mixedgp_atomic_csv(do.call(rbind,histograms),ocean_path(dest,"reference_histograms.csv"))
  mixedgp_atomic_csv(tr[ids,c("cohort_id","C","U")],ocean_path(dest,"reference_profiles.csv"))
  for(ext in c("png","pdf")) {
    ggplot2::ggsave(ocean_path(dest,paste0("class_specific.",ext)),p,width=12,height=11,dpi=160)
    ggplot2::ggsave(ocean_path(dest,paste0("overall.",ext)),pooled,width=12,height=4.5,dpi=160)
  }
  list(plot=p,pooled_plot=pooled,source=sources,common_ids=tr$cohort_id[ids],dest=dest)
}

ocean_paper_registry <- function() {
  data.frame(kind=c(rep("figure",8),rep("table",6)),
    label=c("fig:ocean-si-data","fig:ocean-si-fm","fig:ocean-si-imputation-main","fig:ocean-si-imputation-folds","fig:ocean-si-m-slices","fig:ocean-si-m-folds","fig:ocean-si-fm-folds","fig:ocean-si-mcmc-trace","tab:ocean-si-imputation","tab:ocean-si-prediction","tab:ocean-si-resolution","tab:ocean-si-prediction-full","tab:ocean-si-imputation-class","tab:ocean-si-complete-case"),
    file=c("fig_ocean_cohort_y_u_x.pdf","fig_ocean_f_m_slice.pdf","fig_ocean_u_imputation_all_training_fold1.pdf","fig_ocean_u_imputation_all_training_all_folds.pdf","fig_ocean_m_pressure_fold1.pdf","fig_ocean_m_pressure_all_folds.pdf","fig_ocean_f_m_slice_all_folds.pdf","fig_ocean500_target_traces.pdf","tab_ocean_imputation.tex","tab_ocean_prediction.tex","tab_ocean_resolution.tex","tab_ocean_prediction_full.tex","tab_ocean_imputation_class.tex","tab_ocean_complete_case.tex"),stringsAsFactors=FALSE)
}

ocean_cohort_plot <- function(cohort) {
  stopifnot(all(c("U","Y","G2salinity","log_pressure","C") %in% names(cohort)))
  cohort$oxygen_class<-factor(cohort$C,levels=1:4,labels=c("<50","50-100","100-150","150-180"))
  make<-function(x,label)ggplot2::ggplot(cohort,ggplot2::aes(x=.data[[x]],y=Y,color=oxygen_class))+
    ggplot2::geom_point(alpha=.65,size=1.1)+ggplot2::labs(x=label,y="Silicate (micromol/kg)",color="Oxygen class")+ggplot2::theme_minimal(base_size=11)
  p1<-make("U","Oxygen (micromol/kg)")+ggplot2::geom_vline(xintercept=c(50,100,150),linetype=3,color="grey55")
  patchwork::wrap_plots(p1,make("G2salinity","Salinity"),make("log_pressure","Log pressure"),guides="collect")+
    patchwork::plot_annotation(title=paste(nrow(cohort),"selected ocean profiles")) & ggplot2::theme(legend.position="bottom")
}

ocean_paper_prediction_table <- function(summary,sd_y,full=FALSE) {
  if(!nrow(summary))return(data.frame())
  z<-summary
  out<-data.frame(Method=z$method,Calibration=ifelse(z$cell=="baseline","0%",paste0(sub("O","",z$cell),"%")),
    Test=ifelse(z$prediction_mode=="mixed",paste0(sub("O","",z$cell),"% exact"),"Class only"),Folds=z$n_folds,check.names=FALSE)
  specs<-list(RMSE=c("RMSE",1),`CRPS/sY`=c("CRPS",sd_y),`Coverage (%)`=c("Coverage95",.01),`Width/sY`=c("Width95",sd_y))
  if(full)specs<-c(specs,list(`RMSE/sY`=c("RMSE",sd_y),`MAE/sY`=c("MAE",sd_y),NLPD_Z=c("NLPD",1),`Interval score/sY`=c("IntervalScore95",sd_y)))
  for(n in names(specs)) {
    metric<-specs[[n]][1];scale<-as.numeric(specs[[n]][2]);mu<-z[[paste0(metric,"_mean")]]/scale;se<-z[[paste0(metric,"_sd")]]/scale
    if(n=="NLPD_Z")mu<-mu-log(sd_y)
    out[[n]]<-ocean_mean_sd(mu,se,if(n%in%c("RMSE","Coverage (%)"))2L else 3L)
  }
  groups<-paste(ifelse(z$cell=="baseline","O0",z$cell),z$prediction_mode)
  values<-data.frame(RMSE=round(z$RMSE_mean,2),`CRPS/sY`=round(z$CRPS_mean/sd_y,3),check.names=FALSE)
  attr(out,"best")<-ocean_best_mask(values,groups,z$n_folds==5L)
  out
}

ocean_application_status <- function(output_dir,profile="paper") {
  catalogue<-ocean_task_catalogue(profile);catalogue<-catalogue[catalogue$module!="continuation",]
  rows<-lapply(seq_len(nrow(catalogue)),function(i) {
    t<-catalogue[i,];root<-if(profile=="quick")"smoke_class4" else "full5_class4"
    if(t$module=="eiv") {
      p<-ocean_path(output_dir,root,sprintf("fold%02d",t$fold),t$cell)
      return(data.frame(task=t$task_id,fold=t$fold,cell=t$cell,method="EIV-GP",status=if(file.exists(ocean_path(p,"predictive_metrics.csv")))"scored" else if(file.exists(ocean_path(p,"fit.rds")))"fit_only" else "missing"))
    }
    p<-if(t$module=="original_baselines")ocean_path(output_dir,root,sprintf("fold%02d",t$fold),"baseline_status.csv") else ocean_path(output_dir,"complete_case_baselines/fits",t$task_id,"package_status.csv")
    methods<-if(t$module=="original_baselines")c("UC-GP","LVGP","EzGP") else t$method
    z<-if(file.exists(p))read.csv(p,stringsAsFactors=FALSE) else data.frame(method=character(),status=character());j<-match(methods,z$method)
    data.frame(task=t$task_id,fold=t$fold,cell=t$cell,method=methods,status=ifelse(is.na(j),"missing",z$status[j]))
  });mixedgp_bind(rows)
}

ocean_export_application <- function(output_dir,analysis,settings,data_dir,
  diagnostics_dir=ocean_path(output_dir,"diagnostics"),dest=ocean_path(output_dir,"reports","paper")) {
  inputs<-analysis$inputs;diagnostics<-analysis$diagnostics;prediction<-analysis$prediction
  dir.create(dest,recursive=TRUE,showWarnings=FALSE);registry<-ocean_paper_registry()
  files<-list.files(output_dir,recursive=TRUE,full.names=TRUE)
  files<-files[grepl("(predictive_metrics|baseline_metrics|baseline_status|package_status|convergence_parameters|U_posterior_intervals|m_points|f_points|summary)[.]csv$|fit(_plus20k)?[.]rds$|prediction_distributions[.]rds$",files)]
  files<-files[!startsWith(files,paste0(dest,.Platform$file.sep))]
  settings_id<-settings[setdiff(names(settings),c("action","cores","run_id"))]
  signature<-list(schema=2L,settings=settings_id,files=file.info(files)[,c("size","mtime")],
    input_hashes=tools::md5sum(list.files(data_dir,"[.]csv$",full.names=TRUE)),code=eivgp_module_identity())
  cache<-file.path(dest,"export.rds")
  if(file.exists(cache)) {
    old<-readRDS(cache);m<-old$result$manifest;paths<-file.path(dest,m$file[m$available])
    if(identical(old$signature,signature) && all(file.exists(paths)) &&
       identical(unname(tools::md5sum(paths)),m$md5[m$available]))return(old$result)
  }
  # Remove only stale generated aliases; preserve all fits, raw scores and caches.
  unlink(file.path(dest,registry$file))
  ggplot2::ggsave(file.path(dest,registry$file[1]),ocean_cohort_plot(inputs$cohort),width=11,height=3.8)
  notes<-paste0("Cohort N=",inputs$design$N,"; profile=",settings$profile,". Entries are fold means (sample SD), not independent-replication uncertainty. sY=",format(inputs$response_sd,digits=8),". Missing/failed methods appear in method_status.csv; no values are imputed. Bold marks the smallest displayed RMSE/CRPS among available methods with five folds and identical calibration/test information; it does not imply significance.")
  table<-function(d,name,caption)if(!is.null(d)&&nrow(d)) {
    ocean_export_paper_table(d,name,caption,dest,notes)
    j<-match(paste0(name,".tex"),registry$file)
    if(!is.na(j)) {
      path<-file.path(dest,paste0(name,".tex"));tex<-readLines(path,warn=FALSE)
      tex<-gsub(paste0("\\label{tab:",name,"}"),paste0("\\label{",registry$label[j],"}"),tex,fixed=TRUE)
      ocean_atomic_lines(tex,path)
    }
  }
  if(!is.null(prediction)&&nrow(prediction$summary)) {
    z<-prediction$summary
    table(ocean_paper_prediction_table(z[z$cell %in% c("baseline","O0") & z$prediction_mode=="class_only",],inputs$response_sd),"tab_ocean_prediction","Prediction without exact training or test oxygen")
    table(ocean_paper_prediction_table(z[z$variant=="EIV" & z$cell %in% c("O20","O50"),],inputs$response_sd),"tab_ocean_resolution","Value of exact test oxygen within a fitted EIV-GP")
    full_table<-ocean_paper_prediction_table(z,inputs$response_sd,TRUE)
    table(full_table,"tab_ocean_prediction_full","Complete predictive summaries")
    ocean_export_prediction_panels(full_table,dest,notes)
    evidence<-ocean_prediction_evidence(output_dir,inputs,prediction,settings$profile)
    for(type in c("pooled","subgroups","information","paired"))table(ocean_evidence_table(evidence,type),paste0("tab_ocean_",type),paste("Predictive evidence:",type))
  }
  if(!is.null(analysis$cc_review)) {
    mixedgp_atomic_csv(analysis$cc_review$raw,file.path(dest,"complete_case_fold_metrics.csv"))
    cc<-analysis$cc_review$summary
    mixedgp_atomic_csv(cc,file.path(dest,"complete_case_summary.csv"))
    if(nrow(cc)) {
      display<-data.frame(Method=cc$method,Training=cc$Training,Test=cc$Test,Folds=cc$n_folds,check.names=FALSE)
      for(m in c("RMSE","CRPS","Coverage95","Width95")) {
        scale<-if(m %in% c("CRPS","Width95"))inputs$response_sd else if(m=="Coverage95").01 else 1
        label<-switch(m,RMSE="RMSE",CRPS="CRPS/sY",Coverage95="Coverage (%)",Width95="Width/sY")
        display[[label]]<-ocean_mean_sd(cc[[paste0(m,"_mean")]]/scale,cc[[paste0(m,"_sd")]]/scale,3L)
      }
      table(display,"tab_ocean_complete_case","Complete-case comparisons: matched test subsets and fold counts")
    }
  }
  review<-NULL
  if(!is.null(diagnostics)) {
    u<-ocean_paper_u_tables(diagnostics$u$summary)
    table(u$main,"tab_ocean_imputation","Oxygen imputation on common masked training profiles")
    table(u$class,"tab_ocean_imputation_class","Class-specific oxygen imputation")
    review<-ocean_refinement_review(output_dir,inputs,diagnostics,dest=file.path(diagnostics_dir,"refinement"))
    mcmc<-review$mcmc
    mcmc$iterations_per_chain<-settings$iter+vapply(seq_len(nrow(mcmc)),function(i)if(mcmc$stage[i]=="Initial")0L else ocean_continuation_iterations(output_dir,settings$profile,mcmc$fold[i],mcmc$cell[i]),integer(1))
    for(i in seq_len(nrow(mcmc))) {
      path<-ocean_path(output_dir,if(mcmc$stage[i]=="Initial")"full5_class4"else"continuation_plus20k",
        sprintf("fold%02d",mcmc$fold[i]),mcmc$cell[i],"convergence_summary.csv")
      if(file.exists(path)) {z<-utils::read.csv(path);if("iterations" %in% names(z))mcmc$iterations_per_chain[i]<-z$iterations[1]}
    }
    table(mcmc,"tab_ocean_mcmc","MCMC diagnostics for selected saved fits")
    mapping<-c(fig_ocean_f_m_slice="fm_main",fig_ocean_u_imputation_all_training_fold1="u_main",fig_ocean_u_imputation_all_training_all_folds="u_all_folds",fig_ocean_m_pressure_fold1="m_main",fig_ocean_m_pressure_all_folds="m_all_folds",fig_ocean_f_m_slice_all_folds="fm_all_folds")
    jobs<-lapply(mapping,function(n)diagnostics$plots[[n]])
    mixedgp_save_plot_jobs(jobs,dest,function(n)c(if(grepl("imputation",n))13.2 else 11,if(grepl("all_folds",n))14 else if(grepl("f_m_slice",n))4.9 else 7))
    targets<-ocean_build_target_traces(output_dir,inputs,settings)
    selected<-ocean_target_trace_selection()
    keys<-paste(targets$traces$fold,targets$traces$cell,targets$traces$target)
    if(all(paste(selected$fold,selected$cell,selected$target)%in%keys)) {
      mixedgp_save_plot_jobs(list(fig_ocean500_target_traces=ocean_target_trace_plot(targets$traces)),dest,function(n)c(7.1,6.5))
    }
    audit_dest<-ocean_path(output_dir,"diagnostics","target_traces")
    mixedgp_atomic_csv(selected,file.path(audit_dest,"compact_selection.csv"))
    plots<-lapply(unique(targets$traces$target),function(t) {
      sel<-unique(targets$traces[targets$traces$target==t,c("fold","cell","target")])
      sel<-sel[order(sel$fold,sel$cell),]
      ocean_target_trace_plot(targets$traces,sel,if(t=="m")3L else 2L)
    });names(plots)<-paste0("review_",unique(targets$traces$target))
    mixedgp_save_plot_jobs(plots,audit_dest,function(n)c(7.1,8.3))
    review$targets<-targets
    mixedgp_atomic_csv(targets$summary,file.path(dest,"target_diagnostics.csv"))
  }
  ocean_export_manuscript_tables(dest,inputs)
  status<-ocean_application_status(output_dir,settings$profile);mixedgp_atomic_csv(status,file.path(dest,"method_status.csv"))
  registry$profile<-settings$profile;registry$cohort_N<-inputs$design$N;registry$available<-file.exists(file.path(dest,registry$file))
  registry$md5<-NA_character_;registry$md5[registry$available]<-unname(tools::md5sum(file.path(dest,registry$file[registry$available])))
  registry$status<-ifelse(registry$available,"generated","pending_dependencies")
  if(settings$profile!="paper")registry$status<-paste0("workflow_check_",registry$status)
  mixedgp_atomic_csv(registry,file.path(dest,"artifact_manifest.csv"))
  state<-if(settings$profile!="paper")"WORKFLOW_CHECK" else if(!all(registry$available))"INCOMPLETE_ARTIFACTS" else if(any(!status$status %in% c("scored","success")))"COMPLETE_WITH_METHOD_FAILURES" else "COMPLETE"
  ocean_atomic_lines(c(paste("state:",state),paste("cohort:",inputs$design$N),paste("profile:",settings$profile),paste("generated:",sum(registry$available),"/",nrow(registry)),notes,"Workflow completion does not establish convergence. Nothing is copied to Overleaf.","The revised manuscript uses ocean500. Other cohorts and quick/pilot profiles are distinct analyses; exact historical values require the matching saved fits and environment."),file.path(dest,"SUMMARY.txt"))
  mixedgp_atomic_rds(list(baseline_recovery=if(file.exists(file.path(output_dir,"baseline_recovery.rds")))readRDS(file.path(output_dir,"baseline_recovery.rds")) else NULL,identity=readRDS(ocean_path(output_dir,"workflow_identity.rds")),analysis_code=eivgp_module_identity(),analysis_package=as.character(utils::packageVersion("eivGP")),settings=settings,input_hashes=tools::md5sum(list.files(data_dir,"[.]csv$",full.names=TRUE)),artifact_manifest=registry,method_status=status,time=Sys.time()),file.path(dest,"provenance.rds"))
  result<-list(manifest=registry,status=state,method_status=status,review=review,path=dest)
  mixedgp_atomic_rds(list(signature=signature,result=result),cache)
  result
}

ocean_export_prediction_panels <- function(d,dest,note) {
  model<-paste(d$Method,d$Calibration,ifelse(d$Test=="Class only","OI","MI"),paste0("[",d$Folds," folds]"))
  point<-data.frame(Model=model,d[,c("RMSE","RMSE/sY","MAE/sY")],check.names=FALSE)
  distribution<-data.frame(Model=model,d[,c("CRPS/sY","NLPD_Z","Coverage (%)","Width/sY","Interval score/sY")],check.names=FALSE)
  ocean_export_paper_table(point,"tab_ocean_prediction_full_point","Point prediction",dest,note)
  ocean_export_paper_table(distribution,"tab_ocean_prediction_full_distribution","Predictive distributions",dest,note)
  tex<-c(readLines(file.path(dest,"tab_ocean_prediction_full_point.tex")),readLines(file.path(dest,"tab_ocean_prediction_full_distribution.tex")))
  tex<-gsub("\\label{tab:tab_ocean_prediction_full_point}","\\label{tab:ocean-si-prediction-full}",tex,fixed=TRUE)
  ocean_atomic_lines(tex,file.path(dest,"tab_ocean_prediction_full.tex"))
  ocean_atomic_lines(c("<!doctype html><meta charset='utf-8'>",mixedgp_table_html(point,"Point prediction",note),
    mixedgp_table_html(distribution,"Predictive distributions",note)),file.path(dest,"tab_ocean_prediction_full.html"))
}

ocean_target_trace_settings <- function() list(schema=1L,states=256L,nodes=256L,
  check_states=24L,salinity=34.25,pressure=650,oxygen=75,class=2L)

ocean_target_trace_fit <- function(fit,train,fold,cell,burn=0L,
  spec=ocean_target_trace_settings(),coordinate_diagnostics=NULL) {
  by<-fit$mcmc$samples_by_chain;nc<-length(by$u)
  n<-vapply(by$u,nrow,integer(1));stopifnot(nc>=2L,length(unique(n))==1L)
  n<-n[1];idx<-unique(as.integer(round(seq(1,n,length.out=min(n,spec$states)))))
  ids<-unlist(lapply(0:(nc-1L),function(c)c*n+idx))
  info<-fit$mcmc$mcmc_draw_info
  if(!is.null(info))stopifnot(identical(as.integer(info$chain[ids]),rep(seq_len(nc),each=length(idx))))
  xraw<-matrix(c(spec$salinity,log(spec$pressure)),1,
    dimnames=list(NULL,c('G2salinity','log_pressure')))
  xs<-sweep(sweep(xraw,2,fit$data$x_center,'-'),2,fit$data$x_scale,'/')
  state_mean<-function(j,nodes=spec$nodes) {
    tau<-c(-Inf,fit$mcmc$samples_tau[j,],Inf);c<-spec$class
    lo<-stats::pnorm(tau[c]);hi<-stats::pnorm(tau[c+1L])
    if(!is.finite(hi-lo)||hi<=lo)stop('Degenerate class probability in target trace')
    u<-matrix(stats::qnorm(lo+(seq_len(nodes)-.5)/nodes*(hi-lo)),ncol=1)
    p<-gp_integrated_mean_state_1d(x_train=fit$data$x,u_train=fit$mcmc$samples_u[j,],
      y_train=fit$data$y,x_star=xs,u_mc=u,logtheta=fit$mcmc$samples_logtheta[j,],
      sigma2_eps=fit$mcmc$samples_sigma2[j],kernel='se',matern_nu=NULL,return_cov=FALSE)
    as.numeric(fit$data$y_center+fit$data$y_scale*p$mean)
  }
  # Chunked inner jobs borrow the same pool as outer fold jobs.
  m<-unlist(mixedgp_stage_lapply(split(ids,ceiling(seq_along(ids)/32L)),
    function(js)vapply(js,state_mean,0)),use.names=FALSE)
  check<-unique(round(seq(1,length(ids),length.out=min(length(ids),spec$check_states))))
  alt<-vapply(ids[check],function(j)state_mean(j,2L*spec$nodes),0)
  series<-list(m=matrix(m,nrow=length(idx),ncol=nc));full<-list();individual<-NA_integer_
  if(cell!='O0') {
    f<-predict_eivgp(fit,new_X=xraw,new_U=matrix(spec$oxygen,1),new_U_scale='raw',
      target='surface',draw_ids=ids,include_process_uncertainty=FALSE)
    series$f<-matrix(as.numeric(f),nrow=length(idx),ncol=nc)
    common<-which(!ocean_diag_flag(train$calib_o50)&train$C==spec$class)
    stopifnot(length(common)>0L,all(common %in% fit$data$miss_idx))
    full$Umean<-sapply(by$u,function(x)rowMeans(x[,common,drop=FALSE])*fit$data$U_scale+fit$data$U_center)
    dd<-coordinate_diagnostics
    if(!is.null(dd)&&all(c('parameter','rhat')%in%names(dd))) {
      dd<-dd[grepl('^u\\[[0-9]+\\]$',dd$parameter),]
      dd$id<-as.integer(sub('^u\\[([0-9]+)\\]$','\\1',dd$parameter))
      dd<-dd[dd$id%in%common & is.finite(dd$rhat),]
    } else dd<-data.frame()
    if(!nrow(dd)) {
      rr<-vapply(common,function(j)posterior::rhat(sapply(by$u,function(x)x[,j])),0)
      dd<-data.frame(id=common,rhat=rr);dd<-dd[is.finite(dd$rhat),]
    }
    # Deterministic fallback for constant coordinates; record the chosen row.
    individual<-if(nrow(dd))dd$id[which.max(dd$rhat)] else common[1L]
    full$Uindividual<-sapply(by$u,function(x)x[,individual]*fit$data$U_scale+fit$data$U_center)
    for(nm in names(full))series[[nm]]<-full[[nm]][idx,,drop=FALSE]
  }
  summary<-do.call(rbind,lapply(names(series),function(nm) {
    mat<-if(nm%in%names(full))full[[nm]] else series[[nm]]
    data.frame(fold=fold,cell=cell,target=nm,rhat=posterior::rhat(mat),
      ess=posterior::ess_bulk(mat),mcse_mean=posterior::mcse_mean(mat),
      chain_mean_range=diff(range(colMeans(mat))),pooled_sd=stats::sd(as.numeric(mat)),
      diagnostic_draws=nrow(mat),chains=nc)
  }))
  traces<-do.call(rbind,lapply(names(series),function(nm)do.call(rbind,lapply(seq_len(nc),
    function(c)data.frame(fold=fold,cell=cell,target=nm,chain=c,
      iteration=burn+idx,value=series[[nm]][,c])))))
  list(traces=traces,summary=summary,individual=individual,ids=ids,
    integration_check=data.frame(fold=fold,cell=cell,max_abs_delta=max(abs(alt-m[check])),
      mean_abs_delta=mean(abs(alt-m[check]))))
}

ocean_build_target_traces <- function(output_dir,inputs,settings,
  dest=ocean_path(output_dir,'diagnostics','target_traces'),spec=ocean_target_trace_settings()) {
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  sources<-ocean_final_eiv_sources(output_dir,settings$profile)
  sources<-sources[sources$available,,drop=FALSE]
  if(!nrow(sources))stop('No saved fits for scientific-target traces')
  paths<-mixedgp_stage_lapply(seq_len(nrow(sources)),function(i) {
    s<-sources[i,];path<-ocean_path(output_dir,s$source,s$fit_file)
    dp<-ocean_path(output_dir,s$source,'convergence_parameters.csv')
    train<-inputs$folds[[paste(s$fold,'O50',sep='/')]]$train
    signature<-list(spec=spec,fit=unname(tools::md5sum(path)),burn=settings$burn,
      coordinate_diagnostics=if(file.exists(dp))unname(tools::md5sum(dp)) else NULL,
      train=digest::digest(train),code=digest::digest(body(ocean_target_trace_fit)),
      package_code=eivgp_module_identity())
    cache<-file.path(dest,sprintf('fold%02d_%s.rds',s$fold,s$cell))
    if(file.exists(cache)) {
      old<-readRDS(cache)
      if(identical(old$signature,signature))return(cache)
    }
    fit<-readRDS(path)
    z<-ocean_target_trace_fit(fit,train,s$fold,s$cell,settings$burn,spec,
      if(file.exists(dp))utils::read.csv(dp) else NULL)
    mixedgp_atomic_rds(list(signature=signature,source=s,result=z),cache);cache
  },workers=min(4L,mixedgp_budget_cores()))
  z<-lapply(paths,function(p)readRDS(p)$result)
  out<-list(traces=mixedgp_bind(lapply(z,`[[`,'traces')),summary=mixedgp_bind(lapply(z,`[[`,'summary')),
    integration_check=mixedgp_bind(lapply(z,`[[`,'integration_check')),
    coordinates=data.frame(fold=sources$fold,cell=sources$cell,
      selected_row=vapply(z,`[[`,integer(1),'individual')))
  mixedgp_atomic_csv(out$summary,file.path(dest,'target_diagnostics.csv'))
  mixedgp_atomic_csv(out$integration_check,file.path(dest,'integration_check.csv'))
  mixedgp_atomic_csv(out$coordinates,file.path(dest,'selected_coordinates.csv'))
  mixedgp_atomic_rds(out$traces,file.path(dest,'traces.rds'))
  out
}

ocean_target_trace_selection <- function() {
  # Frozen manuscript display; do not reselect cases from each rerun's R-hat.
  data.frame(fold=c(5L,1L,2L,1L,2L,1L,1L,1L,2L),
    cell=c('O0','O20','O50','O20','O20','O50','O20','O20','O20'),
    target=c('m','m','m','f','f','f','Umean','Uindividual','Uindividual'))
}

ocean_target_trace_plot <- function(traces,selection=ocean_target_trace_selection(),ncol=3L) {
  parts<-lapply(seq_len(nrow(selection)),function(i) {
    a<-selection[i,];z<-traces[traces$fold==a$fold & traces$cell==a$cell & traces$target==a$target,,drop=FALSE]
    if(!nrow(z))return(NULL)
    label<-switch(a$target,m='m(x, class 2)',f='f(x, oxygen = 75)',
      Umean='Mean masked oxygen',Uindividual='Individual masked oxygen')
    z$panel<-paste0('(',letters[i],') ',label,'\n',sub('O','',a$cell),'% calibration, fold ',a$fold)
    z
  });z<-mixedgp_bind(parts)
  if(!nrow(z))stop('No traces for the requested figure panels')
  z$panel<-factor(z$panel,levels=unique(z$panel));z$chain<-factor(z$chain)
  z<-mixedgp_bind(lapply(split(z,interaction(z$panel,z$chain,drop=TRUE)),function(q) {
    q<-q[order(q$iteration),];q$smooth<-if(nrow(q)>=9L)as.numeric(stats::filter(q$value,rep(1/9,9),sides=2)) else q$value;q
  }))
  paper_axis<-identical(as.numeric(range(z$iteration)),c(5001,50000))
  ggplot2::ggplot(z,ggplot2::aes(iteration/1000,value,color=chain,group=chain))+
    ggplot2::geom_line(linewidth=.19,alpha=.28)+
    ggplot2::geom_line(ggplot2::aes(y=smooth),linewidth=.4,na.rm=TRUE)+
    ggplot2::facet_wrap(~panel,ncol=ncol,scales='free_y')+
    ggplot2::scale_color_manual(values=c('#0072B2','#D55E00','#009E73','#CC79A7'))+
    ggplot2::scale_x_continuous(breaks=if(paper_axis)c(10,30,50) else ggplot2::waiver())+
    ggplot2::labs(x=if(paper_axis)'Iteration (thousands; burn-in ends at 5)' else 'Iteration (thousands; retained states)',y='micromol/kg',color='Chain')+
    ggplot2::theme_bw(base_size=9)+ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),
      panel.grid.major=ggplot2::element_line(linewidth=.15,color='grey90'),legend.position='bottom',
      strip.background=ggplot2::element_rect(fill='grey96'),strip.text=ggplot2::element_text(size=8.5),
      panel.spacing=grid::unit(.6,'lines'),legend.margin=ggplot2::margin(0,0,0,0),
      plot.margin=ggplot2::margin(3,4,3,3))
}

ocean_export_manuscript_tables <- function(dest,inputs) {
  read<-function(name) {
    p<-file.path(dest,paste0(name,'.csv'))
    if(file.exists(p))utils::read.csv(p,check.names=FALSE,stringsAsFactors=FALSE) else NULL
  }
  tex<-function(x)gsub('%','\\%',x,fixed=TRUE)
  mark<-function(d,keys,groups=rep('all',nrow(d)),eligible=rep(TRUE,nrow(d))) {
    for(g in unique(groups))for(k in keys) {
      ix<-which(groups==g & eligible);if(!length(ix))next
      v<-suppressWarnings(as.numeric(sub(' .*','',d[[k]][ix])))
      if(!any(is.finite(v)))next
      best<-ix[is.finite(v)&v==min(v,na.rm=TRUE)]
      d[[k]][best]<-paste0('\\textbf{',sub(' .*','',d[[k]][best]),'}',sub('^[^ ]+','',d[[k]][best]))
    };d
  }
  tabular<-function(d,headers,groups=NULL) {
    out<-c(paste0('\\begin{tabular}{',paste0(c('l',rep('r',ncol(d)-1)),collapse=''),'}'),
      '\\toprule',paste0(paste(headers,collapse=' & '),' \\\\'),'\\midrule')
    for(i in seq_len(nrow(d))) {
      if(!is.null(groups)&&(i==1L||groups[i]!=groups[i-1L]))out<-c(out,
        paste0('\\multicolumn{',ncol(d),'}{l}{\\textit{',groups[i],'}}\\\\[2pt]'))
      out<-c(out,paste0(paste(tex(as.character(d[i,])),collapse=' & '),' \\\\'))
    };c(out,'\\bottomrule','\\end{tabular}')
  }
  note<-paste0('Entries are five-fold means (sample SDs), not independent-replication uncertainty. ',
    'RMSE is in $\\mu\\mathrm{mol}\\,\\mathrm{kg}^{-1}$; $s_Y=',sprintf('%.4f',inputs$response_sd),
    '$. Coverage and width refer to 95\\% intervals. Bold marks the lowest displayed loss within each comparison group, not significance. ',
    'Missing or failed competitors and fold counts are retained in the CSV exports and method\\_status.csv.')
  write<-function(name,label,caption,body,extra='')ocean_atomic_lines(c(
    '\\begin{table}[!htbp]','\\centering\\small',paste0('\\caption{',caption,'}'),
    paste0('\\label{',label,'}'),'\\setlength{\\tabcolsep}{3pt}',body,
    '\\par\\smallskip','\\begin{minipage}{\\textwidth}\\footnotesize',note,extra,
    '\\end{minipage}','\\end{table}'),file.path(dest,paste0(name,'.tex')))
  order_methods<-function(d)d[order(match(d$Method,c('UC-GP','LVGP','EzGP','EIV-GP'))),,drop=FALSE]
  d<-read('tab_ocean_prediction')
  if(!is.null(d)&&nrow(d)) {
    d<-order_methods(d);d<-mark(d,c('RMSE','CRPS/sY'),eligible=d$Folds==5L)
    write('tab_ocean_prediction','tab:ocean-si-prediction','Silicate prediction without exact training or test oxygen.',
      tabular(d[,c('Method','RMSE','CRPS/sY','Coverage (%)','Width/sY')],c('Method','RMSE','CRPS/$s_Y$','Cov. (\\%)','Width/$s_Y$')))
  }
  d<-read('tab_ocean_resolution')
  if(!is.null(d)&&nrow(d)) {
    d<-d[order(d$Calibration,d$Test!="Class only"),];d<-mark(d,c('RMSE','CRPS/sY'),d$Calibration,d$Folds==5L)
    write('tab_ocean_resolution','tab:ocean-si-resolution','Value of exact test oxygen within each fitted EIV--GP model.',
      tabular(d[,c('Calibration','Test','RMSE','CRPS/sY','Coverage (%)','Width/sY')],
        c('Training exact','Test oxygen','RMSE','CRPS/$s_Y$','Cov. (\\%)','Width/$s_Y$')))
  }
  for(name in c('tab_ocean_imputation','tab_ocean_imputation_class')) {
    d<-read(name);if(is.null(d)||!nrow(d))next
    class<-name=='tab_ocean_imputation_class'
    d<-mark(d,c('RMSE','MAE'),if(class)d$Class else rep('all',nrow(d)))
    headers<-c('Calibration',if(class)'Class','Bias','RMSE','MAE','Cov. (\\%)','Width')
    # Oxygen, unlike response tables, has physical units in all loss columns.
    old_note<-note;note<-paste0('Entries are five-fold means (sample SDs). Bias, RMSE, MAE and width are in $\\mu\\mathrm{mol}\\,\\mathrm{kg}^{-1}$. ',
      'Bias is posterior mean minus reference oxygen; coverage and width refer to equal-tailed 95\\% individual-oxygen intervals. ',
      'Bold marks the lower displayed RMSE and MAE within each comparison.')
    write(name,if(class)'tab:ocean-si-imputation-class' else 'tab:ocean-si-imputation',
      paste0(if(class)'Class-specific oxygen imputation' else 'Oxygen imputation',
        ' on the common ',inputs$design$N*.4,' masked training profiles per fold.'),tabular(d,headers))
    note<-old_note
  }
  d<-read('tab_ocean_prediction_full')
  if(!is.null(d)&&nrow(d)) {
    cal<-as.numeric(sub('%','',d$Calibration));d<-d[order(cal,match(d$Method,c('UC-GP','LVGP','EzGP','EIV-GP')),d$Test!='Class only'),]
    groups<-ifelse(d$Calibration=='0%','A. No exact training oxygen',paste0(ifelse(d$Calibration=='20%','B. ','C. '),tex(d$Calibration),' exact training oxygen'))
    d<-mark(d,c('RMSE','RMSE/sY','MAE/sY','CRPS/sY','NLPD_Z','Interval score/sY'),groups,d$Folds==5L)
    d$Model<-paste(d$Method,ifelse(d$Method=='EIV-GP',d$Calibration,''),ifelse(d$Test=='Class only','','mixed'))
    b<-c(tabular(d[,c('Model','RMSE','RMSE/sY','MAE/sY')],c('Method','RMSE','RMSE/$s_Y$','MAE/$s_Y$'),groups),
      '\\medskip',tabular(d[,c('Model','CRPS/sY','NLPD_Z','Coverage (%)','Width/sY','Interval score/sY')],
        c('Method','CRPS/$s_Y$','NLPD$_Z$','Cov. (\\%)','Width/$s_Y$','Score/$s_Y$'),groups))
    write('tab_ocean_prediction_full','tab:ocean-si-prediction-full','Complete predictive performance summaries.',b,
      'Unmarked rows use class-only test inputs; mixed rows reveal the indicated calibration fraction at test. NLPD is standardized by subtracting $\\log s_Y$.')
  }
  d<-read('tab_ocean_complete_case')
  if(!is.null(d)&&nrow(d)) {
    # Keep the two common exact-test subsets; other protocols remain in CSV/HTML.
    d<-d[grepl('^[0-9]+ \\(x,u(,c)?\\)$',d$Test),,drop=FALSE]
    if(nrow(d)) {
      cal<-sub('%.*','',d$Training);method<-sub('(CC)','',d$Method,fixed=TRUE)
      ix<-order(as.numeric(cal),match(method,c('UC-GP','LVGP','EzGP','EIV-GP')));d<-d[ix,];cal<-cal[ix]
      groups<-paste0(tex(paste0(cal,'%')),' calibration; ',sub(' .*','',d$Test),' exact-oxygen test profiles per fold')
      d<-mark(d,c('RMSE','CRPS/sY'),groups,d$Folds==5L)
      d[['Training n']]<-as.integer(sub('.*[;:] ([0-9]+) rows.*','\\1',d$Training))
      write('tab_ocean_complete_case','tab:ocean-si-complete-case','Prediction on the common exact-oxygen test subset.',
        tabular(d[,c('Method','Training n','RMSE','CRPS/sY','Coverage (%)','Width/sY')],
          c('Method','Training $n$','RMSE','CRPS/$s_Y$','Cov. (\\%)','Width/$s_Y$'),groups),
        'Complete-case competitors use $(x,u,c)$ and only calibrated training records; EIV--GP retains all training responses and predicts from $(x,u)$.')
    }
  }
  invisible(NULL)
}

# Dataset preparation is separate from fitting and report rendering.
ocean_dataset_dir <- function(repository_dir,dataset="ocean500",data_dir="") {
  dataset<-match.arg(dataset,c("ocean200","ocean500","ocean1000"))
  if(nzchar(data_dir))return(ocean_resolve_path(data_dir,repository_dir))
  if(dataset=="ocean1000") {
    path<-file.path(repository_dir,"reproduction/ocean/prepared/ocean1000")
    if(!file.exists(file.path(path,"fold01_train.csv")))
      stop("Prepare ocean1000 first: Rscript tools/ocean-prepare-folds.R ocean1000; see README_ocean.md.")
    return(path)
  }
  file.path(repository_dir,"data/ocean_silicate",sub("ocean","n",dataset))
}
ocean_default_output <- function(dataset,profile,repository_dir=getwd()) {
  file.path(repository_dir,"results",if(profile=="paper")dataset else paste(dataset,profile,sep="_"))
}
ocean_output_dir <- function(repository_dir,dataset,profile,output_dir="") {
  if(nzchar(output_dir))return(ocean_resolve_path(output_dir,repository_dir))
  base<-ocean_default_output(dataset,profile,repository_dir)
  marker<-paste0(base,".latest-run.rds")
  if(!file.exists(marker))return(base)
  latest<-readRDS(marker)
  if(!is.character(latest)||length(latest)!=1L||!dir.exists(latest))
    stop("Latest run is unavailable. Specify an existing run with --output=...: ",marker)
  latest
}

# Only initial fitting may select a fresh root. Old fits remain untouched.
ocean_choose_fit_output <- function(output_dir,identity,new_base=output_dir) {
  idfile<-ocean_path(output_dir,"workflow_identity.rds")
  if(file.exists(idfile)) {
    if(identical(readRDS(idfile),identity))return(output_dir)
  } else {
    existing<-if(dir.exists(output_dir))list.files(output_dir,all.files=TRUE,no..=TRUE) else character()
    if(!length(setdiff(existing,c("jobs",".workflow.lock","_meta","logs"))))return(output_dir)
  }
  stem<-paste0(new_base,"_run_",format(Sys.time(),"%Y%m%d_%H%M%S"))
  candidate<-stem;i<-1L
  repeat {
    if(dir.create(candidate,showWarnings=FALSE))break
    if(!file.exists(candidate))stop("Cannot create output directory: ",candidate)
    i<-i+1L;candidate<-paste0(stem,"_",i)
  }
  cat("Existing results preserved: ",output_dir,"\nNew run directory: ",candidate,"\n",
    "For subsequent commands targeting this run, use --output=",shQuote(candidate),"\n",sep="")
  candidate
}

# CLI options select the existing workflow; they do not alter model fitting.
ocean_cli_params <- function(args=character()) {
  if(length(args)==1L && !startsWith(args,"--"))return(readRDS(args))
  aliases<-c("data-root"="data_dir","output"="output_dir")
  allowed<-c("action","profile","dataset","cores","data_dir","output_dir",
    "folds","cells","modules","task_ids","run_id","import_dirs","cc_output_dir",
    "diagnostics_dir","representative_fold","continuation_iter","ess_min")
  p<-list()
  for(arg in args) {
    if(!grepl("^--[^=]+=.+$",arg))stop("Use --name=value: ",arg)
    key<-sub("^--([^=]+)=.*$","\\1",arg)
    key<-if(key%in%names(aliases))unname(aliases[[key]]) else gsub("-","_",key)
    if(key%in%c("chains","iter","burn"))stop("Sampling settings are fixed by --profile; remove --",key,". Use --profile=pilot for a short check.")
    if(!key%in%allowed)stop("Unknown option: ",key)
    if(key%in%names(p))stop("Repeated option: ",key)
    p[[key]]<-sub("^--[^=]+=","",arg)
  }
  for(key in intersect(names(p),c("cores","chains","iter","burn","representative_fold","continuation_iter"))) {
    value<-suppressWarnings(as.numeric(p[[key]]))
    if(length(value)!=1L || !is.finite(value) || value<(if(key=="burn")0 else 1) || value>.Machine$integer.max || value!=floor(value))
      stop(key,if(key=="burn")" must be a non-negative integer" else " must be a positive integer")
    p[[key]]<-as.integer(value)
  }
  for(key in intersect(names(p),c("cells","modules","task_ids","import_dirs")))
    p[[key]]<-strsplit(p[[key]],",",fixed=TRUE)[[1]]
  if(!is.null(p$folds)) {
    parts<-strsplit(p$folds,",",fixed=TRUE)[[1]]
    if(any(!grepl("^[1-5](:[1-5])?$",parts)))stop("Use --folds=1:5 or --folds=1,3")
    p$folds<-unique(unlist(lapply(parts,function(part) {
      z<-as.integer(strsplit(part,":",fixed=TRUE)[[1]])
      if(length(z)==2L)seq.int(z[1],z[2]) else z
    })))
  }
  if(length(args) && is.null(p$action))p$action <- "fit"
  p
}

# Existing fits own their original sampling settings; extensions only add steps.
ocean_restore_fit_settings <- function(params,output_dir) {
  path<-ocean_path(output_dir,"workflow_identity.rds")
  if(!file.exists(path))return(params)
  # New fits always use their profile, never inherit old numerical settings.
  if(params$action%in%c("fit","run","plan"))return(params)
  saved<-readRDS(path)$numerics
  if(params$action=="continue" && identical(params$profile,"paper") &&
     !identical(as.integer(unlist(saved[c("chains","iter","burn")])),c(4L,50000L,5000L)))
    stop("This saved run does not use the fixed paper sampling settings; use analyze/report to inspect it or start a new paper fit.")
  for(key in c("chains","iter","burn","ess_min"))params[[key]]<-saved[[key]]
  params
}

ocean_check_report_dependencies <- function() {
  if(!requireNamespace("rmarkdown",quietly=TRUE))stop("HTML reporting requires the R package rmarkdown.")
  if(!rmarkdown::pandoc_available())stop("HTML reporting requires Pandoc on PATH (or RSTUDIO_PANDOC).")
  invisible(TRUE)
}

ocean_render_report <- function(repository_dir,output_dir,settings,dataset,data_dir,
  cc_output_dir="",diagnostics_dir="",representative_fold=1L,analysis=NULL) {
  render_env<-new.env(parent=globalenv())
  if(!is.null(analysis))render_env$.ocean_precomputed_analysis<-list(
    output_dir=normalizePath(output_dir),data_dir=normalizePath(data_dir),
    profile=settings$profile,folds=settings$folds,cells=settings$cells,analysis=analysis)
  rmarkdown::render(file.path(repository_dir,"paper/applications/ocean.Rmd"),
    output_file="ocean.html",output_dir=output_dir,intermediates_dir=output_dir,
    params=list(dataset=dataset,profile=settings$profile,action="report",folds=settings$folds,
      cells=settings$cells,modules=settings$modules,task_ids=settings$task_ids,
      cores=settings$cores,run_id=paste0(settings$run_id,"_report"),data_dir=data_dir,output_dir=output_dir,
      cc_output_dir=cc_output_dir,diagnostics_dir=diagnostics_dir,representative_fold=representative_fold),
    envir=render_env,quiet=TRUE)
}

# Orchestration is separate from the package's fit/analyze implementations.
# A failed later stage leaves verified checkpoints available for resume.
ocean_cli_workflow <- function(settings,execute,render) {
  result<-execute(settings,TRUE)
  if(settings$action%in%c("fit","continue")) {
    analyzed<-settings;analyzed$action<-"analyze"
    analyzed$run_id<-paste0(settings$run_id,"_analyze")
    analysis<-execute(analyzed,FALSE)
    result<-list(fitting=result,analysis=analysis)
  }
  if(settings$action!="plan") {
    completed<-if(settings$action %in% c("fit","continue"))result$analysis else result
    result$html<-if("analysis" %in% names(formals(render)))render(analysis=completed) else render()
  }
  result
}

ocean_run_cli <- function(repository_dir,params=list()) {
  if(any(c("chains","iter","burn")%in%names(params)))
    stop("Sampling settings are fixed by --profile; remove --chains, --iter and --burn.")
  if("diagnostics_workers"%in%names(params))stop("diagnostics_workers was removed; use cores for the shared workflow budget.")
  api<-new.env(parent=environment())
  ocean_load_api(repository_dir,api)
  getp<-function(key,default)if(is.null(params[[key]]))default else params[[key]]
  dataset<-match.arg(getp("dataset","ocean500"),c("ocean200","ocean500","ocean1000"))
  params$run_id<-getp("run_id",paste(getp("profile","paper"),
    if(getp("action","plan")=="run")"fit" else getp("action","plan"),sep="_"))
  profile<-match.arg(getp("profile","paper"),c("paper","pilot","quick"))
  environment(ocean_dataset_dir)<-api
  environment(ocean_output_dir)<-api
  data_dir<-ocean_dataset_dir(repository_dir,dataset,getp("data_dir",""))
  output_dir<-ocean_output_dir(repository_dir,dataset,profile,getp("output_dir",""))
  params$profile<-profile
  params$action<-getp("action","plan")
  params<-ocean_restore_fit_settings(params,output_dir)
  settings<-do.call(api$ocean_settings,
    params[intersect(names(params),names(formals(api$ocean_settings)))])
  resolve<-function(paths)vapply(paths,function(path)
    if(nzchar(path))api$ocean_resolve_path(path,repository_dir) else "",
    character(1),USE.NAMES=FALSE)
  run_identity<-NULL
  if(settings$action=="fit") {
    ocean_check_report_dependencies()
    run_identity<-api$ocean_identity(repository_dir,data_dir,settings)
    output_dir<-ocean_choose_fit_output(output_dir,run_identity,
      if(nzchar(getp("output_dir","")))output_dir else ocean_default_output(dataset,profile,repository_dir))
  }
  cat("Output directory:",output_dir,"\n")
  execute<-function(stage_settings,import=TRUE) {
    if(stage_settings$action=="fit" && !nzchar(getp("output_dir","")))
      on.exit({
        identity_path<-ocean_path(output_dir,"workflow_identity.rds")
        if(file.exists(identity_path) && identical(readRDS(identity_path),run_identity))
          api$mixedgp_atomic_rds(output_dir,paste0(ocean_default_output(dataset,profile,repository_dir),".latest-run.rds"))
      },add=TRUE)
    ocean_analysis(repository_dir,data_dir,output_dir,stage_settings,
    import_dirs=if(import)resolve(getp("import_dirs",character())) else character(),
    cc_output_dir=resolve(getp("cc_output_dir","")),
    diagnostics_dir=resolve(getp("diagnostics_dir","")),
    representative_fold=getp("representative_fold",1L))
  }
  render<-function(analysis=NULL)ocean_render_report(repository_dir,output_dir,settings,dataset,data_dir,
    resolve(getp("cc_output_dir","")),resolve(getp("diagnostics_dir","")),
    getp("representative_fold",1L),analysis=analysis)
  if(settings$action=="plan") {
    result<-execute(settings)
    print(result$tasks,row.names=FALSE)
    print(result$stages,row.names=FALSE)
    print(api$mixedgp_core_plan(settings$cores,sum(result$tasks$module=="eiv"),settings$chains))
    cat("Total shared core budget:",settings$cores,"; up to",settings$chains,
      "chain workers per EIV fit; prediction, diagnostics and plots share the same budget; HTML assembly is serial.\n")
  } else {
    if(settings$action=="fit")api$ocean_storage_init(output_dir)
    # Own the result root through all stages, including reporting and HTML.
    if(dir.exists(output_dir)) {
      application_lock <- file.path(output_dir,".application.lock")
      if(!dir.create(application_lock,showWarnings=FALSE))
        stop("Ocean application already running: ",application_lock,". Inspect its PID before removing a stale lock.")
      writeLines(as.character(Sys.getpid()),file.path(application_lock,"pid"))
      on.exit(unlink(application_lock,recursive=TRUE),add=TRUE)
    }
    event_dir <- file.path(api$ocean_job_dir(output_dir,settings$run_id),"core-events")
    old_options <- options(eivGP.core_events=event_dir)
    on.exit(options(old_options),add=TRUE)
    release_budget <- getFromNamespace("mixedgp_budget_enter","eivGP")(settings$cores)
    on.exit(release_budget(),add=TRUE,after=FALSE)
    # Fail early if HTML cannot be built, before starting an expensive fit.
    ocean_check_report_dependencies()
    result<-api$ocean_logged(file.path(api$ocean_job_dir(output_dir,settings$run_id),"execution.log"),
      ocean_cli_workflow(settings,execute,render))
    cat("Saved results:",output_dir,"\nHTML:",result$html,"\n")
  }
  invisible(result)
}


# Shared CLI orchestration only. All fitting, scoring and plotting use eivGP.
ocean_method_params <- function(args,group) {
  special<-c('reuse-only','plan');flags<-setNames(as.list(rep(FALSE,2)),special)
  for(key in special) {
    ix<-grep(paste0('^--',key,'='),args)
    if(length(ix)>1L)stop('Repeated option: ',key)
    if(length(ix)) {
      value<-sub('^[^=]+=','',args[ix])
      if(!value%in%c('true','false'))stop('--',key,' must be true or false')
      flags[[key]]<-value=='true';args<-args[-ix]
    }
  }
  if(any(grepl('^--(action|modules|import-dirs)=',args)))
    stop('Choose the method script; use --reuse-only=true for saved results or --plan=true to inspect tasks.')
  p<-ocean_cli_params(args)
  p$action<-if(flags$plan)'plan' else if(group=='report'||flags[['reuse-only']])'report' else 'fit'
  p$modules<-switch(group,eivGP=c('eiv','prediction','diagnostics'),
    competitors=c('original_baselines','complete_case_baselines','prediction'),
    report=c('eiv','original_baselines','complete_case_baselines','prediction','diagnostics','surface_comparison'))
  if(is.null(p$run_id))p$run_id<-paste0('ocean_',group,if(p$action=='report')'_saved' else '_run')
  p
}

# Reporting historical fits is permitted; changing their fitting identity is not.
ocean_saved_method_check <- function(output_dir,data_dir,settings) {
  F<-function(n)get(n,envir=environment(ocean_method_run),inherits=TRUE)
  selected<-F('ocean_select_tasks')(settings)
  required<-unlist(lapply(seq_len(nrow(selected)),function(i)
    F('ocean_path')(output_dir,F('ocean_task_files')(selected[i,,drop=FALSE],settings$profile,'analyze'))))
  missing<-required[!file.exists(required)]
  if(length(missing))stop('Saved-only workflow requires completed fits and predictions; no fitting started. Missing: ',paste(head(missing,8),collapse=', '))
  inputs<-list.files(data_dir,'[.]csv$',full.names=TRUE)
  saved<-F('ocean_path')(output_dir,'data',basename(inputs))
  if(!all(file.exists(saved))||!identical(unname(tools::md5sum(inputs)),unname(tools::md5sum(saved))))
    stop('Prepared inputs differ from the saved run; no fitting started.')
  invisible(selected)
}

ocean_method_run <- function(repository_dir,params,group) {
  F<-function(n)get(n,envir=environment(ocean_method_run),inherits=TRUE)
  api<-new.env(parent=environment());ocean_load_api(repository_dir,api)
  getp<-function(n,default)if(is.null(params[[n]]))default else params[[n]]
  profile<-getp('profile','paper');dataset<-getp('dataset','ocean500')
  environment(ocean_dataset_dir)<-api;environment(ocean_output_dir)<-api
  data_dir<-ocean_dataset_dir(repository_dir,dataset,getp('data_dir',''))
  out<-ocean_output_dir(repository_dir,dataset,profile,getp('output_dir',''))
  params$profile<-profile;params<-ocean_restore_fit_settings(params,out)
  settings<-do.call(api$ocean_settings,params[intersect(names(params),names(formals(api$ocean_settings)))])
  if(settings$action=='plan') {
    print(api$ocean_select_tasks(settings),row.names=FALSE)
    cat('Group:',group,'; shared core budget:',settings$cores,'\n');return(invisible(settings))
  }
  if(settings$action=='report')ocean_saved_method_check(out,data_dir,settings)
  if(group=='report')return(ocean_run_cli(repository_dir,params[setdiff(names(params),c('chains','iter','burn'))]))
  if(settings$action=='fit') {
    identity<-api$ocean_identity(repository_dir,data_dir,settings)
    idfile<-F('ocean_path')(out,'workflow_identity.rds')
    # Never silently split the two method groups into different output roots.
    if(file.exists(idfile))api$ocean_assert_identity(readRDS(idfile),identity,out)
    api$ocean_storage_init(out)
  }
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  lock<-file.path(out,'.application.lock')
  if(!dir.create(lock,showWarnings=FALSE))stop('Ocean application already running: ',lock)
  writeLines(as.character(Sys.getpid()),file.path(lock,'pid'));on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  job<-api$ocean_job_dir(out,settings$run_id)
  old<-options(eivGP.core_events=file.path(job,'core-events'));on.exit(options(old),add=TRUE)
  release<-F('mixedgp_budget_enter')(settings$cores);on.exit(release(),add=TRUE,after=FALSE)
  result<-api$ocean_logged(file.path(job,'execution.log'),{
    if(settings$action=='fit') {
      api$ocean_execute(repository_dir,data_dir,out,settings)
      settings$action<-if(group=='eivGP')'analyze' else 'report'
      settings$run_id<-paste0(settings$run_id,'_scores')
    }
    a<-api$ocean_execute(repository_dir,data_dir,out,settings)
    dest<-F('ocean_path')(out,'reports',group)
    dir.create(dest,recursive=TRUE,showWarnings=FALSE)
    keep<-function(d) {
      if(!nrow(d))return(d)
      d[(d$variant=='EIV')==(group=='eivGP'),,drop=FALSE]
    }
    prediction<-a$prediction
    prediction$fold_metrics<-keep(prediction$fold_metrics)
    prediction$fold_metrics_SD_Y<-keep(prediction$fold_metrics_SD_Y)
    prediction$summary<-prediction$summary[(prediction$summary$variant=='EIV')==(group=='eivGP'),,drop=FALSE]
    F('ocean_export_prediction_report')(prediction,dest)
    full_traces<-if(group=='eivGP')F('mixedgp_application_trace_report')(out,file.path(dest,'scientific-traces'))else character()
    sections<-F('mixedgp_table_html')(prediction$summary,paste(group,'predictive performance'))
    for(p in full_traces[grepl('[.]png$',full_traces)])sections<-paste0(sections,
      "<h2>Full-history scientific targets</h2><img style='width:100%' src='scientific-traces/",basename(p),"'>")
    if(group=='competitors') {
      cc<-F('ocean_cc_prediction_tables')(out,a$inputs,a$prediction,
        cc_output_dir=if(nzchar(getp('cc_output_dir','')))api$ocean_resolve_path(params$cc_output_dir,repository_dir) else out,
        profile=profile,dest=file.path(dest,'complete_case'),expected_folds=settings$folds)
      F('mixedgp_atomic_csv')(cc$cc,file.path(dest,'complete_case','competitor_fold_metrics.csv'))
      sections<-c(sections,F('mixedgp_table_html')(cc$cc,'Complete-case competitors: fold scores'))
    } else if(profile!='quick' && all(a$prediction$availability$available)) {
      diagnostics_dir<-if(nzchar(getp('diagnostics_dir','')))api$ocean_resolve_path(params$diagnostics_dir,repository_dir) else F('ocean_path')(out,'diagnostics')
      d<-F('ocean_run_diagnostics')(out,data_dir,report_dir=diagnostics_dir,
        representative_fold=getp('representative_fold',1L),build_prediction=FALSE,
        build_surfaces=FALSE,prediction_report=prediction)
      targets<-F('ocean_build_target_traces')(out,a$inputs,settings)
      F('mixedgp_save_plot_jobs')(list(fig_ocean500_target_traces=F('ocean_target_trace_plot')(targets$traces)),dest,function(n)c(7.1,6.5))
      F('mixedgp_atomic_csv')(targets$summary,file.path(dest,'target_diagnostics.csv'))
      sections<-c(sections,F('mixedgp_table_html')(d$diagnostics$u$summary,'Masked oxygen'),
        '<h2>Scientific quantities: chain traces</h2><img src="fig_ocean500_target_traces.png" style="max-width:100%">',
        F('mixedgp_table_html')(targets$summary,'Target-specific diagnostics'))
    }
    F('ocean_atomic_lines')(c('<!doctype html><meta charset="utf-8">',paste0('<title>Ocean ',group,'</title>'),
      paste0('<h1>Ocean ',group,' — ',dataset,' / ',profile,'</h1>'),
      '<p>Fold means and SDs describe this fixed partition. Workflow completion does not establish convergence.</p>',sections),file.path(dest,'index.html'))
    F('mixedgp_atomic_rds')(list(group=group,settings=settings,no_refit=params$action=='report',
      input_hashes=tools::md5sum(list.files(data_dir,'[.]csv$',full.names=TRUE)),
      fit_identity=readRDS(F('ocean_path')(out,'workflow_identity.rds')),
      analysis_code=F('eivgp_module_identity')(),time=Sys.time()),file.path(dest,'provenance.rds'))
    list(path=dest,html=file.path(dest,'index.html'))
  })
  cat('Saved ',group,' analysis: ',result$html,'\n',sep='');invisible(result)
}

ocean_method_main <- function(group,args=commandArgs(TRUE)) {
  script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
  root<-normalizePath(file.path(dirname(normalizePath(script)),'..'))
  ocean_load_api(root,environment(ocean_method_main))
  if(identical(args,'--help')) {
    cat('Rscript tools/ocean-',group,'.R --cores=16 --profile=paper --dataset=ocean500 --output=/path/to/results\n',
      '  --reuse-only=true: require saved fits and predictions; never fit or continue.\n',
      '  --plan=true: list selected tasks. --folds=1:5 --cells=O0,O20,O50 --diagnostics-dir=...\n',
      '  EIV-GP and competitors each fit and analyze their own group; report only assembles saved results.\n',sep='')
    return(invisible(NULL))
  }
  Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
  ocean_method_run(root,ocean_method_params(args,group),group)
}

# Standalone R data preparation. jsonlite is used only for metadata.

# 1. Input-balanced, class-stratified five-fold assignment.
ocean_fold_permutations <- function(x) {
  if(length(x)==1L)return(matrix(x,nrow=1L))
  do.call(rbind,lapply(seq_along(x),function(i)
    cbind(x[i],ocean_fold_permutations(x[-i]))))
}

ocean_fold_standardize <- function(x, floor_sd=0) {
  centered<-sweep(x,2L,colMeans(x),"-")
  scale<-pmax(sqrt(colMeans(centered^2)),floor_sd)
  stopifnot(all(is.finite(scale)),all(scale>0))
  sweep(centered,2L,scale,"/")
}

ocean_assign_folds <- function(cohort,rule,seed) {
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(seed)
  n<-nrow(cohort);per_class<-n/20L
  stopifnot(n %in% c(500L,1000L),!anyDuplicated(cohort$profile),
    all(table(factor(cohort$C,levels=1:4))==n/4L))
  features<-matrix(0,n,9L);groups<-list()
  quintets<-function(ids)unname(split(ids,rep(seq_len(length(ids)/5L),each=5L)))
  for(class in 1:4) {
    ids<-which(cohort$C==class)
    x<-cbind(cohort$U[ids],log(cohort$G2pressure[ids]),cohort$G2salinity[ids])
    z<-ocean_fold_standardize(x)
    a<-cbind(z,z*z,z[,1]*z[,2],z[,1]*z[,3],z[,2]*z[,3])
    features[ids,]<-ocean_fold_standardize(a,1e-9)
    ordered<-ids[order(cohort$U[ids],cohort$source_row[ids])]
    if(rule=="joint_quintile_pressure") {
      oxygen_groups<-split(ordered,rep(1:5,each=length(ordered)/5L))
      for(group in oxygen_groups) {
        pressure_order<-group[order(cohort$G2pressure[group],cohort$source_row[group])]
        groups<-c(groups,quintets(pressure_order))
      }
    } else if(rule=="oxygen_interleaved")groups<-c(groups,quintets(ordered))
    else stop("Unknown fold rule: ",rule)
  }
  stopifnot(all(lengths(groups)==5L))
  assignment<-integer(n);sums<-array(0,c(5L,9L,4L))
  permutations<-ocean_fold_permutations(1:5)
  for(group in groups) {
    class<-cohort$C[group[1]]
    permutations<-permutations[sample.int(nrow(permutations)),,drop=FALSE]
    scores<-vapply(seq_len(nrow(permutations)),function(i)
      sum((sums[,,class]+features[group[permutations[i,]],,drop=FALSE])^2),numeric(1))
    ordered<-group[permutations[which.min(scores),]]
    assignment[ordered]<-1:5
    sums[,,class]<-sums[,,class]+features[ordered,,drop=FALSE]
  }
  initial<-sum((sums/per_class)^2);accepted<-0L
  for(step in seq_len(20000L)) {
    group<-groups[[sample.int(length(groups),1L)]]
    pair<-group[sample.int(5L,2L)];i<-pair[1];j<-pair[2]
    fi<-assignment[i];fj<-assignment[j];class<-cohort$C[i]
    delta<-features[j,]-features[i,];a<-sums[fi,,class];b<-sums[fj,,class]
    if(sum((a+delta)^2)+sum((b-delta)^2)<sum(a^2)+sum(b^2)-1e-12) {
      sums[fi,,class]<-a+delta;sums[fj,,class]<-b-delta
      assignment[i]<-fj;assignment[j]<-fi;accepted<-accepted+1L
    }
  }
  cohort$fold<-assignment
  stopifnot(all(table(cohort$fold,cohort$C)==per_class))
  list(cohort=cohort,design=list(fold_rule=rule,fold_seed=seed,
    initial_design_score=initial,final_design_score=sum((sums/per_class)^2),
    accepted_input_only_swaps=accepted))
}

# 2. Nested calibration from recorded C and X only. No U, Y or U-derived cell.
ocean_calibration_bins <- function(x,k) {
  # Average empirical ranks keep equal recorded values together. Splitting ties
  # into artificial strata would overstate coverage in the observed X space.
  n<-nrow(x)
  vapply(seq_len(ncol(x)),function(j)
    as.integer(pmin(k,1L+floor(k*(rank(x[,j],ties.method="average")-.5)/n))),integer(n))
}
ocean_calibration_select <- function(x,k,seed) {
  n<-nrow(x);stopifnot(length(k)==2L,all(k>0L),k[1]<k[2],k[2]<=n)
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(seed)
  bins<-lapply(k,function(size)ocean_calibration_bins(x,size))
  included<-function(label,stage)if(stage==1L)label==2L else label>0L
  counts<-function(label)lapply(1:2,function(stage)
    vapply(seq_len(ncol(x)),function(j)tabulate(bins[[stage]][included(label,stage),j],nbins=k[stage]),integer(k[stage])))
  loss<-function(ct)sum(vapply(ct,function(z)sum((z-1L)^2),0))
  # Return changes in marginal occupancy loss for exchanges of membership labels.
  # Labels: 0 hidden at both levels; 1 O50 only; 2 O20 and O50.
  delta<-function(i,j,label,ct) {
    ans<-numeric(length(i))
    for(stage in 1:2) {
      sign<-as.integer(included(label[i],stage))-as.integer(included(label[j],stage))
      for(axis in seq_len(ncol(x))) {
        a<-bins[[stage]][i,axis];b<-bins[[stage]][j,axis]
        ans<-ans+ifelse(a==b,0,2*sign*(ct[[stage]][b,axis]-ct[[stage]][a,axis])+2*sign^2)
      }
    }
    ans
  }
  exchange<-function(i,j,label,ct) {
    for(stage in 1:2) {
      sign<-as.integer(included(label[i],stage))-as.integer(included(label[j],stage))
      if(sign!=0L)for(axis in seq_len(ncol(x))) {
        a<-bins[[stage]][i,axis];b<-bins[[stage]][j,axis]
        ct[[stage]][a,axis]<-ct[[stage]][a,axis]-sign
        ct[[stage]][b,axis]<-ct[[stage]][b,axis]+sign
      }
    }
    old<-label[i];label[i]<-label[j];label[j]<-old
    list(label=label,ct=ct)
  }
  best<-NULL;best_loss<-Inf;attempts<-0L;steps<-max(2000L,30L*n)
  # Fixed design-only search budget, never a search over Y, U or model scores.
  for(restart in 1:3) {
    attempts<-restart;label<-integer(n);order<-sample.int(n)
    label[order[seq_len(k[2])]]<-1L;label[order[seq_len(k[1])]]<-2L
    ct<-counts(label);value<-loss(ct)
    remember<-function() {
      if(value<best_loss){best<<-label;best_loss<<-value}
    }
    remember()
    for(step in seq_len(steps)) {
      if(best_loss==0)break
      pair<-sample.int(n,2L);i<-pair[1];j<-pair[2]
      if(label[i]==label[j])next
      change<-delta(i,j,label,ct);temperature<-2*(.01/2)^((step-1)/max(1,steps-1))
      if(change<=0||runif(1)<exp(-change/temperature)) {
        state<-exchange(i,j,label,ct);label<-state$label;ct<-state$ct;value<-value+change;remember()
      }
    }
    # Polish the best candidate using the best improving membership exchange.
    label<-best;ct<-counts(label);value<-loss(ct)
    repeat {
      if(value==0)break
      pairs<-which(outer(label,label,"!=") & upper.tri(matrix(FALSE,n,n)),arr.ind=TRUE)
      changes<-delta(pairs[,1],pairs[,2],label,ct);minimum<-min(changes)
      if(minimum>=0)break
      tied<-which(changes==minimum);take<-tied[sample.int(length(tied),1L)]
      state<-exchange(pairs[take,1],pairs[take,2],label,ct);label<-state$label;ct<-state$ct;value<-value+minimum;remember()
    }
    if(best_loss==0)break
  }
  ct<-counts(best)
  audit<-do.call(rbind,lapply(1:2,function(stage)do.call(rbind,lapply(seq_len(ncol(x)),function(j) {
    occupancy<-ct[[stage]][,j]
    data.frame(calibration=c("O20","O50")[stage],axis=colnames(x)[j],n_available=n,n_selected=k[stage],
      quantile_bins=k[stage],occupied_bins=sum(occupancy>0L),max_bin_count=max(occupancy),
      squared_occupancy_loss=sum((occupancy-1L)^2),exact_latin_margin=all(occupancy==1L),
      distinct_X_values=length(unique(x[,j])),search_restarts=attempts)
  }))))
  stopifnot(sum(best==2L)==k[1],sum(best>0L)==k[2],identical(loss(ct),best_loss))
  list(o20=best==2L,o50=best>0L,audit=audit)
}
ocean_fold_calibration <- function(train,seed) {
  # Access is deliberately limited to recorded class, salinity/log-pressure,
  # and fixed row identifiers. Hidden oxygen and U-derived `cell` are ignored.
  required<-c("C","G2salinity","G2pressure")
  if(!all(required%in%names(train)))stop("Calibration requires C, G2salinity and G2pressure.")
  if(anyNA(train[,required])||any(!is.finite(as.matrix(train[,required])))||any(train$G2pressure<=0))
    stop("Calibration requires finite recorded C/X and positive pressure.")
  if(!all(train$C%in%1:4)||!all(1:4%in%train$C))stop("Expected four recorded classes.")
  stable<-if("cohort_id"%in%names(train))train$cohort_id else if("source_row"%in%names(train))train$source_row else seq_len(nrow(train))
  if(anyNA(stable)||anyDuplicated(stable))stop("Calibration row identifiers must be unique.")
  o20<-o50<-rep(FALSE,nrow(train));audit<-list()
  for(class in 1:4) {
    ids<-which(train$C==class);ids<-ids[order(stable[ids])];n<-length(ids)
    k<-as.integer(floor(n*c(.2,.5)+.5))
    x<-cbind(salinity=train$G2salinity[ids],log_pressure=log(train$G2pressure[ids]))
    selected<-ocean_calibration_select(x,k,seed+class)
    o20[ids]<-selected$o20;o50[ids]<-selected$o50
    audit[[class]]<-cbind(C=class,selected$audit)
  }
  train$calib_o20<-as.integer(o20);train$calib_o50<-as.integer(o50)
  attr(train,"calibration_design")<-list(method="C-stratified nested empirical-rank X coverage v1",
    selection_variables=c("C","G2salinity","log(G2pressure)"),seed=seed,
    optimizer="three fixed-budget seeded exchange searches; marginal occupancy loss only",
    strict_LHS_guaranteed=FALSE,audit=do.call(rbind,audit))
  stopifnot(all(train$calib_o20<=train$calib_o50))
  train
}

# 3. Nested class-balanced test-input masks, using the existing R seed rule.
ocean_fold_test_masks <- function(test,fold) {
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(202609260L+fold)
  o20<-o50<-rep(FALSE,nrow(test))
  for(class in 1:4) {
    ids<-which(test$C==class);ids<-ids[sample.int(length(ids))];n<-length(ids)
    n20<-round(.2*n)
    n50<-floor(.5*n)+as.integer(n%%2L==1L && (class+fold)%%2L==0L)
    o20[head(ids,n20)]<-TRUE;o50[head(ids,n50)]<-TRUE
  }
  data.frame(cohort_id=test$cohort_id,C=test$C,observed_U_O20=o20,
    observed_U_O50=o50,mask_seed=202609260L+fold)
}

prepare_ocean_folds <- function(dataset="ocean1000",output_dir=NULL,repository_dir=getwd(),calibration_only=FALSE) {
  dataset<-match.arg(dataset,c("ocean500","ocean1000"))
  if(!requireNamespace("jsonlite",quietly=TRUE))stop("Install the R package jsonlite for metadata.")
  root<-normalizePath(repository_dir,mustWork=TRUE)
  input<-file.path(root,"data/ocean_silicate",sub("ocean","n",dataset))
  if(is.null(output_dir))output_dir<-file.path(root,"reproduction/ocean/prepared",
    paste0(dataset,if(calibration_only)"_xcalibration_v1"else""))
  output_dir<-normalizePath(path.expand(output_dir),mustWork=FALSE)
  if(identical(normalizePath(input),output_dir))stop("Choose an output directory separate from the source cohort.")
  if(dir.exists(output_dir)&&length(list.files(output_dir,all.files=TRUE,no..=TRUE)))
    stop("Prepared output already exists; use a new directory: ",output_dir)
  meta_file<-file.path(input,if(dataset=="ocean1000")"cohort_metadata.json" else "freeze_manifest.json")
  metadata<-jsonlite::read_json(meta_file,simplifyVector=FALSE)
  cohort<-read.csv(file.path(input,"cohort_manifest.csv"),stringsAsFactors=FALSE)
  n<-nrow(cohort);seed<-as.integer(metadata$sample_seed)
  stopifnot(n==as.integer(sub("ocean","",dataset)))
  source_files<-file.path(input,"cohort_manifest.csv")
  if(calibration_only) {
    # Validate the existing partitions once, then preserve all rows and test masks.
    frozen<-ocean_prepare_inputs(input)
    source_files<-c(source_files,list.files(input,"^fold[0-9]+_(train|test|test_U_masks)[.]csv$",full.names=TRUE))
    folded<-list(design=list(fold_preparation="Preserved supplied frozen folds; only training calibration replaced"))
  } else {
    cohort$fold<-NULL
    rule<-if(n==500L)"joint_quintile_pressure" else "oxygen_interleaved"
    folded<-ocean_assign_folds(cohort,rule,seed+if(n==500L)100L else 120L)
    cohort<-folded$cohort
  }
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  if(calibration_only)stopifnot(file.copy(file.path(input,"cohort_manifest.csv"),output_dir))
  else write.csv(cohort,file.path(output_dir,"cohort_manifest.csv"),row.names=FALSE)
  audits<-list()
  for(f in 1:5) {
    if(calibration_only) {
      tr<-read.csv(file.path(input,sprintf("fold%02d_train.csv",f)),stringsAsFactors=FALSE)
      test<-read.csv(file.path(input,sprintf("fold%02d_test.csv",f)),stringsAsFactors=FALSE)
      masks<-read.csv(file.path(input,sprintf("fold%02d_test_U_masks.csv",f)),stringsAsFactors=FALSE)
    } else {
      tr<-cohort[cohort$fold!=f,,drop=FALSE];test<-cohort[cohort$fold==f,,drop=FALSE]
      masks<-ocean_fold_test_masks(test,f)
    }
    train<-ocean_fold_calibration(tr,seed+200L+f)
    design<-attr(train,"calibration_design");audits[[f]]<-cbind(fold=f,design$audit)
    stopifnot(!length(intersect(train$profile,test$profile)),
      sum(masks$observed_U_O20)==.2*nrow(test),sum(masks$observed_U_O50)==.5*nrow(test),
      all(masks$observed_U_O20<=masks$observed_U_O50))
    write.csv(train,file.path(output_dir,sprintf("fold%02d_train.csv",f)),row.names=FALSE)
    if(calibration_only) {
      for(suffix in c("test","test_U_masks"))stopifnot(file.copy(file.path(input,sprintf("fold%02d_%s.csv",f,suffix)),output_dir))
    } else {
      write.csv(test,file.path(output_dir,sprintf("fold%02d_test.csv",f)),row.names=FALSE)
      write.csv(masks,file.path(output_dir,sprintf("fold%02d_test_U_masks.csv",f)),row.names=FALSE)
    }
  }
  audit<-do.call(rbind,audits);write.csv(audit,file.path(output_dir,"calibration_X_coverage.csv"),row.names=FALSE)
  metadata[names(folded$design)]<-folded$design
  metadata$data_stage<-"five-fold inputs with C/X-only nested calibration"
  metadata$train<-.8*n;metadata$test<-n/5L
  metadata$calibration_design<-design[!names(design)%in%c("audit","seed")]
  metadata$calibration_design$seed_by_fold<-seed+200L+1:5
  metadata$calibration_uses_U<-FALSE;metadata$calibration_uses_Y<-FALSE;metadata$calibration_uses_cell<-FALSE
  metadata$calibration_seed_rule<-"sample_seed+200+fold; class stream adds C"
  metadata$test_mask_seed_rule<-"202609260+fold"
  metadata$preparation_implementation<-"base-R-X-calibration-v1"
  metadata$R_version<-as.character(getRversion());metadata$RNG_kind<-RNGkind()
  metadata$source_files_md5<-as.list(setNames(unname(tools::md5sum(source_files)),basename(source_files)))
  metadata$preparation_script_md5<-unname(tools::md5sum(file.path(root,"tools/ocean-tools.R")))
  jsonlite::write_json(metadata,file.path(output_dir,"freeze_manifest.json"),auto_unbox=TRUE,pretty=TRUE,digits=16)
  # Run the same reader used by fitting before advertising a prepared directory.
  validated<-ocean_prepare_inputs(output_dir)
  message("Prepared ",n," profiles: ",output_dir,
    "\nCalibration uses C and recorded X only; O20 is nested in O50.",
    "\nExact empirical Latin margins: ",sum(audit$exact_latin_margin),"/",nrow(audit),
    "; see calibration_X_coverage.csv. Strict LHS is not guaranteed.")
  invisible(output_dir)
}
ocean_prepare_main <- function(args,repository_dir) {
  if(identical(args,"--help")) {
    cat("Rscript tools/ocean-prepare-folds.R ocean500 [OUTPUT_DIR] --calibration-only\n",
      "  Preserve the supplied cohort, folds and test masks; replace training calibration using recorded C/X.\n",
      "Rscript tools/ocean-prepare-folds.R ocean1000 [OUTPUT_DIR]\n",
      "  Prepare all folds and C/X-only calibration from the raw selected cohort.\n",sep="")
    return(invisible(NULL))
  }
  flag<-args=="--calibration-only";if(sum(flag)>1L)stop("Repeated --calibration-only")
  positional<-args[!flag]
  if(length(positional)>2L||any(startsWith(positional,"--")))stop("Use --help for Ocean preparation arguments.")
  prepare_ocean_folds(if(length(positional))positional[1]else"ocean1000",
    if(length(positional)>1L)positional[2]else NULL,repository_dir=repository_dir,calibration_only=any(flag))
}



  },e)
  e$ocean_load_api<-function(repository_dir=NULL,envir=parent.frame()) {
    list2env(as.list(e),envir=envir);invisible(NULL)
  }
  e$ocean_workflow_api<-function()as.list(e)
  e$ocean_analysis<-function(repository_dir,data_dir,output_dir,settings,...)
    e$run_ocean_analysis(data_dir,output_dir,repository_dir=repository_dir,settings=settings,...)
  # Application functions may call shared package helpers by their generic names.
  for(n in ls(ns,all.names=TRUE))if(startsWith(n,'mixedgp_'))e[[n]]<-get(n,ns)
  as.list(e)
}
ocean_load_api <- function(repository_dir=.ocean_repository,envir=parent.frame()) {
  list2env(ocean_workflow_api(repository_dir),envir=envir);invisible(NULL)
}
local({
  bindings<-ocean_workflow_api()
  bindings[c('ocean_load_api','ocean_workflow_api')]<-NULL
  list2env(bindings,envir=parent.env(environment()))
})
