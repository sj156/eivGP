ocean_path <- function(...) getFromNamespace("ocean_path","eivGP")(...)
# Ocean application entrypoint: implementation is provided by the installed package.
ocean_load_api <- function(repository_dir=NULL,envir=parent.frame()) {
  if(!is.null(repository_dir)) {
    lib <- Sys.getenv("EIVGP_R_LIBRARY",file.path(repository_dir,".r-library"))
    if(dir.exists(lib)).libPaths(c(normalizePath(lib),.libPaths()))
  }
  if(!requireNamespace("eivGP",quietly=TRUE)||utils::packageVersion("eivGP")<"0.4.0.9008")
    stop("Install this checkout with Rscript tools/build_package.R and R CMD INSTALL eivGP.")
  if(!is.null(repository_dir)) {
    installed<-utils::read.csv(system.file("module-manifest.csv",package="eivGP"))
    native<-grepl("[.]cpp$|^Makevars",installed$module)
    paths<-file.path(repository_dir,"eivGP",ifelse(native,"src","R"),installed$module)
    if(!identical(unname(tools::md5sum(paths)),installed$md5))
      stop("Installed eivGP differs from this checkout; rebuild and reinstall it first.")
  }
  list2env(eivGP::ocean_workflow_api(),envir=envir);invisible(NULL)
}
ocean_analysis <- function(repository_dir,data_dir,output_dir,settings,...) {
  eivGP::run_ocean_analysis(data_dir,output_dir,repository_dir=repository_dir,settings=settings,...)
}
# Dataset preparation is separate from fitting and report rendering.
ocean_dataset_dir <- function(repository_dir,dataset="ocean500",data_dir="") {
  dataset<-match.arg(dataset,c("ocean200","ocean500","ocean1000"))
  if(nzchar(data_dir))return(ocean_resolve_path(data_dir,repository_dir))
  if(dataset=="ocean1000") {
    path<-file.path(repository_dir,"reproduction/ocean/prepared/ocean1000")
    if(!file.exists(file.path(path,"fold01_train.csv")))
      stop("Prepare ocean1000 first: Rscript data/ocean_silicate/prepare_folds.R ocean1000; see README_ocean.md.")
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
    "diagnostics_dir","representative_fold","continuation_iter")
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
  for(key in c("chains","iter","burn"))params[[key]]<-saved[[key]]
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
    output_file="ocean.html",output_dir=output_dir,
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
          api$ocean_atomic_rds(output_dir,paste0(ocean_default_output(dataset,profile,repository_dir),".latest-run.rds"))
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
    print(api$ocean_core_plan(settings$cores,sum(result$tasks$module=="eiv"),settings$chains))
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
    on.exit(release_budget(),add=TRUE)
    # Fail early if HTML cannot be built, before starting an expensive fit.
    ocean_check_report_dependencies()
    result<-api$ocean_logged(file.path(api$ocean_job_dir(output_dir,settings$run_id),"execution.log"),
      ocean_cli_workflow(settings,execute,render))
    cat("Saved results:",output_dir,"\nHTML:",result$html,"\n")
  }
  invisible(result)
}

if(sys.nframe()==0L) {
  args<-commandArgs(TRUE)
  if(identical(args,"--help")) {
    cat("Rscript paper/applications/ocean.R --action=plan|fit|analyze|continue|report --profile=paper|pilot|quick\n",
      "  --dataset=ocean200|ocean500|ocean1000 --cores=16 --output=/path/to/results\n",
      "  Default output: results/<dataset>/ (paper), results/<dataset>_<profile>/ otherwise.\n",
      "  Optional: --data-root=... --folds=1:5 --cells=O0,O20,O50 --modules=...\n",
      "  --run-id=... --diagnostics-dir=...\n",
      "  Sampling is fixed by profile: paper 4 chains/50000 iterations/5000 burn-in; pilot 4/80/30; quick 1/80/30.\n",
      "  Incompatible initial-fit results are preserved; a new folder is selected and printed.\n",
      "  --continuation-iter=20000 (additional iterations per existing chain; pilot default: 80)\n",
      "  --task-ids=... --import-dirs=... --cc-output-dir=... --representative-fold=1\n",
      "No arguments: print the plan. --cores=N alone: fit, analyze, export paper assets and render ocean500.\n",
      "fit (also run) fits, analyzes with the shared core budget and writes HTML; continue extends, analyzes and updates HTML.\n",
      "cores is the shared budget across fitting, prediction, diagnostics and plotting. analyze/report rebuild HTML without fitting. See README_ocean.md.\n",sep="")
  } else {
    script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
    root<-normalizePath(file.path(dirname(normalizePath(script)),"../.."))
    Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",VECLIB_MAXIMUM_THREADS="1")
    lib<-Sys.getenv("EIVGP_R_LIBRARY",file.path(root,".r-library"))
    if(dir.exists(lib)).libPaths(c(normalizePath(lib),.libPaths()))
    ocean_run_cli(root,ocean_cli_params(args))
  }
}
