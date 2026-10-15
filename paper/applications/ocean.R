# Ocean application entrypoint: implementation is provided by the installed package.
ocean_load_api <- function(repository_dir=NULL,envir=parent.frame()) {
  if(!requireNamespace("eivGP",quietly=TRUE)||utils::packageVersion("eivGP")<"0.4.0.9007")
    stop("Install this checkout with Rscript tools/build_package.R and R CMD INSTALL eivGP.")
  if(!is.null(repository_dir)) {
    installed<-utils::read.csv(system.file("module-manifest.csv",package="eivGP"))
    paths<-file.path(repository_dir,"eivGP/R",installed$module)
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
ocean_output_dir <- function(repository_dir,dataset,profile,output_dir="") {
  if(nzchar(output_dir))ocean_resolve_path(output_dir,repository_dir) else
    file.path(repository_dir,"results/ocean",dataset,profile)
}

# CLI options select the existing workflow; they do not alter model fitting.
ocean_cli_params <- function(args=character()) {
  if(length(args)==1L && !startsWith(args,"--"))return(readRDS(args))
  aliases<-c("data-root"="data_dir","output"="output_dir")
  allowed<-c("action","profile","dataset","cores","chains","iter","burn","data_dir","output_dir",
    "folds","cells","modules","task_ids","run_id","import_dirs","cc_output_dir",
    "diagnostics_dir","representative_fold","continuation_iter")
  p<-list()
  for(arg in args) {
    if(!grepl("^--[^=]+=.+$",arg))stop("Use --name=value: ",arg)
    key<-sub("^--([^=]+)=.*$","\\1",arg)
    key<-if(key%in%names(aliases))unname(aliases[[key]]) else gsub("-","_",key)
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
  p
}

# Existing fits own their original sampling settings; extensions only add steps.
ocean_restore_fit_settings <- function(params,output_dir) {
  path<-file.path(output_dir,"workflow_identity.rds")
  if(!file.exists(path))return(params)
  saved<-readRDS(path)$numerics
  for(key in c("chains","iter","burn")) {
    if(is.null(params[[key]]))params[[key]]<-saved[[key]] else
      if(identical(params$action,"continue") && !identical(as.double(params[[key]]),as.double(saved[[key]])))
        stop("Cannot change ",key," during continuation (saved ",saved[[key]],
          "). Use --continuation-iter to add iterations per existing chain.")
  }
  params
}

ocean_check_report_dependencies <- function() {
  if(!requireNamespace("rmarkdown",quietly=TRUE))stop("HTML reporting requires the R package rmarkdown.")
  if(!rmarkdown::pandoc_available())stop("HTML reporting requires Pandoc on PATH (or RSTUDIO_PANDOC).")
  invisible(TRUE)
}

ocean_render_report <- function(repository_dir,output_dir,settings,dataset,data_dir,
  cc_output_dir="",diagnostics_dir="",representative_fold=1L) {
  rmarkdown::render(file.path(repository_dir,"paper/applications/ocean.Rmd"),
    output_file="ocean.html",output_dir=output_dir,
    params=list(dataset=dataset,profile=settings$profile,action="report",folds=settings$folds,
      cells=settings$cells,modules=settings$modules,task_ids=settings$task_ids,
      cores=1L,run_id=paste0(settings$run_id,"_report"),data_dir=data_dir,output_dir=output_dir,
      cc_output_dir=cc_output_dir,diagnostics_dir=diagnostics_dir,representative_fold=representative_fold),
    envir=new.env(parent=globalenv()),quiet=TRUE)
}

# Orchestration is separate from the package's fit/analyze implementations.
# A failed later stage leaves verified checkpoints available for resume.
ocean_cli_workflow <- function(settings,execute,render) {
  result<-execute(settings,TRUE)
  if(settings$action%in%c("fit","continue")) {
    analyzed<-settings;analyzed$action<-"analyze";analyzed$cores<-1L
    analyzed$run_id<-paste0(settings$run_id,"_analyze")
    analysis<-execute(analyzed,FALSE)
    result<-list(fitting=result,analysis=analysis)
  }
  if(settings$action!="plan")result$html<-render()
  result
}

ocean_run_cli <- function(repository_dir,params=list()) {
  if("diagnostics_workers"%in%names(params))stop("diagnostics_workers was removed; analysis always uses one worker.")
  api<-new.env(parent=environment())
  ocean_load_api(repository_dir,api)
  getp<-function(key,default)if(is.null(params[[key]]))default else params[[key]]
  dataset<-match.arg(getp("dataset","ocean500"),c("ocean200","ocean500","ocean1000"))
  params$run_id<-getp("run_id",paste(dataset,getp("profile","paper"),
    if(getp("action","plan")=="run")"fit" else getp("action","plan"),sep="_"))
  profile<-match.arg(getp("profile","paper"),c("paper","pilot","quick"))
  environment(ocean_dataset_dir)<-api
  environment(ocean_output_dir)<-api
  data_dir<-ocean_dataset_dir(repository_dir,dataset,getp("data_dir",""))
  output_dir<-ocean_output_dir(repository_dir,dataset,profile,getp("output_dir",""))
  params<-ocean_restore_fit_settings(params,output_dir)
  settings<-do.call(api$ocean_settings,
    params[intersect(names(params),names(formals(api$ocean_settings)))])
  if(settings$action%in%c("analyze","report"))settings$cores<-1L
  resolve<-function(paths)vapply(paths,function(path)
    if(nzchar(path))api$ocean_resolve_path(path,repository_dir) else "",
    character(1),USE.NAMES=FALSE)
  execute<-function(stage_settings,import=TRUE)ocean_analysis(repository_dir,data_dir,output_dir,stage_settings,
    import_dirs=if(import)resolve(getp("import_dirs",character())) else character(),
    cc_output_dir=resolve(getp("cc_output_dir","")),
    diagnostics_dir=resolve(getp("diagnostics_dir","")),
    representative_fold=getp("representative_fold",1L))
  render<-function()ocean_render_report(repository_dir,output_dir,settings,dataset,data_dir,
    resolve(getp("cc_output_dir","")),resolve(getp("diagnostics_dir","")),
    getp("representative_fold",1L))
  if(settings$action=="plan") {
    result<-execute(settings)
    print(result$tasks,row.names=FALSE)
    print(api$ocean_core_plan(settings$cores,sum(result$tasks$module=="eiv"),settings$chains))
    cat("Total fit core budget:",settings$cores,"; up to",settings$chains,
      "chain workers per EIV fit; analysis and HTML use one worker.\n")
  } else {
    # Fail early if HTML cannot be built, before starting an expensive fit.
    ocean_check_report_dependencies()
    result<-api$ocean_logged(file.path(output_dir,"jobs",settings$run_id,"execution.log"),
      ocean_cli_workflow(settings,execute,render))
    cat("Saved results:",output_dir,"\nHTML:",result$html,"\n")
  }
  invisible(result)
}

if(sys.nframe()==0L) {
  args<-commandArgs(TRUE)
  if(identical(args,"--help")) {
    cat("Rscript paper/applications/ocean.R --action=plan|fit|analyze|continue|report --profile=paper|pilot|quick\n",
      "  --dataset=ocean200|ocean500|ocean1000 --cores=16 --chains=4 --iter=100000 --output=/path/to/results\n",
      "  Optional: --data-root=... --folds=1:5 --cells=O0,O20,O50 --modules=...\n",
      "  --run-id=... --diagnostics-dir=... --burn=5000 (paper default)\n",
      "  --iter: initial total iterations per chain (paper 100000; pilot/quick 80)\n",
      "  --continuation-iter=20000 (additional iterations per existing chain; pilot default: 80)\n",
      "  --task-ids=... --import-dirs=... --cc-output-dir=... --representative-fold=1\n",
      "No arguments: validate inputs and print the paper-profile ocean500 plan.\n",
      "fit (also run) fits, analyzes serially and writes HTML; continue extends, analyzes and updates HTML.\n",
      "cores is the shared fit budget across tasks/chains. analyze/report rebuild HTML without fitting. See README_ocean.md.\n",sep="")
  } else {
    script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
    root<-normalizePath(file.path(dirname(normalizePath(script)),"../.."))
    Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",VECLIB_MAXIMUM_THREADS="1")
    lib<-Sys.getenv("EIVGP_R_LIBRARY",file.path(root,".r-library"))
    if(dir.exists(lib)).libPaths(c(normalizePath(lib),.libPaths()))
    ocean_run_cli(root,ocean_cli_params(args))
  }
}
