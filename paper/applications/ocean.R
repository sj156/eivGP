# Ocean application entrypoint: implementation is provided by the installed package.
ocean_load_api <- function(repository_dir=NULL,envir=parent.frame()) {
  if(!requireNamespace("eivGP",quietly=TRUE)||utils::packageVersion("eivGP")<"0.4.0.9003")
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
  allowed<-c("action","profile","dataset","cores","data_dir","output_dir",
    "folds","cells","modules","task_ids","run_id","import_dirs","cc_output_dir",
    "diagnostics_dir","diagnostics_workers","representative_fold")
  p<-list()
  for(arg in args) {
    if(!grepl("^--[^=]+=.+$",arg))stop("Use --name=value: ",arg)
    key<-sub("^--([^=]+)=.*$","\\1",arg)
    key<-if(key%in%names(aliases))unname(aliases[[key]]) else gsub("-","_",key)
    if(!key%in%allowed)stop("Unknown option: ",key)
    if(key%in%names(p))stop("Repeated option: ",key)
    p[[key]]<-sub("^--[^=]+=","",arg)
  }
  for(key in intersect(names(p),c("cores","diagnostics_workers","representative_fold"))) {
    value<-suppressWarnings(as.numeric(p[[key]]))
    if(length(value)!=1L || !is.finite(value) || value<1 || value!=as.integer(value))
      stop(key," must be a positive integer")
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

ocean_run_cli <- function(repository_dir,params=list()) {
  api<-new.env(parent=environment())
  ocean_load_api(repository_dir,api)
  getp<-function(key,default)if(is.null(params[[key]]))default else params[[key]]
  dataset<-match.arg(getp("dataset","ocean500"),c("ocean200","ocean500","ocean1000"))
  params$run_id<-getp("run_id",paste(dataset,getp("profile","paper"),sep="_"))
  settings<-do.call(api$ocean_settings,
    params[intersect(names(params),names(formals(api$ocean_settings)))])
  # Dataset/path helpers use the package API in this local environment.
  environment(ocean_dataset_dir)<-api
  environment(ocean_output_dir)<-api
  data_dir<-ocean_dataset_dir(repository_dir,dataset,getp("data_dir",""))
  output_dir<-ocean_output_dir(repository_dir,dataset,settings$profile,getp("output_dir",""))
  resolve<-function(paths)vapply(paths,function(path)
    if(nzchar(path))api$ocean_resolve_path(path,repository_dir) else "",
    character(1),USE.NAMES=FALSE)
  execute<-function()ocean_analysis(repository_dir,data_dir,output_dir,settings,
    import_dirs=resolve(getp("import_dirs",character())),
    cc_output_dir=resolve(getp("cc_output_dir","")),
    diagnostics_dir=resolve(getp("diagnostics_dir","")),
    representative_fold=getp("representative_fold",1L),
    diagnostics_workers=getp("diagnostics_workers",1L))
  if(settings$action=="plan") {
    result<-execute()
    print(result$tasks,row.names=FALSE)
  } else {
    result<-api$ocean_logged(file.path(output_dir,"jobs",settings$run_id,"execution.log"),execute())
    cat("Saved results:",output_dir,"\n")
    cat("Render ocean.Rmd with action='report' and the same dataset, profile and output_dir.\n")
  }
  invisible(result)
}

if(sys.nframe()==0L) {
  args<-commandArgs(TRUE)
  if(identical(args,"--help")) {
    cat("Rscript paper/applications/ocean.R --action=plan|run|report --profile=paper|pilot|quick\n",
      "  --dataset=ocean200|ocean500|ocean1000 --cores=4 --output=/path/to/results\n",
      "  Optional: --data-root=... --folds=1:5 --cells=O0,O20,O50 --modules=...\n",
      "  --run-id=... --diagnostics-workers=2 --diagnostics-dir=...\n",
      "  --task-ids=... --import-dirs=... --cc-output-dir=... --representative-fold=1\n",
      "No arguments: validate inputs and print the paper-profile ocean500 plan.\n",
      "run fits/resumes and saves results; report rebuilds summaries without fitting.\n",
      "Render ocean.Rmd separately for HTML. See README_ocean.md.\n",sep="")
  } else {
    script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
    root<-normalizePath(file.path(dirname(normalizePath(script)),"../.."))
    lib<-Sys.getenv("EIVGP_R_LIBRARY",file.path(root,".r-library"))
    if(dir.exists(lib)).libPaths(c(normalizePath(lib),.libPaths()))
    ocean_run_cli(root,ocean_cli_params(args))
  }
}
