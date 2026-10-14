# Standalone application layer; upstream package source is never changed here.
ocean_load_helpers <- function(repository_dir,envir=parent.frame()) {
  helper<-file.path(repository_dir,"paper/applications/ocean_helpers")
  for(name in c("01_setup.R","02_fit_models.R","03_reporting.R","04_diagnostics.R"))
    sys.source(file.path(helper,name),envir)
  invisible(helper)
}
ocean_check_package <- function(repository_dir) {
  if(!requireNamespace("eivGP",quietly=TRUE) || as.character(getNamespaceVersion("eivGP"))!="0.4.0")
    stop("Use the pinned upstream eivGP 0.4.0 library; restart R if an older namespace is loaded.")
  pin<-read.csv(file.path(repository_dir,"paper/applications/ocean_helpers/UPSTREAM_PACKAGE.csv"),stringsAsFactors=FALSE)
  installed<-read.csv(system.file("module-manifest.csv",package="eivGP"),stringsAsFactors=FALSE)
  if(!setequal(pin$module,installed$module) || !identical(pin$md5,installed$md5[match(pin$module,installed$module)]))
    stop("Installed package does not match UPSTREAM_PACKAGE.csv; install the fixed upstream commit.")
  invisible(pin)
}
ocean_analysis <- function(repository_dir,data_dir,output_dir,settings,import_dirs=character(),
  cc_output_dir="",diagnostics_dir="",representative_fold=1L,diagnostics_workers=1L) {
  ocean_check_package(repository_dir)
  a<-ocean_execute(repository_dir,data_dir,output_dir,settings,import_dirs)
  if(settings$action=="plan")return(a)
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
      file.exists(file.path(output_dir,"full5_class4",sprintf("fold%02d",f),"baseline_package_outputs.rds")),logical(1)))
    if(ready) {
      d<-ocean_run_diagnostics(output_dir,data_dir,
        report_dir=if(nzchar(diagnostics_dir))diagnostics_dir else file.path(output_dir,"diagnostics"),
        representative_fold=representative_fold,workers=diagnostics_workers,
        build_prediction=FALSE,prediction_report=available,
        build_surfaces="surface_comparison" %in% settings$modules)
      a$diagnostic_run<-d;a$diagnostics<-d$diagnostics;a$surfaces<-d$surfaces
    }
  }
  a
}
if(sys.nframe()==0L) {
  entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))
  repository_dir<-normalizePath(file.path(dirname(entry),"../.."),mustWork=TRUE)
  ocean_load_helpers(repository_dir,globalenv())
  args<-commandArgs(TRUE);p<-if(length(args))readRDS(args[1]) else list(action="plan")
  settings<-do.call(ocean_settings,p[intersect(names(p),names(formals(ocean_settings)))])
  data_dir<-if(!is.null(p$data_dir) && nzchar(p$data_dir))ocean_resolve_path(p$data_dir,repository_dir) else file.path(repository_dir,"data/ocean_silicate")
  output_dir<-if(!is.null(p$output_dir) && nzchar(p$output_dir))ocean_resolve_path(p$output_dir,repository_dir) else file.path(repository_dir,"reproduction/ocean",settings$profile)
  imports<-as.character(unlist(p$import_dirs));imports<-ocean_resolve_path(imports[nzchar(imports)],repository_dir)
  cc<-if(!is.null(p$cc_output_dir) && nzchar(p$cc_output_dir))ocean_resolve_path(p$cc_output_dir,repository_dir) else ""
  result<-ocean_analysis(repository_dir,data_dir,output_dir,settings,imports,cc_output_dir=cc)
  if(settings$action=="plan")print(result$tasks,row.names=FALSE)
  else if(!is.null(result$prediction))print(result$prediction$summary,row.names=FALSE)
}
