# Ocean application entrypoint: implementation is provided by the installed package.
ocean_load_api <- function(repository_dir=NULL,envir=parent.frame()) {
  if(!requireNamespace("eivGP",quietly=TRUE)||utils::packageVersion("eivGP")<"0.4.0.9002")
    stop("Install this checkout with Rscript tools/build_package.R and R CMD INSTALL eivGP.")
  list2env(eivGP::ocean_workflow_api(),envir=envir);invisible(NULL)
}
ocean_analysis <- function(repository_dir,data_dir,output_dir,settings,...) {
  eivGP::run_ocean_analysis(data_dir,output_dir,repository_dir=repository_dir,settings=settings,...)
}
# Dataset preparation is separate from fitting and report rendering.
ocean_dataset_dir <- function(repository_dir,dataset="ocean200",data_dir="") {
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
if(sys.nframe()==0L) {
  args<-commandArgs(TRUE)
  root<-normalizePath(file.path(dirname(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))),"../.."))
  ocean_load_api(root,globalenv())
  p<-if(length(args))readRDS(args[1]) else list(action="plan")
  settings<-do.call(ocean_settings,p[intersect(names(p),names(formals(ocean_settings)))])
  dataset<-if(is.null(p$dataset))"ocean200" else p$dataset
  data_dir<-ocean_dataset_dir(root,dataset,if(is.null(p$data_dir))"" else p$data_dir)
  output_dir<-if(!is.null(p$output_dir)&&nzchar(p$output_dir))ocean_resolve_path(p$output_dir,root) else
    file.path(root,"reproduction/ocean",dataset,settings$profile)
  result<-ocean_analysis(root,data_dir,output_dir,settings)
  if(settings$action=="plan")print(result$tasks,row.names=FALSE)
}
