#!/usr/bin/env Rscript
# Locate the shared bootstrap beside this script, including when source() is used.
.local_file <- if(sys.nframe()>0L && !is.null(sys.frame(1)$ofile))sys.frame(1)$ofile else sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(normalizePath(.local_file)),"common.R"))
# Rebuild the established studies; the minimal additions live in supplement.R.
run_study_cli <- function(args=commandArgs(TRUE)) {
a<-list();for(x in args){if(!grepl("^--[^=]+=.+$",x))stop("Use --name=value");a[[sub("^--([^=]+)=.*$","\\1",x)]]<-sub("^--[^=]+=","",x)}
get<-function(k,d)if(is.null(a[[k]]))d else a[[k]]
allowed<-c("study","profile","action","data-root","output","cores","competitors")
if(length(setdiff(names(a),allowed)))stop("Unknown option: ",paste(setdiff(names(a),allowed),collapse=","))
numerical_package_check()
study<-get("study","study2");action<-get("action","plan");profile<-get("profile","paper")
output<-get("output",Sys.getenv("EIVGP_RESULTS_ROOT",numerical_default_output(study)));data<-get("data-root",Sys.getenv("EIVGP_DATA_ROOT",numerical_default_data(study)))
if(!nzchar(output)||!nzchar(data))stop("Supply --output (local run root) and --data-root (study-specific frozen data root)")
if(grepl("Dropbox|Overleaf",output,ignore.case=TRUE))stop("Output must be outside the synced manuscript tree")
# The package data stage maintains manifests, so copy synced inputs locally first.
if(action!="plan" && grepl("Dropbox|Overleaf",data,ignore.case=TRUE)) {
  local<-file.path(output,"frozen-data",study)
  paths<-list.files(data,recursive=TRUE,full.names=FALSE)
  if(!length(paths))stop("No frozen input files found")
  for(rel in paths) {
    from<-file.path(data,rel);to<-file.path(local,rel)
    dir.create(dirname(to),recursive=TRUE,showWarnings=FALSE)
    if(file.exists(to)) {
      if(!identical(unname(tools::md5sum(from)),unname(tools::md5sum(to))))stop("Local frozen input differs: ",to)
    } else if(!file.copy(from,to))stop("Cannot copy frozen input: ",from)
  }
  data<-local
}
config<-eivGP::study_config(study,profile,output_root=output,data_root=data,cores=as.integer(get("cores","4")))
competitors<-strsplit(get("competitors",if(profile=="quick")"none"else"UC-GP,LVGP,EzGP"),",",fixed=TRUE)[[1]]
if(identical(competitors,"none"))competitors<-character()
eivGP::run_numerical_study(config,action,competitors=competitors)

}
if(sys.nframe()==0L)run_study_cli()
