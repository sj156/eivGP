#!/usr/bin/env Rscript
# Run from any directory. Use plan first; it never launches fitting.
# Rscript experiments/recover_publication.R plan --archive-root=/path/to/archive --data-root=/path/to/synthetic
# Rscript experiments/recover_publication.R run  --archive-root=/path/to/archive --data-root=/path/to/synthetic
args<-commandArgs(trailingOnly=TRUE)
action<-if(length(args)&&!startsWith(args[1],"--"))args[1]else"plan"
if(any(args%in%c("help","--help"))){
 cat("Usage: Rscript recover_publication.R [plan|run|check] [--archive-root=DIR] [--data-root=DIR] [--output=DIR] [--studies=study1,study2] [--cores=4] [--chain-workers=4] [--competitor-workers=4]\n",
 "archive-root contains study*-publication-* directories and competitor-reports/.\n",
 "data-root contains the original frozen study1/ and study2/ input directories.\n",
 "--cores=N caps both phases (default 4); worker defaults remain 4, even with a larger budget.\n",
 "Defaults: four parallel MCMC chains or four parallel competitor fits, in separate phases.\n",
 "Successful results are reused; run resumes missing tasks; check verifies common IDs and finite scores.\n")
 quit(status=0L)
}
if(!action%in%c("plan","run","check"))stop("Action must be plan, run, or check.")
cli<-grep("^--file=",commandArgs(FALSE),value=TRUE)
repo<-normalizePath(file.path(dirname(sub("^--file=","",cli[1])),".."),mustWork=TRUE)
options<-list()
for(a in args[startsWith(args,"--")]){
 if(!grepl("=",a,fixed=TRUE))stop("Options require --name=value: ",a)
 key<-sub("=.*$","",substring(a,3));value<-sub("^[^=]*=","",a)
 if(!key%in%c("archive-root","data-root","output","studies","cores","chain-workers","competitor-workers")||!is.null(options[[key]]))stop("Unknown or duplicated option: ",key)
 options[[key]]<-value
}
get_option<-function(key,default)if(is.null(options[[key]]))default else options[[key]]
archive<-get_option("archive-root",if(dir.exists("/Volumes/SANDISK USB"))"/Volumes/SANDISK USB"else file.path(repo,"reproduction"))
data<-get_option("data-root",file.path(repo,"reproduction","data","synthetic"))
output<-get_option("output",file.path(repo,"reproduction","recovery-publication"))
cores<-suppressWarnings(as.numeric(get_option("cores","4")))
workers<-suppressWarnings(as.numeric(get_option("chain-workers","4")))
competitor_workers<-suppressWarnings(as.numeric(get_option("competitor-workers","4")))
studies<-strsplit(get_option("studies","study1,study2"),",",fixed=TRUE)[[1]]
Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",VECLIB_MAXIMUM_THREADS="1")
source(file.path(repo,"codes","simulation_helpers.R"))
source(file.path(repo,"codes","publication_recovery.R"))
audit<-mixedgp_recover_publication(repo,archive,data,output,action,workers,studies,competitor_workers,core_budget=cores)
if(action!="plan"&&!all(audit$complete))quit(status=1L)
