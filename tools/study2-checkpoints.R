#!/usr/bin/env Rscript
# Explicit maintenance of a stopped Study II reference directory; never fits.
script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
if(is.na(script))stop("Run with Rscript")
root<-dirname(dirname(normalizePath(script)))
source(file.path(root,"paper/numerical/common.R"))
numerical_package_check()
a<-list()
for(x in commandArgs(TRUE)) {
  if(!grepl("^--[^=]+=.+$",x))stop("Use --name=value")
  k<-sub("^--([^=]+)=.*$","\\1",x)
  if(k%in%names(a))stop("Repeated option: ",k)
  a[[k]]<-sub("^--[^=]+=","",x)
}
if(length(setdiff(names(a),c("output","audit-dir","dry-run","compact-measurement"))))stop("Unknown cleanup option")
if(is.null(a$output)||is.null(a[["audit-dir"]]))stop("Supply --output and --audit-dir")
logical_option<-function(k,default) {
  x<-if(is.null(a[[k]]))default else tolower(a[[k]])
  if(length(x)!=1L||!x%in%c("true","false"))stop(k," must be true or false")
  identical(x,"true")
}
print(eivGP::prune_study_checkpoints(a$output,
  dry_run=logical_option("dry-run","true"),audit_dir=a[["audit-dir"]],
  compact_measurement=logical_option("compact-measurement","false")))
