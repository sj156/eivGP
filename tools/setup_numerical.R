#!/usr/bin/env Rscript
# Restore public versioned dependencies, then install this checkout's package.
args<-commandArgs(TRUE)
if(any(!args%in%c("--check","--install-local")))stop("Options: --check or --install-local (otherwise restore dependencies and install)")
script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
root<-dirname(dirname(normalizePath(script)))
lib<-Sys.getenv("EIVGP_R_LIBRARY",file.path(root,".r-library"))
dir.create(lib,recursive=TRUE,showWarnings=FALSE);.libPaths(c(normalizePath(lib),.libPaths()))
check<-"--check"%in%args
if(!check && !"--install-local"%in%args) {
  options(repos=c(CRAN="https://cloud.r-project.org"))
  if(!requireNamespace("renv",quietly=TRUE))install.packages("renv",lib=lib)
  renv::restore(project=root,library=lib,lockfile=file.path(root,"renv.lock"),prompt=FALSE)
}
if(!requireNamespace("jsonlite",quietly=TRUE))stop("Missing jsonlite; run setup without --check/--install-local")
lock<-jsonlite::read_json(file.path(root,"renv.lock"),simplifyVector=FALSE)
expected<-vapply(lock$Packages,`[[`,character(1),"Version")
actual<-vapply(names(expected),function(p)if(requireNamespace(p,quietly=TRUE))as.character(packageVersion(p))else "missing",character(1))
mismatch<-actual!=as.character(package_version(expected))
if(any(mismatch))stop("Dependency versions differ from renv.lock: ",paste(names(expected)[mismatch],collapse=", "),". Restore with tools/setup_numerical.R")
if(as.character(getRversion())!=lock$R$Version)warning("Lockfile was tested with R ",lock$R$Version,"; this R is ",getRversion(),". Record and validate differences before pooling results")
if(!check) {
  status<-system2(file.path(R.home("bin"),"R"),c("CMD","INSTALL",paste0("--library=",shQuote(lib)),shQuote(file.path(root,"eivGP"))))
  if(status!=0)stop("Local package installation failed")
}
source(file.path(root,"paper/numerical/common.R"));numerical_package_check()
cat("Numerical dependency versions and installed source identity match.\n")
capture.output(sessionInfo(),file=file.path(lib,"numerical-sessionInfo.txt"))
