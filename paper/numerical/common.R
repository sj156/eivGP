# Repository-relative bootstrap: works from any working directory on macOS/Linux.
numerical_repo_root <- function() {
  candidates<-c(Sys.getenv("EIVGP_REPO_ROOT"),getwd())
  args<-grep("^--file=",commandArgs(FALSE),value=TRUE)
  if(length(args))candidates<-c(dirname(normalizePath(sub("^--file=","",args[1]),mustWork=FALSE)),candidates)
  for(fr in sys.frames())if(!is.null(fr$ofile))candidates<-c(dirname(normalizePath(fr$ofile,mustWork=FALSE)),candidates)
  for(p in unique(candidates[nzchar(candidates)]))for(i in 1:8) {
    if(file.exists(file.path(p,"eivGP","DESCRIPTION")) && dir.exists(file.path(p,"paper","numerical")))return(normalizePath(p))
    parent<-dirname(p);if(identical(parent,p))break;p<-parent
  }
  stop("Cannot locate checkout; set EIVGP_REPO_ROOT to the repository root")
}
NUMERICAL_REPO_ROOT<-numerical_repo_root()
lib<-Sys.getenv("EIVGP_R_LIBRARY",file.path(NUMERICAL_REPO_ROOT,".r-library"))
if(dir.exists(lib)).libPaths(c(normalizePath(lib),.libPaths()))
numerical_default_data <- function(study)file.path(NUMERICAL_REPO_ROOT,"results","data",study)
numerical_default_output <- function(name)file.path(NUMERICAL_REPO_ROOT,"results",name)
numerical_package_check <- function() {
  if(!requireNamespace("eivGP",quietly=TRUE)||packageVersion("eivGP")<"0.4.0")stop("Run Rscript tools/setup_numerical.R from the checkout first")
  if(!requireNamespace("jsonlite",quietly=TRUE))stop("Restore the numerical dependencies first")
  lock<-jsonlite::read_json(file.path(NUMERICAL_REPO_ROOT,"renv.lock"))
  expected<-vapply(lock$Packages,`[[`,character(1),"Version")
  actual<-vapply(names(expected),function(p)if(requireNamespace(p,quietly=TRUE))as.character(packageVersion(p))else "missing",character(1))
  mismatch<-actual!=as.character(package_version(expected))
  if(any(mismatch))stop("Dependencies differ from renv.lock: ",paste(names(expected)[mismatch],collapse=", "),"; run tools/setup_numerical.R")
  installed<-getFromNamespace("eivgp_module_identity","eivGP")()
  paths<-file.path(NUMERICAL_REPO_ROOT,"eivGP","R",names(installed))
  if(!identical(unname(tools::md5sum(paths)),unname(installed)))
    stop("Installed eivGP differs from this checkout; rerun tools/setup_numerical.R")
}
