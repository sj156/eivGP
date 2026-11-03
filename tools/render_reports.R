# Resolve the checkout before rendering changes the notebook working directory.
script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
repo<-dirname(dirname(normalizePath(script,mustWork=TRUE)))
# Run from the repository root after installing the current package.
if(!requireNamespace("eivGP",quietly=TRUE)||packageVersion("eivGP")<"0.4.0")
  stop("Install the current eivGP package first.")
root <- Sys.getenv("EIVGP_RESULTS_ROOT")
if(!nzchar(root))stop("Set EIVGP_RESULTS_ROOT to the external results folder.")
root<-normalizePath(path.expand(root),mustWork=TRUE)
for(study in c("study1","study2"))
  rmarkdown::render(file.path(repo,"paper",paste0(study,".Rmd")),
    params=list(results_root=normalizePath(root),action="report"),
    envir=new.env(parent=globalenv()))
for(application in c("ocean","adni"))
  rmarkdown::render(file.path(repo,"paper",paste0(application,".Rmd")),
    params=list(action="plan"),envir=new.env(parent=globalenv()))
