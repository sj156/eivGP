# Run from the repository root after installing the current package.
if(!requireNamespace("eivGP",quietly=TRUE)||packageVersion("eivGP")<"0.4.0")
  stop("Install the current eivGP package first.")
root <- Sys.getenv("EIVGP_RESULTS_ROOT")
if(!nzchar(root))stop("Set EIVGP_RESULTS_ROOT to the external results folder.")
for(study in c("study1","study2"))
  rmarkdown::render(file.path("paper","numerical",paste0(study,".Rmd")),
    params=list(results_root=normalizePath(root),action="report"),
    envir=new.env(parent=globalenv()))
for(application in c("ocean","adni"))
  rmarkdown::render(file.path("paper","applications",paste0(application,".Rmd")),
    params=list(action="plan"),envir=new.env(parent=globalenv()))
