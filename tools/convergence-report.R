# Saved-artifact reporting only. Does not use fitting/provenance resume gates.
# Manifest CSV columns: design,root. Equal design labels assert matching protocols.
convergence_report <- function(manifest,output) {
  absolute<-getFromNamespace("mixedgp_absolute_path","eivGP")
  manifest<-absolute(manifest);output<-absolute(output)
  sources<-read.csv(manifest,stringsAsFactors=FALSE)
  stopifnot(all(c("design","root")%in%names(sources)))
  sources$root<-vapply(sources$root,absolute,character(1))
  eivGP::report_numerical_convergence(setNames(sources$root,sources$design),output)
}
if(sys.nframe()==0L) {
  args<-commandArgs(TRUE)
  if(length(args)!=2L)stop("Usage: Rscript convergence-report.R manifest.csv separate-output-directory")
  lib<-Sys.getenv("EIVGP_R_LIBRARY");if(nzchar(lib)).libPaths(c(lib,.libPaths()))
  convergence_report(args[1],args[2])
}
