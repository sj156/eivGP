# eivGP: shared package implementation, one allocated CPU budget.
if(sys.nframe()==0L) {
  script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
  source(file.path(dirname(normalizePath(script)),"ocean-tools.R"))
  ocean_method_main("eivGP")
}
