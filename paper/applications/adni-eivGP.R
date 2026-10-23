# Author: Xinye Chen. Shared package implementation, one allocated CPU budget.
if(sys.nframe()==0L) {
  script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
  source(file.path(dirname(normalizePath(script)),"adni-workflow.R"))
  adni_method_main("eivGP")
}
