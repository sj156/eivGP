# Application settings only; fitting, comparison, and reporting are package APIs.
adni_analysis <- function(data_dir,output_dir,action=c("plan","report","run"),profile="paper",cores=4L) {
  if(!requireNamespace("eivGP",quietly=TRUE) || packageVersion("eivGP") < "0.4.0")
    stop("Install the repository eivGP package first.")
  eivGP::run_adni_analysis(data_dir,output_dir,match.arg(action),profile,cores)
}
if(sys.nframe()==0L) {
  args<-commandArgs(trailingOnly=TRUE)
  if(length(args)<2L)stop("Usage: Rscript adni.R DATA_DIR OUTPUT_DIR [plan|report|run]")
  print(adni_analysis(args[1],args[2],if(length(args)>2)args[3] else "plan"))
}
