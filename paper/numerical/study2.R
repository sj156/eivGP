# Saved historical paper export only. No fitting or data generation.
study2_analysis <- function(results_root, action="report",
                              profile="paper", cores=4L, tables_dir=NULL, competitor_dir=NULL,
                              mixing_plots=TRUE) {
  if(!identical(action,"report"))stop("This script only reports saved results; use study2-eivGP.R or study2-competitors.R to plan/run jobs")
  if(!requireNamespace("eivGP",quietly=TRUE) || packageVersion("eivGP") < "0.4.0")
    stop("Install the repository eivGP package first.")
  if(action=="report") {
    source<-file.path(results_root,"study2/runs/recovery-publication-linux/study2")
    if(!is.null(tables_dir) && nzchar(tables_dir)) source<-tables_dir
    if(!file.exists(file.path(source,"predictive_metrics.csv"))) stop("Selected paper export not found: ",source)
    results<-eivGP::read_study_results(source,competitor_dir)
    plots<-if(isTRUE(mixing_plots)) eivGP::save_mixing_plots(source) else character()
    return(list(source=source,mixing_plots=plots,summary=eivGP::summarize_comparison(results),diagnostics=results$diagnostics))
  }
}
if(sys.nframe()==0L) {
  args<-commandArgs(TRUE)
  if(length(args))stop("Set EIVGP_RESULTS_ROOT and optionally EIVGP_STUDY2_TABLES; no fitting options here")
  root<-Sys.getenv("EIVGP_RESULTS_ROOT")
  if(!nzchar(root))stop("Set EIVGP_RESULTS_ROOT")
  tables<-Sys.getenv("EIVGP_STUDY2_TABLES")
  print(study2_analysis(root,tables_dir=if(nzchar(tables))tables else NULL))
}
