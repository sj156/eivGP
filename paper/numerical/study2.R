# Study-specific choices; all implementation and reporting come from eivGP.
study2_analysis <- function(results_root, action=c("report","plan","data","fit","run"),
                              profile="paper", cores=4L, tables_dir=NULL, competitor_dir=NULL) {
  action<-match.arg(action)
  if(!requireNamespace("eivGP",quietly=TRUE) || packageVersion("eivGP") < "0.4.0")
    stop("Install the repository eivGP package first.")
  if(action=="report") {
    source<-file.path(results_root,"study2/runs/recovery-publication-linux/study2")
    if(!is.null(tables_dir) && nzchar(tables_dir)) source<-tables_dir
    if(!file.exists(file.path(source,"predictive_metrics.csv"))) stop("Selected paper export not found: ",source)
    results<-eivGP::read_study_results(source,competitor_dir)
    return(list(source=source,summary=eivGP::summarize_comparison(results),diagnostics=results$diagnostics))
  }
  config<-eivGP::study_config("study2",profile,
    output_root=file.path(results_root,"study2/runs"),
    data_root=file.path(results_root,"data/synthetic/study2"),cores=cores)
  eivGP::run_numerical_study(config,action,
    competitors=if(profile=="quick")character() else c("UC-GP","LVGP","EzGP"))
}
if(sys.nframe()==0L) {
  args<-commandArgs(trailingOnly=TRUE)
  if(length(args)&&all(grepl("^--",args))) {
    script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
    env<-new.env(parent=globalenv())
    sys.source(file.path(dirname(normalizePath(script)),"run-study.R"),envir=env)
    env$run_study_cli(c("--study=study2",args))
  } else {
    root<-Sys.getenv("EIVGP_RESULTS_ROOT")
    if(!nzchar(root))stop("Set EIVGP_RESULTS_ROOT, or use --action=plan --data-root=... --output=...")
    print(study2_analysis(root,action=if(length(args))args[1] else "report"))
  }
}
