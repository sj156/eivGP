# Manuscript completion from existing studies: --component=follow-up (defaults to all 50 replications; never refits core).
# Direct follow-up: --action=run --profile=paper --cores=4 --component=publication --stage=pilot
# Publication completion: --component=publication --stage=pilot|integration|full
# Default action is plan; use --action=data, run, or report explicitly.
# Original (X,C) predictors and complete-case (X,U,C) additions share this entry point.
# Study-specific choices; all implementation and reporting come from eivGP.
study1_analysis <- function(results_root, action=c("report","plan","data","fit","run"),
                              profile="paper", cores=4L, tables_dir=NULL, competitor_dir=NULL,
                              mixing_plots=TRUE) {
  action<-match.arg(action)
  if(!requireNamespace("eivGP",quietly=TRUE) || packageVersion("eivGP") < "0.4.0")
    stop("Install the repository eivGP package first.")
  if(action=="report") {
    source<-file.path(results_root,"study1/runs/recovery-publication-mini/study1")
    if(!is.null(tables_dir) && nzchar(tables_dir)) source<-tables_dir
    if(!file.exists(file.path(source,"predictive_metrics.csv"))) stop("Selected paper export not found: ",source)
    results<-eivGP::read_study_results(source,competitor_dir)
    plots<-if(isTRUE(mixing_plots)) eivGP::save_mixing_plots(source) else character()
    return(list(source=source,mixing_plots=plots,summary=eivGP::summarize_comparison(results),diagnostics=results$diagnostics))
  }
  eivGP::run_paper_study("study1",list(action=action,profile=profile,
    `mixing-plots`=mixing_plots,cores=as.character(cores),output=file.path(results_root,"study1")),
    repo_root=Sys.getenv("EIVGP_REPO_ROOT",getwd()))
}
# Notebook and batch commands share the same package orchestration.
study1_workflow <- function(args=list(), repo_root=Sys.getenv("EIVGP_REPO_ROOT",getwd())) {
  env<-new.env(parent=globalenv())
  sys.source(file.path(repo_root,"paper/numerical/common.R"),envir=env)
  env$numerical_package_check()
  eivGP::run_paper_study("study1",args,repo_root)
}
if(sys.nframe()==0L) {
  args<-commandArgs(trailingOnly=TRUE)
  if(!length(args)||all(grepl("^--",args))) {
    script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
    env<-new.env(parent=globalenv())
    sys.source(file.path(dirname(normalizePath(script)),"common.R"),envir=env)
    env$run_study_cli("study1",args)
  } else {
    root<-Sys.getenv("EIVGP_RESULTS_ROOT")
    if(!nzchar(root))stop("Set EIVGP_RESULTS_ROOT, or use --action=plan --data-root=... --output=...")
    print(study1_analysis(root,action=if(length(args))args[1] else "report"))
  }
}
