# Combined compatibility entry point; separate-machine scripts: study1-eivGP.R and study1-competitors.R.
# Complete Study I paper reproduction: original study and calibrated competitors.
# Default action is plan; use --action=data, run, or report explicitly.
# Original (X,C) predictors and complete-case (X,U,C) additions share this entry point.
# Study-specific choices; all implementation and reporting come from eivGP.
study1_analysis <- function(results_root, action=c("report","plan","data","fit","run"),
                              profile="paper", cores=12L, tables_dir=NULL, competitor_dir=NULL,
                              mixing_plots=TRUE, repo_root=Sys.getenv("EIVGP_REPO_ROOT",getwd())) {
  action<-match.arg(action)
  if(!requireNamespace("eivGP",quietly=TRUE) || packageVersion("eivGP") < "0.4.0")
    stop("Install the repository eivGP package first.")
  if(action=="report") {
    if(!is.null(tables_dir) && nzchar(tables_dir)) source<-tables_dir else {
      manifest<-file.path(repo_root,"paper/manifests/study1/manifest.json")
      if(!file.exists(manifest))stop("Set EIVGP_REPO_ROOT to the checkout for the paper result manifest")
      source<-file.path(results_root,jsonlite::read_json(manifest)$current_primary_export$raw_tables)
    }
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
  sys.source(file.path(repo_root,"tools/numerical-common.R"),envir=env)
  env$numerical_package_check()
  if(is.null(args$cores))args$cores<-"12"
  eivGP::run_paper_study("study1",args,repo_root)
}
if(sys.nframe()==0L) {
  args<-commandArgs(trailingOnly=TRUE)
  if(!length(args)||all(grepl("^--",args))) {
    script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
    env<-new.env(parent=globalenv())
    sys.source(file.path(dirname(normalizePath(script)),"numerical-common.R"),envir=env)
    env$run_study_cli("study1",args)
  } else {
    root<-Sys.getenv("EIVGP_RESULTS_ROOT")
    if(!nzchar(root))stop("Set EIVGP_RESULTS_ROOT, or use --action=plan --data-root=... --output=...")
    script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
    print(study1_analysis(root,action=if(length(args))args[1] else "report",
      repo_root=dirname(dirname(normalizePath(script)))))
  }
}
