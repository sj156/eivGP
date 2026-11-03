# Saved historical paper export only. No fitting or data generation.
study2_analysis <- function(results_root, action="report",
                              profile="paper", cores=4L, tables_dir=NULL, competitor_dir=NULL,
                              mixing_plots=FALSE, repo_root=Sys.getenv("EIVGP_REPO_ROOT",getwd())) {
  if(!identical(action,"report"))stop("This script only reports saved results; use study2-eivGP.R or study2-competitors.R to plan/run jobs")
  if(!requireNamespace("eivGP",quietly=TRUE) || packageVersion("eivGP") < "0.4.0")
    stop("Install the repository eivGP package first.")
  if(action=="report") {
    if(!is.null(tables_dir) && nzchar(tables_dir)) source<-tables_dir else {
      manifest<-file.path(repo_root,"paper/manifests/study2/manifest.json")
      if(!file.exists(manifest))stop("Set EIVGP_REPO_ROOT to the checkout for the paper result manifest")
      selection<-jsonlite::read_json(manifest)$current_primary_export$raw_tables
      source<-file.path(results_root,selection)
    }
    if(!file.exists(file.path(source,"predictive_metrics.csv"))) stop("Selected paper export not found: ",source)
    results<-eivGP::read_study_results(source,competitor_dir)
    # Combined CSV exports contain run-level screens, not chain-resolved traces.
    # Historical reporting is CSV-only by default; lane reports own trace plots.
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
  script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
  repo<-normalizePath(file.path(dirname(script),".."))
  print(study2_analysis(root,tables_dir=if(nzchar(tables))tables else NULL,repo_root=repo))
}
