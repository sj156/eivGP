# Author: Xinye Chen. Application entrypoint; implementation is in the installed package.
adni_load_api <- function(repo=NULL) {
  lib<-Sys.getenv('EIVGP_R_LIBRARY',if(is.null(repo))''else file.path(repo,'.r-library'))
  if(nzchar(lib)&&dir.exists(lib)).libPaths(c(normalizePath(lib),.libPaths()))
  if(!requireNamespace('eivGP',quietly=TRUE)||utils::packageVersion('eivGP')<'0.4.0.9014')
    stop('Install this checkout with Rscript tools/build_package.R and R CMD INSTALL eivGP.')
  eivGP::adni_workflow_api(repo)
}
adni_cli_params <- function(...) adni_load_api()$adni_cli_params(...)
adni_settings <- function(...) adni_load_api()$adni_settings(...)
adni_integer <- function(...) adni_load_api()$adni_integer(...)
adni_resolve_path <- function(...) adni_load_api()$adni_resolve_path(...)
adni_method_run <- function(repo,params,group) eivGP::run_adni_application(repo,params,group)
if(sys.nframe()==0L) {
  script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
  source(file.path(dirname(normalizePath(script)),'adni-workflow.R'))
  adni_method_main('report')
}
