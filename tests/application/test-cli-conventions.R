test_that('study-style CLI actions work across all application lanes', {
f<-getFromNamespace('mixedgp_workflow_cli','eivGP')
stopifnot(identical(f(c('--action=plan','--cores=12')),c('--cores=12','--plan=true')),
 identical(f('--action=report'), '--reuse-only=true'),length(f('--action=run'))==0,
 identical(f('--plan=true'),'--plan=true'))
for(a in list(c('--action=run','--plan=true'),c('--action=run','--reuse-only=false'),c('--action=plan','--action=run'),'--action=continue'))stopifnot(inherits(try(f(a),silent=TRUE),'try-error'))
stopifnot(inherits(try(f('--action=run','report'),silent=TRUE),'try-error'))
repo<-normalizePath(testthat::test_path('../..'))
for(name in c('adni','oai','ocean')) {
 env<-new.env(parent=asNamespace('eivGP'));sys.source(file.path(repo,'tools',paste0(name,'-tools.R')),env)
 args<-c('--action=plan','--cores=12','--data-root=/tmp/absent-input','--output=/tmp/example-output')
 if(name=='ocean') {
  api<-env$ocean_workflow_api(repo);p<-api$ocean_method_params(args,'eivGP')
  stopifnot(p$action=='plan',p$data_dir=='/tmp/absent-input')
  stopifnot(api$ocean_method_params('--action=report','eivGP')$action=='report')
  stopifnot(inherits(try(api$ocean_method_params('--action=run','report'),silent=TRUE),'try-error'))
 } else {
  api<-get(paste0(name,'_load_api'),env)(repo)
  parse<-api[[paste0(name,'_cli_params')]]
  p<-parse(args);stopifnot(isTRUE(p$plan),p$data_root=='/tmp/absent-input',p$output=='/tmp/example-output')
  stopifnot(isTRUE(parse('--action=report')$reuse_only),!isTRUE(parse('--action=run')$plan))
 }
}
expect_identical(f('--action=report','eivGP'),'--reuse-only=true')
})
