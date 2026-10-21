method_cli_env <- function() {
  root<-normalizePath(file.path(testthat::test_path(),'../../..'))
  env<-new.env(parent=globalenv())
  sys.source(file.path(root,'paper/applications/ocean.R'),env)
  sys.source(file.path(root,'paper/applications/ocean-workflow.R'),env)
  env
}

test_that('method entrypoints isolate fitting groups and saved-only never selects fitting', {
  e<-method_cli_env()
  for(group in c('eivGP','competitors','report')) {
    p<-e$ocean_method_params('--cores=16',group)
    s<-do.call(ocean_settings,p[intersect(names(p),names(formals(ocean_settings)))])
    tasks<-ocean_select_tasks(s)
    expect_identical(s$cores,16L)
    expect_identical(s$action,if(group=='report')'report' else 'fit')
    if(group=='eivGP')expect_true(all(tasks$module=='eiv'))
    if(group=='competitors') {
      expect_equal(sum(tasks$module=='original_baselines'),5)
      expect_equal(sum(tasks$module=='complete_case_baselines'),30)
      expect_false(any(tasks$module%in%c('eiv','continuation')))
    }
    expect_identical(e$ocean_method_params(c('--cores=16','--reuse-only=true'),group)$action,'report')
    expect_identical(e$ocean_method_params('--plan=true',group)$action,'plan')
  }
  expect_error(e$ocean_method_params('--action=fit','report'),'Choose the method')
  expect_error(e$ocean_method_params('--modules=eiv','competitors'),'Choose the method')
  expect_error(e$ocean_method_params('--reuse-only=maybe','eivGP'),'true or false')
})

test_that('saved-only preflight stops on missing fit outputs before execution', {
  e<-method_cli_env();out<-tempfile();dir.create(out)
  on.exit(unlink(out,recursive=TRUE))
  s<-ocean_settings(profile='quick',action='report',modules='eiv')
  expect_error(e$ocean_saved_method_check(out,out,s),'no fitting started')
})
