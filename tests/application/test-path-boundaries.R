extract_path_function <- function(path,name) {
  find<-function(x) {
    if(missing(x))return(NULL)
    if(is.call(x)&&identical(x[[1L]],as.name('<-'))&&identical(x[[2L]],as.name(name)))return(x)
    if(is.call(x)||is.expression(x)||is.pairlist(x))for(a in as.list(x)) {
      hit<-find(a);if(!is.null(hit))return(hit)
    }
    NULL
  }
  assignment<-find(parse(path));stopifnot(!is.null(assignment))
  env<-new.env(parent=asNamespace('eivGP'));eval(assignment,env);env[[name]]
}
test_that('path normalization is stable for existing and missing destinations', {
  absolute<-getFromNamespace('mixedgp_absolute_path','eivGP')
  root<-tempfile('path-boundaries-');dir.create(root);on.exit(unlink(root,recursive=TRUE))
  dir.create(file.path(root,'checkout'));dir.create(file.path(root,'caller'))
  repo<-normalizePath(file.path(root,'checkout'));caller<-normalizePath(file.path(root,'caller'))
  old<-getwd();on.exit(setwd(old),add=TRUE);setwd(caller)
  expect_identical(absolute('new/results'),file.path(caller,'new/results'))
  expect_identical(absolute('new/../results'),file.path(caller,'results'))
  expect_identical(absolute('./new path/results',base=repo),file.path(repo,'new path/results'))
  expect_identical(absolute(file.path(caller,'new'),base=repo),file.path(caller,'new'))
  expect_identical(absolute('~'),normalizePath(path.expand('~')))
  expect_false(dir.exists(file.path(caller,'new')))
  expect_error(absolute(''),'nonempty');expect_error(absolute(NA_character_),'nonempty')
  for(study in c('study1','study2')) {
    cfg<-eivGP::study_config(study,'quick','new/results','new/data',cores=1L)
    expect_identical(cfg$data_root,file.path(caller,'new/data'))
    setwd(repo);expect_identical(cfg$data_root,file.path(caller,'new/data'));setwd(caller)
  }
  api<-getFromNamespace('eivgp_numerical_api','eivGP')(repo)
  settings<-api$settings(list(profile='quick',output='new/out',`data-root`='new/data',`full-fit-map`='new/map.csv'))
  expect_identical(settings$output,file.path(caller,'new/out'))
  expect_identical(settings$data_root,file.path(caller,'new/data'))
  expect_identical(settings$full_fit_map,file.path(caller,'new/map.csv'))
})
test_that('all application resolvers anchor new relative paths to the checkout', {
  repo<-Sys.getenv('EIVGP_REPO_ROOT',normalizePath(testthat::test_path('../..')))
  root<-tempfile('app-paths-');dir.create(root);on.exit(unlink(root,recursive=TRUE))
  for(name in c('ocean','adni','oai')) {
    resolve<-extract_path_function(file.path(repo,'tools',paste0(name,'-tools.R')),paste0(name,'_resolve_path'))
    expect_identical(resolve('missing/../output with spaces',root),file.path(normalizePath(root),'output with spaces'))
    expect_identical(resolve(file.path(normalizePath(root),'missing/out'),repo),file.path(normalizePath(root),'missing/out'))
    expect_false(dir.exists(file.path(root,'missing')))
  }
})
test_that('publication planning uses absolute roots without preparing or fitting', {
  root<-tempfile('publication-paths-');dir.create(root);on.exit(unlink(root,recursive=TRUE))
  old<-getwd();setwd(root);on.exit(setwd(old),add=TRUE)
  f<-getFromNamespace('mixedgp_publication_workflow','eivGP')
  result<-suppressMessages(capture.output(plan<-f('study2',list(action='plan',stage='full',reps='1',methods='CUT',output='new/out',`data-root`='new/data'),root)))
  expect_true(startsWith(plan$commands[[1]]$output,normalizePath(root)))
  expect_identical(plan$commands[[1]][['data-root']],file.path(normalizePath(root),'new/data'))
  expect_false(dir.exists(file.path(root,'new')))
})
