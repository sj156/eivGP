test_that("Ocean command-line selections preserve the package workflow controls", {
  root<-normalizePath(testthat::test_path("../../.."),mustWork=TRUE)
  testthat::skip_if_not(file.exists(file.path(root,"paper/applications/ocean.R")))
  env<-new.env(parent=globalenv())
  sys.source(file.path(root,"paper/applications/ocean.R"),env)
  p<-env$ocean_cli_params(c("--action=run","--profile=pilot","--dataset=ocean500",
    "--cores=4","--folds=1,3:5","--cells=O20,O50","--output=/external/run",
    "--diagnostics-workers=2","--modules=eiv,prediction"))
  settings<-do.call(ocean_settings,p[intersect(names(p),names(formals(ocean_settings)))])
  expect_identical(settings$folds,c(1L,3L,4L,5L))
  expect_identical(settings$cells,c("O20","O50"))
  expect_identical(settings$modules,c("eiv","prediction"))
  expect_identical(settings$action,"run")
  expect_identical(c(settings$iter,settings$burn,settings$chains),c(80L,30L,4L))
  expect_identical(p$output_dir,"/external/run")
  expect_identical(p$diagnostics_workers,2L)
  expect_error(env$ocean_cli_params("--folds=0:5"),"folds")
  expect_error(env$ocean_cli_params("--cores=1.5"),"positive integer")
  expect_error(env$ocean_cli_params(c("--output=a","--output-dir=b")),"Repeated")
  expect_error(env$ocean_cli_params("--unknown=1"),"Unknown")
  expect_error(env$ocean_cli_params("--action"),"name=value")
  path<-tempfile(fileext=".rds");on.exit(unlink(path),add=TRUE)
  saveRDS(p,path)
  expect_identical(env$ocean_cli_params(path),p)
})

test_that("unspecified Ocean result directories keep their workflow defaults", {
  root<-normalizePath(testthat::test_path("../../.."),mustWork=TRUE)
  testthat::skip_if_not(file.exists(file.path(root,"paper/applications/ocean.R")))
  env<-new.env(parent=globalenv())
  sys.source(file.path(root,"paper/applications/ocean.R"),env)
  env$ocean_load_api<-function(repository_dir,envir)list2env(list(
    ocean_settings=getFromNamespace("ocean_settings","eivGP"),
    ocean_resolve_path=getFromNamespace("ocean_resolve_path","eivGP"),
    ocean_logged=function(path,expr)force(expr)),envir=envir)
  called<-NULL
  env$ocean_analysis<-function(...) {called<<-list(...);list()}
  invisible(capture.output(env$ocean_run_cli(root,list(action="report",
    import_dirs=c("relative-import","")))))
  expect_identical(called$cc_output_dir,"")
  expect_identical(called$diagnostics_dir,"")
  expect_identical(called$import_dirs,c(file.path(root,"relative-import"),""))
  expect_identical(called[[3]],file.path(root,"results/ocean/ocean500/paper"))
})
