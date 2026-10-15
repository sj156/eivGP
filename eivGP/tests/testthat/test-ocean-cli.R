test_that("Ocean command-line selections preserve the package workflow controls", {
  root<-normalizePath(testthat::test_path("../../.."),mustWork=TRUE)
  testthat::skip_if_not(file.exists(file.path(root,"paper/applications/ocean.R")))
  env<-new.env(parent=globalenv())
  sys.source(file.path(root,"paper/applications/ocean.R"),env)
  p<-env$ocean_cli_params(c("--action=run","--profile=pilot","--dataset=ocean500",
    "--cores=4","--folds=1,3:5","--cells=O20,O50","--output=/external/run",
    "--modules=eiv,prediction"))
  settings<-do.call(ocean_settings,p[intersect(names(p),names(formals(ocean_settings)))])
  expect_identical(settings$folds,c(1L,3L,4L,5L))
  expect_identical(settings$cells,c("O20","O50"))
  expect_identical(settings$modules,c("eiv","prediction"))
  expect_identical(settings$action,"fit")
  expect_identical(c(settings$iter,settings$burn,settings$chains),c(80L,30L,4L))
  expect_identical(p$output_dir,"/external/run")
  expect_error(env$ocean_cli_params("--diagnostics-workers=2"),"Unknown")
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
  env$ocean_check_report_dependencies<-function()invisible(TRUE)
  env$ocean_render_report<-function(...)"ocean.html"
  called<-NULL
  env$ocean_analysis<-function(...) {called<<-list(...);list()}
  invisible(capture.output(env$ocean_run_cli(root,list(action="report",
    import_dirs=c("relative-import","")))))
  expect_identical(called$cc_output_dir,"")
  expect_identical(called$diagnostics_dir,"")
  expect_identical(called$import_dirs,c(file.path(root,"relative-import"),""))
  expect_identical(called[[3]],file.path(path.expand("~/eivGP-ocean"),"ocean500","paper"))
})


test_that("Ocean CLI accepts an integer extension length", {
  root<-normalizePath(testthat::test_path("../../.."),mustWork=TRUE)
  testthat::skip_if_not(file.exists(file.path(root,"paper/applications/ocean.R")))
  env<-new.env(parent=globalenv())
  sys.source(file.path(root,"paper/applications/ocean.R"),env)
  p<-env$ocean_cli_params(c("--action=continue","--continuation-iter=12000"))
  expect_identical(do.call(ocean_settings,p)$continuation_iter,12000L)
  for(value in c("0","-1","1.5","Inf","NA","2147483648"))
    expect_error(env$ocean_cli_params(paste0("--continuation-iter=",value)),"positive integer")
})


test_that("one-click fitting and continuation analyze serially before rendering", {
  root<-normalizePath(testthat::test_path("../../.."),mustWork=TRUE)
  env<-new.env(parent=globalenv());sys.source(file.path(root,"paper/applications/ocean.R"),env)
  for(action in c("fit","continue","analyze","report","plan")) {
    stages<-character();budgets<-integer();imports<-logical();chains<-integer()
    execute<-function(settings,import) {
      stages<<-c(stages,settings$action);budgets<<-c(budgets,settings$cores)
      imports<<-c(imports,import);chains<<-c(chains,settings$chains);list()
    }
    render<-function() {stages<<-c(stages,"render");"report.html"}
    settings<-ocean_settings(profile="pilot",action=action,cores=12,chains=8)
    result<-env$ocean_cli_workflow(settings,execute,render)
    expected<-if(action%in%c("fit","continue"))c(action,"analyze","render") else
      if(action=="plan")"plan" else c(action,"render")
    expect_identical(stages,expected)
    if(action%in%c("fit","continue")) {
      expect_identical(budgets,c(12L,1L));expect_identical(imports,c(TRUE,FALSE))
      expect_identical(chains,c(8L,8L))
    }
  }
  fitted<-FALSE;rendered<-FALSE
  expect_error(env$ocean_cli_workflow(ocean_settings(action="fit"),
    function(s,import) {if(s$action=="analyze")stop("analysis failed");fitted<<-TRUE;list()},
    function() {rendered<<-TRUE;"report.html"}),"analysis failed")
  expect_true(fitted);expect_false(rendered)
})


test_that("paper sampling settings are fixed and CLI exposes profiles instead of sampling knobs", {
  expect_identical(c(ocean_settings()$chains,ocean_settings()$iter,ocean_settings()$burn),c(4L,50000L,5000L))
  for(action in c("plan","fit","continue")) {
    expect_error(ocean_settings(action=action,iter=100000),"paper fixes")
    expect_error(ocean_settings(action=action,chains=8),"paper fixes")
    expect_error(ocean_settings(action=action,burn=1000),"paper fixes")
  }
  root<-normalizePath(testthat::test_path("../../.."),mustWork=TRUE)
  env<-new.env(parent=globalenv());sys.source(file.path(root,"paper/applications/ocean.R"),env)
  for(flag in c("--chains=4","--iter=50000","--burn=5000"))
    expect_error(env$ocean_cli_params(flag),"fixed by --profile")
  output<-tempfile();dir.create(output);on.exit(unlink(output,recursive=TRUE),add=TRUE)
  saveRDS(list(numerics=list(iter=80L,burn=30L,chains=2L)),file.path(output,"workflow_identity.rds"))
  p<-env$ocean_restore_fit_settings(list(action="fit",profile="paper"),output)
  expect_null(p$iter);expect_null(p$chains)
  expect_error(env$ocean_restore_fit_settings(list(action="continue",profile="paper"),output),"fixed paper")
  for(action in c("analyze","report")) {
    p<-env$ocean_restore_fit_settings(list(action=action,profile="paper"),output)
    expect_identical(c(p$chains,p$iter,p$burn),c(2L,80L,30L))
  }
  p<-env$ocean_restore_fit_settings(list(action="continue",profile="pilot",continuation_iter=7L),output)
  expect_identical(c(p$chains,p$iter,p$burn),c(2L,80L,30L));expect_identical(p$continuation_iter,7L)
})


test_that("Ocean CLI and Rmd defaults share an external home-directory root", {
  root<-normalizePath(testthat::test_path("../../.."),mustWork=TRUE)
  env<-new.env(parent=globalenv());sys.source(file.path(root,"paper/applications/ocean.R"),env)
  env$ocean_resolve_path<-getFromNamespace("ocean_resolve_path","eivGP")
  for(dataset in c("ocean200","ocean500","ocean1000")) {
    expected<-file.path(path.expand("~/eivGP-ocean"),dataset,"paper")
    expect_identical(env$ocean_output_dir(root,dataset,"paper"),expected)
    lines<-readLines(file.path(root,"paper/config",paste0(dataset,".yml")))
    configured<-strsplit(lines[startsWith(lines,"output_dir:")],'"',fixed=TRUE)[[1]][2]
    expect_identical(env$ocean_output_dir(root,dataset,"paper",configured),expected)
  }
  expect_identical(env$ocean_output_dir(root,"ocean500","pilot"),
    file.path(path.expand("~/eivGP-ocean"),"ocean500","pilot"))
  expect_identical(env$ocean_output_dir(root,"ocean500","paper","~/custom-ocean"),path.expand("~/custom-ocean"))
  expect_identical(env$ocean_output_dir(root,"ocean500","paper","relative-run"),file.path(root,"relative-run"))
  expect_identical(env$ocean_output_dir(root,"ocean500","paper","/external/run"),"/external/run")
})
