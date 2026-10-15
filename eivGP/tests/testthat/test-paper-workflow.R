test_that("paper workflows route both studies and reject misleading selections", {
  calls<-character()
  env<-new.env(parent=environment(run_paper_study))
  env$eivgp_numerical_api<-function(repo_root)list(
    numerical_default_output=function(study)file.path(tempdir(),study),
    numerical_default_data=function(study)file.path(tempdir(),"data",study),
    settings=identity,
    run_mechanism=function(cfg,action){calls<<-c(calls,paste0("mechanism:",action));cfg})
  env$study_config<-function(...)list(...)
  env$run_numerical_study<-function(config,action,competitors){calls<<-c(calls,paste0("core:",action));config}
  run<-run_paper_study;environment(run)<-env
  x<-run("study2",list(action="run",output=tempdir(),"data-root"=tempdir()))
  expect_identical(calls,c("core:run","mechanism:run"))
  expect_named(x,c("core","mechanism"))
  expect_equal(x$mechanism$output,file.path(tempdir(),"mechanism"))
  calls<-character();run("study1",list(action="report",output=tempdir(),"data-root"=tempdir()))
  expect_identical(calls,"core:report")
  calls<-character();run("study2",list(component="mechanism",reps="1",action="report",output=tempdir(),"data-root"=tempdir()))
  expect_identical(calls,"mechanism:report")
  expect_error(run("study2",list(reps="1")),"require component=mechanism")
  expect_error(run("study1",list(component="mechanism")),"belongs to Study II")
  expect_error(run("study2",list(action="fit")),"Use action=run")
  expect_error(run("study2",list(component="mechanism",competitors="none")),"Use methods")
})
