test_that("paper workflows route both studies and reject misleading selections", {
  calls<-character()
  env<-new.env(parent=environment(run_paper_study))
  env$eivgp_numerical_api<-function(repo_root)list(
    numerical_default_output=function(study)file.path(tempdir(),study),
    numerical_default_data=function(study)file.path(tempdir(),"data",study),
    settings=identity,
    run_mechanism=function(cfg,action){calls<<-c(calls,paste0("mechanism:",action));cfg})
  env$mixedgp_study1_cc_workflow<-function(args,repo_root,root,data_root,action){calls<<-c(calls,paste0("complete-case:",action));list(root=root)}
  env$study_config<-function(...)list(...)
  env$run_numerical_study<-function(config,action,competitors){calls<<-c(calls,paste0("core:",action));config}
  run<-run_paper_study;environment(run)<-env
  x<-run("study2",list(action="run",output=tempdir(),"data-root"=tempdir()))
  expect_identical(calls,c("core:run","mechanism:run"))
  expect_named(x,c("core","mechanism"))
  expect_equal(x$mechanism$output,file.path(tempdir(),"mechanism"))
  calls<-character();run("study1",list(action="report",output=tempdir(),"data-root"=tempdir()))
  expect_identical(calls,c("core:report","complete-case:report"))
  calls<-character();run("study1",list(component="complete-case",profile="pilot",action="report",output=tempdir(),"data-root"=tempdir()))
  expect_identical(calls,"complete-case:report")
  calls<-character();run("study1",list(component="core",action="report",output=tempdir(),"data-root"=tempdir()))
  expect_identical(calls,"core:report")
  calls<-character();run("study2",list(component="mechanism",reps="1",action="report",output=tempdir(),"data-root"=tempdir()))
  expect_identical(calls,"mechanism:report")
  expect_error(run("study2",list(reps="1")),"require component=mechanism")
  expect_error(run("study1",list(component="mechanism")),"belongs to Study II")
  expect_error(run("study2",list(action="fit")),"Use action=run")
  calls<-character();expect_error(run("study1",list(action="fit")),"Use action=run")
  expect_length(calls,0)
  expect_error(run("study2",list(component="complete-case")),"Use mechanism")
  expect_error(run("study2",list(component="mechanism",competitors="none")),"Use methods")
})


test_that("mean and response integration budgets are separate and validated", {
  api<-eivgp_numerical_api(tempdir())
  a<-api$settings(list(profile="pilot"))
  expect_equal(a$gp_response_latent,1L)
  expect_equal(a$cc_response_latent,4L)
  b<-api$settings(list(profile="pilot",latent="128",`gp-response-latent`="4",`cc-response-latent`="16"))
  expect_equal(b$latent,128L)
  expect_equal(b$gp_response_latent,4L)
  expect_equal(b$cc_response_latent,16L)
  expect_error(api$settings(list(profile="pilot",latent="8",`gp-response-latent`="9")),"cannot exceed")
  expect_error(api$settings(list(profile="pilot",`cc-response-latent`="0")),"positive")
  expect_equal(api$settings(list(profile="quick",latent="2"))$cc_response_latent,2L)
})
