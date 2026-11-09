# Run with EIVGP_TEST_REPO pointing to this checkout and its package installed.
repo<-Sys.getenv("EIVGP_TEST_REPO");stopifnot(nzchar(repo))
e<-new.env(parent=globalenv());sys.source(file.path(repo,"tools/study2-common.R"),e)
testthat::test_that("resume copies verified source artifacts once and routes only CC by default", {
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  source<-file.path(root,"original");dir.create(source)
  saveRDS(list(value=1),file.path(source,"fit.rds"))
  before<-tools::md5sum(file.path(source,"fit.rds"))
  output<-file.path(root,"resumed");calls<-list()
  runner<-function(lane,args,repo_root) {calls[[length(calls)+1L]]<<-list(lane=lane,args=args)}
  a<-list(`source-root`=source,output=output,action="plan")
  e$study2_resume_competitors(a,repo,runner)
  testthat::expect_false(dir.exists(output))
  a$action<-"run";e$study2_resume_competitors(a,repo,runner)
  testthat::expect_identical(tools::md5sum(names(before)),before)
  testthat::expect_identical(readRDS(file.path(output,"fit.rds")),list(value=1))
  testthat::expect_identical(tail(calls,1)[[1]]$args$methods,"UC-GP-CC,LVGP-CC,EzGP-CC")
  testthat::expect_identical(tail(calls,1)[[1]]$args$cores,"12")
  testthat::expect_null(tail(calls,1)[[1]]$args[["source-root"]])
  # Work done after the copy must survive every subsequent resume.
  saveRDS(list(value=2),file.path(output,"fit.rds"))
  e$study2_resume_competitors(a,repo,runner)
  testthat::expect_identical(readRDS(file.path(output,"fit.rds")),list(value=2))
  testthat::expect_identical(tools::md5sum(names(before)),before)
  bad<-a;bad$output<-source
  testthat::expect_error(e$study2_resume_competitors(bad,repo,runner),"separate")
  # A source with an active owner is not copied.
  lock<-file.path(source,"RUN.lock");dir.create(lock)
  saveRDS(list(pid=Sys.getpid(),host=unname(Sys.info()[["nodename"]])),file.path(lock,"owner.rds"))
  a$output<-file.path(root,"other")
  testthat::expect_error(e$study2_resume_competitors(a,repo,runner),"active")
  testthat::expect_false(dir.exists(a$output))
})

testthat::test_that("duplicate resume exits cleanly without entering the fitting driver", {
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  source<-file.path(root,"original");output<-file.path(root,"resumed")
  dir.create(source);dir.create(output)
  saveRDS(list(source=normalizePath(source)),file.path(output,"resume-source.rds"))
  lock<-file.path(output,"mechanism","RUN.lock");dir.create(lock,recursive=TRUE)
  owner<-file.path(lock,"owner.rds")
  saveRDS(list(pid=Sys.getpid(),host=unname(Sys.info()[["nodename"]]),token="test"),owner)
  before<-tools::md5sum(owner)
  runner<-function(...)stop("Duplicate invocation must not enter fitting")
  result<-e$study2_resume_competitors(list(action="run",`source-root`=source,output=output),repo,runner,
    check_runtime=function()stop("An active job retains its loaded runtime"))
  testthat::expect_identical(result$status,"already_running")
  testthat::expect_identical(result$pid,Sys.getpid())
  testthat::expect_identical(tools::md5sum(owner),before)
  testthat::expect_null(e$study2_active_competitor_run(output,alive=function(pid)FALSE))
  testthat::expect_null(e$study2_active_competitor_run(output,alive=function(pid)NA))
  testthat::expect_null(e$study2_active_competitor_run(output,alive=function(pid) {
    saveRDS(list(pid=pid,host="another-host"),owner);TRUE
  }))
  testthat::expect_null(e$study2_active_competitor_run(output,alive=function(pid)TRUE))
})

testthat::test_that("new work validates its runtime before copying or fitting", {
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  source<-file.path(root,"original");dir.create(source)
  output<-file.path(root,"resumed")
  testthat::expect_error(e$study2_resume_competitors(
    list(action="run",`source-root`=source,output=output),repo,
    runner=function(...)stop("Unexpected fitting"),
    check_runtime=function()stop("Invalid runtime")),"Invalid runtime")
  testthat::expect_false(dir.exists(output))
})
