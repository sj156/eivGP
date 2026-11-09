test_that("follow-up lock recovery preserves live owners and archives dead ones", {
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  lock<-file.path(root,"job.lock");dir.create(lock)
  writeLines("99999999",file.path(lock,"pid"))
  before<-tools::md5sum(file.path(lock,"pid"))
  expect_error(mixedgp_followup_lock(lock,alive=function(pid)TRUE),"active")
  expect_identical(tools::md5sum(file.path(lock,"pid")),before)
  expect_error(mixedgp_followup_lock(lock,alive=function(pid)NA),"cannot be verified")
  release<-mixedgp_followup_lock(lock,alive=function(pid)FALSE)
  archive<-list.dirs(root,recursive=FALSE,full.names=TRUE)
  archive<-archive[grepl("stale-",archive)]
  expect_length(archive,1L)
  expect_identical(readLines(file.path(archive,"pid")),"99999999")
  expect_true(mixedgp_followup_process_alive(Sys.getpid()))
  # Cleanup from a previous owner must not erase a successor's lock.
  owner<-readRDS(file.path(lock,"owner.rds"));owner$token<-"replacement"
  saveRDS(owner,file.path(lock,"owner.rds"));release()
  expect_true(dir.exists(lock))
  owner$host<-paste0(Sys.info()[["nodename"]],"-foreign")
  saveRDS(owner,file.path(lock,"owner.rds"))
  expect_error(mixedgp_followup_lock(lock,alive=function(pid)FALSE),"another host")
  owner$host<-unname(Sys.info()[["nodename"]]);saveRDS(owner,file.path(lock,"owner.rds"))
  expect_error(mixedgp_followup_lock(lock,alive=function(pid) {
    owner$token<-"changed-during-check";saveRDS(owner,file.path(lock,"owner.rds"));FALSE
  }),"ownership changed")
  release<-mixedgp_followup_lock(lock,alive=function(pid)FALSE)
  release();expect_false(dir.exists(lock))
})

test_that("audited aliases retain checkpoint/cache keys and strict scientific identity", {
  compat<-data.frame(field=c("helpers.R","helpers.R","mechanism","mechanism"),
    old=c("h0","h1","m0","m1"),new=c("h1","h2","m1","m2"))
  old<-list(package_manifest=c(helpers.R="h0",mcmc.R="sampler"),code=c(mechanism="m0"),
    settings=list(iter=100,seed=5),versions=c(eivGP="1"))
  new<-old;new$package_manifest["helpers.R"]<-"h2";new$code["mechanism"]<-"m2"
  expect_true(mixedgp_followup_identity_equal(old,new,compat))
  expect_identical(mixedgp_followup_normalize(new,compat),old)
  for(field in c("iter","seed")) {
    bad<-new;bad$settings[[field]]<-0
    expect_false(mixedgp_followup_identity_equal(old,bad,compat))
  }
  bad<-new;bad$package_manifest["mcmc.R"]<-"changed"
  expect_false(mixedgp_followup_identity_equal(old,bad,compat))
  bad<-new;bad$versions["eivGP"]<-"2"
  expect_false(mixedgp_followup_identity_equal(old,bad,compat))
  bad<-new;bad$package_manifest["helpers.R"]<-"unknown"
  expect_false(mixedgp_followup_identity_equal(old,bad,compat))
})

test_that("cached CC predictions resume scoring without fitting or changing predictions", {
  api<-eivgp_numerical_api(tempdir())
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  dat<-list(train=list(C=matrix(1:2)),test=list(C=matrix(1:3),y=c(-1,0,1)))
  pred<-list(y_mean=matrix(0,3,2),y_var=matrix(1,3,2),m_point=rep(0,3))
  saveRDS(pred,file.path(root,"predictions.rds"))
  saveRDS(list(pass=TRUE),file.path(root,"diagnostics.rds"))
  truth<-file.path(root,"truth.rds");saveRDS(list(m=rep(0,3)),truth)
  saveRDS(list(message="interrupted"),file.path(root,"failure.rds"))
  api$measurement_cut<-api$cc_evaluate<-function(...)stop("No fitting allowed")
  before<-tools::md5sum(file.path(root,c("predictions.rds","diagnostics.rds")))
  expect_false(api$score_cached_cc(root,dat,1L,truth,"FULL"))
  expect_true(api$score_cached_cc(root,dat,1L,truth,"UC-GP-CC"))
  expect_identical(tools::md5sum(names(before)),before)
  expect_false(file.exists(file.path(root,"failure.rds")))
  expect_equal(nrow(readRDS(file.path(root,"scores.rds"))),6L)
})
