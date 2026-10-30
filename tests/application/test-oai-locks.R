# Author: Xinye Chen. OAI group ownership and immutable fitting provenance.
repo<-Sys.getenv('EIVGP_TEST_REPO');stopifnot(nzchar(repo))
source(file.path(repo,'tools/oai.R'));api<-oai_load_api(repo)
testthat::test_that('different groups coexist while duplicate writers and combined reports are blocked',{
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  set.seed(19);rng<-.Random.seed
  api$oai_execution_lock(root,'eivGP',function(){
    testthat::expect_error(api$oai_execution_lock(root,'eivGP',function()NULL),'OAI lock is held')
    testthat::expect_error(api$oai_execution_lock(root,'report',function()NULL),'OAI lock is held')
    api$oai_execution_lock(root,'competitors',function(){
      testthat::expect_true(all(file.exists(file.path(root,'_meta/locks',paste0(c('eivGP','competitors'),'.lock/owner.json')))))
    })
  })
  testthat::expect_identical(.Random.seed,rng)
  testthat::expect_false(any(dir.exists(file.path(root,'_meta/locks',paste0(c('eivGP','competitors'),'.lock')))))
  testthat::expect_error(api$oai_execution_lock(root,'eivGP',function()stop('worker error')),'worker error')
  testthat::expect_false(dir.exists(file.path(root,'_meta/locks/eivGP.lock')))
  dir.create(file.path(root,'.application.lock'))
  testthat::expect_error(api$oai_execution_lock(root,'competitors',function()NULL),'earlier exclusive executor')
})
testthat::test_that('execution-only upgrades preserve old receipts and reject numerical changes',{
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  old<-list(schema='adni-upstream-screening-v2',application=list(old='abc'),modules='same',burn=5000L,
    stages=50000L,chains=4L,inputs='fixed',plan='oai')
  path<-file.path(root,'workflow_identity.rds');saveRDS(old,path)
  before<-tools::md5sum(path);current<-old;current$schema<-'oai-upstream-screening-v2';current$application<-list(oai='def')
  testthat::expect_identical(api$oai_workflow_identity(path,current),old)
  testthat::expect_identical(tools::md5sum(path),before)
  current$burn<-500L
  testthat::expect_error(api$oai_workflow_identity(path,current),'settings, package modules or inputs differ')
  testthat::expect_identical(tools::md5sum(path),before)
})
