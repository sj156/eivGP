repo<-Sys.getenv('EIVGP_TEST_REPO');stopifnot(nzchar(repo))
source(file.path(repo,'tools/ocean-tools.R'),local=TRUE)

test_that('Ocean method ownership permits distinct groups and excludes duplicate writers',{
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  set.seed(719);before<-.Random.seed
  ocean_execution_lock(root,'eivGP',function(){
    expect_error(ocean_execution_lock(root,'eivGP',function()NULL),'lock is held')
    expect_error(ocean_execution_lock(root,'report',function()NULL),'lock is held')
    ocean_execution_lock(root,'competitors',function(){
      expect_true(all(file.exists(file.path(root,'_meta/locks',paste0(c('eivGP','competitors'),'.lock/owner.json')))))
    })
  })
  expect_identical(.Random.seed,before)
  expect_false(any(dir.exists(file.path(root,'_meta/locks',paste0(c('eivGP','competitors'),'.lock')))))
  ocean_execution_lock(root,'competitors',function(){
    expect_error(ocean_execution_lock(root,'report',function()NULL),'lock is held')
    expect_false(dir.exists(file.path(root,'_meta/locks/eivGP.lock')))
    expect_error(ocean_migrate_storage(root,paste0(root,'_migrated')),'locked')
  })
  expect_error(ocean_execution_lock(root,'eivGP',function()stop('worker error')),'worker error')
  expect_false(dir.exists(file.path(root,'_meta/locks/eivGP.lock')))
  dir.create(file.path(root,'.application.lock'))
  expect_error(ocean_execution_lock(root,'competitors',function()NULL),'exclusive workflow lock')
})

test_that('Ocean executor updates retain the original fitting identity and reject numerical changes',{
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  path<-file.path(root,'workflow_identity.rds')
  saved<-list(schema='ocean-package-CC-v2',application='previous',data='same',code='same',
    packages='same',numerics=list(chains=4L,iter=50000L,burn=5000L))
  saveRDS(saved,path);hash<-tools::md5sum(path)
  current<-saved;current$application<-'current'
  expect_identical(ocean_workflow_identity(path,current),saved)
  expect_identical(tools::md5sum(path),hash)
  current$numerics$burn<-500L
  expect_error(ocean_workflow_identity(path,current),'incompatible')
  expect_identical(tools::md5sum(path),hash)
})
