test_that("checkpoint storage contains generated rows only and reads legacy arrays", {
  draws<-list(U=array(seq_len(120),c(10,4,3)),tau=matrix(seq_len(20),10,2),
    v=seq_len(10),unused=NULL)
  for(n in c(0L,3L,10L)) {
    prefix<-mixedgp_checkpoint_prefix(draws,n)
    expect_equal(dim(prefix$U)[1],n)
    restored<-mixedgp_checkpoint_restore_draws(prefix,draws)
    expect_identical(restored,draws)
  }
  expect_identical(mixedgp_checkpoint_restore_draws(draws,draws),draws)
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  x<-list(draws=list(U=array(NA_real_,c(10000,40,2))),saved=0L)
  saveRDS(x,file.path(root,'old.rds'),compress=FALSE)
  x$draws<-mixedgp_checkpoint_prefix(x$draws,0L)
  mixedgp_checkpoint_write(x,file.path(root,'new.rds'))
  expect_identical(readRDS(file.path(root,'new.rds')),x)
  expect_lt(file.info(file.path(root,'new.rds'))$size,file.info(file.path(root,'old.rds'))$size/100)
})

test_that("verified full fits retire restart copies, preserving continuation and timing", {
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  old<-options(eivGP.checkpoint_root=file.path(root,'checkpoints'),eivGP.checkpoint_iterations=3L)
  on.exit(options(old),add=TRUE)
  n<-10L;u<-seq(-1,1,length.out=n);X<-matrix(cos(seq_len(n)),ncol=1)
  C<-matrix(as.integer(u>0)+1L,ncol=1);U<-matrix(u,ncol=1);U[3:8,]<-NA_real_
  fit<-suppressWarnings(fit_eivgp(X=X,y=sin(u)+X[,1]/5,C=C,U_obs=U,engine='univariate',
    n_iter=12L,burn=4L,n_chains=2L,parallel=FALSE,seed=913L))
  disk<-attr(fit,'eivgp_restart_storage')$directory
  expect_true(dir.exists(disk))
  chains<-function()list.files(disk,'^chain-[0-9]+[.]rds$',full.names=TRUE)
  expect_length(chains(),2L)
  # Neither prediction-only output, a truncated fit, nor a mismatched file qualifies.
  path<-file.path(root,'final.rds')
  mixedgp_save_artifact(list(states=list(1)),path)
  expect_length(chains(),2L)
  partial<-fit;partial$mcmc$samples_by_chain$u[[1]]<-partial$mcmc$samples_by_chain$u[[1]][1,,drop=FALSE]
  mixedgp_save_artifact(partial,path)
  expect_length(chains(),2L)
  expect_warning(mixedgp_retire_saved_checkpoints(fit,path),'read-back')
  expect_length(chains(),2L)
  dir.create(file.path(disk,'LOCK'))
  mixedgp_save_artifact(fit,path)
  expect_length(chains(),2L)
  unlink(file.path(disk,'LOCK'),recursive=TRUE)
  mixedgp_save_artifact(list(fit=fit),path)
  restored<-mixedgp_read_artifact(path)$fit
  expect_identical(restored,fit)
  expect_length(chains(),0L)
  expect_true(file.exists(file.path(disk,'retired.rds')))
  expect_true(file.exists(file.path(disk,'chain-01-timing.csv')))
  expect_false(file.exists(file.path(disk,'identity.rds')))
  options(eivGP.checkpoint_root=NULL)
  a<-suppressWarnings(continue_eivgp(fit,3L,parallel=FALSE))
  b<-suppressWarnings(continue_eivgp(restored,3L,parallel=FALSE))
  expect_identical(a$mcmc,b$mcmc)
  expect_identical(a$checkpoint,b$checkpoint)
})

test_that("mechanism writer retires multivariate restart copies after full-fit publication", {
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  old<-options(eivGP.checkpoint_root=file.path(root,'checkpoints'))
  on.exit(options(old),add=TRUE)
  n<-10L;u<-seq(-1,1,length.out=n);X<-matrix(cos(seq_len(n)),ncol=1)
  C<-matrix(as.integer(u>0)+1L,ncol=1);U<-matrix(u,ncol=1);U[3:8,]<-NA_real_
  fit<-suppressWarnings(fit_eivgp(X=X,y=sin(u)+X[,1]/5,C=C,U_obs=U,engine='multivariate',
    ident='none',n_iter=12L,burn=4L,n_chains=2L,parallel=FALSE,seed=913L))
  disk<-attr(fit,'eivgp_restart_storage')$directory
  api<-eivgp_numerical_api(tempdir())
  api$atomic_rds(fit,file.path(root,'conditional-001.rds'))
  expect_length(list.files(disk,'^chain-[0-9]+[.]rds$'),0L)
  expect_identical(mixedgp_read_artifact(file.path(root,'conditional-001.rds')),fit)
})
