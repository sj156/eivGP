cut_storage_fixture <- function(root,retain=TRUE) {
  api<-eivgp_numerical_api(tempdir())
  u<-seq(-1,1,length.out=10);U<-cbind(u,rev(u));U[3:8,]<-NA_real_
  C<-cbind(rep(1:4,length.out=10),rep(4:1,length.out=10))
  cfg<-list(iter=16L,warmup=4L,outer=4L,cores=1L,retain_draws=retain)
  object<-suppressWarnings(api$measurement_cut(C,U,cfg,78121L))
  list(api=api,U=U,C=C,cfg=cfg,object=object,
    disk=attr(object,"eivgp_cut_restart_storage")$directory)
}
test_that("CUT compact publication retires exact checkpoints without changing scores", {
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  old<-options(eivGP.checkpoint_root=file.path(root,"checkpoints"),eivGP.checkpoint_iterations=2L,
    eivGP.checkpoint_after_seconds=0);on.exit(options(old),add=TRUE)
  z<-cut_storage_fixture(root);full<-z$object
  compact<-mixedgp_compact_cut_measurement(full,z$U)
  expect_null(compact$trace);expect_identical(compact$states,full$states)
  expect_identical(compact$diagnostics,full$diagnostics)
  expect_identical(compact$checkpoint,full$checkpoint)
  # Compare the previous scoring path against compact summaries using identical predictions.
  mis<-3:8;arr<-full$trace[,,1:12,drop=FALSE]
  U<-array(NA_real_,c(dim(arr)[1]*4,10,2));U[,mis,]<-matrix(arr,ncol=12)
  dat<-list(train=list(X=matrix(1,10,1),C=z$C,U=cbind(seq(-1,1,length.out=10),seq(1,-1,length.out=10))),
    test=list(C=z$C[1:2,,drop=FALSE],y=c(.1,.2)))
  pred<-list(y_mean=matrix(c(.1,.2),2,1),y_var=matrix(1,2,1),m_point=c(0,0),U=U)
  a<-z$api$score_result(pred,dat,c(1L,2L,9L,10L),list(m=c(0,0)),"CUT")
  pred$U<-NULL;pred$U_summary<-compact$imputation
  b<-z$api$score_result(pred,dat,c(1L,2L,9L,10L),list(m=c(0,0)),"CUT")
  expect_identical(a,b)
  path<-file.path(root,"measurement.rds")
  z$api$atomic_rds(compact,path)
  expect_identical(readRDS(path),compact)
  expect_length(list.files(z$disk,"^chain-[0-9]+[.]rds$"),0L)
  expect_equal(readRDS(file.path(z$disk,"retired.rds"))$status,"removed")
  expect_true(file.exists(file.path(z$disk,"chain-01-timing.csv")))
  # Full trace retention is an opt-in output policy, independent of completed checkpoints.
  options(eivGP.checkpoint_root=NULL)
  rerun<-suppressWarnings(z$api$measurement_cut(z$C,z$U,modifyList(z$cfg,list(retain_draws=FALSE)),78121L))
  expect_null(rerun$trace);expect_identical(rerun$imputation,compact$imputation)
  expect_identical(rerun$states,compact$states)
  # A representative expanded trace gives a useful on-disk size regression check.
  big<-full;big$trace<-full$trace[rep(seq_len(dim(full$trace)[1]),100),,,drop=FALSE]
  big$imputation<-NULL
  small<-mixedgp_compact_cut_measurement(big,z$U)
  saveRDS(big,file.path(root,"full.rds"),compress=FALSE)
  saveRDS(small,file.path(root,"compact.rds"),compress=FALSE)
  expect_lt(file.info(file.path(root,"compact.rds"))$size,file.info(file.path(root,"full.rds"))$size/10)
})
test_that("CUT cleanup preserves locked, unfinished, mismatched and partial files", {
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  old<-options(eivGP.checkpoint_root=file.path(root,"checkpoints"),eivGP.checkpoint_iterations=2L,
    eivGP.checkpoint_after_seconds=0);on.exit(options(old),add=TRUE)
  z<-cut_storage_fixture(root,FALSE);path<-file.path(root,"measurement.rds")
  files<-file.path(z$disk,sprintf("chain-%02d.rds",1:4));raw<-readRDS(files[1])
  saveRDS(z$object,path)
  dir.create(file.path(z$disk,"LOCK"));mixedgp_retire_cut_measurement(z$object,path)
  expect_true(all(file.exists(files)));unlink(file.path(z$disk,"LOCK"),recursive=TRUE)
  broken<-raw;broken$iteration<-5L;saveRDS(broken,files[1])
  expect_warning(mixedgp_retire_cut_measurement(z$object,path),"unfinished")
  expect_true(all(file.exists(files)))
  broken<-raw;broken$trace[1,1]<-999;saveRDS(broken,files[1])
  expect_warning(mixedgp_retire_cut_measurement(z$object,path),"evidence_mismatch")
  expect_true(all(file.exists(files)))
  saveRDS(raw,files[1]);extra<-file.path(z$disk,"chain-05.rds");saveRDS(raw,extra)
  expect_warning(mixedgp_retire_cut_measurement(z$object,path),"incomplete_or_extra")
  unlink(extra)
  other<-z$object;other$imputation$mean[1]<-999
  expect_warning(mixedgp_retire_cut_measurement(other,path),"read-back")
  expect_true(all(file.exists(files)))
  identity<-file.path(z$disk,"identity.rds");id<-readRDS(identity)
  badid<-id;badid$model$seed<-0L;saveRDS(badid,identity)
  expect_warning(mixedgp_retire_cut_measurement(z$object,path),"identity changed")
  saveRDS(id,identity)
  # Cached completed jobs retry retirement; existing scores are untouched.
  folder<-file.path(root,"mechanism","rep001-cal20-CUT");dir.create(folder,recursive=TRUE)
  shared<-file.path(dirname(folder),"shared-rep001-cal20");dir.create(shared)
  file.copy(path,file.path(shared,"measurement.rds"))
  saveRDS(data.frame(x=1),file.path(folder,"scores.rds"));hash<-tools::md5sum(file.path(folder,"scores.rds"))
  mixedgp_cleanup_cut_cache(folder,FALSE)
  expect_length(list.files(z$disk,"^chain-[0-9]+[.]rds$"),0L)
  expect_identical(tools::md5sum(file.path(folder,"scores.rds")),hash)
})
test_that("legacy full measurement checkpoints can be audited and compacted explicitly", {
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  old<-options(eivGP.checkpoint_root=file.path(root,"checkpoints"),eivGP.checkpoint_iterations=2L,
    eivGP.checkpoint_after_seconds=0);on.exit(options(old),add=TRUE)
  z<-cut_storage_fixture(root);legacy<-z$object
  legacy$imputation<-legacy$checkpoint<-legacy$provenance<-NULL
  attr(legacy,"eivgp_cut_restart_storage")<-NULL
  mech<-file.path(root,"mechanism");shared<-file.path(mech,"shared-rep001-cal20")
  dir.create(shared,recursive=TRUE);path<-file.path(shared,"measurement.rds");saveRDS(legacy,path)
  saveRDS(list(settings=list(seed=77101L)),file.path(mech,"protocol.rds"))
  audit<-function(dry=TRUE,compact=FALSE) {
    a<-tempfile();on.exit(unlink(a,recursive=TRUE))
    suppressMessages(prune_study_checkpoints(root,dry_run=dry,audit_dir=a,compact_measurement=compact))
  }
  hash<-tools::md5sum(path)
  result<-audit();expect_equal(result$status,"verified");expect_identical(tools::md5sum(path),hash)
  bad<-legacy;bad$trace[1,1,1]<-999;saveRDS(bad,path)
  result<-audit(FALSE);expect_equal(result$status,"preserved");expect_match(result$reason,"trace_mismatch")
  saveRDS(legacy,path)
  result<-audit(FALSE,TRUE)
  expect_equal(result$status,"removed");expect_equal(result$chains_removed,4L)
  saved<-readRDS(path);expect_null(saved$trace);expect_identical(saved$states,legacy$states)
  expect_identical(saved$imputation,z$object$imputation)
  expect_equal(sum(audit(FALSE,TRUE)$removed_bytes),0)
})
test_that("conditional fit retirement requires verified selected states and respects opt-in", {
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  api<-eivgp_numerical_api(tempdir());path<-file.path(root,"conditional-001.rds");sp<-sub(".rds$","-states.rds",path)
  fit<-list(diagnostics=list(x=1));saveRDS(fit,path)
  api$select_gp<-function(fit,number,measurement)list(list(x=1))
  saved<-list(storage_schema=1L,states=list(list(x=1)),pass=FALSE,diagnostics=fit$diagnostics,
    source_md5=unname(tools::md5sum(path)),selection=list(number=1L,measurement=NULL))
  saveRDS(saved,sp)
  api$retire_conditional_fit(path,sp,TRUE);expect_true(file.exists(path))
  bad<-saved;bad$source_md5<-"wrong";saveRDS(bad,sp)
  expect_error(api$retire_conditional_fit(path,sp),"changed");expect_true(file.exists(path))
  bad<-saved;bad$states[[1]]$x<-999;saveRDS(bad,sp)
  expect_error(api$retire_conditional_fit(path,sp),"states or diagnostics differ")
  expect_true(file.exists(path))
  saveRDS(saved,sp);expect_identical(api$retire_conditional_fit(path,sp),saved)
  expect_false(file.exists(path));expect_identical(readRDS(sp),saved)
  expect_identical(api$retire_conditional_fit(path,sp),saved)
})
