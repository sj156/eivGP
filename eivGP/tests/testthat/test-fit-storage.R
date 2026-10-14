test_that("compact storage preserves both engines and exact continuation", {
# Lossless storage, legacy reads, integrity, and exact continuation for both engines.

for(engine in c('univariate','multivariate')) {
 n<-10L;u<-seq(-1,1,length.out=n);X<-matrix(cos(seq_len(n)),ncol=1)
 C<-matrix(as.integer(u>0)+1L,ncol=1);U<-matrix(u,ncol=1);U[3:8,]<-NA_real_
 args<-list(X=X,y=sin(u)+X[,1]/5,C=C,U_obs=U,engine=engine,
   n_iter=16L,burn=4L,n_chains=2L,parallel=FALSE,seed=913L)
 if(engine=='multivariate')args$ident<-'none'
 fit<-suppressWarnings(do.call(fit_eivgp,args))
 compact<-mixedgp_pack_artifact(fit)
 stopifnot(inherits(compact,'mixedgp_compact_fit'),length(compact$pooled_fields)>0,
   identical(mixedgp_unpack_artifact(compact),fit))
 path<-tempfile(fileext='.rds');old<-tempfile(fileext='.rds')
 mixedgp_save_artifact(list(identity=list(seed=913L),fit=fit),path)
 restored<-mixedgp_read_artifact(path)$fit
 saveRDS(fit,old);stopifnot(identical(mixedgp_read_artifact(old),fit),identical(restored,fit))
 a<-suppressWarnings(continue_eivgp(fit,3L,parallel=FALSE))
 b<-suppressWarnings(continue_eivgp(restored,3L,parallel=FALSE))
 stopifnot(identical(a$mcmc,b$mcmc),identical(a$checkpoint,b$checkpoint))
 changed<-fit;key<-paste0('samples_',compact$pooled_fields[1]);changed$mcmc[[key]][1]<-999
 stopifnot(identical(mixedgp_unpack_artifact(mixedgp_pack_artifact(changed)),changed))
 cat(engine,': exact round trip, legacy read, and continuation passed\n')
 unlink(c(path,old))
}

expect_identical(adni_make_C(data.frame(c_ab42_ab40_3=c(1,3),c_gfap_3=c(2,1)))$c_gfap_3, ordered(c(2,1),levels=1:3))
})
