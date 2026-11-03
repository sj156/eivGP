test_that('Rfast squared differences retain ordering and scalar edge cases', {
 set.seed(241)
 for(n in c(1L,2L,8L,20L,100L)) {
  x<-rnorm(n);y<-rnorm(23L)
  expect_identical(mixedgp_squared_differences(x),outer(x,x,'-')^2)
  expect_identical(mixedgp_squared_differences(x,y),outer(x,y,'-')^2)
 }
})
test_that('distance caches and affected rows match the original covariance', {
 set.seed(99);n<-25L
 ctx<-list2env(list(X=matrix(rnorm(n*2),n),n=n,p=2L,d=2L,kernel='se',matern_nu=2.5))
 state<-list(U=matrix(rnorm(n*2),n),logtheta_x=log(c(.4,.8)))
 dec<-list(theta_u=c(.2,.9))
 reference<-function(st,dd){
  Z<-cbind(ctx$X,st$U);D<-matrix(0,n,n);w<-c(exp(st$logtheta_x),dd$theta_u)
  for(j in seq_along(w))D<-D+w[j]*outer(Z[,j],Z[,j],'-')^2
  kernel_from_weighted_sqdist(D,ctx$kernel,ctx$matern_nu)
 }
 for(kernel in c('se','matern'))for(nu in c(.5,1.5,2.5,1.2)) {
  ctx$kernel<-kernel;ctx$matern_nu<-nu
  for(change in c('initial','dictionary','U','theta','return')) {
   st<-state;dd<-dec
   if(change=='dictionary')dd$theta_u[1]<-.6
   if(change=='U')st$U[3,]<-st$U[3,]+.2
   if(change=='theta')st$logtheta_x[1]<-log(.5)
   expected<-reference(st,dd)
   expect_identical(mixedgp_v030_corr(st,dd,ctx),expected)
   for(rows in list(3L,c(7L,2L,19L)))
    expect_identical(mixedgp_v030_corr(st,dd,ctx,rows),expected[rows,,drop=FALSE])
  }
 }
})
