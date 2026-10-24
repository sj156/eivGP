skip_if_not_installed("EzGP")
test_that('EzGP corrected covariance handles exact, rounded and distant inputs',{
 X<-cbind(c(-1,0,1,1.5),c(.2,.8,.4,.6),c(1,1,2,2),c(1,2,1,2));p<-2L;q<-2L;m<-c(2L,2L)
 par<-c(.8,.3,.4,.7,.9,seq(.2,.9,length.out=8));tau<-.1
 fit<-structure(list(param=par,data=list(X=X,Y=c(-.3,.7,.2,.8),p=p,q=q,m=m,tau=tau)),class='EzGP model')
 Xt<-rbind(X[1,],X[1,]+c(1e-7,0,0,0),c(.3,.4,1,2))
 ans<-mixedgp_ezgp_predict(Xt,fit)
 # Independent direct conditioning using the joint latent covariance.
 Kall<-EzGP:::cov_m(rbind(X,Xt),p,q,m,7L,par,0)
 K<-Kall[1:4,1:4]+diag(tau,4);cross<-Kall[5:7,1:4];ki<-solve(K);mu<-sum(ki%*%fit$data$Y)/sum(ki)
 expected_mean<-as.numeric(mu+cross%*%ki%*%(fit$data$Y-mu));expected_var<-diag(Kall[5:7,5:7]-cross%*%ki%*%t(cross))
 expect_equal(ans$Y_hat,expected_mean,tolerance=1e-12);expect_equal(ans$MSE,expected_var,tolerance=1e-12)
 expect_equal(ans$predictive_MSE-ans$MSE,rep(tau,3),tolerance=1e-12);expect_true(all(ans$MSE>0))
 old<-EzGP::EzGP_predict(Xt,fit,1);expect_lt(old$MSE[1],0);expect_lt(old$MSE[2],0)
 expect_equal(ans$Y_hat[3],old$Y_hat[3],tolerance=1e-12);expect_equal(ans$MSE[3],old$MSE[3],tolerance=1e-12)
 expect_equal(mixedgp_ezgp_predict(Xt[1,],fit)$Y_hat,ans$Y_hat[1])
 expect_error(mixedgp_ezgp_predict(c(NA,0,1,1),fit),'finite')
})
test_that('training CV remains finite at duplicate train-validation inputs without optimization',{
 env<-new.env(parent=environment(select_ezgp_tau_cv))
 env$fit_ezgp_at_tau<-function(train_x,y_train,p,q,m_vec,tau,maxeval)
   structure(list(param=c(.8,.3,.7,.4,.9),data=list(X=train_x,Y=y_train,p=p,q=q,m=m_vec,tau=tau)),class='EzGP model')
 cv<-select_ezgp_tau_cv;environment(cv)<-env
 ans<-cv(cbind(rep(.5,8),rep(1,8)),seq_len(8)/10,1L,1L,2L,42L,c(.1,.2),2L,1L)
 expect_true(all(is.finite(ans$cv_loss)));expect_true(ans$tau%in%ans$tau_grid)
})
