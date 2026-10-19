#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
static double ladd(double a,double b){if(a==-INFINITY)return b;if(b==-INFINITY)return a;double m=std::max(a,b);return m+log(exp(a-m)+exp(b-m));}
static double ldiff(double a,double b){return a+log(-expm1(b-a));}
static double interval(double a,double b){return a>=0?ldiff(R::pnorm(a,0,1,false,true),R::pnorm(b,0,1,false,true)):ldiff(R::pnorm(b,0,1,true,true),R::pnorm(a,0,1,true,true));}
struct Target {
 arma::mat X,U,A,lo,hi; arma::vec y,txmean,txsd,tu,lm,lower,upper;arma::uvec mis;
 double av,bv,ar,br,nu;bool threshold,se;int n,p,d,m;
 Target(List x):X(as<arma::mat>(x["X"])),U(as<arma::mat>(x["U"])),A(as<arma::mat>(x["A"])),lo(as<arma::mat>(x["measurement_lo"])),hi(as<arma::mat>(x["measurement_hi"])),y(as<arma::vec>(x["y"])),txmean(as<arma::vec>(x["txmean"])),txsd(as<arma::vec>(x["txsd"])),tu(as<arma::vec>(x["tu"])),lm(as<arma::vec>(x["logmass"])),lower(as<arma::vec>(x["lo"])),upper(as<arma::vec>(x["hi"])),mis(as<arma::uvec>(x["mis"])),av(as<double>(x["av"])),bv(as<double>(x["bv"])),ar(as<double>(x["ar"])),br(as<double>(x["br"])),nu(as<double>(x["nu"])),threshold(as<bool>(x["threshold"])),se(as<bool>(x["se"])) {n=X.n_rows;p=X.n_cols;d=U.n_cols;m=mis.n_elem;}
 bool eval(const arma::vec&q,double&lp,arma::vec&g,arma::mat&u){
  if(!q.is_finite())return false;u=U;arma::vec du(m,arma::fill::ones);
  for(int j=0;j<m;j++){
   double v=q[j];if(threshold){double pl=ladd(R::pnorm(lower[j],0,1,true,true),lm[j]+R::pnorm(v,0,1,true,true));double pr=ladd(R::pnorm(upper[j],0,1,false,true),lm[j]+R::pnorm(v,0,1,false,true));v=pl<=pr?R::qnorm(pl,0,1,true,true):R::qnorm(pr,0,1,false,true);if(!std::isfinite(v)||v<=lower[j]||v>=upper[j])return false;du[j]=exp(lm[j]+.5*(v*v-q[j]*q[j]));}u[mis[j]]=v;
  }
  arma::vec tx=arma::exp(txmean+txsd%q.subvec(m,m+p-1));double r=R::plogis(q[m+p],0,1,true,false);if(!tx.is_finite()||r<=0||r>=1)return false;
  arma::mat D(n,n,arma::fill::zeros);std::vector<arma::mat> xd,ud;
  for(int j=0;j<p;j++){arma::mat z=arma::repmat(X.col(j),1,n)-arma::repmat(X.col(j).t(),n,1);xd.push_back(arma::square(z));D+=tx[j]*xd[j];}
  for(int j=0;j<d;j++){arma::mat z=arma::repmat(u.col(j),1,n)-arma::repmat(u.col(j).t(),n,1);ud.push_back(z);D+=tu[j]*arma::square(z);}
  arma::mat K,F;if(se){K=arma::exp(-D);F=-K;}else if(nu==1.5){arma::mat z=arma::sqrt(3*D);K=(1+z)%arma::exp(-z);F=-1.5*arma::exp(-z);}else{arma::mat z=arma::sqrt(5*D);K=(1+z+arma::square(z)/3)%arma::exp(-z);F=-(5.0/6)*(1+z)%arma::exp(-z);}
  arma::mat B=r*K;B.diag()+=1-r;arma::mat L;if(!arma::chol(L,B,"lower"))return false;
  arma::vec v=arma::solve(arma::trimatl(L),y,arma::solve_opts::fast),alpha=arma::solve(arma::trimatu(L.t()),v,arma::solve_opts::fast);
  double rate=bv+arma::dot(v,v)/2,shape=av+n/2.0;
  lp=-arma::sum(arma::log(L.diag()))-shape*log(rate)-.5*arma::dot(q.head(m+p),q.head(m+p))+ar*R::plogis(q[m+p],0,1,true,true)+br*R::plogis(q[m+p],0,1,false,true);
  arma::mat Li=arma::solve(arma::trimatl(L),arma::eye(n,n),arma::solve_opts::fast),G=.5*(shape/rate*(alpha*alpha.t())-Li.t()*Li),gu(n,d);
  for(int k=0;k<d;k++)gu.col(k)=4*tu[k]*r*arma::sum(G%F%ud[k],1);
  if(!threshold){for(int i=0;i<n;i++)for(unsigned j=0;j<A.n_rows;j++){
   double mu=arma::dot(u.row(i),A.row(j)),a=lo(i,j)-mu,b=hi(i,j)-mu,lz=interval(a,b);if(!std::isfinite(lz))return false;lp+=lz;
   double la=R::dnorm(a,0,1,true),lb=R::dnorm(b,0,1,true),score;
   if(la==lb)score=0;else if(la>lb)score=exp(la-lz)*(-expm1(lb-la));else score=-exp(lb-lz)*(-expm1(la-lb));
   gu.row(i)+=score*A.row(j);
  }}
  g.set_size(m+p+1);for(int j=0;j<m;j++)g[j]=gu[mis[j]]*du[j]-q[j];
  for(int j=0;j<p;j++)g[m+j]=r*tx[j]*txsd[j]*arma::accu(G%F%xd[j])-q[m+j];
  g[m+p]=r*(1-r)*(arma::accu(G%K)-arma::trace(G))+ar*(1-r)-br*r;
  return std::isfinite(lp)&&g.is_finite();
 }
};
// [[Rcpp::export]]
List mixedgp_hmc_target_cpp(arma::vec q,List data){Target t(data);double lp=-INFINITY;arma::vec g;arma::mat u;bool ok=t.eval(q,lp,g,u);return List::create(_["logp"]=ok?lp:-INFINITY,_["gradient"]=g,_["U"]=u,_["ok"]=ok);}
// [[Rcpp::export]]
List mixedgp_hmc_leapfrog_cpp(arma::vec q,arma::vec p,List data,double step,int steps,arma::vec inverse_mass){
 Target t(data);double lp;arma::vec g;arma::mat u;bool ok=t.eval(q,lp,g,u);if(!ok)stop("Nonfinite initial HMC target");double h0=-lp+.5*arma::dot(p%p,inverse_mass),maxerr=0;p+=.5*step*g;
 for(int j=0;j<steps;j++){q+=step*(inverse_mass%p);ok=t.eval(q,lp,g,u);if(!ok)break;arma::vec ps=p+.5*step*g;maxerr=std::max(maxerr,std::abs(-lp+.5*arma::dot(ps%ps,inverse_mass)-h0));p+=(j==steps-1?.5:1)*step*g;}
 double dh=ok?-lp+.5*arma::dot(p%p,inverse_mass)-h0:INFINITY;
 return List::create(_["q"]=q,_["p"]=p,_["delta_h"]=dh,_["max_error"]=ok?maxerr:INFINITY);
}
