#' Title Do Bayesian variable selection via penalized credible region
#'
#' This is a function using the algorithm doing variable selection via penalized credible
#' interval proposed by Bondell et al. (2012). The computation of the proposed sequence is
#' doing matrix computing and using existing LASSO software.
#'
#'
#' @param dlresult Posterior samples of beta. A large matrix (nmc/thin)*p
#'
#' @return \item{betatil}{Variable selection result of beta, a p*1 vector. Most of the values shrinks to 0}
#' @export
#'
#' @examples {
#' p=30
#' n=5
#' #generate x
#' x=matrix(rnorm(n*p),nrow=n)
#' #generate beta
#' beta=c(rep(0,10),runif(n=5,min=-1,max=1),rep(0,10),runif(n=5,min=-1,max=1))
#' #generate y
#' y=x%*%beta+rnorm(n)
#' hyper=dlhyper(x,y)
#' dlresult=dl(x,y,hyper=hyper)
#' dlvs(dlresult)
#' }

dlvs<-function(dlresult){
  #use penalized credible region to do variable selection
  #matrix computation to make the solutions be accomplished by LASSO
  p=ncol(dlresult)
  betamean<-apply(dlresult,2,mean)
  betacov<-cov(dlresult)
  D<-betamean^2
  cov1=expm::sqrtm(betacov)
  covinv<-MASS::ginv(cov1)
  #scale
  xstar=t(as.numeric(D)*covinv)
  ystar=covinv%*%betamean
  xstar<-scale(xstar, scale=F)
  ystar<-ystar-mean(ystar)
  #solve LASSO problem by glmnet package
  model<-glmnet::glmnet(xstar,ystar,standardize=FALSE,alpha=1,family="gaussian")
  lam=glmnet::cv.glmnet(xstar,ystar,standardize=FALSE,alpha=1,family="gaussian")$lambda.1se
  betastar=coef(model,s=lam)[2:(p+1)]
  betatil=D*betastar
  return(betatil)
}
