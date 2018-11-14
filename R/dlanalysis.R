#' dlanalysis
#'
#' This is a function that analyse the MCMC sampling result by computing the posterior mean,
#' median and credible intervals
#'
#'
#' @param dlresult Posterior samples of beta. A large matrix (nmc/thin)*p
#' @param alpha Level for the credible intervals. For example,the default is alpha = 0.05 means 95\% credible intervals
#'
#' @return  \item{betamean}{Posterior mean of beta, a p*1 vector.}
#'            \item{LeftCI}{The left bounds of the credible intervals.}
#'            \item{RightCI}{The right bounds of the credible intervals.}
#'            \item{betamedian}{Posterior median of Beta, a p*1 vector.}
#' @export
#'
#' @examples
#' p=50
#' n=5
#' #generate x
#' x=matrix(rnorm(n*p),nrow=n)
#' #generate beta
#' beta=c(rep(0,10),runif(n=5,min=-1,max=1),rep(0,10),runif(n=5,min=-1,max=1),rep(0,p-30))
#' #generate y
#' y=x%*%beta+rnorm(n)
#' hyper=dlhyper(x,y)
#' dlresult=dl(x,y,hyper=hyper)
#' da=dlanalysis(dlresult,alpha=0.05)
#' da$betamean
#' da$betamedian
#' da$LeftCI
#' da$RightCI
#'
#'


dlanalysis<-function(dlresult,alpha=0.05){
  betamean=apply(dlresult,2,mean)
  betamedian=apply(dlresult,2,median)
  min=apply(dlresult,2,quantile,prob=alpha/2)
  max=apply(dlresult,2,quantile,prob=1-alpha/2)
  result=list("betamean"=betamean,"betamedian"=betamedian,"LeftCI"=min,"RightCI"=max)
  return(result)
}
