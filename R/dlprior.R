#' Title Simulate the dirichlet laplace shrinkage prior
#'
#' This function generates random deviates from dirichlet laplace shrinkage prior and can plot the distribution function.
#'
#'
#' @param hyper important hyperparameter that related to posterior shrinkage scales and prior distribution
#' @param p number of observations
#' @param plt whether to plot the dirichlet laplace prior. default TRUE means plot the distribution
#' @param min left point of the plot graph
#' @param max right point of the plot graph
#' @param sigma the value equals to normal noises' standard deviations
#'
#' @return \item{beta}{A p*1 vector. p observations from the distribution}
#' @export
#'
#' @examples {theta=dlprior(hyper=1/2,p=100000,plt=TRUE,min=-5,max=5,sigma=1)}
#'
#'
#'
dlprior<-function(hyper=1/2,p=100000,plt=TRUE,min=-5,max=5,sigma=1){
  #dirichlet-laplace
  #prior for beta: beta[j] ~ N(0, sigma^2*psi[j]*phi[j]^2*tau^2)
  #phi ~ Dir(a,..., a)
  #tau ~ gamma(p*a, 1/2)
  #psi[j] iid~ exp(1/2)
  a=rep(hyper,p)
  psi=stats::rexp(p,1/2)
  phi=LaplacesDemon::rdirichlet(1,a)
  tau=stats::rgamma(1,shape=p*a,rate=1/2)
  beta=rep(1,p)
  for(j in 1:p){beta[j]=rnorm(1,sd=sigma*sqrt(psi[j]*(phi[j]^2)*(tau^2)))}
  if(plt==TRUE){
  d1=density(beta,from=min,to=max)
  plot(d1)
  }
  return(beta)
}







