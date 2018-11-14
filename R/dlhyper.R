#' Tune the hyperparameter in the prior distribtuion
#'
#'
#' This function is to tune the value of hyperparameter in the prior, which can
#' be [1/max(n,p),1/2]. We use the method proposed by Zhang et al. (2018). This method tune
#'  the hyperparameter by incorporating a prior on R^2. And they give a direct way to
#'  minimize KL directed divergence for special condition.
#'
#'
#' @param x input matrix, each row is an observation vector, dimension n*p. Same as the argument in dlmain
#' @param y Response variable, a n*1 vector. Same as the argument in dlmain
#'
#' @return \item{hyper}{A value that can use in the following posterior computation}
#'
#' @examples
#' p=60
#' n=8
#' #generate x
#'x=matrix(rnorm(n*p),nrow=n)
#' #generate beta
#' beta=c(rep(0,10),runif(n=5,min=-1,max=1),rep(0,20),runif(n=5,min=-1,max=1),rep(0,p-40))
#' #generate y
#' y=x%*%beta+rnorm(n)
#' hyper=dlhyper(x,y)
#'
#'
#' @export
#'
dlhyper<-function(x,y){
  p=ncol(x)
  n=nrow(x)
  #calculate hyperparameter
  xtx=t(x)%*%x
  d=eigen(xtx/n)$values
  P=sum(d)
  Q=4*sum(d^2)-sum(d)^2
  R=-sum(d)^3
  C=P^2/9-Q/3
  A=P*Q/6-P^3/27-R/2
  B=A^2-C^3
  hyper=sqrt(2/((A+sqrt(B))^(1/3)+sign(A-sqrt(B))*abs(A-sqrt(B))^(1/3)-P/3))
  hyper[hyper<1/p]=1/p
  return(hyper)
}



